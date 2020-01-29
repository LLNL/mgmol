// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "HDFrestart.h"
#include "Control.h"
#include "LocalizationRegions.h"
#include "MGmol_MPI.h"
#include "Mesh.h"
#include "hdf_tools.h"
#include "tools.h"

#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

Timer HDFrestart::open_existing_tm_("HDFrestart::open_existing");
Timer HDFrestart::create_file_tm_("HDFrestart::create_file");
Timer HDFrestart::close_file_tm_("HDFrestart::close_file");

const unsigned short MyHDFStrLength = 24;

struct HDF_FixedLengthString
{
    char mystring[MyHDFStrLength];
};

std::string getDatasetName(const std::string& name, const int color)
{
    std::string extension = std::to_string(color);

    std::string datasetname(name);
    if (color < 10) datasetname.append("0");
    if (color < 100) datasetname.append("0");
    if (color < 1000) datasetname.append("0");

    datasetname.append(extension);

    return datasetname;
}

std::string getDatasetName_old(const std::string& name, const int color)
{
    std::string extension = std::to_string(color);

    std::string datasetname(name);
    datasetname.append(extension);

    return datasetname;
}

HDFrestart::~HDFrestart()
{
    if (!closed_) close();

    closeWorkSpace();

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.barrier();

    if (comm_active_ != MPI_COMM_NULL) MPI_Comm_free(&comm_active_);
}

int HDFrestart::close()
{
    if (closed_) return 0;

    close_file_tm_.start();

    Control& ct = *(Control::instance());

    if (ct.verbose > 1 && active_) mgmol_tools::whatisopen(file_id_);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.barrier();

    closed_ = true;

    // Turn off error handling
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

    if (onpe0)
        (*MPIdata::sout) << "~HDFrestart() --- H5Fclose " << filename_
                         << std::endl;
    herr_t err = 0;
    if (active_)
    {
        assert(file_id_ >= 0);
        err = H5Fclose(file_id_);
    }
    else
    {
        assert(file_id_ == -1);
    }

    mmpi.allreduce(&err, 1, MPI_MIN);
    if (err < 0)
    {
        if (onpe0)
            (*MPIdata::serr)
                << "HDFrestart::close(): H5Fclose failed!!!" << std::endl;
        close_file_tm_.stop();
        return err;
    }
    file_id_ = -1;

    close_file_tm_.stop();
    return 0;
}

void HDFrestart::addDateToFilename()
{
    time_t tt;
    time(&tt);

    struct tm* tt1 = gmtime(&tt);

    // year
    filename_.append("_");
    std::string extension = std::to_string(tt1->tm_year - 100);
    filename_.append(extension);

    // day
    filename_.append("_");
    extension = std::to_string(tt1->tm_yday);
    if (tt1->tm_yday < 10)
    {
        filename_.append("00");
    }
    else if (tt1->tm_yday < 100)
    {
        filename_.append("0");
    }
    filename_.append(extension);

    // hour
    filename_.append("_");
    extension = std::to_string(tt1->tm_hour);
    if (tt1->tm_hour < 10)
    {
        filename_.append("0");
    }
    filename_.append(extension);

    // minutes
    filename_.append("_");
    extension = std::to_string(tt1->tm_min);
    if (tt1->tm_min < 10)
    {
        filename_.append("0");
    }
    filename_.append(extension);
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    // make sure all the PEs use the same filename
    std::vector<char> name_buffer(filename_.begin(), filename_.end());
    mmpi.bcast(&name_buffer[0], name_buffer.size());
    const std::string name(name_buffer.begin(), name_buffer.end());
    filename_ = name;
}

void HDFrestart::setActivity()
{
    if (gather_data_x_)
        // activate only PEs on 1st layer in x direction
        active_ = (pes_.my_mpi(0) == 0);
    else
        active_ = true;
}

void HDFrestart::appendTaskNumberToFilename()
{
    int mytask = 0;
    int npes   = 1;
    MPI_Comm_rank(comm_data_, &mytask);
    MPI_Comm_size(comm_data_, &npes);
    filename_.append(".");
    if (mytask < 1000000 && npes > 999999)
    {
        filename_.append("0");
    }
    if (mytask < 100000 && npes > 99999)
    {
        filename_.append("0");
    }
    if (mytask < 10000 && npes > 9999)
    {
        filename_.append("0");
    }
    if (mytask < 1000)
    {
        filename_.append("0");
    }
    if (mytask < 100)
    {
        filename_.append("0");
    }
    if (mytask < 10)
    {
        filename_.append("0");
    }

    std::stringstream oss("");
    oss << mytask;

    filename_.append(oss.str());

    return;
}

void HDFrestart::setupBlocks()
{
    block_[0] = dimsf_[0] / pes_.n_mpi_task(0);
    block_[1] = dimsf_[1] / pes_.n_mpi_task(1);
    block_[2] = dimsf_[2] / pes_.n_mpi_task(2);

    if (gather_data_x_)
    {
        block_[0] = dimsf_[0];
        bsize_    = block_[2] * block_[1] * block_[0] / pes_.n_mpi_task(0);
    }
    else
    {
        bsize_ = block_[2] * block_[1] * block_[0];
    }

    offset_[0] = 0;
    if (use_hdf5p_)
    {
        offset_[1] = pes_.my_mpi(1) * block_[1];
        offset_[2] = pes_.my_mpi(2) * block_[2];
    }
    else
    {
        offset_[1] = 0;
        offset_[2] = 0;
    }
}

void HDFrestart::addMDTime2File(const float run_time)
{
    if (active_)
    {
        std::string attname("MD_time");

        //  Open a dataset attribute.
        hsize_t dims[1]    = { 1 };
        hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
        hid_t attribute_id = H5Acreate2(file_id_, attname.c_str(),
            H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        if (attribute_id < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::addMDTime2File(), Attribute: " << attname
                << " --- H5Acreate failed!!!" << std::endl;
        }
        float attr_data = run_time;
        // if( onpe0 )
        {
            herr_t status
                = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::addMDTime2File(), Attribute: " << attname
                    << " --- H5Awrite failed!!!" << std::endl;
            }
        }
        herr_t status = H5Sclose(dataspace_id);
        if (status < 0)
            (*MPIdata::serr)
                << "HDFrestart::addMDTime2File(), H5Sclose failed!!!"
                << std::endl;
        status = H5Aclose(attribute_id);
        if (status < 0)
            (*MPIdata::serr)
                << "HDFrestart::addMDTime2File(), H5Aclose failed!!!"
                << std::endl;
    }
}

void HDFrestart::add2File(const int data, const std::string& attname)
{
    if (active_)
    {
        //  Open a dataset attribute.
        hsize_t dims[1]    = { 1 };
        hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
        hid_t attribute_id = H5Acreate2(file_id_, attname.c_str(),
            H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        if (attribute_id < 0)
        {
            (*MPIdata::serr) << "HDFrestart::add2File(), Attribute: " << attname
                             << " --- H5Acreate failed!!!" << std::endl;
        }
        int attr_data = data;
        // if( onpe0 )
        {
            herr_t status = H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::add2File(), Attribute: " << attname
                    << " --- H5Awrite failed!!!" << std::endl;
            }
        }
        herr_t status = H5Sclose(dataspace_id);
        if (status < 0)
            (*MPIdata::serr)
                << "HDFrestart::add2File(), H5Sclose failed!!!" << std::endl;
        status = H5Aclose(attribute_id);
        if (status < 0)
            (*MPIdata::serr)
                << "HDFrestart::add2File(), H5Aclose failed!!!" << std::endl;
    }
}

void HDFrestart::addMDstep2File(const int md_step)
{
    add2File(md_step, "MD_step");
}

float HDFrestart::getMDTimeFromFile() const
{
    float run_time            = 0.;
    std::string function_name = "HDFrestart::getMDTimeFromFile()";
    if (onpe0)
    {
        std::string attname("MD_time");
        htri_t exists = H5Aexists(file_id_, attname.c_str());
        if (exists > 0)
        {
            hid_t attribute_id = H5Aopen_name(file_id_, attname.c_str());
            if (attribute_id < 0)
            {
                (*MPIdata::sout)
                    << function_name
                    << " --- WARNING: H5Aopen_name failed!!! --- Attribute: "
                    << attname << std::endl;
            }
            else
            {
                herr_t status
                    = H5Aread(attribute_id, H5T_NATIVE_FLOAT, &run_time);
                if (status < 0)
                {
                    (*MPIdata::serr) << function_name << "--- H5Aread failed!!!"
                                     << std::endl;
                }
                else
                {
                    (*MPIdata::sout) << "MD time for restart file: " << run_time
                                     << std::endl;
                }

                status = H5Aclose(attribute_id);
                if (status < 0)
                {
                    (*MPIdata::serr) << function_name
                                     << " --- H5Aclose failed!!!" << std::endl;
                }
            }
        }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(&run_time, 1);

    return run_time;
}

int HDFrestart::getMDstepFromFile() const { return getFromFile("MD_step"); }

int HDFrestart::getFromFile(const std::string& attname) const
{
    int data                  = 1;
    std::string function_name = "HDFrestart::getFromFile()";
    if (onpe0)
    {
        htri_t exists = H5Aexists(file_id_, attname.c_str());
        if (exists > 0)
        {
            hid_t attribute_id = H5Aopen_name(file_id_, attname.c_str());
            if (attribute_id < 0)
            {
                (*MPIdata::sout)
                    << function_name
                    << " --- WARNING: H5Aopen_name failed!!! --- Attribute: "
                    << attname << std::endl;
            }
            else
            {
                herr_t status = H5Aread(attribute_id, H5T_NATIVE_INT, &data);
                if (status < 0)
                {
                    (*MPIdata::serr) << function_name << "--- H5Aread failed!!!"
                                     << std::endl;
                }
                else
                {
                    (*MPIdata::sout)
                        << attname << " in restart file: " << data << std::endl;
                }

                status = H5Aclose(attribute_id);
                if (status < 0)
                {
                    (*MPIdata::serr) << function_name
                                     << " --- H5Aclose failed!!!" << std::endl;
                }
            }
        }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(&data, 1);

    return data;
}

void HDFrestart::addReleaseNumber2File(const char* release)
{
    if (active_)
    {
        // create type for std::strings of length MyHDFStrLength
        hid_t strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(strtype, MyHDFStrLength);

        std::string attname("MGmol Release");
        //  Open a dataset attribute.
        hsize_t dims[1]    = { 1 };
        hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
        hid_t attribute_id = H5Acreate2(file_id_, attname.c_str(), strtype,
            dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
        if (attribute_id < 0)
        {
            (*MPIdata::serr) << "Attribute: " << attname
                             << " --- H5Acreate failed!!!" << std::endl;
        }
        // if( onpe0 )
        {
            HDF_FixedLengthString t;
            strncpy(t.mystring, release, MyHDFStrLength);
            herr_t status = H5Awrite(attribute_id, strtype, &t);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::addReleaseNumber2File(), Attribute: "
                    << attname << " --- H5Awrite failed!!!" << std::endl;
            }
        }
        herr_t status = H5Sclose(dataspace_id);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::addReleaseNumber2File(), H5Sclose failed!!!"
                << std::endl;
        }
        status = H5Aclose(attribute_id);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::addReleaseNumber2File(), H5Aclose failed!!!"
                << std::endl;
        }
    }
}

// constructor for one layer of PEs writing data
HDFrestart::HDFrestart(const std::string& filename, const pb::PEenv& pes,
    const unsigned gdim[3], const short option_number)
    : pes_(pes), filename_(filename)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    comm_data_ = mmpi.commSameSpin();

    create_file_tm_.start();

    setOptions(option_number);
    verbosity_ = 0;
    closed_    = false;

    //(*MPIdata::sout)<<"HDFrestart::HDFrestart(), filename="<<filename<<endl;
    setActivity();

    count_[0] = count_[1] = count_[2] = 1;

    stride_[0] = stride_[1] = stride_[2] = 1;

    dimsf_[0] = gdim[0];
    dimsf_[1] = gdim[1];
    dimsf_[2] = gdim[2];

    setupBlocks();

    setCommActive();

    // add date to filename if automatic naming strategy
    Control& ct = *(Control::instance());
    if (ct.out_restart_file_naming_strategy)
    {
        addDateToFilename();
        if (use_hdf5p_) filename_.append(".h5");
    }

    if (!use_hdf5p_) // mkdir with filename and create files with task numbers
    {
        // check if dir exists
        if (onpe0)
        {
            struct stat buf;
            int exists  = stat(filename_.c_str(), &buf);
            mode_t mode = (S_IRWXU | S_IRWXG | S_IRWXO);
            if (exists < 0)
            {
                mkdir(filename_.c_str(), mode);
                (*MPIdata::sout) << "Create dir " << filename_ << std::endl;
            }
        }
        filename_.append("/Task");

        appendTaskNumberToFilename();
        // wait to make sure directory is created by task 0
        // before anybody tries to create a file inside that directory
        MPI_Barrier(comm_data_);
        sleep(1);
    }

    if (active_)
    {
        if (onpe0)
            (*MPIdata::sout)
                << "HDFrestart(): create file " << filename_ << std::endl;
        hid_t access_plist
            = H5Pcreate(H5P_FILE_ACCESS); // property list identifier
        if (use_hdf5p_)
        {
            // Set up file access property list with parallel I/O access
            herr_t err_id
                = H5Pset_fapl_mpio(access_plist, comm_active_, MPI_INFO_NULL);
            if (err_id < 0)
            {
                (*MPIdata::serr) << "H5Pset_fapl_mpio failed!!!" << std::endl;
            }
        }
        else
        {
            // The H5FD_CORE driver enables an application to work with a file
            // in memory, speeding reads and writes as no disk access is made.
            // File contents are stored only in memory until the file is closed.
            // The last parameter determines whether file contents are ever
            // written to disk.
            herr_t err_id = H5Pset_fapl_core(access_plist, 1024, 1);
            if (err_id < 0)
                (*MPIdata::serr)
                    << "HDFrestart(): H5Pset_fapl_core failed!!!" << std::endl;
        }
        /* create the file collectively */
        H5Pset_fclose_degree(access_plist, H5F_CLOSE_STRONG);
        file_id_ = H5Fcreate(
            filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access_plist);
        if (file_id_ < 0)
        {
            (*MPIdata::serr) << "HDFrestart(): H5Fcreate failed for file "
                             << filename_ << std::endl;
            ct.global_exit(2);
        }
        herr_t ret = H5Pclose(access_plist);
        if (ret < 0)
        {
            (*MPIdata::serr) << "HDFrestart(): H5Pclose failed!!!" << std::endl;
            ct.global_exit(2);
        }
        assert(ret != -1);
    }
    else
    {
        file_id_ = -1;
    }
    setupWorkSpace();

    verbosity_ = 0;

#ifdef GITHASH
#define xstr(g, x) #g #x
#define LOG(x) addReleaseNumber2File(xstr(git, x));
    LOG(GITHASH);
#endif

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Created HDF file with Data blocks: " << block_[0]
                         << " x " << block_[1] << " x " << block_[2]
                         << std::endl;

    create_file_tm_.stop();
}

// constructor reading data (existing file)
HDFrestart::HDFrestart(const std::string& filename, const pb::PEenv& pes,
    const short option_number)
    : pes_(pes), file_id_(-1)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    comm_data_      = mmpi.commSameSpin();

    open_existing_tm_.start();
    Control& ct = *(Control::instance());

    setOptions(option_number);
    filename_ = filename;
    if (!use_hdf5p_)
    {
        filename_.append("/Task");
        appendTaskNumberToFilename();
    }
    verbosity_ = 0;
    closed_    = false;

    setActivity();

    if (active_)
    {
        if (!fileExists(filename_.c_str()))
        {
            (*MPIdata::serr) << "ERROR MGmol: file " << filename_
                             << " does not exists!!!" << std::endl;
            ct.global_exit(2);
        }
        htri_t ishdf = H5Fis_hdf5(filename_.c_str());
        if (ishdf < 0)
        {
            (*MPIdata::serr)
                << "ERROR MGmol: H5Fis_hdf5() unsuccessful for file "
                << filename_ << std::endl;
            ct.global_exit(2);
        }
    }

    mmpi.barrier();

    count_[0] = count_[1] = count_[2] = 1;
    stride_[0] = stride_[1] = stride_[2] = 1;

    setCommActive();

    if (active_)
    {

        hid_t access_plist
            = H5Pcreate(H5P_FILE_ACCESS); // property list identifier
        if (use_hdf5p_)
        {
            // Set up file access property list with parallel I/O access
            herr_t err_id
                = H5Pset_fapl_mpio(access_plist, comm_active_, MPI_INFO_NULL);
            if (err_id < 0)
            {
                (*MPIdata::serr) << "H5Pset_fapl_mpio failed!!!" << std::endl;
                ct.global_exit(2);
            }
        }
        else
        {
            herr_t err_id = H5Pset_fapl_core(access_plist, 1024, 1);
            if (err_id < 0)
                (*MPIdata::serr)
                    << "HDFrestart(): H5Pset_fapl_core failed!!!" << std::endl;
            // else
            //    if( onpe0 )
            //        (*MPIdata::sout)<<"HDFrestart(): call
            //        H5Pset_fapl_core()"<<endl;
        }
        file_id_ = H5Fopen(filename_.c_str(), H5F_ACC_RDONLY, access_plist);
        H5Pclose(access_plist);

        if (file_id_ < 0)
        {
            (*MPIdata::serr) << "HDFrestart(): open " << filename_
                             << " failed!!!" << std::endl;
            ct.global_exit(2);
        }
        else if (onpe0 && ct.verbose > 0)
        {
            (*MPIdata::sout)
                << "HDFrestart(): open file " << filename_ << std::endl;
        }

        // open one dataset to get dataspace dimension
        hid_t dset_id = open_dset("Vtotal");
        if (dset_id < 0)
        {
            dimsf_[0] = 0;
            dimsf_[1] = 0;
            dimsf_[2] = 0;
        }
        else
        {
            hid_t dspace_id = H5Dget_space(dset_id);

            hsize_t maxdims[3] = { 0, 0, 0 };
            int status
                = H5Sget_simple_extent_dims(dspace_id, dimsf_, &maxdims[0]);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart(): H5Sget_simple_extent_dims failed!!!"
                    << std::endl;
                ct.global_exit(2);
            }

            // get global mesh size from file data
            if (!use_hdf5p_)
            {
                dimsf_[1] *= pes.n_mpi_task(1);
                dimsf_[2] *= pes.n_mpi_task(2);
            }
            if (!gather_data_x_) dimsf_[0] *= pes.n_mpi_task(0);

            // close dataset
            status = H5Dclose(dset_id);
            if (status < 0)
            {
                (*MPIdata::serr) << "H5Dclose failed!!!" << std::endl;
                return;
            }
        }
    }
    // Bcast size of data
    int dimsf[3] = { (int)dimsf_[0], (int)dimsf_[1], (int)dimsf_[2] };
    mmpi.bcast(&dimsf[0], 3);

    dimsf_[0] = dimsf[0];
    dimsf_[1] = dimsf[1];
    dimsf_[2] = dimsf[2];

    setupBlocks();

    setupWorkSpace();

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart() --- Data blocks: " << block_[0]
                         << " x " << block_[1] << " x " << block_[2]
                         << std::endl;
    assert(bsize_ >= 0);

    open_existing_tm_.stop();
}

int writeListCentersAndRadii(
    hid_t dset_id, const unsigned natt, const std::vector<double>& attr_data)
{
    assert(dset_id > -1);
    assert(attr_data.size() < 10000);

    if (natt <= 0) return 0;

    hsize_t attsize = attr_data.size() / natt;

    // Create the data space for the attribute.
    hsize_t dims[2] = { natt, attsize };
    // assert( dims[0]>0 );

    std::string attname("List of centers and radii");
    //  Open a dataset attribute.
    hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);
    hid_t attribute_id = H5Acreate2(dset_id, attname.c_str(), H5T_NATIVE_DOUBLE,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attribute_id < 0)
    {
        (*MPIdata::serr) << "Attribute: " << attname
                         << " --- H5Acreate failed!!!" << std::endl;
        return -1;
    }

    // if( onpe0 && ct.verbose>2 )
    //    (*MPIdata::sout)<<"Write attribute "<<attname<<endl;
    herr_t status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attr_data[0]);
    if (status < 0)
    {
        (*MPIdata::serr) << "Attribute: " << attname
                         << " --- H5Awrite failed!!!" << std::endl;
        return -1;
    }
    status = H5Sclose(dataspace_id);
    if (status < 0)
    {
        (*MPIdata::serr) << "H5Sclose failed!!!" << std::endl;
        return -1;
    }
    status = H5Aclose(attribute_id);
    if (status < 0)
    {
        (*MPIdata::serr) << "Attribute: " << attname
                         << " --- H5Aclose failed!!!" << std::endl;
        return -1;
    }

    return natt;
}

int writeGids(hid_t dset_id, const std::vector<int>& gids)
{
    assert(dset_id >= -1);

    const int natt = gids.size();
    if (natt <= 0) return 0;

    std::string attname("List of gids");

    mgmol_tools::addAttribute2Dataset(dset_id, attname.c_str(), gids);

    return natt;
}

int readListCentersAndRadii(hid_t dset_id, std::vector<double>& attr_data)
{
    std::string attname("List of centers and radii");
    htri_t exists = H5Aexists(dset_id, attname.c_str());
    if (exists < 0)
    {
        (*MPIdata::serr) << "H5Aexists failed for " << attname << "!!!"
                         << std::endl;
        return exists;
    }
    if (exists == 0)
    {
        return 0;
    }

    // Open the data space for the attribute.
    //  Open a dataset attribute.
    hid_t attribute_id = H5Aopen_name(dset_id, attname.c_str());
    if (attribute_id < 0)
    {
        (*MPIdata::serr) << "H5Aopen failed for " << attname << "!!!"
                         << std::endl;
        return -1;
    }

    hid_t attdataspace = H5Aget_space(attribute_id);
    hsize_t dims[2]    = { 0, 0 };
    hsize_t maxdims[2] = { 0, 0 };
    H5Sget_simple_extent_dims(attdataspace, dims, maxdims);
    if (dims[1] != 4)
    {
        (*MPIdata::serr) << "HDFrestart::readAttributes() --- wrong data size"
                         << std::endl;
        return -1;
    }
    assert(dims[0] * dims[1] > 0);
    int natt = dims[0];
    attr_data.resize(dims[0] * dims[1]);
    // Read the attribute data.
#ifdef DEBUG
    (*MPIdata::sout) << "Read attribute " << attname << " of dimension "
                     << dims[0] << "*" << dims[1] << std::endl;
#endif
    herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &attr_data[0]);
    // check validity of data just read
    if (status < 0)
    {
        (*MPIdata::serr) << "H5Aread failed!!!" << std::endl;
        return -1;
    }
    for (unsigned int j = 0; j < dims[0] * dims[1]; j++)
        assert(fabs(attr_data[j]) < 1.e30);
    status = H5Sclose(attdataspace);
    if (status < 0)
    {
        (*MPIdata::serr) << "H5Sclose failed!!!" << std::endl;
        return -1;
    }
    status = H5Aclose(attribute_id);
    if (status < 0)
    {
        (*MPIdata::serr) << "H5Aclose failed!!!" << std::endl;
        return -1;
    }

    return natt;
}

int readGids(hid_t dset_id, std::vector<int>& gids)
{
    // static bool nogids=false;

    // if( nogids ){
    //    return 0;
    //}

    gids.clear();

    std::string attname("List of gids");
    htri_t exists = H5Aexists(dset_id, attname.c_str());
    if (exists < 0)
    {
        (*MPIdata::serr) << "H5Aexists() failed!!!" << std::endl;
        return exists;
    }
    if (exists == 0) // attribute does not exist
    {
        //(*MPIdata::sout)<<"No gids found in dataset "<<endl;
        // nogids=true; // -> don't look for gids for next dataset

        return 0;
    }

    // Open the data space for the attribute.
    //  Open a dataset attribute.
    hid_t attribute_id = H5Aopen_name(dset_id, attname.c_str());
    if (attribute_id < 0)
    {
        (*MPIdata::serr) << "H5Aopen failed for " << attname << "!!!"
                         << std::endl;
        return -1;
    }

    hid_t attdataspace = H5Aget_space(attribute_id);
    hsize_t dims;
    hsize_t maxdims;
    H5Sget_simple_extent_dims(attdataspace, &dims, &maxdims);
    assert(dims > 0);
    gids.resize(dims);
    // Read the attribute data.
#ifdef DEBUG
    (*MPIdata::sout) << "Read attribute " << attname << " of dimension " << dims
                     << std::endl;
#endif
    herr_t status = H5Aread(attribute_id, H5T_NATIVE_INT, &gids[0]);
    if (status < 0)
    {
        (*MPIdata::serr) << "H5Aread failed!!!" << std::endl;
        return -1;
    }
    status = H5Sclose(attdataspace);
    if (status < 0)
    {
        (*MPIdata::serr) << "H5Sclose failed!!!" << std::endl;
        return -1;
    }
    status = H5Aclose(attribute_id);
    if (status < 0)
    {
        (*MPIdata::serr) << "H5Aclose failed!!!" << std::endl;
        return -1;
    }

    return dims;
}

int HDFrestart::getLRCenters(std::multimap<std::string, Vector3D>& centers,
    const int n_max_centers, const std::string& name)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("HDFrestart::getLRCenters()...", (*MPIdata::sout));

    centers.clear();

    int color              = 0;
    int done               = 0;
    int ncenters           = -1;
    short attribute_length = 0;
    while (!done)
    {
        int dim = 0;
        std::vector<double> attr_data;

        std::string datasetname(getDatasetName(name, color));

        int err_id = dset_exists(datasetname);
        if (err_id == 0)
        { // dataset does not exists
            // try older version
            datasetname = getDatasetName_old(name, color);
            err_id      = dset_exists(datasetname);
        }

        if (err_id == 0)
        { // dataset does not exists
            if (onpe0 && ct.verbose > 0)
                (*MPIdata::sout)
                    << "Dataset " << datasetname
                    << " does not exists... Stop reading..." << std::endl;
            done = 1;
        }
        else
        {
            // Open dataset.
            hid_t dset_id = open_dset(datasetname);
            if (dset_id < 0)
            {
                if (verbosity_ > 1)
                    (*MPIdata::sout) << "H5Dopen failed for datasetname "
                                     << datasetname << std::endl;
                return -1;
            }
            else
            {
                if (verbosity_ > 1)
                    (*MPIdata::sout)
                        << "H5Dopen for dataset " << datasetname << std::endl;
            }

            int natt = readListCentersAndRadii(dset_id, attr_data);
            if (natt < 0) return -1;
            attribute_length = natt > 0 ? attr_data.size() / natt : 0;

            dim = natt * attribute_length;

            close_dset(dset_id);
        }

        if (!done)
        {
            if (dim > 0)
            {
                const int nd = dim / attribute_length;
                for (int i = 0; i < nd; i++)
                {
                    Vector3D center(attr_data[4 * i], attr_data[4 * i + 1],
                        attr_data[4 * i + 2]);
                    centers.insert(
                        std::pair<std::string, Vector3D>(datasetname, center));
                    if (verbosity_ > 2 && onpe0)
                    {
                        (*MPIdata::sout) << std::setprecision(8);
                        (*MPIdata::sout) << " HDFrestart::getLRCenters() --- "
                                            "Read function centers:";
                        (*MPIdata::sout) << center[0] << "\t";
                        (*MPIdata::sout) << center[1] << "\t";
                        (*MPIdata::sout) << center[2] << " of color "
                                         << datasetname << std::endl;
                    }
                }
            }

            if (n_max_centers == (int)centers.size())
            {
                done = 1;
            }
        } // !done

        ncenters = (int)centers.size();

        color++;
    }

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Read " << (int)centers.size()
                         << " LR centers in restart file" << std::endl;

    return ncenters;
}

// get distinct function centers and their multiplicities in file
int HDFrestart::getLRs(LocalizationRegions& lrs, const int max_nb_lrs,
    const std::string& name, const bool add)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("Reading localization regions...", (*MPIdata::sout));

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::getLRs() --- Read LRs from dataset "
                         << name << std::endl;

    if (!add) lrs.clear();

    const bool read_radius = (ct.lr_updates_type > 0 || ct.cut_radius < 0.);
    if (!read_radius && onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart: Use uniform truncation radius from "
                            "input file, rcut="
                         << ct.cut_radius << std::endl;

    int nlrs               = 0;
    int color              = 0;
    int done               = 0;
    short attribute_length = 0;
    while (!done)
    {
        int dim = 0;
        std::vector<double> attr_data;
        std::vector<int> gids;

        std::string datasetname(getDatasetName(name, color));

        int err_id = dset_exists(datasetname);
        if (err_id == 0)
        { // dataset does not exists
            // try older version
            datasetname = getDatasetName_old(name, color);
            err_id      = dset_exists(datasetname);
        }
        if (err_id == 0)
        { // dataset does not exists
            if (onpe0 && ct.verbose > 0)
                (*MPIdata::sout)
                    << "HDFrestart::getLRs(), Dataset " << datasetname
                    << " does not exists..." << std::endl;
            done = 1;
        }
        else
        {
            // Open dataset.
            hid_t dset_id = open_dset(datasetname);

            int natt = readListCentersAndRadii(dset_id, attr_data);
            if (natt < 0) return -1;
            attribute_length = natt > 0 ? attr_data.size() / natt : 0;

            dim = natt * attribute_length;

            nlrs += natt;

            readGids(dset_id, gids);
            assert(gids.size() == static_cast<unsigned int>(natt));

            herr_t status = close_dset(dset_id);
            if (status < 0)
            {
                return status;
            }
        }

        if (!done && dim > 0)
        {
            assert(attribute_length > 0);
            const int nd = dim / attribute_length;
            for (int i = 0; i < nd; i++)
            {
                const Vector3D center(attr_data[attribute_length * i],
                    attr_data[attribute_length * i + 1],
                    attr_data[attribute_length * i + 2]);

                // use truncation radius from input file
                // if no adaptation and ct.cut_radius>0.
                double rl = ct.cut_radius;
                if (read_radius) rl = attr_data[attribute_length * i + 3];

                assert(!gids.empty());
                assert(gids.size() > static_cast<unsigned int>(i));
                lrs.push_back_local(center, rl, gids[i]);

#ifdef DEBUG
                (*MPIdata::sout) << setprecision(16);
                (*MPIdata::sout) << "Read LR:";
                (*MPIdata::sout) << center[0] << "\t";
                (*MPIdata::sout) << center[1] << "\t";
                (*MPIdata::sout) << center[2] << "\t, R=";
                (*MPIdata::sout) << rl << std::endl;
#endif

                assert(rl > 0.2);
            }
        }

        if (lrs.hasNcentersReachedNumber(max_nb_lrs)) done = 1;

        color++;

    } // while !done

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::getLRs() --- Read " << nlrs
                         << " LRs from restart file on PE0" << std::endl;

    if (ct.verbose > 0)
        printWithTimeStamp(
            "Done reading localization regions...", (*MPIdata::sout));

    mmpi.allreduce(&nlrs, 1, MPI_MAX);

    return nlrs;
}

/////////////////////////////////////////////////////////////////////////////

int HDFrestart::read_1func_hdf5(double* vv, const std::string& datasetname)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::read_1func_hdf5(). Try to read data "
                         << datasetname << "..." << std::endl;

    assert(block_[0] * block_[1] * block_[2] > 0);

    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        hid_t dset_id  = open_dset(datasetname);
        if (dset_id < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_1func_hdf5() --- open_dset() failed!!!"
                << std::endl;
            return -1;
        }

        // memory dataspace identifier
        hid_t memspace = createMemspace();
        if (memspace < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_1func_hdf5() --- createMemspace failed!!!"
                << std::endl;
            return -1;
        }

        hid_t filespace = H5Dget_space(dset_id);
        if (filespace < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_1func_hdf5() --- H5Dget_space failed!!!"
                << std::endl;
            return -1;
        }

        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            herr_t status = H5Sselect_hyperslab(
                filespace, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                (*MPIdata::serr) << "HDFrestart::read_1func_hdf5() --- "
                                    "H5Sselect_hyperslab failed!!!"
                                 << std::endl;
                return -1;
            }
            // Create property list for collective dataset read.
            if (pes_.n_mpi_tasks() > 1)
            {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
            }
            else
            {
                plist_id = H5P_DEFAULT;
            }
        }
        herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
            plist_id, work_space_double_);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_1func_hdf5() --- H5Dread failed!!!"
                << std::endl;
            return -1;
        }

        // Close/release resources.
        if (use_hdf5p_)
            if (pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);

        status = H5Dclose(dset_id);
        if (status < 0)
        {
            (*MPIdata::serr) << "H5Dclose failed!!!" << std::endl;
            return -1;
        }

        if (use_hdf5p_)
        {
            status = H5Sclose(filespace);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::read_1func_hdf5() --- H5Sclose failed!!!"
                    << std::endl;
                return -1;
            }

            status = H5Sclose(memspace);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::read_1func_hdf5() --- H5Sclose failed!!!"
                    << std::endl;
                return -1;
            }
        }

        memcpy(vv, work_space_double_, bsize_ * sizeof(double));

        if (onpe0 && ct.verbose > 0)
            (*MPIdata::sout)
                << "Dataset " << datasetname << " read" << std::endl;
    }

    // send data to inactive PEs
    if (gather_data_x_)
    {
        int j      = pes_.my_mpi(1);
        int k      = pes_.my_mpi(2);
        int source = pes_.xyz2task(0, j, k);
        MPI_Status mpistatus;
        for (int i = 1; i < pes_.n_mpi_task(0); i++)
        {
            int dest = pes_.xyz2task(i, j, k);
            int tag  = i;
            if (pes_.my_mpi(0) == 0)
            {
                assert(active_);
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send "<<bsize_<<"
                // data to "<<dest<<endl;
                MPI_Send(work_space_double_ + i * bsize_, bsize_, MPI_DOUBLE,
                    dest, tag, comm_data_);
            }
            else if (pes_.my_mpi(0) == i)
            {
                assert(!active_);
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive
                //"<<bsize_<<" data from "<<source<<endl;
                MPI_Recv(vv, bsize_, MPI_DOUBLE, source, tag, comm_data_,
                    &mpistatus);
            }
        }
    }

    return 0;
}

int HDFrestart::read_1func_hdf5(float* vv, const std::string& datasetname)
{
    if (onpe0)
        (*MPIdata::sout) << "HDFrestart::read_1func_hdf5(). Try to read data "
                         << datasetname << "..." << std::endl;

    assert(block_[0] * block_[1] * block_[2] > 0);

    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        hid_t dset_id  = open_dset(datasetname);
        if (dset_id < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_1func_hdf5() --- open_dset() failed!!!"
                << std::endl;
            return -1;
        }

        // memory dataspace identifier
        hid_t memspace = createMemspace();
        if (memspace < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_1func_hdf5() --- createMemspace failed!!!"
                << std::endl;
            return -1;
        }

        hid_t filespace = H5Dget_space(dset_id);
        if (filespace < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_1func_hdf5() --- H5Dget_space failed!!!"
                << std::endl;
            return -1;
        }

        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            herr_t status = H5Sselect_hyperslab(
                filespace, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                (*MPIdata::serr) << "HDFrestart::read_1func_hdf5() --- "
                                    "H5Sselect_hyperslab failed!!!"
                                 << std::endl;
                return -1;
            }
            // Create property list for collective dataset read.
            if (pes_.n_mpi_tasks() > 1)
            {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
            }
            else
            {
                plist_id = H5P_DEFAULT;
            }
        }
        herr_t status = H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
            plist_id, work_space_float_);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_1func_hdf5() --- H5Dread failed!!!"
                << std::endl;
            return -1;
        }

        // Close/release resources.
        if (use_hdf5p_)
            if (pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);

        status = H5Dclose(dset_id);
        if (status < 0)
        {
            (*MPIdata::serr) << "H5Dclose failed!!!" << std::endl;
            return -1;
        }

        if (use_hdf5p_)
        {
            status = H5Sclose(filespace);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::read_1func_hdf5() --- H5Sclose failed!!!"
                    << std::endl;
                return -1;
            }

            status = H5Sclose(memspace);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::read_1func_hdf5() --- H5Sclose failed!!!"
                    << std::endl;
                return -1;
            }
        }

        memcpy(vv, work_space_float_, bsize_ * sizeof(float));

        if (onpe0)
            (*MPIdata::sout)
                << "Dataset " << datasetname << " read" << std::endl;
    }

    // send data to inactive PEs
    if (gather_data_x_)
    {
        int j      = pes_.my_mpi(1);
        int k      = pes_.my_mpi(2);
        int source = pes_.xyz2task(0, j, k);
        MPI_Status mpistatus;
        for (int i = 1; i < pes_.n_mpi_task(0); i++)
        {
            int dest = pes_.xyz2task(i, j, k);
            int tag  = i;
            if (pes_.my_mpi(0) == 0)
            {
                assert(active_);
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send "<<bsize_<<"
                // data to "<<dest<<endl;
                MPI_Send(work_space_float_ + i * bsize_, bsize_, MPI_FLOAT,
                    dest, tag, comm_data_);
            }
            else if (pes_.my_mpi(0) == i)
            {
                assert(!active_);
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive
                //"<<bsize_<<" data from "<<source<<endl;
                MPI_Recv(
                    vv, bsize_, MPI_FLOAT, source, tag, comm_data_, &mpistatus);
            }
        }
    }

    return 0;
}

int HDFrestart::write_1func_hdf5(
    double* vv, const std::string& datasetname, double* ll, double* cell_origin)
{
    assert(ll != nullptr);
    assert(cell_origin != nullptr);

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::write_1func_hdf5(). Try to write data "
                         << datasetname << "..." << std::endl;

    hid_t filespace = -1;
    hid_t dset_id   = -1;
    hid_t memspace  = -1;

    if (active_)
    {
        // filespace identifier
        filespace = createFilespace();

        // memory dataspace identifier
        memspace = createMemspace();

        hid_t plist_id = createPlist();

        dset_id = H5Dcreate2(file_id_, datasetname.c_str(), H5T_NATIVE_DOUBLE,
            filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        if (dset_id < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::write_1func_hdf5: H5Dcreate failed!!!"
                << std::endl;
            return -1;
        }

        releasePlist(plist_id);

        // Write the attribute "Lattice parameters"
        if (ll != nullptr)
        {
            // Create the data space for the attribute "Lattice parameters".
            std::vector<double> attr_data(3);
            attr_data[0] = ll[0];
            attr_data[1] = ll[1];
            attr_data[2] = ll[2];

            std::string attname("Lattice parameters");
            mgmol_tools::addAttribute2Dataset(
                dset_id, attname.c_str(), attr_data);
        }

        if (cell_origin != nullptr)
        {
            std::vector<double> attr_data(3);
            attr_data[0] = cell_origin[0];
            attr_data[1] = cell_origin[1];
            attr_data[2] = cell_origin[2];

            std::string attname("Cell origin");
            mgmol_tools::addAttribute2Dataset(
                dset_id, attname.c_str(), attr_data);
        }
    }

    int ierr = writeData(vv, filespace, memspace, dset_id, 2);
    if (ierr < 0) return ierr;

    if (active_)
    {
        herr_t status = H5Dclose(dset_id);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::write_1func_hdf5: H5Dclose failed!!!"
                << std::endl;
            return -1;
        }
        // if( use_hdf5p_ )
        {
            status = H5Sclose(filespace);
            if (status < 0)
            {
                (*MPIdata::serr) << "HDFrestart::write_1func_hdf5: H5Sclose "
                                    "failed for filespace!!!"
                                 << std::endl;
                return -1;
            }

            status = H5Sclose(memspace);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "H5Sclose failed for memspace!!!" << std::endl;
                return -1;
            }
        }
    }

    return 0;
}

int HDFrestart::write_1func_hdf5(
    float* vv, const std::string& datasetname, double* ll, double* cell_origin)
{
    assert(ll != nullptr);
    assert(cell_origin != nullptr);

    if (onpe0)
        (*MPIdata::sout) << "HDFrestart::write_1func_hdf5(). Try to write data "
                         << datasetname << "..." << std::endl;

    hid_t filespace = -1;
    hid_t dset_id   = -1;
    hid_t memspace  = -1;

    if (active_)
    {
        // filespace identifier
        filespace = createFilespace();

        // memory dataspace identifier
        memspace = createMemspace();

        hid_t plist_id = createPlist();

        dset_id = H5Dcreate2(file_id_, datasetname.c_str(), H5T_NATIVE_FLOAT,
            filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        if (dset_id < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::write_1func_hdf5: H5Dcreate failed!!!"
                << std::endl;
            return -1;
        }

        releasePlist(plist_id);

        // Write the attribute "Lattice parameters"
        if (ll != nullptr)
        {
            // Create the data space for the attribute "Lattice parameters".
            std::vector<double> attr_data(3);
            attr_data[0] = ll[0];
            attr_data[1] = ll[1];
            attr_data[2] = ll[2];

            std::string attname("Lattice parameters");
            mgmol_tools::addAttribute2Dataset(
                dset_id, attname.c_str(), attr_data);
        }

        if (cell_origin != nullptr)
        {
            std::vector<double> attr_data(3);
            attr_data[0] = cell_origin[0];
            attr_data[1] = cell_origin[1];
            attr_data[2] = cell_origin[2];

            std::string attname("Cell origin");
            mgmol_tools::addAttribute2Dataset(
                dset_id, attname.c_str(), attr_data);
        }
    }

    int ierr = writeData(vv, filespace, memspace, dset_id, 1);
    if (ierr < 0) return ierr;

    if (active_)
    {
        herr_t status = H5Dclose(dset_id);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::write_1func_hdf5: H5Dclose failed!!!"
                << std::endl;
            return -1;
        }
        // if( use_hdf5p_ )
        {
            status = H5Sclose(filespace);
            if (status < 0)
            {
                (*MPIdata::serr) << "HDFrestart::write_1func_hdf5: H5Sclose "
                                    "failed for filespace!!!"
                                 << std::endl;
                return -1;
            }

            status = H5Sclose(memspace);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "H5Sclose failed for memspace!!!" << std::endl;
                return -1;
            }
        }
    }

    return 0;
}

int HDFrestart::readData(
    double* data, hid_t memspace, hid_t dset_id, const short precision)
{
    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        herr_t status;
        hid_t filespace = H5Dget_space(dset_id);
        if (filespace < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::readData: H5Dget_space failed !!!" << std::endl;
            return -1;
        }
        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            status = H5Sselect_hyperslab(
                filespace, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::readData: H5Sselect_hyperslab failed!!!"
                    << std::endl;
                return -1;
            }

            // Create property list for collective dataset read.
            if (pes_.n_mpi_tasks() > 1)
            {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
            }
        }
        if (precision == 1)
            status = H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                plist_id, work_space_float_);
        else
            status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                plist_id, work_space_double_);

        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::readData: H5Dread  failed !!!" << std::endl;
            return -1;
        }

        // Close/release resources.
        if (use_hdf5p_ && pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);

        H5Sclose(filespace);
    }
    if (gather_data_x_)
    {
        int j      = pes_.my_mpi(1);
        int k      = pes_.my_mpi(2);
        int source = pes_.xyz2task(0, j, k);
        MPI_Status mpistatus;

        if (precision == 2)
            for (int i = 1; i < pes_.n_mpi_task(0); i++)
            {
                int dest = pes_.xyz2task(i, j, k);
                int tag  = i;
                if (pes_.my_mpi(0) == 0)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send data to
                    //"<<dest<<endl;
                    MPI_Send(work_space_double_ + i * bsize_, bsize_,
                        MPI_DOUBLE, dest, tag, comm_data_);
                }
                else if (pes_.my_mpi(0) == i)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data
                    // from "<<ipe<<endl;
                    MPI_Recv(work_space_double_, bsize_, MPI_DOUBLE, source,
                        tag, comm_data_, &mpistatus);
                }
            }
        else
            for (int i = 1; i < pes_.n_mpi_task(0); i++)
            {
                int dest = pes_.xyz2task(i, j, k);
                int tag  = i;
                if (pes_.my_mpi(0) == 0)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send data to
                    //"<<dest<<endl;
                    MPI_Send(work_space_float_ + i * bsize_, bsize_, MPI_FLOAT,
                        dest, tag, comm_data_);
                }
                else if (pes_.my_mpi(0) == i)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data
                    // from "<<ipe<<endl;
                    MPI_Recv(work_space_float_, bsize_, MPI_FLOAT, source, tag,
                        comm_data_, &mpistatus);
                }
            }
    }
    if (precision == 1)
    {
        for (int i = 0; i < bsize_; i++)
            data[i] = (double)work_space_float_[i];
    }
    else
    {
        memcpy(data, work_space_double_, bsize_ * sizeof(double));
    }
    return 0;
}

// Note precision is ignored. It is left as a argument for convenience
int HDFrestart::readData(
    float* data, hid_t memspace, hid_t dset_id, const short precision)
{
    (void)precision;

    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        herr_t status;
        hid_t filespace = H5Dget_space(dset_id);
        if (filespace < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::readData: H5Dget_space failed !!!" << std::endl;
            return -1;
        }
        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            status = H5Sselect_hyperslab(
                filespace, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::readData: H5Sselect_hyperslab failed!!!"
                    << std::endl;
                return -1;
            }

            // Create property list for collective dataset read.
            if (pes_.n_mpi_tasks() > 1)
            {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
            }
        }
        status = H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
            plist_id, work_space_float_);

        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::readData: H5Dread  failed !!!" << std::endl;
            return -1;
        }

        // Close/release resources.
        if (use_hdf5p_ && pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);

        H5Sclose(filespace);
    }
    if (gather_data_x_)
    {
        int j      = pes_.my_mpi(1);
        int k      = pes_.my_mpi(2);
        int source = pes_.xyz2task(0, j, k);
        MPI_Status mpistatus;

        for (int i = 1; i < pes_.n_mpi_task(0); i++)
        {
            int dest = pes_.xyz2task(i, j, k);
            int tag  = i;
            if (pes_.my_mpi(0) == 0)
            {
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send data to
                //"<<dest<<endl;
                MPI_Send(work_space_float_ + i * bsize_, bsize_, MPI_FLOAT,
                    dest, tag, comm_data_);
            }
            else if (pes_.my_mpi(0) == i)
            {
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data from
                //"<<ipe<<endl;
                MPI_Recv(work_space_float_, bsize_, MPI_FLOAT, source, tag,
                    comm_data_, &mpistatus);
            }
        }
    }
    memcpy(data, work_space_float_, bsize_ * sizeof(float));
    return 0;
}

int HDFrestart::writeData(double* data, hid_t space_id, hid_t memspace,
    hid_t dset_id, const short precision)
{
    if (precision == 1)
    {
        for (int i = 0; i < bsize_; i++)
            work_space_float_[i] = (float)data[i];
    }
    else
    {
        assert(work_space_double_ != nullptr);
        memcpy(work_space_double_, data, bsize_ * sizeof(double));
    }
    // gather data on active PEs
    if (gather_data_x_)
    {
        int j    = pes_.my_mpi(1);
        int k    = pes_.my_mpi(2);
        int dest = pes_.xyz2task(0, j, k);
        MPI_Status mpistatus;

        // blocking send/recv
        if (precision == 2)
            for (int i = 1; i < pes_.n_mpi_task(0); i++)
            {
                int ipe = pes_.xyz2task(i, j, k);
                if (pes_.my_mpi(0) == i)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send data to
                    //"<<dest<<endl;
                    MPI_Send(work_space_double_, bsize_, MPI_DOUBLE, dest, i,
                        comm_data_);
                }
                else if (pes_.my_mpi(0) == 0)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data
                    // from "<<ipe<<endl;
                    MPI_Recv(work_space_double_ + i * bsize_, bsize_,
                        MPI_DOUBLE, ipe, i, comm_data_, &mpistatus);
                }
            }
        else
            for (int i = 1; i < pes_.n_mpi_task(0); i++)
            {
                int ipe = pes_.xyz2task(i, j, k);
                if (pes_.my_mpi(0) == i)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send data to
                    //"<<dest<<endl;
                    MPI_Send(work_space_float_, bsize_, MPI_FLOAT, dest, i,
                        comm_data_);
                }
                else if (pes_.my_mpi(0) == 0)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data
                    // from "<<ipe<<endl;
                    MPI_Recv(work_space_float_ + i * bsize_, bsize_, MPI_FLOAT,
                        ipe, i, comm_data_, &mpistatus);
                }
            }
    }

    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        herr_t status;
        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            status = H5Sselect_hyperslab(
                space_id, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::write_data: H5Sselect_hyperslab failed!!!"
                    << std::endl;
                (*MPIdata::serr)
                    << "HDFrestart::write_data: space_id=" << space_id
                    << std::endl;
                return -1;
            }

            // Create property list for collective dataset write.
            if (pes_.n_mpi_tasks() > 1)
            {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
            }
        }

        assert(dset_id >= 0);
        if (precision == 1)
            status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, space_id,
                plist_id, work_space_float_);
        else
            status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id,
                plist_id, work_space_double_);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::write_data: H5Dwrite failed!!!" << std::endl;
            return -1;
        }

        // Close/release resources.
        if (use_hdf5p_ && pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);
    }

    return 0;
}

// precision is ignored.
int HDFrestart::writeData(float* data, hid_t space_id, hid_t memspace,
    hid_t dset_id, const short precision)
{
    (void)precision;

    memcpy(work_space_float_, data, bsize_ * sizeof(float));

    // gather data on active PEs
    if (gather_data_x_)
    {
        int j    = pes_.my_mpi(1);
        int k    = pes_.my_mpi(2);
        int dest = pes_.xyz2task(0, j, k);
        MPI_Status mpistatus;

        // blocking send/recv
        for (int i = 1; i < pes_.n_mpi_task(0); i++)
        {
            int ipe = pes_.xyz2task(i, j, k);
            if (pes_.my_mpi(0) == i)
            {
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send data to
                //"<<dest<<endl;
                MPI_Send(
                    work_space_float_, bsize_, MPI_FLOAT, dest, i, comm_data_);
            }
            else if (pes_.my_mpi(0) == 0)
            {
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data from
                //"<<ipe<<endl;
                MPI_Recv(work_space_float_ + i * bsize_, bsize_, MPI_FLOAT, ipe,
                    i, comm_data_, &mpistatus);
            }
        }
    }

    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        herr_t status;
        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            status = H5Sselect_hyperslab(
                space_id, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "HDFrestart::write_data: H5Sselect_hyperslab failed!!!"
                    << std::endl;
                (*MPIdata::serr)
                    << "HDFrestart::write_data: space_id=" << space_id
                    << std::endl;
                return -1;
            }

            // Create property list for collective dataset write.
            if (pes_.n_mpi_tasks() > 1)
            {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
            }
        }

        assert(dset_id >= 0);
        status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, space_id,
            plist_id, work_space_float_);

        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::write_data: H5Dwrite failed!!!" << std::endl;
            return -1;
        }

        // Close/release resources.
        if (use_hdf5p_ && pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);
    }

    return 0;
}

int HDFrestart::read_att(const hid_t dset_id, const std::string& attname,
    std::vector<double>& attr_data)
{
    assert(attr_data.size() > 0);

    int dim = (int)attr_data.size();
    if (active_)
    {
        hid_t attribute_id = H5Aopen_name(dset_id, attname.c_str());
        if (attribute_id < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_att() --- H5Aopen_name failed for "
                << attname << std::endl;
            return -1;
        }
        herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &attr_data[0]);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_att(), H5Aread failed!!!" << std::endl;
            return -1;
        }
        status = H5Aclose(attribute_id);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::read_att(), H5Aclose failed!!!" << std::endl;
            return -1;
        }
    }

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(&attr_data[0], dim);

    return 0;
}

void HDFrestart::setCommActive()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int mytask      = 0;
    MPI_Comm_rank(mmpi.commSameSpin(), &mytask);
    if (gather_data_x_)
    {
        int color = (int)active_;
        int mpirc
            = MPI_Comm_split(mmpi.commSameSpin(), color, mytask, &comm_active_);
        if (mpirc != MPI_SUCCESS)
        {
            (*MPIdata::serr) << "ERROR in MPI_Comm_split!!!" << std::endl;
            comm_active_ = MPI_COMM_NULL;
            return;
        }
    }
    else
    {
        comm_active_ = MPI_COMM_NULL;
    }
}

void HDFrestart::setOptions(const short option_number)
{
    switch (option_number)
    {
        case 0: // 1 file/task
            use_hdf5p_     = false;
            gather_data_x_ = false;
            break;
        case 1: // 1 file for all tasks
            use_hdf5p_     = true;
            gather_data_x_ = true;
            break;
        case 2: // 1 file per column in x-direction
            use_hdf5p_     = false;
            gather_data_x_ = true;
            break;
        default:
            (*MPIdata::serr) << "HDFrestart::setOptions() --- type "
                             << option_number << " undefined!!!" << std::endl;
            return;
    }
}

void HDFrestart::setupWorkSpace()
{
    // if( active_ )
    {
        const int n        = block_[0] * block_[1] * block_[2];
        work_space_double_ = new double[n];
        memset(work_space_double_, 0, n * sizeof(double));

        work_space_float_ = new float[n];
        memset(work_space_float_, 0, n * sizeof(float));
    }
}

void HDFrestart::closeWorkSpace()
{
    // if( active_ )
    {
        delete[] work_space_double_;
        delete[] work_space_float_;
    }
}

void HDFrestart::printTimers(std::ostream& os)
{
    open_existing_tm_.print(os);
    create_file_tm_.print(os);
    close_file_tm_.print(os);
}
/*
int HDFrestart::writeRandomState(unsigned short int rand_state[3])
{
    if( active_ ){
        // Create the data space for new datasets
        hsize_t dims[1]={3};
        hid_t    dataspace_id = H5Screate_simple(1, dims, NULL);
        if( dataspace_id<0 ){
            (*MPIdata::serr)<<"HDFrestart::writeRandomState(): H5Screate_simple
failed!!!"<<endl; return -1;
        }
        // Open dataset
        hid_t    dataset_id = H5Dcreate2(file_id_, "/Random_state",
                               H5T_NATIVE_INT,
                               dataspace_id, H5P_DEFAULT, H5P_DEFAULT,
H5P_DEFAULT); if( dataset_id<0 ){
            (*MPIdata::serr)<<"HDFrestart::writeRandomState()::H5Dcreate
/Random_state failed!!!"<<endl; return -1;
        }
        if ( onpe0 ){
            (*MPIdata::sout)<<"HDFrestart::writeRandomState(): State for random
numbers generator: "
                <<rand_state[0]<<","
                <<rand_state[1]<<","
                <<rand_state[2]<<endl;
        }
        int randst[3]={(int)rand_state[0],
                       (int)rand_state[1],
                       (int)rand_state[2]};
        if( onpe0 )
        {
            herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL,
H5S_ALL, H5P_DEFAULT, &randst[0]); if( status<0 ){
                (*MPIdata::serr)<<"HDFrestart::writeRandomState(): H5Dwrite
randst failed!!!"<<endl; return -1; }else{
                (*MPIdata::sout)<<"Random randst written into "
                    <<filename_<<endl;
            }
        }

        herr_t status = H5Dclose(dataset_id);
        if( status<0 ){
            (*MPIdata::serr)<<"HDFrestart::writeRandomState(): H5Dclose
failed!!!"<<endl; return -1;
        }
        status = H5Sclose(dataspace_id);
        if( status<0 ){
            (*MPIdata::serr)<<"HDFrestart::writeRandomState(): H5Sclose
failed!!!"<<endl; return -1;
        }
    }
    return 0;
}
*/

/*
int HDFrestart::readRandomState(unsigned short* rand_state)
{
    std::string function_name("HDFrestart::readRandomState()");

    int randst[3]={0,0,0}; //{(int)_init_rand_state[0],
(int)_init_rand_state[1], (int)_init_rand_state[2]}; // default values

    if( onpe0 )
    {
        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, "/Random_state",H5P_DEFAULT);
        if( dataset_id<0 ){
            if( onpe0 ){
                (*MPIdata::sout)<<function_name<<" --- H5Dopen() failed for
/Random_state "<<endl;
            }
        }else{
            herr_t status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL,
H5S_ALL, H5P_DEFAULT, &randst[0]); if( status<0 ){
                (*MPIdata::serr)<<function_name<<" --- H5Dread() failed for
/Random_state!!!"<<endl; return -1;
            }
            else{
                (*MPIdata::sout)<<"HDFrestart::readRandomState(): State for
random numbers generator: "
                    <<randst[0]<<","
                    <<randst[1]<<","
                    <<randst[2]<<endl;
            }

            // close dataset
            status = H5Dclose(dataset_id);
            if( status<0 ){
                (*MPIdata::serr)<<function_name<<" --- H5Dclose
failed!!!"<<endl; return -1;
            }
        }
    }
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(&randst[0], 3);

    for(short i=0;i<3;i++)rand_state[i]=(unsigned short)randst[i];

    return 0;
}
*/
int HDFrestart::readAtomicNumbers(std::vector<int>& data)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
    {
        (*MPIdata::sout) << "HDFrestart::readAtomicNumbers()..." << std::endl;
    }

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, "/Atomic_numbers", H5P_DEFAULT);
        if (!exists) return 0;

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, "/Atomic_numbers", H5P_DEFAULT);
        if (dataset_id < 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart::readAtomicNumbers() --- H5Dopen2 failed!!!"
                    << std::endl;
            return -1;
        }

        int dim = (int)(H5Dget_storage_size(dataset_id) / sizeof(int));
        if (dim == 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart::readAtomicNumbers() --- No numbers!!!"
                    << std::endl;
            return -1;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart::readAtomicNumbers() --- H5Dread failed!!!"
                << std::endl;
            return -1;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart::readAtomicNumbers() --- H5Dclose failed!!!"
                << std::endl;
            return -1;
        }
    } // if active_
    // send data to inactive PEs
    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

// return -2 means failure
// return -1 means dataset does not exists, and could be from older MGmol
// version
int HDFrestart::readAtomicIDs(std::vector<int>& data)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::readAtomicIDs()..." << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);
        htri_t exists = H5Lexists(file_id_, "/Atomic_IDs", H5P_DEFAULT);
        if (!exists) return -1;

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, "/Atomic_IDs", H5P_DEFAULT);
        if (dataset_id < 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart::readAtomicIDs() --- H5Dopen2 failed!!!"
                    << std::endl;
            return -2;
        }

        int dim = (int)(H5Dget_storage_size(dataset_id) / sizeof(int));
        if (dim == 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart::readAtomicIDs() --- No IDs!!!" << std::endl;
            return -2;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart::readAtomicIDs() --- H5Dread failed!!!"
                << std::endl;
            return -2;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart::readAtomicIDs() --- H5Dclose failed!!!"
                << std::endl;
            return -2;
        }
    }
    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

// return -2 means failure
// return -1 means dataset does not exists, and could be from older MGmol
// version
int HDFrestart::readAtomicNLprojIDs(std::vector<int>& data)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::readAtomicNLprojIDs()..." << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, "/AtomicNLproj_IDs", H5P_DEFAULT);
        if (!exists) return -1;

        hid_t dataset_id = H5Dopen2(file_id_, "/AtomicNLproj_IDs", H5P_DEFAULT);
        if (dataset_id < 0)
        {
            if (onpe0)
                (*MPIdata::sout) << "HDFrestart::readAtomicNLprojIDs() --- "
                                    "H5Dopen2 failed!!!"
                                 << std::endl;
            return -2;
        }

        int dim = (int)(H5Dget_storage_size(dataset_id) / sizeof(int));
        if (dim == 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart::readAtomicNLprojIDs() --- No IDs!!!"
                    << std::endl;
            return -2;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart::readAtomicNLprojIDs() --- H5Dread failed!!!"
                << std::endl;
            return -2;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart::readAtomicNLprojIDs() --- H5Dclose failed!!!"
                << std::endl;
            return -2;
        }
    }

    if (gather_data_x_) gatherDataXdir(data);
    return 0;
}

int HDFrestart::readAtomicPositions(std::vector<double>& data)
{
    if (onpe0)
        (*MPIdata::sout) << "Read ionic positions from hdf5 file" << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, "/Ionic_positions", H5P_DEFAULT);
        if (!exists) return -1;

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, "/Ionic_positions", H5P_DEFAULT);
        if (dataset_id < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart:readAtomicPositions() --- H5Dopen2 failed!!!"
                << std::endl;
            return -2;
        }

        int dim = (int)H5Dget_storage_size(dataset_id) / sizeof(double);
        if (dim == 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart:readAtomicPositions() --- No positions!!!"
                    << std::endl;
            return -2;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart:readAtomicPositions() --- H5Dread failed!!!"
                << std::endl;
            return -2;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout) << "H5Dclose failed!!!" << std::endl;
            return -2;
        }
    }

    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

int HDFrestart::readOldCenterOnMesh(std::vector<double>& data, int i)
{
    if (onpe0)
        (*MPIdata::sout) << "Read old localization centers rounded to mesh "
                            "points from hdf5 file"
                         << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);

        std::stringstream datasetstream;
        datasetstream << "OldCenterOnMesh_" << i;

        std::string datasetname = datasetstream.str();

        htri_t exists = H5Lexists(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (!exists) return -1;

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (dataset_id < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart:readOldCenterOnMesh() --- H5Dopen2 failed!!!"
                << std::endl;
            return -2;
        }

        int dim = (int)H5Dget_storage_size(dataset_id) / sizeof(double);
        if (dim == 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart:readOldCenterOnMesh() --- No old centers!!!"
                    << std::endl;
            return -2;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart:readOldCenterOnMesh() --- H5Dread failed!!!"
                << std::endl;
            return -2;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout) << "H5Dclose failed!!!" << std::endl;
            return -2;
        }
    }

    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

int HDFrestart::readOldCenter(std::vector<double>& data, int i)
{
    if (onpe0)
        (*MPIdata::sout) << "Read old localization centers from hdf5 file"
                         << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);

        std::stringstream datasetstream;
        datasetstream << "OldCenter_" << i;

        std::string datasetname = datasetstream.str();

        htri_t exists = H5Lexists(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (!exists) return -1;

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (dataset_id < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart:readOldCenter() --- H5Dopen2 failed!!!"
                << std::endl;
            return -2;
        }

        int dim = (int)H5Dget_storage_size(dataset_id) / sizeof(double);
        if (dim == 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart:readOldCenter() --- No old centers!!!"
                    << std::endl;
            return -2;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart:readOldCenter() --- H5Dread failed!!!"
                << std::endl;
            return -2;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout) << "H5Dclose failed!!!" << std::endl;
            return -2;
        }
    }

    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

int HDFrestart::readGidsList(std::vector<int>& data)
{
    if (onpe0)
        (*MPIdata::sout) << "Read list of gids from hdf5 file" << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);

        std::string datasetname = "GidsList";

        htri_t exists = H5Lexists(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (!exists) return -1;

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (dataset_id < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart:readGidsList() --- H5Dopen2 failed!!!"
                << std::endl;
            return -2;
        }

        int dim = (int)H5Dget_storage_size(dataset_id) / sizeof(int);
        if (dim == 0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "HDFrestart:readGidsList() --- No GidsList!!!"
                    << std::endl;
            return -2;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart:readGidsList() --- H5Dread failed!!!"
                << std::endl;
            return -2;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout) << "H5Dclose failed!!!" << std::endl;
            return -2;
        }
    }

    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

int HDFrestart::readAtomicVelocities(std::vector<double>& data)
{
    if (onpe0)
        (*MPIdata::sout) << "Read atomic velocities from hdf5 file"
                         << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, "/Ionic_velocities", H5P_DEFAULT);
        if (!exists) return -1;

        // Open an existing dataset
        hid_t dataset_id = H5Dopen2(file_id_, "/Ionic_velocities", H5P_DEFAULT);
        if (dataset_id < 0)
        {
            if (onpe0)
                (*MPIdata::sout) << "HDFrestart::readAtomicVelocities(), "
                                    "H5Dopen failed->no velocities read"
                                 << std::endl;
            data.clear();
        }
        else
        {
            int dim = (int)H5Dget_storage_size(dataset_id) / sizeof(double);
            data.resize(dim);
        }

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            (*MPIdata::serr)
                << "HDFrestart::readAtomicVelocities() --- H5Dread failed!!!"
                << std::endl;
            return -2;
        }

        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout) << "H5Dclose failed!!!" << std::endl;
            return -2;
        }
    }

    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

int HDFrestart::readLockedAtomNames(std::vector<std::string>& data)
{
    if (onpe0)
        (*MPIdata::sout) << "HDFrestart::readLockedAtomNames()..." << std::endl;

    std::vector<char> buffer;
    short name_length = 7; // default, value used before February 2016

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, "/LockedAtomsNames", H5P_DEFAULT);
        if (!exists) return 0;

        hid_t dataset_id = H5Dopen2(file_id_, "/LockedAtomsNames", H5P_DEFAULT);
        if (dataset_id < 0)
        {
            if (onpe0)
                (*MPIdata::sout) << "HDFrestart::readLockedAtomNames(), "
                                    "H5Dopen failed->no locked atoms read"
                                 << std::endl;
            return -1;
        }

        std::string attname("String_Length");
        htri_t existsA = H5Aexists(dataset_id, attname.c_str());
        if (existsA)
        {
            hid_t attribute_id = H5Aopen_name(dataset_id, attname.c_str());
            herr_t status = H5Aread(attribute_id, H5T_NATIVE_INT, &name_length);
            // check validity of data just read
            if (status < 0)
            {
                (*MPIdata::serr) << "H5Aread failed!!!" << std::endl;
                return -1;
            }
        }

        int dim = (int)H5Dget_storage_size(dataset_id) / name_length;

        if (onpe0)
            (*MPIdata::sout)
                << "HDFrestart::readLockedAtomNames(), dataset size=" << dim
                << std::endl;

        if (dim == 0) return 0;

        buffer.resize(dim * name_length);

        // create type for std::strings of length IonData_MaxStrLength
        hid_t strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(strtype, name_length);
        herr_t status = H5Dread(
            dataset_id, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);
        if (status < 0)
        {
            (*MPIdata::sout)
                << "HDFrestart::readLockedAtomNames(), H5Dread failed!!!"
                << std::endl;
            return -1;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            (*MPIdata::sout) << "H5Dclose failed!!!" << std::endl;
            return -1;
        }
    }

    if (gather_data_x_) gatherDataXdir(buffer);

    data.clear();
    for (unsigned short i = 0; i < buffer.size(); i += name_length)
    {
        std::string t(&buffer[i], name_length);
        assert(t.size() > 0);

        stripLeadingAndTrailingBlanks(t);

        assert(t.size() > 0);
        data.push_back(t);
    }

    return 0;
}

int HDFrestart::readAtomicNames(std::vector<std::string>& data)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::readAtomicNames()..." << std::endl;

    std::vector<char> buffer;
    short name_length = 7; // default, value used before February 2016

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, "/Atomic_names", H5P_DEFAULT);
        if (exists)
        {
            // Open the dataset
            hid_t dataset_id = H5Dopen2(file_id_, "/Atomic_names", H5P_DEFAULT);
            if (dataset_id < 0)
            {
                (*MPIdata::sout)
                    << "HDFrestart::readAtomicNames() --- H5Dopen2 failed!!!"
                    << std::endl;
                return -1;
            }

            std::string attname("String_Length");
            htri_t existsA = H5Aexists(dataset_id, attname.c_str());
            if (existsA)
            {
                hid_t attribute_id = H5Aopen_name(dataset_id, attname.c_str());
                herr_t status
                    = H5Aread(attribute_id, H5T_NATIVE_SHORT, &name_length);
                // check validity of data just read
                if (status < 0)
                {
                    (*MPIdata::serr) << "H5Aread failed!!!" << std::endl;
                    return -1;
                }

                status = H5Aclose(attribute_id);
                if (status < 0)
                {
                    (*MPIdata::serr) << "H5Aclose failed!!!" << std::endl;
                    return -1;
                }
            }

            int dim = (int)H5Dget_storage_size(dataset_id) / name_length;
            if (dim == 0)
            {
                if (onpe0)
                    (*MPIdata::sout)
                        << "HDFrestart::readAtomicNames() --- No names!!!"
                        << std::endl;
                return -1;
            }

            buffer.resize(dim * name_length);

            // create type for std::strings of length IonData_MaxStrLength
            hid_t strtype = H5Tcopy(H5T_C_S1);
            H5Tset_size(strtype, name_length);
            herr_t status = H5Dread(
                dataset_id, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);
            if (status < 0)
            {
                (*MPIdata::sout)
                    << "HDFrestart::readAtomicNames() --- H5Dread failed!!!"
                    << std::endl;
                return -1;
            }
            status = H5Dclose(dataset_id);
            if (status < 0)
            {
                (*MPIdata::sout)
                    << "HDFrestart::readAtomicNames() --- H5Dclose failed!!!"
                    << std::endl;
                return -1;
            }
        }
    }

    if (gather_data_x_) gatherDataXdir(buffer);

    data.clear();
    for (unsigned short i = 0; i < buffer.size(); i += name_length)
    {
        std::string t(&buffer[i], name_length);
        assert(t.size() > 0);
        // cout<<"name="<<t<<endl;

        stripLeadingAndTrailingBlanks(t);
        // cout<<"stripped name="<<t<<endl;

        assert(t.size() > 0);
        data.push_back(t);
    }

    return 0;
}

int HDFrestart::readRestartRandomStates(std::vector<unsigned short>& data)
{
    if (onpe0)
        (*MPIdata::sout) << "Read atomic RandomStates from hdf5 file"
                         << std::endl;
    data.clear();

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, "/Ionic_RandomStates", H5P_DEFAULT);
        if (exists)
        {
            int dim = 0;
            // Open an existing dataset
            hid_t dataset_id
                = H5Dopen2(file_id_, "/Ionic_RandomStates", H5P_DEFAULT);
            if (dataset_id < 0)
            {
                if (onpe0)
                    (*MPIdata::sout) << "HDFrestart::readRestartRandomStates(),"
                                        " H5Dopen failed->no random states read"
                                     << std::endl;
                dim = 0;
            }
            else
            {
                dim = (int)H5Dget_storage_size(dataset_id)
                      / sizeof(unsigned short);
            }

            if (dim > 0)
            {
                data.resize(dim);
                herr_t status = H5Dread(dataset_id, H5T_NATIVE_USHORT, H5S_ALL,
                    H5S_ALL, H5P_DEFAULT, &data[0]);
                if (status < 0)
                {
                    (*MPIdata::serr) << "HDFrestart::readRestartRandomStates() "
                                        "--- H5Dread failed!!!"
                                     << std::endl;
                    return -2;
                }

                status = H5Dclose(dataset_id);
                if (status < 0)
                {
                    (*MPIdata::sout) << "H5Dclose failed!!!" << std::endl;
                    return -2;
                }

                if (!data.empty())
                    if (data[0] != data[0])
                    {
                        (*MPIdata::serr)
                            << "ERROR: HDFrestart::readRestartRandomStates() "
                               "--- data[0]="
                            << data[0] << std::endl;
                        return -2;
                    }
            }
        }
    }

    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

void HDFrestart::gatherDataXdir(std::vector<unsigned short>& data)
{

    int j      = pes_.my_mpi(1);
    int k      = pes_.my_mpi(2);
    int source = pes_.xyz2task(0, j, k);
    MPI_Status mpistatus;
    for (int i = 1; i < pes_.n_mpi_task(0); i++)
    {
        int dest = pes_.xyz2task(i, j, k);
        int tag  = i;
        if (pes_.my_mpi(0) == 0)
        {
            assert(active_);
            int datasize = data.size();
            MPI_Send(&datasize, 1, MPI_INT, dest, 2 * tag, comm_data_);
            MPI_Send(&data[0], datasize, MPI_UNSIGNED_SHORT, dest, 2 * tag + 1,
                comm_data_);
        }
        else if (pes_.my_mpi(0) == i)
        {
            assert(!active_);
            int datasize;
            MPI_Recv(
                &datasize, 1, MPI_INT, source, 2 * tag, comm_data_, &mpistatus);
            data.resize(datasize);
            MPI_Recv(&data[0], datasize, MPI_UNSIGNED_SHORT, source,
                2 * tag + 1, comm_data_, &mpistatus);
        }
    }
}

void HDFrestart::gatherDataXdir(std::vector<int>& data)
{

    int j      = pes_.my_mpi(1);
    int k      = pes_.my_mpi(2);
    int source = pes_.xyz2task(0, j, k);
    MPI_Status mpistatus;
    for (int i = 1; i < pes_.n_mpi_task(0); i++)
    {
        int dest = pes_.xyz2task(i, j, k);
        int tag  = i;
        if (pes_.my_mpi(0) == 0)
        {
            assert(active_);
            int datasize = data.size();
            MPI_Send(&datasize, 1, MPI_INT, dest, 2 * tag, comm_data_);
            MPI_Send(
                &data[0], datasize, MPI_INT, dest, 2 * tag + 1, comm_data_);
        }
        else if (pes_.my_mpi(0) == i)
        {
            assert(!active_);
            int datasize;
            MPI_Recv(
                &datasize, 1, MPI_INT, source, 2 * tag, comm_data_, &mpistatus);
            data.resize(datasize);
            MPI_Recv(&data[0], datasize, MPI_INT, source, 2 * tag + 1,
                comm_data_, &mpistatus);
        }
    }
}

void HDFrestart::gatherDataXdir(std::vector<double>& data)
{

    int j      = pes_.my_mpi(1);
    int k      = pes_.my_mpi(2);
    int source = pes_.xyz2task(0, j, k);
    MPI_Status mpistatus;
    for (int i = 1; i < pes_.n_mpi_task(0); i++)
    {
        int dest = pes_.xyz2task(i, j, k);
        int tag  = i;
        if (pes_.my_mpi(0) == 0)
        {
            assert(active_);
            int datasize = data.size();
            MPI_Send(&datasize, 1, MPI_INT, dest, 2 * tag, comm_data_);
            if (datasize > 0)
                MPI_Send(&data[0], datasize, MPI_DOUBLE, dest, 2 * tag + 1,
                    comm_data_);
        }
        else if (pes_.my_mpi(0) == i)
        {
            assert(!active_);
            int datasize;
            MPI_Recv(
                &datasize, 1, MPI_INT, source, 2 * tag, comm_data_, &mpistatus);
            if (datasize > 0)
            {
                data.resize(datasize);
                MPI_Recv(&data[0], datasize, MPI_DOUBLE, source, 2 * tag + 1,
                    comm_data_, &mpistatus);
            }
        }
    }
}

void HDFrestart::gatherDataXdir(std::vector<FixedLengthString>& data)
{
    int j      = pes_.my_mpi(1);
    int k      = pes_.my_mpi(2);
    int source = pes_.xyz2task(0, j, k);
    MPI_Status mpistatus;
    for (int i = 1; i < pes_.n_mpi_task(0); i++)
    {
        int dest = pes_.xyz2task(i, j, k);
        int tag  = i;
        if (pes_.my_mpi(0) == 0)
        {
            assert(active_);
            int datasize = data.size();
            std::vector<char> buffer(datasize * IonData_MaxStrLength);
            for (int i = 0; i < datasize; i++)
                memcpy(&buffer[i * IonData_MaxStrLength], &data[i],
                    IonData_MaxStrLength);

            MPI_Send(&datasize, 1, MPI_INT, dest, 2 * tag, comm_data_);
            MPI_Send(buffer.data(), datasize * IonData_MaxStrLength, MPI_CHAR,
                dest, 2 * tag + 1, comm_data_);
        }
        else if (pes_.my_mpi(0) == i)
        {
            assert(!active_);
            int datasize;
            MPI_Recv(
                &datasize, 1, MPI_INT, source, 2 * tag, comm_data_, &mpistatus);
            std::vector<char> buffer(datasize * IonData_MaxStrLength);
            MPI_Recv(buffer.data(), datasize * IonData_MaxStrLength, MPI_CHAR,
                source, 2 * tag + 1, comm_data_, &mpistatus);
            data.resize(datasize);
            for (int i = 0; i < datasize; i++)
                memcpy(&data[i], &buffer[i * IonData_MaxStrLength],
                    IonData_MaxStrLength);
        }
    }
}

void HDFrestart::gatherDataXdir(std::vector<char>& data)
{
    int j      = pes_.my_mpi(1);
    int k      = pes_.my_mpi(2);
    int source = pes_.xyz2task(0, j, k);
    MPI_Status mpistatus;
    for (int i = 1; i < pes_.n_mpi_task(0); i++)
    {
        int dest = pes_.xyz2task(i, j, k);
        int tag  = i;
        if (pes_.my_mpi(0) == 0)
        {
            assert(active_);
            int datasize = data.size();
            MPI_Send(&datasize, 1, MPI_INT, dest, 2 * tag, comm_data_);
            if (datasize > 0)
                MPI_Send(&data[0], datasize, MPI_CHAR, dest, 2 * tag + 1,
                    comm_data_);
        }
        else if (pes_.my_mpi(0) == i)
        {
            assert(!active_);
            int datasize;
            MPI_Recv(
                &datasize, 1, MPI_INT, source, 2 * tag, comm_data_, &mpistatus);
            if (datasize > 0)
            {
                data.resize(datasize);
                MPI_Recv(&data[0], datasize, MPI_CHAR, source, 2 * tag + 1,
                    comm_data_, &mpistatus);
            }
        }
    }
}
