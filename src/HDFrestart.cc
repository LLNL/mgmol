// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
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

#define MGMOL_HDFRESTART_FAIL(X)                                               \
    {                                                                          \
        (*MPIdata::serr) << "HDFrestart failure in file " << __FILE__          \
                         << " at line " << __LINE__ << std::endl;              \
        (*MPIdata::serr) << "Error Message: " << X << std::endl;               \
    }

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
        if (onpe0) MGMOL_HDFRESTART_FAIL("H5Fclose failed!!!");
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

    // seconds
    filename_.append("_");
    extension = std::to_string(tt1->tm_sec);
    if (tt1->tm_sec < 10)
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
#ifdef MGMOL_USE_HDF5P
    if (use_hdf5p_)
    {
        offset_[1] = pes_.my_mpi(1) * block_[1];
        offset_[2] = pes_.my_mpi(2) * block_[2];
    }
    else
#endif
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
            MGMOL_HDFRESTART_FAIL("H5Acreate failed!!!");
        }
        float attr_data = run_time;
        // if( onpe0 )
        {
            herr_t status
                = H5Awrite(attribute_id, H5T_NATIVE_FLOAT, &attr_data);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Awrite failed!!!");
            }
        }
        herr_t status = H5Sclose(dataspace_id);
        if (status < 0) MGMOL_HDFRESTART_FAIL("H5Sclose failed!!!");
        status = H5Aclose(attribute_id);
        if (status < 0) MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
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
            MGMOL_HDFRESTART_FAIL("H5Acreate failed!!!");
        }
        int attr_data = data;
        // if( onpe0 )
        {
            herr_t status = H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data);
            if (status < 0) MGMOL_HDFRESTART_FAIL("H5Awrite failed!!!");
        }
        herr_t status = H5Sclose(dataspace_id);
        if (status < 0) MGMOL_HDFRESTART_FAIL("H5Sclose failed!!!");
        status = H5Aclose(attribute_id);
        if (status < 0) MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
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
                    MGMOL_HDFRESTART_FAIL("H5Aread failed!!!");
                }
                else
                {
                    (*MPIdata::sout) << "MD time for restart file: " << run_time
                                     << std::endl;
                }

                status = H5Aclose(attribute_id);
                if (status < 0) MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
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
                    MGMOL_HDFRESTART_FAIL("H5Aread failed!!!");
                }
                else
                {
                    (*MPIdata::sout)
                        << attname << " in restart file: " << data << std::endl;
                }

                status = H5Aclose(attribute_id);
                if (status < 0) MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
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
        if (attribute_id < 0) MGMOL_HDFRESTART_FAIL("H5Acreate failed!!!");

        // if( onpe0 )
        {
            HDF_FixedLengthString t;
            strncpy(t.mystring, release, MyHDFStrLength - 1);
            herr_t status = H5Awrite(attribute_id, strtype, &t);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Awrite failed for " << attname);
            }
        }
        herr_t status = H5Sclose(dataspace_id);
        if (status < 0) MGMOL_HDFRESTART_FAIL("H5Sclose failed!!!");

        status = H5Aclose(attribute_id);
        if (status < 0) MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
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

    //(*MPIdata::sout)<<"HDFrestart::HDFrestart(),
    // filename="<<filename<<std::endl;
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
#ifdef MGMOL_USE_HDF5P
        if (use_hdf5p_) filename_.append(".h5");
#endif
    }

#ifdef MGMOL_USE_HDF5P
    if (!use_hdf5p_) // mkdir with filename and create files with task numbers
#endif
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
#ifdef MGMOL_USE_HDF5P
        if (use_hdf5p_)
        {
            // Set up file access property list with parallel I/O access
            herr_t err_id
                = H5Pset_fapl_mpio(access_plist, comm_active_, MPI_INFO_NULL);
            if (err_id < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Pset_fapl_mpio failed!!!");
            }
        }
        else
#endif
        {
            // The H5FD_CORE driver enables an application to work with a file
            // in memory, speeding reads and writes as no disk access is made.
            // File contents are stored only in memory until the file is closed.
            // The last parameter determines whether file contents are ever
            // written to disk.
            herr_t err_id = H5Pset_fapl_core(access_plist, 1024, 1);
            if (err_id < 0) MGMOL_HDFRESTART_FAIL("H5Pset_fapl_core failed!!!");
        }
        /* create the file collectively */
        H5Pset_fclose_degree(access_plist, H5F_CLOSE_STRONG);
        file_id_ = H5Fcreate(
            filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access_plist);
        if (file_id_ < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Fcreate failed for file " << filename_);
            mmpi.abort();
        }
        herr_t ret = H5Pclose(access_plist);
        if (ret < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Pclose failed!!!");
            mmpi.abort();
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
#ifdef MGMOL_USE_HDF5P
    if (!use_hdf5p_)
#endif
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
            MGMOL_HDFRESTART_FAIL(
                "File " << filename_ << " does not exists!!!");
            mmpi.abort();
        }
        htri_t ishdf = H5Fis_hdf5(filename_.c_str());
        if (ishdf < 0)
        {
            MGMOL_HDFRESTART_FAIL(
                "H5Fis_hdf5() unsuccessful for file " << filename_);
            mmpi.abort();
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
#ifdef MGMOL_USE_HDF5P
        if (use_hdf5p_)
        {
            // Set up file access property list with parallel I/O access
            herr_t err_id
                = H5Pset_fapl_mpio(access_plist, comm_active_, MPI_INFO_NULL);
            if (err_id < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Pset_fapl_mpio failed!!!");
                mmpi.abort();
            }
        }
        else
#endif
        {
            herr_t err_id = H5Pset_fapl_core(access_plist, 1024, 1);
            if (err_id < 0) MGMOL_HDFRESTART_FAIL("H5Pset_fapl_core failed!!!");
        }
        file_id_ = H5Fopen(filename_.c_str(), H5F_ACC_RDONLY, access_plist);
        H5Pclose(access_plist);

        if (file_id_ < 0)
        {
            MGMOL_HDFRESTART_FAIL("open " << filename_ << " failed!!!");
            mmpi.abort();
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
                MGMOL_HDFRESTART_FAIL("H5Sget_simple_extent_dims failed!!!");
                mmpi.abort();
            }

            // get global mesh size from file data
#ifdef MGMOL_USE_HDF5P
            if (!use_hdf5p_)
#endif
            {
                dimsf_[1] *= pes.n_mpi_task(1);
                dimsf_[2] *= pes.n_mpi_task(2);
            }
            if (!gather_data_x_) dimsf_[0] *= pes.n_mpi_task(0);

            // close dataset
            status = H5Dclose(dset_id);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Dclose failed!!!");
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
        MGMOL_HDFRESTART_FAIL(
            "H5Acreate failed for attribute " << attname << "!!!");
        return -1;
    }

    // if( onpe0 && ct.verbose>2 )
    //    (*MPIdata::sout)<<"Write attribute "<<attname<<std::endl;
    herr_t status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attr_data[0]);
    if (status < 0)
    {
        MGMOL_HDFRESTART_FAIL(
            "H5Awrite failed for attribute " << attname << "!!!");
        return -1;
    }
    status = H5Sclose(dataspace_id);
    if (status < 0)
    {
        MGMOL_HDFRESTART_FAIL("H5Sclose failed!!!");
        return -1;
    }
    status = H5Aclose(attribute_id);
    if (status < 0)
    {
        MGMOL_HDFRESTART_FAIL(
            "H5Aclose failed for attribute " << attname << "!!!");
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
        MGMOL_HDFRESTART_FAIL(
            "H5Aexists failed for attribute " << attname << "!!!");
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
        MGMOL_HDFRESTART_FAIL(
            "H5Aopen failed for attribute " << attname << "!!!");
        return -1;
    }

    hid_t attdataspace = H5Aget_space(attribute_id);
    hsize_t dims[2]    = { 0, 0 };
    hsize_t maxdims[2] = { 0, 0 };
    H5Sget_simple_extent_dims(attdataspace, dims, maxdims);
    if (dims[1] != 4)
    {
        MGMOL_HDFRESTART_FAIL("Wrong data size");
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
        MGMOL_HDFRESTART_FAIL("H5Aread failed!!!");
        return -1;
    }
    for (unsigned int j = 0; j < dims[0] * dims[1]; j++)
        assert(fabs(attr_data[j]) < 1.e30);
    status = H5Sclose(attdataspace);
    if (status < 0)
    {
        MGMOL_HDFRESTART_FAIL("H5Sclose failed!!!");
        return -1;
    }
    status = H5Aclose(attribute_id);
    if (status < 0)
    {
        MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
        return -1;
    }

    return natt;
}

int readGids(hid_t dset_id, std::vector<int>& gids)
{
    gids.clear();

    std::string attname("List of gids");
    htri_t exists = H5Aexists(dset_id, attname.c_str());
    if (exists < 0)
    {
        MGMOL_HDFRESTART_FAIL("H5Aexists() failed!!!");
        return exists;
    }
    if (exists == 0) // attribute does not exist
    {
        return 0;
    }

    // Open the data space for the attribute.
    //  Open a dataset attribute.
    hid_t attribute_id = H5Aopen_name(dset_id, attname.c_str());
    if (attribute_id < 0)
    {
        MGMOL_HDFRESTART_FAIL("H5Aopen failed for " << attname << "!!!");
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
        MGMOL_HDFRESTART_FAIL("H5Aread failed!!!");
        return -1;
    }
    status = H5Sclose(attdataspace);
    if (status < 0)
    {
        MGMOL_HDFRESTART_FAIL("H5Sclose failed!!!");
        return -1;
    }
    status = H5Aclose(attribute_id);
    if (status < 0)
    {
        MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
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

        int err_id = checkDataExistsLocal(datasetname);
        if (err_id == 0)
        { // dataset does not exists
            // try older version
            datasetname = getDatasetName_old(name, color);
            err_id      = checkDataExistsLocal(datasetname);
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
int HDFrestart::getLRs(std::shared_ptr<LocalizationRegions> lrs,
    const int max_nb_lrs, const std::string& name, const bool add)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("Reading localization regions...", (*MPIdata::sout));

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::getLRs() --- Read LRs from dataset "
                         << name << std::endl;

    if (!add) lrs->clear();

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

        int err_id = checkDataExistsLocal(datasetname);
        if (err_id == 0)
        { // dataset does not exists
            // try older version
            datasetname = getDatasetName_old(name, color);
            err_id      = checkDataExistsLocal(datasetname);
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
                lrs->push_back_local(center, rl, gids[i]);

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

        if (lrs->hasNcentersReachedNumber(max_nb_lrs)) done = 1;

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

template <>
void HDFrestart::getWorkspace<float>(float*& work_space)
{
    work_space = work_space_float_.data();
}

template <>
void HDFrestart::getWorkspace<double>(double*& work_space)
{
    work_space = work_space_double_.data();
}

template <>
int HDFrestart::readDataset(hid_t dset_id, hid_t memspace, hid_t filespace,
    hid_t plist_id, float* work_space)
{
    return H5Dread(
        dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, work_space);
}

template <>
int HDFrestart::readDataset(hid_t dset_id, hid_t memspace, hid_t filespace,
    hid_t plist_id, double* work_space)
{
    return H5Dread(
        dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, work_space);
}

template <>
void HDFrestart::MPI_Send_data(
    float* buffer, const int n, const int dest, const int tag, MPI_Comm comm)
{
    MPI_Send(buffer, n, MPI_FLOAT, dest, tag, comm);
}

template <>
void HDFrestart::MPI_Send_data(
    double* buffer, const int n, const int dest, const int tag, MPI_Comm comm)
{
    MPI_Send(buffer, n, MPI_DOUBLE, dest, tag, comm);
}

template <>
void HDFrestart::MPI_Send_data(
    int* buffer, const int n, const int dest, const int tag, MPI_Comm comm)
{
    MPI_Send(buffer, n, MPI_INT, dest, tag, comm);
}

template <>
void HDFrestart::MPI_Send_data(unsigned short* buffer, const int n,
    const int dest, const int tag, MPI_Comm comm)
{
    MPI_Send(buffer, n, MPI_UNSIGNED_SHORT, dest, tag, comm);
}

template <>
void HDFrestart::MPI_Send_data(
    char* buffer, const int n, const int dest, const int tag, MPI_Comm comm)
{
    MPI_Send(buffer, n, MPI_CHAR, dest, tag, comm);
}

template <>
void HDFrestart::MPI_Recv_data(
    float* buffer, const int n, const int src, const int tag, MPI_Comm comm)
{
    MPI_Status mpistatus;
    MPI_Recv(buffer, n, MPI_FLOAT, src, tag, comm, &mpistatus);
}

template <>
void HDFrestart::MPI_Recv_data(
    double* buffer, const int n, const int src, const int tag, MPI_Comm comm)
{
    MPI_Status mpistatus;
    MPI_Recv(buffer, n, MPI_DOUBLE, src, tag, comm, &mpistatus);
}

template <>
void HDFrestart::MPI_Recv_data(
    int* buffer, const int n, const int src, const int tag, MPI_Comm comm)
{
    MPI_Status mpistatus;
    MPI_Recv(buffer, n, MPI_INT, src, tag, comm, &mpistatus);
}

template <>
void HDFrestart::MPI_Recv_data(unsigned short* buffer, const int n,
    const int src, const int tag, MPI_Comm comm)
{
    MPI_Status mpistatus;
    MPI_Recv(buffer, n, MPI_UNSIGNED_SHORT, src, tag, comm, &mpistatus);
}

template <>
void HDFrestart::MPI_Recv_data(
    char* buffer, const int n, const int src, const int tag, MPI_Comm comm)
{
    MPI_Status mpistatus;
    MPI_Recv(buffer, n, MPI_CHAR, src, tag, comm, &mpistatus);
}

template <>
hid_t HDFrestart::TH5Dcreate2<double>(hid_t file_id,
    const std::string& datasetname, hid_t filespace, hid_t plist_id)
{
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(), H5T_NATIVE_DOUBLE,
        filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if (dset_id < 0)
    {
        MGMOL_HDFRESTART_FAIL("H5Dcreate failedi for " << datasetname << "!!!");
    }
    return dset_id;
}

template <>
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

template <class T>
int HDFrestart::read_1func_hdf5(T* vv, const std::string& datasetname)
{
    if (onpe0)
        (*MPIdata::sout) << "HDFrestart::read_1func_hdf5(). Try to read data "
                         << datasetname << "..." << std::endl;

    assert(block_[0] * block_[1] * block_[2] > 0);

    T* work_space = nullptr;

    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        hid_t dset_id  = open_dset(datasetname);
        if (dset_id < 0)
        {
            MGMOL_HDFRESTART_FAIL("open_dset() failed for " << datasetname);
            return -1;
        }

        // memory dataspace identifier
        hid_t memspace = createMemspace();
        if (memspace < 0)
        {
            MGMOL_HDFRESTART_FAIL("createMemspace failed!!!");
            return -1;
        }

        hid_t filespace = H5Dget_space(dset_id);
        if (filespace < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dget_space failed!!!");
            return -1;
        }

#ifdef MGMOL_USE_HDF5P
        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            herr_t status = H5Sselect_hyperslab(
                filespace, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Sselect_hyperslab failed!!!");
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
#endif
        // read data in work_space
        getWorkspace(work_space);
        herr_t status
            = readDataset(dset_id, memspace, filespace, plist_id, work_space);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dread failed!!!");
            return -1;
        }

#ifdef MGMOL_USE_HDF5P
        // Close/release resources.
        if (use_hdf5p_)
            if (pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);
#endif

        status = H5Dclose(dset_id);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dclose failed!!!");
            return -1;
        }

#ifdef MGMOL_USE_HDF5P
        if (use_hdf5p_)
        {
            status = H5Sclose(filespace);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Sclose failed!!!");
                return -1;
            }

            status = H5Sclose(memspace);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Sclose failed!!!");
                return -1;
            }
        }
#endif

        memcpy(vv, work_space, bsize_ * sizeof(T));

        if (onpe0)
            (*MPIdata::sout)
                << "Dataset " << datasetname << " read" << std::endl;
    } // if active_

    // send data to inactive PEs
    if (gather_data_x_)
    {
        int j      = pes_.my_mpi(1);
        int k      = pes_.my_mpi(2);
        int source = pes_.xyz2task(0, j, k);
        for (int i = 1; i < pes_.n_mpi_task(0); i++)
        {
            int dest = pes_.xyz2task(i, j, k);
            int tag  = i;
            if (pes_.my_mpi(0) == 0)
            {
                assert(active_);
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Send "<<bsize_<<"
                // data to "<<dest<<endl;
                MPI_Send_data(
                    work_space + i * bsize_, bsize_, dest, tag, comm_data_);
            }
            else if (pes_.my_mpi(0) == i)
            {
                assert(!active_);
                //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive
                //"<<bsize_<<" data from "<<source<<endl;
                MPI_Recv_data(vv, bsize_, source, tag, comm_data_);
            }
        }
    }

    return 0;
}

template <class T>
int HDFrestart::write_1func_hdf5(const T* const vv,
    const std::string& datasetname, double* ll, double* cell_origin)
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

        dset_id = TH5Dcreate2<T>(
            file_id_, datasetname.c_str(), filespace, plist_id);
        if (dset_id < 0)
        {
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
            MGMOL_HDFRESTART_FAIL("H5Dclose failed!!!");
            return -1;
        }

        status = H5Sclose(filespace);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Sclose failed for filespace!!!");
            return -1;
        }

        status = H5Sclose(memspace);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Sclose failed for memspace!!!");
            return -1;
        }
    }

    return 0;
}

template <class T>
int HDFrestart::readData(
    T* data, hid_t memspace, hid_t dset_id, const short precision)
{
    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        herr_t status;
        hid_t filespace = H5Dget_space(dset_id);
        if (filespace < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dget_space failed !!!");
            return -1;
        }
#ifdef MGMOL_USE_HDF5P
        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            status = H5Sselect_hyperslab(
                filespace, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Sselect_hyperslab failed!!!");
                return -1;
            }

            // Create property list for collective dataset read.
            if (pes_.n_mpi_tasks() > 1)
            {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
            }
        }
#endif
        if (precision == 1)
        {
            status = H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace,
                plist_id, work_space_float_.data());
        }
        else
        {
            status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                plist_id, work_space_double_.data());
        }

        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dread  failed !!!");
            return -1;
        }

#ifdef MGMOL_USE_HDF5P
        // Close/release resources.
        if (use_hdf5p_ && pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);
#endif

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
                    MPI_Send(work_space_double_.data() + i * bsize_, bsize_,
                        MPI_DOUBLE, dest, tag, comm_data_);
                }
                else if (pes_.my_mpi(0) == i)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data
                    // from "<<ipe<<endl;
                    MPI_Recv(work_space_double_.data(), bsize_, MPI_DOUBLE,
                        source, tag, comm_data_, &mpistatus);
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
                    MPI_Send(work_space_float_.data() + i * bsize_, bsize_,
                        MPI_FLOAT, dest, tag, comm_data_);
                }
                else if (pes_.my_mpi(0) == i)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data
                    // from "<<ipe<<endl;
                    MPI_Recv(work_space_float_.data(), bsize_, MPI_FLOAT,
                        source, tag, comm_data_, &mpistatus);
                }
            }
    }

    // assign data array passed as argument to this function
    if (precision == 1)
    {
        for (int i = 0; i < bsize_; i++)
            data[i] = (T)work_space_float_[i];
    }
    else
    {
        for (int i = 0; i < bsize_; i++)
            data[i] = (T)work_space_double_[i];
    }
    return 0;
}

template <class T>
int HDFrestart::writeData(const T* const data, hid_t space_id, hid_t memspace,
    hid_t dset_id, const short precision)
{
    if (precision == 1)
    {
        assert((int)work_space_float_.size() == bsize_);
        for (int i = 0; i < bsize_; i++)
            work_space_float_[i] = (float)data[i];
    }
    else
    {
        assert((int)work_space_double_.size() == bsize_);
        for (int i = 0; i < bsize_; i++)
            work_space_double_[i] = (double)data[i];
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
                    MPI_Send(work_space_double_.data(), bsize_, MPI_DOUBLE,
                        dest, i, comm_data_);
                }
                else if (pes_.my_mpi(0) == 0)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data
                    // from "<<ipe<<endl;
                    MPI_Recv(work_space_double_.data() + i * bsize_, bsize_,
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
                    MPI_Send(work_space_float_.data(), bsize_, MPI_FLOAT, dest,
                        i, comm_data_);
                }
                else if (pes_.my_mpi(0) == 0)
                {
                    //(*MPIdata::sout)<<"PE: "<<pes_.mytask()<<", Receive data
                    // from "<<ipe<<endl;
                    MPI_Recv(work_space_float_.data() + i * bsize_, bsize_,
                        MPI_FLOAT, ipe, i, comm_data_, &mpistatus);
                }
            }
    }

    if (active_)
    {
        hid_t plist_id = H5P_DEFAULT;
        herr_t status;
#ifdef MGMOL_USE_HDF5P
        if (use_hdf5p_)
        {
            // Select hyperslab in the file.
            status = H5Sselect_hyperslab(
                space_id, H5S_SELECT_SET, offset_, stride_, count_, block_);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Sselect_hyperslab failed for space_id="
                                      << space_id << "!!!");
                return -1;
            }

            // Create property list for collective dataset write.
            if (pes_.n_mpi_tasks() > 1)
            {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
            }
        }
#endif
        assert(dset_id >= 0);
        if (precision == 1)
            status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, space_id,
                plist_id, work_space_float_.data());
        else
            status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id,
                plist_id, work_space_double_.data());
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dwrite failed!!!");
            return -1;
        }

        // Close/release resources.
#ifdef MGMOL_USE_HDF5P
        if (use_hdf5p_ && pes_.n_mpi_tasks() > 1) H5Pclose(plist_id);
#endif
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
            MGMOL_HDFRESTART_FAIL("H5Aopen_name failed for " << attname);
            return -1;
        }
        herr_t status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &attr_data[0]);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Aread failed!!!");
            return -1;
        }
        status = H5Aclose(attribute_id);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
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
            MGMOL_HDFRESTART_FAIL("ERROR in MPI_Comm_split!!!");
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
#ifdef MGMOL_USE_HDF5P
            use_hdf5p_ = false;
#endif
            gather_data_x_ = false;
            break;
#ifdef MGMOL_USE_HDF5P
        case 1: // 1 file for all tasks
            use_hdf5p_     = true;
            gather_data_x_ = true;
            break;
#endif
        case 2: // 1 file per column in x-direction
#ifdef MGMOL_USE_HDF5P
            use_hdf5p_ = false;
#endif
            gather_data_x_ = true;
            break;
        default:
            MGMOL_HDFRESTART_FAIL(
                "setOptions() --- type " << option_number << " undefined!!!");
            return;
    }
}

void HDFrestart::setupWorkSpace()
{
    // if( active_ )
    {
        const int n = block_[0] * block_[1] * block_[2];
        work_space_double_.resize(n);
        memset(work_space_double_.data(), 0, n * sizeof(double));

        work_space_float_.resize(n);
        memset(work_space_float_.data(), 0, n * sizeof(float));
    }
}

void HDFrestart::printTimers(std::ostream& os)
{
    open_existing_tm_.print(os);
    create_file_tm_.print(os);
    close_file_tm_.print(os);
}

int HDFrestart::readAtomicData(std::string datasetname, std::vector<int>& data)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::readAtomicData()..." << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (!exists) return 0;

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (dataset_id < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dopen2 failed for " + datasetname);
            return -1;
        }

        int dim = (int)(H5Dget_storage_size(dataset_id) / sizeof(int));
        if (dim == 0)
        {
            MGMOL_HDFRESTART_FAIL("No " + datasetname);
            return -1;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dread failed for " + datasetname);
            return -1;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dclose failed for " + datasetname);
            return -1;
        }
    } // if active_
    // send data to inactive PEs
    if (gather_data_x_) gatherDataXdir(data);

#ifdef MGMOL_USE_HDF5P
    if (useHdf5p())
    {
        data.erase(std::remove(data.begin(), data.end(), -1), data.end());
    }
#endif

    return 0;
}

int HDFrestart::readAtomicData(
    std::string datasetname, std::vector<double>& data)
{
    if (onpe0)
        (*MPIdata::sout) << "Read atomic data from hdf5 file" << std::endl;

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (!exists) return -1;

        // Open the dataset
        hid_t dataset_id = H5Dopen2(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (dataset_id < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dopen2 failed for " + datasetname);
            return -2;
        }

        int dim = (int)H5Dget_storage_size(dataset_id) / sizeof(double);
        if (dim == 0)
        {
            MGMOL_HDFRESTART_FAIL("readAtomicData() --- No " + datasetname);
            return -2;
        }
        data.resize(dim);

        herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, &data[0]);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dread failed for " + datasetname);
            return -2;
        }
        status = H5Dclose(dataset_id);
        if (status < 0)
        {
            MGMOL_HDFRESTART_FAIL("H5Dclose failed for " + datasetname);
            return -2;
        }
    }

#ifdef MGMOL_USE_HDF5P
    if (useHdf5p())
    {
        data.erase(std::remove(data.begin(), data.end(), 1e+32), data.end());
    }
#endif
    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

int HDFrestart::readOldCenterOnMesh(std::vector<double>& data, int i)
{
    if (onpe0)
        (*MPIdata::sout) << "Read old localization centers rounded to mesh "
                            "points from hdf5 file"
                         << std::endl;

    std::stringstream datasetstream;
    datasetstream << "OldCenterOnMesh_" << i;

    std::string datasetname = datasetstream.str();

    return readAtomicData(datasetname, data);
}

int HDFrestart::readOldCenter(std::vector<double>& data, int i)
{
    if (onpe0)
        (*MPIdata::sout) << "Read old localization centers from hdf5 file"
                         << std::endl;

    std::stringstream datasetstream;
    datasetstream << "OldCenter_" << i;

    std::string datasetname = datasetstream.str();

    return readAtomicData(datasetname, data);
}

int HDFrestart::readAtomicData(
    std::string datasetname, std::vector<std::string>& data)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "HDFrestart::readAtomicData(), dataset = "
                         << datasetname << std::endl;

    std::vector<char> buffer;
    short name_length = 7; // default, value used before February 2016

    if (active_)
    {
        assert(file_id_ >= 0);

        htri_t exists = H5Lexists(file_id_, datasetname.c_str(), H5P_DEFAULT);
        if (exists)
        {
            // Open the dataset
            hid_t dataset_id
                = H5Dopen2(file_id_, datasetname.c_str(), H5P_DEFAULT);
            if (dataset_id < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Dopen2 failed!!!");
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
                    MGMOL_HDFRESTART_FAIL("H5Aread failed!!!");
                    return -1;
                }

                status = H5Aclose(attribute_id);
                if (status < 0)
                {
                    MGMOL_HDFRESTART_FAIL("H5Aclose failed!!!");
                    return -1;
                }
            }

            int dim = (int)H5Dget_storage_size(dataset_id) / name_length;
            if (dim == 0)
            {
                MGMOL_HDFRESTART_FAIL("No names!!!");
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
                MGMOL_HDFRESTART_FAIL("H5Dread failed!!!");
                return -1;
            }
            status = H5Dclose(dataset_id);
            if (status < 0)
            {
                MGMOL_HDFRESTART_FAIL("H5Dclose failed!!!");
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

        stripLeadingAndTrailingBlanks(t);
        // std::cout<<"stripped name="<<t<<std::endl;

        data.push_back(t);
    }

#ifdef MGMOL_USE_HDF5P
    if (useHdf5p())
    {
        data.erase(std::remove(data.begin(), data.end(), ""), data.end());
    }
#endif

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
                    MGMOL_HDFRESTART_FAIL("H5Dread failed!!!");
                    return -2;
                }

                status = H5Dclose(dataset_id);
                if (status < 0)
                {
                    MGMOL_HDFRESTART_FAIL("H5Dclose failed!!!");
                    return -2;
                }

                if (!data.empty())
                    if (std::isnan(data[0]))
                    {
                        MGMOL_HDFRESTART_FAIL(
                            "readRestartRandomStates() is NaN");
                        return -2;
                    }
            }
        }
    }

    if (gather_data_x_) gatherDataXdir(data);

    return 0;
}

template <class T>
void HDFrestart::gatherDataXdir(std::vector<T>& data)
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
                MPI_Send_data(
                    &data[0], datasize, dest, 2 * tag + 1, comm_data_);
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
                MPI_Recv_data(
                    &data[0], datasize, source, 2 * tag + 1, comm_data_);
            }
        }
    }
}

int HDFrestart::countFunctionObjects(std::string& name) const
{
    int count = 0;
    int found = 0;
    do
    {
        std::string datasetname(getDatasetName(name, count));
        // check if dataset exists...
        found = checkDataExists(datasetname);
        if (found) count++;
    } while (found); // dataset exists

    return count;
}

template int HDFrestart::read_1func_hdf5(float*, const std::string&);
template int HDFrestart::read_1func_hdf5(double*, const std::string&);

template int HDFrestart::write_1func_hdf5(
    const double* const, const std::string&, double* ll, double* origin);

template int HDFrestart::readData(
    double*, hid_t memspace, hid_t dset_id, const short precision);
template int HDFrestart::readData(
    float*, hid_t memspace, hid_t dset_id, const short precision);
template int HDFrestart::writeData(const double* const vv, hid_t filespace,
    hid_t memspace, hid_t dset_id, const short precision);
template int HDFrestart::writeData(const float* const vv, hid_t filespace,
    hid_t memspace, hid_t dset_id, const short precision);
