// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef HDFRESTART_H
#define HDFRESTART_H

#include "IonData.h"
#include "MPIdata.h"
#include "PEenv.h"
#include "Timer.h"
#include "Vector3D.h"

#include "hdf5.h"

#ifndef USE_HDF16
#include "H5LTpublic.h"
#endif

#include <iostream>
#include <map>
#include <string>
#include <vector>

class LocalizationRegions;

#ifdef USE_HDF16
hid_t H5Dopen2(hid_t loc_id, const char* name, hid_t dapl_id);
hid_t H5Dcreate2(hid_t loc_id, const char* name, hid_t dtype_id, hid_t space_id,
    hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id);
hid_t H5Acreate2(hid_t loc_id, const char* attr_name, hid_t type_id,
    hid_t space_id, hid_t acpl_id, hid_t aapl_id);
#endif

int readListCentersAndRadii(hid_t dset_id, std::vector<double>& attr_data);
int readGids(hid_t dset_id, std::vector<int>& gids);
int writeListCentersAndRadii(
    hid_t dset_id, const unsigned natt, const std::vector<double>& attr_data);
int writeGids(hid_t dset_id, const std::vector<int>& gids);

class HDFrestart
{
    static Timer open_existing_tm_;
    static Timer create_file_tm_;
    static Timer close_file_tm_;

    const pb::PEenv& pes_;

    std::string filename_;
    hid_t file_id_;

    // flag to denote tasks writing (or not) directly to HDF file
    bool active_;

    hsize_t block_[3];
    hsize_t offset_[3];
    hsize_t count_[3];
    hsize_t stride_[3];

    // global mesh dimensions
    hsize_t dimsf_[3];

    int bsize_;

    double* work_space_double_;
    float* work_space_float_;

    bool use_hdf5p_;
    bool gather_data_x_;

    short verbosity_;

    bool closed_;

    // MPI communicator with all tasks with flag active_=true
    MPI_Comm comm_active_;

    // MPI communicator for all tasks having data to write
    MPI_Comm comm_data_;

    void appendTaskNumberToFilename();
    void setActivity();
    void setupBlocks();
    void setOptions(const short option_number);
    void addReleaseNumber2File(const char* release);
    void gatherDataXdir(std::vector<int>& data);
    void gatherDataXdir(std::vector<unsigned short>& data);
    void gatherDataXdir(std::vector<double>& data);
    void gatherDataXdir(std::vector<FixedLengthString>& data);
    void gatherDataXdir(std::vector<char>& data);

    void closeWorkSpace();
    void setupWorkSpace();

public:
    HDFrestart(const std::string filename, const pb::PEenv& pes,
        const short option_number);
    HDFrestart(const std::string filename, const pb::PEenv& pes,
        const unsigned gdim[3], const short option_number);

    bool gatherDataX() const { return gather_data_x_; }
    bool useHdf5p() const { return use_hdf5p_; }

    hid_t open_dset(const std::string datasetname) const
    {
        return open_dset(datasetname.c_str());
    }
    hid_t open_dset(const char* const datasetname) const
    {
        hid_t dset_id = 0;
        if (active_)
        {
            assert(file_id_ >= 0);
            dset_id = H5Dopen2(file_id_, datasetname, H5P_DEFAULT);
            if (dset_id < 0)
            {
                if (onpe0)
                    (*MPIdata::sout)
                        << "HDFrestart::open_dset() --- Cannot open Dataset "
                        << datasetname << std::endl;
            }
        }
        return dset_id;
    }
    hid_t close_dset(const hid_t dset_id)
    {
        if (active_)
        {
            herr_t status = H5Dclose(dset_id);
            if (status < 0)
            {
                if (onpe0)
                    (*MPIdata::serr) << "H5Dclose failed!!!" << std::endl;
                return -1;
            }
        }
        return 0;
    }

    hid_t dset_exists(const std::string datasetname) const
    {
        if (active_)
        {
            return dset_exists(datasetname.c_str());
        }
        return 0;
    }
    herr_t dset_exists(const char* const datasetname) const
    {
#ifdef USE_HDF16
        herr_t err_id = 1;
#else
        herr_t err_id = 0;
        if (active_)
        {
            err_id = H5LTfind_dataset(file_id_, datasetname);

            if (err_id < 0)
            {
                if (onpe0)
                    (*MPIdata::sout)
                        << "HDFrestart::dset_exists() failed for dataset "
                        << datasetname << std::endl;
            }
        }
#endif
        return err_id;
    }

    void setVerbosity(const short verbosity) { verbosity_ = verbosity; }

    hid_t file_id() const { return file_id_; }
    std::string filename() const { return filename_; }
    bool active() const { return active_; }
#ifdef USE_MPI
    MPI_Comm& comm_active() { return comm_active_; }
#endif

    ~HDFrestart();
    int close();

    void setCommActive();

    hid_t createFilespace() const
    {
        assert(active_);

        hid_t filespace;
        if (use_hdf5p_)
            filespace = H5Screate_simple(3, dimsf_, NULL);
        else
        {
            filespace = H5Screate_simple(3, block_, NULL);
            // hsize_t dims[1]={block_[0]*block_[1]*block_[2]};
            // filespace = H5Screate_simple(1, dims, NULL);
        }
        if (filespace < 0 && onpe0)
            (*MPIdata::serr)
                << "ERROR: createFilespace() failed!!!" << std::endl;
        return filespace;
    }
    hid_t createMemspace() const
    {
        assert(active_);

        hid_t memspace;
        if (use_hdf5p_)
            memspace = H5Screate_simple(3, block_, NULL);
        else
        {
            memspace = H5Screate_simple(3, block_, NULL);
            // hsize_t dims[1]={block_[0]*block_[1]*block_[2]};
            // memspace=H5Screate_simple(1, dims, NULL);
        }
        if (memspace < 0 && onpe0)
            (*MPIdata::serr)
                << "ERROR: createMemspace() failed!!!" << std::endl;
        return memspace;
    }

    int getLRCenters(std::multimap<std::string, Vector3D>&, const int,
        const std::string& name);
    int getLRs(LocalizationRegions& lrs, const int max_nb,
        const std::string& name, const bool add = false);
    int read_1func_hdf5(double*, std::string);
    int read_1func_hdf5(float*, std::string);
    int write_1func_hdf5(
        double*, std::string, double* ll = NULL, double* origin = NULL);
    int write_1func_hdf5(
        float*, std::string, double* ll = NULL, double* origin = NULL);
    int read_att(const hid_t dset_id, std::string attname,
        std::vector<double>& attr_data);

    int writeData(double* vv, hid_t filespace, hid_t memspace, hid_t dset_id,
        const short precision);
    int writeData(float* vv, hid_t filespace, hid_t memspace, hid_t dset_id,
        const short precision);
    int readData(
        double* vv, hid_t memspace, hid_t dset_id, const short precision);
    int readData(
        float* vv, hid_t memspace, hid_t dset_id, const short precision);
    //    int writeRandomState(unsigned short int rand_state[3]);
    //    int readRandomState(unsigned short* rand_state);
    int readAtomicIDs(std::vector<int>& data);
    int readAtomicNLprojIDs(std::vector<int>& data);
    int readAtomicNumbers(std::vector<int>& data);
    int readAtomicNames(std::vector<std::string>& data);
    int readAtomicPositions(std::vector<double>& data);
    int readAtomicVelocities(std::vector<double>& data);
    int readLockedAtomNames(std::vector<std::string>& data);
    int readRestartRandomStates(std::vector<unsigned short>& data);
    int readOldCenter(std::vector<double>& data, int i);
    int readOldCenterOnMesh(std::vector<double>& data, int i);
    int readGidsList(std::vector<int>& data);

    bool olderVersion();

    void addDateToFilename();
    void addMDTime2File(const float run_time);
    void addMDstep2File(const int md_step);
    void add2File(const int data, std::string attname);

    float getMDTimeFromFile() const;
    int getMDstepFromFile() const;
    int getFromFile(std::string attname) const;

    hid_t createPlist()
    {
        hid_t plist_id = H5P_DEFAULT;
        if (use_hdf5p_ && active_)
        {
            // Create chunked dataset.
            plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 3, block_);
        }

        return plist_id;
    }

    void releasePlist(hid_t plist_id)
    {
        if (use_hdf5p_ && active_) H5Pclose(plist_id);
    }

    static void printTimers(std::ostream& os);
};

#endif
