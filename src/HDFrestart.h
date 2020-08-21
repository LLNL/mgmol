// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_HDFRESTART_H
#define MGMOL_HDFRESTART_H

#include "IonData.h"
#include "MPIdata.h"
#include "PEenv.h"
#include "Timer.h"
#include "Vector3D.h"

#include "hdf5.h"

#include "H5LTpublic.h"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

class LocalizationRegions;

herr_t H5LTfind_dataset(hid_t file_id_, const char* datasetname);

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

    template <class T>
    void gatherDataXdir(std::vector<T>& data);

    void closeWorkSpace();
    void setupWorkSpace();

    template <class T>
    void getWorkspace(T*& work_space);

    template <class T>
    int readDataset(hid_t dset_id, hid_t memspace, hid_t filespace,
        hid_t plist_id, T* work_space);

    template <class T>
    void MPI_Send_data(
        T*, const int n, const int dest, const int tag, MPI_Comm);

    template <class T>
    void MPI_Recv_data(T*, const int n, const int src, const int tag, MPI_Comm);

    template <class T>
    hid_t TH5Dcreate2(hid_t file_id, const std::string& datasetname,
        hid_t filespace, hid_t plist_id);

public:
    HDFrestart(const std::string& filename, const pb::PEenv& pes,
        const short option_number);
    HDFrestart(const std::string& filename, const pb::PEenv& pes,
        const unsigned gdim[3], const short option_number);

    bool gatherDataX() const { return gather_data_x_; }
    bool useHdf5p() const { return use_hdf5p_; }

    hid_t open_dset(const std::string& datasetname) const
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

    hid_t dset_exists(const std::string& datasetname) const
    {
        if (active_)
        {
            return dset_exists(datasetname.c_str());
        }
        return 0;
    }
    herr_t dset_exists(const char* const datasetname) const
    {
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
        return err_id;
    }

    void setVerbosity(const short verbosity) { verbosity_ = verbosity; }

    hid_t file_id() const { return file_id_; }
    std::string filename() const { return filename_; }
    bool active() const { return active_; }
    MPI_Comm& comm_active() { return comm_active_; }

    ~HDFrestart();
    int close();

    void setCommActive();

    hid_t createFilespace() const
    {
        assert(active_);

        hid_t filespace;
        if (use_hdf5p_)
            filespace = H5Screate_simple(3, dimsf_, nullptr);
        else
        {
            filespace = H5Screate_simple(3, block_, nullptr);
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
            memspace = H5Screate_simple(3, block_, nullptr);
        else
        {
            memspace = H5Screate_simple(3, block_, nullptr);
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
    int getLRs(std::shared_ptr<LocalizationRegions> lrs, const int max_nb,
        const std::string& name, const bool add = false);

    template <class T>
    int read_1func_hdf5(T*, const std::string&);

    template <class T>
    int write_1func_hdf5(
        T*, const std::string&, double* ll = nullptr, double* origin = nullptr);

    int read_att(const hid_t dset_id, const std::string& attname,
        std::vector<double>& attr_data);

    // write data in file with precision "precision"
    template <class T>
    int writeData(T* vv, hid_t filespace, hid_t memspace, hid_t dset_id,
        const short precision);

    template <class T>
    int readData(T* vv, hid_t memspace, hid_t dset_id, const short precision);

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

    void addDateToFilename();
    void addMDTime2File(const float run_time);
    void addMDstep2File(const int md_step);
    void add2File(const int data, const std::string& attname);

    float getMDTimeFromFile() const;
    int getMDstepFromFile() const;
    int getFromFile(const std::string& attname) const;

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
