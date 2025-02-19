// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "hdf_tools.h"
#include "IonData.h"

#include <cassert>
#include <iostream>
#include <string.h>

#define MGMOL_HDF5_FAIL(X)                                                     \
    {                                                                          \
        std::cerr << "MGMOL_HDF5 failure:" << std::endl;                       \
        std::cerr << "Error Message: " << X << std::endl;                      \
    }

namespace mgmol_tools
{
void string2fixedlength(
    std::vector<std::string>& data, std::vector<FixedLengthString>& tc)
{
    tc.clear();
    for (auto& d : data)
    {
        FixedLengthString t;
        strncpy(t.mystring, d.c_str(), IonData_MaxStrLength - 1);
        tc.push_back(t);
    }
}

void write1d(hid_t file_id, const std::string& datasetname,
    std::vector<int>& data, size_t length)
{
    assert(file_id >= 0);

    if (length == 0) return;

    hsize_t dim = length;

    // Create the data space for the datasets
    hid_t dataspace_id = H5Screate_simple(1, &dim, nullptr);
    if (dataspace_id < 0)
    {
        std::cerr << "write1d(), H5Screate_simple failed!!!" << std::endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(), H5T_NATIVE_INT,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        std::cerr << "write1d(), H5Dcreate2 failed!!!" << std::endl;
        return;
    }
    H5Sclose(dataspace_id);

    herr_t status = H5Dwrite(
        dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    if (status < 0)
    {
        std::cerr << "write1d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    status = H5Dclose(dataset_id);
    if (status < 0)
    {
        std::cerr << "write1d(), H5Dclose failed!!!" << std::endl;
        return;
    }
}

void write2d(hid_t file_id, const std::string& datasetname,
    std::vector<int>& data, size_t* dims)
{
    assert(file_id >= 0);

    // std::cout<<"Write "<<dims[0]<<" atomic numbers..."<<endl;

    if (dims[0] == 0) return;

    hsize_t dimsm[2] = { (hsize_t)dims[0], (hsize_t)dims[1] };

    // Create the data space for the datasets
    hid_t dataspace_id = H5Screate_simple(2, dimsm, nullptr);
    if (dataspace_id < 0)
    {
        std::cerr << "write2d(), H5Screate_simple failed!!!" << std::endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(), H5T_NATIVE_INT,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        std::cerr << "write2d(), H5Dcreate2 failed!!!" << std::endl;
        return;
    }
    H5Sclose(dataspace_id);

    herr_t status = H5Dwrite(
        dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    status = H5Dclose(dataset_id);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Dclose failed!!!" << std::endl;
        return;
    }
}

void write2d(hid_t file_id, const std::string& datasetname,
    std::vector<unsigned short>& data, size_t* dims)
{
    assert(file_id >= 0);

    // std::cout<<"Write "<<dims[0]<<" atomic numbers..."<<endl;

    if (dims[0] == 0) return;

    hsize_t dimsm[2] = { (hsize_t)dims[0], (hsize_t)dims[1] };

    // Create the data space for the datasets
    hid_t dataspace_id = H5Screate_simple(2, dimsm, nullptr);
    if (dataspace_id < 0)
    {
        std::cerr << "write2d(), H5Screate_simple failed!!!" << std::endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(),
        H5T_NATIVE_USHORT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        std::cerr << "write2d(), H5Dcreate2 failed!!!" << std::endl;
        return;
    }
    H5Sclose(dataspace_id);

    herr_t status = H5Dwrite(
        dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    status = H5Dclose(dataset_id);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Dclose failed!!!" << std::endl;
        return;
    }
}

void write2d(hid_t file_id, const std::string& datasetname,
    std::vector<double>& data, size_t* dims)
{
    assert(file_id >= 0);

    // std::cout<<"Write "<<dims[0]<<" atomic numbers..."<<endl;

    if (dims[0] == 0) return;

    hsize_t dimsm[2] = { (hsize_t)dims[0], (hsize_t)dims[1] };

    // Create the data space for the datasets
    hid_t dataspace_id = H5Screate_simple(2, dimsm, nullptr);
    if (dataspace_id < 0)
    {
        std::cerr << "write2d(), H5Screate_simple failed!!!" << std::endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(),
        H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        std::cerr << "write2d(), H5Dcreate2 failed!!!" << std::endl;
        return;
    }
    H5Sclose(dataspace_id);

    herr_t status = H5Dwrite(
        dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    status = H5Dclose(dataset_id);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Dclose failed!!!" << std::endl;
        return;
    }
}

void write2d(hid_t file_id, const std::string& datasetname,
    std::vector<std::string>& data, size_t* dims)
{
    assert(file_id >= 0);

    // create type for std::strings of length IonData_MaxStrLength
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, IonData_MaxStrLength);

    // std::cout<<"Write "<<dims[0]<<" atomic numbers..."<<endl;

    if (dims[0] == 0) return;

    hsize_t dimsm[2] = { (hsize_t)dims[0], (hsize_t)dims[1] };

    // Create the data space for the datasets
    hid_t dataspace_id = H5Screate_simple(2, dimsm, nullptr);
    if (dataspace_id < 0)
    {
        std::cerr << "write2d(), H5Screate_simple failed!!!" << std::endl;
        return;
    }

    // Create dataset
    hid_t dataset_id = H5Dcreate2(file_id, datasetname.c_str(), strtype,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0)
    {
        std::cerr << "write2d(), H5Dcreate2 failed!!!" << std::endl;
        return;
    }
    H5Sclose(dataspace_id);

    // First copy the contents of the vector into a temporary container
    std::vector<FixedLengthString> tc;
    string2fixedlength(data, tc);

    std::string attname("String_Length");
    hsize_t dimsA[1]    = { 1 };
    hid_t dataspaceA_id = H5Screate_simple(1, dimsA, nullptr);
    hid_t attribute_id = H5Acreate2(dataset_id, attname.c_str(), H5T_NATIVE_INT,
        dataspaceA_id, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status
        = H5Awrite(attribute_id, H5T_NATIVE_USHORT, &IonData_MaxStrLength);
    if (status < 0)
    {
        std::cerr << "write2d(), Attribute: " << attname
                  << " --- H5Awrite failed!!!" << std::endl;
    }

    status
        = H5Dwrite(dataset_id, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tc[0]);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    H5Tclose(strtype);

    status = H5Sclose(dataspaceA_id);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Sclose failed!!!" << std::endl;
    }

    status = H5Aclose(attribute_id);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Aclose failed!!!" << std::endl;
    }

    status = H5Dclose(dataset_id);
    if (status < 0)
    {
        std::cerr << "write2d(), H5Dclose failed!!!" << std::endl;
    }
}

#ifdef MGMOL_USE_HDF5P
void parallelWrite2d(hid_t file_id, const std::string& datasetname,
    std::vector<int>& data, size_t* dims, MPI_Comm comm)
{
    assert(file_id >= 0);
    assert(!data.empty());

    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    hsize_t dimsm[2] = { (hsize_t)dims[0], (hsize_t)dims[1] };
    hsize_t dimsf[2] = { (hsize_t)(dimsm[0] * mpi_size), dimsm[1] };

    hid_t filespace = H5Screate_simple(2, dimsf, nullptr);
    if (filespace < 0)
    {
        std::cerr
            << "parallelWrite2d(), H5Screate_simple failed for filespace!!!"
            << std::endl;
        return;
    }
    hid_t memspace = H5Screate_simple(2, dimsm, nullptr);
    if (memspace < 0)
    {
        std::cerr
            << "parallelWrite2d(), H5Screate_simple failed for memspace!!!"
            << std::endl;
        return;
    }

    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, dimsm);
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(), H5T_NATIVE_INT,
        filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if (dset_id < 0)
    {
        std::cerr << "parallelWrite2d() for dataset " << datasetname
                  << ", H5Dcreate2() failed!!!" << std::endl;
        return;
    }
    H5Pclose(plist_id);

    hsize_t offset[2] = { mpi_rank * dimsm[0], 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t count[2]  = { 1, 1 };
    hsize_t block[2]  = { dimsm[0], dimsm[1] };

    /* Select hyperslab in the file. */
    herr_t status = H5Sselect_hyperslab(
        filespace, H5S_SELECT_SET, offset, stride, count, block);
    if (status < 0)
    {
        std::cerr << "parallelWrite2d(), H5Sselect_hyperslab() failed!!!"
                  << std::endl;
        return;
    }

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(
        dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &data[0]);
    if (status < 0)
    {
        std::cerr << "parallelWrite2d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

void parallelWrite2d(hid_t file_id, const std::string& datasetname,
    std::vector<unsigned short>& data, size_t* dims, MPI_Comm comm)
{
    assert(file_id >= 0);
    assert(!data.empty());

    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    hsize_t dimsm[2] = { (hsize_t)dims[0], (hsize_t)dims[1] };
    hsize_t dimsf[2] = { (hsize_t)(dimsm[0] * mpi_size), dimsm[1] };

    hid_t filespace = H5Screate_simple(2, dimsf, nullptr);
    if (filespace < 0)
    {
        std::cerr
            << "parallelWrite2d(), H5Screate_simple failed for filespace!!!"
            << std::endl;
        return;
    }
    hid_t memspace = H5Screate_simple(2, dimsm, nullptr);
    if (memspace < 0)
    {
        std::cerr
            << "parallelWrite2d(), H5Screate_simple failed for memspace!!!"
            << std::endl;
        return;
    }

    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, dimsm);
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(), H5T_NATIVE_USHORT,
        filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if (dset_id < 0)
    {
        std::cerr << "parallelWrite2d() for dataset " << datasetname
                  << ", H5Dcreate2() failed!!!" << std::endl;
        return;
    }
    H5Pclose(plist_id);

    hsize_t offset[2] = { mpi_rank * dimsm[0], 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t count[2]  = { 1, 1 };
    hsize_t block[2]  = { dimsm[0], dimsm[1] };

    /* Select hyperslab in the file. */
    herr_t status = H5Sselect_hyperslab(
        filespace, H5S_SELECT_SET, offset, stride, count, block);
    if (status < 0)
    {
        std::cerr << "parallelWrite2d(), H5Sselect_hyperslab() failed!!!"
                  << std::endl;
        return;
    }

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(
        dset_id, H5T_NATIVE_USHORT, memspace, filespace, plist_id, &data[0]);
    if (status < 0)
    {
        std::cerr << "parallelWrite2d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

void parallelWrite2d(hid_t file_id, const std::string& datasetname,
    std::vector<double>& data, size_t* dims, MPI_Comm comm)
{
    assert(file_id >= 0);
    assert(!data.empty());

    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    hsize_t dimsm[2] = { (hsize_t)dims[0], (hsize_t)dims[1] };
    hsize_t dimsf[2] = { (hsize_t)(dimsm[0] * mpi_size), dimsm[1] };

    hid_t filespace = H5Screate_simple(2, dimsf, nullptr);
    if (filespace < 0)
    {
        std::cerr
            << "parallelWrite2d(), H5Screate_simple failed for filespace!!!"
            << std::endl;
        return;
    }
    hid_t memspace = H5Screate_simple(2, dimsm, nullptr);
    if (memspace < 0)
    {
        std::cerr
            << "parallelWrite2d(), H5Screate_simple failed for memspace!!!"
            << std::endl;
        return;
    }

    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, dimsm);
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(), H5T_NATIVE_DOUBLE,
        filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if (dset_id < 0)
    {
        std::cerr << "parallelWrite2d() for dataset " << datasetname
                  << ", H5Dcreate2() failed!!!" << std::endl;
        return;
    }
    H5Pclose(plist_id);

    hsize_t offset[2] = { mpi_rank * dimsm[0], 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t count[2]  = { 1, 1 };
    hsize_t block[2]  = { dimsm[0], dimsm[1] };

    /* Select hyperslab in the file. */
    herr_t status = H5Sselect_hyperslab(
        filespace, H5S_SELECT_SET, offset, stride, count, block);
    if (status < 0)
    {
        std::cerr << "parallelWrite2d(), H5Sselect_hyperslab() failed!!!"
                  << std::endl;
        return;
    }

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(
        dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &data[0]);
    if (status < 0)
    {
        std::cerr << "parallelWrite2d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

void parallelWrite2d(hid_t file_id, const std::string& datasetname,
    std::vector<std::string>& data, size_t* dims, MPI_Comm comm)
{
    assert(file_id >= 0);
    assert(!data.empty());

    // create type for std::strings of length IonData_MaxStrLength
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, IonData_MaxStrLength);

    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    hsize_t dimsm[2] = { (hsize_t)dims[0], (hsize_t)dims[1] };
    hsize_t dimsf[2] = { (hsize_t)(dimsm[0] * mpi_size), dimsm[1] };

    hid_t filespace = H5Screate_simple(2, dimsf, nullptr);
    if (filespace < 0)
    {
        std::cerr
            << "parallelWrite2d(), H5Screate_simple failed for filespace!!!"
            << std::endl;
        return;
    }
    hid_t memspace = H5Screate_simple(2, dimsm, nullptr);
    if (memspace < 0)
    {
        std::cerr
            << "parallelWrite2d(), H5Screate_simple failed for memspace!!!"
            << std::endl;
        return;
    }

    // Create dataset
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, dimsm);
    hid_t dset_id = H5Dcreate2(file_id, datasetname.c_str(), strtype, filespace,
        H5P_DEFAULT, plist_id, H5P_DEFAULT);
    if (dset_id < 0)
    {
        std::cerr << "parallelWrite2d() for dataset " << datasetname
                  << ", H5Dcreate2() failed!!!" << std::endl;
        return;
    }
    H5Pclose(plist_id);

    hsize_t offset[2] = { mpi_rank * dimsm[0], 0 };
    hsize_t stride[2] = { 1, 1 };
    hsize_t count[2]  = { 1, 1 };
    hsize_t block[2]  = { dimsm[0], dimsm[1] };

    /* Select hyperslab in the file. */
    herr_t status = H5Sselect_hyperslab(
        filespace, H5S_SELECT_SET, offset, stride, count, block);
    if (status < 0)
    {
        std::cerr << "parallelWrite2d(), H5Sselect_hyperslab() failed!!!"
                  << std::endl;
        return;
    }

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // First copy the contents of the vector into a temporary container
    std::vector<FixedLengthString> tc;
    string2fixedlength(data, tc);
    status = H5Dwrite(dset_id, strtype, memspace, filespace, plist_id, &tc[0]);
    if (status < 0)
    {
        std::cerr << "parallelWrite2d(), H5Dwrite failed!!!" << std::endl;
        return;
    }

    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
}

#endif

void addAttribute2Dataset(
    hid_t dset_id, const char* attname, const std::vector<double>& attr_data)
{
    assert(dset_id > -1);

    // Create the data space for the attribute.
    hsize_t dim = (hsize_t)attr_data.size();

    //  Open a dataset attribute.
    hid_t dataspace_id = H5Screate_simple(1, &dim, nullptr);
    if (dataspace_id < 0)
    {
        std::cerr << "H5Screate failed!!!" << std::endl;
        return;
    }
    hid_t attribute_id = H5Acreate2(dset_id, attname, H5T_NATIVE_DOUBLE,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attribute_id < 0)
    {
        std::cerr << "H5Acreate failed!!!" << std::endl;
        return;
    }

    herr_t status = H5Sclose(dataspace_id);
    if (status < 0) std::cerr << "H5Sclose failed!!!" << std::endl;

    //(*MPIdata::sout)<<"Write attribute "<<attname<<endl;
    status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attr_data[0]);
    if (status < 0) std::cerr << "H5Awrite failed!!!" << std::endl;

    status = H5Aclose(attribute_id);
    if (status < 0) std::cerr << "H5Aclose failed!!!" << std::endl;
}

void addAttribute2Dataset(
    hid_t dset_id, const char* attname, const std::vector<int>& attr_data)
{
    assert(dset_id > -1);

    // Create the data space for the attribute.
    hsize_t dim = (hsize_t)attr_data.size();

    //  Open a dataset attribute.
    hid_t dataspace_id = H5Screate_simple(1, &dim, nullptr);
    if (dataspace_id < 0)
    {
        std::cerr << "H5Screate failed!!!" << std::endl;
        return;
    }
    hid_t attribute_id = H5Acreate2(dset_id, attname, H5T_NATIVE_INT,
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attribute_id < 0)
    {
        std::cerr << "H5Acreate failed!!!" << std::endl;
        return;
    }

    herr_t status = H5Sclose(dataspace_id);
    if (status < 0)
    {
        std::cerr << "H5Sclose failed!!!" << std::endl;
        return;
    }

    //(*MPIdata::sout)<<"Write attribute "<<attname<<endl;
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data[0]);
    if (status < 0)
    {
        std::cerr << "H5Awrite failed!!!" << std::endl;
        return;
    }

    status = H5Aclose(attribute_id);
    if (status < 0)
    {
        std::cerr << "H5Aclose failed!!!" << std::endl;
    }
}

// adapted from https://www.hdfgroup.org/HDF5-FAQ.html
int whatisopen(hid_t fid)
{
    char name[1024];

    ssize_t cnt = H5Fget_obj_count(fid, H5F_OBJ_ALL);

    if (cnt <= 0) return cnt;

    if (cnt > 1) std::cout << "HDF5 file: " << cnt << " object(s) open\n";

    // objs = malloc(cnt * sizeof(hid_t));
    hid_t* objs = new hid_t[cnt];

    int howmany = H5Fget_obj_ids(fid, H5F_OBJ_ALL, cnt, objs);

    if (cnt > 1) printf("open objects:\n");

    hid_t* obj = objs;
    for (int i = 0; i < howmany; i++)
    {
        hid_t anobj   = *obj++;
        H5I_type_t ot = H5Iget_type(anobj);
        H5Iget_name(anobj, name, 1024);
        if (ot != H5I_FILE)
            printf("HDF object %d: type %d, name %s\n", i, ot, name);
    }

    delete[] objs;

    return howmany;
}

int write_matrix(
    hid_t file_id, std::string& name, const double* matrix, const int dim)
{
    if (file_id < 0) return 0;

    hsize_t dims[2] = { (hsize_t)dim, (hsize_t)dim };

    // filespace identifier
    hid_t dataspace = H5Screate_simple(2, dims, nullptr);

    hid_t dset_id = H5Dcreate2(file_id, name.c_str(), H5T_NATIVE_DOUBLE,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset_id < 0)
    {
        MGMOL_HDF5_FAIL("H5Dcreate2 failed!!!");
        return -1;
    }

    hid_t memspace  = dataspace;
    hid_t filespace = dataspace;

    herr_t status = H5Dwrite(
        dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, matrix);
    if (status < 0)
    {
        MGMOL_HDF5_FAIL("H5Dwrite failed!!!");
        return -1;
    }

    status = H5Dclose(dset_id);
    if (status < 0)
    {
        MGMOL_HDF5_FAIL("H5Dclose failed!!!");
        return -1;
    }
    status = H5Sclose(dataspace);
    if (status < 0)
    {
        MGMOL_HDF5_FAIL("H5Sclose failed!!!");
        return -1;
    }

    return 0;
}

int read_matrix(hid_t file_id, std::string& name, double* matrix)
{
    int ierr      = 0;
    hid_t dset_id = H5Dopen2(file_id, name.c_str(), H5P_DEFAULT);
    if (dset_id < 0)
    {
        MGMOL_HDF5_FAIL("H5Dopen failed!!");
        ierr = -1;
    }
    else
    {
        herr_t status = H5Dread(
            dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix);
        if (status < 0)
        {
            MGMOL_HDF5_FAIL("H5Dread failed!!");
            ierr = -1;
        }

        status = H5Dclose(dset_id);
        if (status < 0)
        {
            MGMOL_HDF5_FAIL("H5Dclose failed!!!");
            ierr = -1;
        }
    }

    return ierr;
}

} // namespace
