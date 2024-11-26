// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "hdf5.h"

#include <mpi.h>
#include <string>
#include <vector>

namespace mgmol_tools
{
void write1d(hid_t file_id, const std::string& datasetname,
    std::vector<int>& data, size_t length);
void write2d(hid_t file_id, const std::string& datasetname,
    std::vector<int>& data, size_t* dims);
void write2d(hid_t file_id, const std::string& datasetname,
    std::vector<unsigned short>& data, size_t* dims);
void write2d(hid_t file_id, const std::string& datasetname,
    std::vector<double>& data, size_t* dims);
void write2d(hid_t file_id, const std::string& datasetname,
    std::vector<std::string>& data, size_t* dims);

void parallelWrite2d(hid_t file_id, const std::string& datasetname,
    std::vector<int>& data, size_t* dims, MPI_Comm comm);
void parallelWrite2d(hid_t file_id, const std::string& datasetname,
    std::vector<unsigned short>& data, size_t* dims, MPI_Comm comm);
void parallelWrite2d(hid_t file_id, const std::string& datasetname,
    std::vector<double>& data, size_t* dims, MPI_Comm comm);
void parallelWrite2d(hid_t file_id, const std::string& datasetname,
    std::vector<std::string>& data, size_t* dims, MPI_Comm comm);

void addAttribute2Dataset(
    hid_t dset_id, const char* attname, const std::vector<double>& attr_data);
void addAttribute2Dataset(
    hid_t dset_id, const char* attname, const std::vector<int>& attr_data);
int whatisopen(hid_t fid);

int write_matrix(
    hid_t file_id, std::string& name, const double* matrix, const int dim);
int read_matrix(hid_t file_id, std::string& name, double* matrix);
}
