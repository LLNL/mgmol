// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef DistributedIonicData_H
#define DistributedIonicData_H

#include "IonData.h"

#include <iterator>
#include <mpi.h>
#include <string>
#include <vector>

class DistributedIonicData
{
    friend class Ions;
    friend class Ion;

private:
    std::vector<std::string> ion_names_;
    std::vector<double> data_;

    int source_[3];
    int dest_[3];
    MPI_Comm cart_comm_; // MPI cartesian communicator for data distribution

    // unpack data in buffers and add it to local data
    void unpack(char*& cptr, double*& iptr, const short ndata);

    // pack local data into a char buffer and a double buffer
    int pack(char* cbuff, double* ibuff);

public:
    DistributedIonicData(
        const std::vector<std::string>&, const std::vector<double>&);
    DistributedIonicData() {};

    int size() const { return (int)ion_names_.size(); };

    // get data associated with name "ion_name"
    void getData(const std::string& ion_name, double* data)
    {
        short i = 0;
        for (std::vector<std::string>::const_iterator it = ion_names_.begin();
             it != ion_names_.end(); ++it)
        {
            if (it->compare(ion_name) == 0)
            {
                data[0] = data_[3 * i + 0];
                data[1] = data_[3 * i + 1];
                data[2] = data_[3 * i + 2];

                break;
            }
            i++;
        }
    }

    // augment local data by communicating in one direction
    void augmentData(const int nsteps, const int dir, const int disp,
        const int locSize, const int maxLocSize,
        DistributedIonicData& data2send);
};

#endif
