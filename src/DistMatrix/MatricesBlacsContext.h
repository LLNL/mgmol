// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef MatricesBlacsContext_H
#define MatricesBlacsContext_H

#include "BlacsContext.h"
#include <cassert>
#include <iostream>

#include <mpi.h>

class MatricesBlacsContext
{
    MatricesBlacsContext()
        : blactxt_(nullptr), comm_(MPI_COMM_NULL), size_(-1), n_(-1)
    {
    }

    ~MatricesBlacsContext()
    {
        if (blactxt_ != nullptr) delete blactxt_;
    };

    dist_matrix::BlacsContext* blactxt_;
    MPI_Comm comm_;
    int size_;
    int n_;

public:
    static MatricesBlacsContext& instance();

    void setup(const MPI_Comm comm, const int size);

    dist_matrix::BlacsContext* bcxt() const
    {
        assert(blactxt_ != 0);
        return blactxt_;
    }

    void print(std::ostream& os) const
    {
        os << "MatricesBlacsContext: " << n_ << " x " << n_ << std::endl;
    }

    void clear()
    {
        if (blactxt_ != nullptr)
        {
            delete blactxt_;
            blactxt_ = nullptr;
        }
        comm_ = MPI_COMM_NULL;
        size_ = -1;
        n_    = -1;
    }
};

#endif
