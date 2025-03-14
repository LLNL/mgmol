// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef SYMMETRICMATRIX_H
#define SYMMETRICMATRIX_H

#include <cassert>
#include <cstring>
#include <map>
#include <mpi.h>
#include <stdio.h>
#include <vector>

// full lower triangular part of symmetric matrix with compact storage

template <class T>
class SymmetricMatrix
{
private:
    int dimension_; // matrix dimension
    int size_;
    T* data_;

    // global index number for rows and columns
    std::vector<int> gids_;

    std::map<int, int> gid2lid_;

    // MPI
    const MPI_Comm comm_;

    int offset(const int i) const { return (i * (i + 1) >> 1); }

public:
    SymmetricMatrix(const int dimension, const MPI_Comm comm) : comm_(comm)
    {
        dimension_ = dimension;
        size_      = ((dimension * (dimension + 1)) >> 1);
        data_      = new T[size_];
        memset(data_, 0, size_ * sizeof(T));

        gids_.reserve(64);
        for (int i = 0; i < dimension_; i++)
        {
            gids_.push_back(i);
        }

        for (int i = 0; i < dimension_; i++)
        {
            gid2lid_.insert(std::pair<int, int>(i, i));
        }
    }

    SymmetricMatrix(const std::vector<int>& gids, const MPI_Comm comm)
        : gids_(gids), comm_(comm)
    {
        dimension_ = gids_.size();
        size_      = ((dimension_ * (dimension_ + 1)) >> 1);
        data_      = new T[size_];
        memset(data_, 0, size_ * sizeof(T));

        int i = 0;
        for (std::vector<int>::const_iterator it = gids_.begin();
             it != gids_.end(); it++)
        {
            gid2lid_.insert(std::pair<int, int>(*it, i));
            i++;
        }
    }

    ~SymmetricMatrix() { delete[] data_; }

    T val(const int gid1, const int gid2) const
    {
        std::map<int, int>::const_iterator it1 = gid2lid_.find(gid1);
        std::map<int, int>::const_iterator it2 = gid2lid_.find(gid2);

        assert(it1 != gid2lid_.end());
        assert(it2 != gid2lid_.end());
        assert(it1->second < dimension_);
        assert(it2->second < dimension_);

        const int i = it1->second;
        const int j = it2->second;

        const int ii = i > j ? i : j;
        const int jj = i > j ? j : i;
        return data_[offset(ii) + jj];
    }

    void setVal(const int gid1, const int gid2, const T val)
    {
        std::map<int, int>::const_iterator it1 = gid2lid_.find(gid1);
        std::map<int, int>::const_iterator it2 = gid2lid_.find(gid2);

        assert(it1 != gid2lid_.end());
        assert(it2 != gid2lid_.end());
        assert(it1->second < dimension_);
        assert(it2->second < dimension_);

        const int i = it1->second;
        const int j = it2->second;

        const int ii           = i > j ? i : j;
        const int jj           = i > j ? j : i;
        data_[offset(ii) + jj] = val;
    }
    void mpiAllOr();

    int size() const { return size_; }
    int dimension() const { return dimension_; }

    const std::vector<int>& gids() const { return gids_; }

    /* print non-lval pattern in coordinate format */
    void printSparsePattern(int lval, char* fname) const
    {
        FILE* fout = fopen(fname, "w");
        for (int i = 0; i < dimension_; i++)
        {
            for (int j = 0; j < dimension_; j++)
            {
                if (val(i, j) == (T)lval) continue;
                fprintf(fout, "%d %d %d \n", i, j, val(i, j));
            }
        }
    }

    /* convert from dense matrix structure to csr matrix structure */
    void getNNZPattern(
        std::vector<int>& ia, std::vector<int>& ja, int* const nnzrow) const
    {
        ia.push_back(0);
        int nzmax = 0;

        /* get data for first row and reserve memory */
        nnzrow[0] = 0;
        for (int j = 0; j < dimension_; j++)
        {
            if (val(0, j) == 0) continue;
            ja.push_back(j);
            nnzrow[0]++;
        }
        ia.push_back(nnzrow[0]);
        /* done with first row. Reserve memory storage */
        nzmax = (nnzrow[0] + 1)
                * dimension_; // assumes matrix has similar nonzero row count
        ja.reserve(nzmax);

        /* gather remaining matrix data */
        for (int i = 1; i < dimension_; i++)
        {
            nnzrow[i] = 0;
            for (int j = 0; j < dimension_; j++)
            {
                if (val(i, j) == 0) continue;
                ja.push_back(j);
                nnzrow[i]++;
            }
            ia.push_back(ia[i] + nnzrow[i]);
        }
    }
};

#endif
