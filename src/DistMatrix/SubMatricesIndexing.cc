// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SubMatricesIndexing.h"
#include "MPI_AllToSome.h"

#include <iostream>
#include <limits.h>

namespace dist_matrix
{

// constructor based on vectors of indexes corresponding
// to "local" matrix elements
template <class T>
SubMatricesIndexing<T>::SubMatricesIndexing(
    const std::vector<std::vector<int>>&
        indexes, // functions ids for each subdomain on PE
    MPI_Comm comm, const DistMatrix<T>& empty_mat)
    : comm_(comm), local_indexes_(indexes)
{
    MPI_Comm_size(comm_, &npes_);
    MPI_Comm_rank(comm_, &mype_);

    npes_distmat_ = empty_mat.nprow() * empty_mat.npcol();
    // npes_distmat_=npes_;
    assert(npes_distmat_ > 0);

#ifdef DEBUG
    if (mype_ == 0)
        cout << "SubMatrices with " << indexes[0].size() << " indexes"
             << std::endl;
#endif

    nb_local_matrices_ = (int)local_indexes_.size();
    const int nlocal   = (int)local_indexes_[0].size();
    assert(nlocal < SHRT_MAX);
    assert(nlocal <= empty_mat.m());
    assert(nlocal <= empty_mat.n());
    assert(nlocal >= 0);
    assert(nb_local_matrices_ > 0);
    for (int im = 0; im < nb_local_matrices_; im++)
        assert((int)local_indexes_[im].size() == nlocal);

    // mapping between index and position in local matrix
    map_indexes_.resize(nb_local_matrices_);
    for (int im = 0; im < nb_local_matrices_; im++)
    {
        const std::vector<int>& local_indexes_im = local_indexes_[im];
        // cout<<"Matrix "<<im;
        // for(int j=0;j<nlocal;j++)
        //  cout<<" index "<<local_indexes_im[j];
        // cout<<endl;
        for (short color = 0; color < nlocal; color++)
            map_indexes_[im].insert(
                std::pair<int, short>(local_indexes_im[color], color));
        assert((int)map_indexes_[im].size() <= nlocal);
    }

    // initialize arrays of all indexes (i,j) needed from other PEs
    const int npcol = empty_mat.npcol();

    // for each active pe (for DistMatrix), store local associated functions
    // indexes
    std::vector<std::vector<int>> pe_double_local_indexes;
    pe_double_local_indexes.resize(npes_distmat_);
    for (int i = 0; i < npes_distmat_; i++)
        pe_double_local_indexes[i].reserve(32);

    for (int im = 0; im < nb_local_matrices_; im++)
    {
        assert((int)local_indexes_[im].size() == nlocal);

        std::vector<int>::const_iterator p0 = local_indexes_[im].begin();
        while (p0 != local_indexes_[im].end())
        {
            const int i = (*p0);
            assert(i < empty_mat.m());
            if (i >= 0)
            {
                const int pr = empty_mat.pr(i);

                std::vector<int>::const_iterator p1
                    = local_indexes_[im].begin();
                while (p1 != local_indexes_[im].end())
                {
                    const int j = (*p1);
                    if (j >= 0)
                    {
                        const int pc = empty_mat.pc(j);

                        // index (i,j) value should be received from task
                        // pe_distmat
                        const int pe_distmat = pc + pr * npcol;
                        assert(pe_distmat >= 0);
                        assert(pe_distmat < npes_distmat_);
                        // cout<<"(i,j)=("<<i<<","<<j<<"),
                        // pe_distmat="<<pe_distmat<<endl;
                        pe_double_local_indexes[pe_distmat].push_back(i);
                        pe_double_local_indexes[pe_distmat].push_back(j);
                    }
                    p1++;
                }
            }
            p0++;
        }
    }

    // store pe_double_local_indexes into single array for communications
    double_local_indexes_.reserve(256);

    my_sizes_.resize(npes_distmat_, 0);
    for (int pe_distmat = 0; pe_distmat < npes_distmat_; pe_distmat++)
    {
        const std::vector<int>& pe_double_local_indexes_pe
            = pe_double_local_indexes[pe_distmat];
        my_sizes_[pe_distmat] = (int)pe_double_local_indexes_pe.size();
        // cout<<"mype="<<mype_<<",
        // my_sizes_[pe="<<pe_distmat<<"]="<<my_sizes_[pe_distmat]<<endl;
        assert(my_sizes_[pe_distmat] % 2 == 0);
        assert(my_sizes_[pe_distmat]
               <= 2 * nb_local_matrices_ * empty_mat.m() * empty_mat.n());
        for (int ip = 0; ip < my_sizes_[pe_distmat]; ip += 2)
        {
            assert(pe_double_local_indexes_pe[ip] < empty_mat.m());
            assert(pe_double_local_indexes_pe[ip + 1] < empty_mat.n());
            double_local_indexes_.push_back(pe_double_local_indexes_pe[ip]);
            double_local_indexes_.push_back(pe_double_local_indexes_pe[ip + 1]);
        }
    }

    pe_double_local_indexes.clear(); // don't need it anymore, uses memory!

    int my_size = 0;
    for (int p = 0; p < npes_distmat_; p++)
        my_size += my_sizes_[p];

    // tell all the (DistMatrix)active PEs how many indexes I will need from
    // them
    // if( mype_<npes_distmat_ )
    remote_sizes_.resize(npes_, 0);
#if 0
   MPI_Alltoall(&my_sizes_[0],1,MPI_INT,&remote_sizes_[0],1,MPI_INT,comm_);
#else
    MPI_AlltofirstN(&my_sizes_[0], 1, &remote_sizes_[0], npes_distmat_, comm_);
#endif

#ifdef DEBUG
    cout << "mype=" << mype_ << ", my_sizes_[0]=" << my_sizes_[0] << std::endl;
#endif

    int tot_remote_size = 0;
    if (mype_ < npes_distmat_)
        for (int p = 0; p < npes_; p++)
            tot_remote_size += remote_sizes_[p];

    // check data to send
    for (int i = 0; i < (int)double_local_indexes_.size(); i += 2)
    {
        // cout<<"mype="<<mype_<<", my_index="
        //    <<double_local_indexes_[i]<<","<<double_local_indexes_[i+1]<<endl;
        assert(double_local_indexes_[i] < empty_mat.m());
        assert(double_local_indexes_[i + 1] < empty_mat.n());
    }

    // tell remote PEs which indexes I need to exchange
    if (mype_ < npes_distmat_)
        remote_double_indexes_.resize(tot_remote_size, -1);
    std::vector<int> my_displ(npes_distmat_, 0);
    for (int pe = 1; pe < npes_distmat_; pe++)
        my_displ[pe] = my_sizes_[pe - 1] + my_displ[pe - 1];
    std::vector<int> remote_displ;
    if (mype_ < npes_distmat_)
    {
        remote_displ.resize(npes_, 0);
        for (int pe = 1; pe < npes_; pe++)
            remote_displ[pe] = remote_sizes_[pe - 1] + remote_displ[pe - 1];
    }
#if 0
   MPI_Alltoallv(&double_local_indexes_[0],    &my_sizes_[0],    &my_displ[0],    MPI_INT,
                 &remote_double_indexes_[0],&remote_sizes_[0],&remote_displ[0],MPI_INT,
                 comm_);
#else
    MPI_AlltofirstNv(&double_local_indexes_[0], &my_sizes_[0], &my_displ[0],
        &remote_double_indexes_[0], &remote_sizes_[0], &remote_displ[0],
        npes_distmat_, comm_);
#endif
    remote_displ.clear();

    // check received data
    for (int i = 0; i < (int)remote_double_indexes_.size(); i += 2)
    {
        // cout<<"mype="<<mype_<<", remote_index="
        //    <<remote_double_indexes_[i]<<","<<remote_double_indexes_[i+1]<<endl;
        assert(remote_double_indexes_[i] < empty_mat.m());
        assert(remote_double_indexes_[i + 1] < empty_mat.n());
    }

    // adapt sizes and to deal with double value
    // instead of indexes (2 int)
    for (int p = 0; p < npes_distmat_; p++)
    {
        my_sizes_[p] /= 2;
        assert(
            my_sizes_[p] <= nb_local_matrices_ * empty_mat.m() * empty_mat.n());
    }
    if (mype_ < npes_distmat_)
        for (int p = 0; p < npes_; p++)
        {
            remote_sizes_[p] /= 2;
            assert(remote_sizes_[p]
                   <= nb_local_matrices_ * empty_mat.m() * empty_mat.n());
        }
#ifdef DEBUG
    std::cout << "mype=" << mype_ << ", remote_sizes_[0]=" << remote_sizes_[0]
              << std::endl;
#endif

    // precompute mapping between index and position in local matrix
    vector_indexes_.resize(nb_local_matrices_);
    const int maxii = (int)double_local_indexes_.size() / 2;
    for (int im = 0; im < nb_local_matrices_; im++)
    {
        std::map<int, short>& mapind = map_indexes_[im];
        vector_indexes_[im].resize(nlocal * nlocal, -1);
        // loop over indexes (2 int)
        // index ii and ii+1
        for (int ii = 0; ii < maxii; ii++)
        {
            std::map<int, short>::const_iterator p0
                = mapind.find(double_local_indexes_[2 * ii]);
            if (p0 != mapind.end())
            {
                std::map<int, short>::const_iterator p1
                    = mapind.find(double_local_indexes_[2 * ii + 1]);
                if (p1 != mapind.end())
                {
                    short colori = p0->second;
                    short colorj = p1->second;
                    assert(colori >= 0);
                    assert(colorj >= 0);
                    // cout<<"PE:"<<mype_<<", index
                    // "<<double_local_indexes_[ii]<<","<<double_local_indexes_[ii+1]<<endl;
                    // cout<<"PE:"<<mype_<<", (i,j)=("<<i<<","<<j<<")"<<endl;
                    // cout<<"PE:"<<mype_<<",
                    // buf_my_val["<<ii/2<<"]="<<buf_my_val[ii/2]<<endl;
                    assert(static_cast<unsigned int>(colori + nlocal * colorj)
                           < vector_indexes_[im].size());
                    vector_indexes_[im][colori + nlocal * colorj] = ii;
                }
            }
        }
    }
}

template class SubMatricesIndexing<DISTMATDTYPE>;

} // namespace
