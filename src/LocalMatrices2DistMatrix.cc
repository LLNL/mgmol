// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "LocalMatrices2DistMatrix.h"
#include "MGmol_MPI.h"

LocalMatrices2DistMatrix* LocalMatrices2DistMatrix::pinstance_ = nullptr;
MPI_Comm LocalMatrices2DistMatrix::comm_   = MPI_COMM_NULL;
std::vector<std::vector<int>> LocalMatrices2DistMatrix::global_indexes_;
double LocalMatrices2DistMatrix::tol_mat_elements = 1.e-14;
dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>*
    LocalMatrices2DistMatrix::remote_tasks_DistMatrix_ = nullptr;
short LocalMatrices2DistMatrix::sparse_distmatrix_nb_tasks_per_partitions_ = 256;

template <class T>
void LocalMatrices2DistMatrix::convert(
    const LocalMatrices<T>& lmat,
    dist_matrix::SparseDistMatrix<T>& dst,
    const int numst,
    const double tol) const
{
    const int subdiv = (int)global_indexes_.size();

    int* ist = new int[2 * subdiv];
    T* val   = new T[subdiv];

    const short chromatic_number = (short)global_indexes_[0].size();

    // double loop over colors
    for (short icolor = 0; icolor < chromatic_number; icolor++)
    {
        for (short jcolor = 0; jcolor < chromatic_number; jcolor++)
        {
            int pst_old  = -1;
            int valindex = -1;

            // loop over subdomains
            for (short iloc = 0; iloc < subdiv; iloc++)
            {
                const int st1 = global_indexes_[iloc][icolor];
                //(*MPIdata::sout)<<"icolor="<<icolor<<", pst1="<<pst1<<endl;

                if (st1 != -1)
                {
                    const int st2 = global_indexes_[iloc][jcolor];
                    if (st2 != -1)
                    {
                        // unique id for current pair
                        const int pst         = st2 * numst + st1;
                        const T* const ssiloc = lmat.getSubMatrix(iloc);
                        const T tmp
                            = ssiloc[jcolor * chromatic_number + icolor];

                        // accumulate values
                        if (fabs(tmp) > tol)
                        {
                            if (pst == pst_old)
                            {
                                assert(valindex >= 0);
                                val[valindex] += tmp;
                            }
                            else
                            {
                                valindex++;
                                val[valindex]         = tmp;
                                ist[2 * valindex]     = st1;
                                ist[2 * valindex + 1] = st2;
                                pst_old               = pst;
                            }
                        }
                    }
                }

            } // iloc

            // assign values
            for (int i = 0; i <= valindex; i++)
                dst.push_back(ist[2 * i], ist[2 * i + 1], (double)val[i]);

        } // jcolor
    } // icolor

    delete[] val;
    delete[] ist;
}

template<class T>
void LocalMatrices2DistMatrix::convert(const LocalMatrices<T>& lmat,
        dist_matrix::DistMatrix<T>& dst,
        const int numst,
        const double tol)const
{
    assert(remote_tasks_DistMatrix_ != nullptr);
#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sm(comm, dst,
        remote_tasks_DistMatrix_, sparse_distmatrix_nb_tasks_per_partitions_);
#else
    dist_matrix::SparseDistMatrix<DISTMATDTYPE> sm(
        0, dst, remote_tasks_DistMatrix_);
#endif

    convert(lmat, sm, numst, tol);

    sm.parallelSumToDistMatrix();
}

template void LocalMatrices2DistMatrix::convert(
    const LocalMatrices<DISTMATDTYPE>& lmat, dist_matrix::SparseDistMatrix<DISTMATDTYPE>& dst,
    const int numst, const double tol)const;
template void LocalMatrices2DistMatrix::convert(
    const LocalMatrices<DISTMATDTYPE>& lmat, dist_matrix::DistMatrix<DISTMATDTYPE>& dst,
    const int numst, const double tol)const;
