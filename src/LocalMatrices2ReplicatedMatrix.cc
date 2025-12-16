// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "LocalMatrices2ReplicatedMatrix.h"
#include "MGmol_MPI.h"

LocalMatrices2ReplicatedMatrix* LocalMatrices2ReplicatedMatrix::pinstance_
    = nullptr;
std::vector<std::vector<int>> LocalMatrices2ReplicatedMatrix::global_indexes_;
double LocalMatrices2ReplicatedMatrix::tol_mat_elements = 1.e-14;

Timer LocalMatrices2ReplicatedMatrix::convert_tm_(
    "LocalMatrices2ReplicatedMatrix::convert");

void LocalMatrices2ReplicatedMatrix::convert(
    const LocalMatrices<double, MemorySpace::Host>& src, ReplicatedMatrix& dst,
    const int numst, const double tol) const
{
    (void)tol;

    assert(!global_indexes_.empty());

    convert_tm_.start();

    const int subdiv = static_cast<int>(global_indexes_.size());

    std::vector<double> val(subdiv);

    const short chromatic_number
        = static_cast<short>(global_indexes_[0].size());

    std::vector<double> data(numst * numst);

    // double loop over colors
    for (short icolor = 0; icolor < chromatic_number; icolor++)
    {
        for (short jcolor = 0; jcolor < chromatic_number; jcolor++)
        {
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
                        const int idst             = st2 * numst + st1;
                        const double* const ssiloc = src.getSubMatrix(iloc);
                        const double tmp
                            = ssiloc[jcolor * chromatic_number + icolor];

                        // accumulate values
                        data[idst] += tmp;
                    }
                }

            } // iloc

        } // jcolor
    } // icolor

    dst.assign(data.data(), numst);

    convert_tm_.stop();
}

// Sum up all the local contributions (in LocalMatrices) into
// one ReplicatedMatrix
void LocalMatrices2ReplicatedMatrix::accumulate(
    const LocalMatrices<double, MemorySpace::Host>& src, ReplicatedMatrix& dst,
    const double tol) const
{
    // convert into a ReplicatedMatrix
    convert(src, dst, dst.m(), tol);

    // accumulate into a ReplicatedMatrix
    dst.consolidate();
}
