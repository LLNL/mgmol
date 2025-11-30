// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SubspaceProjector.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "ProjectedMatricesInterface.h"

#define NPRINT_ROWS_AND_COLS 10

template <class T>
SubspaceProjector<T>::SubspaceProjector(T& subspace)
    : subspace_(subspace), proj_matrices_(*subspace_.getProjMatrices())
{
    chromatic_number_ = subspace_.chromatic_number();
    subdivx_          = subspace_.subdivx();
    lda_              = subspace_.getLda();
    loc_numpt_        = subspace_.getLocNumpt();
}

// compute [I-P*(S^-1)*P^T]*orbitals
template <class T>
void SubspaceProjector<T>::projectOut(
    T& orbitals, SquareLocalMatrices<MATDTYPE, MemorySpace::Host>* mask)
{
    assert(chromatic_number_ >= 0);
    assert(lda_ >= loc_numpt_);

    // assumption for now, to be removed in the future
    assert(chromatic_number_ == orbitals.chromatic_number());
#if 0
    {
    // test if projection is now 0
    dist_matrix::DistMatrix<DISTMATDTYPE> tmatrix(subspace_.product(orbitals));
    if( onpe0 )
        (*MPIdata::sout)<<"SubspaceProjector: Product before projection"<<endl;
    tmatrix.print((*MPIdata::sout),0,0,NPRINT_ROWS_AND_COLS,NPRINT_ROWS_AND_COLS);
    }
#endif

    SquareLocalMatrices<MATDTYPE, MemorySpace::Host> pmatrix(
        subdivx_, chromatic_number_);

    // compute <xi,orbitals>
    if (chromatic_number_ != 0)
        subspace_.computeLocalProduct(orbitals, pmatrix, false);

    // compute <xi,xi>^{-1}*<xi,orbitals>
    proj_matrices_.applyInvS(pmatrix);

    if (mask != nullptr) pmatrix.applyMask(*mask);

    ORBDTYPE* tproduct = new ORBDTYPE[loc_numpt_ * chromatic_number_];
    memset(tproduct, 0, loc_numpt_ * chromatic_number_ * sizeof(ORBDTYPE));

    // loop over subdomains
    for (short iloc = 0; iloc < subdivx_; iloc++)
    {
        ORBDTYPE* xi            = subspace_.getPsi(0, iloc);
        MATDTYPE* localMat_iloc = pmatrix.getRawPtr(iloc);

        // Compute loc_numpt_ rows (for subdomain iloc)
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(loc_numpt_,
            chromatic_number_, chromatic_number_, 1., xi, lda_, localMat_iloc,
            chromatic_number_, 0., tproduct, loc_numpt_);

        double minus     = -1.;
        ORBDTYPE* parray = orbitals.getPsi(0, iloc);
        for (int j = 0; j < chromatic_number_; j++)
            LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(loc_numpt_, minus,
                tproduct + j * loc_numpt_, parray + j * lda_);
    }

    delete[] tproduct;

    orbitals.incrementIterativeIndex();
}

template class SubspaceProjector<LocGridOrbitals<ORBDTYPE>>;
template class SubspaceProjector<ExtendedGridOrbitals<ORBDTYPE>>;
