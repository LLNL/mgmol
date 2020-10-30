// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Rho.h"
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"
#include "SquareLocalMatrices.h"
#include "SubMatrices.h"
#include "magma_singleton.h"
#include "memory_space.h"
#include "mputils.h"
#include "numerical_kernels.h"

template <class OrbitalsType>
Timer Rho<OrbitalsType>::update_tm_("Rho::update");
template <class OrbitalsType>
Timer Rho<OrbitalsType>::compute_tm_("Rho::compute");
template <class OrbitalsType>
Timer Rho<OrbitalsType>::compute_blas_tm_("Rho::compute_usingBlas");
// template <class OrbitalsType>
// Timer Rho<OrbitalsType>::compute_offdiag_tm_("Rho::compute_offdiag");

#ifdef HAVE_MAGMA
template <typename ScalarType>
using MemoryDev = MemorySpace::Memory<ScalarType, MemorySpace::Device>;
#endif

template <class OrbitalsType>
Rho<OrbitalsType>::Rho()
    : orbitals_type_(OrthoType::UNDEFINED),
      iterative_index_(-10),
      verbosity_level_(0) // default value
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    myspin_         = mmpi.myspin();
    nspin_          = mmpi.nspin();

    // default values for block sizes
    //    block_functions_ = 8;
    //    block_space_     = 256;

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    np_                    = mygrid.size();

    rho_.resize(nspin_);
    for (short is = 0; is < nspin_; is++)
    {
        rho_[is].resize(np_);
        memset(&rho_[is][0], 0, np_ * sizeof(RHODTYPE));
    }

    assert(nspin_ == 1 || nspin_ == 2);
    assert(myspin_ == 0 || myspin_ == 1);
}

template <class OrbitalsType>
void Rho<OrbitalsType>::setup(const OrthoType orbitals_type,
    const std::vector<std::vector<int>>& orbitals_indexes)
{
    if (verbosity_level_ > 2 && onpe0)
        (*MPIdata::sout) << " Rho<OrbitalsType>::setup()" << std::endl;

    orbitals_type_ = orbitals_type;

    orbitals_indexes_ = orbitals_indexes;
}

template <class OrbitalsType>
void Rho<OrbitalsType>::extrapolate()
{
    double minus = -1;
    double two   = 2.;
    if (rho_minus1_.empty())
    {
        rho_minus1_.resize(nspin_);
        rho_minus1_[myspin_].resize(np_);
        memcpy(&rho_minus1_[myspin_][0], &rho_[myspin_][0],
            np_ * sizeof(RHODTYPE));
        return;
    }
    RHODTYPE* tmp = new RHODTYPE[np_];
    memcpy(tmp, &rho_[myspin_][0], np_ * sizeof(RHODTYPE));

    LinearAlgebraUtils<MemorySpace::Host>::MPscal(np_, two, &rho_[myspin_][0]);
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        np_, minus, &rho_minus1_[myspin_][0], &rho_[myspin_][0]);

    memcpy(&rho_minus1_[myspin_][0], tmp, np_ * sizeof(RHODTYPE));

    delete[] tmp;
}

template <class OrbitalsType>
void Rho<OrbitalsType>::axpyRhoc(const double alpha, RHODTYPE* rhoc)
{
    double factor = (nspin_ > 1) ? 0.5 * alpha : alpha;
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        np_, factor, &rhoc[0], &rho_[myspin_][0]);
}

template <class OrbitalsType>
void Rho<OrbitalsType>::update(OrbitalsType& current_orbitals)
{
    const ProjectedMatricesInterface& proj_matrices(
        *(current_orbitals.getProjMatrices()));

    assert(current_orbitals.getIterativeIndex() >= 0);
    assert(proj_matrices.getDMMatrixIndex() >= 0);

    update_tm_.start();

    if (verbosity_level_ > 2 && onpe0)
        (*MPIdata::sout) << "Rho<OrbitalsType>::update()" << std::endl;

    const int new_iterative_index
        = ((1 + current_orbitals.getIterativeIndex()) % 100)
          + (proj_matrices.getDMMatrixIndex() % 100) * 100;

    if (iterative_index_ == new_iterative_index)
    {
        if (onpe0 && verbosity_level_ > 2)
            (*MPIdata::sout) << "Rho already up to date, iterative_index_="
                             << iterative_index_ << std::endl;
        return;
    }
    iterative_index_ = new_iterative_index;
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "Rho<OrbitalsType>::update(), iterative_index_="
                         << iterative_index_ << std::endl;
#endif

    computeRho(current_orbitals);

    rescaleTotalCharge();

    assert(iterative_index_ >= 0);

    update_tm_.stop();
}

// note: rho can be negative because of added background charge
template <class OrbitalsType>
double Rho<OrbitalsType>::computeTotalCharge()
{
    const int nspin = (int)rho_.size();

    double tcharge = 0.;
    for (short ispin = 0; ispin < nspin; ispin++)
    {
        const RHODTYPE* const prho = &rho_[ispin][0];
        for (int idx = 0; idx < np_; idx++)
        {
            assert(prho[idx] < 1.e6);
            tcharge += (double)prho[idx];
        }
    }

    // reduce over spin communicator since charge from other spin is already
    // included
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&tcharge, 1, MPI_SUM);

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    tcharge *= mygrid.vel();
#ifdef DEBUG
    if (mmpi.instancePE0())
        (*MPIdata::sout) << fixed << setprecision(8) << " myspin=" << myspin_
                         << ", Total charge = " << tcharge << std::endl;
#endif

    return tcharge;
}

template <class OrbitalsType>
void Rho<OrbitalsType>::rescaleTotalCharge()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    // Check total charge
    const double tcharge = computeTotalCharge();

    const int nspin = (int)rho_.size();

    Control& ct = *(Control::instance());
    double nel  = ct.getNel();
    if (tcharge > 0.)
    {
        double t1 = nel / tcharge;

        if (mmpi.PE0() && fabs(t1 - 1.) > 0.001)
        {
            (*MPIdata::sout)
                << " Rho<OrbitalsType>::rescaleTotalCharge(), charge = "
                << tcharge << std::endl;
            (*MPIdata::sout) << " Rescaling factor: " << t1 << std::endl;
            (*MPIdata::sout) << " Num. electrons: " << nel << std::endl;
        }

        for (int ispin = 0; ispin < nspin; ispin++)
            LinearAlgebraUtils<MemorySpace::Host>::MPscal(
                np_, t1, &rho_[ispin][0]);
    }
#ifdef DEBUG
    Mesh* mymesh = Mesh::instance();
    // Check total charge again
    double ttcharge = 0.;
    for (int ispin = 0; ispin < nspin; ispin++)
        for (int idx = 0; idx < np_; idx++)
        {
            ttcharge += (double)rho_[ispin][idx];
        }

    mmpi.allreduce(&ttcharge, 1, MPI_SUM);
    const pb::Grid& mygrid = mymesh->grid();
    ttcharge *= mygrid.vel();

    if (mmpi.instancePE0())
        (*MPIdata::sout) << std::fixed << std::setprecision(8)
                         << " Total charge after rescaling: " << ttcharge
                         << std::endl;
#endif
}

// template <class OrbitalsType>
// int Rho<OrbitalsType>::setupSubdomainData(const int iloc,
//    const vector<const T*>& vorbitals,
//    const ProjectedMatricesInterface* const projmatrices,
//    vector<MATDTYPE>& melements, vector<vector<const ORBDTYPE*>>& vmpsi)
//{
//    // printWithTimeStamp("Rho<OrbitalsType>::setupSubdomainData()...",cout);
//
//    const short norb = (short)vorbitals.size();
//    vmpsi.resize(norb);
//
//    const vector<int>& loc_indexes(orbitals_indexes_[iloc]);
//    const int n_colors = (int)loc_indexes.size();
//
//    if (n_colors == 0) return 0;
//
//    vector<int> mycolors;
//    mycolors.reserve(n_colors);
//    for (int icolor = 0; icolor < n_colors; icolor++)
//        if (loc_indexes[icolor] != -1) mycolors.push_back(icolor);
//    const int nmycolors = (int)mycolors.size();
//
//    for (short j = 0; j < norb; j++)
//        vmpsi[j].resize(nmycolors);
//
//    short j = 0;
//    for (typename vector<const T*>::const_iterator it = vorbitals.begin();
//         it != vorbitals.end(); ++it)
//    {
//        for (int color = 0; color < nmycolors; color++)
//        {
//            vmpsi[j][color] = (*it)->getPsi(mycolors[color], iloc);
//            assert(vmpsi[j][color] != NULL);
//        }
//        j++;
//    }
//
//    SquareLocalMatrices<MATDTYPE>& localX(projmatrices->getLocalX());
//    const MATDTYPE* const localX_iloc = localX.getSubMatrix(iloc);
//    melements.clear();
//    melements.resize(nmycolors * nmycolors);
//
//    for (int i = 0; i < nmycolors; i++)
//    {
//        const int icolor = mycolors[i];
//        if (norb == 1)
//        {
//            melements[i * nmycolors + i]
//                = localX_iloc[icolor + n_colors * icolor];
//        }
//        const int jmax = (norb == 1) ? i : nmycolors;
//        for (int j = 0; j < jmax; j++)
//        {
//            int jcolor = mycolors[j];
//            melements[j * nmycolors + i]
//                = localX_iloc[icolor + n_colors * jcolor];
//        }
//    }
//
//    return nmycolors;
//}

template <class OrbitalsType>
void Rho<OrbitalsType>::accumulateCharge(const double alpha, const short ix_max,
    const ORBDTYPE* const psii, const ORBDTYPE* const psij,
    RHODTYPE* const plrho)
{
    for (int ix = 0; ix < ix_max; ix++)
        plrho[ix] += (RHODTYPE)(alpha * (double)psii[ix] * (double)psij[ix]);
}

// template <class OrbitalsType>
// void Rho<OrbitalsType>::computeRhoSubdomain(
//    const int iloc_init, const int iloc_end, const OrbitalsType& orbitals)
//{
//    assert(orbitals_type_ == OrthoType::Eigenfunctions
//        || orbitals_type_ == OrthoType::Nonorthogonal);
//
//    compute_tm_.start();
//
//    Mesh* mymesh = Mesh::instance();
//
//    const int loc_numpt = mymesh->locNumpt();
//
//    RHODTYPE* const prho = &rho_[myspin_][0];
//
//    vector<const T*> vorbitals;
//    vorbitals.push_back(&orbitals);
//
//    const double nondiagfactor = 2.;
//
//    for (int iloc = iloc_init; iloc < iloc_end; iloc++)
//    {
//        const int istart = iloc * loc_numpt;
//
//        RHODTYPE* const lrho = &prho[istart];
//
//        vector<MATDTYPE> melements;
//        vector<vector<const ORBDTYPE*>> vmpsi;
//
//        const int nmycolors = setupSubdomainData(
//            iloc, vorbitals, orbitals.projMatrices(), melements, vmpsi);
//        assert(vmpsi.size() == 1);
//        vector<const ORBDTYPE*> mpsi = vmpsi[0];
//
//        const int nblocks_color = nmycolors / block_functions_;
//        const int max_icolor    = (nmycolors % block_functions_ == 0)
//                                   ? nmycolors
//                                   : nblocks_color * block_functions_;
//        const int missed_rows = nmycolors - max_icolor;
//        //(*MPIdata::sout)<<"max_icolor="<<max_icolor<<endl;
//        //(*MPIdata::sout)<<"Rho<OrbitalsType>::computeRhoSubdomain:
//        //missed_rows="<<missed_rows<<endl;
//        //(*MPIdata::sout)<<"nblocks_color="<<nblocks_color<<endl;
//
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//        for (int idx = 0; idx < loc_numpt; idx += block_space_)
//        {
//            const short ix_max = min(block_space_, loc_numpt - idx);
//            // RHODTYPE* const plrho=lrho+idx;
//
//            // non-diagonal blocks
//            for (int icolor = 0; icolor < max_icolor;
//                 icolor += block_functions_)
//                for (int jcolor = 0; jcolor < icolor;
//                     jcolor += block_functions_)
//                {
//
//                    nonOrthoRhoKernel(icolor, block_functions_, jcolor,
//                        block_functions_, idx, ix_max, &melements[0],
//                        nmycolors, mpsi, nondiagfactor, lrho);
//                }
//            // finish missing rows ()
//            if (missed_rows)
//                for (int jcolor = 0; jcolor < max_icolor;
//                     jcolor += block_functions_)
//                {
//                    nonOrthoRhoKernel(max_icolor, missed_rows, jcolor,
//                        block_functions_, idx, ix_max, &melements[0],
//                        nmycolors, mpsi, nondiagfactor, lrho);
//                }
//
//            // diagonal blocks (jcolor=icolor)
//            for (int icolor = 0; icolor < max_icolor;
//                 icolor += block_functions_)
//            {
//                nonOrthoRhoKernelDiagonalBlock(icolor, block_functions_, idx,
//                    ix_max, &melements[0], nmycolors, mpsi, lrho);
//            }
//            // finish missing rows () for diagonal blocks
//            if (missed_rows)
//            {
//                nonOrthoRhoKernelDiagonalBlock(max_icolor, missed_rows, idx,
//                    ix_max, &melements[0], nmycolors, mpsi, lrho);
//            }
//        }
//    }
//
//    compute_tm_.stop();
//}

// template <class OrbitalsType>
// void Rho<OrbitalsType>::computeRhoSubdomainOffDiagBlock(const int iloc_init,
//    const int iloc_end, const vector<const T*>& vorbitals,
//    const ProjectedMatricesInterface* const projmatrices)
//{
//    assert(orbitals_type_ == OrthoType::Eigenfunctions
//        || orbitals_type_ == OrthoType::Nonorthogonal);
//    assert(vorbitals.size() == 2);
//
//    compute_offdiag_tm_.start();
//
//    //
//    printWithTimeStamp("Rho<OrbitalsType>::computeRhoSubdomainOffDiagBlock()...",cout);
//
//    Mesh* mymesh = Mesh::instance();
//
//    const int loc_numpt = mymesh->locNumpt();
//
//    RHODTYPE* const prho = &rho_[myspin_][0];
//
//    for (int iloc = iloc_init; iloc < iloc_end; iloc++)
//    {
//        const int istart = iloc * loc_numpt;
//
//        RHODTYPE* const lrho = &prho[istart];
//
//        vector<MATDTYPE> melements;
//        vector<vector<const ORBDTYPE*>> vmpsi;
//
//        const int nmycolors = setupSubdomainData(
//            iloc, vorbitals, projmatrices, melements, vmpsi);
//        const int nblocks_color = nmycolors / block_functions_;
//        const int max_icolor    = (nmycolors % block_functions_ == 0)
//                                   ? nmycolors
//                                   : nblocks_color * block_functions_;
//        const int missed_rows = nmycolors - max_icolor;
//        //(*MPIdata::sout)<<"max_icolor="<<max_icolor<<endl;
//        //(*MPIdata::sout)<<"Rho<OrbitalsType>::computeRhoSubdomainOffDiagBlock:
//        //missed_rows="<<missed_rows<<endl;
//        //(*MPIdata::sout)<<"nblocks_color="<<nblocks_color<<endl;
//
//        for (int idx = 0; idx < loc_numpt; idx += block_space_)
//        {
//            const short ix_max    = min(block_space_, loc_numpt - idx);
//            RHODTYPE* const plrho = lrho + idx;
//
//            for (int istart = 0; istart < max_icolor;
//                 istart += block_functions_)
//                for (int jstart = 0; jstart < max_icolor;
//                     jstart += block_functions_)
//                {
//                    for (short jcolor = jstart;
//                         jcolor < jstart + block_functions_; jcolor++)
//                    {
//                        const int jld              = jcolor * nmycolors;
//                        const ORBDTYPE* const psij = &vmpsi[1][jcolor][idx];
//                        for (short icolor = istart;
//                             icolor < istart + block_functions_; icolor++)
//                        {
//                            const ORBDTYPE* const psii =
//                            &vmpsi[0][icolor][idx]; const double alpha
//                                = (double)melements[jld + icolor];
//                            accumulateCharge(alpha, ix_max, psii, psij,
//                            plrho);
//                        }
//                    }
//                }
//            // finish missing rows and cols
//            if (missed_rows)
//            {
//                int icolor = max_icolor;
//                int imax   = nmycolors % block_functions_;
//                for (int jcolor = 0; jcolor < max_icolor;
//                     jcolor += block_functions_)
//                {
//                    for (short j = 0; j < block_functions_; j++)
//                    {
//                        const int jstart = (jcolor + j) * nmycolors + icolor;
//                        const ORBDTYPE* const psij = &vmpsi[1][jcolor +
//                        j][idx]; for (short i = 0; i < imax; i++)
//                        {
//                            const ORBDTYPE* const psii
//                                = &vmpsi[0][icolor + i][idx];
//                            const double alpha = (double)melements[jstart +
//                            i]; accumulateCharge(alpha, ix_max, psii, psij,
//                            plrho);
//                        }
//                    }
//                }
//            }
//            if (missed_rows)
//            {
//                int jcolor = max_icolor;
//                int jmax   = nmycolors % block_functions_;
//                for (int icolor = 0; icolor < max_icolor;
//                     icolor += block_functions_)
//                {
//                    for (short j = 0; j < jmax; j++)
//                    {
//                        const int jstart = (jcolor + j) * nmycolors + icolor;
//                        const ORBDTYPE* const psij = &vmpsi[1][jcolor +
//                        j][idx]; for (short i = 0; i < block_functions_; i++)
//                        {
//                            const ORBDTYPE* const psii
//                                = &vmpsi[0][icolor + i][idx];
//                            const double alpha = (double)melements[jstart +
//                            i]; accumulateCharge(alpha, ix_max, psii, psij,
//                            plrho);
//                        }
//                    }
//                }
//            }
//            if (missed_rows)
//            {
//                int icolor = max_icolor;
//                int imax   = nmycolors % block_functions_;
//                int jcolor = max_icolor;
//                int jmax   = nmycolors % block_functions_;
//                for (short j = 0; j < jmax; j++)
//                {
//                    const int jstart = (jcolor + j) * nmycolors + icolor;
//                    const ORBDTYPE* const psij = &vmpsi[1][jcolor + j][idx];
//                    for (short i = 0; i < imax; i++)
//                    {
//                        const ORBDTYPE* const psii = &vmpsi[0][icolor +
//                        i][idx]; const double alpha = (double)melements[jstart
//                        + i]; accumulateCharge(alpha, ix_max, psii, psij,
//                        plrho);
//                    }
//                }
//            }
//        }
//    }
//
//    compute_offdiag_tm_.stop();
//}

template <class OrbitalsType>
void Rho<OrbitalsType>::computeRhoSubdomain(const int iloc_init,
    const int iloc_end, const OrbitalsType& orbitals,
    const std::vector<double>& occ)
{
    assert(orbitals_type_ != OrthoType::UNDEFINED);
    if (verbosity_level_ > 2 && onpe0)
        (*MPIdata::sout)
            << "Rho<OrbitalsType>::computeRhoSubdomain, diagonal case..."
            << std::endl;

    Mesh* mymesh = Mesh::instance();

    const int loc_numpt = mymesh->locNumpt();
    const int n_colors  = orbitals_indexes_[0].size();

    RHODTYPE* const prho = &rho_[myspin_][0];

    for (int iloc = iloc_init; iloc < iloc_end; iloc++)
    {
        const int istart              = iloc * loc_numpt;
        std::vector<int>& loc_indexes = orbitals_indexes_[iloc];

        RHODTYPE* const lrho = prho + istart;

        // Loop over states and accumulate charge
        for (int icolor = 0; icolor < n_colors; icolor++)
        {
            const ORBDTYPE* const psi = orbitals.getPsi(icolor, iloc);
            const int st              = loc_indexes[icolor];
            assert(st < (int)occ.size());
            if (st >= 0)
            {
                const double t1 = 2. * (double)occ[st];
                for (int idx = 0; idx < loc_numpt; idx++)
                {
                    const double alpha = (double)psi[idx];
                    lrho[idx] += (RHODTYPE)(t1 * alpha * alpha);
                }
            }
        }
    }
}

template <class OrbitalsType>
void Rho<OrbitalsType>::computeRho(OrbitalsType& orbitals)
{
    ProjectedMatricesInterface& proj_matrices(*(orbitals.getProjMatrices()));

    computeRho(orbitals, proj_matrices);
}

template <class OrbitalsType>
void Rho<OrbitalsType>::computeRho(
    OrbitalsType& orbitals, ProjectedMatricesInterface& proj_matrices)
{
    assert(rho_.size() > 0);
    assert(rho_[myspin_].size() > 0);

    Mesh* mymesh        = Mesh::instance();
    const int subdivx   = mymesh->subdivx();
    const int loc_numpt = mymesh->locNumpt();
    Control& ct         = *(Control::instance());

    memset(&rho_[myspin_][0], 0, subdivx * loc_numpt * sizeof(RHODTYPE));

    if (orbitals_type_ == OrthoType::Eigenfunctions
        || (orbitals_type_ == OrthoType::Orthonormal && ct.fullyOccupied()))
    {
        std::vector<double> occ(orbitals.numst());
        proj_matrices.getOccupations(occ);
        computeRhoSubdomain(0, subdivx, orbitals, occ);
    }
    else
    {
        proj_matrices.updateSubMatX();

        // if (dynamic_cast<LocGridOrbitals*>(&orbitals)) but it
        if (std::is_same<OrbitalsType, LocGridOrbitals>::value)
        {
            SquareLocalMatrices<MATDTYPE>& localX(
                (orbitals.projMatrices())->getLocalX());

            // if (verbosity_level_ > 1 && onpe0)
            //    (*MPIdata::sout) << "Mask DM..." << endl;
            localX.applySymmetricMask(orbitals_indexes_);
        }
        computeRhoSubdomainUsingBlas3(0, subdivx, orbitals);
    }

    gatherSpin();
}

template <class OrbitalsType>
template <class MatrixType>
void Rho<OrbitalsType>::computeRho(OrbitalsType& orbitals1,
    OrbitalsType& orbitals2, const MatrixType& dm11, const MatrixType& dm12,
    const MatrixType& /*dm21*/, const MatrixType& dm22)
{
    assert(orbitals_type_ == OrthoType::Nonorthogonal);

    Mesh* mymesh        = Mesh::instance();
    const int subdivx   = mymesh->subdivx();
    const int loc_numpt = mymesh->locNumpt();

    memset(&rho_[myspin_][0], 0, subdivx * loc_numpt * sizeof(RHODTYPE));

    // 11 diagonal block
    ProjectedMatrices<MatrixType>* projmatrices1
        = dynamic_cast<ProjectedMatrices<MatrixType>*>(
            orbitals1.getProjMatrices());
    projmatrices1->updateSubMatX(dm11);
    computeRhoSubdomainUsingBlas3(0, subdivx, orbitals1);

    // 22 diagonal block
    ProjectedMatrices<MatrixType>* projmatrices2
        = dynamic_cast<ProjectedMatrices<MatrixType>*>(
            orbitals2.getProjMatrices());
    projmatrices2->updateSubMatX(dm22);
    computeRhoSubdomainUsingBlas3(0, subdivx, orbitals2);

    // non diagonal blocks
    // vector<const T*> vorbitals;
    // vorbitals.push_back(&orbitals1);
    // vorbitals.push_back(&orbitals2);

    MatrixType dm(dm12);
    dm.scal(2.); // use symmetry to reduce work
    projmatrices1->updateSubMatX(dm);
    //    computeRhoSubdomainOffDiagBlock(0, subdivx, vorbitals, projmatrices1);
    computeRhoSubdomainUsingBlas3(0, subdivx, orbitals1, orbitals2);

    gatherSpin();
}

template <class OrbitalsType>
void Rho<OrbitalsType>::computeRhoSubdomainUsingBlas3(const int iloc_init,
    const int iloc_end, const OrbitalsType& orbitals1,
    const OrbitalsType& orbitals2)
{
    assert(orbitals1.getLda() == orbitals2.getLda());
    assert(orbitals1.chromatic_number() == orbitals2.chromatic_number());

    compute_blas_tm_.start();

    Mesh* mymesh    = Mesh::instance();
    const int nrows = mymesh->locNumpt();

    RHODTYPE* const prho = &rho_[myspin_][0];

    SquareLocalMatrices<MATDTYPE>& localX(
        (orbitals1.projMatrices())->getLocalX());

    const int ld      = orbitals1.getLda();
    const int ncols   = orbitals1.chromatic_number();
    ORBDTYPE* product = new ORBDTYPE[nrows * ncols];

    for (int iloc = iloc_init; iloc < iloc_end; iloc++)
    {
        const int istart = iloc * nrows;

        RHODTYPE* const lrho = &prho[istart];

        const MATDTYPE* const mat = localX.getSubMatrix(iloc);
        ORBDTYPE* phi1            = orbitals1.getPsi(0, iloc);

        // O(N^3) part
#ifdef HAVE_MAGMA
        // If we have magma, we move all the data on the device
        std::unique_ptr<MATDTYPE[], void (*)(MATDTYPE*)> mat_dev(
            MemoryDev<MATDTYPE>::allocate(ncols * ncols),
            MemoryDev<MATDTYPE>::free);
        MemorySpace::copy_to_dev(mat, ncols * ncols, mat_dev);

        std::unique_ptr<ORBDTYPE[], void (*)(ORBDTYPE*)> product_dev(
            MemoryDev<ORBDTYPE>::allocate(nrows * ncols),
            MemoryDev<ORBDTYPE>::free);

        LinearAlgebraUtils<MemorySpace::Device>::MPgemmNN(nrows, ncols, ncols,
            1., phi1, ld, mat_dev.get(), ncols, 0., product_dev.get(), nrows);

        // Move the data back on the host
        MemorySpace::copy_to_host(product_dev, nrows * ncols, product);
#else // MAGMA is not installed
        LinearAlgebraUtils<MemorySpace::Host>::MPgemmNN(
            nrows, ncols, ncols, 1., phi1, ld, mat, ncols, 0., product, nrows);
#endif

        // O(N^2) part: if we have MAGMA and we can offload with OpenMP, we
        // perform the operation on the device. Otherwise, we do it on the host.
#ifdef HAVE_OPENMP_OFFLOAD
        ORBDTYPE* product_alias = product_dev.get();
        ORBDTYPE* phi2_alias    = orbitals2.getPsi(0, iloc);
        // Copy lrho to the device
        std::unique_ptr<RHODTYPE[], void (*)(RHODTYPE*)> lrho_dev(
            MemoryDev<RHODTYPE>::allocate(nrows), MemoryDev<RHODTYPE>::free);
        MemorySpace::copy_to_dev(lrho, nrows, lrho_dev);
        RHODTYPE* const lrho_alias = lrho_dev.get();
#else
        ORBDTYPE* product_alias = product;
        ORBDTYPE* phi2          = orbitals2.getPsi(0, iloc);
        using memory_space_type = typename OrbitalsType::memory_space_type;
        ORBDTYPE* phi2_alias    = MemorySpace::Memory<ORBDTYPE,
            memory_space_type>::allocate_host_view(ld * ncols + nrows);
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::copy_view_to_host(
            phi2, ld * ncols + nrows, phi2_alias);
        RHODTYPE* const lrho_alias = lrho;
#endif
        MGMOL_PARALLEL_FOR_COLLAPSE(2, lrho_alias, product_alias, phi2_alias)
        for (int j = 0; j < ncols; ++j)
        {
            for (int i = 0; i < nrows; ++i)
            {
#pragma omp atomic update
                lrho_alias[i]
                    += product_alias[j * nrows + i] * phi2_alias[j * ld + i];
            }
        }
#if defined(HAVE_MAGMA) && defined(HAVE_OPENMP_OFFLOAD)
        // Move lrho back to the host
        MemorySpace::copy_to_host(lrho_alias, nrows, lrho);
#else
        MemorySpace::Memory<ORBDTYPE, memory_space_type>::free_host_view(
            phi2_alias);
#endif
    }

    delete[] product;

    compute_blas_tm_.stop();
}

template <class OrbitalsType>
template <class MatrixType>
void Rho<OrbitalsType>::computeRho(OrbitalsType& orbitals, const MatrixType& dm)
{
    assert(orbitals_type_ == OrthoType::Nonorthogonal);

    iterative_index_++;

    Mesh* mymesh        = Mesh::instance();
    const int subdivx   = mymesh->subdivx();
    const int loc_numpt = mymesh->locNumpt();

    memset(&rho_[myspin_][0], 0, subdivx * loc_numpt * sizeof(RHODTYPE));

    ProjectedMatrices<MatrixType>* projmatrices
        = dynamic_cast<ProjectedMatrices<MatrixType>*>(
            orbitals.getProjMatrices());
    projmatrices->updateSubMatX(dm);

    computeRhoSubdomainUsingBlas3(0, subdivx, orbitals);

    gatherSpin();
}

template <class OrbitalsType>
void Rho<OrbitalsType>::init(const RHODTYPE* const rhoc)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int ione        = 1;

    // Initialize the charge density
    if (verbosity_level_ > 2 && mmpi.instancePE0())
        (*MPIdata::sout) << "Rho: Initialize electronic density value with rhoc"
                         << std::endl;
    for (unsigned int i = 0; i < rho_.size(); i++)
    {
        Tcopy(&np_, rhoc, &ione, &rho_[i][0], &ione);
        if (rho_.size() == 2)
            LinearAlgebraUtils<MemorySpace::Host>::MPscal(
                np_, 0.5, &rho_[i][0]);
    }
    iterative_index_ = 0;

    rescaleTotalCharge();
}

template <class OrbitalsType>
void Rho<OrbitalsType>::initUniform()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    // Initialize the charge density
    if (mmpi.instancePE0())
        (*MPIdata::sout) << " Initialize electronic density with uniform value"
                         << std::endl;
    for (unsigned int i = 0; i < rho_.size(); i++)
    {
        for (int j = 0; j < np_; j++)
            rho_[i][j] = 1.;
        if (rho_.size() == 2)
            LinearAlgebraUtils<MemorySpace::Host>::MPscal(
                np_, 0.5, &rho_[i][0]);
    }
    iterative_index_ = 0;

    rescaleTotalCharge();
}

// read rho and potentials form a hdf5 file
template <class OrbitalsType>
int Rho<OrbitalsType>::readRestart(HDFrestart& file)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Try to read density" << std::endl;

    // Read the Density
    file.read_1func_hdf5(&rho_[myspin_][0], "Density");
    iterative_index_ = 0;

    return 0;
}

template <class OrbitalsType>
template <typename T2>
double Rho<OrbitalsType>::dotWithRho(const T2* const func) const
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    double val = LinearAlgebraUtils<MemorySpace::Host>::MPdot(
        np_, &rho_[mmpi.myspin()][0], func);

    double esum = 0.;
    mmpi.allreduce(&val, &esum, 1, MPI_SUM);
    val = esum;

    mmpi.allreduceSpin(&val, &esum, 1, MPI_SUM);
    val = esum;

    return val;
}

template <class OrbitalsType>
void Rho<OrbitalsType>::gatherSpin()
{
    if (nspin_ < 2) return;

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    // if(mmpi.instancePE0())cout<<"Rho<OrbitalsType>::gatherSpin(),
    // myspin_="<<myspin_<<endl;
    mmpi.exchangeDataSpin(&rho_[myspin_][0], &rho_[(myspin_ + 1) % 2][0], np_);
}

template <class OrbitalsType>
void Rho<OrbitalsType>::printTimers(std::ostream& os)
{
    update_tm_.print(os);
    compute_tm_.print(os);
    compute_blas_tm_.print(os);
    // compute_offdiag_tm_.print(os);
}

template class Rho<LocGridOrbitals>;
template class Rho<ExtendedGridOrbitals>;

template double Rho<LocGridOrbitals>::dotWithRho<double>(
    const double* const func) const;
template double Rho<ExtendedGridOrbitals>::dotWithRho<double>(
    const double* const func) const;
template void
Rho<ExtendedGridOrbitals>::computeRho<dist_matrix::DistMatrix<double>>(
    ExtendedGridOrbitals&, ExtendedGridOrbitals&,
    const dist_matrix::DistMatrix<double>&,
    const dist_matrix::DistMatrix<double>&,
    const dist_matrix::DistMatrix<double>&,
    const dist_matrix::DistMatrix<double>&);
template void
Rho<ExtendedGridOrbitals>::computeRho<dist_matrix::DistMatrix<double>>(
    ExtendedGridOrbitals&, const dist_matrix::DistMatrix<double>&);
template void Rho<LocGridOrbitals>::computeRho<dist_matrix::DistMatrix<double>>(
    LocGridOrbitals&, const dist_matrix::DistMatrix<double>&);
#ifdef USE_MP
template double Rho<LocGridOrbitals>::dotWithRho<float>(
    const float* const func) const;
#endif
