// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "KBPsiMatrixSparse.h"

#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol_MPI.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"
#include "ReplicatedMatrix.h"

#ifdef MGMOL_USE_SCALAPACK
#include "SquareSubMatrix2DistMatrix.h"
#endif

#include <limits.h>

#define Ry2Ha 0.5;

Timer KBPsiMatrixSparse::global_sum_tm_("KBPsiMatrixSparse::global_sum");
Timer KBPsiMatrixSparse::compute_kbpsi_tm_("KBPsiMatrixSparse::compute_kbpsi");
Timer KBPsiMatrixSparse::computeHvnlMatrix_tm_(
    "KBPsiMatrixSparse::computeHvnlMatrix");
Timer KBPsiMatrixSparse::setup_tm_("KBPsiMatrixSparse::setup");
Timer KBPsiMatrixSparse::trace_tm_("KBPsiMatrixSparse::trace");

static const double tolKBpsi = 1.e-12;

KBPsiMatrixSparse::KBPsiMatrixSparse(
    pb::Lap<ORBDTYPE>* lapop, const bool need2radius)
    : lapop_(lapop), need2radius_(need2radius)
{
    kbpsimat_    = nullptr;
    kbBpsimat_   = nullptr;
    distributor_ = nullptr;

    isDataSetup_ = false;
}

KBPsiMatrixSparse::~KBPsiMatrixSparse() { clearData(); }

void KBPsiMatrixSparse::clearData()
{
    if (isDataSetup_)
    {
        delete kbpsimat_;
        kbpsimat_ = nullptr;
        if (lapop_)
        {
            delete kbBpsimat_;
        }
        kbBpsimat_ = nullptr;
        delete distributor_;
        distributor_ = nullptr;

        isDataSetup_ = false;
    }
}

void KBPsiMatrixSparse::setup(const Ions& ions)
{
    setup_tm_.start();

    setOutdated();

    // clear old data
    clearData();

    // number of projectors overlapping with subdomain
    count_proj_subdomain_ = ions.countProjectorsSubdomain();
    /* setup matrices */
    assert(kbpsimat_ == nullptr);
    assert(kbBpsimat_ == nullptr);
    kbpsimat_
        = new VariableSizeMatrix<sparserow>("KBpsi", count_proj_subdomain_);
    if (lapop_)
        kbBpsimat_ = new VariableSizeMatrix<sparserow>(
            "KBBpsi", count_proj_subdomain_);
    else
        kbBpsimat_ = kbpsimat_;

    spread_radius_ = ions.getMaxVnlRadius();

    // we need to know locally all the info corresponding to projectors
    // overlapping with local subdomain
    if (need2radius_) spread_radius_ *= 2.0;

    //    if( onpe0 )
    //        (*MPIdata::sout)<<"KBPsiMatrixSparse using radius
    //        "<<spread_radius_<<endl;

    /* construct data distribution object */
    assert(distributor_ == nullptr);
    Mesh* mymesh             = Mesh::instance();
    const pb::Grid& mygrid   = mymesh->grid();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double domain[3]         = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    distributor_ = new DataDistribution("KB", spread_radius_, myPEenv, domain);

    isDataSetup_ = true;

    setup_tm_.stop();
}

// consolidate values with overlap of local projectors on other subdomains
void KBPsiMatrixSparse::globalSumKBpsi()
{
    bool append = false; // assuming projectors have rectangular domain
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "KBPsiMatrixSparse::globalSumKBpsi()" << endl;
#endif

    global_sum_tm_.start();

    /* perform data distribution */
    (*distributor_).augmentLocalData((*kbpsimat_), append);
    if (lapop_)
    {
        (*distributor_).augmentLocalData((*kbBpsimat_), append);
    }
    global_sum_tm_.stop();

    return;
}

// Loop over the ions with projectors overlapping with local subdomain
// and evaluate <KB|psi> for some state.
template <class T>
void KBPsiMatrixSparse::computeKBpsi(const Ions& ions, T& orbitals,
    const int first_color, const int nb_colors, const bool flag)
{
    assert(first_color >= 0);
    assert(first_color < 100000);
    assert(nb_colors > 0);
    assert(nb_colors < 100000);
    assert(orbitals.getIterativeIndex() >= 0);

    compute_kbpsi_tm_.start();

    ORBDTYPE* ppsi   = nullptr;
    const int ldsize = orbitals.getLda();
    assert(ldsize > 0);
    assert(ldsize < 1e8);

    if (flag)
    {
        assert(lapop_ != nullptr);
        ppsi = new ORBDTYPE[nb_colors * ldsize];
        orbitals.setDataWithGhosts();
        orbitals.trade_boundaries();
#pragma omp parallel for
        for (int color = 0; color < nb_colors; color++)
        {
            lapop_->rhs(orbitals.getFuncWithGhosts(first_color + color),
                ppsi + color * ldsize);
        }
    }
    else
    {
        ORBDTYPE* psi_first_color    = orbitals.getPsi(first_color);
        unsigned int const ppsi_size = nb_colors * ldsize;
        ppsi                         = MemorySpace::Memory<ORBDTYPE,
            typename T::memory_space_type>::allocate_host_view(ppsi_size);

        MemorySpace::Memory<ORBDTYPE,
            typename T::memory_space_type>::copy_view_to_host(psi_first_color,
            ppsi_size, ppsi);
    }

    // Loop over functions, subdomains and ions
    Mesh* mymesh      = Mesh::instance();
    const int subdivx = mymesh->subdivx();
    for (int iloc = 0; iloc < subdivx; iloc++)
    {
        // Threading here leads to results slightly dependent on number of
        // threads (jlf, 07/15/2016)
        // #pragma omp parallel for
        for (int color = 0; color < nb_colors; color++)
        {
            const int gid = orbitals.getGlobalIndex(iloc, first_color + color);

            // Loop over the ions
            if (gid != -1)
            {
                for (const auto& ion : ions.overlappingNL_ions())
                {
                    computeLocalElement(
                        *ion, gid, iloc, ppsi + color * ldsize, flag);
                }
            }
        }
    }

    if (flag)
    {
        delete[] ppsi;
    }
    else
    {
        MemorySpace::Memory<ORBDTYPE,
            typename T::memory_space_type>::free_host_view(ppsi);
    }

    compute_kbpsi_tm_.stop();
}

void KBPsiMatrixSparse::computeKBpsi(const Ions& ions,
    pb::GridFunc<ORBDTYPE>* phi, const int istate, const bool flag)
{
    assert(lapop_ != nullptr);
    compute_kbpsi_tm_.start();

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const int ldsize = mygrid.size();

    ORBDTYPE* ppsi = new ORBDTYPE[ldsize];

    if (flag)
    {
        lapop_->rhs(*phi, ppsi);
    }
    else
    {
        phi->init_vect(ppsi, 'd');
    }

    // Loop over states, subdomains and ions
    const int subdivx = mymesh->subdivx();
    for (int iloc = 0; iloc < subdivx; iloc++)
    {
        // Loop over the ions
        for (const auto& ion : ions.overlappingNL_ions())
        {
            computeLocalElement(*ion, istate, iloc, ppsi, flag);
        }
    }

    delete[] ppsi;

    compute_kbpsi_tm_.stop();
}

void KBPsiMatrixSparse::scaleWithKBcoeff(const Ions& ions)
{
    for (const auto& ion : ions.overlappingNL_ions())
    {
        std::vector<int> gids;
        ion->getGidsNLprojs(gids);
        std::vector<double> kbcoeffs;
        ion->getKBcoeffs(kbcoeffs);

        const short nprojs = (short)gids.size();
        for (short i = 0; i < nprojs; i++)
        {
            const int gid = gids[i];
            double coeff  = kbcoeffs[i];

            // loop over states to multiply kbpsi_[st][gid] and kbBpsi_[st][gid]
            // by coeff
            kbpsimat_->scaleRow(gid, coeff);
            if (lapop_) kbBpsimat_->scaleRow(gid, coeff);
        }
    }
}

// computeHvnlMatrix:
//    Get the elements of the Hamiltonian matrix due to the non-local
//    potential, and add them into Aij.
// Note: neglecting the small matrix elements reduces the size of hnlij and thus
//       reduces the size of communications later on.
void KBPsiMatrixSparse::computeHvnlElementsIon(
    const KBPsiMatrixSparse* const kbpsi2, const Ion& ion,
    SquareSubMatrix<double>& hnlij) const
{
    assert(ion.here());

    // get projectors ids
    std::vector<int> pgids;
    ion.getGidsNLprojs(pgids);

    const short nprojs = (short)pgids.size();
    if (nprojs == 0) return;

    std::vector<short> kbsigns;
    ion.getKBsigns(kbsigns);

    VariableSizeMatrix<sparserow>* kbBpsimat2 = kbpsi2->kbpsimat_;

    // loop over projectors associated with ion

    int* rindex0 = (int*)(kbBpsimat_->getTableValue(pgids[0]));
    assert(rindex0 != nullptr);
    const int lrindex0 = *rindex0;

    // get number of functions that overlaps with projector
    const int nnzrow = kbBpsimat_->nnzrow(lrindex0);
    assert(kbBpsimat2->nnzrow(lrindex0) == nnzrow);

    std::vector<int> gids1;
    kbBpsimat_->getColumnIndexes(lrindex0, gids1);

    std::vector<int> gids2;
    kbBpsimat2->getColumnIndexes(lrindex0, gids2);

    for (short i = 0; i < nprojs; i++)
    {
        const int gid = pgids[i];

        int* rindex = (int*)(kbBpsimat_->getTableValue(gid));
        assert(rindex != nullptr);

        // double loop over states to fill hnlij[st1][st2] (in general not
        // symmetric... )
        const double coeff = (double)kbsigns[i];
        const int lrindex  = *rindex;

        // assume all rows associated with this atoms have the same
        // number of non-zero elements
        assert(nnzrow == kbBpsimat_->nnzrow(lrindex));

        for (int p1 = 0; p1 < nnzrow; p1++)
        {
            const double kbpsi1 = coeff * kbBpsimat_->getRowEntry(lrindex, p1);
            if (fabs(kbpsi1) > tolKBpsi)
            {
                for (int p2 = 0; p2 < nnzrow; p2++)
                {
                    const double alpha
                        = kbpsi1 * kbBpsimat2->getRowEntry(lrindex, p2);
                    // set hnlij
                    if (fabs(alpha) > tolKBpsi)
                    {
                        // hnlij.push_back(gids1[p1], gids2[p2], alpha);
                        hnlij.addValue(gids1[p1], gids2[p2], alpha);
                    }
                }
            }
        }
    }
}

void KBPsiMatrixSparse::computeHvnlElementsIon(
    const KBPsiMatrixSparse* const kbpsi2, const Ion& ion,
    VariableSizeMatrix<sparserow>& mat) const
{
    assert(ion.here());

    std::vector<int> gids;
    ion.getGidsNLprojs(gids);

    std::vector<short> coeffs;
    ion.getKBsigns(coeffs);

    const short nprojs = (short)gids.size();
    for (short i = 0; i < nprojs; i++)
    {
        const int gid      = gids[i];
        const double coeff = (double)coeffs[i];

        // double loop over states to fill hnlij[st1][st2] (in general not
        // symmetric... )
        int* rindex = (int*)(*kbBpsimat_).getTableValue(gid);
        if (rindex == nullptr) continue;
        const int lrindex = *rindex;
        const int nnzrow1 = kbBpsimat_->nnzrow(lrindex);
        for (int p1 = 0; p1 < nnzrow1; p1++)
        {
            const int st1 = kbBpsimat_->getColumnIndex(lrindex, p1);
            const double kbpsi1
                = coeff * (*kbBpsimat_).getRowEntry(lrindex, p1);
            if (fabs(kbpsi1) > tolKBpsi)
            {
                const int nnzrow2 = (*kbpsi2->kbpsimat_).nnzrow(lrindex);
                for (int p2 = 0; p2 < nnzrow2; p2++)
                {
                    const double alpha
                        = kbpsi1
                          * (*kbpsi2->kbpsimat_).getRowEntry(lrindex, p2);
                    /* set hnlij */
                    if (fabs(alpha) > tolKBpsi)
                    {
                        const int st2
                            = (*kbpsi2->kbpsimat_).getColumnIndex(lrindex, p2);
                        mat.insertMatrixElement(st1, st2, alpha, ADD, true);
                    }
                }
            }
        }
    }
}

SquareSubMatrix<double> KBPsiMatrixSparse::computeHvnlMatrix(
    const Ions& ions) const
{
    return computeHvnlMatrix(this, ions);
}

void KBPsiMatrixSparse::computeHvnlMatrix(
    const Ions& ions, VariableSizeMatrix<sparserow>& mat) const
{
    computeHvnlMatrix(this, ions, mat);
}

void KBPsiMatrixSparse::computeHvnlMatrix(
    const Ions& ions, ProjectedMatricesInterface* proj_matrices) const
{
    computeHvnlMatrix(this, ions, proj_matrices);
}

#ifdef MGMOL_USE_SCALAPACK
template <>
void KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions,
    dist_matrix::DistMatrix<DISTMATDTYPE>& hij) const
{
    SquareSubMatrix<double> submat(computeHvnlMatrix(kbpsi2, ions));

    SquareSubMatrix2DistMatrix* ss2dm = SquareSubMatrix2DistMatrix::instance();
    ss2dm->accumulate(submat, hij, 0.);
}
#endif

template <>
void KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions,
    ReplicatedMatrix& hij) const
{
    SquareSubMatrix<double> submat(computeHvnlMatrix(kbpsi2, ions));

    hij.init(submat.data(), submat.ld());

    hij.consolidate();
}

// build <P|phi> elements, one atom at a time
SquareSubMatrix<double> KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions) const
{
    computeHvnlMatrix_tm_.start();

    std::vector<int> indexes;
    kbpsimat_->getAllColumnIndexes(indexes);

    SquareSubMatrix<double> Aij(indexes);

    // Loop over ions centered on current PE only
    // (distribution of work AND Hvnlij contributions)
    for (const auto& ion : ions.local_ions())
    {
        computeHvnlElementsIon((KBPsiMatrixSparse*)kbpsi2, *ion, Aij);
    }

    computeHvnlMatrix_tm_.stop();

    return Aij;
}

// build <P|phi> elements, one atom at a time
void KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions,
    VariableSizeMatrix<sparserow>& mat) const
{
    computeHvnlMatrix_tm_.start();

    // Loop over ions centered on current PE only
    // (distribution of work AND Hvnlij contributions)
    for (const auto& ion : ions.local_ions())
    {
        computeHvnlElementsIon((KBPsiMatrixSparse*)kbpsi2, *ion, mat);
    }

    computeHvnlMatrix_tm_.stop();
}

void KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions,
    ProjectedMatricesInterface* proj_matrices) const
{
    SquareSubMatrix<double> hnlij(computeHvnlMatrix(kbpsi2, ions));
    proj_matrices->setLocalMatrixElementsHnl(hnlij);
}

template <class T>
void KBPsiMatrixSparse::computeAll(const Ions& ions, T& orbitals)
{
    assert(count_proj_subdomain_ == ions.countProjectorsSubdomain());

    if (getIterativeIndex() >= 0)
        if (orbitals.getIterativeIndex() == getIterativeIndex())
        {
            MGmol_MPI& mmpi = *(MGmol_MPI::instance());
            mmpi.barrier();
#ifdef PRINT_OPERATIONS
            if (onpe0)
                (*MPIdata::sout) << "KBPsi coeff. up to date, "
                                    "KBPsiMatrixSparse::computeAll() skipped"
                                 << endl;
#endif
            return;
        }

#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "KBPsiMatrixSparse::computeAll()" << endl;
#endif
    reset();

    const int iinit = 0;
    const int iend  = orbitals.chromatic_number();
    const int bsize = 32;

    setIterativeIndex(orbitals.getIterativeIndex());

    for (int color = iinit; color < iend; color += bsize)
    {
        int ncolors = std::min(bsize, iend - color);
        computeKBpsi(ions, orbitals, color, ncolors, 0);
        if (lapop_) computeKBpsi(ions, orbitals, color, ncolors, 1);
    }

    // consolidate values with overlap of local projectors on other subdomains
    globalSumKBpsi();

    scaleWithKBcoeff(ions);
}

void KBPsiMatrixSparse::printTimers(std::ostream& os)
{
    computeHvnlMatrix_tm_.print(os);
    global_sum_tm_.print(os);
    computeLocalElement_tm_.print(os);
    compute_kbpsi_tm_.print(os);
    setup_tm_.print(os);
    trace_tm_.print(os);
}

double KBPsiMatrixSparse::getEvnl(
    const Ions& ions, ProjectedMatricesSparse* proj_matrices)
{
    double trace = 0.0;
    // loop over all the ions
    // parallelization over ions by including only those centered in subdomain
    for (const auto& ion : ions.local_ions())
    {
        std::vector<int> gids;
        ion->getGidsNLprojs(gids);

        const short nprojs = (short)gids.size();
        for (short i = 0; i < nprojs; i++)
        {
            const int gid = gids[i];
            trace += getTraceDM(gid, proj_matrices->getDM());
        }
    }

    /* gather trace result */
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    double evnl = 0.0;
    mmpi.allreduce(&trace, &evnl, 1, MPI_SUM);

    return evnl * Ry2Ha;
}

template <class MatrixType>
double KBPsiMatrixSparse::getEvnl(
    const Ions& ions, ProjectedMatrices<MatrixType>* proj_matrices)
{
    SquareLocalMatrices<double, MemorySpace::Host> dm(
        proj_matrices->getReplicatedDM());
    double* replicated_dm = dm.getRawPtr();

    double trace = 0.0;
    // loop over all the ions
    // parallelization over ions by including only those centered in subdomain
    for (const auto& ion : ions.local_ions())
    {
        std::vector<int> gids;
        ion->getGidsNLprojs(gids);

        const short nprojs = (short)gids.size();
        for (short i = 0; i < nprojs; i++)
        {
            const int gid = gids[i];
            trace += getTraceDM(gid, replicated_dm, proj_matrices->dim());
        }
    }

    /* gather trace result */
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    double evnl = 0.0;
    mmpi.allreduce(&trace, &evnl, 1, MPI_SUM);

    return evnl * Ry2Ha;
}

double KBPsiMatrixSparse::getTraceDM(
    const int gid, const double* const mat_X, const int numst) const
{
    trace_tm_.start();

    double trace = 0.;

    int* rindex = (int*)(*kbpsimat_).getTableValue(gid);
    if (rindex == nullptr)
    {
        trace_tm_.stop();
        return trace;
    }

    const int lrindex = *rindex;
    const int nnzrow1 = kbpsimat_->nnzrow(lrindex);
    for (int p1 = 0; p1 < nnzrow1; p1++)
    {
        const int st1            = kbpsimat_->getColumnIndex(lrindex, p1);
        const double t1          = (*kbpsimat_).getRowEntry(lrindex, p1);
        const double* const pmat = &mat_X[st1 * numst];

        for (int p2 = 0; p2 < nnzrow1; p2++)
        {
            const int st2 = kbpsimat_->getColumnIndex(lrindex, p2);

            trace += t1 * (*kbpsimat_).getRowEntry(lrindex, p2) * pmat[st2];
        }
    }

    trace_tm_.stop();

    return trace;
}

double KBPsiMatrixSparse::getTraceDM(
    const int gid, const DensityMatrixSparse& dm) const
{
    trace_tm_.start();

    double trace = 0.;

    int* rindex = (int*)(*kbpsimat_).getTableValue(gid);
    if (rindex == nullptr)
    {
        trace_tm_.stop();
        return trace;
    }

    const int lrindex = *rindex;
    const int nnzrow1 = (*kbpsimat_).nnzrow(lrindex);

    std::vector<int> cols;
    cols.reserve(nnzrow1);
    kbpsimat_->getColumnIndexes(lrindex, cols);
    assert(static_cast<int>(cols.size()) == nnzrow1);

    // get values in row of kbpsimat_
    std::vector<double> vval;
    vval.reserve(nnzrow1);
    kbpsimat_->getRowEntries(lrindex, vval);
    assert(vval.size() == cols.size());

    for (int p1 = 0; p1 < nnzrow1; p1++)
    {
        const int st1   = cols[p1];
        const double t1 = vval[p1];

        // get values in row of DM for specific columns
        std::vector<double> row_values;
        dm.getEntries(st1, cols, row_values);
        assert(row_values.size() == cols.size());

        for (int p2 = 0; p2 < nnzrow1; p2++)
        {
            trace += t1 * vval[p2] * row_values[p2];
        }
    }

    trace_tm_.stop();

    return trace;
}

template void KBPsiMatrixSparse::computeKBpsi(const Ions& ions,
    LocGridOrbitals<ORBDTYPE>& orbitals, const int first_color,
    const int nb_colors, const bool flag);
template void KBPsiMatrixSparse::computeAll(
    const Ions&, LocGridOrbitals<ORBDTYPE>&);

template void KBPsiMatrixSparse::computeKBpsi(const Ions& ions,
    ExtendedGridOrbitals<ORBDTYPE>& orbitals, const int first_color,
    const int nb_colors, const bool flag);
template void KBPsiMatrixSparse::computeAll(
    const Ions&, ExtendedGridOrbitals<ORBDTYPE>&);

#ifdef MGMOL_USE_SCALAPACK
template double KBPsiMatrixSparse::getEvnl(const Ions& ions,
    ProjectedMatrices<dist_matrix::DistMatrix<double>>* proj_matrices);
#endif
template double KBPsiMatrixSparse::getEvnl(
    const Ions& ions, ProjectedMatrices<ReplicatedMatrix>* proj_matrices);
