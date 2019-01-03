// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "KBPsiMatrixSparse.h"
#include "Control.h"
#include "Ions.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol_MPI.h"
#include "Mesh.h"
#include "ProjectedMatrices.h"
#include "ProjectedMatricesSparse.h"

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
    kbpsimat_    = 0;
    kbBpsimat_   = 0;
    distributor_ = 0;

    isDataSetup_ = false;
};

KBPsiMatrixSparse::~KBPsiMatrixSparse() { clearData(); }

void KBPsiMatrixSparse::clearData()
{
    if (isDataSetup_)
    {
        delete kbpsimat_;
        kbpsimat_ = 0;
        if (lapop_)
        {
            delete kbBpsimat_;
        }
        kbBpsimat_ = 0;
        delete distributor_;
        distributor_ = 0;

        isDataSetup_ = false;
    }
}

template <class T>
void KBPsiMatrixSparse::setup(const Ions& ions, const T& orbitals)
{
    setup_tm_.start();

    setOutdated();

    numst_ = orbitals.numst();

    // clear old data
    clearData();

    // number of projectors overlapping with subdomain
    count_proj_subdomain_ = ions.countProjectorsSubdomain();
    /* setup matrices */
    assert(kbpsimat_ == 0);
    assert(kbBpsimat_ == 0);
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
    assert(distributor_ == 0);
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
void KBPsiMatrixSparse::computeKBpsi(Ions& ions, T& orbitals,
    const int first_color, const int nb_colors, const bool flag)
{
    assert(first_color >= 0);
    assert(first_color < 100000);
    assert(nb_colors > 0);
    assert(nb_colors < 100000);
    assert(orbitals.getIterativeIndex() > 0);

    compute_kbpsi_tm_.start();

    ORBDTYPE* ppsi   = 0;
    const int ldsize = orbitals.getLda();
    assert(ldsize > 0);
    assert(ldsize < 1e8);

    if (flag)
    {
        assert(lapop_ != NULL);
        ppsi = new ORBDTYPE[nb_colors * ldsize];
        orbitals.setDataWithGhosts();
        orbitals.trade_boundaries();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int color = 0; color < nb_colors; color++)
        {
            lapop_->rhs(orbitals.getFuncWithGhosts(first_color + color),
                ppsi + color * ldsize);
        }
    }
    else
    {
        ppsi = orbitals.getPsi(first_color);
    }

    // Loop over functions, subdomains and ions
    const vector<Ion*>::const_iterator iend = ions.overlappingNL_ions().end();
    Mesh* mymesh                            = Mesh::instance();
    const int subdivx                       = mymesh->subdivx();
    for (int iloc = 0; iloc < subdivx; iloc++)
    {
// Threading here leads to results slightly dependent on number of threads (jlf,
// 07/15/2016)
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int color = 0; color < nb_colors; color++)
        {
            const int gid = orbitals.getGlobalIndex(iloc, first_color + color);

            // Loop over the ions
            if (gid != -1)
            {
                vector<Ion*>::const_iterator ion
                    = ions.overlappingNL_ions().begin();
                while (ion != iend)
                {
                    computeLocalElement(
                        **ion, gid, iloc, ppsi + color * ldsize, flag);
                    ion++;
                }
            }
        }
    }

    if (flag) delete[] ppsi;

    compute_kbpsi_tm_.stop();
}

void KBPsiMatrixSparse::computeKBpsi(
    Ions& ions, pb::GridFunc<ORBDTYPE>* phi, const int istate, const bool flag)
{
    assert(lapop_ != NULL);
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
    const vector<Ion*>::const_iterator iend = ions.overlappingNL_ions().end();
    const int subdivx                       = mymesh->subdivx();
    for (int iloc = 0; iloc < subdivx; iloc++)
    {
        // Loop over the ions
        vector<Ion*>::const_iterator ion = ions.overlappingNL_ions().begin();
        while (ion != iend)
        {
            computeLocalElement(**ion, istate, iloc, ppsi, flag);
            ion++;
        }
    }

    delete[] ppsi;

    compute_kbpsi_tm_.stop();
}

void KBPsiMatrixSparse::scaleWithKBcoeff(const Ions& ions)
{
    //    int ione=1;
    vector<Ion*>::const_iterator ion  = ions.overlappingNL_ions().begin();
    vector<Ion*>::const_iterator iend = ions.overlappingNL_ions().end();
    while (ion != iend)
    {
        vector<int> gids;
        (*ion)->getGidsNLprojs(gids);
        vector<double> kbcoeffs;
        (*ion)->getKBcoeffs(kbcoeffs);

        const short nprojs = (short)gids.size();
        for (short i = 0; i < nprojs; i++)
        {
            const int gid = gids[i];
            double coeff  = kbcoeffs[i];

            // loop over states to multiply kbpsi_[st][gid] and kbBpsi_[st][gid]
            // by coeff
            (*kbpsimat_).scaleRow(gid, coeff);
            if (lapop_) (*kbBpsimat_).scaleRow(gid, coeff);
        }

        ion++;
    }
}

// computeHvnlMatrix:
//    Get the elements of the Hamiltonian matrix due to the non-local
//    potential, and add them into Aij.
// Note: neglecting the small matrix elements reduces the size of hnlij and thus
//       reduces the size of communications later on.
void KBPsiMatrixSparse::computeHvnlMatrix(const KBPsiMatrixSparse* const kbpsi2,
    const Ion& ion, dist_matrix::SparseDistMatrix<DISTMATDTYPE>& hnlij) const
{
    assert(ion.here());

    vector<int> gids;
    ion.getGidsNLprojs(gids);
    vector<short> kbsigns;
    ion.getKBsigns(kbsigns);
    const short nprojs = (short)gids.size();
    for (short i = 0; i < nprojs; i++)
    {
        const int gid      = gids[i];
        const double coeff = (double)kbsigns[i];

        // double loop over states to fill hnlij[st1][st2] (in general not
        // symmetric... )
        int* rindex = (int*)(*kbBpsimat_).getTableValue(gid);
        if (rindex == NULL) continue;
        const int lrindex = *rindex;
        const int nnzrow1 = kbBpsimat_->nnzrow(lrindex);
        for (int p1 = 0; p1 < nnzrow1; p1++)
        {
            const int st1 = kbBpsimat_->getColumnIndex(lrindex, p1);
            const double kbpsi1
                = coeff * kbBpsimat_->getRowEntry(lrindex, p1);
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
                        hnlij.push_back(st1, st2, alpha);
                    }
                }
            }
        }
    }
}

void KBPsiMatrixSparse::computeHvnlMatrix(const KBPsiMatrixSparse* const kbpsi2,
    const Ion& ion, VariableSizeMatrix<sparserow>& mat) const
{
    assert(ion.here());

    vector<int> gids;
    ion.getGidsNLprojs(gids);

    vector<short> coeffs;
    ion.getKBsigns(coeffs);

    const short nprojs = (short)gids.size();
    for (short i = 0; i < nprojs; i++)
    {
        const int gid      = gids[i];
        const double coeff = (double)coeffs[i];

        // double loop over states to fill hnlij[st1][st2] (in general not
        // symmetric... )
        int* rindex = (int*)(*kbBpsimat_).getTableValue(gid);
        if (rindex == NULL) continue;
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

void KBPsiMatrixSparse::computeHvnlMatrix(const KBPsiMatrixSparse* const kbpsi2,
    const Ion& ion, ProjectedMatricesInterface* proj_matrices) const
{
    assert(ion.here());
    assert(proj_matrices != 0);

    vector<int> gids;
    ion.getGidsNLprojs(gids);

    vector<short> coeffs;
    ion.getKBsigns(coeffs);

    const short nprojs = (short)gids.size();
    for (short i = 0; i < nprojs; i++)
    {
        const int gid      = gids[i];
        const double coeff = (double)coeffs[i];

        // double loop over states to fill hnlij[st1][st2] (in general not
        // symmetric... )
        int* rindex = (int*)(*kbBpsimat_).getTableValue(gid);
        if (rindex == NULL) continue;
        const int lrindex = *rindex;
        const int nnzrow1 = (*kbBpsimat_).nnzrow(lrindex);
        for (int p1 = 0; p1 < nnzrow1; p1++)
        {
            const int st1 = (*kbBpsimat_).getColumnIndex(lrindex, p1);
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
                        proj_matrices->addMatrixElementSparseH(st1, st2, alpha);
                    }
                }
            }
        }
    }
}

void KBPsiMatrixSparse::computeHvnlMatrix(
    const Ions& ions, dist_matrix::SparseDistMatrix<DISTMATDTYPE>& Aij) const
{
    computeHvnlMatrix(this, ions, Aij);
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

void KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions,
    dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE>& Aij) const
{
    computeHvnlMatrix(kbpsi2, ions, Aij.sparse());
}

// build <P|phi> elements, one atom at a time
void KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions,
    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& Aij) const
{
    computeHvnlMatrix_tm_.start();

    // Loop over ions centered on current PE only
    // (distribution of work AND Hvnlij contributions)
    vector<Ion*>::const_iterator ion = ions.local_ions().begin();
    while (ion != ions.local_ions().end())
    {
        computeHvnlMatrix((KBPsiMatrixSparse*)kbpsi2, **ion, Aij);
        ion++;
    }

    computeHvnlMatrix_tm_.stop();
}

// build <P|phi> elements, one atom at a time
void KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions,
    VariableSizeMatrix<sparserow>& mat) const
{
    computeHvnlMatrix_tm_.start();

    // Loop over ions centered on current PE only
    // (distribution of work AND Hvnlij contributions)
    vector<Ion*>::const_iterator ion = ions.local_ions().begin();
    while (ion != ions.local_ions().end())
    {
        computeHvnlMatrix((KBPsiMatrixSparse*)kbpsi2, **ion, mat);
        ion++;
    }

    computeHvnlMatrix_tm_.stop();
}

void KBPsiMatrixSparse::computeHvnlMatrix(
    const KBPsiMatrixInterface* const kbpsi2, const Ions& ions,
    ProjectedMatricesInterface* proj_matrices) const
{
    computeHvnlMatrix_tm_.start();

    // Loop over ions centered on current PE only
    // (distribution of work AND Hvnlij contributions)
    vector<Ion*>::const_iterator ion = ions.local_ions().begin();
    while (ion != ions.local_ions().end())
    {
        computeHvnlMatrix((KBPsiMatrixSparse*)kbpsi2, **ion, proj_matrices);
        ion++;
    }

    computeHvnlMatrix_tm_.stop();
}

// build elements of matrix <phi_i|Vnl|phi_j> (assumed to be symmetric)
// assemble resulting matrix in variable sparse matrix format
void KBPsiMatrixSparse::getPsiKBPsiSym(
    const Ion& ion, VariableSizeMatrix<sparserow>& sm)
{
    vector<int> gids;
    ion.getGidsNLprojs(gids);
    vector<short> kbsigns;
    ion.getKBsigns(kbsigns);

    const short nprojs = (short)gids.size();
    for (short i = 0; i < nprojs; i++)
    {
        const int gid      = gids[i];
        const double coeff = (double)kbsigns[i];
        int* rindex        = (int*)(kbpsimat_->getTableValue(gid));
        if (rindex == NULL) continue;
        const int lrindex = *rindex;
        const int nnzrow1 = kbpsimat_->nnzrow(lrindex);
        for (int p1 = 0; p1 < nnzrow1; p1++)
        {
            double kbpsielement1 = kbpsimat_->getRowEntry(lrindex, p1);
            if (fabs(kbpsielement1) <= tolKBpsi)
                continue;
            const int st1 = kbpsimat_->getColumnIndex(lrindex, p1);
            for (int p2 = 0; p2 < nnzrow1; p2++)
            {
                double kbpsielement2 = kbpsimat_->getRowEntry(lrindex, p2);
                if (fabs(kbpsielement2) <= tolKBpsi)
                    continue;
                const double alpha = coeff * kbpsielement1 * kbpsielement2;
                /* set hnlij */
                if (fabs(alpha) > tolKBpsi)
                {
                    const int st2 = kbpsimat_->getColumnIndex(lrindex, p2);
                    sm.insertMatrixElement(st1, st2, alpha, ADD, true);
                }
            }
        }
    }
}

void KBPsiMatrixSparse::getPsiKBPsiSym(
    const Ions& ions, VariableSizeMatrix<sparserow>& sm)
{
    // loop over all the ions
    // parallelization over ions by including only those centered in subdomain
    for (auto& ion : ions.local_ions())
    {
        getPsiKBPsiSym(*ion, sm);
    }
}

template <class T>
void KBPsiMatrixSparse::computeAll(Ions& ions, T& orbitals)
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
        int ncolors = min(bsize, iend - color);
        computeKBpsi(ions, orbitals, color, ncolors, 0);
        if (lapop_) computeKBpsi(ions, orbitals, color, ncolors, 1);
    }

    // consolidate values with overlap of local projectors on other subdomains
    globalSumKBpsi();

    scaleWithKBcoeff(ions);
}

void KBPsiMatrixSparse::printTimers(ostream& os)
{
    computeHvnlMatrix_tm_.print(os);
    global_sum_tm_.print(os);
    computeLocalElement_tm_.print(os);
    compute_kbpsi_tm_.print(os);
    setup_tm_.print(os);
    trace_tm_.print(os);
}

template <class T>
double KBPsiMatrixSparse::getEvnl(const Ions& ions, T& orbitals,
    ProjectedMatricesSparse* proj_matrices)
{
    const int numst = orbitals.numst();
    if (numst == 0) return 0.;

    double trace = 0.0;
    // loop over all the ions
    // parallelization over ions by including only those centered in subdomain
    for (auto& ion : ions.local_ions())
    {
        vector<int> gids;
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
    MPI_Comm comm   = mmpi.commSpin();

    double evnl = 0.0;
    MPI_Allreduce(&trace, &evnl, 1, MPI_DOUBLE, MPI_SUM, comm);

    return evnl * Ry2Ha;
}

template <class T>
double KBPsiMatrixSparse::getEvnl(const Ions& ions, T& orbitals,
    ProjectedMatrices* proj_matrices)
{
    const int numst = orbitals.numst();
    if (numst == 0) return 0.;

    DISTMATDTYPE* replicated_dm;
    proj_matrices->getReplicatedDM(replicated_dm);

    double trace = 0.0;
    // loop over all the ions
    // parallelization over ions by including only those centered in subdomain
    for (auto& ion : ions.local_ions())
    {
        vector<int> gids;
        ion->getGidsNLprojs(gids);

        const short nprojs = (short)gids.size();
        for (short i = 0; i < nprojs; i++)
        {
            const int gid = gids[i];
            trace += getTraceDM(gid, replicated_dm, numst);
        }
    }

    /* gather trace result */
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSpin();

    double evnl = 0.0;
    MPI_Allreduce(&trace, &evnl, 1, MPI_DOUBLE, MPI_SUM, comm);

    return evnl * Ry2Ha;
}

double KBPsiMatrixSparse::getTraceDM(
    const int gid, const DISTMATDTYPE* const mat_X, const int numst) const
{
    double trace = 0.;

    int* rindex = (int*)(*kbpsimat_).getTableValue(gid);
    if (rindex == NULL) return trace;

    const int lrindex = *rindex;
    const int nnzrow1 = kbpsimat_->nnzrow(lrindex);
    for (int p1 = 0; p1 < nnzrow1; p1++)
    {
        const int st1 = kbpsimat_->getColumnIndex(lrindex, p1);
        const double t1
            = (*kbpsimat_).getRowEntry(lrindex, p1);
        const DISTMATDTYPE* const pmat = &mat_X[st1 * numst];

        for (int p2 = 0; p2 < nnzrow1; p2++)
        {
            const int st2 = kbpsimat_->getColumnIndex(lrindex, p2);

            trace += t1 * (*kbpsimat_).getRowEntry(lrindex, p2)
                     * pmat[st2];
        }
    }

    return trace;
}

double KBPsiMatrixSparse::getTraceDM(
    const int gid, const DensityMatrixSparse& dm) const
{
    trace_tm_.start();

    double trace = 0.;

    int* rindex = (int*)(*kbpsimat_).getTableValue(gid);
    if (rindex == NULL)
    {
        trace_tm_.stop();
        return trace;
    }

    const int lrindex = *rindex;
    const int nnzrow1 = (*kbpsimat_).nnzrow(lrindex);

    vector<int> cols;
    cols.reserve(nnzrow1);
    kbpsimat_->getColumnIndexes(lrindex, cols);
    assert(cols.size() == nnzrow1);

    //get values in row of kbpsimat_
    vector<double> vval;
    vval.reserve(nnzrow1);
    kbpsimat_->getRowEntries(lrindex, vval);
    assert(vval.size() == cols.size());

    for (int p1 = 0; p1 < nnzrow1; p1++)
    {
        const int st1 = cols[p1];
        const double t1 = vval[p1];

        //get values in row of DM for specific columns
        vector<double> row_values;
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

template void KBPsiMatrixSparse::computeKBpsi(Ions& ions,
        LocGridOrbitals& orbitals,
        const int first_color, const int nb_colors, const bool flag);
template double KBPsiMatrixSparse::getEvnl(const Ions& ions,
        LocGridOrbitals& orbitals,
        ProjectedMatricesSparse* proj_matrices);
template double KBPsiMatrixSparse::getEvnl(const Ions& ions,
        LocGridOrbitals& orbitals,
        ProjectedMatrices* proj_matrices);
template void KBPsiMatrixSparse::setup(const Ions&, const LocGridOrbitals&);
template void KBPsiMatrixSparse::computeAll(Ions&, LocGridOrbitals&);

template void KBPsiMatrixSparse::computeKBpsi(Ions& ions,
        ExtendedGridOrbitals& orbitals,
        const int first_color, const int nb_colors, const bool flag);
template double KBPsiMatrixSparse::getEvnl(const Ions& ions,
        ExtendedGridOrbitals& orbitals,
        ProjectedMatricesSparse* proj_matrices);
template double KBPsiMatrixSparse::getEvnl(const Ions& ions,
        ExtendedGridOrbitals& orbitals,
        ProjectedMatrices* proj_matrices);
template void KBPsiMatrixSparse::setup(const Ions&, const ExtendedGridOrbitals&);
template void KBPsiMatrixSparse::computeAll(Ions&, ExtendedGridOrbitals&);

