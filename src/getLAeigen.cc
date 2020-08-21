// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Control.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrixSparse.h"
#include "LocGridOrbitals.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "Potentials.h"
#include "Preconditioning.h"
#include "ProjectedMatrices.h"

#include "MGmol_prototypes.h"
#include <float.h>

#ifdef HAVE_ARPACK
#include "arpack_mangle.h"

typedef int LOGICAL;

extern "C"
{
    void pdsaupd(int*, int*, const char* bmat, const int* const,
        const char* which, const int* const, const double* const, double* resid,
        const int* const, double* V, const int* const, int* iparam, int* ipntr,
        double* workd, double* workl, const int* const, int* info);
    void pdseupd(int*, LOGICAL&, const char* const, LOGICAL*, double*, double*,
        const int* const, const double* const, const char* bmat,
        const int* const, const char* which, const int* const,
        const double* const, double* resid, const int* const, double* V,
        const int* const, int* iparam, int* ipntr, double* workd, double* workl,
        const int* const, int* info);
}

extern Hamiltonian* hamiltonian;

void get_kbpsi(
    KBPsiMatrixSparse& kbpsi, Ions& ions, pb::GridFunc<ORBDTYPE>* phi)
{
    kbpsi.reset();

    kbpsi.computeKBpsi(ions, phi, 0, 0);
    kbpsi.computeKBpsi(ions, phi, 0, 1);

    kbpsi.globalSumKBpsi();

    kbpsi.scaleWithKBcoeff(ions);
}

template <typename T>
void matvec(pb::GridFunc<ORBDTYPE>& gfpsi, double* hpsi,
    KBPsiMatrixSparse& kbpsi, Ions& ions, Preconditioning<T>* precond,
    const double shift)
{
    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    Potentials& pot        = hamiltonian->potential();

    const int numpt = mymesh->numpt();
    double* work    = new double[numpt];

    gfpsi.init_vect(work, 'd');

    // gf_work = -Lap*psi
    pb::GridFunc<ORBDTYPE> gf_work(mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    hamiltonian->lapOper()->apply(gfpsi, gf_work);

    gf_work.init_vect(hpsi, 'd');

    my_daxpy(numpt, shift, work, hpsi);
#if 1
    // gf_work_v = Vtot*psi
    const double* const vtot = pot.vtot();
    pb::GridFunc<ORBDTYPE> gfpot(mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    gfpot.assign(vtot);
    pb::GridFunc<ORBDTYPE> gfvw1(mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    pb::GridFunc<ORBDTYPE> gf_work_v(
        mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    gfvw1.prod(gfpsi, gfpot);
    gfvw1.trade_boundaries();

    hamiltonian->lapOper()->rhs(gfvw1, gf_work_v);
    gf_work_v.init_vect(work, 'd');

    my_daxpy(numpt, 1., work, hpsi);
#endif

#if 1
    // Vnl*psi
    get_kbpsi(kbpsi, ions, &gfpsi);

    vector<int> ptr_func(mymesh->subdivx(), 0);
    get_vnlpsi(ions, ptr_func, kbpsi, work);
    pb::GridFunc<ORBDTYPE> gf_worknl(
        mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);
    gf_worknl.assign(work);

    gf_worknl.trade_boundaries();
    hamiltonian->lapOper()->rhs(gf_worknl, work);

    // Add the contribution of the non-local potential to H phi
    my_daxpy(numpt, 2., work, hpsi);
#endif
    delete[] work;

    // convert to Hartree units
    int ione    = 1;
    double half = 0.5;
    DSCAL(&numpt, &half, hpsi, &ione);

    if (precond != NULL)
    {
        precond->mg(hpsi, hpsi);
    }
}

template <typename T>
double getLAeigen(const double tol, const int maxit, Ions& ions)
{
    if (onpe0)
    {
        (*MPIdata::sout) << "ARPACK to compute largest algebraic eigenvalue"
                         << endl;
        (*MPIdata::sout) << "maxit=" << maxit << endl;
    }

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const int numpt = mymesh->numpt();

    const int nev = 1; // Number of eigenvalues to be computed
    const int ncv = 25; // number of Lanczos vectors. Should be > nev
    char bmat     = 'I';
    char which[3];
    strcpy(which, "LA");
    const int mode = 1; // mode 1: A*x = lambda*x, A symmetric

    int info       = 0;
    int ishift     = 1;
    int iparam[11] = { ishift, 0, maxit, 1, 0, 0, mode, 0, 0, 0, 0 };

    int ipntr[11]; // output

    double* workd    = new double[3 * numpt];
    const int lworkl = ncv * (ncv + 8);
    double* workl    = new double[lworkl];
    double* v        = new double[numpt * ncv];
    double* resid    = new double[numpt];

    int ldv = numpt;

    Preconditioning<T>* precond
        = new Preconditioning<T>(ct.lap_type, ct.mg_levels_, mygrid, ct.bc);
    vector<vector<int>> ptr_func;
    vector<GridMask*> st2mask;
    extern LocGridOrbitals* current_orbitals;
    precond->setup(st2mask, ptr_func);
    precond->setGamma(current_orbitals->get_gamma());

    KBPsiMatrixSparse kbpsi(hamiltonian->lapOper());
    kbpsi.allocate(ions, nev);

    pb::GridFunc<ORBDTYPE> gfpsi(mygrid, ct.bcWF[0], ct.bcWF[1], ct.bcWF[2]);

    // main loop
    int ido         = 0;
    int it          = 1;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.comm();
    do
    {

        // if ( onpe0 )(*MPIdata::sout)<<"Iteration "<<it<<endl;
        pdsaupd(&comm, &ido, &bmat, &numpt, &which[0], &nev, &tol, resid, &ncv,
            v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

        assert(ido != 2);
        assert(ido != 3);
#if 0
        (*MPIdata::sout)<<"ido="<<ido<<endl;
        
        if ( onpe0 ){
            (*MPIdata::sout)<<"# converged values: "<<iparam[4]<<endl;
            double* ritz_val=workl+ipntr[5]-1;
            for(int i=0;i<ncv;i++){
                (*MPIdata::sout)<<"Ritz value   ["<<i<<"]="<<ritz_val[i]<<endl;
            }
        }
#endif
        if (ido == -1 || ido == 1)
        {
            gfpsi.assign(&workd[ipntr[0] - 1]);

            matvec(gfpsi, &workd[ipntr[1] - 1], kbpsi, ions, precond, 0.);
        }
        it++;
    } while (ido != 99);

    double val = DBL_MAX;
    if (info < 0)
    {
        (*MPIdata::sout) << "info=" << info << endl;
    }
    else
    {
        LOGICAL rvec    = 0;
        char howmny     = 'S';
        LOGICAL* select = new LOGICAL[ncv];
        for (int i = 0; i < ncv; i++)
            select[i] = 1;
        double* d = new double[2 * ncv];
        double* z = new double[1];
        int ldz   = 1;
        double sigma;
        pdseupd(&comm, rvec, &howmny, select, d, z, &ldz, &sigma, &bmat, &numpt,
            &which[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd,
            workl, &lworkl, &info);
        if (onpe0)
        {
            for (int i = 0; i < nev; i++)
            {
                (*MPIdata::sout)
                    << "Computed Ritz value   [" << i << "]=" << d[i] << endl;
            }
            (*MPIdata::sout)
                << "The number of Ritz values requested is " << nev << endl;
            (*MPIdata::sout)
                << "The number of Arnoldi vectors generated is " << ncv << endl;
            (*MPIdata::sout)
                << "What portion of the spectrum: " << which << endl;
            (*MPIdata::sout) << "# iterations = " << iparam[2] << endl;
            (*MPIdata::sout) << "The number of converged Ritz values is "
                             << iparam[4] << endl;
            (*MPIdata::sout) << "The number of OP*x is " << iparam[8] << endl;
            (*MPIdata::sout)
                << "The convergence criterion is " << scientific << tol << endl;
            (*MPIdata::sout) << "info=" << info << endl;
        }

        val = d[0];

        delete[] d;
        delete[] z;
        delete[] select;
    }

    if (info == 1 && onpe0)
        (*MPIdata::sout) << "Maximum number of iterations reached." << endl;

    delete[] v;
    delete[] workl;
    delete[] workd;
    delete[] resid;
    delete precond;

    return val;
}

#endif
