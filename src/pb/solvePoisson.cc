// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <iostream>

#include <mpi.h>

#include "GridFunc.h"
#include "Laph4M.h"
#include "PBh4M.h"
#include "PEenv.h"
#include "SolverLap.h"
#include "SolverPB.h"

// parameters for continuum solvation model
double e0    = 78.36;
double rho0  = 0.0004;
double drho0 = 1.3;

#ifdef ADD_ // LINUX, SGI, SUN
#define solvepoisson solvepoisson_
#define mymain mymain_
#endif

// Function solvepoisson
//
// solves Poisson problem: -nabla (epsilon nabla phi)=4*PI*rhs
//
// Input:
//   int pbtype       problem type:
//                      0 for Poisson
//                      1 for Poisson with dielectric
//   int nx, ny, nz:  global grid size assuming storage iz+nz*iy+nz*ny*ix
//   double lx,ly,lz: computation cell dimensions
//   int bcx,bcy,bcz: boundary conditions:
//                      0 for 0 Dirichlet
//                      1 for periodic
//   double* rhs:     right hand side of Poisson equation
//   double* rhoe:    electronic density used to build dielectric epsilon
//                    (not referred if pbtype==0)
//   int niter:       number of V-cycle
//   double* phi:     initial guess of solution of Poisson problem
// Output:
//   double* phi:     solution of Poisson problem
//   double* phi_eps: additional KS potential due to dielectric
extern "C"
{
    void solvepoisson(const int* pbtype, const int* const nx,
        const int* const ny, const int* const nz, const double* const lx,
        const double* const ly, const double* const lz, const int* const bcx,
        const int* const bcy, const int* const bcz, const double* rhs,
        const double* rhoe, const int* const niter, double* phi,
        double* phi_eps)
    {
        if (*pbtype != 0 && *pbtype != 1)
        {
            cout << "solvePoisson: pbtype should be 0 or 1" << endl;
            return;
        }
        if (*bcx != 0 && *bcx != 1)
        {
            cout << "solvePoisson: bcx should be 0 or 1" << endl;
            return;
        }
        if (*bcy != 0 && *bcy != 1)
        {
            cout << "solvePoisson: bcy should be 0 or 1" << endl;
            return;
        }
        if (*bcz != 0 && *bcz != 1)
        {
            cout << "solvePoisson: bcz should be 0 or 1" << endl;
            return;
        }
        if (*pbtype == 1 && *bcx == 1 && *bcy == 1 && *bcz == 1)
        {
            cout << "solvePoisson: periodic BC not implemented for dielectric "
                    "model"
                 << endl;
            return;
        }

        // leading dimension
        const int ldz = *nz;

        GridFunc<GFDTYPE>*
            gf_veps; // Additional potential in KS for solvation model
        GridFunc<GFDTYPE>* gf_rhoe;

        // Parallel environment
        PEenv myPEenv(*nx, *ny, *nz);
        // cout<<"My color is "<< myPEenv.color()<<endl;

        Grid* myGrid;
        SolverPB<PBh4M>* solverPB;
        SolverLap<Laph4M>* solver;
        if (myPEenv.color() == 0)
        {
            if (*pbtype)
            {
                myGrid = new Grid(*lx, *ly, *lz, *nx, *ny, *nz, myPEenv, 2);
                PBh4M oper(*myGrid, e0, rho0, drho0);
                solverPB = new SolverPB<PBh4M>(oper, *bcx, *bcy, *bcz);
                gf_veps  = new GridFunc<GFDTYPE>(*myGrid, 1);
                gf_rhoe
                    = new GridFunc<GFDTYPE>(rhoe, *myGrid, 1, 1, 1, 'g', ldz);
            }
            else
            {
                myGrid = new Grid(*lx, *ly, *lz, *nx, *ny, *nz, myPEenv, 1);
                Laph4M oper(*myGrid);
                solver = new SolverLap<Laph4M>(oper, *bcx, *bcy, *bcz);
            }

            GridFunc<GFDTYPE> gf_phi(
                phi, *myGrid, *bcx, *bcy, *bcz, 'g', ldz); // Coulomb potential
            GridFunc<GFDTYPE> gf_rhs(rhs, *myGrid, 1, 1, 1, 'g', ldz); // r.h.s

            double nv = norm(gf_phi);
            if (myPEenv.onpe0()) cout << " Norm of gf_phi=" << nv << endl;
            nv = norm(gf_rhs);
            if (myPEenv.onpe0()) cout << " Norm of gf_rhs=" << nv << endl;

            // cout<<"Call solver on PE "<<myPEenv.mytask()<<endl;
            double res;
            if (*pbtype)
            {
                if (myPEenv.onpe0()) cout << "PB Solver" << endl;
                (*gf_veps) = 0.;
                res        = solverPB->solve(
                    gf_phi, gf_rhs, *gf_rhoe, *gf_veps, *niter);
                gf_veps->init_vect(phi_eps, 'g', ldz);
            }
            else
            {
                if (myPEenv.onpe0()) cout << "Laplacian Solver" << endl;
                res = solver->solve(gf_phi, gf_rhs, *niter);
            }
            gf_phi.init_vect(phi, 'g', ldz);
            if (myPEenv.onpe0())
                cout << "Residual after solve = " << res << endl;

            // delete objects
            if (*pbtype)
            {
                delete gf_veps;
                delete gf_rhoe;
                delete solverPB;
            }
            else
            {
                delete solver;
            }
            delete myGrid;
        }
        int npes;
        MPI_Comm_size(MPI_COMM_WORLD, &npes);
        // Send data to PE outside of solver partition
        if (npes > myPEenv.n_mpi_tasks())
        {
            int nn         = npes - myPEenv.n_mpi_tasks();
            const int size = (*nx) * (*ny) * ldz;
            if (myPEenv.mytask() < nn && myPEenv.color() == 0)
            {
                // cout<<"Send data from "<<myPEenv.mytask()
                //    <<" to "<<npes-1-myPEenv.mytask()<<endl;
                MPI_Send(phi, size, MPI_DOUBLE, npes - 1 - myPEenv.mytask(), 0,
                    MPI_COMM_WORLD);
                MPI_Send(phi_eps, size, MPI_DOUBLE, npes - 1 - myPEenv.mytask(),
                    1, MPI_COMM_WORLD);
            }
            else if (myPEenv.color() != 0)
            {
                MPI_Status status;
                // cout<<"Receive data from "<<npes-1-myPEenv.mytask()
                //    <<" to "<<myPEenv.mytask()<<endl;
                MPI_Recv(phi, size, MPI_DOUBLE, myPEenv.mytask(), 0,
                    MPI_COMM_WORLD, &status);
                MPI_Recv(phi_eps, size, MPI_DOUBLE, myPEenv.mytask(), 1,
                    MPI_COMM_WORLD, &status);
            }
        }
    }
}

// Need a small wrapper for the fortran main program, since fortran code should
// not be used as the main program
extern "C" int mymain();

int main(int argc, char** argv)
{
    mymain();

    return 0;
}
