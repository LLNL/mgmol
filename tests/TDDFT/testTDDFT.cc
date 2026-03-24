// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE


/*********************************************************
 * Evolves dX/dt = -i(HX - XH)
**********************************************************/ 
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "mgmol_run.h"
#include "DensityMatrix.h"
#include "ProjectedMatrices.h"

#include <cassert>
#include <iostream>
#include <time.h>
#include <vector>
#include <iomanip> 
#include <fstream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define LOGGING 1
#define PRINT_INTERVAL 1
/* timestep size (in a.u.) and number of timesteps */
#define	DELTA_T 0.2
#define NUM_TIME_STEPS 20000

int propagate_density_matrix(ReplicatedMatrix& X_real, ReplicatedMatrix& X_imag, ReplicatedMatrix& H, double dt, double tol, int maxits=10);
void compute_commutator(const ReplicatedMatrix& H, const ReplicatedMatrix& X, ReplicatedMatrix& commHX);

int main(int argc, char** argv)
{
    int mpirc = MPI_Init(&argc, &argv);
    if (mpirc != MPI_SUCCESS)
    {
        std::cerr << "MPI Initialization failed!!!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    /*
     * Initialize general things, like magma, openmp, IO, ...
     */
    mgmol_init(comm);

    /*
     * read runtime parameters
     */
    std::string input_filename("");
    std::string lrs_filename;
    std::string constraints_filename("");

    float total_spin = 0.;
    bool with_spin   = false;

    po::variables_map vm;

    // read from PE0 only
    if (MPIdata::onpe0)
    {
        read_config(argc, argv, vm, input_filename, lrs_filename,
            constraints_filename, total_spin, with_spin);
    }

    MGmol_MPI::setup(comm, std::cout, with_spin);
    MGmol_MPI& mmpi      = *(MGmol_MPI::instance());
    MPI_Comm global_comm = mmpi.commGlobal();

    /*
     * Setup control struct with run time parameters
     */
    Control::setup(global_comm, with_spin, total_spin);
    Control& ct = *(Control::instance());

    ct.setOptions(vm);

    int ret = ct.checkOptions();
    if (ret < 0) return ret;

    mmpi.bcastGlobal(input_filename);
    mmpi.bcastGlobal(lrs_filename);

    // Enter main scope
    {
        if (MPIdata::onpe0)
        {
            std::cout << "-------------------------" << std::endl;
            std::cout << "Construct MGmol object..." << std::endl;
            std::cout << "-------------------------" << std::endl;
        }

        MGmolInterface* mgmol;
        if (ct.isLocMode())
            mgmol = new MGmol<LocGridOrbitals<ORBDTYPE>>(global_comm, *MPIdata::sout,
                input_filename, lrs_filename, constraints_filename);
        else
            mgmol = new MGmol<ExtendedGridOrbitals<ORBDTYPE>>(global_comm, *MPIdata::sout,
                input_filename, lrs_filename, constraints_filename);

        if (MPIdata::onpe0)
        {
            std::cout << "-------------------------" << std::endl;
            std::cout << "MGmol setup..." << std::endl;
            std::cout << "-------------------------" << std::endl;
        }
        mgmol->setup();

        if (MPIdata::onpe0)
        {
            std::cout << "-------------------------" << std::endl;
            std::cout << "Setup done..." << std::endl;
            std::cout << "-------------------------" << std::endl;
        }

        // here we just use the atomic positions read in and used
        // to initialize MGmol
        std::vector<double> positions;
        mgmol->getAtomicPositions(positions);
        std::vector<short> anumbers;
        mgmol->getAtomicNumbers(anumbers);
        if (MPIdata::onpe0)
        {
            std::cout << "Positions:" << std::endl;
            std::vector<short>::iterator ita = anumbers.begin();
            for (std::vector<double>::iterator it = positions.begin();
                 it != positions.end(); it += 3)
            {
                std::cout << *ita;
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
                ita++;
            }
        }

        // compute energy and forces using all MPI tasks
        // expect positions to be replicated on all MPI tasks
        std::vector<double> forces;
        double eks
            = mgmol->evaluateEnergyAndForces(positions, anumbers, forces);
        mgmol->dumpRestart();

        // print out results
        if (MPIdata::onpe0)
        {
            std::cout << "Ground state Energy, Eks : " << eks << std::endl;
            std::cout << "Grounde state Forces :" << std::endl;
            for (std::vector<double>::iterator it = forces.begin();
                 it != forces.end(); it += 3)
            {
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
            }
        }
 
        // Initialize TDDFT:
        // Compute energy and forces again using wavefunctions
        // from previous call
        Mesh* mymesh           = Mesh::instance();
        const pb::Grid& mygrid = mymesh->grid();

        /* cast directly to ProjectedMatrices class to use derived class functionality */
        std::shared_ptr<ProjectedMatrices<ReplicatedMatrix>> dprojmatrices = 
            std::dynamic_pointer_cast<ProjectedMatrices<ReplicatedMatrix>>(mgmol->getProjectedMatrices());            

        ExtendedGridOrbitals<ORBDTYPE> orbitals("new_orbitals", mygrid,
            mymesh->subdivx(), ct.numst, ct.bcWF, dprojmatrices.get(), nullptr,
            nullptr, nullptr, nullptr);

        const pb::PEenv& myPEenv = mymesh->peenv();
        HDFrestart h5file("WF", myPEenv, ct.out_restart_file_type);
        orbitals.read_hdf5(h5file);

        /* TDDFT propagation */
        /* propagation params */
        double t = 0.0;
        double dt = DELTA_T;
        int iters = 0, avg_iters = 0.;

        ReplicatedMatrix X_real("X_real", ct.numst);
        X_real = dprojmatrices->dm();
        ReplicatedMatrix X_imag("X_imag", ct.numst);

        /* Hamiltonian object */
        ReplicatedMatrix H_tot("H_tot", ct.numst);
        
        /*  Compute and store static (nonlocal) static part of Hamiltonian */
        ReplicatedMatrix H_static("H_static", ct.numst);
        mgmol->computeHnl(&orbitals, H_static);
 
        /* Initialize observables */
        double trace = 0., avg_energy = 0.;
        double trace_0 = X_real.trace();
        ReplicatedMatrix temp("temp", ct.numst); 


        std::ofstream outputFile("outfile.txt");
        // Check if the file opened successfully
        if (!outputFile.is_open()) 
        {
            std::cerr << "Error opening logging file!" << std::endl;
            return 1;
        }
//        outputFile<<"I am open"<<std::endl;
        if (MPIdata::onpe0)
        {
            std::cout << "\nInitial State, t = 0.:" << std::endl;
            std::cout << "Initial trace: " << trace_0 << std::endl;
            /* Initial (ground state) energy before density matrix evolution */
            std::cout << "Initial KS Energy: " << eks << std::endl;

            if(LOGGING)
            {
                outputFile << std::left
                         << std::setw(15) <<"t" <<' '
                         << std::setw(12) <<"CN_iters" <<' '
                         << std::setw(18) <<"Tr(X)" <<' '
                         << std::setw(18) <<"delta_Tr(X)" <<' '
                         << std::setw(18) <<"E[XH]" <<' '
                         << std::setw(18) <<"eks"
                         << '\n';
                
                outputFile << std::right;
            }
        }
//        outputFile.flush();
        /* begin time-stepping */
        for(int step=0; step<NUM_TIME_STEPS; step++)
        {
            t = step * dt;

            /* Update Hamiltonian */
            /* Here we update only dynamic part of Hamiltonian */
            mgmol->updateHFromHnl(&orbitals, H_static, H_tot);

            /* Propagate - Crank-Nicolson propagation */
            iters = propagate_density_matrix(X_real, X_imag, H_tot, dt, 1.0e-12, 100);
            avg_iters += iters;
            
            /* symmetrize to enforce Hermiticity */
            X_real.transpose(0.5, X_real, 0.5); /* (X_real + X_real^T)/2 */
            X_imag.transpose(-0.5, X_imag, 0.5); /* X_imag - X_imag^T)/2 */

            /* Update the density matrix (Only real part is needed)*/
            dprojmatrices->setDM(X_real);
            /* Compute Obersevables:
             * 1. trace of density matrix
             * 2. Energy
            */
            /* update total Kohn-Sham energy */
            mgmol->finalEnergy();          
            /* print observables */
            if(MPIdata::onpe0 && (step % PRINT_INTERVAL == 0))
            {
                trace = X_real.trace();
                double trace_diff = std::abs(trace - trace_0);
            
                /* compute average energy (expected value of energy) */
                temp.gemm('n','n',1.0, X_real, H_tot, 0.);
                avg_energy = temp.trace();
                
                /* update total Kohn-Sham energy */
                eks = mgmol->getTotalEnergy();
                
                if(LOGGING)
                {    
                    outputFile << std::scientific << std::setprecision(4)
                         <<  std::setw(5) << (t + dt) <<' '
                         <<  std::setw(10) << iters <<' ';
                    outputFile << std::scientific << std::setprecision(8)
                         <<  std::setw(18) << trace <<' '
                         <<  std::setw(18) << trace_diff <<' '
                         <<  std::setw(18) << avg_energy <<' '
                         << std::setw(18) << eks
                         << '\n';                
                }             
            }        
        }
        
        /* Final summary */
        /* trace of Density matrix -- conserved value */
        trace = X_real.trace();
        double trace_diff = std::abs(trace - trace_0);
        /* compute average energy (expected value of energy) -- conserved value */
        temp.gemm('n','n',1.0, X_real, H_tot, 0.);
        avg_energy = temp.trace();

        /* update total Kohn-Sham energy */
        mgmol->finalEnergy();
        eks = mgmol->getTotalEnergy();
            
        // print out results
        if (MPIdata::onpe0)
        {
            std::cout<<"Final Statistics"<<std::endl;
            std::cout << "\nPropagation complete: t = " << (DELTA_T * NUM_TIME_STEPS)<<" a.u."<< std::endl;
            std::cout << "Average CN_iters: " << (avg_iters/NUM_TIME_STEPS)<< std::endl;
            std::cout << "Initial trace: " << trace_0 << std::endl;
            std::cout << "Final trace:   " << trace << std::endl;
            std::cout << "Total drift:   " << trace_diff << std::endl;
            std::cout << "E[XH] : " << avg_energy << std::endl;
            std::cout << "Final KS Energy : " << eks << std::endl;
            std::cout << "Forces :" << std::endl;
            for (std::vector<double>::iterator it = forces.begin();
                 it != forces.end(); it += 3)
            {
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
            }
        }

/* Print final aromic positions after run */
        mgmol->getAtomicPositions(positions);
        if (MPIdata::onpe0)
        {
            std::cout << "Positions:" << std::endl;
            std::vector<short>::iterator ita = anumbers.begin();
            for (std::vector<double>::iterator it = positions.begin();
                 it != positions.end(); it += 3)
            {
                std::cout << *ita;
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
                ita++;
            }
        }

        outputFile.close();

        delete mgmol;

    } // close main scope

    mgmol_finalize();

    mpirc = MPI_Finalize();
    if (mpirc != MPI_SUCCESS)
    {
        std::cerr << "MPI Finalize failed!!!" << std::endl;
    }

    time_t tt;
    time(&tt);
    if (onpe0) std::cout << " Run ended at " << ctime(&tt) << std::endl;

    return 0;
}

/* Crank-Nicolson evolution of the Density Matrix (real and imaginary parts)
 * X_real' - 0.5*dt*(H*X_imag' - X_imag'*H) =  X_real + 0.5*dt*(H*X_imag - X_imag*H)
 * X_imag' + 0.5*dt*(H*X_real' - X_real'*H) =  X_imag - 0.5*dt*(H*X_real - X_real*H)
 */
int propagate_density_matrix(ReplicatedMatrix& X_real, ReplicatedMatrix& X_imag, ReplicatedMatrix& H, double dt,  double tol, int maxits)
{
    int iters = 0;
    double alpha = 0.5 * dt;
    double m_alpha = -alpha;
    double beta = 1.;
    double error_norm = 0.;
    double error_norm_real, error_norm_imag;

    /* commutator operators */
    ReplicatedMatrix comm_HX_r("comm_HX_r", H.m());
    ReplicatedMatrix comm_HX_i("comm_HX_i", H.m());   

    /* Iteration work matrices */
    /* initialize as previous X */
    ReplicatedMatrix work_r(X_real);
    ReplicatedMatrix work_i(X_imag);     

 //   if (MPIdata::onpe0) std::cout<<"X_real ... "<<std::endl;
 //   X_real.print((*MPIdata::sout), 0, 0, 10, 10);

    /* Precompute right-hand-sides 
     * rhs_r = X_real + 0.5*dt*(H*X_imag - X_imag*H)
     * rhs_i = X_imag - 0.5*dt*(H*X_real - X_real*H)
    */
    ReplicatedMatrix rhs_r(X_real);
    ReplicatedMatrix rhs_i(X_imag);
    /* compute rhs via two gemm operations */
    rhs_r.gemm('n', 'n', alpha, H, X_imag, beta);
    rhs_i.gemm('n', 'n', alpha, X_real, H, beta);
    
    rhs_r.gemm('n', 'n', m_alpha, X_imag, H, beta);
    rhs_i.gemm('n', 'n', m_alpha, H, X_real, beta);

    /* begin Crank-Nicolson loop */
    for(iters=0; iters< maxits; iters++)
    {
        /* Compute new commutators with current guess */
        compute_commutator(H, X_real, comm_HX_r);
        compute_commutator(H, X_imag, comm_HX_i);

        /* update density matrix */
        /* update real part */
        X_real = rhs_r;
        X_real.axpy(alpha, comm_HX_i);
        /* update imaginary part */
        X_imag = rhs_i;
        X_imag.axpy(m_alpha, comm_HX_r);

        /* Check for convergence */
        work_r -= X_real;
        work_i -= X_imag;
        error_norm_real = work_r.norm('F');
        error_norm_imag = work_i.norm('F');
        
        error_norm = std::sqrt(pow(error_norm_real,2.) + pow(error_norm_imag,2.));

        if (MPIdata::onpe0) std::cout<<"iters = "<<iters<<", error norm = "<<error_norm<<std::endl;

        if(error_norm < tol)
            break;
    
        /* update work matrices for next iterate (previous solution)*/
        work_r = X_real;
        work_i = X_imag;
    }
    return iters;
}

/* Compute Commutator operator (H*X - X*H) */
void compute_commutator(const ReplicatedMatrix& H, const ReplicatedMatrix& X, ReplicatedMatrix& op)
{
    op.clear();
    op.gemm('n', 'n', 1., H, X, 0.);
    op.gemm('n', 'n', -1., X, H, 1.);
}
