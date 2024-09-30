// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// This driver uses the header-only LBFGS++ library
// https://github.com/yixuan/LBFGSpp

#include <Eigen/Core>
#include <LBFGS.h>
#include <iostream>

#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "mgmol_run.h"

#include <cassert>
#include <iostream>
#include <time.h>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using Eigen::VectorXd;

class MGmolEnergyAndForces
{
private:
    MGmolInterface* mgmol_;
    int n_;
    std::vector<short>& anumbers_;

public:
    MGmolEnergyAndForces(
        MGmolInterface* mgmol, const int n, std::vector<short>& anumbers)
        : mgmol_(mgmol), n_(n), anumbers_(anumbers)
    {
    }

    double operator()(const VectorXd& x, VectorXd& grad)
    {
        std::vector<double> positions(n_);
        for (int i = 0; i < n_; i++)
        {
            positions[i] = x[i];
        }
        if (MPIdata::onpe0)
        {
            std::cout << "Positions:" << std::endl;
            for (std::vector<double>::iterator it = positions.begin();
                 it != positions.end(); it += 3)
            {
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
            }
        }

        // compute energy and forces using all MPI tasks
        // expect positions to be replicated on all MPI tasks
        std::vector<double> forces(n_);
        double fx
            = mgmol_->evaluateEnergyAndForces(positions, anumbers_, forces);
        // print out results
        if (MPIdata::onpe0)
        {
            std::cout << "Energy: " << fx << std::endl;
            std::cout << "Forces:" << std::endl;
            for (std::vector<double>::iterator it = forces.begin();
                 it != forces.end(); it += 3)
            {
                double norm = 0.;
                for (int i = 0; i < 3; i++)
                {
                    double val = *(it + i);
                    std::cout << "    " << val;
                    norm += val * val;
                }
                std::cout << "  norm: " << std::sqrt(norm);
                std::cout << std::endl;
            }
        }

        // set gradient to negative forces
        for (int i = 0; i < n_; i++)
        {
            grad[i] = -1. * forces[i];
        }
        return fx;
    }
};

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
            mgmol = new MGmol<LocGridOrbitals>(global_comm, *MPIdata::sout,
                input_filename, lrs_filename, constraints_filename);
        else
            mgmol = new MGmol<ExtendedGridOrbitals>(global_comm, *MPIdata::sout,
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

        // Set up parameters
        LBFGSpp::LBFGSParam<double> param;
        param.epsilon        = 4e-4;
        param.max_iterations = 100;

        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver(param);
        const int n = positions.size();
        if (MPIdata::onpe0) std::cout << "n = " << n << std::endl;
        MGmolEnergyAndForces fun(mgmol, n, anumbers);

        // initial guess
        VectorXd x = VectorXd::Zero(n);
        int i      = 0;
        for (auto& pos : positions)
        {
            x[i++] = pos;
        }

        double eks;
        int niter = solver.minimize(fun, x, eks);

        std::cout << niter << " iterations" << std::endl;

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
