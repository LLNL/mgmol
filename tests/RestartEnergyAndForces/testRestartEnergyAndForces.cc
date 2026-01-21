// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

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

        MGmolInterface* mgmol = new MGmol<ExtendedGridOrbitals<ORBDTYPE>>(
            global_comm, *MPIdata::sout, input_filename, lrs_filename,
            constraints_filename);

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

        Mesh* mymesh             = Mesh::instance();
        const pb::Grid& mygrid   = mymesh->grid();
        const pb::PEenv& myPEenv = mymesh->peenv();

        HDFrestart h5file(ct.restart_file, myPEenv, ct.restart_file_type);
        std::string name = "Function";
        int count        = h5file.countFunctionObjects(name);
        if (count != ct.numst)
        {
            std::cerr << "The number of functions in the restart file, "
                      << count << " is not equal to ct.numst, " << ct.numst
                      << std::endl;
            MPI_Abort(mmpi.commSameSpin(), 0);
        }

        std::shared_ptr<ProjectedMatricesInterface> projmatrices
            = mgmol->getProjectedMatrices();

        ExtendedGridOrbitals<ORBDTYPE> orbitals("new_orbitals", mygrid,
            mymesh->subdivx(), ct.numst, ct.bcWF, projmatrices.get(), nullptr,
            nullptr, nullptr, nullptr);

        // read numst_ wavefunction
        int nread = orbitals.read_func_hdf5(h5file, name);
        if (nread != ct.numst)
        {
            std::cerr << "The number of functions read from the restart file, "
                      << nread << " is not equal to ct.numst, " << ct.numst
                      << std::endl;
            MPI_Abort(mmpi.commSameSpin(), 0);
        }

        // set the iterative index to 1 to differentiate it from first instance
        // in MGmol initial() function. This is not very clean and could be
        // better designed, but works for now
        orbitals.setIterativeIndex(1);

        // set initial DM with uniform occupations
        projmatrices->setDMuniform(ct.getNelSpin());
        projmatrices->printDM(std::cout);

        // swap H and O to make sure order of atoms in list does not matter
        double x     = positions[0];
        double y     = positions[1];
        double z     = positions[2];
        positions[0] = positions[3];
        positions[1] = positions[4];
        positions[2] = positions[5];
        positions[3] = x;
        positions[4] = y;
        positions[5] = z;
        short tmp    = anumbers[0];
        anumbers[0]  = anumbers[1];
        anumbers[1]  = tmp;
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

        //
        // evaluate energy and forces with wavefunctions just read
        //
        std::vector<double> forces;
        double eks = mgmol->evaluateDMandEnergyAndForces(
            &orbitals, positions, anumbers, forces);

        // print out results
        if (MPIdata::onpe0)
        {
            std::cout << "Eks : " << eks << std::endl;
            std::cout << "Forces :" << std::endl;
            for (std::vector<double>::iterator it = forces.begin();
                 it != forces.end(); it += 3)
            {
                for (int i = 0; i < 3; i++)
                    std::cout << "    " << *(it + i);
                std::cout << std::endl;
            }
        }

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
