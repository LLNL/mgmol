// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

//
//    Functions to read and write restart data to files.
//
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Control.h"
#include "HDFrestart.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MGmol.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "Potentials.h"
#include "ProjectedMatrices.h"
#include "Rho.h"
#include "hdf5.h"
#include "tools.h"

// read rho and potentials form a hdf5 file
template <class T>
int MGmol<T>::read_rho_and_pot_hdf5(HDFrestart& file, Rho<T>& rho)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        os_ << "Try to read density and potentials" << endl;

    Potentials& pot = hamiltonian_->potential();

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    POTDTYPE* tmp          = new POTDTYPE[mygrid.size()];

    // Read total potential
    file.read_1func_hdf5(pot.vtot(), "Vtotal");

    // Read the hartree potential
    file.read_1func_hdf5(tmp, "Hartree");
    pot.setVh(tmp, 0);

    // Read dielectric potential
    if (ct.diel)
    {
        file.read_1func_hdf5(pot.vepsilon(), "VDielectric");
    }
    // Read the Density
    rho.readRestart(file);

    delete[] tmp;
    return 0;
}

// Writes restart information in a file.
template <class T>
int MGmol<T>::write_hdf5(const string filename, vector<vector<RHODTYPE>>& rho,
    Ions& ions, T& orbitals, LocalizationRegions& lrs)
{
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    const pb::Grid& mygrid   = mymesh->grid();
    Control& ct              = *(Control::instance());

    char zero = '0';
    if (strcmp(&zero, &filename[0]) == 0) return 0;

    unsigned gdim[3] = { mygrid.gdim(0), mygrid.gdim(1), mygrid.gdim(2) };

    // create restart file
    HDFrestart h5f_file(filename, myPEenv, gdim, ct.out_restart_file_type);

    int status = write_hdf5(h5f_file, rho, ions, orbitals, *lrs_);
    if (status < 0 && onpe0)
        (*MPIdata::serr) << "restart.cc: write_hdf5 failed!!!" << endl;

    return status;
}

// Writes restart information in a HDF5 file.
template <class T>
int MGmol<T>::write_hdf5(HDFrestart& h5f_file, vector<vector<RHODTYPE>>& rho,
    Ions& ions, T& orbitals, LocalizationRegions& lrs)
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    MGmol_MPI& mmpi        = *(MGmol_MPI::instance());

    Potentials& pot = hamiltonian_->potential();

    Control& ct = *(Control::instance());

    Timer timer("WriteRestart");
    timer.start();

    if (ct.out_restart_info > 0)
    {
        ions.writeAtomNames(h5f_file);
        ions.writeLockedAtomNames(h5f_file);
        ions.writeAtomicNumbers(h5f_file);
        ions.writeAtomicIDs(h5f_file);
        ions.writeAtomicNLprojIDs(h5f_file);
        ions.writePositions(h5f_file);
        ions.writeRandomStates(h5f_file);
        ions.writeVelocities(h5f_file);
        ions.writeForces(h5f_file);

        h5f_file.addMDTime2File(md_time_);
        h5f_file.addMDstep2File(md_iteration_);
    }

    double ll[3] = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };

    double origin[3] = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };

    if (ct.out_restart_info > 1)
    {
        // Write total potential
        int ierr = h5f_file.write_1func_hdf5(
            pot.vtot(), "Vtotal", &ll[0], &origin[0]);
        if (ierr < 0) return ierr;

        // Write the hartree potential
        ierr = h5f_file.write_1func_hdf5(
            pot.vh_rho(), "Hartree", &ll[0], &origin[0]);
        if (ierr < 0) return ierr;

        // Write
        if (ct.diel)
        {
            ierr = h5f_file.write_1func_hdf5(
                pot.vepsilon(), "VDielectric", &ll[0], &origin[0]);
        }
        if (ierr < 0) return ierr;

        // Write the Density
        ierr = h5f_file.write_1func_hdf5(
            &rho[0][0], "Density", &ll[0], &origin[0]);
        if (ierr < 0) return ierr;

        // Write external potential
        ierr
            = h5f_file.write_1func_hdf5(pot.vext(), "Vext", &ll[0], &origin[0]);
        if (ierr < 0) return ierr;
    }

    // Write wavefunctions and old centers.
    if (ct.out_restart_info > 2)
    {
        int ierr = orbitals.write_hdf5(h5f_file);
        if (ierr < 0) return ierr;

        if (ct.isLocMode() &&
            ct.WFExtrapolation() == WFExtrapolationType::Reversible)
        {
            lrs.writeOldCenters(h5f_file);
        }
    }

    mmpi.barrier();

    timer.stop();
    if (onpe0)
    {
        os_ << "Wrote restart data --- timing: " << endl;
    }
    timer.print(os_);
    return 0;
}

template <class T>
int MGmol<T>::read_restart_lrs(HDFrestart& h5f_file, const string& dset_name)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("MGmol<T>::read_restart_lrs()...", (*MPIdata::sout));

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    int n = 0;
    if (ct.restart_info >= 3)
    {
        n = h5f_file.getLRs(*lrs_, ct.numst, dset_name);
    }
    mmpi.bcast(&n, 1);

    return n;
}

// Reads the restart information from restart files.
template <class T>
int MGmol<T>::read_restart_data(
    HDFrestart& h5f_file, Rho<T>& rho, T& orbitals)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("MGmol<T>::read_restart_data()...", (*MPIdata::sout));

    Timer timer("ReadRestart");
    timer.start();

    if (ct.restart_info > 1)
    {
        read_rho_and_pot_hdf5(h5f_file, rho);
    }

    if (ct.restart_info > 2)
    {
        int ierr = orbitals.read_hdf5(h5f_file);
        if (ierr < 0)
        {
            (*MPIdata::serr) << "MGmol<T>::read_restart_data(): error in reading "
                             << endl;
            return ierr;
        }
    }

    timer.stop();

    if (onpe0) os_ << "Read restart data --- timing: " << endl;
    timer.print(os_);

    return 0;
}

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;
