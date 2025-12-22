// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "Control.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MGmol.h"
#include "Mesh.h"
#include "Species.h"
#include "tools.h"

#include <cassert>
#include <iostream>
#include <string>

#include <mpi.h>

// #define DEBUG 1

template <class OrbitalsType>
int MGmol<OrbitalsType>::readLRsFromInput(std::ifstream* tfile)
{
    Control& ct(*(Control::instance()));

    assert(lrs_);
    assert(ct.restart_info < 3 || !ct.isLocMode());

    if (ct.verbose > 0) printWithTimeStamp("readLRsFromInput", os_);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const double lattice[3] = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };
    const double origin[3]
        = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };
    const double end[3] = { mygrid.origin(0) + mygrid.ll(0),
        mygrid.origin(1) + mygrid.ll(1), mygrid.origin(2) + mygrid.ll(2) };

    const bool read_radius = ct.adaptiveLRsizes() ? true : false;

    if (mmpi.instancePE0())
    {
        if (ct.verbose > 0)
            os_ << "Trying to read " << ct.numst
                << " orbital centers from input file..." << std::endl;
        const std::vector<Ion*>& list_ions(ions_->list_ions());
        std::vector<Ion*>::const_iterator ion = list_ions.begin();

        if (tfile != nullptr) read_comments(*tfile);
        int count = 0;

        // read as many centers as there are functions
        for (int i = 0; i < ct.numst; i++)
        {
            double crds[3] = { 0., 0., 0. };
            bool flag      = false;

            if (tfile != nullptr)
                for (short j = 0; j < 3; j++)
                {
                    std::string sread;
                    (*tfile) >> sread;
                    if (tfile->fail())
                    {
                        os_ << "WARNING: Failed reading Localization center... "
                            << std::endl;
                        flag = false;
                        break;
                    }
                    else
                    {
                        flag    = true;
                        crds[j] = atof(sread.c_str());
                    }
                }

            // set center to ionic position if not read from input file
            if (!flag)
            {
                if (ct.verbose > 1)
                    os_ << "Use atomic position for center " << i << std::endl;

                for (int k = 0; k < 3; k++)
                    crds[k] = (*ion)->position(k);
                ion++;
                if (ion == list_ions.end()) ion = list_ions.begin();
            }

            // move coordinates inside domain boundaries
            for (int j = 0; j < 3; j++)
            {
                while (crds[j] < origin[j])
                    crds[j] += lattice[j];
                while (crds[j] >= end[j])
                    crds[j] -= lattice[j];

                assert(crds[j] > -1000.);
                assert(crds[j] < 1000.);
            }

            float radius = ct.cut_radius;
            if (flag && ct.isLocMode() && read_radius)
            {
                float tmp;
                (*tfile) >> tmp;
                if (tfile->fail())
                {
                    os_ << "WARNING: Failed reading Localization radius!!! Use "
                           "uniform radius..."
                        << std::endl;
                }
                else
                {
                    radius = tmp;
                }
            }

            // this is where a new function with a new gid is created
            Vector3D tmp_center(crds[0], crds[1], crds[2]);
            lrs_->push_back_global(tmp_center, radius);
            if (ct.verbose > 2)
                (*MPIdata::sout) << "Added LR with center " << tmp_center
                                 << " and radius " << radius << std::endl;
            count++;

            // finish reading line
            if (flag)
                while (tfile->get() != '\n')
                    ;
#ifdef DEBUG
            os_ << " Read orbital center (" << crds[0] << "," << crds[1] << ","
                << crds[2] << ")" << std::endl;
#endif
        }

        if (ct.verbose > 1 && onpe0) os_ << "Randomize gids..." << std::endl;
        lrs_->randomizeGids();

        lrs_->fillWithZeroCenters(ct.numst);

        if (onpe0 && ct.verbose > 0)
            os_ << "readInput: Read " << count
                << " orbital centers from input file" << std::endl;
    } // mmpi.instancePE0()

    if (ct.verbose > 0) printWithTimeStamp("setup LRs...", os_);
    lrs_->setup();
    lrs_->printInfo(os_);

#ifdef DEBUG
    lrs_->print(os_);
#endif

    return ct.numst;
}

template <class OrbitalsType>
int MGmol<OrbitalsType>::readCoordinates(
    const std::string& coords_filename, const bool cell_relative)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0) printWithTimeStamp("Read atomic coordinates...", os_);
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const double lattice[3] = { mygrid.ll(0), mygrid.ll(1), mygrid.ll(2) };

    // setup ions
    const std::vector<Species>& sp(ct.getSpecies());
    ions_.reset(new Ions(lattice, sp));

    if (ct.restart_info > 0 && coords_filename.empty())
    {
        // read restart atomic positions
        if (onpe0 && ct.verbose > 0)
        {
            os_ << "Initialize atomic positions from restart file "
                << ct.restart_file << std::endl;
        }
        ions_->initFromRestartFile(*h5f_file_);
    }
    else
    {
        // Coordinates and species type for each ion.
        int info = ions_->readAtoms(coords_filename, cell_relative);

        return info;
    }

    const int num_ions = ions_->getNumIons();
    if (onpe0) os_ << num_ions << " atoms in simulation" << std::endl;

    return 0;
}

template class MGmol<LocGridOrbitals<ORBDTYPE>>;
template class MGmol<ExtendedGridOrbitals<ORBDTYPE>>;
