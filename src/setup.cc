// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ConstraintSet.h"
#include "ExtendedGridOrbitals.h"
#include "Hamiltonian.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "Potentials.h"

template <class T>
int MGmol<T>::setupFromInput(const std::string filename)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("MGmol<T>::setupFromInput()...", std::cout);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    hamiltonian_    = new Hamiltonian<T>();
    Potentials& pot = hamiltonian_->potential();

    ct.registerPotentials(pot);
    ct.setSpecies(pot);

    if (ct.diel) pot.turnOnDiel();

    ct.adjust();

    Mesh* mymesh = Mesh::instance();
    if (ct.isLocMode()) mymesh->subdivGridx(ct.getMGlevels());

    const pb::PEenv& myPEenv = mymesh->peenv();
    if (ct.restart_info > 0)
        h5f_file_
            = new HDFrestart(ct.restart_file, myPEenv, ct.restart_file_type);

    int status = readCoordinates(filename, false);
    if (status == -1) return -1;

    const short myspin = mmpi.myspin();
    const int nval     = ions_->getNValenceElectrons();
    ct.setNumst(myspin, nval);
    ct.setTolEnergy();
    ct.setSpreadRadius();

    // create localization regions
    const pb::Grid& mygrid = mymesh->grid();
    Vector3D vcell(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));
    lrs_ = new LocalizationRegions(vcell, ct.tol_orb_centers_move);

    return 0;
}

template <class T>
int MGmol<T>::setupLRsFromInput(const std::string filename)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    //    Mesh* mymesh           = Mesh::instance();
    //    const pb::Grid& mygrid = mymesh->grid();
    //    Control& ct            = *(Control::instance());

    std::ifstream* tfile = nullptr;
    if (mmpi.instancePE0() && !filename.empty())
    {
        os_ << "Read LRs from file " << filename << std::endl;
        tfile = new std::ifstream(filename.data(), std::ios::in);
        if (!tfile->is_open())
        {
            std::cerr << " Unable to open file " << filename << std::endl;
            global_exit(0);
        }
        else
        {
            os_ << "Open " << filename << std::endl;
        }
    }

    readLRsFromInput(tfile);

    if (!(tfile == nullptr))
    {
        os_ << "Close " << filename << std::endl;
        tfile->close();
        delete tfile;
    }

    return 0;
}

template <class T>
int MGmol<T>::setupConstraintsFromInput(const std::string filename)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    std::ifstream* tfile = nullptr;
    if (mmpi.instancePE0() && !filename.empty())
    {
        os_ << "Read constraints from file " << filename << std::endl;
        tfile = new std::ifstream(filename.data(), std::ios::in);
        if (!tfile->is_open())
        {
            std::cerr << " Unable to open file " << filename << std::endl;
            global_exit(0);
        }
        else
        {
            os_ << "Open " << filename << std::endl;
        }
    }

    constraints_->readConstraints(tfile);

    if (!(tfile == nullptr))
    {
        os_ << "Close " << filename.data() << std::endl;
        tfile->close();
        delete tfile;
    }

    return 0;
}

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;
