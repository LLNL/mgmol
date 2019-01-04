// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ConstraintSet.h"
#include "Hamiltonian.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol.h"
#include "Potentials.h"

template <class T>
int MGmol<T>::setupFromInput(const string filename)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0) printWithTimeStamp("MGmol<T>::setupFromInput()...", cout);

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
int MGmol<T>::setupLRsFromInput(const string filename)
{
    MGmol_MPI& mmpi        = *(MGmol_MPI::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    Control& ct            = *(Control::instance());

    ifstream* tfile = 0;
    if (mmpi.instancePE0() && !filename.empty())
    {
        os_ << "Read LRs from file " << filename << endl;
        tfile = new ifstream(filename.data(), ios::in);
        if (!tfile->is_open())
        {
            cerr << " Unable to open file " << filename << endl;
            global_exit(0);
        }
        else
        {
            os_ << "Open " << filename << endl;
        }
    }

    readLRsFromInput(tfile);

    if (!(tfile == 0))
    {
        os_ << "Close " << filename << endl;
        tfile->close();
        delete tfile;
    }

    return 0;
}

template <class T>
int MGmol<T>::setupConstraintsFromInput(const string filename)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    ifstream* tfile = 0;
    if (mmpi.instancePE0() && !filename.empty())
    {
        os_ << "Read constraints from file " << filename << endl;
        tfile = new ifstream(filename.data(), ios::in);
        if (!tfile->is_open())
        {
            cerr << " Unable to open file " << filename << endl;
            global_exit(0);
        }
        else
        {
            os_ << "Open " << filename << endl;
        }
    }

    constraints_->readConstraints(tfile);

    if (!(tfile == 0))
    {
        os_ << "Close " << filename.data() << endl;
        tfile->close();
        delete tfile;
    }

    return 0;
}

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;

