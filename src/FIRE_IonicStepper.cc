// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

////////////////////////////////////////////////////////////////////////////////
//
// FIRE_IonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id:$

#include "FIRE_IonicStepper.h"
#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "hdf5.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <string>
using namespace std;

// const double scmass=1822.89;

// kb_au=8.617343e-5 [eV/K] / 27.211608 [eV/Ha]
// const double kb_au=3.16678939e-06; // [Ha/K]

FIRE_IonicStepper::FIRE_IonicStepper(const double dt,
    const vector<short>& atmove, vector<double>& tau0, vector<double>& taup,
    vector<double>& fion, const vector<double>& pmass)
    : IonicStepper(dt, atmove, tau0, taup), pmass_(pmass), fion_(fion)
{
    assert(3 * pmass.size() == tau0.size());
    assert(3 * atmove.size() == tau0.size());
    assert(taup.size() == tau0.size());

    nmin_ = 5;
    npp_  = 0;

    alpha_init_ = 0.1;
    alpha_      = alpha_init_;
    falpha_     = 0.99;
    finc_       = 1.1;
    fdec_       = 0.5;
    dt_         = dt;
    dtmax_      = 100.;

    for (int i = 0; i < (int)pmass.size(); i++)
        assert(pmass[i] > 0.);
};

int FIRE_IonicStepper::init(HDFrestart& h5f_file)
{
    assert(tau0_.size() > 0);
    assert(taup_.size() > 0);
    assert(taup_.size() == tau0_.size());

    hid_t file_id            = h5f_file.file_id();
    const string error_msg   = "Error: FIRE_IonicStepper::init(): ";
    const string warning_msg = "Warning: FIRE_IonicStepper::init(): ";

    // Open the dataset
    if (file_id >= 0 && onpe0)
    {
        (*MPIdata::sout) << "Initialize FIRE_IonicStepper with data from "
                         << h5f_file.filename() << endl;

        // read positions into tau0_
        string string_name("/Ionic_positions");
        readPositions_hdf5(h5f_file, string_name);

        // Read velocities equal to (taup-tau0)/dt
        hid_t dataset_id = H5Dopen2(file_id, "/Ionic_velocities", H5P_DEFAULT);
        if (dataset_id < 0)
        {
            if (onpe0)
            {
                (*MPIdata::sout)
                    << warning_msg << "H5Dopen failed for /Ionic_velocities"
                    << endl;
                (*MPIdata::sout) << "Set velocities to zero" << endl;
            }
            memset(&taup_[0], 0, taup_.size() * sizeof(double));
        }
        else
        {
            if (onpe0)
                (*MPIdata::sout) << "Read Ionic velocities from "
                                 << h5f_file.filename() << endl;
            herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                H5S_ALL, H5P_DEFAULT, &taup_[0]);
            if (status < 0)
            {
                (*MPIdata::serr) << error_msg << "H5Dread failed!!!" << endl;
                return -1;
            }
            // close dataset
            status = H5Dclose(dataset_id);
            if (status < 0)
            {
                (*MPIdata::serr) << "H5Dclose failed!!!" << endl;
                return -1;
            }
            if (taup_[0] != taup_[0])
            {
                (*MPIdata::serr)
                    << error_msg << "taup_[0]=" << taup_[0] << endl;
                return -1;
            }
        }
    }

    int n = (int)tau0_.size(), ione = 1;
#ifdef USE_MPI
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(&tau0_[0], n);
    mmpi.bcast(&taup_[0], n);
#endif
    DAXPY(&n, &dt_, &taup_[0], &ione, &tau0_[0], &ione);

    return 0;
}

int FIRE_IonicStepper::write_hdf5(HDFrestart& h5f_file)
{
    hid_t file_id = h5f_file.file_id();
    if (file_id < 0) return 0;

    const bool write_flag = (onpe0 || h5f_file.useHdf5p());

    if (write_flag)
    {
        writePositions(h5f_file);

        writeVelocities(h5f_file);
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////

int FIRE_IonicStepper::run()
{
    //(*MPIdata::sout)<<"FIRE_IonicStepper::run()"<<endl;
    const int na = (int)atmove_.size();
    assert((int)fion_.size() == 3 * na);

    if (taum_.empty()) taum_ = tau0_;

    double params[3] = { 0., 0., 0. };
    double& pp(params[0]);
    double& normf2(params[1]);
    double& normv2(params[2]);

    for (int ia = 0; ia < na; ia++)
    {
        if (atmove_[ia])
        {
            for (int j = 0; j < 3; j++)
            {
                pp += (tau0_[3 * ia + j] - taum_[3 * ia + j])
                      * fion_[3 * ia + j];
                normf2 += fion_[3 * ia + j] * fion_[3 * ia + j];
                normv2 += (tau0_[3 * ia + j] - taum_[3 * ia + j])
                          * (tau0_[3 * ia + j] - taum_[3 * ia + j]);
            }
        }
    }
    assert(pp == pp);
    assert(normf2 == normf2);
    assert(normv2 == normv2);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&params[0], 3, MPI_SUM);

    double inv_norm_f = 1. / sqrt(normf2);
    double normv      = sqrt(normv2);

    if (pp <= 0.)
    {
        dt_ *= fdec_;
        npp_   = 0;
        alpha_ = alpha_init_;
        // freeze system
        if (onpe0)
            (*MPIdata::sout) << "FIRE_IonicStepper: freeze system..." << endl;
        if (onpe0)
            (*MPIdata::sout) << "FIRE_IonicStepper: new dt   =" << dt_ << endl;
        for (int ia = 0; ia < na; ia++)
        {
            for (int j = 0; j < 3; j++)
            {
                taup_[3 * ia + j] = 0.;
            }
        }
    }
    else
    {
        // calculate new velocity
        for (int ia = 0; ia < na; ia++)
        {
            if (atmove_[ia])
            {
                for (int j = 0; j < 3; j++)
                {
                    taup_[3 * ia + j]
                        = (1. - alpha_)
                              * (tau0_[3 * ia + j] - taum_[3 * ia + j])
                          + alpha_ * fion_[3 * ia + j] * inv_norm_f * normv;
                }
            }
        }

        npp_++;
        if (npp_ > nmin_)
        {
            dt_ = min(dt_ * finc_, dtmax_);
            alpha_ *= falpha_;
            if (onpe0)
                (*MPIdata::sout)
                    << "FIRE_IonicStepper: new dt   =" << dt_ << endl;
            if (onpe0)
                (*MPIdata::sout)
                    << "FIRE_IonicStepper: new alpha=" << alpha_ << endl;
        }
    }

    for (int ia = 0; ia < na; ia++)
    {
        if (atmove_[ia])
        {
            // const double invmass=1./ pmass_[ia];
            const double invmass = 1.;
            for (int j = 0; j < 3; j++)
            {
                taup_[3 * ia + j]
                    = tau0_[3 * ia + j] + dt_ * taup_[3 * ia + j]
                      + 0.5 * dt_ * dt_ * invmass * fion_[3 * ia + j];
            }
        }
        else
        {
            for (int j = 0; j < 3; j++)
            {
                taup_[3 * ia + j] = tau0_[3 * ia + j];
            }
        }
    }

    taum_ = tau0_;

    return 0;
}

double FIRE_IonicStepper::etol(void) const { return 0.0; }
