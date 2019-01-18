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
// MD_IonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id$

#ifndef MD_IONICSTEPPER_H
#define MD_IONICSTEPPER_H

#include "IonicStepper.h"
#include <cassert>

class MD_IonicStepper : public IonicStepper
{
private:
    // taum_ stores velocities if dt_<=0
    std::vector<double>& taum_; // taum_[3*na_]

    // mass of each atom
    const std::vector<double>& pmass_;

    std::vector<double>& fion_;

    // thermostat type
    // 0=off
    // 1=Berendsen (Berendsen et al., J. Chem. Phys. 81 (8), 1984)
    // 2=Langevin (Brunger, Brooks, Karplus, Chem. Phys. Lett. 105, 495, 1982)
    // 3=SCALING (according to Qbox implementation)
    int ttherm_;

    // target temperature
    double tkel_;

    double thtime_; // temperature control
    double gamma_; // damping coefficient

    // Temperature control for SCALING algorithm
    double thwidth_;

    std::vector<unsigned short>& rand_states_;

    double getTemperatureAtTimeMinusHalf();
    double getTemperatureAtHalfStep();

    void updateWithVelocityScaling(const double lambda);

public:
    MD_IonicStepper(const double dt, std::vector<short>& atmove,
        std::vector<double>& tau0, std::vector<double>& taup,
        std::vector<double>& taum, std::vector<double>& fion,
        const std::vector<double>& pmass,
        std::vector<unsigned short>& rand_states);

    void setThermostat(const int ttherm, const double tkel, const double thtime,
        const double thwidth, const int nconstraints);

    void addLangevinForces();

    double temperature();
    double kineticEnergy();
    int run();
    void updateTau();
    double etol(void) const;
    int write_hdf5(HDFrestart&);
    int writeForces(HDFrestart& h5f_file);
    int writeTaum(HDFrestart& h5f_file);
    int init(HDFrestart&);
    void printVelocities(std::ostream& os) const;

    void getVelocities(std::vector<float>& velocities) const;
};

#endif
