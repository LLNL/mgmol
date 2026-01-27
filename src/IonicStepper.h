// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef IONICSTEPPER_H
#define IONICSTEPPER_H

#include "HDFrestart.h"
#include "MGmol_MPI.h"

#include <vector>

class IonicStepper
{
protected:
    const std::vector<short>& atmove_; // atmove_[na_]
    int ndofs_;
    double dt_;
    std::vector<double>& tau0_; // tau0_[3*na_]
    std::vector<double>& taup_; // taup_[3*na_]

public:
    IonicStepper(const double dt, const std::vector<short>& atmove,
        std::vector<double>& tau0, std::vector<double>& taup);

    int writeAtomicFields(HDFrestart& h5f_file, const std::vector<double>&,
        const std::string& name, const bool) const;
    int writePositions(HDFrestart& h5f_file,
        const std::string& name = "/Ionic_positions") const;
    int readPositions_hdf5(HDFrestart& h5f_file, const std::string& name);
    int readRandomStates(HDFrestart& h5f_file,
        std::vector<unsigned short>& data,
        const std::string& name = "/Ionic_RandomStates");
    int writeRandomStates(HDFrestart& h5f_file,
        std::vector<unsigned short>& data,
        const std::string& name = "/Ionic_RandomStates") const;
    int writeVelocities(HDFrestart& h5f_file) const;
    int readAtomicFields(
        HDFrestart& h5f_file, std::vector<double>&, const std::string& name);

    virtual bool check_last_step_accepted() const { return true; }
    virtual int run()                   = 0;
    virtual double etol(void) const     = 0;
    virtual int write_hdf5(HDFrestart&) = 0;
    virtual int init(HDFrestart&)       = 0;

    virtual ~IonicStepper() {};
};
#endif
