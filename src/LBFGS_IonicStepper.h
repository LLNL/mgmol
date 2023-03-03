// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifndef LBFGS_IONICSTEPPER_H
#define LBFGS_IONICSTEPPER_H

#include "MPIdata.h"

#include "IonicStepper.h"
#include "hdf5.h"

#include <cassert>
#include <iostream>
#include <vector>

class LBFGS_IonicStepper : public IonicStepper
{
private:
    std::vector<double>& fion_;
    const std::vector<int>& gids_;

    // working arrays (duplicated on every processor)
    std::vector<double> xcurrent_;
    std::vector<double> g_;

    std::vector<double> work_;
    std::vector<double> xref_;
    std::vector<double> diag_;
    int iflag_;

    // global number of variables
    int gndofs_;
    int m_;
    double* etot_;
    int iter_;
    int point_;
    double stp_;
    int info_;
    int infoc_;
    int nfev_;
    int npt_;

    double finit_;
    double dginit_;
    bool brackt_;
    bool stage1_;

    bool freeze_geom_center_;

    //  the variables stx, fx, dgx contain the values of the step,
    //  function, and directional derivative at the best step.
    //  the variables sty, fy, dgy contain the value of the step,
    //  function, and derivative at the other endpoint of
    //  the interval of uncertainty.
    double fx_;
    double fy_;
    double stx_;
    double sty_;
    double dgx_;
    double dgy_;

    double stmin_;
    double stmax_;

    double width_;
    double width1_;

    double eps_; // tolerance on gradient
    // Termination occurs when the relative width of the interval
    // of uncertainty is at most xtol.
    double xtol_;

    int last_step_accepted_;

    int lbfgs(const double, const bool, const int);
    void minushg(const int, const double);
    void mcsrch(double f, double* s, double& stp);
    int writeDoubleAtt(hid_t dataset_id) const;
    int writeIntAtt(hid_t dataset_id) const;
    int read_lbfgs(HDFrestart&);
    void setDataFromLocalData(
        std::vector<double>& x, const std::vector<double>& tau);

public:
    LBFGS_IonicStepper(const double dt, const std::vector<short>& atmove,
        std::vector<double>& tau0, std::vector<double>& taup,
        std::vector<double>& fion, const std::vector<int>& gid, const int m,
        double* etot);

    int run() override;
    double etol(void) const override;
    int write_hdf5(HDFrestart&) override;
    int init(HDFrestart&) override;
    void restart();
    int writeLBFGSinfo(HDFrestart&);

    bool check_last_step_accepted() const override
    {
        return (bool)last_step_accepted_;
    }
};

#endif
