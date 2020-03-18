// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$

#ifndef DISTANCECONSTRAINT_H
#define DISTANCECONSTRAINT_H

#include "Constraint.h"
#include <cassert>

class DistanceConstraint : public Constraint
{
    const std::string name1_;
    const std::string name2_;
    const double distance_;

    double tol_;

    double* tau1_;
    double* tau1p_;
    double m1_;
    bool locked1_;

    double* tau2_;
    double* tau2p_;
    double m2_;
    bool locked2_;

    double f1_, f2_; // mu/m_i

    double* fion1_;
    double* fion2_;

    // is constraint active on set of local ions?
    bool locally_active_;

    // is this task responsible for this constraint?
    // (like printing it)
    bool locally_owned_;

    double projected_force(void) const;
    std::ostream& printConstraint(std::ostream& os) const;

public:
    DistanceConstraint(const std::string& name1, const std::string& name2,
        const double distance, const double tolerance)
        : name1_(name1),
          name2_(name2),
          distance_(distance),
          tol_(tolerance),
          m1_(0.0),
          locked1_(false),
          m2_(0.0),
          locked2_(false)
    {
        names_.push_back(name1);
        names_.push_back(name2);
        tau1_ = tau1p_ = tau2_ = tau2p_ = fion1_ = fion2_ = nullptr;
        locally_active_                                   = false;
        locally_owned_                                    = false;
    };

    bool enforce(void) override;
    bool project_out_forces(void) override;

    void setup(const std::vector<std::string>& names,
        const std::vector<double>& pmass, const std::vector<short>& imove,
        std::vector<double*>& tau, std::vector<double*>& taup,
        std::vector<double*>& fion,
        const std::vector<std::string>& local_names) override;

    std::string type(void) const override { return "distance"; };
    double distance(void) const { return distance_; }
    void print(std::ostream& os) const override;
    void printForce(std::ostream& os) const override;
};
#endif
