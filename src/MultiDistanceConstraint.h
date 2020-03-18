// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$

#ifndef MULTIDISTANCECONSTRAINT_H
#define MULTIDISTANCECONSTRAINT_H

#include "Constraint.h"

#include <cassert>
#include <vector>

class MultiDistanceConstraint : public Constraint
{
    // coefficients
    std::vector<double> m_alpha_;

    // two names associated with each coefficient
    std::vector<std::string> m_name1_;
    std::vector<std::string> m_name2_;

    // value of the constraint (alpha1*d1+alpha2*d2+...)
    double distance_;
    double tol_;

    std::vector<double*> m_tau1_;
    std::vector<double*> m_tau1p_;
    std::vector<double*> u_taup_;
    std::vector<double*> m_tau2_;
    std::vector<double*> m_tau2p_;
    std::vector<double*> m_fion1_;
    std::vector<double*> m_fion2_;
    std::vector<double*> u_fion_;

    // list of unique names
    std::vector<std::string> u_name_;
    std::vector<double> u_mass_;
    std::vector<bool> u_locked_;

    std::vector<double> m_f1_, m_f2_;

    // is constraint active on set of local ions?
    bool locally_active_;

    // is this task responsible for this constraint?
    // (like printing it)
    // there should be one and only one
    bool locally_owned_;

    double projected_force(void) const;
    std::ostream& printConstraint(std::ostream& os) const;

public:
    MultiDistanceConstraint(const std::vector<double>& alpha,
        const std::vector<std::string>& name1,
        const std::vector<std::string>& name2, const double distance,
        const double tolerance)
        : m_alpha_(alpha),
          m_name1_(name1),
          m_name2_(name2),
          distance_(distance),
          tol_(tolerance)
    {
        for (int i = 0; i < (int)alpha.size(); i++)
        {
            names_.push_back(name1[i]);
            names_.push_back(name2[i]);
        }
        locally_active_ = false;
        locally_owned_  = false;
    }
    bool enforce(void) override;
    bool project_out_forces(void) override;

    void setup(const std::vector<std::string>& names,
        const std::vector<double>& pmass, const std::vector<short>& imove,
        std::vector<double*>& tau, std::vector<double*>& taup,
        std::vector<double*>& fion,
        const std::vector<std::string>& local_names) override;

    std::string type(void) const override { return "multidistance"; };
    double distance(void) const { return distance_; };
    void print(std::ostream& os) const override;
    void printForce(std::ostream& os) const override;
};
#endif
