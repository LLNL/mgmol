// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$

#ifndef CONSTRAINTSET_H
#define CONSTRAINTSET_H

#include <iostream>
#include <string>
#include <vector>

class Ions;
class Constraint;

class ConstraintSet
{
private:
    std::vector<Constraint*> vconstraints_;

    // constraints read from PE 0
    std::vector<std::vector<std::string>> msav_;

    void setup(const std::vector<std::string>&,
        const std::vector<double>& pmass, const std::vector<short>&,
        std::vector<double*>& tau, std::vector<double*>& taup,
        std::vector<double*>& fion,
        const std::vector<std::string>& local_names);
    bool addConstraint(Ions&, const std::vector<std::string>& argv);

public:
    ConstraintSet() {};

    void clear();

    void printConstraints(std::ostream& os) const;
    void printConstraintsForces(std::ostream& os);
    int size(void) const { return vconstraints_.size(); }
    void enforceConstraints(const int maxiter);
    void updateConstraints(const double dt);
    void projectOutForces(const int maxiter = 50);
    void setup(Ions&);

    // read all constraints from PE 0
    int readConstraints(const std::string& filename);
    int readConstraints(std::ifstream* tfile);

    // add all constraints to the list of constraints,
    // even if atoms invloved are not known locally
    int addConstraints(Ions& ions);
};
#endif
