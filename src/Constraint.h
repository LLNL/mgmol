// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <cassert>
#include <string>
#include <vector>

class Constraint
{
protected:
    std::vector<std::string>
        names_; // names of atoms involved in the constraint

public:
    virtual ~Constraint() {};

    virtual bool enforce(void)            = 0;
    virtual bool project_out_forces(void) = 0;
    // virtual double projected_force(void) = 0;
    virtual void setup(const std::vector<std::string>&,
        const std::vector<double>& pmass, const std::vector<short>&,
        std::vector<double*>& tau, std::vector<double*>& taup,
        std::vector<double*>& fion, const std::vector<std::string>& local_names)
        = 0;
    virtual std::string type(void) const = 0;
    std::string names(const int i) const
    {
        assert((i >= 0) && (i < (int)names_.size()));
        return names_[i];
    }

    virtual void print(std::ostream& os) const      = 0;
    virtual void printForce(std::ostream& os) const = 0;
};

#endif
