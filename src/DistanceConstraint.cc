// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "DistanceConstraint.h"
#include "MPIdata.h"
#include "tools.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void DistanceConstraint::setup(const vector<string>& names,
    const vector<double>& pmass, const vector<short>& atmove,
    vector<double*>& tau, vector<double*>& taup, vector<double*>& fion,
    const vector<string>& local_names)
{
    assert(names.size() == pmass.size());
    assert(names.size() == tau.size());
    assert(names.size() == fion.size());

    short found = 0;

    // find position in tau array corresponding to atom name1
    for (int ia = 0; ia < (int)names.size(); ia++)
    {
        if (name1_.compare(names[ia]) == 0)
        {
            tau1_    = tau[ia];
            tau1p_   = taup[ia];
            fion1_   = fion[ia];
            m1_      = pmass[ia];
            locked1_ = (atmove[ia] != 1);
            found++;
        }
        if (name2_.compare(names[ia]) == 0)
        {
            tau2_    = tau[ia];
            tau2p_   = taup[ia];
            fion2_   = fion[ia];
            m2_      = pmass[ia];
            locked2_ = (atmove[ia] != 1);
            found++;
        }
    }

    for (vector<string>::const_iterator it = local_names.begin();
         it != local_names.end(); ++it)
    {
        if (name1_.compare(*it) == 0) locally_owned_ = true;
    }

    if (found == 2) locally_active_ = true;

    if (locked1_) m1_ = 1.e15;
    if (locked2_) m2_ = 1.e15;

    // mu/m_i
    f1_ = m2_ / (m1_ + m2_);
    f2_ = m1_ / (m1_ + m2_);

    if (locally_active_)
    {
        assert(tau1_ != nullptr);
        assert(tau1p_ != nullptr);
        assert(fion1_ != nullptr);
        assert(m1_ > 0.0);
        assert(tau2_ != nullptr);
        assert(tau2p_ != nullptr);
        assert(fion2_ != nullptr);
        assert(tau1p_ != tau2p_);
        assert(m2_ > 0.0);
    }
}

////////////////////////////////////////////////////////////////////////////////
bool DistanceConstraint::enforce(void)
{
    if (!locally_active_) return true;

    // if( locally_owned_ )
    //    (*MPIdata::sout)<<"DistanceConstraint::enforce()"<<endl;

    double r12p[3] = { tau1p_[0] - tau2p_[0], tau1p_[1] - tau2p_[1],
        tau1p_[2] - tau2p_[2] };
    double d2      = r12p[0] * r12p[0] + r12p[1] * r12p[1] + r12p[2] * r12p[2];
    double d       = sqrt(d2);

    if (d != d)
    {
        cerr << "mype=" << mype << ", tau1p_[0]=" << tau1p_[0] << endl;
        cerr << "mype=" << mype << ", tau2p_[0]=" << tau2p_[0] << endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (locally_owned_)
        (*MPIdata::sout) << setprecision(8) << "DistanceConstraint, d=" << d
                         << ", constraint = " << distance_ << endl;

    if (fabs(d - distance_) < tol_) return true;

    const double sigma = distance_ * distance_ - d2;

    // make one shake iteration

    double lambda = 0.5 * sigma / d2;

    const double z0 = lambda * (tau1_[0] - tau2_[0]);
    const double z1 = lambda * (tau1_[1] - tau2_[1]);
    const double z2 = lambda * (tau1_[2] - tau2_[2]);

    if (!locked1_)
    {
        tau1p_[0] += f1_ * z0;
        tau1p_[1] += f1_ * z1;
        tau1p_[2] += f1_ * z2;
    }

    if (!locked2_)
    {
        tau2p_[0] -= f2_ * z0;
        tau2p_[1] -= f2_ * z1;
        tau2p_[2] -= f2_ * z2;
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////
bool DistanceConstraint::project_out_forces(void)
{
    if (!locally_active_) return true;

    double dx = tau1_[0] - tau2_[0];
    double dy = tau1_[1] - tau2_[1];
    double dz = tau1_[2] - tau2_[2];
    double d2 = dx * dx + dy * dy + dz * dz;
    d2 *= 2.;
    assert(d2 > 0.);
    double d2inv = 1. / d2;
    double proj1 = (dx * fion1_[0] + dy * fion1_[1] + dz * fion1_[2]);
    double proj2 = (dx * fion2_[0] + dy * fion2_[1] + dz * fion2_[2]);
    double proj  = (proj1 - proj2) * d2inv;
#if 0
    if( locally_owned_ ){
        (*MPIdata::sout)<< "distance constraint between ";
        (*MPIdata::sout)<< setw(4) << name1_ << " ";
        (*MPIdata::sout)<< setw(4) << name2_ << " ";
        (*MPIdata::sout)<<", force on each atom="<<proj*sqrt(d2inv)<<endl;  
    }
#endif
    if (fabs(proj) < tol_) return true;

    fion1_[0] -= proj * dx;
    fion1_[1] -= proj * dy;
    fion1_[2] -= proj * dz;

    fion2_[0] += proj * dx;
    fion2_[1] += proj * dy;
    fion2_[2] += proj * dz;

    return false;
}

////////////////////////////////////////////////////////////////////////////////
double DistanceConstraint::projected_force(void) const
{
    if (!locally_active_) return 0.;

    // if( locally_owned_ )
    //    (*MPIdata::sout)<<"DistanceConstraint::projected_force()"<<endl;
    const double dx    = tau1_[0] - tau2_[0];
    const double dy    = tau1_[1] - tau2_[1];
    const double dz    = tau1_[2] - tau2_[2];
    const double d2    = dx * dx + dy * dy + dz * dz;
    const double proj1 = (dx * fion1_[0] + dy * fion1_[1] + dz * fion1_[2]);
    const double proj2 = (dx * fion2_[0] + dy * fion2_[1] + dz * fion2_[2]);

    return (proj2 - proj1) / (2. * sqrt(d2));
}

////////////////////////////////////////////////////////////////////////////////
ostream& DistanceConstraint::printConstraint(ostream& os) const
{
    if (!locally_owned_) return os;

    os.setf(ios::left, ios::adjustfield);
    os << "constraint ";
    os << setw(4) << name1_ << " ";
    os << setw(4) << name2_ << " ";
    os.setf(ios::fixed, ios::floatfield);
    os.setf(ios::right, ios::adjustfield);
    os << setw(10) << setprecision(6) << distance_;
    return os;
}

////////////////////////////////////////////////////////////////////////////////
void DistanceConstraint::print(ostream& os) const
{
    if (!locally_owned_) return;

    printConstraint(os);
    os << endl;
}

////////////////////////////////////////////////////////////////////////////////
void DistanceConstraint::printForce(ostream& os) const
{
    if (!locally_owned_) return;

    printConstraint(os);
    os << ", force = " << setw(10) << setprecision(6) << projected_force()
       << endl;
}
