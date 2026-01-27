// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$

#include "MultiDistanceConstraint.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <unistd.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void MultiDistanceConstraint::setup(const vector<string>& names_known_locally,
    const vector<double>& pmass, const vector<short>& atmove,
    vector<double*>& tau, vector<double*>& taup, vector<double*>& fion,
    const vector<string>& local_names)
{
    const int nc = m_alpha_.size();

    m_tau1_.resize(nc);
    m_tau1p_.resize(nc);
    m_fion1_.resize(nc);
    m_tau2_.resize(nc);
    m_tau2p_.resize(nc);
    m_fion2_.resize(nc);

    u_name_.clear();
    u_mass_.clear();
    u_taup_.clear();
    u_fion_.clear();
    u_locked_.clear();

    short nfound = 0;

    // loop over all the distances making up this constraint
    for (int i = 0; i < nc; i++)
    {
        m_tau1_[i]  = nullptr;
        m_tau1p_[i] = nullptr;
        m_fion1_[i] = nullptr;

        m_tau2_[i]  = nullptr;
        m_tau2p_[i] = nullptr;
        m_fion2_[i] = nullptr;

        for (int ia = 0; ia < (int)names_known_locally.size(); ia++)
        {
            if (m_name1_[i].compare(names_known_locally[ia]) == 0)
            {
                bool found = false;
                // is m_name1_[i] already involved in a constraint
                for (int j = 0; j < (int)u_name_.size(); j++)
                    if (u_name_[j].compare(m_name1_[i]) == 0) found = true;
                if (!found)
                {
                    u_name_.push_back(names_known_locally[ia]);
                    bool locked = (atmove[ia] != 1);
                    double mass = locked ? 1.e15 : pmass[ia];
                    u_mass_.push_back(mass);
                    u_taup_.push_back(taup[ia]);
                    u_fion_.push_back(fion[ia]);
                    u_locked_.push_back(atmove[ia] != 1);
                }
                m_tau1_[i]  = tau[ia];
                m_tau1p_[i] = taup[ia];
                m_fion1_[i] = fion[ia];
                nfound++;
            }
            else if (m_name2_[i].compare(names_known_locally[ia]) == 0)
            {
                bool found = false;
                for (int j = 0; j < (int)u_name_.size(); j++)
                    if (u_name_[j].compare(m_name2_[i]) == 0) found = true;
                if (!found)
                {
                    u_name_.push_back(names_known_locally[ia]);
                    bool locked = (atmove[ia] != 1);
                    double mass = locked ? 1.e15 : pmass[ia];
                    u_mass_.push_back(mass);
                    u_taup_.push_back(taup[ia]);
                    u_fion_.push_back(fion[ia]);
                    u_locked_.push_back(atmove[ia] != 1);
                }
                m_tau2_[i]  = tau[ia];
                m_tau2p_[i] = taup[ia];
                m_fion2_[i] = fion[ia];
                nfound++;
            }
        }
    }

    // attributes ownership to task where m_name1_[0] is located
    if (!m_name1_.empty())
    {
        for (vector<string>::const_iterator it = local_names.begin();
             it != local_names.end(); ++it)
        {
            if (m_name1_[0].compare(*it) == 0) locally_owned_ = true;
        }
    }
    else
    {
        locally_owned_ = false;
    }

    short lowned    = (short)locally_owned_;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    short check     = 0;
    mmpi.allreduce(&lowned, &check, 1, MPI_MAX);
    if (check != 1)
    {
        int natoms  = local_names.size();
        int tnatoms = 0;
        mmpi.allreduce(&natoms, &tnatoms, 1, MPI_SUM);

        cerr << "MGmol ERROR: no task owns constraint" << endl;
        cerr << "Cannot find " << m_name1_[0]
             << " in any list of local atoms..." << endl;
        cerr << "Total number of local atoms: " << tnatoms << endl;
        sleep(5);
        mmpi.abort();
    }

    locally_active_ = (nfound == 2 * nc);

    if (locally_active_)
        for (int i = 0; i < nc; i++)
        {
            assert(m_tau1_[i] != nullptr);
            assert(m_tau2_[i] != nullptr);
            assert(m_tau1p_[i] != nullptr);
            assert(m_tau2p_[i] != nullptr);
            assert(m_fion1_[i] != nullptr);
            assert(m_fion2_[i] != nullptr);
        }

    short lactive = (short)locally_active_;
    mmpi.allreduce(&lactive, &check, 1, MPI_MAX);
    if (check != 1)
    {
        (*MPIdata::sout) << "ERROR: no task enforces constraint" << endl;
        sleep(5);
        mmpi.abort();
    }
}

////////////////////////////////////////////////////////////////////////////////

bool MultiDistanceConstraint::enforce(void)
{
    if (!locally_active_) return true;

    const int nc = m_alpha_.size();

    double sigma = -distance_;
    for (int k = 0; k < nc; k++)
    {
        double dx = m_tau1p_[k][0] - m_tau2p_[k][0];
        double dy = m_tau1p_[k][1] - m_tau2p_[k][1];
        double dz = m_tau1p_[k][2] - m_tau2p_[k][2];
        double d  = sqrt(dx * dx + dy * dy + dz * dz);
        // if( onpe0 )
        //    (*MPIdata::sout)<<"alpha="<<m_alpha_[k]<<", d="<<d<<"\t";
        sigma += m_alpha_[k] * d;
    }
    if (locally_owned_)
    {
        (*MPIdata::sout) << setprecision(8)
                         << "MultiDistanceConstraint, sigma=" << sigma
                         << ", constraint = " << distance_ << endl;
    }

    if (fabs(sigma) < tol_) return true;

    vector<vector<double>> a_;
    const int na = (int)u_name_.size();
    a_.resize(na);
    for (int i = 0; i < na; i++)
        a_[i].resize(3);

    // make one shake iteration
    for (int p = 0; p < na; p++)
    {
        for (int e = 0; e < 3; e++)
        {
            a_[p][e] = 0.;
            for (int k = 0; k < nc; k++)
            {
                double dx = m_tau1p_[k][0] - m_tau2p_[k][0];
                double dy = m_tau1p_[k][1] - m_tau2p_[k][1];
                double dz = m_tau1p_[k][2] - m_tau2p_[k][2];
                double d  = sqrt(dx * dx + dy * dy + dz * dz);

                int d1 = 0;
                int d2 = 0;
                if (u_name_[p].compare(m_name1_[k]) == 0) d1 = 1;
                if (u_name_[p].compare(m_name2_[k]) == 0) d2 = 1;

                a_[p][e] += m_alpha_[k] * (d1 - d2)
                            * (m_tau1_[k][e] - m_tau2_[k][e]) / d;
            }
        }
    }

    double c = 0.;
    for (int p = 0; p < na; p++)
    {
        for (int e = 0; e < 3; e++)
        {
            c += a_[p][e] * a_[p][e] / u_mass_[p];
        }
    }

    const double lambda = -sigma / c;
    for (int p = 0; p < na; p++)
    {
        if (!u_locked_[p])
        {
            for (int e = 0; e < 3; e++)
            {
                u_taup_[p][e] += lambda * a_[p][e] / u_mass_[p];
            }
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////
bool MultiDistanceConstraint::project_out_forces(void)
{
    if (!locally_active_) return true;

    const int nc = m_alpha_.size();
    const int na = u_name_.size();
    double proj  = 0.;

    vector<vector<double>> a_;
    a_.resize(na);
    for (int i = 0; i < na; i++)
        a_[i].resize(3);

    // loop over (unique) atoms
    for (int p = 0; p < na; p++)
    {
        a_[p][0] = 0.;
        a_[p][1] = 0.;
        a_[p][2] = 0.;

        // loop over constraint distances
        for (int k = 0; k < nc; k++)
        {
            short d1 = 0;
            short d2 = 0;
            if (u_name_[p] == m_name1_[k]) d1 = 1;
            if (u_name_[p] == m_name2_[k]) d2 = 1;
            if (d1 || d2)
            {
                double dx  = m_tau1_[k][0] - m_tau2_[k][0];
                double dy  = m_tau1_[k][1] - m_tau2_[k][1];
                double dz  = m_tau1_[k][2] - m_tau2_[k][2];
                double dr2 = dx * dx + dy * dy + dz * dz;
                assert(dr2 > 0.);
                double dd = (double)(d1 - d2);
                dd *= m_alpha_[k];
                dd /= sqrt(dr2);
                // projection
                proj += dd
                        * (dx * u_fion_[p][0] + dy * u_fion_[p][1]
                            + dz * u_fion_[p][2]);
                // gradient
                a_[p][0] += dd * dx;
                a_[p][1] += dd * dy;
                a_[p][2] += dd * dz;
            }
        }
    }

    // get square norm of gradient
    // (may be 0. if atoms are not local and forces are not computed for these
    // atoms)
    double md2 = 0.;
    for (int p = 0; p < na; p++)
    {
        md2 += a_[p][0] * a_[p][0] + a_[p][1] * a_[p][1] + a_[p][2] * a_[p][2];
    }

    if (fabs(sqrt(md2)) < tol_) return true;
    //(*MPIdata::sout)<<"proj/sqrt(md2)="<<proj/sqrt(md2)<<endl;
    if (fabs(proj / sqrt(md2)) < tol_) return true;

    proj /= md2;

    // substract projection
    for (int p = 0; p < na; p++)
    {
        u_fion_[p][0] -= proj * a_[p][0];
        u_fion_[p][1] -= proj * a_[p][1];
        u_fion_[p][2] -= proj * a_[p][2];
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

double MultiDistanceConstraint::projected_force(void) const
{
    if (!locally_active_) return 0.;

    int nc    = m_alpha_.size();
    int na    = u_name_.size();
    double de = 0.;

    vector<vector<double>> a_;
    a_.resize(na);
    for (int i = 0; i < na; i++)
        a_[i].resize(3);

    // loop over atoms
    for (int p = 0; p < na; p++)
    {
        a_[p][0] = 0.;
        a_[p][1] = 0.;
        a_[p][2] = 0.;

        // loop over constraint distances
        for (int k = 0; k < nc; k++)
        {
            int d1 = 0;
            int d2 = 0;
            if (u_name_[p] == m_name1_[k]) d1 = 1;
            if (u_name_[p] == m_name2_[k]) d2 = 1;
            if (d1 || d2)
            {
                double dd  = (double)(d1 - d2);
                double dx  = m_tau1_[k][0] - m_tau2_[k][0];
                double dy  = m_tau1_[k][1] - m_tau2_[k][1];
                double dz  = m_tau1_[k][2] - m_tau2_[k][2];
                double dr2 = dx * dx + dy * dy + dz * dz;
                assert(dr2 > 0.);
                dd *= m_alpha_[k];
                dd /= sqrt(dr2);
                // projection
                de -= dd
                      * (dx * u_fion_[p][0] + dy * u_fion_[p][1]
                          + dz * u_fion_[p][2]);
                // gradient
                a_[p][0] += dd * dx;
                a_[p][1] += dd * dy;
                a_[p][2] += dd * dz;
            }
        }
    }

    // get square norm of gradient
    double md2 = 0.;
    for (int p = 0; p < na; p++)
    {
        md2 += (a_[p][0] * a_[p][0] + a_[p][1] * a_[p][1]
                + a_[p][2] * a_[p][2]);
    }
    assert(md2 > 0.);
    de /= md2;

    return de;
}

////////////////////////////////////////////////////////////////////////////////

ostream& MultiDistanceConstraint::printConstraint(ostream& os) const
{
    if (!locally_owned_) return os;

    os.setf(ios::left, ios::adjustfield);
    os << "constraint ";
    os << type() << " ";
    for (int ic = 0; ic < (int)m_alpha_.size(); ic++)
        os << setw(10) << setprecision(6) << m_alpha_[ic] << " " << m_name1_[ic]
           << " " << m_name2_[ic] << " ";
    os.setf(ios::fixed, ios::floatfield);
    os.setf(ios::right, ios::adjustfield);
    os << setw(10) << setprecision(6) << distance_;
    return os;
}

////////////////////////////////////////////////////////////////////////////////
void MultiDistanceConstraint::print(ostream& os) const
{
    if (!locally_owned_) return;

    printConstraint(os);
    os << endl;
}

////////////////////////////////////////////////////////////////////////////////
void MultiDistanceConstraint::printForce(ostream& os) const
{
    if (!locally_owned_) return;

    printConstraint(os);
    os << ", force = " << setw(10) << setprecision(6) << projected_force()
       << endl;
}
