// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#include "ConstraintSet.h"
#include "Control.h"
#include "DistanceConstraint.h"
#include "Ions.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "MultiDistanceConstraint.h"
#include "tools.h"

#include <iomanip>
#include <iostream>
#include <string.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::addConstraint(Ions& ions, const vector<string>& argv)
{
    int argc = (int)argv.size();

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(&argc, 1);

    if (onpe0)
    {
        (*MPIdata::sout) << "ConstraintSet::addConstraint" << endl;
        for (int i = 0; i < argc; i++)
            (*MPIdata::sout) << argv[i] << "\t";
        (*MPIdata::sout) << endl;
    }
    enum constraint_type
    {
        unknown,
        distance_type,
        multidistance_type,
        lock_type
    } type
        = unknown;
    double distance;
    const double tolerance = 1.e-6;

    vector<int> len(argc);
    const int size_buf = 256;
    char char_buf[size_buf];
    int shift = 0;
    if (onpe0)
        for (int i = 0; i < argc; i++)
        {
            len[i] = argv[i].size();
            memcpy(
                char_buf + shift, argv[i].c_str(), sizeof(char) * (len[i] + 1));
            shift += len[i] + 1;
        }
    mmpi.bcast(&len[0], argc);
    mmpi.bcast(char_buf, size_buf);
    shift = 0;
    vector<char*> largv(argc);
    for (int i = 0; i < argc; i++)
    {
        largv[i] = new char[len[i]];
        memcpy(largv[i], char_buf + shift, sizeof(char) * (len[i] + 1));
        shift += len[i] + 1;
    }

    if (!strcmp(largv[1], "distance"))
    {
        type = distance_type;
    }
    else if (!strcmp(largv[1], "multidistance"))
    {
        type = multidistance_type;
    }
    else if (!strcmp(largv[1], "lock"))
    {
        type = lock_type;
    }
    else
    {
        (*MPIdata::sout) << " Incorrect constraint type: " << largv[1] << endl;
        return false;
    }

    if (type == distance_type)
    {
        // constraint distance A B dist

        if (argc != 5)
        {
            (*MPIdata::sout)
                << " Incorrect number of arguments for distance constraint"
                << endl;
            return false;
        }
        string name1 = largv[2];
        string name2 = largv[3];
        if (name1 == name2)
        {
            if (onpe0)
                (*MPIdata::sout) << "ERROR: ConstraintSet: cannot define "
                                    "distance constraint between "
                                 << name1 << " and " << name2 << endl;
            return false;
        }

        Ion* ion1 = ions.findIon(name1);
        Ion* ion2 = ions.findIon(name2);

        short found_ion_local = 1;
        if (ion1 == nullptr) found_ion_local = 0;
        if (ion2 == nullptr) found_ion_local = 0;

        short found_ion = found_ion_local;
        mmpi.allreduce(&found_ion, 1, MPI_MAX);

        if (found_ion == 0)
        {
            if (onpe0)
                (*MPIdata::sout) << "ERROR: ConstraintSet: cannot define "
                                    "distance constraint between "
                                 << name1 << " and " << name2 << endl;
            return false;
        }

        distance = atof(largv[4]);
        if (distance <= 0.0)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << " ConstraintSet: distance must be positive " << endl
                    << " ConstraintSet: could not define constraint " << endl;
            return false;
        }

        // check if equivalent constraint is already defined
        bool found = false;
        for (int i = 0; i < (int)vconstraints_.size(); i++)
        {
            Constraint* pc = vconstraints_[i];
            assert(pc != nullptr);
            // check if a constraint (name1,name2) or (name2,name1) is defined
            if (pc->type() == "distance"
                && ((pc->names(0) == name1 && pc->names(1) == name2)
                    || (pc->names(1) == name1 && pc->names(0) == name2)))
                found = true;
        }

        if (found)
        {
            if (onpe0)
                (*MPIdata::sout)
                    << "ERROR: ConstraintSet:addConstraint: a distance "
                       "constraint "
                    << name1 << " " << name2 << " is already defined" << endl
                    << " ConstraintSet: could not add constraint " << endl;
            return false;
        }

        DistanceConstraint* c
            = new DistanceConstraint(name1, name2, distance, tolerance);

        vconstraints_.push_back(c);
    }
    else if (type == multidistance_type)
    {
        // constraint multidistance alpha1 A1 B1 alpha2 A2 B2 ... d
        int nc = 0;
        if (argc < 6)
        {
            if (onpe0)
                (*MPIdata::sout) << "ERROR:Incorrect number of arguments for "
                                    "multidistance constraint"
                                 << endl;
            return false;
        }
        if (argc % 3 == 0)
        {
            distance = atof(largv[argc - 1]);
            nc       = argc - 1;
        }
        else
        {
            if (onpe0)
                (*MPIdata::sout) << "ERROR: Incorrect number of arguments for "
                                    "multidistance constraint"
                                 << endl;
            return false;
        }

        vector<double> m_alpha;
        vector<string> m_name1, m_name2;
        // bool found_all_ions_in_constraint=true;
        for (int ic = 2; ic < nc; ic = ic + 3)
        {
            double alpha12 = atof(largv[ic]);
            string name1   = largv[ic + 1];
            string name2   = largv[ic + 2];
            if (name1 == name2)
            {
                if (onpe0)
                    (*MPIdata::sout) << "ERROR: ConstraintSet: cannot define "
                                        "distance constraint between "
                                     << name1 << " and " << name2 << endl;
                return false;
            }
            //(*MPIdata::sout)<<"name1="<<name1<<", name2="<<name2<<endl;

            Ion* ion1 = ions.findIon(name1);
            Ion* ion2 = ions.findIon(name2);

            short found_ion_local = 1;
            if (ion1 == nullptr) found_ion_local = 0;
            if (ion2 == nullptr) found_ion_local = 0;

            // check if the two names have at least been found by one MPI task
            short found_ion = found_ion_local;
            mmpi.allreduce(&found_ion, 1, MPI_MAX);

            if (found_ion == 0)
            {
                if (onpe0)
                    (*MPIdata::sout)
                        << "ERROR: ConstraintSet: could not find atom " << name1
                        << " or atom " << name2 << endl;
                return false;
            }

            // put the atom names in alphabetical order
            if (name1 < name2)
            {
                m_name1.push_back(name1);
                m_name2.push_back(name2);
            }
            else
            {
                m_name1.push_back(name2);
                m_name2.push_back(name1);
            }
            m_alpha.push_back(alpha12);

            // if( found_ion_local==0 )
            //{
            //    found_all_ions_in_constraint=false;
            //}
        }

        // register constraint on every MPI task, since atoms may move
        // we are thus assuming the number of constraints is small
        // and not growing with system size
        MultiDistanceConstraint* c = new MultiDistanceConstraint(
            m_alpha, m_name1, m_name2, distance, tolerance);

        vconstraints_.push_back(c);
    }
    else if (type == lock_type)
    {
        string name1 = largv[2];
        ions.lockAtom(name1);
    }
    else
    {
        (*MPIdata::sout) << "ConstraintSet::addConstraint: internal error"
                         << endl;
        return false;
    }

    for (int i = 0; i < argc; i++)
    {
        delete[] largv[i];
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::printConstraints(ostream& os) const
{
    for (vector<Constraint*>::const_iterator it = vconstraints_.begin();
         it != vconstraints_.end(); ++it)
    {
        (*it)->print(os);
    }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::printConstraintsForces(ostream& os)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        os << " ConstraintSet: constraints and forces" << endl;
    for (vector<Constraint*>::iterator it = vconstraints_.begin();
         it != vconstraints_.end(); ++it)
    {
        (*it)->printForce(os);
    }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::enforceConstraints(const int maxiter)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << " ConstraintSet::enforceConstraints()" << endl;

    int iter  = 0;
    bool done = false;
    while (!done && (iter < maxiter))
    {
        done = true;
        for (int i = 0; i < (int)vconstraints_.size(); i++)
        {
            bool b = vconstraints_[i]->enforce();
            done &= b;
        }
        iter++;
    }

    // need a reduce to make sure all constraints are enforced
    short local_done  = done ? 1 : 0;
    short global_done = 0;
    mmpi.allreduce(&local_done, &global_done, 1, MPI_MIN);
    done = (global_done == 1);

    if (!done && onpe0)
    {
        (*MPIdata::sout)
            << " WARNING, ConstraintSet: could not enforce constraints in "
            << maxiter << " iterations" << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::projectOutForces(const int maxiter)
{
    int iter  = 0;
    bool done = false;
    while (!done && (iter < maxiter))
    {
        done = true;
        for (int i = 0; i < (int)vconstraints_.size(); i++)
        {
            bool b = vconstraints_[i]->project_out_forces();
            done &= b;
        }
        iter++;
    }

    // need a reduce to make sure all constraints are enforced
    short local_done  = done ? 1 : 0;
    short global_done = 0;
    MGmol_MPI& mmpi   = *(MGmol_MPI::instance());
    mmpi.allreduce(&local_done, &global_done, 1, MPI_MIN);
    done = (global_done == 1);

    if (!done && onpe0)
    {
        (*MPIdata::sout)
            << "WARNING, ConstraintSet: could not project out forces in "
            << maxiter << " iterations" << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::setup(const vector<string>& names,
    const vector<double>& pmass, const vector<short>& imove,
    vector<double*>& tau, vector<double*>& taup, vector<double*>& fion,
    const vector<string>& local_names)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        cout << "ConstraintSet::setup()..." << vconstraints_.size()
             << " constraints" << endl;
    for (vector<Constraint*>::iterator it = vconstraints_.begin();
         it != vconstraints_.end(); ++it)
    {
        (*it)->setup(names, pmass, imove, tau, taup, fion, local_names);
    }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::setup(Ions& ions)
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0) printWithTimeStamp("ConstraintSet::setup()...", cout);

    const vector<double>& pmass(ions.getMassesInteractingIons());

    const vector<string>& names(ions.getNamesInteractingIons());
    const vector<string>& local_names(ions.getLocalNames());

    vector<double*>& tau0(ions.getPositionsInteractingIons());
    vector<double*>& taup(ions.getTaupInteractingIons());
    // cout<<"ConstraintSet::setup(), Size of tau="<<tau0.size()<<endl;

    vector<double*>& fion(ions.getForcesInteractingIons());

    const vector<short>& atmove(ions.getMovesInteractingIons());

    setup(names, pmass, atmove, tau0, taup, fion, local_names);
}

int ConstraintSet::readConstraints(const string& filename)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    ifstream* tfile = nullptr;
    if (mmpi.instancePE0())
    {
        tfile = new ifstream(filename.data(), ios::in);
        if (!tfile->is_open())
        {
            cerr << "MGmol::readConstraints --- Unable to open file "
                 << filename.data() << endl;
            return -1;
        }
        else
        {
            cout << "Open " << filename.data() << endl;
        }
    }

    int nconstraints = readConstraints(tfile);

    if (mmpi.instancePE0())
    {
        delete tfile;
    }

    return nconstraints;
}

int ConstraintSet::readConstraints(ifstream* tfile)
{
    Control& ct = *(Control::instance());

    if (ct.verbose > 0)
        printWithTimeStamp("ConstraintSet::readConstraints()...", cout);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    // read constraints into msav_ from PE 0
    msav_.clear();
    if (mmpi.instancePE0() && tfile != nullptr)
    {
        // if(onpe0)cout<<"readInput: Read constraints..."<<endl;
        read_comments(*tfile);
        while (!tfile->eof())
        {
            string query;
            if (!getline(*tfile, query))
            {
                break;
            }
            stringstream ss(query);
            string command_string;
            ss >> command_string;
            if (command_string.compare("constraint") == 0)
            {
                if (onpe0) cout << "Reading 1 constraint..." << endl;
                string type_string;
                ss >> type_string;
                vector<string> sav;
                if (type_string.compare("distance") == 0)
                {
                    string name1;
                    string name2;
                    ss >> name1;
                    ss >> name2;
                    string distance;
                    ss >> distance;
                    read_comments(*tfile);
                    sav.push_back(command_string);
                    sav.push_back(type_string);
                    sav.push_back(name1);
                    sav.push_back(name2);
                    sav.push_back(distance);
                }
                else if (type_string.compare("multidistance") == 0)
                {
                    string alpha;
                    string name1;
                    string name2;
                    string distance;
                    vector<string> valpha;
                    vector<string> vname1;
                    vector<string> vname2;
                    int n = 0;
                    while (ss >> alpha)
                    {
                        if (ss >> name1)
                        {
                            ss >> name2;
                            valpha.push_back(alpha);
                            vname1.push_back(name1);
                            vname2.push_back(name2);
                            n++;
                        }
                        else
                        {
                            distance = alpha;
                            break;
                        }
                    }
                    read_comments(*tfile);

                    sav.push_back(command_string);
                    sav.push_back(type_string);
                    for (int i = 0; i < n; i++)
                    {
                        sav.push_back(valpha[i]);
                        sav.push_back(vname1[i]);
                        sav.push_back(vname2[i]);
                    }
                    sav.push_back(distance.c_str());
                }
                else if (type_string.compare("lock") == 0)
                {
                    string name;
                    ss >> name;
                    // ions->lockAtom(name);
                    read_comments(*tfile);
                    sav.push_back(command_string);
                    sav.push_back(type_string);
                    sav.push_back(name);
                }
                else
                {
                    cerr << "ERROR reading constraint type=" << type_string
                         << endl;
                    return -1;
                }
                msav_.push_back(sav);
                //}else{
                //    break;
            }
        }
    }

    // Bcast constraints
    int nconstraints = (int)msav_.size();
    mmpi.bcast(&nconstraints);

    if (onpe0 && ct.verbose > 0)
    {
        cout << "Read " << nconstraints << " nconstraints" << endl;
    }

    return nconstraints;
}

int ConstraintSet::addConstraints(Ions& ions)
{
    clear();

    MGmol_MPI& mmpi  = *(MGmol_MPI::instance());
    int nconstraints = (int)msav_.size();
    mmpi.bcast(&nconstraints);

    for (int i = 0; i < nconstraints; i++)
    {
        vector<string> sav;
        if (mmpi.instancePE0())
        {
            sav = msav_[i];
        }
        if (!addConstraint(ions, sav)) return -1;
    }

    mmpi.barrier();

    return nconstraints;
}

void ConstraintSet::clear()
{
    for (vector<Constraint*>::iterator it = vconstraints_.begin();
         it != vconstraints_.end(); ++it)
    {
        delete *it;
    }

    vconstraints_.clear();
}
