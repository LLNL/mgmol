// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cmath>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
using namespace std;

//#include <H5LT.h>
#include "Control.h"
#include "HDFrestart.h"
#include "Ions.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "Species.h"
#include "Timer.h"
#include "hdf_tools.h"
#include "mgmol_mpi_tools.h"
#include "tools.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#define max(a, b) (((a) < (b)) ? (b) : (a))

Timer ions_setupInteractingIons_tm("ions_setupInteractingIons");
Timer ions_setup_tm("ions::setup");

const double ang2bohr = 1.8897269;
// const double rmax = 8.0;

map<string, short> Ions::map_species_;
int Ions::num_ions_          = -1;
short Ions::max_num_proj_    = -1;
double Ions::max_Vl_radius_  = -1.;
double Ions::max_Vnl_radius_ = -1.;

Ions::Ions(const double lat[3], const vector<Species>& sp) : species_(sp)
{
    for (short i = 0; i < 3; i++)
    {
        assert(lat[i] > 0.);
        if (lat[i] > 10000.)
        {
            (*MPIdata::serr) << "Ions constructor: lattice[" << i
                             << "]=" << lat[i] << "!!!" << endl;
            exit(2);
        }
        lattice_[i] = lat[i];
    }
    setup_            = false;
    has_locked_atoms_ = false;

    map_species_["H"]  = 1;
    map_species_["Li"] = 3;
    map_species_["Be"] = 4;
    map_species_["B"]  = 5;
    map_species_["C"]  = 6;
    map_species_["N"]  = 7;
    map_species_["O"]  = 8;
    map_species_["F"]  = 9;
    map_species_["Na"] = 11;
    map_species_["Mg"] = 12;
    map_species_["Al"] = 13;
    map_species_["Si"] = 14;
    map_species_["P"]  = 15;
    map_species_["S"]  = 16;
    map_species_["Cl"] = 17;
    map_species_["K"]  = 19;
    map_species_["Ca"] = 20;
    map_species_["Cr"] = 24;
    map_species_["Mn"] = 25;
    map_species_["Fe"] = 26;
    map_species_["Co"] = 27;
    map_species_["Ni"] = 28;
    map_species_["Cu"] = 29;
    map_species_["Zn"] = 30;
    map_species_["Ga"] = 31;
    map_species_["Ge"] = 32;
    map_species_["Au"] = 79;

    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    for (short i = 0; i < 3; i++)
        div_lattice_[i] = lattice_[i] / (double)(myPEenv.n_mpi_task(i));

    double offset[3] = { (double)myPEenv.my_mpi(0) * div_lattice_[0],
        (double)myPEenv.my_mpi(1) * div_lattice_[1],
        (double)myPEenv.my_mpi(2) * div_lattice_[2] };

    for (short i = 0; i < 3; i++)
        if (myPEenv.n_mpi_task(i) == (myPEenv.my_mpi(i) + 1))
            div_lattice_[i] = lattice_[i] - offset[i];

    // save cartesian communicator info
    cart_comm_     = myPEenv.cart_comm();
    const int disp = -1;
    for (int dir = 0; dir < 3; dir++)
        MPI_Cart_shift(cart_comm_, dir, disp, &source_[dir], &dest_[dir]);
}

Ions::Ions(const Ions& ions, const double shift[3]) : species_(ions.species_)
{
    vector<Ion*>::const_iterator ion = ions.list_ions_.begin();
    while (ion != ions.list_ions_.end())
    {
        Ion* newion = new Ion(**ion);
        newion->shiftPosition(shift);
        newion->setup();
        list_ions_.push_back(newion);
        ion++;
    }
    for (short i = 0; i < 3; ++i)
        lattice_[i] = ions.lattice_[i];

    vector<Ion*>::iterator iion       = list_ions_.begin();
    vector<Ion*>::const_iterator cion = ions.list_ions_.begin();
    while (iion != list_ions_.end())
    {
        (*iion)->set_here((*cion)->here());
        if ((*cion)->here())
        {
            local_ions_.push_back(*iion);
        }
        iion++;
        cion++;
    }
    local_names_      = ions.local_names_;
    atmove_           = ions.atmove_;
    pmass_            = ions.pmass_;
    taum_             = ions.taum_;
    tau0_             = ions.tau0_;
    taup_             = ions.taup_;
    fion_             = ions.fion_;
    velocity_         = ions.velocity_;
    rand_states_      = ions.rand_states_;
    has_locked_atoms_ = ions.has_locked_atoms_;

    gids_ = ions.gids_;

    setupListOverlappingIons();

    setupInteractingIons();

    lstep_[0] = ions.lstep_[0];
    lstep_[1] = ions.lstep_[1];
    lstep_[2] = ions.lstep_[2];
    rstep_[0] = ions.rstep_[0];
    rstep_[1] = ions.rstep_[1];
    rstep_[2] = ions.rstep_[2];
}

void Ions::computeMaxNumProjs()
{
#ifdef DEBUG
    if (onpe0)
        (*MPIdata::sout) << " list_ions of size " << list_ions_.size()
                         << " initialized" << endl;
#endif

    vector<Ion*>::iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        max_num_proj_ = max(max_num_proj_, (*ion)->nProjectors());

        ion++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&max_num_proj_, 1, MPI_MAX);
}

void Ions::setup()
{
    Control& ct = *(Control::instance());

    if (ct.verbose > 0) printWithTimeStamp("Ions::setup()...", cout);

    ions_setup_tm.start();

    updateListIons();

    //#ifndef NDEBUG
    //    checkUnicityLocalIons();
    //#endif

    setupInteractingIons();

    // initialize ionic stepper data
    initStepperData();

    has_locked_atoms_ = hasLockedAtoms();

    // initialize data for constraints
    setupContraintsData(interacting_ions_);

    computeMaxNumProjs();

    if (ct.verbose > 0)
        printWithTimeStamp("Ions::setup()... individual ions...", cout);

    vector<Ion*>::iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        (*ion)->setup();
        ion++;
    }

    setMapVL();
    setupListOverlappingIons();

    max_Vl_radius_  = computeMaxVlRadius();
    max_Vnl_radius_ = computeMaxNLprojRadius();

    setup_ = true;

    assert(max_Vnl_radius_ >= 0.);

    ions_setup_tm.stop();
}

Ions::~Ions()
{
    vector<Ion*>::iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        delete *ion;
        ion++;
    }
}

void Ions::setupListOverlappingIons()
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("Ions::setupListOverlappingIons()...", cout);

    overlappingNL_ions_.clear();
    overlappingVL_ions_.clear();

    vector<Ion*>::const_iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        if ((*ion)->map_nl())
        {
            overlappingNL_ions_.push_back(*ion);
        }
        ion++;
    }

    ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        if ((*ion)->map_l())
        {
            overlappingVL_ions_.push_back(*ion);
        }
        ion++;
    }
}

void Ions::setupInteractingIons()
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("Ions::setupInteractingIons()...", cout);

    ions_setupInteractingIons_tm.start();
    const double rmax = ct.maxDistanceAtomicInfo();
    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::setupInteractingIons() with radius " << rmax
                         << endl;
    }

    interacting_ions_.clear();

    vector<Ion*>::const_iterator ion1 = list_ions_.begin();
    while (ion1 != list_ions_.end())
    {
        // is ion1 interacting with any local ions?
        vector<Ion*>::const_iterator ion2 = local_ions_.begin();
        while (ion2 != local_ions_.end())
        {
            const double r12 = (*ion1)->minimage(**ion2, lattice_, ct.bc);

            if (r12 < rmax)
            {
                // ion1 is interacting with local ions
                interacting_ions_.push_back(*ion1);
                break;
            }

            ion2++;
        }
        ion1++;
    }

    //(*MPIdata::sout)<<"Number of interacting ions =
    //"<<interacting_ions_.size()<<endl;
    ions_setupInteractingIons_tm.stop();
}

// setup arrays to be used in constraints enforcement
// using references to local_ions and extra "dummy" data
void Ions::setupContraintsData(vector<Ion*>& ions_for_constraints)
{
    Control& ct = *(Control::instance());

    if (ct.verbose > 0)
        printWithTimeStamp("Ions::setupContraintsData()...", cout);
    const int nnloc = ions_for_constraints.size() - local_ions_.size();
    // cout<<"interacting_ions_.size()="<<interacting_ions_.size()<<endl;
    // cout<<"local_ions_.size()="<<local_ions_.size()<<endl;
    assert(nnloc >= 0);

    tau0_dummy_.resize(3 * nnloc);
    taup_dummy_.resize(3 * nnloc);
    fion_dummy_.resize(3 * nnloc);
    names_dummy_.clear();
    pmass_dummy_.clear();
    atmove_dummy_.clear();

    interacting_tau0_.clear();
    interacting_taup_.clear();
    interacting_fion_.clear();
    interacting_names_.clear();
    interacting_atmove_.clear();
    interacting_pmass_.clear();

    int ia                           = 0; // count local ions
    int ib                           = 0; // count non-local ions
    vector<Ion*>::const_iterator ion = ions_for_constraints.begin();
    while (ion != ions_for_constraints.end())
    {
        if (isLocal((*ion)->name()))
        {
            assert(3 * ia < tau0_.size());
            assert(3 * ia < taup_.size());
            assert(3 * ia < fion_.size());
            assert(ia < pmass_.size());
            assert(ia < atmove_.size());

            // set reference to local data
            interacting_tau0_.push_back(&tau0_[3 * ia]);

            interacting_taup_.push_back(&taup_[3 * ia]);

            interacting_fion_.push_back(&fion_[3 * ia]);

            interacting_names_.push_back(local_names_[ia]);
            interacting_pmass_.push_back(pmass_[ia]);
            interacting_atmove_.push_back(atmove_[ia]);

            ++ia;
        }
        else
        {
            assert(3 * ib < tau0_dummy_.size());
            assert(3 * ib < taup_dummy_.size());
            assert(3 * ib < fion_dummy_.size());

            // fill up dummy data
            (*ion)->getPosition(&tau0_dummy_[3 * ib]);
            (*ion)->getPosition(&taup_dummy_[3 * ib]);
            (*ion)->getForce(&fion_dummy_[3 * ib]);
            names_dummy_.push_back((*ion)->name());
            pmass_dummy_.push_back((*ion)->getMass());
            atmove_dummy_.push_back(!(*ion)->locked());

            // set references to dummy data
            interacting_tau0_.push_back(&tau0_dummy_[3 * ib]);

            interacting_taup_.push_back(&taup_dummy_[3 * ib]);

            interacting_fion_.push_back(&fion_dummy_[3 * ib]);

            assert(ib < names_dummy_.size());
            assert(ib < pmass_dummy_.size());
            assert(ib < atmove_dummy_.size());
            interacting_names_.push_back(names_dummy_[ib]);
            interacting_pmass_.push_back(pmass_dummy_[ib]);
            interacting_atmove_.push_back(atmove_dummy_[ib]);

            ++ib;
        }

        ++ion;
    }
}

void Ions::iiforce(const short bc[3])
{
    const int nlions = local_ions_.size();
    vector<double> forces(3 * nlions, 0.);

    vector<Ion*>::const_iterator ion1 = local_ions_.begin();
    int ion1_index                    = 0;
    ;
    while (ion1 != local_ions_.end())
    {
        const double z1 = (*ion1)->getZion();
        assert(z1 >= 0.);

        const double rc1 = (*ion1)->getRC();

        vector<Ion*>::const_iterator ion2 = interacting_ions_.begin();
        while (ion2 != interacting_ions_.end())
        {
            if (*ion1 != *ion2)
            {
                // Minimum image convention for r
                double dr[3];
                const double r12 = (*ion1)->minimage(**ion2, lattice_, bc, dr);

                const double invr = 1. / r12;

                const double z2 = (*ion2)->getZion();
                assert(z2 >= 0.);

                const double rc2 = (*ion2)->getRC();

                const double t        = rc1 * rc1 + rc2 * rc2;
                const double invt     = 1. / t;
                const double sqrtinvt = sqrt(invt);

                const double s1 = z1 * z2 * invr * invr;
                const double s2 = erfc(r12 * sqrtinvt) * invr;
                const double s3
                    = M_2_SQRTPI * exp(-r12 * r12 * invt) * sqrtinvt;
                const double alpha = s1 * (s2 + s3);

                for (short i = 0; i < 3; i++)
                    forces[3 * ion1_index + i] += dr[i] * alpha;
            }

            ion2++;
        }
        ion1_index++;
        ion1++;
    }

    vector<Ion*>::iterator lion = local_ions_.begin();
    int ion_index               = 0;
    while (lion != local_ions_.end())
    {
        (*lion)->add_force(forces[3 * ion_index + 0], forces[3 * ion_index + 1],
            forces[3 * ion_index + 2]);
        ion_index++;
        lion++;
    }
}

double Ions::energySelf() const
{
    double energy = 0.;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        energy += (*ion)->eself();
        ion++;
    }
    // multiply by 0.5, 2/sqrt(pi) and 1/sqrt(2) to get 1./sqrt(2*pi)
    energy *= 0.5 * M_2_SQRTPI * M_SQRT1_2;

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&energy, 1, MPI_SUM);

    assert(energy == energy);
    return energy; // Hartree
}

// right-looking approach avoids duplicate computation of energy on different
// processors. This works since ion lists are in ascending order of unique ion
// index
double Ions::energyDiff(const short bc[3]) const
{
    double energy = 0.;

    assert(lattice_[0] > 0.);
    assert(lattice_[1] > 0.);
    assert(lattice_[2] > 0.);

    vector<Ion*>::const_iterator ion1 = interacting_ions_.begin();
    while (ion1 != interacting_ions_.end())
    {
        if ((*ion1)->here())
        {
            vector<Ion*>::const_iterator ion2 = interacting_ions_.begin();
            ;
            while (ion2 != ion1)
            {
                energy += (*ion1)->energyDiff(**ion2, lattice_, bc);

                ion2++;
            }
            // increment ion2 past ion1
            ion2++;
            while (ion2 != interacting_ions_.end())
            {
                energy += (*ion1)->energyDiff(**ion2, lattice_, bc);

                ion2++;
            }
        }

        ion1++;
    }
    assert(energy == energy);

#ifdef USE_MPI
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double tenergy           = 0.;
    MPI_Allreduce(&energy, &tenergy, 1, MPI_DOUBLE, MPI_SUM, myPEenv.comm());

    // take half the result to account for double counting loop over
    // interacting_ions
    energy = 0.5 * tenergy;
#endif

    return energy; // Hartree
}

// Associate each ion with the PE where they are centered
// by setting a flag "here"
// void Ions::associate2PE()
//{
//    Control& ct(*(Control::instance()));
//    Mesh* mymesh             = Mesh::instance();
//    const pb::PEenv& myPEenv = mymesh->peenv();
//    const pb::Grid& mygrid   = mymesh->grid();
//
//    double div_lattice[3];
//    for (short i = 0; i < 3; i++)
//        div_lattice[i] = lattice_[i] / (double)(myPEenv.n_mpi_task(i));
//
//    double offset[3];
//    for (short i = 0; i < 3; i++)
//        offset[i] = (double)myPEenv.my_mpi(i) * div_lattice[i];
//
//#ifdef DEBUG
//    int isum2 = 0;
//    (*MPIdata::sout) << " offset=(" << offset[0] << "," << offset[1] << ","
//                     << offset[2] << ")" << endl;
//#endif
//
//    for (short i = 0; i < 3; i++)
//        if (myPEenv.n_mpi_task(i) == (myPEenv.my_mpi(i) + 1))
//            div_lattice[i] = lattice_[i] - offset[i];
//
//    const double origin[3]
//        = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };
//    const double end[3] = { origin[0] + lattice_[0], origin[1] + lattice_[1],
//        origin[2] + lattice_[2] };
//
//    int isum = 0;
//
//    // Loop over ions
//    local_ions_.clear();
//    vector<Ion*>::iterator ion = list_ions_.begin();
//    while (ion != list_ions_.end())
//    {
//        double t[3];
//        for (short i = 0; i < 3; i++)
//        {
//            t[i] = (*ion)->position(i);
//            while (t[i] >= end[i])
//                t[i] -= lattice_[i];
//            while (t[i] < origin[i])
//                t[i] += lattice_[i];
//            t[i] -= origin[i];
//            t[i] -= offset[i];
//        }
//
//#if DEBUG
//        (*MPIdata::sout) << " t=(" << t[0] << "," << t[1] << "," << t[2] <<
//        ")"
//                         << endl;
//#endif
//
//        if ((t[0] >= 0. && t[0] < (div_lattice[0]))
//            && (t[1] >= 0. && t[1] < (div_lattice[1]))
//            && (t[2] >= 0. && t[2] < (div_lattice[2])))
//        {
//
//            (*ion)->set_here(true);
//
//            local_ions_.push_back(*ion);
//
//            isum++;
//        }
//        else
//        {
//            (*ion)->set_here(false);
//        }
//#ifndef NDEBUG
//        if ((*ion)->here())
//        {
//            (*MPIdata::sout) << " Ion " << (*ion)->name() << " centered on PE
//            "
//                             << myPEenv.mytask() << endl;
//        }
//#endif
//        ion++;
//    }
//    // Test all atoms associated to each processors sum up to
//    // total number of atoms
//
//    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
//    mmpi.allreduce(&isum, 1, MPI_SUM);
//    if (onpe0)
//    {
//        (*MPIdata::sout) << " num_ions=" << num_ions_ << endl;
//        (*MPIdata::sout) << " isum=" << isum << endl;
//    }
//    if (isum != num_ions_)
//    {
//        (*MPIdata::sout) << " num_ions != isum !!!!" << endl;
//        ct.global_exit(2);
//    }
//
//#ifdef DEBUG
//    mmpi.allreduce(&isum2, 1, MPI_SUM);
//    (*MPIdata::sout) << " Isum2=" << isum2 << endl;
//    assert(isum2 == num_ions_);
//#endif
//
//    // if( onpe0 )
//    //    (*MPIdata::sout)<<"Ions::associate2PE() done..."<<endl;
//}

// Writes out the positions of all the ions
// PE root does the writing
void Ions::printPositionsGlobal(ostream& os, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    map<int, string> ion_names;
    vector<double> positions;
    vector<int> islocked;

    gatherNames(ion_names, root, mmpi.commSameSpin());
    gatherPositions(positions, root);
    gatherLockedData(islocked, root);

    if (mmpi.mypeGlobal() == root)
    {

        os << endl << " IONIC POSITIONS:" << endl;

        os << setw(8) << "Atoms" << setw(8) << "X" << setw(10) << "Y"
           << setw(10) << "Z" << endl;

        os.setf(ios::right, ios::adjustfield);
        os.setf(ios::fixed, ios::floatfield);

        int ion_index                           = 0;
        map<int, string>::const_iterator ion_id = ion_names.begin();
        while (ion_id != ion_names.end())
        {
            const int pos = 3 * ion_index;

            os << "$$ ";
            if (islocked[ion_index])
                os << "*";
            else
                os << " ";
            os << setw(4) << ion_id->second << setw(10) << setprecision(4)
               << fixed << positions[pos] << setw(10) << positions[pos + 1]
               << setw(10) << positions[pos + 2] << endl;

            ion_index++;
            ion_id++;
        }

        os << endl;
    }
}

// Writes out the postions of the local ions and their displacements from their
// initial postions.
void Ions::printPositionsLocal(ostream& os, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    if (mmpi.mypeGlobal() == root)
    {

        os << endl
           << " IONIC POSITIONS AND DISPLACEMENTS ON PE" << root << ":" << endl;

        os << setw(8) << "Atoms" << setw(8) << "X" << setw(10) << "Y"
           << setw(10) << "Z" << setw(10) << "dX" << setw(10) << "dY"
           << setw(10) << "dZ" << endl;

        os.setf(ios::right, ios::adjustfield);
        os.setf(ios::fixed, ios::floatfield);

        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            (*ion)->printPosition(os);
            ion++;
        }

        os << endl;
    }
}
// Writes out the postions of the ions and their displacements from their
// initial postions.
void Ions::printPositions(ostream& os, const int root) const
{
    Control& ct(*(Control::instance()));
    if (ct.verbose > 2)
    {
        if (num_ions_ < 256)
            printPositionsGlobal(os, root);
        else
            printPositionsLocal(os, root);
    }
}

void Ions::writeAtomicNumbers(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::writeAtomicNumbers()..." << endl;
    }
    vector<int> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherAtomicNumbers(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            assert((*ion)->atomic_number() > 0);
            assert((*ion)->atomic_number() < 200);

            data.push_back((*ion)->atomic_number());
            ion++;
        }
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short i = s; i < ms; i++)
                data.push_back(-1);
        }
        // for(vector<int>::iterator it =data.begin();
        //                          it!=data.end();
        //                        ++it)
        //    cout<<"Number="<<*it<<endl;

        // cout<<"ms="<<ms<<endl;

        // for(short i=0;i<ms;i++)data.push_back(-100*(mytask+1));
        // vector<int> old(data);
        // data.clear();
        // for(short i=0;i<ms;i++)
        //{
        //    data.push_back(old[i]);
        //    data.push_back(-1);
        //}

        size_t dims[2] = { data.size(), 1 };

        string datasetname("/Atomic_numbers");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

void Ions::writeAtomNames(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::writeAtomNames" << endl;
    }

    // gather data to print locally into vector "data"
    vector<string> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherNames(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            data.push_back((*ion)->name());
            ion++;
        }
    }

    // write data
    hid_t file_id = h5f_file.file_id();

    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            string empty_string;
            for (short i = s; i < ms; i++)
                data.push_back(empty_string);
        }

        size_t dims[2] = { data.size(), 1 };

        string datasetname("/Atomic_names");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

void Ions::lockAtom(const string& name)
{
    vector<Ion*>::iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        string name_ion((*ion)->name());
        if (name.compare(name_ion) == 0)
        {
            (*ion)->lock();
            if (onpe0) (*MPIdata::sout) << "Lock atom " << name << endl;
            break;
        }
        ion++;
    }
}

void Ions::readLockedAtomNames(HDFrestart& h5f_file)
{
    int dim = 0;

    if (dim == 0) return;

    vector<string> data;
    h5f_file.readLockedAtomNames(data);

    for (vector<string>::const_iterator i = data.begin(), end = data.end();
         i != end; ++i)
    {
        lockAtom(*i);
    }
}

void Ions::writeLockedAtomNames(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::writeLockedAtomsNames" << endl;
    }

    // gather data to print locally into vector "data"
    vector<string> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherLockedNames(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            if ((*ion)->locked()) data.push_back((*ion)->name());
            ion++;
        }
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            string empty_string;
            for (short i = s; i < ms; i++)
                data.push_back(empty_string);
            if (ms == 0) return;
        }

        size_t dims[2] = { data.size(), 1 };

        string datasetname("/LockedAtomsNames");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

void Ions::writeAtomicIDs(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));
    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::writeAtomicIDs()..." << endl;
    }

    // gather data to print locally into vector "data"
    vector<int> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherIndexes(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            data.push_back((*ion)->index());
            ion++;
        }
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short i = s; i < ms; i++)
                data.push_back(-1);
        }

        size_t dims[2] = { data.size(), 1 };

        string datasetname("/Atomic_IDs");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

void Ions::writeAtomicNLprojIDs(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::writeAtomicNLprojIDs()..." << endl;
    }

    // gather data to print locally into vector "data"
    vector<int> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherNLprojIds(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            data.push_back((*ion)->nlprojid());
            ion++;
        }
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short i = s; i < ms; i++)
                data.push_back(-1);
        }

        size_t dims[2] = { data.size(), 1 };

        string datasetname("/AtomicNLproj_IDs");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

void Ions::writePositions(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::writePositions" << endl;
    }

    vector<double> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherPositions(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            data.push_back((*ion)->position(0));
            data.push_back((*ion)->position(1));
            data.push_back((*ion)->position(2));
            ion++;
        }
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short i = s; i < ms; i++)
                data.push_back(1.e32);
        }

        size_t dims[2] = { data.size() / 3, 3 };

        string datasetname("/Ionic_positions");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

void Ions::initFromRestartFile(HDFrestart& h5_file)
{
    assert(list_ions_.empty());
    assert(local_ions_.empty());

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    Control& ct(*(Control::instance()));
    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::setFromRestartFile()..." << endl;
    }

    // set up list boundaries
    // get radius of projectors
    double rmax = getMaxListRadius();
    assert(rmax > 0.);
    setupListIonsBoundaries(rmax);

    vector<int> at_numbers;
    h5_file.readAtomicNumbers(at_numbers);
    vector<int> at_indexes;
    int nidxs = h5_file.readAtomicIDs(at_indexes);
    vector<int> at_nlprojIds;
    int npids = h5_file.readAtomicNLprojIDs(at_nlprojIds);
    vector<string> at_names;
    h5_file.readAtomicNames(at_names);
    if (onpe0 && ct.verbose > 2)
    {
        cout << "HDF file: at nb=" << at_numbers.size() << endl;
        cout << "HDF file: names nb=" << at_names.size() << endl;
        cout << "HDF file: indexes nb=" << at_indexes.size() << endl;
        cout << "HDF file: at_nlprojIds nb=" << at_nlprojIds.size() << endl;
    }

    // if reading "old" format with replicated atoms, fill up empty arrays
    // (atomic indexes were not saved before)
    if (nidxs == -1 && at_numbers.size() > 0)
    {
        for (int i = 0; i < at_numbers.size(); i++)
            at_indexes.push_back(i);
    }
    if (npids == -1 && at_numbers.size() > 0)
    {
        for (int i = 0; i < at_numbers.size(); i++)
            at_nlprojIds.push_back(i);
    }

    assert(at_numbers.size() == at_names.size());
    assert(at_numbers.size() == at_indexes.size());
    assert(at_numbers.size() == at_nlprojIds.size());

    num_ions_ = at_names.size();
    mmpi.allreduce(&num_ions_, 1, MPI_SUM);

    if (onpe0 && ct.verbose > 0)
    {
        (*MPIdata::sout) << "Ions::setFromRestartFile(), read " << num_ions_
                         << " names... on PE0" << flush << endl;
    }

    assert(at_numbers.size() == at_names.size());

    double default_vel[3]    = { 0., 0., 0. };
    double default_coords[3] = { 0., 0., 0. };

    for (int i = 0; i < at_numbers.size(); i++)
    {
        const int atnum                     = at_numbers[i];
        vector<Species>::const_iterator spi = species_.begin();
        while (spi != species_.end())
        {
            if (atnum == spi->getAtomicNumber()) break;
            spi++;
        }
        if (onpe0 && ct.verbose > 3)
            (*MPIdata::sout) << "New Ion with name " << at_names[i] << endl;
        if (at_indexes[i] >= 0)
        {
            Ion* new_ion = new Ion(*spi, at_names[i], default_coords,
                default_vel, at_indexes[i], at_nlprojIds[i]);

            list_ions_.push_back(new_ion);
            local_ions_.push_back(new_ion);
        }
    }
    readRestartPositions(h5_file);
    readRestartVelocities(h5_file);
    readRestartRandomStates(h5_file);
    readLockedAtomNames(h5_file);

    // rescale all velocities by factor specified in input
    rescaleVelocities(ct.VelocityScalingFactor());

    setup_ = false;

    // update list ions
    updateListIons();

    //#ifndef NDEBUG
    //    checkUnicityLocalIons();
    //#endif
}

void Ions::readRestartPositions(HDFrestart& h5_file)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Read ionic positions from hdf5 file" << endl;

    vector<double> data;
    h5_file.readAtomicPositions(data);

    int i                      = 0;
    vector<Ion*>::iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        (*ion)->setPosition(data[3 * i], data[3 * i + 1], data[3 * i + 2]);
        ion++;
        i++;
    }
}

void Ions::writeVelocities(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::writeVelocities" << endl;
    }

    vector<double> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherVelocities(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            data.push_back((*ion)->velocity(0));
            data.push_back((*ion)->velocity(1));
            data.push_back((*ion)->velocity(2));
            ion++;
        }
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short i = s; i < ms; i++)
                data.push_back(1.e32);
        }

        size_t dims[2] = { data.size() / 3, 3 };

        string datasetname("/Ionic_velocities");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

void Ions::writeRandomStates(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::writeRandomStates()..." << endl;
    }

    vector<unsigned short> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherRandStates(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            data.push_back((*ion)->randomState(0));
            data.push_back((*ion)->randomState(1));
            data.push_back((*ion)->randomState(2));
            ion++;
        }
    }

    if (!data.empty())
        if (data[0] != data[0])
        {
            (*MPIdata::sout)
                << "WARNING: Ions::writeRandomStates: data[0]=" << data[0]
                << endl;
        }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        if (h5f_file.useHdf5p())
        {
            short s = data.size();
            short ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short i = s; i < ms; i++)
                data.push_back(0);
        }

        size_t dims[2] = { data.size() / 3, 3 };

        string datasetname("/Ionic_RandomStates");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

void Ions::removeMassCenterMotion()
{
    // don't do it if some atoms are locked
    if (has_locked_atoms_) return;

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Remove mass center motion" << endl;

    vector<double> mass(local_ions_.size());
    vector<double> velocities(3 * local_ions_.size());

    vector<Ion*>::iterator ion = local_ions_.begin();
    int i                      = 0;
    double tmp[4]              = { 0., 0., 0., 0. };
    while (ion != local_ions_.end())
    {
        const int threei       = 3 * i;
        velocities[threei]     = (*ion)->velocity(0);
        velocities[threei + 1] = (*ion)->velocity(1);
        velocities[threei + 2] = (*ion)->velocity(2);
        mass[i]                = (*ion)->getMass();
        tmp[0] += mass[i] * velocities[threei];
        tmp[1] += mass[i] * velocities[threei + 1];
        tmp[2] += mass[i] * velocities[threei + 2];
        tmp[3] += mass[i];

        ion++;
        i++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&tmp[0], 4, MPI_SUM);

    const double tmass = 1. / tmp[3];
    for (short j = 0; j < 3; j++)
    {
        tmp[j] *= tmass;
#ifdef DEBUG
        (*MPIdata::sout) << setprecision(12);
        if (onpe0) (*MPIdata::sout) << "V[" << j << "]=" << tmp[j] << endl;
#endif
    }

    ion = local_ions_.begin();
    i   = 0;
    while (ion != local_ions_.end())
    {
        const int threei = 3 * i;
        (*ion)->setVelocity(velocities[threei] - tmp[0],
            velocities[threei + 1] - tmp[1], velocities[threei + 2] - tmp[2]);
        ion++;
        i++;
    }

#ifdef DEBUG // check velocity mass center is 0
    mv[0] = 0.;
    mv[1] = 0.;
    mv[2] = 0.;
    while (ion != local_ions_.end())
    {
        velocities[3 * i]     = (*ion)->velocity(0);
        velocities[3 * i + 1] = (*ion)->velocity(1);
        velocities[3 * i + 2] = (*ion)->velocity(2);
        mv[0] += mass[i] * velocities[3 * i];
        mv[1] += mass[i] * velocities[3 * i + 1];
        mv[2] += mass[i] * velocities[3 * i + 2];

        ion++;
        i++;
    }

    for (short j = 0; j < 3; j++)
    {
        mv[j] *= tmass;
        if (onpe0) (*MPIdata::sout) << "V[" << j << "]=" << mv[j] << endl;
    }
#endif
}

void Ions::readRestartVelocities(HDFrestart& h5_file)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Read atomic velocities from hdf5 file" << endl;

    vector<double> data;
    h5_file.readAtomicVelocities(data);

    vector<Ion*>::iterator ion = local_ions_.begin();
    int i                      = 0;
    while (ion != local_ions_.end())
    {
        (*ion)->setVelocity(data[3 * i], data[3 * i + 1], data[3 * i + 2]);
        ion++;
        i++;
    }
}

void Ions::readRestartRandomStates(HDFrestart& h5f_file)
{
    if (onpe0)
        (*MPIdata::sout) << "Read atomic RandomStates from hdf5 file" << endl;

    vector<unsigned short> data;
    h5f_file.readRestartRandomStates(data);

    vector<Ion*>::iterator ion = local_ions_.begin();
    int i                      = 0;
    while (ion != local_ions_.end())
    {
        (*ion)->setRandomState(data[3 * i], data[3 * i + 1], data[3 * i + 2]);
        ion++;
        i++;
    }
}

void Ions::writeForces(HDFrestart& h5f_file)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "Write ionic forces in hdf5 file" << endl;

    vector<double> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherForces(data, 0, comm);
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions().end())
        {
            // get position of local ion
            double force[3];
            (*ion)->getForce(&force[0]);
            data.push_back(force[0]);
            data.push_back(force[1]);
            data.push_back(force[2]);

            ++ion;
        }
    }

    hid_t file_id = h5f_file.file_id();
    if (file_id >= 0)
    {
        // fill up data array to dimension common to all tasks
        short ms = data.size();
        if (h5f_file.useHdf5p())
        {
            short s = ms;
            mgmol_tools::allreduce(&s, &ms, 1, MPI_MAX, h5f_file.comm_active());
            for (short i = s; i < ms; i++)
                data.push_back(1.e32);
        }

        size_t dims[2] = { data.size() / 3, 3 };

        string datasetname("/Ionic_forces");
        if (h5f_file.useHdf5p())
        {
            mgmol_tools::parallelWrite2d(
                file_id, datasetname, data, dims, h5f_file.comm_active());
        }
        else
        {
            mgmol_tools::write2d(file_id, datasetname, data, dims);
        }
    }
}

// Writes out the postions of the ions and the current forces on them by root
void Ions::printForcesGlobal(ostream& os, const int root) const
{
    double maxf = 0., avfx = 0., avfy = 0., avfz = 0., maxfx = 0., maxfy = 0.,
           maxfz  = 0.;
    int num_atoms = 0;

    double sum_forces[3] = { 0., 0., 0. };
    int num_movable      = 0;

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    map<int, string> ion_names;
    vector<double> positions;
    vector<double> forces;
    vector<int> islocked;

    gatherNames(ion_names, root, mmpi.commSpin());
    gatherPositions(positions, root);
    gatherForces(forces, root);
    gatherLockedData(islocked, root);

    if (mmpi.mypeGlobal() == root)
    {
        os << "IONIC POSITIONS AND FORCES:" << endl;

        os << setiosflags(ios::left) << setw(8) << "Atoms"
           << resetiosflags(ios::left) << setw(8) << "X" << setw(10) << "Y"
           << setw(10) << "Z" << setw(10) << "FX" << setw(10) << "FY"
           << setw(10) << "FZ" << endl;

        int ion_index                           = 0;
        map<int, string>::const_iterator ion_id = ion_names.begin();
        while (ion_id != ion_names.end())
        {

            const int pos = 3 * ion_index;
            os << "## ";
            if (islocked[ion_index])
                os << "*";
            else
                os << " ";
            os << setw(4) << ion_id->second << setiosflags(ios::right)
               << setw(10) << setprecision(4) << fixed << positions[pos]
               << setw(10) << positions[pos + 1] << setw(10)
               << positions[pos + 2] << setprecision(7) << scientific
               << setw(16) << forces[pos] << setw(16) << forces[pos + 1]
               << setw(16) << forces[pos + 2] << endl;

            if (!islocked[ion_index])
            {

                avfx += fabs(forces[pos]);
                avfy += fabs(forces[pos + 1]);
                avfz += fabs(forces[pos + 2]);

                double ff = forces[pos] * forces[pos]
                            + forces[pos + 1] * forces[pos + 1]
                            + forces[pos + 2] * forces[pos + 2];

                maxf = max(maxf, ff);

                maxfx = max(maxfx, fabs(forces[pos]));
                maxfy = max(maxfy, fabs(forces[pos + 1]));
                maxfz = max(maxfz, fabs(forces[pos + 2]));

                num_movable++;
            }

            for (short ii = 0; ii < 3; ii++)
                sum_forces[ii] += forces[pos + ii];

            num_atoms++;
            ion_index++;
            ion_id++;
        }

        if (num_atoms == 0) return;
        // global statistics
        os << endl << "Global Statistics:" << endl;
        os << "==========================" << endl;

        const double inv_num_atoms = 1. / (double)num_atoms;
        for (short ii = 0; ii < 3; ii++)
            sum_forces[ii] *= inv_num_atoms;

        if (num_movable)
        {
            double inv_num_movable = 1. / (double)num_movable;
            avfx                   = avfx * inv_num_movable;
            avfy                   = avfy * inv_num_movable;
            avfz                   = avfz * inv_num_movable;

            os << scientific;
            os << " mean F on movable ions  = (" << avfx << "," << avfy << ","
               << avfz << ")" << endl;
            os << " max F on movable ions   = (" << maxfx << "," << maxfy << ","
               << maxfz << ")" << endl;
            os << " max |F| on movable ions = " << sqrt(maxf) << endl;
        }
        os << " Sum forces on all ions  = (" << sum_forces[0] << ","
           << sum_forces[1] << "," << sum_forces[2] << ")" << endl;
        os << fixed;
    }
}

// Writes out the postions of the local ions and the current forces on them
void Ions::printForcesLocal(ostream& os, const int root) const
{
    int num_atoms   = 0;
    int num_movable = 0;

    double buf[10];
    for (int i = 0; i < 10; i++)
        buf[i] = 0.;

    double* sum_forces = &buf[0];
    double* avg_forces = sum_forces + 3;
    double* max_forces = avg_forces + 3;
    double* maxf       = max_forces + 3;

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.barrier();

    if (mmpi.mypeGlobal() == root)
    {
        os << endl
           << endl
           << "IONIC POSITIONS AND FORCES ON CENTERED ON PE" << root << ":"
           << endl;

        os << setiosflags(ios::left) << setw(8) << "Atoms"
           << resetiosflags(ios::left) << setw(8) << "X" << setw(10) << "Y"
           << setw(10) << "Z" << setw(10) << "FX" << setw(10) << "FY"
           << setw(10) << "FZ" << endl;

        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {

            (*ion)->printPositionAndForce(os);

            if (!(*ion)->locked())
            {

                avg_forces[0] += fabs((*ion)->force(0));
                avg_forces[1] += fabs((*ion)->force(1));
                avg_forces[2] += fabs((*ion)->force(2));

                double ff = (*ion)->norm2F();
                maxf[0]   = max(maxf[0], ff);

                max_forces[0] = max(max_forces[0], fabs((*ion)->force(0)));
                max_forces[1] = max(max_forces[1], fabs((*ion)->force(1)));
                max_forces[2] = max(max_forces[2], fabs((*ion)->force(2)));

                num_movable++;
            }

            for (short ii = 0; ii < 3; ii++)
                sum_forces[ii] += (*ion)->force(ii);
            num_atoms++;

            ion++;
        }
    }
    else
    {
        vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {

            if (!(*ion)->locked())
            {

                avg_forces[0] += fabs((*ion)->force(0));
                avg_forces[1] += fabs((*ion)->force(1));
                avg_forces[2] += fabs((*ion)->force(2));

                double ff = (*ion)->norm2F();
                maxf[0]   = max(maxf[0], ff);

                max_forces[0] = max(max_forces[0], fabs((*ion)->force(0)));
                max_forces[1] = max(max_forces[1], fabs((*ion)->force(1)));
                max_forces[2] = max(max_forces[2], fabs((*ion)->force(2)));

                num_movable++;
            }

            for (short ii = 0; ii < 3; ii++)
                sum_forces[ii] += (*ion)->force(ii);
            num_atoms++;

            ion++;
        }
    }
    // global statistics
    mmpi.allreduce(&num_atoms, 1, MPI_SUM);
    if (num_atoms == 0) return;
    if (onpe0)
    {
        os << endl << "Global Statistics:" << endl;
        os << "==========================" << endl;
    }

    mmpi.allreduce(&num_movable, 1, MPI_SUM);
    mmpi.allreduce(&buf[0], 6, MPI_SUM);
    mmpi.allreduce(&buf[6], 4, MPI_MAX);

    const double inv_num_atoms = 1. / (double)num_atoms;
    for (short ii = 0; ii < 3; ii++)
        sum_forces[ii] *= inv_num_atoms;

    if (mmpi.mypeGlobal() == root)
    {
        if (num_movable)
        {
            double inv_num_movable = 1. / (double)num_movable;
            avg_forces[0]          = avg_forces[0] * inv_num_movable;
            avg_forces[1]          = avg_forces[1] * inv_num_movable;
            avg_forces[2]          = avg_forces[2] * inv_num_movable;

            os << scientific;
            os << " mean F on movable ions  = (" << avg_forces[0] << ","
               << avg_forces[1] << "," << avg_forces[2] << ")" << endl;
            os << " max F on movable ions   = (" << max_forces[0] << ","
               << max_forces[1] << "," << max_forces[2] << ")" << endl;
            os << " max |F| on movable ions = " << sqrt(maxf[0]) << endl;
        }
        os << " Sum forces on all ions  = (" << sum_forces[0] << ","
           << sum_forces[1] << "," << sum_forces[2] << ")" << endl;
        os << fixed;
    }
}

// Writes out the postions of the ions and the current forces on them
void Ions::printForces(ostream& os, const int root) const
{

    Control& ct(*(Control::instance()));

    if (num_ions_ < 512 || ct.verbose > 2)
        printForcesGlobal(os, root);
    else
        printForcesLocal(os, root);
}

int Ions::countIonsHere() const { return (int)local_ions_.size(); }

int Ions::countProjectorsHere() const
{
    int count                        = 0;
    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        count += (*ion)->nProjectors();
        ion++;
    }
    return count;
}

int Ions::countProjectors() const
{
    //    assert( setup_ );
    // cout<<"Num. local ions: "<<local_ions_.size()<<endl;

    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    int nproj                        = 0;
    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        nproj += (*ion)->nProjectors();
        ion++;
    }
    int tmp = nproj;
    MPI_Allreduce(&tmp, &nproj, 1, MPI_INT, MPI_SUM, myPEenv.comm());

    return nproj;
}

int Ions::countProjectorsSubdomain() const
{
    int nproj                        = 0;
    vector<Ion*>::const_iterator ion = overlappingNL_ions_.begin();
    while (ion != overlappingNL_ions_.end())
    {
        nproj += (*ion)->nProjectorsSubdomain();
        ion++;
    }
    return nproj;
}

Ion* Ions::findIon(const string& name) const
{
    vector<Ion*>::const_iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        bool same_name = (*ion)->compareName(name);
        if (same_name) return (*ion);

        ion++;
    }

    // return 0 if name not found in list_ions_
    return 0;
}

bool Ions::isLocal(const string& name) const
{
    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        bool same_name = (*ion)->compareName(name);
        if (same_name) return true;

        ion++;
    }
    return false;
}

Ion* Ions::findIon(const int index) const
{
    vector<Ion*>::const_iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        if ((*ion)->compareIndex(index)) return (*ion);

        ion++;
    }
    return NULL;
}

Ion* Ions::findLocalIon(const int index) const
{
    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        if ((*ion)->compareIndex(index)) return (*ion);

        ion++;
    }
    return NULL;
}

void Ions::setPositions(const vector<double>& tau)
{
    assert(tau.size() == 3 * local_ions_.size());

    vector<Ion*>::iterator ion = local_ions_.begin();
    int ia                     = 0;
    while (ion != local_ions_.end())
    {
        (*ion)->setPosition(tau[3 * ia + 0], tau[3 * ia + 1], tau[3 * ia + 2]);

        ion++;
        ia++;
    }

    setup_ = false;
}

void Ions::get_positions(vector<vector<double>>& rr) const
{
    assert(rr.size() == local_ions_.size());
    if (local_ions_.empty()) return;
    vector<double> tau(3);
    int i                            = 0;
    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        tau[0] = (*ion)->position(0);
        tau[1] = (*ion)->position(1);
        tau[2] = (*ion)->position(2);
        rr[i]  = tau;
        ion++;
        i++;
    }
}
void Ions::set_positions(const vector<vector<double>>& rr)
{
    assert(rr.size() == local_ions_.size());

    if (local_ions_.empty()) return;
    vector<Ion*>::iterator ion = local_ions_.begin();
    int i                      = 0;
    while (ion != local_ions_.end())
    {
        assert(rr[i].size() == 3);
        (*ion)->setPosition(rr[i][0], rr[i][1], rr[i][2]);
        ion++;
        i++;
    }
}
void Ions::get_forces(vector<vector<double>>& ff) const
{
    assert(ff.size() == local_ions_.size());

    if (local_ions_.empty()) return;
    int i                            = 0;
    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        ff[i][0] = (*ion)->force(0);
        ff[i][1] = (*ion)->force(1);
        ff[i][2] = (*ion)->force(2);
        ion++;
        i++;
    }
}
void Ions::set_forces(const vector<vector<double>>& ff)
{
    assert(ff.size() == local_ions_.size());

    if (local_ions_.empty()) return;
    vector<Ion*>::iterator ion = local_ions_.begin();
    int i                      = 0;
    while (ion != local_ions_.end())
    {
        (*ion)->setForce(ff[i][0], ff[i][1], ff[i][2]);
        ion++;
        i++;
    }
}

int Ions::readAtoms(const string& filename, const bool cell_relative)
{

    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Ions::readAtoms() --- Read "
                         << " atomic positions from file " << filename << endl;

    string strxyz(".xyz");
    size_t found = filename.find(strxyz);
    if (found != string::npos)
    {
        num_ions_ = readAtomsFromXYZ(filename, cell_relative);
    }
    else
    {
        num_ions_ = readNatoms(filename, cell_relative);
    }

    return num_ions_;
}

int Ions::readAtoms(ifstream* tfile, const bool cell_relative)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    string first_string_read = "";
    if (mmpi.PE0())
    {
        (*tfile) >> first_string_read;
    }
    mmpi.bcastGlobal(first_string_read);

    string strxyz(".xyz");
    size_t found = first_string_read.find(strxyz);
    if (found != string::npos)
    {
        if (mmpi.PE0())
            (*MPIdata::sout) << "Read atomic positions from file "
                             << first_string_read << endl;
        num_ions_ = readAtomsFromXYZ(first_string_read, cell_relative);
    }
    else
    {
        if (mmpi.PE0())
        {
            for (short i = 0; i < first_string_read.size(); i++)
                tfile->unget();
        }
        num_ions_ = readNatoms(tfile, cell_relative);
    }

    return num_ions_;
}

int Ions::readAtomsFromXYZ(const string& filename, const bool cell_relative)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    Control& ct(*(Control::instance()));

    // set up list boundaries
    // get radius of projectors
    double rmax = getMaxListRadius();
    if (onpe0)
        (*MPIdata::sout) << "Max. radius for species in XYZ file: " << rmax
                         << endl;
    assert(rmax > 0.);

    setupListIonsBoundaries(rmax);

    int natoms = -1;

    ifstream* tfile = 0;
    if (mmpi.PE0())
    {
        tfile = new ifstream(filename.data(), ios::in);
        if (!(*tfile))
        {
            (*MPIdata::serr)
                << "Ions::readAtomsFromXYZ() --- ERROR: cannot open "
                << filename << endl;
        }
        else
        {
            string query;
            if (getline(*tfile, query))
            {
                stringstream na(query);
                na >> natoms;
            }
        }
        if (natoms < 1)
            (*MPIdata::sout)
                << "WARNING: Ions::readAtomsFromXYZ(), number of atoms read = "
                << natoms << endl;
    }
    mmpi.bcastGlobal(&natoms);
    if (natoms < 0) return natoms;

    vector<double> crds(3 * natoms);
    vector<short> spec(natoms);

    int count = 0;
    if (mmpi.PE0())
    {
        string query;
        getline(*tfile, query); // read comment line
        // read atomic species and positions
        for (int ia = 0; ia < natoms; ++ia)
        {
            if (!getline(*tfile, query))
            {
                break;
            }
            stringstream ss(query);
            string name_read;
            ss >> name_read;
            spec[ia] = map_species_.find(name_read)->second;
            for (int j = 0; j < 3; j++)
            {
                ss >> crds[3 * ia + j];
                crds[3 * ia + j] *= ang2bohr;
                assert(crds[3 * ia + j] > -100000.);
                assert(crds[3 * ia + j] < 100000.);
            }
            if (cell_relative)
            {
                crds[3 * ia + 0] *= lattice_[0];
                crds[3 * ia + 1] *= lattice_[1];
                crds[3 * ia + 2] *= lattice_[2];
            }
            ++count;
        }
    }

    mmpi.bcastGlobal(&count, 1);
    if (count < natoms) return count;

    mmpi.bcastGlobal(&crds[0], 3 * natoms);
    mmpi.bcastGlobal(&spec[0], natoms);

    double velocity[3] = { 0., 0., 0. };
    bool locked        = false;
    for (int ia = 0; ia < natoms; ++ia)
    {
        vector<Species>::const_iterator it = species_.begin();
        int isp                            = -1;
        while (it != species_.end())
        {
            ++isp;
            if (it->getAtomicNumber() == spec[ia])
            {
                break;
            }
            ++it;
        }
        string spname("");
        for (map<string, short>::iterator itr = map_species_.begin();
             itr != map_species_.end(); ++itr)
        {
            if (itr->second == spec[ia])
            {
                spname = itr->first;
                break;
            }
        }
        if (spname.compare("") == 0)
        {
            (*MPIdata::serr) << "Ions::readAtomsFromXYZ() --- ERROR: unknown "
                                "species for atomic number "
                             << spec[ia] << endl;
            return -1;
        }

        // make a name for atom based on species and order of reading in
        string aname(spname);
        stringstream ss;
        ss << ia;
        if (ia < 10) aname.append("0");
        if (ia < 100) aname.append("0");
        if (ia < 1000) aname.append("0");
        if (ia < 10000) aname.append("0");
        aname.append(ss.str());
        Ion* new_ion
            = new Ion(species_[isp], aname, &crds[3 * ia], velocity, locked);
#ifdef USE_MPI
        new_ion->bcast(mmpi.commGlobal());
#endif

        // Populate list_ions_ list
        // cout<<"crds: "<<crds[3*ia+0]<<", "<<crds[3*ia+1]<<",
        // "<<crds[3*ia+2]<<endl;
        if (inListIons(crds[3 * ia + 0], crds[3 * ia + 1], crds[3 * ia + 2]))
        {
            list_ions_.push_back(new_ion);
            if (ct.verbose > 2)
                (*MPIdata::sout)
                    << "Ion " << aname << " at position " << crds[3 * ia + 0]
                    << "," << crds[3 * ia + 1] << "," << crds[3 * ia + 2]
                    << " added to the list... on PE" << mmpi.mypeGlobal()
                    << endl;
            // populate local_ions_ list
            if (inLocalIons(
                    crds[3 * ia + 0], crds[3 * ia + 1], crds[3 * ia + 2]))
            {
                (new_ion)->set_here(true);
                local_ions_.push_back(new_ion);
                if (onpe0 && ct.verbose > 2)
                    (*MPIdata::sout)
                        << "Ion " << aname << " at position "
                        << crds[3 * ia + 0] << "," << crds[3 * ia + 1] << ","
                        << crds[3 * ia + 2]
                        << " added to the list of local ions... on PE"
                        << mmpi.mypeGlobal() << endl;
            }
            else
                (new_ion)->set_here(false);
        }
        else
        {
            //(*MPIdata::sout)<<"Ion "<<aname<<" NOT added to the list... on
            // PE"<<mmpi.mype() <<endl;
            delete new_ion;
        }
    }
    //    cout<<mmpi.mype()<<"...list size = "<<list_ions_.size()<<" local ions
    //    size = "<<local_ions_.size()<<endl;

    if (mmpi.PE0()) delete tfile;

    return natoms;
}

int Ions::readNatoms(const string& filename, const bool cell_relative)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Ions::readNAtoms() --- Read "
                         << " atomic positions from file " << filename << endl;

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    ifstream* tfile = 0;
    if (mmpi.instancePE0())
    {
        tfile = new ifstream(filename.data(), ios::in);
        if (!tfile->is_open())
        {
            (*MPIdata::serr)
                << " Unable to open file " << filename.data() << endl;
            return -1;
        }
        else
        {
            (*MPIdata::sout) << "Open " << filename.data() << endl;
        }
    }

    int nread = readNatoms(tfile, cell_relative);

    if (mmpi.instancePE0())
    {
        if (ct.verbose > 0)
            (*MPIdata::sout) << "Close " << filename.data() << endl;
        tfile->close();
        delete tfile;
    }

    return nread;
}

int Ions::readNatoms(ifstream* tfile, const bool cell_relative)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    if (mmpi.PE0()) assert(tfile != 0);

    // set pointer to species object
    // set up list boundaries
    // get radius of projectors
    double rmax = getMaxListRadius();
    assert(rmax > 0.);
    setupListIonsBoundaries(rmax);

    // make sure all processors call this function
    mmpi.barrier();

    Control& ct(*(Control::instance()));
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Ions::readNatoms() --- Try to read "
                         << " atomic positions..." << endl;
    int nread = 0;
    int iread = 1;
    while (iread > 0)
    {
        iread = read1atom(tfile, cell_relative);
        // if(onpe0)cout<<"iread="<<iread<<endl;
        if (iread < 0) return -100 - nread;
        nread += iread;
        if (onpe0 && ct.verbose > 0 && (nread % 1000 == 0))
            (*MPIdata::sout) << "Ions::readNatoms() --- read " << nread
                             << " atomic positions..." << endl;
    }
    if (onpe0) read_comments(*tfile);

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Ions::readNatoms() --- read " << nread
                         << " atomic positions..." << endl;

    int siz = local_ions_.size();
    mmpi.allreduce(&siz, 1, MPI_SUM);
    assert(siz == nread);

    return nread;
}

// return number of atoms read, 0 if end of file, or -1 if failure occurs
int Ions::read1atom(ifstream* tfile, const bool cell_relative)
{
    string name_read = "";
    // short  isp=0;
    double crds[3];
    double velocity[3] = { 0., 0., 0. };

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    Control& ct(*(Control::instance()));

    short movable = 0;
    string query;
    short count = 1;
    if (mmpi.PE0())
    {
        if (getline(*tfile, query))
        {
            if (query[0] == '#') count = 0; // comment
            if (query.empty()) count = 0; // end of line only
            if (query.size() < 3) count = 0; // white spaces only
            // cout<<"length="<<query.size()<<endl;
        }
        else
        {
            if (tfile->eof())
            {
                //(*MPIdata::serr)<<"Ions::read1atom() --- eof..."<<endl;
                count = 0;
            }
            else
            {
                (*MPIdata::serr)
                    << "ERROR: Ions::read1atom() --- getline() failed..."
                    << endl;
                count = -1;
            }
        }
        // cout<<"query="<<query<<endl;
        // cout<<"count="<<count<<endl;
    }

    mmpi.bcastGlobal(&count, 1);
    if (count < 1) return count;

    if (mmpi.PE0())
    {
        stringstream ss(query);
        ss >> name_read;
        if (!checkValidName(name_read))
        {
            cerr << "ERROR: Invalid name read in input file: " << name_read
                 << endl;
            ct.global_exit(2);
        }
        short dummy;
        ss >> dummy; // not used anymore (was species index)

        for (short j = 0; j < 3; j++)
        {
            ss >> crds[j];
            assert(crds[j] > -10000.);
            assert(crds[j] < 10000.);
        }

        if (cell_relative)
        {
            crds[0] *= lattice_[0];
            crds[1] *= lattice_[1];
            crds[2] *= lattice_[2];
        }

        if (!(ss >> movable)) movable = 1; // default value
        if (movable != 0 && movable != 1)
        {
            (*MPIdata::serr) << "Atom " << name_read
                             << ", should be movable (1) or not(0)" << endl;
            return -1;
        }
#ifdef DEBUG
        (*MPIdata::sout) << "movable=" << movable << endl;
#endif
        int j = 0;
        if (movable)
        {
            while (ss >> velocity[j++])
                ;
        }

#ifndef NDEBUG
        (*MPIdata::sout) << "Ions::read1atom() --- Read Ion in position ("
                         << crds[0] << "," << crds[1] << "," << crds[2] << ")"
                         << " and velocity :(" << velocity[0] << ","
                         << velocity[1] << "," << velocity[2] << ")" << endl;
#endif

    } // if onpe0

    double tmp[7] = { crds[0], crds[1], crds[2], velocity[0], velocity[1],
        velocity[2], (double)movable };
    mmpi.bcastGlobal(tmp, 7);

    crds[0]     = tmp[0];
    crds[1]     = tmp[1];
    crds[2]     = tmp[2];
    velocity[0] = tmp[3];
    velocity[1] = tmp[4];
    velocity[2] = tmp[5];

    movable = (int)tmp[6];

    mmpi.bcastGlobal(name_read);

    // find species based on name
    string name(name_read);
    stripName(name);

    short spec_nb                      = map_species_.find(name)->second;
    vector<Species>::const_iterator it = species_.begin();
    int isp                            = -1;
    while (it != species_.end())
    {
        ++isp;
        if (it->getAtomicNumber() == spec_nb)
        {
            break;
        }
        ++it;
    }

    bool locked = (!movable);
    // if(isp>=(int)species_.size())
    //{
    //    (*MPIdata::serr)<<"name="<<name<<", isp="<<isp<<" >= num. species ("
    //                    <<species_.size()<<")"
    //        <<endl;
    //    return -1;
    //}

#ifdef DEBUG
    if (onpe0) (*MPIdata::sout) << "Create new Ion..." << endl;
#endif

    // Ion needs to be created by each MPI task to set global ids
    Ion* new_ion = new Ion(species_[isp], name_read, crds, velocity, locked);

    // assert( (new_ion->locked() == false) || (new_ion->locked() == true) );

#ifdef DEBUG
    if (onpe0) (*MPIdata::sout) << "Ion read..." << endl;
#endif

    // Populate list_ions_ list
    if (inListIons(crds[0], crds[1], crds[2]))
    {
        list_ions_.push_back(new_ion);
        // populate local_ions_ list
        if (inLocalIons(crds[0], crds[1], crds[2]))
        {
            (new_ion)->set_here(true);
            local_ions_.push_back(new_ion);
        }
        else
            (new_ion)->set_here(false);
    }
    else
    {
        delete new_ion;
    }

    return 1;
}

int Ions::getNValenceElectrons() const
{
    vector<Ion*>::const_iterator ion = local_ions_.begin();
    double val                       = 0.;
    while (ion != local_ions_.end())
    {
        val += (*ion)->getZion();
        ion++;
    }
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&val, 1, MPI_SUM);

    return (int)val;
}

double Ions::computeIonicCharge() const
{
    double ionic_charge              = 0.;
    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        ionic_charge += (double)(*ion)->getZion();
        ion++;
    }
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&ionic_charge, 1, MPI_SUM);
    return ionic_charge;
}

void Ions::setVelocities(
    const vector<double>& tau0, const vector<double>& taup, const double dt)
{
    assert(tau0.size() == 3 * local_ions_.size());
    assert(taup.size() == 3 * local_ions_.size());

    int ia                      = 0;
    vector<Ion*>::iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        double v[3];
        for (short i = 0; i < 3; i++)
        {
            v[i] = taup[3 * ia + i];
            v[i] -= tau0[3 * ia + i];
            v[i] /= dt;
        }
        (*iion)->setVelocity(v[0], v[1], v[2]);
        iion++;
        ia++;
    }
}

void Ions::getPositions(vector<double>& tau) const
{
    assert(tau.size() == 3 * local_ions_.size());

    int ia                            = 0;
    vector<Ion*>::const_iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        (*iion)->getPosition(&tau[3 * ia]);
        iion++;
        ia++;
    }
}

void Ions::setTau0()
{
    assert(tau0_.size() == 3 * local_ions_.size());

    int ia                            = 0;
    vector<Ion*>::const_iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        (*iion)->getPosition(&tau0_[3 * ia]);
        iion++;
        ia++;
    }
}

void Ions::setPositionsToTau0()
{
    assert(tau0_.size() == 3 * local_ions_.size());

    int ia                            = 0;
    vector<Ion*>::const_iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        (*iion)->setPosition(
            tau0_[3 * ia + 0], tau0_[3 * ia + 1], tau0_[3 * ia + 2]);
        iion++;
        ia++;
    }
}

void Ions::getVelocities(vector<double>& tau) const
{
    assert(tau.size() == 3 * local_ions_.size());

    int ia                            = 0;
    vector<Ion*>::const_iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        for (short i = 0; i < 3; i++)
        {
            tau[3 * ia + i] = (*iion)->velocity(i);
        }
        iion++;
        ia++;
    }
}

void Ions::getForces(vector<double>& tau) const
{
    assert(tau.size() == 3 * local_ions_.size());

    int ia                            = 0;
    vector<Ion*>::const_iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        assert(3 * ia + 2 < (int)tau.size());
        (*iion)->getForce(&tau[3 * ia]);
        iion++;
        ia++;
    }
}

void Ions::setMapVL()
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const double h0 = mygrid.hgrid(0);
    const double h1 = mygrid.hgrid(1);
    const double h2 = mygrid.hgrid(2);

    const double origin[3]
        = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };

    const int dim0  = mygrid.dim(0);
    const int dim1  = mygrid.dim(1);
    const int dim2  = mygrid.dim(2);
    const int gdim0 = mygrid.gdim(0);
    const int gdim1 = mygrid.gdim(1);
    const int gdim2 = mygrid.gdim(2);

    const int ilow = mygrid.istart(0);
    const int jlow = mygrid.istart(1);
    const int klow = mygrid.istart(2);
    const int ihi  = ilow + dim0 - 1;
    const int jhi  = jlow + dim1 - 1;
    const int khi  = klow + dim2 - 1;

    // Loop over list ions
    vector<Ion*>::iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        /* Generate range of indices over which the short-range difference */
        /* potential will be mapped onto the global grid.                  */
        (*ion)->set_lstart(0, origin[0], h0);
        (*ion)->set_lstart(1, origin[1], h1);
        (*ion)->set_lstart(2, origin[2], h2);

        // Generate indices
        vector<int> Ai0;
        (*ion)->get_Ai(Ai0, gdim0, 0);
        const int dimlx = Ai0.size();

        (*ion)->set_map_l(false);

        bool map0 = false;
        for (int idx = 0; idx < dimlx; idx++)
        {
            if ((Ai0[idx] >= ilow) && (Ai0[idx] <= ihi))
            {
                map0 = true;
                break;
            }
        }

        if (map0)
        {
            bool map1 = false;
            vector<int> Ai1;
            (*ion)->get_Ai(Ai1, gdim1, 1);
            const int dimly = Ai1.size();
            for (int idx = 0; idx < dimly; idx++)
            {
                if ((Ai1[idx] >= jlow) && (Ai1[idx] <= jhi))
                {
                    map1 = true;
                    break;
                }
            }
            if (map1)
            {
                bool map2 = false;
                vector<int> Ai2;
                (*ion)->get_Ai(Ai2, gdim2, 2);
                const int dimlz = Ai2.size();
                for (int idx = 0; idx < dimlz; idx++)
                {
                    if ((Ai2[idx] >= klow) && (Ai2[idx] <= khi))
                    {
                        map2 = true;
                        break;
                    }
                }

                if (map2)
                {
                    (*ion)->set_map_l(true);
                }
            }
        }
        ion++;
    }
}

double Ions::computeMaxVlRadius() const
{
    double radius = 0.;

    vector<Ion*>::const_iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        double r = (*iion)->computeRadiusVl();
        radius   = r > radius ? r : radius;
        iion++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&radius, 1, MPI_MAX);

    assert(radius >= 0.);

    return radius;
}

double Ions::computeMaxNLprojRadius() const
{
    double radius = 0.;

    vector<Ion*>::const_iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        double r = (*iion)->radiusNLproj();
        radius   = r > radius ? r : radius;
        iion++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&radius, 1, MPI_MAX);

    assert(radius >= 0.);

    return radius;
}

void Ions::gatherNames(
    map<int, string>& names, const int root, const MPI_Comm comm) const
{
    assert(comm != 0);

    vector<int> indexes(num_ions_, 0);
    vector<int> local_indexes;
    vector<string> data;
    vector<string> local_names;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        // get local name and index
        local_names.push_back((*ion)->name());
        local_indexes.push_back((*ion)->index());
        ++ion;
    }

    // gather data to PE root
    mgmol_tools::gatherV(local_names, data, root, comm);
    mgmol_tools::gatherV(local_indexes, indexes, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    if (mype == root)
    {
        int num_ions = data.size();
        assert(num_ions = indexes.size());
        for (int i = 0; i < num_ions; i++)
        {
            const int ion_index   = indexes[i];
            const string ion_name = data[i];
            names.insert(std::pair<int, string>(ion_index, ion_name));
        }
    }
}

void Ions::gatherNames(
    vector<string>& names, const int root, const MPI_Comm comm) const
{
    assert(comm != 0);

    vector<string> local_names;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        local_names.push_back((*ion)->name());
        ++ion;
    }

    // gather data to PE root
    vector<string> data;
    mgmol_tools::gatherV(local_names, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    names.clear();
    if (mype == root) names = data;
}

void Ions::gatherLockedNames(
    vector<string>& names, const int root, const MPI_Comm comm) const
{
    vector<string> local_names;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        if ((*ion)->locked()) local_names.push_back((*ion)->name());
        ++ion;
    }

    // gather data to PE root
    vector<string> data;
    mgmol_tools::gatherV(local_names, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    names.clear();
    if (mype == root) names = data;
}

void Ions::gatherIndexes(
    vector<int>& indexes, const int root, const MPI_Comm comm) const
{
    vector<int> local_indexes;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        local_indexes.push_back((*ion)->index());
        ++ion;
    }

    // gather data to PE root
    vector<int> data;
    mgmol_tools::gatherV(local_indexes, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    indexes.clear();
    if (mype == root) indexes = data;
}

void Ions::gatherNLprojIds(
    vector<int>& nlprojids, const int root, const MPI_Comm comm) const
{
    vector<int> local_nlprojids;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        local_nlprojids.push_back((*ion)->nlprojid());
        ++ion;
    }

    // gather data to PE root
    vector<int> data;
    mgmol_tools::gatherV(local_nlprojids, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    nlprojids.clear();
    if (mype == root) nlprojids = data;
}

void Ions::gatherAtomicNumbers(
    vector<int>& atnumbers, const int root, const MPI_Comm comm) const
{
    assert(comm != 0);

    vector<int> local_atnumbers;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        local_atnumbers.push_back((*ion)->atomic_number());
        ++ion;
    }

    // gather data to PE root
    vector<int> data;
    mgmol_tools::gatherV(local_atnumbers, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    atnumbers.clear();
    if (mype == root) atnumbers = data;
}

void Ions::gatherRandStates(
    vector<unsigned short>& rstates, const int root, const MPI_Comm comm) const
{
    vector<unsigned short> local_rstates;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        local_rstates.push_back((*ion)->randomState(0));
        local_rstates.push_back((*ion)->randomState(1));
        local_rstates.push_back((*ion)->randomState(2));
        ++ion;
    }

    // gather data to PE root
    vector<unsigned short> data;
    mgmol_tools::gatherV(local_rstates, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    rstates.clear();
    if (mype == root) rstates = data;
}

void Ions::gatherPositions(
    vector<double>& positions, const int root, const MPI_Comm comm) const
{
    vector<double> local_positions;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        // get position of local ion
        double position[3];
        (*ion)->getPosition(&position[0]);
        local_positions.push_back(position[0]);
        local_positions.push_back(position[1]);
        local_positions.push_back(position[2]);

        ++ion;
    }

    // gather data to PE root
    vector<double> data;
    mgmol_tools::gatherV(local_positions, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    positions.clear();
    if (mype == root) positions = data;
}

void Ions::gatherForces(
    vector<double>& forces, const int root, const MPI_Comm comm) const
{
    vector<double> local_forces;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        // get position of local ion
        double force[3];
        (*ion)->getForce(&force[0]);
        local_forces.push_back(force[0]);
        local_forces.push_back(force[1]);
        local_forces.push_back(force[2]);

        ++ion;
    }

    // gather data to PE root
    vector<double> data;
    mgmol_tools::gatherV(local_forces, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    forces.clear();
    if (mype == root) forces = data;
}

void Ions::gatherVelocities(
    vector<double>& velocities, const int root, const MPI_Comm comm) const
{
    vector<double> local_velocities;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        local_velocities.push_back((*ion)->velocity(0));
        local_velocities.push_back((*ion)->velocity(1));
        local_velocities.push_back((*ion)->velocity(2));

        ++ion;
    }

    // gather data to PE root
    vector<double> data;
    mgmol_tools::gatherV(local_velocities, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    velocities.clear();
    if (mype == root) velocities = data;
}

void Ions::gatherPositions(vector<double>& positions, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    positions.resize(3 * num_ions_, 0.);

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        // get local positions
        const int index = (*ion)->index();
        (*ion)->getPosition(&positions[3 * index]);
        ++ion;
    }

    // gather data to PE root
    mmpi.reduce(&positions[0], 3 * num_ions_, MPI_SUM, root);
}

void Ions::gatherForces(vector<double>& forces, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    forces.resize(3 * num_ions_, 0.);

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        // get local forces
        const int index = (*ion)->index();
        (*ion)->getForce(&forces[3 * index]);
        ion++;
    }

    // gather data to PE root
    const int size = 3 * num_ions_;
    mmpi.reduce(&forces[0], size, MPI_SUM, root);
}

void Ions::gatherLockedData(vector<int>& locked_data, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    locked_data.resize(num_ions_, 0);

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        // get local ion index
        const int index = (*ion)->index();
        if ((*ion)->locked()) locked_data[index] = 1;

        ++ion;
    }

    // gather data to PE root
    mmpi.reduce(&locked_data[0], num_ions_, MPI_SUM, root);
}

bool Ions::hasLockedAtoms() const
{
    short flag = 0;

    vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        if ((*ion)->locked())
        {
            flag = 1;
            break;
        }
        ++ion;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&flag, 1, MPI_MAX);

    return (flag == 1);
}

#if 0
void Ions::syncNames(const int nions, vector<string>& local_names, vector<string>& names)
{
    MGmol_MPI& mmpi ( *(MGmol_MPI::instance()) );
    MPI_Comm& comm=mmpi.comm();  
    // prepare to gather data
    int npes = mmpi.size();
    int num_ions = nions;
    vector<int>recvcnts(npes);
    vector<int>::iterator rcv;
    vector<int>disp(npes);   
    vector<int>::iterator disp_it; 
    int tot = 0;
    int totchars;
    // Gather data for atom names
    // first get length of each name
    vector<int>loc_name_len;
    vector<int>nameLen(num_ions,0);
    tot=0;
    for(vector<string>::iterator str=local_names.begin(); str!=local_names.end(); ++str){
       //tot+= *str.size();
       string s = *str;
       tot+= s.length();
       loc_name_len.push_back(s.length());
    }
    mmpi.allGatherV(loc_name_len, nameLen);

    vector<char>char_loc_names(tot);
    int idx = 0;
    for(vector<string>::iterator str=local_names.begin(); str!=local_names.end(); ++str){
       //tot+= *str.size();
       string s = *str;
       memcpy(&char_loc_names[idx], s.c_str(), s.size());
       idx+=s.size();
    }    
    assert(idx == tot);
    
    // gather data for atom names
    //get recv counts
    rcv = recvcnts.begin();
    MPI_Allgather(&tot, 1, MPI_INT, &(*rcv), 1, MPI_INT, comm);
    
    //get displacements
    totchars = recvcnts[npes-1];
    disp[0] = 0;
    for(int pos = 1; pos < npes; pos++){
       disp[pos] = disp[pos-1] + recvcnts[pos-1];
       totchars+=recvcnts[pos-1];
    }
    // gather atomic names
    vector<char> names_data(totchars);      
    vector<char>::iterator locnames_it=char_loc_names.begin();
    vector<char>::iterator names_it=names_data.begin();  
    disp_it=disp.begin();
    rcv = recvcnts.begin();  
    MPI_Allgatherv(&(*locnames_it), tot, MPI_CHAR, &(*names_it), &(*rcv), &(*disp_it), MPI_CHAR, comm);
      
//    vector<string> names;
    int pos = 0;
    for(int i=0; i<num_ions; i++)
    {
       char str[nameLen[i]+1];
       str[nameLen[i]] = '\0';
       memcpy(str, &names_data[pos], nameLen[i]*sizeof(char));
       string cstr;
       cstr.assign(str);
//       if(onpe0)puts(str);
       names.push_back(cstr);
       pos+=nameLen[i]*sizeof(char);       
    
    }
    return;
}
#endif

double Ions::getMaxNLradius() const
{
    double radius                       = 0;
    vector<Species>::const_iterator spi = species_.begin();
    while (spi != species_.end())
    {
        const double nlradius = spi->nlradius();
        radius                = radius > nlradius ? radius : nlradius;
        spi++;
    }
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&radius, 1, MPI_MAX);

    return radius;
}

double Ions::getMaxLradius() const
{
    double radius                       = 0;
    vector<Species>::const_iterator spi = species_.begin();
    while (spi != species_.end())
    {
        const double lradius = spi->lradius();
        radius               = radius > lradius ? radius : lradius;
        spi++;
    }
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&radius, 1, MPI_MAX);

    return radius;
}

void Ions::setupListIonsBoundaries(const double rmax)
{
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    const pb::Grid& mygrid   = mymesh->grid();

    double offset[3];
    for (short i = 0; i < 3; i++)
        offset[i] = (double)myPEenv.my_mpi(i) * div_lattice_[i];

    const double origin[3]
        = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };
    // const double end[3]    =
    // {origin[0]+lattice_[0],origin[1]+lattice_[1],origin[2]+lattice_[2]};

    // define extents for gathering ion data
    // define extents to be at least size of local domain and at most the size
    // of global domain
    double myrmax[3] = { rmax, rmax, rmax };
    for (short i = 0; i < 3; i++)
    {
        myrmax[i] = myrmax[i] >= div_lattice_[i] ? myrmax[i] : div_lattice_[i];
        myrmax[i] = myrmax[i] <= lattice_[i] ? myrmax[i] : lattice_[i];
    }
    // left / lower extent
    for (short i = 0; i < 3; i++)
    {
        list_boundary_left_[i] = origin[i] + offset[i] - myrmax[i];
    }

    // right/ upper extent
    for (short i = 0; i < 3; i++)
    {
        list_boundary_right_[i]
            = origin[i] + offset[i] + div_lattice_[i] + myrmax[i];
    }

    // setup and save data communication info
    /* get processor distribution */
    int nproc_xyz[3];
    nproc_xyz[0] = myPEenv.n_mpi_task(0);
    nproc_xyz[1] = myPEenv.n_mpi_task(1);
    nproc_xyz[2] = myPEenv.n_mpi_task(2);

    /* get domain info */
    double domain[3];
    domain[0] = mygrid.ll(0);
    domain[1] = mygrid.ll(1);
    domain[2] = mygrid.ll(2);

    /* compute processor width info */
    double proc_width[3];
    proc_width[0] = domain[0] / nproc_xyz[0];
    proc_width[1] = domain[1] / nproc_xyz[1];
    proc_width[2] = domain[2] / nproc_xyz[2];

    /* compute left and right steps in xyz directions */
    /* x-direction */
    lstep_[0] = min((int)(ceil(rmax / (proc_width[0]))), (nproc_xyz[0] - 1));
    rstep_[0] = min(lstep_[0], nproc_xyz[0] - lstep_[0] - 1);

    /* y-direction */
    lstep_[1] = min((int)(ceil(rmax / (proc_width[1]))), (nproc_xyz[1] - 1));
    rstep_[1] = min(lstep_[1], nproc_xyz[1] - lstep_[1] - 1);
    /* z-direction */
    lstep_[2] = min((int)(ceil(rmax / (proc_width[2]))), (nproc_xyz[2] - 1));
    rstep_[2] = min(lstep_[2], nproc_xyz[2] - lstep_[2] - 1);
}

// This function determines if an ion at position x,y,z is within range of
// the list of ions on this PE. This range is obtained from predetermined
// extents in list_boundary_left[x,y,z] and list_boundary_right[x,y,z].
bool Ions::inListIons(const double x, const double y, const double z)
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const double origin[3]
        = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };
    const double end[3] = { origin[0] + lattice_[0], origin[1] + lattice_[1],
        origin[2] + lattice_[2] };

    // define unique position of ion in periodic domain
    double t[3] = { x, y, z };
    for (short i = 0; i < 3; i++)
    {
        while (t[i] >= end[i])
            t[i] -= lattice_[i];
        while (t[i] < origin[i])
            t[i] += lattice_[i];
    }

    bool inList = false;

    // check to see if ion is in list
    if (((t[0] >= list_boundary_left_[0] && t[0] <= list_boundary_right_[0])
            || ((t[0] - lattice_[0]) >= list_boundary_left_[0]
                   && (t[0] - lattice_[0]) <= list_boundary_right_[0])
            || ((t[0] + lattice_[0]) >= list_boundary_left_[0]
                   && (t[0] + lattice_[0]) <= list_boundary_right_[0]))
        && ((t[1] >= list_boundary_left_[1] && t[1] <= list_boundary_right_[1])
               || ((t[1] - lattice_[1]) >= list_boundary_left_[1]
                      && (t[1] - lattice_[1]) <= list_boundary_right_[1])
               || ((t[1] + lattice_[1]) >= list_boundary_left_[1]
                      && (t[1] + lattice_[1]) <= list_boundary_right_[1]))
        && ((t[2] >= list_boundary_left_[2] && t[2] <= list_boundary_right_[2])
               || ((t[2] - lattice_[2]) >= list_boundary_left_[2]
                      && (t[2] - lattice_[2]) <= list_boundary_right_[2])
               || ((t[2] + lattice_[2]) >= list_boundary_left_[2]
                      && (t[2] + lattice_[2]) <= list_boundary_right_[2])))
        inList = true;

    return inList;
}

bool Ions::inLocalIons(const double x, const double y, const double z)
{
    // cout<<"inLocalIons..."<<endl;
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    const pb::Grid& mygrid   = mymesh->grid();

    double offset[3];
    for (short i = 0; i < 3; i++)
        offset[i] = (double)myPEenv.my_mpi(i) * div_lattice_[i];

    const double origin[3]
        = { mygrid.origin(0), mygrid.origin(1), mygrid.origin(2) };
    const double end[3] = { origin[0] + lattice_[0], origin[1] + lattice_[1],
        origin[2] + lattice_[2] };

    double t[3] = { x, y, z };
    for (short i = 0; i < 3; i++)
    {
        while (t[i] >= end[i])
            t[i] -= lattice_[i];
        while (t[i] < origin[i])
            t[i] += lattice_[i];
        t[i] -= origin[i];
        t[i] -= offset[i];
    }

    bool inList = false;

    if ((t[0] >= 0. && t[0] < (div_lattice_[0]))
        && (t[1] >= 0. && t[1] < (div_lattice_[1]))
        && (t[2] >= 0. && t[2] < (div_lattice_[2])))
    {
        inList = true;
    }

    return inList;
}

// void Ions::checkUnicityLocalIons()
//{
//    //cout<<"Ions::checkUnicityLocalIons()..."<<endl;
//    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
//    //cout<<"size list = "<<list_ions_.size()<<endl;
//    for(auto& ion : list_ions_)
//    {
//        int here = ion->here();
//        mmpi.allreduce(&here, 1, MPI_SUM);
//        if(here != 1)
//        {
//            cout<<"Ion "<<ion->name()<<" is here on multiple tasks"<<endl;
//        }
//    }
//}

int Ions::getNumIons(void)
{
    if (num_ions_ < 0) computeNumIons();

    return num_ions_;
}

void Ions::computeNumIons(void)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    num_ions_ = (int)local_ions_.size();
    mmpi.allreduce(&num_ions_, 1, MPI_SUM);
}

double Ions::getMaxListRadius() const
{
    // get radius of projectors
    const double nlradius = getMaxNLradius();
    const double lradius  = getMaxLradius();

    double rmax = nlradius > lradius ? nlradius : lradius;

    Control& ct(*(Control::instance()));
    rmax
        = rmax > ct.maxDistanceAtomicInfo() ? rmax : ct.maxDistanceAtomicInfo();

    return rmax;
}

// augment Ions Data
void Ions::augmentIonsData(const int nsteps, const int dir, const int disp,
    const int locSize, const int maxLocSize, vector<IonData>& iondata,
    int* offset)
{
    // setup buffers for data transfer
    int csize              = IonData_MaxStrLength * locSize;
    const int cbuff_size_r = IonData_MaxStrLength * maxLocSize;
    char* cbuff            = new char[2 * cbuff_size_r];
    char* cbuff_r          = &cbuff[cbuff_size_r];

    // 7 int/atom in IonData
    int isize              = 7 * locSize + 1;
    const int ibuff_size_r = 7 * maxLocSize + 1;
    int* ibuff             = new int[2 * ibuff_size_r];
    int* ibuff_r           = &ibuff[ibuff_size_r];

    // 16 double/atom in IonData
    int dsize              = 16 * locSize;
    const int dbuff_size_r = 16 * maxLocSize;
    double* dbuff          = new double[2 * dbuff_size_r];
    double* dbuff_r        = &dbuff[dbuff_size_r];

    // pack data into buffers
    IonData::packIonData(cbuff, ibuff, dbuff, iondata);

    // transfer data in multiple directions
    MPI_Request request1[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Request request2[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Request request3[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };

    int src;
    int dst;
    if (disp == -1)
    {
        src = source_[dir];
        dst = dest_[dir];
    }
    else
    {
        src = dest_[dir];
        dst = source_[dir];
    }

    int step = 0;
    while (step < nsteps)
    {
        // chars data
        MPI_Irecv(
            cbuff_r, cbuff_size_r, MPI_CHAR, src, 0, cart_comm_, &request1[0]);
        MPI_Isend(cbuff, csize, MPI_CHAR, dst, 0, cart_comm_, &request1[1]);
        // integers data
        MPI_Irecv(
            ibuff_r, ibuff_size_r, MPI_INT, src, 0, cart_comm_, &request2[0]);
        MPI_Isend(ibuff, isize, MPI_INT, dst, 0, cart_comm_, &request2[1]);
        // doubles data
        MPI_Irecv(dbuff_r, dbuff_size_r, MPI_DOUBLE, src, 0, cart_comm_,
            &request3[0]);
        MPI_Isend(dbuff, dsize, MPI_DOUBLE, dst, 0, cart_comm_, &request3[1]);

        // update lists with most recent data
        vector<IonData>::iterator idata = ions_data_.begin() + *offset;
        while (idata != ions_data_.end())
        {
            addIon2Lists(*idata);

            idata++;
        }
        *offset = (int)ions_data_.size();

        // wait to complete data transfer
        MPI_Waitall(2, request1, MPI_STATUS_IGNORE);
        MPI_Waitall(2, request2, MPI_STATUS_IGNORE);
        MPI_Waitall(2, request3, MPI_STATUS_IGNORE);

        // unpack: update ions_data for next call to this function
        char* cptr   = &cbuff_r[0];
        int* iptr    = &ibuff_r[0];
        double* dptr = &dbuff_r[0];
        // first get data size
        const int rsize = *(iptr++);
        // unpack
        for (int i = 0; i < rsize; i++)
        {
            IonData newdata;

            newdata.unpack(cptr, iptr, dptr);

            // update ions_data if ion is in list_ions range
            if (inListIons(newdata.current_position[0],
                    newdata.current_position[1], newdata.current_position[2]))
            {
                ions_data_.push_back(newdata);
            }
        }
        // done with this phase of data transfer

        // update sizes and copy recv buffer to send buffer
        if (step < nsteps - 1)
        {
            csize = rsize * IonData_MaxStrLength;
            isize = 7 * rsize + 1;
            dsize = 16 * rsize;
            memcpy(cbuff, cbuff_r, csize * sizeof(char));
            memcpy(ibuff, ibuff_r, isize * sizeof(int));
            memcpy(dbuff, dbuff_r, dsize * sizeof(double));
        }
        step++;
    }

    delete[] cbuff;
    delete[] ibuff;
    delete[] dbuff;
}

void Ions::updateForcesInteractingIons()
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    // get computed forces into fion_
    getForces(fion_);

    // initialize with local names and forces
    DistributedIonicData forces_data(local_names_, fion_);

    for (short dir = 0; dir < 3; dir++)
    {
        // initial data
        const int lsize = forces_data.size();
        int maxsize     = lsize;
        mmpi.allreduce(&maxsize, 1, MPI_MAX);

        // send local data
        DistributedIonicData data2send(forces_data);

        // send right to left
        int disp = -1;
        forces_data.augmentData(
            lstep_[dir], dir, disp, lsize, maxsize, data2send);

        // send left to right
        disp = 1;
        forces_data.augmentData(
            rstep_[dir], dir, disp, lsize, maxsize, data2send);
    }

    int ia = 0;
    for (vector<string>::const_iterator it = interacting_names_.begin();
         it != interacting_names_.end(); ++it)
    {
        double* force = interacting_fion_[ia];
        forces_data.getData(*it, force);
        ++ia;
    }
}

void Ions::updateTaupInteractingIons()
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    // initialize with local names and forces
    DistributedIonicData taup_data(local_names_, taup_);

    for (short dir = 0; dir < 3; dir++)
    {
        // initial data
        const int lsize = taup_data.size();
        int maxsize     = lsize;
        mmpi.allreduce(&maxsize, 1, MPI_MAX);

        DistributedIonicData data2send(taup_data);

        // send right to left
        int disp = -1;
        taup_data.augmentData(
            lstep_[dir], dir, disp, lsize, maxsize, data2send);

        // send left to right
        disp = 1;
        taup_data.augmentData(
            rstep_[dir], dir, disp, lsize, maxsize, data2send);
    }

    int ia = 0;
    for (vector<string>::const_iterator it = interacting_names_.begin();
         it != interacting_names_.end(); ++it)
    {
        double* taup = interacting_taup_[ia];
        taup_data.getData(*it, taup);
        ++ia;
    }
}

// update list of local ions
void Ions::updateListIons()
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    // collect local_ions data
    // assume local_ions size is same as size of ions_names vector
    ions_data_.clear();

    {
        IonData idata;
        for (auto& lion : local_ions_)
        {
            lion->getIonData(idata);

            // populate ions_data_ list
            ions_data_.push_back(idata);
        }
    }

    // update list ions from current ion positions
    // Note: this is based on data from MD vectors

    // First cleanup list_ions_
    local_ions_.clear();
    // delete current ions from list
    vector<Ion*>::iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        delete *ion;
        ion++;
    }
    list_ions_.clear();

    // Update list starting with local ion data.
    // This enables overlapping data accumulation with communication.
    // Most recent/ final data is accumulated later
    int offset = 0;
    for (short dir = 0; dir < 3; dir++)
    {
        // initial data
        const int lsize = (int)ions_data_.size();
        int maxsize     = lsize;
        mmpi.allreduce(&maxsize, 1, MPI_MAX);
        // if(onpe0)cout<<"lsize = "<<ions_data_.size()<<" list size =
        // "<<list_ions_.size()<<endl;
        vector<IonData> data2send(ions_data_);

        // send right to left
        int disp = -1;
        augmentIonsData(
            lstep_[dir], dir, disp, lsize, maxsize, data2send, &offset);
        // if(onpe0)cout<<"lsize2 = "<<ions_data_.size()<<" list size2 =
        // "<<list_ions_.size()<<endl;
        // send left to right
        disp = 1;
        augmentIonsData(
            rstep_[dir], dir, disp, lsize, maxsize, data2send, &offset);
        // if(onpe0)cout<<"lsize3 = "<<ions_data_.size()<<" list size3 =
        // "<<list_ions_.size()<<endl;
    }

    // perform last data accumulation here
    vector<IonData>::const_iterator idata = ions_data_.begin() + offset;
    while (idata != ions_data_.end())
    {
        addIon2Lists(*idata);

        idata++;
    }

    // clear ions_data_
    ions_data_.clear();
}

void Ions::addIon2Lists(const IonData& data)
{
    if (inListIons(data.current_position[0], data.current_position[1],
            data.current_position[2]))
    {
        Ion* newion = new Ion(getSpecies(data.atomic_num), data);

        // update list_ions_
        list_ions_.push_back(newion);

        // update local_ions_ list
        if (inLocalIons(data.current_position[0], data.current_position[1],
                data.current_position[2]))
        {
            newion->set_here(true);

            newion->setFromIonData(data);

            local_ions_.push_back(newion);
        }
        else
        {
            newion->set_here(false);
        }
    }
}

void Ions::clearStepperData()
{
    local_names_.clear();
    tau0_.clear();
    atmove_.clear();
    taum_.clear();
    taup_.clear();
    fion_.clear();
    velocity_.clear();
    rand_states_.clear();
    pmass_.clear();
    gids_.clear();
}

void Ions::initStepperData()
{
    clearStepperData();

    vector<Ion*>::iterator lion = local_ions_.begin();
    while (lion != local_ions_.end())
    {
        local_names_.push_back((*lion)->name());
        atmove_.push_back(!(*lion)->locked());
        pmass_.push_back((*lion)->getMass());
        gids_.push_back((*lion)->index());

        for (short i = 0; i < 3; i++)
        {
            taum_.push_back((*lion)->old_position(i));
            tau0_.push_back((*lion)->position(i));
            fion_.push_back((*lion)->force(i));
            velocity_.push_back((*lion)->velocity(i));
            rand_states_.push_back((*lion)->randomState(i));
        }

        lion++;
    }
    // initialize taup to enable computing velocities
    int size_tau = (int)tau0_.size();
    taup_        = tau0_;

    int ione     = 1;
    double alpha = 1.;
    DAXPY(&size_tau, &alpha, &tau0_[0], &ione, &taup_[0], &ione);
    alpha = -1.;
    DAXPY(&size_tau, &alpha, &taum_[0], &ione, &taup_[0], &ione);
}

// update local ion data with ionic stepper data
// Note: forces have already been updated.
// fion may include Langevin forces (so ignore??)
void Ions::updateIons()
{
    assert(tau0_.size() == 3 * local_ions_.size());

    Control& ct(*(Control::instance()));

    // update local_ions data
    vector<Ion*>::iterator ion = local_ions_.begin();
    int ia                     = 0;
    while (ion != local_ions_.end())
    {
        (*ion)->setPosition(
            tau0_[3 * ia + 0], tau0_[3 * ia + 1], tau0_[3 * ia + 2]);

        //        (*ion)->setForce(fion_[3*ia+0],
        //                          fion_[3*ia+1],
        //                          fion_[3*ia+2]);

        (*ion)->setRandomState(rand_states_[3 * ia + 0],
            rand_states_[3 * ia + 1], rand_states_[3 * ia + 2]);

        ion++;
        ia++;
    }

    // compute and set velocities if needed
    if (ct.dt > 0.) setVelocities(taum_, tau0_, ct.dt);

    // update various list of ions
    updateListIons();

    setup_ = false;
}

void Ions::shiftIons(const Vector3D& shift)
{
    // update local_ions data
    vector<Ion*>::iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        (*ion)->shiftPositionXLBOMDTest(shift);
        ion++;
    }

    // update various list of ions
    updateListIons();

    setup_ = false;
}

void Ions::rescaleVelocities(const double factor)
{
    if (fabs(factor - 1.) < 1.e-12) return;

    if (onpe0)
    {
        cout << "Ions::rescaleVelocities() with factor " << factor << endl;
    }
    vector<Ion*>::iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        (*ion)->rescaleVelocity(factor);

        ion++;
    }
}

double Ions::computeMinLocalSpacing() const
{
    Control& ct     = *(Control::instance());
    double distance = 1.e6;

    vector<Ion*>::const_iterator ion1 = local_ions_.begin();
    vector<Ion*>::const_iterator ion2 = local_ions_.begin();
    while (ion1 != local_ions_.end())
    {
        while (ion2 != local_ions_.end())
        {
            if (ion1 != ion2)
            {
                double d = (*ion1)->minimage(**ion2, lattice_, ct.bc);
                if (d < distance) distance = d;
            }
            ion2++;
        }
        ion1++;
    }

    return distance;
}
