// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Ions.h"
#include "Control.h"
#include "HDFrestart.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "Species.h"
#include "Timer.h"
#include "hdf_tools.h"
#include "mgmol_mpi_tools.h"
#include "tools.h"

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <iterator>
#include <list>
#include <map>

Timer ions_setupInteractingIons_tm("ions_setupInteractingIons");
Timer ions_setup_tm("ions::setup");

const double ang2bohr = 1.8897269;
// const double rmax = 8.0;

std::map<std::string, short> Ions::map_species_
    = { { "H", 1 }, { "D", 1 }, { "Li", 3 }, { "Be", 4 }, { "B", 5 },
          { "C", 6 }, { "N", 7 }, { "O", 8 }, { "F", 9 }, { "Na", 11 },
          { "Mg", 12 }, { "Al", 13 }, { "Si", 14 }, { "P", 15 }, { "S", 16 },
          { "Cl", 17 }, { "K", 19 }, { "Ca", 20 }, { "Cr", 24 }, { "Mn", 25 },
          { "Fe", 26 }, { "Co", 27 }, { "Ni", 28 }, { "Cu", 29 }, { "Zn", 30 },
          { "Ga", 31 }, { "Ge", 32 }, { "La", 57 }, { "Au", 79 } };

int Ions::num_ions_          = -1;
short Ions::max_num_proj_    = -1;
double Ions::max_Vl_radius_  = -1.;
double Ions::max_Vnl_radius_ = -1.;

Ions::Ions(const double lat[3], const std::vector<Species>& sp) : species_(sp)
{
    for (short i = 0; i < 3; i++)
    {
        assert(lat[i] > 0.);
        if (lat[i] > 10000.)
        {
            (*MPIdata::serr) << "Ions constructor: lattice[" << i
                             << "]=" << lat[i] << "!!!" << std::endl;
            exit(2);
        }
        lattice_[i] = lat[i];
    }
    setup_            = false;
    has_locked_atoms_ = false;

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
    std::vector<Ion*>::const_iterator ion = ions.list_ions_.begin();
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

    std::vector<Ion*>::iterator iion       = list_ions_.begin();
    std::vector<Ion*>::const_iterator cion = ions.list_ions_.begin();
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
                         << " initialized" << std::endl;
#endif

    std::vector<Ion*>::iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        max_num_proj_ = std::max(max_num_proj_, (*ion)->nProjectors());

        ion++;
    }

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&max_num_proj_, 1, MPI_MAX);
}

void Ions::setup()
{
    Control& ct = *(Control::instance());

    if (ct.verbose > 0) printWithTimeStamp("Ions::setup()...", std::cout);

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
        printWithTimeStamp("Ions::setup()... individual ions...", std::cout);

    std::vector<Ion*>::iterator ion = list_ions_.begin();
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
    std::vector<Ion*>::iterator ion = list_ions_.begin();
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
        printWithTimeStamp("Ions::setupListOverlappingIons()...", std::cout);

    overlappingNL_ions_.clear();
    overlappingVL_ions_.clear();

    for (auto ion : list_ions_)
    {
        if (ion->map_nl())
        {
            overlappingNL_ions_.push_back(ion);
        }
    }

    for (auto ion : list_ions_)
    {
        if (ion->map_l())
        {
            overlappingVL_ions_.push_back(ion);
        }
    }
}

void Ions::setupInteractingIons()
{
    Control& ct = *(Control::instance());
    if (ct.verbose > 0)
        printWithTimeStamp("Ions::setupInteractingIons()...", std::cout);

    ions_setupInteractingIons_tm.start();
    const double rmax = ct.maxDistanceAtomicInfo();
    if (onpe0 && ct.verbose > 1)
    {
        (*MPIdata::sout) << "Ions::setupInteractingIons() with radius " << rmax
                         << std::endl;
    }

    interacting_ions_.clear();

    for (auto ion1 : list_ions_)
    {
        // is ion1 interacting with any local ions?
        for (auto ion2 : local_ions_)
        {
            const double r12 = ion1->minimage(*ion2, lattice_, ct.bcPoisson);

            if (r12 < rmax)
            {
                // ion1 is interacting with local ions
                interacting_ions_.push_back(ion1);
                break;
            }
        }
    }

    //(*MPIdata::sout)<<"Number of interacting ions =
    //"<<interacting_ions_.size()<<endl;
    ions_setupInteractingIons_tm.stop();
}

// setup arrays to be used in constraints enforcement
// using references to local_ions and extra "dummy" data
void Ions::setupContraintsData(std::vector<Ion*>& ions_for_constraints)
{
    Control& ct = *(Control::instance());

    if (ct.verbose > 0)
        printWithTimeStamp("Ions::setupContraintsData()...", std::cout);
    const int nnloc = ions_for_constraints.size() - local_ions_.size();
    // std::cout<<"interacting_ions_.size()="<<interacting_ions_.size()<<endl;
    // std::cout<<"local_ions_.size()="<<local_ions_.size()<<endl;
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

    int ia                                = 0; // count local ions
    int ib                                = 0; // count non-local ions
    std::vector<Ion*>::const_iterator ion = ions_for_constraints.begin();
    while (ion != ions_for_constraints.end())
    {
        if (isLocal((*ion)->name()))
        {
            assert(3 * ia < static_cast<int>(tau0_.size()));
            assert(3 * ia < static_cast<int>(taup_.size()));
            assert(3 * ia < static_cast<int>(fion_.size()));
            assert(ia < static_cast<int>(pmass_.size()));
            assert(ia < static_cast<int>(atmove_.size()));

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
            assert(3 * ib < static_cast<int>(tau0_dummy_.size()));
            assert(3 * ib < static_cast<int>(taup_dummy_.size()));
            assert(3 * ib < static_cast<int>(fion_dummy_.size()));

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

            assert(ib < static_cast<int>(names_dummy_.size()));
            assert(ib < static_cast<int>(pmass_dummy_.size()));
            assert(ib < static_cast<int>(atmove_dummy_.size()));
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
    std::vector<double> forces(3 * nlions, 0.);

    std::vector<Ion*>::const_iterator ion1 = local_ions_.begin();
    int ion1_index                         = 0;
    ;
    while (ion1 != local_ions_.end())
    {
        const double z1 = (*ion1)->getZion();
        assert(z1 >= 0.);

        const double rc1 = (*ion1)->getRC();

        std::vector<Ion*>::const_iterator ion2 = interacting_ions_.begin();
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

    std::vector<Ion*>::iterator lion = local_ions_.begin();
    int ion_index                    = 0;
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

    for (auto ion : local_ions_)
    {
        energy += ion->eself();
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

    std::vector<Ion*>::const_iterator ion1 = interacting_ions_.begin();
    while (ion1 != interacting_ions_.end())
    {
        if ((*ion1)->here())
        {
            std::vector<Ion*>::const_iterator ion2 = interacting_ions_.begin();
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

    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double tenergy           = 0.;
    MPI_Allreduce(&energy, &tenergy, 1, MPI_DOUBLE, MPI_SUM, myPEenv.comm());

    // take half the result to account for double counting loop over
    // interacting_ions
    energy = 0.5 * tenergy;

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
//                     << offset[2] << ")" << std::endl;
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
//    std::vector<Ion*>::iterator ion = list_ions_.begin();
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
//                         << std::endl;
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
//                             << myPEenv.mytask() << std::endl;
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
//        (*MPIdata::sout) << " num_ions=" << num_ions_ << std::endl;
//        (*MPIdata::sout) << " isum=" << isum << std::endl;
//    }
//    if (isum != num_ions_)
//    {
//        (*MPIdata::sout) << " num_ions != isum !!!!" << std::endl;
//        ct.global_exit(2);
//    }
//
//#ifdef DEBUG
//    mmpi.allreduce(&isum2, 1, MPI_SUM);
//    (*MPIdata::sout) << " Isum2=" << isum2 << std::endl;
//    assert(isum2 == num_ions_);
//#endif
//
//    // if( onpe0 )
//    //    (*MPIdata::sout)<<"Ions::associate2PE() done..."<<endl;
//}

// Writes out the positions of all the ions
// PE root does the writing
void Ions::printPositionsGlobal(std::ostream& os, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    std::map<int, std::string> ion_names;
    std::vector<double> positions;
    std::vector<int> islocked;

    gatherNames(ion_names, root, mmpi.commSameSpin());
    gatherPositions(positions, root);
    gatherLockedData(islocked, root);

    if (mmpi.mypeGlobal() == root)
    {

        os << std::endl << " IONIC POSITIONS:" << std::endl;

        os << std::setw(8) << "Atoms" << std::setw(8) << "X" << std::setw(10)
           << "Y" << std::setw(10) << "Z" << std::endl;

        os.setf(std::ios::right, std::ios::adjustfield);
        os.setf(std::ios::fixed, std::ios::floatfield);

        int ion_index                                     = 0;
        std::map<int, std::string>::const_iterator ion_id = ion_names.begin();
        while (ion_id != ion_names.end())
        {
            const int pos = 3 * ion_index;

            os << "$$ ";
            if (islocked[ion_index])
                os << "*";
            else
                os << " ";
            os << std::setw(4) << ion_id->second << std::setw(10)
               << std::setprecision(4) << std::fixed << positions[pos]
               << std::setw(10) << positions[pos + 1] << std::setw(10)
               << positions[pos + 2] << std::endl;

            ion_index++;
            ion_id++;
        }

        os << std::endl;
    }
}

// Writes out the postions of the local ions and their displacements from their
// initial postions.
void Ions::printPositionsLocal(std::ostream& os, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    if (mmpi.mypeGlobal() == root)
    {

        os << std::endl
           << " IONIC POSITIONS AND DISPLACEMENTS ON PE" << root << ":"
           << std::endl;

        os << std::setw(8) << "Atoms" << std::setw(8) << "X" << std::setw(10)
           << "Y" << std::setw(10) << "Z" << std::setw(10) << "dX"
           << std::setw(10) << "dY" << std::setw(10) << "dZ" << std::endl;

        os.setf(std::ios::right, std::ios::adjustfield);
        os.setf(std::ios::fixed, std::ios::floatfield);

        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            (*ion)->printPosition(os);
            ion++;
        }

        os << std::endl;
    }
}
// Writes out the postions of the ions and their displacements from their
// initial postions.
void Ions::printPositions(std::ostream& os, const int root) const
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
        (*MPIdata::sout) << "Ions::writeAtomicNumbers()..." << std::endl;
    }
    std::vector<int> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherAtomicNumbers(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
        // for(std::vector<int>::iterator it =data.begin();
        //                          it!=data.end();
        //                        ++it)
        //    std::cout<<"Number="<<*it<<endl;

        // std::cout<<"ms="<<ms<<endl;

        // for(short i=0;i<ms;i++)data.push_back(-100*(mytask+1));
        // std::vector<int> old(data);
        // data.clear();
        // for(short i=0;i<ms;i++)
        //{
        //    data.push_back(old[i]);
        //    data.push_back(-1);
        //}

        size_t dims[2] = { data.size(), 1 };

        std::string datasetname("/Atomic_numbers");
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
        (*MPIdata::sout) << "Ions::writeAtomNames" << std::endl;
    }

    // gather data to print locally into std::vector "data"
    std::vector<std::string> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherNames(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
            std::string empty_string;
            for (short i = s; i < ms; i++)
                data.push_back(empty_string);
        }

        size_t dims[2] = { data.size(), 1 };

        std::string datasetname("/Atomic_names");
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

void Ions::lockAtom(const std::string& name)
{
    std::vector<Ion*>::iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        std::string name_ion((*ion)->name());
        if (name.compare(name_ion) == 0)
        {
            (*ion)->lock();
            if (onpe0) (*MPIdata::sout) << "Lock atom " << name << std::endl;
            break;
        }
        ion++;
    }
}

void Ions::readLockedAtomNames(HDFrestart& h5f_file)
{
    int dim = 0;

    if (dim == 0) return;

    std::vector<std::string> data;
    h5f_file.readLockedAtomNames(data);

    for (std::vector<std::string>::const_iterator i   = data.begin(),
                                                  end = data.end();
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
        (*MPIdata::sout) << "Ions::writeLockedAtomsNames" << std::endl;
    }

    // gather data to print locally into std::vector "data"
    std::vector<std::string> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherLockedNames(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
            std::string empty_string;
            for (short i = s; i < ms; i++)
                data.push_back(empty_string);
            if (ms == 0) return;
        }

        size_t dims[2] = { data.size(), 1 };

        std::string datasetname("/LockedAtomsNames");
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
        (*MPIdata::sout) << "Ions::writeAtomicIDs()..." << std::endl;
    }

    // gather data to print locally into std::vector "data"
    std::vector<int> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherIndexes(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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

        std::string datasetname("/Atomic_IDs");
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
        (*MPIdata::sout) << "Ions::writeAtomicNLprojIDs()..." << std::endl;
    }

    // gather data to print locally into std::vector "data"
    std::vector<int> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherNLprojIds(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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

        std::string datasetname("/AtomicNLproj_IDs");
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
        (*MPIdata::sout) << "Ions::writePositions" << std::endl;
    }

    std::vector<double> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherPositions(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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

        std::string datasetname("/Ionic_positions");
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
        (*MPIdata::sout) << "Ions::setFromRestartFile()..." << std::endl;
    }

    // set up list boundaries
    // get radius of projectors
    double rmax = getMaxListRadius();
    assert(rmax > 0.);
    setupListIonsBoundaries(rmax);

    std::vector<int> at_numbers;
    h5_file.readAtomicNumbers(at_numbers);
    std::vector<int> at_indexes;
    int nidxs = h5_file.readAtomicIDs(at_indexes);
    std::vector<int> at_nlprojIds;
    int npids = h5_file.readAtomicNLprojIDs(at_nlprojIds);
    std::vector<std::string> at_names;
    h5_file.readAtomicNames(at_names);
    if (onpe0 && ct.verbose > 2)
    {
        std::cout << "HDF file: at nb=" << at_numbers.size() << std::endl;
        std::cout << "HDF file: names nb=" << at_names.size() << std::endl;
        std::cout << "HDF file: indexes nb=" << at_indexes.size() << std::endl;
        std::cout << "HDF file: at_nlprojIds nb=" << at_nlprojIds.size()
                  << std::endl;
    }

    // if reading "old" format with replicated atoms, fill up empty arrays
    // (atomic indexes were not saved before)
    if (nidxs == -1 && at_numbers.size() > 0)
    {
        for (unsigned int i = 0; i < at_numbers.size(); i++)
            at_indexes.push_back(i);
    }
    if (npids == -1 && at_numbers.size() > 0)
    {
        for (unsigned int i = 0; i < at_numbers.size(); i++)
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
                         << " names... on PE0" << std::endl;
    }

    assert(at_numbers.size() == at_names.size());

    double default_vel[3]    = { 0., 0., 0. };
    double default_coords[3] = { 0., 0., 0. };

    for (unsigned int i = 0; i < at_numbers.size(); i++)
    {
        const int atnum                          = at_numbers[i];
        std::vector<Species>::const_iterator spi = species_.begin();
        while (spi != species_.end())
        {
            if (atnum == spi->getAtomicNumber()) break;
            spi++;
        }
        if (onpe0 && ct.verbose > 3)
            (*MPIdata::sout)
                << "New Ion with name " << at_names[i] << std::endl;
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
        (*MPIdata::sout) << "Read ionic positions from hdf5 file" << std::endl;

    std::vector<double> data;
    h5_file.readAtomicPositions(data);

    int i                           = 0;
    std::vector<Ion*>::iterator ion = local_ions_.begin();
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
        (*MPIdata::sout) << "Ions::writeVelocities" << std::endl;
    }

    std::vector<double> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherVelocities(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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

        std::string datasetname("/Ionic_velocities");
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
        (*MPIdata::sout) << "Ions::writeRandomStates()..." << std::endl;
    }

    std::vector<unsigned short> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherRandStates(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
                << std::endl;
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

        std::string datasetname("/Ionic_RandomStates");
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
        (*MPIdata::sout) << "Remove mass center motion" << std::endl;

    std::vector<double> mass(local_ions_.size());
    std::vector<double> velocities(3 * local_ions_.size());

    std::vector<Ion*>::iterator ion = local_ions_.begin();
    int i                           = 0;
    double tmp[4]                   = { 0., 0., 0., 0. };
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
        (*MPIdata::sout) << std::setprecision(12);
        if (onpe0) (*MPIdata::sout) << "V[" << j << "]=" << tmp[j] << std::endl;
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
        if (onpe0) (*MPIdata::sout) << "V[" << j << "]=" << mv[j] << std::endl;
    }
#endif
}

void Ions::readRestartVelocities(HDFrestart& h5_file)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Read atomic velocities from hdf5 file"
                         << std::endl;

    std::vector<double> data;
    h5_file.readAtomicVelocities(data);

    std::vector<Ion*>::iterator ion = local_ions_.begin();
    int i                           = 0;
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
        (*MPIdata::sout) << "Read atomic RandomStates from hdf5 file"
                         << std::endl;

    std::vector<unsigned short> data;
    h5f_file.readRestartRandomStates(data);

    std::vector<Ion*>::iterator ion = local_ions_.begin();
    int i                           = 0;
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
        (*MPIdata::sout) << "Write ionic forces in hdf5 file" << std::endl;

    std::vector<double> data;
    if (h5f_file.gatherDataX())
    {
        Mesh* mymesh             = Mesh::instance();
        const pb::PEenv& myPEenv = mymesh->peenv();
        MPI_Comm comm            = myPEenv.comm_x();

        gatherForces(data, 0, comm);
    }
    else
    {
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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

        std::string datasetname("/Ionic_forces");
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
void Ions::printForcesGlobal(std::ostream& os, const int root) const
{
    double maxf = 0., avfx = 0., avfy = 0., avfz = 0., maxfx = 0., maxfy = 0.,
           maxfz  = 0.;
    int num_atoms = 0;

    double sum_forces[3] = { 0., 0., 0. };
    int num_movable      = 0;

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    std::map<int, std::string> ion_names;
    std::vector<double> positions;
    std::vector<double> forces;
    std::vector<int> islocked;

    gatherNames(ion_names, root, mmpi.commSpin());
    gatherPositions(positions, root);
    gatherForces(forces, root);
    gatherLockedData(islocked, root);

    if (mmpi.mypeGlobal() == root)
    {
        os << "IONIC POSITIONS AND FORCES:" << std::endl;

        os << std::setiosflags(std::ios::left) << std::setw(8) << "Atoms"
           << resetiosflags(std::ios::left) << std::setw(8) << "X"
           << std::setw(10) << "Y" << std::setw(10) << "Z" << std::setw(10)
           << "FX" << std::setw(10) << "FY" << std::setw(10) << "FZ"
           << std::endl;

        int ion_index                                     = 0;
        std::map<int, std::string>::const_iterator ion_id = ion_names.begin();
        while (ion_id != ion_names.end())
        {

            const int pos = 3 * ion_index;
            os << "## ";
            if (islocked[ion_index])
                os << "*";
            else
                os << " ";
            os << std::setw(4) << ion_id->second
               << std::setiosflags(std::ios::right) << std::setw(10)
               << std::setprecision(4) << std::fixed << positions[pos]
               << std::setw(10) << positions[pos + 1] << std::setw(10)
               << positions[pos + 2] << std::setprecision(7) << std::scientific
               << std::setw(16) << forces[pos] << std::setw(16)
               << forces[pos + 1] << std::setw(16) << forces[pos + 2]
               << std::endl;

            if (!islocked[ion_index])
            {

                avfx += fabs(forces[pos]);
                avfy += fabs(forces[pos + 1]);
                avfz += fabs(forces[pos + 2]);

                double ff = forces[pos] * forces[pos]
                            + forces[pos + 1] * forces[pos + 1]
                            + forces[pos + 2] * forces[pos + 2];

                maxf = std::max(maxf, ff);

                maxfx = std::max(maxfx, fabs(forces[pos]));
                maxfy = std::max(maxfy, fabs(forces[pos + 1]));
                maxfz = std::max(maxfz, fabs(forces[pos + 2]));

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
        os << std::endl << "Global Statistics:" << std::endl;
        os << "==========================" << std::endl;

        const double inv_num_atoms = 1. / (double)num_atoms;
        for (short ii = 0; ii < 3; ii++)
            sum_forces[ii] *= inv_num_atoms;

        if (num_movable)
        {
            double inv_num_movable = 1. / (double)num_movable;
            avfx                   = avfx * inv_num_movable;
            avfy                   = avfy * inv_num_movable;
            avfz                   = avfz * inv_num_movable;

            os << std::scientific;
            os << " mean F on movable ions  = (" << avfx << "," << avfy << ","
               << avfz << ")" << std::endl;
            os << " max F on movable ions   = (" << maxfx << "," << maxfy << ","
               << maxfz << ")" << std::endl;
            os << " max |F| on movable ions = " << sqrt(maxf) << std::endl;
        }
        os << " Sum forces on all ions  = (" << sum_forces[0] << ","
           << sum_forces[1] << "," << sum_forces[2] << ")" << std::endl;
        os << std::fixed;
    }
}

// Writes out the postions of the local ions and the current forces on them
void Ions::printForcesLocal(std::ostream& os, const int root) const
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
        os << std::endl
           << std::endl
           << "IONIC POSITIONS AND FORCES ON CENTERED ON PE" << root << ":"
           << std::endl;

        os << std::setiosflags(std::ios::left) << std::setw(8) << "Atoms"
           << std::resetiosflags(std::ios::left) << std::setw(8) << "X"
           << std::setw(10) << "Y" << std::setw(10) << "Z" << std::setw(10)
           << "FX" << std::setw(10) << "FY" << std::setw(10) << "FZ"
           << std::endl;

        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {

            (*ion)->printPositionAndForce(os);

            if (!(*ion)->locked())
            {

                avg_forces[0] += fabs((*ion)->force(0));
                avg_forces[1] += fabs((*ion)->force(1));
                avg_forces[2] += fabs((*ion)->force(2));

                double ff = (*ion)->norm2F();
                maxf[0]   = std::max(maxf[0], ff);

                max_forces[0] = std::max(max_forces[0], fabs((*ion)->force(0)));
                max_forces[1] = std::max(max_forces[1], fabs((*ion)->force(1)));
                max_forces[2] = std::max(max_forces[2], fabs((*ion)->force(2)));

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
        std::vector<Ion*>::const_iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {

            if (!(*ion)->locked())
            {

                avg_forces[0] += fabs((*ion)->force(0));
                avg_forces[1] += fabs((*ion)->force(1));
                avg_forces[2] += fabs((*ion)->force(2));

                double ff = (*ion)->norm2F();
                maxf[0]   = std::max(maxf[0], ff);

                max_forces[0] = std::max(max_forces[0], fabs((*ion)->force(0)));
                max_forces[1] = std::max(max_forces[1], fabs((*ion)->force(1)));
                max_forces[2] = std::max(max_forces[2], fabs((*ion)->force(2)));

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
        os << std::endl << "Global Statistics:" << std::endl;
        os << "==========================" << std::endl;
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

            os << std::scientific;
            os << " mean F on movable ions  = (" << avg_forces[0] << ","
               << avg_forces[1] << "," << avg_forces[2] << ")" << std::endl;
            os << " max F on movable ions   = (" << max_forces[0] << ","
               << max_forces[1] << "," << max_forces[2] << ")" << std::endl;
            os << " max |F| on movable ions = " << sqrt(maxf[0]) << std::endl;
        }
        os << " Sum forces on all ions  = (" << sum_forces[0] << ","
           << sum_forces[1] << "," << sum_forces[2] << ")" << std::endl;
        os << std::fixed;
    }
}

// Writes out the postions of the ions and the current forces on them
void Ions::printForces(std::ostream& os, const int root) const
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
    int count                             = 0;
    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
    // std::cout<<"Num. local ions: "<<local_ions_.size()<<endl;

    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    int nproj                             = 0;
    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
    int nproj                             = 0;
    std::vector<Ion*>::const_iterator ion = overlappingNL_ions_.begin();
    while (ion != overlappingNL_ions_.end())
    {
        nproj += (*ion)->nProjectorsSubdomain();
        ion++;
    }
    return nproj;
}

Ion* Ions::findIon(const std::string& name) const
{
    std::vector<Ion*>::const_iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        bool same_name = (*ion)->compareName(name);
        if (same_name) return (*ion);

        ion++;
    }

    // return 0 if name not found in list_ions_
    return nullptr;
}

bool Ions::isLocal(const std::string& name) const
{
    for (auto ion : local_ions_)
    {
        bool same_name = ion->compareName(name);
        if (same_name) return true;
    }
    return false;
}

Ion* Ions::findIon(const int index) const
{
    for (auto ion : list_ions_)
    {
        if (ion->compareIndex(index)) return ion;
    }
    return nullptr;
}

Ion* Ions::findLocalIon(const int index) const
{
    for (auto ion : local_ions_)
    {
        if (ion->compareIndex(index)) return ion;
    }
    return nullptr;
}

void Ions::setPositions(const std::vector<double>& tau)
{
    assert(tau.size() == 3 * local_ions_.size());

    std::vector<Ion*>::iterator ion = local_ions_.begin();
    int ia                          = 0;
    while (ion != local_ions_.end())
    {
        (*ion)->setPosition(tau[3 * ia + 0], tau[3 * ia + 1], tau[3 * ia + 2]);

        ion++;
        ia++;
    }

    setup_ = false;
}

void Ions::get_positions(std::vector<std::vector<double>>& rr) const
{
    assert(rr.size() == local_ions_.size());
    if (local_ions_.empty()) return;
    std::vector<double> tau(3);
    int i                                 = 0;
    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
void Ions::set_positions(const std::vector<std::vector<double>>& rr)
{
    assert(rr.size() == local_ions_.size());

    if (local_ions_.empty()) return;
    std::vector<Ion*>::iterator ion = local_ions_.begin();
    int i                           = 0;
    while (ion != local_ions_.end())
    {
        assert(rr[i].size() == 3);
        (*ion)->setPosition(rr[i][0], rr[i][1], rr[i][2]);
        ion++;
        i++;
    }
}
void Ions::get_forces(std::vector<std::vector<double>>& ff) const
{
    assert(ff.size() == local_ions_.size());

    if (local_ions_.empty()) return;
    int i                                 = 0;
    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        ff[i][0] = (*ion)->force(0);
        ff[i][1] = (*ion)->force(1);
        ff[i][2] = (*ion)->force(2);
        ion++;
        i++;
    }
}
void Ions::set_forces(const std::vector<std::vector<double>>& ff)
{
    assert(ff.size() == local_ions_.size());

    if (local_ions_.empty()) return;
    std::vector<Ion*>::iterator ion = local_ions_.begin();
    int i                           = 0;
    while (ion != local_ions_.end())
    {
        (*ion)->setForce(ff[i][0], ff[i][1], ff[i][2]);
        ion++;
        i++;
    }
}

int Ions::readAtoms(const std::string& filename, const bool cell_relative)
{

    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Ions::readAtoms() --- Read "
                         << " atomic positions from file " << filename
                         << std::endl;

    std::string strxyz(".xyz");
    size_t found = filename.find(strxyz);
    if (found != std::string::npos)
    {
        num_ions_ = readAtomsFromXYZ(filename, cell_relative);
    }
    else
    {
        num_ions_ = readNatoms(filename, cell_relative);
    }

    return num_ions_;
}

int Ions::readAtoms(std::ifstream* tfile, const bool cell_relative)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    std::string first_string_read = "";
    if (mmpi.PE0())
    {
        (*tfile) >> first_string_read;
    }
    mmpi.bcastGlobal(first_string_read);

    std::string strxyz(".xyz");
    size_t found = first_string_read.find(strxyz);
    if (found != std::string::npos)
    {
        if (mmpi.PE0())
            (*MPIdata::sout) << "Read atomic positions from file "
                             << first_string_read << std::endl;
        num_ions_ = readAtomsFromXYZ(first_string_read, cell_relative);
    }
    else
    {
        if (mmpi.PE0())
        {
            for (unsigned short i = 0; i < first_string_read.size(); i++)
                tfile->unget();
        }
        num_ions_ = readNatoms(tfile, cell_relative);
    }

    return num_ions_;
}

int Ions::readAtomsFromXYZ(
    const std::string& filename, const bool cell_relative)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    Control& ct(*(Control::instance()));

    // set up list boundaries
    // get radius of projectors
    double rmax = getMaxListRadius();
    if (onpe0)
        (*MPIdata::sout) << "Max. radius for species in XYZ file: " << rmax
                         << std::endl;
    assert(rmax > 0.);

    setupListIonsBoundaries(rmax);

    int natoms = -1;

    std::ifstream* tfile = nullptr;
    if (mmpi.PE0())
    {
        tfile = new std::ifstream(filename.data(), std::ios::in);
        if (!(*tfile))
        {
            (*MPIdata::serr)
                << "Ions::readAtomsFromXYZ() --- ERROR: cannot open "
                << filename << std::endl;
        }
        else
        {
            std::string query;
            if (getline(*tfile, query))
            {
                std::stringstream na(query);
                na >> natoms;
            }
        }
        if (natoms < 1)
            (*MPIdata::sout)
                << "WARNING: Ions::readAtomsFromXYZ(), number of atoms read = "
                << natoms << std::endl;
    }
    mmpi.bcastGlobal(&natoms);
    if (natoms < 0) return natoms;

    std::vector<double> crds(3 * natoms);
    std::vector<short> spec(natoms);

    int count = 0;
    if (mmpi.PE0())
    {
        std::string query;
        getline(*tfile, query); // read comment line
        // read atomic species and positions
        for (int ia = 0; ia < natoms; ++ia)
        {
            if (!getline(*tfile, query))
            {
                break;
            }
            std::stringstream ss(query);
            std::string name_read;
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
        std::vector<Species>::const_iterator it = species_.begin();
        int isp                                 = -1;
        while (it != species_.end())
        {
            ++isp;
            if (it->getAtomicNumber() == spec[ia])
            {
                break;
            }
            ++it;
        }
        std::string spname("");
        for (std::map<std::string, short>::iterator itr = map_species_.begin();
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
                             << spec[ia] << std::endl;
            return -1;
        }

        // make a name for atom based on species and order of reading in
        std::string aname(spname);
        std::stringstream ss;
        ss << ia;
        if (ia < 10) aname.append("0");
        if (ia < 100) aname.append("0");
        if (ia < 1000) aname.append("0");
        if (ia < 10000) aname.append("0");
        aname.append(ss.str());
        Ion* new_ion
            = new Ion(species_[isp], aname, &crds[3 * ia], velocity, locked);
        new_ion->bcast(mmpi.commGlobal());

        // Populate list_ions_ list
        // std::cout<<"crds: "<<crds[3*ia+0]<<", "<<crds[3*ia+1]<<",
        // "<<crds[3*ia+2]<<endl;
        if (inListIons(crds[3 * ia + 0], crds[3 * ia + 1], crds[3 * ia + 2]))
        {
            list_ions_.push_back(new_ion);
            if (ct.verbose > 2)
                (*MPIdata::sout)
                    << "Ion " << aname << " at position " << crds[3 * ia + 0]
                    << "," << crds[3 * ia + 1] << "," << crds[3 * ia + 2]
                    << " added to the list... on PE" << mmpi.mypeGlobal()
                    << std::endl;
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
                        << mmpi.mypeGlobal() << std::endl;
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
    //    std::cout<<mmpi.mype()<<"...list size = "<<list_ions_.size()<<" local
    //    ions size = "<<local_ions_.size()<<endl;

    if (mmpi.PE0()) delete tfile;

    return natoms;
}

int Ions::readNatoms(const std::string& filename, const bool cell_relative)
{
    Control& ct(*(Control::instance()));

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Ions::readNAtoms() --- Read "
                         << " atomic positions from file " << filename
                         << std::endl;

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    std::ifstream* tfile = nullptr;
    if (mmpi.instancePE0())
    {
        tfile = new std::ifstream(filename.data(), std::ios::in);
        if (!tfile->is_open())
        {
            (*MPIdata::serr)
                << " Unable to open file " << filename.data() << std::endl;
            return -1;
        }
        else
        {
            (*MPIdata::sout) << "Open " << filename.data() << std::endl;
        }
    }

    int nread = readNatoms(tfile, cell_relative);

    if (mmpi.instancePE0())
    {
        if (ct.verbose > 0)
            (*MPIdata::sout) << "Close " << filename.data() << std::endl;
        tfile->close();
        delete tfile;
    }

    return nread;
}

int Ions::readNatoms(std::ifstream* tfile, const bool cell_relative)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    if (mmpi.PE0()) assert(tfile != nullptr);

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
                         << " atomic positions..." << std::endl;
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
                             << " atomic positions..." << std::endl;
    }
    if (onpe0) read_comments(*tfile);

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout) << "Ions::readNatoms() --- read " << nread
                         << " atomic positions..." << std::endl;

    int siz = local_ions_.size();
    mmpi.allreduce(&siz, 1, MPI_SUM);
    assert(siz == nread);

    return nread;
}

// return number of atoms read, 0 if end of file, or -1 if failure occurs
int Ions::read1atom(std::ifstream* tfile, const bool cell_relative)
{
    std::string name_read = "";
    // short  isp=0;
    double crds[3];
    double velocity[3] = { 0., 0., 0. };

    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    Control& ct(*(Control::instance()));

    short movable = 0;
    std::string query;
    short count = 1;
    if (mmpi.PE0())
    {
        if (getline(*tfile, query))
        {
            if (query[0] == '#') count = 0; // comment
            if (query.empty()) count = 0; // end of line only
            if (query.size() < 3) count = 0; // white spaces only
            // std::cout<<"length="<<query.size()<<endl;
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
                    << std::endl;
                count = -1;
            }
        }
        // std::cout<<"query="<<query<<endl;
        // std::cout<<"count="<<count<<endl;
    }

    mmpi.bcastGlobal(&count, 1);
    if (count < 1) return count;

    if (mmpi.PE0())
    {
        std::stringstream ss(query);
        ss >> name_read;
        if (!checkValidName(name_read))
        {
            std::cerr << "ERROR: Invalid name read in input file: " << name_read
                      << std::endl;
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
            (*MPIdata::serr)
                << "Atom " << name_read << ", should be movable (1) or not(0)"
                << std::endl;
            return -1;
        }
#ifdef DEBUG
        (*MPIdata::sout) << "movable=" << movable << std::endl;
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
                         << velocity[1] << "," << velocity[2] << ")"
                         << std::endl;
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
    std::string name(name_read);
    stripName(name);

    auto search = map_species_.find(name);
    assert(search != map_species_.end());

    short spec_nb                           = search->second;
    std::vector<Species>::const_iterator it = species_.begin();
    int isp                                 = -1;
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
    if (onpe0) (*MPIdata::sout) << "Create new Ion..." << std::endl;
#endif

    // Ion needs to be created by each MPI task to set global ids
    Ion* new_ion = new Ion(species_[isp], name_read, crds, velocity, locked);

    // assert( (new_ion->locked() == false) || (new_ion->locked() == true) );

#ifdef DEBUG
    if (onpe0) (*MPIdata::sout) << "Ion read..." << std::endl;
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
    double val = 0.;
    for (auto ion : local_ions_)
    {
        val += ion->getZion();
    }
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&val, 1, MPI_SUM);

    return (int)val;
}

double Ions::computeIonicCharge() const
{
    double ionic_charge = 0.;
    for (auto ion : local_ions_)
    {
        ionic_charge += (double)ion->getZion();
    }
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    mmpi.allreduce(&ionic_charge, 1, MPI_SUM);
    return ionic_charge;
}

void Ions::setVelocities(const std::vector<double>& tau0,
    const std::vector<double>& taup, const double dt)
{
    assert(tau0.size() == 3 * local_ions_.size());
    assert(taup.size() == 3 * local_ions_.size());

    int ia                           = 0;
    std::vector<Ion*>::iterator iion = local_ions_.begin();
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

void Ions::getPositions(std::vector<double>& tau) const
{
    assert(tau.size() == 3 * local_ions_.size());

    int ia                                 = 0;
    std::vector<Ion*>::const_iterator iion = local_ions_.begin();
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

    int ia                                 = 0;
    std::vector<Ion*>::const_iterator iion = local_ions_.begin();
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

    int ia                                 = 0;
    std::vector<Ion*>::const_iterator iion = local_ions_.begin();
    while (iion != local_ions_.end())
    {
        (*iion)->setPosition(
            tau0_[3 * ia + 0], tau0_[3 * ia + 1], tau0_[3 * ia + 2]);
        iion++;
        ia++;
    }
}

void Ions::getVelocities(std::vector<double>& tau) const
{
    assert(tau.size() == 3 * local_ions_.size());

    int ia                                 = 0;
    std::vector<Ion*>::const_iterator iion = local_ions_.begin();
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

void Ions::getForces(std::vector<double>& tau) const
{
    assert(tau.size() == 3 * local_ions_.size());

    int ia                                 = 0;
    std::vector<Ion*>::const_iterator iion = local_ions_.begin();
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
    std::vector<Ion*>::iterator ion = list_ions_.begin();
    while (ion != list_ions_.end())
    {
        /* Generate range of indices over which the short-range difference */
        /* potential will be mapped onto the global grid.                  */
        (*ion)->set_lstart(0, origin[0], h0);
        (*ion)->set_lstart(1, origin[1], h1);
        (*ion)->set_lstart(2, origin[2], h2);

        // Generate indices
        std::vector<int> Ai0;
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
            std::vector<int> Ai1;
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
                std::vector<int> Ai2;
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

    std::vector<Ion*>::const_iterator iion = local_ions_.begin();
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

    std::vector<Ion*>::const_iterator iion = local_ions_.begin();
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

void Ions::gatherNames(std::map<int, std::string>& names, const int root,
    const MPI_Comm comm) const
{
    std::vector<int> indexes(num_ions_, 0);
    std::vector<int> local_indexes;
    std::vector<std::string> data;
    std::vector<std::string> local_names;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
        const unsigned int num_ions = data.size();
        assert(num_ions == indexes.size());
        for (unsigned int i = 0; i < num_ions; i++)
        {
            const int ion_index        = indexes[i];
            const std::string ion_name = data[i];
            names.insert(std::pair<int, std::string>(ion_index, ion_name));
        }
    }
}

void Ions::gatherNames(
    std::vector<std::string>& names, const int root, const MPI_Comm comm) const
{
    std::vector<std::string> local_names;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        local_names.push_back((*ion)->name());
        ++ion;
    }

    // gather data to PE root
    std::vector<std::string> data;
    mgmol_tools::gatherV(local_names, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    names.clear();
    if (mype == root) names = data;
}

void Ions::gatherLockedNames(
    std::vector<std::string>& names, const int root, const MPI_Comm comm) const
{
    std::vector<std::string> local_names;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        if ((*ion)->locked()) local_names.push_back((*ion)->name());
        ++ion;
    }

    // gather data to PE root
    std::vector<std::string> data;
    mgmol_tools::gatherV(local_names, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    names.clear();
    if (mype == root) names = data;
}

void Ions::gatherIndexes(
    std::vector<int>& indexes, const int root, const MPI_Comm comm) const
{
    std::vector<int> local_indexes;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        local_indexes.push_back((*ion)->index());
        ++ion;
    }

    // gather data to PE root
    std::vector<int> data;
    mgmol_tools::gatherV(local_indexes, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    indexes.clear();
    if (mype == root) indexes = data;
}

void Ions::gatherNLprojIds(
    std::vector<int>& nlprojids, const int root, const MPI_Comm comm) const
{
    std::vector<int> local_nlprojids;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        local_nlprojids.push_back((*ion)->nlprojid());
        ++ion;
    }

    // gather data to PE root
    std::vector<int> data;
    mgmol_tools::gatherV(local_nlprojids, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    nlprojids.clear();
    if (mype == root) nlprojids = data;
}

void Ions::gatherAtomicNumbers(
    std::vector<int>& atnumbers, const int root, const MPI_Comm comm) const
{
    std::vector<int> local_atnumbers;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions_.end())
    {
        local_atnumbers.push_back((*ion)->atomic_number());
        ++ion;
    }

    // gather data to PE root
    std::vector<int> data;
    mgmol_tools::gatherV(local_atnumbers, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    atnumbers.clear();
    if (mype == root) atnumbers = data;
}

void Ions::gatherRandStates(std::vector<unsigned short>& rstates,
    const int root, const MPI_Comm comm) const
{
    std::vector<unsigned short> local_rstates;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        local_rstates.push_back((*ion)->randomState(0));
        local_rstates.push_back((*ion)->randomState(1));
        local_rstates.push_back((*ion)->randomState(2));
        ++ion;
    }

    // gather data to PE root
    std::vector<unsigned short> data;
    mgmol_tools::gatherV(local_rstates, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    rstates.clear();
    if (mype == root) rstates = data;
}

void Ions::gatherPositions(
    std::vector<double>& positions, const int root, const MPI_Comm comm) const
{
    std::vector<double> local_positions;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
    std::vector<double> data;
    mgmol_tools::gatherV(local_positions, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    positions.clear();
    if (mype == root) positions = data;
}

void Ions::gatherForces(
    std::vector<double>& forces, const int root, const MPI_Comm comm) const
{
    std::vector<double> local_forces;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
    std::vector<double> data;
    mgmol_tools::gatherV(local_forces, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    forces.clear();
    if (mype == root) forces = data;
}

void Ions::gatherVelocities(
    std::vector<double>& velocities, const int root, const MPI_Comm comm) const
{
    std::vector<double> local_velocities;

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
    while (ion != local_ions().end())
    {
        local_velocities.push_back((*ion)->velocity(0));
        local_velocities.push_back((*ion)->velocity(1));
        local_velocities.push_back((*ion)->velocity(2));

        ++ion;
    }

    // gather data to PE root
    std::vector<double> data;
    mgmol_tools::gatherV(local_velocities, data, root, comm);

    int mype = 0;
    MPI_Comm_rank(comm, &mype);
    velocities.clear();
    if (mype == root) velocities = data;
}

void Ions::gatherPositions(std::vector<double>& positions, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    positions.resize(3 * num_ions_, 0.);

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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

void Ions::gatherForces(std::vector<double>& forces, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    forces.resize(3 * num_ions_, 0.);

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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

void Ions::gatherLockedData(std::vector<int>& locked_data, const int root) const
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));

    locked_data.resize(num_ions_, 0);

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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

    std::vector<Ion*>::const_iterator ion = local_ions_.begin();
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
void Ions::syncNames(const int nions, std::vector<std::string>& local_names, vector<string>& names)
{
    MGmol_MPI& mmpi ( *(MGmol_MPI::instance()) );
    MPI_Comm& comm=mmpi.comm();  
    // prepare to gather data
    int npes = mmpi.size();
    int num_ions = nions;
    std::vector<int>recvcnts(npes);
    std::vector<int>::iterator rcv;
    std::vector<int>disp(npes);   
    std::vector<int>::iterator disp_it; 
    int tot = 0;
    int totchars;
    // Gather data for atom names
    // first get length of each name
    std::vector<int>loc_name_len;
    std::vector<int>nameLen(num_ions,0);
    tot=0;
    for(std::vector<std::string>::iterator str=local_names.begin(); str!=local_names.end(); ++str){
       //tot+= *str.size();
       std::string s = *str;
       tot+= s.length();
       loc_name_len.push_back(s.length());
    }
    mmpi.allGatherV(loc_name_len, nameLen);

    std::vector<char>char_loc_names(tot);
    int idx = 0;
    for(std::vector<std::string>::iterator str=local_names.begin(); str!=local_names.end(); ++str){
       //tot+= *str.size();
       std::string s = *str;
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
    std::vector<char> names_data(totchars);      
    std::vector<char>::iterator locnames_it=char_loc_names.begin();
    std::vector<char>::iterator names_it=names_data.begin();  
    disp_it=disp.begin();
    rcv = recvcnts.begin();  
    MPI_Allgatherv(&(*locnames_it), tot, MPI_CHAR, &(*names_it), &(*rcv), &(*disp_it), MPI_CHAR, comm);
      
//    std::vector<std::string> names;
    int pos = 0;
    for(int i=0; i<num_ions; i++)
    {
       char str[nameLen[i]+1];
       str[nameLen[i]] = '\0';
       memcpy(str, &names_data[pos], nameLen[i]*sizeof(char));
       std::string cstr;
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
    double radius                            = 0;
    std::vector<Species>::const_iterator spi = species_.begin();
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
    double radius                            = 0;
    std::vector<Species>::const_iterator spi = species_.begin();
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
    lstep_[0]
        = std::min((int)(ceil(rmax / (proc_width[0]))), (nproc_xyz[0] - 1));
    rstep_[0] = std::min(lstep_[0], nproc_xyz[0] - lstep_[0] - 1);

    /* y-direction */
    lstep_[1]
        = std::min((int)(ceil(rmax / (proc_width[1]))), (nproc_xyz[1] - 1));
    rstep_[1] = std::min(lstep_[1], nproc_xyz[1] - lstep_[1] - 1);
    /* z-direction */
    lstep_[2]
        = std::min((int)(ceil(rmax / (proc_width[2]))), (nproc_xyz[2] - 1));
    rstep_[2] = std::min(lstep_[2], nproc_xyz[2] - lstep_[2] - 1);
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
    // std::cout<<"inLocalIons..."<<endl;
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
//            std::cout<<"Ion "<<ion->name()<<" is here on multiple
//            tasks"<<endl;
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
    const int locSize, const int maxLocSize, std::vector<IonData>& iondata,
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
        std::vector<IonData>::iterator idata = ions_data_.begin() + *offset;
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
    for (std::vector<std::string>::const_iterator it
         = interacting_names_.begin();
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
    for (std::vector<std::string>::const_iterator it
         = interacting_names_.begin();
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
    // assume local_ions size is same as size of ions_names std::vector
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
    // Note: this is based on data from MD std::vectors

    // First cleanup list_ions_
    local_ions_.clear();
    // delete current ions from list
    std::vector<Ion*>::iterator ion = list_ions_.begin();
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
        std::vector<IonData> data2send(ions_data_);

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
    std::vector<IonData>::const_iterator idata = ions_data_.begin() + offset;
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

    std::vector<Ion*>::iterator lion = local_ions_.begin();
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
    std::vector<Ion*>::iterator ion = local_ions_.begin();
    int ia                          = 0;
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
    std::vector<Ion*>::iterator ion = local_ions_.begin();
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
        std::cout << "Ions::rescaleVelocities() with factor " << factor
                  << std::endl;
    }
    std::vector<Ion*>::iterator ion = local_ions_.begin();
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

    std::vector<Ion*>::const_iterator ion1 = local_ions_.begin();
    std::vector<Ion*>::const_iterator ion2 = local_ions_.begin();
    while (ion1 != local_ions_.end())
    {
        while (ion2 != local_ions_.end())
        {
            if (ion1 != ion2)
            {
                double d = (*ion1)->minimage(**ion2, lattice_, ct.bcPoisson);
                if (d < distance) distance = d;
            }
            ion2++;
        }
        ion1++;
    }

    return distance;
}
