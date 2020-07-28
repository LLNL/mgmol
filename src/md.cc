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
#include "DFTsolver.h"
#include "DMStrategy.h"
#include "Energy.h"
#include "ExtendedGridOrbitals.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "LocalizationRegions.h"
#include "MD_IonicStepper.h"
#include "MDfiles.h"
#include "MGmol.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "MasksSet.h"
#include "Mesh.h"
#include "OrbitalsExtrapolation.h"
#include "OrbitalsExtrapolationFactory.h"
#include "OrbitalsPreconditioning.h"
#include "Potentials.h"
#include "ProjectedMatricesMehrstellen.h"
#include "ProjectedMatricesSparse.h"
#include "Rho.h"
#include "Signal.h"
#include "SpreadsAndCenters.h"
#include "tools.h"

#include <sstream>
#include <string>
#include <vector>

Timer md_iterations_tm("md_iterations");
Timer md_tau_tm("md_tau");
Timer md_moveVnuc_tm("md_moveVnuc");
Timer md_updateMasks_tm("md_updateMasks");
Timer md_extrapolateOrbitals_tm("md_extrapolateOrbitals");

#define DUMP_MAX_NUM_TRY 5

template <class T>
void MGmol<T>::moveVnuc(Ions& ions)
{
    md_moveVnuc_tm.start();

    Potentials& pot = hamiltonian_->potential();

    // Update items that change when the ionic coordinates change
    pot.axpVcompToVh(1.);
    initNuc(ions);
    pot.axpVcompToVh(-1.);

    proj_matrices_->setHiterativeIndex(-1, -1);
    hamiltonian_->setHlOutdated();

    md_moveVnuc_tm.stop();
}

template <class T>
void MGmol<T>::preWFextrapolation()
{
    Control& ct = *(Control::instance());
    if (ct.OuterSolver() == OuterSolverType::ABPG
        || ct.OuterSolver() == OuterSolverType::NLCG)
    {
        if (ct.dm_mix < 1.) proj_matrices_->stripDM();
    }
    else
    {
#if EXTRAPOLATE_H
        ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>* projmat
            = dynamic_cast<
                ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>*>(
                proj_matrices_);
        assert(projmat);
        orbitals_extrapol_->initExtrapolationH(projmat->getMatHB());
#endif
    }
}

template <class T>
void MGmol<T>::postWFextrapolation(T* orbitals)
{
    Control& ct = *(Control::instance());
    if (ct.isLocMode())
    {
        orbitals->computeGramAndInvS();
    }

    if (ct.OuterSolver() == OuterSolverType::ABPG
        || ct.OuterSolver() == OuterSolverType::NLCG)
    {
        if (ct.dm_mix < 1.)
            proj_matrices_->dressupDM();
        else
            proj_matrices_->setDMto2InvS();
    }
    else
    {
        if (orbitals_extrapol_->extrapolatedH()) dm_strategy_->update();
    }
}

template <class T>
void MGmol<T>::extrapolate_centers(bool small_move)
{
    assert(lrs_);

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0)
        os_ << "Extrapolate LR centers..." << std::endl;

    if (small_move && ct.lrs_extrapolation == 2)
    {
        lrs_->extrapolateCentersQuadratic();
    }
    else
    {
        // for large moves, extrapolation does nothing
        // but save current centers
        if (ct.lrs_extrapolation == 1 || ct.lrs_extrapolation == 2)
            lrs_->extrapolateCentersLinear(small_move);

        else if (ct.lrs_extrapolation == 10)
            lrs_->extrapolateCentersVerlet(small_move);
    }
}

template <class T>
T* MGmol<T>::new_orbitals_with_current_LRs(bool setup)
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    Control& ct            = *(Control::instance());

    if (ct.verbose > 1 && onpe0)
        os_ << "Build new orbitals with new masks..." << std::endl;

    if (ct.lr_update) update_masks();

    // need to build new orbitals as masks have changed
    T* new_orbitals = new T("NewMasks", mygrid, mymesh->subdivx(), ct.numst,
        ct.bc, proj_matrices_, lrs_, currentMasks_, corrMasks_, local_cluster_,
        setup);

    return new_orbitals;
}

template <class T>
void MGmol<T>::update_orbitals_LRs(T** orbitals)
{
    T* new_orbitals = new_orbitals_with_current_LRs();
    new_orbitals->assign(**orbitals);
    delete (*orbitals);
    (*orbitals) = new_orbitals;

    (*orbitals)->applyMask();
}

template <class T>
void MGmol<T>::extrapolate_orbitals(T** orbitals)
{
    md_extrapolateOrbitals_tm.start();

    T* new_orbitals = new_orbitals_with_current_LRs();

    orbitals_extrapol_->extrapolate_orbitals(orbitals, new_orbitals);

    (*orbitals)->incrementIterativeIndex();

    md_extrapolateOrbitals_tm.stop();
}

template <class T>
void MGmol<T>::move_orbitals(T** orbitals)
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 1) os_ << "Move orbitals..." << std::endl;

    T* new_orbitals = new_orbitals_with_current_LRs();

    // copy old data
    new_orbitals->assign(**orbitals);

    delete *orbitals;
    *orbitals = new_orbitals;

    if (ct.isLocMode())
    {
        (*orbitals)->applyMask();
        (*orbitals)->normalize();
    }
    (*orbitals)->incrementIterativeIndex();
}

template <class T>
int MGmol<T>::update_masks()
{
    assert(lrs_);

    Control& ct = *(Control::instance());
    if (!ct.isLocMode()) return 0;

    md_updateMasks_tm.start();

    if (onpe0 && ct.verbose > 1) os_ << "update_masks()" << std::endl;

    currentMasks_->update(lrs_);
    corrMasks_->update(lrs_);

    md_updateMasks_tm.stop();

    return 0;
}

void checkMaxForces(const std::vector<double>& fion,
    const std::vector<short>& atmove, std::ostream& os)
{
    assert(3 * atmove.size() == fion.size());

    const int na = (int)atmove.size();
    double f2    = 0.;
    for (int i = 0; i < na; i++)
    {
        if (atmove[i])
            f2 = std::max(f2, fion[3 * i + 0] * fion[3 * i + 0]
                                  + fion[3 * i + 1] * fion[3 * i + 1]
                                  + fion[3 * i + 2] * fion[3 * i + 2]);
    }
    if (onpe0)
        os << "Max. |force| = " << std::setprecision(4) << std::scientific
           << std::sqrt(f2) << std::endl;
}

template <class T>
int MGmol<T>::dumprestartFile(T** orbitals, Ions& ions, Rho<T>& rho,
    const bool write_extrapolated_wf, const short count)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    Control& ct              = *(Control::instance());
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    const pb::Grid& mygrid   = mymesh->grid();
    const unsigned gdim[3] = { mygrid.gdim(0), mygrid.gdim(1), mygrid.gdim(2) };

    std::string filename(std::string(ct.out_restart_file));
    // add an integer corresponding to attempt number/count
    // to allow several attempts at creating and writing file
    std::stringstream s;
    s << count;
    filename += s.str();

    HDFrestart h5file(filename, myPEenv, gdim, ct.out_restart_file_type);

    T previous_orbitals("ForDumping", **orbitals, false);
    if (!orbitals_extrapol_->getRestartData(previous_orbitals))
        previous_orbitals.assign(**orbitals);
    int ierr = write_hdf5(h5file, rho.rho_, ions, previous_orbitals, lrs_);
    mmpi.allreduce(&ierr, 1, MPI_MIN);

    if (ierr < 0)
    {
        if (onpe0)
            (*MPIdata::serr)
                << "dumprestartFile: cannot write ...previous_orbitals..."
                << std::endl;
        return ierr;
    }
    // write_hdf5(h5file, rho.rho_, ions, *orbitals_minus1);
    // stepper->write_hdf5(h5file);

    if (write_extrapolated_wf && ct.out_restart_info > 2)
    {
        ierr = (*orbitals)->write_func_hdf5(h5file, "ExtrapolatedFunction");
        mmpi.allreduce(&ierr, 1, MPI_MIN);
        if (ierr < 0)
        {
            if (onpe0)
                (*MPIdata::serr) << "dumprestartFile: cannot write "
                                    "...ExtrapolatedFunction..."
                                 << std::endl;
            return ierr;
        }
    }

    ierr = h5file.close();
    mmpi.allreduce(&ierr, 1, MPI_MIN);
    if (ierr < 0)
    {
        if (onpe0)
            (*MPIdata::serr)
                << "dumprestartFile: cannot close file..." << std::endl;
        return ierr;
    }

    return 0;
}

template <class T>
void MGmol<T>::md(T** orbitals, Ions& ions)
{
    Control& ct = *(Control::instance());

    // initialize stepper data
    std::vector<std::string>& ions_names(ions.getLocalNames());
    std::vector<double>& tau0(ions.getTau0()); // tau0[3*ia+j]
    std::vector<double>& taup(ions.getTaup()); // taup[3*ia+j]
    std::vector<double>& taum(ions.getTaum()); // taum[3*ia+j]
    std::vector<double>& fion(ions.getFion()); // fion[3*ia+j]
    std::vector<double>& pmass(ions.getPmass()); // pmass[ia]
    std::vector<short>& atmove(ions.getAtmove()); // atmove[ia]
    std::vector<unsigned short>& rand_states(
        ions.getRandStates()); // initial random states

    // initialize taum with velocities
    std::vector<double>& vel(ions.getVelocities());
    taum = vel;

    int size_tau = (int)tau0.size();
    DFTsolver<T>::resetItCount();

    orbitals_extrapol_
        = OrbitalsExtrapolationFactory<T>::create(ct.WFExtrapolation());

    MD_IonicStepper* stepper = new MD_IonicStepper(
        ct.dt, atmove, tau0, taup, taum, fion, pmass, rand_states);
    stepper->setThermostat(ct.thermostat_type, ct.tkel, ct.thtime, ct.thwidth,
        constraints_->size());

    constraints_->printConstraints(os_);

    if (ct.restart_info > 0 && !ct.override_restart)
    {
        if (onpe0) os_ << "Use restart file to initialize MD..." << std::endl;
        stepper->init(*h5f_file_);
    }
    else
    {
        if (onpe0) os_ << "Use input file to initialize MD..." << std::endl;

        // tau0: velocities->displacements
        double dt = -ct.dt;
        int ione  = 1;
        DSCAL(&size_tau, &dt, &tau0[0], &ione);

        // tau0: previous positions given displacement
        double one = 1.;
        DAXPY(&size_tau, &one, &taup[0], &ione, &tau0[0], &ione);

        // enforce constraints before 1st step
        constraints_->enforceConstraints(20);

        stepper->updateTau();
        ions.setPositions(tau0);

        ions.setup();
    }

    ions.printPositions(os_);

    if (ct.restart_info > 1)
    {
        int flag_extrapolated_data = 0;
        if (onpe0)
        {
            flag_extrapolated_data
                = h5f_file_->dset_exists("ExtrapolatedFunction0000");
            if (flag_extrapolated_data == 0)
                flag_extrapolated_data
                    = h5f_file_->dset_exists("ExtrapolatedFunction0");
        }
        MPI_Bcast(&flag_extrapolated_data, 1, MPI_INT, 0, comm_);

        if (ct.restart_info > 2)
        {
            if (flag_extrapolated_data)
            {
                if (onpe0) os_ << "Create new orbitals_minus1..." << std::endl;

                orbitals_extrapol_->setupPreviousOrbitals(&current_orbitals_,
                    proj_matrices_, lrs_, local_cluster_, currentMasks_,
                    corrMasks_, *h5f_file_);

                // need to reset a few things as we just read new orbitals
                (*orbitals)->computeGramAndInvS();
                dm_strategy_->update();
            }

            DFTsolver<T>::setItCountLarge();
        }

        // check if we are restarting from an MD dump
        // no extrapolated functions -> atomic positions were extrapolated
        // in stepper->init()
        if (!flag_extrapolated_data)
        {
            moveVnuc(ions);
        }

        delete h5f_file_;
        h5f_file_ = nullptr;
    }

    // additional SC steps to compensate random start
    if (ct.restart_info < 3)
    {
        double eks = 0.;
        quench(*orbitals, ions, ct.max_electronic_steps, 20, eks);
    }

    ct.max_changes_pot = 0;

    bool extrapolated_flag = true;
    if (ct.dt <= 0.) extrapolated_flag = false;

    MDfiles md_files;

    // main MD iteration loop
    for (int mdstep = 1; mdstep <= ct.num_MD_steps; mdstep++)
    {
        if (ct.checkTimeout())
        {
            if (onpe0)
                os_ << "Signal caught: No more MD iterations" << std::endl;
            break;
        }

        md_iterations_tm.start();

        double eks              = 0.;
        int retval              = 0;
        bool small_move         = true;
        bool last_move_is_small = true;
        do
        {
            retval = quench(*orbitals, ions, ct.max_electronic_steps, 0, eks);

            // update localization regions
            if (ct.adaptiveLRs())
            {
                assert(lrs_);
                adaptLR(spreadf_, nullptr);

                last_move_is_small = lrs_->moveIsSmall();

                // update cluster for load balancing
                if (ct.load_balancing_alpha > 0.0
                    && mdstep % ct.load_balancing_modulo == 0)
                {
                    local_cluster_->computeClusters(
                        ct.load_balancing_max_iterations);
                }

                // printWithTimeStamp("quench done...",cout);

                // do extra cycles if centers resulting from quench
                // are "significantly" different from initial centers
                if (!last_move_is_small)
                {
                    printWithTimeStamp(
                        "WARNING: large move->extra inner cycle...", std::cout);
                    small_move = false;
                    move_orbitals(orbitals);

                    (*orbitals)->computeGramAndInvS();
                    dm_strategy_->update();

                    // reduce number of steps to keep total run time about the
                    // same
                    ct.num_MD_steps--;
                    // printWithTimeStamp("extra work done...",cout);
                }
            }
        } while (!last_move_is_small);

        if (retval < 0)
        {
            extrapolated_flag = false;
            if (ct.dt > 0.)
            {
                if (onpe0)
                {
                    os_ << "WARNING md(): quench returned value " << retval
                        << std::endl;
                    os_ << "STOP md..." << std::endl;
                }
                break;
            }
        }

        if (onpe0)
            os_ << std::setprecision(12) << std::fixed << "%%  "
                << md_iteration_ << "  IONIC CONFIGURATION ENERGY = " << eks
                << std::endl;

        // Compute forces
        force(**orbitals, ions);

        // set fion
        ions.getForces(fion);

        // constraints need to be added again and setup as atoms
        // may have moved and local atoms are not the same anymore
        constraints_->addConstraints(ions);

        constraints_->setup(ions);

        constraints_->printConstraintsForces(os_);

        constraints_->projectOutForces();

        // Print atomic coordinates and forces
        if (md_iteration_ % ct.md_print_freq == 0)
        {
            std::string zero = "0";
            if (zero.compare(ct.md_print_filename) != 0)
            {
                std::vector<float> spreads;
                std::vector<Vector3D> centers;
                std::vector<int> lgids;
                if (ct.isLocMode())
                {
                    spreadf_->computeLocalSpreads(spreads);
                    spreadf_->getLocalCenters(centers);
                    spreadf_->getLocalGids(lgids);
                }

                md_files.printDataInFiles(ions_names, tau0, fion, taum, centers,
                    spreads, lgids, md_iteration_, ct.dt);
            }

            // also print in stdout for small problems or high verbosity
            if (ions_->getNumIons() < 256 || ct.verbose > 2)
            {
                if (ct.verbose > 0) ions.printForcesGlobal(os_);
            }
            else if (zero.compare(ct.md_print_filename) == 0)
            {
                ions.printForcesLocal(os_);
            }
        }

        if (ct.dt <= 0.) checkMaxForces(fion, atmove, os_);

        // add random forces if Langevin dynamics
        stepper->addLangevinForces();
        constraints_->projectOutForces(50);

        // Compute taup
        stepper->run();
        constraints_->enforceConstraints(20);
        if (ct.enforceVmass0) ions.removeMassCenterMotion();
        if (ct.verbose > 2 && ct.dt > 0.) stepper->printVelocities(os_);

        // compute velocities and temperatures at t+1/2
        double ekin = 0.;
        if (ct.dt > 0.)
        {
            // T=Ek/(kB*3/2*N)
            const double temperature = stepper->temperature();
            ekin                     = stepper->kineticEnergy();
            if (onpe0)
            {
                os_ << std::setprecision(6);
                os_ << "Kinetic   Energy= " << ekin;
                os_ << std::setprecision(5);
                os_ << " (T= " << temperature << " )" << std::endl;
            }
        }

        if (onpe0)
        {
            os_ << std::setprecision(10);
            os_ << "Total     Energy= " << ekin + eks << std::endl;
        }

        stepper->updateTau();

        // update local information about ions
        ions.updateIons();

        ions.setup();

        ct.steps = md_iteration_;

#if EXTRAPOLATE_RHO
        if (onpe0) os_ << "Extrapolate rho..." << std::endl;
        rho_->axpyRhoc(-1., rhoc_);
        rho_->extrapolate();
#endif

        moveVnuc(ions);

#if EXTRAPOLATE_RHO
        rho_->axpyRhoc(1., rhoc_);
#endif

        if (!small_move)
        {
            orbitals_extrapol_->clearOldOrbitals();
            lrs_->clearOldCenters();
        }

        preWFextrapolation();

        if (ct.dt > 0.
            || ct.WFExtrapolation() == WFExtrapolationType::Reversible)
        {
            if (ct.lrs_extrapolation > 0) extrapolate_centers(small_move);

            extrapolate_orbitals(orbitals);
        }
        else
            move_orbitals(orbitals);

        postWFextrapolation(*orbitals);

        // rho_->setup(ct.orbital_type,(*orbitals)->getGlobalIndexes());
        // rho_->update(**orbitals);

        md_time_ += ct.dt;
        if (ct.dt > 0.) md_iteration_++;

        if (ct.checkpoint && ct.out_restart_file != "0")
            if (ct.out_restart_info > 0)
                if ((md_iteration_ % ct.checkpoint) == 0
                    && mdstep < ct.num_MD_steps)
                {
                    int ierr    = -1;
                    short count = 0;
                    while (ierr < 0 && count < DUMP_MAX_NUM_TRY)
                    {
                        dump_tm_.start();
                        ierr = dumprestartFile(
                            orbitals, ions, *rho_, extrapolated_flag, count);
                        dump_tm_.stop();
                        if (onpe0 && ierr < 0 && count < (DUMP_MAX_NUM_TRY - 1))
                            std::cout
                                << "dumprestartFile() failed... try again..."
                                << std::endl;
                        if (ierr < 0) sleep(1.);
                        count++;
                    }

                    printWithTimeStamp("dumped restart file...", std::cout);
                }

        md_iterations_tm.stop();

    } // md loop

    // final dump
    if (ct.out_restart_info > 0)
    {
        int ierr    = -1;
        short count = 0;
        while (ierr < 0 && count < DUMP_MAX_NUM_TRY)
        {
            dump_tm_.start();
            ierr = dumprestartFile(
                orbitals, ions, *rho_, extrapolated_flag, count);
            dump_tm_.stop();

            if (onpe0 && ierr < 0 && count < (DUMP_MAX_NUM_TRY - 1))
                std::cout << "dumprestartFile() failed... try again..."
                          << std::endl;
            if (ierr < 0) sleep(1.);
            count++;
        }

        printWithTimeStamp("dumped last restart file...", std::cout);
    }

    delete stepper;
    delete orbitals_extrapol_;
}

template class MGmol<LocGridOrbitals>;
template class MGmol<ExtendedGridOrbitals>;
