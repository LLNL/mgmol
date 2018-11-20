// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ConstraintSet.h"
#include "Control.h"
#include "DFTsolver.h"
#include "DMStrategy.h"
#include "Energy.h"
#include "Hamiltonian.h"
#include "Ions.h"
#include "KBPsiMatrix.h"
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
using namespace std;

static OrbitalsExtrapolation* orbitals_extrapol;

Timer md_iterations_tm("md_iterations");
Timer md_tau_tm("md_tau");
Timer md_moveVnuc_tm("md_moveVnuc");
Timer md_updateMasks_tm("md_updateMasks");
Timer md_extrapolateOrbitals_tm("md_extrapolateOrbitals");

#define DUMP_MAX_NUM_TRY 5

void MGmol::moveVnuc(Ions& ions)
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

void MGmol::preWFextrapolation()
{
    Control& ct = *(Control::instance());
    if (ct.it_algo_type <= 1)
    {
        if (ct.dm_mix < 1.) proj_matrices_->stripDM();
    }
    else
    {
#if EXTRAPOLATE_H
        ProjectedMatrices* projmat
            = dynamic_cast<ProjectedMatrices*>(proj_matrices_);
        assert(projmat);
        projmat->initExtrapolationH();
#endif
    }
}

void MGmol::postWFextrapolation(LocGridOrbitals* orbitals)
{
    Control& ct = *(Control::instance());
    if (ct.isLocMode())
    {
        orbitals->computeGramAndInvS();
    }

    if (ct.it_algo_type > 1)
    {
        if (orbitals_extrapol->extrapolatedH()) dm_strategy_->update();
    }
    else
    {
        if (ct.dm_mix < 1.)
            proj_matrices_->dressupDM();
        else
            proj_matrices_->setDMto2InvS();
    }
}

void MGmol::extrapolate_centers(bool small_move)
{
    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 0) os_ << "Extrapolate LR centers..." << endl;

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

LocGridOrbitals* MGmol::new_orbitals_with_current_LRs(bool setup)
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    Control& ct            = *(Control::instance());

    if (ct.verbose > 1 && onpe0)
        os_ << "Build new orbitals with new masks..." << endl;

    if (ct.lr_update) update_masks();

    // need to build new orbitals as masks have changed
    LocGridOrbitals* new_orbitals = new LocGridOrbitals("NewMasks", mygrid,
        mymesh->subdivx(), ct.numst, ct.bc, proj_matrices_, lrs_, currentMasks_,
        corrMasks_, local_cluster_, setup);

    return new_orbitals;
}

void MGmol::update_orbitals_LRs(LocGridOrbitals** orbitals)
{
    LocGridOrbitals* new_orbitals = new_orbitals_with_current_LRs();
    new_orbitals->assign(**orbitals);
    delete (*orbitals);
    (*orbitals) = new_orbitals;

    (*orbitals)->applyMask();
}

void MGmol::extrapolate_orbitals(LocGridOrbitals** orbitals)
{
    md_extrapolateOrbitals_tm.start();

    LocGridOrbitals* new_orbitals = new_orbitals_with_current_LRs();

    orbitals_extrapol->extrapolate_orbitals(orbitals, new_orbitals);

    (*orbitals)->incrementIterativeIndex();

    md_extrapolateOrbitals_tm.stop();
}

void MGmol::move_orbitals(LocGridOrbitals** orbitals)
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 1) os_ << "Move orbitals..." << endl;

    LocGridOrbitals* new_orbitals = new_orbitals_with_current_LRs();

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

int MGmol::update_masks()
{
    Control& ct = *(Control::instance());
    if (!ct.isLocMode()) return 0;

    md_updateMasks_tm.start();

    if (onpe0 && ct.verbose > 1) os_ << "update_masks()" << endl;

    currentMasks_->update(*lrs_);
    corrMasks_->update(*lrs_);

    md_updateMasks_tm.stop();

    return 0;
}

void checkMaxForces(
    const vector<double>& fion, const vector<short>& atmove, ostream& os)
{
    assert(3 * atmove.size() == fion.size());

    const int na = (int)atmove.size();
    double f2    = 0.;
    for (int i = 0; i < na; i++)
    {
        if (atmove[i])
            f2 = max(f2, fion[3 * i + 0] * fion[3 * i + 0]
                             + fion[3 * i + 1] * fion[3 * i + 1]
                             + fion[3 * i + 2] * fion[3 * i + 2]);
    }
    if (onpe0)
        os << "Max. |force| = " << setprecision(4) << scientific << sqrt(f2)
           << endl;
}

int MGmol::dumprestartFile(LocGridOrbitals** orbitals, Ions& ions, Rho& rho,
    const bool write_extrapolated_wf, const short count)
{
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    Control& ct              = *(Control::instance());
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    const pb::Grid& mygrid   = mymesh->grid();
    const unsigned gdim[3] = { mygrid.gdim(0), mygrid.gdim(1), mygrid.gdim(2) };

    string filename(string(ct.out_restart_file));
    // add an integer corresponding to attempt number/count
    // to allow several attempts at creating and writing file
    stringstream s;
    s << count;
    filename += s.str();

    HDFrestart h5file(filename, myPEenv, gdim, ct.out_restart_file_type);

    LocGridOrbitals previous_orbitals("ForDumping", **orbitals, false);
    if (!orbitals_extrapol->getRestartData(previous_orbitals))
        previous_orbitals.assign(**orbitals);
    int ierr = write_hdf5(h5file, rho.rho_, ions, previous_orbitals, *lrs_);
    mmpi.allreduce(&ierr, 1, MPI_MIN);

    if (ierr < 0)
    {
        if (onpe0)
            (*MPIdata::serr)
                << "dumprestartFile: cannot write ...previous_orbitals..."
                << endl;
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
                                 << endl;
            return ierr;
        }
    }

    ierr = h5file.close();
    mmpi.allreduce(&ierr, 1, MPI_MIN);
    if (ierr < 0)
    {
        if (onpe0)
            (*MPIdata::serr) << "dumprestartFile: cannot close file..." << endl;
        return ierr;
    }

    return 0;
}

void MGmol::md(LocGridOrbitals** orbitals, Ions& ions)
{
    Control& ct = *(Control::instance());

    // initialize stepper data
    vector<string>& ions_names(ions.getLocalNames());
    vector<double>& tau0(ions.getTau0()); // tau0[3*ia+j]
    vector<double>& taup(ions.getTaup()); // taup[3*ia+j]
    vector<double>& taum(ions.getTaum()); // taum[3*ia+j]
    vector<double>& fion(ions.getFion()); // fion[3*ia+j]
    vector<double>& pmass(ions.getPmass()); // pmass[ia]
    vector<short>& atmove(ions.getAtmove()); // atmove[ia]
    vector<unsigned short>& rand_states(
        ions.getRandStates()); // initial random states

    // initialize taum with velocities
    vector<double>& vel(ions.getVelocities());
    taum = vel;

    int size_tau = (int)tau0.size();
    DFTsolver::resetItCount();

    orbitals_extrapol
        = OrbitalsExtrapolationFactory::create(ct.wf_extrapolation, NULL);

    MD_IonicStepper* stepper = new MD_IonicStepper(
        ct.dt, atmove, tau0, taup, taum, fion, pmass, rand_states);
    stepper->setThermostat(ct.thermostat_type, ct.tkel, ct.thtime, ct.thwidth,
        constraints_->size());

    constraints_->printConstraints(os_);

    if (ct.restart_info > 0 && !ct.override_restart)
    {
        if (onpe0) os_ << "Use restart file to initialize MD..." << endl;
        stepper->init(*h5f_file_);
    }
    else
    {
        if (onpe0) os_ << "Use input file to initialize MD..." << endl;

        // tau0: velocities->displacements
        double dt = -ct.dt;
        int ione  = 1;
        dscal(&size_tau, &dt, &tau0[0], &ione);

        // tau0: previous positions given displacement
        double one = 1.;
        daxpy(&size_tau, &one, &taup[0], &ione, &tau0[0], &ione);

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
#ifdef USE_MPI
        MPI_Bcast(&flag_extrapolated_data, 1, MPI_INT, 0, comm_);
#endif

        if (ct.restart_info > 2)
        {
            if (flag_extrapolated_data)
            {
                if (onpe0) os_ << "Create new orbitals_minus1..." << endl;

                orbitals_extrapol->setupPreviousOrbitals(&current_orbitals_,
                    proj_matrices_, lrs_, local_cluster_, currentMasks_,
                    corrMasks_, *h5f_file_);

                // need to reset a few things as we just read new orbitals
                (*orbitals)->computeGramAndInvS();
                dm_strategy_->update();
            }

            DFTsolver::setItCountLarge();
        }

        // check if we are restarting from an MD dump
        // no extrapolated functions -> atomic positions were extrapolated
        // in stepper->init()
        if (!flag_extrapolated_data)
        {
            moveVnuc(ions);
        }

        delete h5f_file_;
        h5f_file_ = 0;
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
            if (onpe0) os_ << "Signal caught: No more MD iterations" << endl;
            break;
        }

        md_iterations_tm.start();

        double eks      = 0.;
        int retval      = 0;
        bool small_move = true;
        do
        {
            retval = quench(*orbitals, ions, ct.max_electronic_steps, 0, eks);

            // update localization regions
            if (ct.adaptiveLRs())
            {
                adaptLR(spreadf_, 0);

                // update cluster for load balancing
                if (ct.load_balancing_alpha > 0.0
                    && mdstep % ct.load_balancing_modulo == 0)
                {
                    local_cluster_->computeClusters(
                        ct.load_balancing_max_iterations);
                }
            }

            // printWithTimeStamp("quench done...",cout);

            // do extra cycles if centers resulting from quench
            // are "significantly" different from initial centers
            if (!lrs_->moveIsSmall())
            {
                printWithTimeStamp(
                    "WARNING: large move->extra inner cycle...", cout);
                small_move = false;
                move_orbitals(orbitals);

                (*orbitals)->computeGramAndInvS();
                dm_strategy_->update();

                // reduce number of steps to keep total run time about the same
                ct.num_MD_steps--;
                // printWithTimeStamp("extra work done...",cout);
            }
        } while (!lrs_->moveIsSmall());

        if (retval < 0)
        {
            extrapolated_flag = false;
            if (ct.dt > 0.)
            {
                if (onpe0)
                {
                    os_ << "WARNING md(): quench returned value " << retval
                        << endl;
                    os_ << "STOP md..." << endl;
                }
                break;
            }
        }

        if (onpe0)
            os_ << setprecision(12) << fixed << "%%  " << md_iteration_
                << "  IONIC CONFIGURATION ENERGY = " << eks << endl;

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
            string zero = "0";
            if (zero.compare(ct.md_print_filename) != 0)
            {
                vector<float> spreads;
                vector<Vector3D> centers;
                vector<int> lgids;
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
                os_ << setprecision(6);
                os_ << "Kinetic   Energy= " << ekin;
                os_ << setprecision(5);
                os_ << " (T= " << temperature << " )" << endl;
            }
        }

        if (onpe0)
        {
            os_ << setprecision(10);
            os_ << "Total     Energy= " << ekin + eks << endl;
        }

        stepper->updateTau();

        // update local information about ions
        ions.updateIons();

        ions.setup();

        ct.steps = md_iteration_;

#if EXTRAPOLATE_RHO
        if (onpe0) os_ << "Extrapolate rho..." << endl;
        rho_->axpyRhoc(-1., rhoc_);
        rho_->extrapolate();
#endif

        moveVnuc(ions);

#if EXTRAPOLATE_RHO
        rho_->axpyRhoc(1., rhoc_);
#endif

        if (!small_move)
        {
            orbitals_extrapol->clearOldOrbitals();
            lrs_->clearOldCenters();
        }

        preWFextrapolation();

        if (ct.dt > 0. || ct.wf_extrapolation == 0)
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
                            cout << "dumprestartFile() failed... try again..."
                                 << endl;
                        if (ierr < 0) sleep(1.);
                        count++;
                    }

                    printWithTimeStamp("dumped restart file...", cout);
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
                cout << "dumprestartFile() failed... try again..." << endl;
            if (ierr < 0) sleep(1.);
            count++;
        }

        printWithTimeStamp("dumped last restart file...", cout);
    }

    delete stepper;
    delete orbitals_extrapol;
}
