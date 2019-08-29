// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_IONS_H
#define MGMOL_IONS_H

#include <fstream>
#include <iomanip>
#include <list>
#include <map>
#include <vector>

#include "DataDistribution.h"
#include "DistributedIonicData.h"
#include "Ion.h"
#include "hdf5.h"

class HDFrestart;

class Ions
{
    static std::map<std::string, short> map_species_;
    static int num_ions_;

    // max l over all the ions
    static short max_num_proj_;

    static double max_Vl_radius_;
    static double max_Vnl_radius_;

    const std::vector<Species>& species_;

    std::vector<Ion*> list_ions_;
    std::vector<Ion*> local_ions_; // centered in local sub-domain
    std::vector<Ion*> interacting_ions_; // for ion-ion interactions
    std::vector<Ion*>
        overlappingNL_ions_; // with projectors overlapping local sub-domain
    std::vector<Ion*> overlappingVL_ions_; // with local potential overlapping
                                           // local sub-domain

    double lattice_[3];
    double div_lattice_[3];

    // boundaries to determine which ions are known on local processor
    double list_boundary_left_[3];
    double list_boundary_right_[3];

    int source_[3];
    int dest_[3];
    MPI_Comm cart_comm_; // MPI cartesian communicator for data distribution

    bool has_locked_atoms_;

    void readRestartVelocities(HDFrestart& h5_file);
    void readRestartRandomStates(HDFrestart& h5_file);
    void readRestartPositions(HDFrestart& h5_file);
    int read1atom(std::ifstream* tfile, const bool cell_relative);

    // void associate2PE();
    void setupInteractingIons();
    void setupListOverlappingIons();
    void setMapVL();
    double computeMaxVlRadius() const;
    double computeMaxNLprojRadius() const;
    double getMaxNLradius() const;
    double getMaxLradius() const;
    void updateListIons();

    // compute boundary for box containing all atoms to be known on local
    // processor that is values for list_boundary_left_, list_boundary_right_
    void setupListIonsBoundaries(const double rmax);

    void augmentIonsData(const int nsteps, const int dir, const int disp,
        const int locSize, const int maxLocSize, std::vector<IonData>& data,
        int* offset);
    bool inListIons(const double x, const double y, const double z);
    bool inLocalIons(const double x, const double y, const double z);

    /* local arrays for controlling ionic positions */
    std::vector<std::string> local_names_; // local ion names
    std::vector<double> taum_; // previous ion positions
    std::vector<double> tau0_; // current ion positions
    std::vector<double> taup_; // next ion positions
    std::vector<double> fion_; // ionic forces
    std::vector<double> velocity_; // ionic forces
    std::vector<double> pmass_;
    std::vector<short> atmove_;
    std::vector<unsigned short> rand_states_;
    std::vector<IonData> ions_data_;

    std::vector<int> gids_;

    // arrays with extra (beyond local) atomic data needed for enforcing
    // constraints.
    // made from list of interacting ions
    std::vector<double*> interacting_tau0_;
    std::vector<double*> interacting_taup_;
    std::vector<double*> interacting_fion_;
    std::vector<std::string> interacting_names_;
    std::vector<double> interacting_pmass_;
    std::vector<short> interacting_atmove_;
    // extra data involved in constraints, but not local
    std::vector<double> tau0_dummy_;
    std::vector<double> taup_dummy_;
    std::vector<double> fion_dummy_;
    std::vector<std::string> names_dummy_;
    std::vector<double> pmass_dummy_;
    std::vector<short> atmove_dummy_;

    int lstep_[3]; // number of steps to the left for each dimension to gather
                   // data
    int rstep_[3]; // number of steps to the right for each dimension to gather
                   // data

    bool setup_;

    void gatherLockedData(std::vector<int>& locked_data, const int root) const;
    void computeNumIons(void);
    int readNatoms(const std::string& input_file, const bool cell_relative);
    int readNatoms(std::ifstream* tfile, const bool cell_relative);
    int readAtomsFromXYZ(const std::string& filename, const bool cell_relative);
    void setupContraintsData(std::vector<Ion*>&);
    void clearStepperData();
    void initStepperData();
    void computeMaxNumProjs();
    void augmentIonicData(DistributedIonicData&, const int nsteps,
        const int dir, const int disp, const int locSize, const int maxLocSize,
        DistributedIonicData& data);

    void gatherNames(std::map<int, std::string>& names, const int root,
        const MPI_Comm comm) const;
    void gatherPositions(
        std::vector<double>& positions, const int root = 0) const;
    void gatherForces(std::vector<double>& forces, const int root = 0) const;

    void gatherNames(std::vector<std::string>& names, const int root,
        const MPI_Comm comm) const;
    void gatherPositions(std::vector<double>& positions, const int root,
        const MPI_Comm comm) const;
    void gatherLockedNames(std::vector<std::string>& names, const int root,
        const MPI_Comm comm) const;
    void gatherIndexes(
        std::vector<int>& indexes, const int root, const MPI_Comm comm) const;
    void gatherNLprojIds(
        std::vector<int>& nlprojids, const int root, const MPI_Comm comm) const;
    void gatherVelocities(std::vector<double>& velocities, const int root,
        const MPI_Comm comm) const;
    void gatherAtomicNumbers(
        std::vector<int>& atnumbers, const int root, const MPI_Comm comm) const;
    void gatherRandStates(std::vector<unsigned short>& rstates, const int root,
        const MPI_Comm comm) const;
    void gatherForces(
        std::vector<double>& forces, const int root, const MPI_Comm comm) const;
    bool hasLockedAtoms() const;

public:
    Ions(const double lat[3], const std::vector<Species>& sp);

    Ions(const Ions&, const double shift[3]);

    ~Ions();

    void setup();

    std::vector<std::string>& getLocalNames() { return local_names_; }
    std::vector<double>& getTau0() { return tau0_; }
    std::vector<double>& getTaup() { return taup_; }
    std::vector<double>& getTaum() { return taum_; }
    std::vector<double>& getVelocities() { return velocity_; }
    std::vector<double>& getFion() { return fion_; }
    std::vector<double>& getPmass() { return pmass_; }
    std::vector<short>& getAtmove() { return atmove_; }
    std::vector<unsigned short>& getRandStates() { return rand_states_; }
    std::vector<int>& getGids() { return gids_; }
    void resetForces()
    {
        std::vector<Ion*>::iterator ion = local_ions_.begin();
        while (ion != local_ions_.end())
        {
            (*ion)->resetForce();
            ion++;
        }
    }
    void removeMassCenterMotion();

    bool hasNLprojectors()
    {
        assert(max_num_proj_ >= 0);
        return (max_num_proj_ > 0);
    }

    double energySelf() const;
    double energyDiff(const short bc[3]) const;

    void printPositions(std::ostream& os, const int root = 0) const;
    void printPositionsLocal(std::ostream& os, const int root = 0) const;
    void printPositionsGlobal(std::ostream& os, const int root = 0) const;
    void printForces(std::ostream& os, const int root = 0) const;
    void printForcesLocal(std::ostream& os, const int root = 0) const;
    void printForcesGlobal(std::ostream& os, const int root = 0) const;
    int getNumIons(void);
    int getNumLocIons(void) const { return local_ions_.size(); }
    int getNumListIons(void) const { return list_ions_.size(); }
    std::vector<Ion*>& local_ions() { return local_ions_; }

    const std::vector<Ion*>& local_ions() const { return local_ions_; }
    const std::vector<Ion*>& list_ions() const { return list_ions_; }
    const std::vector<Ion*>& overlappingNL_ions() const
    {
        return overlappingNL_ions_;
    }
    const std::vector<Ion*>& overlappingVL_ions() const
    {
        return overlappingVL_ions_;
    }

    double computeIonicCharge() const;

    void iiforce(const short periodic[3]);
    double kinetic_E(void) const;

    void writePositions(HDFrestart& h5f_file);
    void writeVelocities(HDFrestart& h5f_file);
    void writeRandomStates(HDFrestart& h5f_file);
    void writeForces(HDFrestart& h5f_file);
    void writeAtomicIDs(HDFrestart& h5f_file);
    void writeAtomicNLprojIDs(HDFrestart& h5f_file);
    void writeAtomicNumbers(HDFrestart& h5f_file);
    void writeAtomNames(HDFrestart& h5f_file);
    void readLockedAtomNames(HDFrestart& h5f_file);
    void writeLockedAtomNames(HDFrestart& h5f_file);

    int countProjectors() const;
    int countProjectorsSubdomain() const;
    int countProjectorsHere() const;
    int countIonsHere() const;

    Ion* findIon(const std::string& name) const;
    Ion* findIon(const int index) const;
    Ion* findLocalIon(const int index) const;
    // check if ion is in list of local ions
    bool isLocal(const std::string& ion_name) const;

    void setPositions(const std::vector<double>& tau);
    void get_positions(std::vector<std::vector<double>>& r) const;
    void set_positions(const std::vector<std::vector<double>>& r);
    void get_forces(std::vector<std::vector<double>>& f) const;
    void set_forces(const std::vector<std::vector<double>>& f);

    void lockAtom(const std::string& name);

    // functions to initialize atomic coordinates
    // either from input file, or restart file
    int readAtoms(const std::string& filename, const bool cell_relative);
    int readAtoms(std::ifstream* tfile, const bool cell_relative);
    void initFromRestartFile(HDFrestart& h5_file);

    int getNValenceElectrons() const;
    void syncForces();
    void setVelocities(
        const std::vector<double>&, const std::vector<double>&, const double);
    void setTau0();
    void setPositionsToTau0();

    void getPositions(std::vector<double>& tau) const;
    void getVelocities(std::vector<double>& tau) const;
    void getForces(std::vector<double>& tau) const;
    void syncData(const std::vector<Species>& sp);
    // void syncNames(const int nions, std::vector<std::string>& local_names,
    // std::vector<std::string>& names);
    double getMaxVlRadius() const { return max_Vl_radius_; }
    double getMaxVnlRadius() const
    {
        assert(max_Vnl_radius_ >= 0.);
        return max_Vnl_radius_;
    }
    double getMaxListRadius() const;
    void updateIons();
    void shiftIons(const Vector3D& shift);

    const std::vector<double>& getMassesInteractingIons() const
    {
        return interacting_pmass_;
    }
    const std::vector<std::string>& getNamesInteractingIons() const
    {
        return interacting_names_;
    }
    std::vector<double*>& getPositionsInteractingIons()
    {
        return interacting_tau0_;
    }
    std::vector<double*>& getTaupInteractingIons() { return interacting_taup_; }
    std::vector<double*>& getForcesInteractingIons()
    {
        return interacting_fion_;
    }

    const std::vector<short>& getMovesInteractingIons() const
    {
        return interacting_atmove_;
    }

    void addIon2Lists(const IonData& data);

    const Species& getSpecies(short atomic_num) const
    {
        std::vector<Species>::const_iterator spi = species_.begin();
        while (spi != species_.end())
        {
            if (atomic_num == spi->getAtomicNumber())
            {
                break;
            }
            spi++;
        }

        return *spi;
    }

    void updateForcesInteractingIons();
    void updateTaupInteractingIons();
    void rescaleVelocities(const double factor);

    /*!
     * Calculate minimum distance between local pairs
     */
    double computeMinLocalSpacing() const;

    // void checkUnicityLocalIons();
};

#endif
