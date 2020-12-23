// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_H
#define MGMOL_H

#include "Energy.h"
#include "GridFuncVector.h"
#include "Hamiltonian.h"
#include "Lap.h"

#include "Control.h"
#include "MGmol_prototypes.h"
#include <iostream>

// inline double one(const double r){ return 1.; }

inline double linear(const double r) { return 1. - r; }

class XConGrid;
class Electrostatic;
class ConstraintSet;
class MLWFTransform;
class NOLMOTransform;
class OrbitalsTransform;
template <class OrbitalsType>
class ProjectedMatrices;
class ProjectedMatricesInterface;
class KBPsiMatrix;
class KBPsiMatrixSparse;
class MasksSet;
class DMStrategy;

template <class OrbitalsType>
class IonicAlgorithm;

#include "AOMMprojector.h"
#include "ClusterOrbitals.h"
#include "ExtendedGridOrbitals.h"
#include "Forces.h"
#include "Ions.h"
#include "LocGridOrbitals.h"
#include "MGmolInterface.h"
#include "OrbitalsExtrapolation.h"
#include "OrbitalsPreconditioning.h"
#include "Rho.h"
#include "SpreadPenaltyInterface.h"
#include "SpreadsAndCenters.h"

template <class OrbitalsType>
class MGmol : public MGmolInterface
{
private:
    std::ostream& os_;

    MPI_Comm comm_;

    XConGrid* xcongrid_;

    OrbitalsType* current_orbitals_;

    AOMMprojector* aomm_;

    Ions* ions_;

    Rho<OrbitalsType>* rho_;

    Energy<OrbitalsType>* energy_;

    Hamiltonian<OrbitalsType>* hamiltonian_;

    Forces<OrbitalsType>* forces_;

    MasksSet* currentMasks_;
    MasksSet* corrMasks_;

    // ProjectedMatrices* proj_matrices_;
    ProjectedMatricesInterface* proj_matrices_;

    IonicAlgorithm<OrbitalsType>* geom_optimizer_;

    std::shared_ptr<LocalizationRegions> lrs_;

    ClusterOrbitals* local_cluster_;

    KBPsiMatrixSparse* g_kbpsi_;

    SpreadsAndCenters<OrbitalsType>* spreadf_;

    SpreadPenaltyInterface<OrbitalsType>* spread_penalty_;

    DMStrategy* dm_strategy_;

    HDFrestart* h5f_file_;

    double total_energy_;
    ConstraintSet* constraints_;

    OrbitalsExtrapolation<OrbitalsType>* orbitals_extrapol_;

    float md_time_;
    int md_iteration_;

    // private functions
    void check_anisotropy();
    double get_charge(RHODTYPE* rho);
    void printTimers();
    int read_rho_and_pot_hdf5(HDFrestart& file, Rho<OrbitalsType>& rho);
    int read_restart_lrs(HDFrestart& h5f_file, const std::string& dset_name);
    int read_restart_data(
        HDFrestart& h5f_file, Rho<OrbitalsType>& rho, OrbitalsType& orbitals);
    void write_header();
    void getKBPsiAndHij(OrbitalsType& orbitals_i, OrbitalsType& orbitals_j,
        Ions& ions, KBPsiMatrixSparse* kbpsi,
        ProjectedMatricesInterface* projmatrices);
    void getKBPsiAndHij(OrbitalsType& orbitals_i, OrbitalsType& orbitals_j,
        Ions& ions, KBPsiMatrixSparse* kbpsi,
        ProjectedMatricesInterface* projmatrices,
        dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void getKBPsiAndHij(OrbitalsType& orbitals, Ions& ions,
        KBPsiMatrixSparse* kbpsi, dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void computeHnlPhiAndAdd2HPhi(Ions& ions, OrbitalsType& phi,
        OrbitalsType& hphi, const KBPsiMatrixSparse* const kbpsi);
    int dumprestartFile(OrbitalsType** orbitals, Ions& ions,
        Rho<OrbitalsType>& rho, const bool write_extrapolated_wf,
        const short count);

    void swapColumnsVect(dist_matrix::DistMatrix<DISTMATDTYPE>& evect,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& hb2N,
        const std::vector<double>& eval,
        dist_matrix::DistMatrix<DISTMATDTYPE>& work);
    void wftransform(OrbitalsType*, OrbitalsType*, Ions&);
    int readLRsFromInput(std::ifstream* tfile);
    void preWFextrapolation();
    void postWFextrapolation(OrbitalsType* orbitals);
    void computeResidualUsingHPhi(OrbitalsType& psi, const OrbitalsType& hphi,
        OrbitalsType& res, const bool applyB = true);

    template <typename MemorySpaceType>
    int initial();
    void initialMasks();
    int setupLRsFromInput(const std::string input_file);

    // timers
    static Timer total_tm_;
    static Timer setup_tm_;
    static Timer closing_tm_;
    static Timer init_tm_;
    static Timer dump_tm_;
    static Timer adaptLR_tm_;
    static Timer evnl_tm_;
    static Timer get_res_tm_;
    static Timer computeHij_tm_;
    static Timer get_Hpsi_and_Hij_tm_;
    static Timer comp_res_tm_;
    static Timer init_nuc_tm_;

    OrbitalsPreconditioning<OrbitalsType>* orbitals_precond_;

public:
    Electrostatic* electrostat_;

    MGmol(MPI_Comm comm, std::ostream& os);

    ~MGmol() override;

    void run() override;
    void initNuc(Ions& ions);
    void initKBR();

    void global_exit(int i);
    void printEigAndOcc();

    int readCoordinates(std::ifstream* tfile, const bool cell_relative);
    int readCoordinates(const std::string& filename, const bool cell_relative);
    double computeConstraintResidual(OrbitalsType& orbitals,
        const OrbitalsType& hphi, OrbitalsType& res, const bool print_residual,
        const bool norm_res);

    void md(OrbitalsType** orbitals, Ions& ions);
    void lbfgsrlx(OrbitalsType** orbitals, Ions& ions);

    template <class MatrixType>
    void computeHij(OrbitalsType& orbitals_i, OrbitalsType& orbitals_j,
        const Ions& ions, const KBPsiMatrixSparse* const kbpsi_i,
        const KBPsiMatrixSparse* const kbpsi_j, MatrixType& mat,
        const bool consolidate);

    void computeHij_private(OrbitalsType& orbitals_i, OrbitalsType& orbitals_j,
        const Ions& ions, const KBPsiMatrixSparse* const kbpsi_i,
        const KBPsiMatrixSparse* const kbpsi_j,
        dist_matrix::DistMatrix<DISTMATDTYPE>& mat);

    void computeHij(OrbitalsType& orbitals_i, OrbitalsType& orbitals_j,
        const Ions& ions, const KBPsiMatrixSparse* const kbpsi,
        dist_matrix::DistMatrix<double>& mat, const bool consolidate);

    void computeHij_private(OrbitalsType& orbitals_i, OrbitalsType& orbitals_j,
        const Ions& ions, const KBPsiMatrixSparse* const kbpsi_i,
        dist_matrix::DistMatrix<DISTMATDTYPE>& mat);

    void computeHij(OrbitalsType& orbitals_i, OrbitalsType& orbitals_j,
        const Ions& ions, const KBPsiMatrixSparse* const kbpsi,
        VariableSizeMatrix<sparserow>& mat, const bool consolidate);

    void computeHij(OrbitalsType& orbitals_i, OrbitalsType& orbitals_j,
        const Ions& ions, const KBPsiMatrixSparse* const kbpsi,
        ProjectedMatricesInterface*);

    template <class MatrixType>
    void addHlocal2matrix(
        OrbitalsType& orbitalsi, OrbitalsType& orbitalsj, MatrixType& mat);
    void addHlocal2matrix(OrbitalsType& orbitalsi, OrbitalsType& orbitalsj,
        VariableSizeMatrix<SparseRow>& mat);

    void update_pot(const pb::GridFunc<POTDTYPE>& vh_init, const Ions& ions);
    void update_pot(const Ions& ions);
    int quench(OrbitalsType* orbitals, Ions& ions, const int max_steps,
        const int iprint, double& last_eks);
    int outerSolve(OrbitalsType& orbitals, OrbitalsType& work_orbitals,
        Ions& ions, const int max_steps, const int iprint, double& last_eks);
    void runfire(OrbitalsType** orbitals, Ions& ions);
    void moveVnuc(Ions& ions);
    void resetProjectedMatricesAndDM(OrbitalsType& orbitals, Ions& ions);
    int getMLWF(MLWFTransform& mlwft, OrbitalsType& orbitals,
        OrbitalsType& work_orbitals, const double dd, const bool apply_flag);
    bool rotateStatesPairsCommonCenter(
        OrbitalsType& orbitals, OrbitalsType& work_orbitals);
    bool rotateStatesPairsOverlap(
        OrbitalsType& orbitals, OrbitalsType& work_orbitals, const double);
    void disentangleOrbitals(
        OrbitalsType& orbitals, OrbitalsType& work_orbitals, Ions&, int&);
    void updateHmatrix(OrbitalsType& orbitals, Ions& ions);
    void getHpsiAndTheta(Ions& ions, OrbitalsType& phi, OrbitalsType& hphi,
        const KBPsiMatrixSparse* const kbpsi);
    void getHpsiAndTheta(Ions& ions, OrbitalsType& phi, OrbitalsType& hphi);
    double computePrecondResidual(OrbitalsType& phi, OrbitalsType& hphi,
        OrbitalsType& res, Ions& ions, KBPsiMatrixSparse* kbpsi,
        const bool print_residual, const bool norm_res);
    void addResidualSpreadPenalty(OrbitalsType& phi, OrbitalsType& res);
    int get_NOLMO(NOLMOTransform& noot, OrbitalsType& orbitals,
        OrbitalsType& work_orbitals, const double dd, const bool apply_flag);
    void adaptLR(const SpreadsAndCenters<OrbitalsType>* spreadf,
        const OrbitalsTransform* ot);
    int update_masks();
    void move_orbitals(OrbitalsType** orbitals);
    int getMLWF2states(const int st1, const int st2, OrbitalsType& orbitals,
        OrbitalsType& work_orbitals);
    void extrapolate_centers(bool small_move);
    void compute_centers(bool small_move, OrbitalsType** orbitals);
    void extrapolate_orbitals(OrbitalsType** orbitals);
    OrbitalsType* new_orbitals_with_current_LRs(bool setup = true);
    void update_orbitals_LRs(OrbitalsType** orbitals);
    void clearOldOrbitals();
    void getKBPsiAndHij(OrbitalsType& orbitals, Ions& ions);
    int write_hdf5(const std::string& filename,
        std::vector<std::vector<RHODTYPE>>& rho, Ions& ions,
        OrbitalsType& orbitals);
    int write_hdf5(HDFrestart& h5f_file,
        std::vector<std::vector<RHODTYPE>>& rho, Ions& ions,
        OrbitalsType& orbitals, std::shared_ptr<LocalizationRegions> lrs);
    double get_evnl(const Ions& ions);
    void sebprintPositions();
    void sebprintForces();
    void get_positions(std::vector<std::vector<double>>& r);
    void set_positions(std::vector<std::vector<double>>& r);
    void get_forces(std::vector<std::vector<double>>& f);
    void set_forces(std::vector<std::vector<double>>& f);
    int nions() { return ions_->getNumIons(); }
    double getTotalEnergy();
    void setup();
    int setupLRs(const std::string input_file) override;
    int setupFromInput(const std::string input_file) override;
    int setupConstraintsFromInput(const std::string input_file) override;
    void cleanup();
    void geomOptimSetup();
    void geomOptimQuench();
    void geomOptimComputeForces();
    int geomOptimRun1Step();
    void geomOptimDumpRestart();
    void geomOptimSetForces(const std::vector<std::vector<double>>& f);
    short geomOptimCheckTolForces(const double tol_force);

    void finalEnergy();
    void printMM();

    void projectOutKernel(OrbitalsType& phi);

    void precond_mg(OrbitalsType& orbitals);
    void setGamma(const pb::Lap<ORBDTYPE>& lapOper, const Potentials& pot);
    double computeResidual(OrbitalsType& orbitals, OrbitalsType& work_orbitals,
        OrbitalsType& res, const bool print_residual, const bool norm_res);
    void applyAOMMprojection(OrbitalsType&);
    void force(OrbitalsType& orbitals, Ions& ions)
    {
        forces_->force(orbitals, ions);
    }
};
// Instantiate static variables here to avoid clang warnings
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::adaptLR_tm_("MGmol::adaptLR");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::dump_tm_("MGmol::dump");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::total_tm_("MGmol::total");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::setup_tm_("MGmol::setup");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::closing_tm_("MGmol::closing");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::init_tm_("MGmol::init");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::evnl_tm_("MGmol::evnl");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::get_res_tm_("MGmol::comp_res_from_Hphi");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::computeHij_tm_("MGmol::computeHij");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::get_Hpsi_and_Hij_tm_("MGmol::get_Hpsi_and_Hij");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::comp_res_tm_("MGmol::comp_res");
template <class OrbitalsType>
Timer MGmol<OrbitalsType>::init_nuc_tm_("MGmol::init_nuc");
#endif
