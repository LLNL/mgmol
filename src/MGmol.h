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
class ProjectedMatrices;
class ProjectedMatricesInterface;
class KBPsiMatrix;
class KBPsiMatrixSparse;
class MasksSet;
class DMStrategy;

template <class T>
class IonicAlgorithm;

#include "AOMMprojector.h"
#include "BasicDataDistributors.h"
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

template <class T>
class MGmol : public MGmolInterface
{
private:
    std::ostream& os_;

    MPI_Comm comm_;

    XConGrid* xcongrid_;

    T* current_orbitals_;

    AOMMprojector* aomm_;

    Ions* ions_;

    Rho<T>* rho_;

    Energy<T>* energy_;

    Hamiltonian<T>* hamiltonian_;

    Forces<T>* forces_;

    MasksSet* currentMasks_;
    MasksSet* corrMasks_;

    // ProjectedMatrices* proj_matrices_;
    ProjectedMatricesInterface* proj_matrices_;

    IonicAlgorithm<T>* geom_optimizer_;

    LocalizationRegions* lrs_;

    ClusterOrbitals* local_cluster_;

    KBPsiMatrixSparse* g_kbpsi_;

    SpreadsAndCenters<T>* spreadf_;

    SpreadPenaltyInterface<T>* spread_penalty_;

    DMStrategy* dm_strategy_;

    HDFrestart* h5f_file_;

    double total_energy_;
    ConstraintSet* constraints_;

    OrbitalsExtrapolation<T>* orbitals_extrapol_;

    float md_time_;
    int md_iteration_;

    // private functions
    void check_anisotropy();
    double get_charge(RHODTYPE* rho);
    void printTimers();
    int read_rho_and_pot_hdf5(HDFrestart& file, Rho<T>& rho);
    int read_restart_lrs(HDFrestart& h5f_file, const std::string& dset_name);
    int read_restart_data(HDFrestart& h5f_file, Rho<T>& rho, T& orbitals);
    void write_header();
    void getKBPsiAndHij(T& orbitals_i, T& orbitals_j, Ions& ions,
        KBPsiMatrixSparse* kbpsi, ProjectedMatricesInterface* projmatrices);
    void getKBPsiAndHij(T& orbitals_i, T& orbitals_j, Ions& ions,
        KBPsiMatrixSparse* kbpsi, ProjectedMatricesInterface* projmatrices,
        dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void getKBPsiAndHij(T& orbitals, Ions& ions, KBPsiMatrixSparse* kbpsi,
        dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void computeHnlPhiAndAdd2HPhi(
        Ions& ions, T& phi, T& hphi, const KBPsiMatrixSparse* const kbpsi);
    void addHlocalij(
        T& orbitalsi, T& orbitalsj, ProjectedMatricesInterface* projmatrices);
    int dumprestartFile(T** orbitals, Ions& ions, Rho<T>& rho,
        const bool write_extrapolated_wf, const short count);

    void swapColumnsVect(dist_matrix::DistMatrix<DISTMATDTYPE>& evect,
        const dist_matrix::DistMatrix<DISTMATDTYPE>& hb2N,
        const std::vector<double>& eval,
        dist_matrix::DistMatrix<DISTMATDTYPE>& work);
    void wftransform(T*, T*, Ions&);
    int readLRsFromInput(std::ifstream* tfile);
    void preWFextrapolation();
    void postWFextrapolation(T* orbitals);
    void computeResidualUsingHPhi(
        T& psi, const T& hphi, T& res, const bool applyB = true);

    template <typename MemorySpaceType>
    int initial();
    void initialMasks();

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

    /* Data distribution objects */
    BasicDataDistributors* data_distributor_;

    OrbitalsPreconditioning<T>* orbitals_precond_;

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
    double computeConstraintResidual(T& orbitals, const T& hphi, T& res,
        const bool print_residual, const bool norm_res);

    void md(T** orbitals, Ions& ions);
    void lbfgsrlx(T** orbitals, Ions& ions);

    template <class T2>
    void computeHij(T& orbitals_i, T& orbitals_j, const Ions& ions,
        const KBPsiMatrixSparse* const kbpsi_i,
        const KBPsiMatrixSparse* const kbpsi_j, T2& mat,
        const bool consolidate);

    void computeHij_private(T& orbitals_i, T& orbitals_j, const Ions& ions,
        const KBPsiMatrixSparse* const kbpsi_i,
        const KBPsiMatrixSparse* const kbpsi_j,
        dist_matrix::DistMatrix<DISTMATDTYPE>& mat);

    void computeHij(T& orbitals_i, T& orbitals_j, const Ions& ions,
        const KBPsiMatrixSparse* const kbpsi,
        dist_matrix::DistMatrix<double>& mat, const bool consolidate);

    void computeHij_private(T& orbitals_i, T& orbitals_j, const Ions& ions,
        const KBPsiMatrixSparse* const kbpsi_i,
        dist_matrix::DistMatrix<DISTMATDTYPE>& mat);

    void computeHij(T& orbitals_i, T& orbitals_j, const Ions& ions,
        const KBPsiMatrixSparse* const kbpsi,
        VariableSizeMatrix<sparserow>& mat, const bool consolidate);

    void computeHij(T& orbitals_i, T& orbitals_j, const Ions& ions,
        const KBPsiMatrixSparse* const kbpsi, ProjectedMatricesInterface*);

    void addHlocal2matrix(
        T& orbitalsi, T& orbitalsj, dist_matrix::DistMatrix<double>& mat);
    void addHlocal2matrix(
        T& orbitalsi, T& orbitalsj, VariableSizeMatrix<SparseRow>& mat);

    void update_pot(const pb::GridFunc<POTDTYPE>& vh_init, const Ions& ions);
    void update_pot(const Ions& ions);
    int quench(T* orbitals, Ions& ions, const int max_steps, const int iprint,
        double& last_eks);
    int outerSolve(T& orbitals, T& work_orbitals, Ions& ions,
        const int max_steps, const int iprint, double& last_eks);
    void runfire(T** orbitals, Ions& ions);
    void moveVnuc(Ions& ions);
    void resetProjectedMatricesAndDM(T& orbitals, Ions& ions);
    int getMLWF(MLWFTransform& mlwft, T& orbitals, T& work_orbitals,
        const double dd, const bool apply_flag);
    bool rotateStatesPairsCommonCenter(T& orbitals, T& work_orbitals);
    bool rotateStatesPairsOverlap(T& orbitals, T& work_orbitals, const double);
    void disentangleOrbitals(T& orbitals, T& work_orbitals, Ions&, int&);
    void updateHmatrix(T& orbitals, Ions& ions);
    void getHpsiAndTheta(
        Ions& ions, T& phi, T& hphi, const KBPsiMatrixSparse* const kbpsi);
    void getHpsiAndTheta(Ions& ions, T& phi, T& hphi);
    double computePrecondResidual(T& phi, T& hphi, T& res, Ions& ions,
        KBPsiMatrixSparse* kbpsi, const bool print_residual,
        const bool norm_res);
    void addResidualSpreadPenalty(T& phi, T& res);
    int get_NOLMO(NOLMOTransform& noot, T& orbitals, T& work_orbitals,
        const double dd, const bool apply_flag);
    void adaptLR(
        const SpreadsAndCenters<T>* spreadf, const OrbitalsTransform* ot);
    int update_masks();
    void move_orbitals(T** orbitals);
    int getMLWF2states(
        const int st1, const int st2, T& orbitals, T& work_orbitals);
    void extrapolate_centers(bool small_move);
    void compute_centers(bool small_move, T** orbitals);
    void extrapolate_orbitals(T** orbitals);
    T* new_orbitals_with_current_LRs(bool setup = true);
    void update_orbitals_LRs(T** orbitals);
    void clearOldOrbitals();
    void getKBPsiAndHij(T& orbitals, Ions& ions);
    int write_hdf5(const std::string& filename,
        std::vector<std::vector<RHODTYPE>>& rho, Ions& ions, T& orbitals);
    int write_hdf5(HDFrestart& h5f_file,
        std::vector<std::vector<RHODTYPE>>& rho, Ions& ions, T& orbitals,
        LocalizationRegions& lrs);
    double get_evnl(const Ions& ions, T& orbitals);
    void sebprintPositions();
    void sebprintForces();
    void get_positions(std::vector<std::vector<double>>& r);
    void set_positions(std::vector<std::vector<double>>& r);
    void get_forces(std::vector<std::vector<double>>& f);
    void set_forces(std::vector<std::vector<double>>& f);
    int nions() { return ions_->getNumIons(); }
    double getTotalEnergy();
    void setup();
    int setupFromInput(const std::string input_file) override;
    int setupLRsFromInput(const std::string input_file) override;
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

    void projectOutKernel(T& phi);

    void precond_mg(T& orbitals);
    void setGamma(const pb::Lap<ORBDTYPE>& lapOper, const Potentials& pot);
    double computeResidual(T& orbitals, T& work_orbitals, T& res,
        const bool print_residual, const bool norm_res);
    void applyAOMMprojection(T&);
    void force(T& orbitals, Ions& ions) { forces_->force(orbitals, ions); }
};
// Instantiate static variables here to avoid clang warnings
template <class T>
Timer MGmol<T>::adaptLR_tm_("MGmol::adaptLR");
template <class T>
Timer MGmol<T>::dump_tm_("MGmol::dump");
template <class T>
Timer MGmol<T>::total_tm_("MGmol::total");
template <class T>
Timer MGmol<T>::setup_tm_("MGmol::setup");
template <class T>
Timer MGmol<T>::closing_tm_("MGmol::closing");
template <class T>
Timer MGmol<T>::init_tm_("MGmol::init");
template <class T>
Timer MGmol<T>::evnl_tm_("MGmol::evnl");
template <class T>
Timer MGmol<T>::get_res_tm_("MGmol::comp_res_from_Hphi");
template <class T>
Timer MGmol<T>::computeHij_tm_("MGmol::computeHij");
template <class T>
Timer MGmol<T>::get_Hpsi_and_Hij_tm_("MGmol::get_Hpsi_and_Hij");
template <class T>
Timer MGmol<T>::comp_res_tm_("MGmol::comp_res");
template <class T>
Timer MGmol<T>::init_nuc_tm_("MGmol::init_nuc");
#endif
