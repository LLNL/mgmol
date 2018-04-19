// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id:$
#ifndef MGMOL_H
#define MGMOL_H

#include "GridFuncVector.h"
#include "Lap.h"

#include "MGmol_prototypes.h"
#include "Control.h"
#include <iostream>

//inline double one(const double r){ return 1.; }

inline double linear(const double r){ return 1.-r; }

//class Rho;
class XConGrid;
class Electrostatic;
class ConstraintSet;
class Energy;
class MLWFTransform;
class NOLMOTransform;
class SpreadsAndCenters;
class OrbitalsTransform;
class Hamiltonian;
class ProjectedMatrices;
class ProjectedMatricesInterface;
class KBPsiMatrix;
class KBPsiMatrixInterface;
class KBPsiMatrixSparse;
class MasksSet;
class IonicAlgorithm;
class AOMMprojector;
class DMStrategy;
class OrbitalsPreconditioning;
class SpreadPenaltyInterface;
//class LocGridOrbitals;

#include "Ions.h"
#include "SparseDistMatrix.h"
#include "RemoteTasksDistMatrix.h"
#include "BasicDataDistributors.h"
#include "Rho.h"
#include "LocGridOrbitals.h"
#include "Forces.h"
#include "ClusterOrbitals.h"
#include "DistMatrixWithSparseComponent.h"

class MGmol
{
private:
    std::ostream& os_;
    
    MPI_Comm comm_;

    dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>* remote_tasks_DistMatrix_;
    static dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>* remote_tasks_DistMatrix_ptr_;
    
    // Compensating charge density
    RHODTYPE* rhoc_;
    
    XConGrid* xcongrid_;

    LocGridOrbitals* current_orbitals_;
    
    AOMMprojector* aomm_;
    
    Ions* ions_;

    Rho* rho_;

    Energy* energy_;

    Hamiltonian* hamiltonian_;
    
    Forces* forces_;

    MasksSet* currentMasks_;
    MasksSet* corrMasks_;
        
    //ProjectedMatrices* proj_matrices_;
    ProjectedMatricesInterface* proj_matrices_;
    
    IonicAlgorithm* geom_optimizer_;
    
    LocalizationRegions* lrs_;
    
    ClusterOrbitals* local_cluster_;
    
    KBPsiMatrixInterface* g_kbpsi_;
        
    SpreadsAndCenters* spreadf_;
    
    SpreadPenaltyInterface* spread_penalty_;
    
    DMStrategy* dm_strategy_;

    HDFrestart* h5f_file_;
        
    double background_charge_;
    double charge_in_cell_;
    double ionic_charge_;
    double total_energy_;
    ConstraintSet* constraints_;
    
    float md_time_;
    int md_iteration_;

    // private functions
    void check_anisotropy();
    void initBackground();
    double get_charge(RHODTYPE *rho);
    void initVcomp(Ions& ions);
    void printTimers();
    void computeLocKBPsiIonProjL(LocGridOrbitals& orbitals,
                       Ion& ion, const int ion_index,
                       dist_matrix::SparseDistMatrix<DISTMATDTYPE>*** prjsum,
                       const int l);
    void computeLocKBPsiIon(LocGridOrbitals& orbitals,
                 Ion& ion, const int ion_index, 
                 dist_matrix::SparseDistMatrix<DISTMATDTYPE>*** loc_kbpsi);
    void computeLocKBPsi(LocGridOrbitals& orbitals,
                  vector<Ion*>& ions_nl, 
                  dist_matrix::SparseDistMatrix<DISTMATDTYPE>*** prjsum);
    int read_rho_and_pot_hdf5(HDFrestart& file, 
                   Rho& rho);
    int read_restart_lrs(HDFrestart& h5f_file, 
              const string& dset_name);
    int read_restart_data(HDFrestart& h5f_file, 
               Rho& rho, 
               LocGridOrbitals& orbitals);
    void write_header();
    void getKBPsiAndHij(LocGridOrbitals& orbitals_i,
                    LocGridOrbitals& orbitals_j,
                    Ions& ions,
                    KBPsiMatrix*kbpsi_i,
                    KBPsiMatrix*kbpsi_j,
                    ProjectedMatricesInterface* projmatrices,
                    dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void getKBPsiAndHij(LocGridOrbitals& orbitals_i,
                    LocGridOrbitals& orbitals_j,
                    Ions& ions,
                    KBPsiMatrixInterface* kbpsi,
                    ProjectedMatricesInterface* projmatrices,
                    dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sh);
    void computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixInterface* const kbpsi,
            dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sparseH);
    void getKBPsiAndHij(LocGridOrbitals& orbitals_i,
                    LocGridOrbitals& orbitals_j,
                    Ions& ions,
                    KBPsiMatrixInterface* kbpsi,
                    ProjectedMatricesInterface* projmatrices);
    void getKBPsiAndHij(LocGridOrbitals& orbitals_i,
                    LocGridOrbitals& orbitals_j,
                    Ions& ions,
                    KBPsiMatrixInterface* kbpsi,
                    ProjectedMatricesInterface* projmatrices,
                    dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void getKBPsiAndHij(LocGridOrbitals& orbitals,
                    Ions& ions,
                    KBPsiMatrixInterface* kbpsi,
                    dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void computeHnlPhiAndAdd2HPhi(Ions& ions,
                LocGridOrbitals& phi,
                LocGridOrbitals& hphi,
                const KBPsiMatrixInterface* const kbpsi);                           
    void addHlocalij(LocGridOrbitals& orbitalsi,
                     LocGridOrbitals& orbitalsj);
    void addHlocal2matrix(LocGridOrbitals& orbitalsi,
                     LocGridOrbitals& orbitalsj, 
                     dist_matrix::SparseDistMatrix<DISTMATDTYPE>& sparseH);
    int dumprestartFile(LocGridOrbitals** orbitals, Ions& ions, Rho& rho,
                     const bool write_extrapolated_wf, const short count);

    void swapColumnsVect(dist_matrix::DistMatrix<DISTMATDTYPE>& evect,
                         const dist_matrix::DistMatrix<DISTMATDTYPE>& hb2N,
                         const std::vector<double>& eval,
                         dist_matrix::DistMatrix<DISTMATDTYPE>& work);
    void wftransform(LocGridOrbitals*,LocGridOrbitals*,Ions&);
    int readLRsFromInput(ifstream* tfile);
    void preWFextrapolation();
    void postWFextrapolation(LocGridOrbitals* orbitals);
    void computeResidualUsingHPhi(LocGridOrbitals& psi,
             const LocGridOrbitals& hphi,
             LocGridOrbitals& res,
             const bool applyB=true);
    
    int initial();

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
    BasicDataDistributors *data_distributor_;    

    OrbitalsPreconditioning* orbitals_precond_;

public:
    Electrostatic* electrostat_;
    
    MGmol(MPI_Comm comm, std::ostream& os);
    
    ~MGmol();
    
    void run();
    void initNuc(Ions& ions);
    void initKBR();
        
    dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>* getRemoteTasksDistMatrix()
    {
        assert( remote_tasks_DistMatrix_!=0 );
        return remote_tasks_DistMatrix_;
    }

    static dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>** getRemoteTasksDistMatrixPtr()
    {
        return &remote_tasks_DistMatrix_ptr_;
    }
    
    void global_exit(int i);
    void printEigAndOcc();
    
    int readInput(const string input_file);
    int readInput(const string input_file1,const string input_file2);

    int readParameters(ifstream* tfile,
                       bool& cell_relative);
    int readCoordinates(ifstream* tfile,
                        const bool cell_relative);
    int readCoordinates(const string filename,
                        const bool cell_relative);
    double computeConstraintResidual(LocGridOrbitals& orbitals, 
                         const LocGridOrbitals& hphi, 
                         LocGridOrbitals& res, 
                         const bool print_residual, 
                         const bool norm_res);
    
    void md(LocGridOrbitals** orbitals, Ions& ions);
    void lbfgsrlx(LocGridOrbitals** orbitals, Ions& ions);
    void computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixSparse* const kbpsi_i,
            const KBPsiMatrixSparse* const kbpsi_j,
            VariableSizeMatrix<sparserow>& mat, const bool consolidate);
    void computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixInterface* const kbpsi,
            VariableSizeMatrix<sparserow>& mat, const bool consolidate);
    void computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixSparse* const kbpsi_i,
            const KBPsiMatrixSparse* const kbpsi_j,
            dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrix* const kbpsi_i,
            const KBPsiMatrix* const kbpsi_j,
            dist_matrix::DistMatrix<DISTMATDTYPE>& hij);            
    void computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixInterface* const kbpsi,
            dist_matrix::DistMatrix<DISTMATDTYPE>& hij);
    void computeHij(LocGridOrbitals& orbitals_i,
            LocGridOrbitals& orbitals_j,
            const Ions& ions,
            const KBPsiMatrixInterface* const kbpsi,
            ProjectedMatricesInterface*);
    void addHlocal2matrix(LocGridOrbitals& orbitalsi,
                 LocGridOrbitals& orbitalsj, 
                 dist_matrix::DistMatrixWithSparseComponent<DISTMATDTYPE>& mat);
    void addHlocal2matrix(LocGridOrbitals& orbitalsi,
                     LocGridOrbitals& orbitalsj, 
                     VariableSizeMatrix<sparserow>& mat);    
    void update_pot(const pb::GridFunc<POTDTYPE>& vh_init,
                const Ions& ions);
    void update_pot(const Ions& ions);
    int mvp(LocGridOrbitals& orbitals,
        Ions& ions,
        Potentials& pot,
        const double etol);
    int quench(LocGridOrbitals* orbitals,
           Ions& ions, 
           const int max_steps,
           const int iprint,
           double& last_eks);
    void runfire(LocGridOrbitals** orbitals, Ions& ions);
    void moveVnuc(Ions& ions);
    void resetProjectedMatricesAndDM(LocGridOrbitals& orbitals, Ions& ions);
    int getMLWF(MLWFTransform& mlwft,
            LocGridOrbitals& orbitals, 
            LocGridOrbitals& work_orbitals,
            const double dd,
            const bool apply_flag);
    bool rotateStatesPairsCommonCenter(LocGridOrbitals& orbitals,
                                   LocGridOrbitals& work_orbitals);
    bool rotateStatesPairsOverlap(LocGridOrbitals& orbitals,
                              LocGridOrbitals& work_orbitals,
                              const double);
    void disentangleOrbitals(LocGridOrbitals& orbitals,
                             LocGridOrbitals& work_orbitals,
                             Ions&, int&);
    void updateHmatrix(LocGridOrbitals& orbitals, Ions& ions);
    void getHpsiAndTheta(Ions& ions, 
                      LocGridOrbitals& phi, 
                      LocGridOrbitals& hphi,
                      const KBPsiMatrixInterface* const kbpsi);
    void getHpsiAndTheta(Ions& ions, 
                      LocGridOrbitals& phi, 
                      LocGridOrbitals& hphi);
    double computePrecondResidual(LocGridOrbitals& phi,
                        LocGridOrbitals& hphi,
                        LocGridOrbitals& res,
                        Ions& ions,
                        KBPsiMatrixInterface* kbpsi, 
                        const bool print_residual,
                        const bool norm_res);
    void addResidualSpreadPenalty(LocGridOrbitals& phi,
                                  LocGridOrbitals& res);
    int get_NOLMO(NOLMOTransform& noot,
              LocGridOrbitals& orbitals, 
              LocGridOrbitals& work_orbitals,
              const double dd,
              const bool apply_flag);
    int read_params1(ifstream* tfile);
    void adaptLR(const SpreadsAndCenters* spreadf,
             const OrbitalsTransform* ot);
    int update_masks();
    void move_orbitals(LocGridOrbitals** orbitals);
    int getMLWF2states(const int st1, const int st2,
            LocGridOrbitals& orbitals, 
            LocGridOrbitals& work_orbitals);
    void extrapolate_centers(bool small_move);
    void compute_centers(bool small_move, LocGridOrbitals** orbitals);
    void extrapolate_orbitals(LocGridOrbitals** orbitals);
    LocGridOrbitals* new_orbitals_with_current_LRs(bool setup=true);
    void update_orbitals_LRs(LocGridOrbitals** orbitals);
    void clearOldOrbitals();
    void getKBPsiAndHij(LocGridOrbitals& orbitals,
                        Ions& ions);
    int write_hdf5(const string filename, 
        vector< vector<RHODTYPE> >& rho, 
        Ions& ions,
        LocGridOrbitals& orbitals,
        LocalizationRegions& lrs);
    int write_hdf5(HDFrestart& h5f_file, 
         vector< vector<RHODTYPE> >& rho, 
         Ions& ions,
         LocGridOrbitals& orbitals,
         LocalizationRegions& lrs);
    double get_evnl(const Ions& ions, LocGridOrbitals& orbitals);
    void sebprintPositions();
    void sebprintForces();
    void get_positions(vector<vector<double> > &r);
    void set_positions(vector<vector<double> > &r);
    void get_forces(vector<vector<double> > &f);
    void set_forces(vector<vector<double> > &f);
    int nions(){ return ions_->getNumIons();}
    double getTotalEnergy();
    void setup();
    int setupFromInput(const string input_file);
    int setupLRsFromInput(const string input_file);
    int setupConstraintsFromInput(const string input_file);
    void cleanup();
    void geomOptimSetup();
    void geomOptimQuench();
    void geomOptimComputeForces();
    int geomOptimRun1Step();
    void geomOptimDumpRestart();
    void geomOptimSetForces(const vector<vector<double> >& f);
    short geomOptimCheckTolForces(const double tol_force);

    void finalEnergy();
    void printMM();

    void projectOutKernel(LocGridOrbitals& phi);

    void precond_mg(LocGridOrbitals& orbitals);
    void setGamma(const pb::Lap<ORBDTYPE>& lapOper, const Potentials& pot);
    double computeResidual(LocGridOrbitals& orbitals,
            LocGridOrbitals& work_orbitals,
            LocGridOrbitals& res,
            const bool print_residual,
            const bool norm_res);

    void force(LocGridOrbitals& orbitals, Ions& ions)
    {
        forces_->force(orbitals, ions);
    }
};

#endif
