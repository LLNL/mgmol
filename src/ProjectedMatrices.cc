// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ProjectedMatrices.h"

#ifdef PROCRUSTES
#include "munkres.h"
void procrustes(dist_matrix::DistMatrix<DISTMATDTYPE>& a,
    dist_matrix::DistMatrix<DISTMATDTYPE>& b,
    dist_matrix::DistMatrix<DISTMATDTYPE>& p);
#endif

#include "Control.h"
#include "DensityMatrix.h"
#include "GramMatrix.h"
#include "HDFrestart.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "ReplicatedWorkSpace.h"
#include "SparseDistMatrix.h"
#include "SubMatrices.h"
#include "fermi.h"
#include "tools.h"

#include <fstream>
#include <iomanip>

#define RY2EV 13.605804

Timer ProjectedMatrices::sygv_tm_("ProjectedMatrices::sygv");
Timer ProjectedMatrices::compute_inverse_tm_(
    "ProjectedMatrices::computeInverse");
Timer ProjectedMatrices::compute_invB_tm_("ProjectedMatrices::computeInvB");
Timer ProjectedMatrices::init_gram_matrix_tm_(
    "ProjectedMatrices::initialize_Gram_Matrix");
Timer ProjectedMatrices::update_theta_tm_("ProjectedMatrices::updateTheta");
Timer ProjectedMatrices::update_submatT_tm_("ProjectedMatrices::updateSubmatT");
Timer ProjectedMatrices::update_submatX_tm_("ProjectedMatrices::updateSubmatX");
Timer ProjectedMatrices::eigsum_tm_("ProjectedMatrices::eigsum");
Timer ProjectedMatrices::consolidate_H_tm_("ProjectedMatrices::consolidate_sH");

// const short LocGridOrbitals::sparse_distmatrix_nb_tasks_per_partitions_=96;
// const short LocGridOrbitals::sparse_distmatrix_nb_tasks_per_partitions_=256;
short ProjectedMatrices::sparse_distmatrix_nb_tasks_per_partitions_ = 256;
short ProjectedMatrices::n_instances_                               = 0;

dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>*
    ProjectedMatrices::remote_tasks_DistMatrix_
    = 0;
GramMatrix* ProjectedMatrices::gram_4dotProducts_  = 0;
DensityMatrix* ProjectedMatrices::dm_4dot_product_ = 0;

static int sparse_distmatrix_nb_partitions = 128;

ProjectedMatrices::ProjectedMatrices(const int ndim, const bool with_spin)
    : with_spin_(with_spin)
{
    dim_             = ndim;
    width_           = 0.;
    min_val_         = 0.25;
    submat_indexing_ = 0;

    dm_        = new DensityMatrix(ndim);
    gm_        = new GramMatrix(ndim);
    mat_X_old_ = 0;
    mat_L_old_ = 0;

#ifdef USE_DIS_MAT
    submatWork_ = 0;
    submatLS_   = 0;
#endif
    localX_ = 0;
    localT_ = 0;

    sH_ = 0;

    if (dim_ > 0)
    {
        eigenvalues_.resize(dim_);
    }

    if (dim_ > 0)
    {
        matH_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("H", ndim, ndim);

        matHB_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("HB", ndim, ndim);
        theta_ = new dist_matrix::DistMatrix<DISTMATDTYPE>("Theta", ndim, ndim);
        u_     = new dist_matrix::DistMatrix<DISTMATDTYPE>("U", ndim, ndim);
        work_  = new dist_matrix::DistMatrix<DISTMATDTYPE>("work", ndim, ndim);
#ifdef PROCRUSTES
        occupation0_.resize(dim_);
        u0_ = 0;
        p_  = new dist_matrix::DistMatrix<DISTMATDTYPE>("P", ndim, ndim);

        u0_->identity();
#endif
    }

    if (onpe0)
        (*MPIdata::sout)
            << "ProjectedMatrices: sparse_distmatrix_nb_tasks_per_partitions_="
            << sparse_distmatrix_nb_tasks_per_partitions_ << endl;
    n_instances_++;
}

ProjectedMatrices::~ProjectedMatrices()
{
    assert(dm_ != 0);
    assert(gm_ != 0);

    delete dm_;
    dm_ = 0;
    delete gm_;
    gm_ = 0;
    if (mat_X_old_)
    {
        delete mat_X_old_;
        mat_X_old_ = 0;
    }
    if (mat_L_old_)
    {
        delete mat_L_old_;
        mat_L_old_ = 0;
    }

    if (dim_ > 0)
    {
        assert(matH_ != 0);
        assert(theta_ != 0);
        assert(u_ != 0);
        assert(work_ != 0);

        delete matHB_;
        matHB_ = 0;
        delete matH_;
        delete theta_;
        theta_ = 0;
        delete u_;
        delete work_;

        delete localX_;
        localX_ = 0;
        delete localT_;
        localT_ = 0;

        delete sH_;
        sH_ = 0;

        delete submatWork_;
        submatWork_ = 0;
        delete submatLS_;
        submatLS_ = 0;

#ifdef PROCRUSTES
        delete u0_;
        delete p_;
#endif
        delete submat_indexing_;
        submat_indexing_ = 0;
    }

    if (n_instances_ == 1)
    {
        if (gram_4dotProducts_ != 0)
        {
            delete gram_4dotProducts_;
            gram_4dotProducts_ = 0;
        }
    }

    n_instances_--;
};

void ProjectedMatrices::setup(
    const double kbt, const int nel, const vector<vector<int>>& global_indexes)
{
    assert(global_indexes.size() > 0);

    ProjectedMatricesInterface::setup(kbt, nel, global_indexes);

    global_indexes_ = global_indexes;

#ifdef USE_DIS_MAT
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSpin();
    if (submat_indexing_ != 0) delete submat_indexing_;
    submat_indexing_ = new dist_matrix::SubMatricesIndexing<DISTMATDTYPE>(
        global_indexes, comm, gm_->getMatrix());

    if (localX_ != NULL) delete localX_;
    if (localT_ != NULL) delete localT_;

    if (sH_ != NULL) delete sH_;

    if (submatWork_ != NULL) delete submatWork_;
    if (submatLS_ != NULL) delete submatLS_;
    // if( orbitals_type_!=0 )
    {

        if (dim_ > 0)
        {
            submatWork_ = new dist_matrix::SubMatrices<DISTMATDTYPE>("Work",
                global_indexes, comm, dm_->getMatrix(), *submat_indexing_);
            submatLS_   = new dist_matrix::SubMatrices<DISTMATDTYPE>("LS",
                global_indexes, comm, gm_->getMatrix(), *submat_indexing_);

            localX_
                = new SquareLocalMatrices<MATDTYPE>(subdiv_, chromatic_number_);
            localT_
                = new SquareLocalMatrices<MATDTYPE>(subdiv_, chromatic_number_);
        }
    }
#endif

#ifdef USE_MPI
    sH_ = new dist_matrix::SparseDistMatrix<DISTMATDTYPE>(comm, *matH_,
        remote_tasks_DistMatrix_, sparse_distmatrix_nb_partitions);
#else
    sH_ = new dist_matrix::SparseDistMatrix<DISTMATDTYPE>(
        0, *matH_, remote_tasks_DistMatrix_);
#endif
}

void ProjectedMatrices::computeInvS()
{
    compute_inverse_tm_.start();
#ifdef PRINT_OPERATIONS
    if (onpe0) (*MPIdata::sout) << "ProjectedMatrices::computeInvS()" << endl;
#endif
    gm_->computeInverse();
    compute_inverse_tm_.stop();
}

void ProjectedMatrices::rotateAll(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& rotation_matrix,
    const bool flag_eigen)
{
    // S -> U^T S U
    // rotate overlap and l_s
    if (flag_eigen)
    {
        gm_->set2Id(-1);
    }
    else
    {
        gm_->rotateAll(rotation_matrix);
    }
    //(*MPIdata::sout)<<"matS"<<endl;
    // matS_->print((*MPIdata::sout),0,0,5,5);

    // rotate matH_
    rotateSym(*matH_, rotation_matrix, *work_);

    computeInvS();

    // theta = invB * matH_
    updateTheta();

    updateHB();

    dm_->rotate(rotation_matrix, flag_eigen);
}

void ProjectedMatrices::applyInvS(SquareLocalMatrices<MATDTYPE>& mat)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> pmatrix("pmatrix", dim_, dim_);

    mat.fillDistMatrix(pmatrix, global_indexes_);

    gm_->applyInv(pmatrix);

    mat.init(pmatrix, global_indexes_, *submat_indexing_);
}

void ProjectedMatrices::setDMto2InvS()
{
    dm_->setto2InvS(gm_->getInverse(), gm_->getAssociatedOrbitalsIndex());
}

void ProjectedMatrices::solveGenEigenProblem(
    dist_matrix::DistMatrix<DISTMATDTYPE>& z, vector<DISTMATDTYPE>& val,
    char job)
{
    assert(val.size() == eigenvalues_.size());

    sygv_tm_.start();

    dist_matrix::DistMatrix<DISTMATDTYPE> mat(*matHB_);

    // Transform the generalized eigenvalue problem to a standard form
    gm_->sygst(mat);

    // solve a standard symmetric eigenvalue problem
    mat.syev(job, 'l', eigenvalues_, *u_);

#ifdef PROCRUSTES
    if (u0_ != 0)
    {
        work_->gemm('t', 'n', 1., *u_, *u0_, 0.);
        for (int j = 0; j < dim_; j++)
        {
            double val = 0.;
            int i      = work_->iamax(j, val);
            if (val < 0.) u_->scalColumn(j, -1.);
        }
        procrustes(*u_, *u0_, *p_);
    }
#endif

    // Get the eigenvectors Z of the generalized eigenvalue problem
    // Solve Z=L**(-T)*U
    z = *u_;
    gm_->solveLST(z);

    val = eigenvalues_;

    sygv_tm_.stop();
}

void ProjectedMatrices::buildDM(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& z, const int orbitals_index)
{
    dm_->build(z, orbitals_index);
}

void ProjectedMatrices::buildDM(const dist_matrix::DistMatrix<DISTMATDTYPE>& z,
    const vector<DISTMATDTYPE>& occ, const int orbitals_index)
{
    dm_->build(z, occ, orbitals_index);
}

void ProjectedMatrices::buildDM(
    const vector<DISTMATDTYPE>& occ, const int orbitals_index)
{
    dm_->build(occ, orbitals_index);
}

void ProjectedMatrices::updateDMwithEigenstates(const int iterative_index)
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "ProjectedMatrices: Compute DM using eigenstates\n";

    dist_matrix::DistMatrix<DISTMATDTYPE> zz("Z", dim_, dim_);
    vector<DISTMATDTYPE> val(dim_);

    // solves generalized eigenvalue problem
    // and return solution in zz and val
    solveGenEigenProblem(zz, val);
    setAuxilliaryEnergiesFromEigenenergies();
    double final_mu = computeChemicalPotentialAndOccupations();
    if (onpe0 && ct.verbose > 1) cout << "Final mu_ = " << final_mu << endl;

    // Build the density matrix X
    // X = Z * gamma * Z^T
    buildDM(zz, iterative_index);
}

void ProjectedMatrices::updateDM(const int iterative_index)
{
    Control& ct = *(Control::instance());

    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver)
        updateDMwithEigenstates(iterative_index);
    else
    {
        cerr<<"Eigensolver not available in ProjectedMatrices::updateDM()\n";
        ct.global_exit(2);
    }
}

void ProjectedMatrices::updateDMwithEigenstatesAndRotate(
    const int iterative_index, dist_matrix::DistMatrix<DISTMATDTYPE>& zz)
{
    vector<DISTMATDTYPE> val(dim_);

    // solves generalized eigenvalue problem
    // and return solution in zz and val
    solveGenEigenProblem(zz, val);
    setAuxilliaryEnergiesFromEigenenergies();
    computeChemicalPotentialAndOccupations();

    rotateAll(zz, true);

    dm_->build(iterative_index);
}

#if 0
void ProjectedMatrices::rotateBackDM()
{
    computeOccupations(occupation0_);
    *u0_=*work_;

#ifdef PROCRUSTES
    procrustes(*u_,*u0_,*p_);
#endif

    dist_matrix::DistMatrix<DISTMATDTYPE> rot("R", dim_, dim_);
    rot.gemm('n','t',1.,*u_,*u0_,0.);
//int nst=rot.m();
//(*MPIdata::sout)<<"u0"<<endl;
//u0_->print((*MPIdata::sout),0,0,nst,nst);
//(*MPIdata::sout)<<"u"<<endl;
//u_->print((*MPIdata::sout),0,0,nst,nst);
//(*MPIdata::sout)<<"rot"<<endl;
//rot.print((*MPIdata::sout),0,0,nst,nst);

#if 0
    (*MPIdata::sout)<<"Test unitarity of u0_"<<endl;
    dist_matrix::DistMatrix<DISTMATDTYPE> u0(*u0_);
    dist_matrix::DistMatrix<DISTMATDTYPE> test("t",nst,nst);
    test.gemm('t','n',1.,u0,*u0_,0.);
    for(int i=0;i<test.mloc();i++)
    for(int j=0;j<i;j++){
        if( fabs( test.val(i+j*test.mloc()) )>1.e-8 )
            (*MPIdata::sout)<<"test["<<i<<"]["<<j<<"]="
                <<test.val(i+j*test.mloc())<<endl;
    }
    for(int i=0;i<test.mloc();i++){
        if( fabs( test.val(i+i*test.mloc())-1. )>1.e-8 )
            (*MPIdata::sout)<<"test["<<i<<"]["<<i<<"]="
                <<test.val(i+i*test.mloc())<<endl;
    }
#endif
#if 0
    (*MPIdata::sout)<<"Test unitarity of u_"<<endl;
    dist_matrix::DistMatrix<DISTMATDTYPE> u(*u_);
    test.gemm('t','n',1.,u,*u_,0.);
    for(int i=0;i<test.mloc();i++)
    for(int j=0;j<i;j++){
        if( fabs( test.val(i+j*test.mloc()) )>1.e-8 )
            (*MPIdata::sout)<<"test["<<i<<"]["<<j<<"]="
                <<test.val(i+j*test.mloc())<<endl;
    }
    for(int i=0;i<test.mloc();i++){
        if( fabs( test.val(i+i*test.mloc())-1. )>1.e-8 )
            (*MPIdata::sout)<<"test["<<i<<"]["<<i<<"]="
                <<test.val(i+i*test.mloc())<<endl;
    }
#endif

    sqrtDistMatrix(rot);
    for(int i=0;i<dim_;i++)
        occupation_[i]=0.5*( occupation0_[i]+occupation_[i] );

    u_->gemm('n','n',1.,rot,*u0_,0.);

    dist_matrix::DistMatrix<DISTMATDTYPE> z(*u_);
    ls_->trtrs('l', 't', 'n', z);

    dist_matrix::DistMatrix<DISTMATDTYPE> gamma("Gamma", &occupation_[0], dim_, dim_);
    const double orbital_occupation = with_spin_ ? 1. : 2.;
    gamma.scal(orbital_occupation); // rescale for spin
 
    work_->symm('r', 'l', 1., gamma, z, 0.);
    dm_->gemm('n', 't', 1., *work_, z, 0.);
}
#endif

void ProjectedMatrices::computeOccupationsFromDM()
{
    Control& ct = *(Control::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver && dim_ > 0)
    {
        assert(dm_ != 0);
        dm_->computeOccupations(gm_->getCholeskyL());
    }
}

void ProjectedMatrices::getOccupations(vector<DISTMATDTYPE>& occ) const
{
    dm_->getOccupations(occ);
}
void ProjectedMatrices::setOccupations(const vector<DISTMATDTYPE>& occ)
{
    dm_->setOccupations(occ);
}
void ProjectedMatrices::printDM(ostream& os) const { dm_->print(os); }

const dist_matrix::DistMatrix<DISTMATDTYPE>& ProjectedMatrices::dm() const
{
    assert(dm_ != 0);
    return dm_->getMatrix();
}
const dist_matrix::DistMatrix<DISTMATDTYPE>&
ProjectedMatrices::kernel4dot() const
{
    return dm_->kernel4dot();
}

double ProjectedMatrices::getNel() const
{
    double val = dm_->dot(gm_->getMatrix());
    if (with_spin_)
    {
        double tmp      = 0.;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduceSpin(&val, &tmp, 1, MPI_SUM);
        val = tmp;
    }

    return val;
}

double ProjectedMatrices::getEigSum()
{
    if (dim_ == 0) return 0.;

    eigsum_tm_.start();
    work_->symm('l', 'l', 1., *matHB_, gm_->getInverse(), 0.);

    // return sum in Ry
    double val = work_->trace();
    eigsum_tm_.stop();

    return val;
}

double ProjectedMatrices::getExpectationH()
{
    if (dim_ == 0) return 0.;

    return getExpectation(*matHB_);
}

double ProjectedMatrices::getExpectation(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& A)
{
    return dm_->getExpectation(A);
}

// strip dm from the overlap contribution
// dm <- Ls**T * dm * Ls
void ProjectedMatrices::stripDM()
{
#ifdef PRINT_OPERATIONS
    if (onpe0) cout << "ProjectedMatrices::stripDM()" << endl;
#endif
#ifdef DEBUG // TEST
    double dd = dm_->getMatrix().trace();
    if (onpe0) cout << "test:  Trace DM = " << dd << endl;
    if (dm_->getMatrix().active()) assert(dd > 0.);
#endif
    dm_->stripS(gm_->getCholeskyL(), gm_->getAssociatedOrbitalsIndex());
}

void ProjectedMatrices::dressupDM()
{
#ifdef PRINT_OPERATIONS
    if (onpe0) cout << "ProjectedMatrices::dressupDM()" << endl;
#endif
    dm_->dressUpS(gm_->getCholeskyL(), gm_->getAssociatedOrbitalsIndex());
}

double ProjectedMatrices::computeEntropy(const double kbt)
{
    return kbt * dm_->computeEntropy();
}

double ProjectedMatrices::computeEntropy()
{
    // if(onpe0)(*MPIdata::sout)<<"ProjectedMatrices::computeEntropy()"<<endl;
    // if(onpe0)(*MPIdata::sout)<<"width_="<<width_<<endl;

    Control& ct = *(Control::instance());
    double entropy;

    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver
        || dm_->fromUniformOccupations())
    {
        if (!occupationsUptodate())
        {
            computeOccupationsFromDM();
        }
        else
        {
            if (onpe0 && ct.verbose > 1)
                (*MPIdata::sout)
                    << "occupations uptodate, skip computation..." << endl;
        }
        entropy = computeEntropy(width_);
    }
    return entropy;
}

void ProjectedMatrices::printOccupations(ostream& os) const
{
    if (dm_->occupationsUptodate()) dm_->printOccupations(os);
}

double ProjectedMatrices::checkCond(const double tol, const bool flag)
{
    double rcond = computeCond();

    if (rcond > tol)
    {
        // ofstream tfile("s.mm", ios::out);
        // gm_->printMM(tfile);
        // tfile.close();
#ifdef USE_MPI
        MGmol_MPI& mgmolmpi = *(MGmol_MPI::instance());
        mgmolmpi.barrier();
#endif
        if (onpe0)
            (*MPIdata::sout)
                << " CONDITION NUMBER OF THE OVERLAP MATRIX EXCEEDS TOL: "
                << rcond << "!!!" << endl;
        Control& ct = *(Control::instance());
        if (flag) ct.global_exit(2);
    }
    return rcond;
}
////// TEMPLATE THIS FOR FLOAT OPTION ??
int ProjectedMatrices::writeDM_hdf5(HDFrestart& h5f_file)
{
    hid_t file_id = h5f_file.file_id();

    ReplicatedWorkSpace<double>& wspace(ReplicatedWorkSpace<double>::instance());

    wspace.initSquareMatrix(dm_->getMatrix());

    if (file_id < 0) return 0;

    hsize_t dims[2] = { dim_, dim_ };

    // filespace identifier
    hid_t dataspace = H5Screate_simple(2, dims, NULL);

    hid_t dset_id = H5Dcreate2(file_id, "/Density_Matrix", H5T_NATIVE_DOUBLE,
        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset_id < 0)
    {
        (*MPIdata::serr)
            << "ProjectedMatrices::write_dm_hdf5: H5Dcreate2 failed!!!" << endl;
        return -1;
    }

    hid_t memspace  = dataspace;
    hid_t filespace = dataspace;

    DISTMATDTYPE* work_matrix = wspace.square_matrix();
    herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
        H5P_DEFAULT, work_matrix);
    if (status < 0)
    {
        (*MPIdata::serr) << "Orbitals: H5Dwrite failed!!!" << endl;
        return -1;
    }

    status = H5Dclose(dset_id);
    if (status < 0)
    {
        (*MPIdata::serr)
            << "ProjectedMatrices::write_dm_hdf5(), H5Dclose failed!!!" << endl;
        return -1;
    }
    status = H5Sclose(dataspace);
    if (status < 0)
    {
        (*MPIdata::serr)
            << "ProjectedMatrices::write_dm_hdf5(), H5Sclose failed!!!" << endl;
        return -1;
    }

    return 0;
}
////// TEMPLATE THIS FOR FLOAT OPTION ??
int ProjectedMatrices::read_dm_hdf5(hid_t file_id)
{
    ReplicatedWorkSpace<double>& wspace(ReplicatedWorkSpace<double>::instance());
    DISTMATDTYPE* work_matrix = wspace.square_matrix();

    int ierr        = 0;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (mmpi.instancePE0())
    {
        hid_t dset_id = H5Dopen2(file_id, "/Density_Matrix", H5P_DEFAULT);
        if (dset_id < 0)
        {
            (*MPIdata::serr) << "H5Dopen failed for /Density_Matrix!!!" << endl;
        }
        else
        {
            ierr          = 1;
            herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                H5S_ALL, H5P_DEFAULT, work_matrix);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "H5Dread failed for /Density_Matrix!!!" << endl;
                return -1;
            }

            status = H5Dclose(dset_id);
            if (status < 0)
            {
                (*MPIdata::serr) << "H5Dclose failed!!!" << endl;
                return -1;
            }
        }
    }
    mmpi.bcast(&ierr, 1);
    if (ierr >= 0) wspace.mpiBcastSquareMatrix();
    if (ierr >= 0) dm_->initMatrix(work_matrix);

    return ierr;
}

void ProjectedMatrices::printEigenvalues(ostream& os) const
{
    printEigenvaluesHa(os);
}

void ProjectedMatrices::printEigenvaluesEV(ostream& os) const
{
    Control& ct = *(Control::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver && onpe0)
    {
        os << endl << " Eigenvalues [eV]:";

        // Print ten to a row.
        os.setf(ios::right, ios::adjustfield);
        os.setf(ios::fixed, ios::floatfield);
        os << setprecision(3);
        for (int i = 0; i < dim_; i++)
        {
            if ((i % 10) == 0) os << endl;
            os << setw(7) << RY2EV * eigenvalues_[i] << " ";
        }
        os << endl;

        if (width_ > 1.e-10)
            os << " FERMI ENERGY   = " << RY2EV * mu_ << "[eV]" << endl;
    }
}

void ProjectedMatrices::printEigenvaluesHa(ostream& os) const
{
    Control& ct = *(Control::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver && onpe0)
    {
        os << endl << " Eigenvalues [Ha]:";

        // Print ten to a row.
        os.setf(ios::right, ios::adjustfield);
        os.setf(ios::fixed, ios::floatfield);
        os << setprecision(3);
        for (int i = 0; i < dim_; i++)
        {
            if ((i % 10) == 0) os << endl;
            os << setw(7) << 0.5 * eigenvalues_[i] << " ";
        }
        os << endl;

        if (width_ > 1.e-10)
            os << " FERMI ENERGY   = " << 0.5 * mu_ << "[Ha]" << endl;
    }
}

void ProjectedMatrices::setAuxilliaryEnergiesFromEigenenergies()
{
    Control& ct = *(Control::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver)
    {
        assert((int)eigenvalues_.size() == dim_);

        if (aux_energies_.size() == 0) aux_energies_.resize(dim_);
        aux_energies_ = eigenvalues_;

        assert((int)aux_energies_.size() == dim_);
    }
}

void ProjectedMatrices::setAuxilliaryEnergiesFromOccupations()
{
    if (aux_energies_.size() == 0) aux_energies_.resize(dim_);

    const double eps   = 1.e-15;
    const double shift = width_ > 1.e-8 ? -log(eps) / width_ : 1000.;

    vector<DISTMATDTYPE> occ(dim_);
    dm_->getOccupations(occ);
    for (int i = 0; i < dim_; i++)
    {
        if (occ[i] < eps)
        {
            aux_energies_[i] = mu_ + shift;
        }
        else if ((1. - occ[i]) < eps)
        {
            aux_energies_[i] = mu_ - shift;
        }
        else
        {
            aux_energies_[i] = mu_ + width_ * log(1. / occ[i] - 1.);
        }
    }
    assert((int)aux_energies_.size() == dim_);
}

// find the Fermi level by a bisection
// algorithm adapted from numerical recipes, 2nd edition
// and fill orbitals accordingly (in fermi_distribution)
double ProjectedMatrices::computeChemicalPotentialAndOccupations(
    const std::vector<DISTMATDTYPE>& energies, const double width,
    const int nel, const int max_numst)
{
    assert(energies.size() > 0);
    assert(nel >= 0);

    Control& ct = *(Control::instance());
    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout)
            << "computeChemicalPotentialAndOccupations() with width=" << width
            << ", for " << nel << " electrons" << endl;

    const int maxit         = 100;
    const double charge_tol = 1.0e-12;

    const double emin = *std::min_element(energies.begin(), energies.end());
    const double emax = *std::max_element(energies.begin(), energies.end());

    double mu1 = emin - 0.001;
    double mu2 = emax + 10. * width;
    assert(mu1 < mu2);
    bool done = false;

    if (nel <= 0)
    {
        // if( onpe0 )
        //    (*MPIdata::sout)<<"computeChemicalPotentialAndOccupations() with
        //    nel="<<nel<<endl;
        mu1 = -10000.;
        mu2 = 10000.;
    }

    vector<DISTMATDTYPE> occ(dim_, 0.);

    if (2 * dim_ <= nel)
    {
        done = true;
        mu_  = mu2;
        for (int i = 0; i < dim_; i++)
            occ[i] = 1.;
    }

    double f2 = 0.;
    if (!done)
    {
        f2 = 2. * fermi_distribution(mu2, max_numst, width, energies, occ)
             - (double)nel;
        // no unoccupied states
        if (fabs(f2) < charge_tol)
        {
            done = true;
            mu_  = mu2;
        }
    }
    double f1 = 0.;
    if (!done)
    {
        f1 = 2. * fermi_distribution(mu1, max_numst, width, energies, occ)
             - (double)nel;
        if (fabs(f1) < charge_tol)
        {
            if (onpe0) (*MPIdata::sout) << "only unoccupied states" << endl;
            done = true;
            mu_  = mu1;
        }
    }

    if (!done)
    {
        if (f1 * f2 > 0.)
        {
            (*MPIdata::sout) << "ERROR: mu1=" << mu1 << ", mu2=" << mu2 << endl;
            (*MPIdata::sout) << "ERROR: f1=" << f1 << ", f2=" << f2 << endl;
            (*MPIdata::sout) << "nel=" << nel << ", width=" << width << endl;
            Control& ct = *(Control::instance());
            ct.global_exit(2);
        }

        double dmu;
        if (f1 < 0.)
        {
            mu_ = mu1;
            dmu = mu2 - mu1;
        }
        else
        {
            mu_ = mu2;
            dmu = mu1 - mu2;
        }

        // main loop
        int iter = 0;
        double f = 0.;
        do
        {
            iter++;

            dmu *= 0.5;
            mu1 = mu_ + dmu;
            f   = 2. * fermi_distribution(mu1, max_numst, width, energies, occ)
                - (double)nel;

            if (f <= 0.)
            {
                mu_ = mu1;
                f   = -f;
            }

        } while ((iter < maxit) && (f > charge_tol));

        if (f > charge_tol)
        {
            if (onpe0)
            {
                (*MPIdata::sout) << "WARNING: "
                                    "ProjectedMatrices::"
                                    "computeChemicalPotentialAndOccupations()"
                                 << endl;
                (*MPIdata::sout) << "Iterations did not converge to tolerance "
                                 << scientific << charge_tol << endl;
                (*MPIdata::sout) << "f= " << f << ", nel=" << nel
                                 << ", max_numst=" << max_numst << endl;
            }
        }
    }

    // if( onpe0 )
    //    (*MPIdata::sout)<<"computeChemicalPotentialAndOccupations() with mu="
    //        <<mu_<<endl;

    dm_->setOccupations(occ);

    return mu_;
}

void ProjectedMatrices::getLoewdinTransform(
    SquareLocalMatrices<MATDTYPE>& localP)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> mat(gm_->getMatrix());
    dist_matrix::DistMatrix<DISTMATDTYPE> vect("eigenvectors", dim_, dim_);
    vector<DISTMATDTYPE> eigenvalues(dim_);
    mat.syev('v', 'l', eigenvalues, vect);

    vector<DISTMATDTYPE> diag_values(dim_);
    for (int i = 0; i < dim_; i++)
        diag_values[i] = (DISTMATDTYPE)(1. / sqrt(eigenvalues[i]));

    dist_matrix::DistMatrix<DISTMATDTYPE> matP("P", dim_, dim_);
    matP.clear();
    matP.setDiagonal(diag_values);
    // if( onpe0 )
    //  (*MPIdata::sout)<<"LocGridOrbitals::orthonormalizeLoewdin() ---
    //  matP"<<endl;
    // matP.print((*MPIdata::sout),0,0,5,5);
    mat.symm('r', 'l', 1., matP, vect, 0.);
    matP.gemm('n', 't', 1., mat, vect, 0.);

    localP.init(matP, global_indexes_, *submat_indexing_);
}

double ProjectedMatrices::getTraceDiagProductWithInvS(
    vector<DISTMATDTYPE>& ddiag)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> diag("diag", dim_, dim_);
    diag.setDiagonal(ddiag);

    work_->gemm('n', 'n', 1., diag, gm_->getInverse(), 0.);

    return work_->trace();
}

void ProjectedMatrices::resetDotProductMatrices()
{
    if (gram_4dotProducts_ != 0) delete gram_4dotProducts_;
    gram_4dotProducts_ = new GramMatrix(*gm_);
}

double ProjectedMatrices::dotProductWithInvS(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    assert(gram_4dotProducts_ != 0);

    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", dim_, dim_);
    local_product.fillDistMatrix(ds, global_indexes_);

    dist_matrix::DistMatrix<DISTMATDTYPE> work("work", dim_, dim_);
    work.gemm('n', 'n', 1., ds, gram_4dotProducts_->getInverse(), 0.);

    return work.trace();
}

double ProjectedMatrices::dotProductWithDM(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", dim_, dim_);
    local_product.fillDistMatrix(ds, global_indexes_);

    dist_matrix::DistMatrix<DISTMATDTYPE> work("work", dim_, dim_);
    work.gemm('n', 'n', 0.5, ds, kernel4dot(), 0.);

    return work.trace();
}

double ProjectedMatrices::dotProductSimple(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    assert(dm_4dot_product_ != 0);

    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", dim_, dim_);
    local_product.fillDistMatrix(ds, global_indexes_);

    return ds.trace();
}

void ProjectedMatrices::printTimers(ostream& os)
{
    sygv_tm_.print(os);
    compute_inverse_tm_.print(os);
    compute_invB_tm_.print(os);
    update_theta_tm_.print(os);
    update_submatX_tm_.print(os);
    update_submatT_tm_.print(os);
    init_gram_matrix_tm_.print(os);
    eigsum_tm_.print(os);
    consolidate_H_tm_.print(os);
}

// Assumes SquareLocalMatrix object contains partial contributions
double ProjectedMatrices::computeTraceInvSmultMat(
    const SquareLocalMatrices<MATDTYPE>& mat)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> pmatrix("pmatrix", dim_, dim_);

    mat.fillDistMatrix(pmatrix, global_indexes_);

    gm_->applyInv(pmatrix);
    return pmatrix.trace();
    /*
        dist_matrix::DistMatrix<DISTMATDTYPE>
       work_matrix("work_matrix",dim_,dim_); work_matrix.symm('l', 'l', 1.,
       gm_->getInverse(), pmatrix, 0.);

        dist_matrix::DistMatrix<DISTMATDTYPE> pm =
       gm_->getInverse();//("pmatrix",dim_,dim_); double itrace = pm.trace();
        if(onpe0)cout<<"INVERSE trace = "<<itrace<<endl;

        return work_matrix.trace();
    */
}

double ProjectedMatrices::computeTraceInvSmultMatMultTheta(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
{
    assert(theta_ != 0);
    dist_matrix::DistMatrix<DISTMATDTYPE> pmat("pmat", dim_, dim_);

    // compute mat*theta_
    pmat.gemm('n', 'n', 1.0, mat, *theta_, 0.);

    // compute invS*pmat = invS*(mat*theta)
    gm_->applyInv(pmat);

    return pmat.trace();
}
