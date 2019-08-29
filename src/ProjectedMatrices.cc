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

#include "Control.h"
#include "DensityMatrix.h"
#include "GramMatrix.h"
#include "HDFrestart.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "Power.h"
#include "PowerGen.h"
#include "ReplicatedWorkSpace.h"
#include "SP2.h"
#include "SparseDistMatrix.h"
#include "fermi.h"
#include "tools.h"

#include <fstream>
#include <iomanip>

#define RY2EV 13.605804

using namespace std;

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
    = nullptr;
GramMatrix* ProjectedMatrices::gram_4dotProducts_  = nullptr;
DensityMatrix* ProjectedMatrices::dm_4dot_product_ = nullptr;

static int sparse_distmatrix_nb_partitions = 128;

ProjectedMatrices::ProjectedMatrices(const int ndim, const bool with_spin)
    : with_spin_(with_spin),
      dim_(ndim),
      dm_(new DensityMatrix(ndim)),
      gm_(new GramMatrix(ndim))
{
    width_   = 0.;
    min_val_ = 0.25;

    sH_ = nullptr;

    if (dim_ > 0)
    {
        eigenvalues_.resize(dim_);
    }

    if (dim_ > 0)
    {
        matH_.reset(new dist_matrix::DistMatrix<DISTMATDTYPE>("H", ndim, ndim));
        matHB_.reset(
            new dist_matrix::DistMatrix<DISTMATDTYPE>("HB", ndim, ndim));
        theta_.reset(
            new dist_matrix::DistMatrix<DISTMATDTYPE>("Theta", ndim, ndim));
        work_.reset(
            new dist_matrix::DistMatrix<DISTMATDTYPE>("work", ndim, ndim));
    }

    if (onpe0)
        (*MPIdata::sout)
            << "ProjectedMatrices: sparse_distmatrix_nb_tasks_per_partitions_="
            << sparse_distmatrix_nb_tasks_per_partitions_ << endl;
    n_instances_++;
}

ProjectedMatrices::~ProjectedMatrices()
{
    if (dim_ > 0)
    {
        delete sH_;
        sH_ = nullptr;
    }

    if (n_instances_ == 1)
    {
        if (gram_4dotProducts_ != nullptr)
        {
            delete gram_4dotProducts_;
            gram_4dotProducts_ = nullptr;
        }
    }

    n_instances_--;
};

void ProjectedMatrices::setup(
    const double kbt, const int nel, const vector<vector<int>>& global_indexes)
{
    assert(global_indexes.size() > 0);

    setupBase(kbt, nel, global_indexes.size(), global_indexes[0].size());

    global_indexes_ = global_indexes;

#ifdef USE_DIS_MAT
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSpin();

    if (sH_ != nullptr) delete sH_;

    if (dim_ > 0)
    {
        DistMatrix2SquareLocalMatrices::setup(
            comm, global_indexes, dm_->getMatrix());
        LocalMatrices2DistMatrix::setup(
            comm, global_indexes, remote_tasks_DistMatrix_);

        localX_.reset(
            new SquareLocalMatrices<MATDTYPE>(subdiv_, chromatic_number_));
        localT_.reset(
            new SquareLocalMatrices<MATDTYPE>(subdiv_, chromatic_number_));
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
    // build DistMatrix from SquareLocalMatrices
    dist_matrix::DistMatrix<DISTMATDTYPE> pmatrix("pmatrix", dim_, dim_);

    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();

    sl2dm->accumulate(mat, pmatrix, dim_);

    gm_->applyInv(pmatrix);

    // convert result back into a SquareLocalMatrices
    DistMatrix2SquareLocalMatrices* dm2sl
        = DistMatrix2SquareLocalMatrices::instance();
    dm2sl->convert(pmatrix, mat);
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
    mat.syev(job, 'l', eigenvalues_, z);

    // Get the eigenvectors Z of the generalized eigenvalue problem
    // Solve Z=L**(-T)*U
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
    double final_mu = computeChemicalPotentialAndOccupations();
    if (onpe0 && ct.verbose > 1) cout << "Final mu_ = " << final_mu << endl;

    // Build the density matrix X
    // X = Z * gamma * Z^T
    buildDM(zz, iterative_index);
}

//"replicated" implementation of SP2.
// Theta is replicated on each MPI task, and SP2 solve run independently
// by each MPI task
void ProjectedMatrices::updateDMwithSP2(const int iterative_index)
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "ProjectedMatrices: Compute DM using SP2\n";

    updateThetaAndHB();

    // generate replicated copy of theta_
    SquareLocalMatrices<double> theta(1, dim_);
    double* work_matrix = theta.getSubMatrix();
    theta_->allgather(work_matrix, dim_);

    double emin;
    double emax;
    double epsilon = 1.e-2;

    static Power<LocalVector<double>, SquareLocalMatrices<double>> power(dim_);

    power.computeEigenInterval(
        theta, emin, emax, epsilon, (onpe0 && ct.verbose > 1));
    if (onpe0 && ct.verbose > 1)
        cout << "emin=" << emin << ", emax=" << emax << endl;

    const bool distributed = false;
    SP2 sp2(ct.dm_tol, distributed);
    {
        // include all the indexes so that traces are computed for the whole
        // replicated matrix
        std::vector<int> ids(dim_);
        for (int i = 0; i < dim_; i++)
            ids[i] = i;
        double buffer = 0.1;
        sp2.initializeLocalMat(theta, emin - buffer, emax + buffer, ids);
    }

    sp2.solve(nel_, (ct.verbose > 1));

    dist_matrix::DistMatrix<DISTMATDTYPE> dm("dm", dim_, dim_);
#ifdef HAVE_BML
    bml_matrix_t* dummy = bml_zero_matrix(
        dense, double_real, thetaSP2.n(), thetaSP2.n(), sequential);
#else
    SquareLocalMatrices<MATDTYPE> dummy(1, theta.n());
#endif

    sp2.getDM(dm, gm_->getInverse(), dummy);
    dm_->setMatrix(dm, iterative_index);

#ifdef HAVE_BML
    bml_deallocate(&dummy);
#endif
}

void ProjectedMatrices::updateDM(const int iterative_index)
{
    Control& ct = *(Control::instance());

    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver)
        updateDMwithEigenstates(iterative_index);
    else if (ct.DMEigensolver() == DMEigensolverType::SP2)
        updateDMwithSP2(iterative_index);
    else
    {
        cerr << "Eigensolver not available in ProjectedMatrices::updateDM()\n";
        ct.global_exit(2);
    }

#ifndef NDEBUG
    double nel = getNel();
    cout << "ProjectedMatrices::updateDM(), nel = " << nel << std::endl;
    assert(fabs(nel - nel) < 1.e-2);
    double energy = getExpectationH();
    cout << "ProjectedMatrices::updateDM(), energy = " << energy << std::endl;
#endif
}

void ProjectedMatrices::updateDMwithEigenstatesAndRotate(
    const int iterative_index, dist_matrix::DistMatrix<DISTMATDTYPE>& zz)
{
    vector<DISTMATDTYPE> val(dim_);

    // solves generalized eigenvalue problem
    // and return solution in zz and val
    solveGenEigenProblem(zz, val);
    computeChemicalPotentialAndOccupations();

    rotateAll(zz, true);

    dm_->build(iterative_index);
}

void ProjectedMatrices::computeOccupationsFromDM()
{
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "ProjectedMatrices::computeOccupationsFromDM()"
                         << std::endl;
#endif
    if (dim_ > 0)
    {
        assert(dm_);
        dm_->computeOccupations(gm_->getCholeskyL());
    }
}

void ProjectedMatrices::getOccupations(vector<DISTMATDTYPE>& occ) const
{
    dm_->getOccupations(occ);
}
void ProjectedMatrices::setOccupations(const vector<DISTMATDTYPE>& occ)
{
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "ProjectedMatrices::setOccupations()" << endl;
#endif
    dm_->setOccupations(occ);
}
void ProjectedMatrices::printDM(ostream& os) const { dm_->print(os); }

const dist_matrix::DistMatrix<DISTMATDTYPE>& ProjectedMatrices::dm() const
{
    assert(dm_);
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
        || ct.DMEigensolver() == DMEigensolverType::SP2
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

    ReplicatedWorkSpace<double>& wspace(
        ReplicatedWorkSpace<double>::instance());

    wspace.initSquareMatrix(dm_->getMatrix());

    if (file_id < 0) return 0;

    hsize_t dims[2] = { dim_, dim_ };

    // filespace identifier
    hid_t dataspace = H5Screate_simple(2, dims, nullptr);

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
    ReplicatedWorkSpace<double>& wspace(
        ReplicatedWorkSpace<double>::instance());
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

void ProjectedMatrices::computeLoewdinTransform(
    SquareLocalMatrices<MATDTYPE>& localP, const int orb_index)
{
    // dm_->computeOccupations(gm_->getCholeskyL());

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
    //  (*MPIdata::sout)<<"ProjectedMatrices::getLoewdinTransform() ---
    //  matP"<<endl;
    // matP.print((*MPIdata::sout),0,0,5,5);
    mat.symm('r', 'l', 1., matP, vect, 0.);
    matP.gemm('n', 't', 1., mat, vect, 0.);

    // new Gram matrix is Identity
    setGram2Id(orb_index);

    // transform DM to reflect Loewdin orthonormalization
    for (int i = 0; i < dim_; i++)
        diag_values[i] = sqrt(eigenvalues[i]);
    dist_matrix::DistMatrix<DISTMATDTYPE> invLoewdin("invLoewdin", dim_, dim_);
    invLoewdin.clear();
    invLoewdin.setDiagonal(diag_values);
    mat.symm('r', 'l', 1., invLoewdin, vect, 0.);
    invLoewdin.gemm('n', 't', 1., mat, vect, 0.);

    dm_->transform(invLoewdin);

    // transform matHB_ to reflect Loewdin orthonormalization
    mat.symm('r', 'l', 1., *matHB_, vect, 0.);
    matHB_->gemm('n', 't', 1., mat, vect, 0.);

    DistMatrix2SquareLocalMatrices* dm2sl
        = DistMatrix2SquareLocalMatrices::instance();
    dm2sl->convert(matP, localP);
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
    if (gram_4dotProducts_ != nullptr) delete gram_4dotProducts_;
    gram_4dotProducts_ = new GramMatrix(*gm_);
}

double ProjectedMatrices::dotProductWithInvS(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    assert(gram_4dotProducts_ != 0);

    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", dim_, dim_);
    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();

    sl2dm->accumulate(local_product, ds, dim_);

    dist_matrix::DistMatrix<DISTMATDTYPE> work("work", dim_, dim_);
    work.gemm('n', 'n', 1., ds, gram_4dotProducts_->getInverse(), 0.);

    return work.trace();
}

double ProjectedMatrices::dotProductWithDM(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", dim_, dim_);
    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();

    sl2dm->accumulate(local_product, ds, dim_);

    dist_matrix::DistMatrix<DISTMATDTYPE> work("work", dim_, dim_);
    work.gemm('n', 'n', 0.5, ds, kernel4dot(), 0.);

    return work.trace();
}

double ProjectedMatrices::dotProductSimple(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    assert(dm_4dot_product_ != 0);

    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", dim_, dim_);
    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();

    sl2dm->accumulate(local_product, ds, dim_);

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

    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();
    sl2dm->accumulate(mat, pmatrix, dim_);

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

/* Use the power method to compute the extents of the spectrum of the
 * generalized eigenproblem.
 */
void ProjectedMatrices::computeGenEigenInterval(
    std::vector<double>& interval, const int maxits, const double pad)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> mat(*matHB_);

    static PowerGen power(dim_);

    power.computeGenEigenInterval(mat, *gm_, interval, maxits, pad);
}
