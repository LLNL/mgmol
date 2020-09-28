// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "ProjectedMatrices.h"

#include "Control.h"
#include "DensityMatrix.h"
#include "DistMatrix2SquareLocalMatrices.h"
#include "DistMatrixTools.h"
#include "HDFrestart.h"
#include "LocalMatrices2DistMatrix.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"
#include "Power.h"
#include "PowerGen.h"
#include "ReplicatedMatrix.h"
#include "ReplicatedVector.h"
#include "ReplicatedWorkSpace.h"
#include "SP2.h"
#include "SparseDistMatrix.h"
#include "SquareSubMatrix2DistMatrix.h"
#include "fermi.h"

#include <fstream>
#include <iomanip>

#define RY2EV 13.605804

template <class MatrixType>
short ProjectedMatrices<MatrixType>::n_instances_ = 0;

template <class MatrixType>
GramMatrix<MatrixType>* ProjectedMatrices<MatrixType>::gram_4dotProducts_
    = nullptr;
template <class MatrixType>
DensityMatrix<MatrixType>* ProjectedMatrices<MatrixType>::dm_4dot_product_
    = nullptr;

static int sparse_distmatrix_nb_partitions = 128;

template <class MatrixType>
ProjectedMatrices<MatrixType>::ProjectedMatrices(
    const int ndim, const bool with_spin, const double width)
    : ProjectedMatricesInterface(with_spin, width),
      with_spin_(with_spin),
      dim_(ndim),
      dm_(new DensityMatrix<MatrixType>(ndim)),
      gm_(new GramMatrix<MatrixType>(ndim))
{
    eigenvalues_.resize(dim_);

    matH_.reset(new MatrixType("H", ndim, ndim));
    matHB_.reset(new MatrixType("HB", ndim, ndim));
    theta_.reset(new MatrixType("Theta", ndim, ndim));
    work_.reset(new MatrixType("work", ndim, ndim));

    n_instances_++;
}

template <class MatrixType>
ProjectedMatrices<MatrixType>::~ProjectedMatrices()
{
    if (n_instances_ == 1)
    {
        if (gram_4dotProducts_ != nullptr)
        {
            delete gram_4dotProducts_;
            gram_4dotProducts_ = nullptr;
        }
    }

    n_instances_--;
}

template <>
void ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>::convert(
    const SquareLocalMatrices<MATDTYPE>& src,
    dist_matrix::DistMatrix<DISTMATDTYPE>& dst)
{
    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();
    sl2dm->accumulate(src, dst);
}

template <>
void ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>::convert(
    const dist_matrix::DistMatrix<DISTMATDTYPE>& src,
    SquareLocalMatrices<MATDTYPE>& dst)
{
    DistMatrix2SquareLocalMatrices* dm2sl
        = DistMatrix2SquareLocalMatrices::instance();
    dm2sl->convert(src, dst);
}

#ifdef HAVE_MAGMA
template <>
void ProjectedMatrices<ReplicatedMatrix>::convert(
    const SquareLocalMatrices<MATDTYPE>& src, ReplicatedMatrix& dst)
{
    dst.init(src.getSubMatrix(), dim_);
}

template <>
void ProjectedMatrices<ReplicatedMatrix>::convert(
    const ReplicatedMatrix& src, SquareLocalMatrices<MATDTYPE>& dst)
{
    src.get(dst.getSubMatrix(), dim_);
}
#endif

template <>
void ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>::setupMPI(
    const std::vector<std::vector<int>>& global_indexes)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSpin();

    DistMatrix2SquareLocalMatrices::setup(
        comm, global_indexes, dm_->getMatrix());
    LocalMatrices2DistMatrix::setup(comm, global_indexes);
}

#ifdef HAVE_MAGMA
template <>
void ProjectedMatrices<ReplicatedMatrix>::setupMPI(
    const std::vector<std::vector<int>>& global_indexes)
{
}
#endif

template <class MatrixType>
void ProjectedMatrices<MatrixType>::setup(
    const std::vector<std::vector<int>>& global_indexes)
{
    assert(global_indexes.size() > 0);

    setupBase(global_indexes.size(), global_indexes[0].size());

    global_indexes_ = global_indexes;

    setupMPI(global_indexes);

    localX_.reset(
        new SquareLocalMatrices<MATDTYPE>(subdiv_, chromatic_number_));
    localT_.reset(
        new SquareLocalMatrices<MATDTYPE>(subdiv_, chromatic_number_));

    localHl_.reset(new SquareLocalMatrices<MATDTYPE>(
        global_indexes.size(), global_indexes[0].size()));
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::computeInvS()
{
    compute_inverse_tm_.start();
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "ProjectedMatrices<MatrixType>::computeInvS()"
                         << std::endl;
#endif
    gm_->computeInverse();
    compute_inverse_tm_.stop();
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::rotateAll(
    const MatrixType& rotation_matrix, const bool flag_eigen)
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
    //(*MPIdata::sout)<<"matS"<<std::endl;
    // matS_->print((*MPIdata::sout),0,0,5,5);

    // rotate matH_
    rotateSym(*matH_, rotation_matrix, *work_);

    computeInvS();

    // theta = invB * matH_
    updateTheta();

    updateHB();

    dm_->rotate(rotation_matrix, flag_eigen);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::applyInvS(
    SquareLocalMatrices<MATDTYPE>& mat)
{
    // build DistMatrix from SquareLocalMatrices
    convert(mat, *work_);

    gm_->applyInv(*work_);

    // convert result back into a SquareLocalMatrices
    convert(*work_, mat);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::setDMto2InvS()
{
    dm_->setto2InvS(gm_->getInverse(), gm_->getAssociatedOrbitalsIndex());
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::solveGenEigenProblem(
    MatrixType& z, char job)
{
    sygv_tm_.start();

    MatrixType mat(*matHB_);

    // Transform the generalized eigenvalue problem to a standard form
    gm_->sygst(mat);

    // solve a standard symmetric eigenvalue problem
    mat.syev(job, 'l', eigenvalues_, z);

    // Get the eigenvectors Z of the generalized eigenvalue problem
    // Solve Z=L**(-T)*U
    gm_->solveLST(z);

    sygv_tm_.stop();
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::buildDM(
    const MatrixType& z, const int orbitals_index)
{
    dm_->build(z, orbitals_index);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::buildDM(const MatrixType& z,
    const std::vector<DISTMATDTYPE>& occ, const int orbitals_index)
{
    dm_->build(z, occ, orbitals_index);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::buildDM(
    const std::vector<DISTMATDTYPE>& occ, const int orbitals_index)
{
    dm_->build(occ, orbitals_index);
}

// Use Chebyshev approximation to compute chemical potential and density matrix
template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDMwithChebApproximation(
    const int iterative_index)
{
    Control& ct = *(Control::instance());
    if (onpe0)
        (*MPIdata::sout)
            << "ProjectedMatrices: Compute DM using Chebyshev approximation"
            << std::endl;

    // CHEBYSHEV APPROXIMATION
    // set pointer to member function to evaluate fermi distribution function
    funcptr_ = &ProjectedMatricesInterface::chebfunDM;

    // Compute interval for Chebyshev approximation
    computeGenEigenInterval(cheb_interval_, ct.dm_approx_power_maxits, 0.05);
    double emin = cheb_interval_[0];
    double emax = cheb_interval_[1];
    //    if (onpe0 && ct.verbose > 1) cout<<"emin ="<<emin<<", emax
    //    ="<<emax<<endl;

    // compute approximation order
    int order         = ct.dm_approx_order;
    const int ndigits = ct.dm_approx_ndigits;
    if (ndigits)
    {
        const double delE
            = (emax - emin) / 2; // scaling factor into range [-1,1]
        const double beta_s = delE / width_; // scale beta = 1/kbt into [-1, 1]
        const double dp_order = 2 * (ndigits - 1) * beta_s / 3;
        order = std::ceil(dp_order) < 2000 ? std::ceil(dp_order) : 2000;
    }
    // compute chemical potential and density matrix with Chebyshev
    // approximation.
    double final_mu = computeChemicalPotentialAndDMwithChebyshev(
        order, emin, emax, iterative_index);
    if (onpe0 && ct.verbose > 1)
        std::cout << "Final mu_ = " << final_mu << " [Ha]" << std::endl;
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDMwithEigenstates(
    const int iterative_index)
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "ProjectedMatrices: Compute DM using eigenstates\n";

    MatrixType zz("Z", dim_, dim_);

    // solves generalized eigenvalue problem
    // and return solution in zz and val
    solveGenEigenProblem(zz);
    computeChemicalPotentialAndOccupations();
    if (onpe0 && ct.verbose > 1)
        std::cout << "Final mu_ = " << 0.5 * mu_ << " [Ha]" << std::endl;

    // Build the density matrix X
    // X = Z * gamma * Z^T
    buildDM(zz, iterative_index);
}

//"replicated" implementation of SP2.
// Theta is replicated on each MPI task, and SP2 solve run independently
// by each MPI task
template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDMwithSP2(const int iterative_index)
{
    Control& ct = *(Control::instance());

    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "ProjectedMatrices: Compute DM using SP2\n";

    updateThetaAndHB();

    // generate replicated copy of theta_
    SquareLocalMatrices<double> theta(1, dim_);
    convert(*theta_, theta);

    double emin;
    double emax;
    double epsilon = 1.e-2;

    static Power<LocalVector<double>, SquareLocalMatrices<double>> power(dim_);

    power.computeEigenInterval(
        theta, emin, emax, epsilon, (onpe0 && ct.verbose > 1));
    if (onpe0 && ct.verbose > 1)
        std::cout << "emin=" << emin << ", emax=" << emax << std::endl;

    const bool distributed = false;
    SP2 sp2(ct.dm_tol, distributed);
    {
        // include all the indexes so that traces are computed for the whole
        // replicated matrix
        std::vector<int> ids(dim_);
        for (unsigned int i = 0; i < dim_; i++)
            ids[i] = i;
        double buffer = 0.1;
        sp2.initializeLocalMat(theta, emin - buffer, emax + buffer, ids);
    }

    const double nel = with_spin_ ? nel_ : 2. * nel_;
    sp2.solve(nel, (ct.verbose > 1));

    MatrixType dm("dm", dim_, dim_);

    sp2.getDM(dm, gm_->getInverse());
    dm_->setMatrix(dm, iterative_index);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDM(const int iterative_index)
{
    Control& ct = *(Control::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver)
        updateDMwithEigenstates(iterative_index);
    else if (ct.DMEigensolver() == DMEigensolverType::Chebyshev)
        updateDMwithChebApproximation(iterative_index);
    else if (ct.DMEigensolver() == DMEigensolverType::SP2)
        updateDMwithSP2(iterative_index);
    else
    {
        std::cerr << "Eigensolver not available in "
                     "ProjectedMatrices<MatrixType>::updateDM()\n";
        ct.global_exit(2);
    }

#ifndef NDEBUG
    double nel = getNel();
    std::cout << "ProjectedMatrices<MatrixType>::updateDM(), nel = " << nel
              << std::endl;
    assert(std::isfinite(nel));
    double energy = getExpectationH();
    std::cout << "ProjectedMatrices<MatrixType>::updateDM(), energy = "
              << energy << std::endl;
#endif
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDMwithEigenstatesAndRotate(
    const int iterative_index, MatrixType& zz)
{
    // solves generalized eigenvalue problem
    // and return solution in zz
    solveGenEigenProblem(zz);
    computeChemicalPotentialAndOccupations();

    rotateAll(zz, true);

    dm_->build(iterative_index);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::computeOccupationsFromDM()
{
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout)
            << "ProjectedMatrices<MatrixType>::computeOccupationsFromDM()"
            << std::endl;
#endif
    Control& ct = *(Control::instance());
    if (ct.DMEigensolver() != DMEigensolverType::Chebyshev)
    {
        assert(dm_);
        dm_->computeOccupations(gm_->getCholeskyL());
    }
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::getOccupations(
    std::vector<DISTMATDTYPE>& occ) const
{
    dm_->getOccupations(occ);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::setOccupations(
    const std::vector<DISTMATDTYPE>& occ)
{
#ifdef PRINT_OPERATIONS
    if (onpe0)
        (*MPIdata::sout) << "ProjectedMatrices<MatrixType>::setOccupations()"
                         << std::endl;
#endif
    dm_->setOccupations(occ);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::printDM(std::ostream& os) const
{
    dm_->print(os);
}

template <class MatrixType>
const MatrixType& ProjectedMatrices<MatrixType>::dm() const
{
    assert(dm_);
    return dm_->getMatrix();
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::getNel() const
{
    double val = dm_->dot(gm_->getMatrix());
    if (with_spin_)
    {
        std::cout << "nel for 1 spin = " << val << std::endl;
        double tmp      = 0.;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduceSpin(&val, &tmp, 1, MPI_SUM);
        val = tmp;
    }

    return val;
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::getEigSum()
{
    eigsum_tm_.start();
    work_->symm('l', 'l', 1., *matHB_, gm_->getInverse(), 0.);

    // return sum in Ry
    double val = work_->trace();
    eigsum_tm_.stop();

    return val;
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::getExpectationH()
{
    return getExpectation(*matHB_);
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::getExpectation(const MatrixType& A)
{
    return dm_->getExpectation(A);
}

// strip dm from the overlap contribution
// dm <- Ls**T * dm * Ls
template <class MatrixType>
void ProjectedMatrices<MatrixType>::stripDM()
{
#ifdef PRINT_OPERATIONS
    if (onpe0)
        std::cout << "ProjectedMatrices<MatrixType>::stripDM()" << std::endl;
#endif
#ifdef DEBUG // TEST
    double dd = dm_->getMatrix().trace();
    if (onpe0) std::cout << "test:  Trace DM = " << dd << std::endl;
    if (dm_->getMatrix().active()) assert(dd > 0.);
#endif
    dm_->stripS(gm_->getCholeskyL());
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::dressupDM()
{
#ifdef PRINT_OPERATIONS
    if (onpe0)
        std::cout << "ProjectedMatrices<MatrixType>::dressupDM()" << std::endl;
#endif
    dm_->dressUpS(gm_->getCholeskyL(), gm_->getAssociatedOrbitalsIndex());
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::computeEntropy(const double kbt)
{
    return kbt * dm_->computeEntropy();
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::computeEntropy()
{
    // if(onpe0)(*MPIdata::sout)<<"ProjectedMatrices<MatrixType>::computeEntropy()"<<std::endl;
    // if(onpe0)(*MPIdata::sout)<<"width_="<<width_<<std::endl;
    compute_entropy_tm_.start();

    Control& ct    = *(Control::instance());
    double entropy = 0.;

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
                    << "occupations uptodate, skip computation..." << std::endl;
        }
        entropy = computeEntropy(width_);
    }
    else
    {
        entropy = computeEntropyWithCheb(width_);
    }

    compute_entropy_tm_.stop();

    return entropy;
}

// compute entropy using Chebyshev Approximation
template <class MatrixType>
double ProjectedMatrices<MatrixType>::computeEntropyWithCheb(const double kbt)
{
    Control& ct = *(Control::instance());

    // compute matrix variable X.S for Chebyshev
    // scale with 1/spin
    MGmol_MPI& mmpi           = *(MGmol_MPI::instance());
    double orbital_occupation = mmpi.nspin() > 1 ? 1. : 2.;
    const double scal         = 1 / orbital_occupation;
    MatrixType pmat("DM-Gram", dim_, dim_);
    pmat.gemm('N', 'N', scal, dm_->getMatrix(), gm_->getMatrix(), 0.);

    const double emin = 0.;
    const double emax = 1.;

    if (onpe0 && ct.verbose > 1)
        (*MPIdata::sout) << "computeEntropyWithChebyshev "
                         << "emin = " << emin << " emax = " << emax
                         << std::endl;

    // set pointer to member function to evaluate entropy function
    funcptr_ = &ProjectedMatricesInterface::chebfunEntropyFromOcc;
    // construct ChebyshevApproximation object
    const int order = 1000;
    static ChebyshevApproximation<MatrixType> chebapp(emin, emax, order, this);
    static bool recompute_entropy_coeffs = true;

    // compute Chebyshev approximation
    MatrixType mat
        = chebapp.computeChebyshevApproximation(pmat, recompute_entropy_coeffs);

    recompute_entropy_coeffs = false;
    // compute trace
    const double ts = mat.trace();
    //    if(onpe0 && ct.verbose > 1)(*MPIdata::sout)<<"entropy =
    //    "<<orbital_occupation*kbt*entropy<<std::endl;

    return -orbital_occupation * kbt * ts;
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::printOccupations(std::ostream& os) const
{
    if (dm_->occupationsUptodate()) dm_->printOccupations(os);
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::checkCond(
    const double tol, const bool flag)
{
    double rcond = computeCond();

    if (rcond > tol)
    {
        // ofstream tfile("s.mm", ios::out);
        // gm_->printMM(tfile);
        // tfile.close();
        MGmol_MPI& mgmolmpi = *(MGmol_MPI::instance());
        mgmolmpi.barrier();
        if (onpe0)
            (*MPIdata::sout)
                << " CONDITION NUMBER OF THE OVERLAP MATRIX EXCEEDS TOL: "
                << rcond << "!!!" << std::endl;
        Control& ct = *(Control::instance());
        if (flag) ct.global_exit(2);
    }
    return rcond;
}

////// TEMPLATE THIS FOR FLOAT OPTION ??
template <class MatrixType>
int ProjectedMatrices<MatrixType>::writeDM_hdf5(HDFrestart& h5f_file)
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
        (*MPIdata::serr) << "ProjectedMatrices<MatrixType>::write_dm_hdf5: "
                            "H5Dcreate2 failed!!!"
                         << std::endl;
        return -1;
    }

    hid_t memspace  = dataspace;
    hid_t filespace = dataspace;

    DISTMATDTYPE* work_matrix = wspace.square_matrix();
    herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
        H5P_DEFAULT, work_matrix);
    if (status < 0)
    {
        (*MPIdata::serr) << "Orbitals: H5Dwrite failed!!!" << std::endl;
        return -1;
    }

    status = H5Dclose(dset_id);
    if (status < 0)
    {
        (*MPIdata::serr) << "ProjectedMatrices<MatrixType>::write_dm_hdf5(), "
                            "H5Dclose failed!!!"
                         << std::endl;
        return -1;
    }
    status = H5Sclose(dataspace);
    if (status < 0)
    {
        (*MPIdata::serr) << "ProjectedMatrices<MatrixType>::write_dm_hdf5(), "
                            "H5Sclose failed!!!"
                         << std::endl;
        return -1;
    }

    return 0;
}
////// TEMPLATE THIS FOR FLOAT OPTION ??
template <class MatrixType>
int ProjectedMatrices<MatrixType>::read_dm_hdf5(hid_t file_id)
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
            (*MPIdata::serr)
                << "H5Dopen failed for /Density_Matrix!!!" << std::endl;
        }
        else
        {
            ierr          = 1;
            herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                H5S_ALL, H5P_DEFAULT, work_matrix);
            if (status < 0)
            {
                (*MPIdata::serr)
                    << "H5Dread failed for /Density_Matrix!!!" << std::endl;
                return -1;
            }

            status = H5Dclose(dset_id);
            if (status < 0)
            {
                (*MPIdata::serr) << "H5Dclose failed!!!" << std::endl;
                return -1;
            }
        }
    }
    mmpi.bcast(&ierr, 1);
    if (ierr >= 0) wspace.mpiBcastSquareMatrix();
    if (ierr >= 0) dm_->initMatrix(work_matrix);

    return ierr;
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::printEigenvalues(std::ostream& os) const
{
    printEigenvaluesHa(os);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::printEigenvaluesEV(std::ostream& os) const
{
    Control& ct = *(Control::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver && onpe0)
    {
        os << std::endl << " Eigenvalues [eV]:";

        // Print ten to a row.
        os.setf(std::ios::right, std::ios::adjustfield);
        os.setf(std::ios::fixed, std::ios::floatfield);
        os << std::setprecision(3);
        for (unsigned int i = 0; i < dim_; i++)
        {
            if ((i % 10) == 0) os << std::endl;
            os << std::setw(7) << RY2EV * eigenvalues_[i] << " ";
        }
        os << std::endl;

        if (width_ > 1.e-10)
            os << " FERMI ENERGY   = " << RY2EV * mu_ << "[eV]" << std::endl;
    }
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::printEigenvaluesHa(std::ostream& os) const
{
    Control& ct = *(Control::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver && onpe0)
    {
        os << std::endl << " Eigenvalues [Ha]:";

        // Print ten to a row.
        os.setf(std::ios::right, std::ios::adjustfield);
        os.setf(std::ios::fixed, std::ios::floatfield);
        os << std::setprecision(3);
        for (unsigned int i = 0; i < dim_; i++)
        {
            if ((i % 10) == 0) os << std::endl;
            os << std::setw(7) << 0.5 * eigenvalues_[i] << " ";
        }
        os << std::endl;

        if (width_ > 1.e-10)
            os << " FERMI ENERGY   = " << 0.5 * mu_ << "[Ha]" << std::endl;
    }
}

// find the Fermi level
// and fill orbitals accordingly (in fermi_distribution)
template <class MatrixType>
void ProjectedMatrices<MatrixType>::computeChemicalPotentialAndOccupations(
    const std::vector<DISTMATDTYPE>& energies, const double width,
    const int max_numst)
{
    assert(energies.size() > 0);
    assert(nel_ >= 0);

    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (mmpi.instancePE0() && ct.verbose > 1)
        (*MPIdata::sout)
            << "computeChemicalPotentialAndOccupations() with width=" << width
            << ", for " << nel_ << " electrons" << std::endl;

    std::vector<DISTMATDTYPE> occ(dim_, 0.);

    mu_ = compute_chemical_potential_and_occupations(
        energies, width, nel_, max_numst, onpe0, occ);
    // if( onpe0 )
    //    (*MPIdata::sout)<<"computeChemicalPotentialAndOccupations() with mu="
    //        <<mu<<std::endl;

    dm_->setOccupations(occ);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::computeLoewdinTransform(
    SquareLocalMatrices<MATDTYPE>& localP, const int orb_index,
    const bool transform_matrices)
{
    assert(gm_ != nullptr);

    MatrixType invSqrtMat("invSqrtMat", dim_, dim_);

    std::shared_ptr<MatrixType> sqrtMat;
    if (transform_matrices)
    {
        sqrtMat.reset(new MatrixType("sqrtMat", dim_, dim_));
    }

    gm_->computeLoewdinTransform(invSqrtMat, sqrtMat, orb_index);

    if (transform_matrices)
    {
        assert(sqrtMat);

        // transform DM to reflect Loewdin orthonormalization
        dm_->transform(*sqrtMat);

        // transform matHB_ to reflect Loewdin orthonormalization
        // (we reuse sqrtMat since we are done with it)
        MatrixType& mat(*sqrtMat);
        mat.symm('r', 'l', 1., *matHB_, invSqrtMat, 0.);
        matHB_->gemm('n', 't', 1., mat, invSqrtMat, 0.);
    }

    convert(invSqrtMat, localP);
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::getTraceDiagProductWithInvS(
    std::vector<DISTMATDTYPE>& ddiag)
{
    return gm_->getTraceDiagProductWithInvS(ddiag);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::resetDotProductMatrices()
{
    if (gram_4dotProducts_ != nullptr) delete gram_4dotProducts_;
    gram_4dotProducts_ = new GramMatrix<MatrixType>(*gm_);
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::dotProductWithInvS(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    assert(gram_4dotProducts_ != nullptr);

    MatrixType ds("ds", dim_, dim_);

    convert(local_product, ds);

    work_->gemm('n', 'n', 1., ds, gram_4dotProducts_->getInverse(), 0.);

    return work_->trace();
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::dotProductWithDM(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    MatrixType ds("ds", dim_, dim_);

    convert(local_product, ds);

    work_->gemm('n', 'n', 0.5, ds, dm_->kernel4dot(), 0.);

    return work_->trace();
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::dotProductSimple(
    const SquareLocalMatrices<MATDTYPE>& local_product)
{
    convert(local_product, *work_);

    return work_->trace();
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::printTimers(std::ostream& os)
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
    compute_entropy_tm_.print(os);
}

// Assumes SquareLocalMatrix object contains partial contributions
template <class MatrixType>
double ProjectedMatrices<MatrixType>::computeTraceInvSmultMat(
    const SquareLocalMatrices<MATDTYPE>& mat)
{
    convert(mat, *work_);

    gm_->applyInv(*work_);
    return work_->trace();
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::computeTraceInvSmultMatMultTheta(
    const MatrixType& mat)
{
    assert(theta_ != nullptr);

    // compute mat*theta_
    work_->gemm('n', 'n', 1.0, mat, *theta_, 0.);

    // compute invS*pmat = invS*(mat*theta)
    gm_->applyInv(*work_);

    return work_->trace();
}

template <class MatrixType>
double
ProjectedMatrices<MatrixType>::computeChemicalPotentialAndDMwithChebyshev(
    const int order, const double emin, const double emax,
    const int iterative_index)
{
    assert(emax > emin);
    assert(nel_ >= 0.);

    Control& ct = *(Control::instance());

    // create Chebyshev approximation object
    ChebyshevApproximation<MatrixType> chebapp(emin, emax, order, this);

    if (onpe0 && ct.verbose > 0)
        (*MPIdata::sout)
            << "computeChemicalPotentialAndDMWithChebyshev(), order = "
            << chebapp.order() << " with width=" << width_ << ", for " << nel_
            << " electrons" << std::endl;

    const int maxit         = 100;
    const double charge_tol = 1.0e-12;

    double mu1 = emin - 0.001;
    double mu2 = emax + 10. * width_;

    assert(mu1 < mu2);
    bool done = false;

    if (nel_ <= 0.)
    {
        mu1 = -10000.;
        mu2 = 10000.;
    }

    if (static_cast<double>(dim_) <= nel_)
    {
        done = true;
        mu_  = mu2;
    }

    // begin
    // compute matrix variable S^{-1}H for Chebyshev
    MatrixType mat(*matHB_);
    gm_->applyInv(mat);
    // build array of Chebyshev polynomials (Chebyshev Nodes)

    /// print matrices
    /*
             ofstream tfile("s.mm", ios::out);
             ofstream tfile2("h.mm", ios::out);
             gm_->printMM(tfile);
             matHB_->printMM(tfile2);
             tfile.close();
             tfile2.close();
    */
    //// end print

    chebapp.buildChebyshevNodes(emin, emax, mat);

    MatrixType tmp("TMP", dim_, dim_);
    MatrixType dm("DM", dim_, dim_);

    if (onpe0 && ct.verbose > 0)
        std::cout << "emin = " << emin << " emax = " << emax << std::endl;

    double f2 = 0.;
    if (!done)
    {
        mu_ = mu2;
        chebapp.computeChebyshevCoeffs();

        // compute Chebyshev approximation
        dm.gemm('N', 'N', 1., chebapp.computeChebyshevApproximation(),
            gm_->getInverse(), 0.);
        tmp.gemm('N', 'N', 1., dm, gm_->getMatrix(), 0.);
        // compute trace and check convergence
        f2 = tmp.trace() - nel_;

        // no unoccupied states
        if (std::abs(f2) < charge_tol)
        {
            done = true;
        }
    }
    double f1 = 0.;
    if (!done)
    {
        mu_ = mu1;
        chebapp.computeChebyshevCoeffs();

        // compute Chebyshev approximation
        dm.gemm('N', 'N', 1., chebapp.computeChebyshevApproximation(),
            gm_->getInverse(), 0.);
        tmp.gemm('N', 'N', 1., dm, gm_->getMatrix(), 0.);
        // compute trace and check convergence
        f1 = tmp.trace() - nel_;

        // no unoccupied states
        if (std::abs(f1) < charge_tol)
        {
            done = true;
        }
    }

    if (!done)
    {
        if (f1 * f2 > 0.)
        {
            (*MPIdata::sout)
                << "ERROR: mu1=" << mu1 << ", mu2=" << mu2 << std::endl;
            (*MPIdata::sout)
                << "ERROR: f1=" << f1 << ", f2=" << f2 << std::endl;
            (*MPIdata::sout)
                << "nel=" << nel_ << ", width=" << width_ << std::endl;
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
        int iter      = 0;
        double f      = 0.;
        double mu_old = mu_;
        do
        {
            iter++;

            dmu *= 0.5;
            mu_ = mu_old + dmu;

            chebapp.computeChebyshevCoeffs();
            // compute Chebyshev approximation
            tmp = chebapp.computeChebyshevApproximation();
            // compute trace and check convergence
            f = tmp.trace() - nel_;
            if (f <= 0.)
            {
                mu_old = mu_;
                f      = -f;
            }

        } while ((iter < maxit) && (f > charge_tol));

        // compute DM and occupations

        if (f > charge_tol)
        {
            if (onpe0)
            {
                (*MPIdata::sout)
                    << "WARNING: "
                       "ProjectedMatrices<MatrixType>::"
                       "computeChemicalPotentialAndDMwithChebyshev()"
                    << std::endl;
                (*MPIdata::sout) << "Iterations did not converge to tolerance "
                                 << std::scientific << charge_tol << std::endl;
            }
        }
    }
    // update mu1 and mu2
    mu1 = mu_ - 10. * width_;
    mu2 = mu_ + 10. * width_;

    // set density matrix
    dm.gemm('N', 'N', 1., tmp, gm_->getInverse(), 0.);
    MGmol_MPI& mmpi           = *(MGmol_MPI::instance());
    double orbital_occupation = mmpi.nspin() > 1 ? 1. : 2.;
    dm.scal(orbital_occupation);
    dm_->setMatrix(dm, iterative_index);

    return mu_;
}

/* Use the power method to compute the extents of the spectrum of the
 * generalized eigenproblem.
 */
template <>
void ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>::computeGenEigenInterval(
    std::vector<double>& interval, const int maxits, const double pad)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> mat(*matHB_);

    static PowerGen<dist_matrix::DistMatrix<DISTMATDTYPE>, dist_matrix::DistVector<DISTMATDTYPE>> power(dim_);

    power.computeGenEigenInterval(mat, *gm_, interval, maxits, pad);
}

#ifdef HAVE_MAGMA
template <>
void ProjectedMatrices<ReplicatedMatrix>::computeGenEigenInterval(
    std::vector<double>& interval, const int maxits, const double pad)
{
    ReplicatedMatrix mat(*matHB_);

    static PowerGen<ReplicatedMatrix, ReplicatedVector> power(dim_);

    power.computeGenEigenInterval(mat, *gm_, interval, maxits, pad);
}
#endif

template <>
void ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>::consolidateH()
{
    consolidate_H_tm_.start();

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSpin();

    dist_matrix::SparseDistMatrix<DISTMATDTYPE> slH(
        comm, *matH_, sparse_distmatrix_nb_partitions);

    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();
    sl2dm->convert(*localHl_, slH, dim_);

    SquareSubMatrix2DistMatrix* ss2dm = SquareSubMatrix2DistMatrix::instance();
    ss2dm->convert(*localHnl_, slH);

    slH.parallelSumToDistMatrix();

    consolidate_H_tm_.stop();
}

#ifdef HAVE_MAGMA
template <>
void ProjectedMatrices<ReplicatedMatrix>::consolidateH()
{
    consolidate_H_tm_.start();

    matH_->assign(*localHl_);
    matH_->add(*localHnl_);

    // sum up across MPI tasks
    matH_->consolidate();

    consolidate_H_tm_.stop();
}
#endif

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateSubMatX(const MatrixType& dm)
{
    convert(dm, *localX_);
}

template <>
SquareLocalMatrices<double>
ProjectedMatrices<dist_matrix::DistMatrix<double>>::getReplicatedDM()
{
    SquareLocalMatrices<double> sldm(1, dim_);
    const dist_matrix::DistMatrix<double>& dm(dm_->getMatrix());
    dm.allgather(sldm.getSubMatrix(), dim_);

    return sldm;
}

#ifdef HAVE_MAGMA
template <>
SquareLocalMatrices<double>
ProjectedMatrices<ReplicatedMatrix>::getReplicatedDM()
{
    SquareLocalMatrices<double> sldm(1, dim_);
    const ReplicatedMatrix& dm(dm_->getMatrix());
    dm.get(sldm.getSubMatrix(), dim_);

    return sldm;
}

#endif

template class ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>;
#ifdef HAVE_MAGMA
template class ProjectedMatrices<ReplicatedMatrix>;
#endif
