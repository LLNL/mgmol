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
#include "HDFrestart.h"
#include "LocalMatrices2ReplicatedMatrix.h"
#include "MGmol_MPI.h"
#include "Orbitals.h"
#include "Power.h"
#include "PowerGen.h"
#include "ReplicatedMatrix.h"
#include "ReplicatedMatrix2SquareLocalMatrices.h"
#include "ReplicatedVector.h"
#include "ReplicatedWorkSpace.h"
#include "SP2.h"
#include "fermi.h"
#include "hdf_tools.h"

#ifdef MGMOL_USE_SCALAPACK
#include "DistMatrix2SquareLocalMatrices.h"
#include "DistMatrixTools.h"
#include "DistVector.h"
#include "LocalMatrices2DistMatrix.h"
#include "SparseDistMatrix.h"
#include "SquareSubMatrix2DistMatrix.h"
#else
typedef double DISTMATDTYPE;
#endif

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

#ifdef MGMOL_USE_SCALAPACK
static int sparse_distmatrix_nb_partitions = 128;
#endif

template <>
std::string ProjectedMatrices<ReplicatedMatrix>::getMatrixType()
{
    return "ReplicatedMatrix";
}

#ifdef MGMOL_USE_SCALAPACK
template <>
std::string ProjectedMatrices<dist_matrix::DistMatrix<double>>::getMatrixType()
{
    return "DistMatrix<double>";
}
#endif

#ifdef MGMOL_USE_SCALAPACK
//
// conversion functions from one matrix format into another
//
#ifndef HAVE_MAGMA
void convert_matrix(const dist_matrix::DistMatrix<double>& src,
    SquareLocalMatrices<double, MemorySpace::Host>& dst)
{
    DistMatrix2SquareLocalMatrices* dm2sl
        = DistMatrix2SquareLocalMatrices::instance();
    dm2sl->convert(src, dst);
}
#else
void convert_matrix(const dist_matrix::DistMatrix<double>& src,
    SquareLocalMatrices<double, MemorySpace::Device>& dst)
{
    DistMatrix2SquareLocalMatrices* dm2sl
        = DistMatrix2SquareLocalMatrices::instance();
    SquareLocalMatrices<double, MemorySpace::Host> tmp(dst.nmat(), dst.m());
    dm2sl->convert(src, tmp);

    dst.assign(tmp);
}
#endif
#endif

#ifndef HAVE_MAGMA
void convert_matrix(const ReplicatedMatrix& src,
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& dst)
{
    assert(dst.m() > 0);

    ReplicatedMatrix2SquareLocalMatrices* r2l
        = ReplicatedMatrix2SquareLocalMatrices::instance();
    r2l->convert(src, dst);
}
#else
void convert_matrix(const ReplicatedMatrix& src,
    SquareLocalMatrices<MATDTYPE, MemorySpace::Device>& dst)
{
    dst.assign(src);
}
#endif

//=====================================================================//

template <class MatrixType>
ProjectedMatrices<MatrixType>::ProjectedMatrices(
    const int ndim, const bool with_spin, const double width)
    : ProjectedMatricesInterface(with_spin, width),
      with_spin_(with_spin),
      dim_(ndim),
      dm_(new DensityMatrix<MatrixType>(ndim)),
      gm_(new GramMatrix<MatrixType>(ndim))
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());

    if (mmpi.instancePE0() && ct.verbose > 1)
    {
        std::cout << "New ProjectedMatrices with MatrixType: "
                  << getMatrixType() << std::endl;
    }

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

#ifdef MGMOL_USE_SCALAPACK
template <>
void ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>::convert(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& src,
    dist_matrix::DistMatrix<DISTMATDTYPE>& dst)
{
    LocalMatrices2DistMatrix* sl2dm = LocalMatrices2DistMatrix::instance();
    sl2dm->accumulate(src, dst);
}
#endif

template <>
void ProjectedMatrices<ReplicatedMatrix>::convert(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& src,
    ReplicatedMatrix& dst)
{
    LocalMatrices2ReplicatedMatrix* sl2rm
        = LocalMatrices2ReplicatedMatrix::instance();

    sl2rm->accumulate(src, dst);
}

#ifdef MGMOL_USE_SCALAPACK
template <>
void ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>::
    setupGlobalIndexes(const std::vector<std::vector<int>>& global_indexes)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSpin();

    DistMatrix2SquareLocalMatrices::setup(
        comm, global_indexes, gm_->getMatrix());
    LocalMatrices2DistMatrix::setup(comm, global_indexes);
}
#endif

template <>
void ProjectedMatrices<ReplicatedMatrix>::setupGlobalIndexes(
    const std::vector<std::vector<int>>& global_indexes)
{
    LocalMatrices2ReplicatedMatrix::setup(global_indexes);

    ReplicatedMatrix2SquareLocalMatrices::setup(global_indexes);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::setup(
    const std::vector<std::vector<int>>& global_indexes)
{
    assert(global_indexes.size() > 0);

    setupBase(global_indexes.size(), global_indexes[0].size());

    global_indexes_ = global_indexes;

    setupGlobalIndexes(global_indexes);

    localX_.reset(new SquareLocalMatrices<MATDTYPE, memory_space_type>(
        subdiv_, chromatic_number_));
    localT_.reset(new SquareLocalMatrices<MATDTYPE, MemorySpace::Host>(
        subdiv_, chromatic_number_));

    localHl_.reset(new SquareLocalMatrices<MATDTYPE, MemorySpace::Host>(
        global_indexes.size(), global_indexes[0].size()));
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateSubMatT()
{
    convert_matrix(*theta_, *localT_);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::computeInvS()
{
    compute_inverse_tm_.start();
#ifdef PRINT_OPERATIONS
    if (mmpi.instancePE0())
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
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& mat)
{
    // build Matrix from SquareLocalMatrices
    convert(mat, *work_);

    gm_->applyInv(*work_);

    // convert result back into a SquareLocalMatrices
    convert_matrix(*work_, mat);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::setDMto2InvS()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());

    if (mmpi.instancePE0() && ct.verbose > 1)
        std::cout << "ProjectedMatrices::setDMto2InvS()..." << std::endl;

    dm_->setto2InvS(gm_->getInverse());
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
void ProjectedMatrices<MatrixType>::buildDM(const MatrixType& z)
{
    dm_->build(z);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::buildDM(
    const MatrixType& z, const std::vector<DISTMATDTYPE>& occ)
{
    dm_->build(z, occ);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::buildDM(
    const std::vector<DISTMATDTYPE>& occ)
{
    dm_->build(occ);
}

// Use Chebyshev approximation to compute chemical potential and density matrix
template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDMwithChebApproximation()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());

    if (mmpi.instancePE0())
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
    //    if (mmpi.instancePE0() && ct.verbose > 1) cout<<"emin ="<<emin<<",
    //    emax
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
    double final_mu
        = computeChemicalPotentialAndDMwithChebyshev(order, emin, emax);
    if (mmpi.instancePE0() && ct.verbose > 1)
        std::cout << "Final mu_ = " << final_mu << " [Ha]" << std::endl;
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDMwithEigenstates()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());

    if (mmpi.PE0() && ct.verbose > 1)
        (*MPIdata::sout) << "ProjectedMatrices: Compute DM using eigenstates\n";

    MatrixType zz("Z", dim_, dim_);

    // solves generalized eigenvalue problem
    // and return solution in zz and val
    solveGenEigenProblem(zz);
    computeChemicalPotentialAndOccupations();
    if (mmpi.instancePE0() && ct.verbose > 1)
        std::cout << "Final mu_ = " << 0.5 * mu_ << " [Ha]" << std::endl;

    // Build the density matrix X
    // X = Z * gamma * Z^T
    buildDM(zz);
}

//"replicated" implementation of SP2.
// Theta is replicated on each MPI task, and SP2 solve run independently
// by each MPI task
template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDMwithSP2()
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    Control& ct     = *(Control::instance());

    if (mmpi.instancePE0() && ct.verbose > 1)
        (*MPIdata::sout) << "ProjectedMatrices: Compute DM using SP2\n";

    updateThetaAndHB();

    // generate replicated copy of theta_
    SquareLocalMatrices<double, MemorySpace::Host> theta(1, dim_);
    convert_matrix(*theta_, theta);

    double emin;
    double emax;
    double epsilon = 1.e-2;

    static Power<LocalVector<double, MemorySpace::Host>,
        SquareLocalMatrices<double, MemorySpace::Host>>
        power(dim_);

    power.computeEigenInterval(
        theta, emin, emax, epsilon, (mmpi.instancePE0() && ct.verbose > 1));
    if (mmpi.instancePE0() && ct.verbose > 1)
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
    dm_->setMatrix(dm);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDM()
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver)
        updateDMwithEigenstates();
    else if (ct.DMEigensolver() == DMEigensolverType::Chebyshev)
        updateDMwithChebApproximation();
    else if (ct.DMEigensolver() == DMEigensolverType::SP2)
        updateDMwithSP2();
    else
    {
        std::cerr << "Eigensolver not available in "
                     "ProjectedMatrices<MatrixType>::updateDM()\n";
        mmpi.abort();
    }

#ifndef NDEBUG
    double nel = getNel();
    if (mmpi.instancePE0())
        std::cout << "ProjectedMatrices<MatrixType>::updateDM(), nel = " << nel
                  << std::endl;
    assert(std::isfinite(nel));
    double energy = getExpectationH();
    if (mmpi.instancePE0())
        std::cout << "ProjectedMatrices<MatrixType>::updateDM(), energy = "
                  << energy << std::endl;
#endif
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateDMwithEigenstatesAndRotate(
    MatrixType& zz)
{
    // solves generalized eigenvalue problem
    // and return solution in zz
    solveGenEigenProblem(zz);
    computeChemicalPotentialAndOccupations();

    rotateAll(zz, true);

    dm_->build(zz);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::computeOccupationsFromDM()
{
#ifdef PRINT_OPERATIONS
    if (mmpi.instancePE0())
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
    assert(!occ.empty());
#ifdef PRINT_OPERATIONS
    if (mmpi.instancePE0())
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
    if (with_spin_)
    {
        double tmp      = 0.;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduceSpin(&val, &tmp, 1, MPI_SUM);
        val = tmp;
    }

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
    double expectation = dm_->getExpectation(A);
    if (with_spin_)
    {
        double tmp      = 0.;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduceSpin(&expectation, &tmp, 1, MPI_SUM);
        expectation = tmp;
    }
    return expectation;
}

// strip dm from the overlap contribution
// dm <- Ls**T * dm * Ls
template <class MatrixType>
void ProjectedMatrices<MatrixType>::stripDM()
{
#ifdef PRINT_OPERATIONS
    if (mmpi.instancePE0())
        std::cout << "ProjectedMatrices<MatrixType>::stripDM()" << std::endl;
#endif
#ifdef DEBUG // TEST
    double dd = dm_->getMatrix().trace();
    if (mmpi.instancePE0())
        std::cout << "test:  Trace DM = " << dd << std::endl;
    if (dm_->getMatrix().active()) assert(dd > 0.);
#endif
    dm_->stripS(gm_->getCholeskyL());
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::dressupDM()
{
#ifdef PRINT_OPERATIONS
    if (mmpi.instancePE0())
        std::cout << "ProjectedMatrices<MatrixType>::dressupDM()" << std::endl;
#endif
    dm_->dressUpS(gm_->getCholeskyL());
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::computeEntropy(const double kbt)
{
    double entropy = dm_->computeEntropy();
    if (with_spin_)
    {
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        double tmp      = 0.;
        mmpi.allreduceSpin(&entropy, &tmp, 1, MPI_SUM);
        entropy = tmp;
    }
    return kbt * entropy;
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::computeEntropy()
{
    compute_entropy_tm_.start();

    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

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
            if (mmpi.PE0() && ct.verbose > 1)
                (*MPIdata::sout) << "computeEntropy: occupations uptodate, "
                                    "skip computation..."
                                 << std::endl;
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

    if (mmpi.PE0() && ct.verbose > 1)
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
    //    if(mmpi.PE0() && ct.verbose > 1)(*MPIdata::sout)<<"entropy =
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
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.barrier();
        if (mmpi.PE0())
            (*MPIdata::sout)
                << " CONDITION NUMBER OF THE OVERLAP MATRIX EXCEEDS TOL: "
                << rcond << "!!!" << std::endl;
        if (flag) mmpi.abort();
    }
    return rcond;
}

template <class MatrixType>
int ProjectedMatrices<MatrixType>::writeDM(HDFrestart& h5f_file)
{
    // std::cout << "ProjectedMatrices<MatrixType>::writeDM()..." << std::endl;
    std::string name("/Density_Matrix");
    return dm_->write(h5f_file, name);
}

template <class MatrixType>
int ProjectedMatrices<MatrixType>::writeSavedDM(HDFrestart& h5f_file)
{
    std::string name("/Density_Matrix_WF");

    ReplicatedWorkSpace<double>& wspace(
        ReplicatedWorkSpace<double>::instance());

    const MatrixType* matrix = mat_X_old_.get();
    wspace.initSquareMatrix(*matrix);

    DISTMATDTYPE* work_matrix = wspace.square_matrix();

    hid_t file_id = h5f_file.file_id();
    return mgmol_tools::write_matrix(file_id, name, work_matrix, dim_);
}

template <class MatrixType>
int ProjectedMatrices<MatrixType>::readDM(HDFrestart& h5f_file)
{
    std::string name("/Density_Matrix");
    return dm_->read(h5f_file, name);
}

template <class MatrixType>
int ProjectedMatrices<MatrixType>::readWFDM(HDFrestart& h5f_file)
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.barrier();
    if (mmpi.PE0()) std::cout << "ProjectedMatrices::readWFDM..." << std::endl;
    std::string name("/Density_Matrix_WF");
    return dm_->read(h5f_file, name);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::printEigenvalues(std::ostream& os) const
{
    printEigenvaluesHa(os);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::printEigenvaluesEV(std::ostream& os) const
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver
        && mmpi.instancePE0())
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
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (ct.DMEigensolver() == DMEigensolverType::Eigensolver
        && mmpi.instancePE0())
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
        energies, width, nel_, max_numst, mmpi.instancePE0(), occ);
    // if( mmpi.instancePE0() )
    //    (*MPIdata::sout)<<"computeChemicalPotentialAndOccupations() with mu="
    //        <<mu<<std::endl;

    dm_->setOccupations(occ);
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::computeLoewdinTransform(
    SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& localP,
    const int orb_index, const bool transform_matrices)
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
        Control& ct     = *(Control::instance());
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        if (mmpi.instancePE0() && ct.verbose > 1)
            std::cout << "Transform DM to reflect Loewdin orthonormalization"
                      << std::endl;
        assert(sqrtMat);

        // transform DM to reflect Loewdin orthonormalization
        dm_->transform(*sqrtMat);

        // transform matHB_ to reflect Loewdin orthonormalization
        // (we reuse sqrtMat since we are done with it)
        MatrixType& mat(*sqrtMat);
        mat.symm('r', 'l', 1., *matHB_, invSqrtMat, 0.);
        matHB_->gemm('n', 't', 1., mat, invSqrtMat, 0.);
    }

    convert_matrix(invSqrtMat, localP);
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
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& local_product)
{
    assert(gram_4dotProducts_ != nullptr);

    MatrixType ds("ds", dim_, dim_);

    convert(local_product, ds);

    work_->gemm('n', 'n', 1., ds, gram_4dotProducts_->getInverse(), 0.);

    return work_->trace();
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::dotProductWithDM(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& local_product)
{
    MatrixType ds("ds", dim_, dim_);

    convert(local_product, ds);

    work_->gemm('n', 'n', 0.5, ds, dm_->kernel4dot(), 0.);

    return work_->trace();
}

template <class MatrixType>
double ProjectedMatrices<MatrixType>::dotProductSimple(
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& local_product)
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
    const SquareLocalMatrices<MATDTYPE, MemorySpace::Host>& mat)
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
    const int order, const double emin, const double emax)
{
    assert(emax > emin);
    assert(nel_ >= 0.);

    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    // create Chebyshev approximation object
    ChebyshevApproximation<MatrixType> chebapp(emin, emax, order, this);

    if (mmpi.instancePE0() && ct.verbose > 0)
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

    if (mmpi.instancePE0() && ct.verbose > 0)
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
            mmpi.abort();
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
            if (mmpi.instancePE0())
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
    double orbital_occupation = mmpi.nspin() > 1 ? 1. : 2.;
    dm.scal(orbital_occupation);
    dm_->setMatrix(dm);

    return mu_;
}

/* Use the power method to compute the extents of the spectrum of the
 * generalized eigenproblem.
 */
#ifdef MGMOL_USE_SCALAPACK
template <>
void ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>::
    computeGenEigenInterval(
        std::vector<double>& interval, const int maxits, const double pad)
{
    dist_matrix::DistMatrix<DISTMATDTYPE> mat(*matHB_);

    static PowerGen<dist_matrix::DistMatrix<DISTMATDTYPE>,
        dist_matrix::DistVector<DISTMATDTYPE>>
        power(dim_);

    power.computeGenEigenInterval(mat, *gm_, interval, maxits, pad);
}
#endif

template <>
void ProjectedMatrices<ReplicatedMatrix>::computeGenEigenInterval(
    std::vector<double>& interval, const int maxits, const double pad)
{
    ReplicatedMatrix mat(*matHB_);

    static PowerGen<ReplicatedMatrix, ReplicatedVector> power(dim_);

    power.computeGenEigenInterval(mat, *gm_, interval, maxits, pad);
}

#ifdef MGMOL_USE_SCALAPACK
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
#endif

template <>
void ProjectedMatrices<ReplicatedMatrix>::consolidateH()
{
    consolidate_H_tm_.start();

    // assign SquareLocalMatrices to matH_
    matH_->assign(*localHl_);
    matH_->add(*localHnl_);

    // sum up across MPI tasks
    matH_->consolidate();

    consolidate_H_tm_.stop();
}

template <class MatrixType>
void ProjectedMatrices<MatrixType>::updateSubMatX(const MatrixType& dm)
{
    convert_matrix(dm, *localX_);
}

#ifdef MGMOL_USE_SCALAPACK
template <>
SquareLocalMatrices<double, MemorySpace::Host>
ProjectedMatrices<dist_matrix::DistMatrix<double>>::getReplicatedDM()
{
    SquareLocalMatrices<double, MemorySpace::Host> sldm(1, dim_);
    const dist_matrix::DistMatrix<double>& dm(dm_->getMatrix());
    dm.allgather(sldm.getRawPtr(), dim_);

    return sldm;
}
#endif

template <>
SquareLocalMatrices<double, MemorySpace::Host>
ProjectedMatrices<ReplicatedMatrix>::getReplicatedDM()
{
    SquareLocalMatrices<double, MemorySpace::Host> sldm(1, dim_);
    const ReplicatedMatrix& dm(dm_->getMatrix());
    dm.get(sldm.getRawPtr(), dim_);

    return sldm;
}

#ifdef MGMOL_USE_SCALAPACK
template class ProjectedMatrices<dist_matrix::DistMatrix<DISTMATDTYPE>>;
#endif
template class ProjectedMatrices<ReplicatedMatrix>;
