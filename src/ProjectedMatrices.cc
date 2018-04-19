// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. 
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved. 
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id: ProjectedMatrices.cc 1750 2017-04-20 00:01:44Z jeanluc $
#include "ProjectedMatrices.h"

#if PROCRUSTES
#include "munkres.h"
#endif

#include "DensityMatrix.h"
#include "GramMatrix.h"
#include "ReplicatedWorkSpace.h"
#include "MPIdata.h"
#include "HDFrestart.h"
#include "Control.h"
#include "MGmol_MPI.h"
#include "SubMatrices.h"
#include "SparseDistMatrix.h"
#include "fermi.h"
#include "tools.h"
#include "random.h"

#include <iomanip>
#include <fstream>

#define RY2EV                13.605804

Timer ProjectedMatrices::sygv_tm_("ProjectedMatrices::sygv");
Timer ProjectedMatrices::compute_inverse_tm_("ProjectedMatrices::computeInverse");
Timer ProjectedMatrices::compute_invB_tm_("ProjectedMatrices::computeInvB");
Timer ProjectedMatrices::init_gram_matrix_tm_("ProjectedMatrices::initialize_Gram_Matrix");
Timer ProjectedMatrices::update_theta_tm_("ProjectedMatrices::updateTheta");
Timer ProjectedMatrices::update_submatT_tm_("ProjectedMatrices::updateSubmatT");
Timer ProjectedMatrices::update_submatX_tm_("ProjectedMatrices::updateSubmatX");
Timer ProjectedMatrices::eigsum_tm_("ProjectedMatrices::eigsum");
Timer ProjectedMatrices::consolidate_H_tm_("ProjectedMatrices::consolidate_sH");
Timer ProjectedMatrices::eig_interval_tm_("ProjectedMatrices::computeEigenInterval");

//const short LocGridOrbitals::sparse_distmatrix_nb_tasks_per_partitions_=96;
//const short LocGridOrbitals::sparse_distmatrix_nb_tasks_per_partitions_=256;
short ProjectedMatrices::sparse_distmatrix_nb_tasks_per_partitions_=256;
short ProjectedMatrices::n_instances_=0;

dist_matrix::RemoteTasksDistMatrix<DISTMATDTYPE>* ProjectedMatrices::remote_tasks_DistMatrix_=0;
GramMatrix* ProjectedMatrices::gram_4dotProducts_=0;
DensityMatrix* ProjectedMatrices::dm_4dot_product_=0;

static int sparse_distmatrix_nb_partitions=128;

ProjectedMatrices::ProjectedMatrices(const int ndim, 
                                     const bool with_spin)
                                     : with_spin_(with_spin)
{
    dim_  =ndim;
    width_=0.;
    min_val_=0.25;
    submat_indexing_=0;
 
    dm_    =new DensityMatrix(ndim);
    gm_    =new GramMatrix(ndim);
    mat_X_old_=0;
    mat_L_old_=0;
        
#ifdef USE_DIS_MAT
    submatWork_=0;
    submatLS_=0;
#endif
    localX_=0;
    localT_=0;
    
    sH_=0;
   
    if( dim_>0 )
    {
        eigenvalues_.resize(dim_);
    }

    if( dim_>0 )
    {
        MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
        const dist_matrix::BlacsContext& bc=*mbc.bcxt();
        matH_  =new dist_matrix::DistMatrix<DISTMATDTYPE>("H",     bc, ndim, ndim);
        
        matHB_ =new dist_matrix::DistMatrix<DISTMATDTYPE>("HB",    bc, ndim, ndim);
        theta_ =new dist_matrix::DistMatrix<DISTMATDTYPE>("Theta", bc, ndim, ndim);
        u_     =new dist_matrix::DistMatrix<DISTMATDTYPE>("U",     bc, ndim, ndim);
        work_  =new dist_matrix::DistMatrix<DISTMATDTYPE>("work",  bc, ndim, ndim);
#if PROCRUSTES
        occupation0_.resize(dim_);
        u0_    =0;
        p_     =new dist_matrix::DistMatrix<DISTMATDTYPE>("P",     bc, ndim, ndim);

        u0_->identity();
#endif
    }

    if(onpe0)(*MPIdata::sout)<<"ProjectedMatrices: sparse_distmatrix_nb_tasks_per_partitions_="
                             <<sparse_distmatrix_nb_tasks_per_partitions_<<endl;
    n_instances_++;
}
    
ProjectedMatrices::~ProjectedMatrices()
{
    assert( dm_!=0 );
    assert( gm_!=0 );
    
    delete dm_; dm_=0;
    delete gm_; gm_=0;
    if( mat_X_old_ )
    {
        delete mat_X_old_;
        mat_X_old_=0;
    }
    if( mat_L_old_ )
    {
        delete mat_L_old_;
        mat_L_old_=0;
    }

    if( dim_>0 )
    {
        assert( matH_!=0 );
        assert( theta_!=0 );
        assert( u_!=0 );
        assert( work_!=0 );
    
        delete matHB_;matHB_=0;
        delete matH_;
        delete theta_;theta_=0;
        delete u_;
        delete work_;
        
        delete localX_;localX_=0;
        delete localT_;localT_=0;
        
        delete sH_;sH_=0;
        
        delete submatWork_;submatWork_=0;
        delete submatLS_;submatLS_=0;

#if PROCRUSTES
        delete u0_;
        delete p_;
#endif
        delete submat_indexing_;submat_indexing_=0;
    }

    if( n_instances_==1 ){
        if( gram_4dotProducts_!=0 )
        {
            delete gram_4dotProducts_;
            gram_4dotProducts_=0;
        }
    }

    n_instances_--;
};

void ProjectedMatrices::setup(const double kbt, const int nel, 
                              const vector<vector<int> >& global_indexes)
{
    assert( global_indexes.size()>0 );
    
    ProjectedMatricesInterface::setup(kbt,nel,global_indexes);
    
    global_indexes_=global_indexes;
    
#ifdef USE_DIS_MAT
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm=mmpi.commSpin();
    if( submat_indexing_!=0 ) delete submat_indexing_;
    submat_indexing_=new dist_matrix::SubMatricesIndexing<DISTMATDTYPE>(
        global_indexes,comm,gm_->getMatrix());
    
    if( localX_!=NULL )delete localX_;
    if( localT_!=NULL )delete localT_;

    if( sH_!=NULL )delete sH_;

    if( submatWork_ !=NULL )delete submatWork_;
    if( submatLS_!=NULL )delete submatLS_;
    //if( orbitals_type_!=0 )
    {
        
        if( dim_>0 )
        {
            submatWork_=new dist_matrix::SubMatrices<DISTMATDTYPE>("Work",global_indexes,
                                         comm,
                                         dm_->getMatrix(),
                                         *submat_indexing_);
            submatLS_=new dist_matrix::SubMatrices<DISTMATDTYPE>("LS",global_indexes,
                                         comm,
                                         gm_->getMatrix(),
                                         *submat_indexing_);
                                         
            localX_=new SquareLocalMatrices<MATDTYPE>(subdiv_,chromatic_number_);
            localT_=new SquareLocalMatrices<MATDTYPE>(subdiv_,chromatic_number_);
        }
    }
#endif

#ifdef USE_MPI
    sH_=new dist_matrix::SparseDistMatrix<DISTMATDTYPE>(comm,*matH_,
                                remote_tasks_DistMatrix_,
                                sparse_distmatrix_nb_partitions);
#else
    sH_=new dist_matrix::SparseDistMatrix<DISTMATDTYPE>(0,*matH_,remote_tasks_DistMatrix_);
#endif
}

#if PROCRUSTES
void procrustes(dist_matrix::DistMatrix<DISTMATDTYPE>& a, dist_matrix::DistMatrix<DISTMATDTYPE>& b,
                dist_matrix::DistMatrix<DISTMATDTYPE>& p)
{
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = mbc. bcxt();
    const int nst    =a.m();

    p.gemm('t','n',1.,b,a,0.);
    
    munkres::Matrix<double> matrix(nst, nst);
    
    double maxval=0.;
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            maxval = max( maxval, p.val(row+nst*col) );
        }
    }
    
    // Initialize matrix
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            matrix(row,col) = maxval-fabs( p.val(row+nst*col) );
        }
    }

    // Apply Munkres algorithm to matrix.
    munkres::Munkres m;
    m.solve(matrix);

    // 0 -> 1, -1 -> 0
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            matrix(row,col)+=1.;
        }
    }
    
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            if( matrix(row,col)>0.5 ){
                if( p.val(row+nst*col)<0. )
                    matrix(row,col)=-1.;
                break;
            }
        }
    }
    
    for ( int row = 0 ; row < nst ; row++ ) {
        for ( int col = 0 ; col < nst ; col++ ) {
            p.setval(row+nst*col, matrix(row,col));
        }
    }

    //(*MPIdata::sout)<<"Procrustes matrix..."<<endl;
    //p.print((*MPIdata::sout),0,0,nst,nst); 

    //w.gemm('n','n',1.,b,p,0.);
}
#endif

       
void ProjectedMatrices::computeInvS()
{
    compute_inverse_tm_.start();
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        (*MPIdata::sout)<<"ProjectedMatrices::computeInvS()"<<endl;
#endif
    gm_->computeInverse();
    compute_inverse_tm_.stop();
}

void ProjectedMatrices::rotateAll(const dist_matrix::DistMatrix<DISTMATDTYPE>&  rotation_matrix,
                                  const bool flag_eigen)
{
    // S -> U^T S U
    // rotate overlap and l_s
    if( flag_eigen ){
        gm_->set2Id(-1);
    }else{
        gm_->rotateAll(rotation_matrix);
    }
    //(*MPIdata::sout)<<"matS"<<endl;
    //matS_->print((*MPIdata::sout),0,0,5,5);

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
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();
    dist_matrix::DistMatrix<DISTMATDTYPE> pmatrix("pmatrix",bc,dim_,dim_);
 
    mat.fillDistMatrix(pmatrix,global_indexes_);

    gm_->applyInv(pmatrix);
    
    mat.init(pmatrix, global_indexes_, *submat_indexing_);    
}

void ProjectedMatrices::setDMto2InvS()
{
    dm_->setto2InvS(gm_->getInverse(),gm_->getAssociatedOrbitalsIndex());
}

void ProjectedMatrices::solveGenEigenProblem(dist_matrix::DistMatrix<DISTMATDTYPE>& z,
                                             vector<DISTMATDTYPE>& val, char job)
{
    assert( val.size()==eigenvalues_.size() );
    
    sygv_tm_.start();

    dist_matrix::DistMatrix<DISTMATDTYPE>  mat(*matHB_);
    
    // Transform the generalized eigenvalue problem to a standard form
    gm_->sygst(mat);
 
    // solve a standard symmetric eigenvalue problem
    mat.syev(job, 'l', eigenvalues_, *u_);

#if PROCRUSTES
    if( u0_!=0 ){
        work_->gemm('t','n',1.,*u_,*u0_,0.);
        for(int j=0;j<dim_;j++){
            double val=0.;
            int i=work_->iamax(j,val);
            if( val<0. )u_->scalColumn(j,-1.);
        }
        procrustes(*u_,*u0_,*p_);
    }
#endif
 
    // Get the eigenvectors Z of the generalized eigenvalue problem
    // Solve Z=L**(-T)*U 
    z = *u_;
    gm_->solveLST(z);

    val = eigenvalues_;

    sygv_tm_.stop();
}

void ProjectedMatrices::buildDM(const dist_matrix::DistMatrix<DISTMATDTYPE>& z, const int orbitals_index)
{
    dm_->build(z,orbitals_index);
}

void ProjectedMatrices::buildDM(const dist_matrix::DistMatrix<DISTMATDTYPE>& z,
                                const vector<DISTMATDTYPE>& occ,
                                const int orbitals_index)
{
    dm_->build(z, occ, orbitals_index);
}
void ProjectedMatrices::buildDM(const vector<DISTMATDTYPE>& occ,const int orbitals_index)
{
    dm_->build(occ,orbitals_index);
}

void ProjectedMatrices::updateDMwithEigenstates(const int iterative_index)
{
    if( onpe0 )(*MPIdata::sout)<<"ProjectedMatrices: Compute DM using eigenstates"<<endl;
    
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();
 
    dist_matrix::DistMatrix<DISTMATDTYPE>  zz("Z",bc, dim_, dim_);
    vector<DISTMATDTYPE> val(dim_);

    // solves generalized eigenvalue problem
    // and return solution in zz and val
    solveGenEigenProblem(zz,val);
    setAuxilliaryEnergiesFromEigenenergies();
    double final_mu = computeChemicalPotentialAndOccupations();
    if(onpe0)cout<<"Final mu_ = "<<final_mu<<endl;
    
    // Build the density matrix X 
    // X = Z * gamma * Z^T
    buildDM(zz,iterative_index); 
}

void ProjectedMatrices::updateDM(const int iterative_index)
{
    Control& ct = *(Control::instance());
    
    if(ct.dm_algo == 0)
       updateDMwithEigenstates(iterative_index); 
}

void ProjectedMatrices::updateDMwithEigenstatesAndRotate(const int iterative_index,
                                                         dist_matrix::DistMatrix<DISTMATDTYPE>&  zz)
{
    vector<DISTMATDTYPE> val(dim_);

    // solves generalized eigenvalue problem
    // and return solution in zz and val
    solveGenEigenProblem(zz,val);
    setAuxilliaryEnergiesFromEigenenergies();
    computeChemicalPotentialAndOccupations();
    
    rotateAll(zz,true);
    
    dm_->build(iterative_index); 
}

#if 0
void ProjectedMatrices::rotateBackDM()
{
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();
 
    computeOccupations(occupation0_);
    *u0_=*work_;

#if PROCRUSTES
    procrustes(*u_,*u0_,*p_);
#endif

    dist_matrix::DistMatrix<DISTMATDTYPE> rot("R",bc, dim_, dim_);
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
    dist_matrix::DistMatrix<DISTMATDTYPE> test("t",bc,nst,nst);
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

    dist_matrix::DistMatrix<DISTMATDTYPE> gamma("Gamma",bc, &occupation_[0], dim_, dim_);
    const double orbital_occupation = with_spin_ ? 1. : 2.;
    gamma.scal(orbital_occupation); // rescale for spin
 
    work_->symm('r', 'l', 1., gamma, z, 0.);
    dm_->gemm('n', 't', 1., *work_, z, 0.);
}
#endif

void ProjectedMatrices::computeOccupationsFromDM()
{
    Control& ct = *(Control::instance());
    if( ct.dm_algo==0 && dim_>0 ){
        assert( dm_!=0 );
        dm_->computeOccupations(gm_->getCholeskyL());    
    }
}

void ProjectedMatrices::getOccupations(vector<DISTMATDTYPE>& occ)const
{
    dm_->getOccupations(occ);
}
void ProjectedMatrices::setOccupations(const vector<DISTMATDTYPE>& occ)
{
    dm_->setOccupations(occ);
}
void ProjectedMatrices::printDM(ostream& os)const
{
    dm_->print(os);
}

const dist_matrix::DistMatrix<DISTMATDTYPE>& ProjectedMatrices::dm()const
{
    assert( dm_!=0 );
    return dm_->getMatrix();
}
const dist_matrix::DistMatrix<DISTMATDTYPE>& ProjectedMatrices::kernel4dot()const
{
    return dm_->kernel4dot();
}

double ProjectedMatrices::getNel()const
{
    double val=dm_->dot(gm_->getMatrix());
    if( with_spin_ ){
        double tmp=0.;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduceSpin(&val, &tmp, 1, MPI_SUM);
        val=tmp;
    }
    
    return val;
}
    
double ProjectedMatrices::getEigSum()
{
    if( dim_==0 )return 0.;
    
    eigsum_tm_.start();
    work_->symm('l','l',1.,*matHB_,gm_->getInverse(),0.);
    
    // return sum in Ry
    double val = work_->trace();
    eigsum_tm_.stop();
    
    return val;
}

double ProjectedMatrices::getExpectationH()
{
    if( dim_==0 )return 0.;
    
    return getExpectation(*matHB_);
}

double ProjectedMatrices::getExpectation(const dist_matrix::DistMatrix<DISTMATDTYPE>& A)
{
    work_->gemm('n', 'n', 1., A, dm_->getMatrix(), 0.);
    double val=work_->trace();
    
    if( with_spin_ ){
        double tmp=0.;
        MGmol_MPI& mmpi = *(MGmol_MPI::instance());
        mmpi.allreduceSpin(&val, &tmp, 1, MPI_SUM);
        val=tmp;
    }
    
    return val;
}

// strip dm from the overlap contribution
// dm <- Ls**T * dm * Ls
void ProjectedMatrices::stripDM()
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        cout<<"ProjectedMatrices::stripDM()"<<endl;
#endif
#ifdef DEBUG // TEST
    double dd=dm_->getMatrix().trace();
    if( onpe0 )
        cout<<"test:  Trace DM = "<<dd<<endl;
    if( dm_->getMatrix().active() )assert( dd>0. );
#endif
    dm_->stripS(gm_->getCholeskyL(), gm_->getAssociatedOrbitalsIndex());
}

void ProjectedMatrices::dressupDM()
{
#ifdef PRINT_OPERATIONS
    if( onpe0 )
        cout<<"ProjectedMatrices::dressupDM()"<<endl;
#endif
    dm_->dressUpS(gm_->getCholeskyL(), gm_->getAssociatedOrbitalsIndex());
}

double ProjectedMatrices::computeEntropy(const double kbt)
{
    return kbt*dm_->computeEntropy();
}

double ProjectedMatrices::computeEntropy()
{
    //if(onpe0)(*MPIdata::sout)<<"ProjectedMatrices::computeEntropy()"<<endl;
    //if(onpe0)(*MPIdata::sout)<<"width_="<<width_<<endl;

    Control& ct = *(Control::instance());
    double entropy;

    if(ct.dm_algo==0 || dm_->fromUniformOccupations())
    {    
       if(!occupationsUptodate() )
       {
          computeOccupationsFromDM();
       }
       else
       {
           if( onpe0 && ct.verbose>1 )
               (*MPIdata::sout)<<"occupations uptodate, skip computation..."<<endl;
       }
       entropy = computeEntropy(width_);
     }
    return entropy;
}


void ProjectedMatrices::printOccupations(ostream& os)const
{
    if(dm_->occupationsUptodate())dm_->printOccupations(os);
}

double ProjectedMatrices::checkCond(const double tol, const bool flag)
{
    double rcond = computeCond();
    
    if( rcond>tol ){
        //ofstream tfile("s.mm", ios::out);
        //gm_->printMM(tfile);
        //tfile.close();
#ifdef USE_MPI
        MGmol_MPI& mgmolmpi = *(MGmol_MPI::instance());
        mgmolmpi.barrier();
#endif
        if( onpe0 )
            (*MPIdata::sout)<<" CONDITION NUMBER OF THE OVERLAP MATRIX EXCEEDS TOL: "
                <<rcond<<"!!!"<<endl;
        Control& ct = *(Control::instance());    
        if( flag )ct.global_exit(2);
    }
    return rcond;
}
////// TEMPLATE THIS FOR FLOAT OPTION ??
int ProjectedMatrices::writeDM_hdf5(HDFrestart& h5f_file)
{
    hid_t file_id=h5f_file.file_id();
    
    ReplicatedWorkSpace& wspace( ReplicatedWorkSpace::instance() );

    wspace.initSquareMatrix(dm_->getMatrix());

    if( file_id<0 )return 0;

    hsize_t dims[2]={dim_, dim_};   

    // filespace identifier
    hid_t  dataspace = H5Screate_simple(2, dims, NULL); 

    hid_t  dset_id = H5Dcreate2(file_id, "/Density_Matrix", H5T_NATIVE_DOUBLE, 
                        dataspace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if( dset_id<0 ){
        (*MPIdata::serr)<<"ProjectedMatrices::write_dm_hdf5: H5Dcreate2 failed!!!"<<endl;
        return -1;
    }

    hid_t memspace=dataspace;
    hid_t filespace=dataspace;

    DISTMATDTYPE* work_matrix=wspace.square_matrix();
    herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                             H5P_DEFAULT, work_matrix);
    if( status<0 )
    {
        (*MPIdata::serr)<<"Orbitals: H5Dwrite failed!!!"<<endl;
        return -1;
    }
    
    status = H5Dclose(dset_id);
    if( status<0 )
    {
        (*MPIdata::serr)<<"ProjectedMatrices::write_dm_hdf5(), H5Dclose failed!!!"<<endl;
        return -1;
    }
    status = H5Sclose(dataspace);
    if( status<0 )
    {
        (*MPIdata::serr)<<"ProjectedMatrices::write_dm_hdf5(), H5Sclose failed!!!"<<endl;
        return -1;
    }

    return 0;
}
////// TEMPLATE THIS FOR FLOAT OPTION ??
int ProjectedMatrices::read_dm_hdf5(hid_t file_id)
{
    ReplicatedWorkSpace& wspace( ReplicatedWorkSpace::instance() );
    DISTMATDTYPE* work_matrix=wspace.square_matrix();

    int ierr=0;
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if( mmpi.instancePE0() )
    {
        hid_t  dset_id = H5Dopen2(file_id, "/Density_Matrix", H5P_DEFAULT);
        if( dset_id<0 )
        {
            (*MPIdata::serr)<<"H5Dopen failed for /Density_Matrix!!!"<<endl;
        }
        else
        {
            ierr=1;
            herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, work_matrix);
            if( status<0 ){
                (*MPIdata::serr)<<"H5Dread failed for /Density_Matrix!!!"<<endl;
                return -1;
            }
            
            status = H5Dclose(dset_id);
            if( status<0 ){
                (*MPIdata::serr)<<"H5Dclose failed!!!"<<endl;
            return -1;
            }
        }
    }
    mmpi.bcast(&ierr, 1);
    if( ierr>=0 )
        wspace.mpiBcastSquareMatrix();
    if( ierr>=0 )dm_->initMatrix(work_matrix);

    return ierr;
}

void ProjectedMatrices::printEigenvalues(ostream& os)const
{
    printEigenvaluesHa(os);
}

void ProjectedMatrices::printEigenvaluesEV(ostream& os)const
{
    Control& ct = *(Control::instance());
    if( ct.dm_algo==0 && onpe0 )
    {
        os<<endl<<" Eigenvalues [eV]:";
        
        // Print ten to a row.
        os.setf(ios::right,ios::adjustfield);
        os.setf(ios::fixed,ios::floatfield);
        os<<setprecision(3);
        for(int i = 0;i < dim_;i++)
        {
            if ( (i%10) == 0 ) os << endl;
            os<<setw(7)<<RY2EV*eigenvalues_[i]<<" ";
        }
        os<<endl;

        if(width_>1.e-10)
            os<<" FERMI ENERGY   = "<<RY2EV * mu_ <<"[eV]"<<endl;

    }
}

void ProjectedMatrices::printEigenvaluesHa(ostream& os)const
{
    Control& ct = *(Control::instance());
    if( ct.dm_algo==0 && onpe0 )
    {
        os<<endl<<" Eigenvalues [Ha]:";
        
        // Print ten to a row.
        os.setf(ios::right,ios::adjustfield);
        os.setf(ios::fixed,ios::floatfield);
        os<<setprecision(3);
        for(int i = 0;i < dim_;i++)
        {
            if ( (i%10) == 0 ) os << endl;
            os<<setw(7)<<0.5*eigenvalues_[i]<<" ";
        }
        os<<endl;

        if(width_>1.e-10)
            os<<" FERMI ENERGY   = "<<0.5 * mu_ <<"[Ha]"<<endl;

    }
}

void ProjectedMatrices::setAuxilliaryEnergiesFromEigenenergies()
{
    Control& ct = *(Control::instance());
    if( ct.dm_algo==0 )
    {
       assert( (int)eigenvalues_.size()==dim_ );

       if( aux_energies_.size()==0 )
           aux_energies_.resize(dim_);
       aux_energies_=eigenvalues_;

       assert( (int)aux_energies_.size()==dim_ );
    }
}

void ProjectedMatrices::setAuxilliaryEnergiesFromOccupations()
{
    if( aux_energies_.size()==0 )
        aux_energies_.resize(dim_);

    const double eps=1.e-15;
    const double shift = width_>1.e-8 ? -log(eps)/width_ : 1000.;
    
    vector<DISTMATDTYPE> occ(dim_);
    dm_->getOccupations(occ);
    for(int i=0;i<dim_;i++){
        if( occ[i]<eps){
            aux_energies_[i]=mu_+shift;
        }else if( (1.-occ[i])<eps){
            aux_energies_[i]=mu_-shift;
        }else{
            aux_energies_[i]=mu_+width_*log(1./occ[i]-1.);
        }
    }
    assert( (int)aux_energies_.size()==dim_ );
}

// find the Fermi level by a bisection
// algorithm adapted from numerical recipes, 2nd edition
// and fill orbitals accordingly (in fermi_distribution)
double ProjectedMatrices::computeChemicalPotentialAndOccupations(const std::vector<DISTMATDTYPE>& energies,
                                                                 const double width,
                                                                 const int nel,
                                                                 const int max_numst)
{
    assert( energies.size()>0 );
    assert( nel>=0 );
    
    Control& ct = *(Control::instance());
    if( onpe0 && ct.verbose>0 )
        (*MPIdata::sout)<<"computeChemicalPotentialAndOccupations() with width="
                        <<width<<", for "<<nel<<" electrons"<<endl;

    const int    maxit = 100;
    const double charge_tol = 1.0e-12;

    const double emin = *std::min_element(energies.begin(),energies.end());
    const double emax = *std::max_element(energies.begin(),energies.end());

    double  mu1=emin     -0.001;
    double  mu2=emax     +10.*width;
    assert(mu1<mu2);
    bool done=false;

    if( nel<=0 )
    {
        //if( onpe0 )
        //    (*MPIdata::sout)<<"computeChemicalPotentialAndOccupations() with nel="<<nel<<endl;
        mu1=-10000.;
        mu2= 10000.;
    }
    
    vector<DISTMATDTYPE> occ(dim_,0.);

    if( 2*dim_<=nel )
    {
        done=true;
        mu_=mu2;
        for(int i=0;i<dim_;i++)occ[i]=1.;
    }

    double f2=0.;
    if( !done )
    {
        f2 = 2.*fermi_distribution(mu2,max_numst,width,energies,occ)-(double)nel;
        // no unoccupied states
        if(fabs(f2)<charge_tol)
        {
            done=true;
            mu_=mu2;
        }
    }
    double f1=0.;
    if( !done )
    {
        f1 = 2.*fermi_distribution(mu1,max_numst,width,energies,occ)-(double)nel;
        if(fabs(f1)<charge_tol)
        {
            if( onpe0 )(*MPIdata::sout)<<"only unoccupied states"<<endl;
            done=true;
            mu_=mu1;
        }
    }

    if( !done )
    {
        if (f1*f2 > 0.)
        {
            (*MPIdata::sout)<<"ERROR: mu1="<<mu1<<", mu2="<<mu2<<endl;
            (*MPIdata::sout)<<"ERROR: f1="<<f1<<", f2="<<f2<<endl;
            (*MPIdata::sout)<<"nel="<<nel<<", width="<<width<<endl;
            Control& ct = *(Control::instance());    
            ct.global_exit(2);
        }

        double  dmu;
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
        int  iter = 0;
        double f =0.;
        do
        {
            iter++;

            dmu *= 0.5;
            mu1 = mu_ + dmu;
            f = 2.*fermi_distribution(mu1,max_numst,width,energies,occ)-(double)nel;

            if (f <= 0.) {
                mu_ = mu1;
                f = - f;
            }

        } while ( (iter < maxit) && (f > charge_tol) );
        
        if ( f > charge_tol )
        {
            if( onpe0 )
            {
                (*MPIdata::sout)<<"WARNING: ProjectedMatrices::computeChemicalPotentialAndOccupations()"<<endl;
                (*MPIdata::sout)<<"Iterations did not converge to tolerance "<<scientific<<charge_tol<<endl;
                (*MPIdata::sout)<<"f= "<<f<<", nel="<<nel<<", max_numst="<<max_numst<<endl;
            }
        }        
    }

    //if( onpe0 )
    //    (*MPIdata::sout)<<"computeChemicalPotentialAndOccupations() with mu="
    //        <<mu_<<endl;

    dm_->setOccupations(occ);
    
    return mu_;
}

void ProjectedMatrices::getLoewdinTransform(SquareLocalMatrices<MATDTYPE>& localP)
{
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();

    dist_matrix::DistMatrix<DISTMATDTYPE>  mat(gm_->getMatrix());
    dist_matrix::DistMatrix<DISTMATDTYPE>  vect("eigenvectors",bc,dim_,dim_);
    vector<DISTMATDTYPE> eigenvalues(dim_);
    mat.syev('v', 'l', eigenvalues, vect);

    vector<DISTMATDTYPE> diag_values(dim_);
    for(int i=0;i<dim_;i++)diag_values[i]=(DISTMATDTYPE)(1./sqrt(eigenvalues[i]));

    dist_matrix::DistMatrix<DISTMATDTYPE>  matP("P",bc,dim_,dim_);
    matP.clear();
    matP.setDiagonal(diag_values);
    //if( onpe0 )
    //  (*MPIdata::sout)<<"LocGridOrbitals::orthonormalizeLoewdin() --- matP"<<endl;
    //matP.print((*MPIdata::sout),0,0,5,5);
    mat.symm('r', 'l', 1., matP, vect, 0.);
    matP.gemm('n', 't', 1., mat, vect, 0.);
    
    localP.init(matP, global_indexes_, *submat_indexing_);
}

double ProjectedMatrices::getTraceDiagProductWithInvS(vector<DISTMATDTYPE>& ddiag)
{
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();

    dist_matrix::DistMatrix<DISTMATDTYPE> diag("diag", bc, dim_, dim_);    
    diag.setDiagonal(ddiag);

    work_->gemm('n', 'n', 1., diag, gm_->getInverse(), 0.);
    
    return work_->trace();    
}

void ProjectedMatrices::resetDotProductMatrices()
{
    if( gram_4dotProducts_!=0 )delete gram_4dotProducts_;
    gram_4dotProducts_=new GramMatrix(*gm_);
}

double ProjectedMatrices::dotProductWithInvS(const SquareLocalMatrices<MATDTYPE>& local_product)
{
    assert( gram_4dotProducts_!=0 );
    
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();
    
    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", bc, dim_, dim_);
    local_product.fillDistMatrix(ds,global_indexes_);

    dist_matrix::DistMatrix<DISTMATDTYPE> work("work", bc, dim_, dim_);
    work.gemm('n', 'n', 1., ds, gram_4dotProducts_->getInverse(), 0.);

    return work.trace();
}

double ProjectedMatrices::dotProductWithDM(const SquareLocalMatrices<MATDTYPE>& local_product)
{
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();
    
    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", bc, dim_, dim_);
    local_product.fillDistMatrix(ds,global_indexes_);

    dist_matrix::DistMatrix<DISTMATDTYPE> work("work", bc, dim_, dim_);
    work.gemm('n', 'n', 0.5, ds, kernel4dot(), 0.);

    return work.trace();
}

double ProjectedMatrices::dotProductSimple(const SquareLocalMatrices<MATDTYPE>& local_product)
{
    assert( dm_4dot_product_!=0 );
    
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();
    
    dist_matrix::DistMatrix<DISTMATDTYPE> ds("ds", bc, dim_, dim_);
    local_product.fillDistMatrix(ds,global_indexes_);

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

//Assumes SquareLocalMatrix object contains partial contributions
double ProjectedMatrices::computeTraceInvSmultMat(const SquareLocalMatrices<MATDTYPE>& mat)
{
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();
    dist_matrix::DistMatrix<DISTMATDTYPE> pmatrix("pmatrix",bc,dim_,dim_);
 
    mat.fillDistMatrix(pmatrix,global_indexes_);

    gm_->applyInv(pmatrix);
    return pmatrix.trace();
/*
    dist_matrix::DistMatrix<DISTMATDTYPE> work_matrix("work_matrix",bc,dim_,dim_); 
    work_matrix.symm('l', 'l', 1., gm_->getInverse(), pmatrix, 0.);
    
    dist_matrix::DistMatrix<DISTMATDTYPE> pm = gm_->getInverse();//("pmatrix",bc,dim_,dim_);
    double itrace = pm.trace();
    if(onpe0)cout<<"INVERSE trace = "<<itrace<<endl;
    
    return work_matrix.trace(); 
*/
}

double ProjectedMatrices::computeTraceInvSmultMatMultTheta(const dist_matrix::DistMatrix<DISTMATDTYPE>& mat)
{
    assert(theta_ != 0);
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();
    dist_matrix::DistMatrix<DISTMATDTYPE> pmat("pmat",bc,dim_,dim_);    
    
    //compute mat*theta_
    pmat.gemm('n', 'n', 1.0, mat, *theta_, 0.); 
    
    //compute invS*pmat = invS*(mat*theta)
    gm_->applyInv(pmat);
    
    return pmat.trace(); 
}


/* Use the power method to compute the extents of the spectrum of the generalized eigenproblem.
 * In order to use a residual-based convergence criterion in an efficient way, we delay 
 * normalization of the vectors to avoid multiple matvecs.
 * NOTE: We are only interested in the eigenvalues, so the final eigenvector may not be normalized.
*/

void ProjectedMatrices::computeGenEigenInterval(std::vector<double>& interval, const int maxits, const double pad)
{   
    srand(13579);

    eig_interval_tm_.start();
    
    interval.clear();
    
    MatricesBlacsContext& mbc( MatricesBlacsContext::instance() );
    const dist_matrix::BlacsContext& bc = *mbc. bcxt();    

    dist_matrix::DistMatrix<DISTMATDTYPE>  mat(*matHB_);
    dist_matrix::DistMatrix<DISTMATDTYPE>  smat(gm_->getMatrix());
    
    // use the power method to get the eigenvalue interval
    const int m = mat.m(); // number of global rows
    const int mloc = mat.mloc(); // number of local rows
    const double one = 1., zero = 0.;

    // define static variables and shift for matrices to speed up convergence
    static double shft = 0.; //shft is initially zero
    // shift
    mat.axpy(shft, gm_->getMatrix());     

    // initialize random vectors for power method
    static std::vector<DISTMATDTYPE>vec1(generate_rand(mloc)); // initial random vector.
    static std::vector<DISTMATDTYPE>vec2(vec1); // initial random vector.

   // initialize solution data
   // initial guess
   dist_matrix::DistMatrix<DISTMATDTYPE> sol("sol",bc, m, 1);
   sol.assignColumn(&vec1[0], 0); // initialize local solution data
   // new solution
   dist_matrix::DistMatrix<DISTMATDTYPE> new_sol("new_sol",bc, m, 1);   
   std::vector<DISTMATDTYPE>vec(mloc,0.);
   new_sol.assignColumn(&vec[0], 0);
   
   // get norm of initial sol
   double alpha = sol.nrm2();
   double gamma = 1./alpha;
   if(onpe0)cout<<"e1:: ITER 0:: = "<< alpha<<" shft = "<<shft<<endl;   
      
   // residual
   dist_matrix::DistMatrix<DISTMATDTYPE> res(new_sol);
   // initial eigenvalue estimate (for shifted system)
   double beta = sol.dot(new_sol);   

   // compute first extent
   int iter1 = 0;
   // loop
   for(int i=0; i<maxits; i++)
   {
      iter1++;
      
      // First compute residual for convergence check
      res.gemv('N', one, mat, sol, zero);
      // store matvec result appropriately scaled for later reuse
      new_sol.clear();
      new_sol.axpy(gamma,res);       
      // Compute residual: res = beta*S*x - mat*x
      res.gemm('N', 'N', beta, gm_->getMatrix(), sol, -1.);      
      // compute residual norm
      double resnorm = res.nrm2();
      // check for convergence
      if(resnorm < 1.0e-2) break;
      
      // apply inverse to new_sol to update solution
      // No need to do matvec with scaled copy of sol.
      // Reuse previously stored matvec from residual calculation
      gm_->applyInv(new_sol); // can also do gemv with gm_->getInverse()
      
      // compute 'shifted' eigenvalue
      beta = sol.dot(new_sol);
      // scale beta by gamma to account for normalizing sol
      beta *= gamma;
      // update solution data
      sol = new_sol;
      // compute norm and update gamma
      alpha = sol.nrm2();
      gamma = 1./alpha;   
   }
   // compute first extent (eigenvalue)
   double e1 = beta - shft;   
   sol.copyDataToVector(vec1);
  
   // shift matrix by beta and compute second extent
   // store shift
   double shft_e1 = -beta;
   mat.axpy(shft_e1, smat);

   // reset data and begin loop
   sol.assignColumn(&vec2[0], 0);
   new_sol.assignColumn(&vec[0], 0);
   alpha = sol.nrm2();
   gamma = 1./alpha;
   beta = sol.dot(new_sol);      
   
   // loop   
   if(onpe0)cout<<"e2:: ITER 0:: = "<< beta<<endl;   
   int iter2=0;
   for(int i=0; i<maxits; i++)
   {
      iter2++;
      
      // First compute residual for convergence check
      res.gemv('N', one, mat, sol, zero);
      // store matvec result appropriately scaled for later reuse
      new_sol.clear();
      new_sol.axpy(gamma,res);       
      // Compute residual: res = beta*S*x - mat*x
      res.gemm('N', 'N', beta, gm_->getMatrix(), sol, -1.);      
      // compute residual norm
      double resnorm = res.nrm2();
      // check for convergence
      if(resnorm < 1.0e-2) break;
      
      // apply inverse to new_sol to update solution
      // No need to do matvec with scaled copy of sol.
      // Reuse previously stored matvec from residual calculation
      gm_->applyInv(new_sol); // can also do gemv with gm_->getInverse()
      
      // compute 'shifted' eigenvalue
      beta = sol.dot(new_sol);
      // scale beta by gamma to account for not normalizing sol
      beta *= gamma;
      // update solution data
      sol = new_sol;
      // compute norm and update gamma
      alpha = sol.nrm2();
      gamma = 1./alpha; 
   }
   // compute second extent
   double e2 = beta - shft_e1 - shft;
   sol.copyDataToVector(vec2);

   // save results
   double tmp = e1;
   e1 = min(tmp, e2);
   e2 = max(tmp, e2);
   double padding = pad*(e2-e1);

if(onpe0)cout<<"Power method Eigen intervals********************  = ( "<<e1<<", "<<e2<<")"<<"iter1 = "<<iter1<<", iter2 = "<<iter2<<endl;   

   e1 -= padding;
   e2 += padding; 
   interval.push_back(e1);
   interval.push_back(e2);   
   
   // update shft
   shft = max(fabs(e1), fabs(e2));

   eig_interval_tm_.stop();
}
