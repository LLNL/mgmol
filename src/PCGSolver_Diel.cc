// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "PCGSolver_Diel.h"

template <class T, typename ScalarType>
void PCGSolver_Diel<T, ScalarType>::clear()
{
    for (short i = 0; i < (short)pc_oper_.size(); i++)
    {
        delete pc_oper_[i];
    }
    for (short i = 0; i < (short)gf_work_.size(); i++)
    {
        assert(gf_work_[i] != nullptr);
        delete gf_work_[i];
    }
    for (short i = 0; i < (short)gf_rcoarse_.size(); i++)
    {
        assert(gf_rcoarse_[i] != nullptr);
        delete gf_rcoarse_[i];
    }
    for (short i = 0; i < (short)gf_newv_.size(); i++)
    {
        assert(gf_newv_[i] != nullptr);
        delete gf_newv_[i];
    }
    // delete grids after pb::GridFunc<ScalarType> objects since those
    // have data members references to grids
    for (short i = 0; i < (short)grid_.size(); i++)
    {
        delete grid_[i];
    }
    pc_oper_.clear();
    grid_.clear();
    gf_work_.clear();
    gf_rcoarse_.clear();
    gf_newv_.clear();
}

template <class T, typename ScalarType>
void PCGSolver_Diel<T, ScalarType>::setupPrecon()
{
    // fine level
    pb::Grid* mygrid = new pb::Grid(oper_.grid());
    grid_.push_back(mygrid);
    const short nghosts = mygrid->ghost_pt();

    T* myoper = new T(oper_);
    pc_oper_.push_back(myoper);

    pb::GridFunc<ScalarType>* gf_work
        = new pb::GridFunc<ScalarType>(*grid_[0], bc_[0], bc_[1], bc_[2]);
    gf_work_.push_back(gf_work);

    // coarse levels
    nlevels_ = max_nlevels_;
    for (short ln = 1; ln <= max_nlevels_; ln++)
    {

        const bool flag_coarsen
            = ((!(mygrid->dim(0)
                   & 1)) // cannot coarsen if mesh not divisible by 2
                && (!(mygrid->dim(1) & 1)) && (!(mygrid->dim(2) & 1))
                && (static_cast<int>(mygrid->dim(0)) >= 2 * nghosts)
                && (static_cast<int>(mygrid->dim(1)) >= 2 * nghosts)
                && (static_cast<int>(mygrid->dim(2)) >= 2 * nghosts));

        if (!flag_coarsen)
        {
            nlevels_ = ln - 1;
            break;
        }

        pb::Grid* coarse_grid = new pb::Grid(mygrid->coarse_grid());
        grid_.push_back(coarse_grid);

        T* myoper = new T(pc_oper_[ln - 1]->coarseOp(*mygrid));
        pc_oper_.push_back(myoper);

        gf_work = new pb::GridFunc<ScalarType>(
            *coarse_grid, bc_[0], bc_[1], bc_[2]);
        gf_work_.push_back(gf_work);

        pb::GridFunc<ScalarType>* gf_rcoarse = new pb::GridFunc<ScalarType>(
            *coarse_grid, bc_[0], bc_[1], bc_[2]);
        gf_rcoarse_.push_back(gf_rcoarse);
        pb::GridFunc<ScalarType>* gf_newv = new pb::GridFunc<ScalarType>(
            *coarse_grid, bc_[0], bc_[1], bc_[2]);
        gf_newv_.push_back(gf_newv);

        mygrid = coarse_grid;
    }
}

template <class T, typename ScalarType>
// MG V-cycle with no mask
void PCGSolver_Diel<T, ScalarType>::preconSolve(pb::GridFunc<ScalarType>& gf_v,
    const pb::GridFunc<ScalarType>& gf_f, const short level)
{
    //(*MPIdata::sout)<<"Preconditioning::mg() at level "<<level<<endl;
    short ncycl = nu1_;
    if (level == nlevels_)
    {
        ncycl = 4 > (nu1_ + nu2_) ? 4 : (nu1_ + nu2_);
    }

    //    pb::Lap* myoper=pc_oper_[level];
    T* myoper = pc_oper_[level];

    // SMOOTHING
    for (short it = 0; it < ncycl; it++)
    {
        myoper->jacobi(gf_v, gf_f, *gf_work_[level]);
    }

    if (level == nlevels_) return;

    // COARSE GRID CORRECTION

    // restrictions
    pb::GridFunc<ScalarType>* rcoarse = gf_rcoarse_[level];
    gf_work_[level]->restrict3D(*rcoarse);

    // storage functions for coarse grid
    pb::GridFunc<ScalarType>* newv = gf_newv_[level];

    // call mgrid solver on a coarser level
    newv->resetData();
    preconSolve(*newv, *rcoarse, level + 1);

    gf_work_[level]->extend3D(*newv);

    gf_v -= (*gf_work_[level]);

    // post-smoothing
    for (short it = 0; it < nu2_; it++)
    {
        myoper->jacobi(gf_v, gf_f, *gf_work_[level]);
    }

    if (bc_[0] != 1 || bc_[2] != 1 || bc_[2] != 1) gf_v.trade_boundaries();
}

template <class T, typename ScalarType>
// Left Preconditioned CG
bool PCGSolver_Diel<T, ScalarType>::solve(
    pb::GridFunc<ScalarType>& gf_phi, pb::GridFunc<ScalarType>& gf_rhs)
{
    if (!oper_.initialized())
    {
        std::cout
            << "Error in PCGSolver_Diel<T>::solve: operator not initialized"
            << std::endl;
        return 0.;
    }

    bool converged           = false;
    const pb::Grid& finegrid = gf_phi.grid();

    // initial data and residual - We assume a nonzero initial guess
    pb::GridFunc<ScalarType> lhs(finegrid, bc_[0], bc_[1], bc_[2]);
    pb::GridFunc<ScalarType> res(finegrid, bc_[0], bc_[1], bc_[2]);
    // scale initial guess with epsilon
    oper_.inv_transform(gf_phi);
    // compute initial residual
    oper_.apply(gf_phi, lhs);
    pb::GridFunc<ScalarType> rhs(gf_rhs);
    oper_.transform(rhs);
    // Hartree units
    rhs *= (4. * M_PI);
    res.diff(rhs, lhs);
    double init_rnorm = res.norm2();
    double rnorm      = init_rnorm;

    // preconditioned residual
    pb::GridFunc<ScalarType> z(finegrid, bc_[0], bc_[1], bc_[2]);
    // preconditioning step
    z = 0.;
    preconSolve(z, res, 0);
    // conjugate vectors
    pb::GridFunc<ScalarType> p(z);
    pb::GridFunc<ScalarType> ap(p.grid(), bc_[0], bc_[1], bc_[2]);

    double rtz = res.gdot(z);

    // main loop
    for (int k = 0; k < maxiters_; k++)
    {
        // matvec: ap = A*p
        oper_.apply(p, ap);
        double ptap = p.gdot(ap);
        double alp  = rtz / ptap;
        //      if(onpe0)printf("rtz = %f, ptap = %f, alp = %f \n",rtz, ptap,
        //      alp);
        // update solution
        gf_phi += p * alp;
        res -= ap * alp;

        // check for convergence
        rnorm = res.norm2();
        if (rnorm <= tol_ * init_rnorm)
        {
            converged = true;
            break;
        }
        z = 0.;
        preconSolve(z, res, 0);
        double rtz_new = res.gdot(z);
        double bet     = rtz_new / rtz;
        p              = z + p * bet;
        rtz            = rtz_new;
        //      if(onpe0)printf("rnorm = %f, bet = %f \n", rnorm, bet);
    }
    oper_.transform(gf_phi);
    final_residual_     = rnorm;
    residual_reduction_ = rnorm / init_rnorm;

    return converged;
}

template <class T, typename ScalarType>
// Left Preconditioned CG
bool PCGSolver_Diel<T, ScalarType>::solve(pb::GridFunc<ScalarType>& gf_phi,
    pb::GridFunc<ScalarType>& gf_rhs, pb::GridFunc<ScalarType>& gf_rhod,
    pb::GridFunc<ScalarType>& gf_vks)
{
    // initialize the linear system operator and the preconditioner
    oper_.init(gf_rhod);

    // setup precon data
    // We assume operator has changed hence new precon is constructed
    setupPrecon();

    bool converged = solve(gf_phi, gf_rhs);

    // compute additional KS potential due to epsilon(rho)
    oper_.get_vepsilon(gf_phi, gf_rhod, gf_vks);

    // release precon data
    clear();
    return converged;
}

template class PCGSolver_Diel<pb::PBh4MP<double>, double>;
template class PCGSolver_Diel<pb::PBh4MP<float>, float>;
template class PCGSolver_Diel<pb::PBh4M<double>, double>;
template class PCGSolver_Diel<pb::PBh4M<float>, float>;
template class PCGSolver_Diel<pb::PBh4<double>, double>;
template class PCGSolver_Diel<pb::PBh4<float>, float>;
template class PCGSolver_Diel<pb::PBh2<double>, double>;
template class PCGSolver_Diel<pb::PBh2<float>, float>;
template class PCGSolver_Diel<pb::PBh6<double>, double>;
template class PCGSolver_Diel<pb::PBh6<float>, float>;
template class PCGSolver_Diel<pb::PBh8<double>, double>;
template class PCGSolver_Diel<pb::PBh8<float>, float>;
