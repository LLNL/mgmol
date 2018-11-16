// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#include "Potentials.h"
#include "Control.h"
#include "Delh4.h"
#include "Grid.h"
#include "GridFunc.h"
#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "Species.h"
#include "tools.h"

#include "TriCubic.h"
#include "mputils.h"

#include <fstream>
using namespace std;

const double ha2ry = 2.;

Potentials::~Potentials()
{
#ifdef HAVE_TRICUBIC
    if (vext_tricubic_ != NULL) delete vext_tricubic_;
#endif
}

Potentials::Potentials(const bool vh_frozen)
{
    //(*MPIdata::sout)<<"Potentials::setup()"<<endl;
    diel_            = false; // default: no dielectric
    itindex_vxc_     = -1;
    itindex_vh_      = -1;
    verbosity_level_ = 0;

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    gdim_[0] = mygrid.gdim(0);
    gdim_[1] = mygrid.gdim(1);
    gdim_[2] = mygrid.gdim(2);

    dim_[0] = mygrid.dim(0);
    dim_[1] = mygrid.dim(1);
    dim_[2] = mygrid.dim(2);

    size_ = dim_[0] * dim_[1] * dim_[2];

    mix_ = 1.;

    scf_dvrho_ = 1000.;
    scf_dv_    = 1000.;

    vh_frozen_ = vh_frozen;

    vtot_.resize(size_);
    vtot_old_.resize(size_);

    vepsilon_.resize(size_);
    vh_rho_.resize(size_);
    vxc_rho_.resize(size_);

    rho_comp_.resize(size_);
    v_comp_.resize(size_);

    v_nuc_.resize(size_);
    v_ext_.resize(size_);

    dv_.resize(size_);

    memset(&vepsilon_[0], 0, size_ * sizeof(POTDTYPE));
    memset(&vh_rho_[0], 0, size_ * sizeof(POTDTYPE));
    memset(&vxc_rho_[0], 0, size_ * sizeof(POTDTYPE));
    memset(&v_ext_[0], 0, size_ * sizeof(POTDTYPE));

#ifdef HAVE_TRICUBIC
    vext_tricubic_ = NULL;
#endif
}

void Potentials::initWithVnuc()
{
    assert(size_ > 0);
    if (verbosity_level_ > 2 && onpe0)
        (*MPIdata::sout) << "Potentials::initWithVnuc()" << endl;
    itindex_vxc_ = 0;
    itindex_vh_  = 0;
    int ione     = 1;
    Tcopy(&size_, &v_nuc_[0], &ione, &vtot_[0], &ione);
    double one = 1.;
    MPaxpy(size_, one, &v_ext_[0], &vtot_[0]);
    // factor 2 to get total potential in [Ry] for calculations
    MPscal(size_, ha2ry, &vtot_[0]);
}

double Potentials::max() const
{
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double vmax              = (*max_element(vtot_.begin(), vtot_.end()));
    vmax                     = myPEenv.double_max_all(vmax);
    return vmax;
}

double Potentials::min() const
{
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();
    double vmin              = -(*min_element(vtot_.begin(), vtot_.end()));
    vmin                     = -myPEenv.double_max_all(vmin);
    return vmin;
}

void Potentials::evalNormDeltaVtotRho(const vector<vector<RHODTYPE>>& rho)
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    scf_dvrho_ = 0.;

    for (short is = 0; is < rho.size(); is++)
        for (int idx = 0; idx < size_; idx++)
        {
            scf_dvrho_ += fabs((double)dv_[idx] * (double)rho[is][idx]);
        }
    scf_dvrho_ *= mygrid.vel();
    scf_dvrho_ *= 0.5;

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.allreduce(&scf_dvrho_, 1, MPI_SUM);
}

double Potentials::update(const vector<vector<RHODTYPE>>& rho)
{
    assert(itindex_vxc_ >= 0);
    assert(itindex_vh_ >= 0);
    assert(itindex_vxc_ == itindex_vh_);

    if (verbosity_level_ > 2 && onpe0)
        (*MPIdata::sout) << "Potentials::update(rho)" << endl;
    int ione                 = 1;
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    // save old potentials
    Tcopy(&size_, &vtot_[0], &ione, &vtot_old_[0], &ione);
    Tcopy(&size_, &vtot_[0], &ione, &dv_[0], &ione);

    // update vtot_ (factor 2. to get Rydbergs)
    for (int idx = 0; idx < size_; idx++)
    {
        vtot_[idx]
            = (POTDTYPE)(ha2ry
                         * ((double)v_nuc_[idx] + (double)v_ext_[idx]
                               + (double)vh_rho_[idx] + (double)vxc_rho_[idx]));
    }
    double two = ha2ry;
    if (diel_) MPaxpy(size_, two, &vepsilon_[0], &vtot_[0]);

    // evaluate correction of vtot
    double minus = -1.;
    MPaxpy(size_, minus, &vtot_[0], &dv_[0]);

    evalNormDeltaVtotRho(rho);

    double dvdot = MPdot(size_, &dv_[0], &dv_[0]);

#ifdef USE_MPI
    double sum = 0.;
    int rc
        = MPI_Allreduce(&dvdot, &sum, 1, MPI_DOUBLE, MPI_SUM, myPEenv.comm());
    if (rc != MPI_SUCCESS)
    {
        cout << "MPI_Allreduce double sum failed!!!" << endl;
        MPI_Abort(myPEenv.comm(), 2);
    }
    dvdot = sum;
#endif

    scf_dv_            = 0.5 * sqrt(dvdot);
    const double gsize = (double)size_ * (double)myPEenv.n_mpi_tasks();
    scf_dv_ /= gsize;

    return scf_dv_;
}

void Potentials::update(const double mix)
{
    assert(itindex_vxc_ == itindex_vh_);

#ifdef DEBUG
    if (onpe0) (*MPIdata::sout) << "Potentials::update(mix)" << endl;
#endif
    //    int ione=1;
    double potmix = mix;
    MPaxpy(size_, potmix, &dv_[0], &vtot_[0]);
}

double Potentials::delta_v(const vector<vector<RHODTYPE>>& rho)
{
    assert(itindex_vxc_ == itindex_vh_);
    assert(size_ > 0);

    if (verbosity_level_ > 2 && onpe0)
        (*MPIdata::sout) << "Potentials::delta_v()" << endl;

    int ione                 = 1;
    Mesh* mymesh             = Mesh::instance();
    const pb::PEenv& myPEenv = mymesh->peenv();

    // save old potentials
    Tcopy(&size_, &vtot_[0], &ione, &vtot_old_[0], &ione);

    // update vtot_ (factor 2. to get Rydbergs)
    for (int idx = 0; idx < size_; idx++)
    {
        dv_[idx]
            = (POTDTYPE)(ha2ry
                         * ((double)v_nuc_[idx] + (double)v_ext_[idx]
                               + (double)vh_rho_[idx] + (double)vxc_rho_[idx]));
    }
    double two = ha2ry;
    if (diel_) MPaxpy(size_, two, &vepsilon_[0], &dv_[0]);

    // evaluate correction of vtot
    double minus = -1.;
    MPaxpy(size_, minus, &vtot_old_[0], &dv_[0]);

    evalNormDeltaVtotRho(rho);

    double dvdot = MPdot(size_, &dv_[0], &dv_[0]);

#ifdef USE_MPI
    double sum = 0.;
    int rc
        = MPI_Allreduce(&dvdot, &sum, 1, MPI_DOUBLE, MPI_SUM, myPEenv.comm());
    if (rc != MPI_SUCCESS)
    {
        cout << "MPI_Allreduce double sum failed!!!" << endl;
        MPI_Abort(myPEenv.comm(), 2);
    }
    dvdot = sum;
#endif

    scf_dv_            = 0.5 * sqrt(dvdot);
    const double gsize = (double)size_ * (double)myPEenv.n_mpi_tasks();
    scf_dv_ /= gsize;

    return scf_dv_;
}

// in Ry
void Potentials::getVofRho(vector<POTDTYPE>& vrho) const
{
    vrho.resize(size_);
    int ione        = 1;
    double minustwo = -2.;

    Tcopy(&size_, &vtot_[0], &ione, &vrho[0], &ione);
    MPaxpy(size_, minustwo, &v_nuc_[0], &vrho[0]);
    MPaxpy(size_, minustwo, &v_ext_[0], &vrho[0]);
}

#ifdef HAVE_TRICUBIC
// type:
// 2->text
// 3->binary
void Potentials::readExternalPot(const string filename, const short type)
{
    assert(type == 2 || type == 3);

    Control& ct            = *(Control::instance());
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    MGmol_MPI& mmpi        = *(MGmol_MPI::instance());

    ifstream* from;

    if (type == 2)
        from = new ifstream(filename.data(), ios::in);
    else if (type == 3)
    {
        from = new ifstream(filename.data(), ios::in | ifstream::binary);
        // get length of file:
        if (mmpi.instancePE0())
        {
            from->seekg(0, ios::end);
            const int length = from->tellg();
            (*MPIdata::sout) << "Length file = " << length << endl;
            from->seekg(0, ios::beg);
            if (length <= 0)
            {
                (*MPIdata::serr)
                    << "ERROR Potential: file length <=0!!!!" << endl;
                ct.global_exit(2);
            }
        }
    }
    if (!from)
    {
        (*MPIdata::serr) << " Cannot open file " << filename << endl;
        ct.global_exit(2);
    }
    if (onpe0)
    {
        (*MPIdata::sout) << "Potentials::read_ExternalPot(), filename="
                         << filename << endl;
        if (type == 2) (*MPIdata::sout) << "text file..." << endl;
        if (type == 3) (*MPIdata::sout) << "binary file..." << endl;
    }

    // read origin and end of cell (to check compatibility)
    float origin[3] = { -1., -1., -1. };
    float end[3]    = { -1., -1., -1. };
    if (type == 3)
    {
        from->read((char*)(&origin[0]), 3 * sizeof(float));
        from->read((char*)end, 3 * sizeof(float));
    }
    else
    {
        read_comments(*from);
        (*from) >> origin[0] >> origin[1] >> origin[2];
        (*from) >> end[0] >> end[1] >> end[2];
    }
    double ll[3]
        = { end[0] - origin[0], end[1] - origin[1], end[2] - origin[2] };

    // check compatibility
    if (onpe0)
        for (short d = 0; d < 3; d++)
        {
            (*MPIdata::sout) << setprecision(8);
            if (fabs(origin[d] - mygrid.origin(d)) > 1.e-3)
            {
                (*MPIdata::serr)
                    << "ERROR Potential: Incompatible cell origin in direction "
                    << d << endl;
                (*MPIdata::serr) << "Potential origin=" << origin[d] << endl;
                (*MPIdata::serr) << "MGmol origin=" << mygrid.origin(d) << endl;
                (*MPIdata::serr)
                    << "Difference=" << fabs(origin[d] - mygrid.origin(d))
                    << endl;
                ct.global_exit(2);
            }
            if (fabs(ll[d] - mygrid.ll(d)) > 1.e-3)
            {
                (*MPIdata::serr) << "ERROR Potential: Incompatible cell "
                                    "dimension in direction "
                                 << d << endl;
                (*MPIdata::serr) << "Potential cell end=" << end[d] << endl;
                (*MPIdata::serr)
                    << "Potential cell dimension=" << ll[d] << endl;
                (*MPIdata::serr)
                    << "MGmol cell dimension=" << mygrid.ll(d) << endl;
                ct.global_exit(2);
            }
        }

    // read mesh size
    int nxyz[3];
    if (type == 3)
    {
        from->read((char*)nxyz, 3 * sizeof(int));
    }
    else
    {
        read_comments(*from);
        (*from) >> nxyz[0] >> nxyz[1] >> nxyz[2];
    }

    // check grid compatibility
    for (short i = 0; i < 3; i++)
        if (nxyz[i] != gdim_[i])
        {
            (*MPIdata::serr) << "Potentials::read_ExternalPot(): dimension "
                             << i << " incompatible with Grid!!!" << endl;
            (*MPIdata::serr)
                << "n=" << nxyz[i] << ", gdim_=" << gdim_[i] << endl;
            ct.global_exit(2);
        }

    const int incx = gdim_[2] * gdim_[1];
    const int incy = gdim_[2];

    // get starting point in file
    int start = mygrid.istart(0) * incx + mygrid.istart(1) * incy;

    const int startz = mygrid.istart(2);
    const int endz   = gdim_[2];

    int index      = 0;
    int file_index = 0;

    if (type == 2)
    {
        double scratch;
        read_comments(*from);
        for (int i = 0; i < dim_[0]; i++)
        {
            // advance (start-file_index) positions
            while (file_index < start)
            {
                (*from) >> scratch;
                file_index++;
            }

            for (int j = 0; j < dim_[1]; j++)
            {
                // advance startz positions
                for (int m = 0; m < startz; m++)
                {
                    (*from) >> scratch;
                }
                for (int k = 0; k < dim_[2]; k++)
                {
                    assert(index < size_);
                    (*from) >> v_ext_[index];
                    //(*MPIdata::sout)<<myPEenv.mytask();
                    //(*MPIdata::sout)<<",
                    //v_ext_["<<index<<"]="<<v_ext_[index]<<endl;
                    index++;
                }

                // advance endz-startz-dim_[2] positions
                for (int m = startz + dim_[2]; m < endz; m++)
                {
                    (*from) >> scratch;
                }
            }
            file_index += gdim_[2] * dim_[1];

            start += incx;
        }
    }

    if (type == 3)
    {
        vector<float> tmp(dim_[2]);
        for (int i = 0; i < dim_[0]; i++)
        {
            // advance (start-file_index) positions
            from->seekg((start - file_index) * sizeof(float), ios::cur);
            file_index = start;

            for (int j = 0; j < dim_[1]; j++)
            {
                // advance startz positions
                from->seekg(startz * sizeof(float), ios::cur);
                from->read((char*)(&tmp[0]), dim_[2] * sizeof(float));
                for (int k = 0; k < dim_[2]; k++)
                {
                    assert(index < size_);
                    v_ext_[index] = tmp[k];
                    index++;
                }

                // advance endz-startz-dim_[2] positions
                from->seekg(
                    (endz - startz - dim_[2]) * sizeof(float), ios::cur);
            }
            file_index += gdim_[2] * dim_[1];

            start += incx;
        }
    }

    delete from;

    setupVextTricubic();
}

void Potentials::setupVextTricubic()
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    Control& ct            = *(Control::instance());

    const short bc[3] = { ct.bcPoisson[0], ct.bcPoisson[1], ct.bcPoisson[2] };

    vext_tricubic_ = new pb::TriCubic<POTDTYPE>(mygrid, bc);

    vext_tricubic_->computeSplineCoeffs(&v_ext_[0]);
}

bool Potentials::withVext() const { return (vext_tricubic_ != NULL); }

void Potentials::getGradVext(const double r[3], double dfdr[3]) const
{
    assert(vext_tricubic_ != NULL);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();
    vext_tricubic_->getGradient(r, dfdr, comm);
}

void Potentials::getValVext(const vector<double>& r, vector<double>& val) const
{
    assert(vext_tricubic_ != NULL);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    MPI_Comm comm   = mmpi.commSameSpin();
    vext_tricubic_->getValues(r, val, comm);
}
#endif

void Potentials::readAll(vector<Species>& sp)
{
    assert(sp.size() <= pot_filenames_.size());

    if (verbosity_level_ > 2 && onpe0)
        (*MPIdata::sout) << "Potentials::readAll() for " << pot_types_.size()
                         << " potentials" << endl;
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    double hmin            = mygrid.hmin();
    if (verbosity_level_ > 2 && onpe0)
        (*MPIdata::sout) << "hmin= " << hmin << endl;

    vector<string>::const_iterator it_filename = pot_filenames_.begin();
    int isp                                    = 0;
    while (it_filename != pot_filenames_.end())
    {
        if (pot_types_[isp] == 0 || pot_types_[isp] == 1)
        {
            assert(isp < (int)sp.size());

            sp[isp].read_1species(*it_filename);

            sp[isp].set_dim_nl(hmin);
            sp[isp].set_dim_l(hmin);
        }
        else
        {
#ifdef HAVE_TRICUBIC
            readExternalPot(*it_filename, pot_types_[isp]);
#else
            (*MPIdata::sout)
                << "ERROR: cannot read external potential "
                << " -> need to compile with Tricubic library" << endl;
#endif
        }
        it_filename++;
        isp++;
    }
}
template <typename T>
void Potentials::setVxc(const T* const vxc, const int iterativeIndex)
{
    assert(iterativeIndex >= 0);
    //    int ione=1;
    itindex_vxc_ = iterativeIndex;
    //    Tcopy(&size_, vxc,  &ione, &vxc_rho_[0], &ione);
    MPcpy(&vxc_rho_[0], vxc, size_);
}
void Potentials::setVh(const POTDTYPE* const vh, const int iterativeIndex)
{
    assert(iterativeIndex >= 0);
    int ione    = 1;
    itindex_vh_ = iterativeIndex;
    Tcopy(&size_, vh, &ione, &vh_rho_[0], &ione);
}

void Potentials::setVh(
    const pb::GridFunc<POTDTYPE>& vh, const int iterativeIndex)
{
    assert(iterativeIndex >= 0);

    itindex_vh_ = iterativeIndex;
    vh.init_vect(&vh_rho_[0], 'd');
}

void Potentials::axpVcompToVh(const double alpha)
{
    //    int ione=1;
    MPaxpy(size_, alpha, &v_comp_[0], &vh_rho_[0]);
    //    daxpy(&size_, &alpha, &v_comp_[0], &ione, &vh_rho_[0], &ione);
}

void Potentials::axpVcomp(POTDTYPE* v, const double alpha)
{
    //    int ione=1;
    MPaxpy(size_, alpha, &v_comp_[0], v);
    //    daxpy(&size_, &alpha, &v_comp_[0], &ione, v, &ione);
}
template void Potentials::setVxc<double>(
    const double* const vxc, const int iterativeIndex);
template void Potentials::setVxc<float>(
    const float* const vxc, const int iterativeIndex);
