// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Potentials.h"
#include "Control.h"
#include "Delh4.h"
#include "Grid.h"
#include "GridFunc.h"
#include "Ions.h"
#include "MGmol_MPI.h"
#include "MGmol_blas1.h"
#include "MPIdata.h"
#include "Mesh.h"
#include "Species.h"
#include "tools.h"

#include "SuperSampling.h"
#include "TriCubic.h"
#include "mputils.h"

#include <fstream>
using namespace std;

// unit conversion factor Ha -> Ry
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

    scf_dvrho_ = 1000.;
    scf_dv_    = 1000.;

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
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        size_, one, &v_ext_[0], &vtot_[0]);
    // factor ha2ry to get total potential in [Ry] for calculations
    LinearAlgebraUtils<MemorySpace::Host>::MPscal(size_, ha2ry, &vtot_[0]);
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

    for (unsigned short is = 0; is < rho.size(); is++)
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
    if (diel_)
        LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
            size_, two, &vepsilon_[0], &vtot_[0]);

    // evaluate correction of vtot
    double minus = -1.;
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        size_, minus, &vtot_[0], &dv_[0]);

    evalNormDeltaVtotRho(rho);

    double dvdot
        = LinearAlgebraUtils<MemorySpace::Host>::MPdot(size_, &dv_[0], &dv_[0]);

    double sum = 0.;
    int rc
        = MPI_Allreduce(&dvdot, &sum, 1, MPI_DOUBLE, MPI_SUM, myPEenv.comm());
    if (rc != MPI_SUCCESS)
    {
        cout << "MPI_Allreduce double sum failed!!!" << endl;
        MPI_Abort(myPEenv.comm(), 2);
    }
    dvdot = sum;

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
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        size_, potmix, &dv_[0], &vtot_[0]);
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
    if (diel_)
        LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
            size_, two, &vepsilon_[0], &dv_[0]);

    // evaluate correction of vtot
    double minus = -1.;
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        size_, minus, &vtot_old_[0], &dv_[0]);

    evalNormDeltaVtotRho(rho);

    double dvdot
        = LinearAlgebraUtils<MemorySpace::Host>::MPdot(size_, &dv_[0], &dv_[0]);

    double sum = 0.;
    int rc
        = MPI_Allreduce(&dvdot, &sum, 1, MPI_DOUBLE, MPI_SUM, myPEenv.comm());
    if (rc != MPI_SUCCESS)
    {
        cout << "MPI_Allreduce double sum failed!!!" << endl;
        MPI_Abort(myPEenv.comm(), 2);
    }
    dvdot = sum;

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
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        size_, minustwo, &v_nuc_[0], &vrho[0]);
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        size_, minustwo, &v_ext_[0], &vrho[0]);
}

#ifdef HAVE_TRICUBIC
// type:
// 2->text
// 3->binary
void Potentials::readExternalPot(const string filename, const char type)
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
                mmpi.abort();
            }
        }
    }
    if (!from)
    {
        (*MPIdata::serr) << " Cannot open file " << filename << endl;
        mmpi.abort();
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
                mmpi.abort();
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
                mmpi.abort();
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
            mmpi.abort();
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
                    // v_ext_["<<index<<"]="<<v_ext_[index]<<endl;
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
        if (pot_types_[isp] == 'n' || pot_types_[isp] == 's'
            || pot_types_[isp] == 'f')
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

    itindex_vxc_ = iterativeIndex;
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
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(
        size_, alpha, &v_comp_[0], &vh_rho_[0]);
}

void Potentials::axpVcomp(POTDTYPE* v, const double alpha)
{
    LinearAlgebraUtils<MemorySpace::Host>::MPaxpy(size_, alpha, &v_comp_[0], v);
}

void Potentials::initializeSupersampledRadialDataOnMesh(
    const Vector3D& position, const Species& sp)
{
    Control& ct = *(Control::instance());

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const int dim0 = mygrid.dim(0);
    const int dim1 = mygrid.dim(1);
    const int dim2 = mygrid.dim(2);

    const double start0 = mygrid.start(0);
    const double start1 = mygrid.start(1);
    const double start2 = mygrid.start(2);

    const double h0 = mygrid.hgrid(0);
    const double h1 = mygrid.hgrid(1);
    const double h2 = mygrid.hgrid(2);

    Vector3D point(0., 0., 0.);

    const double lrad = sp.lradius();

    const RadialInter& lpot(sp.local_pot());
    auto lambda_lpot = [&lpot](double radius) { return lpot.cubint(radius); };
    const Vector3D lattice(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));

    // Debug stuff
    // ofstream myfile;
    // myfile.open("SuperSamplingDebug7-25.txt");
    // Debug stuff

    // Construct subdomain containing molecule
    std::array<double, 3> atomicCenter
        = { position[0], position[1], position[2] };
    std::array<double, 3> botMeshCorner          = { 0, 0, 0 };
    std::array<double, 3> topMeshCorner          = { 0, 0, 0 };
    std::array<double, 3> subDomainBotMeshCorner = { start0, start1, start2 };
    std::array<double, 3> subDomainTopMeshCorner
        = { start0 + dim0 * h0, start1 + dim1 * h1, start2 + dim2 * h2 };
    std::array<int, 3> SSLRad = { 0, 0, 0 };
    SSLRad[0]
        = std::ceil(lrad / h0)
          + 1; // +1 just to be safe and make sure subdomain gets everything
    SSLRad[1]        = std::ceil(lrad / h1) + 1;
    SSLRad[2]        = std::ceil(lrad / h2) + 1;
    botMeshCorner[0] = std::max(std::round((atomicCenter[0] - start0) / h0) * h0
                                    + start0 - SSLRad[0] * h0,
        subDomainBotMeshCorner[0]);
    botMeshCorner[1] = std::max(std::round((atomicCenter[1] - start1) / h1) * h1
                                    + start1 - SSLRad[1] * h1,
        subDomainBotMeshCorner[1]);
    botMeshCorner[2] = std::max(std::round((atomicCenter[2] - start2) / h2) * h2
                                    + start2 - SSLRad[2] * h2,
        subDomainBotMeshCorner[2]);
    topMeshCorner[0] = std::min(
        botMeshCorner[0] + 2 * h0 * SSLRad[0], subDomainTopMeshCorner[0]);
    topMeshCorner[1] = std::min(
        botMeshCorner[1] + 2 * h1 * SSLRad[1], subDomainTopMeshCorner[1]);
    topMeshCorner[2] = std::min(
        botMeshCorner[2] + 2 * h2 * SSLRad[2], subDomainTopMeshCorner[2]);
    const bool harmonics = false;

    SuperSampling<0> current(
        atomicCenter, botMeshCorner, topMeshCorner, harmonics, lambda_lpot);
    int xlimits = std::round((topMeshCorner[0] - botMeshCorner[0]) / h0);
    int ylimits = std::round((topMeshCorner[1] - botMeshCorner[1]) / h1);
    int zlimits = std::round((topMeshCorner[2] - botMeshCorner[2]) / h2);
    int xoffset = std::round((botMeshCorner[0] - start0) / h0);
    int yoffset = std::round((botMeshCorner[1] - start1) / h1);
    int zoffset = std::round((botMeshCorner[2] - start2) / h2);
    int offset  = 0;

    for (int ix = xoffset; ix <= xoffset + xlimits; ix++)
    {
        int istart = ix * dim1 * dim2;
        for (int iy = yoffset; iy <= yoffset + ylimits; iy++)
        {
            int jstart = istart + iy * dim2;
            for (int iz = zoffset; iz <= zoffset + zlimits; iz++)
            {
                const double r
                    = position.minimage(point, lattice, ct.bcPoisson);
                v_nuc_[jstart + iz] += current.values_[0][offset];
                rho_comp_[jstart + iz] += sp.getRhoComp(r);
                v_comp_[jstart + iz] += sp.getVcomp(r);
                offset++;
            }
        }
    }

    /*
    offset = 0;
    for (int ix = 0; ix < dim0; ix++)
    {
        point[0]   = start0 + ix * h0;

        for (int iy = 0; iy < dim1; iy++)
        {
            point[1]   = start1 + iy * h1;

            for (int iz = 0; iz < dim2; iz++)
            {
            point[2] = start2 + iz * h2;

                const double r
                    = position.minimage(point, lattice, ct.bcPoisson);

                if (r < lrad)
                {
                    rho_comp_[offset] += sp.getRhoComp(r);
                    v_comp_[offset]   += sp.getVcomp(r);
                    //myfile << v_nuc_[offset] << ';' << lpot.cubint(r) << ',';
                    //v_nuc_[offset]    += lpot.cubint(r);
                }

                offset++;
            }
        }
    }
    //myfile.close();
    */
}

void Potentials::initializeRadialDataOnMesh(
    const Vector3D& position, const Species& sp)
{
    Control& ct = *(Control::instance());

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    const int dim0 = mygrid.dim(0);
    const int dim1 = mygrid.dim(1);
    const int dim2 = mygrid.dim(2);

    const double start0 = mygrid.start(0);
    const double start1 = mygrid.start(1);
    const double start2 = mygrid.start(2);

    const double h0 = mygrid.hgrid(0);
    const double h1 = mygrid.hgrid(1);
    const double h2 = mygrid.hgrid(2);

    Vector3D point(0., 0., 0.);

    const double lrad = sp.lradius();

    const RadialInter& lpot(sp.local_pot());
    const Vector3D lattice(mygrid.ll(0), mygrid.ll(1), mygrid.ll(2));

    int offset = 0;
    for (int ix = 0; ix < dim0; ix++)
    {
        point[0] = start0 + ix * h0;

        for (int iy = 0; iy < dim1; iy++)
        {
            point[1] = start1 + iy * h1;

            for (int iz = 0; iz < dim2; iz++)
            {
                point[2] = start2 + iz * h2;

                const double r
                    = position.minimage(point, lattice, ct.bcPoisson);

                if (r < lrad)
                {
                    rho_comp_[offset] += sp.getRhoComp(r);
                    v_comp_[offset] += sp.getVcomp(r);
                    v_nuc_[offset] += lpot.cubint(r);
                }

                offset++;
            }
        }
    }
}

// Initialization of the compensating potential
void Potentials::initialize(Ions& ions)
{
    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();
    const int numpt        = mygrid.size();

    memset(&v_comp_[0], 0, numpt * sizeof(POTDTYPE));
    memset(&rho_comp_[0], 0, numpt * sizeof(RHODTYPE));
    memset(&v_nuc_[0], 0, numpt * sizeof(RHODTYPE));

    char flag_filter = pot_type(0);

    // Loop over ions
    for (auto& ion : ions.overlappingVL_ions())
    {
        const Species& sp(ion->getSpecies());

        Vector3D position(ion->position(0), ion->position(1), ion->position(2));

        if (flag_filter == 's')
        {
            const int sampleRate  = 3;
            const int numExtraPts = 40;
            const std::array<double, 3> coarGridSpace
                = { mygrid.hgrid(0), mygrid.hgrid(1), mygrid.hgrid(2) };
            SuperSampling<0>::setup(sampleRate, numExtraPts, coarGridSpace);

            initializeSupersampledRadialDataOnMesh(position, sp);
        }
        else
        {
            initializeRadialDataOnMesh(position, sp);
        }
    }
}

void Potentials::rescaleRhoComp()
{
    Control& ct     = *(Control::instance());
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());

    Mesh* mymesh           = Mesh::instance();
    const pb::Grid& mygrid = mymesh->grid();

    // Check compensating charges
    double comp_rho = getCharge(&rho_comp_[0]);

    if (onpe0 && ct.verbose > 1)
    {
        cout << " Rescaling rhoc" << endl;
    }
    if (ionic_charge_ > 0.)
    {
        const int numpt = mygrid.size();
        double t        = ionic_charge_ / comp_rho;
        LinearAlgebraUtils<MemorySpace::Host>::MPscal(numpt, t, &rho_comp_[0]);

        // Check new compensating charges
        comp_rho = getCharge(&rho_comp_[0]);
    }
    if (onpe0 && ct.verbose > 1)
        cout << " Rescaled compensating charges: " << setprecision(8) << fixed
             << comp_rho << endl;
    if (comp_rho < 0.) mmpi.abort();
}

void Potentials::addBackgroundToRhoComp()
{
    if (fabs(background_charge_) > 0.)
    {
        Control& ct            = *(Control::instance());
        Mesh* mymesh           = Mesh::instance();
        const pb::Grid& mygrid = mymesh->grid();
        const int numpt        = mygrid.size();

        double background
            = background_charge_ / (mygrid.gsize() * mygrid.vel());
        if (ct.bcPoisson[0] == 1 && ct.bcPoisson[1] == 1
            && ct.bcPoisson[2] == 1)
        {
            if (onpe0)
            {
                cout << setprecision(12) << scientific
                     << "Add background charge " << background << " to rhoc "
                     << endl;
            }
            for (int i = 0; i < numpt; i++)
                rho_comp_[i] += background;

            // Check new compensating charges
            getCharge(&rho_comp_[0]);
        }
    }
}

void Potentials::initBackground(Ions& ions)
{
    Control& ct = *(Control::instance());

    // Count up the total ionic charge
    ionic_charge_ = ions.computeIonicCharge();

    // calculation the compensating background charge
    //   for charged supercell calculations
    background_charge_ = 0.;
    charge_in_cell_    = ionic_charge_ - ct.getNel();
    if (ct.bcPoisson[0] != 2 && ct.bcPoisson[1] != 2 && ct.bcPoisson[2] != 2)
    {
        background_charge_ = (-1.) * charge_in_cell_;
    }
    if (onpe0 && ct.verbose > 0)
    {
        cout << "N electrons=      " << ct.getNel() << endl;
        cout << "ionic charge=     " << ionic_charge_ << endl;
        cout << "background charge=" << background_charge_ << endl;
    }

    if (fabs(background_charge_) < 1.e-10) background_charge_ = 0.;
}

template void Potentials::setVxc<double>(
    const double* const vxc, const int iterativeIndex);
template void Potentials::setVxc<float>(
    const float* const vxc, const int iterativeIndex);
