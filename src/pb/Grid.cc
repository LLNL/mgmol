// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "Grid.h"
#include <cmath>
#include <vector>

namespace pb
{

// constructor
Grid::Grid(const double origin[3], const double lattice[3],
    const unsigned ngpts[3], const PEenv& mype_env, const short nghosts,
    const short level)
    : mype_env_(mype_env)
{
    assert(lattice[0] > 1.e-8);
    assert(lattice[1] > 1.e-8);
    assert(lattice[2] > 1.e-8);

    assert(ngpts[0] > 0);
    assert(ngpts[1] > 0);
    assert(ngpts[2] > 0);

    assert(ngpts[0] < 10000);
    assert(ngpts[1] < 10000);
    assert(ngpts[2] < 10000);

    origin_[0] = origin[0];
    origin_[1] = origin[1];
    origin_[2] = origin[2];

    ll_[0] = lattice[0];
    ll_[1] = lattice[1];
    ll_[2] = lattice[2];

    level_ = level;

    assert(mype_env.n_mpi_task(0) > 0);
    assert(mype_env.n_mpi_task(1) > 0);
    assert(mype_env.n_mpi_task(2) > 0);

    dim_[0] = ngpts[0] / mype_env.n_mpi_task(0);
    dim_[1] = ngpts[1] / mype_env.n_mpi_task(1);
    dim_[2] = ngpts[2] / mype_env.n_mpi_task(2);
    assert(dim_[0] * mype_env.n_mpi_task(0) == ngpts[0]);
    assert(dim_[1] * mype_env.n_mpi_task(1) == ngpts[1]);
    assert(dim_[2] * mype_env.n_mpi_task(2) == ngpts[2]);

    gdim_[0] = ngpts[0];
    gdim_[1] = ngpts[1];
    gdim_[2] = ngpts[2];

    ghost_pt_ = nghosts;

    // cout<< "Construct grid of
    // dim"<<dim(0)<<","<<dim(1)<<","<<dim(2)<<endl;

    size_  = dim_[0] * dim_[1] * dim_[2];
    sizeg_ = (dim_[0] + 2 * nghosts) * (dim_[1] + 2 * nghosts)
             * (dim_[2] + 2 * nghosts);

    gsize_ = gdim_[0] * gdim_[1] * gdim_[2];

    for (short i = 0; i < 3; i++)
    {
        hgrid_[i] = ll_[i] / (double)gdim_[i];
        assert(hgrid_[i] > 1.e-8);
    }
    // cout<<"h="<<hgrid(0)<<","<<hgrid(1)<<","<<hgrid(2)<<endl;

    vel_ = hgrid_[0] * hgrid_[1] * hgrid_[2];

    inc_[0] = (dim_[1] + 2 * ghost_pt_) * (dim_[2] + 2 * ghost_pt_);
    inc_[1] = (dim_[2] + 2 * ghost_pt_);
    inc_[2] = 1;

    for (short i = 0; i < 3; i++)
    {
        start_[i]  = mype_env.my_mpi(i) * dim_[i] * hgrid_[i] + origin_[i];
        istart_[i] = mype_env.my_mpi(i) * dim_[i];
    }

    for (short i = 0; i < 3; i++)
    {
        assert(dim_[i] > 0);
        assert(dim_[i] < 10000);
    }
}

// copy constructor
Grid::Grid(const Grid& my_grid, const short nghosts)
    : mype_env_(my_grid.mype_env_)
{
    dim_[0] = my_grid.dim_[0];
    dim_[1] = my_grid.dim_[1];
    dim_[2] = my_grid.dim_[2];

    if (nghosts == -1)
    {
        assert(static_cast<unsigned int>(my_grid.ghost_pt_) <= dim_[0]);
        assert(static_cast<unsigned int>(my_grid.ghost_pt_) <= dim_[1]);
        assert(static_cast<unsigned int>(my_grid.ghost_pt_) <= dim_[2]);

        ghost_pt_ = my_grid.ghost_pt_; // default
        inc_[0]   = my_grid.inc_[0];
        inc_[1]   = my_grid.inc_[1];
        inc_[2]   = my_grid.inc_[2];
        sizeg_    = my_grid.sizeg_;
    }
    else
    {
        assert(static_cast<unsigned int>(nghosts) <= dim_[0]);
        assert(static_cast<unsigned int>(nghosts) <= dim_[1]);
        assert(static_cast<unsigned int>(nghosts) <= dim_[2]);

        ghost_pt_ = nghosts;
        inc_[0]   = (dim_[1] + 2 * ghost_pt_) * (dim_[2] + 2 * ghost_pt_);
        inc_[1]   = (dim_[2] + 2 * ghost_pt_);
        inc_[2]   = 1;
        sizeg_    = (dim_[0] + 2 * ghost_pt_) * (dim_[1] + 2 * ghost_pt_)
                 * (dim_[2] + 2 * ghost_pt_);
    }

    // cout<<"Copy const. for grid\n";
    gdim_[0]   = my_grid.gdim(0);
    gdim_[1]   = my_grid.gdim(1);
    gdim_[2]   = my_grid.gdim(2);
    origin_[0] = my_grid.origin_[0];
    origin_[1] = my_grid.origin_[1];
    origin_[2] = my_grid.origin_[2];
    ll_[0]     = my_grid.ll_[0];
    ll_[1]     = my_grid.ll_[1];
    ll_[2]     = my_grid.ll_[2];
    hgrid_[0]  = my_grid.hgrid(0);
    hgrid_[1]  = my_grid.hgrid(1);
    hgrid_[2]  = my_grid.hgrid(2);
    start_[0]  = my_grid.start(0);
    start_[1]  = my_grid.start(1);
    start_[2]  = my_grid.start(2);
    istart_[0] = my_grid.istart(0);
    istart_[1] = my_grid.istart(1);
    istart_[2] = my_grid.istart(2);

    size_  = my_grid.size_;
    gsize_ = my_grid.gsize_;
    vel_   = my_grid.vel_;
    level_ = my_grid.level_;

    for (short i = 0; i < 3; i++)
    {
        assert(dim_[i] > 0);
        assert(dim_[i] < 10000);
    }
}

Grid& Grid::operator=(const Grid& my_grid)
{
    if (this != &my_grid)
    {
        // cout<<" operator= for grid\n";
        assert(my_grid.dim(0) > 0);
        assert(my_grid.dim(1) > 0);
        assert(my_grid.dim(2) > 0);

        sizeg_ = my_grid.sizeg();
        size_  = my_grid.size();
        gsize_ = my_grid.gsize();

        dim_[0]    = my_grid.dim(0);
        dim_[1]    = my_grid.dim(1);
        dim_[2]    = my_grid.dim(2);
        gdim_[0]   = my_grid.gdim(0);
        gdim_[1]   = my_grid.gdim(1);
        gdim_[2]   = my_grid.gdim(2);
        origin_[0] = my_grid.origin(0);
        origin_[1] = my_grid.origin(1);
        origin_[2] = my_grid.origin(2);
        ll_[0]     = my_grid.ll(0);
        ll_[1]     = my_grid.ll(1);
        ll_[2]     = my_grid.ll(2);
        hgrid_[0]  = my_grid.hgrid(0);
        hgrid_[1]  = my_grid.hgrid(1);
        hgrid_[2]  = my_grid.hgrid(2);
        inc_[0]    = my_grid.inc(0);
        inc_[1]    = my_grid.inc(1);
        inc_[2]    = my_grid.inc(2);
        start_[0]  = my_grid.start(0);
        start_[1]  = my_grid.start(1);
        start_[2]  = my_grid.start(2);
        istart_[0] = my_grid.istart(0);
        istart_[1] = my_grid.istart(1);
        istart_[2] = my_grid.istart(2);

        ghost_pt_ = my_grid.ghost_pt_;
        vel_      = my_grid.vel_;
        level_    = my_grid.level_;
    }
    for (short i = 0; i < 3; i++)
    {
        assert(dim_[i] > 0);
        assert(dim_[i] < 10000);
    }

    return *this;
}

const Grid Grid::coarse_grid() const
{
    // cout<<" Build coarse Grid\n";
    assert((gdim(0) % 2) == 0);
    assert((gdim(1) % 2) == 0);
    assert((gdim(2) % 2) == 0);

    unsigned dim[3] = { gdim_[0] >> 1, gdim_[1] >> 1, gdim_[2] >> 1 };
    Grid coarse_G(origin_, ll_, dim, mype_env_, ghost_pt_, level_ - 1);
    // cout<<"gsize="<<gsize()<<endl;
    // cout<<"coarse_G.gsize="<<coarse_G.gsize()<<endl;
    // cout<<"size="<<size()<<endl;
    // cout<<"coarse_G.size="<<coarse_G.size()<<endl;
    assert(gsize() == 8 * coarse_G.gsize());
    assert(size() == 8 * coarse_G.size());

    return coarse_G;
}

const Grid Grid::replicated_grid(const PEenv& replicated_peenv) const
{
    // cout<<" Build replicated Grid\n";

    Grid replicated_grid(
        origin_, ll_, gdim_, replicated_peenv, ghost_pt_, level_);

    return replicated_grid;
}

template <typename T>
void Grid::getSinCosFunctions(std::vector<T>& sinx, std::vector<T>& siny,
    std::vector<T>& sinz, std::vector<T>& cosx, std::vector<T>& cosy,
    std::vector<T>& cosz) const
{
    sinx.resize(dim_[0]);
    siny.resize(dim_[1]);
    sinz.resize(dim_[2]);
    cosx.resize(dim_[0]);
    cosy.resize(dim_[1]);
    cosz.resize(dim_[2]);

    const T inv_2pi = 0.5 * M_1_PI;
    const T alphax  = ll_[0] * inv_2pi;
    const T alphay  = ll_[1] * inv_2pi;
    const T alphaz  = ll_[2] * inv_2pi;

    // use an hh so that the arguments of trigonometric functions go from
    // 0 to 2*pi from origin of mesh to end of mesh
    const T hhx = 2. * M_PI / (T)gdim_[0];
    const T hhy = 2. * M_PI / (T)gdim_[1];
    const T hhz = 2. * M_PI / (T)gdim_[2];

    const int xoff = istart_[0];
    const int yoff = istart_[1];
    const int zoff = istart_[2];

    for (unsigned int i = 0; i < dim_[0]; i++)
        sinx[i] = std::sin(T(xoff + i) * hhx) * alphax;
    for (unsigned int i = 0; i < dim_[0]; i++)
        cosx[i] = std::cos(T(xoff + i) * hhx) * alphax;
    for (unsigned int i = 0; i < dim_[1]; i++)
        siny[i] = std::sin(T(yoff + i) * hhy) * alphay;
    for (unsigned int i = 0; i < dim_[1]; i++)
        cosy[i] = std::cos(T(yoff + i) * hhy) * alphay;
    for (unsigned int i = 0; i < dim_[2]; i++)
        sinz[i] = std::sin(T(zoff + i) * hhz) * alphaz;
    for (unsigned int i = 0; i < dim_[2]; i++)
        cosz[i] = std::cos(T(zoff + i) * hhz) * alphaz;
}

Vector3D Grid::closestGridPt(Vector3D coords) const
{
    for (int i = 0; i < 3; i++)
    {
        double mod_term = floor(coords[i] / hgrid_[i]) * hgrid_[i];
        double mod      = coords[i] - mod_term;

        if (mod < 0.5 * hgrid_[i])
            mod = 0.;
        else
            mod = hgrid_[i];

        coords[i] = mod_term + mod;
    }

    return coords;
}

double Grid::integralOverMesh(const double* const func) const
{
    assert(vel_ > 0.000001);
    assert(vel_ < 1000.);

    double value = 0.;

    for (unsigned int idx = 0; idx < size_; idx++)
    {
        value += func[idx];
    }
    value = vel_ * mype_env_.double_sum_all(value);

    return value;
}

template void Grid::getSinCosFunctions(std::vector<float>& sinx,
    std::vector<float>& siny, std::vector<float>& sinz,
    std::vector<float>& cosx, std::vector<float>& cosy,
    std::vector<float>& cosz) const;
template void Grid::getSinCosFunctions(std::vector<double>& sinx,
    std::vector<double>& siny, std::vector<double>& sinz,
    std::vector<double>& cosx, std::vector<double>& cosy,
    std::vector<double>& cosz) const;

} // namespace pb
