// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef PB_GRIDFUNC_H
#define PB_GRIDFUNC_H

#include "Grid.h"
#include "Timer.h"

#include <complex>
#include <fstream>
#include <list>
#include <memory>
#include <string.h>
#include <vector>

namespace pb
{

class PEenv;

template <typename ScalarType>
class GridFunc
{
    bool directionDirichlet_[3];
    bool directionPeriodic_[3];
    bool directionMultipole_[3];
    bool directionNeumann_[3];

    int dim_[3];

    int incx_;
    int incy_;

    int mytask_;
    MPI_Request ns_mpireq_[4];
    MPI_Request ud_mpireq_[4];
    MPI_Request ew_mpireq_[4];

    bool north_;
    bool south_;
    bool up_;
    bool down_;
    bool east_;
    bool west_;

    // private functions
    bool def_const() const;

    void setBoundaryValuesNeumannX(const ScalarType alpha);
    void setBoundaryValuesNeumannY(const ScalarType alpha);
    void setBoundaryValuesNeumannZ(const ScalarType alpha);

    ScalarType one(const double r, const double a)
    {
        (void)r; // unused
        (void)a; // unused

        return 1.;
    }
    ScalarType ran0();
    ScalarType radial_func(
        const double r, const double a, const short ftype = 0);

    // allocate memory for function
    void alloc();

    void setup();

    void resizeBuffers();

    static Timer trade_bc_tm_;
    static Timer extend3D_tm_;
    static Timer restrict3D_tm_;
    static Timer prod_tm_;
    static Timer gather_tm_;
    static Timer scatter_tm_;
    static Timer all_gather_tm_;
    static Timer finishExchangeNorthSouth_tm_;
    static Timer finishExchangeUpDown_tm_;
    static Timer finishExchangeEastWest_tm_;

protected:
    const Grid& grid_;

    // data storage
    std::unique_ptr<ScalarType> memory_;

    // raw pointer to data
    ScalarType* uu_;

    bool updated_boundaries_;
    short bc_[3];
    ScalarType valuesNeumann_[3];

    static GridFunc<ScalarType>* bc_func_;

    static std::vector<ScalarType> buf1_;
    static std::vector<ScalarType> buf2_;
    static std::vector<ScalarType> buf3_;
    static std::vector<ScalarType> buf4_;

    void setValues(const int n, const ScalarType* src);

public:
    // Constructors
    GridFunc(const Grid&, const short, const short, const short);

    // constructor with pointer to data allocation
    GridFunc(const Grid&, const short, const short, const short, ScalarType*,
        const bool updated_boundaries = false);

    // copy constructor
    GridFunc(const GridFunc<double>& A);
    GridFunc(const GridFunc<float>& A);

    // copy constructor on different grid
    GridFunc(const GridFunc<ScalarType>&, const Grid&);

    void setValues(const GridFunc<ScalarType>& src);
    void setValues(const ScalarType val);
    void setZero();

    int inc(const short dir) const { return grid_.inc(dir); }

    void assign(const GridFunc<ScalarType>& src, const char dis);
    void scatterFrom(const GridFunc<ScalarType>& src);

    double fmax();

    template <typename ScalarType2, typename MemorySpaceType>
    void getValues(ScalarType2*) const;

    void init_vect(ScalarType*, const char) const;
    void init_vect_shift(ScalarType* global_func) const;
    void gather(ScalarType*) const;
    void allGather(ScalarType*) const;

    short bc(const short dir) const
    {
        assert(bc_[dir] == 0 || bc_[dir] == 1 || bc_[dir] == 2);
        return bc_[dir];
    }
    short fully_periodic() const
    {
        return ((bc_[0] == 1) && (bc_[1] == 1) && (bc_[2] == 1));
    }
    void set_bc(const short px, const short py, const short pz)
    {
        assert(px == 0 || px == 1 || px == 2);
        assert(py == 0 || py == 1 || py == 2);
        assert(pz == 0 || pz == 1 || pz == 2);
        bc_[0]              = px;
        bc_[1]              = py;
        bc_[2]              = pz;
        updated_boundaries_ = false;
    }
    void setValuesNeumann(const ScalarType valuesNeumann[3])
    {
        valuesNeumann_[0]   = valuesNeumann[0];
        valuesNeumann_[1]   = valuesNeumann[1];
        valuesNeumann_[2]   = valuesNeumann[2];
        updated_boundaries_ = false;
    }
    ScalarType* uu(const int i = 0) const
    {
        assert(uu_ != 0);
        return uu_ + i;
    }
    bool updated_boundaries() const { return updated_boundaries_; }

    void set_updated_boundaries(const bool i) { updated_boundaries_ = i; }

    void copyFrom(GridFunc<ScalarType>* src)
    {
        memcpy(uu_, src->uu_, grid_.sizeg() * sizeof(ScalarType));
    }

    void copyFrom(GridFunc<ScalarType>& src)
    {
        memcpy(uu_, src.uu_, grid_.sizeg() * sizeof(ScalarType));
    }

    GridFunc<ScalarType>& operator=(const GridFunc<ScalarType>& func);
    GridFunc<ScalarType>& operator=(const ScalarType val);

    GridFunc<ScalarType> operator+(const GridFunc<ScalarType>& A);
    GridFunc<ScalarType> operator-(const GridFunc<ScalarType>& A);
    GridFunc<ScalarType> operator*(const double val);

    GridFunc<ScalarType>& operator+=(const GridFunc<ScalarType>& func);
    GridFunc<ScalarType>& operator+=(ScalarType alpha);

    GridFunc<ScalarType>& operator-=(const ScalarType alpha);

    GridFunc<ScalarType>& operator*=(const double alpha);
    GridFunc<ScalarType>& operator*=(const GridFunc<ScalarType>& B);

    GridFunc<ScalarType>& operator/=(const GridFunc<ScalarType>& B);

    void axpy(const ScalarType alpha, const GridFunc<ScalarType>& vv);
    void scal(const double alpha);
    void prod(const GridFunc<ScalarType>& A, const GridFunc<ScalarType>& B);
    void diff(const GridFunc<ScalarType>& A, const GridFunc<ScalarType>& B);

    void set_max(const ScalarType val);

    virtual ~GridFunc();

    int count_threshold(const ScalarType);

    void jacobi(const GridFunc<ScalarType>& v,
        const GridFunc<ScalarType>& epsilon, const double c0);
    void add_prod(
        const GridFunc<ScalarType>& v1, const GridFunc<ScalarType>& v2);
    void substract_prod(
        const GridFunc<ScalarType>& v1, const GridFunc<ScalarType>& v2);
    void init_radial(const double[3], const double, const short ftype = 0);
    void initTrigo3d(const short bc[3], const int n[3]);
    void initCos3d(const ScalarType k[3]);
    void print_radial(const char[]);
    void defaultTrade_boundaries();
    virtual void trade_boundaries() { defaultTrade_boundaries(); }
    void set_bc_func(GridFunc<ScalarType>* bc_func) { bc_func_ = bc_func; }

    void setBoundaryValues(const ScalarType, const bool direction[3]);
    void setBoundaryValues(
        const GridFunc<ScalarType>&, const bool direction[3]);
    void setBoundaryValuesNeumann(
        const ScalarType alpha[3], const bool direction[3]);

    bool isZero(const double tol = 1.e-16, const bool wghosts = false);
    void test_setBoundaryValues();
    template <typename ScalarType2>
    double gdot(const GridFunc<ScalarType2>&) const;
    double norm2() const;
    void extend3D(GridFunc<ScalarType>&);
    void restrict3D(GridFunc<ScalarType>&);
    void test_grid_transfer();
    void init_rand();
    void initTriLin(const ScalarType a[4], const bool wghosts = true);
    double integral() const;
    double get_average();
    double average0();

    template <typename ScalarType2>
    void assign(const ScalarType2* const, const char dis = 'd');

    void test_newgrid();
    void print(std::ostream&);
    void write_plt(const char[]) const;
    void write_global_x(const char str[]);
    void global_xyz_task0(ScalarType*);
    void write_xyz(std::ofstream&) const;
    void write_zyx(std::ofstream&) const;
    void write_global_xyz(std::ofstream&);
    void inv(void);
    void inv_sqrt(void);
    void sqrt_func(void);
    void smooth_by_coarsening(int);
    void add_bias(const double bias);
    double get_bias();

    short ghost_pt() const { return grid_.ghost_pt(); }
    int dim(const short i) const { return grid_.dim(i); }
    const PEenv& mype_env() const { return grid_.mype_env(); }
    int sizeg() const { return grid_.sizeg(); }
    int size() const { return grid_.size(); }
    const Grid& grid() const { return grid_; }

    GridFunc<ScalarType>& operator-=(const GridFunc<ScalarType>& func);
    void resetData()
    {
        memset(uu_, 0, grid_.sizeg() * sizeof(ScalarType));
        updated_boundaries_ = true;
    }
    void getCellCornersValues(
        const int i, const int j, const int k, double val[8]) const;

    void initiateExchangeNorthSouth();
    void finishExchangeNorthSouth();
    void initiateExchangeUpDown();
    void finishExchangeUpDown();
    void initiateExchangeEastWest();
    void finishExchangeEastWest();

    void setBoundaryValuesBeforeTrade();

    static void printTimers(std::ostream& os)
    {
        trade_bc_tm_.print(os);
        finishExchangeNorthSouth_tm_.print(os);
        finishExchangeUpDown_tm_.print(os);
        finishExchangeEastWest_tm_.print(os);
        extend3D_tm_.print(os);
        restrict3D_tm_.print(os);
        prod_tm_.print(os);
        gather_tm_.print(os);
        scatter_tm_.print(os);
        all_gather_tm_.print(os);
    }
};

template <typename ScalarType>
double dot(const GridFunc<ScalarType>& A, const GridFunc<ScalarType>& B)
{
    return A.gdot(B) * A.grid().vel();
}
template <typename ScalarType>
double norm(const GridFunc<ScalarType>& A)
{
    return A.norm2();
}

template <typename ScalarType>
GridFunc<ScalarType>* GridFunc<ScalarType>::bc_func_ = nullptr;

template <typename ScalarType>
Timer GridFunc<ScalarType>::trade_bc_tm_(
    "GridFunc::trade_bc_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::restrict3D_tm_(
    "GridFunc::restrict3D_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::extend3D_tm_(
    "GridFunc::extend3D_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::prod_tm_(
    "GridFunc::prod_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::gather_tm_(
    "GridFunc::gather_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::scatter_tm_(
    "GridFunc::scatter_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::all_gather_tm_(
    "GridFunc::all_gather_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::finishExchangeNorthSouth_tm_(
    "GridFunc::finishExNorthSouth_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::finishExchangeUpDown_tm_(
    "GridFunc::finishExUpDown_" + std::to_string(sizeof(ScalarType) * 8));
template <typename ScalarType>
Timer GridFunc<ScalarType>::finishExchangeEastWest_tm_(
    "GridFunc::finishExEastWest_" + std::to_string(sizeof(ScalarType) * 8));

} // namespace pb

#endif
