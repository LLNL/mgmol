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
#include "GridFuncInterface.h"
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

template <typename T>
class GridFunc : public GridFuncInterface
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

    void setBoundaryValuesNeumannX(const T alpha);
    void setBoundaryValuesNeumannY(const T alpha);
    void setBoundaryValuesNeumannZ(const T alpha);

    T one(const double r, const double a)
    {
        (void)r; // unused
        (void)a; // unused

        return 1.;
    }
    T ran0();
    T radial_func(const double r, const double a, const short ftype = 0);

    // allocate memory for function
    void alloc();

    void setup();

    void resizeBuffers();

protected:
    const Grid& grid_;

    // data storage
    std::unique_ptr<T> memory_;

    // raw pointer to data
    T* uu_;

    bool updated_boundaries_;
    short bc_[3];
    T valuesNeumann_[3];

    static GridFunc<T>* bc_func_;

    static std::vector<T> buf1_;
    static std::vector<T> buf2_;
    static std::vector<T> buf3_;
    static std::vector<T> buf4_;

public:
    // Constructors
    GridFunc(const Grid&, const short, const short, const short);

    // constructor with pointer to data allocation
    GridFunc(const Grid&, const short, const short, const short, T*,
        const bool updated_boundaries = false);

    // copy constructor
    GridFunc(const GridFunc<double>& A);
    GridFunc(const GridFunc<float>& A);

    // copy constructor on different grid
    GridFunc(const GridFunc<T>&, const Grid&);

    void setValues(const GridFunc<T>& src);
    void setValues(const T val);
    void setValues(const int n, const T* src, const int pos = 0);

    int inc(const short dir) const { return grid_.inc(dir); }

    void assign(const GridFunc<T>& src, const char dis);
    void scatterFrom(const GridFunc<T>& src);

    double fmax();

    template <typename T2, typename MemorySpaceType>
    void getValues(T2*) const;

    void init_vect(T*, const char) const;
    void init_vect_shift(T* global_func) const;
    void gather(T*) const;
    void allGather(T*) const;

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
    void setValuesNeumann(const T valuesNeumann[3])
    {
        valuesNeumann_[0]   = valuesNeumann[0];
        valuesNeumann_[1]   = valuesNeumann[1];
        valuesNeumann_[2]   = valuesNeumann[2];
        updated_boundaries_ = false;
    }
    T* uu(const int i = 0) const
    {
        assert(uu_ != 0);
        return uu_ + i;
    }
    bool updated_boundaries() const { return updated_boundaries_; }

    void set_updated_boundaries(const bool i) { updated_boundaries_ = i; }

    void copyFrom(GridFunc<T>* src)
    {
        memcpy(uu_, src->uu_, grid_.sizeg() * sizeof(T));
    }

    void copyFrom(GridFunc<T>& src)
    {
        memcpy(uu_, src.uu_, grid_.sizeg() * sizeof(T));
    }

    GridFunc<T>& operator=(const GridFunc<T>& func);
    GridFunc<T>& operator=(const T val);

    GridFunc<T> operator+(const GridFunc<T>& A);
    GridFunc<T> operator-(const GridFunc<T>& A);
    GridFunc<T> operator*(const double val);

    GridFunc<T>& operator+=(const GridFunc<T>& func);
    GridFunc<T>& operator+=(T alpha);

    GridFunc<T>& operator-=(const T alpha);

    GridFunc<T>& operator*=(const double alpha);
    GridFunc<T>& operator*=(const GridFunc<T>& B);

    GridFunc<T>& operator/=(const GridFunc<T>& B);

    void axpy(const double alpha, const GridFunc<T>& vv);
    void scal(const double alpha);
    void prod(const GridFunc<T>& A, const GridFunc<T>& B);
    void diff(const GridFunc<T>& A, const GridFunc<T>& B);

    void set_max(const T val);

    ~GridFunc() override;

    int count_threshold(const T);

    void jacobi(
        const GridFunc<T>& v, const GridFunc<T>& epsilon, const double c0);
    void add_prod(const GridFunc<T>& v1, const GridFunc<T>& v2);
    void substract_prod(const GridFunc<T>& v1, const GridFunc<T>& v2);
    void init_radial(const double[3], const double, const short ftype = 0);
    void initTrigo3d(const short bc[3], const int n[3]);
    void initCos3d(const T k[3]);
    void print_radial(const char[]);
    void defaultTrade_boundaries();
    virtual void trade_boundaries() { defaultTrade_boundaries(); }
    void set_bc_func(GridFunc<T>* bc_func) { bc_func_ = bc_func; }

    void setBoundaryValues(const T, const bool direction[3]);
    void setBoundaryValues(const GridFunc<T>&, const bool direction[3]);
    void setBoundaryValuesNeumann(const T alpha[3], const bool direction[3]);

    bool isZero(const double tol = 1.e-16, const bool wghosts = false);
    void test_setBoundaryValues();
    /* gdot is split into two for consistency between float-double and
     * double-float calls to gdot. See function definition for more details.
     */
    double gdot(const GridFunc<double>&) const;
    double gdot(const GridFunc<float>&) const;
    double norm2() const;
    void extend3D(GridFunc<T>&);
    void restrict3D(GridFunc<T>&);
    void test_grid_transfer();
    void init_rand();
    void initTriLin(const T a[4], const bool wghosts = true);
    double integral() const;
    double get_average();
    double average0();

    template <typename T2>
    void assign(const T2* const, const char dis = 'd');

    void test_newgrid();
    void print(std::ostream&);
    void write_plt(const char[]) const;
    void write_global_x(const char str[]);
    void global_xyz_task0(T*);
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

    GridFunc<T>& operator-=(const GridFunc<T>& func);
    void resetData()
    {
        memset(uu_, 0, grid_.sizeg() * sizeof(T));
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
};

template <typename T>
double dot(const GridFunc<T>& A, const GridFunc<T>& B)
{
    return A.gdot(B) * A.grid().vel();
}
template <typename T>
double norm(const GridFunc<T>& A)
{
    return A.norm2();
}

template <typename T>
GridFunc<T>* GridFunc<T>::bc_func_ = nullptr;

} // namespace pb

#endif
