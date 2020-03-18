// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef GRIDMASK_H
#define GRIDMASK_H

#include "GridFunc.h"
#include "Timer.h"
#include "Vector3D.h"
#include "global.h"

#include <cassert>
#include <vector>

class SubCell;

class GridMask
{
private:
    static Timer init_tm_;

    static int ninstances_;
    static std::vector<SubCell*> sub_cell_; // SubCell for each level

    GridMask(const GridMask&);

    const pb::Grid& grid_;

    unsigned short nlevels_;
    unsigned short subdivx_;
    double delta_;
    std::vector<int> loc_numpt_; // for each level
    std::vector<int> subdim0_; // for each level

    Vector3D center_; // center of mask
    double radius_; // radius of mask

    // Localization mask for all the levels and subdomains
    std::vector<std::vector<std::vector<lmasktype>>> lmask_;

    // Values of mask_not_zero_:
    // -1 -> mask not needed (always 0 everywhere)
    // 0 -> mask 0 everywhere
    // 1 -> mask 1 everywhere
    // 2 -> 0 or 1 according to lmask_
    std::vector<std::vector<short>> mask_not_zero_; // levels, subdivx

protected:
    lmasktype* lmask(const unsigned short level, const unsigned short iloc)
    {
        return &lmask_[level][iloc][0];
    }

public:
    static void printTimers(std::ostream& os) { return init_tm_.print(os); }

    // Constructor
    GridMask(const unsigned short nclevels, const unsigned short subdivx,
        const pb::Grid& mygrid);

    GridMask& operator=(const GridMask& grid_mask);

    // Destructor
    virtual ~GridMask();

    int subdim0(const unsigned short ln) const { return subdim0_[ln]; }

    double center(const int i) const { return center_[i]; }
    double radius() const { return radius_; }
    int loc_numpt(const unsigned short ln) const { return loc_numpt_[ln]; }
    unsigned short subdivx() const { return subdivx_; }

    const Vector3D& center() const { return center_; }

    int init(const Vector3D&, const double rcut, const unsigned short iloc,
        const unsigned short level, lmasktype (*func)(const double));
    int init(const Vector3D&, const double rcut, const unsigned short level,
        lmasktype (*func)(const double));

    void assign(const unsigned short iloc, const unsigned short level,
        const int num, const int xnum, const std::vector<lmasktype>& val);
    void assign1(const unsigned short iloc, const unsigned short level);
    double delta() const { return delta_; }

    void clear(const unsigned short level)
    {
        assert(level < nlevels_);
        for (unsigned short i = 0; i < subdivx_; i++)
        {
            if (mask_not_zero_[level][i] == 2)
            {
                lmask_[level][i].resize(0);
            }
            mask_not_zero_[level][i] = 1;
        }
    }

    void clear()
    {
        for (unsigned short l = 0; l < nlevels_; l++)
            clear(l);
    }

    short mask_not_zero(
        const unsigned short level, const unsigned short iloc) const
    {
        assert(iloc < subdivx_);
        assert(iloc < (unsigned short)mask_not_zero_[level].size());
        assert(level < (unsigned short)mask_not_zero_.size());
        return mask_not_zero_[level][iloc];
    }

    template <typename T>
    void apply(pb::GridFunc<T>&, const unsigned short);

    // pure virtual functions
    virtual void apply(
        pb::GridFunc<float>&, const unsigned short, const unsigned short)
        = 0;
    virtual void apply(
        pb::GridFunc<double>&, const unsigned short, const unsigned short)
        = 0;

    virtual void apply(float*, const unsigned short, const unsigned short,
        const bool first_application = false)
        = 0;
    virtual void apply(double*, const unsigned short, const unsigned short,
        const bool first_application = false)
        = 0;

    void plot(const unsigned short, const int);
    bool overlap_on_pe(const GridMask& gm, const unsigned short level) const;

    bool maskIs0(const unsigned short level, const unsigned short iloc) const
    {
        return (mask_not_zero_[level][iloc] <= 0);
    }

    bool maskIs1(const unsigned short level, const unsigned short iloc) const
    {
        assert((unsigned short)mask_not_zero_.size() > level);
        assert((unsigned short)mask_not_zero_[level].size() > iloc);
        assert(mask_not_zero_[level][iloc] < 3);
        assert(mask_not_zero_[level][iloc] > -2);
        return (mask_not_zero_[level][iloc] == 1);
    }

    bool maskIs2(const unsigned short level, const unsigned short iloc) const
    {
        return (mask_not_zero_[level][iloc] == 2);
    }

    bool maskIsVoid(const unsigned short level, const unsigned short iloc) const
    {
        return (mask_not_zero_[level][iloc] == -1);
    }

    int offset(const unsigned short level, const unsigned short iloc) const
    {
        return loc_numpt_[level] * iloc;
    }

    template <typename T>
    void setZero(
        T* u, const unsigned short level, const unsigned short iloc) const
    {
        memset(u + loc_numpt_[level] * iloc, 0, loc_numpt_[level] * sizeof(T));
    }

    template <typename T>
    void setZero(pb::GridFunc<T>& gu, const unsigned short level,
        const unsigned short iloc) const
    {
        const int shift   = gu.ghost_pt();
        const int incx1   = gu.grid().inc(0);
        const int subdim0 = subdim0_[level];
        const int offset  = (shift + subdim0 * iloc) * incx1;
        assert(offset + incx1 * subdim0 < static_cast<int>(gu.grid().sizeg()));

        memset(gu.uu() + offset, 0, incx1 * subdim0 * sizeof(T));
    }

    template <typename T>
    T limitAbsValue(const T u, const lmasktype umask) const
    {
        if (u > umask) return (T)umask;
        if (u < -umask) return -(T)umask;
        return u;
    }

    template <typename T>
    void multiplyByMask(
        T* u, const unsigned short level, const unsigned short iloc) const;
    template <typename T>
    void cutWithMask(
        T* u, const unsigned short level, const unsigned short iloc) const;
};

#endif
