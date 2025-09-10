// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_KBPROJECTORSPARSE_H
#define MGMOL_KBPROJECTORSPARSE_H

#include "KBprojector.h"
#include "global.h"
#include "mputils.h"

#include <cassert>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

#define CHECK_NORM 0

// one KB projector
class KBprojectorSparse : public KBprojector
{
    // work arrays (1 for each thread)
    static std::vector<std::vector<KBPROJDTYPE>> work_nlindex_;

    static std::vector<std::vector<KBPROJDTYPE>> work_proj_;

    // pointers to projectors for each iloc, l, p, m
    std::vector<std::vector<std::vector<std::vector<KBPROJDTYPE*>>>>
        ptr_projector_;

    // storage for each 'iloc' i sstored in a std::vector
    std::vector<std::vector<KBPROJDTYPE>> projectors_storage_;

    std::vector<int> size_nl_;

    short range_kbproj_;

    // projector start
    short kb_proj_start_index_[3];
    double kb_proj_start_[3];

    std::vector<short> proj_indices_[3];

    // list of nodes where the non-local projector is non-zero on this PE
    // for each subdomain "iloc"
    std::vector<std::vector<int>> nlindex_;

    bool** is_in_domain_;

#if CHECK_NORM
    std::vector<std::vector<double>> norm2_;
#endif

    // private functions
    void allocateProjectors(const short iloc, const int icount);

    void setProjectors(const short iloc, const int icount);

    void setSProjector(const short iloc, const int icount);
    void setPProjector(const short iloc, const int icount);
    void setDProjector(const short iloc, const int icount);
    void setFProjector(const short iloc, const int icount);

    void setNLindex(
        const short iloc, const int size, const std::vector<int>& pvec);

    void setKBProjStart();
    void setProjIndices(const short dir);
    int get_index_array(std::vector<int>& pvec, const short iloc,
        const short index_low[3], const short index_high[3]);
    bool overlapWithBox(const short index_low[3], const short index_high[3]);

    // get projector l,m in pieces corresponding to subdivisions
    const KBPROJDTYPE* getProjector(
        const short iloc, const short l, const short p, const short m) const
    {
        assert(iloc < subdivx_);
        assert(l < static_cast<int>(ptr_projector_[iloc].size()));
        assert(p < static_cast<int>(ptr_projector_[iloc][l].size()));
        assert(m < static_cast<int>(ptr_projector_[iloc][l][p].size()));

        return ptr_projector_[iloc][l][p][m];
    }
    double dotPsi(
        const short iloc, const short l, const short p, const short m) const
    {
        assert(iloc < subdivx_);
        assert(l < static_cast<int>(ptr_projector_[iloc].size()));
        assert(p < static_cast<int>(ptr_projector_[iloc][l].size()));
        assert(m < static_cast<int>(ptr_projector_[iloc][l][p].size()));
        assert(static_cast<unsigned int>(omp_get_thread_num())
               < work_nlindex_.size());

        return LinearAlgebraUtils<MemorySpace::Host>::MPdot(size_nl_[iloc],
            &work_nlindex_[omp_get_thread_num()][0],
            (ptr_projector_[iloc][l][p][m]));
    }
    bool setIndexesAndProjectors();
    void setPtrProjectors();

    void getProjectors(
        const short iloc, std::vector<const KBPROJDTYPE*>& projectors) const;

    template <typename T>
    void axpySKetT(const short iloc, const double alpha, T* const dst) const;

    template <typename T>
    void axpyKetT(
        const short iloc, const std::vector<double>& alpha, T* const dst) const;

public:
    KBprojectorSparse(const Species& sp);
    KBprojectorSparse(const KBprojectorSparse& kbp);

    virtual ~KBprojectorSparse() { clear(); }

    void clear() override
    {
        if (is_in_domain_ != nullptr)
        {
            for (short iloc = 0; iloc < subdivx_; iloc++)
                delete[] is_in_domain_[iloc];
            delete[] is_in_domain_;
            is_in_domain_ = nullptr;
        }

        for (short dir = 0; dir < 3; dir++)
            proj_indices_[dir].clear();
    }

    // setup data that depends on atomic position
    void setup(const double center[3]) override;

    double maxRadius() const override;

    bool overlapPE() const override;
    void registerPsi(const short iloc, const ORBDTYPE* const psi) override;

    bool overlaps(const short iloc) const override
    {
        return (size_nl_[iloc] > 0);
    }
    double dotPsi(const short iloc, const short index) const override;

    // axpySket for templated destination type
    void axpySKet(
        const short iloc, const double alpha, double* const dst) const override
    {
        axpySKetT(iloc, alpha, dst);
    }
    void axpySKet(
        const short iloc, const double alpha, float* const dst) const override
    {
        axpySKetT(iloc, alpha, dst);
    }

    void axpyKet(const short iloc, const std::vector<double>& alpha,
        double* const dst) const override
    {
        axpyKetT(iloc, alpha, dst);
    }
    void axpyKet(const short iloc, const std::vector<double>& alpha,
        float* const dst) const override
    {
        axpyKetT(iloc, alpha, dst);
    }

    void getKBsigns(std::vector<short>& kbsigns) const override;
    void getKBcoeffs(std::vector<double>& coeffs) const override;
};

#endif
