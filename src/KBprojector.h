// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_KBPROJECTOR_H
#define MGMOL_KBPROJECTOR_H

#include "Mesh.h"
#include "Species.h"
#include "global.h"

#include <cassert>
#include <vector>

// one KB projector
class KBprojector
{
protected:
    // reference to species associated to projector
    const Species& species_;

    short subdivx_;
    const short maxl_;
    const short llocal_;

    double center_[3];

    // multiplicity of projector for each "l"
    std::vector<short> multiplicity_;

    double h_[3];

    void initCenter(const double center[3])
    {
        for (short i = 0; i < 3; i++)
            center_[i] = center[i];
    }

public:
    KBprojector(const Species& sp)
        : species_(sp), maxl_(species_.max_l()), llocal_(species_.llocal())
    {
        for (short l = 0; l <= maxl_; l++)
            multiplicity_.push_back(sp.getMultiplicity(l));

        for (short l = 0; l <= maxl_; l++)
            assert(multiplicity_[l] > 0 || l == llocal_);

        Mesh* mymesh           = Mesh::instance();
        const pb::Grid& mygrid = mymesh->grid();
        for (short i = 0; i < 3; i++)
            h_[i] = mygrid.hgrid(i);
        subdivx_ = mymesh->subdivx();

        assert(llocal_ <= maxl_);
        assert(h_[0] > 0.);
        assert(h_[1] > 0.);
        assert(h_[2] > 0.);
        assert(species_.dim_nl() >= 0);
        assert(species_.dim_nl() < 1000);
        for (short l = 0; l <= maxl_; l++)
            assert(multiplicity_[l] > 0 || l == llocal_);
    }

    KBprojector(const KBprojector& kb)
        : species_(kb.species_),
          subdivx_(kb.subdivx_),
          maxl_(kb.maxl_),
          llocal_(kb.llocal_),
          multiplicity_(kb.multiplicity_)
    {
        for (short i = 0; i < 3; i++)
        {
            h_[i] = kb.h_[i];
        }

        assert(llocal_ <= maxl_);
        assert(species_.dim_nl() >= 0);
        assert(species_.dim_nl() < 1000);
    }

    ~KBprojector() { }

    virtual void clear() = 0;

    // setup data that depends on atomic position
    virtual void setup(const double center[3]) { initCenter(center); }

    virtual double maxRadius() const = 0;

    virtual bool overlapPE() const                                        = 0;
    virtual void registerPsi(const short iloc, const ORBDTYPE* const psi) = 0;

    virtual bool overlaps(const short iloc) const                    = 0;
    virtual double dotPsi(const short iloc, const short index) const = 0;

    // axpySket for templated destination type
    virtual void axpySKet(
        const short iloc, const double alpha, double* const) const
        = 0;
    virtual void axpySKet(
        const short iloc, const double alpha, float* const) const
        = 0;

    virtual void axpyKet(const short iloc, const std::vector<double>& alpha,
        double* const dst) const
        = 0;
    virtual void axpyKet(const short iloc, const std::vector<double>& alpha,
        float* const dst) const
        = 0;

    bool onlyOneProjector() const
    {
        return ((maxl_ == 1) && (llocal_ == 1) && (multiplicity_[0] == 1));
    }

    short nProjectors() const
    {
        short nproj = 0;
        assert(llocal_ < 5);
        for (short l = 0; l <= maxl_; l++)
            if (llocal_ != l) nproj += (2 * l + 1) * multiplicity_[l];
        return nproj;
    }

    short nProjectorsSubdomain() const
    {
        short nproj = 0;
        assert(llocal_ < 5);
        if (overlapPE()) nproj = nProjectors();
        return nproj;
    }

    virtual void getKBsigns(std::vector<short>& kbsigns) const  = 0;
    virtual void getKBcoeffs(std::vector<double>& coeffs) const = 0;
};

#endif
