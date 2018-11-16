// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

// $Id$
#ifndef SPECIES_H
#define SPECIES_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "MPIdata.h"
#include "RadialInter.h"

#define SCMASS 1822.89

// Species structure
class Species
{
private:
    // Name of the atom
    std::string name_;

    // Atomic number
    short atomic_number_;

    short zion_;

    // Number of radial grid points in input file
    int n_rad_points_;

    // Gaussian charge parameter
    double rc_;

    // type of pseudopotential
    // 0: KB
    // 1: KB with different convention for radial functions
    // 2: ONCVPSP
    short type_flag_;

    // max. radius for which v_loc != 0
    double lradius_;

    // max. radius for which v_nl != 0
    double nlradius_;

    // size of the grid to represent nl potentials
    short dim_nl_;

    // size of the grid to represent local potentials
    short dim_l_;

    // Unsigned normalization coefficients for the KB projectors
    std::vector<std::vector<double>> ekb_;

    std::vector<std::vector<double>> kb_coeff_;

    // Signs of kb_norm
    std::vector<std::vector<short>> kb_sign_;

    // Max. angular momentum l for non-local projectors
    short max_l_;

    // l-value for local pseudopotential
    short llocal_;

    // multiplicity for each l
    std::vector<short> multiplicity_;

    // Mass
    double mass_;

    // Radial pseudopotentials from input file
    std::vector<RadialInter> input_kbp_;

    // Number of potentials
    short num_potentials_;

    // Filtered Kleinman-Bylander Projectors on linear interpolation grid
    std::vector<std::vector<RadialInter>> kbp_;

    // Linear interpolation for the compensated filtered local potential
    RadialInter local_pot_;

    // parameters for Goedecker pseudopotentials
    double h1s_;
    double h2s_;
    double h1p_;

    // MPI communicator
    MPI_Comm comm_;

    void assignLocalPot(
        const std::vector<double>& x, const std::vector<double>& v)
    {
        local_pot_.assign(x, v);
    }
    std::vector<double>& ref_psi(const int l)
    {
        assert(l < num_potentials_);
        return input_kbp_[l].y(1);
    }
    void setKBcoeff(const short, const short, const double);
    void gauss_filter_local_pot(const double, const bool printFlag);
    void gauss_filter_kbp(
        const short, const short, const double, const bool printFlag);

    void assign_kbp(const short l, const short p, const std::vector<double>& x,
        const std::vector<double>& kbproj)
    {
        assert(l < kbp_.size());
        assert(p < kbp_[l].size());
        kbp_[l][p].assign(x, kbproj);
    }

    void initLocalPotential(const bool flag_filter, const double hmax,
        std::ofstream* tfile, const bool printFlag);
    void initNonlocalKBPotentials(const bool flag_filter, const double hmax,
        std::ofstream* tfile, const bool printFlag);
    void initNonlocalGTHPotentials(const bool flag_filter, const double hmax,
        std::ofstream* tfile, const bool printFlag);
    void initNonlocalMultiProjectorsPotentials(const bool flag_filter,
        const double hmax, std::ofstream* tfile, const bool divide_by_r,
        const bool printFlag);

    void readRadialKBpotentials(std::ifstream* tfile);
    void readRadialGTHpotentials(std::ifstream* tfile);

    void print_kbp(std::ofstream& tfile, const short l = 0, const short p = 0)
    {
        assert(l < (short)kbp_.size());
        assert(p < (short)kbp_[l].size());
        kbp_[l][p].print(tfile, 0);
        tfile << std::endl;
    }

public:
    Species(MPI_Comm comm)
    {
#ifdef DEBUG
        (*MPIdata::sout) << " Construct 1 undefined Species" << endl;
#endif

        name_           = "undefined";
        zion_           = 0;
        n_rad_points_   = 0;
        rc_             = 1.;
        lradius_        = 0.;
        nlradius_       = 0.;
        dim_nl_         = -1;
        dim_l_          = -1;
        max_l_          = -1;
        llocal_         = 0;
        mass_           = 0.;
        num_potentials_ = 0;
        comm_           = comm;
    }

    ~Species() {}

    unsigned short getAtomicNumber() const { return atomic_number_; }
    double getMass() const { return mass_ * SCMASS; } // SCMASS = 1822.89
    bool isMass28() const { return (fabs(mass_ - 28.) < 0.1); }
    bool isMassLargerThan1() const { return (mass_ > 1.1); }
    const RadialInter& getRadialKBP(const short l, const short p) const
    {
        assert(l < (short)kbp_.size());
        assert(p < (short)kbp_[l].size());
        return kbp_[l][p];
    }

    const RadialInter& local_pot() const { return local_pot_; }

    short zion() const { return zion_; }
    short llocal() const { return llocal_; }
    double rc() const { return rc_; }

    double lradius() const { return lradius_; }
    double nlradius() const { return nlradius_; }

    short max_l() const { return max_l_; }
    short dim_nl() const { return dim_nl_; }
    short dim_l() const { return dim_l_; }

    std::string name() const { return name_; }
    void read_1species(const std::string&);
    void print(std::ostream&) const;

    void set_dim_nl(const double hgrid)
    {
        assert(nlradius_ > 0.);
        assert(hgrid > 0.);

        dim_nl_ = (short)(2. * nlradius_ / hgrid) + 1;
        if ((dim_nl_ + 1) % 2) dim_nl_++;
    }

    void set_dim_l(const double hgrid)
    {
        assert(lradius_ > 0.);

        dim_l_ = (short)(2. * lradius_ / hgrid) + 1;
        if ((dim_l_ + 1) % 2) dim_l_++;
    }
    void initPotentials(const bool, const double, const bool);

    void getKBsigns(std::vector<short>& kbsigns) const;
    void getKBcoeffs(std::vector<double>& coeffs) const;

    short getMultiplicity(const short l) const { return multiplicity_[l]; }

    short getKBsign(const short l, const short p) const
    {
        return kb_sign_[l][p];
    }

    double getKBcoeff(const short l, const short p) const
    {
        return kb_coeff_[l][p];
    }
    // void getNLprojIndexes(std::vector<short>& indexes)const;

    void syncKBP(const int root);
};

#endif
