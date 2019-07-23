// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_POTENTIALS_H
#define MGMOL_POTENTIALS_H

#include "Rho.h"
#include "TriCubic.h"

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

class Ions;
class Species;
template <class T>
class GridFunc;

class Potentials
{

    int size_;
    int gdim_[3];
    int dim_[3];
    bool diel_;
    double mix_;

    double scf_dvrho_;
    double scf_dv_;

    double background_charge_;
    double charge_in_cell_;
    double ionic_charge_;

    bool vh_frozen_;

    std::vector<POTDTYPE> vtot_;
    std::vector<POTDTYPE> vtot_old_;

    std::vector<POTDTYPE> vepsilon_;
    std::vector<POTDTYPE> vh_rho_;
    std::vector<POTDTYPE> vxc_rho_;

    // nuclei local potential
    std::vector<POTDTYPE> v_nuc_;

    // external potential (read from input)
    std::vector<POTDTYPE> v_ext_;
#ifdef HAVE_TRICUBIC
    pb::TriCubic<POTDTYPE>* vext_tricubic_;
#endif

    std::vector<POTDTYPE> v_comp_;
    std::vector<RHODTYPE> rho_comp_;

    std::vector<POTDTYPE> dv_;

    int itindex_vxc_;
    int itindex_vh_;

    short verbosity_level_;

    // filenames for various input potentials and their type
    // 0 = radial potential
    // 1 = filtered radial potential
    // 2 = xyz potential
    std::vector<std::string> pot_filenames_;
    std::vector<short> pot_types_;

    void evalNormDeltaVtotRho(const vector<vector<RHODTYPE>>& rho);

    void initializeRadialDataOnMesh(const Vector3D& position, const Species& sp);

public:
    Potentials(const bool vh_frozen = false);

    ~Potentials();

    void setVerbosity(const short vlevel) { verbosity_level_ = vlevel; }

    void registerName(const std::string filename, const short flag)
    {
        pot_filenames_.push_back(filename);
        pot_types_.push_back(flag);
    }

    void writeNames(ostream& os)
    {
        vector<string>::iterator p = pot_filenames_.begin();
        while (p != pot_filenames_.end())
        {
            os << " Potential file: " << *p << endl;
            p++;
        }
    }

    short pot_type(const int isp) const
    {
        assert(isp < (int)pot_types_.size());
        return pot_types_[isp];
    }

    int getIterativeIndex() const
    {
        assert(itindex_vxc_ >= 0);
        assert(itindex_vh_ >= 0);
        assert(itindex_vxc_ == itindex_vh_);
        return itindex_vh_;
    }

    void turnOnDiel() { diel_ = true; }

    int size() const { return size_; }
    bool vh_frozen() const { return vh_frozen_; }
    void freeze_vh() { vh_frozen_ = true; }

    double scf_dvrho(void) const { return scf_dvrho_; }
    double scf_dv(void) const { return scf_dv_; }
    POTDTYPE* vtot() { return &vtot_[0]; }
    POTDTYPE* vh_rho() { return &vh_rho_[0]; }
    POTDTYPE* vxc_rho() { return &vxc_rho_[0]; }
    RHODTYPE* rho_comp() { return &rho_comp_[0]; }

    const vector<POTDTYPE>& vnuc() const { return v_nuc_; }
    POTDTYPE* vnuc() { return &v_nuc_[0]; }
    POTDTYPE* vext() { return &v_ext_[0]; }
    POTDTYPE* vepsilon() { return &vepsilon_[0]; }

    void set_vcomp(const POTDTYPE val)
    {
        const int n = (int)v_comp_.size();
        for (int i = 0; i < n; i++)
            v_comp_[i] = val;
    }
    void axpVcompToVh(const double alpha);
    void axpVcomp(POTDTYPE* v, const double alpha);

    POTDTYPE vtot(const int i) { return vtot_[i]; }
    POTDTYPE vh_rho(const int i) { return vh_rho_[i]; }
    POTDTYPE vxc_rho(const int i) { return vxc_rho_[i]; }
    POTDTYPE vepsilon(const int i) { return vepsilon_[i]; }

    bool diel() const { return diel_; }

    double getChargeInCell() const { return charge_in_cell_; }

    void initWithVnuc();

    void getVofRho(vector<POTDTYPE>& vrho) const;

    double delta_v(const vector<vector<RHODTYPE>>&);
    double update(const vector<vector<RHODTYPE>>&);
    void update(const double);
    double max() const;
    double min() const;
    void readAll(vector<Species>& sp);

    template <typename T>
    void setVxc(const T* const vxc, const int iterativeIndex);
    void setVh(const POTDTYPE* const vh, const int iterativeIndex);
    void setVh(const pb::GridFunc<POTDTYPE>& vh, const int iterativeIndex);

    void initialize(Ions& ions);
    void rescaleRhoComp();
    double getCharge(RHODTYPE* rho);
    void initBackground(Ions& ions);
    void addBackgroundToRhoComp();

#ifdef HAVE_TRICUBIC
    void readExternalPot(const string filename, const short type);
    void setupVextTricubic();
    bool withVext() const;
    void getGradVext(const double[3], double[3]) const;
    void getValVext(const vector<double>& r, vector<double>& val) const;
#endif
};

#endif
