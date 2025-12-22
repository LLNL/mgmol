// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_RADIALMESHFUNCTION_H
#define MGMOL_RADIALMESHFUNCTION_H

#include <mpi.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "MGmol_blas1.h"

class RadialMeshFunction
{
protected:
    std::vector<double> x_;
    std::vector<std::vector<double>> y_; // values at x_
    double drlin_;
    double invdr_;

public:
    RadialMeshFunction()
    {
        // initialize with invalid values
        drlin_ = -1.;
        invdr_ = -1.;
    };

    virtual ~RadialMeshFunction() {};

    RadialMeshFunction(const std::vector<double>& x)
    {
        drlin_ = x[1] - x[0];
        assert(drlin_ > 1.e-15);
        invdr_ = 1. / drlin_;

        x_ = x;
        y_.resize(1);
        y_[0].resize(x_.size());
    };

    RadialMeshFunction(
        const std::vector<double>& x, const std::vector<std::vector<double>>& y)
    {
        drlin_ = x[1] - x[0];
        assert(drlin_ > 1.e-15);
        invdr_ = 1. / drlin_;

        assert(x.size() == y[0].size());
        x_ = x;
        y_ = y;
    }

    RadialMeshFunction(
        const std::vector<double>& x, const std::vector<double>& y)
    {
        drlin_ = x[1] - x[0];
        assert(drlin_ > 1.e-15);
        invdr_ = 1. / drlin_;

        assert(x.size() == y.size());
        x_ = x;
        y_.resize(1);
        y_[0] = y;
    }

    void resize(const int n, const int ny = 1)
    {
        x_.resize(n);
        y_.resize(ny);
        for (int j = 0; j < ny; j++)
            y_[j].resize(n);
        x_[0]    = 1.e152;
        y_[0][0] = 1.e152;
    }

    std::vector<double>& y(const int i = 0)
    {
        assert(i < (int)y_.size());
        assert(y_[i].size() > 0);
        return y_[i];
    }
    double y(const int j, const int i) const
    {
        assert(j < (int)y_.size());
        assert(i < (int)y_[j].size());
        return y_[j][i];
    }
    const std::vector<double>& x() const { return x_; }
    double x(const int i) const
    {
        assert(i < (int)x_.size());
        return x_[i];
    }

    void scale(const double alpha, const int j = 0)
    {
        int n    = (int)y_[j].size();
        int ione = 1;
        DSCAL(&n, &alpha, &y_[j][0], &ione);
    }

    double norm2(const int j = 0) { return radintf2(j); }

    void bcast(MPI_Comm comm, const int root);

    void print(std::ofstream& tfile, const int j = 0) const
    {
        assert(j < (int)y_.size());
        const int n = (int)x_.size();
        for (int i = 0; i < n; i += 3)
        {
            tfile << x_[i];
            tfile << "\t" << y_[j][i];
            tfile << std::endl;
        }
    }

    void assign(
        const std::vector<double>& x, const std::vector<std::vector<double>>& y)
    {
        assert(x.size() == y[0].size());
        x_ = x;
        y_ = y;

        drlin_ = x_[1] - x_[0];
        assert(drlin_ > 1.e-15);
        invdr_ = 1. / drlin_;
    }

    void assign(const std::vector<double>& x, const std::vector<double>& y)
    {
        assert(x.size() == y.size());
        x_ = x;
        y_.resize(1);
        y_[0] = y;

        drlin_ = x_[1] - x_[0];
        assert(drlin_ > 1.e-15);
        invdr_ = 1. / drlin_;
    }

    void assign(const std::vector<std::vector<double>>& y, const double drlin)
    {
        drlin_ = drlin;
        assert(drlin_ > 1.e-15);
        invdr_ = 1. / drlin_;
        y_     = y;
        int n  = (int)y_[0].size();
        x_.resize(n);
        for (int i = 0; i < n; i++)
            x_[i] = i * drlin_;

        assert(drlin_ > 1.e-15);
    }

    void assign(
        const std::vector<double>& y, const double drlin, const int j = 0)
    {
        drlin_ = drlin;
        assert(drlin_ > 1.e-15);
        invdr_ = 1. / drlin_;
        if (y_.size() == 0) y_.resize(1);
        y_[j] = y;
        int n = (int)y_[0].size();
        if (j > 0) assert(n == (int)y.size());
        x_.resize(n);
        for (int i = 0; i < n; i++)
            x_[i] = i * drlin_;
    }

    double gauss_cutoff(const double r, const double rcut);
    double gcutoff(const double g, const double gcut);
    void gauss_filter(const double, const int j = 0);
    void read(const int, const int, std::ifstream* tfile);
    double radint(const int j = 0);
    double radintf2(const int j = 0);
    void inv_ft(const std::vector<double>& rgrid, std::vector<double>& rfunc,
        const int lval, const int j = 0);
    void compute_ft(const std::vector<double>& ggrid, std::vector<double>& gcof,
        const int lval, std::ofstream* tfile, const int j = 0);
    void rft(RadialMeshFunction& filt_func, const int lval, const double hmax,
        const bool, const int j = 0, const bool printFlag = false);

    void printLocalPot(
        const std::string& name, const short ll, std::ofstream* tfile);
};

#endif
