// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef SpreadsAndCenters_H
#define SpreadsAndCenters_H

#include "VariableSizeMatrix.h"
#include "Vector3D.h"

#include <vector>

template <class T>
class SpreadsAndCenters
{
private:
    // number of global ids locally known
    int ngids_;

    // computational domain
    Vector3D origin_;
    Vector3D cell_;

    std::vector<int> localRows_;
    std::vector<int> localRowGid_;

    // cosine and sine matrices for x,y,z directions
    // for gids associated with overlapping regions
    std::vector<std::vector<double>> r_;

    // cos^2 and sin^2 matrices for x,y,z directions
    // for gids associated with overlapping regions
    std::vector<std::vector<double>> r2_;

    // gids of overlapping data
    std::vector<int> gids_;

    std::vector<int> localRowsSq_;
    std::vector<int> localRowSqGid_;

    // private functions
    double computeSpread2(int i, int j) const;
    double computeSpread2(int i) const;
    double computeSpread2(void) const;

    void setData(VariableSizeMatrix<sparserow>& mat,
        const std::vector<int>& gids, const std::vector<int>& localRowGid,
        std::vector<std::vector<double>>& matr);

public:
    virtual ~SpreadsAndCenters() { }

    void computeCenters(std::vector<Vector3D>& centers) const
    {
        assert(ngids_ >= 0);

        centers.resize(ngids_);
        for (int i = 0; i < ngids_; i++)
            centers[i] = computeCenter(i);
    }

    void computeSpreads(std::vector<double>& spreads) const
    {
        assert(ngids_ >= 0);

        spreads.resize(ngids_);
        for (int i = 0; i < ngids_; i++)
            spreads[i] = computeSpread(i);
    }

    void computeSpreads2(std::vector<double>& spreads2) const
    {
        assert(ngids_ >= 0);

        spreads2.resize(ngids_);
        for (int i = 0; i < ngids_; i++)
            spreads2[i] = computeSpread2(i);
    }

    void computeSpreads2(std::vector<float>& spreads2) const
    {
        assert(ngids_ >= 0);

        spreads2.resize(ngids_);
        for (int i = 0; i < ngids_; i++)
            spreads2[i] = computeSpread2(i);
    }

    double computeSpread(int i) const;
    double computeSpread(void) const;
    Vector3D computeCenter(const int gid) const;

    void computeLocalSpreads(std::vector<float>& spreads);
    void computeLocalSpreads2(std::vector<float>& spreads);
    void getLocalCenters(std::vector<Vector3D>& centers);
    void getLocalGids(std::vector<int>& lindex);

    // total volume of spheres of radius "spread"
    double volume() const;

    void print(std::ostream& os, const int root = 0) const;
    void printGlobal(std::ostream& os, const int root = 0) const;
    void printLocal(std::ostream& os, const int root = 0) const;
    void printStats(std::ostream& os) const;

    double computeDistance(const int st1, const int st2) const;

    // set data for cosine and sine matrices
    void setSinCosData(VariableSizeMatrix<sparserow>& mat,
        const std::vector<int>&, const std::vector<int>& localRowGid);
    void setSinCosData(std::vector<std::vector<double>>& a, const int n);

    SpreadsAndCenters(const Vector3D& orig, const Vector3D& ll)
        : origin_(orig), cell_(ll)
    {
        ngids_ = -1; // not set yet

        r_.clear();
        r2_.clear();
        localRows_.clear();
        localRowGid_.clear();
        localRowsSq_.clear();
        localRowSqGid_.clear();
    };

    const std::vector<int>& getGids() const { return gids_; }

    void computePositionMatrix(T& orbitals, T& work_orbitals);
    void computePositionMatrix(const T& orbitals);

    void computeSinCos(const T& orbitals);
    void computeSinCosSquare(const T& orbitals);
    void computeSinCosSquare1D(const T& orbitals, const int dir);
    void computeSinCos2states(const T& orbitals, const int st1, const int st2);
    void computeSinCosDiag2states(
        const T& orbitals, const int st1, const int st2);

    void computeSinCos1D(const T& orbitals, const int dir);
    void computeSinCos(const T& orbitals1, const T& orbitals2);
    void computeSinCosDiag(const T& orbitals, const bool normalized_functions);
};

#endif
