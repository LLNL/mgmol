// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include "SpreadsAndCenters.h"
#include "Control.h"
#include "ExtendedGridOrbitals.h"
#include "LocGridOrbitals.h"
#include "MGmol_MPI.h"
#include "MPIdata.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace std;

static double fourthirdpi = 4. * M_PI / 3.;

////////////////////////////////////////////////////////////////////////////////
template <class T>
Vector3D SpreadsAndCenters<T>::computeCenter(const int index) const
{
    assert(index >= 0 && index < ngids_);
    assert(index < (int)r_[0].size());

    const double itwopi = 0.5 * M_1_PI;
    double cx           = r_[0][index];
    double sx           = r_[1][index];
    // get x between -0.5*cell and +0.5*cell
    double x = cell_[0] * itwopi * atan2(sx, cx);
    // get x 0. and cell
    if (x < 0.) x += cell_[0];
    // shift by origin of cell
    x += origin_[0];

    double cy = r_[2][index];
    double sy = r_[3][index];
    double y  = cell_[1] * itwopi * atan2(sy, cy);
    if (y < 0.) y += cell_[1];
    y += origin_[1];

    double cz = r_[4][index];
    double sz = r_[5][index];
    double z  = cell_[2] * itwopi * atan2(sz, cz);
    if (z < 0.) z += cell_[2];
    z += origin_[2];

    return Vector3D(x, y, z);
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double SpreadsAndCenters<T>::computeSpread2(int i) const
{
    assert((i >= 0) && (i < ngids_));
    return computeSpread2(i, 0) + computeSpread2(i, 1) + computeSpread2(i, 2);
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double SpreadsAndCenters<T>::computeSpread(int i) const
{
    assert(computeSpread2(i) >= 0.);
    return sqrt(computeSpread2(i));
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double SpreadsAndCenters<T>::computeSpread2(void) const
{
    double sum = 0.0;
    for (int i = 0; i < ngids_; i++)
        sum += computeSpread2(i);
    return sum;
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double SpreadsAndCenters<T>::computeSpread(void) const
{
    assert(computeSpread2() >= 0.);
    return sqrt(computeSpread2());
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double SpreadsAndCenters<T>::computeSpread2(int i, int j) const
{
    assert(i >= 0 && i < ngids_);
    assert(j >= 0 && j < 3);
    assert(cell_[j] > 0.);

    const double lby2pi = 0.5 * M_1_PI * cell_[j];
    double c            = r_[2 * j][i];
    double s            = r_[2 * j + 1][i];
    return lby2pi * lby2pi - c * c - s * s;
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
double SpreadsAndCenters<T>::volume() const
{
    double vol = 0.;
    for (int i = 0; i < ngids_; i++)
    {
        double spreadi = computeSpread(i);
        vol += spreadi * spreadi * spreadi;
    }
    return vol * fourthirdpi;
}

////////////////////////////////////////////////////////////////////////////////
template <class T>
void SpreadsAndCenters<T>::printGlobal(ostream& os, const int root) const
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    int ngids       = localRows_.size();
    mmpi.allreduce(&ngids, 1, MPI_SUM);

    vector<double> lspreads;
    vector<double> gspreads(ngids);
    vector<int> lindex;
    vector<int> gindex(ngids);
    vector<double> lcenters;
    vector<double> gcenters(3 * ngids);

    vector<int>::const_iterator it = localRowGid_.begin();
    int i                          = 0;
    while (it != localRowGid_.end())
    {
        Vector3D centeri(computeCenter(localRows_[i]));
        lcenters.push_back(centeri[0]);
        lcenters.push_back(centeri[1]);
        lcenters.push_back(centeri[2]);

        lspreads.push_back(computeSpread(localRows_[i]));
        lindex.push_back(*it);

        it++;
        i++;
    }
    // gather data
    mmpi.gatherV(lspreads, gspreads, root);
    mmpi.gatherV(lcenters, gcenters, root);
    mmpi.gatherV(lindex, gindex, root);
    if (gindex.size() == 0) return;

    if (mmpi.mypeSpin() == root)
    {
        map<int, double> spread;
        map<int, Vector3D> center;

        vector<int>::iterator gid = gindex.begin();
        int i                     = 0;
        while (gid != gindex.end())
        {
            Vector3D centeri(
                gcenters[3 * i], gcenters[3 * i + 1], gcenters[3 * i + 2]);
            center.insert(std::pair<int, Vector3D>(*gid, centeri));
            spread.insert(std::pair<int, double>(*gid, gspreads[i]));

            gid++;
            i++;
        }

        os << endl << " Orbitals centers and spreads " << endl << endl;

        map<int, double>::const_iterator spread_id   = spread.begin();
        map<int, Vector3D>::const_iterator center_id = center.begin();
        while (spread_id != spread.end())
        {
            os << "&& " << setw(4) << (spread_id->first + 1) << "   ";
            os.setf(ios::fixed, ios::floatfield);
            os.setf(ios::right, ios::adjustfield);
            os << setw(10) << setprecision(3) << (center_id->second)[0] << " "
               << setw(10) << setprecision(3) << (center_id->second)[1] << " "
               << setw(10) << setprecision(3) << (center_id->second)[2]
               << "         ";
            os << setw(10) << setprecision(3) << spread_id->second << endl;

            spread_id++;
            center_id++;
        }
    }
}

template <class T>
void SpreadsAndCenters<T>::printLocal(ostream& os, const int root) const
{
    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    if (mmpi.mypeSpin() == root)
    {
        os << endl
           << " Overlapping Orbitals centers and spreads " << endl
           << endl;
        for (int i = 0; i < ngids_; i++)
        {
            Vector3D centeri(computeCenter(i));
            const double spreadi = computeSpread(i);
            if (onpe0)
            {
                os << "&& " << setw(4) << i + 1 << "   ";
                os.setf(ios::fixed, ios::floatfield);
                os.setf(ios::right, ios::adjustfield);
                os << setw(10) << setprecision(3) << centeri[0] << " "
                   << setw(10) << setprecision(3) << centeri[1] << " "
                   << setw(10) << setprecision(3) << centeri[2] << "         ";
                os << setw(10) << setprecision(3) << spreadi << endl;
            }
        }
    }
}

template <class T>
void SpreadsAndCenters<T>::print(ostream& os, const int root) const
{
    Control& ct = *(Control::instance());
    if (ct.numst < 512 || ct.verbose > 2)
        printGlobal(os, root);
    else
        printLocal(os, root);
}
////////////////////////////////////////////////////////////////////////////////
template <class T>
void SpreadsAndCenters<T>::printStats(ostream& os) const
{
    if (onpe0 && ngids_ > 0)
        os << endl << " Orbitals spreads (statistics)" << endl << endl;
    double min_spread = 100000.;
    double max_spread = 0.;
    double avg_spread = 0.;

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    const int nids  = localRows_.size();
    for (int igid = 0; igid < nids; igid++)
    {
        const double spreadi = computeSpread(localRows_[igid]);
        avg_spread += spreadi;
        min_spread = min(spreadi, min_spread);
        max_spread = max(spreadi, max_spread);
    }

    double tmp;
    mmpi.allreduce(&avg_spread, &tmp, 1, MPI_SUM);
    avg_spread = tmp;

    mmpi.allreduce(&min_spread, &tmp, 1, MPI_MIN);
    min_spread = tmp;

    mmpi.allreduce(&max_spread, &tmp, 1, MPI_MAX);
    max_spread = tmp;

    int sum_nids = nids;
    mmpi.allreduce(&sum_nids, 1, MPI_SUM);
    avg_spread /= sum_nids;

    if (onpe0)
    {
        os << "&& ";
        os.setf(ios::fixed, ios::floatfield);
        os.setf(ios::right, ios::adjustfield);
        os << setw(10) << setprecision(3) << "Min. " << min_spread << ", Max. "
           << max_spread << ", Avg. " << avg_spread << endl;
    }
}

// get spreads for functions centered in subdomain
template <class T>
void SpreadsAndCenters<T>::computeLocalSpreads(vector<float>& spreads)
{
    spreads.clear();

    const int nids = localRows_.size();
    for (int igid = 0; igid < nids; igid++)
    {
        spreads.push_back(computeSpread(localRows_[igid]));
    }
}

template <class T>
void SpreadsAndCenters<T>::computeLocalSpreads2(vector<float>& spreads)
{
    spreads.clear();

    const int nids = localRows_.size();
    for (int igid = 0; igid < nids; igid++)
    {
        spreads.push_back(computeSpread2(localRows_[igid]));
    }
}

// get centers for functions centered in subdomain
template <class T>
void SpreadsAndCenters<T>::getLocalCenters(vector<Vector3D>& centers)
{
    centers.clear();

    const int nids = localRows_.size();
    for (int igid = 0; igid < nids; igid++)
    {
        centers.push_back(computeCenter(localRows_[igid]));
    }
}

template <class T>
void SpreadsAndCenters<T>::getLocalGids(vector<int>& lindex)
{
    lindex.clear();

    vector<int>::const_iterator it = localRowGid_.begin();
    while (it != localRowGid_.end())
    {
        lindex.push_back(*it);
        it++;
    }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
double SpreadsAndCenters<T>::computeDistance(const int st1, const int st2) const
{
    assert(st1 < ngids_);
    assert(st2 < ngids_);

    Vector3D center1(computeCenter(st1));
    Vector3D center2(computeCenter(st2));

    const short bc[3] = { 1, 1, 1 };
    double d          = center1.minimage(center2, cell_, bc);

    MGmol_MPI& mmpi = *(MGmol_MPI::instance());
    mmpi.bcast(&d, 1);

    return d;
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
void SpreadsAndCenters<T>::setSinCosData(VariableSizeMatrix<sparserow>& mat,
    const vector<int>& gids, const vector<int>& localRowGid)
{
    setData(mat, gids, localRowGid, r_);
}

///////////////////////////////////////////////////////////////////////////////
template <class T>
void SpreadsAndCenters<T>::setData(VariableSizeMatrix<sparserow>& mat,
    const vector<int>& gids, const vector<int>& localRowGid,
    vector<vector<double>>& matr)
{
    gids_ = gids;

    matr.clear();
    localRows_.clear();
    localRowGid_.clear();

    ngids_         = gids.size();
    const int ndir = mat.nnzrow(0);
    matr.resize(ndir);
    for (int col = 0; col < ndir; col++)
    {
        for (int i = 0; i < ngids_; i++)
        {
            matr[col].push_back(mat.getRowEntry(i, col));
        }
    }
    // setup the localRows vector
    for (vector<int>::const_iterator it = localRowGid.begin();
         it != localRowGid.end(); ++it)
    {
        int* lrindex = (int*)mat.getTableValue(*it);
        localRows_.push_back(*lrindex);
        localRowGid_.push_back(*it);
    }
}

template <class T>
void SpreadsAndCenters<T>::setSinCosData(vector<vector<double>>& a, const int n)
{
    r_.clear();
    localRows_.clear();

    ngids_         = n;
    const int ndir = a.size();
    r_.resize(ndir);
    for (int col = 0; col < ndir; col++)
    {
        for (int i = 0; i < ngids_; i++)
        {
            r_[col].push_back(a[col][i]);
        }
    }
    // setup the localRows vector
    for (int i = 0; i < ngids_; i++)
        localRows_.push_back(i);

    localRowGid_ = localRows_;
}

template <class T>
void SpreadsAndCenters<T>::computePositionMatrix(T& orbitals, T& work_orbitals)
{
    // copy into work array before normalizing
    work_orbitals.copyDataFrom(orbitals);

    // apply extra mask...
    //...

    work_orbitals.normalize();

    computeSinCosDiag(work_orbitals, true);
}

template <class T>
void SpreadsAndCenters<T>::computePositionMatrix(const T& orbitals)
{
    computeSinCosDiag(orbitals, false);
}

template <class T>
void SpreadsAndCenters<T>::computeSinCos(const T& orbitals)
{
    vector<vector<double>> a;
    a.resize(6);

    int n2 = orbitals.numst() * orbitals.numst();
    for (int k = 0; k < 6; k++)
    {
        a[k].resize(n2, 0.);
    }
    SinCosOps<T>::compute(orbitals, a);
    setSinCosData(a, n2);
}

template <class T>
void SpreadsAndCenters<T>::computeSinCosSquare(const T& orbitals)
{
    vector<vector<double>> a;
    a.resize(6);

    int n2 = orbitals.numst() * orbitals.numst();
    for (short k = 0; k < 6; k++)
    {
        a[k].resize(n2, 0.);
    }
    SinCosOps<T>::computeSquare(orbitals, a);
    setSinCosData(a, n2);
}

template <class T>
void SpreadsAndCenters<T>::computeSinCosSquare1D(
    const T& orbitals, const int dir)
{
    vector<vector<double>> a;
    a.resize(2);

    int n2 = orbitals.numst() * orbitals.numst();
    for (short k = 0; k < 2; k++)
    {
        a[k].resize(n2, 0.);
    }
    SinCosOps<T>::computeSquare1D(orbitals, a, dir);
    setSinCosData(a, n2);
}

template <class T>
void SpreadsAndCenters<T>::computeSinCos2states(
    const T& orbitals, const int st1, const int st2)
{
    vector<vector<double>> a;
    a.resize(6);

    int n2 = 4;
    for (int k = 0; k < 6; k++)
    {
        a[k].resize(n2, 0.);
    }
    SinCosOps<T>::compute2states(orbitals, a, st1, st2);
    setSinCosData(a, n2);
}

template <class T>
void SpreadsAndCenters<T>::computeSinCosDiag2states(
    const T& orbitals, const int st1, const int st2)
{
    vector<vector<double>> a;
    a.resize(6);

    for (int k = 0; k < 6; k++)
    {
        a[k].resize(2, 0.);
    }

    SinCosOps<T>::computeDiag2states(orbitals, a, st1, st2);
    setSinCosData(a, 2);
}

template <class T>
void SpreadsAndCenters<T>::computeSinCos1D(const T& orbitals, const int dir)
{
    vector<vector<double>> a;
    a.resize(2);

    int n2 = orbitals.numst() * orbitals.numst();
    for (short k = 0; k < 2; k++)
    {
        a[k].resize(n2, 0.);
    }

    SinCosOps<T>::compute1D(orbitals, a, dir);
    setSinCosData(a, n2);
}

template <class T>
void SpreadsAndCenters<T>::computeSinCos(const T& orbitals1, const T& orbitals2)
{
    vector<vector<double>> a;
    a.resize(6);

    int n2 = orbitals1.numst() * orbitals1.numst();
    for (short k = 0; k < 6; k++)
    {
        a[k].resize(n2, 0.);
    }
    SinCosOps<T>::compute(orbitals1, orbitals2, a);
    setSinCosData(a, n2);
}

template <class T>
void SpreadsAndCenters<T>::computeSinCosDiag(
    const T& orbitals, const bool normalized_functions)
{
    const int initTabSize = 4096;
    VariableSizeMatrix<sparserow> mat("SinCos", initTabSize);

    SinCosOps<T>::computeDiag(orbitals, mat, normalized_functions);
    setSinCosData(
        mat, orbitals.getAllOverlappingGids(), orbitals.getLocalGids());
}

template class SpreadsAndCenters<LocGridOrbitals<ORBDTYPE>>;
template class SpreadsAndCenters<ExtendedGridOrbitals<ORBDTYPE>>;
