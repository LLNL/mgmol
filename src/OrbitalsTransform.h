// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef OrbitalsTRANSFORM_H
#define OrbitalsTRANSFORM_H

#include "BlacsContext.h"
#include "DistMatrix.h"
#include "Vector3D.h"

#include <string.h>
#include <vector>

#define NDIM 3

class OrbitalsTransform
{
protected:
    // row dist_matrix::BlacsContext
    dist_matrix::BlacsContext* bcr_;

    // MPI communicator
    MPI_Comm comm_;

    // max. number columns on each processor
    int bsize_;

    // total number of valid states
    int nst_;

    // number of valid states/columns treated on processor
    int lnst_;

    // cell dimensions
    Vector3D cell_;

    // cell origin
    Vector3D origin_;

    // index first column on processor
    int offset_;

    // transformation matrix
    std::vector<DISTMATDTYPE> mat_;

    dist_matrix::DistMatrix<DISTMATDTYPE>* a_; // transformation

    std::vector<int> iu_; // allowed transformation

    std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>*>
        r_; // cosine and sine matrices for x,y,z directions

public:
    virtual ~OrbitalsTransform() = 0;

    void printTransformationMatrix() const;

    std::vector<dist_matrix::DistMatrix<DISTMATDTYPE>*>& r(void) { return r_; }

    dist_matrix::DistMatrix<DISTMATDTYPE>& r(const int i)
    {
        assert(i < 6 && i >= 0);
        assert(r_[i] != NULL);
        return *r_[i];
    }

    std::vector<DISTMATDTYPE>& mat(void) { return mat_; }

    dist_matrix::DistMatrix<DISTMATDTYPE>& a(void) { return *a_; }

    // compute Orbitals transform
    virtual void compute_transform(const int maxsweep, const double tol) = 0;

    int nst() const { return nst_; }

    virtual double spread2(int i, int j) const = 0;
    double spread2(int i) const;
    double spread2(void) const;
    double spread(int i) const;
    double spread(void) const;
    Vector3D center(int i) const;
    Vector3D center(void) const;

    // total volume of spheres of radius "spread"
    double volume() const;

    void printCentersAndSpreads(std::ostream& os);

    OrbitalsTransform(
        const int nst, const Vector3D& origin, const Vector3D& ll);

    void distributeColumnsR(std::vector<std::vector<DISTMATDTYPE>>& vmm);
};
#endif
