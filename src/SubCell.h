// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_SUBCELL_H
#define MGMOL_SUBCELL_H

#include "Grid.h"
#include "Vector3D.h"

class SubCell
{
    double inner_radius_;
    double outer_radius_;
    std::vector<Vector3D> centers_;
    const Vector3D cell_dimensions_;
    Vector3D subcell_dimensions_;
    short subdivx_;
    /*
     * lower left corner for each subdivision
     */
    std::vector<std::vector<double>> ll_;
    std::vector<std::vector<double>> ur_;

public:
    SubCell(const pb::Grid& fine_grid, const short subdivx, const short level);

    double outerRadius() const { return outer_radius_; }

    // does a sphere of radius "radius" centered at "center" overlap with the
    // subdomain?
    bool spherePossibleOverlap(const Vector3D& center, const double radius,
        const short iloc, const short bcPoisson[3]) const;

    double distance(const Vector3D& center, const short iloc,
        const short bcPoisson[3]) const
    {
        return centers_[iloc].minimage(center, cell_dimensions_, bcPoisson);
    }

    bool includes(const Vector3D& point, const short iloc) const
    {
        const double tol = 1.e-8;
        bool result      = true;
        for (short i = 0; i < 3; i++)
        {
            result = result && (point[i] >= ll_[iloc][i] - tol);
            result = result && (point[i] <= ur_[iloc][i] + tol);
        }
        return result;
    }

    void printCenter(const short iloc, std::ostream& os)
    {
        os << centers_[iloc];
        os << std::endl;
    }

    void print(const short iloc, std::ostream& os)
    {
        os << "Subdomain (" << ll_[iloc][0] << "," << ll_[iloc][1] << ","
           << ll_[iloc][2] << ")-(" << ur_[iloc][0] << "," << ur_[iloc][1]
           << "," << ur_[iloc][2] << ")" << std::endl;
    }
};

#endif
