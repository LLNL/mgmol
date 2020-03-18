// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE
#ifdef HAVE_TRICUBIC

#include "TriCubic.h"
#include "Delh4.h"
#include "Grid.h"
#include "tricubic.h"

namespace pb
{

template <class T>
TriCubic<T>::TriCubic(const Grid& grid, const short bc[3]) : grid_(grid)
{
    for (short i = 0; i < 3; i++)
    {
        h_[i]  = grid_.hgrid(i);
        bc_[i] = bc[i];
    }
}
template <class T>
void TriCubic<T>::computeSplineCoeffs(const T* const f)
{
    Grid tmp_grid(grid_, 2);
    FDoper<T>* myoper_del[3];
    myoper_del[0] = new Delxh4<T>(tmp_grid);
    myoper_del[1] = new Delyh4<T>(tmp_grid);
    myoper_del[2] = new Delzh4<T>(tmp_grid);

    GridFunc<T> gf_f(f, tmp_grid, 1, 1, 1, 'd');

    // compute 1st derivatives
    GridFunc<T>* gf_df[3];
    GridFunc<T>* gf_d2f[3];
    for (short i = 0; i < 3; i++)
    {
        gf_df[i]  = new GridFunc<T>(tmp_grid, 1, 1, 1);
        gf_d2f[i] = new GridFunc<T>(tmp_grid, 1, 1, 1);
    }
    GridFunc<T> gf_d3f(tmp_grid, 1, 1, 1);
    ;

    for (short i = 0; i < 3; i++)
    {
        myoper_del[i]->apply(gf_f, *gf_df[i]);

        // compute 2nd derivatives
        myoper_del[(i + 1) % 3]->apply(*gf_df[i], *gf_d2f[i]);
    }

    // compute 3rd derivative
    myoper_del[2]->apply(*gf_d2f[0], gf_d3f);

    // fill up ghosts when not done by Delxh4
    gf_d2f[1]->trade_boundaries();
    gf_d2f[2]->trade_boundaries();
    gf_d3f.trade_boundaries();

    computeSplineCoeffs(gf_f, *gf_df[0], *gf_df[1], *gf_df[2], *gf_d2f[0],
        *gf_d2f[2], *gf_d2f[1], gf_d3f);

    for (short i = 0; i < 3; i++)
    {
        delete gf_df[i];
        delete gf_d2f[i];
        delete myoper_del[i];
    }
}

template <class T>
void TriCubic<T>::computeSplineCoeffs(const GridFunc<T>& f,
    const GridFunc<T>& fx, const GridFunc<T>& fy, const GridFunc<T>& fz,
    const GridFunc<T>& fxy, const GridFunc<T>& fxz, const GridFunc<T>& fyz,
    const GridFunc<T>& fxyz)
{
    assert(f.ghost_pt() >= 1);

    double fval[8];
    double dfdxval[8];
    double dfdyval[8];
    double dfdzval[8];
    double d2fdxdyval[8];
    double d2fdxdzval[8];
    double d2fdydzval[8];
    double d3fdxdydzval[8];

    const int dim[3] = { grid_.dim(0), grid_.dim(1), grid_.dim(2) };

    const int incx = dim[2] * dim[1];
    const int incy = dim[2];
    const int n3   = dim[2] * dim[1] * dim[0];

    spline_coeff_.resize(n3);
    for (int i = 0; i < n3; i++)
        spline_coeff_[i].resize(64);

    // loop over cells
    for (int i = 0; i < dim[0]; i++)
        for (int j = 0; j < dim[1]; j++)
            for (int k = 0; k < dim[2]; k++)
            {
                f.getCellCornersValues(i, j, k, fval);

                fx.getCellCornersValues(i, j, k, dfdxval);
                fy.getCellCornersValues(i, j, k, dfdyval);
                fz.getCellCornersValues(i, j, k, dfdzval);

                fxy.getCellCornersValues(i, j, k, d2fdxdyval);
                fxz.getCellCornersValues(i, j, k, d2fdxdzval);
                fyz.getCellCornersValues(i, j, k, d2fdydzval);

                fxyz.getCellCornersValues(i, j, k, d3fdxdydzval);

                for (short q = 0; q < 8; q++)
                {
                    dfdxval[q] *= h_[0];
                    dfdyval[q] *= h_[1];
                    dfdzval[q] *= h_[2];

                    d2fdxdyval[q] *= h_[0] * h_[1];
                    d2fdxdzval[q] *= h_[0] * h_[2];
                    d2fdydzval[q] *= h_[1] * h_[2];

                    d3fdxdydzval[q] *= h_[0] * h_[1] * h_[2];
                }

                int index = i * incx + j * incy + k;
                tricubic_get_coeff(&spline_coeff_[index][0], fval, dfdxval,
                    dfdyval, dfdzval, d2fdxdyval, d2fdxdzval, d2fdydzval,
                    d3fdxdydzval);
            }
}
template <class T>
void TriCubic<T>::getValues(
    const vector<double>& vr, vector<double>& val, MPI_Comm comm)
{
    assert(vr.size() / 3 == val.size());
    assert(vr.size() % 3 == 0);

    const int dim[3] = { grid_.dim(0), grid_.dim(1), grid_.dim(2) };

    const int incx = dim[2] * dim[1];
    const int incy = dim[2];

    const int n = (int)val.size();

    const double hinv0 = 1. / h_[0];
    const double hinv1 = 1. / h_[1];
    const double hinv2 = 1. / h_[2];

    // cout<<"TriCubic::getValues(), n="<<n<<endl;
    for (int ir = 0; ir < n; ir++)
    {
        double rincell[3] = { vr[3 * ir], vr[3 * ir + 1], vr[3 * ir + 2] };
        const double origin[3]
            = { grid_.origin(0), grid_.origin(1), grid_.origin(2) };
        const double lattice[3] = { grid_.ll(0), grid_.ll(1), grid_.ll(2) };
        for (short d = 0; d < 3; d++)
        {
            // periodic bc
            if (bc_[d] == 1)
            {
                while (rincell[d] < origin[d])
                    rincell[d] += lattice[d];
                while (rincell[d] >= origin[d] + lattice[d])
                    rincell[d] -= lattice[d];
            }
        }

        const double rlocal[3] = { rincell[0] - grid_.start(0),
            rincell[1] - grid_.start(1), rincell[2] - grid_.start(2) };

        const double sx = floor(rlocal[0] * hinv0);
        const double sy = floor(rlocal[1] * hinv1);
        const double sz = floor(rlocal[2] * hinv2);

        double x = (rlocal[0] * hinv0 - sx);
        double y = (rlocal[1] * hinv1 - sy);
        double z = (rlocal[2] * hinv2 - sz);

        assert(x >= 0.);
        assert(y >= 0.);
        assert(z >= 0.);
        assert(x <= 1.);
        assert(y <= 1.);
        assert(z <= 1.);

        if (sx >= 0. && sy >= 0. && sz >= 0. && (int)sx < dim[0]
            && (int)sy < dim[1] && (int)sz < dim[2])
        {
            const int index = (int)sx * incx + (int)sy * incy + (int)sz;
            assert(index < (int)spline_coeff_.size());
            vector<double>& sp_ref(spline_coeff_[index]);
            double* psp = &sp_ref[0];
            val[ir]     = tricubic_eval(psp, x, y, z);
        }
        else
        {
            val[ir] = 0.;
            cout << "WARNING: TriCubic::getValues(), position outside domain..."
                 << endl;
        }
    }
}

template <class T>
void TriCubic<T>::getGradient(const double r[3], double dfdr[3], MPI_Comm comm)
{
    const int dim[3] = { grid_.dim(0), grid_.dim(1), grid_.dim(2) };

    const int incx = dim[2] * dim[1];
    const int incy = dim[2];

    const double origin[3]
        = { grid_.origin(0), grid_.origin(1), grid_.origin(2) };
    const double lattice[3] = { grid_.ll(0), grid_.ll(1), grid_.ll(2) };

    double rincell[3] = { r[0], r[1], r[2] };
    for (short i = 0; i < 3; i++)
    {
        // periodic bc
        if (bc_[i] == 1)
        {
            while (rincell[i] < origin[i])
                rincell[i] += lattice[i];
            while (rincell[i] >= origin[i] + lattice[i])
                rincell[i] -= lattice[i];
        }
    }

    const double rlocal[3] = { rincell[0] - grid_.start(0),
        rincell[1] - grid_.start(1), rincell[2] - grid_.start(2) };

    const double sx = floor(rlocal[0] / h_[0]);
    const double sy = floor(rlocal[1] / h_[1]);
    const double sz = floor(rlocal[2] / h_[2]);

    double x = (rlocal[0] - h_[0] * sx) / h_[0];
    double y = (rlocal[1] - h_[1] * sy) / h_[1];
    double z = (rlocal[2] - h_[2] * sz) / h_[2];

    assert(x >= 0.);
    assert(y >= 0.);
    assert(z >= 0.);
    assert(x <= 1.);
    assert(y <= 1.);
    assert(z <= 1.);

    if (sx >= 0. && sy >= 0. && sz >= 0. && (int)sx < dim[0] && (int)sy < dim[1]
        && (int)sz < dim[2])
    {
        const int index = (int)sx * incx + (int)sy * incy + (int)sz;
        dfdr[0]
            = tricubic_eval(&spline_coeff_[index][0], x, y, z, 1, 0, 0) / h_[0];
        dfdr[1]
            = tricubic_eval(&spline_coeff_[index][0], x, y, z, 0, 1, 0) / h_[1];
        dfdr[2]
            = tricubic_eval(&spline_coeff_[index][0], x, y, z, 0, 0, 1) / h_[2];
    }
    else
    {
        for (short i = 0; i < 3; i++)
            dfdr[i] = 0.;
        cout << "WARNING: TriCubic::getGradient(), position outside domain..."
             << endl;
    }
}
template class TriCubic<double>;
template class TriCubic<float>;
}
#endif
