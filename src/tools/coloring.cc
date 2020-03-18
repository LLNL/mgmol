// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#include <cassert>
#include <list>
#include <math.h>
#include <string.h>

#include "coloring.h"

/* C version of greedy multicoloring algorithm adapted from SPARSKIT (LGPL) */
void greedyMC(int n, int* ja, int* ia, int* num_colors, int* iord, int* colors)
{
    int i, j, k, maxcol, mycol, ncols;

    std::memset(colors, -1, n * sizeof(int));

    maxcol   = n;
    ncols    = 0;
    int* tmp = new int[n];
    std::memset(tmp, 0, n * sizeof(int));

    /* scan all nodes */
    for (i = 0; i < maxcol; i++)
    {
        /* look at adjacent nodes for previously assigned colors */
        int mcol = -1;
        int ii   = iord[i];
        int k1   = ia[ii];
        int k2   = ia[ii + 1];
        for (k = k1; k < k2; k++)
        {
            j        = ja[k];
            int icol = colors[j];
            if (icol != -1)
            {
                mcol = mcol > icol ? mcol : icol;
                /* mark neighbor's color */
                tmp[icol] = -1;
            }
        }

        /* scan marked colors to assign color */
        for (mycol = 0; mycol <= mcol; mycol++)
            if (tmp[mycol] != -1) break;
        /* reset tmp array */
        std::memset(tmp, 0, (mcol + 1) * sizeof(int));

        /* assign color and update num_colors */
        colors[ii] = mycol;
        ncols      = ncols > mycol ? ncols : mycol;
    }
    /* done with coloring. Assign values */
    *num_colors = ncols + 1;

    /* check result */
#ifndef NDEBUG
    for (i = 0; i < n; i++)
    {
        int k1 = ia[i];
        int k2 = ia[i + 1];
        for (k = k1; k < k2; k++)
        {
            j = ja[k];
            if (i == j) continue;
            assert(colors[i] != colors[j]);
        }
    }
#endif

    delete[] tmp;
    return;
}

void colorRLF(const SymmetricMatrix<char>& overlaps,
    std::list<std::list<int>>& colored_gids, const bool verbose,
    std::ostream& os)
{
    const int dim = overlaps.dimension();
    if (verbose)
    {
        os << "Uses Recursive Largest First (RLF) algorithm"
           << " for problem of size " << dim << std::endl;
    }

    colored_gids.clear();

    // initialize list of "uncolored" gids
    std::list<int> uncolored;
    const std::vector<int>& gids(overlaps.gids());
    for (int i = 0; i < dim; i++)
        uncolored.push_back(gids[i]);

    while (uncolored.size() > 0)
    {

        // compute degrees of vertices
        std::multimap<int, int> degrees;
        std::list<int>::const_iterator pj = uncolored.begin();
        while (pj != uncolored.end())
        {
            const int j0                      = (*pj);
            int degree                        = 0;
            std::list<int>::const_iterator pi = uncolored.begin();
            while (pi != uncolored.end())
            {
                degree += overlaps.val(*pi, j0);
                pi++;
            }
            if (degree) degrees.insert(std::pair<int, int>(degree, j0));
            pj++;
        }

        // find the 1st color: largest degree
        std::multimap<int, int>::const_iterator p0 = degrees.end();
        p0--;
        const int v0 = p0->second;
        uncolored.remove(v0);

        // start list of functions of same color
        std::list<int> vi;
        vi.push_back(v0);

        // compute set u of vertices adjacents to v0
        // compute set v of vertices non-adjacents to v0
        std::list<int> u;
        std::list<int> v;
        pj = uncolored.begin();
        while (pj != uncolored.end())
        {
            int i0 = (*pj);
            if (overlaps.val(i0, v0))
            {
                u.push_back(i0);
            }
            else
            {
                v.push_back(i0);
            }
            pj++;
        }

        // try to color another vertex with color v0
        while (v.size() > 0)
        {

            // compute degrees in u for functions in v
            std::multimap<int, int> udegrees;
            std::list<int>::const_iterator pv = v.begin();
            while (pv != v.end())
            {
                const int j0                      = (*pv);
                int degree                        = 0;
                std::list<int>::const_iterator pu = u.begin();
                while (pu != u.end())
                {
                    int i0 = (*pu);
                    degree += overlaps.val(i0, j0);
                    pu++;
                }
                udegrees.insert(std::pair<int, int>(degree, j0));
                pv++;
            }

            // find the largest degree in udegrees
            assert(udegrees.size() > 0);
            std::multimap<int, int>::const_iterator pu0 = udegrees.end();
            pu0--;
            const int u0 = pu0->second;

            // use same color for v0 and u0
            vi.push_back(u0);
            uncolored.remove(u0);

#ifndef NDEBUG
            if (verbose)
            {
                os << "PE 0: Same color for index " << *vi.begin() << " and "
                   << u0 << std::endl;
                assert(!overlaps.val(*vi.begin(), u0));
            }
#endif

            // add to list u adjacents to u0
            std::list<int> vremove;
            pj = v.begin();
            while (pj != v.end())
            {
                const int i0 = (*pj);
                if (overlaps.val(i0, u0))
                {
                    u.push_back(i0);
                    vremove.push_back(i0);
                }
                pj++;
            }

            // remove from list v adjacents to u0
            pj = vremove.begin();
            while (pj != vremove.end())
            {
                v.remove(*pj);
                pj++;
            }

        } // done with a color

        colored_gids.push_back(vi);
    }
}

/*
   quicksort of the elements in a. elements in b are permuted
   according to the resulting sorted a.
*/
void quicksortI_h2l(int* a, int* b, int lo, int hi)
{
    //  lo is the lower index, hi is the upper index
    //  of the region of array a that is to be sorted
    int i = lo, j = hi;
    int v;
    int x = ceil(a[(lo + hi) / 2]);
    int q;

    //  partition
    do
    {
        while (a[i] > x)
            i++;
        while (a[j] < x)
            j--;
        if (i <= j)
        {
            v    = a[i];
            a[i] = a[j];
            a[j] = v;
            q    = b[i];
            b[i] = b[j];
            b[j] = q;
            i++;
            j--;
        }
    } while (i <= j);

    //  recursion
    if (lo < j) quicksortI_h2l(a, b, lo, j);
    if (i < hi) quicksortI_h2l(a, b, i, hi);
}

void greedyColor(const SymmetricMatrix<char>& overlaps,
    std::list<std::list<int>>& colored_gids, const bool verbose,
    std::ostream& os)
{
    if (verbose)
    {
        os << "Uses greedy algorithm" << std::endl;
    }

    const int dim = overlaps.dimension();
    std::vector<int> vecia;
    std::vector<int> vecja;
    vecia.reserve(dim + 1);
    vecja.reserve(dim);

    int* nnzrow = new int[dim];
    int* iord   = new int[dim];

    /* nonzero pattern of overlap matrix */
    overlaps.getNNZPattern(vecia, vecja, nnzrow);

    for (int i = 0; i < dim; i++)
        iord[i] = i;

    /* sort degree of nodes from hi to lo */
    quicksortI_h2l(nnzrow, iord, 0, dim - 1);
    /* perform greedy multicoloring */
    int num_colors                = 0;
    std::vector<int>::iterator ia = vecia.begin();
    std::vector<int>::iterator ja = vecja.begin();
    /* overwrite nnzrow array -- no longer needed */
    int* colors = nnzrow;
    greedyMC(dim, &(*ja), &(*ia), &num_colors, iord, colors);

    /* populate list of colored_grids */
    for (int color = 0; color < num_colors; color++)
    {
        std::list<int> llist;
        for (int i = 0; i < dim; i++)
        {
            if (colors[i] == color) llist.push_back(i);
        }
        colored_gids.push_back(llist);
    }

    delete[] iord;
    delete[] colors;

    return;
}
