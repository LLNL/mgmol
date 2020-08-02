// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef GRIDFUNCVECTOR_H
#define GRIDFUNCVECTOR_H

#include "GridFunc.h"
#include "GridFuncVectorInterface.h"

#include <map>
#include <vector>

namespace pb
{
template <typename ScalarType>
class GridFuncVector : public GridFuncVectorInterface
{
    // global id for functions in each subdivision
    const std::vector<std::vector<int>>& gid_;

    // gids to communicate to neighbors in each of 6 directions
    std::vector<std::vector<int>> remote_gids_;

    const Grid& grid_;
    const MPI_Comm comm_;

    std::vector<GridFunc<ScalarType>*> functions_;
    std::map<int, short> gid2lid_;

    // number of functions in functions_ (functions_.size())
    int nfunc_;

    int nfunc_max_global_;
    short nsubdivx_;

    bool up_;
    bool down_;
    bool east_;
    bool west_;
    bool north_;
    bool south_;

    bool updated_boundaries_;
    int bc_[3];
    short nghosts_;

    int dimx_;
    int dimy_;
    int dimz_;

    int incx_;
    int incy_;
    int dimxy_;
    int incxy_;
    size_t east_west_size_;
    size_t south_north_size_;
    size_t up_down_size_;

    const bool skinny_stencil_;

    static std::vector<std::vector<ScalarType>> comm_buf1_;
    static std::vector<std::vector<ScalarType>> comm_buf2_;
    static std::vector<std::vector<ScalarType>> comm_buf3_;
    static std::vector<std::vector<ScalarType>> comm_buf4_;

    MPI_Request req_east_west_[4];
    MPI_Request req_north_south_[4];
    MPI_Request req_up_down_[4];

    void wait_north_south();
    void wait_east_west();
    void wait_up_down();

    void allocate_buffers(const int nfunc);

    void initiateNorthSouthComm(const int begin_color, const int end_color);
    void finishNorthSouthComm();
    void initiateUpDownComm(const int begin_color, const int end_color);
    void finishUpDownComm();
    void initiateEastWestComm(const int begin_color, const int end_color);
    void finishEastWestComm();

    void communicateRemoteGids(const int begin_color, const int end_color);

    void allocate(const int n)
    {
        functions_.resize(n);
        for (int i = 0; i < n; i++)
        {
            functions_[i]
                = new GridFunc<ScalarType>(grid_, bc_[0], bc_[1], bc_[2]);
        }
    }

public:
    GridFuncVector(const Grid& my_grid, const int px, const int py,
        const int pz, const std::vector<std::vector<int>>& gid,
        const bool skinny_stencil = false)
        : gid_(gid),
          grid_(my_grid),
          comm_(my_grid.mype_env().comm()),
          skinny_stencil_(skinny_stencil)
    {
        bc_[0] = px;
        bc_[1] = py;
        bc_[2] = pz;

        updated_boundaries_ = false; // boundaries not initialized

        allocate(gid[0].size());

        setup();
    }

    ~GridFuncVector() override
    {
        assert(static_cast<int>(functions_.size()) == nfunc_);
        for (int i = 0; i < nfunc_; i++)
        {
            assert(functions_[i] != 0);
            delete functions_[i];
        }
    }

    void setup();

    template <typename ScalarType2>
    void assign(const int i, const ScalarType2* const v, const char dis = 'd')
    {
        assert(i < static_cast<int>(functions_.size()));

        functions_[i]->assign(v, dis);
        updated_boundaries_ = false;
    }
    GridFunc<ScalarType>& func(const int k) { return *functions_[k]; }
    const GridFunc<ScalarType>& func(const int k) const
    {
        return *functions_[k];
    }
    void trade_boundaries();
    void trade_boundaries_colors(const short, const short);

    size_t size() const { return functions_.size(); }

    // pointwise products this=A*B for each vector in this
    void prod(GridFuncVector<ScalarType>& A, const GridFunc<double>& B);
    void prod(GridFuncVector<ScalarType>& A, const GridFunc<float>& B);

    void extend3D(GridFuncVector<ScalarType>&);
    void restrict3D(GridFuncVector<ScalarType>&);
    void resetData()
    {
        assert(nfunc_ == static_cast<int>(functions_.size()));

        for (short k = 0; k < nfunc_; k++)
            functions_[k]->resetData();
        updated_boundaries_ = true;
    }
    void set_updated_boundaries(const bool flag) { updated_boundaries_ = flag; }
    GridFuncVector<ScalarType>& operator-=(
        const GridFuncVector<ScalarType>& func);
    void axpy(const double alpha, const GridFuncVector<ScalarType>& func);

    void init_vect(const int k, ScalarType* vv, const char dis) const;
    void getValues(const int k, double* vv) const;
    void getValues(const int k, float* vv) const;
};

} // namespace pb

#endif
