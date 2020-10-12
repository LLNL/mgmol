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
#include "memory_space.h"

#include <map>
#include <memory>
#include <vector>

namespace pb
{
#ifdef HAVE_OPENMP_OFFLOAD
template <typename ScalarType, typename MemorySpaceType = MemorySpace::Device>
#else
template <typename ScalarType, typename MemorySpaceType = MemorySpace::Host>
#endif
class GridFuncVector
{
    static Timer trade_bc_tm_;
    static Timer trade_bc_colors_tm_;
    static Timer prod_tm_;
    static Timer finishExchangeNorthSouth_tm_;
    static Timer finishExchangeUpDown_tm_;
    static Timer finishExchangeEastWest_tm_;
    static Timer wait_north_south_tm_;
    static Timer wait_up_down_tm_;
    static Timer wait_east_west_tm_;

    // block of memory for all GridFunc
    std::unique_ptr<ScalarType> memory_;

#ifdef HAVE_OPENMP_OFFLOAD
    // block of memory for device memory
    std::unique_ptr<ScalarType, void (*)(ScalarType*)> functions_dev_;
#endif

    using MemoryST = MemorySpace::Memory<ScalarType, MemorySpaceType>;

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

    // number of functions associated with buffers sizes
    int nfunc4buffers_;

    std::vector<std::vector<ScalarType>> comm_buf1_;
    std::vector<std::vector<ScalarType>> comm_buf2_;
    std::vector<std::vector<ScalarType>> comm_buf3_;
    std::vector<std::vector<ScalarType>> comm_buf4_;

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

    void allocate(const int n);

public:
    GridFuncVector(const Grid& my_grid, const int px, const int py,
        const int pz, const std::vector<std::vector<int>>& gid,
        const bool skinny_stencil = false)
        : gid_(gid),
          grid_(my_grid),
          comm_(my_grid.mype_env().comm()),
          skinny_stencil_(skinny_stencil),
          nfunc4buffers_(0)
#ifdef HAVE_OPENMP_OFFLOAD
          ,
          functions_dev_(MemoryST::allocate(gid[0].size() * my_grid.sizeg()),
              MemoryST::free)
#endif
    {
        bc_[0] = px;
        bc_[1] = py;
        bc_[2] = pz;

        updated_boundaries_ = false; // boundaries not initialized

        allocate(gid[0].size());

        setup();
    }

    ~GridFuncVector()
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
        assert(functions_[i] != nullptr);

        functions_[i]->assign(v, dis);
        updated_boundaries_ = false;
    }

#ifdef HAVE_OPENMP_OFFLOAD
    void copyHtoD(int size)
    {
        MemorySpace::copy_to_dev(memory_.get(), size, functions_dev_.get());
    }

    void copyDtoH(int size)
    {
        MemorySpace::copy_to_host(functions_dev_.get(), size, memory_.get());
    }
#endif

    GridFunc<ScalarType>& getGridFunc(const int k)
    {
        assert(k < static_cast<int>(functions_.size()));
        assert(functions_[k] != nullptr);

        return *functions_[k];
    }

    const GridFunc<ScalarType>& getGridFunc(const int k) const
    {
        assert(k < static_cast<int>(functions_.size()));
        assert(functions_[k] != nullptr);

        return *functions_[k];
    }
    void trade_boundaries();
    void trade_boundaries_colors(const short, const short);

    size_t size() const { return functions_.size(); }

    // pointwise products this=A*B for each vector in this
    void prod(GridFuncVector<ScalarType, MemorySpaceType>& A,
        const GridFunc<double>& B);
    void prod(GridFuncVector<ScalarType, MemorySpaceType>& A,
        const GridFunc<float>& B);

    void extend3D(GridFuncVector<ScalarType, MemorySpaceType>&);
    void restrict3D(GridFuncVector<ScalarType, MemorySpaceType>&);
    void resetData()
    {
        assert(nfunc_ == static_cast<int>(functions_.size()));

        for (short k = 0; k < nfunc_; k++)
            functions_[k]->resetData();
        updated_boundaries_ = true;
    }
    void set_updated_boundaries(const bool flag) { updated_boundaries_ = flag; }
    GridFuncVector<ScalarType, MemorySpaceType>& operator-=(
        const GridFuncVector<ScalarType, MemorySpaceType>& func);
    void axpy(const double alpha,
        const GridFuncVector<ScalarType, MemorySpaceType>& func);

    void init_vect(const int k, ScalarType* vv, const char dis) const;

    template <typename InputScalarType>
    void getValues(const int k, InputScalarType* vv) const;
    // template <typename MemorySpaceType>
    // void getValues(const int k, float* vv) const;

    static void printTimers(std::ostream& os)
    {
        trade_bc_tm_.print(os);
        trade_bc_colors_tm_.print(os);
        prod_tm_.print(os);
        wait_north_south_tm_.print(os);
        wait_up_down_tm_.print(os);
        wait_east_west_tm_.print(os);
        finishExchangeNorthSouth_tm_.print(os);
        finishExchangeUpDown_tm_.print(os);
        finishExchangeEastWest_tm_.print(os);
    }
};

template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::trade_bc_tm_(
    "GridFuncVector::trade_bc");
template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::trade_bc_colors_tm_(
    "GridFuncVector::trade_bc_colors");
template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::prod_tm_(
    "GridFuncVector::prod");
template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::finishExchangeNorthSouth_tm_(
    "GridFuncVector::finishExNorthSouth");
template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::finishExchangeUpDown_tm_(
    "GridFuncVector::finishExUpDown");
template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::finishExchangeEastWest_tm_(
    "GridFuncVector::finishExEastWest");
template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::wait_north_south_tm_(
    "GridFuncVector::waitNS");
template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::wait_up_down_tm_(
    "GridFuncVector::waitUD");
template <typename ScalarType, typename MemorySpaceType>
Timer GridFuncVector<ScalarType, MemorySpaceType>::wait_east_west_tm_(
    "GridFuncVector::waitEW");

} // namespace pb

#endif
