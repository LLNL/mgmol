// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_BLACSCONTEXT_H
#define MGMOL_BLACSCONTEXT_H

#include <mpi.h>

namespace dist_matrix
{

class BlacsContext
{
private:
    int ictxt_;
    int myrow_;
    int mycol_;
    int nprow_;
    int npcol_;
    int nprocs_;
    int size_;
    int myproc_;
    int mype_;
    bool onpe0_;
    bool active_;

    // MPI communicator
    MPI_Comm comm_active_; // corresponding to rectangular context
    const MPI_Comm comm_global_;

    void buildCommunicator();

    // keep assignment and copy constructors private
    BlacsContext& operator=(const BlacsContext&);
    BlacsContext(const BlacsContext&);

    void setup();

public:
    int ictxt() const { return ictxt_; }
    int myrow() const { return myrow_; }
    int mycol() const { return mycol_; }
    int nprow() const { return nprow_; }
    int npcol() const { return npcol_; }

    // number of processes in the context
    // int size() const { return size_; }
    int nprocs() const { return nprocs_; }
    // position of current process in row-major order
    // returns -1 if current process is not part of this context
    int myproc() const { return myproc_; }
    int mype() const { return mype_; }

    bool onpe0(void) const { return onpe0_; }
    bool active(void) const { return active_; }

    bool operator==(const BlacsContext& bc) const { return (this == &bc); }

    // MPI communicator for this context. Returns MPI_COMM_NULL if
    // this process is not part of the context
    MPI_Comm comm_active(void) const { return comm_active_; }
    MPI_Comm comm_global(void) const { return comm_global_; }

    // Constructors

    // default global context: construct a single-row global BlacsContext
    explicit BlacsContext(MPI_Comm comm);

    // specialized global BlacsContext
    // BlacsContext('r'): single-row global context (same as BlacsContext())
    // BlacsContext('c'): single-column global context
    // BlacsContext('s'): largest possible square context
    explicit BlacsContext(
        MPI_Comm comm, const char type, const int max_cpus = 10000);

    // construct a BlacsContext of size nprow * npcol
    explicit BlacsContext(MPI_Comm comm, const int nprow, const int npcol);

    // construct a BlacsContext of size nprow*npcol starting at process ipe
    explicit BlacsContext(
        MPI_Comm, const int ipe, const int nprow, const int npcol);

    // construct a BlacsContext of size nprow*npcol from the processes
    // in context bc lying in the rectangle of size nr * nc starting
    // at process (irow,icol)
    BlacsContext(BlacsContext& bc, const int irow, const int icol, const int nr,
        const int nc, const int nprow, const int npcol);

    // construct a BlacsContext corresponding to "my row" or "my column"
    // of a given context bc
    // use: BlacsContext(bc,'r') or BlacsContext(bc,'c')
    BlacsContext(const BlacsContext& bc, const char type);

    ~BlacsContext();
};

} // namespace
#endif
