// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
 * Variable size csr/csc matrix used for data transfer operations
 */
#ifndef MGMOL_VARIABLESIZEMATRIXINTERFACE_H_
#define MGMOL_VARIABLESIZEMATRIXINTERFACE_H_

// typedef enum INSERTMODE  {INSERT, ADD} INSERTMODE;

/* define maximum and minimum local matrix size */
#define MAX_MAT_SIZE 10000
#define MIN_MAT_SIZE 10
/* define default tolerance for pruning matrix entries */
#define MAT_TOL 1.0e-14
/* define maximum number of print rows */
#define MAX_PRINT_ROWS 100
/* define default number of print rows for diagnostics */
#define NUM_PRINT_ROWS 5

class VariableSizeMatrixInterface
{
protected:
    static Timer insert_tm_;
    static Timer updateRow_tm_;
    static Timer insertRow_tm_;
    static Timer sort_col_tm_;
    static Timer AmultSymBdiag_tm_;
    static Timer AmultSymB_ij_tm_;
    static Timer AmultSymBLocal_tm_;
    static Timer AmultSymB_tm_;

public:
    virtual ~VariableSizeMatrixInterface() { }

    static void printTimers(std::ostream& os)
    {
        AmultSymBdiag_tm_.print(os);
        AmultSymB_ij_tm_.print(os);
        AmultSymBLocal_tm_.print(os);
        AmultSymB_tm_.print(os);
        insert_tm_.print(os);
        updateRow_tm_.print(os);
        insertRow_tm_.print(os);
        sort_col_tm_.print(os);
    }
};

#endif
