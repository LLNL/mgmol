//
// test VariableSizeMatrix
//
// Perform a series of operations to test variableSizeMatrix class.
// Apart from print operations, the operations performed by this class
// are all local operations, hence this unit test is designed to use
// only one processor (without loss of generality).
//
#include "MGmol_MPI.h"
#include "VariableSizeMatrix.h"

#include <boost/test/unit_test.hpp>

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace tt = boost::test_tools;
namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(variable_size_matrix_insert, *utf::tolerance(1.e-12))
{
    const int lsize          = 15;

    // start with a matrix with only half of the rows we need
    VariableSizeMatrix<sparserow> mat("A", lsize / 2);

#ifdef _OPENMP
    int numthreads = omp_get_max_threads();
#else
    int numthreads = 1;
#endif
    std::cout << "Number of threads = " << numthreads << std::endl;

    // repeat insertions 10 times
    for (int n = 0; n < 10; n++)
    {

#pragma omp parallel for
        for (int j = 0; j < lsize; j++)
        {
            for (int i = 0; i < lsize; i++)
                mat.insertMatrixElement(i, j, 1., ADD, true);
        }
    }

    // verify values in matrix
    for (int i = 0; i < lsize; i++)
    {
        for (int j = 0; j < lsize; j++)
        {
            double val = mat.get_value(i, j);
            BOOST_TEST(val == 10.);
        }
    }

}

BOOST_AUTO_TEST_CASE(variable_size_matrix)
{
    int mype;
    int npes;
    const int lsize          = 15;
    const int num_print_rows = 5;
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    BOOST_TEST(
        npes == 1, "ERROR: This example is set up to use only 1 process");

    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    VariableSizeMatrix<sparserow> matA("matA", lsize);
    VariableSizeMatrix<sparserow> matB("matB", lsize);
    VariableSizeMatrix<sparserow> matC("matC", lsize);

    /* generate matrix data */
    const double vals[] = { 1., 2., 1. };
    std::vector<int> cols(3, -1);
    std::vector<int> rows;
    for (int i = 0; i < lsize; i++)
    {
        rows.push_back(i);
    }

    matA.setupSparseRows(rows);
    // row 0:
    cols[0] = 0;
    cols[1] = 1;
    matA.initializeLocalRow(2, rows[0], cols.data(), &vals[1]);
    // row n-1:
    cols[0] = lsize - 2;
    cols[1] = lsize - 1;
    matA.initializeLocalRow(2, rows[lsize - 1], cols.data(), &vals[0]);
    // rows 1 to (n-2)
    for (int i = 1; i < lsize - 1; i++)
    {
        cols[0] = i - 1;
        cols[1] = i;
        cols[2] = i + 1;
        matA.initializeLocalRow(3, rows[i], cols.data(), vals);
    }

    // print matrix
    std::ostream& sout = std::cout;
    sout << " " << num_print_rows << " x " << num_print_rows
         << " principal submatrix of matrix: " << matA.name() << std::endl;
    matA.print(sout, rows, num_print_rows);
    // print some matrix stats
    std::cout << "Check sizes ..." << std::endl;
    std::cout << "matA.n() = " << matA.n() << std::endl;
    std::cout << "matA.nnzmat() = " << matA.nnzmat() << std::endl;
    std::cout << "matA: max row nnz = " << matA.nzmax() << std::endl;
    std::cout << "matA: min row nnz = " << matA.nzmin() << std::endl;
    // Insert new row
    std::cout << "Insert new row ... " << std::endl;
    cols[0] = lsize - 1;
    cols[1] = lsize;
    matA.insertNewRow(2, lsize, cols.data(), &vals[0], true);
    // Insert entry
    std::cout << "Insert entry to existing ... " << std::endl;
    matA.insertMatrixElement((lsize - 1), lsize, 1, INSERT, true);
    // print some matrix stats
    std::cout << "Check sizes ..." << std::endl;
    std::cout << "matA.n() = " << matA.n() << std::endl;
    std::cout << "matA.nnzmat() = " << matA.nnzmat() << std::endl;
    std::cout << "matA: max row nnz = " << matA.nzmax() << std::endl;
    std::cout << "matA: min row nnz = " << matA.nzmin() << std::endl;
    // Initialize matrix matB
    std::cout << "Initialize matB: matB = matA ..." << std::endl;
    matB = matA;
    // Axpy
    std::cout << "Check axpy: matB.axpy(-1., matA) ... " << std::endl;
    matB.axpy(-1., matA);
    double sum = 0.;
    for (int i = 0; i <= lsize; i++)
    {
        sum += matB.pnorm(i, 0);
    }
    BOOST_TEST(sum == 0.0, tt::tolerance(1.e-8));

    // set2Identity
    std::cout << "set to identity ..." << std::endl;
    matB.set2Identity();
    // gemv
    std::cout << "Check gemv ... " << std::endl;
    std::vector<double> x((lsize + 1), 1.0);
    std::vector<double> y((lsize + 1), 0.);
    matB.gemv(1.0, x, 0., y);
    sum = std::accumulate(y.begin(), y.end(), 0.);
    BOOST_TEST(sum == (double)(lsize + 1), tt::tolerance(1.e-8));

    // scale
    double scal = 2.0;
    std::cout << "Check scale ... " << std::endl;
    matB.scale(scal);
    // matrix multiplication
    std::cout << "Check matrix multiplication ... " << std::endl;
    rows.push_back(lsize);
    VariableSizeMatrix<sparserowtab> pattern(matA, true);
    matC.setupSparseRows(rows);
    matA.AmultSymB(&matB, matC, pattern);
    double sumc = 0.;
    sum         = 0.;
    for (int i = 0; i <= lsize; i++)
    {
        sum += matA.pnorm(i, 0);
        sumc += matC.pnorm(i, 0);
    }
    sum *= scal;
    BOOST_TEST(sum == sumc, tt::tolerance(1.e-8));
    // Check trace
    std::cout << "Check trace ..." << std::endl;
    BOOST_TEST(matC.trace() == (scal * matA.trace()), tt::tolerance(1.e-8));
}
