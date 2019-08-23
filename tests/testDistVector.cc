//
// test DistMatrix
//
// multiply a matrix a(m,k) by b(k,n) to get c(m,n)
// using blocks of size (mb,nb) on a process grid (nprow,npcol)
//
// use: testDistMatrix [-check] input_file
// input:
// nprow npcol
// m_a n_a mb_a nb_a transa
// m_b n_b mb_b nb_b transb
// m_c n_c mb_c nb_c
//
#include "BlacsContext.h"
#include "DistMatrix.h"
#include "DistVector.h"
#include "MGmol_MPI.h"

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv)
{
    int mype;
    int npes;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    if (npes != 4)
    {
        std::cout << "Number of processors in MPI: " << npes << std::endl;
        std::cout << "This example to set up to use only 4 processes \n";
        std::cout << "Quitting..." << std::endl;
        return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    int nprow = 2;
    int npcol = 2;

    const int m  = 21;
    const int nb = 6;

    {
        MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

        dist_matrix::BlacsContext bc(MPI_COMM_WORLD, nprow, npcol);
        dist_matrix::DistMatrix<DISTMATDTYPE>::setBlockSize(nb);

        dist_matrix::DistMatrix<DISTMATDTYPE>::setDefaultBlacsContext(&bc);

        dist_matrix::DistVector<double> distv("v", m);

        std::vector<double> stdv;
        for (int i = 0; i < m; i++)
            stdv.push_back(1.);

        distv.assign(stdv);

        double norm2          = distv.nrm2();
        double expected_norm2 = std::sqrt((double)m);

        if (mype == 0) std::cout << "Norm(v)=" << norm2 << std::endl;
        if (fabs(norm2 - expected_norm2) > 0.000001)
        {
            std::cout << "DistVector: problem with nrm2()" << std::endl;
            return 1;
        }

        distv.normalize();

        double normn = distv.nrm2();
        if (mype == 0) std::cout << "Norm(v)=" << normn << std::endl;
        if (fabs(normn - 1.) > 0.000001)
        {
            std::cout << "DistVector: problem with normalize()" << std::endl;
            return 1;
        }

        dist_matrix::DistVector<double> distv2(stdv);

        double diff2 = distv2.scaledDiff2(distv, norm2);
        if (mype == 0) std::cout << "diff2=" << diff2 << std::endl;
        if (fabs(diff2) > 0.000001)
        {
            std::cout << "DistVector: problem with scaledDiff2()" << std::endl;
            return 1;
        }
    }

    MPI_Finalize();

    if (mype == 0) std::cout << "Test SUCCESSFUL!" << std::endl;
    return 0;
}
