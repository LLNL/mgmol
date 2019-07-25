//
// test Replicated2DistMatrix
//
#include "BlacsContext.h"
#include "DistMatrix.h"
#include "SquareLocalMatrices.h"
#include "MGmol_MPI.h"

#include <mpi.h>

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>


int aa(int i, int j) { return i+2*j; }
int bb(int i, int j) { return i-j-3; }

int main(int argc, char **argv)
{
    MPI_Init(&argc,&argv);

    // MPI scope
    {
    int mype;
    int npes;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    if (npes != 4) {
       std::cout<<"Number of processors in MPI: "<<npes<<std::endl;
       std::cout<<"This example to set up to use only 4 processes \n";
       std::cout<<"Quitting..."<<std::endl;
       return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    int nprow=2;
    int npcol=2;

    const int n=11;
    const int nb=3;

    int m=n;
    int mb_a=nb;
    int nb_a=nb;

    MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

    dist_matrix::BlacsContext bc(MPI_COMM_WORLD, nprow,npcol);

    dist_matrix::DistMatrix<double> distm("a",bc,m,n,mb_a,nb_a);

    if ( mype == 0 )
    {
      std::cout << " m x n / mb x nb = "
           << distm.m() << "x" << distm.n() << " / "
           << distm.mb() << "x" << distm.nb() << " / " << std::endl;
    }

    // setup an nxn local matrix
    SquareLocalMatrices<double> replicated(1, n);
    for ( int i = 0; i < m; i++ )
    for ( int j = 0; j < n; j++ )
        replicated.setVal(i, j, 10.*(i+1)+j+1);

    // distribute replicated matrix
    distm.initFromReplicated(replicated.getSubMatrix(), n);

    // convert back to a replicated matrix
    SquareLocalMatrices<double> result(1, n);
    distm.allgather(result.getSubMatrix(), n);

    if(mype==0)
    {
        std::cout << std::setprecision(2);
        result.print(std::cout);
    }

    const double tol = 1.e-8;
    for ( int i = 0; i < m; i++ )
    for ( int j = 0; j < n; j++ )
    {
        double valbefore = replicated.getVal(i, j);
        double valafter = result.getVal(i, j);
        if ( abs(valbefore-valafter) > tol )
        {
            std::cerr<<"i="<<i<<", j="<<j
                     <<", val before: "<<valbefore<<", val. after: "<<valafter
                     <<std::endl;
            return 1;
        }
        // make sure values are not 0
        if ( abs(valafter) < tol )
        {
            std::cerr<<"i="<<i<<", j="<<j
                     <<", val before: "<<valbefore<<", val. after: "<<valafter
                     <<std::endl;
            return 1;
        }
    }

    }
    MPI_Finalize();

    std::cout<<"Successful test!"<<std::endl;
    return 0;
}
