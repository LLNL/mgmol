#include "PEenv.h"
#include "DirectionalReduce.h"

#include <iostream>
using namespace std;

int main(int argc, char **argv)
{
    int mpirc=MPI_Init(&argc, &argv);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if( myrank == 0 )
        cout<<"Test DirectionalReduce"<<endl;

    {
    //mesh
    int ngpts[3]={128,256,196};

    //mesh spacing
    double hh=0.1;

    //distance for data gathering
    double radius=4.;

    pb::PEenv myPEenv(MPI_COMM_WORLD, ngpts[0],ngpts[1],ngpts[2], 1);

    double domain[3]={ngpts[0]*hh,ngpts[1]*hh,ngpts[2]*hh};

    DirectionalReduce dir_reduce(myPEenv.cart_comm(), radius, domain);

    int nprocs[3];
    int periodic[3];
    int coords[3];
    int error = MPI_Cart_get( myPEenv.cart_comm(), 3, nprocs, periodic, coords);

    //direction for test
    const short dir=1;

    if( myrank == 0 )
        cout<<"nprocs["<<dir<<"]="<<nprocs[dir]<<endl;
    if( myrank == 0 )
        cout<<"rstep["<<dir<<"]="<<dir_reduce.rstep(dir)<<endl;

    //compute max. coord in direction dir with radius 'radius'
    int data[2]={coords[dir],1};
    dir_reduce.computeDirMax(dir,data);

    int result = std::min(coords[dir]+dir_reduce.rstep(dir),nprocs[dir]-1);
    if( dir_reduce.lstep(dir) > coords[dir] )result=nprocs[dir]-1;

    std::cout<<"myrank = "<<myrank <<", coords[dir] = "<<coords[dir]<<", data[0] = "<<data[0]<<endl;

    if( result != data[0] )
    {
        std::cout<<"myrank="<<myrank <<", (data[0]="<<data[0]<<") !=  (result="<<result<<")\n";
        return 1;
    }

    }

    mpirc=MPI_Finalize();

    //return 0 for SUCCESS
    return 0;
}

