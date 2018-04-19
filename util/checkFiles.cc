// mpiicpc -I/usr/local/tools/hdf5-intel-parallel-1.8.8/include -O2 -DNDEBUG -g -fp-model precise -fp-model source -o m.o -c checkFiles.cc
// mpiicpc checkFiles.o /usr/local/tools/hdf5-intel-parallel-1.8.8/lib/libhdf5_hl.a /usr/local/tools/hdf5-intel-parallel-1.8.8/lib/libhdf5.a -lz -lifcore -openmp  
//
// BGQ:
// 
// /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r -I/usr/local/tools/hdf5/hdf5-1.8.5/parallel/include -O2 -DNDEBUG -g -o checkFiles.o -c checkFiles.cc
// /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r -o checkFiles checkFiles.o /usr/local/tools/hdf5/hdf5-1.8.5/parallel/lib/libhdf5_hl.a /usr/local/tools/hdf5/hdf5-1.8.5/parallel/lib/libhdf5.a -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 -lxlsmp -L/opt/ibmcmp/xlmass/bg/7.3/bglib64 -L/opt/ibmcmp/xlf/bg/14.1/bglib64 -lxlf90_r -lxlopt -lxlomp_ser -lxl -lxlfmath -qsmp=omp -lm 
// 
// run:
// srun -ppdebug -n128 checkFiles /p/lscratchv/jeanluc/MGmol/1.5K_10_10_14/wave.out

#include "hdf5.h"
#include <mpi.h>

#include <sys/stat.h>
#include <sstream>
#include <string>
#include <iostream>

bool fileExists(const char* file)
{
    struct stat buf;
    return (stat(file, &buf) == 0);
}

void appendTaskNumberToFilename(std::string& filename)
{
    int mytask=0;
    int npes=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mytask);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    filename.append(".");
    if( mytask<1000000 && npes>999999 ){
        filename.append("0");
    }
    if( mytask<100000 && npes>99999 ){
        filename.append("0");
    }
    if( mytask<10000 && npes>9999 ){
        filename.append("0");
    }
    if( mytask<1000 ){
        filename.append("0");
    }
    if( mytask<100  ){
        filename.append("0");
    }
    if( mytask<10   ){
        filename.append("0");
    }

    std::stringstream oss("");
    oss << mytask;
    
    filename.append(oss.str());

    return;
}

int main(int argc, char **argv)
{
    int mpirc=MPI_Init(&argc, &argv);
    int mytask=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mytask);
    
    std::string filename(argv[1]);
    filename.append("/Task");
    appendTaskNumberToFilename(filename);
    
    int found=0;
    int foundHDF=0;
    int total=0;
    if( !fileExists(filename.c_str()) )
        std::cout<<"cannot see "<<filename<<std::endl;
    else
        found=1;
    if( htri_t ishdf=H5Fis_hdf5(filename.c_str())>=0 )
        foundHDF=1;
    
    MPI_Allreduce(&found, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if( mytask==0 )std::cout<<"Found "<<total<<" files."<<std::endl;
    
    MPI_Allreduce(&foundHDF, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if( mytask==0 )std::cout<<"Found "<<total<<" HDF files."<<std::endl;
    
    mpirc=MPI_Finalize();
    return 0;  
} 
