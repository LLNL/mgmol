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
#include "MGmol_MPI.h"

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

int aa(int i, int j) { return i + 2 * j; }
int bb(int i, int j) { return i - j - 3; }

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

    const int n  = 21;
    const int nb = 6;

    int m_a  = n;
    int n_a  = n;
    int mb_a = nb;
    int nb_a = nb;

    int m_b  = n;
    int n_b  = n;
    int mb_b = nb;
    int nb_b = nb;

    int m_c  = n;
    int n_c  = n;
    int mb_c = nb;
    int nb_c = nb;

    char ta = 'n';
    char tb = 'n';

    {
        MGmol_MPI::setup(MPI_COMM_WORLD, std::cout);

        dist_matrix::BlacsContext bc(MPI_COMM_WORLD, nprow, npcol);

        dist_matrix::DistMatrix<double> a("a", bc, m_a, n_a, mb_a, nb_a);
        dist_matrix::DistMatrix<double> b("b", bc, m_b, n_b, mb_b, nb_b);
        dist_matrix::DistMatrix<double> c("c", bc, m_c, n_c, mb_c, nb_c);

        if (mype == 0)
        {
            std::cout << " m_a x n_a / mb_a x nb_a / ta = " << a.m() << "x"
                      << a.n() << " / " << a.mb() << "x" << a.nb() << " / "
                      << ta << std::endl;
            std::cout << " m_b x n_b / mb_b x nb_b / tb = " << b.m() << "x"
                      << b.n() << " / " << b.mb() << "x" << b.nb() << " / "
                      << tb << std::endl;
            std::cout << " m_c x n_c / mb_c x nb_c      = " << c.m() << "x"
                      << c.n() << " / " << c.mb() << "x" << c.nb() << std::endl;
        }

        for (int m = 0; m < a.nblocks(); m++)
            for (int l = 0; l < a.mblocks(); l++)
                for (int y = 0; y < a.nbs(m); y++)
                    for (int x = 0; x < a.mbs(l); x++)
                    {
                        int i      = a.i(l, x);
                        int j      = a.j(m, y);
                        double aij = aa(i, j);
                        int iii    = x + l * a.mb();
                        int jjj    = y + m * a.nb();
                        int ival   = iii + jjj * a.mloc();
                        a.setval(ival, aij);
                    }

        for (int m = 0; m < b.nblocks(); m++)
            for (int l = 0; l < b.mblocks(); l++)
                for (int y = 0; y < b.nbs(m); y++)
                    for (int x = 0; x < b.mbs(l); x++)
                    {
                        int i      = b.i(l, x);
                        int j      = b.j(m, y);
                        double bij = bb(i, j);
                        int iii    = x + l * b.mb();
                        int jjj    = y + m * b.nb();
                        int ival   = iii + jjj * b.mloc();
                        b.setval(ival, bij);
                    }

        c.gemm(ta, tb, 1.0, a, b, 0.0);

        // check result
        std::cout << " checking results..." << std::endl;
        for (int m = 0; m < c.nblocks(); m++)
            for (int l = 0; l < c.mblocks(); l++)
                for (int y = 0; y < c.nbs(m); y++)
                    for (int x = 0; x < c.mbs(l); x++)
                    {
                        int i      = c.i(l, x);
                        int j      = c.j(m, y);
                        double sum = 0.0;
                        int kmax   = (ta == 'n') ? a.n() : a.m();

                        if ((ta == 'n') && (tb == 'n'))
                        {
                            for (int k = 0; k < kmax; k++)
                                sum += aa(i, k) * bb(k, j);
                        }
                        else if ((ta != 'n') && (tb == 'n'))
                        {
                            for (int k = 0; k < kmax; k++)
                                sum += aa(k, i) * bb(k, j);
                        }
                        else if ((ta == 'n') && (tb != 'n'))
                        {
                            for (int k = 0; k < kmax; k++)
                                sum += aa(i, k) * bb(j, k);
                        }
                        else if ((ta != 'n') && (tb != 'n'))
                        {
                            for (int k = 0; k < kmax; k++)
                                sum += aa(k, i) * bb(j, k);
                        }

                        int iii  = x + l * c.mb();
                        int jjj  = y + m * c.nb();
                        int ival = iii + jjj * c.mloc();
                        if (fabs(c.val(ival) - sum) > 1.e-8)
                        {
                            std::cout << " error at element (" << i << "," << j
                                      << ") " << c.val(ival) << " " << sum
                                      << std::endl;
                            return 1;
                        }
                    }

        std::cout << " results checked" << std::endl;

        double norma = a.norm('F');
        if (mype == 0) std::cout << "Norm(a)=" << norma << std::endl;
        if (mype == 0) std::cout << "DistMatrix::allgather..." << std::endl;
        double* aa = new double[a.m() * a.n()];
        a.allgather(aa, a.m());
        if (mype == 0) std::cout << "DistMatrix::init..." << std::endl;
        b.init(aa, a.m());
        double norm = b.norm('F');
        if (mype == 0) std::cout << "Norm(b)=" << norm << std::endl;
        if (fabs(norm - norma) > 0.000001)
        {
            std::cout << "DistMatrix: problem with allgather/init" << std::endl;
            return -1;
        }
        if (mype == 0) std::cout << "DistMatrix::transpose..." << std::endl;
        c.transpose(b);
        norm = c.norm('F');
        if (mype == 0) std::cout << "Norm(c)=" << norm << std::endl;
        if (mype == 0) std::cout << "DistMatrix::scal..." << std::endl;
        c.scal(0.5);
        if (mype == 0) std::cout << "DistMatrix::transpose..." << std::endl;
        b.transpose(c);
        if (mype == 0) std::cout << "DistMatrix::axpy..." << std::endl;
        a.axpy(-2., b);
        if (mype == 0) std::cout << "DistMatrix::operator=..." << std::endl;
        c = a;
        if (mype == 0) std::cout << "DistMatrix::nrm2..." << std::endl;
        norm = c.norm('F');
        if (mype == 0) std::cout << "Norm=" << norm << std::endl;
    }

    MPI_Finalize();
}
