#include "ReplicatedMatrix.h"

#include "catch.hpp"
#include <mpi.h>

TEST_CASE("Check ReplicatedMatrix", "[replicated_matrix]")
{
    int npes;
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    ReplicatedMatrix::setMPIcomm(MPI_COMM_WORLD);

    int n = 10;

    // diagonalize diagonal matrix
    {
        std::vector<double> diag(n);
        for (int i = 0; i < n; i++)
            diag[i] = static_cast<double>(i);

        ReplicatedMatrix mat("diag", diag.data(), n, n);
        ReplicatedMatrix evects("evects", n, n);
        std::vector<double> evals(n);
        mat.syev('T', 'U', evals, evects);

        for (int i = 0; i < n; i++)
            CHECK(evals[i] == Approx(diag[i]).epsilon(1.e-8));

        double tr = mat.trace();
        CHECK(
            tr == Approx(static_cast<double>(n * (n - 1) / 2)).epsilon(1.e-8));
    }

    // set values using SquareLocalMatrices
    // then consolidate matrix
    {
        double value = 0.1;
        SquareLocalMatrices<double, MemorySpace::Host> slm(1, n);
        std::vector<double> tmp(n * n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                tmp[i + n * j] = value;

        slm.setValues(tmp.data(), n);

        ReplicatedMatrix mat("mat", n, n);
        mat.assign(slm);
        mat.consolidate();
        std::vector<double> values(n * n);
        mat.get(values.data(), n);
        for (int i = 0; i < n * n; i++)
            CHECK(values[i]
                  == Approx(value * static_cast<double>(npes)).epsilon(1.e-8));
    }

    {
        assert(0 == 1);
        using MemoryDev = MemorySpace::Memory<double, MemorySpace::Device>;

        double value = 0.1;
        SquareLocalMatrices<double, MemorySpace::Device> slm(1, n);
        std::vector<double> tmp(n * n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                tmp[i + n * j] = value;

        std::unique_ptr<double, void (*)(double*)> tmp_dev(
            MemoryDev::allocate(n * n), MemoryDev::free);

        auto& magma_singleton = MagmaSingleton::get_magma_singleton();
        magma_dsetmatrix(
            n, n, tmp.data(), n, tmp_dev.get(), n, magma_singleton.queue_);

        slm.setValues(tmp_dev.get(), n);

        ReplicatedMatrix mat("mat", n, n);
        mat.assign(slm);
        mat.consolidate();
        std::vector<double> values(n * n);
        mat.get(values.data(), n);
        for (int i = 0; i < n * n; i++)
            CHECK(values[i]
                  == Approx(value * static_cast<double>(npes)).epsilon(1.e-8));
    }

    // set diagonal, then retrieve it
    {
        ReplicatedMatrix mat("mat", n, n);
        std::vector<double> values(n);
        for (int i = 0; i < n; i++)
            values[i] = i;
        mat.setDiagonal(values);

        std::vector<double> newval(n);
        mat.getDiagonalValues(newval.data());
        for (int i = 0; i < n; i++)
            CHECK(newval[i] == Approx(static_cast<double>(i)).epsilon(1.e-15));
    }

    // set matrix to identity, then retrieve diagonal
    {
        ReplicatedMatrix mat("mat", n, n);
        mat.identity();
        std::vector<double> newval(n);
        mat.getDiagonalValues(newval.data());
        for (int i = 0; i < n; i++)
            CHECK(newval[i] == Approx(1.).epsilon(1.e-15));
    }
}
