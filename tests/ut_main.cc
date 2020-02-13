#define CATCH_CONFIG_RUNNER
// Disable colored output because it doesn't play well with ctest
#define CATCH_CONFIG_COLOUR_NONE
#include "catch.hpp"

#include <mpi.h>

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int result = Catch::Session().run(argc, argv);

    MPI_Finalize();

    return result;
}
