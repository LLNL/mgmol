#define CATCH_CONFIG_RUNNER
// Disable colored output because it doesn't play well with ctest
#define CATCH_CONFIG_COLOUR_NONE
#include "catch.hpp"

#include <mpi.h>

#ifdef HAVE_MAGMA
#include "magma_v2.h"
#endif

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

#ifdef HAVE_MAGMA
    magma_int_t magmalog = magma_init();
    if (magmalog == MAGMA_SUCCESS)
    {
        std::cout << "MAGMA Initialization: success" << std::endl;
    }
    else
    {
        if (magmalog == MAGMA_ERR_UNKNOWN)
            std::cout << "MAGMA Initialization: unknown error" << std::endl;
        if (magmalog == MAGMA_ERR_HOST_ALLOC)
            std::cout << "MAGMA Initialization: fails host alloc" << std::endl;
        return 1;
    }
#endif

    int result = Catch::Session().run(argc, argv);

#ifdef HAVE_MAGMA
    magmalog = magma_finalize();

    if (magmalog == MAGMA_SUCCESS)
    {
        std::cout << "MAGMA Finalize: success" << std::endl;
    }
    else
    {
        if (magmalog == MAGMA_ERR_UNKNOWN)
            std::cout << "MAGMA Finalize: unknown error" << std::endl;
        if (magmalog == MAGMA_ERR_HOST_ALLOC)
            std::cout << "MAGMA FINALIZE: fails host alloc" << std::endl;
        return 1;
    }
#endif

    MPI_Finalize();

    return result;
}
