#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>

#include <mpi.h>

struct ExecutionEnvironmentScopeGuard
{
    ExecutionEnvironmentScopeGuard(int argc, char* argv[])
    {
        MPI_Init(&argc, &argv);
    }

    ~ExecutionEnvironmentScopeGuard() { MPI_Finalize(); }
};

bool init_function() { return true; }

int main(int argc, char* argv[])
{
    ExecutionEnvironmentScopeGuard scope_guard(argc, argv);

    return boost::unit_test::unit_test_main(&init_function, argc, argv);
}
