#define CATCH_CONFIG_RUNNER
// Disable colored output because it doesn't play well with ctest
#define CATCH_CONFIG_COLOUR_NONE
#include "catch.hpp"

int main(int argc, char* argv[])
{
    int result = Catch::Session().run(argc, argv);

    return result;
}
