#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "testmatgenunit.h"

int main(int argc, char **argv)
{
    unsigned seed;
    if( argc==2 )
        seed = (unsigned)atoi(argv[1]);
    else
    {
        time_t t;
        seed = (unsigned)time(&t);
    }
    srand(seed);
    amp::mpfr_storage::seedRandState(seed);
    try
    {
        if(!testmatgenunit::testmatgenunit_test_silent<128>())
            throw 0;
    }
    catch(...)
    {
        printf("%-32s FAILED(seed=%ld)\n", "matgen", (long)seed);
        return 1;
    }
    printf("%-32s OK\n", "matgen");
    return 0;
}

