#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "testdensesolverunit.h"

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
        if(!testdensesolverunit::testdensesolverunit_test_silent<128>())
            throw 0;
    }
    catch(...)
    {
        printf("%-32s FAILED(seed=%ld)\n", "densesolver", (long)seed);
        return 1;
    }
    printf("%-32s OK\n", "densesolver");
    return 0;
}

