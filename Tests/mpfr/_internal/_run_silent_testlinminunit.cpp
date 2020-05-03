#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "testlinminunit.h"

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
        if(!testlinminunit::testlinminunit_test_silent<128>())
            throw 0;
    }
    catch(...)
    {
        printf("SEED %9ld    UNIT %s\n", (long)seed, "linmin");
        return 1;
    }
    return 0;
}

