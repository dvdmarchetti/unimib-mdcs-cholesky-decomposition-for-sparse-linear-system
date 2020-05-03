#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "testinvldltunit.h"

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
        if(!testinvldltunit::testinvldltunit_test_silent<128>())
            throw 0;
    }
    catch(...)
    {
        printf("SEED %9ld    UNIT %s\n", (long)seed, "sinverse");
        return 1;
    }
    return 0;
}

