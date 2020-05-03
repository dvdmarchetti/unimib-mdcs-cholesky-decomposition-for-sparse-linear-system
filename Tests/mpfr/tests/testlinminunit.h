
#ifndef _testlinminunit_h
#define _testlinminunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "linmin.h"
namespace testlinminunit
{
    template<unsigned int Precision>
    bool testlinmin(bool silent);
    template<unsigned int Precision>
    bool testlinminunit_test_silent();
    template<unsigned int Precision>
    bool testlinminunit_test();


    template<unsigned int Precision>
    bool testlinmin(bool silent)
    {
        bool result;
        bool waserrors;


        waserrors = false;
        if( !silent )
        {
            printf("TESTING LINMIN\n");
            if( waserrors )
            {
                printf("TEST FAILED\n");
            }
            else
            {
                printf("TEST PASSED\n");
            }
            printf("\n\n");
        }
        result = !waserrors;
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testlinminunit_test_silent()
    {
        bool result;


        result = testlinmin<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testlinminunit_test()
    {
        bool result;


        result = testlinmin<Precision>(false);
        return result;
    }
} // namespace

#endif
