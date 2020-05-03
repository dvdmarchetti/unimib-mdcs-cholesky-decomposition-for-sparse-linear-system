
#ifndef _testgammaunit_h
#define _testgammaunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "gammafunc.h"
namespace testgammaunit
{
    template<unsigned int Precision>
    bool testgamma(bool silent);
    template<unsigned int Precision>
    bool testgammaunit_test_silent();
    template<unsigned int Precision>
    bool testgammaunit_test();


    template<unsigned int Precision>
    bool testgamma(bool silent)
    {
        bool result;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s;
        bool waserrors;
        bool gammaerrors;
        bool lngammaerrors;


        gammaerrors = false;
        lngammaerrors = false;
        waserrors = false;
        threshold = 100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        //
        //
        gammaerrors = gammaerrors || amp::abs<Precision>(gammafunc::gamma<Precision>(amp::ampf<Precision>("0.5"))-amp::sqrt<Precision>(amp::pi<Precision>()))>threshold;
        gammaerrors = gammaerrors || amp::abs<Precision>(gammafunc::gamma<Precision>(amp::ampf<Precision>("1.5"))-amp::ampf<Precision>("0.5")*amp::sqrt<Precision>(amp::pi<Precision>()))>threshold;
        v = gammafunc::lngamma<Precision>(amp::ampf<Precision>("0.5"), s);
        lngammaerrors = lngammaerrors || amp::abs<Precision>(v-amp::log<Precision>(amp::sqrt<Precision>(amp::pi<Precision>())))>threshold || s!=1;
        v = gammafunc::lngamma<Precision>(amp::ampf<Precision>("1.5"), s);
        lngammaerrors = lngammaerrors || amp::abs<Precision>(v-amp::log<Precision>(amp::ampf<Precision>("0.5")*amp::sqrt<Precision>(amp::pi<Precision>())))>threshold || s!=1;
        
        //
        // report
        //
        waserrors = gammaerrors || lngammaerrors;
        if( !silent )
        {
            printf("TESTING GAMMA FUNCTION\n");
            printf("GAMMA:                                   ");
            if( gammaerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("LN GAMMA:                                ");
            if( lngammaerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
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
        
        //
        // end
        //
        result = !waserrors;
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testgammaunit_test_silent()
    {
        bool result;


        result = testgamma<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testgammaunit_test()
    {
        bool result;


        result = testgamma<Precision>(false);
        return result;
    }
} // namespace

#endif
