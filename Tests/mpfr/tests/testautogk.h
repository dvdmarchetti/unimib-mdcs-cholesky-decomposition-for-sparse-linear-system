
#ifndef _testautogk_h
#define _testautogk_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "tsort.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "gammafunc.h"
#include "gq.h"
#include "gkq.h"
#include "autogk.h"
namespace testautogk
{
    template<unsigned int Precision>
    bool testautogkunit(bool silent);
    template<unsigned int Precision>
    bool testautogk_test_silent();
    template<unsigned int Precision>
    bool testautogk_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testautogkunit(bool silent)
    {
        bool result;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        autogk::autogkstate<Precision> state;
        autogk::autogkreport<Precision> rep;
        amp::ampf<Precision> v;
        amp::ampf<Precision> exact;
        amp::ampf<Precision> eabs;
        amp::ampf<Precision> alpha;
        int pkind;
        amp::ampf<Precision> errtol;
        bool simpleerrors;
        bool sngenderrors;
        bool waserrors;


        simpleerrors = false;
        sngenderrors = false;
        waserrors = false;
        errtol = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Simple test: integral(exp(x),+-1,+-2), no maximum width requirements
        //
        a = (2*ap::randominteger(2)-1)*amp::ampf<Precision>("1.0");
        b = (2*ap::randominteger(2)-1)*amp::ampf<Precision>("2.0");
        autogk::autogksmooth<Precision>(a, b, state);
        while( autogk::autogkiteration<Precision>(state) )
        {
            state.f = amp::exp<Precision>(state.x);
        }
        autogk::autogkresults<Precision>(state, v, rep);
        exact = amp::exp<Precision>(b)-amp::exp<Precision>(a);
        eabs = amp::abs<Precision>(amp::exp<Precision>(b)-amp::exp<Precision>(a));
        if( rep.terminationtype<=0 )
        {
            simpleerrors = true;
        }
        else
        {
            simpleerrors = simpleerrors || amp::abs<Precision>(exact-v)>errtol*eabs;
        }
        
        //
        // Simple test: integral(exp(x),+-1,+-2), XWidth=0.1
        //
        a = (2*ap::randominteger(2)-1)*amp::ampf<Precision>("1.0");
        b = (2*ap::randominteger(2)-1)*amp::ampf<Precision>("2.0");
        autogk::autogksmoothw<Precision>(a, b, amp::ampf<Precision>("0.1"), state);
        while( autogk::autogkiteration<Precision>(state) )
        {
            state.f = amp::exp<Precision>(state.x);
        }
        autogk::autogkresults<Precision>(state, v, rep);
        exact = amp::exp<Precision>(b)-amp::exp<Precision>(a);
        eabs = amp::abs<Precision>(amp::exp<Precision>(b)-amp::exp<Precision>(a));
        if( rep.terminationtype<=0 )
        {
            simpleerrors = true;
        }
        else
        {
            simpleerrors = simpleerrors || amp::abs<Precision>(exact-v)>errtol*eabs;
        }
        
        //
        // Simple test: integral(cos(100*x),0,2*pi), no maximum width requirements
        //
        a = 0;
        b = 2*amp::pi<Precision>();
        autogk::autogksmooth<Precision>(a, b, state);
        while( autogk::autogkiteration<Precision>(state) )
        {
            state.f = amp::cos<Precision>(100*state.x);
        }
        autogk::autogkresults<Precision>(state, v, rep);
        exact = 0;
        eabs = 4;
        if( rep.terminationtype<=0 )
        {
            simpleerrors = true;
        }
        else
        {
            simpleerrors = simpleerrors || amp::abs<Precision>(exact-v)>errtol*eabs;
        }
        
        //
        // Simple test: integral(cos(100*x),0,2*pi), XWidth=0.3
        //
        a = 0;
        b = 2*amp::pi<Precision>();
        autogk::autogksmoothw<Precision>(a, b, amp::ampf<Precision>("0.3"), state);
        while( autogk::autogkiteration<Precision>(state) )
        {
            state.f = amp::cos<Precision>(100*state.x);
        }
        autogk::autogkresults<Precision>(state, v, rep);
        exact = 0;
        eabs = 4;
        if( rep.terminationtype<=0 )
        {
            simpleerrors = true;
        }
        else
        {
            simpleerrors = simpleerrors || amp::abs<Precision>(exact-v)>errtol*eabs;
        }
        
        //
        // singular problem on [a,b] = [0.1, 0.5]
        //     f2(x) = (1+x)*(b-x)^alpha, -1 < alpha < 1
        //
        for(pkind=0; pkind<=6; pkind++)
        {
            a = amp::ampf<Precision>("0.1");
            b = amp::ampf<Precision>("0.5");
            if( pkind==0 )
            {
                alpha = -amp::ampf<Precision>("0.9");
            }
            if( pkind==1 )
            {
                alpha = -amp::ampf<Precision>("0.5");
            }
            if( pkind==2 )
            {
                alpha = -amp::ampf<Precision>("0.1");
            }
            if( pkind==3 )
            {
                alpha = amp::ampf<Precision>("0.0");
            }
            if( pkind==4 )
            {
                alpha = amp::ampf<Precision>("0.1");
            }
            if( pkind==5 )
            {
                alpha = amp::ampf<Precision>("0.5");
            }
            if( pkind==6 )
            {
                alpha = amp::ampf<Precision>("0.9");
            }
            
            //
            // f1(x) = (1+x)*(x-a)^alpha, -1 < alpha < 1
            // 1. use singular integrator for [a,b]
            // 2. use singular integrator for [b,a]
            //
            exact = amp::pow<Precision>(b-a, alpha+2)/(alpha+2)+(1+a)*amp::pow<Precision>(b-a, alpha+1)/(alpha+1);
            eabs = amp::abs<Precision>(exact);
            autogk::autogksingular<Precision>(a, b, alpha, amp::ampf<Precision>("0.0"), state);
            while( autogk::autogkiteration<Precision>(state) )
            {
                if( state.xminusa<amp::ampf<Precision>("0.01") )
                {
                    state.f = amp::pow<Precision>(state.xminusa, alpha)*(1+state.x);
                }
                else
                {
                    state.f = amp::pow<Precision>(state.x-a, alpha)*(1+state.x);
                }
            }
            autogk::autogkresults<Precision>(state, v, rep);
            if( rep.terminationtype<=0 )
            {
                sngenderrors = true;
            }
            else
            {
                sngenderrors = sngenderrors || amp::abs<Precision>(v-exact)>errtol*eabs;
            }
            autogk::autogksingular<Precision>(b, a, amp::ampf<Precision>("0.0"), alpha, state);
            while( autogk::autogkiteration<Precision>(state) )
            {
                if( state.bminusx>-amp::ampf<Precision>("0.01") )
                {
                    state.f = amp::pow<Precision>(-state.bminusx, alpha)*(1+state.x);
                }
                else
                {
                    state.f = amp::pow<Precision>(state.x-a, alpha)*(1+state.x);
                }
            }
            autogk::autogkresults<Precision>(state, v, rep);
            if( rep.terminationtype<=0 )
            {
                sngenderrors = true;
            }
            else
            {
                sngenderrors = sngenderrors || amp::abs<Precision>(-v-exact)>errtol*eabs;
            }
            
            //
            // f1(x) = (1+x)*(b-x)^alpha, -1 < alpha < 1
            // 1. use singular integrator for [a,b]
            // 2. use singular integrator for [b,a]
            //
            exact = (1+b)*amp::pow<Precision>(b-a, alpha+1)/(alpha+1)-amp::pow<Precision>(b-a, alpha+2)/(alpha+2);
            eabs = amp::abs<Precision>(exact);
            autogk::autogksingular<Precision>(a, b, amp::ampf<Precision>("0.0"), alpha, state);
            while( autogk::autogkiteration<Precision>(state) )
            {
                if( state.bminusx<amp::ampf<Precision>("0.01") )
                {
                    state.f = amp::pow<Precision>(state.bminusx, alpha)*(1+state.x);
                }
                else
                {
                    state.f = amp::pow<Precision>(b-state.x, alpha)*(1+state.x);
                }
            }
            autogk::autogkresults<Precision>(state, v, rep);
            if( rep.terminationtype<=0 )
            {
                sngenderrors = true;
            }
            else
            {
                sngenderrors = sngenderrors || amp::abs<Precision>(v-exact)>errtol*eabs;
            }
            autogk::autogksingular<Precision>(b, a, alpha, amp::ampf<Precision>("0.0"), state);
            while( autogk::autogkiteration<Precision>(state) )
            {
                if( state.xminusa>-amp::ampf<Precision>("0.01") )
                {
                    state.f = amp::pow<Precision>(-state.xminusa, alpha)*(1+state.x);
                }
                else
                {
                    state.f = amp::pow<Precision>(b-state.x, alpha)*(1+state.x);
                }
            }
            autogk::autogkresults<Precision>(state, v, rep);
            if( rep.terminationtype<=0 )
            {
                sngenderrors = true;
            }
            else
            {
                sngenderrors = sngenderrors || amp::abs<Precision>(-v-exact)>errtol*eabs;
            }
        }
        
        //
        // end
        //
        waserrors = simpleerrors || sngenderrors;
        if( !silent )
        {
            printf("TESTING AUTOGK\n");
            printf("INTEGRATION WITH GIVEN ACCURACY:          ");
            if( simpleerrors || sngenderrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* SIMPLE PROBLEMS:                        ");
            if( simpleerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* SINGULAR PROBLEMS (ENDS OF INTERVAL):   ");
            if( sngenderrors )
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
        result = !waserrors;
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testautogk_test_silent()
    {
        bool result;


        result = testautogkunit<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testautogk_test()
    {
        bool result;


        result = testautogkunit<Precision>(false);
        return result;
    }
} // namespace

#endif
