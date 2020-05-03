
#ifndef _testminlmunit_h
#define _testminlmunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "blas.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "linmin.h"
#include "minlbfgs.h"
#include "minlm.h"
namespace testminlmunit
{
    template<unsigned int Precision>
    bool testminlm(bool silent);
    template<unsigned int Precision>
    bool rkindvsstatecheck(int rkind,
        const minlm::minlmstate<Precision>& state);
    template<unsigned int Precision>
    bool testminlmunit_test_silent();
    template<unsigned int Precision>
    bool testminlmunit_test();


    template<unsigned int Precision>
    bool testminlm(bool silent)
    {
        bool result;
        bool waserrors;
        bool referror;
        bool lin1error;
        bool lin2error;
        bool eqerror;
        bool converror;
        bool scerror;
        bool othererrors;
        int rkind;
        int ckind;
        amp::ampf<Precision> epsf;
        amp::ampf<Precision> epsx;
        amp::ampf<Precision> epsg;
        int maxits;
        int n;
        int m;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > xe;
        ap::template_1d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > xlast;
        int i;
        int j;
        int k;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s;
        amp::ampf<Precision> stpmax;
        ap::template_2d_array< amp::ampf<Precision> > a;
        amp::ampf<Precision> fprev;
        amp::ampf<Precision> xprev;
        minlm::minlmstate<Precision> state;
        minlm::minlmreport<Precision> rep;


        waserrors = false;
        referror = false;
        lin1error = false;
        lin2error = false;
        eqerror = false;
        converror = false;
        scerror = false;
        othererrors = false;
        
        //
        // Reference problem.
        // RKind is a algorithm selector:
        // * 0 = FJ
        // * 1 = FGJ
        // * 2 = FGH
        //
        x.setbounds(0, 2);
        n = 3;
        m = 3;
        for(rkind=0; rkind<=2; rkind++)
        {
            x(0) = 100*amp::ampf<Precision>::getRandom()-50;
            x(1) = 100*amp::ampf<Precision>::getRandom()-50;
            x(2) = 100*amp::ampf<Precision>::getRandom()-50;
            if( rkind==0 )
            {
                minlm::minlmcreatefj<Precision>(n, m, x, state);
            }
            if( rkind==1 )
            {
                minlm::minlmcreatefgj<Precision>(n, m, x, state);
            }
            if( rkind==2 )
            {
                minlm::minlmcreatefgh<Precision>(n, x, state);
            }
            while( minlm::minlmiteration<Precision>(state) )
            {
                
                //
                // (x-2)^2 + y^2 + (z-x)^2
                //
                state.f = amp::sqr<Precision>(state.x(0)-2)+amp::sqr<Precision>(state.x(1))+amp::sqr<Precision>(state.x(2)-state.x(0));
                if( state.needfg || state.needfgh )
                {
                    state.g(0) = 2*(state.x(0)-2)+2*(state.x(0)-state.x(2));
                    state.g(1) = 2*state.x(1);
                    state.g(2) = 2*(state.x(2)-state.x(0));
                }
                if( state.needfij )
                {
                    state.fi(0) = state.x(0)-2;
                    state.fi(1) = state.x(1);
                    state.fi(2) = state.x(2)-state.x(0);
                    state.j(0,0) = 1;
                    state.j(0,1) = 0;
                    state.j(0,2) = 0;
                    state.j(1,0) = 0;
                    state.j(1,1) = 1;
                    state.j(1,2) = 0;
                    state.j(2,0) = -1;
                    state.j(2,1) = 0;
                    state.j(2,2) = 1;
                }
                if( state.needfgh )
                {
                    state.h(0,0) = 4;
                    state.h(0,1) = 0;
                    state.h(0,2) = -2;
                    state.h(1,0) = 0;
                    state.h(1,1) = 2;
                    state.h(1,2) = 0;
                    state.h(2,0) = -2;
                    state.h(2,1) = 0;
                    state.h(2,2) = 2;
                }
                scerror = scerror || !rkindvsstatecheck<Precision>(rkind, state);
            }
            minlm::minlmresults<Precision>(state, x, rep);
            referror = referror || rep.terminationtype<=0 || amp::abs<Precision>(x(0)-2)>amp::ampf<Precision>("0.001") || amp::abs<Precision>(x(1))>amp::ampf<Precision>("0.001") || amp::abs<Precision>(x(2)-2)>amp::ampf<Precision>("0.001");
        }
        
        //
        // 1D problem #1
        //
        for(rkind=0; rkind<=2; rkind++)
        {
            x.setlength(1);
            n = 1;
            m = 1;
            x(0) = 100*amp::ampf<Precision>::getRandom()-50;
            if( rkind==0 )
            {
                minlm::minlmcreatefj<Precision>(n, m, x, state);
            }
            if( rkind==1 )
            {
                minlm::minlmcreatefgj<Precision>(n, m, x, state);
            }
            if( rkind==2 )
            {
                minlm::minlmcreatefgh<Precision>(n, x, state);
            }
            while( minlm::minlmiteration<Precision>(state) )
            {
                state.f = amp::sqr<Precision>(amp::sin<Precision>(state.x(0)));
                if( state.needfg || state.needfgh )
                {
                    state.g(0) = 2*amp::sin<Precision>(state.x(0))*amp::cos<Precision>(state.x(0));
                }
                if( state.needfij )
                {
                    state.fi(0) = amp::sin<Precision>(state.x(0));
                    state.j(0,0) = amp::cos<Precision>(state.x(0));
                }
                if( state.needfgh )
                {
                    state.h(0,0) = 2*(amp::cos<Precision>(state.x(0))*amp::cos<Precision>(state.x(0))-amp::sin<Precision>(state.x(0))*amp::sin<Precision>(state.x(0)));
                }
                scerror = scerror || !rkindvsstatecheck<Precision>(rkind, state);
            }
            minlm::minlmresults<Precision>(state, x, rep);
            lin1error = rep.terminationtype<=0 || amp::abs<Precision>(x(0)/amp::pi<Precision>()-amp::round<Precision>(x(0)/amp::pi<Precision>()))>amp::ampf<Precision>("0.001");
        }
        
        //
        // Linear equations
        //
        for(n=1; n<=10; n++)
        {
            
            //
            // Prepare task
            //
            matgen::rmatrixrndcond<Precision>(n, amp::ampf<Precision>(100), a);
            x.setlength(n);
            xe.setlength(n);
            b.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                xe(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            for(i=0; i<=n-1; i++)
            {
                v = amp::vdotproduct(a.getrow(i, 0, n-1), xe.getvector(0, n-1));
                b(i) = v;
            }
            
            //
            // Test different RKind
            //
            for(rkind=0; rkind<=2; rkind++)
            {
                
                //
                // Solve task
                //
                for(i=0; i<=n-1; i++)
                {
                    x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                if( rkind==0 )
                {
                    minlm::minlmcreatefj<Precision>(n, n, x, state);
                }
                if( rkind==1 )
                {
                    minlm::minlmcreatefgj<Precision>(n, n, x, state);
                }
                if( rkind==2 )
                {
                    minlm::minlmcreatefgh<Precision>(n, x, state);
                }
                while( minlm::minlmiteration<Precision>(state) )
                {
                    if( state.needf || state.needfg || state.needfgh )
                    {
                        state.f = 0;
                    }
                    if( state.needfg || state.needfgh )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            state.g(i) = 0;
                        }
                    }
                    if( state.needfgh )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                state.h(i,j) = 0;
                            }
                        }
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = amp::vdotproduct(a.getrow(i, 0, n-1), state.x.getvector(0, n-1));
                        if( state.needf || state.needfg || state.needfgh )
                        {
                            state.f = state.f+amp::sqr<Precision>(v-b(i));
                        }
                        if( state.needfg || state.needfgh )
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                state.g(j) = state.g(j)+2*(v-b(i))*a(i,j);
                            }
                        }
                        if( state.needfgh )
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                for(k=0; k<=n-1; k++)
                                {
                                    state.h(j,k) = state.h(j,k)+2*a(i,j)*a(i,k);
                                }
                            }
                        }
                        if( state.needfij )
                        {
                            state.fi(i) = v-b(i);
                            amp::vmove(state.j.getrow(i, 0, n-1), a.getrow(i, 0, n-1));
                        }
                    }
                    scerror = scerror || !rkindvsstatecheck<Precision>(rkind, state);
                }
                minlm::minlmresults<Precision>(state, x, rep);
                eqerror = eqerror || rep.terminationtype<=0;
                for(i=0; i<=n-1; i++)
                {
                    eqerror = eqerror || amp::abs<Precision>(x(i)-xe(i))>amp::ampf<Precision>("0.001");
                }
            }
        }
        
        //
        // Testing convergence properties using
        // different optimizer types and different conditions
        //
        s = 100;
        for(rkind=0; rkind<=2; rkind++)
        {
            for(ckind=0; ckind<=3; ckind++)
            {
                epsg = 0;
                epsf = 0;
                epsx = 0;
                maxits = 0;
                if( ckind==0 )
                {
                    epsf = amp::ampf<Precision>("0.0001");
                }
                if( ckind==1 )
                {
                    epsx = amp::ampf<Precision>("0.0001");
                }
                if( ckind==2 )
                {
                    maxits = 2;
                }
                if( ckind==3 )
                {
                    epsg = amp::ampf<Precision>("0.0001");
                }
                x.setlength(3);
                n = 3;
                m = 3;
                for(i=0; i<=2; i++)
                {
                    x(i) = 6;
                }
                if( rkind==0 )
                {
                    minlm::minlmcreatefj<Precision>(n, m, x, state);
                }
                if( rkind==1 )
                {
                    minlm::minlmcreatefgj<Precision>(n, m, x, state);
                }
                if( rkind==2 )
                {
                    minlm::minlmcreatefgh<Precision>(n, x, state);
                }
                minlm::minlmsetcond<Precision>(state, epsg, epsf, epsx, maxits);
                while( minlm::minlmiteration<Precision>(state) )
                {
                    if( state.needf || state.needfg || state.needfgh )
                    {
                        state.f = s*amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(amp::sqr<Precision>(state.x(1))+1)+amp::sqr<Precision>(state.x(2)-state.x(0));
                    }
                    if( state.needfg || state.needfgh )
                    {
                        state.g(0) = s*2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
                        state.g(1) = 2*(amp::sqr<Precision>(state.x(1))+1)*2*state.x(1);
                        state.g(2) = 2*(state.x(2)-state.x(0));
                    }
                    if( state.needfgh )
                    {
                        state.h(0,0) = s*(4*amp::sqr<Precision>(amp::exp<Precision>(state.x(0)))-4*amp::exp<Precision>(state.x(0)))+2;
                        state.h(0,1) = 0;
                        state.h(0,2) = -2;
                        state.h(1,0) = 0;
                        state.h(1,1) = 12*amp::sqr<Precision>(state.x(1))+4;
                        state.h(1,2) = 0;
                        state.h(2,0) = -2;
                        state.h(2,1) = 0;
                        state.h(2,2) = 2;
                    }
                    if( state.needfij )
                    {
                        state.fi(0) = s*(amp::exp<Precision>(state.x(0))-2);
                        state.j(0,0) = s*amp::exp<Precision>(state.x(0));
                        state.j(0,1) = 0;
                        state.j(0,2) = 0;
                        state.fi(1) = amp::sqr<Precision>(state.x(1))+1;
                        state.j(1,0) = 0;
                        state.j(1,1) = 2*state.x(1);
                        state.j(1,2) = 0;
                        state.fi(2) = state.x(2)-state.x(0);
                        state.j(2,0) = -1;
                        state.j(2,1) = 0;
                        state.j(2,2) = 1;
                    }
                    scerror = scerror || !rkindvsstatecheck<Precision>(rkind, state);
                }
                minlm::minlmresults<Precision>(state, x, rep);
                if( ckind==0 )
                {
                    converror = converror || amp::abs<Precision>(x(0)-amp::log<Precision>(amp::ampf<Precision>(2)))>amp::ampf<Precision>("0.05");
                    converror = converror || amp::abs<Precision>(x(1))>amp::ampf<Precision>("0.05");
                    converror = converror || amp::abs<Precision>(x(2)-amp::log<Precision>(amp::ampf<Precision>(2)))>amp::ampf<Precision>("0.05");
                    converror = converror || rep.terminationtype!=1;
                }
                if( ckind==1 )
                {
                    converror = converror || amp::abs<Precision>(x(0)-amp::log<Precision>(amp::ampf<Precision>(2)))>amp::ampf<Precision>("0.05");
                    converror = converror || amp::abs<Precision>(x(1))>amp::ampf<Precision>("0.05");
                    converror = converror || amp::abs<Precision>(x(2)-amp::log<Precision>(amp::ampf<Precision>(2)))>amp::ampf<Precision>("0.05");
                    converror = converror || rep.terminationtype!=2;
                }
                if( ckind==2 )
                {
                    converror = converror || rep.terminationtype!=5 || rep.iterationscount!=maxits;
                }
                if( ckind==3 )
                {
                    converror = converror || amp::abs<Precision>(x(0)-amp::log<Precision>(amp::ampf<Precision>(2)))>amp::ampf<Precision>("0.05");
                    converror = converror || amp::abs<Precision>(x(1))>amp::ampf<Precision>("0.05");
                    converror = converror || amp::abs<Precision>(x(2)-amp::log<Precision>(amp::ampf<Precision>(2)))>amp::ampf<Precision>("0.05");
                    converror = converror || rep.terminationtype!=4;
                }
            }
        }
        
        //
        // Other properties:
        // 1. test reports (F should form monotone sequence)
        // 2. test maximum step
        //
        for(rkind=0; rkind<=2; rkind++)
        {
            
            //
            // reports:
            // * check that first report is initial point
            // * check that F is monotone decreasing
            // * check that last report is final result
            //
            n = 3;
            m = 3;
            s = 100;
            x.setlength(n);
            xlast.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                x(i) = 6;
            }
            if( rkind==0 )
            {
                minlm::minlmcreatefj<Precision>(n, m, x, state);
            }
            if( rkind==1 )
            {
                minlm::minlmcreatefgj<Precision>(n, m, x, state);
            }
            if( rkind==2 )
            {
                minlm::minlmcreatefgh<Precision>(n, x, state);
            }
            minlm::minlmsetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 4);
            minlm::minlmsetxrep<Precision>(state, true);
            fprev = amp::ampf<Precision>::getAlgoPascalMaxNumber();
            while( minlm::minlmiteration<Precision>(state) )
            {
                if( state.needf || state.needfg || state.needfgh )
                {
                    state.f = s*amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(state.x(1))+amp::sqr<Precision>(state.x(2)-state.x(0));
                }
                if( state.needfg || state.needfgh )
                {
                    state.g(0) = s*2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
                    state.g(1) = 2*state.x(1);
                    state.g(2) = 2*(state.x(2)-state.x(0));
                }
                if( state.needfgh )
                {
                    state.h(0,0) = s*(4*amp::sqr<Precision>(amp::exp<Precision>(state.x(0)))-4*amp::exp<Precision>(state.x(0)))+2;
                    state.h(0,1) = 0;
                    state.h(0,2) = -2;
                    state.h(1,0) = 0;
                    state.h(1,1) = 2;
                    state.h(1,2) = 0;
                    state.h(2,0) = -2;
                    state.h(2,1) = 0;
                    state.h(2,2) = 2;
                }
                if( state.needfij )
                {
                    state.fi(0) = amp::sqrt<Precision>(s)*(amp::exp<Precision>(state.x(0))-2);
                    state.j(0,0) = amp::sqrt<Precision>(s)*amp::exp<Precision>(state.x(0));
                    state.j(0,1) = 0;
                    state.j(0,2) = 0;
                    state.fi(1) = state.x(1);
                    state.j(1,0) = 0;
                    state.j(1,1) = 1;
                    state.j(1,2) = 0;
                    state.fi(2) = state.x(2)-state.x(0);
                    state.j(2,0) = -1;
                    state.j(2,1) = 0;
                    state.j(2,2) = 1;
                }
                scerror = scerror || !rkindvsstatecheck<Precision>(rkind, state);
                if( state.xupdated )
                {
                    othererrors = othererrors || state.f>fprev;
                    if( fprev==amp::ampf<Precision>::getAlgoPascalMaxNumber() )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            othererrors = othererrors || state.x(i)!=x(i);
                        }
                    }
                    fprev = state.f;
                    amp::vmove(xlast.getvector(0, n-1), state.x.getvector(0, n-1));
                }
            }
            minlm::minlmresults<Precision>(state, x, rep);
            for(i=0; i<=n-1; i++)
            {
                othererrors = othererrors || x(i)!=xlast(i);
            }
        }
        n = 1;
        x.setlength(n);
        x(0) = 100;
        stpmax = amp::ampf<Precision>("0.05")+amp::ampf<Precision>("0.05")*amp::ampf<Precision>::getRandom();
        minlm::minlmcreatefgh<Precision>(n, x, state);
        minlm::minlmsetcond<Precision>(state, amp::ampf<Precision>("1.0E-9"), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
        minlm::minlmsetstpmax<Precision>(state, stpmax);
        minlm::minlmsetxrep<Precision>(state, true);
        xprev = x(0);
        while( minlm::minlmiteration<Precision>(state) )
        {
            if( state.needf || state.needfg || state.needfgh )
            {
                state.f = amp::exp<Precision>(state.x(0))+amp::exp<Precision>(-state.x(0));
            }
            if( state.needfg || state.needfgh )
            {
                state.g(0) = amp::exp<Precision>(state.x(0))-amp::exp<Precision>(-state.x(0));
            }
            if( state.needfgh )
            {
                state.h(0,0) = amp::exp<Precision>(state.x(0))+amp::exp<Precision>(-state.x(0));
            }
            othererrors = othererrors || amp::abs<Precision>(state.x(0)-xprev)>(1+amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalEpsilon()))*stpmax;
            if( state.xupdated )
            {
                xprev = state.x(0);
            }
        }
        
        //
        // end
        //
        waserrors = referror || lin1error || lin2error || eqerror || converror || scerror || othererrors;
        if( !silent )
        {
            printf("TESTING LEVENBERG-MARQUARDT OPTIMIZATION\n");
            printf("REFERENCE PROBLEM:                        ");
            if( referror )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("1-D PROBLEM #1:                           ");
            if( lin1error )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("1-D PROBLEM #2:                           ");
            if( lin2error )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("LINEAR EQUATIONS:                         ");
            if( eqerror )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("CONVERGENCE PROPERTIES:                   ");
            if( converror )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("STATE FIELDS CONSISTENCY:                 ");
            if( scerror )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("OTHER PROPERTIES:                         ");
            if( othererrors )
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
    Asserts that State fields are consistent with RKind.
    Returns False otherwise.
    *************************************************************************/
    template<unsigned int Precision>
    bool rkindvsstatecheck(int rkind,
        const minlm::minlmstate<Precision>& state)
    {
        bool result;
        int nset;


        nset = 0;
        if( state.needf )
        {
            nset = nset+1;
        }
        if( state.needfg )
        {
            nset = nset+1;
        }
        if( state.needfij )
        {
            nset = nset+1;
        }
        if( state.needfgh )
        {
            nset = nset+1;
        }
        if( state.xupdated )
        {
            nset = nset+1;
        }
        if( nset!=1 )
        {
            result = false;
            return result;
        }
        if( rkind==0 && (state.needfg || state.needfgh) )
        {
            result = false;
            return result;
        }
        if( rkind==1 && state.needfgh )
        {
            result = false;
            return result;
        }
        if( rkind==2 && state.needfij )
        {
            result = false;
            return result;
        }
        result = true;
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testminlmunit_test_silent()
    {
        bool result;


        result = testminlm<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testminlmunit_test()
    {
        bool result;


        result = testminlm<Precision>(false);
        return result;
    }
} // namespace

#endif
