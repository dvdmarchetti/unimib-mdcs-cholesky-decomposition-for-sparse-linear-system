
#ifndef _testminlbfgsunit_h
#define _testminlbfgsunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "linmin.h"
#include "minlbfgs.h"
namespace testminlbfgsunit
{
    template<unsigned int Precision>
    bool testminlbfgs(bool silent);
    template<unsigned int Precision>
    void testfunc1(minlbfgs::minlbfgsstate<Precision>& state);
    template<unsigned int Precision>
    void testfunc2(minlbfgs::minlbfgsstate<Precision>& state);
    template<unsigned int Precision>
    void testfunc3(minlbfgs::minlbfgsstate<Precision>& state);
    template<unsigned int Precision>
    bool testminlbfgsunit_test_silent();
    template<unsigned int Precision>
    bool testminlbfgsunit_test();


    template<unsigned int Precision>
    bool testminlbfgs(bool silent)
    {
        bool result;
        bool waserrors;
        bool referror;
        bool nonconverror;
        bool eqerror;
        bool converror;
        bool crashtest;
        bool othererrors;
        int n;
        int m;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > xe;
        ap::template_1d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > xlast;
        int i;
        int j;
        amp::ampf<Precision> v;
        ap::template_2d_array< amp::ampf<Precision> > a;
        int maxits;
        minlbfgs::minlbfgsstate<Precision> state;
        minlbfgs::minlbfgsreport<Precision> rep;
        amp::ampf<Precision> fprev;
        amp::ampf<Precision> xprev;
        amp::ampf<Precision> stpmax;


        waserrors = false;
        
        //
        // Reference problem
        //
        x.setbounds(0, 2);
        n = 3;
        m = 2;
        x(0) = 100*amp::ampf<Precision>::getRandom()-50;
        x(1) = 100*amp::ampf<Precision>::getRandom()-50;
        x(2) = 100*amp::ampf<Precision>::getRandom()-50;
        minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
        minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
        while( minlbfgs::minlbfgsiteration<Precision>(state) )
        {
            state.f = amp::sqr<Precision>(state.x(0)-2)+amp::sqr<Precision>(state.x(1))+amp::sqr<Precision>(state.x(2)-state.x(0));
            state.g(0) = 2*(state.x(0)-2)+2*(state.x(0)-state.x(2));
            state.g(1) = 2*state.x(1);
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        minlbfgs::minlbfgsresults<Precision>(state, x, rep);
        referror = rep.terminationtype<=0 || amp::abs<Precision>(x(0)-2)>amp::ampf<Precision>("0.001") || amp::abs<Precision>(x(1))>amp::ampf<Precision>("0.001") || amp::abs<Precision>(x(2)-2)>amp::ampf<Precision>("0.001");
        
        //
        // nonconvex problems with hard relief: we start from point with very small
        // gradient, but we need ever smaller gradient in the next step due to
        // Wolfe conditions.
        //
        nonconverror = false;
        x.setlength(1);
        n = 1;
        m = 1;
        v = -100;
        while( v<amp::ampf<Precision>("0.1") )
        {
            x(0) = v;
            minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
            minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>("1.0E-9"), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
            while( minlbfgs::minlbfgsiteration<Precision>(state) )
            {
                state.f = amp::sqr<Precision>(state.x(0))/(1+amp::sqr<Precision>(state.x(0)));
                state.g(0) = (2*state.x(0)*(1+amp::sqr<Precision>(state.x(0)))-amp::sqr<Precision>(state.x(0))*2*state.x(0))/amp::sqr<Precision>(1+amp::sqr<Precision>(state.x(0)));
            }
            minlbfgs::minlbfgsresults<Precision>(state, x, rep);
            nonconverror = nonconverror || rep.terminationtype<=0 || amp::abs<Precision>(x(0))>amp::ampf<Precision>("0.001");
            v = v+amp::ampf<Precision>("0.1");
        }
        
        //
        // Linear equations
        //
        eqerror = false;
        for(n=1; n<=10; n++)
        {
            
            //
            // Prepare task
            //
            a.setbounds(0, n-1, 0, n-1);
            x.setbounds(0, n-1);
            xe.setbounds(0, n-1);
            b.setbounds(0, n-1);
            for(i=0; i<=n-1; i++)
            {
                xe(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                a(i,i) = a(i,i)+3*amp::sign<Precision>(a(i,i));
            }
            for(i=0; i<=n-1; i++)
            {
                v = amp::vdotproduct(a.getrow(i, 0, n-1), xe.getvector(0, n-1));
                b(i) = v;
            }
            
            //
            // Test different M
            //
            for(m=1; m<=n; m++)
            {
                
                //
                // Solve task
                //
                for(i=0; i<=n-1; i++)
                {
                    x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
                minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
                while( minlbfgs::minlbfgsiteration<Precision>(state) )
                {
                    state.f = 0;
                    for(i=0; i<=n-1; i++)
                    {
                        state.g(i) = 0;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = amp::vdotproduct(a.getrow(i, 0, n-1), state.x.getvector(0, n-1));
                        state.f = state.f+amp::sqr<Precision>(v-b(i));
                        for(j=0; j<=n-1; j++)
                        {
                            state.g(j) = state.g(j)+2*(v-b(i))*a(i,j);
                        }
                    }
                }
                minlbfgs::minlbfgsresults<Precision>(state, x, rep);
                eqerror = eqerror || rep.terminationtype<=0;
                for(i=0; i<=n-1; i++)
                {
                    eqerror = eqerror || amp::abs<Precision>(x(i)-xe(i))>amp::ampf<Precision>("0.001");
                }
            }
        }
        
        //
        // Testing convergence properties
        //
        converror = false;
        x.setbounds(0, 2);
        n = 3;
        m = 2;
        for(i=0; i<=2; i++)
        {
            x(i) = 6*amp::ampf<Precision>::getRandom()-3;
        }
        minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
        minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>("0.001"), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
        while( minlbfgs::minlbfgsiteration<Precision>(state) )
        {
            testfunc3<Precision>(state);
        }
        minlbfgs::minlbfgsresults<Precision>(state, x, rep);
        converror = converror || rep.terminationtype!=4;
        for(i=0; i<=2; i++)
        {
            x(i) = 6*amp::ampf<Precision>::getRandom()-3;
        }
        minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
        minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>("0.001"), amp::ampf<Precision>(0), 0);
        while( minlbfgs::minlbfgsiteration<Precision>(state) )
        {
            testfunc3<Precision>(state);
        }
        minlbfgs::minlbfgsresults<Precision>(state, x, rep);
        converror = converror || rep.terminationtype!=1;
        for(i=0; i<=2; i++)
        {
            x(i) = 6*amp::ampf<Precision>::getRandom()-3;
        }
        minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
        minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>("0.001"), 0);
        while( minlbfgs::minlbfgsiteration<Precision>(state) )
        {
            testfunc3<Precision>(state);
        }
        minlbfgs::minlbfgsresults<Precision>(state, x, rep);
        converror = converror || rep.terminationtype!=2;
        for(i=0; i<=2; i++)
        {
            x(i) = 2*amp::ampf<Precision>::getRandom()-1;
        }
        minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
        minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 10);
        while( minlbfgs::minlbfgsiteration<Precision>(state) )
        {
            testfunc3<Precision>(state);
        }
        minlbfgs::minlbfgsresults<Precision>(state, x, rep);
        converror = converror || rep.terminationtype!=5 || rep.iterationscount!=10;
        
        //
        // Crash test: too many iterations on a simple tasks
        // May fail when encounter zero step, underflow or something like that
        //
        crashtest = false;
        x.setbounds(0, 2);
        n = 3;
        m = 2;
        maxits = 10000;
        for(i=0; i<=2; i++)
        {
            x(i) = 6*amp::ampf<Precision>::getRandom()-3;
        }
        minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
        minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), maxits);
        while( minlbfgs::minlbfgsiteration<Precision>(state) )
        {
            state.f = amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(state.x(1))+amp::sqr<Precision>(state.x(2)-state.x(0));
            state.g(0) = 2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
            state.g(1) = 2*state.x(1);
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        minlbfgs::minlbfgsresults<Precision>(state, x, rep);
        crashtest = crashtest || rep.terminationtype<=0;
        
        //
        // Other properties:
        // 1. test reports (F should form monotone sequence)
        // 2. test maximum step
        //
        othererrors = false;
        n = 50;
        m = 2;
        x.setlength(n);
        xlast.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            x(i) = 1;
        }
        minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
        minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 100);
        minlbfgs::minlbfgssetxrep<Precision>(state, true);
        fprev = amp::ampf<Precision>::getAlgoPascalMaxNumber();
        while( minlbfgs::minlbfgsiteration<Precision>(state) )
        {
            if( state.needfg )
            {
                state.f = 0;
                for(i=0; i<=n-1; i++)
                {
                    state.f = state.f+amp::sqr<Precision>((1+i)*state.x(i));
                    state.g(i) = 2*(1+i)*state.x(i);
                }
            }
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
        minlbfgs::minlbfgsresults<Precision>(state, x, rep);
        for(i=0; i<=n-1; i++)
        {
            othererrors = othererrors || x(i)!=xlast(i);
        }
        n = 1;
        m = 1;
        x.setlength(n);
        x(0) = 100;
        stpmax = amp::ampf<Precision>("0.05")+amp::ampf<Precision>("0.05")*amp::ampf<Precision>::getRandom();
        minlbfgs::minlbfgscreate<Precision>(n, m, x, state);
        minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>("1.0E-9"), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
        minlbfgs::minlbfgssetstpmax<Precision>(state, stpmax);
        minlbfgs::minlbfgssetxrep<Precision>(state, true);
        xprev = x(0);
        while( minlbfgs::minlbfgsiteration<Precision>(state) )
        {
            if( state.needfg )
            {
                state.f = amp::exp<Precision>(state.x(0))+amp::exp<Precision>(-state.x(0));
                state.g(0) = amp::exp<Precision>(state.x(0))-amp::exp<Precision>(-state.x(0));
                othererrors = othererrors || amp::abs<Precision>(state.x(0)-xprev)>(1+amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalEpsilon()))*stpmax;
            }
            if( state.xupdated )
            {
                othererrors = othererrors || amp::abs<Precision>(state.x(0)-xprev)>(1+amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalEpsilon()))*stpmax;
                xprev = state.x(0);
            }
        }
        
        //
        // end
        //
        waserrors = referror || nonconverror || eqerror || converror || crashtest || othererrors;
        if( !silent )
        {
            printf("TESTING L-BFGS OPTIMIZATION\n");
            printf("REFERENCE PROBLEM:                        ");
            if( referror )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("NON-CONVEX PROBLEM:                       ");
            if( nonconverror )
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
            printf("CRASH TEST:                               ");
            if( crashtest )
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
    Calculate test function #1

    It may show very interesting behavior when optimized with 'x[0]>=ln(2)'
    constraint.
    *************************************************************************/
    template<unsigned int Precision>
    void testfunc1(minlbfgs::minlbfgsstate<Precision>& state)
    {
        if( state.x(0)<100 )
        {
            state.f = amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(state.x(1))+amp::sqr<Precision>(state.x(2)-state.x(0));
            state.g(0) = 2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
            state.g(1) = 2*state.x(1);
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        else
        {
            state.f = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(0) = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(1) = 0;
            state.g(2) = 0;
        }
    }


    /*************************************************************************
    Calculate test function #2

    Simple variation of #1, much more nonlinear, which makes unlikely premature
    convergence of algorithm .
    *************************************************************************/
    template<unsigned int Precision>
    void testfunc2(minlbfgs::minlbfgsstate<Precision>& state)
    {
        if( state.x(0)<100 )
        {
            state.f = amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(amp::sqr<Precision>(state.x(1)))+amp::sqr<Precision>(state.x(2)-state.x(0));
            state.g(0) = 2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
            state.g(1) = 4*state.x(1)*amp::sqr<Precision>(state.x(1));
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        else
        {
            state.f = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(0) = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(1) = 0;
            state.g(2) = 0;
        }
    }


    /*************************************************************************
    Calculate test function #3

    Simple variation of #1, much more nonlinear, with non-zero value at minimum.
    It achieve two goals:
    * makes unlikely premature convergence of algorithm .
    * solves some issues with EpsF stopping condition which arise when
      F(minimum) is zero

    *************************************************************************/
    template<unsigned int Precision>
    void testfunc3(minlbfgs::minlbfgsstate<Precision>& state)
    {
        amp::ampf<Precision> s;


        s = amp::ampf<Precision>("0.001");
        if( state.x(0)<100 )
        {
            state.f = amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(amp::sqr<Precision>(state.x(1))+s)+amp::sqr<Precision>(state.x(2)-state.x(0));
            state.g(0) = 2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
            state.g(1) = 2*(amp::sqr<Precision>(state.x(1))+s)*2*state.x(1);
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        else
        {
            state.f = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(0) = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(1) = 0;
            state.g(2) = 0;
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testminlbfgsunit_test_silent()
    {
        bool result;


        result = testminlbfgs<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testminlbfgsunit_test()
    {
        bool result;


        result = testminlbfgs<Precision>(false);
        return result;
    }
} // namespace

#endif
