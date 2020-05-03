
#ifndef _testasa_h
#define _testasa_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "linmin.h"
#include "minasa.h"
namespace testasa
{
    template<unsigned int Precision>
    bool testminasa(bool silent);
    template<unsigned int Precision>
    void testfunc1(minasa::minasastate<Precision>& state);
    template<unsigned int Precision>
    void testfunc2(minasa::minasastate<Precision>& state);
    template<unsigned int Precision>
    void testfunc3(minasa::minasastate<Precision>& state);
    template<unsigned int Precision>
    void checkbounds(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& bndl,
        const ap::template_1d_array< amp::ampf<Precision> >& bndu,
        int n,
        bool& err);
    template<unsigned int Precision>
    amp::ampf<Precision> asaboundval(amp::ampf<Precision> x,
        amp::ampf<Precision> b1,
        amp::ampf<Precision> b2);
    template<unsigned int Precision>
    bool testasa_test_silent();
    template<unsigned int Precision>
    bool testasa_test();


    template<unsigned int Precision>
    bool testminasa(bool silent)
    {
        bool result;
        bool waserrors;
        bool referror;
        bool converror;
        bool othererrors;
        int n;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > xe;
        ap::template_1d_array< amp::ampf<Precision> > c;
        ap::template_1d_array< amp::ampf<Precision> > bndl;
        ap::template_1d_array< amp::ampf<Precision> > bndu;
        ap::template_1d_array< amp::ampf<Precision> > xlast;
        amp::ampf<Precision> fprev;
        amp::ampf<Precision> xprev;
        amp::ampf<Precision> stpmax;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s;
        amp::ampf<Precision> tol;
        int algotype;
        ap::template_2d_array< amp::ampf<Precision> > a;
        minasa::minasastate<Precision> state;
        minasa::minasareport<Precision> rep;


        waserrors = false;
        referror = false;
        converror = false;
        othererrors = false;
        
        //
        // Different algorithms
        //
        for(algotype=-1; algotype<=1; algotype++)
        {
            
            //
            // reference problem, simple convex optimization
            //
            for(n=1; n<=5; n++)
            {
                
                //
                // min(x'*diag(c)*x) on a random box
                //
                x.setlength(n);
                xe.setlength(n);
                c.setlength(n);
                bndl.setlength(n);
                bndu.setlength(n);
                for(i=0; i<=n-1; i++)
                {
                    c(i) = 1+amp::ampf<Precision>::getRandom();
                    xe(i) = 4*amp::ampf<Precision>::getRandom()-2;
                    bndl(i) = -amp::maximum<Precision>(amp::ampf<Precision>::getRandom(), amp::ampf<Precision>("0.2"));
                    bndu(i) = +amp::maximum<Precision>(amp::ampf<Precision>::getRandom(), amp::ampf<Precision>("0.2"));
                    x(i) = amp::ampf<Precision>("0.5")*(bndl(i)+bndu(i));
                }
                tol = amp::ampf<Precision>("0.001");
                minasa::minasacreate<Precision>(n, x, bndl, bndu, state);
                minasa::minasasetcond<Precision>(state, tol, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
                minasa::minasasetalgorithm<Precision>(state, algotype);
                while( minasa::minasaiteration<Precision>(state) )
                {
                    checkbounds<Precision>(state.x, bndl, bndu, n, othererrors);
                    state.f = 0;
                    for(i=0; i<=n-1; i++)
                    {
                        state.f = state.f+c(i)*amp::sqr<Precision>(state.x(i)-xe(i));
                        state.g(i) = 2*c(i)*(state.x(i)-xe(i));
                    }
                }
                minasa::minasaresults<Precision>(state, x, rep);
                referror = referror || rep.terminationtype<=0;
                for(i=0; i<=n-1; i++)
                {
                    referror = referror || amp::abs<Precision>(asaboundval<Precision>(xe(i), bndl(i), bndu(i))-x(i))>amp::ampf<Precision>("0.01");
                }
            }
            
            //
            // reference problem 2: non-convex optimization on [-2,2] x [1,2]
            //
            // A saddle function is minimized:
            // * stationary point [0,0] (non-feasible)
            // * constrained minimum [-2,2].
            // * starting point [+2,2]
            //
            // Path from start to end may be very complex, with multiple changes
            // in active constraints, so it is interesting task for our method.
            //
            // Scale parameter is used to make optimization more interesting
            // during GPA runs.
            //
            x.setlength(2);
            bndl.setlength(2);
            bndu.setlength(2);
            bndl(0) = -2;
            bndu(0) = 2;
            x(0) = 2;
            bndl(1) = 1;
            bndu(1) = 2;
            x(1) = 2;
            tol = amp::ampf<Precision>("0.001");
            s = amp::ampf<Precision>("0.01");
            minasa::minasacreate<Precision>(2, x, bndl, bndu, state);
            minasa::minasasetcond<Precision>(state, tol, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
            minasa::minasasetalgorithm<Precision>(state, algotype);
            while( minasa::minasaiteration<Precision>(state) )
            {
                checkbounds<Precision>(state.x, bndl, bndu, 2, othererrors);
                state.f = s*(amp::sqr<Precision>(state.x(0)+state.x(1))-amp::sqr<Precision>(state.x(0)-state.x(1)));
                state.g(0) = s*(2*(state.x(0)+state.x(1))-2*(state.x(0)-state.x(1)));
                state.g(1) = s*(2*(state.x(0)+state.x(1))+2*(state.x(0)-state.x(1)));
            }
            minasa::minasaresults<Precision>(state, x, rep);
            referror = referror || rep.terminationtype<=0 || amp::abs<Precision>(state.x(0)+2)>amp::ampf<Precision>("0.01") || amp::abs<Precision>(state.x(1)-2)>amp::ampf<Precision>("0.01");
            
            //
            // function #1 with 'x[0]>=ln(2)' constraint.
            // may show very interesting behavior.
            //
            x.setlength(3);
            bndl.setlength(3);
            bndu.setlength(3);
            n = 3;
            for(i=0; i<=2; i++)
            {
                bndl(i) = -10000;
                bndu(i) = +10000;
            }
            bndl(0) = amp::log<Precision>(amp::ampf<Precision>(2));
            for(i=0; i<=2; i++)
            {
                x(i) = 3*amp::ampf<Precision>::getRandom()+3;
            }
            minasa::minasacreate<Precision>(n, x, bndl, bndu, state);
            minasa::minasasetcond<Precision>(state, amp::ampf<Precision>("0.0000001"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
            minasa::minasasetalgorithm<Precision>(state, algotype);
            while( minasa::minasaiteration<Precision>(state) )
            {
                checkbounds<Precision>(state.x, bndl, bndu, n, othererrors);
                testfunc1<Precision>(state);
            }
            minasa::minasaresults<Precision>(state, x, rep);
            referror = referror || rep.terminationtype<=0;
            referror = referror || amp::abs<Precision>(x(0)-amp::log<Precision>(amp::ampf<Precision>(2)))>amp::ampf<Precision>("0.05");
            referror = referror || amp::abs<Precision>(x(1))>amp::ampf<Precision>("0.05");
            referror = referror || amp::abs<Precision>(x(2)-amp::log<Precision>(amp::ampf<Precision>(2)))>amp::ampf<Precision>("0.05");
            
            //
            // Testing convergence properties
            //
            x.setlength(3);
            bndl.setlength(3);
            bndu.setlength(3);
            n = 3;
            for(i=0; i<=2; i++)
            {
                bndl(i) = -10000;
                bndu(i) = +10000;
            }
            bndl(0) = amp::log<Precision>(amp::ampf<Precision>(2));
            for(i=0; i<=2; i++)
            {
                x(i) = 3*amp::ampf<Precision>::getRandom()+3;
            }
            minasa::minasacreate<Precision>(n, x, bndl, bndu, state);
            minasa::minasasetcond<Precision>(state, amp::ampf<Precision>("0.001"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
            minasa::minasasetalgorithm<Precision>(state, algotype);
            while( minasa::minasaiteration<Precision>(state) )
            {
                checkbounds<Precision>(state.x, bndl, bndu, n, othererrors);
                testfunc3<Precision>(state);
            }
            minasa::minasaresults<Precision>(state, x, rep);
            converror = converror || rep.terminationtype!=4;
            for(i=0; i<=2; i++)
            {
                x(i) = 3*amp::ampf<Precision>::getRandom()+3;
            }
            minasa::minasacreate<Precision>(n, x, bndl, bndu, state);
            minasa::minasasetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.001"), amp::ampf<Precision>("0.0"), 0);
            minasa::minasasetalgorithm<Precision>(state, algotype);
            while( minasa::minasaiteration<Precision>(state) )
            {
                checkbounds<Precision>(state.x, bndl, bndu, n, othererrors);
                testfunc3<Precision>(state);
            }
            minasa::minasaresults<Precision>(state, x, rep);
            converror = converror || rep.terminationtype!=1;
            for(i=0; i<=2; i++)
            {
                x(i) = 3*amp::ampf<Precision>::getRandom()+3;
            }
            minasa::minasacreate<Precision>(n, x, bndl, bndu, state);
            minasa::minasasetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.001"), 0);
            minasa::minasasetalgorithm<Precision>(state, algotype);
            while( minasa::minasaiteration<Precision>(state) )
            {
                checkbounds<Precision>(state.x, bndl, bndu, n, othererrors);
                testfunc3<Precision>(state);
            }
            minasa::minasaresults<Precision>(state, x, rep);
            converror = converror || rep.terminationtype!=2;
            for(i=0; i<=2; i++)
            {
                x(i) = 3*amp::ampf<Precision>::getRandom()+3;
            }
            minasa::minasacreate<Precision>(n, x, bndl, bndu, state);
            minasa::minasasetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 3);
            minasa::minasasetalgorithm<Precision>(state, algotype);
            while( minasa::minasaiteration<Precision>(state) )
            {
                checkbounds<Precision>(state.x, bndl, bndu, n, othererrors);
                testfunc3<Precision>(state);
            }
            minasa::minasaresults<Precision>(state, x, rep);
            converror = converror || !(rep.terminationtype==5 && rep.iterationscount==3 || rep.terminationtype==7);
            
            //
            // Other properties
            //
            //
            // Other properties:
            // 1. test reports (F should form monotone sequence)
            // 2. test maximum step
            //
            n = 50;
            x.setlength(n);
            xlast.setlength(n);
            bndl.setlength(n);
            bndu.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                x(i) = 1;
                xlast(i) = amp::ampf<Precision>::getRandom();
                bndl(i) = -100000;
                bndu(i) = +100000;
            }
            minasa::minasacreate<Precision>(n, x, bndl, bndu, state);
            minasa::minasasetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 100);
            minasa::minasasetxrep<Precision>(state, true);
            fprev = amp::ampf<Precision>::getAlgoPascalMaxNumber();
            while( minasa::minasaiteration<Precision>(state) )
            {
                checkbounds<Precision>(state.x, bndl, bndu, n, othererrors);
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
            minasa::minasaresults<Precision>(state, x, rep);
            for(i=0; i<=n-1; i++)
            {
                othererrors = othererrors || x(i)!=xlast(i);
            }
            n = 1;
            x.setlength(n);
            bndl.setlength(n);
            bndu.setlength(n);
            x(0) = 100;
            bndl(0) = -1000000;
            bndu(0) = +1000000;
            stpmax = amp::ampf<Precision>("0.05")+amp::ampf<Precision>("0.05")*amp::ampf<Precision>::getRandom();
            minasa::minasacreate<Precision>(n, x, bndl, bndu, state);
            minasa::minasasetcond<Precision>(state, amp::ampf<Precision>("1.0E-9"), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
            minasa::minasasetstpmax<Precision>(state, stpmax);
            minasa::minasasetxrep<Precision>(state, true);
            xprev = x(0);
            while( minasa::minasaiteration<Precision>(state) )
            {
                checkbounds<Precision>(state.x, bndl, bndu, n, othererrors);
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
        }
        
        //
        // end
        //
        waserrors = referror || converror || othererrors;
        if( !silent )
        {
            printf("TESTING ASA OPTIMIZATION\n");
            printf("REFERENCE PROBLEMS:                       ");
            if( referror )
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
    void testfunc1(minasa::minasastate<Precision>& state)
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
    void testfunc2(minasa::minasastate<Precision>& state)
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
    void testfunc3(minasa::minasastate<Precision>& state)
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
    Checks that X is bounded with respect to BndL/BndU.

    If it is not, True is assigned to the Err variable (which is not changed
    otherwise).
    *************************************************************************/
    template<unsigned int Precision>
    void checkbounds(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& bndl,
        const ap::template_1d_array< amp::ampf<Precision> >& bndu,
        int n,
        bool& err)
    {
        int i;


        for(i=0; i<=n-1; i++)
        {
            if( x(i)<bndl(i) || x(i)>bndu(i) )
            {
                err = true;
            }
        }
    }


    /*************************************************************************
    'bound' value: map X to [B1,B2]
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> asaboundval(amp::ampf<Precision> x,
        amp::ampf<Precision> b1,
        amp::ampf<Precision> b2)
    {
        amp::ampf<Precision> result;


        if( x<=b1 )
        {
            result = b1;
            return result;
        }
        if( x>=b2 )
        {
            result = b2;
            return result;
        }
        result = x;
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testasa_test_silent()
    {
        bool result;


        result = testminasa<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testasa_test()
    {
        bool result;


        result = testminasa<Precision>(false);
        return result;
    }
} // namespace

#endif
