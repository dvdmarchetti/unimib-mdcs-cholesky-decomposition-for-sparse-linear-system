
#ifndef _testodesolverunit_h
#define _testodesolverunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "odesolver.h"
namespace testodesolverunit
{
    template<unsigned int Precision>
    bool testodesolver(bool silent);
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void unsetrep(odesolver::odesolverreport<Precision>& rep);
    template<unsigned int Precision>
    bool testodesolverunit_test_silent();
    template<unsigned int Precision>
    bool testodesolverunit_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testodesolver(bool silent)
    {
        bool result;
        int passcount;
        bool curerrors;
        bool rkckerrors;
        bool waserrors;
        ap::template_1d_array< amp::ampf<Precision> > xtbl;
        ap::template_2d_array< amp::ampf<Precision> > ytbl;
        odesolver::odesolverreport<Precision> rep;
        ap::template_1d_array< amp::ampf<Precision> > xg;
        ap::template_1d_array< amp::ampf<Precision> > y;
        amp::ampf<Precision> h;
        amp::ampf<Precision> eps;
        int solver;
        int pass;
        int mynfev;
        amp::ampf<Precision> v;
        int n;
        int m;
        int m2;
        int i;
        amp::ampf<Precision> err;
        odesolver::odesolverstate<Precision> state;


        rkckerrors = false;
        waserrors = false;
        passcount = 10;
        
        //
        // simple test: just A*sin(x)+B*cos(x)
        //
        ap::ap_error::make_assertion(passcount>=2);
        for(pass=0; pass<=passcount-1; pass++)
        {
            for(solver=0; solver<=0; solver++)
            {
                
                //
                // prepare
                //
                h = amp::ampf<Precision>("1.0E-2");
                eps = amp::ampf<Precision>("1.0E-5");
                if( pass%2==0 )
                {
                    eps = -eps;
                }
                y.setlength(2);
                for(i=0; i<=1; i++)
                {
                    y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                m = 2+ap::randominteger(10);
                xg.setlength(m);
                xg(0) = (m-1)*amp::ampf<Precision>::getRandom();
                for(i=1; i<=m-1; i++)
                {
                    xg(i) = xg(i-1)+amp::ampf<Precision>::getRandom();
                }
                v = 2*amp::pi<Precision>()/(xg(m-1)-xg(0));
                amp::vmul(xg.getvector(0, m-1), v);
                if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
                {
                    amp::vmul(xg.getvector(0, m-1), -1);
                }
                mynfev = 0;
                
                //
                // choose solver
                //
                if( solver==0 )
                {
                    odesolver::odesolverrkck<Precision>(y, 2, xg, m, eps, h, state);
                }
                
                //
                // solve
                //
                while( odesolver::odesolveriteration<Precision>(state) )
                {
                    state.dy(0) = state.y(1);
                    state.dy(1) = -state.y(0);
                    mynfev = mynfev+1;
                }
                odesolver::odesolverresults<Precision>(state, m2, xtbl, ytbl, rep);
                
                //
                // check results
                //
                curerrors = false;
                if( rep.terminationtype<=0 )
                {
                    curerrors = true;
                }
                else
                {
                    curerrors = curerrors || m2!=m;
                    err = 0;
                    for(i=0; i<=m-1; i++)
                    {
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(ytbl(i,0)-(+y(0)*amp::cos<Precision>(xtbl(i)-xtbl(0))+y(1)*amp::sin<Precision>(xtbl(i)-xtbl(0)))));
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(ytbl(i,1)-(-y(0)*amp::sin<Precision>(xtbl(i)-xtbl(0))+y(1)*amp::cos<Precision>(xtbl(i)-xtbl(0)))));
                    }
                    curerrors = curerrors || err>10*amp::abs<Precision>(eps);
                    curerrors = curerrors || mynfev!=rep.nfev;
                }
                if( solver==0 )
                {
                    rkckerrors = rkckerrors || curerrors;
                }
            }
        }
        
        //
        // another test:
        //
        //     y(0)   = 0
        //     dy/dx  = f(x,y)
        //     f(x,y) = 0,   x<1
        //              x-1, x>=1
        //
        // with BOTH absolute and fractional tolerances.
        // Starting from zero will be real challenge for
        // fractional tolerance.
        //
        ap::ap_error::make_assertion(passcount>=2);
        for(pass=0; pass<=passcount-1; pass++)
        {
            h = amp::ampf<Precision>("1.0E-4");
            eps = amp::ampf<Precision>("1.0E-4");
            if( pass%2==0 )
            {
                eps = -eps;
            }
            y.setlength(1);
            y(0) = 0;
            m = 21;
            xg.setlength(m);
            for(i=0; i<=m-1; i++)
            {
                xg(i) = amp::ampf<Precision>(2*i)/(amp::ampf<Precision>(m-1));
            }
            mynfev = 0;
            odesolver::odesolverrkck<Precision>(y, 1, xg, m, eps, h, state);
            while( odesolver::odesolveriteration<Precision>(state) )
            {
                state.dy(0) = amp::maximum<Precision>(state.x-1, amp::ampf<Precision>(0));
                mynfev = mynfev+1;
            }
            odesolver::odesolverresults<Precision>(state, m2, xtbl, ytbl, rep);
            if( rep.terminationtype<=0 )
            {
                rkckerrors = true;
            }
            else
            {
                rkckerrors = rkckerrors || m2!=m;
                err = 0;
                for(i=0; i<=m-1; i++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(ytbl(i,0)-amp::sqr<Precision>(amp::maximum<Precision>(xg(i)-1, amp::ampf<Precision>(0)))/2));
                }
                rkckerrors = rkckerrors || err>amp::abs<Precision>(eps);
                rkckerrors = rkckerrors || mynfev!=rep.nfev;
            }
        }
        
        //
        // end
        //
        waserrors = rkckerrors;
        if( !silent )
        {
            printf("TESTING ODE SOLVER\n");
            printf("* RK CASH-KARP:                           ");
            if( rkckerrors )
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
        }
        result = !waserrors;
        return result;
    }


    /*************************************************************************
    Unsets real matrix
    *************************************************************************/
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        x.setlength(1, 1);
        x(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets real vector
    *************************************************************************/
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        x.setlength(1);
        x(0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets report
    *************************************************************************/
    template<unsigned int Precision>
    void unsetrep(odesolver::odesolverreport<Precision>& rep)
    {
        rep.nfev = 0;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testodesolverunit_test_silent()
    {
        bool result;


        result = testodesolver<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testodesolverunit_test()
    {
        bool result;


        result = testodesolver<Precision>(false);
        return result;
    }
} // namespace

#endif
