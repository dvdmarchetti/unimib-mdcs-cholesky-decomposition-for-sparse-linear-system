
#ifndef _testpolintunit_h
#define _testpolintunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "tsort.h"
#include "ratinterpolation.h"
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
#include "lsfit.h"
#include "ratint.h"
#include "apserv.h"
#include "polint.h"
namespace testpolintunit
{
    template<unsigned int Precision>
    bool testpolint(bool silent);
    template<unsigned int Precision>
    amp::ampf<Precision> internalpolint(const ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> > f,
        int n,
        amp::ampf<Precision> t);
    template<unsigned int Precision>
    void brcunset(ratint::barycentricinterpolant<Precision>& b);
    template<unsigned int Precision>
    bool testpolintunit_test_silent();
    template<unsigned int Precision>
    bool testpolintunit_test();


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testpolint(bool silent)
    {
        bool result;
        bool waserrors;
        bool interrors;
        bool fiterrors;
        amp::ampf<Precision> threshold;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > x2;
        ap::template_1d_array< amp::ampf<Precision> > y2;
        ap::template_1d_array< amp::ampf<Precision> > w2;
        ap::template_1d_array< amp::ampf<Precision> > xfull;
        ap::template_1d_array< amp::ampf<Precision> > yfull;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        amp::ampf<Precision> t;
        int i;
        int k;
        ap::template_1d_array< amp::ampf<Precision> > xc;
        ap::template_1d_array< amp::ampf<Precision> > yc;
        ap::template_1d_array< int > dc;
        int info;
        int info2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> v0;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> s;
        amp::ampf<Precision> xmin;
        amp::ampf<Precision> xmax;
        amp::ampf<Precision> refrms;
        amp::ampf<Precision> refavg;
        amp::ampf<Precision> refavgrel;
        amp::ampf<Precision> refmax;
        ratint::barycentricinterpolant<Precision> p;
        ratint::barycentricinterpolant<Precision> p1;
        ratint::barycentricinterpolant<Precision> p2;
        polint::polynomialfitreport<Precision> rep;
        polint::polynomialfitreport<Precision> rep2;
        int n;
        int m;
        int maxn;
        int pass;
        int passcount;


        waserrors = false;
        interrors = false;
        fiterrors = false;
        maxn = 5;
        passcount = 20;
        threshold = amp::ampf<Precision>("1.0E8")*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Test equidistant interpolation
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                
                //
                // prepare task:
                // * equidistant points
                // * random Y
                // * T in [A,B] or near (within 10% of its width)
                //
                do
                {
                    a = 2*amp::ampf<Precision>::getRandom()-1;
                    b = 2*amp::ampf<Precision>::getRandom()-1;
                }
                while( amp::abs<Precision>(a-b)<=amp::ampf<Precision>("0.2") );
                t = a+(amp::ampf<Precision>("1.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1"))*(b-a);
                apserv::taskgenint1dequidist<Precision>(a, b, n, x, y);
                
                //
                // test "fast" equidistant interpolation (no barycentric model)
                //
                interrors = interrors || amp::abs<Precision>(polint::polynomialcalceqdist<Precision>(a, b, y, n, t)-internalpolint<Precision>(x, y, n, t))>threshold;
                
                //
                // test "slow" equidistant interpolation (create barycentric model)
                //
                brcunset<Precision>(p);
                polint::polynomialbuild<Precision>(x, y, n, p);
                interrors = interrors || amp::abs<Precision>(ratint::barycentriccalc<Precision>(p, t)-internalpolint<Precision>(x, y, n, t))>threshold;
                
                //
                // test "fast" interpolation (create "fast" barycentric model)
                //
                brcunset<Precision>(p);
                polint::polynomialbuildeqdist<Precision>(a, b, y, n, p);
                interrors = interrors || amp::abs<Precision>(ratint::barycentriccalc<Precision>(p, t)-internalpolint<Precision>(x, y, n, t))>threshold;
            }
        }
        
        //
        // Test Chebyshev-1 interpolation
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                
                //
                // prepare task:
                // * equidistant points
                // * random Y
                // * T in [A,B] or near (within 10% of its width)
                //
                do
                {
                    a = 2*amp::ampf<Precision>::getRandom()-1;
                    b = 2*amp::ampf<Precision>::getRandom()-1;
                }
                while( amp::abs<Precision>(a-b)<=amp::ampf<Precision>("0.2") );
                t = a+(amp::ampf<Precision>("1.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1"))*(b-a);
                apserv::taskgenint1dcheb1<Precision>(a, b, n, x, y);
                
                //
                // test "fast" interpolation (no barycentric model)
                //
                interrors = interrors || amp::abs<Precision>(polint::polynomialcalccheb1<Precision>(a, b, y, n, t)-internalpolint<Precision>(x, y, n, t))>threshold;
                
                //
                // test "slow" interpolation (create barycentric model)
                //
                brcunset<Precision>(p);
                polint::polynomialbuild<Precision>(x, y, n, p);
                interrors = interrors || amp::abs<Precision>(ratint::barycentriccalc<Precision>(p, t)-internalpolint<Precision>(x, y, n, t))>threshold;
                
                //
                // test "fast" interpolation (create "fast" barycentric model)
                //
                brcunset<Precision>(p);
                polint::polynomialbuildcheb1<Precision>(a, b, y, n, p);
                interrors = interrors || amp::abs<Precision>(ratint::barycentriccalc<Precision>(p, t)-internalpolint<Precision>(x, y, n, t))>threshold;
            }
        }
        
        //
        // Test Chebyshev-2 interpolation
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                
                //
                // prepare task:
                // * equidistant points
                // * random Y
                // * T in [A,B] or near (within 10% of its width)
                //
                do
                {
                    a = 2*amp::ampf<Precision>::getRandom()-1;
                    b = 2*amp::ampf<Precision>::getRandom()-1;
                }
                while( amp::abs<Precision>(a-b)<=amp::ampf<Precision>("0.2") );
                t = a+(amp::ampf<Precision>("1.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1"))*(b-a);
                apserv::taskgenint1dcheb2<Precision>(a, b, n, x, y);
                
                //
                // test "fast" interpolation (no barycentric model)
                //
                interrors = interrors || amp::abs<Precision>(polint::polynomialcalccheb2<Precision>(a, b, y, n, t)-internalpolint<Precision>(x, y, n, t))>threshold;
                
                //
                // test "slow" interpolation (create barycentric model)
                //
                brcunset<Precision>(p);
                polint::polynomialbuild<Precision>(x, y, n, p);
                interrors = interrors || amp::abs<Precision>(ratint::barycentriccalc<Precision>(p, t)-internalpolint<Precision>(x, y, n, t))>threshold;
                
                //
                // test "fast" interpolation (create "fast" barycentric model)
                //
                brcunset<Precision>(p);
                polint::polynomialbuildcheb2<Precision>(a, b, y, n, p);
                interrors = interrors || amp::abs<Precision>(ratint::barycentriccalc<Precision>(p, t)-internalpolint<Precision>(x, y, n, t))>threshold;
            }
        }
        
        //
        // crash-test: ability to solve tasks which will overflow/underflow
        // weights with straightforward implementation
        //
        for(n=1; n<=20; n++)
        {
            a = -amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
            b = +amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
            apserv::taskgenint1dequidist<Precision>(a, b, n, x, y);
            polint::polynomialbuild<Precision>(x, y, n, p);
            for(i=0; i<=n-1; i++)
            {
                interrors = interrors || p.w(i)==0;
            }
        }
        
        //
        // Test rational fitting:
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                
                //
                // N=M+K fitting (i.e. interpolation)
                //
                for(k=0; k<=n-1; k++)
                {
                    apserv::taskgenint1d<Precision>(amp::ampf<Precision>(-1), amp::ampf<Precision>(1), n, xfull, yfull);
                    x.setlength(n-k);
                    y.setlength(n-k);
                    w.setlength(n-k);
                    if( k>0 )
                    {
                        xc.setlength(k);
                        yc.setlength(k);
                        dc.setlength(k);
                    }
                    for(i=0; i<=n-k-1; i++)
                    {
                        x(i) = xfull(i);
                        y(i) = yfull(i);
                        w(i) = 1+amp::ampf<Precision>::getRandom();
                    }
                    for(i=0; i<=k-1; i++)
                    {
                        xc(i) = xfull(n-k+i);
                        yc(i) = yfull(n-k+i);
                        dc(i) = 0;
                    }
                    polint::polynomialfitwc<Precision>(x, y, w, n-k, xc, yc, dc, k, n, info, p1, rep);
                    if( info<=0 )
                    {
                        fiterrors = true;
                    }
                    else
                    {
                        for(i=0; i<=n-k-1; i++)
                        {
                            fiterrors = fiterrors || amp::abs<Precision>(ratint::barycentriccalc<Precision>(p1, x(i))-y(i))>threshold;
                        }
                        for(i=0; i<=k-1; i++)
                        {
                            fiterrors = fiterrors || amp::abs<Precision>(ratint::barycentriccalc<Precision>(p1, xc(i))-yc(i))>threshold;
                        }
                    }
                }
                
                //
                // Testing constraints on derivatives.
                // Special tasks which will always have solution:
                // 1. P(0)=YC[0]
                // 2. P(0)=YC[0], P'(0)=YC[1]
                //
                if( n>1 )
                {
                    for(m=3; m<=5; m++)
                    {
                        for(k=1; k<=2; k++)
                        {
                            apserv::taskgenint1d<Precision>(amp::ampf<Precision>(-1), amp::ampf<Precision>(1), n, x, y);
                            w.setlength(n);
                            xc.setlength(2);
                            yc.setlength(2);
                            dc.setlength(2);
                            for(i=0; i<=n-1; i++)
                            {
                                w(i) = 1+amp::ampf<Precision>::getRandom();
                            }
                            xc(0) = 0;
                            yc(0) = 2*amp::ampf<Precision>::getRandom()-1;
                            dc(0) = 0;
                            xc(1) = 0;
                            yc(1) = 2*amp::ampf<Precision>::getRandom()-1;
                            dc(1) = 1;
                            polint::polynomialfitwc<Precision>(x, y, w, n, xc, yc, dc, k, m, info, p1, rep);
                            if( info<=0 )
                            {
                                fiterrors = true;
                            }
                            else
                            {
                                ratint::barycentricdiff1<Precision>(p1, amp::ampf<Precision>("0.0"), v0, v1);
                                fiterrors = fiterrors || amp::abs<Precision>(v0-yc(0))>threshold;
                                if( k==2 )
                                {
                                    fiterrors = fiterrors || amp::abs<Precision>(v1-yc(1))>threshold;
                                }
                            }
                        }
                    }
                }
            }
        }
        for(m=2; m<=8; m++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // General fitting
                //
                // interpolating function through M nodes should have
                // greater RMS error than fitting it through the same M nodes
                //
                n = 100;
                x2.setlength(n);
                y2.setlength(n);
                w2.setlength(n);
                xmin = 0;
                xmax = 2*amp::pi<Precision>();
                for(i=0; i<=n-1; i++)
                {
                    x2(i) = 2*amp::pi<Precision>()*amp::ampf<Precision>::getRandom();
                    y2(i) = amp::sin<Precision>(x2(i));
                    w2(i) = 1;
                }
                x.setlength(m);
                y.setlength(m);
                for(i=0; i<=m-1; i++)
                {
                    x(i) = xmin+(xmax-xmin)*i/(m-1);
                    y(i) = amp::sin<Precision>(x(i));
                }
                polint::polynomialbuild<Precision>(x, y, m, p1);
                polint::polynomialfitwc<Precision>(x2, y2, w2, n, xc, yc, dc, 0, m, info, p2, rep);
                if( info<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    
                    //
                    // calculate P1 (interpolant) RMS error, compare with P2 error
                    //
                    v1 = 0;
                    v2 = 0;
                    for(i=0; i<=n-1; i++)
                    {
                        v1 = v1+amp::sqr<Precision>(ratint::barycentriccalc<Precision>(p1, x2(i))-y2(i));
                        v2 = v2+amp::sqr<Precision>(ratint::barycentriccalc<Precision>(p2, x2(i))-y2(i));
                    }
                    v1 = amp::sqrt<Precision>(v1/n);
                    v2 = amp::sqrt<Precision>(v2/n);
                    fiterrors = fiterrors || v2>v1;
                    fiterrors = fiterrors || amp::abs<Precision>(v2-rep.rmserror)>threshold;
                }
                
                //
                // compare weighted and non-weighted
                //
                n = 20;
                x.setlength(n);
                y.setlength(n);
                w.setlength(n);
                for(i=0; i<=n-1; i++)
                {
                    x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    w(i) = 1;
                }
                polint::polynomialfitwc<Precision>(x, y, w, n, xc, yc, dc, 0, m, info, p1, rep);
                polint::polynomialfit<Precision>(x, y, n, m, info2, p2, rep2);
                if( info<=0 || info2<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    
                    //
                    // calculate P1 (interpolant), compare with P2 error
                    // compare RMS errors
                    //
                    t = 2*amp::ampf<Precision>::getRandom()-1;
                    v1 = ratint::barycentriccalc<Precision>(p1, t);
                    v2 = ratint::barycentriccalc<Precision>(p2, t);
                    fiterrors = fiterrors || v2!=v1;
                    fiterrors = fiterrors || rep.rmserror!=rep2.rmserror;
                    fiterrors = fiterrors || rep.avgerror!=rep2.avgerror;
                    fiterrors = fiterrors || rep.avgrelerror!=rep2.avgrelerror;
                    fiterrors = fiterrors || rep.maxerror!=rep2.maxerror;
                }
            }
        }
        for(pass=1; pass<=passcount; pass++)
        {
            ap::ap_error::make_assertion(passcount>=2);
            
            //
            // solve simple task (all X[] are the same, Y[] are specially
            // calculated to ensure simple form of all types of errors)
            // and check correctness of the errors calculated by subroutines
            //
            // First pass is done with zero Y[], other passes - with random Y[].
            // It should test both ability to correctly calculate errors and
            // ability to not fail while working with zeros :)
            //
            n = 4;
            if( pass==1 )
            {
                v1 = 0;
                v2 = 0;
                v = 0;
            }
            else
            {
                v1 = amp::ampf<Precision>::getRandom();
                v2 = amp::ampf<Precision>::getRandom();
                v = 1+amp::ampf<Precision>::getRandom();
            }
            x.setlength(4);
            y.setlength(4);
            w.setlength(4);
            x(0) = 0;
            y(0) = v-v2;
            w(0) = 1;
            x(1) = 0;
            y(1) = v-v1;
            w(1) = 1;
            x(2) = 0;
            y(2) = v+v1;
            w(2) = 1;
            x(3) = 0;
            y(3) = v+v2;
            w(3) = 1;
            refrms = amp::sqrt<Precision>((amp::sqr<Precision>(v1)+amp::sqr<Precision>(v2))/2);
            refavg = (amp::abs<Precision>(v1)+amp::abs<Precision>(v2))/2;
            if( pass==1 )
            {
                refavgrel = 0;
            }
            else
            {
                refavgrel = amp::ampf<Precision>("0.25")*(amp::abs<Precision>(v2)/amp::abs<Precision>(v-v2)+amp::abs<Precision>(v1)/amp::abs<Precision>(v-v1)+amp::abs<Precision>(v1)/amp::abs<Precision>(v+v1)+amp::abs<Precision>(v2)/amp::abs<Precision>(v+v2));
            }
            refmax = amp::maximum<Precision>(v1, v2);
            
            //
            // Test errors correctness
            //
            polint::polynomialfit<Precision>(x, y, 4, 1, info, p, rep);
            if( info<=0 )
            {
                fiterrors = true;
            }
            else
            {
                s = ratint::barycentriccalc<Precision>(p, amp::ampf<Precision>(0));
                fiterrors = fiterrors || amp::abs<Precision>(s-v)>threshold;
                fiterrors = fiterrors || amp::abs<Precision>(rep.rmserror-refrms)>threshold;
                fiterrors = fiterrors || amp::abs<Precision>(rep.avgerror-refavg)>threshold;
                fiterrors = fiterrors || amp::abs<Precision>(rep.avgrelerror-refavgrel)>threshold;
                fiterrors = fiterrors || amp::abs<Precision>(rep.maxerror-refmax)>threshold;
            }
        }
        
        //
        // report
        //
        waserrors = interrors || fiterrors;
        if( !silent )
        {
            printf("TESTING POLYNOMIAL INTERPOLATION AND FITTING\n");
            
            //
            // Normal tests
            //
            printf("INTERPOLATION TEST:                      ");
            if( interrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("FITTING TEST:                            ");
            if( fiterrors )
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


    template<unsigned int Precision>
    amp::ampf<Precision> internalpolint(const ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> > f,
        int n,
        amp::ampf<Precision> t)
    {
        amp::ampf<Precision> result;
        int i;
        int j;


        n = n-1;
        for(j=0; j<=n-1; j++)
        {
            for(i=j+1; i<=n; i++)
            {
                f(i) = ((t-x(j))*f(i)-(t-x(i))*f(j))/(x(i)-x(j));
            }
        }
        result = f(n);
        return result;
    }


    template<unsigned int Precision>
    void brcunset(ratint::barycentricinterpolant<Precision>& b)
    {
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > w;


        x.setlength(1);
        y.setlength(1);
        w.setlength(1);
        x(0) = 0;
        y(0) = 0;
        w(0) = 1;
        ratint::barycentricbuildxyw<Precision>(x, y, w, 1, b);
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testpolintunit_test_silent()
    {
        bool result;


        result = testpolint<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testpolintunit_test()
    {
        bool result;


        result = testpolint<Precision>(false);
        return result;
    }
} // namespace

#endif
