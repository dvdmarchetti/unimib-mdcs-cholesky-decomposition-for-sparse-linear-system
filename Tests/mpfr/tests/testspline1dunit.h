
#ifndef _testspline1dunit_h
#define _testspline1dunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "spline3.h"
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
#include "apserv.h"
#include "spline1d.h"
namespace testspline1dunit
{
    template<unsigned int Precision>
    bool testsplineinterpolation(bool silent);
    template<unsigned int Precision>
    void lconst(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const spline1d::spline1dinterpolant<Precision>& c,
        amp::ampf<Precision> lstep,
        amp::ampf<Precision>& l0,
        amp::ampf<Precision>& l1,
        amp::ampf<Precision>& l2);
    template<unsigned int Precision>
    bool testunpack(const spline1d::spline1dinterpolant<Precision>& c,
        const ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void unsetspline1d(spline1d::spline1dinterpolant<Precision>& c);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    bool is1dsolution(int n,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        amp::ampf<Precision> c);
    template<unsigned int Precision>
    bool testspline1dunit_test_silent();
    template<unsigned int Precision>
    bool testspline1dunit_test();


    template<unsigned int Precision>
    bool testsplineinterpolation(bool silent)
    {
        bool result;
        bool waserrors;
        bool crserrors;
        bool cserrors;
        bool hserrors;
        bool aserrors;
        bool lserrors;
        bool dserrors;
        bool uperrors;
        bool cperrors;
        bool lterrors;
        bool ierrors;
        bool fiterrors;
        amp::ampf<Precision> nonstrictthreshold;
        amp::ampf<Precision> threshold;
        int passcount;
        amp::ampf<Precision> lstep;
        amp::ampf<Precision> h;
        int maxn;
        int bltype;
        int brtype;
        bool periodiccond;
        int n;
        int m;
        int i;
        int k;
        int pass;
        int stype;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > yp;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > w2;
        ap::template_1d_array< amp::ampf<Precision> > y2;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > xc;
        ap::template_1d_array< amp::ampf<Precision> > yc;
        ap::template_1d_array< int > dc;
        spline1d::spline1dinterpolant<Precision> c;
        spline1d::spline1dinterpolant<Precision> c2;
        int info;
        int info1;
        int info2;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        amp::ampf<Precision> bl;
        amp::ampf<Precision> br;
        amp::ampf<Precision> t;
        amp::ampf<Precision> sa;
        amp::ampf<Precision> sb;
        amp::ampf<Precision> v;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> l10;
        amp::ampf<Precision> l11;
        amp::ampf<Precision> l12;
        amp::ampf<Precision> l20;
        amp::ampf<Precision> l21;
        amp::ampf<Precision> l22;
        amp::ampf<Precision> p0;
        amp::ampf<Precision> p1;
        amp::ampf<Precision> p2;
        amp::ampf<Precision> s;
        amp::ampf<Precision> ds;
        amp::ampf<Precision> d2s;
        amp::ampf<Precision> s2;
        amp::ampf<Precision> ds2;
        amp::ampf<Precision> d2s2;
        amp::ampf<Precision> vl;
        amp::ampf<Precision> vm;
        amp::ampf<Precision> vr;
        amp::ampf<Precision> err;
        amp::ampf<Precision> tension;
        amp::ampf<Precision> intab;
        spline1d::spline1dfitreport<Precision> rep;
        spline1d::spline1dfitreport<Precision> rep2;
        amp::ampf<Precision> refrms;
        amp::ampf<Precision> refavg;
        amp::ampf<Precision> refavgrel;
        amp::ampf<Precision> refmax;


        waserrors = false;
        passcount = 20;
        lstep = amp::ampf<Precision>("0.005");
        h = amp::ampf<Precision>("0.00001");
        maxn = 10;
        threshold = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        nonstrictthreshold = amp::ampf<Precision>("0.00001");
        lserrors = false;
        cserrors = false;
        crserrors = false;
        hserrors = false;
        aserrors = false;
        dserrors = false;
        cperrors = false;
        uperrors = false;
        lterrors = false;
        ierrors = false;
        fiterrors = false;
        
        //
        // General test: linear, cubic, Hermite, Akima
        //
        for(n=2; n<=maxn; n++)
        {
            x.setbounds(0, n-1);
            y.setbounds(0, n-1);
            yp.setbounds(0, n-1);
            d.setbounds(0, n-1);
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // Prepare task:
                // * X contains abscissas from [A,B]
                // * Y contains function values
                // * YP contains periodic function values
                //
                a = -1-amp::ampf<Precision>::getRandom();
                b = +1+amp::ampf<Precision>::getRandom();
                bl = 2*amp::ampf<Precision>::getRandom()-1;
                br = 2*amp::ampf<Precision>::getRandom()-1;
                for(i=0; i<=n-1; i++)
                {
                    x(i) = amp::ampf<Precision>("0.5")*(b+a)+amp::ampf<Precision>("0.5")*(b-a)*amp::cos<Precision>(amp::pi<Precision>()*(2*i+1)/(2*n));
                    if( i==0 )
                    {
                        x(i) = a;
                    }
                    if( i==n-1 )
                    {
                        x(i) = b;
                    }
                    y(i) = amp::cos<Precision>(amp::ampf<Precision>("1.3")*amp::pi<Precision>()*x(i)+amp::ampf<Precision>("0.4"));
                    yp(i) = y(i);
                    d(i) = -amp::ampf<Precision>("1.3")*amp::pi<Precision>()*amp::sin<Precision>(amp::ampf<Precision>("1.3")*amp::pi<Precision>()*x(i)+amp::ampf<Precision>("0.4"));
                }
                yp(n-1) = yp(0);
                for(i=0; i<=n-1; i++)
                {
                    k = ap::randominteger(n);
                    if( k!=i )
                    {
                        t = x(i);
                        x(i) = x(k);
                        x(k) = t;
                        t = y(i);
                        y(i) = y(k);
                        y(k) = t;
                        t = yp(i);
                        yp(i) = yp(k);
                        yp(k) = t;
                        t = d(i);
                        d(i) = d(k);
                        d(k) = t;
                    }
                }
                
                //
                // Build linear spline
                // Test for general interpolation scheme properties:
                // * values at nodes
                // * continuous function
                // Test for specific properties is implemented below.
                //
                spline1d::spline1dbuildlinear<Precision>(x, y, n, c);
                err = 0;
                for(i=0; i<=n-1; i++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(y(i)-spline1d::spline1dcalc<Precision>(c, x(i))));
                }
                lserrors = lserrors || err>threshold;
                lconst<Precision>(a, b, c, lstep, l10, l11, l12);
                lconst<Precision>(a, b, c, lstep/3, l20, l21, l22);
                lserrors = lserrors || l20/l10>amp::ampf<Precision>("1.2");
                
                //
                // Build cubic spline.
                // Test for interpolation scheme properties:
                // * values at nodes
                // * boundary conditions
                // * continuous function
                // * continuous first derivative
                // * continuous second derivative
                // * periodicity properties
                //
                for(bltype=-1; bltype<=2; bltype++)
                {
                    for(brtype=-1; brtype<=2; brtype++)
                    {
                        
                        //
                        // skip meaningless combination of boundary conditions
                        // (one condition is periodic, another is not)
                        //
                        periodiccond = bltype==-1 || brtype==-1;
                        if( periodiccond && bltype!=brtype )
                        {
                            continue;
                        }
                        
                        //
                        // build
                        //
                        if( periodiccond )
                        {
                            spline1d::spline1dbuildcubic<Precision>(x, yp, n, bltype, bl, brtype, br, c);
                        }
                        else
                        {
                            spline1d::spline1dbuildcubic<Precision>(x, y, n, bltype, bl, brtype, br, c);
                        }
                        
                        //
                        // interpolation properties
                        //
                        err = 0;
                        if( periodiccond )
                        {
                            
                            //
                            // * check values at nodes; spline is periodic so
                            //   we add random number of periods to nodes
                            // * we also test for periodicity of derivatives
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                v = x(i);
                                vm = v+(b-a)*(ap::randominteger(5)-2);
                                t = yp(i)-spline1d::spline1dcalc<Precision>(c, vm);
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(t));
                                spline1d::spline1ddiff<Precision>(c, v, s, ds, d2s);
                                spline1d::spline1ddiff<Precision>(c, vm, s2, ds2, d2s2);
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(s-s2));
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(ds-ds2));
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(d2s-d2s2));
                            }
                            
                            //
                            // periodicity between nodes
                            //
                            v = a+(b-a)*amp::ampf<Precision>::getRandom();
                            vm = v+(b-a)*(ap::randominteger(5)-2);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, v)-spline1d::spline1dcalc<Precision>(c, vm)));
                            spline1d::spline1ddiff<Precision>(c, v, s, ds, d2s);
                            spline1d::spline1ddiff<Precision>(c, vm, s2, ds2, d2s2);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(s-s2));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(ds-ds2));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(d2s-d2s2));
                        }
                        else
                        {
                            
                            //
                            // * check values at nodes
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(y(i)-spline1d::spline1dcalc<Precision>(c, x(i))));
                            }
                        }
                        cserrors = cserrors || err>threshold;
                        
                        //
                        // check boundary conditions
                        //
                        err = 0;
                        if( bltype==0 )
                        {
                            spline1d::spline1ddiff<Precision>(c, a-h, s, ds, d2s);
                            spline1d::spline1ddiff<Precision>(c, a+h, s2, ds2, d2s2);
                            t = (d2s2-d2s)/(2*h);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(t));
                        }
                        if( bltype==1 )
                        {
                            t = (spline1d::spline1dcalc<Precision>(c, a+h)-spline1d::spline1dcalc<Precision>(c, a-h))/(2*h);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(bl-t));
                        }
                        if( bltype==2 )
                        {
                            t = (spline1d::spline1dcalc<Precision>(c, a+h)-2*spline1d::spline1dcalc<Precision>(c, a)+spline1d::spline1dcalc<Precision>(c, a-h))/amp::sqr<Precision>(h);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(bl-t));
                        }
                        if( brtype==0 )
                        {
                            spline1d::spline1ddiff<Precision>(c, b-h, s, ds, d2s);
                            spline1d::spline1ddiff<Precision>(c, b+h, s2, ds2, d2s2);
                            t = (d2s2-d2s)/(2*h);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(t));
                        }
                        if( brtype==1 )
                        {
                            t = (spline1d::spline1dcalc<Precision>(c, b+h)-spline1d::spline1dcalc<Precision>(c, b-h))/(2*h);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(br-t));
                        }
                        if( brtype==2 )
                        {
                            t = (spline1d::spline1dcalc<Precision>(c, b+h)-2*spline1d::spline1dcalc<Precision>(c, b)+spline1d::spline1dcalc<Precision>(c, b-h))/amp::sqr<Precision>(h);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(br-t));
                        }
                        if( bltype==-1 || brtype==-1 )
                        {
                            spline1d::spline1ddiff<Precision>(c, a+100*amp::ampf<Precision>::getAlgoPascalEpsilon(), s, ds, d2s);
                            spline1d::spline1ddiff<Precision>(c, b-100*amp::ampf<Precision>::getAlgoPascalEpsilon(), s2, ds2, d2s2);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(s-s2));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(ds-ds2));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(d2s-d2s2));
                        }
                        cserrors = cserrors || err>amp::ampf<Precision>("1.0E-3");
                        
                        //
                        // Check Lipschitz continuity
                        //
                        lconst<Precision>(a, b, c, lstep, l10, l11, l12);
                        lconst<Precision>(a, b, c, lstep/3, l20, l21, l22);
                        if( l10>amp::ampf<Precision>("1.0E-6") )
                        {
                            cserrors = cserrors || l20/l10>amp::ampf<Precision>("1.2");
                        }
                        if( l11>amp::ampf<Precision>("1.0E-6") )
                        {
                            cserrors = cserrors || l21/l11>amp::ampf<Precision>("1.2");
                        }
                        if( l12>amp::ampf<Precision>("1.0E-6") )
                        {
                            cserrors = cserrors || l22/l12>amp::ampf<Precision>("1.2");
                        }
                    }
                }
                
                //
                // Build Catmull-Rom spline.
                // Test for interpolation scheme properties:
                // * values at nodes
                // * boundary conditions
                // * continuous function
                // * continuous first derivative
                // * periodicity properties
                //
                for(bltype=-1; bltype<=0; bltype++)
                {
                    periodiccond = bltype==-1;
                    
                    //
                    // select random tension value, then build
                    //
                    if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
                    {
                        if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
                        {
                            tension = 0;
                        }
                        else
                        {
                            tension = 1;
                        }
                    }
                    else
                    {
                        tension = amp::ampf<Precision>::getRandom();
                    }
                    if( periodiccond )
                    {
                        spline1d::spline1dbuildcatmullrom<Precision>(x, yp, n, bltype, tension, c);
                    }
                    else
                    {
                        spline1d::spline1dbuildcatmullrom<Precision>(x, y, n, bltype, tension, c);
                    }
                    
                    //
                    // interpolation properties
                    //
                    err = 0;
                    if( periodiccond )
                    {
                        
                        //
                        // * check values at nodes; spline is periodic so
                        //   we add random number of periods to nodes
                        // * we also test for periodicity of first derivative
                        //
                        for(i=0; i<=n-1; i++)
                        {
                            v = x(i);
                            vm = v+(b-a)*(ap::randominteger(5)-2);
                            t = yp(i)-spline1d::spline1dcalc<Precision>(c, vm);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(t));
                            spline1d::spline1ddiff<Precision>(c, v, s, ds, d2s);
                            spline1d::spline1ddiff<Precision>(c, vm, s2, ds2, d2s2);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(s-s2));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(ds-ds2));
                        }
                        
                        //
                        // periodicity between nodes
                        //
                        v = a+(b-a)*amp::ampf<Precision>::getRandom();
                        vm = v+(b-a)*(ap::randominteger(5)-2);
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, v)-spline1d::spline1dcalc<Precision>(c, vm)));
                        spline1d::spline1ddiff<Precision>(c, v, s, ds, d2s);
                        spline1d::spline1ddiff<Precision>(c, vm, s2, ds2, d2s2);
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(s-s2));
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(ds-ds2));
                    }
                    else
                    {
                        
                        //
                        // * check values at nodes
                        //
                        for(i=0; i<=n-1; i++)
                        {
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(y(i)-spline1d::spline1dcalc<Precision>(c, x(i))));
                        }
                    }
                    crserrors = crserrors || err>threshold;
                    
                    //
                    // check boundary conditions
                    //
                    err = 0;
                    if( bltype==0 )
                    {
                        spline1d::spline1ddiff<Precision>(c, a-h, s, ds, d2s);
                        spline1d::spline1ddiff<Precision>(c, a+h, s2, ds2, d2s2);
                        t = (d2s2-d2s)/(2*h);
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(t));
                        spline1d::spline1ddiff<Precision>(c, b-h, s, ds, d2s);
                        spline1d::spline1ddiff<Precision>(c, b+h, s2, ds2, d2s2);
                        t = (d2s2-d2s)/(2*h);
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(t));
                    }
                    if( bltype==-1 )
                    {
                        spline1d::spline1ddiff<Precision>(c, a+100*amp::ampf<Precision>::getAlgoPascalEpsilon(), s, ds, d2s);
                        spline1d::spline1ddiff<Precision>(c, b-100*amp::ampf<Precision>::getAlgoPascalEpsilon(), s2, ds2, d2s2);
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(s-s2));
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(ds-ds2));
                    }
                    crserrors = crserrors || err>amp::ampf<Precision>("1.0E-3");
                    
                    //
                    // Check Lipschitz continuity
                    //
                    lconst<Precision>(a, b, c, lstep, l10, l11, l12);
                    lconst<Precision>(a, b, c, lstep/3, l20, l21, l22);
                    if( l10>amp::ampf<Precision>("1.0E-6") )
                    {
                        crserrors = crserrors || l20/l10>amp::ampf<Precision>("1.2");
                    }
                    if( l11>amp::ampf<Precision>("1.0E-6") )
                    {
                        crserrors = crserrors || l21/l11>amp::ampf<Precision>("1.2");
                    }
                }
                
                //
                // Build Hermite spline.
                // Test for interpolation scheme properties:
                // * values and derivatives at nodes
                // * continuous function
                // * continuous first derivative
                //
                spline1d::spline1dbuildhermite<Precision>(x, y, d, n, c);
                err = 0;
                for(i=0; i<=n-1; i++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(y(i)-spline1d::spline1dcalc<Precision>(c, x(i))));
                }
                hserrors = hserrors || err>threshold;
                err = 0;
                for(i=0; i<=n-1; i++)
                {
                    t = (spline1d::spline1dcalc<Precision>(c, x(i)+h)-spline1d::spline1dcalc<Precision>(c, x(i)-h))/(2*h);
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(d(i)-t));
                }
                hserrors = hserrors || err>amp::ampf<Precision>("1.0E-3");
                lconst<Precision>(a, b, c, lstep, l10, l11, l12);
                lconst<Precision>(a, b, c, lstep/3, l20, l21, l22);
                hserrors = hserrors || l20/l10>amp::ampf<Precision>("1.2");
                hserrors = hserrors || l21/l11>amp::ampf<Precision>("1.2");
                
                //
                // Build Akima spline
                // Test for general interpolation scheme properties:
                // * values at nodes
                // * continuous function
                // * continuous first derivative
                // Test for specific properties is implemented below.
                //
                if( n>=5 )
                {
                    spline1d::spline1dbuildakima<Precision>(x, y, n, c);
                    err = 0;
                    for(i=0; i<=n-1; i++)
                    {
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(y(i)-spline1d::spline1dcalc<Precision>(c, x(i))));
                    }
                    aserrors = aserrors || err>threshold;
                    lconst<Precision>(a, b, c, lstep, l10, l11, l12);
                    lconst<Precision>(a, b, c, lstep/3, l20, l21, l22);
                    hserrors = hserrors || l20/l10>amp::ampf<Precision>("1.2");
                    hserrors = hserrors || l21/l11>amp::ampf<Precision>("1.2");
                }
            }
        }
        
        //
        // Special linear spline test:
        // test for linearity between x[i] and x[i+1]
        //
        for(n=2; n<=maxn; n++)
        {
            x.setbounds(0, n-1);
            y.setbounds(0, n-1);
            
            //
            // Prepare task
            //
            a = -1;
            b = +1;
            for(i=0; i<=n-1; i++)
            {
                x(i) = a+(b-a)*i/(n-1);
                y(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            spline1d::spline1dbuildlinear<Precision>(x, y, n, c);
            
            //
            // Test
            //
            err = 0;
            for(k=0; k<=n-2; k++)
            {
                a = x(k);
                b = x(k+1);
                for(pass=1; pass<=passcount; pass++)
                {
                    t = a+(b-a)*amp::ampf<Precision>::getRandom();
                    v = y(k)+(t-a)/(b-a)*(y(k+1)-y(k));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, t)-v));
                }
            }
            lserrors = lserrors || err>threshold;
        }
        
        //
        // Special Akima test: test outlier sensitivity
        // Spline value at (x[i], x[i+1]) should depend from
        // f[i-2], f[i-1], f[i], f[i+1], f[i+2], f[i+3] only.
        //
        for(n=5; n<=maxn; n++)
        {
            x.setbounds(0, n-1);
            y.setbounds(0, n-1);
            y2.setbounds(0, n-1);
            
            //
            // Prepare unperturbed Akima spline
            //
            a = -1;
            b = +1;
            for(i=0; i<=n-1; i++)
            {
                x(i) = a+(b-a)*i/(n-1);
                y(i) = amp::cos<Precision>(amp::ampf<Precision>("1.3")*amp::pi<Precision>()*x(i)+amp::ampf<Precision>("0.4"));
            }
            spline1d::spline1dbuildakima<Precision>(x, y, n, c);
            
            //
            // Process perturbed tasks
            //
            err = 0;
            for(k=0; k<=n-1; k++)
            {
                amp::vmove(y2.getvector(0, n-1), y.getvector(0, n-1));
                y2(k) = 5;
                spline1d::spline1dbuildakima<Precision>(x, y2, n, c2);
                
                //
                // Test left part independence
                //
                if( k-3>=1 )
                {
                    a = -1;
                    b = x(k-3);
                    for(pass=1; pass<=passcount; pass++)
                    {
                        t = a+(b-a)*amp::ampf<Precision>::getRandom();
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, t)-spline1d::spline1dcalc<Precision>(c2, t)));
                    }
                }
                
                //
                // Test right part independence
                //
                if( k+3<=n-2 )
                {
                    a = x(k+3);
                    b = +1;
                    for(pass=1; pass<=passcount; pass++)
                    {
                        t = a+(b-a)*amp::ampf<Precision>::getRandom();
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, t)-spline1d::spline1dcalc<Precision>(c2, t)));
                    }
                }
            }
            aserrors = aserrors || err>threshold;
        }
        
        //
        // Differentiation, copy/unpack test
        //
        for(n=2; n<=maxn; n++)
        {
            x.setbounds(0, n-1);
            y.setbounds(0, n-1);
            
            //
            // Prepare cubic spline
            //
            a = -1-amp::ampf<Precision>::getRandom();
            b = +1+amp::ampf<Precision>::getRandom();
            for(i=0; i<=n-1; i++)
            {
                x(i) = a+(b-a)*i/(n-1);
                y(i) = amp::cos<Precision>(amp::ampf<Precision>("1.3")*amp::pi<Precision>()*x(i)+amp::ampf<Precision>("0.4"));
            }
            spline1d::spline1dbuildcubic<Precision>(x, y, n, 2, amp::ampf<Precision>("0.0"), 2, amp::ampf<Precision>("0.0"), c);
            
            //
            // Test diff
            //
            err = 0;
            for(pass=1; pass<=passcount; pass++)
            {
                t = a+(b-a)*amp::ampf<Precision>::getRandom();
                spline1d::spline1ddiff<Precision>(c, t, s, ds, d2s);
                vl = spline1d::spline1dcalc<Precision>(c, t-h);
                vm = spline1d::spline1dcalc<Precision>(c, t);
                vr = spline1d::spline1dcalc<Precision>(c, t+h);
                err = amp::maximum<Precision>(err, amp::abs<Precision>(s-vm));
                err = amp::maximum<Precision>(err, amp::abs<Precision>(ds-(vr-vl)/(2*h)));
                err = amp::maximum<Precision>(err, amp::abs<Precision>(d2s-(vr-2*vm+vl)/amp::sqr<Precision>(h)));
            }
            dserrors = dserrors || err>amp::ampf<Precision>("0.001");
            
            //
            // Test copy
            //
            unsetspline1d<Precision>(c2);
            spline1d::spline1dcopy<Precision>(c, c2);
            err = 0;
            for(pass=1; pass<=passcount; pass++)
            {
                t = a+(b-a)*amp::ampf<Precision>::getRandom();
                err = amp::maximum<Precision>(err, amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, t)-spline1d::spline1dcalc<Precision>(c2, t)));
            }
            cperrors = cperrors || err>threshold;
            
            //
            // Test unpack
            //
            uperrors = uperrors || !testunpack<Precision>(c, x);
            
            //
            // Test lin.trans.
            //
            err = 0;
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // LinTransX, general A
                //
                sa = 4*amp::ampf<Precision>::getRandom()-2;
                sb = 2*amp::ampf<Precision>::getRandom()-1;
                t = a+(b-a)*amp::ampf<Precision>::getRandom();
                spline1d::spline1dcopy<Precision>(c, c2);
                spline1d::spline1dlintransx<Precision>(c2, sa, sb);
                err = amp::maximum<Precision>(err, amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, t)-spline1d::spline1dcalc<Precision>(c2, (t-sb)/sa)));
                
                //
                // LinTransX, special case: A=0
                //
                sb = 2*amp::ampf<Precision>::getRandom()-1;
                t = a+(b-a)*amp::ampf<Precision>::getRandom();
                spline1d::spline1dcopy<Precision>(c, c2);
                spline1d::spline1dlintransx<Precision>(c2, amp::ampf<Precision>(0), sb);
                err = amp::maximum<Precision>(err, amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, sb)-spline1d::spline1dcalc<Precision>(c2, t)));
                
                //
                // LinTransY
                //
                sa = 2*amp::ampf<Precision>::getRandom()-1;
                sb = 2*amp::ampf<Precision>::getRandom()-1;
                t = a+(b-a)*amp::ampf<Precision>::getRandom();
                spline1d::spline1dcopy<Precision>(c, c2);
                spline1d::spline1dlintransy<Precision>(c2, sa, sb);
                err = amp::maximum<Precision>(err, amp::abs<Precision>(sa*spline1d::spline1dcalc<Precision>(c, t)+sb-spline1d::spline1dcalc<Precision>(c2, t)));
            }
            lterrors = lterrors || err>threshold;
        }
        
        //
        // Testing integration.
        // Three tests are performed:
        //
        // * approximate test (well behaved smooth function, many points,
        //   integration inside [a,b]), non-periodic spline
        //
        // * exact test (integration of parabola, outside of [a,b], non-periodic spline
        //
        // * approximate test for periodic splines. F(x)=cos(2*pi*x)+1.
        //   Period length is equals to 1.0, so all operations with
        //   multiples of period are done exactly. For each value of PERIOD
        //   we calculate and test integral at four points:
        //   -   0 < t0 < PERIOD
        //   -   t1 = PERIOD-eps
        //   -   t2 = PERIOD
        //   -   t3 = PERIOD+eps
        //
        err = 0;
        for(n=20; n<=35; n++)
        {
            x.setbounds(0, n-1);
            y.setbounds(0, n-1);
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // Prepare cubic spline
                //
                a = -1-amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom();
                b = +1+amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom();
                for(i=0; i<=n-1; i++)
                {
                    x(i) = a+(b-a)*i/(n-1);
                    y(i) = amp::sin<Precision>(amp::pi<Precision>()*x(i)+amp::ampf<Precision>("0.4"))+amp::exp<Precision>(x(i));
                }
                bl = amp::pi<Precision>()*amp::cos<Precision>(amp::pi<Precision>()*a+amp::ampf<Precision>("0.4"))+amp::exp<Precision>(a);
                br = amp::pi<Precision>()*amp::cos<Precision>(amp::pi<Precision>()*b+amp::ampf<Precision>("0.4"))+amp::exp<Precision>(b);
                spline1d::spline1dbuildcubic<Precision>(x, y, n, 1, bl, 1, br, c);
                
                //
                // Test
                //
                t = a+(b-a)*amp::ampf<Precision>::getRandom();
                v = -amp::cos<Precision>(amp::pi<Precision>()*a+amp::ampf<Precision>("0.4"))/amp::pi<Precision>()+amp::exp<Precision>(a);
                v = -amp::cos<Precision>(amp::pi<Precision>()*t+amp::ampf<Precision>("0.4"))/amp::pi<Precision>()+amp::exp<Precision>(t)-v;
                v = v-spline1d::spline1dintegrate<Precision>(c, t);
                err = amp::maximum<Precision>(err, amp::abs<Precision>(v));
            }
        }
        ierrors = ierrors || err>amp::ampf<Precision>("0.001");
        p0 = 2*amp::ampf<Precision>::getRandom()-1;
        p1 = 2*amp::ampf<Precision>::getRandom()-1;
        p2 = 2*amp::ampf<Precision>::getRandom()-1;
        a = -amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
        b = +amp::ampf<Precision>::getRandom()+amp::ampf<Precision>("0.5");
        n = 2;
        x.setlength(n);
        y.setlength(n);
        d.setlength(n);
        x(0) = a;
        y(0) = p0+p1*a+p2*amp::sqr<Precision>(a);
        d(0) = p1+2*p2*a;
        x(1) = b;
        y(1) = p0+p1*b+p2*amp::sqr<Precision>(b);
        d(1) = p1+2*p2*b;
        spline1d::spline1dbuildhermite<Precision>(x, y, d, n, c);
        bl = amp::minimum<Precision>(a, b)-amp::abs<Precision>(b-a);
        br = amp::minimum<Precision>(a, b)+amp::abs<Precision>(b-a);
        err = 0;
        for(pass=1; pass<=100; pass++)
        {
            t = bl+(br-bl)*amp::ampf<Precision>::getRandom();
            v = p0*t+p1*amp::sqr<Precision>(t)/2+p2*amp::sqr<Precision>(t)*t/3-(p0*a+p1*amp::sqr<Precision>(a)/2+p2*amp::sqr<Precision>(a)*a/3);
            v = v-spline1d::spline1dintegrate<Precision>(c, t);
            err = amp::maximum<Precision>(err, amp::abs<Precision>(v));
        }
        ierrors = ierrors || err>threshold;
        n = 100;
        x.setlength(n);
        y.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            x(i) = amp::ampf<Precision>(i)/(amp::ampf<Precision>(n-1));
            y(i) = amp::cos<Precision>(2*amp::pi<Precision>()*x(i))+1;
        }
        y(0) = 2;
        y(n-1) = 2;
        spline1d::spline1dbuildcubic<Precision>(x, y, n, -1, amp::ampf<Precision>("0.0"), -1, amp::ampf<Precision>("0.0"), c);
        intab = spline1d::spline1dintegrate<Precision>(c, amp::ampf<Precision>("1.0"));
        v = amp::ampf<Precision>::getRandom();
        vr = spline1d::spline1dintegrate<Precision>(c, v);
        ierrors = ierrors || amp::abs<Precision>(intab-1)>amp::ampf<Precision>("0.001");
        for(i=-10; i<=10; i++)
        {
            ierrors = ierrors || amp::abs<Precision>(spline1d::spline1dintegrate<Precision>(c, i+v)-(i*intab+vr))>amp::ampf<Precision>("0.001");
            ierrors = ierrors || amp::abs<Precision>(spline1d::spline1dintegrate<Precision>(c, i-1000*amp::ampf<Precision>::getAlgoPascalEpsilon())-i*intab)>amp::ampf<Precision>("0.001");
            ierrors = ierrors || amp::abs<Precision>(spline1d::spline1dintegrate<Precision>(c, amp::ampf<Precision>(i))-i*intab)>amp::ampf<Precision>("0.001");
            ierrors = ierrors || amp::abs<Precision>(spline1d::spline1dintegrate<Precision>(c, i+1000*amp::ampf<Precision>::getAlgoPascalEpsilon())-i*intab)>amp::ampf<Precision>("0.001");
        }
        
        //
        // Test fitting.
        //
        for(pass=1; pass<=passcount; pass++)
        {
            
            //
            // Cubic splines
            // Ability to handle boundary constraints (1-4 constraints on F, dF/dx).
            //
            for(m=4; m<=8; m++)
            {
                for(k=1; k<=4; k++)
                {
                    if( k>=m )
                    {
                        continue;
                    }
                    n = 100;
                    x.setlength(n);
                    y.setlength(n);
                    w.setlength(n);
                    xc.setlength(4);
                    yc.setlength(4);
                    dc.setlength(4);
                    sa = 1+amp::ampf<Precision>::getRandom();
                    sb = 2*amp::ampf<Precision>::getRandom()-1;
                    for(i=0; i<=n-1; i++)
                    {
                        x(i) = sa*amp::ampf<Precision>::getRandom()+sb;
                        y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        w(i) = 1+amp::ampf<Precision>::getRandom();
                    }
                    xc(0) = sb;
                    yc(0) = 2*amp::ampf<Precision>::getRandom()-1;
                    dc(0) = 0;
                    xc(1) = sb;
                    yc(1) = 2*amp::ampf<Precision>::getRandom()-1;
                    dc(1) = 1;
                    xc(2) = sa+sb;
                    yc(2) = 2*amp::ampf<Precision>::getRandom()-1;
                    dc(2) = 0;
                    xc(3) = sa+sb;
                    yc(3) = 2*amp::ampf<Precision>::getRandom()-1;
                    dc(3) = 1;
                    spline1d::spline1dfitcubicwc<Precision>(x, y, w, n, xc, yc, dc, k, m, info, c, rep);
                    if( info<=0 )
                    {
                        fiterrors = true;
                    }
                    else
                    {
                        
                        //
                        // Check that constraints are satisfied
                        //
                        for(i=0; i<=k-1; i++)
                        {
                            spline1d::spline1ddiff<Precision>(c, xc(i), s, ds, d2s);
                            if( dc(i)==0 )
                            {
                                fiterrors = fiterrors || amp::abs<Precision>(s-yc(i))>threshold;
                            }
                            if( dc(i)==1 )
                            {
                                fiterrors = fiterrors || amp::abs<Precision>(ds-yc(i))>threshold;
                            }
                            if( dc(i)==2 )
                            {
                                fiterrors = fiterrors || amp::abs<Precision>(d2s-yc(i))>threshold;
                            }
                        }
                    }
                }
            }
            
            //
            // Cubic splines
            // Ability to handle one internal constraint
            //
            for(m=4; m<=8; m++)
            {
                n = 100;
                x.setlength(n);
                y.setlength(n);
                w.setlength(n);
                xc.setlength(1);
                yc.setlength(1);
                dc.setlength(1);
                sa = 1+amp::ampf<Precision>::getRandom();
                sb = 2*amp::ampf<Precision>::getRandom()-1;
                for(i=0; i<=n-1; i++)
                {
                    x(i) = sa*amp::ampf<Precision>::getRandom()+sb;
                    y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    w(i) = 1+amp::ampf<Precision>::getRandom();
                }
                xc(0) = sa*amp::ampf<Precision>::getRandom()+sb;
                yc(0) = 2*amp::ampf<Precision>::getRandom()-1;
                dc(0) = ap::randominteger(2);
                spline1d::spline1dfitcubicwc<Precision>(x, y, w, n, xc, yc, dc, 1, m, info, c, rep);
                if( info<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    
                    //
                    // Check that constraints are satisfied
                    //
                    spline1d::spline1ddiff<Precision>(c, xc(0), s, ds, d2s);
                    if( dc(0)==0 )
                    {
                        fiterrors = fiterrors || amp::abs<Precision>(s-yc(0))>threshold;
                    }
                    if( dc(0)==1 )
                    {
                        fiterrors = fiterrors || amp::abs<Precision>(ds-yc(0))>threshold;
                    }
                    if( dc(0)==2 )
                    {
                        fiterrors = fiterrors || amp::abs<Precision>(d2s-yc(0))>threshold;
                    }
                }
            }
            
            //
            // Hermite splines
            // Ability to handle boundary constraints (1-4 constraints on F, dF/dx).
            //
            for(m=4; m<=8; m++)
            {
                for(k=1; k<=4; k++)
                {
                    if( k>=m )
                    {
                        continue;
                    }
                    if( m%2!=0 )
                    {
                        continue;
                    }
                    n = 100;
                    x.setlength(n);
                    y.setlength(n);
                    w.setlength(n);
                    xc.setlength(4);
                    yc.setlength(4);
                    dc.setlength(4);
                    sa = 1+amp::ampf<Precision>::getRandom();
                    sb = 2*amp::ampf<Precision>::getRandom()-1;
                    for(i=0; i<=n-1; i++)
                    {
                        x(i) = sa*amp::ampf<Precision>::getRandom()+sb;
                        y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        w(i) = 1+amp::ampf<Precision>::getRandom();
                    }
                    xc(0) = sb;
                    yc(0) = 2*amp::ampf<Precision>::getRandom()-1;
                    dc(0) = 0;
                    xc(1) = sb;
                    yc(1) = 2*amp::ampf<Precision>::getRandom()-1;
                    dc(1) = 1;
                    xc(2) = sa+sb;
                    yc(2) = 2*amp::ampf<Precision>::getRandom()-1;
                    dc(2) = 0;
                    xc(3) = sa+sb;
                    yc(3) = 2*amp::ampf<Precision>::getRandom()-1;
                    dc(3) = 1;
                    spline1d::spline1dfithermitewc<Precision>(x, y, w, n, xc, yc, dc, k, m, info, c, rep);
                    if( info<=0 )
                    {
                        fiterrors = true;
                    }
                    else
                    {
                        
                        //
                        // Check that constraints are satisfied
                        //
                        for(i=0; i<=k-1; i++)
                        {
                            spline1d::spline1ddiff<Precision>(c, xc(i), s, ds, d2s);
                            if( dc(i)==0 )
                            {
                                fiterrors = fiterrors || amp::abs<Precision>(s-yc(i))>threshold;
                            }
                            if( dc(i)==1 )
                            {
                                fiterrors = fiterrors || amp::abs<Precision>(ds-yc(i))>threshold;
                            }
                            if( dc(i)==2 )
                            {
                                fiterrors = fiterrors || amp::abs<Precision>(d2s-yc(i))>threshold;
                            }
                        }
                    }
                }
            }
            
            //
            // Hermite splines
            // Ability to handle one internal constraint
            //
            for(m=4; m<=8; m++)
            {
                if( m%2!=0 )
                {
                    continue;
                }
                n = 100;
                x.setlength(n);
                y.setlength(n);
                w.setlength(n);
                xc.setlength(1);
                yc.setlength(1);
                dc.setlength(1);
                sa = 1+amp::ampf<Precision>::getRandom();
                sb = 2*amp::ampf<Precision>::getRandom()-1;
                for(i=0; i<=n-1; i++)
                {
                    x(i) = sa*amp::ampf<Precision>::getRandom()+sb;
                    y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    w(i) = 1+amp::ampf<Precision>::getRandom();
                }
                xc(0) = sa*amp::ampf<Precision>::getRandom()+sb;
                yc(0) = 2*amp::ampf<Precision>::getRandom()-1;
                dc(0) = ap::randominteger(2);
                spline1d::spline1dfithermitewc<Precision>(x, y, w, n, xc, yc, dc, 1, m, info, c, rep);
                if( info<=0 )
                {
                    fiterrors = true;
                }
                else
                {
                    
                    //
                    // Check that constraints are satisfied
                    //
                    spline1d::spline1ddiff<Precision>(c, xc(0), s, ds, d2s);
                    if( dc(0)==0 )
                    {
                        fiterrors = fiterrors || amp::abs<Precision>(s-yc(0))>threshold;
                    }
                    if( dc(0)==1 )
                    {
                        fiterrors = fiterrors || amp::abs<Precision>(ds-yc(0))>threshold;
                    }
                    if( dc(0)==2 )
                    {
                        fiterrors = fiterrors || amp::abs<Precision>(d2s-yc(0))>threshold;
                    }
                }
            }
        }
        for(m=4; m<=8; m++)
        {
            for(stype=0; stype<=1; stype++)
            {
                for(pass=1; pass<=passcount; pass++)
                {
                    if( stype==1 && m%2!=0 )
                    {
                        continue;
                    }
                    
                    //
                    // cubic/Hermite spline fitting:
                    // * generate "template spline" C2
                    // * generate 2*N points from C2, such that result of
                    //   ideal fit should be equal to C2
                    // * fit, store in C
                    // * compare C and C2
                    //
                    sa = 1+amp::ampf<Precision>::getRandom();
                    sb = 2*amp::ampf<Precision>::getRandom()-1;
                    if( stype==0 )
                    {
                        x.setlength(m-2);
                        y.setlength(m-2);
                        for(i=0; i<=m-2-1; i++)
                        {
                            x(i) = sa*i/(m-2-1)+sb;
                            y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                        spline1d::spline1dbuildcubic<Precision>(x, y, m-2, 1, 2*amp::ampf<Precision>::getRandom()-1, 1, 2*amp::ampf<Precision>::getRandom()-1, c2);
                    }
                    if( stype==1 )
                    {
                        x.setlength(m/2);
                        y.setlength(m/2);
                        d.setlength(m/2);
                        for(i=0; i<=m/2-1; i++)
                        {
                            x(i) = sa*i/(m/2-1)+sb;
                            y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                            d(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                        spline1d::spline1dbuildhermite<Precision>(x, y, d, m/2, c2);
                    }
                    n = 50;
                    x.setlength(2*n);
                    y.setlength(2*n);
                    w.setlength(2*n);
                    for(i=0; i<=n-1; i++)
                    {
                        
                        //
                        // "if i=0" and "if i=1" are needed to
                        // synchronize interval size for C2 and
                        // spline being fitted (i.e. C).
                        //
                        t = amp::ampf<Precision>::getRandom();
                        x(i) = sa*amp::ampf<Precision>::getRandom()+sb;
                        if( i==0 )
                        {
                            x(i) = sb;
                        }
                        if( i==1 )
                        {
                            x(i) = sa+sb;
                        }
                        v = spline1d::spline1dcalc<Precision>(c2, x(i));
                        y(i) = v+t;
                        w(i) = 1+amp::ampf<Precision>::getRandom();
                        x(n+i) = x(i);
                        y(n+i) = v-t;
                        w(n+i) = w(i);
                    }
                    if( stype==0 )
                    {
                        spline1d::spline1dfitcubicwc<Precision>(x, y, w, 2*n, xc, yc, dc, 0, m, info, c, rep);
                    }
                    if( stype==1 )
                    {
                        spline1d::spline1dfithermitewc<Precision>(x, y, w, 2*n, xc, yc, dc, 0, m, info, c, rep);
                    }
                    if( info<=0 )
                    {
                        fiterrors = true;
                    }
                    else
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            v = sa*amp::ampf<Precision>::getRandom()+sb;
                            fiterrors = fiterrors || amp::abs<Precision>(spline1d::spline1dcalc<Precision>(c, v)-spline1d::spline1dcalc<Precision>(c2, v))>threshold;
                        }
                    }
                }
            }
        }
        for(m=4; m<=8; m++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // prepare points/weights
                //
                sa = 1+amp::ampf<Precision>::getRandom();
                sb = 2*amp::ampf<Precision>::getRandom()-1;
                n = 10+ap::randominteger(10);
                x.setlength(n);
                y.setlength(n);
                w.setlength(n);
                for(i=0; i<=n-1; i++)
                {
                    x(i) = sa*amp::ampf<Precision>::getRandom()+sb;
                    y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    w(i) = 1;
                }
                
                //
                // Fit cubic with unity weights, without weights, then compare
                //
                if( m>=4 )
                {
                    spline1d::spline1dfitcubicwc<Precision>(x, y, w, n, xc, yc, dc, 0, m, info1, c, rep);
                    spline1d::spline1dfitcubic<Precision>(x, y, n, m, info2, c2, rep2);
                    if( info1<=0 || info2<=0 )
                    {
                        fiterrors = true;
                    }
                    else
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            v = sa*amp::ampf<Precision>::getRandom()+sb;
                            fiterrors = fiterrors || spline1d::spline1dcalc<Precision>(c, v)!=spline1d::spline1dcalc<Precision>(c2, v);
                            fiterrors = fiterrors || rep.taskrcond!=rep2.taskrcond;
                            fiterrors = fiterrors || rep.rmserror!=rep2.rmserror;
                            fiterrors = fiterrors || rep.avgerror!=rep2.avgerror;
                            fiterrors = fiterrors || rep.avgrelerror!=rep2.avgrelerror;
                            fiterrors = fiterrors || rep.maxerror!=rep2.maxerror;
                        }
                    }
                }
                
                //
                // Fit Hermite with unity weights, without weights, then compare
                //
                if( m>=4 && m%2==0 )
                {
                    spline1d::spline1dfithermitewc<Precision>(x, y, w, n, xc, yc, dc, 0, m, info1, c, rep);
                    spline1d::spline1dfithermite<Precision>(x, y, n, m, info2, c2, rep2);
                    if( info1<=0 || info2<=0 )
                    {
                        fiterrors = true;
                    }
                    else
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            v = sa*amp::ampf<Precision>::getRandom()+sb;
                            fiterrors = fiterrors || spline1d::spline1dcalc<Precision>(c, v)!=spline1d::spline1dcalc<Precision>(c2, v);
                            fiterrors = fiterrors || rep.taskrcond!=rep2.taskrcond;
                            fiterrors = fiterrors || rep.rmserror!=rep2.rmserror;
                            fiterrors = fiterrors || rep.avgerror!=rep2.avgerror;
                            fiterrors = fiterrors || rep.avgrelerror!=rep2.avgrelerror;
                            fiterrors = fiterrors || rep.maxerror!=rep2.maxerror;
                        }
                    }
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
            // Test cubic fitting
            //
            spline1d::spline1dfitcubic<Precision>(x, y, 4, 4, info, c, rep);
            if( info<=0 )
            {
                fiterrors = true;
            }
            else
            {
                s = spline1d::spline1dcalc<Precision>(c, amp::ampf<Precision>(0));
                fiterrors = fiterrors || amp::abs<Precision>(s-v)>threshold;
                fiterrors = fiterrors || amp::abs<Precision>(rep.rmserror-refrms)>threshold;
                fiterrors = fiterrors || amp::abs<Precision>(rep.avgerror-refavg)>threshold;
                fiterrors = fiterrors || amp::abs<Precision>(rep.avgrelerror-refavgrel)>threshold;
                fiterrors = fiterrors || amp::abs<Precision>(rep.maxerror-refmax)>threshold;
            }
            
            //
            // Test cubic fitting
            //
            spline1d::spline1dfithermite<Precision>(x, y, 4, 4, info, c, rep);
            if( info<=0 )
            {
                fiterrors = true;
            }
            else
            {
                s = spline1d::spline1dcalc<Precision>(c, amp::ampf<Precision>(0));
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
        waserrors = lserrors || cserrors || crserrors || hserrors || aserrors || dserrors || cperrors || uperrors || lterrors || ierrors || fiterrors;
        if( !silent )
        {
            printf("TESTING SPLINE INTERPOLATION\n");
            
            //
            // Normal tests
            //
            printf("LINEAR SPLINE TEST:                      ");
            if( lserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("CUBIC SPLINE TEST:                       ");
            if( cserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("CATMULL-ROM SPLINE TEST:                 ");
            if( crserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("HERMITE SPLINE TEST:                     ");
            if( hserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("AKIMA SPLINE TEST:                       ");
            if( aserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("DIFFERENTIATION TEST:                    ");
            if( dserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("COPY/SERIALIZATION TEST:                 ");
            if( cperrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("UNPACK TEST:                             ");
            if( uperrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("LIN.TRANS. TEST:                         ");
            if( lterrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("INTEGRATION TEST:                        ");
            if( ierrors )
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


    /*************************************************************************
    Lipschitz constants for spline inself, first and second derivatives.
    *************************************************************************/
    template<unsigned int Precision>
    void lconst(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const spline1d::spline1dinterpolant<Precision>& c,
        amp::ampf<Precision> lstep,
        amp::ampf<Precision>& l0,
        amp::ampf<Precision>& l1,
        amp::ampf<Precision>& l2)
    {
        amp::ampf<Precision> t;
        amp::ampf<Precision> vl;
        amp::ampf<Precision> vm;
        amp::ampf<Precision> vr;
        amp::ampf<Precision> prevf;
        amp::ampf<Precision> prevd;
        amp::ampf<Precision> prevd2;
        amp::ampf<Precision> f;
        amp::ampf<Precision> d;
        amp::ampf<Precision> d2;


        l0 = 0;
        l1 = 0;
        l2 = 0;
        t = a-amp::ampf<Precision>("0.1");
        vl = spline1d::spline1dcalc<Precision>(c, t-2*lstep);
        vm = spline1d::spline1dcalc<Precision>(c, t-lstep);
        vr = spline1d::spline1dcalc<Precision>(c, t);
        f = vm;
        d = (vr-vl)/(2*lstep);
        d2 = (vr-2*vm+vl)/amp::sqr<Precision>(lstep);
        while( t<=b+amp::ampf<Precision>("0.1") )
        {
            prevf = f;
            prevd = d;
            prevd2 = d2;
            vl = vm;
            vm = vr;
            vr = spline1d::spline1dcalc<Precision>(c, t+lstep);
            f = vm;
            d = (vr-vl)/(2*lstep);
            d2 = (vr-2*vm+vl)/amp::sqr<Precision>(lstep);
            l0 = amp::maximum<Precision>(l0, amp::abs<Precision>((f-prevf)/lstep));
            l1 = amp::maximum<Precision>(l1, amp::abs<Precision>((d-prevd)/lstep));
            l2 = amp::maximum<Precision>(l2, amp::abs<Precision>((d2-prevd2)/lstep));
            t = t+lstep;
        }
    }


    /*************************************************************************
    Unpack testing
    *************************************************************************/
    template<unsigned int Precision>
    bool testunpack(const spline1d::spline1dinterpolant<Precision>& c,
        const ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        bool result;
        int i;
        int n;
        amp::ampf<Precision> err;
        amp::ampf<Precision> t;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        int pass;
        int passcount;
        ap::template_2d_array< amp::ampf<Precision> > tbl;


        passcount = 20;
        err = 0;
        spline1d::spline1dunpack<Precision>(c, n, tbl);
        for(i=0; i<=n-2; i++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                t = amp::ampf<Precision>::getRandom()*(tbl(i,1)-tbl(i,0));
                v1 = tbl(i,2)+t*tbl(i,3)+amp::sqr<Precision>(t)*tbl(i,4)+t*amp::sqr<Precision>(t)*tbl(i,5);
                v2 = spline1d::spline1dcalc<Precision>(c, tbl(i,0)+t);
                err = amp::maximum<Precision>(err, amp::abs<Precision>(v1-v2));
            }
        }
        for(i=0; i<=n-2; i++)
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(i)-tbl(i,0)));
        }
        for(i=0; i<=n-2; i++)
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(i+1)-tbl(i,1)));
        }
        result = err<100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        return result;
    }


    /*************************************************************************
    Unset spline, i.e. initialize it with random garbage
    *************************************************************************/
    template<unsigned int Precision>
    void unsetspline1d(spline1d::spline1dinterpolant<Precision>& c)
    {
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > d;


        x.setlength(2);
        y.setlength(2);
        d.setlength(2);
        x(0) = -1;
        y(0) = amp::ampf<Precision>::getRandom();
        d(0) = amp::ampf<Precision>::getRandom();
        x(1) = 1;
        y(1) = amp::ampf<Precision>::getRandom();
        d(1) = amp::ampf<Precision>::getRandom();
        spline1d::spline1dbuildhermite<Precision>(x, y, d, 2, c);
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
    Tests whether constant C is solution of 1D LLS problem
    *************************************************************************/
    template<unsigned int Precision>
    bool is1dsolution(int n,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        amp::ampf<Precision> c)
    {
        bool result;
        int i;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s1;
        amp::ampf<Precision> s2;
        amp::ampf<Precision> s3;
        amp::ampf<Precision> delta;


        delta = amp::ampf<Precision>("0.001");
        
        //
        // Test result
        //
        s1 = 0;
        for(i=0; i<=n-1; i++)
        {
            s1 = s1+amp::sqr<Precision>(w(i)*(c-y(i)));
        }
        s2 = 0;
        s3 = 0;
        for(i=0; i<=n-1; i++)
        {
            s2 = s2+amp::sqr<Precision>(w(i)*(c+delta-y(i)));
            s3 = s3+amp::sqr<Precision>(w(i)*(c-delta-y(i)));
        }
        result = s2>=s1 && s3>=s1;
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testspline1dunit_test_silent()
    {
        bool result;


        result = testsplineinterpolation<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testspline1dunit_test()
    {
        bool result;


        result = testsplineinterpolation<Precision>(false);
        return result;
    }
} // namespace

#endif
