
#ifndef _testspline2dunit_h
#define _testspline2dunit_h

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
#include "spline2d.h"
namespace testspline2dunit
{
    template<unsigned int Precision>
    bool test2dinterpolation(bool silent);
    template<unsigned int Precision>
    void lconst(const spline2d::spline2dinterpolant<Precision>& c,
        const ap::template_1d_array< amp::ampf<Precision> >& lx,
        const ap::template_1d_array< amp::ampf<Precision> >& ly,
        int m,
        int n,
        amp::ampf<Precision> lstep,
        amp::ampf<Precision>& lc,
        amp::ampf<Precision>& lcx,
        amp::ampf<Precision>& lcy,
        amp::ampf<Precision>& lcxy);
    template<unsigned int Precision>
    void twodnumder(const spline2d::spline2dinterpolant<Precision>& c,
        amp::ampf<Precision> x,
        amp::ampf<Precision> y,
        amp::ampf<Precision> h,
        amp::ampf<Precision>& f,
        amp::ampf<Precision>& fx,
        amp::ampf<Precision>& fy,
        amp::ampf<Precision>& fxy);
    template<unsigned int Precision>
    bool testunpack(const spline2d::spline2dinterpolant<Precision>& c,
        const ap::template_1d_array< amp::ampf<Precision> >& lx,
        const ap::template_1d_array< amp::ampf<Precision> >& ly);
    template<unsigned int Precision>
    bool testlintrans(const spline2d::spline2dinterpolant<Precision>& c,
        amp::ampf<Precision> ax,
        amp::ampf<Precision> bx,
        amp::ampf<Precision> ay,
        amp::ampf<Precision> by);
    template<unsigned int Precision>
    void unsetspline2d(spline2d::spline2dinterpolant<Precision>& c);
    template<unsigned int Precision>
    bool testspline2dunit_test_silent();
    template<unsigned int Precision>
    bool testspline2dunit_test();


    template<unsigned int Precision>
    bool test2dinterpolation(bool silent)
    {
        bool result;
        bool waserrors;
        bool blerrors;
        bool bcerrors;
        bool dserrors;
        bool cperrors;
        bool uperrors;
        bool lterrors;
        bool syerrors;
        bool rlerrors;
        bool rcerrors;
        int pass;
        int passcount;
        int jobtype;
        amp::ampf<Precision> lstep;
        amp::ampf<Precision> h;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        spline2d::spline2dinterpolant<Precision> c;
        spline2d::spline2dinterpolant<Precision> c2;
        ap::template_1d_array< amp::ampf<Precision> > lx;
        ap::template_1d_array< amp::ampf<Precision> > ly;
        ap::template_2d_array< amp::ampf<Precision> > f;
        ap::template_2d_array< amp::ampf<Precision> > fr;
        ap::template_2d_array< amp::ampf<Precision> > ft;
        amp::ampf<Precision> ax;
        amp::ampf<Precision> ay;
        amp::ampf<Precision> bx;
        amp::ampf<Precision> by;
        int i;
        int j;
        int k;
        int n;
        int m;
        int n2;
        int m2;
        amp::ampf<Precision> err;
        amp::ampf<Precision> t;
        amp::ampf<Precision> t1;
        amp::ampf<Precision> t2;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l1x;
        amp::ampf<Precision> l1y;
        amp::ampf<Precision> l1xy;
        amp::ampf<Precision> l2;
        amp::ampf<Precision> l2x;
        amp::ampf<Precision> l2y;
        amp::ampf<Precision> l2xy;
        amp::ampf<Precision> fm;
        amp::ampf<Precision> f1;
        amp::ampf<Precision> f2;
        amp::ampf<Precision> f3;
        amp::ampf<Precision> f4;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v1x;
        amp::ampf<Precision> v1y;
        amp::ampf<Precision> v1xy;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> v2x;
        amp::ampf<Precision> v2y;
        amp::ampf<Precision> v2xy;
        amp::ampf<Precision> mf;
        ap::template_1d_array< amp::ampf<Precision> > ra;
        ap::template_1d_array< amp::ampf<Precision> > ra2;
        int ralen;


        waserrors = false;
        passcount = 10;
        h = amp::ampf<Precision>("0.00001");
        lstep = amp::ampf<Precision>("0.001");
        blerrors = false;
        bcerrors = false;
        dserrors = false;
        cperrors = false;
        uperrors = false;
        lterrors = false;
        syerrors = false;
        rlerrors = false;
        rcerrors = false;
        
        //
        // Test: bilinear, bicubic
        //
        for(n=2; n<=7; n++)
        {
            for(m=2; m<=7; m++)
            {
                x.setbounds(0, n-1);
                y.setbounds(0, m-1);
                lx.setbounds(0, 2*n-2);
                ly.setbounds(0, 2*m-2);
                f.setbounds(0, m-1, 0, n-1);
                ft.setbounds(0, n-1, 0, m-1);
                for(pass=1; pass<=passcount; pass++)
                {
                    
                    //
                    // Prepare task:
                    // * X and Y stores grid
                    // * F stores function values
                    // * LX and LY stores twice dense grid (for Lipschitz testing)
                    //
                    ax = -1-amp::ampf<Precision>::getRandom();
                    bx = +1+amp::ampf<Precision>::getRandom();
                    ay = -1-amp::ampf<Precision>::getRandom();
                    by = +1+amp::ampf<Precision>::getRandom();
                    for(j=0; j<=n-1; j++)
                    {
                        x(j) = amp::ampf<Precision>("0.5")*(bx+ax)-amp::ampf<Precision>("0.5")*(bx-ax)*amp::cos<Precision>(amp::pi<Precision>()*(2*j+1)/(2*n));
                        if( j==0 )
                        {
                            x(j) = ax;
                        }
                        if( j==n-1 )
                        {
                            x(j) = bx;
                        }
                        lx(2*j) = x(j);
                        if( j>0 )
                        {
                            lx(2*j-1) = amp::ampf<Precision>("0.5")*(x(j)+x(j-1));
                        }
                    }
                    for(j=0; j<=n-1; j++)
                    {
                        k = ap::randominteger(n);
                        if( k!=j )
                        {
                            t = x(j);
                            x(j) = x(k);
                            x(k) = t;
                        }
                    }
                    for(i=0; i<=m-1; i++)
                    {
                        y(i) = amp::ampf<Precision>("0.5")*(by+ay)-amp::ampf<Precision>("0.5")*(by-ay)*amp::cos<Precision>(amp::pi<Precision>()*(2*i+1)/(2*m));
                        if( i==0 )
                        {
                            y(i) = ay;
                        }
                        if( i==m-1 )
                        {
                            y(i) = by;
                        }
                        ly(2*i) = y(i);
                        if( i>0 )
                        {
                            ly(2*i-1) = amp::ampf<Precision>("0.5")*(y(i)+y(i-1));
                        }
                    }
                    for(i=0; i<=m-1; i++)
                    {
                        k = ap::randominteger(m);
                        if( k!=i )
                        {
                            t = y(i);
                            y(i) = y(k);
                            y(k) = t;
                        }
                    }
                    for(i=0; i<=m-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            f(i,j) = amp::exp<Precision>(amp::ampf<Precision>("0.6")*x(j))-amp::exp<Precision>(-amp::ampf<Precision>("0.3")*y(i)+amp::ampf<Precision>("0.08")*x(j))+2*amp::cos<Precision>(amp::pi<Precision>()*(x(j)+amp::ampf<Precision>("1.2")*y(i)))+amp::ampf<Precision>("0.1")*amp::cos<Precision>(20*x(j)+15*y(i));
                        }
                    }
                    
                    //
                    // Test bilinear interpolation:
                    // * interpolation at the nodes
                    // * linearity
                    // * continuity
                    // * differentiation in the inner points
                    //
                    spline2d::spline2dbuildbilinear<Precision>(x, y, f, m, n, c);
                    err = 0;
                    for(i=0; i<=m-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(f(i,j)-spline2d::spline2dcalc<Precision>(c, x(j), y(i))));
                        }
                    }
                    blerrors = blerrors || err>10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                    err = 0;
                    for(i=0; i<=m-2; i++)
                    {
                        for(j=0; j<=n-2; j++)
                        {
                            
                            //
                            // Test for linearity between grid points
                            // (test point - geometric center of the cell)
                            //
                            fm = spline2d::spline2dcalc<Precision>(c, lx(2*j+1), ly(2*i+1));
                            f1 = spline2d::spline2dcalc<Precision>(c, lx(2*j), ly(2*i));
                            f2 = spline2d::spline2dcalc<Precision>(c, lx(2*j+2), ly(2*i));
                            f3 = spline2d::spline2dcalc<Precision>(c, lx(2*j+2), ly(2*i+2));
                            f4 = spline2d::spline2dcalc<Precision>(c, lx(2*j), ly(2*i+2));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(amp::ampf<Precision>("0.25")*(f1+f2+f3+f4)-fm));
                        }
                    }
                    blerrors = blerrors || err>10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                    lconst<Precision>(c, lx, ly, m, n, lstep, l1, l1x, l1y, l1xy);
                    lconst<Precision>(c, lx, ly, m, n, lstep/3, l2, l2x, l2y, l2xy);
                    blerrors = blerrors || l2/l1>amp::ampf<Precision>("1.2");
                    err = 0;
                    for(i=0; i<=m-2; i++)
                    {
                        for(j=0; j<=n-2; j++)
                        {
                            spline2d::spline2ddiff<Precision>(c, lx(2*j+1), ly(2*i+1), v1, v1x, v1y, v1xy);
                            twodnumder<Precision>(c, lx(2*j+1), ly(2*i+1), h, v2, v2x, v2y, v2xy);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(v1-v2));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(v1x-v2x));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(v1y-v2y));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(v1xy-v2xy));
                        }
                    }
                    dserrors = dserrors || err>amp::ampf<Precision>("1.0E-3");
                    uperrors = uperrors || !testunpack<Precision>(c, lx, ly);
                    lterrors = lterrors || !testlintrans<Precision>(c, ax, bx, ay, by);
                    
                    //
                    // Test bicubic interpolation.
                    // * interpolation at the nodes
                    // * smoothness
                    // * differentiation
                    //
                    spline2d::spline2dbuildbicubic<Precision>(x, y, f, m, n, c);
                    err = 0;
                    for(i=0; i<=m-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(f(i,j)-spline2d::spline2dcalc<Precision>(c, x(j), y(i))));
                        }
                    }
                    bcerrors = bcerrors || err>10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                    lconst<Precision>(c, lx, ly, m, n, lstep, l1, l1x, l1y, l1xy);
                    lconst<Precision>(c, lx, ly, m, n, lstep/3, l2, l2x, l2y, l2xy);
                    bcerrors = bcerrors || l2/l1>amp::ampf<Precision>("1.2");
                    bcerrors = bcerrors || l2x/l1x>amp::ampf<Precision>("1.2");
                    bcerrors = bcerrors || l2y/l1y>amp::ampf<Precision>("1.2");
                    if( l2xy>amp::ampf<Precision>("0.01") && l1xy>amp::ampf<Precision>("0.01") )
                    {
                        
                        //
                        // Cross-derivative continuity is tested only when
                        // bigger than 0.01. When the task size is too
                        // small, the d2F/dXdY is nearly zero and Lipschitz
                        // constant ratio is ill-conditioned.
                        //
                        bcerrors = bcerrors || l2xy/l1xy>amp::ampf<Precision>("1.2");
                    }
                    err = 0;
                    for(i=0; i<=2*m-2; i++)
                    {
                        for(j=0; j<=2*n-2; j++)
                        {
                            spline2d::spline2ddiff<Precision>(c, lx(j), ly(i), v1, v1x, v1y, v1xy);
                            twodnumder<Precision>(c, lx(j), ly(i), h, v2, v2x, v2y, v2xy);
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(v1-v2));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(v1x-v2x));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(v1y-v2y));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(v1xy-v2xy));
                        }
                    }
                    dserrors = dserrors || err>amp::ampf<Precision>("1.0E-3");
                    uperrors = uperrors || !testunpack<Precision>(c, lx, ly);
                    lterrors = lterrors || !testlintrans<Precision>(c, ax, bx, ay, by);
                    
                    //
                    // Copy/Serialise test
                    //
                    if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
                    {
                        spline2d::spline2dbuildbicubic<Precision>(x, y, f, m, n, c);
                    }
                    else
                    {
                        spline2d::spline2dbuildbilinear<Precision>(x, y, f, m, n, c);
                    }
                    unsetspline2d<Precision>(c2);
                    spline2d::spline2dcopy<Precision>(c, c2);
                    err = 0;
                    for(i=1; i<=5; i++)
                    {
                        t1 = ax+(bx-ax)*amp::ampf<Precision>::getRandom();
                        t2 = ay+(by-ay)*amp::ampf<Precision>::getRandom();
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(spline2d::spline2dcalc<Precision>(c, t1, t2)-spline2d::spline2dcalc<Precision>(c2, t1, t2)));
                    }
                    cperrors = cperrors || err>10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                    unsetspline2d<Precision>(c2);
                    spline2d::spline2dserialize<Precision>(c, ra, ralen);
                    ra2.setlength(ralen);
                    amp::vmove(ra2.getvector(0, ralen-1), ra.getvector(0, ralen-1));
                    spline2d::spline2dunserialize<Precision>(ra2, c2);
                    err = 0;
                    for(i=1; i<=5; i++)
                    {
                        t1 = ax+(bx-ax)*amp::ampf<Precision>::getRandom();
                        t2 = ay+(by-ay)*amp::ampf<Precision>::getRandom();
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(spline2d::spline2dcalc<Precision>(c, t1, t2)-spline2d::spline2dcalc<Precision>(c2, t1, t2)));
                    }
                    cperrors = cperrors || err>10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                    
                    //
                    // Special symmetry test
                    //
                    err = 0;
                    for(jobtype=0; jobtype<=1; jobtype++)
                    {
                        
                        //
                        // Prepare
                        //
                        for(i=0; i<=m-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                ft(j,i) = f(i,j);
                            }
                        }
                        if( jobtype==0 )
                        {
                            spline2d::spline2dbuildbilinear<Precision>(x, y, f, m, n, c);
                            spline2d::spline2dbuildbilinear<Precision>(y, x, ft, n, m, c2);
                        }
                        else
                        {
                            spline2d::spline2dbuildbicubic<Precision>(x, y, f, m, n, c);
                            spline2d::spline2dbuildbicubic<Precision>(y, x, ft, n, m, c2);
                        }
                        
                        //
                        // Test
                        //
                        for(i=1; i<=10; i++)
                        {
                            t1 = ax+(bx-ax)*amp::ampf<Precision>::getRandom();
                            t2 = ay+(by-ay)*amp::ampf<Precision>::getRandom();
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(spline2d::spline2dcalc<Precision>(c, t1, t2)-spline2d::spline2dcalc<Precision>(c2, t2, t1)));
                        }
                    }
                    syerrors = syerrors || err>10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                }
            }
        }
        
        //
        // Test resample
        //
        for(m=2; m<=6; m++)
        {
            for(n=2; n<=6; n++)
            {
                f.setbounds(0, m-1, 0, n-1);
                x.setbounds(0, n-1);
                y.setbounds(0, m-1);
                for(j=0; j<=n-1; j++)
                {
                    x(j) = amp::ampf<Precision>(j)/(amp::ampf<Precision>(n-1));
                }
                for(i=0; i<=m-1; i++)
                {
                    y(i) = amp::ampf<Precision>(i)/(amp::ampf<Precision>(m-1));
                }
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        f(i,j) = amp::exp<Precision>(amp::ampf<Precision>("0.6")*x(j))-amp::exp<Precision>(-amp::ampf<Precision>("0.3")*y(i)+amp::ampf<Precision>("0.08")*x(j))+2*amp::cos<Precision>(amp::pi<Precision>()*(x(j)+amp::ampf<Precision>("1.2")*y(i)))+amp::ampf<Precision>("0.1")*amp::cos<Precision>(20*x(j)+15*y(i));
                    }
                }
                for(m2=2; m2<=6; m2++)
                {
                    for(n2=2; n2<=6; n2++)
                    {
                        for(pass=1; pass<=passcount; pass++)
                        {
                            for(jobtype=0; jobtype<=1; jobtype++)
                            {
                                if( jobtype==0 )
                                {
                                    spline2d::spline2dresamplebilinear<Precision>(f, m, n, fr, m2, n2);
                                    spline2d::spline2dbuildbilinear<Precision>(x, y, f, m, n, c);
                                }
                                if( jobtype==1 )
                                {
                                    spline2d::spline2dresamplebicubic<Precision>(f, m, n, fr, m2, n2);
                                    spline2d::spline2dbuildbicubic<Precision>(x, y, f, m, n, c);
                                }
                                err = 0;
                                mf = 0;
                                for(i=0; i<=m2-1; i++)
                                {
                                    for(j=0; j<=n2-1; j++)
                                    {
                                        v1 = spline2d::spline2dcalc<Precision>(c, amp::ampf<Precision>(j)/(amp::ampf<Precision>(n2-1)), amp::ampf<Precision>(i)/(amp::ampf<Precision>(m2-1)));
                                        v2 = fr(i,j);
                                        err = amp::maximum<Precision>(err, amp::abs<Precision>(v1-v2));
                                        mf = amp::maximum<Precision>(mf, amp::abs<Precision>(v1));
                                    }
                                }
                                if( jobtype==0 )
                                {
                                    rlerrors = rlerrors || err/mf>10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                                }
                                if( jobtype==1 )
                                {
                                    rcerrors = rcerrors || err/mf>10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //
        // report
        //
        waserrors = blerrors || bcerrors || dserrors || cperrors || uperrors || lterrors || syerrors || rlerrors || rcerrors;
        if( !silent )
        {
            printf("TESTING 2D INTERPOLATION\n");
            
            //
            // Normal tests
            //
            printf("BILINEAR TEST:                           ");
            if( blerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("BICUBIC TEST:                            ");
            if( bcerrors )
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
            printf("COPY/SERIALIZE TEST:                     ");
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
            printf("SPECIAL SYMMETRY TEST:                   ");
            if( syerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("BILINEAR RESAMPLING TEST:                ");
            if( rlerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("BICUBIC RESAMPLING TEST:                 ");
            if( rcerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            
            //
            // Summary
            //
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
    void lconst(const spline2d::spline2dinterpolant<Precision>& c,
        const ap::template_1d_array< amp::ampf<Precision> >& lx,
        const ap::template_1d_array< amp::ampf<Precision> >& ly,
        int m,
        int n,
        amp::ampf<Precision> lstep,
        amp::ampf<Precision>& lc,
        amp::ampf<Precision>& lcx,
        amp::ampf<Precision>& lcy,
        amp::ampf<Precision>& lcxy)
    {
        int i;
        int j;
        amp::ampf<Precision> f1;
        amp::ampf<Precision> f2;
        amp::ampf<Precision> f3;
        amp::ampf<Precision> f4;
        amp::ampf<Precision> fx1;
        amp::ampf<Precision> fx2;
        amp::ampf<Precision> fx3;
        amp::ampf<Precision> fx4;
        amp::ampf<Precision> fy1;
        amp::ampf<Precision> fy2;
        amp::ampf<Precision> fy3;
        amp::ampf<Precision> fy4;
        amp::ampf<Precision> fxy1;
        amp::ampf<Precision> fxy2;
        amp::ampf<Precision> fxy3;
        amp::ampf<Precision> fxy4;
        amp::ampf<Precision> s2lstep;


        lc = 0;
        lcx = 0;
        lcy = 0;
        lcxy = 0;
        s2lstep = amp::sqrt<Precision>(amp::ampf<Precision>(2))*lstep;
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Calculate
                //
                twodnumder<Precision>(c, lx(j)-lstep/2, ly(i)-lstep/2, lstep/4, f1, fx1, fy1, fxy1);
                twodnumder<Precision>(c, lx(j)+lstep/2, ly(i)-lstep/2, lstep/4, f2, fx2, fy2, fxy2);
                twodnumder<Precision>(c, lx(j)+lstep/2, ly(i)+lstep/2, lstep/4, f3, fx3, fy3, fxy3);
                twodnumder<Precision>(c, lx(j)-lstep/2, ly(i)+lstep/2, lstep/4, f4, fx4, fy4, fxy4);
                
                //
                // Lipschitz constant for the function itself
                //
                lc = amp::maximum<Precision>(lc, amp::abs<Precision>((f1-f2)/lstep));
                lc = amp::maximum<Precision>(lc, amp::abs<Precision>((f2-f3)/lstep));
                lc = amp::maximum<Precision>(lc, amp::abs<Precision>((f3-f4)/lstep));
                lc = amp::maximum<Precision>(lc, amp::abs<Precision>((f4-f1)/lstep));
                lc = amp::maximum<Precision>(lc, amp::abs<Precision>((f1-f3)/s2lstep));
                lc = amp::maximum<Precision>(lc, amp::abs<Precision>((f2-f4)/s2lstep));
                
                //
                // Lipschitz constant for the first derivative
                //
                lcx = amp::maximum<Precision>(lcx, amp::abs<Precision>((fx1-fx2)/lstep));
                lcx = amp::maximum<Precision>(lcx, amp::abs<Precision>((fx2-fx3)/lstep));
                lcx = amp::maximum<Precision>(lcx, amp::abs<Precision>((fx3-fx4)/lstep));
                lcx = amp::maximum<Precision>(lcx, amp::abs<Precision>((fx4-fx1)/lstep));
                lcx = amp::maximum<Precision>(lcx, amp::abs<Precision>((fx1-fx3)/s2lstep));
                lcx = amp::maximum<Precision>(lcx, amp::abs<Precision>((fx2-fx4)/s2lstep));
                
                //
                // Lipschitz constant for the first derivative
                //
                lcy = amp::maximum<Precision>(lcy, amp::abs<Precision>((fy1-fy2)/lstep));
                lcy = amp::maximum<Precision>(lcy, amp::abs<Precision>((fy2-fy3)/lstep));
                lcy = amp::maximum<Precision>(lcy, amp::abs<Precision>((fy3-fy4)/lstep));
                lcy = amp::maximum<Precision>(lcy, amp::abs<Precision>((fy4-fy1)/lstep));
                lcy = amp::maximum<Precision>(lcy, amp::abs<Precision>((fy1-fy3)/s2lstep));
                lcy = amp::maximum<Precision>(lcy, amp::abs<Precision>((fy2-fy4)/s2lstep));
                
                //
                // Lipschitz constant for the cross-derivative
                //
                lcxy = amp::maximum<Precision>(lcxy, amp::abs<Precision>((fxy1-fxy2)/lstep));
                lcxy = amp::maximum<Precision>(lcxy, amp::abs<Precision>((fxy2-fxy3)/lstep));
                lcxy = amp::maximum<Precision>(lcxy, amp::abs<Precision>((fxy3-fxy4)/lstep));
                lcxy = amp::maximum<Precision>(lcxy, amp::abs<Precision>((fxy4-fxy1)/lstep));
                lcxy = amp::maximum<Precision>(lcxy, amp::abs<Precision>((fxy1-fxy3)/s2lstep));
                lcxy = amp::maximum<Precision>(lcxy, amp::abs<Precision>((fxy2-fxy4)/s2lstep));
            }
        }
    }


    /*************************************************************************
    Numerical differentiation.
    *************************************************************************/
    template<unsigned int Precision>
    void twodnumder(const spline2d::spline2dinterpolant<Precision>& c,
        amp::ampf<Precision> x,
        amp::ampf<Precision> y,
        amp::ampf<Precision> h,
        amp::ampf<Precision>& f,
        amp::ampf<Precision>& fx,
        amp::ampf<Precision>& fy,
        amp::ampf<Precision>& fxy)
    {
        f = spline2d::spline2dcalc<Precision>(c, x, y);
        fx = (spline2d::spline2dcalc<Precision>(c, x+h, y)-spline2d::spline2dcalc<Precision>(c, x-h, y))/(2*h);
        fy = (spline2d::spline2dcalc<Precision>(c, x, y+h)-spline2d::spline2dcalc<Precision>(c, x, y-h))/(2*h);
        fxy = (spline2d::spline2dcalc<Precision>(c, x+h, y+h)-spline2d::spline2dcalc<Precision>(c, x-h, y+h)-spline2d::spline2dcalc<Precision>(c, x+h, y-h)+spline2d::spline2dcalc<Precision>(c, x-h, y-h))/amp::sqr<Precision>(2*h);
    }


    /*************************************************************************
    Unpack test
    *************************************************************************/
    template<unsigned int Precision>
    bool testunpack(const spline2d::spline2dinterpolant<Precision>& c,
        const ap::template_1d_array< amp::ampf<Precision> >& lx,
        const ap::template_1d_array< amp::ampf<Precision> >& ly)
    {
        bool result;
        int i;
        int j;
        int n;
        int m;
        int ci;
        int cj;
        int p;
        amp::ampf<Precision> err;
        amp::ampf<Precision> tx;
        amp::ampf<Precision> ty;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        int pass;
        int passcount;
        ap::template_2d_array< amp::ampf<Precision> > tbl;


        passcount = 20;
        err = 0;
        spline2d::spline2dunpack<Precision>(c, m, n, tbl);
        for(i=0; i<=m-2; i++)
        {
            for(j=0; j<=n-2; j++)
            {
                for(pass=1; pass<=passcount; pass++)
                {
                    p = (n-1)*i+j;
                    tx = (amp::ampf<Precision>("0.001")+amp::ampf<Precision>("0.999")*amp::ampf<Precision>::getRandom())*(tbl(p,1)-tbl(p,0));
                    ty = (amp::ampf<Precision>("0.001")+amp::ampf<Precision>("0.999")*amp::ampf<Precision>::getRandom())*(tbl(p,3)-tbl(p,2));
                    
                    //
                    // Interpolation properties
                    //
                    v1 = 0;
                    for(ci=0; ci<=3; ci++)
                    {
                        for(cj=0; cj<=3; cj++)
                        {
                            v1 = v1+tbl(p,4+ci*4+cj)*amp::pow<Precision>(tx, amp::ampf<Precision>(ci))*amp::pow<Precision>(ty, amp::ampf<Precision>(cj));
                        }
                    }
                    v2 = spline2d::spline2dcalc<Precision>(c, tbl(p,0)+tx, tbl(p,2)+ty);
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(v1-v2));
                    
                    //
                    // Grid correctness
                    //
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(lx(2*j)-tbl(p,0)));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(lx(2*(j+1))-tbl(p,1)));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(ly(2*i)-tbl(p,2)));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(ly(2*(i+1))-tbl(p,3)));
                }
            }
        }
        result = err<10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        return result;
    }


    /*************************************************************************
    LinTrans test
    *************************************************************************/
    template<unsigned int Precision>
    bool testlintrans(const spline2d::spline2dinterpolant<Precision>& c,
        amp::ampf<Precision> ax,
        amp::ampf<Precision> bx,
        amp::ampf<Precision> ay,
        amp::ampf<Precision> by)
    {
        bool result;
        amp::ampf<Precision> err;
        amp::ampf<Precision> a1;
        amp::ampf<Precision> a2;
        amp::ampf<Precision> b1;
        amp::ampf<Precision> b2;
        amp::ampf<Precision> tx;
        amp::ampf<Precision> ty;
        amp::ampf<Precision> vx;
        amp::ampf<Precision> vy;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        int pass;
        int passcount;
        int xjob;
        int yjob;
        spline2d::spline2dinterpolant<Precision> c2;


        passcount = 5;
        err = 0;
        for(xjob=0; xjob<=1; xjob++)
        {
            for(yjob=0; yjob<=1; yjob++)
            {
                for(pass=1; pass<=passcount; pass++)
                {
                    
                    //
                    // Prepare
                    //
                    do
                    {
                        a1 = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    while( a1==0 );
                    a1 = a1*xjob;
                    b1 = 2*amp::ampf<Precision>::getRandom()-1;
                    do
                    {
                        a2 = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    while( a2==0 );
                    a2 = a2*yjob;
                    b2 = 2*amp::ampf<Precision>::getRandom()-1;
                    
                    //
                    // Test XY
                    //
                    spline2d::spline2dcopy<Precision>(c, c2);
                    spline2d::spline2dlintransxy<Precision>(c2, a1, b1, a2, b2);
                    tx = ax+amp::ampf<Precision>::getRandom()*(bx-ax);
                    ty = ay+amp::ampf<Precision>::getRandom()*(by-ay);
                    if( xjob==0 )
                    {
                        tx = b1;
                        vx = ax+amp::ampf<Precision>::getRandom()*(bx-ax);
                    }
                    else
                    {
                        vx = (tx-b1)/a1;
                    }
                    if( yjob==0 )
                    {
                        ty = b2;
                        vy = ay+amp::ampf<Precision>::getRandom()*(by-ay);
                    }
                    else
                    {
                        vy = (ty-b2)/a2;
                    }
                    v1 = spline2d::spline2dcalc<Precision>(c, tx, ty);
                    v2 = spline2d::spline2dcalc<Precision>(c2, vx, vy);
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(v1-v2));
                    
                    //
                    // Test F
                    //
                    spline2d::spline2dcopy<Precision>(c, c2);
                    spline2d::spline2dlintransf<Precision>(c2, a1, b1);
                    tx = ax+amp::ampf<Precision>::getRandom()*(bx-ax);
                    ty = ay+amp::ampf<Precision>::getRandom()*(by-ay);
                    v1 = spline2d::spline2dcalc<Precision>(c, tx, ty);
                    v2 = spline2d::spline2dcalc<Precision>(c2, tx, ty);
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(a1*v1+b1-v2));
                }
            }
        }
        result = err<10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        return result;
    }


    /*************************************************************************
    Unset spline, i.e. initialize it with random garbage
    *************************************************************************/
    template<unsigned int Precision>
    void unsetspline2d(spline2d::spline2dinterpolant<Precision>& c)
    {
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_2d_array< amp::ampf<Precision> > f;


        x.setlength(2);
        y.setlength(2);
        f.setlength(2, 2);
        x(0) = -1;
        x(1) = +1;
        y(0) = -1;
        y(1) = +1;
        f(0,0) = 0;
        f(0,1) = 0;
        f(1,0) = 0;
        f(1,1) = 0;
        spline2d::spline2dbuildbilinear<Precision>(x, y, f, 2, 2, c);
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testspline2dunit_test_silent()
    {
        bool result;


        result = test2dinterpolation<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testspline2dunit_test()
    {
        bool result;


        result = test2dinterpolation<Precision>(false);
        return result;
    }
} // namespace

#endif
