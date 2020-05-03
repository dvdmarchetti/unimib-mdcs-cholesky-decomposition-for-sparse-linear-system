
#ifndef _llstestunit_h
#define _llstestunit_h

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
#include "lsfit.h"
namespace llstestunit
{
    template<unsigned int Precision>
    bool testlls(bool silent);
    template<unsigned int Precision>
    bool isglssolution(int n,
        int m,
        int k,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        const ap::template_2d_array< amp::ampf<Precision> >& cmatrix,
        ap::template_1d_array< amp::ampf<Precision> > c);
    template<unsigned int Precision>
    amp::ampf<Precision> getglserror(int n,
        int m,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        const ap::template_1d_array< amp::ampf<Precision> >& c);
    template<unsigned int Precision>
    void fitlinearnonlinear(int m,
        bool gradonly,
        const ap::template_2d_array< amp::ampf<Precision> >& xy,
        lsfit::lsfitstate<Precision>& state,
        bool& nlserrors);
    template<unsigned int Precision>
    bool llstestunit_test_silent();
    template<unsigned int Precision>
    bool llstestunit_test();


    template<unsigned int Precision>
    bool testlls(bool silent)
    {
        bool result;
        bool waserrors;
        bool llserrors;
        bool nlserrors;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> nlthreshold;
        int maxn;
        int maxm;
        int passcount;
        int n;
        int m;
        int i;
        int j;
        int k;
        int pass;
        amp::ampf<Precision> xscale;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > w2;
        ap::template_1d_array< amp::ampf<Precision> > c;
        ap::template_1d_array< amp::ampf<Precision> > c2;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        ap::template_2d_array< amp::ampf<Precision> > cm;
        amp::ampf<Precision> v;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        lsfit::lsfitreport<Precision> rep;
        lsfit::lsfitreport<Precision> rep2;
        int info;
        int info2;
        amp::ampf<Precision> refrms;
        amp::ampf<Precision> refavg;
        amp::ampf<Precision> refavgrel;
        amp::ampf<Precision> refmax;
        lsfit::lsfitstate<Precision> state;


        waserrors = false;
        llserrors = false;
        nlserrors = false;
        threshold = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        nlthreshold = amp::ampf<Precision>("0.00001");
        maxn = 6;
        maxm = 6;
        passcount = 4;
        
        //
        // Testing unconstrained least squares (linear/nonlinear)
        //
        for(n=1; n<=maxn; n++)
        {
            for(m=1; m<=maxm; m++)
            {
                for(pass=1; pass<=passcount; pass++)
                {
                    
                    //
                    // Solve non-degenerate linear least squares task
                    // Use Chebyshev basis. Its condition number is very good.
                    //
                    a.setlength(n, m);
                    x.setlength(n);
                    y.setlength(n);
                    w.setlength(n);
                    xscale = amp::ampf<Precision>("0.9")+amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getRandom();
                    for(i=0; i<=n-1; i++)
                    {
                        if( n==1 )
                        {
                            x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                        else
                        {
                            x(i) = xscale*(amp::ampf<Precision>(2*i)/(amp::ampf<Precision>(n-1))-1);
                        }
                        y(i) = 3*x(i)+amp::exp<Precision>(x(i));
                        w(i) = 1+amp::ampf<Precision>::getRandom();
                        a(i,0) = 1;
                        if( m>1 )
                        {
                            a(i,1) = x(i);
                        }
                        for(j=2; j<=m-1; j++)
                        {
                            a(i,j) = 2*x(i)*a(i,j-1)-a(i,j-2);
                        }
                    }
                    
                    //
                    // 1. test weighted fitting (optimality)
                    // 2. Solve degenerate least squares task built on the basis
                    //    of previous task
                    //
                    lsfit::lsfitlinearw<Precision>(y, w, a, n, m, info, c, rep);
                    if( info<=0 )
                    {
                        llserrors = true;
                    }
                    else
                    {
                        llserrors = llserrors || !isglssolution<Precision>(n, m, 0, y, w, a, cm, c);
                    }
                    a2.setlength(n, 2*m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            a2(i,2*j+0) = a(i,j);
                            a2(i,2*j+1) = a(i,j);
                        }
                    }
                    lsfit::lsfitlinearw<Precision>(y, w, a2, n, 2*m, info, c2, rep);
                    if( info<=0 )
                    {
                        llserrors = true;
                    }
                    else
                    {
                        
                        //
                        // test answer correctness using design matrix properties
                        // and previous task solution
                        //
                        for(j=0; j<=m-1; j++)
                        {
                            llserrors = llserrors || amp::abs<Precision>(c2(2*j+0)+c2(2*j+1)-c(j))>threshold;
                        }
                    }
                    
                    //
                    // test non-weighted fitting
                    //
                    w2.setlength(n);
                    for(i=0; i<=n-1; i++)
                    {
                        w2(i) = 1;
                    }
                    lsfit::lsfitlinearw<Precision>(y, w2, a, n, m, info, c, rep);
                    lsfit::lsfitlinear<Precision>(y, a, n, m, info2, c2, rep2);
                    if( info<=0 || info2<=0 )
                    {
                        llserrors = true;
                    }
                    else
                    {
                        
                        //
                        // test answer correctness
                        //
                        for(j=0; j<=m-1; j++)
                        {
                            llserrors = llserrors || amp::abs<Precision>(c(j)-c2(j))>threshold;
                        }
                        llserrors = llserrors || amp::abs<Precision>(rep.taskrcond-rep2.taskrcond)>threshold;
                    }
                    
                    //
                    // test nonlinear fitting on the linear task
                    // (only non-degenerate task are tested)
                    // and compare with answer from linear fitting subroutine
                    //
                    if( n>=m )
                    {
                        c2.setlength(m);
                        
                        //
                        // test gradient-only or Hessian-based weighted fitting
                        //
                        lsfit::lsfitlinearw<Precision>(y, w, a, n, m, info, c, rep);
                        for(i=0; i<=m-1; i++)
                        {
                            c2(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                        lsfit::lsfitnonlinearwfg<Precision>(a, y, w, c2, n, m, m, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), state);
                        lsfit::lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), nlthreshold, 0);
                        fitlinearnonlinear<Precision>(m, true, a, state, nlserrors);
                        lsfit::lsfitnonlinearresults<Precision>(state, info, c2, rep2);
                        if( info<=0 )
                        {
                            nlserrors = true;
                        }
                        else
                        {
                            for(i=0; i<=m-1; i++)
                            {
                                nlserrors = nlserrors || amp::abs<Precision>(c(i)-c2(i))>100*nlthreshold;
                            }
                        }
                        for(i=0; i<=m-1; i++)
                        {
                            c2(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                        lsfit::lsfitnonlinearwfgh<Precision>(a, y, w, c2, n, m, m, state);
                        lsfit::lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), nlthreshold, 0);
                        fitlinearnonlinear<Precision>(m, false, a, state, nlserrors);
                        lsfit::lsfitnonlinearresults<Precision>(state, info, c2, rep2);
                        if( info<=0 )
                        {
                            nlserrors = true;
                        }
                        else
                        {
                            for(i=0; i<=m-1; i++)
                            {
                                nlserrors = nlserrors || amp::abs<Precision>(c(i)-c2(i))>100*nlthreshold;
                            }
                        }
                        
                        //
                        // test gradient-only or Hessian-based fitting without weights
                        //
                        lsfit::lsfitlinear<Precision>(y, a, n, m, info, c, rep);
                        for(i=0; i<=m-1; i++)
                        {
                            c2(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                        lsfit::lsfitnonlinearfg<Precision>(a, y, c2, n, m, m, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), state);
                        lsfit::lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), nlthreshold, 0);
                        fitlinearnonlinear<Precision>(m, true, a, state, nlserrors);
                        lsfit::lsfitnonlinearresults<Precision>(state, info, c2, rep2);
                        if( info<=0 )
                        {
                            nlserrors = true;
                        }
                        else
                        {
                            for(i=0; i<=m-1; i++)
                            {
                                nlserrors = nlserrors || amp::abs<Precision>(c(i)-c2(i))>100*nlthreshold;
                            }
                        }
                        for(i=0; i<=m-1; i++)
                        {
                            c2(i) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                        lsfit::lsfitnonlinearfgh<Precision>(a, y, c2, n, m, m, state);
                        lsfit::lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), nlthreshold, 0);
                        fitlinearnonlinear<Precision>(m, false, a, state, nlserrors);
                        lsfit::lsfitnonlinearresults<Precision>(state, info, c2, rep2);
                        if( info<=0 )
                        {
                            nlserrors = true;
                        }
                        else
                        {
                            for(i=0; i<=m-1; i++)
                            {
                                nlserrors = nlserrors || amp::abs<Precision>(c(i)-c2(i))>100*nlthreshold;
                            }
                        }
                    }
                }
            }
            
            //
            // test correctness of the RCond field
            //
            a.setbounds(0, n-1, 0, n-1);
            x.setbounds(0, n-1);
            y.setbounds(0, n-1);
            w.setbounds(0, n-1);
            v1 = amp::ampf<Precision>::getAlgoPascalMaxNumber();
            v2 = amp::ampf<Precision>::getAlgoPascalMinNumber();
            for(i=0; i<=n-1; i++)
            {
                x(i) = amp::ampf<Precision>("0.1")+amp::ampf<Precision>("0.9")*amp::ampf<Precision>::getRandom();
                y(i) = amp::ampf<Precision>("0.1")+amp::ampf<Precision>("0.9")*amp::ampf<Precision>::getRandom();
                w(i) = 1;
                for(j=0; j<=n-1; j++)
                {
                    if( i==j )
                    {
                        a(i,i) = amp::ampf<Precision>("0.1")+amp::ampf<Precision>("0.9")*amp::ampf<Precision>::getRandom();
                        v1 = amp::minimum<Precision>(v1, a(i,i));
                        v2 = amp::maximum<Precision>(v2, a(i,i));
                    }
                    else
                    {
                        a(i,j) = 0;
                    }
                }
            }
            lsfit::lsfitlinearw<Precision>(y, w, a, n, n, info, c, rep);
            if( info<=0 )
            {
                llserrors = true;
            }
            else
            {
                llserrors = llserrors || amp::abs<Precision>(rep.taskrcond-v1/v2)>threshold;
            }
        }
        
        //
        // Test constrained least squares
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                for(m=1; m<=maxm; m++)
                {
                    
                    //
                    // test for K<>0
                    //
                    for(k=1; k<=m-1; k++)
                    {
                        
                        //
                        // Prepare Chebyshev basis. Its condition number is very good.
                        // Prepare constraints (random numbers)
                        //
                        a.setlength(n, m);
                        x.setlength(n);
                        y.setlength(n);
                        w.setlength(n);
                        xscale = amp::ampf<Precision>("0.9")+amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getRandom();
                        for(i=0; i<=n-1; i++)
                        {
                            if( n==1 )
                            {
                                x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                            else
                            {
                                x(i) = xscale*(amp::ampf<Precision>(2*i)/(amp::ampf<Precision>(n-1))-1);
                            }
                            y(i) = 3*x(i)+amp::exp<Precision>(x(i));
                            w(i) = 1+amp::ampf<Precision>::getRandom();
                            a(i,0) = 1;
                            if( m>1 )
                            {
                                a(i,1) = x(i);
                            }
                            for(j=2; j<=m-1; j++)
                            {
                                a(i,j) = 2*x(i)*a(i,j-1)-a(i,j-2);
                            }
                        }
                        cm.setlength(k, m+1);
                        for(i=0; i<=k-1; i++)
                        {
                            for(j=0; j<=m; j++)
                            {
                                cm(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        
                        //
                        // Solve constrained task
                        //
                        lsfit::lsfitlinearwc<Precision>(y, w, a, cm, n, m, k, info, c, rep);
                        if( info<=0 )
                        {
                            llserrors = true;
                        }
                        else
                        {
                            llserrors = llserrors || !isglssolution<Precision>(n, m, k, y, w, a, cm, c);
                        }
                        
                        //
                        // test non-weighted fitting
                        //
                        w2.setlength(n);
                        for(i=0; i<=n-1; i++)
                        {
                            w2(i) = 1;
                        }
                        lsfit::lsfitlinearwc<Precision>(y, w2, a, cm, n, m, k, info, c, rep);
                        lsfit::lsfitlinearc<Precision>(y, a, cm, n, m, k, info2, c2, rep2);
                        if( info<=0 || info2<=0 )
                        {
                            llserrors = true;
                        }
                        else
                        {
                            
                            //
                            // test answer correctness
                            //
                            for(j=0; j<=m-1; j++)
                            {
                                llserrors = llserrors || amp::abs<Precision>(c(j)-c2(j))>threshold;
                            }
                            llserrors = llserrors || amp::abs<Precision>(rep.taskrcond-rep2.taskrcond)>threshold;
                        }
                    }
                }
            }
        }
        
        //
        // nonlinear task for nonlinear fitting:
        //
        //     f(X,C) = 1/(1+C*X^2),
        //     C(true) = 2.
        //
        n = 100;
        c.setlength(1);
        c(0) = 1+2*amp::ampf<Precision>::getRandom();
        a.setlength(n, 1);
        y.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            a(i,0) = 4*amp::ampf<Precision>::getRandom()-2;
            y(i) = 1/(1+2*amp::sqr<Precision>(a(i,0)));
        }
        lsfit::lsfitnonlinearfg<Precision>(a, y, c, n, 1, 1, true, state);
        lsfit::lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), nlthreshold, 0);
        while( lsfit::lsfitnonlineariteration<Precision>(state) )
        {
            if( state.needf )
            {
                state.f = 1/(1+state.c(0)*amp::sqr<Precision>(state.x(0)));
            }
            if( state.needfg )
            {
                state.f = 1/(1+state.c(0)*amp::sqr<Precision>(state.x(0)));
                state.g(0) = -amp::sqr<Precision>(state.x(0))/amp::sqr<Precision>(1+state.c(0)*amp::sqr<Precision>(state.x(0)));
            }
        }
        lsfit::lsfitnonlinearresults<Precision>(state, info, c, rep);
        if( info<=0 )
        {
            nlserrors = true;
        }
        else
        {
            nlserrors = nlserrors || amp::abs<Precision>(c(0)-2)>100*nlthreshold;
        }
        
        //
        // solve simple task (fitting by constant function) and check
        // correctness of the errors calculated by subroutines
        //
        for(pass=1; pass<=passcount; pass++)
        {
            
            //
            // test on task with non-zero Yi
            //
            n = 4;
            v1 = amp::ampf<Precision>::getRandom();
            v2 = amp::ampf<Precision>::getRandom();
            v = 1+amp::ampf<Precision>::getRandom();
            c.setlength(1);
            c(0) = 1+2*amp::ampf<Precision>::getRandom();
            a.setlength(4, 1);
            y.setlength(4);
            a(0,0) = 1;
            y(0) = v-v2;
            a(1,0) = 1;
            y(1) = v-v1;
            a(2,0) = 1;
            y(2) = v+v1;
            a(3,0) = 1;
            y(3) = v+v2;
            refrms = amp::sqrt<Precision>((amp::sqr<Precision>(v1)+amp::sqr<Precision>(v2))/2);
            refavg = (amp::abs<Precision>(v1)+amp::abs<Precision>(v2))/2;
            refavgrel = amp::ampf<Precision>("0.25")*(amp::abs<Precision>(v2)/amp::abs<Precision>(v-v2)+amp::abs<Precision>(v1)/amp::abs<Precision>(v-v1)+amp::abs<Precision>(v1)/amp::abs<Precision>(v+v1)+amp::abs<Precision>(v2)/amp::abs<Precision>(v+v2));
            refmax = amp::maximum<Precision>(v1, v2);
            
            //
            // Test LLS
            //
            lsfit::lsfitlinear<Precision>(y, a, 4, 1, info, c, rep);
            if( info<=0 )
            {
                llserrors = true;
            }
            else
            {
                llserrors = llserrors || amp::abs<Precision>(c(0)-v)>threshold;
                llserrors = llserrors || amp::abs<Precision>(rep.rmserror-refrms)>threshold;
                llserrors = llserrors || amp::abs<Precision>(rep.avgerror-refavg)>threshold;
                llserrors = llserrors || amp::abs<Precision>(rep.avgrelerror-refavgrel)>threshold;
                llserrors = llserrors || amp::abs<Precision>(rep.maxerror-refmax)>threshold;
            }
            
            //
            // Test NLS
            //
            lsfit::lsfitnonlinearfg<Precision>(a, y, c, 4, 1, 1, true, state);
            lsfit::lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), nlthreshold, 0);
            while( lsfit::lsfitnonlineariteration<Precision>(state) )
            {
                if( state.needf )
                {
                    state.f = state.c(0);
                }
                if( state.needfg )
                {
                    state.f = state.c(0);
                    state.g(0) = 1;
                }
            }
            lsfit::lsfitnonlinearresults<Precision>(state, info, c, rep);
            if( info<=0 )
            {
                nlserrors = true;
            }
            else
            {
                nlserrors = nlserrors || amp::abs<Precision>(c(0)-v)>threshold;
                nlserrors = nlserrors || amp::abs<Precision>(rep.rmserror-refrms)>threshold;
                nlserrors = nlserrors || amp::abs<Precision>(rep.avgerror-refavg)>threshold;
                nlserrors = nlserrors || amp::abs<Precision>(rep.avgrelerror-refavgrel)>threshold;
                nlserrors = nlserrors || amp::abs<Precision>(rep.maxerror-refmax)>threshold;
            }
        }
        
        //
        // report
        //
        waserrors = llserrors || nlserrors;
        if( !silent )
        {
            printf("TESTING LEAST SQUARES\n");
            printf("LINEAR LEAST SQUARES:                    ");
            if( llserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("NON-LINEAR LEAST SQUARES:                ");
            if( nlserrors )
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
    Tests whether C is solution of (possibly) constrained LLS problem
    *************************************************************************/
    template<unsigned int Precision>
    bool isglssolution(int n,
        int m,
        int k,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        const ap::template_2d_array< amp::ampf<Precision> >& cmatrix,
        ap::template_1d_array< amp::ampf<Precision> > c)
    {
        bool result;
        int i;
        int j;
        ap::template_1d_array< amp::ampf<Precision> > c2;
        ap::template_1d_array< amp::ampf<Precision> > sv;
        ap::template_1d_array< amp::ampf<Precision> > deltac;
        ap::template_1d_array< amp::ampf<Precision> > deltaproj;
        ap::template_2d_array< amp::ampf<Precision> > u;
        ap::template_2d_array< amp::ampf<Precision> > vt;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s1;
        amp::ampf<Precision> s2;
        amp::ampf<Precision> s3;
        amp::ampf<Precision> delta;
        amp::ampf<Precision> threshold;


        
        //
        // Setup.
        // Threshold is small because CMatrix may be ill-conditioned
        //
        delta = amp::ampf<Precision>("0.001");
        threshold = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalEpsilon());
        c2.setlength(m);
        deltac.setlength(m);
        deltaproj.setlength(m);
        
        //
        // test whether C is feasible point or not (projC must be close to C)
        //
        for(i=0; i<=k-1; i++)
        {
            v = amp::vdotproduct(cmatrix.getrow(i, 0, m-1), c.getvector(0, m-1));
            if( amp::abs<Precision>(v-cmatrix(i,m))>threshold )
            {
                result = false;
                return result;
            }
        }
        
        //
        // find orthogonal basis of Null(CMatrix) (stored in rows from K to M-1)
        //
        if( k>0 )
        {
            svd::rmatrixsvd<Precision>(cmatrix, k, m, 0, 2, 2, sv, u, vt);
        }
        
        //
        // Test result
        //
        result = true;
        s1 = getglserror<Precision>(n, m, y, w, fmatrix, c);
        for(j=0; j<=m-1; j++)
        {
            
            //
            // prepare modification of C which leave us in the feasible set.
            //
            // let deltaC be increment on Jth coordinate, then project
            // deltaC in the Null(CMatrix) and store result in DeltaProj
            //
            amp::vmove(c2.getvector(0, m-1), c.getvector(0, m-1));
            for(i=0; i<=m-1; i++)
            {
                if( i==j )
                {
                    deltac(i) = delta;
                }
                else
                {
                    deltac(i) = 0;
                }
            }
            if( k==0 )
            {
                amp::vmove(deltaproj.getvector(0, m-1), deltac.getvector(0, m-1));
            }
            else
            {
                for(i=0; i<=m-1; i++)
                {
                    deltaproj(i) = 0;
                }
                for(i=k; i<=m-1; i++)
                {
                    v = amp::vdotproduct(vt.getrow(i, 0, m-1), deltac.getvector(0, m-1));
                    amp::vadd(deltaproj.getvector(0, m-1), vt.getrow(i, 0, m-1), v);
                }
            }
            
            //
            // now we have DeltaProj such that if C is feasible,
            // then C+DeltaProj is feasible too
            //
            amp::vmove(c2.getvector(0, m-1), c.getvector(0, m-1));
            amp::vadd(c2.getvector(0, m-1), deltaproj.getvector(0, m-1));
            s2 = getglserror<Precision>(n, m, y, w, fmatrix, c2);
            amp::vmove(c2.getvector(0, m-1), c.getvector(0, m-1));
            amp::vsub(c2.getvector(0, m-1), deltaproj.getvector(0, m-1));
            s3 = getglserror<Precision>(n, m, y, w, fmatrix, c2);
            result = result && s2>=s1/(1+threshold) && s3>=s1/(1+threshold);
        }
        return result;
    }


    /*************************************************************************
    Tests whether C is solution of LLS problem
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> getglserror(int n,
        int m,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        const ap::template_1d_array< amp::ampf<Precision> >& c)
    {
        amp::ampf<Precision> result;
        int i;
        amp::ampf<Precision> v;


        result = 0;
        for(i=0; i<=n-1; i++)
        {
            v = amp::vdotproduct(fmatrix.getrow(i, 0, m-1), c.getvector(0, m-1));
            result = result+amp::sqr<Precision>(w(i)*(v-y(i)));
        }
        return result;
    }


    /*************************************************************************
    Subroutine for nonlinear fitting of linear problem
    *************************************************************************/
    template<unsigned int Precision>
    void fitlinearnonlinear(int m,
        bool gradonly,
        const ap::template_2d_array< amp::ampf<Precision> >& xy,
        lsfit::lsfitstate<Precision>& state,
        bool& nlserrors)
    {
        int i;
        int j;
        amp::ampf<Precision> v;


        while( lsfit::lsfitnonlineariteration<Precision>(state) )
        {
            
            //
            // assume that one and only one of flags is set
            // test that we didn't request hessian in hessian-free setting
            //
            if( gradonly && state.needfgh )
            {
                nlserrors = true;
            }
            i = 0;
            if( state.needf )
            {
                i = i+1;
            }
            if( state.needfg )
            {
                i = i+1;
            }
            if( state.needfgh )
            {
                i = i+1;
            }
            if( i!=1 )
            {
                nlserrors = true;
            }
            
            //
            // test that PointIndex is consistent with actual point passed
            //
            for(i=0; i<=m-1; i++)
            {
                nlserrors = nlserrors || xy(state.pointindex,i)!=state.x(i);
            }
            
            //
            // calculate
            //
            if( state.needf )
            {
                v = amp::vdotproduct(state.x.getvector(0, m-1), state.c.getvector(0, m-1));
                state.f = v;
                continue;
            }
            if( state.needfg )
            {
                v = amp::vdotproduct(state.x.getvector(0, m-1), state.c.getvector(0, m-1));
                state.f = v;
                amp::vmove(state.g.getvector(0, m-1), state.x.getvector(0, m-1));
                continue;
            }
            if( state.needfgh )
            {
                v = amp::vdotproduct(state.x.getvector(0, m-1), state.c.getvector(0, m-1));
                state.f = v;
                amp::vmove(state.g.getvector(0, m-1), state.x.getvector(0, m-1));
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        state.h(i,j) = 0;
                    }
                }
                continue;
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool llstestunit_test_silent()
    {
        bool result;


        result = testlls<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool llstestunit_test()
    {
        bool result;


        result = testlls<Precision>(false);
        return result;
    }
} // namespace

#endif
