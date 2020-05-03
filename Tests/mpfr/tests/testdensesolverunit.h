
#ifndef _testdensesolverunit_h
#define _testdensesolverunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "hqrnd.h"
#include "matgen.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "xblas.h"
#include "densesolver.h"
namespace testdensesolverunit
{
    template<unsigned int Precision>
    bool testdensesolver(bool silent);
    template<unsigned int Precision>
    bool rmatrixchecksolutionm(const ap::template_2d_array< amp::ampf<Precision> >& xe,
        int n,
        int m,
        amp::ampf<Precision> threshold,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_2d_array< amp::ampf<Precision> >& xs);
    template<unsigned int Precision>
    bool rmatrixchecksolution(const ap::template_2d_array< amp::ampf<Precision> >& xe,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_1d_array< amp::ampf<Precision> >& xs);
    template<unsigned int Precision>
    bool rmatrixchecksingularm(int n,
        int m,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_2d_array< amp::ampf<Precision> >& xs);
    template<unsigned int Precision>
    bool rmatrixchecksingular(int n,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_1d_array< amp::ampf<Precision> >& xs);
    template<unsigned int Precision>
    bool cmatrixchecksolutionm(const ap::template_2d_array< amp::campf<Precision> >& xe,
        int n,
        int m,
        amp::ampf<Precision> threshold,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_2d_array< amp::campf<Precision> >& xs);
    template<unsigned int Precision>
    bool cmatrixchecksolution(const ap::template_2d_array< amp::campf<Precision> >& xe,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_1d_array< amp::campf<Precision> >& xs);
    template<unsigned int Precision>
    bool cmatrixchecksingularm(int n,
        int m,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_2d_array< amp::campf<Precision> >& xs);
    template<unsigned int Precision>
    bool cmatrixchecksingular(int n,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_1d_array< amp::campf<Precision> >& xs);
    template<unsigned int Precision>
    void rmatrixmakeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& b);
    template<unsigned int Precision>
    void cmatrixmakeacopy(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& b);
    template<unsigned int Precision>
    void rmatrixdrophalf(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool droplower);
    template<unsigned int Precision>
    void cmatrixdrophalf(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool droplower);
    template<unsigned int Precision>
    void testrsolver(int maxn,
        int maxm,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& rerrors,
        bool& rfserrors);
    template<unsigned int Precision>
    void testspdsolver(int maxn,
        int maxm,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& spderrors,
        bool& rfserrors);
    template<unsigned int Precision>
    void testcsolver(int maxn,
        int maxm,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& cerrors,
        bool& rfserrors);
    template<unsigned int Precision>
    void testhpdsolver(int maxn,
        int maxm,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& hpderrors,
        bool& rfserrors);
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void cunset2d(ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void cunset1d(ap::template_1d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void unsetrep(densesolver::densesolverreport<Precision>& r);
    template<unsigned int Precision>
    void unsetlsrep(densesolver::densesolverlsreport<Precision>& r);
    template<unsigned int Precision>
    bool testdensesolverunit_test_silent();
    template<unsigned int Precision>
    bool testdensesolverunit_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testdensesolver(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > lua;
        ap::template_2d_array< amp::ampf<Precision> > atmp;
        ap::template_1d_array< int > p;
        ap::template_2d_array< amp::ampf<Precision> > xe;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > bv;
        int i;
        int j;
        int k;
        int n;
        int m;
        int pass;
        int taskkind;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> v;
        amp::ampf<Precision> verr;
        int info;
        densesolver::densesolverreport<Precision> rep;
        densesolver::densesolverlsreport<Precision> repls;
        ap::template_2d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > xv;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > tx;
        int maxn;
        int maxm;
        int passcount;
        amp::ampf<Precision> threshold;
        bool rerrors;
        bool cerrors;
        bool spderrors;
        bool hpderrors;
        bool rfserrors;
        bool waserrors;


        maxn = 10;
        maxm = 5;
        passcount = 5;
        threshold = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        rfserrors = false;
        rerrors = false;
        cerrors = false;
        spderrors = false;
        hpderrors = false;
        testrsolver<Precision>(maxn, maxm, passcount, threshold, rerrors, rfserrors);
        testspdsolver<Precision>(maxn, maxm, passcount, threshold, spderrors, rfserrors);
        testcsolver<Precision>(maxn, maxm, passcount, threshold, cerrors, rfserrors);
        testhpdsolver<Precision>(maxn, maxm, passcount, threshold, hpderrors, rfserrors);
        waserrors = rerrors || cerrors || spderrors || hpderrors || rfserrors;
        if( !silent )
        {
            printf("TESTING DENSE SOLVER\n");
            printf("* REAL:                                   ");
            if( rerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* COMPLEX:                                ");
            if( cerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* SPD:                                    ");
            if( spderrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* HPD:                                    ");
            if( hpderrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* ITERATIVE IMPROVEMENT:                  ");
            if( rfserrors )
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
    Checks whether solver results are correct solution.
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixchecksolutionm(const ap::template_2d_array< amp::ampf<Precision> >& xe,
        int n,
        int m,
        amp::ampf<Precision> threshold,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_2d_array< amp::ampf<Precision> >& xs)
    {
        bool result;
        int i;
        int j;


        result = true;
        if( info<=0 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.r1>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.rinf>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    result = result && amp::abs<Precision>(xe(i,j)-xs(i,j))<=threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether solver results are correct solution.
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixchecksolution(const ap::template_2d_array< amp::ampf<Precision> >& xe,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_1d_array< amp::ampf<Precision> >& xs)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > xsm;


        xsm.setlength(n, 1);
        amp::vmove(xsm.getcolumn(0, 0, n-1), xs.getvector(0, n-1));
        result = rmatrixchecksolutionm<Precision>(xe, n, 1, threshold, info, rep, xsm);
        return result;
    }


    /*************************************************************************
    Checks whether solver results indicate singular matrix.
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixchecksingularm(int n,
        int m,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_2d_array< amp::ampf<Precision> >& xs)
    {
        bool result;
        int i;
        int j;


        result = true;
        if( info!=-3 && info!=1 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<0 || rep.r1>1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<0 || rep.rinf>1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            if( info==-3 )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        result = result && xs(i,j)==0;
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether solver results indicate singular matrix.
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixchecksingular(int n,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_1d_array< amp::ampf<Precision> >& xs)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > xsm;


        xsm.setlength(n, 1);
        amp::vmove(xsm.getcolumn(0, 0, n-1), xs.getvector(0, n-1));
        result = rmatrixchecksingularm<Precision>(n, 1, info, rep, xsm);
        return result;
    }


    /*************************************************************************
    Checks whether solver results are correct solution.
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixchecksolutionm(const ap::template_2d_array< amp::campf<Precision> >& xe,
        int n,
        int m,
        amp::ampf<Precision> threshold,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_2d_array< amp::campf<Precision> >& xs)
    {
        bool result;
        int i;
        int j;


        result = true;
        if( info<=0 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.r1>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.rinf>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    result = result && amp::abscomplex<Precision>(xe(i,j)-xs(i,j))<=threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether solver results are correct solution.
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixchecksolution(const ap::template_2d_array< amp::campf<Precision> >& xe,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_1d_array< amp::campf<Precision> >& xs)
    {
        bool result;
        ap::template_2d_array< amp::campf<Precision> > xsm;
        int i_;


        xsm.setlength(n, 1);
        for(i_=0; i_<=n-1;i_++)
        {
            xsm(i_,0) = xs(i_);
        }
        result = cmatrixchecksolutionm<Precision>(xe, n, 1, threshold, info, rep, xsm);
        return result;
    }


    /*************************************************************************
    Checks whether solver results indicate singular matrix.
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixchecksingularm(int n,
        int m,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_2d_array< amp::campf<Precision> >& xs)
    {
        bool result;
        int i;
        int j;


        result = true;
        if( info!=-3 && info!=1 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<0 || rep.r1>1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<0 || rep.rinf>1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            if( info==-3 )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        result = result && xs(i,j)==0;
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether solver results indicate singular matrix.
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixchecksingular(int n,
        int info,
        const densesolver::densesolverreport<Precision>& rep,
        const ap::template_1d_array< amp::campf<Precision> >& xs)
    {
        bool result;
        ap::template_2d_array< amp::campf<Precision> > xsm;
        int i_;


        xsm.setlength(n, 1);
        for(i_=0; i_<=n-1;i_++)
        {
            xsm(i_,0) = xs(i_);
        }
        result = cmatrixchecksingularm<Precision>(n, 1, info, rep, xsm);
        return result;
    }


    /*************************************************************************
    Copy
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixmakeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& b)
    {
        int i;
        int j;


        b.setbounds(0, m-1, 0, n-1);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                b(i,j) = a(i,j);
            }
        }
    }


    /*************************************************************************
    Copy
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixmakeacopy(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& b)
    {
        int i;
        int j;


        b.setbounds(0, m-1, 0, n-1);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                b(i,j) = a(i,j);
            }
        }
    }


    /*************************************************************************
    Drops upper or lower half of the matrix - fills it by special pattern
    which may be used later to ensure that this part wasn't changed
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixdrophalf(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool droplower)
    {
        int i;
        int j;


        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( droplower && i>j || !droplower && i<j )
                {
                    a(i,j) = 1+2*i+3*j;
                }
            }
        }
    }


    /*************************************************************************
    Drops upper or lower half of the matrix - fills it by special pattern
    which may be used later to ensure that this part wasn't changed
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixdrophalf(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool droplower)
    {
        int i;
        int j;


        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( droplower && i>j || !droplower && i<j )
                {
                    a(i,j) = 1+2*i+3*j;
                }
            }
        }
    }


    /*************************************************************************
    Real test
    *************************************************************************/
    template<unsigned int Precision>
    void testrsolver(int maxn,
        int maxm,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& rerrors,
        bool& rfserrors)
    {
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > lua;
        ap::template_2d_array< amp::ampf<Precision> > atmp;
        ap::template_1d_array< int > p;
        ap::template_2d_array< amp::ampf<Precision> > xe;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > bv;
        int i;
        int j;
        int k;
        int n;
        int m;
        int pass;
        int taskkind;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> v;
        amp::ampf<Precision> verr;
        int info;
        densesolver::densesolverreport<Precision> rep;
        densesolver::densesolverlsreport<Precision> repls;
        ap::template_2d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > xv;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > tx;


        
        //
        // General square matrices:
        // * test general solvers
        // * test least squares solver
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                for(m=1; m<=maxm; m++)
                {
                    
                    //
                    // ********************************************************
                    // WELL CONDITIONED TASKS
                    // ability to find correct solution is tested
                    // ********************************************************
                    //
                    // 1. generate random well conditioned matrix A.
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods on original A
                    //
                    matgen::rmatrixrndcond<Precision>(n, amp::ampf<Precision>(1000), a);
                    rmatrixmakeacopy<Precision>(a, n, n, lua);
                    trfac::rmatrixlu<Precision>(lua, n, n, p);
                    xe.setlength(n, m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            xe(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                    }
                    b.setlength(n, m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            v = amp::vdotproduct(a.getrow(i, 0, n-1), xe.getcolumn(j, 0, n-1));
                            b(i,j) = v;
                        }
                    }
                    
                    //
                    // Test solvers
                    //
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset2d<Precision>(x);
                    densesolver::rmatrixsolvem<Precision>(a, n, b, m, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), info, rep, x);
                    rerrors = rerrors || !rmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset1d<Precision>(xv);
                    bv.setlength(n);
                    amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                    densesolver::rmatrixsolve<Precision>(a, n, bv, info, rep, xv);
                    rerrors = rerrors || !rmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset2d<Precision>(x);
                    densesolver::rmatrixlusolvem<Precision>(lua, p, n, b, m, info, rep, x);
                    rerrors = rerrors || !rmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset1d<Precision>(xv);
                    bv.setlength(n);
                    amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                    densesolver::rmatrixlusolve<Precision>(lua, p, n, bv, info, rep, xv);
                    rerrors = rerrors || !rmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset2d<Precision>(x);
                    densesolver::rmatrixmixedsolvem<Precision>(a, lua, p, n, b, m, info, rep, x);
                    rerrors = rerrors || !rmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset1d<Precision>(xv);
                    bv.setlength(n);
                    amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                    densesolver::rmatrixmixedsolve<Precision>(a, lua, p, n, bv, info, rep, xv);
                    rerrors = rerrors || !rmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    
                    //
                    // Test DenseSolverRLS():
                    // * test on original system A*x = b
                    // * test on overdetermined system with the same solution: (A' A')'*x = (b' b')'
                    // * test on underdetermined system with the same solution: (A 0 0 0 ) * z = b
                    //
                    info = 0;
                    unsetlsrep<Precision>(repls);
                    unset1d<Precision>(xv);
                    bv.setlength(n);
                    amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                    densesolver::rmatrixsolvels<Precision>(a, n, n, bv, amp::ampf<Precision>("0.0"), info, repls, xv);
                    if( info<=0 )
                    {
                        rerrors = true;
                    }
                    else
                    {
                        rerrors = rerrors || repls.r2<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || repls.r2>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                        rerrors = rerrors || repls.n!=n || repls.k!=0;
                        for(i=0; i<=n-1; i++)
                        {
                            rerrors = rerrors || amp::abs<Precision>(xe(i,0)-xv(i))>threshold;
                        }
                    }
                    info = 0;
                    unsetlsrep<Precision>(repls);
                    unset1d<Precision>(xv);
                    bv.setlength(2*n);
                    amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                    amp::vmove(bv.getvector(n, 2*n-1), b.getcolumn(0, 0, n-1));
                    atmp.setlength(2*n, n);
                    blas::copymatrix<Precision>(a, 0, n-1, 0, n-1, atmp, 0, n-1, 0, n-1);
                    blas::copymatrix<Precision>(a, 0, n-1, 0, n-1, atmp, n, 2*n-1, 0, n-1);
                    densesolver::rmatrixsolvels<Precision>(atmp, 2*n, n, bv, amp::ampf<Precision>("0.0"), info, repls, xv);
                    if( info<=0 )
                    {
                        rerrors = true;
                    }
                    else
                    {
                        rerrors = rerrors || repls.r2<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || repls.r2>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                        rerrors = rerrors || repls.n!=n || repls.k!=0;
                        for(i=0; i<=n-1; i++)
                        {
                            rerrors = rerrors || amp::abs<Precision>(xe(i,0)-xv(i))>threshold;
                        }
                    }
                    info = 0;
                    unsetlsrep<Precision>(repls);
                    unset1d<Precision>(xv);
                    bv.setlength(n);
                    amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                    atmp.setlength(n, 2*n);
                    blas::copymatrix<Precision>(a, 0, n-1, 0, n-1, atmp, 0, n-1, 0, n-1);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=n; j<=2*n-1; j++)
                        {
                            atmp(i,j) = 0;
                        }
                    }
                    densesolver::rmatrixsolvels<Precision>(atmp, n, 2*n, bv, amp::ampf<Precision>("0.0"), info, repls, xv);
                    if( info<=0 )
                    {
                        rerrors = true;
                    }
                    else
                    {
                        rerrors = rerrors || repls.r2!=0;
                        rerrors = rerrors || repls.n!=2*n || repls.k!=n;
                        for(i=0; i<=n-1; i++)
                        {
                            rerrors = rerrors || amp::abs<Precision>(xe(i,0)-xv(i))>threshold;
                        }
                        for(i=n; i<=2*n-1; i++)
                        {
                            rerrors = rerrors || amp::abs<Precision>(xv(i))>threshold;
                        }
                    }
                    
                    //
                    // ********************************************************
                    // EXACTLY SINGULAR MATRICES
                    // ability to detect singularity is tested
                    // ********************************************************
                    //
                    // 1. generate different types of singular matrices:
                    //    * zero
                    //    * with zero columns
                    //    * with zero rows
                    //    * with equal rows/columns
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods
                    //
                    for(taskkind=0; taskkind<=4; taskkind++)
                    {
                        unset2d<Precision>(a);
                        if( taskkind==0 )
                        {
                            
                            //
                            // all zeros
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j) = 0;
                                }
                            }
                        }
                        if( taskkind==1 )
                        {
                            
                            //
                            // there is zero column
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                            }
                            k = ap::randominteger(n);
                            amp::vmul(a.getcolumn(k, 0, n-1), 0);
                        }
                        if( taskkind==2 )
                        {
                            
                            //
                            // there is zero row
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                            }
                            k = ap::randominteger(n);
                            amp::vmul(a.getrow(k, 0, n-1), 0);
                        }
                        if( taskkind==3 )
                        {
                            
                            //
                            // equal columns
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                            }
                            k = 1+ap::randominteger(n-1);
                            amp::vmove(a.getcolumn(0, 0, n-1), a.getcolumn(k, 0, n-1));
                        }
                        if( taskkind==4 )
                        {
                            
                            //
                            // equal rows
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                            }
                            k = 1+ap::randominteger(n-1);
                            amp::vmove(a.getrow(0, 0, n-1), a.getrow(k, 0, n-1));
                        }
                        xe.setlength(n, m);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=m-1; j++)
                            {
                                xe(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        b.setlength(n, m);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=m-1; j++)
                            {
                                v = amp::vdotproduct(a.getrow(i, 0, n-1), xe.getcolumn(j, 0, n-1));
                                b(i,j) = v;
                            }
                        }
                        rmatrixmakeacopy<Precision>(a, n, n, lua);
                        trfac::rmatrixlu<Precision>(lua, n, n, p);
                        
                        //
                        // Test RMatrixSolveM()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        unset2d<Precision>(x);
                        densesolver::rmatrixsolvem<Precision>(a, n, b, m, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), info, rep, x);
                        rerrors = rerrors || !rmatrixchecksingularm<Precision>(n, m, info, rep, x);
                        
                        //
                        // Test RMatrixSolve()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        unset2d<Precision>(x);
                        bv.setlength(n);
                        amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                        densesolver::rmatrixsolve<Precision>(a, n, bv, info, rep, xv);
                        rerrors = rerrors || !rmatrixchecksingular<Precision>(n, info, rep, xv);
                        
                        //
                        // Test RMatrixLUSolveM()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        unset2d<Precision>(x);
                        densesolver::rmatrixlusolvem<Precision>(lua, p, n, b, m, info, rep, x);
                        rerrors = rerrors || !rmatrixchecksingularm<Precision>(n, m, info, rep, x);
                        
                        //
                        // Test RMatrixLUSolve()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        unset2d<Precision>(x);
                        bv.setlength(n);
                        amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                        densesolver::rmatrixlusolve<Precision>(lua, p, n, bv, info, rep, xv);
                        rerrors = rerrors || !rmatrixchecksingular<Precision>(n, info, rep, xv);
                        
                        //
                        // Test RMatrixMixedSolveM()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        unset2d<Precision>(x);
                        densesolver::rmatrixmixedsolvem<Precision>(a, lua, p, n, b, m, info, rep, x);
                        rerrors = rerrors || !rmatrixchecksingularm<Precision>(n, m, info, rep, x);
                        
                        //
                        // Test RMatrixMixedSolve()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        unset2d<Precision>(x);
                        bv.setlength(n);
                        amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                        densesolver::rmatrixmixedsolve<Precision>(a, lua, p, n, bv, info, rep, xv);
                        rerrors = rerrors || !rmatrixchecksingular<Precision>(n, info, rep, xv);
                    }
                }
            }
        }
        
        //
        // test iterative improvement
        //
        for(pass=1; pass<=passcount; pass++)
        {
            
            //
            // Test iterative improvement matrices
            //
            // A matrix/right part are constructed such that both matrix
            // and solution components are within (-1,+1). Such matrix/right part
            // have nice properties - system can be solved using iterative
            // improvement with |A*x-b| about several ulps of max(1,|b|).
            //
            n = 100;
            a.setlength(n, n);
            b.setlength(n, 1);
            bv.setlength(n);
            tx.setlength(n);
            xv.setlength(n);
            y.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                xv(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                amp::vmove(y.getvector(0, n-1), a.getrow(i, 0, n-1));
                xblas::xdot<Precision>(y, xv, n, tx, v, verr);
                bv(i) = v;
            }
            amp::vmove(b.getcolumn(0, 0, n-1), bv.getvector(0, n-1));
            
            //
            // Test RMatrixSolveM()
            //
            unset2d<Precision>(x);
            densesolver::rmatrixsolvem<Precision>(a, n, b, 1, true, info, rep, x);
            if( info<=0 )
            {
                rfserrors = true;
            }
            else
            {
                xv.setlength(n);
                amp::vmove(xv.getvector(0, n-1), x.getcolumn(0, 0, n-1));
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(y.getvector(0, n-1), a.getrow(i, 0, n-1));
                    xblas::xdot<Precision>(y, xv, n, tx, v, verr);
                    rfserrors = rfserrors || amp::abs<Precision>(v-b(i,0))>8*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::maximum<Precision>(amp::ampf<Precision>(1), amp::abs<Precision>(b(i,0)));
                }
            }
            
            //
            // Test RMatrixSolve()
            //
            unset1d<Precision>(xv);
            densesolver::rmatrixsolve<Precision>(a, n, bv, info, rep, xv);
            if( info<=0 )
            {
                rfserrors = true;
            }
            else
            {
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(y.getvector(0, n-1), a.getrow(i, 0, n-1));
                    xblas::xdot<Precision>(y, xv, n, tx, v, verr);
                    rfserrors = rfserrors || amp::abs<Precision>(v-bv(i))>8*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::maximum<Precision>(amp::ampf<Precision>(1), amp::abs<Precision>(bv(i)));
                }
            }
            
            //
            // Test LS-solver on the same matrix
            //
            densesolver::rmatrixsolvels<Precision>(a, n, n, bv, amp::ampf<Precision>("0.0"), info, repls, xv);
            if( info<=0 )
            {
                rfserrors = true;
            }
            else
            {
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(y.getvector(0, n-1), a.getrow(i, 0, n-1));
                    xblas::xdot<Precision>(y, xv, n, tx, v, verr);
                    rfserrors = rfserrors || amp::abs<Precision>(v-bv(i))>8*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::maximum<Precision>(amp::ampf<Precision>(1), amp::abs<Precision>(bv(i)));
                }
            }
        }
    }


    /*************************************************************************
    SPD test
    *************************************************************************/
    template<unsigned int Precision>
    void testspdsolver(int maxn,
        int maxm,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& spderrors,
        bool& rfserrors)
    {
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > cha;
        ap::template_2d_array< amp::ampf<Precision> > atmp;
        ap::template_1d_array< int > p;
        ap::template_2d_array< amp::ampf<Precision> > xe;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > bv;
        int i;
        int j;
        int k;
        int n;
        int m;
        int pass;
        int taskkind;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> v;
        amp::ampf<Precision> verr;
        bool isupper;
        int info;
        densesolver::densesolverreport<Precision> rep;
        densesolver::densesolverlsreport<Precision> repls;
        ap::template_2d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > xv;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > tx;


        
        //
        // General square matrices:
        // * test general solvers
        // * test least squares solver
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                for(m=1; m<=maxm; m++)
                {
                    
                    //
                    // ********************************************************
                    // WELL CONDITIONED TASKS
                    // ability to find correct solution is tested
                    // ********************************************************
                    //
                    // 1. generate random well conditioned matrix A.
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods on original A
                    //
                    isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
                    matgen::spdmatrixrndcond<Precision>(n, amp::ampf<Precision>(1000), a);
                    rmatrixmakeacopy<Precision>(a, n, n, cha);
                    if( !trfac::spdmatrixcholesky<Precision>(cha, n, isupper) )
                    {
                        spderrors = true;
                        return;
                    }
                    xe.setlength(n, m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            xe(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                    }
                    b.setlength(n, m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            v = amp::vdotproduct(a.getrow(i, 0, n-1), xe.getcolumn(j, 0, n-1));
                            b(i,j) = v;
                        }
                    }
                    rmatrixdrophalf<Precision>(a, n, isupper);
                    rmatrixdrophalf<Precision>(cha, n, isupper);
                    
                    //
                    // Test solvers
                    //
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset2d<Precision>(x);
                    densesolver::spdmatrixsolvem<Precision>(a, n, isupper, b, m, info, rep, x);
                    spderrors = spderrors || !rmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset1d<Precision>(xv);
                    bv.setlength(n);
                    amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                    densesolver::spdmatrixsolve<Precision>(a, n, isupper, bv, info, rep, xv);
                    spderrors = spderrors || !rmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset2d<Precision>(x);
                    densesolver::spdmatrixcholeskysolvem<Precision>(cha, n, isupper, b, m, info, rep, x);
                    spderrors = spderrors || !rmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    unset1d<Precision>(xv);
                    bv.setlength(n);
                    amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                    densesolver::spdmatrixcholeskysolve<Precision>(cha, n, isupper, bv, info, rep, xv);
                    spderrors = spderrors || !rmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    
                    //
                    // ********************************************************
                    // EXACTLY SINGULAR MATRICES
                    // ability to detect singularity is tested
                    // ********************************************************
                    //
                    // 1. generate different types of singular matrices:
                    //    * zero
                    //    * with zero columns
                    //    * with zero rows
                    //    * with equal rows/columns
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods
                    //
                    for(taskkind=0; taskkind<=3; taskkind++)
                    {
                        unset2d<Precision>(a);
                        if( taskkind==0 )
                        {
                            
                            //
                            // all zeros
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j) = 0;
                                }
                            }
                        }
                        if( taskkind==1 )
                        {
                            
                            //
                            // there is zero column
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=i; j<=n-1; j++)
                                {
                                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(j,i) = a(i,j);
                                }
                            }
                            k = ap::randominteger(n);
                            amp::vmul(a.getcolumn(k, 0, n-1), 0);
                            amp::vmul(a.getrow(k, 0, n-1), 0);
                        }
                        if( taskkind==2 )
                        {
                            
                            //
                            // there is zero row
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=i; j<=n-1; j++)
                                {
                                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(j,i) = a(i,j);
                                }
                            }
                            k = ap::randominteger(n);
                            amp::vmul(a.getrow(k, 0, n-1), 0);
                            amp::vmul(a.getcolumn(k, 0, n-1), 0);
                        }
                        if( taskkind==3 )
                        {
                            
                            //
                            // equal columns/rows
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=i; j<=n-1; j++)
                                {
                                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(j,i) = a(i,j);
                                }
                            }
                            k = 1+ap::randominteger(n-1);
                            amp::vmove(a.getcolumn(0, 0, n-1), a.getcolumn(k, 0, n-1));
                            amp::vmove(a.getrow(0, 0, n-1), a.getrow(k, 0, n-1));
                        }
                        xe.setlength(n, m);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=m-1; j++)
                            {
                                xe(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        b.setlength(n, m);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=m-1; j++)
                            {
                                v = amp::vdotproduct(a.getrow(i, 0, n-1), xe.getcolumn(j, 0, n-1));
                                b(i,j) = v;
                            }
                        }
                        rmatrixmakeacopy<Precision>(a, n, n, cha);
                        rmatrixdrophalf<Precision>(a, n, isupper);
                        rmatrixdrophalf<Precision>(cha, n, isupper);
                        
                        //
                        // Test SPDMatrixSolveM()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        unset2d<Precision>(x);
                        densesolver::spdmatrixsolvem<Precision>(a, n, isupper, b, m, info, rep, x);
                        spderrors = spderrors || !rmatrixchecksingularm<Precision>(n, m, info, rep, x);
                        
                        //
                        // Test SPDMatrixSolve()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        unset2d<Precision>(x);
                        bv.setlength(n);
                        amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                        densesolver::spdmatrixsolve<Precision>(a, n, isupper, bv, info, rep, xv);
                        spderrors = spderrors || !rmatrixchecksingular<Precision>(n, info, rep, xv);
                        
                        //
                        // 'equal columns/rows' are degenerate, but
                        // Cholesky matrix with equal columns/rows IS NOT degenerate,
                        // so it is not used for testing purposes.
                        //
                        if( taskkind!=3 )
                        {
                            
                            //
                            // Test SPDMatrixLUSolveM()
                            //
                            info = 0;
                            unsetrep<Precision>(rep);
                            unset2d<Precision>(x);
                            densesolver::spdmatrixcholeskysolvem<Precision>(cha, n, isupper, b, m, info, rep, x);
                            spderrors = spderrors || !rmatrixchecksingularm<Precision>(n, m, info, rep, x);
                            
                            //
                            // Test SPDMatrixLUSolve()
                            //
                            info = 0;
                            unsetrep<Precision>(rep);
                            unset2d<Precision>(x);
                            bv.setlength(n);
                            amp::vmove(bv.getvector(0, n-1), b.getcolumn(0, 0, n-1));
                            densesolver::spdmatrixcholeskysolve<Precision>(cha, n, isupper, bv, info, rep, xv);
                            spderrors = spderrors || !rmatrixchecksingular<Precision>(n, info, rep, xv);
                        }
                    }
                }
            }
        }
    }


    /*************************************************************************
    Real test
    *************************************************************************/
    template<unsigned int Precision>
    void testcsolver(int maxn,
        int maxm,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& cerrors,
        bool& rfserrors)
    {
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > lua;
        ap::template_2d_array< amp::campf<Precision> > atmp;
        ap::template_1d_array< int > p;
        ap::template_2d_array< amp::campf<Precision> > xe;
        ap::template_2d_array< amp::campf<Precision> > b;
        ap::template_1d_array< amp::campf<Precision> > bv;
        int i;
        int j;
        int k;
        int n;
        int m;
        int pass;
        int taskkind;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> verr;
        amp::campf<Precision> v;
        int info;
        densesolver::densesolverreport<Precision> rep;
        densesolver::densesolverlsreport<Precision> repls;
        ap::template_2d_array< amp::campf<Precision> > x;
        ap::template_1d_array< amp::campf<Precision> > xv;
        ap::template_1d_array< amp::campf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > tx;
        int i_;


        
        //
        // General square matrices:
        // * test general solvers
        // * test least squares solver
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                for(m=1; m<=maxm; m++)
                {
                    
                    //
                    // ********************************************************
                    // WELL CONDITIONED TASKS
                    // ability to find correct solution is tested
                    // ********************************************************
                    //
                    // 1. generate random well conditioned matrix A.
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods on original A
                    //
                    matgen::cmatrixrndcond<Precision>(n, amp::ampf<Precision>(1000), a);
                    cmatrixmakeacopy<Precision>(a, n, n, lua);
                    trfac::cmatrixlu<Precision>(lua, n, n, p);
                    xe.setlength(n, m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            xe(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                            xe(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                    }
                    b.setlength(n, m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            v = 0.0;
                            for(i_=0; i_<=n-1;i_++)
                            {
                                v += a(i,i_)*xe(i_,j);
                            }
                            b(i,j) = v;
                        }
                    }
                    
                    //
                    // Test solvers
                    //
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset2d<Precision>(x);
                    densesolver::cmatrixsolvem<Precision>(a, n, b, m, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), info, rep, x);
                    cerrors = cerrors || !cmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset1d<Precision>(xv);
                    bv.setlength(n);
                    for(i_=0; i_<=n-1;i_++)
                    {
                        bv(i_) = b(i_,0);
                    }
                    densesolver::cmatrixsolve<Precision>(a, n, bv, info, rep, xv);
                    cerrors = cerrors || !cmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset2d<Precision>(x);
                    densesolver::cmatrixlusolvem<Precision>(lua, p, n, b, m, info, rep, x);
                    cerrors = cerrors || !cmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset1d<Precision>(xv);
                    bv.setlength(n);
                    for(i_=0; i_<=n-1;i_++)
                    {
                        bv(i_) = b(i_,0);
                    }
                    densesolver::cmatrixlusolve<Precision>(lua, p, n, bv, info, rep, xv);
                    cerrors = cerrors || !cmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset2d<Precision>(x);
                    densesolver::cmatrixmixedsolvem<Precision>(a, lua, p, n, b, m, info, rep, x);
                    cerrors = cerrors || !cmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset1d<Precision>(xv);
                    bv.setlength(n);
                    for(i_=0; i_<=n-1;i_++)
                    {
                        bv(i_) = b(i_,0);
                    }
                    densesolver::cmatrixmixedsolve<Precision>(a, lua, p, n, bv, info, rep, xv);
                    cerrors = cerrors || !cmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    
                    //
                    // ********************************************************
                    // EXACTLY SINGULAR MATRICES
                    // ability to detect singularity is tested
                    // ********************************************************
                    //
                    // 1. generate different types of singular matrices:
                    //    * zero
                    //    * with zero columns
                    //    * with zero rows
                    //    * with equal rows/columns
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods
                    //
                    for(taskkind=0; taskkind<=4; taskkind++)
                    {
                        cunset2d<Precision>(a);
                        if( taskkind==0 )
                        {
                            
                            //
                            // all zeros
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j) = 0;
                                }
                            }
                        }
                        if( taskkind==1 )
                        {
                            
                            //
                            // there is zero column
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                            }
                            k = ap::randominteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(i_,k) = 0*a(i_,k);
                            }
                        }
                        if( taskkind==2 )
                        {
                            
                            //
                            // there is zero row
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                            }
                            k = ap::randominteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(k,i_) = 0*a(k,i_);
                            }
                        }
                        if( taskkind==3 )
                        {
                            
                            //
                            // equal columns
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                            }
                            k = 1+ap::randominteger(n-1);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(i_,0) = a(i_,k);
                            }
                        }
                        if( taskkind==4 )
                        {
                            
                            //
                            // equal rows
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                            }
                            k = 1+ap::randominteger(n-1);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(0,i_) = a(k,i_);
                            }
                        }
                        xe.setlength(n, m);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=m-1; j++)
                            {
                                xe(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        b.setlength(n, m);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=m-1; j++)
                            {
                                v = 0.0;
                                for(i_=0; i_<=n-1;i_++)
                                {
                                    v += a(i,i_)*xe(i_,j);
                                }
                                b(i,j) = v;
                            }
                        }
                        cmatrixmakeacopy<Precision>(a, n, n, lua);
                        trfac::cmatrixlu<Precision>(lua, n, n, p);
                        
                        //
                        // Test CMatrixSolveM()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        cunset2d<Precision>(x);
                        densesolver::cmatrixsolvem<Precision>(a, n, b, m, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), info, rep, x);
                        cerrors = cerrors || !cmatrixchecksingularm<Precision>(n, m, info, rep, x);
                        
                        //
                        // Test CMatrixSolve()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        cunset2d<Precision>(x);
                        bv.setlength(n);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            bv(i_) = b(i_,0);
                        }
                        densesolver::cmatrixsolve<Precision>(a, n, bv, info, rep, xv);
                        cerrors = cerrors || !cmatrixchecksingular<Precision>(n, info, rep, xv);
                        
                        //
                        // Test CMatrixLUSolveM()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        cunset2d<Precision>(x);
                        densesolver::cmatrixlusolvem<Precision>(lua, p, n, b, m, info, rep, x);
                        cerrors = cerrors || !cmatrixchecksingularm<Precision>(n, m, info, rep, x);
                        
                        //
                        // Test CMatrixLUSolve()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        cunset2d<Precision>(x);
                        bv.setlength(n);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            bv(i_) = b(i_,0);
                        }
                        densesolver::cmatrixlusolve<Precision>(lua, p, n, bv, info, rep, xv);
                        cerrors = cerrors || !cmatrixchecksingular<Precision>(n, info, rep, xv);
                        
                        //
                        // Test CMatrixMixedSolveM()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        cunset2d<Precision>(x);
                        densesolver::cmatrixmixedsolvem<Precision>(a, lua, p, n, b, m, info, rep, x);
                        cerrors = cerrors || !cmatrixchecksingularm<Precision>(n, m, info, rep, x);
                        
                        //
                        // Test CMatrixMixedSolve()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        cunset2d<Precision>(x);
                        bv.setlength(n);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            bv(i_) = b(i_,0);
                        }
                        densesolver::cmatrixmixedsolve<Precision>(a, lua, p, n, bv, info, rep, xv);
                        cerrors = cerrors || !cmatrixchecksingular<Precision>(n, info, rep, xv);
                    }
                }
            }
        }
        
        //
        // test iterative improvement
        //
        for(pass=1; pass<=passcount; pass++)
        {
            
            //
            // Test iterative improvement matrices
            //
            // A matrix/right part are constructed such that both matrix
            // and solution components magnitudes are within (-1,+1).
            // Such matrix/right part have nice properties - system can
            // be solved using iterative improvement with |A*x-b| about
            // several ulps of max(1,|b|).
            //
            n = 100;
            a.setlength(n, n);
            b.setlength(n, 1);
            bv.setlength(n);
            tx.setlength(2*n);
            xv.setlength(n);
            y.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                xv(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                xv(i).y = 2*amp::ampf<Precision>::getRandom()-1;
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    y(i_) = a(i,i_);
                }
                xblas::xcdot<Precision>(y, xv, n, tx, v, verr);
                bv(i) = v;
            }
            for(i_=0; i_<=n-1;i_++)
            {
                b(i_,0) = bv(i_);
            }
            
            //
            // Test CMatrixSolveM()
            //
            cunset2d<Precision>(x);
            densesolver::cmatrixsolvem<Precision>(a, n, b, 1, true, info, rep, x);
            if( info<=0 )
            {
                rfserrors = true;
            }
            else
            {
                xv.setlength(n);
                for(i_=0; i_<=n-1;i_++)
                {
                    xv(i_) = x(i_,0);
                }
                for(i=0; i<=n-1; i++)
                {
                    for(i_=0; i_<=n-1;i_++)
                    {
                        y(i_) = a(i,i_);
                    }
                    xblas::xcdot<Precision>(y, xv, n, tx, v, verr);
                    rfserrors = rfserrors || amp::abscomplex<Precision>(v-b(i,0))>8*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::maximum<Precision>(amp::ampf<Precision>(1), amp::abscomplex<Precision>(b(i,0)));
                }
            }
            
            //
            // Test CMatrixSolve()
            //
            cunset1d<Precision>(xv);
            densesolver::cmatrixsolve<Precision>(a, n, bv, info, rep, xv);
            if( info<=0 )
            {
                rfserrors = true;
            }
            else
            {
                for(i=0; i<=n-1; i++)
                {
                    for(i_=0; i_<=n-1;i_++)
                    {
                        y(i_) = a(i,i_);
                    }
                    xblas::xcdot<Precision>(y, xv, n, tx, v, verr);
                    rfserrors = rfserrors || amp::abscomplex<Precision>(v-bv(i))>8*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::maximum<Precision>(amp::ampf<Precision>(1), amp::abscomplex<Precision>(bv(i)));
                }
            }
            
            //
            // TODO: Test LS-solver on the same matrix
            //
        }
    }


    /*************************************************************************
    HPD test
    *************************************************************************/
    template<unsigned int Precision>
    void testhpdsolver(int maxn,
        int maxm,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& hpderrors,
        bool& rfserrors)
    {
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > cha;
        ap::template_2d_array< amp::campf<Precision> > atmp;
        ap::template_1d_array< int > p;
        ap::template_2d_array< amp::campf<Precision> > xe;
        ap::template_2d_array< amp::campf<Precision> > b;
        ap::template_1d_array< amp::campf<Precision> > bv;
        int i;
        int j;
        int k;
        int n;
        int m;
        int pass;
        int taskkind;
        amp::ampf<Precision> mx;
        amp::campf<Precision> v;
        bool isupper;
        int info;
        densesolver::densesolverreport<Precision> rep;
        densesolver::densesolverlsreport<Precision> repls;
        ap::template_2d_array< amp::campf<Precision> > x;
        ap::template_1d_array< amp::campf<Precision> > xv;
        ap::template_1d_array< amp::campf<Precision> > y;
        ap::template_1d_array< amp::campf<Precision> > tx;
        int i_;


        
        //
        // General square matrices:
        // * test general solvers
        // * test least squares solver
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                for(m=1; m<=maxm; m++)
                {
                    
                    //
                    // ********************************************************
                    // WELL CONDITIONED TASKS
                    // ability to find correct solution is tested
                    // ********************************************************
                    //
                    // 1. generate random well conditioned matrix A.
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods on original A
                    //
                    isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
                    matgen::hpdmatrixrndcond<Precision>(n, amp::ampf<Precision>(1000), a);
                    cmatrixmakeacopy<Precision>(a, n, n, cha);
                    if( !trfac::hpdmatrixcholesky<Precision>(cha, n, isupper) )
                    {
                        hpderrors = true;
                        return;
                    }
                    xe.setlength(n, m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            xe(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                            xe(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                        }
                    }
                    b.setlength(n, m);
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=m-1; j++)
                        {
                            v = 0.0;
                            for(i_=0; i_<=n-1;i_++)
                            {
                                v += a(i,i_)*xe(i_,j);
                            }
                            b(i,j) = v;
                        }
                    }
                    cmatrixdrophalf<Precision>(a, n, isupper);
                    cmatrixdrophalf<Precision>(cha, n, isupper);
                    
                    //
                    // Test solvers
                    //
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset2d<Precision>(x);
                    densesolver::hpdmatrixsolvem<Precision>(a, n, isupper, b, m, info, rep, x);
                    hpderrors = hpderrors || !cmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset1d<Precision>(xv);
                    bv.setlength(n);
                    for(i_=0; i_<=n-1;i_++)
                    {
                        bv(i_) = b(i_,0);
                    }
                    densesolver::hpdmatrixsolve<Precision>(a, n, isupper, bv, info, rep, xv);
                    hpderrors = hpderrors || !cmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset2d<Precision>(x);
                    densesolver::hpdmatrixcholeskysolvem<Precision>(cha, n, isupper, b, m, info, rep, x);
                    hpderrors = hpderrors || !cmatrixchecksolutionm<Precision>(xe, n, m, threshold, info, rep, x);
                    info = 0;
                    unsetrep<Precision>(rep);
                    cunset1d<Precision>(xv);
                    bv.setlength(n);
                    for(i_=0; i_<=n-1;i_++)
                    {
                        bv(i_) = b(i_,0);
                    }
                    densesolver::hpdmatrixcholeskysolve<Precision>(cha, n, isupper, bv, info, rep, xv);
                    hpderrors = hpderrors || !cmatrixchecksolution<Precision>(xe, n, threshold, info, rep, xv);
                    
                    //
                    // ********************************************************
                    // EXACTLY SINGULAR MATRICES
                    // ability to detect singularity is tested
                    // ********************************************************
                    //
                    // 1. generate different types of singular matrices:
                    //    * zero
                    //    * with zero columns
                    //    * with zero rows
                    //    * with equal rows/columns
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods
                    //
                    for(taskkind=0; taskkind<=3; taskkind++)
                    {
                        cunset2d<Precision>(a);
                        if( taskkind==0 )
                        {
                            
                            //
                            // all zeros
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a(i,j) = 0;
                                }
                            }
                        }
                        if( taskkind==1 )
                        {
                            
                            //
                            // there is zero column
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=i; j<=n-1; j++)
                                {
                                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                                    if( i==j )
                                    {
                                        a(i,j).y = 0;
                                    }
                                    a(j,i) = a(i,j);
                                }
                            }
                            k = ap::randominteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(i_,k) = 0*a(i_,k);
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(k,i_) = 0*a(k,i_);
                            }
                        }
                        if( taskkind==2 )
                        {
                            
                            //
                            // there is zero row
                            //
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=i; j<=n-1; j++)
                                {
                                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                                    if( i==j )
                                    {
                                        a(i,j).y = 0;
                                    }
                                    a(j,i) = a(i,j);
                                }
                            }
                            k = ap::randominteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(k,i_) = 0*a(k,i_);
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(i_,k) = 0*a(i_,k);
                            }
                        }
                        if( taskkind==3 )
                        {
                            
                            //
                            // equal columns/rows
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a.setlength(n, n);
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=i; j<=n-1; j++)
                                {
                                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                                    if( i==j )
                                    {
                                        a(i,j).y = 0;
                                    }
                                    a(j,i) = a(i,j);
                                }
                            }
                            k = 1+ap::randominteger(n-1);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(i_,0) = a(i_,k);
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a(0,i_) = a(k,i_);
                            }
                        }
                        xe.setlength(n, m);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=m-1; j++)
                            {
                                xe(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        b.setlength(n, m);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=m-1; j++)
                            {
                                v = 0.0;
                                for(i_=0; i_<=n-1;i_++)
                                {
                                    v += a(i,i_)*xe(i_,j);
                                }
                                b(i,j) = v;
                            }
                        }
                        cmatrixmakeacopy<Precision>(a, n, n, cha);
                        cmatrixdrophalf<Precision>(a, n, isupper);
                        cmatrixdrophalf<Precision>(cha, n, isupper);
                        
                        //
                        // Test SPDMatrixSolveM()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        cunset2d<Precision>(x);
                        densesolver::hpdmatrixsolvem<Precision>(a, n, isupper, b, m, info, rep, x);
                        hpderrors = hpderrors || !cmatrixchecksingularm<Precision>(n, m, info, rep, x);
                        
                        //
                        // Test SPDMatrixSolve()
                        //
                        info = 0;
                        unsetrep<Precision>(rep);
                        cunset2d<Precision>(x);
                        bv.setlength(n);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            bv(i_) = b(i_,0);
                        }
                        densesolver::hpdmatrixsolve<Precision>(a, n, isupper, bv, info, rep, xv);
                        hpderrors = hpderrors || !cmatrixchecksingular<Precision>(n, info, rep, xv);
                        
                        //
                        // 'equal columns/rows' are degenerate, but
                        // Cholesky matrix with equal columns/rows IS NOT degenerate,
                        // so it is not used for testing purposes.
                        //
                        if( taskkind!=3 )
                        {
                            
                            //
                            // Test SPDMatrixLUSolveM()
                            //
                            info = 0;
                            unsetrep<Precision>(rep);
                            cunset2d<Precision>(x);
                            densesolver::hpdmatrixcholeskysolvem<Precision>(cha, n, isupper, b, m, info, rep, x);
                            hpderrors = hpderrors || !cmatrixchecksingularm<Precision>(n, m, info, rep, x);
                            
                            //
                            // Test SPDMatrixLUSolve()
                            //
                            info = 0;
                            unsetrep<Precision>(rep);
                            cunset2d<Precision>(x);
                            bv.setlength(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                bv(i_) = b(i_,0);
                            }
                            densesolver::hpdmatrixcholeskysolve<Precision>(cha, n, isupper, bv, info, rep, xv);
                            hpderrors = hpderrors || !cmatrixchecksingular<Precision>(n, info, rep, xv);
                        }
                    }
                }
            }
        }
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
    Unsets real matrix
    *************************************************************************/
    template<unsigned int Precision>
    void cunset2d(ap::template_2d_array< amp::campf<Precision> >& x)
    {
        x.setlength(1, 1);
        x(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets real vector
    *************************************************************************/
    template<unsigned int Precision>
    void cunset1d(ap::template_1d_array< amp::campf<Precision> >& x)
    {
        x.setlength(1);
        x(0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets report
    *************************************************************************/
    template<unsigned int Precision>
    void unsetrep(densesolver::densesolverreport<Precision>& r)
    {
        r.r1 = -1;
        r.rinf = -1;
    }


    /*************************************************************************
    Unsets report
    *************************************************************************/
    template<unsigned int Precision>
    void unsetlsrep(densesolver::densesolverlsreport<Precision>& r)
    {
        r.r2 = -1;
        r.n = -1;
        r.k = -1;
        unset2d<Precision>(r.cx);
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testdensesolverunit_test_silent()
    {
        bool result;


        result = testdensesolver<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testdensesolverunit_test()
    {
        bool result;


        result = testdensesolver<Precision>(false);
        return result;
    }
} // namespace

#endif
