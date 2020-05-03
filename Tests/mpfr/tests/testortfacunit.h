
#ifndef _testortfacunit_h
#define _testortfacunit_h

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
namespace testortfacunit
{
    template<unsigned int Precision>
    bool testortfac(bool silent);
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixdiff(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int n);
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
    void rmatrixfillsparsea(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> sparcity);
    template<unsigned int Precision>
    void cmatrixfillsparsea(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> sparcity);
    template<unsigned int Precision>
    void internalmatrixmatrixmultiply(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ai1,
        int ai2,
        int aj1,
        int aj2,
        bool transa,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int bi1,
        int bi2,
        int bj1,
        int bj2,
        bool transb,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ci1,
        int ci2,
        int cj1,
        int cj2);
    template<unsigned int Precision>
    void testrqrproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& qrerrors);
    template<unsigned int Precision>
    void testcqrproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& qrerrors);
    template<unsigned int Precision>
    void testrlqproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& lqerrors);
    template<unsigned int Precision>
    void testclqproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& lqerrors);
    template<unsigned int Precision>
    void testrbdproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& bderrors);
    template<unsigned int Precision>
    void testrhessproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision> threshold,
        bool& hesserrors);
    template<unsigned int Precision>
    void testrtdproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision> threshold,
        bool& tderrors);
    template<unsigned int Precision>
    void testctdproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        amp::ampf<Precision> threshold,
        bool& tderrors);
    template<unsigned int Precision>
    bool testortfacunit_test_silent();
    template<unsigned int Precision>
    bool testortfacunit_test();


    /*************************************************************************
    Main unittest subroutine
    *************************************************************************/
    template<unsigned int Precision>
    bool testortfac(bool silent)
    {
        bool result;
        int maxmn;
        amp::ampf<Precision> threshold;
        int passcount;
        int mx;
        ap::template_2d_array< amp::ampf<Precision> > ra;
        ap::template_2d_array< amp::campf<Precision> > ca;
        int m;
        int n;
        int pass;
        int i;
        int j;
        bool rqrerrors;
        bool rlqerrors;
        bool cqrerrors;
        bool clqerrors;
        bool rbderrors;
        bool rhesserrors;
        bool rtderrors;
        bool ctderrors;
        bool waserrors;


        waserrors = false;
        rqrerrors = false;
        rlqerrors = false;
        cqrerrors = false;
        clqerrors = false;
        rbderrors = false;
        rhesserrors = false;
        rtderrors = false;
        ctderrors = false;
        maxmn = 3*ablas::ablasblocksize<Precision>(ra)+1;
        passcount = 1;
        threshold = 5*1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Different problems
        //
        for(mx=1; mx<=maxmn; mx++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // Rectangular factorizations: QR, LQ, bidiagonal
                // Matrix types: zero, dense, sparse
                //
                n = 1+ap::randominteger(mx);
                m = 1+ap::randominteger(mx);
                if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
                {
                    n = mx;
                }
                else
                {
                    m = mx;
                }
                ra.setlength(m, n);
                ca.setlength(m, n);
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        ra(i,j) = 0;
                        ca(i,j) = 0;
                    }
                }
                testrqrproblem<Precision>(ra, m, n, threshold, rqrerrors);
                testrlqproblem<Precision>(ra, m, n, threshold, rlqerrors);
                testcqrproblem<Precision>(ca, m, n, threshold, cqrerrors);
                testclqproblem<Precision>(ca, m, n, threshold, clqerrors);
                testrbdproblem<Precision>(ra, m, n, threshold, rbderrors);
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        ra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        ca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                        ca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                }
                testrqrproblem<Precision>(ra, m, n, threshold, rqrerrors);
                testrlqproblem<Precision>(ra, m, n, threshold, rlqerrors);
                testcqrproblem<Precision>(ca, m, n, threshold, cqrerrors);
                testclqproblem<Precision>(ca, m, n, threshold, clqerrors);
                testrbdproblem<Precision>(ra, m, n, threshold, rbderrors);
                rmatrixfillsparsea<Precision>(ra, m, n, amp::ampf<Precision>("0.95"));
                cmatrixfillsparsea<Precision>(ca, m, n, amp::ampf<Precision>("0.95"));
                testrqrproblem<Precision>(ra, m, n, threshold, rqrerrors);
                testrlqproblem<Precision>(ra, m, n, threshold, rlqerrors);
                testcqrproblem<Precision>(ca, m, n, threshold, cqrerrors);
                testclqproblem<Precision>(ca, m, n, threshold, clqerrors);
                testrbdproblem<Precision>(ra, m, n, threshold, rbderrors);
                
                //
                // Square factorizations: Hessenberg, tridiagonal
                // Matrix types: zero, dense, sparse
                //
                ra.setlength(mx, mx);
                ca.setlength(mx, mx);
                for(i=0; i<=mx-1; i++)
                {
                    for(j=0; j<=mx-1; j++)
                    {
                        ra(i,j) = 0;
                        ca(i,j) = 0;
                    }
                }
                testrhessproblem<Precision>(ra, mx, threshold, rhesserrors);
                for(i=0; i<=mx-1; i++)
                {
                    for(j=0; j<=mx-1; j++)
                    {
                        ra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        ca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                        ca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                }
                testrhessproblem<Precision>(ra, mx, threshold, rhesserrors);
                rmatrixfillsparsea<Precision>(ra, mx, mx, amp::ampf<Precision>("0.95"));
                cmatrixfillsparsea<Precision>(ca, mx, mx, amp::ampf<Precision>("0.95"));
                testrhessproblem<Precision>(ra, mx, threshold, rhesserrors);
                
                //
                // Symetric factorizations: tridiagonal
                // Matrix types: zero, dense, sparse
                //
                ra.setlength(mx, mx);
                ca.setlength(mx, mx);
                for(i=0; i<=mx-1; i++)
                {
                    for(j=0; j<=mx-1; j++)
                    {
                        ra(i,j) = 0;
                        ca(i,j) = 0;
                    }
                }
                testrtdproblem<Precision>(ra, mx, threshold, rtderrors);
                testctdproblem<Precision>(ca, mx, threshold, ctderrors);
                for(i=0; i<=mx-1; i++)
                {
                    for(j=i; j<=mx-1; j++)
                    {
                        ra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        ca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                        ca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                        ra(j,i) = ra(i,j);
                        ca(j,i) = amp::conj<Precision>(ca(i,j));
                    }
                }
                for(i=0; i<=mx-1; i++)
                {
                    ca(i,i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                testrtdproblem<Precision>(ra, mx, threshold, rtderrors);
                testctdproblem<Precision>(ca, mx, threshold, ctderrors);
                rmatrixfillsparsea<Precision>(ra, mx, mx, amp::ampf<Precision>("0.95"));
                cmatrixfillsparsea<Precision>(ca, mx, mx, amp::ampf<Precision>("0.95"));
                for(i=0; i<=mx-1; i++)
                {
                    for(j=i; j<=mx-1; j++)
                    {
                        ra(j,i) = ra(i,j);
                        ca(j,i) = amp::conj<Precision>(ca(i,j));
                    }
                }
                for(i=0; i<=mx-1; i++)
                {
                    ca(i,i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                testrtdproblem<Precision>(ra, mx, threshold, rtderrors);
                testctdproblem<Precision>(ca, mx, threshold, ctderrors);
            }
        }
        
        //
        // report
        //
        waserrors = rqrerrors || rlqerrors || cqrerrors || clqerrors || rbderrors || rhesserrors || rtderrors || ctderrors;
        if( !silent )
        {
            printf("TESTING ORTFAC UNIT\n");
            printf("RQR ERRORS:                              ");
            if( !rqrerrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("RLQ ERRORS:                              ");
            if( !rlqerrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("CQR ERRORS:                              ");
            if( !cqrerrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("CLQ ERRORS:                              ");
            if( !clqerrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("RBD ERRORS:                              ");
            if( !rbderrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("RHESS ERRORS:                            ");
            if( !rhesserrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("RTD ERRORS:                              ");
            if( !rtderrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("CTD ERRORS:                              ");
            if( !ctderrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
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
    Diff
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixdiff(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int j;


        result = 0;
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                result = amp::maximum<Precision>(result, amp::abs<Precision>(b(i,j)-a(i,j)));
            }
        }
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
    Sparse fill
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixfillsparsea(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> sparcity)
    {
        int i;
        int j;


        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( amp::ampf<Precision>::getRandom()>=sparcity )
                {
                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                else
                {
                    a(i,j) = 0;
                }
            }
        }
    }


    /*************************************************************************
    Sparse fill
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixfillsparsea(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> sparcity)
    {
        int i;
        int j;


        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( amp::ampf<Precision>::getRandom()>=sparcity )
                {
                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
                else
                {
                    a(i,j) = 0;
                }
            }
        }
    }


    /*************************************************************************
    Matrix multiplication
    *************************************************************************/
    template<unsigned int Precision>
    void internalmatrixmatrixmultiply(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ai1,
        int ai2,
        int aj1,
        int aj2,
        bool transa,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int bi1,
        int bi2,
        int bj1,
        int bj2,
        bool transb,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ci1,
        int ci2,
        int cj1,
        int cj2)
    {
        int arows;
        int acols;
        int brows;
        int bcols;
        int crows;
        int ccols;
        int i;
        int j;
        int k;
        int l;
        int r;
        amp::ampf<Precision> v;
        ap::template_1d_array< amp::ampf<Precision> > work;
        amp::ampf<Precision> beta;
        amp::ampf<Precision> alpha;


        
        //
        // Pre-setup
        //
        k = ap::maxint(ai2-ai1+1, aj2-aj1+1);
        k = ap::maxint(k, bi2-bi1+1);
        k = ap::maxint(k, bj2-bj1+1);
        work.setbounds(1, k);
        beta = 0;
        alpha = 1;
        
        //
        // Setup
        //
        if( !transa )
        {
            arows = ai2-ai1+1;
            acols = aj2-aj1+1;
        }
        else
        {
            arows = aj2-aj1+1;
            acols = ai2-ai1+1;
        }
        if( !transb )
        {
            brows = bi2-bi1+1;
            bcols = bj2-bj1+1;
        }
        else
        {
            brows = bj2-bj1+1;
            bcols = bi2-bi1+1;
        }
        ap::ap_error::make_assertion(acols==brows);
        if( arows<=0 || acols<=0 || brows<=0 || bcols<=0 )
        {
            return;
        }
        crows = arows;
        ccols = bcols;
        
        //
        // Test WORK
        //
        i = ap::maxint(arows, acols);
        i = ap::maxint(brows, i);
        i = ap::maxint(i, bcols);
        work(1) = 0;
        work(i) = 0;
        
        //
        // Prepare C
        //
        if( beta==0 )
        {
            for(i=ci1; i<=ci2; i++)
            {
                for(j=cj1; j<=cj2; j++)
                {
                    c(i,j) = 0;
                }
            }
        }
        else
        {
            for(i=ci1; i<=ci2; i++)
            {
                amp::vmul(c.getrow(i, cj1, cj2), beta);
            }
        }
        
        //
        // A*B
        //
        if( !transa && !transb )
        {
            for(l=ai1; l<=ai2; l++)
            {
                for(r=bi1; r<=bi2; r++)
                {
                    v = alpha*a(l,aj1+r-bi1);
                    k = ci1+l-ai1;
                    amp::vadd(c.getrow(k, cj1, cj2), b.getrow(r, bj1, bj2), v);
                }
            }
            return;
        }
        
        //
        // A*B'
        //
        if( !transa && transb )
        {
            if( arows*acols<brows*bcols )
            {
                for(r=bi1; r<=bi2; r++)
                {
                    for(l=ai1; l<=ai2; l++)
                    {
                        v = amp::vdotproduct(a.getrow(l, aj1, aj2), b.getrow(r, bj1, bj2));
                        c(ci1+l-ai1,cj1+r-bi1) = c(ci1+l-ai1,cj1+r-bi1)+alpha*v;
                    }
                }
                return;
            }
            else
            {
                for(l=ai1; l<=ai2; l++)
                {
                    for(r=bi1; r<=bi2; r++)
                    {
                        v = amp::vdotproduct(a.getrow(l, aj1, aj2), b.getrow(r, bj1, bj2));
                        c(ci1+l-ai1,cj1+r-bi1) = c(ci1+l-ai1,cj1+r-bi1)+alpha*v;
                    }
                }
                return;
            }
        }
        
        //
        // A'*B
        //
        if( transa && !transb )
        {
            for(l=aj1; l<=aj2; l++)
            {
                for(r=bi1; r<=bi2; r++)
                {
                    v = alpha*a(ai1+r-bi1,l);
                    k = ci1+l-aj1;
                    amp::vadd(c.getrow(k, cj1, cj2), b.getrow(r, bj1, bj2), v);
                }
            }
            return;
        }
        
        //
        // A'*B'
        //
        if( transa && transb )
        {
            if( arows*acols<brows*bcols )
            {
                for(r=bi1; r<=bi2; r++)
                {
                    for(i=1; i<=crows; i++)
                    {
                        work(i) = amp::ampf<Precision>("0.0");
                    }
                    for(l=ai1; l<=ai2; l++)
                    {
                        v = alpha*b(r,bj1+l-ai1);
                        k = cj1+r-bi1;
                        amp::vadd(work.getvector(1, crows), a.getrow(l, aj1, aj2), v);
                    }
                    amp::vadd(c.getcolumn(k, ci1, ci2), work.getvector(1, crows));
                }
                return;
            }
            else
            {
                for(l=aj1; l<=aj2; l++)
                {
                    k = ai2-ai1+1;
                    amp::vmove(work.getvector(1, k), a.getcolumn(l, ai1, ai2));
                    for(r=bi1; r<=bi2; r++)
                    {
                        v = amp::vdotproduct(work.getvector(1, k), b.getrow(r, bj1, bj2));
                        c(ci1+l-aj1,cj1+r-bi1) = c(ci1+l-aj1,cj1+r-bi1)+alpha*v;
                    }
                }
                return;
            }
        }
    }


    /*************************************************************************
    Problem testing
    *************************************************************************/
    template<unsigned int Precision>
    void testrqrproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& qrerrors)
    {
        int i;
        int j;
        int k;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > taub;
        ap::template_2d_array< amp::ampf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > r;
        ap::template_2d_array< amp::ampf<Precision> > q2;
        amp::ampf<Precision> v;


        
        //
        // Test decompose-and-unpack error
        //
        rmatrixmakeacopy<Precision>(a, m, n, b);
        ortfac::rmatrixqr<Precision>(b, m, n, taub);
        ortfac::rmatrixqrunpackq<Precision>(b, m, n, taub, m, q);
        ortfac::rmatrixqrunpackr<Precision>(b, m, n, r);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(q.getrow(i, 0, m-1), r.getcolumn(j, 0, m-1));
                qrerrors = qrerrors || amp::abs<Precision>(v-a(i,j))>threshold;
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=ap::minint(i, n-1)-1; j++)
            {
                qrerrors = qrerrors || r(i,j)!=0;
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(q.getrow(i, 0, m-1), q.getrow(j, 0, m-1));
                if( i==j )
                {
                    v = v-1;
                }
                qrerrors = qrerrors || amp::abs<Precision>(v)>=threshold;
            }
        }
        
        //
        // Test for other errors
        //
        k = 1+ap::randominteger(m);
        ortfac::rmatrixqrunpackq<Precision>(b, m, n, taub, k, q2);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                qrerrors = qrerrors || amp::abs<Precision>(q2(i,j)-q(i,j))>10*amp::ampf<Precision>::getAlgoPascalEpsilon();
            }
        }
    }


    /*************************************************************************
    Problem testing
    *************************************************************************/
    template<unsigned int Precision>
    void testcqrproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& qrerrors)
    {
        int i;
        int j;
        int k;
        ap::template_2d_array< amp::campf<Precision> > b;
        ap::template_1d_array< amp::campf<Precision> > taub;
        ap::template_2d_array< amp::campf<Precision> > q;
        ap::template_2d_array< amp::campf<Precision> > r;
        ap::template_2d_array< amp::campf<Precision> > q2;
        amp::campf<Precision> v;
        int i_;


        
        //
        // Test decompose-and-unpack error
        //
        cmatrixmakeacopy<Precision>(a, m, n, b);
        ortfac::cmatrixqr<Precision>(b, m, n, taub);
        ortfac::cmatrixqrunpackq<Precision>(b, m, n, taub, m, q);
        ortfac::cmatrixqrunpackr<Precision>(b, m, n, r);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=m-1;i_++)
                {
                    v += q(i,i_)*r(i_,j);
                }
                qrerrors = qrerrors || amp::abscomplex<Precision>(v-a(i,j))>threshold;
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=ap::minint(i, n-1)-1; j++)
            {
                qrerrors = qrerrors || r(i,j)!=0;
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=m-1;i_++)
                {
                    v += q(i,i_)*amp::conj(q(j,i_));
                }
                if( i==j )
                {
                    v = v-1;
                }
                qrerrors = qrerrors || amp::abscomplex<Precision>(v)>=threshold;
            }
        }
        
        //
        // Test for other errors
        //
        k = 1+ap::randominteger(m);
        ortfac::cmatrixqrunpackq<Precision>(b, m, n, taub, k, q2);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                qrerrors = qrerrors || amp::abscomplex<Precision>(q2(i,j)-q(i,j))>10*amp::ampf<Precision>::getAlgoPascalEpsilon();
            }
        }
    }


    /*************************************************************************
    Problem testing
    *************************************************************************/
    template<unsigned int Precision>
    void testrlqproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& lqerrors)
    {
        int i;
        int j;
        int k;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > taub;
        ap::template_2d_array< amp::ampf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > l;
        ap::template_2d_array< amp::ampf<Precision> > q2;
        amp::ampf<Precision> v;


        
        //
        // Test decompose-and-unpack error
        //
        rmatrixmakeacopy<Precision>(a, m, n, b);
        ortfac::rmatrixlq<Precision>(b, m, n, taub);
        ortfac::rmatrixlqunpackq<Precision>(b, m, n, taub, n, q);
        ortfac::rmatrixlqunpackl<Precision>(b, m, n, l);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(l.getrow(i, 0, n-1), q.getcolumn(j, 0, n-1));
                lqerrors = lqerrors || amp::abs<Precision>(v-a(i,j))>=threshold;
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=ap::minint(i, n-1)+1; j<=n-1; j++)
            {
                lqerrors = lqerrors || l(i,j)!=0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(q.getrow(i, 0, n-1), q.getrow(j, 0, n-1));
                if( i==j )
                {
                    v = v-1;
                }
                lqerrors = lqerrors || amp::abs<Precision>(v)>=threshold;
            }
        }
        
        //
        // Test for other errors
        //
        k = 1+ap::randominteger(n);
        ortfac::rmatrixlqunpackq<Precision>(b, m, n, taub, k, q2);
        for(i=0; i<=k-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                lqerrors = lqerrors || amp::abs<Precision>(q2(i,j)-q(i,j))>10*amp::ampf<Precision>::getAlgoPascalEpsilon();
            }
        }
    }


    /*************************************************************************
    Problem testing
    *************************************************************************/
    template<unsigned int Precision>
    void testclqproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& lqerrors)
    {
        int i;
        int j;
        int k;
        ap::template_2d_array< amp::campf<Precision> > b;
        ap::template_1d_array< amp::campf<Precision> > taub;
        ap::template_2d_array< amp::campf<Precision> > q;
        ap::template_2d_array< amp::campf<Precision> > l;
        ap::template_2d_array< amp::campf<Precision> > q2;
        amp::campf<Precision> v;
        int i_;


        
        //
        // Test decompose-and-unpack error
        //
        cmatrixmakeacopy<Precision>(a, m, n, b);
        ortfac::cmatrixlq<Precision>(b, m, n, taub);
        ortfac::cmatrixlqunpackq<Precision>(b, m, n, taub, n, q);
        ortfac::cmatrixlqunpackl<Precision>(b, m, n, l);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += l(i,i_)*q(i_,j);
                }
                lqerrors = lqerrors || amp::abscomplex<Precision>(v-a(i,j))>=threshold;
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=ap::minint(i, n-1)+1; j<=n-1; j++)
            {
                lqerrors = lqerrors || l(i,j)!=0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += q(i,i_)*amp::conj(q(j,i_));
                }
                if( i==j )
                {
                    v = v-1;
                }
                lqerrors = lqerrors || amp::abscomplex<Precision>(v)>=threshold;
            }
        }
        
        //
        // Test for other errors
        //
        k = 1+ap::randominteger(n);
        ortfac::cmatrixlqunpackq<Precision>(b, m, n, taub, k, q2);
        for(i=0; i<=k-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                lqerrors = lqerrors || amp::abscomplex<Precision>(q2(i,j)-q(i,j))>10*amp::ampf<Precision>::getAlgoPascalEpsilon();
            }
        }
    }


    /*************************************************************************
    Problem testing
    *************************************************************************/
    template<unsigned int Precision>
    void testrbdproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& bderrors)
    {
        int i;
        int j;
        int k;
        ap::template_2d_array< amp::ampf<Precision> > t;
        ap::template_2d_array< amp::ampf<Precision> > pt;
        ap::template_2d_array< amp::ampf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > r;
        ap::template_2d_array< amp::ampf<Precision> > bd;
        ap::template_2d_array< amp::ampf<Precision> > x;
        ap::template_2d_array< amp::ampf<Precision> > r1;
        ap::template_2d_array< amp::ampf<Precision> > r2;
        ap::template_1d_array< amp::ampf<Precision> > taup;
        ap::template_1d_array< amp::ampf<Precision> > tauq;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > e;
        bool up;
        amp::ampf<Precision> v;
        int mtsize;


        
        //
        // Bidiagonal decomposition error
        //
        rmatrixmakeacopy<Precision>(a, m, n, t);
        ortfac::rmatrixbd<Precision>(t, m, n, tauq, taup);
        ortfac::rmatrixbdunpackq<Precision>(t, m, n, tauq, m, q);
        ortfac::rmatrixbdunpackpt<Precision>(t, m, n, taup, n, pt);
        ortfac::rmatrixbdunpackdiagonals<Precision>(t, m, n, up, d, e);
        bd.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                bd(i,j) = 0;
            }
        }
        for(i=0; i<=ap::minint(m, n)-1; i++)
        {
            bd(i,i) = d(i);
        }
        if( up )
        {
            for(i=0; i<=ap::minint(m, n)-2; i++)
            {
                bd(i,i+1) = e(i);
            }
        }
        else
        {
            for(i=0; i<=ap::minint(m, n)-2; i++)
            {
                bd(i+1,i) = e(i);
            }
        }
        r.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(q.getrow(i, 0, m-1), bd.getcolumn(j, 0, m-1));
                r(i,j) = v;
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(r.getrow(i, 0, n-1), pt.getcolumn(j, 0, n-1));
                bderrors = bderrors || amp::abs<Precision>(v-a(i,j))>threshold;
            }
        }
        
        //
        // Orthogonality test for Q/PT
        //
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(q.getcolumn(i, 0, m-1), q.getcolumn(j, 0, m-1));
                if( i==j )
                {
                    bderrors = bderrors || amp::abs<Precision>(v-1)>threshold;
                }
                else
                {
                    bderrors = bderrors || amp::abs<Precision>(v)>threshold;
                }
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(pt.getrow(i, 0, n-1), pt.getrow(j, 0, n-1));
                if( i==j )
                {
                    bderrors = bderrors || amp::abs<Precision>(v-1)>threshold;
                }
                else
                {
                    bderrors = bderrors || amp::abs<Precision>(v)>threshold;
                }
            }
        }
        
        //
        // Partial unpacking test
        //
        k = 1+ap::randominteger(m);
        ortfac::rmatrixbdunpackq<Precision>(t, m, n, tauq, k, r);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                bderrors = bderrors || amp::abs<Precision>(r(i,j)-q(i,j))>10*amp::ampf<Precision>::getAlgoPascalEpsilon();
            }
        }
        k = 1+ap::randominteger(n);
        ortfac::rmatrixbdunpackpt<Precision>(t, m, n, taup, k, r);
        for(i=0; i<=k-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                bderrors = bderrors || r(i,j)-pt(i,j)!=0;
            }
        }
        
        //
        // Multiplication test
        //
        x.setbounds(0, ap::maxint(m, n)-1, 0, ap::maxint(m, n)-1);
        r.setbounds(0, ap::maxint(m, n)-1, 0, ap::maxint(m, n)-1);
        r1.setbounds(0, ap::maxint(m, n)-1, 0, ap::maxint(m, n)-1);
        r2.setbounds(0, ap::maxint(m, n)-1, 0, ap::maxint(m, n)-1);
        for(i=0; i<=ap::maxint(m, n)-1; i++)
        {
            for(j=0; j<=ap::maxint(m, n)-1; j++)
            {
                x(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
            }
        }
        mtsize = 1+ap::randominteger(ap::maxint(m, n));
        rmatrixmakeacopy<Precision>(x, mtsize, m, r);
        internalmatrixmatrixmultiply<Precision>(r, 0, mtsize-1, 0, m-1, false, q, 0, m-1, 0, m-1, false, r1, 0, mtsize-1, 0, m-1);
        rmatrixmakeacopy<Precision>(x, mtsize, m, r2);
        ortfac::rmatrixbdmultiplybyq<Precision>(t, m, n, tauq, r2, mtsize, m, true, false);
        bderrors = bderrors || rmatrixdiff<Precision>(r1, r2, mtsize, m)>threshold;
        rmatrixmakeacopy<Precision>(x, mtsize, m, r);
        internalmatrixmatrixmultiply<Precision>(r, 0, mtsize-1, 0, m-1, false, q, 0, m-1, 0, m-1, true, r1, 0, mtsize-1, 0, m-1);
        rmatrixmakeacopy<Precision>(x, mtsize, m, r2);
        ortfac::rmatrixbdmultiplybyq<Precision>(t, m, n, tauq, r2, mtsize, m, true, true);
        bderrors = bderrors || rmatrixdiff<Precision>(r1, r2, mtsize, m)>threshold;
        rmatrixmakeacopy<Precision>(x, m, mtsize, r);
        internalmatrixmatrixmultiply<Precision>(q, 0, m-1, 0, m-1, false, r, 0, m-1, 0, mtsize-1, false, r1, 0, m-1, 0, mtsize-1);
        rmatrixmakeacopy<Precision>(x, m, mtsize, r2);
        ortfac::rmatrixbdmultiplybyq<Precision>(t, m, n, tauq, r2, m, mtsize, false, false);
        bderrors = bderrors || rmatrixdiff<Precision>(r1, r2, m, mtsize)>threshold;
        rmatrixmakeacopy<Precision>(x, m, mtsize, r);
        internalmatrixmatrixmultiply<Precision>(q, 0, m-1, 0, m-1, true, r, 0, m-1, 0, mtsize-1, false, r1, 0, m-1, 0, mtsize-1);
        rmatrixmakeacopy<Precision>(x, m, mtsize, r2);
        ortfac::rmatrixbdmultiplybyq<Precision>(t, m, n, tauq, r2, m, mtsize, false, true);
        bderrors = bderrors || rmatrixdiff<Precision>(r1, r2, m, mtsize)>threshold;
        rmatrixmakeacopy<Precision>(x, mtsize, n, r);
        internalmatrixmatrixmultiply<Precision>(r, 0, mtsize-1, 0, n-1, false, pt, 0, n-1, 0, n-1, true, r1, 0, mtsize-1, 0, n-1);
        rmatrixmakeacopy<Precision>(x, mtsize, n, r2);
        ortfac::rmatrixbdmultiplybyp<Precision>(t, m, n, taup, r2, mtsize, n, true, false);
        bderrors = bderrors || rmatrixdiff<Precision>(r1, r2, mtsize, n)>threshold;
        rmatrixmakeacopy<Precision>(x, mtsize, n, r);
        internalmatrixmatrixmultiply<Precision>(r, 0, mtsize-1, 0, n-1, false, pt, 0, n-1, 0, n-1, false, r1, 0, mtsize-1, 0, n-1);
        rmatrixmakeacopy<Precision>(x, mtsize, n, r2);
        ortfac::rmatrixbdmultiplybyp<Precision>(t, m, n, taup, r2, mtsize, n, true, true);
        bderrors = bderrors || rmatrixdiff<Precision>(r1, r2, mtsize, n)>threshold;
        rmatrixmakeacopy<Precision>(x, n, mtsize, r);
        internalmatrixmatrixmultiply<Precision>(pt, 0, n-1, 0, n-1, true, r, 0, n-1, 0, mtsize-1, false, r1, 0, n-1, 0, mtsize-1);
        rmatrixmakeacopy<Precision>(x, n, mtsize, r2);
        ortfac::rmatrixbdmultiplybyp<Precision>(t, m, n, taup, r2, n, mtsize, false, false);
        bderrors = bderrors || rmatrixdiff<Precision>(r1, r2, n, mtsize)>threshold;
        rmatrixmakeacopy<Precision>(x, n, mtsize, r);
        internalmatrixmatrixmultiply<Precision>(pt, 0, n-1, 0, n-1, false, r, 0, n-1, 0, mtsize-1, false, r1, 0, n-1, 0, mtsize-1);
        rmatrixmakeacopy<Precision>(x, n, mtsize, r2);
        ortfac::rmatrixbdmultiplybyp<Precision>(t, m, n, taup, r2, n, mtsize, false, true);
        bderrors = bderrors || rmatrixdiff<Precision>(r1, r2, n, mtsize)>threshold;
    }


    /*************************************************************************
    Problem testing
    *************************************************************************/
    template<unsigned int Precision>
    void testrhessproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision> threshold,
        bool& hesserrors)
    {
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_2d_array< amp::ampf<Precision> > h;
        ap::template_2d_array< amp::ampf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > t1;
        ap::template_2d_array< amp::ampf<Precision> > t2;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        int i;
        int j;
        amp::ampf<Precision> v;


        rmatrixmakeacopy<Precision>(a, n, n, b);
        
        //
        // Decomposition
        //
        ortfac::rmatrixhessenberg<Precision>(b, n, tau);
        ortfac::rmatrixhessenbergunpackq<Precision>(b, n, tau, q);
        ortfac::rmatrixhessenbergunpackh<Precision>(b, n, h);
        
        //
        // Matrix properties
        //
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(q.getcolumn(i, 0, n-1), q.getcolumn(j, 0, n-1));
                if( i==j )
                {
                    v = v-1;
                }
                hesserrors = hesserrors || amp::abs<Precision>(v)>threshold;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=i-2; j++)
            {
                hesserrors = hesserrors || h(i,j)!=0;
            }
        }
        
        //
        // Decomposition error
        //
        t1.setlength(n, n);
        t2.setlength(n, n);
        internalmatrixmatrixmultiply<Precision>(q, 0, n-1, 0, n-1, false, h, 0, n-1, 0, n-1, false, t1, 0, n-1, 0, n-1);
        internalmatrixmatrixmultiply<Precision>(t1, 0, n-1, 0, n-1, false, q, 0, n-1, 0, n-1, true, t2, 0, n-1, 0, n-1);
        hesserrors = hesserrors || rmatrixdiff<Precision>(t2, a, n, n)>threshold;
    }


    /*************************************************************************
    Tridiagonal tester
    *************************************************************************/
    template<unsigned int Precision>
    void testrtdproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision> threshold,
        bool& tderrors)
    {
        int i;
        int j;
        ap::template_2d_array< amp::ampf<Precision> > ua;
        ap::template_2d_array< amp::ampf<Precision> > la;
        ap::template_2d_array< amp::ampf<Precision> > t;
        ap::template_2d_array< amp::ampf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > t2;
        ap::template_2d_array< amp::ampf<Precision> > t3;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > e;
        amp::ampf<Precision> v;


        ua.setbounds(0, n-1, 0, n-1);
        la.setbounds(0, n-1, 0, n-1);
        t.setbounds(0, n-1, 0, n-1);
        q.setbounds(0, n-1, 0, n-1);
        t2.setbounds(0, n-1, 0, n-1);
        t3.setbounds(0, n-1, 0, n-1);
        
        //
        // fill
        //
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                ua(i,j) = 0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=i; j<=n-1; j++)
            {
                ua(i,j) = a(i,j);
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                la(i,j) = 0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=i; j++)
            {
                la(i,j) = a(i,j);
            }
        }
        
        //
        // Test 2tridiagonal: upper
        //
        ortfac::smatrixtd<Precision>(ua, n, true, tau, d, e);
        ortfac::smatrixtdunpackq<Precision>(ua, n, true, tau, q);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                t(i,j) = 0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            t(i,i) = d(i);
        }
        for(i=0; i<=n-2; i++)
        {
            t(i,i+1) = e(i);
            t(i+1,i) = e(i);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(q.getcolumn(i, 0, n-1), a.getcolumn(j, 0, n-1));
                t2(i,j) = v;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(t2.getrow(i, 0, n-1), q.getcolumn(j, 0, n-1));
                t3(i,j) = v;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                tderrors = tderrors || amp::abs<Precision>(t3(i,j)-t(i,j))>threshold;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(q.getrow(i, 0, n-1), q.getrow(j, 0, n-1));
                if( i==j )
                {
                    v = v-1;
                }
                tderrors = tderrors || amp::abs<Precision>(v)>threshold;
            }
        }
        
        //
        // Test 2tridiagonal: lower
        //
        ortfac::smatrixtd<Precision>(la, n, false, tau, d, e);
        ortfac::smatrixtdunpackq<Precision>(la, n, false, tau, q);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                t(i,j) = 0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            t(i,i) = d(i);
        }
        for(i=0; i<=n-2; i++)
        {
            t(i,i+1) = e(i);
            t(i+1,i) = e(i);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(q.getcolumn(i, 0, n-1), a.getcolumn(j, 0, n-1));
                t2(i,j) = v;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(t2.getrow(i, 0, n-1), q.getcolumn(j, 0, n-1));
                t3(i,j) = v;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                tderrors = tderrors || amp::abs<Precision>(t3(i,j)-t(i,j))>threshold;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(q.getrow(i, 0, n-1), q.getrow(j, 0, n-1));
                if( i==j )
                {
                    v = v-1;
                }
                tderrors = tderrors || amp::abs<Precision>(v)>threshold;
            }
        }
    }


    /*************************************************************************
    Hermitian problem tester
    *************************************************************************/
    template<unsigned int Precision>
    void testctdproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        amp::ampf<Precision> threshold,
        bool& tderrors)
    {
        int i;
        int j;
        ap::template_2d_array< amp::campf<Precision> > ua;
        ap::template_2d_array< amp::campf<Precision> > la;
        ap::template_2d_array< amp::campf<Precision> > t;
        ap::template_2d_array< amp::campf<Precision> > q;
        ap::template_2d_array< amp::campf<Precision> > t2;
        ap::template_2d_array< amp::campf<Precision> > t3;
        ap::template_1d_array< amp::campf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > e;
        amp::campf<Precision> v;
        int i_;


        ua.setbounds(0, n-1, 0, n-1);
        la.setbounds(0, n-1, 0, n-1);
        t.setbounds(0, n-1, 0, n-1);
        q.setbounds(0, n-1, 0, n-1);
        t2.setbounds(0, n-1, 0, n-1);
        t3.setbounds(0, n-1, 0, n-1);
        
        //
        // fill
        //
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                ua(i,j) = 0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=i; j<=n-1; j++)
            {
                ua(i,j) = a(i,j);
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                la(i,j) = 0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=i; j++)
            {
                la(i,j) = a(i,j);
            }
        }
        
        //
        // Test 2tridiagonal: upper
        //
        ortfac::hmatrixtd<Precision>(ua, n, true, tau, d, e);
        ortfac::hmatrixtdunpackq<Precision>(ua, n, true, tau, q);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                t(i,j) = 0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            t(i,i) = d(i);
        }
        for(i=0; i<=n-2; i++)
        {
            t(i,i+1) = e(i);
            t(i+1,i) = e(i);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += amp::conj(q(i_,i))*a(i_,j);
                }
                t2(i,j) = v;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += t2(i,i_)*q(i_,j);
                }
                t3(i,j) = v;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                tderrors = tderrors || amp::abscomplex<Precision>(t3(i,j)-t(i,j))>threshold;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += q(i,i_)*amp::conj(q(j,i_));
                }
                if( i==j )
                {
                    v = v-1;
                }
                tderrors = tderrors || amp::abscomplex<Precision>(v)>threshold;
            }
        }
        
        //
        // Test 2tridiagonal: lower
        //
        ortfac::hmatrixtd<Precision>(la, n, false, tau, d, e);
        ortfac::hmatrixtdunpackq<Precision>(la, n, false, tau, q);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                t(i,j) = 0;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            t(i,i) = d(i);
        }
        for(i=0; i<=n-2; i++)
        {
            t(i,i+1) = e(i);
            t(i+1,i) = e(i);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += amp::conj(q(i_,i))*a(i_,j);
                }
                t2(i,j) = v;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += t2(i,i_)*q(i_,j);
                }
                t3(i,j) = v;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                tderrors = tderrors || amp::abscomplex<Precision>(t3(i,j)-t(i,j))>threshold;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += q(i,i_)*amp::conj(q(j,i_));
                }
                if( i==j )
                {
                    v = v-1;
                }
                tderrors = tderrors || amp::abscomplex<Precision>(v)>threshold;
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testortfacunit_test_silent()
    {
        bool result;


        result = testortfac<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testortfacunit_test()
    {
        bool result;


        result = testortfac<Precision>(false);
        return result;
    }
} // namespace

#endif
