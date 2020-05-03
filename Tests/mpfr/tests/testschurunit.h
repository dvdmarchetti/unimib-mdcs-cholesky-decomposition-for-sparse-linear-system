
#ifndef _testschurunit_h
#define _testschurunit_h

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
#include "hsschur.h"
#include "schur.h"
namespace testschurunit
{
    template<unsigned int Precision>
    bool testschur(bool silent);
    template<unsigned int Precision>
    void fillsparsea(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision> sparcity);
    template<unsigned int Precision>
    void testschurproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision>& materr,
        amp::ampf<Precision>& orterr,
        bool& errstruct,
        bool& wfailed);
    template<unsigned int Precision>
    bool testschurunit_test_silent();
    template<unsigned int Precision>
    bool testschurunit_test();


    /*************************************************************************
    Testing Schur decomposition subroutine
    *************************************************************************/
    template<unsigned int Precision>
    bool testschur(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        int n;
        int maxn;
        int i;
        int j;
        int pass;
        int passcount;
        bool waserrors;
        bool errstruct;
        bool wfailed;
        amp::ampf<Precision> materr;
        amp::ampf<Precision> orterr;
        amp::ampf<Precision> threshold;


        materr = 0;
        orterr = 0;
        errstruct = false;
        wfailed = false;
        waserrors = false;
        maxn = 70;
        passcount = 1;
        threshold = 5*100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        a.setbounds(0, maxn-1, 0, maxn-1);
        
        //
        // zero matrix, several cases
        //
        for(i=0; i<=maxn-1; i++)
        {
            for(j=0; j<=maxn-1; j++)
            {
                a(i,j) = 0;
            }
        }
        for(n=1; n<=maxn; n++)
        {
            if( n>30 && n%2==0 )
            {
                continue;
            }
            testschurproblem<Precision>(a, n, materr, orterr, errstruct, wfailed);
        }
        
        //
        // Dense matrix
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                if( n>30 && n%2==0 )
                {
                    continue;
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                }
                testschurproblem<Precision>(a, n, materr, orterr, errstruct, wfailed);
            }
        }
        
        //
        // Sparse matrices, very sparse matrices, incredible sparse matrices
        //
        for(pass=1; pass<=1; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                if( n>30 && n%3!=0 )
                {
                    continue;
                }
                fillsparsea<Precision>(a, n, amp::ampf<Precision>("0.8"));
                testschurproblem<Precision>(a, n, materr, orterr, errstruct, wfailed);
                fillsparsea<Precision>(a, n, amp::ampf<Precision>("0.9"));
                testschurproblem<Precision>(a, n, materr, orterr, errstruct, wfailed);
                fillsparsea<Precision>(a, n, amp::ampf<Precision>("0.95"));
                testschurproblem<Precision>(a, n, materr, orterr, errstruct, wfailed);
                fillsparsea<Precision>(a, n, amp::ampf<Precision>("0.997"));
                testschurproblem<Precision>(a, n, materr, orterr, errstruct, wfailed);
            }
        }
        
        //
        // report
        //
        waserrors = materr>threshold || orterr>threshold || errstruct || wfailed;
        if( !silent )
        {
            printf("TESTING SCHUR DECOMPOSITION\n");
            printf("Schur decomposition error:               %5.3le\n",
                double(amp::ampf<Precision>(materr).toDouble()));
            printf("Schur orthogonality error:               %5.3le\n",
                double(amp::ampf<Precision>(orterr).toDouble()));
            printf("T matrix structure:                      ");
            if( !errstruct )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("Always converged:                        ");
            if( !wfailed )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("Threshold:                               %5.3le\n",
                double(amp::ampf<Precision>(threshold).toDouble()));
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


    template<unsigned int Precision>
    void fillsparsea(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision> sparcity)
    {
        int i;
        int j;


        for(i=0; i<=n-1; i++)
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


    template<unsigned int Precision>
    void testschurproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision>& materr,
        amp::ampf<Precision>& orterr,
        bool& errstruct,
        bool& wfailed)
    {
        ap::template_2d_array< amp::ampf<Precision> > s;
        ap::template_2d_array< amp::ampf<Precision> > t;
        ap::template_1d_array< amp::ampf<Precision> > sr;
        ap::template_1d_array< amp::ampf<Precision> > astc;
        ap::template_1d_array< amp::ampf<Precision> > sastc;
        int i;
        int j;
        int k;
        amp::ampf<Precision> v;
        amp::ampf<Precision> locerr;


        sr.setbounds(0, n-1);
        astc.setbounds(0, n-1);
        sastc.setbounds(0, n-1);
        
        //
        // Schur decomposition, convergence test
        //
        t.setbounds(0, n-1, 0, n-1);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                t(i,j) = a(i,j);
            }
        }
        if( !schur::rmatrixschur<Precision>(t, n, s) )
        {
            wfailed = true;
            return;
        }
        
        //
        // decomposition error
        //
        locerr = 0;
        for(j=0; j<=n-1; j++)
        {
            amp::vmove(sr.getvector(0, n-1), s.getrow(j, 0, n-1));
            for(k=0; k<=n-1; k++)
            {
                v = amp::vdotproduct(t.getrow(k, 0, n-1), sr.getvector(0, n-1));
                astc(k) = v;
            }
            for(k=0; k<=n-1; k++)
            {
                v = amp::vdotproduct(s.getrow(k, 0, n-1), astc.getvector(0, n-1));
                sastc(k) = v;
            }
            for(k=0; k<=n-1; k++)
            {
                locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(sastc(k)-a(k,j)));
            }
        }
        materr = amp::maximum<Precision>(materr, locerr);
        
        //
        // orthogonality error
        //
        locerr = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(s.getcolumn(i, 0, n-1), s.getcolumn(j, 0, n-1));
                if( i!=j )
                {
                    locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(v));
                }
                else
                {
                    locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(v-1));
                }
            }
        }
        orterr = amp::maximum<Precision>(orterr, locerr);
        
        //
        // T matrix structure
        //
        for(j=0; j<=n-1; j++)
        {
            for(i=j+2; i<=n-1; i++)
            {
                if( t(i,j)!=0 )
                {
                    errstruct = true;
                }
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testschurunit_test_silent()
    {
        bool result;


        result = testschur<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testschurunit_test()
    {
        bool result;


        result = testschur<Precision>(false);
        return result;
    }
} // namespace

#endif
