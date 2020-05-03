
#ifndef _testevdunit_h
#define _testevdunit_h

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
#include "evd.h"
namespace testevdunit
{
    template<unsigned int Precision>
    bool testevd(bool silent);
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
    void rmatrixsymmetricsplit(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& al,
        ap::template_2d_array< amp::ampf<Precision> >& au);
    template<unsigned int Precision>
    void cmatrixhermitiansplit(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& al,
        ap::template_2d_array< amp::campf<Precision> >& au);
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void cunset2d(ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void cunset1d(ap::template_1d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    amp::ampf<Precision> tdtestproduct(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& z,
        const ap::template_1d_array< amp::ampf<Precision> >& lambda);
    template<unsigned int Precision>
    amp::ampf<Precision> testproduct(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& z,
        const ap::template_1d_array< amp::ampf<Precision> >& lambda);
    template<unsigned int Precision>
    amp::ampf<Precision> testort(const ap::template_2d_array< amp::ampf<Precision> >& z,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> testcproduct(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& z,
        const ap::template_1d_array< amp::ampf<Precision> >& lambda);
    template<unsigned int Precision>
    amp::ampf<Precision> testcort(const ap::template_2d_array< amp::campf<Precision> >& z,
        int n);
    template<unsigned int Precision>
    void testsevdproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& al,
        const ap::template_2d_array< amp::ampf<Precision> >& au,
        int n,
        amp::ampf<Precision> threshold,
        bool& serrors,
        int& failc,
        int& runs);
    template<unsigned int Precision>
    void testhevdproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_2d_array< amp::campf<Precision> >& al,
        const ap::template_2d_array< amp::campf<Precision> >& au,
        int n,
        amp::ampf<Precision> threshold,
        bool& herrors,
        int& failc,
        int& runs);
    template<unsigned int Precision>
    void testsevdbiproblem(const ap::template_2d_array< amp::ampf<Precision> >& afull,
        const ap::template_2d_array< amp::ampf<Precision> >& al,
        const ap::template_2d_array< amp::ampf<Precision> >& au,
        int n,
        bool distvals,
        amp::ampf<Precision> threshold,
        bool& serrors,
        int& failc,
        int& runs);
    template<unsigned int Precision>
    void testhevdbiproblem(const ap::template_2d_array< amp::campf<Precision> >& afull,
        const ap::template_2d_array< amp::campf<Precision> >& al,
        const ap::template_2d_array< amp::campf<Precision> >& au,
        int n,
        bool distvals,
        amp::ampf<Precision> threshold,
        bool& herrors,
        int& failc,
        int& runs);
    template<unsigned int Precision>
    void testtdevdproblem(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        amp::ampf<Precision> threshold,
        bool& tderrors,
        int& failc,
        int& runs);
    template<unsigned int Precision>
    void testtdevdbiproblem(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        bool distvals,
        amp::ampf<Precision> threshold,
        bool& serrors,
        int& failc,
        int& runs);
    template<unsigned int Precision>
    void testnsevdproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision> threshold,
        bool& nserrors,
        int& failc,
        int& runs);
    template<unsigned int Precision>
    void testevdset(const int& n,
        const amp::ampf<Precision>& threshold,
        const amp::ampf<Precision>& bithreshold,
        int& failc,
        int& runs,
        bool& nserrors,
        bool& serrors,
        bool& herrors,
        bool& tderrors,
        bool& sbierrors,
        bool& hbierrors,
        bool& tdbierrors);
    template<unsigned int Precision>
    bool testevdunit_test_silent();
    template<unsigned int Precision>
    bool testevdunit_test();


    /*************************************************************************
    Testing symmetric EVD subroutine
    *************************************************************************/
    template<unsigned int Precision>
    bool testevd(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > ra;
        int n;
        int j;
        int failc;
        int runs;
        amp::ampf<Precision> failr;
        amp::ampf<Precision> failthreshold;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> bithreshold;
        bool waserrors;
        bool nserrors;
        bool serrors;
        bool herrors;
        bool tderrors;
        bool sbierrors;
        bool hbierrors;
        bool tdbierrors;
        bool wfailed;


        failthreshold = amp::ampf<Precision>("0.005");
        threshold = 100000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        bithreshold = amp::ampf<Precision>("1.0E-6");
        nserrors = false;
        serrors = false;
        herrors = false;
        tderrors = false;
        sbierrors = false;
        hbierrors = false;
        tdbierrors = false;
        failc = 0;
        runs = 0;
        
        //
        // Test problems
        //
        for(n=1; n<=ablas::ablasblocksize<Precision>(ra); n++)
        {
            testevdset<Precision>(n, threshold, bithreshold, failc, runs, nserrors, serrors, herrors, tderrors, sbierrors, hbierrors, tdbierrors);
        }
        for(j=2; j<=3; j++)
        {
            for(n=j*ablas::ablasblocksize<Precision>(ra)-1; n<=j*ablas::ablasblocksize<Precision>(ra)+1; n++)
            {
                testevdset<Precision>(n, threshold, bithreshold, failc, runs, nserrors, serrors, herrors, tderrors, sbierrors, hbierrors, tdbierrors);
            }
        }
        
        //
        // report
        //
        wfailed = amp::ampf<Precision>(failc)/amp::ampf<Precision>(runs)>failthreshold;
        waserrors = nserrors || serrors || herrors || tderrors || sbierrors || hbierrors || tdbierrors || wfailed;
        if( !silent )
        {
            printf("TESTING EVD UNIT\n");
            printf("NS ERRORS:                               ");
            if( !nserrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("S ERRORS:                                ");
            if( !serrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("H ERRORS:                                ");
            if( !herrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("TD ERRORS:                               ");
            if( !tderrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("SBI ERRORS:                              ");
            if( !sbierrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("HBI ERRORS:                              ");
            if( !hbierrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("TDBI ERRORS:                             ");
            if( !tdbierrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("FAILURE THRESHOLD:                       ");
            if( !wfailed )
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
    Copies A to AL (lower half) and AU (upper half), filling unused parts by
    random garbage.
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixsymmetricsplit(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& al,
        ap::template_2d_array< amp::ampf<Precision> >& au)
    {
        int i;
        int j;


        for(i=0; i<=n-1; i++)
        {
            for(j=i+1; j<=n-1; j++)
            {
                al(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                al(j,i) = a(i,j);
                au(i,j) = a(i,j);
                au(j,i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            al(i,i) = a(i,i);
            au(i,i) = a(i,i);
        }
    }


    /*************************************************************************
    Copies A to AL (lower half) and AU (upper half), filling unused parts by
    random garbage.
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixhermitiansplit(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& al,
        ap::template_2d_array< amp::campf<Precision> >& au)
    {
        int i;
        int j;


        for(i=0; i<=n-1; i++)
        {
            for(j=i+1; j<=n-1; j++)
            {
                al(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                al(j,i) = amp::conj<Precision>(a(i,j));
                au(i,j) = a(i,j);
                au(j,i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            al(i,i) = a(i,i);
            au(i,i) = a(i,i);
        }
    }


    /*************************************************************************
    Unsets 2D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& a)
    {
        a.setbounds(0, 0, 0, 0);
        a(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets 2D array.
    *************************************************************************/
    template<unsigned int Precision>
    void cunset2d(ap::template_2d_array< amp::campf<Precision> >& a)
    {
        a.setbounds(0, 0, 0, 0);
        a(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets 1D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& a)
    {
        a.setbounds(0, 0);
        a(0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets 1D array.
    *************************************************************************/
    template<unsigned int Precision>
    void cunset1d(ap::template_1d_array< amp::campf<Precision> >& a)
    {
        a.setbounds(0, 0);
        a(0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Tests Z*Lambda*Z' against tridiag(D,E).
    Returns relative error.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> tdtestproduct(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& z,
        const ap::template_1d_array< amp::ampf<Precision> >& lambda)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        int k;
        amp::ampf<Precision> v;
        amp::ampf<Precision> mx;


        result = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Calculate V = A[i,j], A = Z*Lambda*Z'
                //
                v = 0;
                for(k=0; k<=n-1; k++)
                {
                    v = v+z(i,k)*lambda(k)*z(j,k);
                }
                
                //
                // Compare
                //
                if( abs(i-j)==0 )
                {
                    result = amp::maximum<Precision>(result, amp::abs<Precision>(v-d(i)));
                }
                if( abs(i-j)==1 )
                {
                    result = amp::maximum<Precision>(result, amp::abs<Precision>(v-e(ap::minint(i, j))));
                }
                if( abs(i-j)>1 )
                {
                    result = amp::maximum<Precision>(result, amp::abs<Precision>(v));
                }
            }
        }
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            mx = amp::maximum<Precision>(mx, amp::abs<Precision>(d(i)));
        }
        for(i=0; i<=n-2; i++)
        {
            mx = amp::maximum<Precision>(mx, amp::abs<Precision>(e(i)));
        }
        if( mx==0 )
        {
            mx = 1;
        }
        result = result/mx;
        return result;
    }


    /*************************************************************************
    Tests Z*Lambda*Z' against A
    Returns relative error.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> testproduct(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& z,
        const ap::template_1d_array< amp::ampf<Precision> >& lambda)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        int k;
        amp::ampf<Precision> v;
        amp::ampf<Precision> mx;


        result = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Calculate V = A[i,j], A = Z*Lambda*Z'
                //
                v = 0;
                for(k=0; k<=n-1; k++)
                {
                    v = v+z(i,k)*lambda(k)*z(j,k);
                }
                
                //
                // Compare
                //
                result = amp::maximum<Precision>(result, amp::abs<Precision>(v-a(i,j)));
            }
        }
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                mx = amp::maximum<Precision>(mx, amp::abs<Precision>(a(i,j)));
            }
        }
        if( mx==0 )
        {
            mx = 1;
        }
        result = result/mx;
        return result;
    }


    /*************************************************************************
    Tests Z*Z' against diag(1...1)
    Returns absolute error.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> testort(const ap::template_2d_array< amp::ampf<Precision> >& z,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;


        result = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(z.getcolumn(i, 0, n-1), z.getcolumn(j, 0, n-1));
                if( i==j )
                {
                    v = v-1;
                }
                result = amp::maximum<Precision>(result, amp::abs<Precision>(v));
            }
        }
        return result;
    }


    /*************************************************************************
    Tests Z*Lambda*Z' against A
    Returns relative error.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> testcproduct(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& z,
        const ap::template_1d_array< amp::ampf<Precision> >& lambda)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        int k;
        amp::campf<Precision> v;
        amp::ampf<Precision> mx;


        result = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Calculate V = A[i,j], A = Z*Lambda*Z'
                //
                v = 0;
                for(k=0; k<=n-1; k++)
                {
                    v = v+z(i,k)*lambda(k)*amp::conj<Precision>(z(j,k));
                }
                
                //
                // Compare
                //
                result = amp::maximum<Precision>(result, amp::abscomplex<Precision>(v-a(i,j)));
            }
        }
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                mx = amp::maximum<Precision>(mx, amp::abscomplex<Precision>(a(i,j)));
            }
        }
        if( mx==0 )
        {
            mx = 1;
        }
        result = result/mx;
        return result;
    }


    /*************************************************************************
    Tests Z*Z' against diag(1...1)
    Returns absolute error.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> testcort(const ap::template_2d_array< amp::campf<Precision> >& z,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::campf<Precision> v;
        int i_;


        result = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += z(i_,i)*amp::conj(z(i_,j));
                }
                if( i==j )
                {
                    v = v-1;
                }
                result = amp::maximum<Precision>(result, amp::abscomplex<Precision>(v));
            }
        }
        return result;
    }


    /*************************************************************************
    Tests SEVD problem
    *************************************************************************/
    template<unsigned int Precision>
    void testsevdproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& al,
        const ap::template_2d_array< amp::ampf<Precision> >& au,
        int n,
        amp::ampf<Precision> threshold,
        bool& serrors,
        int& failc,
        int& runs)
    {
        ap::template_1d_array< amp::ampf<Precision> > lambda;
        ap::template_1d_array< amp::ampf<Precision> > lambdaref;
        ap::template_2d_array< amp::ampf<Precision> > z;
        int i;
        int j;
        amp::ampf<Precision> v;


        
        //
        // Test simple EVD: values and full vectors, lower A
        //
        unset1d<Precision>(lambdaref);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevd<Precision>(al, n, 1, false, lambdaref, z) )
        {
            failc = failc+1;
            return;
        }
        serrors = serrors || testproduct<Precision>(a, n, z, lambdaref)>threshold;
        serrors = serrors || testort<Precision>(z, n)>threshold;
        for(i=0; i<=n-2; i++)
        {
            if( lambdaref(i+1)<lambdaref(i) )
            {
                serrors = true;
                return;
            }
        }
        
        //
        // Test simple EVD: values and full vectors, upper A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevd<Precision>(au, n, 1, true, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        serrors = serrors || testproduct<Precision>(a, n, z, lambda)>threshold;
        serrors = serrors || testort<Precision>(z, n)>threshold;
        for(i=0; i<=n-2; i++)
        {
            if( lambda(i+1)<lambda(i) )
            {
                serrors = true;
                return;
            }
        }
        
        //
        // Test simple EVD: values only, lower A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevd<Precision>(al, n, 0, false, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(i)-lambdaref(i))>threshold;
        }
        
        //
        // Test simple EVD: values only, upper A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevd<Precision>(au, n, 0, true, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(i)-lambdaref(i))>threshold;
        }
    }


    /*************************************************************************
    Tests SEVD problem
    *************************************************************************/
    template<unsigned int Precision>
    void testhevdproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_2d_array< amp::campf<Precision> >& al,
        const ap::template_2d_array< amp::campf<Precision> >& au,
        int n,
        amp::ampf<Precision> threshold,
        bool& herrors,
        int& failc,
        int& runs)
    {
        ap::template_1d_array< amp::ampf<Precision> > lambda;
        ap::template_1d_array< amp::ampf<Precision> > lambdaref;
        ap::template_2d_array< amp::campf<Precision> > z;
        int i;
        int j;
        amp::campf<Precision> v;


        
        //
        // Test simple EVD: values and full vectors, lower A
        //
        unset1d<Precision>(lambdaref);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevd<Precision>(al, n, 1, false, lambdaref, z) )
        {
            failc = failc+1;
            return;
        }
        herrors = herrors || testcproduct<Precision>(a, n, z, lambdaref)>threshold;
        herrors = herrors || testcort<Precision>(z, n)>threshold;
        for(i=0; i<=n-2; i++)
        {
            if( lambdaref(i+1)<lambdaref(i) )
            {
                herrors = true;
                return;
            }
        }
        
        //
        // Test simple EVD: values and full vectors, upper A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevd<Precision>(au, n, 1, true, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        herrors = herrors || testcproduct<Precision>(a, n, z, lambda)>threshold;
        herrors = herrors || testcort<Precision>(z, n)>threshold;
        for(i=0; i<=n-2; i++)
        {
            if( lambda(i+1)<lambda(i) )
            {
                herrors = true;
                return;
            }
        }
        
        //
        // Test simple EVD: values only, lower A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevd<Precision>(al, n, 0, false, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(i)-lambdaref(i))>threshold;
        }
        
        //
        // Test simple EVD: values only, upper A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevd<Precision>(au, n, 0, true, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(i)-lambdaref(i))>threshold;
        }
    }


    /*************************************************************************
    Tests EVD problem

    DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                    are solving sparse task with  lots  of  zero  eigenvalues.
                    In such cases some tests related to the  eigenvectors  are
                    not performed.
    *************************************************************************/
    template<unsigned int Precision>
    void testsevdbiproblem(const ap::template_2d_array< amp::ampf<Precision> >& afull,
        const ap::template_2d_array< amp::ampf<Precision> >& al,
        const ap::template_2d_array< amp::ampf<Precision> >& au,
        int n,
        bool distvals,
        amp::ampf<Precision> threshold,
        bool& serrors,
        int& failc,
        int& runs)
    {
        ap::template_1d_array< amp::ampf<Precision> > lambda;
        ap::template_1d_array< amp::ampf<Precision> > lambdaref;
        ap::template_2d_array< amp::ampf<Precision> > z;
        ap::template_2d_array< amp::ampf<Precision> > zref;
        ap::template_2d_array< amp::ampf<Precision> > a1;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        ap::template_2d_array< amp::ampf<Precision> > ar;
        bool wsucc;
        int i;
        int j;
        int k;
        int m;
        int i1;
        int i2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;


        lambdaref.setbounds(0, n-1);
        zref.setbounds(0, n-1, 0, n-1);
        a1.setbounds(0, n-1, 0, n-1);
        a2.setbounds(0, n-1, 0, n-1);
        
        //
        // Reference EVD
        //
        runs = runs+1;
        if( !evd::smatrixevd<Precision>(afull, n, 1, true, lambdaref, zref) )
        {
            failc = failc+1;
            return;
        }
        
        //
        // Select random interval boundaries.
        // If there are non-distinct eigenvalues at the boundaries,
        // we move indexes further until values splits. It is done to
        // avoid situations where we can't get definite answer.
        //
        i1 = ap::randominteger(n);
        i2 = i1+ap::randominteger(n-i1);
        while( i1>0 )
        {
            if( amp::abs<Precision>(lambdaref(i1-1)-lambdaref(i1))>10*threshold )
            {
                break;
            }
            i1 = i1-1;
        }
        while( i2<n-1 )
        {
            if( amp::abs<Precision>(lambdaref(i2+1)-lambdaref(i2))>10*threshold )
            {
                break;
            }
            i2 = i2+1;
        }
        
        //
        // Select A, B
        //
        if( i1>0 )
        {
            a = amp::ampf<Precision>("0.5")*(lambdaref(i1)+lambdaref(i1-1));
        }
        else
        {
            a = lambdaref(0)-1;
        }
        if( i2<n-1 )
        {
            b = amp::ampf<Precision>("0.5")*(lambdaref(i2)+lambdaref(i2+1));
        }
        else
        {
            b = lambdaref(n-1)+1;
        }
        
        //
        // Test interval, no vectors, lower A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevdr<Precision>(al, n, 0, false, a, b, m, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test interval, no vectors, upper A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevdr<Precision>(au, n, 0, true, a, b, m, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test indexes, no vectors, lower A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevdi<Precision>(al, n, 0, false, i1, i2, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test indexes, no vectors, upper A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevdi<Precision>(au, n, 0, true, i1, i2, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test interval, vectors, lower A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevdr<Precision>(al, n, 1, false, a, b, m, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            
            //
            // Distinct eigenvalues, test vectors
            //
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(z.getcolumn(j, 0, n-1), zref.getcolumn(i1+j, 0, n-1));
                if( v<0 )
                {
                    amp::vmul(z.getcolumn(j, 0, n-1), -1);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    serrors = serrors || amp::abs<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
        
        //
        // Test interval, vectors, upper A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevdr<Precision>(au, n, 1, true, a, b, m, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            
            //
            // Distinct eigenvalues, test vectors
            //
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(z.getcolumn(j, 0, n-1), zref.getcolumn(i1+j, 0, n-1));
                if( v<0 )
                {
                    amp::vmul(z.getcolumn(j, 0, n-1), -1);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    serrors = serrors || amp::abs<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
        
        //
        // Test indexes, vectors, lower A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevdi<Precision>(al, n, 1, false, i1, i2, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            
            //
            // Distinct eigenvalues, test vectors
            //
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(z.getcolumn(j, 0, n-1), zref.getcolumn(i1+j, 0, n-1));
                if( v<0 )
                {
                    amp::vmul(z.getcolumn(j, 0, n-1), -1);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    serrors = serrors || amp::abs<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
        
        //
        // Test indexes, vectors, upper A
        //
        unset1d<Precision>(lambda);
        unset2d<Precision>(z);
        runs = runs+1;
        if( !evd::smatrixevdi<Precision>(au, n, 1, true, i1, i2, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            
            //
            // Distinct eigenvalues, test vectors
            //
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(z.getcolumn(j, 0, n-1), zref.getcolumn(i1+j, 0, n-1));
                if( v<0 )
                {
                    amp::vmul(z.getcolumn(j, 0, n-1), -1);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    serrors = serrors || amp::abs<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
    }


    /*************************************************************************
    Tests EVD problem

    DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                    are solving sparse task with  lots  of  zero  eigenvalues.
                    In such cases some tests related to the  eigenvectors  are
                    not performed.
    *************************************************************************/
    template<unsigned int Precision>
    void testhevdbiproblem(const ap::template_2d_array< amp::campf<Precision> >& afull,
        const ap::template_2d_array< amp::campf<Precision> >& al,
        const ap::template_2d_array< amp::campf<Precision> >& au,
        int n,
        bool distvals,
        amp::ampf<Precision> threshold,
        bool& herrors,
        int& failc,
        int& runs)
    {
        ap::template_1d_array< amp::ampf<Precision> > lambda;
        ap::template_1d_array< amp::ampf<Precision> > lambdaref;
        ap::template_2d_array< amp::campf<Precision> > z;
        ap::template_2d_array< amp::campf<Precision> > zref;
        ap::template_2d_array< amp::campf<Precision> > a1;
        ap::template_2d_array< amp::campf<Precision> > a2;
        ap::template_2d_array< amp::campf<Precision> > ar;
        bool wsucc;
        int i;
        int j;
        int k;
        int m;
        int i1;
        int i2;
        amp::campf<Precision> v;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        int i_;


        lambdaref.setbounds(0, n-1);
        zref.setbounds(0, n-1, 0, n-1);
        a1.setbounds(0, n-1, 0, n-1);
        a2.setbounds(0, n-1, 0, n-1);
        
        //
        // Reference EVD
        //
        runs = runs+1;
        if( !evd::hmatrixevd<Precision>(afull, n, 1, true, lambdaref, zref) )
        {
            failc = failc+1;
            return;
        }
        
        //
        // Select random interval boundaries.
        // If there are non-distinct eigenvalues at the boundaries,
        // we move indexes further until values splits. It is done to
        // avoid situations where we can't get definite answer.
        //
        i1 = ap::randominteger(n);
        i2 = i1+ap::randominteger(n-i1);
        while( i1>0 )
        {
            if( amp::abs<Precision>(lambdaref(i1-1)-lambdaref(i1))>10*threshold )
            {
                break;
            }
            i1 = i1-1;
        }
        while( i2<n-1 )
        {
            if( amp::abs<Precision>(lambdaref(i2+1)-lambdaref(i2))>10*threshold )
            {
                break;
            }
            i2 = i2+1;
        }
        
        //
        // Select A, B
        //
        if( i1>0 )
        {
            a = amp::ampf<Precision>("0.5")*(lambdaref(i1)+lambdaref(i1-1));
        }
        else
        {
            a = lambdaref(0)-1;
        }
        if( i2<n-1 )
        {
            b = amp::ampf<Precision>("0.5")*(lambdaref(i2)+lambdaref(i2+1));
        }
        else
        {
            b = lambdaref(n-1)+1;
        }
        
        //
        // Test interval, no vectors, lower A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevdr<Precision>(al, n, 0, false, a, b, m, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test interval, no vectors, upper A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevdr<Precision>(au, n, 0, true, a, b, m, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test indexes, no vectors, lower A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevdi<Precision>(al, n, 0, false, i1, i2, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test indexes, no vectors, upper A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevdi<Precision>(au, n, 0, true, i1, i2, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test interval, vectors, lower A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevdr<Precision>(al, n, 1, false, a, b, m, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            
            //
            // Distinct eigenvalues, test vectors
            //
            for(j=0; j<=m-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += z(i_,j)*amp::conj(zref(i_,i1+j));
                }
                v = amp::conj<Precision>(v/amp::abscomplex<Precision>(v));
                for(i_=0; i_<=n-1;i_++)
                {
                    z(i_,j) = v*z(i_,j);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    herrors = herrors || amp::abscomplex<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
        
        //
        // Test interval, vectors, upper A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevdr<Precision>(au, n, 1, true, a, b, m, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            
            //
            // Distinct eigenvalues, test vectors
            //
            for(j=0; j<=m-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += z(i_,j)*amp::conj(zref(i_,i1+j));
                }
                v = amp::conj<Precision>(v/amp::abscomplex<Precision>(v));
                for(i_=0; i_<=n-1;i_++)
                {
                    z(i_,j) = v*z(i_,j);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    herrors = herrors || amp::abscomplex<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
        
        //
        // Test indexes, vectors, lower A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevdi<Precision>(al, n, 1, false, i1, i2, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            
            //
            // Distinct eigenvalues, test vectors
            //
            for(j=0; j<=m-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += z(i_,j)*amp::conj(zref(i_,i1+j));
                }
                v = amp::conj<Precision>(v/amp::abscomplex<Precision>(v));
                for(i_=0; i_<=n-1;i_++)
                {
                    z(i_,j) = v*z(i_,j);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    herrors = herrors || amp::abscomplex<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
        
        //
        // Test indexes, vectors, upper A
        //
        unset1d<Precision>(lambda);
        cunset2d<Precision>(z);
        runs = runs+1;
        if( !evd::hmatrixevdi<Precision>(au, n, 1, true, i1, i2, lambda, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            herrors = herrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            
            //
            // Distinct eigenvalues, test vectors
            //
            for(j=0; j<=m-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    v += z(i_,j)*amp::conj(zref(i_,i1+j));
                }
                v = amp::conj<Precision>(v/amp::abscomplex<Precision>(v));
                for(i_=0; i_<=n-1;i_++)
                {
                    z(i_,j) = v*z(i_,j);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    herrors = herrors || amp::abscomplex<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
    }


    /*************************************************************************
    Tests EVD problem
    *************************************************************************/
    template<unsigned int Precision>
    void testtdevdproblem(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        amp::ampf<Precision> threshold,
        bool& tderrors,
        int& failc,
        int& runs)
    {
        ap::template_1d_array< amp::ampf<Precision> > lambda;
        ap::template_1d_array< amp::ampf<Precision> > ee;
        ap::template_1d_array< amp::ampf<Precision> > lambda2;
        ap::template_2d_array< amp::ampf<Precision> > z;
        ap::template_2d_array< amp::ampf<Precision> > zref;
        ap::template_2d_array< amp::ampf<Precision> > a1;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        bool wsucc;
        int i;
        int j;
        amp::ampf<Precision> v;


        lambda.setbounds(0, n-1);
        lambda2.setbounds(0, n-1);
        zref.setbounds(0, n-1, 0, n-1);
        a1.setbounds(0, n-1, 0, n-1);
        a2.setbounds(0, n-1, 0, n-1);
        if( n>1 )
        {
            ee.setbounds(0, n-2);
        }
        
        //
        // Test simple EVD: values and full vectors
        //
        for(i=0; i<=n-1; i++)
        {
            lambda(i) = d(i);
        }
        for(i=0; i<=n-2; i++)
        {
            ee(i) = e(i);
        }
        runs = runs+1;
        wsucc = evd::smatrixtdevd<Precision>(lambda, ee, n, 2, z);
        if( !wsucc )
        {
            failc = failc+1;
            return;
        }
        tderrors = tderrors || tdtestproduct<Precision>(d, e, n, z, lambda)>threshold;
        tderrors = tderrors || testort<Precision>(z, n)>threshold;
        for(i=0; i<=n-2; i++)
        {
            if( lambda(i+1)<lambda(i) )
            {
                tderrors = true;
                return;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                zref(i,j) = z(i,j);
            }
        }
        
        //
        // Test values only variant
        //
        for(i=0; i<=n-1; i++)
        {
            lambda2(i) = d(i);
        }
        for(i=0; i<=n-2; i++)
        {
            ee(i) = e(i);
        }
        runs = runs+1;
        wsucc = evd::smatrixtdevd<Precision>(lambda2, ee, n, 0, z);
        if( !wsucc )
        {
            failc = failc+1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            tderrors = tderrors || amp::abs<Precision>(lambda2(i)-lambda(i))>threshold;
        }
        
        //
        // Test multiplication variant
        //
        for(i=0; i<=n-1; i++)
        {
            lambda2(i) = d(i);
        }
        for(i=0; i<=n-2; i++)
        {
            ee(i) = e(i);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a1(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                a2(i,j) = a1(i,j);
            }
        }
        runs = runs+1;
        wsucc = evd::smatrixtdevd<Precision>(lambda2, ee, n, 1, a1);
        if( !wsucc )
        {
            failc = failc+1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            tderrors = tderrors || amp::abs<Precision>(lambda2(i)-lambda(i))>threshold;
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(a2.getrow(i, 0, n-1), zref.getcolumn(j, 0, n-1));
                
                //
                // next line is a bit complicated because
                // depending on algorithm used we can get either
                // z or -z as eigenvector. so we compare result
                // with both A*ZRef and -A*ZRef
                //
                tderrors = tderrors || amp::abs<Precision>(v-a1(i,j))>threshold && amp::abs<Precision>(v+a1(i,j))>threshold;
            }
        }
        
        //
        // Test first row variant
        //
        for(i=0; i<=n-1; i++)
        {
            lambda2(i) = d(i);
        }
        for(i=0; i<=n-2; i++)
        {
            ee(i) = e(i);
        }
        runs = runs+1;
        wsucc = evd::smatrixtdevd<Precision>(lambda2, ee, n, 3, z);
        if( !wsucc )
        {
            failc = failc+1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            tderrors = tderrors || amp::abs<Precision>(lambda2(i)-lambda(i))>threshold;
            
            //
            // next line is a bit complicated because
            // depending on algorithm used we can get either
            // z or -z as eigenvector. so we compare result
            // with both z and -z
            //
            tderrors = tderrors || amp::abs<Precision>(z(0,i)-zref(0,i))>threshold && amp::abs<Precision>(z(0,i)+zref(0,i))>threshold;
        }
    }


    /*************************************************************************
    Tests EVD problem

    DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                    are solving sparse task with  lots  of  zero  eigenvalues.
                    In such cases some tests related to the  eigenvectors  are
                    not performed.
    *************************************************************************/
    template<unsigned int Precision>
    void testtdevdbiproblem(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        bool distvals,
        amp::ampf<Precision> threshold,
        bool& serrors,
        int& failc,
        int& runs)
    {
        ap::template_1d_array< amp::ampf<Precision> > lambda;
        ap::template_1d_array< amp::ampf<Precision> > lambdaref;
        ap::template_2d_array< amp::ampf<Precision> > z;
        ap::template_2d_array< amp::ampf<Precision> > zref;
        ap::template_2d_array< amp::ampf<Precision> > a1;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        ap::template_2d_array< amp::ampf<Precision> > ar;
        bool wsucc;
        int i;
        int j;
        int k;
        int m;
        int i1;
        int i2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;


        lambdaref.setbounds(0, n-1);
        zref.setbounds(0, n-1, 0, n-1);
        a1.setbounds(0, n-1, 0, n-1);
        a2.setbounds(0, n-1, 0, n-1);
        
        //
        // Reference EVD
        //
        lambdaref.setlength(n);
        amp::vmove(lambdaref.getvector(0, n-1), d.getvector(0, n-1));
        runs = runs+1;
        if( !evd::smatrixtdevd<Precision>(lambdaref, e, n, 2, zref) )
        {
            failc = failc+1;
            return;
        }
        
        //
        // Select random interval boundaries.
        // If there are non-distinct eigenvalues at the boundaries,
        // we move indexes further until values splits. It is done to
        // avoid situations where we can't get definite answer.
        //
        i1 = ap::randominteger(n);
        i2 = i1+ap::randominteger(n-i1);
        while( i1>0 )
        {
            if( amp::abs<Precision>(lambdaref(i1-1)-lambdaref(i1))>10*threshold )
            {
                break;
            }
            i1 = i1-1;
        }
        while( i2<n-1 )
        {
            if( amp::abs<Precision>(lambdaref(i2+1)-lambdaref(i2))>10*threshold )
            {
                break;
            }
            i2 = i2+1;
        }
        
        //
        // Test different combinations
        //
        
        //
        // Select A, B
        //
        if( i1>0 )
        {
            a = amp::ampf<Precision>("0.5")*(lambdaref(i1)+lambdaref(i1-1));
        }
        else
        {
            a = lambdaref(0)-1;
        }
        if( i2<n-1 )
        {
            b = amp::ampf<Precision>("0.5")*(lambdaref(i2)+lambdaref(i2+1));
        }
        else
        {
            b = lambdaref(n-1)+1;
        }
        
        //
        // Test interval, no vectors
        //
        lambda.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            lambda(i) = d(i);
        }
        runs = runs+1;
        if( !evd::smatrixtdevdr<Precision>(lambda, e, n, 0, a, b, m, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test indexes, no vectors
        //
        lambda.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            lambda(i) = d(i);
        }
        runs = runs+1;
        if( !evd::smatrixtdevdi<Precision>(lambda, e, n, 0, i1, i2, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        
        //
        // Test interval, transform vectors
        //
        lambda.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            lambda(i) = d(i);
        }
        a1.setbounds(0, n-1, 0, n-1);
        a2.setbounds(0, n-1, 0, n-1);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a1(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                a2(i,j) = a1(i,j);
            }
        }
        runs = runs+1;
        if( !evd::smatrixtdevdr<Precision>(lambda, e, n, 1, a, b, m, a1) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            ar.setbounds(0, n-1, 0, m-1);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    v = amp::vdotproduct(a2.getrow(i, 0, n-1), zref.getcolumn(i1+j, 0, n-1));
                    ar(i,j) = v;
                }
            }
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(a1.getcolumn(j, 0, n-1), ar.getcolumn(j, 0, n-1));
                if( v<0 )
                {
                    amp::vmul(ar.getcolumn(j, 0, n-1), -1);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    serrors = serrors || amp::abs<Precision>(a1(i,j)-ar(i,j))>threshold;
                }
            }
        }
        
        //
        // Test indexes, transform vectors
        //
        lambda.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            lambda(i) = d(i);
        }
        a1.setbounds(0, n-1, 0, n-1);
        a2.setbounds(0, n-1, 0, n-1);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a1(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                a2(i,j) = a1(i,j);
            }
        }
        runs = runs+1;
        if( !evd::smatrixtdevdi<Precision>(lambda, e, n, 1, i1, i2, a1) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            ar.setbounds(0, n-1, 0, m-1);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    v = amp::vdotproduct(a2.getrow(i, 0, n-1), zref.getcolumn(i1+j, 0, n-1));
                    ar(i,j) = v;
                }
            }
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(a1.getcolumn(j, 0, n-1), ar.getcolumn(j, 0, n-1));
                if( v<0 )
                {
                    amp::vmul(ar.getcolumn(j, 0, n-1), -1);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    serrors = serrors || amp::abs<Precision>(a1(i,j)-ar(i,j))>threshold;
                }
            }
        }
        
        //
        // Test interval, do not transform vectors
        //
        lambda.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            lambda(i) = d(i);
        }
        z.setbounds(0, 0, 0, 0);
        runs = runs+1;
        if( !evd::smatrixtdevdr<Precision>(lambda, e, n, 2, a, b, m, z) )
        {
            failc = failc+1;
            return;
        }
        if( m!=i2-i1+1 )
        {
            failc = failc+1;
            return;
        }
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(z.getcolumn(j, 0, n-1), zref.getcolumn(i1+j, 0, n-1));
                if( v<0 )
                {
                    amp::vmul(z.getcolumn(j, 0, n-1), -1);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    serrors = serrors || amp::abs<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
        
        //
        // Test indexes, do not transform vectors
        //
        lambda.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            lambda(i) = d(i);
        }
        z.setbounds(0, 0, 0, 0);
        runs = runs+1;
        if( !evd::smatrixtdevdi<Precision>(lambda, e, n, 2, i1, i2, z) )
        {
            failc = failc+1;
            return;
        }
        m = i2-i1+1;
        for(k=0; k<=m-1; k++)
        {
            serrors = serrors || amp::abs<Precision>(lambda(k)-lambdaref(i1+k))>threshold;
        }
        if( distvals )
        {
            for(j=0; j<=m-1; j++)
            {
                v = amp::vdotproduct(z.getcolumn(j, 0, n-1), zref.getcolumn(i1+j, 0, n-1));
                if( v<0 )
                {
                    amp::vmul(z.getcolumn(j, 0, n-1), -1);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    serrors = serrors || amp::abs<Precision>(z(i,j)-zref(i,i1+j))>threshold;
                }
            }
        }
    }


    /*************************************************************************
    Non-symmetric problem
    *************************************************************************/
    template<unsigned int Precision>
    void testnsevdproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision> threshold,
        bool& nserrors,
        int& failc,
        int& runs)
    {
        amp::ampf<Precision> mx;
        int i;
        int j;
        int k;
        int vjob;
        bool needl;
        bool needr;
        ap::template_1d_array< amp::ampf<Precision> > wr0;
        ap::template_1d_array< amp::ampf<Precision> > wi0;
        ap::template_1d_array< amp::ampf<Precision> > wr1;
        ap::template_1d_array< amp::ampf<Precision> > wi1;
        ap::template_1d_array< amp::ampf<Precision> > wr0s;
        ap::template_1d_array< amp::ampf<Precision> > wi0s;
        ap::template_1d_array< amp::ampf<Precision> > wr1s;
        ap::template_1d_array< amp::ampf<Precision> > wi1s;
        ap::template_2d_array< amp::ampf<Precision> > vl;
        ap::template_2d_array< amp::ampf<Precision> > vr;
        ap::template_1d_array< amp::ampf<Precision> > vec1r;
        ap::template_1d_array< amp::ampf<Precision> > vec1i;
        ap::template_1d_array< amp::ampf<Precision> > vec2r;
        ap::template_1d_array< amp::ampf<Precision> > vec2i;
        ap::template_1d_array< amp::ampf<Precision> > vec3r;
        ap::template_1d_array< amp::ampf<Precision> > vec3i;
        amp::ampf<Precision> curwr;
        amp::ampf<Precision> curwi;
        amp::ampf<Precision> vt;
        amp::ampf<Precision> tmp;


        vec1r.setbounds(0, n-1);
        vec2r.setbounds(0, n-1);
        vec3r.setbounds(0, n-1);
        vec1i.setbounds(0, n-1);
        vec2i.setbounds(0, n-1);
        vec3i.setbounds(0, n-1);
        wr0s.setbounds(0, n-1);
        wr1s.setbounds(0, n-1);
        wi0s.setbounds(0, n-1);
        wi1s.setbounds(0, n-1);
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( amp::abs<Precision>(a(i,j))>mx )
                {
                    mx = amp::abs<Precision>(a(i,j));
                }
            }
        }
        if( mx==0 )
        {
            mx = 1;
        }
        
        //
        // Load values-only
        //
        runs = runs+1;
        if( !evd::rmatrixevd<Precision>(a, n, 0, wr0, wi0, vl, vr) )
        {
            failc = failc+1;
            return;
        }
        
        //
        // Test different jobs
        //
        for(vjob=1; vjob<=3; vjob++)
        {
            needr = vjob==1 || vjob==3;
            needl = vjob==2 || vjob==3;
            runs = runs+1;
            if( !evd::rmatrixevd<Precision>(a, n, vjob, wr1, wi1, vl, vr) )
            {
                failc = failc+1;
                return;
            }
            
            //
            // Test values:
            // 1. sort by real part
            // 2. test
            //
            amp::vmove(wr0s.getvector(0, n-1), wr0.getvector(0, n-1));
            amp::vmove(wi0s.getvector(0, n-1), wi0.getvector(0, n-1));
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-2-i; j++)
                {
                    if( wr0s(j)>wr0s(j+1) )
                    {
                        tmp = wr0s(j);
                        wr0s(j) = wr0s(j+1);
                        wr0s(j+1) = tmp;
                        tmp = wi0s(j);
                        wi0s(j) = wi0s(j+1);
                        wi0s(j+1) = tmp;
                    }
                }
            }
            amp::vmove(wr1s.getvector(0, n-1), wr1.getvector(0, n-1));
            amp::vmove(wi1s.getvector(0, n-1), wi1.getvector(0, n-1));
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-2-i; j++)
                {
                    if( wr1s(j)>wr1s(j+1) )
                    {
                        tmp = wr1s(j);
                        wr1s(j) = wr1s(j+1);
                        wr1s(j+1) = tmp;
                        tmp = wi1s(j);
                        wi1s(j) = wi1s(j+1);
                        wi1s(j+1) = tmp;
                    }
                }
            }
            for(i=0; i<=n-1; i++)
            {
                nserrors = nserrors || amp::abs<Precision>(wr0s(i)-wr1s(i))>threshold;
                nserrors = nserrors || amp::abs<Precision>(wi0s(i)-wi1s(i))>threshold;
            }
            
            //
            // Test right vectors
            //
            if( needr )
            {
                k = 0;
                while( k<=n-1 )
                {
                    if( wi1(k)==0 )
                    {
                        amp::vmove(vec1r.getvector(0, n-1), vr.getcolumn(k, 0, n-1));
                        for(i=0; i<=n-1; i++)
                        {
                            vec1i(i) = 0;
                        }
                        curwr = wr1(k);
                        curwi = 0;
                    }
                    if( wi1(k)>0 )
                    {
                        amp::vmove(vec1r.getvector(0, n-1), vr.getcolumn(k, 0, n-1));
                        amp::vmove(vec1i.getvector(0, n-1), vr.getcolumn(k+1, 0, n-1));
                        curwr = wr1(k);
                        curwi = wi1(k);
                    }
                    if( wi1(k)<0 )
                    {
                        amp::vmove(vec1r.getvector(0, n-1), vr.getcolumn(k-1, 0, n-1));
                        amp::vmoveneg(vec1i.getvector(0, n-1), vr.getcolumn(k, 0, n-1));
                        curwr = wr1(k);
                        curwi = wi1(k);
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        vt = amp::vdotproduct(a.getrow(i, 0, n-1), vec1r.getvector(0, n-1));
                        vec2r(i) = vt;
                        vt = amp::vdotproduct(a.getrow(i, 0, n-1), vec1i.getvector(0, n-1));
                        vec2i(i) = vt;
                    }
                    amp::vmove(vec3r.getvector(0, n-1), vec1r.getvector(0, n-1), curwr);
                    amp::vsub(vec3r.getvector(0, n-1), vec1i.getvector(0, n-1), curwi);
                    amp::vmove(vec3i.getvector(0, n-1), vec1r.getvector(0, n-1), curwi);
                    amp::vadd(vec3i.getvector(0, n-1), vec1i.getvector(0, n-1), curwr);
                    for(i=0; i<=n-1; i++)
                    {
                        nserrors = nserrors || amp::abs<Precision>(vec2r(i)-vec3r(i))>threshold;
                        nserrors = nserrors || amp::abs<Precision>(vec2i(i)-vec3i(i))>threshold;
                    }
                    k = k+1;
                }
            }
            
            //
            // Test left vectors
            //
            if( needl )
            {
                k = 0;
                while( k<=n-1 )
                {
                    if( wi1(k)==0 )
                    {
                        amp::vmove(vec1r.getvector(0, n-1), vl.getcolumn(k, 0, n-1));
                        for(i=0; i<=n-1; i++)
                        {
                            vec1i(i) = 0;
                        }
                        curwr = wr1(k);
                        curwi = 0;
                    }
                    if( wi1(k)>0 )
                    {
                        amp::vmove(vec1r.getvector(0, n-1), vl.getcolumn(k, 0, n-1));
                        amp::vmove(vec1i.getvector(0, n-1), vl.getcolumn(k+1, 0, n-1));
                        curwr = wr1(k);
                        curwi = wi1(k);
                    }
                    if( wi1(k)<0 )
                    {
                        amp::vmove(vec1r.getvector(0, n-1), vl.getcolumn(k-1, 0, n-1));
                        amp::vmoveneg(vec1i.getvector(0, n-1), vl.getcolumn(k, 0, n-1));
                        curwr = wr1(k);
                        curwi = wi1(k);
                    }
                    for(j=0; j<=n-1; j++)
                    {
                        vt = amp::vdotproduct(vec1r.getvector(0, n-1), a.getcolumn(j, 0, n-1));
                        vec2r(j) = vt;
                        vt = amp::vdotproduct(vec1i.getvector(0, n-1), a.getcolumn(j, 0, n-1));
                        vec2i(j) = -vt;
                    }
                    amp::vmove(vec3r.getvector(0, n-1), vec1r.getvector(0, n-1), curwr);
                    amp::vadd(vec3r.getvector(0, n-1), vec1i.getvector(0, n-1), curwi);
                    amp::vmove(vec3i.getvector(0, n-1), vec1r.getvector(0, n-1), curwi);
                    amp::vsub(vec3i.getvector(0, n-1), vec1i.getvector(0, n-1), curwr);
                    for(i=0; i<=n-1; i++)
                    {
                        nserrors = nserrors || amp::abs<Precision>(vec2r(i)-vec3r(i))>threshold;
                        nserrors = nserrors || amp::abs<Precision>(vec2i(i)-vec3i(i))>threshold;
                    }
                    k = k+1;
                }
            }
        }
    }


    /*************************************************************************
    Testing EVD subroutines for one N

    NOTES:
    * BIThreshold is a threshold for bisection-and-inverse-iteration subroutines.
      special threshold is needed because these subroutines may have much more
      larger error than QR-based algorithms.
    *************************************************************************/
    template<unsigned int Precision>
    void testevdset(const int& n,
        const amp::ampf<Precision>& threshold,
        const amp::ampf<Precision>& bithreshold,
        int& failc,
        int& runs,
        bool& nserrors,
        bool& serrors,
        bool& herrors,
        bool& tderrors,
        bool& sbierrors,
        bool& hbierrors,
        bool& tdbierrors)
    {
        ap::template_2d_array< amp::ampf<Precision> > ra;
        ap::template_2d_array< amp::ampf<Precision> > ral;
        ap::template_2d_array< amp::ampf<Precision> > rau;
        ap::template_2d_array< amp::campf<Precision> > ca;
        ap::template_2d_array< amp::campf<Precision> > cal;
        ap::template_2d_array< amp::campf<Precision> > cau;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > e;
        int pass;
        int i;
        int j;
        int mkind;


        
        //
        // Test symmetric problems
        //
        
        //
        // Test symmetric problem: zero, random, sparse matrices.
        //
        ra.setlength(n, n);
        ral.setlength(n, n);
        rau.setlength(n, n);
        ca.setlength(n, n);
        cal.setlength(n, n);
        cau.setlength(n, n);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                ra(i,j) = 0;
                ca(i,j) = 0;
            }
        }
        rmatrixsymmetricsplit<Precision>(ra, n, ral, rau);
        cmatrixhermitiansplit<Precision>(ca, n, cal, cau);
        testsevdproblem<Precision>(ra, ral, rau, n, threshold, serrors, failc, runs);
        testhevdproblem<Precision>(ca, cal, cau, n, threshold, herrors, failc, runs);
        testsevdbiproblem<Precision>(ra, ral, rau, n, false, bithreshold, sbierrors, failc, runs);
        testhevdbiproblem<Precision>(ca, cal, cau, n, false, bithreshold, hbierrors, failc, runs);
        for(i=0; i<=n-1; i++)
        {
            for(j=i+1; j<=n-1; j++)
            {
                ra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                ca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                ca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                ra(j,i) = ra(i,j);
                ca(j,i) = amp::conj<Precision>(ca(i,j));
            }
            ra(i,i) = 2*amp::ampf<Precision>::getRandom()-1;
            ca(i,i) = 2*amp::ampf<Precision>::getRandom()-1;
        }
        rmatrixsymmetricsplit<Precision>(ra, n, ral, rau);
        cmatrixhermitiansplit<Precision>(ca, n, cal, cau);
        testsevdproblem<Precision>(ra, ral, rau, n, threshold, serrors, failc, runs);
        testhevdproblem<Precision>(ca, cal, cau, n, threshold, herrors, failc, runs);
        testsevdbiproblem<Precision>(ra, ral, rau, n, true, bithreshold, sbierrors, failc, runs);
        testhevdbiproblem<Precision>(ca, cal, cau, n, true, bithreshold, hbierrors, failc, runs);
        rmatrixfillsparsea<Precision>(ra, n, n, amp::ampf<Precision>("0.995"));
        cmatrixfillsparsea<Precision>(ca, n, n, amp::ampf<Precision>("0.995"));
        for(i=0; i<=n-1; i++)
        {
            for(j=i+1; j<=n-1; j++)
            {
                ra(j,i) = ra(i,j);
                ca(j,i) = amp::conj<Precision>(ca(i,j));
            }
            ca(i,i).y = 0;
        }
        rmatrixsymmetricsplit<Precision>(ra, n, ral, rau);
        cmatrixhermitiansplit<Precision>(ca, n, cal, cau);
        testsevdproblem<Precision>(ra, ral, rau, n, threshold, serrors, failc, runs);
        testhevdproblem<Precision>(ca, cal, cau, n, threshold, herrors, failc, runs);
        testsevdbiproblem<Precision>(ra, ral, rau, n, false, bithreshold, sbierrors, failc, runs);
        testhevdbiproblem<Precision>(ca, cal, cau, n, false, bithreshold, hbierrors, failc, runs);
        
        //
        // testing tridiagonal problems
        //
        for(mkind=0; mkind<=4; mkind++)
        {
            d.setlength(n);
            if( n>1 )
            {
                e.setlength(n-1);
            }
            if( mkind==0 )
            {
                
                //
                // Zero matrix
                //
                for(i=0; i<=n-1; i++)
                {
                    d(i) = 0;
                }
                for(i=0; i<=n-2; i++)
                {
                    e(i) = 0;
                }
            }
            if( mkind==1 )
            {
                
                //
                // Diagonal matrix
                //
                for(i=0; i<=n-1; i++)
                {
                    d(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                for(i=0; i<=n-2; i++)
                {
                    e(i) = 0;
                }
            }
            if( mkind==2 )
            {
                
                //
                // Off-diagonal matrix
                //
                for(i=0; i<=n-1; i++)
                {
                    d(i) = 0;
                }
                for(i=0; i<=n-2; i++)
                {
                    e(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            if( mkind==3 )
            {
                
                //
                // Dense matrix with blocks
                //
                for(i=0; i<=n-1; i++)
                {
                    d(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                for(i=0; i<=n-2; i++)
                {
                    e(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                j = 1;
                i = 2;
                while( j<=n-2 )
                {
                    e(j) = 0;
                    j = j+i;
                    i = i+1;
                }
            }
            if( mkind==4 )
            {
                
                //
                // dense matrix
                //
                for(i=0; i<=n-1; i++)
                {
                    d(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                for(i=0; i<=n-2; i++)
                {
                    e(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            testtdevdproblem<Precision>(d, e, n, threshold, tderrors, failc, runs);
            testtdevdbiproblem<Precision>(d, e, n, mkind==1 || mkind==2 || mkind==4, bithreshold, tdbierrors, failc, runs);
        }
        
        //
        // Test non-symmetric problems
        //
        
        //
        // Test non-symmetric problems: zero, random, sparse matrices.
        //
        ra.setlength(n, n);
        ca.setlength(n, n);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                ra(i,j) = 0;
                ca(i,j) = 0;
            }
        }
        testnsevdproblem<Precision>(ra, n, threshold, nserrors, failc, runs);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                ra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                ca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                ca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
            }
        }
        testnsevdproblem<Precision>(ra, n, threshold, nserrors, failc, runs);
        rmatrixfillsparsea<Precision>(ra, n, n, amp::ampf<Precision>("0.995"));
        cmatrixfillsparsea<Precision>(ca, n, n, amp::ampf<Precision>("0.995"));
        testnsevdproblem<Precision>(ra, n, threshold, nserrors, failc, runs);
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testevdunit_test_silent()
    {
        bool result;


        result = testevd<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testevdunit_test()
    {
        bool result;


        result = testevd<Precision>(false);
        return result;
    }
} // namespace

#endif
