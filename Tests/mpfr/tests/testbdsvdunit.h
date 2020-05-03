
#ifndef _testbdsvdunit_h
#define _testbdsvdunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "rotations.h"
#include "bdsvd.h"
namespace testbdsvdunit
{
    template<unsigned int Precision>
    bool testbdsvd(bool silent);
    template<unsigned int Precision>
    void fillidentity(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    void fillsparsede(ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        amp::ampf<Precision> sparcity);
    template<unsigned int Precision>
    void getbdsvderror(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::ampf<Precision> >& u,
        const ap::template_2d_array< amp::ampf<Precision> >& c,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& vt,
        amp::ampf<Precision>& materr,
        amp::ampf<Precision>& orterr,
        bool& wsorted);
    template<unsigned int Precision>
    void testbdsvdproblem(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        amp::ampf<Precision>& materr,
        amp::ampf<Precision>& orterr,
        bool& wsorted,
        bool& wfailed);
    template<unsigned int Precision>
    bool testbdsvdunit_test_silent();
    template<unsigned int Precision>
    bool testbdsvdunit_test();


    static int failcount;
    static int succcount;


    /*************************************************************************
    Testing bidiagonal SVD decomposition subroutine
    *************************************************************************/
    template<unsigned int Precision>
    bool testbdsvd(bool silent)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > e;
        ap::template_2d_array< amp::ampf<Precision> > mempty;
        int n;
        int maxn;
        int i;
        int j;
        int gpass;
        int pass;
        bool waserrors;
        bool wsorted;
        bool wfailed;
        bool failcase;
        amp::ampf<Precision> materr;
        amp::ampf<Precision> orterr;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> failthreshold;
        amp::ampf<Precision> failr;


        failcount = 0;
        succcount = 0;
        materr = 0;
        orterr = 0;
        wsorted = true;
        wfailed = false;
        waserrors = false;
        maxn = 15;
        threshold = 5*100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        failthreshold = amp::ampf<Precision>("1.0E-2");
        d.setbounds(0, maxn-1);
        e.setbounds(0, maxn-2);
        
        //
        // special case: fail matrix
        //
        n = 5;
        d(0) = -amp::ampf<Precision>("8.27448347422711894000e-01");
        d(1) = -amp::ampf<Precision>("8.16705832087160854600e-01");
        d(2) = -amp::ampf<Precision>("2.53974358904729382800e-17");
        d(3) = -amp::ampf<Precision>("1.24626684881972815700e+00");
        d(4) = -amp::ampf<Precision>("4.64744131545637651000e-01");
        e(0) = -amp::ampf<Precision>("3.25785088656270038800e-01");
        e(1) = -amp::ampf<Precision>("1.03732413708914436580e-01");
        e(2) = -amp::ampf<Precision>("9.57365642262031357700e-02");
        e(3) = -amp::ampf<Precision>("2.71564153973817390400e-01");
        failcase = bdsvd::rmatrixbdsvd<Precision>(d, e, n, true, false, mempty, 0, mempty, 0, mempty, 0);
        
        //
        // special case: zero divide matrix
        // unfixed LAPACK routine should fail on this problem
        //
        n = 7;
        d(0) = -amp::ampf<Precision>("6.96462904751731892700e-01");
        d(1) = amp::ampf<Precision>("0.00000000000000000000e+00");
        d(2) = -amp::ampf<Precision>("5.73827770385971991400e-01");
        d(3) = -amp::ampf<Precision>("6.62562624399371191700e-01");
        d(4) = amp::ampf<Precision>("5.82737148001782223600e-01");
        d(5) = amp::ampf<Precision>("3.84825263580925003300e-01");
        d(6) = amp::ampf<Precision>("9.84087420830525472200e-01");
        e(0) = -amp::ampf<Precision>("7.30307931760612871800e-02");
        e(1) = -amp::ampf<Precision>("2.30079042939542843800e-01");
        e(2) = -amp::ampf<Precision>("6.87824621739351216300e-01");
        e(3) = -amp::ampf<Precision>("1.77306437707837570600e-02");
        e(4) = amp::ampf<Precision>("1.78285126526551632000e-15");
        e(5) = -amp::ampf<Precision>("4.89434737751289969400e-02");
        bdsvd::rmatrixbdsvd<Precision>(d, e, n, true, false, mempty, 0, mempty, 0, mempty, 0);
        
        //
        // zero matrix, several cases
        //
        for(i=0; i<=maxn-1; i++)
        {
            d(i) = 0;
        }
        for(i=0; i<=maxn-2; i++)
        {
            e(i) = 0;
        }
        for(n=1; n<=maxn; n++)
        {
            testbdsvdproblem<Precision>(d, e, n, materr, orterr, wsorted, wfailed);
        }
        
        //
        // Dense matrix
        //
        for(n=1; n<=maxn; n++)
        {
            for(pass=1; pass<=10; pass++)
            {
                for(i=0; i<=maxn-1; i++)
                {
                    d(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                for(i=0; i<=maxn-2; i++)
                {
                    e(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                testbdsvdproblem<Precision>(d, e, n, materr, orterr, wsorted, wfailed);
            }
        }
        
        //
        // Sparse matrices, very sparse matrices, incredible sparse matrices
        //
        for(n=1; n<=maxn; n++)
        {
            for(pass=1; pass<=10; pass++)
            {
                fillsparsede<Precision>(d, e, n, amp::ampf<Precision>("0.5"));
                testbdsvdproblem<Precision>(d, e, n, materr, orterr, wsorted, wfailed);
                fillsparsede<Precision>(d, e, n, amp::ampf<Precision>("0.8"));
                testbdsvdproblem<Precision>(d, e, n, materr, orterr, wsorted, wfailed);
                fillsparsede<Precision>(d, e, n, amp::ampf<Precision>("0.9"));
                testbdsvdproblem<Precision>(d, e, n, materr, orterr, wsorted, wfailed);
                fillsparsede<Precision>(d, e, n, amp::ampf<Precision>("0.95"));
                testbdsvdproblem<Precision>(d, e, n, materr, orterr, wsorted, wfailed);
            }
        }
        
        //
        // report
        //
        failr = amp::ampf<Precision>(failcount)/(amp::ampf<Precision>(succcount+failcount));
        waserrors = materr>threshold || orterr>threshold || !wsorted || failr>failthreshold;
        if( !silent )
        {
            printf("TESTING BIDIAGONAL SVD DECOMPOSITION\n");
            printf("SVD decomposition error:                 %5.3le\n",
                double(amp::ampf<Precision>(materr).toDouble()));
            printf("SVD orthogonality error:                 %5.3le\n",
                double(amp::ampf<Precision>(orterr).toDouble()));
            printf("Singular values order:                   ");
            if( wsorted )
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
                printf("YES\n");
            }
            else
            {
                printf("NO\n");
                printf("Fail ratio:                              %5.3lf\n",
                    double(amp::ampf<Precision>(failr).toDouble()));
            }
            printf("Fail matrix test:                        ");
            if( !failcase )
            {
                printf("AS EXPECTED\n");
            }
            else
            {
                printf("CONVERGED (UNEXPECTED)\n");
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
    void fillidentity(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n)
    {
        int i;
        int j;


        a.setbounds(0, n-1, 0, n-1);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( i==j )
                {
                    a(i,j) = 1;
                }
                else
                {
                    a(i,j) = 0;
                }
            }
        }
    }


    template<unsigned int Precision>
    void fillsparsede(ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        amp::ampf<Precision> sparcity)
    {
        int i;
        int j;


        d.setbounds(0, n-1);
        e.setbounds(0, ap::maxint(0, n-2));
        for(i=0; i<=n-1; i++)
        {
            if( amp::ampf<Precision>::getRandom()>=sparcity )
            {
                d(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            else
            {
                d(i) = 0;
            }
        }
        for(i=0; i<=n-2; i++)
        {
            if( amp::ampf<Precision>::getRandom()>=sparcity )
            {
                e(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            else
            {
                e(i) = 0;
            }
        }
    }


    template<unsigned int Precision>
    void getbdsvderror(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::ampf<Precision> >& u,
        const ap::template_2d_array< amp::ampf<Precision> >& c,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& vt,
        amp::ampf<Precision>& materr,
        amp::ampf<Precision>& orterr,
        bool& wsorted)
    {
        int i;
        int j;
        int k;
        amp::ampf<Precision> locerr;
        amp::ampf<Precision> sm;


        
        //
        // decomposition error
        //
        locerr = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                sm = 0;
                for(k=0; k<=n-1; k++)
                {
                    sm = sm+w(k)*u(i,k)*vt(k,j);
                }
                if( isupper )
                {
                    if( i==j )
                    {
                        locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(d(i)-sm));
                    }
                    else
                    {
                        if( i==j-1 )
                        {
                            locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(e(i)-sm));
                        }
                        else
                        {
                            locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(sm));
                        }
                    }
                }
                else
                {
                    if( i==j )
                    {
                        locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(d(i)-sm));
                    }
                    else
                    {
                        if( i-1==j )
                        {
                            locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(e(j)-sm));
                        }
                        else
                        {
                            locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(sm));
                        }
                    }
                }
            }
        }
        materr = amp::maximum<Precision>(materr, locerr);
        
        //
        // check for C = U'
        // we consider it as decomposition error
        //
        locerr = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(u(i,j)-c(j,i)));
            }
        }
        materr = amp::maximum<Precision>(materr, locerr);
        
        //
        // orthogonality error
        //
        locerr = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=i; j<=n-1; j++)
            {
                sm = amp::vdotproduct(u.getcolumn(i, 0, n-1), u.getcolumn(j, 0, n-1));
                if( i!=j )
                {
                    locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(sm));
                }
                else
                {
                    locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(sm-1));
                }
                sm = amp::vdotproduct(vt.getrow(i, 0, n-1), vt.getrow(j, 0, n-1));
                if( i!=j )
                {
                    locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(sm));
                }
                else
                {
                    locerr = amp::maximum<Precision>(locerr, amp::abs<Precision>(sm-1));
                }
            }
        }
        orterr = amp::maximum<Precision>(orterr, locerr);
        
        //
        // values order error
        //
        for(i=1; i<=n-1; i++)
        {
            if( w(i)>w(i-1) )
            {
                wsorted = false;
            }
        }
    }


    template<unsigned int Precision>
    void testbdsvdproblem(const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        amp::ampf<Precision>& materr,
        amp::ampf<Precision>& orterr,
        bool& wsorted,
        bool& wfailed)
    {
        ap::template_2d_array< amp::ampf<Precision> > u;
        ap::template_2d_array< amp::ampf<Precision> > vt;
        ap::template_2d_array< amp::ampf<Precision> > c;
        ap::template_1d_array< amp::ampf<Precision> > w;
        int i;
        int j;
        int k;
        amp::ampf<Precision> v;
        amp::ampf<Precision> mx;


        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            if( amp::abs<Precision>(d(i))>mx )
            {
                mx = amp::abs<Precision>(d(i));
            }
        }
        for(i=0; i<=n-2; i++)
        {
            if( amp::abs<Precision>(e(i))>mx )
            {
                mx = amp::abs<Precision>(e(i));
            }
        }
        if( mx==0 )
        {
            mx = 1;
        }
        
        //
        // Upper BDSVD tests
        //
        w.setbounds(0, n-1);
        fillidentity<Precision>(u, n);
        fillidentity<Precision>(vt, n);
        fillidentity<Precision>(c, n);
        for(i=0; i<=n-1; i++)
        {
            w(i) = d(i);
        }
        if( !bdsvd::rmatrixbdsvd<Precision>(w, e, n, true, false, u, n, c, n, vt, n) )
        {
            failcount = failcount+1;
            wfailed = true;
            return;
        }
        getbdsvderror<Precision>(d, e, n, true, u, c, w, vt, materr, orterr, wsorted);
        fillidentity<Precision>(u, n);
        fillidentity<Precision>(vt, n);
        fillidentity<Precision>(c, n);
        for(i=0; i<=n-1; i++)
        {
            w(i) = d(i);
        }
        if( !bdsvd::rmatrixbdsvd<Precision>(w, e, n, true, true, u, n, c, n, vt, n) )
        {
            failcount = failcount+1;
            wfailed = true;
            return;
        }
        getbdsvderror<Precision>(d, e, n, true, u, c, w, vt, materr, orterr, wsorted);
        
        //
        // Lower BDSVD tests
        //
        w.setbounds(0, n-1);
        fillidentity<Precision>(u, n);
        fillidentity<Precision>(vt, n);
        fillidentity<Precision>(c, n);
        for(i=0; i<=n-1; i++)
        {
            w(i) = d(i);
        }
        if( !bdsvd::rmatrixbdsvd<Precision>(w, e, n, false, false, u, n, c, n, vt, n) )
        {
            failcount = failcount+1;
            wfailed = true;
            return;
        }
        getbdsvderror<Precision>(d, e, n, false, u, c, w, vt, materr, orterr, wsorted);
        fillidentity<Precision>(u, n);
        fillidentity<Precision>(vt, n);
        fillidentity<Precision>(c, n);
        for(i=0; i<=n-1; i++)
        {
            w(i) = d(i);
        }
        if( !bdsvd::rmatrixbdsvd<Precision>(w, e, n, false, true, u, n, c, n, vt, n) )
        {
            failcount = failcount+1;
            wfailed = true;
            return;
        }
        getbdsvderror<Precision>(d, e, n, false, u, c, w, vt, materr, orterr, wsorted);
        
        //
        // update counter
        //
        succcount = succcount+1;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testbdsvdunit_test_silent()
    {
        bool result;


        result = testbdsvd<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testbdsvdunit_test()
    {
        bool result;


        result = testbdsvd<Precision>(false);
        return result;
    }
} // namespace

#endif
