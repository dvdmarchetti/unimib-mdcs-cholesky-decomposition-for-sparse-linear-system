
#ifndef _testcorrunit_h
#define _testcorrunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "ftbase.h"
#include "fft.h"
#include "conv.h"
#include "corr.h"
namespace testcorrunit
{
    template<unsigned int Precision>
    bool testcorr(bool silent);
    template<unsigned int Precision>
    void refcorrc1d(const ap::template_1d_array< amp::campf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::campf<Precision> >& r);
    template<unsigned int Precision>
    void refcorrc1dcircular(const ap::template_1d_array< amp::campf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::campf<Precision> >& r);
    template<unsigned int Precision>
    void refcorrr1d(const ap::template_1d_array< amp::ampf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::ampf<Precision> >& r);
    template<unsigned int Precision>
    void refcorrr1dcircular(const ap::template_1d_array< amp::ampf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::ampf<Precision> >& r);
    template<unsigned int Precision>
    void refconvc1d(const ap::template_1d_array< amp::campf<Precision> >& a,
        int m,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& r);
    template<unsigned int Precision>
    void refconvc1dcircular(const ap::template_1d_array< amp::campf<Precision> >& a,
        int m,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& r);
    template<unsigned int Precision>
    void refconvr1d(const ap::template_1d_array< amp::ampf<Precision> >& a,
        int m,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& r);
    template<unsigned int Precision>
    void refconvr1dcircular(const ap::template_1d_array< amp::ampf<Precision> >& a,
        int m,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& r);
    template<unsigned int Precision>
    bool testcorrunit_test_silent();
    template<unsigned int Precision>
    bool testcorrunit_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testcorr(bool silent)
    {
        bool result;
        int m;
        int n;
        int i;
        ap::template_1d_array< amp::ampf<Precision> > ra;
        ap::template_1d_array< amp::ampf<Precision> > rb;
        ap::template_1d_array< amp::ampf<Precision> > rr1;
        ap::template_1d_array< amp::ampf<Precision> > rr2;
        ap::template_1d_array< amp::campf<Precision> > ca;
        ap::template_1d_array< amp::campf<Precision> > cb;
        ap::template_1d_array< amp::campf<Precision> > cr1;
        ap::template_1d_array< amp::campf<Precision> > cr2;
        int maxn;
        amp::ampf<Precision> referr;
        amp::ampf<Precision> refrerr;
        amp::ampf<Precision> inverr;
        amp::ampf<Precision> invrerr;
        amp::ampf<Precision> errtol;
        bool referrors;
        bool refrerrors;
        bool inverrors;
        bool invrerrors;
        bool waserrors;


        maxn = 32;
        errtol = 100000*amp::pow<Precision>(amp::ampf<Precision>(maxn), amp::ampf<Precision>(3)/amp::ampf<Precision>(2))*amp::ampf<Precision>::getAlgoPascalEpsilon();
        referrors = false;
        refrerrors = false;
        inverrors = false;
        invrerrors = false;
        waserrors = false;
        
        //
        // Test against reference O(N^2) implementation.
        //
        referr = 0;
        refrerr = 0;
        for(m=1; m<=maxn; m++)
        {
            for(n=1; n<=maxn; n++)
            {
                
                //
                // Complex correlation
                //
                ca.setlength(m);
                for(i=0; i<=m-1; i++)
                {
                    ca(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                    ca(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
                cb.setlength(n);
                for(i=0; i<=n-1; i++)
                {
                    cb(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                    cb(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
                cr1.setlength(1);
                corr::corrc1d<Precision>(ca, m, cb, n, cr1);
                refcorrc1d<Precision>(ca, m, cb, n, cr2);
                for(i=0; i<=m+n-2; i++)
                {
                    referr = amp::maximum<Precision>(referr, amp::abscomplex<Precision>(cr1(i)-cr2(i)));
                }
                cr1.setlength(1);
                corr::corrc1dcircular<Precision>(ca, m, cb, n, cr1);
                refcorrc1dcircular<Precision>(ca, m, cb, n, cr2);
                for(i=0; i<=m-1; i++)
                {
                    referr = amp::maximum<Precision>(referr, amp::abscomplex<Precision>(cr1(i)-cr2(i)));
                }
                
                //
                // Real correlation
                //
                ra.setlength(m);
                for(i=0; i<=m-1; i++)
                {
                    ra(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                rb.setlength(n);
                for(i=0; i<=n-1; i++)
                {
                    rb(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                rr1.setlength(1);
                corr::corrr1d<Precision>(ra, m, rb, n, rr1);
                refcorrr1d<Precision>(ra, m, rb, n, rr2);
                for(i=0; i<=m+n-2; i++)
                {
                    refrerr = amp::maximum<Precision>(refrerr, amp::abs<Precision>(rr1(i)-rr2(i)));
                }
                rr1.setlength(1);
                corr::corrr1dcircular<Precision>(ra, m, rb, n, rr1);
                refcorrr1dcircular<Precision>(ra, m, rb, n, rr2);
                for(i=0; i<=m-1; i++)
                {
                    refrerr = amp::maximum<Precision>(refrerr, amp::abs<Precision>(rr1(i)-rr2(i)));
                }
            }
        }
        referrors = referrors || referr>errtol;
        refrerrors = refrerrors || refrerr>errtol;
        
        //
        // end
        //
        waserrors = referrors || refrerrors;
        if( !silent )
        {
            printf("TESTING CORRELATION\n");
            printf("FINAL RESULT:                             ");
            if( waserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* AGAINST REFERENCE COMPLEX CORR:         ");
            if( referrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* AGAINST REFERENCE REAL CORR:            ");
            if( refrerrors )
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
    Reference implementation
    *************************************************************************/
    template<unsigned int Precision>
    void refcorrc1d(const ap::template_1d_array< amp::campf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::campf<Precision> >& r)
    {
        int i;
        int j;
        amp::campf<Precision> v;
        ap::template_1d_array< amp::campf<Precision> > s;
        int i_;


        s.setlength(m+n-1);
        for(i_=0; i_<=n-1;i_++)
        {
            s(i_) = signal(i_);
        }
        for(i=n; i<=m+n-2; i++)
        {
            s(i) = 0;
        }
        r.setlength(m+n-1);
        for(i=0; i<=n-1; i++)
        {
            v = 0;
            for(j=0; j<=m-1; j++)
            {
                if( i+j>=n )
                {
                    break;
                }
                v = v+amp::conj<Precision>(pattern(j))*s(i+j);
            }
            r(i) = v;
        }
        for(i=1; i<=m-1; i++)
        {
            v = 0;
            for(j=i; j<=m-1; j++)
            {
                v = v+amp::conj<Precision>(pattern(j))*s(j-i);
            }
            r(m+n-1-i) = v;
        }
    }


    /*************************************************************************
    Reference implementation
    *************************************************************************/
    template<unsigned int Precision>
    void refcorrc1dcircular(const ap::template_1d_array< amp::campf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::campf<Precision> >& r)
    {
        int i;
        int j;
        amp::campf<Precision> v;


        r.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            v = 0;
            for(j=0; j<=m-1; j++)
            {
                v = v+amp::conj<Precision>(pattern(j))*signal((i+j)%n);
            }
            r(i) = v;
        }
    }


    /*************************************************************************
    Reference implementation
    *************************************************************************/
    template<unsigned int Precision>
    void refcorrr1d(const ap::template_1d_array< amp::ampf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::ampf<Precision> >& r)
    {
        int i;
        int j;
        amp::ampf<Precision> v;
        ap::template_1d_array< amp::ampf<Precision> > s;


        s.setlength(m+n-1);
        amp::vmove(s.getvector(0, n-1), signal.getvector(0, n-1));
        for(i=n; i<=m+n-2; i++)
        {
            s(i) = 0;
        }
        r.setlength(m+n-1);
        for(i=0; i<=n-1; i++)
        {
            v = 0;
            for(j=0; j<=m-1; j++)
            {
                if( i+j>=n )
                {
                    break;
                }
                v = v+pattern(j)*s(i+j);
            }
            r(i) = v;
        }
        for(i=1; i<=m-1; i++)
        {
            v = 0;
            for(j=i; j<=m-1; j++)
            {
                v = v+pattern(j)*s(-i+j);
            }
            r(m+n-1-i) = v;
        }
    }


    /*************************************************************************
    Reference implementation
    *************************************************************************/
    template<unsigned int Precision>
    void refcorrr1dcircular(const ap::template_1d_array< amp::ampf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::ampf<Precision> >& r)
    {
        int i;
        int j;
        amp::ampf<Precision> v;


        r.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            v = 0;
            for(j=0; j<=m-1; j++)
            {
                v = v+pattern(j)*signal((i+j)%n);
            }
            r(i) = v;
        }
    }


    /*************************************************************************
    Reference implementation
    *************************************************************************/
    template<unsigned int Precision>
    void refconvc1d(const ap::template_1d_array< amp::campf<Precision> >& a,
        int m,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& r)
    {
        int i;
        amp::campf<Precision> v;
        int i_;
        int i1_;


        r.setlength(m+n-1);
        for(i=0; i<=m+n-2; i++)
        {
            r(i) = 0;
        }
        for(i=0; i<=m-1; i++)
        {
            v = a(i);
            i1_ = (0) - (i);
            for(i_=i; i_<=i+n-1;i_++)
            {
                r(i_) = r(i_) + v*b(i_+i1_);
            }
        }
    }


    /*************************************************************************
    Reference implementation
    *************************************************************************/
    template<unsigned int Precision>
    void refconvc1dcircular(const ap::template_1d_array< amp::campf<Precision> >& a,
        int m,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& r)
    {
        int i1;
        int i2;
        int j2;
        ap::template_1d_array< amp::campf<Precision> > buf;
        int i_;
        int i1_;


        refconvc1d<Precision>(a, m, b, n, buf);
        r.setlength(m);
        for(i_=0; i_<=m-1;i_++)
        {
            r(i_) = buf(i_);
        }
        i1 = m;
        while( i1<=m+n-2 )
        {
            i2 = ap::minint(i1+m-1, m+n-2);
            j2 = i2-i1;
            i1_ = (i1) - (0);
            for(i_=0; i_<=j2;i_++)
            {
                r(i_) = r(i_) + buf(i_+i1_);
            }
            i1 = i1+m;
        }
    }


    /*************************************************************************
    Reference FFT
    *************************************************************************/
    template<unsigned int Precision>
    void refconvr1d(const ap::template_1d_array< amp::ampf<Precision> >& a,
        int m,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& r)
    {
        int i;
        amp::ampf<Precision> v;


        r.setlength(m+n-1);
        for(i=0; i<=m+n-2; i++)
        {
            r(i) = 0;
        }
        for(i=0; i<=m-1; i++)
        {
            v = a(i);
            amp::vadd(r.getvector(i, i+n-1), b.getvector(0, n-1), v);
        }
    }


    /*************************************************************************
    Reference implementation
    *************************************************************************/
    template<unsigned int Precision>
    void refconvr1dcircular(const ap::template_1d_array< amp::ampf<Precision> >& a,
        int m,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& r)
    {
        int i1;
        int i2;
        int j2;
        ap::template_1d_array< amp::ampf<Precision> > buf;


        refconvr1d<Precision>(a, m, b, n, buf);
        r.setlength(m);
        amp::vmove(r.getvector(0, m-1), buf.getvector(0, m-1));
        i1 = m;
        while( i1<=m+n-2 )
        {
            i2 = ap::minint(i1+m-1, m+n-2);
            j2 = i2-i1;
            amp::vadd(r.getvector(0, j2), buf.getvector(i1, i2));
            i1 = i1+m;
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testcorrunit_test_silent()
    {
        bool result;


        result = testcorr<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testcorrunit_test()
    {
        bool result;


        result = testcorr<Precision>(false);
        return result;
    }
} // namespace

#endif
