
#ifndef _testconvunit_h
#define _testconvunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "ftbase.h"
#include "fft.h"
#include "conv.h"
namespace testconvunit
{
    template<unsigned int Precision>
    bool testconv(bool silent);
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
    bool testconvunit_test_silent();
    template<unsigned int Precision>
    bool testconvunit_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testconv(bool silent)
    {
        bool result;
        int m;
        int n;
        int i;
        int rkind;
        int circkind;
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
        // Automatic ConvC1D() and different algorithms of ConvC1DX() are tested.
        //
        referr = 0;
        refrerr = 0;
        for(m=1; m<=maxn; m++)
        {
            for(n=1; n<=maxn; n++)
            {
                for(circkind=0; circkind<=1; circkind++)
                {
                    for(rkind=-3; rkind<=1; rkind++)
                    {
                        
                        //
                        // skip impossible combinations of parameters:
                        // * circular convolution, M<N, RKind<>-3 - internal subroutine does not support M<N.
                        //
                        if( circkind!=0 && m<n && rkind!=-3 )
                        {
                            continue;
                        }
                        
                        //
                        // Complex convolution
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
                        if( rkind==-3 )
                        {
                            
                            //
                            // test wrapper subroutine:
                            // * circular/non-circular
                            //
                            if( circkind==0 )
                            {
                                conv::convc1d<Precision>(ca, m, cb, n, cr1);
                            }
                            else
                            {
                                conv::convc1dcircular<Precision>(ca, m, cb, n, cr1);
                            }
                        }
                        else
                        {
                            
                            //
                            // test internal subroutine
                            //
                            if( m>=n )
                            {
                                
                                //
                                // test internal subroutine:
                                // * circular/non-circular mode
                                //
                                conv::convc1dx<Precision>(ca, m, cb, n, circkind!=0, rkind, 0, cr1);
                            }
                            else
                            {
                                
                                //
                                // test internal subroutine - circular mode only
                                //
                                ap::ap_error::make_assertion(circkind==0);
                                conv::convc1dx<Precision>(cb, n, ca, m, false, rkind, 0, cr1);
                            }
                        }
                        if( circkind==0 )
                        {
                            refconvc1d<Precision>(ca, m, cb, n, cr2);
                        }
                        else
                        {
                            refconvc1dcircular<Precision>(ca, m, cb, n, cr2);
                        }
                        if( circkind==0 )
                        {
                            for(i=0; i<=m+n-2; i++)
                            {
                                referr = amp::maximum<Precision>(referr, amp::abscomplex<Precision>(cr1(i)-cr2(i)));
                            }
                        }
                        else
                        {
                            for(i=0; i<=m-1; i++)
                            {
                                referr = amp::maximum<Precision>(referr, amp::abscomplex<Precision>(cr1(i)-cr2(i)));
                            }
                        }
                        
                        //
                        // Real convolution
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
                        if( rkind==-3 )
                        {
                            
                            //
                            // test wrapper subroutine:
                            // * circular/non-circular
                            //
                            if( circkind==0 )
                            {
                                conv::convr1d<Precision>(ra, m, rb, n, rr1);
                            }
                            else
                            {
                                conv::convr1dcircular<Precision>(ra, m, rb, n, rr1);
                            }
                        }
                        else
                        {
                            if( m>=n )
                            {
                                
                                //
                                // test internal subroutine:
                                // * circular/non-circular mode
                                //
                                conv::convr1dx<Precision>(ra, m, rb, n, circkind!=0, rkind, 0, rr1);
                            }
                            else
                            {
                                
                                //
                                // test internal subroutine - non-circular mode only
                                //
                                conv::convr1dx<Precision>(rb, n, ra, m, circkind!=0, rkind, 0, rr1);
                            }
                        }
                        if( circkind==0 )
                        {
                            refconvr1d<Precision>(ra, m, rb, n, rr2);
                        }
                        else
                        {
                            refconvr1dcircular<Precision>(ra, m, rb, n, rr2);
                        }
                        if( circkind==0 )
                        {
                            for(i=0; i<=m+n-2; i++)
                            {
                                refrerr = amp::maximum<Precision>(refrerr, amp::abs<Precision>(rr1(i)-rr2(i)));
                            }
                        }
                        else
                        {
                            for(i=0; i<=m-1; i++)
                            {
                                refrerr = amp::maximum<Precision>(refrerr, amp::abs<Precision>(rr1(i)-rr2(i)));
                            }
                        }
                    }
                }
            }
        }
        referrors = referrors || referr>errtol;
        refrerrors = refrerrors || refrerr>errtol;
        
        //
        // Test inverse convolution
        //
        inverr = 0;
        invrerr = 0;
        for(m=1; m<=maxn; m++)
        {
            for(n=1; n<=maxn; n++)
            {
                
                //
                // Complex circilar and non-circular
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
                cr2.setlength(1);
                conv::convc1d<Precision>(ca, m, cb, n, cr2);
                conv::convc1dinv<Precision>(cr2, m+n-1, cb, n, cr1);
                for(i=0; i<=m-1; i++)
                {
                    inverr = amp::maximum<Precision>(inverr, amp::abscomplex<Precision>(cr1(i)-ca(i)));
                }
                cr1.setlength(1);
                cr2.setlength(1);
                conv::convc1dcircular<Precision>(ca, m, cb, n, cr2);
                conv::convc1dcircularinv<Precision>(cr2, m, cb, n, cr1);
                for(i=0; i<=m-1; i++)
                {
                    inverr = amp::maximum<Precision>(inverr, amp::abscomplex<Precision>(cr1(i)-ca(i)));
                }
                
                //
                // Real circilar and non-circular
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
                rr2.setlength(1);
                conv::convr1d<Precision>(ra, m, rb, n, rr2);
                conv::convr1dinv<Precision>(rr2, m+n-1, rb, n, rr1);
                for(i=0; i<=m-1; i++)
                {
                    invrerr = amp::maximum<Precision>(invrerr, amp::abs<Precision>(rr1(i)-ra(i)));
                }
                rr1.setlength(1);
                rr2.setlength(1);
                conv::convr1dcircular<Precision>(ra, m, rb, n, rr2);
                conv::convr1dcircularinv<Precision>(rr2, m, rb, n, rr1);
                for(i=0; i<=m-1; i++)
                {
                    invrerr = amp::maximum<Precision>(invrerr, amp::abs<Precision>(rr1(i)-ra(i)));
                }
            }
        }
        inverrors = inverrors || inverr>errtol;
        invrerrors = invrerrors || invrerr>errtol;
        
        //
        // end
        //
        waserrors = referrors || refrerrors || inverrors || invrerrors;
        if( !silent )
        {
            printf("TESTING CONVOLUTION\n");
            printf("FINAL RESULT:                             ");
            if( waserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* AGAINST REFERENCE COMPLEX CONV:         ");
            if( referrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* AGAINST REFERENCE REAL CONV:            ");
            if( refrerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* COMPLEX INVERSE:                        ");
            if( inverrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* REAL INVERSE:                           ");
            if( invrerrors )
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
    bool testconvunit_test_silent()
    {
        bool result;


        result = testconv<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testconvunit_test()
    {
        bool result;


        result = testconv<Precision>(false);
        return result;
    }
} // namespace

#endif
