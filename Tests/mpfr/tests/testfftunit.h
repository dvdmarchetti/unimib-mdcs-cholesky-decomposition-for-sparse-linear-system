
#ifndef _testfftunit_h
#define _testfftunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "ftbase.h"
#include "fft.h"
namespace testfftunit
{
    template<unsigned int Precision>
    bool testfft(bool silent);
    template<unsigned int Precision>
    void reffftc1d(ap::template_1d_array< amp::campf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    void reffftc1dinv(ap::template_1d_array< amp::campf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    void refinternalcfft(ap::template_1d_array< amp::ampf<Precision> >& a,
        int nn,
        bool inversefft);
    template<unsigned int Precision>
    void refinternalrfft(const ap::template_1d_array< amp::ampf<Precision> >& a,
        int nn,
        ap::template_1d_array< amp::campf<Precision> >& f);
    template<unsigned int Precision>
    bool testfftunit_test_silent();
    template<unsigned int Precision>
    bool testfftunit_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testfft(bool silent)
    {
        bool result;
        int n;
        int i;
        int k;
        ap::template_1d_array< amp::campf<Precision> > a1;
        ap::template_1d_array< amp::campf<Precision> > a2;
        ap::template_1d_array< amp::campf<Precision> > a3;
        ap::template_1d_array< amp::ampf<Precision> > r1;
        ap::template_1d_array< amp::ampf<Precision> > r2;
        ap::template_1d_array< amp::ampf<Precision> > buf;
        ftbase::ftplan<Precision> plan;
        int maxn;
        amp::ampf<Precision> bidierr;
        amp::ampf<Precision> bidirerr;
        amp::ampf<Precision> referr;
        amp::ampf<Precision> refrerr;
        amp::ampf<Precision> reinterr;
        amp::ampf<Precision> errtol;
        bool referrors;
        bool bidierrors;
        bool refrerrors;
        bool bidirerrors;
        bool reinterrors;
        bool waserrors;


        maxn = 128;
        errtol = 100000*amp::pow<Precision>(amp::ampf<Precision>(maxn), amp::ampf<Precision>(3)/amp::ampf<Precision>(2))*amp::ampf<Precision>::getAlgoPascalEpsilon();
        bidierrors = false;
        referrors = false;
        bidirerrors = false;
        refrerrors = false;
        reinterrors = false;
        waserrors = false;
        
        //
        // Test bi-directional error: norm(x-invFFT(FFT(x)))
        //
        bidierr = 0;
        bidirerr = 0;
        for(n=1; n<=maxn; n++)
        {
            
            //
            // Complex FFT/invFFT
            //
            a1.setlength(n);
            a2.setlength(n);
            a3.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                a1(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                a1(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                a2(i) = a1(i);
                a3(i) = a1(i);
            }
            fft::fftc1d<Precision>(a2, n);
            fft::fftc1dinv<Precision>(a2, n);
            fft::fftc1dinv<Precision>(a3, n);
            fft::fftc1d<Precision>(a3, n);
            for(i=0; i<=n-1; i++)
            {
                bidierr = amp::maximum<Precision>(bidierr, amp::abscomplex<Precision>(a1(i)-a2(i)));
                bidierr = amp::maximum<Precision>(bidierr, amp::abscomplex<Precision>(a1(i)-a3(i)));
            }
            
            //
            // Real
            //
            r1.setlength(n);
            r2.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                r1(i) = 2*amp::ampf<Precision>::getRandom()-1;
                r2(i) = r1(i);
            }
            fft::fftr1d<Precision>(r2, n, a1);
            amp::vmul(r2.getvector(0, n-1), 0);
            fft::fftr1dinv<Precision>(a1, n, r2);
            for(i=0; i<=n-1; i++)
            {
                bidirerr = amp::maximum<Precision>(bidirerr, amp::abscomplex<Precision>(r1(i)-r2(i)));
            }
        }
        bidierrors = bidierrors || bidierr>errtol;
        bidirerrors = bidirerrors || bidirerr>errtol;
        
        //
        // Test against reference O(N^2) implementation
        //
        referr = 0;
        refrerr = 0;
        for(n=1; n<=maxn; n++)
        {
            
            //
            // Complex FFT
            //
            a1.setlength(n);
            a2.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                a1(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                a1(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                a2(i) = a1(i);
            }
            fft::fftc1d<Precision>(a1, n);
            reffftc1d<Precision>(a2, n);
            for(i=0; i<=n-1; i++)
            {
                referr = amp::maximum<Precision>(referr, amp::abscomplex<Precision>(a1(i)-a2(i)));
            }
            
            //
            // Complex inverse FFT
            //
            a1.setlength(n);
            a2.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                a1(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                a1(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                a2(i) = a1(i);
            }
            fft::fftc1dinv<Precision>(a1, n);
            reffftc1dinv<Precision>(a2, n);
            for(i=0; i<=n-1; i++)
            {
                referr = amp::maximum<Precision>(referr, amp::abscomplex<Precision>(a1(i)-a2(i)));
            }
            
            //
            // Real forward/inverse FFT:
            // * calculate and check forward FFT
            // * use precalculated FFT to check backward FFT
            //   fill unused parts of frequencies array with random numbers
            //   to ensure that they are not really used
            //
            r1.setlength(n);
            r2.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                r1(i) = 2*amp::ampf<Precision>::getRandom()-1;
                r2(i) = r1(i);
            }
            fft::fftr1d<Precision>(r1, n, a1);
            refinternalrfft<Precision>(r2, n, a2);
            for(i=0; i<=n-1; i++)
            {
                refrerr = amp::maximum<Precision>(refrerr, amp::abscomplex<Precision>(a1(i)-a2(i)));
            }
            a3.setlength(amp::floor<Precision>(amp::ampf<Precision>(n)/amp::ampf<Precision>(2))+1);
            for(i=0; i<=amp::floor<Precision>(amp::ampf<Precision>(n)/amp::ampf<Precision>(2)); i++)
            {
                a3(i) = a2(i);
            }
            a3(0).y = 2*amp::ampf<Precision>::getRandom()-1;
            if( n%2==0 )
            {
                a3(amp::floor<Precision>(amp::ampf<Precision>(n)/amp::ampf<Precision>(2))).y = 2*amp::ampf<Precision>::getRandom()-1;
            }
            for(i=0; i<=n-1; i++)
            {
                r1(i) = 0;
            }
            fft::fftr1dinv<Precision>(a3, n, r1);
            for(i=0; i<=n-1; i++)
            {
                refrerr = amp::maximum<Precision>(refrerr, amp::abs<Precision>(r2(i)-r1(i)));
            }
        }
        referrors = referrors || referr>errtol;
        refrerrors = refrerrors || refrerr>errtol;
        
        //
        // test internal real even FFT
        //
        reinterr = 0;
        for(k=1; k<=maxn/2; k++)
        {
            n = 2*k;
            
            //
            // Real forward FFT
            //
            r1.setlength(n);
            r2.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                r1(i) = 2*amp::ampf<Precision>::getRandom()-1;
                r2(i) = r1(i);
            }
            ftbase::ftbasegeneratecomplexfftplan<Precision>(n/2, plan);
            buf.setlength(n);
            fft::fftr1dinternaleven<Precision>(r1, n, buf, plan);
            refinternalrfft<Precision>(r2, n, a2);
            reinterr = amp::maximum<Precision>(reinterr, amp::abs<Precision>(r1(0)-a2(0).x));
            reinterr = amp::maximum<Precision>(reinterr, amp::abs<Precision>(r1(1)-a2(n/2).x));
            for(i=1; i<=n/2-1; i++)
            {
                reinterr = amp::maximum<Precision>(reinterr, amp::abs<Precision>(r1(2*i+0)-a2(i).x));
                reinterr = amp::maximum<Precision>(reinterr, amp::abs<Precision>(r1(2*i+1)-a2(i).y));
            }
            
            //
            // Real backward FFT
            //
            r1.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                r1(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            a2.setlength(amp::floor<Precision>(amp::ampf<Precision>(n)/amp::ampf<Precision>(2))+1);
            a2(0) = r1(0);
            for(i=1; i<=amp::floor<Precision>(amp::ampf<Precision>(n)/amp::ampf<Precision>(2))-1; i++)
            {
                a2(i).x = r1(2*i+0);
                a2(i).y = r1(2*i+1);
            }
            a2(amp::floor<Precision>(amp::ampf<Precision>(n)/amp::ampf<Precision>(2))) = r1(1);
            ftbase::ftbasegeneratecomplexfftplan<Precision>(n/2, plan);
            buf.setlength(n);
            fft::fftr1dinvinternaleven<Precision>(r1, n, buf, plan);
            fft::fftr1dinv<Precision>(a2, n, r2);
            for(i=0; i<=n-1; i++)
            {
                reinterr = amp::maximum<Precision>(reinterr, amp::abs<Precision>(r1(i)-r2(i)));
            }
        }
        reinterrors = reinterrors || reinterr>errtol;
        
        //
        // end
        //
        waserrors = bidierrors || bidirerrors || referrors || refrerrors || reinterrors;
        if( !silent )
        {
            printf("TESTING FFT\n");
            printf("FINAL RESULT:                             ");
            if( waserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* BI-DIRECTIONAL COMPLEX TEST:            ");
            if( bidierrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* AGAINST REFERENCE COMPLEX FFT:          ");
            if( referrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* BI-DIRECTIONAL REAL TEST:               ");
            if( bidirerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* AGAINST REFERENCE REAL FFT:             ");
            if( refrerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* INTERNAL EVEN FFT:                      ");
            if( reinterrors )
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
    Reference FFT
    *************************************************************************/
    template<unsigned int Precision>
    void reffftc1d(ap::template_1d_array< amp::campf<Precision> >& a,
        int n)
    {
        ap::template_1d_array< amp::ampf<Precision> > buf;
        int i;


        ap::ap_error::make_assertion(n>0);
        buf.setlength(2*n);
        for(i=0; i<=n-1; i++)
        {
            buf(2*i+0) = a(i).x;
            buf(2*i+1) = a(i).y;
        }
        refinternalcfft<Precision>(buf, n, false);
        for(i=0; i<=n-1; i++)
        {
            a(i).x = buf(2*i+0);
            a(i).y = buf(2*i+1);
        }
    }


    /*************************************************************************
    Reference inverse FFT
    *************************************************************************/
    template<unsigned int Precision>
    void reffftc1dinv(ap::template_1d_array< amp::campf<Precision> >& a,
        int n)
    {
        ap::template_1d_array< amp::ampf<Precision> > buf;
        int i;


        ap::ap_error::make_assertion(n>0);
        buf.setlength(2*n);
        for(i=0; i<=n-1; i++)
        {
            buf(2*i+0) = a(i).x;
            buf(2*i+1) = a(i).y;
        }
        refinternalcfft<Precision>(buf, n, true);
        for(i=0; i<=n-1; i++)
        {
            a(i).x = buf(2*i+0);
            a(i).y = buf(2*i+1);
        }
    }


    /*************************************************************************
    Internal complex FFT stub.
    Uses straightforward formula with O(N^2) complexity.
    *************************************************************************/
    template<unsigned int Precision>
    void refinternalcfft(ap::template_1d_array< amp::ampf<Precision> >& a,
        int nn,
        bool inversefft)
    {
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        int i;
        int j;
        int k;
        amp::ampf<Precision> hre;
        amp::ampf<Precision> him;
        amp::ampf<Precision> c;
        amp::ampf<Precision> s;
        amp::ampf<Precision> re;
        amp::ampf<Precision> im;


        tmp.setbounds(0, 2*nn-1);
        if( !inversefft )
        {
            for(i=0; i<=nn-1; i++)
            {
                hre = 0;
                him = 0;
                for(k=0; k<=nn-1; k++)
                {
                    re = a(2*k);
                    im = a(2*k+1);
                    c = amp::cos<Precision>(-2*amp::pi<Precision>()*k*i/nn);
                    s = amp::sin<Precision>(-2*amp::pi<Precision>()*k*i/nn);
                    hre = hre+c*re-s*im;
                    him = him+c*im+s*re;
                }
                tmp(2*i) = hre;
                tmp(2*i+1) = him;
            }
            for(i=0; i<=2*nn-1; i++)
            {
                a(i) = tmp(i);
            }
        }
        else
        {
            for(k=0; k<=nn-1; k++)
            {
                hre = 0;
                him = 0;
                for(i=0; i<=nn-1; i++)
                {
                    re = a(2*i);
                    im = a(2*i+1);
                    c = amp::cos<Precision>(2*amp::pi<Precision>()*k*i/nn);
                    s = amp::sin<Precision>(2*amp::pi<Precision>()*k*i/nn);
                    hre = hre+c*re-s*im;
                    him = him+c*im+s*re;
                }
                tmp(2*k) = hre/nn;
                tmp(2*k+1) = him/nn;
            }
            for(i=0; i<=2*nn-1; i++)
            {
                a(i) = tmp(i);
            }
        }
    }


    /*************************************************************************
    Internal real FFT stub.
    Uses straightforward formula with O(N^2) complexity.
    *************************************************************************/
    template<unsigned int Precision>
    void refinternalrfft(const ap::template_1d_array< amp::ampf<Precision> >& a,
        int nn,
        ap::template_1d_array< amp::campf<Precision> >& f)
    {
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        int i;


        tmp.setbounds(0, 2*nn-1);
        for(i=0; i<=nn-1; i++)
        {
            tmp(2*i) = a(i);
            tmp(2*i+1) = 0;
        }
        refinternalcfft<Precision>(tmp, nn, false);
        f.setlength(nn);
        for(i=0; i<=nn-1; i++)
        {
            f(i).x = tmp(2*i+0);
            f(i).y = tmp(2*i+1);
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testfftunit_test_silent()
    {
        bool result;


        result = testfft<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testfftunit_test()
    {
        bool result;


        result = testfft<Precision>(false);
        return result;
    }
} // namespace

#endif
