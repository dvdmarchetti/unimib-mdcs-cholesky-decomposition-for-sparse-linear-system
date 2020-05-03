
#ifndef _testfhtunit_h
#define _testfhtunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "ftbase.h"
#include "fft.h"
#include "fht.h"
namespace testfhtunit
{
    template<unsigned int Precision>
    bool testfht(bool silent);
    template<unsigned int Precision>
    void reffhtr1d(ap::template_1d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    void reffhtr1dinv(ap::template_1d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    bool testfhtunit_test_silent();
    template<unsigned int Precision>
    bool testfhtunit_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testfht(bool silent)
    {
        bool result;
        int n;
        int i;
        ap::template_1d_array< amp::ampf<Precision> > r1;
        ap::template_1d_array< amp::ampf<Precision> > r2;
        ap::template_1d_array< amp::ampf<Precision> > r3;
        int maxn;
        amp::ampf<Precision> bidierr;
        amp::ampf<Precision> referr;
        amp::ampf<Precision> errtol;
        bool referrors;
        bool bidierrors;
        bool waserrors;


        maxn = 128;
        errtol = 100000*amp::pow<Precision>(amp::ampf<Precision>(maxn), amp::ampf<Precision>(3)/amp::ampf<Precision>(2))*amp::ampf<Precision>::getAlgoPascalEpsilon();
        bidierrors = false;
        referrors = false;
        waserrors = false;
        
        //
        // Test bi-directional error: norm(x-invFHT(FHT(x)))
        //
        bidierr = 0;
        for(n=1; n<=maxn; n++)
        {
            
            //
            // FHT/invFHT
            //
            r1.setlength(n);
            r2.setlength(n);
            r3.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                r1(i) = 2*amp::ampf<Precision>::getRandom()-1;
                r2(i) = r1(i);
                r3(i) = r1(i);
            }
            fht::fhtr1d<Precision>(r2, n);
            fht::fhtr1dinv<Precision>(r2, n);
            fht::fhtr1dinv<Precision>(r3, n);
            fht::fhtr1d<Precision>(r3, n);
            for(i=0; i<=n-1; i++)
            {
                bidierr = amp::maximum<Precision>(bidierr, amp::abs<Precision>(r1(i)-r2(i)));
                bidierr = amp::maximum<Precision>(bidierr, amp::abs<Precision>(r1(i)-r3(i)));
            }
        }
        bidierrors = bidierrors || bidierr>errtol;
        
        //
        // Test against reference O(N^2) implementation
        //
        referr = 0;
        for(n=1; n<=maxn; n++)
        {
            
            //
            // FHT
            //
            r1.setlength(n);
            r2.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                r1(i) = 2*amp::ampf<Precision>::getRandom()-1;
                r2(i) = r1(i);
            }
            fht::fhtr1d<Precision>(r1, n);
            reffhtr1d<Precision>(r2, n);
            for(i=0; i<=n-1; i++)
            {
                referr = amp::maximum<Precision>(referr, amp::abs<Precision>(r1(i)-r2(i)));
            }
            
            //
            // inverse FHT
            //
            r1.setlength(n);
            r2.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                r1(i) = 2*amp::ampf<Precision>::getRandom()-1;
                r2(i) = r1(i);
            }
            fht::fhtr1dinv<Precision>(r1, n);
            reffhtr1dinv<Precision>(r2, n);
            for(i=0; i<=n-1; i++)
            {
                referr = amp::maximum<Precision>(referr, amp::abs<Precision>(r1(i)-r2(i)));
            }
        }
        referrors = referrors || referr>errtol;
        
        //
        // end
        //
        waserrors = bidierrors || referrors;
        if( !silent )
        {
            printf("TESTING FHT\n");
            printf("FINAL RESULT:                             ");
            if( waserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* BI-DIRECTIONAL TEST:                    ");
            if( bidierrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* AGAINST REFERENCE FHT:                  ");
            if( referrors )
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
    Reference FHT
    *************************************************************************/
    template<unsigned int Precision>
    void reffhtr1d(ap::template_1d_array< amp::ampf<Precision> >& a,
        int n)
    {
        ap::template_1d_array< amp::ampf<Precision> > buf;
        int i;
        int j;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(n>0);
        buf.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            v = 0;
            for(j=0; j<=n-1; j++)
            {
                v = v+a(j)*(amp::cos<Precision>(2*amp::pi<Precision>()*i*j/n)+amp::sin<Precision>(2*amp::pi<Precision>()*i*j/n));
            }
            buf(i) = v;
        }
        for(i=0; i<=n-1; i++)
        {
            a(i) = buf(i);
        }
    }


    /*************************************************************************
    Reference inverse FHT
    *************************************************************************/
    template<unsigned int Precision>
    void reffhtr1dinv(ap::template_1d_array< amp::ampf<Precision> >& a,
        int n)
    {
        int i;


        ap::ap_error::make_assertion(n>0);
        reffhtr1d<Precision>(a, n);
        for(i=0; i<=n-1; i++)
        {
            a(i) = a(i)/n;
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testfhtunit_test_silent()
    {
        bool result;


        result = testfht<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testfhtunit_test()
    {
        bool result;


        result = testfht<Precision>(false);
        return result;
    }
} // namespace

#endif
