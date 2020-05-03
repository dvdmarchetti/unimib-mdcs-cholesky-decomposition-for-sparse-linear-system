
#ifndef _testhqrndunit_h
#define _testhqrndunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "hqrnd.h"
namespace testhqrndunit
{
    template<unsigned int Precision>
    void calculatemv(const ap::template_1d_array< amp::ampf<Precision> >& x,
        int n,
        amp::ampf<Precision>& mean,
        amp::ampf<Precision>& means,
        amp::ampf<Precision>& stddev,
        amp::ampf<Precision>& stddevs);
    template<unsigned int Precision>
    bool testhqrnd(bool silent);
    template<unsigned int Precision>
    void unsetstate(hqrnd::hqrndstate<Precision>& state);
    template<unsigned int Precision>
    bool testhqrndunit_test_silent();
    template<unsigned int Precision>
    bool testhqrndunit_test();


    template<unsigned int Precision>
    void calculatemv(const ap::template_1d_array< amp::ampf<Precision> >& x,
        int n,
        amp::ampf<Precision>& mean,
        amp::ampf<Precision>& means,
        amp::ampf<Precision>& stddev,
        amp::ampf<Precision>& stddevs)
    {
        int i;
        amp::ampf<Precision> v;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> variance;


        mean = 0;
        means = 1;
        stddev = 0;
        stddevs = 1;
        variance = 0;
        if( n<=1 )
        {
            return;
        }
        
        //
        // Mean
        //
        for(i=0; i<=n-1; i++)
        {
            mean = mean+x(i);
        }
        mean = mean/n;
        
        //
        // Variance (using corrected two-pass algorithm)
        //
        if( n!=1 )
        {
            v1 = 0;
            for(i=0; i<=n-1; i++)
            {
                v1 = v1+amp::sqr<Precision>(x(i)-mean);
            }
            v2 = 0;
            for(i=0; i<=n-1; i++)
            {
                v2 = v2+(x(i)-mean);
            }
            v2 = amp::sqr<Precision>(v2)/n;
            variance = (v1-v2)/(n-1);
            if( variance<0 )
            {
                variance = 0;
            }
            stddev = amp::sqrt<Precision>(variance);
        }
        
        //
        // Errors
        //
        means = stddev/amp::sqrt<Precision>(amp::ampf<Precision>(n));
        stddevs = stddev*amp::sqrt<Precision>(amp::ampf<Precision>(2))/amp::sqrt<Precision>(amp::ampf<Precision>(n-1));
    }


    template<unsigned int Precision>
    bool testhqrnd(bool silent)
    {
        bool result;
        bool waserrors;
        int samplesize;
        amp::ampf<Precision> sigmathreshold;
        int passcount;
        int n;
        int i;
        int pass;
        int s1;
        int s2;
        int i1;
        int i2;
        amp::ampf<Precision> r1;
        amp::ampf<Precision> r2;
        ap::template_1d_array< amp::ampf<Precision> > x;
        amp::ampf<Precision> mean;
        amp::ampf<Precision> means;
        amp::ampf<Precision> stddev;
        amp::ampf<Precision> stddevs;
        amp::ampf<Precision> lambda;
        bool seederrors;
        bool urerrors;
        amp::ampf<Precision> ursigmaerr;
        bool uierrors;
        amp::ampf<Precision> uisigmaerr;
        bool normerrors;
        amp::ampf<Precision> normsigmaerr;
        bool experrors;
        amp::ampf<Precision> expsigmaerr;
        hqrnd::hqrndstate<Precision> state;


        waserrors = false;
        sigmathreshold = 7;
        samplesize = 100000;
        passcount = 50;
        x.setbounds(0, samplesize-1);
        
        //
        // Test seed errors
        //
        seederrors = false;
        for(pass=1; pass<=passcount; pass++)
        {
            s1 = 1+ap::randominteger(32000);
            s2 = 1+ap::randominteger(32000);
            unsetstate<Precision>(state);
            hqrnd::hqrndseed<Precision>(s1, s2, state);
            i1 = hqrnd::hqrnduniformi<Precision>(100, state);
            unsetstate<Precision>(state);
            hqrnd::hqrndseed<Precision>(s1, s2, state);
            i2 = hqrnd::hqrnduniformi<Precision>(100, state);
            seederrors = seederrors || i1!=i2;
            unsetstate<Precision>(state);
            hqrnd::hqrndseed<Precision>(s1, s2, state);
            r1 = hqrnd::hqrnduniformr<Precision>(state);
            unsetstate<Precision>(state);
            hqrnd::hqrndseed<Precision>(s1, s2, state);
            r2 = hqrnd::hqrnduniformr<Precision>(state);
            seederrors = seederrors || r1!=r2;
        }
        
        //
        // Test HQRNDRandomize() and real uniform generator
        //
        unsetstate<Precision>(state);
        hqrnd::hqrndrandomize<Precision>(state);
        urerrors = false;
        ursigmaerr = 0;
        for(i=0; i<=samplesize-1; i++)
        {
            x(i) = hqrnd::hqrnduniformr<Precision>(state);
        }
        for(i=0; i<=samplesize-1; i++)
        {
            urerrors = urerrors || x(i)<=0 || x(i)>=1;
        }
        calculatemv<Precision>(x, samplesize, mean, means, stddev, stddevs);
        if( means!=0 )
        {
            ursigmaerr = amp::maximum<Precision>(ursigmaerr, amp::abs<Precision>((mean-amp::ampf<Precision>("0.5"))/means));
        }
        else
        {
            urerrors = true;
        }
        if( stddevs!=0 )
        {
            ursigmaerr = amp::maximum<Precision>(ursigmaerr, amp::abs<Precision>((stddev-amp::sqrt<Precision>(amp::ampf<Precision>(1)/amp::ampf<Precision>(12)))/stddevs));
        }
        else
        {
            urerrors = true;
        }
        urerrors = urerrors || ursigmaerr>sigmathreshold;
        
        //
        // Test HQRNDRandomize() and integer uniform
        //
        unsetstate<Precision>(state);
        hqrnd::hqrndrandomize<Precision>(state);
        uierrors = false;
        uisigmaerr = 0;
        for(n=2; n<=10; n++)
        {
            for(i=0; i<=samplesize-1; i++)
            {
                x(i) = hqrnd::hqrnduniformi<Precision>(n, state);
            }
            for(i=0; i<=samplesize-1; i++)
            {
                uierrors = uierrors || x(i)<0 || x(i)>=n;
            }
            calculatemv<Precision>(x, samplesize, mean, means, stddev, stddevs);
            if( means!=0 )
            {
                uisigmaerr = amp::maximum<Precision>(uisigmaerr, amp::abs<Precision>((mean-amp::ampf<Precision>("0.5")*(n-1))/means));
            }
            else
            {
                uierrors = true;
            }
            if( stddevs!=0 )
            {
                uisigmaerr = amp::maximum<Precision>(uisigmaerr, amp::abs<Precision>((stddev-amp::sqrt<Precision>((amp::sqr<Precision>(amp::ampf<Precision>(n))-1)/12))/stddevs));
            }
            else
            {
                uierrors = true;
            }
        }
        uierrors = uierrors || uisigmaerr>sigmathreshold;
        
        //
        // Special 'close-to-limit' test on uniformity of integers
        // (straightforward implementation like 'RND mod N' will return
        //  non-uniform numbers for N=2/3*LIMIT)
        //
        unsetstate<Precision>(state);
        hqrnd::hqrndrandomize<Precision>(state);
        uierrors = false;
        uisigmaerr = 0;
        n = amp::round<Precision>(amp::ampf<Precision>("2.0")/amp::ampf<Precision>("3.0")*amp::ampf<Precision>("2147483563.0"));
        for(i=0; i<=samplesize-1; i++)
        {
            x(i) = hqrnd::hqrnduniformi<Precision>(n, state);
        }
        for(i=0; i<=samplesize-1; i++)
        {
            uierrors = uierrors || x(i)<0 || x(i)>=n;
        }
        calculatemv<Precision>(x, samplesize, mean, means, stddev, stddevs);
        if( means!=0 )
        {
            uisigmaerr = amp::maximum<Precision>(uisigmaerr, amp::abs<Precision>((mean-amp::ampf<Precision>("0.5")*(n-1))/means));
        }
        else
        {
            uierrors = true;
        }
        if( stddevs!=0 )
        {
            uisigmaerr = amp::maximum<Precision>(uisigmaerr, amp::abs<Precision>((stddev-amp::sqrt<Precision>((amp::sqr<Precision>(amp::ampf<Precision>(n))-1)/12))/stddevs));
        }
        else
        {
            uierrors = true;
        }
        uierrors = uierrors || uisigmaerr>sigmathreshold;
        
        //
        // Test normal
        //
        unsetstate<Precision>(state);
        hqrnd::hqrndrandomize<Precision>(state);
        normerrors = false;
        normsigmaerr = 0;
        i = 0;
        while( i<samplesize )
        {
            hqrnd::hqrndnormal2<Precision>(state, r1, r2);
            x(i) = r1;
            if( i+1<samplesize )
            {
                x(i+1) = r2;
            }
            i = i+2;
        }
        calculatemv<Precision>(x, samplesize, mean, means, stddev, stddevs);
        if( means!=0 )
        {
            normsigmaerr = amp::maximum<Precision>(normsigmaerr, amp::abs<Precision>((mean-0)/means));
        }
        else
        {
            normerrors = true;
        }
        if( stddevs!=0 )
        {
            normsigmaerr = amp::maximum<Precision>(normsigmaerr, amp::abs<Precision>((stddev-1)/stddevs));
        }
        else
        {
            normerrors = true;
        }
        normerrors = normerrors || normsigmaerr>sigmathreshold;
        
        //
        // Test exponential
        //
        unsetstate<Precision>(state);
        hqrnd::hqrndrandomize<Precision>(state);
        experrors = false;
        expsigmaerr = 0;
        lambda = 2+5*amp::ampf<Precision>::getRandom();
        for(i=0; i<=samplesize-1; i++)
        {
            x(i) = hqrnd::hqrndexponential<Precision>(lambda, state);
        }
        for(i=0; i<=samplesize-1; i++)
        {
            uierrors = uierrors || x(i)<0;
        }
        calculatemv<Precision>(x, samplesize, mean, means, stddev, stddevs);
        if( means!=0 )
        {
            expsigmaerr = amp::maximum<Precision>(expsigmaerr, amp::abs<Precision>((mean-amp::ampf<Precision>("1.0")/lambda)/means));
        }
        else
        {
            experrors = true;
        }
        if( stddevs!=0 )
        {
            expsigmaerr = amp::maximum<Precision>(expsigmaerr, amp::abs<Precision>((stddev-amp::ampf<Precision>("1.0")/lambda)/stddevs));
        }
        else
        {
            experrors = true;
        }
        experrors = experrors || expsigmaerr>sigmathreshold;
        
        //
        // Final report
        //
        waserrors = seederrors || urerrors || uierrors || normerrors || experrors;
        if( !silent )
        {
            printf("RNG TEST\n");
            printf("SEED TEST:                               ");
            if( !seederrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("UNIFORM CONTINUOUS:                      ");
            if( !urerrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("UNIFORM INTEGER:                         ");
            if( !uierrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("NORMAL:                                  ");
            if( !normerrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("EXPONENTIAL:                             ");
            if( !experrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            if( waserrors )
            {
                printf("TEST SUMMARY: FAILED\n");
            }
            else
            {
                printf("TEST SUMMARY: PASSED\n");
            }
            printf("\n\n");
        }
        result = !waserrors;
        return result;
    }


    /*************************************************************************
    Unsets HQRNDState structure
    *************************************************************************/
    template<unsigned int Precision>
    void unsetstate(hqrnd::hqrndstate<Precision>& state)
    {
        state.s1 = 0;
        state.s2 = 0;
        state.v = 0;
        state.magicv = 0;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testhqrndunit_test_silent()
    {
        bool result;


        result = testhqrnd<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testhqrndunit_test()
    {
        bool result;


        result = testhqrnd<Precision>(false);
        return result;
    }
} // namespace

#endif
