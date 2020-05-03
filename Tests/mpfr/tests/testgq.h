
#ifndef _testgq_h
#define _testgq_h

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
#include "gammafunc.h"
#include "gq.h"
namespace testgq
{
    template<unsigned int Precision>
    bool testgqunit(bool silent);
    template<unsigned int Precision>
    void buildgausshermitequadrature(int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    amp::ampf<Precision> mapkind(int k);
    template<unsigned int Precision>
    void buildgausslegendrequadrature(int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void buildgaussjacobiquadrature(int n,
        amp::ampf<Precision> alpha,
        amp::ampf<Precision> beta,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void buildgausslaguerrequadrature(int n,
        amp::ampf<Precision> alpha,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    bool testgq_test_silent();
    template<unsigned int Precision>
    bool testgq_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testgqunit(bool silent)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > alpha;
        ap::template_1d_array< amp::ampf<Precision> > beta;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > x2;
        ap::template_1d_array< amp::ampf<Precision> > w2;
        amp::ampf<Precision> err;
        int n;
        int i;
        int info;
        int akind;
        int bkind;
        amp::ampf<Precision> alphac;
        amp::ampf<Precision> betac;
        amp::ampf<Precision> errtol;
        amp::ampf<Precision> nonstricterrtol;
        amp::ampf<Precision> stricterrtol;
        bool recerrors;
        bool specerrors;
        bool waserrors;


        recerrors = false;
        specerrors = false;
        waserrors = false;
        errtol = amp::ampf<Precision>("1.0E-12");
        nonstricterrtol = amp::ampf<Precision>("1.0E-6");
        stricterrtol = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Three tests for rec-based Gauss quadratures with known weights/nodes:
        // 1. Gauss-Legendre with N=2
        // 2. Gauss-Legendre with N=5
        // 3. Gauss-Chebyshev with N=1, 2, 4, 8, ..., 512
        //
        err = 0;
        alpha.setlength(2);
        beta.setlength(2);
        alpha(0) = 0;
        alpha(1) = 0;
        beta(1) = amp::ampf<Precision>(1)/(amp::ampf<Precision>(4*1*1-1));
        gq::gqgeneraterec<Precision>(alpha, beta, amp::ampf<Precision>("2.0"), 2, info, x, w);
        if( info>0 )
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+amp::sqrt<Precision>(amp::ampf<Precision>(3))/3));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)-amp::sqrt<Precision>(amp::ampf<Precision>(3))/3));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-1));
            for(i=0; i<=0; i++)
            {
                recerrors = recerrors || x(i)>=x(i+1);
            }
        }
        else
        {
            recerrors = true;
        }
        alpha.setlength(5);
        beta.setlength(5);
        alpha(0) = 0;
        for(i=1; i<=4; i++)
        {
            alpha(i) = 0;
            beta(i) = amp::sqr<Precision>(amp::ampf<Precision>(i))/(4*amp::sqr<Precision>(amp::ampf<Precision>(i))-1);
        }
        gq::gqgeneraterec<Precision>(alpha, beta, amp::ampf<Precision>("2.0"), 5, info, x, w);
        if( info>0 )
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+amp::sqrt<Precision>(245+14*amp::sqrt<Precision>(amp::ampf<Precision>(70)))/21));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+x(4)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)+amp::sqrt<Precision>(245-14*amp::sqrt<Precision>(amp::ampf<Precision>(70)))/21));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)+x(3)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(2)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-(322-13*amp::sqrt<Precision>(amp::ampf<Precision>(70)))/900));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-w(4)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-(322+13*amp::sqrt<Precision>(amp::ampf<Precision>(70)))/900));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-w(3)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(2)-amp::ampf<Precision>(128)/amp::ampf<Precision>(225)));
            for(i=0; i<=3; i++)
            {
                recerrors = recerrors || x(i)>=x(i+1);
            }
        }
        else
        {
            recerrors = true;
        }
        n = 1;
        while( n<=512 )
        {
            alpha.setlength(n);
            beta.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                alpha(i) = 0;
                if( i==0 )
                {
                    beta(i) = 0;
                }
                if( i==1 )
                {
                    beta(i) = amp::ampf<Precision>(1)/amp::ampf<Precision>(2);
                }
                if( i>1 )
                {
                    beta(i) = amp::ampf<Precision>(1)/amp::ampf<Precision>(4);
                }
            }
            gq::gqgeneraterec<Precision>(alpha, beta, amp::pi<Precision>(), n, info, x, w);
            if( info>0 )
            {
                for(i=0; i<=n-1; i++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(x(i)-amp::cos<Precision>(amp::pi<Precision>()*(n-i-amp::ampf<Precision>("0.5"))/n)));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(w(i)-amp::pi<Precision>()/n));
                }
                for(i=0; i<=n-2; i++)
                {
                    recerrors = recerrors || x(i)>=x(i+1);
                }
            }
            else
            {
                recerrors = true;
            }
            n = n*2;
        }
        recerrors = recerrors || err>errtol;
        
        //
        // Three tests for rec-based Gauss-Lobatto quadratures with known weights/nodes:
        // 1. Gauss-Lobatto with N=3
        // 2. Gauss-Lobatto with N=4
        // 3. Gauss-Lobatto with N=6
        //
        err = 0;
        alpha.setlength(2);
        beta.setlength(2);
        alpha(0) = 0;
        alpha(1) = 0;
        beta(0) = 0;
        beta(1) = amp::ampf<Precision>(1*1)/(amp::ampf<Precision>(4*1*1-1));
        gq::gqgenerategausslobattorec<Precision>(alpha, beta, amp::ampf<Precision>("2.0"), amp::ampf<Precision>(-1), amp::ampf<Precision>(+1), 3, info, x, w);
        if( info>0 )
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(2)-1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-amp::ampf<Precision>(1)/amp::ampf<Precision>(3)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-amp::ampf<Precision>(4)/amp::ampf<Precision>(3)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(2)-amp::ampf<Precision>(1)/amp::ampf<Precision>(3)));
            for(i=0; i<=1; i++)
            {
                recerrors = recerrors || x(i)>=x(i+1);
            }
        }
        else
        {
            recerrors = true;
        }
        alpha.setlength(3);
        beta.setlength(3);
        alpha(0) = 0;
        alpha(1) = 0;
        alpha(2) = 0;
        beta(0) = 0;
        beta(1) = amp::ampf<Precision>(1*1)/(amp::ampf<Precision>(4*1*1-1));
        beta(2) = amp::ampf<Precision>(2*2)/(amp::ampf<Precision>(4*2*2-1));
        gq::gqgenerategausslobattorec<Precision>(alpha, beta, amp::ampf<Precision>("2.0"), amp::ampf<Precision>(-1), amp::ampf<Precision>(+1), 4, info, x, w);
        if( info>0 )
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)+amp::sqrt<Precision>(amp::ampf<Precision>(5))/5));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(2)-amp::sqrt<Precision>(amp::ampf<Precision>(5))/5));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(3)-1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-amp::ampf<Precision>(1)/amp::ampf<Precision>(6)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-amp::ampf<Precision>(5)/amp::ampf<Precision>(6)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(2)-amp::ampf<Precision>(5)/amp::ampf<Precision>(6)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(3)-amp::ampf<Precision>(1)/amp::ampf<Precision>(6)));
            for(i=0; i<=2; i++)
            {
                recerrors = recerrors || x(i)>=x(i+1);
            }
        }
        else
        {
            recerrors = true;
        }
        alpha.setlength(5);
        beta.setlength(5);
        alpha(0) = 0;
        alpha(1) = 0;
        alpha(2) = 0;
        alpha(3) = 0;
        alpha(4) = 0;
        beta(0) = 0;
        beta(1) = amp::ampf<Precision>(1*1)/(amp::ampf<Precision>(4*1*1-1));
        beta(2) = amp::ampf<Precision>(2*2)/(amp::ampf<Precision>(4*2*2-1));
        beta(3) = amp::ampf<Precision>(3*3)/(amp::ampf<Precision>(4*3*3-1));
        beta(4) = amp::ampf<Precision>(4*4)/(amp::ampf<Precision>(4*4*4-1));
        gq::gqgenerategausslobattorec<Precision>(alpha, beta, amp::ampf<Precision>("2.0"), amp::ampf<Precision>(-1), amp::ampf<Precision>(+1), 6, info, x, w);
        if( info>0 )
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)+amp::sqrt<Precision>((7+2*amp::sqrt<Precision>(amp::ampf<Precision>(7)))/21)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(2)+amp::sqrt<Precision>((7-2*amp::sqrt<Precision>(amp::ampf<Precision>(7)))/21)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(3)-amp::sqrt<Precision>((7-2*amp::sqrt<Precision>(amp::ampf<Precision>(7)))/21)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(4)-amp::sqrt<Precision>((7+2*amp::sqrt<Precision>(amp::ampf<Precision>(7)))/21)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(5)-1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-amp::ampf<Precision>(1)/amp::ampf<Precision>(15)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-(14-amp::sqrt<Precision>(amp::ampf<Precision>(7)))/30));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(2)-(14+amp::sqrt<Precision>(amp::ampf<Precision>(7)))/30));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(3)-(14+amp::sqrt<Precision>(amp::ampf<Precision>(7)))/30));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(4)-(14-amp::sqrt<Precision>(amp::ampf<Precision>(7)))/30));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(5)-amp::ampf<Precision>(1)/amp::ampf<Precision>(15)));
            for(i=0; i<=4; i++)
            {
                recerrors = recerrors || x(i)>=x(i+1);
            }
        }
        else
        {
            recerrors = true;
        }
        recerrors = recerrors || err>errtol;
        
        //
        // Three tests for rec-based Gauss-Radau quadratures with known weights/nodes:
        // 1. Gauss-Radau with N=2
        // 2. Gauss-Radau with N=3
        // 3. Gauss-Radau with N=3 (another case)
        //
        err = 0;
        alpha.setlength(1);
        beta.setlength(2);
        alpha(0) = 0;
        beta(0) = 0;
        beta(1) = amp::ampf<Precision>(1*1)/(amp::ampf<Precision>(4*1*1-1));
        gq::gqgenerategaussradaurec<Precision>(alpha, beta, amp::ampf<Precision>("2.0"), amp::ampf<Precision>(-1), 2, info, x, w);
        if( info>0 )
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)-amp::ampf<Precision>(1)/amp::ampf<Precision>(3)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-amp::ampf<Precision>("0.5")));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-amp::ampf<Precision>("1.5")));
            for(i=0; i<=0; i++)
            {
                recerrors = recerrors || x(i)>=x(i+1);
            }
        }
        else
        {
            recerrors = true;
        }
        alpha.setlength(2);
        beta.setlength(3);
        alpha(0) = 0;
        alpha(1) = 0;
        for(i=0; i<=2; i++)
        {
            beta(i) = amp::sqr<Precision>(amp::ampf<Precision>(i))/(4*amp::sqr<Precision>(amp::ampf<Precision>(i))-1);
        }
        gq::gqgenerategaussradaurec<Precision>(alpha, beta, amp::ampf<Precision>("2.0"), amp::ampf<Precision>(-1), 3, info, x, w);
        if( info>0 )
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)-(1-amp::sqrt<Precision>(amp::ampf<Precision>(6)))/5));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(2)-(1+amp::sqrt<Precision>(amp::ampf<Precision>(6)))/5));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-amp::ampf<Precision>(2)/amp::ampf<Precision>(9)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-(16+amp::sqrt<Precision>(amp::ampf<Precision>(6)))/18));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(2)-(16-amp::sqrt<Precision>(amp::ampf<Precision>(6)))/18));
            for(i=0; i<=1; i++)
            {
                recerrors = recerrors || x(i)>=x(i+1);
            }
        }
        else
        {
            recerrors = true;
        }
        alpha.setlength(2);
        beta.setlength(3);
        alpha(0) = 0;
        alpha(1) = 0;
        for(i=0; i<=2; i++)
        {
            beta(i) = amp::sqr<Precision>(amp::ampf<Precision>(i))/(4*amp::sqr<Precision>(amp::ampf<Precision>(i))-1);
        }
        gq::gqgenerategaussradaurec<Precision>(alpha, beta, amp::ampf<Precision>("2.0"), amp::ampf<Precision>(+1), 3, info, x, w);
        if( info>0 )
        {
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(2)-1));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(1)+(1-amp::sqrt<Precision>(amp::ampf<Precision>(6)))/5));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(0)+(1+amp::sqrt<Precision>(amp::ampf<Precision>(6)))/5));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(2)-amp::ampf<Precision>(2)/amp::ampf<Precision>(9)));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(1)-(16+amp::sqrt<Precision>(amp::ampf<Precision>(6)))/18));
            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(0)-(16-amp::sqrt<Precision>(amp::ampf<Precision>(6)))/18));
            for(i=0; i<=1; i++)
            {
                recerrors = recerrors || x(i)>=x(i+1);
            }
        }
        else
        {
            recerrors = true;
        }
        recerrors = recerrors || err>errtol;
        
        //
        // test recurrence-based special cases (Legendre, Jacobi, Hermite, ...)
        // against another implementation (polynomial root-finder)
        //
        for(n=1; n<=20; n++)
        {
            
            //
            // test gauss-legendre
            //
            err = 0;
            gq::gqgenerategausslegendre<Precision>(n, info, x, w);
            if( info>0 )
            {
                buildgausslegendrequadrature<Precision>(n, x2, w2);
                for(i=0; i<=n-1; i++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(x(i)-x2(i)));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(w(i)-w2(i)));
                }
            }
            else
            {
                specerrors = true;
            }
            specerrors = specerrors || err>errtol;
            
            //
            // Test Gauss-Jacobi.
            // Since task is much more difficult we will use less strict
            // threshold.
            //
            err = 0;
            for(akind=0; akind<=9; akind++)
            {
                for(bkind=0; bkind<=9; bkind++)
                {
                    alphac = mapkind<Precision>(akind);
                    betac = mapkind<Precision>(bkind);
                    gq::gqgenerategaussjacobi<Precision>(n, alphac, betac, info, x, w);
                    if( info>0 )
                    {
                        buildgaussjacobiquadrature<Precision>(n, alphac, betac, x2, w2);
                        for(i=0; i<=n-1; i++)
                        {
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(x(i)-x2(i)));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(w(i)-w2(i)));
                        }
                    }
                    else
                    {
                        specerrors = true;
                    }
                }
            }
            specerrors = specerrors || err>nonstricterrtol;
            
            //
            // special test for Gauss-Jacobi (Chebyshev weight
            // function with analytically known nodes/weights)
            //
            err = 0;
            gq::gqgenerategaussjacobi<Precision>(n, -amp::ampf<Precision>("0.5"), -amp::ampf<Precision>("0.5"), info, x, w);
            if( info>0 )
            {
                for(i=0; i<=n-1; i++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(x(i)+amp::cos<Precision>(amp::pi<Precision>()*(i+amp::ampf<Precision>("0.5"))/n)));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(w(i)-amp::pi<Precision>()/n));
                }
            }
            else
            {
                specerrors = true;
            }
            specerrors = specerrors || err>stricterrtol;
            
            //
            // Test Gauss-Laguerre
            //
            err = 0;
            for(akind=0; akind<=9; akind++)
            {
                alphac = mapkind<Precision>(akind);
                gq::gqgenerategausslaguerre<Precision>(n, alphac, info, x, w);
                if( info>0 )
                {
                    buildgausslaguerrequadrature<Precision>(n, alphac, x2, w2);
                    for(i=0; i<=n-1; i++)
                    {
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(x(i)-x2(i)));
                        err = amp::maximum<Precision>(err, amp::abs<Precision>(w(i)-w2(i)));
                    }
                }
                else
                {
                    specerrors = true;
                }
            }
            specerrors = specerrors || err>nonstricterrtol;
            
            //
            // Test Gauss-Hermite
            //
            err = 0;
            gq::gqgenerategausshermite<Precision>(n, info, x, w);
            if( info>0 )
            {
                buildgausshermitequadrature<Precision>(n, x2, w2);
                for(i=0; i<=n-1; i++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(x(i)-x2(i)));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(w(i)-w2(i)));
                }
            }
            else
            {
                specerrors = true;
            }
            specerrors = specerrors || err>nonstricterrtol;
        }
        
        //
        // end
        //
        waserrors = recerrors || specerrors;
        if( !silent )
        {
            printf("TESTING GAUSS QUADRATURES\n");
            printf("FINAL RESULT:                             ");
            if( waserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* SPECIAL CASES (LEGENDRE/JACOBI/..)      ");
            if( specerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* RECURRENCE-BASED:                       ");
            if( recerrors )
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
            printf("\n\n");
        }
        result = !waserrors;
        return result;
    }


    /*************************************************************************
    Gauss-Hermite, another variant
    *************************************************************************/
    template<unsigned int Precision>
    void buildgausshermitequadrature(int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        int i;
        int j;
        amp::ampf<Precision> r;
        amp::ampf<Precision> r1;
        amp::ampf<Precision> p1;
        amp::ampf<Precision> p2;
        amp::ampf<Precision> p3;
        amp::ampf<Precision> dp3;
        amp::ampf<Precision> pipm4;
        amp::ampf<Precision> tmp;


        x.setbounds(0, n-1);
        w.setbounds(0, n-1);
        pipm4 = amp::pow<Precision>(amp::pi<Precision>(), -amp::ampf<Precision>("0.25"));
        for(i=0; i<=(n+1)/2-1; i++)
        {
            if( i==0 )
            {
                r = amp::sqrt<Precision>(amp::ampf<Precision>(2*n+1))-amp::ampf<Precision>("1.85575")*amp::pow<Precision>(amp::ampf<Precision>(2*n+1), -amp::ampf<Precision>(1)/amp::ampf<Precision>(6));
            }
            else
            {
                if( i==1 )
                {
                    r = r-amp::ampf<Precision>("1.14")*amp::pow<Precision>(amp::ampf<Precision>(n), amp::ampf<Precision>("0.426"))/r;
                }
                else
                {
                    if( i==2 )
                    {
                        r = amp::ampf<Precision>("1.86")*r-amp::ampf<Precision>("0.86")*x(0);
                    }
                    else
                    {
                        if( i==3 )
                        {
                            r = amp::ampf<Precision>("1.91")*r-amp::ampf<Precision>("0.91")*x(1);
                        }
                        else
                        {
                            r = 2*r-x(i-2);
                        }
                    }
                }
            }
            do
            {
                p2 = 0;
                p3 = pipm4;
                for(j=0; j<=n-1; j++)
                {
                    p1 = p2;
                    p2 = p3;
                    p3 = p2*r*amp::sqrt<Precision>(amp::ampf<Precision>(2)/(amp::ampf<Precision>(j+1)))-p1*amp::sqrt<Precision>(amp::ampf<Precision>(j)/(amp::ampf<Precision>(j+1)));
                }
                dp3 = amp::sqrt<Precision>(amp::ampf<Precision>(2*j))*p2;
                r1 = r;
                r = r-p3/dp3;
            }
            while( amp::abs<Precision>(r-r1)>=amp::ampf<Precision>::getAlgoPascalEpsilon()*(1+amp::abs<Precision>(r))*100 );
            x(i) = r;
            w(i) = 2/(dp3*dp3);
            x(n-1-i) = -x(i);
            w(n-1-i) = w(i);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-2-i; j++)
            {
                if( x(j)>=x(j+1) )
                {
                    tmp = x(j);
                    x(j) = x(j+1);
                    x(j+1) = tmp;
                    tmp = w(j);
                    w(j) = w(j+1);
                    w(j+1) = tmp;
                }
            }
        }
    }


    /*************************************************************************
    Maps:
        0   =>  -0.9
        1   =>  -0.5
        2   =>  -0.1
        3   =>   0.0
        4   =>  +0.1
        5   =>  +0.5
        6   =>  +0.9
        7   =>  +1.0
        8   =>  +1.5
        9   =>  +2.0
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> mapkind(int k)
    {
        amp::ampf<Precision> result;


        result = 0;
        if( k==0 )
        {
            result = -amp::ampf<Precision>("0.9");
        }
        if( k==1 )
        {
            result = -amp::ampf<Precision>("0.5");
        }
        if( k==2 )
        {
            result = -amp::ampf<Precision>("0.1");
        }
        if( k==3 )
        {
            result = amp::ampf<Precision>("0.0");
        }
        if( k==4 )
        {
            result = +amp::ampf<Precision>("0.1");
        }
        if( k==5 )
        {
            result = +amp::ampf<Precision>("0.5");
        }
        if( k==6 )
        {
            result = +amp::ampf<Precision>("0.9");
        }
        if( k==7 )
        {
            result = +amp::ampf<Precision>("1.0");
        }
        if( k==8 )
        {
            result = +amp::ampf<Precision>("1.5");
        }
        if( k==9 )
        {
            result = +amp::ampf<Precision>("2.0");
        }
        return result;
    }


    /*************************************************************************
    Gauss-Legendre, another variant
    *************************************************************************/
    template<unsigned int Precision>
    void buildgausslegendrequadrature(int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        int i;
        int j;
        amp::ampf<Precision> r;
        amp::ampf<Precision> r1;
        amp::ampf<Precision> p1;
        amp::ampf<Precision> p2;
        amp::ampf<Precision> p3;
        amp::ampf<Precision> dp3;
        amp::ampf<Precision> tmp;


        x.setbounds(0, n-1);
        w.setbounds(0, n-1);
        for(i=0; i<=(n+1)/2-1; i++)
        {
            r = amp::cos<Precision>(amp::pi<Precision>()*(4*i+3)/(4*n+2));
            do
            {
                p2 = 0;
                p3 = 1;
                for(j=0; j<=n-1; j++)
                {
                    p1 = p2;
                    p2 = p3;
                    p3 = ((2*j+1)*r*p2-j*p1)/(j+1);
                }
                dp3 = n*(r*p3-p2)/(r*r-1);
                r1 = r;
                r = r-p3/dp3;
            }
            while( amp::abs<Precision>(r-r1)>=amp::ampf<Precision>::getAlgoPascalEpsilon()*(1+amp::abs<Precision>(r))*100 );
            x(i) = r;
            x(n-1-i) = -r;
            w(i) = 2/((1-r*r)*dp3*dp3);
            w(n-1-i) = 2/((1-r*r)*dp3*dp3);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-2-i; j++)
            {
                if( x(j)>=x(j+1) )
                {
                    tmp = x(j);
                    x(j) = x(j+1);
                    x(j+1) = tmp;
                    tmp = w(j);
                    w(j) = w(j+1);
                    w(j+1) = tmp;
                }
            }
        }
    }


    /*************************************************************************
    Gauss-Jacobi, another variant
    *************************************************************************/
    template<unsigned int Precision>
    void buildgaussjacobiquadrature(int n,
        amp::ampf<Precision> alpha,
        amp::ampf<Precision> beta,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        int i;
        int j;
        amp::ampf<Precision> r;
        amp::ampf<Precision> r1;
        amp::ampf<Precision> t1;
        amp::ampf<Precision> t2;
        amp::ampf<Precision> t3;
        amp::ampf<Precision> p1;
        amp::ampf<Precision> p2;
        amp::ampf<Precision> p3;
        amp::ampf<Precision> pp;
        amp::ampf<Precision> an;
        amp::ampf<Precision> bn;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        amp::ampf<Precision> c;
        amp::ampf<Precision> tmpsgn;
        amp::ampf<Precision> tmp;
        amp::ampf<Precision> alfbet;
        amp::ampf<Precision> temp;
        int its;


        x.setbounds(0, n-1);
        w.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            if( i==0 )
            {
                an = alpha/n;
                bn = beta/n;
                t1 = (1+alpha)*(amp::ampf<Precision>("2.78")/(4+n*n)+amp::ampf<Precision>("0.768")*an/n);
                t2 = 1+amp::ampf<Precision>("1.48")*an+amp::ampf<Precision>("0.96")*bn+amp::ampf<Precision>("0.452")*an*an+amp::ampf<Precision>("0.83")*an*bn;
                r = (t2-t1)/t2;
            }
            else
            {
                if( i==1 )
                {
                    t1 = (amp::ampf<Precision>("4.1")+alpha)/((1+alpha)*(1+amp::ampf<Precision>("0.156")*alpha));
                    t2 = 1+amp::ampf<Precision>("0.06")*(n-8)*(1+amp::ampf<Precision>("0.12")*alpha)/n;
                    t3 = 1+amp::ampf<Precision>("0.012")*beta*(1+amp::ampf<Precision>("0.25")*amp::abs<Precision>(alpha))/n;
                    r = r-t1*t2*t3*(1-r);
                }
                else
                {
                    if( i==2 )
                    {
                        t1 = (amp::ampf<Precision>("1.67")+amp::ampf<Precision>("0.28")*alpha)/(1+amp::ampf<Precision>("0.37")*alpha);
                        t2 = 1+amp::ampf<Precision>("0.22")*(n-8)/n;
                        t3 = 1+8*beta/((amp::ampf<Precision>("6.28")+beta)*n*n);
                        r = r-t1*t2*t3*(x(0)-r);
                    }
                    else
                    {
                        if( i<n-2 )
                        {
                            r = 3*x(i-1)-3*x(i-2)+x(i-3);
                        }
                        else
                        {
                            if( i==n-2 )
                            {
                                t1 = (1+amp::ampf<Precision>("0.235")*beta)/(amp::ampf<Precision>("0.766")+amp::ampf<Precision>("0.119")*beta);
                                t2 = 1/(1+amp::ampf<Precision>("0.639")*(n-4)/(1+amp::ampf<Precision>("0.71")*(n-4)));
                                t3 = 1/(1+20*alpha/((amp::ampf<Precision>("7.5")+alpha)*n*n));
                                r = r+t1*t2*t3*(r-x(i-2));
                            }
                            else
                            {
                                if( i==n-1 )
                                {
                                    t1 = (1+amp::ampf<Precision>("0.37")*beta)/(amp::ampf<Precision>("1.67")+amp::ampf<Precision>("0.28")*beta);
                                    t2 = 1/(1+amp::ampf<Precision>("0.22")*(n-8)/n);
                                    t3 = 1/(1+8*alpha/((amp::ampf<Precision>("6.28")+alpha)*n*n));
                                    r = r+t1*t2*t3*(r-x(i-2));
                                }
                            }
                        }
                    }
                }
            }
            alfbet = alpha+beta;
            do
            {
                temp = 2+alfbet;
                p1 = (alpha-beta+temp*r)*amp::ampf<Precision>("0.5");
                p2 = 1;
                for(j=2; j<=n; j++)
                {
                    p3 = p2;
                    p2 = p1;
                    temp = 2*j+alfbet;
                    a = 2*j*(j+alfbet)*(temp-2);
                    b = (temp-1)*(alpha*alpha-beta*beta+temp*(temp-2)*r);
                    c = 2*(j-1+alpha)*(j-1+beta)*temp;
                    p1 = (b*p2-c*p3)/a;
                }
                pp = (n*(alpha-beta-temp*r)*p1+2*(n+alpha)*(n+beta)*p2)/(temp*(1-r*r));
                r1 = r;
                r = r1-p1/pp;
            }
            while( amp::abs<Precision>(r-r1)>=amp::ampf<Precision>::getAlgoPascalEpsilon()*(1+amp::abs<Precision>(r))*100 );
            x(i) = r;
            w(i) = amp::exp<Precision>(gammafunc::lngamma<Precision>(alpha+n, tmpsgn)+gammafunc::lngamma<Precision>(beta+n, tmpsgn)-gammafunc::lngamma<Precision>(amp::ampf<Precision>(n+1), tmpsgn)-gammafunc::lngamma<Precision>(n+alfbet+1, tmpsgn))*temp*amp::pow<Precision>(amp::ampf<Precision>(2), alfbet)/(pp*p2);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-2-i; j++)
            {
                if( x(j)>=x(j+1) )
                {
                    tmp = x(j);
                    x(j) = x(j+1);
                    x(j+1) = tmp;
                    tmp = w(j);
                    w(j) = w(j+1);
                    w(j+1) = tmp;
                }
            }
        }
    }


    /*************************************************************************
    Gauss-Laguerre, another variant
    *************************************************************************/
    template<unsigned int Precision>
    void buildgausslaguerrequadrature(int n,
        amp::ampf<Precision> alpha,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        int i;
        int j;
        amp::ampf<Precision> r;
        amp::ampf<Precision> r1;
        amp::ampf<Precision> p1;
        amp::ampf<Precision> p2;
        amp::ampf<Precision> p3;
        amp::ampf<Precision> dp3;
        amp::ampf<Precision> tsg;
        amp::ampf<Precision> tmp;


        x.setbounds(0, n-1);
        w.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            if( i==0 )
            {
                r = (1+alpha)*(3+amp::ampf<Precision>("0.92")*alpha)/(1+amp::ampf<Precision>("2.4")*n+amp::ampf<Precision>("1.8")*alpha);
            }
            else
            {
                if( i==1 )
                {
                    r = r+(15+amp::ampf<Precision>("6.25")*alpha)/(1+amp::ampf<Precision>("0.9")*alpha+amp::ampf<Precision>("2.5")*n);
                }
                else
                {
                    r = r+((1+amp::ampf<Precision>("2.55")*(i-1))/(amp::ampf<Precision>("1.9")*(i-1))+amp::ampf<Precision>("1.26")*(i-1)*alpha/(1+amp::ampf<Precision>("3.5")*(i-1)))/(1+amp::ampf<Precision>("0.3")*alpha)*(r-x(i-2));
                }
            }
            do
            {
                p2 = 0;
                p3 = 1;
                for(j=0; j<=n-1; j++)
                {
                    p1 = p2;
                    p2 = p3;
                    p3 = ((-r+2*j+alpha+1)*p2-(j+alpha)*p1)/(j+1);
                }
                dp3 = (n*p3-(n+alpha)*p2)/r;
                r1 = r;
                r = r-p3/dp3;
            }
            while( amp::abs<Precision>(r-r1)>=amp::ampf<Precision>::getAlgoPascalEpsilon()*(1+amp::abs<Precision>(r))*100 );
            x(i) = r;
            w(i) = -amp::exp<Precision>(gammafunc::lngamma<Precision>(alpha+n, tsg)-gammafunc::lngamma<Precision>(amp::ampf<Precision>(n), tsg))/(dp3*n*p2);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-2-i; j++)
            {
                if( x(j)>=x(j+1) )
                {
                    tmp = x(j);
                    x(j) = x(j+1);
                    x(j+1) = tmp;
                    tmp = w(j);
                    w(j) = w(j+1);
                    w(j+1) = tmp;
                }
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testgq_test_silent()
    {
        bool result;


        result = testgqunit<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testgq_test()
    {
        bool result;


        result = testgqunit<Precision>(false);
        return result;
    }
} // namespace

#endif
