
#ifndef _testsblasunit_h
#define _testsblasunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "sblas.h"
namespace testsblasunit
{
    template<unsigned int Precision>
    bool testsblas(bool silent);
    template<unsigned int Precision>
    bool testsblasunit_test_silent();
    template<unsigned int Precision>
    bool testsblasunit_test();


    template<unsigned int Precision>
    bool testsblas(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > ua;
        ap::template_2d_array< amp::ampf<Precision> > la;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y1;
        ap::template_1d_array< amp::ampf<Precision> > y2;
        ap::template_1d_array< amp::ampf<Precision> > y3;
        int n;
        int maxn;
        int i;
        int j;
        int i1;
        int i2;
        int gpass;
        bool waserrors;
        amp::ampf<Precision> mverr;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> alpha;
        amp::ampf<Precision> v;


        mverr = 0;
        waserrors = false;
        maxn = 10;
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Test MV
        //
        for(n=2; n<=maxn; n++)
        {
            a.setbounds(1, n, 1, n);
            ua.setbounds(1, n, 1, n);
            la.setbounds(1, n, 1, n);
            x.setbounds(1, n);
            y1.setbounds(1, n);
            y2.setbounds(1, n);
            y3.setbounds(1, n);
            
            //
            // fill A, UA, LA
            //
            for(i=1; i<=n; i++)
            {
                a(i,i) = 2*amp::ampf<Precision>::getRandom()-1;
                for(j=i+1; j<=n; j++)
                {
                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    a(j,i) = a(i,j);
                }
            }
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=n; j++)
                {
                    ua(i,j) = 0;
                }
            }
            for(i=1; i<=n; i++)
            {
                for(j=i; j<=n; j++)
                {
                    ua(i,j) = a(i,j);
                }
            }
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=n; j++)
                {
                    la(i,j) = 0;
                }
            }
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=i; j++)
                {
                    la(i,j) = a(i,j);
                }
            }
            
            //
            // test on different I1, I2
            //
            for(i1=1; i1<=n; i1++)
            {
                for(i2=i1; i2<=n; i2++)
                {
                    
                    //
                    // Fill X, choose Alpha
                    //
                    for(i=1; i<=i2-i1+1; i++)
                    {
                        x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    alpha = 2*amp::ampf<Precision>::getRandom()-1;
                    
                    //
                    // calculate A*x, UA*x, LA*x
                    //
                    for(i=i1; i<=i2; i++)
                    {
                        v = amp::vdotproduct(a.getrow(i, i1, i2), x.getvector(1, i2-i1+1));
                        y1(i-i1+1) = alpha*v;
                    }
                    sblas::symmetricmatrixvectormultiply<Precision>(ua, true, i1, i2, x, alpha, y2);
                    sblas::symmetricmatrixvectormultiply<Precision>(la, false, i1, i2, x, alpha, y3);
                    
                    //
                    // Calculate error
                    //
                    amp::vsub(y2.getvector(1, i2-i1+1), y1.getvector(1, i2-i1+1));
                    v = amp::vdotproduct(y2.getvector(1, i2-i1+1), y2.getvector(1, i2-i1+1));
                    mverr = amp::maximum<Precision>(mverr, amp::sqrt<Precision>(v));
                    amp::vsub(y3.getvector(1, i2-i1+1), y1.getvector(1, i2-i1+1));
                    v = amp::vdotproduct(y3.getvector(1, i2-i1+1), y3.getvector(1, i2-i1+1));
                    mverr = amp::maximum<Precision>(mverr, amp::sqrt<Precision>(v));
                }
            }
        }
        
        //
        // report
        //
        waserrors = mverr>threshold;
        if( !silent )
        {
            printf("TESTING SYMMETRIC BLAS\n");
            printf("MV error:                                %5.3le\n",
                double(amp::ampf<Precision>(mverr).toDouble()));
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


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testsblasunit_test_silent()
    {
        bool result;


        result = testsblas<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testsblasunit_test()
    {
        bool result;


        result = testsblas<Precision>(false);
        return result;
    }
} // namespace

#endif
