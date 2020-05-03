
#ifndef _testhblasunit_h
#define _testhblasunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "hblas.h"
namespace testhblasunit
{
    template<unsigned int Precision>
    bool testhblas(bool silent);
    template<unsigned int Precision>
    bool testhblasunit_test_silent();
    template<unsigned int Precision>
    bool testhblasunit_test();


    template<unsigned int Precision>
    bool testhblas(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > ua;
        ap::template_2d_array< amp::campf<Precision> > la;
        ap::template_1d_array< amp::campf<Precision> > x;
        ap::template_1d_array< amp::campf<Precision> > y1;
        ap::template_1d_array< amp::campf<Precision> > y2;
        ap::template_1d_array< amp::campf<Precision> > y3;
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
        amp::campf<Precision> alpha;
        amp::campf<Precision> v;
        int i_;
        int i1_;


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
                a(i,i).x = 2*amp::ampf<Precision>::getRandom()-1;
                a(i,i).y = 0;
                for(j=i+1; j<=n; j++)
                {
                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    a(j,i) = amp::conj<Precision>(a(i,j));
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
                        x(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                        x(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    alpha.x = 2*amp::ampf<Precision>::getRandom()-1;
                    alpha.y = 2*amp::ampf<Precision>::getRandom()-1;
                    
                    //
                    // calculate A*x, UA*x, LA*x
                    //
                    for(i=i1; i<=i2; i++)
                    {
                        i1_ = (1)-(i1);
                        v = 0.0;
                        for(i_=i1; i_<=i2;i_++)
                        {
                            v += a(i,i_)*x(i_+i1_);
                        }
                        y1(i-i1+1) = alpha*v;
                    }
                    hblas::hermitianmatrixvectormultiply<Precision>(ua, true, i1, i2, x, alpha, y2);
                    hblas::hermitianmatrixvectormultiply<Precision>(la, false, i1, i2, x, alpha, y3);
                    
                    //
                    // Calculate error
                    //
                    for(i_=1; i_<=i2-i1+1;i_++)
                    {
                        y2(i_) = y2(i_) - y1(i_);
                    }
                    v = 0.0;
                    for(i_=1; i_<=i2-i1+1;i_++)
                    {
                        v += y2(i_)*amp::conj(y2(i_));
                    }
                    mverr = amp::maximum<Precision>(mverr, amp::sqrt<Precision>(amp::abscomplex<Precision>(v)));
                    for(i_=1; i_<=i2-i1+1;i_++)
                    {
                        y3(i_) = y3(i_) - y1(i_);
                    }
                    v = 0.0;
                    for(i_=1; i_<=i2-i1+1;i_++)
                    {
                        v += y3(i_)*amp::conj(y3(i_));
                    }
                    mverr = amp::maximum<Precision>(mverr, amp::sqrt<Precision>(amp::abscomplex<Precision>(v)));
                }
            }
        }
        
        //
        // report
        //
        waserrors = mverr>threshold;
        if( !silent )
        {
            printf("TESTING HERMITIAN BLAS\n");
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
    bool testhblasunit_test_silent()
    {
        bool result;


        result = testhblas<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testhblasunit_test()
    {
        bool result;


        result = testhblas<Precision>(false);
        return result;
    }
} // namespace

#endif
