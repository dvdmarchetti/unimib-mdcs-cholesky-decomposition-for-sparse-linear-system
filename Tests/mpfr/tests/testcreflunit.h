
#ifndef _testcreflunit_h
#define _testcreflunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "creflections.h"
namespace testcreflunit
{
    template<unsigned int Precision>
    bool testcrefl(bool silent);
    template<unsigned int Precision>
    bool testcreflunit_test_silent();
    template<unsigned int Precision>
    bool testcreflunit_test();


    template<unsigned int Precision>
    bool testcrefl(bool silent)
    {
        bool result;
        int i;
        int j;
        int k;
        int n;
        int m;
        int maxmn;
        ap::template_1d_array< amp::campf<Precision> > x;
        ap::template_1d_array< amp::campf<Precision> > v;
        ap::template_1d_array< amp::campf<Precision> > work;
        ap::template_2d_array< amp::campf<Precision> > h;
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > b;
        ap::template_2d_array< amp::campf<Precision> > c;
        amp::campf<Precision> tmp;
        amp::campf<Precision> beta;
        amp::campf<Precision> tau;
        amp::ampf<Precision> err;
        amp::ampf<Precision> mer;
        amp::ampf<Precision> mel;
        amp::ampf<Precision> meg;
        int pass;
        int passcount;
        bool waserrors;
        amp::ampf<Precision> threshold;
        int i_;


        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        passcount = 1000;
        mer = 0;
        mel = 0;
        meg = 0;
        for(pass=1; pass<=passcount; pass++)
        {
            
            //
            // Task
            //
            n = 1+ap::randominteger(10);
            m = 1+ap::randominteger(10);
            maxmn = ap::maxint(m, n);
            
            //
            // Initialize
            //
            x.setbounds(1, maxmn);
            v.setbounds(1, maxmn);
            work.setbounds(1, maxmn);
            h.setbounds(1, maxmn, 1, maxmn);
            a.setbounds(1, maxmn, 1, maxmn);
            b.setbounds(1, maxmn, 1, maxmn);
            c.setbounds(1, maxmn, 1, maxmn);
            
            //
            // GenerateReflection
            //
            for(i=1; i<=n; i++)
            {
                x(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                x(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                v(i) = x(i);
            }
            creflections::complexgeneratereflection<Precision>(v, n, tau);
            beta = v(1);
            v(1) = 1;
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=n; j++)
                {
                    if( i==j )
                    {
                        h(i,j) = 1-tau*v(i)*amp::conj<Precision>(v(j));
                    }
                    else
                    {
                        h(i,j) = -tau*v(i)*amp::conj<Precision>(v(j));
                    }
                }
            }
            err = 0;
            for(i=1; i<=n; i++)
            {
                tmp = 0.0;
                for(i_=1; i_<=n;i_++)
                {
                    tmp += amp::conj(h(i_,i))*x(i_);
                }
                if( i==1 )
                {
                    err = amp::maximum<Precision>(err, amp::abscomplex<Precision>(tmp-beta));
                }
                else
                {
                    err = amp::maximum<Precision>(err, amp::abscomplex<Precision>(tmp));
                }
            }
            err = amp::maximum<Precision>(err, amp::abs<Precision>(beta.y));
            meg = amp::maximum<Precision>(meg, err);
            
            //
            // ApplyReflectionFromTheLeft
            //
            for(i=1; i<=m; i++)
            {
                x(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                x(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                v(i) = x(i);
            }
            for(i=1; i<=m; i++)
            {
                for(j=1; j<=n; j++)
                {
                    a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    b(i,j) = a(i,j);
                }
            }
            creflections::complexgeneratereflection<Precision>(v, m, tau);
            beta = v(1);
            v(1) = 1;
            creflections::complexapplyreflectionfromtheleft<Precision>(b, tau, v, 1, m, 1, n, work);
            for(i=1; i<=m; i++)
            {
                for(j=1; j<=m; j++)
                {
                    if( i==j )
                    {
                        h(i,j) = 1-tau*v(i)*amp::conj<Precision>(v(j));
                    }
                    else
                    {
                        h(i,j) = -tau*v(i)*amp::conj<Precision>(v(j));
                    }
                }
            }
            for(i=1; i<=m; i++)
            {
                for(j=1; j<=n; j++)
                {
                    tmp = 0.0;
                    for(i_=1; i_<=m;i_++)
                    {
                        tmp += h(i,i_)*a(i_,j);
                    }
                    c(i,j) = tmp;
                }
            }
            err = 0;
            for(i=1; i<=m; i++)
            {
                for(j=1; j<=n; j++)
                {
                    err = amp::maximum<Precision>(err, amp::abscomplex<Precision>(b(i,j)-c(i,j)));
                }
            }
            mel = amp::maximum<Precision>(mel, err);
            
            //
            // ApplyReflectionFromTheRight
            //
            for(i=1; i<=n; i++)
            {
                x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                v(i) = x(i);
            }
            for(i=1; i<=m; i++)
            {
                for(j=1; j<=n; j++)
                {
                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    b(i,j) = a(i,j);
                }
            }
            creflections::complexgeneratereflection<Precision>(v, n, tau);
            beta = v(1);
            v(1) = 1;
            creflections::complexapplyreflectionfromtheright<Precision>(b, tau, v, 1, m, 1, n, work);
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=n; j++)
                {
                    if( i==j )
                    {
                        h(i,j) = 1-tau*v(i)*amp::conj<Precision>(v(j));
                    }
                    else
                    {
                        h(i,j) = -tau*v(i)*amp::conj<Precision>(v(j));
                    }
                }
            }
            for(i=1; i<=m; i++)
            {
                for(j=1; j<=n; j++)
                {
                    tmp = 0.0;
                    for(i_=1; i_<=n;i_++)
                    {
                        tmp += a(i,i_)*h(i_,j);
                    }
                    c(i,j) = tmp;
                }
            }
            err = 0;
            for(i=1; i<=m; i++)
            {
                for(j=1; j<=n; j++)
                {
                    err = amp::maximum<Precision>(err, amp::abscomplex<Precision>(b(i,j)-c(i,j)));
                }
            }
            mer = amp::maximum<Precision>(mer, err);
        }
        
        //
        // Overflow crash test
        //
        x.setbounds(1, 10);
        v.setbounds(1, 10);
        for(i=1; i<=10; i++)
        {
            v(i) = amp::ampf<Precision>::getAlgoPascalMaxNumber()*amp::ampf<Precision>("0.01")*(2*amp::ampf<Precision>::getRandom()-1);
        }
        creflections::complexgeneratereflection<Precision>(v, 10, tau);
        
        //
        // report
        //
        waserrors = meg>threshold || mel>threshold || mer>threshold;
        if( !silent )
        {
            printf("TESTING COMPLEX REFLECTIONS\n");
            printf("Generate error:                          %5.3le\n",
                double(amp::ampf<Precision>(meg).toDouble()));
            printf("Apply(L) error:                          %5.3le\n",
                double(amp::ampf<Precision>(mel).toDouble()));
            printf("Apply(R) error:                          %5.3le\n",
                double(amp::ampf<Precision>(mer).toDouble()));
            printf("Threshold:                               %5.3le\n",
                double(amp::ampf<Precision>(threshold).toDouble()));
            printf("Overflow crash test:                     PASSED\n");
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
    bool testcreflunit_test_silent()
    {
        bool result;


        result = testcrefl<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testcreflunit_test()
    {
        bool result;


        result = testcrefl<Precision>(false);
        return result;
    }
} // namespace

#endif
