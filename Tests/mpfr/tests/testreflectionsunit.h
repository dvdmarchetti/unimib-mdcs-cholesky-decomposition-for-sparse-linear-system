
#ifndef _testreflectionsunit_h
#define _testreflectionsunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "reflections.h"
namespace testreflectionsunit
{
    template<unsigned int Precision>
    bool testreflections(bool silent);
    template<unsigned int Precision>
    bool testreflectionsunit_test_silent();
    template<unsigned int Precision>
    bool testreflectionsunit_test();


    template<unsigned int Precision>
    bool testreflections(bool silent)
    {
        bool result;
        int i;
        int j;
        int n;
        int m;
        int maxmn;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > v;
        ap::template_1d_array< amp::ampf<Precision> > work;
        ap::template_2d_array< amp::ampf<Precision> > h;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_2d_array< amp::ampf<Precision> > c;
        amp::ampf<Precision> tmp;
        amp::ampf<Precision> beta;
        amp::ampf<Precision> tau;
        amp::ampf<Precision> err;
        amp::ampf<Precision> mer;
        amp::ampf<Precision> mel;
        amp::ampf<Precision> meg;
        int pass;
        int passcount;
        amp::ampf<Precision> threshold;
        int tasktype;
        amp::ampf<Precision> xscale;


        passcount = 10;
        threshold = 100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        mer = 0;
        mel = 0;
        meg = 0;
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=10; n++)
            {
                for(m=1; m<=10; m++)
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
                    // GenerateReflection, three tasks are possible:
                    // * random X
                    // * zero X
                    // * non-zero X[1], all other are zeros
                    // * random X, near underflow scale
                    // * random X, near overflow scale
                    //
                    for(tasktype=0; tasktype<=4; tasktype++)
                    {
                        xscale = 1;
                        if( tasktype==0 )
                        {
                            for(i=1; i<=n; i++)
                            {
                                x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        if( tasktype==1 )
                        {
                            for(i=1; i<=n; i++)
                            {
                                x(i) = 0;
                            }
                        }
                        if( tasktype==2 )
                        {
                            x(1) = 2*amp::ampf<Precision>::getRandom()-1;
                            for(i=2; i<=n; i++)
                            {
                                x(i) = 0;
                            }
                        }
                        if( tasktype==3 )
                        {
                            for(i=1; i<=n; i++)
                            {
                                x(i) = (ap::randominteger(21)-10)*amp::ampf<Precision>::getAlgoPascalMinNumber();
                            }
                            xscale = 10*amp::ampf<Precision>::getAlgoPascalMinNumber();
                        }
                        if( tasktype==4 )
                        {
                            for(i=1; i<=n; i++)
                            {
                                x(i) = (2*amp::ampf<Precision>::getRandom()-1)*amp::ampf<Precision>::getAlgoPascalMaxNumber();
                            }
                            xscale = amp::ampf<Precision>::getAlgoPascalMaxNumber();
                        }
                        amp::vmove(v.getvector(1, n), x.getvector(1, n));
                        reflections::generatereflection<Precision>(v, n, tau);
                        beta = v(1);
                        v(1) = 1;
                        for(i=1; i<=n; i++)
                        {
                            for(j=1; j<=n; j++)
                            {
                                if( i==j )
                                {
                                    h(i,j) = 1-tau*v(i)*v(j);
                                }
                                else
                                {
                                    h(i,j) = -tau*v(i)*v(j);
                                }
                            }
                        }
                        err = 0;
                        for(i=1; i<=n; i++)
                        {
                            tmp = amp::vdotproduct(h.getrow(i, 1, n), x.getvector(1, n));
                            if( i==1 )
                            {
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(tmp-beta));
                            }
                            else
                            {
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(tmp));
                            }
                        }
                        meg = amp::maximum<Precision>(meg, err/xscale);
                    }
                    
                    //
                    // ApplyReflectionFromTheLeft
                    //
                    for(i=1; i<=m; i++)
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
                    reflections::generatereflection<Precision>(v, m, tau);
                    beta = v(1);
                    v(1) = 1;
                    reflections::applyreflectionfromtheleft<Precision>(b, tau, v, 1, m, 1, n, work);
                    for(i=1; i<=m; i++)
                    {
                        for(j=1; j<=m; j++)
                        {
                            if( i==j )
                            {
                                h(i,j) = 1-tau*v(i)*v(j);
                            }
                            else
                            {
                                h(i,j) = -tau*v(i)*v(j);
                            }
                        }
                    }
                    for(i=1; i<=m; i++)
                    {
                        for(j=1; j<=n; j++)
                        {
                            tmp = amp::vdotproduct(h.getrow(i, 1, m), a.getcolumn(j, 1, m));
                            c(i,j) = tmp;
                        }
                    }
                    err = 0;
                    for(i=1; i<=m; i++)
                    {
                        for(j=1; j<=n; j++)
                        {
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(b(i,j)-c(i,j)));
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
                    reflections::generatereflection<Precision>(v, n, tau);
                    beta = v(1);
                    v(1) = 1;
                    reflections::applyreflectionfromtheright<Precision>(b, tau, v, 1, m, 1, n, work);
                    for(i=1; i<=n; i++)
                    {
                        for(j=1; j<=n; j++)
                        {
                            if( i==j )
                            {
                                h(i,j) = 1-tau*v(i)*v(j);
                            }
                            else
                            {
                                h(i,j) = -tau*v(i)*v(j);
                            }
                        }
                    }
                    for(i=1; i<=m; i++)
                    {
                        for(j=1; j<=n; j++)
                        {
                            tmp = amp::vdotproduct(a.getrow(i, 1, n), h.getcolumn(j, 1, n));
                            c(i,j) = tmp;
                        }
                    }
                    err = 0;
                    for(i=1; i<=m; i++)
                    {
                        for(j=1; j<=n; j++)
                        {
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(b(i,j)-c(i,j)));
                        }
                    }
                    mer = amp::maximum<Precision>(mer, err);
                }
            }
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
        reflections::generatereflection<Precision>(v, 10, tau);
        result = meg<=threshold && mel<=threshold && mer<=threshold;
        if( !silent )
        {
            printf("TESTING REFLECTIONS\n");
            printf("Pass count is %0ld\n",
                long(passcount));
            printf("Generate     absolute error is       %5.3le\n",
                double(amp::ampf<Precision>(meg).toDouble()));
            printf("Apply(Left)  absolute error is       %5.3le\n",
                double(amp::ampf<Precision>(mel).toDouble()));
            printf("Apply(Right) absolute error is       %5.3le\n",
                double(amp::ampf<Precision>(mer).toDouble()));
            printf("Overflow crash test passed\n");
            if( result )
            {
                printf("TEST PASSED\n");
            }
            else
            {
                printf("TEST FAILED\n");
            }
        }
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testreflectionsunit_test_silent()
    {
        bool result;


        result = testreflections<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testreflectionsunit_test()
    {
        bool result;


        result = testreflections<Precision>(false);
        return result;
    }
} // namespace

#endif
