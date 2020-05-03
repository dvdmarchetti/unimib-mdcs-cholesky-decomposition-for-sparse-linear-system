
#ifndef _testspdgevdunit_h
#define _testspdgevdunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "sblas.h"
#include "blas.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "hblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "spdgevd.h"
namespace testspdgevdunit
{
    template<unsigned int Precision>
    bool testspdgevd(bool silent);
    template<unsigned int Precision>
    bool testspdgevdunit_test_silent();
    template<unsigned int Precision>
    bool testspdgevdunit_test();


    /*************************************************************************
    Testing bidiagonal SVD decomposition subroutine
    *************************************************************************/
    template<unsigned int Precision>
    bool testspdgevd(bool silent)
    {
        bool result;
        int pass;
        int n;
        int passcount;
        int maxn;
        int atask;
        int btask;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_2d_array< amp::ampf<Precision> > afull;
        ap::template_2d_array< amp::ampf<Precision> > bfull;
        ap::template_2d_array< amp::ampf<Precision> > l;
        ap::template_2d_array< amp::ampf<Precision> > z;
        bool isuppera;
        bool isupperb;
        int i;
        int j;
        int minij;
        amp::ampf<Precision> v;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        bool cw;
        amp::ampf<Precision> err;
        amp::ampf<Precision> valerr;
        amp::ampf<Precision> threshold;
        bool waserrors;
        bool wfailed;
        bool wnsorted;


        threshold = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        valerr = 0;
        wfailed = false;
        wnsorted = false;
        maxn = 20;
        passcount = 5;
        
        //
        // Main cycle
        //
        for(n=1; n<=maxn; n++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                for(atask=0; atask<=1; atask++)
                {
                    for(btask=0; btask<=1; btask++)
                    {
                        isuppera = atask==0;
                        isupperb = btask==0;
                        
                        //
                        // Initialize A, B, AFull, BFull
                        //
                        t1.setbounds(0, n-1);
                        a.setbounds(0, n-1, 0, n-1);
                        b.setbounds(0, n-1, 0, n-1);
                        afull.setbounds(0, n-1, 0, n-1);
                        bfull.setbounds(0, n-1, 0, n-1);
                        l.setbounds(0, n-1, 0, n-1);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                a(j,i) = a(i,j);
                                afull(i,j) = a(i,j);
                                afull(j,i) = a(i,j);
                            }
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=i+1; j<=n-1; j++)
                            {
                                l(i,j) = amp::ampf<Precision>::getRandom();
                                l(j,i) = l(i,j);
                            }
                            l(i,i) = amp::ampf<Precision>("1.5")+amp::ampf<Precision>::getRandom();
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                minij = ap::minint(i, j);
                                v = amp::vdotproduct(l.getrow(i, 0, minij), l.getcolumn(j, 0, minij));
                                b(i,j) = v;
                                b(j,i) = v;
                                bfull(i,j) = v;
                                bfull(j,i) = v;
                            }
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                if( isuppera )
                                {
                                    if( j<i )
                                    {
                                        a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                    }
                                }
                                else
                                {
                                    if( i<j )
                                    {
                                        a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                    }
                                }
                                if( isupperb )
                                {
                                    if( j<i )
                                    {
                                        b(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                    }
                                }
                                else
                                {
                                    if( i<j )
                                    {
                                        b(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                                    }
                                }
                            }
                        }
                        
                        //
                        // Problem 1
                        //
                        if( !spdgevd::smatrixgevd<Precision>(a, n, isuppera, b, isupperb, 1, 1, d, z) )
                        {
                            wfailed = true;
                            continue;
                        }
                        err = 0;
                        for(j=0; j<=n-1; j++)
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                v1 = amp::vdotproduct(afull.getrow(i, 0, n-1), z.getcolumn(j, 0, n-1));
                                v2 = amp::vdotproduct(bfull.getrow(i, 0, n-1), z.getcolumn(j, 0, n-1));
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(v1-d(j)*v2));
                            }
                        }
                        valerr = amp::maximum<Precision>(err, valerr);
                        
                        //
                        // Problem 2
                        //
                        if( !spdgevd::smatrixgevd<Precision>(a, n, isuppera, b, isupperb, 1, 2, d, z) )
                        {
                            wfailed = true;
                            continue;
                        }
                        err = 0;
                        for(j=0; j<=n-1; j++)
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                v1 = amp::vdotproduct(bfull.getrow(i, 0, n-1), z.getcolumn(j, 0, n-1));
                                t1(i) = v1;
                            }
                            for(i=0; i<=n-1; i++)
                            {
                                v2 = amp::vdotproduct(afull.getrow(i, 0, n-1), t1.getvector(0, n-1));
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(v2-d(j)*z(i,j)));
                            }
                        }
                        valerr = amp::maximum<Precision>(err, valerr);
                        
                        //
                        // Test problem 3
                        //
                        if( !spdgevd::smatrixgevd<Precision>(a, n, isuppera, b, isupperb, 1, 3, d, z) )
                        {
                            wfailed = true;
                            continue;
                        }
                        err = 0;
                        for(j=0; j<=n-1; j++)
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                v1 = amp::vdotproduct(afull.getrow(i, 0, n-1), z.getcolumn(j, 0, n-1));
                                t1(i) = v1;
                            }
                            for(i=0; i<=n-1; i++)
                            {
                                v2 = amp::vdotproduct(bfull.getrow(i, 0, n-1), t1.getvector(0, n-1));
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(v2-d(j)*z(i,j)));
                            }
                        }
                        valerr = amp::maximum<Precision>(err, valerr);
                    }
                }
            }
        }
        
        //
        // report
        //
        waserrors = valerr>threshold || wfailed || wnsorted;
        if( !silent )
        {
            printf("TESTING SYMMETRIC GEVD\n");
            printf("Av-lambdav error (generalized):          %5.3le\n",
                double(amp::ampf<Precision>(valerr).toDouble()));
            printf("Eigen values order:                      ");
            if( !wnsorted )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("Always converged:                        ");
            if( !wfailed )
            {
                printf("YES\n");
            }
            else
            {
                printf("NO\n");
            }
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
    bool testspdgevdunit_test_silent()
    {
        bool result;


        result = testspdgevd<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testspdgevdunit_test()
    {
        bool result;


        result = testspdgevd<Precision>(false);
        return result;
    }
} // namespace

#endif
