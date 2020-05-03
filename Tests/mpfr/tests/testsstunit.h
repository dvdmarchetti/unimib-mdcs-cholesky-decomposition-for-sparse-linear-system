
#ifndef _testsstunit_h
#define _testsstunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "trlinsolve.h"
namespace testsstunit
{
    template<unsigned int Precision>
    bool testsst(bool silent);
    template<unsigned int Precision>
    void makeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& b);
    template<unsigned int Precision>
    bool testsstunit_test_silent();
    template<unsigned int Precision>
    bool testsstunit_test();


    /*************************************************************************
    Main unittest subroutine
    *************************************************************************/
    template<unsigned int Precision>
    bool testsst(bool silent)
    {
        bool result;
        int maxmn;
        int passcount;
        amp::ampf<Precision> threshold;
        ap::template_2d_array< amp::ampf<Precision> > aeffective;
        ap::template_2d_array< amp::ampf<Precision> > aparam;
        ap::template_1d_array< amp::ampf<Precision> > xe;
        ap::template_1d_array< amp::ampf<Precision> > b;
        int n;
        int pass;
        int i;
        int j;
        int cnts;
        int cntu;
        int cntt;
        int cntm;
        bool waserrors;
        bool isupper;
        bool istrans;
        bool isunit;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s;


        waserrors = false;
        maxmn = 15;
        passcount = 15;
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Different problems
        //
        for(n=1; n<=maxmn; n++)
        {
            aeffective.setbounds(0, n-1, 0, n-1);
            aparam.setbounds(0, n-1, 0, n-1);
            xe.setbounds(0, n-1);
            b.setbounds(0, n-1);
            for(pass=1; pass<=passcount; pass++)
            {
                for(cnts=0; cnts<=1; cnts++)
                {
                    for(cntu=0; cntu<=1; cntu++)
                    {
                        for(cntt=0; cntt<=1; cntt++)
                        {
                            for(cntm=0; cntm<=2; cntm++)
                            {
                                isupper = cnts==0;
                                isunit = cntu==0;
                                istrans = cntt==0;
                                
                                //
                                // Skip meaningless combinations of parameters:
                                // (matrix is singular) AND (matrix is unit diagonal)
                                //
                                if( cntm==2 && isunit )
                                {
                                    continue;
                                }
                                
                                //
                                // Clear matrices
                                //
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        aeffective(i,j) = 0;
                                        aparam(i,j) = 0;
                                    }
                                }
                                
                                //
                                // Prepare matrices
                                //
                                if( isupper )
                                {
                                    for(i=0; i<=n-1; i++)
                                    {
                                        for(j=i; j<=n-1; j++)
                                        {
                                            aeffective(i,j) = amp::ampf<Precision>("0.9")*(2*amp::ampf<Precision>::getRandom()-1);
                                            aparam(i,j) = aeffective(i,j);
                                        }
                                        aeffective(i,i) = (2*ap::randominteger(2)-1)*(amp::ampf<Precision>("0.8")+amp::ampf<Precision>::getRandom());
                                        aparam(i,i) = aeffective(i,i);
                                    }
                                }
                                else
                                {
                                    for(i=0; i<=n-1; i++)
                                    {
                                        for(j=0; j<=i; j++)
                                        {
                                            aeffective(i,j) = amp::ampf<Precision>("0.9")*(2*amp::ampf<Precision>::getRandom()-1);
                                            aparam(i,j) = aeffective(i,j);
                                        }
                                        aeffective(i,i) = (2*ap::randominteger(2)-1)*(amp::ampf<Precision>("0.8")+amp::ampf<Precision>::getRandom());
                                        aparam(i,i) = aeffective(i,i);
                                    }
                                }
                                if( isunit )
                                {
                                    for(i=0; i<=n-1; i++)
                                    {
                                        aeffective(i,i) = 1;
                                        aparam(i,i) = 0;
                                    }
                                }
                                if( istrans )
                                {
                                    if( isupper )
                                    {
                                        for(i=0; i<=n-1; i++)
                                        {
                                            for(j=i+1; j<=n-1; j++)
                                            {
                                                aeffective(j,i) = aeffective(i,j);
                                                aeffective(i,j) = 0;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        for(i=0; i<=n-1; i++)
                                        {
                                            for(j=i+1; j<=n-1; j++)
                                            {
                                                aeffective(i,j) = aeffective(j,i);
                                                aeffective(j,i) = 0;
                                            }
                                        }
                                    }
                                }
                                
                                //
                                // Prepare task, solve, compare
                                //
                                for(i=0; i<=n-1; i++)
                                {
                                    xe(i) = 2*amp::ampf<Precision>::getRandom()-1;
                                }
                                for(i=0; i<=n-1; i++)
                                {
                                    v = amp::vdotproduct(aeffective.getrow(i, 0, n-1), xe.getvector(0, n-1));
                                    b(i) = v;
                                }
                                trlinsolve::rmatrixtrsafesolve<Precision>(aparam, n, b, s, isupper, istrans, isunit);
                                amp::vmul(xe.getvector(0, n-1), s);
                                amp::vsub(xe.getvector(0, n-1), b.getvector(0, n-1));
                                v = amp::vdotproduct(xe.getvector(0, n-1), xe.getvector(0, n-1));
                                v = amp::sqrt<Precision>(v);
                                waserrors = waserrors || v>threshold;
                            }
                        }
                    }
                }
            }
        }
        
        //
        // report
        //
        if( !silent )
        {
            printf("TESTING RMatrixTRSafeSolve\n");
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
    Copy
    *************************************************************************/
    template<unsigned int Precision>
    void makeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& b)
    {
        int i;
        int j;


        b.setbounds(0, m-1, 0, n-1);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                b(i,j) = a(i,j);
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testsstunit_test_silent()
    {
        bool result;


        result = testsst<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testsstunit_test()
    {
        bool result;


        result = testsst<Precision>(false);
        return result;
    }
} // namespace

#endif
