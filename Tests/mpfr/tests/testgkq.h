
#ifndef _testgkq_h
#define _testgkq_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "tsort.h"
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
#include "gkq.h"
namespace testgkq
{
    template<unsigned int Precision>
    bool testgkqunit(bool silent);
    template<unsigned int Precision>
    amp::ampf<Precision> mapkind(int k);
    template<unsigned int Precision>
    bool testgkq_test_silent();
    template<unsigned int Precision>
    bool testgkq_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testgkqunit(bool silent)
    {
        bool result;
        int pkind;
        amp::ampf<Precision> errtol;
        amp::ampf<Precision> eps;
        amp::ampf<Precision> nonstricterrtol;
        int n;
        int i;
        int k;
        int info;
        amp::ampf<Precision> err;
        int akind;
        int bkind;
        amp::ampf<Precision> alphac;
        amp::ampf<Precision> betac;
        ap::template_1d_array< amp::ampf<Precision> > x1;
        ap::template_1d_array< amp::ampf<Precision> > wg1;
        ap::template_1d_array< amp::ampf<Precision> > wk1;
        ap::template_1d_array< amp::ampf<Precision> > x2;
        ap::template_1d_array< amp::ampf<Precision> > wg2;
        ap::template_1d_array< amp::ampf<Precision> > wk2;
        int info1;
        int info2;
        bool successatleastonce;
        bool intblerrors;
        bool vstblerrors;
        bool generrors;
        bool waserrors;


        intblerrors = false;
        vstblerrors = false;
        generrors = false;
        waserrors = false;
        errtol = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        nonstricterrtol = 1000*errtol;
        
        //
        // test recurrence-based Legendre nodes against the precalculated table
        //
        for(pkind=0; pkind<=5; pkind++)
        {
            n = 0;
            if( pkind==0 )
            {
                n = 15;
            }
            if( pkind==1 )
            {
                n = 21;
            }
            if( pkind==2 )
            {
                n = 31;
            }
            if( pkind==3 )
            {
                n = 41;
            }
            if( pkind==4 )
            {
                n = 51;
            }
            if( pkind==5 )
            {
                n = 61;
            }
            gkq::gkqlegendrecalc<Precision>(n, info, x1, wk1, wg1);
            gkq::gkqlegendretbl<Precision>(n, x2, wk2, wg2, eps);
            if( info<=0 )
            {
                generrors = true;
                break;
            }
            for(i=0; i<=n-1; i++)
            {
                vstblerrors = vstblerrors || amp::abs<Precision>(x1(i)-x2(i))>errtol;
                vstblerrors = vstblerrors || amp::abs<Precision>(wk1(i)-wk2(i))>errtol;
                vstblerrors = vstblerrors || amp::abs<Precision>(wg1(i)-wg2(i))>errtol;
            }
        }
        
        //
        // Test recurrence-baced Gauss-Kronrod nodes against Gauss-only nodes
        // calculated with subroutines from GQ unit.
        //
        for(k=1; k<=30; k++)
        {
            n = 2*k+1;
            
            //
            // Gauss-Legendre
            //
            err = 0;
            gkq::gkqgenerategausslegendre<Precision>(n, info1, x1, wk1, wg1);
            gq::gqgenerategausslegendre<Precision>(k, info2, x2, wg2);
            if( info1>0 && info2>0 )
            {
                for(i=0; i<=k-1; i++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(x1(2*i+1)-x2(i)));
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(wg1(2*i+1)-wg2(i)));
                }
            }
            else
            {
                generrors = true;
            }
            generrors = generrors || err>errtol;
        }
        for(k=1; k<=15; k++)
        {
            n = 2*k+1;
            
            //
            // Gauss-Jacobi
            //
            successatleastonce = false;
            err = 0;
            for(akind=0; akind<=9; akind++)
            {
                for(bkind=0; bkind<=9; bkind++)
                {
                    alphac = mapkind<Precision>(akind);
                    betac = mapkind<Precision>(bkind);
                    gkq::gkqgenerategaussjacobi<Precision>(n, alphac, betac, info1, x1, wk1, wg1);
                    gq::gqgenerategaussjacobi<Precision>(k, alphac, betac, info2, x2, wg2);
                    if( info1>0 && info2>0 )
                    {
                        successatleastonce = true;
                        for(i=0; i<=k-1; i++)
                        {
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(x1(2*i+1)-x2(i)));
                            err = amp::maximum<Precision>(err, amp::abs<Precision>(wg1(2*i+1)-wg2(i)));
                        }
                    }
                    else
                    {
                        generrors = generrors || info1!=-5;
                    }
                }
            }
            generrors = generrors || err>errtol || !successatleastonce;
        }
        
        //
        // end
        //
        waserrors = intblerrors || vstblerrors || generrors;
        if( !silent )
        {
            printf("TESTING GAUSS-KRONROD QUADRATURES\n");
            printf("FINAL RESULT:                             ");
            if( waserrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* PRE-CALCULATED TABLE:                   ");
            if( intblerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* CALCULATED AGAINST THE TABLE:           ");
            if( vstblerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* GENERAL PROPERTIES:                     ");
            if( generrors )
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
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testgkq_test_silent()
    {
        bool result;


        result = testgkqunit<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testgkq_test()
    {
        bool result;


        result = testgkqunit<Precision>(false);
        return result;
    }
} // namespace

#endif
