
#ifndef _testinvldltunit_h
#define _testinvldltunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "sblas.h"
#include "ldlt.h"
#include "sinverse.h"
namespace testinvldltunit
{
    template<unsigned int Precision>
    bool testinvldlt(bool silent);
    template<unsigned int Precision>
    void generatematrix(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        int task);
    template<unsigned int Precision>
    void restorematrix(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool upperin);
    template<unsigned int Precision>
    bool testinvldltunit_test_silent();
    template<unsigned int Precision>
    bool testinvldltunit_test();


    template<unsigned int Precision>
    bool testinvldlt(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        ap::template_2d_array< amp::ampf<Precision> > a3;
        ap::template_1d_array< int > p;
        int n;
        int pass;
        int mtask;
        int i;
        int j;
        int k;
        int minij;
        bool upperin;
        bool cr;
        amp::ampf<Precision> v;
        amp::ampf<Precision> err;
        bool waserrors;
        int passcount;
        int maxn;
        int htask;
        amp::ampf<Precision> threshold;


        err = 0;
        passcount = 10;
        maxn = 20;
        threshold = 10000000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        waserrors = false;
        
        //
        // Test
        //
        for(n=1; n<=maxn; n++)
        {
            a.setbounds(0, n-1, 0, n-1);
            a2.setbounds(0, n-1, 0, n-1);
            a3.setbounds(0, n-1, 0, n-1);
            for(mtask=2; mtask<=2; mtask++)
            {
                for(htask=0; htask<=1; htask++)
                {
                    for(pass=1; pass<=passcount; pass++)
                    {
                        upperin = htask==0;
                        
                        //
                        // Prepare task:
                        // * A contains symmetric matrix
                        // * A2, A3 contains its upper (or lower) half
                        //
                        generatematrix<Precision>(a, n, mtask);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a2(i,j) = a(i,j);
                                a3(i,j) = a(i,j);
                            }
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                if( upperin )
                                {
                                    if( j<i )
                                    {
                                        a2(i,j) = 0;
                                        a3(i,j) = 0;
                                    }
                                }
                                else
                                {
                                    if( i<j )
                                    {
                                        a2(i,j) = 0;
                                        a3(i,j) = 0;
                                    }
                                }
                            }
                        }
                        
                        //
                        // Test 1: inv(A2)
                        //
                        sinverse::smatrixinverse<Precision>(a2, n, upperin);
                        restorematrix<Precision>(a2, n, upperin);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                v = amp::vdotproduct(a.getrow(i, 0, n-1), a2.getcolumn(j, 0, n-1));
                                if( i==j )
                                {
                                    v = v-1;
                                }
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(v));
                            }
                        }
                        
                        //
                        // Test 2: inv(LDLt(A3))
                        //
                        ldlt::smatrixldlt<Precision>(a3, n, upperin, p);
                        sinverse::smatrixldltinverse<Precision>(a3, p, n, upperin);
                        restorematrix<Precision>(a3, n, upperin);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                v = amp::vdotproduct(a.getrow(i, 0, n-1), a3.getcolumn(j, 0, n-1));
                                if( i==j )
                                {
                                    v = v-1;
                                }
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(v));
                            }
                        }
                    }
                }
            }
        }
        
        //
        // report
        //
        waserrors = err>threshold;
        if( !silent )
        {
            printf("TESTING LDLT INVERSE\n");
            printf("ERROR:                                   %5.3le\n",
                double(amp::ampf<Precision>(err).toDouble()));
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


    template<unsigned int Precision>
    void generatematrix(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        int task)
    {
        int i;
        int j;


        if( task==0 )
        {
            
            //
            // Zero matrix
            //
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a(i,j) = 0;
                }
            }
        }
        if( task==1 )
        {
            
            //
            // Sparse matrix
            //
            for(i=0; i<=n-1; i++)
            {
                for(j=i+1; j<=n-1; j++)
                {
                    if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.95") )
                    {
                        a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    else
                    {
                        a(i,j) = 0;
                    }
                    a(j,i) = a(i,j);
                }
                if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.95") )
                {
                    a(i,i) = (2*ap::randominteger(2)-1)*(amp::ampf<Precision>("0.8")+amp::ampf<Precision>::getRandom());
                }
                else
                {
                    a(i,i) = 0;
                }
            }
        }
        if( task==2 )
        {
            
            //
            // Dense matrix
            //
            for(i=0; i<=n-1; i++)
            {
                for(j=i+1; j<=n-1; j++)
                {
                    a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    a(j,i) = a(i,j);
                }
                a(i,i) = (2*ap::randominteger(2)-1)*(amp::ampf<Precision>("0.7")+amp::ampf<Precision>::getRandom());
            }
        }
    }


    template<unsigned int Precision>
    void restorematrix(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool upperin)
    {
        int i;
        int j;


        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( upperin )
                {
                    if( j<i )
                    {
                        a(i,j) = a(j,i);
                    }
                }
                else
                {
                    if( i<j )
                    {
                        a(i,j) = a(j,i);
                    }
                }
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testinvldltunit_test_silent()
    {
        bool result;


        result = testinvldlt<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testinvldltunit_test()
    {
        bool result;


        result = testinvldlt<Precision>(false);
        return result;
    }
} // namespace

#endif
