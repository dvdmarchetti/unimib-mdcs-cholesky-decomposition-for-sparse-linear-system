
#ifndef _testldltunit_h
#define _testldltunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "ldlt.h"
namespace testldltunit
{
    template<unsigned int Precision>
    bool testldlt(bool silent);
    template<unsigned int Precision>
    void generatematrix(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        int task);
    template<unsigned int Precision>
    bool testldltunit_test_silent();
    template<unsigned int Precision>
    bool testldltunit_test();


    template<unsigned int Precision>
    bool testldlt(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        ap::template_2d_array< amp::ampf<Precision> > l;
        ap::template_2d_array< amp::ampf<Precision> > d;
        ap::template_2d_array< amp::ampf<Precision> > u;
        ap::template_2d_array< amp::ampf<Precision> > t;
        ap::template_2d_array< amp::ampf<Precision> > t2;
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
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        waserrors = false;
        
        //
        // Test
        //
        for(n=1; n<=maxn; n++)
        {
            a.setbounds(0, n-1, 0, n-1);
            a2.setbounds(0, n-1, 0, n-1);
            l.setbounds(0, n-1, 0, n-1);
            u.setbounds(0, n-1, 0, n-1);
            d.setbounds(0, n-1, 0, n-1);
            t.setbounds(0, n-1, 0, n-1);
            t2.setbounds(0, n-1, 0, n-1);
            for(mtask=0; mtask<=2; mtask++)
            {
                for(htask=0; htask<=1; htask++)
                {
                    for(pass=1; pass<=passcount; pass++)
                    {
                        upperin = htask==0;
                        
                        //
                        // Prepare task:
                        // * A contains symmetric matrix
                        // * A2 contains its upper (or lower) half
                        //
                        generatematrix<Precision>(a, n, mtask);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a2(i,j) = a(i,j);
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
                                    }
                                }
                                else
                                {
                                    if( i<j )
                                    {
                                        a2(i,j) = 0;
                                    }
                                }
                            }
                        }
                        
                        //
                        // LDLt
                        //
                        ldlt::smatrixldlt<Precision>(a2, n, upperin, p);
                        
                        //
                        // Test (upper or lower)
                        //
                        if( upperin )
                        {
                            
                            //
                            // Unpack D
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    d(i,j) = 0;
                                }
                            }
                            k = 0;
                            while( k<n )
                            {
                                if( p(k)>=0 )
                                {
                                    d(k,k) = a2(k,k);
                                    k = k+1;
                                }
                                else
                                {
                                    d(k,k) = a2(k,k);
                                    d(k,k+1) = a2(k,k+1);
                                    d(k+1,k) = a2(k,k+1);
                                    d(k+1,k+1) = a2(k+1,k+1);
                                    k = k+2;
                                }
                            }
                            
                            //
                            // Unpack U
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    u(i,j) = 0;
                                }
                                u(i,i) = 1;
                            }
                            k = 0;
                            while( k<n )
                            {
                                
                                //
                                // unpack Uk
                                //
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        t(i,j) = 0;
                                    }
                                    t(i,i) = 1;
                                }
                                if( p(k)>=0 )
                                {
                                    for(i=0; i<=k-1; i++)
                                    {
                                        t(i,k) = a2(i,k);
                                    }
                                }
                                else
                                {
                                    for(i=0; i<=k-1; i++)
                                    {
                                        t(i,k) = a2(i,k);
                                        t(i,k+1) = a2(i,k+1);
                                    }
                                }
                                
                                //
                                // multiply U
                                //
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        v = amp::vdotproduct(t.getrow(i, 0, n-1), u.getcolumn(j, 0, n-1));
                                        t2(i,j) = v;
                                    }
                                }
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        u(i,j) = t2(i,j);
                                    }
                                }
                                
                                //
                                // permutations
                                //
                                if( p(k)>=0 )
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        v = u(k,j);
                                        u(k,j) = u(p(k),j);
                                        u(p(k),j) = v;
                                    }
                                    k = k+1;
                                }
                                else
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        v = u(k,j);
                                        u(k,j) = u(n+p(k),j);
                                        u(n+p(k),j) = v;
                                    }
                                    k = k+2;
                                }
                            }
                            
                            //
                            // Calculate U*D*U', store result in T2
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    v = amp::vdotproduct(u.getrow(i, 0, n-1), d.getcolumn(j, 0, n-1));
                                    t(i,j) = v;
                                }
                            }
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    v = amp::vdotproduct(t.getrow(i, 0, n-1), u.getrow(j, 0, n-1));
                                    t2(i,j) = v;
                                }
                            }
                        }
                        else
                        {
                            
                            //
                            // Unpack D
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    d(i,j) = 0;
                                }
                            }
                            k = 0;
                            while( k<n )
                            {
                                if( p(k)>=0 )
                                {
                                    d(k,k) = a2(k,k);
                                    k = k+1;
                                }
                                else
                                {
                                    d(k,k) = a2(k,k);
                                    d(k,k+1) = a2(k+1,k);
                                    d(k+1,k) = a2(k+1,k);
                                    d(k+1,k+1) = a2(k+1,k+1);
                                    k = k+2;
                                }
                            }
                            
                            //
                            // Unpack L
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    l(i,j) = 0;
                                }
                                l(i,i) = 1;
                            }
                            k = 0;
                            while( k<n )
                            {
                                
                                //
                                // permutations
                                //
                                if( p(k)>=0 )
                                {
                                    for(i=0; i<=n-1; i++)
                                    {
                                        v = l(i,k);
                                        l(i,k) = l(i,p(k));
                                        l(i,p(k)) = v;
                                    }
                                }
                                else
                                {
                                    for(i=0; i<=n-1; i++)
                                    {
                                        v = l(i,k+1);
                                        l(i,k+1) = l(i,n+p(k+1));
                                        l(i,n+p(k+1)) = v;
                                    }
                                }
                                
                                //
                                // unpack Lk
                                //
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        t(i,j) = 0;
                                    }
                                    t(i,i) = 1;
                                }
                                if( p(k)>=0 )
                                {
                                    for(i=k+1; i<=n-1; i++)
                                    {
                                        t(i,k) = a2(i,k);
                                    }
                                }
                                else
                                {
                                    for(i=k+2; i<=n-1; i++)
                                    {
                                        t(i,k) = a2(i,k);
                                        t(i,k+1) = a2(i,k+1);
                                    }
                                }
                                
                                //
                                // multiply L
                                //
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        v = amp::vdotproduct(l.getrow(i, 0, n-1), t.getcolumn(j, 0, n-1));
                                        t2(i,j) = v;
                                    }
                                }
                                for(i=0; i<=n-1; i++)
                                {
                                    for(j=0; j<=n-1; j++)
                                    {
                                        l(i,j) = t2(i,j);
                                    }
                                }
                                
                                //
                                // Next K
                                //
                                if( p(k)>=0 )
                                {
                                    k = k+1;
                                }
                                else
                                {
                                    k = k+2;
                                }
                            }
                            
                            //
                            // Calculate L*D*L', store result in T2
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    v = amp::vdotproduct(l.getrow(i, 0, n-1), d.getcolumn(j, 0, n-1));
                                    t(i,j) = v;
                                }
                            }
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    v = amp::vdotproduct(t.getrow(i, 0, n-1), l.getrow(j, 0, n-1));
                                    t2(i,j) = v;
                                }
                            }
                        }
                        
                        //
                        // Test
                        //
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                err = amp::maximum<Precision>(err, amp::abs<Precision>(a(i,j)-t2(i,j)));
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
            printf("TESTING LDLT DECOMPOSITION\n");
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
                a(i,i) = (2*ap::randominteger(2)-1)*(amp::ampf<Precision>("0.8")+amp::ampf<Precision>::getRandom());
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testldltunit_test_silent()
    {
        bool result;


        result = testldlt<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testldltunit_test()
    {
        bool result;


        result = testldlt<Precision>(false);
        return result;
    }
} // namespace

#endif
