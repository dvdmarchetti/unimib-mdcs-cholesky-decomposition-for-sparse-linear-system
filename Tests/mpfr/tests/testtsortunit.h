
#ifndef _testtsortunit_h
#define _testtsortunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "tsort.h"
namespace testtsortunit
{
    template<unsigned int Precision>
    bool testsort(bool silent);
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void unset1di(ap::template_1d_array< int >& a);
    template<unsigned int Precision>
    void testsortresults(const ap::template_1d_array< amp::ampf<Precision> >& asorted,
        const ap::template_1d_array< int >& p1,
        const ap::template_1d_array< int >& p2,
        const ap::template_1d_array< amp::ampf<Precision> >& aoriginal,
        int n,
        bool& waserrors);
    template<unsigned int Precision>
    bool testtsortunit_test_silent();
    template<unsigned int Precision>
    bool testtsortunit_test();


    /*************************************************************************
    Testing tag sort
    *************************************************************************/
    template<unsigned int Precision>
    bool testsort(bool silent)
    {
        bool result;
        bool waserrors;
        int n;
        int i;
        int pass;
        int passcount;
        int maxn;
        ap::template_1d_array< amp::ampf<Precision> > a;
        ap::template_1d_array< amp::ampf<Precision> > a0;
        ap::template_1d_array< amp::ampf<Precision> > a2;
        ap::template_1d_array< int > p1;
        ap::template_1d_array< int > p2;


        waserrors = false;
        maxn = 100;
        passcount = 10;
        
        //
        // Test
        //
        for(n=1; n<=maxn; n++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // (probably) distinct sort
                //
                unset1di<Precision>(p1);
                unset1di<Precision>(p2);
                a.setbounds(0, n-1);
                a0.setbounds(0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    a(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    a0(i) = a(i);
                }
                tsort::tagsort<Precision>(a0, n, p1, p2);
                testsortresults<Precision>(a0, p1, p2, a, n, waserrors);
                
                //
                // non-distinct sort
                //
                unset1di<Precision>(p1);
                unset1di<Precision>(p2);
                a.setbounds(0, n-1);
                a0.setbounds(0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    a(i) = i/2;
                    a0(i) = a(i);
                }
                tsort::tagsort<Precision>(a0, n, p1, p2);
                testsortresults<Precision>(a0, p1, p2, a, n, waserrors);
            }
        }
        
        //
        // report
        //
        if( !silent )
        {
            printf("TESTING TAGSORT\n");
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
    Unsets 2D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::campf<Precision> >& a)
    {
        a.setbounds(0, 0, 0, 0);
        a(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets 1D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& a)
    {
        a.setbounds(0, 0);
        a(0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets 1D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset1di(ap::template_1d_array< int >& a)
    {
        a.setbounds(0, 0);
        a(0) = ap::randominteger(3)-1;
    }


    template<unsigned int Precision>
    void testsortresults(const ap::template_1d_array< amp::ampf<Precision> >& asorted,
        const ap::template_1d_array< int >& p1,
        const ap::template_1d_array< int >& p2,
        const ap::template_1d_array< amp::ampf<Precision> >& aoriginal,
        int n,
        bool& waserrors)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > a2;
        amp::ampf<Precision> t;
        ap::template_1d_array< int > f;


        a2.setbounds(0, n-1);
        f.setbounds(0, n-1);
        
        //
        // is set ordered?
        //
        for(i=0; i<=n-2; i++)
        {
            waserrors = waserrors || asorted(i)>asorted(i+1);
        }
        
        //
        // P1 correctness
        //
        for(i=0; i<=n-1; i++)
        {
            waserrors = waserrors || asorted(i)!=aoriginal(p1(i));
        }
        for(i=0; i<=n-1; i++)
        {
            f(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            f(p1(i)) = f(p1(i))+1;
        }
        for(i=0; i<=n-1; i++)
        {
            waserrors = waserrors || f(i)!=1;
        }
        
        //
        // P2 correctness
        //
        for(i=0; i<=n-1; i++)
        {
            a2(i) = aoriginal(i);
        }
        for(i=0; i<=n-1; i++)
        {
            if( p2(i)!=i )
            {
                t = a2(i);
                a2(i) = a2(p2(i));
                a2(p2(i)) = t;
            }
        }
        for(i=0; i<=n-1; i++)
        {
            waserrors = waserrors || asorted(i)!=a2(i);
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testtsortunit_test_silent()
    {
        bool result;


        result = testsort<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testtsortunit_test()
    {
        bool result;


        result = testsort<Precision>(false);
        return result;
    }
} // namespace

#endif
