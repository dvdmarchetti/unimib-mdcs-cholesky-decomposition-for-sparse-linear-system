
#ifndef _testsafesolveunit_h
#define _testsafesolveunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "safesolve.h"
namespace testsafesolveunit
{
    template<unsigned int Precision>
    bool testsafesolve(bool silent);
    template<unsigned int Precision>
    void rmatrixmakeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& b);
    template<unsigned int Precision>
    void cmatrixmakeacopy(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& b);
    template<unsigned int Precision>
    bool testsafesolveunit_test_silent();
    template<unsigned int Precision>
    bool testsafesolveunit_test();


    /*************************************************************************
    Main unittest subroutine
    *************************************************************************/
    template<unsigned int Precision>
    bool testsafesolve(bool silent)
    {
        bool result;
        int maxmn;
        amp::ampf<Precision> threshold;
        bool rerrors;
        bool cerrors;
        bool waserrors;
        bool isupper;
        int trans;
        bool isunit;
        amp::ampf<Precision> scalea;
        amp::ampf<Precision> growth;
        int i;
        int j;
        int n;
        int j1;
        int j2;
        amp::campf<Precision> cv;
        ap::template_2d_array< amp::campf<Precision> > ca;
        ap::template_2d_array< amp::campf<Precision> > cea;
        ap::template_2d_array< amp::campf<Precision> > ctmpa;
        ap::template_1d_array< amp::campf<Precision> > cxs;
        ap::template_1d_array< amp::campf<Precision> > cxe;
        amp::ampf<Precision> rv;
        ap::template_2d_array< amp::ampf<Precision> > ra;
        ap::template_2d_array< amp::ampf<Precision> > rea;
        ap::template_2d_array< amp::ampf<Precision> > rtmpa;
        ap::template_1d_array< amp::ampf<Precision> > rxs;
        ap::template_1d_array< amp::ampf<Precision> > rxe;
        int i_;


        maxmn = 30;
        threshold = 100000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        rerrors = false;
        cerrors = false;
        waserrors = false;
        
        //
        // Different problems: general tests
        //
        for(n=1; n<=maxmn; n++)
        {
            
            //
            // test complex solver with well-conditioned matrix:
            // 1. generate A: fill off-diagonal elements with small values,
            //    diagonal elements are filled with larger values
            // 2. generate 'effective' A
            // 3. prepare task (exact X is stored in CXE, right part - in CXS),
            //    solve and compare CXS and CXE
            //
            isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
            trans = ap::randominteger(3);
            isunit = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
            scalea = amp::ampf<Precision>::getRandom()+amp::ampf<Precision>("0.5");
            ca.setlength(n, n);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( i==j )
                    {
                        ca(i,j).x = (2*ap::randominteger(2)-1)*(5+amp::ampf<Precision>::getRandom());
                        ca(i,j).y = (2*ap::randominteger(2)-1)*(5+amp::ampf<Precision>::getRandom());
                    }
                    else
                    {
                        ca(i,j).x = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                        ca(i,j).y = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                    }
                }
            }
            cmatrixmakeacopy<Precision>(ca, n, n, ctmpa);
            for(i=0; i<=n-1; i++)
            {
                if( isupper )
                {
                    j1 = 0;
                    j2 = i-1;
                }
                else
                {
                    j1 = i+1;
                    j2 = n-1;
                }
                for(j=j1; j<=j2; j++)
                {
                    ctmpa(i,j) = 0;
                }
                if( isunit )
                {
                    ctmpa(i,i) = 1;
                }
            }
            cea.setlength(n, n);
            for(i=0; i<=n-1; i++)
            {
                if( trans==0 )
                {
                    for(i_=0; i_<=n-1;i_++)
                    {
                        cea(i,i_) = scalea*ctmpa(i,i_);
                    }
                }
                if( trans==1 )
                {
                    for(i_=0; i_<=n-1;i_++)
                    {
                        cea(i_,i) = scalea*ctmpa(i,i_);
                    }
                }
                if( trans==2 )
                {
                    for(i_=0; i_<=n-1;i_++)
                    {
                        cea(i_,i) = scalea*amp::conj(ctmpa(i,i_));
                    }
                }
            }
            cxe.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                cxe(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                cxe(i).y = 2*amp::ampf<Precision>::getRandom()-1;
            }
            cxs.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                cv = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    cv += cea(i,i_)*cxe(i_);
                }
                cxs(i) = cv;
            }
            if( safesolve::cmatrixscaledtrsafesolve<Precision>(ca, scalea, n, cxs, isupper, trans, isunit, amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber())) )
            {
                for(i=0; i<=n-1; i++)
                {
                    cerrors = cerrors || amp::abscomplex<Precision>(cxs(i)-cxe(i))>threshold;
                }
            }
            else
            {
                cerrors = true;
            }
            
            //
            // same with real
            //
            isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
            trans = ap::randominteger(2);
            isunit = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
            scalea = amp::ampf<Precision>::getRandom()+amp::ampf<Precision>("0.5");
            ra.setlength(n, n);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( i==j )
                    {
                        ra(i,j) = (2*ap::randominteger(2)-1)*(5+amp::ampf<Precision>::getRandom());
                    }
                    else
                    {
                        ra(i,j) = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                    }
                }
            }
            rmatrixmakeacopy<Precision>(ra, n, n, rtmpa);
            for(i=0; i<=n-1; i++)
            {
                if( isupper )
                {
                    j1 = 0;
                    j2 = i-1;
                }
                else
                {
                    j1 = i+1;
                    j2 = n-1;
                }
                for(j=j1; j<=j2; j++)
                {
                    rtmpa(i,j) = 0;
                }
                if( isunit )
                {
                    rtmpa(i,i) = 1;
                }
            }
            rea.setlength(n, n);
            for(i=0; i<=n-1; i++)
            {
                if( trans==0 )
                {
                    amp::vmove(rea.getrow(i, 0, n-1), rtmpa.getrow(i, 0, n-1), scalea);
                }
                if( trans==1 )
                {
                    amp::vmove(rea.getcolumn(i, 0, n-1), rtmpa.getrow(i, 0, n-1), scalea);
                }
            }
            rxe.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                rxe(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            rxs.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                rv = amp::vdotproduct(rea.getrow(i, 0, n-1), rxe.getvector(0, n-1));
                rxs(i) = rv;
            }
            if( safesolve::rmatrixscaledtrsafesolve<Precision>(ra, scalea, n, rxs, isupper, trans, isunit, amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber())) )
            {
                for(i=0; i<=n-1; i++)
                {
                    rerrors = rerrors || amp::abs<Precision>(rxs(i)-rxe(i))>threshold;
                }
            }
            else
            {
                rerrors = true;
            }
        }
        
        //
        // Special test with diagonal ill-conditioned matrix:
        // * ability to solve it when resulting growth is less than threshold
        // * ability to stop solve when resulting growth is greater than threshold
        //
        // A = diag(1, 1/growth)
        // b = (1, 0.5)
        //
        n = 2;
        growth = 10;
        ca.setlength(n, n);
        ca(0,0) = 1;
        ca(0,1) = 0;
        ca(1,0) = 0;
        ca(1,1) = 1/growth;
        cxs.setlength(n);
        cxs(0) = amp::ampf<Precision>("1.0");
        cxs(1) = amp::ampf<Precision>("0.5");
        cerrors = cerrors || !safesolve::cmatrixscaledtrsafesolve<Precision>(ca, amp::ampf<Precision>("1.0"), n, cxs, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), ap::randominteger(3), false, amp::ampf<Precision>("1.05")*amp::maximum<Precision>(amp::abscomplex<Precision>(cxs(1))*growth, amp::ampf<Precision>("1.0")));
        cerrors = cerrors || !safesolve::cmatrixscaledtrsafesolve<Precision>(ca, amp::ampf<Precision>("1.0"), n, cxs, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), ap::randominteger(3), false, amp::ampf<Precision>("0.95")*amp::maximum<Precision>(amp::abscomplex<Precision>(cxs(1))*growth, amp::ampf<Precision>("1.0")));
        ra.setlength(n, n);
        ra(0,0) = 1;
        ra(0,1) = 0;
        ra(1,0) = 0;
        ra(1,1) = 1/growth;
        rxs.setlength(n);
        rxs(0) = amp::ampf<Precision>("1.0");
        rxs(1) = amp::ampf<Precision>("0.5");
        rerrors = rerrors || !safesolve::rmatrixscaledtrsafesolve<Precision>(ra, amp::ampf<Precision>("1.0"), n, rxs, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), ap::randominteger(2), false, amp::ampf<Precision>("1.05")*amp::maximum<Precision>(amp::abs<Precision>(rxs(1))*growth, amp::ampf<Precision>("1.0")));
        rerrors = rerrors || !safesolve::rmatrixscaledtrsafesolve<Precision>(ra, amp::ampf<Precision>("1.0"), n, rxs, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), ap::randominteger(2), false, amp::ampf<Precision>("0.95")*amp::maximum<Precision>(amp::abs<Precision>(rxs(1))*growth, amp::ampf<Precision>("1.0")));
        
        //
        // Special test with diagonal degenerate matrix:
        // * ability to solve it when resulting growth is less than threshold
        // * ability to stop solve when resulting growth is greater than threshold
        //
        // A = diag(1, 0)
        // b = (1, 0.5)
        //
        n = 2;
        ca.setlength(n, n);
        ca(0,0) = 1;
        ca(0,1) = 0;
        ca(1,0) = 0;
        ca(1,1) = 0;
        cxs.setlength(n);
        cxs(0) = amp::ampf<Precision>("1.0");
        cxs(1) = amp::ampf<Precision>("0.5");
        cerrors = cerrors || safesolve::cmatrixscaledtrsafesolve<Precision>(ca, amp::ampf<Precision>("1.0"), n, cxs, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), ap::randominteger(3), false, amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber()));
        ra.setlength(n, n);
        ra(0,0) = 1;
        ra(0,1) = 0;
        ra(1,0) = 0;
        ra(1,1) = 0;
        rxs.setlength(n);
        rxs(0) = amp::ampf<Precision>("1.0");
        rxs(1) = amp::ampf<Precision>("0.5");
        rerrors = rerrors || safesolve::rmatrixscaledtrsafesolve<Precision>(ra, amp::ampf<Precision>("1.0"), n, rxs, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), ap::randominteger(2), false, amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber()));
        
        //
        // report
        //
        waserrors = rerrors || cerrors;
        if( !silent )
        {
            printf("TESTING SAFE TR SOLVER\n");
            printf("REAL:                                    ");
            if( !rerrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("COMPLEX:                                 ");
            if( !cerrors )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
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
    Copy
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixmakeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
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
    Copy
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixmakeacopy(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& b)
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
    bool testsafesolveunit_test_silent()
    {
        bool result;


        result = testsafesolve<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testsafesolveunit_test()
    {
        bool result;


        result = testsafesolve<Precision>(false);
        return result;
    }
} // namespace

#endif
