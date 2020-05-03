
#ifndef _testmatinvunit_h
#define _testmatinvunit_h

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
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
namespace testmatinvunit
{
    template<unsigned int Precision>
    bool testmatinv(bool silent);
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
    bool rmatrixcheckinverse(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep);
    template<unsigned int Precision>
    bool spdmatrixcheckinverse(ap::template_2d_array< amp::ampf<Precision> > a,
        ap::template_2d_array< amp::ampf<Precision> > inva,
        bool isupper,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep);
    template<unsigned int Precision>
    bool hpdmatrixcheckinverse(ap::template_2d_array< amp::campf<Precision> > a,
        ap::template_2d_array< amp::campf<Precision> > inva,
        bool isupper,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep);
    template<unsigned int Precision>
    bool rmatrixcheckinversesingular(const ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep);
    template<unsigned int Precision>
    bool cmatrixcheckinverse(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_2d_array< amp::campf<Precision> >& inva,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep);
    template<unsigned int Precision>
    bool cmatrixcheckinversesingular(const ap::template_2d_array< amp::campf<Precision> >& inva,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void rmatrixdrophalf(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool droplower);
    template<unsigned int Precision>
    void cmatrixdrophalf(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool droplower);
    template<unsigned int Precision>
    void testrtrinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& rtrerrors);
    template<unsigned int Precision>
    void testctrinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& ctrerrors);
    template<unsigned int Precision>
    void testrinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& rerrors);
    template<unsigned int Precision>
    void testcinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& cerrors);
    template<unsigned int Precision>
    void testspdinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& spderrors);
    template<unsigned int Precision>
    void testhpdinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& hpderrors);
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void cunset2d(ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void cunset1d(ap::template_1d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void unsetrep(matinv::matinvreport<Precision>& r);
    template<unsigned int Precision>
    bool testmatinvunit_test_silent();
    template<unsigned int Precision>
    bool testmatinvunit_test();


    /*************************************************************************
    Test
    *************************************************************************/
    template<unsigned int Precision>
    bool testmatinv(bool silent)
    {
        bool result;
        int maxrn;
        int maxcn;
        int passcount;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> rcondtol;
        bool rtrerrors;
        bool ctrerrors;
        bool rerrors;
        bool cerrors;
        bool spderrors;
        bool hpderrors;
        bool waserrors;
        ap::template_2d_array< amp::ampf<Precision> > emptyra;
        ap::template_2d_array< amp::ampf<Precision> > emptyca;


        maxrn = 3*ablas::ablasblocksize<Precision>(emptyra)+1;
        maxcn = 3*ablas::ablasblocksize<Precision>(emptyca)+1;
        passcount = 1;
        threshold = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        rcondtol = amp::ampf<Precision>("0.01");
        rtrerrors = false;
        ctrerrors = false;
        rerrors = false;
        cerrors = false;
        spderrors = false;
        hpderrors = false;
        testrtrinv<Precision>(maxrn, passcount, threshold, rtrerrors);
        testctrinv<Precision>(maxcn, passcount, threshold, ctrerrors);
        testrinv<Precision>(maxrn, passcount, threshold, rerrors);
        testspdinv<Precision>(maxrn, passcount, threshold, spderrors);
        testcinv<Precision>(maxcn, passcount, threshold, cerrors);
        testhpdinv<Precision>(maxcn, passcount, threshold, hpderrors);
        waserrors = rtrerrors || ctrerrors || rerrors || cerrors || spderrors || hpderrors;
        if( !silent )
        {
            printf("TESTING MATINV\n");
            printf("* REAL TRIANGULAR:                        ");
            if( rtrerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* COMPLEX TRIANGULAR:                     ");
            if( ctrerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* REAL:                                   ");
            if( rerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* COMPLEX:                                ");
            if( cerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* SPD:                                    ");
            if( spderrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* HPD:                                    ");
            if( hpderrors )
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
    Checks whether inverse is correct
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixcheckinverse(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep)
    {
        bool result;
        int i;
        int j;
        amp::ampf<Precision> v;


        result = true;
        if( info<=0 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.r1>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.rinf>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = amp::vdotproduct(a.getrow(i, 0, n-1), inva.getcolumn(j, 0, n-1));
                    if( i==j )
                    {
                        v = v-1;
                    }
                    result = result && amp::abs<Precision>(v)<=threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether inverse is correct
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool spdmatrixcheckinverse(ap::template_2d_array< amp::ampf<Precision> > a,
        ap::template_2d_array< amp::ampf<Precision> > inva,
        bool isupper,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep)
    {
        bool result;
        int i;
        int j;
        amp::ampf<Precision> v;


        for(i=0; i<=n-2; i++)
        {
            if( isupper )
            {
                amp::vmove(a.getcolumn(i, i+1, n-1), a.getrow(i, i+1, n-1));
                amp::vmove(inva.getcolumn(i, i+1, n-1), inva.getrow(i, i+1, n-1));
            }
            else
            {
                amp::vmove(a.getrow(i, i+1, n-1), a.getcolumn(i, i+1, n-1));
                amp::vmove(inva.getrow(i, i+1, n-1), inva.getcolumn(i, i+1, n-1));
            }
        }
        result = true;
        if( info<=0 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.r1>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.rinf>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = amp::vdotproduct(a.getrow(i, 0, n-1), inva.getcolumn(j, 0, n-1));
                    if( i==j )
                    {
                        v = v-1;
                    }
                    result = result && amp::abs<Precision>(v)<=threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether inverse is correct
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool hpdmatrixcheckinverse(ap::template_2d_array< amp::campf<Precision> > a,
        ap::template_2d_array< amp::campf<Precision> > inva,
        bool isupper,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep)
    {
        bool result;
        int i;
        int j;
        amp::campf<Precision> v;
        int i_;


        for(i=0; i<=n-2; i++)
        {
            if( isupper )
            {
                for(i_=i+1; i_<=n-1;i_++)
                {
                    a(i_,i) = amp::conj(a(i,i_));
                }
                for(i_=i+1; i_<=n-1;i_++)
                {
                    inva(i_,i) = amp::conj(inva(i,i_));
                }
            }
            else
            {
                for(i_=i+1; i_<=n-1;i_++)
                {
                    a(i,i_) = amp::conj(a(i_,i));
                }
                for(i_=i+1; i_<=n-1;i_++)
                {
                    inva(i,i_) = amp::conj(inva(i_,i));
                }
            }
        }
        result = true;
        if( info<=0 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.r1>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.rinf>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += a(i,i_)*inva(i_,j);
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    result = result && amp::abscomplex<Precision>(v)<=threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether inversion result indicate singular matrix
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixcheckinversesingular(const ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep)
    {
        bool result;
        int i;
        int j;


        result = true;
        if( info!=-3 && info!=1 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<0 || rep.r1>1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<0 || rep.rinf>1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            if( info==-3 )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        result = result && inva(i,j)==0;
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether inverse is correct
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixcheckinverse(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_2d_array< amp::campf<Precision> >& inva,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep)
    {
        bool result;
        int i;
        int j;
        amp::campf<Precision> v;
        int i_;


        result = true;
        if( info<=0 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.r1>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<100*amp::ampf<Precision>::getAlgoPascalEpsilon() || rep.rinf>1+1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += a(i,i_)*inva(i_,j);
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    result = result && amp::abscomplex<Precision>(v)<=threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Checks whether inversion result indicate singular matrix
    Returns True on success.
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixcheckinversesingular(const ap::template_2d_array< amp::campf<Precision> >& inva,
        int n,
        amp::ampf<Precision> threshold,
        int info,
        const matinv::matinvreport<Precision>& rep)
    {
        bool result;
        int i;
        int j;


        result = true;
        if( info!=-3 && info!=1 )
        {
            result = false;
        }
        else
        {
            result = result && !(rep.r1<0 || rep.r1>1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            result = result && !(rep.rinf<0 || rep.rinf>1000*amp::ampf<Precision>::getAlgoPascalEpsilon());
            if( info==-3 )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        result = result && inva(i,j)==0;
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Drops upper or lower half of the matrix - fills it by special pattern
    which may be used later to ensure that this part wasn't changed
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixdrophalf(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool droplower)
    {
        int i;
        int j;


        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( droplower && i>j || !droplower && i<j )
                {
                    a(i,j) = 1+2*i+3*j;
                }
            }
        }
    }


    /*************************************************************************
    Drops upper or lower half of the matrix - fills it by special pattern
    which may be used later to ensure that this part wasn't changed
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixdrophalf(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool droplower)
    {
        int i;
        int j;


        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( droplower && i>j || !droplower && i<j )
                {
                    a(i,j) = 1+2*i+3*j;
                }
            }
        }
    }


    /*************************************************************************
    Real TR inverse
    *************************************************************************/
    template<unsigned int Precision>
    void testrtrinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& rtrerrors)
    {
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > b;
        int n;
        int pass;
        int i;
        int j;
        int task;
        bool isupper;
        bool isunit;
        amp::ampf<Precision> v;
        bool waserrors;
        int info;
        matinv::matinvreport<Precision> rep;


        waserrors = false;
        
        //
        // Test
        //
        for(n=1; n<=maxn; n++)
        {
            a.setlength(n, n);
            b.setlength(n, n);
            for(task=0; task<=3; task++)
            {
                for(pass=1; pass<=passcount; pass++)
                {
                    
                    //
                    // Determine task
                    //
                    isupper = task%2==0;
                    isunit = task/2%2==0;
                    
                    //
                    // Generate matrix
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            if( i==j )
                            {
                                a(i,i) = 1+amp::ampf<Precision>::getRandom();
                            }
                            else
                            {
                                a(i,j) = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                            }
                            b(i,j) = a(i,j);
                        }
                    }
                    
                    //
                    // Inverse
                    //
                    matinv::rmatrixtrinverse<Precision>(b, n, isupper, isunit, info, rep);
                    if( info<=0 )
                    {
                        rtrerrors = true;
                        return;
                    }
                    
                    //
                    // Structural test
                    //
                    if( isunit )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            rtrerrors = rtrerrors || a(i,i)!=b(i,i);
                        }
                    }
                    if( isupper )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=i-1; j++)
                            {
                                rtrerrors = rtrerrors || a(i,j)!=b(i,j);
                            }
                        }
                    }
                    else
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=i+1; j<=n-1; j++)
                            {
                                rtrerrors = rtrerrors || a(i,j)!=b(i,j);
                            }
                        }
                    }
                    
                    //
                    // Inverse test
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            if( j<i && isupper || j>i && !isupper )
                            {
                                a(i,j) = 0;
                                b(i,j) = 0;
                            }
                        }
                    }
                    if( isunit )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            a(i,i) = 1;
                            b(i,i) = 1;
                        }
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            v = amp::vdotproduct(a.getrow(i, 0, n-1), b.getcolumn(j, 0, n-1));
                            if( j!=i )
                            {
                                rtrerrors = rtrerrors || amp::abs<Precision>(v)>threshold;
                            }
                            else
                            {
                                rtrerrors = rtrerrors || amp::abs<Precision>(v-1)>threshold;
                            }
                        }
                    }
                }
            }
        }
    }


    /*************************************************************************
    Complex TR inverse
    *************************************************************************/
    template<unsigned int Precision>
    void testctrinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& ctrerrors)
    {
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > b;
        int n;
        int pass;
        int i;
        int j;
        int task;
        bool isupper;
        bool isunit;
        amp::campf<Precision> v;
        bool waserrors;
        int info;
        matinv::matinvreport<Precision> rep;
        int i_;


        waserrors = false;
        
        //
        // Test
        //
        for(n=1; n<=maxn; n++)
        {
            a.setlength(n, n);
            b.setlength(n, n);
            for(task=0; task<=3; task++)
            {
                for(pass=1; pass<=passcount; pass++)
                {
                    
                    //
                    // Determine task
                    //
                    isupper = task%2==0;
                    isunit = task/2%2==0;
                    
                    //
                    // Generate matrix
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            if( i==j )
                            {
                                a(i,i).x = 1+amp::ampf<Precision>::getRandom();
                                a(i,i).y = 1+amp::ampf<Precision>::getRandom();
                            }
                            else
                            {
                                a(i,j).x = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                                a(i,j).y = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                            }
                            b(i,j) = a(i,j);
                        }
                    }
                    
                    //
                    // Inverse
                    //
                    matinv::cmatrixtrinverse<Precision>(b, n, isupper, isunit, info, rep);
                    if( info<=0 )
                    {
                        ctrerrors = true;
                        return;
                    }
                    
                    //
                    // Structural test
                    //
                    if( isunit )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            ctrerrors = ctrerrors || a(i,i)!=b(i,i);
                        }
                    }
                    if( isupper )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=i-1; j++)
                            {
                                ctrerrors = ctrerrors || a(i,j)!=b(i,j);
                            }
                        }
                    }
                    else
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=i+1; j<=n-1; j++)
                            {
                                ctrerrors = ctrerrors || a(i,j)!=b(i,j);
                            }
                        }
                    }
                    
                    //
                    // Inverse test
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            if( j<i && isupper || j>i && !isupper )
                            {
                                a(i,j) = 0;
                                b(i,j) = 0;
                            }
                        }
                    }
                    if( isunit )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            a(i,i) = 1;
                            b(i,i) = 1;
                        }
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            v = 0.0;
                            for(i_=0; i_<=n-1;i_++)
                            {
                                v += a(i,i_)*b(i_,j);
                            }
                            if( j!=i )
                            {
                                ctrerrors = ctrerrors || amp::abscomplex<Precision>(v)>threshold;
                            }
                            else
                            {
                                ctrerrors = ctrerrors || amp::abscomplex<Precision>(v-1)>threshold;
                            }
                        }
                    }
                }
            }
        }
    }


    /*************************************************************************
    Real test
    *************************************************************************/
    template<unsigned int Precision>
    void testrinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& rerrors)
    {
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > lua;
        ap::template_2d_array< amp::ampf<Precision> > inva;
        ap::template_2d_array< amp::ampf<Precision> > invlua;
        ap::template_1d_array< int > p;
        int i;
        int j;
        int k;
        int n;
        int pass;
        int taskkind;
        amp::ampf<Precision> v;
        int info;
        matinv::matinvreport<Precision> rep;


        
        //
        // General square matrices:
        // * test general solvers
        // * test least squares solver
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                
                //
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                matgen::rmatrixrndcond<Precision>(n, amp::ampf<Precision>(1000), a);
                rmatrixmakeacopy<Precision>(a, n, n, lua);
                trfac::rmatrixlu<Precision>(lua, n, n, p);
                rmatrixmakeacopy<Precision>(a, n, n, inva);
                rmatrixmakeacopy<Precision>(lua, n, n, invlua);
                info = 0;
                unsetrep<Precision>(rep);
                matinv::rmatrixinverse<Precision>(inva, n, info, rep);
                rerrors = rerrors || !rmatrixcheckinverse<Precision>(a, inva, n, threshold, info, rep);
                info = 0;
                unsetrep<Precision>(rep);
                matinv::rmatrixluinverse<Precision>(invlua, p, n, info, rep);
                rerrors = rerrors || !rmatrixcheckinverse<Precision>(a, invlua, n, threshold, info, rep);
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                //    * with equal rows/columns
                // 2. test different methods
                //
                for(taskkind=0; taskkind<=4; taskkind++)
                {
                    unset2d<Precision>(a);
                    if( taskkind==0 )
                    {
                        
                        //
                        // all zeros
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 0;
                            }
                        }
                    }
                    if( taskkind==1 )
                    {
                        
                        //
                        // there is zero column
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        amp::vmul(a.getcolumn(k, 0, n-1), 0);
                    }
                    if( taskkind==2 )
                    {
                        
                        //
                        // there is zero row
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        amp::vmul(a.getrow(k, 0, n-1), 0);
                    }
                    if( taskkind==3 )
                    {
                        
                        //
                        // equal columns
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        amp::vmove(a.getcolumn(0, 0, n-1), a.getcolumn(k, 0, n-1));
                    }
                    if( taskkind==4 )
                    {
                        
                        //
                        // equal rows
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        amp::vmove(a.getrow(0, 0, n-1), a.getrow(k, 0, n-1));
                    }
                    rmatrixmakeacopy<Precision>(a, n, n, lua);
                    trfac::rmatrixlu<Precision>(lua, n, n, p);
                    info = 0;
                    unsetrep<Precision>(rep);
                    matinv::rmatrixinverse<Precision>(a, n, info, rep);
                    rerrors = rerrors || !rmatrixcheckinversesingular<Precision>(a, n, threshold, info, rep);
                    info = 0;
                    unsetrep<Precision>(rep);
                    matinv::rmatrixluinverse<Precision>(lua, p, n, info, rep);
                    rerrors = rerrors || !rmatrixcheckinversesingular<Precision>(lua, n, threshold, info, rep);
                }
            }
        }
    }


    /*************************************************************************
    Complex test
    *************************************************************************/
    template<unsigned int Precision>
    void testcinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& cerrors)
    {
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > lua;
        ap::template_2d_array< amp::campf<Precision> > inva;
        ap::template_2d_array< amp::campf<Precision> > invlua;
        ap::template_1d_array< int > p;
        int i;
        int j;
        int k;
        int n;
        int pass;
        int taskkind;
        amp::ampf<Precision> v;
        int info;
        matinv::matinvreport<Precision> rep;
        int i_;


        
        //
        // General square matrices:
        // * test general solvers
        // * test least squares solver
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                
                //
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                matgen::cmatrixrndcond<Precision>(n, amp::ampf<Precision>(1000), a);
                cmatrixmakeacopy<Precision>(a, n, n, lua);
                trfac::cmatrixlu<Precision>(lua, n, n, p);
                cmatrixmakeacopy<Precision>(a, n, n, inva);
                cmatrixmakeacopy<Precision>(lua, n, n, invlua);
                info = 0;
                unsetrep<Precision>(rep);
                matinv::cmatrixinverse<Precision>(inva, n, info, rep);
                cerrors = cerrors || !cmatrixcheckinverse<Precision>(a, inva, n, threshold, info, rep);
                info = 0;
                unsetrep<Precision>(rep);
                matinv::cmatrixluinverse<Precision>(invlua, p, n, info, rep);
                cerrors = cerrors || !cmatrixcheckinverse<Precision>(a, invlua, n, threshold, info, rep);
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                //    * with equal rows/columns
                // 2. test different methods
                //
                for(taskkind=0; taskkind<=4; taskkind++)
                {
                    cunset2d<Precision>(a);
                    if( taskkind==0 )
                    {
                        
                        //
                        // all zeros
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 0;
                            }
                        }
                    }
                    if( taskkind==1 )
                    {
                        
                        //
                        // there is zero column
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            a(i_,k) = 0*a(i_,k);
                        }
                    }
                    if( taskkind==2 )
                    {
                        
                        //
                        // there is zero row
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            a(k,i_) = 0*a(k,i_);
                        }
                    }
                    if( taskkind==3 )
                    {
                        
                        //
                        // equal columns
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            a(i_,0) = a(i_,k);
                        }
                    }
                    if( taskkind==4 )
                    {
                        
                        //
                        // equal rows
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            a(0,i_) = a(k,i_);
                        }
                    }
                    cmatrixmakeacopy<Precision>(a, n, n, lua);
                    trfac::cmatrixlu<Precision>(lua, n, n, p);
                    info = 0;
                    unsetrep<Precision>(rep);
                    matinv::cmatrixinverse<Precision>(a, n, info, rep);
                    cerrors = cerrors || !cmatrixcheckinversesingular<Precision>(a, n, threshold, info, rep);
                    info = 0;
                    unsetrep<Precision>(rep);
                    matinv::cmatrixluinverse<Precision>(lua, p, n, info, rep);
                    cerrors = cerrors || !cmatrixcheckinversesingular<Precision>(lua, n, threshold, info, rep);
                }
            }
        }
    }


    /*************************************************************************
    SPD test
    *************************************************************************/
    template<unsigned int Precision>
    void testspdinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& spderrors)
    {
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > cha;
        ap::template_2d_array< amp::ampf<Precision> > inva;
        ap::template_2d_array< amp::ampf<Precision> > invcha;
        bool isupper;
        int i;
        int j;
        int k;
        int n;
        int pass;
        int taskkind;
        amp::ampf<Precision> v;
        int info;
        matinv::matinvreport<Precision> rep;


        
        //
        // General square matrices:
        // * test general solvers
        // * test least squares solver
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
                
                //
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                matgen::spdmatrixrndcond<Precision>(n, amp::ampf<Precision>(1000), a);
                rmatrixdrophalf<Precision>(a, n, isupper);
                rmatrixmakeacopy<Precision>(a, n, n, cha);
                if( !trfac::spdmatrixcholesky<Precision>(cha, n, isupper) )
                {
                    continue;
                }
                rmatrixmakeacopy<Precision>(a, n, n, inva);
                rmatrixmakeacopy<Precision>(cha, n, n, invcha);
                info = 0;
                unsetrep<Precision>(rep);
                matinv::spdmatrixinverse<Precision>(inva, n, isupper, info, rep);
                spderrors = spderrors || !spdmatrixcheckinverse<Precision>(a, inva, isupper, n, threshold, info, rep);
                info = 0;
                unsetrep<Precision>(rep);
                matinv::spdmatrixcholeskyinverse<Precision>(invcha, n, isupper, info, rep);
                spderrors = spderrors || !spdmatrixcheckinverse<Precision>(a, invcha, isupper, n, threshold, info, rep);
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                // 2. test different methods
                //
                for(taskkind=0; taskkind<=2; taskkind++)
                {
                    unset2d<Precision>(a);
                    if( taskkind==0 )
                    {
                        
                        //
                        // all zeros
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 0;
                            }
                        }
                    }
                    if( taskkind==1 )
                    {
                        
                        //
                        // there is zero column
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        amp::vmul(a.getcolumn(k, 0, n-1), 0);
                    }
                    if( taskkind==2 )
                    {
                        
                        //
                        // there is zero row
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        amp::vmul(a.getrow(k, 0, n-1), 0);
                    }
                    info = 0;
                    unsetrep<Precision>(rep);
                    matinv::spdmatrixcholeskyinverse<Precision>(a, n, isupper, info, rep);
                    if( info!=-3 && info!=1 )
                    {
                        spderrors = true;
                    }
                    else
                    {
                        spderrors = spderrors || rep.r1<0 || rep.r1>1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                        spderrors = spderrors || rep.rinf<0 || rep.rinf>1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                    }
                }
            }
        }
    }


    /*************************************************************************
    HPD test
    *************************************************************************/
    template<unsigned int Precision>
    void testhpdinv(int maxn,
        int passcount,
        amp::ampf<Precision> threshold,
        bool& hpderrors)
    {
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > cha;
        ap::template_2d_array< amp::campf<Precision> > inva;
        ap::template_2d_array< amp::campf<Precision> > invcha;
        bool isupper;
        int i;
        int j;
        int k;
        int n;
        int pass;
        int taskkind;
        amp::campf<Precision> v;
        int info;
        matinv::matinvreport<Precision> rep;
        int i_;


        
        //
        // General square matrices:
        // * test general solvers
        // * test least squares solver
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(n=1; n<=maxn; n++)
            {
                isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
                
                //
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                matgen::hpdmatrixrndcond<Precision>(n, amp::ampf<Precision>(1000), a);
                cmatrixdrophalf<Precision>(a, n, isupper);
                cmatrixmakeacopy<Precision>(a, n, n, cha);
                if( !trfac::hpdmatrixcholesky<Precision>(cha, n, isupper) )
                {
                    continue;
                }
                cmatrixmakeacopy<Precision>(a, n, n, inva);
                cmatrixmakeacopy<Precision>(cha, n, n, invcha);
                info = 0;
                unsetrep<Precision>(rep);
                matinv::hpdmatrixinverse<Precision>(inva, n, isupper, info, rep);
                hpderrors = hpderrors || !hpdmatrixcheckinverse<Precision>(a, inva, isupper, n, threshold, info, rep);
                info = 0;
                unsetrep<Precision>(rep);
                matinv::hpdmatrixcholeskyinverse<Precision>(invcha, n, isupper, info, rep);
                hpderrors = hpderrors || !hpdmatrixcheckinverse<Precision>(a, invcha, isupper, n, threshold, info, rep);
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                // 2. test different methods
                //
                for(taskkind=0; taskkind<=2; taskkind++)
                {
                    cunset2d<Precision>(a);
                    if( taskkind==0 )
                    {
                        
                        //
                        // all zeros
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j) = 0;
                            }
                        }
                    }
                    if( taskkind==1 )
                    {
                        
                        //
                        // there is zero column
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            a(i_,k) = 0*a(i_,k);
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            a(k,i_) = 0*a(k,i_);
                        }
                    }
                    if( taskkind==2 )
                    {
                        
                        //
                        // there is zero row
                        //
                        a.setlength(n, n);
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                a(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                                a(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        for(i_=0; i_<=n-1;i_++)
                        {
                            a(k,i_) = 0*a(k,i_);
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            a(i_,k) = 0*a(i_,k);
                        }
                    }
                    info = 0;
                    unsetrep<Precision>(rep);
                    matinv::hpdmatrixcholeskyinverse<Precision>(a, n, isupper, info, rep);
                    if( info!=-3 && info!=1 )
                    {
                        hpderrors = true;
                    }
                    else
                    {
                        hpderrors = hpderrors || rep.r1<0 || rep.r1>1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                        hpderrors = hpderrors || rep.rinf<0 || rep.rinf>1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
                    }
                }
            }
        }
    }


    /*************************************************************************
    Unsets real matrix
    *************************************************************************/
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        x.setlength(1, 1);
        x(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets real vector
    *************************************************************************/
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        x.setlength(1);
        x(0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets real matrix
    *************************************************************************/
    template<unsigned int Precision>
    void cunset2d(ap::template_2d_array< amp::campf<Precision> >& x)
    {
        x.setlength(1, 1);
        x(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets real vector
    *************************************************************************/
    template<unsigned int Precision>
    void cunset1d(ap::template_1d_array< amp::campf<Precision> >& x)
    {
        x.setlength(1);
        x(0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets report
    *************************************************************************/
    template<unsigned int Precision>
    void unsetrep(matinv::matinvreport<Precision>& r)
    {
        r.r1 = -1;
        r.rinf = -1;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testmatinvunit_test_silent()
    {
        bool result;


        result = testmatinv<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testmatinvunit_test()
    {
        bool result;


        result = testmatinv<Precision>(false);
        return result;
    }
} // namespace

#endif
