
#ifndef _testrcondunit_h
#define _testrcondunit_h

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
namespace testrcondunit
{
    template<unsigned int Precision>
    bool testrcond(bool silent);
    template<unsigned int Precision>
    void rmatrixmakeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& b);
    template<unsigned int Precision>
    void rmatrixdrophalf(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool droplower);
    template<unsigned int Precision>
    void cmatrixdrophalf(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool droplower);
    template<unsigned int Precision>
    void rmatrixgenzero(ap::template_2d_array< amp::ampf<Precision> >& a0,
        int n);
    template<unsigned int Precision>
    bool rmatrixinvmattr(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunittriangular);
    template<unsigned int Precision>
    bool rmatrixinvmatlu(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n);
    template<unsigned int Precision>
    bool rmatrixinvmat(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    void rmatrixrefrcond(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision>& rc1,
        amp::ampf<Precision>& rcinf);
    template<unsigned int Precision>
    void cmatrixmakeacopy(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::campf<Precision> >& b);
    template<unsigned int Precision>
    void cmatrixgenzero(ap::template_2d_array< amp::campf<Precision> >& a0,
        int n);
    template<unsigned int Precision>
    bool cmatrixinvmattr(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunittriangular);
    template<unsigned int Precision>
    bool cmatrixinvmatlu(ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n);
    template<unsigned int Precision>
    bool cmatrixinvmat(ap::template_2d_array< amp::campf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    void cmatrixrefrcond(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        amp::ampf<Precision>& rc1,
        amp::ampf<Precision>& rcinf);
    template<unsigned int Precision>
    bool testrmatrixtrrcond(int maxn,
        int passcount);
    template<unsigned int Precision>
    bool testcmatrixtrrcond(int maxn,
        int passcount);
    template<unsigned int Precision>
    bool testrmatrixrcond(int maxn,
        int passcount);
    template<unsigned int Precision>
    bool testspdmatrixrcond(int maxn,
        int passcount);
    template<unsigned int Precision>
    bool testcmatrixrcond(int maxn,
        int passcount);
    template<unsigned int Precision>
    bool testhpdmatrixrcond(int maxn,
        int passcount);
    template<unsigned int Precision>
    bool testrcondunit_test_silent();
    template<unsigned int Precision>
    bool testrcondunit_test();


    template<unsigned int Precision>
    const amp::ampf<Precision>& threshold50()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("0.25");
        return v;
    }
    template<unsigned int Precision>
    const amp::ampf<Precision>& threshold90()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("0.10");
        return v;
    }


    template<unsigned int Precision>
    bool testrcond(bool silent)
    {
        bool result;
        int maxn;
        int passcount;
        bool waserrors;
        bool rtrerr;
        bool ctrerr;
        bool rerr;
        bool cerr;
        bool spderr;
        bool hpderr;


        maxn = 10;
        passcount = 100;
        
        //
        // report
        //
        rtrerr = !testrmatrixtrrcond<Precision>(maxn, passcount);
        ctrerr = !testcmatrixtrrcond<Precision>(maxn, passcount);
        rerr = !testrmatrixrcond<Precision>(maxn, passcount);
        cerr = !testcmatrixrcond<Precision>(maxn, passcount);
        spderr = !testspdmatrixrcond<Precision>(maxn, passcount);
        hpderr = !testhpdmatrixrcond<Precision>(maxn, passcount);
        waserrors = rtrerr || ctrerr || rerr || cerr || spderr || hpderr;
        if( !silent )
        {
            printf("TESTING RCOND\n");
            printf("REAL TRIANGULAR:                         ");
            if( !rtrerr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("COMPLEX TRIANGULAR:                      ");
            if( !ctrerr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("REAL:                                    ");
            if( !rerr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("SPD:                                     ");
            if( !spderr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("HPD:                                     ");
            if( !hpderr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("COMPLEX:                                 ");
            if( !cerr )
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
    Generate matrix with given condition number C (2-norm)
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixgenzero(ap::template_2d_array< amp::ampf<Precision> >& a0,
        int n)
    {
        int i;
        int j;


        a0.setlength(n, n);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a0(i,j) = 0;
            }
        }
    }


    /*************************************************************************
    triangular inverse
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixinvmattr(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunittriangular)
    {
        bool result;
        bool nounit;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> ajj;
        ap::template_1d_array< amp::ampf<Precision> > t;


        result = true;
        t.setbounds(0, n-1);
        
        //
        // Test the input parameters.
        //
        nounit = !isunittriangular;
        if( isupper )
        {
            
            //
            // Compute inverse of upper triangular matrix.
            //
            for(j=0; j<=n-1; j++)
            {
                if( nounit )
                {
                    if( a(j,j)==0 )
                    {
                        result = false;
                        return result;
                    }
                    a(j,j) = 1/a(j,j);
                    ajj = -a(j,j);
                }
                else
                {
                    ajj = -1;
                }
                
                //
                // Compute elements 1:j-1 of j-th column.
                //
                if( j>0 )
                {
                    amp::vmove(t.getvector(0, j-1), a.getcolumn(j, 0, j-1));
                    for(i=0; i<=j-1; i++)
                    {
                        if( i<j-1 )
                        {
                            v = amp::vdotproduct(a.getrow(i, i+1, j-1), t.getvector(i+1, j-1));
                        }
                        else
                        {
                            v = 0;
                        }
                        if( nounit )
                        {
                            a(i,j) = v+a(i,i)*t(i);
                        }
                        else
                        {
                            a(i,j) = v+t(i);
                        }
                    }
                    amp::vmul(a.getcolumn(j, 0, j-1), ajj);
                }
            }
        }
        else
        {
            
            //
            // Compute inverse of lower triangular matrix.
            //
            for(j=n-1; j>=0; j--)
            {
                if( nounit )
                {
                    if( a(j,j)==0 )
                    {
                        result = false;
                        return result;
                    }
                    a(j,j) = 1/a(j,j);
                    ajj = -a(j,j);
                }
                else
                {
                    ajj = -1;
                }
                if( j<n-1 )
                {
                    
                    //
                    // Compute elements j+1:n of j-th column.
                    //
                    amp::vmove(t.getvector(j+1, n-1), a.getcolumn(j, j+1, n-1));
                    for(i=j+1; i<=n-1; i++)
                    {
                        if( i>j+1 )
                        {
                            v = amp::vdotproduct(a.getrow(i, j+1, i-1), t.getvector(j+1, i-1));
                        }
                        else
                        {
                            v = 0;
                        }
                        if( nounit )
                        {
                            a(i,j) = v+a(i,i)*t(i);
                        }
                        else
                        {
                            a(i,j) = v+t(i);
                        }
                    }
                    amp::vmul(a.getcolumn(j, j+1, n-1), ajj);
                }
            }
        }
        return result;
    }


    /*************************************************************************
    LU inverse
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixinvmatlu(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > work;
        int i;
        int iws;
        int j;
        int jb;
        int jj;
        int jp;
        amp::ampf<Precision> v;


        result = true;
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            return result;
        }
        work.setbounds(0, n-1);
        
        //
        // Form inv(U)
        //
        if( !rmatrixinvmattr<Precision>(a, n, true, false) )
        {
            result = false;
            return result;
        }
        
        //
        // Solve the equation inv(A)*L = inv(U) for inv(A).
        //
        for(j=n-1; j>=0; j--)
        {
            
            //
            // Copy current column of L to WORK and replace with zeros.
            //
            for(i=j+1; i<=n-1; i++)
            {
                work(i) = a(i,j);
                a(i,j) = 0;
            }
            
            //
            // Compute current column of inv(A).
            //
            if( j<n-1 )
            {
                for(i=0; i<=n-1; i++)
                {
                    v = amp::vdotproduct(a.getrow(i, j+1, n-1), work.getvector(j+1, n-1));
                    a(i,j) = a(i,j)-v;
                }
            }
        }
        
        //
        // Apply column interchanges.
        //
        for(j=n-2; j>=0; j--)
        {
            jp = pivots(j);
            if( jp!=j )
            {
                amp::vmove(work.getvector(0, n-1), a.getcolumn(j, 0, n-1));
                amp::vmove(a.getcolumn(j, 0, n-1), a.getcolumn(jp, 0, n-1));
                amp::vmove(a.getcolumn(jp, 0, n-1), work.getvector(0, n-1));
            }
        }
        return result;
    }


    /*************************************************************************
    Matrix inverse
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixinvmat(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n)
    {
        bool result;
        ap::template_1d_array< int > pivots;


        trfac::rmatrixlu<Precision>(a, n, n, pivots);
        result = rmatrixinvmatlu<Precision>(a, pivots, n);
        return result;
    }


    /*************************************************************************
    reference RCond
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrefrcond(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision>& rc1,
        amp::ampf<Precision>& rcinf)
    {
        ap::template_2d_array< amp::ampf<Precision> > inva;
        amp::ampf<Precision> nrm1a;
        amp::ampf<Precision> nrminfa;
        amp::ampf<Precision> nrm1inva;
        amp::ampf<Precision> nrminfinva;
        amp::ampf<Precision> v;
        int k;
        int i;


        
        //
        // inv A
        //
        rmatrixmakeacopy<Precision>(a, n, n, inva);
        if( !rmatrixinvmat<Precision>(inva, n) )
        {
            rc1 = 0;
            rcinf = 0;
            return;
        }
        
        //
        // norm A
        //
        nrm1a = 0;
        nrminfa = 0;
        for(k=0; k<=n-1; k++)
        {
            v = 0;
            for(i=0; i<=n-1; i++)
            {
                v = v+amp::abs<Precision>(a(i,k));
            }
            nrm1a = amp::maximum<Precision>(nrm1a, v);
            v = 0;
            for(i=0; i<=n-1; i++)
            {
                v = v+amp::abs<Precision>(a(k,i));
            }
            nrminfa = amp::maximum<Precision>(nrminfa, v);
        }
        
        //
        // norm inv A
        //
        nrm1inva = 0;
        nrminfinva = 0;
        for(k=0; k<=n-1; k++)
        {
            v = 0;
            for(i=0; i<=n-1; i++)
            {
                v = v+amp::abs<Precision>(inva(i,k));
            }
            nrm1inva = amp::maximum<Precision>(nrm1inva, v);
            v = 0;
            for(i=0; i<=n-1; i++)
            {
                v = v+amp::abs<Precision>(inva(k,i));
            }
            nrminfinva = amp::maximum<Precision>(nrminfinva, v);
        }
        
        //
        // result
        //
        rc1 = nrm1inva*nrm1a;
        rcinf = nrminfinva*nrminfa;
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
    Generate matrix with given condition number C (2-norm)
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixgenzero(ap::template_2d_array< amp::campf<Precision> >& a0,
        int n)
    {
        int i;
        int j;


        a0.setlength(n, n);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a0(i,j) = 0;
            }
        }
    }


    /*************************************************************************
    triangular inverse
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixinvmattr(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunittriangular)
    {
        bool result;
        bool nounit;
        int i;
        int j;
        amp::campf<Precision> v;
        amp::campf<Precision> ajj;
        ap::template_1d_array< amp::campf<Precision> > t;
        int i_;


        result = true;
        t.setbounds(0, n-1);
        
        //
        // Test the input parameters.
        //
        nounit = !isunittriangular;
        if( isupper )
        {
            
            //
            // Compute inverse of upper triangular matrix.
            //
            for(j=0; j<=n-1; j++)
            {
                if( nounit )
                {
                    if( a(j,j)==0 )
                    {
                        result = false;
                        return result;
                    }
                    a(j,j) = 1/a(j,j);
                    ajj = -a(j,j);
                }
                else
                {
                    ajj = -1;
                }
                
                //
                // Compute elements 1:j-1 of j-th column.
                //
                if( j>0 )
                {
                    for(i_=0; i_<=j-1;i_++)
                    {
                        t(i_) = a(i_,j);
                    }
                    for(i=0; i<=j-1; i++)
                    {
                        if( i<j-1 )
                        {
                            v = 0.0;
                            for(i_=i+1; i_<=j-1;i_++)
                            {
                                v += a(i,i_)*t(i_);
                            }
                        }
                        else
                        {
                            v = 0;
                        }
                        if( nounit )
                        {
                            a(i,j) = v+a(i,i)*t(i);
                        }
                        else
                        {
                            a(i,j) = v+t(i);
                        }
                    }
                    for(i_=0; i_<=j-1;i_++)
                    {
                        a(i_,j) = ajj*a(i_,j);
                    }
                }
            }
        }
        else
        {
            
            //
            // Compute inverse of lower triangular matrix.
            //
            for(j=n-1; j>=0; j--)
            {
                if( nounit )
                {
                    if( a(j,j)==0 )
                    {
                        result = false;
                        return result;
                    }
                    a(j,j) = 1/a(j,j);
                    ajj = -a(j,j);
                }
                else
                {
                    ajj = -1;
                }
                if( j<n-1 )
                {
                    
                    //
                    // Compute elements j+1:n of j-th column.
                    //
                    for(i_=j+1; i_<=n-1;i_++)
                    {
                        t(i_) = a(i_,j);
                    }
                    for(i=j+1; i<=n-1; i++)
                    {
                        if( i>j+1 )
                        {
                            v = 0.0;
                            for(i_=j+1; i_<=i-1;i_++)
                            {
                                v += a(i,i_)*t(i_);
                            }
                        }
                        else
                        {
                            v = 0;
                        }
                        if( nounit )
                        {
                            a(i,j) = v+a(i,i)*t(i);
                        }
                        else
                        {
                            a(i,j) = v+t(i);
                        }
                    }
                    for(i_=j+1; i_<=n-1;i_++)
                    {
                        a(i_,j) = ajj*a(i_,j);
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    LU inverse
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixinvmatlu(ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n)
    {
        bool result;
        ap::template_1d_array< amp::campf<Precision> > work;
        int i;
        int iws;
        int j;
        int jb;
        int jj;
        int jp;
        amp::campf<Precision> v;
        int i_;


        result = true;
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            return result;
        }
        work.setbounds(0, n-1);
        
        //
        // Form inv(U)
        //
        if( !cmatrixinvmattr<Precision>(a, n, true, false) )
        {
            result = false;
            return result;
        }
        
        //
        // Solve the equation inv(A)*L = inv(U) for inv(A).
        //
        for(j=n-1; j>=0; j--)
        {
            
            //
            // Copy current column of L to WORK and replace with zeros.
            //
            for(i=j+1; i<=n-1; i++)
            {
                work(i) = a(i,j);
                a(i,j) = 0;
            }
            
            //
            // Compute current column of inv(A).
            //
            if( j<n-1 )
            {
                for(i=0; i<=n-1; i++)
                {
                    v = 0.0;
                    for(i_=j+1; i_<=n-1;i_++)
                    {
                        v += a(i,i_)*work(i_);
                    }
                    a(i,j) = a(i,j)-v;
                }
            }
        }
        
        //
        // Apply column interchanges.
        //
        for(j=n-2; j>=0; j--)
        {
            jp = pivots(j);
            if( jp!=j )
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    work(i_) = a(i_,j);
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    a(i_,j) = a(i_,jp);
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    a(i_,jp) = work(i_);
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Matrix inverse
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixinvmat(ap::template_2d_array< amp::campf<Precision> >& a,
        int n)
    {
        bool result;
        ap::template_1d_array< int > pivots;


        trfac::cmatrixlu<Precision>(a, n, n, pivots);
        result = cmatrixinvmatlu<Precision>(a, pivots, n);
        return result;
    }


    /*************************************************************************
    reference RCond
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrefrcond(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        amp::ampf<Precision>& rc1,
        amp::ampf<Precision>& rcinf)
    {
        ap::template_2d_array< amp::campf<Precision> > inva;
        amp::ampf<Precision> nrm1a;
        amp::ampf<Precision> nrminfa;
        amp::ampf<Precision> nrm1inva;
        amp::ampf<Precision> nrminfinva;
        amp::ampf<Precision> v;
        int k;
        int i;


        
        //
        // inv A
        //
        cmatrixmakeacopy<Precision>(a, n, n, inva);
        if( !cmatrixinvmat<Precision>(inva, n) )
        {
            rc1 = 0;
            rcinf = 0;
            return;
        }
        
        //
        // norm A
        //
        nrm1a = 0;
        nrminfa = 0;
        for(k=0; k<=n-1; k++)
        {
            v = 0;
            for(i=0; i<=n-1; i++)
            {
                v = v+amp::abscomplex<Precision>(a(i,k));
            }
            nrm1a = amp::maximum<Precision>(nrm1a, v);
            v = 0;
            for(i=0; i<=n-1; i++)
            {
                v = v+amp::abscomplex<Precision>(a(k,i));
            }
            nrminfa = amp::maximum<Precision>(nrminfa, v);
        }
        
        //
        // norm inv A
        //
        nrm1inva = 0;
        nrminfinva = 0;
        for(k=0; k<=n-1; k++)
        {
            v = 0;
            for(i=0; i<=n-1; i++)
            {
                v = v+amp::abscomplex<Precision>(inva(i,k));
            }
            nrm1inva = amp::maximum<Precision>(nrm1inva, v);
            v = 0;
            for(i=0; i<=n-1; i++)
            {
                v = v+amp::abscomplex<Precision>(inva(k,i));
            }
            nrminfinva = amp::maximum<Precision>(nrminfinva, v);
        }
        
        //
        // result
        //
        rc1 = nrm1inva*nrm1a;
        rcinf = nrminfinva*nrminfa;
    }


    /*************************************************************************
    Returns True for successful test, False - for failed test
    *************************************************************************/
    template<unsigned int Precision>
    bool testrmatrixtrrcond(int maxn,
        int passcount)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > ea;
        ap::template_1d_array< int > p;
        int n;
        int i;
        int j;
        int j1;
        int j2;
        int pass;
        bool err50;
        bool err90;
        bool errspec;
        bool errless;
        amp::ampf<Precision> erc1;
        amp::ampf<Precision> ercinf;
        ap::template_1d_array< amp::ampf<Precision> > q50;
        ap::template_1d_array< amp::ampf<Precision> > q90;
        amp::ampf<Precision> v;
        bool isupper;
        bool isunit;


        err50 = false;
        err90 = false;
        errless = false;
        errspec = false;
        q50.setlength(2);
        q90.setlength(2);
        for(n=1; n<=maxn; n++)
        {
            
            //
            // special test for zero matrix
            //
            rmatrixgenzero<Precision>(a, n);
            errspec = errspec || rcond::rmatrixtrrcond1<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
            errspec = errspec || rcond::rmatrixtrrcondinf<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
            
            //
            // general test
            //
            a.setlength(n, n);
            for(i=0; i<=1; i++)
            {
                q50(i) = 0;
                q90(i) = 0;
            }
            for(pass=1; pass<=passcount; pass++)
            {
                isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
                isunit = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    a(i,i) = 1+amp::ampf<Precision>::getRandom();
                }
                rmatrixmakeacopy<Precision>(a, n, n, ea);
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
                        ea(i,j) = 0;
                    }
                    if( isunit )
                    {
                        ea(i,i) = 1;
                    }
                }
                rmatrixrefrcond<Precision>(ea, n, erc1, ercinf);
                
                //
                // 1-norm
                //
                v = 1/rcond::rmatrixtrrcond1<Precision>(a, n, isupper, isunit);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(0) = q50(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(0) = q90(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
                
                //
                // Inf-norm
                //
                v = 1/rcond::rmatrixtrrcondinf<Precision>(a, n, isupper, isunit);
                if( v>=threshold50<Precision>()*ercinf )
                {
                    q50(1) = q50(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*ercinf )
                {
                    q90(1) = q90(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>ercinf*amp::ampf<Precision>("1.001");
            }
            for(i=0; i<=1; i++)
            {
                err50 = err50 || q50(i)<amp::ampf<Precision>("0.50");
                err90 = err90 || q90(i)<amp::ampf<Precision>("0.90");
            }
            
            //
            // degenerate matrix test
            //
            if( n>=3 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                a(0,0) = 1;
                a(n-1,n-1) = 1;
                errspec = errspec || rcond::rmatrixtrrcond1<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
                errspec = errspec || rcond::rmatrixtrrcondinf<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
            }
            
            //
            // near-degenerate matrix test
            //
            if( n>=2 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    a(i,i) = 1;
                }
                i = ap::randominteger(n);
                a(i,i) = amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
                errspec = errspec || rcond::rmatrixtrrcond1<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
                errspec = errspec || rcond::rmatrixtrrcondinf<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
            }
        }
        
        //
        // report
        //
        result = !(err50 || err90 || errless || errspec);
        return result;
    }


    /*************************************************************************
    Returns True for successful test, False - for failed test
    *************************************************************************/
    template<unsigned int Precision>
    bool testcmatrixtrrcond(int maxn,
        int passcount)
    {
        bool result;
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > ea;
        ap::template_1d_array< int > p;
        int n;
        int i;
        int j;
        int j1;
        int j2;
        int pass;
        bool err50;
        bool err90;
        bool errspec;
        bool errless;
        amp::ampf<Precision> erc1;
        amp::ampf<Precision> ercinf;
        ap::template_1d_array< amp::ampf<Precision> > q50;
        ap::template_1d_array< amp::ampf<Precision> > q90;
        amp::ampf<Precision> v;
        bool isupper;
        bool isunit;


        err50 = false;
        err90 = false;
        errless = false;
        errspec = false;
        q50.setlength(2);
        q90.setlength(2);
        for(n=1; n<=maxn; n++)
        {
            
            //
            // special test for zero matrix
            //
            cmatrixgenzero<Precision>(a, n);
            errspec = errspec || rcond::cmatrixtrrcond1<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
            errspec = errspec || rcond::cmatrixtrrcondinf<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
            
            //
            // general test
            //
            a.setlength(n, n);
            for(i=0; i<=1; i++)
            {
                q50(i) = 0;
                q90(i) = 0;
            }
            for(pass=1; pass<=passcount; pass++)
            {
                isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
                isunit = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j).x = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
                        a(i,j).y = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    a(i,i).x = 1+amp::ampf<Precision>::getRandom();
                    a(i,i).y = 1+amp::ampf<Precision>::getRandom();
                }
                cmatrixmakeacopy<Precision>(a, n, n, ea);
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
                        ea(i,j) = 0;
                    }
                    if( isunit )
                    {
                        ea(i,i) = 1;
                    }
                }
                cmatrixrefrcond<Precision>(ea, n, erc1, ercinf);
                
                //
                // 1-norm
                //
                v = 1/rcond::cmatrixtrrcond1<Precision>(a, n, isupper, isunit);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(0) = q50(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(0) = q90(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
                
                //
                // Inf-norm
                //
                v = 1/rcond::cmatrixtrrcondinf<Precision>(a, n, isupper, isunit);
                if( v>=threshold50<Precision>()*ercinf )
                {
                    q50(1) = q50(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*ercinf )
                {
                    q90(1) = q90(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>ercinf*amp::ampf<Precision>("1.001");
            }
            for(i=0; i<=1; i++)
            {
                err50 = err50 || q50(i)<amp::ampf<Precision>("0.50");
                err90 = err90 || q90(i)<amp::ampf<Precision>("0.90");
            }
            
            //
            // degenerate matrix test
            //
            if( n>=3 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                a(0,0) = 1;
                a(n-1,n-1) = 1;
                errspec = errspec || rcond::cmatrixtrrcond1<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
                errspec = errspec || rcond::cmatrixtrrcondinf<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
            }
            
            //
            // near-degenerate matrix test
            //
            if( n>=2 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    a(i,i) = 1;
                }
                i = ap::randominteger(n);
                a(i,i) = amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
                errspec = errspec || rcond::cmatrixtrrcond1<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
                errspec = errspec || rcond::cmatrixtrrcondinf<Precision>(a, n, amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5"), false)!=0;
            }
        }
        
        //
        // report
        //
        result = !(err50 || err90 || errless || errspec);
        return result;
    }


    /*************************************************************************
    Returns True for successful test, False - for failed test
    *************************************************************************/
    template<unsigned int Precision>
    bool testrmatrixrcond(int maxn,
        int passcount)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > lua;
        ap::template_1d_array< int > p;
        int n;
        int i;
        int j;
        int pass;
        bool err50;
        bool err90;
        bool errspec;
        bool errless;
        amp::ampf<Precision> erc1;
        amp::ampf<Precision> ercinf;
        ap::template_1d_array< amp::ampf<Precision> > q50;
        ap::template_1d_array< amp::ampf<Precision> > q90;
        amp::ampf<Precision> v;


        err50 = false;
        err90 = false;
        errless = false;
        errspec = false;
        q50.setbounds(0, 3);
        q90.setbounds(0, 3);
        for(n=1; n<=maxn; n++)
        {
            
            //
            // special test for zero matrix
            //
            rmatrixgenzero<Precision>(a, n);
            rmatrixmakeacopy<Precision>(a, n, n, lua);
            trfac::rmatrixlu<Precision>(lua, n, n, p);
            errspec = errspec || rcond::rmatrixrcond1<Precision>(a, n)!=0;
            errspec = errspec || rcond::rmatrixrcondinf<Precision>(a, n)!=0;
            errspec = errspec || rcond::rmatrixlurcond1<Precision>(lua, n)!=0;
            errspec = errspec || rcond::rmatrixlurcondinf<Precision>(lua, n)!=0;
            
            //
            // general test
            //
            a.setbounds(0, n-1, 0, n-1);
            for(i=0; i<=3; i++)
            {
                q50(i) = 0;
                q90(i) = 0;
            }
            for(pass=1; pass<=passcount; pass++)
            {
                matgen::rmatrixrndcond<Precision>(n, amp::exp<Precision>(amp::ampf<Precision>::getRandom()*amp::log<Precision>(amp::ampf<Precision>(1000))), a);
                rmatrixmakeacopy<Precision>(a, n, n, lua);
                trfac::rmatrixlu<Precision>(lua, n, n, p);
                rmatrixrefrcond<Precision>(a, n, erc1, ercinf);
                
                //
                // 1-norm, normal
                //
                v = 1/rcond::rmatrixrcond1<Precision>(a, n);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(0) = q50(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(0) = q90(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
                
                //
                // 1-norm, LU
                //
                v = 1/rcond::rmatrixlurcond1<Precision>(lua, n);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(1) = q50(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(1) = q90(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
                
                //
                // Inf-norm, normal
                //
                v = 1/rcond::rmatrixrcondinf<Precision>(a, n);
                if( v>=threshold50<Precision>()*ercinf )
                {
                    q50(2) = q50(2)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*ercinf )
                {
                    q90(2) = q90(2)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>ercinf*amp::ampf<Precision>("1.001");
                
                //
                // Inf-norm, LU
                //
                v = 1/rcond::rmatrixlurcondinf<Precision>(lua, n);
                if( v>=threshold50<Precision>()*ercinf )
                {
                    q50(3) = q50(3)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*ercinf )
                {
                    q90(3) = q90(3)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>ercinf*amp::ampf<Precision>("1.001");
            }
            for(i=0; i<=3; i++)
            {
                err50 = err50 || q50(i)<amp::ampf<Precision>("0.50");
                err90 = err90 || q90(i)<amp::ampf<Precision>("0.90");
            }
            
            //
            // degenerate matrix test
            //
            if( n>=3 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                a(0,0) = 1;
                a(n-1,n-1) = 1;
                errspec = errspec || rcond::rmatrixrcond1<Precision>(a, n)!=0;
                errspec = errspec || rcond::rmatrixrcondinf<Precision>(a, n)!=0;
                errspec = errspec || rcond::rmatrixlurcond1<Precision>(a, n)!=0;
                errspec = errspec || rcond::rmatrixlurcondinf<Precision>(a, n)!=0;
            }
            
            //
            // near-degenerate matrix test
            //
            if( n>=2 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    a(i,i) = 1;
                }
                i = ap::randominteger(n);
                a(i,i) = amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
                errspec = errspec || rcond::rmatrixrcond1<Precision>(a, n)!=0;
                errspec = errspec || rcond::rmatrixrcondinf<Precision>(a, n)!=0;
                errspec = errspec || rcond::rmatrixlurcond1<Precision>(a, n)!=0;
                errspec = errspec || rcond::rmatrixlurcondinf<Precision>(a, n)!=0;
            }
        }
        
        //
        // report
        //
        result = !(err50 || err90 || errless || errspec);
        return result;
    }


    /*************************************************************************
    Returns True for successful test, False - for failed test
    *************************************************************************/
    template<unsigned int Precision>
    bool testspdmatrixrcond(int maxn,
        int passcount)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > cha;
        ap::template_1d_array< int > p;
        int n;
        int i;
        int j;
        int pass;
        bool err50;
        bool err90;
        bool errspec;
        bool errless;
        bool isupper;
        amp::ampf<Precision> erc1;
        amp::ampf<Precision> ercinf;
        ap::template_1d_array< amp::ampf<Precision> > q50;
        ap::template_1d_array< amp::ampf<Precision> > q90;
        amp::ampf<Precision> v;


        err50 = false;
        err90 = false;
        errless = false;
        errspec = false;
        q50.setlength(2);
        q90.setlength(2);
        for(n=1; n<=maxn; n++)
        {
            isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
            
            //
            // general test
            //
            a.setlength(n, n);
            for(i=0; i<=1; i++)
            {
                q50(i) = 0;
                q90(i) = 0;
            }
            for(pass=1; pass<=passcount; pass++)
            {
                matgen::spdmatrixrndcond<Precision>(n, amp::exp<Precision>(amp::ampf<Precision>::getRandom()*amp::log<Precision>(amp::ampf<Precision>(1000))), a);
                rmatrixrefrcond<Precision>(a, n, erc1, ercinf);
                rmatrixdrophalf<Precision>(a, n, isupper);
                rmatrixmakeacopy<Precision>(a, n, n, cha);
                trfac::spdmatrixcholesky<Precision>(cha, n, isupper);
                
                //
                // normal
                //
                v = 1/rcond::spdmatrixrcond<Precision>(a, n, isupper);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(0) = q50(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(0) = q90(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
                
                //
                // Cholesky
                //
                v = 1/rcond::spdmatrixcholeskyrcond<Precision>(cha, n, isupper);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(1) = q50(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(1) = q90(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
            }
            for(i=0; i<=1; i++)
            {
                err50 = err50 || q50(i)<amp::ampf<Precision>("0.50");
                err90 = err90 || q90(i)<amp::ampf<Precision>("0.90");
            }
            
            //
            // degenerate matrix test
            //
            if( n>=3 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                a(0,0) = 1;
                a(n-1,n-1) = 1;
                errspec = errspec || rcond::spdmatrixrcond<Precision>(a, n, isupper)!=-1;
                errspec = errspec || rcond::spdmatrixcholeskyrcond<Precision>(a, n, isupper)!=0;
            }
            
            //
            // near-degenerate matrix test
            //
            if( n>=2 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    a(i,i) = 1;
                }
                i = ap::randominteger(n);
                a(i,i) = amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
                errspec = errspec || rcond::spdmatrixrcond<Precision>(a, n, isupper)!=0;
                errspec = errspec || rcond::spdmatrixcholeskyrcond<Precision>(a, n, isupper)!=0;
            }
        }
        
        //
        // report
        //
        result = !(err50 || err90 || errless || errspec);
        return result;
    }


    /*************************************************************************
    Returns True for successful test, False - for failed test
    *************************************************************************/
    template<unsigned int Precision>
    bool testcmatrixrcond(int maxn,
        int passcount)
    {
        bool result;
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > lua;
        ap::template_1d_array< int > p;
        int n;
        int i;
        int j;
        int pass;
        bool err50;
        bool err90;
        bool errless;
        bool errspec;
        amp::ampf<Precision> erc1;
        amp::ampf<Precision> ercinf;
        ap::template_1d_array< amp::ampf<Precision> > q50;
        ap::template_1d_array< amp::ampf<Precision> > q90;
        amp::ampf<Precision> v;


        q50.setbounds(0, 3);
        q90.setbounds(0, 3);
        err50 = false;
        err90 = false;
        errless = false;
        errspec = false;
        
        //
        // process
        //
        for(n=1; n<=maxn; n++)
        {
            
            //
            // special test for zero matrix
            //
            cmatrixgenzero<Precision>(a, n);
            cmatrixmakeacopy<Precision>(a, n, n, lua);
            trfac::cmatrixlu<Precision>(lua, n, n, p);
            errspec = errspec || rcond::cmatrixrcond1<Precision>(a, n)!=0;
            errspec = errspec || rcond::cmatrixrcondinf<Precision>(a, n)!=0;
            errspec = errspec || rcond::cmatrixlurcond1<Precision>(lua, n)!=0;
            errspec = errspec || rcond::cmatrixlurcondinf<Precision>(lua, n)!=0;
            
            //
            // general test
            //
            a.setbounds(0, n-1, 0, n-1);
            for(i=0; i<=3; i++)
            {
                q50(i) = 0;
                q90(i) = 0;
            }
            for(pass=1; pass<=passcount; pass++)
            {
                matgen::cmatrixrndcond<Precision>(n, amp::exp<Precision>(amp::ampf<Precision>::getRandom()*amp::log<Precision>(amp::ampf<Precision>(1000))), a);
                cmatrixmakeacopy<Precision>(a, n, n, lua);
                trfac::cmatrixlu<Precision>(lua, n, n, p);
                cmatrixrefrcond<Precision>(a, n, erc1, ercinf);
                
                //
                // 1-norm, normal
                //
                v = 1/rcond::cmatrixrcond1<Precision>(a, n);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(0) = q50(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(0) = q90(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
                
                //
                // 1-norm, LU
                //
                v = 1/rcond::cmatrixlurcond1<Precision>(lua, n);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(1) = q50(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(1) = q90(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
                
                //
                // Inf-norm, normal
                //
                v = 1/rcond::cmatrixrcondinf<Precision>(a, n);
                if( v>=threshold50<Precision>()*ercinf )
                {
                    q50(2) = q50(2)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*ercinf )
                {
                    q90(2) = q90(2)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>ercinf*amp::ampf<Precision>("1.001");
                
                //
                // Inf-norm, LU
                //
                v = 1/rcond::cmatrixlurcondinf<Precision>(lua, n);
                if( v>=threshold50<Precision>()*ercinf )
                {
                    q50(3) = q50(3)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*ercinf )
                {
                    q90(3) = q90(3)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>ercinf*amp::ampf<Precision>("1.001");
            }
            for(i=0; i<=3; i++)
            {
                err50 = err50 || q50(i)<amp::ampf<Precision>("0.50");
                err90 = err90 || q90(i)<amp::ampf<Precision>("0.90");
            }
            
            //
            // degenerate matrix test
            //
            if( n>=3 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                a(0,0) = 1;
                a(n-1,n-1) = 1;
                errspec = errspec || rcond::cmatrixrcond1<Precision>(a, n)!=0;
                errspec = errspec || rcond::cmatrixrcondinf<Precision>(a, n)!=0;
                errspec = errspec || rcond::cmatrixlurcond1<Precision>(a, n)!=0;
                errspec = errspec || rcond::cmatrixlurcondinf<Precision>(a, n)!=0;
            }
            
            //
            // near-degenerate matrix test
            //
            if( n>=2 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    a(i,i) = 1;
                }
                i = ap::randominteger(n);
                a(i,i) = amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
                errspec = errspec || rcond::cmatrixrcond1<Precision>(a, n)!=0;
                errspec = errspec || rcond::cmatrixrcondinf<Precision>(a, n)!=0;
                errspec = errspec || rcond::cmatrixlurcond1<Precision>(a, n)!=0;
                errspec = errspec || rcond::cmatrixlurcondinf<Precision>(a, n)!=0;
            }
        }
        
        //
        // report
        //
        result = !(err50 || err90 || errless || errspec);
        return result;
    }


    /*************************************************************************
    Returns True for successful test, False - for failed test
    *************************************************************************/
    template<unsigned int Precision>
    bool testhpdmatrixrcond(int maxn,
        int passcount)
    {
        bool result;
        ap::template_2d_array< amp::campf<Precision> > a;
        ap::template_2d_array< amp::campf<Precision> > cha;
        ap::template_1d_array< int > p;
        int n;
        int i;
        int j;
        int pass;
        bool err50;
        bool err90;
        bool errspec;
        bool errless;
        bool isupper;
        amp::ampf<Precision> erc1;
        amp::ampf<Precision> ercinf;
        ap::template_1d_array< amp::ampf<Precision> > q50;
        ap::template_1d_array< amp::ampf<Precision> > q90;
        amp::ampf<Precision> v;


        err50 = false;
        err90 = false;
        errless = false;
        errspec = false;
        q50.setlength(2);
        q90.setlength(2);
        for(n=1; n<=maxn; n++)
        {
            isupper = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
            
            //
            // general test
            //
            a.setlength(n, n);
            for(i=0; i<=1; i++)
            {
                q50(i) = 0;
                q90(i) = 0;
            }
            for(pass=1; pass<=passcount; pass++)
            {
                matgen::hpdmatrixrndcond<Precision>(n, amp::exp<Precision>(amp::ampf<Precision>::getRandom()*amp::log<Precision>(amp::ampf<Precision>(1000))), a);
                cmatrixrefrcond<Precision>(a, n, erc1, ercinf);
                cmatrixdrophalf<Precision>(a, n, isupper);
                cmatrixmakeacopy<Precision>(a, n, n, cha);
                trfac::hpdmatrixcholesky<Precision>(cha, n, isupper);
                
                //
                // normal
                //
                v = 1/rcond::hpdmatrixrcond<Precision>(a, n, isupper);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(0) = q50(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(0) = q90(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
                
                //
                // Cholesky
                //
                v = 1/rcond::hpdmatrixcholeskyrcond<Precision>(cha, n, isupper);
                if( v>=threshold50<Precision>()*erc1 )
                {
                    q50(1) = q50(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                if( v>=threshold90<Precision>()*erc1 )
                {
                    q90(1) = q90(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                }
                errless = errless || v>erc1*amp::ampf<Precision>("1.001");
            }
            for(i=0; i<=1; i++)
            {
                err50 = err50 || q50(i)<amp::ampf<Precision>("0.50");
                err90 = err90 || q90(i)<amp::ampf<Precision>("0.90");
            }
            
            //
            // degenerate matrix test
            //
            if( n>=3 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                a(0,0) = 1;
                a(n-1,n-1) = 1;
                errspec = errspec || rcond::hpdmatrixrcond<Precision>(a, n, isupper)!=-1;
                errspec = errspec || rcond::hpdmatrixcholeskyrcond<Precision>(a, n, isupper)!=0;
            }
            
            //
            // near-degenerate matrix test
            //
            if( n>=2 )
            {
                a.setlength(n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = amp::ampf<Precision>("0.0");
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    a(i,i) = 1;
                }
                i = ap::randominteger(n);
                a(i,i) = amp::ampf<Precision>("0.1")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
                errspec = errspec || rcond::hpdmatrixrcond<Precision>(a, n, isupper)!=0;
                errspec = errspec || rcond::hpdmatrixcholeskyrcond<Precision>(a, n, isupper)!=0;
            }
        }
        
        //
        // report
        //
        result = !(err50 || err90 || errless || errspec);
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testrcondunit_test_silent()
    {
        bool result;


        result = testrcond<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testrcondunit_test()
    {
        bool result;


        result = testrcond<Precision>(false);
        return result;
    }
} // namespace

#endif
