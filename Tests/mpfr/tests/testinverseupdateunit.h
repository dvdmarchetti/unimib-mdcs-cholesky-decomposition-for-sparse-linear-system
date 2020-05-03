
#ifndef _testinverseupdateunit_h
#define _testinverseupdateunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "inverseupdate.h"
namespace testinverseupdateunit
{
    template<unsigned int Precision>
    bool testinverseupdate(bool silent);
    template<unsigned int Precision>
    void makeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& b);
    template<unsigned int Precision>
    void matlu(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots);
    template<unsigned int Precision>
    void generaterandomorthogonalmatrix(ap::template_2d_array< amp::ampf<Precision> >& a0,
        int n);
    template<unsigned int Precision>
    void generaterandommatrixcond(ap::template_2d_array< amp::ampf<Precision> >& a0,
        int n,
        amp::ampf<Precision> c);
    template<unsigned int Precision>
    bool invmattr(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunittriangular);
    template<unsigned int Precision>
    bool invmatlu(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n);
    template<unsigned int Precision>
    bool invmat(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> matrixdiff(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int n);
    template<unsigned int Precision>
    bool updandinv(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< amp::ampf<Precision> >& u,
        const ap::template_1d_array< amp::ampf<Precision> >& v,
        int n);
    template<unsigned int Precision>
    bool testinverseupdateunit_test_silent();
    template<unsigned int Precision>
    bool testinverseupdateunit_test();


    template<unsigned int Precision>
    bool testinverseupdate(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > inva;
        ap::template_2d_array< amp::ampf<Precision> > b1;
        ap::template_2d_array< amp::ampf<Precision> > b2;
        ap::template_1d_array< amp::ampf<Precision> > u;
        ap::template_1d_array< amp::ampf<Precision> > v;
        int n;
        int maxn;
        int i;
        int j;
        int updrow;
        int updcol;
        amp::ampf<Precision> val;
        int pass;
        int passcount;
        bool waserrors;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> c;


        waserrors = false;
        maxn = 10;
        passcount = 100;
        threshold = amp::ampf<Precision>("1.0E-6");
        
        //
        // process
        //
        for(n=1; n<=maxn; n++)
        {
            a.setbounds(0, n-1, 0, n-1);
            b1.setbounds(0, n-1, 0, n-1);
            b2.setbounds(0, n-1, 0, n-1);
            u.setbounds(0, n-1);
            v.setbounds(0, n-1);
            for(pass=1; pass<=passcount; pass++)
            {
                c = amp::exp<Precision>(amp::ampf<Precision>::getRandom()*amp::log<Precision>(amp::ampf<Precision>(10)));
                generaterandommatrixcond<Precision>(a, n, c);
                makeacopy<Precision>(a, n, n, inva);
                if( !invmat<Precision>(inva, n) )
                {
                    waserrors = true;
                    break;
                }
                
                //
                // Test simple update
                //
                updrow = ap::randominteger(n);
                updcol = ap::randominteger(n);
                val = amp::ampf<Precision>("0.1")*(2*amp::ampf<Precision>::getRandom()-1);
                for(i=0; i<=n-1; i++)
                {
                    if( i==updrow )
                    {
                        u(i) = val;
                    }
                    else
                    {
                        u(i) = 0;
                    }
                    if( i==updcol )
                    {
                        v(i) = 1;
                    }
                    else
                    {
                        v(i) = 0;
                    }
                }
                makeacopy<Precision>(a, n, n, b1);
                if( !updandinv<Precision>(b1, u, v, n) )
                {
                    waserrors = true;
                    break;
                }
                makeacopy<Precision>(inva, n, n, b2);
                inverseupdate::rmatrixinvupdatesimple<Precision>(b2, n, updrow, updcol, val);
                waserrors = waserrors || matrixdiff<Precision>(b1, b2, n, n)>threshold;
                
                //
                // Test row update
                //
                updrow = ap::randominteger(n);
                for(i=0; i<=n-1; i++)
                {
                    if( i==updrow )
                    {
                        u(i) = 1;
                    }
                    else
                    {
                        u(i) = 0;
                    }
                    v(i) = amp::ampf<Precision>("0.1")*(2*amp::ampf<Precision>::getRandom()-1);
                }
                makeacopy<Precision>(a, n, n, b1);
                if( !updandinv<Precision>(b1, u, v, n) )
                {
                    waserrors = true;
                    break;
                }
                makeacopy<Precision>(inva, n, n, b2);
                inverseupdate::rmatrixinvupdaterow<Precision>(b2, n, updrow, v);
                waserrors = waserrors || matrixdiff<Precision>(b1, b2, n, n)>threshold;
                
                //
                // Test column update
                //
                updcol = ap::randominteger(n);
                for(i=0; i<=n-1; i++)
                {
                    if( i==updcol )
                    {
                        v(i) = 1;
                    }
                    else
                    {
                        v(i) = 0;
                    }
                    u(i) = amp::ampf<Precision>("0.1")*(2*amp::ampf<Precision>::getRandom()-1);
                }
                makeacopy<Precision>(a, n, n, b1);
                if( !updandinv<Precision>(b1, u, v, n) )
                {
                    waserrors = true;
                    break;
                }
                makeacopy<Precision>(inva, n, n, b2);
                inverseupdate::rmatrixinvupdatecolumn<Precision>(b2, n, updcol, u);
                waserrors = waserrors || matrixdiff<Precision>(b1, b2, n, n)>threshold;
                
                //
                // Test full update
                //
                for(i=0; i<=n-1; i++)
                {
                    v(i) = amp::ampf<Precision>("0.1")*(2*amp::ampf<Precision>::getRandom()-1);
                    u(i) = amp::ampf<Precision>("0.1")*(2*amp::ampf<Precision>::getRandom()-1);
                }
                makeacopy<Precision>(a, n, n, b1);
                if( !updandinv<Precision>(b1, u, v, n) )
                {
                    waserrors = true;
                    break;
                }
                makeacopy<Precision>(inva, n, n, b2);
                inverseupdate::rmatrixinvupdateuv<Precision>(b2, n, u, v);
                waserrors = waserrors || matrixdiff<Precision>(b1, b2, n, n)>threshold;
            }
        }
        
        //
        // report
        //
        if( !silent )
        {
            printf("TESTING INVERSE UPDATE (REAL)\n");
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
    LU decomposition
    *************************************************************************/
    template<unsigned int Precision>
    void matlu(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots)
    {
        int i;
        int j;
        int jp;
        ap::template_1d_array< amp::ampf<Precision> > t1;
        amp::ampf<Precision> s;


        pivots.setbounds(0, ap::minint(m-1, n-1));
        t1.setbounds(0, ap::maxint(m-1, n-1));
        ap::ap_error::make_assertion(m>=0 && n>=0);
        
        //
        // Quick return if possible
        //
        if( m==0 || n==0 )
        {
            return;
        }
        for(j=0; j<=ap::minint(m-1, n-1); j++)
        {
            
            //
            // Find pivot and test for singularity.
            //
            jp = j;
            for(i=j+1; i<=m-1; i++)
            {
                if( amp::abs<Precision>(a(i,j))>amp::abs<Precision>(a(jp,j)) )
                {
                    jp = i;
                }
            }
            pivots(j) = jp;
            if( a(jp,j)!=0 )
            {
                
                //
                //Apply the interchange to rows
                //
                if( jp!=j )
                {
                    amp::vmove(t1.getvector(0, n-1), a.getrow(j, 0, n-1));
                    amp::vmove(a.getrow(j, 0, n-1), a.getrow(jp, 0, n-1));
                    amp::vmove(a.getrow(jp, 0, n-1), t1.getvector(0, n-1));
                }
                
                //
                //Compute elements J+1:M of J-th column.
                //
                if( j<m )
                {
                    jp = j+1;
                    s = 1/a(j,j);
                    amp::vmul(a.getcolumn(j, jp, m-1), s);
                }
            }
            if( j<ap::minint(m, n)-1 )
            {
                
                //
                //Update trailing submatrix.
                //
                jp = j+1;
                for(i=j+1; i<=m-1; i++)
                {
                    s = a(i,j);
                    amp::vsub(a.getrow(i, jp, n-1), a.getrow(j, jp, n-1), s);
                }
            }
        }
    }


    /*************************************************************************
    Generate matrix with given condition number C (2-norm)
    *************************************************************************/
    template<unsigned int Precision>
    void generaterandomorthogonalmatrix(ap::template_2d_array< amp::ampf<Precision> >& a0,
        int n)
    {
        amp::ampf<Precision> t;
        amp::ampf<Precision> lambda;
        int s;
        int i;
        int j;
        amp::ampf<Precision> u1;
        amp::ampf<Precision> u2;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > v;
        ap::template_2d_array< amp::ampf<Precision> > a;
        amp::ampf<Precision> sm;


        if( n<=0 )
        {
            return;
        }
        w.setbounds(1, n);
        v.setbounds(1, n);
        a.setbounds(1, n, 1, n);
        a0.setbounds(0, n-1, 0, n-1);
        
        //
        // Prepare A
        //
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                if( i==j )
                {
                    a(i,j) = 1;
                }
                else
                {
                    a(i,j) = 0;
                }
            }
        }
        
        //
        // Calculate A using Stewart algorithm
        //
        for(s=2; s<=n; s++)
        {
            
            //
            // Prepare v and Lambda = v'*v
            //
            do
            {
                i = 1;
                while( i<=s )
                {
                    u1 = 2*amp::ampf<Precision>::getRandom()-1;
                    u2 = 2*amp::ampf<Precision>::getRandom()-1;
                    sm = u1*u1+u2*u2;
                    if( sm==0 || sm>1 )
                    {
                        continue;
                    }
                    sm = amp::sqrt<Precision>(-2*amp::log<Precision>(sm)/sm);
                    v(i) = u1*sm;
                    if( i+1<=s )
                    {
                        v(i+1) = u2*sm;
                    }
                    i = i+2;
                }
                lambda = amp::vdotproduct(v.getvector(1, s), v.getvector(1, s));
            }
            while( lambda==0 );
            lambda = 2/lambda;
            
            //
            // A * (I - 2 vv'/v'v ) =
            //   = A - (2/v'v) * A * v * v' =
            //   = A - (2/v'v) * w * v'
            //  where w = Av
            //
            for(i=1; i<=s; i++)
            {
                t = amp::vdotproduct(a.getrow(i, 1, s), v.getvector(1, s));
                w(i) = t;
            }
            for(i=1; i<=s; i++)
            {
                t = w(i)*lambda;
                amp::vsub(a.getrow(i, 1, s), v.getvector(1, s), t);
            }
        }
        
        //
        //
        //
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                a0(i-1,j-1) = a(i,j);
            }
        }
    }


    template<unsigned int Precision>
    void generaterandommatrixcond(ap::template_2d_array< amp::ampf<Precision> >& a0,
        int n,
        amp::ampf<Precision> c)
    {
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;
        ap::template_2d_array< amp::ampf<Precision> > q1;
        ap::template_2d_array< amp::ampf<Precision> > q2;
        ap::template_1d_array< amp::ampf<Precision> > cc;
        int i;
        int j;
        int k;


        generaterandomorthogonalmatrix<Precision>(q1, n);
        generaterandomorthogonalmatrix<Precision>(q2, n);
        cc.setbounds(0, n-1);
        l1 = 0;
        l2 = amp::log<Precision>(1/c);
        cc(0) = amp::exp<Precision>(l1);
        for(i=1; i<=n-2; i++)
        {
            cc(i) = amp::exp<Precision>(amp::ampf<Precision>::getRandom()*(l2-l1)+l1);
        }
        cc(n-1) = amp::exp<Precision>(l2);
        a0.setbounds(0, n-1, 0, n-1);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a0(i,j) = 0;
                for(k=0; k<=n-1; k++)
                {
                    a0(i,j) = a0(i,j)+q1(i,k)*cc(k)*q2(j,k);
                }
            }
        }
    }


    /*************************************************************************
    triangular inverse
    *************************************************************************/
    template<unsigned int Precision>
    bool invmattr(ap::template_2d_array< amp::ampf<Precision> >& a,
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
    bool invmatlu(ap::template_2d_array< amp::ampf<Precision> >& a,
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
        if( !invmattr<Precision>(a, n, true, false) )
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
    bool invmat(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n)
    {
        bool result;
        ap::template_1d_array< int > pivots;


        matlu<Precision>(a, n, n, pivots);
        result = invmatlu<Precision>(a, pivots, n);
        return result;
    }


    /*************************************************************************
    Diff
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> matrixdiff(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int j;


        result = 0;
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                result = amp::maximum<Precision>(result, amp::abs<Precision>(b(i,j)-a(i,j)));
            }
        }
        return result;
    }


    /*************************************************************************
    Update and inverse
    *************************************************************************/
    template<unsigned int Precision>
    bool updandinv(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< amp::ampf<Precision> >& u,
        const ap::template_1d_array< amp::ampf<Precision> >& v,
        int n)
    {
        bool result;
        ap::template_1d_array< int > pivots;
        int i;
        amp::ampf<Precision> r;


        for(i=0; i<=n-1; i++)
        {
            r = u(i);
            amp::vadd(a.getrow(i, 0, n-1), v.getvector(0, n-1), r);
        }
        matlu<Precision>(a, n, n, pivots);
        result = invmatlu<Precision>(a, pivots, n);
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testinverseupdateunit_test_silent()
    {
        bool result;


        result = testinverseupdate<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testinverseupdateunit_test()
    {
        bool result;


        result = testinverseupdate<Precision>(false);
        return result;
    }
} // namespace

#endif
