
#ifndef _testrcondldltunit_h
#define _testrcondldltunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "ldlt.h"
#include "ssolve.h"
#include "estnorm.h"
#include "srcond.h"
namespace testrcondldltunit
{
    template<unsigned int Precision>
    bool testrcondldlt(bool silent);
    template<unsigned int Precision>
    void generatematrix(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        int task);
    template<unsigned int Precision>
    void makeacopy(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_2d_array< amp::ampf<Precision> >& b);
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
    void matlu(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots);
    template<unsigned int Precision>
    bool invmat(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    void refrcond(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        amp::ampf<Precision>& rc1,
        amp::ampf<Precision>& rcinf);
    template<unsigned int Precision>
    bool testrcondldltunit_test_silent();
    template<unsigned int Precision>
    bool testrcondldltunit_test();


    template<unsigned int Precision>
    bool testrcondldlt(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        ap::template_2d_array< amp::ampf<Precision> > a3;
        int minij;
        ap::template_1d_array< int > p;
        int n;
        int maxn;
        int i;
        int j;
        int pass;
        int passcount;
        bool waserrors;
        bool err50;
        bool err90;
        bool errless;
        amp::ampf<Precision> erc1;
        amp::ampf<Precision> ercinf;
        ap::template_1d_array< amp::ampf<Precision> > q50;
        ap::template_1d_array< amp::ampf<Precision> > q90;
        amp::ampf<Precision> v;
        amp::ampf<Precision> threshold50;
        amp::ampf<Precision> threshold90;
        int htask;
        int mtask;
        bool upperin;


        waserrors = false;
        err50 = false;
        err90 = false;
        errless = false;
        maxn = 10;
        passcount = 100;
        threshold50 = amp::ampf<Precision>("0.5");
        threshold90 = amp::ampf<Precision>("0.1");
        q50.setbounds(0, 1);
        q90.setbounds(0, 1);
        
        //
        // process
        //
        for(n=1; n<=maxn; n++)
        {
            a.setbounds(0, n-1, 0, n-1);
            a2.setbounds(0, n-1, 0, n-1);
            a3.setbounds(0, n-1, 0, n-1);
            for(htask=0; htask<=1; htask++)
            {
                for(mtask=2; mtask<=2; mtask++)
                {
                    for(i=0; i<=1; i++)
                    {
                        q50(i) = 0;
                        q90(i) = 0;
                    }
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
                        refrcond<Precision>(a, n, erc1, ercinf);
                        
                        //
                        // normal
                        //
                        v = srcond::smatrixrcond<Precision>(a2, n, upperin);
                        if( v<=erc1/threshold50 )
                        {
                            q50(0) = q50(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                        }
                        if( v<=erc1/threshold90 )
                        {
                            q90(0) = q90(0)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                        }
                        errless = errless || v<erc1/amp::ampf<Precision>("1.001");
                        
                        //
                        // LDLT
                        //
                        ldlt::smatrixldlt<Precision>(a3, n, upperin, p);
                        v = srcond::smatrixldltrcond<Precision>(a3, p, n, upperin);
                        if( v<=erc1/threshold50 )
                        {
                            q50(1) = q50(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                        }
                        if( v<=erc1/threshold90 )
                        {
                            q90(1) = q90(1)+amp::ampf<Precision>(1)/amp::ampf<Precision>(passcount);
                        }
                        errless = errless || v<erc1/amp::ampf<Precision>("1.001");
                    }
                    for(i=0; i<=1; i++)
                    {
                        err50 = err50 || q50(i)<amp::ampf<Precision>("0.50");
                        err90 = err90 || q90(i)<amp::ampf<Precision>("0.90");
                    }
                }
            }
        }
        
        //
        // report
        //
        waserrors = err50 || err90 || errless;
        if( !silent )
        {
            printf("TESTING RCOND (LDLT)\n");
            printf("50%% within [0.5*cond,cond]:              ");
            if( err50 || errless )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("90%% within [0.1*cond,cond]               ");
            if( err90 || errless )
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
    reference RCond
    *************************************************************************/
    template<unsigned int Precision>
    void refrcond(const ap::template_2d_array< amp::ampf<Precision> >& a,
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
        makeacopy<Precision>(a, n, n, inva);
        if( !invmat<Precision>(inva, n) )
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
        rc1 = amp::ampf<Precision>("1.0")/(nrm1inva*nrm1a);
        rcinf = amp::ampf<Precision>("1.0")/(nrminfinva*nrminfa);
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testrcondldltunit_test_silent()
    {
        bool result;


        result = testrcondldlt<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testrcondldltunit_test()
    {
        bool result;


        result = testrcondldlt<Precision>(false);
        return result;
    }
} // namespace

#endif
