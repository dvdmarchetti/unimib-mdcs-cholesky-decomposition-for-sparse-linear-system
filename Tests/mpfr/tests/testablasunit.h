
#ifndef _testablasunit_h
#define _testablasunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "ablasf.h"
#include "ablas.h"
namespace testablasunit
{
    template<unsigned int Precision>
    bool testablas(bool silent);
    template<unsigned int Precision>
    void refcmatrixrighttrsm(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::campf<Precision> >& x,
        int i2,
        int j2);
    template<unsigned int Precision>
    void refcmatrixlefttrsm(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::campf<Precision> >& x,
        int i2,
        int j2);
    template<unsigned int Precision>
    void refrmatrixrighttrsm(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int i2,
        int j2);
    template<unsigned int Precision>
    void refrmatrixlefttrsm(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int i2,
        int j2);
    template<unsigned int Precision>
    bool internalcmatrixtrinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunittriangular);
    template<unsigned int Precision>
    bool internalrmatrixtrinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunittriangular);
    template<unsigned int Precision>
    void refcmatrixsyrk(int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::campf<Precision> >& c,
        int ic,
        int jc,
        bool isupper);
    template<unsigned int Precision>
    void refrmatrixsyrk(int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ic,
        int jc,
        bool isupper);
    template<unsigned int Precision>
    void refcmatrixgemm(int m,
        int n,
        int k,
        amp::campf<Precision> alpha,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int ib,
        int jb,
        int optypeb,
        amp::campf<Precision> beta,
        ap::template_2d_array< amp::campf<Precision> >& c,
        int ic,
        int jc);
    template<unsigned int Precision>
    void refrmatrixgemm(int m,
        int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int ib,
        int jb,
        int optypeb,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ic,
        int jc);
    template<unsigned int Precision>
    void naivematrixmatrixmultiply(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ai1,
        int ai2,
        int aj1,
        int aj2,
        bool transa,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int bi1,
        int bi2,
        int bj1,
        int bj2,
        bool transb,
        amp::ampf<Precision> alpha,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ci1,
        int ci2,
        int cj1,
        int cj2,
        amp::ampf<Precision> beta);
    template<unsigned int Precision>
    bool testtrsm(int minn,
        int maxn);
    template<unsigned int Precision>
    bool testsyrk(int minn,
        int maxn);
    template<unsigned int Precision>
    bool testgemm(int minn,
        int maxn);
    template<unsigned int Precision>
    bool testtrans(int minn,
        int maxn);
    template<unsigned int Precision>
    bool testrank1(int minn,
        int maxn);
    template<unsigned int Precision>
    bool testmv(int minn,
        int maxn);
    template<unsigned int Precision>
    bool testcopy(int minn,
        int maxn);
    template<unsigned int Precision>
    bool testablasunit_test_silent();
    template<unsigned int Precision>
    bool testablasunit_test();


    template<unsigned int Precision>
    bool testablas(bool silent)
    {
        bool result;
        amp::ampf<Precision> threshold;
        bool trsmerrors;
        bool syrkerrors;
        bool gemmerrors;
        bool transerrors;
        bool rank1errors;
        bool mverrors;
        bool copyerrors;
        bool waserrors;
        ap::template_2d_array< amp::ampf<Precision> > ra;


        trsmerrors = false;
        syrkerrors = false;
        gemmerrors = false;
        transerrors = false;
        rank1errors = false;
        mverrors = false;
        copyerrors = false;
        waserrors = false;
        threshold = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        trsmerrors = trsmerrors || testtrsm<Precision>(1, 3*ablas::ablasblocksize<Precision>(ra)+1);
        syrkerrors = syrkerrors || testsyrk<Precision>(1, 3*ablas::ablasblocksize<Precision>(ra)+1);
        gemmerrors = gemmerrors || testgemm<Precision>(1, 3*ablas::ablasblocksize<Precision>(ra)+1);
        transerrors = transerrors || testtrans<Precision>(1, 3*ablas::ablasblocksize<Precision>(ra)+1);
        rank1errors = rank1errors || testrank1<Precision>(1, 3*ablas::ablasblocksize<Precision>(ra)+1);
        mverrors = mverrors || testmv<Precision>(1, 3*ablas::ablasblocksize<Precision>(ra)+1);
        copyerrors = copyerrors || testcopy<Precision>(1, 3*ablas::ablasblocksize<Precision>(ra)+1);
        
        //
        // report
        //
        waserrors = trsmerrors || syrkerrors || gemmerrors || transerrors || rank1errors || mverrors || copyerrors;
        if( !silent )
        {
            printf("TESTING ABLAS\n");
            printf("* TRSM:                                  ");
            if( trsmerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* SYRK:                                  ");
            if( syrkerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* GEMM:                                  ");
            if( gemmerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* TRANS:                                 ");
            if( transerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* RANK1:                                 ");
            if( rank1errors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* MV:                                    ");
            if( mverrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* COPY:                                  ");
            if( copyerrors )
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
    Reference implementation

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void refcmatrixrighttrsm(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::campf<Precision> >& x,
        int i2,
        int j2)
    {
        ap::template_2d_array< amp::campf<Precision> > a1;
        ap::template_2d_array< amp::campf<Precision> > a2;
        ap::template_1d_array< amp::campf<Precision> > tx;
        int i;
        int j;
        amp::campf<Precision> vc;
        bool rupper;
        int i_;
        int i1_;


        if( n*m==0 )
        {
            return;
        }
        a1.setlength(n, n);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a1(i,j) = 0;
            }
        }
        if( isupper )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=i; j<=n-1; j++)
                {
                    a1(i,j) = a(i1+i,j1+j);
                }
            }
        }
        else
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=i; j++)
                {
                    a1(i,j) = a(i1+i,j1+j);
                }
            }
        }
        rupper = isupper;
        if( isunit )
        {
            for(i=0; i<=n-1; i++)
            {
                a1(i,i) = 1;
            }
        }
        a2.setlength(n, n);
        if( optype==0 )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a2(i,j) = a1(i,j);
                }
            }
        }
        if( optype==1 )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a2(i,j) = a1(j,i);
                }
            }
            rupper = !rupper;
        }
        if( optype==2 )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a2(i,j) = amp::conj<Precision>(a1(j,i));
                }
            }
            rupper = !rupper;
        }
        internalcmatrixtrinverse<Precision>(a2, n, rupper, false);
        tx.setlength(n);
        for(i=0; i<=m-1; i++)
        {
            i1_ = (j2) - (0);
            for(i_=0; i_<=n-1;i_++)
            {
                tx(i_) = x(i2+i,i_+i1_);
            }
            for(j=0; j<=n-1; j++)
            {
                vc = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    vc += tx(i_)*a2(i_,j);
                }
                x(i2+i,j2+j) = vc;
            }
        }
    }


    /*************************************************************************
    Reference implementation

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void refcmatrixlefttrsm(int m,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::campf<Precision> >& x,
        int i2,
        int j2)
    {
        ap::template_2d_array< amp::campf<Precision> > a1;
        ap::template_2d_array< amp::campf<Precision> > a2;
        ap::template_1d_array< amp::campf<Precision> > tx;
        int i;
        int j;
        amp::campf<Precision> vc;
        bool rupper;
        int i_;
        int i1_;


        if( n*m==0 )
        {
            return;
        }
        a1.setlength(m, m);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                a1(i,j) = 0;
            }
        }
        if( isupper )
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=i; j<=m-1; j++)
                {
                    a1(i,j) = a(i1+i,j1+j);
                }
            }
        }
        else
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=i; j++)
                {
                    a1(i,j) = a(i1+i,j1+j);
                }
            }
        }
        rupper = isupper;
        if( isunit )
        {
            for(i=0; i<=m-1; i++)
            {
                a1(i,i) = 1;
            }
        }
        a2.setlength(m, m);
        if( optype==0 )
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    a2(i,j) = a1(i,j);
                }
            }
        }
        if( optype==1 )
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    a2(i,j) = a1(j,i);
                }
            }
            rupper = !rupper;
        }
        if( optype==2 )
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    a2(i,j) = amp::conj<Precision>(a1(j,i));
                }
            }
            rupper = !rupper;
        }
        internalcmatrixtrinverse<Precision>(a2, m, rupper, false);
        tx.setlength(m);
        for(j=0; j<=n-1; j++)
        {
            i1_ = (i2) - (0);
            for(i_=0; i_<=m-1;i_++)
            {
                tx(i_) = x(i_+i1_,j2+j);
            }
            for(i=0; i<=m-1; i++)
            {
                vc = 0.0;
                for(i_=0; i_<=m-1;i_++)
                {
                    vc += a2(i,i_)*tx(i_);
                }
                x(i2+i,j2+j) = vc;
            }
        }
    }


    /*************************************************************************
    Reference implementation

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void refrmatrixrighttrsm(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int i2,
        int j2)
    {
        ap::template_2d_array< amp::ampf<Precision> > a1;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        ap::template_1d_array< amp::ampf<Precision> > tx;
        int i;
        int j;
        amp::ampf<Precision> vr;
        bool rupper;


        if( n*m==0 )
        {
            return;
        }
        a1.setlength(n, n);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a1(i,j) = 0;
            }
        }
        if( isupper )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=i; j<=n-1; j++)
                {
                    a1(i,j) = a(i1+i,j1+j);
                }
            }
        }
        else
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=i; j++)
                {
                    a1(i,j) = a(i1+i,j1+j);
                }
            }
        }
        rupper = isupper;
        if( isunit )
        {
            for(i=0; i<=n-1; i++)
            {
                a1(i,i) = 1;
            }
        }
        a2.setlength(n, n);
        if( optype==0 )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a2(i,j) = a1(i,j);
                }
            }
        }
        if( optype==1 )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a2(i,j) = a1(j,i);
                }
            }
            rupper = !rupper;
        }
        internalrmatrixtrinverse<Precision>(a2, n, rupper, false);
        tx.setlength(n);
        for(i=0; i<=m-1; i++)
        {
            amp::vmove(tx.getvector(0, n-1), x.getrow(i2+i, j2, j2+n-1));
            for(j=0; j<=n-1; j++)
            {
                vr = amp::vdotproduct(tx.getvector(0, n-1), a2.getcolumn(j, 0, n-1));
                x(i2+i,j2+j) = vr;
            }
        }
    }


    /*************************************************************************
    Reference implementation

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void refrmatrixlefttrsm(int m,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int i1,
        int j1,
        bool isupper,
        bool isunit,
        int optype,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int i2,
        int j2)
    {
        ap::template_2d_array< amp::ampf<Precision> > a1;
        ap::template_2d_array< amp::ampf<Precision> > a2;
        ap::template_1d_array< amp::ampf<Precision> > tx;
        int i;
        int j;
        amp::ampf<Precision> vr;
        bool rupper;


        if( n*m==0 )
        {
            return;
        }
        a1.setlength(m, m);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                a1(i,j) = 0;
            }
        }
        if( isupper )
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=i; j<=m-1; j++)
                {
                    a1(i,j) = a(i1+i,j1+j);
                }
            }
        }
        else
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=i; j++)
                {
                    a1(i,j) = a(i1+i,j1+j);
                }
            }
        }
        rupper = isupper;
        if( isunit )
        {
            for(i=0; i<=m-1; i++)
            {
                a1(i,i) = 1;
            }
        }
        a2.setlength(m, m);
        if( optype==0 )
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    a2(i,j) = a1(i,j);
                }
            }
        }
        if( optype==1 )
        {
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    a2(i,j) = a1(j,i);
                }
            }
            rupper = !rupper;
        }
        internalrmatrixtrinverse<Precision>(a2, m, rupper, false);
        tx.setlength(m);
        for(j=0; j<=n-1; j++)
        {
            amp::vmove(tx.getvector(0, m-1), x.getcolumn(j2+j, i2, i2+m-1));
            for(i=0; i<=m-1; i++)
            {
                vr = amp::vdotproduct(a2.getrow(i, 0, m-1), tx.getvector(0, m-1));
                x(i2+i,j2+j) = vr;
            }
        }
    }


    /*************************************************************************
    Internal subroutine.
    Triangular matrix inversion

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    bool internalcmatrixtrinverse(ap::template_2d_array< amp::campf<Precision> >& a,
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
                        if( i+1<j )
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
                if( j+1<n )
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
    Internal subroutine.
    Triangular matrix inversion

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    bool internalrmatrixtrinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
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
    Reference SYRK subroutine.

      -- ALGLIB routine --
         16.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void refcmatrixsyrk(int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::campf<Precision> >& c,
        int ic,
        int jc,
        bool isupper)
    {
        ap::template_2d_array< amp::campf<Precision> > ae;
        int i;
        int j;
        amp::campf<Precision> vc;
        int i_;


        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( isupper && j>=i || !isupper && j<=i )
                {
                    if( beta==0 )
                    {
                        c(i+ic,j+jc) = 0;
                    }
                    else
                    {
                        c(i+ic,j+jc) = c(i+ic,j+jc)*beta;
                    }
                }
            }
        }
        if( alpha==0 )
        {
            return;
        }
        if( n*k>0 )
        {
            ae.setlength(n, k);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                if( optypea==0 )
                {
                    ae(i,j) = a(ia+i,ja+j);
                }
                if( optypea==2 )
                {
                    ae(i,j) = amp::conj<Precision>(a(ia+j,ja+i));
                }
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                vc = 0;
                if( k>0 )
                {
                    vc = 0.0;
                    for(i_=0; i_<=k-1;i_++)
                    {
                        vc += ae(i,i_)*amp::conj(ae(j,i_));
                    }
                }
                vc = alpha*vc;
                if( isupper && j>=i )
                {
                    c(ic+i,jc+j) = vc+c(ic+i,jc+j);
                }
                if( !isupper && j<=i )
                {
                    c(ic+i,jc+j) = vc+c(ic+i,jc+j);
                }
            }
        }
    }


    /*************************************************************************
    Reference SYRK subroutine.

      -- ALGLIB routine --
         16.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void refrmatrixsyrk(int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ic,
        int jc,
        bool isupper)
    {
        ap::template_2d_array< amp::ampf<Precision> > ae;
        int i;
        int j;
        amp::ampf<Precision> vr;


        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( isupper && j>=i || !isupper && j<=i )
                {
                    if( beta==0 )
                    {
                        c(i+ic,j+jc) = 0;
                    }
                    else
                    {
                        c(i+ic,j+jc) = c(i+ic,j+jc)*beta;
                    }
                }
            }
        }
        if( alpha==0 )
        {
            return;
        }
        if( n*k>0 )
        {
            ae.setlength(n, k);
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                if( optypea==0 )
                {
                    ae(i,j) = a(ia+i,ja+j);
                }
                if( optypea==1 )
                {
                    ae(i,j) = a(ia+j,ja+i);
                }
            }
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                vr = 0;
                if( k>0 )
                {
                    vr = amp::vdotproduct(ae.getrow(i, 0, k-1), ae.getrow(j, 0, k-1));
                }
                vr = alpha*vr;
                if( isupper && j>=i )
                {
                    c(ic+i,jc+j) = vr+c(ic+i,jc+j);
                }
                if( !isupper && j<=i )
                {
                    c(ic+i,jc+j) = vr+c(ic+i,jc+j);
                }
            }
        }
    }


    /*************************************************************************
    Reference GEMM,
    ALGLIB subroutine
    *************************************************************************/
    template<unsigned int Precision>
    void refcmatrixgemm(int m,
        int n,
        int k,
        amp::campf<Precision> alpha,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int ib,
        int jb,
        int optypeb,
        amp::campf<Precision> beta,
        ap::template_2d_array< amp::campf<Precision> >& c,
        int ic,
        int jc)
    {
        ap::template_2d_array< amp::campf<Precision> > ae;
        ap::template_2d_array< amp::campf<Precision> > be;
        int i;
        int j;
        amp::campf<Precision> vc;
        int i_;


        ae.setlength(m, k);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                if( optypea==0 )
                {
                    ae(i,j) = a(ia+i,ja+j);
                }
                if( optypea==1 )
                {
                    ae(i,j) = a(ia+j,ja+i);
                }
                if( optypea==2 )
                {
                    ae(i,j) = amp::conj<Precision>(a(ia+j,ja+i));
                }
            }
        }
        be.setlength(k, n);
        for(i=0; i<=k-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( optypeb==0 )
                {
                    be(i,j) = b(ib+i,jb+j);
                }
                if( optypeb==1 )
                {
                    be(i,j) = b(ib+j,jb+i);
                }
                if( optypeb==2 )
                {
                    be(i,j) = amp::conj<Precision>(b(ib+j,jb+i));
                }
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                vc = 0.0;
                for(i_=0; i_<=k-1;i_++)
                {
                    vc += ae(i,i_)*be(i_,j);
                }
                vc = alpha*vc;
                if( beta!=0 )
                {
                    vc = vc+beta*c(ic+i,jc+j);
                }
                c(ic+i,jc+j) = vc;
            }
        }
    }


    /*************************************************************************
    Reference GEMM,
    ALGLIB subroutine
    *************************************************************************/
    template<unsigned int Precision>
    void refrmatrixgemm(int m,
        int n,
        int k,
        amp::ampf<Precision> alpha,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ia,
        int ja,
        int optypea,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int ib,
        int jb,
        int optypeb,
        amp::ampf<Precision> beta,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ic,
        int jc)
    {
        ap::template_2d_array< amp::ampf<Precision> > ae;
        ap::template_2d_array< amp::ampf<Precision> > be;
        int i;
        int j;
        amp::ampf<Precision> vc;


        ae.setlength(m, k);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                if( optypea==0 )
                {
                    ae(i,j) = a(ia+i,ja+j);
                }
                if( optypea==1 )
                {
                    ae(i,j) = a(ia+j,ja+i);
                }
            }
        }
        be.setlength(k, n);
        for(i=0; i<=k-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( optypeb==0 )
                {
                    be(i,j) = b(ib+i,jb+j);
                }
                if( optypeb==1 )
                {
                    be(i,j) = b(ib+j,jb+i);
                }
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                vc = amp::vdotproduct(ae.getrow(i, 0, k-1), be.getcolumn(j, 0, k-1));
                vc = alpha*vc;
                if( beta!=0 )
                {
                    vc = vc+beta*c(ic+i,jc+j);
                }
                c(ic+i,jc+j) = vc;
            }
        }
    }


    template<unsigned int Precision>
    void naivematrixmatrixmultiply(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int ai1,
        int ai2,
        int aj1,
        int aj2,
        bool transa,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int bi1,
        int bi2,
        int bj1,
        int bj2,
        bool transb,
        amp::ampf<Precision> alpha,
        ap::template_2d_array< amp::ampf<Precision> >& c,
        int ci1,
        int ci2,
        int cj1,
        int cj2,
        amp::ampf<Precision> beta)
    {
        int arows;
        int acols;
        int brows;
        int bcols;
        int i;
        int j;
        int k;
        int l;
        int r;
        amp::ampf<Precision> v;
        ap::template_1d_array< amp::ampf<Precision> > x1;
        ap::template_1d_array< amp::ampf<Precision> > x2;


        
        //
        // Setup
        //
        if( !transa )
        {
            arows = ai2-ai1+1;
            acols = aj2-aj1+1;
        }
        else
        {
            arows = aj2-aj1+1;
            acols = ai2-ai1+1;
        }
        if( !transb )
        {
            brows = bi2-bi1+1;
            bcols = bj2-bj1+1;
        }
        else
        {
            brows = bj2-bj1+1;
            bcols = bi2-bi1+1;
        }
        ap::ap_error::make_assertion(acols==brows);
        if( arows<=0 || acols<=0 || brows<=0 || bcols<=0 )
        {
            return;
        }
        l = arows;
        r = bcols;
        k = acols;
        x1.setbounds(1, k);
        x2.setbounds(1, k);
        for(i=1; i<=l; i++)
        {
            for(j=1; j<=r; j++)
            {
                if( !transa )
                {
                    if( !transb )
                    {
                        v = amp::vdotproduct(b.getcolumn(bj1+j-1, bi1, bi2), a.getrow(ai1+i-1, aj1, aj2));
                    }
                    else
                    {
                        v = amp::vdotproduct(b.getrow(bi1+j-1, bj1, bj2), a.getrow(ai1+i-1, aj1, aj2));
                    }
                }
                else
                {
                    if( !transb )
                    {
                        v = amp::vdotproduct(b.getcolumn(bj1+j-1, bi1, bi2), a.getcolumn(aj1+i-1, ai1, ai2));
                    }
                    else
                    {
                        v = amp::vdotproduct(b.getrow(bi1+j-1, bj1, bj2), a.getcolumn(aj1+i-1, ai1, ai2));
                    }
                }
                if( beta==0 )
                {
                    c(ci1+i-1,cj1+j-1) = alpha*v;
                }
                else
                {
                    c(ci1+i-1,cj1+j-1) = beta*c(ci1+i-1,cj1+j-1)+alpha*v;
                }
            }
        }
    }


    /*************************************************************************
    ?Matrix????TRSM tests

    Returns False for passed test, True - for failed
    *************************************************************************/
    template<unsigned int Precision>
    bool testtrsm(int minn,
        int maxn)
    {
        bool result;
        int n;
        int m;
        int mx;
        int i;
        int j;
        int optype;
        int uppertype;
        int unittype;
        int xoffsi;
        int xoffsj;
        int aoffsitype;
        int aoffsjtype;
        int aoffsi;
        int aoffsj;
        ap::template_2d_array< amp::ampf<Precision> > refra;
        ap::template_2d_array< amp::ampf<Precision> > refrxl;
        ap::template_2d_array< amp::ampf<Precision> > refrxr;
        ap::template_2d_array< amp::campf<Precision> > refca;
        ap::template_2d_array< amp::campf<Precision> > refcxl;
        ap::template_2d_array< amp::campf<Precision> > refcxr;
        ap::template_2d_array< amp::ampf<Precision> > ra;
        ap::template_2d_array< amp::campf<Precision> > ca;
        ap::template_2d_array< amp::ampf<Precision> > rxr1;
        ap::template_2d_array< amp::ampf<Precision> > rxl1;
        ap::template_2d_array< amp::campf<Precision> > cxr1;
        ap::template_2d_array< amp::campf<Precision> > cxl1;
        ap::template_2d_array< amp::ampf<Precision> > rxr2;
        ap::template_2d_array< amp::ampf<Precision> > rxl2;
        ap::template_2d_array< amp::campf<Precision> > cxr2;
        ap::template_2d_array< amp::campf<Precision> > cxl2;
        amp::ampf<Precision> threshold;


        threshold = amp::sqr<Precision>(amp::ampf<Precision>(maxn))*100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        result = false;
        for(mx=minn; mx<=maxn; mx++)
        {
            
            //
            // Select random M/N in [1,MX] such that max(M,N)=MX
            //
            m = 1+ap::randominteger(mx);
            n = 1+ap::randominteger(mx);
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                m = mx;
            }
            else
            {
                n = mx;
            }
            
            //
            // Initialize RefRA/RefCA by random matrices whose upper
            // and lower triangle submatrices are non-degenerate
            // well-conditioned matrices.
            //
            // Matrix size is 2Mx2M (four copies of same MxM matrix
            // to test different offsets)
            //
            refra.setlength(2*m, 2*m);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    refra(i,j) = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                }
            }
            for(i=0; i<=m-1; i++)
            {
                refra(i,i) = (2*ap::randominteger(1)-1)*(2*m+amp::ampf<Precision>::getRandom());
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    refra(i+m,j) = refra(i,j);
                    refra(i,j+m) = refra(i,j);
                    refra(i+m,j+m) = refra(i,j);
                }
            }
            refca.setlength(2*m, 2*m);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    refca(i,j).x = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                    refca(i,j).y = amp::ampf<Precision>("0.2")*amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.1");
                }
            }
            for(i=0; i<=m-1; i++)
            {
                refca(i,i).x = (2*ap::randominteger(2)-1)*(2*m+amp::ampf<Precision>::getRandom());
                refca(i,i).y = (2*ap::randominteger(2)-1)*(2*m+amp::ampf<Precision>::getRandom());
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    refca(i+m,j) = refca(i,j);
                    refca(i,j+m) = refca(i,j);
                    refca(i+m,j+m) = refca(i,j);
                }
            }
            
            //
            // Generate random XL/XR.
            //
            // XR is NxM matrix (matrix for 'Right' subroutines)
            // XL is MxN matrix (matrix for 'Left' subroutines)
            //
            refrxr.setlength(n, m);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    refrxr(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            refrxl.setlength(m, n);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    refrxl(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            refcxr.setlength(n, m);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    refcxr(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refcxr(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            refcxl.setlength(m, n);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    refcxl(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refcxl(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            
            //
            // test different types of operations, offsets, and so on...
            //
            // to avoid unnecessary slowdown we don't test ALL possible
            // combinations of operation types. We just generate one random
            // set of parameters and test it.
            //
            ra.setlength(2*m, 2*m);
            rxr1.setlength(n, m);
            rxr2.setlength(n, m);
            rxl1.setlength(m, n);
            rxl2.setlength(m, n);
            ca.setlength(2*m, 2*m);
            cxr1.setlength(n, m);
            cxr2.setlength(n, m);
            cxl1.setlength(m, n);
            cxl2.setlength(m, n);
            optype = ap::randominteger(3);
            uppertype = ap::randominteger(2);
            unittype = ap::randominteger(2);
            xoffsi = ap::randominteger(2);
            xoffsj = ap::randominteger(2);
            aoffsitype = ap::randominteger(2);
            aoffsjtype = ap::randominteger(2);
            aoffsi = m*aoffsitype;
            aoffsj = m*aoffsjtype;
            
            //
            // copy A, XR, XL (fill unused parts with random garbage)
            //
            for(i=0; i<=2*m-1; i++)
            {
                for(j=0; j<=2*m-1; j++)
                {
                    if( i>=aoffsi && i<aoffsi+m && j>=aoffsj && j<aoffsj+m )
                    {
                        ca(i,j) = refca(i,j);
                        ra(i,j) = refra(i,j);
                    }
                    else
                    {
                        ca(i,j) = amp::ampf<Precision>::getRandom();
                        ra(i,j) = amp::ampf<Precision>::getRandom();
                    }
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    if( i>=xoffsi && j>=xoffsj )
                    {
                        cxr1(i,j) = refcxr(i,j);
                        cxr2(i,j) = refcxr(i,j);
                        rxr1(i,j) = refrxr(i,j);
                        rxr2(i,j) = refrxr(i,j);
                    }
                    else
                    {
                        cxr1(i,j) = amp::ampf<Precision>::getRandom();
                        cxr2(i,j) = cxr1(i,j);
                        rxr1(i,j) = amp::ampf<Precision>::getRandom();
                        rxr2(i,j) = rxr1(i,j);
                    }
                }
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( i>=xoffsi && j>=xoffsj )
                    {
                        cxl1(i,j) = refcxl(i,j);
                        cxl2(i,j) = refcxl(i,j);
                        rxl1(i,j) = refrxl(i,j);
                        rxl2(i,j) = refrxl(i,j);
                    }
                    else
                    {
                        cxl1(i,j) = amp::ampf<Precision>::getRandom();
                        cxl2(i,j) = cxl1(i,j);
                        rxl1(i,j) = amp::ampf<Precision>::getRandom();
                        rxl2(i,j) = rxl1(i,j);
                    }
                }
            }
            
            //
            // Test CXR
            //
            ablas::cmatrixrighttrsm<Precision>(n-xoffsi, m-xoffsj, ca, aoffsi, aoffsj, uppertype==0, unittype==0, optype, cxr1, xoffsi, xoffsj);
            refcmatrixrighttrsm<Precision>(n-xoffsi, m-xoffsj, ca, aoffsi, aoffsj, uppertype==0, unittype==0, optype, cxr2, xoffsi, xoffsj);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    result = result || amp::abscomplex<Precision>(cxr1(i,j)-cxr2(i,j))>threshold;
                }
            }
            
            //
            // Test CXL
            //
            ablas::cmatrixlefttrsm<Precision>(m-xoffsi, n-xoffsj, ca, aoffsi, aoffsj, uppertype==0, unittype==0, optype, cxl1, xoffsi, xoffsj);
            refcmatrixlefttrsm<Precision>(m-xoffsi, n-xoffsj, ca, aoffsi, aoffsj, uppertype==0, unittype==0, optype, cxl2, xoffsi, xoffsj);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    result = result || amp::abscomplex<Precision>(cxl1(i,j)-cxl2(i,j))>threshold;
                }
            }
            if( optype<2 )
            {
                
                //
                // Test RXR
                //
                ablas::rmatrixrighttrsm<Precision>(n-xoffsi, m-xoffsj, ra, aoffsi, aoffsj, uppertype==0, unittype==0, optype, rxr1, xoffsi, xoffsj);
                refrmatrixrighttrsm<Precision>(n-xoffsi, m-xoffsj, ra, aoffsi, aoffsj, uppertype==0, unittype==0, optype, rxr2, xoffsi, xoffsj);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        result = result || amp::abs<Precision>(rxr1(i,j)-rxr2(i,j))>threshold;
                    }
                }
                
                //
                // Test RXL
                //
                ablas::rmatrixlefttrsm<Precision>(m-xoffsi, n-xoffsj, ra, aoffsi, aoffsj, uppertype==0, unittype==0, optype, rxl1, xoffsi, xoffsj);
                refrmatrixlefttrsm<Precision>(m-xoffsi, n-xoffsj, ra, aoffsi, aoffsj, uppertype==0, unittype==0, optype, rxl2, xoffsi, xoffsj);
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        result = result || amp::abs<Precision>(rxl1(i,j)-rxl2(i,j))>threshold;
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    SYRK tests

    Returns False for passed test, True - for failed
    *************************************************************************/
    template<unsigned int Precision>
    bool testsyrk(int minn,
        int maxn)
    {
        bool result;
        int n;
        int k;
        int mx;
        int i;
        int j;
        int uppertype;
        int xoffsi;
        int xoffsj;
        int aoffsitype;
        int aoffsjtype;
        int aoffsi;
        int aoffsj;
        int alphatype;
        int betatype;
        ap::template_2d_array< amp::ampf<Precision> > refra;
        ap::template_2d_array< amp::ampf<Precision> > refrc;
        ap::template_2d_array< amp::campf<Precision> > refca;
        ap::template_2d_array< amp::campf<Precision> > refcc;
        amp::ampf<Precision> alpha;
        amp::ampf<Precision> beta;
        ap::template_2d_array< amp::ampf<Precision> > ra1;
        ap::template_2d_array< amp::ampf<Precision> > ra2;
        ap::template_2d_array< amp::campf<Precision> > ca1;
        ap::template_2d_array< amp::campf<Precision> > ca2;
        ap::template_2d_array< amp::ampf<Precision> > rc;
        ap::template_2d_array< amp::ampf<Precision> > rct;
        ap::template_2d_array< amp::campf<Precision> > cc;
        ap::template_2d_array< amp::campf<Precision> > cct;
        amp::ampf<Precision> threshold;


        threshold = maxn*100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        result = false;
        for(mx=minn; mx<=maxn; mx++)
        {
            
            //
            // Select random M/N in [1,MX] such that max(M,N)=MX
            //
            k = 1+ap::randominteger(mx);
            n = 1+ap::randominteger(mx);
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                k = mx;
            }
            else
            {
                n = mx;
            }
            
            //
            // Initialize RefRA/RefCA by random Hermitian matrices,
            // RefRC/RefCC by random matrices
            //
            // RA/CA size is 2Nx2N (four copies of same NxN matrix
            // to test different offsets)
            //
            refra.setlength(2*n, 2*n);
            refca.setlength(2*n, 2*n);
            for(i=0; i<=n-1; i++)
            {
                refra(i,i) = 2*amp::ampf<Precision>::getRandom()-1;
                refca(i,i) = 2*amp::ampf<Precision>::getRandom()-1;
                for(j=i+1; j<=n-1; j++)
                {
                    refra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    refra(j,i) = refra(i,j);
                    refca(j,i) = amp::conj<Precision>(refca(i,j));
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    refra(i+n,j) = refra(i,j);
                    refra(i,j+n) = refra(i,j);
                    refra(i+n,j+n) = refra(i,j);
                    refca(i+n,j) = refca(i,j);
                    refca(i,j+n) = refca(i,j);
                    refca(i+n,j+n) = refca(i,j);
                }
            }
            refrc.setlength(n, k);
            refcc.setlength(n, k);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=k-1; j++)
                {
                    refrc(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    refcc(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refcc(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            
            //
            // test different types of operations, offsets, and so on...
            //
            // to avoid unnecessary slowdown we don't test ALL possible
            // combinations of operation types. We just generate one random
            // set of parameters and test it.
            //
            ra1.setlength(2*n, 2*n);
            ra2.setlength(2*n, 2*n);
            ca1.setlength(2*n, 2*n);
            ca2.setlength(2*n, 2*n);
            rc.setlength(n, k);
            rct.setlength(k, n);
            cc.setlength(n, k);
            cct.setlength(k, n);
            uppertype = ap::randominteger(2);
            xoffsi = ap::randominteger(2);
            xoffsj = ap::randominteger(2);
            aoffsitype = ap::randominteger(2);
            aoffsjtype = ap::randominteger(2);
            alphatype = ap::randominteger(2);
            betatype = ap::randominteger(2);
            aoffsi = n*aoffsitype;
            aoffsj = n*aoffsjtype;
            alpha = alphatype*(2*amp::ampf<Precision>::getRandom()-1);
            beta = betatype*(2*amp::ampf<Precision>::getRandom()-1);
            
            //
            // copy A, C (fill unused parts with random garbage)
            //
            for(i=0; i<=2*n-1; i++)
            {
                for(j=0; j<=2*n-1; j++)
                {
                    if( i>=aoffsi && i<aoffsi+n && j>=aoffsj && j<aoffsj+n )
                    {
                        ca1(i,j) = refca(i,j);
                        ca2(i,j) = refca(i,j);
                        ra1(i,j) = refra(i,j);
                        ra2(i,j) = refra(i,j);
                    }
                    else
                    {
                        ca1(i,j) = amp::ampf<Precision>::getRandom();
                        ca2(i,j) = ca1(i,j);
                        ra1(i,j) = amp::ampf<Precision>::getRandom();
                        ra2(i,j) = ra1(i,j);
                    }
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=k-1; j++)
                {
                    if( i>=xoffsi && j>=xoffsj )
                    {
                        rc(i,j) = refrc(i,j);
                        rct(j,i) = refrc(i,j);
                        cc(i,j) = refcc(i,j);
                        cct(j,i) = refcc(i,j);
                    }
                    else
                    {
                        rc(i,j) = amp::ampf<Precision>::getRandom();
                        rct(j,i) = rc(i,j);
                        cc(i,j) = amp::ampf<Precision>::getRandom();
                        cct(j,i) = cct(j,i);
                    }
                }
            }
            
            //
            // Test complex
            // Only one of transform types is selected and tested
            //
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                ablas::cmatrixsyrk<Precision>(n-xoffsi, k-xoffsj, alpha, cc, xoffsi, xoffsj, 0, beta, ca1, aoffsi, aoffsj, uppertype==0);
                refcmatrixsyrk<Precision>(n-xoffsi, k-xoffsj, alpha, cc, xoffsi, xoffsj, 0, beta, ca2, aoffsi, aoffsj, uppertype==0);
            }
            else
            {
                ablas::cmatrixsyrk<Precision>(n-xoffsi, k-xoffsj, alpha, cct, xoffsj, xoffsi, 2, beta, ca1, aoffsi, aoffsj, uppertype==0);
                refcmatrixsyrk<Precision>(n-xoffsi, k-xoffsj, alpha, cct, xoffsj, xoffsi, 2, beta, ca2, aoffsi, aoffsj, uppertype==0);
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    result = result || amp::abscomplex<Precision>(ca1(i,j)-ca2(i,j))>threshold;
                }
            }
            
            //
            // Test real
            // Only one of transform types is selected and tested
            //
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                ablas::rmatrixsyrk<Precision>(n-xoffsi, k-xoffsj, alpha, rc, xoffsi, xoffsj, 0, beta, ra1, aoffsi, aoffsj, uppertype==0);
                refrmatrixsyrk<Precision>(n-xoffsi, k-xoffsj, alpha, rc, xoffsi, xoffsj, 0, beta, ra2, aoffsi, aoffsj, uppertype==0);
            }
            else
            {
                ablas::rmatrixsyrk<Precision>(n-xoffsi, k-xoffsj, alpha, rct, xoffsj, xoffsi, 1, beta, ra1, aoffsi, aoffsj, uppertype==0);
                refrmatrixsyrk<Precision>(n-xoffsi, k-xoffsj, alpha, rct, xoffsj, xoffsi, 1, beta, ra2, aoffsi, aoffsj, uppertype==0);
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    result = result || amp::abs<Precision>(ra1(i,j)-ra2(i,j))>threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    GEMM tests

    Returns False for passed test, True - for failed
    *************************************************************************/
    template<unsigned int Precision>
    bool testgemm(int minn,
        int maxn)
    {
        bool result;
        int m;
        int n;
        int k;
        int mx;
        int i;
        int j;
        int aoffsi;
        int aoffsj;
        int aoptype;
        int aoptyper;
        int boffsi;
        int boffsj;
        int boptype;
        int boptyper;
        int coffsi;
        int coffsj;
        ap::template_2d_array< amp::ampf<Precision> > refra;
        ap::template_2d_array< amp::ampf<Precision> > refrb;
        ap::template_2d_array< amp::ampf<Precision> > refrc;
        ap::template_2d_array< amp::campf<Precision> > refca;
        ap::template_2d_array< amp::campf<Precision> > refcb;
        ap::template_2d_array< amp::campf<Precision> > refcc;
        amp::ampf<Precision> alphar;
        amp::ampf<Precision> betar;
        amp::campf<Precision> alphac;
        amp::campf<Precision> betac;
        ap::template_2d_array< amp::ampf<Precision> > rc1;
        ap::template_2d_array< amp::ampf<Precision> > rc2;
        ap::template_2d_array< amp::campf<Precision> > cc1;
        ap::template_2d_array< amp::campf<Precision> > cc2;
        amp::ampf<Precision> threshold;


        threshold = maxn*100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        result = false;
        for(mx=minn; mx<=maxn; mx++)
        {
            
            //
            // Select random M/N/K in [1,MX] such that max(M,N,K)=MX
            //
            m = 1+ap::randominteger(mx);
            n = 1+ap::randominteger(mx);
            k = 1+ap::randominteger(mx);
            i = ap::randominteger(3);
            if( i==0 )
            {
                m = mx;
            }
            if( i==1 )
            {
                n = mx;
            }
            if( i==2 )
            {
                k = mx;
            }
            
            //
            // Initialize A/B/C by random matrices with size (MaxN+1)*(MaxN+1)
            //
            refra.setlength(maxn+1, maxn+1);
            refrb.setlength(maxn+1, maxn+1);
            refrc.setlength(maxn+1, maxn+1);
            refca.setlength(maxn+1, maxn+1);
            refcb.setlength(maxn+1, maxn+1);
            refcc.setlength(maxn+1, maxn+1);
            for(i=0; i<=maxn; i++)
            {
                for(j=0; j<=maxn; j++)
                {
                    refra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    refrb(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    refrc(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    refcb(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refcb(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    refcc(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refcc(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            
            //
            // test different types of operations, offsets, and so on...
            //
            // to avoid unnecessary slowdown we don't test ALL possible
            // combinations of operation types. We just generate one random
            // set of parameters and test it.
            //
            rc1.setlength(maxn+1, maxn+1);
            rc2.setlength(maxn+1, maxn+1);
            cc1.setlength(maxn+1, maxn+1);
            cc2.setlength(maxn+1, maxn+1);
            aoffsi = ap::randominteger(2);
            aoffsj = ap::randominteger(2);
            aoptype = ap::randominteger(3);
            aoptyper = ap::randominteger(2);
            boffsi = ap::randominteger(2);
            boffsj = ap::randominteger(2);
            boptype = ap::randominteger(3);
            boptyper = ap::randominteger(2);
            coffsi = ap::randominteger(2);
            coffsj = ap::randominteger(2);
            alphar = ap::randominteger(2)*(2*amp::ampf<Precision>::getRandom()-1);
            betar = ap::randominteger(2)*(2*amp::ampf<Precision>::getRandom()-1);
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                alphac.x = 2*amp::ampf<Precision>::getRandom()-1;
                alphac.y = 2*amp::ampf<Precision>::getRandom()-1;
            }
            else
            {
                alphac = 0;
            }
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                betac.x = 2*amp::ampf<Precision>::getRandom()-1;
                betac.y = 2*amp::ampf<Precision>::getRandom()-1;
            }
            else
            {
                betac = 0;
            }
            
            //
            // copy C
            //
            for(i=0; i<=maxn; i++)
            {
                for(j=0; j<=maxn; j++)
                {
                    rc1(i,j) = refrc(i,j);
                    rc2(i,j) = refrc(i,j);
                    cc1(i,j) = refcc(i,j);
                    cc2(i,j) = refcc(i,j);
                }
            }
            
            //
            // Test complex
            //
            ablas::cmatrixgemm<Precision>(m, n, k, alphac, refca, aoffsi, aoffsj, aoptype, refcb, boffsi, boffsj, boptype, betac, cc1, coffsi, coffsj);
            refcmatrixgemm<Precision>(m, n, k, alphac, refca, aoffsi, aoffsj, aoptype, refcb, boffsi, boffsj, boptype, betac, cc2, coffsi, coffsj);
            for(i=0; i<=maxn; i++)
            {
                for(j=0; j<=maxn; j++)
                {
                    result = result || amp::abscomplex<Precision>(cc1(i,j)-cc2(i,j))>threshold;
                }
            }
            
            //
            // Test real
            //
            ablas::rmatrixgemm<Precision>(m, n, k, alphar, refra, aoffsi, aoffsj, aoptyper, refrb, boffsi, boffsj, boptyper, betar, rc1, coffsi, coffsj);
            refrmatrixgemm<Precision>(m, n, k, alphar, refra, aoffsi, aoffsj, aoptyper, refrb, boffsi, boffsj, boptyper, betar, rc2, coffsi, coffsj);
            for(i=0; i<=maxn; i++)
            {
                for(j=0; j<=maxn; j++)
                {
                    result = result || amp::abs<Precision>(rc1(i,j)-rc2(i,j))>threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    transpose tests

    Returns False for passed test, True - for failed
    *************************************************************************/
    template<unsigned int Precision>
    bool testtrans(int minn,
        int maxn)
    {
        bool result;
        int m;
        int n;
        int mx;
        int i;
        int j;
        int aoffsi;
        int aoffsj;
        int boffsi;
        int boffsj;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> threshold;
        ap::template_2d_array< amp::ampf<Precision> > refra;
        ap::template_2d_array< amp::ampf<Precision> > refrb;
        ap::template_2d_array< amp::campf<Precision> > refca;
        ap::template_2d_array< amp::campf<Precision> > refcb;


        result = false;
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        for(mx=minn; mx<=maxn; mx++)
        {
            
            //
            // Select random M/N in [1,MX] such that max(M,N)=MX
            // Generate random V1 and V2 which are used to fill
            // RefRB/RefCB with control values.
            //
            m = 1+ap::randominteger(mx);
            n = 1+ap::randominteger(mx);
            if( ap::randominteger(2)==0 )
            {
                m = mx;
            }
            else
            {
                n = mx;
            }
            v1 = amp::ampf<Precision>::getRandom();
            v2 = amp::ampf<Precision>::getRandom();
            
            //
            // Initialize A by random matrix with size (MaxN+1)*(MaxN+1)
            // Fill B with control values
            //
            refra.setlength(maxn+1, maxn+1);
            refrb.setlength(maxn+1, maxn+1);
            refca.setlength(maxn+1, maxn+1);
            refcb.setlength(maxn+1, maxn+1);
            for(i=0; i<=maxn; i++)
            {
                for(j=0; j<=maxn; j++)
                {
                    refra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    refrb(i,j) = i*v1+j*v2;
                    refcb(i,j) = i*v1+j*v2;
                }
            }
            
            //
            // test different offsets (zero or one)
            //
            // to avoid unnecessary slowdown we don't test ALL possible
            // combinations of operation types. We just generate one random
            // set of parameters and test it.
            //
            aoffsi = ap::randominteger(2);
            aoffsj = ap::randominteger(2);
            boffsi = ap::randominteger(2);
            boffsj = ap::randominteger(2);
            ablas::rmatrixtranspose<Precision>(m, n, refra, aoffsi, aoffsj, refrb, boffsi, boffsj);
            for(i=0; i<=maxn; i++)
            {
                for(j=0; j<=maxn; j++)
                {
                    if( i<boffsi || i>=boffsi+n || j<boffsj || j>=boffsj+m )
                    {
                        result = result || amp::abs<Precision>(refrb(i,j)-(v1*i+v2*j))>threshold;
                    }
                    else
                    {
                        result = result || amp::abs<Precision>(refrb(i,j)-refra(aoffsi+j-boffsj,aoffsj+i-boffsi))>threshold;
                    }
                }
            }
            ablas::cmatrixtranspose<Precision>(m, n, refca, aoffsi, aoffsj, refcb, boffsi, boffsj);
            for(i=0; i<=maxn; i++)
            {
                for(j=0; j<=maxn; j++)
                {
                    if( i<boffsi || i>=boffsi+n || j<boffsj || j>=boffsj+m )
                    {
                        result = result || amp::abscomplex<Precision>(refcb(i,j)-(v1*i+v2*j))>threshold;
                    }
                    else
                    {
                        result = result || amp::abscomplex<Precision>(refcb(i,j)-refca(aoffsi+j-boffsj,aoffsj+i-boffsi))>threshold;
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    rank-1tests

    Returns False for passed test, True - for failed
    *************************************************************************/
    template<unsigned int Precision>
    bool testrank1(int minn,
        int maxn)
    {
        bool result;
        int m;
        int n;
        int mx;
        int i;
        int j;
        int aoffsi;
        int aoffsj;
        int uoffs;
        int voffs;
        amp::ampf<Precision> threshold;
        ap::template_2d_array< amp::ampf<Precision> > refra;
        ap::template_2d_array< amp::ampf<Precision> > refrb;
        ap::template_2d_array< amp::campf<Precision> > refca;
        ap::template_2d_array< amp::campf<Precision> > refcb;
        ap::template_1d_array< amp::ampf<Precision> > ru;
        ap::template_1d_array< amp::ampf<Precision> > rv;
        ap::template_1d_array< amp::campf<Precision> > cu;
        ap::template_1d_array< amp::campf<Precision> > cv;


        result = false;
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        for(mx=minn; mx<=maxn; mx++)
        {
            
            //
            // Select random M/N in [1,MX] such that max(M,N)=MX
            //
            m = 1+ap::randominteger(mx);
            n = 1+ap::randominteger(mx);
            if( ap::randominteger(2)==0 )
            {
                m = mx;
            }
            else
            {
                n = mx;
            }
            
            //
            // Initialize A by random matrix with size (MaxN+1)*(MaxN+1)
            // Fill B with control values
            //
            refra.setlength(maxn+maxn, maxn+maxn);
            refrb.setlength(maxn+maxn, maxn+maxn);
            refca.setlength(maxn+maxn, maxn+maxn);
            refcb.setlength(maxn+maxn, maxn+maxn);
            for(i=0; i<=2*maxn-1; i++)
            {
                for(j=0; j<=2*maxn-1; j++)
                {
                    refra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    refrb(i,j) = refra(i,j);
                    refcb(i,j) = refca(i,j);
                }
            }
            ru.setlength(2*m);
            cu.setlength(2*m);
            for(i=0; i<=2*m-1; i++)
            {
                ru(i) = 2*amp::ampf<Precision>::getRandom()-1;
                cu(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                cu(i).y = 2*amp::ampf<Precision>::getRandom()-1;
            }
            rv.setlength(2*n);
            cv.setlength(2*n);
            for(i=0; i<=2*n-1; i++)
            {
                rv(i) = 2*amp::ampf<Precision>::getRandom()-1;
                cv(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                cv(i).y = 2*amp::ampf<Precision>::getRandom()-1;
            }
            
            //
            // test different offsets (zero or one)
            //
            // to avoid unnecessary slowdown we don't test ALL possible
            // combinations of operation types. We just generate one random
            // set of parameters and test it.
            //
            aoffsi = ap::randominteger(maxn);
            aoffsj = ap::randominteger(maxn);
            uoffs = ap::randominteger(m);
            voffs = ap::randominteger(n);
            ablas::cmatrixrank1<Precision>(m, n, refca, aoffsi, aoffsj, cu, uoffs, cv, voffs);
            for(i=0; i<=2*maxn-1; i++)
            {
                for(j=0; j<=2*maxn-1; j++)
                {
                    if( i<aoffsi || i>=aoffsi+m || j<aoffsj || j>=aoffsj+n )
                    {
                        result = result || amp::abscomplex<Precision>(refca(i,j)-refcb(i,j))>threshold;
                    }
                    else
                    {
                        result = result || amp::abscomplex<Precision>(refca(i,j)-(refcb(i,j)+cu(i-aoffsi+uoffs)*cv(j-aoffsj+voffs)))>threshold;
                    }
                }
            }
            ablas::rmatrixrank1<Precision>(m, n, refra, aoffsi, aoffsj, ru, uoffs, rv, voffs);
            for(i=0; i<=2*maxn-1; i++)
            {
                for(j=0; j<=2*maxn-1; j++)
                {
                    if( i<aoffsi || i>=aoffsi+m || j<aoffsj || j>=aoffsj+n )
                    {
                        result = result || amp::abs<Precision>(refra(i,j)-refrb(i,j))>threshold;
                    }
                    else
                    {
                        result = result || amp::abs<Precision>(refra(i,j)-(refrb(i,j)+ru(i-aoffsi+uoffs)*rv(j-aoffsj+voffs)))>threshold;
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    MV tests

    Returns False for passed test, True - for failed
    *************************************************************************/
    template<unsigned int Precision>
    bool testmv(int minn,
        int maxn)
    {
        bool result;
        int m;
        int n;
        int mx;
        int i;
        int j;
        int aoffsi;
        int aoffsj;
        int xoffs;
        int yoffs;
        int opca;
        int opra;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> rv1;
        amp::ampf<Precision> rv2;
        amp::campf<Precision> cv1;
        amp::campf<Precision> cv2;
        ap::template_2d_array< amp::ampf<Precision> > refra;
        ap::template_2d_array< amp::campf<Precision> > refca;
        ap::template_1d_array< amp::ampf<Precision> > rx;
        ap::template_1d_array< amp::ampf<Precision> > ry;
        ap::template_1d_array< amp::campf<Precision> > cx;
        ap::template_1d_array< amp::campf<Precision> > cy;
        int i_;
        int i1_;


        result = false;
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        for(mx=minn; mx<=maxn; mx++)
        {
            
            //
            // Select random M/N in [1,MX] such that max(M,N)=MX
            //
            m = 1+ap::randominteger(mx);
            n = 1+ap::randominteger(mx);
            if( ap::randominteger(2)==0 )
            {
                m = mx;
            }
            else
            {
                n = mx;
            }
            
            //
            // Initialize A by random matrix with size (MaxN+MaxN)*(MaxN+MaxN)
            // Initialize X by random vector with size (MaxN+MaxN)
            // Fill Y by control values
            //
            refra.setlength(maxn+maxn, maxn+maxn);
            refca.setlength(maxn+maxn, maxn+maxn);
            for(i=0; i<=2*maxn-1; i++)
            {
                for(j=0; j<=2*maxn-1; j++)
                {
                    refra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    refca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            rx.setlength(2*maxn);
            cx.setlength(2*maxn);
            ry.setlength(2*maxn);
            cy.setlength(2*maxn);
            for(i=0; i<=2*maxn-1; i++)
            {
                rx(i) = 2*amp::ampf<Precision>::getRandom()-1;
                cx(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                cx(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                ry(i) = i;
                cy(i) = i;
            }
            
            //
            // test different offsets (zero or one)
            //
            // to avoid unnecessary slowdown we don't test ALL possible
            // combinations of operation types. We just generate one random
            // set of parameters and test it.
            //
            aoffsi = ap::randominteger(maxn);
            aoffsj = ap::randominteger(maxn);
            xoffs = ap::randominteger(maxn);
            yoffs = ap::randominteger(maxn);
            opca = ap::randominteger(3);
            opra = ap::randominteger(2);
            ablas::cmatrixmv<Precision>(m, n, refca, aoffsi, aoffsj, opca, cx, xoffs, cy, yoffs);
            for(i=0; i<=2*maxn-1; i++)
            {
                if( i<yoffs || i>=yoffs+m )
                {
                    result = result || cy(i)!=i;
                }
                else
                {
                    cv1 = cy(i);
                    if( opca==0 )
                    {
                        i1_ = (xoffs)-(aoffsj);
                        cv2 = 0.0;
                        for(i_=aoffsj; i_<=aoffsj+n-1;i_++)
                        {
                            cv2 += refca(aoffsi+i-yoffs,i_)*cx(i_+i1_);
                        }
                    }
                    if( opca==1 )
                    {
                        i1_ = (xoffs)-(aoffsi);
                        cv2 = 0.0;
                        for(i_=aoffsi; i_<=aoffsi+n-1;i_++)
                        {
                            cv2 += refca(i_,aoffsj+i-yoffs)*cx(i_+i1_);
                        }
                    }
                    if( opca==2 )
                    {
                        i1_ = (xoffs)-(aoffsi);
                        cv2 = 0.0;
                        for(i_=aoffsi; i_<=aoffsi+n-1;i_++)
                        {
                            cv2 += amp::conj(refca(i_,aoffsj+i-yoffs))*cx(i_+i1_);
                        }
                    }
                    result = result || amp::abscomplex<Precision>(cv1-cv2)>threshold;
                }
            }
            ablas::rmatrixmv<Precision>(m, n, refra, aoffsi, aoffsj, opra, rx, xoffs, ry, yoffs);
            for(i=0; i<=2*maxn-1; i++)
            {
                if( i<yoffs || i>=yoffs+m )
                {
                    result = result || ry(i)!=i;
                }
                else
                {
                    rv1 = ry(i);
                    if( opra==0 )
                    {
                        rv2 = amp::vdotproduct(refra.getrow(aoffsi+i-yoffs, aoffsj, aoffsj+n-1), rx.getvector(xoffs, xoffs+n-1));
                    }
                    if( opra==1 )
                    {
                        rv2 = amp::vdotproduct(refra.getcolumn(aoffsj+i-yoffs, aoffsi, aoffsi+n-1), rx.getvector(xoffs, xoffs+n-1));
                    }
                    result = result || amp::abs<Precision>(rv1-rv2)>threshold;
                }
            }
        }
        return result;
    }


    /*************************************************************************
    COPY tests

    Returns False for passed test, True - for failed
    *************************************************************************/
    template<unsigned int Precision>
    bool testcopy(int minn,
        int maxn)
    {
        bool result;
        int m;
        int n;
        int mx;
        int i;
        int j;
        int aoffsi;
        int aoffsj;
        int boffsi;
        int boffsj;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> rv1;
        amp::ampf<Precision> rv2;
        amp::campf<Precision> cv1;
        amp::campf<Precision> cv2;
        ap::template_2d_array< amp::ampf<Precision> > ra;
        ap::template_2d_array< amp::ampf<Precision> > rb;
        ap::template_2d_array< amp::campf<Precision> > ca;
        ap::template_2d_array< amp::campf<Precision> > cb;


        result = false;
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        for(mx=minn; mx<=maxn; mx++)
        {
            
            //
            // Select random M/N in [1,MX] such that max(M,N)=MX
            //
            m = 1+ap::randominteger(mx);
            n = 1+ap::randominteger(mx);
            if( ap::randominteger(2)==0 )
            {
                m = mx;
            }
            else
            {
                n = mx;
            }
            
            //
            // Initialize A by random matrix with size (MaxN+MaxN)*(MaxN+MaxN)
            // Initialize X by random vector with size (MaxN+MaxN)
            // Fill Y by control values
            //
            ra.setlength(maxn+maxn, maxn+maxn);
            ca.setlength(maxn+maxn, maxn+maxn);
            rb.setlength(maxn+maxn, maxn+maxn);
            cb.setlength(maxn+maxn, maxn+maxn);
            for(i=0; i<=2*maxn-1; i++)
            {
                for(j=0; j<=2*maxn-1; j++)
                {
                    ra(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    ca(i,j).x = 2*amp::ampf<Precision>::getRandom()-1;
                    ca(i,j).y = 2*amp::ampf<Precision>::getRandom()-1;
                    rb(i,j) = 1+2*i+3*j;
                    cb(i,j) = 1+2*i+3*j;
                }
            }
            
            //
            // test different offsets (zero or one)
            //
            // to avoid unnecessary slowdown we don't test ALL possible
            // combinations of operation types. We just generate one random
            // set of parameters and test it.
            //
            aoffsi = ap::randominteger(maxn);
            aoffsj = ap::randominteger(maxn);
            boffsi = ap::randominteger(maxn);
            boffsj = ap::randominteger(maxn);
            ablas::cmatrixcopy<Precision>(m, n, ca, aoffsi, aoffsj, cb, boffsi, boffsj);
            for(i=0; i<=2*maxn-1; i++)
            {
                for(j=0; j<=2*maxn-1; j++)
                {
                    if( i<boffsi || i>=boffsi+m || j<boffsj || j>=boffsj+n )
                    {
                        result = result || cb(i,j)!=1+2*i+3*j;
                    }
                    else
                    {
                        result = result || amp::abscomplex<Precision>(ca(aoffsi+i-boffsi,aoffsj+j-boffsj)-cb(i,j))>threshold;
                    }
                }
            }
            ablas::rmatrixcopy<Precision>(m, n, ra, aoffsi, aoffsj, rb, boffsi, boffsj);
            for(i=0; i<=2*maxn-1; i++)
            {
                for(j=0; j<=2*maxn-1; j++)
                {
                    if( i<boffsi || i>=boffsi+m || j<boffsj || j>=boffsj+n )
                    {
                        result = result || rb(i,j)!=1+2*i+3*j;
                    }
                    else
                    {
                        result = result || amp::abs<Precision>(ra(aoffsi+i-boffsi,aoffsj+j-boffsj)-rb(i,j))>threshold;
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testablasunit_test_silent()
    {
        bool result;


        result = testablas<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testablasunit_test()
    {
        bool result;


        result = testablas<Precision>(false);
        return result;
    }
} // namespace

#endif
