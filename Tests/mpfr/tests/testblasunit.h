
#ifndef _testblasunit_h
#define _testblasunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "blas.h"
namespace testblasunit
{
    template<unsigned int Precision>
    bool testblas(bool silent);
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
    bool testblasunit_test_silent();
    template<unsigned int Precision>
    bool testblasunit_test();


    template<unsigned int Precision>
    bool testblas(bool silent)
    {
        bool result;
        int pass;
        int passcount;
        int n;
        int i;
        int i1;
        int i2;
        int j;
        int j1;
        int j2;
        int l;
        int k;
        int r;
        int i3;
        int j3;
        int col1;
        int col2;
        int row1;
        int row2;
        ap::template_1d_array< amp::ampf<Precision> > x1;
        ap::template_1d_array< amp::ampf<Precision> > x2;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_2d_array< amp::ampf<Precision> > c1;
        ap::template_2d_array< amp::ampf<Precision> > c2;
        amp::ampf<Precision> err;
        amp::ampf<Precision> e1;
        amp::ampf<Precision> e2;
        amp::ampf<Precision> e3;
        amp::ampf<Precision> v;
        amp::ampf<Precision> scl1;
        amp::ampf<Precision> scl2;
        amp::ampf<Precision> scl3;
        bool was1;
        bool was2;
        bool trans1;
        bool trans2;
        amp::ampf<Precision> threshold;
        bool n2errors;
        bool hsnerrors;
        bool amaxerrors;
        bool mverrors;
        bool iterrors;
        bool cterrors;
        bool mmerrors;
        bool waserrors;


        n2errors = false;
        amaxerrors = false;
        hsnerrors = false;
        mverrors = false;
        iterrors = false;
        cterrors = false;
        mmerrors = false;
        waserrors = false;
        threshold = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Test Norm2
        //
        passcount = 1000;
        e1 = 0;
        e2 = 0;
        e3 = 0;
        scl2 = amp::ampf<Precision>("0.5")*amp::ampf<Precision>::getAlgoPascalMaxNumber();
        scl3 = 2*amp::ampf<Precision>::getAlgoPascalMinNumber();
        for(pass=1; pass<=passcount; pass++)
        {
            n = 1+ap::randominteger(1000);
            i1 = ap::randominteger(10);
            i2 = n+i1-1;
            x1.setbounds(i1, i2);
            x2.setbounds(i1, i2);
            for(i=i1; i<=i2; i++)
            {
                x1(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            v = 0;
            for(i=i1; i<=i2; i++)
            {
                v = v+amp::sqr<Precision>(x1(i));
            }
            v = amp::sqrt<Precision>(v);
            e1 = amp::maximum<Precision>(e1, amp::abs<Precision>(v-blas::vectornorm2<Precision>(x1, i1, i2)));
            for(i=i1; i<=i2; i++)
            {
                x2(i) = scl2*x1(i);
            }
            e2 = amp::maximum<Precision>(e2, amp::abs<Precision>(v*scl2-blas::vectornorm2<Precision>(x2, i1, i2)));
            for(i=i1; i<=i2; i++)
            {
                x2(i) = scl3*x1(i);
            }
            e3 = amp::maximum<Precision>(e3, amp::abs<Precision>(v*scl3-blas::vectornorm2<Precision>(x2, i1, i2)));
        }
        e2 = e2/scl2;
        e3 = e3/scl3;
        n2errors = e1>=threshold || e2>=threshold || e3>=threshold;
        
        //
        // Testing VectorAbsMax, Column/Row AbsMax
        //
        x1.setbounds(1, 5);
        x1(1) = amp::ampf<Precision>("2.0");
        x1(2) = amp::ampf<Precision>("0.2");
        x1(3) = -amp::ampf<Precision>("1.3");
        x1(4) = amp::ampf<Precision>("0.7");
        x1(5) = -amp::ampf<Precision>("3.0");
        amaxerrors = blas::vectoridxabsmax<Precision>(x1, 1, 5)!=5 || blas::vectoridxabsmax<Precision>(x1, 1, 4)!=1 || blas::vectoridxabsmax<Precision>(x1, 2, 4)!=3;
        n = 30;
        x1.setbounds(1, n);
        a.setbounds(1, n, 1, n);
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
            }
        }
        was1 = false;
        was2 = false;
        for(pass=1; pass<=1000; pass++)
        {
            j = 1+ap::randominteger(n);
            i1 = 1+ap::randominteger(n);
            i2 = i1+ap::randominteger(n+1-i1);
            amp::vmove(x1.getvector(i1, i2), a.getcolumn(j, i1, i2));
            if( blas::vectoridxabsmax<Precision>(x1, i1, i2)!=blas::columnidxabsmax<Precision>(a, i1, i2, j) )
            {
                was1 = true;
            }
            i = 1+ap::randominteger(n);
            j1 = 1+ap::randominteger(n);
            j2 = j1+ap::randominteger(n+1-j1);
            amp::vmove(x1.getvector(j1, j2), a.getrow(i, j1, j2));
            if( blas::vectoridxabsmax<Precision>(x1, j1, j2)!=blas::rowidxabsmax<Precision>(a, j1, j2, i) )
            {
                was2 = true;
            }
        }
        amaxerrors = amaxerrors || was1 || was2;
        
        //
        // Testing upper Hessenberg 1-norm
        //
        a.setbounds(1, 3, 1, 3);
        x1.setbounds(1, 3);
        a(1,1) = 2;
        a(1,2) = 3;
        a(1,3) = 1;
        a(2,1) = 4;
        a(2,2) = -5;
        a(2,3) = 8;
        a(3,1) = 99;
        a(3,2) = 3;
        a(3,3) = 1;
        hsnerrors = amp::abs<Precision>(blas::upperhessenberg1norm<Precision>(a, 1, 3, 1, 3, x1)-11)>threshold;
        
        //
        // Testing MatrixVectorMultiply
        //
        a.setbounds(2, 3, 3, 5);
        x1.setbounds(1, 3);
        x2.setbounds(1, 2);
        a(2,3) = 2;
        a(2,4) = -1;
        a(2,5) = -1;
        a(3,3) = 1;
        a(3,4) = -2;
        a(3,5) = 2;
        x1(1) = 1;
        x1(2) = 2;
        x1(3) = 1;
        x2(1) = -1;
        x2(2) = -1;
        blas::matrixvectormultiply<Precision>(a, 2, 3, 3, 5, false, x1, 1, 3, amp::ampf<Precision>("1.0"), x2, 1, 2, amp::ampf<Precision>("1.0"));
        blas::matrixvectormultiply<Precision>(a, 2, 3, 3, 5, true, x2, 1, 2, amp::ampf<Precision>("1.0"), x1, 1, 3, amp::ampf<Precision>("1.0"));
        e1 = amp::abs<Precision>(x1(1)+5)+amp::abs<Precision>(x1(2)-8)+amp::abs<Precision>(x1(3)+1)+amp::abs<Precision>(x2(1)+2)+amp::abs<Precision>(x2(2)+2);
        x1(1) = 1;
        x1(2) = 2;
        x1(3) = 1;
        x2(1) = -1;
        x2(2) = -1;
        blas::matrixvectormultiply<Precision>(a, 2, 3, 3, 5, false, x1, 1, 3, amp::ampf<Precision>("1.0"), x2, 1, 2, amp::ampf<Precision>("0.0"));
        blas::matrixvectormultiply<Precision>(a, 2, 3, 3, 5, true, x2, 1, 2, amp::ampf<Precision>("1.0"), x1, 1, 3, amp::ampf<Precision>("0.0"));
        e2 = amp::abs<Precision>(x1(1)+3)+amp::abs<Precision>(x1(2)-3)+amp::abs<Precision>(x1(3)+1)+amp::abs<Precision>(x2(1)+1)+amp::abs<Precision>(x2(2)+1);
        mverrors = e1+e2>=threshold;
        
        //
        // testing inplace transpose
        //
        n = 10;
        a.setbounds(1, n, 1, n);
        b.setbounds(1, n, 1, n);
        x1.setbounds(1, n-1);
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                a(i,j) = amp::ampf<Precision>::getRandom();
            }
        }
        passcount = 10000;
        was1 = false;
        for(pass=1; pass<=passcount; pass++)
        {
            i1 = 1+ap::randominteger(n);
            i2 = i1+ap::randominteger(n-i1+1);
            j1 = 1+ap::randominteger(n-(i2-i1));
            j2 = j1+(i2-i1);
            blas::copymatrix<Precision>(a, i1, i2, j1, j2, b, i1, i2, j1, j2);
            blas::inplacetranspose<Precision>(b, i1, i2, j1, j2, x1);
            for(i=i1; i<=i2; i++)
            {
                for(j=j1; j<=j2; j++)
                {
                    if( a(i,j)!=b(i1+(j-j1),j1+(i-i1)) )
                    {
                        was1 = true;
                    }
                }
            }
        }
        iterrors = was1;
        
        //
        // testing copy and transpose
        //
        n = 10;
        a.setbounds(1, n, 1, n);
        b.setbounds(1, n, 1, n);
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                a(i,j) = amp::ampf<Precision>::getRandom();
            }
        }
        passcount = 10000;
        was1 = false;
        for(pass=1; pass<=passcount; pass++)
        {
            i1 = 1+ap::randominteger(n);
            i2 = i1+ap::randominteger(n-i1+1);
            j1 = 1+ap::randominteger(n);
            j2 = j1+ap::randominteger(n-j1+1);
            blas::copyandtranspose<Precision>(a, i1, i2, j1, j2, b, j1, j2, i1, i2);
            for(i=i1; i<=i2; i++)
            {
                for(j=j1; j<=j2; j++)
                {
                    if( a(i,j)!=b(j,i) )
                    {
                        was1 = true;
                    }
                }
            }
        }
        cterrors = was1;
        
        //
        // Testing MatrixMatrixMultiply
        //
        n = 10;
        a.setbounds(1, 2*n, 1, 2*n);
        b.setbounds(1, 2*n, 1, 2*n);
        c1.setbounds(1, 2*n, 1, 2*n);
        c2.setbounds(1, 2*n, 1, 2*n);
        x1.setbounds(1, n);
        x2.setbounds(1, n);
        for(i=1; i<=2*n; i++)
        {
            for(j=1; j<=2*n; j++)
            {
                a(i,j) = amp::ampf<Precision>::getRandom();
                b(i,j) = amp::ampf<Precision>::getRandom();
            }
        }
        passcount = 1000;
        was1 = false;
        for(pass=1; pass<=passcount; pass++)
        {
            for(i=1; i<=2*n; i++)
            {
                for(j=1; j<=2*n; j++)
                {
                    c1(i,j) = amp::ampf<Precision>("2.1")*i+amp::ampf<Precision>("3.1")*j;
                    c2(i,j) = c1(i,j);
                }
            }
            l = 1+ap::randominteger(n);
            k = 1+ap::randominteger(n);
            r = 1+ap::randominteger(n);
            i1 = 1+ap::randominteger(n);
            j1 = 1+ap::randominteger(n);
            i2 = 1+ap::randominteger(n);
            j2 = 1+ap::randominteger(n);
            i3 = 1+ap::randominteger(n);
            j3 = 1+ap::randominteger(n);
            trans1 = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
            trans2 = amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5");
            if( trans1 )
            {
                col1 = l;
                row1 = k;
            }
            else
            {
                col1 = k;
                row1 = l;
            }
            if( trans2 )
            {
                col2 = k;
                row2 = r;
            }
            else
            {
                col2 = r;
                row2 = k;
            }
            scl1 = amp::ampf<Precision>::getRandom();
            scl2 = amp::ampf<Precision>::getRandom();
            blas::matrixmatrixmultiply<Precision>(a, i1, i1+row1-1, j1, j1+col1-1, trans1, b, i2, i2+row2-1, j2, j2+col2-1, trans2, scl1, c1, i3, i3+l-1, j3, j3+r-1, scl2, x1);
            naivematrixmatrixmultiply<Precision>(a, i1, i1+row1-1, j1, j1+col1-1, trans1, b, i2, i2+row2-1, j2, j2+col2-1, trans2, scl1, c2, i3, i3+l-1, j3, j3+r-1, scl2);
            err = 0;
            for(i=1; i<=l; i++)
            {
                for(j=1; j<=r; j++)
                {
                    err = amp::maximum<Precision>(err, amp::abs<Precision>(c1(i3+i-1,j3+j-1)-c2(i3+i-1,j3+j-1)));
                }
            }
            if( err>threshold )
            {
                was1 = true;
                break;
            }
        }
        mmerrors = was1;
        
        //
        // report
        //
        waserrors = n2errors || amaxerrors || hsnerrors || mverrors || iterrors || cterrors || mmerrors;
        if( !silent )
        {
            printf("TESTING BLAS\n");
            printf("VectorNorm2:                             ");
            if( n2errors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("AbsMax (vector/row/column):              ");
            if( amaxerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("UpperHessenberg1Norm:                    ");
            if( hsnerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("MatrixVectorMultiply:                    ");
            if( mverrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("InplaceTranspose:                        ");
            if( iterrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("CopyAndTranspose:                        ");
            if( cterrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("MatrixMatrixMultiply:                    ");
            if( mmerrors )
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
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testblasunit_test_silent()
    {
        bool result;


        result = testblas<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testblasunit_test()
    {
        bool result;


        result = testblas<Precision>(false);
        return result;
    }
} // namespace

#endif
