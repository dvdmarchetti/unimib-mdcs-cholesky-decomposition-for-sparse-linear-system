
#ifndef _testxblasunit_h
#define _testxblasunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "xblas.h"
namespace testxblasunit
{
    template<unsigned int Precision>
    bool testxblas(bool silent);
    template<unsigned int Precision>
    bool testxblasunit_test_silent();
    template<unsigned int Precision>
    bool testxblasunit_test();


    template<unsigned int Precision>
    const amp::ampf<Precision>& xchunk()
    {
        static amp::ampf<Precision> v = 1048576;
        return v;
    }
    static const int xchunkcount = 4;


    template<unsigned int Precision>
    bool testxblas(bool silent)
    {
        bool result;
        bool approxerrors;
        bool exactnesserrors;
        bool waserrors;
        amp::ampf<Precision> approxthreshold;
        int maxn;
        int passcount;
        int n;
        int i;
        int pass;
        amp::ampf<Precision> rv1;
        amp::ampf<Precision> rv2;
        amp::ampf<Precision> rv2err;
        amp::campf<Precision> cv1;
        amp::campf<Precision> cv2;
        amp::ampf<Precision> cv2err;
        amp::ampf<Precision> cv2errx;
        amp::ampf<Precision> cv2erry;
        ap::template_1d_array< amp::ampf<Precision> > rx;
        ap::template_1d_array< amp::ampf<Precision> > ry;
        ap::template_1d_array< amp::campf<Precision> > cx;
        ap::template_1d_array< amp::campf<Precision> > cy;
        ap::template_1d_array< amp::ampf<Precision> > temp;
        amp::ampf<Precision> b;
        amp::ampf<Precision> s;
        int i_;


        approxerrors = false;
        exactnesserrors = false;
        waserrors = false;
        approxthreshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        maxn = 1000;
        passcount = 10;
        
        //
        // tests:
        // 1. ability to calculate dot product
        // 2. higher precision
        //
        for(n=1; n<=maxn; n++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                //  ability to approximately calculate real dot product
                //
                rx.setlength(n);
                ry.setlength(n);
                temp.setlength(n);
                for(i=0; i<=n-1; i++)
                {
                    if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.2") )
                    {
                        rx(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    else
                    {
                        rx(i) = 0;
                    }
                    if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.2") )
                    {
                        ry(i) = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    else
                    {
                        ry(i) = 0;
                    }
                }
                rv1 = amp::vdotproduct(rx.getvector(0, n-1), ry.getvector(0, n-1));
                xblas::xdot<Precision>(rx, ry, n, temp, rv2, rv2err);
                approxerrors = approxerrors || amp::abs<Precision>(rv1-rv2)>approxthreshold;
                
                //
                //  ability to approximately calculate complex dot product
                //
                cx.setlength(n);
                cy.setlength(n);
                temp.setlength(2*n);
                for(i=0; i<=n-1; i++)
                {
                    if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.2") )
                    {
                        cx(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                        cx(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    else
                    {
                        cx(i) = 0;
                    }
                    if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.2") )
                    {
                        cy(i).x = 2*amp::ampf<Precision>::getRandom()-1;
                        cy(i).y = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    else
                    {
                        cy(i) = 0;
                    }
                }
                cv1 = 0.0;
                for(i_=0; i_<=n-1;i_++)
                {
                    cv1 += cx(i_)*cy(i_);
                }
                xblas::xcdot<Precision>(cx, cy, n, temp, cv2, cv2err);
                approxerrors = approxerrors || amp::abscomplex<Precision>(cv1-cv2)>approxthreshold;
            }
        }
        
        //
        // test of precision: real
        //
        n = 50000;
        rx.setlength(n);
        ry.setlength(n);
        temp.setlength(n);
        for(pass=0; pass<=passcount-1; pass++)
        {
            ap::ap_error::make_assertion(n%2==0);
            
            //
            // First test: X + X + ... + X - X - X - ... - X = 1*X
            //
            s = amp::exp<Precision>(amp::ampf<Precision>(ap::maxint(pass, 50)));
            if( pass==passcount-1 && pass>1 )
            {
                s = amp::ampf<Precision>::getAlgoPascalMaxNumber();
            }
            ry(0) = (2*amp::ampf<Precision>::getRandom()-1)*s*amp::sqrt<Precision>(2*amp::ampf<Precision>::getRandom());
            for(i=1; i<=n-1; i++)
            {
                ry(i) = ry(0);
            }
            for(i=0; i<=n/2-1; i++)
            {
                rx(i) = 1;
            }
            for(i=n/2; i<=n-2; i++)
            {
                rx(i) = -1;
            }
            rx(n-1) = 0;
            xblas::xdot<Precision>(rx, ry, n, temp, rv2, rv2err);
            exactnesserrors = exactnesserrors || rv2err<0;
            exactnesserrors = exactnesserrors || rv2err>4*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::abs<Precision>(ry(0));
            exactnesserrors = exactnesserrors || amp::abs<Precision>(rv2-ry(0))>rv2err;
            
            //
            // First test: X + X + ... + X = N*X
            //
            s = amp::exp<Precision>(amp::ampf<Precision>(ap::maxint(pass, 50)));
            if( pass==passcount-1 && pass>1 )
            {
                s = amp::ampf<Precision>::getAlgoPascalMaxNumber();
            }
            ry(0) = (2*amp::ampf<Precision>::getRandom()-1)*s*amp::sqrt<Precision>(2*amp::ampf<Precision>::getRandom());
            for(i=1; i<=n-1; i++)
            {
                ry(i) = ry(0);
            }
            for(i=0; i<=n-1; i++)
            {
                rx(i) = 1;
            }
            xblas::xdot<Precision>(rx, ry, n, temp, rv2, rv2err);
            exactnesserrors = exactnesserrors || rv2err<0;
            exactnesserrors = exactnesserrors || rv2err>4*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::abs<Precision>(ry(0))*n;
            exactnesserrors = exactnesserrors || amp::abs<Precision>(rv2-n*ry(0))>rv2err;
        }
        
        //
        // test of precision: complex
        //
        n = 50000;
        cx.setlength(n);
        cy.setlength(n);
        temp.setlength(2*n);
        for(pass=0; pass<=passcount-1; pass++)
        {
            ap::ap_error::make_assertion(n%2==0);
            
            //
            // First test: X + X + ... + X - X - X - ... - X = 1*X
            //
            s = amp::exp<Precision>(amp::ampf<Precision>(ap::maxint(pass, 50)));
            if( pass==passcount-1 && pass>1 )
            {
                s = amp::ampf<Precision>::getAlgoPascalMaxNumber();
            }
            cy(0).x = (2*amp::ampf<Precision>::getRandom()-1)*s*amp::sqrt<Precision>(2*amp::ampf<Precision>::getRandom());
            cy(0).y = (2*amp::ampf<Precision>::getRandom()-1)*s*amp::sqrt<Precision>(2*amp::ampf<Precision>::getRandom());
            for(i=1; i<=n-1; i++)
            {
                cy(i) = cy(0);
            }
            for(i=0; i<=n/2-1; i++)
            {
                cx(i) = 1;
            }
            for(i=n/2; i<=n-2; i++)
            {
                cx(i) = -1;
            }
            cx(n-1) = 0;
            xblas::xcdot<Precision>(cx, cy, n, temp, cv2, cv2err);
            exactnesserrors = exactnesserrors || cv2err<0;
            exactnesserrors = exactnesserrors || cv2err>4*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::abscomplex<Precision>(cy(0));
            exactnesserrors = exactnesserrors || amp::abscomplex<Precision>(cv2-cy(0))>cv2err;
            
            //
            // First test: X + X + ... + X = N*X
            //
            s = amp::exp<Precision>(amp::ampf<Precision>(ap::maxint(pass, 50)));
            if( pass==passcount-1 && pass>1 )
            {
                s = amp::ampf<Precision>::getAlgoPascalMaxNumber();
            }
            cy(0) = (2*amp::ampf<Precision>::getRandom()-1)*s*amp::sqrt<Precision>(2*amp::ampf<Precision>::getRandom());
            for(i=1; i<=n-1; i++)
            {
                cy(i) = cy(0);
            }
            for(i=0; i<=n-1; i++)
            {
                cx(i) = 1;
            }
            xblas::xcdot<Precision>(cx, cy, n, temp, cv2, cv2err);
            exactnesserrors = exactnesserrors || cv2err<0;
            exactnesserrors = exactnesserrors || cv2err>4*amp::ampf<Precision>::getAlgoPascalEpsilon()*amp::abscomplex<Precision>(cy(0))*n;
            exactnesserrors = exactnesserrors || amp::abscomplex<Precision>(cv2-n*cy(0))>cv2err;
        }
        
        //
        // report
        //
        waserrors = approxerrors || exactnesserrors;
        if( !silent )
        {
            printf("TESTING XBLAS\n");
            printf("APPROX.TESTS:                            ");
            if( approxerrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("EXACT TESTS:                             ");
            if( exactnesserrors )
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
        
        //
        // end
        //
        result = !waserrors;
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testxblasunit_test_silent()
    {
        bool result;


        result = testxblas<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testxblasunit_test()
    {
        bool result;


        result = testxblas<Precision>(false);
        return result;
    }
} // namespace

#endif
