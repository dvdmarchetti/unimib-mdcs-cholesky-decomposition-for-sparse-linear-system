
#ifndef _testmatgenunit_h
#define _testmatgenunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
namespace testmatgenunit
{
    template<unsigned int Precision>
    bool testmatgen(bool silent);
    template<unsigned int Precision>
    bool isspd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        bool isupper);
    template<unsigned int Precision>
    bool obsoletesvddecomposition(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::ampf<Precision> >& v);
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void unset2dc(ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    bool ishpd(ap::template_2d_array< amp::campf<Precision> > a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> svdcond(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> extsign(amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    amp::ampf<Precision> mymax(amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    amp::ampf<Precision> pythag(amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    bool testmatgenunit_test_silent();
    template<unsigned int Precision>
    bool testmatgenunit_test();


    static const int maxsvditerations = 60;


    template<unsigned int Precision>
    bool testmatgen(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a;
        ap::template_2d_array< amp::ampf<Precision> > b;
        ap::template_2d_array< amp::ampf<Precision> > u;
        ap::template_2d_array< amp::ampf<Precision> > v;
        ap::template_2d_array< amp::campf<Precision> > ca;
        ap::template_2d_array< amp::campf<Precision> > cb;
        ap::template_2d_array< amp::ampf<Precision> > r1;
        ap::template_2d_array< amp::ampf<Precision> > r2;
        ap::template_2d_array< amp::campf<Precision> > c1;
        ap::template_2d_array< amp::campf<Precision> > c2;
        ap::template_1d_array< amp::ampf<Precision> > w;
        int n;
        int maxn;
        int i;
        int j;
        int pass;
        int passcount;
        bool waserrors;
        amp::ampf<Precision> cond;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> vt;
        amp::campf<Precision> ct;
        amp::ampf<Precision> minw;
        amp::ampf<Precision> maxw;
        bool serr;
        bool herr;
        bool spderr;
        bool hpderr;
        bool rerr;
        bool cerr;
        int i_;


        rerr = false;
        cerr = false;
        serr = false;
        herr = false;
        spderr = false;
        hpderr = false;
        waserrors = false;
        maxn = 20;
        passcount = 15;
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Testing orthogonal
        //
        for(n=1; n<=maxn; n++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                r1.setbounds(0, n-1, 0, 2*n-1);
                r2.setbounds(0, 2*n-1, 0, n-1);
                c1.setbounds(0, n-1, 0, 2*n-1);
                c2.setbounds(0, 2*n-1, 0, n-1);
                
                //
                // Random orthogonal, real
                //
                unset2d<Precision>(a);
                unset2d<Precision>(b);
                matgen::rmatrixrndorthogonal<Precision>(n, a);
                matgen::rmatrixrndorthogonal<Precision>(n, b);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        
                        //
                        // orthogonality test
                        //
                        vt = amp::vdotproduct(a.getrow(i, 0, n-1), a.getrow(j, 0, n-1));
                        if( i==j )
                        {
                            rerr = rerr || amp::abs<Precision>(vt-1)>threshold;
                        }
                        else
                        {
                            rerr = rerr || amp::abs<Precision>(vt)>threshold;
                        }
                        vt = amp::vdotproduct(b.getrow(i, 0, n-1), b.getrow(j, 0, n-1));
                        if( i==j )
                        {
                            rerr = rerr || amp::abs<Precision>(vt-1)>threshold;
                        }
                        else
                        {
                            rerr = rerr || amp::abs<Precision>(vt)>threshold;
                        }
                        
                        //
                        // test for difference in A and B
                        //
                        if( n>=2 )
                        {
                            rerr = rerr || a(i,j)==b(i,j);
                        }
                    }
                }
                
                //
                // Random orthogonal, complex
                //
                unset2dc<Precision>(ca);
                unset2dc<Precision>(cb);
                matgen::cmatrixrndorthogonal<Precision>(n, ca);
                matgen::cmatrixrndorthogonal<Precision>(n, cb);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        
                        //
                        // orthogonality test
                        //
                        ct = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            ct += ca(i,i_)*amp::conj(ca(j,i_));
                        }
                        if( i==j )
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct-1)>threshold;
                        }
                        else
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct)>threshold;
                        }
                        ct = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            ct += cb(i,i_)*amp::conj(cb(j,i_));
                        }
                        if( i==j )
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct-1)>threshold;
                        }
                        else
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct)>threshold;
                        }
                        
                        //
                        // test for difference in A and B
                        //
                        if( n>=2 )
                        {
                            cerr = cerr || ca(i,j)==cb(i,j);
                        }
                    }
                }
                
                //
                // From the right real tests:
                // 1. E*Q is orthogonal
                // 2. Q1<>Q2 (routine result is changing)
                // 3. (E E)'*Q = (Q' Q')' (correct handling of non-square matrices)
                //
                unset2d<Precision>(a);
                unset2d<Precision>(b);
                a.setbounds(0, n-1, 0, n-1);
                b.setbounds(0, n-1, 0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = 0;
                        b(i,j) = 0;
                    }
                    a(i,i) = 1;
                    b(i,i) = 1;
                }
                matgen::rmatrixrndorthogonalfromtheright<Precision>(a, n, n);
                matgen::rmatrixrndorthogonalfromtheright<Precision>(b, n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        
                        //
                        // orthogonality test
                        //
                        vt = amp::vdotproduct(a.getrow(i, 0, n-1), a.getrow(j, 0, n-1));
                        if( i==j )
                        {
                            rerr = rerr || amp::abs<Precision>(vt-1)>threshold;
                        }
                        else
                        {
                            rerr = rerr || amp::abs<Precision>(vt)>threshold;
                        }
                        vt = amp::vdotproduct(b.getrow(i, 0, n-1), b.getrow(j, 0, n-1));
                        if( i==j )
                        {
                            rerr = rerr || amp::abs<Precision>(vt-1)>threshold;
                        }
                        else
                        {
                            rerr = rerr || amp::abs<Precision>(vt)>threshold;
                        }
                        
                        //
                        // test for difference in A and B
                        //
                        if( n>=2 )
                        {
                            rerr = rerr || a(i,j)==b(i,j);
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        r2(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        r2(i+n,j) = r2(i,j);
                    }
                }
                matgen::rmatrixrndorthogonalfromtheright<Precision>(r2, 2*n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        rerr = rerr || amp::abs<Precision>(r2(i+n,j)-r2(i,j))>threshold;
                    }
                }
                
                //
                // From the left real tests:
                // 1. Q*E is orthogonal
                // 2. Q1<>Q2 (routine result is changing)
                // 3. Q*(E E) = (Q Q) (correct handling of non-square matrices)
                //
                unset2d<Precision>(a);
                unset2d<Precision>(b);
                a.setbounds(0, n-1, 0, n-1);
                b.setbounds(0, n-1, 0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = 0;
                        b(i,j) = 0;
                    }
                    a(i,i) = 1;
                    b(i,i) = 1;
                }
                matgen::rmatrixrndorthogonalfromtheleft<Precision>(a, n, n);
                matgen::rmatrixrndorthogonalfromtheleft<Precision>(b, n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        
                        //
                        // orthogonality test
                        //
                        vt = amp::vdotproduct(a.getrow(i, 0, n-1), a.getrow(j, 0, n-1));
                        if( i==j )
                        {
                            rerr = rerr || amp::abs<Precision>(vt-1)>threshold;
                        }
                        else
                        {
                            rerr = rerr || amp::abs<Precision>(vt)>threshold;
                        }
                        vt = amp::vdotproduct(b.getrow(i, 0, n-1), b.getrow(j, 0, n-1));
                        if( i==j )
                        {
                            rerr = rerr || amp::abs<Precision>(vt-1)>threshold;
                        }
                        else
                        {
                            rerr = rerr || amp::abs<Precision>(vt)>threshold;
                        }
                        
                        //
                        // test for difference in A and B
                        //
                        if( n>=2 )
                        {
                            rerr = rerr || a(i,j)==b(i,j);
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        r1(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        r1(i,j+n) = r1(i,j);
                    }
                }
                matgen::rmatrixrndorthogonalfromtheleft<Precision>(r1, n, 2*n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        rerr = rerr || amp::abs<Precision>(r1(i,j)-r1(i,j+n))>threshold;
                    }
                }
                
                //
                // From the right complex tests:
                // 1. E*Q is orthogonal
                // 2. Q1<>Q2 (routine result is changing)
                // 3. (E E)'*Q = (Q' Q')' (correct handling of non-square matrices)
                //
                unset2dc<Precision>(ca);
                unset2dc<Precision>(cb);
                ca.setbounds(0, n-1, 0, n-1);
                cb.setbounds(0, n-1, 0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        ca(i,j) = 0;
                        cb(i,j) = 0;
                    }
                    ca(i,i) = 1;
                    cb(i,i) = 1;
                }
                matgen::cmatrixrndorthogonalfromtheright<Precision>(ca, n, n);
                matgen::cmatrixrndorthogonalfromtheright<Precision>(cb, n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        
                        //
                        // orthogonality test
                        //
                        ct = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            ct += ca(i,i_)*amp::conj(ca(j,i_));
                        }
                        if( i==j )
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct-1)>threshold;
                        }
                        else
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct)>threshold;
                        }
                        ct = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            ct += cb(i,i_)*amp::conj(cb(j,i_));
                        }
                        if( i==j )
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct-1)>threshold;
                        }
                        else
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct)>threshold;
                        }
                        
                        //
                        // test for difference in A and B
                        //
                        cerr = cerr || ca(i,j)==cb(i,j);
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c2(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        c2(i+n,j) = c2(i,j);
                    }
                }
                matgen::cmatrixrndorthogonalfromtheright<Precision>(c2, 2*n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        cerr = cerr || amp::abscomplex<Precision>(c2(i+n,j)-c2(i,j))>threshold;
                    }
                }
                
                //
                // From the left complex tests:
                // 1. Q*E is orthogonal
                // 2. Q1<>Q2 (routine result is changing)
                // 3. Q*(E E) = (Q Q) (correct handling of non-square matrices)
                //
                unset2dc<Precision>(ca);
                unset2dc<Precision>(cb);
                ca.setbounds(0, n-1, 0, n-1);
                cb.setbounds(0, n-1, 0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        ca(i,j) = 0;
                        cb(i,j) = 0;
                    }
                    ca(i,i) = 1;
                    cb(i,i) = 1;
                }
                matgen::cmatrixrndorthogonalfromtheleft<Precision>(ca, n, n);
                matgen::cmatrixrndorthogonalfromtheleft<Precision>(cb, n, n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        
                        //
                        // orthogonality test
                        //
                        ct = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            ct += ca(i,i_)*amp::conj(ca(j,i_));
                        }
                        if( i==j )
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct-1)>threshold;
                        }
                        else
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct)>threshold;
                        }
                        ct = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            ct += cb(i,i_)*amp::conj(cb(j,i_));
                        }
                        if( i==j )
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct-1)>threshold;
                        }
                        else
                        {
                            cerr = cerr || amp::abscomplex<Precision>(ct)>threshold;
                        }
                        
                        //
                        // test for difference in A and B
                        //
                        cerr = cerr || ca(i,j)==cb(i,j);
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c1(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                        c1(i,j+n) = c1(i,j);
                    }
                }
                matgen::cmatrixrndorthogonalfromtheleft<Precision>(c1, n, 2*n);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        cerr = cerr || amp::abscomplex<Precision>(c1(i,j)-c1(i,j+n))>threshold;
                    }
                }
            }
        }
        
        //
        // Testing GCond
        //
        for(n=2; n<=maxn; n++)
        {
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // real test
                //
                unset2d<Precision>(a);
                cond = amp::exp<Precision>(amp::log<Precision>(amp::ampf<Precision>(1000))*amp::ampf<Precision>::getRandom());
                matgen::rmatrixrndcond<Precision>(n, cond, a);
                b.setbounds(1, n, 1, n);
                for(i=1; i<=n; i++)
                {
                    for(j=1; j<=n; j++)
                    {
                        b(i,j) = a(i-1,j-1);
                    }
                }
                if( obsoletesvddecomposition<Precision>(b, n, n, w, v) )
                {
                    maxw = w(1);
                    minw = w(1);
                    for(i=2; i<=n; i++)
                    {
                        if( w(i)>maxw )
                        {
                            maxw = w(i);
                        }
                        if( w(i)<minw )
                        {
                            minw = w(i);
                        }
                    }
                    vt = maxw/minw/cond;
                    if( amp::abs<Precision>(amp::log<Precision>(vt))>amp::log<Precision>(1+threshold) )
                    {
                        rerr = true;
                    }
                }
            }
        }
        
        //
        // Symmetric/SPD
        // N = 2 .. 30
        //
        for(n=2; n<=maxn; n++)
        {
            
            //
            // SPD matrices
            //
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // Generate A
                //
                unset2d<Precision>(a);
                cond = amp::exp<Precision>(amp::log<Precision>(amp::ampf<Precision>(1000))*amp::ampf<Precision>::getRandom());
                matgen::spdmatrixrndcond<Precision>(n, cond, a);
                
                //
                // test condition number
                //
                spderr = spderr || svdcond<Precision>(a, n)/cond-1>threshold;
                
                //
                // test SPD
                //
                spderr = spderr || !isspd<Precision>(a, n, true);
                
                //
                // test that A is symmetic
                //
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        spderr = spderr || amp::abs<Precision>(a(i,j)-a(j,i))>threshold;
                    }
                }
                
                //
                // test for difference between A and B (subsequent matrix)
                //
                unset2d<Precision>(b);
                matgen::spdmatrixrndcond<Precision>(n, cond, b);
                if( n>=2 )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            spderr = spderr || a(i,j)==b(i,j);
                        }
                    }
                }
            }
            
            //
            // HPD matrices
            //
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // Generate A
                //
                unset2dc<Precision>(ca);
                cond = amp::exp<Precision>(amp::log<Precision>(amp::ampf<Precision>(1000))*amp::ampf<Precision>::getRandom());
                matgen::hpdmatrixrndcond<Precision>(n, cond, ca);
                
                //
                // test HPD
                //
                hpderr = hpderr || !ishpd<Precision>(ca, n);
                
                //
                // test that A is Hermitian
                //
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        hpderr = hpderr || amp::abscomplex<Precision>(ca(i,j)-amp::conj<Precision>(ca(j,i)))>threshold;
                    }
                }
                
                //
                // test for difference between A and B (subsequent matrix)
                //
                unset2dc<Precision>(cb);
                matgen::hpdmatrixrndcond<Precision>(n, cond, cb);
                if( n>=2 )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            hpderr = hpderr || ca(i,j)==cb(i,j);
                        }
                    }
                }
            }
            
            //
            // Symmetric matrices
            //
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // test condition number
                //
                unset2d<Precision>(a);
                cond = amp::exp<Precision>(amp::log<Precision>(amp::ampf<Precision>(1000))*amp::ampf<Precision>::getRandom());
                matgen::smatrixrndcond<Precision>(n, cond, a);
                serr = serr || svdcond<Precision>(a, n)/cond-1>threshold;
                
                //
                // test for difference between A and B
                //
                unset2d<Precision>(b);
                matgen::smatrixrndcond<Precision>(n, cond, b);
                if( n>=2 )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            serr = serr || a(i,j)==b(i,j);
                        }
                    }
                }
            }
            
            //
            // Hermitian matrices
            //
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // Generate A
                //
                unset2dc<Precision>(ca);
                cond = amp::exp<Precision>(amp::log<Precision>(amp::ampf<Precision>(1000))*amp::ampf<Precision>::getRandom());
                matgen::hmatrixrndcond<Precision>(n, cond, ca);
                
                //
                // test that A is Hermitian
                //
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        herr = herr || amp::abscomplex<Precision>(ca(i,j)-amp::conj<Precision>(ca(j,i)))>threshold;
                    }
                }
                
                //
                // test for difference between A and B (subsequent matrix)
                //
                unset2dc<Precision>(cb);
                matgen::hmatrixrndcond<Precision>(n, cond, cb);
                if( n>=2 )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            herr = herr || ca(i,j)==cb(i,j);
                        }
                    }
                }
            }
        }
        
        //
        // report
        //
        waserrors = rerr || cerr || serr || spderr || herr || hpderr;
        if( !silent )
        {
            printf("TESTING MATRIX GENERATOR\n");
            printf("REAL TEST:                               ");
            if( !rerr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("COMPLEX TEST:                            ");
            if( !cerr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("SYMMETRIC TEST:                          ");
            if( !serr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("HERMITIAN TEST:                          ");
            if( !herr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("SPD TEST:                                ");
            if( !spderr )
            {
                printf("OK\n");
            }
            else
            {
                printf("FAILED\n");
            }
            printf("HPD TEST:                                ");
            if( !hpderr )
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
    Test whether matrix is SPD
    *************************************************************************/
    template<unsigned int Precision>
    bool isspd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        bool isupper)
    {
        bool result;
        int i;
        int j;
        amp::ampf<Precision> ajj;
        amp::ampf<Precision> v;


        
        //
        //     Test the input parameters.
        //
        ap::ap_error::make_assertion(n>=0);
        
        //
        //     Quick return if possible
        //
        result = true;
        if( n<=0 )
        {
            return result;
        }
        if( isupper )
        {
            
            //
            // Compute the Cholesky factorization A = U'*U.
            //
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Compute U(J,J) and test for non-positive-definiteness.
                //
                v = amp::vdotproduct(a.getcolumn(j, 0, j-1), a.getcolumn(j, 0, j-1));
                ajj = a(j,j)-v;
                if( ajj<=0 )
                {
                    result = false;
                    return result;
                }
                ajj = amp::sqrt<Precision>(ajj);
                a(j,j) = ajj;
                
                //
                // Compute elements J+1:N of row J.
                //
                if( j<n-1 )
                {
                    for(i=j+1; i<=n-1; i++)
                    {
                        v = amp::vdotproduct(a.getcolumn(i, 0, j-1), a.getcolumn(j, 0, j-1));
                        a(j,i) = a(j,i)-v;
                    }
                    v = 1/ajj;
                    amp::vmul(a.getrow(j, j+1, n-1), v);
                }
            }
        }
        else
        {
            
            //
            // Compute the Cholesky factorization A = L*L'.
            //
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Compute L(J,J) and test for non-positive-definiteness.
                //
                v = amp::vdotproduct(a.getrow(j, 0, j-1), a.getrow(j, 0, j-1));
                ajj = a(j,j)-v;
                if( ajj<=0 )
                {
                    result = false;
                    return result;
                }
                ajj = amp::sqrt<Precision>(ajj);
                a(j,j) = ajj;
                
                //
                // Compute elements J+1:N of column J.
                //
                if( j<n-1 )
                {
                    for(i=j+1; i<=n-1; i++)
                    {
                        v = amp::vdotproduct(a.getrow(i, 0, j-1), a.getrow(j, 0, j-1));
                        a(i,j) = a(i,j)-v;
                    }
                    v = 1/ajj;
                    amp::vmul(a.getcolumn(j, j+1, n-1), v);
                }
            }
        }
        return result;
    }


    template<unsigned int Precision>
    bool obsoletesvddecomposition(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::ampf<Precision> >& v)
    {
        bool result;
        int nm;
        int minmn;
        int l;
        int k;
        int j;
        int jj;
        int its;
        int i;
        amp::ampf<Precision> z;
        amp::ampf<Precision> y;
        amp::ampf<Precision> x;
        amp::ampf<Precision> vscale;
        amp::ampf<Precision> s;
        amp::ampf<Precision> h;
        amp::ampf<Precision> g;
        amp::ampf<Precision> f;
        amp::ampf<Precision> c;
        amp::ampf<Precision> anorm;
        ap::template_1d_array< amp::ampf<Precision> > rv1;
        bool flag;


        rv1.setbounds(1, n);
        w.setbounds(1, n);
        v.setbounds(1, n, 1, n);
        result = true;
        if( m<n )
        {
            minmn = m;
        }
        else
        {
            minmn = n;
        }
        g = amp::ampf<Precision>("0.0");
        vscale = amp::ampf<Precision>("0.0");
        anorm = amp::ampf<Precision>("0.0");
        for(i=1; i<=n; i++)
        {
            l = i+1;
            rv1(i) = vscale*g;
            g = 0;
            s = 0;
            vscale = 0;
            if( i<=m )
            {
                for(k=i; k<=m; k++)
                {
                    vscale = vscale+amp::abs<Precision>(a(k,i));
                }
                if( vscale!=amp::ampf<Precision>("0.0") )
                {
                    for(k=i; k<=m; k++)
                    {
                        a(k,i) = a(k,i)/vscale;
                        s = s+a(k,i)*a(k,i);
                    }
                    f = a(i,i);
                    g = -extsign<Precision>(amp::sqrt<Precision>(s), f);
                    h = f*g-s;
                    a(i,i) = f-g;
                    if( i!=n )
                    {
                        for(j=l; j<=n; j++)
                        {
                            s = amp::ampf<Precision>("0.0");
                            for(k=i; k<=m; k++)
                            {
                                s = s+a(k,i)*a(k,j);
                            }
                            f = s/h;
                            for(k=i; k<=m; k++)
                            {
                                a(k,j) = a(k,j)+f*a(k,i);
                            }
                        }
                    }
                    for(k=i; k<=m; k++)
                    {
                        a(k,i) = vscale*a(k,i);
                    }
                }
            }
            w(i) = vscale*g;
            g = amp::ampf<Precision>("0.0");
            s = amp::ampf<Precision>("0.0");
            vscale = amp::ampf<Precision>("0.0");
            if( i<=m && i!=n )
            {
                for(k=l; k<=n; k++)
                {
                    vscale = vscale+amp::abs<Precision>(a(i,k));
                }
                if( vscale!=amp::ampf<Precision>("0.0") )
                {
                    for(k=l; k<=n; k++)
                    {
                        a(i,k) = a(i,k)/vscale;
                        s = s+a(i,k)*a(i,k);
                    }
                    f = a(i,l);
                    g = -extsign<Precision>(amp::sqrt<Precision>(s), f);
                    h = f*g-s;
                    a(i,l) = f-g;
                    for(k=l; k<=n; k++)
                    {
                        rv1(k) = a(i,k)/h;
                    }
                    if( i!=m )
                    {
                        for(j=l; j<=m; j++)
                        {
                            s = amp::ampf<Precision>("0.0");
                            for(k=l; k<=n; k++)
                            {
                                s = s+a(j,k)*a(i,k);
                            }
                            for(k=l; k<=n; k++)
                            {
                                a(j,k) = a(j,k)+s*rv1(k);
                            }
                        }
                    }
                    for(k=l; k<=n; k++)
                    {
                        a(i,k) = vscale*a(i,k);
                    }
                }
            }
            anorm = mymax<Precision>(anorm, amp::abs<Precision>(w(i))+amp::abs<Precision>(rv1(i)));
        }
        for(i=n; i>=1; i--)
        {
            if( i<n )
            {
                if( g!=amp::ampf<Precision>("0.0") )
                {
                    for(j=l; j<=n; j++)
                    {
                        v(j,i) = a(i,j)/a(i,l)/g;
                    }
                    for(j=l; j<=n; j++)
                    {
                        s = amp::ampf<Precision>("0.0");
                        for(k=l; k<=n; k++)
                        {
                            s = s+a(i,k)*v(k,j);
                        }
                        for(k=l; k<=n; k++)
                        {
                            v(k,j) = v(k,j)+s*v(k,i);
                        }
                    }
                }
                for(j=l; j<=n; j++)
                {
                    v(i,j) = amp::ampf<Precision>("0.0");
                    v(j,i) = amp::ampf<Precision>("0.0");
                }
            }
            v(i,i) = amp::ampf<Precision>("1.0");
            g = rv1(i);
            l = i;
        }
        for(i=minmn; i>=1; i--)
        {
            l = i+1;
            g = w(i);
            if( i<n )
            {
                for(j=l; j<=n; j++)
                {
                    a(i,j) = amp::ampf<Precision>("0.0");
                }
            }
            if( g!=amp::ampf<Precision>("0.0") )
            {
                g = amp::ampf<Precision>("1.0")/g;
                if( i!=n )
                {
                    for(j=l; j<=n; j++)
                    {
                        s = amp::ampf<Precision>("0.0");
                        for(k=l; k<=m; k++)
                        {
                            s = s+a(k,i)*a(k,j);
                        }
                        f = s/a(i,i)*g;
                        for(k=i; k<=m; k++)
                        {
                            a(k,j) = a(k,j)+f*a(k,i);
                        }
                    }
                }
                for(j=i; j<=m; j++)
                {
                    a(j,i) = a(j,i)*g;
                }
            }
            else
            {
                for(j=i; j<=m; j++)
                {
                    a(j,i) = amp::ampf<Precision>("0.0");
                }
            }
            a(i,i) = a(i,i)+amp::ampf<Precision>("1.0");
        }
        for(k=n; k>=1; k--)
        {
            for(its=1; its<=maxsvditerations; its++)
            {
                flag = true;
                for(l=k; l>=1; l--)
                {
                    nm = l-1;
                    if( amp::abs<Precision>(rv1(l))+anorm==anorm )
                    {
                        flag = false;
                        break;
                    }
                    if( amp::abs<Precision>(w(nm))+anorm==anorm )
                    {
                        break;
                    }
                }
                if( flag )
                {
                    c = amp::ampf<Precision>("0.0");
                    s = amp::ampf<Precision>("1.0");
                    for(i=l; i<=k; i++)
                    {
                        f = s*rv1(i);
                        if( amp::abs<Precision>(f)+anorm!=anorm )
                        {
                            g = w(i);
                            h = pythag<Precision>(f, g);
                            w(i) = h;
                            h = amp::ampf<Precision>("1.0")/h;
                            c = g*h;
                            s = -f*h;
                            for(j=1; j<=m; j++)
                            {
                                y = a(j,nm);
                                z = a(j,i);
                                a(j,nm) = y*c+z*s;
                                a(j,i) = -y*s+z*c;
                            }
                        }
                    }
                }
                z = w(k);
                if( l==k )
                {
                    if( z<amp::ampf<Precision>("0.0") )
                    {
                        w(k) = -z;
                        for(j=1; j<=n; j++)
                        {
                            v(j,k) = -v(j,k);
                        }
                    }
                    break;
                }
                if( its==maxsvditerations )
                {
                    result = false;
                    return result;
                }
                x = w(l);
                nm = k-1;
                y = w(nm);
                g = rv1(nm);
                h = rv1(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(amp::ampf<Precision>("2.0")*h*y);
                g = pythag<Precision>(f, amp::ampf<Precision>(1));
                f = ((x-z)*(x+z)+h*(y/(f+extsign<Precision>(g, f))-h))/x;
                c = amp::ampf<Precision>("1.0");
                s = amp::ampf<Precision>("1.0");
                for(j=l; j<=nm; j++)
                {
                    i = j+1;
                    g = rv1(i);
                    y = w(i);
                    h = s*g;
                    g = c*g;
                    z = pythag<Precision>(f, h);
                    rv1(j) = z;
                    c = f/z;
                    s = h/z;
                    f = x*c+g*s;
                    g = -x*s+g*c;
                    h = y*s;
                    y = y*c;
                    for(jj=1; jj<=n; jj++)
                    {
                        x = v(jj,j);
                        z = v(jj,i);
                        v(jj,j) = x*c+z*s;
                        v(jj,i) = -x*s+z*c;
                    }
                    z = pythag<Precision>(f, h);
                    w(j) = z;
                    if( z!=amp::ampf<Precision>("0.0") )
                    {
                        z = amp::ampf<Precision>("1.0")/z;
                        c = f*z;
                        s = h*z;
                    }
                    f = c*g+s*y;
                    x = -s*g+c*y;
                    for(jj=1; jj<=m; jj++)
                    {
                        y = a(jj,j);
                        z = a(jj,i);
                        a(jj,j) = y*c+z*s;
                        a(jj,i) = -y*s+z*c;
                    }
                }
                rv1(l) = amp::ampf<Precision>("0.0");
                rv1(k) = f;
                w(k) = x;
            }
        }
        return result;
    }


    /*************************************************************************
    Unsets 2D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::ampf<Precision> >& a)
    {
        a.setbounds(0, 0, 0, 0);
        a(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets 2D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset2dc(ap::template_2d_array< amp::campf<Precision> >& a)
    {
        a.setbounds(0, 0, 0, 0);
        a(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Tests whether A is HPD
    *************************************************************************/
    template<unsigned int Precision>
    bool ishpd(ap::template_2d_array< amp::campf<Precision> > a,
        int n)
    {
        bool result;
        int j;
        amp::ampf<Precision> ajj;
        amp::campf<Precision> v;
        amp::ampf<Precision> r;
        ap::template_1d_array< amp::campf<Precision> > t;
        ap::template_1d_array< amp::campf<Precision> > t2;
        ap::template_1d_array< amp::campf<Precision> > t3;
        int i;
        ap::template_2d_array< amp::campf<Precision> > a1;
        int i_;


        t.setbounds(0, n-1);
        t2.setbounds(0, n-1);
        t3.setbounds(0, n-1);
        result = true;
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        for(j=0; j<=n-1; j++)
        {
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            v = 0.0;
            for(i_=0; i_<=j-1;i_++)
            {
                v += amp::conj(a(i_,j))*a(i_,j);
            }
            ajj = (a(j,j)-v).x;
            if( ajj<=0 )
            {
                a(j,j) = ajj;
                result = false;
                return result;
            }
            ajj = amp::sqrt<Precision>(ajj);
            a(j,j) = ajj;
            
            //
            // Compute elements J+1:N-1 of row J.
            //
            if( j<n-1 )
            {
                for(i_=0; i_<=j-1;i_++)
                {
                    t2(i_) = amp::conj(a(i_,j));
                }
                for(i_=j+1; i_<=n-1;i_++)
                {
                    t3(i_) = a(j,i_);
                }
                for(i=j+1; i<=n-1; i++)
                {
                    v = 0.0;
                    for(i_=0; i_<=j-1;i_++)
                    {
                        v += a(i_,i)*t2(i_);
                    }
                    t3(i) = t3(i)-v;
                }
                for(i_=j+1; i_<=n-1;i_++)
                {
                    a(j,i_) = t3(i_);
                }
                r = 1/ajj;
                for(i_=j+1; i_<=n-1;i_++)
                {
                    a(j,i_) = r*a(j,i_);
                }
            }
        }
        return result;
    }


    /*************************************************************************
    SVD condition number
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> svdcond(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n)
    {
        amp::ampf<Precision> result;
        ap::template_2d_array< amp::ampf<Precision> > a1;
        ap::template_2d_array< amp::ampf<Precision> > v;
        ap::template_1d_array< amp::ampf<Precision> > w;
        int i;
        int j;
        amp::ampf<Precision> minw;
        amp::ampf<Precision> maxw;


        a1.setbounds(1, n, 1, n);
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                a1(i,j) = a(i-1,j-1);
            }
        }
        if( !obsoletesvddecomposition<Precision>(a1, n, n, w, v) )
        {
            result = 0;
            return result;
        }
        minw = w(1);
        maxw = w(1);
        for(i=2; i<=n; i++)
        {
            if( w(i)<minw )
            {
                minw = w(i);
            }
            if( w(i)>maxw )
            {
                maxw = w(i);
            }
        }
        result = maxw/minw;
        return result;
    }


    template<unsigned int Precision>
    amp::ampf<Precision> extsign(amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        amp::ampf<Precision> result;


        if( b>=0 )
        {
            result = amp::abs<Precision>(a);
        }
        else
        {
            result = -amp::abs<Precision>(a);
        }
        return result;
    }


    template<unsigned int Precision>
    amp::ampf<Precision> mymax(amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        amp::ampf<Precision> result;


        if( a>b )
        {
            result = a;
        }
        else
        {
            result = b;
        }
        return result;
    }


    template<unsigned int Precision>
    amp::ampf<Precision> pythag(amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        amp::ampf<Precision> result;


        if( amp::abs<Precision>(a)<amp::abs<Precision>(b) )
        {
            result = amp::abs<Precision>(b)*amp::sqrt<Precision>(1+amp::sqr<Precision>(a/b));
        }
        else
        {
            result = amp::abs<Precision>(a)*amp::sqrt<Precision>(1+amp::sqr<Precision>(b/a));
        }
        return result;
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testmatgenunit_test_silent()
    {
        bool result;


        result = testmatgen<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testmatgenunit_test()
    {
        bool result;


        result = testmatgen<Precision>(false);
        return result;
    }
} // namespace

#endif
