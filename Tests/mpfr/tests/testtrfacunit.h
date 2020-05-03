
#ifndef _testtrfacunit_h
#define _testtrfacunit_h

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
namespace testtrfacunit
{
    template<unsigned int Precision>
    bool testtrfac(bool silent);
    template<unsigned int Precision>
    void testcluproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& err,
        bool& properr);
    template<unsigned int Precision>
    void testrluproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& err,
        bool& properr);
    template<unsigned int Precision>
    bool testtrfacunit_test_silent();
    template<unsigned int Precision>
    bool testtrfacunit_test();


    template<unsigned int Precision>
    bool testtrfac(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > ra;
        ap::template_2d_array< amp::ampf<Precision> > ral;
        ap::template_2d_array< amp::ampf<Precision> > rau;
        ap::template_2d_array< amp::campf<Precision> > ca;
        ap::template_2d_array< amp::campf<Precision> > cal;
        ap::template_2d_array< amp::campf<Precision> > cau;
        int m;
        int n;
        int mx;
        int maxmn;
        int i;
        int j;
        int minij;
        int pass;
        amp::campf<Precision> vc;
        amp::ampf<Precision> vr;
        bool waserrors;
        bool spderr;
        bool hpderr;
        bool rerr;
        bool cerr;
        bool properr;
        amp::ampf<Precision> threshold;
        int i_;


        rerr = false;
        spderr = false;
        cerr = false;
        hpderr = false;
        properr = false;
        waserrors = false;
        maxmn = 4*ablas::ablasblocksize<Precision>(ra)+1;
        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon()*maxmn;
        
        //
        // test LU
        //
        for(mx=1; mx<=maxmn; mx++)
        {
            
            //
            // Initialize N/M, both are <=MX,
            // at least one of them is exactly equal to MX
            //
            n = 1+ap::randominteger(mx);
            m = 1+ap::randominteger(mx);
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                n = mx;
            }
            else
            {
                m = mx;
            }
            
            //
            // First, test on zero matrix
            //
            ra.setlength(m, n);
            ca.setlength(m, n);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    ra(i,j) = 0;
                    ca(i,j) = 0;
                }
            }
            testcluproblem<Precision>(ca, m, n, threshold, cerr, properr);
            testrluproblem<Precision>(ra, m, n, threshold, rerr, properr);
            
            //
            // Second, random matrix with moderate condition number
            //
            ra.setlength(m, n);
            ca.setlength(m, n);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    ra(i,j) = 0;
                    ca(i,j) = 0;
                }
            }
            for(i=0; i<=ap::minint(m, n)-1; i++)
            {
                ra(i,i) = 1+10*amp::ampf<Precision>::getRandom();
                ca(i,i) = 1+10*amp::ampf<Precision>::getRandom();
            }
            matgen::cmatrixrndorthogonalfromtheleft<Precision>(ca, m, n);
            matgen::cmatrixrndorthogonalfromtheright<Precision>(ca, m, n);
            matgen::rmatrixrndorthogonalfromtheleft<Precision>(ra, m, n);
            matgen::rmatrixrndorthogonalfromtheright<Precision>(ra, m, n);
            testcluproblem<Precision>(ca, m, n, threshold, cerr, properr);
            testrluproblem<Precision>(ra, m, n, threshold, rerr, properr);
        }
        
        //
        // Test Cholesky
        //
        for(n=1; n<=maxmn; n++)
        {
            
            //
            // Load CA (HPD matrix with low condition number),
            //      CAL and CAU - its lower and upper triangles
            //
            matgen::hpdmatrixrndcond<Precision>(n, 1+50*amp::ampf<Precision>::getRandom(), ca);
            cal.setlength(n, n);
            cau.setlength(n, n);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    cal(i,j) = i;
                    cau(i,j) = j;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(i_=0; i_<=i;i_++)
                {
                    cal(i,i_) = ca(i,i_);
                }
                for(i_=i; i_<=n-1;i_++)
                {
                    cau(i,i_) = ca(i,i_);
                }
            }
            
            //
            // Test HPDMatrixCholesky:
            // 1. it must leave upper (lower) part unchanged
            // 2. max(A-L*L^H) must be small
            //
            if( trfac::hpdmatrixcholesky<Precision>(cal, n, false) )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        if( j>i )
                        {
                            hpderr = hpderr || cal(i,j)!=i;
                        }
                        else
                        {
                            vc = 0.0;
                            for(i_=0; i_<=j;i_++)
                            {
                                vc += cal(i,i_)*amp::conj(cal(j,i_));
                            }
                            hpderr = hpderr || amp::abscomplex<Precision>(ca(i,j)-vc)>threshold;
                        }
                    }
                }
            }
            else
            {
                hpderr = true;
            }
            if( trfac::hpdmatrixcholesky<Precision>(cau, n, true) )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        if( j<i )
                        {
                            hpderr = hpderr || cau(i,j)!=j;
                        }
                        else
                        {
                            vc = 0.0;
                            for(i_=0; i_<=i;i_++)
                            {
                                vc += amp::conj(cau(i_,i))*cau(i_,j);
                            }
                            hpderr = hpderr || amp::abscomplex<Precision>(ca(i,j)-vc)>threshold;
                        }
                    }
                }
            }
            else
            {
                hpderr = true;
            }
            
            //
            // Load RA (SPD matrix with low condition number),
            //      RAL and RAU - its lower and upper triangles
            //
            matgen::spdmatrixrndcond<Precision>(n, 1+50*amp::ampf<Precision>::getRandom(), ra);
            ral.setlength(n, n);
            rau.setlength(n, n);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    ral(i,j) = i;
                    rau(i,j) = j;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(ral.getrow(i, 0, i), ra.getrow(i, 0, i));
                amp::vmove(rau.getrow(i, i, n-1), ra.getrow(i, i, n-1));
            }
            
            //
            // Test SPDMatrixCholesky:
            // 1. it must leave upper (lower) part unchanged
            // 2. max(A-L*L^H) must be small
            //
            if( trfac::spdmatrixcholesky<Precision>(ral, n, false) )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        if( j>i )
                        {
                            spderr = spderr || ral(i,j)!=i;
                        }
                        else
                        {
                            vr = amp::vdotproduct(ral.getrow(i, 0, j), ral.getrow(j, 0, j));
                            spderr = spderr || amp::abs<Precision>(ra(i,j)-vr)>threshold;
                        }
                    }
                }
            }
            else
            {
                spderr = true;
            }
            if( trfac::spdmatrixcholesky<Precision>(rau, n, true) )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        if( j<i )
                        {
                            spderr = spderr || rau(i,j)!=j;
                        }
                        else
                        {
                            vr = amp::vdotproduct(rau.getcolumn(i, 0, i), rau.getcolumn(j, 0, i));
                            spderr = spderr || amp::abs<Precision>(ra(i,j)-vr)>threshold;
                        }
                    }
                }
            }
            else
            {
                spderr = true;
            }
        }
        
        //
        // report
        //
        waserrors = rerr || spderr || cerr || hpderr || properr;
        if( !silent )
        {
            printf("TESTING TRIANGULAR FACTORIZATIONS\n");
            printf("* REAL:                                  ");
            if( rerr )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* SPD:                                   ");
            if( spderr )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* COMPLEX:                               ");
            if( cerr )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* HPD:                                   ");
            if( hpderr )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("* OTHER PROPERTIES:                      ");
            if( properr )
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
    void testcluproblem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& err,
        bool& properr)
    {
        ap::template_2d_array< amp::campf<Precision> > ca;
        ap::template_2d_array< amp::campf<Precision> > cl;
        ap::template_2d_array< amp::campf<Precision> > cu;
        ap::template_2d_array< amp::campf<Precision> > ca2;
        ap::template_1d_array< amp::campf<Precision> > ct;
        int i;
        int j;
        int minmn;
        amp::campf<Precision> v;
        ap::template_1d_array< int > p;
        int i_;


        minmn = ap::minint(m, n);
        
        //
        // PLU test
        //
        ca.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            for(i_=0; i_<=n-1;i_++)
            {
                ca(i,i_) = a(i,i_);
            }
        }
        trfac::cmatrixplu<Precision>(ca, m, n, p);
        for(i=0; i<=minmn-1; i++)
        {
            if( p(i)<i || p(i)>=m )
            {
                properr = false;
                return;
            }
        }
        cl.setlength(m, minmn);
        for(j=0; j<=minmn-1; j++)
        {
            for(i=0; i<=j-1; i++)
            {
                cl(i,j) = amp::ampf<Precision>("0.0");
            }
            cl(j,j) = amp::ampf<Precision>("1.0");
            for(i=j+1; i<=m-1; i++)
            {
                cl(i,j) = ca(i,j);
            }
        }
        cu.setlength(minmn, n);
        for(i=0; i<=minmn-1; i++)
        {
            for(j=0; j<=i-1; j++)
            {
                cu(i,j) = amp::ampf<Precision>("0.0");
            }
            for(j=i; j<=n-1; j++)
            {
                cu(i,j) = ca(i,j);
            }
        }
        ca2.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=minmn-1;i_++)
                {
                    v += cl(i,i_)*cu(i_,j);
                }
                ca2(i,j) = v;
            }
        }
        ct.setlength(n);
        for(i=minmn-1; i>=0; i--)
        {
            if( i!=p(i) )
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    ct(i_) = ca2(i,i_);
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    ca2(i,i_) = ca2(p(i),i_);
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    ca2(p(i),i_) = ct(i_);
                }
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                err = err || amp::abscomplex<Precision>(a(i,j)-ca2(i,j))>threshold;
            }
        }
        
        //
        // LUP test
        //
        ca.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            for(i_=0; i_<=n-1;i_++)
            {
                ca(i,i_) = a(i,i_);
            }
        }
        trfac::cmatrixlup<Precision>(ca, m, n, p);
        for(i=0; i<=minmn-1; i++)
        {
            if( p(i)<i || p(i)>=n )
            {
                properr = false;
                return;
            }
        }
        cl.setlength(m, minmn);
        for(j=0; j<=minmn-1; j++)
        {
            for(i=0; i<=j-1; i++)
            {
                cl(i,j) = amp::ampf<Precision>("0.0");
            }
            for(i=j; i<=m-1; i++)
            {
                cl(i,j) = ca(i,j);
            }
        }
        cu.setlength(minmn, n);
        for(i=0; i<=minmn-1; i++)
        {
            for(j=0; j<=i-1; j++)
            {
                cu(i,j) = amp::ampf<Precision>("0.0");
            }
            cu(i,i) = amp::ampf<Precision>("1.0");
            for(j=i+1; j<=n-1; j++)
            {
                cu(i,j) = ca(i,j);
            }
        }
        ca2.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = 0.0;
                for(i_=0; i_<=minmn-1;i_++)
                {
                    v += cl(i,i_)*cu(i_,j);
                }
                ca2(i,j) = v;
            }
        }
        ct.setlength(m);
        for(i=minmn-1; i>=0; i--)
        {
            if( i!=p(i) )
            {
                for(i_=0; i_<=m-1;i_++)
                {
                    ct(i_) = ca2(i_,i);
                }
                for(i_=0; i_<=m-1;i_++)
                {
                    ca2(i_,i) = ca2(i_,p(i));
                }
                for(i_=0; i_<=m-1;i_++)
                {
                    ca2(i_,p(i)) = ct(i_);
                }
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                err = err || amp::abscomplex<Precision>(a(i,j)-ca2(i,j))>threshold;
            }
        }
    }


    template<unsigned int Precision>
    void testrluproblem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        amp::ampf<Precision> threshold,
        bool& err,
        bool& properr)
    {
        ap::template_2d_array< amp::ampf<Precision> > ca;
        ap::template_2d_array< amp::ampf<Precision> > cl;
        ap::template_2d_array< amp::ampf<Precision> > cu;
        ap::template_2d_array< amp::ampf<Precision> > ca2;
        ap::template_1d_array< amp::ampf<Precision> > ct;
        int i;
        int j;
        int minmn;
        amp::ampf<Precision> v;
        ap::template_1d_array< int > p;


        minmn = ap::minint(m, n);
        
        //
        // PLU test
        //
        ca.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            amp::vmove(ca.getrow(i, 0, n-1), a.getrow(i, 0, n-1));
        }
        trfac::rmatrixplu<Precision>(ca, m, n, p);
        for(i=0; i<=minmn-1; i++)
        {
            if( p(i)<i || p(i)>=m )
            {
                properr = false;
                return;
            }
        }
        cl.setlength(m, minmn);
        for(j=0; j<=minmn-1; j++)
        {
            for(i=0; i<=j-1; i++)
            {
                cl(i,j) = amp::ampf<Precision>("0.0");
            }
            cl(j,j) = amp::ampf<Precision>("1.0");
            for(i=j+1; i<=m-1; i++)
            {
                cl(i,j) = ca(i,j);
            }
        }
        cu.setlength(minmn, n);
        for(i=0; i<=minmn-1; i++)
        {
            for(j=0; j<=i-1; j++)
            {
                cu(i,j) = amp::ampf<Precision>("0.0");
            }
            for(j=i; j<=n-1; j++)
            {
                cu(i,j) = ca(i,j);
            }
        }
        ca2.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(cl.getrow(i, 0, minmn-1), cu.getcolumn(j, 0, minmn-1));
                ca2(i,j) = v;
            }
        }
        ct.setlength(n);
        for(i=minmn-1; i>=0; i--)
        {
            if( i!=p(i) )
            {
                amp::vmove(ct.getvector(0, n-1), ca2.getrow(i, 0, n-1));
                amp::vmove(ca2.getrow(i, 0, n-1), ca2.getrow(p(i), 0, n-1));
                amp::vmove(ca2.getrow(p(i), 0, n-1), ct.getvector(0, n-1));
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                err = err || amp::abs<Precision>(a(i,j)-ca2(i,j))>threshold;
            }
        }
        
        //
        // LUP test
        //
        ca.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            amp::vmove(ca.getrow(i, 0, n-1), a.getrow(i, 0, n-1));
        }
        trfac::rmatrixlup<Precision>(ca, m, n, p);
        for(i=0; i<=minmn-1; i++)
        {
            if( p(i)<i || p(i)>=n )
            {
                properr = false;
                return;
            }
        }
        cl.setlength(m, minmn);
        for(j=0; j<=minmn-1; j++)
        {
            for(i=0; i<=j-1; i++)
            {
                cl(i,j) = amp::ampf<Precision>("0.0");
            }
            for(i=j; i<=m-1; i++)
            {
                cl(i,j) = ca(i,j);
            }
        }
        cu.setlength(minmn, n);
        for(i=0; i<=minmn-1; i++)
        {
            for(j=0; j<=i-1; j++)
            {
                cu(i,j) = amp::ampf<Precision>("0.0");
            }
            cu(i,i) = amp::ampf<Precision>("1.0");
            for(j=i+1; j<=n-1; j++)
            {
                cu(i,j) = ca(i,j);
            }
        }
        ca2.setlength(m, n);
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = amp::vdotproduct(cl.getrow(i, 0, minmn-1), cu.getcolumn(j, 0, minmn-1));
                ca2(i,j) = v;
            }
        }
        ct.setlength(m);
        for(i=minmn-1; i>=0; i--)
        {
            if( i!=p(i) )
            {
                amp::vmove(ct.getvector(0, m-1), ca2.getcolumn(i, 0, m-1));
                amp::vmove(ca2.getcolumn(i, 0, m-1), ca2.getcolumn(p(i), 0, m-1));
                amp::vmove(ca2.getcolumn(p(i), 0, m-1), ct.getvector(0, m-1));
            }
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                err = err || amp::abs<Precision>(a(i,j)-ca2(i,j))>threshold;
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testtrfacunit_test_silent()
    {
        bool result;


        result = testtrfac<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testtrfacunit_test()
    {
        bool result;


        result = testtrfac<Precision>(false);
        return result;
    }
} // namespace

#endif
