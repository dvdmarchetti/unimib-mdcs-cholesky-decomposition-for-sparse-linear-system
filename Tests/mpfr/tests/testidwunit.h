
#ifndef _testidwunit_h
#define _testidwunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "tsort.h"
#include "nearestneighbor.h"
#include "reflections.h"
#include "hblas.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "hqrnd.h"
#include "matgen.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "xblas.h"
#include "densesolver.h"
#include "idwint.h"
namespace testidwunit
{
    template<unsigned int Precision>
    bool testidw(bool silent);
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void testxy(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        const int& n,
        const int& nx,
        const int& d,
        const int& nq,
        const int& nw,
        bool& idwerrors);
    template<unsigned int Precision>
    void testrxy(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        const int& n,
        const int& nx,
        const amp::ampf<Precision>& r,
        bool& idwerrors);
    template<unsigned int Precision>
    void testdegree(const int& n,
        const int& nx,
        const int& d,
        const int& dtask,
        bool& idwerrors);
    template<unsigned int Precision>
    void testnoisy(bool& idwerrors);
    template<unsigned int Precision>
    bool testidwunit_test_silent();
    template<unsigned int Precision>
    bool testidwunit_test();


    /*************************************************************************
    Testing IDW interpolation
    *************************************************************************/
    template<unsigned int Precision>
    bool testidw(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > xy;
        int i;
        int j;
        amp::ampf<Precision> vx;
        amp::ampf<Precision> vy;
        amp::ampf<Precision> vz;
        int d;
        int dtask;
        int nx;
        int n;
        int nq;
        int nw;
        int smalln;
        int largen;
        bool waserrors;
        bool idwerrors;


        idwerrors = false;
        smalln = 256;
        largen = 1024;
        nq = 10;
        nw = 18;
        
        //
        // Simple test:
        // * F = x^3 + sin(pi*y)*z^2 - (x+y)^2
        // * space is either R1=[-1,+1] (other dimensions are
        //   fixed at 0), R1^2 or R1^3.
        //* D = -1, 0, 1, 2
        //
        for(nx=1; nx<=2; nx++)
        {
            xy.setlength(largen, nx+1);
            for(i=0; i<=largen-1; i++)
            {
                for(j=0; j<=nx-1; j++)
                {
                    xy(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                if( nx>=1 )
                {
                    vx = xy(i,0);
                }
                else
                {
                    vx = 0;
                }
                if( nx>=2 )
                {
                    vy = xy(i,1);
                }
                else
                {
                    vy = 0;
                }
                if( nx>=3 )
                {
                    vz = xy(i,2);
                }
                else
                {
                    vz = 0;
                }
                xy(i,nx) = vx*vx*vx+amp::sin<Precision>(amp::pi<Precision>()*vy)*amp::sqr<Precision>(vz)-amp::sqr<Precision>(vx+vy);
            }
            for(d=-1; d<=2; d++)
            {
                testxy<Precision>(xy, largen, nx, d, nq, nw, idwerrors);
            }
        }
        
        //
        // Another simple test:
        // * five points in 2D - (0,0), (0,1), (1,0), (-1,0) (0,-1)
        // * F is random
        // * D = -1, 0, 1, 2
        //
        nx = 2;
        xy.setlength(5, nx+1);
        xy(0,0) = 0;
        xy(0,1) = 0;
        xy(0,2) = 2*amp::ampf<Precision>::getRandom()-1;
        xy(1,0) = 1;
        xy(1,1) = 0;
        xy(1,2) = 2*amp::ampf<Precision>::getRandom()-1;
        xy(2,0) = 0;
        xy(2,1) = 1;
        xy(2,2) = 2*amp::ampf<Precision>::getRandom()-1;
        xy(3,0) = -1;
        xy(3,1) = 0;
        xy(3,2) = 2*amp::ampf<Precision>::getRandom()-1;
        xy(4,0) = 0;
        xy(4,1) = -1;
        xy(4,2) = 2*amp::ampf<Precision>::getRandom()-1;
        for(d=-1; d<=2; d++)
        {
            testxy<Precision>(xy, 5, nx, d, nq, nw, idwerrors);
        }
        
        //
        // Degree test.
        //
        // F is either:
        // * constant (DTask=0)
        // * linear (DTask=1)
        // * quadratic (DTask=2)
        //
        // Nodal functions are either
        // * constant (D=0)
        // * linear (D=1)
        // * quadratic (D=2)
        //
        // When DTask<=D, we can interpolate without errors.
        // When DTask>D, we MUST have errors.
        //
        for(nx=1; nx<=3; nx++)
        {
            for(d=0; d<=2; d++)
            {
                for(dtask=0; dtask<=2; dtask++)
                {
                    testdegree<Precision>(smalln, nx, d, dtask, idwerrors);
                }
            }
        }
        
        //
        // Noisy test
        //
        testnoisy<Precision>(idwerrors);
        
        //
        // report
        //
        waserrors = idwerrors;
        if( !silent )
        {
            printf("TESTING INVERSE DISTANCE WEIGHTING\n");
            printf("* IDW:                                   ");
            if( !idwerrors )
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
    Unsets 2D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::campf<Precision> >& a)
    {
        a.setbounds(0, 0, 0, 0);
        a(0,0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Unsets 1D array.
    *************************************************************************/
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& a)
    {
        a.setbounds(0, 0);
        a(0) = 2*amp::ampf<Precision>::getRandom()-1;
    }


    /*************************************************************************
    Testing IDW:
    * generate model using N/NX/D/NQ/NW
    * test basic properties
    *************************************************************************/
    template<unsigned int Precision>
    void testxy(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        const int& n,
        const int& nx,
        const int& d,
        const int& nq,
        const int& nw,
        bool& idwerrors)
    {
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> lipschitzstep;
        int i;
        int j;
        int i1;
        int i2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> t;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;
        idwint::idwinterpolant<Precision> z1;
        ap::template_1d_array< amp::ampf<Precision> > x;


        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        lipschitzstep = amp::ampf<Precision>("0.001");
        x.setlength(nx);
        
        //
        // build
        //
        idwint::idwbuildmodifiedshepard<Precision>(xy, n, nx, d, nq, nw, z1);
        
        //
        // first, test interpolation properties at nodes
        //
        for(i=0; i<=n-1; i++)
        {
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
            idwerrors = idwerrors || idwint::idwcalc<Precision>(z1, x)!=xy(i,nx);
        }
        
        //
        // test Lipschitz continuity
        //
        i1 = ap::randominteger(n);
        do
        {
            i2 = ap::randominteger(n);
        }
        while( i2==i1 );
        l1 = 0;
        t = 0;
        while( t<1 )
        {
            v = 1-t;
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i1, 0, nx-1), v);
            v = t;
            amp::vadd(x.getvector(0, nx-1), xy.getrow(i2, 0, nx-1), v);
            v1 = idwint::idwcalc<Precision>(z1, x);
            v = 1-(t+lipschitzstep);
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i1, 0, nx-1), v);
            v = t+lipschitzstep;
            amp::vadd(x.getvector(0, nx-1), xy.getrow(i2, 0, nx-1), v);
            v2 = idwint::idwcalc<Precision>(z1, x);
            l1 = amp::maximum<Precision>(l1, amp::abs<Precision>(v2-v1)/lipschitzstep);
            t = t+lipschitzstep;
        }
        l2 = 0;
        t = 0;
        while( t<1 )
        {
            v = 1-t;
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i1, 0, nx-1), v);
            v = t;
            amp::vadd(x.getvector(0, nx-1), xy.getrow(i2, 0, nx-1), v);
            v1 = idwint::idwcalc<Precision>(z1, x);
            v = 1-(t+lipschitzstep/3);
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i1, 0, nx-1), v);
            v = t+lipschitzstep/3;
            amp::vadd(x.getvector(0, nx-1), xy.getrow(i2, 0, nx-1), v);
            v2 = idwint::idwcalc<Precision>(z1, x);
            l2 = amp::maximum<Precision>(l2, amp::abs<Precision>(v2-v1)/(lipschitzstep/3));
            t = t+lipschitzstep/3;
        }
        idwerrors = idwerrors || l2>amp::ampf<Precision>("2.0")*l1;
    }


    /*************************************************************************
    Testing IDW:
    * generate model using R-based model
    * test basic properties
    *************************************************************************/
    template<unsigned int Precision>
    void testrxy(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        const int& n,
        const int& nx,
        const amp::ampf<Precision>& r,
        bool& idwerrors)
    {
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> lipschitzstep;
        int i;
        int j;
        int i1;
        int i2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> t;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;
        idwint::idwinterpolant<Precision> z1;
        ap::template_1d_array< amp::ampf<Precision> > x;


        threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        lipschitzstep = amp::ampf<Precision>("0.001");
        x.setlength(nx);
        
        //
        // build
        //
        idwint::idwbuildmodifiedshepardr<Precision>(xy, n, nx, r, z1);
        
        //
        // first, test interpolation properties at nodes
        //
        for(i=0; i<=n-1; i++)
        {
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
            idwerrors = idwerrors || idwint::idwcalc<Precision>(z1, x)!=xy(i,nx);
        }
        
        //
        // test Lipschitz continuity
        //
        i1 = ap::randominteger(n);
        do
        {
            i2 = ap::randominteger(n);
        }
        while( i2==i1 );
        l1 = 0;
        t = 0;
        while( t<1 )
        {
            v = 1-t;
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i1, 0, nx-1), v);
            v = t;
            amp::vadd(x.getvector(0, nx-1), xy.getrow(i2, 0, nx-1), v);
            v1 = idwint::idwcalc<Precision>(z1, x);
            v = 1-(t+lipschitzstep);
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i1, 0, nx-1), v);
            v = t+lipschitzstep;
            amp::vadd(x.getvector(0, nx-1), xy.getrow(i2, 0, nx-1), v);
            v2 = idwint::idwcalc<Precision>(z1, x);
            l1 = amp::maximum<Precision>(l1, amp::abs<Precision>(v2-v1)/lipschitzstep);
            t = t+lipschitzstep;
        }
        l2 = 0;
        t = 0;
        while( t<1 )
        {
            v = 1-t;
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i1, 0, nx-1), v);
            v = t;
            amp::vadd(x.getvector(0, nx-1), xy.getrow(i2, 0, nx-1), v);
            v1 = idwint::idwcalc<Precision>(z1, x);
            v = 1-(t+lipschitzstep/3);
            amp::vmove(x.getvector(0, nx-1), xy.getrow(i1, 0, nx-1), v);
            v = t+lipschitzstep/3;
            amp::vadd(x.getvector(0, nx-1), xy.getrow(i2, 0, nx-1), v);
            v2 = idwint::idwcalc<Precision>(z1, x);
            l2 = amp::maximum<Precision>(l2, amp::abs<Precision>(v2-v1)/(lipschitzstep/3));
            t = t+lipschitzstep/3;
        }
        idwerrors = idwerrors || l2>amp::ampf<Precision>("2.0")*l1;
    }


    /*************************************************************************
    Testing degree properties

    F is either:
    * constant (DTask=0)
    * linear (DTask=1)
    * quadratic (DTask=2)

    Nodal functions are either
    * constant (D=0)
    * linear (D=1)
    * quadratic (D=2)

    When DTask<=D, we can interpolate without errors.
    When DTask>D, we MUST have errors.
    *************************************************************************/
    template<unsigned int Precision>
    void testdegree(const int& n,
        const int& nx,
        const int& d,
        const int& dtask,
        bool& idwerrors)
    {
        amp::ampf<Precision> threshold;
        int nq;
        int nw;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> c0;
        ap::template_1d_array< amp::ampf<Precision> > c1;
        ap::template_2d_array< amp::ampf<Precision> > c2;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_2d_array< amp::ampf<Precision> > xy;
        idwint::idwinterpolant<Precision> z1;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;


        threshold = amp::ampf<Precision>("1.0E6")*amp::ampf<Precision>::getAlgoPascalEpsilon();
        nq = 2*(nx*nx+nx+1);
        nw = 10;
        ap::ap_error::make_assertion(nq<=n);
        
        //
        // prepare model
        //
        c0 = 2*amp::ampf<Precision>::getRandom()-1;
        c1.setlength(nx);
        for(i=0; i<=nx-1; i++)
        {
            c1(i) = 2*amp::ampf<Precision>::getRandom()-1;
        }
        c2.setlength(nx, nx);
        for(i=0; i<=nx-1; i++)
        {
            for(j=i+1; j<=nx-1; j++)
            {
                c2(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                c2(j,i) = c2(i,j);
            }
            do
            {
                c2(i,i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            while( amp::abs<Precision>(c2(i,i))<=amp::ampf<Precision>("0.3") );
        }
        
        //
        // prepare points
        //
        xy.setlength(n, nx+1);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=nx-1; j++)
            {
                xy(i,j) = 4*amp::ampf<Precision>::getRandom()-2;
            }
            xy(i,nx) = c0;
            if( dtask>=1 )
            {
                v = amp::vdotproduct(c1.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
                xy(i,nx) = xy(i,nx)+v;
            }
            if( dtask==2 )
            {
                for(j=0; j<=nx-1; j++)
                {
                    v = amp::vdotproduct(c2.getrow(j, 0, nx-1), xy.getrow(i, 0, nx-1));
                    xy(i,nx) = xy(i,nx)+xy(i,j)*v;
                }
            }
        }
        
        //
        // build interpolant, calculate value at random point
        //
        idwint::idwbuildmodifiedshepard<Precision>(xy, n, nx, d, nq, nw, z1);
        x.setlength(nx);
        for(i=0; i<=nx-1; i++)
        {
            x(i) = 4*amp::ampf<Precision>::getRandom()-2;
        }
        v1 = idwint::idwcalc<Precision>(z1, x);
        
        //
        // calculate model value at the same point
        //
        v2 = c0;
        if( dtask>=1 )
        {
            v = amp::vdotproduct(c1.getvector(0, nx-1), x.getvector(0, nx-1));
            v2 = v2+v;
        }
        if( dtask==2 )
        {
            for(j=0; j<=nx-1; j++)
            {
                v = amp::vdotproduct(c2.getrow(j, 0, nx-1), x.getvector(0, nx-1));
                v2 = v2+x(j)*v;
            }
        }
        
        //
        // Compare
        //
        if( dtask<=d )
        {
            idwerrors = idwerrors || amp::abs<Precision>(v2-v1)>threshold;
        }
        else
        {
            idwerrors = idwerrors || amp::abs<Precision>(v2-v1)<threshold;
        }
    }


    /*************************************************************************
    Noisy test:
     * F = x^2 + y^2 + z^2 + noise on [-1,+1]^3
     * space is either R1=[-1,+1] (other dimensions are
       fixed at 0), R1^2 or R1^3.
     * D = 1, 2
     * 4096 points is used for function generation,
       4096 points - for testing
     * RMS error of "noisy" model on test set must be
       lower than RMS error of interpolation model.
    *************************************************************************/
    template<unsigned int Precision>
    void testnoisy(bool& idwerrors)
    {
        amp::ampf<Precision> noiselevel;
        int nq;
        int nw;
        int d;
        int nx;
        int ntrn;
        int ntst;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> t;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> ve;
        ap::template_2d_array< amp::ampf<Precision> > xy;
        ap::template_1d_array< amp::ampf<Precision> > x;
        idwint::idwinterpolant<Precision> z1;
        idwint::idwinterpolant<Precision> z2;
        amp::ampf<Precision> rms1;
        amp::ampf<Precision> rms2;


        nq = 20;
        nw = 40;
        noiselevel = amp::ampf<Precision>("0.2");
        ntrn = 256;
        ntst = 1024;
        for(d=1; d<=2; d++)
        {
            for(nx=1; nx<=2; nx++)
            {
                
                //
                // prepare dataset
                //
                xy.setlength(ntrn, nx+1);
                for(i=0; i<=ntrn-1; i++)
                {
                    v = noiselevel*(2*amp::ampf<Precision>::getRandom()-1);
                    for(j=0; j<=nx-1; j++)
                    {
                        t = 2*amp::ampf<Precision>::getRandom()-1;
                        v = v+amp::sqr<Precision>(t);
                        xy(i,j) = t;
                    }
                    xy(i,nx) = v;
                }
                
                //
                // build interpolants
                //
                idwint::idwbuildmodifiedshepard<Precision>(xy, ntrn, nx, d, nq, nw, z1);
                idwint::idwbuildnoisy<Precision>(xy, ntrn, nx, d, nq, nw, z2);
                
                //
                // calculate RMS errors
                //
                x.setlength(nx);
                rms1 = 0;
                rms2 = 0;
                for(i=0; i<=ntst-1; i++)
                {
                    ve = 0;
                    for(j=0; j<=nx-1; j++)
                    {
                        t = 2*amp::ampf<Precision>::getRandom()-1;
                        x(j) = t;
                        ve = ve+amp::sqr<Precision>(t);
                    }
                    v1 = idwint::idwcalc<Precision>(z1, x);
                    v2 = idwint::idwcalc<Precision>(z2, x);
                    rms1 = rms1+amp::sqr<Precision>(v1-ve);
                    rms2 = rms2+amp::sqr<Precision>(v2-ve);
                }
                idwerrors = idwerrors || rms2>rms1;
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testidwunit_test_silent()
    {
        bool result;


        result = testidw<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testidwunit_test()
    {
        bool result;


        result = testidw<Precision>(false);
        return result;
    }
} // namespace

#endif
