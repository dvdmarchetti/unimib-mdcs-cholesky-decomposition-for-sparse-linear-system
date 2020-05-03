
#ifndef _testpsplineunit_h
#define _testpsplineunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "spline3.h"
#include "blas.h"
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
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "linmin.h"
#include "minlbfgs.h"
#include "minlm.h"
#include "lsfit.h"
#include "apserv.h"
#include "spline1d.h"
#include "tsort.h"
#include "hsschur.h"
#include "evd.h"
#include "gammafunc.h"
#include "gq.h"
#include "gkq.h"
#include "autogk.h"
#include "pspline.h"
namespace testpsplineunit
{
    template<unsigned int Precision>
    bool testpsplineinterpolation(bool silent);
    template<unsigned int Precision>
    void unsetp2(pspline::pspline2interpolant<Precision>& p);
    template<unsigned int Precision>
    void unsetp3(pspline::pspline3interpolant<Precision>& p);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    bool testpsplineunit_test_silent();
    template<unsigned int Precision>
    bool testpsplineunit_test();


    template<unsigned int Precision>
    bool testpsplineinterpolation(bool silent)
    {
        bool result;
        bool waserrors;
        bool p2errors;
        bool p3errors;
        amp::ampf<Precision> nonstrictthreshold;
        amp::ampf<Precision> threshold;
        int passcount;
        amp::ampf<Precision> lstep;
        amp::ampf<Precision> h;
        int maxn;
        int periodicity;
        int skind;
        int pkind;
        bool periodic;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        int n;
        int tmpn;
        int i;
        int k;
        amp::ampf<Precision> vx;
        amp::ampf<Precision> vy;
        amp::ampf<Precision> vz;
        amp::ampf<Precision> vx2;
        amp::ampf<Precision> vy2;
        amp::ampf<Precision> vz2;
        amp::ampf<Precision> vdx;
        amp::ampf<Precision> vdy;
        amp::ampf<Precision> vdz;
        amp::ampf<Precision> vdx2;
        amp::ampf<Precision> vdy2;
        amp::ampf<Precision> vdz2;
        amp::ampf<Precision> vd2x;
        amp::ampf<Precision> vd2y;
        amp::ampf<Precision> vd2z;
        amp::ampf<Precision> vd2x2;
        amp::ampf<Precision> vd2y2;
        amp::ampf<Precision> vd2z2;
        amp::ampf<Precision> v0;
        amp::ampf<Precision> v1;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > z;
        ap::template_1d_array< amp::ampf<Precision> > t;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        ap::template_1d_array< amp::ampf<Precision> > t3;
        ap::template_2d_array< amp::ampf<Precision> > xy;
        ap::template_2d_array< amp::ampf<Precision> > xyz;
        pspline::pspline2interpolant<Precision> p2;
        pspline::pspline3interpolant<Precision> p3;
        spline1d::spline1dinterpolant<Precision> s;


        waserrors = false;
        passcount = 20;
        lstep = amp::ampf<Precision>("0.005");
        h = amp::ampf<Precision>("0.00001");
        maxn = 10;
        threshold = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        nonstrictthreshold = amp::ampf<Precision>("0.00001");
        p2errors = false;
        p3errors = false;
        
        //
        // Test basic properties of 2- and 3-dimensional splines:
        // * PSpline2ParameterValues() properties
        // * values at nodes
        // * for periodic splines - periodicity properties
        //
        // Variables used:
        // * N              points count
        // * SKind          spline
        // * PKind          parameterization
        // * Periodicity    whether we have periodic spline or not
        //
        for(n=2; n<=maxn; n++)
        {
            for(skind=0; skind<=2; skind++)
            {
                for(pkind=0; pkind<=2; pkind++)
                {
                    for(periodicity=0; periodicity<=1; periodicity++)
                    {
                        periodic = periodicity==1;
                        
                        //
                        // skip unsupported combinations of parameters
                        //
                        if( periodic && n<3 )
                        {
                            continue;
                        }
                        if( periodic && skind==0 )
                        {
                            continue;
                        }
                        if( n<5 && skind==0 )
                        {
                            continue;
                        }
                        
                        //
                        // init
                        //
                        xy.setlength(n, 2);
                        xyz.setlength(n, 3);
                        apserv::taskgenint1dequidist<Precision>(amp::ampf<Precision>(-1), amp::ampf<Precision>(+1), n, t2, x);
                        amp::vmove(xy.getcolumn(0, 0, n-1), x.getvector(0, n-1));
                        amp::vmove(xyz.getcolumn(0, 0, n-1), x.getvector(0, n-1));
                        apserv::taskgenint1dequidist<Precision>(amp::ampf<Precision>(-1), amp::ampf<Precision>(+1), n, t2, y);
                        amp::vmove(xy.getcolumn(1, 0, n-1), y.getvector(0, n-1));
                        amp::vmove(xyz.getcolumn(1, 0, n-1), y.getvector(0, n-1));
                        apserv::taskgenint1dequidist<Precision>(amp::ampf<Precision>(-1), amp::ampf<Precision>(+1), n, t2, z);
                        amp::vmove(xyz.getcolumn(2, 0, n-1), z.getvector(0, n-1));
                        unsetp2<Precision>(p2);
                        unsetp3<Precision>(p3);
                        if( periodic )
                        {
                            pspline::pspline2buildperiodic<Precision>(xy, n, skind, pkind, p2);
                            pspline::pspline3buildperiodic<Precision>(xyz, n, skind, pkind, p3);
                        }
                        else
                        {
                            pspline::pspline2build<Precision>(xy, n, skind, pkind, p2);
                            pspline::pspline3build<Precision>(xyz, n, skind, pkind, p3);
                        }
                        
                        //
                        // PSpline2ParameterValues() properties
                        //
                        pspline::pspline2parametervalues<Precision>(p2, tmpn, t2);
                        if( tmpn!=n )
                        {
                            p2errors = true;
                            continue;
                        }
                        pspline::pspline3parametervalues<Precision>(p3, tmpn, t3);
                        if( tmpn!=n )
                        {
                            p3errors = true;
                            continue;
                        }
                        p2errors = p2errors || t2(0)!=0;
                        p3errors = p3errors || t3(0)!=0;
                        for(i=1; i<=n-1; i++)
                        {
                            p2errors = p2errors || t2(i)<=t2(i-1);
                            p3errors = p3errors || t3(i)<=t3(i-1);
                        }
                        if( periodic )
                        {
                            p2errors = p2errors || t2(n-1)>=1;
                            p3errors = p3errors || t3(n-1)>=1;
                        }
                        else
                        {
                            p2errors = p2errors || t2(n-1)!=1;
                            p3errors = p3errors || t3(n-1)!=1;
                        }
                        
                        //
                        // Now we have parameter values stored at T,
                        // and want to test whether the actully correspond to
                        // points
                        //
                        for(i=0; i<=n-1; i++)
                        {
                            
                            //
                            // 2-dimensional test
                            //
                            pspline::pspline2calc<Precision>(p2, t2(i), vx, vy);
                            p2errors = p2errors || amp::abs<Precision>(vx-x(i))>threshold;
                            p2errors = p2errors || amp::abs<Precision>(vy-y(i))>threshold;
                            
                            //
                            // 3-dimensional test
                            //
                            pspline::pspline3calc<Precision>(p3, t3(i), vx, vy, vz);
                            p3errors = p3errors || amp::abs<Precision>(vx-x(i))>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vy-y(i))>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vz-z(i))>threshold;
                        }
                        
                        //
                        // Test periodicity (if needed)
                        //
                        if( periodic )
                        {
                            
                            //
                            // periodicity at nodes
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                
                                //
                                // 2-dimensional test
                                //
                                pspline::pspline2calc<Precision>(p2, t2(i)+ap::randominteger(10)-5, vx, vy);
                                p2errors = p2errors || amp::abs<Precision>(vx-x(i))>threshold;
                                p2errors = p2errors || amp::abs<Precision>(vy-y(i))>threshold;
                                pspline::pspline2diff<Precision>(p2, t2(i)+ap::randominteger(10)-5, vx, vdx, vy, vdy);
                                p2errors = p2errors || amp::abs<Precision>(vx-x(i))>threshold;
                                p2errors = p2errors || amp::abs<Precision>(vy-y(i))>threshold;
                                pspline::pspline2diff2<Precision>(p2, t2(i)+ap::randominteger(10)-5, vx, vdx, vd2x, vy, vdy, vd2y);
                                p2errors = p2errors || amp::abs<Precision>(vx-x(i))>threshold;
                                p2errors = p2errors || amp::abs<Precision>(vy-y(i))>threshold;
                                
                                //
                                // 3-dimensional test
                                //
                                pspline::pspline3calc<Precision>(p3, t3(i)+ap::randominteger(10)-5, vx, vy, vz);
                                p3errors = p3errors || amp::abs<Precision>(vx-x(i))>threshold;
                                p3errors = p3errors || amp::abs<Precision>(vy-y(i))>threshold;
                                p3errors = p3errors || amp::abs<Precision>(vz-z(i))>threshold;
                                pspline::pspline3diff<Precision>(p3, t3(i)+ap::randominteger(10)-5, vx, vdx, vy, vdy, vz, vdz);
                                p3errors = p3errors || amp::abs<Precision>(vx-x(i))>threshold;
                                p3errors = p3errors || amp::abs<Precision>(vy-y(i))>threshold;
                                p3errors = p3errors || amp::abs<Precision>(vz-z(i))>threshold;
                                pspline::pspline3diff2<Precision>(p3, t3(i)+ap::randominteger(10)-5, vx, vdx, vd2x, vy, vdy, vd2y, vz, vdz, vd2z);
                                p3errors = p3errors || amp::abs<Precision>(vx-x(i))>threshold;
                                p3errors = p3errors || amp::abs<Precision>(vy-y(i))>threshold;
                                p3errors = p3errors || amp::abs<Precision>(vz-z(i))>threshold;
                            }
                            
                            //
                            // periodicity between nodes
                            //
                            v0 = amp::ampf<Precision>::getRandom();
                            pspline::pspline2calc<Precision>(p2, v0, vx, vy);
                            pspline::pspline2calc<Precision>(p2, v0+ap::randominteger(10)-5, vx2, vy2);
                            p2errors = p2errors || amp::abs<Precision>(vx-vx2)>threshold;
                            p2errors = p2errors || amp::abs<Precision>(vy-vy2)>threshold;
                            pspline::pspline3calc<Precision>(p3, v0, vx, vy, vz);
                            pspline::pspline3calc<Precision>(p3, v0+ap::randominteger(10)-5, vx2, vy2, vz2);
                            p3errors = p3errors || amp::abs<Precision>(vx-vx2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vy-vy2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vz-vz2)>threshold;
                            
                            //
                            // near-boundary test for continuity of function values and derivatives:
                            // 2-dimensional curve
                            //
                            ap::ap_error::make_assertion(skind==1 || skind==2);
                            v0 = 100*amp::ampf<Precision>::getAlgoPascalEpsilon();
                            v1 = 1-v0;
                            pspline::pspline2calc<Precision>(p2, v0, vx, vy);
                            pspline::pspline2calc<Precision>(p2, v1, vx2, vy2);
                            p2errors = p2errors || amp::abs<Precision>(vx-vx2)>threshold;
                            p2errors = p2errors || amp::abs<Precision>(vy-vy2)>threshold;
                            pspline::pspline2diff<Precision>(p2, v0, vx, vdx, vy, vdy);
                            pspline::pspline2diff<Precision>(p2, v1, vx2, vdx2, vy2, vdy2);
                            p2errors = p2errors || amp::abs<Precision>(vx-vx2)>threshold;
                            p2errors = p2errors || amp::abs<Precision>(vy-vy2)>threshold;
                            p2errors = p2errors || amp::abs<Precision>(vdx-vdx2)>nonstrictthreshold;
                            p2errors = p2errors || amp::abs<Precision>(vdy-vdy2)>nonstrictthreshold;
                            pspline::pspline2diff2<Precision>(p2, v0, vx, vdx, vd2x, vy, vdy, vd2y);
                            pspline::pspline2diff2<Precision>(p2, v1, vx2, vdx2, vd2x2, vy2, vdy2, vd2y2);
                            p2errors = p2errors || amp::abs<Precision>(vx-vx2)>threshold;
                            p2errors = p2errors || amp::abs<Precision>(vy-vy2)>threshold;
                            p2errors = p2errors || amp::abs<Precision>(vdx-vdx2)>nonstrictthreshold;
                            p2errors = p2errors || amp::abs<Precision>(vdy-vdy2)>nonstrictthreshold;
                            if( skind==2 )
                            {
                                
                                //
                                // second derivative test only for cubic splines
                                //
                                p2errors = p2errors || amp::abs<Precision>(vd2x-vd2x2)>nonstrictthreshold;
                                p2errors = p2errors || amp::abs<Precision>(vd2y-vd2y2)>nonstrictthreshold;
                            }
                            
                            //
                            // near-boundary test for continuity of function values and derivatives:
                            // 3-dimensional curve
                            //
                            ap::ap_error::make_assertion(skind==1 || skind==2);
                            v0 = 100*amp::ampf<Precision>::getAlgoPascalEpsilon();
                            v1 = 1-v0;
                            pspline::pspline3calc<Precision>(p3, v0, vx, vy, vz);
                            pspline::pspline3calc<Precision>(p3, v1, vx2, vy2, vz2);
                            p3errors = p3errors || amp::abs<Precision>(vx-vx2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vy-vy2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vz-vz2)>threshold;
                            pspline::pspline3diff<Precision>(p3, v0, vx, vdx, vy, vdy, vz, vdz);
                            pspline::pspline3diff<Precision>(p3, v1, vx2, vdx2, vy2, vdy2, vz2, vdz2);
                            p3errors = p3errors || amp::abs<Precision>(vx-vx2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vy-vy2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vz-vz2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vdx-vdx2)>nonstrictthreshold;
                            p3errors = p3errors || amp::abs<Precision>(vdy-vdy2)>nonstrictthreshold;
                            p3errors = p3errors || amp::abs<Precision>(vdz-vdz2)>nonstrictthreshold;
                            pspline::pspline3diff2<Precision>(p3, v0, vx, vdx, vd2x, vy, vdy, vd2y, vz, vdz, vd2z);
                            pspline::pspline3diff2<Precision>(p3, v1, vx2, vdx2, vd2x2, vy2, vdy2, vd2y2, vz2, vdz2, vd2z2);
                            p3errors = p3errors || amp::abs<Precision>(vx-vx2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vy-vy2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vz-vz2)>threshold;
                            p3errors = p3errors || amp::abs<Precision>(vdx-vdx2)>nonstrictthreshold;
                            p3errors = p3errors || amp::abs<Precision>(vdy-vdy2)>nonstrictthreshold;
                            p3errors = p3errors || amp::abs<Precision>(vdz-vdz2)>nonstrictthreshold;
                            if( skind==2 )
                            {
                                
                                //
                                // second derivative test only for cubic splines
                                //
                                p3errors = p3errors || amp::abs<Precision>(vd2x-vd2x2)>nonstrictthreshold;
                                p3errors = p3errors || amp::abs<Precision>(vd2y-vd2y2)>nonstrictthreshold;
                                p3errors = p3errors || amp::abs<Precision>(vd2z-vd2z2)>nonstrictthreshold;
                            }
                        }
                    }
                }
            }
        }
        
        //
        // Test differentiation, tangents, calculation between nodes.
        //
        // Because differentiation is done in parameterization/spline/periodicity
        // oblivious manner, we don't have to test all possible combinations
        // of spline types and parameterizations.
        //
        // Actually we test special combination with properties which allow us
        // to easily solve this problem:
        // * 2 (3) variables
        // * first variable is sampled from equidistant grid on [0,1]
        // * other variables are random
        // * uniform parameterization is used
        // * periodicity - none
        // * spline type - any (we use cubic splines)
        // Same problem allows us to test calculation BETWEEN nodes.
        //
        for(n=2; n<=maxn; n++)
        {
            
            //
            // init
            //
            xy.setlength(n, 2);
            xyz.setlength(n, 3);
            apserv::taskgenint1dequidist<Precision>(amp::ampf<Precision>(0), amp::ampf<Precision>(+1), n, t, x);
            amp::vmove(xy.getcolumn(0, 0, n-1), x.getvector(0, n-1));
            amp::vmove(xyz.getcolumn(0, 0, n-1), x.getvector(0, n-1));
            apserv::taskgenint1dequidist<Precision>(amp::ampf<Precision>(0), amp::ampf<Precision>(+1), n, t, y);
            amp::vmove(xy.getcolumn(1, 0, n-1), y.getvector(0, n-1));
            amp::vmove(xyz.getcolumn(1, 0, n-1), y.getvector(0, n-1));
            apserv::taskgenint1dequidist<Precision>(amp::ampf<Precision>(0), amp::ampf<Precision>(+1), n, t, z);
            amp::vmove(xyz.getcolumn(2, 0, n-1), z.getvector(0, n-1));
            unsetp2<Precision>(p2);
            unsetp3<Precision>(p3);
            pspline::pspline2build<Precision>(xy, n, 2, 0, p2);
            pspline::pspline3build<Precision>(xyz, n, 2, 0, p3);
            
            //
            // Test 2D/3D spline:
            // * build non-parametric cubic spline from T and X/Y
            // * calculate its value and derivatives at V0
            // * compare with Spline2Calc/Spline2Diff/Spline2Diff2
            // Because of task properties both variants should
            // return same answer.
            //
            v0 = amp::ampf<Precision>::getRandom();
            spline1d::spline1dbuildcubic<Precision>(t, x, n, 0, amp::ampf<Precision>("0.0"), 0, amp::ampf<Precision>("0.0"), s);
            spline1d::spline1ddiff<Precision>(s, v0, vx2, vdx2, vd2x2);
            spline1d::spline1dbuildcubic<Precision>(t, y, n, 0, amp::ampf<Precision>("0.0"), 0, amp::ampf<Precision>("0.0"), s);
            spline1d::spline1ddiff<Precision>(s, v0, vy2, vdy2, vd2y2);
            spline1d::spline1dbuildcubic<Precision>(t, z, n, 0, amp::ampf<Precision>("0.0"), 0, amp::ampf<Precision>("0.0"), s);
            spline1d::spline1ddiff<Precision>(s, v0, vz2, vdz2, vd2z2);
            
            //
            // 2D test
            //
            pspline::pspline2calc<Precision>(p2, v0, vx, vy);
            p2errors = p2errors || amp::abs<Precision>(vx-vx2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vy-vy2)>threshold;
            pspline::pspline2diff<Precision>(p2, v0, vx, vdx, vy, vdy);
            p2errors = p2errors || amp::abs<Precision>(vx-vx2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vy-vy2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vdx-vdx2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vdy-vdy2)>threshold;
            pspline::pspline2diff2<Precision>(p2, v0, vx, vdx, vd2x, vy, vdy, vd2y);
            p2errors = p2errors || amp::abs<Precision>(vx-vx2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vy-vy2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vdx-vdx2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vdy-vdy2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vd2x-vd2x2)>threshold;
            p2errors = p2errors || amp::abs<Precision>(vd2y-vd2y2)>threshold;
            
            //
            // 3D test
            //
            pspline::pspline3calc<Precision>(p3, v0, vx, vy, vz);
            p3errors = p3errors || amp::abs<Precision>(vx-vx2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vy-vy2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vz-vz2)>threshold;
            pspline::pspline3diff<Precision>(p3, v0, vx, vdx, vy, vdy, vz, vdz);
            p3errors = p3errors || amp::abs<Precision>(vx-vx2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vy-vy2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vz-vz2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vdx-vdx2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vdy-vdy2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vdz-vdz2)>threshold;
            pspline::pspline3diff2<Precision>(p3, v0, vx, vdx, vd2x, vy, vdy, vd2y, vz, vdz, vd2z);
            p3errors = p3errors || amp::abs<Precision>(vx-vx2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vy-vy2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vz-vz2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vdx-vdx2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vdy-vdy2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vdz-vdz2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vd2x-vd2x2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vd2y-vd2y2)>threshold;
            p3errors = p3errors || amp::abs<Precision>(vd2z-vd2z2)>threshold;
            
            //
            // Test tangents for 2D/3D
            //
            pspline::pspline2tangent<Precision>(p2, v0, vx, vy);
            p2errors = p2errors || amp::abs<Precision>(vx-vdx2/apserv::safepythag2<Precision>(vdx2, vdy2))>threshold;
            p2errors = p2errors || amp::abs<Precision>(vy-vdy2/apserv::safepythag2<Precision>(vdx2, vdy2))>threshold;
            pspline::pspline3tangent<Precision>(p3, v0, vx, vy, vz);
            p3errors = p3errors || amp::abs<Precision>(vx-vdx2/apserv::safepythag3<Precision>(vdx2, vdy2, vdz2))>threshold;
            p3errors = p3errors || amp::abs<Precision>(vy-vdy2/apserv::safepythag3<Precision>(vdx2, vdy2, vdz2))>threshold;
            p3errors = p3errors || amp::abs<Precision>(vz-vdz2/apserv::safepythag3<Precision>(vdx2, vdy2, vdz2))>threshold;
        }
        
        //
        // Arc length test.
        //
        // Simple problem with easy solution (points on a straight line with
        // uniform parameterization).
        //
        for(n=2; n<=maxn; n++)
        {
            xy.setlength(n, 2);
            xyz.setlength(n, 3);
            for(i=0; i<=n-1; i++)
            {
                xy(i,0) = i;
                xy(i,1) = i;
                xyz(i,0) = i;
                xyz(i,1) = i;
                xyz(i,2) = i;
            }
            pspline::pspline2build<Precision>(xy, n, 1, 0, p2);
            pspline::pspline3build<Precision>(xyz, n, 1, 0, p3);
            a = amp::ampf<Precision>::getRandom();
            b = amp::ampf<Precision>::getRandom();
            p2errors = p2errors || amp::abs<Precision>(pspline::pspline2arclength<Precision>(p2, a, b)-(b-a)*amp::sqrt<Precision>(amp::ampf<Precision>(2))*(n-1))>nonstrictthreshold;
            p3errors = p3errors || amp::abs<Precision>(pspline::pspline3arclength<Precision>(p3, a, b)-(b-a)*amp::sqrt<Precision>(amp::ampf<Precision>(3))*(n-1))>nonstrictthreshold;
        }
        
        //
        // report
        //
        waserrors = p2errors || p3errors;
        if( !silent )
        {
            printf("TESTING SPLINE INTERPOLATION\n");
            
            //
            // Normal tests
            //
            printf("2D TEST:                                 ");
            if( p2errors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("3D TEST:                                 ");
            if( p3errors )
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
    Unset spline, i.e. initialize it with random garbage
    *************************************************************************/
    template<unsigned int Precision>
    void unsetp2(pspline::pspline2interpolant<Precision>& p)
    {
        ap::template_2d_array< amp::ampf<Precision> > xy;


        xy.setlength(2, 2);
        xy(0,0) = -1;
        xy(0,1) = -1;
        xy(1,0) = +1;
        xy(1,1) = +1;
        pspline::pspline2build<Precision>(xy, 2, 1, 0, p);
    }


    /*************************************************************************
    Unset spline, i.e. initialize it with random garbage
    *************************************************************************/
    template<unsigned int Precision>
    void unsetp3(pspline::pspline3interpolant<Precision>& p)
    {
        ap::template_2d_array< amp::ampf<Precision> > xy;


        xy.setlength(2, 3);
        xy(0,0) = -1;
        xy(0,1) = -1;
        xy(0,2) = -1;
        xy(1,0) = +1;
        xy(1,1) = +1;
        xy(1,2) = +1;
        pspline::pspline3build<Precision>(xy, 2, 1, 0, p);
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
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testpsplineunit_test_silent()
    {
        bool result;


        result = testpsplineinterpolation<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testpsplineunit_test()
    {
        bool result;


        result = testpsplineinterpolation<Precision>(false);
        return result;
    }
} // namespace

#endif
