
#ifndef _testnearestneighborunit_h
#define _testnearestneighborunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "tsort.h"
#include "nearestneighbor.h"
namespace testnearestneighborunit
{
    template<unsigned int Precision>
    bool testnearestneighbor(bool silent);
    template<unsigned int Precision>
    void unset2d(ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    void unset1d(ap::template_1d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    bool kdtresultsdifferent(const ap::template_2d_array< amp::ampf<Precision> >& refxy,
        int ntotal,
        const ap::template_2d_array< amp::ampf<Precision> >& qx,
        const ap::template_2d_array< amp::ampf<Precision> >& qxy,
        const ap::template_1d_array< int >& qt,
        int n,
        int nx,
        int ny);
    template<unsigned int Precision>
    amp::ampf<Precision> vnorm(const ap::template_1d_array< amp::ampf<Precision> >& x,
        int n,
        int normtype);
    template<unsigned int Precision>
    void testkdtuniform(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        const int& n,
        const int& nx,
        const int& ny,
        const int& normtype,
        bool& kdterrors);
    template<unsigned int Precision>
    bool testnearestneighborunit_test_silent();
    template<unsigned int Precision>
    bool testnearestneighborunit_test();


    /*************************************************************************
    Testing Nearest Neighbor Search
    *************************************************************************/
    template<unsigned int Precision>
    bool testnearestneighbor(bool silent)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > xy;
        int i;
        int j;
        amp::ampf<Precision> v;
        int normtype;
        int nx;
        int ny;
        int n;
        int smalln;
        int largen;
        int passcount;
        int pass;
        bool waserrors;
        bool kdterrors;


        kdterrors = false;
        passcount = 2;
        smalln = 256;
        largen = 2048;
        ny = 3;
        
        //
        //
        //
        for(pass=1; pass<=passcount; pass++)
        {
            for(normtype=0; normtype<=2; normtype++)
            {
                for(nx=1; nx<=3; nx++)
                {
                    
                    //
                    // Test in hypercube
                    //
                    xy.setlength(largen, nx+ny);
                    for(i=0; i<=largen-1; i++)
                    {
                        for(j=0; j<=nx+ny-1; j++)
                        {
                            xy(i,j) = 10*amp::ampf<Precision>::getRandom()-5;
                        }
                    }
                    for(n=1; n<=10; n++)
                    {
                        testkdtuniform<Precision>(xy, n, nx, ap::randominteger(ny+1), normtype, kdterrors);
                    }
                    testkdtuniform<Precision>(xy, largen, nx, ap::randominteger(ny+1), normtype, kdterrors);
                    
                    //
                    // Test clustered (2*N points, pairs of equal points)
                    //
                    xy.setlength(2*smalln, nx+ny);
                    for(i=0; i<=smalln-1; i++)
                    {
                        for(j=0; j<=nx+ny-1; j++)
                        {
                            xy(2*i+0,j) = 10*amp::ampf<Precision>::getRandom()-5;
                            xy(2*i+1,j) = xy(2*i+0,j);
                        }
                    }
                    testkdtuniform<Precision>(xy, 2*smalln, nx, ap::randominteger(ny+1), normtype, kdterrors);
                    
                    //
                    // Test degenerate case: all points are same except for one
                    //
                    xy.setlength(smalln, nx+ny);
                    v = amp::ampf<Precision>::getRandom();
                    for(i=0; i<=smalln-2; i++)
                    {
                        for(j=0; j<=nx+ny-1; j++)
                        {
                            xy(i,j) = v;
                        }
                    }
                    for(j=0; j<=nx+ny-1; j++)
                    {
                        xy(smalln-1,j) = 10*amp::ampf<Precision>::getRandom()-5;
                    }
                    testkdtuniform<Precision>(xy, smalln, nx, ap::randominteger(ny+1), normtype, kdterrors);
                }
            }
        }
        
        //
        // report
        //
        waserrors = kdterrors;
        if( !silent )
        {
            printf("TESTING NEAREST NEIGHBOR SEARCH\n");
            printf("* KD TREES:                              ");
            if( !kdterrors )
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
    Compare results from different queries:
    * X     just X-values
    * XY    X-values and Y-values
    * XT    X-values and tag values
    *************************************************************************/
    template<unsigned int Precision>
    bool kdtresultsdifferent(const ap::template_2d_array< amp::ampf<Precision> >& refxy,
        int ntotal,
        const ap::template_2d_array< amp::ampf<Precision> >& qx,
        const ap::template_2d_array< amp::ampf<Precision> >& qxy,
        const ap::template_1d_array< int >& qt,
        int n,
        int nx,
        int ny)
    {
        bool result;
        int i;
        int j;


        result = false;
        for(i=0; i<=n-1; i++)
        {
            if( qt(i)<0 || qt(i)>=ntotal )
            {
                result = true;
                return result;
            }
            for(j=0; j<=nx-1; j++)
            {
                result = result || qx(i,j)!=refxy(qt(i),j);
                result = result || qxy(i,j)!=refxy(qt(i),j);
            }
            for(j=0; j<=ny-1; j++)
            {
                result = result || qxy(i,nx+j)!=refxy(qt(i),nx+j);
            }
        }
        return result;
    }


    /*************************************************************************
    Returns norm
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> vnorm(const ap::template_1d_array< amp::ampf<Precision> >& x,
        int n,
        int normtype)
    {
        amp::ampf<Precision> result;
        int i;


        result = amp::ampf<Precision>::getRandom();
        if( normtype==0 )
        {
            result = 0;
            for(i=0; i<=n-1; i++)
            {
                result = amp::maximum<Precision>(result, amp::abs<Precision>(x(i)));
            }
            return result;
        }
        if( normtype==1 )
        {
            result = 0;
            for(i=0; i<=n-1; i++)
            {
                result = result+amp::abs<Precision>(x(i));
            }
            return result;
        }
        if( normtype==2 )
        {
            result = 0;
            for(i=0; i<=n-1; i++)
            {
                result = result+amp::sqr<Precision>(x(i));
            }
            result = amp::sqrt<Precision>(result);
            return result;
        }
        return result;
    }


    /*************************************************************************
    Testing Nearest Neighbor Search on uniformly distributed hypercube

    NormType: 0, 1, 2
    D: space dimension
    N: points count
    *************************************************************************/
    template<unsigned int Precision>
    void testkdtuniform(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        const int& n,
        const int& nx,
        const int& ny,
        const int& normtype,
        bool& kdterrors)
    {
        amp::ampf<Precision> errtol;
        ap::template_1d_array< int > tags;
        ap::template_1d_array< amp::ampf<Precision> > ptx;
        ap::template_1d_array< amp::ampf<Precision> > tmpx;
        ap::template_1d_array< bool > tmpb;
        nearestneighbor::kdtree<Precision> treex;
        nearestneighbor::kdtree<Precision> treexy;
        nearestneighbor::kdtree<Precision> treext;
        ap::template_2d_array< amp::ampf<Precision> > qx;
        ap::template_2d_array< amp::ampf<Precision> > qxy;
        ap::template_1d_array< int > qtags;
        ap::template_1d_array< amp::ampf<Precision> > qr;
        int kx;
        int kxy;
        int kt;
        int kr;
        amp::ampf<Precision> eps;
        int i;
        int j;
        int k;
        int task;
        bool isequal;
        amp::ampf<Precision> r;
        int q;
        int qcount;


        qcount = 10;
        
        //
        // Tol - roundoff error tolerance (for '>=' comparisons)
        //
        errtol = 100000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // fill tags
        //
        tags.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            tags(i) = i;
        }
        
        //
        // build trees
        //
        nearestneighbor::kdtreebuild<Precision>(xy, n, nx, 0, normtype, treex);
        nearestneighbor::kdtreebuild<Precision>(xy, n, nx, ny, normtype, treexy);
        nearestneighbor::kdtreebuildtagged<Precision>(xy, tags, n, nx, 0, normtype, treext);
        
        //
        // allocate arrays
        //
        tmpx.setlength(nx);
        tmpb.setlength(n);
        qx.setlength(n, nx);
        qxy.setlength(n, nx+ny);
        qtags.setlength(n);
        qr.setlength(n);
        ptx.setlength(nx);
        
        //
        // test general K-NN queries (with self-matches):
        // * compare results from different trees (must be equal) and
        //   check that correct (value,tag) pairs are returned
        // * test results from XT tree - let R be radius of query result.
        //   then all points not in result must be not closer than R.
        //
        for(q=1; q<=qcount; q++)
        {
            
            //
            // Select K: 1..N
            //
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                k = 1+ap::randominteger(n);
            }
            else
            {
                k = 1;
            }
            
            //
            // Select point (either one of the points, or random)
            //
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                i = ap::randominteger(n);
                amp::vmove(ptx.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
            }
            else
            {
                for(i=0; i<=nx-1; i++)
                {
                    ptx(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            
            //
            // Test:
            // * consistency of results from different queries
            // * points in query are IN the R-sphere (or at the boundary),
            //   and points not in query are outside of the R-sphere (or at the boundary)
            // * distances are correct and are ordered
            //
            kx = nearestneighbor::kdtreequeryknn<Precision>(treex, ptx, k, true);
            kxy = nearestneighbor::kdtreequeryknn<Precision>(treexy, ptx, k, true);
            kt = nearestneighbor::kdtreequeryknn<Precision>(treext, ptx, k, true);
            if( kx!=k || kxy!=k || kt!=k )
            {
                kdterrors = true;
                return;
            }
            kx = 0;
            kxy = 0;
            kt = 0;
            nearestneighbor::kdtreequeryresultsx<Precision>(treex, qx, kx);
            nearestneighbor::kdtreequeryresultsxy<Precision>(treexy, qxy, kxy);
            nearestneighbor::kdtreequeryresultstags<Precision>(treext, qtags, kt);
            nearestneighbor::kdtreequeryresultsdistances<Precision>(treext, qr, kr);
            if( kx!=k || kxy!=k || kt!=k || kr!=k )
            {
                kdterrors = true;
                return;
            }
            kdterrors = kdterrors || kdtresultsdifferent<Precision>(xy, n, qx, qxy, qtags, k, nx, ny);
            for(i=0; i<=n-1; i++)
            {
                tmpb(i) = true;
            }
            r = 0;
            for(i=0; i<=k-1; i++)
            {
                tmpb(qtags(i)) = false;
                amp::vmove(tmpx.getvector(0, nx-1), ptx.getvector(0, nx-1));
                amp::vsub(tmpx.getvector(0, nx-1), qx.getrow(i, 0, nx-1));
                r = amp::maximum<Precision>(r, vnorm<Precision>(tmpx, nx, normtype));
            }
            for(i=0; i<=n-1; i++)
            {
                if( tmpb(i) )
                {
                    amp::vmove(tmpx.getvector(0, nx-1), ptx.getvector(0, nx-1));
                    amp::vsub(tmpx.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
                    kdterrors = kdterrors || vnorm<Precision>(tmpx, nx, normtype)<r*(1-errtol);
                }
            }
            for(i=0; i<=k-2; i++)
            {
                kdterrors = kdterrors || qr(i)>qr(i+1);
            }
            for(i=0; i<=k-1; i++)
            {
                amp::vmove(tmpx.getvector(0, nx-1), ptx.getvector(0, nx-1));
                amp::vsub(tmpx.getvector(0, nx-1), xy.getrow(qtags(i), 0, nx-1));
                kdterrors = kdterrors || amp::abs<Precision>(vnorm<Precision>(tmpx, nx, normtype)-qr(i))>errtol;
            }
        }
        
        //
        // test general approximate K-NN queries (with self-matches):
        // * compare results from different trees (must be equal) and
        //   check that correct (value,tag) pairs are returned
        // * test results from XT tree - let R be radius of query result.
        //   then all points not in result must be not closer than R/(1+Eps).
        //
        for(q=1; q<=qcount; q++)
        {
            
            //
            // Select K: 1..N
            //
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                k = 1+ap::randominteger(n);
            }
            else
            {
                k = 1;
            }
            
            //
            // Select Eps
            //
            eps = amp::ampf<Precision>("0.5")+amp::ampf<Precision>::getRandom();
            
            //
            // Select point (either one of the points, or random)
            //
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                i = ap::randominteger(n);
                amp::vmove(ptx.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
            }
            else
            {
                for(i=0; i<=nx-1; i++)
                {
                    ptx(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            
            //
            // Test:
            // * consistency of results from different queries
            // * points in query are IN the R-sphere (or at the boundary),
            //   and points not in query are outside of the R-sphere (or at the boundary)
            // * distances are correct and are ordered
            //
            kx = nearestneighbor::kdtreequeryaknn<Precision>(treex, ptx, k, true, eps);
            kxy = nearestneighbor::kdtreequeryaknn<Precision>(treexy, ptx, k, true, eps);
            kt = nearestneighbor::kdtreequeryaknn<Precision>(treext, ptx, k, true, eps);
            if( kx!=k || kxy!=k || kt!=k )
            {
                kdterrors = true;
                return;
            }
            kx = 0;
            kxy = 0;
            kt = 0;
            nearestneighbor::kdtreequeryresultsx<Precision>(treex, qx, kx);
            nearestneighbor::kdtreequeryresultsxy<Precision>(treexy, qxy, kxy);
            nearestneighbor::kdtreequeryresultstags<Precision>(treext, qtags, kt);
            nearestneighbor::kdtreequeryresultsdistances<Precision>(treext, qr, kr);
            if( kx!=k || kxy!=k || kt!=k || kr!=k )
            {
                kdterrors = true;
                return;
            }
            kdterrors = kdterrors || kdtresultsdifferent<Precision>(xy, n, qx, qxy, qtags, k, nx, ny);
            for(i=0; i<=n-1; i++)
            {
                tmpb(i) = true;
            }
            r = 0;
            for(i=0; i<=k-1; i++)
            {
                tmpb(qtags(i)) = false;
                amp::vmove(tmpx.getvector(0, nx-1), ptx.getvector(0, nx-1));
                amp::vsub(tmpx.getvector(0, nx-1), qx.getrow(i, 0, nx-1));
                r = amp::maximum<Precision>(r, vnorm<Precision>(tmpx, nx, normtype));
            }
            for(i=0; i<=n-1; i++)
            {
                if( tmpb(i) )
                {
                    amp::vmove(tmpx.getvector(0, nx-1), ptx.getvector(0, nx-1));
                    amp::vsub(tmpx.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
                    kdterrors = kdterrors || vnorm<Precision>(tmpx, nx, normtype)<r*(1-errtol)/(1+eps);
                }
            }
            for(i=0; i<=k-2; i++)
            {
                kdterrors = kdterrors || qr(i)>qr(i+1);
            }
            for(i=0; i<=k-1; i++)
            {
                amp::vmove(tmpx.getvector(0, nx-1), ptx.getvector(0, nx-1));
                amp::vsub(tmpx.getvector(0, nx-1), xy.getrow(qtags(i), 0, nx-1));
                kdterrors = kdterrors || amp::abs<Precision>(vnorm<Precision>(tmpx, nx, normtype)-qr(i))>errtol;
            }
        }
        
        //
        // test general R-NN queries  (with self-matches):
        // * compare results from different trees (must be equal) and
        //   check that correct (value,tag) pairs are returned
        // * test results from XT tree - let R be radius of query result.
        //   then all points not in result must be not closer than R.
        //
        for(q=1; q<=qcount; q++)
        {
            
            //
            // Select R
            //
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.3") )
            {
                r = amp::maximum<Precision>(amp::ampf<Precision>::getRandom(), amp::ampf<Precision>::getAlgoPascalEpsilon());
            }
            else
            {
                r = amp::ampf<Precision>::getAlgoPascalEpsilon();
            }
            
            //
            // Select point (either one of the points, or random)
            //
            if( amp::ampf<Precision>::getRandom()>amp::ampf<Precision>("0.5") )
            {
                i = ap::randominteger(n);
                amp::vmove(ptx.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
            }
            else
            {
                for(i=0; i<=nx-1; i++)
                {
                    ptx(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
            }
            
            //
            // Test:
            // * consistency of results from different queries
            // * points in query are IN the R-sphere (or at the boundary),
            //   and points not in query are outside of the R-sphere (or at the boundary)
            // * distances are correct and are ordered
            //
            kx = nearestneighbor::kdtreequeryrnn<Precision>(treex, ptx, r, true);
            kxy = nearestneighbor::kdtreequeryrnn<Precision>(treexy, ptx, r, true);
            kt = nearestneighbor::kdtreequeryrnn<Precision>(treext, ptx, r, true);
            if( kxy!=kx || kt!=kx )
            {
                kdterrors = true;
                return;
            }
            kx = 0;
            kxy = 0;
            kt = 0;
            nearestneighbor::kdtreequeryresultsx<Precision>(treex, qx, kx);
            nearestneighbor::kdtreequeryresultsxy<Precision>(treexy, qxy, kxy);
            nearestneighbor::kdtreequeryresultstags<Precision>(treext, qtags, kt);
            nearestneighbor::kdtreequeryresultsdistances<Precision>(treext, qr, kr);
            if( kxy!=kx || kt!=kx || kr!=kx )
            {
                kdterrors = true;
                return;
            }
            kdterrors = kdterrors || kdtresultsdifferent<Precision>(xy, n, qx, qxy, qtags, kx, nx, ny);
            for(i=0; i<=n-1; i++)
            {
                tmpb(i) = true;
            }
            for(i=0; i<=kx-1; i++)
            {
                tmpb(qtags(i)) = false;
            }
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(tmpx.getvector(0, nx-1), ptx.getvector(0, nx-1));
                amp::vsub(tmpx.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
                if( tmpb(i) )
                {
                    kdterrors = kdterrors || vnorm<Precision>(tmpx, nx, normtype)<r*(1-errtol);
                }
                else
                {
                    kdterrors = kdterrors || vnorm<Precision>(tmpx, nx, normtype)>r*(1+errtol);
                }
            }
            for(i=0; i<=kx-2; i++)
            {
                kdterrors = kdterrors || qr(i)>qr(i+1);
            }
        }
        
        //
        // Test self-matching:
        // * self-match - nearest neighbor of each point in XY is the point itself
        // * no self-match - nearest neighbor is NOT the point itself
        //
        if( n>1 )
        {
            
            //
            // test for N=1 have non-general form, but it is not really needed
            //
            for(task=0; task<=1; task++)
            {
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(ptx.getvector(0, nx-1), xy.getrow(i, 0, nx-1));
                    kx = nearestneighbor::kdtreequeryknn<Precision>(treex, ptx, 1, task==0);
                    nearestneighbor::kdtreequeryresultsx<Precision>(treex, qx, kx);
                    if( kx!=1 )
                    {
                        kdterrors = true;
                        return;
                    }
                    isequal = true;
                    for(j=0; j<=nx-1; j++)
                    {
                        isequal = isequal && qx(0,j)==ptx(j);
                    }
                    if( task==0 )
                    {
                        kdterrors = kdterrors || !isequal;
                    }
                    else
                    {
                        kdterrors = kdterrors || isequal;
                    }
                }
            }
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testnearestneighborunit_test_silent()
    {
        bool result;


        result = testnearestneighbor<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testnearestneighborunit_test()
    {
        bool result;


        result = testnearestneighbor<Precision>(false);
        return result;
    }
} // namespace

#endif
