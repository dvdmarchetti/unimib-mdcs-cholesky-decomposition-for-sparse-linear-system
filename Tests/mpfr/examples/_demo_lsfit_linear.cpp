#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "lsfit.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    int m;
    int n;
    ap::template_1d_array< amp::ampf<Precision> > y;
    ap::template_2d_array< amp::ampf<Precision> > fmatrix;
    ap::template_2d_array< amp::ampf<Precision> > cmatrix;
    lsfit::lsfitreport<Precision> rep;
    int info;
    ap::template_1d_array< amp::ampf<Precision> > c;
    int i;
    int j;
    amp::ampf<Precision> x;
    amp::ampf<Precision> a;
    amp::ampf<Precision> b;


    printf("\n\nFitting tan(x) by third degree polynomial\n\n");
    printf("Fit type             rms.err max.err    p(0)   dp(0)\n");
    
    //
    // Fitting tan(x) at [0, 0.4*pi] by third degree polynomial:
    // a) without constraints
    // b) constrained at x=0: p(0)=0
    // c) constrained at x=0: p'(0)=1
    // c) constrained at x=0: p(0)=0, p'(0)=1
    //
    m = 4;
    n = 100;
    a = 0;
    b = amp::ampf<Precision>("0.4")*amp::pi<Precision>();
    
    //
    // Prepare task matrix
    //
    y.setlength(n);
    fmatrix.setlength(n, m);
    for(i=0; i<=n-1; i++)
    {
        x = a+(b-a)*i/(n-1);
        y(i) = amp::tan<Precision>(x);
        fmatrix(i,0) = amp::ampf<Precision>("1.0");
        for(j=1; j<=m-1; j++)
        {
            fmatrix(i,j) = x*fmatrix(i,j-1);
        }
    }
    
    //
    // Solve unconstrained task
    //
    lsfit::lsfitlinear<Precision>(y, fmatrix, n, m, info, c, rep);
    printf("Unconstrained        %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()),
        double(amp::ampf<Precision>(rep.maxerror).toDouble()),
        double(amp::ampf<Precision>(c(0)).toDouble()),
        double(amp::ampf<Precision>(c(1)).toDouble()));
    
    //
    // Solve constrained task, p(0)=0
    // Prepare constraints matrix:
    // * first M columns store values of basis functions at X=0
    // * last column stores zero (desired value at X=0)
    //
    cmatrix.setlength(1, m+1);
    cmatrix(0,0) = 1;
    for(i=1; i<=m-1; i++)
    {
        cmatrix(0,i) = 0;
    }
    cmatrix(0,m) = 0;
    lsfit::lsfitlinearc<Precision>(y, fmatrix, cmatrix, n, m, 1, info, c, rep);
    printf("Constrained, p(0)=0  %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()),
        double(amp::ampf<Precision>(rep.maxerror).toDouble()),
        double(amp::ampf<Precision>(c(0)).toDouble()),
        double(amp::ampf<Precision>(c(1)).toDouble()));
    
    //
    // Solve constrained task, p'(0)=0
    // Prepare constraints matrix:
    // * first M columns store derivatives of basis functions at X=0
    // * last column stores 1.0 (desired derivative at X=0)
    //
    cmatrix.setlength(1, m+1);
    for(i=0; i<=m-1; i++)
    {
        cmatrix(0,i) = 0;
    }
    cmatrix(0,1) = 1;
    cmatrix(0,m) = 1;
    lsfit::lsfitlinearc<Precision>(y, fmatrix, cmatrix, n, m, 1, info, c, rep);
    printf("Constrained, dp(0)=1 %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()),
        double(amp::ampf<Precision>(rep.maxerror).toDouble()),
        double(amp::ampf<Precision>(c(0)).toDouble()),
        double(amp::ampf<Precision>(c(1)).toDouble()));
    
    //
    // Solve constrained task, p(0)=0, p'(0)=0
    // Prepare constraints matrix:
    // * first M columns store values/derivatives of basis functions at X=0
    // * last column stores desired values/derivative at X=0
    //
    cmatrix.setlength(2, m+1);
    cmatrix(0,0) = 1;
    for(i=1; i<=m-1; i++)
    {
        cmatrix(0,i) = 0;
    }
    cmatrix(0,m) = 0;
    for(i=0; i<=m-1; i++)
    {
        cmatrix(1,i) = 0;
    }
    cmatrix(1,1) = 1;
    cmatrix(1,m) = 1;
    lsfit::lsfitlinearc<Precision>(y, fmatrix, cmatrix, n, m, 2, info, c, rep);
    printf("Constrained, both    %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()),
        double(amp::ampf<Precision>(rep.maxerror).toDouble()),
        double(amp::ampf<Precision>(c(0)).toDouble()),
        double(amp::ampf<Precision>(c(1)).toDouble()));
    printf("\n\n");    
    return 0;
}