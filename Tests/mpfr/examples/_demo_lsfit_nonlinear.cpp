#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "lsfit.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    int m;
    int n;
    int k;
    ap::template_1d_array< amp::ampf<Precision> > y;
    ap::template_2d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > c;
    lsfit::lsfitreport<Precision> rep;
    lsfit::lsfitstate<Precision> state;
    int info;
    amp::ampf<Precision> epsf;
    amp::ampf<Precision> epsx;
    int maxits;
    int i;
    int j;
    amp::ampf<Precision> a;
    amp::ampf<Precision> b;


    printf("Fitting 0.5(1+cos(x)) on [-pi,+pi] with exp(-alpha*x^2)\n");
    
    //
    // Fitting 0.5(1+cos(x)) on [-pi,+pi] with Gaussian exp(-alpha*x^2):
    // * without Hessian (gradient only)
    // * using alpha=1 as initial value
    // * using 1000 uniformly distributed points to fit to
    //
    // Notes:
    // * N - number of points
    // * M - dimension of space where points reside
    // * K - number of parameters being fitted
    //
    n = 1000;
    m = 1;
    k = 1;
    a = -amp::pi<Precision>();
    b = +amp::pi<Precision>();
    
    //
    // Prepare task matrix
    //
    y.setlength(n);
    x.setlength(n, m);
    c.setlength(k);
    for(i=0; i<=n-1; i++)
    {
        x(i,0) = a+(b-a)*i/(n-1);
        y(i) = amp::ampf<Precision>("0.5")*(1+amp::cos<Precision>(x(i,0)));
    }
    c(0) = amp::ampf<Precision>("1.0");
    epsf = amp::ampf<Precision>("0.0");
    epsx = amp::ampf<Precision>("0.0001");
    maxits = 0;
    
    //
    // Solve
    //
    lsfit::lsfitnonlinearfg<Precision>(x, y, c, n, m, k, true, state);
    lsfit::lsfitnonlinearsetcond<Precision>(state, epsf, epsx, maxits);
    while( lsfit::lsfitnonlineariteration<Precision>(state) )
    {
        if( state.needf )
        {
            
            //
            // F(x) = Exp(-alpha*x^2)
            //
            state.f = amp::exp<Precision>(-state.c(0)*amp::sqr<Precision>(state.x(0)));
        }
        if( state.needfg )
        {
            
            //
            // F(x)      = Exp(-alpha*x^2)
            // dF/dAlpha = (-x^2)*Exp(-alpha*x^2)
            //
            state.f = amp::exp<Precision>(-state.c(0)*amp::sqr<Precision>(state.x(0)));
            state.g(0) = -amp::sqr<Precision>(state.x(0))*state.f;
        }
    }
    lsfit::lsfitnonlinearresults<Precision>(state, info, c, rep);
    printf("alpha:   %0.3lf\n",
        double(amp::ampf<Precision>(c(0)).toDouble()));
    printf("rms.err: %0.3lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()));
    printf("max.err: %0.3lf\n",
        double(amp::ampf<Precision>(rep.maxerror).toDouble()));
    printf("Termination type: %0ld\n",
        long(info));
    printf("\n\n");    
    return 0;
}