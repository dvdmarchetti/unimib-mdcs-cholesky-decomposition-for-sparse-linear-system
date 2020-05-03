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


    printf("Fitting 1-x^2 on [-1,+1] with cos(alpha*pi*x)+beta\n");
    
    //
    // Fitting 1-x^2 on [-1,+1] with cos(alpha*pi*x)+beta:
    // * using Hessian
    // * using alpha=1 and beta=0 as initial values
    // * using 1000 uniformly distributed points to fit to
    //
    // Notes:
    // * N - number of points
    // * M - dimension of space where points reside
    // * K - number of parameters being fitted
    //
    n = 1000;
    m = 1;
    k = 2;
    a = -1;
    b = +1;
    
    //
    // Prepare task matrix
    //
    y.setlength(n);
    x.setlength(n, m);
    c.setlength(k);
    for(i=0; i<=n-1; i++)
    {
        x(i,0) = a+(b-a)*i/(n-1);
        y(i) = 1-amp::sqr<Precision>(x(i,0));
    }
    c(0) = amp::ampf<Precision>("1.0");
    c(1) = amp::ampf<Precision>("0.0");
    epsf = amp::ampf<Precision>("0.0");
    epsx = amp::ampf<Precision>("0.0001");
    maxits = 0;
    
    //
    // Solve
    //
    lsfit::lsfitnonlinearfgh<Precision>(x, y, c, n, m, k, state);
    lsfit::lsfitnonlinearsetcond<Precision>(state, epsf, epsx, maxits);
    while( lsfit::lsfitnonlineariteration<Precision>(state) )
    {
        
        //
        // F(x) = Cos(alpha*pi*x)+beta
        //
        state.f = amp::cos<Precision>(state.c(0)*amp::pi<Precision>()*state.x(0))+state.c(1);
        
        //
        // F(x)      = Cos(alpha*pi*x)+beta
        // dF/dAlpha = -pi*x*Sin(alpha*pi*x)
        // dF/dBeta  = 1.0
        //
        if( state.needfg || state.needfgh )
        {
            state.g(0) = -amp::pi<Precision>()*state.x(0)*amp::sin<Precision>(state.c(0)*amp::pi<Precision>()*state.x(0));
            state.g(1) = amp::ampf<Precision>("1.0");
        }
        
        //
        // F(x)            = Cos(alpha*pi*x)+beta
        // d2F/dAlpha2     = -(pi*x)^2*Cos(alpha*pi*x)
        // d2F/dAlphadBeta = 0
        // d2F/dBeta2     =  0
        //
        if( state.needfgh )
        {
            state.h(0,0) = -amp::sqr<Precision>(amp::pi<Precision>()*state.x(0))*amp::cos<Precision>(state.c(0)*amp::pi<Precision>()*state.x(0));
            state.h(0,1) = amp::ampf<Precision>("0.0");
            state.h(1,0) = amp::ampf<Precision>("0.0");
            state.h(1,1) = amp::ampf<Precision>("0.0");
        }
    }
    lsfit::lsfitnonlinearresults<Precision>(state, info, c, rep);
    printf("alpha:   %0.3lf\n",
        double(amp::ampf<Precision>(c(0)).toDouble()));
    printf("beta:    %0.3lf\n",
        double(amp::ampf<Precision>(c(1)).toDouble()));
    printf("rms.err: %0.3lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()));
    printf("max.err: %0.3lf\n",
        double(amp::ampf<Precision>(rep.maxerror).toDouble()));
    printf("Termination type: %0ld\n",
        long(info));
    printf("\n\n");    
    return 0;
}