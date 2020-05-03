#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "spline1d.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    ap::template_1d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > y;
    ap::template_1d_array< amp::ampf<Precision> > w;
    ap::template_1d_array< amp::ampf<Precision> > xc;
    ap::template_1d_array< amp::ampf<Precision> > yc;
    ap::template_1d_array< int > dc;
    int n;
    int i;
    int info;
    spline1d::spline1dinterpolant<Precision> s;
    amp::ampf<Precision> t;
    spline1d::spline1dfitreport<Precision> rep;


    
    //
    // Fitting by constrained Hermite spline
    //
    printf("FITTING BY CONSTRAINED HERMITE SPLINE\n\n");
    printf("F(x)=sin(x)      function being fitted\n");
    printf("[0, pi]          interval\n");
    printf("M=6              number of basis functions to use\n");
    printf("S(0)=0           first constraint\n");
    printf("S(pi)=0          second constraint\n");
    printf("N=100            number of points to fit\n");
    
    //
    // Create and fit:
    // * X  contains points
    // * Y  contains values
    // * W  contains weights
    // * XC contains constraints locations
    // * YC contains constraints values
    // * DC contains derivative indexes (0 = constrained function value)
    //
    n = 100;
    x.setlength(n);
    y.setlength(n);
    w.setlength(n);
    for(i=0; i<=n-1; i++)
    {
        x(i) = amp::pi<Precision>()*i/(n-1);
        y(i) = amp::sin<Precision>(x(i));
        w(i) = 1;
    }
    xc.setlength(2);
    yc.setlength(2);
    dc.setlength(2);
    xc(0) = 0;
    yc(0) = 0;
    dc(0) = 0;
    xc(0) = amp::pi<Precision>();
    yc(0) = 0;
    dc(0) = 0;
    spline1d::spline1dfithermitewc<Precision>(x, y, w, n, xc, yc, dc, 2, 6, info, s, rep);
    
    //
    // Output results
    //
    if( info>0 )
    {
        printf("\nOK, we have finished\n\n");
        printf("     x   F(x)   S(x)  Error\n");
        t = 0;
        while( t<amp::ampf<Precision>("0.999999")*amp::pi<Precision>() )
        {
            printf("%6.3lf %6.3lf %6.3lf %6.3lf\n",
                double(amp::ampf<Precision>(t).toDouble()),
                double(amp::ampf<Precision>(amp::sin<Precision>(t)).toDouble()),
                double(amp::ampf<Precision>(spline1d::spline1dcalc<Precision>(s, t)).toDouble()),
                double(amp::ampf<Precision>(amp::abs<Precision>(spline1d::spline1dcalc<Precision>(s, t)-amp::sin<Precision>(t))).toDouble()));
            t = amp::minimum<Precision>(amp::pi<Precision>(), t+amp::ampf<Precision>("0.25"));
        }
        printf("%6.3lf %6.3lf %6.3lf %6.3lf\n\n",
            double(amp::ampf<Precision>(t).toDouble()),
            double(amp::ampf<Precision>(amp::sin<Precision>(t)).toDouble()),
            double(amp::ampf<Precision>(spline1d::spline1dcalc<Precision>(s, t)).toDouble()),
            double(amp::ampf<Precision>(amp::abs<Precision>(spline1d::spline1dcalc<Precision>(s, t)-amp::sin<Precision>(t))).toDouble()));
        printf("rms error is %6.3lf\n",
            double(amp::ampf<Precision>(rep.rmserror).toDouble()));
        printf("max error is %6.3lf\n",
            double(amp::ampf<Precision>(rep.maxerror).toDouble()));
        printf("S(0) = S(pi) = 0 (exactly)\n\n");
    }
    else
    {
        printf("\nSomething wrong, Info=%0ld",
            long(info));
    }    
    return 0;
}