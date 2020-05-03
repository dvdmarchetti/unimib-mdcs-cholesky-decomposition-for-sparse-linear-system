#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "spline1d.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    ap::template_1d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > y;
    int n;
    int i;
    int info;
    spline1d::spline1dinterpolant<Precision> s;
    amp::ampf<Precision> t;
    spline1d::spline1dfitreport<Precision> rep;


    
    //
    // Fitting by unconstrained natural cubic spline
    //
    printf("FITTING BY UNCONSTRAINED NATURAL CUBIC SPLINE\n\n");
    printf("F(x)=sin(x)      function being fitted\n");
    printf("[0, pi]          interval\n");
    printf("M=4              number of basis functions to use\n");
    printf("N=100            number of points to fit\n");
    
    //
    // Create and fit
    //
    n = 100;
    x.setlength(n);
    y.setlength(n);
    for(i=0; i<=n-1; i++)
    {
        x(i) = amp::pi<Precision>()*i/(n-1);
        y(i) = amp::sin<Precision>(x(i));
    }
    spline1d::spline1dfitcubic<Precision>(x, y, n, 4, info, s, rep);
    
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
    }
    else
    {
        printf("\nSomething wrong, Info=%0ld",
            long(info));
    }    
    return 0;
}