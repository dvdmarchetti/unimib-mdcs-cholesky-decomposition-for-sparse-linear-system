#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "spline1d.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    ap::template_1d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > y;
    ap::template_1d_array< amp::ampf<Precision> > d;
    int n;
    int i;
    amp::ampf<Precision> t;
    spline1d::spline1dinterpolant<Precision> s;
    amp::ampf<Precision> err;
    amp::ampf<Precision> maxerr;


    
    //
    // Interpolation by natural Cubic spline.
    //
    printf("INTERPOLATION BY HERMITE SPLINE\n\n");
    printf("F(x)=sin(x), [0, pi], 3 nodes\n\n");
    printf("     x   F(x)   S(x)  Error\n");
    
    //
    // Create spline
    //
    n = 3;
    x.setlength(n);
    y.setlength(n);
    d.setlength(n);
    for(i=0; i<=n-1; i++)
    {
        x(i) = amp::pi<Precision>()*i/(n-1);
        y(i) = amp::sin<Precision>(x(i));
        d(i) = amp::cos<Precision>(x(i));
    }
    spline1d::spline1dbuildhermite<Precision>(x, y, d, n, s);
    
    //
    // Output results
    //
    t = 0;
    maxerr = 0;
    while( t<amp::ampf<Precision>("0.999999")*amp::pi<Precision>() )
    {
        err = amp::abs<Precision>(spline1d::spline1dcalc<Precision>(s, t)-amp::sin<Precision>(t));
        maxerr = amp::maximum<Precision>(err, maxerr);
        printf("%6.3lf %6.3lf %6.3lf %6.3lf\n",
            double(amp::ampf<Precision>(t).toDouble()),
            double(amp::ampf<Precision>(amp::sin<Precision>(t)).toDouble()),
            double(amp::ampf<Precision>(spline1d::spline1dcalc<Precision>(s, t)).toDouble()),
            double(amp::ampf<Precision>(err).toDouble()));
        t = amp::minimum<Precision>(amp::pi<Precision>(), t+amp::ampf<Precision>("0.25"));
    }
    err = amp::abs<Precision>(spline1d::spline1dcalc<Precision>(s, amp::pi<Precision>())-amp::sin<Precision>(amp::pi<Precision>()));
    maxerr = amp::maximum<Precision>(err, maxerr);
    printf("%6.3lf %6.3lf %6.3lf %6.3lf\n\n",
        double(amp::ampf<Precision>(amp::pi<Precision>()).toDouble()),
        double(amp::ampf<Precision>(amp::sin<Precision>(amp::pi<Precision>())).toDouble()),
        double(amp::ampf<Precision>(spline1d::spline1dcalc<Precision>(s, amp::pi<Precision>())).toDouble()),
        double(amp::ampf<Precision>(err).toDouble()));
    printf("max|error| = %0.3lf\n",
        double(amp::ampf<Precision>(maxerr).toDouble()));
    printf("Try other demos (spline1d_linear, spline1d_cubic) and compare errors...\n\n\n");    
    return 0;
}