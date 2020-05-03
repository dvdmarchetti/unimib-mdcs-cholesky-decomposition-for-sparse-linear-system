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
    amp::ampf<Precision> t;
    spline1d::spline1dinterpolant<Precision> s;
    amp::ampf<Precision> v;
    amp::ampf<Precision> dv;
    amp::ampf<Precision> d2v;
    amp::ampf<Precision> err;
    amp::ampf<Precision> maxerr;


    
    //
    // Demonstration of Spline1DCalc(), Spline1DDiff(), Spline1DIntegrate()
    //
    printf("DEMONSTRATION OF Spline1DCalc(), Spline1DDiff(), Spline1DIntegrate()\n\n");
    printf("F(x)=sin(x), [0, pi]\n");
    printf("Natural cubic spline with 3 nodes is used\n\n");
    
    //
    // Create spline
    //
    n = 3;
    x.setlength(n);
    y.setlength(n);
    for(i=0; i<=n-1; i++)
    {
        x(i) = amp::pi<Precision>()*i/(n-1);
        y(i) = amp::sin<Precision>(x(i));
    }
    spline1d::spline1dbuildcubic<Precision>(x, y, n, 2, amp::ampf<Precision>("0.0"), 2, amp::ampf<Precision>("0.0"), s);
    
    //
    // Output results
    //
    spline1d::spline1ddiff<Precision>(s, amp::ampf<Precision>(0), v, dv, d2v);
    printf("                 S(x)    F(x) \n");
    printf("function       %6.3lf  %6.3lf \n",
        double(amp::ampf<Precision>(spline1d::spline1dcalc<Precision>(s, amp::ampf<Precision>(0))).toDouble()),
        double(amp::ampf<Precision>(0).toDouble()));
    printf("d/dx(0)        %6.3lf  %6.3lf \n",
        double(amp::ampf<Precision>(dv).toDouble()),
        double(amp::ampf<Precision>(1).toDouble()));
    printf("d2/dx2(0)      %6.3lf  %6.3lf \n",
        double(amp::ampf<Precision>(d2v).toDouble()),
        double(amp::ampf<Precision>(0).toDouble()));
    printf("integral(0,pi) %6.3lf  %6.3lf \n",
        double(amp::ampf<Precision>(spline1d::spline1dintegrate<Precision>(s, amp::pi<Precision>())).toDouble()),
        double(amp::ampf<Precision>(2).toDouble()));
    printf("\n\n");    
    return 0;
}