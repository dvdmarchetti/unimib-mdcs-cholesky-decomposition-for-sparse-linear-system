#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "polint.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    ap::template_1d_array< amp::ampf<Precision> > y;
    int n;
    int i;
    amp::ampf<Precision> t;
    ratint::barycentricinterpolant<Precision> p;
    amp::ampf<Precision> v;
    amp::ampf<Precision> dv;
    amp::ampf<Precision> d2v;
    amp::ampf<Precision> err;
    amp::ampf<Precision> maxerr;


    
    //
    // Demonstration
    //
    printf("POLYNOMIAL INTERPOLATION\n\n");
    printf("F(x)=sin(x), [0, pi]\n");
    printf("Second degree polynomial is used\n\n");
    
    //
    // Create polynomial interpolant
    //
    n = 3;
    y.setlength(n);
    for(i=0; i<=n-1; i++)
    {
        y(i) = amp::sin<Precision>(amp::ampf<Precision>("0.5")*amp::pi<Precision>()*(amp::ampf<Precision>("1.0")+amp::cos<Precision>(amp::pi<Precision>()*i/(n-1))));
    }
    polint::polynomialbuildcheb2<Precision>(amp::ampf<Precision>(0), amp::pi<Precision>(), y, n, p);
    
    //
    // Output results
    //
    ratint::barycentricdiff2<Precision>(p, amp::ampf<Precision>(0), v, dv, d2v);
    printf("                 P(x)    F(x) \n");
    printf("function       %6.3lf  %6.3lf \n",
        double(amp::ampf<Precision>(ratint::barycentriccalc<Precision>(p, amp::ampf<Precision>(0))).toDouble()),
        double(amp::ampf<Precision>(0).toDouble()));
    printf("d/dx(0)        %6.3lf  %6.3lf \n",
        double(amp::ampf<Precision>(dv).toDouble()),
        double(amp::ampf<Precision>(1).toDouble()));
    printf("d2/dx2(0)      %6.3lf  %6.3lf \n",
        double(amp::ampf<Precision>(d2v).toDouble()),
        double(amp::ampf<Precision>(0).toDouble()));
    printf("\n\n");    
    return 0;
}