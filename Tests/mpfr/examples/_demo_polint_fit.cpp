#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "polint.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    int m;
    int n;
    ap::template_1d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > y;
    ap::template_1d_array< amp::ampf<Precision> > w;
    ap::template_1d_array< amp::ampf<Precision> > xc;
    ap::template_1d_array< amp::ampf<Precision> > yc;
    ap::template_1d_array< int > dc;
    polint::polynomialfitreport<Precision> rep;
    int info;
    ratint::barycentricinterpolant<Precision> p;
    int i;
    int j;
    amp::ampf<Precision> a;
    amp::ampf<Precision> b;
    amp::ampf<Precision> v;
    amp::ampf<Precision> dv;


    printf("\n\nFitting exp(2*x) at [-1,+1] by polinomial\n\n");
    printf("Fit type             rms.err max.err    p(0)   dp(0)\n");
    
    //
    // Prepare points
    //
    m = 5;
    a = -1;
    b = +1;
    n = 1000;
    x.setlength(n);
    y.setlength(n);
    w.setlength(n);
    for(i=0; i<=n-1; i++)
    {
        x(i) = a+(b-a)*i/(n-1);
        y(i) = amp::exp<Precision>(2*x(i));
        w(i) = amp::ampf<Precision>("1.0");
    }
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) without constraints
    //
    polint::polynomialfit<Precision>(x, y, n, m, info, p, rep);
    ratint::barycentricdiff1<Precision>(p, amp::ampf<Precision>("0.0"), v, dv);
    printf("Unconstrained        %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()),
        double(amp::ampf<Precision>(rep.maxerror).toDouble()),
        double(amp::ampf<Precision>(v).toDouble()),
        double(amp::ampf<Precision>(dv).toDouble()));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: p(0)=1
    //
    xc.setlength(1);
    yc.setlength(1);
    dc.setlength(1);
    xc(0) = 0;
    yc(0) = 1;
    dc(0) = 0;
    polint::polynomialfitwc<Precision>(x, y, w, n, xc, yc, dc, 1, m, info, p, rep);
    ratint::barycentricdiff1<Precision>(p, amp::ampf<Precision>("0.0"), v, dv);
    printf("Constrained, p(0)=1  %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()),
        double(amp::ampf<Precision>(rep.maxerror).toDouble()),
        double(amp::ampf<Precision>(v).toDouble()),
        double(amp::ampf<Precision>(dv).toDouble()));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: dp(0)=2
    //
    xc.setlength(1);
    yc.setlength(1);
    dc.setlength(1);
    xc(0) = 0;
    yc(0) = 2;
    dc(0) = 1;
    polint::polynomialfitwc<Precision>(x, y, w, n, xc, yc, dc, 1, m, info, p, rep);
    ratint::barycentricdiff1<Precision>(p, amp::ampf<Precision>("0.0"), v, dv);
    printf("Constrained, dp(0)=2 %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()),
        double(amp::ampf<Precision>(rep.maxerror).toDouble()),
        double(amp::ampf<Precision>(v).toDouble()),
        double(amp::ampf<Precision>(dv).toDouble()));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: p(0)=1, dp(0)=2
    //
    xc.setlength(2);
    yc.setlength(2);
    dc.setlength(2);
    xc(0) = 0;
    yc(0) = 1;
    dc(0) = 0;
    xc(1) = 0;
    yc(1) = 2;
    dc(1) = 1;
    polint::polynomialfitwc<Precision>(x, y, w, n, xc, yc, dc, 2, m, info, p, rep);
    ratint::barycentricdiff1<Precision>(p, amp::ampf<Precision>("0.0"), v, dv);
    printf("Constrained, both    %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(amp::ampf<Precision>(rep.rmserror).toDouble()),
        double(amp::ampf<Precision>(rep.maxerror).toDouble()),
        double(amp::ampf<Precision>(v).toDouble()),
        double(amp::ampf<Precision>(dv).toDouble()));
    printf("\n\n");    
    return 0;
}