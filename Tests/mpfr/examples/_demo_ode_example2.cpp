#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "odesolver.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    ap::template_1d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > y;
    ap::template_2d_array< amp::ampf<Precision> > ytbl;
    amp::ampf<Precision> eps;
    amp::ampf<Precision> h;
    int m;
    int i;
    odesolver::odesolverstate<Precision> state;
    odesolver::odesolverreport<Precision> rep;


    
    //
    // ODESolver unit is used to solve simple ODE:
    // y'' = -y, y(0) = 0, y'(0)=1.
    //
    // This ODE may be written as first-order system:
    // y' =  z
    // z' = -y
    //
    // Its solution is well known in academic circles :)
    //
    // Three intermediate values are calculated,
    // plus starting and final points.
    //
    y.setlength(2);
    y(0) = 0;
    y(1) = 1;
    x.setlength(5);
    x(0) = amp::pi<Precision>()*0/4;
    x(1) = amp::pi<Precision>()*1/4;
    x(2) = amp::pi<Precision>()*2/4;
    x(3) = amp::pi<Precision>()*3/4;
    x(4) = amp::pi<Precision>()*4/4;
    eps = amp::ampf<Precision>("1.0E-8");
    h = amp::ampf<Precision>("0.01");
    odesolver::odesolverrkck<Precision>(y, 2, x, 5, eps, h, state);
    while( odesolver::odesolveriteration<Precision>(state) )
    {
        state.dy(0) = state.y(1);
        state.dy(1) = -state.y(0);
    }
    odesolver::odesolverresults<Precision>(state, m, x, ytbl, rep);
    printf("     X   Y(X)     Error\n");
    for(i=0; i<=m-1; i++)
    {
        printf("%6.3lf %6.3lf  %8.1le\n",
            double(amp::ampf<Precision>(x(i)).toDouble()),
            double(amp::ampf<Precision>(ytbl(i,0)).toDouble()),
            double(amp::ampf<Precision>(amp::abs<Precision>(ytbl(i,0)-amp::sin<Precision>(x(i)))).toDouble()));
    }    
    return 0;
}