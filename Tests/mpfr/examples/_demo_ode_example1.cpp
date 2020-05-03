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
    // y' = y, y(0) = 1.
    //
    // Its solution is well known in academic circles :)
    //
    // No intermediate values are calculated,
    // just starting and final points.
    //
    y.setlength(1);
    y(0) = 1;
    x.setlength(2);
    x(0) = 0;
    x(1) = 1;
    eps = amp::ampf<Precision>("1.0E-4");
    h = amp::ampf<Precision>("0.01");
    odesolver::odesolverrkck<Precision>(y, 1, x, 2, eps, h, state);
    while( odesolver::odesolveriteration<Precision>(state) )
    {
        state.dy(0) = state.y(0);
    }
    odesolver::odesolverresults<Precision>(state, m, x, ytbl, rep);
    printf("    X  Y(X)\n");
    for(i=0; i<=m-1; i++)
    {
        printf("%5.3lf %5.3lf\n",
            double(amp::ampf<Precision>(x(i)).toDouble()),
            double(amp::ampf<Precision>(ytbl(i,0)).toDouble()));
    }    
    return 0;
}