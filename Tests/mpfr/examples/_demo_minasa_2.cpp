#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "minasa.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    int n;
    int i;
    minasa::minasastate<Precision> state;
    minasa::minasareport<Precision> rep;
    ap::template_1d_array< amp::ampf<Precision> > s;
    ap::template_1d_array< amp::ampf<Precision> > bndl;
    ap::template_1d_array< amp::ampf<Precision> > bndu;
    amp::ampf<Precision> x;
    amp::ampf<Precision> y;
    amp::ampf<Precision> z;


    
    //
    // Function being minimized:
    //     F = x+4y+9z subject to 0<=x<=1, 0<=y<=1, 0<=z<=1.
    //
    // Take a look at MinASASetStpMax() - it restricts step length by
    // a small value, so we can see the current point traveling through
    // a feasible set, sticking to its bounds.
    //
    n = 3;
    s.setlength(n);
    bndl.setlength(n);
    bndu.setlength(n);
    for(i=0; i<=n-1; i++)
    {
        s(i) = 1;
        bndl(i) = 0;
        bndu(i) = 1;
    }
    minasa::minasacreate<Precision>(n, s, bndl, bndu, state);
    minasa::minasasetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.00001"), 0);
    minasa::minasasetxrep<Precision>(state, true);
    minasa::minasasetstpmax<Precision>(state, amp::ampf<Precision>("0.2"));
    printf("\n\nF = x+4y+9z subject to 0<=x<=1, 0<=y<=1, 0<=z<=1\n");
    printf("OPTIMIZATION STARTED\n");
    while( minasa::minasaiteration<Precision>(state) )
    {
        if( state.needfg )
        {
            x = state.x(0);
            y = state.x(1);
            z = state.x(2);
            state.f = x+4*y+9*z;
            state.g(0) = 1;
            state.g(1) = 4;
            state.g(2) = 9;
        }
        if( state.xupdated )
        {
            printf("    F(%4.2lf, %4.2lf, %4.2lf) = %0.3lf\n",
                double(amp::ampf<Precision>(state.x(0)).toDouble()),
                double(amp::ampf<Precision>(state.x(1)).toDouble()),
                double(amp::ampf<Precision>(state.x(2)).toDouble()),
                double(amp::ampf<Precision>(state.f).toDouble()));
        }
    }
    printf("OPTIMIZATION STOPPED\n");
    minasa::minasaresults<Precision>(state, s, rep);
    
    //
    // output results
    //
    printf("X = %4.2lf (should be 0.00)\n",
        double(amp::ampf<Precision>(s(0)).toDouble()));
    printf("Y = %4.2lf (should be 0.00)\n",
        double(amp::ampf<Precision>(s(1)).toDouble()));
    printf("Z = %4.2lf (should be 0.00)\n\n\n",
        double(amp::ampf<Precision>(s(2)).toDouble()));    
    return 0;
}