#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mincg.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    int n;
    mincg::mincgstate<Precision> state;
    mincg::mincgreport<Precision> rep;
    ap::template_1d_array< amp::ampf<Precision> > s;
    amp::ampf<Precision> x;
    amp::ampf<Precision> y;


    
    //
    // Function minimized:
    //     F = (x-1)^4 + (y-x)^2
    // N = 2 - task dimension.
    //
    n = 2;
    s.setlength(2);
    s(0) = 10;
    s(1) = 11;
    mincg::mincgcreate<Precision>(n, s, state);
    mincg::mincgsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.00001"), 0);
    mincg::mincgsetxrep<Precision>(state, true);
    printf("\n\nF = (x-1)^4 + (y-x)^2\n");
    printf("OPTIMIZATION STARTED\n");
    while( mincg::mincgiteration<Precision>(state) )
    {
        if( state.needfg )
        {
            x = state.x(0);
            y = state.x(1);
            state.f = amp::sqr<Precision>(amp::sqr<Precision>(x-1))+amp::sqr<Precision>(y-x);
            state.g(0) = 4*amp::sqr<Precision>(x-1)*(x-1)+2*(x-y);
            state.g(1) = 2*(y-x);
        }
        if( state.xupdated )
        {
            printf("    F(%8.5lf,%8.5lf)=%0.5lf\n",
                double(amp::ampf<Precision>(state.x(0)).toDouble()),
                double(amp::ampf<Precision>(state.x(1)).toDouble()),
                double(amp::ampf<Precision>(state.f).toDouble()));
        }
    }
    printf("OPTIMIZATION STOPPED\n");
    mincg::mincgresults<Precision>(state, s, rep);
    
    //
    // output results
    //
    printf("X = %4.2lf (should be 1.00)\n",
        double(amp::ampf<Precision>(s(0)).toDouble()));
    printf("Y = %4.2lf (should be 1.00)\n\n\n",
        double(amp::ampf<Precision>(s(1)).toDouble()));    
    return 0;
}