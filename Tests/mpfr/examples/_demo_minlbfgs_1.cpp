#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "minlbfgs.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    int n;
    int m;
    minlbfgs::minlbfgsstate<Precision> state;
    minlbfgs::minlbfgsreport<Precision> rep;
    ap::template_1d_array< amp::ampf<Precision> > s;
    amp::ampf<Precision> x;
    amp::ampf<Precision> y;


    
    //
    // Function minimized:
    //     F = exp(x-1) + exp(1-x) + (y-x)^2
    // N = 2 - task dimension
    // M = 1 - build tank-1 model
    //
    n = 2;
    m = 1;
    s.setlength(2);
    s(0) = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
    s(1) = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
    minlbfgs::minlbfgscreate<Precision>(n, m, s, state);
    minlbfgs::minlbfgssetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0001"), 0);
    while( minlbfgs::minlbfgsiteration<Precision>(state) )
    {
        if( state.needfg )
        {
            x = state.x(0);
            y = state.x(1);
            state.f = amp::exp<Precision>(x-1)+amp::exp<Precision>(1-x)+amp::sqr<Precision>(y-x);
            state.g(0) = amp::exp<Precision>(x-1)-amp::exp<Precision>(1-x)+2*(x-y);
            state.g(1) = 2*(y-x);
        }
    }
    minlbfgs::minlbfgsresults<Precision>(state, s, rep);
    
    //
    // output results
    //
    printf("\n\nF = exp(x-1) + exp(1-x) + (y-x)^2\n");
    printf("X = %4.2lf (should be 1.00)\n",
        double(amp::ampf<Precision>(s(0)).toDouble()));
    printf("Y = %4.2lf (should be 1.00)\n\n\n",
        double(amp::ampf<Precision>(s(1)).toDouble()));    
    return 0;
}