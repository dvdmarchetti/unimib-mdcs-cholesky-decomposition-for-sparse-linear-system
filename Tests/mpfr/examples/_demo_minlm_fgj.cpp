#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "minlm.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    minlm::minlmstate<Precision> state;
    minlm::minlmreport<Precision> rep;
    ap::template_1d_array< amp::ampf<Precision> > s;
    amp::ampf<Precision> x;
    amp::ampf<Precision> y;


    
    //
    // Example of solving simple task using FGJ scheme.
    //
    // Function minimized:
    //     F = (x-2*y)^2 + (x-2)^2 + (y-1)^2
    // exact solution is (2,1).
    //
    s.setlength(2);
    s(0) = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
    s(1) = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
    minlm::minlmcreatefgj<Precision>(2, 3, s, state);
    minlm::minlmsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.001"), 0);
    while( minlm::minlmiteration<Precision>(state) )
    {
        x = state.x(0);
        y = state.x(1);
        if( state.needf )
        {
            state.f = amp::sqr<Precision>(x-2*y)+amp::sqr<Precision>(x-2)+amp::sqr<Precision>(y-1);
        }
        if( state.needfg )
        {
            state.f = amp::sqr<Precision>(x-2*y)+amp::sqr<Precision>(x-2)+amp::sqr<Precision>(y-1);
            state.g(0) = 2*(x-2*y)+2*(x-2)+0;
            state.g(1) = -4*(x-2*y)+0+2*(y-1);
        }
        if( state.needfij )
        {
            state.fi(0) = x-2*y;
            state.fi(1) = x-2;
            state.fi(2) = y-1;
            state.j(0,0) = 1;
            state.j(0,1) = -2;
            state.j(1,0) = 1;
            state.j(1,1) = 0;
            state.j(2,0) = 0;
            state.j(2,1) = 1;
        }
    }
    minlm::minlmresults<Precision>(state, s, rep);
    
    //
    // output results
    //
    printf("X = %4.2lf (correct value - 2.00)\n",
        double(amp::ampf<Precision>(s(0)).toDouble()));
    printf("Y = %4.2lf (correct value - 1.00)\n",
        double(amp::ampf<Precision>(s(1)).toDouble()));
    printf("TerminationType = %0ld (should be 2 - stopping when step is small enough)\n",
        long(rep.terminationtype));
    printf("NFunc = %0ld\n",
        long(rep.nfunc));
    printf("NJac  = %0ld\n",
        long(rep.njac));
    printf("NGrad = %0ld\n",
        long(rep.ngrad));
    printf("NHess = %0ld\n",
        long(rep.nhess));    
    return 0;
}