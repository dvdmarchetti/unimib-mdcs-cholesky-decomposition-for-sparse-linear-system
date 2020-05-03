#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "autogk.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    autogk::autogkstate<Precision> state;
    amp::ampf<Precision> v;
    autogk::autogkreport<Precision> rep;
    amp::ampf<Precision> a;
    amp::ampf<Precision> b;
    amp::ampf<Precision> alpha;


    
    //
    // f1(x) = (1+x)*(x-a)^alpha, alpha=-0.3
    // Exact answer is (B-A)^(Alpha+2)/(Alpha+2) + (1+A)*(B-A)^(Alpha+1)/(Alpha+1)
    //
    // This code demonstrates use of the State.XMinusA (State.BMinusX) field.
    //
    // If we try to use State.X instead of State.XMinusA,
    // we will end up dividing by zero! (in 64-bit precision)
    //
    a = amp::ampf<Precision>("1.0");
    b = amp::ampf<Precision>("5.0");
    alpha = -amp::ampf<Precision>("0.9");
    autogk::autogksingular<Precision>(a, b, alpha, amp::ampf<Precision>("0.0"), state);
    while( autogk::autogkiteration<Precision>(state) )
    {
        state.f = amp::pow<Precision>(state.xminusa, alpha)*(1+state.x);
    }
    autogk::autogkresults<Precision>(state, v, rep);
    printf("integral((1+x)*(x-a)^alpha) on [%0.1lf; %0.1lf] = %0.2lf\nExact answer is %0.2lf\n",
        double(amp::ampf<Precision>(a).toDouble()),
        double(amp::ampf<Precision>(b).toDouble()),
        double(amp::ampf<Precision>(v).toDouble()),
        double(amp::ampf<Precision>(amp::pow<Precision>(b-a, alpha+2)/(alpha+2)+(1+a)*amp::pow<Precision>(b-a, alpha+1)/(alpha+1)).toDouble()));    
    return 0;
}