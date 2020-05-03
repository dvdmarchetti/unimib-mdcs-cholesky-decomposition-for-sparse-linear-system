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


    
    //
    // f(x) = x*sin(x), integrated at [-pi, pi].
    // Exact answer is 2*pi
    //
    autogk::autogksmooth<Precision>(-amp::pi<Precision>(), +amp::pi<Precision>(), state);
    while( autogk::autogkiteration<Precision>(state) )
    {
        state.f = state.x*amp::sin<Precision>(state.x);
    }
    autogk::autogkresults<Precision>(state, v, rep);
    printf("integral(x*sin(x),-pi,+pi) = %0.2lf\nExact answer is %0.2lf\n",
        double(amp::ampf<Precision>(v).toDouble()),
        double(amp::ampf<Precision>(2*amp::pi<Precision>()).toDouble()));    
    return 0;
}