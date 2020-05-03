#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "minlm.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    minlm::minlmstate<Precision> state;
    minlm::minlmreport<Precision> rep;
    int i;
    ap::template_1d_array< amp::ampf<Precision> > s;
    ap::template_1d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > y;
    amp::ampf<Precision> fi;
    int n;
    int m;


    
    //
    // Example of solving polynomial approximation task using FJ scheme.
    //
    // Data points:
    //     xi are random numbers from [-1,+1],
    //
    // Function being fitted:
    //     yi = exp(xi) - sin(xi) - x^3/3
    //
    // Function being minimized:
    //     F(a,b,c) =
    //         (a + b*x0 + c*x0^2 - y0)^2 +
    //         (a + b*x1 + c*x1^2 - y1)^2 + ...
    //
    n = 3;
    s.setlength(n);
    for(i=0; i<=n-1; i++)
    {
        s(i) = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
    }
    m = 100;
    x.setlength(m);
    y.setlength(m);
    for(i=0; i<=m-1; i++)
    {
        x(i) = amp::ampf<Precision>(2*i)/(amp::ampf<Precision>(m-1))-1;
        y(i) = amp::exp<Precision>(x(i))-amp::sin<Precision>(x(i))-x(i)*x(i)*x(i)/3;
    }
    
    //
    // Now S stores starting point, X and Y store points being fitted.
    //
    minlm::minlmcreatefj<Precision>(n, m, s, state);
    minlm::minlmsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.001"), 0);
    while( minlm::minlmiteration<Precision>(state) )
    {
        if( state.needf )
        {
            state.f = 0;
        }
        for(i=0; i<=m-1; i++)
        {
            
            //
            // "a" is stored in State.X[0]
            // "b" - State.X[1]
            // "c" - State.X[2]
            //
            fi = state.x(0)+state.x(1)*x(i)+state.x(2)*amp::sqr<Precision>(x(i))-y(i);
            if( state.needf )
            {
                
                //
                // F is equal to sum of fi squared.
                //
                state.f = state.f+amp::sqr<Precision>(fi);
            }
            if( state.needfij )
            {
                
                //
                // Fi
                //
                state.fi(i) = fi;
                
                //
                // dFi/da
                //
                state.j(i,0) = 1;
                
                //
                // dFi/db
                //
                state.j(i,1) = x(i);
                
                //
                // dFi/dc
                //
                state.j(i,2) = amp::sqr<Precision>(x(i));
            }
        }
    }
    minlm::minlmresults<Precision>(state, s, rep);
    
    //
    // output results
    //
    printf("A = %4.2lf\n",
        double(amp::ampf<Precision>(s(0)).toDouble()));
    printf("B = %4.2lf\n",
        double(amp::ampf<Precision>(s(1)).toDouble()));
    printf("C = %4.2lf\n",
        double(amp::ampf<Precision>(s(2)).toDouble()));
    printf("TerminationType = %0ld (should be 2 - stopping when step is small enough)\n",
        long(rep.terminationtype));    
    return 0;
}