#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "rcond.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    int n;
    int i;
    int j;
    amp::ampf<Precision> c1;
    amp::ampf<Precision> x;
    ap::template_2d_array< amp::ampf<Precision> > a;


    printf("                 CONDITION NUMBERS\n");
    printf("OF VANDERMONDE AND CHEBYSHEV INTERPOLATION MATRICES\n\n");
    printf("    VANDERMONDE   CHEBYSHEV\n");
    printf("  N      1-norm      1-norm\n");
    for(n=2; n<=14; n++)
    {
        a.setlength(n, n);
        printf("%3ld",
            long(n));
        
        //
        // Vandermone matrix
        //
        for(i=0; i<=n-1; i++)
        {
            x = amp::ampf<Precision>(2*i)/(amp::ampf<Precision>(n-1))-1;
            a(i,0) = 1;
            for(j=1; j<=n-1; j++)
            {
                a(i,j) = a(i,j-1)*x;
            }
        }
        c1 = 1/rcond::rmatrixrcond1<Precision>(a, n);
        printf(" %11.1lf",
            double(amp::ampf<Precision>(c1).toDouble()));
        
        //
        // Chebyshev interpolation matrix
        //
        for(i=0; i<=n-1; i++)
        {
            x = amp::ampf<Precision>(2*i)/(amp::ampf<Precision>(n-1))-1;
            a(i,0) = 1;
            if( n>=2 )
            {
                a(i,1) = x;
            }
            for(j=2; j<=n-1; j++)
            {
                a(i,j) = 2*x*a(i,j-1)-a(i,j-2);
            }
        }
        c1 = 1/rcond::rmatrixrcond1<Precision>(a, n);
        printf(" %11.1lf\n",
            double(amp::ampf<Precision>(c1).toDouble()));
    }    
    return 0;
}