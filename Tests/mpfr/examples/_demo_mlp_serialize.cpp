#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mlpbase.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    mlpbase::multilayerperceptron<Precision> network1;
    mlpbase::multilayerperceptron<Precision> network2;
    mlpbase::multilayerperceptron<Precision> network3;
    ap::template_1d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > y;
    ap::template_1d_array< amp::ampf<Precision> > r;
    int rlen;
    amp::ampf<Precision> v1;
    amp::ampf<Precision> v2;


    
    //
    // Generate two networks filled with small random values.
    // Use MLPSerialize/MLPUnserialize to make network copy.
    //
    mlpbase::mlpcreate0<Precision>(1, 1, network1);
    mlpbase::mlpcreate0<Precision>(1, 1, network2);
    mlpbase::mlpserialize<Precision>(network1, r, rlen);
    mlpbase::mlpunserialize<Precision>(r, network2);
    
    //
    // Now Network1 and Network2 should be identical.
    // Let's demonstrate it.
    //
    printf("Test serialization/unserialization\n");
    x.setlength(1);
    y.setlength(1);
    x(0) = 2*amp::ampf<Precision>::getRandom()-1;
    mlpbase::mlpprocess<Precision>(network1, x, y);
    v1 = y(0);
    printf("Network1(X) = %0.2lf\n",
        double(amp::ampf<Precision>(y(0)).toDouble()));
    mlpbase::mlpprocess<Precision>(network2, x, y);
    v2 = y(0);
    printf("Network2(X) = %0.2lf\n",
        double(amp::ampf<Precision>(y(0)).toDouble()));
    if( v1==v2 )
    {
        printf("Results are equal, OK.\n");
    }
    else
    {
        printf("Results are not equal... Strange...");
    }    
    return 0;
}