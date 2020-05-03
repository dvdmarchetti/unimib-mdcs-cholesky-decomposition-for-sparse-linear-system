#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mlpbase.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    mlpbase::multilayerperceptron<Precision> net;
    ap::template_1d_array< amp::ampf<Precision> > x;
    ap::template_1d_array< amp::ampf<Precision> > y;


    
    //
    // regression task with 2 inputs (independent variables)
    // and 2 outputs (dependent variables).
    //
    // network weights are initialized with small random values.
    //
    mlpbase::mlpcreate0<Precision>(2, 2, net);
    x.setlength(2);
    y.setlength(2);
    x(0) = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
    x(1) = amp::ampf<Precision>::getRandom()-amp::ampf<Precision>("0.5");
    mlpbase::mlpprocess<Precision>(net, x, y);
    printf("Regression task\n");
    printf("IN[0]  = %5.2lf\n",
        double(amp::ampf<Precision>(x(0)).toDouble()));
    printf("IN[1]  = %5.2lf\n",
        double(amp::ampf<Precision>(x(1)).toDouble()));
    printf("OUT[0] = %5.2lf\n",
        double(amp::ampf<Precision>(y(0)).toDouble()));
    printf("OUT[1] = %5.2lf\n",
        double(amp::ampf<Precision>(y(1)).toDouble()));    
    return 0;
}