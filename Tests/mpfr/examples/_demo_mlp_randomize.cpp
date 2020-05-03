#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mlpbase.h"

int main(int argc, char **argv)
{
    const int Precision = 128;
    
    mlpbase::multilayerperceptron<Precision> net;


    mlpbase::mlpcreate0<Precision>(2, 1, net);
    mlpbase::mlprandomize<Precision>(net);    
    return 0;
}