
#include "func_NN.h"
#include <math.h>

double funk(int i, double x, double* y) 
{
    double result;

    switch (i) 
    {
        case 0:
            result = sin(x) * y[0] + cos(x) * y[1]; 
            break;
        case 1:
            result = sin(x) * y[1] + cos(x) * y[0];
            break;
        default:
            return NAN;
    }

    return result;
}