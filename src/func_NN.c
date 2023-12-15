
#include "func_NN.h"
#include <math.h>

double funk(int i, double x, double* y) 
{
    double result;

    switch (i) 
    {
        case 0:
            // 1 result = sin(x) * y[0] + cos(x) * y[1]; 
            result = 4 * y[0] - 3 * y[1] + 2 * sin(x); 
            break;
        case 1:
            // 1 result = sin(x) * y[1] + cos(x) * y[0];
            result = 2 * y[0] - y[1] - 2 * cos(x);
            break;
        default:
            return NAN;
    }

    return result;
}