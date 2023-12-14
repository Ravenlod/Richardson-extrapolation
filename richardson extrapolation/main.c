#include <stdio.h>
#include <stdlib.h>
#include "odu.h"
#include "funk.h"

int main() 
{
    FILE *input, *output;
    input = fopen("input.txt", "r");
    int n;
    double a,b; //= 0, b = 20;
    double e;
    int k, *k1 = (int*)malloc(sizeof(int));//k=50
    fscanf(input, "%d %lf %lf %lf %d", &n, &a, &b, &e, &k);
    double *y0 = (double* )malloc(sizeof(double) * n);
    double *y0Line = (double* )malloc(sizeof(double) * n);//double y0Line[] = { 1.0, -2.0, -2.0 };
    
    for(int i = 0; i < n; i ++)
    {
        fscanf(input, "%lf", &y0Line[i]);
    }
    for(int i = 0; i < n; i ++)
    {
        y0[i] = y0Line[i];
    }
    printf("####");
    printf("n:%d a:%lf b:%lf e:%lf k:%d\n", n, a, b, e, k);
    output = fopen("output.txt", "w");
    fprintf(output, "n:%d a:%lf b:%lf e:%lf k:%d ", n, a, b, e, k);
    
    double** result, **viewOutput;

    result = (double**)malloc((k + 1) * sizeof(double*));
    for (int i = 0; i <= k; i++) 
    {
        result[i] = (double*)malloc((n + 1) * sizeof(double));
    }

    solveODE(n, a, b, e, k, y0Line, result, k1);

    double h = (b - a) / k;
    double temp_a = a;

    printf("n:%d a:%lf b:%lf e:%lf k:%d K1:%d\n", n, a, b, e, k, *k1);
    fprintf(output, "k1:%d\n", *k1);
    for(int i = 0; i < k + 1; i ++)
    {
        printf("%lf: ",temp_a);
        fprintf(output, "%lf: ",temp_a);
        for(int j = 0; j < n; j ++)
        {
            printf("%lf ", result[i][j]);
            fprintf(output, "%lf ", result[i][j]);
        }
        temp_a += h;
        printf("\n");
        fprintf(output, "\n");
    }

    viewOutput = (double**)malloc((*k1 + 1) * sizeof(double*));
    for (int i = 0; i <= *k1; i++) 
    {
        viewOutput[i] = (double*)malloc((n + 1) * sizeof(double));
    }
    graphViewer(n, a, b, *k1, y0Line, viewOutput);
    
    // printf("*******************************\n");

    // for (int i = 0; i < (*k1 + 1); i++) {
    //     for (int j = 0; j < (n + 1); j++) {
    //         printf("%lf ", viewOutput[i][j]);
    //         fprintf(file, "%lf ", viewOutput[i][j]);
    //     }
    //     printf("\n");
    //     fprintf(file, "\n");
    // }



    printf("*****\n%lf\n", viewOutput[0][0]);
    for(int i = 0; i < *k1 + 1; i ++)
    {
        free(viewOutput[i]);
    }
    free(viewOutput);

    for(int i = 0; i < k + 1; i ++)
    {
        free(result[i]);
    }
    free(result);
    free(k1);
    scanf("%*d");
    return 0;
}

