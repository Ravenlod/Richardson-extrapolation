// main_NN.c
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "odu_NN.h"

void graphViewer(int n, double a, double b, int k, double* y0, double** result, int k_old);

int main() 
{
    int n, k, *k1 = (int*)malloc(sizeof(int));
    double a, b, e, *y0;
    FILE *inputVars = fopen("var/input.txt", "r");
    FILE *output = fopen("var/output.txt", "w");

    fscanf(inputVars, "%d %lf %lf %lf %d", &n, &a, &b, &e, &k);

    y0 = (double*)malloc(n * sizeof(double));
    for(int i = 0; i < n; i ++)
    {
        fscanf(inputVars, "%lf", &y0[i]);
    }

    double** result, **viewOutput;

    result = (double**)malloc((k + 1) * sizeof(double*));
    for (int i = 0; i <= k; i++) 
    {
        result[i] = (double*)malloc((n + 1) * sizeof(double));
    }

    solveODE(n, a, b, e, k, y0, result, k1);

    double h = (b - a) / k;
    double temp_a = a;

    printf("n:%d a:%lf b:%lf e:%lf k:%d K1:%d\n", n, a, b, e, k, *k1);
    fprintf(output, "n:%d a:%lf b:%lf e:%lf k:%d K1:%d\n", n, a, b, e, k, *k1);
    // for(int i = 0; i < k + 1; i ++)
    // {
    //     printf("%lf: ",temp_a);
    //     for(int j = 0; j < n; j ++)
    //     {
    //         printf("%lf ", result[i][j]);
    //     }
    //     temp_a += h;
    //     printf("\n");
    // }

    viewOutput = (double**)malloc((*k1 + 1) * sizeof(double*));
    for (int i = 0; i <= *k1; i++) 
    {
        viewOutput[i] = (double*)malloc((n + 1) * sizeof(double));
    }
    graphViewer(n, a, b, *k1, y0, viewOutput, k);
    //printf("%lf|%lf", y0[0], y0[1]);


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
    free(y0);
    return 0;
}

void graphViewer(int n, double a, double b, int k, double* y0, double** result, int k_old)
{
    double h = (b - a) / k, segmentBegin = a, segmentEnd, temp_a = a;
    double *iterValue = (double*)malloc(sizeof(double) * n), *nextValue;
    
    for(int i = 0; i < n; i++)
    {
        iterValue[i] = y0[i];
        result[0][i] = y0[i];
    }
    
    for(int i = 0; i < k; i++)
    {
        segmentEnd = segmentBegin + h;
        nextValue = RichardsonExtrapolation(n, segmentBegin, segmentEnd, 
        ELEMENTS_COUNT_PER_SEGMENT, iterValue);
        segmentBegin = segmentEnd;

        for(int j = 0; j < n; j++)
        {
            result[i + 1][j] = nextValue[j];
        }
        free(iterValue);
        iterValue = nextValue;

    }


    double *max = (double*)malloc(n * sizeof(double));
    double *min = (double*)malloc(n * sizeof(double));
    int width = k_old * 16, height = k_old * 9;
    
    double xAxisRelativePosition = 0.1, yAxisRelativePosition = 0.05, indentRelative = 0.01;
    int xAxisAbsolutePosition = xAxisRelativePosition * width;
    int yAxisAbsolutePosition = height - yAxisRelativePosition * height;
    int indentXAbsolutePosition = indentRelative * width;
    int indentYAbsolutePosition = indentRelative * height;

    double minimalStepLength = (b - a) * 2 / width;
    int step = floor(minimalStepLength / h);
    if(step == 0)
    {
        step = 1;
    }

    for(int i = 0; i < n; i++)
    {
        max[i] = result[0][i];
        min[i] = result[0][i];
    }

    for (int i = 0; i < k + 1; i++) 
    {  
        for (int j = 0; j < n; j++) 
        {
            if (result[i][j] > max[j]) 
            {
                max[j] = result[i][j];
            }
            if (result[i][j] < min[j])
            {
                min[j] = result[i][j];
            }
        }    
    }
    
    // printf("h = %lf\n", h);
    // printf("hCount = %d\n", step);
    // printf("minimalStepLength = %lf\n", minimalStepLength);
    
    for(int m = 0; m < n; m++)
    {
        const char fileNameConst[] = {'o', 'p', 't', '/','l', 'i', 'n', 'k', (char)(m + 48), '.', 's', 'v', 'g', '\0'};
        FILE *svgFile = fopen(fileNameConst, "w");

        int triangleUp[][2]  = {{0 - 5, 15}, {5 - 5, 0}, {10 - 5, 15}};
        int triangleRight[][2]  = {{0, 0 - 5}, {15, 5 - 5}, {0, 10 - 5}};

        for(int i = 0; i < 3; i++)
        {
            triangleUp[i][0] = triangleUp[i][0] + xAxisAbsolutePosition;
            triangleUp[i][1] = triangleUp[i][1] + height - yAxisAbsolutePosition - 15;
            triangleRight[i][0] = triangleRight[i][0] + width - xAxisAbsolutePosition;
            triangleRight[i][1] = triangleRight[i][1] + yAxisAbsolutePosition;
        }

        fprintf(svgFile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(svgFile, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n", width, height);

        //X Axis
        fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" />\n", 
                xAxisAbsolutePosition, yAxisAbsolutePosition, width - xAxisAbsolutePosition, yAxisAbsolutePosition);
        //Y Axis
        fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" />\n", 
                xAxisAbsolutePosition, height - yAxisAbsolutePosition, xAxisAbsolutePosition, yAxisAbsolutePosition);
        
        fprintf(svgFile, "<polygon points=\"%d,%d %d,%d %d,%d\"  fill=\"black\" /> \n", 
            triangleUp[0][0], triangleUp[0][1], triangleUp[1][0], triangleUp[1][1], triangleUp[2][0], triangleUp[2][1]); 
        
        fprintf(svgFile, "<polygon points=\"%d,%d %d,%d %d,%d\"  fill=\"black\" /> \n", 
                triangleRight[0][0], triangleRight[0][1], triangleRight[1][0], triangleRight[1][1], triangleRight[2][0], triangleRight[2][1]);

        //Y max
        fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" />\n",
                xAxisAbsolutePosition - 5, (int)(yAxisRelativePosition * height) + indentYAbsolutePosition,
                xAxisAbsolutePosition + 5, (int)(yAxisRelativePosition * height) + indentYAbsolutePosition);
        fprintf(svgFile, "<text x=\"%d\" y=\"%d\" font-size=\"4\" fill=\"black\">%.1lf</text>\n",
                xAxisAbsolutePosition - 15, (int)(yAxisRelativePosition * height) + indentYAbsolutePosition + 1,
                max[m]);
        //Y min
        fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" />\n",
                xAxisAbsolutePosition - 5, yAxisAbsolutePosition - indentYAbsolutePosition,
                xAxisAbsolutePosition + 5, yAxisAbsolutePosition - indentYAbsolutePosition);
        fprintf(svgFile, "<text x=\"%d\" y=\"%d\" font-size=\"4\" fill=\"black\">%.1lf</text>\n",
            xAxisAbsolutePosition - 15, yAxisAbsolutePosition - indentYAbsolutePosition + 1,
            min[m]);
        
        for(int i = 0; i < k_old + 1; i ++)
        {
            fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" stroke-width=\"1\" />\n",
                    xAxisAbsolutePosition + indentXAbsolutePosition + (width - 2 * (xAxisAbsolutePosition + indentXAbsolutePosition) ) * (i) / k_old, 
                    yAxisAbsolutePosition - 5,
                    xAxisAbsolutePosition + indentXAbsolutePosition + (width - 2 * (xAxisAbsolutePosition + indentXAbsolutePosition) ) * (i) / k_old, 
                    yAxisAbsolutePosition + 5);
            fprintf(svgFile, "<text x=\"%d\" y=\"%d\" fill=\"black\" font-size=\"4\">%.1lf</text>\n",
                    xAxisAbsolutePosition + indentXAbsolutePosition - 3 + (width - 2 * (xAxisAbsolutePosition + indentXAbsolutePosition) ) * (i) / k_old, 
                    yAxisAbsolutePosition + 15,
                    a + i * (b - a) / k_old);
        }
        for(int i = 0; i < k / step; i ++)
        {
            fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"salmon\" stroke-width=\"1\" />\n",
                    xAxisAbsolutePosition + indentXAbsolutePosition + (width - 2 * (xAxisAbsolutePosition + indentXAbsolutePosition) ) * (i) * step / k, 
                    (int)floor(yAxisAbsolutePosition - indentYAbsolutePosition - (result[i * step][m] - min[m]) / 
                    fabs(max[m] - min[m]) * (yAxisAbsolutePosition - yAxisRelativePosition * height - 2 * indentYAbsolutePosition)),
                    xAxisAbsolutePosition + indentXAbsolutePosition + (width - 2 * (xAxisAbsolutePosition + indentXAbsolutePosition) ) * (i + 1) * step / k, 
                    (int)floor(yAxisAbsolutePosition - indentYAbsolutePosition - (result[(i + 1) * step][m] - min[m]) / 
                    fabs(max[m] - min[m]) * (yAxisAbsolutePosition - yAxisRelativePosition * height - 2 * indentYAbsolutePosition)));
        }
        fprintf(svgFile, "</svg>\n");



        // for(int i = 0; i < k + 1; i++)
        // {
        //     //printf("%lf ",temp_a);
        //     //fprintf(graphOutput, "%lf ", temp_a);
        //     temp_a += h;

        //     for(int j = 0; j < n; j++)
        //     {
        //         //printf("%lf ", result[i][j]);
        //         fprintf(graphOutput, "%lf ", result[i][j]);
        //     }
        //     //printf("\n");
        //     fprintf(graphOutput, "\n");
        // }
        fclose(svgFile);
    }
    
    free(min);
    free(max);
    free(iterValue);
}
