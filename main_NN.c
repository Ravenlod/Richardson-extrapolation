// main_NN.c
// #define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define EXTRAPOLATION_MATRIX_SIZE 5
#define ELEMENTS_COUNT_PER_SEGMENT 25
//#include "odu_NN.h"
double funk(int i, double x, double* y);
void solveODE(int n, double a, double b, double e, int k, double* y0, double** result, int *k1);
double* solveRunge(int n, double a, double b,  int k, double* y0);
double* GaussElimination(int pos, int n, double **matrix);
double* RichardsonExtrapolation(int n, double segmentBegin, double segmentEnd, int k, double* y0);
void graphViewer(int n, double a, double b, int k, double* y0, double** result, int k_old);


int main() 
{
    int n, k, *k1 = (int*)malloc(sizeof(int));
    double a, b, e, *y0;
    FILE *inputVars = fopen("input.txt", "r");
    fscanf(inputVars, "%d %lf %lf %lf %d", &n, &a, &b, &e, &k);

    y0 = (double*)malloc(n * sizeof(double));
    for(int i = 0; i < n; i ++)
    {
        fscanf(inputVars, "%lf", &y0[i]);
    }

    // double* out = solveRunge(n, a, b, k, y0);
    // printf("Runge: %lf %lf\n", out[0], out[1]);

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

double* RichardsonExtrapolation(int n, double segmentBegin, double segmentEnd, int k, double* y0) 
{
    double h = (segmentEnd - segmentBegin) / k;
    double **leftLine = (double**)malloc(EXTRAPOLATION_MATRIX_SIZE * sizeof(double*));
    double **matrix = (double**)malloc(EXTRAPOLATION_MATRIX_SIZE * sizeof(double*));

    for(int i = 0; i < EXTRAPOLATION_MATRIX_SIZE; i++)
    {
        matrix[i] = (double *)malloc((EXTRAPOLATION_MATRIX_SIZE + 1) * sizeof(double));

        for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j ++)
        {
            matrix[i][j] = exp( pow(2, j) * h);
            //printf("%lf ", matrix[i][j]);
        }

        //printf("\n");
        h /= 2;
    }

    for(int i = 0, k_change = k; i < EXTRAPOLATION_MATRIX_SIZE; i++, k_change *= 2)
    {
        leftLine[i] = solveRunge(n, segmentBegin, segmentEnd, k_change, y0);
        //printf("RE:%lf:\n", y0[0]);
    }

    double *extrapolationResult = (double*)malloc(sizeof(double) * n);

    for(int i = 0; i < n; i ++)
    {
        for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j ++)
        {
            matrix[j][EXTRAPOLATION_MATRIX_SIZE] = leftLine[j][i];
        }
        
        double**tempMatrix = (double**)malloc(EXTRAPOLATION_MATRIX_SIZE * sizeof(double*));
        for (int i = 0; i < EXTRAPOLATION_MATRIX_SIZE; i++) 
        {
            tempMatrix[i] = (double*)malloc((EXTRAPOLATION_MATRIX_SIZE+1) * sizeof(double));
        }
        for(int i = 0; i < EXTRAPOLATION_MATRIX_SIZE; i ++)
        {
            for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE+1; j ++)
            {
                tempMatrix[i][j] = matrix[i][j];
            }
        }

        double *rootLine = GaussElimination(0,EXTRAPOLATION_MATRIX_SIZE,tempMatrix);
        double feedbackCheck = 0;
        for(int j = 0;j < EXTRAPOLATION_MATRIX_SIZE; j++){
            feedbackCheck += matrix[0][j] * rootLine[j];
        }
        
        //printf("feedback check: |%lf|\n" , feedbackCheck - leftLine[0][i]);

        extrapolationResult[i] = 0;
        for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j ++)
        {
            extrapolationResult[i] += rootLine[j];
            
        }
        free(rootLine);
        // printf("RESULT: %lf\n", extrapolationResult[i]);
        
        for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j ++)
        {
            free(tempMatrix[j]);
        }
        free(tempMatrix);
    }
    for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j ++)
    {
        free(matrix[j]);
    }
        free(matrix);
    for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j ++)
    {
        free(leftLine[j]);
    }
        free(leftLine);

    return extrapolationResult;
}

double* GaussElimination(int pos, int n, double **matrix)
{
    double sin, cos;
    double** copyMatrix;

    copyMatrix = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) 
    {
        copyMatrix[i] = (double*)malloc((n + 1) * sizeof(double));
    }


    cos = matrix[pos][pos] / (sqrt(matrix[pos][pos] * matrix[pos][pos] + matrix[pos + 1][pos] * matrix[pos + 1][pos]));
    sin = matrix[pos + 1][pos] * cos / matrix[pos][pos];

    //основной цикл; количество итераций уменьшается с ростом вложености
    for (int i = 0; i < n - 1 - pos; i++) 
    {
        
        //Создаем копию матрицу для преобразования
        for (int m = 0; m < n; m ++) 
        {
            for (int j = 0; j < n + 1; j ++) 
            {
                copyMatrix[m][j] = matrix[m][j];
            }

        }

        //расчет матрицы поворотов
        for (int j = 0; j < n + 1; j++) 
        {
            copyMatrix[pos][j] = matrix[pos][j] * cos + matrix[pos + 1 + i][j] * sin;
            copyMatrix[pos + i + 1][j] = -matrix[pos][j] * sin + matrix[pos + 1 + i][j] * cos;
        }
        
        //Перезаписываем изменения в основную матрицу
        for (int m = 0; m < n; m ++) 
        {
            for (int j = 0; j < n + 1; j++) 
            {
                matrix[m][j] = copyMatrix[m][j];
            }
        }

        //Костыль для ограничения расчётов косинуса и синуса
        if (i + 2 != n - pos) 
        {
            cos = matrix[pos][pos] / (sqrt(matrix[pos][pos] * matrix[pos][pos] + matrix[pos + 2 + i][pos] * matrix[pos + 2 + i][pos]));
            sin = matrix[pos + 2 + i][pos] * cos / matrix[pos][pos];
        }

    }
    for (int i = 0; i < n; i++) 
    {
        free(copyMatrix[i]);
    }
    free(copyMatrix);
    pos += 1;

    //выход из функции
    if (pos == n - 1) 
    {

        double *roots = (double*)malloc(sizeof(double) * n);
        double numerator = matrix[n - 1][n];
        for (int i = n; i > 0; i--)
        {
            roots[i - 1] = numerator/matrix[i - 1][i - 1];
            if(i == 1)
            {
                break;
            }
            numerator = matrix[i - 2][n];
            for(int j = n; j > i - 1; j --)
            {
                numerator -= matrix[i - 2][j - 1] * roots[j - 1]; 
            }
        }
        return roots;
    }

    //вход в рекурсию со сдвигом pos + 1
    else
    {
        double* result = GaussElimination(pos, n, matrix);
        return result;
    }
}

void solveODE(int n, double a, double b, double e, int k, double* y0, double** result, int *k1) 
{
    double  h = (b - a) / k, segmentBegin = a, segmentEnd, accuracy, 
    kNext = EXTRAPOLATION_MATRIX_SIZE, localSegmentBegin = segmentBegin, 
    localSegmentEnd;
    int iterationsCount = 1, localIterationCount, counter = 1;
    double *nextValueLine = (double*)malloc(sizeof(double) * n);
    double *iterationValue;
    int accuracyFlag = 1;
    double **accuracyMatrix = (double**)malloc(sizeof(double*) * k);

    for(int i = 0; i < k; i ++)
    {
        accuracyMatrix[i] = (double*)malloc(sizeof(double) * n);
    }
    for(int i = 0; i < n; i ++)
    {
        result[0][i] = y0[i];
    }

    do
    {
        
        localIterationCount = iterationsCount;

        for(int m = 0; m < k; m ++)
        {
            for(int s = 0; s < 2; s ++)
            {   
                for(int j = 0; j < n; j ++ )
                {
                    nextValueLine[j] = result[m][j];
                }
                localSegmentBegin = segmentBegin;
                for(int i = 0; i < localIterationCount; i ++)
                {
                    
                    localSegmentEnd = localSegmentBegin + h / localIterationCount;
                    iterationValue = RichardsonExtrapolation(n, localSegmentBegin, 
                                        localSegmentEnd, kNext , nextValueLine);
                    localSegmentBegin = localSegmentEnd;
                    for(int j = 0; j < n; j ++ )
                    {
                        nextValueLine[j] = iterationValue[j];
                    }
                    if(i < localIterationCount - 1)
                    {
                        free(iterationValue);
                    } 
                }
                
                if(s == 0)
                {
                    for(int i = 0; i < n; i ++)
                    {
                        accuracyMatrix[m][i] = iterationValue[i];
                    }
                    free(iterationValue);
                }

                localIterationCount *= 2;
            }
            localIterationCount = iterationsCount;
            for(int s = 0; s < n; s ++)
            {
                result[m + 1][s] = iterationValue[s];
                accuracyMatrix[m][s] -= iterationValue[s]; 

            }
            free(iterationValue);
            segmentBegin += h;
        }

        for(int m = 0; m < n; m ++)
        {
            for(int i = 0; i < k; i ++)
            {   
                accuracy += accuracyMatrix[i][m] * accuracyMatrix[i][m];
            }
            accuracy = sqrt(accuracy) / k;
            if(accuracy > e)
            {
                accuracyFlag = 0;
                accuracy = 0;
                iterationsCount *= 2;
            }
            else
            {
                accuracyFlag = 1;
            }
        }
        
    segmentBegin = a;
    } while(accuracyFlag == 0);

    *k1 = iterationsCount * k;
    free(nextValueLine);
    for(int i = 0; i < k; i ++)
    {
        free(accuracyMatrix[i]);
    }
    free(accuracyMatrix);

}

double* solveRunge(int n, double a, double b, int k, double* y0) 
{
    double **fButcherTableau = (double**)malloc(5 * sizeof(double*));

    fButcherTableau[0] = (double*)malloc(3  * sizeof(double));
    fButcherTableau[1] = (double*)malloc(1  * sizeof(double));
    fButcherTableau[2] = (double*)malloc(2  * sizeof(double));
    fButcherTableau[3] = (double*)malloc(3  * sizeof(double));
    fButcherTableau[4] = (double*)malloc(4  * sizeof(double));

    fButcherTableau[0][0] = 0.262;
    fButcherTableau[0][1] = 0.160457849033965;
    fButcherTableau[0][2] = 1;
    fButcherTableau[1][0] = fButcherTableau[0][0];
    fButcherTableau[2][1] = (fButcherTableau[0][1] * (fButcherTableau[0][1] - fButcherTableau[0][0])) / (2 * fButcherTableau[0][0] * (1 - 2 * fButcherTableau[0][0]));
    fButcherTableau[2][0] = fButcherTableau[0][1] - fButcherTableau[2][1];
    fButcherTableau[4][1] = (2 * fButcherTableau[0][1] - 1) / (12 * fButcherTableau[0][0] * (fButcherTableau[0][1] - fButcherTableau[0][0]) * (1 - fButcherTableau[0][0]));
    fButcherTableau[4][2] = (2 * fButcherTableau[0][0] - 1) / (12 * fButcherTableau[0][1] * (fButcherTableau[0][0] - fButcherTableau[0][1]) * (1 - fButcherTableau[0][1]));
    fButcherTableau[4][3] = (6 * fButcherTableau[0][0] *fButcherTableau[0][1] - 4 * fButcherTableau[0][0] - 4 * fButcherTableau[0][1] + 3) / (12 * (1 - fButcherTableau[0][0]) * (1 - fButcherTableau[0][1]));
    fButcherTableau[3][1] = (fButcherTableau[4][1] * (1 - fButcherTableau[0][0]) - fButcherTableau[4][2] * fButcherTableau[2][1]) / fButcherTableau[4][3];
    fButcherTableau[3][2] = (fButcherTableau[4][2] * (1 - fButcherTableau[0][1]) / fButcherTableau[4][3]);
    fButcherTableau[3][0] = 1 - fButcherTableau[3][1] - fButcherTableau[3][2];
    fButcherTableau[4][0] = 1 - fButcherTableau[4][1] - fButcherTableau[4][2] - fButcherTableau[4][3];
    
    //printf("b:%lf %lf %lf %lf %lf %lf\n", fButcherTableau[1][0], fButcherTableau[2][0], 
    //fButcherTableau[2][1], fButcherTableau[3][0], fButcherTableau[3][1], fButcherTableau[3][2]);
    //printf("c: %lf %lf %lf %lf\n", fButcherTableau[4][0], fButcherTableau[4][1], 
    //fButcherTableau[4][2], fButcherTableau[4][3]);

    double currentX = a, *localY = (double*)malloc(n * sizeof(double)), h = (b - a) / ELEMENTS_COUNT_PER_SEGMENT, localX;
    double **kMatrix = (double**)malloc( n * sizeof(double*)), *currentY = (double*)malloc(n * sizeof(double));

    for(int i = 0; i < n; i ++)
    {
        localY[i] = y0[i];
        currentY[i] = y0[i];
        kMatrix[i] = (double*)malloc(4 * sizeof(double));
    }


    for(int i = 0; i < ELEMENTS_COUNT_PER_SEGMENT; i ++)
    {
        localX = currentX;
        for(int j = 0; j < n; j ++)
        {
            localY[j] = currentY[j];
        }
        for(int j = 0; j < 4; j ++)
        {
            for(int s = 0; s < n; s ++)
            {
                kMatrix[s][j] = funk(s, localX, localY) * h;
            }
            if(j < 3)
            {
                for(int m = 0; m < n; m ++)
                {
                    localY[m] = currentY[m];
                    for(int s = 0; s <= j; s ++)
                    {
                        localY[m] += fButcherTableau[j + 1][s] * kMatrix[m][s];
                    }
                }
                localX = currentX + fButcherTableau[0][j] * h;
            }
           
        }

        for(int j = 0; j < n; j ++)
        {

            currentY[j] += fButcherTableau[4][0] * kMatrix[j][0] + fButcherTableau[4][1] * kMatrix[j][1] + 
            fButcherTableau[4][2] * kMatrix[j][2] + fButcherTableau[4][3] * kMatrix[j][3];         
        }
        currentX += h;
    }
    //printf("%lf %lf\n", currentY[0], currentY[1]);
    //printf("||%lf||\n", currentX);
    //printf("%lf, %lf %lf a:%lf b:%lf\n", h, y0[0], y0[1]);
    return currentY;
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
    
    printf("h = %lf\n", h);
    printf("hCount = %d\n", step);
    printf("minimalStepLength = %lf\n", minimalStepLength);
    
    for(int m = 0; m < n; m++)
    {
        const char fileNameConst[] = {'l', 'i', 'n', 'k', (char)(m + 48), '.', 's', 'v', 'g', '\0'};
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
