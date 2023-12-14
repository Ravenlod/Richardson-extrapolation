#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "odu_NN.h"
#include "func_NN.h"

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
    return currentY;
}
