// main_NN.c
// #define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define EXTRAPOLATION_MATRIX_SIZE 5
#define ELEMENTS_COUNT_PER_SEGMENT 25
#define K_MULTIPLIER 2

//#include "odu_NN.h"
double funk(int i, double x, double* y);
double** solveODE(int n, double a, double b, double e, int k, double* y0, double** result);
double* solveRunge(int n, double a, double b,  int k, double* y0);
double* GaussElimination(int pos, int n, double **matrix);
double* RichardsonExtrapolation(int n, double segmentBegin, double segmentEnd, int k, double* y0);



struct Runge
{
    int k;
    double a, b, a2, a3, a4, b21, b31, b32, b41, b42, b43, c1, c2, c3, c4;
};

int main() 
{
    int n = 1;
    double a = 0, b = 5;
    double e = 0.0000001;
    int k = 10;
    double y0Line[] = { 1 };;
    double** result;

    result = (double**)malloc((k + 1) * sizeof(double*));

    for (int i = 0; i <= k; i++) 
    {
        result[i] = (double*)malloc((n + 1) * sizeof(double));
    }



    // Чтение входных данных из файла
    FILE* input = fopen("input.txt", "r");
    if (input == NULL) 
    {
        fprintf(stderr, "Ошибка открытия файла входных данных.\n");
        return 1;
    }

    /* fscanf(input, "%d %lf %lf %lf %d", &n, &a, &b, &e, &k);

     y0 = (double *)malloc(n * sizeof(double));
     if (y0 == NULL) {
         fprintf(stderr, "Ошибка выделения памяти.\n");
         fclose(input);
         return 1;
     }
     fscanf(input, "%lf", &y0[0]);
    */
    double **output = solveODE(n, a, b, e, k, y0Line, result);

    for(int i = 0; i < k + 1; i ++)
    {
        for(int j = 0; j < n; j ++)
        {
            printf("%lf ", output[i][j]);
        }
        printf("\n");
    }

    for(int i = 0; i < k + 1; i ++)
    {
        free(output[i]);
    }
    free(output);
    return 0;
}

double funk(int i, double x, double* y) 
{
    double result;

    switch (i) 
    {
        // case 0:
        //     result = x * y[0] + x * y[1];  // Пример: y'[0] = x * y[0]
        //     break;
        // case 1:
        //     result = x * y[1] + x * y[0];  // Пример: y'[1] = sin(x) * y[1]
        //     break;

         case 0:
            //result = sin(x) * y[0] + cos(x) * y[1]; 
            result = ((3 * x + y[0] * y[0] * y[0] - 1)/y[0])*((3 * x + y[0] * y[0] * y[0] - 1)/y[0]);
            break;
        /* case 1:
            result = sin(x) * y[1] + cos(x) * y[0];  // Пример: y'[1] = sin(x) * y[1]
            break; */
        default:
            return NAN;
    }

    return result;
}

double* RichardsonExtrapolation(int n, double segmentBegin, double segmentEnd, int k, double* y0) 
{
    double h = (segmentEnd - segmentBegin) / k;
    double **leftLine = (double**)malloc(k * sizeof(double*));
    double **transportedLeftLine = (double**)malloc(n * sizeof(double*));
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
    }


    // for(int j = 0; j < n; j++)
    // {
    // printf("%lf ", leftLine[EXTRAPOLATION_MATRIX_SIZE - 1][j]);
        
    // }
    // printf("|\n");

    


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
        for(int j=0;j<EXTRAPOLATION_MATRIX_SIZE;j++){
            feedbackCheck += matrix[0][j] * rootLine[j];
        }
        
        //printf("feedback check: |%lf|\n" , feedbackCheck - leftLine[0][i]);

        extrapolationResult[i] = 0;
        for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j ++)
        {
            extrapolationResult[i] += rootLine[j];
            
        }

        // printf("RESULT: %lf\n", extrapolationResult[i]);
        
        for(int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j ++)
        {
            free(tempMatrix[j]);
        }
        free(tempMatrix);
    }
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

        // for(int i = 0; i < n; i++)
        // {
        //     for(int j = 0; j < n + 1; j ++)
        //     {
        //         printf("%lf ", matrix[i][j]);
        //     }
        //     printf("\n");

        // }

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

double** solveODE(int n, double a, double b, double e, int k, double* y0, double** result) 
{
    double  h = (b - a) / k, segmentBegin = a, segmentEnd, accuracy, 
    kNext = EXTRAPOLATION_MATRIX_SIZE, localSegmentBegin = segmentBegin, 
    localSegmentEnd;
    int iterationsCount = 1, localIterationCount;
    double *nextValueLine = (double*)malloc(sizeof(double) * n);
    double *iterationValue = (double*)malloc(sizeof(double) * n);

    int accuracyFlag = 1;
    double **accuracyMatrix = (double**)malloc(sizeof(double*) * k);

    double **resultedOutput = (double**)malloc(sizeof(double*) * (k + 1));

    resultedOutput[0] = (double* )malloc(sizeof(double) * n);
    for(int i = 0; i < k; i ++)
    {
        resultedOutput[i + 1] = (double* )malloc(sizeof(double) * n);
        accuracyMatrix[i] = (double*)malloc(sizeof(double) * n);
    }
    for(int i = 0; i < n; i ++)
    {
        resultedOutput[0][i] = y0[i];
    }

    do
    {
        printf("NEW ITERATION\n");
        localIterationCount = iterationsCount;

        for(int m = 0; m < k; m ++)
        {
            for(int s = 0; s < 2; s ++)
            {   
                for(int j = 0; j < n; j ++ )
                {
                    nextValueLine[j] = resultedOutput[m][j];
                    //printf("%lf ", resultedOutput[m][j]);
                }
                //printf("\n___START___\n");
                localSegmentBegin = segmentBegin;
                for(int i = 0; i < localIterationCount; i ++)
                {
                    
                    localSegmentEnd = localSegmentBegin + h / localIterationCount;
                    iterationValue = RichardsonExtrapolation(n, localSegmentBegin, 
                                        localSegmentEnd, kNext , nextValueLine);
                    //kNext *= 2;
                    localSegmentBegin = localSegmentEnd;
                    for(int j = 0; j < n; j ++ )
                    {
                        nextValueLine[j] = iterationValue[j];
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
            //kNext = EXTRAPOLATION_MATRIX_SIZE;
            for(int s = 0; s < n; s ++)
            {
                resultedOutput[m + 1][s] = iterationValue[s];
                //printf("%lf %lf\n", accuracyMatrix[m][s], resultedOutput[m+1][s]);
                accuracyMatrix[m][s] -= iterationValue[s]; 

            }
            free(iterationValue);
            segmentBegin += h;
        }

        // printf("RESULT MATRIX\n");
        // for(int m = 0; m < k + 1; m ++)
        // {
        //     for(int s = 0; s < n; s ++)
        //     {
        //         printf("%lf ", resultedOutput[m][s]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        for(int m = 0; m < n; m ++)
        {
            for(int i = 0; i < k; i ++)
            {   
                accuracy += accuracyMatrix[i][m] * accuracyMatrix[i][m];
                
                //printf("SEGMENT %d ACCURACY: %lf\n", i + 1, accuracyMatrix[i][m]);

            }
            accuracy = sqrt(accuracy) / k;
            printf("|%.10lf\n", accuracy);
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

    free(nextValueLine);
    for(int i = 0; i < k; i ++)
    {
        free(accuracyMatrix[i]);
    }
    free(accuracyMatrix);
    // for(int i = 0; i < k + 1; k ++)
    // {
    //     free(resultedOutput[i]);
    // }
    // free(resultedOutput);
    return resultedOutput;
}


double* solveRunge(int n, double a, double b, int k, double* y0) 
{
    struct Runge Runge;
    Runge.a2 = 0.2;
    Runge.a3 = 0.6;
    Runge.a4 = 1;
    Runge.b21 = Runge.a2;

    Runge.b32 = (Runge.a3 * (Runge.a3 - Runge.a2)) / (2 * Runge.a2 * (1 - 2 * Runge.a2));
    Runge.b31 = Runge.a3 - Runge.b32;

    Runge.c2 = (2 * Runge.a3 - 1) / (12 * Runge.a2 * (Runge.a3 - Runge.a2) * (1 - Runge.a2));
    Runge.c3 = (2 * Runge.a2 - 1) / (12 * Runge.a3 * (Runge.a2 - Runge.a3) * (1 - Runge.a3));
    Runge.c4 = (6 * Runge.a2 * Runge.a3 - 4 * Runge.a2 - 4 * Runge.a3 + 3) / (12 * (1 - Runge.a2) * (1 - Runge.a3));
    Runge.c1 = 1 - Runge.c2 - Runge.c3 - Runge.c4;
    Runge.b42 = (Runge.c2 * (1 - Runge.a2) - Runge.c3 * Runge.b32) / Runge.c4;

    Runge.b43 = (Runge.c3 * (1 - Runge.a3)) / Runge.c4;
    Runge.b41 = 1 - Runge.b42 - Runge.b43;

    double h = (b - a) / k;

    double* f_result;
    double x_temp;

    double* y_temp, * y_new_temp, * y_super_new_temp;
    double* k1_line, * k2_line, * k3_line, * k4_line;
    double* result_line;

    double* output;

    y_temp = (double*)malloc(n * sizeof(double));
    y_new_temp = (double*)malloc(n * sizeof(double));
    y_super_new_temp = (double*)malloc(n * sizeof(double));
    k1_line = (double*)malloc(n * sizeof(double));//массив из к1 для каждой функции из системы
    k2_line = (double*)malloc(n * sizeof(double));
    k3_line = (double*)malloc(n * sizeof(double));
    k4_line = (double*)malloc(n * sizeof(double));

    // result_line = (double*)malloc((n + 1) * sizeof(double));
    result_line = (double*)malloc(n * sizeof(double));
    output = (double*)malloc((n + 1) * sizeof(double));

    for (int i = 0; i < n; i++) 
    {
        y_temp[i] = y0[i];
    }

    x_temp = a;

    f_result = (double*)malloc(k * sizeof(double));
    for (int i = 0; i < k + 1; i++) 
    {

        double x_1 = x_temp + Runge.a2 * h;
        double x_2 = x_temp + Runge.a3 * h;
        double x_3 = x_temp + Runge.a4 * h;
        //result_line[0] = x_temp;

        for (int j = 0; j < n; j++) 
        {
            f_result[j] = funk(j, x_temp, y_temp);
            k1_line[j] = f_result[j] * h;

            y_new_temp[j] = y_temp[j] + Runge.b21 * k1_line[j];
            //printf("%lf| %lf| ", y_new_temp[j], k1_line[j]);
        }

        for (int j = 0; j < n; j++) 
        {
            y_super_new_temp[j] = y_new_temp[j];
        }
        //printf("\n");
        for (int j = 0; j < n; j++) 
        {
            f_result[j] = funk(j, x_1, y_super_new_temp);
            k2_line[j] = f_result[j] * h;
            y_new_temp[j] = y_temp[j] + Runge.b31 * k1_line[j] + Runge.b32 * k2_line[j];
            //printf("%lf| %lf| ", y_new_temp[j], k2_line[j]);
        }
        for (int j = 0; j < n; j++) 
        {
            y_super_new_temp[j] = y_new_temp[j];
        }

        //printf("\n");
        for (int j = 0; j < n; j++) 
        {
            f_result[j] = funk(j, x_2, y_super_new_temp);
            k3_line[j] = f_result[j] * h;
            y_new_temp[j] = y_temp[j] + Runge.b41 * k1_line[j] + Runge.b42 * k2_line[j] + Runge.b43 * k3_line[j];
            //printf("%lf| %lf| ", y_new_temp[j], k3_line[j]);

        }

        for (int j = 0; j < n; j++) 
        {
            y_super_new_temp[j] = y_new_temp[j];
        }
        //printf("\n");
        for (int j = 0; j < n; j++) 
        {
            f_result[j] = funk(j, x_3, y_super_new_temp);
            k4_line[j] = f_result[j] * h;

            y_temp[j] += Runge.c1 * k1_line[j] + Runge.c2 * k2_line[j] + Runge.c3 * k3_line[j] + Runge.c4 * k4_line[j];
            result_line[j] = y_temp[j];
        }
        //printf("\n");
        x_temp = x_temp + h;
    }

    free(y_temp);
    free(y_new_temp);
    free(y_super_new_temp);
    free(k1_line);
    free(k2_line);
    free(k3_line);
    free(k4_line);
    free(output);
    free(f_result);
    return result_line;
}