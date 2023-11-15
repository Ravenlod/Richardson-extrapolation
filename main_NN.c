// main_NN.c
// #define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define EXTRAPOLATION_MATRIX_SIZE 4
#define ELEMENTS_COUNT_PER_SEGMENT 20

//#include "odu_NN.h"
double funk(int i, double x, double* y);
void solveODE(int n, double a, double b, double e, int k, double* y0, double** result);
double* solveRunge(int n, double a, double b,  int k, double* y0);
int recursiveSearch(int pos, int n, double matrix[][6]);
void systemOfLinearFunctions(int n, double* leftLine, int* hLine);



struct Runge
{
    int k;
    double a, b, a2, a3, a4, b21, b31, b32, b41, b42, b43, c1, c2, c3, c4;
};

int main() 
{



    int n = 2;
    double a, b;
    double e;
    int k = 10;
    double* y0;
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
    double y0Line[] = { 1,-2 };
    solveODE(n, 0, 20, 0.01, k, y0Line, result);

    /*double** matrix;
    matrix = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        matrix[i] = (double*)malloc(n * sizeof(double));
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = exp(ratio * (j + 1) * hLine[i]);
        }
    }*/
    

    return 0;
}
double funk(int i, double x, double* y) 
{
    double result;

    switch (i) 
    {
        case 0:
            result = sin(x) * y[0] + cos(x) * y[1];  // Пример: y'[0] = x * y[0]
            break;
        case 1:
            result = sin(x) * y[1] + cos(x) * y[0];  // Пример: y'[1] = sin(x) * y[1]
            break;
        default:
            return NAN;
    }

    return result;
}

void systemOfLinearFunctions(int n, double* leftLine, int* hLine) 
{
    int ratio=2;
    double** matrix;
    matrix = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) 
    {
        matrix[i] = (double*)malloc(n * sizeof(double));
    }
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            matrix[i][j] = exp(ratio * (j + 1) * hLine[i]);
        }
    }
    
}

//Рекурсивный проход по матрице. Неявно возвращает матрицу в треугольном виде.
int recursiveSearch(int pos, int n, double matrix[][6]) 
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
                //printf("%lf ", copyMatrix[m][j]);
            }
            //printf("\n");
        }
        //printf("\n");
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
                //printf("%lf ", matrix[m][j]);
            }
            //printf("\n ");
        }
        //printf("\n");

        //Костыль для ограничения расчётов косинуса и синуса
        if (i + 2 != n - pos) 
        {
            cos = matrix[pos][pos] / (sqrt(matrix[pos][pos] * matrix[pos][pos] + matrix[pos + 2 + i][pos] * matrix[pos + 2 + i][pos]));
            sin = matrix[pos + 2 + i][pos] * cos / matrix[pos][pos];
        }

    }
    free(copyMatrix);
    pos += 1;

    //выход из функции
    if (pos == n - 1) 
    {
        return 1;
    }

    //вход в рекурсию со сдвигом pos + 1
    else
    {
        recursiveSearch(pos, n, matrix);
        return 1;
    }
}

void solveODE(int n, double a, double b, double e, int k, double* y0, double** result) 
{
    //Allocate memory to store the result of the Runge-Kutta method
    double *myResult, **leftLine;
    leftLine = (double**)malloc(EXTRAPOLATION_MATRIX_SIZE * sizeof(double*));
    double** transportedLeftLine = (double**)malloc(n * sizeof(double*));
    /* for(int i=0;i<EXTRAPOLATION_MATRIX_SIZE;i++){
        leftLine[i]=(double*)malloc((n + 1) * sizeof(double));
    } */
    //Call the Runge-Kutta method to solve the ODE

    double h = (b-a)/k;
    double x0 = a;
    double x1 = x0+h;
    myResult = solveRunge(n, a, b, k, y0);
    for(int i=0, k_change=ELEMENTS_COUNT_PER_SEGMENT; i<EXTRAPOLATION_MATRIX_SIZE; i++, k_change*=2){
        leftLine[i] = solveRunge(n, x0, x1, k_change, y0);
    }
    printf("#################START\n");

    for(int i=0, k_change=2;i<EXTRAPOLATION_MATRIX_SIZE;i++, k_change*=2){
        for(int j=0;j<n;j++){
        printf("%lf ", leftLine[i][j]);
            
        }
        printf("\n");

    }
    printf("#################MiDDLE\n");

    

    for (int i = 0; i < n; i++) {
        transportedLeftLine[i] = (double*)malloc(EXTRAPOLATION_MATRIX_SIZE * sizeof(double));
        for (int j = 0; j < EXTRAPOLATION_MATRIX_SIZE; j++) {
            transportedLeftLine[i][j] = leftLine[j][i];
        }
    }

    for(int i=0, k_change=2;i<n;i++, k_change*=2){
        for(int j=0;j<EXTRAPOLATION_MATRIX_SIZE;j++){
        printf("%lf ", transportedLeftLine[i][j]);    
        }
        printf("\n");

    }
    printf("#################STOP\n");


    double matrix[5][6] = { {4,	3,	-2,	5,	-7,73}, {-3,2,	4,	-5,	2,-40}, {5,	2,	5,	-3,	6,-77},
                            {-2, 9,	-7,	3,	2,66}, {-6,	2,	4,	-1,	8,-54}};

    recursiveSearch(0,n,matrix);
    printf("*********************************\n");
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n + 1; j++) 
        {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }

    double *roots = (double*)malloc(sizeof(double) * n);
    double numerator = matrix[n - 1][n];
    for (int i = n; i > 0; i--)
    {
        roots[i - 1] = numerator/matrix[i - 1][i - 1];
        numerator = matrix[i - 2][n];
        for(int j = n; j > i - 1; j --)
        {
            numerator -= matrix[i - 2][j - 1] * roots[j - 1]; 
        }
    }

    printf("*********************************\n");
    for (int i = 0; i < n; i++) 
    {
        printf("%lf ", roots[i]);
    }
    

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

    return result_line;
}