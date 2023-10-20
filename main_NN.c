// main_NN.c
#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#include "odu_NN.h"
double funk(int i, double x, double* y);
void solveODE(int n, double a, double b, double e, int k, double* y0, double** result);


struct Runge
{
    int k;
    double a, b, a2, a3, a4, b21, b31, b32, b41, b42, b43, c1, c2, c3, c4;
};

int main() {



    int n = 2;
    double a, b;
    double e;
    int k = 50;
    double* y0;
    double** result;

    result = (double**)malloc((k + 1) * sizeof(double*));

    for (int i = 0; i <= k; i++) {
        result[i] = (double*)malloc((n + 1) * sizeof(double));
    }



    // Чтение входных данных из файла
    FILE* input = fopen("input.txt", "r");
    if (input == NULL) {
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
    double temp[] = { 1,-2 };
    solveODE(2, 0, 20, 0.01, k, temp, result);
    for (int i = 0; i <= k; i++) {
        for (int j = 0; j < n + 1; j++) {
            printf("%lf ", result[i][j]);
        }
        printf("\n");
    }


    return 0;
}
double funk(int i, double x, double* y) {
    double result;

    switch (i) {
    case 0:
        result = sin(x) * y[0] + cos(x) * y[1];  // Пример: y'[0] = x * y[0]
        break;
    case 1:
        result = sin(x) * y[1] + cos(x) * y[0];  // Пример: y'[1] = sin(x) * y[1]
        break;
    default:
        return NAN;
        // Добавьте обработку других функций, если есть
    }

    return result;
}

void solveODE(int n, double a, double b, double e, int k, double* y0, double** result) {}
double solveRunge(int n, double a, double b, double e, int k, double* y0, double** result) {
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

    result_line = (double*)malloc((n + 1) * sizeof(double));
    output = (double*)malloc((n + 1) * sizeof(double));

    for (int i = 0; i < n; i++) {
        y_temp[i] = y0[i];
    }

    x_temp = a;

    f_result = (double*)malloc(k * sizeof(double));
    for (int i = 0; i < k + 1; i++) {

        double x_1 = x_temp + Runge.a2 * h;
        double x_2 = x_temp + Runge.a3 * h;
        double x_3 = x_temp + Runge.a4 * h;
        result_line[0] = x_temp;

        for (int j = 0; j < n; j++) {
            f_result[j] = funk(j, x_temp, y_temp);
            k1_line[j] = f_result[j] * h;

            y_new_temp[j] = y_temp[j] + Runge.b21 * k1_line[j];
            //printf("%lf| %lf| ", y_new_temp[j], k1_line[j]);
        }

        for (int j = 0; j < n; j++) {
            y_super_new_temp[j] = y_new_temp[j];
        }
        //printf("\n");
        for (int j = 0; j < n; j++) {
            f_result[j] = funk(j, x_1, y_super_new_temp);
            k2_line[j] = f_result[j] * h;
            y_new_temp[j] = y_temp[j] + Runge.b31 * k1_line[j] + Runge.b32 * k2_line[j];
            //printf("%lf| %lf| ", y_new_temp[j], k2_line[j]);
        }
        for (int j = 0; j < n; j++) {
            y_super_new_temp[j] = y_new_temp[j];
        }

        //printf("\n");
        for (int j = 0; j < n; j++) {
            f_result[j] = funk(j, x_2, y_super_new_temp);
            k3_line[j] = f_result[j] * h;
            y_new_temp[j] = y_temp[j] + Runge.b41 * k1_line[j] + Runge.b42 * k2_line[j] + Runge.b43 * k3_line[j];
            //printf("%lf| %lf| ", y_new_temp[j], k3_line[j]);

        }

        for (int j = 0; j < n; j++) {
            y_super_new_temp[j] = y_new_temp[j];
        }
        //printf("\n");
        for (int j = 0; j < n; j++) {
            f_result[j] = funk(j, x_3, y_super_new_temp);
            k4_line[j] = f_result[j] * h;

            y_temp[j] += Runge.c1 * k1_line[j] + Runge.c2 * k2_line[j] + Runge.c3 * k3_line[j] + Runge.c4 * k4_line[j];
            result_line[j + 1] = y_temp[j];
            //printf("%lf|", y_temp[j]);

            //printf("%lf| %lf| ", y_temp[j], k4_line[j]);
        }
        //printf("\n");

        for (int j = 0; j < n + 1; j++) {
            result[i][j] = result_line[j];
        }
        x_temp = x_temp + h;
    }

    for (int i = 0; i < n + 1; i++) {
        output[i];
    }

}