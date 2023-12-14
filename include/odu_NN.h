#pragma once

#define EXTRAPOLATION_MATRIX_SIZE 5
#define ELEMENTS_COUNT_PER_SEGMENT 25

void solveODE(int n, double a, double b, double e, int k, double* y0, double** result, int *k1);
double* solveRunge(int n, double a, double b,  int k, double* y0);
double* GaussElimination(int pos, int n, double **matrix);
double* RichardsonExtrapolation(int n, double segmentBegin, double segmentEnd, int k, double* y0);