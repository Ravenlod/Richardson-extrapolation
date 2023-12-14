// odu.h
#ifndef ODU_H
#define ODU_H

void solveODE(int n, double a, double b, double e, int k, double* y0, double** result, int *k1);
double* solveRunge(int n, double a, double b, int k, double* y0);
double* RichardsonExtrapolation(int n, double segmentBegin, double segmentEnd, int k, double* y0);
double* GaussElimination(int pos, int n, double **matrix);
void graphViewer(int n, double a, double b, int k, double* y0, double** result);


#endif  // ODU_H
