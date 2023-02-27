# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#include "spkmeans.h"



double distance(double *xi, double *xj, int len){
    int i;
    double sum =0;
    for (i = 0; i < len; i++){
        sum += pow(xi[i] - xj[i], 2);
    }
    sum = - sum / 2;
    return exp(sum);
}

double **wam(double **dataPoints, int len){
    int i, j;
    double **w = (double**)malloc(len * sizeof(double*));
    for (i = 0; i < len; i++){
        w[i] = (double*)malloc(len * sizeof(double*));
    }

    for (i = 0; i < len; i++){
        for (j = i; j < len; j++){
            if (i == j){
                w[i][j] = 0;
            }
            else {
                w[i][j] = distance(dataPoints[i], dataPoints[j], len);
                w[j][i] = w[i][j];
             }
        }
    }

    return w;
}

void main(){
    double *x1 = malloc(3 * sizeof(double));
    double *x2 = malloc(3 * sizeof(double));
    x1[0] = 1.0, x1[1] = 2., x1[2] = 3.;
    x2[0] = 2., x2[1] = 3., x2[2] = 4.;
    double res = distance(x1, x2, 3);
    printf("%lf\n", res);
}