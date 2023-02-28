# include <stdio.h>
# include <stdlib.h>
# include <math.h>
#include "spkmeans.h"

void print_1D_Array(double *arr, int len){
    int i;
    for (i = 0; i < len; i++){
        printf("%lf ", arr[i]);
    }
}

void print_2D_Array(double **arr, int len){
    int i, j;
    for (i = 0; i < len; i++){
        for (j = 0; j < len; j++){
            printf("%lf ", arr[i][j]);
        }
        printf("\n");
    }
}

void free_arr(double **arr, int len){
    int i, j;
    for (i = 0; i < len; i++){
        free(arr[i]);
    }
    free(arr);
}

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

double **ddg(double **dataPoints, int len){
    int i, j;
    double **w = wam(dataPoints, len);
    double **diagonal_Degree_Matrix = (double**)malloc(len * sizeof(double*));
    for (int i =0; i < len; i++){
        diagonal_Degree_Matrix[i] = (double*)calloc(len, sizeof(double));
    }
    for (int i =0; i < len; i++){
        double sum = 0;
        for (int j =0; j < len; j++){
            sum += w[i][j];
        }

        diagonal_Degree_Matrix[i][i] = sum;
    }
    return diagonal_Degree_Matrix;
}

double **gl(double **dataPoints, int len){
    int i, j;
    double **w, **dd_Matrix, **gl_Matrix;
    gl_Matrix = (double**)malloc(len * sizeof(double*));
    for (i = 0; i < len; i++){
        gl_Matrix[i] = (double*)malloc(len * sizeof(double));
    }
    w = wam(dataPoints, len);
    dd_Matrix = ddg(dataPoints, len);
    for (i = 0; i < len; i++){
        for (j = 0; j < len; j++){
            gl_Matrix[i][j] = dd_Matrix[i][j] - w[i][j];
        }
    }
    return gl_Matrix;
}

int main(){
    return 0;
}
