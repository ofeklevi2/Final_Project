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

void print_2D_Array(double **arr, int dim_1, int dim_2){
    int i, j;
    for (i = 0; i < dim_1; i++){
        for (j = 0; j < dim_2; j++){
            printf("%lf ", arr[i][j]);
        }
        printf("\n");
    }
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

double **wam_c(double **dataPoints, int len){
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


double **ddg_c(double **dataPoints, int len){
    int i, j;
    double **w = wam_c(dataPoints, len);
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

double **gl_c(double **dataPoints, int len){
    int i, j;
    double **w, **dd_Matrix, **gl_Matrix;
    gl_Matrix = (double**)malloc(len * sizeof(double*));
    for (i = 0; i < len; i++){
        gl_Matrix[i] = (double*)malloc(len * sizeof(double));
    }
    w = wam_c(dataPoints, len);
    dd_Matrix = ddg_c(dataPoints, len);
    for (i = 0; i < len; i++){
        for (j = 0; j < len; j++){
            gl_Matrix[i][j] = dd_Matrix[i][j] - w[i][j];
        }
    }
    return gl_Matrix;
}

void free_arr(double **arr, int len){ //frees 2d array 
    int i;
    for (i=0; i<len; i++){
        free(arr[i]);
    }
    free(arr);
}

void main(){
    double *x1 = malloc(3 * sizeof(double));
    double *x2 = malloc(3 * sizeof(double));
    x1[0] = 1.0, x1[1] = 2., x1[2] = 3.;
    x2[0] = 2., x2[1] = 3., x2[2] = 4.;
    double res = distance(x1, x2, 3);
    printf("%lf\n", res);
}
