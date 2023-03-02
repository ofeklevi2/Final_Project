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

void free_arr(double **arr, int len){ //frees 2d array 
    int i;
    for (i=0; i<len; i++){
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

int *find_Indexes_Of_Max_Element(double **L, int len){
    int i, j;
    double maximum = fabs(L[0][1]);
    int *ij = (int*)malloc(2 * sizeof(int));
    if (ij == NULL){
        return NULL;
    }
    ij[0] = 0, ij[1] = 1;
    for (i = 0; i < len; i++){
        for (j = i + 1; j < len; j++){
            if (fabs(L[i][j]) > maximum){
                maximum = L[i][j];
                ij[0] = i;
                ij[1] = j;
            }
        }
    }
    return ij;
}
double calc_t(int i, int j, double **A){
    double jj = A[j][j];
    double ii = A[i][i];
    double theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    int sign = (theta >= 0) ? 1 : -1; //if theta >= 0 then sign = 1, else sign = -1
    double t = sign / (fabs(theta) + sqrt((theta * theta) + 1));
    return t;
}
double calc_c(int i, int j, double **A){
    double t = calc_t(i, j, A);


    double c = 1 / sqrt((t * t) + 1); 
    return c;
}

double calc_s(int i, int j, double **A){
    double t = calc_t(i, j, A);
    double c = 1 / sqrt((t * t) + 1); 
    double s = t * c;
    return s;

}

double **build_Rotation_Matrix_P(double **A, int len){
    int k,i,j;
    double **P = (double**)malloc(len * sizeof(double*));
    if (P == NULL){
        return NULL;
    }
    int *ij = find_Indexes_Of_Max_Element(A, len);
    if (ij == NULL){
        return NULL;
    }
    i = ij[0];
    j = ij[1];
        for (k = 0; k < len; k++){
        P[k] = (double*)calloc(len, sizeof(double*));
        if (k == i){
            P[i][i] = calc_c(i, j, A);
            P[i][j] = calc_s(i, j, A);
        }
        else if(k == j){
            P[j][j] =  P[i][i];
            P[j][i] = -P[i][j];
        }
        else{
            P[k][k] = 1;
        }
       
    }
    free(ij);
    return P;

}

void get_A_Prime(int i, int j, double **A, int len, double c, double s){
    int row, col;
    double ii, jj, ij;
    for (row =0; row < len; row++){
        if (row != i && row != j){
            double ri = c * A[row][i] - s * A[row][j];
            double rj = c *A[row][j] + s * A[row][i];
            A[row][i] = ri, A[i][row] = ri;
            A[row][j] = rj, A[j][row] = rj;
        }
    }
    ii = (c * c) * A[i][i] + (s * s) * A[j][j] - 2 * s * c * A[i][j]; 
    jj = (c * c) * A[i][i] + (s * s) * A[j][j] + 2 * s * c * A[i][j]; 
    ij = 0;

    A[i][i] = ii;
    A[j][j] = jj;
    A[i][j] = ij, A[j][i] = ij;
}

double off(double **A, int len){
    int i, j;
    double sum = 0;
    for (i = 0; i < len; i++){
        for (j =0; j < len; j++){
            if (i != j){
                sum += pow(A[i][j], 2);
            } 
        }
    }
    return sum;
}


double **wam_c(double **dataPoints, int len){
    int i, j;
    double **W = (double**)malloc(len * sizeof(double*));
    if (W == NULL){
        return NULL;
    }
    for (i = 0; i < len; i++){
       W[i] = (double*)malloc(len * sizeof(double*));
    }

    for (i = 0; i < len; i++){
        for (j = i; j < len; j++){
            if (i == j){
                W[i][j] = 0;
            }
            else {
                W[i][j] = distance(dataPoints[i], dataPoints[j], len);
                W[j][i] = W[i][j];
             }
        }
    }

    return W;
}


double **ddg_c(double **dataPoints, int len){
    int i, j;
    double **w = wam_c(dataPoints, len);
    double **diagonal_Degree_Matrix = (double**)malloc(len * sizeof(double*));
    double **W = wam_c(dataPoints, len);
    if (W == NULL){
        return NULL;
    }
    double **D = (double**)malloc(len * sizeof(double*));
    if (D == NULL){
        return NULL;
    }

    for (int i =0; i < len; i++){
        D[i] = (double*)calloc(len, sizeof(double));
    }
    for (int i =0; i < len; i++){
        double sum = 0;
        for (int j =0; j < len; j++){
            sum += W[i][j];
        }
        D[i][i] = sum;
    }
    free_arr(W, len);
    return D;
}

double **gl_c(double **dataPoints, int len){
    int i, j;
    double **W, **D, **L;
    L = (double**)malloc(len * sizeof(double*));
    if (L == NULL){
        return NULL;
    }
    for (i = 0; i < len; i++){
        L[i] = (double*)malloc(len * sizeof(double));
    }
    W = wam_c(dataPoints, len);
    if (W == NULL){
        return NULL;
    }
    D = ddg_c(dataPoints, len);
    if (W == NULL){
        return NULL;
    }

    w = wam_c(dataPoints, len);
    dd_Matrix = ddg_c(dataPoints, len);
    for (i = 0; i < len; i++){
        for (j = 0; j < len; j++){
           L[i][j] = D[i][j] - W[i][j];
        }
    }
    free_arr(W, len);
    free_arr(D, len);
    return L;
}

void get_A_Converence(double **A, int len){
    int i, j;
    double c, s, sum1, sum2;
    int *ij;
    double eps = 1.0 * pow(10, -5);

    ij = find_Indexes_Of_Max_Element(A, len);
    i = ij[0], j =ij[1];
    c = calc_c(i, j, A);
    s = calc_s(i, j, A);
    sum1 = off(A, len);
    printf("A:\n");
    print_2D_Array(A, len);
    printf("\n");
    get_A_Prime(i, j, A, len, c, s);
    sum2 = off(A, len);
    printf("A':\n");
    print_2D_Array(A, len);
    printf("\n");
    printf("sum1 = %lf, sum2 = %lf, Diff = %lf ", sum1, sum2, sum1 - sum2);
    printf("\n\n\n");

    while (sum1 - sum2 > eps){
        sum1 = sum2;
        ij = find_Indexes_Of_Max_Element(A, len);
        i = ij[0], j =ij[1];
        c = calc_c(i, j, A);
        s = calc_s(i, j, A);
        get_A_Prime(i, j, A, len, c, s);
        sum2 = off(A, len);
        printf("A':\n");
        print_2D_Array(A, len);
        printf("\n");
        printf("sum1 = %lf, sum2 = %lf, Diff = %lf ", sum1, sum2, sum1 - sum2);
        printf("\n\n\n");
    }
}

//##################### Tests_Section #################################

void test_1(){
    int i, j, len;
    int *ij;
    double c, s, sum1, sum2;
    double **W, **D, **L, **P, **A_Prime;
    double **dataPoint = malloc(3 * sizeof(double*));
    double *x1 = malloc(3 * sizeof(double));
    double *x2 = malloc(3 * sizeof(double));
    double *x3 = malloc(3 * sizeof(double));

    len = 3;
    x1[0] = 1.0, x1[1] = 2.0, x1[2] = 3.0;
    x2[0] = 1.1, x2[1] = 2.1, x2[2] = 3.1;
    x3[0] = 1.2, x3[1] = 2.2, x3[2] = 3.2;

    dataPoint[0] = x1,   dataPoint[1] = x2,   dataPoint[2] = x3;

    W = wam(dataPoint,len);
    D = ddg(dataPoint, len);
    L = gl(dataPoint, len);
    P = build_Rotation_Matrix_P(L, len);

    // printf("W:\n");
    // print_2D_Array(W, len);
    // printf("\nD:\n");
    // print_2D_Array(D, len);
    // printf("\nL:\n");
    // print_2D_Array(L, len);
    // printf("\nP:\n");
    // print_2D_Array(P, len);
    // printf("\n");
    ij = find_Indexes_Of_Max_Element(L, len);
    i = ij[0], j =ij[1];
    c = calc_c(i, j, L);
    s = calc_s(i, j, L);
    // sum1 = off(L, len);
    // get_A_Prime(i, j, L, len, c, s);
    // print_2D_Array(L, 3);
    // printf("\n");
    // sum2 = off(L, len);
    // printf("sum1 = %lf, sum2 = %lf, eps = %lf ", sum1, sum2, sum1 - sum2);
    // printf("\n");

    get_A_Converence(L, len);
    // print_2D_Array(L, len);

    free(ij);
    free_arr(W, len);
    free_arr(D, len);
    free_arr(L, len);
    free_arr(P, len);
}

void test_2(){
    int i, j, len;
    int *ij;
    double c, s, sum1, sum2;
    double **W, **D, **L, **P, **A_Prime;
    double **dataPoint = malloc(4 * sizeof(double*));
    double *x1 = malloc(4 * sizeof(double));
    double *x2 = malloc(4 * sizeof(double));
    double *x3 = malloc(4 * sizeof(double));
    double *x4 = malloc(4 * sizeof(double));
    len = 4;
    x1[0] = 1.0,  x1[1] = 2.0, x1[2] = 3.0,  x1[3] = 4.0;
    x2[0] = 1.1, x2[1] = 2.1,  x2[2] = 3.1, x2[3] = 4.1;
    x3[0] = 1.2,  x3[1] = 2.2, x3[2] = 3.2,  x3[3] = 10.3;
    x4[0] = 1.3, x4[1] = 2.3,  x4[2] = 3.3, x4[3] = 4.3;

    dataPoint[0] = x1,   dataPoint[1] = x2,   dataPoint[2] = x3, dataPoint[3] = x4;

    W = wam(dataPoint,len);
    D = ddg(dataPoint, len);
    L = gl(dataPoint, len);
    P = build_Rotation_Matrix_P(L, len);

    // printf("W:\n");
    // print_2D_Array(W, len);
    // printf("\nD:\n");
    // print_2D_Array(D, len);
    // printf("\nL:\n");
    // print_2D_Array(L, len);
    // printf("\nP:\n");
    // print_2D_Array(P, len);
    // printf("\n");
    ij = find_Indexes_Of_Max_Element(L, len);
    i = ij[0], j =ij[1];
    c = calc_c(i, j, L);
    s = calc_s(i, j, L);
    // sum1 = off(L, len);
    // get_A_Prime(i, j, L, len, c, s);
    // print_2D_Array(L, 3);
    // printf("\n");
    // sum2 = off(L, len);
    // printf("sum1 = %lf, sum2 = %lf, eps = %lf ", sum1, sum2, sum1 - sum2);
    // printf("\n");

    get_A_Converence(L, len);
    // print_2D_Array(L, len);

    free(ij);
    free_arr(W, len);
    free_arr(D, len);
    free_arr(L, len);
    free_arr(P, len);
}

//##################### End_Of_Tests_Section #################################


int main(){
    //test_1();
    test_2();
    return 0;
}

