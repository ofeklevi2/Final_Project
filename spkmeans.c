# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
#include "spkmeans.h"

void delete_cords(struct cord *head){
  if (head != NULL){
    delete_cords(head->next);
    free(head);
  }
}

void delete_vectors(struct vector *head){
  if (head != NULL){
    delete_cords(head->cords);
    delete_vectors(head->next);
    free(head);
  }
}

double **allocate_Memory(int rows, int cols){
    int i;
    double **A = (double**)malloc(rows * sizeof(double*));
    if (A == NULL){
        return NULL;
    }
    for (i = 0; i < rows; i++){
        A[i] = (double*)calloc(cols, sizeof(double));
        if (A[i] == NULL){
        return NULL;
        }
    }
    return A;
}

void print_1D_Array(double *arr, int len){
    int i;
    for (i = 0; i < len; i++){
        if (arr[i] == 0) printf("0.0000");
        else printf("%.4f", arr[i]);
        if (i != len - 1) printf(",");
    }
}

void print_2D_Array(double **arr, int dim_1, int dim_2){
    int i, j;
    for (i = 0; i < dim_1; i++){
        for (j = 0; j < dim_2; j++){
            if (arr[i][j] == 0) printf("0.0000");
            else printf("%.4f", arr[i][j]);
            if (j != dim_2 - 1) printf(",");
        }
        printf("\n");
    }
}

double *vector_to_arr(vector *head, int dim){
    double *res;
    cord *curr_cord;
    int i;
    curr_cord = head->cords;
    res = malloc(dim * sizeof(double));
    if (res == NULL){
        return NULL;
    }
    for (i=0; i<dim; i++){
        res[i] = curr_cord->value;
        curr_cord = curr_cord->next;
    }
    return res;
}

double **linked_list_to_arr(vector *head, int dim_1, int dim_2){
    double **res;
    vector *curr_vec; 
    int i;
    curr_vec = head;
    res = malloc(dim_1*sizeof(double*));
    if (res == NULL){
        return NULL;
    }
    for (i=0; i<dim_1; i++){
        res[i] = vector_to_arr(curr_vec,dim_2);
        curr_vec = curr_vec->next;
    }
    return res;
}


int linked_list_len(struct vector *head_vec){
    int cnt;
    cnt = 0;
    for(; head_vec->next != NULL; head_vec = head_vec->next){
      ++cnt;
  }
  return cnt;
}

int vector_len(vector *vec){
    cord *curr_cord;
    int cnt = 0;
    curr_cord = vec->cords;
    while(curr_cord != NULL){
        cnt++;
        curr_cord = curr_cord->next;
    }
    return cnt;
}
void free_arr(double **arr, int rows){ //frees 2d array 
    int i;
    for (i=0; i < rows; i++){
        free(arr[i]);
    }
    free(arr);
}

double distance(double *xi, double *xj, int len, int vector_len){
    int i;
    double sum =0;
    for (i = 0; i < vector_len; i++){
        sum += pow(xi[i] - xj[i], 2);
    }
    // printf("sum = %lf \n", sum);
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
            if (fabs(L[i][j]) > fabs(maximum)){
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
        if (P[k] == NULL){
            return NULL;
        }
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
    jj = (s * s) * A[i][i] + (c * c) * A[j][j] + 2 * s * c * A[i][j]; 
    ij = 0.0;

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

double **matrix_Multiplication(double **A, double **B, int len){
    int i, j, k;
    double **res = allocate_Memory(len, len);
    if (res == NULL){
        return NULL;
    }
    for (i = 0; i < len; i++){
        for (j = 0; j < len; j++){
            for (k = 0; k < len; k++){
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

double **I(int len){
    int i, j;
    double **res = allocate_Memory(len, len);
    if (res == NULL){
        return NULL;
    }
    for (i = 0; i < len; i++){
        for (j = 0; j < len; j++){
            res[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    return res;
}


double **wam_c(double **dataPoints, int len, int vector_len){
    int i, j;
    double **W = (double**)malloc(len * sizeof(double*));
    if (W == NULL){
        return NULL;
    }
    for (i = 0; i < len; i++){
       W[i] = (double*)malloc(len * sizeof(double*));
        if (W[i] == NULL){
            return NULL;
        }   
    }

    for (i = 0; i < len; i++){
        for (j = i; j < len; j++){
            if (i == j){
                W[i][j] = 0;
            }
            else {
                W[i][j] = distance(dataPoints[i], dataPoints[j], len, vector_len);
                W[j][i] = W[i][j];
             }
        }
    }

    return W;
}


double **ddg_c(double **dataPoints, int len, int vector_len){
    int i, j;
    double **W = wam_c(dataPoints, len, vector_len);
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

double **gl_c(double **dataPoints, int len, int vector_len){
    int i, j;
    double **W, **D, **L;
    L = (double**)malloc(len * sizeof(double*));
    if (L == NULL){
        return NULL;
    }
    for (i = 0; i < len; i++){
        L[i] = (double*)malloc(len * sizeof(double));
        if (L[i] == NULL){
            return NULL;
        }    
    }
    W = wam_c(dataPoints, len, vector_len);
    if (W == NULL){
        return NULL;
    }
    D = ddg_c(dataPoints, len, vector_len);
    if (W == NULL){
        return NULL;
    }

    for (i = 0; i < len; i++){
        for (j = 0; j < len; j++){
           L[i][j] = D[i][j] - W[i][j];
        }
    }
    free_arr(W, len);
    free_arr(D, len);
    return L;
}

double **transpose(double **J, int rows, int cols){ //rows is the #rows of J, cols is ther #cols of J
    int i, j;
    double **J_Transpose = allocate_Memory(cols, rows);
    if (J_Transpose == NULL){
        return NULL;
    }
    for (i = 0; i < cols; i++){
        for (j = 0; j < rows; j++){
            J_Transpose[i][j] = J[j][i];
        }
    }
    return J_Transpose;
}

int comparator (const void *x1, const void *x2){ // x1 and x2 are eigenvalue struct (rows of J_Transpose)
    eigenvalue *a =  (eigenvalue*) x1;
    eigenvalue *b =  (eigenvalue*) x2;

    if (a->value < b->value) return -1;
    else if (a->value > b->value) return 1;
    else if (a->index < b->index) return -1;
    else if (a->index > b->index) return 1;
    else return 0;
}

double **sort_Rows(double **J_Transpose, int len){ // sort matrix by increasing order by its first entry of each row
                                                   // len is number of rows == number of eigenvalues
    int i, j;
    double val;
    eigenvalue *arr = (eigenvalue*)malloc(len * sizeof(eigenvalue));
    if (arr == NULL){
        return NULL;
    }
    double **res = allocate_Memory(len, len + 1);
    if (res == NULL){
        return NULL;
    }
    for (i = 0 ; i < len; i++){
        val = J_Transpose[i][0];
        eigenvalue a = {.index = i, .value = val, .row = J_Transpose[i]};
        arr[i] = a;
    }
    qsort(arr, len, sizeof(eigenvalue), comparator);
    for (i = 0; i < len; i++){
        for (j = 0; j < len + 1; j++){
            res[i][j] = arr[i].row[j];
        }
    }
    free(arr);

    return res;
}


double **jacobi_c(double **A, int len, int sort){

    int i, j, k, iter = 0;
    double c, s, sum1, sum2, **tmp, ***V_saver;
    int *ij;
    double **J = allocate_Memory(len + 1, len); // J for jacobi
    double **J_Transpose, **sorted_J_Transpose; 
    if (J == NULL){
        return NULL;
    }
    double **V = I(len);
    if (V == NULL){
        return NULL;
    }
    double eps = 1.0 * pow(10, -5);
    double **P = build_Rotation_Matrix_P(A, len);
    if (P == NULL){
        return NULL;
    }
    V_saver = (double***) malloc (102 * sizeof(double**));
    tmp = V;  
    V = matrix_Multiplication(V, P, len);
    if (V == NULL){
        return NULL;
    }  
    V_saver[0] = V;
    free_arr(tmp,len); 
    ij = find_Indexes_Of_Max_Element(A, len);
    if (ij == NULL){
        return NULL;
    }  
    i = ij[0], j =ij[1];
    c = calc_c(i, j, A);
    s = calc_s(i, j, A);
    sum1 = off(A, len);
    get_A_Prime(i, j, A, len, c, s);
    iter++;
    sum2 = off(A, len);
    free(ij);
    free_arr(P, len);      
    while (sum1 - sum2 > eps && iter < 100){
        sum1 = sum2;
        P = build_Rotation_Matrix_P(A, len);
        //tmp = V;
        V = matrix_Multiplication(V, P, len);
        V_saver[iter] = V;
        //free_arr(tmp,len);
        if (P == NULL){
            return NULL;
        } 

        ij = find_Indexes_Of_Max_Element(A, len);
        if (ij == NULL){
            return NULL;
        } 
        i = ij[0], j =ij[1];
        c = calc_c(i, j, A);
        s = calc_s(i, j, A);       
        get_A_Prime(i, j, A, len, c, s);
        iter++;
        sum2 = off(A, len);
        free(ij);
        free_arr(P, len); 

        
    }

    for(j = 0; j < len; j++){ //first row of J contains eigenvalues
        if (A[j][j] < 0 && A[j][j] > -0.0001){
            A[j][j] = 0;
            for (k = 0; k < len; k++){
                V[k][j] *= -1;
            }
        }
        J[0][j] = A[j][j];
    }
    
    for (i = 1; i < len + 1; i++){ //The other rows are the corresponding eigenvectors of the first rows (which exactly idencial to V's rows)
         for (j = 0; j < len; j++){
             J[i][j] = V[i - 1][j];
       }
     }


    //############################## Start spk() sort J code here #################################################
    if (sort == 1){
        J_Transpose = transpose(J, len + 1, len);  // (len + 1) : #rows of J, (len) - #cols of J
        if (J_Transpose == NULL){
            return NULL;
        }  
        sorted_J_Transpose = sort_Rows(J_Transpose, len); 
        
        //now transpose sorted_J_Transpose to get sorted_J
        for (i = 0; i < len + 1; i++){
            for(j = 0; j < len; j++){
                J[i][j] = sorted_J_Transpose[j][i];
            }
        }
        
        free_arr(J_Transpose,len);
        free_arr(sorted_J_Transpose,len); 
    }

        //############################## End spk() sort J code here #################################################
    for (i = 0; i < iter ; i++){
        free_arr(V_saver[i], len);
    }
    free(V_saver);
    return J;
}


//################ Testers #########################
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
    x1[0] = 1.0,       x1[1] = sqrt(2.0), x1[2] = 2.0;
    x2[0] = sqrt(2.0), x2[1] = 3.0,       x2[2] = sqrt(2.0);
    x3[0] = 2.0,       x3[1] = sqrt(2.0), x3[2] = 1.0;

    dataPoint[0] = x1,   dataPoint[1] = x2,   dataPoint[2] = x3;

    W = wam_c(dataPoint,len,len);
    D = ddg_c(dataPoint, len,len);
    L = gl_c(dataPoint, len,len);
    P = build_Rotation_Matrix_P(L, len);

    // printf("W:\n");
    // print_2D_Array(W, len, len);
    // printf("\nD:\n");
    // print_2D_Array(D, len, len);
    // printf("\nL:\n");
    // print_2D_Array(L, len, len);
    // printf("\nP:\n");
    // print_2D_Array(P, len, len);
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

    jacobi_c(L, len, 0);
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
    double **dataPoint = malloc(3 * sizeof(double*));
    double *x1 = malloc(3 * sizeof(double));
    double *x2 = malloc(3 * sizeof(double));
    double *x3 = malloc(3 * sizeof(double));

    len = 3;
    x1[0] = 1.0, x1[1] = 2.0, x1[2] = 3.0;
    x2[0] = 1.1, x2[1] = 2.1, x2[2] = 3.1;
    x3[0] = 1.2, x3[1] = 2.2, x3[2] = 3.2;

    dataPoint[0] = x1,   dataPoint[1] = x2,   dataPoint[2] = x3;

    W = wam_c(dataPoint,len,len);
    D = ddg_c(dataPoint, len,len);
    L = gl_c(dataPoint, len,len);
    P = build_Rotation_Matrix_P(L, len);

    // printf("W:\n");
    // print_2D_Array(W, len);
    // printf("\nD:\n");
    // print_2D_Array(D, len);
    // printf("\nL:\n");
    // print_2D_Array(L, len, len);
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

    jacobi_c(L, len, 0);
    // print_2D_Array(L, len);

    free(ij);
    free_arr(W, len);
    free_arr(D, len);
    free_arr(L, len);
    free_arr(P, len);
}

// ############### End of testers ##################


int main(int argc, char** argv){

    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    double n, **res;
    char c, *goal;
    int dim1,dim2;
    if (argc != 3){
        printf("An Error Has Occurred");
        return 1;
    }
    goal = argv[1];


    head_cord = calloc(1, sizeof(struct cord));
    if (head_cord == NULL){
        printf("An Error Has Occurred\n");
        return 1; 
    }
    
    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = calloc(1, sizeof(struct vector));
    if (head_vec == NULL){  
        printf("An Error Has Occurred\n");
        return 1; 
    }
    curr_vec = head_vec;
    curr_vec->next = NULL;

    FILE *input_data = fopen(argv[2], "r");
    if (input_data == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    while (fscanf(input_data,"%lf%c", &n, &c) == 2)
        {

            if (c == '\n')
            {
                curr_cord->value = n;
                curr_vec->cords = head_cord;
                curr_vec->next = calloc(1, sizeof(struct vector));
                if (curr_vec->next == NULL){
                printf("An Error Has Occurred\n");
                return 1; 
                }
                curr_vec = curr_vec->next;
                curr_vec->next = NULL;
                head_cord = malloc(sizeof(struct cord));
                if (head_cord == NULL){
                printf("An Error Has Occurred\n");
                return 1; 
                }
                curr_cord = head_cord;
                curr_cord->next = NULL;
                continue;
            }

            curr_cord->value = n;
            curr_cord->next = calloc(1, sizeof(struct cord));
            if (curr_cord->next == NULL){
            printf("An Error Has Occurred\n");
            return 1; 
            }
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
    }
    fclose(input_data);
    
    dim1 = linked_list_len(head_vec);
    dim2 = vector_len(head_vec);
    double** dataPoints = linked_list_to_arr(head_vec, dim1,dim2);
    if (strcmp(goal, "wam") == 0){
        res = wam_c(dataPoints,dim1,dim2);
        
    }
    else if (strcmp(goal, "ddg") == 0){
        res = ddg_c(dataPoints,dim1,dim2);
    }
    else if(strcmp(goal, "gl") == 0){
        res = gl_c(dataPoints,dim1,dim2);
    }
    else if(strcmp(goal, "jacobi") == 0){
        res = jacobi_c(dataPoints,dim1, 0);
    }

    else if(strcmp(goal, "ofek") == 0){
        res = jacobi_c(dataPoints,dim1, 1);
    }

    if (strcmp(goal, "jacobi") == 0){
        print_2D_Array(res,dim1 + 1,dim2);
        free_arr(res,dim1+1);
    }
    else if (strcmp(goal, "ofek") == 0){
        print_2D_Array(res,dim1 + 1,dim2);
        free_arr(res,dim1+1);
    }
    else{
        print_2D_Array(res,dim1,dim1);
        free_arr(res,dim1);   
    }

    delete_vectors(head_vec);
    free(head_cord);
    free_arr(dataPoints,dim1);

   return 0;
}


