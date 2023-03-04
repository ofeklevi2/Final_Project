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

int main(int argc, char** argv){
    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    double n, **res;
    char c, *goal;
    int dim1,dim2;
    if (argc != 3){
        printf("An Error Has Occurred");
        return NULL;
    }
    goal = argv[1];


    head_cord = calloc(1, sizeof(struct cord));
    if (head_cord == NULL){
        printf("An Error Has Occurred\n");
        return NULL; 
    }
    
    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = calloc(1, sizeof(struct vector));
    if (head_vec == NULL){  
        printf("An Error Has Occurred\n");
        return NULL; 
    }
    curr_vec = head_vec;
    curr_vec->next = NULL;

    FILE *input_data = fopen(argv[2], "r");
    if (input_data == NULL){
        printf("An Error Has Occurred\n");
        return NULL;
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
                return NULL; 
                }
                curr_vec = curr_vec->next;
                curr_vec->next = NULL;
                head_cord = malloc(sizeof(struct cord));
                if (head_cord == NULL){
                printf("An Error Has Occurred\n");
                return NULL; 
                }
                curr_cord = head_cord;
                curr_cord->next = NULL;
                continue;
            }

            curr_cord->value = n;
            curr_cord->next = calloc(1, sizeof(struct cord));
            if (curr_cord->next == NULL){
            printf("An Error Has Occurred\n");
            return NULL; 
            }
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
    }
    dim1 = linked_list_len(head_vec);
    dim2 = vector_len(head_vec);
    double** dataPoints = linked_list_to_arr(head_vec, dim1,dim2);
    if (strcmp(goal, "wam") == 0){
        res = wam_c(dataPoints,dim1);
    }
    else if (strcmp(goal, "ddg") == 0){
        res = ddg_c(dataPoints,dim1);
    }
    else if(strcmp(goal, "gl") == 0){
        res = gl_c(dataPoints,dim1);
    }
    print_2D_Array(res,dim1,dim2);
    delete_vectors(head_vec);
}
