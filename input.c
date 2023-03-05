# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

typedef struct cord
{
    double value;
    struct cord *next;
}cord;
typedef struct vector
{
    struct vector *next;
    struct cord *cords;
}vector;

double *vector_to_arr(vector *head, int dim){
    double *res;
    cord *curr_cord;
    int i;
    curr_cord = head->cords;
    res = malloc(dim * sizeof(double));
    for (i=0; i<dim; i++){
        res[i] = curr_cord->value;
        curr_cord = curr_cord->next;
    }
    return res;
}

double** centroids_to_arr(struct vector **centroids, int dim_1, int dim_2){
    double **res;
    int i,j;
    struct cord *curr_cord;
    res = malloc(dim_1*sizeof(double*));
    if (res == NULL){
        return NULL;
    }
    for (i=0; i<dim_1; i++){
        res[i] = malloc(dim_2 * sizeof(double));
        if (res[i] == NULL){
            return NULL;
        }
        curr_cord = centroids[i]->cords;
        for(j=0; j<dim_2; j++){
            res[i][j] = curr_cord->value;
            curr_cord = curr_cord->next;
        }
    }
    return res;
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

int vec_list_len(struct vector *head_vec){
    int cnt;
    cnt = 0;
    for(; head_vec->next != NULL; head_vec = head_vec->next){
      ++cnt;
    }   
    return cnt;
}

  
int main(int argc, char **argv){

    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    double n, *res;
    char c;


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


    while (scanf("%lf%c", &n, &c) == 2)
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

    res = vector_to_arr(head_vec, 3);
    print_1D_Array(res,3);
    return 1;
}