# ifndef SPKMEANS_H_
# define SPKMEANS_H_

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

typedef struct eigenvalue{
    double value;
    int index;
    double *row;

}eigenvalue;

double** wam_c(double **dataPoints, int len, int vector_len);

double** ddg_c(double **dataPoints, int len, int vector_len);

double** gl_c(double **dataPoints, int len, int vector_len);

void free_arr(double **arr, int len);

int linked_list_len(struct vector *head_vec);

void delete_cords(struct cord *head);

void delete_vectors(struct vector *head);

# endif