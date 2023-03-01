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

double** wam_c(double **dataPoints, int len);

double** ddg_c(double **dataPoints);

double** gl_c(double **dataPoints);

void free_arr(double **arr, int len);

# endif