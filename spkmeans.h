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

double** wam(double **dataPoints);

# endif