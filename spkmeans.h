# ifndef SPKMEANS_H_
# define SPKMEANS_H_

struct cord
{
    double value;
    struct cord *next;
};
struct vector
{
    struct vector *next;
    struct cord *cords;
};

public **double wam(**double dataPoints)

# endif