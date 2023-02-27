# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>


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

void delete_vectors_arr(struct vector **arr, int size){
  int i;
  for (i=0; i<size; i++){
    delete_vectors(arr[i]);
  }
  free(arr);
}

int vec_list_len(struct vector *head_vec){
    int cnt;
    cnt = 0;
    for(; head_vec->next != NULL; head_vec = head_vec->next){
      ++cnt;
}
return cnt;
}

struct vector* compute_centroid_by_cluster(struct vector *head_vec){
  int cluster_len;
  struct vector *curr_vec, *centroid;
  struct cord *curr_cent_cord, *curr_vec_cord, *head_cord; 
  centroid = calloc(1, sizeof(struct vector));
  if (centroid == NULL){
    return NULL; 
  }
  cluster_len = vec_list_len(head_vec);
  head_cord = calloc(1, sizeof(struct cord));
  if (head_cord == NULL){
    return NULL; 
  }
  curr_vec = head_vec;

  while(curr_vec->next != NULL){
    curr_cent_cord = head_cord;
    curr_vec_cord = curr_vec->cords;
    while (curr_vec_cord != NULL){
      curr_cent_cord->value = (curr_cent_cord->value) + ((curr_vec_cord->value)/cluster_len);
      curr_vec_cord = curr_vec_cord -> next;
      if ((curr_vec_cord != NULL) & (curr_cent_cord->next == NULL)){
        curr_cent_cord->next = calloc(1, sizeof(struct cord));
        if (curr_cent_cord->next == NULL){
          return NULL; 
        }
        }
      
      curr_cent_cord = curr_cent_cord->next;
    }
    curr_vec = curr_vec->next;
  }
  centroid->cords = head_cord;
  return centroid;
}

double d(struct vector *v1, struct vector *v2){
  struct cord *curr_v1_cord, *curr_v2_cord;
  double res;
  res = 0; 
  curr_v1_cord = v1->cords;
  curr_v2_cord = v2->cords;
  while (curr_v1_cord != NULL){
    res = res + pow(((curr_v1_cord->value) - (curr_v2_cord->value)),2);
    curr_v1_cord = curr_v1_cord->next;
    curr_v2_cord = curr_v2_cord->next;
  }
  res = pow(res,0.5);
  return res;
}

int findClosest(struct vector *vec, struct vector **centroids, int K){
  double minDist, currDist;
  int minNum, i; 
  minDist = d(vec,centroids[0]);
  minNum=0;
  for (i = 0 ; i < K; i++){
    currDist = d(vec, centroids[i]);
    if (currDist < minDist){
      minDist = currDist; 
      minNum = i;
    }
  }
  return minNum; 
}

double compute_max_delta(struct vector** updated_centroids,struct vector** centroids, int K){
  double currDelta, maxDelta;
  int i;
  currDelta = 0;
  maxDelta = 0; 
  for (i=0; i<K; i++){
    currDelta = d(centroids[i], updated_centroids[i]);
    if (currDelta > maxDelta){
      maxDelta = currDelta;
    }
  }
  return maxDelta;
}

int gather_up(struct vector *head_vec, struct vector** clusters, int K){
  struct vector *curr,*prev, *first;
  int i;
  first = clusters[0];
  curr = first;
  for (i=1; i<K; i++){
    while(curr->next != NULL){
      curr = curr->next;
    }
    *curr = *clusters[i];
  }
  
  if(head_vec != clusters[0]){
    prev = first;
    curr = prev->next;
    while (curr->next != head_vec->next){
      prev = curr;
      curr = curr->next;
    }
    prev->next = curr->next;
    head_vec->next = first;
  }
  return 0; 
}

struct vector** compute_new_centroids(int K, struct vector *head_vec, struct vector **centroids){
    struct vector **clusters, *curr_vec, *next_vec, *tmp, **updated_centroids;
    int closest,i, check_no_error;
    curr_vec = head_vec;
    clusters = calloc(K, sizeof(struct vector*));
    if (clusters == NULL){
      return NULL; 
    }
    for (i = 0; i<K; i++){
      clusters[i] = calloc(1, sizeof(struct vector));
      if (clusters[i] == NULL){
        return NULL; 
      }
    }
    while(curr_vec -> next != NULL){
      closest = findClosest(curr_vec,centroids,K); 
      next_vec = curr_vec->next; 
      tmp = clusters[closest]; 
      clusters[closest] = curr_vec; 
      clusters[closest]->next = tmp;
      curr_vec = next_vec; 
    }
    free(curr_vec);
    updated_centroids = calloc(K, sizeof(struct vector*));
    if (updated_centroids == NULL){
      return NULL; 
    }
    for (i=0; i<K; i++){
      updated_centroids[i] = compute_centroid_by_cluster(clusters[i]);
      if (updated_centroids[i] == NULL){
        return NULL;
      }
    }
    check_no_error = gather_up(head_vec, clusters, K);
    if (check_no_error == 1){
      return NULL;
    }
    for (i=1; i<K; i++){
        free(clusters[i]);
    }
    free(clusters);
    return updated_centroids;
}
double** centroids_to_arr(struct vector **centroids, int K, int vec_len){
    double **res;
    int i,j;
    struct cord *curr_cord;
    res = malloc(K*sizeof(sizeof(double*)));
    if (res == NULL){
        return NULL;
    }
    for (i=0; i<K; i++){
        res[i] = malloc(vec_len * sizeof(double));
        if (res[i] == NULL){
            return NULL;
        }
        curr_cord = centroids[i]->cords;
        for(j=0; j<vec_len; j++){
            res[i][j] = curr_cord->value;
            curr_cord = curr_cord->next;
        }
    }
    return res;
}

double** kmeans(int K, int iter,int eps, struct vector *head_vec,int vec_len){
    int i, iteration_number;
    double delta, **res;
    struct vector **tmp, *curr_vec;
    struct vector **centroids, **updated_centroids; 
    struct cord *curr_cord, *curr_cent_cord;

    centroids = calloc(K, sizeof(struct vector*));
    if (centroids == NULL){
      return NULL; 
    }
    curr_vec = head_vec;
    for (i=0; i<K; i++){
      centroids[i] = calloc(1,sizeof(struct vector));
      if (centroids[i] == NULL){
        return NULL; 
      }
      centroids[i]->cords = calloc(1,sizeof(struct cord));
      if (centroids[i]->cords == NULL){
        return NULL; 
      }
      curr_cord = curr_vec->cords;
      curr_cent_cord = centroids[i]->cords;
      while(curr_cord != NULL){
        curr_cent_cord->value = curr_cord->value;
        if (curr_cord->next != NULL){
          curr_cent_cord->next = calloc(1,sizeof(struct cord));
          if (curr_cent_cord->next == NULL){
            return NULL; 
          }
        }
        curr_cord = curr_cord->next;
        curr_cent_cord = curr_cent_cord->next;
      }
      curr_vec = curr_vec->next;
    }


    iteration_number = 0;
    delta = eps + 1;
      while ((delta >= eps) & (iteration_number < iter)){ 
        updated_centroids = compute_new_centroids(K, head_vec, centroids);
        if (updated_centroids == NULL){
          return NULL;
        }
        delta = compute_max_delta(updated_centroids,centroids, K);
        tmp = centroids;
        centroids = updated_centroids;
        delete_vectors_arr(tmp, K);
        iteration_number++;
      }
    res = centroids_to_arr(centroids,K, vec_len);
    delete_vectors_arr(centroids,K);
    delete_vectors(head_vec);
    return res;
}


struct vector* initialize_head_vec (double **data_arr, int dim_1, int dim_2){
    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    int i,j; 
    
    head_vec = calloc(1, sizeof(struct vector));
    if (head_vec == NULL){
        return NULL; 
    }
    curr_vec = head_vec;

    for (i = 0; i < dim_1; i++){
      head_cord = calloc(1, sizeof(struct cord));
      if (head_cord == NULL){
        return NULL; 
      }
      curr_cord = head_cord;
      for (j=0; j < dim_2; j++){
        curr_cord->value = data_arr[i][j];
        if (j != dim_2 - 1){
          curr_cord->next = malloc(sizeof(struct cord));
          if (curr_cord == NULL){
            return NULL; 
          }
        }
        else{
          curr_cord->next = NULL;
        }
         curr_cord = curr_cord->next;
      }

      curr_vec->cords = head_cord;
      curr_vec ->next = calloc(1, sizeof(struct vector));
      if (curr_vec->next == NULL){
        return NULL; 
      }

      curr_vec = curr_vec->next;
    }
    return head_vec;
}

void fill_arr(double **arr, PyObject *list, int dim_1, int dim_2){
    int i,j;
    PyObject *item;
    
    for (i=0; i<dim_1; i++){
        item = PyList_GetItem(list,i);
        for (j=0; j<dim_2; j++){
            arr[i][j] = PyFloat_AsDouble(PyList_GetItem(item,j));
    }
  } 
}

static PyObject* fit(PyObject *self, PyObject *args){
    int K, iter, data_length, vec_len, i, j;
    double eps, **final_centroids, **data_arr;
    PyObject *data_list, *res, *python_float;
    struct vector *head_vec;
    if(!PyArg_ParseTuple(args, "iidOii", &K, &iter, &eps, &data_list,&data_length, &vec_len)) {
        return NULL; 
    }
    data_arr = malloc(data_length * sizeof(double*));
    if (data_arr == NULL){
      return NULL;
    }
    for (i=0; i < data_length; i++){
      data_arr[i] = malloc(vec_len * sizeof(double));
      if (data_arr[i] == NULL){
        return NULL;
    }
    }
    fill_arr(data_arr, data_list, data_length, vec_len);

    head_vec = initialize_head_vec (data_arr, data_length, vec_len);
    if (head_vec == NULL){
        return NULL;
    }

    final_centroids = kmeans(K, iter, eps, head_vec, vec_len);
    if (final_centroids == NULL){
      return NULL;
    }
    
    for (i=0; i < data_length; i++){
      free(data_arr[i]);
    }
    free(data_arr);

    res = PyList_New(K);
    for(i = 0; i<K; i++){
        PyList_SetItem(res, i, PyList_New(vec_len));
        for (j=0; j < vec_len; j++){
            python_float = Py_BuildValue("f",final_centroids[i][j]);
            PyList_SetItem(PyList_GetItem(res,i), j, python_float);
        }
    }

    for (i=0; i<K; i++){
      free(final_centroids[i]);
    }
    free(final_centroids);
    return res;
}

static PyMethodDef fitMethod[] = {
    {"fit",                   
      (PyCFunction) fit, 
      METH_VARARGS,          
      PyDoc_STR("The fit method expects 6 arguments : K (number of clusters), iter (max number of iterations), eps (the delta for converence), data_list (list of data points), data_length (the number of data points in data_list), vec_len (number of cordinates of each data point)")}, 
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", 
    NULL, 
    -1, 
    fitMethod
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}





    

