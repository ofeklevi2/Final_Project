# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "spkmeans.h"


double** py_to_c_arr(PyObject *list, int dim_1, int dim_2){//converts python list to c array 
    int i,j;
    double** data_arr;
    PyObject *item;
    data_arr = malloc(dim_1 * sizeof(double*));
    if (data_arr == NULL){
      return NULL;
    }
    for (i=0; i < dim_1; i++){
      data_arr[i] = malloc(dim_2 * sizeof(double));
      if (data_arr[i] == NULL){
        return NULL;
    }
    }
    for (i=0; i<dim_1; i++){
        item = PyList_GetItem(list,i);
        for (j=0; j<dim_2; j++){
            data_arr[i][j] = PyFloat_AsDouble(PyList_GetItem(item,j));
    }
  }
  return data_arr; 
}
PyObject* c_to_py_list(double **data_arr, int dim_1, int dim_2){//converts c array to python list 
    PyObject *res, *python_float;
    int i,j;
    res = PyList_New(dim_1);
    for(i=0; i < dim_1; i++){
        PyList_SetItem(res, i, PyList_New(dim_2));
        for (j=0; j < dim_2; j++){
            python_float = Py_BuildValue("f",data_arr[i][j]);
            PyList_SetItem(PyList_GetItem(res,i), j, python_float);
        }
    }
    return res;
}

void delete_vectors_arr(struct vector **arr, int size){
  int i;
  for (i=0; i<size; i++){
    delete_vectors(arr[i]);
  }
  free(arr);
}

struct vector* compute_centroid_by_cluster(struct vector *head_vec){
  int cluster_len;
  struct vector *curr_vec, *centroid;
  struct cord *curr_cent_cord, *curr_vec_cord, *head_cord; 
  centroid = calloc(1, sizeof(struct vector));
  if (centroid == NULL){
    return NULL; 
  }
  cluster_len = linked_list_len(head_vec);
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

double** kmeans(int K, struct vector *head_vec,int vec_len){
    int i, iteration_number;
    double  **res;
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
      while (iteration_number < 300){ 
        updated_centroids = compute_new_centroids(K, head_vec, centroids);
        if (updated_centroids == NULL){
          return NULL;
        }
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


static PyObject* wam(PyObject *self, PyObject *args){
    int len;
    double **data_arr,**res_as_arr;
    PyObject *data_list, *res_as_list;
    if(!PyArg_ParseTuple(args, "Oi", &data_list,&len)) {
        return NULL; 
    }
    data_arr = py_to_c_arr(data_list,len,len);
    res_as_arr = wam_c(data_arr,len,len);
    res_as_list = c_to_py_list(res_as_arr, len, len);
    free_arr(data_arr,len);
    free_arr(res_as_arr,len);
    return res_as_list;
}

static PyObject* ddg(PyObject *self, PyObject *args){
    int len;
    double **data_arr,**res_as_arr;
    PyObject *data_list, *res_as_list;
    if(!PyArg_ParseTuple(args, "Oi", &data_list,&len)) {
        return NULL; 
    }
    data_arr = py_to_c_arr(data_list,len,len);
    res_as_arr = ddg_c(data_arr,len, len);
    res_as_list = c_to_py_list(res_as_arr, len, len);
    free_arr(data_arr,len);
    free_arr(res_as_arr,len);
    return res_as_list;
}

static PyObject* gl(PyObject *self, PyObject *args){
    int len;
    double **data_arr,**res_as_arr;
    PyObject *data_list, *res_as_list;
    if(!PyArg_ParseTuple(args, "Oi", &data_list,&len)) {
        return NULL; 
    }
    data_arr = py_to_c_arr(data_list,len,len);
    res_as_arr = gl_c(data_arr,len,len);
    res_as_list = c_to_py_list(res_as_arr, len, len);
    free_arr(data_arr,len);
    free_arr(res_as_arr,len);
    return res_as_list;
}

static PyObject* jacobi(PyObject *self, PyObject *args){
    int len,sort;
    double **data_arr,**res_as_arr;
    PyObject *data_list, *res_as_list;
    if(!PyArg_ParseTuple(args, "Oii", &data_list,&len,&sort)) {
        return NULL; 
    }
    data_arr = py_to_c_arr(data_list,len,len);
    res_as_arr = jacobi_c(data_arr,len,sort);
    res_as_list = c_to_py_list(res_as_arr, len+1, len);
    free_arr(data_arr,len);
    free_arr(res_as_arr,len+1);
    return res_as_list;
}

static PyObject* spk(PyObject *self, PyObject *args){
    int K, data_length, vec_len;
    double **final_centroids, **data_arr;
    PyObject *data_list, *res;
    struct vector *head_vec;
    if(!PyArg_ParseTuple(args, "iOii", &K, &data_list,&data_length, &vec_len)) {
        return NULL; 
    }
    data_arr = py_to_c_arr(data_list, data_length, vec_len);

    head_vec = initialize_head_vec (data_arr, data_length, vec_len);
    if (head_vec == NULL){
        return NULL;
    }

    final_centroids = kmeans(K, head_vec, vec_len);
    if (final_centroids == NULL){
      return NULL;
    }
    
    free_arr(data_arr,data_length);
    res = c_to_py_list (final_centroids, K, vec_len);
    free_arr(final_centroids,K);
    return res;
}


static PyMethodDef spkmeansMethods[] = {
    {"wam",  
        (PyCFunction) wam, 
      METH_VARARGS, 
      PyDoc_STR("wam doc")},
    {"ddg",  
        (PyCFunction) ddg, 
      METH_VARARGS, 
      PyDoc_STR("ddg doc")},
    {"gl",  
        (PyCFunction) gl, 
      METH_VARARGS, 
      PyDoc_STR("gl doc")},
    {"jacobi",  
        (PyCFunction) jacobi, 
      METH_VARARGS, 
      PyDoc_STR("gl doc")},
    {"spk",  
        (PyCFunction) spk, 
      METH_VARARGS, 
      PyDoc_STR("spk doc")},
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

static struct PyModuleDef spkmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", 
    NULL, 
    -1, 
    spkmeansMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&spkmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}
    
    