# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "spkmeans.h"


double** py_to_c_arr(PyObject *list, int dim_1, int dim_2){//converts python list to c array 
    int i,j;
    PyObject *item;
    data_arr = malloc(dim_1 * sizeof(double*));
    if (data_arr == NULL){
      return NULL;
    }
    for (i=0; i < dim_2; i++){
      data_arr[i] = malloc(vec_len * sizeof(double));
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
            python_float = Py_BuildValue("f",final_centroids[i][j]);
            PyList_SetItem(PyList_GetItem(res,i), j, python_float);
        }
    }
    return res;
}

static PyObject* wam(PyObject *self, PyObject *args){
    int **data_arr,**res_as_arr, len;
    PyObject *data_list, *res_as_list;
    if(!PyArg_ParseTuple(args, "Oi", &data_list,&len)) {
        return NULL; 
    }
    data_arr = py_to_c_arr(data_list,len,len);
    free_arr(data_arr,len,len);
    res_as_arr = wam_c(data_arr,len);
    res_as_list = c_to_py_list();
    free_arr(res_as_arr,len,len);
    return res_as_list;
}

static PyMethodDef spkmeansMethods[] = {
    {"wam",  
        (PyCFunction) wam, 
      METH_VARARGS, 
      PyDoc_STR("wam doc")},
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
    
    