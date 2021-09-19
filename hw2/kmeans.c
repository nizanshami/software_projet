#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
/*
 * Helper functions prototypes
 */
static int kmeans(int , float*, float*, float*, int, int, int );
static void assign(float *, float *, int, int, int);
static short re_estimate(float *, float*,float *, int, int, int);
static void vec_sum(float *, float *, int, int, int);
static void zero_mat(float *, int, int);
static float distance(float *, float *, int,int, int);
static float get(float*, int, int, int);
static void set(float* , int, int, int, float);
static PyObject *Convert_Big_Array(float *, int);

/*
* kmeans code start here
*/
static int kmeans(int k, float* data_points, float* centroids, float* utl ,int max_iter, int dim, int n){
    short convergece = 1;
    int i;

    for (i=0; i<max_iter; i++){
        assign(data_points, centroids, dim, n, k);
        convergece = re_estimate(data_points, centroids, utl, dim, n, k);
        if (convergece == 1) {
            return 0;
        }
    }
    return 0;
}


/*
 * assigns data points to their closest cluster (measure distance from the centroid)
 * updates the number of cluster for each data point
 */
static void assign(float* data_points, float* clusters, int dim, int n, int k){
    int int_max = 2147483647;
    int cluster;
    int v,c;
    float min_dis, dis;

    min_dis = int_max;
    for(v = 0; v < n; v++){
        for(c = 0;c < k; c++){
            dis = distance(data_points, clusters, dim, v, c);
            if( dis <= min_dis){
                min_dis = dis;
                cluster = c;
            }
        }
        set(data_points, v, dim, dim + 1, cluster);
        min_dis = int_max;
    }
}


/*
 * re-estimates a centroid for each cluster:
 * for each cluster calculate the average of the points assign to it,
 * updates centroids to be the average vector,
 * returns 1 if the old centroids are equal to the new ones.
 */
static short re_estimate(float* data_points, float* clusters,float *utl, int dim, int n, int k) {
    short isEqual = 1;
    int i, j;
    float x;

    zero_mat(utl, dim + 1, k);

    /* sum all vectors for each cluster */
    for (i = 0; i < n; i++) {
        j = get(data_points, i, dim, dim+1);
        vec_sum(utl, data_points, dim, j, i);
        x = get(utl, j, dim, dim + 1) + 1;
        set(utl, j, dim, dim + 1, x);
    }

    /* Divides each sum by the number of vectors to get average */
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            x = get(utl, i, j, dim+1); 
            set(utl, i, j, dim + 1, (x / get(utl, i, dim, dim+1)));
        }
    }

    /* Compare the old centroids to the new ones */
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            if (!(((get(clusters, i, j, dim) - get(utl, i, j, dim+1)) < 0.000001) && ((get(clusters, i, j, dim) - get(utl, i, j, dim+1)) > -0.000001))) {
                isEqual = 0;
                break;
            }
        }
        if (!isEqual){
            break;
        }
    }

    /* if there is no change in centroid- we have reach convergence */
    if (isEqual == 1) {
        return 1;
    }

    /* copy the new centroids to the old ones place */
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            x = get(utl, i, j, dim+1);
            set(clusters, i, j, dim, x);
        }
    }
    return isEqual;
}


/*
* calculate distance between two vectors (v1-v2)^2
*/
static float distance(float *v1, float *v2, int dim,int row_v1, int row_v2){
    int i;
    float result = 0;
    float x;
    
    for(i = 0;i < dim;i++){
        x = (get(v1, row_v1, i, dim + 1)-get(v2, row_v2, i, dim));
        x *= x;
        result += x;
    }
    return result;
}


/*
 * adds vec2 to vec1 coordinate wise
 */
static void vec_sum(float* vec1, float* vec2, int dim, int row_vec1, int row_vec2){
    int i;
    float sum;
    
    for(i = 0;i < dim;i++){
        sum = get(vec1, row_vec1, i, dim+1) + get(vec2, row_vec2, i, dim+1);
        set(vec1, row_vec1, i, dim + 1, sum);
    }
}


/*
 * zeros a given matrix from row start to row end
 */
static void zero_mat(float* clusters , int dim, int n){
    int i,j;

    for(i = 0; i < n; i++){
        for(j=0; j < dim; j++){
            set(clusters, i, j, dim, 0);
        }
    }
}



static float get(float* arr, int i, int j, int dim){
    int index;

    index = (i*dim + j);
    return arr[index];
}


static void set(float* arr, int i, int j, int dim, float item){
    int index;

    index = (i*dim + j);
    arr[index] = item;
}
/*
* python c extension start here
*/

static PyObject *Convert_Big_Array(float *array, int length){
    PyObject *pylist, *item;
    int i;
    pylist = PyList_New(length);
    for (i=0; i<length; i++) {
      item = PyFloat_FromDouble((double) array[i]);
      PyList_SetItem(pylist, i, item);
    }
    return pylist;
  }


static PyObject * fit_c(int k, PyObject *pyData_points, PyObject *pyCentroid, int max_iter, int dim, int n){
    float *data_points, *centroid, *utl,x;
    Py_ssize_t index;
    int i,j;
    PyObject *item, *pylist;

    data_points = malloc(sizeof(float) *((dim+1) *n));
    centroid = malloc(sizeof(float) * (k * dim));
    utl = malloc(sizeof(float) * (k * (dim+1)));

    /*convert python lists : pyData_points and pyCentroid to c arrays*/
    for(i = 0;i < n;i++){
        for(j = 0;j < dim;j++){
            index = i*dim+j;
            item = PyList_GetItem(pyData_points,index);
            x = (float) PyFloat_AsDouble(item);
            if(x == -1.0 && PyErr_Occurred()){
                puts("float error");
            }
            set(data_points, i, j, dim+1, x);

            if(i < k){
                item = PyList_GetItem(pyCentroid, index);
                x = (float) PyFloat_AsDouble(item);
                if(x == -1.0 && PyErr_Occurred()){
                puts("float error");
                }
                set(centroid, i, j, dim, x);
            }
        }
        set(data_points, i, dim, dim+1, 0);
    }


    kmeans(k,data_points, centroid, utl, max_iter, dim, n);

    pylist = Convert_Big_Array(centroid, k*dim);

    free(data_points);
    free(centroid);
    free(utl);

    return pylist;
} 

static PyObject *fit_capi(PyObject* self, PyObject* args){
    PyObject *pyData_points, *pyCentroid;
    int k, dim, n, max_iter;
    if(!PyArg_ParseTuple(args, "iOOiii", &k,
                                        &pyData_points,
                                        &pyCentroid,
                                        &max_iter,
                                        &dim,
                                        &n)){
        return NULL;
    }
    return fit_c(k, pyData_points, pyCentroid, max_iter, dim, n);
    
}


/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef mykmeansMethods[] = {
    {"fit",                   /* the Python method name that will be used */
      (PyCFunction) fit_capi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("compute centroid useing the givan initial points")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    mykmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};


PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
