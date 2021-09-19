#include <stdio.h> 
extern const int MAX_ITER;

double get(double* arr, int i, int j, int dim);

/**
 * Sets arr[j][i]=item
 */
void set(double* arr, int i, int j, int dim, double item);

/**
 * Return target matrix for goal == "wam"
 */
double *wam(double *data_points, int n, int dim);

/**
 * Return an array that contains the diagonal of the target matrix for goal == "ddg"
 */
double *ddg(double *data_points,int n, int dim);

/**
 * Return target matrix for goal == "lnorm"
 */
double *lnorm(double *data_points, int n, int dim);

/**
 * Return target matrix in shape (n+1,n)
 * the first row of the returned matrix will be the eigenvalues,
 * the other n rows will be the eigenvectors matrix when the goal == "jacobi"
 * */
double *jacobi(double *data_points, int n);

/**
 * Return target matrix T in shape (n,k) when goal == "spk"
 */
double *spk(double *data_points, int n , int dim, int *k);

/**
 * K-means algorithm implementation
 */
int kmeans(int k, double* data_points, double* centroids, double* utl, int dim, int n);

typedef enum {WAM, DDG, LNORM, JACOBI, SPK} Goal;

/**
 * Prints matrix in the template required
 */
Goal get_enum(char* goal_string);
