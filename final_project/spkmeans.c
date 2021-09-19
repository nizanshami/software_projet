#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/**
 * Normalized Spectral Clustering implementation
 */

const int MAX_ITER = 300;

/**
 * find the eigengap and return the number of cluster
 */
int calc_eigenvalue_gap(double *mat, int *sorted_eigenvalues_indexes, int n) {
    double *deltas;
    double get(double *, int, int, int), max = 0;
    int i, index, next_index, half_n, result = 0;

    deltas = malloc(sizeof(double) * n);
    if(deltas == NULL){
        printf("An Error Has Occured");
        exit(1);
    }

    for (i = 0; i < n - 1; i++) {
        index = sorted_eigenvalues_indexes[i];
        next_index = sorted_eigenvalues_indexes[i + 1];
        deltas[i] = fabs(get(mat, index, index, n) -
                                get(mat, next_index, next_index, n));
    }
    half_n = (int) (n / 2);
    for (i = 0; i < half_n; i++) {
        if (max < deltas[i]) {
            max = deltas[i];
            result = i + 1;
        }
    }
    free(deltas);
    return result;
}

/**
 * given A real symmetric matrix A and a unit matrix V, updates A and V
 * with A to be an eigenvalues diagonal matrix,
 * and V's columns to be the corresponding eigenvectors.
 */
void jacobi_algorithm_for_eigenvalues(double *A, double *V, int n) {
    double get(double *, int, int, int), get_c(double), get_t(double *, int, int, int), calc_off_square(double *, int);
    double off_of_A, off_of_Atag, epsilon, s, c, t, val_row, val_col, a_ii, a_jj, a_ij;
    void indexes_of_max_off_diag(double *, int *, int *, int), set(double *, int, int, int, double);
    int row, col, i, counter = 0;

    /* if A diagonal matrix we skip the while loop with this starting settings*/
    off_of_A = 0;
    off_of_Atag = calc_off_square(A, n);
    epsilon = pow(10, -15);

    while (fabs(off_of_A - off_of_Atag) > epsilon && counter < 100) {
        counter++;
        off_of_A = off_of_Atag;
        indexes_of_max_off_diag(A, &row, &col, n);
        t = get_t(A, row, col, n);
        c = get_c(t);
        s = t * c;

        /*update A */
        for (i = 0; i < n; i++) {
            if (i != row && i != col) {
                val_row = c * get(A, i, row, n) - s * get(A, i, col, n);
                val_col = c * get(A, i, col, n) + s * get(A, i, row, n);
                set(A, i, row, n, val_row);
                set(A, row, i, n, val_row);
                set(A, i, col, n, val_col);
                set(A, col, i, n, val_col);
            }
        }
        a_ii = get(A, row, row, n);
        a_jj = get(A, col, col, n);
        a_ij = get(A, row, col, n);
        set(A, row, row, n, ((c * c) * a_ii + (s * s) * a_jj - 2 * s * c * a_ij));
        set(A, col, col, n, ((s * s) * a_ii + (c * c) * a_jj + 2 * s * c * a_ij));
        set(A, row, col, n, (((c * c) - (s * s)) * a_ij + s * c * (a_ii - a_jj)));
        set(A, col, row, n, (((c * c) - (s * s)) * a_ij + s * c * (a_ii - a_jj)));

        off_of_Atag = calc_off_square(A, n);

        /*update V */
        for (i = 0; i < n; i++) {
            val_row = get(V, i, row, n);
            val_col = get(V, i, col, n);

            set(V, i, row, n, c * val_row - s * val_col);
            set(V, i, col, n, c * val_col + s * val_row);
        }
    }
}

/**
 * form a unit matrix in mat
 * */
void form_unit_matrix(double *mat, int n) {
    void set(double *, int, int, int, double);
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                set(mat, i, j, n, 1.0);
            } else {
                set(mat, i, j, n, 0.0);
            }
        }
    }
}

/**
 * given the weighted adjacency matrix,
 * form the diagonal degree matrix into diagonal_mat
 */
void form_diagonal_mat(double *diagonal_mat, double *weighted_adj_mat, int n) {
    int i, j;
    double d = 0;
    for (i = 0; i < n; i++) {
        d = 0;
        for (j = 0; j < n; j++) {
            d += weighted_adj_mat[0];
            weighted_adj_mat++;
        }
        diagonal_mat[i] = d;
    }
}

/**
 * given data points, form it's weighted adjacency matrix into mat
 */
void form_weighted_adj_mat(double *mat, double *data_points, int dim, int n) {
    double distance(double *, double *, int, int, int);
    double get(double *, int, int, int), w;
    void set(double *, int, int, int, double);
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                set(mat, i, j, n, 0);
            } else {
                w = distance(data_points, data_points, dim, i, j) / 2;
                w = exp(-w);
                set(mat, i, j, n, w);
            }
        }
    }
}

/**
 * given D=diag(diagonal_mat) and W=weights_mat
 * calculates I-DWD into normalized_laplacian
 */
void calc_normalized_laplacian(double *normalized_laplacian, double *diagonal_mat, double *weights_mat, int dim) {
    void set(double *, int, int, int, double);
    double get(double *, int, int, int), result = 0;
    int i, j;

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if (i == j) {
                result = 1 - (1 / diagonal_mat[i]) * get(weights_mat, i, i, dim);
            } else {
                result = (-1)
                         * (1 / sqrt(diagonal_mat[i]))
                         * get(weights_mat, i, j, dim)
                         * (1 / sqrt(diagonal_mat[j]));
            }
            set(normalized_laplacian, i, j, dim, result);
        }
    }
}

/**
 * given a square matrix (mat), sets row and col to be
 * the indexes of the off-diagonal largest absolute value
 */
void indexes_of_max_off_diag(double *mat, int *row, int *col, int dim) {
    double max = 0, val = 0, get(double *, int, int, int);
    int i, j;

    for (i = 0; i < dim; i++) {
        for (j = i + 1; j < dim; j++) {
            val = get(mat, i, j, dim);
            if (fabs(val) > max) {
                max = fabs(val);
                *row = i;
                *col = j;
            }
        }
    }
}

/**
 * given a square matrix (mat) ,
 * the indexes of its off-diagonal largest absolute value (i,j)
 * and it's dimension (dim)
 * returns t
 */
double get_t(double *mat, int i, int j, int dim) {
    double theta, t, get(double *, int, int, int);

    theta = (get(mat, j, j, dim) - get(mat, i, i, dim)) / (2 * get(mat, i, j, dim));
    t = 1 / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    if (theta < 0) {
        t = (-1) * t;
    }
    return t;
}

/**
 * given t, returns c
 */
double get_c(double t) {
    return 1 / (sqrt(pow(t, 2) + 1));
}

/**
 * given matrix mat calculates the off(mat)^2 = notm_f(mat) - sum_(i=1)^n a_ii
 */
double calc_off_square(double *mat, int n) {
    double result = 0;
    double get(double *, int, int, int);
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                result += pow(get(mat, i, j, n), 2);
            }
        }
    }
    return result;
}

/**
 * given diagonal matrix [n*n] and an empty array [n], fills the array with the diagonal values
 */
void diag_to_array(double *diag, double *arr, int n) {
    double get(double *, int, int, int);
    int i;
    for (i = 0; i < n; i++) {
        arr[i] = get(diag, i, i, n);
    }
}

/**
 * given array of length dim, and its dimention, fills the array with integers
 * such that arr[i]=i
 */
void fill_ordered_ints(int *arr, int dim) {
    int i;
    for (i = 0; i < dim; i++) {
        arr[i] = i;
    }
}

/**
 * given an array- arr, and another array- indexes, both of length dim,
 * fills indexes with integers [0,dim-1] such that their order dictates a sorted order on arr values.
 * (using stable QuickSort algorithm)
 */
void quickSort_indexes(double *arr, int *indexes, int dim) {
    void fill_ordered_ints(int *, int);
    void quickSort_rec(double *, int *, int, int);
    void stable(double *, int *, int);

    fill_ordered_ints(indexes, dim);
    quickSort_rec(arr, indexes, 0, dim - 1);
    stable(arr, indexes, dim);
}

/**
 * Recursive QuickSort
 */
void quickSort_rec(double *arr, int *indexes, int low, int high) {
    void quickSort_rec(double *, int *, int, int);
    int partition(double *, int *, int, int);
    int mid = 0;

    if (low < high) {
        mid = partition(arr, indexes, low, high);
        quickSort_rec(arr, indexes, low, mid - 1);
        quickSort_rec(arr, indexes, mid + 1, high);
    }
}

/**
 * Implementation of the partition algorithm that is used for QS
 */
int partition(double *arr, int *indexes, int low, int high) {
    double pivot = arr[indexes[high]];
    int i = low - 1, j, temp;

    for (j = low; j <= high - 1; j++) {
        if (arr[indexes[j]] < pivot) {
            i++;
            temp = indexes[i];
            indexes[i] = indexes[j];
            indexes[j] = temp;
        }
    }
    temp = indexes[i + 1];
    indexes[i + 1] = indexes[high];
    indexes[high] = temp;

    return i + 1;
}

/**
 * Stables the sort
 */
void stable(double *arr, int *indexes, int dim) {
    void arr_int_to_double(double *, int *, int);
    void sorted_double_to_int(double *, int *, int *, int);
    void fill_ordered_ints(int *, int);
    double *double_indexes;
    int *ordered_ints;
    int low = 0, high = 1;

    double_indexes = malloc(sizeof(double) * (dim));
    if(double_indexes == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    arr_int_to_double(double_indexes, indexes, dim);

    ordered_ints = malloc(sizeof(int) * (dim));
    if(ordered_ints == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    fill_ordered_ints(ordered_ints, dim);
    if (dim <= 1) {
        return;
    }
    while (high < dim) {
        while (high + 1 < dim
               && arr[(int) double_indexes[high + 1]] == arr[(int) double_indexes[low]]) {
            high += 1;
        }
        if (arr[(int) double_indexes[high]] != arr[(int) double_indexes[low]]) {
            low += 1;
            high = low + 1;
            continue;
        }
        quickSort_rec(double_indexes, ordered_ints, low, high);
        low = high + 1;
        high = low + 1;
    }
    sorted_double_to_int(double_indexes, indexes, ordered_ints, dim);

    free(ordered_ints);
    free(double_indexes);
}

/**
 * Convert integers array into doubles array
 */
void arr_int_to_double(double *d_arr, int *i_arr, int dim) {
    int i;
    for (i = 0; i < dim; i++) {
        d_arr[i] = (double) i_arr[i];
    }
}

/**
 * Convert doubles array into integers array,
 * while sorting it according to the sorted indexes array
 */
void sorted_double_to_int(double *d_arr, int *i_arr, int *sorted_indexes, int dim) {
    int i;
    for (i = 0; i < dim; i++) {
        i_arr[i] = (int) d_arr[sorted_indexes[i]];
    }
}

/**
 * Fill given matrix with the data from input file
 */
double * read_data(FILE *fp, int *n, int *dim) {
    void set(double *, int, int, int, double);
    int num_of_lines(FILE*);
    int num_of_columns(FILE*);
    char *p, *line;
    long bOfFile;
    int size, i, j;
    double *data_points;

    bOfFile = ftell(fp);/*save the address of the beginning of the file */
    *n = num_of_lines(fp);
    fseek(fp, bOfFile, SEEK_SET);/*set the file position back to the beginning */
    *dim = num_of_columns(fp);
    fseek(fp, bOfFile, SEEK_SET);/*set the file position back to the beginning */

    line = malloc(sizeof(char) * (30 * (*dim)));
    if(line == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    data_points = malloc(sizeof(double) * (*dim) * (*n));
    if(data_points == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    size = 50 * (*dim);
    for (i = 0; i < *n; i++) {

        /* p is a pointer to a beginning of a line in the file */
        p = fgets(line, size, fp);

        for (j = 0; j < *dim; j++) {

            /* extract double  form the line */
            set(data_points, i, j, *dim, strtod(p, &p));

            p += 1; /* skip comma */
        }
    }
    free(line);
    return data_points;
}

/**
 * Counts the lines in input file
 */
int num_of_lines(FILE *fp) {
    int ch;
    int lines = 1;

    while (!feof(fp)) {
        ch = fgetc(fp);
        if (ch == '\n') {
            lines++;
        }
    }
    return lines;
}

/**
 * Count the columns in input file
 * */
int num_of_columns(FILE *fp) {
    int ch;
    int columns = 1;

    while (!feof(fp)) {
        ch = fgetc(fp);
        if (ch == ',') {
            columns++;
        }
        if (ch == '\n') {
            break;
        }
    }
    return columns;
}

/**
 * Prints matrix in the template required,
 * if row == -1 that means mat is an array containing the diagonal elements
 * of the diagonal matrix to print
 */
void print_matrix(double *mat, int row, int col) {
    int i, j, diag = 0;
    double item;
    double get(double *, int, int, int);
    if (row == -1) {
        diag = 1;
        row = col;
    }

    for (i = 0; i < row; i++) {
        for (j = 0; j < col ; j++) {
            if (diag == 1) {
                item = mat[i];
                if (i == j && fabs(item) >= 0.00005) {
                    printf("%0.4f", item);
                } else {
                    printf("%0.4f", 0.0);
                }
            } else {
                item = get(mat, i, j, col);
                if(fabs(item) >= 0.00005){
                    printf("%0.4f", item);
                } else {
                    printf("%0.4f", 0.0);
                }
            }

            /*adding comma when needed*/
            if (j < col-1 ){
                printf(",");
            }
        }
        /*adding line break when needed*/
        if (i < row-1){
            printf("\n");
        }
    }
}

/**
 * Given the eigenvectors matrix [n*n] (as columns),
 * an array [n] that dictates an order on the eigenvalues,
 * n (the mat dimention), k, and a result matrix [n*k],
 * fills the result matrix with the normalized first k eigenvectors as columns
 */
void normalized_k_eigenvectors(double *eigenvectors, int *indexes, int n, int k, double *result) {
    double get(double *, int, int, int);
    void set(double *, int, int, int, double);
    int i, j, t;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        sum = 0.0;
        for (t = 0; t < k; t++) {
            j = indexes[t];
            sum += get(eigenvectors, i, j, n) * get(eigenvectors, i, j, n);
        }
        sum = sqrt(sum);
        for (t = 0; t < k; t++) {
            j = indexes[t];
            set(result, i, t, k, get(eigenvectors, i, j, n) / sum);
        }
    }
}

/**
 * Calculate distance between two vectors sqrt(sum_(i=1)^n (v1_i-v2_i)^2) (L2 norm)
 */
double distance(double *v1, double *v2, int dim, int row_v1, int row_v2) {
    int i;
    double result = 0;
    double x;
    double get(double *, int, int, int);

    for (i = 0; i < dim; i++) {
        x = (get(v1, row_v1, i, dim) - get(v2, row_v2, i, dim));
        x *= x;
        result += x;
    }
    result = sqrt(result);
    return result;
}

/**
 * Calculate distance between two vectors sqrt(sum_(i=1)^n (v1_i-v2_i)^2) (L2 norm)
 */
double kmeans_distance(double *v1, double *v2, int dim, int row_v1, int row_v2) {
    int i;
    double result = 0;
    double x;
    double get(double *, int, int, int);

    for (i = 0; i < dim; i++) {
        x = (get(v1, row_v1, i, dim + 1) - get(v2, row_v2, i, dim));
        x *= x;
        result += x;
    }
    return result;
}

/**
 * Returns arr[i][j]
 */
double get(double *arr, int i, int j, int dim) {
    int index;

    index = (i * dim + j);
    return arr[index];
}

/**
 * Sets arr[j][i]=item
 */
void set(double *arr, int i, int j, int dim, double item) {
    int index;

    index = (i * dim + j);
    arr[index] = item;
}

/**
 * K-means algorithm implementation
 */
int kmeans(int k, double *data_points, double *centroids, double *utl, int dim, int n) {
    void assign(double *, double *, int, int, int);
    short re_estimate(double *, double *, double *, int, int, int);
    short convergence;
    int i;

    for (i = 0; i < MAX_ITER; i++) {
        assign(data_points, centroids, dim, n, k);
        convergence = re_estimate(data_points, centroids, utl, dim, n, k);
        if (convergence == 1) {
            return 0;
        }
    }
    return 0;
}

/**
 * Assigns data points to their closest cluster (measure distance from the centroid)
 * updates the number of cluster for each data point
 */
void assign(double *data_points, double *clusters, int dim, int n, int k) {
    int int_max = 2147483647;
    int cluster = 0;
    int v, c;
    void set(double *, int, int, int, double);
    double min_dis, dis, kmeans_distance(double *, double *, int, int, int);

    min_dis = int_max;
    for (v = 0; v < n; v++) {
        for (c = 0; c < k; c++) {
            dis = kmeans_distance(data_points, clusters, dim, v, c);
            if (dis <= min_dis) {
                min_dis = dis;
                cluster = c;
            }
        }
        set(data_points, v, dim, dim + 1, cluster);
        min_dis = int_max;
    }
}

/**
 * Re-estimates a centroid for each cluster:
 * for each cluster calculate the average of the points assign to it,
 * updates centroids to be the average vector,
 * returns 1 if the old centroids are equal to the new ones.
 */
short re_estimate(double *data_points, double *clusters, double *utl, int dim, int n, int k) {
    short isEqual = 1;
    int i, j;
    double x, get(double *, int, int, int);
    void zero_mat(double *, int, int),
            set(double *, int, int, int, double),
            vec_sum(double *, double *, int, int, int);

    zero_mat(utl, dim + 1, k);

    /* sum all vectors for each cluster */
    for (i = 0; i < n; i++) {
        j = get(data_points, i, dim, dim + 1);
        vec_sum(utl, data_points, dim, j, i);
        x = get(utl, j, dim, dim + 1) + 1;
        set(utl, j, dim, dim + 1, x);
    }

    /* Divides each sum by the number of vectors to get average */
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            x = get(utl, i, j, dim + 1);
            set(utl, i, j, dim + 1, (x / get(utl, i, dim, dim + 1)));
        }
    }

    /* Compare the old centroids to the new ones */
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            if (!(((get(clusters, i, j, dim) - get(utl, i, j, dim + 1)) < 0.000001) &&
                  ((get(clusters, i, j, dim) - get(utl, i, j, dim + 1)) > -0.000001))) {
                isEqual = 0;
                break;
            }
        }
        if (!isEqual) {
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
            x = get(utl, i, j, dim + 1);
            set(clusters, i, j, dim, x);
        }
    }
    return isEqual;
}

/**
 * Adds vec2 to vec1 coordinate wise
 */
void vec_sum(double *vec1, double *vec2, int dim, int row_vec1, int row_vec2) {
    int i;
    void set(double *, int, int, int, double);
    double sum, get(double *, int, int, int);

    for (i = 0; i < dim; i++) {
        sum = get(vec1, row_vec1, i, dim + 1) + get(vec2, row_vec2, i, dim + 1);
        set(vec1, row_vec1, i, dim + 1, sum);
    }
}

/**
 * Zeros a given matrix
 */
void zero_mat(double *mat, int dim, int n) {
    int i, j;
    void set(double *, int, int, int, double);

    for (i = 0; i < n; i++) {
        for (j = 0; j < dim; j++) {
            set(mat, i, j, dim, 0);
        }
    }
}

/**
 * Return target matrix for goal == "wam"
 */
double *wam(double *data_points, int n, int dim) {
    void form_weighted_adj_mat(double *, double *, int, int);
    double *target_matrix;

    target_matrix = malloc(sizeof(double) * n * n);
    if(target_matrix == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    form_weighted_adj_mat(target_matrix, data_points, dim, n);

    return target_matrix;
}

/**
 * Return an array that contains the diagonal of the target matrix for goal == "ddg"
 */
double *ddg(double *data_points, int n, int dim) {
    void form_diagonal_mat(double *, double *, int);
    double *target_diagonal, *weighted_adj_mat, *wam(double *, int, int);

    target_diagonal = malloc(sizeof(double) * n);
    if(target_diagonal == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    weighted_adj_mat = wam(data_points, n, dim);
    form_diagonal_mat(target_diagonal, weighted_adj_mat, n);

    free(weighted_adj_mat);

    return target_diagonal;
}

/**
 * Return target matrix for goal == "lnorm"
 */
double *lnorm(double *data_points, int n, int dim) {
    void form_diagonal_mat(double *, double *, int),
            form_weighted_adj_mat(double *, double *, int, int);
    void calc_normalized_laplacian(double *, double *, double *, int);
    double *target_matrix, *weighted_adj_mat, *diagonal_mat;

    target_matrix = malloc(sizeof(double) * n * n);
    weighted_adj_mat = malloc(sizeof(double) * n * n);
    diagonal_mat = malloc(sizeof(double) * n);
    if(target_matrix == NULL || weighted_adj_mat == NULL || diagonal_mat == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    form_weighted_adj_mat(weighted_adj_mat, data_points, dim, n);
    form_diagonal_mat(diagonal_mat, weighted_adj_mat, n);
    calc_normalized_laplacian(target_matrix, diagonal_mat, weighted_adj_mat, n);

    free(weighted_adj_mat);
    free(diagonal_mat);

    return target_matrix;
}

/**
 * Return target matrix in shape (n+1,n)
 * the first row of the returned matrix will be the eigenvalues,
 * the other n rows will be the eigenvectors matrix when the goal == "jacobi"
 * */
double *jacobi(double *data_points, int n) {
    double *target_matrix, *V;
    void jacobi_algorithm_for_eigenvalues(double *, double *, int),
            form_unit_matrix(double *, int);
    int i, j;

    target_matrix = malloc(sizeof(double) * n * (n + 1));
    V = malloc(sizeof(double) * n * n);
    if(target_matrix == NULL || V == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    form_unit_matrix(V, n);
    jacobi_algorithm_for_eigenvalues(data_points, V, n);

    for (i = 0; i <= n; i++) {
        for (j = 0; j < n; j++) {
            if (i == 0) {
                set(target_matrix, i, j, n, get(data_points, j, j, n));
            } else {
                set(target_matrix, i, j, n, get(V, j, (i - 1), n));
            }
        }
    }
    free(V);

    return target_matrix;
}

/**
 * Return target matrix T in shape (n,k) when goal == "spk"
 */
double *spk(double *data_points, int n, int dim, int *k) {
    double *target_matrix, *lnorm(double *, int, int), *V,
            *normalized_laplacian, *eigenvalues;
    void jacobi_algorithm_for_eigenvalues(double *, double *, int),
            form_unit_matrix(double *, int);
    void diag_to_array(double *, double *, int), quickSort_indexes(double *, int *, int)
    , normalized_k_eigenvectors(double *, int *, int, int, double *);
    int *indexes;

    V = malloc(sizeof(double) * n * n);
    indexes = malloc(sizeof(int) * n);
    eigenvalues = malloc(sizeof(double) * n);
    if(V == NULL || indexes == NULL || eigenvalues == NULL){
        printf("An Error Has Occured");
        exit(1);
    }

    normalized_laplacian = lnorm(data_points, n, dim);
    form_unit_matrix(V, n);
    jacobi_algorithm_for_eigenvalues(normalized_laplacian, V, n);
    diag_to_array(normalized_laplacian, eigenvalues, n);
    quickSort_indexes(eigenvalues, indexes, n);
    if (*k == 0) {
        *k = calc_eigenvalue_gap(normalized_laplacian, indexes, n);
    }
    if (*k <= 0 || *k >= n) {
        printf("invalid input!");
        exit(1);
    }
    target_matrix = malloc(sizeof(double) * n * (*k));
    if(target_matrix == NULL){
        printf("An Error Has Occured");
        exit(1);
    }
    normalized_k_eigenvectors(V, indexes, n, *k, target_matrix);

    free(normalized_laplacian);
    free(indexes);
    free(eigenvalues);
    free(V);

    return target_matrix;
}

typedef enum {WAM, DDG, LNORM, JACOBI, SPK} Goal;

/**
 * Given string, returns the respective enum (default enum is spk)
 */
Goal get_enum(char* goal_string){
    if (!strcmp(goal_string, "wam")){
        return WAM;
    } else if (!strcmp(goal_string, "ddg")){
        return DDG;
    } else if (!strcmp(goal_string, "lnorm")){
        return LNORM;
    } else if (!strcmp(goal_string, "jacobi")){
        return JACOBI;
    } else if (!strcmp(goal_string, "spk")) {
        return SPK;
    } else {
        return SPK;
    }
}

/**
 * Copies the fist last_row_to_copy rows of mat_to_copy to result_mat
 */
void copy_rows(double* result_mat, double* mat_to_copy, int last_row_to_copy, int dim1, int dim2){
    int i, j;
    for (i = 0; i < last_row_to_copy; i++) {
        for (j = 0; j < dim2; j++) {
            set(result_mat, i, j, dim1, get(mat_to_copy, i, j, dim2));
        }
    }
}

int main(int argc, char *argv[]) {
    void print_matrix(double *, int, int);
    double* read_data(FILE *, int *, int *);
    void copy_rows(double*, double*, int, int, int);
    int kmeans(int, double *, double *, double *, int, int);
    double *target_matrix, *data_points, *centroids, *util, *transform_points,
            *wam(double *, int, int), *lnorm(double *, int, int), *ddg(double *, int, int),
            *spk(double *, int, int, int *), *jacobi(double *, int);
    Goal get_enum(char*), goal;
    int dim=0, k, n=0, rows, cols;
    FILE *f;
    target_matrix = 0;

    if (argc != 4) {
        printf("invalid input!");
        return 1;
    }
    /* reading arguments*/
    k = strtol(argv[1], NULL, 10);
    goal = get_enum(argv[2]);
    f = fopen(argv[3], "r");

    /* build matrix that contains all the points */
    data_points = read_data(f, &n, &dim);
    fclose(f);

    rows = k;
    cols = k;

    /*calculate the goal matrix*/
    switch (goal){
        case WAM:
            target_matrix = wam(data_points, n, dim);
            rows = n;
            cols = n;
            break;
        case DDG:
            target_matrix = ddg(data_points, n, dim);
            rows = -1;
            cols = n;
            break;
        case LNORM:
            target_matrix = lnorm(data_points, n, dim);
            rows = n;
            cols = n;
            break;
        case JACOBI:
            target_matrix = jacobi(data_points, n);
            rows = (n+1);
            cols = n;
            break;
        case SPK:
            target_matrix = spk(data_points, n, dim, &k);
            rows = k;
            cols = k;
            break;
    };

    free(data_points);
    if (goal != SPK) {
        print_matrix(target_matrix, rows, cols);
        free(target_matrix);
        return 0;
    }
    transform_points = malloc(sizeof(double) * n * (k + 1));
    centroids = malloc(sizeof(double) * k * k);
    util = malloc(sizeof(double) * k * (k + 1));
    if(transform_points == NULL || centroids == NULL || util == NULL){
        printf("An Error Has Occured");
        exit(1);
    }

    copy_rows(centroids, target_matrix ,k , k, k);
    copy_rows(transform_points, target_matrix ,n ,k+1 , k);

    kmeans(k, transform_points, centroids, util, k, n);
    print_matrix(centroids, k, k);

    /* free the memory used */
    free(target_matrix);
    free(transform_points);
    free(centroids);
    free(util);

    return 0;
}
