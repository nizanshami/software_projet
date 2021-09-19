/*
 * Kmeans implementation in C
 * Software Project ex. 1
 *
 * Amit Elyasi
 * 316291434
 *
 * Nizan Shemi
 * 206962912
 *
 */

#include <stdio.h>
#include <stdlib.h>


/*
 * implementations of kmeans algorithm:
 * given n data point (in R^d), group the data into k clusters,
 * each data point is assigned to exacly one cluster (note that k<n)
 */
int kmeans(int k, float* data_points, float* centroids, float* utl ,int max_iter, int dim, int n){
    short convergece = 1;
    void assign(float*, float*, int, int, int);
    short re_estimate(float*, float* , float*, int, int, int);
    void build_clusters(int, int, float*, float*);
    int i;

    build_clusters(k, dim, data_points, centroids);
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
void assign(float* data_points, float* clusters, int dim, int n, int k){
    int INT_MAX = 2147483647;
    int cluster;
    float distance(float *, float *, int, int, int);
    int v,c;
    float min_dis, dis;
    void set(float *, int, int , int,float);

    min_dis = INT_MAX;
    for(v = 0; v < n; v++){
        for(c = 0;c < k; c++){
            dis = distance(data_points, clusters, dim, v, c);
            if( dis <= min_dis){
                min_dis = dis;
                cluster = c;
            }
        }
        set(data_points, v, dim, dim + 1, cluster);
        min_dis = INT_MAX;
    }
}


/*
 * re-estimates a centroid for each cluster:
 * for each cluster calculate the average of the points assign to it,
 * updates centroids to be the average vector,
 * returns 1 if the old centroids are equal to the new ones.
 */
short re_estimate(float* data_points, float* clusters,float *utl, int dim, int n, int k) {
    void vec_sum(float* , float* , int, int, int);
    void zero_mat(float *, int, int);
    float get(float *, int, int, int);
    void set(float *, int, int, int, float);
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
float distance(float *v1, float *v2, int dim,int row_v1, int row_v2){
    int i;
    float result = 0;
    float x;
    float get(float *, int, int, int);

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
void vec_sum(float* vec1, float* vec2, int dim, int row_vec1, int row_vec2){
    float get(float *, int, int, int);
    void set(float *, int ,int, int, float);
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
void zero_mat(float* clusters , int dim, int n){
    int i,j;
    void set(float *, int, int, int, float);

    for(i = 0; i < n; i++){
        for(j=0; j < dim; j++){
            set(clusters, i, j, dim, 0);
        }
    }
}


/* fill given matrix with the data from input file */
void read_data(FILE* fp, float *data_points, char *line, int n, int dim){
    void set(float *, int, int, int, float);
    char *p;
    int size, i, j;

    size = 50 * dim;
    for(i = 0; i < n; i++){

        /* p is a pointer to a beging of a line in the file */
        p = fgets(line, size, fp);

        for(j = 0; j < dim; j++){

            /* extract float form the line */
            set(data_points, i, j, dim + 1, strtod(p, &p));

            p += 1; /* skip comma */
        }
        set(data_points, i, j, dim + 1, 0);
    }
}


/*
 * counts the lines in input file
 */
int num_of_lines(FILE *fp){
    int ch;
    int lines = 0;

    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
        {
            lines++;
        }
    }
    return lines;
}


/* count the columns in input file */
int num_of_columns(FILE *fp){
    int ch;
    int columns = 1;

    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == ','){
            columns++;
        }
        if(ch == '\n'){
            break;
        }
    }
    return columns;
}


void build_clusters(int k, int dim, float *vectors, float *centroids){
    void set(float *, int, int, int, float);
    float get(float *, int, int ,int);
    int i, j;

    for(i = 0; i < k; i++){
        for(j = 0;j < dim;j++){
            set(centroids, i, j, dim, get(vectors, i, j, dim + 1));
        }
    }
}


/*
 * prints the centroids in the template requierd
 */
void print_centroids(float* clusters, int k, int dim){
    int i,j;
    float get(float *, int, int, int);

    for(i = 0; i < k;i++){
        for(j = 0;j < dim-1;j++){
            printf("%0.4f,", get(clusters, i, j, dim));
        }
        printf("%.4f\n", get(clusters, i, dim-1, dim));
    }
}


float get(float* arr, int i, int j, int dim){
    int index;

    index = (i*dim + j);
    return arr[index];
}


void set(float* arr, int i, int j, int dim, float item){
    int index;

    index = (i*dim + j);
    arr[index] = item;
}


int main( int argc, char* argv[]) {
    int kmeans(int, float* , float*, float*, int, int, int);
    void read_data(FILE*, float *, char *, int, int );
    void print_centroids(float*, int, int);
    int max_iter, dim, k, n;
    long bOfFile;
    float * utl;
    float *data_points;
    float *centroids;
    char *line;

    /* reading arguments */
    k = strtol(argv[1], NULL, 10);
    if(argc == 3){
        max_iter = strtol(argv[2], NULL, 10);

        if (max_iter<=0){
            printf("INPUT ERROR:\nmaximum iterations is invalid");
            return 1;
        }
    }

    else{
        max_iter = 200;
    }
    bOfFile = ftell(stdin);/*save the address of the beginning of the file */
    n = num_of_lines(stdin);

    if(k<=0){
        printf("INPUT ERROR:\nk is invalid");
        return 1;
    }
    if(n < k){
        printf("INPUT ERROR:\nthere are less then k=%d data points",k);
        return 1;
    }

    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */
    dim = num_of_columns(stdin);
    fseek(stdin, bOfFile, SEEK_SET);/*set the file position back to the beginning */

    line = malloc(sizeof(char) * (30*dim));
    /* build matrix that contins all the points */;
    data_points = malloc(sizeof(float) * ((dim+1) * n));
    /* build matrix that contins all the centroids */;
    centroids = malloc(sizeof(float) * (dim*k));
    /* build matrix that use for calclations */;
    utl = malloc(sizeof(float) * ((dim+1)*k));
    read_data(stdin, data_points, line, n, dim);
    kmeans(k, data_points, centroids, utl, max_iter, dim, n);
    print_centroids(centroids, k, dim);


    /* free the memory used */
    free(line);
    free(data_points);
    free(centroids);
    free(utl);

    return 0;
}