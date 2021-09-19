"""
Kmeans_pp implementation in Python
Software Project ex. 2

Amit Elyasi
316291434

Nizan Shemi
206962912
"""

import mykmeanssp as kmeans
import numpy as np
import pandas as pd
import sys


"""
implementation of kmeans++, finding a good initialization of centroids for the "kmeans" algorithm
"""
def Kmeans_pp(datapoints,k):
    n = len(datapoints)
    datapoints = np.array([np.array(point) for point in datapoints])
    centroids = [0]
    centroids[0] = np.random.choice(n,1)[0]
    dists = np.zeros((n,k))
    min_dists = np.zeros(n)
    for z in range(k-1):
        for i in range(n):
            dists[i][z] = np.inner(datapoints[i]-datapoints[centroids[z]], datapoints[i]-datapoints[centroids[z]])
            min_dists[i] = min(dists[i][0:z+1])
        dists_sum = np.sum(min_dists)
        probabilities = [dist/dists_sum for dist in min_dists]

        new_centroid_indx = np.random.choice(n, 1, p=probabilities)[0]
        centroids.append(new_centroid_indx)
    return centroids


"""
given 2 file names, where the files containes data points with the first column being the index, 
reads the files and joins them with inner join (key=the first column) 
returns pandas DataFrame containing the joint data 
"""
def readAndJoin(file1, file2):
    data1 = pd.read_csv(file1, header=None)
    data2 = pd.read_csv(file2, header=None)
    data = pd.merge(data1, data2, on=0, sort=True).drop(0,1)
    return data


"""
given an array of arrays, returns a 1D python array with all of the entrys of the matrix in one long vector
"""
def arrToSeq(arr):
    return [item for vec in arr for item in vec]


def seqRoundAndPrint(seq, d):
    n = len(seq)//d
    seq = np.round(seq, 4)
    for i in range(n):
        for j in range(d-1):
            print(f"{seq[i*d+j]},",end="")
        print(seq[i * d + d-1])


def main():
    np.random.seed(0)
    args = sys.argv
    if len(args)<4 or len(args)>5:
        print("INPUT ERROR:\nthe number of argument is incorrect")
        return 1
    if len(args) == 5:
        k, max_iter, file1, file2 = args[1], args[2], args[3], args[4]
    else:
        k, file1, file2 = args[1], args[2], args[3]
        max_iter = 300

    try:
        k = int(k)
    except ValueError:
        print("INPUT ERROR:\nk has to be an integer")
        return 1

    if (k <= 0):
        print("INPUT ERROR:\nk can't be <= 0")
        return 1
    try:
        max_iter = int(max_iter)
    except ValueError:
        print("INPUT ERROR:\nk can't be a letter")
        return 1
    if max_iter <= 0:
        print("INPUT ERROR:\nmax iteration can't be <= 0")
        return 1

    datapoints_table = readAndJoin(file1, file2)
    datapoints_matrix = datapoints_table.values.tolist()
    if len(datapoints_matrix) <= k:
        print(f"INPUT ERROR:\nthere are less or equal to {k} data points")
        return False
    centroids_indexes = Kmeans_pp(datapoints_matrix, k)
    datapoints = arrToSeq(datapoints_matrix)
    centroids = arrToSeq([datapoints_matrix[i] for i in centroids_indexes])
    print(",".join(str(indx) for indx in centroids_indexes))
    centroids = kmeans.fit(k, datapoints, centroids, max_iter, len(datapoints_matrix[0]), len(datapoints_matrix))
    seqRoundAndPrint(centroids, len(datapoints_matrix[0]))

if __name__ == "__main__":
    main()



