"""
Kmeans implementation in Python
Software Project ex. 1

Amit Elyasi
316291434

Nizan Shemi
206962912
"""

import sys


def Kmeans(k, max_iter=200):

    data_points = read_data()
    if not good_input(data_points, k):
        return

    centroids = data_points[:k]
    for i in range(max_iter):
        clusters = assign(data_points, centroids)
        new_centroids = re_estimate(clusters)
        if sorted(new_centroids) == sorted(centroids):
            break
        centroids = new_centroids
    return centroids


def assign(data_points, centroids):
    clusters = [[centroid] for centroid in centroids]

    for point in data_points:
        min_dst = distance(point, centroids[0])
        closest_centroid_num = 0

        for i,centroid in enumerate(centroids):
            dst = distance(point, centroid)
            if dst < min_dst:
                min_dst = dst
                closest_centroid_num = i

        clusters[closest_centroid_num].append(point)

    return clusters


def re_estimate(clusters):
    new_centroids = []
    for cluster in clusters:
        cluster = cluster[0:]
        new_centroids.append(average(cluster))

    return new_centroids


def read_data():
    vectors = []
    while(True):
        try:
            for line in input().split("\n"):
                vectors.append(toFloat(line.split(",")))
        except EOFError:
            break
    return vectors


def toFloat(arr):
    for i,str in enumerate(arr):
        arr[i] = float(str)
    return arr


def good_input(vectors, k):
    dim = len(vectors[0])
    for vector in vectors:
        if len(vector) != dim:
            print(f"INPUT ERROR:\nthe vector {vector} is {len(vector)} dimentional, different from the first vector dimention: {dim}")
            return False
    if len(vectors) < k:
        print(f"INPUT ERROR:\nthere are less then k={k} data points")
        return False
    return True


def average(cluster):
    dim = len(cluster[0])
    avg = [0 for i in range(dim)]
    cluster_size = len(cluster)
    for point in cluster:
        avg = vectors_sum(point, avg)
    avg = scalar_div(avg, cluster_size)
    return avg


def distance(vec1, vec2):
    # if vec is scalar:
    if type(vec2) in (int, float):
        return abs(vec1-vec2)

    if len(vec1) != len(vec2):
        print("The vectors are not in the same length")
        return -1
    sum = 0
    for i in range(len(vec1)):
        sum += (vec1[i]-vec2[i])**2
    return sum


def vectors_sum(vec1,vec2):
    # if vec is scalar:
    if type(vec1) in (int, float):
        return vec1+vec2

    if len(vec1) != len(vec2):
        print("The vectors are not in the same length")
        return -1

    sum_vec = [0 for i in range(len(vec1))]
    for i in range(len(vec1)):
        sum_vec[i] = vec1[i] + vec2[i]
    return sum_vec


def scalar_div(vec, x):
    if x==0 :
        print("you can't divide by zero")
        return -1

    # if vec is scalar:
    if type(vec) in (int, float) :
        return vec/x

    for i in range(len(vec)):
        vec[i] = vec[i]/x
    return vec


def main():

    try:
        k = int(sys.argv[1])    
        if (k <= 0):
            print("INPUT ERROR:\nk can't be <= 0")
            return 1
    except ValueError:
        print("INPUT ERROR:\nk can't be a letter")
        return 1

    # algorithm
    if len(sys.argv) >= 3:
        try:
            max_iter = int(sys.argv[2])    

        except ValueError:
            print("INPUT ERROR:\nk can't be a letter")
            return 1
    
        if max_iter <= 0:
            print("INPUT ERROR:\nmax iteration can't be <= 0")
            return 1

        centroids = Kmeans(k, max_iter=max_iter)
    else:
        centroids = Kmeans(k)

    # print in the right format
    for j,centroid in enumerate(centroids):
        if j != 0:
            print("\n", end="")
        for i,coordinate in enumerate(centroid):
            if i!=0 :
                print(",", end="")
            print("%.4f" % coordinate, end="")


if __name__ == "__main__":
    main()
