import numpy as np
import pandas as pd
import mykmeanssp
import math
import sys
np.random.seed(0)

def eigengapHuristic(eigenvalues):
    n = len(eigenvalues)
    if (n <= 1):
      return 0
    eigenvalues.sort()
    deltas = [abs(eigenvalues[i]-eigenvalues[i+1]) for i in range(len(eigenvalues)-1)]
    relevant_deltas = np.array(deltas[:math.floor(n/2)])
    return np.argmax(relevant_deltas)

def draw_initial_centroids(K, dataPoints):
    pass


def format_print(centroids):
    for centroid in centroids:
        c = ""
        for cord in range(len(centroid)-1):
            x = "%.4f"%centroid[cord]
            c += x +','
        c += "%.4f"%centroid[cord+1]
        print(c)

def main():
    print("main")
    args = sys.argv[1:]
    if (len(args) == 3):
        K = args[0]
        goal = args[1]
        input_data = args[2]
    else:
        goal = args[0]
        input_data = args[1]
        
    dataPoints = np.loadtxt(input_data,delimiter=",",dtype=float)
    dataPoints = dataPoints.tolist()
    if (goal == "wam"):
        res = mykmeanssp.wam(dataPoints, len(dataPoints))
    elif (goal == "ddg"):
        res = mykmeanssp.ddg(dataPoints, len(dataPoints))
    elif (goal == "gl"):
        res = mykmeanssp.gl(dataPoints, len(dataPoints))
    elif (goal == "jacobi"){
        res = mykmeanssp.jacobi(dataPoints, len(dataPoints))
    }
    elif (goal == "spk"):
        gl = mykmeanssp.gl(dataPoints,len(dataPoints))
        jacobi = mykmeanssp.jacobi(gl, len(gl))
        if (len(args) == 2):
            K = eigengapHuristic(jacobi[0])
        spectral_kmeans_vectors = np.array(jacobi[1:K+1])
        spectral_kmeans_vectors = np.transpose(spectral_kmeans_vectors)
        chosen_centroids, initial_centroids = draw_initial_centroids(K, spectral_kmeans_vectors)
        #res = spk(initial_centroids, kmeans_vectors, len(kmeans_vectors), len(kmeans_vectors[0]))
    format_print(res)


if __name__ == '__main__':
    main()
