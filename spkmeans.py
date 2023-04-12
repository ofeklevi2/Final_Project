import numpy as np
import pandas as pd
import mykmeanssp
import math
import sys
np.random.seed(0)

def eigengapHuristic(eigenvalues):
    n = len(eigenvalues)
    deltas = [abs(eigenvalues[i]-eigenvalues[i+1]) for i in range(len(eigenvalues)-1)]
    relevant_deltas = np.array(deltas[:math.floor(n/2)])
    return np.argmax(relevant_deltas) + 1

def draw_new_index(data,centroids,indexes):
    def D(vec,centroids):
        distances = [np.linalg.norm(vec - centroids[i]) for i in range(len(centroids))]
        return np.min(distances)
    
    D_arr = [D(vec, centroids) for vec in data]
    D_sum = np.sum(D_arr)
    probabilities = [D_arr[i]/D_sum for i in range(len(data))]
    return np.random.choice(indexes, p=probabilities)

def draw_initial_centroids(K,dataPoints):
  indexes = range(len(dataPoints))
  selected_indexes = []
  centroids = []
  new_index = np.random.choice(indexes)
  selected_indexes.append(new_index)
  centroids.append(dataPoints[new_index].tolist())
  while len(centroids) < K:
    new_index = draw_new_index(dataPoints, centroids, indexes)
    selected_indexes.append(new_index)
    centroids.append(dataPoints[new_index].tolist())
  return selected_indexes,centroids

def format_print(elements, isInts = False):
    if (isInts):
        c=""
        for elem in elements:
            c += str(elem)+","
        c = c[:-1]
        print(c)
        return 

    for elem in elements:
        c = ""
        for cord in range(len(elem)-1):
            x = "%.4f"%elem[cord]
            c += x +','
        c += "%.4f"%elem[cord+1]
        print(c)

def main():
    args = sys.argv[1:]
    if (len(args) == 3):
        K = int(args[0])
        goal = args[1]
        input_data = args[2]
    else:
        goal = args[0]
        input_data = args[1]
        
    dataPoints = np.loadtxt(input_data,delimiter=",",dtype=float) 
    dataPoints = dataPoints.tolist()
    if (goal == "wam"):
        print(len(dataPoints), len(dataPoints[0]))
        res = mykmeanssp.wam(dataPoints, len(dataPoints),len(dataPoints[0]))
    elif (goal == "ddg"):
        res = mykmeanssp.ddg(dataPoints, len(dataPoints),len(dataPoints[0]))
    elif (goal == "gl"):
        res = mykmeanssp.gl(dataPoints, len(dataPoints),len(dataPoints[0]))
    elif (goal == "jacobi"):
        res = mykmeanssp.jacobi(dataPoints, len(dataPoints),0)
    elif (goal == "spk"):
        gl = mykmeanssp.gl(dataPoints,len(dataPoints), len(dataPoints[0]))
        jacobi = mykmeanssp.jacobi(gl, len(gl),1)
        print(jacobi)
        if (len(args) == 2):
            K = eigengapHuristic(jacobi[0])
        spectral_kmeans_vectors = [lst[:K] for lst in jacobi]
        spectral_kmeans_vectors = spectral_kmeans_vectors [1:]
        spectral_kmeans_vectors = np.array(spectral_kmeans_vectors)
        selected_indexes, initial_centroids = draw_initial_centroids(K, spectral_kmeans_vectors)
        spectral_kmeans_vectors = spectral_kmeans_vectors.tolist()
        for centroid in initial_centroids:
            spectral_kmeans_vectors.remove(centroid)
        dataPoints = initial_centroids + spectral_kmeans_vectors
        format_print(selected_indexes, isInts =True)
        res = mykmeanssp.spk(K, dataPoints ,len(dataPoints), len(dataPoints[0]))
    format_print(res)


if __name__ == '__main__':
    main()
