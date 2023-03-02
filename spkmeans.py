import numpy as np
import pandas as pd
import mykmeanssp
import sys
np.random.seed(0)

def eigengapHuristic():
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
        K = eigengapHuristic()
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
    format_print(res)


if __name__ == '__main__':
    main()
