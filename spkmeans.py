import numpy as np
import pandas as pd
import mykmeanssp
import sys
np.random.seed(0)

def eigengapHuristic():
    pass


def main():
    args =sys.argv[1:]
    goalValues = ["spk", "wam", "ddg", "gl", "jacobi"]
    length = len(args)
    if (length < 2 or length > 3 ):
        return False
    elif (length == 2):
        k = eigengapHuristic()
        goal = args[0]
        file_name = args[1]
    else:
        k = args[0]
        goal = args[1]
        file_name = args[2]
    


    N = getNumberOfDataPoints()