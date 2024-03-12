# Spectral K-Means Clustering Implementation

This project implements the unnormalized spectral k-means clustering algorithm, which uses a hybrid approach integrating C for computational efficiency and Python for ease of interface and data handling.

## Files Overview

- `makefile`: Contains the compilation instructions for the C components of the project, ensuring that the code adheres to strict standards (`-ansi`, `-Wall`, etc.) and linking necessary libraries (e.g., `-lm` for math).

- `spkmeans.c`: The core C file that contains functions for performing weighted adjacency matrix calculation, diagonal degree matrix generation, graph Laplacian, and Jacobi transformations.

- `spkmeans.h`: The header file for `spkmeans.c`, which includes the declarations of the functions and the data structures used throughout the C code.

- `input.c`: An example file that demonstrates how to structure the input data and utilize the functions defined in `spkmeans.c` to convert vector data to arrays and print them.

- `setup.py`: A Python setup file used for creating the Python extension module, allowing the C functions to be called from Python scripts.

- `spkmeans.py`: The main Python script that utilizes the `mykmeanssp` module to perform spectral clustering. It includes functions for the eigengap heuristic, initializing centroids, and executing the k-means algorithm.

- `spkmeansmodule.c`: Defines the methods of the Python extension module, which allows the translation of Python objects to C data types and vice versa, enabling the seamless execution of the underlying C algorithms.



