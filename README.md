# CVRPFeatureExtractors
A large collection of implemented feature extractors for capacitated vehicle routing problems

# Usage

Configure and run `main.py` using Python 2.7 and with all of the dependencies installed.

# Input and Output

Reads [TSPLIB formatted](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/) formatted `.vrp` files and outputs `.csv`files with around 400 columns (features), with each line describing a single CVRP instance.

# Configuring

* Change the folders path in `main.py/main()` to point to the actual folders the `.vrp` files reside in
* Change the `ACOTSP_EXE_PATH` in `tsp_solver.py` to point to the ACOTSP executable
* Compile the VRPH binaries and SYMPHONY binaries and place them into `solvers`folder

# Dependencies

* Python 2.7
* Numpy - for vector and matrix operations
* sklearn - for min-max scaling, PCA, clustering, and MDS
* scipy - for convenient distance matrix calculation and statistics
* networkx - for calculating digraph nearest neighbor features
* matplotlib and PIL (OPTIONAL) - for plotting clusters (analyze_clusters.py only)
* VRPH (init routine)

# TODO

* Add support for XML CVRP instances (use helpers.cvrp_rkm16_io/read_TSPLIB_CVRP() as a template)
* Add a command line user interface to main.py, e.g., by using argparse
* Write unittests to verify correct operation of the extractors  