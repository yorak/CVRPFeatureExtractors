# CVRPFeatureExtractors
A large collection of implemented feature extractors for capacitated vehicle routing problems

![Clusters](results/clusters.png)

_A set of classical CVRP instances clustered using the computed features (Rasku et al. 2016)_

# Usage

[Configure](https://github.com/yorak/CVRPFeatureExtractors#configuring) and run `main.py` using Python 2.7 and with [all of the other dependencies installed](https://github.com/yorak/CVRPFeatureExtractors#dependencies).

# Input and Output

Reads [TSPLIB formatted](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/) formatted `.vrp` files and outputs `.csv`files with around 400 columns (features), with each line describing a single CVRP instance.

# Configuring

* Set the PYTHONPATH to point to the `CVRPFeatureExtractors` folder
* Change the `ACOTSP_EXE_PATH` in `tsp_solver.py` to point to the ACOTSP executable
* Compile the VRPH binaries and SYMPHONY binaries and place them into `solvers`folder
* Change the folders path in `main.py/main()` to point to the actual folders the `.vrp` files reside in

# Dependencies

* Python 2.7
* Numpy - for vector and matrix operations
* sklearn - for min-max scaling, PCA, clustering, and MDS
* scipy - for convenient distance matrix calculation and statistics
* networkx - for calculating digraph nearest neighbor features
* matplotlib and PIL (OPTIONAL) - for plotting clusters (analyze_clusters.py only)
* [custom VRPH](https://github.com/yorak/VRPH/tree/local_search_stats) init executable - used for local search probing
* SYMPHONY branch-and-cut solver - used for MIP probing
* [custom ACOTSP](https://github.com/juherask/ACOTSP) TSP solver - has an option to disable ant system, used to solve TSPs fast and accurately in estimation of constraint tightness

# TODO

* Add support for XML CVRP instances (use helpers.cvrp_rkm16_io/read_TSPLIB_CVRP() as a template)
* Add a command line user interface to main.py, e.g., by using argparse
* Try to find the ACOTSP from the `solvers` folder before resorting to hard coded path for the exe
* Write unit tests to verify correct operation of the extractors  

# Citing 

Rasku, J., Kärkkäinen, T., and Musliu, N. (2016). *Feature extractors for describing vehicle routing problem instances.* In OASIcs, SCOR'16, volume 50. Dagstuhl Publishing.