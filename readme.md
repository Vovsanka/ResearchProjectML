# Research Project Machine Learning 2025
## Partially Optimal Cubic Subspace Clustering

C++ implementation of the partial optimality algorithm for cubic clique partitioning problem.

### Requires:
- Boost
- Eigen3

## Build and run: 
mkdir build
cd build
cmake ..
cmake --build .
./CubicClustering

- The seeds.txt contains 15 random seeds that can be used in the experiments (see main.cpp)
- The evaluation of the cost function (optional) is written to same.txt and diff.txt in the build folder;
- The generated planes and points are written to planes.csv and points.csv respectively in the build folder;
- The problem solving logs are written to logs.txt in the build folder;
- The computation time, partial optimality and accuracy are written to data.csv (configurable) in the build folder;

There are some Python scripts (require numpy, pandas, matplotlib) that can process the log files 
mentioned before.
- run plotting.py to plot the cost evaluation and the generated planes and points afterwards;
- run process.py to process data.csv and get the computation time, partial optimality and accuracy for the 1., 2., 3. quartiles (after sorting).
Modify the script output as needed!
(15 experiments: q1=4, q2=8, q3=12) 
(7 experiments: q1=2, q2=4, q3=6) 

