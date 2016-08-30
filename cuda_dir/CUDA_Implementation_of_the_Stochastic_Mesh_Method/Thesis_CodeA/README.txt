Code Description:
The files main.cu, cuda_mesh_generation.cu, convert_to_array.cpp, cuda_mesh_estimator.cu, cuda_path_estimator.cu, index.cpp, one_dim_results.cpp, cuda_mesh_weights.cu collectively implement the stochastic mesh method for pricing multidimensional Bermudan options in a CUDA framework. The input parameters are specified in settings.txt. The cuda_mesh_generation.cu file generates the nodes of the stochastic mesh and the file cuda_mesh_weights.cu generates the mesh weights. The file cuda_mesh_estimator.cu generates the high bias option value and the file cuda_path_estimator.cu generates the low bias option value. The file main.cu contains the main function which reads in the input parameters, initialises the random seeds on the GPU threads, and performs the final calculation of confidence intervals.

Parameters:

Compilation:

Execution:

Output:

Matlab Code:



