all: executable

executable: convert_to_array.o cuda_mesh_estimator.o main.o cuda_path_estimator.o index.o cuda_mesh_weights.o one_dim_results.o cuda_mesh_generation.o 
	nvcc convert_to_array.o cuda_mesh_estimator.o main.o cuda_path_estimator.o index.o cuda_mesh_weights.o one_dim_results.o cuda_mesh_generation.o -arch=sm_35 -o executable 


main.o: main.cu
	nvcc -c main.cu -arch=sm_35

cuda_mesh_generation.o: cuda_mesh_generation.cu
	nvcc -c cuda_mesh_generation.cu -arch=sm_35 

convert_to_array.o: convert_to_array.cpp
	g++ -c convert_to_array.cpp

cuda_mesh_estimator.o: cuda_mesh_estimator.cu
	nvcc -c cuda_mesh_estimator.cu -arch=sm_35 

cuda_path_estimator.o: cuda_path_estimator.cu
	nvcc -c cuda_path_estimator.cu -arch=sm_35

index.o: index.cpp
	g++ -c index.cpp 

one_dim_results.o: one_dim_results.cpp
	g++ -c one_dim_results.cpp

cuda_mesh_weights.o: cuda_mesh_weights.cu
	nvcc -c cuda_mesh_weights.cu -arch=sm_35

clean:
	rm -rf *o mesh


