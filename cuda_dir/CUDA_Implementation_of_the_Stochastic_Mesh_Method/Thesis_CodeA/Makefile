all: executable

executable: convert_to_array.o cuda_meshestimatortwo.o main.o cudapathestimator.o index.o cuda_meshweights.o plot_high_bias.o cuda_mesh_generation.o 
	nvcc convert_to_array.o cuda_meshestimatortwo.o main.o cudapathestimator.o index.o cuda_meshweights.o plot_high_bias.o cuda_mesh_generation.o -arch=sm_35 -o executable 


main.o: main.cu
	nvcc -c main.cu -arch=sm_35

cuda_mesh_generation.o: cuda_mesh_generation.cu
	nvcc -c cuda_mesh_generation.cu -arch=sm_35 

convert_to_array.o: convert_to_array.cpp
	g++ -c convert_to_array.cpp

cuda_meshestimatortwo.o: cuda_meshestimatortwo.cu
	nvcc -c cuda_meshestimatortwo.cu -arch=sm_35 

cudapathestimator.o: cudapathestimator.cu
	nvcc -c cudapathestimator.cu -arch=sm_35

index.o: index.cpp
	g++ -c index.cpp 

plot_high_bias.o: plot_high_bias.cpp
	g++ -c plot_high_bias.cpp

cuda_meshweights.o: cuda_meshweights.cu
	nvcc -c cuda_meshweights.cu -arch=sm_35

clean:
	rm -rf *o mesh

