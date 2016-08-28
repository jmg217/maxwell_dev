all: executable

executable: convert_to_array.o cuda_meshestimatortwo.o main.o cudapathestimator.o index.o cuda_meshweights.o Payoff.o plot_high_bias.o cuda_mesh_generation.o density.o payoff.o innerMeshControl.o
	nvcc convert_to_array.o cuda_meshestimatortwo.o main.o cudapathestimator.o index.o cuda_meshweights.o Payoff.o plot_high_bias.o cuda_mesh_generation.o density.o payoff.o innerMeshControl.o -arch=sm_35 -o executable 

main.o: main.cu
	nvcc -dc main.cu -arch=sm_35

cuda_mesh_generation.o: cuda_mesh_generation.cu
	nvcc -dc cuda_mesh_generation.cu -arch=sm_35 

convert_to_array.o: convert_to_array.cpp
	g++ -c convert_to_array.cpp

cuda_meshestimatortwo.o: cuda_meshestimatortwo.cu
	nvcc -dc cuda_meshestimatortwo.cu -arch=sm_35 

cudapathestimator.o: cudapathestimator.cu
	nvcc -dc cudapathestimator.cu -arch=sm_35

index.o: index.cu
	nvcc -dc index.cu -arch=sm_35

plot_high_bias.o: plot_high_bias.cpp
	g++ -c plot_high_bias.cpp

cuda_meshweights.o: cuda_meshweights.cu
	nvcc -dc cuda_meshweights.cu -arch=sm_35

density.o: density.cu
	nvcc -dc density.cu -arch=sm_35

payoff.o: payoff.cu
	nvcc -dc payoff.cu -arch=sm_35

innerMeshControl.o: innerMeshControl.cu
	nvcc -dc innerMeshControl.cu -arch=sm_35

clean:
	rm -rf *o mesh

