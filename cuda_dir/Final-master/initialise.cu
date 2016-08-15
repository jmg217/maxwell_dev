#include <cuda.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include "enum_header.h"
#include <unistd.h>
#include <stdio.h>

/* we need these includes for CUDA's random number stuff */
#include <curand.h>
#include <curand_kernel.h>

__global__ void init(unsigned int seed, curandState_t* states) {
int idx=blockDim.x*blockIdx.x + threadIdx.x;

  /* we have to initialize the state */
  curand_init(seed, /* the seed can be the same for each core, here we pass the time in from the CPU */
              idx, /* the sequence number should be different for each core (unless you want all
                             cores to get the same sequence of numbers for some reason - use thread id! */
              0, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
              &states[idx]);
}

void cuda_init_random(int N){

dim3 gridDim((int)ceil(N/512.0));
dim3 blockDim(512.0);

curandState_t* states;

cudaMalloc((void**) &states, N * sizeof(curandState_t));

init<<<gridDim, blockDim>>>(time(0), states);

cudaDeviceSynchronize();

cudaError_t error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }

}
