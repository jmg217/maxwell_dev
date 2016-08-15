#include <unistd.h>
#include <stdio.h>

 /* we need these includes for CUDA's random number stuff */
#include <curand.h>
#include <curand_kernel.h>

__global__ void init(unsigned int seed, curandState_t* states) {

   /* we have to initialize the state */
   curand_init(seed, /* the seed can be the same for each core, here we pass the time in from the CPU */
               threadIdx.x, /* the sequence number should be different for each core (unless you want all
 +                             cores to get the same sequence of numbers for some reason - use thread id! */
               0, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
               &states[threadIdx.x]);
 }


void cuda_init_random(int N, curandState_t* states){

//curandState_t* states;
cudaMalloc((void**) &states, N * sizeof(curandState_t));

   /* invoke the GPU to initialize all of the random states */
   init<<<1, N>>>(time(0), states);

}
