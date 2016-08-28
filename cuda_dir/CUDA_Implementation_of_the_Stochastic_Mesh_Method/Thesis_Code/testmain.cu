#include <unistd.h>
#include <stdio.h>

 /* we need these includes for CUDA's random number stuff */
 #include <curand.h>
 #include <curand_kernel.h>
#define N 5

 #define MAX 100

__global__ void init(unsigned int seed, curandState_t* states) {

   /* we have to initialize the state */
   curand_init(seed, /* the seed can be the same for each core, here we pass the time in from the CPU */
               threadIdx.x, /* the sequence number should be different for each core (unless you want all
 +                             cores to get the same sequence of numbers for some reason - use thread id! */
               0, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
               &states[threadIdx.x]);
 }


void cuda_init_random(curandState_t* States, curandState_t* states){

//curandState_t* states;
//curandState_t* states;
//cudaMalloc((void**) &states, N * sizeof(curandState_t));

   /* invoke the GPU to initialize all of the random states */
   init<<<1, N>>>(time(0), states);
cudaDeviceSynchronize();
cudaMemcpy(States, states, sizeof(curandState_t)*N, cudaMemcpyDeviceToHost);

}

void cppfunc(int it, int final, curandState_t* States, curandState_t* states);
//void cppfunc(int it, int final, curandState_t* states);

void hello_world(int it);

//void cuda_init_random(int N, curandState_t* states);

int main(){

curandState_t* States;
States= new curandState_t[N];
//curandState_t* states;
//curandState_t* states;
//cudaMalloc((void**) &states, N * sizeof(curandState_t));
//cuda_init_random(States, states);
int it;
for( it=0; it<10; it++){

curandState_t* states;
//curandState_t* states;
cudaMalloc((void**) &states, N * sizeof(curandState_t));
if(it==0){cuda_init_random(States, states);}
else{cudaMemcpy(states, States, N*sizeof(curandState_t), cudaMemcpyHostToDevice);}

cppfunc(it, 10, States, states);

cudaDeviceReset();
}
hello_world(it);

delete[] States;
return 0;

}
