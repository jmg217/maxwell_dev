#include <unistd.h>
#include <stdio.h>
 
 /* we need these includes for CUDA's random number stuff */
 #include <curand.h>
 #include <curand_kernel.h>
 
 #define N 5
 
 #define MAX 100
 
 /* this GPU kernel function is used to initialize the random states */
// __global__ void init(unsigned int seed, curandState_t* states) {
 
   /* we have to initialize the state */
//   curand_init(seed, /* the seed can be the same for each core, here we pass the time in from the CPU */
//               threadIdx.x, /* the sequence number should be different for each core (unless you want all
// +                             cores to get the same sequence of numbers for some reason - use thread id! */
//               0, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
//               &states[threadIdx.x]);
// }
 
 /* this GPU kernel takes an array of states, and an array of ints, and puts a random int into each */
 __global__ void randoms(curandState_t* states, double* numbers) {
   /* curand works like rand - except that it takes a state as a parameter */
   numbers[threadIdx.x] = curand_normal_double(&states[threadIdx.x]);
 }
 
void cppfunc(int it, int final, curandState_t* States, curandState_t* states) {
   /* CUDA's random number library uses curandState_t to keep track of the seed value
 +     we will store a random state for every thread  */
//curandState_t* states;
/*
curandState_t* states;
//states=&pointer;
if(it==0){  
   //curandState_t* states;
//   allocate space on the GPU for the random states 
   cudaMalloc((void**) &states, N * sizeof(curandState_t));
 //    invoke the GPU to initialize all of the random states 
   init<<<1, N>>>(time(0), states);
} 
*/  
 /* allocate an array of unsigned ints on the CPU and GPU */
   double cpu_nums[N];
   double* gpu_nums;
   cudaMalloc((void**) &gpu_nums, N * sizeof(double));
// for(int it=0; it<10;it++){
printf("NEW BLOCK!!!!!!!!!!!!\n");
   /* invoke the kernel to get some random numbers */
 cudaMemcpy(states, States, N * sizeof(double), cudaMemcpyHostToDevice);
   randoms<<<1, N>>>(states, gpu_nums);
 cudaDeviceSynchronize();
cudaMemcpy(States, states, sizeof(curandState_t)*N, cudaMemcpyDeviceToHost);
   /* copy the random numbers back */
   cudaMemcpy(cpu_nums, gpu_nums, N * sizeof(double), cudaMemcpyDeviceToHost);
 
   /* print them out */
   for (int i = 0; i < N; i++) {
      printf("%lf\n", cpu_nums[i]);
   }
 //}
   /* free the memory we allocated for the states and numbers */
cudaFree(gpu_nums);

if(it==final){
   cudaFree(states);
 } // cudaFree(gpu_nums);
//cudaDeviceReset(); 
//return states;
   
 }
