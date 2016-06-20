/* This code will multiply two vectors and
   check the result.
*/

#include <cuda.h>
#include <iostream>

/* Fill in your dotProduct kernel here...
 */
#define THREADS_PER_BLOCK 256
__device__ float result;

__global__ void calcDotProductKern(float *x, float *y, int N)
{
  __shared__ float product[THREADS_PER_BLOCK]; //All threads in a block must be able 
                                               //to access this array

    int index = threadIdx.x + blockIdx.x * blockDim.x; //index

    product[threadIdx.x] = x[index] * y[index]; //result of elementwise
                                                //multiplication goes into product

    //Make sure every thread has finished
    __syncthreads();

    //Sum the elements serially to obtain dot product
    if( 0 == threadIdx.x ) //Pick one thread to sum, otherwise all will execute
    {
        float sum = 0;
        for(int j=0; j < THREADS_PER_BLOCK; j++) 
		sum += product[j];
        atomicAdd(&result,sum);
    }
}

float calcDotProduct(float *x, float *y, int N)
{
  int threads = THREADS_PER_BLOCK;
  int blocks = (N + threads - 1)/ threads;

  float result = 0;
  cudaMemcpyToSymbol(result, &result, sizeof(result), 0, cudaMemcpyHostToDevice);

  calcDotProductKern<<<blocks,threads>>>(x, y, N);

  cudaMemcpyFromSymbol(&result, result, sizeof(result), 0, cudaMemcpyDeviceToHost);

  return result;
}


int main(void)
{
  const int N = 10000;

  float *x_host = new float[N];
  float *y_host = new float[N];
 

 

  // Fill matrix and vector on host
  for(int i=0 ; i < N ; i++)
  {
    x_host[i] = sin(i*0.1f);
    y_host[i] = cos(i*0.23f);
  }
  
  float *x;
  float *y;
 

  cudaMalloc(&x, N*sizeof(float));
  cudaMalloc(&y, N*sizeof(float));
 

  // Copy x and y to device
  cudaMemcpy(x, x_host, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(y, y_host, N*sizeof(float), cudaMemcpyHostToDevice);  

  cudaEvent_t start, end;
  cudaEventCreate(&start);
  cudaEventCreate(&end);

  cudaEventRecord(start, 0);

  float prodGPU = calcDotProduct(x, y, N);

  cudaEventRecord(end, 0);
  cudaEventSynchronize(end);

  cudaMemcpy(y_host, y, N*sizeof(float), cudaMemcpyDeviceToHost);

  // Check result
  float prod = 0;
  for(int i=0 ; i < N ; i++)
  {
    prod += y_host[i] * x_host[i];
  }
  

  if( fabs(prod - prodGPU) / prod < 1e-4 )
  {
    std::cout << "Multiplication correct!" << std::endl;
    float timeInMs;
    cudaEventElapsedTime(&timeInMs, start, end);
    
    std::cout << "Time: " << timeInMs << "ms" << std::endl;
    std::cout << "Bandwidth: " << (2*N*sizeof(float)) / 1.0e9 / (timeInMs/1000) << "Gbps" << std::endl;
  }
  else
  {
    std::cout << "Multiplication wrong!" << std::endl;
  }
  
  cudaFree(x);
  cudaFree(y);

  delete[] x_host;
  delete[] y_host;
}
