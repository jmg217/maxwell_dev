/* This code will multiply two vectors and
   check the result.
*/

#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <iostream>
#include <stdio.h>

#define CUDA_CHECK {cudaThreadSynchronize();	\
  cudaError_t err = cudaGetLastError();\
  if(err){\
    std::cout << "Error: " << cudaGetErrorString(err) << " line " << __LINE__ << std::endl; \
    exit(1);\
  }}

//#define CUDA_CHECK

/* Fill in your dotProduct kernel here...
 */

// Naive (and wrong) way
__global__ void calcDotProductKern1(float *x, float *y, float *res, int N)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i < N){
    (*res) += x[i] * y[i];
  }
}

// Naive, correct, but slow way
__global__ void calcDotProductKern2(float *x, float *y, float *res, int N)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i < N){
    atomicAdd(res, x[i] * y[i]);
  }
}

// Better reduction first pass
// res needs to point to at least 'blocks' floats.
__global__ void calcDotProductKern3(float *x, float *y, float *res, int N)
{
  __shared__ float product[512];
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int locI = threadIdx.x;
  if(i < N)
    {
      product[locI] = x[i] * y[i];
    }
  else{
    product[locI] = 0;
  }
  __syncthreads();

  int blockSize = blockDim.x;

  if( blockSize >= 1024 && locI < 512){
    product[locI] = product[locI] + product[locI+512];
    __syncthreads();
  }
  if( blockSize >= 512 && locI < 256){
    product[locI] = product[locI] + product[locI+256];
    __syncthreads();
  }
  if( blockSize >= 256 && locI < 128){
    product[locI] = product[locI] + product[locI+128];
    __syncthreads();
  }
  if( blockSize >= 128 && locI < 64){
    product[locI] = product[locI] + product[locI+64];
    __syncthreads();
  }
  if( blockSize >= 64 && locI < 32){
    product[locI] = product[locI] + product[locI+32];
    __syncthreads();
  }
  if( blockSize >= 32 && locI < 16){
    product[locI] = product[locI] + product[locI+16];
    __syncthreads();
  }
  if( blockSize >= 16 && locI < 8){
    product[locI] = product[locI] + product[locI+8];
  }
  if( blockSize >= 8 && locI < 4){
    product[locI] = product[locI] + product[locI+4];
  }
  if( blockSize >= 4 && locI < 2){
    product[locI] = product[locI] + product[locI+2];
  }
  if( blockSize >= 2 && locI < 1){
    product[locI] = product[locI] + product[locI+1];
  }
  if( locI == 0){
    res[blockIdx.x] = product[0];
  }
}

// Generic reduction
// x[] is of size N, and y[] is of size N/blockDim.x
__global__ void reduce(float *x, float *y, int N)
{
  __shared__ float result[512];
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int locI = threadIdx.x;
  int blockSize = blockDim.x;

  if(i < N){
    result[locI] = x[i];
  }else{
    result[locI] = 0;
  }
  __syncthreads();
  if( blockSize >= 1024 && locI < 512){
    result[locI] = result[locI] + result[locI+512];
    __syncthreads();
  }
  if( blockSize >= 512 && locI < 256){
    result[locI] = result[locI] + result[locI+256];
    __syncthreads();
  }
  if( blockSize >= 256 && locI < 128){
    result[locI] = result[locI] + result[locI+128];
    __syncthreads();
  }
  if( blockSize >= 128 && locI < 64){
    result[locI] = result[locI] + result[locI+64];
    __syncthreads();
  }
  if( blockSize >= 64 && locI < 32){
    result[locI] = result[locI] + result[locI+32];
    __syncthreads();
  }
  if( blockSize >= 32 && locI < 16){
    result[locI] = result[locI] + result[locI+16];
    __syncthreads();
  }
  if( blockSize >= 16 && locI < 8){
    result[locI] = result[locI] + result[locI+8];
  }
  if( blockSize >= 8 && locI < 4){
    result[locI] = result[locI] + result[locI+4];
  }
  if( blockSize >= 4 && locI < 2){
    result[locI] = result[locI] + result[locI+2];
  }
  if( blockSize >= 2 && locI < 1){
    result[locI] = result[locI] + result[locI+1];
  }
  if(locI == 0){
    y[blockIdx.x] = result[0];
  }
}

float calcDotProduct1(float* x, float* y, int N){

  int threads = 512;
  int blocks = (N + threads - 1)/ threads;

  float* res;
  cudaMalloc(&res, sizeof(float));
  float resHost = 0;
  cudaMemcpy(res, &resHost, sizeof(float), cudaMemcpyHostToDevice);
  CUDA_CHECK;
  calcDotProductKern1<<<blocks, threads>>>(x, y, res, N);
  CUDA_CHECK;
  cudaMemcpy(&resHost, res, sizeof(float), cudaMemcpyDeviceToHost);
  cudaFree(res);
  CUDA_CHECK;
  return resHost;
}

float calcDotProduct2(float* x, float* y, int N){

  int threads = 512;
  int blocks = (N + threads - 1)/ threads;

  float* res;
  cudaMalloc(&res, sizeof(float));
  float resHost = 0;
  cudaMemcpy(res, &resHost, sizeof(float), cudaMemcpyHostToDevice);
  CUDA_CHECK;
  calcDotProductKern2<<<blocks, threads>>>(x, y, res, N);
  CUDA_CHECK;
  cudaMemcpy(&resHost, res, sizeof(float), cudaMemcpyDeviceToHost);
  cudaFree(res);
  CUDA_CHECK;
  return resHost;
}

float calcDotProduct3(float* x, float* y, int N){
  int threads = 512;
  int blocks = (N + threads - 1)/ threads;

  float* res;
  cudaMalloc(&res, blocks*sizeof(float));
  CUDA_CHECK;
  calcDotProductKern3<<<blocks, threads>>>(x, y, res, N);
  CUDA_CHECK;
  float* resHost = new float[blocks];
  cudaMemcpy(resHost, res, sizeof(float) * blocks, cudaMemcpyDeviceToHost);
  CUDA_CHECK;
  float p=0;
  for(int i=0 ; i < blocks ; i++){
    p += resHost[i];
  }

  delete[] resHost;
  cudaFree(res);
  CUDA_CHECK;
  return p;
}


float calcDotProduct3Reduce(float* x, float* y, int N){
  int threads = 512;
  int blocks = (N + threads - 1)/ threads;

  float* res;
  cudaMalloc(&res, blocks*sizeof(float));
  CUDA_CHECK;
  calcDotProductKern3<<<blocks, threads>>>(x, y, res, N);
  CUDA_CHECK;
  while(blocks > 1){
    int blocksOrig = blocks;
    blocks = ceil((float)blocks / threads);
    float* resOrig = res;
    cudaMalloc(&res, blocks*sizeof(float));
    reduce<<<blocks, threads>>>(resOrig, res, blocksOrig);

    cudaFree(resOrig);
  }
  CUDA_CHECK;
  float resHost;
  cudaMemcpy(&resHost, res, sizeof(float), cudaMemcpyDeviceToHost);
  CUDA_CHECK;
  cudaFree(res);
  return resHost;
}

float calcDotProductThrust(float* x, float* y, int N){
  thrust::device_ptr<float> xThStart(x);
  thrust::device_ptr<float> yThStart(y);
  thrust::device_ptr<float> xThEnd(x + N);
  thrust::device_ptr<float> yThEnd(y + N);

  return thrust::inner_product(xThStart, xThEnd, yThStart, 0.0f);
}

float timeDotProduct(float (*kernel)(float*, float*, int), float *x, float *y, int N, float ans)
{
  CUDA_CHECK;
  cudaEvent_t start, end;
  cudaEventCreate(&start);
  cudaEventCreate(&end);
  CUDA_CHECK;
  cudaEventRecord(start, 0);
  
  double prod = kernel(x, y, N);

  cudaEventRecord(end, 0);
  CUDA_CHECK;
  cudaEventSynchronize(end);
  CUDA_CHECK;
  cudaError_t err = cudaGetLastError();

  if(err){
    std::cout << "Error: " << cudaGetErrorString(err) << std::endl;
  }

  if( fabs(prod - ans) / fabs(ans) < 1e-4 )
  {
    std::cout << "Multiplication correct! " << prod << " = " << ans << std::endl;
  }
  else
  {
    std::cout << "Multiplication wrong! " << prod << " != " << ans << std::endl;
  }

  float timeInMs;
  cudaEventElapsedTime(&timeInMs, start, end);
  std::cout << "Time: " << timeInMs << "ms" << std::endl << std::endl;

  CUDA_CHECK;

  cudaEventDestroy(start);
  cudaEventDestroy(end);

  CUDA_CHECK;
  
  return 0;
}


int main(void)
{
  const int N = 20000000;

  float *x_host = new float[N];
  float *y_host = new float[N];
  

  // Fill matrix and vector on host
  for(int i=0 ; i < N ; i++)
  {
    x_host[i] = sin(i*0.013);
    y_host[i] = cos(i*0.019);
  }
  
  float *x;
  float *y;

  cudaMalloc(&x, N*sizeof(float));
  cudaMalloc(&y, N*sizeof(float));

  CUDA_CHECK;
  
  // Copy x and y to device
  cudaMemcpy(x, x_host, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(y, y_host, N*sizeof(float), cudaMemcpyHostToDevice);  

  CUDA_CHECK;

  //cudaMemcpy(y_host, y, N*sizeof(float), cudaMemcpyDeviceToHost);

  // Check result

  clock_t st = clock();
  float prod = 0;
  for(int i=0 ; i < N ; i++)
  {
    prod += y_host[i] * x_host[i];
  }
  clock_t end = clock();

  std::cout << "CPU time = " << (end - st) / (float)CLOCKS_PER_SEC * 1000 << " ms" << std::endl;

  std::cout << "Naive approach - wrong" << std::endl;
  timeDotProduct(calcDotProduct1, x, y, N, prod);
  std::cout << "Using atomic operations" << std::endl;
  timeDotProduct(calcDotProduct2, x, y, N, prod);
  std::cout << "Reduction across one thread block only" << std::endl;
  timeDotProduct(calcDotProduct3, x, y, N, prod);
  std::cout << "Repeated reduction" << std::endl;
  timeDotProduct(calcDotProduct3Reduce, x, y, N, prod);
  std::cout << "Thrust" << std::endl;
  timeDotProduct(calcDotProductThrust, x, y, N, prod);

  cudaFree(x);
  cudaFree(y);

  delete[] x_host;
  delete[] y_host;
}
