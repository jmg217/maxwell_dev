#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <iostream>
#include <stdio.h>


//NOTE: COMPILE WITH -arch=sm_20

#define CUDA_CHECK {cudaThreadSynchronize();	\
  cudaError_t err = cudaGetLastError();\
  if(err){\
    std::cout << "Error: " << cudaGetErrorString(err) << " line " << __LINE__ << std::endl; \
    exit(1);\
  }}


double calcDotProductThrust(double* x, double* y, int N){
  thrust::device_ptr<double> xThStart(x);
  thrust::device_ptr<double> yThStart(y);
  thrust::device_ptr<double> xThEnd(x + N);
  thrust::device_ptr<double> yThEnd(y + N);

  return thrust::inner_product(xThStart, xThEnd, yThStart, 0.0);
}

double timeDotProduct(double (*kernel)(double*, double*, int), double *x, double *y, int N, double ans)
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


//int main(void)

double inner_product(int N, std::vector<double>& first_vector, std::vector<double>& second_vector)
{
  //const int N = 20000000;
/*
  float *x_host = new float[N];
  float *y_host = new float[N];
  */
  
  double *x_host = new double[N];
  double *y_host = new double[N];

  std::copy(first_vector.begin(), first_vector.end(), x_host);
  std::copy(second_vector.begin(), second_vector.end(), y_host);
  // Fill matrix and vector on host
  /*for(int i=0 ; i < N ; i++)
  {
   
    //x_host[i]=first_vector[i];
    //y_host[i]=second_vector[i];
    x_host[i] = sin(i*0.013);
    y_host[i] = cos(i*0.019);
  }*/


 /* 
  float *x;
  float *y;
*/

  double *x;
  double *y;

  cudaMalloc(&x, N*sizeof(double));
  cudaMalloc(&y, N*sizeof(double));

  CUDA_CHECK;
  
  // Copy x and y to device
  cudaMemcpy(x, x_host, N*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(y, y_host, N*sizeof(double), cudaMemcpyHostToDevice);  

  CUDA_CHECK;

  //cudaMemcpy(y_host, y, N*sizeof(float), cudaMemcpyDeviceToHost);

  // Check result

  clock_t st = clock();
  double prod = 0;
  for(int i=0 ; i < N ; i++)
  {
    prod += y_host[i] * x_host[i];
  }
  clock_t end = clock();

  std::cout << "CPU time = " << (end - st) / (float)CLOCKS_PER_SEC * 1000 << " ms" << std::endl;

/*  std::cout << "Naive approach - wrong" << std::endl;
  timeDotProduct(calcDotProduct1, x, y, N, prod);
  std::cout << "Using atomic operations" << std::endl;
  timeDotProduct(calcDotProduct2, x, y, N, prod);
  std::cout << "Reduction across one thread block only" << std::endl;
  timeDotProduct(calcDotProduct3, x, y, N, prod);
  std::cout << "Repeated reduction" << std::endl;
  timeDotProduct(calcDotProduct3Reduce, x, y, N, prod);
 */ 
  std::cout << "Thrust" << std::endl;
  timeDotProduct(calcDotProductThrust, x, y, N, prod);

  cudaFree(x);
  cudaFree(y);

  delete[] x_host;
  delete[] y_host;
}
