#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#define PI 3.14159265358979323846

__device__ double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t);


__host__ __device__ double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);


__host__ __device__ double* two_dim_index(double* vector, int i, int j, double m, int b);


__device__ double kahansum(double* sortvector, int b){
double sum=0, c=0, y, t;

        for(int i=0; i<b; i++){
                y=sortvector[i]-c;
                t=sum+y;
                c=(t-sum)-y;
                sum=t;
        }

return sum;
}


__global__ void valuesKernel(double* tempW_device ,double m,int  b, double* sigma_device,double* delta_device,double r, double delta_t,double* X_device,int num_assets){


int idx =blockDim.x*blockIdx.x + threadIdx.x;


int m_int=(int)m;
if(idx<(m_int-1)*b*b){

double w;

int i=idx/(b*b);
int j=idx/b;
if(j>(b-1)){
j=j%b;
}
int k=idx%b;


w=1;
for(int jjj=0; jjj<num_assets; jjj++){
	w = w * density(*three_dim_index(X_device, (i), k, jjj, m, b, num_assets), *three_dim_index(X_device, i+1, j, jjj, m, b, num_assets), sigma_device[jjj], r, delta_device[jjj], delta_t);
        }


tempW_device[idx]=w;

}
}


__global__ void sumweightsKernel(double* tempW_device , int b, double* weight_denominator_device, double m){

int idx =blockDim.x*blockIdx.x + threadIdx.x;
int m_int=(int)m;
if(idx<(m_int-1)*b){

double sum=0, c=0, y, t;

int start=idx*b;

for(int i=start; i<start+b; i++){
                y=tempW_device[i]-c;
                t=sum+y;
                c=(t-sum)-y;
                sum=t;
}

weight_denominator_device[idx]=sum;
}
}


__global__ void meshweightsKernel(double* W_device, double m, int b, double* sigma_device, double* delta_device, double r, double delta_t, double* X_device, int num_assets, double* weight_denominator_device, double* tempW_device){
double wdenominator;

int idx =blockDim.x*blockIdx.x + threadIdx.x;
int m_int=(int)m;

if(idx<b*b*m_int){

int i=idx/(b*b);
int k=idx/b;
if(k>(b-1)){
k=k%b;
}
int j=idx%b;


if(i==0){

                                if(j==0){

                                        *three_dim_index(W_device, i, k, j, m, b, b)=1;
                                }// all weights from the starting node are equal to 1

                                else{

                                        *three_dim_index(W_device, i, k, j, m, b, b)=0;
                                }
}	
if(i>0){		//sortvector.clear();


			wdenominator= *two_dim_index(weight_denominator_device, i-1, k, m-1, b);
			*three_dim_index(W_device, (i), k, j, m, b, b)=(((double)b) * (*three_dim_index(tempW_device, i-1, k, j, m-1, b, b)))/wdenominator;
		}

}
}

void meshweights(double* W, double m, int b, double sigma[], double delta[], double r, double delta_t, double* X, int num_assets, double* weight_denominator){
int m_int=(int)m;

int temp_N=(m_int-1) * b*b;

double* sigma_host;
sigma_host =sigma;
double* delta_host;
delta_host=delta;
double* tempW;
tempW= new double[temp_N];


int X_N=(m_int) * b * (num_assets);
int W_N=(m_int) * b*b;
int w_N=(m_int-1)*b;
int sigma_N=num_assets;
int delta_N=num_assets;


double* X_device;
double* W_device;
double* weight_denominator_device;
double* sigma_device;
double* delta_device;
double* tempW_device;

cudaMalloc((void**) &X_device, X_N*sizeof(double) );
cudaMemcpy(X_device, X, X_N*sizeof(double), cudaMemcpyHostToDevice);

cudaError_t error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


cudaMalloc((void**) &W_device, W_N*sizeof(double) );
cudaMemcpy(W_device, W, W_N*sizeof(double), cudaMemcpyHostToDevice);

error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


cudaMalloc((void**) &sigma_device, sigma_N*sizeof(double) );
cudaMemcpy(sigma_device, sigma_host, sigma_N*sizeof(double), cudaMemcpyHostToDevice);

error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


cudaMalloc((void**) &delta_device, delta_N*sizeof(double) );
cudaMemcpy(delta_device, delta_host, delta_N*sizeof(double), cudaMemcpyHostToDevice);

error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


cudaMalloc((void**) &weight_denominator_device, w_N*sizeof(double) );
cudaMemcpy(weight_denominator_device, weight_denominator, w_N*sizeof(double), cudaMemcpyHostToDevice);

error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


cudaMalloc((void**) &tempW_device, temp_N*sizeof(double) );
cudaMemcpy(tempW_device, tempW, temp_N*sizeof(double), cudaMemcpyHostToDevice);


error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }

cudaDeviceSetLimit(cudaLimitMallocHeapSize, 80000000*sizeof(double));

dim3 VgridDim((int)ceil(temp_N/512.0));
dim3 VblockDim(512.0);

valuesKernel<<<VgridDim,VblockDim>>>(tempW_device , m, b, sigma_device, delta_device, r, delta_t, X_device, num_assets);
cudaDeviceSynchronize();

cudaMemcpy(tempW, tempW_device, sizeof(double)*temp_N, cudaMemcpyDeviceToHost);
cudaMemcpy(tempW_device, tempW, temp_N*sizeof(double), cudaMemcpyHostToDevice);

error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


dim3 sgridDim((int)ceil(w_N/512.0));
dim3 sblockDim(512.0);

error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


sumweightsKernel<<<sgridDim, sblockDim>>>(tempW_device , b, weight_denominator_device, m);

cudaDeviceSynchronize();

cudaMemcpy(weight_denominator, weight_denominator_device, sizeof(double)*w_N, cudaMemcpyDeviceToHost);
cudaMemcpy(weight_denominator_device, weight_denominator, w_N*sizeof(double), cudaMemcpyHostToDevice);


error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }

dim3 mgridDim((int)ceil(W_N/512.0));
dim3 mblockDim(512.0);
meshweightsKernel<<<mgridDim, mblockDim>>>(W_device , m, b, sigma_device, delta_device, r, delta_t, X_device, num_assets, weight_denominator_device, tempW_device);

cudaDeviceSynchronize();



 error = cudaGetLastError();


  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }

cudaMemcpy(W, W_device, sizeof(double)*W_N, cudaMemcpyDeviceToHost);
cudaMemcpy(weight_denominator, weight_denominator_device, sizeof(double)*w_N, cudaMemcpyDeviceToHost);


cudaFree(X_device);
cudaFree(sigma_device);
cudaFree(delta_device);
cudaFree(W_device);
cudaFree(weight_denominator_device);
cudaFree(tempW_device);

delete[] tempW;
}


