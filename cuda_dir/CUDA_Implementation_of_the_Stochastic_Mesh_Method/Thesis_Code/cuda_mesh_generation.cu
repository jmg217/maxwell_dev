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


__host__ __device__ double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);


__host__ __device__ double* two_dim_index(double* vector, int i, int j, double m, int b);


//this kernel generates the stochastic mesh.
__global__ void MeshGenKernel(double* X_device, double* delta_device, double* sigma_device,double* X0_device, int N, double strike, double r, double delta_t, int b, double m, int num_assets, curandState_t* states, double* asset_amount_device){
//thread index
int idx =blockDim.x*blockIdx.x + threadIdx.x;

if(idx<N){
int m_int=(int)m;
double Xi, Z;

//we set the values of each mesh node in this loop
for(int i=0; i<m_int; i++){
        if(i==0){
        
                for(int ll=0; ll<num_assets; ll++){
                       
			Z=curand_normal_double(&states[idx]);
			
        		*three_dim_index(X_device, i, idx, ll, m, b, num_assets) = X0_device[ll] +  (r-delta_device[ll]-0.5*pow(sigma_device[ll], 2))*delta_t + sigma_device[ll]*sqrt(delta_t)*Z;

                }
        
        }

        if(i>0){
                for(int jj=0; jj<num_assets; jj++){
                        
			Z=curand_normal_double(&states[idx]);
                        Xi=*three_dim_index(X_device, (i-1), idx, jj, m, b, num_assets);
                        *three_dim_index(X_device, i, idx, jj, m, b, num_assets)=Xi +  (r-delta_device[jj]-0.5*pow(sigma_device[jj], 2))*delta_t + sigma_device[jj]*sqrt(delta_t)*Z;
                }
        
        }
}

}

}
 
void mesh_generation(int b, int num_assets, double m, double X0[], double sigma[], double delta[], double asset_amount[], double* X, double strike, double r, double delta_t, curandState_t* States, curandState_t* states, int threads){

int m_int=(int)m;

int N= b;

double* X0_host;
X0_host =X0;

double* sigma_host;
sigma_host =sigma;

double* delta_host;
delta_host =delta;

double* asset_amount_host;
asset_amount_host =asset_amount;

int X_N=(m_int) * b * (num_assets);
int delta_N= num_assets;
int sigma_N=num_assets;
int X0_N=num_assets;
int asset_amount_N = num_assets;

double* X_device;
double* sigma_device;
double* delta_device;
double* X0_device;
double* asset_amount_device;

//allocate memory on the device and copy from cpu to device

cudaMalloc((void**) &X_device, X_N*sizeof(double) );
cudaMemcpy(X_device, X, X_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &X0_device, X0_N*sizeof(double) );
cudaMemcpy(X0_device, X0_host, X0_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &sigma_device, sigma_N*sizeof(double) );
cudaMemcpy(sigma_device, sigma_host, sigma_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &delta_device, delta_N*sizeof(double) );
cudaMemcpy(delta_device, delta_host, delta_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &asset_amount_device, asset_amount_N*sizeof(double) );
cudaMemcpy(asset_amount_device, asset_amount_host, asset_amount_N*sizeof(double), cudaMemcpyHostToDevice);


cudaMemcpy(states, States, threads*sizeof(curandState_t*), cudaMemcpyHostToDevice);

//set the number of threads
dim3 gridDim((int)ceil(N/512.0));
dim3 blockDim(512.0);

//launch kernel
MeshGenKernel<<<gridDim, blockDim>>>(X_device, delta_device, sigma_device, X0_device, N, strike, r, delta_t, b,  m, num_assets, states, asset_amount_device);

cudaDeviceSynchronize();

cudaMemcpy(X, X_device, sizeof(double)*X_N, cudaMemcpyDeviceToHost);
cudaError_t error = cudaGetLastError(); 
 if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("Found at line %d\n", __LINE__);
    exit(1);
  }
cudaMemcpy(States, states, sizeof(curandState_t)*threads, cudaMemcpyDeviceToHost);

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("Found at line %d\n", __LINE__);
    exit(1);
  }

//free the gpu memory

cudaFree(X_device);
cudaFree(sigma_device);
cudaFree(delta_device);
cudaFree(X0_device);
cudaFree(asset_amount_device);


}

