#include <cuda.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include "Payoff.h"
#include "enum_header.h"

__host__ __device__ double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);

__host__ __device__ double* two_dim_index(double* vector, int i, int j, double m, int b);

__device__ double GeometricPayOffCallM(double* X, int i, int j, double m, int b, int num_assets, double Strike);

__device__ double GeometricPayOffPutM(double* X, int i, int j, double m, int b, int num_assets, double Strike);

__device__ double inner_control_meshME(int i, int j, int b, double r, double delta_t, double m, double* W_device, double* X_device, double* V_device, int num_assets );

//this function returns the high bias mesh price
__global__ void MeshEstimatorKernel(double strike, double r, double delta_t, int b, double m, double* X_device, double* W_device, double* V_device, double* asset_amount_device, int num_assets, int ker){

double H; //payoff variable 
double C; //continuation value variable
int idx =blockDim.x*blockIdx.x + threadIdx.x;

if(idx<b){

	if(ker==0){
			
		H=GeometricPayOffCallM( X_device, m-1-ker, idx, m, b, num_assets, strike)*exp(-r*delta_t*(m-ker));
	
		*two_dim_index(V_device, ker, idx, m, b)=H;
	}
	
	else{
		
		//the inner control function calculate the continuation value using a control variate
		C=inner_control_meshME(ker, idx, b, r, delta_t, m, W_device, X_device, V_device, num_assets);


		H=GeometricPayOffCallM( X_device, m-1-ker, idx, m, b, num_assets, strike)*exp(-r*delta_t*(m-ker));	
		
		if(H>=C){
			*two_dim_index(V_device, ker, idx, m, b)=H;
			
		}

		else{
			*two_dim_index(V_device, ker, idx, m, b)=C;
		}	
	}
	
}
//end of if statement
}

//this function allocates the gpu memory copies data to the gpu
double MeshEstimator(double strike, double r, double delta_t, int b, double m, double* X, double* W, double* V, double asset_amount[], int num_assets ){

double V_0;
int m_int=(int)m;


double* asset_amount_host;
asset_amount_host =asset_amount;


int X_N=(m_int) * b * (num_assets);
int asset_amount_N=num_assets;
int W_N=(m_int) * b*b;
int V_N=(m_int) * b;
double* X_device;
double* V_device;
double* asset_amount_device;
double* W_device;


cudaMalloc((void**) &X_device, X_N*sizeof(double) );
cudaMemcpy(X_device, X, X_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &asset_amount_device, asset_amount_N*sizeof(double) );
cudaMemcpy(asset_amount_device, asset_amount_host, asset_amount_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &V_device, V_N*sizeof(double) );
cudaMemcpy(V_device, V, V_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &W_device, W_N*sizeof(double) );
cudaMemcpy(W_device, W, W_N*sizeof(double), cudaMemcpyHostToDevice);

//set the number of threads
dim3 gridDim((int)ceil(b/512.0));
dim3 blockDim(512.0);


cudaError_t error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }

//we loop of the high bias kernel for each time step
for(int ker=0; ker<m; ker++){
MeshEstimatorKernel<<<gridDim, blockDim>>>(strike, r, delta_t, b, m,  X_device,  W_device, V_device, asset_amount_device, num_assets, ker);

cudaDeviceSynchronize();


error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


cudaMemcpy(V, V_device, sizeof(double)*V_N, cudaMemcpyDeviceToHost);
if(ker<m-1){
cudaMemcpy(V_device, V, V_N*sizeof(double), cudaMemcpyHostToDevice);
}
}

cudaFree(X_device);
cudaFree(asset_amount_device);
cudaFree(V_device);
cudaFree(W_device);

double sum=0;

for(int k=0; k<b; k++){
sum+=*two_dim_index(V, (m_int-1), k, m, b);
}

//this is the high bias option value at time 0
V_0=(1/((double)b))*sum;


return V_0;

}
