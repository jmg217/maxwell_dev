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

 
#define PI 3.14159265358979323846

__host__ __device__ double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);

__host__ __device__ double* two_dim_index(double* vector, int i, int j, double m, int b);

__device__ double GeometricPayOffCallV(double* X, double m, int b, int num_assets, double Strike);

__device__ double GeometricPayOffPutV(double* X, double m, int b, int num_assets, double Strike);

__device__ double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t);

//this function calculates the weights at each time step for each sub optimal stopping routine
__device__ void S_weights(double* S_Weights, double* X_device, double* S_new, int m, int b, double* sigma_device, double* delta_device, double delta_t, int num_assets, double r , int i, double* weight_denominator_device ){

double sum, w_s;

	for(int h=0; h<b; h++){   

	sum=0;
	w_s=1;
		for(int kk=0; kk<num_assets; kk++){
			w_s*=density(S_new[kk], *three_dim_index(X_device, (i+1), h, kk, m, b, num_assets), sigma_device[kk], r, delta_device[kk], delta_t);

		}
	sum = *two_dim_index(weight_denominator_device, i, h, m-1, b);
	if(sum==0){printf("division by zero in weights function of path estimator\n");}
	w_s = (((double)b)*w_s)/sum;	
	S_Weights[h]=w_s;
	}

}

//this kernel calculates the low bias price estimate
__global__ void PathEstimatorKernel(double* X_device, double* weight_denominator_device, double* V_device, double* delta_device, double* sigma_device, double* X0_device, int N, double strike, double r, double delta_t, int b, int m, int num_assets, curandState_t* states, double* results_dev, double* asset_amount_device){

//thread number
int idx =blockDim.x*blockIdx.x + threadIdx.x;
if(idx<N){

double v_0, S_i, Z, C, H, sum, weight; //, w_s, sum_Z;
const int S_N= num_assets;
const int S_W_N= b; 

//dynamic allocation of memory on the device heap memory.
double* S_new;
S_new= new double[S_N]; 
double* S_Weights;
S_Weights=new double[S_W_N];

int i=0;

//on each thread perform a do-while loop for the simulated basket of stock. the do_while loop exits when the option is exercised.
do {

	if(i==0){
		for(int ll=0; ll<num_assets; ll++){
			Z=curand_normal_double(&states[idx]);
			S_i=X0_device[ll] +  (r-delta_device[ll]-0.5*pow(sigma_device[ll], 2))*delta_t + sigma_device[ll]*sqrt(delta_t)*Z;
			S_new[ll]=S_i;			
		}
	}

	else{
		for(int jj=0; jj<num_assets; jj++){
			Z=curand_normal_double(&states[idx]);
			S_i=S_new[jj] + (r-delta_device[jj]-0.5*pow(sigma_device[jj], 2))*delta_t + sigma_device[jj]*sqrt(delta_t)*Z;
			S_new[jj]=S_i;
		}
	}

if(i<m-1){

S_weights(S_Weights, X_device, S_new, m, b, sigma_device, delta_device, delta_t, num_assets, r, i, weight_denominator_device);

}

double con_val=0; //continuation value variable
	sum=0;

	if(i==m-1){
	C=0;//continuation value at the last time step
	}
	
	else{
		for(int k=0; k<b; k++){	
			weight= S_Weights[k];
			
			con_val= *two_dim_index(V_device, (m-1-i-1), k, m, b);
			sum+=(weight) * (con_val); 			
		}
	
       
    
        C=(1/(double)b)*sum; //continuation value
	}	
	
H= GeometricPayOffCallV(S_new, m, num_assets, num_assets, strike)*exp(-r*delta_t*((i+1)));


i=i+1;
}while(H<C);//this will stop once H is less then the continuation value. at m-1, c=0 therefore m-1 is the max amount of loops.

v_0=H;

results_dev[idx]=v_0;


delete[] S_new;
delete[] S_Weights;
}



}

//this function allocates memory on the device for the path estimator.
double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma[], double delta[], double X0[], double* X, double* weight_denominator, double* V, double asset_amount[], int num_assets, int Path_estimator_iterations, int iterator, int Final_iteration, curandState_t* States, curandState_t* states, int threads ){


cudaError_t error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


int N= Path_estimator_iterations;

double* sigma_host;
sigma_host =sigma;

double* delta_host;
delta_host =delta;

double* X0_host;
X0_host =X0;

double* asset_amount_host;
asset_amount_host =asset_amount;

int m_int=(int)m;

int X_N=(m_int) * b * (num_assets);
int W_N=(m_int-1) * b;
int V_N=(m_int) * b;
int delta_N= num_assets;
int sigma_N=num_assets;
int X0_N=num_assets;
int asset_amount_N = num_assets;

double* X_device;
double* V_device;
double* weight_denominator_device;
double* sigma_device;
double* delta_device;
double* X0_device;
double* asset_amount_device;



error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


cudaMalloc((void**) &X_device, X_N*sizeof(double) );
cudaMemcpy(X_device, X, X_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &V_device, V_N*sizeof(double) );
cudaMemcpy(V_device, V, V_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &weight_denominator_device, W_N*sizeof(double) );
cudaMemcpy(weight_denominator_device, weight_denominator, W_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &X0_device, X0_N*sizeof(double) );
cudaMemcpy(X0_device, X0_host, X0_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &sigma_device, sigma_N*sizeof(double) );
cudaMemcpy(sigma_device, sigma_host, sigma_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &delta_device, delta_N*sizeof(double) );
cudaMemcpy(delta_device, delta_host, delta_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &asset_amount_device, asset_amount_N*sizeof(double) );
cudaMemcpy(asset_amount_device, asset_amount_host, asset_amount_N*sizeof(double), cudaMemcpyHostToDevice);


cudaMemcpy(states, States, threads*sizeof(curandState_t*), cudaMemcpyHostToDevice);

// set the number of threads
dim3 gridDim((int)ceil(N/512.0));
dim3 blockDim(512.0);


error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


double* results;
results = new double[N];
double* results_dev;
cudaMalloc((void**) &results_dev, N*sizeof(double) );

// CALL RANDOM SEEDING KERNEL HERE
error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


 
error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }

PathEstimatorKernel<<<gridDim, blockDim>>>(X_device, weight_denominator_device, V_device, delta_device, sigma_device, X0_device, N, strike, r, delta_t, b,  m_int, num_assets, states, results_dev, asset_amount_device);


cudaDeviceSynchronize();


error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
	printf("found at line %d\n", __LINE__);
    exit(1);
  }
cudaMemcpy(results, results_dev, sizeof(double)*N, cudaMemcpyDeviceToHost);

error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("Found at line %d\n", __LINE__);
    exit(1);
  }

cudaMemcpy(States, states, sizeof(curandState_t)*threads, cudaMemcpyDeviceToHost);

error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("Found at line %d\n", __LINE__);
    exit(1);
  }


double result=0;
for(int f=0; f<Path_estimator_iterations; f++){
result+=results[f];
}
result=(1/double(N))*result;

delete[] results;

error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


cudaFree(X_device);
cudaFree(V_device);
cudaFree(weight_denominator_device);
cudaFree(sigma_device);
cudaFree(delta_device);
cudaFree(X0_device);
cudaFree(results_dev);
cudaFree(asset_amount_device);

if(iterator==Final_iteration-1){
cudaFree(states);
}

return result;
}



