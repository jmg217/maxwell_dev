#include <cuda.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>


double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);


double* two_dim_index(double* vector, int i, int j, double m, int b);


__device__ double* two_dim_indexME(double* vector, int i, int j, double m, int b){


double* p;

//specify index layout here
p=&vector[b*(i)+(j)];

return p;


}

__device__ double* three_dim_indexME(double* matrix, int i, int j, int k, double m, int b, int num_assets){

double* p;

p=&matrix[i*b*num_assets+j*num_assets+k];
return p;

}

//this function returns the payoff of a geometric average call option.
__device__ double GeometricPayOffCallM(double* X, int i, int j, double m, int b, int num_assets, double Strike){
double h;
h=1;
for(int l=0; l<num_assets; l++){

                h*= exp(*three_dim_indexME(X, i, j, l, m, b, num_assets));
        }
        h=pow(h,1.0/(num_assets));

        if(h-Strike>0){
                h=h-Strike;
        }
        else{
                h=0;
        }
return h;
}

//this function returns the payoff of a geometric average put option.
__device__ double GeometricPayOffPutM(double* X, int i, int j, double m, int b, int num_assets, double Strike){
double h;
h=1;
for(int l=0; l<num_assets; l++){

                h*=exp(*three_dim_indexME(X, i, j, l, m, b, num_assets));
        }
        h=pow(h,1.0/(num_assets));

        if(Strike-h>0){
                h=Strike-h;
        }
        else{
                h=0;
        }
return h;
}

//this function returns the continuation value given by the inner control.
__device__ double inner_control_meshME(int i, int j, int b, double r, double delta_t, double m, double* W_device, double* X_device, double* V_device, int num_assets ){

int m_int= (int)m;
double ControlMeshContinVal=0;

double sum=0, ContinVal=0, stock_expectation=0, true_stock_expectation=0, numerator=0, denominator=0, beta=0;

for(int k=0; k<b; k++){


        sum+= (*three_dim_indexME(W_device, (m_int-i), k, j, m, b, b)) * (*two_dim_indexME(V_device, (i-1), k, m, b));
}

ContinVal=(1/((double)b))*sum;

sum=0;

for(int l=0; l<b; l++){


        sum+=(*three_dim_indexME(W_device, (m_int-i), l, j, m, b, b)) * exp((*three_dim_indexME(X_device, (m_int-i), l,0 , m, b, num_assets))); 

}
stock_expectation=(1/((double)b))*sum; 

true_stock_expectation=exp(*three_dim_indexME(X_device, (m_int-i-1), j, 0, m, b, num_assets)) * exp(r*delta_t);//old


for(int p=0; p<b; p++){


numerator += ( (*three_dim_indexME(W_device, (m_int-i), p, j, m, b, b)) * exp( *three_dim_indexME(X_device, (m_int-i), p, 0, m, b, num_assets) ) - stock_expectation ) * ( (*three_dim_indexME(W_device, (m_int-i), p, j, m, b, b)) * (*two_dim_indexME(V_device, (i-1), p, m, b)) - ContinVal );


denominator += pow( ( (*three_dim_indexME(W_device, (m_int-i), p, j, m, b, b))  * exp( *three_dim_indexME(X_device, (m_int-i), p, 0, m, b, num_assets) ) - stock_expectation ) , 2 );

}

beta=numerator/denominator;

ControlMeshContinVal= ContinVal-beta*(stock_expectation-true_stock_expectation);

return ControlMeshContinVal;
}



//this kernel updates the V matrix of high bias mesh prices
__global__ void MeshEstimatorKernel(double strike, double r, double delta_t, int b, double m, double* X_device, double* W_device, double* V_device, double* asset_amount_device, int num_assets, int ker){


double H; //payoff variable 
double C; //continuation value variable
int idx =blockDim.x*blockIdx.x + threadIdx.x;

if(idx<b){


	if(ker==0){
			

		H=GeometricPayOffCallM( X_device, m-1-ker, idx, m, b, num_assets, strike)*exp(-r*delta_t*(m-ker));

		*two_dim_indexME(V_device, ker, idx, m, b)=H;
	}
	
	else{


		C=inner_control_meshME(ker, idx, b, r, delta_t, m, W_device, X_device, V_device, num_assets);



		H=GeometricPayOffCallM( X_device, m-1-ker, idx, m, b, num_assets, strike)*exp(-r*delta_t*(m-ker));	
		if(H>=C){

			*two_dim_indexME(V_device, ker, idx, m, b)=H;
			
		}

		else{
			*two_dim_indexME(V_device, ker, idx, m, b)=C;
		}	
	}
	

}

}


//This function allocates memory on the device and reutrns the high bias estimate to the main function.
double MeshEstimator(double strike, double r, double delta_t, int b, double m, double* X, double* W, double* V, double asset_amount[], int num_assets){

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


dim3 gridDim((int)ceil(b/512.0));
dim3 blockDim(512.0);


cudaError_t error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


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
