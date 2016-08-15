#include <cuda.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <curand.h>
#include <curand_kernel.h>
#define PI 3.14159265358979323846


__device__ double densityMW(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t){
double f=0, x=0;
//x=(1/(sigma*sqrt(delta_t)))*(log(Xnew)-log(Xold)-(r-delta-0.5*sigma*sigma)*delta_t);
x=(1/(sigma*sqrt(delta_t)))*(Xnew-Xold-(r-delta-0.5*sigma*sigma)*delta_t);
//f= (1/(sigma*sqrt(delta_t)*Xnew))*(1/(sqrt(2*PI)))*exp(-0.5*x*x); // this is the transition density
f= (1/(sigma*sqrt(delta_t)))*(1/(sqrt(2*PI)))*exp(-0.5*x*x);
return f;

}



__device__ double* three_dim_indexMW(double* matrix, int i, int j, int k, double m, int b, int num_assets){

int m_int = (int)m;
double* p;

//specify index layout here
//p=&matrix[(m_int)*b*(k)+(m_int)*(j)+(i)];
p=&matrix[i*b*num_assets+j*num_assets+k];
return p;

}


__device__ double* two_dim_indexMW(double* vector, int i, int j, double m, int b){

//int m_int= (int)m;
double* p;

//specify index layout here
p=&vector[b*(i)+(j)];

return p;


}


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


__global__ void meshweightsKernel(double* W_device, double  m, int b, double* sigma_device, double* delta_device, double r, double delta_t, double* X_device, int num_assets, double* weight_denominator_device){
double wdenominator, w;
int N=(int)m;

int idx =blockDim.x*blockIdx.x + threadIdx.x;

if(idx<N){

//printf("m=%f",m);
double* sortvector;
sortvector=new double[b];

	
	if(idx==0){
		for(int k=0; k<b; k++){
        	
			for(int j=0; j<b; j++){

				if(j==0){
					
					*three_dim_indexMW(W_device, idx, k, j, m, b, b)=1;
				}// all weights from the starting node are equal to 1

				else{
					
					*three_dim_indexMW(W_device, idx, k, j, m, b, b)=0;
				}
			}

		}
	}


	if(idx>0){

		for(int k=0; k<b; k++){	
		//dim1temp.clear();
		//sortvector.clear();

		wdenominator=0;

			for(int j=0; j<b; j++){
			//std::cout<<j<<std::endl;	
				w=1;
				//w=0; //set w to 1 since it will be equal to a product
				for(int jj=0; jj<num_assets; jj++){
				w = w * densityMW(*three_dim_indexMW(X_device, (idx-1), j, jj, m, b, num_assets), *three_dim_indexMW(X_device, idx, k, jj, m, b, num_assets), sigma_device[jj], r, delta_device[jj], delta_t);

				}
			//w = exp(w);
			//dim1temp.push_back(w);
			//sortvector.push_back(w);   
			sortvector[j]=w;                                                                
			}
			//std::sort(sortvector.begin(), sortvector.end());
			wdenominator=kahansum(sortvector, b);
				//devide each element by the denominator
//			std::cout<<"before"<<std::endl;
			*two_dim_indexMW(weight_denominator_device, idx-1, k, m-1, b)=wdenominator;
//			std::cout<<"after and I= "<<I<<std::endl;
//			std::cout<<*two_dim_index(weight_denominator, I, k, m-1, b)<<std::endl;
			//printf("wdenom at idx=%i and k=%i is=%f\n",idx,k,wdenominator);
			for(int t=0; t<b; t++){
				*three_dim_indexMW(W_device, (idx), k, t, m, b, b)=(((double)b)*(sortvector[t]))/wdenominator;
			}
			//std::cout<<"after"<<std::endl;
		}
	//std::cout<<"outside 1"<<std::endl;	
	
	}

//std::cout<<"outside 2"<<std::endl;

//std::cout<<"outside 3"<<std::endl;
delete[] sortvector;
}//std::cout<<"end"<<std::endl;
}



void meshweights(double* W, double m, int b, double sigma[], double delta[], double r, double delta_t, double* X, int num_assets, double* weight_denominator, int threads){

int m_int=(int)m;


double* sigma_host;
sigma_host =sigma;

double* delta_host;
delta_host =delta;

//double* asset_amount_host;
//asset_amount_host =asset_amount;

int X_N=(m_int) * b * (num_assets);
int delta_N= num_assets;
int sigma_N=num_assets;
int W_N=(m_int) * b * b;
int w_N=(m_int-1) * b ;
//int weight_denominator_N=(N-1) * b;

double* X_device;
double* sigma_device;
double* delta_device;
double* W_device;
double* weight_denominator_device;

cudaMalloc((void**) &X_device, X_N*sizeof(double) );
cudaMemcpy(X_device, X, X_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &sigma_device, sigma_N*sizeof(double) );
cudaMemcpy(sigma_device, sigma_host, sigma_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &delta_device, delta_N*sizeof(double) );
cudaMemcpy(delta_device, delta_host, delta_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &W_device, W_N*sizeof(double) );
cudaMemcpy(W_device, W, W_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &weight_denominator_device, w_N*sizeof(double) );
cudaMemcpy(weight_denominator_device, weight_denominator, w_N*sizeof(double), cudaMemcpyHostToDevice);


cudaError_t error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


//printf("here");

dim3 gridDim((int)ceil((m_int)/512.0));
dim3 blockDim(512.0);

meshweightsKernel<<<gridDim, blockDim>>>(W_device, m, b, sigma_device, delta_device, r, delta_t,  X_device, num_assets, weight_denominator_device);


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
cudaFree(weight_denominator_device);
cudaFree(delta_device);
cudaFree(W_device);

//for(int test=0; test<W_N; test++){
//printf("weights=%f\n",W[test]);
//}

error = cudaGetLastError();
if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }

}
