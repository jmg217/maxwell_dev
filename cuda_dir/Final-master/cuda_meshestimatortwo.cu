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
//double inner_control_mesh(int i, int j, int b, double r, double delta_t, double m, std::vector< std::vector< std::vector<double> > >& W, std::vector<std::vector< std::vector<double> > >& X, std::vector< std::vector<double> >& V);
double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);

//__host__ __device__ double* two_dim_index(double* vector, int i, int j, double m, int b);
double* two_dim_index(double* vector, int i, int j, double m, int b);


__device__ double* two_dim_indexME(double* vector, int i, int j, double m, int b){

//int m_int= (int)m;
double* p;

//specify index layout here
p=&vector[b*(i)+(j)];

return p;


}

__device__ double* three_dim_indexME(double* matrix, int i, int j, int k, double m, int b, int num_assets){

int m_int = (int)m;
double* p;

//specify index layout here
//p=&matrix[(m_int)*b*(k)+(m_int)*(j)+(i)];
p=&matrix[i*b*num_assets+j*num_assets+k];
return p;

}


__device__ double GeometricPayOffCallM(double* X, int i, int j, double m, int b, int num_assets, double Strike){
double h;
h=1;
for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);
                h*= exp(*three_dim_indexME(X, i, j, l, m, b, num_assets));
        }
        h=pow(h,1.0/(num_assets));
//h=std::max(h-Strike,0.0);
        if(h-Strike>0){
                h=h-Strike;
        }
        else{
                h=0;
        }
return h;
}

__device__ double GeometricPayOffPutM(double* X, int i, int j, double m, int b, int num_assets, double Strike){
double h;
h=1;
for(int l=0; l<num_assets; l++){
        //h*=exp(X[i][j][l]);
                h*=exp(*three_dim_indexME(X, i, j, l, m, b, num_assets));
        }
        h=pow(h,1.0/(num_assets));
//h=std::max(Strike-h,0.0);
        if(Strike-h>0){
                h=Strike-h;
        }
        else{
                h=0;
        }
return h;
}


__device__ double inner_control_meshME(int i, int j, int b, double r, double delta_t, double m, double* W_device, double* X_device, double* V_device, int num_assets ){

int m_int= (int)m;
double ControlMeshContinVal=0;
//double outersum=0;//new
//double true_stock_expectation_sum=0;//new
double sum=0, ContinVal=0, stock_expectation=0, true_stock_expectation=0, numerator=0, denominator=0, beta=0;

for(int k=0; k<b; k++){

        //sum+=(W[(m-i)][k][j])*V[i-1][k]; //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]                         
        sum+= (*three_dim_indexME(W_device, (m_int-i), k, j, m, b, b)) * (*two_dim_indexME(V_device, (i-1), k, m, b));
}

ContinVal=(1/((double)b))*sum;

sum=0;
//true_stock_expectation_sum=0;
//outersum=0;
//for(int ll=0; ll< num_assets; ll++){  //new 
//sum=0;
for(int l=0; l<b; l++){

        //sum+=(W[(m-i)][l][j])*exp(X[m-i][l][0]); //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]                         
        sum+=(*three_dim_indexME(W_device, (m_int-i), l, j, m, b, b)) * exp((*three_dim_indexME(X_device, (m_int-i), l,0 , m, b, num_assets))); //old
        //sum+=(*three_dim_index(W, (m_int-i), l, j, m, b, b)) * exp((*three_dim_index(X, (m_int-i), l,ll , m, b, num_assets)));//new

}
stock_expectation=(1/((double)b))*sum; //old
//stock_expectation=(1/((double)num_assets))*outersum;//new

//true_stock_expectation=exp(X[m-i-1][j][0])*exp(r*delta_t);
true_stock_expectation=exp(*three_dim_indexME(X_device, (m_int-i-1), j, 0, m, b, num_assets)) * exp(r*delta_t);//old
//true_stock_expectation=(1/((double)num_assets))*true_stock_expectation_sum;

for(int p=0; p<b; p++){

//numerator += ( W[(m-i)][p][j] * exp( X[m-i][p][0] ) - stock_expectation ) * ( W[(m-i)][p][j] * V[i-1][p] - ContinVal );
numerator += ( (*three_dim_indexME(W_device, (m_int-i), p, j, m, b, b)) * exp( *three_dim_indexME(X_device, (m_int-i), p, 0, m, b, num_assets) ) - stock_expectation ) * ( (*three_dim_indexME(W_device, (m_int-i), p, j, m, b, b)) * (*two_dim_indexME(V_device, (i-1), p, m, b)) - ContinVal );

//denominator += pow( ( W[(m-i)][p][j]  * exp( X[m-i][p][0] ) - stock_expectation ) , 2 );
denominator += pow( ( (*three_dim_indexME(W_device, (m_int-i), p, j, m, b, b))  * exp( *three_dim_indexME(X_device, (m_int-i), p, 0, m, b, num_assets) ) - stock_expectation ) , 2 );

}

beta=numerator/denominator;

ControlMeshContinVal= ContinVal-beta*(stock_expectation-true_stock_expectation);

return ControlMeshContinVal;
}


/*
//this function returns the payoff value
double payoff(std::vector<std::vector< std::vector<double> > >& X, double k, std::vector<double>& asset_amount, int i, int j){
double h;
//h=k-x;
//h=0;
h=1;
//for(int l=0; l<asset_amount.size(); l++){
//	h+=asset_amount[l]*exp(X[i][j][l]);
//}
for(int l=0; l<asset_amount.size(); l++){
	h*=exp(X[i][j][l]);
}
h=pow(h,1.0/(asset_amount.size()));
h=h-k;
	if(h<0){
	h=0;
	}

	return h;
	}
*/
//this function returns the high bias mesh price
__global__ void MeshEstimatorKernel(double strike, double r, double delta_t, int b, double m, double* X_device, double* W_device, double* V_device, double* asset_amount_device, int num_assets, int ker){
//Containers Matrix = matrix;
//Containers Vector = vector;



double H; //payoff variable 
double C; //continuation value variable
//double sum;

// temp vector in Estimator loop
//std::vector< double > tempvec;
int idx =blockDim.x*blockIdx.x + threadIdx.x;

if(idx<b){

//Mesh Estimator loop

//for(int i=0; i<m; i++){
//tempvec.clear();

	
		if(ker==0){
			
		//H=payoff(X, strike, asset_amount, m-1-i, j)*exp(-r*delta_t*(m-i));	
		//H=thePayOff(X, asset_amount_device, m-1-i,idx, m, b, num_assets, Strike)*exp(-r*delta_t*(m-i));
		H=GeometricPayOffCallM( X_device, m-1-ker, idx, m, b, num_assets, strike)*exp(-r*delta_t*(m-ker));
		//tempvec.push_back(H);
		//V[m*(j)+(i)]=H;
		*two_dim_indexME(V_device, ker, idx, m, b)=H;
		}
	
		else{
//continue to develope continuation vale			

//NEW CONTROL VARIATE VERSION OF CONTIN VAL

		C=inner_control_meshME(ker, idx, b, r, delta_t, m, W_device, X_device, V_device, num_assets);


////////////PREVIOUS CALCULATION of CONTIN VAL. 
	/*	sum=0;
			for(int k=0; k<b; k++){
//std::cout<< sum<<std::endl;			
			//sum+=(W[(m-i)][k][j])*V[i-1][k]; //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]	
			sum+=(*two_dim_index(V, i-1, k, m, b)) * (*three_dim_index(W,(m-i),  k, j, m, b, b));	
			}
	
		C=(1/((double)b))*sum; //continuation value
	*/		
	
////////////END OF CALCULATION OF PREVIOUS CONTIN VAL. (WITHOUT CONTROL VARIATES)
	
		/*if(m-i==2 && j == 50){
		std::cout<<"contin value="<< C<<std::endl; 
		}*/
		//H=payoff(X, strike, asset_amount, m-1-i, j)*exp(-r*delta_t*(m-i));
		//H=thePayOff(X_device, asset_amount_device, m-1-i,idx, m, b,  Matrix, num_assets_device)*exp(-r*delta_t*(m-i));
		H=GeometricPayOffCallM( X_device, m-1-ker, idx, m, b, num_assets, strike)*exp(-r*delta_t*(m-ker));	
			if(H>=C){
				//tempvec.push_back(H);
				*two_dim_indexME(V_device, ker, idx, m, b)=H;
			
			}

			else{
				//tempvec.push_back(C);
				*two_dim_indexME(V_device, ker, idx, m, b)=C;
			}	
		}
	
//V.push_back(tempvec); //high bias option value matrix
//}

}
//sum over option values at the second time step

}

double MeshEstimator(double strike, double r, double delta_t, int b, double m, double* X, double* W, double* V, double asset_amount[], int num_assets){

double V_0;
int m_int=(int)m;


double* asset_amount_host;
asset_amount_host =asset_amount;


//double* asset_amount_host;
//asset_amount_host =asset_amount;

int X_N=(m_int) * b * (num_assets);
int asset_amount_N=num_assets;
int W_N=(m_int) * b*b;
int V_N=(m_int) * b;
//int weight_denominator_N=(N-1) * b;

double* X_device;
double* V_device;
double* asset_amount_device;
double* W_device;


cudaMalloc((void**) &X_device, X_N*sizeof(double) );
cudaMemcpy(X_device, X, X_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &asset_amount_device, asset_amount_N*sizeof(double) );
cudaMemcpy(asset_amount_device, asset_amount, asset_amount_N*sizeof(double), cudaMemcpyHostToDevice);

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
//cudaDeviceSynchronize();
if(ker<m-1){
cudaMemcpy(V_device, V, V_N*sizeof(double), cudaMemcpyHostToDevice);
}
//cudaDeviceSynchronize();
}

cudaFree(X_device);
cudaFree(asset_amount_device);
cudaFree(V_device);
cudaFree(W_device);

double sum=0;

for(int k=0; k<b; k++){
//sum+=V[m-1][k];
sum+=*two_dim_index(V, (m_int-1), k, m, b);
}
//this is the high bias option value at time 0
V_0=(1/((double)b))*sum;

//Calculate estimate european value
/*sum=0;
for(int hh=0; hh<b; hh++){
//sum+=V[m-1][k];
sum+=thePayOff(X, asset_amount, m-1, hh, m, b,  Matrix, num_assets)*exp(-r*delta_t*(m));
}
double euroval=(1/((double)b))*sum;
EuroVals.push_back(euroval);
*/
return V_0;

}
