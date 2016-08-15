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

// REMEMBER TO PUT __host__ __device__ IN FRONT OF CLASS METHODS

#define PI 3.14159265358979323846
#define CHECK(x) do {\
cudaError_t err =(x);\
if (err !=cudaSuccess){\
	fprintf(stderr, "API error"\
	"%s:%d Returned:%d\n", \
	__FILE__, __LINE__, err);\
	exit(1);\
} while(0)


__device__ double* two_dim_indexGPU(double* vector, int i, int j, double m, int b){

//int m_int= (int)m;
double* p;

//specify index layout here
p=&vector[b*(i)+(j)];

return p;


}

__device__ double* three_dim_indexGPU(double* matrix, int i, int j, int k, double m, int b){

int m_int = (int)m;
double* p;

//specify index layout here
p=&matrix[(m_int)*b*(k)+(m_int)*(j)+(i)];

return p;

}

__device__ double densityGPU(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t){

double f=0, x=0;
//x=(1/(sigma*sqrt(delta_t)))*(log(Xnew)-log(Xold)-(r-delta-0.5*sigma*sigma)*delta_t);
x=(1/(sigma*sqrt(delta_t)))*(Xnew-Xold-(r-delta-0.5*sigma*sigma)*delta_t);
//f= (1/(sigma*sqrt(delta_t)*Xnew))*(1/(sqrt(2*PI)))*exp(-0.5*x*x); // this is the transition density
f= (1/(sigma*sqrt(delta_t)))*(1/(sqrt(2*PI)))*exp(-0.5*x*x);
return f;

}

__global__ void init(unsigned int seed, curandState_t* states) {
int idx=blockDim.x*blockIdx.x + threadIdx.x;

  /* we have to initialize the state */
  curand_init(seed, /* the seed can be the same for each core, here we pass the time in from the CPU */
              idx, /* the sequence number should be different for each core (unless you want all
                             cores to get the same sequence of numbers for some reason - use thread id! */
              0, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
              &states[idx]);
}


__device__ double GeometricPayOffCallV(double* X, int i, double m, int b, int num_assets, double Strike){
double h;
h=1;
for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);

                h*= exp(*two_dim_indexGPU(X, i, l, m, b));
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

__device__ double GeometricPayOffPutV(double* X, int i, double m, int b, int num_assets, double Strike){
double h;
h=1; 
for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);
                h*= exp(*two_dim_indexGPU(X, i, l, m, b));
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

__device__ void S_weights(double* S_Weights, double* X_device, double* S, double m, int b, double* sigma_device, double* delta_device, double delta_t, int num_assets, double r , int i ){

double density_product, sum, w_s;

	for(int h=0; h<b; h++){   //h=k
	sum=0;
	w_s=1;
		for(int kk=0; kk<num_assets; kk++){
			w_s*=densityGPU(*two_dim_indexGPU(S, i, kk, m, num_assets), *three_dim_indexGPU(X_device, (i+1), h, kk, m, b), sigma_device[kk], r, delta_device[kk], delta_t);
		}

	density_product=1;
	
		for(int g=0; g<b; g++){   //g=l
			for(int gg=0; gg<num_assets; gg++){
			//density_product*=density(X[i][g][gg], X[i+1][h][gg], sigma[gg], r, delta[gg], delta_t);
			density_product*=densityGPU(*three_dim_indexGPU(X_device, i, g, gg, m, b), *three_dim_indexGPU(X_device, (i+1), h, gg, m, b), sigma_device[gg], r, delta_device[gg], delta_t);

			}
		sum+=(1/((double)b))*density_product;
		}
	w_s=w_s/sum;	
	
	*two_dim_indexGPU(S_Weights, i, h, m, b)=w_s;
	}

}


__global__ void PathEstimatorKernel(double* X_device, double* W_device, double* V_device, double* delta_device, double* sigma_device, double* X0_device, int N, double strike, double r, double delta_t, int b, double m, int num_assets, curandState_t* states, double* results_dev, double* asset_amount_device){

//GeometricPayOffPut thePayOff(strike);
//GeometricPayOffPut payoff(strike);

//enum Containers { vector, matrix };

//Containers Vector = vector;
//Containers Matrix = matrix;


double v_0, S_i, Z, C, H, sum, weight; //, w_s, sum_Z;
//srand((unsigned)time(NULL));
//std::random_device rd;
//std::default_random_engine generator;
//generator.seed( rd() );
//std::normal_distribution<double> distribution(0.0,1.0);
/// ARRAY CODE
int m_int=(int)m;
const int S_N=(m_int)*num_assets;
const int S_W_N=(m_int)*b;

double* S;
S= new double[S_N]; 


double* S_Weights;
S_Weights=new double[S_W_N];

int idx =blockDim.x*blockIdx.x + threadIdx.x;



int i=0;
do {
	if(i==0){
		for(int ll=0; ll<num_assets; ll++){
			//Z=boxmuller();
			// NEED TO CHANGE THE RANDOM NUMBER GENERATOR	
			//Z=distribution(generator);
			Z=curand_normal_double(&states[idx]);
			S_i=X0_device[ll] +  (r-delta_device[ll]-0.5*pow(sigma_device[ll], 2))*delta_t + sigma_device[ll]*sqrt(delta_t)*Z;
			//tempnodevector.push_back(S_i);
			*two_dim_indexGPU(S, i, ll, m, num_assets)=S_i;			
		}
	}

	else{
		for(int jj=0; jj<num_assets; jj++){
			//Z=boxmuller();
			//Z=distribution(generator);
			Z=curand_normal_double(&states[idx]);
			S_i=(*two_dim_indexGPU(S, (i-1), jj, m, num_assets)) +  (r-delta_device[jj]-0.5*pow(sigma_device[jj], 2))*delta_t + sigma_device[jj]*sqrt(delta_t)*Z;
			//tempnodevector.push_back(S_i);
			*two_dim_indexGPU(S, i, jj, m, num_assets)=S_i;
		}
	}


if(i<m-1){
//S_weights(tempvec, S_Weights, X, S, m, b, sigma, delta, delta_t, asset_amount, r, i  );
S_weights(S_Weights, X_device, S, m, b, sigma_device, delta_device, delta_t, num_assets, r, i );
}

double con_val=0; //continuation value variable
	sum=0;

	if(i==m-1){
	C=0;//continuation value at the last time step
	}
	
	else{
		for(int k=0; k<b; k++){	
			weight=*two_dim_indexGPU(S_Weights, i, k, m, b);
			//con_val=V[(m-1)-i-1][k];
			con_val=*two_dim_indexGPU(V_device, (m-1-i-1), k, m, b);
			sum+=(weight) * (con_val); 			
		}
	
        //con_val=inner_product(b, first_vector, second_vector);
	
    
        C=(1/(double)b)*sum; //continuation value
//	C=(1/(double)b)*con_val;
	}	
	

//H=Payoff(S, strike, asset_amount, i)*exp(-r*delta_t*((i+1)));
//H=thePayOff(S, i, 0, m, num_assets, Vector, num_assets)*exp(-r*delta_t*((i+1)));
H= GeometricPayOffPutV(S, i, m, num_assets, num_assets, strike)*exp(-r*delta_t*((i+1)));

i=i+1;
}while(H<C);

v_0=H;

results_dev[idx]=v_0;


delete[] S;
delete[] S_Weights;
//return v_0;

}

double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma[], double delta[], double X0[], double* X, double* W, double* V, double asset_amount[], int num_assets, int Path_estimator_iterations ){

cudaError_t error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }
;

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
int W_N=(m_int) * b * b;
int V_N=(m_int) * b;
int delta_N= num_assets;
int sigma_N=num_assets;
int X0_N=num_assets;
int asset_amount_N = num_assets;

double* X_device;
double* V_device;
double* W_device;
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

cudaMalloc((void**) &W_device, W_N*sizeof(double) );
cudaMemcpy(W_device, W, W_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &X0_device, X0_N*sizeof(double) );
cudaMemcpy(X0_device, X0_host, X0_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &sigma_device, sigma_N*sizeof(double) );
cudaMemcpy(sigma_device, sigma_host, sigma_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &delta_device, delta_N*sizeof(double) );
cudaMemcpy(delta_device, delta_host, delta_N*sizeof(double), cudaMemcpyHostToDevice);

cudaMalloc((void**) &asset_amount_device, asset_amount_N*sizeof(double) );
cudaMemcpy(asset_amount_device, asset_amount_host, asset_amount_N*sizeof(double), cudaMemcpyHostToDevice);

dim3 gridDim(1);
dim3 blockDim(N);


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


curandState_t* states;

cudaMalloc((void**) &states, N * sizeof(curandState_t));

init<<<gridDim, blockDim>>>(time(0), states);

cudaDeviceSynchronize();

 
error = cudaGetLastError();

  if( error != cudaSuccess )
  {
    std::cout << cudaGetErrorString(error) << std::endl;
    printf("found at line %d\n", __LINE__);
    exit(1);
  }

//printf("after");
PathEstimatorKernel<<<gridDim, blockDim>>>(X_device, W_device, V_device, delta_device, sigma_device, X0_device, N, strike, r, delta_t, b,  m, num_assets, states, results_dev, asset_amount_device);

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
    printf("found at line %d\n", __LINE__);
    exit(1);
  }


double result=0;
for(int f=0; f<Path_estimator_iterations; f++){
result+=results[f];
//printf()
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
cudaFree(W_device);
cudaFree(sigma_device);
cudaFree(delta_device);
cudaFree(X0_device);
cudaFree(results_dev);
cudaFree(asset_amount_device);

//cudaDeviceReset();

return result;
}
