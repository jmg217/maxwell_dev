#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>


#define PI 3.14159265358979323846
//this function initialises the kernel which constructs the mesh according to the independent path construction.
void mesh_generation(int b, int num_assets, double m, double X0[], double sigma[], double delta[], double asset_amount[], double* X, double strike, double r, double delta_t, curandState_t* States, curandState_t* statesi, int threads);

//this function initialises all the mesh weights kernels
void meshweights(double* W, double m, int b, double sigma[], double delta[], double r, double delta_t, double* X, int num_assets, double* weight_denominator);

//this function converts one dim vectors to one dim arrays
void one_dim_array(std::vector< double > &vals, double array[], int N);

//this function provides an indexing interface for 3-d matrices stored in 1 dim arrays
double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);

//this function provides an indexing interface for 2-d matrices stored in 1 dim arrays
double* two_dim_index(double* vector, int i, int j, double m, int b);

//this function initialises the mesh estimator kernel
double MeshEstimator(double strike, double r, double delta_t, int b, double m, double* X, double* W, double* V, double asset_amount[], int num_assets);

//this function initialises the pathestimator kernel.
double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma[], double delta[], double X0[], double* X, double* weight_denominator, double* V, double asset_amount[], int num_assets, int Path_estimator_iterations, int iterator, int Final_iteration, curandState_t* States, curandState_t* states, int threads);

//this function prints results for 1-dim options to 'OneDimHighPayoff.txt'. These results can be plotted using the program 'OneDimPayoff.py'
void print_high_payoff(int b, double m, double* X, double* V, double asset_amount[], double* W );

//this function prints results for 1-dim options to 'stocksimulation.txt'. These results can be plotted using the program 'simulation.py'
void SimulationPaths(int b, double m, double* X, double* V, double asset_amount[], double* W, double x );

//this kernel initialises the random seeds on each thread
__global__ void initialise(unsigned int seed, curandState_t* states) {
int idx=blockDim.x*blockIdx.x + threadIdx.x;

  // we have to initialize the state 
  curand_init(seed,  idx, 0, &states[idx]);
}


int main(){

srand((unsigned)time(NULL));



//begin timer
clock_t begin=clock();

//read in parameters from setting.txt
std::ifstream setting( "settings.txt" );
std::string line;
std::vector<std::string> settings;
int linenumber=0;
while(std::getline( setting, line))
    {
          if(linenumber%2==1)
              settings.push_back(line);
          linenumber++;
    }
setting.close();

int integer; 

std::vector < double > X0V;
std::vector < double > deltaV;
std::vector <double> sigmaV;
std::vector <double> asset_amountV;



std::istringstream ss(settings[0]);
std::string token;
while(std::getline(ss, token, ',')) 
    {
          X0V.push_back(atof(token.c_str()));
    }

double T = atof(settings[1].c_str());
double m = atof(settings[2].c_str());
double delta_t=T/m;

double v_0, V_0, vtotal_sum=0, Vtotal_sum=0;

double r= atof(settings[3].c_str());

std::istringstream ss2(settings[4]);
while(std::getline(ss2, token, ','))
    {
        deltaV.push_back(atof(token.c_str()));
    }

std::istringstream ss3(settings[5]);
while(std::getline(ss3, token, ','))
    {
        sigmaV.push_back(atof(token.c_str()));
    }

int Path_estimator_iterations=atof(settings[6].c_str());
double strike=atof(settings[7].c_str());
int b=atoi(settings[8].c_str());
int N=atoi(settings[9].c_str());
double quantile=atof(settings[10].c_str());
int num_assets=atof(settings[11].c_str());

std::istringstream ss4(settings[12]);
while(std::getline(ss4, token, ','))
    {
        asset_amountV.push_back(atof(token.c_str()));
    }


if(X0V.size() != num_assets || sigmaV.size() != num_assets || deltaV.size() !=num_assets || asset_amountV.size() !=num_assets){
          std::cout<<"Either the starting price, volatility, number of assets or dividend yield was not specified for all assets"<<std::endl;
           exit (EXIT_FAILURE);
       }

std::cout<<"The parameters of this simulation are:"<<std::endl;
       //Print these values to screen 
for(integer=0; integer<X0V.size(); integer++){
std::cout<<"Starting Price="<<X0V[integer]<<std::endl;
}

std::cout<<"Time to expiry="<<T<<"\n"<<"Number of time steps="<<m<<"\n"<<"interest rate="<<r<<std::endl;

for(integer=0; integer<sigmaV.size(); integer++){
    std::cout<<"volatility="<<sigmaV[integer]<<std::endl;
}

for(integer=0; integer<deltaV.size(); integer++){
    std::cout<<"dividend yield="<<deltaV[integer]<<std::endl;
}
std::cout<<"number of iterations over path estimator="<<Path_estimator_iterations<<"\n"<<"strike  price="<<strike<<"\n"<<"number of nodes per time step="<<b<<"\n"<<"number mesh generations="<<N<<"\n"<<"Number of Assets="<<num_assets<<std::endl; 

for(integer=0; integer<asset_amountV.size(); integer++){
    std::cout<<"asset amount="<<asset_amountV[integer]<<std::endl;
}

// CONVERT TO ARRAYS

double X0 [num_assets];
double delta [num_assets];
double sigma [num_assets];
double asset_amount [num_assets];

one_dim_array(X0V, X0, num_assets);
one_dim_array(deltaV, delta, num_assets);
one_dim_array(sigmaV, sigma, num_assets);
one_dim_array(asset_amountV, asset_amount, num_assets);


//V values from each iteration over meshes
std::vector< double > Vvector;
//v values from each iteration over meshes
std::vector< double > vvector;
//asset vector
std::vector< double > assets;
//1 d vector in Weightsgen for-loop
std:: vector<double> dim1temp;

std::vector<double> sortvector;
//std::cout<<"Before loop"<<std::endl;


// ARRAY CODE
int m_int= (int)m;

double* X;
int X_dim = (m_int) * b * (num_assets);
X= new double[X_dim];

double* W;
int W_dim = (m_int) * b * b;
W= new double[W_dim];

double* V;
int V_dim = (m_int) * b;
V = new double[V_dim];

 
double* weight_denominator;
int denom_dim = (m_int-1) * b;
weight_denominator =new double[denom_dim];

for(int init=0; init<num_assets; init++){
	X0[init]=log(X0[init]);
}	

int threads=Path_estimator_iterations;
if(b>Path_estimator_iterations){
threads=b;
}
curandState_t* States;
States= new curandState_t[threads];

//for-loop over different meshes
for(int iterator=0; iterator<N; iterator++){

	curandState_t* states;

	cudaMalloc((void**) &states, threads * sizeof(curandState_t));
	if(iterator==0){
		dim3 gridDim((int)ceil(threads/512.0));
		dim3 blockDim(512.0);
		initialise<<<gridDim, blockDim>>>(time(0), states);
		cudaDeviceSynchronize();
		cudaMemcpy(States, states, sizeof(curandState_t)*threads, cudaMemcpyDeviceToHost);
	}
	else{cudaMemcpy(states, States, threads*sizeof(curandState_t), cudaMemcpyHostToDevice);}

	mesh_generation(b, num_assets, m, X0, sigma, delta, asset_amount, X, strike, r, delta_t, States, states, threads);

	meshweights(W,  m, b, sigma, delta, r, delta_t, X, num_assets, weight_denominator);


	double check=0;
	//check all the weights from X0 are 1
	for(int e=0; e<b; e++){
		if(*three_dim_index(W, 0, e, 0, m, b, b)!=1){
			std::cout<<"there is an error with the weights. check that W[0][k][0]'s =1"<<std::endl;
		}
	}
	//check that the weights going into a node sum to 1
	for(int q=1; q<m; q++){ 
		for(int a=0; a<b; a++){
			check=0;
			for(int E=0; E<b; E++){
			
				check+=*three_dim_index(W, (q), a, E, m, b, num_assets);
			}
		}
	}

	V_0=MeshEstimator(strike, r, delta_t, b, m, X, W, V, asset_amount, num_assets);
	Vvector.push_back(V_0);//vector containing high bias option prices
	Vtotal_sum+=V_0;

	std::cout<<"High Bias price (V_0) for mesh iteration "<<iterator<<" is "<<V_0<<std::endl;

	v_0=PathEstimator(strike, r, delta_t, b,  m, sigma, delta, X0, X, weight_denominator, V, asset_amount, num_assets, Path_estimator_iterations, iterator, N, States, states, threads);

	cudaDeviceReset();

	vvector.push_back(v_0);
	vtotal_sum+=v_0;

	std::cout<<"Low Bias price (v_0) for mesh iteration "<<iterator<<" is "<<v_0<<std::endl;


}//this is the end of the loop over the whole process.

if(num_assets==1){
	print_high_payoff(b, m, X, V, asset_amount,W);
	SimulationPaths(b, m, X, V, asset_amount, W, X0[0] );
}

//Calculate V(N) and v(N)
V_0=(1/double(N))*Vtotal_sum;
v_0=(1/double(N))*vtotal_sum;


//calculate errors
double std_div_V=0, std_div_v=0, squaresumV=0, squaresumv=0, Verror=0, verror=0;

for(int h=0; h<N; h++){
	squaresumV+=(Vvector[h]-V_0)*(Vvector[h]-V_0);
	squaresumv+=(vvector[h]-v_0)*(vvector[h]-v_0);
}
std_div_V=sqrt((1/double(N))*squaresumV); //standard deviation of V
std_div_v=sqrt((1/double(N))*squaresumv); //standard deviation of v

double standardErrorV=std_div_V*(1/sqrt(double(N)));
double standardErrorv=std_div_v*(1/sqrt(double(N)));

Verror=quantile*standardErrorV;
verror=quantile*standardErrorv;

std::cout<<"V(N)_0="<<V_0<<"\t"<<"V error="<<Verror<<std::endl;
std::cout<<"v(N)_0="<<v_0<<"\t"<<"v error="<<verror<<std::endl;

double pointEst=(V_0+v_0)/2;
double EstimatedError=((Verror+V_0)-(v_0-verror))/(2*pointEst);
clock_t end =clock();
double elapsedtime=double(end-begin) / CLOCKS_PER_SEC;

std::ofstream outFile("results.txt", std::ios_base::app | std::ios_base::out);

outFile << N <<"\t"<< b <<"\t"<< Path_estimator_iterations<<"\t"<<exp(X0[0])<<"\t" << v_0 <<"\t"<< standardErrorv <<"\t"<< V_0 <<"\t"<< standardErrorV <<"\t"<< v_0-verror<<"\t"<<Verror+V_0 <<"\t"<<pointEst<<"\t"<<EstimatedError<<"\t" <<elapsedtime<< std::endl;

outFile.close();

std::ofstream LoutFile("latexresults.txt", std::ios_base::app | std::ios_base::out);

LoutFile <<std::fixed<<std::setprecision(3) << N <<"&"<< b <<"&"<< Path_estimator_iterations<<"&"<<exp(X0[0])<<"&" << v_0 <<"&"<< standardErrorv <<"&"<< V_0 <<"&"<< standardErrorV <<"&"<< v_0-verror<<"&"<<Verror+V_0 <<"&"<<pointEst<<"&"<<EstimatedError<<"&" <<elapsedtime<< std::endl;

outFile.close();

delete[] X;
delete[] W;
delete[] V;
delete[] weight_denominator;
delete[] States;
return 0;


}


