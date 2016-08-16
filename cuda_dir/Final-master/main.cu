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
//#include <random>
#include "Payoff.h"
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>


#define PI 3.14159265358979323846
double outer_control_mesh(std::vector<double>& Vvector, std::vector<double>& EuroVals, int N, int b, double EuroGeoOption );

void mesh_generation(int b, int num_assets, double m, double X0[], double sigma[], double delta[], double asset_amount[], double* X, double strike, double r, double delta_t, curandState_t* States, curandState_t* statesi, int threads);

void meshweights(double* W, double m, int b, double sigma[], double delta[], double r, double delta_t, double* X, int num_assets, double* weight_denominator);

//void meshweights(double* W, double m, int b, double sigma[], double delta[], double r, double delta_t, double* X, int num_assets, double* weight_denominator, int threads);

void one_dim_array(std::vector< double > &vals, double array[], int N);

double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);

double* two_dim_index(double* vector, int i, int j, double m, int b);

//void cuda_init_random(int Path_estimator_iterations);

/*
//This function is contained within meshestimator.cpp
double MeshEstimator(double strike, double r, double delta_t, int b, double m,  std::vector<std::vector< std::vector<double> > >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V, std::vector<double>& asset_amount, const PayOff& thePayOff);
*/

//double MeshEstimator(double strike, double r, double delta_t, int b, double m,  double* X, double* W, double* V, double asset_amount[], const PayOff& thePayOff, int num_assets, std::vector <double>& EuroVals);

double MeshEstimator(double strike, double r, double delta_t, int b, double m, double* X, double* W, double* V, double asset_amount[], int num_assets);

/*
//This function is contained within patestimator.cpp.
double PathEstimator(double strike, double r, double delta_t, int b, double m, std::vector<double>& sigma, std::vector<double>& delta, std::vector<double>& X0, std::vector<std::vector< std::vector<double> > >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V, std::vector<double>& asset_amount, const PayOff& thePayOff);
*/

double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma[], double delta[], double X0[], double* X, double* weight_denominator, double* V, double asset_amount[], int num_assets, int Path_estimator_iterations, int iterator, int Final_iteration, curandState_t* States, curandState_t* states, int threads);

void print_high_payoff(int b, double m, double* X, double* V, double asset_amount[], double* W );

__global__ void initialise(unsigned int seed, curandState_t* states) {
int idx=blockDim.x*blockIdx.x + threadIdx.x;

  // we have to initialize the state 
  curand_init(seed, // the seed can be the same for each core, here we pass the time in from the CPU 
              idx, // the sequence number should be different for each core (unless you want all
                            // cores to get the same sequence of numbers for some reason - use thread id! 
              0, // the offset is how much extra we advance in the sequence for each call, can be 0 
              &states[idx]);
}


/*
__host__ __device__ double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t)// This function returns the transition density between node values.
{
double f=0, x=0;
//x=(1/(sigma*sqrt(delta_t)))*(log(Xnew)-log(Xold)-(r-delta-0.5*sigma*sigma)*delta_t);
x=(1/(sigma*sqrt(delta_t)))*(Xnew-Xold-(r-delta-0.5*sigma*sigma)*delta_t);
//f= (1/(sigma*sqrt(delta_t)*Xnew))*(1/(sqrt(2*PI)))*exp(-0.5*x*x); // this is the transition density
f= (1/(sigma*sqrt(delta_t)))*(1/(sqrt(2*PI)))*exp(-0.5*x*x);
return f;
}
*/

int main(){

srand((unsigned)time(NULL));
//std::random_device rd;
//std::default_random_engine generator;
//generator.seed( rd() );
//std::normal_distribution<double> distribution(0.0,1.0);

//Now we read in the parameters from settings.txt
//new code

clock_t begin=clock();

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

//double double_num;
int integer; 

std::vector < double > X0V;
std::vector < double > deltaV;
std::vector <double> sigmaV;
std::vector <double> asset_amountV;

std::vector <double> EuroVals;


std::istringstream ss(settings[0]);
std::string token;
while(std::getline(ss, token, ',')) 
    {
          X0V.push_back(atof(token.c_str()));
    }

double T = atof(settings[1].c_str());
double m = atof(settings[2].c_str());
double delta_t=T/m;
//int Rn;
double v_0, V_0, vtotal_sum=0, Vtotal_sum=0;
 //Z, Xi, Xj, v_sum, sum_Z=0, vtotal_sum=0, Vtotal_sum=0;
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
double EuroGeoOption=atof(settings[13].c_str());

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
std::cout<<"number of iterations over path estimator="<<Path_estimator_iterations<<"\n"<<"strike  price="<<strike<<"\n"<<"number of nodes per time step="<<b<<"\n"<<"number mesh generations="<<N<<"\n"<<"Number of Assets="<<num_assets<<"\n"<<"European Option="<<EuroGeoOption<<std::endl; 

for(integer=0; integer<asset_amountV.size(); integer++){
    std::cout<<"asset amount="<<asset_amountV[integer]<<std::endl;
}

//////////////////////////////////////////
//// DECLARE AN INSTANCE OF THE PAYOFF////
//////////////////////////////////////////

//GeometricPayOffCall payoff(strike);
//GeometricPayOffPut *payoff_dev;  

// CONVERT TO ARRAYS

double X0 [num_assets];
double delta [num_assets];
double sigma [num_assets];
double asset_amount [num_assets];

one_dim_array(X0V, X0, num_assets);
one_dim_array(deltaV, delta, num_assets);
one_dim_array(sigmaV, sigma, num_assets);
one_dim_array(asset_amountV, asset_amount, num_assets);


// VECTOR CODE
/*
//Mesh matrix
std::vector< std::vector< std::vector<double> > > X;
//WEIGHTS 3-dimensional matrix for step 1 and beyond
std::vector< std::vector< std::vector<double> > > W;
//2-d temp vector in MeshGen for-loop
std::vector < std::vector< double > > myvector;
//temp vecotr in MeshGen for-loop
std::vector<double > nodevector;
//2 d temp vector in WeightsGen for-loop
std::vector< std::vector<double> > dim2temp;
//1 d vector in Weightsgen for-loop
std:: vector<double> dim1temp;
//mesh estimator high bias 2-d matrix
std::vector< std::vector<double> > V;
*/
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

//std::cout<<"before"<<std::endl; 
double* weight_denominator;
int denom_dim = (m_int-1) * b;
weight_denominator =new double[denom_dim];
//std::cout<<"after"<<std::endl;
for(int init=0; init<num_assets; init++){
	X0[init]=log(X0[init]);
}	

int threads=Path_estimator_iterations;
if(b>Path_estimator_iterations){
threads=b;
}
curandState_t* States;
States= new curandState_t[threads];
//curandState_t* states;
//curandState_t* states;
//cudaMalloc((void**) &states, threads * sizeof(curandState_t));
//dim3 gridDim((int)ceil(threads/512.0));
//dim3 blockDim(512.0);
//initialise<<<gridDim, blockDim>>>(time(0), states);
//cudaDeviceSynchronize();
//cudaMemcpy(States, states, sizeof(curandState_t)*threads, cudaMemcpyDeviceToHost);

//for-loop over different meshes
for(int iterator=0; iterator<N; iterator++){
//X.clear();
//W.clear();
//V.clear();
//for-loop to generate the mesh

curandState_t* states;
//curandState_t* states;
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

//THIS IS THE SERIAL MESH GENERATION CODE
/*
for(int i=0; i<m; i++){


	//myvector.clear();

	if(i==0){



	for(int l=0; l<b; l++){

		//nodevector.clear();

		for(int ll=0; ll<num_assets; ll++){
			//Z=boxmuller();//standard normally distributed variable 
			Z=distribution(generator);
		//std::cout<<"meshgen="<<Z<<std::endl;
			//Xi=X0[ll] * (exp ((r-delta[ll]-0.5*pow(sigma[ll], 2))*delta_t + sigma[ll]*sqrt(delta_t)*Z));//node value at the second time step
		//	Xi=X0[ll] +  (r-delta[ll]-0.5*pow(sigma[ll], 2))*delta_t + sigma[ll]*sqrt(delta_t)*Z;
			//nodevector.push_back(Xi);
 
	*three_dim_index(X, i, l, ll, m, b, num_assets) = X0[ll] +  (r-delta[ll]-0.5*pow(sigma[ll], 2))*delta_t + sigma[ll]*sqrt(delta_t)*Z;
		
	//X[m*b*(ll)+m*(l)+(i)]=X0[ll] +  (r-delta[ll]-0.5*pow(sigma[ll], 2))*delta_t + sigma[ll]*sqrt(delta_t)*Z;
		//three_dim_index(X, i, l, ll, m, b, Xi);

		}
	//	myvector.push_back(nodevector);	//store the value in a temp vector
		
	}
	}

	if(i>0){
	
	for(int j=0; j<b; j++){
	
	//	nodevector.clear();
	//	Rn=UniRandom(b);
//std::cout<<Rn<<std::endl;
		for(int jj=0; jj<num_assets; jj++){
			//std::cout<<"in loop"<<"\t"<<Rn<<std::endl;
			//Z=boxmuller();
			Z=distribution(generator); 
			//std::cout<<"meshgen="<<Z<<std::endl;
			//Xi=X[i-1][j][jj];
			//Xi=X[m*b*(jj)+m*(j)+(i-1)];
			Xi=*three_dim_index(X, (i-1), j, jj, m, b, num_assets);
			//Xj=Xi * (exp ((r-delta[jj]-0.5*pow(sigma[jj], 2))*delta_t + sigma[jj]*sqrt(delta_t)*Z));
			//X[m*b*(jj)+m*(j)+(i)]=Xi +  (r-delta[jj]-0.5*pow(sigma[jj], 2))*delta_t + sigma[jj]*sqrt(delta_t)*Z; 
			*three_dim_index(X, i, j, jj, m, b, num_assets)=Xi +  (r-delta[jj]-0.5*pow(sigma[jj], 2))*delta_t + sigma[jj]*sqrt(delta_t)*Z;
			//nodevector.push_back(Xj);
		}	
	//	myvector.push_back(nodevector);
	}
	}

//X.push_back(myvector);

}
*/
//std::cout<<"after loop"<<std::endl;
meshweights(W,  m, b, sigma, delta, r, delta_t, X, num_assets, weight_denominator);
//std::cout<<"after"<<std::endl;

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
			//check+=W[m*b*(E)+m*(a)+(q)];
			check+=*three_dim_index(W, (q), a, E, m, b, num_assets);
		}
	}
}

//V_0=MeshEstimator(strike, r, delta_t, b, m, X, W, V, asset_amount, payoff, num_assets, EuroVals);//high bias option price

//V_0= MeshEstimator(strike, r, delta_t, b, m, X, W, V, asset_amount[], num_assets);

V_0=MeshEstimator(strike, r, delta_t, b, m, X, W, V, asset_amount, num_assets);
Vvector.push_back(V_0);//vector containing high bias option prices
Vtotal_sum+=V_0;

std::cout<<"High Bias price (V_0) for mesh iteration "<<iterator<<" is "<<V_0<<std::endl;

//average over path estimators

/*
v_sum=0;
for(int f=0; f<Path_estimator_iterations; f++){
v_sum=PathEstimator(strike, r, delta_t, b,  m, sigma, delta, X0, X, W, V, asset_amount, payoff, num_assets);
}
*/
//v_0=(1/double(Path_estimator_iterations))*v_sum;
v_0=PathEstimator(strike, r, delta_t, b,  m, sigma, delta, X0, X, weight_denominator, V, asset_amount, num_assets, Path_estimator_iterations, iterator, N, States, states, threads);

cudaDeviceReset();

vvector.push_back(v_0);
vtotal_sum+=v_0;

std::cout<<"Low Bias price (v_0) for mesh iteration "<<iterator<<" is "<<v_0<<std::endl;


}//this is the end of the loop over the whole process.
print_high_payoff(b, m, X, V, asset_amount,W);
//Calculate V(N) and v(N)
V_0=(1/double(N))*Vtotal_sum;
v_0=(1/double(N))*vtotal_sum;

if(EuroGeoOption>0){
std::cout<<"here"<<std::endl;
V_0= outer_control_mesh(Vvector, EuroVals, N, b, EuroGeoOption );
}
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


