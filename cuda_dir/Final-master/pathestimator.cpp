#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <random>
#include "Payoff.h"
#include "enum_header.h" 



//declare transition density function. This is contained in meshgen.cpp

double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t);
//declare the box muller function. This is contained in meshgen.cpp
//double boxmuller();

double inner_product(int N, std::vector<double>& first_vector, std::vector<double>& second_vector);
//this returns the payoff value

double* two_dim_index(double* vector, int i, int j, double m, int b);

double* three_dim_index(double* matrix, int i, int j, int k, double m, int b);


/*double Payoffder.h (std::vector<std::vector<double> >& S, double k, std::vector<double>& asset_amount, int i){
double h;
//h=k-x;
//h=0;
h=1;

//for(int l=0; l<asset_amount.size(); l++){
//	h+=asset_amount[l]*exp(S[i][l]);
//}
for(int l=0; l<asset_amount.size(); l++){
         h*=exp(S[i][l]);
}
h=pow(h,1.0/(asset_amount.size()));
h=h-k; 
//h=h-k;
	if(h<0){
	h=0;
	}

return h;
}
*/


void S_weights(double* S_Weights, double* X, double* S, double m, int b, double sigma[], double delta[], double delta_t, int num_assets, double r , int i ){

double density_product, sum, w_s;

//this for-loop generates the S weights 

//std::cout<<"t="<<t<<std::endl;
//tempvec.clear();
	for(int h=0; h<b; h++){   //h=k
	sum=0;
	w_s=1;
		for(int kk=0; kk<num_assets; kk++){
		//w_s*=density(S[i][0][kk], X[i+1][h][kk], sigma[kk], r, delta[kk], delta_t);
		w_s*=density(*two_dim_index(S, i, kk, m, num_assets), *three_dim_index(X, (i+1), h, kk, m, b), sigma[kk], r, delta[kk], delta_t);
		}

	density_product=1;
	
		for(int g=0; g<b; g++){   //g=l
			for(int gg=0; gg<num_assets; gg++){
			//density_product*=density(X[i][g][gg], X[i+1][h][gg], sigma[gg], r, delta[gg], delta_t);
			density_product*=density(*three_dim_index(X, i, g, gg, m, b), *three_dim_index(X, (i+1), h, gg, m, b), sigma[gg], r, delta[gg], delta_t);

			}
		sum+=(1/((double)b))*density_product;
		}
	w_s=w_s/sum;	
	//tempvec.push_back(w_s);
*two_dim_index(S_Weights, i, h, m, b)=w_s;
	}

//S_Weights.push_back(tempvec); //vector storing S weights


}




//This function returns a low bias price of an option.
/*double PathEstimator(double strike, double r, double delta_t, int b, double m, std::vector<double>& sigma, std::vector<double>& delta, std::vector<double>& X0, std::vector<std::vector< std::vector<double> > >& X , std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V, std::vector<double>& asset_amount, const PayOff& thePayOff)
*/
double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma[], double delta[], double X0[], double* X, double* W, double* V, double asset_amount[], const PayOff& thePayOff, int num_assets){

Containers Vector = vector;
Containers Matrix = matrix;

double v_0, S_i, Z, C, H, sum, weight, w_s, sum_Z;
//srand((unsigned)time(NULL));
std::random_device rd;
std::default_random_engine generator;
generator.seed( rd() );
std::normal_distribution<double> distribution(0.0,1.0);
//// VECTOR CODE
/*
//Simulated path for sub optimal stoppage
std::vector< std::vector< std::vector< double > > > S;
//temp vector of S weights
std::vector< double > tempvec;
//weights matrix for S
std::vector< std::vector<double> > S_Weights; 
//temp vector in simulated path loop
std::vector< std::vector< double > > nodevector;

std::vector<double> tempnodevector;
*/

/// ARRAY CODE
int m_int=(int)m;

double* S;
S = new double[(m_int)*num_assets];


double* S_Weights;
S_Weights = new double[(m_int)*b];

int i=0;
//simulated path loop
//for(int i=0; i<m; i++){
do {
//nodevector.clear();
//tempnodevector.clear();
//std::cout<<"i="<<i<<std::endl;
	if(i==0){
		for(int ll=0; ll<num_assets; ll++){
			//Z=boxmuller();
			Z=distribution(generator);
		//	std::cout<<Z<<std::endl;
			//S_i=X0[ll] * (exp((r-delta[ll]-0.5*sigma[ll]*sigma[ll])*delta_t + sigma[ll]*sqrt(delta_t)*Z));//the second value in the simulated path
			S_i=X0[ll] +  (r-delta[ll]-0.5*pow(sigma[ll], 2))*delta_t + sigma[ll]*sqrt(delta_t)*Z;
			//tempnodevector.push_back(S_i);
			*two_dim_index(S, i, ll, m, num_assets)=S_i;			
		}
	}

	else{
		for(int jj=0; jj<num_assets; jj++){
			//Z=boxmuller();
			Z=distribution(generator);
			//S_i=S[i-1][jj]*(exp((r-delta[jj]-0.5*sigma[jj]*sigma[jj])*delta_t + sigma[jj]*sqrt(delta_t)*Z));//the simulate path values
			//S_i=S[i-1][0][jj] +  (r-delta[jj]-0.5*pow(sigma[jj], 2))*delta_t + sigma[jj]*sqrt(delta_t)*Z;
			S_i=(*two_dim_index(S, (i-1), jj, m, num_assets)) +  (r-delta[jj]-0.5*pow(sigma[jj], 2))*delta_t + sigma[jj]*sqrt(delta_t)*Z;
			//tempnodevector.push_back(S_i);
			*two_dim_index(S, i, jj, m, num_assets)=S_i;
		}
	}
//nodevector.push_back(tempnodevector);
//S.push_back(nodevector);//simulated path is stored in this vector 
//}end of simulated path loop
if(i<m-1){
//S_weights(tempvec, S_Weights, X, S, m, b, sigma, delta, delta_t, asset_amount, r, i  );
S_weights(S_Weights, X, S, m, b, sigma, delta, delta_t, num_assets, r, i );
}
/*
double density_product;

//this for-loop generates the S weights 
for(int t=0; t<(m-1); t++){
//std::cout<<"t="<<t<<std::endl;
tempvec.clear();
	for(int h=0; h<b; h++){   //h=k
	sum=0;
	w_s=1;
		for(int kk=0; kk<asset_amount.size(); kk++){
		w_s*=density(S[t][0][kk], X[t+1][h][kk], sigma[kk], r, delta[kk], delta_t);
		}

	density_product=1;
	
		for(int g=0; g<b; g++){   //g=l
			for(int gg=0; gg<asset_amount.size(); gg++){
			density_product*=density(X[t][g][gg], X[t+1][h][gg], sigma[gg], r, delta[gg], delta_t);
			}
		sum+=(1/((double)b))*density_product;
		}
	w_s=w_s/sum;	
	tempvec.push_back(w_s);
	}

S_Weights.push_back(tempvec); //vector storing S weights
}
*/
double con_val=0; //continuation value variable
//sub optimal path loop
//for(int i=0; i<m; i++){
	sum=0;

//std::vector<double> first_vector;
//std::vector<double> second_vector;

	if(i==m-1){
	C=0;//continuation value at the last time step
	}
	
	else{
		for(int k=0; k<b; k++){	
			//first_vector.push_back(S_Weights[i][k]);
                        //second_vector.push_back(V[(m-1)-i-1][k]);
			//weight=S_Weights[i][k];
			weight=*two_dim_index(S_Weights, i, k, m, b);
			//con_val=V[(m-1)-i-1][k];
			con_val=*two_dim_index(V, (m-1-i-1), k, m, b);
			sum+=(weight) * (con_val); 			
		}
	
        //con_val=inner_product(b, first_vector, second_vector);
	
    
        C=(1/(double)b)*sum; //continuation value
//	C=(1/(double)b)*con_val;
	}	
	

//H=Payoff(S, strike, asset_amount, i)*exp(-r*delta_t*((i+1)));
H=thePayOff(S, asset_amount, i, 0, m, num_assets, Vector, num_assets)*exp(-r*delta_t*((i+1)));
//check if continuation value is greater then the immediate payoff
/*	if(H>=C || i==m-1){
	v_0=H; 
	break;
	}*/

i=i+1;
}while(H<C);

v_0=H;

delete[] S;
delete[] S_Weights;



return v_0;

}
