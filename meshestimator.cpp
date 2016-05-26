#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>
#include "Payoff.h"
//this function returns the payoff value
double payoff(std::vector<std::vector< std::vector<double> > >& X, double k, std::vector<double>& asset_amount, int i, int j){
double h;
//h=k-x;
//h=0;
h=1;
/*for(int l=0; l<asset_amount.size(); l++){
	h+=asset_amount[l]*exp(X[i][j][l]);
}*/
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

//this function returns the high bias mesh price
double MeshEstimator(double strike, double r, double delta_t, int b, double m,  std::vector<std::vector< std::vector<double> > >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V, std::vector<double>& asset_amount, const PayOff& thePayOff){

double H; //payoff variable 
double C; //continuation value variable
double sum, V_0;

// temp vector in Estimator loop
std::vector< double > tempvec;


//Mesh Estimator loop

for(int i=0; i<m; i++){
tempvec.clear();

	for(int j=0; j<b; j++){
		if(i==0){

		//H=payoff(X, strike, asset_amount, m-1-i, j)*exp(-r*delta_t*(m-i));	
		H=thePayOff(X, asset_amount, m-1-i,j)*exp(-r*delta_t*(m-i));
		tempvec.push_back(H);
		}
	
		else{
//continue to develope continuation vale			
		sum=0;
			for(int k=0; k<b; k++){
//std::cout<< sum<<std::endl;			
			sum+=(W[(m-i)][k][j])*V[i-1][k]; //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]	
			/*if(m-i==2 && j == 50){
			std::cout<<"k="<<k<<"\t"<<"weight="<<W[(m-i)][k][j]<<"\t"<<"V_i+1="<<V[i-1][k]<<std::endl;
			}*/	
			}
	
		C=(1/((double)b))*sum; //continuation value
			
		
		/*if(m-i==2 && j == 50){
		std::cout<<"contin value="<< C<<std::endl; 
		}*/
		//H=payoff(X, strike, asset_amount, m-1-i, j)*exp(-r*delta_t*(m-i));
		H=thePayOff(X, asset_amount, m-1-i,j)*exp(-r*delta_t*(m-i));
			if(H>=C){
				tempvec.push_back(H);
			
			}

			else{
				tempvec.push_back(C);
			
			}	
		}
	}
V.push_back(tempvec); //high bias option value matrix
}


//sum over option values at the second time step
sum=0;
for(int k=0; k<b; k++){
sum+=V[m-1][k];        
}
//this is the high bias option value at time 0
V_0=(1/((double)b))*sum;


return V_0;
}
