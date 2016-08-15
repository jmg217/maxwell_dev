#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>
#include "Payoff.h"
#include "enum_header.h"
//double inner_control_mesh(int i, int j, int b, double r, double delta_t, double m, std::vector< std::vector< std::vector<double> > >& W, std::vector<std::vector< std::vector<double> > >& X, std::vector< std::vector<double> >& V);


double inner_control_mesh(int i, int j, int b, double r, double delta_t, double m, double* W, double* X, double* V);

double* three_dim_index(double* matrix, int i, int j, int k, double m, int b);

double* two_dim_index(double* vector, int i, int j, double m, int b);

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
double MeshEstimator(double strike, double r, double delta_t, int b, double m, double* X, double* W, double* V, double asset_amount[], const PayOff& thePayOff, int num_assets){
Containers Matrix = matrix;
Containers Vector = vector;
double H; //payoff variable 
double C; //continuation value variable
double sum, V_0;

// temp vector in Estimator loop
//std::vector< double > tempvec;


//Mesh Estimator loop

for(int i=0; i<m; i++){
//tempvec.clear();

	for(int j=0; j<b; j++){
		if(i==0){

		//H=payoff(X, strike, asset_amount, m-1-i, j)*exp(-r*delta_t*(m-i));	
		H=thePayOff(X, asset_amount, m-1-i,j, m, b, Matrix, num_assets)*exp(-r*delta_t*(m-i));
		//tempvec.push_back(H);
		//V[m*(j)+(i)]=H;
		*two_dim_index(V, i, j, m, b)=H;
		}
	
		else{
//continue to develope continuation vale			

//NEW CONTROL VARIATE VERSION OF CONTIN VAL

		C=inner_control_mesh(i, j, b, r, delta_t, m, W, X, V);


////////////PREVIOUS CALCULATION of CONTIN VAL. 
	/*	sum=0;
			for(int k=0; k<b; k++){
//std::cout<< sum<<std::endl;			
			sum+=(W[(m-i)][k][j])*V[i-1][k]; //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]	
				
			}
	
		C=(1/((double)b))*sum; //continuation value
			
	*/
////////////END OF CALCULATION OF PREVIOUS CONTIN VAL. (WITHOUT CONTROL VARIATES)
	
		/*if(m-i==2 && j == 50){
		std::cout<<"contin value="<< C<<std::endl; 
		}*/
		//H=payoff(X, strike, asset_amount, m-1-i, j)*exp(-r*delta_t*(m-i));
		H=thePayOff(X, asset_amount, m-1-i,j, m, b,  Matrix, num_assets)*exp(-r*delta_t*(m-i));
			if(H>=C){
				//tempvec.push_back(H);
				*two_dim_index(V, i, j, m, b)=H;
			
			}

			else{
				//tempvec.push_back(C);
				*two_dim_index(V, i, j, m, b)=C;
			}	
		}
	}
//V.push_back(tempvec); //high bias option value matrix
}


//sum over option values at the second time step
sum=0;
int m_int =(int)m;
for(int k=0; k<b; k++){
//sum+=V[m-1][k];
sum+=*two_dim_index(V, (m_int-1), k, m, b);        
}
//this is the high bias option value at time 0
V_0=(1/((double)b))*sum;


return V_0;
}
