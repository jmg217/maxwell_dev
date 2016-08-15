#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);

double* two_dim_index(double* vector, int i, int j, double m, int b);

void print_high_payoff(int b, double m, double* X, double* V, double asset_amount[], double* W ){
double x=0;
/*
std::ofstream outFile("highpayoff.txt", std::ios_base::app | std::ios_base::out);

for (int t=0; t<m; t++){

	for (int i=0; i<b; i++){
	x=0;
		for(int tt=0; tt<asset_amount.size(); tt++){
		//x+=asset_amount[tt]*X[m-1-t][i][tt];
		x+=asset_amount* (*three_dim_index(X, m-1-t, i, tt));
		}
	//outFile << m-t <<"\t"<< exp(x) <<"\t"<< V[t][i]<<"\t"<<X[m-1-t][i][0]<<"\t"<<X[m-1-t][i][1]<< std::endl;

	outFile << m-t <<"\t"<< exp(x) <<"\t"<< *two_dim_index(V, t, i, m, b)<<"\t"<< *three_dim_index(X, m-1-t, i, 1) <<"\t"<< *three_dim_index(X, m-1-t, i, 1) << std::endl;	
//	if(t>0){
//	for(int ttt=0; ttt<b; ttt++){
//	outFile <<ttt<<"\t"<<V[t-1][ttt]<<"\t"<<W[m-t][ttt][i]<<std::endl;
//	}
//	}
	}
}

outFile.close();
*/
}


