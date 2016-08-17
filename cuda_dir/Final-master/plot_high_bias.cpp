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
//double x=0;

std::ofstream outFile("OneDimHighPayoff.txt", std::ios_base::app | std::ios_base::out);

for (int t=0; t<m; t++){


	for (int i=0; i<b; i++){
	outFile << (m-t)/10 <<"\t"<< *two_dim_index(V, t, i, m, b)<<"\t"<< exp(*three_dim_index(X, m-1-t, i, 0, m, b, 1)) << std::endl;	


	}

}

outFile.close();

}

void SimulationPaths(int b, double m, double* X, double* V, double asset_amount[], double* W, double x ){


std::ofstream outFile("stocksimulation.txt", std::ios_base::app | std::ios_base::out);

for (int t=0; t<m+1; t++){

outFile<<t/m<<"\t";
        for (int i=0; i<b; i++){
		if(t==0){outFile <<exp(x)<<"\t";}
		else{ outFile <<exp(*three_dim_index(X, t-1, i, 0, m, b, 1))<<"\t";}
        }
outFile<<std::endl;
}

outFile.close();

}

