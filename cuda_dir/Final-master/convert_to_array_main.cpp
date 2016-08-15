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
#include "Payoff.h"

double** two_dim_array(std::vector< std::vector<double> > &vals, int N, int M);
void one_dim_array(std::vector< double > &vals, double array[], int N);
double*** three_dim_array(std::vector<std::vector< std::vector<double> > > &vals, int N, int M, int L);

int main(){

std::vector< std::vector<double> > V;
std::vector<double> tempv;

int N=5;
int M=5;
int L=5;
/*
for(int i=0; i<M; i++){
tempv.clear();
	for(int j=0; j<N; j++){
	tempv.push_back(i*j);
	}
V.push_back(tempv);
}
double **p;
  

std::cout<<"before"<<std::endl;
p=two_dim_array(V, M, N);
std::cout<<"after"<<std::endl;



for(int l=0; l<M; l++){

	for(int k=0; k<N; k++){
	std::cout<<V[l][k]<<"\t";
	}
std::cout<<std::endl;
}

for(int t=0; t<M; t++){

        for(int z=0; z<N; z++){
        std::cout<<p[t][z]<<"\t";
        }
std::cout<<std::endl;
}
*/

std::vector<double > v;
double g [N];

for(int y=0; y<N; y++){

v.push_back(y);
}

one_dim_array(v, g, N);

for(int iter=0; iter<N; iter++){
std::cout<<g[iter]<<"\t"<<v[iter]<<std::endl;
}
/*
std::vector<std::vector<std::vector<double> > > Matrix;
double*** Arr;

std::vector<double> tempvec;
std::vector<std::vector<double> > tempmatrix;

M=5;
N=5;
L=5;

for(int xx=0; xx<M; xx++){
tempmatrix.clear();
std::cout<<"xx="<<xx<<std::endl;
	for(int dd=0; dd<N; dd++){
	tempvec.clear();
std::cout<<"dd="<<dd<<std::endl;

		for(int zz=0; zz<L; zz++){
std::cout<<"zz="<<zz<<std::endl;

		tempvec.push_back(zz*dd*xx);	
		}

	tempmatrix.push_back(tempvec);
	}
Matrix.push_back(tempmatrix);
}

std::cout<<"mat "<<Matrix[0][1][0];
Arr=three_dim_array(Matrix, M, N, L);
*/
return 0;
}
