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

__host__ __device__ double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);

__host__ __device__ double* two_dim_index(double* vector, int i, int j, double m, int b);


__device__ double GeometricPayOffCallM(double* X, int i, int j, double m, int b, int num_assets, double Strike){
double h;
h=1;
for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);
                h*= exp(*three_dim_index(X, i, j, l, m, b, num_assets));
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
                h*=exp(*three_dim_index(X, i, j, l, m, b, num_assets));
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

__device__ double GeometricPayOffCallV(double* X, double m, int b, int num_assets, double Strike){
double h;
h=1;
for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);

                //h*= exp(*two_dim_indexGPU(X, i, l, m, b));
                h*=exp(X[l]);
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

__device__ double GeometricPayOffPutV(double* X, double m, int b, int num_assets, double Strike){
double h;
h=1;
for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);
                //h*= exp(*two_dim_indexGPU(X, i, l, m, b));
                h*=exp(X[l]);
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

