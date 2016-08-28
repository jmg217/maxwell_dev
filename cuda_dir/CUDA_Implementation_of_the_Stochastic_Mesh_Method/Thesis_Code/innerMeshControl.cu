#include <cuda.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include "enum_header.h"
#include <unistd.h>
#include <stdio.h>

__host__ __device__ double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);

__host__ __device__ double* two_dim_index(double* vector, int i, int j, double m, int b);


__device__ double inner_control_meshME(int i, int j, int b, double r, double delta_t, double m, double* W_device, double* X_device, double* V_device, int num_assets ){

int m_int= (int)m;
double ControlMeshContinVal=0;
double sum=0, ContinVal=0, stock_expectation=0, true_stock_expectation=0, numerator=0, denominator=0, beta=0;

for(int k=0; k<b; k++){


        sum+= (*three_dim_index(W_device, (m_int-i), k, j, m, b, b)) * (*two_dim_index(V_device, (i-1), k, m, b));
}

ContinVal=(1/((double)b))*sum;

sum=0;
for(int l=0; l<b; l++){

        sum+=(*three_dim_index(W_device, (m_int-i), l, j, m, b, b)) * exp((*three_dim_index(X_device, (m_int-i), l,0 , m, b, num_assets)));


}
stock_expectation=(1/((double)b))*sum;
true_stock_expectation=exp(*three_dim_index(X_device, (m_int-i-1), j, 0, m, b, num_assets)) * exp(r*delta_t);//old

for(int p=0; p<b; p++){

numerator += ( (*three_dim_index(W_device, (m_int-i), p, j, m, b, b)) * exp( *three_dim_index(X_device, (m_int-i), p, 0, m, b, num_assets) ) - stock_expectation ) * ( (*three_dim_index(W_device, (m_int-i), p, j, m, b, b)) * (*two_dim_index(V_device, (i-1), p, m, b)) - ContinVal );

denominator += pow( ( (*three_dim_index(W_device, (m_int-i), p, j, m, b, b))  * exp( *three_dim_index(X_device, (m_int-i), p, 0, m, b, num_assets) ) - stock_expectation ) , 2 );

}

beta=numerator/denominator;

ControlMeshContinVal= ContinVal-beta*(stock_expectation-true_stock_expectation);

return ControlMeshContinVal;
}

