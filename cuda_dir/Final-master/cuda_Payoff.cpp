#include "Payoff.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include "enum_header.h"


__host__ __device__ double* three_dim_index(double* matrix, int i, int j, int k, double m, int b);


__host__ __device__ double* two_dim_index(double* vector, int i, int j, double m, int b);


__host__ __device__ GeometricPayOffCall::GeometricPayOffCall(double k) : Strike(k)
{
}

__host__ __device__ GeometricPayOffPut::GeometricPayOffPut(double k) : Strike(k)
{
}

__host__ __device__ double GeometricPayOffCall::operator () (double* X, int i, int j, double m, int b, Containers container, int num_assets ) const
{
double h;
h=1;


if(container==matrix){
	for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);
		h*= exp(*three_dim_index(X, i, j, l, m, b));
	}
	h=pow(h,1.0/(num_assets));
//h=std::max(h-Strike,0.0);
	if(h-Strike>0){
		h=h-Strike;
	}
	else{
		h=0;
	}

}

else if(container==vector){
	for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);
        	h*= exp(*two_dim_index(X, i, l, m, b));
	}
	h=pow(h,1.0/(num_assets));
	if(h-Strike>0){
		h=h-Strike;
	}
	else{
		h=0;
	}

}

else{
printf("%s \n", "The container you choose when calling the payoff functor was incorrect");
return 0;
}

return h;

}

__host__ __device__ double GeometricPayOffPut::operator () (double* X, int i, int j, double m, int b, Containers container, int num_assets) const
{
double h;
h=1;
if(container==matrix){
	for(int l=0; l<num_assets; l++){
        //h*=exp(X[i][j][l]);
		h*=exp(*three_dim_index(X, i, j, l, m, b));
	}
	h=pow(h,1.0/(num_assets));
//h=std::max(Strike-h,0.0);
	if(Strike-h>0){
		h=Strike-h;
	}
	else{
		h=0;
	}

}

else if(container==vector){
	for(int l=0; l<num_assets; l++){
       // h*=exp(X[i][j][l]);
        	h*= exp(*two_dim_index(X, i, l, m, b));
	}
	h=pow(h,1.0/(num_assets));
	if(Strike-h>0){
		h=Strike-h;
	}
	else{
		h=0;
	}
}

else{
printf("%s \n", "The container you choose when calling the payoff functor was incorrect");
return 0;
}



return h;

}


