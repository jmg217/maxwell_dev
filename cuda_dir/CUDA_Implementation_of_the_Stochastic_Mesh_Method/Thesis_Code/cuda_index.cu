#include<iostream>

//think about inlining this

__device__ double* three_dim_indexGPU(double* matrix, int i, int j, int k, double m, int b){
int m_int = (int)m;
double* p;

//specify index layout here
p=&matrix[(m_int)*b*(k)+(m_int)*(j)+(i)];

return p;
}


__device__ double* two_dim_indexGPU(double* vector, int i, int j, double m, int b){
//int m_int= (int)m;
double* p;

//specify index layout here
p=&vector[b*(i)+(j)];

return p;
}


