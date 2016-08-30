#include<iostream>



double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets){
int m_int = (int)m;
double* p;

p=&matrix[i*b*num_assets+j*num_assets+k];
return p;
}


double* two_dim_index(double* vector, int i, int j, double m, int b){
int m_int= (int)m;
double* p;

p=&vector[b*(i)+(j)];

return p;
}


