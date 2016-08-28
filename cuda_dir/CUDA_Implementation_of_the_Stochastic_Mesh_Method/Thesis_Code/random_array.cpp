#include <random>
#include <iostream>

double* random_array(int N){

std::random_device rd;
std::default_random_engine generator;
generator.seed( rd() );
std::normal_distribution<double> distribution(0.0,1.0);

double* p; 
p = new double[N];
double z;
for(int i=0; i<N; i++ ){
z=distribution(generator);
p[i]=z;
}

return p;

}
