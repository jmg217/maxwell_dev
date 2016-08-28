#include <iostream>
#include <random>

double* random_array(int N);

int main(){
srand((unsigned)time(NULL));
std::random_device rd;
std::default_random_engine generator;
generator.seed( rd() );
std::normal_distribution<double> distribution(0.0,1.0);

double* p;

p=random_array(100);

for(int i=0; i<100; i++){
std::cout<<p[i]<<std::endl;
}

return 0;
}
