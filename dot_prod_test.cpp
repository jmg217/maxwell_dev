#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>

double inner_product(int N, std::vector<double>& first_vector, std::vector<double>& second_vector);

int main(){

int N=10000000;

std::vector<double> first_vector(N);
std::vector<double> second_vector(N);

 for(int i=0 ; i < N ; i++)
  {
    
    //x_host[i]=first_vector[i];
    //y_host[i]=second_vector[i];
    first_vector[i] = sin(i*0.013);
    second_vector[i] = cos(i*0.019);
  }

double dot=inner_product(N, first_vector, second_vector);
std::cout<<dot<<std::endl;

return 0;
}
