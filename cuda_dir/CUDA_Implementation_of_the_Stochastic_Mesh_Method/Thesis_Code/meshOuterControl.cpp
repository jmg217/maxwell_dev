#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>

double outer_control_mesh(std::vector<double>& Vvector, std::vector<double>& EuroVals, int N, int b, double EuroGeoOption ){

double ControlMeshVal=0, sum=0;

double beta=0;
double numerator, denominator, Average_V, EstEuroGeoOption;


for(int f=0;f<N;f++){
sum+=Vvector[f];
}

Average_V=(1/(double)N)*sum;
sum=0;

for(int l=0; l<N; l++){
sum+=EuroVals[l];
//std::cout<<EuroVals[l]<<std::endl;
}

EstEuroGeoOption=(1/(double)N)*sum;
sum=0;

////////////////YOUR HERE. GET BETA NEXT
for(int w=0; w<N; w++){
sum+=(EuroVals[w]-EstEuroGeoOption)*(Vvector[w]-Average_V);

}
numerator=sum;
sum=0;
for(int e=0; e<N; e++){
sum+=(EuroVals[e]-EstEuroGeoOption)*(EuroVals[e]-EstEuroGeoOption);

}
denominator=sum;
sum=0;

beta=numerator/denominator;

ControlMeshVal=Average_V-beta*(EstEuroGeoOption-EuroGeoOption);

return ControlMeshVal;
}

