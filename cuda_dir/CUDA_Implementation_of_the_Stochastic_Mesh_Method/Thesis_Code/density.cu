#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#define PI 3.14159265358979323846

__device__ double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t){

double f=0, x=0;
//x=(1/(sigma*sqrt(delta_t)))*(log(Xnew)-log(Xold)-(r-delta-0.5*sigma*sigma)*delta_t);
x=(1/(sigma*sqrt(delta_t)))*(Xnew-Xold-(r-delta-0.5*sigma*sigma)*delta_t);
//f= (1/(sigma*sqrt(delta_t)*Xnew))*(1/(sqrt(2*PI)))*exp(-0.5*x*x); // this is the transition density
f= (1/(sigma*sqrt(delta_t)))*(1/(sqrt(2*PI)))*exp(-0.5*x*x);
return f;
}

