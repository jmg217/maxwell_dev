#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include "enum_header.h"
#include <unistd.h>
#include <stdio.h>

#define PI 3.14159265358979323846


/* we need these includes for CUDA's random number stuff */

double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t)// This function returns the transition density between node values.
{
double f=0, x=0;
//x=(1/(sigma*sqrt(delta_t)))*(log(Xnew)-log(Xold)-(r-delta-0.5*sigma*sigma)*delta_t);
x=(1/(sigma*sqrt(delta_t)))*(Xnew-Xold-(r-delta-0.5*sigma*sigma)*delta_t);
//f= (1/(sigma*sqrt(delta_t)*Xnew))*(1/(sqrt(2*PI)))*exp(-0.5*x*x); // this is the transition density
f= (1/(sigma*sqrt(delta_t)))*(1/(sqrt(2*PI)))*exp(-0.5*x*x);
return f;
}

