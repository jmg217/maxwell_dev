#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "Payoff.h"

////////Remeber to use delete when using new

void one_dim_array(std::vector< double > &vals, double array[], int N)
{
//double* temp;
//temp =new double[N];

for(int i=0; i<N; i++){
array[i]=vals[i];

}

//return temp;

}


double** two_dim_array(std::vector< std::vector<double> > &vals, int N, int M)
{

   double** temp;
   temp = new double*[N];


   for(int i=0; (i < N); i++)
   { 
      temp[i] = new double[M];
      for(int j=0; (j < M); j++)
      {
          temp[i][j] = vals[i][j];
      } 
   }

return temp;
 }


double*** three_dim_array(std::vector<std::vector< std::vector<double> > > &vals, int N, int M, int L)
{


double*** temp;
temp =new double**[N];


for(int i=0; (i < N); i++)
   {
      temp[i] = new double*[M];
      for(int j=0; (j < M); j++)
      {
	temp[i][j]= new double[L];
	for(int k=0; k<L; k++){
          	temp[i][j][k] = vals[i][j][k];
      		}
	}
   }

return temp;
}

