#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <sstream>


double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t);

double* three_dim_index(double* matrix, int i, int j, int k, double m, int b, int num_assets);
double* two_dim_index(double* vector, int i, int j, double m, int b);

double kahansum(double* sortvector, int b){
double sum=0, c=0, y, t;

        for(int i=0; i<b; i++){
                y=sortvector[i]-c;
                t=sum+y;
                c=(t-sum)-y;
                sum=t;
        }

return sum;
}


void meshweights(double* W, double m, int b, double sigma[], double delta[], double r, double delta_t, double* X, int num_assets, double* weight_denominator){
double wdenominator, w;

double* sortvector;
sortvector=new double[b];

for(int I=0; I<m; I++){
	
	if(I==0){
		for(int k=0; k<b; k++){
        	
			for(int j=0; j<b; j++){

				if(j==0){
					
					*three_dim_index(W, I, k, j, m, b, b)=1;
				}// all weights from the starting node are equal to 1

				else{
					
					*three_dim_index(W, I, k, j, m, b, b)=0;
				}
			}

		}
	}


	if(I>0){

		for(int k=0; k<b; k++){	
		//dim1temp.clear();
		//sortvector.clear();

		wdenominator=0;

			for(int j=0; j<b; j++){
			//std::cout<<j<<std::endl;	
				w=1;
				//w=0; //set w to 1 since it will be equal to a product
				for(int jj=0; jj<num_assets; jj++){
				w = w * density(*three_dim_index(X, (I-1), j, jj, m, b, num_assets), *three_dim_index(X, I, k, jj, m, b, num_assets), sigma[jj], r, delta[jj], delta_t);

				}
			//w = exp(w);
			//dim1temp.push_back(w);
			//sortvector.push_back(w);   
			sortvector[j]=w;                                                                
			}
			//std::sort(sortvector.begin(), sortvector.end());
			wdenominator=kahansum(sortvector, b);
				//devide each element by the denominator
//			std::cout<<"before"<<std::endl;
			*two_dim_index(weight_denominator, I-1, k, m-1, b)=wdenominator;
//			std::cout<<"after and I= "<<I<<std::endl;
//			std::cout<<*two_dim_index(weight_denominator, I, k, m-1, b)<<std::endl;
			for(int t=0; t<b; t++){
				*three_dim_index(W, (I), k, t, m, b, b)=(((double)b)*(sortvector[t]))/wdenominator;
			}
			//std::cout<<"after"<<std::endl;
		}
	//std::cout<<"outside 1"<<std::endl;	
	
	}

//std::cout<<"outside 2"<<std::endl;
}
//int W_N=b*(m-1);
//for(int test=0; test<W_N; test++){
//printf("weights=%f\n",weight_denominator[test]);
//}

//std::cout<<"outside 3"<<std::endl;
delete[] sortvector;
//std::cout<<"end"<<std::endl;
}
