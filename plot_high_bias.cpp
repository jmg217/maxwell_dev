#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

void print_high_payoff(int b, double m, std::vector<std::vector< std::vector<double> > >& X, std::vector< std::vector<double> >& V, std::vector<double>& asset_amount, std::vector<std::vector< std::vector<double> > >& W ){
double x=0;
std::ofstream outFile("highpayoff.txt", std::ios_base::app | std::ios_base::out);

for (int t=0; t<m; t++){

	for (int i=0; i<b; i++){
	x=0;
		for(int tt=0; tt<asset_amount.size(); tt++){
		x+=asset_amount[tt]*X[m-1-t][i][tt];
		}
	outFile << m-t <<"\t"<< exp(x) <<"\t"<< V[t][i]<<"\t"<<X[m-1-t][i][0]<<"\t"<<X[m-1-t][i][1]<< std::endl;
/*	if(t>0){
	for(int ttt=0; ttt<b; ttt++){
	outFile <<ttt<<"\t"<<V[t-1][ttt]<<"\t"<<W[m-t][ttt][i]<<std::endl;
	}
	}*/
	}
}

outFile.close();

}


