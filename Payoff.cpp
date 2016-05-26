#include "Payoff.h"
#include <algorithm>
#include <vector>
#include <cmath>

GeometricPayOffCall::GeometricPayOffCall(double k) : Strike(k)
{
}

GeometricPayOffPut::GeometricPayOffPut(double k) : Strike(k)
{
}

double GeometricPayOffCall::operator () (const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const
{
double h;
h=1;

for(int l=0; l<asset_amount.size(); l++){
        h*=exp(X[i][j][l]);
}
h=pow(h,1.0/(asset_amount.size()));
h=std::max(h-Strike,0.0);

return h;

}

double GeometricPayOffPut::operator () (const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const
{
double h;
h=1;

for(int l=0; l<asset_amount.size(); l++){
        h*=exp(X[i][j][l]);
}
h=pow(h,1.0/(asset_amount.size()));
h=std::max(Strike-h,0.0);

return h;

}


