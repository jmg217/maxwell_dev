#include "Payoff.h"
#include <algorithm>
#include <vector>
#include <cmath>

MeshPayOffCall::MeshPayOffCall(double k) : Strike(k)
{
}

MeshPayOffPut::MeshPayOffPut(double k) : Strike(k)
{
}

double MeshPayOffCall::operator () (const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const
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

double MeshPayOffPut::operator () (const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const
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


double PathPayOffPut::operator () (const std::vector<std::vector<double> >& S, std::vector<double>& asset_amount, int i) const
{
double h;

h=1;

for(int l=0; l<asset_amount.size(); l++){
         h*=exp(S[i][l]);
}
h=pow(h,1.0/(asset_amount.size()));
h=std::max(Strike-h,0.0);
return h;

}

double PathPayOffCall::operator () (const std::vector<std::vector<double> >& S, std::vector<double>& asset_amount, int i) const
{
double h;

h=1;

for(int l=0; l<asset_amount.size(); l++){
         h*=exp(S[i][l]);
}
h=pow(h,1.0/(asset_amount.size()));
h=std::max(h-Strike,0.0);
return h;

}


PathPayOffPut::PathPayOffPut(double k) : Strike(k)
{
}

PathPayOffCall::PathPayOffCall(double k) : Strike(k)
{
}
