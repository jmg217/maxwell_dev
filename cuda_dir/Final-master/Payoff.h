#include <vector>
#include "enum_header.h" 
#ifndef PAYOFF2_H
#define PAYOFF2_H


class PayOff
{
public:
	PayOff(){};
	virtual double operator()( double* X, double asset_amount[], int i, int j, double m, int b, Containers container, int num_assets ) const=0;
	virtual ~PayOff(){}
private:
};


class GeometricPayOffCall : public PayOff
{
public:
	GeometricPayOffCall(double k);
	virtual double operator()(double* X, double asset_amount[], int i, int j, double m, int b, Containers container, int num_assets) const;
	virtual ~GeometricPayOffCall(){}
private:
	double Strike;
};


class GeometricPayOffPut : public PayOff
{
public:
	GeometricPayOffPut(double k);
	virtual double operator()(double* X, double asset_amount[], int i, int j, double m, int b, Containers container, int num_assets) const;
	
	virtual ~GeometricPayOffPut(){}
private:
	double Strike;
};
#endif
