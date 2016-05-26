#include <vector>
#ifndef PAYOFF2_H
#define PAYOFF2_H


class PayOff
{
public:
	PayOff(){};
	virtual double operator()(const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const=0;
	virtual ~PayOff(){}
private:
};


class GeometricPayOffCall : public PayOff
{
public:
	GeometricPayOffCall(double k);
	virtual double operator()(const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const;
	virtual ~GeometricPayOffCall(){}
private:
	double Strike;
};


class GeometricPayOffPut : public PayOff
{
public:
	GeometricPayOffPut(double k);
	virtual double operator()(const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const;
	
	virtual ~GeometricPayOffPut(){}
private:
	double Strike;
};
#endif
