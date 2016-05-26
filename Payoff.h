#include <vector>
#ifndef PAYOFF2_H
#define PAYOFF2_H


class PayOff
{
public:
	PayOff(){};
	virtual double operator()(double Spot) const=0;
	virtual ~PayOff(){}
private:
};


class MeshPayOffCall : public PayOff
{
public:
	MeshPayOffCall(double k);
	virtual double operator()(const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const;
	virtual ~MeshPayOffCall(){}
private:
	double Strike;
};


class MeshPayOffPut : public PayOff
{
public:
	MeshPayOffPut(double k);
	virtual double operator()(const std::vector<std::vector< std::vector<double> > >& X, std::vector<double>& asset_amount, int i, int j) const;
	
	virtual ~MeshPayOffPut(){}
private:
	double Strike;
};

class PathPayOffCall : public PayOff
{
public:
        PathPayOffCall(double k);
        virtual double operator()(const std::vector<std::vector<double> >& S, std::vector<double>& asset_amount, int i) const;
        virtual ~PathPayOffCall(){}
private:
        double Strike;
};


class PathPayOffPut : public PayOff
{
public:
        PathPayOffPut(double k);
        virtual double operator()(const std::vector<std::vector<double> >& S, std::vector<double>& asset_amount, int i) const;
        virtual ~PathPayOffPut(){}
private:
        double Strike;
};

#endif
