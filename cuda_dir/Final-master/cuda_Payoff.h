#include <vector>
#include "enum_header.h" 
#ifndef PAYOFF2_H
#define PAYOFF2_H


class PayOff
{
public:
	__host__ __device__ PayOff(){};
	__host__ __device__ virtual double operator()( double* X, int i, int j, double m, int b, Containers container, int num_assets ) const=0;
	__host__ __device__ virtual ~PayOff(){}
private:
};


class GeometricPayOffCall : public PayOff
{
public:
	__host__ __device__ GeometricPayOffCall(double k);
	__host__ __device__ virtual double operator()(double* X, int i, int j, double m, int b, Containers container, int num_assets) const;
	__host__ __device__ virtual ~GeometricPayOffCall(){}
private:
	__host__ __device__ double Strike;
};


class GeometricPayOffPut : public PayOff
{
public:
	__host__ __device__ GeometricPayOffPut(double k);
	__host__ __device__ virtual double operator()(double* X, int i, int j, double m, int b, Containers container, int num_assets) const;
	
	__host__ __device__ virtual ~GeometricPayOffPut(){}
private:
	__host__ __device__ double Strike;
};
#endif
