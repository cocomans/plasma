#ifndef PROBLEM_INITIALIZER_H
#define PROBLEM_INITIALIZER_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "PlasmaData.h"
#include <stdint.h>


class ParticleList;
class HOMoments;
class FieldData;
class ParallelInfo;

extern const realkind pi_const;

class ProblemInitializer
{
public:

	virtual ~ProblemInitializer(){};

	virtual void initialize_particles(ParticleList* particles, HOMoments* moments,ParallelInfo* myinfo){};

	virtual void initialize_fields(FieldData** fields, HOMoments** moments,ParallelInfo* myinfo){};

	virtual void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,int ispecies,int iptcl){};

	virtual void init_velocities(realkind& vx, realkind& vy, realkind& vz,
								int ispecies,int iptcl){};

	virtual void check_step(ParticleList* particles,
			HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo){}

	virtual void finish(ParticleList* particles,
			HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo){}

	PlasmaData* pdata;

	char* title;


};




extern double drandom();


__inline__
double box_muller(double mean = 0.0,double std = 1.0)
{
	return std*sqrt(-2.0*log(drandom()))*cos(2.0*3.14159265359*drandom()) + mean;
}

class raised_cos
{
public:

	raised_cos(float alpha_in,float k_in){alpha = alpha_in;k = k_in;}

	float operator()(float xp)
	{
		return xp*k/(2.0*pi_const) + alpha*sin(xp*k)/(2.0*pi_const);
	}


	float alpha,k;
};

class raised_sin
{
public:

	raised_sin(double alpha_in,double k_in){alpha = alpha_in;k = k_in;}

	double operator()(double xp)
	{
		return xp*k/(2.0*pi_const) - alpha*cos(xp*k)/(2.0*pi_const)+alpha/(k*2.0*pi_const);
	}


	float alpha,k;
};

class island_coalescence
{
public:
	island_coalescence(const double _epsilon,
					   const double _lambda,
					   const double _Lx,
					   const double _Ly,
					   const double _x0,
					   const double _y0,
					   const int 	_nx,
					   const int 	_ny) : epsilon(_epsilon),
	lambda(_lambda),nx(_nx),ny(_ny),Lx(_Lx),Ly(_Ly),x0(_x0),y0(_y0)
	{
		pdf = (double*)malloc((nx*ny+1)*sizeof(double));
		pdf[0] = 0;

		for(int j=0;j<ny;j++)
			for(int i=0;i<nx;i++)
			{
				double x,y;
				x = i*Lx/nx + x0;
				y = j*Ly/ny + y0;

				pdf[i+nx*j + 1] = pdf[i+nx*j]
				 + ((1.0-epsilon*epsilon)/pow((cosh(y/lambda)+epsilon*cos(x/lambda)),2.0));
			}

		for(int i=0;i<nx*ny+1;i++)
			pdf[i] /= pdf[nx*ny];
	}

	void getxy(double xp, float &x, float &y)
	{
		uint64_t* xp_i;
		uint hi,lo;

		// cast xp as a uint64
		xp_i = (uint64_t*)&xp;
		// split xp into a hi and lo component
		hi = (*xp_i >> 24) & 0x0FFF;
		lo = (*xp_i) & 0x00000FFF;

		uint scale = 0x0FFF;

		// x is the scaled lo component, y is the scaled hi component
		x = Lx*lo/((float)scale) + x0;
		y = Ly*hi/((float)scale) + y0;
	}

	double operator()(double xp)
	{
		int i = floor(xp);
		double frac = xp-i;

		return pdf[i]*(frac) + pdf[i+1]*(1-frac);
	}

	const double epsilon;
	const double lambda;
	const double x0,y0;
	const double Lx,Ly;
	double* pdf;
	const int nx,ny;


};

template<class Op>__inline__
double distribution_intrp(double pin,double x_low0,double x_high0,Op cdf)
{
	// solves pin = cdf(x) for x
	double tol = 1.0e-10*fabs(x_high0+x_low0);
	double diff,x_high,x_low,x_half;
	double result;
	int iter_max = 100;
	int iter =0;

	pin = fmodf(pin,1.0f);

	diff = 2.0*tol;

	x_low = x_low0;
	x_high = x_high0;


	while(fabs(diff) > tol)
	{
		x_half = 0.5*(x_low + x_high);

		diff = pin - cdf(x_half);

		if(fabs(diff) <= tol)
		{
			result = x_half;
			break;
		}
		else if(diff > 0)
		{
			x_low = x_half;
		}
		else
		{
			x_high = x_half;
		}

		iter++;

		if(iter > iter_max)
		{
			result = x_half;
			break;
		}

	}

	return result;


}


#endif /* PROBLEM_INITIALIZER_H */
