/*-------------------------------------------------------------------------*/
/**
	@file	PExitCheck.h
	@author	J. Payne
	@date		3/07/2013
	@brief	Declares the virtual PExitCheck class and its variants


*/
/*--------------------------------------------------------------------------*/
#ifndef PEXIT_CHECHER_H
#define PEXIT_CHECHER_H



#include "PlasmaData.h"


#ifdef DOUBLE_PRECISION
#define RK_EPSILON 1.0e-14
#else
#define RK_EPSILON 1.0e-6
#endif


class PExitCheck
{
public:
	__host__ __device__
	virtual ~PExitCheck(){};

	__host__ __device__
	virtual int check_particle(const realkind dt_finished,
			const int ix, const int iy, const int iz,
					   	   	    const int iter){}
	__host__ __device__
	virtual int check_particle(const realkind dt_finished,
								const int ix, const int iy,
								const int iter){}
	__host__ __device__
	virtual int check_particle(const realkind dt_finished,
			const int ix,
								const int iter){}

	realkind dt;

	int nSubcycle_max;


};


class PExitCheckCPU : public PExitCheck
{
public:
	__host__ __device__
	PExitCheckCPU(realkind _dt,int _nSubcycle_max){dt=_dt;nSubcycle_max=_nSubcycle_max;}

	__host__ __device__
	~PExitCheckCPU(){};

	__host__ __device__
	int check_particle(const realkind dt_finished,
			const int ix, const int iy, const int iz,
							    const int iter)
	{	//printf("dt check: %f - %f = %g \n",dt_finished,dt*(1.0-RK_EPSILON),dt_finished-dt*(1.0-RK_EPSILON));

		return  (dt_finished >= dt)||(iter >= nSubcycle_max);}
	__host__ __device__
	int check_particle(const realkind dt_finished,
			const int ix, const int iy,
								const int iter)
	{//printf("dt check: %f - %f = %g \n",dt_finished,dt*(1.0-RK_EPSILON),dt_finished-dt*(1.0-RK_EPSILON));
		return  (dt_finished >= (dt))||(iter >= nSubcycle_max);}
	__host__ __device__
	int check_particle(const realkind dt_finished,
			const int ix,
								const int iter)
	{//	printf("dt check: %f - %f = %g \n",dt_finished,dt*(1.0-RK_EPSILON),dt_finished-dt*(1.0-RK_EPSILON));

		return  (dt_finished >= (dt))||(iter >= nSubcycle_max);}

};

class PExitCheckGPU : public PExitCheck
{
public:
	__host__ __device__
	PExitCheckGPU():
	ix0(0),iy0(0),iz0(0),
	nx(0),ny(0),nz(0){};

	__host__ __device__
	PExitCheckGPU(realkind _dt,int _nSubcycle_max,
				int _ix0,int _iy0,int _iz0,
				int _nx,int _ny,int _nz):
				ix0(_ix0),iy0(_iy0),iz0(_iz0),
				nx(_nx),ny(_ny),nz(_nz){dt = _dt; nSubcycle_max = _nSubcycle_max;};
	__host__ __device__
	~PExitCheckGPU(){};

	__host__ __device__
	int check_particle(const realkind dt_finished,
			const int ix, const int iy, const int iz,
							    const int iter)
	{

		// Time step and iter count check
		bool qtime = (dt_finished >= dt)||(iter >= nSubcycle_max);

		// Position Check
		int xadj = ix - ix0;
		qtime = qtime || (xadj < 0) || (xadj > nx);
		xadj = iy - iy0;
		qtime = qtime || (xadj < 0) || (xadj > ny);
		xadj = iz - iz0;
		qtime = qtime || (xadj < 0) || (xadj > nz);

		return qtime;

	}

	__host__ __device__
	int check_particle(const realkind dt_finished,
			const int ix, const int iy,
								const int iter)
	{
	// Time step and iter count check
	bool qtime = (dt_finished >= dt)||(iter >= nSubcycle_max);

	// Position Check
	int xadj = ix - ix0;
	qtime = qtime || (xadj < 0) || (xadj > nx);
	xadj = iy - iy0;
	qtime = qtime || (xadj < 0) || (xadj > ny);

	return qtime;
	}

	__host__ __device__
	int check_particle(const realkind dt_finished,
								const int ix,
								const int iter)
	{return  (dt_finished >= dt)||(iter >= nSubcycle_max);}


	realkind dt;

	int nSubcycle_max;

	int ix0,iy0,iz0;
	char nx,ny,nz;

};














#endif
