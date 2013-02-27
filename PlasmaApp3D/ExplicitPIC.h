#ifndef EXPLICIT_PIC_H
#define EXPLICIT_PIC_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include <omp.h>


class PlasmaData;
class ProblemInitializer;
class ParticleList;
class FieldData;
class ParallelInfo;
class FieldSolver;
class HOMoments;

class ExplicitPIC
{
public:

	ExplicitPIC(PlasmaData* pdata_in, FieldSolver* LOsolver_in, ProblemInitializer* initializer_in,
			ParallelInfo** pinfo,int ndevices_in);

	ExplicitPIC(void){};

	void init(void);

	void run(void);

	void simulate(void);





	FieldData** fields;
	ParticleList** particles;
	HOMoments**	moments;
	HOMoments* moments_old;


	PlasmaData* pdata;
	ParallelInfo** pll_info;


	FieldSolver* LOsolver;
	ProblemInitializer* initializer;

	int myid_mpi;
	int ndevices;
	int nthreads;
	int n_nodes;
};











#endif /* EXPLICIT_PIC_H */
