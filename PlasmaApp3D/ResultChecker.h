#ifndef FIELD_SOLVER_H
#define FIELD_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include <omp.h>


class PlasmaData;
class FieldData;
class ParallelInfo;
class HOMoments;


class FieldSolver
{
public:

	virtual ~FieldSolver(){};

	virtual void solve(PlasmaData* pdata,FieldData* fields,HOMoments* moments){};
};


#endif /* FIELD_SOLVER_H */
