#ifndef AMPERE_SOLVER_H
#define AMPERE_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include "../FieldSolver.h"


class FieldData;
class ParallelInfo;
class HOMoments;


class AmpereSolver : public FieldSolver
{
public:

	~AmpereSolver();

	void solve(PlasmaData* pdata,FieldData* fields,HOMoments* moments);

	realkind calc_residual(PlasmaData* pdata,
			FieldData* fields_next,
			FieldData* fields_old,
			HOMoments* moments);

};






#endif /* AMPERE_SOLVER_H */
