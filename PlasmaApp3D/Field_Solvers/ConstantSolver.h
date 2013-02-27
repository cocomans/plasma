#ifndef CONSTANT_SOLVER_H
#define CONSTANT_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include "../FieldSolver.h"

class PlasmaData;
class FieldData;
class ParallelInfo;
class HOMoments;


class ConstantSolver : public FieldSolver
{
public:

	~ConstantSolver();

	void solve(PlasmaData* pdata,FieldData* fields,HOMoments* moments);
};






#endif /* CONSTANT_SOLVER_H */
