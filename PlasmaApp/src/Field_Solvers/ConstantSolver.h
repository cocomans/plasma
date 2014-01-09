#ifndef CONSTANT_SOLVER_H
#define CONSTANT_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "../cuda_defines.h"
#include <gnuplot_i.h>
#include <mpi.h>
#include "../FieldSolver.h"

class PlasmaData;
class FieldDataCPU;
class ParallelInfo;
class NodeHOMoments;


class ConstantSolver : public FieldSolver
{
public:

	~ConstantSolver();

	void solve(PlasmaData* pdata,
			FieldDataCPU* fields, //output
			NodeHOMoments* moments);
};






#endif /* CONSTANT_SOLVER_H */
