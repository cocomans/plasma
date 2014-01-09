/*
 * LowOrder_Solver.h
 *
 *  Created on: Dec 11, 2013
 *      Author: payne
 */

#ifndef LOWORDER_SOLVER_H_
#define LOWORDER_SOLVER_H_
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "NOX_Epetra.H"


class StateVector;
namespace MultiScale
{class SysStateVector;}

namespace HoLo
{
namespace LowOrder
{
class Problem;
class NonLinearSolver;

class Solver
{
public:

	SimParams* params;
	Problem* problem;
	NonLinearSolver* nlsolver;

	Solver(SimParams* _params);

	void solve(MultiScale::SysStateVector** Q_in, HoLo::StateVector* gamma);

	void ComputeConsistancy(MultiScale::SysStateVector** Q_in,
			HoLo::StateVector* consTerm);

	void ComputeResidual(MultiScale::SysStateVector** Q_in,
				HoLo::StateVector* resid);



};



} /* namespace LowOrder */
}
#endif /* LOWORDER_SOLVER_H_ */
