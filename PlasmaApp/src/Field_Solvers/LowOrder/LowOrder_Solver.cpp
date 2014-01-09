/*
 * LowOrder_Solver.cpp
 *
 *  Created on: Dec 11, 2013
 *      Author: payne
 */

#include "LowOrder_Solver.h"
#include "SysStateVector.h"
namespace HoLo
{
namespace LowOrder
{

void Solver::solve(MultiScale::SysStateVector** Q_in, HoLo::StateVector* gamma)
{
	// Set the latest problem states
	problem->Q_lo[0] -> CopyFrom(Q_in[0]->moments);
	problem->gamma -> CopyFrom(gamma);

	// Solve the nonlinear system
	nlsolver -> Solve();


}

void Solver::ComputeConsistancy(MultiScale::SysStateVector** Q_in,
		HoLo::StateVector* consTerm)
{
	problem->ComputeConsTerms(consTerm);
}

void Solver::ComputeResidual(MultiScale::SysStateVector** Q_in,
			HoLo::StateVector* resid)
{
	problem->ComputeHoResidual(resid);
}

} /* namespace LowOrder */
}
