/*
 * DiscretePDE.h
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#ifndef DISCRETEPDE_H_
#define DISCRETEPDE_H_

#include "SimParams.h"

namespace HoLo
{
class StateVector;
class StressTensor;
namespace LowOrder
{

class DiscretePDE
{
public:
	int nx, ny, nz;
	real dt,qe,qi,me,mi,Te;



	virtual ~DiscretePDE();

	virtual void UpdateAuxillary(StateVector** Q);

	virtual void EvaluatePDE(StateVector** Q,
			StressTensor** stress,
			StateVector* Gamma,
			StateVector* Res);

	virtual double FieldResidual(StateVector* Res);
};

} /* namespace LowOrder */
} /* namespace HoLo */
#endif /* DISCRETEPDE_H_ */
