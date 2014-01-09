/*
 * DiscretePDE1D3V.h
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#ifndef DISCRETEPDE1D3V_H_
#define DISCRETEPDE1D3V_H_

#include "DiscretePDE.h"

namespace HoLo
{
namespace LowOrder
{

class DiscretePDE1D3V: public DiscretePDE
{
public:

	class Dep{public:enum {ne,pe_u,pe_v,pe_w,ni,pi_u,pi_v,pi_w,Av,Aw,phi,LF};};
	class Aux{public:enum {Eu,Ev,Ew,Bu,Bv,Bw,Au,LF};};
	class Prec{public:enum {ne,ni,LF};};

	// This must be defined for all DiscretePDE's
	// Negative signs are used to indicate auxiliary variables
	static const int pdeVar2genVar[SimParams::Vars::LF];

	void UpdateAuxillary(StateVector** Q);

	void EvaluatePDE(StateVector** Q,
				StressTensor** stress,
				StateVector* Gamma,
				StateVector* Res);

	double FieldResidual(StateVector* Res);


};



} /* namespace LowOrder */
} /* namespace HoLo */
#endif /* DISCRETEPDE1D3V_H_ */
