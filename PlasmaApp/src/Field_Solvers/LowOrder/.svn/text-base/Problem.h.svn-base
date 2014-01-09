/*
 * Problem.h
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#ifndef PROBLEM_H_
#define PROBLEM_H_

#include "SimParams.h"

namespace HoLo
{
class StateVector;
class StressTensor;
namespace LowOrder
{
class BoundaryValues;
class DiscretePDE;
class Preconditioner;
class Problem :
		public NOX::Epetra::Interface::Required,
		public NOX::Epetra::Interface::Preconditioner
{
public:

	// Basically tied to the SysStateVectors
	StateVector** Q_hi;
	// Tied to Q_lo in LowOrder::Solver
	StateVector** Q_lo;
	StateVector*  gamma;
	// Basically tied to the SysStateVectors
	StressTensor** stress;
	StateVector* Res;


	Problem (SimParams* _params, Epetra_Comm* comm );


	~Problem (  );

	/// Returns Q_lo[0]->Dep
	Epetra_Vector* GetX();

	void ComputeConsTerms(HoLo::StateVector* consTerm);

	void ComputeHoResidual(HoLo::StateVector* Resid);


	// inherited functions
	bool computeF( const Epetra_Vector &x,
		 Epetra_Vector &F,
		 const FillType fillFlag );

	bool computePreconditioner( const Epetra_Vector &x,
				  Epetra_Operator &M,
				  Teuchos::ParameterList *precParams );

	HoLo::LowOrder::Preconditioner* GetPreconditioner();

	Epetra_Comm*             comm;
	SimParams* simParams;
	BoundaryValues* BCs;
	int nSpatial;
	int nVel;
	int nStates;


	DiscretePDE*   pde;
	HoLo::LowOrder::Preconditioner* precond;


};

} /* namespace LowOrder */
} /* namespace HoLo */
#endif /* PROBLEM_H_ */
