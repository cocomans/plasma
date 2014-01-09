/*
 * Problem.cpp
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#include "Problem.h"
#include "StateVector.h"
#include "DiscretePDE.h"
#include "BoundaryValues.h"
#include "DiscretePDE1D3V.h"

namespace HoLo
{
namespace LowOrder
{

Problem::Problem ( SimParams* _params, Epetra_Comm* comm )
{
	simParams = _params;
	nSpatial = simParams->nSpatial;
	nVel = simParams->nVel;

	switch(nVel)
	{
	case 1:
		printf("Warning no pdes for this case\n");
		exit(1);
	case 2:
		printf("Warning no pdes for this case\n");
		exit(1);
	case 3:
		switch(nSpatial)
		{
		case 1:
			pde = new DiscretePDE1D3V();
			break;
		case 2:
			printf("Warning no pdes for this case\n");
			exit(1);
		case 3:
			printf("Warning no pdes for this case\n");
			exit(1);
		default:
			printf("Warning no pdes for this case\n");
			exit(1);
		}
		break;
	default:
		printf("Warning no pdes for this case\n");
		exit(1);
	}


}

Problem::~Problem()
{
	// TODO Auto-generated destructor stub



}

Epetra_Vector* Problem::GetX()
{
	return Q_lo[0]->Dep_vector;
}

void Problem::ComputeConsTerms()
{
	gamma->Set(0.0);

	pde->EvaluatePDE(Q_hi,stress,gamma,Res);

	gamma->CopyFrom(Res);

}

bool Problem::computeF( const Epetra_Vector &x,
	 Epetra_Vector &F,
	 const FillType fillFlag )
{
	// Copy stuf into StateVectors
	Res -> CopyDepFrom(&F);
	Q_lo[0] ->  CopyDepFrom(&x);

	// Apply boundary conditions
	Q_lo[0] -> ApplyBC(BCs);

	pde->EvaluatePDE(Q_lo,stress,gamma,Res);

	Res -> CopyDepTo(&F);

}

bool Problem::computePreconditioner( const Epetra_Vector &x,
			  Epetra_Operator &M,
			  Teuchos::ParameterList *precParams )
{

}

HoLo::LowOrder::Preconditioner* Problem::GetPreconditioner()
{
	return precond;
}


} /* namespace LowOrder */
} /* namespace HoLo */
