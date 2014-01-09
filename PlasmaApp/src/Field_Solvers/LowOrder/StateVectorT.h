/*
 * StateVectorT.h
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#ifndef STATEVECTORT_H_
#define STATEVECTORT_H_

#include "StateVector.h"
#include "DiscretePDE1D3V.h"

class SimParams;
namespace HoLo
{


template<class PDE>
class StateVector_T: public HoLo::StateVector
{

	StateVector_T(SimParams* _params);
	StateVector_T(StateVector vec_in);

	~StateVector_T();

	void Allocate(DeviceNS::CPU::Mesh* _mesh);

	void InitFiller();

	real& GetDep(const int iType,const int& i,const int& j,const int& k);


	real& GetAux(const int iType,const int& i,const int& j,const int& k);

	real& Get(const int iType,const int& i,const int& j=0,const int& k=0);

	void ApplyBC(HoLo::LowOrder::BoundaryValues* BCs);

	void CopyDepFrom(const Epetra_Vector* _src);

	void CopyAuxFrom(const Epetra_Vector* _src);

	void CopyDepTo(Epetra_Vector* _src);

	void CopyAuxTo(Epetra_Vector* _src);

	void CopyFrom(StateVector* _src);

};

typedef StateVector_T<LowOrder::DiscretePDE1D3V> StateVector1D3V;

HoLo::StateVector* StateVectorGen(SimParams* params)
{
	HoLo::StateVector* result = NULL;

	switch(params->nSpatial)
	{
	case 1:
		switch(params->nVelocity)
		{
		case 3:
			result = new StateVector1D3V(params);
			break;
		default:
			break;
		}
		break;
		default:
			break;
	}

	if(result == NULL)
	{
		printf("Error no matching StateVector for specified case\n");
		exit(0);
	}
}



} /* namespace HoLo */
#endif /* STATEVECTORT_H_ */
