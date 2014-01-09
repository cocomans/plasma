/*
 * SysStateVector.h
 *
 *  Created on: Dec 19, 2013
 *      Author: payne
 */

#ifndef SYSSTATEVECTOR_H_
#define SYSSTATEVECTOR_H_
#include "SimParams.h"

namespace HoLo
{
	class StateVector;
	namespace HighOrder
	{
	class NodeParticleList;
	}
}

namespace MultiScale
{

class SysStateVector
{
public:
	HoLo::StateVector* moments;
	HoLo::HighOrder::NodeParticleList* particles;
	HoLo::StressTensor* stress;

	SysStateVector(SimParams* _params);

	void CopyDepFrom(const Epetra_Vector* _src);

	void CopyAuxFrom(const Epetra_Vector* _src);

	void CopyDepTo(Epetra_Vector* _src);

	void CopyAuxTo(Epetra_Vector* _src);

	void CopyFrom(StateVector* _src);

	void CopyFrom(SysStateVector* _src);


	void Set(realkind val);

};

} /* namespace MultiScale */
#endif /* SYSSTATEVECTOR_H_ */
