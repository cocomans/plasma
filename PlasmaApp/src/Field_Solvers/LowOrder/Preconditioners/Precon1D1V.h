/*
 * Precon1D1V.h
 *
 *  Created on: Dec 18, 2013
 *      Author: payne
 */

#ifndef PRECON1D1V_H_
#define PRECON1D1V_H_

#include "Preconditioner.h"

namespace HoLo
{
namespace LowOrder
{

class Precon1D1V: public HoLo::LowOrder::Preconditioner
{
public:
	~Precon1D1V();

	Precon1D1V(SimParams* _params,
			DeviceNS::CPU::Mesh* _mesh,
				Epetra_Comm* _comm);

	void UpdateMatrix(const Epetra_Vector &xg);
	int ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const;;

};

} /* namespace LowOrder */
} /* namespace HoLo */
#endif /* PRECON1D1V_H_ */
