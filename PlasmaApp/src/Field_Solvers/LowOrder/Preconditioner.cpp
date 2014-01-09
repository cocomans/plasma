/*
 * Preconditioner.cpp
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#include "Preconditioner.h"
#include "MLSolver.h"

namespace HoLo
{
namespace LowOrder
{

int Preconditioner::SetUseTranspose(bool UseTranspose)
{
  return 1;
}

int Preconditioner::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return 1;
}


double Preconditioner::NormInf() const
{
  return 1;
}

const char* Preconditioner::Label() const
{
  return "User Preconditioner";
}

bool Preconditioner::UseTranspose() const
{
  return false;
}

bool Preconditioner::HasNormInf() const
{
  return false;
}

const Epetra_Comm& Preconditioner::Comm () const
{
	return * comm;
}
const Epetra_Map& Preconditioner::OperatorDomainMap () const
{
	return *red_map;
}
const Epetra_Map& Preconditioner::OperatorRangeMap () const
{
	return *red_map;
}

} /* namespace LowOrder */
} /* namespace HoLo */
