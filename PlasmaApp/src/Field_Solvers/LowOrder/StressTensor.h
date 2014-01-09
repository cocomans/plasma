/*
 * StressTensor.h
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#ifndef STRESSTENSOR_H_
#define STRESSTENSOR_H_
#include "SimParams.h"

namespace DeviceNS
{
namespace CPU
{
class Mesh;
}
}

namespace HoLo
{


class StressTensor
{
public:
	enum {Suu,Suv,Suw,Svv,Svw,Sww,LF};
	enum {electron, ion};


	void Allocate(SimParams* params);

	void Allocate(SimParams* params,StressTensor* _in);

	void PartialAllocate(SimParams* params);



	double& Get(const int& iType,const int& iSpecies,args_ijk);

	void InitFiller();


	int nVec;
	int nx,ny,nz;
	int nBufx,nBufy,nBufz;
	int nTx, nTy, nTz;

	long long int nElements;
	SimParams* params;
	DeviceNS::CPU::Mesh* mesh;
	Epetra_BlockMap* map;
	Epetra_Vector* vec;
};

} /* namespace HoLo */
#endif /* STRESSTENSOR_H_ */
