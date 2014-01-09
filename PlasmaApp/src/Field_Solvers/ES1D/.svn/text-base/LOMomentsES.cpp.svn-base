//============================================================================================
//
//	AUTHOR: WILLIAM T. TAITANO
//	This is the source code for the LOMomentsES class
//
//============================================================================================
#include "LOMomentsES.h"

//============================================================================================
//	Constructor
//============================================================================================
LOMomentsES::LOMomentsES(PlasmaData* pDataIn, HOMoments* hoMoments)
{
	// copy pData pointer from pDataIn
	pData 	= pDataIn;
	// extract some grid parameters
	int nx 			= pData->nx;
	int ny 			= pData->ny;
	int nz 			= pData->nz;
	int nspecies	= pData->nspecies;
	int ntotal 		= nx*ny*nz*nspecies;
	// allocate the size of member variables
	charge 			= (realkind*)malloc(ntotal*sizeof(realkind));
	currentx 		= (realkind*)malloc(ntotal*sizeof(realkind));
	currenty 		= (realkind*)malloc(ntotal*sizeof(realkind));
	currentz 		= (realkind*)malloc(ntotal*sizeof(realkind));
	for (int s = 0;s < nspecies; s++) {
		for (int i = 0;i < nx;i++) {
			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < nz; k++) {
					this->get_val(i,j,k,s,LOMoments_charge) 	= hoMoments->get_val(i,j,k,s,HOMoments_charge);
					this->get_val(i,j,k,s,LOMoments_currentx)	= hoMoments->get_val(i,j,k,s,HOMoments_currentx);
					this->get_val(i,j,k,s,LOMoments_currenty) 	= hoMoments->get_val(i,j,k,s,HOMoments_currenty);
					this->get_val(i,j,k,s,LOMoments_currentz) 	= hoMoments->get_val(i,j,k,s,HOMoments_currentz);
				}
			}
		}
	}
}
//============================================================================================
//	Destructor
//============================================================================================
LOMomentsES::~LOMomentsES()
{
	free(charge);
	free(currentx);
	free(currenty);
	free(currentz);
}
//============================================================================================
//	Member Function: get_val function
//============================================================================================
realkind& LOMomentsES::get_val(	const int ix, const int iy, const int iz,
					const int ispecies,enum LOMoments_moment moment)
{
	realkind* result;

	int ix2,iy2,iz2;


	ix2 = ((ix%pData->nx)+pData->nx)%pData->nx;
	iy2 = ((iy%pData->ny)+pData->ny)%pData->ny;
	iz2 = ((iz%pData->nz)+pData->nz)%pData->nz;

	int iout = ix2 + pData->nx * (iy2 + pData->ny * (iz2 + pData->nz * ispecies));
	switch(moment)
	{
	case LOMoments_charge:
		result = charge + iout;
		break;
	case LOMoments_currentx:
		result = currentx + iout;
		break;
	case LOMoments_currenty:
		result = currenty + iout;
		break;
	case LOMoments_currentz:
		result = currentz + iout;
		break;
	default:
		break;
	}

	return *result;
}
