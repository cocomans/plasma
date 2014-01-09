//================================================================================================
//
//	AUTHOR: WILLIAM T. TAITANO
//	This is the source code for the ConsistencyTerm class
//
//================================================================================================
#include "ConsistencyTermES.h"

//=======================================================================
//	Constructor
//=======================================================================

ConsistencyTermES::ConsistencyTermES(PlasmaData* pDataIn)
{
	// Assign pData pointer to pDataIn
	pData 			= pDataIn;
	// Extract some input parameters
	int nx 			= pData->nx;
	int ny 			= pData->ny;
	int nz 			= pData->nz;
	int nspecies 	= pData->nspecies;
	int ntotal 		= nx*ny*nz*nspecies;
	// Allocate the array size
	gammaN 			= (realkind*)malloc(ntotal*sizeof(realkind));
	gammaNU 		= (realkind*)malloc(ntotal*sizeof(realkind));
}
//=======================================================================
//	Destructor
//=======================================================================
ConsistencyTermES::~ConsistencyTermES()
{
	// Free the allocated memory for the consistency term
	free(gammaN);
	free(gammaNU);
	free(pData);
}
//=======================================================================
//	Member Function: get_val
//=======================================================================
realkind& ConsistencyTermES::get_val(const int ix, const int iy, const int iz,
		const int ispecies,enum consistencyTerm consistencyterm)
{
	realkind* result;

	int ix2,iy2,iz2;


	ix2 = ((ix%pData->nx)+pData->nx)%pData->nx;
	iy2 = ((iy%pData->ny)+pData->ny)%pData->ny;
	iz2 = ((iz%pData->nz)+pData->nz)%pData->nz;

	int iout = ix2 + pData->nx * (iy2 + pData->ny * (iz2 + pData->nz * ispecies));
	switch(consistencyterm)
	{
	case consistencyterm_continuity:
		result = gammaN + iout;
		break;
	case consistencyterm_momentum:
		result = gammaNU + iout;
		break;
	default:
		break;
	}

	return *result;
}
//=======================================================================
//	Member Function: consistencyCalcContinuity
//=======================================================================
void ConsistencyTermES::consistencyCalcContinuity(	int species, HOMoments* curHOMoments,
													HOMoments* oldHOMoments)
{
	// Extract some input parameters
	int nx 			= pData->nx;
	int ny 			= pData->ny;
	int nz 			= pData->nz;
	realkind dt 	= pData->dt;
	realkind dx		= pData->dxdi;
	realkind dy 	= pData->dydi;
	realkind dz 	= pData->dzdi;
	// Declare and calculate some quantities
	realkind dtRecip= 1.0/dt;
	realkind dxRecip= 1.0/dx;
	realkind dyRecip= 1.0/dy;
	realkind dzRecip= 1.0/dz;
	realkind nNew, nOld, nuHalf_ip1, nuHalf_i;
	realkind m 		= pData->mspecies[species];
	realkind tempVal= 0.;
	// Loop through cell to assign consistency term
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				// Extract the relevant moments at the current cell index
				nNew 	= curHOMoments->get_val(i,j,k,species,HOMoments_charge);
				nOld 	= oldHOMoments->get_val(i,j,k,species,HOMoments_charge);
				if (i == nx - 1) {
					nuHalf_i= curHOMoments->get_val(i,j,k,species,HOMoments_currentx);
					nuHalf_ip1= curHOMoments->get_val(0,j,k,species,HOMoments_currentx);
				}
				else {
					nuHalf_i= curHOMoments->get_val(i,j,k,species,HOMoments_currentx);
					nuHalf_ip1= curHOMoments->get_val(i+1,j,k,species,HOMoments_currentx);
				}
				// Calculate consistency term
				tempVal = 	m*dtRecip*(nNew - nOld) +
							m*dxRecip*(nuHalf_ip1 - nuHalf_i);
				tempVal /= 	nNew;		// scaling by dividing by the number density.
										// we may want to be careful for very large scale simualtion
										// when we have zero particles. Perhaps we can normalize by the
										// half time density to assure an increased possibility in avoiding
										// divide by zeros.
				// Assign consistency term to relevant cell index
				this->get_val(i,j,k,species,consistencyterm_continuity) = tempVal;
			}
		}
	}
}
//=======================================================================
//	Member Function: consistencyCalcMomentum
//=======================================================================
void ConsistencyTermES::consistencyCalcMomentum(	int species, HOMoments* curHOMoments,
													HOMoments* oldHOMoments, FieldData* curFieldData,
													FieldData* oldFieldData)
{
	// Extract some input parameters
	int nx 			= pData->nx;
	int ny 			= pData->ny;
	int nz 			= pData->nz;
	realkind dt 	= pData->dt;
	realkind dx		= pData->dxdi;
	realkind dy 	= pData->dydi;
	realkind dz 	= pData->dzdi;
	// Declare and calculate some quantities
	realkind dtRecip= 1.0/dt;
	realkind dxRecip= 1.0/dx;

	realkind 	nuHalf_i, nuHalfOld_i,
				nNew_i, nNew_im1, nOld_i, nOld_im1,
				S2xxNew_i, S2xxOld_i, S2xxNew_im1, S2xxOld_im1,
				ms, qs;
	ms 			= pData->mspecies[species];
	qs 			= qe*pData->qspecies[species];
	realkind 	nFaceHalf, nFaceNew, nFaceOld;
	realkind 	Ehalf;
	realkind qms 	= qs/ms;
	realkind tempVal_x = 0.;
	// Loop through cell to assign consistency term
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				// Initialize tempVal_x to 0
				tempVal_x 	= 0.;
				// Extract the relevant moments at the current cell index
				nuHalf_i	= curHOMoments->get_val(i,j,k,species,HOMoments_currentx);
				nuHalfOld_i = oldHOMoments->get_val(i,j,k,species,HOMoments_currentx);

				if(i == 0 || i == nx) {
					nNew_i 	= curHOMoments->get_val(0,j,k,species,HOMoments_charge);
					nOld_i 	= oldHOMoments->get_val(0,j,k,species,HOMoments_charge);
					nNew_im1= curHOMoments->get_val(nx-1,j,k,species,HOMoments_charge);
					nOld_im1= oldHOMoments->get_val(nx-1,j,k,species,HOMoments_charge);

					S2xxNew_i 	= curHOMoments->get_val(0,j,k,species,HOMoments_S2xx);
					S2xxOld_i 	= oldHOMoments->get_val(0,j,k,species,HOMoments_S2xx);
					S2xxNew_im1 = curHOMoments->get_val(nx-1,j,k,species,HOMoments_S2xx);
					S2xxOld_im1 = oldHOMoments->get_val(nx-1,j,k,species,HOMoments_S2xx);

				}
				else {
					nNew_i 		= curHOMoments->get_val(i,j,k,species,HOMoments_charge);
					nOld_i 		= oldHOMoments->get_val(i,j,k,species,HOMoments_charge);
					nNew_im1	= curHOMoments->get_val(i-1,j,k,species,HOMoments_charge);
					nOld_im1	= oldHOMoments->get_val(i-1,j,k,species,HOMoments_charge);

					S2xxNew_i 	= curHOMoments->get_val(i,j,k,species,HOMoments_S2xx);
					S2xxOld_i 	= oldHOMoments->get_val(i,j,k,species,HOMoments_S2xx);
					S2xxNew_im1 = curHOMoments->get_val(i-1,j,k,species,HOMoments_S2xx);
					S2xxOld_im1 = oldHOMoments->get_val(i-1,j,k,species,HOMoments_S2xx);
				}
				nFaceNew	= 	0.5*(nNew_i + nNew_im1);
				nFaceOld 	= 	0.5*(nOld_i + nOld_im1);
				nFaceHalf 	= 	0.5*(nFaceNew + nFaceOld);
				Ehalf 		= 	0.5*(curFieldData->getE(i,j,k,0) + oldFieldData->getE(i,j,k,0));
				// Calculate the consistency term
				tempVal_x 	= 	ms*dtRecip*(nuHalf_i - nuHalfOld_i) +
								ms*0.5*dxRecip*(S2xxNew_i + S2xxOld_i - S2xxNew_im1 - S2xxOld_im1) -
								qs*nFaceHalf*Ehalf;
				tempVal_x 	/= 	nFaceNew;
				// Store the consistency term
				this->get_val(i,j,k,species,consistencyterm_momentum) = tempVal_x;

			}
		}
	}
}
