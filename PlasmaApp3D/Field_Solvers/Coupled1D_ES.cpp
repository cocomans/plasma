#include "Coupled1D_ES.h"
#include "../FieldData.h"
#include "../HOMoments.h"
#include "LOSolverSI.h"
#include "NLResidual.h"


void Coupled1D_ES::init(PlasmaData* pDataIn,
						FieldData* fields_old_in,
						FieldData* fields_half_in,
						FieldData* fields_next_in,
						HOMoments* moments_old_in,
						HOMoments* moments_next_in)
{

	// Setup framework members members
	pData 			= pDataIn;
	fields_old 		= fields_old_in;
	fields_half 	= fields_half_in;
	fields_next 	= fields_next_in;
	moments_old 	= moments_old_in;
	moments_next 	= moments_next_in;

	solver 			= new LOSolverSI(	pData,
										moments_next_in, moments_old_in,
										fields_next_in, fields_old_in);
	residual 		= new NLResidual(pData);

}

Coupled1D_ES::~Coupled1D_ES()
{

}

void Coupled1D_ES::solve(	PlasmaData* pData,
							FieldData* fields,
							HOMoments* moments_next_in)
{
	int nx = pData->nx;
	int ny = pData->ny;
	int nz = pData->nz;
	// Print out the solve execution
//	printf("solve() is called under the LO solver");
	// Solve the low order system
	solver -> solve(	moments_next_in, moments_old,
						fields, fields_old);
	// Update field data
	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int i=0;i<nx+1;i++)
			{
				fields->getE(i,j,k,0) = fields->getE(i,0,0,0);
				fields_half->getE(i,j,k,0) = 0.5*(fields_next->getE(i,0,0,0) + fields_old->getE(i,0,0,0));
			}
		}
	}
	// Calculate the residual for Ampere equation using LO field and HO current
	calc_residual( 	pData,
					fields, fields_old,
					moments_next);

}


realkind Coupled1D_ES::calc_residual(	PlasmaData* pData,
										FieldData* fields_next,
										FieldData* fields_old,
										HOMoments* curHOMoments)
{
	//=============================================================
	// 	Extract some simulation parameters
	//=============================================================
	int		nx 			= pData->nx;
	int 	nspecies 	= pData->nspecies;
	double 	dt 			= pData->dt;
	double 	dtRecip 	= 1.0/dt;
	double 	dx 			= pData->dxdi;
	double 	eps0 		= epsilon_naught;

	//=============================================================
	//	Calculate the average HO current
	//=============================================================
	double 	j_avg 		= 0.;
	double 	qs;
	for (int s = 0; s < nspecies; s++) {
		qs 	= pData->qspecies[s]*qe;
		for (int i = 0; i < nx; i++) {
			j_avg += qs*curHOMoments->get_val(i,0,0,s,HOMoments_currentx);
		}
	}
	j_avg 	/= double(nx);
	//=============================================================
	//	Calculate the L2 norm of Ampere equation with LO field and
	//	HO currents
	//=============================================================
	// Initialize the L2 norm
	double L2_FE 		= 0.;	// L2 of Ampere equation
	double FE_i;				// residual at cell i
	// Calculate the L2 norm
	for (int i = 0; i < nx; i++) {
		FE_i 		= eps0*dtRecip*(fields_next->getE(i,0,0,0) - fields_old->getE(i,0,0,0));
		for (int s = 0; s < nspecies; s++) {
			qs 		= 	pData->qspecies[s]*qe;
			FE_i	+=	qs*curHOMoments->get_val(i,0,0,s,HOMoments_currentx);
		}
		FE_i 		-= 	j_avg;
	}
	// Scale the L2 norm of nonlinear residual for the Ampere equation
	L2_FE 			= sqrt(dx*FE_i*FE_i);
	// Return the L2 norm of nonlinea residual for the Ampere equation
	return realkind(L2_FE);
}
