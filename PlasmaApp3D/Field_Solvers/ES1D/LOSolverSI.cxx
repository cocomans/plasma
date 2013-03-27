/*=========================================================================
                                                                                
Copyright (c) 2011, Los Alamos National Security, LLC

All rights reserved.

Copyright 2011. Los Alamos National Security, LLC. 
This software was produced under U.S. Government contract DE-AC52-06NA25396 
for Los Alamos National Laboratory (LANL), which is operated by 
Los Alamos National Security, LLC for the U.S. Department of Energy. 
The U.S. Government has rights to use, reproduce, and distribute this software. 
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  
If software is modified to produce derivative works, such modified software 
should be clearly marked, so as not to confuse it with the version available 
from LANL.
 
Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
-   Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer. 
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution. 
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software 
    without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
                                                                                
=========================================================================*/

#include "PlasmaUtility.h"
#include "LOSolverSI.h"

//=====================================================================
//	Constructor
//=====================================================================

LOSolverSI::LOSolverSI(	PlasmaData* pDataIn,
						HOMoments* curHOMoments, HOMoments* oldHOMoments,
						FieldData* curFieldData, FieldData* oldFieldData)
{
	this->pData 		= pDataIn;
	curLOMomentsES 		= new LOMomentsES(pDataIn, curHOMoments);
	oldLOMomentsES 		= new LOMomentsES(pDataIn, oldHOMoments);
	consistencyTermES	= new ConsistencyTermES(pDataIn);
	int nspecies 		= pData->nspecies;
	for (int s = 0; s < nspecies; s ++) {
		consistencyTermES->consistencyCalcContinuity(s,curHOMoments,oldHOMoments);
		consistencyTermES->consistencyCalcMomentum(s,curHOMoments,oldHOMoments,curFieldData,oldFieldData);
	}
//	printf("Right before initializing the curIntField() \n");
	curIntField = new FieldDataCPU();
//	printf("Right before initializing the curIntField() 2 \n");
	curIntField -> allocate(pDataIn);
//	printf("Right before initializing the curIntField() 3 \n");
	for(int i = 0; i < pData->nx;i++) {
		curIntField->getE(i,0,0,0) 	= curFieldData->getE(i,0,0,0);
	}
//	printf("Initialization of LOSolverSI was successful() \n");
}
//=====================================================================
//	Destructor
//=====================================================================
LOSolverSI::~LOSolverSI()
{
	delete(curLOMomentsES);
	delete(oldLOMomentsES);
	delete(consistencyTermES);
}

//=====================================================================
// Solver high level
//=====================================================================

int LOSolverSI::solve(
						HOMoments* curHOMoments,         // updates values
						HOMoments* oldHOMoments,
						FieldData* curFieldData,
						FieldData* oldFieldData)
{
	//===============================================
	//	Declare variable
	//===============================================
	int k_inner;
	//===============================================
	//	Choose which solver based on the number of species specified by user
	//===============================================
	if (pData->nspecies == 1)
		k_inner = siSolverElectron(curHOMoments, oldHOMoments, curFieldData, oldFieldData);
	else
		k_inner = siSolverTwoSpecies(curHOMoments, oldHOMoments, curFieldData, oldFieldData);
	//===============================================
	//	A return statement
	//===============================================
	return k_inner;
}

//=====================================================================
// Semi-implicit low order solver for electron only
//=====================================================================

int LOSolverSI::siSolverElectron( 	HOMoments* curHOMoments, HOMoments* oldHOMoments,
									FieldData* curFieldData, FieldData* oldFieldData)
{
	// Declare some variable
	int 		nx 		= pData->nx;
	double 		dt 		= pData->dt;
	double 		dx 		= pData->dxdi;
	double 		dxRecip = 1.0/dx;
	double 		eps0 	= epsilon_naught;
	double 		tol_lo 	= 1.0e-10;
	int 		loTrunc = 100;
	double 		me 		= pData->mspecies[0];

	// Pre-allocate the memory for variables
	double* 	s2e_pb_half 	= (double*)malloc((nx)*sizeof(double));
	double* 	gammaNe 		= (double*)malloc((nx)*sizeof(double));
	double* 	gammaNUe 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	ne 				= (double*)malloc((nx)*sizeof(double));
	double* 	ne_old 			= (double*)malloc((nx)*sizeof(double));
	double* 	ne_half 		= (double*)malloc((nx)*sizeof(double));
	double* 	nue_half 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	nue_half_old 	= (double*)malloc((nx+1)*sizeof(double));
	double* 	E_old 			= (double*)malloc((nx+1)*sizeof(double));
	double* 	E 				= (double*)malloc((nx+1)*sizeof(double));

	double* 	ne_f 			= (double*)malloc((nx+1)*sizeof(double));
	double* 	ne_f_half 		= (double*)malloc((nx+1)*sizeof(double));

	// Calculate the HO density normalized stress tensor
	for(int i = 0; i < nx; i++) {
		s2e_pb_half[i] 	= 0.5*(	curHOMoments->get_val(i,0,0,0,HOMoments_S2xx) +
								oldHOMoments->get_val(i,0,0,0,HOMoments_S2xx))/
								curHOMoments->get_val(i,0,0,0,HOMoments_charge);
	}
	// Copy pointer of LOMomentsES to HOMoments
	oldLOMomentsES->charge 		= oldHOMoments->charge;
	oldLOMomentsES->currentx	= oldHOMoments->currentx;
	// Calculate the consistency term
	consistencyTermES->consistencyCalcContinuity(0,curHOMoments,oldHOMoments);
	consistencyTermES->consistencyCalcMomentum(0,curHOMoments,oldHOMoments,curIntField,oldFieldData);
	// Copy values to allocated array
	for(int i = 0; i < nx; i++) {
		// For face values
		if (i == 0) {
			gammaNUe[i] 		= consistencyTermES->get_val(i,0,0,0,consistencyterm_momentum);
//			gammaNUe[i] 		= 0.;
			gammaNUe[nx] 		= gammaNUe[0];
			nue_half[i] 		= curLOMomentsES->get_val(i,0,0,0,LOMoments_currentx);
			nue_half_old[i] 	= oldLOMomentsES->get_val(i,0,0,0,LOMoments_currentx);
			nue_half[nx] 		= nue_half[0];
			nue_half_old[nx] 	= nue_half_old[0];
			E[i] 				= curFieldData->getE(i,0,0,0);
			E_old[i] 			= oldFieldData->getE(i,0,0,0);
			E[nx] 				= E[0];
			E_old[nx] 			= E_old[0];
		}
		else {
			gammaNUe[i] 		= consistencyTermES->get_val(i,0,0,0,consistencyterm_momentum);
//			gammaNUe[i] 		= 0.;
			nue_half[i] 		= curLOMomentsES->get_val(i,0,0,0,LOMoments_currentx);
			nue_half_old[i] 	= oldLOMomentsES->get_val(i,0,0,0,LOMoments_currentx);
			E[i]	 			= curFieldData->getE(i,0,0,0);
			E_old[i] 			= oldFieldData->getE(i,0,0,0);
		}
		// For cell centered values
		ne[i] 			= curLOMomentsES->get_val(i,0,0,0,LOMoments_charge);
		ne_old[i] 		= oldLOMomentsES->get_val(i,0,0,0,LOMoments_charge);
//		gammaNe[i] 		= consistencyTermES->get_val(i,0,0,0,consistencyterm_continuity);
		gammaNe[i] 		= 0.0;
	}
	// Calculate initial non-linear residual for low order system for tolerance
	NLResidual* residual = new NLResidual(this->pData);
	double maxResidual0, maxResidual;
	maxResidual0 = residual->calculateNLResidual(	curHOMoments, oldHOMoments,
													curLOMomentsES, oldLOMomentsES,
													curFieldData, oldFieldData,
													consistencyTermES);
//  cout << endl <<  "   The initial inner: " << maxResidual0 << endl;
//	residual->printResidual();

	// Declare average current
	double 	j_avg;
	// Convergence loop
	int 	flag_conv = 0;
	int 	k_inner = 0;
	while (flag_conv == 0) {
		// Increment the iteration counter
		k_inner++;
//		cout << "This is the " << k_inner <<"th inner iteration." << endl;

		// Calculate electron density face quantity for n+1/2, p+1
		for (int i = 0; i < nx; i++) {
			ne_half[i] 	= 0.5*(ne[i] + ne_old[i]);
		}
		makeFaceFlux(ne_half,ne_f_half,nx);
		makeFaceFlux(ne,ne_f,nx);
		// Calculate the new average current at the current Picard iteration
		j_avg 	= 0.0;
		for (int i = 0; i < nx; i++) {
			j_avg	+= nue_half[i]*pData->qspecies[0];
//			j_avg 	+= curHOMoments->get_val(i,0,0,0,HOMoments_currentx)*pData->qspecies[0];
		}
		j_avg	/= 	double(nx);
		// Back calculate to compute number density
		for (int i = 0; i < nx; i++) {
			ne[i] 	= (ne_old[i] - dt*dxRecip*(nue_half[i+1] - nue_half[i]))/(1 - dt*gammaNe[i]/me);
		}
		// Back calculate to compute electric field
		for (int i = 0; i < nx; i++) {
			E[i] 	= E_old[i] - dt*qe*(nue_half[i]*pData->qspecies[0] - j_avg)/eps0;
//			E[i] 	= E_old[i] - dt*qe*(curHOMoments->get_val(i,0,0,0,HOMoments_currentx) - j_avg);
		}
		E[nx] 	= E[0];
		// Linear solver returning rue
		solveLinearSystemElectron(	dt,
									dx,
									s2e_pb_half,
									pData->qspecies[0]*qe,
									me,
									eps0,
									gammaNe,
									gammaNUe,
									ne,
									ne_old,
									ne_f,
									ne_f_half,
									E_old,
									nue_half_old,
									j_avg,
									nue_half);

		// Store the updated values to the LOMomentsES class (damn why did I decide to make a LOMomentsES class... :( )
		for (int i = 0; i < nx; i++) {
			curLOMomentsES->get_val(i,0,0,0,LOMoments_charge)	= ne[i];
			curLOMomentsES->get_val(i,0,0,0,LOMoments_currentx) = nue_half[i];
			curFieldData->getE(i,0,0,0) 						= E[i];
		}
		// Check for convergence by calculating the L2 of non-linear residual
		maxResidual = residual->calculateNLResidual(	curHOMoments, oldHOMoments,
    													curLOMomentsES, oldLOMomentsES,
    													curFieldData, oldFieldData,
    													consistencyTermES);
//		residual->printResidual();

		// Check for convergence in the Picard loop (particle solution)
		if (maxResidual <= (1.0e-6 * maxResidual0) || k_inner == loTrunc) {
			flag_conv = 1;
		}
//    getchar();
	}
	for (int i = 0;i < nx; i++) {
		for (int j = 0; j < pData->ny; j++) {
			for (int k = 0; k < pData->nz; k++) {
				curFieldData->getE(i,j,k,0) 	= E[i];
				curIntField->getE(i,j,k,0) 		= E[i];
			}
		}
	}
//	// Calculate the consistency term
//	consistencyTermES->consistencyCalcContinuity(0,curHOMoments,oldHOMoments);
//	consistencyTermES->consistencyCalcMomentum(0,curHOMoments,oldHOMoments,curFieldData,oldFieldData);

//	// check the L2 norm of relative difference between ne_LO, ne_HO and nue_LO, nue_HO
//	double* relDiffne 	= (double*)malloc(nx*sizeof(double));
//	double* relDiffnue 	= (double*)malloc((nx+1)*sizeof(double));
//	double L2_ne = 0.;
//	double L2_nue = 0.;
//	for (int i = 0;i < nx; i++) {
//		relDiffne[i] 		= 	fabs((curHOMoments->get_val(i,0,0,0,HOMoments_charge) - ne[i] )/
//								(curHOMoments->get_val(i,0,0,0,HOMoments_charge)));
//		L2_ne 				+= 	sqrt(relDiffne[i]*relDiffne[i]*dx);
//		printf("relDiffne[%i] = %e \n",i,relDiffne[i]);
//	}
//	for (int i = 0;i < nx+1; i++){
//		relDiffnue[i] 		= 	fabs((curHOMoments->get_val(i,0,0,0,HOMoments_currentx) - nue_half[i])/
//								(curHOMoments->get_val(i,0,0,0,HOMoments_currentx)));
//		L2_nue 				+= 	sqrt(relDiffnue[i]*relDiffnue[i]*dx);
//		printf("relDiffnue[%i] = %e \n",i,relDiffnue[i]);
//	}
//	printf("L2_reldiff_ne = %e \n",L2_ne);
//	printf("L2_reldiff_nue = %e \n",L2_nue);

	return k_inner;
	// Free-up allocated memorys
	free(s2e_pb_half);
	free(gammaNe);
	free(gammaNUe);
	free(ne);
	free(ne_old);
	free(ne_half);
	free(nue_half);
	free(nue_half_old);
	free(E_old);
	free(E);
	free(ne_f);
	free(ne_f_half);
}

//=====================================================================
// Semi-implicit low order solver for two-species case
//=====================================================================

int LOSolverSI::siSolverTwoSpecies( 	HOMoments* curHOMoments, HOMoments* oldHOMoments,
										FieldData* curFieldData, FieldData* oldFieldData)
{
	// Declare some variable
	int 		nx 		= pData->nx;
	double 		dt 		= pData->dt;
	double 		dx 		= pData->dxdi;
	double 		dxRecip = 1.0/dx;
	double 		eps0 	= epsilon_naught;
	double 		tol_lo 	= 1.0e-12;
	int 		loTrunc = 500;
	double 		me 		= pData->mspecies[0];
	double 		mi 		= pData->mspecies[1];
	int 		nspecies= pData->nspecies;
	double 		qes 	= qe*pData->qspecies[0];
	double 		qis 	= qe*pData->qspecies[1];
	double 		last_residual = 0;

	// Pre-allocate the memory for variables
	double* 	s2e_pb_half 	= (double*)malloc((nx)*sizeof(double));
	double* 	s2i_pb_half 	= (double*)malloc((nx)*sizeof(double));
	double* 	gammaNe 		= (double*)malloc((nx)*sizeof(double));
	double* 	gammaNi 		= (double*)malloc((nx)*sizeof(double));
	double* 	gammaNUe 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	gammaNUi 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	ne 				= (double*)malloc((nx)*sizeof(double));
	double* 	ni 				= (double*)malloc((nx)*sizeof(double));
	double* 	ne_old 			= (double*)malloc((nx)*sizeof(double));
	double* 	ni_old 			= (double*)malloc((nx)*sizeof(double));
	double* 	ne_half 		= (double*)malloc((nx)*sizeof(double));
	double* 	ni_half 		= (double*)malloc((nx)*sizeof(double));
	double* 	nue_half 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	nui_half 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	nue_half_old 	= (double*)malloc((nx+1)*sizeof(double));
	double* 	nui_half_old 	= (double*)malloc((nx+1)*sizeof(double));
	double* 	E_old 			= (double*)malloc((nx+1)*sizeof(double));
	double* 	E 				= (double*)malloc((nx+1)*sizeof(double));

	double* 	ne_f 			= (double*)malloc((nx+1)*sizeof(double));
	double* 	ni_f 			= (double*)malloc((nx+1)*sizeof(double));
	double* 	ne_f_half 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	ni_f_half 		= (double*)malloc((nx+1)*sizeof(double));

	double* 	E_source 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	E_factor 		= (double*)malloc((nx+1)*sizeof(double));
	double* 	nui_source 		= (double*)malloc((nx+1)*sizeof(double));

	// Calculate the HO density normalized stress tensor
	for(int i = 0; i < nx; i++) {
		s2e_pb_half[i] 	= 0.5*(	curHOMoments->get_val(i,0,0,0,HOMoments_S2xx) +
								oldHOMoments->get_val(i,0,0,0,HOMoments_S2xx))/
								curHOMoments->get_val(i,0,0,0,HOMoments_charge);

		s2i_pb_half[i] 	= 0.5*(	curHOMoments->get_val(i,0,0,1,HOMoments_S2xx) +
								oldHOMoments->get_val(i,0,0,1,HOMoments_S2xx))/
								curHOMoments->get_val(i,0,0,1,HOMoments_charge);
	}
	// Copy pointer of LOMomentsES to HOMoments
	oldLOMomentsES->charge 		= oldHOMoments->charge;
	oldLOMomentsES->currentx	= oldHOMoments->currentx;
	// Calculate the consistency term
	for (int s = 0;s < nspecies; s++) {
		consistencyTermES->consistencyCalcContinuity(s,curHOMoments,oldHOMoments);
		consistencyTermES->consistencyCalcMomentum(s,curHOMoments,oldHOMoments,curIntField,oldFieldData);
	}
	// Copy values to allocated array
	for(int i = 0; i < nx; i++) {
		// For face values
		if (i == 0) {
			gammaNUe[i] 		= consistencyTermES->get_val(i,0,0,0,consistencyterm_momentum);
			gammaNUi[i] 		= consistencyTermES->get_val(i,0,0,1,consistencyterm_momentum);
			gammaNUe[nx] 		= gammaNUe[0];
			gammaNUi[nx] 		= gammaNUi[0];
			nue_half[i] 		= curLOMomentsES->get_val(i,0,0,0,LOMoments_currentx);
			nue_half_old[i] 	= oldLOMomentsES->get_val(i,0,0,0,LOMoments_currentx);
			nue_half[nx] 		= nue_half[0];
			nue_half_old[nx] 	= nue_half_old[0];
			nui_half[i] 		= curLOMomentsES->get_val(i,0,0,1,LOMoments_currentx);
			nui_half_old[i] 	= oldLOMomentsES->get_val(i,0,0,1,LOMoments_currentx);
			nui_half[nx] 		= nui_half[0];
			nui_half_old[nx] 	= nui_half_old[0];
			E[i] 				= curFieldData->getE(i,0,0,0);
			E_old[i] 			= oldFieldData->getE(i,0,0,0);
			E[nx] 				= E[0];
			E_old[nx] 			= E_old[0];
		}
		else {
			gammaNUe[i] 		= consistencyTermES->get_val(i,0,0,0,consistencyterm_momentum);
			gammaNUi[i] 		= consistencyTermES->get_val(i,0,0,1,consistencyterm_momentum);
			nue_half[i] 		= curLOMomentsES->get_val(i,0,0,0,LOMoments_currentx);
			nue_half_old[i] 	= oldLOMomentsES->get_val(i,0,0,0,LOMoments_currentx);
			nui_half[i] 		= curLOMomentsES->get_val(i,0,0,1,LOMoments_currentx);
			nui_half_old[i] 	= oldLOMomentsES->get_val(i,0,0,1,LOMoments_currentx);
			E[i]	 			= curFieldData->getE(i,0,0,0);
			E_old[i] 			= oldFieldData->getE(i,0,0,0);
		}
		// For cell centered values
		ne[i] 			= curLOMomentsES->get_val(i,0,0,0,LOMoments_charge);
		ne_old[i] 		= oldLOMomentsES->get_val(i,0,0,0,LOMoments_charge);
		ni[i] 			= curLOMomentsES->get_val(i,0,0,1,LOMoments_charge);
		ni_old[i] 		= oldLOMomentsES->get_val(i,0,0,1,LOMoments_charge);
		gammaNe[i] 		= consistencyTermES->get_val(i,0,0,0,consistencyterm_continuity);
		gammaNi[i] 		= consistencyTermES->get_val(i,0,0,1,consistencyterm_continuity);
//		gammaNe[i] 		= 0.;
//		gammaNi[i] 		= 0.;
	}
	// Calculate initial non-linear residual for low order system for tolerance
	NLResidual* residual = new NLResidual(this->pData);
	double maxResidual0, maxResidual;
	maxResidual0 = residual->calculateNLResidual(	curHOMoments, oldHOMoments,
													curLOMomentsES, oldLOMomentsES,
												   	curIntField, oldFieldData,
													consistencyTermES);
//  cout << endl <<  "   The initial inner: " << maxResidual0 << endl;
//	residual->printResidual();

	// Declare average current
	double 	j_avg;
	// Convergence loop
	int 	flag_conv = 0;
	int 	k_inner = 0;
	while (flag_conv == 0) {
		// Increment the iteration counter
		k_inner++;
//		cout << "This is the " << k_inner <<"th inner iteration." << endl;

		// Calculate the new average current at the current Picard iteration
		j_avg 	= 0.0;
		for (int i = 0; i < nx; i++) {
			j_avg	= j_avg + nue_half[i]*qes + nui_half[i]*qis;
		}
		j_avg	= 	j_avg/double(nx);
		// Calculate the new ion number density
		for (int i = 0; i < nx ; i++) {
			ni[i] 	= (ni_old[i] - (dt/dx)*(nui_half[i+1] - nui_half[i]))/(1 - dt*gammaNi[i]/mi);
		}
		// Calculate density face quantity for n+1/2, p+1
		for (int i = 0; i < nx; i++) {
			ne_half[i] 	= 0.5*(ne[i] + ne_old[i]);
			ni_half[i] 	= 0.5*(ni[i] + ni_old[i]);
		}
		makeFaceFlux(ne_half,ne_f_half,nx);
		makeFaceFlux(ne,ne_f,nx);
		makeFaceFlux(ni_half,ni_f_half,nx);
		makeFaceFlux(ni,ni_f,nx);

		// Calculate the nui_source
		for (int i = 0; i < nx+1; i++) {
			if (i == 0 || i == nx) {
				nui_source[i]	= 	nui_half_old[i] - (dt/dx)*(ni[0]*s2i_pb_half[0] - ni[nx-1]*s2i_pb_half[nx-1]) +
									(qis/mi)*dt*ni_f_half[i]*E_old[i]*0.5 + (dt/mi)*gammaNUi[i]*ni_f[i];
			}
			else {
				nui_source[i]	= 	nui_half_old[i] - (dt/dx)*(ni[i]*s2i_pb_half[i] - ni[i-1]*s2i_pb_half[i-1]) +
									(qis/mi)*dt*ni_f_half[i]*E_old[i]*0.5 + (dt/mi)*gammaNUi[i]*ni_f[i];
			}
		}
		// Calculate the source term from Ampere equation
		for (int i = 0; i < nx+1; i++) {
			E_factor[i] 	= 1.0/(1.0 + dt*dt*qis*qis*ni_f_half[i]*0.5/mi/eps0);
			E_source[i] 	= E_old[i] - (dt/eps0)*(qis*nui_source[i] - j_avg);
		}
		// Linear solver returning rue
		solveLinearSystemTwoSpecies(
									dt,
									dx,
									s2e_pb_half,
									qes,
									me,
									eps0,
									gammaNe,
									gammaNUe,
									ne,
									ne_old,
									ne_f,
									ne_f_half,
									E_old,
									E_factor,
									E_source,
									nue_half_old,
									j_avg,
									nue_half);           // output

		// Back calculate to compute number density
		for (int i = 0; i < nx; i++) {
			ne[i] 			= (ne_old[i] - dt*dxRecip*(nue_half[i+1] - nue_half[i]))/(1 - dt*gammaNe[i]/me);
		}
		// Back calculate to compute electric field
		for (int i = 0; i < nx+1; i++) {
			E[i] 			= E_factor[i]*E_source[i] - dt*qes*nue_half[i]*E_factor[i]/eps0;
		}
		// Back calculate to compute ion momentum
		for (int i = 0; i < nx+1; i++) {
			nui_half[i] 	= nui_source[i] + dt*qis*ni_f_half[i]*E[i]*0.5/mi;
		}

		// Store the updated values to the LOMomentsES class (damn why did I decide to make a LOMomentsES class... :( )
		for (int i = 0; i < nx; i++) {
			curLOMomentsES->get_val(i,0,0,0,LOMoments_charge)	= ne[i];
			curLOMomentsES->get_val(i,0,0,0,LOMoments_currentx) = nue_half[i];
			curLOMomentsES->get_val(i,0,0,1,LOMoments_charge) 	= ni[i];
			curLOMomentsES->get_val(i,0,0,1,LOMoments_currentx) = nui_half[i];
			curFieldData->getE(i,0,0,0) 						= E[i];
		}
		last_residual = maxResidual;
		// Check for convergence by calculating the L2 of non-linear residual
		maxResidual = residual->calculateNLResidual(	curHOMoments, oldHOMoments,
    													curLOMomentsES, oldLOMomentsES,
    													curFieldData, oldFieldData,
    													consistencyTermES);
//		residual->printResidual();

		double residual_change = fabs(maxResidual-last_residual)/(maxResidual+last_residual);

		// Check for convergence in the Picard loop (particle solution)
		if (maxResidual <= (tol_lo) || k_inner == loTrunc || residual_change <= 100*tol_lo) {
			flag_conv = 1;
		}
//    getchar();
	}
	for (int i = 0;i < nx; i++) {
		for (int j = 0; j < pData->ny; j++) {
			for (int k = 0; k < pData->nz; k++) {
				curFieldData->getE(i,j,k,0) 	= E[i];
				curIntField->getE(i,j,k,0) 		= E[i];
			}
		}
	}
	// Calculate the consistency term
	for (int s = 0;s < nspecies;s++) {
		consistencyTermES->consistencyCalcContinuity(s,curHOMoments,oldHOMoments);
		consistencyTermES->consistencyCalcMomentum(s,curHOMoments,oldHOMoments,curFieldData,oldFieldData);
	}
//	for (int i = 0;i < nx+1; i++) {
//		printf("gammaNUe[%i] = %e \n",i,gammaNUe[i]);
//		printf("gammaNUi[%i] = %e \n",i,gammaNUi[i]);
//	}
//	getchar();
//		// check the L2 norm of relative difference between ne_LO, ne_HO and nue_LO, nue_HO
//		double* relDiffne 	= (double*)malloc(nx*sizeof(double));
//		double* relDiffnue 	= (double*)malloc((nx+1)*sizeof(double));
//		double L2_ne = 0.;
//		double L2_nue = 0.;
//		for (int i = 0;i < nx; i++) {
//			relDiffne[i] 		= 	fabs((curHOMoments->get_val(i,0,0,0,HOMoments_charge) - ne[i] )/
//									(curHOMoments->get_val(i,0,0,0,HOMoments_charge)));
//			L2_ne 				+= 	sqrt(relDiffne[i]*relDiffne[i]*dx);
//			printf("relDiffne[%i] = %e \n",i,relDiffne[i]);
//		}
//		for (int i = 0;i < nx+1; i++){
//			relDiffnue[i] 		= 	fabs((curHOMoments->get_val(i,0,0,0,HOMoments_currentx) - nue_half[i])/
//									(curHOMoments->get_val(i,0,0,0,HOMoments_currentx)));
//			L2_nue 				+= 	sqrt(relDiffnue[i]*relDiffnue[i]*dx);
//			printf("relDiffnue[%i] = %e \n",i,relDiffnue[i]);
//		}
//		printf("L2_reldiff_ne = %e \n",L2_ne);
//		printf("L2_reldiff_nue = %e \n",L2_nue);
//	getchar();
	return k_inner;
	// Free-up allocated memorys
	free(s2e_pb_half);
	free(gammaNe);
	free(gammaNUe);
	free(ne);
	free(ne_old);
	free(ne_half);
	free(ne_f);
	free(ne_f_half);
	free(nue_half);
	free(nue_half_old);

	free(s2i_pb_half);
	free(gammaNi);
	free(gammaNUi);
	free(ni);
	free(ni_old);
	free(ni_half);
	free(ni_f);
	free(ni_f_half);
	free(nui_half);
	free(nui_half_old);

	free(E_old);
	free(E);
}

//=====================================================================
// Triadiagonal linear solver for single species electron case
//=====================================================================

void LOSolverSI::solveLinearSystemElectron(
											double 		dt,
											double 		dx,
											double* 	s2e_pb_half,
											double 		qes,
											double 		me,
											double 		eps0,
											double* 	gammaN,
											double* 	gammaNU,
											double* 	ne,
											double* 	ne_old,
											double* 	ne_f,
											double* 	ne_f_half,
											double* 	E_old,
											double* 	nue_half_old,
											double 		j_avg,
											double* 	nue_half)           // output
{
	// Declare some constants for optimization
	int		nx 				= pData->nx;
	double 	mdtRecip		= me/dt;
	double 	dxRecip 		= 1.0/dx;
	double 	mdxRecip 		= me*dxRecip;
	double 	mdxsqRecip 		= me*dxRecip*dxRecip;
	double 	dteps0Recip 	= dt/eps0;
	double 	mdtdxsqRecip 	= me*dt*dxRecip*dxRecip;
	double 	qesq 			= qes*qes;
	// Declare matrix and rhs
	double** 	A 			= new double*[this->pData->nx+1];
	double* 	rhs 		= new double[this->pData->nx+1];
	// Allocate memory for E_source
	double* E_source 		= (double*)malloc((nx+1)*sizeof(double));
	// Initialize the matrix components
	for (int i = 0; i < this->pData->nx+1; i++) {
		A[i] = new double[this->pData->nx+1];
		for (int j = 0; j < this->pData->nx+1; j++) {
			A[i][j] = 0.0;
		}
	}
	// Fill E_source value
	for (int i = 0; i < nx+1; i++) {
		E_source[i] 	= E_old[i] + dt*qes*j_avg/eps0;
	}
	// Form the coefficient matrix for electron momentum equation
	for (int fIndx = 0; fIndx < this->pData->nx+1; fIndx++) {
		if (fIndx == 0 || fIndx == this->pData->nx) {
			A[fIndx][fIndx]		= 	mdtRecip +
									mdtdxsqRecip*(s2e_pb_half[0] + s2e_pb_half[nx-1]) +
									qes*qes*0.5*dt*ne_f_half[fIndx]/eps0;

			A[fIndx][1] 		=	-mdtdxsqRecip*s2e_pb_half[0];

			A[fIndx][nx-1] 		=	-mdtdxsqRecip*s2e_pb_half[nx-1];
		}
		else {
			A[fIndx][fIndx] 	= 	mdtRecip +
									mdtdxsqRecip*(s2e_pb_half[fIndx] + s2e_pb_half[fIndx-1]) +
									qes*qes*0.5*dt*ne_f_half[fIndx]/eps0;

			A[fIndx][fIndx + 1] =	-mdtdxsqRecip*s2e_pb_half[fIndx];

			A[fIndx][fIndx - 1] =	-mdtdxsqRecip*s2e_pb_half[fIndx-1];
		}
	}
	// Right hand side
//	cout << endl;
	for (int fIndx = 0; fIndx < this->pData->nx+1; fIndx++) {
		if (fIndx == 0 || fIndx == this->pData->nx) {
			rhs[fIndx] 			=	mdtRecip*nue_half_old[fIndx] -
									mdxRecip*(	ne_old[0]*s2e_pb_half[0] -
												ne_old[nx-1]*s2e_pb_half[nx-1]
											 ) +
									0.5*qes*ne_f_half[fIndx]*E_source[fIndx] +
									0.5*qes*ne_f_half[fIndx]*E_old[fIndx] +
									gammaNU[fIndx]*ne_f[fIndx];
		}
		else {
			rhs[fIndx] 			=	mdtRecip*nue_half_old[fIndx] -
									mdxRecip*(	ne_old[fIndx]*s2e_pb_half[fIndx] -
												ne_old[fIndx-1]*s2e_pb_half[fIndx-1]
											) +
									0.5*qes*ne_f_half[fIndx]*E_source[fIndx] +
									0.5*qes*ne_f_half[fIndx]*E_old[fIndx] +
									gammaNU[fIndx]*ne_f[fIndx];


		}
	}
	// Calculate the new time electron momentum from direct inversion
	double residual;
	int flag;
	GaussElim(A, this->pData->nx+1, rhs, this->pData->nx+1, nue_half, residual, flag);
	// Free allocated memory
	for (int fIndx = 0; fIndx < this->pData->nx+1; fIndx++) {
	    delete [] A[fIndx];
	}
	delete [] A;
	delete [] rhs;
	free(E_source);
}
void LOSolverSI::solveLinearSystemTwoSpecies(
											double 		dt,
											double 		dx,
											double* 	s2e_pb_half,
											double 		qes,
											double 		me,
											double 		eps0,
											double* 	gammaN,
											double* 	gammaNU,
											double* 	ne,
											double* 	ne_old,
											double* 	ne_f,
											double* 	ne_f_half,
											double* 	E_old,
											double* 	E_factor,
											double* 	E_source,
											double* 	nue_half_old,
											double 		j_avg,
											double* 	nue_half)		// output
{
	// Declare some constants for optimization
	int		nx 				= pData->nx;
	double 	mdtRecip		= me/dt;
	double 	dxRecip 		= 1.0/dx;
	double 	mdxRecip 		= me*dxRecip;
	double 	mdxsqRecip 		= me*dxRecip*dxRecip;
	double 	dteps0Recip 	= dt/eps0;
	double 	mdtdxsqRecip 	= me*dt*dxRecip*dxRecip;
	double 	qesq 			= qes*qes;
	// Declare matrix and rhs
	double** 	A 			= new double*[this->pData->nx+1];
	double* 	rhs 		= new double[this->pData->nx+1];
	// Initialize the matrix components
	for (int i = 0; i < this->pData->nx+1; i++) {
		A[i] = new double[this->pData->nx+1];
		for (int j = 0; j < this->pData->nx+1; j++) {
			A[i][j] = 0.0;
		}
	}
	// Form the coefficient matrix for electron momentum equation
	for (int fIndx = 0; fIndx < this->pData->nx+1; fIndx++) {
		if (fIndx == 0 || fIndx == this->pData->nx) {
			A[fIndx][fIndx]		= 	mdtRecip +
									mdtdxsqRecip*(s2e_pb_half[0] + s2e_pb_half[nx-1]) +
									qes*qes*0.5*dt*ne_f_half[fIndx]*E_factor[fIndx]/eps0;

			A[fIndx][1] 		=	-mdtdxsqRecip*s2e_pb_half[0];

			A[fIndx][nx-1] 		=	-mdtdxsqRecip*s2e_pb_half[nx-1];
		}
		else {
			A[fIndx][fIndx] 	= 	mdtRecip +
									mdtdxsqRecip*(s2e_pb_half[fIndx] + s2e_pb_half[fIndx-1]) +
									qes*qes*0.5*dt*ne_f_half[fIndx]*E_factor[fIndx]/eps0;

			A[fIndx][fIndx + 1] =	-mdtdxsqRecip*s2e_pb_half[fIndx];

			A[fIndx][fIndx - 1] =	-mdtdxsqRecip*s2e_pb_half[fIndx-1];
		}
	}
	// Right hand side
//	cout << endl;
	for (int fIndx = 0; fIndx < this->pData->nx+1; fIndx++) {
		if (fIndx == 0 || fIndx == this->pData->nx) {
			rhs[fIndx] 			=	mdtRecip*nue_half_old[fIndx] -
									mdxRecip*(	ne_old[0]*s2e_pb_half[0] -
												ne_old[nx-1]*s2e_pb_half[nx-1]
											 ) +
									0.5*qes*ne_f_half[fIndx]*E_source[fIndx]*E_factor[fIndx] +
									0.5*qes*ne_f_half[fIndx]*E_old[fIndx] +
									gammaNU[fIndx]*ne_f[fIndx] -
									(dt/me)*(gammaN[0]*s2e_pb_half[0]*ne[0] -
											 gammaN[nx-1]*s2e_pb_half[nx-1]*ne[nx-1]);
		}
		else {
			rhs[fIndx] 			=	mdtRecip*nue_half_old[fIndx] -
									mdxRecip*(	ne_old[fIndx]*s2e_pb_half[fIndx] -
												ne_old[fIndx-1]*s2e_pb_half[fIndx-1]
											) +
									0.5*qes*ne_f_half[fIndx]*E_source[fIndx]*E_factor[fIndx] +
									0.5*qes*ne_f_half[fIndx]*E_old[fIndx] +
									gammaNU[fIndx]*ne_f[fIndx]  -
									(dt/me)*(gammaN[fIndx]*s2e_pb_half[fIndx]*ne[fIndx] -
											gammaN[fIndx-1]*s2e_pb_half[fIndx-1]*ne[fIndx-1]);


		}
	}
	// Calculate the new time electron momentum from direct inversion
	double residual;
	int flag;
	GaussElim(A, this->pData->nx+1, rhs, this->pData->nx+1, nue_half, residual, flag);
	// Free allocated memory
	for (int fIndx = 0; fIndx < this->pData->nx+1; fIndx++) {
	    delete [] A[fIndx];
	}
	delete [] A;
	delete [] rhs;
}
