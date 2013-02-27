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

#include "NLResidual.h"
#include "PlasmaUtility.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

NLResidual::NLResidual(PlasmaData* pDataIn)
{
  pData = pDataIn;

  this->F_r 	= new double*[pData->nspecies];
  this->F_ru 	= new double*[pData->nspecies];

  for (int sIndx = 0; sIndx < this->pData->nspecies; sIndx++) {
    // Cell based (differs from old version)
    this->F_r[sIndx] = new double[this->pData->nx];

    // Face based (differs from old version)
    this->F_ru[sIndx] = new double[this->pData->nx+1];
  }
  this->F_E = new double[this->pData->nx+1];

  this->maxResidual_r = new double[this->pData->nspecies];
  this->maxResidual_ru = new double[this->pData->nspecies];
}

NLResidual::~NLResidual()
{
  for (int sIndx = 0; sIndx < this->pData->nspecies; sIndx++) {
    delete [] this->F_r[sIndx];
    delete [] this->F_ru[sIndx];
  }
  delete this->F_r;
  delete this->F_ru;
  delete this->F_E;
  delete this->maxResidual_r;
  delete this->maxResidual_ru;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate the initial non-linear residual
// nlres_init.m
//
///////////////////////////////////////////////////////////////////////

double NLResidual::calculateNLResidual(
										HOMoments* 		curHOMoments,
										HOMoments* 		oldHOMoments,
										LOMomentsES* 	curLOMomentsES,
										LOMomentsES* 	oldLOMomentsES,
										FieldData* 		curFieldData,
										FieldData* 		oldFieldData,
										ConsistencyTermES* consistencyTermES)
{
  nlres_E(curHOMoments, oldHOMoments, curLOMomentsES, oldLOMomentsES, curFieldData, oldFieldData);

  this->maxResidual_E = maxL2(this->F_E, this->pData->nx+1);
  this->totMaxResidual = this->maxResidual_E;

  for (int sIndx = 0; sIndx < this->pData->nspecies; sIndx++) {

    nlres_r(sIndx,curHOMoments, oldHOMoments, curLOMomentsES, oldLOMomentsES, curFieldData, oldFieldData, consistencyTermES);
    nlres_ru(sIndx,curHOMoments, oldHOMoments, curLOMomentsES, oldLOMomentsES, curFieldData, oldFieldData, consistencyTermES);

    this->maxResidual_r[sIndx] = maxL2(this->F_r[sIndx], this->pData->nx);
    if (this->totMaxResidual < this->maxResidual_r[sIndx])
      this->totMaxResidual = this->maxResidual_r[sIndx];
    
    this->maxResidual_ru[sIndx] = maxL2(this->F_ru[sIndx], this->pData->nx+1);
    if (this->totMaxResidual < this->maxResidual_ru[sIndx])
      this->totMaxResidual = this->maxResidual_ru[sIndx];
  }
  return this->totMaxResidual;
}

///////////////////////////////////////////////////////////////////////
//
// Construct the non-linear residual function for Ampere Law equation
// Values are returned in parameter, function returns the maximum value
//
///////////////////////////////////////////////////////////////////////

void NLResidual::nlres_E(	HOMoments* 		curHOMoments,
							HOMoments* 		oldHOMoments,
							LOMomentsES* 	curLOMoments,
							LOMomentsES* 	oldLOMoments,
							FieldData* 		curFieldData,
							FieldData* 		oldFieldData)
{
	realkind* 	E 		= (realkind*)malloc((pData->nx+1)*sizeof(realkind));
	realkind* 	E_old 	= (realkind*)malloc((pData->nx+1)*sizeof(realkind));
	double 		dt 		= pData->dt;
	double 		dtRecip	= 1.0/dt;
	for	(int i = 0;i < pData->nx+1; i++) {
		E[i] 		= curFieldData->getE(i,0,0,0);
		E_old[i] 	= oldFieldData->getE(i,0,0,0);
	}

	// Sum all the currents from different species
	realkind* 	j_tot 	= (realkind*)malloc((pData->nx+1)*sizeof(realkind));
	for (int fIndx = 0; fIndx < this->pData->nx+1; fIndx++)
	{
		j_tot[fIndx] = 0.0;
		for (int sIndx = 0; sIndx < pData->nspecies; sIndx++) {
			j_tot[fIndx] += pData->qspecies[sIndx]*curLOMoments->get_val(fIndx,0,0,sIndx,LOMoments_currentx);
		}
	}
	realkind j_avg = 0.0;
	for (int cIndx = 0; cIndx < this->pData->nx; cIndx++) {
		j_avg += j_tot[cIndx];
	}
	j_avg /= realkind(pData->nx);

	// Construct the non-linear residual for the Ampere Law equation on node
	for (int fIndx = 0; fIndx < pData->nx+1; fIndx++)
	this->F_E[fIndx] = epsilon_naught * (E[fIndx] - E_old[fIndx])*dtRecip + qe*(j_tot[fIndx] - j_avg);

	free(j_tot);
	free(E);
	free(E_old);
}

///////////////////////////////////////////////////////////////////////
//
// Construct the non-linear residual function for continuity equation
//
///////////////////////////////////////////////////////////////////////

void NLResidual::nlres_r(	int sIndx,
							HOMoments* curHOMoments,
							HOMoments* oldHOMoments,
							LOMomentsES* curLOMomentsES,
							LOMomentsES* oldLOMomentsES,
							FieldData* curFieldData,
							FieldData* oldFieldData,
							ConsistencyTermES* consistencyTermES)
{
  realkind 	m 			= 	this->pData->mspecies[sIndx];
  realkind 	mdt_recip 	= 	m / this->pData->dt;
  realkind 	mdx_recip 	= 	m / this->pData->dxdi;
  int nx 				=	pData->nx;
  // Allocate memory
  realkind* r 		= (realkind*)malloc(nx*sizeof(realkind));
  realkind* r_old	= (realkind*)malloc(nx*sizeof(realkind));
  realkind* ru 		= (realkind*)malloc((nx+1)*sizeof(realkind));
  realkind* gammaN 	= (realkind*)malloc(nx*sizeof(realkind));
  for (int i = 0; i < nx; i++) {
	  r[i]  	= curLOMomentsES->get_val(i,0,0,sIndx,LOMoments_charge);
	  r_old[i] 	= oldLOMomentsES->get_val(i,0,0,sIndx,LOMoments_charge);
	  gammaN[i] = consistencyTermES->get_val(i,0,0,sIndx,consistencyterm_continuity);
//	  gammaN[i] = 0.;
	  if (i == 0) {
		  ru[0] 	= curLOMomentsES->get_val(i,0,0,sIndx,LOMoments_currentx);
		  ru[nx] 	= ru[0];
	  }
	  else {
		  ru[i] 	= curLOMomentsES->get_val(i,0,0,sIndx,LOMoments_currentx);
	  }
  }
  // Construct the non-linear residual for the continuity equation on node
  for (int cIndx = 0; cIndx < this->pData->nx; cIndx++) {
	  this->F_r[sIndx][cIndx] =
			   (mdt_recip * (r[cIndx] - r_old[cIndx]) +
				mdx_recip * (ru[cIndx + 1] - ru[cIndx]) -
				gammaN[cIndx]*r[cIndx]);
//	  printf("gammaN[%i] = %e \n",cIndx,gammaN[cIndx]);
  }
  // Free the allocated memory
  free(r);
  free(r_old);
  free(ru);
  free(gammaN);
}

///////////////////////////////////////////////////////////////////////
//
// Construct the non-linear residual function for momentum equation
//
///////////////////////////////////////////////////////////////////////

void NLResidual::nlres_ru(
							int sIndx,
							HOMoments* curHOMoments,
							HOMoments* oldHOMoments,
							LOMomentsES* curLOMomentsES,
							LOMomentsES* oldLOMomentsES,
							FieldData* curFieldData,
							FieldData* oldFieldData,
							ConsistencyTermES* consistencyTermES)
{
  realkind 	q 			= this->pData->qspecies[sIndx]*qe;
  realkind 	m 			= this->pData->mspecies[sIndx];
  realkind 	mdt_recip 	= m / this->pData->dt;
  realkind 	mdx_recip 	= m / this->pData->dxdi;
  int 		nx 			= pData->nx;
  // Extract moment quantities at time n+1 and n+1/2
  realkind* r 	= (realkind*)malloc(nx*sizeof(realkind));
  realkind* ru 	= (realkind*)malloc((nx+1)*sizeof(realkind));
  realkind* E 	= (realkind*)malloc((nx+1)*sizeof(realkind));

  // Extract moment quantities at time n
  realkind* r_old	= (realkind*)malloc(nx*sizeof(realkind));
  realkind* ru_old 	= (realkind*)malloc((nx+1)*sizeof(realkind));
  realkind* E_old 	= (realkind*)malloc((nx+1)*sizeof(realkind));
  // Extract HO moment quantities
  realkind*	r_p 		= (realkind*)malloc(nx*sizeof(realkind));
  realkind* r_p_old 	= (realkind*)malloc(nx*sizeof(realkind));
  realkind* s2_p 		= (realkind*)malloc(nx*sizeof(realkind));
  realkind* s2_p_half	= (realkind*)malloc(nx*sizeof(realkind));
  realkind* s2_pb_half 	= (realkind*)malloc(nx*sizeof(realkind));
  for (int i = 0; i < nx; i++) {
	  r[i] 			=	curLOMomentsES->get_val(i,0,0,sIndx,LOMoments_charge);
	  r_old[i]		=	oldLOMomentsES->get_val(i,0,0,sIndx,LOMoments_charge);
	  r_p[i] 		= 	curHOMoments->get_val(i,0,0,sIndx,HOMoments_charge);
	  r_p_old[i]	= 	oldHOMoments->get_val(i,0,0,sIndx,HOMoments_charge);
	  s2_p[i] 		= 	curHOMoments->get_val(i,0,0,sIndx,HOMoments_S2xx);
	  s2_p_half[i]	= 	0.5*(s2_p[i] + oldHOMoments->get_val(i,0,0,sIndx,HOMoments_S2xx));
	  s2_pb_half[i]	=	s2_p_half[i]/r_p[i];
	  if (i == 0) {
		  ru[0] 	= curLOMomentsES->get_val(0,0,0,sIndx,LOMoments_currentx);
		  ru[nx] 	= ru[0];
		  ru_old[0] = oldLOMomentsES->get_val(0,0,0,sIndx,LOMoments_currentx);
		  ru_old[nx]= ru_old[0];
		  E[0] 		= curFieldData->getE(0,0,0,0);
		  E[nx] 	= E[0];
		  E_old[0] 	= oldFieldData->getE(0,0,0,0);
		  E_old[nx] = E_old[0];
	  }
	  else {
		  ru[i] 	= curLOMomentsES->get_val(i,0,0,sIndx,LOMoments_currentx);
		  ru_old[i] = oldLOMomentsES->get_val(i,0,0,sIndx,LOMoments_currentx);
		  E[i] 		= curFieldData->getE(i,0,0,0);
		  E_old[i] 	= oldFieldData->getE(i,0,0,0);
	  }
  }

  // Construct half time face quantities
  realkind* r_f 	= (realkind*)malloc((nx+1)*sizeof(realkind));
  realkind* r_half 	= (realkind*)malloc(nx*sizeof(realkind));
  realkind* r_f_half= (realkind*)malloc((nx+1)*sizeof(realkind));
  realkind* E_half 	= (realkind*)malloc((nx+1)*sizeof(realkind));
  realkind* gammaNU = (realkind*)malloc((nx+1)*sizeof(realkind));
  for (int i = 0; i < pData->nx; i++) {
	  r_half[i] 	= (r[i] + r_old[i])*0.5;
  }
  for (int i = 0; i < pData->nx+1; i++) {
	  if (i == 0 || i == nx) {
		  r_f_half[0] 	= 0.5*(r_half[0] + r_half[nx-1]);
		  r_f_half[nx] 	= r_f_half[0];
		  r_f[0] 		= 0.5*(r[0] + r[nx-1]);
		  r_f[nx] 		= r_f[0];
		  gammaNU[0] 	= consistencyTermES->get_val(0,0,0,sIndx,consistencyterm_momentum);
		  gammaNU[nx] 	= gammaNU[0];
	  }
	  else {
		  r_f_half[i] 	= 0.5*(r_half[i] + r_half[i-1]);
		  r_f[i] 		= 0.5*(r[i] + r[i-1]);
		  gammaNU[i] 	= consistencyTermES->get_val(i,0,0,sIndx,consistencyterm_momentum);
	  }
  }
  for (int i = 0; i < pData->nx+1; i++) {
	  E_half[i] 	= (E[i] + E_old[i])*0.5;
  }
//  cout << endl;
  // Construct the non-linear residual
    for (int fIndx = 0; fIndx < this->pData->nx+1; fIndx++) {

      if (fIndx == 0 || fIndx == this->pData->nx) {
        this->F_ru[sIndx][fIndx] =
            (mdt_recip * (ru[fIndx] - ru_old[fIndx]) +
             mdx_recip * ((r[0] * s2_pb_half[0]) - 
                          (r[nx-1] * s2_pb_half[nx-1])) -
             (q * r_f_half[fIndx] * E_half[fIndx]) - 
             (gammaNU[fIndx] * r_f[fIndx]));
      }
      else {
        this->F_ru[sIndx][fIndx] =
            (mdt_recip * (ru[fIndx] - ru_old[fIndx]) +
             mdx_recip * ((r[fIndx] * s2_pb_half[fIndx]) - 
                          (r[fIndx - 1] * s2_pb_half[fIndx - 1])) -
             (q * r_f_half[fIndx] * E_half[fIndx]) - 
             (gammaNU[fIndx] * r_f[fIndx]));
      }
//      printf("F_ru[%i] = %e \n",fIndx,F_ru[sIndx][fIndx]);
//      printf("r_f_half[%i] = %f \n",fIndx,r_f_half[fIndx]);
//      printf("r_f[%i] = %f \n",fIndx,r_f[fIndx]);
//      printf("ru[%i] = %f \n",fIndx,ru[fIndx]);
//      printf("E_half[%i] = %f \n",fIndx,E_half[fIndx]);
//      printf("gammaNU[%i] = %f \n",fIndx,gammaNU[fIndx]);
//      printf("gammaNU[%i] = %f \n",fIndx,gammaNU[fIndx]);
    }
    free(r);
    free(ru);
    free(E);
    free(r_old);
    free(ru_old);
    free(E_old);
    free(r_p);
    free(r_p_old);
    free(s2_p);
    free(s2_p_half);
    free(s2_pb_half);
    free(r_f);
    free(r_f_half);
    free(r_half);
    free(E_half);
    free(gammaNU);
}

///////////////////////////////////////////////////////////////////////
//
// Return the maximum L2 of the array
//
///////////////////////////////////////////////////////////////////////

double NLResidual::maxL2(double* F, int size)
{
  double sumVal = 0.0;
  for (int i = 0; i < size; i++) {
    sumVal += (F[i] * F[i]);
  }
  double maxVal = sqrt(this->pData->dxdi * sumVal);

  return maxVal;
}

///////////////////////////////////////////////////////////////////////
//
// Return the absolute value maximum of the array
//
///////////////////////////////////////////////////////////////////////

double NLResidual::absoluteMax(double* F, int size)
{
  double maxVal = 0.0;
  for (int i = 0; i < size; i++) {
    if (maxVal < fabs(F[i]))
      maxVal = fabs(F[i]);
  }

  return maxVal;
}

///////////////////////////////////////////////////////////////////////
//
// Print the residual maxima
//
///////////////////////////////////////////////////////////////////////

void NLResidual::printResidual()
{
   cout << "Maximum L2 of NL residual is " << this->totMaxResidual << endl;
  for (int sIndx = 0; sIndx < this->pData->nspecies; sIndx++)
    cout << "Fr" << sIndx+1 << " = " << this->maxResidual_r[sIndx] << endl;
  for (int sIndx = 0; sIndx < this->pData->nspecies; sIndx++)
    cout << "Fru" << sIndx+1 << " = " << this->maxResidual_ru[sIndx] << endl;
  cout << "FE = " << this->maxResidual_E << endl;
}
