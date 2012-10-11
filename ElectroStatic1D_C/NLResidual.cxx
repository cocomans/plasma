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
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"
#include "DhatConsistency.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

NLResidual::NLResidual(PlasmaData* data)
{
  this->data = data;

  this->F_r = new double*[this->data->p_size];
  this->F_ru = new double*[this->data->p_size];

  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    // Cell based (differs from old version)
    this->F_r[sIndx] = new double[this->data->nx];

    // Face based (differs from old version)
    this->F_ru[sIndx] = new double[this->data->nfx];
  }
  this->F_E = new double[this->data->nfx];

  this->maxResidual_r = new double[this->data->p_size];
  this->maxResidual_ru = new double[this->data->p_size];
}

NLResidual::~NLResidual()
{
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
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
                        MomentSolution* curMoment,
                        MomentSolution* oldMoment,
                        ParticleSolution* curParticle,
                        ParticleSolution* oldParticle,
                        DhatConsistency* Dhat)
{
  nlres_E(curMoment, oldMoment);
  this->maxResidual_E = maxL2(this->F_E, this->data->nfx);
  this->totMaxResidual = this->maxResidual_E;

  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    nlres_r(sIndx, curMoment, oldMoment, Dhat);
    nlres_ru(sIndx, curMoment, oldMoment, curParticle, oldParticle, Dhat);

    this->maxResidual_r[sIndx] = maxL2(this->F_r[sIndx], this->data->nx);
    if (this->totMaxResidual < this->maxResidual_r[sIndx])
      this->totMaxResidual = this->maxResidual_r[sIndx];
    
    this->maxResidual_ru[sIndx] = maxL2(this->F_ru[sIndx], this->data->nfx);
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

void NLResidual::nlres_E(
                        MomentSolution* curMoment,
                        MomentSolution* oldMoment)
{
  double* E = curMoment->getE();
  double* E_old = oldMoment->getE();
  double wp2dt = 1.0 / (this->data->wp * this->data->wp * this->data->dt);

  // Sum all the currents from different species
  double* j_tot = new double[this->data->nfx];
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    j_tot[fIndx] = 0.0;

  if (this->data->temp_order == ORDER_1_BACK_EULER) {
    for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
      double* ru = curMoment->getRU(sIndx);
      for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
        j_tot[fIndx] += this->data->q[sIndx] * ru[fIndx];
    }
  }

  else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
    for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
      double* ru_half = curMoment->getRU(sIndx);
      for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
        j_tot[fIndx] += this->data->q[sIndx] * ru_half[fIndx];
    }
  }

  double j_avg = 0.0;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    j_avg += j_tot[cIndx];
  j_avg /= this->data->nx;

  // Construct the non-linear residual for the Ampere Law equation on node
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    this->F_E[fIndx] = wp2dt * (E[fIndx] - E_old[fIndx]) + j_tot[fIndx] - j_avg;

/*
  cout << "wp2dt " << wp2dt << " javg " << this->data->j_avg << endl;
  cout << "E from nlres_E" << endl;
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    cout << E[fIndx] << endl;
  cout << "E_old from nlres_E" << endl;
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    cout << E_old[fIndx] << endl;
  cout << "ru_tot from nlres_E" << endl;
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    cout << ru_tot[fIndx] << endl;
  cout << "F_E from nlres_E" << endl;
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    cout << this->F_E[fIndx] << endl;
*/
  delete [] j_tot;
}

///////////////////////////////////////////////////////////////////////
//
// Construct the non-linear residual function for continuity equation
//
///////////////////////////////////////////////////////////////////////

void NLResidual::nlres_r(
                        int sIndx,
                        MomentSolution* curMoment,
                        MomentSolution* oldMoment,
                        DhatConsistency* Dhat)
{
  double m = this->data->m[sIndx];
  double mdt_recip = m / this->data->dt;
  double mdx_recip = m / this->data->dx;

  // Cell based
  double* r = curMoment->getR(sIndx);
  double* r_old = oldMoment->getR(sIndx);

  // Face based
  double* ru = curMoment->getRU(sIndx);

  // Consistency
  double* Dhat_r = Dhat->getR(sIndx);

/*
cout << "m " << m << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_r r " << r[cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_r r_old " << r_old[cIndx] << endl;
for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
cout << "nlres_r ru " << ru[fIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_r Dhat_r " << Dhat_r[cIndx] << endl;
*/

  // Construct the non-linear residual for the continuity equation on node
  // Cell based and no periodic condition
  if (this->data->temp_order == ORDER_1_BACK_EULER) {
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
      this->F_r[sIndx][cIndx] =
               (mdt_recip * (r[cIndx] - r_old[cIndx]) +
                mdx_recip * (ru[cIndx + 1] - ru[cIndx])) / m - 
               Dhat_r[cIndx] * r[cIndx] / m;
    }
  }
  else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
      this->F_r[sIndx][cIndx] =
               (mdt_recip * (r[cIndx] - r_old[cIndx]) +
                mdx_recip * (ru[cIndx + 1] - ru[cIndx])) / m;
    }
    if (this->data->cc_flag == CHARGE_CONSERVE_OFF) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        this->F_r[sIndx][cIndx] -= Dhat_r[cIndx] * r[cIndx] / m;
      }
    }
  }
/*
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_r F_r " << this->F_r[sIndx][cIndx] << endl;
*/
}

///////////////////////////////////////////////////////////////////////
//
// Construct the non-linear residual function for momentum equation
//
///////////////////////////////////////////////////////////////////////

void NLResidual::nlres_ru(
                        int sIndx,
                        MomentSolution* curMoment,
                        MomentSolution* oldMoment,
                        ParticleSolution* curParticle,
                        ParticleSolution* oldParticle,
                        DhatConsistency* Dhat)
{
  double q = this->data->q[sIndx];
  double m = this->data->m[sIndx];
  double mdt_recip = m / this->data->dt;
  double mdx_recip = m / this->data->dx;

  // Extract moment quantities at time n+1 and n+1/2
  double* r = curMoment->getR(sIndx);
  double* ru = curMoment->getRU(sIndx);
  double* E = curMoment->getE();

  // Extract moment quantities at time n
  double* r_old = oldMoment->getR(sIndx);
  double* ru_old = oldMoment->getRU(sIndx);
  double* E_old = oldMoment->getE();

  // Extract particle quantities at time n+1 and n+1/2
  double* r_p = curParticle->getR(sIndx);
  double* s2_p = curParticle->getS2(sIndx);
  double* s2_p_half = curParticle->getS2Half(sIndx);

  // Extract particle quantities at time n
  double* r_p_old = oldParticle->getR(sIndx);

  // Construct density normalized stress tensor
  double* s2_pb = new double[this->data->nx];
  double* s2_pb_half = new double[this->data->nx];
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    s2_pb[cIndx ] = s2_p[cIndx] / r_p[cIndx];
    s2_pb_half[cIndx ] = s2_p_half[cIndx] / r_p[cIndx];
  }

  // Construct half time face quantities
  double* r_half = new double[this->data->nx];
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    r_half[cIndx] = (r[cIndx] + r_old[cIndx]) / 2.0;
  }
  double* r_f_half = new double[this->data->nfx];
  makeFaceFlux(r_half, r_f_half, this->data->nx);

  // Construct full time face quantities
  double* r_f = new double[this->data->nfx];
  makeFaceFlux(r, r_f, this->data->nx);

  // Construct half time electric field
  double* E_half = new double[this->data->nfx];
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    E_half[fIndx] = (E[fIndx] + E_old[fIndx]) / 2.0;
  }

  // Consistency
  double* Dhat_ru = Dhat->getRU(sIndx);

/*
for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
cout << "nlres_ru ru " << ru[fIndx] << endl;
for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
cout << "nlres_ru ru_old " << ru_old[fIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_ru r " << r[cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_ru r_p " << r_p[cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_ru s2_p " << s2_p[cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_ru s2_p_half " << s2_p_half[cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "nlres_ru s2_pb_half " << s2_pb_half[cIndx] << endl;
for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
cout << "nlres_ru r_f_half " << r_f_half[fIndx] << endl;
for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
cout << "nlres_ru E_half " << E_half[fIndx] << endl;
for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
cout << "nlres_ru Dhat_ru " << Dhat_ru[fIndx] << endl;
for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
cout << "nlres_ru r_f " << r_f[fIndx] << endl;
*/
//Will Taitano edit: 06/25/2012 for debugging purpose
/*    for (int fIndx = 0; fIndx < this->data->nx; fIndx++) {
        Dhat_ru[fIndx] = 0.0;
    }*/
  // Construct the non-linear residual
  int last = this->data->nx - 1;
  if (this->data->temp_order == ORDER_1_BACK_EULER) {
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {

      if (fIndx == 0 || fIndx == this->data->nx) {
        this->F_ru[sIndx][fIndx] =
            (mdt_recip * (ru[fIndx] - ru_old[fIndx]) +
             mdx_recip * ((r[0] * s2_pb[0]) - 
                          (r[last] * s2_pb[last])) -
             (q * r_f[fIndx] * E[fIndx]) - 
             (Dhat_ru[fIndx] * r_f[fIndx])) / m;
      }
      else {
        this->F_ru[sIndx][fIndx] =
            (mdt_recip * (ru[fIndx] - ru_old[fIndx]) +
             mdx_recip * ((r[fIndx] * s2_pb[fIndx]) - 
                          (r[fIndx - 1] * s2_pb[fIndx - 1])) -
             (q * r_f[fIndx] * E[fIndx]) - 
             (Dhat_ru[fIndx] * r_f[fIndx])) / m;
      }
    }
  }
  else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {

      if (fIndx == 0 || fIndx == this->data->nx) {
        this->F_ru[sIndx][fIndx] =
            (2.0 * mdt_recip * (ru[fIndx] - ru_old[fIndx]) +
             mdx_recip * ((r[0] * s2_pb_half[0]) - 
                          (r[last] * s2_pb_half[last])) -
             (q * r_f_half[fIndx] * E_half[fIndx]) - 
             (Dhat_ru[fIndx] * r_f[fIndx])) / m;
      }
      else {
        this->F_ru[sIndx][fIndx] =
            (2.0 * mdt_recip * (ru[fIndx] - ru_old[fIndx]) +
             mdx_recip * ((r[fIndx] * s2_pb_half[fIndx]) - 
                          (r[fIndx - 1] * s2_pb_half[fIndx - 1])) -
             (q * r_f_half[fIndx] * E_half[fIndx]) - 
             (Dhat_ru[fIndx] * r_f[fIndx])) / m;
      }
    }
  }
//Will Taitano edit: 06/24/2012 for debugging purpose
/*
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    cout << "nlres_ru " << this->F_ru[sIndx][fIndx] << endl;
*/

  delete [] r_half;
  delete [] E_half;
  delete [] s2_pb;
  delete [] s2_pb_half;
  delete [] r_f;
  delete [] r_f_half;
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
  double maxVal = sqrt(this->data->dx * sumVal);

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
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    cout << "Fr" << sIndx+1 << " = " << this->maxResidual_r[sIndx] << endl;
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    cout << "Fru" << sIndx+1 << " = " << this->maxResidual_ru[sIndx] << endl;
  cout << "FE = " << this->maxResidual_E << endl;
}
