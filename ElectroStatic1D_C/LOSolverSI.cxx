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

#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"
#include "NLResidual.h"
#include "LOSolverSI.h"
#include "DhatConsistency.h"
#include "stdio.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

LOSolverSI::LOSolverSI(PlasmaData* data) : LOSolver(data)
{
}

LOSolverSI::~LOSolverSI()
{
}

///////////////////////////////////////////////////////////////////////
//
// Solver high level
//
///////////////////////////////////////////////////////////////////////

int LOSolverSI::solve(
                     MomentSolution* curMoment,
                     MomentSolution* oldMoment,
                     ParticleSolution* curParticle,
                     ParticleSolution* oldParticle,
                     DhatConsistency* Dhat)
{
  this->curParticle = curParticle;
  this->oldParticle = oldParticle;
  this->curMoment = curMoment;
  this->oldMoment = oldMoment;
  this->Dhat = Dhat;

  int k_inner;
  if (this->data->p_size == 1)
    k_inner = siSolverElectron();
  else
    k_inner = siSolverTwoSpecies();

  return k_inner;
}

///////////////////////////////////////////////////////////////////////
//
// Semi-implicit low order solver for electron only
//
///////////////////////////////////////////////////////////////////////

int LOSolverSI::siSolverElectron()
{
  cout << "LOSolverSI::siSolverElectron()" << endl;
  double* E_source = new double[this->data->nfx];
  double* E_factor = new double[this->data->nfx];
  double* E_factor_half = new double[this->data->nfx];
  double* s2e_pb = new double[this->data->nx];
  double* s2e_pb_half = new double[this->data->nx];
  double* re_f = new double[this->data->nfx];
  double* re_f_old = new double[this->data->nfx];
  double* re_f_half = new double[this->data->nfx];

  double* Dhat_re = this->Dhat->getR(ELECTRON);
  double* Dhat_rue = this->Dhat->getRU(ELECTRON);

  double dt = this->data->dt;
  double wps = this->data->wp * this->data->wp;
  double j_avg = this->data->j_avg;
  double dx_recip = 1.0 / this->data->dx;
  double dt2_recip = this->data->dt / 2.0;
  int lo_trunc = 100;

  // Electron charge and mass
  double qe = this->data->q[ELECTRON];
  double me = this->data->m[ELECTRON];

  // Calculate E factor for either order solution
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
      E_factor[fIndx] = 1.0;
      E_factor_half[fIndx] = 0.5;
  }

  // Electric fields
  double* E = this->curMoment->getE();
  double* E_old = this->oldMoment->getE();

  // Moment quantities at n+1 and n+1/2 time
  double* re = this->curMoment->getR(ELECTRON);        // electron density
  double* rue = this->curMoment->getRU(ELECTRON);      // electron momentum

  // Moment quantities at n time
  double* re_old = this->oldMoment->getR(ELECTRON);    // electron density
  double* rue_old = this->oldMoment->getRU(ELECTRON);  // electron momentum

  // Kinetic moments at n+1 and n+1/2 time
  double* re_p = this->curParticle->getR(ELECTRON);
  double* rue_half_p = this->curParticle->getRUHalf(ELECTRON);
  double* s2e_p_half = this->curParticle->getS2Half(ELECTRON);
  double* s2e_p = this->curParticle->getS2(ELECTRON);
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    s2e_pb[cIndx] = s2e_p[cIndx] / re_p[cIndx];
    s2e_pb_half[cIndx] = s2e_p_half[cIndx] / re_p[cIndx];
  }

  // Calculate initial non-linear residual for low order system for tolerance
  NLResidual* residual = new NLResidual(this->data);
  double maxResidual0, maxResidual;
  maxResidual0 = residual->calculateNLResidual(
                                this->curMoment, this->oldMoment,
                                this->curParticle, this->oldParticle,
                                this->Dhat);
  cout << endl <<  "   The initial inner: " << maxResidual0 << endl;
  residual->printResidual();

  // Convergence loop
  int flag_conv = 0;
  int k_inner = 0;
  while (flag_conv == 0) {
    double j_avg = 0.0;
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      j_avg += qe * rue[cIndx];
    j_avg /= this->data->nx;

    k_inner++;

    // Calculate E source
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      E_source[fIndx] = E_old[fIndx] + (dt * wps * j_avg);

    // Calculate electron density face quantity for n+1/2, p+1
    makeFaceFlux(re_old, re_f_old, this->data->nx);
    makeFaceFlux(re, re_f, this->data->nx);
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      re_f_half[fIndx] = (re_f[fIndx] + re_f_old[fIndx]) / 2.0;

    // Linear solver returning rue
    if (this->data->temp_order == ORDER_1_BACK_EULER) {
      solveLinearSystem(dt, s2e_pb, re_f, re_f, E_factor,
                        re, re_old, rue_old, E_old,
                        E_source, Dhat_re, Dhat_rue, rue);
    }
    else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
      solveLinearSystem(dt2_recip, s2e_pb_half, re_f_half, re_f, E_factor_half,
                        re, re_old, rue_old, E_old,
                        E_source, Dhat_re, Dhat_rue, rue);
    }

    // Back substitute to calculate electron continuity
    if (this->data->cc_flag == CHARGE_CONSERVE_ON) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
        re[cIndx] = re_old[cIndx] - 
                    (dt * dx_recip * (rue[cIndx + 1] - rue[cIndx]));
    }
    else if (this->data->cc_flag == CHARGE_CONSERVE_OFF) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
        re[cIndx] = (re_old[cIndx] - 
                     (dt * dx_recip * (rue[cIndx + 1] - rue[cIndx]))) /
                    (1.0 - (dt * Dhat_re[cIndx] / me));
    }

    // Back substitute to calculate electric field
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      E[fIndx] = E_source[fIndx] - (wps * dt * qe * rue[fIndx]);

    // Check for convergence by calculating the L2 of non-linear residual
    maxResidual = residual->calculateNLResidual(
                                this->curMoment, this->oldMoment,
                                this->curParticle, this->oldParticle,
                                this->Dhat);
    cout << endl << "   k_inner = " << k_inner << " Picard iteration: " 
         << maxResidual << endl;
    residual->printResidual();

    // Check for convergence in the Picard loop (particle solution)
    if (maxResidual <= (this->data->tol_nlres_inner * maxResidual0) || 
        k_inner == lo_trunc)
        flag_conv = 1;
  }
  delete residual;
  delete [] E_source;
  delete [] E_factor;
  delete [] E_factor_half;
  delete [] s2e_pb;
  delete [] s2e_pb_half;
  delete [] re_f;
  delete [] re_f_old;
  delete [] re_f_half;

  return k_inner;
}

///////////////////////////////////////////////////////////////////////
//
// Semi-implicit low order solver for two species
//
///////////////////////////////////////////////////////////////////////

int LOSolverSI::siSolverTwoSpecies()
{
  cout << "LOSolverSI::siSolverTwoSpecies()" << endl;
  double* E_source = new double[this->data->nfx];
  double* E_source_half = new double[this->data->nfx];
  double* E_factor = new double[this->data->nfx];
  double* E_factor_half = new double[this->data->nfx];

  double* re_f = new double[this->data->nfx];
  double* re_f_old = new double[this->data->nfx];
  double* re_f_half = new double[this->data->nfx];
  double* s2e_pb = new double[this->data->nx];
  double* s2e_pb_half = new double[this->data->nx];

  double* ri_f = new double[this->data->nfx];
  double* ri_half = new double[this->data->nfx];
  double* ri_f_old = new double[this->data->nfx];
  double* ri_f_half = new double[this->data->nfx];
  double* rui_source = new double[this->data->nfx];
  double* rui_source_half = new double[this->data->nfx];
  double* s2i_pb = new double[this->data->nx];
  double* s2i_pb_half = new double[this->data->nx];

  double* Dhat_re = this->Dhat->getR(ELECTRON);
  double* Dhat_rue = this->Dhat->getRU(ELECTRON);
  double* Dhat_ri = this->Dhat->getR(ION);
  double* Dhat_rui = this->Dhat->getRU(ION);
//Will Taitano edit: 06/24/2012 The below is for debugging purpose
/*    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        Dhat_re[cIndx] = 0.0;
        Dhat_ri[cIndx] = 0.0;
    }
    for (int fIndx = 0; fIndx < this->data->nx+1; fIndx++) {
        Dhat_rue[fIndx] = 0.0;
        Dhat_rui[fIndx] = 0.0;
    }*/
  double dt = this->data->dt;
  double wps = this->data->wp * this->data->wp;
  double j_avg = this->data->j_avg;
  double dx_recip = 1.0 / this->data->dx;
  double dt2_recip = this->data->dt / 2.0;
  double dts = this->data->dt * this->data->dt;
  int lo_trunc = 100;

  // Electron and ion charge and mass
  double qe = this->data->q[ELECTRON];
  double qi = this->data->q[ION];
  double me = this->data->m[ELECTRON];
  double mi = this->data->m[ION];

  // Full time and old time electric fields
  double* E = this->curMoment->getE();
  double* E_old = this->oldMoment->getE();

  // Moment quantities at n+1 and n+1/2 time
  double* re = this->curMoment->getR(ELECTRON);        // electron density
  double* rue = this->curMoment->getRU(ELECTRON);      // electron momentum
  double* ri = this->curMoment->getR(ION);             // ion density
  double* rui = this->curMoment->getRU(ION);           // ion momentum

  // Moment quantities at n time
  double* re_old = this->oldMoment->getR(ELECTRON);    // electron density
  double* rue_old = this->oldMoment->getRU(ELECTRON);  // electron momentum
  double* ri_old = this->oldMoment->getR(ION);         // ion density
  double* rui_old = this->oldMoment->getRU(ION);       // ion momentum

  // Kinetic moments at n+1 and n+1/2 time
  double* re_p = this->curParticle->getR(ELECTRON);
  double* rue_half_p = this->curParticle->getRUHalf(ELECTRON);
  double* s2e_p_half = this->curParticle->getS2Half(ELECTRON);
  double* s2e_p = this->curParticle->getS2(ELECTRON);
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    s2e_pb[cIndx] = s2e_p[cIndx] / re_p[cIndx];
    s2e_pb_half[cIndx] = s2e_p_half[cIndx] / re_p[cIndx];
  }

  double* ri_p = this->curParticle->getR(ION);
  double* rui_half_p = this->curParticle->getRUHalf(ION);
  double* s2i_p_half = this->curParticle->getS2Half(ION);
  double* s2i_p = this->curParticle->getS2(ION);
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    s2i_pb[cIndx] = s2i_p[cIndx] / ri_p[cIndx];
    s2i_pb_half[cIndx] = s2i_p_half[cIndx] / ri_p[cIndx];
  }

  // Calculate initial non-linear residual for low order system for tolerance
  NLResidual* residual = new NLResidual(this->data);
  double maxResidual0, maxResidual;
  maxResidual0 = residual->calculateNLResidual(
                                this->curMoment, this->oldMoment,
                                this->curParticle, this->oldParticle,
                                this->Dhat);
  cout << endl <<  "   The initial inner: " << maxResidual0 << endl;
  residual->printResidual();

  // Convergence loop
  int flag_conv = 0;
  int k_inner = 0;
  while (flag_conv == 0) 
  {
    double j_avg = 0.0;
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      j_avg += qe * rue[cIndx] + qi * rui[cIndx];
    j_avg /= this->data->nx;

    k_inner++;

    // Calculate ion n+1, p+1 ion density
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      ri[cIndx] = (ri_old[cIndx] -
                   dt * dx_recip * (rui[cIndx + 1] - rui[cIndx]));//
//                  / (1.0 - dt * Dhat_ri[cIndx] / mi);

//    if (this->data->temp_order == ORDER_2_CRANK_NICOLSON)
      if (this->data->cc_flag == CHARGE_CONSERVE_OFF) {
          for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
              ri[cIndx] /= (1.0 - dt * Dhat_ri[cIndx] / mi);          
      }

    // Calculate face quantity for n+1, p+1 and n+1/2, p+1 ion density
    makeFaceFlux(ri, ri_f, this->data->nx);
    makeFaceFlux(ri_old, ri_f_old, this->data->nx);
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      ri_f_half[fIndx] = (ri_f[fIndx] + ri_f_old[fIndx]) / 2.0;

    // Source term from Ampere equation
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
      E_factor[fIndx] = 
               1.0 / (1.0 + wps * dts * qi *qi * ri_f[fIndx] / mi);
      E_factor_half[fIndx] = 
               1.0 / (1.0 + wps * dts * qi *qi * ri_f_half[fIndx] / mi / 4.0);
    }

    // Source term from ion momentum equation 
    int last = this->data->nx - 1;
    if (this->data->temp_order == ORDER_1_BACK_EULER) {
      for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
        if (fIndx == 0 || fIndx == this->data->nx) {
          rui_source[fIndx] = rui_old[fIndx] -
               dt * dx_recip * 
               (ri[0] * s2i_pb[0] - ri[last] * s2i_pb[last]) +
               dt * Dhat_rui[fIndx] * ri_f[fIndx] / mi;
        }
        else {
          rui_source[fIndx] = rui_old[fIndx] -
               dt * dx_recip * 
               (ri[fIndx] * s2i_pb[fIndx] - ri[fIndx - 1] * s2i_pb[fIndx - 1]) +
               dt * Dhat_rui[fIndx] * ri_f[fIndx] / mi;
        }
      }
      for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
        E_source[fIndx] = (E_old[fIndx] -
                         (wps * dt * qi * rui_source[fIndx]) +
                         (wps * dt * qi * j_avg))*E_factor[fIndx];
//                         (wps * dt * qi * this->data->j_avg)) *
//                        E_factor[fIndx];
      }
    }
    else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
      for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
        if (fIndx == 0 || fIndx == this->data->nx) {
          rui_source_half[fIndx] = rui_old[fIndx] -
               dt2_recip * dx_recip * 
               (ri[0] * s2i_pb_half[0] - ri[last] * s2i_pb_half[last]) +
               dt2_recip * Dhat_rui[fIndx] * ri_f[fIndx] / mi +
               dt2_recip * qi * ri_f_half[fIndx] * E_old[fIndx] / 2.0 / mi;
        }
        else {
          rui_source_half[fIndx] = rui_old[fIndx] -
               dt2_recip * dx_recip * 
               (ri[fIndx] * s2i_pb_half[fIndx] - ri[fIndx - 1] * s2i_pb_half[fIndx - 1]) +
               dt2_recip * Dhat_rui[fIndx] * ri_f[fIndx] / mi +
               dt2_recip * qi * ri_f_half[fIndx] * E_old[fIndx] / 2.0 / mi;
        }
      }
      // Source term from ampere equation
      for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
        E_source_half[fIndx] = (E_old[fIndx] -
                         (wps * dt * qi * rui_source_half[fIndx]) +
                         (wps * dt * qi * j_avg)) * E_factor_half[fIndx];
//                         (wps * dt * qi * this->data->j_avg)) *
//                        E_factor_half[fIndx];
      }
    }
    // Calculate face quantity for n+1/2, p+1 electron density
    makeFaceFlux(re_old, re_f_old, this->data->nx);
    makeFaceFlux(re, re_f, this->data->nx);

    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      re_f_half[fIndx] = (re_f[fIndx] + re_f_old[fIndx]) / 2.0;
/*      for (int fIndx = 0; fIndx < this->data->nx+1; fIndx++) {
          cout << "re_f[" << fIndx << "] = " << re_f[fIndx] << endl;
      }    */
    // Solve the linear system
    if (this->data->temp_order == ORDER_1_BACK_EULER) {
      solveLinearSystem(dt, s2e_pb, re_f, re_f, E_factor,
                        re, re_old, rue_old, E_old,
                        E_source, Dhat_re, Dhat_rue, rue);
    }
    else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
      solveLinearSystem(dt2_recip, s2e_pb_half, re_f_half, re_f, E_factor_half,
                        re, re_old, rue_old, E_old,
                        E_source_half, Dhat_re, Dhat_rue, rue);

    }

    // Back substitute to calculate electron continuity
    if (this->data->cc_flag == CHARGE_CONSERVE_ON) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
          re[cIndx] = re_old[cIndx] - 
                      (dt * dx_recip * (rue[cIndx + 1] - rue[cIndx]));
    }
    else if (this->data->cc_flag == CHARGE_CONSERVE_OFF) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
        re[cIndx] = (re_old[cIndx] - 
                     dt * dx_recip * (rue[cIndx + 1] - rue[cIndx])) /
                    (1.0 - (dt * Dhat_re[cIndx] / me));
    }
      if (this->data->temp_order == ORDER_1_BACK_EULER) {
          // Back substitute to calculate electric field
          for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
              E[fIndx] = E_source[fIndx] - 
              (wps * dt * qe * rue[fIndx] * E_factor[fIndx]);
          
          for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
              rui[fIndx] = rui_source[fIndx] +
              (dt * qi * ri_f_half[fIndx] * E[fIndx]) /4.0 / mi;          
      }
      else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
          // Back substitute to calculate electric field
          for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
              E[fIndx] = E_source_half[fIndx] - 
              (wps * dt * qe * rue[fIndx] * E_factor_half[fIndx]);
          
          for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
              rui[fIndx] = rui_source_half[fIndx] +
              (dt * qi * ri_f_half[fIndx] * E[fIndx]) /4.0/ mi;
      }
      


    // Check for convergence by calculating the L2 of non-linear residual
    maxResidual = residual->calculateNLResidual(
                                this->curMoment, this->oldMoment,
                                this->curParticle, this->oldParticle,
                                this->Dhat);
    cout << endl << "   k_inner = " << k_inner << " Picard iteration: " 
         << maxResidual << endl;
    residual->printResidual();

    // Check for convergence in the Picard loop (particle solution)
    if (maxResidual <= this->data->tol_nlres_inner * maxResidual0 || 
        k_inner == lo_trunc)
          flag_conv = 1;
  }
  delete residual;
  delete [] E_source; 
  delete [] E_source_half; 
  delete [] E_factor; 
  delete [] E_factor_half; 

  delete [] re_f;
  delete [] re_f_old;
  delete [] re_f_half;
  delete [] s2e_pb;
  delete [] s2e_pb_half;

  delete [] ri_f;
  delete [] ri_f_old;
  delete [] ri_f_half;
  delete [] rui_source;
  delete [] rui_source_half;
  delete [] s2i_pb;
  delete [] s2i_pb_half;

  return k_inner;
}

///////////////////////////////////////////////////////////////////////
//
// Triadiagonal linear solver for two species electron and ion
// First order has dt_factor = dt and uses s2e_pb and re_f
// Second order has dt_factor = dt/2 and uses s2e_pb_half and re_f_half
//
///////////////////////////////////////////////////////////////////////

void LOSolverSI::solveLinearSystem(
                     double dt_factor,   // solver order dependent input
                     double* s2e,        // solver order dependent input
                     double* re_f,       // solver order dependent input
                     double* re_f_full,  // fixed input, full time electron number density
                     double* E_factor,   // solver order dependent input
                     double* re,         // input
                     double* re_old,     // input
                     double* rue_old,    // input
                     double* E_old,      // input
                     double* E_source,   // input
                     double* Dhat_re,    // input
                     double* Dhat_rue,   // input
                     double* rue)        // output
{
  double me = this->data->m[ELECTRON];
  double qe = this->data->q[ELECTRON];
  double dt = this->data->dt;
  double dx_recip = 1.0 / this->data->dx;
  double dxs_recip = 1.0 / this->data->dx / this->data->dx;
  double wps = this->data->wp * this->data->wp;
  double t_recip = 1.0/2.0;
  double** A = new double*[this->data->nfx];
  double* rhs = new double[this->data->nfx];
  for (int i = 0; i < this->data->nfx; i++) {
    A[i] = new double[this->data->nfx];
    for (int j = 0; j < this->data->nfx; j++)
      A[i][j] = 0.0;
  }
  // Form the coefficient matrix for electron momentum equation
  int last = this->data->nx - 1;
    if (this->data->temp_order == ORDER_1_BACK_EULER) {
        for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
            if (fIndx == 0 || fIndx == this->data->nx) {
                A[fIndx][fIndx] = 
                (me / dt_factor) +
                (me * dxs_recip * dt * s2e[0]) + 
                (me * dxs_recip * dt * s2e[last]) + 
                (qe * qe * wps * dt * re_f[fIndx] * E_factor[fIndx]);
                A[fIndx][1] = 
                -me * dt * dxs_recip * s2e[0];
                A[fIndx][last] = 
                -me * dt * dxs_recip * s2e[last];
            }
            else {
                A[fIndx][fIndx] = 
                (me / dt_factor) +
                (me * dxs_recip * dt * s2e[fIndx]) + 
                (me * dxs_recip * dt * s2e[fIndx - 1]) + 
                (qe * qe * wps * dt * re_f[fIndx] * E_factor[fIndx]);
                A[fIndx][fIndx + 1] = 
                -me * dt * dxs_recip * s2e[fIndx];
                A[fIndx][fIndx - 1] = 
                -me * dt * dxs_recip * s2e[fIndx - 1];
            }
        }
    }
    else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
/*solveLinearSystem(dt2_recip, s2e_pb_half, re_f_half, re_f, E_factor_half,
 re, re_old, rue_old, E_old,
 E_source_half, Dhat_re, Dhat_rue, rue);*/       
        for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
            if (fIndx == 0 || fIndx == this->data->nx) {
                A[fIndx][fIndx] = 
                (me / dt_factor) +
                (me * dxs_recip * dt * s2e[0]) + 
                (me * dxs_recip * dt * s2e[last]) + 
                t_recip*(qe * qe * wps * dt * re_f[fIndx] * E_factor[fIndx]);
                
                A[fIndx][1] = 
                -me * dt * dxs_recip * s2e[0];
                
                A[fIndx][last] = 
                -me * dt * dxs_recip * s2e[last];
            }
            else {
                A[fIndx][fIndx] = 
                (me / dt_factor) +
                (me * dxs_recip * dt * s2e[fIndx]) + 
                (me * dxs_recip * dt * s2e[fIndx - 1]) + 
                t_recip*(qe * qe * wps * dt * re_f[fIndx] * E_factor[fIndx]);
                
                A[fIndx][fIndx + 1] = 
                -me * dt * dxs_recip * s2e[fIndx];
                
                A[fIndx][fIndx - 1] = 
                -me * dt * dxs_recip * s2e[fIndx - 1];
            }
        }
    }
  // Right hand side
    if (this->data->temp_order == ORDER_1_BACK_EULER) {
        for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
            if (fIndx == 0 || fIndx == this->data->nx) {
                rhs[fIndx] = 
                (rue_old[fIndx] * (me / dt_factor)) -
                (me * dx_recip *
                 (s2e[0] * re_old[0] - s2e[last] * re_old[last])) +
                (qe * re_f[fIndx] * E_source[fIndx]) +
                (Dhat_rue[fIndx] * re_f_full[fIndx]);
            }
            else {
                rhs[fIndx] = 
                (rue_old[fIndx] * (me / dt_factor)) -
                (me * dx_recip *
                 (s2e[fIndx] * re_old[fIndx] - 
                  s2e[fIndx - 1] * re_old[fIndx - 1])) +
                (qe * re_f[fIndx] * E_source[fIndx]) +
                (Dhat_rue[fIndx] * re_f_full[fIndx]);
            }
        }
    }
    else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
        for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
            if (fIndx == 0 || fIndx == this->data->nx) {
                rhs[fIndx] = 
                (rue_old[fIndx] * (me / dt_factor)) -
                (me * dx_recip *
                 (s2e[0] * re_old[0] - s2e[last] * re_old[last])) +
                t_recip*(qe * re_f[fIndx] * E_source[fIndx]) +
                t_recip*(qe * re_f[fIndx] * E_old[fIndx]) +
                (Dhat_rue[fIndx] * re_f_full[fIndx]);

/*                    cout << "re_f_full[" << fIndx << "] = " << re_f_full[fIndx] << endl;
                cout << "Dhat_rue[" << fIndx << "] = " << Dhat_rue[fIndx] << endl;                */
            }
            else {
                rhs[fIndx] = 
                (rue_old[fIndx] * (me / dt_factor)) -
                (me * dx_recip *
                 (s2e[fIndx] * re_old[fIndx] - 
                  s2e[fIndx - 1] * re_old[fIndx - 1])) +
                t_recip*(qe * re_f[fIndx] * E_source[fIndx]) +
                t_recip*(qe * re_f[fIndx] * E_old[fIndx]) + 
                (Dhat_rue[fIndx] * re_f_full[fIndx]);

/*                    cout << "re_f_full[" << fIndx << "] = " << re_f_full[fIndx] << endl;
                cout << "Dhat_rue[" << fIndx << "] = " << Dhat_rue[fIndx] << endl;
 */
            }
        }        
    }

  // Charge conserving scheme is not used
  if (this->data->cc_flag == CHARGE_CONSERVE_OFF) {
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
      if (fIndx == 0 || fIndx == this->data->nx) {
        rhs[fIndx] -= me * dx_recip *
                      (Dhat_re[0] * s2e[0] * re[0] -
                       Dhat_re[last] * s2e[last] * re[last]) *
                      dt / me;
      }
      else {
        rhs[fIndx] -= me * dx_recip *
                      (Dhat_re[fIndx] * s2e[fIndx] * re[fIndx] -
                       Dhat_re[fIndx - 1] * s2e[fIndx - 1] * re[fIndx - 1]) *
                      dt / me;
      }
    }
  }

  // Calculate the new time electron momentum from direct inversion
  double residual;
  int flag;
  GaussElim(A, this->data->nfx, rhs, this->data->nfx, rue, residual, flag);

  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    delete [] A[fIndx];
  delete [] A;
  delete [] rhs;
}
