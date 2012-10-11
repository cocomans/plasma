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

#include "PushSimple.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"
// Will Taitano edit 06/25/2012: For debugging purpose
#include "stdio.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

PushSimple::PushSimple(PlasmaData* data) : ParticlePush(data)
{
}

PushSimple::~PushSimple()
{
}

///////////////////////////////////////////////////////////////////////
//
// Calculate the new position and velocity of particles after push
// Based on the old Picard point electric field, using subcycling of particles
// Moves the current ParticleSolution and MomentSolution and calculates new
//
///////////////////////////////////////////////////////////////////////

void PushSimple::push(
                         MomentSolution* curMoment,       // unchanged
                         MomentSolution* oldMoment,       // unchanged
                         ParticleSolution* curParticle,   // updated
                         ParticleSolution* oldParticle)   // unchanged
{
  cout << "PushSimple::particlePush()" << endl;
  double tol_picard = 1e-10;
  double dx_recip = 1.0 / this->data->dx;

  double* E = new double[this->data->nfx];
  double* E_old = new double[this->data->nfx];
  double* E_half = new double[this->data->nfx];

  // Copy the electric field from the solver because it gets smoothed
  // Calculate the average E face value between time steps
  double* E_from_moment = curMoment->getE();
  double* E_old_from_moment = oldMoment->getE();

  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    E[fIndx] = E_from_moment[fIndx];
    E_old[fIndx] = E_old_from_moment[fIndx];
    E_half[fIndx] = (E[fIndx] + E_old[fIndx]) / 2.0; 
  }

  // Smooth electric field
  for (int i = 0; i < this->data->fil_num; i++) {
    smoothFaceFilter(this->data->nfx, E_half);
    smoothFaceFilter(this->data->nfx, E);
  }

  // Calculate full time and half time moment solutions
  MomentSolution* full = new MomentSolution(this->data);
  MomentSolution* half = new MomentSolution(this->data);

  /////////////////////////////////////////////////////////////////
  //
  // Loop through all species
  //
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    double dt = this->data->sub_dt[sIndx];
    double qm = this->data->q[sIndx] / this->data->m[sIndx];
    int numberOfParticles = this->data->NP0_tot[sIndx];

    // Get the particle information from the solutions
    double* Xold = oldParticle->getX(sIndx);
    double* Vold = oldParticle->getV(sIndx);
    double* X = curParticle->getX(sIndx);
    double* V = curParticle->getV(sIndx);

    /////////////////////////////////////////////////////////////////
    //
    // Loop through all particles for the species
    //
    for (int pIndx = 0; pIndx < numberOfParticles; pIndx++) {
      double x_old_st = Xold[pIndx];
      double v_old_st = Vold[pIndx];
      double x_half, v_half;
      double x, v;

      /////////////////////////////////////////////////////////////////
      //
      // Particle subcycling loop
      //
      for (int st = 0; st < this->data->sub_nt[sIndx]; st++) {
        v = v_old_st;               // full time velocity
        v_half = v_old_st;          // half time velocity
        double x_k = x_old_st;      // Picard iteration position
        double v_k = v_old_st;      // Picard iteration velocity

        /////////////////////////////////////////////////////////////////
        //
        // Picard Crank-Nicholson convergence loop
        // 
        int numIter = 0;
        int maxIter = 10;
        double rel_v = 1.0;
        double rel_x = 1.0;

        while ((fabs(rel_v) >= tol_picard && 
                fabs(rel_x) >= tol_picard) &&
               numIter < maxIter) {
          numIter++;
  
          // Push position based on old position and current velocity
          // Push velocity based on current position and old velocity
          if (this->data->temp_order == ORDER_1_BACK_EULER) {
            xSolve(x_old_st, v, &x, &x_half, dt);
            vSolve(x, v_old_st, E, &v, qm, dt);
          }
          else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
            xSolve(x_old_st, v_half, &x, &x_half, dt);
            vSolve(x_half, v_old_st, E_half, &v, qm, dt);
          }
          v_half = (v + v_old_st) / 2.0;
  
          // Relative difference between last Picard iteration and this one
          rel_x = fabs(x - x_k) / fabs(x);
          rel_v = fabs(v - v_k) / fabs(v);

          // Save position and velocity from this Picard iteration
          x_k = x;
          v_k = v;
        }
        // Picard convergence loop
        /////////////////////////////////////////////////////////////////
  
        // Set the starting value for the next subcycle
        x_old_st = x;
        v_old_st = v;
      }
      // Subcycling loop
      /////////////////////////////////////////////////////////////////

      // Store final particle position and velocity for time and half time
      X[pIndx] = x;
      V[pIndx] = v;

      // Tally half time moment quantities
      half->tallyMomentCalculation(sIndx, x_half, v_half, 1.0);

      // Tally full time moment quantities
      full->tallyMomentCalculation(sIndx, x, v, 1.0);
    }
    // Every particle loop
    /////////////////////////////////////////////////////////////////

    // Scale the half and full moment quantities after accumulating
    full->scaleMoment(sIndx);
    half->scaleMoment(sIndx);

    // Smooth results
    full->smoothMoment(sIndx);
    half->smoothMoment(sIndx);
 
    // Calculate charge conservation properties
    double* r_p_old = oldParticle->getR(sIndx);
    double cc_rms = 0.0;
    double dt_recip = 1.0 / this->data->dt;

    double* r = full->getR(sIndx);
    double* ru = full->getRU(sIndx);
    double* ru_half = half->getRU(sIndx);
// Will Taitano edit 06/25/2012 for debugging purpose
/*      cout << "For specie " << sIndx << endl;
      for (int fIndx = 0; fIndx < this->data->nx+1; fIndx++) {
          cout << "ru_half[" << fIndx << "] = " << ru_half[fIndx] << endl;
      }*/
    if (this->data->temp_order == ORDER_1_BACK_EULER) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        double cc = dt_recip * (r[cIndx] - r_p_old[cIndx]) +
                    dx_recip * (ru[cIndx+1] - ru[cIndx]);
        cc_rms += cc * cc;
      }
    }
    else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        double cc = dt_recip * (r[cIndx] - r_p_old[cIndx]) +
                    dx_recip * (ru_half[cIndx+1] - ru_half[cIndx]);
        cc_rms += cc * cc;
      }
    }

    // Calculate the root mean square of conservation equation
    cc_rms = sqrt(cc_rms / this->data->nx);
    cout << endl <<  "***   CC_RMS " << cc_rms << endl << endl;

    // Store full and half time moment quantities into current particle solution
    curParticle->storeR(sIndx, full->getR(sIndx));
    curParticle->storeS2(sIndx, full->getS2(sIndx));
    curParticle->storeRU(sIndx, full->getRU(sIndx));

    curParticle->storeS2Half(sIndx, half->getS2(sIndx));
    curParticle->storeRUHalf(sIndx, half->getRU(sIndx));

    curParticle->storeCC_RMS(sIndx, cc_rms);
  }
  // Every species loop 
  /////////////////////////////////////////////////////////////////

  delete full;
  delete half;
  delete [] E;
  delete [] E_old;
  delete [] E_half;
}

///////////////////////////////////////////////////////////////////////
//
// Solves for the new particle position via Picard iteration by
// non-linear updates to the velocity
//
///////////////////////////////////////////////////////////////////////

void PushSimple::xSolve(
                        double x_old,          // Input current location
                        double v,              // Input velocity update
                        double* x,             // Output location update
                        double* x_half,        // Output location update
                        double dt)
{
  // Calculate a new location and an update for the velocity push
  *x = x_old + dt * v;
  *x_half = (*x + x_old) / 2.0;

  // Make sure the new locations haven't escaped the system
  xBoundaryCheck(x);
  xBoundaryCheck(x_half);
}

///////////////////////////////////////////////////////////////////////
//
// Reassign a new position if a particle escapes the system
// Updates in place
//
///////////////////////////////////////////////////////////////////////

void PushSimple::xBoundaryCheck(double* x)
{
  if (*x >= this->data->lx) {
    double dx_bound = *x - this->data->lx;
    double num_div = ceil(dx_bound / this->data->lx);
    dx_bound = *x - num_div * this->data->lx;
    *x = dx_bound;
  }
  else if (*x <= 0.0) {
    double dx_bound = fabs(*x);
    double num_div = floor(fabs(dx_bound / this->data->lx));
    dx_bound = fabs(*x) - num_div * this->data->lx;
    *x = this->data->lx - fabs(dx_bound);
  }
}

///////////////////////////////////////////////////////////////////////
//
// Solves for the new particle velocity via Picard iteration by
// non-linear updates to the position and electric field
//
///////////////////////////////////////////////////////////////////////

void PushSimple::vSolve(
                        double x,              // Input location update
                        double v_old,          // Input current velocity
                        double* E,             // Input electric field
                        double* v,             // Output new velocity
                        double qm,
                        double dt)
{
  // Electric interpolation for particle
  // Get the lower face index for this particle
  int fIndx = (int) (x / this->data->dx);

  // Distance to the lower face
  double ddx = x - this->data->xpos_face[fIndx]; 
  double dEdx = (E[fIndx+1] - E[fIndx]) / this->data->dx;
  double E_p = E[fIndx] + dEdx * ddx;

  // New velocity
  *v = v_old + dt * qm * E_p;
}
