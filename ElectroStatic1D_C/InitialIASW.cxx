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

#include "InitialIASW.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

InitialIASW::InitialIASW(PlasmaData* data) : ParticleInitialize(data)
{
}

InitialIASW::~InitialIASW()
{
}

///////////////////////////////////////////////////////////////////////
//
// Calculate spatial density, velocity and temperature profile for ELECTRON
//
///////////////////////////////////////////////////////////////////////

void InitialIASW::calculateProfileElectron()
{
  cout << "InitialIASW::calculateProfileElectron()" << endl;
  int sIndx = ELECTRON;
  this->mean_T[sIndx] = 0.0;
  this->mean_u[sIndx] = 0.0;
  this-> mean_rho[sIndx] = 0.0;

  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {

    // Wave factor times cell center in position space
    double value = sin(this->data->k_wave * this->data->xpos_node[cIndx]);

    // Profile is unperturbed plus perturbed times distance factor
    if (this->data->prob_case == 1) {
      this->rho[sIndx][cIndx] = 
                  this->data->rho0[sIndx] + (this->data->alp_r[sIndx] * value);
      this->u[sIndx][cIndx] = 
                  this->data->u0[sIndx] + this->data->alp_u[sIndx] * value;
      this->T[sIndx][cIndx] =
                  this->data->T0[sIndx] + this->data->alp_T[sIndx] * value;
    }
    else if (this->data->prob_case == 2 || this->data->prob_case == 3) {
        this->rho[sIndx][cIndx] = 
        this->data->rho0[sIndx] + (this->data->alp_r[sIndx] * value);
        this->u[sIndx][cIndx] = 
        this->data->u0[sIndx] + this->data->alp_u[sIndx] * value;
        this->T[sIndx][cIndx] =
        this->data->T0[sIndx] + this->data->alp_T[sIndx] * value;
        /*      this->rho[sIndx][cIndx] = 
                  this->data->rho0[sIndx] +
                  (this->data->alp_r[sIndx] * value) - 
                  this->data->T0[sIndx] * this->data->alp_r[sIndx] *
                  this->data->k_wave * this->data->k_wave * value /
                  (this->data->w_scale * this->data->w_scale);
      this->u[sIndx][cIndx] = 
                  this->data->u0[sIndx] + this->data->alp_u[sIndx] * value;
      this->T[sIndx][cIndx] =
                  this->data->T0[sIndx] + this->data->alp_T[sIndx] * value;*/
    }

    this->mean_T[sIndx] += this->T[sIndx][cIndx];
    this->mean_u[sIndx] += this->u[sIndx][cIndx];
    this->mean_rho[sIndx] += this->rho[sIndx][cIndx];
  }
  this->mean_T[sIndx] /= this->data->nx;
  this->mean_u[sIndx] /= this->data->nx;
  this->mean_rho[sIndx] /= this->data->nx;
/*
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "initIASW rho " << rho[sIndx][cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "initIASW u " << u[sIndx][cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "initIASW T " << T[sIndx][cIndx] << endl;
*/

  // Calculate the velocity space parameters for species
  this->data->v_th[sIndx] =
                  sqrt((2.0 * this->mean_T[sIndx]) / this->data->m[sIndx]);
  this->data->vmax[sIndx] =
                  5.0 * this->data->v_th[sIndx] - this->mean_u[sIndx];
  this->data->vmin[sIndx] =
                 -5.0 * this->data->v_th[sIndx] - this->mean_u[sIndx];
  this->data->dv[sIndx] = 
                  (this->data->vmax[sIndx] - this->data->vmin[sIndx]) /
                  (2.0 * this->data->nv + 1);

  // Calculate the velocity space node center
  for (int vIndx = 0; vIndx < this->data->nvc; vIndx++)
    this->data->v_vec[sIndx][vIndx] =
                  this->data->vmin[sIndx] + (vIndx * this->data->dv[sIndx]);
/*
for (int vIndx = 0; vIndx < this->data->nvc; vIndx++)
cout << "v_vec " << data->v_vec[sIndx][vIndx] << endl;
*/

  // Calculate the velocity space face
  for (int vIndx = 0; vIndx < this->data->nvf; vIndx++)
    this->data->v_vec_face[sIndx][vIndx] =
                   this->data->vmin[sIndx] +
                   ((vIndx - 1.5) * this->data->dv[sIndx]);
}

///////////////////////////////////////////////////////////////////////
//
// Calculate spatial density, velocity and temperature profile for ION
//
///////////////////////////////////////////////////////////////////////

void InitialIASW::calculateProfileIon()
{
  cout << "InitialIASW::calculateProfileIon()" << endl;
  int sIndx = ION;
  this->mean_T[sIndx] = 0.0;
  this->mean_u[sIndx] = 0.0;
  this-> mean_rho[sIndx] = 0.0;

  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {

    // Wave factor times cell center in position space
    double value = sin(this->data->k_wave * this->data->xpos_node[cIndx]);

    // Profile is unperturbed plus perturbed times distance factor
    this->rho[sIndx][cIndx] = this->data->rho0[sIndx] +
                              this->data->alp_r[sIndx] * value;
    this->u[sIndx][cIndx] = this->data->u0[sIndx] +
                            this->data->alp_u[sIndx] * value;
    this->T[sIndx][cIndx] = this->data->T0[sIndx] +
                            this->data->alp_T[sIndx] * value;

    this->mean_T[sIndx] += this->T[sIndx][cIndx];
    this->mean_u[sIndx] += this->u[sIndx][cIndx];
    this->mean_rho[sIndx] += this->rho[sIndx][cIndx];
  }
  this->mean_T[sIndx] /= this->data->nx;
  this->mean_u[sIndx] /= this->data->nx;
  this->mean_rho[sIndx] /= this->data->nx;
/*
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "initIASW ION rho " << rho[sIndx][cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "initIASW ION u " << u[sIndx][cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "initIASW ION T " << T[sIndx][cIndx] << endl;
*/

  // Calculate the velocity space parameters for species
  this->data->v_th[sIndx] =
                  sqrt((2.0 * this->mean_T[sIndx]) / this->data->m[sIndx]);
  this->data->vmax[sIndx] = 10.0 * this->data->v_th[sIndx];
  this->data->vmin[sIndx] = -1.0 * this->data->vmax[sIndx];
  this->data->dv[sIndx] =
                  (this->data->vmax[sIndx] - this->data->vmin[sIndx]) /
                  (2.0 * this->data->nv + 1);

  // Calculate the velocity space node center
  for (int vIndx = 0; vIndx < this->data->nvc; vIndx++)
    this->data->v_vec[sIndx][vIndx] =
                  this->data->vmin[sIndx] + (vIndx * this->data->dv[sIndx]);
/*
for (int vIndx = 0; vIndx < this->data->nvc; vIndx++)
cout << "v_vec " << data->v_vec[sIndx][vIndx] << endl;
*/

  // Calculate the velocity space face
  for (int vIndx = 0; vIndx < this->data->nvf; vIndx++)
    this->data->v_vec_face[sIndx][vIndx] =
                   this->data->vmin[sIndx] +
                   ((vIndx - 1.5) * this->data->dv[sIndx]);
}

///////////////////////////////////////////////////////////////////////
//
// Ionic acoustic shock wave initialization
//
///////////////////////////////////////////////////////////////////////

void InitialIASW::initialize(ParticleSolution* curParticle)
{
  cout << "InitialIASW::initialize()" << endl;
  // Calculate positions and velocities for particles
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    if (sIndx == ELECTRON)
      calculateProfileElectron();
    else if (sIndx == ION)
      calculateProfileIon();

    // Random placement of particles
    if (this->data->quiet_start_flag == RANDOM) {
      initializeRandom(sIndx, curParticle);
/*        cout << "The number for e is: " << this->data->NP0_tot[0] << endl;
        cout << "The number for i is: " << this->data->NP0_tot[1] << endl;
      createDistributionFunction(sIndx);
        cout << "This is a random distribution \n" << endl;
        cout << "The quiet flag is: " << this->data->quiet_start_flag << " \n" << endl;
        cout << "The number for e is: " << this->data->NP0_tot[0] << endl;
        cout << "The number for i is: " << this->data->NP0_tot[1] << endl;*/
/*
      if (this->data->wo == 0 || this->data->wo == 1) {
        for (int k = 0; k < 1; k++) {
          initialConditionCleanup(sIndx, curParticle);
          momentMatching(sIndx, curParticle);
        }
      }
      else if (this->data->wo == 2) {
        for (int k = 0; k < 5; k++) {
          initialConditionCleanup(sIndx, curParticle);
          momentMatching(sIndx, curParticle);
        }
      }
*/
    }

    // Deterministic placement of particles
    else {
      initializeDeterministic(sIndx, curParticle);
    }
  }
}
