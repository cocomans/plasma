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

#include "InitialTSI.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

InitialTSI::InitialTSI(PlasmaData* data) : ParticleInitialize(data)
{
}

InitialTSI::~InitialTSI()
{
}

///////////////////////////////////////////////////////////////////////
//
// Two stream instability particle position
//
///////////////////////////////////////////////////////////////////////

void InitialTSI::initialize(ParticleSolution* curParticle)
{
  // Calculate the spatial density, velocity and temperature by species
  calculateProfile();

  // Calculate positions and velocities for particles
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    // Place ELECTRON particles using two stream delta function
    if (sIndx == ELECTRON) {
      initializeElectron(curParticle);
    }
    else {
      // Random placement of particles
      if (this->data->quiet_start_flag == RANDOM) {
        initializeRandom(sIndx, curParticle);
      }
      // Deterministic placement of particles
      else {
        initializeDeterministic(sIndx, curParticle);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
// Electron two stream instability particle position
// Based on the spatial profile, place ELECTRON particles
// in position space and velocity space using delta function placement
//
///////////////////////////////////////////////////////////////////////

void InitialTSI::initializeElectron(ParticleSolution* curParticle)
{
  cout << "InitialTSI::initializeElectron()" << endl;
  int* rho_NP = new int[this->data->nx];
  double* ddx = new double[this->data->nx];
  this->data->NP0_tot[ELECTRON] = 0;
  double dxd2 = this->data->dx / 2.0;

  // Calculate number of particles per cell and separation between particles
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    rho_NP[cIndx] = (int) round(this->data->NP_ref[ELECTRON] * 
                                this->rho[ELECTRON][cIndx] * this->data->dx);
    ddx[cIndx] = (double) this->data->dx / (double) rho_NP[cIndx];
    this->data->NP0_tot[ELECTRON] += rho_NP[cIndx];
  }
  this->data->mean_rho_NP[ELECTRON] = (double) this->data->NP0_tot[ELECTRON] / 
                                      (double) this->data->nx;

  // Double the number of particles for two beams
  this->data->NP0_tot[ELECTRON] *= 2;

  // Allocate space for particles inside the solution and get pointers
  // Write directly into the arrays for efficiency
  curParticle->allocateParticles(ELECTRON, this->data->NP0_tot[ELECTRON]);
  double* X = curParticle->getX(ELECTRON);
  double* V = curParticle->getV(ELECTRON);

  // Particles with positivite velocity
  int pIndx = 0;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int p = 0; p < rho_NP[cIndx]; p++) {
      X[pIndx] = this->data->xpos_node[cIndx] - dxd2 + ddx[cIndx] * (p + 0.5);
      V[pIndx] = this->data->u0[ELECTRON]; 
      pIndx++;
    }
  }

  // Particles with negative velocity
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int p = 0; p < rho_NP[cIndx]; p++) {
      X[pIndx] = this->data->xpos_node[cIndx] - dxd2 + ddx[cIndx] * (p + 0.5);
      V[pIndx] = -(this->data->u0[ELECTRON]); 
      pIndx++;
    }
  }

  // Mass per particle per species
  // Division by 2 comes from the symmetry of beams
  this->data->mpp[ELECTRON] = (this->mean_rho[ELECTRON] * this->data->dx) / 
                              this->data->mean_rho_NP[ELECTRON] / 2.0;
  delete [] rho_NP;
  delete [] ddx;
}
