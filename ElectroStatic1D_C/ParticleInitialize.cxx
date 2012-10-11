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

#include "ParticleInitialize.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"
#include <vector>

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

ParticleInitialize::ParticleInitialize(PlasmaData* data)
{
  this->data = data;

  // Position space arrays
  this->rho = new double*[this->data->p_size];
  this->u = new double*[this->data->p_size];
  this->T = new double*[this->data->p_size];

  for (int sIndx = 0; sIndx < this->data->nx; sIndx++) {
    this->rho[sIndx] = new double[this->data->nx];
    this->u[sIndx] = new double[this->data->nx];
    this->T[sIndx] = new double[this->data->nx];
  }

  this->mean_rho = new double[this->data->p_size];
  this->mean_u = new double[this->data->p_size];
  this->mean_T = new double[this->data->p_size];

  // Distribution function for deterministic placement of particles
  this->fDist = new double*[this->data->nx];
  this->fCount = new int*[this->data->nx];
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    fDist[cIndx] = new double[2 * this->data->nv + 1];
    fCount[cIndx] = new int[2 * this->data->nv + 1];
  }
}

ParticleInitialize::~ParticleInitialize()
{
  for (int sIndx = 0; sIndx < this->data->nx; sIndx++) {
    delete [] this->rho[sIndx];
    delete [] this->u[sIndx];
    delete [] this->T[sIndx];
    delete [] this->fDist[sIndx];
    delete [] this->fCount[sIndx];
  }
  delete [] this->rho;
  delete [] this->u;
  delete [] this->T;
  delete [] this->mean_rho;
  delete [] this->mean_u;
  delete [] this->mean_T;
  delete [] this->fDist;
  delete [] this->fCount;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate spatial density, velocity and temperature profile
//
///////////////////////////////////////////////////////////////////////

void ParticleInitialize::calculateProfile()
{
  cout << "ParticleInitialize::calculateProfile()" << endl;
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    this->mean_T[sIndx] = 0.0;
    this->mean_u[sIndx] = 0.0;
    this->mean_rho[sIndx] = 0.0;

    // Calculate spatial density, velocity, temperature profiles
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {

      // Wave factor times cell center in position space
      double value = cos(this->data->k_wave * this->data->xpos_node[cIndx]);

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
cout << "profile rho " << rho[sIndx][cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "profile u " << u[sIndx][cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "profile T " << T[sIndx][cIndx] << endl;
*/

    // Calculate the velocity space parameters for species
    this->data->v_th[sIndx] = sqrt((2.0 * this->mean_T[sIndx]) / 
                              this->data->m[sIndx]);
    this->data->vmax[sIndx] = 5.0 * this->data->v_th[sIndx];
    this->data->vmin[sIndx] = -1.0 * this->data->vmax[sIndx];
    this->data->dv[sIndx] = this->data->vmax[sIndx] / this->data->nv;

    // Calculate the velocity space node center
    for (int vIndx = 0; vIndx < this->data->nvc; vIndx++) {
      this->data->v_vec[sIndx][vIndx] = this->data->vmin[sIndx] +
                                        (vIndx * this->data->dv[sIndx]);
    }

    // Calculate the velocity space face
    for (int vIndx = 0; vIndx < this->data->nvf; vIndx++) {
      this->data->v_vec_face[sIndx][vIndx] = this->data->vmin[sIndx] +
                                       ((vIndx - 0.5) * this->data->dv[sIndx]);
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
// Particle initialization random placement of particles
// For Landau all species and two-stream non-electron
//
///////////////////////////////////////////////////////////////////////

void ParticleInitialize::initializeRandom(
                              int sIndx,
                              ParticleSolution* curParticle)
{
  double dxd2 = this->data->dx / 2.0;
  int* rho_NP = new int[this->data->nx];

  // Calculate number of particles per cell
  this->data->NP0_tot[sIndx] = 0;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    rho_NP[cIndx] = (int) round(this->data->NP_ref[sIndx] * 
                                this->rho[sIndx][cIndx] * this->data->dx);
    this->data->NP0_tot[sIndx] += rho_NP[cIndx];
  }
  this->data->mean_rho_NP[sIndx] = (double) this->data->NP0_tot[sIndx] / 
                                   (double) this->data->nx;

  // Allocate space for particles inside the solution and get pointers
  // Write directly into the arrays for efficiency
  curParticle->allocateParticles(sIndx, this->data->NP0_tot[sIndx]);
  double* X = curParticle->getX(sIndx);
  double* V = curParticle->getV(sIndx);

  // Initialize particle position and velocity
  int pIndx = 0;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    // Separation between particles in the cell
    double ddx = this->data->dx / (double) rho_NP[cIndx];

    for (int p = 0; p < rho_NP[cIndx]; p++) {
      X[pIndx] = this->data->xpos_node[cIndx] - dxd2 + ddx * (p + 0.5);
      V[pIndx] = this->u[sIndx][cIndx] + 
                 sqrt(-2.0 * this->T[sIndx][cIndx] * 
                   log(genRandom()) / this->data->m[sIndx]) *
                 cos(2.0 * M_PI * genRandom());
      pIndx++;
    }
  }
  // Mass per particle per species
  this->data->mpp[sIndx] = (this->mean_rho[sIndx] * this->data->dx) / 
                           this->data->mean_rho_NP[sIndx];
  delete [] rho_NP;
}

///////////////////////////////////////////////////////////////////////
//
// Particle initialization deterministic placement of particles
// For Landau all species and two-stream non-electron
//
///////////////////////////////////////////////////////////////////////

void ParticleInitialize::initializeDeterministic(
                              int sIndx,
                              ParticleSolution* curParticle)
{
  cout << "ParticleInitialize::initializeDeterministic()" << endl;
  // Initialization of distribution function and velocity array
  createDistributionFunction(sIndx);

  // Allocate space for particles inside the solution and get pointers
  // Write directly into the arrays for efficiency
  curParticle->allocateParticles(sIndx, this->data->NP0_tot[sIndx]);
  double* X = curParticle->getX(sIndx);
  double* V = curParticle->getV(sIndx);

  // Initialize particle position and velocity vector
  int pIndx = 0;
  double dxd2 = this->data->dx / 2.0;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 0; vIndx < 2*this->data->nv+1; vIndx++) {

      // Number of particles in the cell (cIndx,vIndx)
      int size_v_part = this->fCount[cIndx][vIndx];

      if (size_v_part != 0) {
        // Separation between particles in the cell
        double ddx = (double) this->data->dx / (double) size_v_part;

        for (int p = 0; p < size_v_part; p++) {
          X[pIndx] = this->data->xpos_node[cIndx] - dxd2 + ddx * (p + 0.5);
          V[pIndx] = this->data->v_vec[sIndx][vIndx];
          pIndx++;
        }
      }
    }
  }
  // Mass per particle per species
  this->data->mpp[sIndx] = (this->mean_rho[sIndx] * this->data->dx) /
                            this->data->mean_rho_NP[sIndx];

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
}

///////////////////////////////////////////////////////////////////////
//
// Create distribution function from determinism
//
///////////////////////////////////////////////////////////////////////

void ParticleInitialize::createDistributionFunction(int sIndx)
{
  cout << "ParticleInitialize::createDistributionFunction()" << endl;
  // Initialization of distribution function and velocity array
  this->data->NP0_tot[sIndx] = 0;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 0; vIndx < 2*this->data->nv+1; vIndx++) {
      this->fDist[cIndx][vIndx] =
               sqrt(this->data->m[sIndx] /
                    (2.0 * M_PI * this->T[sIndx][cIndx])) *
               exp(-1.0 * this->data->m[sIndx] *
                   ((this->data->v_vec[sIndx][vIndx] - this->u[sIndx][cIndx]) *
                    (this->data->v_vec[sIndx][vIndx] - this->u[sIndx][cIndx])) /
                   (2.0 * this->T[sIndx][cIndx])) *
               this->rho[sIndx][cIndx];
/*
if (sIndx == 1) {
cout << "T " << T[sIndx][cIndx] << endl;
cout << "Fdist " << fDist[cIndx][vIndx] << endl;
cout << "sqrt part " << sqrt(this->data->m[sIndx] /
                    (2.0 * M_PI * this->T[sIndx][cIndx])) << endl;
cout << "exp innards " << (-1.0 * this->data->m[sIndx] *
                   ((this->data->v_vec[sIndx][vIndx] - this->u[sIndx][cIndx]) *
                    (this->data->v_vec[sIndx][vIndx] - this->u[sIndx][cIndx])) /
                   (2.0 * this->T[sIndx][cIndx])) << endl;
cout << "square part " << ((this->data->v_vec[sIndx][vIndx] - this->u[sIndx][cIndx]) * (this->data->v_vec[sIndx][vIndx] - this->u[sIndx][cIndx])) << endl;
cout << "square part w divide " << ((this->data->v_vec[sIndx][vIndx] - this->u[sIndx][cIndx]) * (this->data->v_vec[sIndx][vIndx] - this->u[sIndx][cIndx])) << endl;
cout << "v_vec " << data->v_vec[sIndx][vIndx] << endl;
cout << "u " << u[sIndx][cIndx] << endl;
}
*/
      this->fCount[cIndx][vIndx] = (int)
                   round(this->data->NP_ref[sIndx] * fDist[cIndx][vIndx] *
                         this->data->dx * this->data->dv[sIndx]);
      this->data->NP0_tot[sIndx] += fCount[cIndx][vIndx];
    }
  }

/*
  cout << "fCount 1 through 15" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 0; vIndx < 15; vIndx++)
      cout << "   " << fCount[cIndx][vIndx];
    cout << endl;
  }
  cout << endl;
  cout << "fCount 16 through 30" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 15; vIndx < 30; vIndx++)
      cout << "   " << fCount[cIndx][vIndx];
    cout << endl;
  }
  cout << endl;
  cout << "fCount 31 through 45" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 30; vIndx < 45; vIndx++)
      cout << "   " << fCount[cIndx][vIndx];
    cout << endl;
  }
  cout << endl;
  cout << "fCount 46 through 60" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 45; vIndx < 60; vIndx++)
      cout << "   " << fCount[cIndx][vIndx];
    cout << endl;
  }
  cout << endl;
  cout << "fCount 61 through 75" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 60; vIndx < 75; vIndx++)
      cout << "   " << fCount[cIndx][vIndx];
    cout << endl;
  }
  cout << endl;
  cout << "fCount 76 through 90" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 75; vIndx < 90; vIndx++)
      cout << "   " << fCount[cIndx][vIndx];
    cout << endl;
  }
  cout << endl;
  cout << "fCount 91 through 105" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 90; vIndx < 105; vIndx++)
      cout << "   " << fCount[cIndx][vIndx];
    cout << endl;
  }
  cout << endl;
  cout << "fCount 106 through 121" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    for (int vIndx = 105; vIndx < 121; vIndx++)
      cout << "   " << fCount[cIndx][vIndx];
    cout << endl;
  }
  cout << endl;
cout << "TOTAL COUNT " << data->NP0_tot[sIndx] << endl;
*/

  this->data->mean_rho_NP[sIndx] = (double) this->data->NP0_tot[sIndx] /
                                   (double) this->data->nx;
}

///////////////////////////////////////////////////////////////////////
//
// Clean up the distribution function for quiet start
//
///////////////////////////////////////////////////////////////////////

void ParticleInitialize::initialConditionCleanup(
                              int sIndx,
                              ParticleSolution* curParticle)
{
  cout << "ParticleInitialize::initialConditionCleanup()" << endl;
  int numberOfParticles = this->data->NP0_tot[sIndx];
  double* X = curParticle->getX(sIndx);
  double* V = curParticle->getV(sIndx);

  MomentSolution* initMoment = new MomentSolution(this->data);
  initMoment->momentCalculation(sIndx, numberOfParticles, X, V);

  // Calculate the number of particles to be added or subtracted from cell
  double* rho_p = initMoment->getR(sIndx);
  int* rho_diff = new int[this->data->nx];
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
/* PKF possible bug in matlab */
    rho_diff[cIndx] = (int) round(this->data->dx * 
                            (this->rho[sIndx][cIndx] - rho_p[cIndx]) /
                            this->data->mpp[sIndx]);
/*
cout << "mpp " << data->mpp[sIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "init_cond rho " << rho[sIndx][cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "init_cond rho_p " << rho_p[cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "init_cond rho_diff " << rho_diff[cIndx] << endl;
*/

  // Number of new particles that are added or subtracted
  int numberToAdd = 0;
  int numberToSub = 0;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    if (rho_diff[cIndx] > 0)
      numberToAdd += rho_diff[cIndx];
    else if (rho_diff[cIndx] < 0)
      numberToSub += rho_diff[cIndx];
  }

  // Allocate array of arrays of particles to add
  double* x = new double[numberToAdd];
  double* v = new double[numberToAdd];

  // Adjust the number of particles per cell by randomly choosing a
  // particle to delete or add.  When adding use Maxwellian
  int tal_as = 0;
  int p_add = 0;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {

    // Number of particles in cell
    int num_cell_i = 0;
    for (int vIndx = 0; vIndx < 2*this->data->nv + 1; vIndx++)
      num_cell_i += this->fCount[cIndx][vIndx];

    // Delete particles from cell
    if (rho_diff[cIndx] < 0) {
      int numberToDelete = -rho_diff[cIndx];
      vector<int> rand_vec;
      for (int i = 0; i < numberToDelete; i++)
        rand_vec.push_back(i);

      // Randomly permute 
      random_shuffle(rand_vec.begin(), rand_vec.end());

      // Delete the random particles
      for (int i = 0; i < numberToDelete; i++) {
        X[tal_as + rand_vec[i]] = -1;
      }
    }
    // Add particles to cell
    else if (rho_diff[cIndx] > 0) {
      for (int i = 0; i < rho_diff[cIndx]; i++) {
        x[p_add] = this->data->xpos_node[cIndx] + 
                   this->data->dx * (genRandom() - 0.5);
        v[p_add] = this->u[sIndx][cIndx] + 
                   sqrt(-2.0 * this->T[sIndx][cIndx] * log(genRandom())) *
                   cos(2.0 * M_PI * genRandom());
        p_add++;
      }
    }
    tal_as += num_cell_i;
  }

  // Update the x,v with added and subtracted particles
  this->data->NP0_tot[sIndx] = numberOfParticles + numberToAdd + numberToSub;
cout << "numberToAdd " << numberToAdd << endl;
cout << "numberToSub " << numberToSub << endl;
cout << "Updated number of particles " << data->NP0_tot[sIndx] << endl;
  double* X_new = new double[this->data->NP0_tot[sIndx]];
  double* V_new = new double[this->data->NP0_tot[sIndx]];
  int pIndxNew = 0;

  // Copy original particles with subtracted particles skipped
  for (int pIndx = 0; pIndx < numberOfParticles; pIndx++) {
    if (X[pIndx] != -1) {
      X_new[pIndxNew] = X[pIndx];
      V_new[pIndxNew] = V[pIndx];
      pIndxNew++;
    }
  }

  // Add in new particles
  for (int pIndx = 0; pIndx < numberToAdd; pIndx++) {
    X_new[pIndxNew] = x[pIndx];
    V_new[pIndxNew] = v[pIndx];
    pIndxNew++;
  }

  // Delete old particles and return new particles
  curParticle->freeParticles(sIndx);
  curParticle->allocateParticles(sIndx, this->data->NP0_tot[sIndx]);
  curParticle->copyParticles(sIndx, X_new, V_new);

  delete [] rho_diff;
  delete [] v;
  delete [] x;
  delete [] X_new;
  delete [] V_new;
}

///////////////////////////////////////////////////////////////////////
//
// Match the fluid velocity moment
// Velocity array is altered in place
//
///////////////////////////////////////////////////////////////////////

void ParticleInitialize::momentMatching(
                              int sIndx,
                              ParticleSolution* curParticle)
{
  cout << "ParticleInitialize::momentMatching()" << endl;
  int numberOfParticles = this->data->NP0_tot[sIndx];
  double* X = curParticle->getX(sIndx);
  double* V = curParticle->getV(sIndx);

  MomentSolution* initMoment = new MomentSolution(this->data);
  initMoment->momentCalculation(sIndx, numberOfParticles, X, V);

  double* u_p = initMoment->getU(sIndx);
  double* T_p = initMoment->getT(sIndx);

  double* du = new double[this->data->nx];
  double* alp = new double[this->data->nx];

  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    du[cIndx] = u_p[cIndx] - this->u[sIndx][cIndx];
    alp[cIndx] = sqrt(this->T[sIndx][cIndx] / T_p[cIndx]);
  }
  for (int pIndx = 0; pIndx < numberOfParticles; pIndx++) {
    int cl = (int) (X[pIndx] / this->data->dx);
    if (cl == this->data->nx)
      cl = this->data->nx - 1;
    V[pIndx] = alp[cl] * (V[pIndx] - du[cl]);
  }

  delete initMoment;
  delete [] du;
  delete [] alp;
}
