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

#include "AmpereVlasov.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "Particles.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"
#include "NLResidual.h"
#include "LOSolverJFNK.h"
#include "LOSolverSI.h"
#include "DhatConsistency.h"
#include "OutGnuplot.h"
#include "stdio.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

AmpereVlasov::AmpereVlasov(string plasmaInputFile)
{
  // Set system parameters
  this->data = new PlasmaData(plasmaInputFile);
  // Create particle solutions
  this->curParticle = new ParticleSolution(this->data);
  this->oldParticle = new ParticleSolution(this->data);

  // Create moment solutions
  this->curMoment = new MomentSolution(this->data);
  this->oldMoment = new MomentSolution(this->data);
  // Create the particle initializer and pusher
  this->particles = new Particles(this->data);

  // Create the low order solver
  if (this->data->si_lo_flag == JFNK_LO_SOLVER) {
    this->loSolver = new LOSolverJFNK(this->data);
  }
  else if (this->data->si_lo_flag == SI_LO_SOLVER) {
    this->loSolver = new LOSolverSI(this->data);
  }

  // Create the consistency calculator
  this->Dhat = new DhatConsistency(this->data);

  // Create the output processor
  this->plasmaOut = new OutGnuplot(this->data);

  // Charge conservation and kinetic energy by species
  this->cc_rms = new double*[this->data->p_size];
  this->KE = new double*[this->data->p_size];

  // Initialization of time and summary arrays
  // Store initial state in time index 0
  this->numberOfTimeSteps = (int) (this->data->tmax / this->data->dt) + 1;
  this->time_vec = new double[this->numberOfTimeSteps];

  this->E_tot = new double[this->numberOfTimeSteps];
  this->mv_tot = new double[this->numberOfTimeSteps];

  this->E_diff = new double[this->numberOfTimeSteps];
  this->mv_diff = new double[this->numberOfTimeSteps];

  this->EE = new double[this->numberOfTimeSteps];
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    this->KE[sIndx] = new double[this->numberOfTimeSteps];
    this->cc_rms[sIndx] = new double[this->numberOfTimeSteps];
  }
}

AmpereVlasov::~AmpereVlasov()
{
  delete this->data;
  delete this->particles;
  delete this->curParticle;
  delete this->oldParticle;
  delete this->curMoment;
  delete this->oldMoment;
  delete this->plasmaOut;

  delete [] this->time_vec;
  delete [] this->E_tot;
  delete [] this->mv_tot;
  delete [] this->E_diff;
  delete [] this->mv_diff;

  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    delete [] this->cc_rms[sIndx];
    delete [] this->KE[sIndx];
  }
  delete [] this->cc_rms;
  delete [] this->KE;
  delete [] this->EE;
}

///////////////////////////////////////////////////////////////////////
//
// Preprocessing
//
///////////////////////////////////////////////////////////////////////
void AmpereVlasov::preprocess()
{
  cout << "AmpereVlasov::preprocess()" << endl;
  cout << setprecision(16) << endl;
    
  // Calculate the cell and face positions and connectivity
  this->data->setCellFaceInformation();

  // Initialize particle positions and velocities (decides number to create)

  this->particles->initialCondition(this->curParticle);
  // Allocate same number of particles in old solution
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    this->oldParticle->allocateParticles(sIndx, 
                            this->curParticle->getNumberOfParticles(sIndx));
    
  // Calculate moments on initial particle positions
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    this->curMoment->momentCalculation(sIndx,
                      this->curParticle->getNumberOfParticles(sIndx),
                      this->curParticle->getX(sIndx),
                      this->curParticle->getV(sIndx));
    //Print moment values
    //this->curMoment->printMoment(sIndx);
  }

  // Poisson solve uses spatial cell density from moment calculation
  this->curMoment->electricFieldCalculation();

  // Smooth moment quantities
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    this->curMoment->smoothMoment(sIndx);

  // Copy parts of initial moment solution (low order)
  // to initial particle solution (high order)
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    this->curParticle->storeR(sIndx, this->curMoment->getR(sIndx));
    this->curParticle->storeRU(sIndx, this->curMoment->getRU(sIndx));
    this->curParticle->storeS2(sIndx, this->curMoment->getS2(sIndx));

    this->curParticle->storeRUHalf(sIndx, this->curMoment->getRU(sIndx));
    this->curParticle->storeS2Half(sIndx, this->curMoment->getS2(sIndx));
  }

  // Calculate average initial current
  this->data->j_avg = 0.0;
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    double* ru = curMoment->getRU(sIndx);
    double mean_ru = 0.0;
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      mean_ru += ru[fIndx];
    mean_ru /= this->data->nfx;
    this->data->j_avg += this->data->q[sIndx] * mean_ru;
  }

  // Copy the current solutions to the old solutions
  this->oldParticle->copyAllParticle(this->curParticle);
  this->oldMoment->copyAllMoment(this->curMoment);

    //Will Taitano edit: 06/25/2012 for debugging purposes
    for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
/*        double* ru_old = this->oldParticle->getRU(sIndx);
        //        for (int fIndx = 0; fIndx this->data->nx; fIndx++) {
        cout << "old time ru at left boundary is: " << ru_old[0] << endl;
        cout << "old time ru at right boundary is: " << ru_old[this->data->nx] << endl;*/
        //        }
    }    
    
  // Calculate total energy of the system
  totalEnergyCalculation(0);

  // Print initial condition
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    cout << "The total number of particles in the system for "
         << sIndx << " th specie is " << this->data->NP0_tot[sIndx] << endl;
    cout << "The plasma wave frequency for the "
        << sIndx << " th specie is " << this->data->w[sIndx] << endl;
    cout << "The plasma wave time-scale for the "
         << sIndx << " th specie is " << this->data->tau[sIndx] << endl;
  }
  cout << "Total initial energy: " << this->E_tot[0] << endl;
  cout << "Total initial momentum: " << this->mv_tot[0] << endl;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate the total energy of the system
// Total energy conservation definition:
//   KE^n+1 + EE^n+1  M = KE^n + EE^n
// Initial value stored in 0, time steps numbered starting with 1
//
///////////////////////////////////////////////////////////////////////

void AmpereVlasov::totalEnergyCalculation(int tIndx)
{
  // Kinetic energy summation
  this->mv_tot[tIndx] = 0.0;

  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    this->KE[sIndx][tIndx] = 0.0;
    double mass = this->data->m[sIndx];
    double mpp = this->data->mpp[sIndx];
    double mp_mass = mpp * mass;

    double* x = this->curParticle->getX(sIndx);
    double* v = this->curParticle->getV(sIndx);
    int numberOfParticles = this->curParticle->getNumberOfParticles(sIndx);

    for (int pIndx = 0; pIndx < numberOfParticles; pIndx++) {
      this->mv_tot[tIndx] += v[pIndx];
      this->KE[sIndx][tIndx] += v[pIndx] * v[pIndx];
    }
    this->mv_tot[tIndx] *= mp_mass;
    this->KE[sIndx][tIndx] *= (0.5 * mp_mass);
  }

  // Electric field energy summation
  double* E = this->curMoment->getE();
  this->EE[tIndx] = 0.0;
  for (int fIndx  = 0; fIndx < this->data->nfx; fIndx++)
    this->EE[tIndx] += E[fIndx] * E[fIndx];
  this->EE[tIndx] *= this->data->dx * 0.5 / this->data->wp / this->data->wp;

  // Total energy
  this->E_tot[tIndx] = this->EE[tIndx];
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    this->E_tot[tIndx] += this->KE[sIndx][tIndx];

  this->E_diff[tIndx] = this->E_tot[0] - this->E_tot[tIndx];
  this->mv_diff[tIndx] = this->mv_tot[0] - this->mv_tot[tIndx];
}

///////////////////////////////////////////////////////////////////////
//
// Simulation loop
//
///////////////////////////////////////////////////////////////////////

void AmpereVlasov::simulation()
{
  // Create non-linear residual containers for outer and inner loops
  NLResidual* residual = new NLResidual(this->data);
  double maxResidual0, maxResidual;
  int ho_trunc = 30;
  double time = 0.0;
  int tIndx = 0;
  this->time_vec[tIndx] = 0.0;
  int picard_outer_tot = 0;             // Total outer Picard iteration count
  int newton_inner_tot = 0;             // Total inner Newton iteration count

  /////////////////////////////////////////////////////////////////////////
  //``
  // Iterate through time
  //
  while (time < this->data->tmax) {
    time += this->data->dt;
    tIndx++;

    // Normalize time to electron plasma wave scale
    this->time_vec[tIndx] = time * this->data->w[0];
    cout << "******************************************************" << endl;
    cout << "The current time = " << time << endl;
    cout << "The current normalized time = " << this->time_vec[tIndx] << endl;
    cout << "******************************************************" << endl;

    // Initialize consistency term calculation for first Newton
    this->Dhat->initialize();

    // Calculate the initial outer loop residual on E and all species r, ru
    maxResidual0 = residual->calculateNLResidual(
                               this->curMoment, this->oldMoment,
                               this->curParticle, this->oldParticle,
                               this->Dhat);
    cout << endl << "The initial outer: " << maxResidual0 << endl;
    residual->printResidual();

    ///////////////////////////////////////////////////////////////////////
    //
    // Outer Picard iteration (high order solver)
    //
    int k_outer = 0;                    // Outer Newton iteration count
    int flag_conv_outer = 0;            // Convergence for outer Newton loop
    while (flag_conv_outer == 0) {
      k_outer++;
// Edit Will Taitano 06/25/2012: For debugging purpose
/*        for (int sIndx = 0 ; sIndx < this->data->p_size; sIndx++) {
            
            double* ru_half = this->curMoment->getRU(sIndx);
            cout << "For " << k_outer << " th outer iteration, before LO solve" << endl;
            cout << "For specie " << sIndx << endl;
            for (int fIndx = 0; fIndx < this->data->nx+1; fIndx++) {
                cout << "ru_half[" << fIndx << "] = " << ru_half[fIndx] << endl;
            }        
        }      */
      // Inner Newton iteration (low order solver)
      int k_inner = this->loSolver->solve(
                               this->curMoment, this->oldMoment,
                               this->curParticle, this->oldParticle,
                               this->Dhat);
      newton_inner_tot += k_inner;

      // Push the particles (updates curParticle solution)
      this->particles->pushParticles(
                               this->curMoment, this->oldMoment,
                               this->curParticle, this->oldParticle);

//Will Taitano edit: 06/25/2012: For debugging purpose only.
/*        for (int sIndx = 0 ; sIndx < this->data->p_size; sIndx++) {
            
            double* ru_half = this->curMoment->getRU(sIndx);
            cout << "For " << k_outer << " th outer iteration, after HO solve" << endl;
            cout << "For specie " << sIndx << endl;
            for (int fIndx = 0; fIndx < this->data->nx+1; fIndx++) {
                cout << "ru_half[" << fIndx << "] = " << ru_half[fIndx] << endl;
            }        
        }      */
      // Calculate the consistency term for both continuity and momentum
      this->Dhat->DhatMake(
                               this->curMoment, this->oldMoment,
                               this->curParticle, this->oldParticle);
        //Will Taitano edit: 06/25/2012: For debugging purpose only.
/*        for (int sIndx = 0 ; sIndx < this->data->p_size; sIndx++) {
            
            double* ru_half = this->curMoment->getRU(sIndx);
            cout << "For " << k_outer << " th outer iteration, after Dhat calc" << endl;
            cout << "For specie " << sIndx << endl;
            for (int fIndx = 0; fIndx < this->data->nx+1; fIndx++) {
                cout << "ru_half[" << fIndx << "] = " << ru_half[fIndx] << endl;
            }        
        }*/
      // Calculation of the new non-linear residual function for convergence
      maxResidual = residual->calculateNLResidual(
                               this->curMoment, this->oldMoment,
                               this->curParticle, this->oldParticle,
                               this->Dhat);
      cout << endl << "k_outer = " << k_outer << " Newton iteration: " 
           << maxResidual << endl;
      residual->printResidual();
// Will Taitano edit: 06/24/2012 accurate commenting. 
      //Check for convergence in the outer HO Picard loop (particle solution)
//      // Check for convergence in the Newton loop (particle solution)
      if (this->data->tol_flag == 0) {
        if (maxResidual <= this->data->tol_nlres * maxResidual0 || 
            k_outer == ho_trunc) {
          flag_conv_outer = 1;
        }
      }
      else if (this->data->tol_flag == 1) {
        if (maxResidual <= this->data->tol_nlres || k_outer == ho_trunc) {
          flag_conv_outer = 1;
        }
      }
    }
    // Close of outer loop
    //
    ///////////////////////////////////////////////////////////////////////
    
    picard_outer_tot += k_outer;

    // Calculate total energy and momentum of system
    totalEnergyCalculation(tIndx);

    // Store charge conservation
    for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
      this->cc_rms[sIndx][tIndx] = this->curParticle->getCC_RMS(sIndx);

    // Output for this time step
    this->plasmaOut->outputPlots(
                       tIndx, this->curMoment, this->curParticle,
                       this->time_vec, this->cc_rms, 
                       this->E_tot, this->E_diff, this->mv_diff);

    // Copy current moment and particle solutions to the old moment
    this->oldMoment->copyParticleSolution(curParticle);
    this->oldMoment->copyMomentSolution(curMoment);
  
    // Swap the particle solution buffers to save copying
    ParticleSolution* ptemp = this->oldParticle;
    this->oldParticle = this->curParticle;
    this->curParticle = ptemp;
  }
  // Close of time loop
  //
  /////////////////////////////////////////////////////////////////////////

  //CPU_t_tot = toc;
  double picard_avg = (double) picard_outer_tot / (double) tIndx;
  double newton_avg = (double) newton_inner_tot / (double) picard_outer_tot;
  cout << "The average Picard iteration was " << picard_avg << endl;
  cout << "The average Newton iteration was " << newton_avg << endl;

  for (int tIndx = 0; tIndx < numberOfTimeSteps; tIndx++) {
    cout << "tIndx: " << tIndx << "  time: " << this->time_vec[tIndx] << endl;
    cout << "    E_tot   " << this->E_tot[tIndx] << endl;
    cout << "    E_diff  " << this->E_diff[tIndx] << endl;
    cout << "    mv_diff " << this->mv_diff[tIndx] << endl;
    for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
      cout << "  cc_rms " << sIndx << " " << this->cc_rms[sIndx][tIndx] << endl;
  }
    
  delete residual;
}
