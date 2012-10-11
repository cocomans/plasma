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

#include "DhatConsistency.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

DhatConsistency::DhatConsistency(PlasmaData* data)
{
  this->data = data;

  this->Dhat_r = new double*[this->data->p_size];
  this->Dhat_ru = new double*[this->data->p_size];
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    this->Dhat_r[sIndx] = new double[this->data->nx];
    this->Dhat_ru[sIndx] = new double[this->data->nfx]; 
  } 

  initialize();
}

DhatConsistency::~DhatConsistency()
{
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    delete [] this->Dhat_r[sIndx];
    delete [] this->Dhat_ru[sIndx];
  }
  delete [] this->Dhat_r;
  delete [] this->Dhat_ru;
}

///////////////////////////////////////////////////////////////////////
//
// Initialize consistency terms
//
///////////////////////////////////////////////////////////////////////

void DhatConsistency::initialize()
{
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      this->Dhat_r[sIndx][cIndx] = 0.0;
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      this->Dhat_ru[sIndx][fIndx] = 0.0;
  }
}

///////////////////////////////////////////////////////////////////////
//
// Create the consistency enforcing term for both momentum and continuity
// equations using the higher order particle solution
//
///////////////////////////////////////////////////////////////////////

void DhatConsistency::DhatMake(
                     MomentSolution* curMoment,
                     MomentSolution* oldMoment,
                     ParticleSolution* curParticle,
                     ParticleSolution* oldParticle)

{
  double* r_p_half = new double[this->data->nx];
  double* r_p_f_half = new double[this->data->nx];
  double* s2_pb_half = new double[this->data->nx];
  double* s2_pb = new double[this->data->nx];
  double* r_p_f = new double[this->data->nx];
  double* E_half = new double[this->data->nfx];

  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    // Electric field at n+1 and n time step
    double* E = curMoment->getE();
    double* E_old = oldMoment->getE();

    // Particle quantities at time n+1 and n+1/2
    double* r_p = curParticle->getR(sIndx);
    double* ru_p = curParticle->getRU(sIndx);
    double* ru_p_half = curParticle->getRUHalf(sIndx);
    double* s2_p = curParticle->getS2(sIndx);
    double* s2_p_half = curParticle->getS2Half(sIndx);

    // Particle quantities at time n
    double* r_p_old = oldParticle->getR(sIndx);
    double* ru_p_old = oldParticle->getRU(sIndx);

    // Construct density normalized stress tensor
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
      s2_pb_half[cIndx] = s2_p_half[cIndx] / r_p[cIndx];
      s2_pb[cIndx] = s2_p[cIndx] / r_p[cIndx];
    }

    // Construct half time face quantities
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      r_p_half[cIndx] = (r_p[cIndx] + r_p_old[cIndx]) / 2.0;
    makeFaceFlux(r_p_half, r_p_f_half, this->data->nx);

    // Construct full time face quantities
    makeFaceFlux(r_p, r_p_f, this->data->nx);

    // Construct half time electric field
    for (int fIndx = 0; fIndx < (this->data->nfx); fIndx++)
      E_half[fIndx] = (E[fIndx] + E_old[fIndx]) / 2.0;

    // Construct consistency terms
    if (this->data->temp_order == ORDER_1_BACK_EULER) {
      consistencyTerms(sIndx, 1.0, ru_p, s2_pb, r_p_f, E,
                       r_p, r_p_old, ru_p_old);
    }
    else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
      consistencyTerms(sIndx, 2.0, ru_p_half, s2_pb_half, r_p_f_half, E_half,
                       r_p, r_p_old, ru_p_old);
    }
/*
cout << endl << endl;
for (int fIndx = 0; fIndx < (this->data->nfx); fIndx++)
cout << "Dhat E " << E[fIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "Dhat r_p " << r_p[cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "Dhat r_p_old " << r_p_old[cIndx] << endl;
for (int fIndx = 0; fIndx < (this->data->nfx); fIndx++)
cout << "Dhat ru_p " << ru_p[fIndx] << endl;
for (int fIndx = 0; fIndx < (this->data->nfx); fIndx++)
cout << "Dhat ru_p_old " << ru_p_old[fIndx] << endl;
for (int fIndx = 0; fIndx < (this->data->nfx); fIndx++)
cout << "Dhat r_p_f " << r_p_f[fIndx] << endl;
for (int cIndx = 0; cIndx < (this->data->nx); cIndx++)
cout << "Dhat s2_p " << s2_p[cIndx] << endl;
for (int cIndx = 0; cIndx < (this->data->nx); cIndx++)
cout << "Dhat s2_pb " << s2_pb[cIndx] << endl;
for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
cout << "Dhat_r " << Dhat_r[sIndx][cIndx] << endl;
for (int fIndx = 0; fIndx < (this->data->nfx); fIndx++)
cout << "Dhat_ru " << Dhat_ru[sIndx][fIndx] << endl;
*/
  }

  delete[] r_p_half;
  delete[] r_p_f;
  delete[] r_p_f_half;
  delete[] s2_pb;
  delete[] s2_pb_half;
  delete[] E_half;
}

///////////////////////////////////////////////////////////////////////
//
// Construct consistency terms
//
///////////////////////////////////////////////////////////////////////

void DhatConsistency::consistencyTerms(
                        int sIndx,          // input species
                        double factor,      // solver order dependent input
                        double* ru_p,       // solver order dependent input
                        double* s2_pb,      // solver order dependent input
                        double* r_p_f,      // solver order dependent input
                        double* E,          // solver order dependent input
                        double* r_p,
                        double* r_p_old,
                        double* ru_p_old)
{
  double q = this->data->q[sIndx];
  double m = this->data->m[sIndx];
  double dtdx = this->data->dt / this->data->dx;
  double mdx_recip = m / this->data->dx;
  double mdt_recip = m / this->data->dt;

  // Construct consistency term for continuity equation
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    this->Dhat_r[sIndx][cIndx] =
                    ((mdt_recip * (r_p[cIndx] - r_p_old[cIndx])) +
                     (mdx_recip * (ru_p[cIndx+1] - ru_p[cIndx]))) /
                    r_p[cIndx];
  }
  // Construct consistency term for continuity equation
  int last = this->data->nx - 1;
//Will Taitano edit: 06/25/2012 for debugging purpose
/*    cout << "For specie " << sIndx << ":" << endl;*/
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    // Fixed the consistency term calculation
    if (fIndx == 0 || fIndx == this->data->nx)
    {
      this->Dhat_ru[sIndx][fIndx] =
                      (factor * mdt_recip * (ru_p[fIndx] - ru_p_old[fIndx]) +
                       (mdx_recip * 
                        (r_p[0] * s2_pb[0] - r_p[last] * s2_pb[last])) -
                       q * r_p_f[fIndx] * E[fIndx]) / r_p_f[fIndx];
/*      cout << "Dhat_ru[" << fIndx << "] = " << Dhat_ru[sIndx][fIndx] << endl;
        cout << "ru_half[" << fIndx << "] = " << ru_p[fIndx] << endl;
        cout << "ru_old[" << fIndx << "] = " << ru_p_old[fIndx] << endl;*/
    }
    else
    {
      this->Dhat_ru[sIndx][fIndx] =
                      (factor * mdt_recip * (ru_p[fIndx] - ru_p_old[fIndx]) +
                       (mdx_recip * 
                         (r_p[fIndx] * s2_pb[fIndx] - 
                          r_p[fIndx - 1] * s2_pb[fIndx - 1])) -
                       q * r_p_f[fIndx] * E[fIndx]) / r_p_f[fIndx];
/*      cout << "Dhat_ru[" << fIndx << "] = " << Dhat_ru[sIndx][fIndx] << endl;*/
    }
  }
}
