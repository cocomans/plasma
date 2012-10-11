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

#include "ParticleSolution.h"
#include "PlasmaData.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

ParticleSolution::ParticleSolution(PlasmaData* data)
{
  this->data = data;

  this->numberOfParticles = new int[this->data->p_size];
  this->X = new double*[this->data->p_size];
  this->V = new double*[this->data->p_size];

  this->R = new double*[this->data->p_size];
  this->RU = new double*[this->data->p_size];
  this->RU_half = new double*[this->data->p_size];
  this->S2 = new double*[this->data->p_size];
  this->S2_half = new double*[this->data->p_size];
  this->P = new double*[this->data->p_size];

  this->cc_rms = new double[this->data->p_size];

  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    this->R[sIndx] = new double[this->data->nx];
    this->S2[sIndx] = new double[this->data->nx];
    this->S2_half[sIndx] = new double[this->data->nx];
    this->P[sIndx] = new double[this->data->nx];

    this->RU[sIndx] = new double[this->data->nfx];
    this->RU_half[sIndx] = new double[this->data->nfx];
  }
}
/////////////////////////////////////////////////////////////////////////
// 
// Is this the Desctructor? (Will Taitano 06/21/2012 12:55 PM)
//
/////////////////////////////////////////////////////////////////////////
ParticleSolution::~ParticleSolution()
{
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    delete [] this->X[sIndx];
    delete [] this->V[sIndx];
  }
  delete [] this->X;
  delete [] this->V;

  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    delete [] this->R[sIndx];
    delete [] this->RU[sIndx];
    delete [] this->RU_half[sIndx];
    delete [] this->S2[sIndx];
    delete [] this->S2_half[sIndx];
    delete [] this->P[sIndx];
  }
  delete [] this->R;
  delete [] this->RU;
  delete [] this->RU_half;
  delete [] this->S2;
  delete [] this->S2_half;
  delete [] this->P;
  delete [] this->cc_rms;
}

///////////////////////////////////////////////////////////////////////
//
// Allocate space for the particles
//
///////////////////////////////////////////////////////////////////////

void ParticleSolution::allocateParticles(int sIndx, int size)
{
  this->numberOfParticles[sIndx] = size;
  this->X[sIndx] = new double[size];
  this->V[sIndx] = new double[size];
}

///////////////////////////////////////////////////////////////////////
//
// Free space for the particles
//
///////////////////////////////////////////////////////////////////////

void ParticleSolution::freeParticles(int sIndx)
{
  this->numberOfParticles[sIndx] = 0;
  delete [] this->X[sIndx];
  delete [] this->V[sIndx];
}

///////////////////////////////////////////////////////////////////////
//
// Copy particles
//
///////////////////////////////////////////////////////////////////////

void ParticleSolution::copyParticles(int sIndx, double* x, double* v)
{
  this->X[sIndx] = new double[this->numberOfParticles[sIndx]];
  this->V[sIndx] = new double[this->numberOfParticles[sIndx]];
  for (int pIndx = 0; pIndx < this->numberOfParticles[sIndx]; pIndx++) {
    this->X[sIndx][pIndx] = x[pIndx];
    this->V[sIndx][pIndx] = v[pIndx];
  }
}

///////////////////////////////////////////////////////////////////////
//
// Copy Particles
// Needed for first time step and so must size the particles
//
///////////////////////////////////////////////////////////////////////

void ParticleSolution::copyAllParticle(ParticleSolution* p)
{
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    // Copy other parts of solution
    double* copyrho = p->getR(sIndx);
    double* copyru = p->getRU(sIndx);
    double* copyruhalf = p->getRUHalf(sIndx);
    double* copys2 = p->getS2(sIndx);
    double* copys2half = p->getS2Half(sIndx);
    double* copyP = p->getP(sIndx);

    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
      this->R[sIndx][cIndx] = copyrho[cIndx];
      this->S2[sIndx][cIndx] = copys2[cIndx];
      this->S2_half[sIndx][cIndx] = copys2half[cIndx];
      this->P[sIndx][cIndx] = copyP[cIndx];
    }
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
      this->RU[sIndx][fIndx] = copyru[fIndx];
      this->RU_half[sIndx][fIndx] = copyruhalf[fIndx];
    }

    double* copyX = p->getX(sIndx);
    double* copyV = p->getV(sIndx);

    for (int pIndx = 0; pIndx < this->numberOfParticles[sIndx]; pIndx++) {
      this->X[sIndx][pIndx] = copyX[pIndx];
      this->V[sIndx][pIndx] = copyV[pIndx];
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
// Store ParticleSolution arrays
//
///////////////////////////////////////////////////////////////////////

void ParticleSolution::storeR(int sIndx, double* r)
{
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    this->R[sIndx][cIndx] = r[cIndx];
}

void ParticleSolution::storeRU(int sIndx, double* ru)
{
// Will Taitano edit: 06/25/2012
// The maximum index is suppose to be cell face array size
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
        this->RU[sIndx][fIndx] = ru[fIndx];
//  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
//    this->RU[sIndx][cIndx] = ru[cIndx];
}

void ParticleSolution::storeS2(int sIndx, double* s2)
{
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    this->S2[sIndx][cIndx] = s2[cIndx];
}

void ParticleSolution::storeRUHalf(int sIndx, double* ru_half)
{
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    this->RU_half[sIndx][fIndx] = ru_half[fIndx];
}

void ParticleSolution::storeS2Half(int sIndx, double* s2_half)
{
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    this->S2_half[sIndx][fIndx] = s2_half[fIndx];
}
