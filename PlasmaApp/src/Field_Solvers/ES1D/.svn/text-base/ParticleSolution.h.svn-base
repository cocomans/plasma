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

#ifndef ParticleSolution_h
#define ParticleSolution_h

using namespace std;

class PlasmaData2;

class ParticleSolution {
public:
  ParticleSolution(PlasmaData2* data);   // Initialization parameters
  ~ParticleSolution();

  // Manage space for particles
  void allocateParticles(int sIndx, int size);
  void freeParticles(int sIndx);
  void copyParticles(int sIndx, double* X, double* V);

  // Complete copy for initialization
  void copyAllParticle(ParticleSolution* p);

  // Store updated arrays in the solution
  void storeR(int sIndx, double* r);
  void storeRU(int sIndx, double* ru);
  void storeS2(int sIndx, double* s2);

  void storeRUHalf(int sIndx, double* ru_half);
  void storeS2Half(int sIndx, double* s2_half);

  void storeCC_RMS(int sIndx, double ccrms)  { this->cc_rms[sIndx] = ccrms; }

  // Retrieve particle solution
  int getNumberOfParticles(int sIndx) { return this->numberOfParticles[sIndx]; }
  double* getX(int sIndx)             { return this->X[sIndx]; }
  double* getV(int sIndx)             { return this->V[sIndx]; }

  double* getR(int sIndx)             { return this->R[sIndx]; }
  double* getRU(int sIndx)            { return this->RU[sIndx]; }
  double* getS2(int sIndx)            { return this->S2[sIndx]; }
  double* getP(int sIndx)             { return this->P[sIndx]; }

  double* getRUHalf(int sIndx)        { return this->RU_half[sIndx]; }
  double* getS2Half(int sIndx)        { return this->S2_half[sIndx]; }

  double  getCC_RMS(int sIndx)        { return this->cc_rms[sIndx]; }

  void ReduceSolution();

private:
  PlasmaData2* data;            // Input parameters

  // Full moment quantities
  double** R;                  // Cell density
  double** RU;                 // Face first moment
  double** S2;                 // Cell second moment
  double** P;                  // Pressure 

  // Half moment quantities
  double** RU_half;            // Face first moment
  double** S2_half;            // Cell second moment

  // Particles
  int* numberOfParticles;      // Number of particles per species for this thread
  int* numberOfParticlesTot; // Total number of particles for all threads
  double** X;                  // Particles in position space by species
  double** V;                  // Particles in velocity space by species

  double* cc_rms;              // Charge conservation by species
};

#endif
