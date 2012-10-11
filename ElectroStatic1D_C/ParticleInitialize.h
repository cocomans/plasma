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

#ifndef ParticleInitialize_h
#define ParticleInitialize_h

using namespace std;

#include "PlasmaData.h"

class ParticleSolution;

class ParticleInitialize {
public:
  ParticleInitialize(PlasmaData* data);
  ~ParticleInitialize();

  friend class InitialLandau;
  friend class InitialTSI;
  friend class InitialIASW;

  virtual void initialize(ParticleSolution* curParticle) = 0;

  // Calculate the density, velocity and temperature profiles
  void calculateProfile();

  // Random placement of particles
  void initializeRandom(
        int sIndx,
        ParticleSolution* curParticle);

  // Deterministic placement of particles using distribution function
  void initializeDeterministic(
        int sIndx,
        ParticleSolution* curParticle);
  void createDistributionFunction(int sIndx);
  void initialConditionCleanup(
        int sIndx,
        ParticleSolution* curParticle);
  void momentMatching(
        int sIndx,
        ParticleSolution* curParticle);

private:
  PlasmaData* data;            // Input parameters

  double** rho;                // Spatial density profile
  double** u;                  // Velocity profile
  double** T;                  // Temperature profile

  double* mean_rho;
  double* mean_u;
  double* mean_T;

  double** fDist;              // Distribution function for initial placement
  int** fCount;                // Particles per spatial-veloctiy cell
};

#endif
