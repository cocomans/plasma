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

#ifndef AmpereVlasov_h
#define AmpereVlasov_h

#include <vector>
#include <string>

class PlasmaData;
class Particles;
class ParticleSolution;
class MomentSolution;
class LOSolver;
class DhatConsistency;
class PlasmaOut;

using namespace std;

class AmpereVlasov {
public:
  AmpereVlasov(string filename);
  ~AmpereVlasov();

  // Preprocessing initial conditions
  void preprocess();

  // Simulation loop over time
  void simulation();

  // Summary of simulation time step
  void totalEnergyCalculation(int tIndx);

private:
  // Input object
  PlasmaData* data;

  // Particle initializer and pusher
  Particles* particles;

  // Particle double buffer for maintaining old and new particle information
  ParticleSolution* curParticle;
  ParticleSolution* oldParticle;

  // Moment double buffer for maintaining old and new moment information
  MomentSolution* curMoment;
  MomentSolution* oldMoment;

  // Low order solver
  LOSolver* loSolver;

  // Consistency calculator
  DhatConsistency* Dhat;

  // Output
  PlasmaOut* plasmaOut;

  // Summary information over time
  int numberOfTimeSteps;
  double* time_vec;                      // Normalized time
  double* E_tot;                         // Total electric
  double* mv_tot;                        // Total momentum
  double* E_diff;                        // Difference in total electric field
  double* mv_diff;                       // Difference in total momentum

  double* EE;                            // Electric field energy
  double** KE;                           // Kinetic energy by species
  double** cc_rms;                       // Charge conservation by species
};

#endif
