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

#ifndef MomentSolution_h
#define MomentSolution_h

using namespace std;

class PlasmaData;
class ParticleSolution;

//
// Container for the moment solution for double buffering through time steps
//
class MomentSolution {
public:
  MomentSolution(PlasmaData2* data);
  ~MomentSolution();

  // Calculate moments on all particles
  void momentCalculation(
        int sIndx,               // Species index
        int numberOfParticles,
        double* x,               // Positions
        double* v);              // Velocities

  // Calculate moment increment and accumulate for one particle
  void tallyMomentCalculation(
        int sIndx,               // Species index
        double x,                // Particle location
        double v,                // Particle velocity
        double factor);          // Factor on velocity

  // Not getting the same answers as the Matlab code in low digits
  void electricFieldCalculation();
  void electricFieldTridiagonal();

  void scaleMoment(int sIndx);
  void scaleMoment(int sIndx, double factor);
  void smoothMoment(int sIndx);
  void printMoment(int sIndx);

  // Copy complete moment for initialization
  void copyAllMoment(MomentSolution* m);

  // Copy only information needed for time step pass
  void copyParticleSolution(ParticleSolution* m);
  void copyMomentSolution(MomentSolution* m);

  double* getR(int sIndx)       { return this->R[sIndx]; }
  double* getRU(int sIndx)      { return this->RU[sIndx]; }
  double* getS2(int sIndx)      { return this->S2[sIndx]; }
  double* getS3(int sIndx)      { return this->S3[sIndx]; }

  double* getU(int sIndx)       { return this->U[sIndx]; }
  double* getT(int sIndx)       { return this->T[sIndx]; }
  double* getP(int sIndx)       { return this->P[sIndx]; }
  double* getETh(int sIndx)     { return this->E_TH[sIndx]; }
  double* getETot(int sIndx)    { return this->E_TOT[sIndx]; }
  double* getE()                { return this->E; }

  double* getRFace(int sIndx)   { return this->R_f[sIndx]; }
  double* getRUCell(int sIndx)  { return this->RU[sIndx]; }

  // Reduce Moment data accross MPI nodes
  void  ReduceMoment();

  // Broadcast Moment Solution to all MPI nodes
  void BCastMoment();

private:
  PlasmaData2* data;     // Parameters and data

  double** R;           // Cell density (Oth moment)
  double** RU;          // Face momentum (1st moment)
  double** S2;          // Cell stress (2nd moment)
  double** S3;          // Cell third moment

  double** U;           // Cell velocity
  double** T;           // Cell temperature
  double** P;           // Cell pressure

  double** E_TH;        // Cell thermal energy moment
  double** E_TOT;
  double* E;            // Face electric field

  double** R_f;         // Face density
  double** RU_c;        // Cell momentum
};

#endif
