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

#include "Particles.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"
#include "InitialLandau.h"
#include "InitialTSI.h"
#include "InitialIASW.h"
#include "PushSimple.h"
#include "PushRayTrace.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

Particles::Particles(PlasmaData* data)
{
  this->data = data;

  // Create the particle initializer
  if (this->data->prob_type == LANDAU_DAMPING) {
    this->particleInitialize = new InitialLandau(this->data);
    this->data->wp = 1;
  }
  else if (this->data->prob_type == TWO_STREAM_INSTABILITY) {
    this->particleInitialize = new InitialTSI(this->data);
    this->data->wp = 1;
  }
  else {
    this->particleInitialize = new InitialIASW(this->data);
    this->data->wp = 1;
  }

  // Create the particle pusher
  if (this->data->cc_flag == CHARGE_CONSERVE_OFF) {
    this->particlePush = new PushSimple(this->data);
  }
  else if (this->data->cc_flag == CHARGE_CONSERVE_ON) {
    this->particlePush = new PushRayTrace(this->data);
  }
}

Particles::~Particles()
{
  delete this->particleInitialize;
  delete this->particlePush;
}

///////////////////////////////////////////////////////////////////////
//
// Initialize particle position and velocity
//
///////////////////////////////////////////////////////////////////////

void Particles::initialCondition(ParticleSolution* curParticle)
{
  cout << "Particles::initialCondition()" << endl;

  // Reset weighting to 0th order for initialization
  int wo_sav = this->data->wo;
  this->data->wo = 0;

  // Initialize particle position and velocity
  this->particleInitialize->initialize(curParticle);

  // Redefine the plasma wave time scale by scaling with w_p
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    this->data->tau[sIndx] /= this->data->wp;
  }
  this->data->wo = wo_sav;
}

///////////////////////////////////////////////////////////////////////
//
// Push particle position and velocity
//
///////////////////////////////////////////////////////////////////////

void Particles::pushParticles(
                      MomentSolution* curMoment,
                      MomentSolution* oldMoment,
                      ParticleSolution* curParticle,    // Updates values
                      ParticleSolution* oldParticle)

{
  cout << "Particles::particlePush()" << endl;
  this->particlePush->push(curMoment, oldMoment,
                           curParticle, oldParticle);
}
