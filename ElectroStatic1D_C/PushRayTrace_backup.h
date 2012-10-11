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

#ifndef PushRayTrace_h
#define PushRayTrace_h

using namespace std;

#include "PlasmaData.h"
#include "ParticlePush.h"
#include "stdio.h"

class MomentSolution;
class ParticleSolution;

// Type of ray trace particles
const int BAD_CROSS = 0;
const int CROSSED   = 1;
const int TRAPPED   = 2;
const int REVERSED  = 3;

class PushRayTrace : public ParticlePush {
public:
  PushRayTrace(PlasmaData* data);   // Initialization parameters
  ~PushRayTrace();

  ///////////////////////////////////////////////////////////////////
  //
  // Advance particle move
  //
  void push(
        MomentSolution* curMoment,
        MomentSolution* oldMoment,
        ParticleSolution* curParticle,    // Updates values
        ParticleSolution* oldParticle);

  double DtauEstimation(
                                double* E_E,               // E or E_half depending on solve order
                                double x_E,                // particle position
                                double v_E,                // particle velocity
                                int cl_E);                 // particle cell    
/*  double EFieldInterpolation(
        double* E_E,               // E or E_half depending on solve order
        double x_E,                // particle position
        double v_E,                // particle velocity
        int cl_E);                 // particle cell*/

  double DtauCalculationCross(
                              double* E_E,               // E or E_half depending on solve order
                              double x_E,                // particle position
                              int cl_E);                 // particle cell    
/*  double EFieldInterpolationCross(
        double* E_E,               // E or E_half depending on solve order
        double x_E,                // particle position
        int cl_E);                 // particle cell*/

  int particleRayTrace(double dx_half_p);

  int particleReverse(
        int cf_i,
        double v0);

private:
  // Ray trace variables and shortcuts
  int* cell_r;                 // Cell to the right of each cell
  int* cell_l;                 // Cell to the left of each cell
  int* face_c_r;               // Face to the right of each cell
  int* face_c_l;               // Face to the left of each cell
  double* cell_f_c;            // Face positions
  int** c_f_id;
  int** c_f_n;

  // Species loop variables
  double m;                    // Mass of the species
  double q;                    // Charge of the species
  double qm;                   // Charge to mass ratio of species

  double* x_old;               // Old positions of every particle
  double* v_old;               // Old velocities of every particle

  double* x_p;                 // New positions of every particle
  double* v_p;                 // New velocities of every particle

  // Particle loop variables
  double x_sub_f;
  double v_sub_f;
  int cl_sub_f;

  // Subcycling loop variables
  double x_sub_old;            // Old subcycle full time final position 
  double v_sub_old;            // Old subcycle full time final velocity
  int cl_sub_old;              // Old subcycle full time final cell

  double x_sub_half;           // Old subcycle half time final position 
  double v_sub_half;           // Old subcycle half time final velocity
  int cl_sub_half;             // Old subcycle half time final cell

  double dtau;                 // Delta for subcycle moments
  double dt_res;               // Convergence

  // Picard solver loop variables
  double x_k;                  // Starting position for each picard loop
  double v_k;                  // Starting velocity for each picard loop

  int part_cf_xf0;             // Particle reversed or stayed across face
  double x_px;                 // Final position for substep
  double v_px;                 // Final velocity for substep
  int cl_px;                   // Final cell for substep
  int cl_xf;                   // Cell particle crossed into (or face crossed?)
  double c_f_n_xf;             // Normal vector of first cell crossed
  double x_xf;                 // Face position at which particle crossed
  double xf_real;              // Final position of particle
  double E_p;

  // Ray trace loop variables
  double x0;                   // Starting position for ray trace
  int cl;                      // Starting cell for ray trace
  double dx_p;                 // Displacement distance

  // Flags
  int flag_sc;                 // Subcycling termination flag
  int flag_conv;               // Picard termination flag
  int flag_hs;                 // Half step position flag
  int flag_fx;                 // Particle crossed a face
  int flag_rt;                 // Ray trace termination flag
};

#endif
