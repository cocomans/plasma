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

#ifndef PlasmaData2_h
#define PlasmaData2_h

#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
class PlasmaData;

using namespace std;

// Species type
const int ELECTRON  = 0;
const int ION       = 1;
const int SPECIES   = 2;

// Face type
const int PERIODIC  = 0;
const int INTERFACE = 1;
const int OUTFLOW   = 2;

// Problem type
const int LANDAU_DAMPING = 0;
const int TWO_STREAM_INSTABILITY = 1;
const int ION_ACOUSTIC_SHOCK_WAVE = 2;

// Particle initialization type
const int RANDOM = 0;
const int DETERMINISTIC = 1;

// Electric field interpolation methods
const int NGP_INTERPOLATION = 0;                // Nearest grid point
const int LINEAR_INTERPOLATION = 1;
const int BSPLINE_INTERPOLATION = 2;

// Time order of discretization being used
const int ORDER_1_BACK_EULER = 1;
const int ORDER_2_CRANK_NICOLSON = 2;

// Low order solver
const int JFNK_LO_SOLVER = 0;
const int SI_LO_SOLVER   = 1;

// Particle pushing
const int CHARGE_CONSERVE_OFF = 0;      // Non charge conserving
const int CHARGE_CONSERVE_ON  = 1;      // Charge conserving

class PlasmaData2 {
public:
  PlasmaData2(const string& inFile);

  PlasmaData2(PlasmaData* pdata);

  ~PlasmaData2();

  void   readInput(const string& inFile);
  void   getKeyword(char* inBuf, string& keyword, string& rest);

  void   setCellFaceInformation();

  //////////////////////////////////////////////////////////////////////////
  //
  // Problem parameters
  //
  int prob_type;           // Landau, Two stream, Ion acoustic shock wave
  int prob_case;           // Cases of the problem type
  double lx;               // System length
  int nx;                  // Number of cells in system
  int nfx;                 // Number of position space faces
  int nv;                  // Number of velocity space cells each direction
  int nvc;                 // Number of velocity space cells total
  int nvf;                 // Number of velocity space face cells total

  double tmax;             // Maximum duration of simulation
  double tol_nlres;        // Tolerance for non-linear residual convergence
  double tol_nlres_inner;  // Tolerance for inner Newton residual
  double tol_face;         // Distance tolerance to face
  int tol_flag;            // 0 relative tolerance, 1 absolute tolerance

  int fil_num;             // Number of filtering operations on moment quantity
  double mp;               // Perturbation of non-linear residual function
  int cc_flag;             // Charge conserving adaptive orbit integrator flag
  int si_lo_flag;          // Low order solver flag
  int temp_order;          // Discretization order of time flag
  double j_avg;            // Initial average current
  int stress_flag;         // Stress tensor or pressure closure flag
  int wo;                  // Order of particle weighting method
  int wo_E;                // Order of electric field weighting method
  double eps_r;            // Relative tolerance for adaptive orbit integrator
  double eps_a;            // Absolute tolerance for adaptive orbit integrator
  double k_wave;           // Wave vector for the perturbation

  //////////////////////////////////////////////////////////////////////////
  //
  // Particle parameters
  //
  int p_size;              // Number of particle species
  int NP_ref[SPECIES];     // Reference number of particle array
  int p_type[SPECIES];     // ELECTRON or ION
  string p_name[SPECIES];  // ELECTRON or ION
  string p_color[SPECIES]; // ELECTRON or ION
  double q[SPECIES];       // Charge array
  double m[SPECIES];       // Mass array
  int sub_nt[SPECIES];     // Number of subcycle array
                           // for non-charge conserving scheme
  double rho0[SPECIES];    // Unperturbed density quantity
  double u0[SPECIES];      // Unperturbed fluid velocity quantity
  double T0[SPECIES];      // Unperturbed fluid temperature quantity
  double alp_r[SPECIES];   // Perturbed density quantity
  double alp_u[SPECIES];   // Perturbed fluid velocity quantity
  double alp_T[SPECIES];   // Perturbed fluid temperature quantity
  double w[SPECIES];       // Plasma frequency array
  double tau[SPECIES];     // Plasma wave time scale array

  double mpp[SPECIES];     // Mass per particle
  int NP0_tot[SPECIES];    // Number of particles in system per species
  double mean_rho_NP[SPECIES]; // Average density per species

  double v_th[SPECIES];    // Average density per species
  double vmax[SPECIES];    // Average density per species
  double vmin[SPECIES];    // Average density per species
  double dv[SPECIES];      // Average density per species
  double* v_vec[SPECIES];  // Average density per species
  double* v_vec_face[SPECIES];// Average density per species
  
  //////////////////////////////////////////////////////////////////////////
  //
  // Simulation parameters
  //
  int quiet_start_flag;    // Placement of particles noisy or quiet
  double dt;               // Global time step size
  double sub_dt[SPECIES];  // Global subcycle time step
  double dx;               // Spatial cell size
  int jacobi_flag;         // JACOBI or GAUSS-SEIDEL
  int ns;                  // Number of Gauss-Seidel or Jacobi sweep
  int wp;                  // Plasma frequency for scaling

  double w_scale;
  int picard_flag;         // Use outer and inner picard
  int precon_flag;         // Use physics based preconditioner for jfnk
  int nle_ke_flag;         // Use nonlinear eliminated kinetic for jfnk
  double ient_fact;        // Constant for inexact Newton tolerance

  //////////////////////////////////////////////////////////////////////////
  //
  // Cell and face information
  //
  double* xpos_node;       // Cell centroid
  double* xpos_face;       // Face location
  double* cell_f_c;        // Face centroid

  int* cell_l;             // Left face of cell
  int* cell_r;             // Right face of cell

  int* face_type;          // Periodic, interface or outflow
  int* face_c_l;           // Left cell of face
  int* face_c_r;           // Right cell of face
  int* face_n_l;           // Left cell of face normalized
  int* face_n_r;           // Right cell of face normalized

  int** face_cell;         // Face-cell connectivity (left and right)
  int** c_f_id;            // Face list (left and right)
  int** c_f_n;             // Face normal (left and right)

  //////////////////////////////////////////////////////////////////////////
  //
  // MPI Information
  //
  int myid; // MPI Thread ID
  int numprocs; // Number of MPI Threads
};

#endif
