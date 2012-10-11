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

#include "PlasmaData.h"
#include "stdio.h"

const int LINESIZE = 512;

///////////////////////////////////////////////////////////////////////////
//
// Read parameter input file and store results
//
///////////////////////////////////////////////////////////////////////////

PlasmaData::PlasmaData(const string& inFile)
{
  // Problem parameters
  this->prob_type = 1;
  this->prob_case = 1;
  this->lx = 1.0;

  this->nx = 32;
  this->nv = 60;

  this->tmax = 60.0;
  this->fil_num = 1;
  this->mp = 1.0e-6;
  this->j_avg = 0.0;
  this->stress_flag = 1;
  this->wo = 1;
  this->tol_nlres = 1.0e-6;
  this->tol_nlres_inner = 1.0e-6;
  this->tol_face = 1e-6;
  this->tol_flag = 0;

  this->si_lo_flag = SI_LO_SOLVER;
  this->temp_order = ORDER_1_BACK_EULER;
  this->cc_flag = CHARGE_CONSERVE_ON;
  this->wo_E = LINEAR_INTERPOLATION;

  // Particle parameters
  this->eps_r = 1e-4;
  this->eps_a = 1e-4;
  this->p_size = 1;
  this->NP_ref[ELECTRON] = 100000;  this->NP_ref[ION] = 409600; 
  this->p_type[ELECTRON] = 0;       this->p_type[ION] = 1;
  this->p_name[ELECTRON] = "ELECTRON"; this->p_name[ION] = "ION";
  this->p_color[ELECTRON] = "red";  this->p_color[ION] = "green";
  this->q[ELECTRON] = -1;           this->q[ION] = 1;
  this->m[ELECTRON] = 1;            this->m[ION] = 1834;
  this->sub_nt[ELECTRON] = 1;       this->sub_nt[ION] = 1;
  this->rho0[ELECTRON] = 1.0;       this->rho0[ION] = 1.0;
  this->u0[ELECTRON] = 0.1;         this->u0[ION] = 0.0;
  this->T0[ELECTRON] = 1.0;         this->T0[ION] = 0.0;
  this->alp_r[ELECTRON] = 0.001;    this->alp_r[ION] = 0.0;
  this->alp_u[ELECTRON] = 0.0;      this->alp_u[ION] = 0.0;
  this->alp_T[ELECTRON] = 0.0;      this->alp_T[ION] = 0.0;

  // Simulation parameters
  this->quiet_start_flag = DETERMINISTIC;
  this->dt = 0.5;
  this->jacobi_flag = 0;
  this->ns = 25;
  this->wp = 1;

  this->w_scale = 1;
  this->precon_flag = 0;
  this->picard_flag = 1;
  this->nle_ke_flag = 0;
  this->ient_fact = 1e-2;

  ///////////////////////////////////////////////////////////////////////////
  //
  // Parameters changed in the input deck
  //
  if (inFile != "")
    readInput(inFile);

  // Specialized derivations for problem types
  if (this->prob_type == LANDAU_DAMPING) {
    this->k_wave = (2.0 * M_PI) / this->lx;
    for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
      this->w[sIndx] = sqrt(1.0 / this->m[sIndx]);
      this->tau[sIndx] = 1.0 / this->w[sIndx];
        
    }
  }
  else if (this->prob_type == TWO_STREAM_INSTABILITY) {
    if (this->prob_case == 1) {
      this->k_wave = (2.0 * M_PI) / this->lx;
      for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
        this->w[sIndx] = sqrt(1.0 / this->m[sIndx]);
        this->tau[sIndx] = 1.0 / this->w[sIndx];
      }
    }
    else if (this->prob_case == 2) {
      this->k_wave = (2.0 * M_PI) / this->lx;
      for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
        this->w[sIndx] = sqrt(1.0 / this->m[sIndx]);
        this->tau[sIndx] = 1.0 / this->w[sIndx];
      }
    }
  }
  else if (this->prob_type == ION_ACOUSTIC_SHOCK_WAVE) {
    if (this->prob_case == 1) {
      this->k_wave = (2.0 * M_PI) / this->lx;
      for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
        this->w[sIndx] = sqrt(1.0 / this->m[sIndx]);
        this->tau[sIndx] = 1.0 / this->w[sIndx];
        cout << "prob case = " << prob_case << " \n";
        cout << "for species " << sIndx << " the w_p is = " << w[sIndx] << "\n";
      }
    }
    else if (this->prob_case == 2) {
      this->k_wave = (2.0 * M_PI) / this->lx;
      double dtau_omeg_e = sqrt(this->m[ELECTRON]) / this->w_scale;
      this->dt = 1.0 * dtau_omeg_e;
      for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
        this->w[sIndx] = sqrt(1.0 / this->m[sIndx]);
        this->tau[sIndx] = 1.0 / this->w[sIndx];
          cout << "prob case = " << prob_case << " \n";
          cout << "for species " << sIndx << " the w_p is = " << w[sIndx] << "\n";
      }
    }
    else if (this->prob_case == 3) {
      this->k_wave = (2.0 * M_PI) / this->lx;
      for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
        this->w[sIndx] = sqrt(this->rho0[ELECTRON] / this->m[sIndx]);
        this->tau[sIndx] = 1.0 / this->w[sIndx];
          cout << "prob case = " << prob_case << " \n";
          cout << "for species " << sIndx << " the w_p is = " << w[sIndx] << "\n";
      }
    }
 }

  // Parameters derived from others
  this->nfx = this->nx + 1;
  this->nvc = (this->nv * 2) + 1;
  this->nvf = (this->nv * 2) + 2;
  this->dx = this->lx / (double) this->nx;
  for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
    this->sub_dt[sIndx] = this->dt / this->sub_nt[sIndx];
    this->NP0_tot[sIndx] = this->NP_ref[sIndx];
  }

  // Allocate storage for cell and face location arrays used by all
  this->xpos_node = new double[this->nx];
  this->xpos_face = new double[this->nfx];
  this->cell_f_c = new double[this->nfx];

  this->cell_l = new int[this->nx];
  this->cell_r = new int[this->nx];

  this->face_type = new int[this->nfx];
  this->face_c_l = new int[this->nfx];
  this->face_c_r = new int[this->nfx];
  this->face_n_l = new int[this->nfx];
  this->face_n_r = new int[this->nfx];

  this->face_cell = new int*[2];
  this->c_f_id = new int*[2];
  this->c_f_n = new int*[2];
  for (int i = 0; i < 2; i++) {
    this->face_cell[i] = new int[this->nx];
//Will Taitano edit 06/26/2012: Assigning correct array size for c_f_id
    this->c_f_id[i] = new int[this->nx];    
//    this->c_f_id[i] = new int[this->nfx];
    this->c_f_n[i] = new int[this->nfx];
  }

  for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
    this->v_vec[sIndx] = new double[this->nvc];
    this->v_vec_face[sIndx] = new double[this->nvf];
  }
}

PlasmaData::~PlasmaData()
{
  delete [] this->xpos_node;
  delete [] this->xpos_face;

  for (int sIndx = 0; sIndx < SPECIES; sIndx++) {
    delete [] this->v_vec[sIndx];
    delete [] this->v_vec_face[sIndx];
  }
  delete this->v_vec;
  delete this->v_vec_face;

  delete [] this->cell_f_c;
  delete [] this->cell_l;
  delete [] this->cell_r;

  delete [] this->face_type;
  delete [] this->face_c_l;
  delete [] this->face_c_r;
  delete [] this->face_n_l;
  delete [] this->face_n_r;

  for (int i = 0; i < 2; i++) {
    delete [] this->face_cell[i];
    delete [] this->c_f_id[i];
    delete [] this->c_f_n[i];
  }
  delete [] this->face_cell;
  delete [] this->c_f_id;
  delete [] this->c_f_n;
}

////////////////////////////////////////////////////////////////////////////
//
// Set cell and face positions and connectivity
// Lots of redundancy here
// cell_info_calc.m (which is indexed on 1 not 0)
//
////////////////////////////////////////////////////////////////////////////

void PlasmaData::setCellFaceInformation()
{
  // Cell array information
  for (int cell = 0; cell < this->nx; cell++) {
    this->xpos_node[cell] = this->dx * (cell + 0.5);

    this->cell_l[cell] = cell;
    this->cell_r[cell] = cell + 1;
//Will Taitano edit: 06/26/2012
//For some reason c_f_id[0][0] and c_f_id[0][1] gets overwritten and cannot be changed after the final loop for the Face array information* section below so moving the c_f_id assigning process afterwards.
//    this->c_f_id[0][cell] = this->cell_l[cell];
//    this->c_f_id[1][cell] = this->cell_r[cell];          
  }

  // Face array information*
  for (int face = 0; face < this->nfx; face++) 
  {
    this->xpos_face[face] = this->dx * face;
    this->cell_f_c[face] = this->xpos_face[face];

    this->face_type[face] = INTERFACE;
    this->face_n_l[face] = -1.0;
    this->face_n_r[face] = 1.0;
    this->c_f_n[0][face] = this->face_n_l[face];
    this->c_f_n[1][face] = this->face_n_r[face];

    this->face_c_l[face] = (this->nx + face - 1) % this->nx;
    this->face_c_r[face] = face % this->nx;
    this->face_cell[0][face] = this->face_c_l[face];
    this->face_cell[1][face] = this->face_c_r[face];
    //Will Taitano edit: 06/26/2012 
    //For some magical mystic reason which I cannot figure out why, at the end of this for loop, c_f_n[0][0] and c_f_n[1][0] gets overwritten to some other value. Is there a bizzarre pointer used? I don't know. Please help. For now, I'm just assigning c_f_id after this loop
  }
    // Cell array information
    for (int cell = 0; cell < this->nx; cell++) {
        this->c_f_id[0][cell] = this->cell_l[cell];
        this->c_f_id[1][cell] = this->cell_r[cell];          
    }
  // Boundary faces
  this->face_type[0] = PERIODIC;
  this->face_type[this->nx] = PERIODIC;

/*
  cout << "face_cell" << endl;
  for (int fIndx = 0; fIndx < this->nfx; fIndx++)
    cout << face_cell[0][fIndx] << "   " << face_cell[1][fIndx] << endl;

  cout << "face_type" << endl;
  for (int fIndx = 0; fIndx < this->nfx; fIndx++)
    cout << face_type[fIndx] << endl;

  cout << "cell_l, cell_r" << endl;
  for (int cIndx = 0; cIndx < this->nx; cIndx++)
    cout << cell_l[cIndx] << "   " << cell_r[cIndx] << endl;
*/
/*
  cout << "face_n_l, face_n_r" << endl;
  for (int fIndx = 0; fIndx < this->nfx; fIndx++)
    cout << face_n_l[fIndx] << "   " << face_n_r[fIndx] << endl;
*/
/*
  cout << "cell_f_c" << endl;
  for (int fIndx = 0; fIndx < this->nfx; fIndx++)
    cout << cell_f_c[fIndx] << endl;

  cout << "face_c_l, face_c_r" << endl;
  for (int fIndx = 0; fIndx < this->nfx; fIndx++)
    cout << face_c_l[fIndx] << "   " << face_c_r[fIndx] << endl;

  cout << "c_f_id" << endl;
  for (int cIndx = 0; cIndx < this->nx; cIndx++)
    cout << c_f_id[0][cIndx] << "   " << c_f_id[1][cIndx] << endl;

  cout << "c_f_n" << endl;
  for (int fIndx = 0; fIndx < this->nfx; fIndx++)
    cout << c_f_n[0][fIndx] << "   " << c_f_n[1][fIndx] << endl;
*/
}

////////////////////////////////////////////////////////////////////////////
//
// Parameters which can be altered from an input deck
//
////////////////////////////////////////////////////////////////////////////

void PlasmaData::readInput(const string& inFile)
{
  ifstream inStr(inFile.c_str());
  if (!inStr) {
    cout << "Could not open input file " << inFile << endl;
    exit(-1);
  }

  char inBuf[LINESIZE];
  string keyword;
  string rest;

  while (inStr.getline(inBuf, LINESIZE)) {
    if (inBuf[0] != '#' && inStr.gcount() > 1) {

      getKeyword(inBuf, keyword, rest);
      istringstream line(rest.c_str());
      cout << "KEYWORD: " << keyword << "   " << rest << endl;

      if (keyword == "PROBLEM_TYPE")
        line >> this->prob_type;
      else if (keyword == "PROBLEM_CASE")
        line >> this->prob_case;
      else if (keyword == "QUIET_START")
        line >> this->quiet_start_flag;
      else if (keyword == "SYSTEM_LENGTH")
        line >> this->lx;
      else if (keyword == "NUMBER_OF_POS_CELLS")
        line >> this->nx;
      else if (keyword == "NUMBER_OF_VEL_CELLS")
        line >> this->nv;

      else if (keyword == "TIME_MAX")
        line >> this->tmax;
      else if (keyword == "TIME_STEP")
        line >> this->dt;

      else if (keyword == "NUMBER_OF_FILTERING_OPERATIONS")
        line >> this->fil_num;
      else if (keyword == "CHARGE_CONSERVE_FLAG")
        line >> this->cc_flag;
      else if (keyword == "LO_SOLVER")
        line >> this->si_lo_flag;
      else if (keyword == "TEMPORAL_ORDER")
        line >> this->temp_order;
      else if (keyword == "W_SCALE")
        line >> this->w_scale;

      else if (keyword == "NUMBER_OF_PARTICLE_SPECIES")
        line >> this->p_size;
      else if (keyword == "PARTICLE_TYPE")
        line >> this->p_type[0] >> this->p_type[1];
      else if (keyword == "PARTICLE_REFERENCE_NUMBER")
        line >> this->NP_ref[0] >> this->NP_ref[1];
      else if (keyword == "SPECIES_CHARGE")
        line >> this->q[0] >> this->q[1];
      else if (keyword == "SPECIES_MASS")
        line >> this->m[0] >> this->m[1];
      else if (keyword == "SPECIES_SUBCYCLE")
        line >> this->sub_nt[0] >> this->sub_nt[1];

      else if (keyword == "DENSITY_UNPERTURBED")
        line >> this->rho0[0] >> this->rho0[1];
      else if (keyword == "VELOCITY_UNPERTURBED")
        line >> this->u0[0] >> this->u0[1];
      else if (keyword == "TEMPERATURE_UNPERTURBED")
        line >> this->T0[0] >> this->T0[1];
      else if (keyword == "DENSITY_PERTURBED")
        line >> this->alp_r[0] >> this->alp_r[1];
      else if (keyword == "VELOCITY_PERTURBED")
        line >> this->alp_u[0] >> this->alp_u[1];
      else if (keyword == "TEMPERATURE_PERTURBED")
        line >> this->alp_T[0] >> this->alp_T[1];
    }
  }
}


/////////////////////////////////////////////////////////////////////////////
//
// Keywords start in position 0 and are delimited by white space
//
/////////////////////////////////////////////////////////////////////////////

void PlasmaData::getKeyword(char* inBuf, string& keyword, string& rest)
{
  string line(inBuf);
  string::size_type keyPos = line.find(' ');
  keyword = line.substr(0, keyPos);
  rest = line.substr(keyPos + 1);
}
