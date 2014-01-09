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

#include "PlasmaData2.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"
#include "stdio.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

MomentSolution::MomentSolution(PlasmaData2* data)
{
  this->data = data;

  this->R = new double*[this->data->p_size];
  this->R_f = new double*[this->data->p_size];
  this->RU = new double*[this->data->p_size];
  this->RU_c = new double*[this->data->p_size];
  this->S2 = new double*[this->data->p_size];
  this->S3 = new double*[this->data->p_size];

  this->U = new double*[this->data->p_size];
  this->T = new double*[this->data->p_size];
  this->P = new double*[this->data->p_size];
  this->E_TH = new double*[this->data->p_size];
  this->E_TOT = new double*[this->data->p_size];

  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    this->R[sIndx] = new double[this->data->nx];
    this->R_f[sIndx] = new double[this->data->nfx];
    this->RU[sIndx] = new double[this->data->nfx];
    this->RU_c[sIndx] = new double[this->data->nx];

    this->S2[sIndx] = new double[this->data->nx];
    this->S3[sIndx] = new double[this->data->nx];

    this->U[sIndx] = new double[this->data->nx];
    this->T[sIndx] = new double[this->data->nx];
    this->P[sIndx] = new double[this->data->nx];

    this->E_TH[sIndx] = new double[this->data->nx];
    this->E_TOT[sIndx] = new double[this->data->nx];

    // Initialize
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
      this->R[sIndx][cIndx] = 0.0;
      this->RU_c[sIndx][cIndx] = 0.0;
      this->S2[sIndx][cIndx] = 0.0;
      this->S3[sIndx][cIndx] = 0.0;
      this->E_TH[sIndx][cIndx] = 0.0;
    }
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
      this->R_f[sIndx][fIndx] = 0.0;
      this->RU[sIndx][fIndx] = 0.0;
    }
  }
  this->E = new double[this->data->nfx];
}

MomentSolution::~MomentSolution()
{
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    delete [] this->R[sIndx];
    delete [] this->R_f[sIndx];
    delete [] this->RU[sIndx];
    delete [] this->RU_c[sIndx];

    delete [] this->S2[sIndx];
    delete [] this->S3[sIndx];

    delete [] this->U[sIndx];
    delete [] this->T[sIndx];
    delete [] this->P[sIndx];
    delete [] this->E_TH[sIndx];
    delete [] this->E_TOT[sIndx];
  }
  delete [] this->R;
  delete [] this->R_f;
  delete [] this->RU;
  delete [] this->RU_c;
  delete [] this->S2;
  delete [] this->S3;
  delete [] this->U;
  delete [] this->T;
  delete [] this->P;
  delete [] this->E_TH;
  delete [] this->E_TOT;
  delete [] this->E;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate the electric field via Poisson equation for the 0th time step
// Inconsistencies in the formulation of the A matrix (not symmetric)
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::electricFieldCalculation()
{
  double dx2 = this->data->dx * this->data->dx;
  double wp2 = this->data->wp * this->data->wp;

  double** A = new double*[this->data->nx];
  for (int i = 0; i < this->data->nx; i++) {
    A[i] = new double[this->data->nx];
    for (int j = 0; j < this->data->nx; j++)
      A[i][j] = 0.0;
  }
  double* rhs = new double[this->data->nx];
  double* phi = new double[this->data->nx];

  // Tally of total number density
  double* rho_tot = new double[this->data->nx];
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    rho_tot[cIndx] = 0.0;

  // Collect density multiplied by charge per spatial cell per species
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      rho_tot[cIndx] += this->data->q[sIndx] * this->R[sIndx][cIndx];

  // If only one species of particle increment by 1
  if (this->data->p_size == 1)
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      rho_tot[cIndx] += 1.0;

  // Set up vector for right hand side of linear system Ax = b
  if (this->data->prob_type != ION_ACOUSTIC_SHOCK_WAVE) {
    for (int k = 0; k < this->data->fil_num; k++)
      smoothCellFilter(this->data->nx, rho_tot);
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      rhs[cIndx] = rho_tot[cIndx] * dx2;
  } 
  else {
    for (int k = 0; k < this->data->fil_num; k++)
      smoothCellFilter(this->data->nx, rho_tot);
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
      rhs[cIndx] = rho_tot[cIndx] * dx2 * wp2;
  }

  // Form the coefficient matrix for poisson solver
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    if (cIndx == 0) {
      A[cIndx][cIndx] = 3;
      A[cIndx][cIndx + 1] = -1;
    }
    else if (cIndx == (this->data->nx - 1)) {
      A[cIndx][cIndx] = 3;
      A[cIndx][cIndx - 1] = -1;
    }
    else {
      A[cIndx][cIndx - 1] = -1;
      A[cIndx][cIndx] = 2;
      A[cIndx][cIndx + 1] = -1;
    }
  }

  // Calculate the new time electron momentum from direct inversion
  double residual;
  int flag;
  GaussElim(A, this->data->nx, rhs, this->data->nx, phi, residual, flag);

  // Calculate electric field
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    if (fIndx == 0 || fIndx == this->data->nx)
      this->E[fIndx] = -(phi[0] - phi[this->data->nx - 1]) / this->data->dx;
    else
      this->E[fIndx] = -(phi[fIndx] - phi[fIndx-1]) / this->data->dx;
  }

  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    delete [] A[cIndx];
  delete [] A;
  delete [] rhs;
  delete [] phi;
  delete [] rho_tot;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate the electric field via Poisson equation for the 0th time step
// Thomas tridiagonal algorithm
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::electricFieldTridiagonal()
{
  double dx2 = this->data->dx * this->data->dx;
  double wp2 = this->data->wp * this->data->wp;
  double two_dx = 2.0 * this->data->dx;

  // Allocate vectors for Ax = v solver
  // Since sparse matrix A is really tridiagonal we use simplified
  // Gaussian elimination with forward sweep to remove unknowns
  // and then back substitution (Thomas algorithm)
  double* a = new double[this->data->nfx];
  double* b = new double[this->data->nfx];
  double* c = new double[this->data->nfx];
  double* v = new double[this->data->nfx];
  double* phi = new double[this->data->nfx];
  double* rho_tot = new double[this->data->nfx];

  // Tally of total number density
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    rho_tot[fIndx] = 0.0;

  // Collect density multiplied by charge per spatial cell per species
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++)
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      rho_tot[fIndx] += this->data->q[sIndx] * this->R_f[sIndx][fIndx];

  // If only one species of particle increment by 1
  if (this->data->p_size == 1)
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      rho_tot[fIndx] += 1.0;

  // Set up vector for right hand side of linear system Ax = b
  if (this->data->prob_type != ION_ACOUSTIC_SHOCK_WAVE) {
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      v[fIndx] = rho_tot[fIndx] * dx2;
  } 
  else {
    for (int k = 0; k < this->data->fil_num; k++)
      smoothFaceFilter(this->data->nfx, rho_tot);
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      v[fIndx] = rho_tot[fIndx] * dx2 * wp2;
  }

  // Set up tridiagonal matrix
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    a[fIndx] = -1;
    b[fIndx] = 2;
    c[fIndx] = -1;
  }

  // Solution of Ax = b where x is phi
  solveMatrix2(this->data->nfx, a, b, c, v, phi);

  // Calculate electric field
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    if (fIndx == 0 || fIndx == this->data->nx)
      this->E[fIndx] = -(phi[1] - phi[this->data->nx - 1]) / two_dx;
    else
      this->E[fIndx] = -(phi[fIndx+1] - phi[fIndx-1]) / two_dx;
  }

  delete [] a;
  delete [] b;
  delete [] c;
  delete [] v;
  delete [] phi;
  delete [] rho_tot;
}

///////////////////////////////////////////////////////////////////////
//
// Moment calculation
// moment_calc.m
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::momentCalculation(
                        int sIndx,                 // Species index
                        int numberOfParticles,
                        double* x,                 // Particle location
                        double* v)                 // Particle velocity
{
  double dx_recip = 1.0 / this->data->dx;
  double c_part, r_part, l_part;
  int cl, cl_r, cl_l;
  if (numberOfParticles == 0) {
    cout << "No particles for specie " << sIndx << endl;
    return;
  }

  // Tally the moment quantities for all particles
  // Increments R, RU, S2, S3 for cell and face
  for (int pIndx = 0; pIndx < numberOfParticles; pIndx++)
    tallyMomentCalculation(sIndx, x[pIndx], v[pIndx], 1.0);

  // Average velocity of particles in a cell
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    this->U[sIndx][cIndx] = this->RU_c[sIndx][cIndx] / this->R[sIndx][cIndx];

  // Accumulate thermal energy moment
  for (int pIndx = 0; pIndx < numberOfParticles; pIndx++) {
    int cl = (int) (x[pIndx] / this->data->dx);
    if (cl == this->data->nx)
      cl = this->data->nx - 1;
    int cl_r = this->data->face_c_r[this->data->cell_r[cl]];
    int cl_l = this->data->face_c_l[this->data->cell_l[cl]];

    // Accumulate thermal energy moment
    double vu1 = v[pIndx] - this->U[sIndx][cl];
    double vu2 = vu1 * vu1;

    if (x[pIndx] >= this->data->xpos_node[cl]) {
      // Particle location is past cell centroid, share with cell on right
      r_part = (x[pIndx] - this->data->xpos_node[cl]) * dx_recip;
      c_part = 1.0 - r_part;

      this->E_TH[sIndx][cl] += vu2 * c_part;
      this->E_TH[sIndx][cl_r] += vu2 * r_part;
    }
    else {
      // Particle location is before cell centroid, share with cell on left
      l_part = (this->data->xpos_node[cl] - x[pIndx]) * dx_recip;
      c_part = 1.0 - l_part;

      this->E_TH[sIndx][cl] += vu2 * c_part;
      this->E_TH[sIndx][cl_l] += vu2 * l_part;
    }
  }

  // Scale the resulting quantities
  scaleMoment(sIndx);
}

///////////////////////////////////////////////////////////////////////
//
// Moment calculation accumulator taking information for a single particle
// and accumulating the result in arrays
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::tallyMomentCalculation(
                        int sIndx,                 // Species index
                        double x,                  // Particle location
                        double v,                  // Particle velocity
                        double factor)             // Velocity factor
{
  double s = 1.0 / this->data->dx;
  double dx_recip = 1.0 / this->data->dx;
  double dx2_recip = 1.0 / (this->data->dx * 2.0);
  double dxd2 = this->data->dx / 2.0;
  double c_part, r_part, l_part;
  int cl, cl_r, cl_l, cl_f, cl_f_l, cl_f_r;

  double* rho_f = this->R_f[sIndx];
  double* rho = this->R[sIndx];
  double* s1_f = this->RU[sIndx];
  double* s1 = this->RU_c[sIndx];
  double* s2 = this->S2[sIndx];
  double* s3 = this->S3[sIndx];

  // Cell index and cells on either side
  cl = (int) (x / this->data->dx);
  if (cl == this->data->nx)
    cl--;
  cl_r = this->data->face_c_r[this->data->cell_r[cl]];
  cl_l = this->data->face_c_l[this->data->cell_l[cl]];

  double v1 = v * factor;
  double v2 = v * v * factor;
  double v3 = v * v * v * factor;

  double eps_w = x - this->data->xpos_node[cl];
  double s_1 = s * (dxd2 - eps_w) * dx_recip;
  double s_2 = s - s_1;
  double a = s * ((dxd2 - eps_w) * (dxd2 - eps_w)) * dx2_recip;
  double b = s_2 * (dxd2 + eps_w) - 
             s * ((dxd2 + eps_w) * (dxd2 + eps_w)) * dx2_recip;
  double diff = 1.0 - a - b;

  // Tally variables at cell center using m = 2, b-spline function
  rho[cl_l] += a;
  rho[cl_r] += b;
  rho[cl] += diff;

  s2[cl_l] += v2 * a;
  s2[cl_r] += v2 * b;
  s2[cl] += v2 * diff;

  // Tally variables at cell center using m = 1, b-spline function
  if (x >= this->data->xpos_node[cl]) {

    // Particle is located to the right of the cell centroid
    r_part = (x - this->data->xpos_node[cl]) * dx_recip;
    c_part = 1.0 - r_part;

    s1[cl] += v1 * c_part;
    s1[cl_r] += v1 * r_part;
    s3[cl] += v3 * c_part;
    s3[cl_r] += v3 * r_part;
  }
  else {
    // Particle is located to the left of the cell centroid
    l_part = (this->data->xpos_node[cl] - x) * dx_recip;
    c_part = 1.0 - l_part;

    s1[cl] += v1 * c_part;
    s1[cl_l] += v1 * l_part;
    s3[cl] += v3 * c_part;
    s3[cl_l] += v3 * l_part;
  }

  // Set the faces for this particle location
  if (x >= this->data->xpos_node[cl]) {
    cl_f = cl + 1;
    cl_f_r = cl_f + 1;
    cl_f_l = cl_f - 1;
    if (cl_f == this->data->nx)
      cl_f_r = 1;
  }
  else {
    cl_f = cl;
    cl_f_r = cl_f + 1;
    cl_f_l = cl_f - 1;
    if (cl_f == 0)
      cl_f_l = this->data->nx - 1;
  }

  // Tally surface density using 2nd order shape function
  eps_w = x - this->data->xpos_face[cl_f];
  s_1 = s * (dxd2 - eps_w) * dx_recip;
  s_2 = s - s_1;
  a = s * ((dxd2 - eps_w) * (dxd2 - eps_w)) * dx2_recip;
  b = s_2 * (dxd2 + eps_w) - s * ((dxd2 + eps_w) * (dxd2 + eps_w)) * dx2_recip;
  diff = 1.0 - a - b;

  // This extra confusion is because while first cell and last cell are
  // next to each other, the first face and the last face are identical
  // so on the boundary values must be added in two places
  if (cl_f == 0) {                          // First face
    rho_f[cl_f] += diff;
    rho_f[this->data->nx] += diff;
    rho_f[cl_f_r] += b;
    rho_f[cl_f_l] += a;
  }
  else if (cl_f == 1) {                     // Second face
    rho_f[cl_f] += diff;
    rho_f[cl_f_r] += b;
    rho_f[cl_f_l] += a;
    rho_f[this->data->nx] += a;
  }
  else if (cl_f == (this->data->nx - 1)) {  // Next to last face
    rho_f[cl_f] += diff;
    rho_f[cl_f_r] += b;
    rho_f[0] += b;
    rho_f[cl_f_l] += a;
  }
  else if (cl_f == this->data->nx) {        // Last face
    rho_f[cl_f] += diff;
    rho_f[0] += diff;
    rho_f[cl_f_r] += b;
    rho_f[cl_f_l] += a;
  }
  else {                                    // Inner faces
    rho_f[cl_f] += diff;
    rho_f[cl_f_r] += b;
    rho_f[cl_f_l] += a;
  }

  // Tally surface current using 1st order shape function
  if (x >= this->data->xpos_face[cl_f]) {
    // Particle is located to the right of the face
//  Will Taitano Edit 06/25/2012
//      r_part  = (this->data->xpos_face[cl_f] + this->data->dx - x)*dx_recip;
    r_part = (x - this->data->xpos_face[cl_f]) * dx_recip;
    c_part = 1.0 - r_part;

    if (cl_f == 0) {                          // First face
      s1_f[cl_f] += v1 * c_part;
      s1_f[this->data->nx] += v1 * c_part;
      s1_f[cl_f_r] += v1 * r_part;
    }
    else if (cl_f == (this->data->nx - 1)) {  // Next to last face
      s1_f[cl_f] += v1 * c_part;
      s1_f[cl_f_r] += v1 * r_part;
      s1_f[0] += v1 * r_part;
    }
    else if (cl_f == this->data->nx) {        // Last face
      s1_f[cl_f] += v1 * c_part;
      s1_f[0] += v1 * c_part;
      s1_f[cl_f_r] += v1 * r_part;
    }
    else {                                    // Inner faces
      s1_f[cl_f] += v1 * c_part;
      s1_f[cl_f_r] += v1 * r_part;
    }
  }
  else {
    // Particle is located to the left of the face
//  Will Taitano Edit 06/25/2012
//      l_part  = (x + this->data->dx - this->data->xpos_face[cl_f])*dx_recip;
    l_part = (this->data->xpos_face[cl_f] - x) * dx_recip;
    c_part = 1.0 - l_part;

    if (cl_f == this->data->nx) {             // Last face
      s1_f[cl_f] += v1 * c_part;
      s1_f[0] += v1 * c_part;
      s1_f[cl_f_l] += v1 * l_part;
    }
    else if (cl_f == 0) {                     // First face
      s1_f[cl_f] += v1 * c_part;
      s1_f[cl_f_l] += v1 * l_part;
      s1_f[this->data->nx - 1] += v1 * l_part;
    }
    else if (cl_f == 1) {                     // Second face
      s1_f[cl_f] += v1 * c_part;
      s1_f[cl_f_l] += v1 * l_part;
      s1_f[this->data->nx] += v1 * l_part;
    }
    else {                                    // Inner faces
      s1_f[cl_f] += v1 * c_part;
      s1_f[cl_f_l] += v1 * l_part;
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
// Scale the moment quantities
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::scaleMoment(int sIndx)
{
  // Scale the resulting quantities
  double dxmpp = this->data->mpp[sIndx] / this->data->dx;

  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    this->R[sIndx][cIndx] *= dxmpp;
    this->RU_c[sIndx][cIndx] *= dxmpp;
    this->S2[sIndx][cIndx] *= dxmpp;
    this->S3[sIndx][cIndx] = dxmpp * (this->S3[sIndx][cIndx] / 2.0);
    this->E_TH[sIndx][cIndx] = dxmpp * (this->E_TH[sIndx][cIndx] / 2.0);
    this->T[sIndx][cIndx] = 2.0 * this->data->m[sIndx] * 
                            this->E_TH[sIndx][cIndx] / this->R[sIndx][cIndx];
    this->E_TOT[sIndx][cIndx] = this->S2[sIndx][cIndx] / 2.0;
    this->U[sIndx][cIndx] = this->RU_c[sIndx][cIndx] / this->R[sIndx][cIndx];
    this->P[sIndx][cIndx] = this->S2[sIndx][cIndx] - 
                      (2.0 * this->RU_c[sIndx][cIndx] * this->U[sIndx][cIndx]) +
                      (this->U[sIndx][cIndx] * this->U[sIndx][cIndx]);
  }

  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    this->R_f[sIndx][fIndx] *= dxmpp;
    this->RU[sIndx][fIndx] *= dxmpp;
  }
}

void MomentSolution::scaleMoment(int sIndx, double factor)
{
  // Scale the resulting quantities
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
    this->R[sIndx][cIndx] *= factor;
    this->RU_c[sIndx][cIndx] *= factor;
    this->S2[sIndx][cIndx] *= factor;
  }

  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    this->R_f[sIndx][fIndx] *= factor;
    this->RU[sIndx][fIndx] *= factor;
  }
}

///////////////////////////////////////////////////////////////////////
//
// Filter moment quantities to smooth
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::smoothMoment(int sIndx)
{
  for (int fil = 0; fil < this->data->fil_num; fil++) {
    smoothCellFilter(this->data->nx, this->R[sIndx]);
    smoothCellFilter(this->data->nx, this->S2[sIndx]);
    smoothFaceFilter(this->data->nfx, this->RU[sIndx]);
  }
}

///////////////////////////////////////////////////////////////////////
//
// Print moment information
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::printMoment(int sIndx)
{
  cout << setprecision(16);
  cout << endl << "rho =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->R[sIndx][cIndx] << endl;

  cout << endl << "rho_f =" << endl;
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    cout << this->R_f[sIndx][fIndx] << endl;

  cout << endl << "u =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->U[sIndx][cIndx] << endl;

  cout << endl << "s1 =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->RU_c[sIndx][cIndx] << endl;

  cout << endl << "s1_f" << endl;
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    cout << this->RU[sIndx][fIndx] << endl;

  cout << endl << "s2 =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->S2[sIndx][cIndx] << endl;

  cout << endl << "s3 =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->S3[sIndx][cIndx] << endl;

  cout << endl << "P =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->P[sIndx][cIndx] << endl;

  cout << endl << "T =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->T[sIndx][cIndx] << endl;

  cout << endl << "E_th =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->E_TH[sIndx][cIndx] << endl;

  cout << endl << "E_tot =" << endl;
  for (int cIndx = 0; cIndx < this->data->nx; cIndx++)
    cout << this->E_TOT[sIndx][cIndx] << endl;
}

///////////////////////////////////////////////////////////////////////
//
// Copy given MomentSolution into this, basic information is set in constructor
// Needed for first time step
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::copyAllMoment(MomentSolution* m)
{
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    double* copyr = m->getR(sIndx);
    double* copyrface = m->getRFace(sIndx);
    double* copyru = m->getRU(sIndx);
    double* copyrucell = m->getRUCell(sIndx);
    double* copys2 = m->getS2(sIndx);
    double* copys3 = m->getS3(sIndx);
    double* copyE_th = m->getETh(sIndx);
    double* copyE_tot = m->getETot(sIndx);
    double* copyu = m->getU(sIndx);
    double* copyT = m->getT(sIndx);
    double* copyP = m->getP(sIndx);

    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
      this->R[sIndx][cIndx] = copyr[cIndx];
      this->RU_c[sIndx][cIndx] = copyrucell[cIndx];
      this->S2[sIndx][cIndx] = copys2[cIndx];
      this->S3[sIndx][cIndx] = copys3[cIndx];
      this->E_TH[sIndx][cIndx] = copyE_th[cIndx];
      this->E_TOT[sIndx][cIndx] = copyE_tot[cIndx];
      this->U[sIndx][cIndx] = copyu[cIndx];
      this->T[sIndx][cIndx] = copyT[cIndx];
      this->P[sIndx][cIndx] = copyP[cIndx];
    }

    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
      this->R_f[sIndx][fIndx] = copyrface[fIndx];
      this->RU[sIndx][fIndx] = copyru[fIndx];
    }
  }

  double* copyE = m->getE();
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    this->E[fIndx] = copyE[fIndx];
  }
}

///////////////////////////////////////////////////////////////////////
//
// Copy ParticleSolution to this MomentSolution
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::copyParticleSolution(ParticleSolution* p)
{
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {

    // Copy other parts of solution
    double* copyrho = p->getR(sIndx);
    double* copyru = p->getRU(sIndx);
    double* copys2 = p->getS2(sIndx);

    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
      this->R[sIndx][cIndx] = copyrho[cIndx];
      this->S2[sIndx][cIndx] = copys2[cIndx];
    }
    for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
      this->RU[sIndx][fIndx] = copyru[fIndx];
  }
}

///////////////////////////////////////////////////////////////////////
//
// Copy given MomentSolution into this
// Needed for time step iteration
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::copyMomentSolution(MomentSolution* m)
{
  double* copyE = m->getE();
  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++)
    this->E[fIndx] = copyE[fIndx];
}

///////////////////////////////////////////////////////////////////////
//
// Reduce Moment data from all nodes to node 0
//
///////////////////////////////////////////////////////////////////////

void MomentSolution::ReduceMoment()
{
//	const MPI::Datatype datatype(MPI_DOUBLE_INT);
//	const MPI::Op op(MPI_SUM);
//	MomentSolution* recBuf = new MomentSolution(this->data);
//
//	recBuf->copyMomentSolution(this);
//
//
//	  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
//
////		  MPI::Intracomm::Reduce(this->getR(sIndx), recBuf->getR(sIndx),
////				  this->data->nx,  datatype, op,0);
//#ifdef MPI_INC
//		  MPI_Reduce(this->getR(sIndx), recBuf->getR(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//		  MPI_Reduce(this->getRU(sIndx), recBuf->getRU(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//		  MPI_Reduce(this->getS2(sIndx), recBuf->getS2(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//		  MPI_Reduce(this->getS3(sIndx), recBuf->getS3(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//
//		  MPI_Reduce(this->getU(sIndx), recBuf->getU(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//		  MPI_Reduce(this->getT(sIndx), recBuf->getT(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//		  MPI_Reduce(this->getP(sIndx), recBuf->getP(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//
//		  MPI_Reduce(this->getETh(sIndx), recBuf->getETh(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//		  MPI_Reduce(this->getETot(sIndx), recBuf->getETot(sIndx),
//				  this->data->nx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//
//		  MPI_Reduce(this->getRFace(sIndx), recBuf->getRFace(sIndx),
//				  this->data->nfx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//		  MPI_Reduce(this->getRUCell(sIndx), recBuf->getRUCell(sIndx),
//				  this->data->nfx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//#endif
//
//	  }
//
//
//		  this->copyAllMoment(recBuf);


	  //delete [] recBuf;

}

void MomentSolution::BCastMoment()
{
//#ifdef MPI_INC
//	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_Bcast(this->getE(), this->data->nfx, MPI_DOUBLE,
//			0, MPI_COMM_WORLD);
//#endif
}

























