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

#include "PushRayTrace.h"
#include "PlasmaData.h"
#include "PlasmaUtility.h"
#include "ParticleSolution.h"
#include "MomentSolution.h"
#include "stdio.h"

/////////////////////////////////////////////////////////////////////////
// 
// Constructor
//
/////////////////////////////////////////////////////////////////////////

PushRayTrace::PushRayTrace(PlasmaData* data) : ParticlePush(data)
{
}

PushRayTrace::~PushRayTrace()
{
}

///////////////////////////////////////////////////////////////////////
//
// Particle tracking based on the paper by
// Graham B. Macpherson, Niklas Nordin and Henry G. Weller
// "Particle tracking in unstructured, arbitrary polyhedral meshes for use in
// CFD and molecular dynamics", Commun. Numer. Meth. Engng (2008)
// again, similar to the cell_info_calc.m, the implementation is a first cut 
// naive implementation and there's a lot of inefficient implementation.
// There's a lot that can be done here to optimize the particle tracking.
//
///////////////////////////////////////////////////////////////////////

void PushRayTrace::push(
                           MomentSolution* curMoment,       // unchanged
                           MomentSolution* oldMoment,       // unchanged
                           ParticleSolution* curParticle,   // updated
                           ParticleSolution* oldParticle)   // unchanged
{
  cout << "PushRayTrace::push()" << endl;
  double dx     = this->data->dx;
  double dx_recip = 1.0 / this->data->dx;
  double two_recip = 1.0 / 2.0;
  double tol_picard = 1.0e-10;
  double dx2_recip = 1.0 / (this->data->dx * 2.0);
  double dxd2 = this->data->dx / 2.0;
  double s = 1.0 / this->data->dx;
  double dt_res_tol = 0.0;

  cell_r = this->data->cell_r;
  cell_l = this->data->cell_l;
  face_c_r = this->data->face_c_r;
  face_c_l = this->data->face_c_l;
  cell_f_c = this->data->cell_f_c;
  c_f_id = this->data->c_f_id;
  c_f_n = this->data->c_f_n;
  
  double* E = new double[this->data->nfx];
  double* E_old = new double[this->data->nfx];
  double* E_half = new double[this->data->nfx];

  // Copy the electric field from the solver because it gets smoothed
  // Calculate the average E face value between time steps
  double* E_from_moment = curMoment->getE();
  double* E_old_from_moment = oldMoment->getE();

  for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
    E[fIndx] = E_from_moment[fIndx];
    E_old[fIndx] = E_old_from_moment[fIndx];
    E_half[fIndx] = (E[fIndx] + E_old[fIndx]) / 2.0;
  }

  // Smooth electric field
  for (int fil = 0; fil < this->data->fil_num; fil++) {
    if (this->data->temp_order == ORDER_1_BACK_EULER)
      smoothFaceFilter(this->data->nfx, E);
    else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON)
      smoothFaceFilter(this->data->nfx, E_half);
  }

  // Full and half moment solutions
  MomentSolution* full = new MomentSolution(this->data);
  MomentSolution* half = new MomentSolution(this->data);

  /////////////////////////////////////////////////////////////////
  //
  // Loop through all species
  //
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) {
    m = this->data->m[sIndx];            // mass of species
    q = this->data->q[sIndx];            // charge of species
    qm = q / m;                          // charge to mass ratio
    double mpp = this->data->mpp[sIndx]; // weight of species
    double dxmpp = mpp / this->data->dx;
    int numberOfParticles = this->data->NP0_tot[sIndx];
    // Old time particle information
    this->x_old = oldParticle->getX(sIndx);
    this->v_old = oldParticle->getV(sIndx);

    // New time particle information
    this->x_p = curParticle->getX(sIndx);
    this->v_p = curParticle->getV(sIndx);

    /////////////////////////////////////////////////////////////////
    //
    // Loop through all particles for the species
    //
    int tot_count = 0;           // total operation counter
    int tot_picard = 0;          // total picard operation counter
    int tot_sc = 0;              // total sub-cycling counter
    int tot_rt = 0;              // total ray tracing counter
      
    for (int pIndx = 0; pIndx < numberOfParticles; pIndx++) {
      x_sub_f = x_old[pIndx];
      v_sub_f = v_old[pIndx];

      cl_sub_f = (int) (x_sub_f / this->data->dx);
      if (cl_sub_f == this->data->nx)
        cl_sub_f = this->data->nx - 1;

      /////////////////////////////////////////////////////////////////
      //
      // Particle subcycling loop
      //
      dt_res = this->data->dt;         // residual for subcycling
      flag_sc = 0;

      int sc_counter = 0;              // sub-cycle counter
      while (flag_sc == 0) {
        if (sc_counter == 1000) {
          cout << "more than 1000 subcycling was required for particle " << pIndx << ", terminating subcycling for the particle." << endl;
            cout << "The x_sub_f = " << x_sub_f << " and x_sub_half = " << x_sub_half << endl;
            cout << "The last subcycling time step size was, dtau = " << dtau << endl;
          flag_sc = 1;
        }
        // Check if total time of particle is within tolerance
        tot_sc++;
        sc_counter++;

        // Set old time sub-cycle to previous sub-cycle final values
        x_sub_old = x_sub_f;
        v_sub_old = v_sub_f;
        cl_sub_old = cl_sub_f;

        // Set half time sub-cycle to previous sub-cycle final values
        x_sub_half = x_sub_f;
        v_sub_half = v_sub_f;
        cl_sub_half = cl_sub_f;

        // Electric field interpolation is hard coded for CPU speed
        // Uses parameters but does not change them
        if (this->data->temp_order == ORDER_1_BACK_EULER) {
          dtau = DtauEstimation(E, x_sub_old, v_sub_old, cl_sub_old);
//          dtau = EFieldInterpolation(E, x_sub_old, v_sub_old, cl_sub_old);
        }
        else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
          dtau = DtauEstimation(E_half, x_sub_old, v_sub_old, cl_sub_old);
//          dtau = EFieldInterpolation(E_half, x_sub_old, v_sub_old, cl_sub_old);
        }
        // Initialize the kth Picard solution
        x_k = x_sub_old;
        v_k = v_sub_old;
        flag_conv = 0;      // Picard convergence flag
        int k = 0;          // Picard iteration index
        flag_fx = 0;        // Particle crossed a face

        /////////////////////////////////////////////////////////////////
        //
        // Picard Crank-Nicholson convergence loop
        // 
        double rel_v = 1.0;
        double rel_x = 1.0;

        while (flag_conv == 0 && k < 10) {
          tot_picard++;
          k++;
          x0 = x_sub_old;
          cl = cl_sub_old;
          x_sub_half = x_sub_old;

          // Total displacement distance
//          dx_p;
          if (this->data->temp_order == ORDER_1_BACK_EULER)
            dx_p = dtau * v_sub_f;
          else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON)
            dx_p = dtau * v_sub_half;
          double dx_half_p = dx_p * two_recip;
          /////////////////////////////////////////////////////////////////
          //
          // Particle tracking algorithm (ray trace)
          //
          flag_hs = 0;                   // half step position flag
          part_cf_xf0 = 0;               // particle crosses a face

          int rt_count = particleRayTrace(dx_half_p);

          tot_count += rt_count;
          tot_rt += rt_count;
          //
          /////////////////////////////////////////////////////////////////

          // End of ray tracing, calculate convergence of particle trajectory
          // E field interpolation is hard coded for CPU speed
          double* E_E;
          double x_E;
          int cl_E;
          if (this->data->temp_order == ORDER_1_BACK_EULER) {
            E_E = E;
            x_E = x_sub_f;
            cl_E = cl_sub_f;
          }
          else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
            E_E = E_half;
            x_E = x_sub_half;
            cl_E = cl_sub_half;
          }

          // Set parameters for electric field interpolation
          double ddx = x_E - this->data->xpos_face[cl_E];
          double dEdx = (E_E[cl_E + 1] - E_E[cl_E]) * dx_recip;
          E_p = E_E[cl_E] + dEdx * ddx;

          // Particle reversed direction at a face
          if (part_cf_xf0 == 1) {
            // Final position for sub step at the face where particle reversed
            x_sub_f = x_px;
            x_sub_half = x_px;
            v_sub_f = v_px;
            v_sub_half = v_px;
            cl_sub_f = cl_px;
            cl_sub_half = cl_px;
            dtau = 0.0; // subcycle time step size since it reversed
          }
          // Particle passed through surface
          else {
            v_sub_f = v_sub_old + dtau * qm * E_p; // new time velocity
            cl_sub_f = cl;
            v_sub_half = (v_sub_f + v_sub_old) * two_recip; //half time velocity
          }

          // Relative difference between successive Picard iterations
          rel_x = (x_sub_f - x_k) / x_sub_f;
          rel_v = (v_sub_f - v_k) / v_sub_f;

          // Set position and velocity for next Picard
          x_k = x_sub_f;
          v_k = v_sub_f;

          if (fabs(rel_v) <= tol_picard && fabs(rel_x) <= tol_picard)
            flag_conv = 1;
        }
        // Picard convergence loop
        /////////////////////////////////////////////////////////////////

        // Particle crossed a surface
        if (flag_fx == 1) {
          if (part_cf_xf0 == 0) {

            // Particle crossed face, adjust dtau to stop particle at face
            // Assign final particle position to be the face at which the
            // particle crossed, approached from the direction of vel vector
            x_sub_f = x_xf;
            x_sub_half = (x_sub_f + x_sub_old) * two_recip;
            cl_sub_half = cl_sub_old;

// PKF This should be (E,x_sub_f, cl_sub_f) but it doesn't work
            if (this->data->temp_order == ORDER_1_BACK_EULER) {
              dtau = DtauCalculationCross(E, x_sub_half, cl_sub_half);
//              dtau = EFieldInterpolationCross(E, x_sub_half, cl_sub_half);
            }
            else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
              dtau = DtauCalculationCross(E_half, x_sub_half, cl_sub_half);
//              dtau = EFieldInterpolationCross(E_half, x_sub_half, cl_sub_half);
            }

            // With the back calculated subcycle time size
            // back calculate the final velocity and half time velocity
            v_sub_f = v_sub_old + dtau * qm * E_p;
            v_sub_half = (v_sub_f + v_sub_old) * two_recip;
  
            // Now update the final particle position to be the face location
            // approached from the opposite direction of the particle
            // velocity vector.  This is done due to the periodic boundary
            // condition which technically allows 2 positions for the same face.
            x_sub_f = xf_real;
            cl_sub_f = cl_xf;
          }
          // Particle reversed direction at cell face at current subcycle
          else if (part_cf_xf0 == 1) {
            x_sub_f = x_px;
            x_sub_half = x_px;
            v_sub_f = v_px;
            v_sub_half = v_px;
            cl_sub_f = cl_px;
            cl_sub_half = cl_px;
            dtau = 0.0;
          }
        }

        // Tally half time moment quantities (s1_p_half and s2_p_half)
        half->tallyMomentCalculation(sIndx, x_sub_half, v_sub_half, dtau);

        // Update sub step residual to see if loop can terminate
        dt_res -= dtau;
        if (dt_res <= dt_res_tol)
          flag_sc = 1;
      }
      // Subcycling loop
      /////////////////////////////////////////////////////////////////

      // Tally full time moment quantities (r_p and s1_p)
      full->tallyMomentCalculation(sIndx, x_sub_f, v_sub_f, 1.0);

      // Store the new particle position and velocity
      x_p[pIndx] = x_sub_f;
      v_p[pIndx] = v_sub_f;
    }
    // Every particle loop
    /////////////////////////////////////////////////////////////////

    // Scale the moment quantities
    full->scaleMoment(sIndx, dxmpp);
    half->scaleMoment(sIndx, dxmpp / this->data->dt);

    // Smooth moment quantities
    full->smoothMoment(sIndx);
    half->smoothMoment(sIndx);

    cout << "The average ray tracing operation per particle for species "
         << sIndx << " is: " 
         << (double) tot_count / (double) numberOfParticles << endl;
    cout << "The average sub-cycling operation per particle for species "
         << sIndx << " is: " 
         << (double) tot_sc / (double) numberOfParticles << endl;
    cout << "The average Picard ops per particle per sub-cycle for species "
         << sIndx << " is: " 
         << (double) tot_picard / (double) tot_sc 
         << endl;
    cout << "The average ray-tracing ops per particle per sub-cycle "
         << "per picard iteration for species "
         << sIndx << " is: " 
         << (double) tot_rt / (double) tot_sc 
         << endl;

    // Calculate charge conservation properties
    double* r_old = oldParticle->getR(sIndx);
    double cc_rms = 0.0;
    double dt_recip = 1.0 / this->data->dt;

    double* r = full->getR(sIndx);
    double* s1 = full->getRU(sIndx);
    double* s1_half = half->getRU(sIndx);

    if (this->data->temp_order == ORDER_1_BACK_EULER) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        double cc = dt_recip * (r[cIndx] - r_old[cIndx]) +
                    dx_recip * (s1[cIndx + 1] - s1[cIndx]);
        cc_rms += cc * cc;
      }
    } else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        double cc = dt_recip * (r[cIndx] - r_old[cIndx]) +
                    dx_recip * (s1_half[cIndx + 1] - s1_half[cIndx]);
        cc_rms += cc * cc;
      }
    }
/*      cout << "r for species " << sIndx << ": \n";
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++) 
      {
          cout << r[cIndx] << "\n";
      }
      cout << "s1_half for species " << sIndx << ": \n";
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++) 
      {
          cout << s1_half[cIndx] << "\n";
      }*/
    // Calculate root mean square of conservation equation
    cc_rms = sqrt(cc_rms / this->data->nx);
    cout << endl <<  "***   CC_RMS " << cc_rms << endl << endl;

    // Store full and half time moment quantities into current particle solution
    curParticle->storeR(sIndx, full->getR(sIndx));
    curParticle->storeS2(sIndx, full->getS2(sIndx));
    curParticle->storeRU(sIndx, full->getRU(sIndx));

    curParticle->storeS2Half(sIndx, half->getS2(sIndx));
    curParticle->storeRUHalf(sIndx, half->getRU(sIndx));

    curParticle->storeCC_RMS(sIndx, cc_rms);

  }
  // Every species loop 
  /////////////////////////////////////////////////////////////////

  delete full;
  delete half;
  delete [] E;
  delete [] E_old;
  delete [] E_half;
}

///////////////////////////////////////////////////////////////////////////
//
//  Dtau estimation based on old subcycle time particle velocity and cell size
//
///////////////////////////////////////////////////////////////////////////
double PushRayTrace::DtauEstimation(
                                         double* E_E,
                                         double x_E,
                                         double v_E,
                                         int cl_E)
/*double PushRayTrace::EFieldInterpolation(
                             double* E_E,
                             double x_E,
                             double v_E,
                             int cl_E)*/
{
  double dx_recip = 1.0 / this->data->dx;

  // Set parameters for electric field interpolation
  double ddx = x_E - this->data->xpos_face[cl_E];
  double dEdx = (E_E[cl_E + 1] - E_E[cl_E]) * dx_recip;

  // Electric field interpolation is hard coded for CPU speed
  // Set parameters for electric field interpolation
  double E_p0 = E_E[cl_E] + dEdx * ddx;
  double E_p0_sq = E_p0 * E_p0;
  double qm_sq = qm * qm;
  double two_recip = 1.0 / 2.0;
  double eps_r = this->data->eps_r;
  double eps_a = this->data->eps_a;

  // Calculate coefficients for quadratic step for upper limit of time
  // step size.  This will try to calculate dtau from
  // A_ul * dtau^2 + B_ul * dtau + C_ul
  double A_ul = sqrt(qm_sq * E_p0_sq +
                     qm_sq * dEdx * dEdx * v_E * v_E) * two_recip;
  double A_ul_recip = 1.0 / A_ul;
  double B_ul = -eps_r * sqrt(v_E * v_E + E_p0_sq * qm_sq);
  double C_ul = -eps_a;

  // Calculate the upper limit time step size
  double pm_part = sqrt(B_ul * B_ul - 4.0 * A_ul * C_ul);
  double dtau_ul_p = two_recip * A_ul_recip * (-B_ul + pm_part);
  double dtau_ul_m;
  double dtau_ul;
  double dtau;
  if (dtau_ul_p < 0) {
    dtau_ul_m = two_recip * A_ul_recip * (-B_ul - pm_part);
    dtau_ul = dtau_ul_m;
  } else {
    dtau_ul_m = two_recip * A_ul_recip * (-B_ul - pm_part);
    if (dtau_ul_p < dtau_ul_m && dtau_ul_m > 0.0)
      dtau_ul = dtau_ul_m;
    else
      dtau_ul = dtau_ul_p;
/* PKF changed in Will's code but it causes a hang
    if (dtau_ul_p < dtau_ul_m && dtau_ul_m > 0.0)
      dtau_ul = dtau_ul_p;
    else
      dtau_ul = dtau_ul_m;
*/
  }

  // Check if the subcycle time step is larger than the time residual
  if (dtau_ul <= dt_res)
    dtau = dtau_ul;
  else
    dtau = dt_res;

  double dtau_se = this->data->dx / fabs(v_sub_old);
  if (dtau >= dtau_se)
    dtau = dtau_se;

  return dtau;
}

///////////////////////////////////////////////////////////////////////////
//
//  Dtau calculation based on cell face crossing
//
///////////////////////////////////////////////////////////////////////////
double PushRayTrace::DtauCalculationCross(
                                              double* E_E,
                                              double x_E,
                                              int cl_E)
/*double PushRayTrace::EFieldInterpolationCross(
                             double* E_E,
                             double x_E,
                             int cl_E)*/
{
  double dx_recip = 1.0 / this->data->dx;
  double two_recip = 1.0 / 2.0;

  // Set parameters for electric field interpolation
  double ddx = x_E - this->data->xpos_face[cl_E];
  double dEdx = (E_E[cl_E + 1] - E_E[cl_E]) * dx_recip;
  E_p = E_E[cl_E] + dEdx * ddx;
  
  // Back calculate the subcycle time step size based on the face
  // which we know since we stopped particle at the face
  double A_tau = -two_recip * qm * E_p;
  double B_tau = -v_sub_old;
  double C_tau = x_sub_f - x_sub_old;
  double A_tau_recip = 1.0 / A_tau;

  double dtau_pm = sqrt(B_tau * B_tau - 4.0 * A_tau * C_tau);
  double dtau_p = two_recip * A_tau_recip * (-B_tau + dtau_pm);
  double dtau_m = two_recip * A_tau_recip * (-B_tau - dtau_pm);

  double v_sub_f_p = v_sub_old + dtau_p * qm * E_p;
  double dtau;
  
  // If positive, correct sign
  if (c_f_n_xf * v_sub_f_p >= 0.0) {
    if (dtau_p <= dtau_m && dtau_m > 0.0)
      dtau = dtau_p;
    else if (dtau_p > dtau_m && dtau_m > 0.0)
      dtau = dtau_m;
    else if (dtau_m < 0.0)
      dtau = dtau_p;
  } else {
    dtau = dtau_m;
  }
  return dtau;
}

///////////////////////////////////////////////////////////////////////////
//
// Loop through all faces in the cell to check surface crossing
// Choose the face which is closest to starting point
// Hard coded for a 1D problem
//
///////////////////////////////////////////////////////////////////////////

int PushRayTrace::particleRayTrace(double dx_half_p)
{
  int lenf = 2;              // particle can only cross 2 faces in 1D
  double two_recip = 1.0 / 2.0;
//  double dx_half_p = dx_p * two_recip;

  int rt_count = 0;
  flag_rt = 0;
//  while (flag_rt == 0 && rt_count < 100) {
  while (flag_rt == 0) {
    rt_count++;

    // Loop through all faces in the cell to check surface crossing
    int cf_i = 0;                // index of face crossed by particle
    double lambda = ((cell_f_c[c_f_id[0][cl]] - x0) * 
                     c_f_n[0][cl]) / (dx_p * c_f_n[0][cl]);

    for (int fIndx = 1; fIndx < lenf; fIndx++) {
      double lambda_next = ((cell_f_c[c_f_id[fIndx][cl]] - x0) * 
                             c_f_n[fIndx][cl]) / (dx_p * c_f_n[fIndx][cl]);
      if (lambda_next >= 0.0) {
        if (lambda_next > lambda || lambda < 0.0) {
          lambda = lambda_next;
          cf_i = fIndx;
        }
      }
    }
  
    // Find where particle crossed
    int cross_flag = 0;          // type of face cross
    double lambda_a = 0.0;       // distance to the face normalized

    // Particle crossed face
    if (lambda > 0.0 && lambda < 1.0) {
      cross_flag = CROSSED;
      lambda_a = lambda;
    }
    // Particle is trapped
    else if (lambda >= 1.0) {
      cross_flag = TRAPPED;
    }
    // Particle is reversing into previous cell
    else if (lambda == 0.0) {
      lambda_a = 0.0;
      if (this->data->temp_order == ORDER_1_BACK_EULER) {
        cross_flag = particleReverse(cf_i, v_sub_f);
      }
      else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) {
        cross_flag = particleReverse(cf_i, v_sub_half);
      }
    }

    // Particle crossed cell face
    if (cross_flag == CROSSED) {
      double lambdx_p = lambda_a * dx_p;
      dx_p -= lambdx_p;

      // Check if particle traveled half of total displacement
      if (fabs(dx_p) <= fabs(dx_half_p) && flag_hs == 0) {
        // Half time position of particle
        x_sub_half = x0 + (dx_p - dx_half_p + lambdx_p);
        cl_sub_half = cl;
        flag_hs = 1;
      }

      // Determine which cell particle streams into
      // Normal vector of surface is positive, stream to right cell
      if (c_f_n[cf_i][cl] == 1 ) {

        // If first surface that particle crosses
        if (flag_fx == 0) {
          flag_fx = 1;                    // flag for face cross
          cl_xf = face_c_r[cell_r[cl]];   // cell particle crosses into
          c_f_n_xf = c_f_n[cf_i][cl];     // normal vector of that cell
          x_xf = cell_f_c[cell_r[cl]];    // face pos approached from left side
          xf_real = cell_f_c[cell_l[cl_xf]]; // final position of particle
        }
        // Update cell particle is in, set position to face
        cl = face_c_r[cell_r[cl]];
        x0 = cell_f_c[cell_l[cl]];
      }
      // Normal vector of surface is negative, stream to left cell
      else {
        // If first surface that particle crosses
        if (flag_fx == 0) {
          flag_fx = 1;                    // flag for face cross
          cl_xf = face_c_l[cell_l[cl]];   // cell particle crosses into
          c_f_n_xf = c_f_n[cf_i][cl];     // normal vector of that cell
          x_xf = cell_f_c[cell_l[cl]];    // face pos approached from right side
          xf_real = cell_f_c[cell_r[cl_xf]]; // final position of particle
        }
        // Update cell particle is in, set position to face
        cl = face_c_l[cell_l[cl]];
        x0 = cell_f_c[cell_r[cl]];
      }
    }
    // Particle is trapped in cell
    else if (cross_flag == TRAPPED) {
      x_sub_f = x0 + dx_p;                // final position of particle
      cl_sub_f = cl;                      // cell it is in

      // If half distance hasn't been occupied yet
      if (flag_hs == 0) {
        x_sub_half = x0 + (dx_p - dx_half_p);
        cl_sub_half = cl;
        flag_hs = 1;
      }
      dx_p = 0.0;
      flag_rt = 1;
    }
    // Particle reverses direction at face
    else {
      x_sub_f = x_px;                     // final time same starting position
      x_sub_half = x_px;                  // half time same starting position
      v_sub_f = v_px;                     // final velocity
      v_sub_half = (v_sub_old + v_sub_f); // time centered velocity
      dtau = 0.0;                         // particle does not move
      flag_fx = 1;                        // particle crossed a face
    }
  }
  return rt_count;                        // Number of ray trace iterations
}

///////////////////////////////////////////////////////////////////////////
//
// Particle reverse
//
///////////////////////////////////////////////////////////////////////////

int PushRayTrace::particleReverse(
                  int cf_i,
                  double v0)
{
  int cross_flag = BAD_CROSS;

  if (c_f_n[cf_i][cl] * v0 >= 0.0) {
    part_cf_xf0 = 1;
    dtau = 0.0;
    flag_conv = 1;
    flag_rt = 1;
    flag_hs = 1;
  
    // Determine which cell particle streams into
    // Normal vector of surface positive, stream to right cell
    if (c_f_n[cf_i][cl] == 1.0) {
      cl = face_c_r[cell_r[cl]];  // cell particle traverse to
      x0 = cell_f_c[cell_l[cl]];  // assign particle pos to face
      x_px = x0;
      v_px = v0;
      cl_px = cl;
    }
    // Normal vector of surface negative, stream to left cell
    else {
      cl = face_c_l[cell_l[cl]]; //cell particle traversed to
      x0 = cell_f_c[cell_r[cl]]; // assign particle pos to face
      x_px = x0;
      v_px = v0;
      cl_px = cl;
    }
    cross_flag = REVERSED;
  }
  return cross_flag;
}
