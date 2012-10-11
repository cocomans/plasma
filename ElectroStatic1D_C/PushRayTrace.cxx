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
    double x_E;

    // Electric field interpolation is hard coded for CPU speed
    // Set parameters for electric field interpolation
    double eps_r = this->data->eps_r;
    double eps_a = this->data->eps_a;
    
    double E_p0;
    double E_p0_sq;
    double qm_sq;
    
    
    double A_ul;
    double A_ul_recip;
    double B_ul;
    double C_ul;
    
    double A_tau;
    double B_tau;
    double C_tau;
    double A_tau_recip;
    double dtau_pm;
    double dtau_p;
    double dtau_m;
    double v_sub_f_p;
    double x_px, v_px, x_xf, xf_real, x0;
    double dx_p, dx_half_p;
    
    int    cl_px, cl_xf, cl, cl_sub_half, c_f_n_xf;
    int    flag_fx, part_cf_xf0, k, flag_hs, flag_rt, cross_flag;
    double lambda_a, lambda, lambda_next, lambda_first, lambdx_p;
    int    lam_i;
    int    rt_count, tot_count_sc, tot_rt_sc, cf_i, tot_count, tot_picard, tot_rt, tot_sc;
    
    double pm_part;
    double dtau_ul_p;
    double dtau_ul_m;
    double dtau_ul;
    double dtau;
    
    double x_e;
    double ddx;
    double dEdx;
    
    int cl_E;
    
    cell_r = this->data->cell_r;
    cell_l = this->data->cell_l;
    face_c_r = this->data->face_c_r;
    face_c_l = this->data->face_c_l;
    cell_f_c = this->data->cell_f_c;
    c_f_n = this->data->c_f_n;
    c_f_id = this->data->c_f_id;
/*
    cout << "cell_l, cell_r" << endl;
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        cout << cell_l[cIndx] << ", " << cell_r[cIndx] << endl;
    }
    cout << "face_c_l, face_c_r" << endl;
    for (int cIndx = 0; cIndx < this->data->nfx; cIndx++) {
        cout << face_c_l[cIndx] << ", " << face_c_r[cIndx] << endl;
    }
    cout << "cell_f_c[cIndx]" << endl;
    for (int cIndx = 0; cIndx < this->data->nfx; cIndx++) {
        cout << cell_f_c[cIndx] << endl;
    }
    cout << "c_f_id[cIndx][0], c_f_id[cInx][1]" << endl;
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        cout << c_f_id[0][cIndx] << " , " << c_f_id[1][cIndx] << endl;
    }
    cout << "cell_f_n[cIndx][0], cell_f_n[cInx][1]" << endl;
    for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        cout << c_f_n[0][cIndx] << " , " << c_f_n[1][cIndx] << endl;
    }    
*/
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
  for (int sIndx = 0; sIndx < this->data->p_size; sIndx++) 
  {
    m = this->data->m[sIndx];            // mass of species
    q = this->data->q[sIndx];            // charge of species
    qm = q / m;                          // charge to mass ratio
    double mpp = this->data->mpp[sIndx]; // weight of species
    double dxmpp = mpp / this->data->dx;
    int numberOfParticles = this->data->NP0_tot[sIndx];
    // Old time particle information
    this->x_old = oldParticle->getX(sIndx);
    this->v_old = oldParticle->getV(sIndx);
      double* E_m = new double[this->data->nfx];
      if (this->data->temp_order == ORDER_1_BACK_EULER) {
          for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
              E_m[fIndx] = E[fIndx];
          }
      }
      else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON){
          for (int fIndx = 0; fIndx < this->data->nfx; fIndx++) {
              E_m[fIndx] = E_half[fIndx];
          }
      }

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
      
    for (int pIndx = 0; pIndx < numberOfParticles; pIndx++) 
    {
      x_sub_f = x_old[pIndx];
      v_sub_f = v_old[pIndx];

      cl_sub_f = (int) (x_sub_f / this->data->dx);
      if (cl_sub_f == this->data->nx)
        cl_sub_f = this->data->nx - 1;
/*        cl_sub_f = ceil(x_sub_f/this->data->dx);
        cout << cl_sub_f << endl;
*/
      /////////////////////////////////////////////////////////////////
      //
      // Particle subcycling loop
      //
      dt_res = this->data->dt;         // residual for subcycling
      flag_sc = 0;

      int sc_counter = 0;              // sub-cycle counter
      double  absv_old  = fabs(v_sub_f);
      double  dtau_se   = dx/absv_old;
      while (flag_sc == 0) 
      {

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

          //Electric field interpolation is hard-coded for now 
          //set parameters for electric field interpolation
          double x_E         = x_sub_old;                //set position for particle for E-field interpolation
          int cl_E        = cl_sub_old;               //set cell index for E-field interpolation for particle
          double ddx         = x_E - this->data->xpos_face[cl_E];    //calculate the distance from the left face of the cell in which particle lives in
          double dEdx        = (E_m[cl_E+1] - E_m[cl_E])*dx_recip;//calculate gradiet in electric field for the cell in which particle lives in
          double E_p0        = E_m[cl_E] + dEdx*ddx;     //calculation of electric field for particle
          double E_p0_sq     = E_p0*E_p0;                //calculate square of initial E-field for particle
          //calculate the coefficients for the quadratic step for upper limit of dtau for time estimation
          //This part will try to calculate dtau from:
          //A_ul*dtau^2 + B_ul*dtau + C_ul
          A_ul        = sqrt(qm_sq*E_p0_sq + qm_sq*dEdx*dEdx*v_sub_old*v_sub_old)*two_recip;  //Coefficient for dtau calculation for second order term
          A_ul_recip  = 1.0/A_ul;                                                             //Reciprical for A_ul
          B_ul        = -eps_r*sqrt(v_sub_old*v_sub_old + E_p0_sq*qm_sq);                     //Coefficient for dtau calculation for first order term
          C_ul        = -eps_a;                                                               //Coefficient for dtau calculation for 0th order term
          //calculate the upper limit time step size
          pm_part     = sqrt(B_ul*B_ul - 4.0*A_ul*C_ul);            
          dtau_ul_p   = two_recip*A_ul_recip*(-B_ul + pm_part);
          if (dtau_ul_p < 0.0)                                        //if upper limit dtau for position root is less than 0.0, set the dtau_ul to the negative root case
          {
              dtau_ul_m   = two_recip*A_ul_recip*(-B_ul - pm_part);
              dtau_ul     = dtau_ul_m;
          }
          else 
          {
              dtau_ul_m   = two_recip*A_ul_recip*(-B_ul - pm_part);
              if (dtau_ul_p < dtau_ul_m && dtau_ul_m > 0.0)           //check if the upper limit positive root is less than the negative root but the negative root is positive
              {
                  //    dtau_ul = dtau_ul_p;                                //set the dtau_ul to the positive root
                  dtau_ul = dtau_ul_m;                                //set the dtau_ul to the negative root
              }
              else
              {
                  dtau_ul = dtau_ul_p;                                //set the dtau_ul to the positive root
                  //    dtau_ul = dtau_ul_m;                                //set the dtau_ul to the negative root
              }
          }
          //check if the subcycle time-step is larger than the time residual
          if (dtau_ul <= dt_res) 
          {
              dtau = dtau_ul;
          }
          else 
          {
              dtau = dt_res;
          }
          //Now get a SMART estimate of dtau based on global time step size and initial sub-cycle particle velocity 
          //Since we know that the particle will HAVE to stop at cell faces, there's no meaning to take a large 
          //sub-cycle time step size that makes particle traverse multiple cell
          //absv_old = fabs(v_sub_old);
          //dtau_se  = dx/absv_old;
          //check if the estimated upper limit time is larger than the optical depth of particle in the cell
          if (dtau >= dtau_se) 
          {
              dtau = dtau_se;
          }
          if (sc_counter >= 1000) //i.e. if there's more than 100 sub-cycling, somethings wrong, simply recompute dtau to be dt_res
          {
              dtau = dt_res;
              cout << "More than 1000 sub-cycling was required for " << pIndx << "particle, terminating sub-cycling." << endl;
          }                

        // Initialize the kth Picard solution
        x_k = x_sub_old;
        v_k = v_sub_old;
        flag_conv = 0;      // Picard convergence flag
        k = 0;          // Picard iteration index
        flag_fx = 0;        // Particle crossed a face

        /////////////////////////////////////////////////////////////////
        //
        // Picard Crank-Nicholson convergence loop
        // 
        double rel_v = 1.0;
        double rel_x = 1.0;

        while (flag_conv == 0) 
        {
          tot_picard++;
          k++;
          x0 = x_sub_old;
          cl = cl_sub_old;
          x_sub_half = x_sub_old;

          // Total displacement distance
          dx_p;
          if (this->data->temp_order == ORDER_1_BACK_EULER)
          {
            dx_p = dtau * v_sub_f;
          }
          else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON)
          {
            dx_p = dtau * v_sub_half;
          }
          dx_half_p = dx_p*two_recip;
          /////////////////////////////////////////////////////////////////
          //
          // Particle tracking algorithm (ray trace)
          //
          flag_hs = 0;                   // half step position flag
          flag_rt = 0;
          part_cf_xf0 = 0;               // particle crosses a face
          rt_count  = 0;
          flag_fx   = 0;
          cross_flag = 0;
          part_cf_xf0 = 0;
          while (flag_rt == 0) 
          {
              tot_count_sc     = tot_count_sc + 1;    //increment total operation counter
              tot_rt_sc        = tot_rt_sc + 1;       //increment total ray tracing counter
              rt_count         = rt_count + 1;        //count the number of cells the particle traverses
              //lenf_a = 2;//hard coding the # of faces per cell
              //loop through all faces in the cell to check surface x-ing
              lam_i   = 0; //initialize the face index
              while (lam_i < 2) 
              {
                  if (lam_i == 0)     //first surface in cell to check
                  {
                      cf_i    = lam_i;                            
                      lambda = ((cell_f_c[c_f_id[lam_i][cl]]) - x0)*c_f_n[lam_i][cl]
                      /
                      ((dx_p)*c_f_n[lam_i][cl]); //calculate lambda from the reference paper in the top of the code
                      lambda_first = lambda;                      //lambda here just for debugging purposes
                  }
                  else if (lam_i > 0) //next surface in cell to check
                  {
                      lambda_next = ((cell_f_c[c_f_id[lam_i][cl]]) - x0)*c_f_n[lam_i][cl]
                      /
                      ((dx_p)*c_f_n[lam_i][cl]);
                      if (lambda_next >= 0.) 
                      {
                          if (lambda_next > lambda) 
                          {
                              lambda  = lambda_next;
                              cf_i    = lam_i;
                          }
                          else if (lambda < 0.) 
                          {
                              lambda  = lambda_next;
                              cf_i    = lam_i;
                          }
                      }
                  }
                  lam_i++;
              }
              if (lambda > 0.0 && lambda < 1.0)   //check if particle crossed a face
              {
                  cross_flag      = 1;
                  lambda_a        = lambda;
              }
              else if (lambda >= 1.0)             //check if particle is trapped
              {
                  cross_flag      = 2;
              }
              else if (lambda == 0.0)             //check if particle is reversing into previous cell
              {
                  lambda_a        = 0.0;
                  if (this->data->temp_order == 1) 
                  {
                      if (c_f_n[cf_i][cl]*v_sub_f >= 0.0) //if negative, particle crosses into adjacel cell
                      {
                          part_cf_xf0     = 1;                    //turn on flag to indicate the particle direction reversed at cell face
                          dtau            = 0.0;                  //set sub-cycle time step size to 0 since it didn't "move" in the sub-cycle, but just reversed directioin
                          flag_conv       = 1;                    //turn on flag to indicate the sub-cycle converged
                          flag_rt         = 1;                    //turn on ray-trace flag ot indicate the ray trace has been completed for the sub-cycle
                          flag_hs         = 1;                    //turn on the half sub-cycle time flag since the particle didn't move
                          //determine which cell particle streams into
                          if (c_f_n[cf_i][cl] == 1)       //if normal vector of surface is positive, streams into cell to the right
                          {   
                              cl      = face_c_r[cell_r[cl]];     //determine which cell particle traversed into
                              x0      = cell_f_c[cell_l[cl]];     //assign particle position to face
                              x_px    = x0;                       //store final sub-cycle position of particle to initial position since particle essentially just reversed direction at cell face
                              v_px    = v_sub_f;                  //however, store final sub-cycle particle velocity to the final velocity to indicate direction of particle
                              cl_px   = cl;                       //the final sub-cycle cell index for particle will just be the initial cell in which sub-cycle began.
                          }
                          else                                    //if normal vector of surface is negative, stream into the left cell
                          {
                              cl      = face_c_l[cell_l[cl]];     //determine which cell particle traversed into
                              x0      = cell_f_c[cell_r[cl]];     //assign particle position to face
                              x_px    = x0;                       //store final sub-cycle position of particle to initial position since particle essentially just reversed direction at cell face
                              v_px    = v_sub_f;                  //however, store final sub-cycle particle velocity to the final velocity to indicate direction of particle
                              cl_px   = cl;                       //the final sub-cycle cell index for particle will just be the initial cell in which sub-cycle began.
                          }
                          cross_flag  = 3;                        //set cross_flag to 3 to indicate that the particle reversed direction at cell face
                      }
                      else                                        //for debugging purposes
                      {
                          cout << "Something is wrong in ray tracing \n" << endl;
                          lambda_a = 0.;
                      }
                  }
                  else if (this->data->temp_order == 2) 
                  {
                      if (c_f_n[cf_i][cl]*v_sub_half >= 0.0) 
                      {
                          part_cf_xf0     = 1;                    //turn on flag to indicate the particle direction reversed at cell face
                          dtau            = 0.0;                  //set sub-cycle time step size to 0 since it didn't "move" in the sub-cycle, but just reversed directioin
                          flag_conv       = 1;                    //turn on flag to indicate the sub-cycle converged
                          flag_rt         = 1;                    //turn on ray-trace flag ot indicate the ray trace has been completed for the sub-cycle
                          flag_hs         = 1;                    //turn on the half sub-cycle time flag since the particle didn't move
                          //determine which cell particle streams into
                          if (c_f_n[cf_i][cl] == 1)       //if normal vector of surface is positive, streams into cell to the right
                          {   
                              cl      = face_c_r[cell_r[cl]];     //determine which cell particle traversed into
                              x0      = cell_f_c[cell_l[cl]];     //assign particle position to face
                              x_px    = x0;                       //store final sub-cycle position of particle to initial position since particle essentially just reversed direction at cell face
                              v_px    = v_sub_f;                  //however, store final sub-cycle particle velocity to the final velocity to indicate direction of particle
                              cl_px   = cl;                       //the final sub-cycle cell index for particle will just be the initial cell in which sub-cycle began.
                          }
                          else                                    //if normal vector of surface is negative, stream into the left cell
                          {
                              cl      = face_c_l[cell_l[cl]];     //determine which cell particle traversed into
                              x0      = cell_f_c[cell_r[cl]];     //assign particle position to face
                              x_px    = x0;                       //store final sub-cycle position of particle to initial position since particle essentially just reversed direction at cell face
                              v_px    = v_sub_f;                  //however, store final sub-cycle particle velocity to the final velocity to indicate direction of particle
                              cl_px   = cl;                       //the final sub-cycle cell index for particle will just be the initial cell in which sub-cycle began.
                          }
                          cross_flag  = 3;                        //set cross_flag to 3 to indicate that the particle reversed direction at cell face
                      }
                      else                                        //for debugging purposes
                      {
                          cout <<"Something is wrong in ray tracing\n" << endl;
                          lambda_a = 0.;
                      }
                  }
              }
              switch (cross_flag) 
              {
                  case 1:                                                     //if particle crosses surface
                      lambdx_p    = lambda_a*dx_p;                            //calculate fraction of total displacement for current sub-cycle at current cell
                      dx_p        = dx_p - lambdx_p;                          //update the remaining distance for particle to travel for current sub-cycle
                      if (fabs(dx_p) <= fabs(dx_half_p) && flag_hs == 0)            //check if the absolute value of dx_p is less than or equal to the absolute value of dx_half
                      {
                          x_sub_half  = x0 + (dx_p - dx_half_p + lambdx_p);   //set half sub-cycle time position
                          cl_sub_half = cl;                                   //set half sub-cycle time cell index
                          flag_hs     = 1;                                    //turn on that half distance have been reached
                      }
                      if (c_f_n[cf_i][cl] ==1)                        //check if particle crosses right surface of cell
                      {   
                          if (flag_fx == 0)                                   //if particle has not yet crossed a face yet within a sub-cycle
                          {
                              flag_fx     = 1;                                //turn on flag saying that the particle has already crossed a face
                              cl_xf       = face_c_r[cell_r[cl]];             //determine which cell is to the right (automatically takes into account of periodic condition
                              c_f_n_xf    = c_f_n[cf_i][cl];          //store the normal vector of face that particle crossed
                              x_xf        = cell_f_c[cell_r[cl]];             //store the particle position at face from the left limit
                              xf_real     = cell_f_c[cell_l[cl_xf]];          //store the "real" particle position "after" it crossed a face (i.e. particle position on face from right limit) this automatically takes care of the periodic condition
                          }
                          cl      = face_c_r[cell_r[cl]];                     //update the cell index of particle
                          x0      = cell_f_c[cell_l[cl]];                     //update the new "initial" sub-cycle particle position
                      }
                      else                                                    //check if particle crosses left surface of cell
                      {
                          if (flag_fx == 0) 
                          {
                              flag_fx     = 1;                                //turn on flag saying that the particle has already crossed a face
                              cl_xf       = face_c_l[cell_l[cl]];             //determine which cell is to the right (automatically takes into account of periodic condition
                              c_f_n_xf    = c_f_n[cf_i][cl];          //store the normal vector of face that particle crossed
                              x_xf        = cell_f_c[cell_l[cl]];             //store the particle position at face from the left limit
                              xf_real     = cell_f_c[cell_r[cl_xf]];          //store the "real" particle position "after" it crossed a face (i.e. particle position on face from right
                          }
                          cl      = face_c_l[cell_l[cl]];                     //update the cell index of particle
                          x0      = cell_f_c[cell_r[cl]];                     //update the new "initial" sub-cycle particle position
                      }
                      break;                            
                  case 2:                                                     //if particle is trapped in cell
                      x_sub_f     = x0 + dx_p;                                //set final position for sub-cycle
                      cl_sub_f    = cl;                                       //set final cell index for sub-cycle
                      if (flag_hs == 0)                                       //check if half displacement has been achieved or not
                      {
                          x_sub_half  = x0 + (dx_p - dx_half_p);              //set the half sub-cycle time particle position
                          cl_sub_half = cl;                                   //set the half sub-cycle time particle cell index
                          flag_hs     = 1;                                    //turn on flag indicating that half particle displacement has been achieved
                      }
                      dx_p    = 0.0;                                            //set the total remaining displacement to zero since particle finished traversing within sub-cycle
                      flag_rt = 1;                                            //turn on flag indicating that ray tracing has completed
                      break;                            
                  default:                                                    //if particle reverses direction at face
                      x_sub_f     = x_px;                                     //set final sub-cycle time particle position to the face in which the particle reversed direction
                      x_sub_half  = x_px;                                     //set half sub-cycle time particle position to the face in which the particle reversed direction
                      v_sub_f     = v_px;                                     //set final sub-cycle time particle velocity to the actual value in which the non-linear picard iteration converged
                      v_sub_half  = (v_sub_old + v_sub_f);                    //set half sub-cycle time particle velocity to the half value in which the non-linear picard iteration converged
                      dtau        = 0.0;                                      //set sub-cycle time to zero since particle just reversed direction but didn't move
                      flag_fx     = 1;                                        //turn on flag indicating that the particle has crossed a face
                      break;
              }
          }
          //end of ray tracing, now calculate convergence of particle traj.
          //E-field interpolation is hard-coded for CPU speed optimization purpose
          if (this->data->temp_order == 1) 
          {
              x_E     = x_sub_f;                              //set particle position for E-field interpolation
              cl_E    = cl_sub_f;                             //set particle cell index for E-field interpolation
          }
          else if (this->data->temp_order == 2) 
          {
              x_E     = x_sub_half;                           //set particle position for E-field interpolation
              cl_E    = cl_sub_half;                          //set particle cell index for E-field interpolation
          }
          ddx     = x_E - this->data->xpos_face[cl_E];                    //calculate distance of particle from the left face of cell
          dEdx    = (E_m[cl_E+1] - E_m[cl_E])*dx_recip;       //calculate the gradient of electric field in which particle is in
          E_p    = E_m[cl_E] + dEdx*ddx;                      //actually interpolate the electric field to the particle
          if (part_cf_xf0 == 1)                               //check if the particle reversed direction at the cell face
          {
              x_sub_f     = x_px;                             //store final sub-cycle time particle position to the cell face in which particle reversed direction
              x_sub_half  = x_px;                             //store half sub-cycle time particle position to the cell face in which particle reversed direction
              v_sub_f     = v_px;                             //store final sub-cycle time particle velocity to the value in which the non-linear picard iteration converged
              v_sub_half  = v_px;                             //store half sub-cycle time particle velocity to the value in which the non-linear picard iteration converged and average it
              cl_sub_f    = cl_px;                            //store the final sub-cycle time particle cell index to the cell in which face is associated with
              cl_sub_half = cl_px;                            //store the half sub-cycle time particle cell index to the cell in which face is associated with
              dtau        = 0.0;                              //set the sub-cycle time to zero since particle did not "move" but just reversed direction
          }
          else                                                //check if the particle didn't reverse direction at the cell face
          {
              v_sub_f     = v_sub_old + dtau*qm*E_p;          //calculate final velocity based on position of particle and the interpolated electric field
              cl_sub_f    = cl;                               //set the final sub-cycle time particle cell index
              v_sub_half  = (v_sub_f + v_sub_old)*two_recip;  //calculate the half sub-cycle time particle velocity as an average of final and old time value
          }
          rel_x       = (x_sub_f - x_k)/(x_sub_f);            //calculate the relative difference of particle position for convergence purposes
          rel_v       = (v_sub_f - v_k)/(v_sub_f);            //calculat ethe relative difference of particle velocity for convergence purposes
          x_k         = x_sub_f;                              //update reference particle position for convergence purposes
          v_k         = v_sub_f;                              //update reference particle velocity for convergence purposes
          if (fabs(rel_x) <= tol_picard && fabs(rel_v) <= tol_picard ) //check for convergence condition
          {
              flag_conv = 1;                                  //turn flag on to indicate that the non-linear picard iteration converged for both particle position and velocity
          }
          if (k == 10) //terminate if more than 10 picard iteration is required. Something is wrong
          {
              flag_conv = 1;
              cout <<"More than 10 picard iteration was required for the" << pIndx << " th particle, terminating picard iteration\n" << endl;
          }
        }            
        //if particle crossed a surface, stop particle and recalculate the position and velocity            
        if (flag_fx == 1 && part_cf_xf0 == 0)                               //if particle crossed a face, adjust dtau to stop particle at cell face
        {
            x_sub_f     = x_xf;                                             //set final sub-cycle time particle position to be the cell face in which particle crossed approached from the direction of particle
            x_sub_half  = (x_sub_f + x_sub_old)*two_recip;                  //set half sub-cycle time particle position to be the average of full time and old time position
            cl_sub_half = cl_sub_old;                                       //set half sub-cycle time particle cell index to be the initial cell index since we're stopping the particle at the cell face
            //E-field interpolation is hard coded for now for optimization
            if (this->data->temp_order == 1) 
            {
                x_E     = x_sub_half;
                cl_E    = cl_sub_half;
            }
            else if (this->data->temp_order == 2) 
            {
                x_E     = x_sub_half;
                cl_E    = cl_sub_half;
            }
            ddx     = x_E - this->data->xpos_face[cl_E];                                //calculate the distance of particle from left face of cell in which it is in
            dEdx    = (E_m[cl_E + 1] - E_m[cl_E])*dx_recip;                 //calculate the gradient of electric field in which the particle is in
            E_p     = E_m[cl_E] + dEdx*ddx;                                 //actually interpolate the electric field to the particle
            //back calculate the sub-cycle time step size based on the fact that we KNOW where the particle is since we're stopping them at cell surface
            //based on:A_tau*dtau^2 + B_tau*dtau^1 + C_tau = 0
            A_tau       = -two_recip*qm*E_p;                                //calculate the coefficient to calculate the new sub-cycle time to stop the particle at cell face for 2nd order term
            B_tau       = -v_sub_old;                                       //calculate the coefficient to calculate the new sub-cycle time to stop the particle at cell face for 1st order term
            C_tau       = x_sub_f - x_sub_old;                              //calculate the coefficient to caluclate the new sub-cycle time to stop the particle at cell face for 0th order term
            A_tau_recip = 1.0/A_tau;                                        //calculate reciprical of A_tau
            dtau_pm     = sqrt(B_tau*B_tau - 4.0*A_tau*C_tau);              //calculate the the sqrt{} for the root for quadratic equation
            dtau_p      = two_recip*A_tau_recip*(-B_tau + dtau_pm);         //calculate the positive of the sqrt{} for the root
            dtau_m      = two_recip*A_tau_recip*(-B_tau - dtau_pm);         //calculate the negative of the sqrt{} for the root
            v_sub_f_p   = v_sub_old + dtau_p*qm*E_p;                        //calculate the final sub-cycle velocity of particle for positive root of dtau
            if (c_f_n_xf*v_sub_f_p >= 0.0)                                  //check if the normal vector of face in which particle crossed is parallel to the velocity calculated by the positive root of dtau
            {
                //check if dtau_p is less than dtau_m and other conditions
                if (dtau_p <= dtau_m && dtau_m > 0.0) 
                {
                    dtau    = dtau_p;
                }
                else if (dtau_p > dtau_m && dtau_m > 0.0) 
                {
                    dtau    = dtau_m;
                }
                else if (dtau_m < 0.0) 
                {
                    dtau    = dtau_p;
                }
            }
            else                                                            //if the normal vector of face in which particle crossed is anti-parallel to the velocity calculated by the positive root of dtau, then dtau must be dtau_m for physical reasons
            {
                dtau    = dtau_m;
            }
            //now with the back calculated sub-cycle time size, back calculate the velocities
            v_sub_f     = v_sub_old + dtau*qm*E_p;                          //calculate the new final sub-cycle time velocity that stops particle at cell face
            v_sub_half  = (v_sub_f + v_sub_old)*two_recip;                  //calculate time aeraged sub-cycle half time vleocity 
            x_sub_f     = xf_real;                                          //set the particle position at the face position approached form the opposite direction as the particle velocity
            cl_sub_f    = cl_xf;                                            //set the final sub-cycle particle cell index to the cell index in which the face is associated with when approached form the opposite direction of the particle velocity
        }
        //The following is a redundancy since the particle tracking algorithm is a bit more complicated...
        else if (part_cf_xf0 == 1)                                          //if particle reversed direction at cell face in sub-cycle
        {
            x_sub_f     = x_px;                                             //set the final sub-cycle time particle position at the cell face
            x_sub_half  = x_px;                                             //set the half sub-cycle time particle position at the cell face
            v_sub_f     = v_px;                                             //set the final sub-cycle time particle velocity
            v_sub_half  = v_px;                                             //set the half sub-cycle time particle velocity
            cl_sub_f    = cl_px;                                            //set the final sub-cycle time particle cell index
            cl_sub_half = cl_px;                                            //set the half sub-cycle time particle cell index
            dtau        = 0.0;                                              //set the sub-cycle time to be equal to zero
        }            
//Will Taitano edit: 06/28/2012 for debugging purpose only
          cout << "pIndx = " << pIndx << endl;
        cout << "sc_count = " << sc_counter << endl;
        cout << "x_sub_old = " << x_sub_old << ", x_sub_half = " << x_sub_half << ", x_sub_f = " << x_sub_f << endl;
        cout << "v_sub_old = " << v_sub_old << ", v_sub_half = " << v_sub_half << ", v_sub_f = " << v_sub_f << endl;
        cout << "cl_sub_old = " << cl_sub_old << ", cl_sub_half = " << cl_sub_half << ", cl_sub_f = " << cl_sub_f << endl;
          cout << "dtau = " << dtau << endl;
        // Tally half time moment quantities (s1_p_half and s2_p_half)
        half->tallyMomentCalculation(sIndx, x_sub_half, v_sub_half, dtau);
            
        // Update sub step residual to see if loop can terminate
        dt_res -= dtau;
        if (dt_res <= dt_res_tol)
          flag_sc = 1;
      }
//        cout << "The dtau_res for the " << pIndx << " th particle is, dt_res = " << dt_res << endl;
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

    if (this->data->temp_order == ORDER_1_BACK_EULER) 
    {
      for (int cIndx = 0; cIndx < this->data->nx; cIndx++) {
        double cc = dt_recip * (r[cIndx] - r_old[cIndx]) +
                    dx_recip * (s1[cIndx + 1] - s1[cIndx]);
        cc_rms += cc * cc;
      }
    } else if (this->data->temp_order == ORDER_2_CRANK_NICOLSON) 
    {
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