#include "DiscretePDE_1D3Vem.h"
#include <omp.h>

DiscretePDE_1D3Vem::DiscretePDE_1D3Vem(const Teuchos::RCP<MapManager> &map_manager,
	      const Teuchos::RCP<SimParams> &_params,
	      Epetra_Comm* _comm)
{
  comm       = _comm;
  simParams  = _params;
  mapManager = map_manager;

  stagMesh   = Teuchos::rcp(new StagMesh(mapManager, simParams, comm) );
  OvX        = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );
  OvXold     = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );

  ImportParams();
};

//DiscretePDE::~DiscretePDE(){}

////////////////////////////////////////////////////////////////////
// USER: Edit any constants you may need in SimParams.cpp  //
///////////////////////////////////////////////////////////////////
void DiscretePDE_1D3Vem::ImportParams()
{

  me_h   = simParams->GetParamsPhys()->get<double>("me_h");
  mi_h   = simParams->GetParamsPhys()->get<double>("mi_h");
  qe_h   = simParams->GetParamsPhys()->get<double>("qe_h");
  qi_h   = simParams->GetParamsPhys()->get<double>("qi_h");

  zeta = simParams->GetParamsPhys()->get<double>("zeta");
  xi    = simParams->GetParamsPhys()->get<double>("xi");
  dx    = simParams->GetParamsPhys()->get<double>("dx");
  dt    = simParams->GetParamsPhys()->get<double>("dt");
}



////////////////////////////////////////////////////////////////////
// USER: This is where users must build the main residual calc.   //
// Propogate the enum below with names for each of your unknowns. //
// The number of these unknowns must match "Number Unknowns"      //
///////////////////////////////////////////////////////////////////
void DiscretePDE_1D3Vem::EvaluateResidual( const Epetra_Vector& x,
				    const Epetra_Vector& xold,
				    const Epetra_Vector& gamma,
				    Epetra_Vector& res )
{
  // Build an enum of variable names

	using namespace EM1D3V;

//  std::cout << "gamma = " << std::endl;
//  std::cout << gamma;

  // Load x into a discrete operator wrapper
  OvX->SetVector(x);

  // Load xold into a discrte operator wrapper
  OvXold->SetVector(xold);

  //Global Operations (ie. averages)
  double pe_avg = 0.0;
  double pi_avg = 0.0;

  double pe_avg_old = 0.0;
  double pi_avg_old = 0.0;

  // double pe_avg = OvX->Mean(p_e);
  // double pi_avg = OvX->Mean(p_i);

  // double pe_avg_old = OvXold->Mean(p_e);
  // double pi_avg_old = OvXold->Mean(p_i);

  double jmeanx1 = 0;
  double jmeany1 = 0;
  double jmeanz1 = 0;

  // Loop over elements
  int N = mapManager->LocSdNumElemX;
  for(int i = 0; i<N; i++)
  {
      // Convert elem index to first point index
      int j = mapManager->LocElem_to_LocPtSd(i,0);

      jmeanx1 += qe_h*OvX->GetF2C(px_e, i) + qi_h*OvX->GetF2C(px_i,i);
      jmeany1 += qe_h*OvX->Get(py_e, i) + qi_h*OvX->Get(py_i,i);
      jmeanz1 += qe_h*OvX->Get(pz_e, i) + qi_h*OvX->Get(pz_i,i);


  }

  jmeanx1 /= (double)(N);
  jmeany1 /= (double)(N);
  jmeanz1 /= (double)(N);

//  for(int i = 0; i<N; i++)
//  {
//
//	  OvX->Get(Ey,i) = -2.0*(OvX->Get(Ay,i) - OvXold->Get(Ay,i))/dt - OvXold->Get(Ey,i);
//	  OvX->Get(Ez,i) = -2.0*(OvX->Get(Az,i) - OvXold->Get(Az,i))/dt - OvXold->Get(Ez,i);
//
//
//  }

omp_set_num_threads(32);
#pragma omp parallel for
for(int i = 0; i<N; i++)
    {   
      // Convert elem index to first point index
      int j = mapManager->LocElem_to_LocPtSd(i,0);

//      printf("B[%i] = %e %e %e\n",i,OvX->Get(Bx,i),OvX->Get(By,i),OvX->Get(Bz,i));
//      printf("E[%i] = %e %e %e\n",i,OvX->Get(Ex,i),OvX->Get(Ey,i),OvX->Get(Ez,i));
//      printf("A[%i] = %e %e %e\n",i,OvX->Get(Ax,i),OvX->Get(Ay,i),OvX->Get(Az,i));
//      printf("pe[%i] = %e %e %e\n",i,OvX->Get(px_e,i),OvX->Get(py_e,i),OvX->Get(pz_e,i));
//      printf("pi[%i] = %e %e %e\n",i,OvX->Get(px_i,i),OvX->Get(py_i,i),OvX->Get(pz_i,i));


//	  double Ext = OvXold->Get(Ex,i) - dt*(qe_h*OvX->Get(px_e, i) + qi_h*OvX->Get(px_i,i) - jmeanx1)/(xi*xi);
	  double Eyt = -2.0*(OvX->Get(Ay,i) - OvXold->Get(Ay,i))/dt - OvXold->Get(Ey,i);
	  double Ezt = -2.0*(OvX->Get(Az,i) - OvXold->Get(Az,i))/dt - OvXold->Get(Ez,i);

	  double Ext = OvX->Get(Ex,i);
//	  double Eyt = OvX->Get(Ey,i);
//	  double Ezt = OvX->Get(Ez,i);

//      double tmp_W = (OvX->Get(Sxx_e, i-1)+OvXold->Get(Sxx_e, i-1))*(OvX->Get(n_e,i-1)+OvXold->Get(n_e,i-1));
//      double tmp_C = (OvX->Get(Sxx_e, i)+OvXold->Get(Sxx_e, i)) * (OvX->Get(n_e,i)+OvXold->Get(n_e,i));
//      double tmp_E = (OvX->Get(Sxx_e, i+1)+OvXold->Get(Sxx_e, i+1)) * (OvX->Get(n_e,i+1)+OvXold->Get(n_e,i+1));
//      double Dx_Sxx_ne = 0.25*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = 0.25*(OvX->Get(Sxy_e, i-1)+OvXold->Get(Sxy_e, i-1)) * (OvX->Get(n_e,i-1)+OvXold->Get(n_e,i-1));
//      tmp_C = 0.25*(OvX->Get(Sxy_e, i)+OvXold->Get(Sxy_e, i)) * (OvX->Get(n_e,i)+OvXold->Get(n_e,i));
//      tmp_E = 0.25*(OvX->Get(Sxy_e, i+1)+OvXold->Get(Sxy_e, i+1)) * (OvX->Get(n_e,i+1)+OvXold->Get(n_e,i+1));
//      double Dx_Sxy_ne = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = 0.25*(OvX->Get(Sxz_e, i-1)+OvXold->Get(Sxz_e, i-1)) * (OvX->Get(n_e,i-1)+OvXold->Get(n_e,i-1));
//      tmp_C = 0.25*(OvX->Get(Sxz_e, i)+OvXold->Get(Sxz_e, i)) * (OvX->Get(n_e,i)+OvXold->Get(n_e,i));
//      tmp_E = 0.25*(OvX->Get(Sxz_e, i+1)+OvXold->Get(Sxz_e, i+1)) * (OvX->Get(n_e,i+1)+OvXold->Get(n_e,i+1));
//      double Dx_Sxz_ne = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      // Ion
//      tmp_W = 0.25*(OvX->Get(Sxx_i, i-1)+OvXold->Get(Sxx_i, i-1)) * (OvX->Get(n_i,i-1)+OvXold->Get(n_i,i-1));
//      tmp_C = 0.25*(OvX->Get(Sxx_i, i)+OvXold->Get(Sxx_i, i)) * (OvX->Get(n_i,i)+OvXold->Get(n_i,i));
//      tmp_E = 0.25*(OvX->Get(Sxx_i, i+1)+OvXold->Get(Sxx_i, i+1)) * (OvX->Get(n_i,i+1)+OvXold->Get(n_i,i+1));
//      double Dx_Sxx_ni = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = 0.25*(OvX->Get(Sxy_i, i-1)+OvXold->Get(Sxy_i, i-1)) * (OvX->Get(n_i,i-1)+OvXold->Get(n_i,i-1));
//      tmp_C = 0.25*(OvX->Get(Sxy_i, i)+OvXold->Get(Sxy_i, i)) * (OvX->Get(n_i,i)+OvXold->Get(n_i,i));
//      tmp_E = 0.25*(OvX->Get(Sxy_i, i+1)+OvXold->Get(Sxy_i, i+1)) * (OvX->Get(n_i,i+1)+OvXold->Get(n_i,i+1));
//      double Dx_Sxy_ni = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = 0.25*(OvX->Get(Sxz_i, i-1)+OvXold->Get(Sxz_i, i-1)) * (OvX->Get(n_i,i-1)+OvXold->Get(n_i,i-1));
//      tmp_C = 0.25*(OvX->Get(Sxz_i, i)+OvXold->Get(Sxz_i, i)) * (OvX->Get(n_i,i)+OvXold->Get(n_i,i));
//      tmp_E = 0.25*(OvX->Get(Sxz_i, i+1)+OvXold->Get(Sxz_i, i+1)) * (OvX->Get(n_i,i+1)+OvXold->Get(n_i,i+1));
//      double Dx_Sxz_ni = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//

      double tmp_W = (OvX->Get(Sxx_e, i-1)+3.0*OvXold->Get(Sxx_e, i-1))*(OvX->Get(n_e,i-1)+3.0*OvXold->Get(n_e,i-1));
      double tmp_C = (OvX->Get(Sxx_e, i)+3.0*OvXold->Get(Sxx_e, i)) * (OvX->Get(n_e,i)+3.0*OvXold->Get(n_e,i));
      double tmp_E = (OvX->Get(Sxx_e, i+1)+3.0*OvXold->Get(Sxx_e, i+1)) * (OvX->Get(n_e,i+1)+3.0*OvXold->Get(n_e,i+1));
      double Dx_Sxx_ne = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);

      tmp_W = (OvX->Get(Sxy_e, i-1)+3.*OvXold->Get(Sxy_e, i-1)) * (OvX->Get(n_e,i-1)+3.*OvXold->Get(n_e,i-1));
      tmp_C = (OvX->Get(Sxy_e, i)+3.*OvXold->Get(Sxy_e, i)) * (OvX->Get(n_e,i)+3.*OvXold->Get(n_e,i));
      tmp_E = (OvX->Get(Sxy_e, i+1)+3.*OvXold->Get(Sxy_e, i+1)) * (OvX->Get(n_e,i+1)+3.*OvXold->Get(n_e,i+1));
      double Dx_Sxy_ne = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);

      tmp_W = (OvX->Get(Sxz_e, i-1)+3.*OvXold->Get(Sxz_e, i-1)) * (OvX->Get(n_e,i-1)+3.*OvXold->Get(n_e,i-1));
      tmp_C = (OvX->Get(Sxz_e, i)+3.*OvXold->Get(Sxz_e, i)) * (OvX->Get(n_e,i)+3.*OvXold->Get(n_e,i));
      tmp_E = (OvX->Get(Sxz_e, i+1)+3.*OvXold->Get(Sxz_e, i+1)) * (OvX->Get(n_e,i+1)+3.*OvXold->Get(n_e,i+1));
      double Dx_Sxz_ne = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);

      // Ion
      tmp_W = (OvX->Get(Sxx_i, i-1)+3.*OvXold->Get(Sxx_i, i-1)) * (OvX->Get(n_i,i-1)+3.*OvXold->Get(n_i,i-1));
      tmp_C = (OvX->Get(Sxx_i, i)+3.*OvXold->Get(Sxx_i, i)) * (OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i));
      tmp_E = (OvX->Get(Sxx_i, i+1)+3.*OvXold->Get(Sxx_i, i+1)) * (OvX->Get(n_i,i+1)+3.*OvXold->Get(n_i,i+1));
      double Dx_Sxx_ni = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);

      tmp_W = (OvX->Get(Sxy_i, i-1)+3.*OvXold->Get(Sxy_i, i-1)) * (OvX->Get(n_i,i-1)+3.*OvXold->Get(n_i,i-1));
      tmp_C = (OvX->Get(Sxy_i, i)+3.*OvXold->Get(Sxy_i, i)) * (OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i));
      tmp_E = (OvX->Get(Sxy_i, i+1)+3.*OvXold->Get(Sxy_i, i+1)) * (OvX->Get(n_i,i+1)+3.*OvXold->Get(n_i,i+1));
      double Dx_Sxy_ni = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);

      tmp_W = (OvX->Get(Sxz_i, i-1)+3.*OvXold->Get(Sxz_i, i-1)) * (OvX->Get(n_i,i-1)+3.*OvXold->Get(n_i,i-1));
      tmp_C = (OvX->Get(Sxz_i, i)+3.*OvXold->Get(Sxz_i, i)) * (OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i));
      tmp_E = (OvX->Get(Sxz_i, i+1)+3.*OvXold->Get(Sxz_i, i+1)) * (OvX->Get(n_i,i+1)+3.*OvXold->Get(n_i,i+1));
      double Dx_Sxz_ni = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);



      // Fill each unknown
      //    Electron continuity equation
      res[j+n_e] = (OvX->Get(n_e, i) - OvXold->Get(n_e, i)) + dt*( OvX->DxF2C(px_e, i, dx)) - gamma[j + n_e];
      //    Ion continuity equation
      res[j+n_i] = (OvX->Get(n_i, i) - OvXold->Get(n_i, i)) + dt*( OvX->DxF2C(px_i, i, dx) ) - gamma[j + n_i];

      //    Electron momentum equation
      //    momentum electron-x
      res[j+px_e] = (OvX->Get(px_e, i) - OvXold->Get(px_e, i)) +
                          (dt/(2.0*me_h))*(
                                  me_h*Dx_Sxx_ne -
                                  0.0625*qe_h*(OvX->GetC2F(n_e,i)+3.*OvXold->GetC2F(n_e,i))*(3.*OvXold->Get(Ex,i)+Ext) -
                                  (
                                          0.125*qe_h*(OvX->GetC2F(py_e,i)+OvXold->GetC2F(py_e,i))*(OvX->DxC2F(Ay,i,dx) + 3.*OvXold->DxC2F(Ay,i,dx))+
                                          0.125*qe_h*(OvX->GetC2F(pz_e,i)+OvXold->GetC2F(pz_e,i))*(OvX->DxC2F(Az,i,dx) + 3.*OvXold->DxC2F(Az,i,dx))
                                  	  )
                                  )
                                   - gamma[j + px_e];

      //    momentum electron-y
      res[j+py_e] = (OvX->Get(py_e, i) - OvXold->Get(py_e, i)) +
                          (dt/(2.0*me_h))*(
                                  me_h*Dx_Sxy_ne -
                                  0.0625*qe_h*(OvX->Get(n_e,i)+3.*OvXold->Get(n_e,i))*(3.*OvXold->Get(Ey,i)+Eyt) -
                                  (
                                		  qe_h*OvX->Get(pz_e,i)*OvXold->Get(Bx,i)-
                                          0.125*qe_h*(OvX->GetF2C(px_e,i)+OvXold->GetF2C(px_e,i))*(OvX->Dx(Ay,i,dx) + 3.*OvXold->Dx(Ay,i,dx))
                                  	  )
                                  )
                                   - gamma[j + py_e];

      //    momentum electron-z
      res[j+pz_e] = (OvX->Get(pz_e, i) - OvXold->Get(pz_e, i)) +
                          (dt/(2.0*me_h))*(
                                  me_h*Dx_Sxz_ne -
                                  0.0625*qe_h*(OvX->Get(n_e,i)+3.*OvXold->Get(n_e,i))*(3.*OvXold->Get(Ez,i)+Ezt) -
                                  (
                                		  -qe_h*OvX->Get(py_e,i)*OvXold->Get(Bx,i)-
                                          0.125*qe_h*(OvX->GetF2C(px_e,i)+OvXold->GetF2C(px_e,i))*(OvX->Dx(Az,i,dx) + 3.*OvXold->Dx(Az,i,dx))
                                  	  )
                                  )
                                   - gamma[j + pz_e];



      //    Ion momentum equation
      //    momentum ion-x
      res[j+px_i] = (OvX->Get(px_i, i) - OvXold->Get(px_i, i)) +
                          (dt/(2.0*mi_h))*(
                                  mi_h*Dx_Sxx_ni -
                                  0.0625*qi_h*(OvX->GetC2F(n_i,i)+3.*OvXold->GetC2F(n_i,i))*(3.*OvXold->Get(Ex,i)+Ext) -
                                  (
                                          0.125*qi_h*(OvX->GetC2F(py_i,i)+OvXold->GetC2F(py_i,i))*(OvX->DxC2F(Ay,i,dx) + 3.*OvXold->DxC2F(Ay,i,dx)) +
                                          0.125*qi_h*(OvX->GetC2F(pz_i,i)+OvXold->GetC2F(pz_i,i))*(OvX->DxC2F(Az,i,dx) + 3.*OvXold->DxC2F(Az,i,dx))
                                  	  )
                                  )
                                   - gamma[j + px_i];

      //    momentum electron-y
      res[j+py_i] = (OvX->Get(py_i, i) - OvXold->Get(py_i, i)) +
                          (dt/(2.0*mi_h))*(
                                  mi_h*Dx_Sxy_ni -
                                  0.0625*qi_h*(OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i))*(3.*OvXold->Get(Ey,i)+Eyt) -
                                  (
                                		  qi_h*OvX->Get(pz_i,i)*OvXold->Get(Bx,i)-
                                          0.125*qi_h*(OvX->GetF2C(px_i,i)+OvXold->GetF2C(px_i,i))*(OvX->Dx(Ay,i,dx) + 3.*OvXold->Dx(Ay,i,dx))
                                  	  )
                                  )
                                   - gamma[j + py_i];

      //    momentum electron-z
      res[j+pz_i] = (OvX->Get(pz_i, i) - OvXold->Get(pz_i, i)) +
                          (dt/(2.0*mi_h))*(
                                  mi_h*Dx_Sxz_ni -
                                  0.0625*qi_h*(OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i))*(3.*OvXold->Get(Ez,i)+Ezt) -
                                  (
                                		  -qi_h*OvX->Get(py_i,i)*OvXold->Get(Bx,i)-
                                          0.125*qi_h*(OvX->GetF2C(px_i,i)+OvXold->GetF2C(px_i,i))*(OvX->Dx(Az,i,dx) + 3.*OvXold->Dx(Az,i,dx))
                                  	  )
                                  )
                                   - gamma[j + pz_i];


//
//      // Fill each unknown
//      //    Electron continuity equation
//      res[j+n_e] = (OvX->Get(n_e, i) - OvXold->Get(n_e, i)) + dt*( OvX->DxF2C(px_e, i, dx)) - gamma[j + n_e];
//      //    Ion continuity equation
//      res[j+n_i] = (OvX->Get(n_i, i) - OvXold->Get(n_i, i)) + dt*( OvX->DxF2C(px_i, i, dx) ) - gamma[j + n_i];
//
//      //    Electron momentum equation
//      //    momentum electron-x
//      res[j+px_e] = (OvX->Get(px_e, i) - OvXold->Get(px_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxx_ne -
//                                  0.25*(qe_h*OvX->GetC2F(n_e,i)+qe_h*OvXold->GetC2F(n_e,i))*(OvXold->Get(Ex,i)+OvX->Get(Ex,i)) -
//                                  (
//                                          0.5*qe_h*OvX->GetC2F(py_e,i)*(OvX->DxC2F(Ay,i,dx) + OvXold->DxC2F(Ay,i,dx))+
//                                          0.5*qe_h*OvX->GetC2F(pz_e,i)*(OvX->DxC2F(Az,i,dx) + OvXold->DxC2F(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + px_e];
//
//      //    momentum electron-y
//      res[j+py_e] = (OvX->Get(py_e, i) - OvXold->Get(py_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxy_ne -
//                                  0.25*(qe_h*OvX->Get(n_e,i)+qe_h*OvXold->Get(n_e,i))*(OvXold->Get(Ey,i)+OvX->Get(Ey,i)) -
//                                  (
//                                		  qe_h*OvX->Get(pz_e,i)*OvXold->Get(Bx,i)-
//                                          0.5*qe_h*OvX->GetF2C(px_e,i)*(OvX->Dx(Ay,i,dx) + OvXold->Dx(Ay,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + py_e];
//
//      //    momentum electron-z
//      res[j+pz_e] = (OvX->Get(pz_e, i) - OvXold->Get(pz_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxz_ne -
//                                  0.25*(qe_h*OvX->Get(n_e,i)+qe_h*OvXold->Get(n_e,i))*(OvXold->Get(Ez,i)+OvX->Get(Ez,i)) -
//                                  (
//                                		  -qe_h*OvX->Get(py_e,i)*OvXold->Get(Bx,i)-
//                                          0.5*qe_h*OvX->GetF2C(px_e,i)*(OvX->Dx(Az,i,dx) + OvXold->Dx(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + pz_e];
//
//
//
//      //    Ion momentum equation
//      //    momentum ion-x
//      res[j+px_i] = (OvX->Get(px_i, i) - OvXold->Get(px_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxx_ni -
//                                  0.25*(qi_h*OvX->GetC2F(n_i,i)+qi_h*OvXold->GetC2F(n_i,i))*(OvXold->Get(Ex,i)+OvX->Get(Ex,i)) -
//                                  (
//                                          0.5*qi_h*OvX->GetC2F(py_i,i)*(OvX->DxC2F(Ay,i,dx) + OvXold->DxC2F(Ay,i,dx)) +
//                                          0.5*qi_h*OvX->GetC2F(pz_i,i)*(OvX->DxC2F(Az,i,dx) + OvXold->DxC2F(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + px_i];
//
//      //    momentum electron-y
//      res[j+py_i] = (OvX->Get(py_i, i) - OvXold->Get(py_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxy_ni -
//                                  0.25*(qi_h*OvX->Get(n_i,i)+qi_h*OvXold->Get(n_i,i))*(OvXold->Get(Ey,i)+OvX->Get(Ey,i)) -
//                                  (
//                                		  qi_h*OvX->Get(pz_i,i)*OvXold->Get(Bx,i)-
//                                          0.5*qi_h*OvX->GetF2C(px_i,i)*(OvX->Dx(Ay,i,dx) + OvXold->Dx(Ay,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + py_i];
//
//      //    momentum electron-z
//      res[j+pz_i] = (OvX->Get(pz_i, i) - OvXold->Get(pz_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxz_ni -
//                                  0.25*(qi_h*OvX->Get(n_i,i)+qi_h*OvXold->Get(n_i,i))*(OvXold->Get(Ez,i)+OvX->Get(Ez,i)) -
//                                  (
//                                		  -qi_h*OvX->Get(py_i,i)*OvXold->Get(Bx,i)-
//                                          0.5*qi_h*OvX->GetF2C(px_i,i)*(OvX->Dx(Az,i,dx) + OvXold->Dx(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + pz_i];


//      //    Construct nonlinear terms
//      //    Electron
//      //    Compute advection for momentum x
//      double tmp_W = OvX->Get(Sxx_e, i-1) * OvX->Get(n_e,i-1);
//      double tmp_C = OvX->Get(Sxx_e, i) * OvX->Get(n_e,i);
//      double tmp_E = OvX->Get(Sxx_e, i+1) * OvX->Get(n_e,i+1);
//      double Dx_Sxx_ne = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = OvXold->Get(Sxx_e, i-1) * OvXold->Get(n_e,i-1);
//      tmp_C = OvXold->Get(Sxx_e, i) * OvXold->Get(n_e,i);
//      tmp_E = OvXold->Get(Sxx_e, i+1) * OvXold->Get(n_e,i+1);
//      double Dx_Sxx_ne_old = OvXold->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//      //    Compute advection for momentum y
//      tmp_W = OvX->Get(Sxy_e, i-1) * OvX->Get(n_e,i-1);
//      tmp_C = OvX->Get(Sxy_e, i) * OvX->Get(n_e,i);
//      tmp_E = OvX->Get(Sxy_e, i+1) * OvX->Get(n_e,i+1);
//      double Dx_Sxy_ne = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = OvXold->Get(Sxy_e, i-1) * OvXold->Get(n_e,i-1);
//      tmp_C = OvXold->Get(Sxy_e, i) * OvXold->Get(n_e,i);
//      tmp_E = OvXold->Get(Sxy_e, i+1) * OvXold->Get(n_e,i+1);
//      double Dx_Sxy_ne_old = OvXold->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//      //    Compute advection for momentum z
//      tmp_W = OvX->Get(Sxz_e, i-1) * OvX->Get(n_e,i-1);
//      tmp_C = OvX->Get(Sxz_e, i) * OvX->Get(n_e,i);
//      tmp_E = OvX->Get(Sxz_e, i+1) * OvX->Get(n_e,i+1);
//      double Dx_Sxz_ne = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = OvXold->Get(Sxz_e, i-1) * OvXold->Get(n_e,i-1);
//      tmp_C = OvXold->Get(Sxz_e, i) * OvXold->Get(n_e,i);
//      tmp_E = OvXold->Get(Sxz_e, i+1) * OvXold->Get(n_e,i+1);
//      double Dx_Sxz_ne_old = OvXold->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      //    ion
//      //    Compute advection for momentum x
//      tmp_W = OvX->Get(Sxx_i, i-1) * OvX->Get(n_i,i-1);
//      tmp_C = OvX->Get(Sxx_i, i) * OvX->Get(n_i,i);
//      tmp_E = OvX->Get(Sxx_i, i+1) * OvX->Get(n_i,i+1);
//      double Dx_Sxx_ni = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = OvXold->Get(Sxx_i, i-1) * OvXold->Get(n_i,i-1);
//      tmp_C = OvXold->Get(Sxx_i, i) * OvXold->Get(n_i,i);
//      tmp_E = OvXold->Get(Sxx_i, i+1) * OvXold->Get(n_i,i+1);
//      double Dx_Sxx_ni_old = OvXold->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//      //    Compute advection for momentum y
//      tmp_W = OvX->Get(Sxy_i, i-1) * OvX->Get(n_i,i-1);
//      tmp_C = OvX->Get(Sxy_i, i) * OvX->Get(n_i,i);
//      tmp_E = OvX->Get(Sxy_i, i+1) * OvX->Get(n_i,i+1);
//      double Dx_Sxy_ni = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = OvXold->Get(Sxy_i, i-1) * OvXold->Get(n_i,i-1);
//      tmp_C = OvXold->Get(Sxy_i, i) * OvXold->Get(n_i,i);
//      tmp_E = OvXold->Get(Sxy_i, i+1) * OvXold->Get(n_i,i+1);
//      double Dx_Sxy_ni_old = OvXold->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//      //    Compute advection for momentum z
//      tmp_W = OvX->Get(Sxz_i, i-1) * OvX->Get(n_i,i-1);
//      tmp_C = OvX->Get(Sxz_i, i) * OvX->Get(n_i,i);
//      tmp_E = OvX->Get(Sxz_i, i+1) * OvX->Get(n_i,i+1);
//      double Dx_Sxz_ni = OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = OvXold->Get(Sxz_i, i-1) * OvXold->Get(n_i,i-1);
//      tmp_C = OvXold->Get(Sxz_i, i) * OvXold->Get(n_i,i);
//      tmp_E = OvXold->Get(Sxz_i, i+1) * OvXold->Get(n_i,i+1);
//      double Dx_Sxz_ni_old = OvXold->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
////      // Fill each unknown
////      //    Electron continuity equation
//      res[j+n_e] = (OvX->Get(n_e, i) - OvXold->Get(n_e, i)) + dt*( OvX->DxF2C(px_e, i, dx)) - gamma[j + n_e];
////      //    Ion continuity equation
//      res[j+n_i] = (OvX->Get(n_i, i) - OvXold->Get(n_i, i)) + dt*( OvX->DxF2C(px_i, i, dx) ) - gamma[j + n_i];
//
//      //    Electron momentum equation
//      //    momentum electron-x
//      res[j+px_e] = (OvX->Get(px_e, i) - OvXold->Get(px_e, i)) +
//                          (dt/(2.0*me_h))*(
//                        		  0.5*me_h*Dx_Sxx_ne -
//                                  0.5*(qe_h*OvX->GetC2F(n_e,i)+qe_h*OvXold->GetC2F(n_e,i))*(OvXold->Get(Ex,i)+Ext) -
//                                  (
//                                          qe_h*(OvX->GetC2F(py_e,i))*OvX->DxC2F(Ay,i,dx) +
//                                          qe_h*(OvX->GetC2F(pz_e,i))*OvX->DxC2F(Az,i,dx)
//                                  	  ) +
//                                  	0.5*me_h*Dx_Sxx_ne_old -
//                                  (
//                                          qe_h*(OvX->GetC2F(py_e,i))*OvXold->DxC2F(Ay,i,dx) +
//                                          qe_h*(OvX->GetC2F(pz_e,i))*OvXold->DxC2F(Az,i,dx)
//                                  	  ))
//                                   - gamma[j + px_e];
//
//      //    momentum electron-y
//      res[j+py_e] = (OvX->Get(py_e, i) - OvXold->Get(py_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  0.5*me_h*Dx_Sxy_ne -
//                                  0.5*(qe_h*OvX->Get(n_e,i)+qe_h*OvXold->Get(n_e,i))*(OvXold->Get(Ey,i)+Eyt) -
//                                  (
//                                		  qe_h*OvX->Get(pz_e,i)*OvXold->Get(Bx,i)-
//                                          qe_h*(OvX->GetF2C(px_e,i))*OvX->Dx(Ay,i,dx)
//                                  	  ) +
//                                  	0.5*me_h*Dx_Sxy_ne_old -
//                                  (
//                                		  qe_h*OvXold->Get(pz_e,i)*OvXold->Get(Bx,i)-
//                                		  qe_h*(OvX->GetF2C(px_e,i))*OvXold->Dx(Ay,i,dx)
//                                  	  ))
//                                   - gamma[j + py_e];
//
//      //    momentum electron-z
//      res[j+pz_e] = (OvX->Get(pz_e, i) - OvXold->Get(pz_e, i)) +
//                          (dt/(2.0*me_h))*(
//                        		  0.5*me_h*Dx_Sxz_ne -
//                                  0.5*(qe_h*OvX->Get(n_e,i)+qe_h*OvXold->Get(n_e,i))*(OvXold->Get(Ez,i)+Ezt) -
//                                  (
//                                		  -qe_h*OvX->Get(py_e,i)*OvXold->Get(Bx,i)-
//                                          qe_h*(OvX->GetF2C(px_e,i))*OvX->Dx(Az,i,dx)
//                                  	  ) +
//                                  	0.5*me_h*Dx_Sxz_ne_old -
//                                  (
//                                		  -qe_h*OvXold->Get(py_e,i)*OvXold->Get(Bx,i)-
//                                          qe_h*(OvX->GetF2C(px_e,i))*OvXold->Dx(Az,i,dx)
//                                  	  ))
//                                   - gamma[j + pz_e];
//
//
//
//      //    Ion momentum equation
//      //    momentum ion-x
//      res[j+px_i] = (OvX->Get(px_i, i) - OvXold->Get(px_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                        		  0.5*mi_h*Dx_Sxx_ni -
//                                  0.5*(qi_h*OvX->GetC2F(n_i,i)+qi_h*OvXold->GetC2F(n_i,i))*(OvXold->Get(Ex,i)+Ext) -
//                                  (
//                                		  qi_h*(OvX->GetC2F(py_i,i))*OvX->DxC2F(Ay,i,dx) +
//                                		  qi_h*(OvX->GetC2F(pz_i,i))*OvX->DxC2F(Az,i,dx)
//                                  	  ) +
//                                  	0.5*mi_h*Dx_Sxx_ni_old -
//
//                                  (
//                                          qi_h*(OvX->GetC2F(py_i,i))*OvXold->DxC2F(Ay,i,dx) +
//                                          qi_h*(OvX->GetC2F(pz_i,i))*OvXold->DxC2F(Az,i,dx)
//                                  	  ))
//                                   - gamma[j + px_i];
//
//      //    momentum electron-y
//      res[j+py_i] = (OvX->Get(py_i, i) - OvXold->Get(py_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                        		  0.5*mi_h*Dx_Sxy_ni -
//                                  0.5*(qi_h*OvX->Get(n_i,i)+qi_h*OvXold->Get(n_i,i))*(OvXold->Get(Ey,i)+Eyt) -
//                                  (
//                                		  qi_h*OvX->Get(pz_i,i)*OvXold->Get(Bx,i)-
//                                          qi_h*(OvX->GetF2C(px_i,i))*OvX->Dx(Ay,i,dx)
//                                  	  ) +
//                                  	0.5*mi_h*Dx_Sxy_ni_old -
//                                  (
//                                		  qi_h*OvXold->Get(pz_i,i)*OvXold->Get(Bx,i)-
//                                		  qi_h*(OvX->GetF2C(px_i,i))*OvXold->Dx(Ay,i,dx)
//                                  	  ))
//                                   - gamma[j + py_i];
//
//      //    momentum electron-z
//      res[j+pz_i] = (OvX->Get(pz_i, i) - OvXold->Get(pz_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                        		  0.5*mi_h*Dx_Sxz_ni -
//                                  0.5*(qi_h*OvX->Get(n_i,i)+qi_h*OvXold->Get(n_i,i))*(OvXold->Get(Ez,i)+Ezt) -
//                                  (
//                                		  -qi_h*OvX->Get(py_i,i)*OvXold->Get(Bx,i)-
//                                          qi_h*(OvX->GetF2C(px_i,i))*OvX->Dx(Az,i,dx)
//                                  	  ) +
//                                  	0.5*mi_h*Dx_Sxz_ni_old -
//                                  (
//                                		  -qi_h*OvXold->Get(py_i,i)*OvXold->Get(Bx,i)-
//                                          qi_h*(OvX->GetF2C(px_i,i))*OvXold->Dx(Az,i,dx)
//                                  	  ))
//                                   - gamma[j + pz_i];

//          Electric field equaitons
//          Ex
//      res[j+Ax] = 		(0.5*dt*(OvX->Get(Ex,i) + OvXold->Get(Ex,i)) +
//              	  	  (OvX->Get(Ax,i) - OvXold->Get(Ax,i)));
//    //    Ey
//      res[j+Ey] = (0.5*dt*(OvX->Get(Ey,i) + OvXold->Get(Ey,i)) +
//                          (OvX->Get(Ay,i) - OvXold->Get(Ay,i)));
//      //    Ez
//      res[j+Ez] = (0.5*dt*(OvX->Get(Ez,i) + OvXold->Get(Ez,i)) +
//                          (OvX->Get(Az,i) - OvXold->Get(Az,i)));
//
//      //    Ay
//      res[j+Ex]     = 0.0*dt*zeta/(xi*xi)*(OvX->Dxx(Ax,i,dx)+OvXold->Dxx(Ax,i,dx))
//    				  + (dt*(qe_h*OvX->Get(px_e,i) + qi_h*OvX->Get(px_i,i))/(xi*xi)
//    				  + (OvX->Get(Ex,i) - OvXold->Get(Ex,i)));
//
//      //    Ay
//      res[j+Ay]     = 0.5*dt*zeta/(xi*xi)*(OvX->Dxx(Ay,i,dx)+OvXold->Dxx(Ay,i,dx))
//    				  + (dt*(qe_h*OvX->Get(py_e,i) + qi_h*OvX->Get(py_i,i))/(xi*xi)
//    				  + (OvX->Get(Ey,i) - OvXold->Get(Ey,i)));
//      //    Az
//      res[j+Az]     = 0.5*dt*zeta/(xi*xi)*(OvX->Dxx(Az,i,dx)+OvXold->Dxx(Az,i,dx))
//    				  + (dt*(qe_h*OvX->Get(pz_e,i) + qi_h*OvX->Get(pz_i,i))/(xi*xi)
//    				  + (OvX->Get(Ez,i) - OvXold->Get(Ez,i)));








      //    Ey
      res[j+Ey] = (0.5*dt*(OvX->Get(Ey,i) + OvXold->Get(Ey,i)) +
                          (OvX->Get(Ay,i) - OvXold->Get(Ay,i)));
      //    Ez
      res[j+Ez] = (0.5*dt*(OvX->Get(Ez,i) + OvXold->Get(Ez,i)) +
                          (OvX->Get(Az,i) - OvXold->Get(Az,i)));

      //    Ay
      res[j+Ex]     = (dt*(qe_h*OvX->Get(px_e,i) + qi_h*OvX->Get(px_i,i) - jmeanx1)/(xi*xi)
    				  + (OvX->Get(Ex,i) - OvXold->Get(Ex,i)));

      //    Ay
      res[j+Ay]     = (0.5*zeta*dt/(xi*xi)*(OvX->Dxx(Ay,i,dx)+OvXold->Dxx(Ay,i,dx))
    				  + (dt*(qe_h*OvX->Get(py_e,i) + qi_h*OvX->Get(py_i,i) - jmeany1)/(xi*xi)));
      //    Az
      res[j+Az]     = (0.5*zeta*dt/(xi*xi)*(OvX->Dxx(Az,i,dx)+OvXold->Dxx(Az,i,dx))
    				  + (dt*(qe_h*OvX->Get(pz_e,i) + qi_h*OvX->Get(pz_i,i) - jmeanz1)/(xi*xi)));


//      res[j+Ax] = 0.0;
//      res[j+Ey] = 0.0;
//      res[j+Ez] = 0.0;



//      res[j+n_e] = 0.0;
//      res[j+n_i] = 0.0;
//      res[j+px_e] = 0.0;
//      res[j+py_e] = 0.0;
//      res[j+pz_e] = 0.0;
//      res[j+px_i] = 0.0;
//      res[j+py_i] = 0.0;
//      res[j+pz_i] = 0.0;

//      res[j+Bx] = 0.0;
//      res[j+By] = 0.0;
//      res[j+Bz] = 0.0;




    }


}
