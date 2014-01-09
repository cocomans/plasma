#include "DiscretePDE2D_1D3Vem.h"

DiscretePDE2D_1D3Vem::DiscretePDE2D_1D3Vem(const Teuchos::RCP<SimParams> &_params,
	      Epetra_Comm* _comm)
{
  comm       = _comm;
  simParams  = _params;
//  mapManager = map_manager;
//
//  stagMesh   = Teuchos::rcp(new StagMesh(mapManager, simParams, comm) );
//  OvX        = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );
//  OvXold     = Teuchos::rcp(new DiffMath(mapManager, simParams, comm) );

  ImportParams();
};

//DiscretePDE::~DiscretePDE(){}

////////////////////////////////////////////////////////////////////
// USER: Edit any constants you may need in SimParams.cpp  //
///////////////////////////////////////////////////////////////////
void DiscretePDE2D_1D3Vem::ImportParams()
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
void DiscretePDE2D_1D3Vem::EvaluateResidual(EpVecWrapper& x,
				    const EpVecWrapper& xold,
				    const EpVecWrapper& gamma,
				    EpVecWrapper& res )
{
  // Build an enum of variable names

	using namespace EM1D3V;

//  std::cout << "gamma = " << std::endl;
//  std::cout << gamma;



  //Global Operations (ie. averages)
  double jmeanx1 = qe_h*x.mean(px_e) + qi_h*x.mean(px_i);
  double jmeany1 = qe_h*x.mean(py_e) + qi_h*x.mean(py_i);
  double jmeanz1 = qe_h*x.mean(pz_e) + qi_h*x.mean(pz_i);

//  for(int i = 0; i<N; i++)
//  {
//
//	  x.Get(Ey,i) = -2.0*(x.Get(Ay,i) - xold.Get(Ay,i))/dt - xold.Get(Ey,i);
//	  x.Get(Ez,i) = -2.0*(x.Get(Az,i) - xold.Get(Az,i))/dt - xold.Get(Ez,i);
//
//
//  }


#pragma omp parallel for
for(int i = 0; i<x.nx; i++)
    {   
      // Convert elem index to first point index

//      printf("B[%i] = %e %e %e\n",i,x.Get(Bx,i),x.Get(By,i),x.Get(Bz,i));
//      printf("E[%i] = %e %e %e\n",i,x.Get(Ex,i),x.Get(Ey,i),x.Get(Ez,i));
//      printf("A[%i] = %e %e %e\n",i,x.Get(Ax,i),x.Get(Ay,i),x.Get(Az,i));
//      printf("pe[%i] = %e %e %e\n",i,x.Get(px_e,i),x.Get(py_e,i),x.Get(pz_e,i));
//      printf("pi[%i] = %e %e %e\n",i,x.Get(px_i,i),x.Get(py_i,i),x.Get(pz_i,i));


//	  double Ext = xold.Get(Ex,i) - dt*(qe_h*x.Get(px_e, i) + qi_h*x.Get(px_i,i) - jmeanx1)/(xi*xi);
	  double Eyt = -2.0*(x.Get(Ay,i) - xold.Get(Ay,i))/dt - xold.Get(Ey,i);
	  double Ezt = -2.0*(x.Get(Az,i) - xold.Get(Az,i))/dt - xold.Get(Ez,i);

	  double Ext = x.Get(Ex,i);
//	  double Eyt = x.Get(Ey,i);
//	  double Ezt = x.Get(Ez,i);

//      double tmp_W = (x.Get(Sxx_e, i-1)+xold.Get(Sxx_e, i-1))*(x.Get(n_e,i-1)+xold.Get(n_e,i-1));
//      double tmp_C = (x.Get(Sxx_e, i)+xold.Get(Sxx_e, i)) * (x.Get(n_e,i)+xold.Get(n_e,i));
//      double tmp_E = (x.Get(Sxx_e, i+1)+xold.Get(Sxx_e, i+1)) * (x.Get(n_e,i+1)+xold.Get(n_e,i+1));
//      double Dx_Sxx_ne = 0.25*x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = 0.25*(x.Get(Sxy_e, i-1)+xold.Get(Sxy_e, i-1)) * (x.Get(n_e,i-1)+xold.Get(n_e,i-1));
//      tmp_C = 0.25*(x.Get(Sxy_e, i)+xold.Get(Sxy_e, i)) * (x.Get(n_e,i)+xold.Get(n_e,i));
//      tmp_E = 0.25*(x.Get(Sxy_e, i+1)+xold.Get(Sxy_e, i+1)) * (x.Get(n_e,i+1)+xold.Get(n_e,i+1));
//      double Dx_Sxy_ne = x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = 0.25*(x.Get(Sxz_e, i-1)+xold.Get(Sxz_e, i-1)) * (x.Get(n_e,i-1)+xold.Get(n_e,i-1));
//      tmp_C = 0.25*(x.Get(Sxz_e, i)+xold.Get(Sxz_e, i)) * (x.Get(n_e,i)+xold.Get(n_e,i));
//      tmp_E = 0.25*(x.Get(Sxz_e, i+1)+xold.Get(Sxz_e, i+1)) * (x.Get(n_e,i+1)+xold.Get(n_e,i+1));
//      double Dx_Sxz_ne = x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      // Ion
//      tmp_W = 0.25*(x.Get(Sxx_i, i-1)+xold.Get(Sxx_i, i-1)) * (x.Get(n_i,i-1)+xold.Get(n_i,i-1));
//      tmp_C = 0.25*(x.Get(Sxx_i, i)+xold.Get(Sxx_i, i)) * (x.Get(n_i,i)+xold.Get(n_i,i));
//      tmp_E = 0.25*(x.Get(Sxx_i, i+1)+xold.Get(Sxx_i, i+1)) * (x.Get(n_i,i+1)+xold.Get(n_i,i+1));
//      double Dx_Sxx_ni = x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = 0.25*(x.Get(Sxy_i, i-1)+xold.Get(Sxy_i, i-1)) * (x.Get(n_i,i-1)+xold.Get(n_i,i-1));
//      tmp_C = 0.25*(x.Get(Sxy_i, i)+xold.Get(Sxy_i, i)) * (x.Get(n_i,i)+xold.Get(n_i,i));
//      tmp_E = 0.25*(x.Get(Sxy_i, i+1)+xold.Get(Sxy_i, i+1)) * (x.Get(n_i,i+1)+xold.Get(n_i,i+1));
//      double Dx_Sxy_ni = x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = 0.25*(x.Get(Sxz_i, i-1)+xold.Get(Sxz_i, i-1)) * (x.Get(n_i,i-1)+xold.Get(n_i,i-1));
//      tmp_C = 0.25*(x.Get(Sxz_i, i)+xold.Get(Sxz_i, i)) * (x.Get(n_i,i)+xold.Get(n_i,i));
//      tmp_E = 0.25*(x.Get(Sxz_i, i+1)+xold.Get(Sxz_i, i+1)) * (x.Get(n_i,i+1)+xold.Get(n_i,i+1));
//      double Dx_Sxz_ni = x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//
//      double tmp_W = (x.Get(Sxx_e, i-1)+3.0*xold.Get(Sxx_e, i-1))*(x.Get(n_e,i-1)+3.0*xold.Get(n_e,i-1));
//      double tmp_C = (x.Get(Sxx_e, i)+3.0*xold.Get(Sxx_e, i)) * (x.Get(n_e,i)+3.0*xold.Get(n_e,i));
//      double tmp_E = (x.Get(Sxx_e, i+1)+3.0*xold.Get(Sxx_e, i+1)) * (x.Get(n_e,i+1)+3.0*xold.Get(n_e,i+1));
//      double Dx_Sxx_ne = 0.0625*x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = (x.Get(Sxy_e, i-1)+3.*xold.Get(Sxy_e, i-1)) * (x.Get(n_e,i-1)+3.*xold.Get(n_e,i-1));
//      tmp_C = (x.Get(Sxy_e, i)+3.*xold.Get(Sxy_e, i)) * (x.Get(n_e,i)+3.*xold.Get(n_e,i));
//      tmp_E = (x.Get(Sxy_e, i+1)+3.*xold.Get(Sxy_e, i+1)) * (x.Get(n_e,i+1)+3.*xold.Get(n_e,i+1));
//      double Dx_Sxy_ne = 0.0625*x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = (x.Get(Sxz_e, i-1)+3.*xold.Get(Sxz_e, i-1)) * (x.Get(n_e,i-1)+3.*xold.Get(n_e,i-1));
//      tmp_C = (x.Get(Sxz_e, i)+3.*xold.Get(Sxz_e, i)) * (x.Get(n_e,i)+3.*xold.Get(n_e,i));
//      tmp_E = (x.Get(Sxz_e, i+1)+3.*xold.Get(Sxz_e, i+1)) * (x.Get(n_e,i+1)+3.*xold.Get(n_e,i+1));
//      double Dx_Sxz_ne = 0.0625*x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      // Ion
//      tmp_W = (x.Get(Sxx_i, i-1)+3.*xold.Get(Sxx_i, i-1)) * (x.Get(n_i,i-1)+3.*xold.Get(n_i,i-1));
//      tmp_C = (x.Get(Sxx_i, i)+3.*xold.Get(Sxx_i, i)) * (x.Get(n_i,i)+3.*xold.Get(n_i,i));
//      tmp_E = (x.Get(Sxx_i, i+1)+3.*xold.Get(Sxx_i, i+1)) * (x.Get(n_i,i+1)+3.*xold.Get(n_i,i+1));
//      double Dx_Sxx_ni = 0.0625*x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = (x.Get(Sxy_i, i-1)+3.*xold.Get(Sxy_i, i-1)) * (x.Get(n_i,i-1)+3.*xold.Get(n_i,i-1));
//      tmp_C = (x.Get(Sxy_i, i)+3.*xold.Get(Sxy_i, i)) * (x.Get(n_i,i)+3.*xold.Get(n_i,i));
//      tmp_E = (x.Get(Sxy_i, i+1)+3.*xold.Get(Sxy_i, i+1)) * (x.Get(n_i,i+1)+3.*xold.Get(n_i,i+1));
//      double Dx_Sxy_ni = 0.0625*x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//      tmp_W = (x.Get(Sxz_i, i-1)+3.*xold.Get(Sxz_i, i-1)) * (x.Get(n_i,i-1)+3.*xold.Get(n_i,i-1));
//      tmp_C = (x.Get(Sxz_i, i)+3.*xold.Get(Sxz_i, i)) * (x.Get(n_i,i)+3.*xold.Get(n_i,i));
//      tmp_E = (x.Get(Sxz_i, i+1)+3.*xold.Get(Sxz_i, i+1)) * (x.Get(n_i,i+1)+3.*xold.Get(n_i,i+1));
//      double Dx_Sxz_ni = 0.0625*x.DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//
//
//      // Fill each unknown
//      //    Electron continuity equation
//      res[j+n_e] = (x.Get(n_e, i) - xold.Get(n_e, i)) + dt*( x.DxF2C(px_e, i, dx)) - gamma[j + n_e];
//      //    Ion continuity equation
//      res[j+n_i] = (x.Get(n_i, i) - xold.Get(n_i, i)) + dt*( x.DxF2C(px_i, i, dx) ) - gamma[j + n_i];
//
//      //    Electron momentum equation
//      //    momentum electron-x
//      res[j+px_e] = (x.Get(px_e, i) - xold.Get(px_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxx_ne -
//                                  0.0625*qe_h*(x.GetC2F(n_e,i)+3.*xold.GetC2F(n_e,i))*(3.*xold.Get(Ex,i)+Ext) -
//                                  (
//                                          0.125*qe_h*(x.GetC2F(py_e,i)+xold.GetC2F(py_e,i))*(x.DxC2F(Ay,i,dx) + 3.*xold.DxC2F(Ay,i,dx))+
//                                          0.125*qe_h*(x.GetC2F(pz_e,i)+xold.GetC2F(pz_e,i))*(x.DxC2F(Az,i,dx) + 3.*xold.DxC2F(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + px_e];
//
//      //    momentum electron-y
//      res[j+py_e] = (x.Get(py_e, i) - xold.Get(py_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxy_ne -
//                                  0.0625*qe_h*(x.Get(n_e,i)+3.*xold.Get(n_e,i))*(3.*xold.Get(Ey,i)+Eyt) -
//                                  (
//                                		  qe_h*x.Get(pz_e,i)*xold.Get(Bx,i)-
//                                          0.125*qe_h*(x.GetF2C(px_e,i)+xold.GetF2C(px_e,i))*(x.Dx(Ay,i,dx) + 3.*xold.Dx(Ay,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + py_e];
//
//      //    momentum electron-z
//      res[j+pz_e] = (x.Get(pz_e, i) - xold.Get(pz_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxz_ne -
//                                  0.0625*qe_h*(x.Get(n_e,i)+3.*xold.Get(n_e,i))*(3.*xold.Get(Ez,i)+Ezt) -
//                                  (
//                                		  -qe_h*x.Get(py_e,i)*xold.Get(Bx,i)-
//                                          0.125*qe_h*(x.GetF2C(px_e,i)+xold.GetF2C(px_e,i))*(x.Dx(Az,i,dx) + 3.*xold.Dx(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + pz_e];
//
//
//
//      //    Ion momentum equation
//      //    momentum ion-x
//      res[j+px_i] = (x.Get(px_i, i) - xold.Get(px_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxx_ni -
//                                  0.0625*qi_h*(x.GetC2F(n_i,i)+3.*xold.GetC2F(n_i,i))*(3.*xold.Get(Ex,i)+Ext) -
//                                  (
//                                          0.125*qi_h*(x.GetC2F(py_i,i)+xold.GetC2F(py_i,i))*(x.DxC2F(Ay,i,dx) + 3.*xold.DxC2F(Ay,i,dx)) +
//                                          0.125*qi_h*(x.GetC2F(pz_i,i)+xold.GetC2F(pz_i,i))*(x.DxC2F(Az,i,dx) + 3.*xold.DxC2F(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + px_i];
//
//      //    momentum electron-y
//      res[j+py_i] = (x.Get(py_i, i) - xold.Get(py_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxy_ni -
//                                  0.0625*qi_h*(x.Get(n_i,i)+3.*xold.Get(n_i,i))*(3.*xold.Get(Ey,i)+Eyt) -
//                                  (
//                                		  qi_h*x.Get(pz_i,i)*xold.Get(Bx,i)-
//                                          0.125*qi_h*(x.GetF2C(px_i,i)+xold.GetF2C(px_i,i))*(x.Dx(Ay,i,dx) + 3.*xold.Dx(Ay,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + py_i];
//
//      //    momentum electron-z
//      res[j+pz_i] = (x.Get(pz_i, i) - xold.Get(pz_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxz_ni -
//                                  0.0625*qi_h*(x.Get(n_i,i)+3.*xold.Get(n_i,i))*(3.*xold.Get(Ez,i)+Ezt) -
//                                  (
//                                		  -qi_h*x.Get(py_i,i)*xold.Get(Bx,i)-
//                                          0.125*qi_h*(x.GetF2C(px_i,i)+xold.GetF2C(px_i,i))*(x.Dx(Az,i,dx) + 3.*xold.Dx(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + pz_i];


//
//      // Fill each unknown
//      //    Electron continuity equation
//      res[j+n_e] = (x.Get(n_e, i) - xold.Get(n_e, i)) + dt*( x.DxF2C(px_e, i, dx)) - gamma[j + n_e];
//      //    Ion continuity equation
//      res[j+n_i] = (x.Get(n_i, i) - xold.Get(n_i, i)) + dt*( x.DxF2C(px_i, i, dx) ) - gamma[j + n_i];
//
//      //    Electron momentum equation
//      //    momentum electron-x
//      res[j+px_e] = (x.Get(px_e, i) - xold.Get(px_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxx_ne -
//                                  0.25*(qe_h*x.GetC2F(n_e,i)+qe_h*xold.GetC2F(n_e,i))*(xold.Get(Ex,i)+x.Get(Ex,i)) -
//                                  (
//                                          0.5*qe_h*x.GetC2F(py_e,i)*(x.DxC2F(Ay,i,dx) + xold.DxC2F(Ay,i,dx))+
//                                          0.5*qe_h*x.GetC2F(pz_e,i)*(x.DxC2F(Az,i,dx) + xold.DxC2F(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + px_e];
//
//      //    momentum electron-y
//      res[j+py_e] = (x.Get(py_e, i) - xold.Get(py_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxy_ne -
//                                  0.25*(qe_h*x.Get(n_e,i)+qe_h*xold.Get(n_e,i))*(xold.Get(Ey,i)+x.Get(Ey,i)) -
//                                  (
//                                		  qe_h*x.Get(pz_e,i)*xold.Get(Bx,i)-
//                                          0.5*qe_h*x.GetF2C(px_e,i)*(x.Dx(Ay,i,dx) + xold.Dx(Ay,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + py_e];
//
//      //    momentum electron-z
//      res[j+pz_e] = (x.Get(pz_e, i) - xold.Get(pz_e, i)) +
//                          (dt/(2.0*me_h))*(
//                                  me_h*Dx_Sxz_ne -
//                                  0.25*(qe_h*x.Get(n_e,i)+qe_h*xold.Get(n_e,i))*(xold.Get(Ez,i)+x.Get(Ez,i)) -
//                                  (
//                                		  -qe_h*x.Get(py_e,i)*xold.Get(Bx,i)-
//                                          0.5*qe_h*x.GetF2C(px_e,i)*(x.Dx(Az,i,dx) + xold.Dx(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + pz_e];
//
//
//
//      //    Ion momentum equation
//      //    momentum ion-x
//      res[j+px_i] = (x.Get(px_i, i) - xold.Get(px_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxx_ni -
//                                  0.25*(qi_h*x.GetC2F(n_i,i)+qi_h*xold.GetC2F(n_i,i))*(xold.Get(Ex,i)+x.Get(Ex,i)) -
//                                  (
//                                          0.5*qi_h*x.GetC2F(py_i,i)*(x.DxC2F(Ay,i,dx) + xold.DxC2F(Ay,i,dx)) +
//                                          0.5*qi_h*x.GetC2F(pz_i,i)*(x.DxC2F(Az,i,dx) + xold.DxC2F(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + px_i];
//
//      //    momentum electron-y
//      res[j+py_i] = (x.Get(py_i, i) - xold.Get(py_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxy_ni -
//                                  0.25*(qi_h*x.Get(n_i,i)+qi_h*xold.Get(n_i,i))*(xold.Get(Ey,i)+x.Get(Ey,i)) -
//                                  (
//                                		  qi_h*x.Get(pz_i,i)*xold.Get(Bx,i)-
//                                          0.5*qi_h*x.GetF2C(px_i,i)*(x.Dx(Ay,i,dx) + xold.Dx(Ay,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + py_i];
//
//      //    momentum electron-z
//      res[j+pz_i] = (x.Get(pz_i, i) - xold.Get(pz_i, i)) +
//                          (dt/(2.0*mi_h))*(
//                                  mi_h*Dx_Sxz_ni -
//                                  0.25*(qi_h*x.Get(n_i,i)+qi_h*xold.Get(n_i,i))*(xold.Get(Ez,i)+x.Get(Ez,i)) -
//                                  (
//                                		  -qi_h*x.Get(py_i,i)*xold.Get(Bx,i)-
//                                          0.5*qi_h*x.GetF2C(px_i,i)*(x.Dx(Az,i,dx) + xold.Dx(Az,i,dx))
//                                  	  )
//                                  )
//                                   - gamma[j + pz_i];


      //    Construct nonlinear terms
      //    Electron
      //    Compute advection for momentum x
      double tmp_W = x.Get(Sxx_e, i-1) * x.Get(n_e,i-1);
      double tmp_C = x.Get(Sxx_e, i) * x.Get(n_e,i);
      double tmp_E = x.Get(Sxx_e, i+1) * x.Get(n_e,i+1);
      double Dx_Sxx_ne = x.DxC2F(tmp_W, tmp_C, tmp_E);

      tmp_W = xold.Get(Sxx_e, i-1) * xold.Get(n_e,i-1);
      tmp_C = xold.Get(Sxx_e, i) * xold.Get(n_e,i);
      tmp_E = xold.Get(Sxx_e, i+1) * xold.Get(n_e,i+1);
      double Dx_Sxx_ne_old = xold.DxC2F(tmp_W, tmp_C, tmp_E);
      //    Compute advection for momentum y
      tmp_W = x.Get(Sxy_e, i-1) * x.Get(n_e,i-1);
      tmp_C = x.Get(Sxy_e, i) * x.Get(n_e,i);
      tmp_E = x.Get(Sxy_e, i+1) * x.Get(n_e,i+1);
      double Dx_Sxy_ne = x.DxC2F(tmp_W, tmp_C, tmp_E);

      tmp_W = xold.Get(Sxy_e, i-1) * xold.Get(n_e,i-1);
      tmp_C = xold.Get(Sxy_e, i) * xold.Get(n_e,i);
      tmp_E = xold.Get(Sxy_e, i+1) * xold.Get(n_e,i+1);
      double Dx_Sxy_ne_old = xold.DxC2F(tmp_W, tmp_C, tmp_E);
      //    Compute advection for momentum z
      tmp_W = x.Get(Sxz_e, i-1) * x.Get(n_e,i-1);
      tmp_C = x.Get(Sxz_e, i) * x.Get(n_e,i);
      tmp_E = x.Get(Sxz_e, i+1) * x.Get(n_e,i+1);
      double Dx_Sxz_ne = x.DxC2F(tmp_W, tmp_C, tmp_E);

      tmp_W = xold.Get(Sxz_e, i-1) * xold.Get(n_e,i-1);
      tmp_C = xold.Get(Sxz_e, i) * xold.Get(n_e,i);
      tmp_E = xold.Get(Sxz_e, i+1) * xold.Get(n_e,i+1);
      double Dx_Sxz_ne_old = xold.DxC2F(tmp_W, tmp_C, tmp_E);

      //    ion
      //    Compute advection for momentum x
      tmp_W = x.Get(Sxx_i, i-1) * x.Get(n_i,i-1);
      tmp_C = x.Get(Sxx_i, i) * x.Get(n_i,i);
      tmp_E = x.Get(Sxx_i, i+1) * x.Get(n_i,i+1);
      double Dx_Sxx_ni = x.DxC2F(tmp_W, tmp_C, tmp_E);

      tmp_W = xold.Get(Sxx_i, i-1) * xold.Get(n_i,i-1);
      tmp_C = xold.Get(Sxx_i, i) * xold.Get(n_i,i);
      tmp_E = xold.Get(Sxx_i, i+1) * xold.Get(n_i,i+1);
      double Dx_Sxx_ni_old = xold.DxC2F(tmp_W, tmp_C, tmp_E);
      //    Compute advection for momentum y
      tmp_W = x.Get(Sxy_i, i-1) * x.Get(n_i,i-1);
      tmp_C = x.Get(Sxy_i, i) * x.Get(n_i,i);
      tmp_E = x.Get(Sxy_i, i+1) * x.Get(n_i,i+1);
      double Dx_Sxy_ni = x.DxC2F(tmp_W, tmp_C, tmp_E);

      tmp_W = xold.Get(Sxy_i, i-1) * xold.Get(n_i,i-1);
      tmp_C = xold.Get(Sxy_i, i) * xold.Get(n_i,i);
      tmp_E = xold.Get(Sxy_i, i+1) * xold.Get(n_i,i+1);
      double Dx_Sxy_ni_old = xold.DxC2F(tmp_W, tmp_C, tmp_E);
      //    Compute advection for momentum z
      tmp_W = x.Get(Sxz_i, i-1) * x.Get(n_i,i-1);
      tmp_C = x.Get(Sxz_i, i) * x.Get(n_i,i);
      tmp_E = x.Get(Sxz_i, i+1) * x.Get(n_i,i+1);
      double Dx_Sxz_ni = x.DxC2F(tmp_W, tmp_C, tmp_E);

      tmp_W = xold.Get(Sxz_i, i-1) * xold.Get(n_i,i-1);
      tmp_C = xold.Get(Sxz_i, i) * xold.Get(n_i,i);
      tmp_E = xold.Get(Sxz_i, i+1) * xold.Get(n_i,i+1);
      double Dx_Sxz_ni_old = xold.DxC2F(tmp_W, tmp_C, tmp_E);

//      // Fill each unknown
//      //    Electron continuity equation
      res(n_e,i) = (x.Get(n_e, i) - xold.Get(n_e, i)) + dt*( x.DxF2C(px_e, i)) - gamma(n_e,i);
//      //    Ion continuity equation
      res(n_i,i) = (x.Get(n_i, i) - xold.Get(n_i, i)) + dt*( x.DxF2C(px_i, i)) - gamma(n_i,i);

      //    Electron momentum equation
      //    momentum electron-x
      res(px_e,i) = (x.Get(px_e, i) - xold.Get(px_e, i)) +
                          (dt/(2.0*me_h))*(
                                  me_h*Dx_Sxx_ne -
                                  0.5*(qe_h*x.GetC2F(n_e,i)+qe_h*xold.GetC2F(n_e,i))*(xold.Get(Ex,i)+Ext) -
                                  (
                                          qe_h*(x.GetC2F(py_e,i))*x.DxC2F(Ay,i) +
                                          qe_h*(x.GetC2F(pz_e,i))*x.DxC2F(Az,i)
                                  	  ) +
                                  me_h*Dx_Sxx_ne_old -
                                  (
                                          qe_h*(x.GetC2F(py_e,i))*xold.DxC2F(Ay,i) +
                                          qe_h*(x.GetC2F(pz_e,i))*xold.DxC2F(Az,i)
                                  	  ))
                                   - gamma(px_e,i);

      //    momentum electron-y
      res(py_e,i) = (x.Get(py_e, i) - xold.Get(py_e, i)) +
                          (dt/(2.0*me_h))*(
                                  me_h*Dx_Sxy_ne -
                                  0.5*(qe_h*x.Get(n_e,i)+qe_h*xold.Get(n_e,i))*(xold.Get(Ey,i)+Eyt) -
                                  (
                                		  qe_h*x.Get(pz_e,i)*xold.Get(Bx,i)-
                                          qe_h*(x.GetF2C(px_e,i))*x.Dx(Ay,i)
                                  	  ) +
                                  me_h*Dx_Sxy_ne_old -
                                  (
                                		  qe_h*xold.Get(pz_e,i)*xold.Get(Bx,i)-
                                		  qe_h*(x.GetF2C(px_e,i))*xold.Dx(Ay,i)
                                  	  ))
                                   - gamma(py_e,i);

      //    momentum electron-z
      res(pz_e,i) = (x.Get(pz_e, i) - xold.Get(pz_e, i)) +
                          (dt/(2.0*me_h))*(
                                  me_h*Dx_Sxz_ne -
                                  0.5*(qe_h*x.Get(n_e,i)+qe_h*xold.Get(n_e,i))*(xold.Get(Ez,i)+Ezt) -
                                  (
                                		  -qe_h*x.Get(py_e,i)*xold.Get(Bx,i)-
                                          qe_h*(x.GetF2C(px_e,i))*x.Dx(Az,i)
                                  	  ) +
                                  me_h*Dx_Sxz_ne_old -
                                  (
                                		  -qe_h*xold.Get(py_e,i)*xold.Get(Bx,i)-
                                          qe_h*(x.GetF2C(px_e,i))*xold.Dx(Az,i)
                                  	  ))
                                   - gamma(pz_e,i);



      //    Ion momentum equation
      //    momentum ion-x
      res(px_i,i) = (x.Get(px_i, i) - xold.Get(px_i, i)) +
                          (dt/(2.0*mi_h))*(
                                  mi_h*Dx_Sxx_ni -
                                  0.5*(qi_h*x.GetC2F(n_i,i)+qi_h*xold.GetC2F(n_i,i))*(xold.Get(Ex,i)+Ext) -
                                  (
                                		  qi_h*(x.GetC2F(py_i,i))*x.DxC2F(Ay,i) +
                                		  qi_h*(x.GetC2F(pz_i,i))*x.DxC2F(Az,i)
                                  	  ) +
                                  mi_h*Dx_Sxx_ni_old -

                                  (
                                          qi_h*(x.GetC2F(py_i,i))*xold.DxC2F(Ay,i) +
                                          qi_h*(x.GetC2F(pz_i,i))*xold.DxC2F(Az,i)
                                  	  ))
                                   - gamma(px_i,i);

      //    momentum electron-y
      res(py_i,i) = (x.Get(py_i, i) - xold.Get(py_i, i)) +
                          (dt/(2.0*mi_h))*(
                                  mi_h*Dx_Sxy_ni -
                                  0.5*(qi_h*x.Get(n_i,i)+qi_h*xold.Get(n_i,i))*(xold.Get(Ey,i)+Eyt) -
                                  (
                                		  qi_h*x.Get(pz_i,i)*xold.Get(Bx,i)-
                                          qi_h*(x.GetF2C(px_i,i))*x.Dx(Ay,i)
                                  	  ) +
                                  mi_h*Dx_Sxy_ni_old -
                                  (
                                		  qi_h*xold.Get(pz_i,i)*xold.Get(Bx,i)-
                                		  qi_h*(x.GetF2C(px_i,i))*xold.Dx(Ay,i)
                                  	  ))
                                   - gamma(py_i,i);

      //    momentum electron-z
      res(pz_i,i) = (x.Get(pz_i, i) - xold.Get(pz_i, i)) +
                          (dt/(2.0*mi_h))*(
                                  mi_h*Dx_Sxz_ni -
                                  0.5*(qi_h*x.Get(n_i,i)+qi_h*xold.Get(n_i,i))*(xold.Get(Ez,i)+Ezt) -
                                  (
                                		  -qi_h*x.Get(py_i,i)*xold.Get(Bx,i)-
                                          qi_h*(x.GetF2C(px_i,i))*x.Dx(Az,i)
                                  	  ) +
                                  mi_h*Dx_Sxz_ni_old -
                                  (
                                		  -qi_h*xold.Get(py_i,i)*xold.Get(Bx,i)-
                                          qi_h*(x.GetF2C(px_i,i))*xold.Dx(Az,i)
                                  	  ))
                                   - gamma(pz_i,i);

//          Electric field equaitons
//          Ex
//      res[j+Ax] = 		(0.5*dt*(x.Get(Ex,i) + xold.Get(Ex,i)) +
//              	  	  (x.Get(Ax,i) - xold.Get(Ax,i)));
//    //    Ey
//      res[j+Ey] = (0.5*dt*(x.Get(Ey,i) + xold.Get(Ey,i)) +
//                          (x.Get(Ay,i) - xold.Get(Ay,i)));
//      //    Ez
//      res[j+Ez] = (0.5*dt*(x.Get(Ez,i) + xold.Get(Ez,i)) +
//                          (x.Get(Az,i) - xold.Get(Az,i)));
//
//      //    Ay
//      res[j+Ex]     = 0.0*dt*zeta/(xi*xi)*(x.Dxx(Ax,i,dx)+xold.Dxx(Ax,i,dx))
//    				  + (dt*(qe_h*x.Get(px_e,i) + qi_h*x.Get(px_i,i))/(xi*xi)
//    				  + (x.Get(Ex,i) - xold.Get(Ex,i)));
//
//      //    Ay
//      res[j+Ay]     = 0.5*dt*zeta/(xi*xi)*(x.Dxx(Ay,i,dx)+xold.Dxx(Ay,i,dx))
//    				  + (dt*(qe_h*x.Get(py_e,i) + qi_h*x.Get(py_i,i))/(xi*xi)
//    				  + (x.Get(Ey,i) - xold.Get(Ey,i)));
//      //    Az
//      res[j+Az]     = 0.5*dt*zeta/(xi*xi)*(x.Dxx(Az,i,dx)+xold.Dxx(Az,i,dx))
//    				  + (dt*(qe_h*x.Get(pz_e,i) + qi_h*x.Get(pz_i,i))/(xi*xi)
//    				  + (x.Get(Ez,i) - xold.Get(Ez,i)));








      //    Ey
      res(Ey,i) = (0.5*dt*(x.Get(Ey,i) + xold.Get(Ey,i)) +
                          (x.Get(Ay,i) - xold.Get(Ay,i)));
      //    Ez
      res(Ez,i) = (0.5*dt*(x.Get(Ez,i) + xold.Get(Ez,i)) +
                          (x.Get(Az,i) - xold.Get(Az,i)));

      //    Ay
      res(Ex,i)     = (dt*(qe_h*x.Get(px_e,i) + qi_h*x.Get(px_i,i) - jmeanx1)/(xi*xi)
    				  + (x.Get(Ex,i) - xold.Get(Ex,i)));

      //    Ay
      res(Ay,i)     = (0.5*zeta*dt/(xi*xi)*(x.Dxx(Ay,i)+xold.Dxx(Ay,i))
    				  + (dt*(qe_h*x.Get(py_e,i) + qi_h*x.Get(py_i,i) - jmeany1)/(xi*xi)));
      //    Az
      res(Az,i)     = (0.5*zeta*dt/(xi*xi)*(x.Dxx(Az,i)+xold.Dxx(Az,i))
    				  + (dt*(qe_h*x.Get(pz_e,i) + qi_h*x.Get(pz_i,i) - jmeanz1)/(xi*xi)));


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
