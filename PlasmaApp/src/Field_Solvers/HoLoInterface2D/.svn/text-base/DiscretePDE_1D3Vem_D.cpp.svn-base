#include "DiscretePDE_1D3Vem_D.h"

DiscretePDE_1D3Vem_D::DiscretePDE_1D3Vem_D(const Teuchos::RCP<MapManager> &map_manager,
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
void DiscretePDE_1D3Vem_D::ImportParams()
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

//int did_i_calc_means = 0;


////////////////////////////////////////////////////////////////////
// USER: This is where users must build the main residual calc.   //
// Propogate the enum below with names for each of your unknowns. //
// The number of these unknowns must match "Number Unknowns"      //
///////////////////////////////////////////////////////////////////
void DiscretePDE_1D3Vem_D::EvaluateResidual( const Epetra_Vector& x,
				    const Epetra_Vector& xold,
				    const Epetra_Vector& gamma,
				    Epetra_Vector& res )
{
  // Build an enum of variable names

	using namespace EM1D3V_D;

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

      jmeanx1 += qe_h*OvX->Get(px_e, i) + qi_h*OvX->Get(px_i,i);
      jmeany1 += qe_h*OvX->Get(py_e, i) + qi_h*OvX->Get(py_i,i);
      jmeanz1 += qe_h*OvX->Get(pz_e, i) + qi_h*OvX->Get(pz_i,i);



  }

  jmeanx1 /= (double)(N);
  jmeany1 /= (double)(N);
  jmeanz1 /= (double)(N);

//  printf("\n\n\n#############################################\n");
//  printf("Average Currents = %e %e %e\n",jmeanx1,jmeany1,jmeanz1);
//  printf("#############################################\n\n\n");

//  did_i_calc_means = 0;

  for(int i = 0; i<N; i++)
  {
      (*(OvX->X))[mapManager->LocElem_to_LocPtOv(i, Ex)] = OvXold->Get(Ex,i) - ((qe_h*OvX->Get(px_e,i) + qi_h*OvX->Get(px_i,i) - jmeanx1)*dt/(xi*xi));
      (*(OvX->X))[mapManager->LocElem_to_LocPtOv(i, Ey)] = -2.0*(OvX->Get(Ay,i) - OvXold->Get(Ay,i))/dt - OvXold->Get(Ey,i);
      (*(OvX->X))[mapManager->LocElem_to_LocPtOv(i, Ez)] = -2.0*(OvX->Get(Az,i) - OvXold->Get(Az,i))/dt - OvXold->Get(Ez,i);

  }


#pragma omp parallel for
  for(int i = 0; i<N; i++)
    {   
      // Convert elem index to first point index
      int j = mapManager->LocElem_to_LocPtSd(i,0);

      //    Electric field equaitons
      //    Ex
//      res[j+Ax] = 		(0.5*dt*(OvX->Get(Ex,i) + OvXold->Get(Ex,i)) +
//              	  	  (OvX->Get(Ax,i) - OvXold->Get(Ax,i)));
//      //     Ey
//      res[j+Ey] = (0.5*dt*(OvX->Get(Ey,i) + OvXold->Get(Ey,i)) +
//                          (OvX->Get(Ay,i) - OvXold->Get(Ay,i)));
//      //    Ez
//      res[j+Ez] = (0.5*dt*(OvX->Get(Ez,i) + OvXold->Get(Ez,i)) +
//                          (OvX->Get(Az,i) - OvXold->Get(Az,i)));
//
//      //    Ay
//      res[j+Ex]     = (dt*(qe_h*OvX->Get(px_e,i) + qi_h*OvX->Get(px_i,i) - jmeanx1)/(xi*xi)
//    				  + (OvX->Get(Ex,i) - OvXold->Get(Ex,i)));

      //    Ay
      res[j+Ay]     = dt*(OvX->Dxx(Ay,i,dx)+OvXold->Dxx(Ay,i,dx))
    				  + (2.0*dt/zeta*(qe_h*OvX->Get(py_e,i) + qi_h*OvX->Get(py_i,i) - jmeany1)
    				  + 0.0*(OvX->Get(Ey,i) - OvXold->Get(Ey,i)));
      //    Az
      res[j+Az]     = dt*(OvX->Dxx(Az,i,dx)+OvXold->Dxx(Az,i,dx))
    				  + (2.0*dt/zeta*(qe_h*OvX->Get(pz_e,i) + qi_h*OvX->Get(pz_i,i) - jmeanz1)
    				  + 0.0*(OvX->Get(Ez,i) - OvXold->Get(Ez,i)));


      res[j+Ax] = 0.0;
      res[j+Ex] = (x[j+Ex] - OvX->Get(Ex,i));
      res[j+Ey] = (x[j+Ey] - OvX->Get(Ey,i));
      res[j+Ez] = (x[j+Ez] - OvX->Get(Ez,i));




//      res[j+n_e] = 0.0;
//      res[j+n_i] = 0.0;
//      res[j+px_e] = 0.0;
//      res[j+py_e] = 0.0;
//      res[j+pz_e] = 0.0;
//      res[j+px_i] = 0.0;
//      res[j+py_i] = 0.0;
//      res[j+pz_i] = 0.0;

      res[j+Bx] = 0.0;
      res[j+By] = 0.0;
      res[j+Bz] = 0.0;




    }


}
