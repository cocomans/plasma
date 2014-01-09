#include "SimParams.h"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Ifpack_Utils.h"

////////////////////////////////////////
SimParams::SimParams(Epetra_Comm* Comm,PlasmaData* _pdata)
{

  comm = Comm;
  
  Phys = Teuchos::rcp(new Teuchos::ParameterList);
  Map   = Teuchos::rcp(new Teuchos::ParameterList);
  NLS  = Teuchos::rcp(new Teuchos::ParameterList);
  ML  = Teuchos::rcp(new Teuchos::ParameterList);
  Flags  = Teuchos::rcp(new Teuchos::ParameterList);

  pdata = _pdata;

  ImportParams();

}
////////////////////////////////////////

////////////////////////////////////////
SimParams::~SimParams(){}
////////////////////////////////////////

////////////////////////////////////////
void SimParams::ImportParams()
{

	int argc = pdata->argc_s;
	char** argv = pdata->argv_s;
    
  // Parse Command Line Options
  Trilinos_Util::CommandLineParser CLP(argc, argv);
  bool pre = CLP.Has("-pre");
  pre = false;

  int nx = pdata->nx;
  double x0 = pdata->xmin;
  double xf = pdata->Lx+pdata->xmin;

  int nUknwn;
  if((pdata->ndimensions == 1)&&(pdata->nVelocity == 1))
  {
	  using namespace ES1D1V;
	  nUknwn = Sxx_i + 1;
  }
  else if((pdata->ndimensions == 1)&&(pdata->nVelocity == 3))
  {
	  using namespace EM1D3V;
	  nUknwn = Sxz_i + 1;
  }
  if((pdata->ndimensions == 2)&&(pdata->nVelocity == 3))
  {
	  using namespace EM2D3V;
	  nUknwn = Szz_i + 1;
	  nx = pdata->nx * pdata->ny;
  }
  // Flags
  Flags->set("pre", pre);
  Flags->set("Output Directory", "data");

  // Set physics params

  Phys->set("x0", x0);
  Phys->set("xf", xf);

  Phys->set("dx", pdata->dxdi);
  Phys->set("dt", pdata->dt);
  Phys->set("xi", pdata->xi);
  Phys->set("zeta", pdata->bobzeta);
  Phys->set("me_h",pdata->m_h[0]);
  Phys->set("mi_h",pdata->m_h[1]);
  Phys->set("qe_h",pdata->q_h[0]);
  Phys->set("qi_h",pdata->q_h[1]);

  Phys->set("nSpatial",pdata->ndimensions);
  Phys->set("nVel",pdata->nVelocity);


  Phys->set("Newton Relative Tolerance", 1.0e-6);
  Phys->set("Newton Max Iterations", 20);

  // Map Parameters
  Map->set("Index", 0);
  Map->set("Number Elements", nx);
  Map->set("Element Overlap", 8); // must be >= 2
  Map->set("Number Unknowns", nUknwn); //23
  Map->set("Number Reduced Unknowns", 1);       // Must edit map manager array for

  ML_Epetra::SetDefaults("SA",*ML);
  // ML Parameters
  ML->set("ML output", 0);
  ML->set("cycle applications", 1); 

  ML->set("smoother: damping factor", 2.0/3.0); 
  ML->set("smoother: sweeps", 3); 
  ML->set("smoother: pre or post", "both"); 

  ML->set("coarse: type","Jacobi");
  ML->set("coarse: pre or post", "both"); 
  ML->set("coarse: sweeps", 250); 
  ML->set("coarse: damping factor", 2.0/3.0); 
  

  ////////////////////////////////////////////////////////////
  // Nonlinear Solver Params
  ////////////////////////////////////////////////////////////
  // Printing Params
  NLS->set("Nonlinear Solver", "Line Search Based");
  Teuchos::ParameterList& printParams = NLS->sublist("Printing");
  printParams.set("MyPID", comm->MyPID());
  printParams.set("Output Processor", 0);
  printParams.set("Output Information",
		  0
		  + NOX::Utils::OuterIteration
//		  + NOX::Utils::OuterIterationStatusTest
//		  + NOX::Utils::InnerIteration
//		  + NOX::Utils::Parameters
//		  + NOX::Utils::Details
//		  + NOX::Utils::Warning
		  ); 

  // Search Params
  Teuchos::ParameterList& searchParams = NLS->sublist("Line Search");
  searchParams.set("Method", "Full Step");

  // Direction Params
  Teuchos::ParameterList& dirParams = NLS->sublist("Direction");
  dirParams.set("Method", "Newton");

  // Newton Params
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Type 1");

  int iter_max = sqrt(nUknwn*nx);
  // Linear Solver Params
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 2*nUknwn*nx);
  lsParams.set("Size of Krylov Subspace",nUknwn*nx);
  lsParams.set("Tolerance", 1e-2);
//  lsParams.set("Preconditioner", "AZ_none");


}
////////////////////////////////////////

////////////////////////////////////////
Teuchos::RCP<Teuchos::ParameterList> SimParams::GetParamsPhys(){return Phys;}
Teuchos::RCP<Teuchos::ParameterList> SimParams::GetParamsMap(){return Map;}
Teuchos::RCP<Teuchos::ParameterList> SimParams::GetParamsNLS(){return NLS;}
Teuchos::RCP<Teuchos::ParameterList> SimParams::GetParamsML(){return ML;}
Teuchos::RCP<Teuchos::ParameterList> SimParams::GetParamsFlags(){return Flags;}

Teuchos::ParameterList& SimParams::GetPrintParams()
{
  return NLS->sublist("Printing");
}

Teuchos::ParameterList& SimParams::GetLSParams()
{
  Teuchos::ParameterList& dirParams = NLS->sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  return newtonParams.sublist("Linear Solver");
}

////////////////////////////////////////

////////////////////////////////////////
void SimParams::Print()
{
  LO_log << *NLS;
  LO_log << *Phys;
  LO_log << *Map;
  LO_log << *ML;
  LO_log << *Flags;
}
////////////////////////////////////////
