#ifndef TIMESTEPPER_H_
#define TIMESTEPPER_H_

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_Time.hpp"

// User Includes
#include "SimParams.h"
#include <NOX_Epetra.H>

#include "JFNK.h"
#include "SimParams.h"
#include "LowOrderProblem.h"

#include "DiffMath.h"

//=========================
//Include HO header-files
//========================
#include "../../HOMomentsCPU.h"
#include "../../FieldDataCPU.h"
#include "../../NodeHOMoments.h"
#include "../../PlasmaData.h"
#include "../FieldSolver.h"

class HoLoInterface : public FieldSolver
{
 public:

  HoLoInterface(PlasmaData* _pdata);

  ~HoLoInterface();

  void init(
	    PlasmaData* PData,
		FieldDataCPU* fields_old,
		NodeFieldData* fields_half,
		FieldDataCPU* fields_next,
		NodeHOMoments* moments);

  void solve(PlasmaData* Pdata,
	     FieldDataCPU* Fields, //output
	     NodeHOMoments* Moments);

  void InitStepSolve(PlasmaData* Pdata,
  	     FieldDataCPU* Fields, //output
  	     NodeHOMoments* Moments);

  realkind calc_residual(PlasmaData* pdata,
		  FieldDataCPU* fields_next,
		  FieldDataCPU* fields_old,
		  NodeHOMoments* moments);

  void update_solution(void);

  void ImportParams();

  void HoToLo( FieldDataCPU* field_in, HOMomentsCPU* moments_in, const Teuchos::RCP<Epetra_Vector>& v);
  void LoToFields();

  // For Stand alone run
  void SolveTimeStepping();

  double HO_resid;

 private:
  
  Epetra_Comm* comm;
  Teuchos::RCP<SimParams> simParams;
  Teuchos::RCP<LowOrderProblem> problem;
  Teuchos::RCP<JFNK> jfnk;

  Teuchos::RCP<DiffMath> OvX;


  Teuchos::RCP<Epetra_Vector> X;
  Teuchos::RCP<Epetra_Vector> Xold;
  Teuchos::RCP<Epetra_Vector> ConsTerms;

  PlasmaData* pdata;
  FieldDataCPU* fields_old;
  NodeFieldData* fields_half;
  FieldDataCPU* fields_next;
  HOMomentsCPU* moments_old;
  HOMomentsCPU* moments_next;

  // local storage for const values
  double me_h ;
  double mi_h ;
  double qe_h ;
  double qi_h ;
  double xi   ;
  double dx   ;
  double dt   ;

  double tol;
  double res;
  bool isConv;

  
};

#endif // TIMESTEPPER_H_
