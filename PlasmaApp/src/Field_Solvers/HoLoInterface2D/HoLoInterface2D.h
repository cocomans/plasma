#ifndef HOLO_INTERFACE_2D_H_
#define HOLO_INTERFACE_2D_H_

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_Time.hpp"

// User Includes
#include "../HoLoInterface/SimParams.h"
#include <NOX_Epetra.H>

#include "JFNK2D.h"
#include "LowOrderProblem2D.h"

#include "EpVecWrapper.h"

//=========================
//Include HO header-files
//========================
#include "../../HOMomentsCPU.h"
#include "../../FieldDataCPU.h"
#include "../../NodeHOMoments.h"
#include "../../PlasmaData.h"
#include "../FieldSolver.h"

class HoLoInterface2D : public FieldSolver
{
 public:

  HoLoInterface2D(PlasmaData* _pdata);

  ~HoLoInterface2D();

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

  void HoToLo( FieldDataCPU* field_in, HOMomentsCPU* moments_in,
		  Teuchos::RCP<EpVecWrapper>& v);
  void LoToFields();

  void LoToMoments();


  // For Stand alone run
  void SolveTimeStepping();

  double HO_resid;


  
  Epetra_Comm* comm;
  Teuchos::RCP<SimParams> simParams;
  Teuchos::RCP<LowOrderProblem2D> problem;
  Teuchos::RCP<JFNK2D> jfnk;


  Teuchos::RCP<EpVecWrapper> X;
  Teuchos::RCP<EpVecWrapper> Xold;
  Teuchos::RCP<EpVecWrapper> ConsTerms;
  Teuchos::RCP<EpVecWrapper> ResVec;


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

#endif // HOLO_INTERFACE_2D_H_
