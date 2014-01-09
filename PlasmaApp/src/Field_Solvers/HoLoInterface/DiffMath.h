#ifndef DISCRETEMATH_H_
#define DISCRETEMATH_H_
// Trilinos Includes
#include "NOX_Epetra.H"
#include "MapManager.h"

class DiffMath
{
 public:

  DiffMath( const Teuchos::RCP<MapManager>& map_manager,
	    const Teuchos::RCP<SimParams> &Params,
	    Epetra_Comm* Comm );
  ~DiffMath();

  void SetVector(const Epetra_Vector& x);

  double& Get   (int t, int i);
  double GetSM(int t, int i);
  double GetC2F(int t, int i);
  double GetF2C(int t, int i);

  double Mean(int t);
   
  double Dx   (int t, int i, double dx);
  double Dxx   (int t, int i, double dx);
  double DxC2F(int t, int i, double dx);
  double DxF2C(int t, int i, double dx);

  double Dx   (double W, double C, double E, double dx);
  double DxC2F(double W, double C, double E, double dx);
  double DxF2C(double W, double C, double E, double dx);

 private:

  Epetra_Comm*                comm;
  Teuchos::RCP<SimParams>     simParams;
  Teuchos::RCP<MapManager>    mapManager;
  Teuchos::RCP<Epetra_Vector> X;

  Teuchos::RCP<Epetra_Vector> avg;

};
#endif // DISCRETEMATH_H_
