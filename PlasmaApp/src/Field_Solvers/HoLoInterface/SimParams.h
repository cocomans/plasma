#ifndef SIM_PARAMS_H_
#define SIM_PARAMS_H_

// Trilinos Includes
#include "NOX_Epetra.H"
#include <Epetra_LocalMap.h>
#include "az_aztec.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CommandLineParser.h"
#include "../../PlasmaData.h"
#include "../../Util/webui/LogStream.h"

namespace EM2D3V{
enum {
	n_e,
	n_i,
	px_e,
	py_e,
	pz_e,
	px_i,
	py_i,
	pz_i,
	Ex,
	Ey,
	Ez,
	Bx,
	By,
	Bz,
	Ax,
	Ay,
	Az,
	Chi,
	Phi,
	Sxx_e,
	Sxx_i,
	Sxy_e,
	Sxy_i,
	Sxz_e,
	Sxz_i,
	Syy_e,
	Syy_i,
	Syz_e,
	Syz_i,
	Szz_e,
	Szz_i};
}
namespace EM1D3V{
enum { n_e, n_i,  px_e, py_e, pz_e, px_i,  py_i, pz_i, Ex, Ey, Ez, Bx, By, Bz, Ax, Ay, Az, Sxx_e, Sxx_i, Sxy_e, Sxy_i, Sxz_e, Sxz_i};
}
namespace ES1D1V{
enum { n_e,  p_e, Ex, p_i, n_i,  Sxx_e, Sxx_i};
}


enum { pe_red } ;


class SimParams
{
 public:

  SimParams(Epetra_Comm* Comm,PlasmaData* _pdata);
  ~SimParams();
  
  Teuchos::RCP<Teuchos::ParameterList> GetParamsPhys();
  Teuchos::RCP<Teuchos::ParameterList> GetParamsMap();
  Teuchos::RCP<Teuchos::ParameterList> GetParamsNLS();
  Teuchos::RCP<Teuchos::ParameterList> GetParamsML();
  Teuchos::RCP<Teuchos::ParameterList> GetParamsFlags();

  Teuchos::ParameterList& GetPrintParams();
  Teuchos::ParameterList& GetLSParams();

  void Print();
  PlasmaData* pdata;
 private:

  void ImportParams();

  Epetra_Comm* comm;
  
  Teuchos::RCP<Teuchos::ParameterList> Phys;
  Teuchos::RCP<Teuchos::ParameterList> Flags;
  Teuchos::RCP<Teuchos::ParameterList> Map;
  Teuchos::RCP<Teuchos::ParameterList> NLS;
  Teuchos::RCP<Teuchos::ParameterList> ML;



};

#endif // SIM_PARAMS_H_
