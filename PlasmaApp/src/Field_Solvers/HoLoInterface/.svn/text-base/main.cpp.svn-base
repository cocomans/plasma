// Trilinos Includes
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_Time.hpp"

// User Includes
#include "SimParams.h"
#include "HoLoInterface.h"

int main( int argc, char *argv[] )
{
  //////////////////////////////////////
  // Initialization
  //////////////////////////////////////
  #ifdef HAVE_MPI
  MPI_Init( &argc, &argv );
  Epetra_MpiComm comm( MPI_COMM_WORLD );
  std::cout << comm<<std::endl;
  #else
  Epetra_SerialComm comm;
  std::cout << comm <<std::endl;
  #endif
  //////////////////////////////////////

  // Initialize Simulation Parameters
  Teuchos::RCP<SimParams>  simParams =
    Teuchos::rcp(new SimParams(argc, argv, &comm));

  // Create High Order / Low Order interface object
  Teuchos::RCP<HoLoInterface>  interface=
    Teuchos::rcp(new HoLoInterface(simParams, &comm));
  //////////////////////////////////////
  
  //////////////////////////////////////
  // Solve!
  //////////////////////////////////////
  interface->init();
  interface->solve();
  // interface->GetState();

  //////////////////////////////////////
  // Finalize
  //////////////////////////////////////
  #ifdef HAVE_MPI
  MPI_Finalize();
  #endif

  return 0;
}
