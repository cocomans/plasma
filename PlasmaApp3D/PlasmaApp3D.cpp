#include <mpi.h>







int main(int argc,char* argv[])
{
	int rc,myid,num_nodes;

	rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS)
	{
		printf ("Error starting MPI program. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD,&num_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);


	// Setup plasma data
	PlasmaData pdata(argc,argv);


	// Figure out how many compute devices are on this node

	// Set the number of OMP threads to the number of compute devices

	//





	MPI_Finalize();
}
