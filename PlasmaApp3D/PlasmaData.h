/*-------------------------------------------------------------------------*/
/**
	@file	PlasmaData.h
	@author	J. Payne
	@date		1/04/2012
	@brief	Declares the PlasmaData class. Also defines simulation precision,
	several macro calls, constants, and utility functions.


*/
/*--------------------------------------------------------------------------*/

#ifndef PLASMA_DATA_H
#define PLASMA_DATA_H

#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>


#include <mpi.h>

#  define CUDA_SAFE_KERNEL(call) {                                         \
	call;																					\
	cudaDeviceSynchronize();														\
	cudaError err = cudaGetLastError();										\
    if ( cudaSuccess != err) {                                               \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
                exit(EXIT_FAILURE);                                                  \
    } }

#  define CUDA_SAFE_CALL(call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }


__device__
void atomicAddD(double* address,double value);


#define DOUBLE_PRECISION



#define DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION

#define FLOAT_PRECISION double

#define MPI_REALKIND MPI_DOUBLE
#else

#define FLOAT_PRECISION float

	#define MPI_REALKIND MPI_FLOAT
#endif


/*-------------------------------------------------------------------------*/
/**
	@typedef realkind

	Defines floating precision at runtime.
*/
/*--------------------------------------------------------------------------*/
typedef FLOAT_PRECISION realkind;

const realkind qe2me = -1.0;
const realkind qe = -1.0;
const realkind mass_e = 1.0;
const realkind c_light = 1.0;
extern realkind epsilon_naught;// = 8.854187817620e-12;
extern realkind mu_naught;// = 1.2566370614e-6;
extern const realkind pi_const;

#define NSPECIES_MAX 10
#define NSUBCYCLE_BINS 50



/*-------------------------------------------------------------------------*/
/**
	@class PlasmaData
	@author	J. Payne
	@date		12/20/2012
	@brief	Class to control interpolation of the field information
	to the High Order problem

	This class is a pure virtual class that is overloaded with device and
	optimization specific implementations. This parent class merely declares
	the interface for allocating, plotting, interpolating, and setting the
	field data.

	The members declared in this class may be overwritten by children classes
	in order to encompass things such as non-uniform mesh and other optimizations.

	Every child of this class should hold to the following discritization and
	interpolation rules.

	First the Electric and Magnetic fields are chosen to live on a Yee mesh \cite yee1966.

	** Insert Figures here **

	Second, in order for energy conservation the shape functions for the electric
	field interpolation must depend on the interpolation used for the current tally.

	\f{eqnarray}{
	E^{l}\left(\vec{r}_p\right) = \sum_{i,j,k}E^l_{i_l,j_l,k_l}{\cal S}_{E^l}\left(\vec{r}_{i_l,j_l,k_l} - \vec{r}_{p} \right)
	\f}

	With \f$l\f$ indicating the component of \f$E\f$, \f$i_l\f$, \f$j_l\f$, and \f$k_l\f$ are the
	cell index positions where \f$E^l\f$ is stored and are equivilent to \f$f_l = f+1/2\delta_{fl}\f$.
	Where \f${\cal S}_{E^l}\left(\vec{r}_{i_l,j_l,k_l} - \vec{r}_{p} \right)\f$ is
	\f{eqnarray}{
	{\cal S}_{E^l}\left(\vec{r}_{i_l,j_l,k_l} - \vec{r}_{p} \right) = \prod_{m}S^{1+\delta_{lm}}\left(r^m_{i,j,k} - r^m_p \right)
	\f}

	@todo Add more information on Yee mesh, figures, etc. Also add more information on
	shape functions.

*/
/*--------------------------------------------------------------------------*/
class PlasmaData
{
public:

	PlasmaData(){};

	PlasmaData(int argc, char* argv[]);

	void readInput(const std::string& inFile);

	void getKeyword(char* inBuf, std::string& keyword, std::string& rest);

	void update_group(void);

	void set_cuda_device();

	void setup(void)
	{
		epsilon_naught = n0*qe*qe/(qe/qe2me*omega_pe*omega_pe);

		//epsilon_naught = 1.0;

		L_debye = sqrt(epsilon_naught*Te/(qe*qe*n0));

	//	xmin *= L_debye;
	//	ymin *= L_debye;
	//	zmin *= L_debye;

	//	Lx *= L_debye;
	//	Ly *= L_debye;
	//	Lz *= L_debye;



		dxdi = Lx/((realkind)nx);
		dydi = Ly/((realkind)ny);
		dzdi = Lz/((realkind)nz);

		didx = 1.0/dxdi;
		didy = 1.0/dydi;
		didz = 1.0/dzdi;
	}

	realkind dt;
	realkind tmax;
	realkind time_done;

	realkind qspecies[NSPECIES_MAX]; // qs/qe
	realkind mspecies[NSPECIES_MAX]; // ms/me
	realkind wspecies[NSPECIES_MAX]; // statistical weight per species

	realkind epsilon_r, epsilon_a, tol_outer_piccard;

	realkind dxdi,dydi,dzdi;
	realkind didx,didy,didz;

	realkind xmin, ymin, zmin; // In terms of electron Debye Length
	realkind Lx,Ly,Lz; // In terms of electron Debye Length
	realkind L_debye;

	realkind n0; // Plasma Density
	realkind Te; // Electron Temperature
	realkind omega_pe; // Electron plasma frequency

	realkind Bmag_avg;
	realkind Emag_avg;

	int niter_max;
	int nSubcycle_max;
	int nSubcycle_min;
	int nsteps;
	short int nreturn; // How many particles finished before exiting particleObj

	int iParticleListCPU; // Which CPU particle list to use
	int iFieldDataCPU; // Which CPU field data to use

	int nx,ny,nz;
	int nptcls;
	int nspecies;
	int nptcls_total;
	int nptcls_gpu;
	int nptcls_cpu;
	int my_nptcls;
	int gpu_multiplier; // How many more particles to run on GPU than CPU
	int my_species;

	/// Dimension of sorting cluster (GPU)
	int ClusterSortDim;
	/// Dimension of Storage Cluster (GPU)
	int ClusterStorDim;
	/// Number of thread blocks per cluster
	int TBpCluster;


	int nptcls_species[NSPECIES_MAX];

	int prob_type;
	int prob_case;

	int ndimensions; // Number of spatial dimensions
	int nVelocity; // Number of velocity dimensions
	int iEM; // Electrostatic = 0 Electromagnetic = 1


	short int cpu_vec_length;
	int num_cores;
	int mynode;
	int num_nodes; // total number of mpi tasks
	int num_procs; // Number of unique mpi processors
	int my_proc; //
	int myrank_node;
	int device_type;
	int ndevices;
	int ngpus;

	char output_id[16];
	char* SimName; // Name of simulation e.g. IonAcoustic, TwoStream, etc

	bool plot_flag;
	bool lo_all; // do the lo solve on all nodes, no field broadcast.

	int runid;

	int* random_seeds;


	int npiccard_outer;


};


#endif /* PLASMA_DATA_H */
