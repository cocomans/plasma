/*-------------------------------------------------------------------------*/
/**
	@file	PlasmaData.h
	@author	J. Payne
	@date		1/04/2013
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
#include <string>
#include <string.h>
#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include "cuda_defines.h"
#include <gnuplot_i.h>


#define SIGN_MASK 0x7FFFFFFF

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

#ifndef NO_CUDA
#  define CUDA_SAFE_CALL(call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }
#else
#  define CUDA_SAFE_CALL(call) {}
#endif

__device__
void atomicAddD(double* address,double value);

#ifndef NO_DOUBLE_PREC
#define DOUBLE_PRECISION
#endif


//#define DOUBLE_PRECISION

#ifdef DOUBLE_PRECISION

#define FLOAT_PRECISION double

#define MPI_REALKIND MPI_DOUBLE
#else

#define FLOAT_PRECISION float

	#define MPI_REALKIND MPI_FLOAT
#endif

#ifdef DOUBLE_PRECISION
#define HALF_r 0.5
#define QUARTER_r 0.25
#define ZERO_r 0.0
#define CELL_OFFSET_r 0.5+1.0e-12
#define ONE_r 1.0
#define ONE_5_r 1.5
#define ZERO_75_r 1.75
#else
#define HALF_r 0.5f
#define QUARTER_r 0.25f
#define ZERO_r 0.0f
#define CELL_OFFSET_r 0.5f+1.0e-5f
#define ONE_r 1.0f
#define ONE_5_r 1.5f
#define ZERO_75_r 0.75f
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

class ParallelInfo;
class RunStatus;
class RunData;
class RandomGenerator;
class OutputCtrl;
class Normalization;
class LogStream;
class sseServer;

extern LogStream  LO_log;
extern LogStream  HO_log;
extern LogStream  Outer_log;
extern LogStream  debug_log;



/// Initial Condition Parameters
#define qe_p -1.60217657e-19 // Coulomb per electron
#define qi_p 1.60217657e-19 // Coulomb per proton
const realkind me_p = 9.10938215e-31; // Mass of electron in kg
const realkind mi_p1 = 1.672621777e-27; // Mass of proton in kg
const realkind epsilon_0_p = 8.854187817620e-12; // s^2*C^2/(m^3*kg)
const realkind mu_0_p = 1.2566370614e-6; // m*kg/C^2

const int LINESIZE = 512;


class InputParser
{
public:
	__attribute__((noinline))
	InputParser(int argc, char** argv);


	int argc;
	char** argv;
	std::string* inFile;
	std::ifstream* inStr;
	// Flag that indicates input file exists
	bool iFile;
	bool fileOpen;

	template<typename T> inline
	void GetParam(const char* fileKey, const char* comKey,T& output, T defVal)
	{
		using namespace std;

		output = defVal;

		if(iFile == 1)
		{
			// Open the file if it got closed
			if(fileOpen != 1)
			{
				fileOpen = 1;
				inStr = new ifstream(inFile->c_str());
				if (!inStr) {
					cout << "Could not open input file " << inFile << endl;
					exit(-1);
				}
			}

			// Set the file buffer to the beginning of the file
			inStr->clear();
			inStr->seekg(0, ios::beg);

			char inBuf[LINESIZE];
			string keyword;
			string rest;

			while (inStr->getline(inBuf, LINESIZE)) {
				if (inBuf[0] != '#' && inStr->gcount() > 1) {

					getKeyword(inBuf, keyword, rest);
					istringstream line(rest.c_str());
	//				cout << "KEYWORD: " << keyword << "   " << rest << endl;

					if (keyword == fileKey)
						line >> output;
				}
			}
		}

		// Command Line input overwrites file input
		for(int i=0;i<argc;i++)
		{
			std::string arg(argv[i]);

			if(arg == comKey)
				istringstream(argv[i+1]) >> output;
		}

		cout << string(fileKey) << string(" = ") << (output) << endl;


	}

	__attribute__((noinline))
	void getKeyword(char* inBuf, std::string& keyword, std::string& rest);

};

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


	PlasmaData(int argc,char* argv[]);

	PlasmaData(int argc,char* argv[],Normalization* normalizer);


	void readInput(const std::string& inFile);

	void getKeyword(char* inBuf, std::string& keyword, std::string& rest);

	void set_cuda_device();

	void setup(void)
	{
		epsilon_naught = n0*qe*qe/(qe/qe2me*omega_pe*omega_pe);

		//epsilon_naught = 1.0;

		L_debye = sqrt(epsilon_naught*Te/(qe*qe*n0));



		dxdi = Lx/((realkind)nx);
		dydi = Ly/((realkind)ny);
		dzdi = Lz/((realkind)nz);

		didx = 1.0/dxdi;
		didy = 1.0/dydi;
		didz = 1.0/dzdi;
	}

	void parse_normalization(int argc, char* argv[]);

	realkind dt;
	realkind tmax;
	realkind time_done;

	realkind qspecies[NSPECIES_MAX]; // qs/qe
	realkind mspecies[NSPECIES_MAX]; // ms/me
	realkind wspecies[NSPECIES_MAX]; // statistical weight per species

	realkind epsilon_r, epsilon_a, tol_outer_picard;

	realkind dxdi,dydi,dzdi;
	realkind didx,didy,didz;

	realkind xmin, ymin, zmin; // In terms of electron Debye Length
	realkind L_debye;

	realkind Te; // Electron Temperature
	realkind omega_pe; // Electron plasma frequency

	//realkind Bmag_avg;
	//realkind Emag_avg;

	int niter_max;
	int nSubcycle_max;
	int nSubcycle_min;
	int nOuterPicard_max;
	int nsteps;
	short int nreturn; // How many particles finished before exiting particleObj


	// naught quantities (user inputed in si units)
	/// reference mass in kg
	realkind m0;
	/// Reference charge in C
	realkind q0;
	/// Reference Length in m
	realkind L0;
	/// Reference number density in atom/m^3
	realkind n0;
	/// Reference Temperature in Joules
	realkind T0;




	/// System Lengths in m
	realkind Lx_p;
	realkind Ly_p;
	realkind Lz_p;
	/// Electron temperature in physical units (same units as T0)
	realkind Te_p;
	/// ion temperature in physical units (same units as T0)
	realkind Ti_p;
	/// Electron density in physical units (same units as n0)
	realkind ne_p;
	/// Ion density in physical units (same units as n0)
	realkind ni_p;

	realkind mi_p;

	/// q_alpha / q0
	realkind q_h[NSPECIES_MAX];
	/// m_alpha / m0
	realkind m_h[NSPECIES_MAX];
	/// n_bar / n0
	realkind n_bh[NSPECIES_MAX];
	/// T_alpha / T0
	realkind T_h[NSPECIES_MAX];
	/// System lengths Lx_p / L0, Ly_p / L0 ...
	realkind Lx,Ly,Lz;
	/// Constant normalized parameter
	realkind xi;
	realkind bobzeta;

	realkind spacer[10];

	realkind eps0_h;
	realkind mu0_h;






	int nx,ny,nz;
	/// Number of particles to run on node
	int nptcls;
	int nspecies;

	int nptcls_device[3];
	float device_multiplier[3];
	int nptcls_species[NSPECIES_MAX];

	int prob_type;
	int prob_case;

	int ndimensions; // Number of spatial dimensions
	int nVelocity; // Number of velocity dimensions
	int iEM; // Electrostatic = 0 Electromagnetic = 1


	ParallelInfo* node_info;
	ParallelInfo** thread_info;
	RunData* rdata;
	RunStatus* rstatus;
	RandomGenerator** RandGens;
	OutputCtrl* output_ctrl;

	InputParser* parser;

	int argc_s;
	char** argv_s;

	std::string* pythonScript;




};


#endif /* PLASMA_DATA_H */
