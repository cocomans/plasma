#ifndef PARALLEL_INFO_H
#define PARALLEL_INFO_H


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <gnuplot_i.h>
#include <mpi.h>

class RandGen;

class GPUInfo
{
public:
	/// Dimension of sorting cluster (GPU)
	int ClusterSortDim;
	/// Dimension of Storage Cluster (GPU)
	int ClusterStorDim;
	/// Number of thread blocks per cluster
	int TBpCluster;
	/// Number of particles multiplier
	float multiplier;

	/// GPU ID
	int igpu;
};

class MICInfo
{
public:
	/// CPU vector length
	int vec_length;
	/// number of MIC threads to run
	int nthreads;
	/// Dimension of sorting cluster
	int ClusterSortDim;
	/// Dimension of Storage Cluster
	int ClusterStorDim;
	/// Number of particles multiplier
	float multiplier;
};

class CPUInfo
{
public:
	/// Which CPU particle list to use
	char iParticleListCPU;
	/// Which CPU field data to use
	char iFieldDataCPU;
	/// CPU vector length
	int vec_length;
	/// number of CPU threads to run
	int nthreads;
	/// Number of particles multiplier
	float multiplier;

};

class ParallelInfo
{
public:

	/*-------------------------------------------------------------------------*/
	/**
	 * @brief ParallelInfo Constructor, creates sub-classes
	 *
	 */
	/*--------------------------------------------------------------------------*/
	ParallelInfo()
	{
		cpu_info = new CPUInfo();
		gpu_info = new GPUInfo();
		mic_info = new MICInfo();
	}
	/*-------------------------------------------------------------------------*/
	/**
	 * @brief Populate this information, figure out devices,
	 *
	 * @return array containing pointers for information for each device-level thread
	 * ie one for CPU, one for each GPU, and one for each MIC on node
	 *
	 */
	/*--------------------------------------------------------------------------*/
	ParallelInfo** init_info(int argc,char* argv[]);

	void update_group(void);

	void count_devices();

	void set_affinity(void);



	/// Global MPI rank
	int rank_g;
	/// Node MPI rank
	int rank_n;
	/// MPI node
	int inode;

	/// Number of unique mpi nodes
	int nNodes;
	/// Number of MPI tasks global
	int nTasks_g;
	/// Number of MPI tasks on node
	int nTasks_n;

	/// number of GPUs
	int nGPU;
	/// number of MIC devices
	int nMIC;
	/// number of CPU threads to run
	int nCPU;
	/// Number of physical processors on node
	int nProcs;

	/// Total number of computational devices (1 + nGPU + nMIC)
	int nDevices;
	int nDevices2;
	/// Total number of OpenMP threads (nCPU + nGPU + nMIC)
	int nThreads;
	/// Total number of distinct processing units (nProcs + nGPU + nMIC)
	int nUnits;

	/// Type of device to use 0=CPU, 1=GPU, 2=MIC
	int device_type;
	/// Which GPU / CPU / MIC am I associated with
	int device_id;

	/// Which species to use
	int ispecies;
	int nspecies;


	int random_seed;
	/// Random number generator for this thread
	RandGen* generator;

	CPUInfo* cpu_info;
	GPUInfo* gpu_info;
	MICInfo* mic_info;
};



#endif /* PARALLEL_INFO_H */
