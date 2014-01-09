/*-------------------------------------------------------------------------*/
/**
	@file	ParallelInfo.cpp
	@author	J. Payne
	@date		04/30/2013
	@brief	Defines the ParallelInfo methods


*/
/*--------------------------------------------------------------------------*/
#include "ParallelInfo.h"
#include "PlasmaData.h"
#include <list>
#include <string>
#include <string.h>
#include <cctype>
#include <omp.h>
#include <mpi.h>
#include <errno.h>
#define handle_error_en(en, msg) \
	   do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)

MPI_Group* node_groups;
MPI_Comm* node_comms;



	/*-------------------------------------------------------------------------*/
	/**
	 * @brief Populate this information, figure out devices,
	 *
	 * @return array containing pointers for information for each device-level thread
	 * ie one for CPU, one for each GPU, and one for each MIC on node
	 *
	 */
	/*--------------------------------------------------------------------------*/
ParallelInfo** ParallelInfo::init_info(int argc,char* argv[])
{

	InputParser* parser = new InputParser(argc,argv);

	this->cpu_info->vec_length = 32;
	this->cpu_info->iParticleListCPU = 0;
	this->cpu_info->iFieldDataCPU = 1;
	this->gpu_info->ClusterSortDim = 2;
	this->gpu_info->ClusterStorDim = 32;
	this->gpu_info->TBpCluster = 1;

	this->nCPU = omp_get_num_procs();
	this->nProcs = omp_get_num_procs();
	this->nGPU = 0;
	this->nMIC = 0;

	device_type = 0;
	device_id = 0;
	ispecies = 0;

	nspecies = 1;

	random_seed = rand()%1023789456;
	srand(random_seed);



	MPI_Comm_size(MPI_COMM_WORLD,&(nTasks_g));

	MPI_Comm_rank(MPI_COMM_WORLD,&(rank_g));

	update_group();

	count_devices();

	MPI_Barrier(MPI_COMM_WORLD);
//	system("nvidia-smi");
	int ndevices, ntesla;
	int itesla[20];

	#ifndef NO_CUDA

			// First get the total number of cuda capable devices

			cudaGetDeviceCount(&ndevices);




			// now loop through the devices, we are looking for tesla cards,
			// The main thing that sets all Tesla cards apart is 2 async engines
			ntesla = 0;
			for(int i=0;i<ndevices;i++)
			{
				cudaDeviceProp prop;
				CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop,i));
				printf("Device[%i] name = %s\n",i,prop.name);
				std::string name(prop.name);
				if(name.find("Tesla") != std::string::npos)
				{
					itesla[ntesla] = i;
					ntesla += 1;
				}
			}

			nGPU=ntesla;

//			nGPU = 1;

//			nGPU=2;
	//	device_type = 0;

	#endif

	parser->GetParam("CPU_VECTOR_LENGTH","--vec_length",cpu_info->vec_length,32);


	parser->GetParam("NUMBER_CPU_CORES","--num-cores",nCPU,32);

	parser->GetParam("GPU_NPTCL_MULTIPLIER","--gpu-mult",gpu_info->multiplier,1.0f);

	parser->GetParam("NUMBER_SPECIES","-ns",nspecies,1);




	set_affinity();


	// Command line arguments overwrite input file
	for(int i=0;i<argc;i++)
	{
		std::string arg(argv[i]);


		if(arg == ("--vec_length"))
			cpu_info->vec_length = atoi(argv[i+1]);

		else if(arg == "--num_cores")
			nCPU = atoi(argv[i+1]);

		else if(arg == "--gpu-mult")
			gpu_info->multiplier = atof(argv[i+1]);

		else if(arg == "--plist-cpu")
			cpu_info->iParticleListCPU = atoi(argv[i+1]);
		else if(arg == "--fdata-cpu")
			cpu_info->iFieldDataCPU = atoi(argv[i+1]);
		else if(arg == ("-ns"))
			nspecies = atoi(argv[i+1]);

	}


	// Calculate some generated parameters
	nThreads = nCPU+nGPU+nMIC;
	nUnits = nProcs+nGPU+nMIC;
	nDevices = 1+nGPU+nMIC;
	nDevices2 = nspecies+nGPU+nMIC;

	ParallelInfo** thread_info = (ParallelInfo**)malloc(nDevices2*sizeof(ParallelInfo*));
	// Initialize everything for each device-thread
	for(int i=0;i<nDevices2;i++)
	{
		thread_info[i] = new ParallelInfo();
		CPUInfo* temp_c = thread_info[i] -> cpu_info;
		GPUInfo* temp_g = thread_info[i] -> gpu_info;
		MICInfo* temp_m = thread_info[i] -> mic_info;
		thread_info[i][0] = *this;

		thread_info[i] -> cpu_info = temp_c;
		thread_info[i] -> gpu_info = temp_g;
		thread_info[i] -> mic_info = temp_m;

		if(i < nspecies)
		{
			// CPU thread info
			thread_info[i] -> device_type = 0;
			thread_info[i] -> device_id = i;
			thread_info[i] -> ispecies = i;
		}
		else if(i < (nGPU+nspecies))
		{
			// GPU thread info
			thread_info[i] -> device_type = 1;
			thread_info[i] -> device_id = i-nspecies;
			thread_info[i] -> gpu_info -> igpu = itesla[i-nspecies];
			thread_info[i] -> ispecies = 0;
			printf("thread %i has gpu %i with device_id %i\n",i,thread_info[i]->gpu_info->igpu,thread_info[i]->device_id);
		}
		else if(i < (nGPU+nspecies+nMIC))
		{
			// GPU thread info
			thread_info[i] -> device_type = 2;
			thread_info[i] -> device_id = i-(nGPU+nspecies);
		}
	}




	return thread_info;

}

void ParallelInfo::set_affinity(void)
{
//    int s, j;
//    cpu_set_t cpuset;
//    pthread_t thread;
//
//    thread = pthread_self();
//
//    /* Set affinity mask to include CPUs 0 to nCPU */
//
//    CPU_ZERO(&cpuset);
//    for (j = 0; j < nCPU; j++)
//        CPU_SET(j, &cpuset);
//
//    s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
//    if (s != 0)
//        handle_error_en(s, "pthread_setaffinity_np");
//
//    /* Check the actual affinity mask assigned to the thread */
//
//    s = pthread_getaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
//    if (s != 0)
//        handle_error_en(s, "pthread_getaffinity_np");

//    printf("Set returned by pthread_getaffinity_np() contained:\n");
//    for (j = 0; j < CPU_SETSIZE; j++)
//        if (CPU_ISSET(j, &cpuset))
//            printf("    CPU %d\n", j);
}

using namespace std;
void ParallelInfo::count_devices(void)
{
	this->nGPU = 0;
	this->nMIC = 0;
}

// Function to count the number of unique nodes, and assign a group integer to each node
void ParallelInfo::update_group(void)
{
	// There has got to be a better way of doing this, but well this should work for now...

	char processor_names[nTasks_g][MPI_MAX_PROCESSOR_NAME];
	int processor_name_lengths[nTasks_g];

	char my_processor_name[MPI_MAX_PROCESSOR_NAME];
	int my_processor_name_length;

	// Every node gets its processor name
	printf("Getting Names for %i nodes\n",nTasks_g);
	MPI_Get_processor_name(my_processor_name,&my_processor_name_length);
	printf("My name is %s\n",my_processor_name);
	printf("Reducing Name Lengths\n");
	if(rank_g != 0)
	{
		// send your processor name length to the root node
		MPI_Send(&my_processor_name_length,1,
				MPI_INT,0,rank_g,MPI_COMM_WORLD);
	}
	else
	{
		// recieve the processor name lengths
		for(int i=1;i<nTasks_g;i++)
		{
			MPI_Recv(&processor_name_lengths[i],1,
					MPI_INT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	printf("Reducing Names\n");
	if(rank_g != 0)
	{
		// send your processor name to the root node
		MPI_Send(my_processor_name,MPI_MAX_PROCESSOR_NAME,
				MPI_CHAR,0,rank_g+nTasks_g,MPI_COMM_WORLD);
	}
	else
	{
		strcpy(processor_names[0], my_processor_name);
		processor_name_lengths[0] = my_processor_name_length;
		// recieve the processor names
		for(int i=1;i<nTasks_g;i++)
		{
			MPI_Recv(processor_names[i],MPI_MAX_PROCESSOR_NAME,
					MPI_CHAR,i,i+nTasks_g,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

//			printf("Reciveing name[%i]: %s\n",mynode,processor_names[i]);

		}
	}


	MPI_Barrier(MPI_COMM_WORLD);

	printf("Broadcasting Names\n");

	// Broadcast all of the names
	for(int i=0;i<nTasks_g;i++)
	{
		MPI_Bcast(processor_name_lengths+i,1,MPI_INT,0,MPI_COMM_WORLD);

	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(int i=0;i<nTasks_g;i++)
	{

		MPI_Bcast(processor_names[i],MPI_MAX_PROCESSOR_NAME,
				MPI_CHAR,0,MPI_COMM_WORLD);

	}


	MPI_Barrier(MPI_COMM_WORLD);

	const char* unique_processor_names0[nTasks_g];
	int unique_processor_name_lengths0[nTasks_g];
	int num_tasks_per_processor[nTasks_g];
	int num_unique_procs;
	int k = 0;
	if(rank_g == 0){
		printf("Sorting Names\n");
	list<string> unique_procs;
	list<string> unique_procs2;
	for(int i=0;i<nTasks_g;i++)
	{
		string temp;
		temp = string(processor_names[i]);
		unique_procs.push_back(temp);
	}

	// sort the list
	unique_procs.sort();

	unique_procs2 = unique_procs;

	// remove duplicate values
	unique_procs2.unique();



	for(list<string>::iterator it=unique_procs2.begin();it!=unique_procs2.end();++it)
	{
		unique_processor_name_lengths0[k] = it->size();
		unique_processor_names0[k] = it->c_str();

		//memcpy(unique_processor_names0[k],it->c_str(),it->size()*sizeof(char));

		k++;
	}
	}


	MPI_Barrier(MPI_COMM_WORLD);
	nNodes = k;

	// Broadcast the number of names to everyone

	MPI_Bcast(&nNodes,1,MPI_INT,0,MPI_COMM_WORLD);


	char unique_processor_names[nNodes][MPI_MAX_PROCESSOR_NAME];
	int unique_processor_name_lengths[nNodes];


	for(int i=0;i<nNodes;i++)
	{
		if(rank_g == 0)
//		printf("Broadcasting length of names = %i\n",unique_processor_name_lengths0[i]);
		unique_processor_name_lengths[i] = unique_processor_name_lengths0[i];
		MPI_Bcast(unique_processor_name_lengths+i,1,MPI_INT,0,MPI_COMM_WORLD);

		if(rank_g == 0)
		{
			strcpy(unique_processor_names[i],unique_processor_names0[i]);
//			printf("Broadcasting name: %s\n",unique_processor_names[i]);
		}
		MPI_Bcast(unique_processor_names[i],MPI_MAX_PROCESSOR_NAME,
				MPI_CHAR,0,MPI_COMM_WORLD);

	}

	MPI_Barrier(MPI_COMM_WORLD);

	int* group_members[nNodes];
	int num_group_members[nNodes];

	if(rank_g < nNodes)
	{
//		printf("Grouping names %i\n",mynode);
		// since this is parallelizable we might as well

		// count the number of tasks with a given processor name
		num_group_members[rank_g] = 0;
		string str1,str2;
		str1 = string(unique_processor_names[rank_g]);
		for(int i=0;i<nTasks_g;i++)
		{
			str2 = string(processor_names[i]);

			cout << str1 << " " << str2 << endl;

			if(str1 == str2)
				num_group_members[rank_g]++;

		}

		// allocate the number of groups members
		group_members[rank_g] = (int*)malloc(num_group_members[rank_g]*sizeof(int));

		printf("group %i has %i members\n",rank_g,num_group_members[rank_g]);
		// Setup the list of group members
		str1 = unique_processor_names[rank_g];
		k = 0;
		for(int i=0;i<nTasks_g;i++)
		{
			str2 = processor_names[i];

			if(str1 == str2)
			{
				group_members[rank_g][k] = i;
				k++;
			}

		}




	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Broadcast all of the group members
	printf("Broadcasting group members\n");
	for(int i=0;i<nNodes;i++)
	{
		MPI_Bcast(&num_group_members[i],1,MPI_INT,i,MPI_COMM_WORLD);

		if(rank_g != i)
			group_members[i] = (int*)malloc(num_group_members[i]*sizeof(int));
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(int i=0;i<nNodes;i++)
	{
		MPI_Bcast(group_members[i],num_group_members[i],
				MPI_INT,i,MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	node_groups = (MPI_Group*)malloc(nNodes*sizeof(MPI_Group));
	node_comms = (MPI_Comm*)malloc(nNodes*sizeof(MPI_Comm));

	// Create the MPI Groups
	MPI_Group orig_group;

	// Get the original group
	MPI_Comm_group(MPI_COMM_WORLD,&orig_group);
	MPI_Barrier(MPI_COMM_WORLD);
	// Divide tasks into processor based groups

	for(int i=0;i<nNodes;i++)
	{
		bool ibelong = false;

		for(int j=0;j<num_group_members[i];j++)
		{
			if(rank_g == group_members[i][j])
			{
				ibelong = true;
				inode = i;
			}

		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Group_incl(orig_group,num_group_members[i],group_members[i],node_groups+i);

		if(ibelong)
		{

			if(rank_g == group_members[i][0])
			{
				printf("Node %i Creating group %i\n",rank_g,i);
			}
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);
	// Create MPI Communicators
	for(int i=0;i<nNodes;i++)
	{
		if(rank_g == group_members[i][0])
		printf("Creating comm for group %i\n",i);
		MPI_Comm_create(MPI_COMM_WORLD,node_groups[i],node_comms+i);
	}

	// Get the on node rank
	MPI_Comm_rank(node_comms[inode],&rank_n);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("rank = %i, node = %i, node rank = %i\n",rank_g,inode,rank_n);

//	getchar();
	MPI_Barrier(MPI_COMM_WORLD);

}

