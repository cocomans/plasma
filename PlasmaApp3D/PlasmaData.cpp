#include "PlasmaData.h"
#include <mpi.h>
#include <omp.h>
#include <list>
#include <string>
#include <string.h>
#include <cctype>
#include <tr1/random>
#include <random>


MPI_Group* node_groups;
MPI_Comm* node_comms;

const int LINESIZE = 512;

realkind epsilon_naught;
realkind mu_naught;

const realkind pi_const = 3.14159265359;

std::tr1::uniform_int<int> rand_unif;

std::tr1::ranlux64_base_01 gen1;
std::tr1::mt19937 gen2;
std::tr1::linear_congruential<unsigned long, 33614, 0, 2147483647> seeder;
std::tr1::ranlux64_base_01 generator;


std::tr1::xor_combine<std::tr1::ranlux64_base_01,7,std::tr1::mt19937,13> xorgen;

bool irand = 0;

double drandom()
{
	if(!irand)
	{
		// initialize the seeder
		seeder.seed((rand()%2386157)+(rand()%98463));



		// combined engine
		xorgen = *(new std::tr1::xor_combine<std::tr1::ranlux64_base_01,7,std::tr1::mt19937,13>(gen1, gen2));

		// initialize the combined generator
		xorgen.seed(seeder);



		rand_unif = *(new std::tr1::uniform_int<int>(0,1000000-1));

		irand = 1;
	}


	return (rand_unif(xorgen)+1.0e-6)*1.0e-6;

}

int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}

PlasmaData::PlasmaData(int argc, char* argv[])
{

//	system("numactl --nodebind=0,1 --physcpubind=+0-7");
	system("numactl --show");

	setvbuf( stdout, NULL, _IONBF, 0 );

	// Default Values
	nx = 32;
	ny = 32;
	nz = 32;

	nspecies = 1;
	nptcls = 6400;

	Lx = 4.0*pi_const;
	Ly = 1.0;
	Lz = 1.0;

	xmin = 0;
	ymin = 0;
	zmin = 0;

	n0 = 1;

	omega_pe = 1.0;
	Te = 1.0;

	epsilon_r = 1.0e-10;
	epsilon_a = 1.0e-4;
	niter_max = 20;
	nSubcycle_max = 1000;
	nSubcycle_min = 4;

	cpu_vec_length = 32;

	dt = 0.5;
	nsteps = 20;

	ndimensions = 1;
	nVelocity = 1;
	iEM = 0;

	Bmag_avg = 0;


	num_cores = omp_get_num_procs();

	plot_flag = 0;
	lo_all = 0;
	runid = 0;

	gpu_multiplier = 1;

	nreturn = 1;
	iParticleListCPU = 0;
	iFieldDataCPU = 0;

	SimName = "unknown";



	// Check for an input file
	for(int i=0;i<argc;i++)
	{
		std::string arg(argv[i]);
		if(arg == ("-f"))
		{
			// Need to read the input file first
			std::string inFile(argv[i+1]);
			readInput(inFile);
			break;

		}
	}


	// Command line arguments overwrite input file
	for(int i=0;i<argc;i++)
	{
		std::string arg(argv[i]);

		if(arg == ("-nx"))
			nx = atoi(argv[i+1]);
		else if(arg == ("-ny"))
			ny = atoi(argv[i+1]);
		else if(arg == ("-nz"))
			nz = atoi(argv[i+1]);

		else if(arg == ("-Lx"))
			Lx = atof(argv[i+1]);
		else if(arg == ("-Ly"))
			Ly = atof(argv[i+1]);
		else if(arg == ("-Lz"))
			Lz = atof(argv[i+1]);

		else if(arg == ("-x0"))
			xmin = atof(argv[i+1]);
		else if(arg == ("-y0"))
			ymin = atof(argv[i+1]);
		else if(arg == ("-z0"))
			zmin = atof(argv[i+1]);

		else if(arg == ("--vec_length"))
			cpu_vec_length = atoi(argv[i+1]);
		else if(arg == ("-dt"))
			dt = atof(argv[i+1]);
		else if(arg == ("-s"))
			nsteps = atoi(argv[i+1]);

		else if(arg == ("-np"))
			nptcls = atoi(argv[i+1]);
		else if(arg == ("-ns"))
			nspecies = atoi(argv[i+1]);
		else if(arg == "--min-subcycles")
			nSubcycle_min = atoi(argv[i+1]);

		else if(arg == "--epsilon")
			epsilon_r = atof(argv[i+1]);

		else if(arg == "-n0")
			n0 = atof(argv[i+1]);
		else if(arg == "-Te")
			Te = atof(argv[i+1]);
		else if(arg == "-We")
			omega_pe = atof(argv[i+1]);

		else if(arg == "--num_cores")
			num_cores = atoi(argv[i+1]);

		else if(arg == "--gpu-mult")
			gpu_multiplier = atoi(argv[i+1]);

		else if(arg == "-g")
			plot_flag = 1;
		else if(arg == "--lo-all")
			lo_all = 1;
		else if(arg == "--runid")
			runid = atoi(argv[i+1]);

		else if(arg == "--nSpatial")
			ndimensions = atoi(argv[i+1]);
		else if(arg == "--nVel")
			nVelocity = atoi(argv[i+1]);
		else if(arg == "--plist-cpu")
			iParticleListCPU = atoi(argv[i+1]);

	}

	MPI_Comm_size(MPI_COMM_WORLD,&num_nodes);

	MPI_Comm_rank(MPI_COMM_WORLD,&mynode);

	update_group();

	this->set_cuda_device();


	// Check to make sure that nx, ny, and nz are powers of 2
	if(!isPowerOfTwo(nx))
	{	if(mynode == 0)
		printf("Error!: nx = %i is not a power of 2\n",nx);
		exit(0);
	}
	else if(!isPowerOfTwo(ny))
	{	if(mynode == 0)
		printf("Error!: ny = %i is not a power of 2\n",ny);
		exit(0);
	}
	else if(!isPowerOfTwo(nz))
	{
		if(mynode == 0)
		printf("Error!: nz = %i is not a power of 2\n",nz);
		exit(0);
	}


	this->setup();



	my_nptcls = nptcls*((gpu_multiplier-1)*(1-device_type) +1);

	nptcls_cpu = nptcls;

	nptcls = my_nptcls;
	float mion = 2000.0;

	int device_type_temp = device_type + 1;
	int device_type_temp2 = device_type == 1;
	MPI_Allreduce(&device_type_temp,&ndevices,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&device_type_temp2,&ngpus,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	if(nspecies == 1)// Single Species case
		my_species = 0;
	else // Multi Species Case
	{
		my_species = mynode%2;
		// Electrons on GPU, ions on CPU.
		if(device_type == 1)
			my_species = 1;
		else
		{
//			my_species = (mynode%(nspecies-1)) + 1;
			my_species = myrank_node%nspecies;
			//my_species = 0;

		}

	}



//	mion = 100;

	int nptcls_species_temp[NSPECIES_MAX];
	for(int i=0;i<NSPECIES_MAX;i++)
	{
		nptcls_species_temp[i] = 0;
	}

	nptcls_species_temp[my_species] = my_nptcls;

	// Reduce the number of particles for each species
	MPI_Allreduce(nptcls_species_temp,nptcls_species,
			NSPECIES_MAX,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	//printf("num nodes = %i\n",num_nodes);

	for(int i=0;i<NSPECIES_MAX;i++)
	{
		qspecies[i] = -1;
		mspecies[i] = i*mion;
	}

	qspecies[0] = 1;
	mspecies[0] = 1;

	for(int i=0;i<NSPECIES_MAX;i++)
	{
		wspecies[i] = n0*Lx*Ly*Lz/((double)nptcls_species[i]);

		printf("wspecies[%i] = %e\n",i,wspecies[i]);
	}

	unsigned long long int FNV_prime = 1099511628211ULL;
	unsigned long long int FNV_offset_basis = 14695981039346656037ULL;
	char data_octet[19];
	unsigned long long int hash = FNV_offset_basis;

	int nptcls_t=0;
	for(int i=0;i<nspecies;i++)
		nptcls_t += nptcls_species[i];

	data_octet[0] = nx;
	data_octet[1] = nptcls_t*1.0e-5;
	data_octet[2] = ndimensions+3*nVelocity;
	data_octet[3] = Lx;
	data_octet[4] = nspecies;
	data_octet[5] = num_nodes;
	data_octet[6] = num_cores;
	data_octet[7] = cpu_vec_length;
	data_octet[8] = sizeof(realkind);
	data_octet[9] = dt*10;
	data_octet[10] = nsteps;
	data_octet[11] = lo_all+1;
	data_octet[12] = gpu_multiplier;
	data_octet[13] = iEM+1;
	data_octet[14] = ny;
	data_octet[15] = nz;
	data_octet[16] = Ly;
	data_octet[17] = Lz;
	data_octet[18] = ndevices;

	for(int i=0;i<19;i++)
	{
		hash *= FNV_prime;
		hash = hash ^ data_octet[i];
	}

	sprintf(output_id,"%X",hash);

	timespec res;
	clock_getres(CLOCK_REALTIME,&res);
	tmax = nsteps * dt;

	nreturn = std::min(nreturn,cpu_vec_length);

	if(mynode == 0){

	// Print out the simulation parameters
	printf("############################################\n");
	printf("Simulation will be run with the following parameters:\n");
	// Number of cells
	printf("\n/********* Spatial Parameters *********/\n");
	printf("Grid Dimensions (nx,ny,nz): %i, %i, %i\n",nx,ny,nz);
	printf("System Lengths (Lx,Ly,Lz): %e, %e, %e\n",Lx,Ly,Lz);
	printf("System Origins (x0,y0,z0): %e, %e, %e\n",xmin,ymin,zmin);

	printf("\n/********* Particle Parameters *********/\n");
	printf("Number of particles per species: %i\n",nptcls);
	printf("Number of species: %i\n",nspecies);

	printf("\n/********* Time Step Parameters *********/\n");
	printf("Total Simulation Time: %e\n",tmax);
	printf("Time Step Size: %e\n",dt);

	printf("\n/********* Other Parameters *********/\n");
	printf("Relative Tolerance: %e\n",epsilon_r);
	printf("Absolute Tolerance: %e\n",epsilon_a);
	printf("CPU Vector Length: %i\n",cpu_vec_length);
	printf("Number of CPU Cores to Use: %i\n",num_cores);
	printf("Output ID: %s\n",output_id);
	printf("Clock resolution: %e (ns)\n",1.0e9*res.tv_sec+res.tv_nsec);
	printf("Sizeof or realkind is %i\n",sizeof(realkind));
	printf("############################################\n");

	}


}

void PlasmaData::set_cuda_device()
{
	device_type = 0;
#ifndef NO_CUDA
	if(myrank_node > 0)
	{
		// First get the total number of cuda capable devices
		int ndevices, ntesla;
		CUDA_SAFE_CALL(cudaGetDeviceCount(&ndevices));

		int itesla[ndevices];

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

		if((myrank_node-1) < ntesla)
		{
			device_type = 1;
			printf("node %i setting device %i\n",myrank_node,itesla[myrank_node-1]);
			CUDA_SAFE_CALL(cudaSetDevice(itesla[myrank_node-1]));
		}
	}

	device_type = 0;
#endif
}

using namespace std;

// Function to count the number of unique nodes, and assign a group integer to each node
void PlasmaData::update_group(void)
{
	// There has got to be a better way of doing this, but well this should work for now...

	char processor_names[num_nodes][MPI_MAX_PROCESSOR_NAME];
	int processor_name_lengths[num_nodes];

	char my_processor_name[MPI_MAX_PROCESSOR_NAME];
	int my_processor_name_length;

	// Every node gets its processor name
	printf("Getting Names for %i nodes\n",num_nodes);
	MPI_Get_processor_name(my_processor_name,&my_processor_name_length);
	printf("My name is %s\n",my_processor_name);
	printf("Reducing Name Lengths\n");
	if(mynode != 0)
	{
		// send your processor name length to the root node
		MPI_Send(&my_processor_name_length,1,
				MPI_INT,0,mynode,MPI_COMM_WORLD);
	}
	else
	{
		// recieve the processor name lengths
		for(int i=1;i<num_nodes;i++)
		{
			MPI_Recv(&processor_name_lengths[i],1,
					MPI_INT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	printf("Reducing Names\n");
	if(mynode != 0)
	{
		// send your processor name to the root node
		MPI_Send(my_processor_name,MPI_MAX_PROCESSOR_NAME,
				MPI_CHAR,0,mynode+num_nodes,MPI_COMM_WORLD);
	}
	else
	{
		strcpy(processor_names[0], my_processor_name);
		processor_name_lengths[0] = my_processor_name_length;
		// recieve the processor names
		for(int i=1;i<num_nodes;i++)
		{
			MPI_Recv(processor_names[i],MPI_MAX_PROCESSOR_NAME,
					MPI_CHAR,i,i+num_nodes,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

//			printf("Reciveing name[%i]: %s\n",mynode,processor_names[i]);

		}
	}


	MPI_Barrier(MPI_COMM_WORLD);

	printf("Broadcasting Names\n");

	// Broadcast all of the names
	for(int i=0;i<num_nodes;i++)
	{
		MPI_Bcast(processor_name_lengths+i,1,MPI_INT,0,MPI_COMM_WORLD);

	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(int i=0;i<num_nodes;i++)
	{

		MPI_Bcast(processor_names[i],MPI_MAX_PROCESSOR_NAME,
				MPI_CHAR,0,MPI_COMM_WORLD);

	}


	MPI_Barrier(MPI_COMM_WORLD);

	const char* unique_processor_names0[num_nodes];
	int unique_processor_name_lengths0[num_nodes];
	int num_tasks_per_processor[num_nodes];
	int num_unique_procs;
	int k = 0;
	if(mynode == 0){
		printf("Sorting Names\n");
	list<string> unique_procs;
	list<string> unique_procs2;
	for(int i=0;i<num_nodes;i++)
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
	num_procs = k;

	// Broadcast the number of names to everyone

	MPI_Bcast(&num_procs,1,MPI_INT,0,MPI_COMM_WORLD);


	char unique_processor_names[num_procs][MPI_MAX_PROCESSOR_NAME];
	int unique_processor_name_lengths[num_procs];


	for(int i=0;i<num_procs;i++)
	{
		if(mynode == 0)
//		printf("Broadcasting length of names = %i\n",unique_processor_name_lengths0[i]);
		unique_processor_name_lengths[i] = unique_processor_name_lengths0[i];
		MPI_Bcast(unique_processor_name_lengths+i,1,MPI_INT,0,MPI_COMM_WORLD);

		if(mynode == 0)
		{
			strcpy(unique_processor_names[i],unique_processor_names0[i]);
//			printf("Broadcasting name: %s\n",unique_processor_names[i]);
		}
		MPI_Bcast(unique_processor_names[i],MPI_MAX_PROCESSOR_NAME,
				MPI_CHAR,0,MPI_COMM_WORLD);

	}

	MPI_Barrier(MPI_COMM_WORLD);

	int* group_members[num_procs];
	int num_group_members[num_procs];

	if(mynode < num_procs)
	{
//		printf("Grouping names %i\n",mynode);
		// since this is parallelizable we might as well

		// count the number of tasks with a given processor name
		num_group_members[mynode] = 0;
		string str1,str2;
		str1 = string(unique_processor_names[mynode]);
		for(int i=0;i<num_nodes;i++)
		{
			str2 = string(processor_names[i]);

			cout << str1 << " " << str2 << endl;

			if(str1 == str2)
				num_group_members[mynode]++;

		}

		// allocate the number of groups members
		group_members[mynode] = (int*)malloc(num_group_members[mynode]*sizeof(int));

		printf("group %i has %i members\n",mynode,num_group_members[mynode]);
		// Setup the list of group members
		str1 = unique_processor_names[mynode];
		k = 0;
		for(int i=0;i<num_nodes;i++)
		{
			str2 = processor_names[i];

			if(str1 == str2)
			{
				group_members[mynode][k] = i;
				k++;
			}

		}




	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Broadcast all of the group members
	printf("Broadcasting group members\n");
	for(int i=0;i<num_procs;i++)
	{
		MPI_Bcast(&num_group_members[i],1,MPI_INT,i,MPI_COMM_WORLD);

		if(mynode != i)
			group_members[i] = (int*)malloc(num_group_members[i]*sizeof(int));
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(int i=0;i<num_procs;i++)
	{
		MPI_Bcast(group_members[i],num_group_members[i],
				MPI_INT,i,MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	node_groups = (MPI_Group*)malloc(num_procs*sizeof(MPI_Group));
	node_comms = (MPI_Comm*)malloc(num_procs*sizeof(MPI_Comm));

	// Create the MPI Groups
	MPI_Group orig_group;

	// Get the original group
	MPI_Comm_group(MPI_COMM_WORLD,&orig_group);
	MPI_Barrier(MPI_COMM_WORLD);
	// Divide tasks into processor based groups

	for(int i=0;i<num_procs;i++)
	{
		bool ibelong = false;

		for(int j=0;j<num_group_members[i];j++)
		{
			if(mynode == group_members[i][j])
			{
				ibelong = true;
				my_proc = i;
			}

		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Group_incl(orig_group,num_group_members[i],group_members[i],node_groups+i);

		if(ibelong)
		{

			if(mynode == group_members[i][0])
			{
				printf("Node %i Creating group %i\n",mynode,i);
			}
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);
	// Create MPI Communicators
	for(int i=0;i<num_procs;i++)
	{
		if(mynode == group_members[i][0])
		printf("Creating comm for group %i\n",i);
		MPI_Comm_create(MPI_COMM_WORLD,node_groups[i],node_comms+i);
	}

	// Get the on node rank
	MPI_Comm_rank(node_comms[my_proc],&myrank_node);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("rank = %i, node = %i, node rank = %i\n",mynode,my_proc,myrank_node);

//	getchar();
	MPI_Barrier(MPI_COMM_WORLD);

}




void PlasmaData::readInput(const std::string& inFile)
{
  ifstream inStr(inFile.c_str());
  if (!inStr) {
    cout << "Could not open input file " << inFile << endl;
    exit(-1);
  }

  char inBuf[LINESIZE];
  string keyword;
  string rest;

  while (inStr.getline(inBuf, LINESIZE)) {
    if (inBuf[0] != '#' && inStr.gcount() > 1) {

      getKeyword(inBuf, keyword, rest);
      istringstream line(rest.c_str());
      cout << "KEYWORD: " << keyword << "   " << rest << endl;

      /*
      if (keyword == "PROBLEM_TYPE")
        line >> this->prob_type;
      else if (keyword == "PROBLEM_CASE")
        line >> this->prob_case;
      else if (keyword == "QUIET_START")
        line >> this->quiet_start_flag;
      else if (keyword == "SYSTEM_LENGTH")
        line >> this->lx;
      else if (keyword == "NUMBER_OF_POS_CELLS")
        line >> this->nx;
      else if (keyword == "NUMBER_OF_VEL_CELLS")
        line >> this->nv;

      else if (keyword == "TIME_MAX")
        line >> this->tmax;
      else if (keyword == "TIME_STEP")
        line >> this->dt;

      else if (keyword == "NUMBER_OF_FILTERING_OPERATIONS")
        line >> this->fil_num;
      else if (keyword == "CHARGE_CONSERVE_FLAG")
        line >> this->cc_flag;
      else if (keyword == "LO_SOLVER")
        line >> this->si_lo_flag;
      else if (keyword == "TEMPORAL_ORDER")
        line >> this->temp_order;
      else if (keyword == "W_SCALE")
        line >> this->w_scale;

      else if (keyword == "NUMBER_OF_PARTICLE_SPECIES")
        line >> this->p_size;
      else if (keyword == "PARTICLE_TYPE")
        line >> this->p_type[0] >> this->p_type[1];
      else if (keyword == "PARTICLE_REFERENCE_NUMBER")
        line >> this->NP_ref[0] >> this->NP_ref[1];
      else if (keyword == "SPECIES_CHARGE")
        line >> this->q[0] >> this->q[1];
      else if (keyword == "SPECIES_MASS")
        line >> this->m[0] >> this->m[1];
      else if (keyword == "SPECIES_SUBCYCLE")
        line >> this->sub_nt[0] >> this->sub_nt[1];

      else if (keyword == "DENSITY_UNPERTURBED")
        line >> this->rho0[0] >> this->rho0[1];
      else if (keyword == "VELOCITY_UNPERTURBED")
        line >> this->u0[0] >> this->u0[1];
      else if (keyword == "TEMPERATURE_UNPERTURBED")
        line >> this->T0[0] >> this->T0[1];
      else if (keyword == "DENSITY_PERTURBED")
        line >> this->alp_r[0] >> this->alp_r[1];
      else if (keyword == "VELOCITY_PERTURBED")
        line >> this->alp_u[0] >> this->alp_u[1];
      else if (keyword == "TEMPERATURE_PERTURBED")
        line >> this->alp_T[0] >> this->alp_T[1];
        */



      if (keyword == "SYSTEM_LENGTH_X")
        line >> this->Lx;
      else if (keyword == "SYSTEM_LENGTH_Y")
        line >> this->Ly;
      else if (keyword == "SYSTEM_LENGTH_Z")
        line >> this->Lz;

      else if (keyword == "NUMBER_OF_POS_CELLS_X")
        line >> this->nx;
      else if (keyword == "NUMBER_OF_POS_CELLS_Y")
        line >> this->ny;
      else if (keyword == "NUMBER_OF_POS_CELLS_Z")
        line >> this->nz;

      else if (keyword == "X_ORIGIN")
        line >> this->xmin;
      else if (keyword == "Y_ORIGIN")
        line >> this->ymin;
      else if (keyword == "Z_ORIGIN")
        line >> this->zmin;

      else if (keyword == "NUMBER_OF_TIME_STEPS")
        line >> this->nsteps;
      else if (keyword == "TIME_STEP")
        line >> this->dt;

      else if (keyword == "NUMBER_OF_PARTICLES")
        line >> this->nptcls;
      else if (keyword == "NUMBER_OF_PARTICLE_SPECIES")
        line >> this->nspecies;

      else if (keyword == "CPU_VECTOR_LENGTH")
        line >> this->cpu_vec_length;
      else if (keyword == "PICARD_TOLERANCE")
        line >> this->epsilon_r;


    }
  }
}

void PlasmaData::getKeyword(char* inBuf, string& keyword, string& rest)
{
  string line(inBuf);
  string::size_type keyPos = line.find(' ');
  keyword = line.substr(0, keyPos);
  rest = line.substr(keyPos + 1);
}

