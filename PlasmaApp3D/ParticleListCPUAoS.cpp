#include "ParticleListCPUAoS.h"


void ParticleListCPUAoS::allocate(PlasmaData* pdata,int nptcls_in)
{
	// Allocate memory for particles
	nptcls_allocated = nptcls_in;

	nptcls = nptcls_in;

	particles = (CPUParticleData*)malloc(nptcls_allocated*sizeof(CPUParticleData));

	// allocate short ints for cluster id's
	cluster_id = (int*)malloc(nptcls_allocated*sizeof(short int));
}


void ParticleListCPUAoS::CPUfree()
{
	// Free Particle array
	free(particles);

	// allocate short ints for cluster id's
	free(cluster_id);

	nptcls_allocated = 0;
	nptcls = 0;
}

void ParticleListCPUAoS::copy_from(const ParticleList* list_in)
{
	ispecies = list_in -> ispecies;
	// copy realkind data
	for(int i=0;i<ParticleList_nfloats-1;i++)
	{
		for(int j=0;j<nptcls;j++)
		{
			get_fvalue(j,i) = list_in->get_fvalue(j,i);
		}
	}

	// copy int data
	for(int i=0;i<ParticleList_nints;i++)
	{
		for(int j=0;j<nptcls;j++)
		{
			get_ivalue(j,i) = list_in->get_ivalue(j,i);
		}
	}

	// copy short ints for cluster id's
	memcpy(cluster_id,list_in->cluster_id,nptcls*sizeof(short int));
}

void ParticleListCPUAoS::copy_from(const ParticleListCPUAoS* list_in)
{
	ispecies = list_in -> ispecies;
	// copy realkind data
	memcpy(particles,list_in->particles,nptcls*sizeof(CPUParticleData));

	// copy short ints for cluster id's
	memcpy(cluster_id,list_in->cluster_id,nptcls*sizeof(short int));
}
