#include "ParticleListGPU.cuh"
#include "CurrentTally.h"
#include "HOMoments.h"
#include "PlasmaData.h"
#include "ProblemInitializer.h"
#include "math.h"
#include "omp.h"




ParticleListGPU::ParticleListGPU()
{

}

ParticleListGPU::~ParticleListGPU()
{
	/*
	// Free realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		free(*get_float(i));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		free(*get_int(i));
	}

	// allocate short ints for cluster id's
	free(cluster_id);
	*/
}

void ParticleListGPU::copy_from(const ParticleList* list_in)
{
	ispecies = list_in -> ispecies;
	// Free realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		memcpy(*get_float(i),*(list_in->get_float(i)),nptcls*sizeof(realkind));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		memcpy(*get_int(i),*(list_in->get_int(i)),nptcls*sizeof(int));
	}

	// allocate short ints for cluster id's
	memcpy(cluster_id,list_in->cluster_id,nptcls*sizeof(short int));
}

void ParticleListGPU::allocate(PlasmaData* pdata,int nptcls_in)
{
	// Allocate memory for particles
	nptcls_allocated = nptcls_in;

	nptcls = nptcls_in;

	// Allocate realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		*get_float(i) = (realkind*)malloc(nptcls_allocated*sizeof(realkind));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		*get_int(i) = (int*)malloc(nptcls_allocated*sizeof(int));
	}

	// allocate short ints for cluster id's
	cluster_id = (short int*)malloc(nptcls_allocated*sizeof(short int));

}

void ParticleListGPU::init(ProblemInitializer* initializer, HOMoments* moments)
{

	CurrentTally currents(moments->currentx,
						  moments->currenty,
						  moments->currentz,
						  make_int3(moments->pdata->nx,moments->pdata->ny,moments->pdata->nz),
						  moments->pdata->dxdi,moments->pdata->dydi,moments->pdata->dzdi,
						  moments->pdata->ndimensions);

	moments -> set_vals(0);

	for(int i=0;i<nptcls;i++)
	{
		realkind px,py,pz,vx,vy,vz;
		int ix,iy,iz;


		initializer->init_particle(px,py,pz,ix,iy,iz,vx,vy,vz,ispecies,i);




		currents.tally(px,py,pz,vx,vy,vz,ix,iy,iz,1.0f);


		// Set Position Values, ifloat = 0-2
		this->get_fvalue(i,0) = px;
		this->get_fvalue(i,1) = py;
		this->get_fvalue(i,2) = pz;

		// Set Position Index Values, iint = 0-2
		this->get_ivalue(i,0) = ix;
		this->get_ivalue(i,1) = iy;
		this->get_ivalue(i,2) = iz;

		// Set Velocity Values, ifloat = 3-5
		this->get_fvalue(i,3) = vx;
		this->get_fvalue(i,4) = vy;
		this->get_fvalue(i,5) = vz;
	}

}

realkind ParticleListGPU::evaluate_energy(PlasmaData* pdata)
{
	double etotal = 0.0;
	for(int i=0;i<nptcls;i++)
	{
		etotal += get_fvalue(i,3)* get_fvalue(i,3);
		etotal += get_fvalue(i,4)* get_fvalue(i,4);
		etotal += get_fvalue(i,5)* get_fvalue(i,5);
	}

	etotal = etotal * pdata->mspecies[ispecies] * 0.5/((double)pdata->nptcls_total);

	return etotal;
}

void ParticleListGPU::CPUfree()
{
	// Allocate realkind arrays
	for(int i=0;i<ParticleList_nfloats;i++)
	{
		free(*get_float(i));
	}

	// Allocate int arrays
	for(int i=0;i<ParticleList_nints;i++)
	{
		free(*get_int(i));
	}

	// allocate short ints for cluster id's
	free(cluster_id);
}

long long int ParticleListGPU::push(PlasmaData* pdata, FieldData* fields, HOMoments* moments)
{
	return 0;
}

void ParticleListGPU::plot_particles(PlasmaData* pdata){}
