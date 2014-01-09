/*
 * NodeHOMoments.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: payne
 */

#include "omp.h"
#include "mpi.h"
#include "PlasmaData.h"
#include "NodeHOMoments.h"
#include "ParallelInfo.h"
#include "FieldDataCPU.h"







void NodeHOMoments::allocate(PlasmaData* _pdata)
{
	pdata = _pdata;

	nx = pdata->nx;
	ny = pdata->ny;
	nz = pdata->nz;
	nspecies = pdata->nspecies;

	moments_old = new HOMomentsCPU(pdata);
	moments_next = new HOMomentsCPU(pdata);

	omp_set_num_threads(pdata->node_info->nCPU);

	cpu_moments = (HOMomentsCPU*)malloc(pdata->node_info->nThreads*sizeof(HOMomentsCPU));

// This only sets up the HOMoments for the CPU, the other CPU HOMoments for
// the GPU and MIC must be set up when those device specific versions are set up
	printf("Allocating CPU Moments\n");
	   omp_set_dynamic(0);
	   omp_set_nested(0);
	for(int i=0;i<pdata->node_info->nCPU;i++)
		   cpu_moments[i] = *(new HOMomentsCPU(pdata));


		printf("Allocating GPU Moments\n");
	if(pdata->node_info->nGPU > 0)
	{
		gpu_moments = (HOMomentsGPU*)malloc(pdata->node_info->nGPU*sizeof(HOMomentsGPU));
		for(int i=0;i<pdata->node_info->nGPU;i++)
		{
			CUDA_SAFE_CALL(cudaSetDevice(pdata->thread_info[pdata->node_info->nspecies+i]->gpu_info->igpu));
			gpu_moments[i] = *(new HOMomentsGPU(pdata));
			cpu_moments[pdata->node_info->nCPU+i] = *(new HOMomentsCPU(pdata));
		}
	}

	if(pdata->node_info->nMIC > 0)
	{
		mic_moments = (HOMomentsMIC*)malloc(pdata->node_info->nMIC*sizeof(HOMomentsMIC));
		for(int i=0;i<pdata->node_info->nMIC;i++)
		{
			mic_moments[i] = *(new HOMomentsMIC(pdata));
			cpu_moments[pdata->node_info->nCPU+pdata->node_info->nGPU+i] = *(new HOMomentsCPU(pdata));
		}
	}






}

void NodeHOMoments::reduce()
{

	int alloc_size = nx*ny*nz*nspecies*10;

	omp_set_num_threads(pdata->node_info->nCPU);
	int nthreads = pdata->node_info->nCPU;

//#pragma omp critical
//	printf("Current number of omp threads = %i\n",omp_get_num_threads());
//	for(int i=0;i<nx;i++)
//	{
//		printf("CPU Moments[%i] = %e\n",cpu_moments[0].get_val(i,0,0,0,HOMoments_currentx));
//	}


#pragma omp parallel for
	for(int i=0;i<alloc_size;i++)
	{


		for(int j=1;j<pdata->node_info->nThreads;j++)
		{
//#pragma omp critical
			cpu_moments[0].all_data[i] += cpu_moments[j].all_data[i];
		}
	}



	moments_next -> copy_from(cpu_moments);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Allreduce(cpu_moments->all_data,moments_next->all_data,alloc_size,
			MPI_REALKIND,MPI_SUM,MPI_COMM_WORLD);
}


void NodeHOMoments::apply_weights()
{

	moments_next -> apply_weights();

}

void NodeHOMoments::UpdateMoments()
{

	moments_old -> copy_from(moments_next);

}

realkind& NodeHOMoments::get_val(int i, int j, int k,
		int ispecies,enum HOMoments_moment moment,int itime)
{
	switch(itime)
	{
	case 0:
		return moments_old -> get_val(i,j,k,ispecies,moment);
		break;
	case 1:
		return moments_next -> get_val(i,j,k,ispecies,moment);
		break;
	default:
		break;
	}

	return moments_old -> get_val(i,j,k,ispecies,moment);
}


void NodeHOMoments::set_vals(realkind val_in)
{



#pragma omp parallel for
	for(int i=0;i<pdata->node_info->nThreads;i++)
	{
//		for(int h=0;h<nx;h++)
//		{
//			printf("Currentx[%i][%i] = %e\n",i,h,cpu_moments[i].get_val(h,0,0,0,HOMoments_currentx));
//		}
		cpu_moments[i].set_vals(val_in);
	}

	if(pdata->node_info->nGPU > 0)
	{
#pragma omp parallel for
		for(int i=0;i<pdata->node_info->nGPU;i++)
		{
			CUDA_SAFE_CALL(cudaSetDevice(pdata->thread_info[pdata->node_info->nspecies+i]->gpu_info->igpu));
			gpu_moments[i].set_vals(val_in);
		}
	}

	if(pdata->node_info->nMIC > 0)
	{
//#pragma omp for
		for(int i=0;i<pdata->node_info->nMIC;i++)
		{
			mic_moments[i].set_vals(val_in);
		}
	}
}

double NodeHOMoments::evaluate_energy()
{
	return moments_next->evaluate_energy();
}

__host__
double NodeHOMoments::evaluate_energy(int iS)
{
	return moments_next->evaluate_energy(iS);
}



double NodeHOMoments::check_charge()
{

	double charge_cons = 0;

	int j=0;
	int k=0;

		realkind dv = pdata->didx*pdata->didy*pdata->didz;
	if(pdata->ndimensions == 1)
	{
		for(int i=0;i<pdata->nx;i++)
		{
			double deltaj  = 0;
			double deltan = 0;
			for(int l=0;l<pdata->nspecies;l++)
			{



				deltaj += pdata->qspecies[l]*(get_val(i+1,j,k,l,HOMoments_currentx,1)-get_val(i,j,k,l,HOMoments_currentx,1))*pdata->didx;
				deltan += pdata->qspecies[l]*(get_val(i,j,k,l,HOMoments_charge,1)-get_val(i,j,k,l,HOMoments_charge,0))/pdata->dt;



		//	printf("charge_cons(%i) = %e + %e = %e\n",i,deltaj,deltan,deltaj+deltan);

			}
			charge_cons += pow((deltan + deltaj),2);

		}
		charge_cons = sqrt(charge_cons/((double)pdata->nx));
	}
	else
	{
		for(int j=0;j<pdata->ny;j++)
		{
			for(int i=0;i<pdata->nx;i++)
			{
				realkind deltaj  = 0;
				realkind deltan = 0;

				for(int l=0;l<pdata->nspecies;l++)
				{
					deltaj += pdata->qspecies[l]*(get_val(i+1,j,k,l,HOMoments_currentx)-get_val(i,j,k,l,HOMoments_currentx))*pdata->didx;
					deltan += pdata->qspecies[l]*(get_val(i,j,k,l,HOMoments_charge)-moments_old->get_val(i,j,k,l,HOMoments_charge))/(pdata->dt);

					if(pdata->ndimensions >1)
					deltaj += pdata->qspecies[l]*(get_val(i,j+1,k,l,HOMoments_currenty)-get_val(i,j,k,l,HOMoments_currenty))*pdata->didy;

				}

			//	printf("charge_cons(%i) = %e + %e = %e\n",i,deltaj,deltan,deltaj+deltan);

				charge_cons += pow((deltan + deltaj),2);
			}
		}

		charge_cons = sqrt(charge_cons/((double)pdata->nx*pdata->ny));

	}


	return charge_cons;

}

double NodeHOMoments::check_momentum(FieldDataCPU* fields_next, FieldDataCPU* fields_old)
{

	double mom_cons = 0;

	int j=0;
	int k=0;

		realkind dv = pdata->didx*pdata->didy*pdata->didz;
	if(pdata->ndimensions == 1)
		for(int l=0;l<pdata->nspecies;l++)
		{
		for(int i=0;i<pdata->nx;i++)
		{
			double deltaVy  = 0;
			double deltaAy = 0;

			double deltaVz  = 0;
			double deltaAz = 0;


				deltaVy = pdata->mspecies[l]*(get_val(i,j,k,l,HOMoments_currenty,1)-get_val(i,j,k,l,HOMoments_currenty,0));
				deltaAy = pdata->qspecies[l]*(
						0.5*(get_val(i,j,k,l,HOMoments_charge,1)+get_val(i-1,j,k,l,HOMoments_charge,1))*fields_next->getA(i,j,k,1)
						- 0.5*(get_val(i,j,k,l,HOMoments_charge,0)+get_val(i-1,j,k,l,HOMoments_charge,0))*fields_old->getA(i,j,k,1));

				deltaVz = pdata->mspecies[l]*(get_val(i,j,k,l,HOMoments_currentz,1)-get_val(i,j,k,l,HOMoments_currentz,0));
				deltaAz = pdata->qspecies[l]*(
						0.5*(get_val(i,j,k,l,HOMoments_charge,1)+get_val(i-1,j,k,l,HOMoments_charge,1))*fields_next->getA(i,j,k,2)
						- 0.5*(get_val(i,j,k,l,HOMoments_charge,0)+get_val(i-1,j,k,l,HOMoments_charge,0))*fields_old->getA(i,j,k,2));
				double temp = (abs(deltaVy+deltaAy)+abs(deltaVz + deltaAz))/((double)pdata->nptcls_species[l]);
		//	printf("charge_cons(%i) = %e + %e = %e\n",i,deltaj,deltan,deltaj+deltan);

				mom_cons = fmax(temp,mom_cons);
		}
		}



	return mom_cons;

}





