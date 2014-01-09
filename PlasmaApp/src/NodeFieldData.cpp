/*
 * NodeFieldData.cpp
 *
 *  Created on: May 1, 2013
 *      Author: payne
 */

#include "NodeFieldData.h"
#include "PlasmaData.h"
#include "CPUTimer.h"
#include "ParallelInfo.h"
#include "RunData.h"
#include <omp.h>
#include <mpi.h>

void NodeFieldData::allocate(PlasmaData* _pdata)
{
	pdata = _pdata;

	nx = pdata->nx;
	ny = pdata->ny;
	nz = pdata->nz;


	cpu_fields = new FieldDataCPU();
	cpu_fields -> allocate(pdata);


	if(pdata->node_info->nGPU > 0)
	{
		gpu_fields = (FieldDataGPU*)malloc(pdata->node_info->nGPU * sizeof(FieldDataGPU));
#pragma omp parallel for
		for(int i=0;i<pdata->node_info->nGPU;i++)
		{
			CUDA_SAFE_CALL(cudaSetDevice(pdata->thread_info[pdata->node_info->nspecies+i]->gpu_info->igpu));
			gpu_fields[i] = *(new FieldDataGPU());
			gpu_fields[i].allocate(pdata);
		}
	}

	if(pdata->node_info->nMIC > 0)
	{
		mic_fields = new FieldDataMIC();
		mic_fields -> allocate(pdata);
	}
	bcast_timer = new CPUTimer();
}



void NodeFieldData::allocate(PlasmaData* _pdata, NodeFieldData* _fields)
{
	pdata = _pdata;

	nx = pdata->nx;
	ny = pdata->ny;
	nz = pdata->nz;

	alloc_size = nx*ny*nz;


	cpu_fields = new FieldDataCPU();
	cpu_fields -> allocate(pdata);


	if(pdata->node_info->nGPU > 0)
	{
		gpu_fields = (FieldDataGPU*)malloc(pdata->node_info->nGPU * sizeof(FieldDataGPU));
#pragma omp parallel for
		for(int i=0;i<pdata->node_info->nGPU;i++)
		{
			CUDA_SAFE_CALL(cudaSetDevice(pdata->thread_info[pdata->node_info->nspecies+i]->gpu_info->igpu));
			gpu_fields[i] = *(new FieldDataGPU());
			gpu_fields[i].allocate(pdata);
		}
	}

	if(pdata->node_info->nMIC > 0)
	{
		mic_fields = new FieldDataMIC();
		mic_fields -> allocate(pdata);
	}

	bcast_timer = _fields->bcast_timer;

}

void NodeFieldData::average(FieldDataCPU* field0,FieldDataCPU* field1)
{




	if(pdata->ndimensions == 1)
	{

//		for(int k=0;k<nz;k++)
////#pragma omp parallel for
//		for(int j=0;j<ny;j++)
		int j=0;
		int k=0;
#pragma omp parallel for
		for(int i=0;i<nx;i++)
		{
			for(int l=0;l<3;l++)
			{
				cpu_fields->getE(i,j,k,l) = 0.125*(field0->getE(i+1,j,k,l)
						+ 2.0*field0->getE(i,j,k,l)
						+ field0->getE(i-1,j,k,l)
						+ field1->getE(i+1,j,k,l)
						+ 2.0*field1->getE(i,j,k,l)
						+ field1->getE(i-1,j,k,l));

				cpu_fields->getB(i,j,k,l) = 0.125*(field0->getB(i+1,j,k,l)
						+ 2.0*field0->getB(i,j,k,l)
						+ field0->getB(i-1,j,k,l)
						+ field1->getB(i+1,j,k,l)
						+ 2.0*field1->getB(i,j,k,l)
						+ field1->getB(i-1,j,k,l));

				cpu_fields->getA(i,j,k,l) = 0.25*(field0->getA(i+1,j,k,l)
						+ 2.0*field0->getA(i,j,k,l)
						+ field0->getA(i-1,j,k,l));

				cpu_fields->getAp(i,j,k,l) = 0.25*(
						  field1->getA(i+1,j,k,l)
						+ 2.0*field1->getA(i,j,k,l)
						+ field1->getA(i-1,j,k,l));



//				cpu_fields->getE(i,j,k,l) = 0.5*(
//						+ field0->getE(i,j,k,l)
//						+ field1->getE(i,j,k,l));
//
//				cpu_fields->getB(i,j,k,l) = 0.5*(
//						+ field0->getB(i,j,k,l)
//						+ field1->getB(i,j,k,l));
			}

//			cpu_fields->getE(i,j,k,1) = 0.25*(field0->getA(i+1,j,k,1)
//					+ 2.0*field0->getA(i,j,k,1)
//					+ field0->getA(i-1,j,k,1)
//					- (field1->getA(i+1,j,k,1)
//					+ 2.0*field1->getA(i,j,k,1)
//					+ field1->getA(i-1,j,k,1)))/pdata->dt;
//
//			cpu_fields->getE(i,j,k,2) = 0.25*(field0->getA(i+1,j,k,2)
//								+ 2.0*field0->getA(i,j,k,2)
//								+ field0->getA(i-1,j,k,2)
//								- (field1->getA(i+1,j,k,2)
//								+ 2.0*field1->getA(i,j,k,2)
//								+ field1->getA(i-1,j,k,2)))/pdata->dt;
		}
	}
	else
	{

		for(int k=0;k<nz;k++)
				{
					for(int j=0;j<ny;j++)
					{
						for(int i=0;i<nx;i++)
						{
							for(int l=0;l<3;l++)
							{
								cpu_fields->getE(i,j,k,l) = 0.25*(
										0.125*(field0->getE(i+1,j-1,0,l)
										+ 2.0*field0->getE(i,j-1,0,l)
										+ field0->getE(i-1,j-1,0,l)
										+ field1->getE(i+1,j-1,0,l)
										+ 2.0*field1->getE(i,j-1,0,l)
										+ field1->getE(i-1,j-1,0,l)) +

										2.0*0.125*(field0->getE(i+1,j,0,l)
										+ 2.0*field0->getE(i,j,0,l)
										+ field0->getE(i-1,j,0,l)
										+ field1->getE(i+1,j,0,l)
										+ 2.0*field1->getE(i,j,0,l)
										+ field1->getE(i-1,j,0,l)) +

										0.125*(field0->getE(i+1,j+1,0,l)
										+ 2.0*field0->getE(i,j+1,0,l)
										+ field0->getE(i-1,j+1,0,l)
										+ field1->getE(i+1,j+1,0,l)
										+ 2.0*field1->getE(i,j+1,0,l)
										+ field1->getE(i-1,j+1,0,l))
										);

								cpu_fields->getB(i,j,k,l) = 0.25*(
										0.125*(field0->getB(i+1,j-1,0,l)
										+ 2.0*field0->getB(i,j-1,0,l)
										+ field0->getB(i-1,j-1,0,l)
										+ field1->getB(i+1,j-1,0,l)
										+ 2.0*field1->getB(i,j-1,0,l)
										+ field1->getB(i-1,j-1,0,l)) +

										2.0*0.125*(field0->getB(i+1,j,0,l)
										+ 2.0*field0->getB(i,j,0,l)
										+ field0->getB(i-1,j,0,l)
										+ field1->getB(i+1,j,0,l)
										+ 2.0*field1->getB(i,j,0,l)
										+ field1->getB(i-1,j,0,l)) +

										0.125*(field0->getB(i+1,j+1,0,l)
										+ 2.0*field0->getB(i,j+1,0,l)
										+ field0->getB(i-1,j+1,0,l)
										+ field1->getB(i+1,j+1,0,l)
										+ 2.0*field1->getB(i,j+1,0,l)
										+ field1->getB(i-1,j+1,0,l))
										);

								cpu_fields->getA(i,j,k,l) = 0.25*(field0->getA(i+1,j,0,l)
										+ 2.0*field0->getA(i,j,0,l)
										+ field0->getA(i-1,j,0,l));

								cpu_fields->getAp(i,j,k,l) = 0.25*(
										  field1->getA(i+1,j,0,l)
										+ 2.0*field1->getA(i,j,0,l)
										+ field1->getA(i-1,j,0,l));
							}
						}
					}
				}
//		for(int k=0;k<nz;k++)
//		{
//			for(int j=0;j<ny;j++)
//			{
//				for(int i=0;i<nx;i++)
//				{
//					for(int l=0;l<3;l++)
//					{
//						cpu_fields -> getE(i,j,k,l) = 0.5*(field0 -> getE(i,j,k,l)
//								+ field1 -> getE(i,j,k,l));
//
//						cpu_fields -> getB(i,j,k,l) = 0.5*(field0 -> getB(i,j,k,l)
//													+ field1 -> getB(i,j,k,l));
//					}
//				}
//			}
//		}
	}

}

void NodeFieldData::broadcast()
{

	bcast_timer->start();
	// Cluster level broadcast
	if(!(pdata->rdata->lo_all))
	{
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(cpu_fields->all_data,12*alloc_size,MPI_REALKIND,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

	}





	// Node level broadcast aka send the data to all devices (we can skip the CPU)
	omp_set_num_threads(pdata->node_info->nDevices);

#pragma omp parallel for
for(int i=0;i<pdata->node_info->nDevices2;i++)
	{
		int tid = omp_get_thread_num();

		int device_type = pdata->thread_info[i]->device_type;
		int device_id = pdata->thread_info[i] -> device_id;

		if(device_type == 1)
		{
			CUDA_SAFE_CALL(cudaSetDevice(pdata->thread_info[i]->gpu_info->igpu));
			gpu_fields[device_id].copy_from(cpu_fields);
		}
		else if(device_type == 2)
		{
			mic_fields[device_id].copy_from(cpu_fields);
		}
	}

	bcast_timer->stop();
}


double NodeFieldData::evaluate_energy()
{
	return cpu_fields->evaluate_energy();
}




























