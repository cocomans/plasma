#include "EnergyCons_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../HOMoments.h"
#include "../CurrentTally.h"
#include "../PlasmaData.h"
#include "../FieldData.h"
#include "../ParallelInfo.h"

EnergyCons_Initializer::EnergyCons_Initializer(PlasmaData* pdata_in)
{
		pdata = pdata_in;
		Ey = 2.0;
		Bz = 50;
		kt = 1.0*(pdata->Lx)/(2.0*3.14159265359);
		v_phase = 2.0*sqrt(2.0*pdata->Te/(pdata->mspecies[0]*mass_e));
		printf("Phase Velocity = %f\n",v_phase);

}

void EnergyCons_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
									 int& ix, int& iy, int& iz,
									 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl)
{
	// Set Position Values, ifloat = 0-2
	px = ((2.0*iptcl)/((double)pdata->nptcls))*pdata->Lx+pdata->xmin;
	vx = box_muller(0.0,sqrt(2.0*pdata->Te/(pdata->mspecies[0]*mass_e)));

	py = 0.5;
	pz = 0.5;



	ix = floor(px*pdata->didx);
	iy = floor(py*pdata->didy);
	iz = floor(pz*pdata->didz);

	px = px*pdata->didx - ix;
	py = py*pdata->didy - iy;
	pz = pz*pdata->didz - iz;

	ix = ((ix%pdata->nx)+pdata->nx)%pdata->nx;

	// Set Velocity Values, ifloat = 3-5
	vy = 0;

	vz = 0;




}


void EnergyCons_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{
	// Set Velocity Values, ifloat = 3-5
	if(ispecies == 0)
		vx = 0.1*(1-(2*(rand()%2)));
	else
		vx = 0;
	vy = 0;

	vz = 0;



}

void EnergyCons_Initializer::initialize_particles(ParticleList* particles,
										  HOMoments* moments,
										  ParallelInfo* myinfo)
{

	particles -> ispecies = myinfo->myid_mpi;

	particles -> init(this,moments);

	kE0 = particles->evaluate_energy(pdata);



}

void EnergyCons_Initializer::initialize_fields(FieldData** fields,
										HOMoments** moments,
										ParallelInfo* myinfo)
{

	int tid = myinfo->tid;

	// Setup E-field
	for(int k=0;k<pdata->nz;k++)
	{
		for(int j=0;j<pdata->ny;j++)
		{
			for(int i=0;i<pdata->nx;i++)
			{
				realkind x = (i)*pdata->dxdi+pdata->xmin;
				//realkind y = j*pdata->dydi+pdata->ymin;
				//realkind z = k*pdata->dzdi+pdata->zmin;

				realkind Ex = -0.1*sin(x/kt)*(kt);



				fields[tid+1] -> getE(i,j,k,0) = 0;
				fields[tid+1] -> getE(i,j,k,1) = 0;
				fields[tid+1] -> getE(i,j,k,2) = 0;

				fields[tid+1] -> getB(i,j,k,0) = 0;
				fields[tid+1] -> getB(i,j,k,1) = 0;
				fields[tid+1] -> getB(i,j,k,2) = 0;

				fields[0] -> getE(i,j,k,0) = 0;
				fields[0] -> getE(i,j,k,1) = 0;
				fields[0] -> getE(i,j,k,2) = 0;

				fields[0] -> getB(i,j,k,0) = 0;
				fields[0] -> getB(i,j,k,1) = 0;
				fields[0] -> getB(i,j,k,2) = 0;



			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

	if(myinfo->myid_mpi == 0)
	{
		float gamma;

		gamma = pi_const * powf(pdata->omega_pe,3.0f);

		plot = gnuplot_init();
		kenergy_plot = gnuplot_init();
		Tenergy_plot = gnuplot_init();
		Energy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
		kEnergy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
		TEnergy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));

		max_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
		max_t_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
		nplot = 0;
		nmax_array = 0;


		E0 = fields[0] -> evaluate_energy();
		Energy_array[0] = (E0);

		kEnergy_array[0] = kE0;

		TEnergy_array[0] = E0+kE0;
		nplot++;

		//gnuplot_cmd(plot,"set log y");
		//gnuplot_cmd(kenergy_plot,"set log y");


	}

}

void EnergyCons_Initializer::check_step(ParticleList* particles,
		HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo)
{



	// Print out the sum of the x-component of the current
	if(myinfo->myid_mpi == 0)
	{
		Energy_array[nplot] = ((fields[0] -> evaluate_energy()));
		kEnergy_array[nplot] = particles->evaluate_energy(pdata);

		TEnergy_array[nplot] = kEnergy_array[nplot] + Energy_array[nplot];
		printf("Energy Error= %e\n",(E0+kE0 - (Energy_array[nplot]+kEnergy_array[nplot]))/(E0+kE0));

		if(nplot == 1)
		{
			max_array[0] = log(Energy_array[0]);
			max_t_array[nmax_array] = pdata->dt;
			nmax_array++;
		}
		else if(nplot > 1)
		{
			if(Energy_array[nplot-1] > fmax(Energy_array[nplot-2],Energy_array[nplot]))
			{
				max_array[nmax_array] = log(Energy_array[nplot-1]);
				max_t_array[nmax_array] = (nplot-1);

				float period = (max_t_array[nmax_array]-max_t_array[nmax_array-1])*pdata->dt;
				float vphase = 2.0*pi_const/period* 1.0*pdata->Lx/(2.0*pi_const);
				printf("Energy fluctuation period = %f\n",period);
				printf("Phase velocity = %f\n",vphase);
				printf("Energy decay rate = %f\n",(max_array[nmax_array]-max_array[nmax_array-1])/
						(pdata->dt*(max_t_array[nmax_array]-max_t_array[nmax_array-1])));

				nmax_array++;
			}
		}


		if(nplot > 0)
		{
			gnuplot_resetplot(plot);
			gnuplot_resetplot(kenergy_plot);
			gnuplot_resetplot(Tenergy_plot);
		}

		nplot++;
		gnuplot_plot_x(plot,Energy_array,nplot,"Field Energy");
		gnuplot_plot_x(kenergy_plot,kEnergy_array,nplot,"Particle Energy");
		gnuplot_plot_x(Tenergy_plot,TEnergy_array,nplot,"Total Energy");
		//gnuplot_setstyle(plot,"lines");
		gnuplot_plot_xy(plot,max_t_array,max_array,nmax_array,NULL);
		//gnuplot_setstyle(plot,"points");


	}
















}

