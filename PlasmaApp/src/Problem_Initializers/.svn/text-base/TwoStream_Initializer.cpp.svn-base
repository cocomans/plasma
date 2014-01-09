#include "TwoStream_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../NodeHOMoments.h"
#include "../NodeParticleList.h"
#include "../FieldDataCPU.h"
#include "../PlasmaData.h"

#include "../ParallelInfo.h"

void TwoStream_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
		 int& ix, int& iy, int& iz,
		 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl,
		 int nptcls,int ioffset)
{
	// Set Position Values, ifloat = 0-2

	raised_cos cdf(alpha,1.0/kt);

	int iptcl_new = iptcl;

	int idist = iptcl_new%(pdata->nptcls_device[0]*pdata->node_info->nTasks_g/2);
	int ivel = iptcl_new/(pdata->nptcls_device[0]*pdata->node_info->nTasks_g/2);

//	idist += ioffset;
	idist = iptcl_new%(nptcls/2);
	ivel = iptcl_new/(nptcls/2);
	double temp_x;

	if(ispecies == 0)
		temp_x = distribution_intrp((idist)/(0.5*nptcls),0.0,pdata->Lx,cdf);
	else
		temp_x = (iptcl/100.0)/(0.01*nptcls) * pdata->Lx + pdata->xmin;

//	int iptcl_new = iptcl;
//
//	int idist = iptcl_new%(pdata->nptcls/2);
//	int ivel = iptcl_new/(pdata->nptcls/2);
//
//	double temp_x;
//
//	if(ispecies == 0)
//		temp_x = distribution_intrp((idist+0.5)/(0.5*pdata->nptcls),0.0,pdata->Lx,cdf);
//	else
//		px = iptcl/(1.0*pdata->nptcls) * pdata->Lx + pdata->xmin;
	//px = box_muller(0.5,0.1);

	///px = ((2.0*iptcl)/((double)pdata->nptcls))*pdata->Lx+pdata->xmin;
	//vx = 1.0*(2*(ivel%2)-1)*(exp(-(ivel/2)/100.0));

	//vx = sqrt(2.0*pdata->Te/(pdata->mspecies[0]*mass_e))
	//	*sqrt(-2.0*log((ivel2+0.5)/((realkind)nv2)))*cos(2.0*3.14159265359*(ivel1)/((realkind)nv1));

	//temp_x = drandom();

	//printf("random number = %f\n",px);

	// Set Velocity Values, ifloat = 3-5
	if(ispecies == 0)
		vx = 0.1*(1-(2*(ivel)));
	else
		vx = box_muller(0.0,0.01/sqrt(pdata->mspecies[ispecies]));


	py = 0.0;
	pz = 0.0;

	ix = floor(temp_x*pdata->didx);
	iy = floor(py*pdata->didy);
	iz = floor(pz*pdata->didz);

	px = temp_x*pdata->didx - ix;
	py = py*pdata->didy - iy;
	pz = pz*pdata->didz - iz;

	//px = drandom();

	ix = ((ix%pdata->nx)+pdata->nx)%pdata->nx;


}

void TwoStream_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{
	// Set Velocity Values, ifloat = 3-5
	if(ispecies == 0)
		vx = 0.1*(1-(2*(rand()%2)));
	else
		vx = 0;
	vy = 0;

	vz = 0;



}

void TwoStream_Initializer::initialize_particles(NodeParticleList* particles,
										  NodeHOMoments* moments)
{




	particles -> init(this,moments);
	MPI_Barrier(MPI_COMM_WORLD);


}

void TwoStream_Initializer::initialize_fields(FieldDataCPU* fields,
										NodeHOMoments* moments)
{


	// Setup E-field
	for(int k=0;k<pdata->nz;k++)
	{
		for(int j=0;j<pdata->ny;j++)
		{
			for(int i=0;i<pdata->nx;i++)
			{
				realkind x = i*pdata->dxdi+pdata->xmin;
				//realkind y = j*pdata->dydi+pdata->ymin;
				//realkind z = k*pdata->dzdi+pdata->zmin;

				realkind Ex = alpha*sin(x/kt)*(kt)*pdata->q_h[0]/(pdata->xi*pdata->xi);



				fields -> getE(i,j,k,0) = Ex;
				fields -> getE(i,j,k,1) = 0;
				fields -> getE(i,j,k,2) = 0;

				fields -> getB(i,j,k,0) = 0;
				fields -> getB(i,j,k,1) = 0;
				fields -> getB(i,j,k,2) = 0;

				fields -> getA(i,j,k,0) = 0;
				fields -> getA(i,j,k,1) = 0;
				fields -> getA(i,j,k,2) = 0;

				printf("fields(%i,%i,%i) = %f at %f\n",i,j,k,fields -> getE(i,j,k,0),x);




			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}



//
//	if(myinfo->myid_mpi == 0)
//	{
//		realkind gamma;
//
//		gamma = pi_const * powf(pdata->omega_pe,3.0f);
//
//		if(pdata->plot_flag)
//		{
//			plot = gnuplot_init();
//			kenergy_plot = gnuplot_init();
//			Tenergy_plot = gnuplot_init();
//			charge_cons_plot = gnuplot_init();
//		}
//		Energy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
//		kEnergy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
//		TEnergy_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
//		charge_cons_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
//
//		max_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
//		max_t_array = (float*)malloc(2*pdata->nsteps*sizeof(float));
//		nplot = 0;
//		nmax_array = 0;
//
//
//		kE0 = moments[0] -> evaluate_energy();
//		kEnergy_array[0] = kE0;
//		E0 = fields[0] -> evaluate_energy();
//		Energy_array[0] = (E0);
//
//		TEnergy_array[0] =0;
//		max_t_array[0] = 0;
//
//		charge_cons_array[0] = 0;
//
//		nplot++;
//
//		if(pdata->plot_flag) gnuplot_cmd(plot,"set log y");
//		//gnuplot_cmd(kenergy_plot,"set log y");
//
//
//	}

}

void TwoStream_Initializer::check_step(NodeParticleList* particles,
		NodeHOMoments* moments,	FieldDataCPU* fields_old, FieldDataCPU* fields_next)
{
//
//		if((pdata->mynode < 2)&&(pdata->plot_flag))
//			particles->plot_particles(pdata);
//
//	if(myinfo->myid_mpi == 0)
//	{
//		if((pdata->mynode < 2)&&(pdata->plot_flag))
//		particles->plot_particles(pdata);
//
//
//		double pEnergy = ((fields[0] -> evaluate_energy()));
//		double kEnergy = moments[0] -> evaluate_energy();
//
//		Energy_array[nplot] = pEnergy;
//		kEnergy_array[nplot] = kEnergy;
//
//		TEnergy_array[nplot] = (E0+kE0 - (kEnergy + pEnergy))/(E0+kE0);
//		//printf("Energy Error= %e\n",(E0+kE0 - (Energy_array[nplot]+kEnergy_array[nplot]))/(E0+kE0));
//		max_t_array[nplot] = (nplot-1)*pdata->dt*2;
//		if(nplot == 1)
//		{
//			max_array[0] = log(Energy_array[0]);
//			max_t_array[nmax_array] = pdata->dt;
//			nmax_array++;
//		}
//		else if(nplot > 1)
//		{
//			if(Energy_array[nplot-1] > fmax(Energy_array[nplot-2],Energy_array[nplot]))
//			{
//				//max_array[nmax_array] = log(Energy_array[nplot-1]);
//				//max_t_array[nmax_array] = (nplot-1);
//
//				//realkind period = (max_t_array[nmax_array]-max_t_array[nmax_array-1])*pdata->dt;
//				//realkind vphase = 2.0*pi_const/period* 1.0*pdata->Lx/(2.0*pi_const);
//				//printf("Energy fluctuation period = %f\n",period);
//				//printf("Phase velocity = %f\n",vphase);
//				//printf("Energy decay rate = %f\n",(max_array[nmax_array]-max_array[nmax_array-1])/
//				//		(pdata->dt*(max_t_array[nmax_array]-max_t_array[nmax_array-1])));
//
//				nmax_array++;
//			}
//		}
//
//		realkind charge_cons = moments[0] -> check_charge(moments_old);
//
//		printf("Charge conservation = %e\n",charge_cons);
//		charge_cons_array[nplot] = charge_cons;
//
//
//		if((nplot > 0)&&(nplot%10 == 9))
//		{
//			if(pdata->plot_flag){
//			gnuplot_resetplot(plot);
//			gnuplot_resetplot(kenergy_plot);
//			gnuplot_resetplot(Tenergy_plot);
//			gnuplot_resetplot(charge_cons_plot);
//			}
//		}
//
//		nplot++;
//
//
//		if((pdata->plot_flag)&&(nplot%10 == 9))
//		{
//			gnuplot_plot_xy(plot,max_t_array+1,Energy_array+1,nplot-1,"Field Energy");
//			gnuplot_cmd(plot,"replot \"E_vec.csv\" title \"Matlab Code\" with points");
//			gnuplot_plot_x(kenergy_plot,kEnergy_array,nplot,"Particle Energy");
//			gnuplot_plot_x(Tenergy_plot,TEnergy_array,nplot,"Total Energy");
//			gnuplot_plot_x(charge_cons_plot,charge_cons_array,nplot,"Charge Conservation");
//		}
//
//
//	}

}

void TwoStream_Initializer::finish(NodeParticleList* particles,
		NodeHOMoments* moments,NodeFieldData* fields_half)
{


//	// Print out the sum of the x-component of the current
//	if(myinfo->myid_mpi == 0)
//	{
//		if(pdata->plot_flag)
//		{
//			gnuplot_resetplot(plot);
//			gnuplot_cmd(plot,"set xlabel \"time\"");
//			gnuplot_cmd(plot,"set ylabel \"Field Energy\"");
//			gnuplot_plot_xy(plot,max_t_array,Energy_array,nplot,"PlasmaApp3D");
//			gnuplot_cmd(plot,"replot \"E_vec.csv\" title \"Matlab Code\" with points");
//			//getchar();
//			gnuplot_save_pdf(plot,"TS_field_energy");
//		}
//
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);

}

