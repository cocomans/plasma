#include "../FieldDataCPU.h"
#include "../ParticleListCPU.h"
#include "../PlasmaData.h"
#include "../CPUTimer.h"
#include "../HOMoments.h"
#include "../AmpereSolve.h"
#include "../ShapeFunctions.h"
#include <omp.h>
#include <gnuplot_i.h>
#include <time.h>


int main(int argc, char* argv[])
{


	PlasmaData pdata(argc,argv);

	pdata.nspecies = 2;

	int nspecies = pdata.nspecies;
	int nptcls = pdata.nptcls;
	int steps = pdata.nsteps;

	int nx = pdata.nx;
	int ny = pdata.ny;
	int nz = pdata.nz;

	int iptcl[nptcls];

	float Ey = 2.0;
	float Bz = 100.0;

	double energy_total = 0.0;

	pdata.niter_max = 200;

	pdata.nSubcycle_max = 100;

	pdata.setup();

	FieldDataCPU fields;
	ParticleListCPU electrons;
	ParticleListCPU electrons_next;
	ParticleListCPU ions;
	HOMoments* moments;

	int numprocs = omp_get_num_procs();

	moments = (HOMoments*)malloc(numprocs*sizeof(HOMoments));

	for(int i=0;i<numprocs;i++)
	{
		moments[i] = *new HOMoments(&pdata);
	}

	//float x_plot[nptcls][steps];
	//float y_plot[nptcls][steps];



	fields.allocate(nx,ny,nz,nspecies);


	fields.dx = pdata.dxdi;
	fields.dy = pdata.dydi;
	fields.dz = pdata.dzdi;

	electrons.allocate(nptcls);
	electrons.ispecies = 0;
	electrons_next.allocate(nptcls);
	electrons_next.ispecies = 0;

	ions.allocate(nptcls);
	ions.ispecies = 1;

	pdata.qspecies[1] = 1.0;
	pdata.mspecies[1] = 1.0;






	for(int i=0;i<nptcls;i++)
	{
		iptcl[i] = i;

		electrons.px[i] = rand()%100/100.0;
		electrons.py[i] = rand()%100/100.0;
		electrons.pz[i] = 0.5;

		electrons.ix[i] = (nx/2);
		electrons.iy[i] = rand()%(ny/2)+ny/4;
		electrons.iz[i] = nz/2;

		electrons.vx[i] = 1.5*((2*(rand()%10000))/10000.0 - 1);
		electrons.vy[i] = 1.5*((2*(rand()%10000))/10000.0 - 1);
		electrons.vz[i] = 0.00*(rand()%50 / 50.0f - 0.5);

		energy_total += pdata.n0*0.5*qe/(qe2me*pdata.nptcls)*(powf(electrons.vx[i],2)+powf(electrons.vy[i],2)+powf(electrons.vz[i],2));
	}




	// Setup E-field
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			for(int k=0;k<nz;k++)
			{
				float x = i*pdata.dxdi+pdata.xmin;
				float y = j*pdata.dydi+pdata.ymin;
				float z = k*pdata.dzdi+pdata.zmin;



				float Ex = sgn(x)*pdata.n0*qe*qe2me*(1.0-(x));
				Ey =  0;


				fields.getE(i,j,k,0) = 0;
				fields.getE(i,j,k,1) = 0;
				fields.getE(i,j,k,2) = 0.0;

				fields.getB(i,j,k,0) = 0;
				fields.getB(i,j,k,1) = 0;
				fields.getB(i,j,k,2) = 0;


			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

	fields.q2m[0] = 1.0;
	fields.q2m[1] = -1.0/1.0;

	printf("Efield setup complete\n");

	float time;
	double avg_error = 0.0;
	int n_error = 1;

	CPUTimer timer;


	moments->init_plot();
	fields.init_plot();


	FieldDataCPU fields_old = fields;
	FieldDataCPU fields_next = fields;
	FieldDataCPU fields_half = fields;


	fields_next.allocate(nx,ny,nz,nspecies);
	fields_half.allocate(nx,ny,nz,nspecies);

	timer.start();
	for(int i=0;i<steps;i++)
	{
		//time = dtau0*(i);

		fields_next.copy_from(&fields);


		double residual = 1.0e-6 * 2.0;

		while(residual > 1.0e-6)
		{
			electrons_next.copy_from(&electrons);

			Average_Fields(&fields,&fields_next,&fields_half);

			for(int j=0;j<numprocs;j++)
				moments[j].set_vals(0);


			electrons_next.push(&pdata,&fields_half,moments);
			moments->reset_plot();
			moments->plot(nz/2,0,HOMoments_currentx);


			residual = Calc_Residual(&pdata,&fields_next,&fields,moments);

			printf("Field Residual = %e\n",residual);


			// Do the field solve
			AmpereSolve(&pdata,&fields_next,&fields_old,moments);





		}

		electrons.copy_from(&electrons_next);

		fields.copy_from(&fields_next);

		fields.reset_plot();
		fields.plot(&pdata,nz/2,0,0,0);





		printf("finished step %i\n",i);
		double kenergy_total_step = 0.0;
		double penergy_step;

		for(int j=0;j<nptcls;j++)
		{

			kenergy_total_step += pdata.n0*0.5*qe/(qe2me*pdata.nptcls)
					*(powf(electrons.vx[j],2)
					+powf(electrons.vy[j],2)
					+powf(electrons.vz[j],2));
		}

		penergy_step = fields.evaluate_energy();

		float energy_step = penergy_step + kenergy_total_step;

		printf("Total Energy = %e, -> %e Kinetic + %e Potential = %e\n",
				energy_total,kenergy_total_step,penergy_step,energy_step);
		printf("Energy Error = %e\n",
				fabs(energy_total-energy_step)/energy_total);







	}
	timer.stop();


	printf("average error = %e \n",avg_error/((float)n_error));
	printf("Run did %f electrons per second\n",nptcls*steps/(timer.diff()*1.0e-3));


	moments->plot(nz/2,0,HOMoments_currentx);


	printf("Press 'Enter' to continue\n");
		getchar();

	moments->close_plot();
	fields.close_plot();


}
