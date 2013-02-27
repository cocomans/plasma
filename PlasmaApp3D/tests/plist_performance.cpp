#include "../FieldDataCPU.h"
#include "../ParticleListCPU.h"
#include "../ParticleListCPUAoS.h"
#include "../PlasmaData.h"
#include "../CPUTimer.h"
#include "../HOMoments.h"
#include <omp.h>
#include <gnuplot_i.h>
#include <time.h>

int main(int argc,char* argv[])
{

	PlasmaData pdata(argc,argv);

	pdata.nspecies = 1;

	int nspecies = pdata.nspecies;
	int nptcls = pdata.nptcls;
	int steps = pdata.nsteps;

	int nx = pdata.nx;
	int ny = pdata.ny;
	int nz = pdata.nz;

	pdata.ndimensions = 3;

	int iptcl[nptcls];

	float Ey = 2.0;
	float Bz = 50.0;


	FieldDataCPU fields;
	ParticleListCPU particles;
	ParticleListCPU particles2;
	HOMoments* moments;

	int numprocs = omp_get_num_procs();

	moments = (HOMoments*)malloc(numprocs*sizeof(HOMoments));

	for(int i=0;i<numprocs;i++)
	{
		moments[i] = *new HOMoments(&pdata);
	}


	fields.allocate(&pdata);
	particles.allocate(nptcls);
	particles2.allocate(nptcls);

	fields.dx = pdata.dxdi;
	fields.dy = pdata.dydi;
	fields.dz = pdata.dzdi;

	particles.ispecies = 0;








	for(int i=0;i<nptcls;i++)
	{
		iptcl[i] = i;

		particles.px[i] = rand()%10000/10000.0;
		particles.py[i] = rand()%10000/10000.0;
		particles.pz[i] = 0.5;

		particles.ix[i] = i%nx;
		particles.iy[i] = (i%(nx*ny))/nx;
		particles.iz[i] = 0;

		particles.vx[i] = 0.5*(2*(rand()%10000))/10000.0;
		particles.vy[i] = 0.5*(2*(rand()%10000))/10000.0;
		particles.vz[i] = 0.1* (rand()%50000 / 50000.0f - 0.5);

		particles.dt_finished[i] = 0;


	}

	particles2.copy_from(&particles);


	for(int i=0;i<nptcls;i++)
	{
		float px1,py1,pz1,vx1,vy1,vz1;
		int ix1,iy1,iz1;

		float px2,py2,pz2,vx2,vy2,vz2;
		int ix2,iy2,iz2;

		px1 = particles.get_fvalue(i,0);
		py1 = particles.get_fvalue(i,1);
		pz1 = particles.get_fvalue(i,2);

		vx1 = particles.get_fvalue(i,3);
		vy1 = particles.get_fvalue(i,4);
		vz1 = particles.get_fvalue(i,5);

		ix1 = particles.get_ivalue(i,0);
		iy1 = particles.get_ivalue(i,1);
		iz1 = particles.get_ivalue(i,2);



		px2 = particles2.get_fvalue(i,0);
		py2 = particles2.get_fvalue(i,1);
		pz2 = particles2.get_fvalue(i,2);

		vx2 = particles2.get_fvalue(i,3);
		vy2 = particles2.get_fvalue(i,4);
		vz2 = particles2.get_fvalue(i,5);

		ix2 = particles2.get_ivalue(i,0);
		iy2 = particles2.get_ivalue(i,1);
		iz2 = particles2.get_ivalue(i,2);


		if((px1 != px2)||(py1 != py2)||(pz1 != pz2)||
		   (vx1 != vx2)||(vy1 != vy2)||(vz1 != vz2))
		{
			printf("Data Structure Values Different for particle %i!!!\n",i);
			printf("%f, %f, %f, %f, %f, %f\n",px1,py1,pz1,vx1,vy1,vz1);
			printf("%f, %f, %f, %f, %f, %f\n",px2,py2,pz2,vx2,vy2,vz2);
		}



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

				float Ex = -1.0*x;


				fields.getE(i,j,k,0) = 0;
				fields.getE(i,j,k,1) = Ey;
				fields.getE(i,j,k,2) = 0.0;

				fields.getB(i,j,k,0) = 0;
				fields.getB(i,j,k,1) = 0;
				fields.getB(i,j,k,2) = Bz;


			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

	fields.q2m[0] = 1.0;

	printf("Efield setup complete\n");

	float time;
	double avg_error = 0.0;
	int n_error = 0;

	CPUTimer timer;

	CPUTimer timer2;


	moments->init_plot();

	int nsteps1 = 0;
	int nsteps2 = 0;

	timer.start();
	for(int i=0;i<steps;i++)
	{
		//time = dtau0*(i);

		//moments.set_vals(0);
		nsteps1 += particles.push(&pdata,&fields,moments);
		printf("finished step %i\n",i);


	}
	timer.stop();


	timer2.start();
	for(int i=0;i<steps;i++)
	{
		//time = dtau0*(i);

		//moments.set_vals(0);
		nsteps2 += particles2.push(&pdata,&fields,moments);
		printf("finished step %i\n",i);


	}
	timer2.stop();

	printf("SoA Run did %f particle-steps per second with %i steps\n",nsteps1/(timer.diff()*1.0e-3),nsteps1);
	printf("AoS Run did %f particle-steps per second with %i steps\n",nsteps2/(timer2.diff()*1.0e-3),nsteps2);

	//moments->plot(nz/2,0,HOMoments_currentx);


	printf("Press 'Enter' to continue\n");
		getchar();

	//moments->close_plot();


	//particles2.CPUfree();

}
