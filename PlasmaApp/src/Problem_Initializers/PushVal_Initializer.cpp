#include "PushVal_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../NodeHOMoments.h"
#include "../NodeParticleList.h"
#include "../FieldDataCPU.h"
#include "../PlasmaData.h"
#include "../RunData.h"
#include "../ShapeFunctions.h"

#include "../ParallelInfo.h"

using namespace std;

PushVal_Initializer::PushVal_Initializer(PlasmaData* pdata_in)
{
		pdata = pdata_in;

		const char* val_file = "validation/particle_pos";



		string line;
		ifstream myfile(val_file);

		if(myfile.is_open())
		{

			getline(myfile,line);
			getline(myfile,line,',');
			nt_v = atoi(line.c_str());
			getline(myfile,line,',');
			np_v = atoi(line.c_str());
			getline(myfile,line,',');
			dt_ref = atof(line.c_str());
			getline(myfile,line,',');
			printf("np = %i, nt = %i, dt_ref = %e\n",np_v,nt_v,dt_ref);

			val_particles = (ValParticle**)malloc(nt_v*sizeof(ValParticle*));

			val_particles[0] = (ValParticle*)malloc(np_v*nt_v*sizeof(ValParticle));

			time = (double*)malloc(nt_v*sizeof(double));

			for(int i=1;i<nt_v;i++)
			{
				val_particles[i] = val_particles[0] + np_v*i;
			}


			int counter = 0;

			// skip the first line, don't need it

			getline(myfile,line);

			bool igood = 1;
			while(igood)
			{
				string value;
				igood = myfile.good();
				int icol = counter%7;
				int irow = counter/7;
				int itime = irow/np_v;
				int iptcl = irow%np_v;

				if(itime > nt_v-1)
					break;
				double tempx;
				getline(myfile,value,',');
//printf("%i, %i: %s\n",itime,iptcl,value.c_str());
				switch(icol)
				{
				case 0:

					break;
				case 1:

					//time[itime] = (double)atof(value.c_str());
					break;
				case 2:
					 tempx = (double)atof(value.c_str());;
					val_particles[itime][iptcl].x = tempx;
					break;
				case 3:
					val_particles[itime][iptcl].y = (double)atof(value.c_str());
					break;
				case 4:
					val_particles[itime][iptcl].vx = (double)atof(value.c_str());
					break;
				case 5:
					val_particles[itime][iptcl].vy = (double)atof(value.c_str());
					break;
				case 6:
					val_particles[itime][iptcl].vz = (double)atof(value.c_str());
					break;
				default:
					break;
				}




				counter++;
			}

			printf("finished reading in particles\n");
			myfile.close();
		}









}

void PushVal_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
		 int& ix, int& iy, int& iz,
		 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl,
		 int nptcls,int ioffset)
{

	int iptcl_new = iptcl%np_v;

	px = val_particles[0][iptcl_new].x - pdata->xmin;
	py = val_particles[0][iptcl_new].y - pdata->ymin;


	px *= pdata->didx;
	py *= pdata->didy;

	ix = floor(px);
	iy = floor(py);
	px = px-ix;
	py = py-iy;

	vx = val_particles[0][iptcl_new].vx;
	vy = val_particles[0][iptcl_new].vy;
	vz = val_particles[0][iptcl_new].vz;



}

void PushVal_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{




}

void PushVal_Initializer::initialize_particles(NodeParticleList* particles,
										  NodeHOMoments* moments)
{


	particles -> init(this,moments);

	int nptcls = particles->cpu_particles->nptcls;

	ParticleListCPU* temp_list = particles->cpu_particles;
	for(int i=0;i<nptcls;i++)
	{
		double px = temp_list->get_fvalue(i,0);
		double py = temp_list->get_fvalue(i,1);

		int ix = temp_list->get_ivalue(i,0);
		int iy = temp_list->get_ivalue(i,1);

		double vx = temp_list->get_fvalue(i,3);
		double vy = temp_list->get_fvalue(i,4);
		double vz = temp_list->get_fvalue(i,5);

		int iptcl_ref = i%np_v;

		double x = (px+ix)*pdata->dxdi + pdata->xmin;
		double y = (py+iy)*pdata->dydi + pdata->ymin;

		double ts = pdata->dt*pdata->rstatus->istep;
		int is = floor(ts/0.01);
		ts = ts - is;

		ValParticle rp0 = val_particles[min(is,np_v-1)][iptcl_ref];
		ValParticle rp1 = val_particles[min(is+1,np_v-1)][iptcl_ref];

		ValParticle rp;

		double sm = S1_shape(ts);
		double sp = S1_shape(1.0-ts);

		rp.x = sm*rp0.x+sp*rp1.x;
		rp.y = sm*rp0.y+sp*rp1.y;

		rp.vx = sm*rp0.vx+sp*rp1.vx;
		rp.vy = sm*rp0.vy+sp*rp1.vy;
		rp.vz = sm*rp0.vz+sp*rp1.vz;


		double x_error = fabs(rp.x - x)/(fabs(x)+fabs(rp.x));
		double y_error = fabs(rp.y - y)/(fabs(y)+fabs(rp.y));

		double vx_error = fabs(rp.vx - vx)/(fabs(vx)+fabs(rp.vx));
		double vy_error = fabs(rp.vy - vy)/(fabs(vy)+fabs(rp.vy));
		double vz_error = fabs(rp.vz - vz)/(fabs(vz)+fabs(rp.vz));

		double t_error = sqrt(x_error*x_error
							+y_error*y_error
							+vx_error*vx_error
							+vy_error*vy_error
							+vz_error*vz_error);


			printf("Particle[%i] = %e %e %e %e %e\n",i,
					x,y,vx,vy,vz);


	}




}

void PushVal_Initializer::initialize_fields(FieldDataCPU* fields,
										NodeHOMoments* moments)
{
	// Setup E-field
	for(int k=0;k<pdata->nz;k++)
	{
#pragma omp parallel for
		for(int j=0;j<pdata->ny;j++)
		{
			for(int i=0;i<pdata->nx;i++)
			{
				realkind x = (i+0.5)*pdata->dxdi+pdata->xmin;
				realkind y = (j+0.5)*pdata->dydi+pdata->ymin;
				realkind x2 = (i)*pdata->dxdi+pdata->xmin;
				realkind y2 = (j)*pdata->dydi+pdata->ymin;

				//realkind z = k*pdata->dzdi+pdata->zmin;

				fields -> getE(i,j,k,0) = y*y;
				fields -> getE(i,j,k,1) = x*x;
				fields -> getE(i,j,k,2) = 0.0;

				fields -> getB(i,j,k,0) = cos(x2)+sin(y2);
				fields -> getB(i,j,k,1) = sin(x2)+cos(y2);
				fields -> getB(i,j,k,2) = 10.0;

//				fields -> getE(i,j,k,0) = cos(x+y);
//				fields -> getE(i,j,k,1) = sin(x+y);
//				fields -> getE(i,j,k,2) = 1.0;
//
//				fields -> getB(i,j,k,0) = -sin(x+y);
//				fields -> getB(i,j,k,1) = -cos(x+y);
//				fields -> getB(i,j,k,2) = -10.0*(sin(x+y)+cos(x+y));


			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

}

void PushVal_Initializer::check_step(NodeParticleList* particles_next,
		NodeHOMoments* moments,FieldDataCPU* fields_old,FieldDataCPU* fields_next)
{

	int nptcls = particles_next->cpu_particles->nptcls;

	ParticleListCPU* temp_list = particles_next->cpu_particles;

	double norm_error=0;
	for(int i=0;i<nptcls;i++)
	{
		double px = temp_list->get_fvalue(i,0);
		double py = temp_list->get_fvalue(i,1);

		int ix = temp_list->get_ivalue(i,0);
		int iy = temp_list->get_ivalue(i,1);

		double vx = temp_list->get_fvalue(i,3);
		double vy = temp_list->get_fvalue(i,4);
		double vz = temp_list->get_fvalue(i,5);

		int iptcl_ref = i%np_v;

		double x = (px+ix)*pdata->dxdi + pdata->xmin;
		double y = (py+iy)*pdata->dydi + pdata->ymin;

		double ts0 = pdata->dt*(pdata->rstatus->istep);
		int is = floor(ts0/dt_ref);
		double ts = ts0/dt_ref - is;

		ValParticle rp0 = val_particles[min(is,nt_v-1)][iptcl_ref];
		ValParticle rp1 = val_particles[min(is+1,nt_v-1)][iptcl_ref];

		ValParticle rp;

		double sm = S1_shape(ts);
		double sp = S1_shape(1.0-ts);

		rp.x = sm*rp0.x+sp*rp1.x;
		rp.y = sm*rp0.y+sp*rp1.y;

		rp.vx = sm*rp0.vx+sp*rp1.vx;
		rp.vy = sm*rp0.vy+sp*rp1.vy;
		rp.vz = sm*rp0.vz+sp*rp1.vz;


		double x_error = fabs(rp.x - x)/(fabs(x)+fabs(rp.x));
		double y_error = fabs(rp.y - y)/(fabs(y)+fabs(rp.y));

		double vx_error = fabs(rp.vx - vx)/(fabs(vx)+fabs(rp.vx));
		double vy_error = fabs(rp.vy - vy)/(fabs(vy)+fabs(rp.vy));
		double vz_error = fabs(rp.vz - vz)/(fabs(vz)+fabs(rp.vz));

		double t_error = sqrt(x_error*x_error
							+y_error*y_error
							+vx_error*vx_error
							+vy_error*vy_error
							+vz_error*vz_error)/5.0;

		norm_error += t_error*t_error;

//		printf("error = %e\n",t_error);

		if(t_error > 1.0e-4)
		{
//			printf("Warning Particle position error[%i] at %f(%i) = %e %e %e %e %e\n",i,ts0,is,
//					x_error,y_error,vx_error,vy_error,vz_error);

			printf("Warning Particle position error[%i] at %f(%i) = %e(%e) %e(%e) %e(%e) %e(%e) %e(%e)\n",
					i,ts0,is,
					x,rp.x,y,rp.y,vx,rp.vx,vy,rp.vy,vz,rp.vz);
		}

	}

	printf("norm error = %e\n",sqrt(norm_error)/nptcls);

	//getchar();

}

void PushVal_Initializer::finish(NodeParticleList* particles,
		NodeHOMoments* moments,NodeFieldData* fields_half)
{




}

