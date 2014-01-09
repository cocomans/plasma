#include "Weibel_Initializer.h"
#include "../PlasmaData.h"
#include "../ParticleListCPU.h"
#include "../ParticleListGPU.cuh"
#include "../NodeHOMoments.h"
#include "../NodeParticleList.h"
#include "../FieldDataCPU.h"
#include "../PlasmaData.h"

#include "../ParallelInfo.h"
#include "../RunData.h"
#include"../Util/rand.h"
#include "../Field_Solvers/InitSolve.h"

Weibel_Initializer::Weibel_Initializer(PlasmaData* pdata_in)
{
	printf("Constructing Initializer\n");
	pdata = pdata_in;
	kt = pdata->Lx/(2.0*3.14159265359);
	alpha = 0.01;
	title = "Weibel";

	pdata->rdata->SimName = "Weibel";

	np_in = 256000;

	px_dist[0] = (double*)malloc(np_in*sizeof(double));
	ix_dist[0] = (int*)malloc(np_in*sizeof(int));

	px_dist[1] = (double*)malloc(np_in*sizeof(double));
	ix_dist[1] = (int*)malloc(np_in*sizeof(int));



	vx_dist[0] = (double*)malloc(np_in*sizeof(double));
	vy_dist[0] = (double*)malloc(np_in*sizeof(double));
	vz_dist[0] = (double*)malloc(np_in*sizeof(double));

	vx_dist[1] = (double*)malloc(np_in*sizeof(double));
	vy_dist[1] = (double*)malloc(np_in*sizeof(double));
	vz_dist[1] = (double*)malloc(np_in*sizeof(double));

	random_d = (double*)malloc(pdata->nptcls*sizeof(double));

	dfdn = 2.0/pdata->nptcls;

	std::ofstream ofile("vxp_check.txt");

	const char* vxp_file = "vxp/vxp.txt";
	std::ifstream* input_file = new std::ifstream(vxp_file);
	std::string line;
	std::string delim = " ";
	 if ((input_file->is_open())&&(ofile.is_open()))
	  {
	    while ( std::getline (*input_file,line) )
	    {
			 std::stringstream sstream(line);

	    	if(line.find("#") == std::string::npos)
	    	{
	    		int i = 0;
	    		int ispecies;
	    		int iptcl;
	    		std::string value;

	    		char* delim = " ";
	    		double tmp;

	    		while(sstream >> tmp)
	    		{

	    			switch(i)
	    			{
	    			case 0:
	    				ispecies = ((int)tmp) - 1;
	    				break;
	    			case 1:
	    				iptcl = ((int)tmp) - 1;
	    				break;
	    			case 2:
	    				ix_dist[ispecies][iptcl] = ((int)tmp) - 1;
	    				break;
	    			case 3:
	    				px_dist[ispecies][iptcl] = pdata->didx*tmp;
	    				break;
	    			case 4:
	    				vx_dist[ispecies][iptcl] = tmp;
	    				break;
	    			case 5:
	    				vy_dist[ispecies][iptcl] = tmp;
	    				break;
	    			case 6:
	    				vz_dist[ispecies][iptcl] = tmp;
	    				break;
	    			}

	    			i++;
	    		}



	    	}
	    }

	    input_file->close();
	  }


	for(int l=0;l<pdata->nspecies;l++)
		for(int i=0;i<pdata->nptcls;i++)
		{


//			init_velocities(vx_dist[l][i],vy_dist[l][i],vz_dist[l][i],l,i);
//
//			int ivel = i/(pdata->nptcls/2);
//			int itemp = i%(np_in/2);
//			if(ivel == 1)
//			{
//				vx_dist[l][i] = -vx_dist[l][itemp];
//				vy_dist[l][i] = -vy_dist[l][itemp];
//				vz_dist[l][i] = -vz_dist[l][itemp];
//				px_dist[l][i] =  px_dist[l][itemp];
//				ix_dist[l][i] =  ix_dist[l][itemp];
//
//			}

//			random_d[i] = drandom();
//    		char temp[1024];
//    		sprintf(temp,"Reading in Particle: %i %i %i %e %e %e %e\n",l,i,
//    				ix_dist[l][i],
//    				px_dist[l][i],
//    				vx_dist[l][i],
//    				vy_dist[l][i],
//    				vz_dist[l][i]);
//
//    		ofile << temp;

		}

    ofile.close();


	printf("Finished Constructing Initializer\n");

}



void Weibel_Initializer::init_particle(realkind& px, realkind& py, realkind& pz,
		 int& ix, int& iy, int& iz,
		 realkind& vx, realkind& vy, realkind& vz,int ispecies, int iptcl,
		 int nptcls,int ioffset)
{
	// Set Position Values, ifloat = 0-2

//	raised_cos cdf(alpha,1.0/kt);
//
//	int iptcl_new = iptcl;
//
//	int idist = iptcl_new%(pdata->nptcls_device[0]*pdata->node_info->nTasks_g/2);
//	int ivel = iptcl_new/(pdata->nptcls_device[0]*pdata->node_info->nTasks_g/2);
//
////	idist += ioffset;
//	idist = iptcl_new%(nptcls/4);
//	ivel = iptcl_new/(nptcls/2);
//	double temp_x;
//
//
////	temp_x = distribution_intrp((idist)/(0.5*nptcls),0.0,pdata->Lx,cdf);
//
//	temp_x = (idist)/(nptcls/4.0)*pdata->Lx;
//	temp_x += alpha*cos(temp_x/kt)*pdata->Lx;
////	if(ispecies == 1)
////		temp_x += log(drandom())*pdata->Lx;
//
//	double tempf = (iptcl%(nptcls/2))/(nptcls/2.0);
//	int im = floor(tempf/dfdn);
//	tempf = tempf/dfdn - im;
//	int ip = im + 1;
//
//	im = ((im%(pdata->nptcls/2) + (pdata->nptcls/2))%(pdata->nptcls/2));
//	ip = ((ip%(pdata->nptcls/2) + (pdata->nptcls/2))%(pdata->nptcls/2));
//
//	if(ispecies == 1)
//	{
////		double tempr = random_d[im]*(1.0-tempf)+random_d[ip]*(tempf);
////		temp_x += log(tempr)*pdata->Lx;
//	}
//
//
//	vx = (1-(2*(ivel)))*(vx_dist[ispecies][im]*(1.0-tempf)+vx_dist[ispecies][ip]*(tempf));
//	vy = (1-(2*(ivel)))*(vy_dist[ispecies][im]*(1.0-tempf)+vy_dist[ispecies][ip]*(tempf));
//	vz = (1-(2*(ivel)))*(vz_dist[ispecies][im]*(1.0-tempf)+vz_dist[ispecies][ip]*(tempf));
//
////	if(ispecies == 1)
////	{
////		vx = 0;
////		vy = 0;
////		vz = 0;
////	}
//
//	iptcl = iptcl;
	py = 0.0;
	pz = 0.0;

	iy = floor(py*pdata->didy);
	iz = floor(pz*pdata->didz);

	py = py*pdata->didy - iy;
	pz = pz*pdata->didz - iz;
//
//	//px = drandom();
//

	int ip_t = iptcl%np_in;
	int ivel = iptcl/(nptcls/2);
	ix = ix_dist[ispecies][ip_t];
	px = px_dist[ispecies][ip_t];//+drandom();
	vx = vx_dist[ispecies][ip_t];
	vy = vy_dist[ispecies][ip_t];
	vz = vz_dist[ispecies][ip_t];


	ix = ix + floor(px);

	px = px - floor(px);

//	if(ispecies == 1)
//	{
//		vy /= 4.0;
//		vz /= 4.0;
//	}
//	ix = ((ix%pdata->nx)+pdata->nx)%pdata->nx;


//	iy = 0;
//	iz = 0;



}

void Weibel_Initializer::init_velocities(realkind& vx, realkind& vy, realkind& vz, int ispecies, int iptcl)
{
	double aniso_mult;
//	double v0  = sqrt((2*0.08*0.08+0.02*0.02)/3.0);
	double v0 = 0.1;
	if(ispecies == 0)
	{	aniso_mult = 4.0; v0 = 0.1;}

	else
		aniso_mult = 1.0;


//	vx = box_muller(0.0,v0);
//	vy = box_muller(0.0,(aniso_mult*v0));
//	vz = box_muller(0.0,(aniso_mult*v0));

	randn(vx);
	randn(vy);
	randn(vz);

	vx *= v0;
	vy *= aniso_mult*v0;
	vz *= aniso_mult*v0;

//	if(ispecies != 0)
//	{
//		vx = 0;
//		vy = 0;
//		vz = 0;
//	}



}

void Weibel_Initializer::initialize_particles(NodeParticleList* particles,
										  NodeHOMoments* moments)
{


	printf("Initializing particles\n");

	particles -> init(this,moments);
	MPI_Barrier(MPI_COMM_WORLD);

	printf("Finished init parts\n");

}

void Weibel_Initializer::initialize_fields(FieldDataCPU* fields,
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



				fields -> getE(i,j,k,0) = 0.0;
				fields -> getE(i,j,k,1) = 0;
				fields -> getE(i,j,k,2) = 0;

				fields -> getB(i,j,k,0) = 0.0;
				fields -> getB(i,j,k,1) = 0;
				fields -> getB(i,j,k,2) = 0;

				fields -> getA(i,j,k,0) = 000;
				fields -> getA(i,j,k,1) = 0;
				fields -> getA(i,j,k,2) = 0;

				printf("fields(%i,%i,%i) = %f at %f\n",i,j,k,fields -> getE(i,j,k,0),x);




			//	printf("fields(%i,%i,%i) = %f, %f, %f\n",i,j,k,
				//	fields.getE(i,j,k,0),fields.getE(i,j,k,1),fields.getE(i,j,k,2));
			}
		}
	}

	InitSolve(pdata,fields,moments->moments_next);





}

void Weibel_Initializer::check_step(NodeParticleList* particles,
		NodeHOMoments* moments,	FieldDataCPU* fields_old, FieldDataCPU* fields_next)
{


}

void Weibel_Initializer::finish(NodeParticleList* particles,
		NodeHOMoments* moments,NodeFieldData* fields_half)
{


}

