#ifndef ISLAND_INITIALIZER_H
#define ISLAND_INITIALIZER_H

#include "../ProblemInitializer.h"
#include "../PlasmaData.h"

class ParticleList;
class HOMoments;
class FieldData;


class Island_Initializer : public ProblemInitializer
{
	/**@class Island_Initializer Island_Initializer.h Problem_Initializers/Island_Initializer.h
	 * This initializer is for the 2D3V Magnetic Island Coalescence problem.
	 * A detailed problem description can be found in \cite Pritchett1992 \cite Knoll2006
	 *
	 */
public:

	Island_Initializer(PlasmaData* pdata_in)
	{
		pdata = pdata_in;
		kt = pdata->Lx/(2.0*3.14159265359);
		alpha = 0.001;
		title = "Island Coalescence";

		B0 = 0.3;

		lambda = 1.0/sqrt(pdata->mspecies[0]/pdata->mspecies[1]);

		cdf = new island_coalescence(0.25,
				lambda,
				pdata->Lx,
				pdata->Ly,
				0.0,
				pdata->ymin,
				pdata->nx,
				pdata->ny);



	}
	~Island_Initializer(){};

	void initialize_particles(ParticleList* particles, HOMoments* moments,ParallelInfo* myinfo);

	void initialize_fields(FieldData** fields, HOMoments** moments,ParallelInfo* myinfo);

	void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,int ispecies,int iptcl);


	void init_velocities(realkind& vx, realkind& vy, realkind& vz,
						int ispecies,int iptcl);

	void check_step(ParticleList* particles,
			HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo);

	void finish(ParticleList* particles,
			HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo);


	double E0,kE0;
	double B0,lambda;
	float kt;
	float alpha;
	gnuplot_ctrl* plot;
	gnuplot_ctrl* kenergy_plot;
	gnuplot_ctrl* Tenergy_plot;
	float* Energy_array;
	float* kEnergy_array;
	float* TEnergy_array;
	float* max_array;
	float* max_t_array;
	int nplot;
	int nmax_array;
	island_coalescence* cdf;

};

#endif /* ISLAND_INITIALIZER_H */
