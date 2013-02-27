#ifndef TWOSTREAM_INITIALIZER_H
#define TWOSTREAM_INITIALIZER_H

#include "../ProblemInitializer.h"
#include "../PlasmaData.h"

class ParticleList;
class HOMoments;
class FieldData;


class TwoStream_Initializer : public ProblemInitializer
{
public:

	TwoStream_Initializer(PlasmaData* pdata_in)
	{
		pdata = pdata_in;
		kt = pdata->Lx/(2.0*3.14159265359);
		alpha = 0.001;
		title = "TwoStream";



	}
	~TwoStream_Initializer(){};

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
	float kt;
	float alpha;
	gnuplot_ctrl* plot;
	gnuplot_ctrl* kenergy_plot;
	gnuplot_ctrl* Tenergy_plot;
	gnuplot_ctrl* charge_cons_plot;
	float* Energy_array;
	float* kEnergy_array;
	float* TEnergy_array;
	float* max_array;
	float* max_t_array;
	float* charge_cons_array;
	int nplot;
	int nmax_array;

};

#endif /* TWOSTREAM_INITIALIZER_H */
