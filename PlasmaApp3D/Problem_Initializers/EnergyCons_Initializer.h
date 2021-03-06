#ifndef ENERGYCONS_INITIALIZER_H
#define ENERGYCONS_INITIALIZER_H

#include "../ProblemInitializer.h"


class ParticleList;
class HOMoments;
class FieldData;

class EnergyCons_Initializer : public ProblemInitializer
{
public:

	EnergyCons_Initializer(PlasmaData* pdata_in);

	~EnergyCons_Initializer(){};

	void initialize_particles(ParticleList* particles, HOMoments* moments,ParallelInfo* myinfo);

	void initialize_fields(FieldData** fields, HOMoments** moments,ParallelInfo* myinfo);

	void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,int ispecies,int iptcl);


	void init_velocities(realkind& vx, realkind& vy, realkind& vz,
						int ispecies,int iptcl);

	void check_step(ParticleList* particles,
			HOMoments** moments,HOMoments* moments_old,FieldData** fields,ParallelInfo* myinfo);


	float Ey,Bz;
	float E0,kE0;
	float kt;
	float v_phase;
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

};

#endif /* ENERGYCONS_INITIALIZER_H */
