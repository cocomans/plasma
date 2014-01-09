#ifndef WEIBEL_INITIALIZER_H
#define WEIBEL_INITIALIZER_H

#include "../ProblemInitializer.h"
#include "../PlasmaData.h"



class Weibel_Initializer : public ProblemInitializer
{
public:

	Weibel_Initializer(PlasmaData* pdata_in);

	~Weibel_Initializer(){};

	void initialize_particles(NodeParticleList* particles, NodeHOMoments* moments);

	void initialize_fields(FieldDataCPU* fields, NodeHOMoments* moments);

	void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,int ispecies,int iptcl,
								int nptcls, int ioffset);


	void init_velocities(realkind& vx, realkind& vy, realkind& vz,
						int ispecies,int iptcl);

	void check_step(NodeParticleList* particles_next,
			NodeHOMoments* moments,FieldDataCPU* fields_old,FieldDataCPU* fields_next);

	void finish(NodeParticleList* particles,
			NodeHOMoments* moments,NodeFieldData* fields_half);

	double E0,kE0;
	float kt;
	float alpha;

	int np_in;

	double* px_dist[2];
	int* ix_dist[2];
	double* vx_dist[2];
	double* vy_dist[2];
	double* vz_dist[2];
	double* random_d;

	double dfdn;

	int nplot;
	int nmax_array;

};

#endif /* WEIBEL_INITIALIZER_H */
