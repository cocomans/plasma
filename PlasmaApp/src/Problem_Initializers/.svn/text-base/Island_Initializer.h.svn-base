#ifndef ISLAND_INITIALIZER_H
#define ISLAND_INITIALIZER_H

#include "../ProblemInitializer.h"
#include "../PlasmaData.h"




class Island_Initializer : public ProblemInitializer
{
	/**@class Island_Initializer Island_Initializer.h Problem_Initializers/Island_Initializer.h
	 * This initializer is for the 2D3V Magnetic Island Coalescence problem.
	 * A detailed problem description can be found in \cite Pritchett1992 \cite Knoll2006
	 *
	 */
public:

	Island_Initializer(PlasmaData* pdata_in);

	~Island_Initializer(){};

	void initialize_particles(NodeParticleList* particles, NodeHOMoments* moments);

	void initialize_fields(FieldDataCPU* fields, NodeHOMoments* moments);

	void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,int ispecies,int iptcl,
								int nptcls,int ioffset);


	void init_velocities(realkind& vx, realkind& vy, realkind& vz,
						int ispecies,int iptcl);

	void check_step(NodeParticleList* particles_next,
			NodeHOMoments* moments,FieldDataCPU* fields_old,FieldDataCPU* fields_next);

	void finish(NodeParticleList* particles,
			NodeHOMoments* moments,NodeFieldData* fields_half);



	double B0,lambda;
	float kt;
	float alpha;

	island_coalescence* cdf;

};

#endif /* ISLAND_INITIALIZER_H */
