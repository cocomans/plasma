#ifndef IONACOUSTIC_INITIALIZER_H
#define IONACOUSTIC_INITIALIZER_H

#include "../ProblemInitializer.h"

class PlasmaData;
class ParticleList;
class HOMoments;
class FieldData;


class IonAcoustic_Initializer : public ProblemInitializer
{
public:

	IonAcoustic_Initializer(PlasmaData* pdata_in);

	~IonAcoustic_Initializer(){};

	void initialize_particles(NodeParticleList* particles, NodeHOMoments* moments);

	void initialize_fields(FieldDataCPU* fields, NodeHOMoments* moments);

	void init_particle(realkind& px, realkind& py, realkind& pz, int& ix, int& iy, int& iz,
								realkind& vx,realkind& vy,realkind& vz,
								int ispecies,int iptcl, int nptcls,int ioffset);


	void init_velocities(realkind& vx, realkind& vy, realkind& vz,
						int ispecies,int iptcl);






	float Ey,Bz;
	double E0,kE0;
	float alpha;
	float kt;
	float v_shift;


};

#endif /* IONACOUSTIC_INITIALIZER_H */
