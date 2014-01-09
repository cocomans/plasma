//===============================================================================
//
//	AUTHOR: WILLIAM T. TAITANO
//	This class will store the LO moment quantities
//
//===============================================================================

#ifndef LOMOMENTSES_H_
#define LOMOMENTSES_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <gnuplot_i.h>
#include <mpi.h>
#include "../PlasmaData.h"
#include "../HOMoments.h"

enum LOMoments_moment
{
	LOMoments_charge = 0,
	LOMoments_currentx = 1,
	LOMoments_currenty = 2,
	LOMoments_currentz = 3,
};

class LOMomentsES {
public:
	//================================================
	//	Constructor
	//================================================
	LOMomentsES(PlasmaData* pDataIn, HOMoments* hoMoments);
	//================================================
	//	Destructor
	//================================================
	~LOMomentsES();
	//================================================
	//	MemberFunction: get_val
	//================================================
	realkind& get_val(	const int ix, const int iy, const int iz,
						const int ispecies,enum LOMoments_moment moment);
	//================================================
	//	Members
	//================================================
	/// charge (0th Moment)
	realkind* charge;
	/// x component of current (1st Moment)
	realkind* currentx;
	/// y component of current (1st Moment)
	realkind* currenty;
	/// z component of current (1st Moment)
	realkind* currentz;
	/// xx component of Stress (2nd Moment)
	realkind* S2;
	/// Simulation information
	PlasmaData* pData;
};


#endif /* LOMOMENTSES_H_ */
