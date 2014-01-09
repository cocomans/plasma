/*
 * NodeHOMoments.h
 *
 *  Created on: Apr 22, 2013
 *      Author: payne
 */

#ifndef NODEHOMOMENTS_H_
#define NODEHOMOMENTS_H_

#include "PlasmaData.h"
#include "HOMomentsCPU.h"
#include "HOMomentsGPU.h"
#include "HOMomentsMIC.h"

class NodeFieldData;
class FieldDataCPU;


class NodeHOMoments
{
public:

	/// root HO Moments for time \f$\mathcal{M}^{t}\f$
	HOMomentsCPU*	moments_old;
	/// root HO Moments for time \f$\mathcal{M}^{t+1}\f$
	HOMomentsCPU* 	moments_next;


	// Device Specific Moment tallies
	HOMomentsCPU* cpu_moments;
	HOMomentsGPU* gpu_moments;
	HOMomentsMIC* mic_moments;


	PlasmaData* pdata;

	int nx,ny,nz,nspecies;



	void allocate(PlasmaData* _pdata);

	void apply_weights();

	void reduce();

	void UpdateMoments();

	realkind& get_val(int i, int j, int k, int ispecies,
			enum HOMoments_moment moment,int itime = 1);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Check charge conservation

		@param[in] moments_old Moments from time \f$t\f$, calling moment is at
		time \f$t+1\f$

		This method performs the following charge conservation check.

		\f{eqnarray}{
		CC^{t+1} = \frac{1}{n_x n_y n_z}\sum_{i,j,k,l}\frac{n^{t+1}_{i,j,k,l}-n^{t}_{i,j,k,l}}{\Delta t}
		+ \frac{\overline{nu}^{t+1/2}_{i+1/2,j,k,l}-\overline{nu}^{t+1/2}_{i-1/2,j,k,l}}{\Delta x}
		+ \frac{\overline{nv}^{t+1/2}_{i,j+1/2,k,l}-\overline{nv}^{t+1/2}_{i,j-1/2,k,l}}{\Delta y}
		+ \frac{\overline{nw}^{t+1/2}_{i,j,k+1/2,l}-\overline{nw}^{t+1/2}_{i,j,k-1/2,l}}{\Delta z}
		\f}
	*/
	/*--------------------------------------------------------------------------*/
	__host__
	double check_charge();

	/*-------------------------------------------------------------------------*/
	/**
		@brief Calculate and return the total particle kinetic energy from
		the stress tensor.

		@result Total system kinetic energy

		This function evaluates the following equation in order to return the
		total plasma kinetic energy.

		\f{eqnarray}{
		KE^{t+1} = \frac{1}{\Delta x\Delta y}\sum_{i,j,k,l}\left(S_{xx,i,j,k,l}^{HO,t+1}
					+S_{yy,i,j,k,l}^{HO,t+1}
					+S_{zz,i,j,k,l}^{HO,t+1}\right)
		\f}

	*/
	/*--------------------------------------------------------------------------*/
	__host__
	double evaluate_energy(void);

	__host__
	double evaluate_energy(int iS);

	double check_momentum(FieldDataCPU* fields_next, FieldDataCPU* fields_old);


	void set_vals(realkind val_in);


};


#endif /* NODEHOMOMENTS_H_ */
