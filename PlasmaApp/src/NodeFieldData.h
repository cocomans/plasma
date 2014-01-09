/*
 * NodeFieldData.h
 *
 *  Created on: Apr 22, 2013
 *      Author: payne
 */

#ifndef NODEFIELDDATA_H_
#define NODEFIELDDATA_H_

#include "FieldDataCPU.h"
#include "FieldDataGPU.h"
#include "FieldDataMIC.h"

class CPUTimer;

class NodeFieldData
{
public:


	FieldDataCPU* cpu_fields;
	FieldDataGPU* gpu_fields;
	FieldDataMIC* mic_fields;

	PlasmaData* pdata;
	int nx,ny,nz;
	int alloc_size;

	CPUTimer* bcast_timer;

	void allocate(PlasmaData* _pdata);

	void allocate(PlasmaData* _pdata,NodeFieldData* _fields);

	/*-------------------------------------------------------------------------*/
	/**
		@brief Compute the average of two fields

		@param[in] a First field to average
		@param[in] b Second field to average
		@param[in,out] c Field to store the average of a and b.
	*/
	/*--------------------------------------------------------------------------*/
	void average(FieldDataCPU* field0, FieldDataCPU* field1);

	void broadcast();

	double evaluate_energy();


};



#endif /* NODEFIELDDATA_H_ */
