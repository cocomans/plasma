/*-------------------------------------------------------------------------*/
/**
	@file		NormElectronEM.h
	@author	J. Payne
	@date		08/29/2013
	@brief	Declares the NormElectronEM class, defines the Electron-based
	electro-magnetic normalization.


*/
/*--------------------------------------------------------------------------*/
#ifndef NORM_ELECTRON_EM_H
#define NORM_ELECTRON_EM_H

#include "../Normalization.h"



class NormElectronEM : public Normalization
{
public:

	~NormElectronEM(){};

	void CalcNorms(PlasmaData* pdata,InputParser* parser);
	// naught quantities (user inputed in si units)
//	/// reference mass in kg
//	realkind m0;
//	/// Reference charge in C
//	realkind q0;
//	/// Reference Length in m
//	realkind L0;
//	/// Reference number density in atom/m^3
//	realkind n0;
//	/// Reference Temperature in Joules
//	realkind T0;
//
//	/// System Lengths in m
//	realkind Lx_p;
//	realkind Ly_p;
//	realkind Lz_p;
//	/// Electron temperature in physical units (same units as T0)
//	realkind Te_p;
//	/// ion temperature in physical units (same units as T0)
//	realkind Ti_p;
//	/// Electron density in physical units (same units as n0)
//	realkind ne_p;
//	/// Ion density in physical units (same units as n0)
//	realkind ni_p;
//
//	realkind Lx,Ly,Lz;
//
//	int nptcls_species[10];


};




#endif /* NORM_ELECTRON_EM_H */
