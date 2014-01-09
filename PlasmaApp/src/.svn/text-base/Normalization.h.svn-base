/*-------------------------------------------------------------------------*/
/**
	@file		Normalization.h
	@author	J. Payne
	@date		08/29/2013
	@brief	Declares the Normalization class, a pure virtual class that is overloaded
	with specific normalization constant calculations.


*/
/*--------------------------------------------------------------------------*/
#ifndef NORMALIZATION_H
#define NORMALIZATION_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <gnuplot_i.h>

class PlasmaData;
class InputParser;

class Normalization
{
public:

	virtual ~Normalization(){};

	virtual void CalcNorms(PlasmaData* pdata,InputParser* parser){};
	// naught quantities (user inputed in si units)
	/// reference mass in kg
	realkind m0;
	/// Reference charge in C
	realkind q0;
	/// Reference Length in m
	realkind L0;
	/// Reference number density in atom/m^3
	realkind n0;
	/// Reference Temperature in Joules
	realkind T0;

	/// System Lengths in m
	realkind Lx_p;
	realkind Ly_p;
	realkind Lz_p;
	/// Electron temperature in physical units (same units as T0)
	realkind Te_p;
	/// ion temperature in physical units (same units as T0)
	realkind Ti_p;
	/// Electron density in physical units (same units as n0)
	realkind ne_p;
	/// Ion density in physical units (same units as n0)
	realkind ni_p;

	realkind Lx,Ly,Lz;

	int nptcls_species[10];


};




#endif /* NORMALIZATION_H */
