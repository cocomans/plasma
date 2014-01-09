/*-------------------------------------------------------------------------*/
/**
	@file		NormElectronEM.cpp
	@author	J. Payne
	@date		08/29/2013
	@brief	Defines the NormElectronEM class, defines the Electron-based
	electro-magnetic normalization.


*/
/*--------------------------------------------------------------------------*/

#include "../PlasmaData.h"
#include "NormElectronEM.h"



void NormElectronEM::CalcNorms(PlasmaData* pdata, InputParser* parser)
{

	double Ne = 256000;
	double Ni = pdata->nptcls_species[1];
	double c2 = 1.0/(sqrt(epsilon_0_p*mu_0_p));

	parser->GetParam("NORMALIZED_SYSTEM_LENGTH_X","-Lx",Lx,pi_const);
	parser->GetParam("NORMALIZED_SYSTEM_LENGTH_Y","-Ly",Ly,1.0);
	parser->GetParam("NORMALIZED_SYSTEM_LENGTH_Z","-Lz",Lz,1.0);

	parser->GetParam("LENGTH_NORMALIZER","-L0",L0,c2*c2*me_p*epsilon_0_p*Lx/(qe_p*qe_p*Ne));



	parser->GetParam("PHYSICAL_ELECTRON_DENSITY","-ne_p",ne_p,Ne/(Lx*L0));
	parser->GetParam("PHYSICAL_ION_DENSITY","-ni_p",ni_p,ne_p);

	parser->GetParam("MASS_NORMALIZER","-m0",m0,me_p*L0*ne_p);
	parser->GetParam("CHARGE_NORMALIZER","-q0",q0,-qe_p*L0*ne_p);

	parser->GetParam("DENSITY_NORMALIZER","-n0",n0,1.0/L0);
	parser->GetParam("TEMPERATURE_NORMALIZER","-T0",T0,m0*c2*c2);

	parser->GetParam("PHYSICAL_ELECTRON_TEMPERATURE","-Te_p",Te_p,T0);
	parser->GetParam("PHYSICAL_ION_TEMPERATURE","-Ti_p",Ti_p,T0);


	parser->GetParam("PHYSICAL_SYSTEM_LENGTH_X","-Lx_p",Lx_p,L0*Lx);
	parser->GetParam("PHYSICAL_SYSTEM_LENGTH_Y","-Ly_p",Ly_p,L0*Ly);
	parser->GetParam("PHYSICAL_SYSTEM_LENGTH_Z","-Lz_p",Lz_p,L0*Lz);


//	parser->GetParam("NORMALIZED_SYSTEM_LENGTH_X","-Lx",Lx,pi_const);
//	parser->GetParam("NORMALIZED_SYSTEM_LENGTH_Y","-Ly",Ly,1.0);
//	parser->GetParam("NORMALIZED_SYSTEM_LENGTH_Z","-Lz",Lz,1.0);
//
//	parser->GetParam("LENGTH_NORMALIZER","-L0",L0,3.46514e+08);
//
//
//	parser->GetParam("PHYSICAL_ELECTRON_DENSITY","-ne_p",ne_p,0.000235163);
//	parser->GetParam("PHYSICAL_ION_DENSITY","-ni_p",ni_p,ne_p);
//
//
//	parser->GetParam("MASS_NORMALIZER","-m0",m0,me_p*L0*ne_p);
//	parser->GetParam("CHARGE_NORMALIZER","-q0",q0,qe_p*L0*ne_p);
//
//	parser->GetParam("DENSITY_NORMALIZER","-n0",n0,1.0/L0);
//	parser->GetParam("TEMPERATURE_NORMALIZER","-T0",T0,m0*c2*c2);
//
//	parser->GetParam("PHYSICAL_ELECTRON_TEMPERATURE","-Te_p",Te_p,T0);
//	parser->GetParam("PHYSICAL_ION_TEMPERATURE","-Ti_p",Ti_p,T0);
//
//
//	parser->GetParam("PHYSICAL_SYSTEM_LENGTH_X","-Lx_p",Lx_p,L0*Lx);
//	parser->GetParam("PHYSICAL_SYSTEM_LENGTH_Y","-Ly_p",Ly_p,L0*Ly);
//	parser->GetParam("PHYSICAL_SYSTEM_LENGTH_Z","-Lz_p",Lz_p,L0*Lz);

	pdata->L0 = L0;
	pdata->Lx = Lx;
	pdata->Ly = Ly;
	pdata->ne_p = ne_p;
	pdata->ni_p = ni_p;
	pdata->m0 = m0;
	pdata->q0 = q0;
	pdata->n0 = n0;
	pdata->T0 = T0;
	pdata->Te_p = Te_p;
	pdata->Ti_p = Ti_p;
	pdata->Lx_p = Lx_p;
	pdata->Ly_p = Ly_p;
	pdata->Lz_p = Lz_p;

	printf("\n\n##############NORMS################\n\n");



}




