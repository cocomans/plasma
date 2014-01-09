/*
 * DiscretePDE1D3V.cpp
 *
 *  Created on: Dec 12, 2013
 *      Author: payne
 */

#include "DiscretePDE1D3V.h"
#include "ScalarArray.h"
#include "StateVector.h"
#include "StressTensor.h"

namespace HoLo
{
namespace LowOrder
{
const int DiscretePDE1D3V::pdeVar2genVar[SimParams::Vars::LF] = {
		Dep::ne,Dep::pe_u,Dep::pe_v,Dep::pe_w,
		Dep::ni,Dep::pi_u,Dep::pi_v,Dep::pi_w,
		-Aux::Eu,-Aux::Ev,-Aux::Ew,
		-Aux::Bu,-Aux::Bv,-Aux::Bw,
		-Aux::Au,Dep::Av,Dep::Aw,
		Dep::phi
};
void DiscretePDE1D3V::UpdateAuxillary(StateVector** Q)
{
	using namespace Util::PDE;

	auto Av_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetDep(Dep::Av,i,j,k);};
	auto Aw_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetDep(Dep::Aw,i,j,k);};

	auto phi_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetDep(Dep::phi,i,j,k);};

	auto Av_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetDep(Dep::Av,i,j,k);};
	auto Aw_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetDep(Dep::Aw,i,j,k);};

	auto phi_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetDep(Dep::phi,i,j,k);};

	auto Eu_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetAux(Aux::Eu,i,j,k);};

	auto Ev_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetAux(Aux::Ev,i,j,k);};
	auto Ew_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetAux(Aux::Ew,i,j,k);};

	auto Bv_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetAux(Aux::Bv,i,j,k);};

	auto Bw_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetAux(Aux::Bw,i,j,k);};

	auto Eu_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetAux(Aux::Eu,i,j,k);};

	auto Ev_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetAux(Aux::Ev,i,j,k);};
	auto Ew_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetAux(Aux::Ew,i,j,k);};

	auto Bv_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetAux(Aux::Bv,i,j,k);};

	auto Bw_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetAux(Aux::Bw,i,j,k);};



	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Av_oF)> Av_o(Av_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Aw_oF)> Aw_o(Aw_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(phi_oF)> phi_o(phi_oF);

	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Av_nF)> Av_n(Av_nF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Aw_nF)> Aw_n(Aw_nF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(phi_nF)> phi_n(phi_nF);

	ScalarArray<real,StagMesh::F,StagMesh::C,StagMesh::C,decltype(Eu_oF)> Eu_o(Eu_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Ev_oF)> Ev_o(Ev_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Ew_oF)> Ew_o(Ew_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Bv_oF)> Bv_o(Bv_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Bw_oF)> Bw_o(Bw_oF);


	for(int i=-1;i<nx+1;i++)
	{
		Q[0]->GetAux(Aux::Eu,i,0,0) = Du(phi_n)(i,0,0);
		Q[0]->GetAux(Aux::Ev,i,0,0) = ((-2.0/dt)*(Av_n - Av_o) - Ev_o)(i,0,0);
		Q[0]->GetAux(Aux::Ew,i,0,0) = ((-2.0/dt)*(Aw_n - Aw_o) - Ew_o)(i,0,0);
		Q[0]->GetAux(Aux::Bv,i,0,0) = (-1.0*Du(Aw_n))(i,0,0);
		Q[0]->GetAux(Aux::Bw,i,0,0) = (Du(Av_n))(i,0,0);
	}






}

double DiscretePDE1D3V::FieldResidual(StateVector* Res)
{

	int FieldIDs[3] = {Dep::phi,Dep::Av,Dep::Aw};
	double result = 0;

	for(int l=0;l<3;l++)
	for(int i=0;i<nx;i++)
	{
		result += pow(Res->GetDep(FieldIDs[l],i,0,0),2);
	}

	// Zero out the non-field residuals
	for(int l=0;l<Dep::LF;l++)
	{
		if((l != Dep::phi)||(l != Dep::Av)||(l != Dep::Aw))
			for(int i=0;i<nx;i++)
			{
				Res->GetDep(FieldIDs[l],i,0,0) = 0;
			}
	}

	return sqrt(result)/((double)nx);
}

void DiscretePDE1D3V::EvaluatePDE(StateVector** Q,
			StressTensor** stress,
			StateVector* Gamma,
			StateVector* Res)
{
	using namespace Util::PDE;
	// Setting up lambdas for pde access. DO NOT TOUCH!!!!!
	auto ne_oF = [&](const int& i, const int &j, const int& k){return Q[1]->GetDep(Dep::ne,i,j,k);};
	auto ni_oF = [&](const int& i, const int &j, const int& k){return Q[1]->GetDep(Dep::ni,i,j,k);};

	auto pe_uoF = [&](const int& i, const int &j, const int& k){return Q[1]->GetDep(Dep::pe_u,i,j,k);};
	auto pe_voF = [&](const int& i, const int &j, const int& k){return Q[1]->GetDep(Dep::pe_v,i,j,k);};
	auto pe_woF = [&](const int& i, const int &j, const int& k){return Q[1]->GetDep(Dep::pe_w,i,j,k);};

	auto pi_uoF = [&](const int& i, const int &j, const int& k){return Q[1]->GetDep(Dep::pi_u,i,j,k);};
	auto pi_voF = [&](const int& i, const int &j, const int& k){return Q[1]->GetDep(Dep::pi_v,i,j,k);};
	auto pi_woF = [&](const int& i, const int &j, const int& k){return Q[1]->GetDep(Dep::pi_w,i,j,k);};

	auto ne_nF = [&](const int& i, const int &j, const int& k){return Q[0]->GetDep(Dep::ne,i,j,k);};
	auto ni_nF = [&](const int& i, const int &j, const int& k){return Q[0]->GetDep(Dep::ni,i,j,k);};

	auto pe_unF = [&](const int& i, const int &j, const int& k){return Q[0]->GetDep(Dep::pe_u,i,j,k);};
	auto pe_vnF = [&](const int& i, const int &j, const int& k){return Q[0]->GetDep(Dep::pe_v,i,j,k);};
	auto pe_wnF = [&](const int& i, const int &j, const int& k){return Q[0]->GetDep(Dep::pe_w,i,j,k);};

	auto pi_unF = [&](const int& i, const int &j, const int& k){return Q[0]->GetDep(Dep::pi_u,i,j,k);};
	auto pi_vnF = [&](const int& i, const int &j, const int& k){return Q[0]->GetDep(Dep::pi_v,i,j,k);};
	auto pi_wnF = [&](const int& i, const int &j, const int& k){return Q[0]->GetDep(Dep::pi_w,i,j,k);};

	auto suu_eoF = [&](const int& i, const int &j, const int& k)
			{return stress[1]->Get(StressTensor::Suu,StressTensor::electron,i,j,k)*ne_oF(i,j,k);};
	auto suv_eoF = [&](const int& i, const int &j, const int& k)
			{return stress[1]->Get(StressTensor::Suv,StressTensor::electron,i,j,k)*ne_oF(i,j,k);};
	auto suw_eoF = [&](const int& i, const int &j, const int& k)
			{return stress[1]->Get(StressTensor::Suw,StressTensor::electron,i,j,k)*ne_oF(i,j,k);};

	auto suu_ioF = [&](const int& i, const int &j, const int& k)
			{return stress[1]->Get(StressTensor::Suu,StressTensor::ion,i,j,k)*ni_oF(i,j,k);};
	auto suv_ioF = [&](const int& i, const int &j, const int& k)
			{return stress[1]->Get(StressTensor::Suv,StressTensor::ion,i,j,k)*ni_oF(i,j,k);};
	auto suw_ioF = [&](const int& i, const int &j, const int& k)
			{return stress[1]->Get(StressTensor::Suw,StressTensor::ion,i,j,k)*ni_oF(i,j,k);};

	auto suu_enF = [&](const int& i, const int &j, const int& k)
			{return stress[0]->Get(StressTensor::Suu,StressTensor::electron,i,j,k)*ne_nF(i,j,k);};
	auto suv_enF = [&](const int& i, const int &j, const int& k)
			{return stress[0]->Get(StressTensor::Suv,StressTensor::electron,i,j,k)*ne_nF(i,j,k);};
	auto suw_enF = [&](const int& i, const int &j, const int& k)
			{return stress[0]->Get(StressTensor::Suw,StressTensor::electron,i,j,k)*ne_nF(i,j,k);};

	auto suu_inF = [&](const int& i, const int &j, const int& k)
			{return stress[0]->Get(StressTensor::Suu,StressTensor::ion,i,j,k)*ni_nF(i,j,k);};
	auto suv_inF = [&](const int& i, const int &j, const int& k)
			{return stress[0]->Get(StressTensor::Suv,StressTensor::ion,i,j,k)*ni_nF(i,j,k);};
	auto suw_inF = [&](const int& i, const int &j, const int& k)
			{return stress[0]->Get(StressTensor::Suw,StressTensor::ion,i,j,k)*ni_nF(i,j,k);};



	auto Av_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetDep(Dep::Av,i,j,k);};
	auto Aw_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetDep(Dep::Aw,i,j,k);};

	auto phi_oF = [&](const int& i, const int &j, const int& k)
			{return Q[1]->GetDep(Dep::phi,i,j,k);};

	auto Av_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetDep(Dep::Av,i,j,k);};
	auto Aw_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetDep(Dep::Aw,i,j,k);};

	auto phi_nF = [&](const int& i, const int &j, const int& k)
			{return Q[0]->GetDep(Dep::phi,i,j,k);};



	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(ne_oF)> ne_o(ne_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(ne_nF)> ne_n(ne_nF);

	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(ni_oF)> ni_o(ni_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(ni_nF)> ni_n(ni_nF);

	ScalarArray<real,StagMesh::F,StagMesh::C,StagMesh::C,decltype(pe_uoF)> pe_uo(pe_uoF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(pe_voF)> pe_vo(pe_voF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(pe_woF)> pe_wo(pe_woF);

	ScalarArray<real,StagMesh::F,StagMesh::C,StagMesh::C,decltype(pi_uoF)> pi_uo(pi_uoF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(pi_voF)> pi_vo(pi_voF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(pi_woF)> pi_wo(pi_woF);


	ScalarArray<real,StagMesh::F,StagMesh::C,StagMesh::C,decltype(pe_unF)> pe_un(pe_unF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(pe_vnF)> pe_vn(pe_vnF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(pe_wnF)> pe_wn(pe_wnF);

	ScalarArray<real,StagMesh::F,StagMesh::C,StagMesh::C,decltype(pi_unF)> pi_un(pi_unF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(pi_vnF)> pi_vn(pi_vnF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(pi_wnF)> pi_wn(pi_wnF);


	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Av_oF)> Av_o(Av_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Aw_oF)> Aw_o(Aw_oF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(phi_oF)> phi_o(phi_oF);

	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Av_nF)> Av_n(Av_nF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(Aw_nF)> Aw_n(Aw_nF);
	ScalarArray<real,StagMesh::C,StagMesh::C,StagMesh::C,decltype(phi_nF)> phi_n(phi_nF);



#pragma omp parallel for
	for(int i=0;i<Q[0]->nx;i++)
	{

		Res->GetDep(Dep::ne,i,0,0) = ((ne_n - ne_o) + dt*Du(pe_un))(i,0,0);
		Res->GetDep(Dep::ni,i,0,0) = ((ni_n - ni_o) + dt*Du(pi_un))(i,0,0);

		Res->GetDep(Dep::pe_u,i,0,0) = ((pe_un - pe_uo) + dt*(Te/me * Du(ne_n) + qe/me*Du(phi_n)))(i,0,0);
		Res->GetDep(Dep::pe_v,i,0,0) = ((pe_un - pe_uo) + dt*(Te/me * Du(ne_n) + qe/me*Du(phi_n)))(i,0,0);
		Res->GetDep(Dep::pe_w,i,0,0) = ((pe_un - pe_uo) + dt*(Te/me * Du(ne_n) + qe/me*Du(phi_n)))(i,0,0);

		Res->GetDep(Dep::pi_u,i,0,0) = ((pe_un - pe_uo) + dt*(Te/me * Du(ne_n) + qe/me*Du(phi_n)))(i,0,0);
		Res->GetDep(Dep::pi_v,i,0,0) = ((pe_un - pe_uo) + dt*(Te/me * Du(ne_n) + qe/me*Du(phi_n)))(i,0,0);
		Res->GetDep(Dep::pi_w,i,0,0) = ((pe_un - pe_uo) + dt*(Te/me * Du(ne_n) + qe/me*Du(phi_n)))(i,0,0);

		Res->GetDep(Dep::phi,i,0,0) = (Du(Du((phi_n-phi_o))) + dt*Du(qe*pe_un + qi*pi_un))(i,0,0);
		Res->GetDep(Dep::Av,i,0,0) = (0.5*Du(Du(Av_n+Av_o)) - dt*(qe*pe_vn + qi*pi_vn))(i,0,0);
		Res->GetDep(Dep::Aw,i,0,0) = (0.5*Du(Du(Aw_n+Aw_o)) - dt*(qe*pe_wn + qi*pi_wn))(i,0,0);
	}

// Before....
//    double tmp_W = (OvX->Get(Sxx_e, i-1)+3.0*OvXold->Get(Sxx_e, i-1))*(OvX->Get(n_e,i-1)+3.0*OvXold->Get(n_e,i-1));
//     double tmp_C = (OvX->Get(Sxx_e, i)+3.0*OvXold->Get(Sxx_e, i)) * (OvX->Get(n_e,i)+3.0*OvXold->Get(n_e,i));
//     double tmp_E = (OvX->Get(Sxx_e, i+1)+3.0*OvXold->Get(Sxx_e, i+1)) * (OvX->Get(n_e,i+1)+3.0*OvXold->Get(n_e,i+1));
//     double Dx_Sxx_ne = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//     tmp_W = (OvX->Get(Sxy_e, i-1)+3.*OvXold->Get(Sxy_e, i-1)) * (OvX->Get(n_e,i-1)+3.*OvXold->Get(n_e,i-1));
//     tmp_C = (OvX->Get(Sxy_e, i)+3.*OvXold->Get(Sxy_e, i)) * (OvX->Get(n_e,i)+3.*OvXold->Get(n_e,i));
//     tmp_E = (OvX->Get(Sxy_e, i+1)+3.*OvXold->Get(Sxy_e, i+1)) * (OvX->Get(n_e,i+1)+3.*OvXold->Get(n_e,i+1));
//     double Dx_Sxy_ne = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//     tmp_W = (OvX->Get(Sxz_e, i-1)+3.*OvXold->Get(Sxz_e, i-1)) * (OvX->Get(n_e,i-1)+3.*OvXold->Get(n_e,i-1));
//     tmp_C = (OvX->Get(Sxz_e, i)+3.*OvXold->Get(Sxz_e, i)) * (OvX->Get(n_e,i)+3.*OvXold->Get(n_e,i));
//     tmp_E = (OvX->Get(Sxz_e, i+1)+3.*OvXold->Get(Sxz_e, i+1)) * (OvX->Get(n_e,i+1)+3.*OvXold->Get(n_e,i+1));
//     double Dx_Sxz_ne = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//     // Ion
//     tmp_W = (OvX->Get(Sxx_i, i-1)+3.*OvXold->Get(Sxx_i, i-1)) * (OvX->Get(n_i,i-1)+3.*OvXold->Get(n_i,i-1));
//     tmp_C = (OvX->Get(Sxx_i, i)+3.*OvXold->Get(Sxx_i, i)) * (OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i));
//     tmp_E = (OvX->Get(Sxx_i, i+1)+3.*OvXold->Get(Sxx_i, i+1)) * (OvX->Get(n_i,i+1)+3.*OvXold->Get(n_i,i+1));
//     double Dx_Sxx_ni = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//     tmp_W = (OvX->Get(Sxy_i, i-1)+3.*OvXold->Get(Sxy_i, i-1)) * (OvX->Get(n_i,i-1)+3.*OvXold->Get(n_i,i-1));
//     tmp_C = (OvX->Get(Sxy_i, i)+3.*OvXold->Get(Sxy_i, i)) * (OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i));
//     tmp_E = (OvX->Get(Sxy_i, i+1)+3.*OvXold->Get(Sxy_i, i+1)) * (OvX->Get(n_i,i+1)+3.*OvXold->Get(n_i,i+1));
//     double Dx_Sxy_ni = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//     tmp_W = (OvX->Get(Sxz_i, i-1)+3.*OvXold->Get(Sxz_i, i-1)) * (OvX->Get(n_i,i-1)+3.*OvXold->Get(n_i,i-1));
//     tmp_C = (OvX->Get(Sxz_i, i)+3.*OvXold->Get(Sxz_i, i)) * (OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i));
//     tmp_E = (OvX->Get(Sxz_i, i+1)+3.*OvXold->Get(Sxz_i, i+1)) * (OvX->Get(n_i,i+1)+3.*OvXold->Get(n_i,i+1));
//     double Dx_Sxz_ni = 0.0625*OvX->DxC2F(tmp_W, tmp_C, tmp_E, dx);
//
//
//
//     // Fill each unknown
//     //    Electron continuity equation
//     res[j+n_e] = (OvX->Get(n_e, i) - OvXold->Get(n_e, i)) + dt*( OvX->DxF2C(px_e, i, dx)) - gamma[j + n_e];
//     //    Ion continuity equation
//     res[j+n_i] = (OvX->Get(n_i, i) - OvXold->Get(n_i, i)) + dt*( OvX->DxF2C(px_i, i, dx) ) - gamma[j + n_i];
//
//     //    Electron momentum equation
//     //    momentum electron-x
//     res[j+px_e] = (OvX->Get(px_e, i) - OvXold->Get(px_e, i)) +
//                         (dt/(2.0*me_h))*(
//                                 me_h*Dx_Sxx_ne -
//                                 0.0625*qe_h*(OvX->GetC2F(n_e,i)+3.*OvXold->GetC2F(n_e,i))*(3.*OvXold->Get(Ex,i)+Ext) -
//                                 (
//                                         0.125*qe_h*(OvX->GetC2F(py_e,i)+OvXold->GetC2F(py_e,i))*(OvX->DxC2F(Ay,i,dx) + 3.*OvXold->DxC2F(Ay,i,dx))+
//                                         0.125*qe_h*(OvX->GetC2F(pz_e,i)+OvXold->GetC2F(pz_e,i))*(OvX->DxC2F(Az,i,dx) + 3.*OvXold->DxC2F(Az,i,dx))
//                                 	  )
//                                 )
//                                  - gamma[j + px_e];
//
//     //    momentum electron-y
//     res[j+py_e] = (OvX->Get(py_e, i) - OvXold->Get(py_e, i)) +
//                         (dt/(2.0*me_h))*(
//                                 me_h*Dx_Sxy_ne -
//                                 0.0625*qe_h*(OvX->Get(n_e,i)+3.*OvXold->Get(n_e,i))*(3.*OvXold->Get(Ey,i)+Eyt) -
//                                 (
//                               		  qe_h*OvX->Get(pz_e,i)*OvXold->Get(Bx,i)-
//                                         0.125*qe_h*(OvX->GetF2C(px_e,i)+OvXold->GetF2C(px_e,i))*(OvX->Dx(Ay,i,dx) + 3.*OvXold->Dx(Ay,i,dx))
//                                 	  )
//                                 )
//                                  - gamma[j + py_e];
//
//     //    momentum electron-z
//     res[j+pz_e] = (OvX->Get(pz_e, i) - OvXold->Get(pz_e, i)) +
//                         (dt/(2.0*me_h))*(
//                                 me_h*Dx_Sxz_ne -
//                                 0.0625*qe_h*(OvX->Get(n_e,i)+3.*OvXold->Get(n_e,i))*(3.*OvXold->Get(Ez,i)+Ezt) -
//                                 (
//                               		  -qe_h*OvX->Get(py_e,i)*OvXold->Get(Bx,i)-
//                                         0.125*qe_h*(OvX->GetF2C(px_e,i)+OvXold->GetF2C(px_e,i))*(OvX->Dx(Az,i,dx) + 3.*OvXold->Dx(Az,i,dx))
//                                 	  )
//                                 )
//                                  - gamma[j + pz_e];
//
//
//
//     //    Ion momentum equation
//     //    momentum ion-x
//     res[j+px_i] = (OvX->Get(px_i, i) - OvXold->Get(px_i, i)) +
//                         (dt/(2.0*mi_h))*(
//                                 mi_h*Dx_Sxx_ni -
//                                 0.0625*qi_h*(OvX->GetC2F(n_i,i)+3.*OvXold->GetC2F(n_i,i))*(3.*OvXold->Get(Ex,i)+Ext) -
//                                 (
//                                         0.125*qi_h*(OvX->GetC2F(py_i,i)+OvXold->GetC2F(py_i,i))*(OvX->DxC2F(Ay,i,dx) + 3.*OvXold->DxC2F(Ay,i,dx)) +
//                                         0.125*qi_h*(OvX->GetC2F(pz_i,i)+OvXold->GetC2F(pz_i,i))*(OvX->DxC2F(Az,i,dx) + 3.*OvXold->DxC2F(Az,i,dx))
//                                 	  )
//                                 )
//                                  - gamma[j + px_i];
//
//     //    momentum electron-y
//     res[j+py_i] = (OvX->Get(py_i, i) - OvXold->Get(py_i, i)) +
//                         (dt/(2.0*mi_h))*(
//                                 mi_h*Dx_Sxy_ni -
//                                 0.0625*qi_h*(OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i))*(3.*OvXold->Get(Ey,i)+Eyt) -
//                                 (
//                               		  qi_h*OvX->Get(pz_i,i)*OvXold->Get(Bx,i)-
//                                         0.125*qi_h*(OvX->GetF2C(px_i,i)+OvXold->GetF2C(px_i,i))*(OvX->Dx(Ay,i,dx) + 3.*OvXold->Dx(Ay,i,dx))
//                                 	  )
//                                 )
//                                  - gamma[j + py_i];
//
//     //    momentum electron-z
//     res[j+pz_i] = (OvX->Get(pz_i, i) - OvXold->Get(pz_i, i)) +
//                         (dt/(2.0*mi_h))*(
//                                 mi_h*Dx_Sxz_ni -
//                                 0.0625*qi_h*(OvX->Get(n_i,i)+3.*OvXold->Get(n_i,i))*(3.*OvXold->Get(Ez,i)+Ezt) -
//                                 (
//                               		  -qi_h*OvX->Get(py_i,i)*OvXold->Get(Bx,i)-
//                                         0.125*qi_h*(OvX->GetF2C(px_i,i)+OvXold->GetF2C(px_i,i))*(OvX->Dx(Az,i,dx) + 3.*OvXold->Dx(Az,i,dx))
//                                 	  )
//                                 )
//                                  - gamma[j + pz_i];
//    //    Ey
//    res[j+Ey] = (0.5*dt*(OvX->Get(Ey,i) + OvXold->Get(Ey,i)) +
//                        (OvX->Get(Ay,i) - OvXold->Get(Ay,i)));
//    //    Ez
//    res[j+Ez] = (0.5*dt*(OvX->Get(Ez,i) + OvXold->Get(Ez,i)) +
//                        (OvX->Get(Az,i) - OvXold->Get(Az,i)));
//
//    //    Ay
//    res[j+Ex]     = (dt*(qe_h*OvX->Get(px_e,i) + qi_h*OvX->Get(px_i,i) - jmeanx1)/(xi*xi)
//  				  + (OvX->Get(Ex,i) - OvXold->Get(Ex,i)));
//
//    //    Ay
//    res[j+Ay]     = (0.5*zeta*dt/(xi*xi)*(OvX->Dxx(Ay,i,dx)+OvXold->Dxx(Ay,i,dx))
//  				  + (dt*(qe_h*OvX->Get(py_e,i) + qi_h*OvX->Get(py_i,i) - jmeany1)/(xi*xi)));
//    //    Az
//    res[j+Az]     = (0.5*zeta*dt/(xi*xi)*(OvX->Dxx(Az,i,dx)+OvXold->Dxx(Az,i,dx))
//  				  + (dt*(qe_h*OvX->Get(pz_e,i) + qi_h*OvX->Get(pz_i,i) - jmeanz1)/(xi*xi)));
























}

} /* namespace LowOrder */
} /* namespace HoLo */
