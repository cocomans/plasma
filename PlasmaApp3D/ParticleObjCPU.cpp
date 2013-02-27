
#include "ParticleObjCPU.h"
#include "FieldData.h"
#include "PlasmaData.h"
#include "vec_funcs.h"
#include "CurrentTally.h"

void ParticleObjCPU::push(PlasmaData* pdata, FieldData* fields,
		CurrentTally* currents, const int nSubcycl_max)
{

	dt_finished = 0.0;
	int iter = 0;


	// Begin Subcycle
	while(dt_finished < pdata->dt)
	{

		// Estimate dtau
		float dtau = estimate_dtau(pdata,fields);


		// Begin Picard - Crank-Nicholson
		ParticleObj new_ptcl = PicardIterate(pdata,fields,dtau);


		// Check Cell crossing

		// Handle Cell crossing

		// Accumulate Current
		accumulate_current(pdata,currents);



		dt_finished += dtau;

		iter++;

		if(iter >= nSubcycl_max)
			break;
	}

}

ParticleObjCPU ParticleObjCPU::PicardIterate(PlasmaData* pdata, FieldData* fields, const float dtau)
{

	int iter;

	float3 x_half;
	float3 v_half;

	float3 x_next = make_float3(px,py,pz);
	float3 v_next = make_float3(vx,vy,vz);

	float3 x_last = make_float3(px,py,pz);
	float3 v_last = make_float3(vx,vy,vz);



	float residual = 10.0f * pdata->picard_tolerance;

	int niter_max = pdata -> niter_max;

	iter = 0;

	while(residual > pdata->picard_tolerance)
	{

		v_half = 0.5f * (v_next + v_last);


		x_next = __fmad(dtau, v_next, make_float3(px,py,pz));

		x_half = 0.5f * (x_next + x_last);

		float3 accel = fields->get_accel(x_half,v_half);

		v_next = __fmad(dtau,accel,make_float3(vx,vy,vz));

		residual = l1norm(v_next - v_last)/l1norm(v_next);


		iter++;

		if(iter >= niter_max)
			break;

	}

	ParticleObj new_ptcl = *this;

	new_ptcl.px = x_next.x;
	new_ptcl.py = x_next.y;
	new_ptcl.pz = x_next.z;

	new_ptcl.vx = v_next.x;
	new_ptcl.vy = v_next.y;
	new_ptcl.vz = v_next.z;

	return new_ptcl;

}


float3 ParticleObjCPU::calc_div_a(PlasmaData* pdata, FieldData* fields,const float3& accel)
{
	// Returns divergence of acceleration dotted into velocity

	/*
	 * Tensor form:
	 * [0].x [1].x [2].x
	 * [0].y [1].y [2].y
	 * [0].z [1].z [2].z
	 */

	float3 divB[3];

	float3 div_vxB[3];


	float3 Bvec = make_float3(
				fields->intrpB<0,FieldData_deriv_f>(px,py,pz),
				fields->intrpB<1,FieldData_deriv_f>(px,py,pz),
				fields->intrpB<2,FieldData_deriv_f>(px,py,pz));

	float3 velocity = make_float3(vx,vy,vz);

	float3 axB = cross_product(accel,Bvec);


	// v x B jacobian

	float3 div_a[3];
	div_a[0] = 1.0f/vx * axB;
	div_a[1] = 1.0f/vy * axB;
	div_a[2] = 1.0f/vz * axB;

	divB[0] = make_float3(
				fields->intrpB<0,FieldData_deriv_dfdx>(px,py,pz),
				fields->intrpB<1,FieldData_deriv_dfdx>(px,py,pz),
				fields->intrpB<2,FieldData_deriv_dfdx>(px,py,pz));

	divB[1] = make_float3(
				fields->intrpB<0,FieldData_deriv_dfdy>(px,py,pz),
				fields->intrpB<1,FieldData_deriv_dfdy>(px,py,pz),
				fields->intrpB<2,FieldData_deriv_dfdy>(px,py,pz));

	divB[2] = make_float3(
				fields->intrpB<0,FieldData_deriv_dfdz>(px,py,pz),
				fields->intrpB<1,FieldData_deriv_dfdz>(px,py,pz),
				fields->intrpB<2,FieldData_deriv_dfdz>(px,py,pz));

	div_vxB[0] = cross_product(velocity,divB[0]);
	div_vxB[1] = cross_product(velocity,divB[1]);
	div_vxB[2] = cross_product(velocity,divB[2]);


	div_a[0] += div_vxB[0];
	div_a[1] += div_vxB[1];
	div_a[2] += div_vxB[2];

	// Electric Field Jacobian

	float3 divE[3];

	divE[0] = make_float3(
				fields->intrpE<0,FieldData_deriv_dfdx>(px,py,pz),
				fields->intrpE<1,FieldData_deriv_dfdx>(px,py,pz),
				fields->intrpE<2,FieldData_deriv_dfdx>(px,py,pz));

	divE[1] = make_float3(
				fields->intrpE<0,FieldData_deriv_dfdy>(px,py,pz),
				fields->intrpE<1,FieldData_deriv_dfdy>(px,py,pz),
				fields->intrpE<2,FieldData_deriv_dfdy>(px,py,pz));

	divE[2] = make_float3(
				fields->intrpE<0,FieldData_deriv_dfdz>(px,py,pz),
				fields->intrpE<1,FieldData_deriv_dfdz>(px,py,pz),
				fields->intrpE<2,FieldData_deriv_dfdz>(px,py,pz));

	// dot div_a into velocity
	float3 result;

	result.x = div_a[0].x*velocity.x + div_a[1].x*velocity.x + div_a[2].x*velocity.x;
	result.y = div_a[0].y*velocity.y + div_a[1].y*velocity.y + div_a[2].y*velocity.y;
	result.z = div_a[0].z*velocity.z + div_a[1].z*velocity.z + div_a[2].z*velocity.z;

	return result;




}


float ParticleObjCPU::estimate_dtau(PlasmaData* pdata, FieldData* fields)
{
	float dtau;

	float3 accel;
	float3 div_a;
	float alpha,beta,gamma2;

	// Calculate acceleration
	accel = fields->get_accel(make_float3(px,py,pz),make_float3(vx,vy,vz));

	// Calculate div_a dot v
	div_a = calc_div_a(pdata,fields,accel);

	// Using the l1norm instead of the l2norm for speed reasons
	float norm_diva = l1norm(div_a);
	float norm_a = l1norm(accel);

	alpha = (abs(norm_diva) + abs(norm_a))/2;
	beta = pdata->epsilon_r * (abs(norm_a) + l1norm(make_float3(vx,vy,vz)));
	gamma2 = pdata->epsilon_a;

	float d = 1.0f/sqrt(alpha);

	dtau = fmax(gamma2*d,beta*d*d);

	return dtau;

}

void ParticleObjCPU::accumulate_current(PlasmaData* pdata, CurrentTally* currents)
{
	currents->tally(make_float3(px,py,pz),make_float3(vx,vy,vz),make_int3(ix,iy,iz));
}
















































