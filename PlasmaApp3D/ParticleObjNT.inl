
#include "ParticleObjNT.h"
#include "FieldData.h"
#include "PlasmaData.h"
#include "vec_funcs.h"
#include "typevecN.h"
#include "CurrentTally.h"
#include "ShapeFunctions.h"



template<const int ileft,const int nSpatial,const int nVel,const bool iEM> __attribute__ ((noinline))
void PushNT(PlasmaData* 		pdata,
				 FieldData* 		fields,
				 CurrentTally* 		current,
				 ParticleList* 	plist,
				 int**				iter_array,
				 int**				iptcl,
				 int**				iptcl_new,
				 int&				nptcls_left,
				 int&				nptcl_done,
				 int&				nptcls_process,
				 long long int&				nSubSteps_done)
{
//	printf("Doing a shrink push\n");
	ParticleObjNT<ileft,nSpatial,nVel,iEM> particles2(*iptcl);
	typevecN<int,ileft> iter2;


	iter2 = *iter_array;

	// Copy from main array
	particles2 = *plist;
#ifndef GPU_CODE
	int tid = omp_get_thread_num();
	particles2.piccard_timer = plist->piccard_timer+tid;
	particles2.accel_timer = plist->accel_timer+tid;
	particles2.tally_timer = plist->tally_timer+tid;
	particles2.crossing_timer = plist->crossing_timer+tid;
	particles2.dtau_est_timer = plist->dtau_est_timer+tid;
#endif

	// Push
	particles2.push(pdata,fields,current,iter2,pdata->nSubcycle_max);

	// Write back to main array
	//printf("writing particles back\n");
	for(int j=0;j<ileft;j++)
	{
		//printf("iptcl[%i] = %i\n",j,iptcl[0][j]);
		particles2.write_back(*plist,j);
	}

	//*plist = particles2;

	// Condense iptcl list
	int k = 0;


	for(int j=0;j<ileft;j++)
	{
		bool idone = 0;

		//printf("iter2(%i) = %f\n",j,particles2.dt_finished(j));
		if(particles2.dt_finished(j) >= pdata->dt)
		{
			idone = 1;
		}
		else if(iter2(j) >= pdata->nSubcycle_max)
		{
//			printf("warning particle finished before time step was finished dt_left[%i] = %f\n",iptcl[0][j],pdata->dt-particles2.dt_finished(j));

			idone = 1;
		}
		else
			idone = 0;


		if(idone)
		{
			nSubSteps_done += iter2(j);
			plist->num_subcycles[iptcl[0][j]] += iter2(j);
			iter2(j) = 0;
			plist->dt_finished[particles2.pid[j]] = 0;

			// Accumulate Charge and S2 moment

		}
		else
		{
			iptcl[0][k] = iptcl[0][j];
			iter_array[0][k] = iter2(j);


			k++;
		}


	}




//	printf("nptcl_done, k = %i, %i\n",nptcl_done,k);
	nptcl_done = nptcls_process - k ;
	nptcls_left = k;


	//int* iptcl_temp = *iptcl;
	//*iptcl = *iptcl_new;
	//*iptcl_new = iptcl_temp;
}

template<const int ileft,const int nSpatial,const int nVel,const bool iEM> __attribute__ ((noinline))
void shrink_pushT(PlasmaData* 		pdata,
				 FieldData* 		fields,
				 CurrentTally* 		current,
				 ParticleList* 	plist,
				 int**				iter_array,
				 int**				iptcl,
				 int**				iptcl_new,
				 int&				nptcls_left,
				 int&				nptcl_done,
				 int&				nptcls_process,
				 long long int&				nSubSteps_done)
{
	//printf("shrink_push with %i ptcls\n",ileft);
if(nptcls_left > 0)
{
	if(nptcls_left == ileft)
	{
		PushNT<ileft,nSpatial,nVel,iEM>(pdata,fields,current,plist,
							iter_array,iptcl,iptcl_new,
							nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);

	} /* if(ileft == nptcls_left) */
	else if(ileft > 0)
	{
		shrink_pushT<ileft-1,nSpatial,nVel,iEM>(pdata,fields,current,plist,
							iter_array,iptcl,iptcl_new,
							nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);
	}
}


}


template<const int N,const int nSpatial,const int nVel,const bool iEM> __attribute__((noinline)) __device__
void ParticleObjNT<N,nSpatial,nVel,iEM>::push(PlasmaData* pdata, FieldData* fields,
		CurrentTally* currents, typevecN<int,N>& iter, const int nSubcycl_max)
{

	bool irun = 1;
	typevecN<typevecN<realkind,N>,nVel> accel;
	accel = get_accel(fields,position,velocity);

//	for(int i=0;i<N;i++)
//		printf("accel[%i] = %f(%f,%f)\n",pid[i],accel(0)(i),position(0)(i),velocity(0)(i));


	// Begin Subcycle
	while(irun)
	{

		// Estimate dtau
		//typevecN<realkind,N> dtau = estimate_dtau(pdata,fields);

		typevecN<realkind,N> dtau;

		dtau = estimate_dtau(pdata,fields);
		//dtau = 0.001f;


		// Begin Picard - Crank-Nicholson
		//printf("Picard Iteration %i\n",iter(0));
#ifndef GPU_CODE
		piccard_timer->start();
#endif
		PicardIterate(pdata,fields,currents,dtau);
#ifndef GPU_CODE
		piccard_timer->stop();
#endif

		// Check Cell crossing

		// Handle Cell crossing
/*		if(pdata->ndimensions == 1)
		{
			for(int i=0;i<N;i++)
			{
				ix(i) = ix(i) + icross(i);
				px(i) = px(i) - icross(i);

				if(icross(i) != 0)
				{
					printf("Particle %i at boundary %f, %i\n",i,px(i),ix(i));
				}
			}

		}
*/
		// Accumulate Current
		//printf("Accumulate Current %i\n",iter(0));
		//accumulate_current(pdata,currents,v_half,dtau);


//		for(int i=0;i<N;i++)
//		{
//			if(dtau(i) < 1.0e-6*pdata->dt)
//				printf("Warning dtau(%i) at %i = %e is really small\n",pid[i],iter(i),dtau(i));
//		}

		dt_finished += dtau;

		for(int i=0;i<N;i++)
			iter(i) = iter(i) + (!ifinished(i));

		// If this is a vector list, then we need to check the exit condition
		if(N > 1)
		{
			irun = !check_exit(pdata,iter,nSubcycl_max);

		}
		else
		{

			if(iter(0) >= nSubcycl_max)
				irun = 0;
			else if(dt_finished(0) >= pdata->dt)
				irun = 0;
		}

	}


}



template<const int N,const int nSpatial,const int nVel,const bool iEM> __attribute__((noinline)) __device__
void ParticleObjNT<N,nSpatial,nVel,iEM>::PicardIterate(PlasmaData* pdata, FieldData* fields,
		CurrentTally* currents,typevecN<realkind,N>& dtau0)
{

	int iter;

	typevecN<typevecN<realkind,N>,nSpatial> position_half;
	typevecN<typevecN<realkind,N>,nSpatial> position_next;
	typevecN<typevecN<realkind,N>,nSpatial> position_last;

	typevecN<typevecN<realkind,N>,nVel> velocity_half;
	typevecN<typevecN<realkind,N>,nVel> velocity_next;
	typevecN<typevecN<realkind,N>,nVel> velocity_last;


	typevecN<realkind,N> dtau = dtau0;


	typevecN<typevecN<realkind,N>,nVel> accel;

	//typevecN<float3,N> div_a;

	realkind inv_divs[3];

	if(nSpatial > 0)
		inv_divs[0] = pdata->didx;

	if(nSpatial > 1)
		inv_divs[1] = pdata->didy;

	if(nSpatial > 2)
		inv_divs[2] = pdata->didz;



	typevecN<realkind,N> residual;
//	typevecN<realkind,N> residual_d;
//	typevecN<realkind,N> residual_n;

	double avg_resid  = 10.0f;

	int niter_max = pdata -> niter_max;

/*
	x_next = (vx * inv_divs.x) * dtau + px;
	y_next = (vy * inv_divs.y) * dtau + py;
	z_next = (vz * inv_divs.z) * dtau + pz;

	x_half = (x_next + px) * 0.5f;
	y_half = (y_next + py) * 0.5f;
	z_half = (z_next + pz) * 0.5f;

	accel = get_accel(fields,x_half,y_half,z_half,
			vx,vy,vz);
*/
//	position_next = position;
//	position_last = position;


	for(int i=0;i<nSpatial;i++)
		position_next(i) = (velocity(i) * dtau) * inv_divs[i] + position(i);


	position_half = (position_next + position) * realkind(0.5);


	accel = get_accel(fields,position_half,velocity);



	for(int i=0;i<nVel;i++)
		velocity_next(i) = accel(i) * dtau + velocity(i);

	residual = 0.0;


	iter = 0;
	while(avg_resid > pdata->epsilon_r)
	{

//		for(int i=0;i<N;i++)
//			if(isnan(velocity(0)(i))||isnan(position(0)(i))||isnan(velocity_next(0)(i)))
//			printf("Stuff out(%i) = %e | %e | %e | %e | %e | %e\n",i,position_next(0)(i),position_half(0)(i),
//					velocity_next(0)(i),velocity_half(0)(i),accel(0)(i),dtau(i));


		velocity_half = (velocity_next + velocity) * realkind(0.5);

//		if(nVel == 3)
//			for(int i=0;i<nVel;i++)printf("vhalf[%i]: %e = (%e + %e) * 0.5\n",i,
//					velocity_half(i)(0),velocity_next(i)(0),velocity(i)(0));

		position_last = position_next;


		// Get the delta tau till a cell crossing
		for(int i=0;i<nSpatial;i++)
			dtau = time_till_crossing(velocity_half(i),position(i),dtau0,inv_divs[i]);


		// Dampen solution if it gets stuck
		if(iter%8 == 7)
			dtau0 = dtau * 0.5;

		for(int i=0;i<nSpatial;i++)
			position_next(i) = ((velocity_half(i) * inv_divs[i]) * dtau)  + position(i);



		position_half = (position_next + position) * realkind(0.5);


		accel = get_accel(fields,position_half,velocity_half);
		//accel = get_accel2(fields,position_half,dtau);


//		residual = abs((velocity_next(0) - velocity(0))/dtau - accel(0));

		velocity_last = velocity_next;

		for(int i=0;i<nVel;i++)
			velocity_next(i) = (accel(i) * dtau) + velocity(i);

		// Residual Calculation

//		if(nVel == 1)
//			residual =  abs(velocity_next(0) - velocity_last(0)) /
//			(abs(velocity_next(0)) + abs(velocity_last(0)));
//		else if(nVel == 2)
//			residual =  (abs(velocity_next(0) - velocity_last(0))
//					+ abs(velocity_next(1) - velocity_last(1))) /
//					(abs(velocity_next(0)) + abs(velocity_last(0))
//					+ abs(velocity_next(1)) + abs(velocity_last(1)));
//		else if(nVel == 3)
//			residual =  (abs(velocity_next(0) - velocity_last(0))
//					+ abs(velocity_next(1) - velocity_last(1))
//					+ abs(velocity_next(2) - velocity_last(2))) /
//					(abs(velocity_next(0)) + abs(velocity_last(0))
//					+ abs(velocity_next(1)) + abs(velocity_last(1))
//					+ abs(velocity_next(2)) + abs(velocity_last(2)));
//
//		if(nSpatial == 1)
//			residual +=  abs(position_next(0) - position_last(0)) /
//			(abs(position_next(0)) + abs(position_last(0)));
//		else if(nSpatial == 2)
//			residual +=  (abs(position_next(0) - position_last(0))
//					+ abs(position_next(1) - position_last(1))) /
//					(abs(position_next(0)) + abs(position_last(0))
//					+ abs(position_next(1)) + abs(position_last(1)));
//		else if(nSpatial == 3)
//			residual +=  (abs(position_next(0) - position_last(0))
//					+ abs(position_next(1) - position_last(1))
//					+ abs(position_next(2) - position_last(2))) /
//					(abs(position_next(0)) + abs(position_last(0))
//					+ abs(position_next(1)) + abs(position_last(1))
//					+ abs(position_next(2)) + abs(position_last(2)));


		avg_resid = 0;

		residual = abs(((position_next(0) - position(0)))/(dtau*inv_divs[0])
				- (velocity_next(0) + velocity(0))*0.5);

		for(int i=0;i<N;i++)
		{
			avg_resid += (residual(i))*(!ifinished(i));
		}

		avg_resid /= (realkind)N;

//		printf("avg_resid = %e at %i\n",avg_resid,iter);


		iter++;

		if(iter > niter_max)
			break;


	}

	for(int i=0;i<N;i++)
	{
		npiccard(i) += iter*(!ifinished(i));
//		npiccard2(i) += iter*iter;
	}

	velocity_half = (velocity_next + velocity) * realkind(0.5);

	for(int i=0;i<nSpatial;i++)
		position_next(i) = ((velocity_half(i) * inv_divs[i]) * dtau)  + position(i);

	position_half = (position_next + position) * realkind(0.5);
//#ifndef GPU_CODE
	accumulate_current(pdata,currents,position_half,velocity_half,dtau);
//#endif

	// Check Boundary Conditions and update positions / velocities
	typevecN<int,N> ishift;
	typevecN<int,nSpatial> dims;

	if(nSpatial > 0)
		dims(0) = pdata->nx;
	if(nSpatial > 1)
		dims(1) = pdata->ny;
	if(nSpatial > 2)
		dims(2) = pdata->nz;

	for(int i=0;i<nSpatial;i++)
	{
//		ishift = ifloor((position_next(i)-1.0e-14)/(1.0-2.0e-14));

		ishift = ifloor((position_next(i)));
		position(i) = position_next(i) - ishift;
		iposition(i) += ishift;
		iposition(i) = imod(imod(iposition(i),dims(i))+dims(i),dims(i));

//		for(int j=0;j<N;j++)
//		{
//			if((position_next(i)(j) > 1.0+1.0e-10)||(position_next(i)(j) < -1.0e-10))
//				printf("Cell Crossing Detected: %e %i %i\n",position_next(i)(j),iposition(i)(j),ishift(i));
//		}

	}

	for(int i=0;i<nVel;i++)
		velocity(i) = velocity_next(i);


	dtau0 = dtau;




}


template<int N,const int nSpatial,const int nVel,const bool iEM> __device__
typevecN<typevecN<realkind,N>,3> ParticleObjNT<N,nSpatial,nVel,iEM>::calc_div_a(
		PlasmaData* pdata, FieldData* fields,
		const typevecN<typevecN<realkind,N>,3>& accel)
{
	// Returns divergence of acceleration dotted into velocity

	/*
	 * Tensor form:
	 * [0].x [1].x [2].x
	 * [0].y [1].y [2].y
	 * [0].z [1].z [2].z
	 */
/*
	typevecN<typevecN<realkind,N>,3> divB[3];

	typevecN<typevecN<realkind,N>,3> div_vxB[3];

	//typevecN<typevecN<realkind,N>,3> Bvec = get_B<FieldData_deriv_f>(fields,px,py,pz,ix,iy,iz);


	// v x B jacobian

	typevecN<typevecN<realkind,N>,3> div_a[3];


	divB[0] = get_B<FieldData_deriv_dfdx>(fields,px,py,pz,ix,iy,iz);

	divB[1] = get_B<FieldData_deriv_dfdy>(fields,px,py,pz,ix,iy,iz);

	divB[2] = get_B<FieldData_deriv_dfdz>(fields,px,py,pz,ix,iy,iz);

	div_vxB[0] = cross_productN3(vx,vy,vz,divB[0]);
	div_vxB[1] = cross_productN3(vx,vy,vz,divB[1]);
	div_vxB[2] = cross_productN3(vx,vy,vz,divB[2]);


	div_a[0] += div_vxB[0];
	div_a[1] += div_vxB[1];
	div_a[2] += div_vxB[2];

	// Electric Field Jacobian

	typevecN<typevecN<realkind,N>,3> divE[3];

	divE[0] = get_E<FieldData_deriv_dfdx>(fields,px,py,pz,ix,iy,iz);

	divE[1] = get_E<FieldData_deriv_dfdy>(fields,px,py,pz,ix,iy,iz);

	divE[2] = get_E<FieldData_deriv_dfdz>(fields,px,py,pz,ix,iy,iz);

	// dot div_a into velocity
	typevecN<typevecN<realkind,N>,3> result;

	result = (div_a[0](0) + div_a[1](0) + div_a[2](0)) * vx;
	result = (div_a[0](1) + div_a[1](1) + div_a[2](1)) * vy;
	result = (div_a[0](2) + div_a[1](2) + div_a[2](2)) * vz;

	//for(int i=0;i<N;i++)
	//{
	//	printf("div_a = %e %e %e\n",result(i).x,result(i).y,result(i).z);
	//}

	return result;

*/

	typevecN<typevecN<realkind,N>,3> result;
	return result;

}

template<int N,const int nSpatial,const int nVel,const bool iEM> __device__
typevecN<realkind,N> ParticleObjNT<N,nSpatial,nVel,iEM>::estimate_dtau(PlasmaData* pdata, FieldData* fields)
{
#ifndef GPU_CODE
	dtau_est_timer->start();
#endif
	typevecN<realkind,N> dtau;

	typevecN<typevecN<realkind,N>,nVel> accel;
	typevecN<realkind,N> div_a;
	typevecN<realkind,N> alpha,beta,gamma2;


	dtau = pdata->dt/pdata->nSubcycle_min;

	if(nVel == 1){
	// Calculate acceleration
	accel = get_accel(fields,position,velocity);

	// Calculate div_a dot v
	//div_a = calc_div_a(pdata,fields,accel);


	for(int i=0;i<N;i++)
	{
		div_a(i) = (fields->getE(iposition(0)(i),0,0,0)-fields->getE(iposition(0)(i)-1,0,0,0))
				*fields->q2m[species]*pdata->didx;
	}
	// Using the l1norm instead of the l2norm for speed reasons
	typevecN<realkind,N> norm_diva;

	typevecN<realkind,N> norm_a;

	norm_diva = abs(div_a);
	norm_a = abs(accel(0));

	alpha = (abs(norm_diva) + abs(norm_a)) * realkind(0.5);
	beta =  ((norm_a)  + abs(velocity(0))) * (pdata->epsilon_r);
	gamma2 = pdata->epsilon_a;




	typevecN<realkind,N> d = 1.0/sqrt(alpha);


//	dtau = vmin(dtau,(pdata->dt - dt_finished));

	dtau = vmin(dtau,max(sqrt(gamma2)*d,beta*d*d));

	dtau = vmin(dtau,abs((pdata->dxdi)/velocity(0)));

	//dtau = 0.001f;
	}





	dtau = vmin(dtau,(pdata->dt - dt_finished));
	if(iEM)
	{
		typevecN<typevecN<realkind,N>,nVel> Bfield;

		for(int i=0;i<N;i++)
		{
			Bfield(0)(i)= fields->getB(iposition(0)(i),iposition(1)(i),0,0);
			Bfield(1)(i)= fields->getB(iposition(0)(i),iposition(1)(i),0,1);
			Bfield(2)(i)= fields->getB(iposition(0)(i),iposition(1)(i),0,2);
		}
		typevecN<realkind,N> B_mag = sqrt(Bfield(0)*Bfield(0) + Bfield(1)*Bfield(1) + Bfield(2)*Bfield(2));
		typevecN<realkind,N> omega_c = abs(B_mag * pdata->qspecies[species]/pdata->mspecies[species]);

		dtau = vmin(dtau,0.1f/omega_c);

		dtau = vmin(dtau,abs((pdata->dxdi)/velocity(0)));
		dtau = vmin(dtau,abs((pdata->dydi)/velocity(1)));
	}

//	for(int i=0;i<N;i++)
//		dtau(i) = fminf(dtau(i),pdata->dt/pdata->nSubcycle_min);

#ifndef GPU_CODE
	dtau_est_timer->stop();
#endif
	return dtau;

}

template<int N,const int nSpatial,const int nVel,const bool iEM> __attribute__((noinline)) __device__
typevecN<realkind,N> ParticleObjNT<N,nSpatial,nVel,iEM>::time_till_crossing(typevecN<realkind,N>& vin_half,
								typevecN<realkind,N>& pin,typevecN<realkind,N>& dtau0,const realkind scale)
{
#ifndef GPU_CODE
	crossing_timer->start();
#endif
	const realkind cell_min = -1.0e-6;
	const realkind cell_max = 1.0+1.0e-6;
	typevecN<realkind,N> dtau;
	typevecN<typevecN<realkind,N>,2> dtau_new;
	typevecN<realkind,N> cell_face;
//	typevecN<realkind,N> radical;
//	typevecN<realkind,N> sqrt_radical;

	cell_face = sgn(vin_half) * (0.5+1.0e-6);

	for(int i=0;i<N;i++)
	{
		cell_face(i) += 0.5;
	}

	dtau = (pin-cell_face)/(vin_half*(-1.0f*scale));

	for(int i=0;i<N;i++)
	{

		// Gaurd against negative and imaginary dtau's
		dtau(i) = ((dtau(i) > 0.0)) ? dtau(i):dtau0(i);

	}

	dtau = vmin(dtau0,dtau);

//	dtau_new(0) = (pin-cell_min)/(vin_half*(-scale));
//	dtau_new(1) = (pin-cell_max)/(vin_half*(-scale));
//
//	for(int i=0;i<N;i++)
//	{
//		if(!ifinished(i)){
//		// Gaurd against negative and imaginary dtau's
//		dtau_new(0)(i) = ((dtau_new(0)(i) > 0.0f)) ? dtau_new(0)(i):dtau0(i);
//		dtau_new(1)(i) = ((dtau_new(1)(i) > 0.0f)) ? dtau_new(1)(i):dtau0(i);
//
//		dtau(i) = fmin(dtau0(i),fmin((dtau_new(0))(i),(dtau_new(1))(i)));
//		}
//		else
//			dtau(i) = 0;
//	}


//	// first we get the time to the lower cell boundary
//	typevecN<realkind,N> a;
//	typevecN<realkind,N> b;
//
//	a = accel * (0.5f*scale);
//	b = vin_half * scale;
//
//
//
//	for(int i=0;i<N;i++)
//	{
//		dtau_new(0)(i) = (-b(i) + sgn(b(i))*sqrtf(b(i)*b(i) - 4.0f*a(i)*(pin(i) - cell_min)))/(2.0f*a(i));
//		dtau_new(1)(i) = (-b(i) + sgn(b(i))*sqrtf(b(i)*b(i) - 4.0f*a(i)*(pin(i) - cell_max)))/(2.0f*a(i));
//
//		// Gaurd against negative and imaginary dtau's
//		dtau_new(0)(i) = ((dtau_new(0)(i) > 0.0f)) ? dtau_new(0)(i):dtau0(i);
//		dtau_new(1)(i) = ((dtau_new(1)(i) > 0.0f)) ? dtau_new(1)(i):dtau0(i);
//
//		dtau(i) = fmin(dtau0(i),fmin((dtau_new(0))(i),(dtau_new(1))(i)));
//	}

#ifndef GPU_CODE
	crossing_timer->stop();
#endif
	return dtau;

}

template<int N,const int nSpatial,const int nVel,const bool iEM> __attribute__((noinline)) __device__
void ParticleObjNT<N,nSpatial,nVel,iEM>::accumulate_current(PlasmaData* pdata, CurrentTally* currents,
		const typevecN<typevecN<realkind,N>,nSpatial> position_half,
		const typevecN<typevecN<realkind,N>,nVel> velocity_half,
		const typevecN<realkind,N>& dtau)
{
#ifndef GPU_CODE
	tally_timer->start();
#endif
	typevecN<typevecN<realkind,N>,nSpatial> x_out;
	typevecN<typevecN<int,N>,nSpatial> ix_out;

	typevecN<int,N> ishift;


	for(int i=0;i<nSpatial;i++)
	{
		ishift = ifloor(position_half(i));
		x_out(i) = position_half(i) - ishift;
		ix_out(i) = iposition(i) + ishift;


	}

	for(int i=0;i<N;i++)
	{
		realkind scale;
		scale = dtau(i)/pdata->dt;
		if(!ifinished(i)){
		if(nSpatial == 1)
		{
			if(nVel == 1)
//				currents->tally(x_out(0)(i),0,0,velocity_half(0)(i),0,0,ix_out(0)(i),0,0,dtau(i)/pdata->dt);
				currents->tally1d1v(x_out(0)(i),velocity_half(0)(i),
								ix_out(0)(i),
								scale);
			else if (nVel == 2)
				currents->tally1d2v(x_out(0)(i),
								velocity_half(0)(i),velocity_half(1)(i),
								ix_out(0)(i),
								dtau(i)/pdata->dt);
			else if (nVel == 3)
				currents->tally1d3v(x_out(0)(i),
								velocity_half(0)(i),velocity_half(1)(i),velocity_half(2)(i),
								ix_out(0)(i),
								dtau(i)/pdata->dt);
		}
		else if(nSpatial == 2)
		{

			if (nVel == 2)
				currents->tally2d2v(x_out(0)(i),x_out(1)(i),
								velocity_half(0)(i),velocity_half(1)(i),
								ix_out(0)(i),ix_out(1)(i),
								dtau(i)/pdata->dt);
			else if (nVel == 3)
				currents->tally2d3v(x_out(0)(i),x_out(1)(i),
								velocity_half(0)(i),velocity_half(1)(i),velocity_half(2)(i),
								ix_out(0)(i),ix_out(1)(i),
								dtau(i)/pdata->dt);

		}
		else if(nSpatial == 3)
		{

			if (nVel == 3)
				currents->tally3d3v(x_out(0)(i),x_out(1)(i),x_out(2)(i),
								velocity_half(0)(i),velocity_half(1)(i),velocity_half(2)(i),
								ix_out(0)(i),ix_out(1)(i),ix_out(2)(i),
								dtau(i)/pdata->dt);

		}
		}

	}

#ifndef GPU_CODE
	tally_timer->stop();
#endif
}



template<int N,const int nSpatial,const int nVel,const bool iEM> __device__
bool ParticleObjNT<N,nSpatial,nVel,iEM>::check_exit(PlasmaData* pdata,const typevecN<int,N>& iter, const int nSubcycl_max)
{
	int result = 0;

	for(int i=0;i<N;i++)
	{
		if(dt_finished(i) >= pdata->dt*(1-1.0e-8))
		{
			dt_finished(i) = pdata->dt;
			ifinished(i) = 1;
			result += 1;
		}
		else if(iter(i) >= nSubcycl_max)
		{
			ifinished(i) = 1;
			result += 1;
			break;
		}
	}

	return result >= (N);
}



template<int N,const int nSpatial,const int nVel,const bool iEM>
template<enum FieldData_deriv ideriv> __device__
typevecN<typevecN<realkind,N>,nVel> ParticleObjNT<N,nSpatial,nVel,iEM>::get_B(FieldData* fields,
								const typevecN<typevecN<realkind,N>,nSpatial>& x,
								const typevecN<typevecN<int,N>,nSpatial>& ix)
{

	typevecN<typevecN<realkind,N>,nVel> B_out;

	for(int j=0;j<nVel;j++)
	{

		for(int i=0; i<N; i++)
		{
			if(!ifinished(i)){
			if(nSpatial == 1)
				B_out(j)(i) = fields->intrpB(x(0)(i),0,0,ix(0)(i),0,0,j,ideriv);
			else if(nSpatial == 2)
				B_out(j)(i) = fields->intrpB(x(0)(i),x(1)(i),0,ix(0)(i),ix(1)(i),0,j,ideriv);
			else if(nSpatial == 3)
				B_out(j)(i) = fields->intrpB(x(0)(i),x(1)(i),x(2)(i),ix(0)(i),ix(1)(i),ix(2)(i),j,ideriv);
			}
		}
	}

//	if(nVel == 3)
//	for(int i=0;i<N;i++)
//		printf("B_out(%i) = %e, %e, %e\n",i,B_out(0)(i),B_out(1)(i),B_out(2)(i));

	return B_out;
}


template<int N,const int nSpatial,const int nVel,const bool iEM>
template<enum FieldData_deriv ideriv> __device__
typevecN<typevecN<realkind,N>,nVel> ParticleObjNT<N,nSpatial,nVel,iEM>::get_E(FieldData* fields,
								const typevecN<typevecN<realkind,N>,nSpatial>& x,
								const typevecN<typevecN<int,N>,nSpatial>& ix)
{

	typevecN<typevecN<realkind,N>,nVel> E_out;


		for(int j=0;j<nVel;j++)
		{

			for(int i=0; i<N; i++)
			{
				if(!ifinished(i)){
				if(nSpatial == 1)
					E_out(j)(i) = fields->intrpE(x(0)(i),0,0,ix(0)(i),0,0,j,ideriv);
				else if(nSpatial == 2)
					E_out(j)(i) = fields->intrpE(x(0)(i),x(1)(i),0,ix(0)(i),ix(1)(i),0,j,ideriv);
				else if(nSpatial == 3)
					E_out(j)(i) = fields->intrpE(x(0)(i),x(1)(i),x(2)(i),ix(0)(i),ix(1)(i),ix(2)(i),j,ideriv);
				}
			}
		}

//	if(nVel == 3)
//	for(int i=0;i<N;i++)
//		printf("E_out(%i) = %e, %e, %e\n",i,E_out(0)(i),E_out(1)(i),E_out(2)(i));



	return E_out;
}

template<int N,const int nSpatial,const int nVel,const bool iEM> __device__
typevecN<typevecN<realkind,N>,nVel> ParticleObjNT<N,nSpatial,nVel,iEM>::get_accel(FieldData* fields,
										const typevecN<typevecN<realkind,N>,nSpatial>& x,
										const typevecN<typevecN<realkind,N>,nVel>& v)
{

#ifndef GPU_CODE
		accel_timer->start();
#endif

	typevecN<typevecN<realkind,N>,nSpatial> x_out;
	typevecN<typevecN<int,N>,nSpatial> ix_out;

	typevecN<typevecN<realkind,N>,nVel> accel;





	ix_out = ifloor(x);
	x_out = x - ix_out;
	ix_out += iposition;




	if((nSpatial == 2)&&(nVel == 3))
	{
		for(int i=0;i<N;i++)
		{
			if(!ifinished(i))
			{
				fields->intrpAccel(x_out(0)(i),x_out(1)(i),0,
								v(0)(i),v(1)(i),v(2)(i),
								ix_out(0)(i),ix_out(1)(i),0,
								accel(0)(i),accel(1)(i),accel(2)(i));
			}
		}
	}
	else
	{

		typevecN<typevecN<realkind,N>,nVel> E;

		E = get_E<FieldData_deriv_f>(fields,x_out,ix_out);

		if(iEM)
		{
			if(nVel == 3)
			{
				//main formulation
				typevecN<typevecN<realkind,N>,nVel> B;
				B = get_B<FieldData_deriv_f>(fields,x_out,ix_out);
				accel = cross_productNV(v(0),v(1),v(2),B);

				for(int i=0;i<nVel;i++)
					accel(i) = accel(i) + E(i);


			}
			else
				accel = E;
		}
		else
			accel = E;
	}


	accel = accel * (-(fields->q2m[species]));

//	if(nVel == 3)
//	for(int i=0;i<N;i++)
//		printf("Accel(%i) = %e, %e, %e\n",i,accel(0)(i),accel(1)(i),accel(2)(i));

#ifndef GPU_CODE
		accel_timer->stop();
#endif

	return accel;

}

template<int N,const int nSpatial,const int nVel,const bool iEM> __device__
typevecN<typevecN<realkind,N>,nVel> ParticleObjNT<N,nSpatial,nVel,iEM>::get_accel2(FieldData* fields,
										const typevecN<typevecN<realkind,N>,nSpatial>& x,
										typevecN<realkind,N> dtau)
{

#ifndef GPU_CODE
		accel_timer->start();
#endif

	typevecN<typevecN<realkind,N>,nSpatial> x_out;
	typevecN<typevecN<int,N>,nSpatial> ix_out;

	typevecN<typevecN<realkind,N>,nVel> accel;

	typevecN<typevecN<realkind,N>,nVel> E;


	typevecN<typevecN<realkind,N>,nVel> vtil;
	typevecN<typevecN<realkind,N>,nVel> vbar;

	realkind q2mt = -fields->q2m[species];


	ix_out = ifloor(x);
	x_out = x - ix_out;
	ix_out += iposition;



	E = get_E<FieldData_deriv_f>(fields,x_out,ix_out);

	for(int i=0;i<nVel;i++)
	vtil(i) = velocity(i) + E(i)*dtau*(0.5*q2mt);

	vbar = vtil;

	if(iEM)
	{
		if(nVel == 3)
		{
			//main formulation
			typevecN<typevecN<realkind,N>,nVel> B;
			typevecN<typevecN<realkind,N>,nVel> vtilxB;
			typevecN<realkind,N> vdotB;
			typevecN<realkind,N> Bmag;
			typevecN<realkind,N> dtaut;

			dtaut = dtau*(realkind(0.5)*q2mt);

			B = get_B<FieldData_deriv_f>(fields,x_out,ix_out);

			vtilxB = cross_productNV(vtil(0),vtil(1),vtil(2),B);

			Bmag = B(0)*B(0);
			for(int i=1;i<nVel;i++)
				Bmag = Bmag + B(i)*B(i);


			for(int i=0;i<N;i++)
			Bmag(i) = (Bmag(i)*dtaut(i)*dtaut(i))+1.0;

			vdotB = vtil(0)*B(0);

			for(int i=1;i<nVel;i++)
				vdotB = vdotB + vtil(i)*B(i);

			for(int i=0;i<nVel;i++)
			{
				vbar(i) = vtil(i) + (vtilxB(i) + vdotB*B(i)*dtaut)*dtaut/Bmag;
			}



		}
		else
			accel = E;
	}
	else
		accel = E;


	for(int i=0;i<nVel;i++)
		accel(i) = ((vbar(i)-velocity(i))/dtau)*2.0;

	for(int j=0;j<nVel;j++)
		for(int i=0;i<N;i++)
			accel(j)(i) *= (!ifinished(i));

//	if(nVel == 3)
//	for(int i=0;i<N;i++)
//		printf("Accel(%i) = %e, %e, %e\n",i,accel(0)(i),accel(1)(i),accel(2)(i));

#ifndef GPU_CODE
		accel_timer->stop();
#endif

	return accel;

}




/*

ploop<VEC_LENGTH_MAX> dummy_ploop;

ParticleObjN<1> dummy_p;

ploop<VEC_LENGTH_MAX> DummyFUNCTION(PlasmaData* pdata,FieldData* fields,
		CurrentTally* currents,typevecN<int,VEC_LENGTH_MAX>& iter,const int nSubcycl_max)
{
	ploop<VEC_LENGTH_MAX> ploop_out;
	ploop_out.particle.push(pdata,fields,currents,iter,nSubcycl_max);

	return ploop_out;
}
*/

 /* PARTICLEOBJNT_INL */







