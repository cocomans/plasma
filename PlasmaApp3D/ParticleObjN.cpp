
#include "ParticleObjN.h"
#include "FieldData.h"
#include "PlasmaData.h"
#include "vec_funcs.h"
#include "typevecN.h"
#include "CurrentTally.h"
#include "ShapeFunctions.h"



template<int ileft> __attribute__ ((noinline))
void PushN(PlasmaData* 		pdata,
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
	ParticleObjN<ileft> particles2(*iptcl);
	typevecN<int,ileft> iter2;

	iter2 = *iter_array;

	// Copy from main array
	particles2 = *plist;


	// Push
	particles2.push(pdata,fields,current,iter2,pdata->nSubcycle_max);

	// Write back to main array
	//printf("writing particles back\n");

	//*plist = particles2;
	for(int j=0;j<ileft;j++)
	{
		//printf("iptcl[%i] = %i\n",j,iptcl[0][j]);
		particles2.write_back(*plist,j);
	}

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
			//printf("warning particle finished before time step was finished\n");
			idone = 1;
		}
		else
			idone = 0;


		if(idone)
		{
			nSubSteps_done += iter2(j);
			iter2(j) = 0;
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

template<const int ileft> __attribute__ ((noinline))
void shrink_push(PlasmaData* 		pdata,
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
		PushN<ileft>(pdata,fields,current,plist,
							iter_array,iptcl,iptcl_new,
							nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);

	} /* if(ileft == nptcls_left) */
	else if(ileft > 0)
	{
		shrink_push<ileft-1>(pdata,fields,current,plist,
							iter_array,iptcl,iptcl_new,
							nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);
	}
}


}

template<> __attribute__ ((noinline))
void shrink_push<1>(PlasmaData* 		pdata,
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
	if(nptcls_left == 1)
	{
		PushN<1>(pdata,fields,current,plist,
							iter_array,iptcl,iptcl_new,
							nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);

	} /* if(1 == nptcls_left) */

}

template<> __attribute__ ((noinline))
void shrink_push<VEC_LENGTH_MAX>(PlasmaData* 		pdata,
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

	if(nptcls_left == VEC_LENGTH_MAX)
	{
		PushN<VEC_LENGTH_MAX>(pdata,fields,current,plist,
							iter_array,iptcl,iptcl_new,
							nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);

	} /* if(ileft == nptcls_left) */
	else
	{
		shrink_push<VEC_LENGTH_MAX-1>(pdata,fields,current,plist,
							iter_array,iptcl,iptcl_new,
							nptcls_left,nptcl_done,nptcls_process,nSubSteps_done);
	}


}




template<int N> __attribute__((noinline))
void ParticleObjN<N>::push(PlasmaData* pdata, FieldData* fields,
		CurrentTally* currents, typevecN<int,N>& iter, const int nSubcycl_max)
{

	bool irun = 1;

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
		PicardIterate(pdata,fields,currents,dtau);


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

		iter++;

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



template<int N> __attribute__((noinline))
void ParticleObjN<N>::PicardIterate(PlasmaData* pdata, FieldData* fields,
		CurrentTally* currents,typevecN<realkind,N>& dtau0)
{

	int iter;

	typevecN<realkind,N> x_half;
	typevecN<realkind,N> y_half;
	typevecN<realkind,N> z_half;

	typevecN<realkind,N> vx_half;
	typevecN<realkind,N> vy_half;
	typevecN<realkind,N> vz_half;


	typevecN<realkind,N> x_next;
	typevecN<realkind,N> y_next;
	typevecN<realkind,N> z_next;
	typevecN<realkind,N> x_last;
	typevecN<realkind,N> y_last;
	typevecN<realkind,N> z_last;
	typevecN<realkind,N> vx_next;
	typevecN<realkind,N> vy_next;
	typevecN<realkind,N> vz_next;
	typevecN<realkind,N> vx_last;
	typevecN<realkind,N> vy_last;
	typevecN<realkind,N> vz_last;

	typevecN<realkind,N> dtau = dtau0;




	//typevecN<float3,N> v_last = make_float3N(vx,vy,vz);
	typevecN<typevecN<realkind,N>,3> accel;

	//typevecN<float3,N> div_a;

	float3 inv_divs = make_float3(pdata->didx,pdata->didy,pdata->didz);



	typevecN<realkind,N> residual;

	double avg_resid  = 10.0f;

	int niter_max = pdata -> niter_max;

/*
	accel = get_accel(fields,x_old,v_old);

	div_a = calc_div_a(pdata,fields,accel);

	v_next += (accel + div_a * dtau * 0.5f)* dtau;

	v_half = (v_next + v_old) * 0.5f;

	x_next = (v_half * inv_divs) * dtau + x_old;

	accel = get_accel(fields,x_half,v_half);

	v_next = accel * dtau + v_old;

	v_half = (v_next + v_old) * 0.5f;

	residual = l1normN((x_next - x_old)/dtau - v_half*inv_divs)/l1normN(v_old*inv_divs);

	avg_resid = 0;

	for(int i=0;i<N;i++)
	{
		avg_resid += (residual(i));
	}
	printf("residual = %e\n",avg_resid/(2.0 * N));
	//v_next = v_old;
*/

/*
	dtau_m = (px-cell_min)/(vx*(-1.0f*pdata->didx));
	dtau_p = (px-cell_max)/(vx*(-1.0f*pdata->didx));

	for(int i=0;i<N;i++)
	{
		// Gaurd against negative and imaginary dtau's
		dtau_p(i) = ((dtau_p(i) > 0.0f)) ? dtau_p(i):dtau0(i);
		dtau_m(i) = ((dtau_m(i) > 0.0f)) ? dtau_m(i):dtau0(i);

		dtau(i) = fmin(dtau(i),fmin(dtau_p(i),dtau_m(i)));
	}
*/
	x_next = (vx * inv_divs.x) * dtau + px;
	y_next = (vy * inv_divs.y) * dtau + py;
	z_next = (vz * inv_divs.z) * dtau + pz;

//	for(int i=0;i<N;i++)
//	{
//		x_next(i) = fmin(fmax(x_next(i),cell_min),cell_max);
//	}


	x_half = (x_next + px) * 0.5f;
	y_half = (y_next + py) * 0.5f;
	z_half = (z_next + pz) * 0.5f;

	accel = get_accel(fields,x_half,y_half,z_half,
			vx,vy,vz);




/*
	// first we get the time to the lower cell boundary
	radical = (vx*vx*(pdata->didx*pdata->didx) - accel(0)*(px - cell_min)*(2.0f*pdata->didx));
	dtau_p = (vx*(-pdata->didx) + sqrt(radical))/(accel(0) * (pdata->didx));
	dtau_m = (vx*(-pdata->didx) - sqrt(radical))/(accel(0) * (pdata->didx));

	for(int i=0;i<N;i++)
	{
		// Gaurd against negative and imaginary dtau's
		dtau_p(i) = ((dtau_p(i) > 0.0f)&&(radical(i) >= 0.0f)) ? dtau_p(i):dtau0(i);
		dtau_m(i) = ((dtau_m(i) > 0.0f)&&(radical(i) >= 0.0f)) ? dtau_m(i):dtau0(i);

		dtau(i) = fmin(dtau0(i),fmin(dtau_p(i),dtau_m(i)));
	}

	// now we get the time to the upper cell boundary
	radical = (vx*vx*(pdata->didx*pdata->didx) - accel(0)*(px - cell_max)*(2.0f*pdata->didx));
	dtau_p = (vx*(-pdata->didx) + sqrt(radical))/(accel(0) * (pdata->didx));
	dtau_m = (vx*(-pdata->didx) - sqrt(radical))/(accel(0) * (pdata->didx));

	for(int i=0;i<N;i++)
	{
		// Gaurd against negative and imaginary dtau's
		dtau_p(i) = ((dtau_p(i) > 0.0f)&&(radical(i) >= 0.0f)) ? dtau_p(i):dtau0(i);
		dtau_m(i) = ((dtau_m(i) > 0.0f)&&(radical(i) >= 0.0f)) ? dtau_m(i):dtau0(i);

		dtau(i) = fmin(dtau0(i),fmin(dtau_p(i),dtau_m(i)));
	}
*/
	//v_last = v_next;

//	dtau = time_till_crossing(vx,px,accel(0),dtau0,pdata->didx);
//	if(pdata->ndimensions == 3)
//	{
//		dtau = time_till_crossing(vy,py,accel(1),dtau,pdata->didy);
//		dtau = time_till_crossing(vz,pz,accel(2),dtau,pdata->didz);
//	}

	vx_next = (accel(0)) * dtau + vx;
	vy_next = (accel(1)) * dtau + vy;
	vz_next = (accel(2)) * dtau + vz;

	residual = 0.0f;


	iter = 0;
	while(avg_resid > pdata->epsilon_r)
	{

//		for(int i=0;i<N;i++)
//			printf("Stuff out(%i) = %e | %e | %e | %e | %e | %e\n",i,x_next(i),x_half(i),
//					vx_next(i),vx_half(i),accel(0)(i),dtau(i));

		vx_half = (vx_next + vx) * 0.5f;
		vy_half = (vy_next + vy) * 0.5f;
		vz_half = (vz_next + vz) * 0.5f;


		x_last = x_next;
		y_last = y_next;
		z_last = z_next;

		dtau = time_till_crossing(vx_half,px,accel(0),dtau,pdata->didx);
		if(pdata->ndimensions == 3)
		{
			dtau = time_till_crossing(vy,py,accel(1),dtau,pdata->didy);
			dtau = time_till_crossing(vz,pz,accel(2),dtau,pdata->didz);
		}


		if(iter%6 == 5)
			dtau0 = dtau0 * 0.5;

		x_next = (vx_half * inv_divs.x) * dtau + px;
		y_next = (vy_half * inv_divs.y) * dtau + py;
		z_next = (vz_half * inv_divs.z) * dtau + pz;




		x_half = (x_next + px) * 0.5f;
		y_half = (y_next + py) * 0.5f;
		z_half = (z_next + pz) * 0.5f;


		accel = get_accel(fields,x_half,y_half,z_half,
				vx_half,vy_half,vz_half);




		vx_last = vx_next;
		vy_last = vy_next;
		vz_last = vz_next;

		vx_next = (accel(0)) * dtau + vx;
		vy_next = (accel(1)) * dtau + vy;
		vz_next = (accel(2)) * dtau + vz;


		residual = (abs(vx_next - vx_last))/
					(abs(vx_next)+abs(vx_last));

		residual += (abs(x_next - x_last))/
					(abs(x_next)+abs(x_last));

//		residual = (abs(vx_next - vx_last)+abs(vy_next-vy_last)+abs(vz_next-vz_last))/
//					(abs(vx_next)+abs(vy_next)+abs(vz_next));
//
//		residual += (abs(x_next - x_last)+abs(y_next-y_last)+abs(z_next-z_last))/
//					(abs(x_next)+abs(y_next)+abs(z_next));

//		for(int i=0;i<N;i++)
//		{
//
//			residual(i) = fabs(((x_next(i)-px(i)))/(dtau(i)*pdata->didx)
//					- (vx_next(i) + vx(i))*(0.5f));
//
//			residual(i) += fabs(((y_next(i)-py(i)))/(dtau(i)*pdata->didy)
//					- (vy_next(i) + vy(i))*(0.5f));
//
//			residual(i) += fabs(((z_next(i)-pz(i)))/(dtau(i)*pdata->didz)
//					- (vz_next(i) + vz(i))*(0.5f));
//
//		}
//
//		residual += abs(((vx_last-vx))/(dtau)
//				- accel(0));
//
//		residual += abs(((vy_last-vy))/(dtau)
//				- accel(1));
//
//		residual += abs(((vz_last-vz))/(dtau)
//				- accel(2));

		avg_resid = 0;

		for(int i=0;i<N;i++)
		{
			avg_resid += (residual(i));
			//printf("residual = %e\n",residual(i));
		}

		avg_resid /= (realkind)N;

//		printf("avg_resid = %e at %i\n",avg_resid,iter);
		iter++;

		if(iter > niter_max)
			break;




	}

	//	printf("residual(%i) = %e\n",iter,avg_resid);



	// Check Boundary Conditions

	vx_half = (vx_next + vx) * 0.5f;
	vy_half = (vy_next + vy) * 0.5f;
	vz_half = (vz_next + vz) * 0.5f;

	x_next = (vx_half * inv_divs.x) * dtau + px;
	y_next = (vy_half * inv_divs.y) * dtau + py;
	z_next = (vz_half * inv_divs.z) * dtau + pz;

	x_half = (x_next + px) * 0.5f;
	y_half = (y_next + py) * 0.5f;
	z_half = (z_next + pz) * 0.5f;



//	for(int i=0;i<N;i++)
//	{
//
////		residual(i) = fabs(((x_next(i)-px(i)))/(dtau(i)*pdata->didx)
////				- (vx_next(i) + vx(i))*(0.5f));
////
////		residual(i) += fabs(((y_next(i)-py(i)))/(dtau(i)*pdata->didy)
////				- (vy_next(i) + vy(i))*(0.5f));
////
////		residual(i) += fabs(((z_next(i)-pz(i)))/(dtau(i)*pdata->didz)
////				- (vz_next(i) + vz(i))*(0.5f));
//
//		residual(i) = fabs(((vx_next(i)-vx(i)))/(dtau(i))
//				- accel(0)(i));
//
//		residual(i) += fabs(((vy_next(i)-vy(i)))/(dtau(i))
//				- accel(1)(i));
//
//		residual(i) += fabs(((vz_next(i)-vz(i)))/(dtau(i))
//				- accel(2)(i));
//
//	}

//	residual = abs(((vx_next-vx))
//			- accel(0)*(dtau));
//
//	residual += abs(((vy_next-vy))
//			- accel(1)*(dtau));
//
//	residual += abs(((vz_next-vz))
//			- accel(2)*(dtau));
//
//
//	residual = residual/(dtau*3.0f);
//
//	avg_resid = 0;
//
//	for(int i=0;i<N;i++)
//	{
//		avg_resid += (residual(i));
////		printf("residual(%i) = %e\n",i,residual(i));
//	}
//
//	avg_resid /= 3.0*(realkind)N;
//
//
//
//
//
//	for(int i=0;i<N;i++)
//	{
//		//if(residual(i) > 1.0e-3){
//		//printf("residual(%i) = %e at %i iterations dtau = %e\n",i,residual(i),iter,dtau(i));
//		if((x_next(i) <= -1.0e-9)||(x_next(i) >= (1.0+1.0e-9)))
//			printf("x(%i) clost to crossing, x_old = %e, x_next = %e, resid = %e\n",i,px(i),x_next(i),residual(i));
//		if((y_next(i) <= 0.0f)||(y_next(i) >= 1.0f))
//			printf("x(%i) clost to crossing, y_old = %e, y_next = %e, resid = %e\n",i,py(i),y_next(i),residual(i));
//		if((z_next(i) <= 0.0f)||(z_next(i) >= 1.0f))
//			printf("x(%i) clost to crossing, z_old = %e, z_next = %e, resid = %e\n",i,pz(i),z_next(i),residual(i));
//
//		//}
//	}

	accumulate_current(pdata,currents,x_half,y_half,z_half,
			vx_half,vy_half,vz_half,dtau);

	typevecN<int,N> ishift;

	ishift = ifloor(x_next);
	px = x_next - ishift;
	ix += ishift;

	ix = imod(imod(ix,pdata->nx)+pdata->nx,pdata->nx);

	vx = vx_next;


	if(pdata->ndimensions == 3)
	{
		ishift = ifloor(y_next);
		py = y_next - ishift;
		iy += ishift;

		iy = imod(imod(iy,pdata->ny)+pdata->ny,pdata->ny);

		vy = vy_next;

		ishift = ifloor(z_next);
		pz = z_next - ishift;
		iz += ishift;

		iz = imod(imod(iz,pdata->nz)+pdata->nz,pdata->nz);

		vz = vz_next;
	}

	dtau0 = dtau;




/*
	ishift = ifloor(x_next);
	new_ptcl.px = x_next - ishift;
	new_ptcl.ix = ix + ishift;

	ishift = ifloor(y_next);
	new_ptcl.py = y_next - ishift;
	new_ptcl.iy = iy + ishift;

	ishift = ifloor(z_next);
	new_ptcl.pz = z_next - ishift;
	new_ptcl.iz = iz + ishift;

	new_ptcl.ix = imod(imod(new_ptcl.ix,pdata->nx)+pdata->nx,pdata->nx);
	new_ptcl.iy = imod(imod(new_ptcl.iy,pdata->ny)+pdata->ny,pdata->ny);
	new_ptcl.iz = imod(imod(new_ptcl.iz,pdata->nz)+pdata->nz,pdata->nz);
*/
//	for(int i=0;i<N;i++)
//	{
//		int ishift;
//
//		ishift = floorf(x_next(i));
//		new_ptcl.px(i) = x_next(i) - ishift;
//		new_ptcl.ix(i) = ix(i) + ishift;
//
//		ishift = floorf(y_next(i));
//		new_ptcl.py(i) = y_next(i) - ishift;
//		new_ptcl.iy(i) = iy(i) + ishift;
//
//		ishift = floorf(z_next(i));
//		new_ptcl.pz(i) = z_next(i) - ishift;
//		new_ptcl.iz(i) = iz(i) + ishift;

//		new_ptcl.ix(i) = ((new_ptcl.ix(i)%pdata->nx)+pdata->nx)%pdata->nx;
	//	new_ptcl.iy(i) = ((new_ptcl.iy(i)%pdata->ny)+pdata->ny)%pdata->ny;
	//	new_ptcl.iz(i) = ((new_ptcl.iz(i)%pdata->nz)+pdata->nz)%pdata->nz;


/*
		if((new_ptcl.px(i) < 0.0f)||(new_ptcl.px(i) > 1.0f)
			||(new_ptcl.py(i) < 0.0f)||(new_ptcl.py(i) > 1.0f)
			||(new_ptcl.pz(i) < 0.0f)||(new_ptcl.pz(i) > 1.0f)
			||(new_ptcl.ix(i) >= pdata->nx)||(new_ptcl.iy(i) >= pdata->ny))
		{
			printf("Warning position out of bounds\n");
		}



	//}

	new_ptcl.vx = vx_next;
	new_ptcl.vy = vy_next;
	new_ptcl.vz = vz_next;

	new_ptcl.species = species;
	new_ptcl.pid = pid;
	new_ptcl.cluster_id = cluster_id;
	new_ptcl.dt_finished = dt_finished;

	*/

	//return new_ptcl;

}


template<int N>
typevecN<typevecN<realkind,N>,3> ParticleObjN<N>::calc_div_a(PlasmaData* pdata, FieldData* fields,
		const typevecN<typevecN<realkind,N>,3>& accel)
{
	// Returns divergence of acceleration dotted into velocity

	/*
	 * Tensor form:
	 * [0].x [1].x [2].x
	 * [0].y [1].y [2].y
	 * [0].z [1].z [2].z
	 */

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




}

template<int N>
typevecN<realkind,N> ParticleObjN<N>::estimate_dtau(PlasmaData* pdata, FieldData* fields)
{
	typevecN<realkind,N> dtau;

	typevecN<typevecN<realkind,N>,3> accel;
	//typevecN<typevecN<realkind,N>,3> div_a;
	typevecN<realkind,N> div_a;
	typevecN<realkind,N> alpha,beta,gamma2;




	if(false){
	// Calculate acceleration

/*
	// Calculate div_a dot v
	div_a = calc_div_a(pdata,fields,accel);
*/

	for(int i=0;i<N;i++)
	{
		div_a(i) = (fields->getE(ix(i),iy(i),iz(i),0)-fields->getE(ix(i)-1,iy(i),iz(i),0))
				*fields->q2m[species]*pdata->didx*vx(i);
	}
	// Using the l1norm instead of the l2norm for speed reasons
	typevecN<realkind,N> norm_diva;

	typevecN<realkind,N> norm_a;// = l1normN(accel);

	norm_diva = abs(div_a);
	norm_a = abs(accel(0));

	alpha = (abs(norm_diva) + abs(norm_a)) * 0.5f;
	beta =  ((norm_a)  + abs(vx)) * (pdata->epsilon_r);
	gamma2 = pdata->epsilon_a;




	typevecN<realkind,N> d = 1.0f/sqrt(alpha);

	dtau = pdata->dt/pdata->nSubcycle_min;

	dtau = vmin(dtau,(pdata->dt - dt_finished));

	dtau = vmin(dtau,max(sqrt(gamma2)*d,beta*d*d));

	//dtau = 0.001f;


	}
	dtau = pdata->dt/pdata->nSubcycle_min;

	dtau = vmin(dtau,(pdata->dt - dt_finished));
	if(pdata->Bmag_avg > 0.0f)
	{
		typevecN<typevecN<realkind,N>,3> Bfield = get_B<FieldData_deriv_f>(fields,px,py,pz,ix,iy,iz);
		typevecN<realkind,N> B_mag = sqrt(Bfield(0)*Bfield(0) + Bfield(1)*Bfield(1) + Bfield(2)*Bfield(2));
		typevecN<realkind,N> omega_c = B_mag * pdata->qspecies[species]/pdata->mspecies[species];

		dtau = vmin(dtau,0.1f/omega_c);
	}






//	accel = get_accel(fields,px,py,pz,vx,vy,vz);
//
//	dtau = time_till_crossing(vx,px,accel(0),dtau,pdata->didx);
//if(pdata->ndimensions == 3)
//{
//	dtau = time_till_crossing(vy,py,accel(1),dtau,pdata->didy);
//	dtau = time_till_crossing(vz,pz,accel(2),dtau,pdata->didz);
//}



	//dtau = min(dtau,abs((pdata->dxdi+pdata->dydi+pdata->dzdi)/(abs(vx)+abs(vy)+abs(vz)))*0.25f);



	//dtau = min(dtau,(pdata->dt - dt_finished));



/*
	// If this is 1D electrostatic problem then we have an analytic solution to the time of cell crossing
	if(pdata->ndimensions == 1)
	{
		realkind A,B,C,F;
		realkind omega,x0,E0,E1,alpha,q,dt1,dt2,v0;
		for(int i=0;i<N;i++)
		{
			// First we need to calculate the Electric field strength
			E0 = fields->getE(ix(i),iy(i),iz(i),0)*pdata->didx;
			E1 = fields->getE(ix(i)+1,iy(i),iz(i),0)*pdata->didx;
			alpha = E1-E0;
			q = pdata->qspecies[species];
			x0 = px(i);
			v0 = vx(i) * pdata->didx;

			omega = sqrt(abs(q*alpha/pdata->mspecies[species]));

			// The formula used depends on the signs of the charge and electric field slope
			if(sgn(q) == sgn(alpha))
			{
				B = x0/2.0f + E0/(2.0f*alpha) - v0/(2.0f*omega);
				C = v0/(2.0f*omega) + x0/2.0f + E0/(2.0f*alpha);
				// first check the time till px=0
				A = E0/alpha;
				F = A*A - 4.0f*B*C;
				if(F >= 0.0f)
				{
					F = sqrt(F);
					dt1 = log((A+F)/B)/omega;
					dt2 = log((A-F)/B)/omega;

					// guard against negative dt
					dt1 = dt1 > 0.0f ? dt1:dtau(i);
					dt2 = dt2 > 0.0f ? dt2:dtau(i);

					dt1 = fmin(dt1,dt2);
					if(dt1 < dtau(i))
					{
						dtau(i) = dt1;
						icross(i) = -1;
					}
				}

				// check the time till px=1
				A = E0/alpha+1.0f;
				F = A*A - 4.0f*B*C;
				if(F >= 0.0f)
				{
					F = sqrt(F);
					dt1 = log((A+F)/B)/omega;
					dt2 = log((A-F)/B)/omega;

					// guard against negative dt
					dt1 = dt1 > 0.0f ? dt1:dtau(i);
					dt2 = dt2 > 0.0f ? dt2:dtau(i);

					dt1 = fmin(dt1,dt2);
					if(dt1 < dtau(i))
					{
						dtau(i) = dt1;
						icross(i) = 1;
					}
				}


			}
			else
			{
				B = x0+E0/alpha;
				C = v0/omega;
				// first check the time till px=0
				A = E0/alpha;
				F = C*C*(B*B+C*C-A*A);
				if(F >= 0.0f)
				{
					F = sqrt(F);
					dt1 = acos((A*B+F)/(B*B+C*C))/omega;
					dt2 = acos((A*B-F)/(B*B+C*C))/omega;

					// guard against negative dt
					dt1 = dt1 > 0.0f ? dt1:dtau(i);
					dt2 = dt2 > 0.0f ? dt2:dtau(i);

					dt1 = fmin(dt1,dt2);
					if(dt1 < dtau(i))
					{
						dtau(i) = dt1;
						icross(i) = -1;
					}
				}

				// check the time till px=1
				A = E0/alpha+1.0f;
				F = A*A - 4.0f*B*C;
				if(F >= 0.0f)
				{
					F = sqrt(F);
					dt1 = log((A+F)/B)/omega;
					dt2 = log((A-F)/B)/omega;

					// guard against negative dt
					dt1 = dt1 > 0.0f ? dt1:dtau(i);
					dt2 = dt2 > 0.0f ? dt2:dtau(i);

					dt1 = fmin(dt1,dt2);
					if(dt1 < dtau(i))
					{
						dtau(i) = dt1;
						icross(i) = 1;
					}
				}


			}

		}


	}
*/
	//for(int i=0;i<N;i++)
//	{
		//if(isnan(dtau(i))||(dtau(i) <= 0.0f))
		//	printf("dtau(%i) = %e\n",i,dtau(i));
//	}


	return dtau;

}

template<int N> __attribute__((noinline))
typevecN<realkind,N> ParticleObjN<N>::time_till_crossing(typevecN<realkind,N>& vin_half,typevecN<realkind,N>& pin,
				typevecN<realkind,N>& accel, typevecN<realkind,N>& dtau0,const realkind scale)
{
	const realkind cell_min = -1.0e-6;
	const realkind cell_max = 1.0+1.0e-6;
	typevecN<realkind,N> dtau;
	typevecN<typevecN<realkind,N>,2> dtau_new;
	//typevecN<realkind,N> cell_face;
//	typevecN<realkind,N> radical;
//	typevecN<realkind,N> sqrt_radical;

//	cell_face = sgn(vin_half) * (0.5+1.0e-8);
//
//	for(int i=0;i<N;i++)
//	{
//		cell_face(i) += 0.5;
//	}
//
//	dtau = (pin-cell_face)/(vin_half*(-1.0f*scale));
//
//	for(int i=0;i<N;i++)
//	{
//		// Gaurd against negative and imaginary dtau's
//		dtau(i) = ((dtau(i) > 0.0f)) ? dtau(i):dtau0(i);
//
//
//		dtau(i) = fmin(dtau0(i),dtau(i));
//	}

	dtau_new(0) = (pin-cell_min)/(vin_half*(-1.0f*scale));
	dtau_new(1) = (pin-cell_max)/(vin_half*(-1.0f*scale));

	for(int i=0;i<N;i++)
	{
		// Gaurd against negative and imaginary dtau's
		dtau_new(0)(i) = ((dtau_new(0)(i) > 0.0f)) ? dtau_new(0)(i):dtau0(i);
		dtau_new(1)(i) = ((dtau_new(1)(i) > 0.0f)) ? dtau_new(1)(i):dtau0(i);

		dtau(i) = fmin(dtau0(i),fmin((dtau_new(0))(i),(dtau_new(1))(i)));
	}


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

	return dtau;

}

template<int N> __attribute__((noinline))
void ParticleObjN<N>::accumulate_current(PlasmaData* pdata, CurrentTally* currents,
		const typevecN<realkind,N> x_half,const typevecN<realkind,N> y_half,const typevecN<realkind,N> z_half,
		const typevecN<realkind,N> vx_half,const typevecN<realkind,N> vy_half,const typevecN<realkind,N> vz_half,
		const typevecN<realkind,N>& dtau)
{

	typevecN<realkind,N> x_out,y_out,z_out;
	typevecN<int,N> ix_out,iy_out,iz_out;

	ix_out = ifloor(x_half);
	x_out = x_half - ix_out;
	ix_out += ix;

	iy_out = ifloor(y_half);
	y_out = y_half - iy_out;
	iy_out += iy;

	iz_out = ifloor(z_half);
	z_out = z_half - iz_out;
	iz_out += iz;

	for(int i=0;i<N;i++)
	{
		currents->tally(x_out(i),y_out(i),z_out(i),
						vx_half(i),vy_half(i),vz_half(i),
						ix_out(i),iy_out(i),iz_out(i),
						dtau(i)/pdata->dt);
	}
}



template<int N>
bool ParticleObjN<N>::check_exit(PlasmaData* pdata,const typevecN<int,N>& iter, const int nSubcycl_max)
{
	bool result = 0;

	for(int i=0;i<N;i++)
	{
		if(dt_finished(i) >= pdata->dt)
		{
			result = 1;
			break;
		}
		else if(iter(i) >= nSubcycl_max)
		{
			result = 1;
			break;
		}
	}

	return result;
}

template<int N>
template<enum FieldData_deriv ideriv>
typevecN<typevecN<realkind,N>,3> ParticleObjN<N>::get_B(FieldData* fields, const typevecN<float3,N>& x)
{
	typevecN<typevecN<realkind,N>,3> B_out;

	for(int i=0; i<N; i++)
	{
		float3 pos;
		int3 index;
		int ishift;

		ishift = floorf(x(i).x);
		pos.x = x(i).x - ishift;
		index.x = ix(i) + ishift;

		ishift = floorf(x(i).y);
		pos.y = x(i).y - ishift;
		index.y = iy(i) + ishift;

		ishift = floorf(x(i).z);
		pos.z = x(i).z - ishift;
		index.z = iz(i) + ishift;

		//index.x = ((index.x%fields->nx)+fields->nx)%fields->nx;
		//index.y = ((index.y%fields->ny)+fields->ny)%fields->ny;
		//index.z = ((index.z%fields->nz)+fields->nz)%fields->nz;



		(B_out(0))(i) = fields->intrpB(pos.x,pos.y,pos.z,
											index.x,index.y,index.z,0,ideriv);
		(B_out(1))(i) = fields->intrpB(pos.x,pos.y,pos.z,
											index.x,index.y,index.z,1,ideriv);
		(B_out(2))(i) = fields->intrpB(pos.x,pos.y,pos.z,
											index.x,index.y,index.z,2,ideriv);

	}

	return B_out;
}

template<int N>
template<enum FieldData_deriv ideriv>
typevecN<typevecN<realkind,N>,3> ParticleObjN<N>::get_B(FieldData* fields,
								const typevecN<realkind,N>& x,
								const typevecN<realkind,N>& y,
								const typevecN<realkind,N>& z,
								const typevecN<int,N>& ix_out,
								const typevecN<int,N>& iy_out,
								const typevecN<int,N>& iz_out)
{

	typevecN<typevecN<realkind,N>,3> B_out;

	for(int j=0;j<3;j++)
	{

		for(int i=0; i<N; i++)
		{
			(B_out(j))(i) = fields->intrpB(x(i),y(i),z(i),
												ix_out(i),iy_out(i),iz_out(i),j,ideriv);
		}
	}

	//for(int i=0;i<N;i++)
	//	printf("B_out(%i) = %e, %e, %e\n",i,B_out(0)(i),B_out(1)(i),B_out(2)(i));

	return B_out;
}


template<int N>
template<enum FieldData_deriv ideriv>
typevecN<typevecN<realkind,N>,3> ParticleObjN<N>::get_E(FieldData* fields, const typevecN<float3,N>& x)
{
	typevecN<typevecN<realkind,N>,3> E_out;

	for(int i=0; i<N; i++)
	{
		float3 pos;
		int3 index;
		int ishift;

		//ishift = floorf(x(i).x);
		//pos.x = x(i).x - ishift;
		//index.x = ix(i) + ishift;

		//ishift = floorf(x(i).y);
		//pos.y = x(i).y - ishift;
		//index.y = iy(i) + ishift;


		//ishift = floorf(x(i).z);
		//pos.z = x(i).z - ishift;
		//index.z = iz(i) + ishift;



		//index.x = ((index.x%fields->nx)+fields->nx)%fields->nx;
		//index.y = ((index.y%fields->ny)+fields->ny)%fields->ny;
		//index.z = ((index.z%fields->nz)+fields->nz)%fields->nz;

	//	printf("read E = %f, %f, %f, %i, %i, %i\n",pos.x,pos.y,pos.z,index.x,index.y,index.z);

		E_out(i).x = fields->intrpE(pos.x,pos.y,pos.z,
											index.x,index.y,index.z,0,ideriv);
		E_out(i).y = fields->intrpE(pos.x,pos.y,pos.z,
											index.x,index.y,index.z,1,ideriv);
		E_out(i).z = fields->intrpE(pos.x,pos.y,pos.z,
											index.x,index.y,index.z,2,ideriv);

	//	printf("E_out(%i) = %e, %e, %e\n",i,E_out(i).x,E_out(i).y,E_out(i).z);
	}

	return E_out;
}

template<int N>
template<enum FieldData_deriv ideriv>
typevecN<typevecN<realkind,N>,3> ParticleObjN<N>::get_E(FieldData* fields,
						const typevecN<realkind,N>& x,
						const typevecN<realkind,N>& y,
						const typevecN<realkind,N>& z,
						const typevecN<int,N>& ix_out,
						const typevecN<int,N>& iy_out,
						const typevecN<int,N>& iz_out)
{

	typevecN<typevecN<realkind,N>,3> E_out;

	if(fields->ndimensions == 1)
	{
		for(int i=0; i<N; i++)
		{
			E_out(0)(i) = fields->intrpE(x(i),y(i),z(i),
												ix_out(i),iy_out(i),iz_out(i),0,ideriv);
		}

		E_out(1) = 0.0;
		E_out(2) = 0.0;
	}
	else
	{
		for(int j=0;j<3;j++)
		{

			for(int i=0; i<N; i++)
			{
				E_out(j)(i) = fields->intrpE(x(i),y(i),z(i),
													ix_out(i),iy_out(i),iz_out(i),j,ideriv);
			}
		}
	}

//	for(int i=0;i<N;i++)
//		printf("E_out(%i) = %e @ %e, %i\n",i,E_out(0)(i),x(i),ix_out(i));



	return E_out;
}

template<int N>
typevecN<typevecN<realkind,N>,3> ParticleObjN<N>::get_accel(FieldData* fields,
							typevecN<realkind,N>& x,
							typevecN<realkind,N>& y,
							typevecN<realkind,N>& z,
							const typevecN<realkind,N>& vx,
							const typevecN<realkind,N>& vy,
							const typevecN<realkind,N>& vz)
{


	typevecN<realkind,N> x_out,y_out,z_out;
	typevecN<int,N> ix_out,iy_out,iz_out;

	ix_out = ifloor(x);
	x_out = x - ix_out;
	ix_out += ix;

	iy_out = ifloor(y);
	y_out = y - iy_out;
	iy_out += iy;

	iz_out = ifloor(z);
	z_out = z - iz_out;
	iz_out += iz;

	typevecN<typevecN<realkind,N>,3> accel;

	typevecN<typevecN<realkind,N>,3> E;


	E = get_E<FieldData_deriv_f>(fields,x_out,y_out,z_out,ix_out,iy_out,iz_out);

	if(fields->pdata->Bmag_avg > 0.0f)
	{
		typevecN<typevecN<realkind,N>,3> B;
		B = get_B<FieldData_deriv_f>(fields,x_out,y_out,y_out,ix_out,iy_out,iz_out);
		accel = E + cross_productN3(vx,vy,vz,B);
	}
	else
		accel = E;// + cross_productN3(vx,vy,vz,B);


	accel = accel * (fields->q2m[species]) * qe2me;
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
ploop<1> DummyFUNCTION1(PlasmaData* pdata,FieldData* fields,
		CurrentTally* currents,typevecN<int,1>& iter,const int nSubcycl_max)
{
	ploop<1> ploop_out;
	ploop_out.particle.push(pdata,fields,currents,iter,nSubcycl_max);

	return ploop_out;
}

/*
ploop<2> DummyFUNCTION1(PlasmaData* pdata,FieldData* fields,
		CurrentTally* currents,typevecN<int,2>& iter,const int nSubcycl_max)
{
	ploop<2> ploop_out;
	ploop_out.particle.push(pdata,fields,currents,iter,nSubcycl_max);

	return ploop_out;
}

ploop<4> DummyFUNCTION1(PlasmaData* pdata,FieldData* fields,
		CurrentTally* currents,typevecN<int,4>& iter,const int nSubcycl_max)
{
	ploop<4> ploop_out;
	ploop_out.particle.push(pdata,fields,currents,iter,nSubcycl_max);

	return ploop_out;
}

ploop<6> DummyFUNCTION1(PlasmaData* pdata,FieldData* fields,
		CurrentTally* currents,typevecN<int,6>& iter,const int nSubcycl_max)
{
	ploop<6> ploop_out;
	ploop_out.particle.push(pdata,fields,currents,iter,nSubcycl_max);

	return ploop_out;
}

ploop<8> DummyFUNCTION1(PlasmaData* pdata,FieldData* fields,
		CurrentTally* currents,typevecN<int,8>& iter,const int nSubcycl_max)
{
	ploop<8> ploop_out;
	ploop_out.particle.push(pdata,fields,currents,iter,nSubcycl_max);

	return ploop_out;
}

ploop<16> DummyFUNCTION1(PlasmaData* pdata,FieldData* fields,
		CurrentTally* currents,typevecN<int,16>& iter,const int nSubcycl_max)
{
	ploop<16> ploop_out;
	ploop_out.particle.push(pdata,fields,currents,iter,nSubcycl_max);

	return ploop_out;
}



ploop<32> DummyFUNCTION1(PlasmaData* pdata,FieldData* fields,
		CurrentTally* currents,typevecN<int,32>& iter,const int nSubcycl_max)
{
	ploop<32> ploop_out;
	ploop_out.particle.push(pdata,fields,currents,iter,nSubcycl_max);

	return ploop_out;
}




*/










