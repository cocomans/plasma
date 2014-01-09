/*
 * BoundaryValues.h
 *
 *  Created on: Dec 11, 2013
 *      Author: payne
 */

#ifndef BOUNDARYVALUES_H_
#define BOUNDARYVALUES_H_

#include "StateVector.h"
namespace HoLo
{
namespace LowOrder
{
	class BoundaryValues
	{
	public:

		class BCs{public:enum{north,south,east,west,LF};};
		class BCType{public:enum {Periodic, Dirichlet, Neumann};};

		int nx,ny,nz;
		int nBufx,nBufy,nBufz;
		SimParams* params;

		int north[SimParams::Vars::LF];
		int south[SimParams::Vars::LF];
		int east[SimParams::Vars::LF];
		int west[SimParams::Vars::LF];

		BoundaryValues(SimParams* _params)
	{
			params = _params;

			// Default values for boundary conditions
			for(int iVar=0;iVar<SimParams::Vars::LF;iVar++)
			{
				for(int iBC=0;iBC<BCs::LF;iBC++)
				{
					GetBC(iBC,iVar) = BCType::Periodic;
				}
			}


	}

		int& GetBC(int iBC,int iVar)
		{

			int* result;
			switch(iBC)
			{
			case BCs::north:
				result = north+iVar;
				break;
			case BCs::south:
				result = south+iVar;
				break;
			case BCs::east:
				result = east+iVar;
				break;
			case BCs::west:
				result = west+iVar;
				break;
			default: break;
			}

			return *result;
		}

		template<class Getter>
		void ApplyBCs(Getter op,int iVar)
		{
			for(int l=0;l<BCs::LF;l++)
			{
				if(GetBC(l,iVar) == BCType::Periodic)
				{
					ApplyPeriodic(op,l);
				}
			}
		}


		template<class Getter>
		void ApplyPeriodic(Getter op,int dir)
		{
			int nBuf,n,i,j,k;
			int iget,jget,kget;
			switch(dir)
			{
			case BCs::east: // x direction
				nBuf = nBufx;
				n = nx;
				break;
			default:
				break;
			}
			for(int i=-nBuf;i<n+nBuf;i++)
			{
				int iget = ((i%n)+n)%n;
				op(i,0,0) = op(iget,0,0);
			}
		}

		template<class Getter>
		void ApplyDirchlet(Getter op)
		{
			for(int i=-nBufx;i<nx+nBufx;i++)
			{
				int iget = ((i%nx)+nx)%nx;
				op(i) = op(iget);
			}
		}

		template<class Getter>
		void ApplyNeumann(Getter op)
		{
			for(int i=-nBufx;i<nx+nBufx;i++)
			{
				int iget = ((i%nx)+nx)%nx;
				op(i) = op(iget);
			}
		}
	};



} /* namespace LowOrder */
} /* namespace HoLo */
#endif /* BOUNDARYVALUES_H_ */
