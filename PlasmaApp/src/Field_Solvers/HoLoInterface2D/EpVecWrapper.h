#ifndef EP_VEC_WRAPPER_H
#define EP_VEC_WRAPPER_H
//Trilinos Includes
#include "NOX_Epetra.H"
#include <Teuchos_FancyOStream.hpp>
#include "../HoLoInterface/SimParams.h"
class MapManager2D;


class EpVecWrapper
{
public:
	EpVecWrapper(Teuchos::RCP<MapManager2D> _map,
			const Teuchos::RCP<SimParams> &params);

	EpVecWrapper(Epetra_Vector& _x,const Teuchos::RCP<SimParams> &params);

	EpVecWrapper(const Epetra_Vector& _x,const Teuchos::RCP<SimParams> &params);


	EpVecWrapper& operator=(const Epetra_Vector& _x);

	EpVecWrapper& operator=(Epetra_Vector& _x);

	void PutScalar(double val)
	{
		X->PutScalar(val);
	}

	double& operator()(const int elm, int i, int j)
	{
		i = (i&(nx-1)+nx)&(nx-1);
		j = (j&(ny-1)+ny)&(ny-1);
		return (*X)[elm+nUnk*(i+nx*(j))];
	}

	const double& operator()(const int elm, int i, int j)
	const
	{
		i = (i&(nx-1)+nx)&(nx-1);
		j = (j&(ny-1)+ny)&(ny-1);
		return (*X)[elm+nUnk*(i+nx*(j))];
	}

	double& operator()(const int elm, int i)
	{
		i = (i&(nx-1)+nx)&(nx-1);
		return (*X)[elm+nUnk*(i)];
	}

	const double& operator()(const int elm, int i)
	const
	{
		i = (i&(nx-1)+nx)&(nx-1);
		return (*X)[elm+nUnk*(i)];
	}

	double& Get(const int elm, int i, int j=0)
	{
		i = (i&(nx-1)+nx)&(nx-1);
		j = (j&(ny-1)+ny)&(ny-1);
		return (*X)[elm+nUnk*(i+nx*(j))];
	}

	const double& Get(const int elm, int i, int j=0)
	const
	{
		i = (i&(nx-1)+nx)&(nx-1);
		j = (j&(ny-1)+ny)&(ny-1);
		return (*X)[elm+nUnk*(i+nx*(j))];
	}


	const double GetC2F(const int elm, int i, int j=0)
	const
	{
		return 0.5*(Get(elm,i-1,j) + Get(elm,i,j));
	}

	const double GetF2C(const int elm, int i, int j=0)
	const
	{

		return 0.5*(Get(elm,i+1,j) + Get(elm,i,j));
	}

	const double GetyF2xF(const int elm, int i, int j=0)
	const
	{

		return 0.25*(Get(elm,i-1,j) + Get(elm,i-1,j+1)
				+ Get(elm,i,j) + Get(elm,i,j+1));
	}

	const double GetxF2yF(const int elm, int i, int j=0)
	const
	{

		return 0.25*(Get(elm,i,j) + Get(elm,i,j-1)
				+ Get(elm,i+1,j) + Get(elm,i+1,j-1));
	}

	const double GetC2xF(const int elm, int i, int j=0)
	const
	{
		return 0.5*(Get(elm,i-1,j) + Get(elm,i,j));
	}

	const double GetC2yF(const int elm, int i, int j=0)
	const
	{
		return 0.5*(Get(elm,i,j-1) + Get(elm,i,j));
	}

	const double GetxF2C(const int elm, int i, int j=0)
	const
	{

		return 0.5*(Get(elm,i+1,j) + Get(elm,i,j));
	}
	const double GetyF2C(const int elm, int i, int j=0)
	const
	{

		return 0.5*(Get(elm,i,j+1) + Get(elm,i,j));
	}

	const double GetC2V(const int elm, int i, int j=0)
	const
	{

		return 0.25*(Get(elm,i,j) + Get(elm,i-1,j)
				+ Get(elm,i,j-1) + Get(elm,i-1,j-1));
	}

	const double GetV2C(const int elm, int i, int j=0)
	const
	{

		return 0.25*(Get(elm,i,j) + Get(elm,i+1,j)
				+ Get(elm,i,j+1) + Get(elm,i+1,j+1));
	}


	const double GetV2xF(const int elm, int i, int j=0)
	const
	{

		return 0.5*(Get(elm,i,j) + Get(elm,i,j+1));
	}

	const double GetV2yF(const int elm, int i, int j=0)
	const
	{

		return 0.5*(Get(elm,i,j) + Get(elm,i+1,j));
	}

	const double DxC2V(const int elm, int i, int j=0)
	const
	{

		return 0.5*(DxC2F(elm,i,j) + DxC2F(elm,i,j-1));
	}

	const double DyC2V(const int elm, int i, int j=0)
	const
	{

		return 0.5*(DyC2F(elm,i,j) + DyC2F(elm,i-1,j));
	}

	const double DxV2C(const int elm, int i, int j=0)
	const
	{

		return 0.5*(DxF2C(elm,i,j) + DxF2C(elm,i,j+1));
	}

	const double DyV2C(const int elm, int i, int j=0)
	const
	{

		return 0.5*(DyF2C(elm,i,j) + DyF2C(elm,i+1,j));
	}



	const double Dx(const int elm, int i, int j=0)
	const
	{
		return (Get(elm,i+1,j)-Get(elm,i-1,j))/(2.0*dx);
	}

	const double Dy(const int elm, int i, int j=0)
	const
	{
		return (Get(elm,i,j+1)-Get(elm,i,j-1))/(2.0*dy);
	}

	const double DxC2F(const int elm, int i, int j=0)
	const
	{
		return (Get(elm,i,j)-Get(elm,i-1,j))/(dx);
	}

	const double DyC2F(const int elm, int i, int j=0)
	const
	{
		return ((*this)(elm,i,j)-(*this)(elm,i,j-1))/(dy);
	}

	const double DyC2xF(const int elm, int i, int j=0)
	const
	{
		return 0.5*(
				Dy(elm,i,j)
				+
				Dy(elm,i-1,j)
				);
	}

	const double DxC2yF(const int elm, int i, int j=0)
	const
	{
		return 0.5*(
				Dx(elm,i,j)
				+
				Dx(elm,i,j-1)
				);
	}

	const double DxF2C(const int elm, int i, int j=0)
	const
	{
		return ((*this)(elm,i+1,j)-(*this)(elm,i,j))/(dx);
	}

	const double DyF2C(const int elm, int i, int j=0)
	const
	{
		return ((*this)(elm,i,j+1)-(*this)(elm,i,j))/(dy);
	}

	const double Dxx(const int elm, int i, int j=0)
	const
	{
		return (Get(elm,i+1,j) - 2.0*Get(elm,i,j) + Get(elm,i-1,j))/(dx*dx);
	}

	const double Dyy(const int elm, int i, int j=0)
	const
	{
		return (Get(elm,i,j+1) - 2.0*Get(elm,i,j) + Get(elm,i,j-1))/(dy*dy);
	}

	const double DxC2F(double W, double C, double E)
	const
	{
	  return ( C - W ) / dx;
	}

	double DxF2C(double W, double C, double E)
	const
	{
	  return ( E - C ) / dx;
	}


	const double mean(const int elm)
	const
	{
		double result = 0;
		for(int j=0;j<ny;j++)
			for(int i=0;i<nx;i++)
				result += (*this)(elm,i,j);

		return result/((double)nx*ny);
	}

	double subNorm2(const int n,...)
	{
		va_list args;
		double result = 0;

		for(int j=0;j<ny;j++)
		{
			for(int i=0;i<nx;i++)
			{
				va_start(args,n);
				for(int l=0;l<n;l++)
				{
					int elm = va_arg(args,int);
					result += pow(Get(elm,i,j),2.0);
				}
				va_end(args);

			}
		}

		return sqrt(result/((double)nx*ny*n));
	}



	Teuchos::RCP<SimParams> simParams;
	Teuchos::RCP<Epetra_Vector> X;
	int nx,ny,nUnk;
	double dx,dy;

};


































#endif /* EP_VEC_WRAPPER_H */
