 
#include "EpVecWrapper.h"
#include  "../HoLoInterface/SimParams.h"
#include "MapManager2D.h"

EpVecWrapper::EpVecWrapper(Teuchos::RCP<MapManager2D> _map,
		const Teuchos::RCP<SimParams> &params)
{
	X = Teuchos::rcp(new Epetra_Vector(*(_map->GetPtSdMap() )));
	simParams = params;
	nx = params->pdata->nx;
	ny = params->pdata->ny;
	if(params->pdata->ndimensions == 1)
		ny = 1;

	nUnk = simParams->GetParamsMap()->get<int>("Number Unknowns");
	dx = params->pdata->dxdi;
	dy = params->pdata->dydi;
}

EpVecWrapper::EpVecWrapper(Epetra_Vector& _x,
		const Teuchos::RCP<SimParams> &params)
{
	X = Teuchos::rcp(new Epetra_Vector(View,_x,0));
	simParams = params;
	nx = params->pdata->nx;
	ny = params->pdata->ny;
	if(params->pdata->ndimensions == 1)
		ny = 1;

	nUnk = simParams->GetParamsMap()->get<int>("Number Unknowns");
	dx = params->pdata->dxdi;
	dy = params->pdata->dydi;
}

EpVecWrapper::EpVecWrapper(const Epetra_Vector& _x,
		const Teuchos::RCP<SimParams> &params)
{
	X = Teuchos::rcp(new Epetra_Vector(View,_x,0));
	simParams = params;
	nx = params->pdata->nx;
	ny = params->pdata->ny;
	if(params->pdata->ndimensions == 1)
		ny = 1;

	nUnk = simParams->GetParamsMap()->get<int>("Number Unknowns");
	dx = params->pdata->dxdi;
	dy = params->pdata->dydi;
}


EpVecWrapper& EpVecWrapper::operator=(const Epetra_Vector& _x)
{
	(*X) = _x;

	return *this;
}

EpVecWrapper& EpVecWrapper::operator=(Epetra_Vector& _x)
{
	(*X) = _x;

	return *this;
}
