 
template<typename T>
class ArrayManager;

class FlexBC
{
public:
	virtual realkind CalcBC(ArrayManager<realkind>& array,int i, int j, int k);
};

template<typename T>
class ArrayManager
{
public:

	ArrayManager(T* _data,int _nbuff,int _nd,int _nx,int _ny,int _nz,int _nelm) :
		nbuff(_nbuff),data(_data),nx(_nx),ny(_ny),nz(_nz),nelm(_nelm),
		nxf(2*_nbuff+_nx),nyf(2*_nbuff+_ny),nzf(2*_nbuff+_nz),nD(_nd)
	{
		if(nD == 1)
			off_data = data + hash_function(nbuff);
		else if(nD == 2)
			off_data = data + hash_function(nbuff,nbuff);
		else if(nD == 3)
			off_data = data + hash_function(nbuff,nbuff,nbuff);

		ntotal = hash_function(nxf-1,nyf-1,nzf-1,nelm-1);
	};

	int hash_function(const int i, const int j, const int k,const int Elm)
	const
	{
		return i + nxf*(j+nyf*(k+nzf*Elm));
	}

	int hash_function(const int i, const int j, const int Elm)
	const
	{
		return i + nxf*(j+nyf*(Elm));
	}

	int hash_function(const int i,const int Elm)
	const
	{
		return i + nxf*((Elm));
	}

	int hash_function(const int i)
	const
	{
		return i;
	}

	const T& Get(const int i, const int j, const int k, const int elm)
	const
	{
		return off_data[hash_function(i,j,k,elm)];
	}

	const T& Get(const int i, const int j,  const int elm)
	const
	{
		return off_data[hash_function(i,j,elm)];
	}

	const T& Get(const int i, const int elm)
	const
	{
		return off_data[hash_function(i,elm)];
	}

	const T& Get(const int i)
	const
	{
		return off_data[hash_function(i)];
	}

	T& Get(const int i, const int j, const int k, const int elm)
	{
		return off_data[hash_function(i,j,k,elm)];
	}

	T& Get(const int i, const int j,  const int elm)
	{
		return off_data[hash_function(i,j,elm)];
	}

	T& Get(const int i, const int elm)
	{
		return off_data[hash_function(i,elm)];
	}

	T& Get(const int i)
	{
		return off_data[hash_function(i)];
	}

	T& operator()(const int i, const int j, const int k, const int elm)
	{
		return off_data[hash_function(i,j,k,elm)];
	}

	T& operator()(const int i, const int j,  const int elm)
	{
		return off_data[hash_function(i,j,elm)];
	}

	T& operator()(const int i, const int elm)
	{
		return off_data[hash_function(i,elm)];
	}

	T& operator()(const int i)
	{
		return off_data[hash_function(i)];
	}

	T Mean()
	{
		T result = 0;
		for(int k=0;k<nz;k++)
			for(int j=0;j<ny;j++)
				for(int i=0;i<nx;i++)
				{
					result += Get(i,j,k);
				}

		return result/((double)nx*ny*nz);
	}


	void UpdateBoundaries(FlexBC* Bc)
	{
		for(int elm=0;elm<nelm;elm++)
		{
			if(nD == 1)
			{
				for(int i=-nbuff;i<nx+nbuff;i++)
				{

				}
			}
			else if(nD == 2)
			{
				for(int k=-nbuff;k<nz+nbuff;k++)
				{
					for(int j=-nbuff;j<ny+nbuff;j++)
					{

					}
				}
			}
			else if(nD == 3)
			{
				for(int k=-nbuff;k<nz+nbuff;k++)
				{

				}
			}

		}

	}


private:

	T const* data;
	T* off_data;

	const int nD;
	const int nx,ny,nz,nelm;
	const int nbuff;
	const int nxf,nyf,nzf;
	const int ntotal;
};
