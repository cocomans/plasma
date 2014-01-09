 
#include "LogStream.h"
#include "../../PlasmaData.h"
#include "../../RunData.h"
#include "../../OutputCtrl.h"
#include <algorithm>


LogStream::LogStream(char* _name)
{
	name = _name;

	qNewUpdates=0;
	io_mutex = new std::mutex();
}



bool LogStream::getStream(std::string& result)
{
	bool qUpdated = 0;

	io_mutex->lock();
//		if(qNewUpdates)
//		{
//			sstream->seekg(0);
		char inBuf[512];
		std::stringstream tmp(sstream);
//			std::cout << "stream:" << name <<"\n" << sstream << std::endl;
		int nlines = std::count(std::istreambuf_iterator<char>(tmp),std::istreambuf_iterator<char>(),'\n');
//			printf("there are %i lines in stream\n",nlines);
		tmp.seekg(0);
		int iline = 0;
		while(tmp.getline(inBuf,512))
		{

			if(iline >= (nlines-100)){
			std::string temp = inBuf;

			result +=  "\"" + temp + "\",";
			}

			iline++;
		}

		int lastc = result.rfind(",");
		if(lastc != std::string::npos)
			result.erase(lastc,1);

		qNewUpdates = 0;
		qUpdated = 1;
//		}
	io_mutex->unlock();

	return qUpdated;
}







std::string LogStream::getName(void)
{
	return name;
}

void googleChart::plotLines(int _nplots,int length,float* x,...)
{
	nplots = _nplots;
	int nmax = 500;
	int stride = (length + nmax -1)/nmax;
	io_mutex->lock();
	output = "";
	data = "";
	std::stringstream datas(data);
	datas << "\"" << xname << "\",";
	va_list args;
	va_start(args,nplots);
	for(int j=0;j<nplots;j++)
	{
		datas << "\"" << va_arg(args,const char*) << "\",";
		double* dummy = va_arg(args,double*);

	}


	for(int i=0;i<length;i+=stride)
	{
		va_start(args,nplots);
		datas << x[i] << ",";
		for(int j=0;j<nplots;j++)
		{
			const char* dummy = va_arg(args,const char*);
			datas << va_arg(args,float*)[i] << ",";

		}
		va_end(args);
	}

	data = datas.str();

	int lastc = data.rfind(",");
	if(lastc != std::string::npos)
		data.erase(lastc,1);

	io_mutex->unlock();


}

void googleChart::appendPoints(int _nplots,int length,float* x,...)
{

	int nmax = 500;
	int stride = (length + nmax -1)/nmax;
	io_mutex->lock();
	data += ",";
	std::stringstream datas(data);
	va_list args;
	va_start(args,nplots);



	for(int i=0;i<length;i+=stride)
	{
		va_start(args,nplots);
		datas << x[i] << ",";
		for(int j=0;j<nplots;j++)
		{
			const char* dummy = va_arg(args,const char*);
			datas << va_arg(args,float*)[i] << ",";

		}
		va_end(args);
	}



	data = datas.str();

	int lastc = data.rfind(",");
	if(lastc != std::string::npos)
		data.erase(lastc,1);



	io_mutex->unlock();


}

void googleChart::getStream(std::string& result)
{
	io_mutex->lock();



	std::string tmpname = name;
	int lastc = tmpname.rfind(" ");
	if(lastc != std::string::npos)
		tmpname.erase(lastc,1);
	std::stringstream outputs(output);

	outputs << "\"" << tmpname << "\":{"
			"\"title\": \"" << name <<"\","
			"\"nplots\":" << nplots <<","
			"\"data\": ["<< data << "]}";



	result = outputs.str();
	io_mutex->unlock();
}


