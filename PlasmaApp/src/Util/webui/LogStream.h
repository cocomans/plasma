#ifndef LOG_STREAM_H
#define LOG_STREAM_H
#include <pthread.h>

#include <iostream>
#include <stdlib.h>
#include <mutex>
#include <string>
#include <iostream>
#include <sstream>
class LogStream
{
public:

	LogStream(char* _name);

	template <class ...Args>
	void print(const char* templ,const Args& ...args)
	{
	    io_mutex->lock();
			char temp[1024];
			sprintf(temp,templ,args...);
			printf("%s",temp);
	    	sstream += temp;
	    	qNewUpdates = 1;
	    io_mutex->unlock();
	}

	template<class T>
	void operator<<(T& in)
	{
	    io_mutex->lock();
	    	std::stringstream tmp;
	    	tmp << in;
	    	sstream += tmp.str();
	    	qNewUpdates = 1;
	    io_mutex->unlock();
	}

	bool getStream(std::string& result);

	std::string getName(void);

	std::string sstream;
private:
	std::string name;
	std::mutex* io_mutex;
	int qNewUpdates;


};

typedef struct gChartInfo
{
	std::string name;
	float* ys;


}ChartInfo;

class googleChart
{
public:

	googleChart(const char* _name,const char* _xname):name(_name),xname(_xname)
	{
		io_mutex = new std::mutex();

	};

	void plotLines(int _nplots,int length,float* x,...);

	void appendPoints(int _nplots,int length,float* x,...);


	void getStream(std::string& result);
	std::string data;
	std::string output;
private:
	int nplots;
	std::string name;
	std::string xname;
	std::mutex* io_mutex;

};

//class googleMotionChart : public googleChart
//{
//public:
//
//	googleMotionChart(const char* _name,const char* _xname):name(_name),xname(_xname)
//	{
//		io_mutex = new std::mutex();
//
//	};
//
//	void plotLines(int _nplots,int length,std::string* entitiy,...);
//
//	void appendPoints(int _nplots,int length,std::string* entitiy,...);
//
//
//	std::string data;
//	std::string output;
//private:
//	int nplots;
//	std::string name;
//	std::string xname;
//	std::mutex* io_mutex;
//
//};

#endif
