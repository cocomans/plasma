#ifndef SSE_DATA_STREAM_H
#define SSE_DATA_STREAM_H
#include <pthread.h>

#include <iostream>
#include <stdlib.h>
#include <mutex>

class PlasmaData;
class googleChart;
class LogStream;

class sseDataStream
{
public:

	sseDataStream(){};

	sseDataStream(PlasmaData* _pdata,int _nStreams, int _nCharts);

	sseDataStream(PlasmaData* _pdata);


	bool getMessage(std::string &msg_out, int& updateIndex_in);

//	void getAllMessages(std::string &msg_out);

	PlasmaData* pdata;
	LogStream** dStreams;
	std::string** rawStreams;
	googleChart** charts;
	int nStreams;
	int nCharts;
	int ParamsOnly;


	// All of the messages
	std::string* AllMsg;
	// Last message sent
	std::string* LastMsg;

private:
	std::mutex io_mutex;
	int updateIndex;
};

#endif


