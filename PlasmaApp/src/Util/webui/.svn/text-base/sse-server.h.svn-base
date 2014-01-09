#ifndef SSE_SERVER_H
#define SSE_SERVER_H
/* Sample TCP server */
//#define __GXX_EXPERIMENTAL_CXX0X__
#include <pthread.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <pthread.h>
#include <mutex>
#include <cstdlib>
#include <omp.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <stdio.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <stdlib.h>
#include <string.h>
#include <arpa/inet.h>             // for the use of inet_ntop()
#include <netinet/in.h>
#include<unistd.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<sys/socket.h>
#include<arpa/inet.h>
#include<netdb.h>
#include<signal.h>
#include<fcntl.h>
#include <regex.h>
#include <algorithm>
#include <Teuchos_ParameterList.hpp>
#include "../../PlasmaData.h"


#define CONN_MAX 1000
class PlasmaData;
class sseDataStream;


class InputParser;
class WebControlSystem
{
public:

	WebControlSystem(){}

	WebControlSystem(PlasmaData* _pdata,InputParser* parser);

	void SetiRun(int iStatus)
	{
		io_mutex.lock();
		iRun = iStatus;
		io_mutex.unlock();
	}

	int CheckiRun()
	{
		int result;
		io_mutex.lock();
		result = iRun;
		io_mutex.unlock();
		return result;
	}


	PlasmaData* pdata;
	Teuchos::ParameterList* params;
private:
	std::mutex io_mutex;
	int iRun;


};

bool replace(std::string& str, const std::string& from, const std::string& to);


void* fork_server(void* data);
void* fork_connection(void* data);
void* stream_connection(void* data);

class sseServer
{
public:

	int listenfd,connfd[CONN_MAX];
	struct sockaddr_in servaddr,cliaddr[CONN_MAX];
	socklen_t clilen[CONN_MAX];
	char PORTc[6];
	int PORTi;
	char* ROOT;

	struct addrinfo *resd;

	int* exit_flag;
	std::mutex io_mutex[CONN_MAX];


	pthread_t conn_threads[CONN_MAX];

	sseDataStream* stream;
	WebControlSystem* pgrm_controls;

	sseServer(){};

	sseServer(char* _ROOT,char _PORTc[6],sseDataStream* _stream,
			WebControlSystem* _pgrm_controls,int* _exit_flag);

	// binds the listen socket
	int StartServer(void);

	// Creates a new thread for the server, and starts listening for
	// connections.
	int RunServer(void);

	int HandleResponse(int iCon);

	int ProcessPost(int iCon)
	{

		pgrm_controls->pdata->Lx *= 2.0;
		return 0;
	}


	int EventStream(int iCon);

};

typedef struct conn_data
{
public:
	sseServer* server;
	pthread_t* pid;
	int* iCon;
	int slot;
} strconn_data;









#endif
