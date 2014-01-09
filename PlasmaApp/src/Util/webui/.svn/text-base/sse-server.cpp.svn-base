 
#include "sse-server.h"
#include "sseDataStream.h"
#include "../../PlasmaData.h"

bool replace(std::string& str, const std::string& from, const std::string& to) {

	size_t start_pos = 0;

	while(true)
	{

		start_pos = str.find(from,start_pos);

		if(start_pos == std::string::npos)
			break;


		str.replace(start_pos, from.length(), to);

		start_pos += to.length();

	}
    return true;
}

sseServer::sseServer(char* _ROOT,char _PORTc[6],sseDataStream* _stream,
		WebControlSystem* _pgrm_controls,int* _exit_flag)
{
	stream = _stream;
	pgrm_controls = _pgrm_controls;
	exit_flag = _exit_flag;
    ROOT = getenv("PWD");
    strcpy(PORTc,_PORTc);
    char c;

    int slot=0;


    PORTi = atoi(PORTc);

    printf("Server started at port no. %s%s%s with root directory as %s%s%s\n","\033[92m",PORTc,"\033[0m","\033[92m",ROOT,"\033[0m");

    for(int i=0;i<CONN_MAX;i++)
    	connfd[i] = -1;
}

int sseServer::StartServer(void)
{
	int ierr = 0;

    struct addrinfo hints, *res, *p;

    // getaddrinfo for host
    memset (&hints, 0, sizeof(hints));
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_flags = AI_PASSIVE;
    if (getaddrinfo( NULL, PORTc, &hints, &res) != 0)
    {
//        perror ("getaddrinfo() error");
        exit(1);
    }
    // socket and bind
    for (p = res; p!=NULL; p=p->ai_next)
    {
        listenfd = socket (p->ai_family, p->ai_socktype, 0);
        if (listenfd == -1) continue;
        if (bind(listenfd, p->ai_addr, p->ai_addrlen) == 0) break;
    }
    if (p==NULL)
    {
        perror ("socket() or bind()");
        exit(1);
    }

    freeaddrinfo(res);

	pthread_t* pid = new pthread_t();

	strconn_data* data = new strconn_data();

	data->server = this;
	data->pid = pid;

	// Server runs on its own thread
	ierr = pthread_create(pid,NULL,fork_server,(void*)data);
	if (ierr)
	{
		printf("ERROR; return code from pthread_create() is %d\n", ierr);
		exit(-1);
	}

    return ierr;
}

void* fork_server(void* data)
{
	strconn_data* data2 = (strconn_data*)data;

	// Run the server
	data2->server->RunServer();

	// cleanup
	free(data2->pid);
	free(data2);
	pthread_exit(0);
	return NULL;
}

void* fork_connection(void* data)
{
	strconn_data* data2 = (strconn_data*)data;

	// Run the server
	data2->server->HandleResponse(data2->slot);

	return NULL;
}

void* stream_connection(void* data)
{
	strconn_data* data2 = (strconn_data*)data;

	// Run the server
	data2->server->EventStream(data2->slot);

	// cleanup
	free(data2->pid);
	free(data2);
	pthread_exit(0);
	return NULL;
}


int sseServer::RunServer(void)
{

	int ierr = 0;
	strconn_data data[CONN_MAX];

	printf("listenfd = %i\n",listenfd);
    // listen for incoming connections
    if ( listen (listenfd, 1024) != 0 )
    {
        perror("listen() error");
        exit(1);
    }

	int slot = 0;
	for(;;)
	{


		clilen[slot]=sizeof(cliaddr[slot]);
		connfd[slot] = accept(listenfd,(struct sockaddr *)&(cliaddr[slot]),&(clilen[slot]));
		if(*exit_flag)
		{
			shutdown(connfd[slot],SHUT_RDWR);
			close(connfd[slot]);
			break;
		}
		printf("got connection on fd[%i] = %i\n",slot,connfd[slot]);
		data[slot].pid = &(conn_threads[slot]);;
		data[slot].server = this;
		data[slot].iCon = &(connfd[slot]);
		data[slot].slot = slot;

		int rc = pthread_create(&(conn_threads[slot]),NULL,fork_connection,(void*)&(data[slot]));
		if (rc) {
		printf("ERROR; return code from pthread_create() is %d\n", rc);
		exit(-1);
		}

		while (connfd[slot]!=-1) slot = (slot+1)%100;



	}

	shutdown(listenfd,SHUT_RDWR);
	close(listenfd);


	return ierr;
}

int sseServer::HandleResponse(int iCon)
{
	int streamflag = 0;
	char imesg[99999];
	char mesg[99999], *reqline[3], data_to_send[4096], path[99999];
	int fd, bytes_read;
	int rcvd;

	rcvd = recvfrom(connfd[iCon],imesg,99999,0,(struct sockaddr *)&(cliaddr[iCon]),&(clilen[iCon]));

	if (rcvd<0)    // receive error
		fprintf(stderr,("recv() error\n"));
	else if (rcvd==0)    // receive socket closed
		fprintf(stderr,"Client disconnected upexpectedly.\n");
	else
	{
		std::string sendstuff = "";
		std::string header  = "";
		std::cout << "Recieved:\n" << imesg << std::endl;

		reqline[0] = strtok (imesg, " \t\n");
		if ( strncmp(reqline[0], "GET\0", 4)==0 )
		{
			reqline[1] = strtok (NULL, " \t");
			reqline[2] = strtok (NULL, " \t\n");

			if ( strncmp( reqline[2], "HTTP/1.0", 8)!=0 && strncmp( reqline[2], "HTTP/1.1", 8)!=0 )
			{
				write(connfd[iCon], "HTTP/1.1 400 Bad Request\n", 25);
			}
			else
			{
				if ( strncmp(reqline[1], "/\0", 2)==0 )
				{
					reqline[1] = "/src/Util/webui/index.html";
				}

				std::string req = reqline[1];
				// Check to see if a stream event is being requested
				if((req.find("/stream")) != std::string::npos)
				{

					// A stream request has been made, we need to continously send to
					// the stream in a new thread.
					strconn_data* data = new strconn_data();
					streamflag = 1;
					pthread_t* thread = new pthread_t;
					data->server = this;
					data->pid = thread;
					data->iCon = connfd+iCon;
					data->slot = iCon;


					int rc = pthread_create(thread,NULL,stream_connection,(void*)(data));
					if (rc) { printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
					}

					pthread_exit(0);
				}
				else
				{
					strcpy(path, ROOT);
					strcpy(&path[strlen(ROOT)], reqline[1]);
					printf("file: %s\n", path);

					std::ifstream is (path);


					if ( is )    //FILE FOUND
					{
						std::string temp = "";
						header = "HTTP/1.1 200 OK\nContent-Type: text/html\n\n";
					    is.seekg (0, is.end);
					    int length = is.tellg();
					    is.seekg (0, is.beg);
					    char * buffer = new char [length];
					    is.read (buffer,length);

					    temp = buffer;

						replace(temp,std::string("<!PORT_NUMBER!>"),std::string(PORTc));

						sendstuff = temp + "\n\n\0";
						std::cout << header + sendstuff << std::endl;


					}
					else    header ="HTTP/1.1 404 Not Found\n\n"; //FILE NOT FOUND
					std::string temp = (header+sendstuff);

					io_mutex[connfd[iCon]].lock();
					send(connfd[iCon], temp.c_str(), temp.length(),0);
					io_mutex[connfd[iCon]].unlock();


				}
			}
		}
		else if ( strncmp(reqline[0], "POST\0", 4)==0 )
		{
			// Process a post
			ProcessPost(iCon);
			printf("Got a Post\n %s\n %s/n",reqline[1],reqline[2]);
			if ( strncmp(reqline[1], "/\0", 2)==0 )
			{
				reqline[1] = "/src/Util/webui/index.html";
			}

			strcpy(path, ROOT);
			strcpy(&path[strlen(ROOT)], reqline[1]);
			printf("file: %s\n", path);

			std::ifstream is (path);


			if ( is )    //FILE FOUND
			{
				std::string temp = "";
				header = "HTTP/1.1 200 OK\nContent-Type: text/html\n\n";
			    is.seekg (0, is.end);
			    int length = is.tellg();
			    is.seekg (0, is.beg);
			    char * buffer = new char [length];
			    is.read (buffer,length);

			    temp = buffer;

				replace(temp,std::string("<!PORT_NUMBER!>"),std::string(PORTc));

				sendstuff = temp + "\n\n\0";
				std::cout << header + sendstuff << std::endl;


			}
			else    header ="HTTP/1.1 404 Not Found\n\n"; //FILE NOT FOUND
			std::string temp = (header+sendstuff);

			io_mutex[connfd[iCon]].lock();
			send(connfd[iCon], temp.c_str(), temp.length(),0);
			io_mutex[connfd[iCon]].unlock();


		}
	}

	shutdown (connfd[iCon], SHUT_RDWR);
	close(connfd[iCon]);
	connfd[iCon] = -1;
	pthread_exit(0);
}

int sseServer::EventStream(int iCon)
{
	std::string header = "HTTP/1.1 200 OK\nContent-Type: text/event-stream\nContent-Location: /stream\n\n";
//	send(connfd[iCon], header.c_str(), header.length(),0);
	std::string sendstuff;

	int updateIndex = 0;

	stream->getMessage(sendstuff,updateIndex);
	header += sendstuff;

	std::cout << "sending:\n" << header << std::endl;
	io_mutex[connfd[iCon]].lock();
	int send_err = send(connfd[iCon], header.c_str(), header.length(),MSG_NOSIGNAL);
	io_mutex[connfd[iCon]].unlock();
	if((*exit_flag) || (send_err <= 0))
	{
		printf("Connection %i closed, exiting stream\n",connfd[iCon]);
		shutdown (connfd[iCon], SHUT_RDWR);
		close(connfd[iCon]);
		connfd[iCon] = -1;

		return 0;
	}

	for(;;)
	{

		sendstuff = "";
		// Loop until we have a message to get
		while(!(stream->getMessage(sendstuff,updateIndex)))
		{
			printf("trying to get message\n");
			sleep(1);

		}

		std::cout << "sending:\n" << sendstuff << std::endl;
		io_mutex[connfd[iCon]].lock();
		int send_err = send(connfd[iCon], sendstuff.c_str(), sendstuff.length(),MSG_NOSIGNAL);
		io_mutex[connfd[iCon]].unlock();

		if((*exit_flag) || (send_err <= 0))
		{
			printf("Connection %i closed, exiting stream\n",connfd[iCon]);
			shutdown (connfd[iCon], SHUT_RDWR);
			close(connfd[iCon]);
			connfd[iCon] = -1;

			return 0;
		}
		sleep(1);


	}

	// Something went wrong if we ever make it here
	return 1;


}

WebControlSystem::WebControlSystem(PlasmaData* _pdata,InputParser* parser)
{
	parser->GetParam("ImmediatlyStartRun","--now",iRun,0);
	pdata = _pdata;

}

//int main(int argc, char* argv[])
//{
//	LogStream* logger1 = new LogStream("logger1");
//	LogStream* logger2 = new LogStream("logger2");
//	LogStream* logger3 = new LogStream("logger3");
//	LogStream* logger4 = new LogStream("logger4");
//
//	LogStream* loggers[4] = {logger1, logger2,logger3,logger4};
//
//	sseDataStream* stream = new sseDataStream(loggers,4);
//
//	WebControlSystem* controls = new WebControlSystem();
//
//	int exit_flag = 0;
//
//	sseServer* myserver = new sseServer(argc,argv,stream,controls,&exit_flag);
//
//	myserver->StartServer();
//
//	getchar();
//	omp_set_num_threads(4);
//
//#pragma omp parallel for
//	for(int i=0;i<10000;i++)
//	{
//		int tid = omp_get_thread_num();
//		logger1->print("l1: %i thread[%i]\n",i,tid);
//		logger2->print("l2: %i\n",i);
//		logger3->print("l3: %i\n",i);
//		logger4->print("l4: %i\n",i);
//
//		sleep(1.0+0.1*(rand()%10));
//	}
//
//	exit_flag = 1;
//
//	sleep(5);
//
//
//}


















