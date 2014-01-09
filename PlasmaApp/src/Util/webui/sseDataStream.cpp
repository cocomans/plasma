 
#include "sseDataStream.h"
#include "LogStream.h"
#include "../../PlasmaData.h"
#include "../../RunData.h"
#include "../../OutputCtrl.h"
#include "../CPUTimer.h"

#include <algorithm>

sseDataStream::sseDataStream(PlasmaData* _pdata,int _nStreams, int _nCharts)
{
	pdata = _pdata;


	nStreams = _nStreams;
	nCharts = _nCharts;

	dStreams = (LogStream**)malloc(nStreams*sizeof(LogStream*));
	dStreams[0] = &LO_log;
	dStreams[1] = &HO_log;
	dStreams[2] = &Outer_log;
	dStreams[3] = &debug_log;

	rawStreams = (std::string**)malloc(nStreams*sizeof(std::string*));

	for(int i=0;i<nStreams;i++)
		rawStreams[i] = (new std::string());

	AllMsg = new std::string();
	LastMsg = new std::string();
	updateIndex = 0;
	ParamsOnly = 0;
}

sseDataStream::sseDataStream(PlasmaData* _pdata)
{
	pdata = _pdata;


	nStreams = 4;
	nCharts = 0;

	dStreams = (LogStream**)malloc(nStreams*sizeof(LogStream*));
	dStreams[0] = &LO_log;
	dStreams[1] = &HO_log;
	dStreams[2] = &Outer_log;
	dStreams[3] = &debug_log;

	rawStreams = (std::string**)malloc(nStreams*sizeof(std::string*));

	for(int i=0;i<nStreams;i++)
		rawStreams[i] = (new std::string());

	AllMsg = new std::string();
	LastMsg = new std::string();
	updateIndex = 0;
	ParamsOnly = 1;

}

bool sseDataStream::getMessage(std::string &msg_out, int& updateIndex_in)
{



	bool res = false;
	io_mutex.lock();



//			for(int i=0;i<nStreams;i++)
//				res = res || dStreams[i]->getStream((temp[i]));


			// Something has been updated, need to rebuild data structure
			updateIndex++;

			msg_out = "data: {\ndata: \"container\": {\n";
			std::string temp;

			if(!ParamsOnly){
			for(int i=0;i<nStreams;i++)
			{
				std::string temp;
				dStreams[i]->getStream((temp));
				msg_out += "data: \""+dStreams[i]->getName()+"\": "
						+"[" + (temp) + "]" +  ",\n";
			}

			char tempc[1024];

			double t1=0,t2=0,t3=0,t4=1;

//			t1 = fmax(pdata->output_ctrl->push_timer->get_cummulative(),0);
//			t2 = fmax(pdata->output_ctrl->Comm_timer->get_cummulative(),0);
//			t3 = fmax(pdata->output_ctrl->LOSolve_timer->get_cummulative(),0);
//			t4 = fmax(0.01*(t1+t2+t3),0);



			sprintf(tempc,"data: \"guagues\": ["
					"\"Label\", \"Value\","
					"\"Push Time\", %e,"
					"\"Comm Time\", %e,"
					"\"LOSolve Time\", %e],\n",t1/t4,t2/t4,t3/t4);

			msg_out += tempc;

			sprintf(tempc,"data: \"opicard_data\":{"
					"\"istep\": %i,"
					"\"npiccard_o_avg\": %e,"
					"\"npiccard_o_last\": %i,"
					"\"npiccars_o_total\": %i},\n",pdata->rstatus->istep,
					pdata->rstatus->npiccard_outer_avg,
					pdata->rstatus->npiccard_outer_last,
					pdata->rstatus->npiccard_outer);


			msg_out += tempc;

			msg_out += "data: \"charts\": {\n";
			for(int i=0;i<nCharts;i++)
			{

				charts[i]->getStream((temp));
				msg_out += "data: " + (temp)+  ",\n";
			}

			int lastc = msg_out.rfind(",");
			if(lastc != std::string::npos)
				msg_out.erase(lastc,1);
			}
			temp = "";
			std::stringstream params(temp);

			params << "data: \"params\": {"
					"\"outputid\": \"" << pdata->rdata->output_id  << "\","
					"\"ncells\": [" << pdata->nx << "," << pdata->ny << "," << pdata->nz << "],"
					"\"lengths\": [" << pdata->Lx << "," << pdata->Ly << "," << pdata->Lz << "],"
					"\"species_info\": [";
			for(int i=0;i<pdata->nspecies;i++)
			{
				char sp_nm[128];
				if(i == 0)
					sprintf(sp_nm,"Electrons");
				else
					sprintf(sp_nm,"Ion%i",i);

				params << "{ \"name\": \"" << sp_nm << "\",";
				params << "\"nptcls\":" << pdata->nptcls_species[i] << ",";
				params << "\"charge\":" << pdata->qspecies[i] << ",";
				params << "\"mass\":" << pdata->mspecies[i] << ",";
				params << "\"weight\":" << pdata->wspecies[i] << "}";


				if(i < pdata->nspecies-1)
					params << ",";
			}

			msg_out += "data: },\n";
			params << "] }\n";

			msg_out += params.str();




			msg_out += "data: }\ndata: }\n\n";

//				msg_out = "data: {\n"
//				"data: \"msg\": \"hello world\",\n"
//				"data: \"id\": 12345\n"
//				"data: }\n\n";

//			if(updateIndex_in >= updateIndex)
//				res= false;

		updateIndex_in = updateIndex;


	io_mutex.unlock();



	return true;
}
