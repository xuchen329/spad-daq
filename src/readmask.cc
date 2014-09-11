//system
#include <string.h>
#include <iostream>
#include <fstream>
#include <bitset>
#include <vector>
#include <stdio.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

//ROOT
#include <TFile.h>
#include <TTree.h>

//others
#include "ConfigFile.hh"
#include "readmask.hh"

using namespace std;

//convert mask file to bitstream
const vector<unsigned char> mask2hex(string const maskfname){
    fstream maskfp(maskfname.c_str(),fstream::in);
    vector<unsigned char> obuf,retbuf;
    string ibuf;
    if(!maskfp.is_open()){
	cout<<"Mask File: "<<maskfname<<" not found!"<<endl;
	return retbuf;
    }
    cout<<"Mask file: "<<maskfname<<endl;
/*
  masking info are streamed to obuf

  b1b0<--cluster 0, bottom line
  ...
  until cluster 161, top line
  So, in the mask file, the first line is actually the bottom line in chip
*/
    while(!maskfp.eof()){
	getline(maskfp,ibuf);
	if(ibuf.size()!=16 or ibuf[0]=='#') continue; //strict condition check
	else{ //encode mask data
	    for(size_t i=0;i<2;i++){
		obuf.push_back(static_cast<unsigned char>(bitset<8> (ibuf.substr(8*(1-i),8)).to_ulong()));
	    }
	}
    }
//    cout<<obuf.size()<<endl;
/*
  reorder masking data
  from bottom, line by line cross clusters
*/
    for(size_t cluster_row=0;cluster_row<18;cluster_row++){
	for(size_t row=0;row<26;row++){
	    for(size_t cluster_clnm=0;cluster_clnm<9;cluster_clnm++){
		retbuf.push_back(obuf.at((cluster_clnm*18+cluster_row)*(2*26)+row*2));
		retbuf.push_back(obuf.at((cluster_clnm*18+cluster_row)*(2*26)+row*2+1));
	    }
	}
    }
    obuf.clear();
    return retbuf;
}



int rawdata_daq(const string daqconfigfile){
/*
  read connection config and maskfile name
*/
    ConfigFile config(daqconfigfile);
    string FPGA_ip;
    config.readInto(FPGA_ip,"FPGA_ip",string("192.168.1.10"));
    cout<<"Connecting to :\t"<<FPGA_ip<<endl;
    int port = config.read<int>("port", 7);
    string maskfilename;
    config.readInto(maskfilename,"maskfile");
    
    //create socket to connect
    struct sockaddr_in server;
    int sock;
    sock = socket(AF_INET, SOCK_STREAM, 0);
    server.sin_family = AF_INET;

//    server.sin_port = htons(50007); 
//    server.sin_addr.s_addr = inet_addr("127.0.0.1");

//    server.sin_port = htons(7);
//    server.sin_addr.s_addr = inet_addr("192.168.1.10");
    server.sin_addr.s_addr = inet_addr(FPGA_ip.c_str());
    server.sin_port = htons(port);

    connect(sock, (struct sockaddr *)&server, sizeof(server));

    char sendbuf[1024]; //command buffer
    memset(sendbuf,0,sizeof(sendbuf));
    
    //send the first message to server
    sprintf(sendbuf,"%d\n\r",0);
    cout<<"First message sent to chip: "<<sendbuf<<endl;
    sendto(sock, sendbuf, sizeof(sendbuf), 0, (struct sockaddr *)&server, sizeof(server));


    //constructing masking data
    vector<unsigned char> confdata = mask2hex(maskfilename);
    //cout<<"Size of maskdata: "<<confdata.size()<<endl;

/*
  read chip configurations
*/
    const float CLKPERIOD  = config.read<float>("clock_period");
    const float TIMEWINDOW = config.read<float>("daq_timewindow");
    const float REALTIMEWINDOW = TIMEWINDOW + 20e-9;
    int MODE,ROWCALSEL_ACT,PREAMP,TDCIN,ENETH,ENETEST;
    int SREN,TDCTH,SWBIASTDC,SWBIASVCO,DISEXTON;
    int N_CYCLE;

    
    MODE          = config.read<int>("daq_mode");
    ROWCALSEL_ACT = 0;
    PREAMP        = config.read<int>("preamp");
    TDCIN         = 0;
    ENETH         = config.read<int>("ene_thr");
    ENETEST       = config.read<int>("ene_test");
    //mode 1 force enetest=0
    if(MODE==1) ENETEST=0;

    SREN          = config.read<int>("smart_reset");
    TDCTH         = config.read<int>("tdc_thr");
    SWBIASTDC     = 6;
    SWBIASVCO     = 5;
    DISEXTON      = 1;
    N_CYCLE       = config.read<int>("n_cycle");
    const float DISEXTPERIOD  = 3.0e-6+0.16e-6;
    const float SRINTERVAL    = config.read<float>("sr_interval");
    int ROTIME    = config.read<int>("ro_time"); //readout time in us

    //manually construct the configuration stream
    confdata.push_back(static_cast<unsigned long>(htons((bitset<16> (TIMEWINDOW/CLKPERIOD)).to_ulong())>>8));
    confdata.push_back(static_cast<unsigned long>(htons((bitset<16> (TIMEWINDOW/CLKPERIOD)).to_ulong())));
    confdata.push_back(static_cast<unsigned char>(MODE));
    confdata.push_back(static_cast<unsigned char>(ROWCALSEL_ACT));
    confdata.push_back(static_cast<unsigned char>(PREAMP));
    confdata.push_back(static_cast<unsigned char>(TDCIN));
    confdata.push_back(static_cast<unsigned char>((bitset<16> (ENETH+(ENETEST<<15)).to_ulong())));
    confdata.push_back(static_cast<unsigned char>((bitset<16> (ENETH+(ENETEST<<15)).to_ulong())>>8));
    confdata.push_back(static_cast<unsigned char>(htons((bitset<16> (SRINTERVAL/CLKPERIOD)).to_ulong())>>8));
    confdata.push_back(static_cast<unsigned char>(htons((bitset<16> (SRINTERVAL/CLKPERIOD)).to_ulong())));
    confdata.push_back(static_cast<unsigned char>(SREN));
    confdata.push_back(static_cast<unsigned char>(TDCTH));
    confdata.push_back(static_cast<unsigned char>(SWBIASTDC));
    confdata.push_back(static_cast<unsigned char>(SWBIASVCO));
    confdata.push_back(static_cast<unsigned char>(htons(bitset<16> (DISEXTPERIOD/CLKPERIOD+(DISEXTON<<15)).to_ulong())>>8));
    confdata.push_back(static_cast<unsigned char>(htons(bitset<16> (DISEXTPERIOD/CLKPERIOD+(DISEXTON<<15)).to_ulong())));
    confdata.push_back(static_cast<unsigned char>(0));
    confdata.push_back(static_cast<unsigned char>(ROTIME));

    //send everything in the chank of 1024
    cout<<"Start sending masking and config data to FPGA ..."<<endl;
    for(size_t idx=0;idx<9;idx++){
	memcpy(sendbuf, &confdata[idx*1024],1024);
	int n = write(sock,sendbuf,sizeof(sendbuf));
    }

    //preparing for readout
#define N_ROW_PIXEL 26
#define N_COL_CLUSTER 9
#define N_ROW_CLUSTER 18
#define UNIT_DATA_BYTE 2
#define TDC_PER_CLUSTER 48
    
    int numofdata_mode0 = N_ROW_PIXEL*N_COL_CLUSTER*N_ROW_CLUSTER*UNIT_DATA_BYTE+6+TDC_PER_CLUSTER*N_COL_CLUSTER*UNIT_DATA_BYTE+2; //6 is config info, this is 256 events in mode 0
    int numofdata_mode1 = 1+9+4+48*9*2+2; //some how...

    int maxcyclenum = numofdata_mode0*N_CYCLE; //<--- total number of cycles
    
    unsigned char buf1024[1024];
    
    //preparing root file
    TFile* datafile = new TFile("./data/rawdata.root","recreate");
    TTree *tr = new TTree("data","--raw data from chip--");
    tr->Branch("buf",buf1024,"buf[1024]/b");
    
    cout<<"Start reading data from FPGA ..."<<endl;
    cout<<"Total number of r/o cycles: "<<maxcyclenum<<endl;
    for(int cycle=0;cycle<maxcyclenum;cycle++){ //readout cycles
	int n = read(sock, buf1024, sizeof(buf1024)); //first readout data from socket
	if(n>0) tr->Fill();
	else cout<<"Reading socket error !.."<<endl;
    }
    
    tr->Write();
    datafile->Write();
    datafile->Close();
    
    cout<<"R/O has finished, rawdata.root is generated ..."<<endl;
    return 0;
}
