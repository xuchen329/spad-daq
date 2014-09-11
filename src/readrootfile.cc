#include <stdlib.h>
#include <iostream>
#include <time.h>

//ROOT
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>

#include "readrootfile.hh"

using namespace std;

int thermo2dec_shift(int ph){
    int ret = 16;
    if(ph==255)      ret=0;
    else if(ph==127) ret=1;
    else if(ph==126) ret=2;
    else if(ph==124) ret=3;
    else if(ph==120) ret=4;
    else if(ph==112) ret=5;
    else if(ph==96) ret=6;
    else if(ph==64) ret=7;
    else if(ph==0) ret=8;
    else if(ph==1) ret=9;
    else if(ph==129) ret=10;
    else if(ph==133) ret=11;
    else if(ph==135) ret=12;
    else if(ph==151) ret=13;
    else if(ph==159) ret=14;
    else if(ph==223) ret=15;
    return ret;
}

const string get_now_as_string(){
    time_t timer;
    timer = time(NULL);
    char buf[80];
    memset(buf,0,80);
    strftime(buf,80,"%F_%H-%M.root",localtime(&timer));
    string retbuf = buf;
    return retbuf;
}

//decode mode 0 data
int rawdata_decode_mode0(){
    gSystem->Load("libTree.so");
/*
    if(argc<2){
	cout<<"Decode raw data:./raw2sth <ifile> [outfile]"<<endl;
	return 0;
    }

    TString filename = argv[1];
    cout<<"File name is :"<<filename<<endl;
*/
    cout<<endl<<"Start decoding rawdata.root"<<endl;
    TFile* fpt = new TFile("./data/rawdata.root","read");
    TTree* tr = (TTree*)fpt->Get("data");

    unsigned char buf[1024];
    tr->SetBranchAddress("buf",buf);

    const int nentries  = tr->GetEntries();
    cout<<"Total number of Entries: "<<nentries<<endl;

    const unsigned char fullbuf[4] ={0xff,0xff,0xff,0xff};

    TString ofilename = "./data/Meas0_";
    ofilename+=get_now_as_string();
    //ofilename+=".root";
/*
    if(argc>2){
	ofilename = argv[2];
    }
*/
    TFile* fout= new TFile(ofilename,"recreate");
    TTree* otr = new TTree("data","--decoded dSiPM data--");
    
/*
#define N_ROW_PIXEL 26
    int numofdata1=N_ROW_PIXEL*9*18*2;
    int numofdata2=N_ROW_PIXEL*9*18+6;
    int numofdata3=N_ROW_PIXEL*9*18*2+6+48*9*2+2;
*/  
    int linesum=0;

    int cluster[9*18][26] = {};
    int tdc[9*48+1] = {};
    unsigned char conf[24];

    otr->Branch("cluster",&cluster,"cluster[162][26]/I");
    otr->Branch("tdc",&tdc,"tdc[433]/I");
    otr->Branch("conf",&conf,"conf[24]/b");
    
    int line_idx    = 0;
    bool fb = 1;//firstbyte
    int event_idx =0;
    int coarse_cnt = 0;
    int fb_ph = 0;
    cout<<"Decoding MODE 0 data"<<endl;
    for(int ientry=0;ientry<nentries;ientry++){
	tr->GetEntry(ientry);
	for(int idx=0;idx<1024/4;idx++){ //4 byte (32 bit) a group
	    if (memcmp(buf+idx*4,fullbuf,4)==0){ //all ones, should not happen
		cout<<"!!!"<<endl;
		continue;
	    }
	    //mode 0
	    if(line_idx>=26*9*18+6 && line_idx<26*9*18+6+48*9+1){//TDC
		if(fb){//first byte has coarse c*ount and 3 phase bits
//coarse count is 12 bits
		    coarse_cnt = ((int)((buf[4*idx+2])&(0xF)))<<8; //bit 3-0 from byte 2
		    coarse_cnt+=(int)(buf[4*idx+3]); //all byte 3
//take out bit 6,5,4 from byte2
		    fb_ph = (((int)(buf[4*idx+2]))&(0x70))>>4;
		}
		else {//second byte 
		    int tdc_idx = line_idx-(26*9*18+6);//calculate tdc position
		    if(tdc_idx==0){//dummy TDC
			tdc[tdc_idx]=(4095-coarse_cnt)*32;
			//cout<<"TDC["<<tdc_idx<<"] coarse count: "<<coarse_cnt<<endl;
			//cout<<"TDC["<<tdc_idx<<"] val: "<<(4095-coarse_cnt)*32<<endl;
			
		    }
		    else if(coarse_cnt!=4095){
//continue get phase
			int ph1 = fb_ph+((((int)(buf[4*idx+3]))&0x1F)<<3);//take first 5 bits for phase1
			int ph2 = (((int)(buf[4*idx+3]))>>6)+((((int)(buf[4*idx+2]))&0x3F)<<2);
			int ph1val = thermo2dec_shift(ph1);
			int ph2val = thermo2dec_shift(ph2);
			int edge = (((int)(buf[4*idx+2]))>>6)&0x1;
			int tdccode =0;
			//cout<<"coarse: "<<coarse_cnt<<" ph1: "<<ph1val<<" ph2: "<<ph2val<<" edge: "<<edge<<endl;
			if(edge==0 && ph1val>10 && ph2val<10){
			    tdccode=(4095-coarse_cnt)*16+(4095-coarse_cnt-1)*16+(16-ph2val)+(16-ph1val);
			}
			else if(edge==0 && ph1val<10 && ph2val<10){
			    tdccode=(4095-coarse_cnt-1)*16+(4095-coarse_cnt-1)*16+(16-ph2val)+(16-ph1val);
			}
			else{
			    tdccode=(4095-coarse_cnt)*32+(32-ph1val-ph2val);
			}
			tdc[tdc_idx]=tdccode;
			//cout<<"TDC["<<tdc_idx<<"] code: "<<tdccode<<endl;
		    }
		    coarse_cnt=0;
		    line_idx++;
		}
		fb=!fb;
	    }
	    if(line_idx<26*9*18){//energy
		//sum the two bytes of one line, the later byte (!fb) is shifted by 8 bits
		linesum += ((int)(buf[4*idx+3]))<<((!fb)*8);
		if(!fb){
		    //push in linesum
		    //calculate cluster number, pixels are line by line cross clusters
		    int cluster_idx = (line_idx%9)*18+line_idx/(9*26);
		    int line_nr = (line_idx/9)%26;
		    /*if(linesum){
			cout<<"Cluster: "<<cluster_idx<<" linesum: "<<linesum<<endl;
			}*/
		    cluster[cluster_idx][line_nr]=linesum;
		    linesum=0;
		    line_idx++;
		}
		fb=!fb;//fb flips between 1,0
		continue;
	    }
	    if((line_idx>=26*9*18) && (line_idx<26*9*18+6)){//conf
		conf[(line_idx-26*9*18)*4]=buf[idx*4];
		conf[(line_idx-26*9*18)*4+1]=buf[idx*4+1];
		conf[(line_idx-26*9*18)*4+2]=buf[idx*4+2];
		conf[(line_idx-26*9*18)*4+3]=buf[idx*4+3];
		fb=!fb;
		line_idx++;
	    }
	    if(line_idx==26*9*18+6+48*9+1){//finish one event
		event_idx++;
		line_idx = 0;
		fb=1;
		otr->Fill();
		memset(cluster,0,sizeof(cluster));
		memset(tdc,0,sizeof(tdc));
		memset(conf,0,sizeof(conf));
	    }
	}
    }
    cout<<"Total number of events: "<<event_idx<<endl;
    fpt->Close();
    fout->cd();
    otr->Write();
    fout->Write();
    fout->Close();
    cout<<"Decoding has finished, ROOT file: "<<ofilename<<" is generated."<<endl;
/*
    for(int i=0;i<cluster[10].size();i++){
	cout<<cluster[11].at(i)<<endl;
    }
*/    
    return 0;
}

//decode mode 1 data
int rawdata_decode_mode1(){
    gSystem->Load("libTree.so");
/*
    if(argc<2){
	cout<<"Decode raw data:./raw2sth <ifile> [outfile]"<<endl;
	return 0;
    }

    TString filename = argv[1];
    cout<<"File name is :"<<filename<<endl;
*/
    cout<<endl<<"Start decoding rawdata.root"<<endl;
    TFile* fpt = new TFile("./data/rawdata.root","read");
    TTree* tr = (TTree*)fpt->Get("data");

    unsigned char buf[1024];
    tr->SetBranchAddress("buf",buf);

    const int nentries  = tr->GetEntries();
    cout<<"Total number of Entries: "<<nentries<<endl;

    const unsigned char fullbuf[4] ={0xff,0xff,0xff,0xff};

    TString ofilename = "./data/Meas1_";
    ofilename+=get_now_as_string();
    //ofilename+=".root";
/*
    if(argc>2){
	ofilename = argv[2];
    }
*/
    TFile* fout= new TFile(ofilename,"recreate");
    TTree* otr = new TTree("data","--decoded dSiPM data--");
    
/*
#define N_ROW_PIXEL 26
    int numofdata1=N_ROW_PIXEL*9*18*2;
    int numofdata2=N_ROW_PIXEL*9*18+6;
    int numofdata3=N_ROW_PIXEL*9*18*2+6+48*9*2+2;
*/  
    int linesum=0;

    int cluster[9*18][26] = {};
    int tdc[9*48+1] = {};
    unsigned char conf[24];

    otr->Branch("cluster",&cluster,"cluster[162][26]/I");
    otr->Branch("tdc",&tdc,"tdc[433]/I");
    otr->Branch("conf",&conf,"conf[24]/b");
    
    int line_idx    = 0;
    bool fb = 1;//firstbyte
    int event_idx =0;
    int coarse_cnt = 0;
    int fb_ph = 0;
    cout<<"Decoding MODE 0 data"<<endl;
    for(int ientry=0;ientry<nentries;ientry++){
	tr->GetEntry(ientry);
	for(int idx=0;idx<1024/4;idx++){ //4 byte (32 bit) a group
	    if (memcmp(buf+idx*4,fullbuf,4)==0){ //all ones, should not happen
		cout<<"!!!"<<endl;
		continue;
	    }
	    //mode 0
	    if(line_idx>=26*9*18+6 && line_idx<26*9*18+6+48*9+1){//TDC
		if(fb){//first byte has coarse c*ount and 3 phase bits
//coarse count is 12 bits
		    coarse_cnt = ((int)((buf[4*idx+2])&(0xF)))<<8; //bit 3-0 from byte 2
		    coarse_cnt+=(int)(buf[4*idx+3]); //all byte 3
//take out bit 6,5,4 from byte2
		    fb_ph = (((int)(buf[4*idx+2]))&(0x70))>>4;
		}
		else {//second byte 
		    int tdc_idx = line_idx-(26*9*18+6);//calculate tdc position
		    if(tdc_idx==0){//dummy TDC
			tdc[tdc_idx]=(4095-coarse_cnt)*32;
			//cout<<"TDC["<<tdc_idx<<"] coarse count: "<<coarse_cnt<<endl;
			//cout<<"TDC["<<tdc_idx<<"] val: "<<(4095-coarse_cnt)*32<<endl;
			
		    }
		    else if(coarse_cnt!=4095){
//continue get phase
			int ph1 = fb_ph+((((int)(buf[4*idx+3]))&0x1F)<<3);//take first 5 bits for phase1
			int ph2 = (((int)(buf[4*idx+3]))>>6)+((((int)(buf[4*idx+2]))&0x3F)<<2);
			int ph1val = thermo2dec_shift(ph1);
			int ph2val = thermo2dec_shift(ph2);
			int edge = (((int)(buf[4*idx+2]))>>6)&0x1;
			int tdccode =0;
			//cout<<"coarse: "<<coarse_cnt<<" ph1: "<<ph1val<<" ph2: "<<ph2val<<" edge: "<<edge<<endl;
			if(edge==0 && ph1val>10 && ph2val<10){
			    tdccode=(4095-coarse_cnt)*16+(4095-coarse_cnt-1)*16+(16-ph2val)+(16-ph1val);
			}
			else if(edge==0 && ph1val<10 && ph2val<10){
			    tdccode=(4095-coarse_cnt-1)*16+(4095-coarse_cnt-1)*16+(16-ph2val)+(16-ph1val);
			}
			else{
			    tdccode=(4095-coarse_cnt)*32+(32-ph1val-ph2val);
			}
			tdc[tdc_idx]=tdccode;
			//cout<<"TDC["<<tdc_idx<<"] code: "<<tdccode<<endl;
		    }
		    coarse_cnt=0;
		    line_idx++;
		}
		fb=!fb;
	    }
	    if(line_idx<26*9*18){//energy
		//sum the two bytes of one line, the later byte (!fb) is shifted by 8 bits
		linesum += ((int)(buf[4*idx+3]))<<((!fb)*8);
		if(!fb){
		    //push in linesum
		    //calculate cluster number, pixels are line by line cross clusters
		    int cluster_idx = (line_idx%9)*18+line_idx/(9*26);
		    int line_nr = (line_idx/9)%26;
		    /*if(linesum){
			cout<<"Cluster: "<<cluster_idx<<" linesum: "<<linesum<<endl;
			}*/
		    cluster[cluster_idx][line_nr]=linesum;
		    linesum=0;
		    line_idx++;
		}
		fb=!fb;//fb flips between 1,0
		continue;
	    }
	    if((line_idx>=26*9*18) && (line_idx<26*9*18+6)){//conf
		conf[(line_idx-26*9*18)*4]=buf[idx*4];
		conf[(line_idx-26*9*18)*4+1]=buf[idx*4+1];
		conf[(line_idx-26*9*18)*4+2]=buf[idx*4+2];
		conf[(line_idx-26*9*18)*4+3]=buf[idx*4+3];
		fb=!fb;
		line_idx++;
	    }
	    if(line_idx==26*9*18+6+48*9+1){//finish one event
		event_idx++;
		line_idx = 0;
		fb=1;
		otr->Fill();
		memset(cluster,0,sizeof(cluster));
		memset(tdc,0,sizeof(tdc));
		memset(conf,0,sizeof(conf));
	    }
	}
    }
    cout<<"Total number of events: "<<event_idx<<endl;
    fpt->Close();
    fout->cd();
    otr->Write();
    fout->Write();
    fout->Close();
    cout<<"Decoding has finished, ROOT file: "<<ofilename<<" is generated."<<endl;
/*
    for(int i=0;i<cluster[10].size();i++){
	cout<<cluster[11].at(i)<<endl;
    }
*/    
    return 0;
}
