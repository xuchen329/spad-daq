#include <stdlib.h>
#include <iostream>
#include <vector>
#include <bitset>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TRint.h>
#include <TH1I.h>
#include <TCanvas.h>
#include <TH2I.h>
#include <TStyle.h>

using namespace std;

int main(int argc, char** argv){
    gSystem->Load("libTree.so");

    if(argc<2){
	cout<<"Plot data:./plotdata <infile>"<<endl;
	return 0;
    }
    
    TString infilename = argv[1];
    TFile* fin = new TFile(infilename,"read");
    TTree* tree = (TTree*)fin->Get("data");

    bool map_mode =0;
    Int_t CLUSTER_ID = 10;
    if(argc>2){
      CLUSTER_ID = atoi(argv[2]);
    }
    else map_mode=1;


    TRint *theApp = new TRint("app",&argc,argv);

    int cluster[9*18][26] = {};
    int tdc[9*48+1] = {};
    unsigned char conf[24] = {};

    tree->SetBranchAddress("cluster",&cluster);
    tree->SetBranchAddress("tdc",&tdc);
    tree->SetBranchAddress("conf",&conf);

    const int nentry = tree->GetEntries();
    cout<<"Total number of events: "<<nentry<<endl;

    vector<int> activated_cluster;
    activated_cluster.clear();

    Int_t COL=CLUSTER_ID/18;
    TH1I *histene = new TH1I("energy",Form("Energy_Cluster %d;Pixel;Counts",CLUSTER_ID),416,1,416);
    TH1I *histtdc = new TH1I("tdc",Form("TDC_Group %d;Code;Counts",COL),20000,1,20000);
    TH1I *histntdc = new TH1I("ntdc",Form("nTDC_fired_Group %d;nTDC;Counts",COL),48,1,48);
    TH2I *histmap = new TH2I("map","Pixel firing map;Cluster_ID;Pixel_fired",72,0,71,416,1,416);

    if(!map_mode){
      for(int ientry=0;ientry<nentry;ientry++){
	tree->GetEntry(ientry);
	Int_t nTDC=0;
	for(int idx=0;idx<9*18;idx++){
	  if(CLUSTER_ID!=idx) continue;
	  int npixelfired = 0;
	  for(int row=0;row<26;row++){
	    if(cluster[idx][row]!=0){
	      bitset<16> foo (cluster[idx][row]);
	      npixelfired+=foo.count();
	      //		    activated_cluster.push_back();
	    }
	  }
	  histene->Fill(npixelfired);
	  //   if(npixelfired!=0)
	  //cout<<"cluster: "<<idx<<" fired: "<<npixelfired<<endl;
	}
	for(int tidx=COL*48+1;tidx<=(COL+1)*48;tidx++){
	  if(tdc[tidx]!=0){
	    //cout<<"TDC:"<<tidx<<" value: "<<tdc[tidx]<<endl;
	    histtdc->Fill(-1*(tdc[0]-tdc[tidx]));
	    nTDC++;
	  }
	}
	/*	cout<<"Activated clusters: ";
		for(vector<int>::iterator it=activated_cluster.begin();it!=activated_cluster.end();it++){
		cout<<*it<<"\t";
	}
	cout<<endl;
*/ //	activated_cluster.clear();
//	if(ientry==1) break;
	//cout<<"Dummy TDC: "<<tdc[0]<<endl;
	//cout<<sizeof(conf)<<endl;
	histntdc->Fill(nTDC);
      }
    }
    else{
      for(int ientry=0;ientry<nentry;ientry++){
	tree->GetEntry(ientry);
	for(int idx=0;idx<9*18;idx++){
	  if(idx>71) continue;
	  int npixelfired = 0;
	  for(int row=0;row<26;row++){
	    if(cluster[idx][row]!=0){
	      bitset<16> foo (cluster[idx][row]);
	      npixelfired+=foo.count();
	    }
	  }
	  histmap->Fill(idx,npixelfired);
	}
      }
    }

    gStyle->SetOptFit(1111);
    if(!map_mode){
      TCanvas* can_tdc = new TCanvas("TDC","TDC");
      histtdc->Draw();
      TCanvas* can_ene = new TCanvas("Energy","Energy");
      histene->Draw();
      histene->Fit("gaus");
      can_ene->UseCurrentStyle();
      TCanvas* can_ntdc = new TCanvas("nTDC","nTDC");
      histntdc->Draw();
      histntdc->Fit("gaus");
      can_ntdc->UseCurrentStyle();
    }
    else{
      TCanvas* can_ene = new TCanvas("Energy","Energy");
      histmap->Draw();
    }

    theApp->Run();
    fin->Close();
    return 0;
}
