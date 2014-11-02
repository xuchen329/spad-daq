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
#include <TMath.h>

using namespace std;

int main(int argc, char** argv){
    gSystem->Load("libTree.so");

    if(argc<2){
	cout<<"Plot data:./macrocaldcr <infile>"<<endl;
	return 0;
    }

    double Timewindow = 240e-9; //frame length
    if(argc>2){
	Timewindow = atof(argv[2]);
    }
    cout<<"Frame length is: "<<Timewindow*1e9<<" ns"<<endl;
    
    TString infilename = argv[1];
    TFile* fin = new TFile(infilename,"read");
    TTree* tree = (TTree*)fin->Get("data");
    
    TRint *theApp = new TRint("app",&argc,argv);

    int cluster[9*18][26] = {};
    //int tdc[9*48+1] = {};
    unsigned char conf[24] = {};

    tree->SetBranchAddress("cluster",&cluster);
    //tree->SetBranchAddress("tdc",&tdc);
    tree->SetBranchAddress("conf",&conf);

    const int nentry = tree->GetEntries();
    cout<<"Total number of events: "<<nentry<<endl;

    vector<int> activated_cluster;
    activated_cluster.clear();
/*
    Int_t COL=CLUSTER_ID/18;
    TH1I *histene = new TH1I("energy",Form("Energy_Cluster %d;Pixel;Counts",CLUSTER_ID),416,1,416);
    TH1I *histtdc = new TH1I("tdc",Form("TDC_Group %d;Code;Counts",COL),20000,1,20000);
    TH1I *histntdc = new TH1I("ntdc",Form("nTDC_fired_Group %d;nTDC;Counts",COL),48,1,48);
    TH2I *histmap = new TH2I("map","Pixel firing map;Cluster_ID;Pixel_fired",72,0,71,416,1,416);
*/
    int pixel_fire_counter[9*18][26][16] = {};
    for(int ientry=0;ientry<nentry;ientry++){
	tree->GetEntry(ientry);
	Int_t nTDC=0;
	for(int idx=0;idx<9*18;idx++){ //loop through all clusters
	  for(int row=0;row<26;row++){//loop throught all rows
	      if(cluster[idx][row]==0) continue;
	      bitset<16> foo (cluster[idx][row]);
	      for(std::size_t px=0; px<foo.size(); px++){
		  if(foo.test(px)) pixel_fire_counter[idx][row][px]++;
	      }
	  }
	}
    }
    //calculate dcr
    double pixel_dcr[9*18][26][16] = {};
    for(int idx=0;idx<9*18;idx++){
	for(int row=0;row<26;row++){
	    for(int px=0;px<16;px++){
		double firing_prop = (double)pixel_fire_counter[idx][row][px]/(double)nentry;
		if(abs(1.-firing_prop)<1e-10) pixel_dcr[idx][row][px] = 0;
		else pixel_dcr[idx][row][px] = 1./Timewindow*log(1./(1.-firing_prop));
		//if(pixel_dcr[idx][row][px]>0)
		//    cout<<"Cluster["<<idx<<"] "<<"Pixel["<<px+16*row<<"] DCR: "<<pixel_dcr[idx][row][px]<<endl;
	    }
	}
    }

    //calculate dcr median of every cluster
    double dcr_median[9*18] = {0.};
    double dcr_mean[9*18] = {0.};
    for(int idx=0;idx<9*18;idx++){
	dcr_median[idx] = TMath::Median(416,&(pixel_dcr[idx][0][0]));
	dcr_mean[idx] = TMath::Mean(416,&(pixel_dcr[idx][0][0]));
	if(dcr_median[idx]>0){ //print results for non empty clusters
	    cout<<"Cluster["<<idx<<"] DCR median: "<<dcr_median[idx];
	    cout<<" mean: "<<dcr_mean[idx]<<endl;
	}
    }

    theApp->Run();
    fin->Close();
    return 0;
}
