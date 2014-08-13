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

using namespace std;

int main(int argc, char** argv){
    gSystem->Load("libTree.so");

    if(argc<3){
	cout<<"Plot data:./macroplot <infile> <cluster_id>"<<endl;
	return 0;
    }
    TString infilename = argv[1];
    Int_t CLUSTER_ID = atoi(argv[2]);
    Int_t COL=CLUSTER_ID/18;
    if(CLUSTER_ID>72){
	cout<<"Cluster number: "<<CLUSTER_ID<<endl;
	cout<<"In colunm: "<<COL<<endl;
	cout<<"Out of working cluster range [0,72]"<<endl;
	return 0;
    }
    else {
	cout<<"Cluster number: "<<CLUSTER_ID<<endl;
	cout<<"In colunm: "<<COL<<endl;
    }
    TFile* fin = new TFile(infilename,"read");
    TTree* tree = (TTree*)fin->Get("data");

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
    
    //Int_t CLUSTER_ID = 10;
    // Int_t COL=CLUSTER_ID/18;
    TH1I *histene = new TH1I("energy",Form("Energy_Cluster %d;Pixel fired (per frame);Counts",CLUSTER_ID),416,1,416);
    TH1I *histtdc = new TH1I("tdc",Form("TDC_Group %d;Code;Counts",COL),20000,1,20000);
    TH1I *histntdc = new TH1I("ntdc",Form("TDC of colunm %d;TDC fired (per frame);Counts",COL),48,1,48);

    for(int ientry=0;ientry<nentry;ientry++){
	tree->GetEntry(ientry);
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
	int ntdc_cnter = 0;
	for(int tidx=COL*48+1;tidx<=(COL+1)*48;tidx++){
	    if(tdc[tidx]!=0){
		//cout<<"TDC:"<<tidx<<" value: "<<tdc[tidx]<<endl;
		histtdc->Fill(-1*(tdc[0]-tdc[tidx]));
		ntdc_cnter++;
	    }
	}
	histntdc->Fill(ntdc_cnter);
/*	cout<<"Activated clusters: ";
	for(vector<int>::iterator it=activated_cluster.begin();it!=activated_cluster.end();it++){
	    cout<<*it<<"\t";
	}
	cout<<endl;
*/ //	activated_cluster.clear();
//	if(ientry==1) break;
	//cout<<"Dummy TDC: "<<tdc[0]<<endl;
	//cout<<sizeof(conf)<<endl;
    }

    TCanvas* can_tdc = new TCanvas("TDC","TDC");
    histtdc->Draw();
    TCanvas* can_ene = new TCanvas("Energy","Energy");
    histene->Draw();
    TCanvas* can_ntdc = new TCanvas("NTDC","NTDC");
    histntdc->Draw();
    
    theApp->Run();
    fin->Close();
    return 0;
}
