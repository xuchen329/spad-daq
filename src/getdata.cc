#include <iostream>

#include "readmask.hh"
#include "readrootfile.hh"

using namespace std;

int main(int argc, char** argv){
    if(argc<2){
	cout<<"Start with:./bin/getdata <configfile>"<<endl;
	return 1;
    }
    string daqconfig = argv[1];
    cout<<"Config file: "<<daqconfig<<endl;
    int ret = rawdata_daq(daqconfig);
    if(ret) cout<<"Wrong with communicating to FPGA!"<<endl;
    ret = rawdata_decode_mode0();
    return ret;
}
