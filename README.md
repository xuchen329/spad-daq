SpadDAQ Readme

Has the same functionality as getdata_sum.c written by Shingo.
All configurations are read from etc/dsipm_daq.cfg
data from dSiPM are sorted and stored in data/Meas_<data+time>.root
es
In order for everything to work, the system needs an installation of ROOT.
In SpadDAQ folder, type

   make

everything will be compiled and two executables are generated in bin folder.

In order to get data from dSiPM, follow the same procedure as before, then in stead of using getdata_sum, go to SpadDAQ folder and run

./bin/getdata etc/dsipm_config.cfg

all configurations for dSiPM can be changed in etc/dsipm_config.cfg file. When the data acquisition is finished, a root file will be generated in data folder.

To get the pixel and tdc information from the root file, run

./bin/macroplot data/Meas_<data+time>.root 0

this will show pixel firing histogram of cluster 0, and also tdc firing histogram and all the tdc data.

In case the getdata program complains can not find libdsipm.so, do the following in bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PAT:.

this will set the environment to look for lib files in the current directory