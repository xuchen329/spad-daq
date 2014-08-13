#!/usr/bin/env python
import numpy as np
import os,sys

#require python>2.6 with numpy>1.3.0 to work
#usage: generatemask.py <one/zero/[cluster number]> <maskfilename>
#generatemask.py one mask1.txt will generate a mask file with first 4 colunms are 1s
#cluster number separated by comma: generatemask.py 5,6,7,9 mask1.txt

def allzeromask(filename):
    FULLMASK = np.array([])
    for i in range(162):
        if i is 0:
            FULLMASK = np.reshape(np.zeros(16*26),(26,16)).astype(int)
        else:
            FULLMASK = np.append(FULLMASK,np.reshape(np.zeros(16*26),(26,16)).astype(int),axis=0)
    np.savetxt(filename+"~",FULLMASK,fmt='%1d',delimiter='',newline='\n')
    destination=open(filename,"w")
    source=open(filename+"~","r")
    cnt=0
    for line in source:
        if not cnt%26:
            clusternumber=int(cnt/26)
            comment_sentense="#D"+str(clusternumber)+"\n"
            destination.write(comment_sentense)
            destination.write(line)
        else:
            destination.write(line)
        cnt+=1

def allonemask(filename):
    FULLMASK = np.array([])
    for i in range(162):
        if i is 0:
            FULLMASK = np.reshape(np.ones(16*26),(26,16)).astype(int)
        elif i<72:
            FULLMASK = np.append(FULLMASK,np.reshape(np.ones(16*26),(26,16)).astype(int),axis=0)
        else:
            FULLMASK = np.append(FULLMASK,np.reshape(np.zeros(16*26),(26,16)).astype(int),axis=0)
    np.savetxt(filename+"~",FULLMASK,fmt='%1d',delimiter='',newline='\n')
    destination=open(filename,"w")
    source=open(filename+"~","r")
    cnt=0
    for line in source:
        if not cnt%26:
            clusternumber=int(cnt/26)
            comment_sentense="#D"+str(clusternumber)+"\n"
            destination.write(comment_sentense)
            destination.write(line)
        else:
            destination.write(line)
        cnt+=1

def selectivemask(filename,enablelist):
    FULLMASK = np.array([])
    for i in range(162):
        if i is 0:
            if i in enablelist:
                myfunc=np.ones
            else:
                myfunc=np.zeros
            FULLMASK=np.reshape(myfunc(16*26),(26,16)).astype(int)
        else:
            if i in enablelist:
                myfunc=np.ones
            else:
                myfunc=np.zeros
            FULLMASK=np.append(FULLMASK,np.reshape(myfunc(16*26),(26,16)).astype(int),axis=0)
    np.savetxt(filename+"~",FULLMASK,fmt='%1d',delimiter='',newline='\n')
    destination=open(filename,"w")
    source=open(filename+"~","r")
    cnt=0
    for line in source:
        if not cnt%26:
            clusternumber=int(cnt/26)
            comment_sentense="#D"+str(clusternumber)+"\n"
            destination.write(comment_sentense)
            destination.write(line)
        else:
            destination.write(line)
        cnt+=1

if __name__=="__main__":
    if "zero" in sys.argv[1]:
        print "Generate mask file will all 0s"
        filename = sys.argv[2]
        allzeromask(filename)
    if "one" in sys.argv[1]:
        print "Generate mask file with all 1s"
        filename = sys.argv[2]
        allonemask(filename)
    else:
        mylist = sys.argv[1]
        filename= sys.argv[2]
        enablelist = tuple([int(i) for i in mylist.split(',')])
        selectivemask(filename,enablelist)
