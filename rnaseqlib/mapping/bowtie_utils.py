##
## Utilities for working with Bowtie.
##
import os
import sys
import time

import rnaseqlib

def convertBowtieToJxns(inpath, outpath, jxns, jxnmod, compress=False):
  
    if isinstance(compress,str):
        compress = eval(compress)
 
    if compress==True: 
        out = gzip.open(outpath,'wb')
    else:
        out = open(outpath,'w')

    if inpath.endswith(".gz"):
        infile = gzip.open(inpath)
    else:
        infile = open(inpath)
    
    storedjxns = {}
    print "Converting Bowtie to jxns: ", inpath, outpath
    for line in infile:
        if not line.startswith("@"):
            if 'Jxn' in line:
                vals = line.strip().split("\t")
                name = vals[0]
                coord = int(vals[3])
                jxnnum = int(coord)/jxnmod
                jxnpos = int(coord)%jxnmod
                key = jxnnum/mod
                jxn = jxns[key][jxnnum]
                if jxn not in storedjxns:
                    storedjxns[jxn] = {}
                try:
                    #storedjxns[jxn].append(name+":"+jxnpos)
                    storedjxns[jxn][jxnpos] += 1
                except:
                    #storedjxns[jxn] = [name+":"+jxnpos]
                    storedjxns[jxn][jxnpos] = 1
 
    for jxn in storedjxns:
        #out.write(jxn+"\t"+"\t".join(storedjxns[jxn])+"\n")
        data = [[x, storedjxns[jxn][x]] for x in storedjxns[jxn]]
        data.sort(key=operator.itemgetter(0))
        out.write(jxn+"\t"+";".join([":".join(map(str,x)) \
            for x in data])+"\n")
    out.close() 
