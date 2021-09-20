#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 19 April 2020
# Author : Rafael Borges
# Objective: Given ConSurf query_msa.aln and SLIDER table,
# return percentage of SLIDER residue

# USAGE: ConSurf_SLIDER.py query_msa.aln SLIDER.txt OutputFile
# OutputFile should be a new path

# RETURNS: a table with SLIDER details plus % and n of 100% of given residue

from __future__ import print_function
import sys,os

import sys,os

CSinputAln  = sys.argv[1]
SLIDERi     = sys.argv[2]
outfile     = sys.argv[3]

with open(CSinputAln)  as  f: frCSA=f.readlines()
with open(SLIDERi) as      f: frSL =f.readlines()

#aalist=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

dicCSA={}
k=''

for l in frCSA:
    if l.startswith('>'):
        # print (l)
        if k=='':
            seq=''
            rightkey=l[:-1]
        else:
            dicCSA[k]=seq
            seq=''
            # print (dicCS)
            # exit()
        k=l[:-1]
    else:
        seq+=l[:-1]
dicCSA[k]=seq

dicN={}
dicResTN={}

seqfixed=dicCSA[rightkey]
c=1
for i,aa in enumerate(seqfixed):
    if aa!='-':
        for k in dicCSA:
            #print ( dicCS[k])
            restype=dicCSA[k][i]
            if restype!='-':
                if c not in dicN: dicN[c]=1
                else:             dicN[c]+=1

                if c not in dicResTN:          dicResTN[c]={}
                if restype not in dicResTN[c]: dicResTN[c][restype]=1
                else:                          dicResTN[c][restype]+=1
        c+=1

#print (dicN)

with open(outfile,'w') as fw:
    fw.write(frSL[0][:-1]+'\t#Struct\t%\n')
    for l in frSL[1:]:
        if l.startswith('A') or l.startswith('B') or l.startswith('C') and 'MainCh' not in l:
            l2      = l.split()
            resn    = l2[1]
            restype = l2[2][0]
            fw.write(l[:-1]+'\t')
            try:    countrestype=dicResTN[int(resn)][restype]
            except: countrestype=0
            total=dicN[int(resn)]
            perc='%.1f'%(100.0*countrestype/total)
            fw.write(str(countrestype)+'\t')
            fw.write(perc+'\n')
        else:
            fw.write(l)
