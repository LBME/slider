#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 20 September 2021
# Author : Rafael Borges
# Objective: Given Dubious and Resolved log files from SLIDER_VENOM with residues and their RSCC,
# return variability of sequences

# USAGE: GivenDubiousResolvedResiduesReturnSequenceVariability.py Dubious.log Resolved.log OutputFile
# OutputFile should be a new path

# RETURNS: a file with amino acid possibilities (line) by residue number (column)

from __future__ import print_function
import sys,os #,time

inputAll2  = sys.argv[1]
outfile    = sys.argv[2]

with open (inputAll2) as f: frAll2=f.readlines()

dicall={}
vrscc1=0.0
for l in frAll2[1:]:
    l=l.split()
    if len(l)>1:
        if l[3]!='BaselineMainChain':
            #print(l)
            # try:
            ch=l[0]
            resn=int(l[1])
            restype=l[2][0]
            vrscc2=float(l[3])
            #print (vrscc1,vrscc2)
            if ch not in dicall:       dicall[ch]={}
            if resn not in dicall[ch]: dicall[ch][resn]=[]
            if vrscc1<vrscc2: vrscc1=float(vrscc2)
            if restype not in dicall[ch][resn] and vrscc1-vrscc2<3.0:
                dicall[ch][resn].append(restype)
            #print (dicall[ch])
            #time.sleep(1)
            # except:
            #     pass
            # exit()
    else: vrscc1=0.0


if len(dicall)>1:
    dicall['all']={}
    for ch in dicall:
        for resn in dicall[ch]:
            if resn not in dicall['all']: dicall['all'][resn]=[]
            for restype in dicall[ch][resn]:
                if restype not in dicall['all'][resn]: dicall['all'][resn].append(restype)

with open(outfile,'w') as fo:
    for ch in dicall:
        if len(dicall)>1: fo.write('>'+ch+'\n')
        variabN=0
        for resn in dicall[ch]:
            if len(dicall[ch][resn])>variabN: variabN=len(dicall[ch][resn])
        for i in range(variabN):
            for resn in sorted(dicall[ch]):
                try:    fo.write(dicall[ch][resn][i])
                except: fo.write(' ')
            fo.write('\n')