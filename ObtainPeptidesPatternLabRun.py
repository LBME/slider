#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 08 April 2020
# Author : Rafael Borges
# Objetive: Obtain positive SLIDER peptides from overall run from PATTERNLAB
# It reads a sequence and patternlab fasta file with all values
# It returns each peptide among its PrimaryScore and PPM

# USAGE: ObtainPeptidesPatternLabRun.py seq.seq fasta.fasta OutputFile
# seq.seq has just the sequence without >
# fasta.fasta exactly like PatternLab output
# OutputFile should be a new path

# RETURNS: a table with peptide PrimaryScore PPM
#          a file with each peptide with its position in sequence

from __future__ import print_function
import sys,os

import sys,os

seqinput         = sys.argv[1] #FASTA - single line
PatternLabInput  = sys.argv[2] #PatterLab file
outfile          = sys.argv[3] #where files will be organized

def GivenTwoSeqsReturnIndexMatch (seq1,seq2):
    if len(seq1)>len(seq2): seqA,seqB = seq1 , seq2
    else:                   seqA,seqB = seq2 , seq1
    finalId=0
    for iA in range(len(seqA)-len(seqB)+1):
        cid=0
        # print seqA[iA:iA+len(seqB)]
        # print seqB
        for iB in range(len(seqB)):
            cA=seqA[iA+iB]
            cB=seqB[iB]
            if cA==cB: cid+=1
        # print cid
        if cid>finalId:
            finalId=cid
            indfinalId=iA
    return finalId,indfinalId

with open(seqinput) as f: seq=f.read()

with open(PatternLabInput) as f: fr=f.readlines()

dicall={}
dicline={}

checkk=False
for l in fr:
    if l.startswith('	Locus:'):
        if not l.startswith('	Locus:contaminant'): checkk=True
        else:                                        checkk=False
    else:
        if checkk and len(l)>40:
            l2=l.split()
            pept = l2[-1]
            ipept=pept.index('.')+1
            fpept=pept.rindex('.')
            pept=pept[ipept:fpept]
            pept=pept.replace('(+15.994900)','')
            PScore=float(l2[-7])
            #print (pept,PScore)
            if pept not in dicall: dicall[pept]=PScore
            if pept not in dicline: dicline[pept]=l
            else:
                if PScore>dicall[pept]: dicall[pept]=PScore
                if PScore>dicall[pept]: dicline[pept]=l


print (seq,'\n')

with open(outfile,'w') as fw:
    fw.write('#Scan Line: Unique	FileName	ScanNumber	ChargeState	PrimaryScore	DeltCN	M+H+	CalcM+H+	ZScore	BayesianScore	RedundancyAtPtnLevel	Sequence')
    for pept in dicline:
        fw.write(dicline[pept])


lall=[]
for pept in dicall:
    finalId,indfinalId=GivenTwoSeqsReturnIndexMatch(seq,pept)
    lall.append([indfinalId,pept,dicall[pept]])

#print (lall)
lall=sorted(lall, key=lambda x: x[2], reverse=True)
#print (lall)

lallclean=[]
for i,l in enumerate(lall):
    insert = True
    pept1 = l[1]
    for i2,l2 in enumerate(lall[:i]):
        pept2=l2[1]
        #print (i,pept1,i2,pept2)
        if pept1 in pept2:
            insert=False
            #print ('Insert False',pept1,'in',pept2)
    if insert and pept1 not in lallclean: lallclean.append(l)
    #if i==3: exit()

lallclean=sorted(lallclean, key=lambda x: x[0])#, reverse=True)
for l in lallclean:
    indfinalId=l[0]
    pept=l[1]
    pscore=l[2]
    print (' '*indfinalId+pept,pscore)


