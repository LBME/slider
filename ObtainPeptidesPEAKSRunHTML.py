#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 08 April 2020
# Author : Rafael Borges
# Objetive: Obtain positive SLIDER peptides from overall run from PEAKS HTML output
# It reads a sequence and PEAKS HTML output with all values
# It returns each peptide among its -10lgP

# USAGE: ObtainPeptidesPatternLabRun.py seq.seq PEAKS.html OutputFile
# seq.seq has just the sequence without >
# PEAKS.html exactly like PEAKS html output
# OutputFile should be a new path

# RETURNS: a table with peptide PrimaryScore PPM
#          a file with each peptide with its position in sequence

from __future__ import print_function
import sys,os

import sys,os

seqinput         = sys.argv[1] #FASTA - single line
PEAKSInput       = sys.argv[2] #PEAKSInput file
RSCCInput        = sys.argv[3]
outfile          = sys.argv[4] #where files will be organized
try: exclude     = sys.argv[5].split(',')
except: exclude=False

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

outfile1=open(outfile+'.log','w')
outfile2=open(outfile+'.log2','w')

with open(seqinput) as f: seq=f.read()

with open(PEAKSInput) as f: fr2=f.read()
with open(PEAKSInput) as f: fr1=f.readlines()

fr2=fr2[fr2.index('Supporting Peptides'):]
fr2=fr2.split('\n')

dicall={}

checkk=False
check2=False
check3=False
checkPept=False
dicprotcodes={}

for l in fr1:

    #print (l[:-1])
    #get names (Description) proteins
    if l=='<div id="bar"><div id="tit"><a name="ps">Protein List</a></div></div><br>\n': check2=True
    if check2 and l=='<tbody>\n':                                                        check3=True
    if check2 and check3:
        if l.startswith('<td id="cl">'):   protcode=l[12+l[12:].index('>')+1:l.rindex('</a></td>')] #protcode
        elif l.startswith('<td id="wt">'): protdesc=l[12:-5]                                         #prot description
        #end of a protein section
        elif l=='</tr>\n':
            #check if strings in list is contained in protein description
            ccount=0
            if exclude:
                for strr in exclude:
                    if strr not in protdesc: ccount+=1
                if ccount==len(exclude): dicprotcodes[protcode]=protdesc
                else:                    print ('Excluded peptides of',protdesc)
        #end section of protein names (Description)
        elif l=='</tbody>\n': check2,check3,checkPept=False,False,True

    if checkPept:
        if l.startswith('<div id="bar"><div id="tit">') and 'Summary</div></div><br>' not in l and 'Protein List</a></div></div><br>' not in l:
            if l[28:28+l[28:].index('<')] in dicprotcodes: check4=True
            else:                                          check4=False
        if '<tr id="row">' in l:
            counter=0
            checkk=True
        if '<td id="cc">' in l:
            counter+=1
            if counter==2:
                PScore=l[12:-6] #float(l[12:-5])

        if '<td id="wt">' in l:
            l2=l[12:-6]
            if l2[1] =='.': l2=l2[2:]
            if l2[-2]=='.': l2=l2[:-2]
            l2=l2.replace('.','')
            #print (l2)
            #for r in '0123456789+.()': l2=l2.replace(r,'')
            for i in range(l2.count('(')):
                i1 = l2.index('(')
                i2 = l2.index(')')
                l2=l2[:i1]+l2[i2+1:]
                #print (l2)

            pept=l2

        if '</tr>' in l and checkk and check4:
            #print (pept,PScore)


            if pept not in dicall: dicall[pept]=PScore
            else:
                if PScore>dicall[pept]: dicall[pept]=PScore

if exclude: print (seq,'\n')

# with open(outfile,'w') as fw:
#     fw.write('#Scan Line: Unique	FileName	ScanNumber	ChargeState	PrimaryScore	DeltCN	M+H+	CalcM+H+	ZScore	BayesianScore	RedundancyAtPtnLevel	Sequence')
#     for pept in dicline:
#         fw.write(dicline[pept])


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
lvarseq=[]
for l in lallclean:
    indfinalId=l[0]
    pept=l[1]
    pscore=l[2]
    print (' '*indfinalId+pept,pscore)
    outfile1.write(' '*indfinalId+pept+' '+pscore+'\n')
    lvarseq.append(' '*indfinalId+pept)



print ('\n\n\n'+seq,'\n')
lvarres=[]
for i,res in enumerate(seq):
    lvar=[]
    for seq in lvarseq:
        if len(seq)>i:
            res2=seq[i]
            if res2!=' ' and res2 not in lvar: lvar.append(res2)
    lvarres.append(lvar)

imax=0
for l in lvarres:
    if len(l)>imax: imax=len(l)
#print (lvarres)
strvarres=''
for i in range(imax):
    #print (i)
    for l in lvarres:
        #print (l)
        if len(l)<=i: strvarres+=' '
        else:        strvarres+=l[i]
    strvarres += '\n'
print (strvarres)
outfile2.write(strvarres)

outfile1.close()
outfile2.close()

with open(outfile+'.log') as fscore: fr=fscore.readlines()

dicscore={}
for l in fr:
    score=l[-7:-1]
    try: float(score)
    except: score=score[1:]
    score=score.replace(' ','')
    #print (l)
    ri=l.rindex(' ')
    #print (l[:ri])
    for i,s in enumerate(l[:ri]):
        if s!=' ':
            if i+1 not in dicscore: dicscore[i+1]={}
            if s in dicscore[i+1]:
                if float(score)>float(dicscore[i+1][s]): dicscore[i+1][s]=score
            else: dicscore[i+1][s]=score

# for i in dicscore:
#     for ii in dicscore[i]:
#         print (i,ii,dicscore[i][ii])

with open(RSCCInput) as frCC: frCCl=frCC.readlines()

with open(outfile+'_table.log','w') as fwCC:
    fwCC.write(frCCl[0][:-1]+'\t-10lgP\n')
    for l in frCCl[1:]:
        fwCC.write(l[:-1])
        #print (l[:-1])
        if l!='\n' and 'BaselineMainChain' not in l:
            l=l.split()
            resn=int(l[1])
            restype=l[2][0]
            #print (resn,restype,dicscore[resn][restype])
            #exit()
            if resn in dicscore and restype in dicscore[resn]:
                fwCC.write('\t'+dicscore[resn][restype])
        fwCC.write('\n')
