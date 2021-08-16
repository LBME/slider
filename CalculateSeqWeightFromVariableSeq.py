#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import range
from builtins import bytes, str

import os,sys,time
import multiprocessing
import collections
import numpy
import datetime
import itertools

fseqvar = sys.argv[1]
fMwExp  = sys.argv[2]
NumbSeq = int(sys.argv[3])
fout    = sys.argv[4]

amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
amino_acid_list_3L=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
dicAAweightMonoIso={'A':71.037114 ,'R':156.101111,'N':114.042927,'D':115.026943,'C':103.009185,'E':129.042593,
                    'Q':128.058578,'G':57.021464 ,'H':137.058912,'I':113.084064,'L':113.084064,'K':128.094963,
                    'M':131.040485,'F':147.068414,'P':97.052764 ,'S':87.032028 ,'T':101.047679,'U':150.95363 ,
                    'W':186.079313,'Y':163.06332 ,'V':99.068414}
#extracted from: http://www.matrixscience.com/help/aa_help.html

fo=open(fout,'w')
print ('\nVaried sequence being evaluated:')
fo.write('\nVaried sequence being evaluated:\n')

seqvar=[]
lNumbSeq=len(str(NumbSeq))

with open(fseqvar) as fr: frl=fr.readlines()
for l in frl:
    print (' '*lNumbSeq+'  '+l[:-1])
    fo.write(' '*lNumbSeq+'  '+l)

lMwExp=[]
with open(fMwExp) as fr: frl2=fr.readlines()
for Mw in frl2:
    Mw=Mw.replace(' ','').replace('\n','')
    lMwExp.append(float(Mw))

for i in range(len(frl[0])):
    laa=[]
    for ii in range(len(frl)):
        if len(frl[ii])>i:
            aa=frl[ii][i]
            if aa!=' ' and aa!='\n': laa.append(aa)
    if laa!=[]: seqvar.append(laa)

#print (seqvar)
seqvar2=[]

seqpos=1
for v in seqvar:
    seqpos=seqpos*len(v)

print ('\nNumber of possible sequences is:',seqpos)
fo.write('\nNumber of possible sequences is: '+str(seqpos)+'\n')

lseqvar=[]
counter=0
checkk=True
diffXtalMSMw=1000
for EachPossib in itertools.product(*seqvar):
    counter+=1
    #if counter%1000000==0:
    #if str(counter)[-6:]=='000000':
    #if checkk and str(counter)[-5:] == '00000':
    if str(counter)[-5:] == '00000':
        print(seqpos - counter, 'to go.')
        lseqvar = sorted(lseqvar, key=lambda x: x[3], reverse=False)[:NumbSeq + 1]
        diffXtalMSMw=lseqvar[NumbSeq-1][3]
    #     checkk=False
    # if not checkk and str(counter)[-6:] == '000000':
    #     print(seqpos - counter, 'to go.')
    #     lseqvar = sorted(lseqvar, key=lambda x: x[3], reverse=False)[:NumbSeq + 1]
    #     diffXtalMSMw=lseqvar[NumbSeq][3]

    #SMw=0
    SMw=18.012245 # Additional H in N-term and OH in C-terminal, correspondent to a water
                  # http://www.iapws.org/faq1/isotope.html 18.015268
                  # http://nuclearmasses.org/resources_folder/Wang_2017_Chinese_Phys_C_41_030003.pdf
                  # O   15.994915
                  # H    1.008665
                  # H2O 18.012245

    # for res in EachPossib:
    #     SMw+=dicAAweightMonoIso[res]
    EachPossib2=collections.Counter(EachPossib)
    for res in EachPossib2:
        SMw+=EachPossib2[res]*dicAAweightMonoIso[res]

    #Calculating different Mw from possible oxidized methionines
    nMet=EachPossib2['M']
    lSMw=[]
    for i in range(nMet+1):
        lSMw.append( SMw+(i*15.9949) ) #https://www.ionsource.com/Card/MetOx/metox.htm

    # diffXtalMSMw=1000
    for MSMw in lMwExp:
        for i,SMwvar in enumerate(lSMw):
            diffXtalMSMwVar = MSMw - SMwvar
            if diffXtalMSMwVar<diffXtalMSMw:
                lseqvar.append((''.join(EachPossib), SMwvar, diffXtalMSMwVar, abs(diffXtalMSMwVar), MSMw, i))
    #         diff=abs(MSMw-SMwvar)
    #         if diff<diffXtalMSMw:
    #             diffXtalMSMw = diff
    #             MwExp=MSMw
    #             MwTheo=SMwvar
    #             iM=i
    #
    # lseqvar.append((''.join(EachPossib), MwTheo, MwExp - MwTheo, abs(MwExp - MwTheo), MwExp, iM))

lseqvar = sorted(lseqvar, key=lambda x: x[3], reverse=False)[:NumbSeq + 1]

print ('\nCalculated Mass of sequence:')
fo.write('\nCalculated Mass of sequence:\n')

for i in range(NumbSeq):
    print ( (lNumbSeq-len(str(i+1)))*' ' +str(i+1)+')',lseqvar[i][0],'%.2f'%(lseqvar[i][1]),'%.2f'%(lseqvar[i][4]) ,'%.3f'%(lseqvar[i][2]),lseqvar[i][5])
    fo.write( (lNumbSeq-len(str(i+1)))*' ' +str(i+1)+') '+str(lseqvar[i][0])+' '+'%.2f'%(lseqvar[i][1])+' '+'%.2f'%(lseqvar[i][4]) +' '+'%.3f'%(lseqvar[i][2])+ ' '+str(lseqvar[i][5])+'\n')

fo.close()