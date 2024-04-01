#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 28 October 2021
# Author : Rafael Borges
# Objective: Given ConSurf query_msa.aln, PDB file and chain ID (1 letter),
# return percentage of amino acids by residue position

# USAGE: ConSurf_query_msa_export_percentages.py query_msa.aln pdb.pdb chainID OutputFile
# OutputFile should be a new path

# RETURNS: a table with percentage of amino acids by residue position

from __future__ import print_function
import sys,os
from numpy import std

CSinputAln  = sys.argv[1]
pdbi        = sys.argv[2]
chain       = sys.argv[3]
outfile     = sys.argv[4]
try:    perccutoff=float(sys.argv[5])
except:perccutoff=False



def BestAlignment2StringsReturnIndexScore (st1,st2):
    if '~' in st1:
        print ('~ in string 1, failure, correct code')
        exit()
    st11='~'*(len(st2)-1)+st1+'~'*(len(st2)-1)
    bestscore=0
    bestindex=0
    for i in range( len(st11)-len(st2)+1 ):
        scorevar = 0
        v1=st11[i:]
        for ii in range(len(st2)):
            if st2[ii]==v1[ii]: scorevar+=1
        # print st2
        # print v1
        # print scorevar
        # print '\n'
        if scorevar>bestscore:
            bestscore=scorevar
            bestindex=i
    bestindex-=len(st2)-1
    # print bestscore
    # print ' '*(-1*bestindex)+st1
    # print ' '*bestindex+st2
    # print bestindex
    return bestindex,bestscore

amino_acid_list=           ['A',  'C',  'D',  'E',  'F',  'G',  'H',  'I',  'K',  'L',  'M',  'N',  'P',  'Q',  'R',  'S',  'T',  'V',  'W',  'Y',  'M'  ]
amino_acid_list_3L=        ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR','MSE']



with open(CSinputAln)  as  f: frCSA=f.readlines()
with open(pdbi) as         f: frpdb=f.readlines()

#aalist=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

#Read sequence and residue number from CA atoms given a provided chain and pdb file
seqPDB=''
lres=[]
for l in frpdb:
    #Extracting sequence from PDB is done evaluating ATOM, CA, given chain and absent altConf or equals to A
    if l.startswith('ATOM') and l[13:15]=='CA' and l[21]==chain and l[16] in [' ','A']:
        resn=int(l[22:26])
        rest3L=l[17:20]
        rest1L=amino_acid_list[amino_acid_list_3L.index(rest3L)]
        seqPDB+=rest1L
        lres.append(resn)

# print (seqPDB,len(seqPDB))
# print (lres,len(lres))
# exit()

#Extract contiguous fragments of sequence in lseqPDB
prevres=-999
lseqPDB=[]
vseq=''
for i,s in enumerate (seqPDB):
    #print (i,s,lres[i],vseq)
    res=lres[i]
    if prevres==res-1 or prevres==-999:
        vseq+=s
    else:
        lseqPDB.append(vseq)
        vseq=''
    prevres=int(res)
lseqPDB.append(vseq)

# print (lseqPDB,len(lseqPDB[0]))
# exit()

#Extract sequences from Multiple Sequence Alignment (MSA) from ConSurf (query_msa.aln)
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

#Sequence from input of ConSurf
seqfixed=dicCSA[rightkey]
seqfixed2=seqfixed.replace('-','')

#Calculate Identity and Stdev
lseqID=[]
TotalResSeqCSA=len(seqfixed)
for vseqkey in dicCSA:
    if vseqkey!=rightkey:
        vId=0
        for i,s in enumerate(dicCSA[vseqkey]):
            if seqfixed[i]!='-' and seqfixed[i]==s: vId+=1
            lseqID.append(vId*100/TotalResSeqCSA)
Avg=sum(lseqID)/len(lseqID)
Avg='%.1f'%(Avg)
StDev=std(lseqID)
StDev='%.1f'%(StDev)
with open (outfile+'_AvgStDev.log','w') as fw: fw.write('Avg\tStDev\n'+Avg+'\t'+StDev)

# print ('PDB seq    ',seqPDB)
# print ('ConSurf seq',seqfixed2)

#Extract alignment between Fragments of sequences from PDB in respect to Consurf sequence (seqfixed2)
# print (lseqPDB)
# exit()
indexeConsurfSeqMatchedPDB=[]
for i,seq in enumerate(lseqPDB):
    vbestindex,vbestscore=BestAlignment2StringsReturnIndexScore(seqfixed2,seq)
    #print (range( vbestindex,vbestindex+len(seq)))
    indexeConsurfSeqMatchedPDB+=range(vbestindex,vbestindex+len(seq))
    print ('\nSequence in PDB contiguous Fragment '+str(i+1))
    print('PDBseq:    '+ ' '* vbestindex + seq)
    print('ConsurfSeq '+seqfixed2)
    #print (len(indexeConsurfSeqMatchedPDB))
    if len(seq)>vbestscore:
        print ('Contiguous residues in PDB did not match perfectly any fragment of sequence in Consurf.')
        # print (' '*vbestindex+seq)
        # print (seqfixed2)
        exit()
# exit()


# if seqPDB!=seqfixed2:
#     print ('Exiting program because sequences from PDB and from ConSurf input are different.')
#     exit()

# #Alignment between PDB and ConSurf sequences
# if seqPDB!=seqfixed:
# #if True:
#     from Bio import pairwise2
#     pairwise_tuple = pairwise2.align.localxx(seqPDB, seqfixed)
#     bestalignment=pairwise_tuple[0]
#     seqPDBalign  =bestalignment[0]
#     seqfixedalign=bestalignment[1]
#     #print (pairwise_tuple)
#
#     # for i,s in enumerate(pairwise_tuple):
#     #     for ii,ss in enumerate(s):
#     #         print (i,ii,ss)
#     # #print (pairwise_tuple)
#     #exit()


dicN={}
dicResTN={}

counterPDB=0
counterseqfixed2=-1
#c=1
for i,aa in enumerate(seqfixed):
    if aa!='-':
        counterseqfixed2+=1
        if counterseqfixed2 in indexeConsurfSeqMatchedPDB:
            for k in dicCSA:
                #print ( dicCS[k])
                restype=dicCSA[k][i]
                c=lres[counterPDB]
                if restype!='-':
                    if c not in dicN: dicN[c]=1
                    else:             dicN[c]+=1

                    if c not in dicResTN:          dicResTN[c]={}
                    if restype not in dicResTN[c]: dicResTN[c][restype]=1
                    else:                          dicResTN[c][restype]+=1
            counterPDB+=1

#print (dicN)

with open(outfile,'w') as fw:
    fw.write('Outputting residue number of PDB file ('+pdbi+') and percentages of observed amino acids in ConSurf output ('+CSinputAln+')')
    for c in dicResTN:
        fw.write( '\n\nResidue: '+str(c) )
        fw.write('\nResType\t#Seqs\t%')
        for restype in sorted(dicResTN[c]):
            perc='%.1f'%(100*dicResTN[c][restype]/dicN[c])
            fw.write('\n'+restype+'\t'+str(dicResTN[c][restype])+'\t'+perc)

for c in dicResTN:
    for restype in amino_acid_list:
        if restype not in dicResTN[c]: dicResTN[c][restype]=0

with open(outfile+'_Joao','w') as fw:
    #fw.write('Outputting residue number of PDB file ('+pdbi+') and percentages of all amino acids in ConSurf output ('+CSinputAln+')')
    fw.write('Res#\tResType\t#Seqs\t%')
    for c in dicResTN:
        #fw.write( '\n\nResidue: '+str(c) )
        for restype in sorted(dicResTN[c]):
            perc='%.1f'%(100*dicResTN[c][restype]/dicN[c])
            fw.write('\n'+str(c)+'\t'+restype+'\t'+str(dicResTN[c][restype])+'\t'+perc)

#print (dicResTN)
with open(outfile+'_PDB','w') as fw:
    #fw.write('Outputting residue number of PDB file ('+pdbi+') and percentages of all amino acids in ConSurf output ('+CSinputAln+')')
    fw.write('Res#\tResType\t#Seqs\t%')
    for i,c in enumerate(dicResTN):
        #fw.write( '\n\nResidue: '+str(c) )
        restype=seqPDB[i]
        perc='%.1f'%(100*dicResTN[c][restype]/dicN[c])
        fw.write('\n'+str(c)+'\t'+restype+'\t'+str(dicResTN[c][restype])+'\t'+perc)

#Write sequence variation given a cutoff
if perccutoff:
    lcutoff=[]
    countmaxres=0
    for c in dicResTN:
        vcountmaxres=0
        vaa=[]
        for restype in sorted(dicResTN[c]):
            vperc=100*dicResTN[c][restype]/dicN[c]
            if vperc>=perccutoff:
                vcountmaxres += 1
                vaa.append([restype,vperc])
        vaa=sorted(vaa,key=(lambda item: item[1]), reverse=True)
        #exit()
        vaa2=''
        for s in vaa: vaa2+=s[0]
        lcutoff.append(vaa2)
        if vcountmaxres>countmaxres: countmaxres=int(vcountmaxres)
    # print (lcutoff)
    # print (countmaxres)
    # exit()
    with open(outfile + '_cutoff.seq', 'w') as fw:
        for i in range(countmaxres):
            for s in lcutoff:
                #print (s,i)
                if len(s)>i: fw.write(s[i])
                else:         fw.write(' ')
            fw.write('\n')
