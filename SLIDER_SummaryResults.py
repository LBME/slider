#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import range
from builtins import bytes, str

import os,sys,time

input_folder = sys.argv[1]

res1_5=['1VL9','4QMC','3U8I','2BCH','2BAX','4QEM','3G8F','3FO7','4QF7','1UNE','1HN4','2PYC']
res2 = ['1L8S','3ELO','5VFH','4RFP','4YU7','3G8G','2B96','1MKT','5OWC','3QNL','5Y5E','5G3N','3JTI','1MKU','3U8D','3CXI',
        '5G3M','2B00','6G5J','1FXF','1O3W','1MKV','2ZP4','4AUP','2ZP3','2AZY','4YV5','5WZO','1C74','1GH4','2ZP5','2BD1',
        '1MKS','5OW8','1IRB','1FDK','4KF3','5WZW','1LE6','5VET','1Y6O','2WG9','5WZM','2WG7','1FX9','1KVO']
res2_5=['3JQ5','3GCI','5VFM','5VFJ','1LE7','2AZZ','5WZV','5WZU','2B01','1J1A','4UY1','2QHW','4HMB','1Y6P','2B03','4FGA',
        '3Q4Y','3BJW','5WZS','2GNS','3FVJ','4DBK','1P7O','3U8B','3U8H','2WG8','3L30','4G5I','5WZT','2B04','3HSW','3FG5',
        '4O1Y']
res3 = ['5TFV','6AL3','1O2E','1B4W','4GFY','1OZY','3FVI','1A2A','1GOD','1C1J','4EIX']

amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
amino_acid_list_3L=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']



batch_struct=[res1_5,res2,res2_5,res3]
strbatch=['res1_5','res2','res2_5','res3']

for batch in batch_struct:
    print (strbatch[batch_struct.index(batch)])
    dicaa={}
    for a in amino_acid_list:
        dicaa[a]={'Tp':0,'Fp':0,'Tn':0,'Fn':0,'TogetherNeg':0,'TogetherPos':0}
    for pdb in batch:
        with open(input_folder+'/'+pdb+'/output-'+pdb+'_all.log') as f: fr=f.readlines()
        dicvar={}
        for l in fr[1:]:
            if len(l)>1:
                l=l.split()
                ch=l[0]
                resn=l[1]
                rest=l[2]
                RSCC=float(l[3])
                if ch not in dicvar: dicvar[ch]={}
                if resn not in dicvar[ch]: dicvar[ch][resn]={}
                dicvar[ch][resn][rest]=RSCC
        for ch in dicvar:
            for resn in dicvar[ch]:
                lvar=[]
                for rest in dicvar[ch][resn]:
                    lvar.append([rest,dicvar[ch][resn][rest]])
                #print (lvar)
                lvar=sorted(lvar,key=(lambda item: item[1]), reverse=True)
                #print (lvar)
                #exit()
                bestRSCC=lvar[0][1]
                lvar2=[]
                for l in lvar:
                    rest=l[0]
                    RSCC=l[1]
                    if bestRSCC-RSCC<3: lvar2.append(l)

                if len(lvar2)==1:
                    if '!' in lvar2[0][0]: dicaa[lvar[0][0][0]]['Tp']+=1
                    else:                  dicaa[lvar[0][0][0]]['Fp']+=1
                else:
                    for l in lvar2:
                        rest = l[0]
                        if '!' in rest: dicaa[rest[0]]['TogetherPos']+=1
                        else:           dicaa[rest[0]]['TogetherNeg']+=1

                for l in lvar:
                    if l not in lvar2:
                        rest = l[0]
                        if '!' in rest: dicaa[rest[0]]['Fn']+=1
                        else:           dicaa[rest[0]]['Tn']+=1

    for a in dicaa:
        s=a
        for i in ['Tp','Fp','Tn','Fn','TogetherNeg','TogetherPos']:
            s+=' '+i+' '+str(dicaa[a][i])
        print (s)