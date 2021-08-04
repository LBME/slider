#!/usr/bin/env python
# -*- coding: utf-8 -*-

#to import personal my library
#import RJB_lib
#import #

from __future__ import print_function
from builtins import range
from builtins import bytes, str

import string
import sys
import os
import subprocess
import time
import numpy
import math
from math import acos, asin, atan, sqrt, degrees, atan2
import numpy
def cosdeg (angle_degree):
    return math.cos(math.radians(angle_degree))
def sindeg (angle_degree):
    return math.sin(math.radians(angle_degree))
import itertools
#import ConfigParser
from collections import defaultdict


#sys.path.insert(0, "/cri4/rafael/Git_Scripts/git-tools/tools")
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/git-tools/tools/Brasil")
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/SLIDER/seq_slider")
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/REPO_EXTERNAL/ARCIMBOLDO/borges-arcimboldo/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/REPO_EXTERNAL/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/REPO_EXTERNAL/borges-arcimboldo/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/git-tools/tools")
#sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/git-tools/tools/Brasil")
## sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/REPO_EXTERNAL/ARCIMBOLDO/borges-arcimboldo/ARCIMBOLDO_FULL")
## sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/REPO_EXTERNAL/borges-arcimboldo/ARCIMBOLDO_FULL")
## sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/REPO_EXTERNAL/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/SLIDER/seq_slider")
#sys.path.insert(0, "/home/rborges/REPO-ARCIMBOLDO/borges-arcimboldo/ARCIMBOLDO_FULL")


#import BORGES_MATRIX
from operator import itemgetter, attrgetter, methodcaller
from termcolor import colored
import traceback
# import Grid
from collections import defaultdict
##from alixe_library import generate_fake_ins_for_shelxe , read_cell_and_sg_from_pdb
# import SELSLIB2
import datetime
from math import pi
import Bio.PDB
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as MatList

amino_acid_list=           ['A',  'C',  'D',  'E',  'F',  'G',  'H',  'I',  'K',  'L',  'M',  'N',  'P',  'Q',  'R',  'S',  'T',  'V',  'W',  'Y',  'M'  ]
amino_acid_list_3L=        ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR','MSE']
amino_acid_list_numb_atoms=[   5,    6,    8  ,  9  , 11  ,  4  , 10  ,  8  ,  9  ,  8  ,  8  ,  8  ,  7  ,  9  , 11  ,  6  ,  7  ,  7 ,  14 ,  12 , 8   ]
DicHydrophobic={'ARG':('CB','CG'),'LYS':('CB','CG','CD'),'HIS':('CB'),'ASP':('CB'),'GLU':('CB','CG'),'PRO':('CB','CG'),
                'MET':('CB','CG','CE'),'ALA':('CB'),'VAL':('CB','CG1','CG2'),'LEU':('CB','CG','CD1','CD2'),
                'ILE':('CB','CG1','CG2','CD1'),'TYR':('CB','CG','CD1','CD2','CE1','CE2'),
                'PHE':('CB','CG','CD1','CD2','CE1','CE2','CZ'),'TRP':('CG','CD2','CE3','CZ3','CH2','CZ2'),'CIS':('CB'),
                'GLN':('CB','CG'),'ASN':('CB'),'THR':('CG2'),'GLY':(),'SER':()}
#alphabet=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqstuvwxyz'


def generate_HLC_from_FOM ( mtz_input , mtz_output ):
    os.system ('chltofom -mtzin '+ mtz_input + ' -colin-phifom "/*/*/[PHIC,FOM]" -mtzout ' + mtz_output + ' -colout HL' )

def extract_labels_mtz_to_list ( file_input_mtz ):
    os.system ( 'mtzinfo ' + file_input_mtz + ' | grep LABELS > log_labels_mtz.txt' )
    labels_file=open('log_labels_mtz.txt','r')
    labels_list=labels_file.read()
    labels_file.close()
    labels_list=labels_list[:-1].split()    
    os.system('rm log_labels_mtz.txt')
    return labels_list

def extract_types_mtz_to_list ( file_input_mtz ):
    os.system('mtzinfo '+file_input_mtz+' | grep TYPES > log_types_mtz.txt')
    types_file=open('log_types_mtz.txt','r')
    types_list=types_file.read()
    types_file.close()
    types_list=types_list[:-1].split()
    os.system('rm log_types_mtz.txt')
    return types_list

def extract_unit_cell_space_group_from_mtz_to_list ( file_input_mtz ):
    os.system ( 'mtzinfo ' + file_input_mtz + ' | grep XDATA > log_xdata_mtz.txt' )
    xdata_file=open('log_xdata_mtz.txt','r')
    xdata_list=xdata_file.read()
    xdata_file.close()
    try:
        xdata_list=xdata_list[:-1].split()
        unit_cell=xdata_list[1:7]
        space_group=xdata_list[9]
        resolution=xdata_list[7:9]
        os.system('rm log_xdata_mtz.txt')
        return unit_cell,space_group
    except:
        try:
            os.system ( 'mtzdmp ' + file_input_mtz + ' | grep " * Number of Columns" -B3 > log_cell_mtz.txt' )
            cell_file=open('log_cell_mtz.txt','r')
            unit_cell=cell_file.readlines()[0].split()
            #cell_line=cell_list[0].split()
            cell_file.close()
            os.system ( 'mtzdmp ' + file_input_mtz + ' | grep " * Space group " > log_space_group_mtz.txt' )
            space_group_file=open('log_space_group_mtz.txt','r')
            space_group_list=space_group_file.read().split()
            space_group=space_group_list[-1][:-1]
            space_group_file.close()
            return unit_cell,space_group
        except:
            print ('Unfortunately, problem extracting information from both "mtzinfo/mtzdmp '+file_input_mtz+'"')
        exit()

def phs2mtz ( file_input_phs , file_input_mtz , file_output_mtz , printt=False):
    unit_cell,space_group=extract_unit_cell_space_group_from_mtz_to_list ( file_input_mtz )
    phs2mtz_job_descr=open('phs2mtz_job_instr.txt','w')
    phs2mtz_job_descr.write('CELL\t')
    unit_cell_string=''
    for item in unit_cell:
        unit_cell_string+=item+' '
    phs2mtz_job_descr.write(unit_cell_string+'\n')
    phs2mtz_job_descr.write('SYMM\t'+space_group+'\n')
    phs2mtz_job_descr.write('labout\tH K L F FOM PHI SIGF\n')
    phs2mtz_job_descr.write('CTYPOUT\tH H H F W P Q\n')
    phs2mtz_job_descr.close()
    phs2mtz_job_instr=open('phs2mtz_job_instr.txt','r')
    phs2mtz_job = subprocess.Popen([ 'f2mtz','HKLIN',file_input_phs,'HKLOUT',file_output_mtz], stdin=phs2mtz_job_instr, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    out, err = phs2mtz_job.communicate()
    if printt: print (out)
    phs2mtz_job_instr.close()
    #os.system('phs2mtz_job_instr.txt')

##http://shelx.uni-ac.gwdg.de/~tg/teaching/anl-ccp4/tutorial.pdf
##A .phs-ﬁle can be converted to .mtz-format with the CCP4-program f2mtz and a short script like the
##following:
###!/bin/bash
##f2mtz hklin $1 hklout ${1%phs}mtz << eof
##CELL 92.8 92.8 129.2 90.000 90.000 120.00
##SYMM P6122
##title My short title
##labout H K L F FOM PHI sigF
##CTYPOUT H H H F W P Q
##pname CCP4 Workshop 2010
##dname phs to mtz conversion
##END
##eof
    
def check_add_HLC_coefficients ( file_input_mtz ):
    types_list_mtz=extract_types_mtz_to_list ( file_input_mtz )
    if 'A' not in types_list_mtz:
        os.system ('mv ' + file_input_mtz + ' ' + file_input_mtz+'_bk')
        #time.sleep(5)
        generate_HLC_from_FOM ( file_input_mtz+'_bk' , file_input_mtz )

def mtz2hlc_script ( file_input_mtz , file_output_hlc ):
    labels_list=extract_labels_mtz_to_list ( file_input_mtz )
    types_list=extract_types_mtz_to_list ( file_input_mtz )
    mtz2hlc_job_descr=open('mtz2hlc_job_instr.txt','w')
    mtz2hlc_job_descr.write('LABIN FP=')
    if 'FP' in labels_list :
        mtz2hlc_job_descr.write('FP')
    elif 'F' in labels_list :
        mtz2hlc_job_descr.write('F')
    elif 'FOBS' in labels_list :
        mtz2hlc_job_descr.write('FOBS')
    elif 'F-obs' in labels_list :
        mtz2hlc_job_descr.write('F-obs')
    else:
        print ('mtz2hlc_script was unable to find FP/F from mtz file:'+file_input_mtz)
        quit()
    mtz2hlc_job_descr.write(' DUM1=')
    FOM=labels_list[types_list.index('W')]
    mtz2hlc_job_descr.write(FOM)
    mtz2hlc_job_descr.write(' DUM2=PHIC SIGFP=')
    if 'SIGFP' in labels_list :
        mtz2hlc_job_descr.write('SIGFP \\\n')
    elif 'SIGF' in labels_list :
        mtz2hlc_job_descr.write('SIGF \\\n')
    elif 'SIGFOBS' in labels_list :
        mtz2hlc_job_descr.write('SIGFOBS \\\n')
    elif 'SIGF-obs' in labels_list:
        mtz2hlc_job_descr.write('SIGF-obs')
    else:
        print ('mtz2hlc_script was unable to find SIGFP/SIGF from mtz file:'+file_input_mtz)
        quit()

    ##    HLC=[]
    ##for i, j in enumerate(['A']):
    ##    HLC.append(i)
    ##    for label in range(len(labels_list)):
    ##        if labels_list[label]=='A':
    ##            HLC.append(label)

##    mtz2hlc_job_descr.write('HLA=HL.ABCD.A HLB=HL.ABCD.B \\\n')
##    mtz2hlc_job_descr.write('HLC=HL.ABCD.C HLD=HL.ABCD.D \n')

    mtz2hlc_job_descr.write('HLA=HLA HLB=HLB \\\n')
    mtz2hlc_job_descr.write('HLC=HLC HLD=HLD \n')
    mtz2hlc_job_descr.write("OUTPUT USER '(3I4,F9.2,F6.3,F7.1,F8.2,4F9.4)' \n")
    mtz2hlc_job_descr.write('END')
    mtz2hlc_job_descr.close()
    mtz2hlc_job_instr=open('mtz2hlc_job_instr.txt','r')

    mtz2hlc_job = subprocess.Popen([ 'mtz2various','HKLIN',file_input_mtz,'HKLOUT',file_output_hlc], stdin=mtz2hlc_job_instr, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    out, err = mtz2hlc_job.communicate()
    print (out)
    mtz2hlc_job_instr.close()
    os.system('rm mtz2hlc_job_instr.txt')

    ##mtz2various HKLIN tempppp.mtz HKLOUT zzz.hlc << EOF
    ##LABIN FP=FP DUM1=FOM DUM2=PHIC SIGFP=SIGFP \
    ##HLA=HL.ABCD.A HLB=HL.ABCD.B \
    ##HLC=HL.ABCD.C HLD=HL.ABCD.D 
    ##OUTPUT USER '(3I4,F9.2,F6.3,F7.1,F8.2,4F9.4)'
    ##END
    ##EOF
    
def mtz2phi_script ( file_input_mtz , file_output_phi , mtz_f=False , mtz_sigf=False , mtz_fom=False, mtz_strf=False , count=0, printt=False): #change file to variable

    if printt: print ('Converting',file_input_mtz,'into',file_output_phi)

    mtz2phi_job_descr = open('mtz2phi_job_instr' + str(count) + '.txt', 'w')

    if mtz_f!=False and mtz_sigf!=False and mtz_fom!=False and mtz_strf!=False:
        mtz2phi_job_descr.write('LABIN FP='+mtz_f)
        mtz2phi_job_descr.write(' DUM1='+mtz_fom)
        mtz2phi_job_descr.write(' DUM2='+mtz_strf)
        mtz2phi_job_descr.write(' SIGFP='+mtz_sigf+' \n')

    else:
        labels_list=extract_labels_mtz_to_list ( file_input_mtz )
        types_list=extract_types_mtz_to_list ( file_input_mtz )

        mtz2phi_job_descr.write('LABIN FP=')
        if mtz_f==False:
            if 'FP' in labels_list :
                mtz2phi_job_descr.write('FP')
            elif 'F' in labels_list :
                mtz2phi_job_descr.write('F')
            elif 'FOBS' in labels_list :
                mtz2phi_job_descr.write('FOBS')
            elif 'F-obs' in labels_list:
                mtz2phi_job_descr.write('F-obs')
            elif 'F-obs-filtered' in labels_list:
                mtz2phi_job_descr.write('F-obs-filtered')
            else:
                print ('mtz2phi_script was unable to find F/FP/FOBS/F-obs/F-obs-filtered from mtz file:'+file_input_mtz)
                quit()
        else:
            mtz2phi_job_descr.write(mtz_f)

        mtz2phi_job_descr.write(' DUM1=')
        if mtz_fom==False:
            try:
                FOM=labels_list[types_list.index('W')]
                mtz2phi_job_descr.write(FOM)
            except:
                mtz2phi_job_descr.write('FOM')
        else: mtz2phi_job_descr.write(mtz_fom)

        mtz2phi_job_descr.write(' DUM2=')
        if mtz_fom==False:
            if 'PHIC' in labels_list :
                mtz2phi_job_descr.write('PHIC SIGFP=')
            elif 'PHIFCALC' in labels_list :
                mtz2phi_job_descr.write('PHIFCALC SIGFP=')
            elif 'PHIF-model' in labels_list :
                mtz2phi_job_descr.write('PHIF-model SIGFP=')
            else:
                print ('mtz2phi_script was unable to find PHIC/PHIFCALC from mtz file:' + file_input_mtz)
                quit()
        else: mtz2phi_job_descr.write(mtz_strf)
        if mtz_sigf==False:
            if 'SIGFP' in labels_list :
                mtz2phi_job_descr.write('SIGFP \n')
            elif 'SIGF' in labels_list :
                mtz2phi_job_descr.write('SIGF \n')
            elif 'SIGFOBS' in labels_list :
                mtz2phi_job_descr.write('SIGFOBS \n')
            elif 'SIGF-obs' in labels_list:
                mtz2phi_job_descr.write('SIGF-obs \n')
            elif 'SIGF-obs-filtered' in labels_list:
                mtz2phi_job_descr.write('SIGF-obs-filtered \n')
            else:
                print ('mtz2phi_script was unable to find SIGFP/SIGF/SIGFOBS/SIGF-obs/SIGF-obs-filtered from mtz file:'+file_input_mtz)
                quit()
        else:
            mtz2phi_job_descr.write(mtz_sigf+' \n')

    mtz2phi_job_descr.write("OUTPUT USER '(3I4,F9.2,F6.3,F7.1,F8.2)' \n")
    mtz2phi_job_descr.write('END')
    mtz2phi_job_descr.close()
    mtz2phi_job_instr=open('mtz2phi_job_instr'+str(count)+'.txt','r')
    # print 'mtz2various','HKLIN',file_input_mtz,'HKLOUT',file_output_phi
    # exit()
    mtz2phi_job = subprocess.Popen([ 'mtz2various','HKLIN',file_input_mtz,'HKLOUT',file_output_phi], stdin=mtz2phi_job_instr, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    out, err = mtz2phi_job.communicate()
    #print out
    mtz2phi_job_instr.close()
    if not os.path.isfile(file_output_phi):
        print ('\n\n\n\n\n\nFAILURE IN FUNCTION RJB_LIB.mtz2phi_script CONVERTING',file_input_mtz,'TO',file_output_phi)
        print ('Command used:')
        print ('mtz2various', 'HKLIN', file_input_mtz, 'HKLOUT', file_output_phi, mtz2phi_job_instr)
        print ('Instruction file:,mtz2phi_job_instr'+str(count)+'.txt')
        print ('out:')
        print (out)
        print ('err')
        print (err)
        mtz2phi_job_instr.close()
        quit()
    else:
        mtz2phi_job_instr.close()
        os.system('rm mtz2phi_job_instr'+str(count)+'.txt')



##mtz2various HKLIN ${1} HKLOUT $output.phi > /dev/null << EOF
##LABIN FP=F DUM1=FOM DUM2=PHIC SIGFP=SIGF
##OUTPUT USER '(3I4,F9.2,F6.3,F7.1,F8.2,4F9.4)'
##END
##EOF

def phs2mtz_script ( file_input_phs , instructions , file_output_mtz ):
    if instructions.endswith('.mtz'):
        labels_list=extract_labels_mtz_to_list ( file_input_mtz )
    types_list=extract_types_mtz_to_list ( file_input_mtz )
    mtz2phi_job_descr=open('mtz2phi_job_instr.txt','w')
    mtz2phi_job_descr.write('LABIN FP=')
    if 'FP' in labels_list :
        mtz2phi_job_descr.write('FP')
    elif 'F' in labels_list :
        mtz2phi_job_descr.write('F')
    elif 'FOBS' in labels_list :
        mtz2phi_job_descr.write('FOBS')
    else:
        print ('mtz2phi_script was unable to find F/FP/FOBS from mtz file:'+file_input_mtz)
        quit()
    mtz2phi_job_descr.write(' DUM1=')
    FOM=labels_list[types_list.index('W')]
    mtz2phi_job_descr.write(FOM)
    mtz2phi_job_descr.write(' DUM2=PHIC SIGFP=')

    if 'SIGFP' in labels_list :
        mtz2phi_job_descr.write('SIGFP \n')
    elif 'SIGF' in labels_list :
        mtz2phi_job_descr.write('SIGF \n')
    elif 'SIGFOBS' in labels_list :
        mtz2phi_job_descr.write('SIGFOBS \n')
    else:
        print ('mtz2phi_script was unable to find SIGFP/SIGF from mtz file:'+file_input_mtz )
        quit()

    mtz2phi_job_descr.write("OUTPUT USER '(3I4,F9.2,F6.3,F7.1,F8.2)' \n")
    mtz2phi_job_descr.write('END')
    mtz2phi_job_descr.close()
    mtz2phi_job_instr=open('mtz2phi_job_instr.txt','r')

    mtz2phi_job = subprocess.Popen([ 'mtz2various','HKLIN',file_input_mtz,'HKLOUT',file_output_phi], stdin=mtz2phi_job_instr, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    out, err = mtz2phi_job.communicate()
    print (out)
    mtz2phi_job_instr.close()
    #os.system('rm mtz2phi_job_instr.txt')


def extract_table_for_dictionary ( input_file , label_index_to_write_key_in_output_dictionary ):
    input_fileee=open( input_file , 'r' )
    input_file_list=input_fileee.readlines()
    input_fileee.close()
    labels=input_file_list[0].split()
    dictio_all={}
    for line in input_file_list[1:]:
        dictio_var={}
        line=line.split()
        for i in range(len(labels)):
            dictio_var[labels[i]]=line[i]
        dictio_all[line[label_index_to_write_key_in_output_dictionary]]=dictio_var
    return dictio_all

def extract_table_for_dictionary_of_list ( input_file , label_index_to_write_key_in_output_dictionary ):
    input_fileee=open( input_file , 'r' )
    input_file_list=input_fileee.readlines()
    input_fileee.close()
    labels=input_file_list[0].split()
    dictio_all={}
    for line in input_file_list[1:]:
        dictio_var={}
        line=line.split()
        for i in range(len(labels)):
            dictio_var[labels[i]]=[]
        dictio_all[line[label_index_to_write_key_in_output_dictionary]]=dictio_var
##        dictionary created with number of residues and lists
    for line in input_file_list[1:]:
        line=line.split()
        for i in range(len(labels)):
            dictio_all[line[label_index_to_write_key_in_output_dictionary]][labels[i]].append(line[i])
    return dictio_all


def extract_table_to_list_of_dict_with_first_line_as_key ( input_file ): #script generated in Jul7,2016 based on extract_table_for_dictionary to be more general, I will create in list instead of chosen a dic, thus it has no down side!
    input_fileee=open( input_file , 'r' )
    input_file_list=input_fileee.readlines()
    input_fileee.close()
    labels=input_file_list[0].split()
    list_dic_all=[]
    for line in input_file_list[1:]:
        dictio_var={}
        line=line.split()
        for i in range(len(labels)):
            dictio_var[labels[i]]=line[i]
        list_dic_all.append(dictio_var)
    return list_dic_all


    
def from_list_dic_generate_dic_dic_lists_chosen_index ( list_dic_all , chosen_key  ):
    dic_dic_list_all={}
    if not chosen_key in list_dic_all[0]:
        print ('Chosen key does not exist in given list of dics. Failure in RJB_lib.from_list_dic_generate_dic_dic_lists_chosen_index.')
        print ('\nExisting keys are:')
        for key in sorted(list_dic_all[0]):
            print (key )
        exit()
    for dic in list_dic_all:
        try:
            dic_dic_list_all[dic[chosen_key]]
        except:
            dic_dic_list_all[dic[chosen_key]]=defaultdict(list)
        for key in dic:
            if not key==chosen_key:
                dic_dic_list_all [dic[chosen_key]] [key] .append( dic[key] )
    return dic_dic_list_all


def given_dic_of_dic_of_lists_return_list_dic_mean_stdev ( dic_dic_list_all , list_reject_keys ):
    list_dic_summ=[]
    for chain in sorted(dic_dic_list_all):
        dic_var={}
        dic_var['chain']=chain
        dic=dic_dic_list_all[chain]
        for key in dic:
            if not key in list_reject_keys :
                lista=dic[key]
                while 'n/a' in lista: lista.remove('n/a')
                if len( lista)==0:
                    mean='n/a'
                    std='n/a'
                else:
                    lista=map(float,lista)
                    mean=numpy.mean(lista)
                    std=numpy.std(lista)
                    mean=mean_dev_return_good_value ( mean , 'mean' )
                    std=mean_dev_return_good_value ( std , 'std' )
                dic_var[key+'_m']=mean
                dic_var[key+'_sd']=std
        list_dic_summ.append(dic_var)
    return list_dic_summ


def given_out_edstats_return_list_mean_stdev_by_chain (input_file):
    list_dic_all=extract_table_to_list_of_dict_with_first_line_as_key ( input_file )
    dic_dic_list_all=from_list_dic_generate_dic_dic_lists_chosen_index ( list_dic_all , 'CI'  )
    list_dic_summ=given_dic_of_dic_of_lists_return_list_dic_mean_stdev ( dic_dic_list_all , ['RT','MN','CP','NR'] )
    return list_dic_summ
    
def generate_deviation_by_resid ( out_input_file , amino_acid_list_3L , out_output_file ): #corrected in 2 October 2015
    out_fileee=open( out_input_file , 'r' )
    out_file_list=out_fileee.readlines()
##    print out_file_list
    try:
        columns=out_file_list[0].split() #column of first line
    except:
        print ("\n\n*******" + out_input_file, "should contain more values:",out_file_list , "*******\n\n")
        return None
    #part that writes in a list, columns that should be saved in new table file
    columns_chosen=[] #columns to be inserted in 
    columns_chosen.append(columns[2])
    columns_chosen.append(columns[3])
    for label in columns:
        if label.endswith('s') or label.endswith('a') or label.endswith('m'):
            columns_chosen.append(label)
    #
    dictio_table={}
##    #generate dictionary (containing all) of dictionary of possible amino_acids containing empty lists
##    for amino_acid in amino_acid_list:
##        dictio_var={}
##        for label2 in columns_chosen:
##            dictio_var[label2]=[]
##        dictio_table[amino_acid]=dictio_var
    for amino_acid in amino_acid_list_3L :
##        if amino_acid!='GLY' and amino_acid!='ALA':
        dic_amino_acid={}
        for label in columns_chosen:
            var_list=[]
            for line in out_file_list:
                #print line
                line_list=line.split()
                #print line_list
                if line_list[1]==amino_acid:
                    for index_list_line in range(len(line_list) ) :
                        #print out_file_list[0][index_list_line]
                        if columns[index_list_line]==label:
                            #if (amino_acid!='GLY' and amino_acid!='ALA' and label.endswith("s") ) or not label.endswith("s"):
                            if line_list[index_list_line]!="n/a":
                                var_list.append(line_list[index_list_line])
                            else:
                                var_list.append("NaN")
            dic_amino_acid[label]=var_list
            dictio_table[amino_acid]=dic_amino_acid
    #for item in dictio_table:
   # print dictio_table
    out_output_file2=open(out_output_file,'w')
    out_output_file2.write('RT\t#models\tCI\tRN\t')
    for item in columns_chosen[2:]:
        out_output_file2.write(item+'(mean)\t'+item+'(dev)\t')
    out_output_file2.write('\n')
    for amino_acid in dictio_table:
        if not len(dictio_table[amino_acid]['CI'])==0:
            out_output_file2.write(amino_acid+'\t')
            numb_models=str(len(dictio_table[amino_acid]['CI']))
            out_output_file2.write(numb_models+'\t')
            for key in columns_chosen:
                if key=='CI' or key=='RN':
                    out_output_file2.write(dictio_table[amino_acid][key][0]+'\t')
                else:
                    if key=='NPm':
                        list=[int(i) for i in dictio_table[amino_acid][key]]
                    else:
                        list=[float(i) for i in dictio_table[amino_acid][key]]
                    mean=numpy.mean(list)
                    mean=mean_dev_return_good_value ( mean , 'mean' )
                    std=numpy.std(list)
                    std=mean_dev_return_good_value ( std , 'std' )
                    out_output_file2.write(mean+'\t')
                    out_output_file2.write(std+'\t')
            out_output_file2.write('\n')
    out_fileee.close()
    out_output_file2.close()



def generate_deviation_by_resid_old ( out_input_file , amino_acid_list_3L , out_output_file ):
    out_fileee=open( out_input_file , 'r' )
    out_file_list=out_fileee.readlines()
    columns=out_file_list[0].split()
    columns_chosen=[]
    columns_chosen.append(columns[2])
    columns_chosen.append(columns[3])
    for label in columns:
        if label.endswith('s') or label.endswith('a') or label.endswith('m'):
            columns_chosen.append(label)
    dictio_table={}
##    #generate dictionary (containing all) of dictionary of possible amino_acids containing empty lists
##    for amino_acid in amino_acid_list:
##        dictio_var={}
##        for label2 in columns_chosen:
##            dictio_var[label2]=[]
##        dictio_table[amino_acid]=dictio_var
    for amino_acid in amino_acid_list_3L :
##        if amino_acid!='GLY' and amino_acid!='ALA':
        dic_amino_acid={}
        for label in columns_chosen:
            var_list=[]
            for line in out_file_list:
                line_list=line.split()
                if line_list[1]==amino_acid:
                    for index_list_line in range(len(line_list) ) :
                        #print out_file_list[0][index_list_line]
                        if columns[index_list_line]==label:
                            if (amino_acid!='GLY' and amino_acid!='ALA' and label.endswith("s") ) or not label.endswith("s"):
                                var_list.append(line_list[index_list_line])
                            else:
                                var_list.append("1000")
            dic_amino_acid[label]=var_list
            dictio_table[amino_acid]=dic_amino_acid
    #for item in dictio_table:
   # print dictio_table
    out_output_file2=open(out_output_file,'w')
    out_output_file2.write('RT\t#models\tCI\tRN\t')
    for item in columns_chosen[2:]:
        out_output_file2.write(item+'(mean)\t'+item+'(dev)\t')
    out_output_file2.write('\n')
    for amino_acid in dictio_table:
        if not len(dictio_table[amino_acid]['CI'])==0:
            out_output_file2.write(amino_acid+'\t')
            numb_models=str(len(dictio_table[amino_acid]['CI']))
            out_output_file2.write(numb_models+'\t')
            for key in columns_chosen:
                if key=='CI' or key=='RN':
                    out_output_file2.write(dictio_table[amino_acid][key][0]+'\t')
                else:
                    if key=='NPm':
                        list=[int(i) for i in dictio_table[amino_acid][key]]
                    else:
                        list=[float(i) for i in dictio_table[amino_acid][key]]
                    mean=numpy.mean(list)
                    mean=mean_dev_return_good_value ( mean , 'mean' )
                    std=numpy.std(list)
                    std=mean_dev_return_good_value ( std , 'std' )
                    out_output_file2.write(mean+'\t')
                    out_output_file2.write(std+'\t')
            out_output_file2.write('\n')
    out_output_file2.close()

def mean_dev_return_good_value ( input , type_of_variable ) :
    if type_of_variable=='mean':
        if input<1 and input>-1:
            input='%.3f'%(input)
        else:
            input='%.1f'%(input)
    elif type_of_variable=='std':
        if input<1 and input>-1:
            input='%.4f'%(input)
        else:
            input='%.1f'%(input)
    else:
        print ('exception with value: '+input)
        print ('quitting')
        exit()
    return input

def retrieve_res_statistic_to_table ( out_input_file, min_res , max_res , output_file ):
    out_res_file_new_overall=open(output_file , 'w') # new table overall main chain protein
    out_res_file_new_overall.write('R#\tBAm(mean)\tBAm(dev)\tNPm(mean)\tNPm(dev)\tRm(mean)\tRm(dev)\tRGm(mean)\tRGm(dev)\tSRGm(mean)\tSRGm(dev)\tCCSm(mean)\tCCSm(dev)\tCCPm(mean)\tCCPm(dev)\tZCCm(mean)\tZCCm(dev)\tZOm(mean)\tZOm(dev)\tZDm(mean)\tZDm(dev)\tZD-m(mean)\tZD-m(dev)\tZD+m(mean)\tZD+m(dev)\n')
    for i in range (min_res , max_res+1):
        i_str=str(i)
        out_fileee=open( out_input_file+i_str, 'r' )
        out_file_list=out_fileee.readlines()
        out_res_file_new_overall.write(i_str+'\t')

        list_of_list=[]
        for column in range (4, 16 ):
            list=[]
            for line in range (1, len(out_file_list) ):
                line2=out_file_list[line].split('\t')
                if line2[column]!='n/a':
                    list.append( float(line2[column] ) )
            list_of_list.append(list)
        for listaa in list_of_list:
            mean=numpy.mean(listaa)
            if mean<1 and mean>-1:
                mean='%.3f'%(mean)
            else:
                mean='%.1f'%(mean)
            out_res_file_new_overall.write(mean+'\t')
            std=numpy.std(listaa)
            if std<1 and std>-1:
                std='%.4f'%(std)
            else:
                std='%.1f'%(std)
            out_res_file_new_overall.write(std+'\t')


        out_res_file_new_overall.write('\n')
    out_res_file_new_overall.close()
    
    
def extract_protein_chainID_res_number (pdb_input) :
    file=open(pdb_input,'r')
    file_list=file.readlines()
    file.close()
    DicChResNResType={}
    listChResNCA=[]
    chain_ID=[]
    res=[]
    number=0
    res_check = 'FirstOccurrence'
    chain_check=0
    res_string=''
    for line in file_list:
        if line.startswith('ATOM') or (line.startswith('HETATM') and line[17:20]=='MSE'):
            #line2=line.split()
            if line[12:15]==' CA' and [line[21],int(line[22:26])] not in listChResNCA: listChResNCA.append([line[21],int(line[22:26])]) #Constructing a list of lists [ [chain1,resnumb1] , [chain2,resnumb2] ]
            if res_check=='FirstOccurrence':
                res_check=int(line[22:26])
                chain_check=line[21]
                chain_ID.append(chain_check)
                residue=amino_acid_list[ amino_acid_list_3L.index(line[17:20]) ]
                DicChResNResType[chain_check]={res_check:residue}
                res_string=res_string+residue
            else:
                if not res_check==int(line[22:26]):
                    if chain_check==line[21]:
                        if int(line[22:26])==res_check+1:
                            res_check=res_check+1
                            residue=amino_acid_list[ amino_acid_list_3L.index(line[17:20]) ]
                            DicChResNResType[chain_check][res_check]=residue
                            res_string=res_string+residue
                        else:
                            res_check=int(line[22:26])
                            residue=amino_acid_list[ amino_acid_list_3L.index(line[17:20]) ]
                            DicChResNResType[chain_check][res_check] = residue
                            res_string=res_string+'/'+residue
                    else:
                        res.append(res_string)
                        res_string=''
                        chain_check=line[21]
                        chain_ID.append(chain_check)
                        residue=amino_acid_list[ amino_acid_list_3L.index(line[17:20]) ]
                        res_string=residue
                        res_check=int(line[22:26])
                        if chain_check not in DicChResNResType: DicChResNResType[chain_check] = {res_check: residue}
                        else:                                   DicChResNResType[chain_check][res_check]  = residue
    res.append(res_string)
    count=0
    for chain in res:
        count=count+len(chain.replace('/',''))
##    print 'verifying fn'
##    print 'chain_ID',chain_ID
##    print 'res',res
##    print 'count',count
    return chain_ID , res , count , DicChResNResType , listChResNCA

def extract_minimum_maximum_residue_number (pdb_input) : #script that extracts minimum and maximum residue number from a PDB
    file=open(pdb_input,'r')
    file_list=file.readlines()
    file.close()
    for line in file_list:
        if line.startswith('ATOM'):
            #taking minimum number
            try:
                minimum_res_numb
                if minimum_res_numb>int(line[23:26]):
                    minimum_res_numb=int(line[23:26])
            except:
                minimum_res_numb=int(line[23:26])
            #taking maximum number
            try:
                maximum_res_numb
                if maximum_res_numb<int(line[23:26]):
                    maximum_res_numb=int(line[23:26])
            except:
                maximum_res_numb=int(line[23:26])
    return minimum_res_numb , maximum_res_numb


def rewrite_seq_variables (seq_input, seq_output):  #script that given one sequence in a file containing just the sequence in one letter code, generates variation based on a dictionary. Such a dictionary was made by me looking at similarities of charge and chemical bonds of each residue
	amino_acid_group={"A":["G","S","T","V"],"C":["-"],"D":["T","N","G","A","L"],"E":["D","L","N","A","G"],"F":["Y","W","A","G"],"G":["A","V","T"],"H":["Y","W","I","A","G"],"I":["V","N","L","A","G"],"K":["R","Q","M","G","A"],"L":["I","V","N","A","G"],"M":["K","I","C","A","G"],"N":["Q","T","D","A","G"],"P":["F","V","A","G"],"Q":["N","E","H","G","A"],"R":["K","H","S","G","A"],"S":["T","N","Q","A","G"],"T":["S","N","Q","A","G"],"V":["T","S","A","G"],"W":["F","H","M","A","G"],"Y":["W","F","A","G","H"] }
	seq_inputtt=open(seq_input,"r")
	seq_input_string=seq_inputtt.read()
	seq_input_string=seq_input_string[:-1]
	list_write=["","","","","",""]
	for character in seq_input_string:
		list_write[0]=list_write[0]+character
		for i in range (5):
			try:
				list_write[i+1]=list_write[i+1]+amino_acid_group[character][i]
			except:
				list_write[i+1]=list_write[i+1]+"-"
	file_output=open(seq_output,"w")
	for seq in list_write:
		file_output.write(seq)
		file_output.write("\n")
	file_output.close()


def plot_line_error (input_file , output_file ): #output_file_without_.png
    input_fileee=open(input_file,"r")
    input_file_list=input_fileee.readlines()
    labels=input_file_list[0].split()
    gnuplot_instr=open("gnuplot_instr","w")
    gnuplot_instr.write("set terminal png font \"default\"\n")
    gnuplot_instr.write('set nokey\n')
    if input_file_list[0].startswith("RT"):
        gnuplot_instr.write("set style line 1 lc rgb \'grey30\' ps 0 lt 1 lw 2\n")
        gnuplot_instr.write("set style line 2 lc rgb 'grey70' lt 1 lw 2\n")
        gnuplot_instr.write("set label '*' at 3,0.8 center\n")
        gnuplot_instr.write("set label '*' at 4,0.8 center\n")
        #gnuplot_instr.write("set label '*' at 4,0.8 center\n")
        gnuplot_instr.write("set border 3\n")
        gnuplot_instr.write("set xtics nomirror scale 0\n")
        gnuplot_instr.write("set xtics rotate\n")
        for label in labels:
            if label.endswith("(mean)") and not label.startswith("ZDm"):
                gnuplot_instr.write('set ylabel \"'+label[:-6]+'\" \n')
                residue_number_table=output_file.split('_')[-1]
                gnuplot_instr.write('set xlabel \"Residue #'+residue_number_table+'\" \n')
                gnuplot_instr.write("set output \""+output_file+"_"+label[:-6]+".png\n")
                index_label=labels.index(label)
                index_label+=1
                gnuplot_instr.write("plot \""+input_file+"\" using 0:"+str(index_label)+":"+str(index_label+1)+" with yerrorbars ls 1, \"\" using 0:"+str(index_label)+":(0.7):xtic(1) with boxes ls 2\n")
    else:
        try:
            xrange_f=int( input_file_list[-1].split()[0] )
            xrange_f+=1
            xrange_i=int( input_file_list[1].split()[0] )
            xrange_i-=1
            gnuplot_instr.write("set xrange [")
            gnuplot_instr.write(str(xrange_i)+":"+str(xrange_f)+"]\n")
        except:
            blabla=0
        for label in labels:
            if label.endswith("(mean)") and not label.startswith("ZDm"):
                gnuplot_instr.write('set ylabel \"'+label[:-6]+'\" \n')
                gnuplot_instr.write('set xlabel \"Residues number\" \n')
                gnuplot_instr.write("set output \""+output_file+"_"+label[:-6]+".png\n")
                index_label=labels.index(label)
                index_label+=1
                gnuplot_instr.write("plot \""+input_file+"\" using 1:"+str(index_label)+" with lines lw 3 lt -1 , \"\" using 1:"+str(index_label)+":"+str(index_label+1)+" with yerrorbars \n")
    gnuplot_instr.close()
    os.system("gnuplot gnuplot_instr > /dev/null")


def plot_line (input_file , output_file , column_to_be_labeled , columns_to_be_evaluated): #script generated August 13rd, number given should not account for 0 as first column
    input_fileee=open(input_file,"r")
    input_file_list=input_fileee.readlines()
    labels=input_file_list[0].split()
    gnuplot_instr=open("gnuplot_instr","w")
    gnuplot_instr.write("set terminal png font \"default\"\n")
    gnuplot_instr.write('set nokey\n')
    try:
        xrange_f=int( input_file_list[-1].split()[column_to_be_labeled-1] )
        xrange_f+=1
        xrange_i=int( input_file_list[1].split()[column_to_be_labeled-1] )
        xrange_i-=1
        gnuplot_instr.write("set xrange [")
        gnuplot_instr.write(str(xrange_i)+":"+str(xrange_f)+"]\n")
    except:
        pass
    for label_number in columns_to_be_evaluated:
            gnuplot_instr.write('set ylabel \"'+labels[label_number-1]+'\" \n')
            gnuplot_instr.write('set xlabel \"Residues number\" \n')
            gnuplot_instr.write("set output \""+output_file+"_"+labels[label_number-1]+".png\n")
            gnuplot_instr.write("plot \""+input_file+"\" using "+str(column_to_be_labeled)+":"+str(label_number)+" with lines lw 3 lt -1 \n")#+", \"\" using 1:"+str(index_label)+":"+str(index_label+1)+" with yerrorbars \n")
    gnuplot_instr.close()
    os.system("gnuplot gnuplot_instr > /dev/null")


def generate_mark_true_res_SC_in_txt ( txt_input , true_res , txt_output ):
    txt_input_file=open(txt_input,'r')
    txt_input_list=txt_input_file.readlines()
    txt_input_file.close()
    txt_output_file=open(txt_output,'w')
    for line in txt_input_list:
        if line.startswith(true_res) and line[3]!="!":
            line_post_mortem=line[:3]+'!'+line[3:]
            txt_output_file.write(line_post_mortem)
        else:
            txt_output_file.write(line)
    txt_output_file.close()
    
def return_low_high_value_from_out_file ( statistic , txt_input ):
    txt_input_file=open(txt_input,'r')
    txt_input_list=txt_input_file.readlines()
    txt_input_file.close()
    labels=txt_input_list[0].split()
    index_statistic=labels.index(statistic)
    var_list=[]
    for line_number in range(1,len(txt_input_list)):
        var_list.append(float(txt_input_list[line_number].split()[index_statistic]))
    if statistic.startswith('CC') or statistic.startswith('ZCC') or statistic=='ZOs(mean)' or statistic=='ZD-s(mean)' :
        value=max(var_list)
    elif statistic.startswith('R') or statistic.startswith('BA') or statistic.startswith('ZD+'):
        value=min(var_list)
    value='%.3f'%value
    return value

def mtzfix (file_input,file_output,type_refinement_program):
    if not os.path.isfile( file_output ):
        if 'sigmaa' in type_refinement_program or 'sigmaa' in file_input:
            mtzfix_job = subprocess.Popen([ 'mtzfix','FLABEL F SIGF FC PHIC FC_ALL PHIC_ALL 2FOFCWT PH2FOFCWT FOFCWT PHFOFCWT FOM'.split(),'HKLIN',file_input,'HKLOUT',file_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
            #mtzfix_job = subprocess.Popen([ 'mtzfix','FLABEL F SIGF FC PHIC FC_ALL PHIC_ALL 2FOFCWT PH2FOFCWT FOFCWT PHFOFCWT FOM'.split(),'HKLIN',file_input,'HKLOUT',file_output], stdin=fft_job_instr, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            #mtzfix_job = subprocess.Popen([ 'mtzfix','HKLIN',file_input,'HKLOUT',file_output], stdin=fft_job_instr, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            mtzfix_job = subprocess.Popen([ 'mtzfix','HKLIN',file_input,'HKLOUT',file_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
        out, err = mtzfix_job.communicate()
        #print out[-110:]
        if not os.path.isfile( file_output ):
            print ("mtzfix did not generate mtzfixfile, previous file",file_input," being used.")
            os.system("cp "+file_input+" "+file_output)
        else:
            print ("mtzfix performed correctly in file:",file_input)

def RSZD_calculation (file_mtz, file_pdb , file_out , type_refinement_program): #script generated ?, upgraded for evaluating coot output mtz in May, 4th 2016 
    #type_refinement_program = refmac buster sigmaa (coot)
    file_mtz_fix=file_mtz[:-4]+'_mtzfix.mtz'
    file_map_fo=file_mtz[:-4]+'_fo.map'
    file_map_df=file_mtz[:-4]+'_df.map'
##    print 'mtz to be fixed:',file_mtz_fix
##    print 'map fo:',file_map_fo
##    print 'map fc:',file_map_df
    mtzfix (file_mtz,file_mtz_fix,type_refinement_program)
    #correct mtz through mtzfix line: mtzfix  [FLABEL <string>]  HKLIN in.mtz  HKLOUT out.mtz   [FLABEL <string>]=not necessary
    #generate maps through fft
    if 'BUSTER' in file_mtz or 'buster' in type_refinement_program or 'sigmaa' in type_refinement_program or 'sigmaa' in file_pdb or 'phenix.refine' in type_refinement_program:
        fft_script ( file_mtz_fix , file_map_fo , '2FOFCWT' , 'PH2FOFCWT' )
        fft_script ( file_mtz_fix , file_map_df , 'FOFCWT' , 'PHFOFCWT' )
    elif 'refmac' in file_mtz or 'refmac' in type_refinement_program:
        fft_script ( file_mtz_fix , file_map_fo , 'FWT' , 'PHWT' )
        fft_script ( file_mtz_fix , file_map_df , 'DELFWT' , 'PHDELWT'  )
    # elif 'phenix' in file_mtz or 'phenix' in type_refinement_program:
    #     print 'phenix.refine is under development, it may be advisable to use phenix.model_vs_data'
    #     exit()
    #     fft_script ( file_mtz_fix , file_map_fo , 'FWT' , 'PHWT' )
    #     fft_script ( file_mtz_fix , file_map_df , 'DELFWT' , 'PHDELWT'  )
    else:
        print ('Unable to find if mtz file '+file_mtz+' was coming from BUSTER or REFMAC or sigmaa (coot>get-eds-pdb-and-mtz), thus unable to generate edstats .out file')
    #extract resolution from mtz with mtzdmp
    low_res,high_res=extract_resolution_from_mtz(file_mtz)
    #generate out through edstats
    os.system('echo resl=' + low_res + ',resh=' + high_res + ' | edstats  MAPIN1 ' + file_map_fo + ' MAPIN2 ' + file_map_df + '  XYZIN ' + file_pdb + ' OUT ' + file_out + ' > ' + file_out[:-3] +'log')
    os.system('rm '+file_map_fo+' '+file_map_df+' '+file_mtz_fix)

def RSZD_calculation_delete (file_mtz, file_pdb , file_out , type_refinement_program): #script generated ?, upgraded for evaluating coot output mtz in May, 4th 2016
    #type_refinement_program = refmac buster sigmaa (coot)
    file_mtz_fix=file_mtz[:-4]+'_mtzfix.mtz'
    file_map_fo=file_pdb[:-4]+'_fo.map'
    file_map_df=file_pdb[:-4]+'_df.map'
    if 'BUSTER' in file_mtz or 'buster' in type_refinement_program or 'sigmaa' in type_refinement_program or 'sigmaa' in file_pdb or 'phenix.refine' in type_refinement_program :
        fft_script ( file_mtz_fix , file_map_fo , '2FOFCWT' , 'PH2FOFCWT' )
        fft_script ( file_mtz_fix , file_map_df , 'FOFCWT' , 'PHFOFCWT' )
    elif 'refmac' in file_mtz or 'refmac' in type_refinement_program:
        fft_script ( file_mtz_fix , file_map_fo , 'FWT' , 'PHWT' )
        fft_script ( file_mtz_fix , file_map_df , 'DELFWT' , 'PHDELWT'  )
    else:
        print ('Unable to find if mtz file '+file_mtz+' was coming from BUSTER or REFMAC or sigmaa (coot>get-eds-pdb-and-mtz), thus unable to generate edstats .out file')
    #extract resolution from mtz with mtzdmp
    low_res,high_res=extract_resolution_from_mtz(file_mtz)
    #generate out through edstats
    os.system('echo resl=' + low_res + ',resh=' + high_res + ' | edstats  MAPIN1 ' + file_map_fo + ' MAPIN2 ' + file_map_df + '  XYZIN ' + file_pdb + ' OUT ' + file_out + ' > ' + file_out[:-3] +'log')


def extract_resolution_from_mtz ( file_input_mtz ):
    os.system('mtzdmp ' + file_input_mtz + ' | grep " *  Resolution Range :" -A2 | grep A > '+file_input_mtz[:-4]+'_res.log')
    b=open(file_input_mtz[:-4]+'_res.log', 'r')
    c=b.read().split()
    b.close()
    low_res=c[3]
    high_res=c[5]
    os.system('rm '+file_input_mtz[:-4]+'_res.log')
    return low_res,high_res

def extract_resolution_from_mtz_delete ( file_input_mtz ):
    os.system('mtzdmp ' + file_input_mtz + ' | grep " *  Resolution Range :" -A2 | grep A > '+file_input_mtz[:-4]+'_res.log')
    b=open(file_input_mtz[:-4]+'_res.log', 'r')
    c=b.read().split()
    b.close()
    low_res=c[3]
    high_res=c[5]
    #os.system('rm '+file_input_mtz[:-4]+'_res.log')
    return low_res,high_res


def fft_script ( file_input_mtz , file_output_map , F , PHI ):
    #writting fft_job_description_file
    fft_job_description=open(file_output_map+'_fft_job.txt','w')
    fft_job_description.write('XYZLIM ASU\n')
    fft_job_description.write('GRID SAMPLE 4.5\n')
    fft_job_description.write('LABIN -\n')
    fft_job_description.write('     F1='+F+' PHI='+PHI)
    fft_job_description.write('\nEND')
    fft_job_description.close()
    fft_job_instr=open(file_output_map+'_fft_job.txt', 'r')
    fft_job = subprocess.Popen([ 'fft','HKLIN',file_input_mtz,'MAPOUT',file_output_map], stdin=fft_job_instr, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    out, err = fft_job.communicate()
    fft_job_instr.close()
    os.system('rm '+file_output_map+'_fft_job.txt')

def extract_list_mean_from_outfile ( input_outfile ): #script generated in 6 July 2015 to extract from out_file all the statistics´ label and mean
    input_outfileee=open ( input_outfile , "r")
    input_outfile_list=input_outfileee.readlines()
    input_outfileee.close()
    labels=[]
    final_list=[]
    first_out_line=input_outfile_list[0].split()
    for column in range ( 3,len( first_out_line ) ):
        labels.append ( first_out_line[column] )
        var_list=[]
        for line in range (1,len(input_outfile_list)):
            var_list.append ( float(input_outfile_list[line].split()[column]) )
        var_mean=numpy.mean(var_list)
        final_list.append (var_mean)
    return labels , final_list		


def extract_list_mean_std_from_outfile ( input_outfile ): #script generated in 16 July 2015 to extract from out_file all the statistics´ label and mean
    input_outfileee=open ( input_outfile , "r")
    input_outfile_list=input_outfileee.readlines()
    input_outfileee.close()
    labels=[]
    final_list_mean=[]
    final_list_std=[]
    first_out_line=input_outfile_list[0].split()
    for column in range ( 3,len( first_out_line ) ):
        labels.append ( first_out_line[column] )
        var_list=[]
        for line in range (1,len(input_outfile_list)):
            var_list.append ( float(input_outfile_list[line].split()[column]) )
        var_mean=numpy.mean(var_list)
        var_std=numpy.std(var_list)
        if var_mean<1 and var_mean>-1:
            var_mean='%.3f'%(var_mean)
        else:
            var_mean='%.1f'%(var_mean)
        if var_std<1 and var_std>-1:
            var_std='%.4f'%(var_std)
        else:
            var_std='%.1f'%(var_std)
        final_list_mean.append (var_mean)
        final_list_std.append (var_std)
    return labels , final_list_mean , final_list_std


def extract_chain_statistics_from_edstats_log_file ( edstats_log_file ) :
    file=open( edstats_log_file , 'r' )
    log=file.readlines()
    file.close()
    dic_all={}
    dic_overall=defaultdict(float)
    for ind,line in enumerate(log):
        if line.startswith ('DF map Q-Q difference plot statistics for chain'):
            chain=line[-4]
            dic_var=defaultdict(float)
            value0=log[ind+3].split()
            dic_overall['Biso']=float(value0[-1])
            value1=log[ind+7].split()
            dic_var['Q-Qplot+']=float(value1[-1])
            dic_var['Q-Qplot-']=float(value1[-2])
            value2=log[ind+8].split()
            dic_var['Q-QplotZD+']=float(value2[-1])
            dic_var['Q-QplotZD-']=float(value2[-2])
            dic_all[chain]=dic_var
        elif line.startswith('DF map Q-Q difference plot statistics for all density values'):
            value1=log[ind+5].split()
            dic_overall['Q-Qplot+']=float(value1[-1])
            dic_overall['Q-Qplot-']=float(value1[-2])
            value2=log[ind+6].split()
            dic_overall['Q-QplotZD+']=float(value2[-1])
            dic_overall['Q-QplotZD-']=float(value2[-2])
        elif line.startswith('RMS Z-scores for protein residues'):
            value1=log[ind+3].split()
            dic_overall[value1[0]+'s']=float(value1[-2])
            dic_overall[value1[0]+'all']=float(value1[-1])
            value2=log[ind+4].split()
            dic_overall[value2[0]+'s']=float(value2[-2])
            dic_overall[value2[0]+'all']=float(value2[-1])
        elif line.startswith('RSZD- score, count & percentage of protein residues with worse score:'):
            value1=log[ind+5].split()
            try:
                dic_overall['n_'+value1[0][:2]+'s']=float(value1[3])
                dic_overall['n_'+value1[0][:2]+'all']=float(value1[-2])
            except:
                dic_overall['n_-2s']=0
                dic_overall['n_-2all']=0
            value2=log[ind+6].split()
            try:
                dic_overall['n_'+value2[0][:2]+'s']=float(value2[3])
                dic_overall['n_'+value2[0][:2]+'all']=float(value2[-2])
            except:
                dic_overall['n_-3s']=0
                dic_overall['n_-3all']=0
        elif line.startswith('Overall mean Biso:'):
            dic_overall['Biso']=float(line.split()[-1])
    return dic_overall,dic_all


def add_Bfactors_occup_to_pdb ( pdb_input_file , pdb_output_file , Bf_number_in_string_6characters , occ_number_in_string_3characters=False ): #example Bf_number_in_string_5characters=" 20.00" or "  4.00" and occ_number_in_string_3characters_'1.00' or '0.00'
    if len(Bf_number_in_string_6characters)<6: Bf_number_in_string_6characters=' '*(6-len(Bf_number_in_string_6characters))+Bf_number_in_string_6characters
    with open(pdb_input_file) as f: f2=f.readlines()
    with open(pdb_output_file,'w') as f:
        for line in f2:
            if line.startswith('ATOM'):
                if occ_number_in_string_3characters == False: occ = line[56:60]
                else:                                         occ = occ_number_in_string_3characters
                f.write( line[:56] + occ + Bf_number_in_string_6characters + line[66:] )
            else:f.write(line)
            
            
            
def read_table_file_output_table_file_with_normalized_values ( input_file , list_strings_labels_to_be_normalized , output_file ): #changed the List_Column_to_be_normalized to list of strings
    input_fileee=open( input_file , 'r')
    input_file=input_fileee.readlines()
    input_fileee.close()
    
    number_of_columns=len(input_file[0].split())
    all_columns=input_file[0].split()
    
    mean_0=[]
    stdev_0=[]

##    for i in range (number_of_columns):
##        if i+1 in List_Column_to_be_normalized:
##            normalized_labels.append(all_columns[i])
    
    for i in range( number_of_columns ):
        if all_columns[i] in list_strings_labels_to_be_normalized:
            var_values=[]
            for line in input_file[1:]:
                if not line.startswith('initial'):
                    line=line.split()
                    if line[i]!='n/a': var_values.append(float(line[i]))
            mean_0.append(numpy.mean(var_values))
            stdev_0.append(numpy.std(var_values))
    # print list_strings_labels_to_be_normalized
    # for i in range(len(list_strings_labels_to_be_normalized)):
    #     print "Label normalized:"
    #     print list_strings_labels_to_be_normalized[i]
    #     print "Mean:"
    #     print mean_0[i]
    #     print "STDEV:"
    #     print stdev_0[i]
    #     print "\n\n"

    output_file2=open ( output_file , "w" )
    output_file2.write(input_file[0])
    for line in input_file[1:]:
        if not line.startswith('initial'):
            line=line.split()
            for i in range (number_of_columns):
                if all_columns[i] in list_strings_labels_to_be_normalized and line[i]!='n/a':
                    index_normalized_label=list_strings_labels_to_be_normalized.index(all_columns[i])
                    normalized_value= (  float(line[i]) - mean_0[index_normalized_label])/stdev_0[index_normalized_label]
                    if normalized_value<1 and normalized_value>-1:
                        normalized_value='%.3f'%(normalized_value)
                    else:
                        normalized_value='%.1f'%(normalized_value)
                    output_file2.write(str(normalized_value+"\t"))

                else:
                    output_file2.write(line[i]+"\t")
            output_file2.write("\n")
    output_file2.close()
    
    
def given_two_sequences_return_identity_similarity(seq1_string,seq2_string):
# 1 letter code residues
# ALA/GLY are assumed to be the same
    dictio_similarity={"s":['A',"G","T","S"],"c":['C'],"n":['D','E'],"a":['F',"Y","W"],"p":['H','K','R'],"l":['I','L',"V"],"m":['M'],"o":['N',"Q"],"r":["P"],'x':['X']}
    #groups discription:       [small],       ['C'],  [negative],  [aromatic],        [positive],      [aliphatic],     ['M'],     [polar]]        [rigid]
    #based on: http://imed.med.ucm.es/Tools/sias.html
    identity_check=0
    similarity_check=0
    seq1_sim=""
    seq2_sim=""
##    print "seq1:",seq1_string
##    print "seq2:",seq2_string
    if len (seq1_string)!=len(seq2_string):
        print (seq1_string,seq2_string,"with different lengths (",len (seq1_string),len(seq2_string),")")
        exit()
    for i in seq1_string:
        for key in dictio_similarity:
            if i.upper() in dictio_similarity[key]:
                seq1_sim+=key
    for i in seq2_string:
        for key in dictio_similarity:
            if i.upper() in dictio_similarity[key]:
                seq2_sim+=key
##    print "seq1_sim:",seq1_sim
##    print "seq2_sim:",seq2_sim
    seq1_string=seq1_string.replace("A","G")
    seq1_string=seq1_string.replace("a","g")
    seq2_string=seq2_string.replace("A","G")
    seq2_string=seq2_string.replace("a","g")
    for i in range(len(seq1_string)):
        if seq1_string[i]==seq2_string[i]:
            identity_check+=1
        if seq1_sim[i]==seq2_sim[i]:
            similarity_check+=1
##    print identity_check
##    print similarity_check
    identity=identity_check*100/len(seq1_string)
    similarity=similarity_check*100/len(seq1_string)
    return identity,similarity

def return_indices_character_in_string (string, character): #return list of all indeces given character
    return [i for i, char in enumerate(string) if char == character]

def remove_given_list_characters_of_both_strings_together (seq1_string,seq2_string,list_characters): #ex: it will remove by index from both strings even with one has not that character in given position 
##    print '\n\n\n\n'
##    print list_characters
##    print seq1_string
##    print seq2_string
    if len (seq1_string)!=len(seq2_string):
        print (seq1_string,seq2_string,"with different lengths (",len (seq1_string),len(seq2_string),")")
        exit()
    list_unwanted_indices=[]
    for char in list_characters:
##        print char
        list_unwanted_indices1=return_indices_character_in_string (seq1_string, char)
        list_unwanted_indices2=return_indices_character_in_string (seq2_string, char)
        list_unwanted_indices3=list_unwanted_indices1 + list(set(list_unwanted_indices2) - set(list_unwanted_indices1))
        list_unwanted_indices=list_unwanted_indices + list(set(list_unwanted_indices3) - set(list_unwanted_indices))
    new_seq1=''
    new_seq2=''
    for ind in range(len(seq1_string)):
        if ind not in list_unwanted_indices:
            new_seq1+=seq1_string[ind]
            new_seq2+=seq2_string[ind]
##    print new_seq1
##    print new_seq2
    return new_seq1,new_seq2

def given_two_sequences_return_excluding_list_characters_identity_similarity(seq1_string,seq2_string,list_characters): #excludes X or x before assigning: modified 10Jun2016 from given_two_sequences_return_identity_similarity function to generate possible identity and possible similarity
    # 1 letter code residues
    # ALA/GLY are assumed to be the same
    seq1_string,seq2_string=remove_given_list_characters_of_both_strings_together (seq1_string,seq2_string,list_characters)
    identity,similarity=given_two_sequences_return_identity_similarity(seq1_string,seq2_string)
    return identity,similarity


def given_two_sequences_return_evaluated_identity_similarity(seq1_string,seq2_string): #modified 10Jun2016 from given_two_sequences_return_identity_similarity function to generate possible identity and possible similarity
# 1 letter code residues
# ALA/GLY are assumed to be the same
    dictio_similarity={"s":['A',"G","T","S"],"c":['C'],"n":['D','E'],"a":['F',"Y","W"],"p":['H','K','R'],"l":['I','L',"V"],"m":['M'],"o":['N',"Q"],"r":["P"],'x':["X"]}
    #groups discription:       [small],       ['C'],  [negative],  [aromatic],        [positive],      [aliphatic],     ['M'],     [polar]]        [rigid]
    #based on: http://imed.med.ucm.es/Tools/sias.html
    identity_check=0
    similarity_check=0
    seq1_sim=""
    seq2_sim=""
##    print "seq1:",seq1_string
##    print "seq2:",seq2_string
    total_string=len(seq1_string)
    if len (seq1_string)!=len(seq2_string):
        print (seq1_string,seq2_string,"with different lengths (",len (seq1_string),len(seq2_string),")")
        exit()
    for i in seq1_string:
        for key in dictio_similarity:
            if i in dictio_similarity[key]:
                seq1_sim+=key
    for i in seq2_string:
        for key in dictio_similarity:
            if i in dictio_similarity[key]:
                seq2_sim+=key
##    print "seq1_sim:",seq1_sim
##    print "seq2_sim:",seq2_sim
    seq1_string=seq1_string.replace("A","G")
    seq2_string=seq2_string.replace("A","G")
    for i in range(len(seq1_string)):
        if seq1_string[i]==seq2_string[i] and seq1_string[i]!='X':
            identity_check+=1
        if seq1_sim[i]==seq2_sim[i]  and seq1_string[i]!='x':
            similarity_check+=1
##    print identity_check
##    print similarity_check
    identity=identity_check*100/total_string
    similarity=similarity_check*100/total_string
    return identity,similarity


def given_two_sequences_return_its_identity_in_0_1 (seq1_string,seq2_string):
    seq1_string=seq1_string.replace("A","G")
    seq2_string=seq2_string.replace("A","G")
    seq=""
    for i in range(len(seq1_string)):
        if seq1_string[i]==seq2_string[i]:
            seq+="1"
        else:
            seq+="0"
    return seq


def refine_refmac_multiprocesses ( pdb_input_file , mtz_input_file , pdb_output_file , refmac_tmp , refmac_path ):
    refmacTMP=open(refmac_tmp,'r')
    xyzin=pdb_input_file
    xyzout=pdb_output_file
    hklin=mtz_input_file
    hklout=xyzout[:-4]+'.mtz'
    cifout=xyzout[:-4]+'.cif'
    refmac_output_log = open(xyzout[:-4]+'.log', 'w')
    p = subprocess.Popen([refmac_path,'XYZIN',xyzin,'XYZOUT',xyzout,'HKLIN',hklin,'HKLOUT',hklout,'LIBOUT',cifout], stdin=refmacTMP, stdout=refmac_output_log, stderr=subprocess.PIPE , text=True)
    out, err = p.communicate()
    refmacTMP.close()
    refmac_output_log.close()
    if not os.path.isfile( xyzout ):
        print ('ERROR in refinement.\nFile',xyzout,'not generated.')
        exit()

def sh_refine_refmac_multiprocesses ( pdb_input_file , mtz_input_file , pdb_output_file , refmac_tmp_file , refmac_path , sh_file , distribute_computing ):
    file_sh=open(sh_file,'w')
##    if distribute_computing=='multiprocessing':
##        file_sh.write('#!/bin/sh'+ '\n')
    if distribute_computing=='local_grid':
        file_sh.write('#!/bin/tcsh'+ '\n')
        #file_sh.write('source /xtal/xtalsetup'+ '\n')
        file_sh.write('hostname > '+sh_file.split('/')[-1]+ '_host \n')
        #file_sh.write('source /xtal/ccp4/ccp4/include/ccp4.setup\n') # obsolute after May 10, 2019
        file_sh.write('source /xtal/ccp4/ccp4-7.0/include/ccp4.setup-csh.in\n')  # updated on May  10, 2019
##        file_sh.write('source /xtal/ccp4/ccp4/include/ccp4.setup\n')
##        file_sh.write('source /xtal/buster/buster/setup.csh\n')
##        #file_sh.write('setenv LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu\n')
##        file_sh.write('setenv PATH /usr/local/bin\n')
##        file_sh.write('setenv PATH /usr/bin:$PATH\nsetenv PATH /bin:$PATH\n')
##        #file_sh.write('echo $LD_LIBRARY_PATH > '+sh_file.split('/')[-1]+'_lib'+ '\n')
##        #file_sh.write('echo $PATH > '+sh_file.split('/')[-1]+'_path'+ '\n')
##        #file_sh.write('sleep 60'+ '\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /lib/x86_64-linux-gnu\n')
##        file_sh.write('setenv PATH /usr/local/bin\n')
##        file_sh.write('setenv PATH /usr/bin:$PATH\n')
##        file_sh.write('setenv PATH /bin:$PATH\n')
        
##        file_sh.write('setenv LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MCR/lib/:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/runtime/glnxa64:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/bin/glnxa64:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/sys/os/glnxa64:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/sys/java/jre/glnxa64/jre/lib/amd64/server:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/sys/java/jre/glnxa64/jre/lib/amd64/:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/atsas/atsas/lib/atsas:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/Xmipp/Xmipp-2.4-x64/lib:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/imod/imod-4.5.3-RH6-64/IMOD/lib:$LD_LIBRARY_PATH\n')
        #file_sh.write('setenv LD_LIBRARY_PATH /xtal/imp/imp/lib64:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/XtalView/lib/ibmpcLinux2:$LD_LIBRARY_PATH\n')
##        file_sh.write('source /xtal/xtalsetup\n')
##        file_sh.write('source /xtal/ccp4/ccp4/include/ccp4.setup\n')
##        file_sh.write('source /xtal/buster/buster/setup.csh\n')
    xyzin=pdb_input_file
    xyzout=pdb_output_file
    hklin=mtz_input_file
    hklout=xyzout[:-4]+'.mtz'
    cifout=xyzout[:-4]+'.cif'
    refmac_output_log = xyzout[:-4]+'.log'
    if distribute_computing=='local_grid':
        xyzin =xyzin.split('/')[-1]
        xyzout=xyzout.split('/')[-1]
        hklin =hklin.split('/')[-1]
        hklout=hklout.split('/')[-1]
        cifout=cifout.split('/')[-1]
        refmac_output_log = refmac_output_log.split('/')[-1]
        refmac_tmp_file=refmac_tmp_file.split('/')[-1]
    #p = subprocess.Popen([refmac_path,'XYZIN',xyzin,'XYZOUT',xyzout,'HKLIN',hklin,'HKLOUT',hklout,'LIBOUT',cifout], stdin=refmacTMP, stdout=refmac_output_log, stderr=subprocess.PIPE)
    #os.system( '/xtal/ccp4/ccp4-6.4.0/bin/refmac5 XYZIN "0_Scwrl4/seq' + str(b) + '.pdb" XYZOUT "1_refmac/seq' + str(b) + '_refmac1.pdb" HKLIN "' + mtz_file + '" HKLOUT "1_refmac/seq' + str(b) + '_refmac1.mtz" LIBOUT "1_refmac/seq' + str(b) + '_refmac1.cif" < 1_refmac/refmacTMP.tmp > 1_refmac/seq' + str(b) + '_refmac1.log ')
    file_sh.write(refmac_path+' XYZIN '+xyzin+' XYZOUT '+xyzout+' HKLIN '+hklin+' HKLOUT '+hklout+' LIBOUT '+cifout+" < "+refmac_tmp_file+" > "+refmac_output_log+' \n\n')
    file_sh.close()

def refine_refmac_singleprocess ( pdb_file , mtz_file , refmac_tmp ):
    refmac_path=which_program ("refmac5")    
    xyzin=pdb_file
    xyzout=xyzin[:-4]+'_refmac1.pdb'
    hklin=mtz_file
    hklout=xyzout[:-4]+'.mtz'
    cifout=xyzout[:-4]+'.cif'
    logout=xyzout[:-4]+'.log'
    os.system(refmac_path+' XYZIN '+xyzin+' XYZOUT '+xyzout+' HKLIN '+hklin+' HKLOUT '+hklout+' LIBOUT '+cifout+" < "+refmac_tmp+" > "+logout)



def which_program (name_program):
    os.system("which "+name_program+" > path_program")
    program_path_file=open("path_program","r")
    program_path=program_path_file.read()
    program_path=program_path[:-1]
    program_path_file.close()
    return program_path

def number_of_processor():
    os.system('grep "model name" /proc/cpuinfo | wc -l > nproc')  # if one wants number of threads
    nproc_file=open("nproc","r")
    nproc=int(nproc_file.read())
    nproc_file.close()
    print (nproc,"cores found")
    return nproc

def extract_wWPE_from_lst (lst_file): #Script generated in July 17th 2015 for script EDSTATS_ARCIMBOLDO_evaluation.py
    os.system('grep " model in .ent:" -B3 '+lst_file+' > wMPE.log')
    file=open("wMPE.log","r")
    lst_list=file.readlines()
    file.close()
##    os.system("rm wMPE.log")
    list_wMPE=[]
##    I am taking the wMPE generated by the last cycle of density modification. Each number in this range corresponds to each cycle of auto_tracing, if list_wMPE[-1] last cycle and list_wMPE[0] no cycle of auto_tracing
    for i in range(0,len(lst_list),5):
        line=lst_list[i].split()
        list_wMPE.append(line[-3])
    return list_wMPE

def extract_CC_from_lst (lst_file): #Script generated in July 17th 2015 for script EDSTATS_ARCIMBOLDO_evaluation.py
    os.system('grep "CC" '+lst_file+' > CC.log')
    file=open("CC.log","r")
    lst_list=file.readlines()
    file.close()
##    os.system("rm CC.log")
    list_CC=[]
##    I am taking the wMPE generated by the last cycle of density modification. Each number in this range corresponds to each cycle of auto_tracing, if list_wMPE[-1] last cycle and list_wMPE[0] no cycle of auto_tracing
    for line in lst_list:
        line=line.split()
        if len(list_CC)==0:
            if line[0]=="Overall":
                list_CC.append(line[-1])
        else:
            if line[0]=="Best":
                list_CC.append(line[6][:-1])
    return list_CC



def massimoLIB_readCCValFromSUM (sumPath):
    if not os.path.exists(sumPath):
        print (sumPath,"does not exist")
        exit()
    f = open(sumPath,"r")
    sum_list = f.readlines()
    dictio_all = {}
##    while line != None and line != "":
    for line_index in range(len(sum_list)-1):
        #process line
        if sum_list[line_index].startswith("=======") and sum_list[line_index+1].startswith("MODEL"):
##            print "start reading CC_Val n. "+str(numClus)
            riga1 = sum_list[line_index]
            riga2 = sum_list[line_index+1].split()
            riga3 = sum_list[line_index+2].split()
            riga4 = sum_list[line_index+3].split()

            model = int(riga2[1])
            corresp = riga2[3]
            cluster = riga2[5]
            if cluster == "None":
                cluster = "None"
            else:
                cluster = cluster

            nAtoms = int(riga3[1])
            nER = int(riga3[3])
            initCC = float(riga3[5])
            finalCC = float(riga3[7])
            wMPEa = float(riga3[9])
            wMPEb = float(riga3[11])
            wMPEc = float(riga3[13])
            wMPEd = float(riga3[15])

            shx = float(riga4[1])
            shy = float(riga4[2])
            shz = float(riga4[3])
            contrast = float(riga4[5])
            connect = float(riga4[7])
            dizio = {"model":model,"corresp":corresp,"natoms":nAtoms, "ner":nER,"initcc":initCC,"finalcc":finalCC,"cluster":cluster,"wMPE_init":[wMPEa,wMPEb],"wMPE_end":[wMPEc,wMPEd], "shift_origin":[shx,shy,shz], "contrast":contrast, "connect":connect}
##            mfom = float(riga4[9])
##            sfom = float(riga4[11])
##            dizio = {"model":model,"corresp":corresp,"natoms":nAtoms, "ner":nER,"initcc":initCC,"finalcc":finalCC,"cluster":cluster,"wMPE_init":[wMPEa,wMPEb],"wMPE_end":[wMPEc,wMPEd], "shift_origin":[shx,shy,shz], "contrast":contrast, "connect":connect,"mfom":mfom,"sfom":sfom}
            dictio_all[dizio["corresp"]]=dizio
    return (dictio_all)



def run_Scwrl4_multiprocess (pdb_input, seq_input, pdb_output, log_output , Scwrl4_path ) :
    #Scwrl4_path=which_program ("Scwrl4")
    options_Scwrl4='-h -t -#'
    p = subprocess.Popen ( [Scwrl4_path,'-i',pdb_input,'-o',pdb_output,'-s',seq_input]+options_Scwrl4.split(),stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
##    print "pdb_input is:",pdb_input
##    print "pdb_output is:",pdb_output
##    print "seq_input is:",seq_input
##    print "log_output is:",log_output
    out, err = p.communicate()
    file=open(log_output,'w')
    file.write(out)
    file.close()
    if not os.path.isfile( pdb_output ):
        print ('ERROR in Scwrl4 generation.\nFile',pdb_output,'not generated.')
        exit()
    add_Bfactors_occup_to_pdb(pdb_input_file=pdb_output, pdb_output_file=pdb_output,Bf_number_in_string_6characters=' 20.00', occ_number_in_string_3characters='1.00')


def run_Scwrl4_singleprocess (pdb_input, seq_input, pdb_output, log_output) :
    Scwrl4_path=which_program ("Scwrl4")
    options_Scwrl4=' -h -t -# '
    os.system(Scwrl4_path+' -i '+pdb_input+' -o '+pdb_output+' -s '+seq_input+options_Scwrl4+" > "+log_output)
    print ("pdb_input is:",pdb_input)
    print ("pdb_output is:",pdb_output)
    print ("seq_input is:",seq_input)
    print ("log_output is:",log_output)

def coot_ROTAMER_automatic_refinement ( input_PDB_file , input_map_file , output_file , number_process , coot_path=False ) :
    script_file=open('coot_script_'+str(number_process),'w')
    script_file.write('(turn-off-backup 0)\n')
    script_file.write('(fit-protein 0)\n')
    script_file.write('(stepped-refine-protein-for-rama 0)\n')
    script_file.write('(write-pdb-file 0 "'+output_file+'")\n' )
    script_file.write('(coot-real-exit 0)' )
    script_file.close()
    if coot_path==False:
        coot_path='coot'
    os.system( coot_path +' --pdb ' + input_PDB_file + ' --auto ' + input_map_file + ' --no-guano -s coot_script_'+str(number_process)+' --no-graphics > log_coot_ROTAMER'+str(number_process))
    os.system( 'rm coot_script_'+str(number_process)+" log_coot_ROTAMER"+str(number_process))
    if not os.path.isfile( output_file ):
        print ('ERROR in COOT MODELING.\nFile',output_file,'not generated.')
        exit()

# line: coot --pdb input_PDB_file --auto input_map_file --no-guano -s coot_script --no-graphics > log_coot_ROTAMER    
# --no-guano                        do not save 0... scripts of coot to backup and reinicialize it from where it stoped
# (turn-off-backup 0)			    do not save 'coot-backup'
#(fit-protein 0)				    fit protein of macromolecule '0'
#(write-pdb-file 0 "output_file")	save macromolecule '0' in file-name output-file ("" necessary)
#(coot-real-exit 0)			        quit

def coot_py_script_rotamer_sphere_refine ( pdb , dic , outputCootPy , output_pdb  , radius_sph_ref=5 , removeWAT=False , radius_sph_removeWat=2.5): #Script created Dec5,2016. Dic should be: dic[chain][resnumb]=False if no mutation is desired or dic[chain][resnumb]='aminoacid_1L' if mutation is desired
    if not outputCootPy.endswith('.py'):
        print ('wrong variable outputCootPy',outputCootPy,'given in function RJB_lib.coot_py_script_sphere_refine')
        exit()
    if removeWAT:
        line_remWAT = "other_residues = residues_near_residue(0, centred_residue, " + str(radius_sph_removeWat) + ")\n"
        line_remWAT+="for resi in other_residues:\n    chain=resi[0]\n    resnumb=resi[1]\n    resi_name=resname_from_serial_number (0,chain,resnumb)\n    if resi_name=='WAT' or resi_name=='HOH': delete_atom (0 , chain, resnumb , '' , resi_name , '' ) \n"
    else: line_remWAT=''
    ou=open(outputCootPy,'w')
    ou.write( 'turn_off_backup(0) \n')
    dic_pdb=return_dic_chain_resnumb_restype (pdb)
    for ch,item in dic.items():
        for resn,mut in item.items():
            if mut!=False:
                if len(mut)==1:
                    mut=amino_acid_list_3L[amino_acid_list.index(mut)]
                if amino_acid_list[amino_acid_list_3L.index(mut)]==dic_pdb[ch][resn]:
                    if mut!='ALA':
                        ou.write( 'mutate(0,"'+ch+'",'+str(resn)+',"","'+'ALA'+'")\n')
                    else:
                        ou.write( 'mutate(0,"'+ch+'",'+str(resn)+',"","'+'GLY'+'")\n')
                ou.write( 'mutate(0,"'+ch+'",'+str(resn)+',"","'+mut+'")\n')
                ou.write( 'auto_fit_best_rotamer ( '+str(resn)+', "", "", "'+ch+'", 0 , 1 , 1 , 0.10)\n')
                if float(radius_sph_ref)<50.0:
                    ou.write( 'centred_residue=["'+ch+'",'+str(resn)+',""]\n')
                    ou.write( line_remWAT )
                    ou.write( 'other_residues = residues_near_residue(0, centred_residue, '+str(radius_sph_ref)+')\n')
                    ou.write( 'all_residuess = [centred_residue]\n')
                    ou.write( 'if (type(other_residues) is ListType):\n')
                    ou.write( '    all_residuess += other_residues\n')
                    ou.write( 'refine_residues(0, all_residuess)\n')
                    ou.write( 'accept_regularizement()\n')
    if float(radius_sph_ref)>=50.0:
        ou.write( 'centred_residue=["'+ch+'",'+str(resn)+',""]\n')
        ou.write( 'other_residues = residues_near_residue(0, centred_residue, '+str(radius_sph_ref)+')\n')
        ou.write( 'all_residuess = [centred_residue]\n')
        ou.write( 'if (type(other_residues) is ListType):\n')
        ou.write( '    all_residuess += other_residues\n')
        ou.write( 'refine_residues(0, all_residuess)\n')
        ou.write( 'accept_regularizement()\n')
    ###auto-fit-best-rotamer ( resnumb, "", "", chain-id, 0 , 1 , 1 , 0.10)
    ou.write( 'write_pdb_file(0,"'+output_pdb+'")\n')
    ou.write( 'coot_real_exit (0)\n')
    ou.close()

##turn_off_backup(0)
##mutate(0,"A",1,"","TRP")
##auto_fit_best_rotamer ( 1, "", "", "A", 0 , 1 , 1 , 0.10)
##centred_residue=['A',1,'']
##other_residues = residues_near_residue(0, centred_residue, 5)
##all_residuess = [centred_residue]
##if (type(other_residues) is ListType):
##    all_residuess += other_residues
##refine_residues(0, all_residuess)
##accept_regularizement()
##write_pdb_file(0,"test.pdb")
##coot_real_exit (0)

##adding command to remove water in surrounding of chosing residues

def coot_run_rotamer_sphere_refinement ( input_PDB_file , input_mtz_file , output_pdb , outputCootPy , dic , radius_sph_ref=5 , coot_path=False ) : #Script created Dec5,2016
    coot_py_script_rotamer_sphere_refine ( input_PDB_file , dic , outputCootPy , output_pdb , radius_sph_ref=5 )
    if coot_path==False:
        coot_path='coot'
    os.system( coot_path +' --pdb ' + input_PDB_file + ' --auto ' + input_mtz_file + ' --no-guano -s '+outputCootPy+' --no-graphics > '+outputCootPy+'_log')
    os.system( 'rm '+outputCootPy+' '+outputCootPy+'_log')
    if not os.path.isfile( output_file ):
        print ('ERROR in COOT MODELING.\nFile',output_file,'not generated.')
        exit()

def coot_run_rotamer_sphere_refinement_multi ( input_PDB_file , input_mtz_file , output_pdb , outputCootPy , dic , radius_sph_ref=5  , removeWAT=False , coot_path=False ) : #Script created Dec5,2016
    coot_py_script_rotamer_sphere_refine ( input_PDB_file , dic , outputCootPy , output_pdb , radius_sph_ref , removeWAT )
    if coot_path==False:
        coot_path='coot'
    #os.system( coot_path +' --pdb ' + input_PDB_file + ' --auto ' + input_mtz_file + ' --no-guano -s '+outputCootPy+' --no-graphics > '+outputCootPy+'_log')
    #print coot_path,'--pdb',input_PDB_file,'--auto',input_mtz_file,'--no-guano','-s',outputCootPy,'--no-graphics'
    p = subprocess.Popen([coot_path,'--pdb',input_PDB_file,'--auto',input_mtz_file,'--no-guano','-s',outputCootPy,'--no-graphics'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    out, err = p.communicate()
    #file=open(output_pdb[:-4]+'_coot.log','w')
    #file.write(out)
    #file.close()
    #os.system( 'rm '+outputCootPy+' '+outputCootPy+'_log')
    if not os.path.isfile( output_pdb ):
        print ('ERROR in COOT MODELING.\nFile',output_pdb,'not generated.')
        print (out)
        print (err)
        exit()

def coot_mutate_sph_ref_correcting_files (pdb,mtz_phases,outpdb,dic,dic_pdb=False, radius_sph_ref=5, removeWAT=False , coot_path=False , printtt=False): #Script created Dec5,2016 dic_pdb should come from fn [{ return_dic_chain_resnumb_restype (pdb) }]
    if dic_pdb==False or dic_pdb=='':
        dic_pdb=return_dic_chain_resnumb_restype (pdb)
    dic_new_res={}
    dic_impres=return_impartial_res (pdb)
    for ch,dresnrest in dic.items():
        for i,lmut in dresnrest.items():
            if printtt: print ('SEQ',i,'in chain',ch,'in PDB was originnaly',dic_pdb[ch][i])
            try:
                if dic[ch][i]!=dic_pdb[ch][i] or dic_impres[ch][i]==dic[ch][i]:
                    try:
                        dic_new_res[ch][i]=dic[ch][i]
                    except:
                        dic_new_res[ch]={i:dic[ch][i]}
            except:
                pass
    #print 'dic is',dic
    #print 'dic_new_res',dic_new_res
    if len(dic_new_res)>0:
        for ch,dii in dic_new_res.items():
            for i,rest in dii.items():
                if printtt: print ('New assigned residue in chain',ch,'by coot is:',rest)
        if not os.path.isfile(outpdb):
            coot_run_rotamer_sphere_refinement_multi( input_PDB_file=pdb , input_mtz_file=mtz_phases , output_pdb=outpdb , outputCootPy=outpdb[:-3]+'py' , dic=dic_new_res , radius_sph_ref=radius_sph_ref , removeWAT=removeWAT , coot_path=coot_path )
        if not os.path.isfile(outpdb):
            print ('ERROR in COOT MODELING.\nFile',outpdb,'not generated.')
            exit()
        if printtt: print ('\n\n')
    else:
        if not os.path.isfile(outpdb):
            remove_Bfactor_occ_res_pdb (pdb_input=pdb,pdb_output=outpdb,dic_ch_resn=dic)


##def coot_mutate_sph_ref_polder_multi (pdb,mtz_init,mtz_phases,outpdb,outlog,dic, radius_sph_ref=5, coot_path=False , polder_path=False , output_map=False): #Script created Dec5,2016
##    dic_pdb=return_dic_chain_resnumb_restype (pdb)
##    dic_new_res={}
##    dic_impres=return_impartial_res (pdb)
##    for ch,dresnrest in dic.items():
##        for i,lmut in dresnrest.items():
##            print 'SEQ',i,'in PDB is',dic_pdb[ch][i]
##            try:
##                if dic[ch][i]!=dic_pdb[ch][i] or dic_impres[ch][i]==dic[ch][i]:
##                    try:
##                        dic_new_res[ch][i]=dic[ch][i]
##                    except:
##                        dic_new_res[ch]={i:dic[ch][i]}
##            except:
##                pass
##    #print 'dic is',dic
##    #print 'dic_new_res',dic_new_res
##    if len(dic_new_res)>0:
##        for ch,dii in dic_new_res.items():
##            for i,rest in dii.items():
##                print 'Assigned residue is:',rest
##        if not os.path.isfile(outpdb):
##            coot_run_rotamer_sphere_refinement_multi( input_PDB_file=pdb , input_mtz_file=mtz_phases , output_pdb=output+'.pdb' , outputCootPy=output+'.py' , dic=dic_new_res , radius_sph_ref=radius_sph_ref , coot_path=coot_path )
##            pdb=outpdb
##        if not os.path.isfile(outpdb):
##            print 'ERROR in COOT MODELING.\nFile',output_pdb,'not generated.'
##            exit()
##        print '\n\n'
##    else:
##        if not os.path.isfile(outpdb):
##            remove_Bfactor_occ_res_pdb (pdb_input=pdb,pdb_output=outpdb,dic_ch_resn=dic)
##            pdb=outpdb
##    if output_map:
##        run_phenix_polder_multiproc ( pdb , mtz_init , dic , output , polder_path=polder_path, output_mtz=True) #dic should be dic[chain][resn]
##    else:
##        run_phenix_polder_multiproc ( pdb , mtz_init , dic , output , polder_path=polder_path, output_mtz=False)

##def coot_mutate_sph_ref_overlapmap_multi (pdb,mtz_init,mtz_phases,outpdb,outlog,dic, radius_sph_ref=5, coot_path=False): #, polder_path=False , output_map=False): #Script created Dec5,2016
##    dic_pdb=return_dic_chain_resnumb_restype (pdb)
##    dic_new_res={}
##    dic_impres=return_impartial_res (pdb)
##    for ch,dresnrest in dic.items():
##        for i,lmut in dresnrest.items():
##            print 'SEQ',i,'in PDB is',dic_pdb[ch][i]
##            try:
##                if dic[ch][i]!=dic_pdb[ch][i] or dic_impres[ch][i]==dic[ch][i]:
##                    try:
##                        dic_new_res[ch][i]=dic[ch][i]
##                    except:
##                        dic_new_res[ch]={i:dic[ch][i]}
##            except:
##                pass
##    #print 'dic is',dic
##    #print 'dic_new_res',dic_new_res
##    if len(dic_new_res)>0:
##        for ch,dii in dic_new_res.items():
##            for i,rest in dii.items():
##                print 'Assigned residue is:',rest
##        if not os.path.isfile(output+'.pdb'):
##            coot_run_rotamer_sphere_refinement_multi( input_PDB_file=pdb , input_mtz_file=mtz_phases , output_pdb=output+'.pdb' , outputCootPy=output+'.py' , dic=dic_new_res , radius_sph_ref=radius_sph_ref , coot_path=coot_path )
##            pdb=output+'.pdb'
##        if not os.path.isfile(output+'.pdb'):
##            print 'ERROR in COOT MODELING.\nFile',output_pdb,'not generated.'
##            exit()
##        print '\n\n'
##    else:
##        if not os.path.isfile(output+'.pdb'):
##            remove_Bfactor_occ_res_pdb (pdb_input=pdb,pdb_output=output+'.pdb',dic_ch_resn=dic)
##            pdb=output+'.pdb'
##    sfall_overlapmap ( pdb=pdb , mtz=mtz_init , out=output )
##    #if output_map:
##        #run_phenix_polder_multiproc ( pdb , mtz_init , dic , output , polder_path=polder_path, output_mtz=True) #dic should be dic[chain][resn]
##    #else:
##        #run_phenix_polder_multiproc ( pdb , mtz_init , dic , output , polder_path=polder_path, output_mtz=False)
##

def remove_Bfactor_occ_res_pdb (pdb_input,pdb_output,dic_ch_resn,AtomsExclude=['CA','N ','C ','O '],AtomsInclude=[] ):
    ff=open(pdb_input)
    f=ff.readlines()
    ff.close()
    o=open(pdb_output,'w')
    AtomsExclude=set(AtomsExclude)
    AtomsInclude=set(AtomsInclude)
    #AtomsExclude=['CA','N ','C ','O ']
    if len( AtomsExclude.intersection(AtomsInclude) ) > 0 :
        print ('AtomsExclude and AtomsInclude has(ve)',len( AtomsExclude.intersection(AtomsInclude) ),'identical atoms, lists should contain no identical atoms.')
        print ('AtomsExclude:',AtomsExclude)
        print ('AtomsInclude',AtomsInclude)
        exit()

    for l in f:
        if l.startswith('ATOM'):
            bfChange=False
            occChange=False
            ch=l[21]
            altconf=l[16]
            resn=int(l[22:26])
            atomt=l[13:15]
            for ch2,ddic in dic_ch_resn.items():
                for resn2 in ddic:
                    if resn2==resn and ch==ch2 and atomt not in AtomsExclude and atomt in AtomsInclude:
                        bfChange=True
                        if altconf==' ':
                            occChange=True
            o.write(l[:56])
            if occChange:
                o.write('1.00 ')
            else:
                o.write(l[56:61])
            if bfChange:
                o.write('30.00')
            else:
                o.write(l[61:66])
            o.write(l[66:])
        else:
            o.write(l)
    o.close()

def return_impartial_res (pdb):
    # amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    # amino_acid_list_3L=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
    # amino_acid_list_numb_atoms=[5,6,  8  ,  9  , 11  ,  4  , 10  ,  8  ,  9  ,  8  ,  8  ,  8  ,  7  ,  9  , 11  , 6   ,  7  ,  7  ,  14 ,  12 ]
    f1=open(pdb)
    fi=f1.readlines()
    f1.close()
    resnumb=''
    restype=''
    c=0 # counter of number of atoms per residue
    chain=''
    dic_ch_resn={}
    for l in fi:
        if l.startswith('ATOM') and (l[16]==' ' or l[16]=='A') and not l[12]=='H' and not l[13]=='H': # read each line of ATOM without considering hydrogen atoms
            resnumb_new=int(l[22:26])
            if resnumb_new==resnumb:
                c+=1
            else:
                if l[21]!=chain and not c==0:
                    c-=1
                if not c==0:
                    if not amino_acid_list_numb_atoms[amino_acid_list_3L.index(restype)]==c:
                        try:
                            dic_ch_resn[chain][resnumb]=amino_acid_list[amino_acid_list_3L.index(restype)].upper()
                        except:
                            dic_ch_resn[chain]={resnumb:amino_acid_list[amino_acid_list_3L.index(restype)].upper()}
                c=1
            resnumb=int(l[22:26])
            restype=l[17:20]
            chain=l[21]
    c-=1
    if not amino_acid_list_numb_atoms[amino_acid_list_3L.index(restype)]==c:
        try:
            dic_ch_resn[chain][resnumb]=amino_acid_list[amino_acid_list_3L.index(restype)].upper()
        except:
            dic_ch_resn[chain]={resnumb:amino_acid_list[amino_acid_list_3L.index(restype)].upper()}
    return dic_ch_resn


def run_phenix_polder_multiproc ( pdb , mtz , dic , output , polder_path=False , output_mtz=False , sidechain=True): #dic should be dic[chain][resn] old: resn may be a integer or an string-> could follow phenix syntax, such as 1:5 #Script created Dec5,2016 #dic should be dic[chain][resn]
    if polder_path==False:
        polder_path='phenix.polder'
    #l='phenix.polder '+pdb+' '+mtz+' selection="chain '+ch+' resid '+str(resn)+'"'
    selection='selection="chain '
    for ch,item in dic.items():
        for ch1L in ch:
            for resn,mut in item.items():
                if selection=='selection="chain ':
                    selection+=ch1L+' and resid '+str(resn)
                else:
                    selection+=' or chain '+ch1L+' and resid '+str(resn)
                if sidechain==True:                  selection+=' and not name C and not name N  and not name O '
                elif sidechain=='OnlyMainChain': selection+=' and (name C or name N or name O or name CA)'
    selection+='"'
    if not output_mtz:
        p = subprocess.Popen([polder_path,pdb,mtz,selection], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
        #p = subprocess.Popen([polder_path,pdb,mtz,selection,'data_labels="F_XDSdataset"'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        p = subprocess.Popen([polder_path,pdb,mtz,selection,'output_file_name_prefix='+output], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
        #p = subprocess.Popen([polder_path,pdb,mtz,selection,'output_file_name_prefix='+output[:-4],'data_labels="F_XDSdataset"'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    file=open(output,'w')
    file.write(out)
    file.close()
    #os.system(l)
    if not output_mtz:
        if os.path.isfile(pdb[:-4].split('/')[-1]+'_polder_map_coeffs.mtz'): os.system('rm ./'+pdb[:-4].split('/')[-1]+'_polder_map_coeffs.mtz')
    if not os.path.isfile(output):
        print ('Failure in RJB_lib.run_phenix_polder_multiproc function, polder log file not created')
        exit()
    #os.system('mv polder.log '+output+'.log')
    #os.system('mv '+pdb[:-4]+'_polder_map_coeffs.mtz '+output+'.mtz')


def extract_CC_R_Rfree_from_polder_log (log,printtt=False): #Script created Dec5,2016
    f2=open(log)
    f=f2.readlines()
    cc12=''
    #print log
    for i,l in enumerate(f):
        if l.startswith('Map 3: real Fobs data'):
            cc12=float(f[i+1].split()[1])
            cc13=float(f[i+2].split()[1])
            cc23=float(f[i+3].split()[1])
        elif l.startswith('R factors for unmodified input model and data:'):
            ll=f[i+1]
            #rmodel=float(ll[7:13])
            #rfreemodel=float(ll[21:27])
            rmodel=float(ll.split()[2])
            rfreemodel=float(ll.split()[-1])
        elif l.startswith('R factor when ligand is excluded for mask calculation') or l.startswith('R factor for OMIT map (ligand is excluded for mask calculation):'):
            ll=f[i+1]
            #rexcl=float(ll[7:13])
            rexcl=float(ll.split()[2])
            rfreeexcl=float(ll.split()[-1])
            #rfreeexcl=float(ll[21:27])
    rimpr=(rexcl-rmodel)*100
    rfreeimpr=(rfreeexcl-rfreemodel)*100
    rsimpr=rimpr+rfreeimpr
    d={'cc13':cc13,'rmodel':rmodel,'rfreemodel':rfreemodel,'rexcl':rexcl,'rfreeexcl':rfreeexcl,'rimpr':rimpr,'rfreeimpr':rfreeimpr,'rsimpr':rsimpr}
    return d
    # if cc13>cc12 and cc13>cc23:
    #     return d
    # else:
    #     stcc12='%.2f'%(cc12)
    #     stcc13='%.2f'%(cc13)
    #     stcc23='%.2f'%(cc23)
    #     if printtt: print ('phenix.polder map ('+log+') should not be used, since the density in the OMIT region resembles to bulk solvent density, as seen by CC(1,2){'+stcc12+'} and CC(2,3){'+stcc23+'} are larger or comparable to CC(1,3){'+stcc13+'}.')
    #     return False


def sfall_overlapmap ( pdb , mtz , out ): #script created
    rmf=[]
    #writting sfall_description_file
    rmf.append(out+'_1sfall.ins')
    d1=open(out+'_1sfall.ins','w')
    d1.write('LABIN  FP=F_XDSdataset SIGFP=SIGF_XDSdataset FREE=FreeR_flag\n')
    d1.write('labout -\n   FC=FCalc PHIC=PHICalc\n')
    d1.write('MODE SFCALC -\n    XYZIN -\n    HKLIN\n')
    d1.write("symmetry 'P 65 2 2'\n")
    d1.write('badd 0.0\n')
    d1.write('vdwr 2.5\n')
    d1.write('end')
    d1.close()
    d1=open(out+'_1sfall.ins', 'r')
    #create map1
    #running sfall
    sfall = subprocess.Popen([ 'sfall','HKLIN',mtz,'HKLOUT',out+'_1sfall.mtz','XYZIN',pdb], stdin=d1, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    ou, err = sfall.communicate()
##    print 'sfall1'
##    print ou
##    print '\n\n\n\n\n\n\n'
    d1.close()
    rmf.append(out+'_1sfall.mtz')
    #os.system('rm '+file_output_map+'_fft_job.txt')
    rmf.append(out+'_2fft.ins')
    d2=open(out+'_2fft.ins','w')
    d2.write('labin  F1=F_XDSdataset SIG1=SIGF_XDSdataset PHI=PHICalc\nEND')
    d2.close()
    d2=open(out+'_2fft.ins','r')
    fft = subprocess.Popen([ 'fft','HKLIN',out+'_1sfall.mtz','MAPOUT',out+'_2fft.map'], stdin=d2, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    ou, err = fft.communicate()
    d2.close()
    rmf.append(out+'_2fft.map')
##    print 'fft2'
##    print ou
##    print '\n\n\n\n\n\n\n'
    #map mask
    rmf.append(out+'_3mapmask.ins')
    d3=open(out+'_3mapmask.ins','w')
    d3.write('XYZLIM ASU')
    d3.close()
    d3=open(out+'_3mapmask.ins','r')
    mapmask = subprocess.Popen([ 'mapmask','MAPIN',out+'_2fft.map','MAPOUT',out+'_3mapmask.map'], stdin=d3, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    ou, err = mapmask.communicate()
    d3.close()
    rmf.append(out+'_3mapmask.map')
##    print 'mapmask3'
##    print ou
##    print '\n\n\n\n\n\n\n'
    #create map2
    #writting sfall_description_file
    rmf.append(out+'_4sfall.ins')
    d4=open(out+'_4sfall.ins','w')
    d4.write('MODE ATMMAP\n')
    d4.write('grid 108 108 156\n')
    d4.write("symmetry 'P 65 2 2'\n")
    d4.write('badd 0.0\n')
    d4.write('vdwr 2.5\n')
    d4.write('end')
    d4.close()
    d4=open(out+'_4sfall.ins', 'r')
    sfall2 = subprocess.Popen([ 'sfall','MAPOUT',out+'_4model_sfall.map','XYZIN',pdb], stdin=d4, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    ou, err = sfall2.communicate()
    d4.close()
    rmf.append(out+'_4model_sfall.map')
##    print 'sfall4'
##    print ou
##    print '\n\n\n\n\n\n\n'
    #map mask
    rmf.append(out+'_5mapmask.ins')
    d5=open(out+'_5mapmask.ins','w')
    d5.write('XYZLIM 0 54 0 107 0 51\n')
    d5.write('AXIS Y X Z')
    d5.close()
    d5=open(out+'_5mapmask.ins','r')
    mapmask2 = subprocess.Popen([ 'mapmask','MAPIN',out+'_4model_sfall.map','MAPOUT',out+'_5mapmask.map'], stdin=d5, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    ou, err = mapmask2.communicate()
    d5.close()
    rmf.append(out+'_5mapmask.map')
##    print 'mapmask5'
##    print ou
##    print '\n\n\n\n\n\n\n'
    #sfall map3
    rmf.append(out+'_6sfall.ins')
    d6=open(out+'_6sfall.ins','w')
    d6.write('MODE ATMMAP RESMOD\n')
    d6.write('grid 108 108 156\n')
    d6.write("symmetry 'P 65 2 2'\n")
    d6.write('chain A 1 121\n')
    d6.write('badd 0.0\n')
    d6.write('vdwr 2.5\n')
    d6.write('end')
    d6.close()
    d6=open(out+'_6sfall.ins', 'r')
    sfall3 = subprocess.Popen([ 'sfall','MAPOUT',out+'_6model_sfall.map','XYZIN',pdb], stdin=d6, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    ou, err = sfall3.communicate()
    d6.close()
    rmf.append(out+'_6model_sfall.map')
##    print 'sfall6'
##    print ou
##    print '\n\n\n\n\n\n\n'
    #mapmask map3
    rmf.append(out+'_7mapmask.ins')
    d7=open(out+'_7mapmask.ins','w')
    d7.write('XYZLIM 0 54 0 107 0 51\n')
    d7.write('AXIS Y X Z')
    d7.close()
    d7=open(out+'_7mapmask.ins','r')
    mapmask3 = subprocess.Popen([ 'mapmask','MAPIN',out+'_6model_sfall.map','MAPOUT',out+'_7mapmask.map'], stdin=d7, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    ou, err = mapmask3.communicate()
    d7.close()
    rmf.append(out+'_7mapmask.map')
##    print 'mapmask6'
##    print ou
##    print '\n\n\n\n\n\n\n'
    #correlate 3 maps with overlapmap
    rmf.append(out+'_8overlapmap.ins')
    d8=open(out+'_8overlapmap.ins','w')
    d8.write('correlate residue\nchain A 1 121\nEND')
    d8.close()
    d8=open(out+'_8overlapmap.ins','r')
    overlapmap = subprocess.Popen([ 'overlapmap','MAPIN1',out+'_3mapmask.map','MAPIN2',out+'_5mapmask.map','MAPIN3',out+'_7mapmask.map'], stdin=d8, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    ou, err = overlapmap.communicate()
##    print 'overlapmap'
##    print ou
    log=open(out+'.log','w')
    log.write(ou)
    log.close()
    for f in rmf:
        os.system('rm '+f)
   
def extract_CC_overlapmaplog_to_dic_ch (log):
    ff=open(log)
    f=ff.readlines()
    ff.close()
    d={}
    for i,l in enumerate(f):
        if l.startswith(' $TABLE:  Correlation Residue by Residue (chain'):
            ch=l[l.index(')')-1]
            d[ch]={}
            ii=i+5
            ll=f[ii]
            while ll!='  $$\n':
                #print ll
                lls=ll.split()
                resn=int(lls[0])
                ccmc=float(lls[1])
                ccsc=float(lls[2])
                d[ch][resn]={'ccmc':ccmc,'ccsc':ccsc}
                ii+=1
                ll=f[ii]
    return d
    
    
    
    
    
    
def generate_from_refmac_ROTAMER_file_list_frequency_sequences (input_file_with_list , number_models ): #Script modificated in 30th July 2015 to extract frequency of amino acids of chosen # of models from _refmac_ROTAMER table
    input_file=open(input_file_with_list,"r")
    input_file_list=input_file.readlines()
    evaluate_sequences=[]
    list_all=[]
    for line_number in range(number_models):
        evaluate_sequences.append(input_file_list[line_number+1].split()[7])
    length_seq=len(evaluate_sequences[0])
    for i in range(length_seq):
        list_by_seqnumb=[]
        var_seq=""
        for sequence in evaluate_sequences:
            var_seq+=sequence[i]
        for amino_acid in amino_acid_list:
            var_list_by_aa=[]
            if var_seq.count(amino_acid)>0:
                var_list_by_aa=[amino_acid,var_seq.count(amino_acid)]
                list_by_seqnumb.append(var_list_by_aa)
        list_all.append(list_by_seqnumb)
    return list_all



def given_list_of_possible_res_write_seqfile_for_SLIDER_POINTMUTATION ( list_of_residues , output_filename ):
    output_file=open(output_filename,"w")

    i=0
    for item in list_of_residues:
        if i<len(item):
            i=len(item)

    for line_output in range(i):
        print ("line",i)
        for seq_number in range (len(list_of_residues)):
    ##        print "seq_number:",seq_number
    ##        print "with content",seq[seq_number]
            try:
                residue=list_of_residues[seq_number][line_output]
                if len(residue)==3:
                    residue=amino_acid_list [amino_acid_list_3L.index(residue)]
                output_file.write(residue)
    ##            print "residue_output:",residue
            except:
                output_file.write('-')
        output_file.write('\n')
    output_file.close()
    
def given_seqfile_SLIDER_POINTMUTATION_and_list_from_frequencies_refmac_ROTAMER_generate_star_fix_residue ( list_of_residues_refmac , input_seqfilename , output_seqfilename ):
    input_seqfile=open(input_seqfilename,"r")
    input_seqfile_list=input_seqfile.readlines()
    input_seqfile.close()
    output_seqfile=open(output_seqfilename,"w")
    for line in range(len(input_seqfile_list)):
        if not line==1:
            output_seqfile.write(input_seqfile_list[line])
        else:
            for seq_number in range(len(input_seqfile_list[0])-1):
                if len(list_of_residues_refmac[seq_number])==1 and input_seqfile_list[0][seq_number]==list_of_residues_refmac[seq_number][0][0] and input_seqfile_list[1][seq_number]=="-":
                    output_seqfile.write(".")
                else:
                    output_seqfile.write(input_seqfile_list[1][seq_number])
            if line!=len(input_seqfile_list): #condition to prevent from writting a \n in end of file if seqfile has only two lines...
                output_seqfile.write("\n")
                
                
def retrieve_information_R_CCampl_from_PDB ( dictionary_var , dictionary_overall ) : #copied from SLIDER_POINTMUTATION.py September 13rd 2015
    xyzout=dictionary_var['folder']+'/'+dictionary_var['file']
    refmac_R,refmac_Rfree=retrieve_Rfactor_from_PDB_file ( xyzout )
    refmac_CC_ampl,refmac_CC_ampl_free=retrieve_CC_ampl_from_PDB_file ( xyzout )
    number_res=count_number_of_residues ( xyzout )
    dictionary_var.update({ 'R':refmac_R , 'Rfree':refmac_Rfree , 'CC_ampl':refmac_CC_ampl , 'CC_ampl_free':refmac_CC_ampl_free , '#res':number_res })
    dictionary_overall[ dictionary_var['file'] ]=dictionary_var
    return dictionary_overall
    
def retrieve_CC_from_HTML_buster ( dictionary_var , dictionary_overall ) : #copied from SLIDER_POINTMUTATION.py September 13rd 2015
    htmlout=dictionary_var['folder']+'/'+dictionary_var['file'].split('.pd')[0]+'.html'
    CCmc,CCsc=retrieve_busterCC_from_HTML_file ( htmlout )
    dictionary_overall[ dictionary_var['file'] ]['CCmc']=CCmc
    dictionary_overall[ dictionary_var['file'] ]['CCsc']=CCsc
    return dictionary_overall


def retrieve_Rfactor_from_PDB_file ( PDB_file ): #copied from SLIDER_POINTMUTATION.py September 13rd 2015
    R_factor, Rfree = '',''
    PDB_Rfactor_file=open ( PDB_file , 'r' )
    PDB_Rfactors_list=PDB_Rfactor_file.readlines()
    PDB_Rfactor_file.close()
    skipp=False
    for PDB_file_Rfac in range( len(PDB_Rfactors_list) ):
        if PDB_Rfactors_list[PDB_file_Rfac].startswith('REMARK') and len( PDB_Rfactors_list[PDB_file_Rfac].split() )>=5 and not skipp:
            if PDB_Rfactors_list[PDB_file_Rfac].split()[2]=='R' and PDB_Rfactors_list[PDB_file_Rfac].split()[5]=='SET)' :
                R_factor=float( PDB_Rfactors_list[PDB_file_Rfac].split()[7] )
            elif PDB_Rfactors_list[PDB_file_Rfac].split()[2]=='FREE' and PDB_Rfactors_list[PDB_file_Rfac].split()[5]==':' :
                Rfree=float ( PDB_Rfactors_list[PDB_file_Rfac].split()[6] )
            elif PDB_Rfactors_list[PDB_file_Rfac].startswith('REMARK best refinement for ') and 'R/Rfree' in PDB_Rfactors_list[PDB_file_Rfac]:
                R_factor=float(PDB_Rfactors_list[PDB_file_Rfac][-14:-8])
                Rfree  =float(PDB_Rfactors_list[PDB_file_Rfac][-7:-1])
                skipp=True

    if R_factor=='' and Rfree=='':
        print ('Failure obtaining R/Rfree of PDB',PDB_file,'EXITING.')
        quit()
    else: return R_factor,Rfree


def retrieve_CC_ampl_from_PDB_file ( PDB_file ):
##Tested in BUSTER and REFMAC refinements PDB output
    PDB_CC_ampl_file=open ( PDB_file , 'r' )
    PDB_CC_ampl_list=PDB_CC_ampl_file.readlines()
    PDB_CC_ampl_file.close()
    for PDB_CC_ampl_ in range( len(PDB_CC_ampl_list) ):
        if PDB_CC_ampl_list[PDB_CC_ampl_].startswith('REMARK') and len( PDB_CC_ampl_list[PDB_CC_ampl_].split() )>=6:
            if PDB_CC_ampl_list[PDB_CC_ampl_].split()[2]=='CORRELATION' and PDB_CC_ampl_list[PDB_CC_ampl_].split()[5]==':':
                CC_ampl_float=float( PDB_CC_ampl_list[PDB_CC_ampl_].split()[6] )
            elif PDB_CC_ampl_list[PDB_CC_ampl_].split()[2]=='CORRELATION' and PDB_CC_ampl_list[PDB_CC_ampl_].split()[5]=='FREE':
                CC_ampl_free_float=float ( PDB_CC_ampl_list[PDB_CC_ampl_].split()[7] )
    try: return CC_ampl_float, CC_ampl_free_float
    except:
        print ('Failure retrieving CC_ampl from',PDB_file,'assuming -1 as value to prevent failure')
        CC_ampl_float=-1.0
        CC_ampl_free_float=-1.0
        return CC_ampl_float,CC_ampl_free_float


def count_number_of_residues_not_good ( folder_PDB ) : #including waters
    ProteinPDB = Bio.PDB.PDBParser()
    ProteinPDB2 = ProteinPDB.get_structure('whatever', folder_PDB )
    chain_list = Selection.unfold_entities(ProteinPDB2, 'C')
    number_res = 0
    for b in range(len(chain_list)):
        res_chain=Selection.unfold_entities(chain_list[b], 'R')
        number_res=number_res+len(res_chain)
    return number_res

def refine_buster_multiprocesses ( pdb_input , mtz_input , pdb_output , buster_path , buster_options='-noWAT -nbig 10 -RB -nthread 1  UsePdbchk="no"' , list_keep=['pdb','mtz','tar.gz','log','html']):
##def refine_buster_multiprocesses ( dictionary_buster , mtz_input ):
#(pdb_input , mtz_input , folder_output , log_output):
#dictio_buster = {}
#dictio_buster = {'file':PDB_file , 'folder':diretorio, 'last_folder':diretorio.split('/')[len(diretorio.split('/'))-1 ] }
##    file_name=dictionary_buster['file'].split('.pd')[0]
##    folder_output=dictionary_buster['folder']+'/BUSTER_temp_'+file_name
##    xyzin=dictionary_buster['folder']+'/'+dictionary_buster['file']
##    os.system( 'mkdir ' + folder_output )
    xyzin=pdb_input
    folder=pdb_output[:-4]
##    os.system('mkdir '+folder)
    #buster_options='-noWAT -nbig 10 -RB -nthread 1 UsePdbchk="no"'
     #-w 50 AdjustXrayWeightAutomatically="no" -autoncs 
    p = subprocess.Popen([buster_path,'-p',xyzin,'-m',mtz_input,'-d',folder]+buster_options.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
    out, err = p.communicate()
    if 'log' in list_keep:
        file=open(folder+'.log','w')
        file.write(out)
        file.close()
    #rename_compress_remove_BUSTER_jobs ( folder_output , dictionary_buster['folder'] , file_name )
    rename_chosen_outputname_compress_remove_BUSTER_jobs ( folder , pdb_output[:-4] , list_keep)
    if not os.path.isfile( pdb_output ):
        print ('\n\n\nERROR in refinement.\nFile',pdb_output,'not generated.\n\n\n')
        exit()

def sh_refine_buster_multiprocesses ( pdb_input , mtz_input , pdb_output , buster_path , sh_file , distribute_computing , buster_options='-noWAT -nbig 10 -RB -nthread 1  UsePdbchk="no"' ):
    file_sh=open(sh_file,'w')
    if distribute_computing=='multiprocessing':
        file_sh.write('#!/bin/bash'+ '\n')
    if distribute_computing=='local_grid':
        file_sh.write('#!/bin/tcsh'+ '\n')
        file_sh.write('setenv PATH /usr/local/bin\nsetenv PATH /usr/bin:$PATH\nsetenv PATH /bin:$PATH\n') # suggested by MAX
        file_sh.write('hostname > '+sh_file.split('/')[-1]+ '_host \n')
##        file_sh.write('setenv LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu\n')
##        file_sh.write('setenv PATH /usr/local/bin\n')
##        file_sh.write('setenv PATH /usr/bin:$PATH\n')
##        file_sh.write('setenv PATH /bin:$PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MCR/lib/:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/runtime/glnxa64:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/bin/glnxa64:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/sys/os/glnxa64:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/sys/java/jre/glnxa64/jre/lib/amd64/server:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/saxs/MATLAB/MATLAB_Compiler_Runtime/v82/sys/java/jre/glnxa64/jre/lib/amd64/:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/atsas/atsas/lib/atsas:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/Xmipp/Xmipp-2.4-x64/lib:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/imod/imod-4.5.3-RH6-64/IMOD/lib:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/imp/imp/lib64:$LD_LIBRARY_PATH\n')
##        file_sh.write('setenv LD_LIBRARY_PATH /xtal/XtalView/lib/ibmpcLinux2:$LD_LIBRARY_PATH\n')
##        file_sh.write('source /xtal/xtalsetup\n')
##        file_sh.write('source /xtal/ccp4/ccp4/include/ccp4.setup\n')
        file_sh.write('source /xtal/ccp4/ccp4-7.0/include/ccp4.setup\n')
        file_sh.write('source /xtal/buster/buster/setup.csh\n')
        #file_sh.write('source /xtal/xtalsetup'+ '\n')
        
    xyzin=pdb_input
    folder=pdb_output[:-4]
    if distribute_computing=='local_grid':
        xyzin =xyzin.split('/')[-1]
        folder=folder.split('/')[-1]
        mtz_input=mtz_input.split('/')[-1]
##    os.system('mkdir '+folder)
    if buster_options==False:
    #buster_options=' -noWAT -nbig 10 -RB -nthread 1 "StopOnGellySanityCheckError=no"'
        buster_options=' -noWAT -nbig 10 -RB -nthread 1 UsePdbchk="no" '
    #-w 50 AdjustXrayWeightAutomatically="no" -autoncs
    #check if file has arrived
    #file_sh.write('if (! -f "'+xyzin+'" ) then\n        echo "Files did not arrive yet, waiting 15 seconds."\n        sleep 15\n        endif\n')
    file_sh.write(buster_path+' -p '+xyzin+' -m '+mtz_input+' -d '+folder+' '+buster_options)
    file_sh.write(' > '+folder+'.log\n')

    #p = subprocess.Popen([buster_path,'-p',xyzin,'-m',mtz_input,'-d',folder]+buster_options.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #out, err = p.communicate()
    #file=open(folder+'.log','w')
    #file.write(out)
    #file.close()
    #rename_compress_remove_BUSTER_jobs ( folder_output , dictionary_buster['folder'] , file_name )
    folder_input_BUSTER=folder
    filename_output_without_extension=pdb_output[:-4]
    if distribute_computing=='local_grid':
        filename_output_without_extension=filename_output_without_extension.split('/')[-1]

    #added for new version
    # file:///xtal/buster/BUSTER_snapshot_20170920/docs/autobuster/manual/autoBUSTER7.html#corr
    #general run line: % corr -p refine.pdb -m refine.mtz -F 2FOFCWT -P PH2FOFCWT
    file_sh.write(buster_path[:buster_path.rindex('/')]+'/corr'+' -p '+folder+'/refine.pdb -m '+folder+'/refine.mtz -F 2FOFCWT -P PH2FOFCWT -d '+folder+'_RSCC  > '+folder+'_RSCC.log\n')

    #recent BUSTER refine does not output html file
    #file_sh.write( 'cp ' + folder_input_BUSTER + '/analyse.html ' + filename_output_without_extension + '.html' + '\n')
    file_sh.write( 'cp ' + folder_input_BUSTER + '/refine.pdb ' + filename_output_without_extension + '.pdb' + '\n')
    file_sh.write( 'cp ' + folder_input_BUSTER + '/refine.mtz ' + filename_output_without_extension + '.mtz' + '\n')

    file_sh.write( 'tar -zcf ' + filename_output_without_extension  + '.tar.gz ' + folder_input_BUSTER + ' '+folder_input_BUSTER+'_RSCC \n')
    file_sh.write( 'rm -r ' + folder_input_BUSTER + ' '+folder_input_BUSTER+'_RSCC  \n\n')
    #rename_chosen_outputname_compress_remove_BUSTER_jobs ( folder , pdb_output[:-4] )
    file_sh.close()
    os.system('chmod 755 '+sh_file)
    os.system('chmod -R 755 '+pdb_output[:pdb_output.rindex('/')])
##    f=open(sh_file)
##    f2=f.read()
##    print f2
##    f.close()

def refine_buster_singleprocess ( pdb_input , mtz_input , water):
    folder_output="folder_"+pdb_input[:-4]+"_folder"
    xyzin=pdb_input
    os.system( 'mkdir ' + folder_output )
    if water:
        buster_options=' -nbig 10 -RB -WAT 1 '
    else:
        buster_options=' -noWAT -nbig 10 -RB '
    #-w 50 AdjustXrayWeightAutomatically="no" -autoncs 
    #os.system("/xtal/buster/buster_snapshot_20150316/autoBUSTER/bin/linux64/refine -p "+xyzin+' -m '+mtz_input+' -d '+folder_output+buster_options+" > "+pdb_input[:-4]+"_BUSTER.log")
    os.system("refine -p " + xyzin + ' -m ' + mtz_input + ' -d ' + folder_output + buster_options + " > " + pdb_input[:-4] + "_BUSTER.log")
    rename_compress_remove_BUSTER_jobs ( folder_output , "." , pdb_input[:-4] )


def rename_compress_remove_BUSTER_jobs ( folder_input_BUSTER , folder_output_files , filename_PDB_mtz_log_compression_without_extension ):
    #new version does not output .html anymore
    #os.system( 'cp ' + folder_input_BUSTER + '/analyse.html ' + folder_output_files + '/' + filename_PDB_mtz_log_compression_without_extension + '_BUSTER.html' )
    os.system( 'cp ' + folder_input_BUSTER + '/refine.pdb ' + folder_output_files + '/' + filename_PDB_mtz_log_compression_without_extension + '_BUSTER.pdb' )
    os.system( 'cp ' + folder_input_BUSTER + '/refine.mtz ' + folder_output_files + '/' + filename_PDB_mtz_log_compression_without_extension + '_BUSTER.mtz' )
    os.system( 'tar -zcf ' + folder_output_files + '/' + filename_PDB_mtz_log_compression_without_extension  + '_BUSTER.tar.gz ' + folder_input_BUSTER )
    os.system( 'rm -r ' + folder_input_BUSTER )

def rename_chosen_outputname_compress_remove_BUSTER_jobs ( folder_input_BUSTER , filename_output_without_extension , list_keep=['pdb','mtz','tar.gz','log']): #,'html']):
    #if 'html' in list_keep and os.path.isfile(folder_input_BUSTER + '/analyse.html'):
    #    os.system( 'cp ' + folder_input_BUSTER + '/analyse.html ' + filename_output_without_extension + '.html' )
    if 'pdb' in list_keep:
        os.system( 'cp ' + folder_input_BUSTER + '/refine.pdb ' + filename_output_without_extension + '.pdb' )
    if 'mtz' in list_keep:
        os.system( 'cp ' + folder_input_BUSTER + '/refine.mtz ' + filename_output_without_extension + '.mtz' )
    if 'tar.gz' in list_keep:
        os.system( 'tar -zcf ' + filename_output_without_extension  + '.tar.gz ' + folder_input_BUSTER )
    os.system( 'rm -r ' + folder_input_BUSTER )


#obsolete
# def retrieve_busterCC_from_HTML_file ( HTML_file ):
#     HTML_file_BUSTER_CCs  = open ( HTML_file , 'r' )
#     HTML_file_BUSTER_CCs_list = HTML_file_BUSTER_CCs.readlines()
#     HTML_file_BUSTER_CCs.close()
#     return float ( HTML_file_BUSTER_CCs_list[len(HTML_file_BUSTER_CCs_list)-11].split('&nbsp;')[4] ) , float ( HTML_file_BUSTER_CCs_list[len(HTML_file_BUSTER_CCs_list)-11].split('&nbsp;')[len(HTML_file_BUSTER_CCs_list[len(HTML_file_BUSTER_CCs_list)-11].split('&nbsp;'))-2].split('<')[0] )


def retrieve_busterCC_from_log_file ( log_file ): # old, I dont know why but this line '       Total' wasnt outputting in the run in: /localdata/rafael/SLIDER_phasing/1YZF/5-Article/control-seqs/control-1YZF_stats/1_ref/c0_true_ref1_RSCC.log
    #print log_file
    LOG_file_BUSTER_CCs  = open ( log_file , 'r' )
    LOG_file_BUSTER_CCs_list = LOG_file_BUSTER_CCs.readlines()
    LOG_file_BUSTER_CCs.close()
    CCmc = 0
    CCsc = 0
    i=0
    while CCmc == 0 and i<=len(LOG_file_BUSTER_CCs_list):
        if LOG_file_BUSTER_CCs_list[i].startswith('   CC: main-chain = '):
            l=LOG_file_BUSTER_CCs_list[i].split()
            CCmc = float ( l[-4] )
            CCsc = float(l[-1])
    # while CCmc == 0 or i<=len(LOG_file_BUSTER_CCs_list):
    #     if LOG_file_BUSTER_CCs_list[i].startswith('       Total'):
    #         l=LOG_file_BUSTER_CCs_list[i].split()
    #         CCmc = float ( l[6] )
    #         CCsc = float(l[-1])
        else: i+=1
    if CCmc == 0:
        print ('Condition not found in',log_file)
        print ('Failure in function RJB_lib.retrieve_busterCC_from_log_file')
        exit()
    return CCmc , CCsc

def retrieve_busterRSCC_from_log_file ( log_file ): #update
    LOG_file_BUSTER_RSCCs  = open ( log_file , 'r' )
    LOG_file_BUSTER_RSCCs_list = LOG_file_BUSTER_RSCCs.readlines()
    LOG_file_BUSTER_RSCCs.close()
    RSCCmc = 0
    RSCCsc = 0
    i=0
    while RSCCmc == 0:
        if LOG_file_BUSTER_RSCCs_list[i].startswith('   CC: main-chain =') or LOG_file_BUSTER_RSCCs_list[i].startswith('       Total'):
            l=LOG_file_BUSTER_RSCCs_list[i].split()
            RSCCmc = float ( l[-4] )
            RSCCsc = float(l[-1])
        # if LOG_file_BUSTER_RSCCs_list[i].startswith('       Total'):
            #l=LOG_file_BUSTER_RSCCs_list[i].split()
            # RSCCmc = float ( l[6] )
            # RSCCsc = float(l[-1])
        else: i+=1
    return RSCCmc , RSCCsc



def separate_PDB_file_by_chains (pdb_input):
    chain_list,seq,aa_number=extract_protein_chainID_res_number (pdb_input)
    file=open(pdb_input,'r')
    file_list=file.readlines()
    file.close()
    print ("PDB file:",pdb_input)
    for chain in chain_list:
        output_by_chain=open(pdb_input[:-4]+"_"+chain+".pdb","w")
        for line in file_list:
            if line.startswith('ATOM') or line.startswith("ANISOU") or line.startswith("HETATM"):
                line2=line.split()
                if chain == line2[4]:
                    output_by_chain.write(line)
        else:
            output_by_chain.write(line)
        output_by_chain.close()
        print ("PDB file: ",pdb_input[:-4]+"_"+chain+".pdb","containing chain",chain,"has just been written.")
                    

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def remove_empty_keys_in_dictionary (dictio):
    remove_keys=[]
    for key in dictio:
        if len(dictio[key])==0:
            remove_keys.append(key)
    for key in remove_keys:
        del dictio[key]
        print (key,"removed from dictionary")


def separate_discontinuous_residues_in_chains (pdb_input,pdb_output):
    file=open(pdb_input,'r')
    file_list=file.readlines()
    file.close()
    pdb_output_file=open(pdb_output,"w")
    dict_missing_lines={}
    chains_existing=[]
    #outputing lines that are not related to ATOMS and inserting ATOM lines into lists organized by chain in a dictionary
    for line in file_list:
        if not line.startswith('ATOM') and not line.startswith('ANISOU') and not line.startswith('END') and not line.startswith('TER'):
            pdb_output_file.write(line)
        elif line.startswith('ATOM'):
            chain=line[21]
            #condition to creat list if it does not exist
            try:
                dict_missing_lines[chain].append(line)
            except:
                dict_missing_lines[chain]=[line]
    #outputing consecutive residues in all chains in dictionary, only first consecutive residue will be written, when gap is encontered, different chain will be assigned
##    print "Writting consecutive residues"
    for chain in dict_missing_lines:
        res_numb_check=int(dict_missing_lines[chain][0][22:26])
        for line in dict_missing_lines[chain]:
            if res_numb_check==int(line[22:26]) or res_numb_check+1==int(line[22:26]):
                res_numb_check=int(line[22:26])
                pdb_output_file.write(line)
##                print line
                #dict_missing_lines[chain].remove(line)
                #removing written lines
                index_line=dict_missing_lines[chain].index(line)
                dict_missing_lines[chain]=dict_missing_lines[chain][:index_line]+dict_missing_lines[chain][index_line+1:]
                #add chain into a list to prevent from using same chain letter
                if not chain in chains_existing:
                    chains_existing.append(chain)
##                    print chain,"appended"
            #not consecutive residues, will be further written
    #remove empty list in dictionary[chain], so as program runs, things go empty
    remove_empty_keys_in_dictionary (dict_missing_lines)
    #condition to run until dictio is empty
##    print "Writting not consecutive residues"
    while len(dict_missing_lines)>0:
        for chain in dict_missing_lines:
            #chain_options = chains that were not assigned to something else
            chain_options=[]
            for chain2 in alphabet:
                if chain2 not in chains_existing:
                    chain_options.append(chain2)
            res_numb_check=int(dict_missing_lines[chain][0][22:26])
            for line in dict_missing_lines[chain]:
                if res_numb_check==int(line[22:26]) or res_numb_check+1==int(line[22:26]):
                    res_numb_check=int(line[22:26])
                    pdb_output_file.write(line[:21]+chain_options[0]+line[22:])
##                    print line[:21]+chain_options[0]+line[22:]
                    #removing written lines
                    index_line=dict_missing_lines[chain].index(line)
                    dict_missing_lines[chain]=dict_missing_lines[chain][:index_line]+dict_missing_lines[chain][index_line+1:]
                    #dict_missing_lines[chain].remove(line)
                    if not chain_options[0] in chains_existing:
                        chains_existing.append(chain_options[0])
        remove_empty_keys_in_dictionary (dict_missing_lines)
            





##                res_numb_check
##                chain_check
##                if ( res_numb_check==int(line[22:26]) or res_numb_check==int(line[22:26])+1 ) and chain_check==line[21]:
##                    pdb_output_file.write(line)
##                    res_numb_check=int(line[22:26])
##                    chain_check=line[21]
##                    if not chain_check in chains_existing:
##                        chains_existing.append(chain_check)
##                else:
##                    try:
##                        dict_missing_lines[chain_check].append(line)
##                    except:
##                        dict_missing_lines[chain_check]=[line]
##            except:
##                res_numb_check=int(line[22:26])
##                chain_check=line[21]
##                pdb_output_file.write(line)
##    for chain in dict_missing_lines:
##        if chain not in chains_existing:
##            for line in dict_missing_lines[chain]:
##                
##    for missing_line in missing_lines:
##        if missing_line[21] not in chains:
                
                
#Trying to write something to use chain_list,seq,aa_number=extract_protein_chainID_res_number info to keep same lines that do not need to be changed
##    chain_list,seq,aa_number=extract_protein_chainID_res_number (pdb_input)
##    letter_for_chains=[]
##    for chain in alphabet:
##        if chain not in chain_list:
##            letter_for_chains.append(chain)
##    for line in file_list:
##        if not line.startswith('ATOM') or not line.startswith('ANISOU') or not line.startswith('TER') or not line.startswith('END'):
##            pdb_output_file.write(line)
##    for chain_number in range(len(chain_list)):
##        if not "/" in seq[chain_number]:
##            for line in file_list:
##                if line.startswith('ATOM') and chain_list[chain_number]==line[21]:
##                    pdb_output_file.write(line)
##        else:
##            res_numb_check=None
##            list_seq_no_gaps=seq[chain_number].split("/")
##            list_size_seq_no_gaps=[]
##            total_size=0
##            for i in list_seq_no_gaps:
##                list_size_seq_no_gaps.append(len(i))
##            for line in file_list:
##                if line.startswith('ATOM') and chain_list[chain_number]==line[21]:
##                    if res_numb_check==None:
##                        res_numb_check=int(line[22:26])
##                    if res_numb_check==int(line[22:26]) and total_size<list_size_seq_no_gaps[0]:
##                        pdb_output_file.write(line)
##                    elif res_numb_check==int(line[22:26])+1 and not total_size<list_size_seq_no_gaps[0]:
##                        pdb_output_file.write(line)
##                        total_size+=1
##                        res_numb_check=int(line[22:26])
##                    try:
##                    pdb_output_file.write(line)

            

#Trying to write something that would iterate over line and it would take into account all possibilities
##    chain_existing=[]
##    res=[]
##    number=0
###    res_numb_check=None
##    chain_check=None
###    res_string=''
##    for line in file_list:
##        if not line.startswith('ATOM') or not line.startswith('ANISOU'):
##            pdb_output_file.write(line)
##        elif line.startswith('ATOM') :
###           condition for first line of ATOM to take starting conditions
##            if chain_check==None:
##                res_numb_check=int(line[22:26])
##                chain_check=line[21]
##                chain_existing.append(chain_check)
###                residue=amino_acid_list[ amino_acid_list_3L.index(line[17:20]) ]
###                res_string=res_string+residue
##                pdb_output_file.write(line)
##            else:
##                #condition for ATOMs from same chain
##                if chain_check==line[21]:
##                    #condition for ATOMs from same chain and same res numb
##                    if res_numb_check==int(line[22:26]) or res_numb_check==int(line[22:26])+1:
##                        pdb_output_file.write(line)
##                        chain_check==line[21]
##                        res_numb_check=int(line[22:26])
##                    #condition for ATOMs from same chain and different res numb
##                    else:
##                        chain_new=[]
##                        for chain in alphabet:
##                            if not chain in chain_existing:
##                                chain_new.append(chain)
##                        chain_check=chain_new[0]
##                        chain_existing.append(chain_check)
##                        res_numb_check=int(line[22:26])
##                        pdb_output_file.write(line[:21]+chain_check+line[22:])
##                #condition for ATOMs from different chain
##                else:
##                    pdb_output_file.write(line)
##                    
##                    if res_numb_check==int(line[22:26]) or res_numb_check==int(line[22:26])+1:
##                        pdb_output_file.write(line)
##                    res_numb_check==int(line[22:26]) and :
##                    pdb_output_file.write(line)
##                else:
##                    print "blah"
##                if not res_check==int(line[22:26]):
##                    if chain_check==line[21]:
##                        if int(line[22:26])==res_check+1:
##                            res_check=res_check+1
###                            residue=amino_acid_list[ amino_acid_list_3L.index(line[17:20]) ]
###                            res_string=res_string+residue
##                            pdb_output_file.write(line)
##                        else:
##                            res_check=int(line[22:26])
###                            residue=amino_acid_list[ amino_acid_list_3L.index(line[17:20]) ]
###                            res_string=res_string+'/'+residue
##                            
##                    else:
##                        res.append(res_string)
##                        res_string=''
##                        chain_check=line[21]
##                        chain_ID.append(chain_check)
##                        residue=amino_acid_list[ amino_acid_list_3L.index(line[17:20]) ]
##                        res_string=residue
##                        res_check=int(line[22:26])
##    res.append(res_string)
##    count=0
##    for chain in res:
##        count=count+len(chain.replace('/',''))
##    return chain_ID , res , count



##    chain_list,seq,aa_number=extract_protein_chainID_res_number (pdb_input)
##    file=open(pdb_input,'r')
##    file_list=file.readlines()
##    file.close()
##    pdb_output_file=open(pdb_output,"w")
##    print "PDB file:",pdb_input
##    for line in file_list:
##    for chain in chain_list:
##        output_by_chain=open(pdb_input[:-4]+"_"+chain+".pdb","w")
##        for line in file_list:
##            if line.startswith('ATOM') or line.startswith("ANISOU") or line.startswith("HETATM"):
##                line2=line.split()
##                if chain == line2[4]:
##                    output_by_chain.write(line)
##        else:
##            output_by_chain.write(line)
##        output_by_chain.close()
##        print "PDB file: ",pdb_input[:-4]+"_"+chain+".pdb","containing chain",chain,"has just been written."



def PDB_measure_distance (pdb_input, residue_number1 , atom1 , chain1 , residue_number2 , atom2 , chain2 ):
    file=open(pdb_input,'r')
    file_list=file.readlines()
    file.close()
##    print "PDB file:",pdb_input
##    print "res1 is",residue_number1
##    print "res2 is",residue_number2
##    print "atom1 is",atom1
##    print "atom2 is",atom2
##    print "chain1 is",chain1
##    print "chain2 is",chain2
    for line in file_list:
        if line.startswith('ATOM') and line[13:16].replace(" ","")==atom1 and int(line[22:26])==residue_number1 and line[21]==chain1:
            X1=float(line[31:38])
            Y1=float(line[39:46])
            Z1=float(line[47:54])
        elif line.startswith('ATOM') and line[13:16].replace(" ","")==atom2 and int(line[22:26])==residue_number2 and line[21]==chain2:
            X2=float(line[31:38])
            Y2=float(line[39:46])
            Z2=float(line[47:54])
    try:
        list1=[X1,Y1,Z1]
        list2=[X2,Y2,Z2]
##        print list1
##        print list2
        distance=CalculateEuclidianDistance (list1,list2)
        print ("PDB file: ",pdb_input,"with distance of",distance,"Angstron", "between",atom1+str(residue_number1)+chain1,"and",atom2+str(residue_number2)+chain2)
        return distance
    except:
        print ("Failure:",pdb_input, residue_number1 , atom1 , chain1 , residue_number2 , atom2 , chain2)
        return 0

def CalculateEuclidianDistance (v1,v2):
    #def calculate_distance_XYZ (list1,list2): OLD
    if len(v1)==len(v2):
        distance = math.sqrt (sum( [(a - b) ** 2 for a, b in zip(v1, v2)] ) )
        #distance=numpy.sqrt( (list1[0]-list2[0])**2 + (list1[1]-list2[1])**2 + (list1[2]-list2[2])**2 )
        return distance
    else:
        print ("Different vectors lengths, correct input. Exiting function RJB_lib.CalculateEuclidianDistance.")
        exit()

def read_sequence_from_file_FASTA (seq_file):
    seq_file_read=open(seq_file,'r')
    seq_string_dirty=seq_file_read.read()
    while ">" in seq_string_dirty:
        seq_string_index1=seq_string_dirty.index(">")
        seq_string_index2=seq_string_dirty.index("\n",seq_string_index1)
        seq_string_dirty=seq_string_dirty[:seq_string_index1]+seq_string_dirty[seq_string_index2+2:]
    seq_string=seq_string_dirty.replace("\n","").replace("\t","").replace(" ","")
    return seq_string
    
##def given_file_containing_seq_for_each_line_writes_seqfile_with_star_for_SLIDERPOINTMUTATION ( file_with_seq , output_filename ): #Script written 25sep2015. seqfile should be given containing one seq per line already alligned
##    input_seq_file=open(file_with_seq,"r")
##    list_seq_byseq=input_seq_file.readlines()
##    output_file=open(output_filename,"w")
##    list_seq_by_position=[]
##    for column_number in range(len(file_with_seq[0]):
##        column_aminoacids=[]
##        for line_number in range(len(file_with_seq):
##            column_aminoacids.append(file_with_seq[file_with_seq][column_number])
##        list_seq_by_position.append

##def extract_table_filename_for_dictionary ( input_file ): #Script written in 5th October 2015 to extract info from OUT files 
##    input_fileee=open( input_file , 'r' )
##    input_file_list=input_fileee.readlines()
##    input_fileee.close()
##    labels=input_file_list[0].split()
##    dictio_all={}
##    for line in input_file_list[1:]:
##        dictio_var={}
##        line=line.split()
##        for i in range(len(labels)):
##            dictio_var[labels[i]]=line[i]
##        dictio_all[filename.split("/")[-1]]=dictio_var
##    return dictio_all


def remove_SLIDER_files_keeping_n_models ( number_models_keep , folder_to_be_removed):
    table_files=[]
    if folder_to_be_removed.endswith('/'):
        folder_to_be_removed=folder_to_be_removed[:-1]
    table_files.append(folder_to_be_removed+'_refmac_sorted')
    table_files.append(folder_to_be_removed+'_refmac_ROTAMER_sorted')
    table_files.append(folder_to_be_removed+'_BUSTER_sorted')
    table_files.append(folder_to_be_removed+'_buster_ROTAMER_sorted')
    os.system('rm '+folder_to_be_removed+'/model* > /dev/null')
    for table in table_files:
        try:
            table_file_read=open(table,"r")
            table_file_list=table_file_read.readlines()
            table_file_read.close()
            number_models_keep=int(number_models_keep)
            save_files=[]
            for i in range(len(table_file_list)):
                if i>number_models_keep and i!=0:
                    save_files.append(table_file_list[i].split()[0].split("_")[0]+"_")
            if 'refmac_sorted' in table:
                for item_remove in save_files:
                    os.system("rm "+folder_to_be_removed+"/refine/"+item_remove+"refmac1.*"+" "+folder_to_be_removed+"/refine/"+item_remove[:-1]+".* > /dev/null")
            elif '_refmac_ROTAMER_sorted' in table:
                for item_remove in save_files:
                    os.system("rm "+folder_to_be_removed+"/ROTAMER/"+item_remove+"*ROTAMERcoot.pdb"+" "+folder_to_be_removed+"/ROTAMER/"+item_remove+"*ROTAMERcoot_refmac1.* > /dev/null")
            elif '_BUSTER_sorted' in table:
                for item_remove in save_files:
                    os.system("rm "+folder_to_be_removed+"/refine/"+item_remove+"BUSTER.*"+" "+folder_to_be_removed+"/refine/"+item_remove[:-1]+".* > /dev/null")                    
            elif '_buster_ROTAMER_sorted' in table:
                for item_remove in save_files:
                    os.system("rm "+folder_to_be_removed+"/ROTAMER/"+item_remove+"*ROTAMERcoot.pdb"+" "+folder_to_be_removed+"/ROTAMER/"+item_remove+"*ROTAMERcoot_BUSTER.* > /dev/null")                    
        except:
            print ('\n\n', table ,'not found\n\n')
        


def run_MOLE_tunnel_search ( input_pdbfile , residue , chain , output_folder , path_mole2):
    xml_file=open(input_pdbfile[:-3]+'xml','w')
    xml_file.write('<?xml version="1.0" encoding="utf-8"?>\n<Tunnels>\n')
    xml_file.write('\t<Input>'+input_pdbfile+'</Input>\n')
    xml_file.write('\t<WorkingDirectory>./'+output_folder+'/</WorkingDirectory>\n')
    xml_file.write('\t<Params ProbeRadius="4" InteriorThreshold="1.25" SurfaceCoverRadius="10" OriginRadius="7" IgnoreHETAtoms="1" BottleneckRadius="1.4" />\n')
    xml_file.write('\t<Export Mesh="1" MeshGz="0" Cavities="1" MeshDensity="1.33" PyMol="1" PDB="1" />\n')
    xml_file.write('\t<Origin Auto="0">\n')
    #<Path>
    #xml_file.write('\t\t<Start>\n')
    xml_file.write('\t\t<Residue Chain="'+chain+'" SequenceNumber="'+str(residue)+'" />\n')
      #</Start>
      #<End>
       # <Residue Chain="B" SequenceNumber="111" />
       # <Residue Chain="A" SequenceNumber="10" />
      #</End>
   # </Path>
    xml_file.write('\t</Origin>\n')
    xml_file.write('</Tunnels>')
    xml_file.close()
    mono_path=which_program ("mono")
    os.system(mono_path+' '+path_mole2+' '+input_pdbfile[:-3]+'xml')



def extract_info_cavities_Mole2 (cavities_inputfile,number_of_cavity):
    cavit_file=open(cavities_inputfile,'r')
    cavit_file_list=cavit_file.readlines()
    cavit_file.close()
    hpo_hpa_pol=[]
    for i,line in enumerate(cavit_file_list):
        if line.startswith('  <Cavity Type="Cavity" ') and line.endswith('" Id="'+str(number_of_cavity)+'">\n'):
            print ('Desired cavity found:')
            print (line)
            volume=float(line.split()[2][8:-1])
            ii=i+1
            l=cavit_file_list[ii]
            while not l.startswith('  </Cavity>'):#l.startswith('  <Cavity Type="Cavity" ') and not l.startswith('</Cavities>'):
                #print l
                ii+=1
                l=cavit_file_list[ii+1]
                if l.replace(' ','').startswith('<PropertiesCharge='):
                    kk=['Charge="','NumPositives="','NumNegatives="','Hydrophobicity="','Hydropathy="','Polarity="']
                    print ('cavity function extracting values of line:')
                    print (l)
                    lista=[]
                    for n,k in enumerate(kk):
                        v=l[l.index(k)+len(k): l[l.index(k)+len(k):].index('"')+l.index(k)+len(k) ]
                        lista.append(v)
                    hpo_hpa_pol.append(lista)
    try:
        print ('hpo_hpa_pol',hpo_hpa_pol)
        return volume,hpo_hpa_pol
    except:
        return 0

def extract_chemical_info_tunnel_Mole2 (xmlfile,tun_id):
    ff=open(xmlfile)
    f=ff.readlines()
    ff.close()
    hpo_hpa_pol=[]
    for i,l in enumerate(f):
        if l.startswith('  <Tunnel Id="'+str(tun_id)):
            print ('Desired tunnel found:')
            print (l)
            kk=['Charge="','NumPositives="','NumNegatives="','Hydrophobicity="','Hydropathy="','Polarity="']
            ll=f[i+1]
            print ('tunnel function extracting values of line:')
            print (ll)
            k='Cavity="'
            ncavity=l[l.index(k)+len(k): l[l.index(k)+len(k):].index('"')+l.index(k)+len(k) ]
            for n,k in enumerate(kk):
                v=ll[ll.index(k)+len(k): ll[ll.index(k)+len(k):].index('"')+ll.index(k)+len(k) ]
                hpo_hpa_pol.append(v)
    tunpdb=xmlfile[:-5]+'_'+str(tun_id)+'.pdb'
    tunVol=extract_info_tunnel_Mole2 (tunnel_input_file=tunpdb)
    print ('hpo_hpa_pol',hpo_hpa_pol)
    print ('Tunnel volume',tun_id,'calculated from:',tunpdb,tunVol)
    return hpo_hpa_pol,ncavity,tunVol


def extract_chemical_info_nparts_tunnel_Mole2 (xmlfile,tun_id,nparts):
    ff=open(xmlfile)
    f=ff.readlines()
    ff.close()
    overallvalues=[]
    lvalues=[]
    for i,l in enumerate(f):
        if l.startswith('  <Tunnel Id="'+str(tun_id)):
            print ('Desired tunnel found:')
            print (l)
            kk=['Charge="','NumPositives="','NumNegatives="','Hydrophobicity="','Hydropathy="','Polarity="']
            ll=f[i+1]
            print ('tunnel function extracting values of line:')
            print (ll)
            k='Cavity="'
            ncavity=l[l.index(k)+len(k): l[l.index(k)+len(k):].index('"')+l.index(k)+len(k) ]
            for n,k in enumerate(kk):
                v=ll[ll.index(k)+len(k): ll[ll.index(k)+len(k):].index('"')+ll.index(k)+len(k) ]
                overallvalues.append(v)
            ii=i
            lll=f[ii]
            while not lll.startswith('    </Layers>'):#l.startswith('  <Cavity Type="Cavity" ') and not l.startswith('</Cavities>'):
                #print lll
                if lll.startswith('      <Layer MinRadius'):
                    print ('tunnel function extracting values of line:')
                    print (lll)
                    k='StartDistance="'
                    v=lll[lll.index(k)+len(k): lll[lll.index(k)+len(k):].index('"')+lll.index(k)+len(k) ]
                    k='EndDistance="'
                    enddist=lll[lll.index(k)+len(k): lll[lll.index(k)+len(k):].index('"')+lll.index(k)+len(k) ]
                    lista=[]
                    lista.append(v)
                    ii+=3
                    lll=f[ii]
                    print (lll )
                    kk=['Charge="','NumPositives="','NumNegatives="','Hydrophobicity="','Hydropathy="','Polarity="']
                    for n,k in enumerate(kk):
                        #print k
                        v=lll[lll.index(k)+len(k): lll[lll.index(k)+len(k):].index('"')+lll.index(k)+len(k) ]
                        lista.append(v)
                    lista.append(enddist)
                    lvalues.append(lista)
                    #print lvalues
                ii+=1
                lll=f[ii]
            
            lengthtunnel=lvalues[-1][-1]
            div=float(lengthtunnel)/nparts
            ldiv=[]
            svalues=[]
            print ('Calculating values')
            print ('Size of parts of tunnels are:',div)
            for i in range(1,nparts+1):
                ldiv.append(div*i)
            for FlEndT in ldiv:
                #print FlEndT
                values=[]
                #values will contain:kk=['Charge="','NumPositives="','NumNegatives="','Hydrophobicity="','Hydropathy="','Polarity="']
                for i in range(len(kk)):
                    nsum=0
                    var=[]
                    for lit in lvalues:
                        startv=float(lit[0])
                        endv=float(lit[-1])
                        if endv>FlEndT:
                            endv=FlEndT
                        diff=endv-startv
                        if startv<FlEndT and endv<=FlEndT:
                            #print startv
                            #print endv
                            #print 'list inside float',lit
                            var.append(float(lit[i+1])*diff)
                            nsum+=diff
                    val=sum(var)/nsum
                    val='%.1f'%(val)
                    values.append(val)
                svalues.append(values)
    tunpdb=xmlfile[:-5]+'_'+str(tun_id)+'.pdb'
    tunVol=extract_info_tunnel_Mole2 (tunnel_input_file=tunpdb)
    #print 'overallvalues',overallvalues
    #print 'Tunnel volume',tun_id,'calculated from:',tunpdb,tunVol
    print ('calculated values are:',svalues)
    return overallvalues,ncavity,tunVol,svalues



def extract_info_tunnel_Mole2 (tunnel_input_file):
    try:
        tunnel_file=open(tunnel_input_file,'r')
        tunnel_file_list=tunnel_file.readlines()
        tunnel_file.close()
        dist=[]
        radius=[]
        for line in tunnel_file_list:
            if line.startswith('HETATM'):
                dist.append(float(line[54:60]))
                radius.append(float(line[60:66]))
        volume=numpy.mean(radius)*numpy.mean(radius)*math.pi*dist[-1]
        return volume
    except:
        return 0


def extract_average_final_tunnel_Mole2 (tunnel_input_file):
    try:
        tunnel_file=open(tunnel_input_file,'r')
        tunnel_file_list=tunnel_file.readlines()
        tunnel_file.close()
        dist=[]
        radius=[]
        for line in tunnel_file_list:
            if line.startswith('HETATM'):
                dist.append(float(line[54:60]))
                radius.append(float(line[60:66]))
        volume=numpy.mean(radius)*numpy.mean(radius)*math.pi*dist[-1]
        avg_dist=numpy.mean(dist)
        avg_radius=numpy.mean(radius)
        return [avg_dist,dist[-1]],[avg_radius,radius[-1]]
    except:
        return 0



def extract_coordinates (pdb_path, residue , chain , atom=False , exitt=True):
    if not atom:
        atom='CA'
    atom+=' '*(3-len(atom))
    pdb_fileee=open(pdb_path,'r')
    pdb_list=pdb_fileee.readlines()
    pdb_fileee.close()
    check=False
    for line in pdb_list :
        if line.startswith("ATOM"):
            if line[21]==chain and int(line[22:26])==int(residue) and line[13:16]==atom:
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                check=True
    if check: return [x,y,z]
    else:
        atom=atom.replace(' ','')
        print ("Failure looking for pdb file",pdb_path,'residue', residue , 'chain' ,chain, 'atom','"'+atom+'"')
        if atom=='CB':
            print ("Obtaining theoretical Cb position based on N/C/CA positions using BioPython.\n")
            cb=generateGlyCbBioPython(pdbinput=pdb_path, resn=int(residue), chain=chain)
            return cb
        if exitt:
            print ('Exiting.')
            exit()
        else: return False

def generateGlyCbBioPython (pdbinput , resn, chain ):
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdbinput[:-4], pdbinput)
    residuee=struct[0][chain][resn]
    # get atom coordinates as vectors
    n = residuee['N'].get_vector()
    c = residuee['C'].get_vector()
    ca =residuee['CA'].get_vector()
    # center at origin
    n = n - ca
    c = c - ca
    # find rotation matrix that rotates n -120 degrees along the ca-c vector
    rot = Bio.PDB.rotaxis(-pi * 120.0 / 180.0, c)
    # apply rotation to ca-n vector
    cb_at_origin = n.left_multiply(rot)
    # put on top of ca atom
    cb = cb_at_origin + ca
    cbb=[cb[0],cb[1],cb[2]]
    return cbb

def normal_vector_from_two_3D_vectors (v1,v2): #v1 and v2 should be a list 3 [x,y,z]
    v1x=v1[0]
    v1y=v1[1]
    v1z=v1[2]
    v2x=v2[0]
    v2y=v2[1]
    v2z=v2[2]
    nx=v1y*v2z-v1z*v2y
    ny=v1x*v2z-v1z*v2x
    nz=v1x*v2y-v1y*v2x
    return [nx,ny,nz]

def write_table_from_list_dictio (list_item,dic_all,output_filename): #not tested, generated 4/4/16 to write superimposition tables monomer B with monB
    output_file=open(output_filename,'w')
    for item in list_item[1:]:
        output_file.write('\t'+item[:4])
    output_file.write('\n')
    for n1 in range(len(list_item)):
        output_file.write(list_item[n1][:4]+'\t'*n1)
        for item2 in list_item[n1+1:]:
            angle='%.1f'%(dic_all[list_item[n1]][item2])
            output_file.write('\t'+angle)
        output_file.write('\n')
    output_file.close()

    
def write_table_inverted_from_list_dictio (list_item,dic_all,output_filename):  #generated 4/4/16 to write superimposition tables (inverted)
    output_file=open(output_filename,'w')
    columns=list_item[::-1][:-1]
    lines=list_item[:-1]
    for item in columns: #[::-1] -> invert list
        output_file.write('\t'+item)
    output_file.write('\n')
    for line_item in lines:
        output_file.write(line_item)
        for column_item in columns:
            try:
                angle='%.1f'%(dic_all[line_item][column_item])
                output_file.write('\t'+angle)
                del dic_all[line_item][column_item]
            except:
                try:
                    angle='%.1f'%(dic_all[column_item][line_item])
                    output_file.write('\t'+angle)
                    del dic_all[column_item][line_item]
                except:
                    pass
        output_file.write('\n')
    output_file.close()



def superimpose_dif_seq_chain_rot_angle (pdb_input1_fixed,residue1_1,residue1_2,chain1,pdb_input2_moving,residue2_1,residue2_2, chain2):#,pdb_output):
    superimpose_instructions=open("superimpose_instructions.log","w")
    superimpose_instructions.write("fit res CA "+str(residue1_1)+" to "+str(residue1_2)+' chain '+chain1.upper()+'\n')
    superimpose_instructions.write("match "+str(residue2_1)+" to "+str(residue2_2)+' chain '+chain2.upper()+'\n')
    #superimpose_instructions.write("output -\n")
    #superimpose_instructions.write("    xyz\n")
    superimpose_instructions.write("end\n")
    superimpose_instructions.close()
    superimpose_instructions=open("superimpose_instructions.log","r")
    superimpose_job = subprocess.Popen([ 'lsqkab','XYZIN1',pdb_input1_fixed,'XYZIN2',pdb_input2_moving,], stdin=superimpose_instructions, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)#'XYZOUT',pdb_output], stdin=superimpose_instructions, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = superimpose_job.communicate()
    #print out
    try:
        print (pdb_input1_fixed[:4],pdb_input2_moving[:4],"lsqkab angle:",out[out.index("  SPHERICAL POLARS OMEGA PHI CHI")+60:out.index("  SPHERICAL POLARS OMEGA PHI CHI")+69])
        return float(out[out.index("  SPHERICAL POLARS OMEGA PHI CHI")+60:out.index("  SPHERICAL POLARS OMEGA PHI CHI")+69])
    except:
        print (out)
        return "error"
        exit()

def superimpose_lsqkab_dif_seq_chain_pdb_atomtype_output (pdb_input1_moving,chain1,residue1_1,residue1_2, pdb_input2_reference,chain2,residue2_1,residue2_2,type_superposition,output_name,output_options='xyz',show_rmsd=False):
    #type_superposition should be CA | MAIN | SIDE | ALL
    #output_options should be xyz (coordinates), log , xyzlog , deltas (distances), xyzdeltas order is not essential
    if output_name.endswith('.pdb') or output_name.endswith('.log') or output_name.endswith('.dat'):
        output_name=output_name[:-4]
    if 'xyz' in output_options:
        xyz=True
    else:
        xyz=False
    if 'log' in output_options:
        log=True
    else:
        log=False
    if 'deltas' in output_options:
        deltas=True
    else:
        deltas=False
    lch1=[]
    chain1=chain1.replace(',','').replace('.','').replace(';','')
    for i in chain1:
        lch1.append(i)
    lch2=[]
    chain2=chain2.replace(',','').replace('.','').replace(';','')
    for i in chain2:
        lch2.append(i)
    
    superimpose_instructions=open("superimpose_instructions.log","w")
    for i,s in enumerate(lch1):
        superimpose_instructions.write("fit res "+type_superposition+" "+str(residue1_1)+" to "+str(residue1_2)+' chain '+lch1[i].upper()+'\n')
        superimpose_instructions.write("match "+str(residue2_1)+" to "+str(residue2_2)+' chain '+lch2[i].upper()+'\n')
    #superimpose_instructions.write("output -\n")
    if xyz:
        #superimpose_instructions.write("    xyz\n")
        superimpose_instructions.write("output xyz\n")
    if deltas:
        #superimpose_instructions.write("    deltas\n")
        superimpose_instructions.write("output deltas\n")

    superimpose_instructions.write("end\n")
    superimpose_instructions.close()
    superimpose_instructions=open("superimpose_instructions.log","r")
    if xyz and deltas:
        superimpose_job = subprocess.Popen([ 'lsqkab','XYZINM',pdb_input1_moving,'XYZINF',pdb_input2_reference,'XYZOUT',output_name+'.pdb','DELTAS',output_name+'.dat' ], stdin=superimpose_instructions, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)        # XYZINM as in match XYZINF as in fit, DELTAS (A list of ALL differences between atom pairs is written to a file assigned extracted from http://www.ccp4.ac.uk/html/lsqkab.html#output)
##        if not os.path.isfile( output_name+'.pdb' ) and not os.path.isfile( output_name+'.dat' ):
##            print 'error in function RJB_lib.superimpose_dif_seq_chain_output_pdb'
##            print 'PDB moving:',pdb_input1_moving,'\nPDB reference:',pdb_input2_reference,'\nPDB output:',output_name+'.pdb','\nDELTAS:',output_name+'.dat'
##            exit()
    elif xyz:
        superimpose_job = subprocess.Popen([ 'lsqkab','XYZINM',pdb_input1_moving,'XYZINF',pdb_input2_reference,'XYZOUT',output_name+'.pdb'], stdin=superimpose_instructions, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)        # XYZINM as in match XYZINF as in fit
##        if not os.path.isfile( output_name+'.pdb' ) :
##            print 'error in function RJB_lib.superimpose_dif_seq_chain_output_pdb'
##            print 'PDB moving:',pdb_input1_moving,'\nPDB reference:',pdb_input2_reference,'\nPDB output:',output_name+'.pdb'
##            exit()
    elif deltas:
        superimpose_job = subprocess.Popen([ 'lsqkab','XYZINM',pdb_input1_moving,'XYZINF',pdb_input2_reference,'DELTAS',output_name+'.dat'], stdin=superimpose_instructions, stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)        # XYZINM as in match XYZINF as in fit, DELTAS (A list of ALL differences between atom pairs is written to a file assigned extracted from http://www.ccp4.ac.uk/html/lsqkab.html#output)
##        if not os.path.isfile( output_name+'.dat' ):
##            print 'error in function RJB_lib.superimpose_dif_seq_chain_output_pdb'
##            print 'PDB moving:',pdb_input1_moving,'\nPDB reference:',pdb_input2_reference,'\nDELTAS:',output_name+'.dat'
##            exit()
    else:
        print ('ERROR in RJB_lib.superimpose_dif_seq_chain_output_pdb function verifying output_options',output_options)
    out, err = superimpose_job.communicate()
    if log:
        log_file=open(output_name+'.log','w')
        log_file.write(out)
        log_file.close()
    try:
        if show_rmsd:
            print ('PDB flexible',pdb_input1_moving,'chain',chain1,'residues',residue1_1,'to',residue1_2,'being superimposed (lsqkab) into PDB reference',pdb_input2_reference,'chain',chain2,'residues',residue2_1,'to',residue2_2,'and saved with filename:')
            if xyz:
                print (output_name+'.pdb')
            if deltas:
                print (output_name+'.dat')
            if log:
                print (output_name+'.log')
            print (out[out.index("          RMS     XYZ DISPLACEMENT = ")+36:out.index("          RMS     XYZ DISPLACEMENT = ")+44])
        return out[out.index("          RMS     XYZ DISPLACEMENT = ")+36:out.index("          RMS     XYZ DISPLACEMENT = ")+44]
    except:
        print (out)
        print ("error in function RJB_lib.superimpose_dif_seq_chain_output_pdb")
        return "error in function RJB_lib.superimpose_dif_seq_chain_output_pdb"
        exit()


def write_pdb_Euler_angles (input_table, output_pdbfile):
    table=open(input_table,'r')
    table_lines=table.readlines()[1:]
    table.close()
    pdbfile=open(output_pdbfile,'w')
    originfile=open('origin.pdb','w')
    originfile.write('HETATM    0      ORI0    0       0.000   0.000   0.000\n')
    originfile.write('HETATM    1      ORI1    1     180.000   0.000   0.000\n')
    originfile.write('HETATM   -1      ORI1   -1    -180.000   0.000   0.000\n')
    originfile.write('HETATM    2      ORI2    2       0.000 180.000   0.000\n')
    originfile.write('HETATM   -2      ORI2   -2       0.000-180.000   0.000\n')
    originfile.write('HETATM    3      ORI3    3       0.000   0.000 180.000\n')
    originfile.write('HETATM   -3      ORI3   -3       0.000   0.000-180.000\n')
    originfile.close()
    i=1
    for line in table_lines:
        line=line.split()
        x='%.3f'%(float(line[1]))
        y='%.3f'%(float(line[2]))
        z='%.3f'%(float(line[3]))
        pdbfile.write('HETATM'+' '*(5-len(str(i))) +str(i)+' '*(9-len(line[0][:4]))+line[0][:4]+' '*(6-len(str(i))) +str(i)+' '*4)
        pdbfile.write(' '*(8-len(x))+x )
        pdbfile.write(' '*(8-len(y))+y )
        pdbfile.write(' '*(8-len(z))+z )
        pdbfile.write('\n')
        i+=1
    pdbfile.close

def change_pdb_chain_A_with_B (pdb_input, pdb_output): # used for angle calculation script created 2016/04/28
    pdb_input_read=open(pdb_input,'r')
    list_pdb_input=pdb_input_read.readlines()
    pdb_input_read.close()
    file_out=open(pdb_output,'w')
    for line in list_pdb_input:
        if line.startswith('ATOM'):
            if line[21]=='A':
                file_out.write(line[:21]+'B'+line[22:])
            elif line[21]=='B':
                file_out.write(line[:21]+'A'+line[22:])
            else:
                print ('Script RJB_lib.change_pdb_chain_A_with_B was unable to find chain letter A or B, check following PDB and line if that is the case, check line if that is the case')
                print (pdb_input,line)
        elif line.startswith('REMARK'): file_out.write(line)
    file_out.close()


def renumber_protein_leaving_HETATM_unchanged_continuous_respecting_chain (pdb_input,pdb_output): # used for angle calculation script created 2016/04/28
    pdb_fileee = open(pdb_input,"r")
    pdb_file_list = pdb_fileee.readlines()
    pdb_fileee.close()

    out_file = open(pdb_output,"w")

    new_res_numbering=1

    for line in range (len(pdb_file_list)):
        if not pdb_file_list[line].startswith("ATOM") and not pdb_file_list[line].startswith("ANISOU") and not pdb_file_list[line].startswith("TER"):
            out_file.write(pdb_file_list[line])
        elif pdb_file_list[line].startswith("ATOM"):
            #print pdb_file_list[line]
            residue_number=int(pdb_file_list[line][22:26])
            chain_ID=pdb_file_list[line][21]
            #print residue_number
            #chain_ID=pdb_file_list[line][21]
            #print chain_ID

    #       if pdb_file_list[line-1].startswith('ATOM'):
    #            n=1
    #       else:
    #               if pdb_file_list[line-2].startswith('ATOM'):
    #                       n=2
    #               else:
    #                       if pdb_file_list[line-3].startswith('ATOM'):
    #                               n=3

            try:
                previous_res
                chain_ID_previous
                #residue_number_previous=previous_line[22:26]
                #chain_ID_previous=pdb_file_list[line-n][21]
                #if chain_ID_previous!=chain_ID:
                ##res_numb=int(residue_number)
            except:
                previous_res=residue_number
                chain_ID_previous=chain_ID
                #previous_res=int(pdb_file_list[22:26])
                #residue_number_previous=residue_number
                #chain_ID_previous=chain_ID
                ##res_numb=int(residue_number)

            #if previous_res==residue_number :
                #out_file.write(pdb_file_list[line])
            if previous_res!=residue_number and chain_ID_previous==chain_ID:
                new_res_numbering=new_res_numbering+1
            elif previous_res!=residue_number and not chain_ID_previous==chain_ID:
                new_res_numbering=1

            size_res_numb=len(str (new_res_numbering) )
            out_file.write(pdb_file_list[line][:22])
            out_file.write( (26-22-size_res_numb)*" "  )
            out_file.write(str(new_res_numbering))
            out_file.write(pdb_file_list[line][26:])
            previous_res=residue_number
            chain_ID_previous=chain_ID
    out_file.close()

def convert_file_gro_pdb (gro_input, pdb_output): #PROBLEM!!!
    gro_input_read=open(gro_input)
    list_gro_input=gro_input_read.readlines()
    gro_input_read.close()
    pdb=open(pdb_output,'w')
    for line in list_gro_input:
        if len(line.split())==6:
            try: #exclude different lines
                pdb.write('ATOM  ') #beginning
                pdb.write(line[15:20]+' ') #atom number
                line_split=line.split() #some writting characters are easier separating in list of strings
                if len(line_split[1])==4: #for atom with 4 characters
                    pdb.write(line_split[1]+' ') #atom writting
                else:
                    pdb.write(' '+line_split[1]+' '*(4-len(line_split[1])))
                pdb.write(line[5:8]+' ') #number of residue
                protein_lenth=121 #how many residues do your protein has? I.e., what is the number of last residue in chain A? Rest goes multiplicating, i.e., equal number of residues in each chain is assumed.
                if int(line[:5])<(protein_lenth+1): #write chain A
                    pdb.write('A')
                elif protein_lenth<int(line[:5])<(protein_lenth*2+1): #write chain B
                    pdb.write('B')
                else:
                    pdb.write('C') #rest goes to chain C
                pdb.write(line[1:5]+' '*4)
                pdb.write(line[20:28]) #X coordinate
                pdb.write(line[28:36]) #Y coordinate
                pdb.write(line[36:44]) #Z coordinate
                pdb.write('  1.00 20.00           ')
                pdb.write(line_split[1][0]+' \n')
            except: #show not used lines
                print (gro_input,'line:\n',line,'\nnot written')
        else: #show not used lines
            print (gro_input,'line:\n',line,'\nnot written')


def generate_PDB_with_CarbonCentroid7_from_PDB (input_pdb,output_pdb): #generated April 28, 2016 for automatic angle calculation
    fileee=open(input_pdb)
    list_pdb_file=fileee.readlines()
    fileee.close()
    file_out=open(output_pdb,'w')
    dic_all_CA={}
    list_chains=[]
    for line in list_pdb_file:
        if line.startswith('ATOM') and line[13:15]=='CA':
            x=float(line[31:38])
            y=float(line[39:45])
            z=float(line[47:54])
            chain=line[21]
            res_numb=int(line[22:26])
            res_type=line[17:20]
            list_CA=[x,y,z,res_numb,res_type]
            try:
                dic_all_CA[chain].append(list_CA)
            except:
                dic_all_CA[chain]=[]
                dic_all_CA[chain].append(list_CA)
    for chain in dic_all_CA:
        for i in range(3,len(dic_all_CA[chain])-3):
    ##                list_CA1=dic_all_CA[chain][i-1]
    ##                x1=list_CA1[0]
    ##                y1=list_CA1[1]
    ##                z1=list_CA1[2]
    ##                list_CA2=dic_all_CA[chain][i]
    ##                x2=list_CA2[0]
    ##                y2=list_CA2[1]
    ##                z2=list_CA2[2]
    ##                res_numb=list_CA2[3]
    ##                res_numb=' '*(4-len(str(res_numb)))+str(res_numb)
    ##                res_type=list_CA2[4]
    ##                list_CA3=dic_all_CA[chain][i+1]
    ##                x3=list_CA3[0]
    ##                y3=list_CA3[1]
    ##                z3=list_CA3[2]
    ##                Cent_x=(x1+x2+x3)/3
    ##                Cent_y=(y1+y2+y3)/3
    ##                Cent_z=(z1+z2+z3)/3
            x=0
            y=0
            z=0
            n_sum=0
            for ii in range(-3,4):
                list_CA=dic_all_CA[chain][i-ii]
                x+=list_CA[0]
                y+=list_CA[1]
                z+=list_CA[2]
                if ii==0:
                    res_numb=list_CA[3]
                    res_numb=' '*(4-len(str(res_numb)))+str(res_numb)
                    res_type=list_CA[4]
                n_sum+=1
            Cent_x=x/n_sum
            Cent_y=y/n_sum
            Cent_z=z/n_sum

            Cent_x='%.3f'%(Cent_x)
            Cent_x=' '*(7-len(Cent_x))+Cent_x+' '
            Cent_y='%.3f'%(Cent_y)
            Cent_y=' '*(7-len(Cent_y))+Cent_y+' '
            Cent_z='%.3f'%(Cent_z)
            Cent_z=' '*(7-len(Cent_z))+Cent_z+' '
            file_out.write('ATOM   '+res_numb+'  CA  '+res_type+' '+chain+res_numb+'     '+Cent_x+Cent_y+Cent_z+' 1.00 20.00           C  \n')
    file_out.close()

def generate_CarbonCentroid7_in_Morphing_PDB (input_pdb,output_pdb,list_res): #generated April 28, 2016 for automatic angle calculation #modified Nov1, 2016. For add axes in morph file
    fileee=open(input_pdb)
    list_pdb_file=fileee.readlines()
    fileee.close()
    file_out=open(output_pdb,'w')
    dic_all_CA={}
    list_res2=[]
    for res in list_res:
        for i in range(res-3,res+4):
            if not i in list_res2:
                list_res2.append(i)
    for ind,line in enumerate(list_pdb_file):
        if line.startswith('ATOM') and line[13:15]=='CA' and int(line[22:26]) in list_res2:
            x=float(line[31:38])
            y=float(line[39:45])
            z=float(line[47:54])
            chain=line[21]
            res_numb=int(line[22:26])
            res_type=line[17:20]
            list_CA=[x,y,z,res_numb,res_type]
            try:
                dic_all_CA[chain][res_numb]=list_CA
            except:
                dic_all_CA[chain]={}
                dic_all_CA[chain][res_numb]=list_CA
        elif line.startswith('ENDMDL'):
            last_atom_numb=int(list_pdb_file[ind-1][5:11])
            c=0
            for chain in dic_all_CA:
                for res in list_res:
                    for i in range(3,len(dic_all_CA[chain])-3):
                        x=0
                        y=0
                        z=0
                        n_sum=0
                        for i in range(-3,4):
                            list_CA=dic_all_CA[chain][res+i]
                            x+=list_CA[0]
                            y+=list_CA[1]
                            z+=list_CA[2]
                            if i==0:
                                res_numb=200+res
                                res_numb=' '*(4-len(str(res_numb)))+str(res_numb)
                                res_type='HOH'
                            n_sum+=1
                    Cent_x=x/n_sum
                    Cent_y=y/n_sum
                    Cent_z=z/n_sum

                    Cent_x='%.3f'%(Cent_x)
                    Cent_x=' '*(7-len(Cent_x))+Cent_x+' '
                    Cent_y='%.3f'%(Cent_y)
                    Cent_y=' '*(7-len(Cent_y))+Cent_y+' '
                    Cent_z='%.3f'%(Cent_z)
                    Cent_z=' '*(7-len(Cent_z))+Cent_z+' '
                    atomnumb=last_atom_numb+c
                    atomnumb=' '*(5-len(str(atomnumb)))+str(atomnumb)
                    file_out.write('ATOM  '+atomnumb+'  O   '+res_type+' '+chain+res_numb+'     '+Cent_x+Cent_y+Cent_z+' 1.00 20.00           O  \n')
                    c+=1
            dic_all_CA={}
        file_out.write(line)
    file_out.close()

def lsqkab_delta_read_into_dict (input_filedat): #script reads file.dat of lsqkab (deltas output) and return distances into a dictio like: dic[res#][atomtype]=dist# > res# in int ; atomtype string (like CA) ; dist# in float  , extracted from "~/DATA/MD/debian/Rafael/
    input_filedat_file=open(input_filedat,"r")
    input_filedat_list=input_filedat_file.readlines()
    dic_dat={}
    for line in input_filedat_list:
        #started writting code assuming output was based on characters number, but probably other ccp4 version was different, thus wrote something with split
        #resnumb=int(line[11:15])
        #atomtype=line[15:17]
        #atomtype=atomtype.replace(" ","")
        #distance=float(line[:10])
        line_split=line.split()
        distance=float(line_split[0])
        for i in range(len(line_split[1])):
            if line_split[1][i].isdigit():
                last_number=i+1
        resnumb=int(line_split[1][:last_number])
        atomtype=line_split[1][last_number:]
#        try:
#	    resnumb=int(line_split[1][:-1]
#            atomtype=line_split[1][-1:]
#	except:
#	    try:
#	      resnumb=int(line_split[1][:-2]
#            except:
#	      resnumb=int(line_split[1][:-3]
        try:
            dic_dat[resnumb][atomtype]=distance
        except:
            dic_var={}
            dic_var[atomtype]=distance
            dic_dat[resnumb]=dic_var
    return dic_dat

#
def convert_list_of_dictios_of_two_keys_and_stringfloat_into_one_dictio_of_two_keys_of_lists ( list_dictios ): #used for automatic-angle-calculation to manage various pair of distances
    dic_all={}
    for dic in list_dictios:
        for key1 in dic:
                for key2 in dic[key1]:
                    try:
                        dic_all[key1][key2].append(dic[key1][key2])
                    except:
                        try:
                            dic_all[key1][key2]=[]
                            dic_all[key1][key2].append(dic[key1][key2])
                        except:
                            dic_var={}
                            dic_all[key1]=dic_var
                            dic_all[key1][key2]=[]
                            dic_all[key1][key2].append(dic[key1][key2])
    return dic_all

## UK: 20th May 2016
def read_secstr_pred_psipred (input_file):
    file_read=open(input_file,"r")
    file_list=file_read.readlines()
    file_read.close()
##    confidence=file_list[2].split()[1]
##    sectr=file_list[3].split()[1]
##    sequence=file_list[4].split()[1]
    #considering that 
    for line_index in range(len(file_list)):
        line=file_list[line_index]
        if file_list[line_index].startswith('Conf:'):
            try:
                confidence+=file_list[line_index].split()[1]
                sec_str_pred+=file_list[line_index+1].split()[1]
                sequence+=file_list[line_index+2].split()[1]
            except:
                confidence=file_list[line_index].split()[1]
                sec_str_pred=file_list[line_index+1].split()[1]
                sequence=file_list[line_index+2].split()[1]
    list_secstr_pred=[]
    aa_verify=sec_str_pred[0]
    string_var=''
    for letter in sec_str_pred:
        if aa_verify==letter:
            string_var+=letter
        else:
            list_secstr_pred.append(string_var)
            string_var=letter
            aa_verify=letter 
    if len(string_var)==1:
        list_secstr_pred.append(string_var)
    #return sequence,sec_str_pred,confidence,list_secstr_pred
    dic_secstr_pred_psipred={}
    dic_secstr_pred_psipred['sequence']=sequence
    dic_secstr_pred_psipred['sec_str_pred']=sec_str_pred
    dic_secstr_pred_psipred['confidence']=confidence
    dic_secstr_pred_psipred['list_secstr_pred']=list_secstr_pred
    return dic_secstr_pred_psipred
            

def verify_and_correct_aa_from_shelxepdb_pdbRJBdic (pdbRJBdic): #script to exclude partial aminoacids from shelxe output pdb
    correction=False
    list_exclude=[]
    for chain in pdbRJBdic:
        for res_numb in pdbRJBdic[chain] :
            if len(pdbRJBdic[chain][res_numb])<4:
                correction=True
                list_exclude.append((chain,res_numb))
    for excl in list_exclude:
        del pdbRJBdic[excl[0]][excl[1]]
    return correction,pdbRJBdic

    
##    file_input=open(pdb_input,"r")
##    pdb_file_list=file_read.readlines()
##    file_input=close()
##    file_output=open(pdb_output,'w')
##    pdbRJBdic={} #pdbRJBdic will have key chain and then residue_number and then list of atoms
##    for index_line in range(len(pdb_file_list)):
##        if pdb_file_list[index_line].startswith("ATOM") or pdb_file_list[index_line].startswith("ANISOU"):
##            residue_number=int(pdb_file_list[index_line][22:26])
##            chain_ID=pdb_file_list[index_line][21]
##            atom_type=pdb_file_list[index_line][13:15]
##            try:
##                pdbRJBdic[chain_ID][residue_number].append(atom_type)
##            except:
##                try:
##                    pdbRJBdic[chain_ID][residue_number]=[]
##                    pdbRJBdic[chain_ID][residue_number].append(atom_type)
##                except:
##                    pdbRJBdic[chain_ID]={}
##                    pdbRJBdic[chain_ID][residue_number]=[]
##                    pdbRJBdic[chain_ID][residue_number].append(atom_type)
##    #remove partial residues
##    for chain in pdbRJBdic:
##        for residue in pdbRJBdic[chain]:
##            if not "CA" and not "C" and not "N" and not "O" in pdbRJBdic[chain][residue_number]:
##                del pdbRJBdic[chain][residue_number]
##    #write complete main chain residues
##    for index_line in range(len(pdb_file_list)):
##        if not pdb_file_list[index_line].startswith("ATOM") or pdb_file_list[index_line].startswith("ANISOU"):
##            file_output.write(pdb_file_list[index_line])
##        else:
##            residue_number=int(pdb_file_list[index_line][22:26])
##            chain_ID=pdb_file_list[index_line][21]
##            try:
##                pdbRJBdic[chain][residue_number]
##                file_output.write(pdb_file_list[index_line])
##            except:
##                pass
    file_write.close()
    
    
def convert_pdb_into_pdbRJBdic (pdb_input): #pdbRJBdic[chain][res_number][atom_type/res_type] : res_type=string atom_type keys: [X],[Y],[Z],[occupancy],[Bfactor] floats
    file_input=open(pdb_input,"r")
    pdb_file_list=file_input.readlines()
    file_input.close()
    pdbRJBdic={}
    pdb_header_titles=('HEADER', 'TITLE ', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL  ', 'REMARK', 'DBREF ', 'SEQRES', 'FORMUL', 'HELIX ', 'SHEET ', 'SSBOND', 'CISPEP', 'CRYST1', 'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3')
    pdb_header_list_lines=[]
    for line in pdb_file_list:
        if line.startswith("ATOM"):
            chain_ID=line[21]
            try:
                pdbRJBdic[chain_ID]
            except:
                pdbRJBdic[chain_ID]={}
            residue_number=int(line[22:26])
            res_type=line[17:20]
            try:
                pdbRJBdic[chain_ID][residue_number]
            except:
                pdbRJBdic[chain_ID][residue_number]={}
                pdbRJBdic[chain_ID][residue_number]['res_type']=res_type
            atom_type=line[13:15]
            dic_var={}
            atom_numb=int(line[5:11])
            x=float(line[31:38])
            y=float(line[39:46])
            z=float(line[47:54])
            Bfactor=float(line[61:66])
            occupancy=float(line[56:60])
            dic_var['atom_numb']=atom_numb
            dic_var['X']=x
            dic_var['Y']=y
            dic_var['Z']=z
            dic_var['Bfactor']=Bfactor
            dic_var['occupancy']=occupancy
            pdbRJBdic[chain_ID][residue_number][atom_type]=dic_var
        elif line[:6] in pdb_header_titles:
            pdb_header_list_lines.append(line)
        else:
            print ('\nRJB_lib.pdbRJBdic function >>> Following line was excluded from RJB_lib.pdbRJBdic function')
            print ('RJB_lib.pdbRJBdic function >>> '+line+'\n')
    return pdbRJBdic,pdb_header_list_lines

def write_pdb_from_pdbRJBdic (pdbRJBdic,pdb_header_list_lines,output_pdb):
    #pdbRJBdic[chain][res_number][atom_type/res_type] : res_type=string atom_type keys: [X],[Y],[Z],[occupancy],[Bfactor] floats
    #transform dic into list ordered by atom_numb: list= [atom_numb,atom_type_res_type,chain,res_numb,X,Y,Z,occupancy,Bfactor]
    now = datetime.datetime.now()
    date=now.isoformat()[:10]+' '+now.isoformat()[11:16]
    file_output=open(output_pdb,'w')
    file_output.write('REMARK   PDB written by pdbRJBdic +date+\n')
    for line in pdb_header_list_lines:
        file_output.write(line)
    list_pdbRJBdic=[]
    for chain in pdbRJBdic:
        for res_numb in pdbRJBdic[chain] :
            for atom_type in pdbRJBdic[chain][res_numb]:
                if atom_type!='res_type':
                    atom_list=( pdbRJBdic[chain][res_numb][atom_type]['atom_numb'],atom_type,pdbRJBdic[chain][res_numb]['res_type'],chain,res_numb,pdbRJBdic[chain][res_numb][atom_type]['X'],pdbRJBdic[chain][res_numb][atom_type]['Y'],pdbRJBdic[chain][res_numb][atom_type]['Z'],pdbRJBdic[chain][res_numb][atom_type]['occupancy'],pdbRJBdic[chain][res_numb][atom_type]['Bfactor'] )
                    list_pdbRJBdic.append(atom_list)
    list_pdbRJBdic=sorted(list_pdbRJBdic,key=(lambda item: item[0]), reverse=False)
    for atom in list_pdbRJBdic:
        #list= [atom_numb,atom_type_res_type,chain,res_numb,X,Y,Z,occupancy,Bfactor]
        atom_numb=atom[0]
        atom_type=atom[1]
        res_type=atom[2]
        chain=atom[3]
        res_numb=atom[4]
        X=atom[5]
        Y=atom[6]
        Z=atom[7]
        occupancy=atom[8]
        Bfactor=atom[9]
        atom_numb=' '*(6-len(str(atom_numb)))+str(atom_numb)
        if len(atom_type)<4:
            atom_type='  '+atom_type
        elif len(atom_type)==4:
            atom_type=' '+atom_type
        else:
            print ('Error in function RJB_lib.write_pdb_from_pdbRJBdic')
            print ('Wrong atom_type:',atom_type,'length='+str(len(atom_type)))
            exit()
        atom_type=atom_type+' '*(6-len(atom_type))
        res_numb=' '*(4-len(str(res_numb)))+str(res_numb)
        X='%.3f'%(X)
        X=' '*(7-len(X))+X+' '
        Y='%.3f'%(Y)
        Y=' '*(7-len(Y))+Y+' '
        Z='%.3f'%(Z)
        Z=' '*(7-len(Z))+Z+' '
        occupancy='%.2f'%(occupancy)
        Bfactor='%.2f'%(Bfactor)
        file_output.write('ATOM '+atom_numb+atom_type+res_type+' '+chain+res_numb+'     '+X+Y+Z+' '+occupancy+' '+Bfactor+'           '+atom_type.replace(' ','')[0]+'\n')
    file_output.write('END')
    file_output.close()


def phenix_sec_str_restraints (input_pdb, output_ss , method): # method = [ 'ksdssp' , 'mmtbx_dssp' , 'from_ca' ]
    # general line:  phenix.secondary_structure_restraints model.pdb search_method=from_ca format=pdb extracted from: https://www.phenix-online.org/documentation/reference/secondary_structure_restraints.html
    os.system('phenix.secondary_structure_restraints ' + input_pdb + ' search_method='+method+' format=pdb   > log_ss_'+method+'.log' ) #> /dev/null')
    if os.path.isfile( input_pdb+'_ss.eff'):
        os.system('mv '+input_pdb+'_ss.eff'+' '+output_ss)
    else:
        print ('Failure in RJ_lib.phenix_sec_str_restraints function')
        exit()

def read_pdb_ss_information_return_list_dics (input_ss_eff): #list=[dic(s)] dic keys: 'ss': 'HELIX' or 'SHEET', 'final_res': int, 'initial_res': int, 'length': int, 'chain': letter}
    file=open(input_ss_eff)
    file_list=file.readlines()
    file.close()
    list_file_all=[]
    for line in file_list:
        line=line.split()
        dic_var={}
        dic_var['ss']=line[0]
        if line[0]=='SHEET':
            dic_var['chain']=line[5]
            dic_var['initial_res']=int(line[6])
            dic_var['final_res']=int(line[9])
            dic_var['res_numb_connection_in']=int(line[-5])
            dic_var['chain_connection_in']=line[-6]
            dic_var['res_numb_connection_out']=int(line[-1])
            dic_var['chain_connection_out']=line[-2]
            list_file_all.append(dic_var)
        elif line[0]=='HELIX':
            dic_var['chain']=line[-4]
            dic_var['initial_res']=int(line[5])
            dic_var['final_res']=int(line[-3])
            dic_var['length']=int(line[-1])
            list_file_all.append(dic_var)
    return list_file_all

def read_pdb_ss_information_return_dic_bychain_byelement_ranges (input_ss_eff):
    list_dic=read_pdb_ss_information_return_list_dics (input_ss_eff)
    #print list_dic
    dic_ss_range={}
    for dic in list_dic:
        try:
            dic_ss_range[dic['chain']]
        except:
            dic_ss_range[dic['chain']]={}
            dic_ss_range[dic['chain']]['HELIX']=[]
            dic_ss_range[dic['chain']]['SHEET']=[]
        if len(dic_ss_range[dic['chain']][dic['ss']])>0:
            if dic['initial_res']-1==dic_ss_range[dic['chain']][dic['ss']][-1]:
                dic_ss_range[dic['chain']][dic['ss']].remove(dic_ss_range[dic['chain']][dic['ss']][-1])
        #if dic['initial_res']-1 in dic['chain']][dic['ss']]:
            
        #dic_ss_range[dic['chain']]['COIL']=[]
        for each_number in range(dic['initial_res'],dic['final_res']+1):
            dic_ss_range[dic['chain']][dic['ss']].append( each_number )
    return dic_ss_range

def from_pdbRJBdic_and_dic_ss_range_generate_dic_bychain_sequence (pdbRJBdic,dic_ss_range):
##    dic_ss_range=read_pdb_ss_information_return_dic_bychain_byelement_ranges (input_ss_eff)
##    print dic_ss_range
    dic_ss_seq={}
    for chain in pdbRJBdic:
##        res_numb_old=9999
        try:
            dic_ss_seq[chain]
        except:
            dic_ss_seq[chain]=''
        for res_numb in pdbRJBdic[chain] :
##            if not res_numb==res_numb_old+1:
##                dic_ss_seq[chain]+=' '*(res_numb-res_numb_old)
            try:
                if res_numb in dic_ss_range[chain]['HELIX']:
                    dic_ss_seq[chain]+='H'
                elif res_numb in dic_ss_range[chain]['SHEET']:
                    dic_ss_seq[chain]+='E'
                else:
                    dic_ss_seq[chain]+='C'
            except:
                dic_ss_seq[chain]+='C'
##            res_numb_old=res_numb
        dic_ss_seq[chain]=dic_ss_seq[chain]
    return dic_ss_seq

def from_pdbRJBdic_and_input_ss_eff_generate_dic_bychain_sequence (pdbRJBdic,input_ss_eff):
    dic_ss_range=read_pdb_ss_information_return_dic_bychain_byelement_ranges (input_ss_eff)
    dic_bychain_sequence=from_pdbRJBdic_and_dic_ss_range_generate_dic_bychain_sequence (pdbRJBdic,dic_ss_range)
    return dic_bychain_sequence
#def read_pdb_ss_information_return_list_secstr_seq (input_pdb, input_ss_eff):
    #dic_init=read_pdb_ss_information_return_list_dics (input_ss_eff)
    #pdbRJBdic,header=convert_pdb_into_pdbRJBdic (input_pdb)
    
##    count_S=0
##    count_H=0
##    for line in file_list:
##        line=line.split()
##        dic_var={}
##        dic_var['ss']=line[0]
##        if line[0]=='SHEET':
##            try:
##                dic_all[line[5]]
##            except:
##                dic_all[line[5]]={}
##            chain=line[5]
##            dic_var['initial_res']=line[6]
##            dic_var['final_res']=line[9]
##            dic_var['res_connection']=line[-1]
##            dic_var['chain_connection']=line[-2]
##            dic_all[chain][count_S]=dic_var
##            count_S+=1
##        elif line[0]=='HELIX':
##            try:
##                dic_all[line[-3]]
##            except:
##                dic_all[line[-3]]={}
##            chain=line[-3]
##            dic_var['initial_res']=line[5]
##            dic_var['final_res']=line[-2]
##            dic_var['length']=line[-1]
##            dic_all[chain][count_H]=dic_var
##            count_H+=1
##    print dic_all
##    return dic_all
##


def from_secstr_pred_psipred_generate_trusted_fragments (dic_secstr_pred_psipred,confidence_level,min_size_frag, trust_loop ): #confidence_level_from_1 to 9, should be an integer adding 4Jun16: ['initial_res'] / ['final_res'] #changing in Jul18,2016 to include loops
    min_size_frag=int(min_size_frag)
    confidence_level=int(confidence_level)
    list_dic_ss_frag=[]
    #dic_ss_frag will have same keys as: dic_secstr_pred_psipred keys: ['sequence'] ['sec_str_pred'] ['confidence']
    count_previous=0
    ss_previous=dic_secstr_pred_psipred['sec_str_pred'][0]
    seq=''
    conf=''
    ss=''
    for int_index in range(len(dic_secstr_pred_psipred['sequence'])):
        #print 'index',int_index,'seq',dic_secstr_pred_psipred['sequence'][int_index],'conf',dic_secstr_pred_psipred['confidence'][int_index],'ss',dic_secstr_pred_psipred['sec_str_pred'][int_index]
        #condition to extend fragment of same SS of previous residue but it won't take last residue into consideration
        if (dic_secstr_pred_psipred['sec_str_pred'][int_index]=='H' or dic_secstr_pred_psipred['sec_str_pred'][int_index]=='E' or (trust_loop==True and dic_secstr_pred_psipred['sec_str_pred'][int_index]=='C')) and int(dic_secstr_pred_psipred['confidence'][int_index]) >= confidence_level and ( (dic_secstr_pred_psipred['sec_str_pred'][int_index]==ss_previous and len(seq)>0) or len(seq)==0 ) and int_index!=len(dic_secstr_pred_psipred['sequence'])-1:
            #print 'COND1'
            #dic_secstr_pred_psipred['confidence'][int_index], 'ss', dic_secstr_pred_psipred['sec_str_pred'][int_index]
            seq+=dic_secstr_pred_psipred['sequence'][int_index]
            conf+=dic_secstr_pred_psipred['confidence'][int_index]
            ss+=dic_secstr_pred_psipred['sec_str_pred'][int_index]
        else:
            #when last res is the same ss as previous and should be also counted
            if int_index==len(dic_secstr_pred_psipred['sequence'])-1 and dic_secstr_pred_psipred['sec_str_pred'][int_index]==ss_previous:
                #print 'COND2','index',int_index,'seq',dic_secstr_pred_psipred['sequence'][int_index],'conf',dic_secstr_pred_psipred['confidence'][int_index],'ss',dic_secstr_pred_psipred['sec_str_pred'][int_index]
                seq+=dic_secstr_pred_psipred['sequence'][int_index]
                conf+=dic_secstr_pred_psipred['confidence'][int_index]
                ss+=dic_secstr_pred_psipred['sec_str_pred'][int_index]
            #when other res is found that below 
            if len(seq)>=min_size_frag:
                #print 'COND3'
                dic_ss_frag={}
                dic_ss_frag['sequence']=seq
                dic_ss_frag['sec_str_pred']=ss
                dic_ss_frag['confidence']=conf
                dic_ss_frag['initial_res']=int_index+1-len(seq)
                dic_ss_frag['final_res']=int_index
                list_dic_ss_frag.append(dic_ss_frag)
            seq=dic_secstr_pred_psipred['sequence'][int_index]
            conf=dic_secstr_pred_psipred['confidence'][int_index]
            ss=dic_secstr_pred_psipred['sec_str_pred'][int_index]
        ss_previous=dic_secstr_pred_psipred['sec_str_pred'][int_index]
    #print list_dic_ss_frag
    return list_dic_ss_frag

def SLIDER_read_borfile (bor_file):
    Config = ConfigParser.ConfigParser()
    Config.readfp(open(bor_file))

    #extracting values from Config Parser
    #input files
    pdb_file      = Config.get("GENERAL", "pdb_file")
    hkl_file      = Config.get("GENERAL", "hkl_file")
    mtz_file      = Config.get("GENERAL", "mtz_file")
    input_folder  = Config.get("GENERAL", "working_directory")
    output_folder = Config.get("GENERAL", "output_directory")
    try:    secstr_file    = Config.get("GENERAL", "secstr_file")
    except: secstr_file    = False
    try:    seq_file       = Config.get("GENERAL", "seq_file")
    except: seq_file       = False
    try:
        align_file     = Config.get("GENERAL", "align_file")
        RemoteHomologous=True
    except:
        align_file     = False
        RemoteHomologous= False
    try:    shelxe_ins_file= Config.get("GENERAL", 'shelxe_ins_file')
    except: shelxe_ins_file= False

    try:    ent_file      = Config.get("GENERAL", "ent_file")
    except:
        ent_file=False
        #ent_mtz_file=False
    # if ent_file!=False:
    #     try: ent_mtz_file = Config.get("GENERAL", "ent_mtz_file")
    #     except: ent_mtz_file=False

    try:
        sspdb_file      = Config.get("GENERAL", "sspdb_file")
        sliding_by_ss   = True
    # should be:
    # A: CCHHHHCCEEE
    # B: HHHHHHHHHHH
    except:
        sspdb_file      = False
        sliding_by_ss   = False

    if not RemoteHomologous and not sliding_by_ss:
        print ('Either alignment (pir format) or secondary structure prediction file (PSIPRED) should be provided.')
        exit()
    elif RemoteHomologous==True and sliding_by_ss==True:
        print ('Choose one prior information only, either alignment (PIR format) or secondary structure prediction file (PSIPRED).')

    try:
        refmac_parameters_file = Config.get("GENERAL", "refmac_parameters_file")
    except:
        refmac_parameters_file = ''

    distribute_computing = Config.get("CONNECTION", "distribute_computing")

    #extracting values from Config Parser

    try:
        use_coot = Config.get("SLIDER", "use_coot")
        if use_coot.lower() =='true' or int(use_coot) == 1:   use_coot = True
        else:                                                 use_coot = False
    #parameters_other_programs:
    except:
        use_coot = False

    #input files
    shelxe_path                               = Config.get("LOCAL", "path_local_shelxe")
    scwrl4_path                               = Config.get("LOCAL", "path_local_scwrl4")
    #if ent_mtz_file!=False:
    #      MtzMtzCC_path                         = Config.get("LOCAL", 'path_local_getMtzMtzCC')
    #     PhenixCutMap_path                     = Config.get("LOCAL", 'path_local_PhenixCutMap')
    #else: MtzMtzCC_path    = False
    try:    MtzMtzCC_path                     = Config.get("LOCAL", 'path_local_getMtzMtzCC')
    except: MtzMtzCC_path    = False
    #     PhenixCutMap_path= False
    # MtzMtzCC_path    = False
    PhenixCutMap_path= False
    if use_coot==True:             coot_path  = Config.get("LOCAL", "path_local_coot")
    else:                          coot_path  = False
    edstats_path                              = Config.get("LOCAL", "path_local_edstats")

    #extracting values from Config Parser
    #SLIDER parameters
    # if shelxe_path.endswith('magic'):
    #     try:
    #         fcf_path=Config.get("GENERAL", "fcf_path")
    #     except:
    #         print 'Check if fcf_path was given in GENERAL option in instruction bor file.'
    #         exit()
    # else:
    #     fcf_path=False

    mtz_f    = Config.get("SLIDER", 'f_label')
    mtz_sf   = Config.get("SLIDER", 'sigf_label')
    try:     mtz_i = Config.get("SLIDER", 'i_label')
    except:  mtz_i = False
    try:     mtz_si = Config.get("SLIDER", 'sigi_label')
    except:  mtz_si = False


    mtz_free = Config.get("SLIDER", 'rfree_label')

    try:
        Mw = Config.get("SLIDER", 'molecular_weight')
    except:
        print ('molecular_weight variable not set, use protparam and sequence to obtain protein molecular weight, exiting.')
        exit()
    try:
        nASU = Config.get("SLIDER", 'number_of_component')
    except:
        print ("number_of_component variable not set, use Matthew's Coefficient to obtain most probable number of molecules per Asymmetric Unit, exiting.")
        exit()

    try:     shelxe_line = Config.get("SLIDER", 'shelxe_line')
    except:  shelxe_line = False

    # try:
    #     sliding_by_ss = Config.getint("SLIDER", "sliding_only_ss")
    #     if sliding_by_ss == 1 or sliding_by_ss.lower() == 'true':
    #         sliding_by_ss = True
    #     else:                    sliding_by_ss = False
    #parameters_other_programs:
    # except:
    #     sliding_by_ss = False

    try: RandomModels = Config.getint("SLIDER", "RandomModels")
    except: RandomModels=0
    if RandomModels!=0:
        if RandomModels>1000:
            print ('Random models given exceeds limit, value setting to 1000.')
            RandomModels=1000
        try:
            RandomOnlyEvalSequence=Config.get("SLIDER", "RandomOnlyEvalSequence")
            if RandomOnlyEvalSequence == '1' or RandomOnlyEvalSequence.lower() == 'true' : RandomOnlyEvalSequence = True
            else:                                                                    RandomOnlyEvalSequence = False
        except:                                                                      RandomOnlyEvalSequence = False
    else: RandomOnlyEvalSequence = None

    try:
        ReduceComplexity = Config.get("SLIDER", "ReduceComplexity").lower()
        if ReduceComplexity == '1' or ReduceComplexity == 'true': ReduceComplexity=True
        # if ReduceComplexity!='hydpho' and ReduceComplexity!='lowrot':
        #     print 'ReduceComplexity option in SLIDER should be either hydpho or lowrot'
        #     quit()
    except: ReduceComplexity=False

    # try:
    #     RemoteHomologous = Config.get("SLIDER", "RemoteHomologous")
    #     if RemoteHomologous == '1' or RemoteHomologous.lower() == 'true' : RemoteHomologous = True
    #     else:                                                              RemoteHomologous = False
    # #parameters_other_programs:
    # except:
    #     RemoteHomologous = False

    try:
        Recover0occup= Config.get("SLIDER", "Recover0occup")
        if Recover0occup == '1' or Recover0occup.lower() == 'true': Recover0occup = True
    except:                                                         Recover0occup = False

    try:
        refinement_program = Config.get("SLIDER", "refinement_program") # should be buster / phenix.refine / refmac
        if refinement_program == 'buster':
            try:
                buster_path = Config.get("LOCAL", "path_local_buster")
                refmac_path = False
                phenixrefine_path = False
            except:
                print ('buster path not given in bor file')
                exit()
        elif refinement_program == 'phenix.refine':
            try:
                phenixrefine_path = Config.get("LOCAL", "path_local_phenix.refine")
                refmac_path = False
                buster_path = False
            except:
                print ('phenix.refine path not given in bor file')
                exit()
        elif refinement_program == 'refmac':
            try:
                refmac_path = Config.get("LOCAL", "path_local_refmac")
                phenixrefine_path = False
                buster_path = False
            except:
                print ('refmac path not given in bor file')
                exit()
        else:
            print ('refinement_program should be either buster / phenix.refine / refmac.')
            exit()
    except:
        print ('refinement_program not given in bor file')
        exit()
    try:
        pp_conf=Config.getint("SLIDER", "psipred_confidence_level")
    except:
        pp_conf=0
    try:
        pp_frag_size=Config.getint("SLIDER", "psipred_min_frag_size")
    except:    
        pp_frag_size=4
    try:
        ss_method=Config.get("SLIDER", "ss_eval_pdb_method")
    except:
        ss_method='borgesmatrix'
    
    try:
        minimum_ss_frag=Config.getint("SLIDER", "minimum_ss_frag")
    except:
        minimum_ss_frag=5
    try:
        sliding_tolerance=Config.getint("SLIDER", "sliding_tolerance")
    except:
        sliding_tolerance=3
    # try:
    #     size_frag_tolerance=Config.getint("SLIDER", "size_frag_tolerance")
    # except:
    #     size_frag_tolerance=3
    #size_frag_tolerance=None
    try:
        models_by_chain=Config.getint("SLIDER", "models_by_chain")
        if models_by_chain>1000:
            print ('models_by_chain given exceeds limit, value setting to 1000.')
            models_by_chain=1000
    except:
        models_by_chain=100

    try:
        seq_pushed_refinement=Config.getint("SLIDER", "seq_pushed_refinement")
        if seq_pushed_refinement>1000:
            print ('seq_pushed_refinement given exceeds limit, value setting to 1000.')
            seq_pushed_refinement=1000
    except:
        seq_pushed_refinement=models_by_chain
##    try:
##        models_combination=Config.getint("SLIDER", "models_combination")
##    except:
##        models_combination=3

    try:
        Config.get("SLIDER", "chosen_chains_independent")
        print ('Obsolete variable "chosen_chains_independent" was used in borfile, please use "chosen_chains" and separate by , dependent chains and ; independent chains.')
        exitt = True
    except: exitt = False
    if exitt==True: exit()

    try:
        Config.get("SLIDER", "chosen_chains_dependent")
        print ('Obsolete variable "chosen_chains_dependent" was used in borfile, please use "chosen_chains" and separate by , dependent chains and ; independent chains.')
        exitt = True
    except: exitt = False
    if exitt==True: exit()

    # chosen_chains_independent=[]
    # try:
    #     for i in Config.get("SLIDER", ":q
    # :q:qains_independent").replace(' ','').split(','):
    #         chosen_chains_independent.append( [i] )
    # except:
    #     pass
    # chosen_chains_dependent=[]
    # try:
    #     chosen_chains_dependent.append ( Config.get("SLIDER", "chosen_chains_dependent").replace(' ','').split(',') )
    # except:
    #     for i in range(1,10):
    #         try:
    #             chosen_chains_dependent.append ( Config.get("SLIDER", "chosen_chains_dependent"+str(i)).replace(' ','').split(',') )
    #         except:
    #             pass

    chosen_chains=[]
    try:
        for i in Config.get("SLIDER", "chosen_chains").replace(' ','').split(';'):
            chosen_chains.append( i.split(',') )
    except:
        pass

    try:
        fixed_residues_modelled=Config.get("SLIDER", "fixed_residues_modelled") # should be given as chain:residueN1-residueN2,residueN3-residueN4 such as: A:1-10,20-30 B:1-5,30-35
        #DicFixResMod={}
        DicFixResMod =defaultdict(list)
        FixRes2=fixed_residues_modelled.split()
        for s in FixRes2:
            #DicFixRes[s[0]]=[]
            for ss in s[2:].split(','):
                resn1= int(ss.split('-')[0])
                resn2 =int(ss.split('-')[1])
                for i in range( resn1 , resn2+1 ):
                    DicFixResMod[s[0]].append( i )
    except: DicFixResMod=defaultdict(list)

    try:
        fixed_residues_notmodelled=Config.get("SLIDER", "fixed_residues_notmodelled") # should be given as chain:residueN1-residueN2,residueN3-residueN4 such as: A:1-10,20-30 B:1-5,30-35
        #DicFixResNotMod={}
        DicFixResNotMod =defaultdict(list)
        FixRes3=fixed_residues_notmodelled.split()
        for s in FixRes3:
            #DicFixRes[s[0]]=[]
            for ss in s[2:].split(','):
                resn1= int(ss.split('-')[0])
                resn2 =int(ss.split('-')[1])
                for i in range( resn1 , resn2+1 ):
                    DicFixResNotMod[s[0]].append( i )
    except: DicFixResNotMod=defaultdict(list)

    for ch in DicFixResMod:
        for i in DicFixResMod:
            if i in DicFixResNotMod[ch]:
                print ('Chosen fixed_residues_modelled and fixed_residues_notmodelled cannot have overlap. They have same residue',i,'in chain',ch,'.')
                exit()

    # if DicFixRes!=False:
    #     try:
    #         ModelFixedResidues = Config.get("SLIDER", "ModelFixedResidues")
    #         if 'true' == ModelFixedResidues.lower() or '1' == ModelFixedResidues: ModelFixedResidues = True
    #         else:                                                                 ModelFixedResidues = False
    #     except: ModelFixedResidues=False
    #
    # else: ModelFixedResidues=False

    try:
        ncschains=Config.get("SLIDER", "ncschains")
        if 'true'==ncschains.lower() or '1'==ncschains: ncschains=True
        else: ncschains=False
    except:   ncschains=False

    try:
        merge_chosen_chains=Config.get("SLIDER", "merge_chosen_chains")
        merge_chosen_chains=merge_chosen_chains.replace(' ','').split(',')
    except:
        merge_chosen_chains=[]
    try:
        models_by_merge_chain=Config.getint("SLIDER", 'models_by_merge_chain')
    except:
        models_by_merge_chain=int(models_by_chain)
    models_by_merge_chain=int(math.sqrt(models_by_merge_chain))


    try:
        number_shelxe_trials=Config.getint("SLIDER", 'number_shelxe_trials')
    except:
        number_shelxe_trials=15

    try:
        expand_from_map = Config.get("SLIDER", "expand_from_map")
        if 'true'==expand_from_map.lower() or '1'==expand_from_map:                                 expand_from_map=True
        else:                                                                                       expand_from_map=False
    except:
        if refinement_program=='buster': expand_from_map=True
        else:                            expand_from_map=False


    try:
        trust_loop=Config.get("SLIDER", "trust_loop")
        if trust_loop=='1' or trust_loop.lower()=='true':
            trust_loop=True
        elif trust_loop=='-1' or trust_loop=='0' or trust_loop.lower()=='false':
            trust_loop=False
    except:
        trust_loop=False

    try:
        ModelEdge=Config.get("SLIDER", "ModelEdge")
        if    ModelEdge == '1' or ModelEdge.lower() == 'true':  ModelEdge = True
        elif  ModelEdge == '0' or ModelEdge.lower() == 'false': ModelEdge = False
        else:
            print ('Invalid option given for ModelEdge',ModelEdge)
            exit()
    except: ModelEdge = False

    try:
        buster_parameters=Config.get("SLIDER", "buster_parameters")
    except:
        buster_parameters='-noWAT -nbig 10 -RB -nthread 1 UsePdbchk="no"'

    try:
        PhenixRefineParameters=Config.get("SLIDER", "PhenixRefineParameters")
    except:
        PhenixRefineParameters=False

    try:
        merge_FOM=Config.get("SLIDER", "merge_FOM")
    except:
        merge_FOM='refine1_LLG'

    try:
        LLG_testing=Config.get("SLIDER", "LLG_testing")
        if    LLG_testing == '1' or LLG_testing.lower() == 'true':  LLG_testing = True
        elif  LLG_testing == '0' or LLG_testing.lower() == 'false': LLG_testing = False
    except: LLG_testing=False

    # try:
    #     control_folder = Config.get("SLIDER", "control_folder")
    # #parameters_other_programs:
    # except:
    #     control_folder = False

    dic_borfile={ 'pdb_file': pdb_file,'hkl_file': hkl_file, 'mtz_file': mtz_file, 'shelxe_ins_file': shelxe_ins_file, 'fcf_path':fcf_path,"ent_path":ent_file,'ent_mtz_file':ent_mtz_file,'sspdb_file':sspdb_file,'output_folder': output_folder, 'secstr_file': secstr_file, 'seq_file':seq_file, 'align_file':align_file,'psipred_min_frag_size': pp_frag_size , 'ss_eval_pdb_method': ss_method, 'psipred_confidence_level': pp_conf , 'refinement_program' : refinement_program , 'refmac_parameters_file': refmac_parameters_file,'sliding_by_ss': sliding_by_ss,'RandomModels':RandomModels,'RandomOnlyEvalSequence':RandomOnlyEvalSequence,'ReduceComplexity':ReduceComplexity,'RemoteHomologous':RemoteHomologous,'Recover0occup':Recover0occup,'minimum_ss_frag':minimum_ss_frag, 'size_frag_tolerance':size_frag_tolerance, "models_by_chain":models_by_chain, 'seq_pushed_refinement':seq_pushed_refinement,'chosen_chains':chosen_chains, 'fixed_residues':DicFixRes,'ModelFixedResidues':ModelFixedResidues,"ncschains":ncschains,"shelxe_path": shelxe_path, "phenixrefine_path": phenixrefine_path , "refmac_path":refmac_path , 'buster_path':buster_path , 'scwrl4_path':scwrl4_path , 'coot_path':coot_path , 'edstats_path':edstats_path , 'MtzMtzCC_path':MtzMtzCC_path,'PhenixCutMap_path':PhenixCutMap_path, 'distribute_computing':distribute_computing, 'shelxe_line':shelxe_line , 'number_shelxe_trials':number_shelxe_trials ,'expand_from_map':expand_from_map,'models_by_merge_chain':models_by_merge_chain, 'merge_chosen_chains':merge_chosen_chains, 'trust_loop':trust_loop, 'ModelEdge':ModelEdge, 'sliding_tolerance':sliding_tolerance  ,     'mtz_f':mtz_f , 'mtz_sf':mtz_sf , 'mtz_i':mtz_i , 'mtz_si':mtz_si , 'mtz_free':mtz_free  , 'buster_parameters':buster_parameters , 'PhenixRefineParameters':PhenixRefineParameters,'merge_FOM': merge_FOM , 'use_coot':use_coot , 'nASU':nASU , 'Mw':Mw , 'LLG_testing':LLG_testing }
    return dic_borfile

def extract_chains_from_pdb(pdb_input):
    list_ordered_chains=[]
    file_input=open(pdb_input,"r")
    pdb_file_list=file_input.readlines()
    file_input.close()
    for line in pdb_file_list:
        if line.startswith("ATOM") and line[12:16]==' CA ' and line[21] not in list_ordered_chains:
            list_ordered_chains.append(line[21])
    return list_ordered_chains


def counter_consecutive_letters(string):
    previous_letter = ''
    counter = []
    for letter in string:
        if letter == previous_letter:
            counter[-1] = (letter, counter[-1][1] +1)
        else:
            counter.append((letter, 1))
            previous_letter = letter
    return counter

def separate_consecutive_letters_string_into_list(string):
    a=counter_consecutive_letters(string)
    list=[]
    for item in a:
        l=item[0]*item[1]
        list.append(l)
    return list

def convert_ss_borgesmatrix_dic_into_dic_bychain_byelement_ranges (ss_borgesmatrix_dic): #generate Jun2,2016 to convert ss_borgesmatrix_dic into my dic
    dic_ss_range={}
    for dic in ss_borgesmatrix_dic:
        for tuple_residue in dic['reslist']:
            chain=tuple_residue[2]
            res_numb=tuple_residue[3][1]
            try:
                dic_ss_range[chain]
            except:
                dic_ss_range[chain]={}
                dic_ss_range[chain]['HELIX']=[]
                dic_ss_range[chain]['SHEET']=[]
            if dic['sstype']=='ah':
                dic_ss_range[chain]['HELIX'].append(res_numb)
            elif dic['sstype']=='bs':
                dic_ss_range[chain]['SHEET'].append(res_numb)
    return dic_ss_range

def from_ss_borgesmatrix_and_pdbRJBdic_return_dic_bychain_sequence (ss_borgesmatrix_dic,pdbRJBdic):
    dic_ss_range=convert_ss_borgesmatrix_dic_into_dic_bychain_byelement_ranges (ss_borgesmatrix_dic)
    dic_bychain_sequence=from_pdbRJBdic_and_dic_ss_range_generate_dic_bychain_sequence (pdbRJBdic,dic_ss_range)
    return dic_bychain_sequence

def given_string_return_list_same_characters (string):
    list=[]
    prev=string[0]
    var=''
    for char in string:
        if char==prev:
            var+=char
        else:
            list.append(var)
            var=char
        prev=char
    list.append(var)
    return list











def generate_fragments_restricting_by_ss ( pdb_file , pdb_dic_seq_ss , list_dics_sspred_conf , minimum_ss_frag , size_frag_tolerance , chosen_chains ):
    list_seqfrags_H=[]
    list_dicfrags_H=[]
    list_seqfrags_S=[]
    list_dicfrags_S=[]
    for dic in list_dics_sspred_conf:
        if dic['sec_str_pred'][0]=='H':
            list_dicfrags_H.append(dic)
            list_seqfrags_H.append(dic['sequence'])
        elif dic['sec_str_pred'][0]=='S':
            list_dicfrags_S.append(dic)
            list_seqfrags_S.append(dic['sequence'])
    #chain_list=return_chain_list (pdb_file)
    resnumb_dic=return_dic_resnumb_list (pdb_file)
    seq_dic=return_dic_sequence (pdb_file)
    dic_variation={}
    print ('sspred')
    for i in list_dics_sspred_conf:
        print (i)
##    list_seq_keep=[]
##    list_chain_keep=[]
##    super_list_seq=[]
    count_chain=0
    #for chain in chain_list:
    for chain in chosen_chains:
##        print ('\n\nYEY\n\n', chain)
        res_count=0
        seq=seq_dic[chain]
        ss_string=pdb_dic_seq_ss[chain]
        ss_list=separate_consecutive_letters_string_into_list(ss_string)
        list_chain=[]
        for piece_ss in ss_list:
##            print ('piece_ss',piece_ss,'length',len(piece_ss))
##            print ('rescount',res_count)
            list_var=[]
            res_count_end=res_count+len(piece_ss)
##            list_chain_keep.append(chain)
##            list_seq_keep_var=chain_seq[res_count :res_count_end].lower()
##                print (list_var)
            if piece_ss[0]=='H' and len(piece_ss)>=minimum_ss_frag:
##                print 'H'
                count_H_frag=0
                bigger_fragH=max(list_seqfrags_H, key=len)
                for fragH in list_seqfrags_H:
##                    print ('fragH',fragH)
##                    print ('fragH',len(fragH),'piece_ss',len(piece_ss))
                    #prob= 1 /  float (  (abs( len(fragH)-len(piece_ss))+1 ))
                    prob= 1 /  float (  (abs( len(fragH)-len(piece_ss))+1 ))
##                    print ('prob',prob)
                    score=math.log( prob )
                    #score=score / len(piece_ss)
                    if len(fragH)>= (len(piece_ss)-size_frag_tolerance):
                        if len(fragH)<len(piece_ss):
                            begin=len(fragH)-len(piece_ss)
                            end=0
                        else:
                            begin=0
                            end=len(fragH)-len(piece_ss)+1
                        for c in range( begin,end ):
                            if c==begin: #this is the condition starting sequence... So size_frag_tolerance may generate A in beginning of seq
##                                print ('cond1')
                                for i in range(1,size_frag_tolerance):
                                    seq='a'*i+fragH[:len(piece_ss)-i]
                                    if len(seq)<len(piece_ss):
                                        seq+='a'*(len(piece_ss)-len(seq))
                                    init_res=list_dicfrags_H[count_H_frag]['initial_res']-c-i
                                    final_res=init_res+len(piece_ss)-1
                                    score2=score+math.log( 1/float( (i+1)**2 ) )
                                    list_var.append( (seq,init_res,final_res,score2) )
##                                    print (seq)
                            if c==(end-1): #this is the condition ending sequence... So size_frag_tolerance may generate A in ending of seq
##                                print ('cond2')
                                for i in range(1,size_frag_tolerance):
                                    seq=fragH[-len(piece_ss)+i:]
                                    seq+='a'*i
                                    #seq+='a'*(len(piece_ss)-len(seq))
                                    init_res=list_dicfrags_H[count_H_frag]['initial_res']+c-i
                                    final_res=init_res+len(piece_ss)-1
                                    score2=score+math.log( 1/float((i+1)**2) )
                                    if len(seq)==len(piece_ss):
                                        list_var.append( (seq,init_res,final_res,score2) )
##                                        print (seq)
##                                    else:
##                                        print ('\n\n\nREJECTED!!!!\n',seq,'\n\n\n')
                            if c>=0:
##                                print ('cond3')
                                seq=fragH[c:len(piece_ss)+c]
                                init_res=list_dicfrags_H[count_H_frag]['initial_res']+c
                                final_res=init_res+len(piece_ss)-1
                                list_var.append( (seq,init_res,final_res,score) )
##                                print (seq)
##                        print ('list_var',list_var)
                    elif len(bigger_fragH)<len(piece_ss):
                        if len(bigger_fragH)-size_frag_tolerance<len(fragH):
                            print ('\nFailure relating SS psipred prediction and SS pdb evaluation.\nA helix with length',len(piece_ss),'was found in chain',chain,'and it is bigger than biggest predicted helix (',len(bigger_fragH),')by PSIPRED. Thus, polyAla will be completed in edges in all possibilities of the helix:',fragH,'with length',len(fragH))
                            for c in range( -len(fragH)+len(piece_ss)+1 ):
                                d=-len(fragH)+len(piece_ss)-c
                                seq='a'*c+fragH+'a'*d
                                init_res=list_dicfrags_H[count_H_frag]['initial_res']-c
                                final_res=init_res+len(piece_ss)-1
                                list_var.append( (seq,init_res,final_res,score) )
    ##                        print ('list_var special',list_var)
                    count_H_frag+=1
##                    if kill:
##                        exit()
            elif piece_ss[0]=='E' and len(piece_ss)>=minimum_ss_frag:
                print ('In preparation')
                exit()
            else:   #this would be kind of equal to 'if piece_ss[0]=='C':', although len(frag)>=len(minimum_frag) has to be taken into account also
##                print ('C')
                seq_var=seq_dic[chain][res_count :res_count_end].lower()
                list_var.append((seq_var,0,0,0))

            list_chain.append(list_var)
##            dic_var['possibilities']=list_var
##            dic_var['frag_numb']=frag_numb
##            list_dic_variation.append(dic_var)
##            super_list_seq.append(list_var)
##            list_seq_keep.append(list_seq_keep_var)
            res_count+=len(piece_ss)
        dic_variation[chain]=list_chain
    return dic_variation

##def scoring_function_logprobseq (pred_frag,obs_frag,prev_score):
##    try:
##        score=prev_score+math.log(1/abs(len(pred_frag)-len(obs_frag)))
##    except:
##        score=math.log(1/abs(len(pred_frag)-len(obs_frag)))
##    return score


def generate_fragments_restricting_by_ss_new ( pdb_file , pdb_dic_seq_ss , list_dics_sspred_conf , minimum_ss_frag , size_frag_tolerance , chosen_chains , show_state=False ): # modified 14 Jul 2016 to go in a more general approach
##    sspred_list_H=[]
##    sspred_dic_H=[]
##    sspred_list_S=[]
##    sspred_dic_S=[]
    #example list_dics_sspred_conf=[{'confidence': '06688524489999999999998230', 'sec_str_pred': 'HHHHHHHHHHHHHHHHHHHHHHHHHH', 'initial_res': 9, 'final_res': 34, 'sequence': 'NEAWVKDTNGFDILMGQFAHNIENIW'}, {'confidence': '4752', 'sec_str_pred': 'EEEE', 'initial_res': 47, 'final_res': 50, 'sequence': 'YVKY'}, {'confidence': '0543', 'sec_str_pred': 'EEEE', 'initial_res': 58, 'final_res': 61, 'sequence': 'SHIN'}, {'confidence': '08998', 'sec_str_pred': 'EEEEE', 'initial_res': 66, 'final_res': 70, 'sequence': 'TITIE'}, {'confidence': '2899999888210', 'sec_str_pred': 'HHHHHHHHHHHHH', 'initial_res': 77, 'final_res': 89, 'sequence': 'PAAHLRRAIIKTL'}, {'confidence': '56750045899999984211', 'sec_str_pred': 'HHHHHHHHHHHHHHHHHHHH', 'initial_res': 130, 'final_res': 149, 'sequence': 'EGRASNFADYLLKNRLKSRS'}, {'confidence': '103442', 'sec_str_pred': 'HHHHHH', 'initial_res': 155, 'final_res': 160, 'sequence': 'IYSVTI'}]
 
    #pdb_resnumb_dic=return_dic_resnumb_list (pdb_file)
    pdb_seq_dic=return_dic_sequence (pdb_file)
    dic_list_all_seq_variation={} #key chain letter
##    list_seq_keep=[]
##    list_chain_keep=[]
##    super_list_seq=[]
    count_chain=0
    #for chain in chain_list:
    if show_state:
        print ('sspred is')
    for i,dic_sspred in enumerate(list_dics_sspred_conf):
        if show_state:
            print (dic_sspred['sec_str_pred']  )
        list_dics_sspred_conf[i]['index']=i
    if show_state:
        print ('\n' )
    for list_chain in chosen_chains:
        for chain in list_chain:
            if show_state:
                print ('\n\nchain is',chain)
    ##        print '\n\nYEY\n\n', chain
            res_count=0
            pdb_seq=pdb_seq_dic[chain]
            pdb_ss_string=pdb_dic_seq_ss[chain]
            pdb_ss_list=separate_consecutive_letters_string_into_list(pdb_ss_string)
            pdb_seq_list=given_string_convert_in_list_string_based_on_other_list_strings ( pdb_seq , pdb_ss_list )
            if show_state:
                print ('\npdb_ss_list is',pdb_ss_list)
            list_var_chain=[]
            list_clean=[]
            if show_state:
                print ('pdb_seq_list',pdb_seq_list)
                print ('pdb_ss_list',pdb_ss_list)
            for ind_piece_pdb,piece_pdb in enumerate(pdb_ss_list):
                if show_state:
                    print ('\n\n\npiece_pdb is',piece_pdb)
                list_by_piece=[]
                if show_state:
                    print ('accepted minimum legth',len(piece_pdb)>=minimum_ss_frag, len(piece_pdb),'>=',minimum_ss_frag)
                if len(piece_pdb)>=minimum_ss_frag:
                    for dic_sspred in list_dics_sspred_conf:
                        if show_state:
                            print ('\neval',dic_sspred['sec_str_pred'],len(dic_sspred['sec_str_pred']))
                            print ('same ss',dic_sspred['sec_str_pred'][0]==piece_pdb[0])
                            print ('size in accordance',(len(dic_sspred['sec_str_pred'])+size_frag_tolerance)>=len(piece_pdb))
                        if dic_sspred['sec_str_pred'][0]==piece_pdb[0] and (len(dic_sspred['sec_str_pred'])+size_frag_tolerance)>=len(piece_pdb):
                            dic_new=dict(dic_sspred)
                            dic_new['pdb_piece_seq']=pdb_seq_list[ind_piece_pdb]
                            list_by_piece.append(dic_new)
                            if show_state:
                                print ('frag inserted',dic_sspred['sec_str_pred'])
                    if len(list_by_piece)==0: #a loop a bigger fragment than prediction
                        list_by_piece=[{'sec_str_pred':piece_pdb,'pdb_piece_seq':pdb_seq_list[ind_piece_pdb]}]
                        if show_state:
                            print ('loop',piece_pdb,'placed')
                else:
                    list_by_piece=[{'sec_str_pred':piece_pdb,'pdb_piece_seq':pdb_seq_list[ind_piece_pdb]}]
                    if show_state:
                        print ('number of residues not long enough',piece_pdb,'placed')
                if len(list_by_piece)==0: # condition with no pdbFrag has the size of any ss prediction
                    chosen_pdb_frags=[]
                    bigger_frag=0
                    for ind_piece_pdb,piece_pdb in enumerate(pdb_ss_list):
                        if dic_sspred['sec_str_pred'][0]==piece_pdb[0]:
                            if bigger_frag<len(dic_sspred['sec_str_pred']): bigger_frag=len(dic_sspred['sec_str_pred'])
                    for ind_piece_pdb,piece_pdb in enumerate(pdb_ss_list):
                        if dic_sspred['sec_str_pred'][0]==piece_pdb[0] and (len(dic_sspred['sec_str_pred'])+size_frag_tolerance)>= bigger_frag:
                            dic_new=dict(dic_sspred)
                            dic_new['pdb_piece_seq']=pdb_seq_list[ind_piece_pdb]
                            list_by_piece.append(dic_new)
                list_var_chain.append(list_by_piece)
            if show_state:
                print ('\n\n\npdb_ss_list is',pdb_ss_list)
                print ('\nPossibilities')
                for ind,i in enumerate(list_var_chain):
                    for ii in i:
                        try:
                            print (ind,ii['sec_str_pred'])
                        except:
                            print (ind,'no assignment',ii)
            for each_possib in itertools.product(*list_var_chain):
                prev_i=-1
                ind_check=True
                for i,frag in enumerate(each_possib):
                    if ind_check:
                        try:
                            if not frag['index']>prev_i:
                                ind_check=False
                            prev_i=frag['index']
                        except:
                            pass
                if ind_check:
                    list_clean.append(each_possib)
            if show_state:
                print ('\n\n\nshowing results')
                for pos in list_clean:
                    for p in pos:
                        print (p)
                    print ('\n\n')
            dic_list_all_seq_variation[chain]=list_clean
    return dic_list_all_seq_variation

def given_string_convert_in_list_string_based_on_other_list_strings ( string , list_string ):
    l=[]
    c=0
    for s in list_string:
        l.append(string[c:(c+len(s))])
        c+=len(s)
    if not len(string)==c:
        print ('RJB_lib.given_string_convert_in_list_string_based_on_other_list_strings function given with different number variables: #characters in list = ',c,'#characters in string',len(string))
        exit()
    return l

def given_dic_variation_sspred_groups_generated_sequences ( dic_list_all_seq_variation , minimum_ss_frag , size_frag_tolerance , diag=False ):
    dic_all_possib={}
    #dic_list_all_seq_variation=[{'sec_str_pred':piece_pdb,'pdb_piece_seq':pdb_seq_list[ind_piece_pdb]}] in case loop or smaller pdb_frag than given tolerance
    #dic_list_all_seq_variation=[{'sec_str_pred':'HHHHHHHH','pdb_piece_seq':'ACYRL', 'confidence': '06688524489999999999998230', 'sec_str_pred': 'HHHHHHHHHHHHHHHHHHHHHHHHHH', 'initial_res': 9, 'final_res': 34, 'sequence': 'NEAWVKDTNGFDILMGQFAHNIENIW'}, {'confidence': '4752', 'sec_str_pred': 'EEEE', 'initial_res': 47, 'final_res': 50, 'sequence': 'YVKY'}, {'confidence': '0543', 'sec_str_pred': 'EEEE', 'initial_res': 58, 'final_res': 61, 'sequence': 'SHIN'}, {'confidence': '08998', 'sec_str_pred': 'EEEEE', 'initial_res': 66, 'final_res': 70, 'sequence': 'TITIE'}, {'confidence': '2899999888210', 'sec_str_pred': 'HHHHHHHHHHHHH', 'initial_res': 77, 'final_res': 89, 'sequence': 'PAAHLRRAIIKTL'}, {'confidence': '56750045899999984211', 'sec_str_pred': 'HHHHHHHHHHHHHHHHHHHH', 'initial_res': 130, 'final_res': 149, 'sequence': 'EGRASNFADYLLKNRLKSRS'}, {'confidence': '103442', 'sec_str_pred': 'HHHHHH', 'initial_res': 155, 'final_res': 160, 'sequence': 'IYSVTI'}]
    for chain in dic_list_all_seq_variation:
        if diag:
            print (chain  )
        pos_all_ch=[]
        for list_dic in dic_list_all_seq_variation[chain]:
            if diag:
                print (list_dic)
            l=[]
            for dic in list_dic:
                list_var=[]
                piece_ss=dic['pdb_piece_seq']
                sspred=dic['sec_str_pred']
                if len(dic)==2: #if no frag will be assigned
                    list_var.append( ( piece_ss.lower() ,0,0) ) #0,0 should be changed by initial and final res numbers from pdb
                elif len(sspred)>=len(piece_ss)-size_frag_tolerance and len(dic)>2:
                    assign_seq=dic['sequence']
                    if len(assign_seq)<len(piece_ss):
                        begin=len(assign_seq)-len(piece_ss)
                        end=0
                    else:
                        begin=0
                        end=len(assign_seq)-len(piece_ss)+1
                    for c in range( begin,end ):
                        if c>=0: #normal condition
    ##                                print 'cond3'
                            seq=assign_seq[c:len(piece_ss)+c]
                            init_res=dic['initial_res']+c
                            final_res=init_res+len(piece_ss)-1
                            list_var.append( (seq,init_res,final_res) )
                        if c==begin: #this is the condition starting sequence... So size_frag_tolerance may generate A in beginning of seq
    ##                                print 'cond1'
                            for i in range(1,size_frag_tolerance+1):
                                seq='a'*i+assign_seq[:len(piece_ss)-i]
                                if len(seq)<len(piece_ss):
                                    seq+='a'*(len(piece_ss)-len(seq))
                                init_res=dic['initial_res']-c-i
                                final_res=init_res+len(piece_ss)-1
                                list_var.append( (seq,init_res,final_res) )
    ##                                    print seq
                        if c==(end-1): #this is the condition ending sequence... So size_frag_tolerance may generate A in ending of seq
    ##                                print 'cond2'
                            for i in range(1,size_frag_tolerance+1):
                                seq=assign_seq[-len(piece_ss)+i:]
                                seq+='a'*i
                                #seq+='a'*(len(piece_ss)-len(seq))
                                init_res=dic['initial_res']+c-i
                                final_res=init_res+len(piece_ss)-1
                                if len(seq)==len(piece_ss):
                                    list_var.append( (seq,init_res,final_res) )
                elif len(sspred)<len(piece_ss)-size_frag_tolerance and len(dic)>2:
                    assign_seq=dic['sequence']
                    print ('\nFailure relating SS psipred prediction and SS pdb evaluation.\nA helix with length',len(piece_ss),'was found in chain',chain,'and it is bigger than biggest predicted helix by PSIPRED. Thus, polyAla will be completed in edges in all possibilities of the helix:',assign_seq,'with length',len(assign_seq))
                    for c in range( -len(assign_seq)+len(piece_ss)+1 ):
                        d=-len(assign_seq)+len(piece_ss)-c
                        seq='a'*c+assign_seq+'a'*d
                        init_res=dic['initial_res']-c
                        final_res=init_res+len(piece_ss)-1
                        list_var.append( (seq,init_res,final_res) )
    ##                        print 'list_var special',list_var
                else:
                    print ('This should not have passed, please evaluate cause.')
                    exit()
                l.append(list_var)
            pos=itertools.product(*l)
            for p in pos:
                pos_all_ch.append(p)
        dic_all_possib[chain]=pos_all_ch

        if diag:
            print (chain)
            for each in pos_all_ch:
                print (each)
                seq=''
                for i in each:
                    seq+=i[0]
                print (seq )
                        
    if diag:
        for l in dic_list_all_seq_variation:
            for ll in l:
                for key in ll:
                    print (ll[key],key)
                print ('\n\n')

    return dic_all_possib



def generate_sequence_combination_restricted_by_tolerance ( dic_variation , sliding_tolerance , length_seq ):
    dic_possib={}
##    print '\n\n\n\nNOW IT COMES'
    for chain in dic_variation:
        possib_by_chain=dic_variation[chain]
        list_var=[]
##        print possib_by_chain
        for each_possib in itertools.product(*possib_by_chain):
            seq=each_possib[0][0]
            score=each_possib[0][-1]
            length=len(seq)
            init_res1=each_possib[0][1]
            final_res1=each_possib[0][2]
##            print each_possib
##            print score
            if final_res1==0: #correcting final residue for loop
                final_res1=length
##            print init_res1,final_res1
            for possib in each_possib[1:]:
##                print possib
                score+=possib[-1]
##                print score
                seq2=possib[0]
                length+=len(seq2)
                init_res2=possib[1]
                final_res2=possib[2]
                if init_res2==0 and final_res2==0: #correcting final residue for loop
                    final_res2=len(seq2)+final_res1
                    init_res2=final_res1+1
##                print init_res2,final_res2
                if final_res1<=init_res2+sliding_tolerance and final_res2>=length-sliding_tolerance:
                    if not init_res1==0:
##                        print 'final_res1',final_res1,'init_res2',init_res2
                        prob= 1 /  float (  (abs( final_res1-init_res2+1)+1 ))
##                        print 'prob',prob
                        score+=math.log( prob )
                        #score+=math.log(1/((abs(final_res1-init_res2+1))+1))
                    seq+=seq2
                    init_res1=init_res2
                    final_res1=final_res2
##                    print score
            if len(seq)==length:
                list_var.append((seq,score))
##                print score
##                print 'YEAH!!\n\n\n\n',seq
        dic_possib[chain]=list_var
    print ('\n\nConsidering initial and final residues of each fragment and its sequence assignment, the following was generated by chain:')
    for chain in sorted(dic_possib):
        print ('chain',chain)
##        print list_possib[i][0][0]
        print ('number of assigned residues',len(filter(lambda x: x in string.uppercase, dic_possib[chain][-1][0])))
        print ('number of possible trials:',len(dic_possib[chain]),'\n')
##        print 'Seq\tScore'
##        for ii in list_possib[i]:
##            print ii[0],'\t\t',ii[1]
##        print '\n\n\n\n'
##        for each in list_possib[i]:
##            print each[1]
    #print list_possib
    return dic_possib

def organized_possibilities_list_restrict_by_given_numbers ( dic_possib , models_by_chain ): #changed in 27Jun2016 to input dic_possib
    short_dic_possib={}
    for chain in dic_possib:
        chain_list=dic_possib[chain]
        #chain_list=sorted(chain_list,key=(lambda item: item[1]), reverse=True)
        cut=models_by_chain
        try:
            possib_ind=models_by_chain
            while chain_list[possib_ind][1]==chain_list[possib_ind-1][1] and chain_list[possib_ind][2]==chain_list[possib_ind-1][2] and possib_ind<len(chain_list):
                cut=possib_ind
                possib_ind+=1
        except:
            pass
        short_dic_possib[chain]=chain_list[:cut]
    print ('\n\nConsidering the scoring function and the chosen number of models to be pushed through refinement, the following was generated by chain:')
    for chain,lista in short_dic_possib.items():
        print ('chain',chain)
        print ('number of assigned residues',max ([ len( filter(lambda y: y in string.uppercase , x[0]) ) for x in lista]))
        print ('number of possible trials:',len(lista),'\n')
    return short_dic_possib

def generate_sequences_from_seqcombination_chosenchains ( dic_possib , pdb_file , chosen_chains , models_by_chain , BORGES_MATRIX_listSSpdb , ncs=False ): #script generated for seqslider_uk1.2
    seq_list=[]
    pdb_seq_list=return_list_sequence(pdb_file)
    pdb_chain_list=return_chain_list(pdb_file)
    #create list containing possibilities
    for list_chain in chosen_chains:
        l=[]
        l3=[]
        #account more than one chain
        if len(list_chain)>1:
            #print len(list_chain)
            n=0
            if not ncs:
                while abs( (n**len(list_chain) - models_by_chain) ) > abs( ((n+1)**len(list_chain) - models_by_chain) ) :
                    #print n,n**len(list_chain),abs( (n**len(list_chain) - models_by_chain) ),'>',abs( ((n+1)**len(list_chain) - models_by_chain) )
                    n+=1
            #no combinations for NCS
            else:
                if models_by_chain>len(dic_possib[chosen_chains[0][0]]): n=len(dic_possib[chosen_chains[0][0]])
                else:                                                    n=models_by_chain
            for chain in list_chain:
                l2=[]
                for tuple_eval in dic_possib[chain][:n]:
                    t2=list(tuple_eval)
                    t2.append(chain)
                    l2.append(tuple(t2))
                l3.append(l2)
            if not ncs: l=itertools.product(*l3)
            else:
                #print l3
                for i in range(n):
                    #print i
                    l4=[]
                    for ii in range(len(l3)):
                        #print ii
                        #print l3[ii][i]
                        l4.append(l3[ii][i])
                    l.append( l4 )
        #account for only one chain
        else:
            for i in dic_possib[ list_chain[0] ]:
                t2=list(i)
                t2.append(list_chain[0])
                l.append( [tuple(t2)] )
        # generate tuple
        for l_tuple_eval in l:
            seq=''
            ll_score=0
            alig_score=0
            list_seq_chain=[]
            if len(l_tuple_eval[0])==5:
                id=0
                idres=0
                count=0
            for pdb_ind in range(len(pdb_seq_list)):
                if pdb_chain_list[pdb_ind] not in list_chain:
                    seq+=pdb_seq_list[pdb_ind].lower()
                else:
                    for tuple_eval in l_tuple_eval:
                        #t=s(seq,list_chain, ll_score , alig_score , list_seq_chain , id ) id only in post_mortem
                        if pdb_chain_list[pdb_ind]==tuple_eval[-1]:
                            seq+=tuple_eval[0]
                            list_seq_chain.append(tuple_eval[0])
                            ll_score+=tuple_eval[1]
                            alig_score+=tuple_eval[2]
                            if len(tuple_eval)==5:
                                assign_res=len(filter(lambda x: x in string.uppercase, tuple_eval[0]))
                                id+=(tuple_eval[3]* assign_res /100)
                                idres+=tuple_eval[3]
                                count+=assign_res
            #added in 10May2018 to separate chain and resnumber
            diclistresnumb = defaultdict(list)
            for i,s in enumerate(seq):
                if s==s.upper():
                    diclistresnumb[BORGES_MATRIX_listSSpdb[i][0]].append(BORGES_MATRIX_listSSpdb[i][1])
            # print seq
            # exit()
            if len(l_tuple_eval[0])==5:
                id=100*id/count
                seq_list.append((seq, diclistresnumb, ll_score, alig_score, list_seq_chain, id))
            else:
                seq_list.append( (seq,diclistresnumb, ll_score , alig_score , list_seq_chain ) )
    seq_list_clean=[]
    
    return seq_list

##    for chain_eval,list_eval in  dic_possib.items():
##        for tuple_eval in list_eval:
##            seq_eval=tuple_eval[0]
##            seq=''
##            for pdb_ind in range(len(pdb_seq_list)):
##                if pdb_chain_list[pdb_ind]!=chain_eval:
##                    seq+=pdb_seq_list[pdb_ind].lower()
##                else:
##                    seq+=tuple_eval[0]
##            seq_list.append((seq,chain_eval, tuple_eval[1] , tuple_eval[2]))
##    count=0
##    for i in seq_list:
##        count+=1
##        print count
##        print i
##    exit()
##    return seq_list


def return_chain_list (pdb_input): #['A', 'B', 'W']
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdb_input[:-4], pdb_input)
    list_chains=[]
    for chain in struct[0].get_list():
        list_chains.append(chain.get_id())
    return list_chains
#
# def return_list_resnumb_list (pdb_input): #[[1,2,3],[1,2]] 1st list: res of 1st chain, 2nd list: res of 2nd chain
#     p = Bio.PDB.PDBParser(PERMISSIVE=1)
#     struct = p.get_structure(pdb_input[:-4], pdb_input)
#     list_chains=return_chain_list (pdb_input)
#     list_res_numb=[]
#     for chain in list_chains:
#         list_var=[]
#         for res in struct[0][chain].get_list():
#             if res.get_id()[0]==' ':
#                 list_var.append(res.get_id()[1])
#         list_res_numb.append(list_var)
#     return list_res_numb
#
#
def return_dic_resnumb_list (pdb_input): #{'A': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121], 'W': []}
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdb_input[:-4], pdb_input)
    list_chains=return_chain_list (pdb_input)
    dic_res_numb={}
    for chain in list_chains:
        list_var=[]
        for res in struct[0][chain].get_list():
            if res.get_id()[0]==' ':
                list_var.append(res.get_id()[1])
        dic_res_numb[chain]=list_var
    return dic_res_numb


    
def total_numb_res (pdb_input): #it also counts water as residues return integer
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdb_input[:-4], pdb_input)
    residues = struct[0].get_residues()
    rescount=0
    for r in residues:
        rescount+=1
    return rescount

def return_tuple_chain_resnumb_restype (pdb_input):#[[('A', 1, 'S'), ('A', 2, 'A')],[('B', 1, 'S')]]
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdb_input[:-4], pdb_input)
    #list_chains=return_chain_list (pdb_input)
    list_chains=[]
    for chain in struct[0].get_list():
         list_chains.append(chain.get_id())
    list_tuples=[]
    for chain in list_chains:
        list_var=[]
        for res in struct[0][chain].get_list():
            if res.get_id()[0]==' ':
                if res.get_resname() in amino_acid_list_3L:
                    seq=amino_acid_list[amino_acid_list_3L.index(res.get_resname())]
                    resnumb=res.get_id()[1]
                    tuple=(chain,resnumb,seq)
                    list_var.append(tuple)
                else:
                    print ('Compound',res.get_resname() , 'numbered',res.get_id()[1] , 'chained',chain, 'was rejected by RJB_lib.return_tuple_chain_resnumb_restype!')
        list_tuples.append(list_var)
    return list_tuples

def return_restype_list (pdb_input,printe=False): #['SAVQFEVSIIKIAGKSGVEEYGSLGCYCGSGGASRPLLASDACCAVHDCCFGLVSSCDPGAGVYTYAKLLGVATCGLGNPCAVQICECDKVAATCFRKNAAVFDIKRQFKPAKICAEAQPC', '']
    if isinstance(pdb_input, basestring):
        p = Bio.PDB.PDBParser(PERMISSIVE=1)
        struct = p.get_structure(pdb_input[:-4], pdb_input)
    else:
        struct=pdb_input
    list_chains=return_chain_list (pdb_input)
    list_restype=[]
    for chain in list_chains:
        string_var=''
        for res in struct[0][chain].get_list():
            if res.get_id()[0]==' ':
                try:
                    seq=amino_acid_list[amino_acid_list_3L.index(res.get_resname())]
                    string_var+=seq
                except:
                    if printe:
                        print ('Restype',res.get_resname(),'not a residue, therefore it was not considered.')
        list_restype.append(string_var)
    return list_restype
#
# def return_list_of_number_residues_by_chain (pdb_input):
#     l=return_restype_list (pdb_input)
#     ll=[len(i) for i in l]
#     return ll
#
def return_dic_sequence (pdb_input) : #{'A': 'SAVQFEVSIIKIAGKSGVEEYGSLGCYCGSGGASRPLLASDACCAVHDCCFGLVSSCDPGAGVYTYAKLLGVATCGLGNPCAVQICECDKVAATCFRKNAAVFDIKRQFKPAKICAEAQPC', 'W': ''}
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdb_input[:-4], pdb_input)
    list_chains=return_chain_list (pdb_input)
    dic_seq={}
    for chain in list_chains:
        string_var=''
        for res in struct[0][chain].get_list():
            if res.get_id()[0]==' ':
                if res.get_resname() in amino_acid_list_3L:
                    string_var+=amino_acid_list[amino_acid_list_3L.index(res.get_resname())]
                else:
                    print ('Compound',res.get_resname() , 'numbered',res.get_id()[1] , 'chained',chain, 'was rejected by RJB_lib.return_dic_sequence!')
        dic_seq[chain]=string_var
    return dic_seq


def return_list_sequence (pdb_input) : #['SAVQFEVSIIKIAGKSGVEEYGSLGCYCGSGGASRPLLASDACCAVHDCCFGLVSSCDPGAGVYTYAKLLGVATCGLGNPCAVQICECDKVAATCFRKNAAVFDIKRQFKPAKICAEAQPC', '']
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdb_input[:-4], pdb_input)
    list_chains=return_chain_list (pdb_input)
    list_seq=[]
    for chain in list_chains:
        string_var=''
        for res in struct[0][chain].get_list():
            if res.get_id()[0]==' ':
                string_var+=amino_acid_list[amino_acid_list_3L.index(res.get_resname())]
        list_seq.append(string_var)
    return list_seq



def return_chainlist_seqlist_numbres (pdb_input): #(['A', 'W'], [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121], []], 196)
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdb_input[:-4], pdb_input)

    list_chains=return_chain_list (pdb_input)

    list_res_numb=[]
    for chain in list_chains:
        list_var=[]
        for res in struct[0][chain].get_list():
            if res.get_id()[0]==' ':
                list_var.append(res.get_id()[1])
        list_res_numb.append(list_var)

    rescount=total_numb_res (pdb_input)

    return list_chains , list_res_numb , rescount

def return_dic_chain_resnumb_restype (pdb):
    ltup=return_tuple_chain_resnumb_restype (pdb) #[[('A', 1, 'S'), ('A', 2, 'A')],[('B', 1, 'S')]]
    dall={}
    for i in ltup:
        for ii in i:
            ch=ii[0]
            resn=ii[1]
            rest=ii[2]
            try:
                dall[ch][resn]=rest
            except:
                dall[ch]={}
                dall[ch][resn]=rest
    return dall

def remove_partial_res_from_PDB_BioPython ( pdb_input, pdb_output ):
    p = Bio.PDB.PDBParser(PERMISSIVE=1)
    struct = p.get_structure(pdb_input[:-4], pdb_input)
    residues = struct[0].get_residues()
    verify=False
    while verify==False:
        for res in residues:
            if len(res)>=4:
                verify=True
    if verify:
        class FullResSelect(Bio.PDB.Select):
            def accept_residue(self, residue):
                if len(residue)>=4:
                    return True
                else:
                    return False
        io = Bio.PDB.PDBIO()
        io.set_structure(struct)
        io.save( pdb_output , FullResSelect() )
        print ('Script RJB_lib.remove_partial_res_from_PDB_BioPython corrected pdb',pdb_input,'and named to',pdb_output)
        add_previous_ATOM_lines_to_pdb ( pdb_input, pdb_output , pdb_output )
        print ('Script RJB_lib.add_previous_ATOM_lines_to_pdb added header into same name file.')
    else:
        print ('\n\nPDB without partial residues.\n\n')
    return verify

def SortResUniChainAfromSHREDDER ( pdb_input, pdb_output ):
    with open(pdb_input) as f: f2=f.readlines()
    with open(pdb_output,'w') as f:
        listatoms=[]
        check=True
        for l in f2:
            if not l.startswith('ATOM') and not l.startswith('HETATM'):
                if check: f.write(l)
            else:
                check=False
                listatoms.append(l)
        listatoms=sorted(listatoms,key=lambda x: int(x[22:26]) )
        #sorted(student_tuples, key=lambda student: student[2])
        #(key = (lambda item: item[-1]), reverse = False)
        for l in listatoms:
            ll=l[:21]+'A'+l[22:]
            f.write(ll)


def add_previous_ATOM_lines_to_pdb (pdb_input_header, pdb_input_atoms , pdb_output):
    input_header=open(pdb_input_header,'r')
    header=input_header.read()
    input_header.close()
    input_atoms=open(pdb_input_atoms,'r')
    atoms=input_atoms.read()
    input_atoms.close()
    output=open(pdb_output,'w')
    output.write(header[:header.index('\nATOM')]+'\n')
    output.write(atoms)
    output.close()

def print_list_dics_sspred_conf_in_nice_way (list_dics_sspred_conf):
    for i in range(len(list_dics_sspred_conf)):
        print ('Fragment',i+1,'of',list_dics_sspred_conf[i]['sec_str_pred'][0],'with length',len(list_dics_sspred_conf[i]['sequence']))
        print (list_dics_sspred_conf[i]['initial_res'],list_dics_sspred_conf[i]['sequence'],list_dics_sspred_conf[i]['final_res'],'\n')
        


def given_2_superposed_pdbs_return_list_res_match ( pdb_match , pdb_fit , dist=1.0):#method='one') : #method may be all (get all atoms)
    #Script written 10Jun2016 for SLIDER_phasing to calculate correspondence of res of two structures to be calculated post_morten analysis
    #method=X if distance is above 1.0A between trace and reference, an 'X' is given
    #method=write write correspondent residues
    ProteinPDB = Bio.PDB.PDBParser()
    trace = ProteinPDB.get_structure ( pdb_match[:-4] , pdb_match )
    ent = ProteinPDB.get_structure ( pdb_fit[:-4] , pdb_fit )

    #list_CA_trace=[]
    list_CA_trace_related_ent_BioPDB = []
    for atom in trace.get_atoms():
        if atom.get_name()=='CA':
            a_tr = atom
            r_tr = a_tr.get_parent()
            r_tr_type = amino_acid_list[amino_acid_list_3L.index(r_tr.get_resname())]
            c_tr = r_tr.get_parent()
            #list_CA_trace.append(atom)
            list_CA_trace_related_ent_BioPDB.append([(c_tr.get_id(), r_tr.get_id()[1], r_tr_type, a_tr)])

    list_CA_ent=[]
    for atom in ent.get_atoms():
        if atom.get_name()=='CA':
            list_CA_ent.append(atom)
    ns=Bio.PDB.NeighborSearch(list_CA_ent)
    #
    # listforcomb=[]
    # for match in list_CA_trace_related_ent_BioPDB:
    #     listvar=[]
    #     coord = CAtrace.get_coord()
    #     pair = ns.search(coord, i, "A")
    #     for CAent in pair:
    #         listvar.append(CAent)
    #     listforcomb.append(listvar)



    table_remove=[]
    for i in numpy.arange(0.5,dist+0.1,0.25):
        # print i
        #count = 0
        #for CA in list_CA_trace:
        for listpair in list_CA_trace_related_ent_BioPDB:
            if len(listpair)==1:
                # print 'listpair1',listpair
                CAtrace=listpair[0][-1]
                coord= CAtrace.get_coord()
                pair = ns.search(coord, i, "A")
                # print 'pair',pair,len(pair)
                if len(pair)>0:
                    listvar=[]
                    for CAent in pair:
                        if CAent not in table_remove: listvar.append( (CAent,CAtrace-CAent) )

                    if len(listvar)>0:
                        # print 'listpair2',listvar
                        listvar.sort(key=(lambda item: item[-1]), reverse=False)
                        # print 'listpair2',listvar
                        a_ent= listvar[0][0]
                        r_ent=a_ent.get_parent()
                        r_ent_type=amino_acid_list[amino_acid_list_3L.index(r_ent.get_resname())]
                        c_ent=r_ent.get_parent()
                        distvar=listvar[0][1]

                        table_remove.append(a_ent)
                        listpair.append((c_ent.get_id(),r_ent.get_id()[1],r_ent_type, a_ent))
                        listpair.append(distvar)

    identicalresiduesn = 0
    for i in list_CA_trace_related_ent_BioPDB:
        #print i
        if len(i) == 3 and i[0][2] == i[1][2]: identicalresiduesn += 1
    percentidentity = 100 * float(identicalresiduesn) / len(list_CA_trace_related_ent_BioPDB)

    return list_CA_trace_related_ent_BioPDB , identicalresiduesn , percentidentity
    #output is like list_CA_trace_related_ent_BioPDB[0]=[('A', 2, 'N', <Atom CA>)]       #here no match was found
                   #list_CA_trace_related_ent_BioPDB[1]=[('A', 3, 'T', <Atom CA>), ('A', 3, 'K', <Atom CA>), 1.2654423]


    # print 'Algorithm made to obtain minimum independent distance:'
    #
    # listdist = []
    # for i in list_CA_trace_related_ent_BioPDB:
    #     print i
    #     if len(i) > 1:
    #         listdist.append(i[-1])
    #
    # rmsd = calculateRMSDfromdistances(listdist)
    # rmsd = '%.1f' % (rmsd)
    # print rmsd,'\n\n\n'
    #
    # return list_CA_trace_related_ent_BioPDB , rmsd
    #
    # listbyfrag=[]
    # chprev   = list_CA_trace_related_ent_BioPDB[0][0][0]
    # resnprev = list_CA_trace_related_ent_BioPDB[0][0][1]-1
    # listvar=[]
    # for match in list_CA_trace_related_ent_BioPDB:
    #     ch  =match[0][0]
    #     resn=match[0][1]
    #     if ch==chprev and resn-1==resnprev:
    #         listvar.append(match)
    #     else:
    #         listbyfrag.append(listvar)
    #         listvar=[]
    #     chprev=ch
    #     resnprev=resn
    # listbyfrag.append(listvar)
    #
    # listbyfrag.sort(key=(lambda item: len(item)), reverse=True)
    #
    # list_CA_trace_related_ent_BioPDB2=[]
    #
    # for frag in listbyfrag:
    #     lresnmatch= [x[1] for x in frag]
    #     minm=numpy.min(lresnmatch)
    #     maxm = numpy.max(lresnmatch)
    #     listvar=[]
    #     for i in range(minm-5,maxm+5):
    #         countt=i
    #         count2=0
    #         for aa in frag:
    #             resn=aa[1][1]
    #             if countt==resn: count2+=1
    #         listvar.append(count2)
    #     best=minm-5+listvar.index( max(listvar) )
    #     for i, aa in enumerate(frag):
    #         tupletrace=aa[0]
    #         atomCAentnew = ent[0]['A'][100]['CA']
    #         list_CA_trace_related_ent_BioPDB2
    #
    #
    # return list_CA_trace_related_ent_BioPDB



def GivenListMatchCA2PDBsReturnCorrespondentResString (list_CA_trace_related_ent):
    seq=''
    for item in list_CA_trace_related_ent:
        try:
            seq+=item[1][2]
        except:
            seq+='X'
    return seq #seq='XXKIVLF'

# def given_list2pdbs_superposed_return_correspondent_res_string (pdb_match,pdb_fit,dist=1.0): OBSOLETE NOW USE FIRST given_2_superposed_pdbs_return_list_res_match THEN GivenListMatchCA2PDBsReturnCorrespondentResString

# def given_2pdbs_superposed_return_correspondent_res_list_string (pdb_match,pdb_fit): # OBSOLETE NOW USE FIRST given_2_superposed_pdbs_return_list_res_match THEN GivenListListMatchCA2PDBsReturnListSequenceByChain
#     listt=given_2_superposed_pdbs_return_list_res_match ( pdb_match , pdb_fit )
#     list_seq=given_list_pair_res_return_list_sequence_by_chain (listt)
#     return list_seq

# def given_2pdbs_superposed_return_correspondent_res_dic_string (pdb_match,pdb_fit,dist=1.0): # OBSOLETE NOW USE FIRST given_2_superposed_pdbs_return_list_res_match THEN GivenListListMatchCA2PDBsReturnDicSequenceByChain
#     listt=given_2_superposed_pdbs_return_list_res_match ( pdb_match , pdb_fit , dist)
#     dic_seq=given_list_pair_res_return_dic_sequence_by_chain (listt)
#     return dic_seq #seq['A']='XXKIVLF'


def GivenListListMatchCA2PDBsReturnListSequenceByChain(list_CA_trace_related_ent):
# def given_list_pair_res_return_list_sequence_by_chain (list_CA_trace_related_ent): OBSOLETE
    verify=list_CA_trace_related_ent[0][0][0]
    seq=''
    list_seq=[]
    for item in list_CA_trace_related_ent:
        if not item[0][0]==verify:
            list_seq.append(seq)
            seq=''
            verify=item[0][0]
        try:
            seq+=item[1][2]
        except:
            seq+='X'
    list_seq.append(seq)
    print (list_seq)
    return list_seq


def GivenListListMatchCA2PDBsReturnDicSequenceByChain(list_CA_trace_related_ent):
# def given_list_pair_res_return_dic_sequence_by_chain (list_CA_trace_related_ent): # OBSOLETE
    dic_seq={}
    for item in list_CA_trace_related_ent:
        try:
            dic_seq[item[0][0]]
        except:
            dic_seq[item[0][0]]=''
        try:
            dic_seq[item[0][0]]+=item[1][2]
        except:
            dic_seq[item[0][0]]+='X'
    return dic_seq #dic_seq['A']='XXKIVLF'

def GivenListListMatchCA2PDBsReturnDicChResNResT(list_CA_trace_related_ent):
    DicResNResT={}
    for item in list_CA_trace_related_ent:
        ch  =item[0][0]
        resn=item[0][1]
        if len(item)>1: restrue=item[1][2]
        else:           restrue='X'
        if ch not in DicResNResT: DicResNResT[ch]={resn:restrue}
        else:                     DicResNResT[ch][resn]=restrue
    return DicResNResT

def GivenDicChResNResTReturnDicChResNRangeSeq(DicChResNResT,DicChResNResSS=False):
    DicChResNRangeSeq={}
    for ch,DicResNResT in DicChResNResT.items():
        DicChResNRangeSeq[ch]={}
        prevRN = sorted(DicResNResT)[0]
        resni  = str(prevRN) + '-'
        seq    = DicResNResT[prevRN]
        if DicChResNResSS!=False: prevSS=DicChResNResSS[ch][prevRN]

        for resn in sorted(DicResNResT.keys())[1:]:
            if prevRN == resn-1 and (DicChResNResSS==False or DicChResNResSS[ch][resn]==prevSS): #condition to increase segment
                seq   += DicResNResT[resn]
            else:                #condition to save previous segment and start new one
                DicChResNRangeSeq[ch][resni+str(prevRN)]=seq
                resni=str(resn)+'-'
                seq=DicResNResT[resn]
            prevRN = resn
            if DicChResNResSS!=False: prevSS = DicChResNResSS[ch][resn]
        # condition to save last and forgotten segment
        DicChResNRangeSeq[ch][resni + str(prevRN)] = seq
    # for ch in TrueSeqDicResNRangeSeq:
    #     print ch
    #     for resn in sorted(TrueSeqDicResNRangeSeq[ch]):
    #         print resn, TrueSeqDicResNRangeSeq[ch][resn]
    return DicChResNRangeSeq

def GivenPDBlistChResNCAReturnlistChResNRangeCA (listChResNCA): #input = [['A',1],['A',2],['A',3],['A',10],['A',11],['A',12],['A',13]]
    listChResNRangeCA=[]                                        #output= [['A','1-3'],            ['A','10-13']                      ]
    prevresn=listChResNCA[0][-1]
    prevch = listChResNCA[0][0]
    listvar=[prevch,[prevresn]]
    for ch,resn in listChResNCA[1:]:
        if ch==prevch and resn-1==prevresn: listvar[-1].append(resn)
        else:
            listChResNRangeCA.append([listvar[0],str(listvar[-1][0])+'-'+str(listvar[-1][-1])])
            listvar = [ch, [resn]]
        prevresn=resn
        prevch  =ch
    listChResNRangeCA.append([listvar[0], str(listvar[-1][0]) + '-' + str(listvar[-1][-1])])
    return listChResNRangeCA


def GivenListMatchCA2PDBsReturnNumberIdenticalRes (list_CA_trace_related_ent):
    NumberIdenticalRes=0
    for match in list_CA_trace_related_ent:
        if len(match)>1 and match[0][2]==match[1][2]: NumberIdenticalRes+=1
    return NumberIdenticalRes


#def phenix_refine (pdb_file,mtz_file,output_file,F,SIGF,parameters):
    #print 'to be implemented'
    #template extracted from: https://www.phenix-online.org/documentation/reference/refinement.html#running-phenix-refine
    #% phenix.refine <pdb-file(s)> <reflection-file(s)> <monomer-library-file(s)>  <parameter-keyword(s)> <parameter-file(s)>
    #'phenix.refine '+pdb_file+' '+mtz_file+' refinement.input.xray_data.labels="'+F+','+SIGF+'" strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body ramachandran_restraints=True' refinement.output.prefix='+output_file+' main.number_of_macro_cycles=20 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=True optimize_adp_weight=True  secondary_structure.enabled=True nproc=1'
    #phenix.refine n1-polyA.pdb 1YZF.mtz      refinement.input.xray_data.labels="FOBS,SIGFOBS"   strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body ramachandran_restraints=True  refinement.output.prefix=n1-polyA_phenix main.number_of_macro_cycles=20 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=True optimize_adp_weight=True  secondary_structure.enabled=True nproc=1


def sh_refine_phenix_multiprocesses ( pdb_input_file , mtz_input_file , pdb_output_file , phenixrefine_path , F, SIGF, sh_file , distribute_computing , options=False):

    file_sh=open(sh_file,'w')
##    if distribute_computing=='multiprocessing':
##        file_sh.write('#!/bin/sh'+ '\n')
    if distribute_computing=='local_grid':
        file_sh.write('#!/bin/tcsh'+ '\n')
        file_sh.write('setenv PATH /usr/local/bin\nsetenv PATH /usr/bin:$PATH\nsetenv PATH /bin:$PATH\n')  # suggested by MAX
        #file_sh.write('source /xtal/xtalsetup'+ '\n')
        file_sh.write('hostname > '+sh_file.split('/')[-1]+ '_host \n')
        file_sh.write('source /xtal/phenix/phenix/phenix_env.csh\n')
    xyzin=pdb_input_file
    outf=pdb_output_file[:-4].split('/')[-1]
    hklin=mtz_input_file
    #hklout=xyzout[:-4]+'.mtz'
    #phenix_output_log = xyzout[:-4]+'.log'
    if distribute_computing=='local_grid':
        xyzin =xyzin.split('/')[-1]
        pdb_output_file = pdb_output_file.split('/')[-1]
        hklin =hklin.split('/')[-1]
        #hklout=hklout.split('/')[-1]
        #phenix_output_log = phenix_output_log.split('/')[-1]

    #if not options: options = 'strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body main.number_of_macro_cycles=5 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=False optimize_adp_weight=False nproc=1 export_final_f_model=true'

    # fewer restrictions   5LCY#
    if not options:                          options = 'strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body                              main.number_of_macro_cycles=5 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=False optimize_adp_weight=False                                   nproc=1 export_final_f_model=true '
    #with NCS             5LCY/2#
    #if not options:                          options = 'strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body                              main.number_of_macro_cycles=5 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=False optimize_adp_weight=False                                   nproc=1 export_final_f_model=true ncs_search.enabled=True '

    #various restrictions 1YZF#if not options: options= 'strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body ramachandran_restraints=True main.number_of_macro_cycles=5 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=False optimize_adp_weight=False  secondary_structure.enabled=True nproc=1 export_final_f_model=true'
    #1YZF#                     if not options: options= 'strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body ramachandran_restraints=True main.number_of_macro_cycles=5 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=False optimize_adp_weight=False  secondary_structure.enabled=True nproc=1 export_final_f_model=true simulated_annealing=true'
    #if not options: options = 'strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body ramachandran_restraints=True main.number_of_macro_cycles=5 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=False optimize_adp_weight=False  secondary_structure.enabled=True nproc=1 export_final_f_model=true simulated_annealing=true'

    file_sh.write(phenixrefine_path+' '+xyzin+' '+hklin+' refinement.input.xray_data.labels="'+F+','+SIGF+'" prefix='+outf+' '+options+'\n')
    file_sh.write('mv ' + outf + '_001.mtz '+pdb_output_file[:-4]+'.mtz\n')
    file_sh.write('mv '+outf+'_001_f_model.mtz '+pdb_output_file[:-4]+'_FOM.mtz\n')
    file_sh.write('mv ' + outf + '_001.log ' + pdb_output_file[:-4] + '.log\n')
    file_sh.write('mv ' + outf + '_001.pdb ' + pdb_output_file+'\n')
    file_sh.close()

# def phenix_refine_multiprocessing (pdb_input,mtz_input,output,F,SF,options=False): #type is a string
#     type could be < None , strategy=rigid_body , 'optimize_xyz_weight=true optimize_adp_weight=true' , 'strategy=tls tls.find_automatically=True' , ncs_search.enabled=True >
#    file_name=dictionary_buster['file'].split('.pd')[0]
#    folder_output=dictionary_buster['folder']+'/BUSTER_temp_'+file_name
#    xyzin=dictionary_buster['folder']+'/'+dictionary_buster['file']
#    os.system( 'mkdir ' + folder_output )
#    buster_options='-noWAT -nbig 10 -RB -nthread 1'
#     #-w 50 AdjustXrayWeightAutomatically="no" -autoncs
    #'phenix.refine '+pdb_file+' '+mtz_file+' refinement.input.xray_data.labels="'+F+','+SIGF+'" refinement.output.prefix='+output_file+' strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body ramachandran_restraints=True main.number_of_macro_cycles=20 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=True optimize_adp_weight=True  secondary_structure.enabled=True nproc=1'
    #phenix.refine n1-polyA.pdb 1YZF.mtz      refinement.input.xray_data.labels="FOBS,SIGFOBS"  refinement.output.prefix=example          strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body ramachandran_restraints=True main.number_of_macro_cycles=20 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=True optimize_adp_weight=True  secondary_structure.enabled=True nproc=1
    # if not options: options='strategy=individual_sites+individual_adp+individual_sites_real_space+rigid_body ramachandran_restraints=True main.number_of_macro_cycles=20 write_eff_file=false write_geo_file=false  write_def_file=false optimize_xyz_weight=True optimize_adp_weight=True  secondary_structure.enabled=True nproc=1'
    # '=+' export_final_f_model=true output.write_eff_file=False output.write_geo_file=False output.write_def_file=False refinement.main.number_of_macro_cycles=10'
    # p = subprocess.Popen(['phenix.refine',pdb_input,mtz_input,'refinement.input.xray_data.labels="'+F+','+SIGF,'prefix='+output]+options.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # out, err = p.communicate()
    #export_final_f_model=true
#    file=open(folder_output+'/'+file_name+'.log','w')
#    file.write(out)
#    file.close()
    
    #phenix.refine mtz_input pdb_input prefix=output_name
    #< None , strategy=rigid_body , 'optimize_xyz_weight=true optimize_adp_weight=true' , 'strategy=tls tls.find_automatically=True' , ncs_search.enabled=True >
    #export_final_f_model=true output.write_eff_file=False output.write_geo_file=False output.write_def_file=False

def GetPhenixRefineLLG (log):
    with open(log) as f:
        fl=f.read()
        #
        indnormllg=fl.rindex('normalized target function (ml) (work):')+len('normalized target function (ml) (work):')
        strnormllg=fl[indnormllg:indnormllg+28]
        normllg=float(strnormllg)
        #
        indllgwork=fl.rindex('target function (ml) not normalized (work):')+len('target function (ml) not normalized (work):')
        strllgwork=fl[indllgwork:indllgwork+28]
        llgwork=float(strllgwork)
        #
        indllgfree=fl.rindex('target function (ml) not normalized (free):')+len('target function (ml) not normalized (free):')
        strllgfree=fl[indllgfree:indllgfree+28]
        llgfree=float(strllgfree)
    return normllg , llgwork , llgfree

def generate_matrix_SLIDER_alignment (includeX=False):
    matrix_align={}
    for ind1,aa1 in enumerate(amino_acid_list):
        for ind2,aa2 in enumerate(amino_acid_list):
            if aa1==aa2:
                matrix_align[(aa1.lower(),aa2)]=4
                matrix_align[(aa1,aa2)]=6
            else:
                matrix_align[(aa1.lower(),aa2)]=-1
                matrix_align[(aa1,aa2)]=-4
    matrix_align[('a','A')]=0
    if includeX:
        for aa1 in amino_acid_list:
            for k in ['X']:#,'-']:
                matrix_align[(aa1, k)] = 0
                matrix_align[(k, aa1)] = 0
    #from Bio.SubsMat import MatrixInfo
    #matrix_align=MatrixInfo.ident
    return matrix_align

    
def alignment_score (dic_possib , sequence , post_mortem , true_seq_dic_string ): #dic_possib key:chain letter=list of tuples of seq,my_score)
    matrix_align=generate_matrix_SLIDER_alignment ()
    from Bio import pairwise2
    new_dic={}
    for chain in sorted (dic_possib):
        element=dic_possib[chain]
        print ('Aligning and scoring chain',chain,'('+str(len(dic_possib[chain])),'sequences)')
        list_var=[]
#_align(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, penalize_extend_when_opening, penalize_end_gaps, align_globally, gap_char, force_generic, score_only, one_alignment_only)
#pairwise_tuple=pairwise2.align.localdd (sequence, tuple_elem[0] ,matrix_align,-15,-5) #(1st seq,2nd seq,penalization for opening gap,penalization extending gap)
        for tuple_ind , tuple_elem in enumerate (element):
            #seq=tuple_elem[0].replace('a','X')
            #seq=seq.upper()
            seq=tuple_elem[0]
            #pairwise_tuple=pairwise2.align.localdd (sequence, seq ,matrix_align,-50,-25,-15,-5) #(1st seq,2nd seq,penalization for opening gap,penalization extending gap)
            pairwise_tuple=pairwise2.align.localds (sequence, seq ,matrix_align,-2,-1) #(1st seq,2nd seq,penalization for opening gap,penalization extending gap)
            align_score=float(pairwise_tuple[0][2])
            #print(pairwise2.format_alignment(*pairwise_tuple[0])) #print best alignment
##            for a in pairwise_tuple:
##                print(pairwise2.format_alignment(*a))
            if not post_mortem:
                list_var.append( (tuple_elem[0],float(tuple_elem[1]),float(align_score)) )
            else:
                seq1=tuple_elem[0]
                seq2=true_seq_dic_string[chain]
##                print seq1
##                print seq2
                seq11=''
                seq22=''
                for ind in range(len(seq1)):
                    if seq1[ind].upper()==seq1[ind]:
                        seq11+=seq1[ind]
                        seq22+=seq2[ind]
                #ident,simil=given_two_sequences_return_identity_similarity(tuple_elem[0], true_seq_list_string[chain_ind] )
                # print 'seq11',seq11
                # print 'seq22',seq22
                # exit()
                ident,simil=given_two_sequences_return_identity_similarity(seq11, seq22 )
                list_var.append( (tuple_elem[0],float(tuple_elem[1]),float(align_score),ident) )
        new_dic[chain]=list_var
    return new_dic

def AlignTwoSequences (seq1,seq2,matrix_align=False,printt=False):#,MatrixAlign=matlist.blosum62): # this was written in Feb27 2018 and it has proven really bad to find remoto homologous sequence pairing
                                                                   # clustalomega was effective in 1YZF sequence pairing and was replaced by it, its file should be given by user
    #options available on:
    #matrix_align =  MatrixAlign
    #matrix_align = MatList.blosum62
    #matrix_align = MatList.blosum30
    #matrix_align = MatList.gonnet
    #matrix_align = MatList.ident
    if matrix_align==False: matrix_align = generate_matrix_SLIDER_alignment(includeX=True)
    #print 'Aligning and scoring:'
    #print seq1,seq2
    # _align(sequenceA, sequenceB, match_fn, gap_A_fn, gap_B_fn, penalize_extend_when_opening, penalize_end_gaps, align_globally, gap_char, force_generic, score_only, one_alignment_only)
    # pairwise_tuple=pairwise2.align.localdd (sequence, tuple_elem[0] ,matrix_align,-15,-5) #(1st seq,2nd seq,penalization for opening gap,penalization extending gap)
    # pairwise_tuple=pairwise2.align.localdd (sequence, seq ,matrix_align,-50,-25,-15,-5) #(1st seq,2nd seq,penalization for opening gap,penalization extending gap)
    #pairwise_tuple = pairwise2.align.localds(seq1, seq2, matrix_align, -10, -1)  # (1st seq,2nd seq,penalization for opening gap,penalization extending gap)
    pairwise_tuple = pairwise2.align.localds(seq1, seq2, matrix_align, -2, -1)  # (1st seq,2nd seq,penalization for opening gap,penalization extending gap)
    align_score = float(pairwise_tuple[0][2])
    if printt:
        for a in pairwise_tuple:
            print(pairwise2.format_alignment(*a))
    # exit()
    return pairwise_tuple[0] , align_score


def write_refmac_TMP_default_file ( refmacTMP_filename , ncycles , additional_line , mtz_f=False , mtz_sf=False , mtz_free=False ): # if you would like to add additional options not written by default, give a string here, otherwise just use empty string ''
    refmacTMP=open(refmacTMP_filename,'w')
    refmacTMP.write( 'make check NONE\nmake -\n    hydrogen ALL -\n    hout NO -\n    peptide NO -\n    cispeptide YES -\n    ssbridge YES -\n    symmetry YES -\n    sugar YES -\n    connectivity NO -\n    link NO\nrefi -\n    type REST -\n    resi MLKF -\n    meth CGMAT -\n    bref ISOT\n')
    refmacTMP.write( additional_line )
    refmacTMP.write( '\nncyc ' + ncycles + '\nscal -\n    type SIMP -\n    LSSC -\n    ANISO -\n    EXPE\nsolvent YES\nweight -\n    AUTO\nmonitor MEDIUM -\n    torsion 10.0 -\n    distance 10.0 -\n    angle 10.0 -\n    plane 10.0 -\n    chiral 10.0 -\n    bfactor 10.0 -\n    bsphere 10.0 -\n    rbond 10.0 -\n    ncsr 10.0\n')
    refmacTMP.write( 'labin  FP=')
    if mtz_f==False:
        refmacTMP.write( 'F')
    else:
        refmacTMP.write( mtz_f )
    refmacTMP.write(' SIGFP=')
    if mtz_sf==False:
        refmacTMP.write( 'SIGF')
    else:
        refmacTMP.write( mtz_sf )
    refmacTMP.write(' FREE=')
    if mtz_free==False:
        refmacTMP.write( 'FreeR_flag')
    else:
        refmacTMP.write(mtz_free)
    refmacTMP.write('\n')
    #refmacTMP.write('labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM\nPNAME unknown\nDNAME unknown101\nRSIZE 80\nEXTERNAL WEIGHT SCALE 10.0\nEXTERNAL USE MAIN\nEXTERNAL DMAX 4.2\nEND' )
    refmacTMP.write( 'labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM\n')
    #refmacTMP.write( 'PNAME unknown\nDNAME unknown101\nRSIZE 80\nEXTERNAL WEIGHT SCALE 10.0\nEXTERNAL USE MAIN\nEXTERNAL DMAX 4.2\nEND')
    refmacTMP.write('RIDG DIST SIGM 0.02\nEND')
    refmacTMP.close()

def retrieve_energies_from_Scwrl4_log_file ( Scwrl4_log_file ) : #extracted from SLIDER script
    energy_file=open(Scwrl4_log_file,'r')
    energy_list=energy_file.readlines()
    energy_file.close()
    for energ in energy_list:
        if energ.startswith('Total minimal energy of the graph'):
#energ.startswith('Energy of this cluster') or
            return float( energ[energ.index('=')+2:energ.index('\n')] )



def given_list_dic_files_add_factors_calculate_FOM ( list_dics_files , PDB_input_key , inicial_key , refinement_program ):
    list_dics_files=given_list_dic_files_add_factors ( list_dics_files , PDB_input_key , inicial_key , refinement_program )
    dic_FOM=given_list_dic_files_generate_dic_FOM_by_chains ( list_dics_files , inicial_key , refinement_program )
    list_dics_files=given_list_dic_files_and_dic_FOM_calculate_all_FOM ( list_dics_files , dic_FOM , inicial_key , refinement_program )
    return list_dics_files


def given_list_dic_files_add_factors ( list_dics_files , PDB_input_key , inicial_key , refinement_program ):
    for dic in list_dics_files:
        PDB_input_file=dic[PDB_input_key]
        if refinement_program=='refmac' or refinement_program=='buster':
            dic[inicial_key+'_R'],dic[inicial_key+'_Rfree']=retrieve_Rfactor_from_PDB_file ( PDB_input_file )
            dic[inicial_key+'_CC_ampl'],dic[inicial_key+'_CC_ampl_free']=retrieve_CC_ampl_from_PDB_file ( PDB_input_file )
        if refinement_program=='buster':
            #dic[inicial_key+'_CCmc'],dic[inicial_key+'_CCsc']=retrieve_busterCC_from_HTML_file ( PDB_input_file[:-4]+'.html' )
            dic[inicial_key + '_CCmc'], dic[inicial_key + '_CCsc'] = retrieve_busterCC_from_log_file(PDB_input_file[:-4] + '_RSCC.log')
    return list_dics_files

def obtainRfactorsBUSTERMapOnly (pdbfile):
    Rfactor=[]
    Rfree = []
    with open(pdbfile, 'r') as f:
        f2=f.readlines()
        f.close()
        skipp=False
        for l in f2:
            if l.startswith('REMARK best refinement for ') and 'R/Rfree' in l and not skipp:
                # print 'Rfactor',l[-14:-8]
                # print 'Rfree  ',l[-7:-1]
                # print l,f2.index(l)
                Rfactor.append(float(l[-14:-8]))
                Rfree.append(float(l[-7:-1]))
                skipp=True
    if len(Rfactor)==1 and len(Rfree)==1: return Rfactor[0],Rfree[0]
    else:
        print ('More than one line meeting criteria in file',pdbfile,'function RJB_lib.obtainRfactorsBUSTERMapOnly')
        print (Rfactor,Rfree)
        exit()


def given_list_dic_files_generate_dic_FOM_by_chains ( list_dics_files , inicial_key , refinement_program ): #should be refine1 or refine2
##    dic_FOM=defaultdict(dict)
    dic_FOM={}
    for dic in list_dics_files:
        ch='_'.join(dic['chain_eval'])
        if refinement_program=='refmac' or refinement_program=='buster':
            try:
                dic_FOM[ch]
            except:
                dic_FOM[ch]={'Rs':[],'Rsfree':[],'CCs_ampl':[],'CCs_ampl_free':[]}
            dic_FOM[ch]['Rs'].append(dic[inicial_key+'_R'])
            dic_FOM[ch]['Rsfree'].append(dic[inicial_key+'_Rfree'])
            dic_FOM[ch]['CCs_ampl'].append(dic[inicial_key+'_CC_ampl'])
            dic_FOM[ch]['CCs_ampl_free'].append(dic[inicial_key+'_CC_ampl_free'])
        if refinement_program=='buster':
            try:
                dic_FOM[ch]['minCCmc']
            except:
                dic_FOM[ch]['CCsmc']=[]
                dic_FOM[ch]['CCssc']=[]
            dic_FOM[ch]['CCsmc'].append(dic[inicial_key+'_CCmc'])
            dic_FOM[ch]['CCssc'].append(dic[inicial_key+'_CCsc'])
    for chain,dic in dic_FOM.items():
        dic_FOM[chain]['maxR']=max(dic['Rs'])
        dic_FOM[chain]['maxRf']=max(dic['Rsfree'])
        dic_FOM[chain]['minCC']=min(dic['CCs_ampl'])
        dic_FOM[chain]['minCCf']=min(dic['CCs_ampl_free'])
        if refinement_program=='buster':
            dic_FOM[chain]['minCCmc']=min(dic['CCsmc'])
            dic_FOM[chain]['minCCsc']=min(dic['CCssc'])
    return dic_FOM


def given_list_dic_files_add_factors_calculate_myFOM ( list_dics_files , inicial_key ):
    dic_FOM=given_list_dic_files_generate_dic_myFOM_by_chains ( list_dics_files , inicial_key  )
    list_dics_files=given_list_dic_files_and_dic_myFOM_calculate_all_myFOM ( list_dics_files , dic_FOM , inicial_key )
    return list_dics_files

def given_list_dic_files_generate_dic_myFOM_by_chains ( list_dics_files , inicial_key ): #should be refine1 or refine2
##    dic_FOM=defaultdict(dict)
    dic_FOM={}
    for dic in list_dics_files:
            ch='_'.join(dic['chain_eval'])
            try:
                dic_FOM[ch]
            except:
                dic_FOM[ch]={'Rs':[],'Rsfree':[],'CCsmc':[],'CCssc':[]}
            dic_FOM[ch]['Rs'].append(dic[inicial_key+'_R'])
            dic_FOM[ch]['Rsfree'].append(dic[inicial_key+'_Rfree'])
            dic_FOM[ch]['CCsmc'].append(dic[inicial_key+'_CCmc_chain_mean'])
            dic_FOM[ch]['CCssc'].append(dic[inicial_key+'_CCsc_chain_mean'])
    for chain,dic in dic_FOM.items():
        dic_FOM[chain]['maxR']=max(dic['Rs'])
        dic_FOM[chain]['maxRf']=max(dic['Rsfree'])
        dic_FOM[chain]['minCCmc']=min(dic['CCsmc'])
        dic_FOM[chain]['minCCsc']=min(dic['CCssc'])
    return dic_FOM

def given_list_dic_files_and_dic_myFOM_calculate_all_myFOM ( list_dics_files , dic_FOM , inicial_key ):
    for dic in list_dics_files:
        chain='_'.join(dic['chain_eval'])
        maxR=dic_FOM[chain]['maxR']
        maxRfree=dic_FOM[chain]['maxRf']
        minCCmc=dic_FOM[chain]['minCCmc']
        minCCsc=dic_FOM[chain]['minCCsc']
        myFOM=(maxR+maxRfree-minCCmc-minCCsc-dic[inicial_key+'_R']-dic[inicial_key+'_Rfree']+dic[inicial_key+'_CCmc_chain_mean']+dic[inicial_key+'_CCsc_chain_mean'])*10
        dic[inicial_key+'_myFOM']=myFOM
    return list_dics_files

##    if refinement_program=='refmac' or refinement_program=='buster':
##        r1_maxR = max(map(lambda x: x['refine1_R'], list_dics_all_files.values()))
##        r1_maxRfree = max(map(lambda x: x['refine1_Rfree'], list_dics_all_files.values()))
##        r1_minCC_ampl = min(map(lambda x: x['refine1_CC_ampl'], list_dics_all_files.values()))
##        r1_minCC_ampl_free = min(map(lambda x: x['refine1_CC_ampl_free'], list_dics_all_files.values()))
##    if refinement_program=='buster':
##        r1_minCCmc = min(map(lambda x: x['refine1_CCmc'], dictio_buster.values()))
##        r1_minCCsc = min(map(lambda x: x['refine1_CCsc'], dictio_buster.values()))





def given_list_dic_files_and_dic_FOM_calculate_all_FOM ( list_dics_files , dic_FOM , inicial_key , refinement_program ):
    for dic in list_dics_files:
        chain='_'.join(dic['chain_eval'])
        maxR=dic_FOM[chain]['maxR']
        maxRfree=dic_FOM[chain]['maxRf']
        minCC_ampl=dic_FOM[chain]['minCC']
        minCC_ampl_free=dic_FOM[chain]['minCCf']
        FOM=(maxR+maxRfree-minCC_ampl-minCC_ampl_free-dic[inicial_key+'_R']-dic[inicial_key+'_Rfree']+dic[inicial_key+'_CC_ampl']+dic[inicial_key+'_CC_ampl_free'])*10
        dic[inicial_key+'_FOM']=FOM
        if refinement_program=='buster':
            minCCmc=dic_FOM[chain]['minCCmc']
            minCCsc=dic_FOM[chain]['minCCsc']
            FOM_Katherin=(maxR+maxRfree-minCCmc-minCCsc-dic[inicial_key+'_R']-dic[inicial_key+'_Rfree']+dic[inicial_key+'_CCmc']+dic[inicial_key+'_CCsc'])*10
            dic[inicial_key+'_kFOM']=FOM_Katherin
    return list_dics_files

def write_table_from_list_dict_and_keys_sorted_by (list_item,dic_all,key_sorted,reverse_state,output_filename): #not tested, generated 4/4/16 to write superimposition tables monomer B with monB
    output_file=open(output_filename,'w')
    for item in list_item:
        if item.startswith('refine'):
            item=item[8:]
        output_file.write(item[:7]+'\t')
    output_file.write('\n')
    dic_all.sort(key=(lambda item: item[key_sorted]), reverse=reverse_state)
    for lista in dic_all:
        for chosen_key in list_item:
            var=lista[chosen_key]
            if type(var) is float or type(var) is numpy.float64:
                var='%.3f'%(var)
                if len(str(var))>7: var='%.2f'%(float(var))
                if 'Ident%' in chosen_key or 'wMPE' in chosen_key: var='%.1f'%(float(var))
            elif type(var) is int:
                var=str(var)
            output_file.write(var+'\t')
        output_file.write('\n')
    output_file.close()




def read_setup_bor_connection (bor_file):
    Config = ConfigParser.ConfigParser()
    Config.readfp(open(bor_file))
    #extracting values from Config Parser
    #input files
    DicGridConn={}
    distribute_computing = Config.get("CONNECTION", "distribute_computing")
 # Reading the setup.bor file
    if distribute_computing == "remote_grid":
            path_bor = Config.get("CONNECTION", "setup_bor_path")
            if path_bor is None or path_bor == "" or not os.path.exists(path_bor):
                print (colored("ATTENTION: the path given for the setup.bor does not exist.\n Please contact your administrator","red"))
                sys.exit(1)
            try:
                setupbor = ConfigParser.ConfigParser()
                setupbor.readfp(open(path_bor))
                DicGridConn["username"] = setupbor.get("GRID", "remote_frontend_username")
                DicGridConn["host"] = setupbor.get("GRID", "remote_frontend_host")
                DicGridConn["port"] = setupbor.getint("GRID", "remote_frontend_port")
                DicGridConn["passkey"] = Config.get("CONNECTION", "remote_frontend_passkey")
                DicGridConn["promptA"] = (setupbor.get("GRID", "remote_frontend_prompt")).strip()+" "
                DicGridConn["isnfs"] = setupbor.getboolean("GRID", "remote_fylesystem_isnfs")
                try:
                    DicGridConn["remote_submitter_username"] = setupbor.get("GRID", "remote_submitter_username")
                    DicGridConn["remote_submitter_host"] = setupbor.get("GRID", "remote_submitter_host")
                    DicGridConn["remote_submitter_port"] = setupbor.getint("GRID", "remote_submitter_port")
                    DicGridConn["promptB"] = (setupbor.get("GRID", "remote_submitter_prompt")).strip()+" "
                except:
                    pass
                DicGridConn["home_frontend_directory"] = setupbor.get("GRID", "home_frontend_directory")
                #SELSLIB2.PATH_NEW_PHASER = setupbor.get("GRID", "path_remote_phaser")
                #SELSLIB2.PATH_NEW_SHELXE = setupbor.get("GRID", "path_remote_shelxe")
                #SELSLIB2.GRID_TYPE_R = setupbor.get("GRID","type_remote")
                #if SELSLIB2.GRID_TYPE_R == "Condor":
                    #SELSLIB2.SHELXE_REQUIREMENTS = setupbor.get("CONDOR", "requirements_shelxe")
                    #SELSLIB2.PHASER_REQUIREMENTS = setupbor.get("CONDOR", "requirements_phaser")
                    #SELSLIB2.BORGES_REQUIREMENTS = setupbor.get("CONDOR", "requirements_borges")
                #SELSLIB2.LOCAL = False
            except:
                print (colored("ATTENTION: Some keyword in your configuration files are missing. Contact your administrator","red"))
                print ("Path bor given: ",path_bor)
                print (traceback.print_exc(file=sys.stdout))
                sys.exit(1)
    cm=None

     #STARTING THE GRID MANAGER

    GRID_TYPE = ""
    if distribute_computing == "remote_grid":
        GRID_TYPE = setupbor.get("GRID","type_remote")
    elif distribute_computing == "local_grid":
        path_bor = Config.get("CONNECTION", "setup_bor_path")
        if path_bor is None or path_bor == "" or not os.path.exists(path_bor):
            print (colored("ATTENTION: the path given for the setup.bor does not exist.\n Please contact your administrator","red"))
            sys.exit(1)
        setupbor = ConfigParser.ConfigParser()
        setupbor.readfp(open(path_bor))
        GRID_TYPE = setupbor.get("GRID","type_local")

    if cm == None:
        if GRID_TYPE == "Condor":
            cm = Grid.condorManager( should_transfer_files='TRUE' )
        elif GRID_TYPE == "SGE":
            QNAME = setupbor.get("SGE","qname")
            FRACTION = setupbor.getfloat("SGE","fraction")
            cm = Grid.SGEManager(qname=QNAME,fraction=FRACTION)
        elif GRID_TYPE == "MOAB":
            PARTITION = setupbor.get("MOAB","partition")
            #FRACTION = setupbor.getfloat("MOAB","partition")
            cm = Grid.MOABManager(partition=PARTITION)
        elif GRID_TYPE == "SLURM":
            PARTITION = setupbor.get("SLURM","partition")
            if PARTITION != None and PARTITION != '':
                cm = Grid.SLURMManager(partition=PARTITION)
            else:
                cm = Grid.SLURMManager()
        elif GRID_TYPE == "TORQUE":
            QNAME = setupbor.get("TORQUE","qname")
            FRACTION = setupbor.getint("TORQUE","cores_per_node")
            PARALLEL_JOBS = setupbor.getint("TORQUE","number_of_parallel_jobs")
            MAUI = setupbor.getboolean("TORQUE","maui")
            cm = Grid.TORQUEManager(qname=QNAME,cores_per_node=FRACTION,parallel_jobs=PARALLEL_JOBS,maui=MAUI)

    if cm is not None:
        cm.setRank("kflops")
        cm.nice_user = "true"
        #TODO: Eliminate the SGE.py
        PATH_REMOTE_SGEPY = setupbor.get("GRID", "path_remote_sgepy")
        PATH_REMOTE_PYTHON_INTERPRETER = setupbor.get("GRID", "python_remote_interpreter")
        PATH_LOCAL_PYTHON_INTERPRETER = setupbor.get("LOCAL", "python_local_interpreter")

        if PATH_REMOTE_PYTHON_INTERPRETER.strip() in ["", None]:
            PATH_REMOTE_PYTHON_INTERPRETER = "/usr/bin/python"

        if PATH_LOCAL_PYTHON_INTERPRETER.strip() in ["", None]:
            PATH_LOCAL_PYTHON_INTERPRETER = "/usr/bin/python"
    
    #restrict_computers='((Machine == "nuno.ibmb.csic.es")||(Machine == "ripley.ibmb.csic.es")||(Machine == "grievous.ibmb.csic.es")||(Machine == "angora.ibmb.csic.es")||(Machine == "rumpel.ibmb.csic.es")||(Machine == "dooku.ibmb.csic.es")||(Machine == "nuno.ibmb.csic.es"))'
    #restrict_computers='(Machine == "nuno.ibmb.csic.es")||(Machine == "ripley.ibmb.csic.es")||(Machine == "grievous.ibmb.csic.es")||(Machine == "angora.ibmb.csic.es")||(Machine == "rumpel.ibmb.csic.es")||(Machine == "dooku.ibmb.csic.es")'
    #restrict_computers='(Machine == "ripley.ibmb.csic.es")||(Machine == "grievous.ibmb.csic.es")||(Machine == "rumpel.ibmb.csic.es")||(Machine == "dooku.ibmb.csic.es")'
    #cm.setRequirements(restrict_computers)
    
    DicParameters={}
    DicParameters["nameExecution"] = Config.get("SLIDER", "name_job")

    #SystemUtility.open_connection(DicGridConn,DicParameters,cm)
    return DicGridConn,DicParameters,cm


def submit_job_sh ( list_files_input , sh_file , output_directory , distribute_computing, cm ):
#def submit_job_sh ( list_files_input , list_files_output , sh_file , distribute_computing, cm ):
    if distribute_computing=='multiprocessing':
        #os.system('nohup /bin/tcsh '+sh_file+' &')
        p = subprocess.Popen(['nohup','/bin/bash',sh_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)
        out, err = p.communicate()
    else:
        job = Grid.gridJob('nameJob')
        #job.setExecutable( '/xtal/cri4-software/utilities/tcsh_condor' )
        job.setExecutable(sh_file) #suggested by max
        job.addInputFile(sh_file,False)
        #job.setStdIn([sh_file])
        for file in list_files_input:
            job.addInputFile(file,False)
        job.setInitialDir(output_directory)
        #lia = lineStdIn.split()
        sh_file=sh_file.split('/')[-1]
        job.setArguments([sh_file])
##        for file in list_files_output:
##            job.addOutputFile(file,False)
        out=sh_file[:-2]+'out'
        err=outt=sh_file[:-2]+'err'
        job.addOutputFile(out,False)
        job.addOutputFile(err,False)
        #restrict_computers written 24/Oct/2017 : || or && and
        restrict_computers = '(Machine != "cri2-01.ibmb.csic.es")'# )&&(Machine != "ray.ibmb.csic.es"))' #||(Machine == "ripley.ibmb.csic.es")||(Machine == "grievous.ibmb.csic.es")||(Machine == "angora.ibmb.csic.es")||(Machine == "rumpel.ibmb.csic.es")||(Machine == "dooku.ibmb.csic.es")||(Machine == "nuno.ibmb.csic.es"))'
##        restrict_computers='((Machine == "nuno.ibmb.csic.es")||(Machine == "ripley.ibmb.csic.es")||(Machine == "grievous.ibmb.csic.es")||(Machine == "angora.ibmb.csic.es")||(Machine == "rumpel.ibmb.csic.es")||(Machine == "dooku.ibmb.csic.es")||(Machine == "nuno.ibmb.csic.es"))'
        cm.setRequirements(restrict_computers)
##        cm.setRank("kflops")
        cm.submitJob(job)
        #(nc,nq) = cm.submitJob(job)
        #if hasattr(cm,"channel"):
        #    cm.change_remote_dir(current_dir)


def submit_job_sh_argums_execute ( list_files_input , executable , argums , output_directory , distribute_computing, cm , log_name ):
    #def submit_job_sh ( list_files_input , list_files_output , sh_file , distribute_computing, cm ):
    file_name=list_files_input[0][list_files_input[0].rindex('/')+1:]
    argums=file_name+' '+argums
    argums=argums.split()
    if distribute_computing=='multiprocessing':
        #os.system('nohup /bin/tcsh '+sh_file+' &')
        #p = subprocess.Popen([shelxe_path,dictionary_shelxe['file']]+shelxe_options.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=dictionary_shelxe['folder']) extracted from shelxe_ByDir.py
        #print executable , file_name , argums
        #p = subprocess.Popen([executable,file_name ]+argums, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_directory )
        p = subprocess.Popen([executable]+argums, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_directory  , text=True)
        out, err = p.communicate()
    else:
        job = Grid.gridJob('nameJob')
        job.setExecutable( executable )
        for file in list_files_input:
            job.addInputFile(file,False)
        job.setInitialDir(output_directory)
        job.setArguments(argums)
##        for file in list_files_output:
##            job.addOutputFile(file,False)
        out=log_name+'.out'
        err=log_name+'.err'
        job.addOutputFile(out,False)
        job.addOutputFile(err,False)
##        restrict_computers='((Machine == "nuno.ibmb.csic.es")||(Machine == "ripley.ibmb.csic.es")||(Machine == "grievous.ibmb.csic.es")||(Machine == "angora.ibmb.csic.es")||(Machine == "rumpel.ibmb.csic.es")||(Machine == "dooku.ibmb.csic.es")||(Machine == "nuno.ibmb.csic.es"))'
##        cm.setRequirements(restrict_computers)
##        cm.setRank("kflops")
        cm.submitJob(job)
        #(nc,nq) = cm.submitJob(job)
        #if hasattr(cm,"channel"):
        #    cm.change_remote_dir(current_dir)


def given_table_file_return_list_of_labels ( input_table ):
    f=open(input_table,'r')
    f2=f.readlines()
    f.close()
    labels=f2[0].split()
    return labels

def given_table_values_return_table_by_position ( input_table , sorted_by_key , output_table , reverse_state ): #give sorted_by_key=False if no key is desired to be sorted
    list_dic_all=extract_table_to_list_of_dict_with_first_line_as_key ( input_table )
    labels=given_table_file_return_list_of_labels ( input_table )
    f=open(output_table,'w')
    list_dic_all=given_list_dic_all_from_table_return_same_list_with_floats (list_dic_all)
    list_dic_new=[]
    list_dic_var=[]
    for i in list_dic_all:
        list_dic_var.append(i)
    for ind,dic in enumerate(list_dic_all):
        d={}
        d[labels[0]]=dic[labels[0]]
        for l in labels[1:]:
            if l!=sorted_by_key:
                statistic=l
                if 'CC' in statistic or statistic.startswith('ZO') or statistic.startswith('ZD-') or statistic.endswith('FOM') or statistic=='kFOM' or statistic=='align_s' or statistic=='score' or 'LLG' in statistic:
                    list_dic_var.sort(key=(lambda item: item[statistic]), reverse=True)
                    d[statistic]=list_dic_var.index(dic)
                elif statistic.startswith('R') or statistic.startswith('BA') or statistic.startswith('ZD+') or statistic.startswith('SRG') or statistic.startswith('ZDa') or statistic.startswith('ZDm') or statistic.startswith('ZDs') or statistic=='energy' or statistic.endswith('wMPE'):
                    list_dic_var.sort(key=(lambda item: item[statistic]), reverse=False)
                    d[statistic]=list_dic_var.index(dic)
                #for Table_edstats_overall_X_2
                else:
                    print (statistic,'not sorted')
                    d[statistic]='-'
            else:
                d[l]=dic[l]
                #list_dic_new[ind][statistic]=list_dic_var.index(dic)
        list_dic_new.append(d)
    
    if sorted_by_key!=False:
        list_dic_new.sort(key=(lambda item: item[sorted_by_key]), reverse=reverse_state)
    for l in labels[:-1]:
        f.write(l+'\t')
    f.write(labels[-1]+'\n')
    for line in list_dic_new:
        for l in labels[:-1]:
            f.write(str(line[l])+'\t')
        f.write(str(line[labels[-1]])+'\n')
    f.close()
   

def given_list_dic_all_from_table_return_same_list_with_floats (list_dic_all) :
    for dic in list_dic_all:
        for key in dic:
                try:
                    dic[key]=float(dic[key])
                except:
                    pass
    return list_dic_all


def read_table_file_sum_chosen_list_labels_output_newlabel_to_last_column_in_table ( input_table , list_labels , output_table , new_label_string ):
    f=open(input_table,'r')
    f2=f.readlines()
    f.close()
    labels=f2[0].split()
    f2=f2[1:]
    fw=open(output_table,'w')
    for lb in labels:
        fw.write(lb+'\t')
    fw.write( new_label_string +'\n')
    labels.append(new_label_string)
    for i_li,li in enumerate(f2):
        li=li.split()
        f2[i_li]=li
        value=0
        for lbc in list_labels:
            i_lbc = labels.index(lbc)
            if li[i_lbc]!='n/a':
                if lbc.startswith('R'):
                    value-=float(li[i_lbc])
                else:
                    value+=float(li[i_lbc])
        f2[i_li].append(value)
    f2=sorted(f2,key=(lambda item: item[-1]), reverse=True)
    for li in f2:
        for i in li[:-1]:
            fw.write(i+'\t')
        fw.write(str(li[-1])+'\n')
    fw.close()
    
    

def mkdir(folder):
    if not os.path.isdir( folder ):
        os.system('mkdir '+folder)


def shelxe_multiprocesses ( input , shelxe_path , shelxe_options ):
    if shelxe_path=='' or shelxe_path==False:
        shelxe_path=which_program ("shelxe")
    print ('Shelxe expanding',input,'with options',shelxe_options)
    folder=input[:input.rindex('/')]
    input=input[input.rindex('/')+1:]
    p = subprocess.Popen([shelxe_path,input]+shelxe_options.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=folder , text=True)
    out, err = p.communicate()
    
def return_iCC_wMPE ( lst_file ):
    f=open(lst_file)
    f2=f.read()
    f.close()
    try:
        i=f2.index('wMPE')+len('wMPE')
        wMPE=float(f2[i:i+5])
    except:
        wMPE=''
    if 'Overall CC between native Eobs and Ecalc (from fragment) =' in f2:
        i=f2.index('Overall CC between native Eobs and Ecalc (from fragment) =')+len('Overall CC between native Eobs and Ecalc (from fragment) =')
        iCC=float(f2[i:i+6])
    else:
        iCC=False
    return iCC,wMPE

def read_PDB_return_number_molec_averaged_Bfactor_occup  ( input_file ):
    filee=open(input_file)
    f=filee.readlines()
    filee.close()
    d={}
    #d=defaultdict(list)
    for l in f:
        if l.startswith('ATOM') or l.startswith('HET'):
            ch=l[21]
            resn=int(l[22:26])
            rest=l[17:20]
            occ=float(l[56:60])
            wbf=(float(l[61:66])*float(l[56:60]))
            n=1
            if rest in amino_acid_list_3L:
                rest='prot'
            if rest not in d:
                d[rest]=defaultdict(list)
            if ch not in d[rest]['ch']: d[rest]['ch'].append(ch)
            if resn not in d[rest]['resn']: d[rest]['resn'].append(resn)
            d[rest]['occ'].append(occ)
            d[rest]['wbf'].append(wbf)
    li=[]
    for key,val in d.items():
        #li= ( restype , listchain , listresnumb , avocc , wocc )
        li.append( ( key , val['ch'] , val['resn'] , numpy.mean( val['occ'] ) , (sum(val['wbf']))/(sum(val['occ'])) ) )
    return li
    
    
    
def extract_DGdiss0_pisa(log):
    f1=open(log)
    f=f1.readlines()
    f1.close()
    for i,l in enumerate(f):
        if len(l.split())>4 and l.split()[4]=='DGdiss0':
            DGdiss0=float(f[i+2].split()[4])
    return DGdiss0

def pisa_assemblies_multiproc (pdb,log):
    os.system('pisa aa -analyse '+pdb+'  --lig=fixed > delete.log')
    os.system('pisa aa -list assemblies > '+log)
    if not os.path.isfile( log ):
        print ('ERROR in pisa_assemblies_multiproc.\nFile',log,'not generated.')
        exit()



def from_table_write_coot_coordinates (list_ind,file,output_pdb): #list_ind should have three number in float/int
    f1=open(file)
    f=f1.readlines()
    f1.close()
    if not len(f[1].split())>=max(list_ind):
        print ('Table',file,'has less values than chosen indexes:',list_ind)
        exit()
    if not len(list_ind)==3:
        print ('Given numbers for list_ind different than 3')
        exit()
        
    i1=list_ind[0]
    i2=list_ind[1]
    i3=list_ind[2]

    ou=open(output_pdb,'w')

    now = datetime.datetime.now()
    date=now.isoformat()[:10]+' '+now.isoformat()[11:16]
    ou.write('REMARK   PDB written by RJB_lib.from_table_write_coot_coordinates function '+date+'\n')
    ou.write('REMARK   PDB written from table file '+file+'\n')

    at_n=0

    for l in f[1:]:
        if len(f[1].split())>=max(list_ind):
            l=l.split()
            atom_numb=int(at_n)
            #res_type='0'*(3-len(str(atom_numb)))+str(atom_numb)
            #atom_type='O'+res_type
            res_type='HOH'
            atom_type='O'
            chain='A'
            res_numb=atom_numb
            X=float(l[i1])
            Y=float(l[i2])
            Z=float(l[i3])
            occupancy=1.00
            Bfactor=20.00
            atom_numb=' '*(6-len(str(atom_numb)))+str(atom_numb)
            if len(atom_type)<4:
                atom_type='  '+atom_type
            elif len(atom_type)==4:
                atom_type=' '+atom_type
            else:
                print ('Error in function from_table_write_coot_coordinates addapted from RJB_lib.write_pdb_from_pdbRJBdic in Oct19,2016')
                print ('Wrong atom_type:',atom_type,'length='+str(len(atom_type)))
                exit()
            atom_type=atom_type+' '*(6-len(atom_type))
            res_numb=' '*(4-len(str(res_numb)))+str(res_numb)
            X='%.3f'%(X)
            X=' '*(7-len(X))+X+' '
            Y='%.3f'%(Y)
            Y=' '*(7-len(Y))+Y+' '
            Z='%.3f'%(Z)
            Z=' '*(7-len(Z))+Z+' '
            occupancy='%.2f'%(occupancy)
            Bfactor='%.2f'%(Bfactor)
            ou.write('ATOM '+atom_numb+atom_type+res_type+' '+chain+res_numb+'     '+X+Y+Z+' '+occupancy+' '+Bfactor+'           '+atom_type.replace(' ','')[0]+'\n')
            at_n+=1
        else:
            print ('Failure reading line',l,'from file:',file)
            exit()
    ou.write('END')
    ou.close()


def pymol_center_of_mass (pdb_file,selection): #both should be strings
    #insert path of pymol
    sys.path.insert(0, '/usr/lib/python2.7/dist-packages')
    import pymol
    pymol.pymol_argv = ['pymol','-qc']# + sys.argv[1:]
    pymol.finish_launching()
    cmd = pymol.cmd
    cmd.load (pdb_file)
    from pymol import center_of_mass
    com=center_of_mass.get_com(selection)
    pymol.cmd.reinitialize()
    return com

def CenterOfMassByChain (pdbinput,chain,OnlyResidues=True): #Script written 09/01/2018 to make automatic-angle-calculation_article a binary not dependent of anything else
    f2=open(pdbinput)
    f=f2.readlines()
    f2.close()
    listx=[]
    countweight=0
    listy=[]
    listz=[]
    check=True
    weight = {'H': 1, 'C': 12, 'N': 14, 'O': 16, 'S': 32}
    for l in f:
        if l.startswith('ATOM') and l[21]==chain and ( (OnlyResidues==True and l[17:20] in amino_acid_list_3L) or OnlyResidues==False):
            occup=float(l[54:60])
            x=float(l[30:38])
            y=float(l[38:46])
            z=float(l[46:54])
            atom=l[76:78].replace(' ','')
            wvar=weight[atom]
            listx.append(x*wvar*occup)
            listy.append(y*wvar*occup)
            listz.append(z*wvar*occup)
            countweight+=(wvar*occup)
    X = sum(listx) / countweight
    Y = sum(listy) / countweight
    Z = sum(listz) / countweight
    CenterMass=(X,Y,Z)
    return CenterMass


def CheckIfChainHasContinuousResidues (pdbinput,chain):
    f2=open(pdbinput)
    checkresnumb=''
    for l in f2:
        if l.startswith('ATOM') and l[17:20] in amino_acid_list_3L and l[21]==chain:
            if checkresnumb=='': checkresnumb=int(l[22:26])
            if not checkresnumb==int(l[22:26]) and not checkresnumb==int(l[22:26])-1:
                return False
    return True



def SuperimposeBioPython (pdbfix,lchainfix,lresnumbfix,pdbvar,lchainvar,lresnumbvar,atoms='CA',pdbvaroutput=False,distance=False,MatrixRotTra=False,logout=False): #atoms may be all or CA
    #l should be a list of lists related to pairs,
    # example residues 1:3 of chain A should be superimposed to residues 1:3 of chain B and residues 21:23 of chain A should be superimposed to residues 21:23 of chain B, then:
    # lchainfix=['A','A']   lresnumbfix=[[1,2,3],[21,22,23]]
    # lchainvar=['A','A']   lresnumbvar=[[1,2,3],[21,22,23]]
    # Script written 09/01/2018 to make automatic-angle-calculation_article a binary not dependent of anything else

    if isinstance(pdbfix, basestring):
        p = Bio.PDB.PDBParser(PERMISSIVE=1)
        strfix = p.get_structure(pdbfix[:-4], pdbfix)
    else:
        strfix=pdbfix

    #parser = Bio.PDB.PDBParser()
    #strfix = parser.get_structure(pdbfix[:-4], pdbfix)

    listresnfix=[]
    for i,ch in enumerate(lchainfix):
        for resn in lresnumbfix[i]:

            if atoms == 'CA': listresnfix.append( strfix[0][ch][resn]['CA'] )

            elif atoms == 'all':
                atom_list = strfix[0][ch][resn].get_unpacked_list()
                for at0 in atom_list:
                    if not at0.get_id().startswith('H'): listresnfix.append( at0 )

            else:
                print ('Wrong variable given to atoms option in function RJB_lib.SuperimposeBioPython . Exiting.')
                exit()

    if isinstance(pdbvar, basestring):
        p = Bio.PDB.PDBParser(PERMISSIVE=1)
        strvar = p.get_structure(pdbvar[:-4], pdbvar)
    else:
        strvar=pdbvar
    #strvar = parser.get_structure(pdbvar[:-4], pdbvar)

    listresnvar=[]
    for i,ch in enumerate(lchainvar):
        for resn in lresnumbvar[i]:
            if atoms == 'CA': listresnvar.append( strvar[0][ch][resn]['CA'] )

            elif atoms == 'all':
                atom_list = strvar[0][ch][resn].get_unpacked_list()
                for at0 in atom_list:
                    if not at0.get_id().startswith('H'): listresnvar.append( at0 )

            else:
                print ('Wrong variable given to atoms option in function RJB_lib.SuperimposeBioPython . Exiting.')
                exit()


    sup = Bio.PDB.Superimposer()
    sup.set_atoms(listresnfix, listresnvar)
    rmsd=sup.rms

    if pdbvaroutput != False or distance:
        sup.apply(strvar)

    if MatrixRotTra:
        MatrixRotTra=sup.rotran

    if pdbvaroutput!=False:

        io = Bio.PDB.PDBIO()
        io.set_structure(strvar)
        io.save(pdbvaroutput)

    if distance:
        ldist=[]
        if logout!=False:
            fo=open(logout,'w')
            fo.write('Chain1'+'\t'+'Res#1'+'\t'+'Atom1'+'\t'+'Chain2'+'\t'+'Res#2'+'\t'+'Atom2'+'\t'+'Dist(A)')
        for i , at1 in enumerate( listresnfix ):
            at2=listresnvar[i]
            dist=at1-at2
            ldist.append( (at1,at2,dist) )
            #at1.get_full_id()=('tCc7_tCc7', 0, 'A', (' ', 4, ' '), ('CA', ' '))
            if logout != False:
                at11     = at1.get_full_id()
                at11ch   = at11[2]
                at11resn = at11[3][1]
                at11at   = at11[4][0]
                at22     = at2.get_full_id()
                at22ch   = at22[2]
                at22resn = at22[3][1]
                at22at   = at22[4][0]
                dist='%.1f' %(dist)
                fo.write('\n'+ at11ch + '\t' + str(at11resn) + '\t' + at11at )
                fo.write('\t'+ at22ch + '\t' + str(at22resn) + '\t' + at22at )
                fo.write('\t' + dist)
        if logout != False: fo.close()
    else:
        ldist=False

    return rmsd , ldist , MatrixRotTra , strvar

def GivenListAtsDistManyModelsBioPythonReturnRMSFByRes (lista):
    ListResnRMSF=[]
    for i in range(len(lista[0])):
        lresdists = []
        resn = lista[0][i][0].get_full_id()[3][1]
        #print resn
        for ll in lista:
            lres=ll[i]
            dist=lres[-1]
            #print dist
            lresdists.append(dist)
        #print lresdists
        mean=numpy.mean(lresdists)
        stdev = numpy.std(lresdists)
        maxx=numpy.max(lresdists)
        #print stdev which is RMSF
        ListResnRMSF.append((resn, mean , stdev , maxx))

    return ListResnRMSF

def countresnumb (pdb_input):
    fi2 = open(pdb_input)
    fi = fi2.readlines()
    fi2.close()
    checkchain=''
    checkresn=''
    rescount=0
    for l in fi:
        if l.startswith('ATOM'):
            ch=l[21:22]
            resn=int(l[22:26])
            restype=l[17:20]
            if ch!=checkchain or resn!=checkresn and restype in amino_acid_list_3L:
                rescount+=1
            checkchain=ch
            checkresn=resn
    return rescount

def RemoveResidue(pdbinput, ListTupleRejected, pdboutput):
    # Done 9/1/2018
    # ListTupleRejected should be a list of tuples, whereas each tuple should contain chain, resnumb (integer).
    fi2=open(pdbinput)
    fi=fi2.readlines()
    fi2.close()
    fo=open(pdboutput,'w')
    for l in fi:
        if l.startswith('ATOM') or l.startswith('HETATM'):
            if (l[21] , int(l[22:26])) not in ListTupleRejected: fo.write(l)
        else: fo.write(l)
    fo.close()


def generate_dict_count_alignment (alig_file): #done in Dec20,2016
    ff=open(alig_file)
    f=ff.readlines()
    ff.close()
    d={}
    for i,c in enumerate(f[0][:-1]):
        d[i+1]={}
    for i,c in enumerate(f[0][:-1]):
        for ii in f:
            if ii[i]!='-' and ii[i]!='X':
                if ii[i] not in d[i+1]:
                        d[i+1][ii[i]]=1
                else:
                        d[i+1][ii[i]]+=1
    return d

def extract_RotMatrixTransl_from_log_LSQKAB (LSQKAB_log_file):
    log_file=open(LSQKAB_log_file,'r')
    log_list=log_file.readlines()
    for line_index in range (len(log_list)):
        if log_list[line_index]=='      ROTATION MATRIX:\n':
            #print log_list[line_index:]
            RotMatrix=[log_list[line_index+1].split(),log_list[line_index+2].split(),log_list[line_index+3].split()]
            for line_index2 in range (len(RotMatrix)):
                for column_index in range (len(RotMatrix[line_index2])):
                    RotMatrix[line_index2][column_index]=float(RotMatrix[line_index2][column_index])
            Transl=map(lambda n: float(n),log_list[line_index+4].split()[-3:])
            #print LSQKAB_log_file,RotMatrix,Transl
    try:
        return RotMatrix , Transl
    except:
        print ('Failure running extract_RotMatrixTransl_from_log_LSQKAB function, probably "      ROTATION MATRIX:" not found in ',LSQKAB_log_file)



def extract_proper_Tait_Bryan_angles_from_RotMatrix_James_Diebel_PDF_321 ( RotMatrix ): #RotMatrix in a tuple 3x3 Formulas extracted from: James_Diebel_PDF
    phi=degrees ( atan2( -RotMatrix[1][0] , RotMatrix[0][0] ) )
    theta=degrees ( asin( RotMatrix[2][0] ) )
    psi=degrees ( atan2( -RotMatrix[2][1] , RotMatrix[2][2] ) )
    try:
        return [psi,theta,phi]
    except:
        print ("Failure running extract_proper_Tait_Bryan_angles_from_RotMatrix_James_Diebel_PDF_321 function")
        exit()

def extract_proper_Tait_Bryan_angles_from_RotMatrix_James_Diebel_PDF_323 ( RotMatrix ): #RotMatrix in a tuple 3x3 Formulas extracted from: James_Diebel_PDF
    phi=degrees ( atan2( RotMatrix[1][2] , -RotMatrix[0][2] ) )
    theta=degrees ( acos( RotMatrix[2][2] ) )
    psi=degrees ( atan2( -RotMatrix[2][1] , RotMatrix[2][0] ) )
    try:
        return [psi,theta,phi]
    except:
        print ("Failure running extract_proper_Tait_Bryan_angles_from_RotMatrix_James_Diebel_PDF_323 function")
        exit()


def cross_product_from_3P (list_coordinates2A,list_coordinates2B,list_coordinates2C):
    v1x=-list_coordinates2A[0]+list_coordinates2B[0]
    v1y=-list_coordinates2A[1]+list_coordinates2B[1]
    v1z=-list_coordinates2A[2]+list_coordinates2B[2]
    #
    v2x=-list_coordinates2A[0]+list_coordinates2C[0]
    v2y=-list_coordinates2A[1]+list_coordinates2C[1]
    v2z=-list_coordinates2A[2]+list_coordinates2C[2]
    #
    nv=cross_product_from_2V ([v1x,v1y,v1z],[v2x,v2y,v2z])
    #
    return nv

def cross_product_from_2V (v1,v2):
    #http://tutorial.math.lamar.edu/Classes/CalcII/CrossProduct.aspx#Vectors_CrossProd_Ex2
    #https://en.wikipedia.org/wiki/Cross_product
    nx=v1[1]*v2[2]-v1[2]*v2[1]
    ny=v1[2]*v2[0]-v1[0]*v2[2]
    nz=v1[0]*v2[1]-v1[1]*v2[0]
    #
    return [nx, ny, nz]


def dot_product_from_2V (v1,v2):
    #http://tutorial.math.lamar.edu/Classes/CalcII/DotProduct.aspx#OrthogFact
    dp=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
    #print 'dp',dp
    return dp

def magnitude_vector(v):
    mv=0
    for i in v:
        mv+=i**2
    mv=sqrt(mv)
    #print 'mv',mv
    return mv

def extract_coordinates_n_atom (pdb_path, n_atom ):
    pdb_fileee=open(pdb_path,'r')
    pdb_list=pdb_fileee.readlines()
    pdb_fileee.close()
    for line in pdb_list :
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if int(line.split()[1])==n_atom:
                #print line
                x=float(line[31:38])
                y=float(line[39:45])
                z=float(line[47:54])
    try:
        return [x,y,z]
    except:
        print ("Failure looking for pdb file",pdb_path,'and atom number', n_atom)
        exit()
        
def formula_plane_from_nvector_point (vp,p):
    #http://tutorial.math.lamar.edu/Classes/CalcIII/EqnsOfPlanes.aspx
    #vp=(A,B,C)
    #p=(px,py,pz)
    #A(x-px)+B(y-py)+C(z-pz)=0
    D=(vp[0]*-p[0])+(vp[1]*-p[1])+(vp[2]*-p[2])
    return [vp[0] , vp[1] , vp[2] , D]

def distance_from_point_to_plane (point,plane_formula):
    #http://mathinsight.org/distance_point_plane
    #http://mathinsight.org/distance_point_plane_examples
    pl=plane_formula
    p=point
    d=( pl[0]*p[0] + pl[1]*p[1] + pl[2]*p[2] + pl[3] ) / sqrt(pl[0]**2+pl[1]**2+pl[2]**2)
    return d

def angle_between_two_vectors (v1,v2): #v should be v=[float_x,fl_y,fl_Z]
    angl=degrees ( acos  ( abs (  dot_product_from_2V (v1,v2) / ( magnitude_vector(v1)*magnitude_vector(v2) )   )))
    return angl

def angle_between_two_vectors_atan2 (v1,v2): #v should be v=[float_x,fl_y,fl_Z]
    angl=degrees ( atan2  (  magnitude_vector(cross_product_from_2V (v1,v2)),  dot_product_from_2V (v1,v2))   )
    return angl

def write_coord_end_pdb (input,output,ch,resn,restype,atomtype,coord):
    f1=open(input)
    ff1=f1.read()
    f1.close()
    f2=open(output,'w')
    if 'END' in ff1:f2.write(ff1[:ff1.index('END')]+'\n')
    else: f2.write(ff1+'\n' )
    natom=9000
    for i in range(len(ch)):
        chh=ch[i]
        resnn=str(resn[i])
        resnn=' '*(4-len(resnn))+resnn
        restypee=restype[i]
        restypee+=' '*(4-len(restypee))
        atomtypee=atomtype[i]
        if len(atomtypee)==4: atomtypee=' '+atomtypee+'  '
        else: atomtypee='  '+atomtypee+' '*(3-len(atomtypee))+' '
        x=coord[i][0]
        y=coord[i][1]
        z=coord[i][2]
        x='%.3f'%(x)
        x=' '*(7-len(x))+x+' '
        y='%.3f'%(y)
        y=' '*(7-len(y))+y+' '
        z='%.3f'%(z)
        z=' '*(7-len(z))+z+' '
        f2.write('ATOM   '+str(natom)+atomtypee+restypee+chh+resnn+'     '+x+y+z+' 1.00 20.00           '+atomtype[i][0]+'\n')
        natom+=1
    f2.close()
    
def correct_NM_charmm_models (file_input,file_output):
    if '/' in file_output and not os.path.isdir(file_output[:file_output.rindex('/')]): os.system('mkdir -p '+file_output[:file_output.rindex('/')])
    f2=open(file_input)
    f3=f2.readlines()
    f2.close()
    fo=open(file_output,'w')
    for l in f3:
        if l.startswith ('ATOM') and l[13]!='H' and l[12]!='H':
            l=l.replace('\n','')
            if 'HSD' in l:l=l.replace('HSD','HIS')
            if 'HSE' in l:l=l.replace('HSE','HIS')
            ln=l[:21]+l[-1]+l[22:66]+' '*11+l[13]
            fo.write(ln+'\n')
        elif l.startswith('REMARK'): fo.write(l)
    fo.close()

def gnuplot_angles (title,titleXaxis,titleYaxis,rangeXaxis,rangeYaxis,Ylabel_column_type,list_plotline_column_title_color,input_file,outputname):
    #list_columns_plot_title should be [[1,'title1'],[2,title2]]
    #Ylabel_column_type should be [ integer(index column starting in 1) , 'string'/'float' ]
    #list_plotline_column_title_color [2,'Roll','E69F00']
    #example from automatic-anglecalc: RJB_lib.gnuplot_angles (title=output_table+'_rot321_CenterMass_angles',                    titleXaxis='Structures',titleYaxis='Degrees', rangeXaxis=False,rangeYaxis=False,Ylabel_column_type=[1,'string'],list_plotline_column_title_color=[[2,'Roll' ,'CC79A7'],[3,'Twist','0090A8'],[4,'Tilt','D55E00']],input_file=output_table+'_rot321_CenterMass.log',outputname=output_table+'_rot321_CenterMass_angles.png')
    #example from auto_NM.py           RJB_lib.gnuplot_angles (title=output_table+'_angles',                    titleXaxis='Structures among '+output_table[-3:],titleYaxis='Degrees', rangeXaxis=False,rangeYaxis=False,Ylabel_column_type=[1,'float'],list_plotline_column_title_color=[[2,'Roll' ,'CC79A7'],[3,'Twist','0090A8'],[4,'Tilt','D55E00']],input_file=output_table+'_rot321_CenterMass.log',outputname=output_table+'_angles.png')
    if Ylabel_column_type[-1]=='string':
        Ylabel_column='xtic('+str(Ylabel_column_type[0])+')'
        rotateYaxis=True
    if Ylabel_column_type[-1]=='float':
        Ylabel_column=str(Ylabel_column_type[0])
        rotateYaxis=False
    gnuplot_instr=open(outputname[:-3]+"gnu","w")
    gnuplot_instr.write('set terminal png font "Verdana,24" size 1280,960\n')
    #gnuplot_instr.write("set terminal png font \"default\"\n")
    if rotateYaxis: gnuplot_instr.write("set xtics rotate\n")
    if rangeXaxis!=False: gnuplot_instr.write('set  xrange '+rangeXaxis+'\n')
    if rangeYaxis!=False: gnuplot_instr.write('set  yrange '+rangeYaxis+'\n')
    if title!=False: gnuplot_instr.write('set  title "'+title+'"\n')
    if titleYaxis!=False: gnuplot_instr.write('set  ylabel "'+titleYaxis+'"\n')
    if titleXaxis!=False: gnuplot_instr.write('set  xlabel "'+titleXaxis+'"\n')
    #gnuplot_instr.write('set xlabel font "2"\n')
    gnuplot_instr.write('set output "'+outputname+'"\n')
    gnuplot_instr.write('plot ')
    for i in list_plotline_column_title_color:
        #print i
        if i[1]!='': t=' title "'+i[1]+'" '
        else: t=' '
        #print t
        if i[-1]!='': c=' linetype rgb "#'+i[-1]+'" '
        else: c=' '
        #print c
        if rotateYaxis: colplot=str(i[0])+':'+Ylabel_column
        else: colplot=Ylabel_column+':'+str(i[0])
        #print colplot
        gnuplot_instr.write('"'+input_file+'" using '+colplot+t+c+' linewidth 5 with lines,')
    gnuplot_instr.close()
    os.system("gnuplot "+outputname[:-3]+"gnu > /dev/null")
    
    
    
    #gnuplot_instr.write('set nokey\n') #no legend
##    if input_file_list[0].startswith("RT"):
##        gnuplot_instr.write("set style line 1 lc rgb \'grey30\' ps 0 lt 1 lw 2\n")
##        gnuplot_instr.write("set style line 2 lc rgb 'grey70' lt 1 lw 2\n")
##        gnuplot_instr.write("set label '*' at 3,0.8 center\n")
##        gnuplot_instr.write("set label '*' at 4,0.8 center\n")
##        gnuplot_instr.write("set label '*' at 4,0.8 center\n")
        ##gnuplot_instr.write("set border 3\n") #excluir borda em volta do grafio
        ##gnuplot_instr.write("set xtics nomirror scale 0\n") #tracinhos para cada label em Ylabel
        






def gnuplot_matrix (input_file,outputname,title=False,titleXaxis=False,titleYaxis=False):
    gnuplot_instr=open(outputname[:-3]+"gnu","w")
    gnuplot_instr.write("set terminal png font \"default\"\n")
    if title!=False: gnuplot_instr.write('set  title "'+title+'"\n')
    if titleYaxis!=False: gnuplot_instr.write('set  ylabel "'+titleYaxis+'"\n')
    if titleXaxis!=False: gnuplot_instr.write('set  xlabel "'+titleXaxis+'"\n')
    gnuplot_instr.write('set terminal png font "Verdana,24" size 1280,960\n')
    gnuplot_instr.write('set pm3d map\n')
    gnuplot_instr.write('set xtics\n')
    gnuplot_instr.write('set ytics\n')
    gnuplot_instr.write('set cbrange [0:14]\n')
    gnuplot_instr.write('set palette rgb 34,35,36\n')
    gnuplot_instr.write('set output "'+outputname+'"\n')
    gnuplot_instr.write('splot "'+input_file+'"  matrix nonuniform')
    gnuplot_instr.close()
    os.system("gnuplot "+outputname[:-3]+"gnu > /dev/null")

##colorblindnes options
##colourset_color cb1= [0 , 0 , 0]
### black (000000)
##set_color cb2= [230 , 159 , 0]
### orange (E69F00)
##set_color cb3= [86 , 180 , 233]
### sky blue (56B4E9)
##set_color cb4= [0 , 158 , 115]
### bluewish green (009E73)
##set_color cb5= [240 , 228 , 66]
### yellow (F0E442)
##set_color cb6= [0 , 114 , 168]
### blue (0090A8)
##set_color cb7= [213 , 94 , 0]
### vermilion (D55E00)
##set_color cb8= [204 , 121 , 167]
### reddish purple (CC79A7)


def output_runline(output_file,printtt=False):
    runline = open(output_file, 'w')
    p=''
    for i in sys.argv:
        runline.write(i+' ')
        p+=i+' '
    if printtt:
        print (p)
    runline.close()

def runAREAIMOLccp4 (pdbfile,outpdb,instr_file='areaimol-instructions.ins',areaimolpath='areaimol'):
    if not os.path.isfile(instr_file):
        fo=open(instr_file,'w')
        fo.write('DIFFMODE OFF\nMODE -\n    NOHOH\nSMODE OFF\nREPORT -\n    CONTACT -\n    YES -\n    RESAREA -\n    YES\nPNTDEN 10\nPROBE 1.4\nOUTPUT\nEND')
        fo.close()
    instr=open(instr_file)
    AreaimolLog = open(outpdb[:-3]+'log', 'w')
    p = subprocess.Popen(
        [areaimolpath, 'XYZIN', pdbfile, 'XYZOUT', outpdb],stdin=instr, stdout=AreaimolLog, stderr=subprocess.PIPE , text=True)
    out, err = p.communicate()
    instr.close()
    AreaimolLog.close()
    if not os.path.isfile(outpdb):
        print ('ERROR in areaimol calculation.\nFile', outpdb, 'not generated.')
        exit()
    os.system('rm '+instr_file)

# DIFFMODE OFF
# MODE -
#     NOHOH
# SMODE OFF
# REPORT -
#     CONTACT -
#     YES -
#     RESAREA -
#     YES
# PNTDEN 10
# PROBE 1.4
# OUTPUT
# END

def ObtainAREAIMOLASA (pdbfile,DicChResN,ExcludeAtoms=[]):#,IncludeAtoms=False):
    asa=0
    with open(pdbfile) as fp:
        for l in fp.readlines():
            if l.startswith('ATOM'):
                ch = l[21]
                resn = int(l[22:26])
                atomt = l[14:16]
                if atomt not in ExcludeAtoms and (DicChResN==False or (ch in DicChResN and resn in DicChResN[ch])):
                    asa+=float(l[60:66])
    return asa

def remove_pdb_double_occupancy_atoms_pdb (input_pdb , output_pdb): #adapted from http://biopython.org/wiki/Remove_PDB_disordered_atoms
    parser = Bio.PDB.PDBParser()
    s = parser.get_structure(input_pdb[:-4], input_pdb)
    io = Bio.PDB.PDBIO()

    keepAltID = 'A'

    class NMROutputSelector2(Bio.PDB.Select):  # Inherit methods from Select class #taken from http://biopython.org/wiki/Remove_PDB_disordered_atoms
        def accept_atom(self, atom):
            if atom.get_parent().resname not in amino_acid_list_3L: #exclude non protein atoms
                return True
            else:
                if ( (not atom.is_disordered()) or atom.get_altloc() == keepAltID) :
                    atom.set_altloc(' ')  # Eliminate alt location ID before output.
                    return True
                else:  # Alt location was not one to be output.
                    return False
                    # end of accept_atom()

    # end of NMROutputSelector2()
    check=False
    dic_ch_list_resnumb_disordered=defaultdict(list)
    for residue in s.get_residues():
        if residue.is_disordered():
            check=True
            #dic={}
            chain=residue.get_parent().id
            resname=residue.resname
            resnumb=residue.id[1]
            if resname in amino_acid_list_3L: dic_ch_list_resnumb_disordered[chain].append(resnumb)
            #print residue
    if check:
        io = Bio.PDB.PDBIO()
        io.set_structure(s)
        io.save(output_pdb, select=NMROutputSelector2())
        # Note that the code above does not eliminate the alternate location identifier (‘A’ in the example above). It is the programmer’s responsibility to eliminate the identifier when necessary.
    #print check
    if check: return dic_ch_list_resnumb_disordered
    else: return False


def phaser_Ani_tNCS_correction (input_mtz,F,SigF,Amplitudes,sh_file,log_file): #Amplitudes should be True or False
    print ('Running phaser to perform ANISOTROPY and tNCS CORRECTION')
    f=open(sh_file,'w')
    f.write('#!/bin/tcsh\n')
    f.write('phaser << EOF - phaser\n')
    f.write('MODE NCS\n')
    f.write('MR_NCS\n')
    f.write('HKLIN "'+input_mtz+'"\n')
    if Amplitudes: f.write('LABIN F='+F+' SIGF='+SigF+'\n')
    else:          f.write('LABIN I='+F+' SIGI='+SigF+'\n')
    f.write('TITLE Anisotropy and TNCS Correction\n')
    f.write('TNCS EPSFAC WRITE anis.tncs\n')
    f.write('NORM EPSFAC WRITE anis.norm\n')
    f.write('ROOT anis\n')
    f.write('END\n')
    f.write('EOF-phaser\n')
    f.close()
    os.system('chmod 755 ./'+sh_file)
    os.system('./'+sh_file+' > ' + log_file)

def PhaserMR_AUTO_Map (InputDataMTZ,FData,SigFData,InputMapMTZ,extent,center,rms,Mw,outputfile,PHASER_path='phaser'):#,FMap='FWT',SigFMap='PHWT'): extend/center=[x,y,z]
    print ('Running PHASER to perform location of fragment of map',InputMapMTZ,'and saving to',outputfile)
    f=open(outputfile+'.sh','w')
    f.write('#!/bin/tcsh\n')
    f.write(PHASER_path+' << EOF - phaser\n')
    f.write('MODE MR_AUTO\n')
    f.write('HKLIN "'+InputDataMTZ+'"\n')
    f.write('LABIN F='+FData+' SIGF='+SigFData+'\n')
    f.write('ENSEMBLE ense_1 HKLIN "' + InputMapMTZ + '" F=FWT PHI=PHWT EXTENT '+extent[0]+' '+extent[1]+' '+extent[2]+' RMS '+rms+' CENTRE '+center[0]+' '+center[1]+' '+center[2]+' PROTEIN MW '+Mw+' NUCLEIC MW 0 CELL SCALE 1\n')
    f.write('SEARCH ENSEMBLE ense_1 NUM 1\n')
    f.write('ROOT '+outputfile+'\n')
    f.write('END\n')
    f.write('EOF-phaser\n')
    f.close()
    # f.write('TITLE Anisotropy and TNCS Correction\n')
    # f.write('TNCS EPSFAC WRITE anis.tncs\n')
    # f.write('NORM EPSFAC WRITE anis.norm\n')
    # f.write('ROOT anis\n')
    os.system('/bin/tcsh '+outputfile+'.sh > ' + outputfile+'.log')
    #os.system('./'+sh_file+' > ' + log_file)


# /xtal/phaser_nightly/dev-2733/phaser_nightly << eof
# MODE MR_AUTO
# HKLIN "/localdata/rafael/SLIDER_phasing/2018/Map/1YZF/1YZF.mtz"
# LABIN F=FOBS SIGF=SIGFOBS
# ENSEMBLE ense_1 HKLIN "/localdata/rafael/SLIDER_phasing/2018/Map/1YZF/maps/map0out.mtz" F=FWT PHI=PHWT EXTENT 12 12 12 RMS 0.5 CENTRE  PROTEIN MW 10000 NUCLEIC MW 0 CELL SCALE 1
# SEARCH ENSEMBLE ense_1 NUM 1
# ROOT AUTO_ense_1
# eof



def calculate_LLG (input_pdb,input_mtz,SG,F,SigF,Amplitudes,HighRes,Mw,nASU,RMS,sh_file,log_file,VRMS=False,RigidBodyRef=False,Gimble=False):
    if VRMS: VRMS='ON'
    else:    VRMS='OFF'
    if Gimble: Gimble='ON'
    else:      Gimble='OFF'
    if not RigidBodyRef: print ('Running phaser to calculate LLG of model:',input_pdb)
    else:                print ('Running phaser to perform Rigid Body Refinement and LLG calculation of model:',input_pdb)
    f=open(sh_file,'w')
    f.write('#!/bin/tcsh\n')
    f.write('phaser << EOF - phaser\n')
    f.write('MODE MR_RNP\n')
    f.write('HKLIN "'+input_mtz+'"\n')
    f.write('HKLOUT OFF\n')
    if Amplitudes: f.write('LABIN F='+F+' SIGF='+SigF+'\n')
    else:          f.write('LABIN I='+F+' SIGI='+SigF+'\n')
    if not RigidBodyRef: f.write('TITLE Calculate LLG\n')
    else:                f.write('TITLE RBR and LLG Calculation\n')
    f.write('JOBS 1\n')
    f.write('SGALTERNATIVE SELECT NONE\n')
    f.write('COMPOSITION PROTEIN MW '+Mw+' NUMBER '+nASU+'\n')
    #f.write('MACMR PROTOCOL OFF\n')
    if not RigidBodyRef: f.write('MACMR ROT OFF TRA OFF BFAC OFF VRMS '+VRMS+' CELL OFF LAST OFF NCYCLE 1000\n')#NCYCLE <NCYC> MINIMIZER [BFGS|NEWTON|DESCENT]\n')
    elif RigidBodyRef:   f.write('MACMR ROT ON  TRA ON  BFAC OFF VRMS '+VRMS+' CELL OFF LAST ON \nMACMR  CHAINS '+Gimble+'\n') #Changed to LOW Mahan Case
    #elif RigidBodyRef:   f.write('MACMR ROT ON  TRA ON  BFAC OFF VRMS '+VRMS+'\nMACMR  CHAINS '+Gimble+'\n')
    f.write('MACANO PROTOCOL OFF\n')
    f.write('MACTNCS PROTOCOL OFF\n')
    f.write('TNCS EPSFAC READ anis.tncs\n')
    f.write('NORM EPSFAC READ anis.norm\n')
    f.write('RESOLUTION 99 '+HighRes+'\n')
    f.write('SOLPARAMETERS BULK USE OFF\n')
    if not RigidBodyRef: f.write('XYZOUT OFF\n')
    else:                f.write('XYZOUT ON\nTOPFILES 34\n')
    #f.write('TOPFILES 34\n')
    #f.write('ENSEMBLE ensarci0 PDBFILE '+input_pdb+' RMS 0.1\n')
    #f.write('ENSEMBLE ensarci1 PDBFILE marenhigh.pdb RMS '+RMS+'\n')#'#0.1\n')
    #f.write('ENSEMBLE ensarci1 DISABLE CHECK ON\n')
    f.write('ENSEMBLE ensarci0 PDBFILE ' + input_pdb + ' RMS '+RMS+'\n')#'#0.1\n')
    f.write('ENSEMBLE ensarci0 DISABLE CHECK ON\n')
    f.write('SOLU SET\n')
    f.write('SOLU SPAC '+SG+'\n')
    #f.write('SOLU 6DIM ENSE ensarci1 EULER 	0.0 0.0 0.0	FRAC 0.0 0.0 0.0	BFAC 0.0\n')
    f.write('SOLU 6DIM ENSE ensarci0 EULER 	0.0 0.0 0.0	FRAC 0.0 0.0 0.0	BFAC 0.0\n')
    f.write('ROOT "'+input_pdb[:-4]+'_phaser_out"\n')
    f.write('END\n')
    f.write('EOF-phaser\n')
    f.close()
    #os.system('chmod 755 '+sh_file)
    os.system('tcsh '+sh_file+' > ' + log_file)


#def calculate_LLG (input_pdb,input_mtz,SG,F,SigF,Amplitudes,HighRes,Mw,nASU,eRMSD,sh_file,log_file,VRMS=True): #Before Cambridge 2019
    # if VRMS: VRMS='ON'
    # else:    VRMS='OFF'
    # print 'Running phaser to calculate LLG of model:',input_pdb
    # f=open(sh_file,'w')
    # f.write('#!/bin/tcsh\n')
    # f.write('phaser << EOF - phaser\n')
    # f.write('MODE MR_RNP\n')
    # f.write('HKLIN "'+input_mtz+'"\n')
    # f.write('HKLOUT OFF\n')
    # if Amplitudes: f.write('LABIN F='+F+' SIGF='+SigF+'\n')
    # else:          f.write('LABIN I='+F+' SIGI='+SigF+'\n')
    # f.write('TITLE Calculate LLG\n')
    # f.write('JOBS 1\n')
    # f.write('SGALTERNATIVE SELECT NONE\n')
    # f.write('COMPOSITION PROTEIN MW '+Mw+' NUMBER '+nASU+'\n')
    # #f.write('MACMR PROTOCOL OFF\n')
    # f.write('MACMR ROT OFF TRA OFF BFAC OFF VRMS '+VRMS+' CELL OFF LAST OFF NCYCLE 1000\n')#NCYCLE <NCYC> MINIMIZER [BFGS|NEWTON|DESCENT]\n')
    # f.write('MACANO PROTOCOL OFF\n')
    # f.write('MACTNCS PROTOCOL OFF\n')
    # f.write('TNCS EPSFAC READ anis.tncs\n')
    # f.write('NORM EPSFAC READ anis.norm\n')
    # f.write('RESOLUTION 99 '+HighRes+'\n')
    # f.write('SOLPARAMETERS BULK USE OFF\n')
    # f.write('XYZOUT OFF\n')
    # f.write('TOPFILES 34\n')
    # #f.write('ENSEMBLE ensarci0 PDBFILE '+input_pdb+' RMS 0.1\n')
    # f.write('ENSEMBLE ensarci0 PDBFILE ' + input_pdb + ' RMS '+eRMSD+'\n')#'#0.1\n')
    # f.write('ENSEMBLE ensarci0 DISABLE CHECK ON\n')
    # f.write('SOLU SET\n')
    # f.write('SOLU SPAC '+SG+'\n')
    # f.write('SOLU 6DIM ENSE ensarci0 EULER 	0.0 0.0 0.0	FRAC 0.0 0.0 0.0	BFAC 0.0\n')
    # f.write('ROOT "0"\n')
    # f.write('END\n')
    # f.write('EOF-phaser\n')
    # f.close()
    # #os.system('chmod 755 '+sh_file)
    # os.system('tcsh '+sh_file+' > ' + log_file)

def return_LLG(phaser_log_file):
    f2=open(phaser_log_file)
    f=f2.readlines()
    f2.close()
    for i,l in enumerate(f):
        #if l.startswith('   Refined TF/TFZ equivalent ='):
        if l.startswith('$$ loggraph $$'):
            LLG=float( f[i+1].split()[1] )
            #TFZ= float( l[l.rindex('/')+1:l.rindex('(')] )
            #return eLLG,TFZ
            return LLG
    #print 'Failure in function RJB_lib.return_LLG , it did not found eLLG and TFZ, exiting.'
    print ('Function RJB_lib.return_LLG failed to find LLG in file: '+phaser_log_file+' . Exiting.')
    exit()

def return_VRMS(phaser_log_file):
    f2=open(phaser_log_file)
    f=f2.read()
    f2.close()
    i=len('   SOLU ENSEMBLE ensarci0 VRMS DELTA -0.9944 #RMSD  1.00 #VRMS  ')
    ivrms=f.index('   SOLU ENSEMBLE ensarci0 VRMS DELTA')
    vrms=f[ ivrms+i:ivrms+i+4]
    return f[ ivrms+i:ivrms+i+4]
    # print 'Function RJB_lib.return_LLG failed to find LLG in file: '+phaser_log_file+' . Exiting.'
    # exit()


def runPhenixMapMtzMtzCC (mtz1,mtz2,log,MtzMtzCC_path='phenix.get_cc_mtz_mtz',keep=False):
    print ('Running phenix.get_cc_mtz_mtz to calculate MapCC between',mtz1,'and',mtz2)
    if '/' in log:
        folder=log[:log.rindex('/')+1]
        log=log[log.rindex('/')+1:]
    else:  folder=''
    #os.system(MtzMtzCC_path+' '+mtz1+' '+mtz2+' scale=True keep_f_mag=False use_only_refl_present_in_mtz_1=True  offset_mtz=offset_'+log[:-3]+'mtz'+' temp_dir="'+folder+'z_temppp_'+log[:-4]+'" > '+folder+log)
    os.system(
        MtzMtzCC_path + ' ' + mtz1 + ' ' + mtz2 + ' offset_mtz=offset_' + log[:-3] + 'mtz' + ' temp_dir="' + folder + 'z_temppp_' + log[:-4] + '" use_only_refl_present_in_mtz_1=True scale=True keep_f_mag=False  > ' + folder + log)
    #print MtzMtzCC_path+' '+mtz1+' '+mtz2+' scale=True keep_f_mag=False use_only_refl_present_in_mtz_1=True  offset_mtz=offset_'+log[:-3]+'mtz'+' temp_dir="z_temppp_'+log[:-4]+'" > '+log
    if not keep:
        os.system('rm -r offset.log offset_'+log[:-3]+'mtz '+folder+'z_temppp_'+log[:-4])
    else: os.system('mv  offset.log offset_'+log)

def get_PhenixMapMtzMtzCC (log):
    f=open(log)
    f2=f.readlines()
    #f2.reverse()
    f.close()
    for i in range(len(f2)-1,-1,-1):
        l=f2[i]
        if l.startswith('Starting correlation:'):
            cci=float(l.split()[2])
            ccf=float(l.split()[-1])
            return cci,ccf
    print ('Failure obtaining mapCC from file:',log)

def runPhenixCutMapInBoxGivenCoord (pdb,mtz,outmtz,PhenixCutMap_path='phenix.cut_out_density',keep=False):
    print ('Running phenix.cut_out_density on box '+pdb+' and mtzmap '+mtz)
    os.system(PhenixCutMap_path+' pdb_in='+pdb+' mtz_in='+mtz+' mtz_out='+outmtz+' temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out='+outmtz[:-3]+'eff cutout_fix_position=True cutout_type=box cutout_model_radius=0 cutout_sphere_radius=0 padding=0 > '+outmtz[:-3]+'log')
    if not keep: os.system('rm -r '+outmtz[:-4]+'_temp1 '+outmtz[:-4]+'_output '+outmtz[:-3]+'eff')

def runPhenixCutMapSpherePhenix (pdb,mtz,outmtz,SphereRadius,PhenixCutMap_path='phenix.cut_out_density',keep=False):
    print ('Running phenix.cut_out_density using phaser defaults settings around a sphere of '+str(SphereRadius)+' around coordinates on '+pdb+' using mtzmap '+mtz+' and saving map in '+outmtz)
    os.system(PhenixCutMap_path+' pdb_in='+pdb+' mtz_in='+mtz+' mtz_out='+outmtz+' temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out='+outmtz[:-3]+'eff for_phaser=True cutout_model_radius='+str(SphereRadius)+'  > '+outmtz[:-3]+'log')
    #print     PhenixCutMap_path+' pdb_in='+pdb+' mtz_in='+mtz+' mtz_out='+outmtz+' temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out='+outmtz[:-3]+'eff for_phaser=True cutout_model_radius='+str(SphereRadius)+'  > '+outmtz[:-3]+'log'
    #exit()
    #phenix.cut_out_density pdb_in=input.pdb mtz_in=input.mtz mtz_out=output.mtz temp_dir=remove output_dir=remove2 params_out=remove3.eff for_phaser=True cutout_model_radius=10  > output.log
    if not keep: os.system('rm -r '+outmtz[:-4]+'_temp1 '+outmtz[:-4]+'_output '+outmtz[:-3]+'eff')

def runPhenixCutMapGivenAroundCoord (pdb,mtz,outpdb,outmtz,selection=False,PhenixCutMap_path='phenix.cut_out_density',keep=False):
    print ('Running phenix.cut_out_density on atoms in '+pdb+' and mtzmap '+mtz)
    if selection!=False:
        if 'atom_selection' in selection: sel=selection
        else: sel='atom_selection'+'="'+selection+'"'
    else: sel=''
    #print '\n\n\n'+PhenixCutMap_path+' pdb_in="'+pdb+'" mtz_in="'+mtz+'" mtz_out="'+outmtz+'" pdb_out="'+outpdb+'" temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out="'+outmtz[:-3]+'eff" '+sel+' padding=0.5 cutout_type=model cutout_sphere_radius=10 cutout_model_radius=0.5 > "'+outmtz[:-3]+'log"'
    output_cc=open(outmtz[:-3] + 'log','w')
    p = subprocess.Popen (PhenixCutMap_path+' pdb_in="'+pdb+'" mtz_in="'+mtz+'" mtz_out="'+outmtz+'" pdb_out="'+outpdb+'" temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out="'+outmtz[:-3]+'eff" '+sel+' padding=0.5 cutout_type=model cutout_sphere_radius=10 cutout_model_radius=0.5 ' , shell=True , stdin=subprocess.PIPE,stdout=output_cc,stderr=output_cc , text=True)
    outp,errorp=p.communicate()
    output_cc.close()
    #os.system(PhenixCutMap_path+' pdb_in="'+pdb+'" mtz_in="'+mtz+'" mtz_out="'+outmtz+'" pdb_out="'+outpdb+'" temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out="'+outmtz[:-3]+'eff" '+sel+' padding=0.5 cutout_type=model cutout_sphere_radius=10 cutout_model_radius=0.5 > "'+outmtz[:-3]+'log"')
    if not keep: os.system('rm -r '+outmtz[:-4]+'_temp1 '+outmtz[:-4]+'_output '+outmtz[:-3]+'eff')



def runPhenixSuperposeMapsFromModels (pdb1_var,pdb2_fixed,mtz1_var,mtz2_fixed,outpdb,outmtz,selection=False,PhenixCutMap_path='phenix.cut_out_density',keep=False):
    print ('Running phenix.cut_out_density on atoms in '+pdb+' and mtzmap '+mtz)
    if selection!=False:
        if 'atom_selection' in selection: sel=selection
        else: sel='atom_selection'+'="'+selection+'"'
    else: sel=''
    #print '\n\n\n'+PhenixCutMap_path+' pdb_in="'+pdb+'" mtz_in="'+mtz+'" mtz_out="'+outmtz+'" pdb_out="'+outpdb+'" temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out="'+outmtz[:-3]+'eff" '+sel+' padding=0.5 cutout_type=model cutout_sphere_radius=10 cutout_model_radius=0.5 > "'+outmtz[:-3]+'log"'
    output_cc=open(outmtz[:-3] + 'log','w')
    p = subprocess.Popen (PhenixCutMap_path+' pdb_in="'+pdb+'" mtz_in="'+mtz+'" mtz_out="'+outmtz+'" pdb_out="'+outpdb+'" temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out="'+outmtz[:-3]+'eff" '+sel+' padding=0.5 cutout_type=model cutout_sphere_radius=10 cutout_model_radius=0.5 ' , shell=True , stdin=subprocess.PIPE,stdout=output_cc,stderr=output_cc , text=True)
    outp,errorp=p.communicate()
    output_cc.close()
    #os.system(PhenixCutMap_path+' pdb_in="'+pdb+'" mtz_in="'+mtz+'" mtz_out="'+outmtz+'" pdb_out="'+outpdb+'" temp_dir="'+outmtz[:-4]+'_temp1" output_dir="'+outmtz[:-4]+'_output" params_out="'+outmtz[:-3]+'eff" '+sel+' padding=0.5 cutout_type=model cutout_sphere_radius=10 cutout_model_radius=0.5 > "'+outmtz[:-3]+'log"')
    if not keep: os.system('rm -r '+outmtz[:-4]+'_temp1 '+outmtz[:-4]+'_output '+outmtz[:-3]+'eff')


def extract_edges_coordinates_pdb_write_pdb (pdb,atoms_accept,dist_cut,pdbout):
    #dicall=defaultdict(list)
    dist_cut=float(dist_cut)
    dicall={}
    f=open(pdb)
    f2=f.readlines()
    f.close()
    for l in f2:
        if l.startswith('ATOM') or l.startswith('HETATM'):
            chain_ID = l[21]
            atom_type = l[13:15].replace(' ','')
            if atom_type in atoms_accept:
                x = float(l[31:38])
                y = float(l[39:46])
                z = float(l[47:54])
                if chain_ID not in dicall:
                    dicall[chain_ID]=[[x],[y],[z]]
                else:
                    dicall[chain_ID][0].append(x)
                    dicall[chain_ID][1].append(y)
                    dicall[chain_ID][2].append(z)

    dic_lines_bychain_extremities = {}
    for ch,listall in dicall.items():
        lx=[min(listall[0])-dist_cut,max(listall[0])+dist_cut]
        ly=[min(listall[1])-dist_cut,max(listall[1])+dist_cut]
        lz=[min(listall[2])-dist_cut,max(listall[2])+dist_cut]
        select=[lx,ly,lz]
        #print select

        dic_lines_bychain_extremities[ch]=''

        if pdbout.endswith('.pdb'): pdbout=pdbout[:-4]
        fo=open(pdbout+'_extremities_chain'+ch+'.pdb','w')
        i=9001
        for x in lx:
            for y in ly:
                for z in lz:
                    #print type(x)
                    xx='%.3f'%(x)
                    xx=' '*(7-len(xx))+xx+' '
                    yy='%.3f'%(y)
                    yy=' '*(7-len(yy))+yy+' '
                    zz='%.3f'%(z)
                    zz=' '*(7-len(zz))+zz+' '
                    atomnumb=str(i)
                    atomnumb=' '*(5-len(str(atomnumb)))+str(atomnumb)
                    i+=1
                    ll='ATOM  '+atomnumb+'  O   '+'HOH'+' '+ch+atomnumb[1:]+'     '+xx+yy+zz+' 1.00 20.00           O  \n'
                    fo.write(ll)
                    dic_lines_bychain_extremities[ch]+=ll
        fo.close()
    return dic_lines_bychain_extremities

def cutmap_calculate_CC (boxpdb,mtz1_input,mtz1_cut,mtz2_ent,PhenixCutMap_path,MtzMtzCC_path,keep=False):
    runPhenixCutMapGivenCoord (pdb,mtz,outmtz,PhenixCutMap_path='phenix.cut_out_density',keep=False)
    runPhenixMapMtzMtzCC (mtz1,mtz2,log,MtzMtzCC_path='phenix.get_cc_mtz_mtz',keep=False)


def cutmap_coord_calculate_CC (pdb_input,mtz1_input,pdb1_cut,mtz1_cut,mtz2_ent,log,selection=False,PhenixCutMap_path='phenix.cut_out_density',MtzMtzCC_path='phenix.get_cc_mtz_mtz',keep=False):
    runPhenixCutMapGivenAroundCoord(pdb=pdb_input, mtz=mtz1_input, outpdb=pdb1_cut, outmtz=mtz1_cut, selection=selection,PhenixCutMap_path=PhenixCutMap_path, keep=keep)
    runPhenixMapMtzMtzCC (mtz1=mtz2_ent,mtz2=mtz1_cut,log=log,MtzMtzCC_path=MtzMtzCC_path,keep=keep)


def listres0occup(pdb):
    f2=open(pdb)
    f=f2.readlines()
    lres0occ=[]
    for l in f:
        if l.startswith ('ATOM') and l[56:60]=='0.00' and l[13:15]=='CA':
            ch=l[21]
            restype=amino_acid_list [ amino_acid_list_3L.index(l[17:20])]
            resnumb=int(l[22:26])
            tup=(ch, resnumb)
            lres0occ.append(tup)
            #lres0occ.append ( (ch,resnumb,restype) )

    return lres0occ

def GivenListRes0occScwrl4PDBRecover0occ (lres0occ,pdb_input,pdb_output):
    with open(pdb_input) as f: f2=f.readlines()
    with open(pdb_output,'w') as f:
        for l in f2:
            if not l.startswith('ATOM'): f.write(l)
            else:
                ch=l[21]
                resn=int(l[22:26])
                if not (ch,resn) in lres0occ: f.write(l)
                else:                         f.write(l[:56]+'0.00'+l[60:])

def convertBORGES_MATRIXfrags2listSSpdb (pdbinicial,pdbborgesmatrix):
    residuespdbinicial=return_tuple_chain_resnumb_restype(pdbinicial)
    listSSpdb=[]
    for i1 in residuespdbinicial:
        for i2 in i1:
            tupleres=i2
            sstype='coil'
            ch=i2[0]
            resn=i2[1]
            for i3 in pdbborgesmatrix:
                reslist = i3['reslist']
                if i3['sstype']=='coil': continue
                sstypev=i3['sstype']
                for res in reslist:
                    chv=res[2]
                    resnv=res[3][1]
                    if chv==ch and resn==resnv: sstype=sstypev
            listSSpdb.append(tupleres+tuple([sstype]))
    return listSSpdb

#example listSSpdb: [('A', 7, 'A', 'coil'), ('A', 8, 'A', 'ah'), ('A', 9, 'A', 'ah'), ('A', 10, 'A', 'ah'),('F', 32, 'A', 'ah'), ('F', 33, 'A', 'ah'), ('F', 34, 'A', 'ah')]

def countnumberexpectedATOMS (sequence,nASU):#,resolution):
    atomspermonomer=0
    for i,a in enumerate (amino_acid_list[:-1]):
        print (i,a,sequence.count(a),amino_acid_list_numb_atoms[i])
        atomspermonomer+= (sequence.count(a) * amino_acid_list_numb_atoms[i])
    atomspermonomer=atomspermonomer*nASU
    #expectedwaters=nASU
    return atomspermonomer

def calculateRMSDfromdistances (listdist):
    listdistSQ=[x**2 for x in listdist]
    rmsd=(sum(listdistSQ)/len(listdistSQ))**0.5
    return rmsd

def readCLUSTALOoutput (alignfile):
    dic={}
    from Bio import AlignIO
    align = AlignIO.read(alignfile, "clustal")
    # print type(align)
    for a in align:
    #     print a
    # print(align)
        dic[a.id]=a.seq
        # print a.id,a.seq
    #record.seq + " " + record.id
    #exit()
    return dic

def readPIRalignment (alignfile):
    dic={}
    from Bio.SeqIO import PirIO
    with open (alignfile) as f:
        for record in PirIO.PirIterator(f):
            #print record.seq
            dic[record.id]=record.seq
    return dic


#written in 11May2018 to convert dic['A']=[1,2,3,4] to folder string seq_A1-4 for SLIDER1.9.py
def convertDicChResNStr ( DicChResNStr ):
    strr = 'seq'
    for k, v in DicChResNStr.items():
        strr += '_' + k
        for n in v:
            if not n - 1 in v:
                strr += str(n)
            if not n + 1 in v:
                strr += '-' + str(n) + '_'
        strr = strr[:-1]
    return strr

def GivenEdstatsOutDicChResnStats12ReturnList (RSS_output_file,DicChResn,Stat1,Stat2): #Stat usually CCSm or CCSs or CCSa m,s,a=main,side,all atoms
    f2=open(RSS_output_file)
    f=f2.readlines()
    llabels=f[0].split()
    iCh=llabels.index('CI')
    iResn=llabels.index('RN')
    if Stat1 in llabels: iStat1 = llabels.index(Stat1)
    else: print (Stat1,'not found in',RSS_output_file,'"n/a" value being atributed.')
    if Stat2 in llabels: iStat2 = llabels.index(Stat2)
    else: print (Stat2,'not found in',RSS_output_file,'"n/a" value being atributed.')
    lvar1 = []
    lvar2 = []
    for l in f[1:]:
        l=l.split()
        if l[iCh] in DicChResn and int(l[iResn]) in DicChResn[l[iCh]]:
            if Stat1 in llabels and l[iStat1]!='n/a':  lvar1.append(float( l[iStat1] ))
            if Stat2 in llabels and l[iStat2] != 'n/a':lvar2.append(float( l[iStat2] ))
    if lvar1 == []: lvar1 = 0
    if lvar2 == []: lvar2 = 0
    return lvar1,lvar2

def BestAlignment2StringsReturnIndex (st1,st2):
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
    return bestindex

def GivenTwoSequencesIdenticalRes(st1,st2):
    if len (st1)!=len(st2):
        print (st1,st2,"with different lengths (",len (st1),len(st2),")")
        exit()
    id=0
    for i in range(len(st1)):
        if st1[i]==st2[i]: id+=1
    return id

def CheckHeaderPDB (inputpdb):
    checknumb=0
    with open(inputpdb) as f:
        f2 = f.readlines()
        for l in f2:
            if l.startswith('SCALE') or l.startswith('CRYST1'): checknumb+=1
    if checknumb!=4:
        print ('REQUIRED INFORMATION IN PDB HEADER NOT FOUND (CRYST1/SCALE CARDS). CORRECT THE FILE:',inputpdb)
        exit()

def ObtainExtentCentreLogPHENIXCutOutDensity (logfile):
    extent=''
    centre=''
    with open(logfile) as f:
        for l in f.readlines():
            if l.startswith('Center of cutout region will be from model center at:'):
                centre=l[ l.index('[')+1:l.index(']')].replace(',','')
            elif l.startswith('Dimensions of cutout region will be:'):
                extent=l[ l.index('[')+1:l.index(']')].replace(',','')
    if extent!='' and centre!='': return extent , centre
    else:
        print ('Failure finding string:')
        print ('Center of cutout region will be from model center at:')
        print ('and/or')
        print ('Dimensions of cutout region will be:')
        print ('in file:',logfile)
        exit()

def ReturnShelxeCCwMPEPdbLlst ( pdb_file,lst_file ):
    with open(lst_file) as ff:
        f2=ff.readlines()
        f3=ff.read()
        if '** Unable to trace map - giving up **' in f3:
            return 'n/a','n/a','n/a','n/a','n/a' #CC,wMPE,cycle,nres,nchains
        else:
            with open(pdb_file) as f:
                l1stline = f.readlines()[0].split()
                CC = float(l1stline[6][:-1])
                cycle = int(l1stline[3])
                nres = int(l1stline[7])
                nchains = int(l1stline[10])
                # print CC,cycle,nres,nchains
            count = 0
            for l in f2:
                if l.startswith(' CC for partial structure against native data ='):
                    count+=1
                    if count==cycle:
                        wMPE = float ( f2 [f2.index(l) + 3 ] .split()[-3] )
                        break
    return CC,wMPE,cycle,nres,nchains

# def ReturnShelxeCCwMPEPdbLlst ( pdb_file,lst_file ):
#     with open(lst_file) as ff:
#         f2=ff.readlines()
#         f3=ff.read()
#         if '** Unable to trace map - giving up **' in f3:
#             return 'n/a','n/a','n/a','n/a','n/a' #CC,wMPE,cycle,nres,nchains
#         else:
#             with open(pdb_file) as f:
#                 l1stline = f.readlines()[0].split()
#                 CC = float(l1stline[6][:-1])
#                 cycle = int(l1stline[3])
#                 nres = int(l1stline[7])
#                 nchains = int(l1stline[10])
#                 # print CC,cycle,nres,nchains
#             count = 0
#             for l in f2:
#                 if l.startswith(' CC for partial structure against native data ='):
#                     count+=1
#                     if count==cycle:
#                         wMPE = float ( f2 [f2.index(l) + 3 ] .split()[-3] )
#                         break
#     return CC,wMPE,cycle,nres,nchains


def PhenixModelVsData (pdb_file,mtz_data,F,output_file,PhenixModelVsData_path=False):
    if not PhenixModelVsData_path: PhenixModelVsData_path='phenix.model_vs_data'
    print (PhenixModelVsData_path+' '+pdb_file+' '+mtz_data+' f_obs_label="'+F+'" comprehensive=true > '+output_file)
    os.system(PhenixModelVsData_path+' '+pdb_file+' '+mtz_data+' f_obs_label="'+F+'" comprehensive=true > '+output_file)

def ExtractFromPhenixModelVsDataRSCCallList (log):
    with open(log) as f:
        contentfile=f.read()
        ind1=contentfile.index('  Overall map cc(Fc,2mFo-DFc): ')+len('  Overall map cc(Fc,2mFo-DFc): ')
        RSCCall=float(contentfile[ind1+1:ind1+7])
        contentfilelist=contentfile[contentfile.index('id string')+4:].split('\n')
        ccindex=contentfilelist[0].split()[::-1].index('CC')*-1-1
        RSCClist=[]
        for l in contentfilelist[1:]:
            l=l.split()
            #print 'l',l,'RSCC',l[ccindex]
            if len(l)>1: RSCClist.append(float(l[ccindex]))
    return RSCCall,RSCClist

def GivenSeqListChResNCASeqReturnStringRef(seqSeqPushRef , listChResNCA , Sequence ) :

    DicChResNRangeSeqVar={}

    DicChResNSeqVar=seqSeqPushRef[1]
    seqvarall=seqSeqPushRef[0]
    chvar=''
    resnvar=''
    for Ch in DicChResNSeqVar:
        #print Ch
        DicChResNRangeSeqVar[Ch]={}
        for ResN in DicChResNSeqVar[Ch]:
            #print ResN
            if (Ch==chvar and ResN-1!=resnvar) or Ch!=chvar:
                if chvar!='': DicChResNRangeSeqVar[Ch][str(ResIvar)+'-'+str(resnvar)]=seqvar
                    #print 'placing dic',str(ResIvar)+'-'+str(resnvar),seqvar
                seqvar=''
                ResIvar=ResN
            #print seqvar
            seqvar += seqvarall[listChResNCA.index([Ch, ResN])]
            chvar = Ch
            resnvar = ResN
    DicChResNRangeSeqVar[Ch][str(ResIvar) + '-' + str(resnvar)] = seqvar

    strvar=''
    for Ch in DicChResNRangeSeqVar:
        for ResNRange , seqvar in DicChResNRangeSeqVar[Ch].items():
            bestindexalignment=BestAlignment2StringsReturnIndex( Sequence , seqvar )
            if strvar!='': strvar += '&'
            strvar += str(bestindexalignment) + '-' + str(bestindexalignment + len(seqvar) - 1)
            #strvar+=Ch+'_'+ResNRange+'_seq'+str(bestindexalignment)+'-'+str(bestindexalignment+len(seqvar)-1)
    return strvar

def GetFreeMemory():
    linesfm=subprocess.check_output(['free', '-m'], text=True)
    linesfm=linesfm.split('\n')
    freemem=float( linesfm[1].split() [linesfm[0].split().index('free')+1] )
    return freemem
    #result = subprocess.run(['free', '-m'], stdout=subprocess.PIPE) #python 3
    #print result.stdout

# def GivenLSQKABlogReturnRMSD (log):
#     with open(log)as f:
#         out=f.read()
#         rmsd=float ([out.index("          RMS     XYZ DISPLACEMENT = ") + 36:out.index(
#             "          RMS     XYZ DISPLACEMENT = ") + 44])
#     return rmsd

def CalculateRSCCfromMAP (inputmtz,inputpdb,outputt,programm='phenix.get_cc_mtz_pdb'):
    l=programm+' raise_sorry=True debug=True fix_xyz=True '+inputmtz+' '+inputpdb+' temp_dir='+outputt
    ll=l.split()
    #os.system(programm+' raise_sorry=True debug=True fix_xyz=True '+inputmtz+' '+inputpdb+' temp_dir='+output + ' > '+output+'_run.log')
    p = subprocess.Popen(ll, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE , text=True)#, cwd=folder)
    out, err = p.communicate()
    with open(outputt+'_run.log','w') as fw: fw.write(out)
    os.system('cp '+outputt+'/cc.log '+outputt+'_RSCC.log')
    shutil.rmtree(outputt)

def ReadRSCCPhenix (inputlog):
    with open(inputlog) as ff: f=ff.read()
    lf=f[f.index('  RESIDUE           CC       ALL     MAIN     SIDE'):].split('\n')
    i=1
    li=lf[i]
    dicc={}
    irscc_overall=f.index('Overall map correlation:')
    rscc_overall=float(f[irscc_overall+25:irscc_overall+32])
    # print rscc_overall
    # exit()
    dicc['all']=rscc_overall
    for l in lf[1:]:
        if l=='': break
        #li=l
        l=l.split()
        ch=l[2]
        resn=int(l[1])
        rscc=float(l[3])
        if ch not in dicc: dicc[ch]={}
        dicc[ch][resn]=rscc
    return dicc

def readshelxeline(lstfile):
    with open(lstfile) as ff: f=ff.readlines()
    ShelxeLine=f[f.index(' Command line parameters:\n')+1].split()[1:]
    return ShelxeLine


# Functions created 14 May 2021

#
def GenerateSym (pdbin,pdbout,dist=5.5,pymolpath='pymol',pymolins=''): #both should be strings
    pdb=pdbin[pdbin.rindex('/')+1:-4]
    #print(pdbin, pdbout, pymolins)
    if pymolins=='': pymolins=pdbout[:-3] + 'pml'
    with open(pymolins, 'w') as fw:
        fw.write('load '+pdbin+'\n')
        fw.write('symexp sym,'+pdb+',('+pdb+'),'+str(dist)+'\n')
        fw.write('select all\n')
        fw.write('save '+pdbout+', sele\n')
        fw.write('quit')
    os.system(pymolpath+' -c '+pymolins+' > /dev/null') #'+pymolins[:-4]+'_pymol.log')

def ChangeChSym(pdbin,pdbout):
    with open(pdbin) as f: fl=f.readlines()
    CheckStr=[]
    with open(pdbout,'w') as fw:
        for l in fl:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                lin=l[11:26]
                if lin in CheckStr:
                    for aa in alphabet:
                        lin=lin[:10]+aa+lin[11:]
                        if lin not in CheckStr: break
                CheckStr.append(lin)
                l=l[:11]+lin+l[26:]
                fw.write(l)


def RemoveCheckResAboveDist(pdbin,chf,resnf,pdboutshort,pdboutlarge,distout,clashout,distint=4.0,distclash=2.5):
    #dist=20.5): # trash 20.5 because it is the maximum distance of a Ca to Nh in Arg (7.5) * 2 + 5.5 (maximum considered interaction)
    resnf=int(resnf)
    with open(pdbin) as f: fl=f.readlines()
    fw=open(distout,'w')
    fwclash=open(clashout,'w')
    fw2 = open(pdboutshort, 'w')
    fw.write     ('Ch1\tResN1\tResT1\tAtType1\tCh2\tResN2\tResT2\tAtType2\tDist\n')
    fwclash.write('Ch1\tResN1\tResT1\tAtType1\tCh2\tResN2\tResT2\tAtType2\tDist\n')
    dicChKeep=defaultdict(list)
    DicRes={}
    for l in fl:
        if l.startswith('ATOM') and l[21]==chf and int(l[22:26])==resnf and l[12:16] not in (' N  ',' C  ',' O  '):#' CA ',
            atomt1=l[12:16]
            rest1=l[17:20]
            pos1=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
            DicRes[atomt1]=pos1
            fw2.write(l)
    for l in fl:
        if l.startswith('ATOM') and (l[21]!=chf or int(l[22:26])!=resnf) or l.startswith('HETATM'):
            atomt2 = l[12:16]
            rest2 = l[17:20]
            ch2=l[21]
            resn2=int(l[22:26])
            pos2=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
            write2=True
            for at1,pos1 in DicRes.items():
                distvar=CalculateEuclidianDistance(pos1,pos2)
                if distvar<=distint:
                    sdistvar='%.3f'%(distvar)
                    if resn2 not in dicChKeep[ch2]: dicChKeep[ch2].append(resn2)
                    if distvar >= distclash: fw.write     (chf+'\t'+str(resnf)+'\t'+rest1+'\t'+at1+'\t'+ch2+'\t'+str(resn2)+'\t'+rest2+'\t'+atomt2+'\t'+sdistvar+'\n')
                    else:                    fwclash.write(chf+'\t'+str(resnf)+'\t'+rest1+'\t'+at1+'\t'+ch2+'\t'+str(resn2)+'\t'+rest2+'\t'+atomt2+'\t'+sdistvar+'\n')
                    if write2:
                        fw2.write(l)
                        write2=False

    #Add previous and later residue to prevent extra addition of a H by reduce
    dicChKeep2 = defaultdict(list)
    for ch in dicChKeep:
        for resn in dicChKeep[ch]:
            for resnv in [resn-1,resn,resn+1]:
                if resnv not in dicChKeep2: dicChKeep2[ch].append(resnv)
    with open(pdboutlarge,'w') as fw3:
        for l in fl:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                ch2=l[21]
                resn2=int(l[22:26])
                if (l[21]==chf and int(l[22:26])==resnf) or resn2 in dicChKeep2[ch2]: fw3.write(l)
    fw.close()
    fw2.close()

def CheckWatersFromSymmetry(pdbin,pdbsym,chf,resnf,dist=2.5):
    resnf=int(resnf)
    with open(pdbin)  as f: fl1=f.readlines()
    with open(pdbsym) as f: fl2=f.readlines()
    DicRes={}
    for l in fl1:
        if l.startswith('ATOM') and l[21]==chf and int(l[22:26])==resnf and l[12:16] not in (' N  ',' C  ',' O  '):#' CA ',
            atomt1=l[12:16]
            rest1=l[17:20]
            pos1=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
            DicRes[atomt1]=pos1
    #print (DicRes)
    LChResnRemove=[]
    for l in fl2:
        if l.startswith('HETATM') and l[17:20]=='HOH':
            ch2=l[21]
            resn2=int(l[22:26])
            pos2=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
            for at1,pos1 in DicRes.items():
                distvar=CalculateEuclidianDistance(pos1,pos2)
                if distvar>15: break
                #print(ch2, resn2, pos2,distvar)
                #print (at1,ch2,resn2,'%.3f'%(distvar))
                if distvar<dist:
                    LChResnRemove.append( (ch2,resn2) )
    return LChResnRemove

def RemoveWat(LChResnWat,pdbin,pdbout):
    with open(pdbin) as f: fl = f.readlines()
    with open(pdbout,'w') as fw:
        for l in fl:
            if l.startswith('HETATM') and l[17:20]=='HOH' and (l[21],int(l[22:26])) in LChResnWat: pass
            else: fw.write(l)

def CheckWatPDBs(pdb1,pdb2):
    with open(pdb1) as f: fl1 = f.readlines()
    with open(pdb2) as f: fl2 = f.readlines()
    linWat1=[]
    for l in fl1:
        if l.startswith('HETATM') and l[17:20] == 'HOH': linWat1.append(l[11:26])
    linWat2=[]
    for l in fl2:
        if l.startswith('HETATM') and l[17:20] == 'HOH': linWat2.append(l[11:26])
    linWat1=set(linWat1)
    linWat2=set(linWat2)
    ldiff=list( linWat1-linWat2 ) + list( linWat2-linWat1 )
    return ldiff

def RetrieveSumClashSSbond(clashin,maxdist=2.4):
    with open(clashin) as f: fl=f.readlines()
    SumClash=0
    nSSb=0
    for l in fl[1:]:
        l=l.split()
        if l[2]=='CYS' and l[6]=='CYS' and l[3]=='SG' and l[7]=='SG': nSSb+=1
        else: SumClash+=maxdist-float(l[-1])
    return SumClash,nSSb

def readPhenixClashscore(log,chf,resnf):
    with open(log) as f: fr=f.readlines()
    clashscoresum=0.0
    check=False
    for l in fr:
        if l=='Bad Clashes >= 0.4 Angstrom:\n':
            check=True
            #print ('1:',l)
        elif check:
            #print ('2:',l)
            if l.startswith('clashscore'): clashscore=float(l.split()[-1])
            else:
                l=l.split()
                ch1,resn1=l[0],int(l[1])
                ch2,resn2=l[4],int(l[5])
                if chf in [ch1,ch2] and resnf in [resn1,resn2]: clashscoresum+=float(l[-1][1:])
    return clashscore,clashscoresum

def ReturnHSaltbonds(login,chf,resnf):
    with open(login) as f: fr = f.readlines()
    nH=0
    nSalt=0
    check=False
    for l in fr:
        if l=='n    s   type  num  typ     dist DA aas  dist angle  dist       angle   num\n': check=True
        elif check:
            lc=[ (l[0],int(l[1:5])), (l[14],int(l[15:19])) ]
            if (chf,resnf) in lc:
                if l[6:12]  in ('LYS NZ','ARG NH') and l[20:26] in ('GLU OE','ASP OD') or l[6:12] in ('GLU OE','ASP OD') and l[20:26] in ('LYS NZ','ARG NH'): nSalt+=1
                else: nH+=1
    return nH,nSalt


def ReturnHydInt(login):
    with open(login) as f: fr = f.readlines()
    nHydInt=0
    for l in fr[1:]:
        l=l.split()
        atomt1=l[3]
        rest1=l[2]
        atomt2=l[7]
        rest2=l[6]
        if rest1 in DicHydrophobic and atomt1 in DicHydrophobic[rest1] and rest2 in DicHydrophobic and atomt2 in DicHydrophobic[rest2]: nHydInt+=1
    return nHydInt

# def ReadRotamerProb(tablein='/home/rborges/DATA/DunbrackRotamerLibrary/Everything-5/SimpleOpt1-5/ALL.bbdep.rotamers.lib'):
#     with open(tablein) as f: fr.readlines()
#     dicaa={}
#     for l in fr:
