#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 15 March 2024
# Author : Rafael Borges
# Objective: Generate amino acid possibilities by residue, fit best rotamer in electron density and calculate their
# side chain and main chain real-space correlation coefficient

# USAGE: SLIDER_VENOM.py pdb_file reflections.mtz map_coeffs.mtz output TRYALL

# pdb_file should be a protein coordinate file with best possible fit against electron density
# reflections.mtz should contain intensities or amplitudes
# map_coeffs.mtz with calculated map (2FOFCWT/PH2FOFCWT)
# output should be a new path, a folder 'output' will contain all SLIDER output and files starting with output will
#     summarize information:
# type_of_run (which was TRYALL in example USAGE)
#
# Type_of_run should be a string with keywords and can change how SLIDER is run:
# Amino acid possibilities can be generated in three different scenarios:
# Trying all 20 possibilities for each residue (keyword: TRYALL)
# Possibilities restricted by:
#  mass spectrometry (keyword: MASSPEC)
#  alignment (keyword: ALIGN)
# If the last two options are chosen, an additional file containing either mass spectrometry or alignment file should be given.
#
# The mass spectrometry file should be a file containing text of amino acids (aa) in one letter code (FASTA), each column should contain generated aa, per example, if 1st and 3rd residue should be a F and T, and 2nd residue either L,D or N, the file should be:
# -------------------------------------------------------------------------
# FLT
#  D
#  N
# -------------------------------------------------------------------------
#
# Keyword SKIPTEST should be given to skip RAM memory calculation;
#
# Links to the external software of Phenix and Coot should be accessible through the terminal.
#
# Output files:
# output_all.log contains information of chain and residue number, RSCC and delta contrast
# output_Resolved.log same as output_all.log, but only resolved residues
# output_Dubious.log  same as output_all.log, but only dubious residues
# output_runline.log has the line used to run SLIDER_VENOM.py
# output_polder_coot_open_maps contains a script to open omit maps in coot (Go to Calculate -> Run Script -> select file)

from __future__ import print_function
from builtins import range
from builtins import bytes, str

import os,sys,time,shutil
import multiprocessing
from collections import defaultdict
import RJB_lib
import numpy
import datetime
import Bio.PDB

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
# from Bio import BiopythonParserWarning
# warnings.simplefilter('ignore', BiopythonParserWarning)
# from Bio import BiopythonDeprecationWarning
# warnings.simplefilter('ignore', BiopythonDeprecationWarning)

import re
import subprocess

pdb = sys.argv[1]
mtz_init = sys.argv[2]
mtz_phases = sys.argv[3]
output_folder= sys.argv[4]
typeee = sys.argv[5]

#typeee may be:
#MASSPEC: use partial MASS SPECTROMETRY DATA
#TRYALL: try all possibilities
##CONSTRUCT: add word CONSTRUCT to construct the model while it runs... Desactivated
#SINGLE: for structures composed of single protein
#ALIGN: alignment file
#MAINCHAIN: keep all atoms in residue for generating phenix.polder CC
#now default: SIDECHAIN: skip C, N, O in generating phenix.polder CC
#SELCH: do not calculate chains given in sys.argv[6] , should be A,B,C or A
#SKIPTEST: will not use RAM memory calculation / TEST
#BRAGG: use 40 cores
#ML: calculation of residue depth, number and type of side chain interactions and clashes

if 'SELCH' in typeee: RemoveChains=sys.argv[6].split(',')
else:                 RemoveChains=False

if 'ML' in typeee: ML=True
else:              ML=False

#clashes check VDW radii (0.4 A):
#PON-SC           https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1947-7

#Molprobity probe https://www.sciencedirect.com/science/article/abs/pii/S0022283698924007?via%3Dihub
#Table 1. Atomic parameters used in Reduce and Probe
#A. Bond lengths (Å)
# blCH=1.1
# blNH_OH=1.0
# blSH=1.3
# #B. Van der Waals radii (Å)
# vdwH=1.17
# vdwHarom=1.0
# vdwHpol=1.0
# vdwC=1.75
# vdwCcarbonyl=1.65
# vdwN=1.55
# vdwO=1.4
# vdwP=1.8
# vdwS=1.8

#chimera https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/findclash/findclash.html
#overlapij = rVDWi + rVDWj – dij – allowanceij
#For detecting clashes, cutoff values of 0.4-1.0 Å and allowance values of 0.2-0.6 Å are generally reasonable (default clash criteria 0.6 and 0.4 Å, respectively).

#hydrogen bond
# 2.2-2.5 Å as “strong, mostly covalent”,
# 2.5-3.2 Å as “moderate, mostly electrostatic”,
# 3.2-4.0 Å as “weak, electrostatic”
#http://biomodel.uah.es/en/water/hbonds.htm
#https://www.researchgate.net/post/What-are-the-standard-values-of-distance-and-angle-cutoffs-to-be-considered-for-hydrogen-bond-formation

#https://en.wikipedia.org/wiki/Van_der_Waals_force
# 4.0-6.0 Å as “weak, electrostatic”

#LigPlus
#https://www.ebi.ac.uk/thornton-srv/software/LigPlus/manual2/manual.html
# hydrogen bond 2.7 Å (H)
# hydrogen bond 3.35 Å (without H)
# non-bonded    2.9  Å (H)
# non-bonded    3.9  Å (without H)

#PLIP
#https://projects.biotec.tu-dresden.de/plip-web/plip/help
# hydrophobic distance   4.0 Å
# hydrogen bond distance 4.1 Å
# pistack distance       7.5 Å
# pi cation distance     6.0 Å
# salt bridge distance   5.5 Å


#https://www.frontiersin.org/articles/10.3389/fmolb.2015.00056/full
#interaction https://www.frontiersin.org/files/Articles/164719/fmolb-02-00056-HTML/image_m/fmolb-02-00056-g001.jpg

#types of interaction: #https://www.cambridgemedchemconsulting.com/resources/molecular_interactions.html
#Typical Energies Salt Bridge ~2 kcal/mol H-Bond ~1 kcal/mol Hydrophobic ~0.7 kcal/mol Aromatic ~1-3 kcal/mol

#https://www.ncbi.nlm.nih.gov/books/NBK21726/#:~:text=There%20are%20four%20main%20types,1%20to%205%20kcal%2Fmol.
#hydrogen bonds in proteins and nucleic acids are only 1 to 2 kcal/mol
#hydrogen bond in water (≈5 kcal/mol)
#van der Waals interaction is about 1 kcal/mol

#Biomolecular Crystallography Bernhard Rupp
#S-S bond     ≈62 kcal/mol
#ionic       3-10 kcal/mol
#Hydrogen    3-7  kcal/mol
#Polar       1-8  kcal/mol
#Van der Waals ≈5 kcal/mol
#Hydrophobic ~0.7 kcal/mol

# I will find only hydrophobic interaction within:
# minimum distance 2.5 Å
# maximum distance 4.0 Å
# use hbplus for hydrogen bond search https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/

#Rotamers
# DicRotamers={'':}

# multiprocessing implementation of phenix.clashscore, done Mar15,2024
def phenixclashscore(clashscore_path, pdb, log):
    p = subprocess.Popen([clashscore_path, pdb], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = p.communicate()
    # print (out)
    with open(log, 'w') as f: f.write(out)

# multiprocessing implementation of phenix.rotalize, done Mar15,2024
def phenixrotalize(rotalyze_path, pdb, log):
    p = subprocess.Popen([rotalyze_path, pdb], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    out, err = p.communicate()
    with open(log, 'w') as f: f.write(out)


print ('\n\n\n')

now = datetime.datetime.now()
date = now.isoformat()[:10] + ' ' + now.isoformat()[11:16]
print ('Initiating '+sys.argv[0]+' '+date)


RJB_lib.output_runline(output_file=output_folder+'_runline.log',printtt=True)

now = datetime.datetime.now()
date = now.isoformat()[:10] + ' ' + now.isoformat()[11:16]
print ('Running '+sys.argv[0]+' '+date)

nproc=RJB_lib.number_of_processor()
#nproc=100
#nproc=20
#nproc=24
#nproc=48
#nproc=6
if 'BRAGG' in typeee: nproc=40

####Check if files exist and if mtz_phases does not, it generates it with phenix.maps
for i in [pdb,mtz_init,mtz_phases]:
    if not os.path.isfile(i):
        if i==mtz_phases:
            #mtz_phases=pdb[:-4] + '_map_coeffs.mtz'
            print ('\nMtz file with 2Fo-Fc map and Fo-Fc not present, generating with phenix.maps with output: '+mtz_phases+'\n')
            os.system('phenix.maps '+pdb+' '+mtz_init+' > '+mtz_phases[:-4]+'_run.log')
            os.system('mv '+pdb[:-4] + '_map_coeffs.mtz '+mtz_phases)
            if not os.path.isfile(mtz_phases):
                print ('Failure to generate file. Exiting program')
                exit()
        else:
            print (i , 'file does not exist, exiting program.')
            exit()


####Checks if pdb_file contains disordered residues, if so, it removes them.
dic_disorder=RJB_lib.remove_pdb_double_occupancy_atoms_pdb (pdb , pdb[:-4]+'_removed_disordered_atoms.pdb')

#added 20181031 to remove partial occupancy atoms from dictionary
if dic_disorder!=False:
    #print 'Removing residues with double occupancy using BioPython:'
    #print dic_disorder
    print ('The following residues with less than 1.0 occupancy atoms will be kept in model, but will be removed from SLIDER_VENOM evaluation:')
    for ch in dic_disorder:
        print ('in chain',ch, ' ,'.join(str(x) for x in dic_disorder[ch]))
    #exit()

amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
amino_acid_list_3L=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

#dfiles=defaultdict(list)


dic_res=RJB_lib.return_dic_resnumb_list (pdb)
TotResN=0
for ch,listres in dic_res.items():
    TotResN+=len(listres)
if dic_disorder:
    for ch,listres in dic_disorder.items():
        TotResN-=len(listres)

dic_pdb=RJB_lib.return_dic_chain_resnumb_restype (pdb)
#delete empty keys (chains with ligs or waters)
deletekeys=[]
for ch,lres in dic_res.items():
    #print ch,lres
    if len(lres)==0 or (RemoveChains!=False and ch in RemoveChains):
        deletekeys.append(ch)
for delet in deletekeys:
    del dic_res[delet]

RJB_lib.mkdir(output_folder)
RJB_lib.mkdir(output_folder+'/'+'mainchain')

dic_pos_aa={}
for ch in dic_res:
    dic_pos_aa[ch]=defaultdict(list)

alig=False

if 'ALIGN' in typeee:
    alig=True
    ali=sys.argv[6]
    dic_ali=RJB_lib.generate_dict_count_alignment (ali)

##cc=0
##for resn,lmut in dic_ali.items():
##        cc+=len(lmut)
##print cc
##exit()

if 'MAINCHAIN' in typeee: sidechainoption=False
else: sidechainoption=True

if 'MASSPEC' in typeee:
    seq=sys.argv[6]
    with open(seq) as f: lseq=f.readlines()

#print(dic_pos_aa)

#for i in dic_res['A']:
#Modification to consider multiple chains in MASSPEC input, Mar14,2024
for ch,lres in dic_res.items():
    #print ('\n')
    for imasspec,ires in enumerate(lres):
        #print (ires)
        if 'MASSPEC' in typeee:
            #EvalResMASSPEC=''
            for line in lseq:
                if line.startswith('>'): chmasspec=line[-2]
                #print (chmasspec)
                elif ch==chmasspec:
                    if imasspec<len(line) and line[imasspec] not in [' ','-','\n']:
                        #print ('here',line[imasspec])
                        #exit()
                        if line[imasspec] in dic_pos_aa[ch][ires]:
                            print ('Repeated aminoacid',line[imasspec-1],'found in the file\n',seq)
                            print ('Please correct it.')
                            exit()
                        dic_pos_aa[ch][ires].append(line[imasspec])
                        #EvalResMASSPEC+=line[imasspec]
            # if EvalResMASSPEC=='':
            #     for a in amino_acid_list:
            #         dic_pos_aa[ch][i].append(a)
        elif 'TRYALL' in typeee:
            for a in amino_acid_list:
                dic_pos_aa[ch][ires].append(a)

#print(dic_pos_aa)



####CONVERT SINGLE CHAIN TO ALL CHAINS
##chall=''
##lch=[]
##for ch in dic_pos_aa:
##    print ch
##    chall+=ch
##    lch.append(ch)
##dic_pos_aa[chall]=dic_pos_aa['A']
##for delet in lch:
##    del dic_pos_aa[delet]
##for ch in dic_pos_aa:
##    print ch

if 'MASSPEC' not in typeee and 'TRYALL' not in typeee and 'ALIGN' not in typeee:
    print ('Failure. Wrong option for typeee.')
    exit()


if 'SINGLE' in typeee:
    newkey=''
    for ch in dic_pos_aa:
        newkey+=ch
    newd=dic_pos_aa[ch]
    dic_pos_aa={newkey:newd}

if 'ALIGN' in typeee and not 'MASSPEC' in typeee and not 'TRYALL' in typeee:
    for ch in dic_pos_aa:
        for i in dic_ali:
            for a in dic_ali[i]:
                dic_pos_aa[ch][i].append(a)



#check if coot script mutate is going to work:
##dchkcoot={}
##lseq=[]
##for ch,d1 in dic_pdb.items():
##    sseq=''
##    dchkcoot[ch]={}
##    for resn,d2 in d1.items():
##        if not d2=='A':
##            dchkcoot[ch][resn]='A'
##            sseq+='A'
##        else:
##            dchkcoot[ch][resn]='T'
##            sseq+='T'
##    lseq.append(sseq)
###print dchkcoot
##RJB_lib.coot_run_rotamer_sphere_refinement_multi ( input_PDB_file=pdb , input_mtz_file=mtz_phases , output_pdb='testcootmutatefn.pdb' , outputCootPy='testcootmutatefn.py' , dic=dchkcoot , radius_sph_ref=0  , coot_path=False )
###chcheckcoot,seqcheckcoot,tnrescheckcoot=RJB_lib.extract_protein_chainID_res_number('testcootmutatefn.pdb')
##seqcheckcoot=RJB_lib.return_restype_list ('testcootmutatefn.pdb')
##checkcoot=True
##print 'seqcheckcoot',seqcheckcoot
##print 'lseq        ',lseq
##for i,s in enumerate (seqcheckcoot):
##    for ii,aa in enumerate(s):
##        if aa!=lseq[i][ii]:
##            print 'COOT MUTATE FUNCTION DID NOT WORK WITH RESIDUE NUMBER',ii+1,'!!! By my testing, I would say that the residue has not all main chain atoms partially occupied, this should be the reason. Please make sure to change until coot mutate works and change PDB input file.'
##            checkcoot=False
##if not checkcoot:
##    print 'Exitting now'
##    exit()




#all resn that only has one possibility will be first integrated into model
#if 'MASSPEC' in typeee:
dic_seq_initial=RJB_lib.return_dic_sequence(pdb)
d_imp={}
impr=False
listchains=[]
for ch,dires in dic_pos_aa.items():
    listchains.append(ch)
    for i,lmut in dires.items():
        if len(lmut)==1 and dic_seq_initial[ch][ dic_res[ch].index(i) ]!=lmut[0]:
            try:
                d_imp[ch][i]=lmut[0]
                impr=True
            except:
                d_imp[ch]={}
                d_imp[ch][i]=lmut[0]
                impr=True

# if impr:
#     print 'Improving initial model to model residues that were given only one possibility. If the residue were already a component of that position in the initial model, then no improvement will be done in that specific residue.'
#     print 'Initial PDB:',pdb,'will be saved to',output_folder+'/'+pdb[:-4]+'_imp_single_possibilities.pdb'
#     if not os.path.isfile(output_folder+'/'+pdb[:-4]+'_imp_single_possibilities.pdb'):
#         RJB_lib.coot_run_rotamer_sphere_refinement_multi ( input_PDB_file=pdb , input_mtz_file=mtz_phases , output_pdb=output_folder+'/'+pdb[:-4]+'_imp_single_possibilities.pdb' , outputCootPy=output_folder+'/'+pdb[:-4]+'_imp_single_possibilities.py' , dic=d_imp , radius_sph_ref=500  , removeWAT=False, coot_path=False )
#     pdb=output_folder+'/'+pdb[:-4]+'_imp_single_possibilities.pdb'
#            #coot_run_rotamer_sphere_refinement_multi ( input_PDB_file , input_mtz_file , output_pdb , outputCootPy , dic , radius_sph_ref=5  , coot_path=False )

####RUNNING AREAIMOL FOR SIDE CHAIN ATOMS
if not os.path.isfile (pdb[:-4]+'-areaimol.pdb'): RJB_lib.runAREAIMOLccp4_bigMolec (pdbfile=pdb,outpdb=pdb[:-4]+'-areaimol.pdb',listchains=listchains)


#added 20181031 to remove partial occupancy atoms from dictionary
if dic_disorder!=False:
    for ch in dic_disorder:
        for resn in dic_disorder[ch]:
            if resn in dic_pos_aa[ch]: del dic_pos_aa[ch][resn]




################added 20181109 to test maximum RAM memory usage
if 'SKIPTEST' in typeee:
    print ('Test of RAM memory usage for coot and polder jobs was selected to be skipped.')
    NewNProcCoot=nproc
    NewNProcPolder=nproc
else:
    ch=list(dic_pos_aa.keys())[0]
    dires=dic_pos_aa[ch]
    resn=list(dic_pos_aa[ch].keys())[0]
    stresn=str(resn)
    #a=lmut[0]
    a='A'
    # print 'Chain',ch
    # print 'Residue number',resn
    # print 'Residue type',a

    print ('\nTesting RAM Memory usage of coot & phenix.polder processes with chain',ch,'Residue number',resn,'Residue type',a)

    counterr=0
    ifmem=RJB_lib.GetFreeMemory()
    for i in range(9):
        varmem=float( RJB_lib.GetFreeMemory() )
        if varmem<ifmem: ifmem=varmem
        time.sleep(0.1)

    print ('\nFree Memory before runnning external programs',ifmem,'Mb')


    ####RUNNING COOT MODELING
    aa = amino_acid_list_3L[amino_acid_list.index(a)]
    # COOT section
    outf = output_folder + '/test-mem-coot-' + stresn + a + ''
    dic = {}
    for c in ch:
        dic[c] = {resn: a}
    # print dic


    # print 'Running coot mutate, rotamer search and sphere refine in chain',ch,'and residue',stresn,'with side chain',aa,'from file:',pdb,'to save in:',outf+'.pdb'
    process = multiprocessing.Process(target=RJB_lib.coot_mutate_sph_ref_correcting_files, args=(pdb, mtz_phases, outf + '.pdb', dic, dic_pdb, 5, True, False))
    # coot_mutate_sph_ref_correcting_files          (pdb,mtz_phases,outpdb,dic,dic_pdb radius_sph_ref=5, coot_path=False , printtt=False)
    process.start()
    vcootfmem=RJB_lib.GetFreeMemory()
    #time.sleep(0.1)
    # break

    #wait coot jobs finish
    while len(multiprocessing.active_children()) != 0:
        varmem=float( RJB_lib.GetFreeMemory() )
        if varmem<vcootfmem: vcootfmem=varmem
        #print RJB_lib.GetFreeMemory()
        time.sleep(0.1)
        #break


    print ('\nFree memory running 1 process of coot =',vcootfmem)
    cootfmem=ifmem-vcootfmem
    print ('Maximum RAM memory spent on coot',cootfmem,'Mb')

    #print '\nNow with phenix.polder'
    ####RUNNING PHENIX.POLDER
    #polder section
    pdbinput=str(outf)+'.pdb'
    outf = str ( output_folder + '/test-mem-polder' + stresn + a )
    dic={ch:{resn:a}}
    #print dic
    #process = multiprocessing.Process(target= RJB_lib.coot_run_rotamer_sphere_refinement_multi , args= ( pdb=pdb,mtz=mtz,output=outf,ch=ch,resn=resn,rest3L_mutate=aa , coot_path=False , polder_path=False ) )
    output_map=True
    # while 1:
    process = multiprocessing.Process(target= RJB_lib.run_phenix_polder_multiproc , args= ( pdbinput , mtz_init , dic , outf+'.log' ,  False , output_map , sidechainoption) )
                                                     #run_phenix_polder_multiproc ( pdb , mtz_init , dic , output , polder_path=polder_path, output_mtz=True) #dic should be dic[chain][resn]
    print ('Running phenix.polder to calculate CC of refine in chain(s)',ch,'and residue',stresn,'with side chain',aa,'from file:',pdbinput,'to save in:',outf+'.log\n')
    process.start()
    vpolderfmem=RJB_lib.GetFreeMemory()

    #wait phenix.polder jobs finish
    while len(multiprocessing.active_children()) != 0:
        varmem=float( RJB_lib.GetFreeMemory() )
        if varmem<vpolderfmem: vpolderfmem=varmem
        time.sleep(0.1)

    print ('\nFree memory running 1 process of polder =',vpolderfmem)
    polderfmem=ifmem-vpolderfmem
    print ('Maximum RAM memory spent on polder',polderfmem,'Mb')


    # '\nIt was found',nproc,'processors (counting HyperThreading if available).'
    # '600 Mb will be left free as tolerance and to be free for other programs.'
    if cootfmem<1.0: cootfmem=1.0
    if polderfmem < 1.0: polderfmem = 1.0
    NewNProcCoot=   int ( (ifmem-600)/cootfmem   )
    NewNProcPolder= int ( (ifmem-600)/polderfmem )
    print ('There is',ifmem-600,'available RAM memory')
    print ('To be divided to',cootfmem,'for coot jobs')
    print ('To be divided to',polderfmem,'for polder jobs')

    if NewNProcCoot==0 or NewNProcPolder==0:
        print ('Run requires more RAM memory than available.')
        exit()
    else:
        if NewNProcCoot   > nproc: NewNProcCoot   = int(nproc)
        if NewNProcPolder > nproc: NewNProcPolder = int(nproc)

    print ('\nTherefore, it will be used:')
    print (NewNProcCoot  ,'processors for coot jobs.')
    print (NewNProcPolder,'processors for polder jobs.')

################added 20181109 to test maximum RAM memory usage END

#As coot now uses multiple CPUs, I am reducing by 4 the #CPU processors running coot.
NewNProcCoot=int(NewNProcCoot/4)
print ('As coot now uses multiple CPUs, I am reducing by 4 the #CPU processors running coot.')
print (NewNProcCoot  ,'processors for coot jobs.')

#Run Residue Depth from https://biopython.org/docs/1.75/api/Bio.PDB.ResidueDepth.html
if ML:
    if not os.path.isfile(output_folder+'_ResDepth.log'):
        parser = Bio.PDB.PDBParser()
        structure = parser.get_structure("ResDepth", pdb)
        model = structure[0]
        rd = Bio.PDB.ResidueDepth(model)
        resdepth=open(output_folder+'_ResDepth.log','w')
        resdepth.write('Ch\tResN\tResDepth')
        for ch,dires in dic_pos_aa.items():
            for resn,lmut in dires.items():
                #print(ch,resn,lmut)
                resdepthvar=rd[ch, (' ', resn, ' ')]
                resdepthvar='%.1f'%(resdepthvar[1])
                resdepth.write('\n'+ch+'\t'+str(resn)+'\t'+resdepthvar)

RJB_lib.mkdir(output_folder+'/eval')
for ch,dires in dic_pos_aa.items():
    #print ch
    RJB_lib.mkdir(output_folder+'/'+ch)
    RJB_lib.mkdir(output_folder + '/eval/'+ch)
    for resn,lmut in dires.items():
        print ('\nEvaluating chain',ch,'and residue',resn,' with: ',', '.join(lmut),' (total: ',str(len(lmut))+')')
        print ('Running coot mutate, rotamer search and sphere refine')
        #print lmut
        #print resn
        stresn=str(resn)
        outfold=output_folder+'/'+ch+'/'+stresn+'/'
        RJB_lib.mkdir(outfold)
        if (not os.path.isfile(outfold+ch+'_'+stresn+'_polder.log')):
            ####RUNNING COOT MODELING
            #print (lmut)
            for a in lmut:
                #print a
                aa=amino_acid_list_3L[amino_acid_list.index(a)]
                #COOT section
                outf=outfold+stresn+a+''
                dic={}
                for c in ch:
                    dic[c]={resn:a}
                #print dic
                if not os.path.isfile( outf+'.pdb' ):
                    if  NewNProcCoot > -1: #NOTE: PROCESSES es el numero de cores que quieres lanzar, default == numero de cores-1
                #                    print "I found ", sym.REALPROCESSES, "CPUs." #NOTE: REALPROCESSES es el numero de cores de tu ordenador
                        while 1:
                            time.sleep(0.1)
                            if len(multiprocessing.active_children()) < NewNProcCoot:
                                #print 'Running coot mutate, rotamer search and sphere refine in chain',ch,'and residue',stresn,'with side chain',aa,'from file:',pdb,'to save in:',outf+'.pdb'
                                process = multiprocessing.Process(target= RJB_lib.coot_mutate_sph_ref_correcting_files , args= ( pdb , mtz_phases , outf+'.pdb' , dic ,dic_pdb, 5, False , False ) )
                                                                                 #coot_mutate_sph_ref_correcting_files         (pdb,mtz_phases,outpdb,dic,dic_pdb=False, radius_sph_ref=5, removeWAT=False , coot_path=False , printtt=False)
                                process.start()
                                time.sleep(0.1)
                                break
                    else:
                        print ("FATAL ERROR: I cannot load correctly information of CPUs.")
                        exit()

# wait coot jobs finish
while 1:
    time.sleep(0.1)
    if len(multiprocessing.active_children()) == 0:
        break

if ML:
    for ch,dires in dic_pos_aa.items():
        for resn,lmut in dires.items():
            #Evaluation of clashes and interactions (hydrogen / SS-bond / hydrophobic) and their energy
            stresn = str(resn)
            outfold = output_folder + '/' + ch + '/' + stresn + '/'
            outfold2=output_folder + '/eval/'+ch+'/'+stresn
            #print ('Generating files in',outfold2)
            RJB_lib.mkdir(outfold2)
            # fwclashes.write(str(resn))
            #if not os.path.isfile(outfold2+'/'+ch+'_'+stresn+'_clash.log') and not os.path.isfile(outfold2+'/'+ch+'_'+stresn+'_nInt.log'):
            if True:
                print ('Calculating clashes and interactions in folder',outfold2)
            # Found a water molecule (from symmetry) laying in top of 2B04/A1E, need to remove waters
                for a in lmut:
                    #print(ch, stresn,a)
                    outi = output_folder+'/'+ch+'/'+stresn+'/'+stresn+a+'.pdb'
                    outf = outfold2 +'/'+ stresn + a + '.pdb'
                    # t1 = time.time()
                    if not os.path.isfile(outf):
                        #RJB_lib.GenerateSym(pdbin=outi,pdbout=outf,dist=5.5)
                        RJB_lib.GenerateSymKeepNearRes(pdbin=outi,pdbout=outf, ch=ch, NRes=stresn, dist=5.5, pymolpath='pymol', pymolins='')
                    # t2 = time.time()
                    #print (a)
                        LRemoveWatChResnVar=RJB_lib.CheckWatersFromSymmetry(pdbin=outi,pdbsym=outf,chf=ch,resnf=resn,dist=2.5)
                        RJB_lib.RemoveWat(LChResnWat=LRemoveWatChResnVar,pdbin=outf,pdbout=outf)
                    #print (a,LRemoveWatChResnVar)
                    #ldiff=RJB_lib.CheckWatPDBs(pdb1=output_folder+'/'+ch+'/'+stresn+'/'+stresn+a+'.pdb2',pdb2=output_folder+'/'+ch+'/'+stresn+'/'+stresn+a+'.pdb')
                    #if ldiff!=[]: print (resn,a,ldiff)
                #exit()

                #t2 = time.time()
                for a in lmut:
                    outi = output_folder+'/'+ch+'/'+stresn+'/'+stresn+a+'.pdb'
                    outf = outfold2 +'/'+ stresn + a + '.pdb'
                    RJB_lib.ChangeChSym(pdbin=outf,pdbout=outf)
                    outfhbplus=outf[:-4]+'_4hbplus.pdb'
                    loghbplus=outfhbplus[:-3]+'hb2'
                    outfphenixclash=outf[:-4]+'_4phenix.clash.pdb'
                    logfphenixclash=outfphenixclash[:-3]+'log'
                    logfclash=outf[:-4]+'_clash.log'
                    #print('checking files ' + outfhbplus + ' ' + outfphenixclash + ' ' + outf[:-4] + '_dist.log ' + logfclash)
                    if not os.path.isfile(outfhbplus) and not os.path.isfile(outfphenixclash) and not os.path.isfile(outf[:-4]+'_dist.log') and not os.path.isfile(logfclash):
                        #print('checking files '+outfhbplus+' '+outfphenixclash+' '+outf[:-4]+'_dist.log'+logfclash)
                        RJB_lib.RemoveCheckResAboveDist(pdbin=outf,chf=ch,resnf=resn,pdboutshort=outfhbplus,pdboutlarge=outfphenixclash,distout=outf[:-4]+'_dist.log',clashout=logfclash,distint=4.0,distclash=2.4)

                    # Single processor
                    # if not os.path.isfile(logfphenixclash):  # os.system('phenix.rotalyze '+outfold2+res+'.pdb > '+outfold2+res+'.log')
                    #     os.system('phenix.clashscore ' + outfphenixclash + ' > ' + logfphenixclash)

                    #Multiprocessing
                    # print('Checking files ' + outfphenixclash + ' ' + logfphenixclash)
                    # print('Checking file1 '+outfphenixclash , os.path.isfile(outfphenixclash))
                    # print('Checking file2 ' + logfphenixclash , os.path.isfile(logfphenixclash))
                    if not os.path.isfile(logfphenixclash):
                        # print('Checking file2.2 ' + logfphenixclash, os.path.isfile(logfphenixclash))
                        if nproc > -1:  # NOTE: PROCESSES es el numero de cores que quieres lanzar, default == numero de cores-1
                            #                    print "I found ", sym.REALPROCESSES, "CPUs." #NOTE: REALPROCESSES es el numero de cores de tu ordenador
                            while 1:
                                time.sleep(0.1)
                                if len(multiprocessing.active_children()) < nproc:
                                    #print('Running phenix.clashscore', outfphenixclash,logfphenixclash)

                                    process = multiprocessing.Process(target=phenixclashscore, args=(
                                    'phenix.clashscore', outfphenixclash,logfphenixclash))
                                    process.start()
                                    time.sleep(0.1)
                                    break
                        else:
                            print("FATAL ERROR: I cannot load correctly information of CPUs.")
                            exit()
                    # exit()
    while 1:
        time.sleep(0.1)
        if len(multiprocessing.active_children()) == 0:
            break
                #t3 = time.time()

                # for a in lmut:
                #     outfphenixclash=outf[:-4]+'_4phenix.clash.pdb'
                #     logfphenixclash=outfphenixclash[:-3]+'log'
                #
                #     if not os.path.isfile(logfphenixclash):
                #         os.system('phenix.clashscore ' + outfphenixclash + ' > '+logfphenixclash)

                # for a in lmut: fwclash.write()

for ch, dires in dic_pos_aa.items():
    for resn,lmut in dires.items():
        stresn=str(resn)
        if ML:
            outfold2=output_folder + '/eval/'+ch+'/'+stresn

            fwclash = open(outfold2 + '/' + ch + '_' + stresn + '_clash.log', 'w')
            fwclash.write('Aa\tSClDist\tPhClSum\tPhClSc')
            fwInt = open(outfold2 + '/' + ch + '_' + stresn + '_nInt.log', 'w')
            fwInt.write('Aa\tnSSb\tnHb\tnSalt\tnHydInt\tEnergy')

            for a in lmut:
                outi = output_folder+'/'+ch+'/'+stresn+'/'+stresn+a+'.pdb'
                outf = outfold2 +'/'+ stresn + a + '.pdb'
                # t2 = time.time()
                RJB_lib.ChangeChSym(pdbin=outf,pdbout=outf)
                # t3 = time.time()
                outfhbplus=outf[:-4]+'_4hbplus.pdb'
                loghbplus=outfhbplus[:-3]+'hb2'
                outfphenixclash=outf[:-4]+'_4phenix.clash.pdb'
                logfphenixclash=outfphenixclash[:-3]+'log'
                logfclash=outf[:-4]+'_clash.log'
                # if not os.path.isfile(outfhbplus) and not os.path.isfile(outfphenixclash) and not os.path.isfile(outf[:-4]+'_dist.log') and not os.path.isfile(logfclash):
                #     RJB_lib.RemoveCheckResAboveDist(pdbin=outf,chf=ch,resnf=resn,pdboutshort=outfhbplus,pdboutlarge=outfphenixclash,distout=outf[:-4]+'_dist.log',clashout=logfclash,distint=4.0,distclash=2.4)
                SumClashDist,nSSb=RJB_lib.RetrieveSumClashSSbond(clashin=logfclash,maxdist=2.5)
                # if not os.path.isfile(logfphenixclash):
                #     os.system('phenix.clashscore ' + outfphenixclash + ' > '+logfphenixclash)
                #print (logfphenixclash)
                clashscore,clashscoresum=RJB_lib.readPhenixClashscore(log=logfphenixclash,chf=ch,resnf=resn)
                SumClashDist,clashscoresum,clashscore='%.1f'%(SumClashDist),'%.1f'%(clashscoresum),'%.1f'%(clashscore)
                #print (a,SumClashDist,clashscoresum,clashscore )
                fwclash.write('\n'+a+'\t'+SumClashDist+'\t'+clashscoresum+'\t'+clashscore)
                # print (clashscore,clashscoresum)

                #os.remove(outf)
                # t4 = time.time()
                # print('t of GenerateSym = ',t2-t1)
                # print('t of ChangeChSym = ', t3 - t2)
                # print('t of RemoveCheckResAboveDist = ', t4 - t3)
                if not os.path.isfile(loghbplus):
                    os.system('hbplus '+outfhbplus+' > /dev/null')
                    # print (loghbplus[loghbplus.rindex('/')+1:],loghbplus)
                    shutil.move(loghbplus[loghbplus.rindex('/')+1:],loghbplus)
                nH,nSalt=RJB_lib.ReturnHSaltbonds(login=loghbplus,chf=ch,resnf=resn)
                #print (loghbplus,nH)
                nHydInt=RJB_lib.ReturnHydInt(login=outf[:-4]+'_dist.log')
                energySC=62*nSSb+(nH+nSalt)*3+nHydInt*0.7
                senergySC='%.1f'%(energySC)
                fwInt.write('\n'+a+'\t'+str(nSSb)+'\t'+str(nH)+'\t'+str(nSalt)+'\t'+str(nHydInt)+'\t'+senergySC) #'Aa\tnSSb\tnHb\tnSalt\tnHydInt\tEnergy')
                #exit()
                #RJB_lib.Run_hbplus(pdbin=outi,hbplus='/home/rborges/LigPlus/lib/exe_linux64/hbplus')
                # print (outf)
                # exit()
            #exit()
            fwclash.close()
            fwInt.close()

            #t4 = time.time()
            # print('running All evalutation clashscore', t3 - t2)
            # print('running Post evalutation clashscore', t4 - t3)

if ML:
    #Join information from clashes
    with open(output_folder+'_clashes.log','w') as fwclash:
        fwclash.write('Ch\tResN\tResT\tSClDist\tPhClSum\tPhClSc\n\n')
        for ch,dires in dic_pos_aa.items():
            for resn in dires:
                stresn=str(resn)
                outfold2 = output_folder + '/eval/' + ch + '/' + stresn
                with open(outfold2+'/'+ ch+'_'+stresn + '_clash.log') as f: fr=f.readlines()
                for l in fr[1:]: fwclash.write(ch+'\t'+stresn+'\t'+l)
                fwclash.write('\n\n')

    #Join information from interaction
    with open(output_folder+'_interaction.log','w') as fwInt:
        fwInt.write('Ch\tResN\tResT\tnSSb\tnHb\tnSalt\tnHydInt\tEnergy\n\n')
        for ch,dires in dic_pos_aa.items():
            for resn in dires:
                stresn=str(resn)
                outfold2 = output_folder + '/eval/' + ch + '/' + stresn
                with open(outfold2+'/'+ ch+'_'+stresn + '_nInt.log') as f: fr=f.readlines()
                for l in fr[1:]: fwInt.write(ch+'\t'+stresn+'\t'+l)
                fwInt.write('\n\n')



    # ROTAMER EVALUATION
    DicRotResults={}         # dic that contains rotamer results [ch][resn][restype]=Favored/Allowed/OUTLIER
    for res in ['C',  'D',  'E',  'F',  'H',  'I',  'K',  'L',  'M',  'N',  'P',  'Q',  'R',  'S',  'T',  'V',  'W',  'Y'  ]:
        fwrot=open(output_folder+'_rotamers.log','w')
        #output labels in output_folder+'_rotamers.log'
        fwrot.write('Ch\tResN')
        for r in ['C', 'D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']: fwrot.write('\t'+r)
        DicRot=defaultdict(list) # dic that contains what residues were evaluated
        for ch,dires in dic_pos_aa.items():
            if ch not in DicRotResults: DicRotResults[ch]={}
            for resn, lmut in dires.items():
                if res in lmut: DicRot[ch].append(resn)
        outfold2 = output_folder + '/eval/'
        countatom=1
        #for each amino acid type, output pdb file containing all of residues within this aa type
        if not os.path.isfile(outfold2+res+'.pdb'):
            with open(outfold2+res+'.pdb','w') as fwpdb:
                for ch,lresn in DicRot.items():
                    for resn in lresn:
                        stresn=str(resn)
                        with open(output_folder+'/'+ch+'/'+stresn+'/'+stresn+res+'.pdb') as f: fr=f.readlines()
                        for l in fr:
                            if l.startswith('ATOM') and l[21]==ch and int(l[22:26])==resn:
                                #correct atom number
                                scountatom=str(countatom)
                                l='ATOM'+' '*(7-len(scountatom))+scountatom+l[11:]
                                if l[72]==' ': l=l[:72]+l[21]+l[73:]
                                fwpdb.write(l)
                                countatom+=1
        #print (outfold2+res+'.pdb')
        #print ('phenix.rotalyze '+outfold2+res+'.pdb > '+outfold2+res+'.log')
        #exit()

    for res in ['C', 'D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
        if not os.path.isfile(outfold2+res+'.log'): #os.system('phenix.rotalyze '+outfold2+res+'.pdb > '+outfold2+res+'.log')
            if nproc > -1:  # NOTE: PROCESSES es el numero de cores que quieres lanzar, default == numero de cores-1
                #                    print "I found ", sym.REALPROCESSES, "CPUs." #NOTE: REALPROCESSES es el numero de cores de tu ordenador
                while 1:
                    time.sleep(0.1)
                    if len(multiprocessing.active_children()) < nproc:
                        print ('Running phenix.rotalyze',outfold2+res+'.pdb > ',outfold2+res+'.log')
                        process = multiprocessing.Process(target=phenixrotalize, args=('phenix.rotalyze', outfold2+res+'.pdb',outfold2+res+'.log'))
                        process.start()
                        time.sleep(0.1)
                        break
            else:
                print("FATAL ERROR: I cannot load correctly information of CPUs.")
                exit()
    while 1:
        time.sleep(0.1)
        if len(multiprocessing.active_children()) == 0:
            break

    for res in ['C', 'D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
        with open(outfold2+res+'.log') as f: fr=f.readlines()
        for l in fr[1:-1]:
            ch=l[1]
            resn=int(l[2:6])
            l=l.split(':')
            rottype=l[-1]
            rotval=l[-2]
            if resn not in DicRotResults[ch]: DicRotResults[ch][resn]={}
            DicRotResults[ch][resn][res]=rotval
    for ch,dicresnaarot in DicRotResults.items():
        for resn,dicaarot in dicresnaarot.items():
            fwrot.write('\n'+ch+'\t'+str(resn))
            for aa,rotval in dicaarot.items():
                fwrot.write('\t'+rotval)



for ch, dires in dic_pos_aa.items():
    #print (ch)
    RJB_lib.mkdir(output_folder + '/' + ch)
    for resn, lmut in dires.items():
        #print('resn', resn)
        #print('lmut', lmut)
        stresn=str(resn)
        outfmc=output_folder + '/mainchain/' + ch + stresn + '.pdb'
        if not os.path.isfile(outfmc):
            dic = {}
            for c1 in ch:
                #print ('c1',c1)
                for c in c1:
                    #print ('c',c)
                    #print ('dic_pdb[c][resn]',dic_pdb[c][resn])
                    dic[c] = {resn:dic_pdb[c][resn]}
            RJB_lib.remove_Bfactor_occ_res_pdb(pdb_input=pdb,pdb_output=outfmc,dic_ch_resn=dic,AtomsExclude=[],AtomsInclude=['CA','N ','C ','O '])
                                                                                                            # I was writing all atoms names to be excluded,but decided to do the opposite
                                                                                                            # ['CB','CD','CE','CG','CZ','ND','NE','NZ','OD','SD','SG'] )
                                                                                                            #['CB ' ,'SG ' ,'CG1','OD1','OD2','CG ','CD1','CD2','CE1','CE2',        'CZ','ND1','NE2','CG1','CG2',     'SD ','CE ','','','','','','',''] )
                                                                                                            #A/ALA,C/CYS,D/ASP,            F/PHE                           G/GLY0  H/HIS        I/ILE           L/LEU M/MET      N/ASN

#wait coot jobs finish
while 1:
    time.sleep(0.1)
    if len(multiprocessing.active_children()) == 0 :
        break

####RUNNING PHENIX.POLDER
for ch, dires in dic_pos_aa.items():
    # print ch
    for resn, lmut in dires.items():
        print ('\nEvaluating chain', ch, 'and residue', resn, ' with: ', ', '.join(lmut), ' (total: ', str(len(lmut)) + ')')
        print ('Running phenix.polder to calculate CC of refine in chain(s)\n')  # print lmut
        # print resn
        stresn = str(resn)
        outfold = output_folder + '/' + ch + '/' + stresn + '/'
        RJB_lib.mkdir(outfold)
        if (not os.path.isfile(outfold + ch + '_' + stresn + '_polder.log')):
            for a in lmut:
                #print a
                aa=amino_acid_list_3L[amino_acid_list.index(a)]
                #polder section
                outf=outfold+stresn+a+''
                dic={ch:{resn:a}}
                pdbinput=outfold+stresn+a+'.pdb'
                #print dic
                if not os.path.isfile( outf+'.log' ):
                    if  NewNProcPolder > -1: #NOTE: PROCESSES es el numero de cores que quieres lanzar, default == numero de cores-1
                #                    print "I found ", sym.REALPROCESSES, "CPUs." #NOTE: REALPROCESSES es el numero de cores de tu ordenador
                        while 1:
                            time.sleep(0.1)
                            if len(multiprocessing.active_children()) < NewNProcPolder:
                                #process = multiprocessing.Process(target= RJB_lib.coot_run_rotamer_sphere_refinement_multi , args= ( pdb=pdb,mtz=mtz,output=outf,ch=ch,resn=resn,rest3L_mutate=aa , coot_path=False , polder_path=False ) )
                                if 'A' in lmut:
                                    output_map=False
                                    if a=='A':
                                        output_map=True
                                else:
                                    output_map=True
                                process = multiprocessing.Process(target= RJB_lib.run_phenix_polder_multiproc , args= ( pdbinput , mtz_init , dic , outf+'.log' ,  False , output_map , sidechainoption) )
                                                                                 #run_phenix_polder_multiproc ( pdb , mtz_init , dic , output , polder_path=polder_path, output_mtz=True) #dic should be dic[chain][resn]
                                #print 'Running phenix.polder to calculate CC of refine in chain(s)',ch,'and residue',stresn,'with side chain',aa,'from file:',pdbinput,'to save in:',outf+'.log\n'
                                process.start()
                                time.sleep(0.1)
                                break
                    else:
                        print ("FATAL ERROR: I cannot load correctly information of CPUs.")
                        exit()





####RUNNING PHENIX.POLDER MAIN CHAIN
for ch, dires in dic_pos_aa.items():
    # print ch
    for resn, lmut in dires.items():
        print ('Running phenix.polder to calculate RSCC of main chain atoms of residue',resn,'in chain(s)', ch)
        #print '\n'  # print lmut
        # print resn
        stresn = str(resn)
        outfmc = output_folder + '/mainchain/' + ch + stresn + '_polder.log'
        if (not os.path.isfile(outfmc)):
                #polder section
                dic={ch:{resn:a}}
                pdbinput=outfmc[:-11]+'.pdb'
                #print dic
                if not os.path.isfile( outfmc ):
                    if  NewNProcPolder > -1: #NOTE: PROCESSES es el numero de cores que quieres lanzar, default == numero de cores-1
                #                    print "I found ", sym.REALPROCESSES, "CPUs." #NOTE: REALPROCESSES es el numero de cores de tu ordenador
                        while 1:
                            time.sleep(0.1)
                            if len(multiprocessing.active_children()) < NewNProcPolder:
                                #process = multiprocessing.Process(target= RJB_lib.coot_run_rotamer_sphere_refinement_multi , args= ( pdb=pdb,mtz=mtz,output=outf,ch=ch,resn=resn,rest3L_mutate=aa , coot_path=False , polder_path=False ) )
                                output_map=True
                                process = multiprocessing.Process(target= RJB_lib.run_phenix_polder_multiproc , args= ( pdbinput , mtz_init , dic , outfmc ,  False , output_map , 'OnlyMainChain') )
                                                                                 #run_phenix_polder_multiproc ( pdb , mtz_init , dic , output , polder_path=polder_path, output_mtz=True) #dic should be dic[chain][resn]
                                #print 'Running phenix.polder to calculate CC of refine in chain(s)',ch,'and residue',stresn,'with side chain',aa,'from file:',pdbinput,'to save in:',outf+'.log\n'
                                process.start()
                                time.sleep(0.1)
                                break
                    else:
                        print ("FATAL ERROR: I cannot load correctly information of CPUs.")
                        exit()


while 1:
    time.sleep(0.1)
    if len(multiprocessing.active_children()) == 0:
        break

#Dictionary of CC1/3
dicallSC={}

print ('Building Summary of polder RSCC (cc13) of side chain atoms')
###SUMMARY RESULTS PHENIX.POLDER BY CHAIN AND BY RESIDUE NUMBER
quitt = False
for ch, dires in dic_pos_aa.items():
    dicallSC[ch]={}
    # print ch
    for resn, lmut in dires.items():
        dicallSC[ch][resn]=[]
        stresn=str(resn)
        outfold=output_folder+'/'+ch+'/'+stresn+'/'
        #if not os.path.isfile(outfold+ch+'_'+stresn+'_polder.log'):
        ou2=open(outfold+ch+'_'+stresn+'_polder.log','w')
        ou2.write('Residue\tCC1,3\tR\tRfree\tRimp\tRfImp\tRimp')
        lisel=[]
        for a in lmut:
            aa=amino_acid_list_3L[amino_acid_list.index(a)]
            outlog=outfold+stresn+a+'.log'
            try: di=RJB_lib.extract_CC_R_Rfree_from_polder_log (outlog)
            except:
                print ('Failure extracting CC, R, Rfree of file:',outlog)
                di = False
                quitt = True
                #exit()
            if di!=False:
                lisel.append( (a,di['cc13'],di) )
                if dic_pdb[ch[0]][resn]==a: dicallSC[ch][resn].append( [a+'!',di['cc13']] )
                else:                       dicallSC[ch][resn].append( [a,di['cc13']] )

            #d={'cc13':cc13,'rmodel':rmodel,'rfreemodel':rfreemodel,'rexcl':rexcl,'rfreeexcl':rfreeexcl,'rimpr':rimpr,'rfreeimpr':rfreeimpr}
            #'\t'+di['cc13']+'\t'+di[]+'\t'++'\t'+)
        lisel=sorted(lisel,key=(lambda item: item[1]), reverse=True)
        #for i in lisel[0]:
        for i in lisel:
            di=i[2]
            for key,value in di.items():
                di[key]='%.3f'%(value)
            tab=['cc13','rmodel','rfreemodel','rimpr','rfreeimpr','rsimpr']
            ou2.write('\n'+i[0])
            if dic_pdb[ch[0]][resn]==i[0]:
                ou2.write('!')
            for t in tab:
                ou2.write('\t'+di[t])
        ou2.close()

print ('Building Summary of polder RSCC (cc13) of main chain atoms')
#### Building dictionary of polder RSCC (cc13) of main chain atoms
dicmainchain={}
for ch, dires in dic_pos_aa.items():
    dicmainchain[ch] = {}
    n = 0
    for resn, lmut in dires.items():
        n += 1
        stresn = str(resn)
        outfmc = output_folder + '/mainchain/' + ch + stresn + '_polder.log'
        try:
            di = RJB_lib.extract_CC_R_Rfree_from_polder_log(outfmc)
        except:
            print('Failure extracting CC, R, Rfree of file:', outfmc)
            di = False
            quitt=True
        if di==False: di={'cc13':'Null'}
        else:         di['cc13'] = '%.1f' % (di['cc13']*100)
        dicmainchain[ch][resn] = di['cc13']

if quitt: exit()

####writting overall table

dic_impartialres=RJB_lib.return_impartial_res (pdb)

ou=open(output_folder+'_summary.log','w')
ou2=open(output_folder+'_all.log','w')
ou3=open(output_folder+'_all2.log','w')
tab=['Chain','ResN','AAcid','RSCC_sc','DContr','SideCh']
tab3=['Chain','ResN','AAcid','RSCC_sc','DContr','R','Rfree','Rimp','RfImp','SideCh']
if alig:
    tab.append('Align%')
    tab.append('Impartial?')

ou.write(tab[0])
ou2.write(tab[0])
ou3.write(tab3[0])
for t in tab[1:]:
    ou.write('\t'+t)
    ou2.write('\t'+t)
for t in tab3[1:]:
    ou3.write('\t'+t)

dall={'polderCC':{},'CCSa':{},'ZCCa':{},'ZDa':{}}
eval=['CorDist','Within1','Within2','Within3','Within4','Within5']
for e in eval:
    for k,v in dall.items():
        v[e]=0

##for ch,dires in dic_pos_aa.items():
##    n=0
##    for resn,lmut in dires.items():
##        print ch,resn,lmut
##exit()

for ch,dires in dic_pos_aa.items():
    n=0
    for resn,lmut in dires.items():
        #print ch,resn
        n+=1
        stresn=str(resn)
        outfold=output_folder+'/'+ch+'/'+stresn+'/'
        #summary on POLDER STATS
        listadic=RJB_lib.extract_table_to_list_of_dict_with_first_line_as_key ( outfold+ch+'_'+stresn+'_polder.log' )
        listadic=RJB_lib.given_list_dic_all_from_table_return_same_list_with_floats (listadic)
        listadic=sorted(listadic,key=(lambda item: item['CC1,3']), reverse=True)
        #print 'len(listadic)',len(listadic)

        #writting full polder summary
        for fulli,fulll in enumerate(listadic):
                ou2.write('\n'+ch+'\t'+stresn)
                ou2.write('\t'+listadic[fulli]['Residue'])
                ffcc='%.1f'%(listadic[fulli]['CC1,3']*100)
                ou2.write('\t'+ffcc)
                if fulli==len(listadic)-1:
                    ffsdif='last'
                else:
                    ffdif=100*(listadic[fulli]['CC1,3']-listadic[fulli+1]['CC1,3'])
                    ffsdif='%.1f'%( ffdif)
                ou2.write('\t'+ffsdif)
                if alig:
                    if listadic[fulli]['Residue'][0] in dic_ali[resn]:
                        ou2.write('\t'+str(dic_ali[resn][listadic[fulli]['Residue'][0]]))
                    else:
                        ou2.write('\t0')
                if ch in dic_impartialres and resn in dic_impartialres[ch]:
                    ou2.write('\tYes')
                else:
                    ou2.write('\tNo')

        #writting full polder summary 3
        FR=False
        bcc = listadic[0]['CC1,3'] * 100
        for fulli,fulll in enumerate(listadic):
            if '!' in listadic[fulli]['Residue']: FR=True
            if not FR or bcc-listadic[fulli]['CC1,3']*100<3 or '!' in listadic[fulli]['Residue']:
                #if stresn == '64': print ('here')
                ou3.write('\n'+ch+'\t'+stresn)
                ou3.write('\t'+listadic[fulli]['Residue'])
                ffcc='%.1f'%(listadic[fulli]['CC1,3']*100)
                ou3.write('\t'+ffcc)
                if fulli==len(listadic)-1:
                    ffsdif='last'
                else:
                    ffdif=bcc-100*(listadic[fulli+1]['CC1,3'])
                    ffsdif='%.1f'%( ffdif)
                ou3.write('\t'+ffsdif)
                ou3.write('\t'+'%.1f'%(listadic[fulli]['R']*100)+'\t'+'%.1f'%(listadic[fulli]['Rfree']*100)+
                          '\t'+'%.1f'%(listadic[fulli]['Rimp'])+'\t'+'%.1f'%(listadic[fulli]['RfImp']))
                if ch in dic_impartialres and resn in dic_impartialres[ch]: ou3.write('\tYes')
                else:                                                       ou3.write('\tNo')
        ou3.write('\n')
        ou3.write(ch + '\t' + stresn)
        ou3.write('\t\t' + dicmainchain[ch][resn] + '\tBaselineMainChain\n')

        for numb in range(5):
            #print 'numb',numb
            if numb+1<len(listadic):
                ou.write('\n'+ch+'\t'+stresn)
                ou.write('\t'+listadic[numb]['Residue'])
                cc='%.1f'%(listadic[numb]['CC1,3']*100)
                ou.write('\t'+cc)
                dif=100*(listadic[numb]['CC1,3']-listadic[numb+1]['CC1,3'])
                sdif='%.1f'%( dif)
                ou.write('\t'+sdif)
                if alig:
                    if listadic[numb]['Residue'][0] in dic_ali[resn]:
                        ou.write('\t'+str(dic_ali[resn][listadic[numb]['Residue'][0]]))
                    else:
                        ou.write('\t0')
                if ch in dic_impartialres and resn in dic_impartialres[ch]:
                    ou.write('\tYes')
                else:
                    ou.write('\tNo')
            else:
                ou.write('\t\t\t')
        ou.write('\n')
        #to stay in same line! A B & AB
        #ou2.write('\n'*(20-fulli))
        if len(listadic)>0: ou2.write('\n')
        
        #correctness statistics
        check=True
        #for ii in listadic:
        #    print ii
        for i in range(5):
            #POLDER
            #print 'i',i
            #print len(listadic)
            #print len(listadic[i]['Residue'])
            #print listadic[i]['Residue']
            if len(listadic)>=i+1:
                if len(listadic[i]['Residue'])>1 and listadic[i]['Residue'][1]=='!':
                    if len(listadic)>i+1:
                        dif=100*(listadic[0]['CC1,3']-listadic[1]['CC1,3'])
                    elif len(listadic)==i+1:
                        dif=100
                    #print 'resn',resn,'correct','in position',i
                    if i==0 and dif>=1.5:
                        dall['polderCC']['CorDist']+=1
                    for ii in range(i,5):
                        #print 'within',ii
                        dall['polderCC']['Within'+str(ii+1)]+=1
                    check=False
        if check:
            print ('true answer not found in POLDER evaluation',resn)

ou.close()
ou2.close()
ou3.close()

#print dall

ou=open(output_folder+'_correctness.log','w')
lall=['polderCC']
eval=['CorDist','Within1','Within2','Within3','Within4','Within5']

for stati in lall:
    ou.write('Summary of correctness of statistic in number of residues: '+stati)
    for e in eval:
        if e=='CorDist':
            ou.write('\n#Res distinguished:\t')
        else:
            ou.write('\n#Res within best '+e[-1]+'\t')
        ou.write(str(dall[stati][e]))
        #write percentage
        #print
        ou.write('\t'+'%.1f'%(int(dall[stati][e])*100/TotResN)+'%')
    ou.write('\n\n')
ou.close()

maps=open(output_folder+'_polder_coot_open_maps','w')
maps.write( '(set-default-initial-contour-level-for-difference-map 3.0)\n')
lfmaps=[]
pathorig=os.getcwd()

for ch,dires in dic_pos_aa.items():
    for resn,lmut in dires.items():
    #for ch,dires in dic_pos_aa.items():
        #print ch
        #for resn,lmut in dires.items():
        for ch,dires2 in dic_pos_aa.items():
            if resn in dires2:
                for file in os.listdir(output_folder+'/'+ch+'/'+str(resn)):
                    if file.endswith('.mtz') and 'polder' in file:
                        f=pathorig+'/'+output_folder+'/'+ch+'/'+str(resn)+'/'+file
                lfmaps.append(f)
for fm in lfmaps:
    maps.write( '( make-and-draw-map "' + str(fm)+'" "mFo-DFc_polder" "PHImFo-DFc_polder" "" 0 1)\n')
maps.close()

##
##
####FINAL TABLE

ou=open(output_folder+'_final_model.log','w')
tab=['Chain','ResN','AAcid','RSCC_mc','DContr']
if alig:
    tab.append('Align%')
    tab.append('Impartial?')
ou.write(tab[0])
for t in tab[1:]:   ou.write('\t'+t)



if 'SINGLE' not in typeee:

    print ('Generating Final Table')

    for ch, dires in dic_pos_aa.items():
        for resn, lmut in dires.items():
            restype=dic_pdb[ch][resn]
            print ('Chain',ch,'and residue',resn,'is',restype,'in model')
            #n+=1
            stresn=str(resn)
            outfold=output_folder+'/'+ch+'/'+stresn+'/'
            #outfold=output_folder+'/'+'AB'+'/'+stresn+'/'
            #summary on POLDER STATS
            listadic=RJB_lib.extract_table_to_list_of_dict_with_first_line_as_key ( outfold+ch+'_'+stresn+'_polder.log' )
            #listadic=RJB_lib.extract_table_to_list_of_dict_with_first_line_as_key ( outfold+'AB_'+stresn+'_polder.log' )
            listadic=RJB_lib.given_list_dic_all_from_table_return_same_list_with_floats (listadic)
            listadic=sorted(listadic,key=(lambda item: item['CC1,3']), reverse=True)
            #print 'len(listadic)',len(listadic)
            verif=False
            for i,l in enumerate(listadic):
                if not verif:
                    ou.write ('\n' + ch + '\t' + stresn)
                    ou.write ('\t' + l['Residue'])
                    cc='%.1f'%(l['CC1,3']*100)
                    ou.write('\t'+cc)
                    if i+1!=len(listadic):
                        dif=100*(l['CC1,3']-listadic[i+1]['CC1,3'])
                        sdif='%.1f'%( dif)
                    else:
                        sdif='last'
                    ou.write('\t'+sdif)
                    if alig:
                        if l['Residue'][0] in dic_ali[resn]:
                            ou.write('\t'+str(dic_ali[resn][l['Residue'][0]]))
                        else:
                            ou.write('\t0')
                    if ch in dic_impartialres and resn in dic_impartialres[ch]:
                        ou.write('\tYes')
                    else:
                        ou.write('\tNo')
                    # if l['Residue'].startswith(restype):
                    #     verif=True

            ####write base line of main chain polder RSCC
            ou.write('\n'+ch+'\t'+stresn)
            ou.write('\t\t' + dicmainchain[ch][resn] + '\tBaselineMainChain')

            ou.write('\n')

    ou.close()

# if 'ALIGN' in typeee:
#     for i in dic_ali:
#         print (i,dic_ali[i])



now = datetime.datetime.now()
date = now.isoformat()[:10] + ' ' + now.isoformat()[11:16]
print ('Finishing '+sys.argv[0]+' '+date)




####writting table main chain

ou = open(output_folder + '/mainchain/mainchain_summary.log', 'w')
tab = ['Chain', 'ResN', 'RSCC_mc']
ou.write(tab[0])
for t in tab[1:]:
    ou.write('\t' + t)

for ch,dicresn in dicmainchain.items():
    for resn, cc13 in dicresn.items():
        stresn=str(resn)
        ou.write('\n'+ch+'\t'+stresn+'\t')
        ou.write(cc13)
ou.close()

tab = ['Chain', 'ResN', 'AAcid','DContr','RSCC_mc']
outResolved=open(output_folder+'_Resolved.log','w')
outDubious =open(output_folder+'_Dubious.log','w')

for i in tab:
    outResolved.write(i+'\t')


for ch in dicallSC:
    for resn in dicallSC[ch]:
        #print dicallSC[ch][resn]
        lvar=sorted(dicallSC[ch][resn], key=(lambda item: item[1]), reverse=True)
        #print (ch,resn,lvar)
        if len(lvar)>1 and '!' in lvar[0][0] and lvar[0][1]-lvar[1][1]>0.03:
            cc13='%.1f'%(lvar[0][1]*100)
            ccdiff= '%.1f'%((lvar[0][1]-lvar[1][1])*100)
            outResolved.write('\n'+ch+'\t'+str(resn)+'\t'+lvar[0][0]+'\t'+cc13 + '\t' + ccdiff + '\t' + dicmainchain[ch][resn])
        else:
            for i,l in enumerate(lvar):
                cc13 = '%.1f' % (l[1] * 100)
                #print i,l,len(lvar)
                if i+1!=len(lvar):
                    ccdiff = '%.1f' % ((lvar[i][1] - lvar[i+1][1]) * 100)
                else: ccdiff = 'last'
                outDubious.write('\n'+ch+'\t'+str(resn)+'\t'+l[0]+'\t'+cc13 + '\t' + ccdiff )
            outDubious.write('\n'+ch+'\t'+str(resn)+'\tMainCh\t'+dicmainchain[ch][resn]+'\n')


outResolved.close()
outDubious.close()

#Writting Analysis based on PDB/MTZ data statistics

#llabels=['PDB','Source','#Res','#ProtCh','Res(A)','Compl%','ComplHigh%','Redundancy','SG','B-Wilson','B-Wilson-mean','I/Sigma','I/Sigma(High)','Specie','R','Rfree']
llabels=['PDB','Res','RPDB','RfreePDB','ComplPDB','ComplAllX','ComplHRPDB','ComplHRX','ISX','WPX']
amino_acid_list_3L=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR','MSE']

out=output_folder+'_DATA_statistics.log'

fw = open(out, 'w')
fw.write(llabels[0])
for i in llabels[1:]: fw.write(',' + i)

#Getting values from PDB
with open(pdb) as f: frPDB = f.read()

try:        iR=frPDB.index('REMARK   3   R VALUE            (WORKING SET)')+48
except:     iR=frPDB.index('REMARK   3   R VALUE          (WORKING SET, NO CUTOFF) :')+57
RPDB=re.sub(r'[^0-9.]', '',frPDB[iR:iR+7])
try:    iRfree=frPDB.index('REMARK   3   FREE R VALUE                    ')+48
except: iRfree=frPDB.index('REMARK   3   FREE R VALUE                  (NO CUTOFF) :')+57
RfreePDB=re.sub(r'[^0-9.]', '',frPDB[iRfree:iRfree+7])

if 'REMARK   2 RESOLUTION.' in frPDB: iRes=frPDB.index('REMARK   2 RESOLUTION.')+23
elif 'REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :' in frPDB: iRes=frPDB.index('REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :')+48
elif 'REMARK   3   BIN RESOLUTION RANGE HIGH           :' in frPDB: iRes=frPDB.index('REMARK   3   BIN RESOLUTION RANGE HIGH           :')+51
Res=re.sub(r'[^0-9.]', '',frPDB[iRes:iRes+8])


if   'REMARK 200  COMPLETENESS FOR RANGE     (%) :'      in frPDB: iCompl = frPDB.index('REMARK 200  COMPLETENESS FOR RANGE     (%) :') + 44
elif 'REMARK   3   COMPLETENESS FOR RANGE        (%) :' in frPDB: iCompl = frPDB.index('REMARK   3   COMPLETENESS FOR RANGE        (%) :')+49
ComplPDB=re.sub(r'[^0-9.]', '',frPDB[iCompl:iCompl+6])
if 'REMARK 200  COMPLETENESS FOR SHELL     (%) :' in frPDB: iComplHR = frPDB.index('REMARK 200  COMPLETENESS FOR SHELL     (%) :') + 44
elif 'REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) :' in frPDB: iComplHR = frPDB.index('REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) :') + 53
try:
    ComplHRPDB = re.sub(r'[^0-9.]', '',frPDB[iComplHR:iComplHR + 6])
except:
    iPhenixHR=frPDB.index('REMARK   3   BIN  RESOLUTION RANGE  COMPL.    NWORK NFREE   RWORK  RFREE  CCWORK CCFREE')
    fPhenixHR=frPDB[iPhenixHR:].index('REMARK   3  \n')+iPhenixHR
    lPhenixHR=frPDB[iPhenixHR:fPhenixHR].split('\n')
    #print (lPhenixHR)
    for l in lPhenixHR:
        if l.startswith('REMARK   3    '):
            #print(l)
            ComplHRPDB=float(l[38:42])*100
    ComplHRPDB='%.1f'%(ComplHRPDB)
    #print (ComplHRPDB)
    #exit()


# try:
#  iRed=frPDB.index('REMARK 200  DATA REDUNDANCY                :')+44
#  RedPDB=frPDB[iRed:iRed+6].replace(' ', '')
# except: RedPDB='NULL'
#
# try:
#  iWP=frPDB.index('REMARK   3   FROM WILSON PLOT           (A**2) :')+48
#  WPPDB=frPDB[iWP:iWP+6].replace(' ', '')
# except: WPPDB='NULL'
#
# try:
#  iISall=frPDB.index('REMARK 200  <I/SIGMA(I)> FOR THE DATA SET  :')+44
#  ISallPDB=frPDB[iISall:iISall+9].replace(' ', '')
#  iISHR=frPDB.index('REMARK 200  <I/SIGMA(I)> FOR SHELL         :')+44
#  ISHRPDB=frPDB[iISHR:iISHR+9].replace(' ', '')
# except:
#  ISallPDB='NULL'
#  ISHRPDB = 'NULL'

#Running phenix.xtriage
if not os.path.isfile(mtz_init[:-4]+'_xtriage.log'):
 print('Running phenix.xtriage '+mtz_init+' > '+mtz_init[:-4]+'_xtriage.log')
 os.system('phenix.xtriage '+mtz_init+' > '+mtz_init[:-4]+'_xtriage.log')

#Getting values from phenix.xtriage
if os.path.isfile(mtz_init[:-4]+'_xtriage.log'):
 ISX,ComplAllX,ComplHRX,ComplIS1,ComplIS3,ComplIS5,ComplIS10,ComplIS15,WPX,t1i,t2f,t3e,t4l,t5l2,t6z='NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL','NULL'

with open(mtz_init[:-4]+'_xtriage.log') as f: fX=f.readlines()
checkCompl=False
checkComplBinning=False
for i,l in enumerate(fX):

    if l.startswith('Completeness in resolution range: '): ComplAllX=l.split()[-1]
    if l.startswith('Overall <I/sigma> for this dataset is'): ISX=l.split()[-1]
    if l.startswith('  | Completeness and data strength'): checkCompl=True
    if checkCompl and l.startswith('  The completeness of data for which I/sig(I)>3.00'):
     checkCompl=False
     lHR=fX[i-3].split()
     #print (lHR)
     ComplIS1,ComplIS3,ComplIS5,ComplIS10,ComplIS15=lHR[5],lHR[9],lHR[11],lHR[13],lHR[15]
    if l.startswith('                 ----------Completeness (log-binning)----------'): checkComplBinning=True
    #if checkComplBinning: print (l[:-1])
    if checkComplBinning and l.startswith('               ----------Analysis of resolution limits------'):
     checkComplBinning=False
     #print ('Here')
     ComplHRX=fX[i-3].split()[-2][:-1]
    if l.endswith(' A**2\n'): WPX=l.split()[0]

    if l.startswith('  <I^2>/<I>^2 :'): t1i=l[15:21].replace(' ','')
    if l.startswith('  <F>^2/<F^2> :'): t2f=l[15:21].replace(' ', '')
    if l.startswith('  <|E^2-1|>   :'): t3e=l[15:21].replace(' ', '')
    if l.startswith('  <|L|>       :'): t4l=l[15:21].replace(' ', '')
    if l.startswith('  <L^2>       :'): t5l2=l[15:21].replace(' ', '')
    if l.startswith('  Multivariate Z score L-test:'): t6z = l.split()[-1]


ComplAllX=float(ComplAllX)*100
ComplAllX = '%.1f' % (ComplAllX)

# if ISallPDB!='NULL':
#  ISallPDB = float(ISallPDB)
#  ISallPDB = '%.1f' % (ISallPDB)

fw.write('\n'+pdb[:-4]+','+Res+','+RPDB+','+RfreePDB+','+ComplPDB+','+ComplAllX+','+ComplHRPDB+','+ComplHRX+','+ISX+','+WPX)
