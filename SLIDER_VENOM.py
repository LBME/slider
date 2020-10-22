#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import range
from builtins import bytes, str

import os,sys,time
import multiprocessing
from collections import defaultdict
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/git-tools/tools")
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/git-tools/tools/Brasil")
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/SLIDER/seq_slider")
#sys.path.insert(0, "/cri4/rafael/Git_Scripts/REPO_EXTERNAL/ARCIMBOLDO/borges-arcimboldo/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/git-tools/tools")
#sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/git-tools/tools/Brasil")
#sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/REPO_EXTERNAL/ARCIMBOLDO/borges-arcimboldo/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/REPO_EXTERNAL/borges-arcimboldo/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/home/rborges/Dropbox/Git_Scripts/SLIDER/seq_slider")
#sys.path.insert(0, "/home/rjborges/Dropbox/Git_Scripts/git-tools/tools")
#sys.path.insert(0, "/home/rjborges/Dropbox/Git_Scripts/git-tools/tools/Brasil")
#sys.path.insert(0, "/home/rjborges/Dropbox/Git_Scripts/REPO_EXTERNAL/ARCIMBOLDO/borges-arcimboldo/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/home/rjborges/Dropbox/Git_Scripts/REPO_EXTERNAL/borges-arcimboldo/ARCIMBOLDO_FULL")
#sys.path.insert(0, "/home/rjborges/Dropbox/Git_Scripts/SLIDER/seq_slider")
import RJB_lib
import numpy
import datetime

pdb = sys.argv[1]
mtz_init = sys.argv[2]
mtz_phases = sys.argv[3]
output_folder= sys.argv[4]
typeee = sys.argv[5]

#typeee may be:
#MASPEC: use partial MASS SPECTROMETRY DATA
#TRYALL: try all possibilities
#CONSTRUCT: add word CONSTRUCT to construct the model while it runs... Desactivated
#SINGLE: for structures composed of single protein
#ALIGN: alignment file
#MAINCHAIN: keep all atoms in residue for generating phenix.polder CC
#now default: SIDECHAIN: skip C, N, O in generating phenix.polder CC
#SELCH: do not calculate chains given in sys.argv[6] , should be A,B,C or A
#SKIPTEST: will not use RAM memory calculation / TEST

if 'SELCH' in typeee: RemoveChains=sys.argv[6].split(',')
else:                 RemoveChains=False



print ('\n\n\n')

now = datetime.datetime.now()
date = now.isoformat()[:10] + ' ' + now.isoformat()[11:16]
print ('Initiating '+sys.argv[0]+' '+date)


RJB_lib.output_runline(output_file=output_folder+'_runline.log',printtt=True)

now = datetime.datetime.now()
date = now.isoformat()[:10] + ' ' + now.isoformat()[11:16]
print ('Running '+sys.argv[0]+' '+date)

nproc=RJB_lib.number_of_processor()
#nproc=4
#nproc=20
#nproc=24
#nproc=48
#nproc=6

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

if 'MASPEC' in typeee:
    seq=sys.argv[6]
    llseq=open(seq)
    lseq=llseq.readlines()
#for i in dic_res['A']:
for ch,lres in dic_res.items():
    for i in lres:
        if 'MASPEC' in typeee:
            check=''
            for line in lseq:
                if line[i-1]!=' ' and line[i-1]!='-':
                    dic_pos_aa[ch][i].append(line[i-1])
                    check+=line[i-1]
            if check=='':
                for a in amino_acid_list:
                    dic_pos_aa[ch][i].append(a)
        elif 'TRYALL' in typeee:
            for a in amino_acid_list:
                dic_pos_aa[ch][i].append(a)




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

if 'MASPEC' not in typeee and 'TRYALL' not in typeee and 'ALIGN' not in typeee:
    print ('Failure. Wrong option for typeee.')
    exit()


if 'SINGLE' in typeee:
    newkey=''
    for ch in dic_pos_aa:
        newkey+=ch
    newd=dic_pos_aa[ch]
    dic_pos_aa={newkey:newd}

if 'ALIGN' in typeee and not 'MASPEC' in typeee and not 'TRYALL' in typeee:
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
#if 'MASPEC' in typeee:
dic_seq_initial=RJB_lib.return_dic_sequence(pdb)
d_imp={}
impr=False
for ch,dires in dic_pos_aa.items():
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
if not os.path.isfile (pdb[:-4]+'-areaimol.pdb'): RJB_lib.runAREAIMOLccp4 (pdbfile=pdb,outpdb=pdb[:-4]+'-areaimol.pdb')


#added 20181031 to remove partial occupancy atoms from dictionary
if dic_disorder!=False:
    for ch in dic_disorder:
        for resn in dic_disorder[ch]:
            del dic_pos_aa[ch][resn]




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





for ch,dires in dic_pos_aa.items():
    #print ch
    RJB_lib.mkdir(output_folder+'/'+ch)
    for resn,lmut in dires.items():
        print ('\nEvaluating chain',ch,'and residue',resn,' with: ',', '.join(lmut),' (total: ',str(len(lmut))+')')
        print ('Running coot mutate, rotamer search and sphere refine')
        #print lmut
        #print resn
        stresn=str(resn)
        outfold=output_folder+'/'+ch+'/'+stresn+'/'
        RJB_lib.mkdir(outfold)
        if (not os.path.isfile(outfold+ch+'_'+stresn+'_polder.log') ):
            ####RUNNING COOT MODELING
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
                                process = multiprocessing.Process(target= RJB_lib.coot_mutate_sph_ref_correcting_files , args= ( pdb , mtz_phases , outf+'.pdb' , dic ,dic_pdb, 5, True , False ) )
                                                                                 #coot_mutate_sph_ref_correcting_files          (pdb,mtz_phases,outpdb,dic,dic_pdb radius_sph_ref=5, coot_path=False , printtt=False)
                                process.start()
                                time.sleep(0.1)
                                break
                    else:
                        print ("FATAL ERROR: I cannot load correctly information of CPUs.")
                        exit()


for ch, dires in dic_pos_aa.items():
    # print ch
    RJB_lib.mkdir(output_folder + '/' + ch)
    for resn, lmut in dires.items():
        stresn=str(resn)
        outfmc=output_folder + '/mainchain/' + ch + stresn + '.pdb'
        if not os.path.isfile(outfmc):
            dic = {}
            for c1 in ch:
                for c in c1:
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

###SUMMARY RESULTS PHENIX.POLDER BY CHAIN AND BY RESIDUE NUMBER
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
            di=RJB_lib.extract_CC_R_Rfree_from_polder_log (outlog)
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


#### Building dictionary of polder RSCC (cc13) of main chain atoms
dicmainchain={}
for ch, dires in dic_pos_aa.items():
    dicmainchain[ch] = {}
    n = 0
    for resn, lmut in dires.items():
        n += 1
        stresn = str(resn)
        outfmc = output_folder + '/mainchain/' + ch + stresn + '_polder.log'
        di = RJB_lib.extract_CC_R_Rfree_from_polder_log(outfmc)
        if di==False: di={'cc13':'Null'}
        else:         di['cc13'] = '%.1f' % (di['cc13']*100)
        dicmainchain[ch][resn] = di['cc13']


####writting overall table

dic_impartialres=RJB_lib.return_impartial_res (pdb)

ou=open(output_folder+'_summary.log','w')
ou2=open(output_folder+'_all.log','w')
tab=['Chain','ResN','PCCres','PCC','PCCdif%']
if alig:
    tab.append('Align%')
    tab.append('Impartial?')

ou.write(tab[0])
ou2.write(tab[0])
for t in tab[1:]:
    ou.write('\t'+t)
    ou2.write('\t'+t)

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
        ou.write('\t'+str(int(dall[stati][e])*100/TotResN)+'%')
    ou.write('\n\n')
ou.close()

maps=open(output_folder+'_polder_coot_open_maps','w')
maps.write( '(set-default-initial-contour-level-for-difference-map 3.0)\n')
lfmaps=[]
pathorig=os.getcwd()
dires=dic_pos_aa[ch]
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
    maps.write( '( make-and-draw-map "' + fm+'" "mFo-DFc_polder" "PHImFo-DFc_polder" "" 0 1)\n')
maps.close()

##
##
####FINAL TABLE

ou=open(output_folder+'_final_model.log','w')
tab=['Chain','ResN','PCCres','PCC','PCCdif%']
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

if 'ALIGN' in typeee:
    for i in dic_ali:
        print (i,dic_ali[i])



now = datetime.datetime.now()
date = now.isoformat()[:10] + ' ' + now.isoformat()[11:16]
print ('Finishing '+sys.argv[0]+' '+date)




####writting table main chain

ou = open(output_folder + '/mainchain/mainchain_summary.log', 'w')
tab = ['Chain', 'ResN', 'PCC']
ou.write(tab[0])
for t in tab[1:]:
    ou.write('\t' + t)

for ch,dicresn in dicmainchain.items():
    for resn, cc13 in dicresn.items():
        stresn=str(resn)
        ou.write('\n'+ch+'\t'+stresn+'\t')
        ou.write(cc13)
ou.close()

tab = ['Chain', 'ResN', 'PCC','PCCdiff','PCCmainCh']
outResolved=open(output_folder+'_Resolved.log','w')
outDubious =open(output_folder+'_Dubious.log','w')

for i in tab:
    outResolved.write(i+'\t')


for ch in dicallSC:
    for resn in dicallSC[ch]:
        #print dicallSC[ch][resn]
        lvar=sorted(dicallSC[ch][resn], key=(lambda item: item[1]), reverse=True)
        #print (resn,lvar)
        if '!' in lvar[0][0] and len(lvar)>1 and lvar[0][1]-lvar[1][1]>0.03:
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

