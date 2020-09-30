#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, sys, time
##import multiprocessing
##from collections import defaultdict
##import RJB_lib
import numpy as np
import math
import datetime
from Bio.PDB import *

'''
parser = PDBParser()
structure = parser.get_structure('1clp', 'C:/Users/Usuario/PycharmProjects/SLIDER/Testes/1clpTESTE.pdb')
atomresn_list = []
for model in structure:
    for chain in model:
        if chain == structure[0]['C']:
            for residue in chain:
                if residue == structure[0]['C'][1]:
                    for atom in residue:
                        if atom.get_id() == 'CA':
                            atomresn_list.append(structure[0]['C'][1]['CA'])
                        elif atom.get_id() == 'C':
                            atomresn_list.append(structure[0]['C'][1]['C'])
                        elif atom.get_id() == 'N':
                            atomresn_list.append(structure[0]['C'][1]['N'])
print(atomresn_list)
sup = Superimposer()
sup.set_atoms([structure[0]['A'][1]['CA'],structure[0]['A'][1]['N'],structure[0]['A'][1]['C']],[structure[0]['C'][1]['CA'],structure[0]['C'][1]['N'],structure[0]['C'][1]['C']])
#sup.set_atoms([structure[0]['A'][1]['CA'], structure[0]['A'][1]['N'], structure[0]['A'][1]['C']],atomresn_list)
#print sup.rotran
#print sup.rms
sup.apply([structure[0]['C'][1]['CA'],structure[0]['C'][1]['N'],structure[0]['C'][1]['C']])
oiCaArg = structure[0]['A'][1]['CA'].get_coord()
oiCaAlanina = structure[0]['C'][1]['CA'].get_coord()
oiNArg = structure[0]['A'][1]['N'].get_coord()
oiNAlanina = structure[0]['C'][1]['N'].get_coord()
print(oiCaArg,oiCaAlanina)


#def cootAutomaticRotamersSave (A3L):
'''
'''
dicionarioResidues = {'A3L': ['ARG','ASN','ASP','CYS','GLN','GLY','GLU','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']}
dicionarioResiduesRotamers = {  'ARG': ['ptp180', 'ptt85','ptt180', 'ptt-85', 'ptm180', 'ptm-85', 'tpp85', 'tpp180', 'tpt85', 'tpt180', 'ttp85', 'ttp180', 'ttp-105', 'ttt85', 'ttt180', 'ttt-85', 'ttm105', 'ttm180', 'ttm-85', 'mtp85', 'mtp180', 'mtp-105', 'mtt85', 'mtt180', 'mtt-85', 'mtm105', 'mtm180', 'mtm-85', 'mmt85', 'mmt180', 'mmt-85', 'mmm180', 'mmm-85'],
                                'ASN': ['p-10', 'p30', 't-20', 't30', 'm-20', 'm-80', 'm120'],
                                'ASP': ['p-10', 'p30', 't0', 't70', 'm-20'],
                                'CYS': ['p', 't', 'm'],
                                'GLN': ['mm100', 'mm-40','mt-30', 'mp0', 'tt0', 'tp60', 'tp-100', 'pm0', 'pt20'],
                                'GLY': [],
                                'GLU': ['pt-20', 'pm0', 'tp10', 'tt0', 'tm-20', 'mp0', 'mt-10', 'mm-40'],
                                'HIS': ['p-80', 'p80', 't-160', 't-80', 't60', 'm-70', 'm170', 'm80'],
                                'ILE': ['pp', 'pt', 'tp', 'tt', 'mp', 'mt', 'mm'],
                                'LEU': ['pp', 'tp', 'tt', 'mp', 'mt'],
                                'LYS': ['ptpt', 'pttp', 'pttt', 'pttm', 'tptp', 'tptt', 'tptm', 'ttpp', 'ttpt', 'tttp', 'tttt', 'tttm', 'ttmt', 'mtpp', 'mtpt', 'mttp', 'mttt', 'mttm', 'mtmt', 'mtmm', 'mmtp', 'mmtt', 'mmtm', 'mmmt'],
                                'MET': ['ptp', 'ptm', 'tpp', 'tpt', 'ttp', 'ttt', 'ttm', 'mtp', 'mtt', 'mtm', 'mmp', 'mmt', 'mmm'],
                                'PHE': ['p90', 't80', 'm-85', 'm-30'],
                                'PRO': ['Cg endo', 'Cg exo', 'cis Cg endo'],
                                'SER': ['p', 't', 'm'],
                                'THR': ['p', 't', 'm'],
                                'TRP': ['p-90', 'p90', 't-105', 't90', 'm-90', 'm0', 'm95'],
                                'TYR': ['p90', 't80', 'm-85', 'm-30'],
                                'VAL': ['p', 't', 'm']
                              }

value4321 = dicionarioResidues['A3L']
for v in value4321:
    value1234 = dicionarioResiduesRotamers[''+v+'']
    coot_path = 'coot'
    if not os.path.exists("/home/jbruno/PycharmProjects/SLIDER/aa/"+v+""):
        os.mkdir("/home/jbruno/PycharmProjects/SLIDER/aa/"+v+"")
        print("Directory", "/home/jbruno/PycharmProjects/SLIDER/aa/"+v+"", "Created")
    else:
        print("Directory", "/home/jbruno/PycharmProjects/SLIDER/aa/"+v+"", "already exists!")
    script_file = open('/home/jbruno/PycharmProjects/SLIDER/aa/'+v+'/'+v+'.py', 'w')
    script_file.write('mutate (0,"A",1,"","'+v+'")\n')
    for k in value1234:
        script_file.write('set_residue_to_rotamer_name (0,"A",1,"","","'+k+'")\n')
        script_file.write('write_pdb_file (0,"/home/jbruno/PycharmProjects/SLIDER/aa/'+v+'/'+v+'_'+k+'.pdb")\n')
    script_file.close()
    os.system( coot_path +' --pdb "/home/jbruno/PycharmProjects/SLIDER/aa/ALA-trial1.pdb" -s "/home/jbruno/PycharmProjects/SLIDER/aa/'+v+'/'+v+'.py"')
 '''

'''
parser = PDBParser()
structure = parser.get_structure('1clp', '/home/jbruno/Database2018/natural/1clp/1clp.pdb')
for model in structure:
     for chain in model:
        for residue in model.get_residues():
                print(chain,residue,residue['CA'])
                if residue.get_resname() == 'GLY':
                    print 'oi sou uma glicina'
                    n = residue['N'].get_vector()
                    c = residue['C'].get_vector()
                    ca = residue['CA'].get_vector()
                    n = n - ca
                    c = c - ca
                    rot = rotaxis(-math.pi*128/180.0, c)
                    cb_at_origin = n.left_multiply(rot)
                    cb = cb_at_origin + ca
                    print(cb)

'''

# now = datetime.datetime.now()
# date = now.isoformat()[:10] + ' ' + now.isoformat()[11:16]
# print 'Running '+sys.argv[0]+' '+date

# nproc=RJB_lib.number_of_processor()
# nproc=4
# nproc=20
# nproc=24
# nproc=48


# amino_acid_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
# amino_acid_list_3L=['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']


############# CODE START HERE
# Precisa de dois arquivos para abrir; o input(base) e o model(alvo), comentei a parte que ele salva a translação para não editar seus arquivos.
pdbinput = sys.argv[1]
pdbmodel = sys.argv[2]


class InfoFile:
    def __init__(self, type, atomindex, atomtype, rest, ch, resn, v1, v2, v3, occ, wbf, specificAtom, startIn):
        self.type = type
        self.atomindex = atomindex
        self.atomtype = atomtype
        self.rest = rest
        self.ch = ch
        self.resn = resn
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.occ = occ
        self.wbf = wbf
        self.specificAtom = specificAtom
        self.startIn = startIn

    def getValues(self):
        return (self.v1, self.v2, self.v3)

    def getOffset(self):
        return self.startIn

    def setValues(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

    def save(self):
        with open(pdbmodel, 'r+') as mr:
            mr.seek(self.getOffset(), 0)
            mr.write("{0:80s}".format(self.print()))
            mr.flush()
            print(self.print())

    def print(self):
        return "ATOM {0:6d} {1:4s} {2:3s} {3} {4:3d} {5:>11.3f} {6:>7.3f} {7:>7.3f} {8:>5.2f} {9:>5.2f} {10:>11s}".format(
            self.atomindex, self.atomtype, self.rest, self.ch, self.resn, self.v1, self.v2, self.v3, self.occ, self.wbf,
            self.specificAtom)


def addV(v1, v2):
    if len(v1) == len(v2):
        addvector = []
        for i in range(len(v1)):
            addvector.append(v1[i] + v2[i])
        return tuple(i for i in addvector)
    else:
        print('Vetores diferentes')
        exit()


with open(pdbinput, 'r') as fr:
    print(fr)
    frl = fr.readlines()

dictpdb = {}
for l in frl:

    if l.startswith('ATOM'):
        ch = l[21]
        atomindex = l[6:11]
        atomtype = l[12:16].replace(' ', '')
        resn = int(l[22:26])
        rest = l[17:20]
        occ = float(l[56:60])
        wbf = (float(l[61:66]) * float(l[56:60]))

        if ch not in dictpdb:
            dictpdb[ch] = {}
        if resn not in dictpdb[ch]:
            dictpdb[ch][resn] = {}
            dictpdb[ch][resn]['restype'] = rest
        dictpdb[ch][resn][atomtype] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))

CApos = dictpdb['A'][1]['CA']
Npos = dictpdb['A'][1]['N']
Cpos = dictpdb['A'][1]['C']
Opos = dictpdb['A'][1]['O']

with open(pdbmodel, 'r+') as mr:
    print(mr)
    mrl = mr.readlines()

dictpdb_model = {}
siz = 0
for m in mrl:
    if m.startswith('ATOM'):
        ch = m[21]
        atomindex = m[6:11]
        atomtype = m[12:16].replace(' ', '')
        resn = int(m[22:26])
        rest = m[17:20]
        occ = float(m[56:60])
        wbf = (float(m[61:66]) * float(m[56:60]))
        if ch not in dictpdb_model:
            dictpdb_model[ch] = {}
        if resn not in dictpdb_model[ch]:
            dictpdb_model[ch][resn] = {}
            dictpdb_model[ch][resn]['restype'] = rest
        dictpdb_model[ch][resn][atomtype] = InfoFile('ATOM', int(m[5:11]), m[12:15], rest, ch, resn, float(m[30:38]),
                                                     float(m[38:46]), float(m[46:54]), float(m[56:60]),
                                                     (float(m[61:66]) * float(m[56:60])), str(m[77:78]), siz)
        siz = siz + len(m)
        if atomtype == 'CA' or atomtype == 'C' or atomtype == 'N' or atomtype == 'O':
            model = dictpdb_model['A'][1][atomtype]
            v1 = float(m[30:38])
            v2 = float(m[38:46])
            v3 = float(m[46:54])
            # print(model.print()) ## Código para salvar a translação no modelo.
            # print(siz)
            # model.setValues(CApos[0] + v1, CApos[1] + v2, CApos[2] + v3)
            # print(CApos[0] + v1, CApos[1] + v2, CApos[2] + v3)
            # model.save()

import numpy as np

# Pos. dos átomos; a = Nalvo; b = CApos; c = Nbase. (Não fiz direto para testar, caso queira testar outros átomos como C, trocar "a" e "c".)
# a = np.array([20.646, 18.358, -36.909])
# b = np.array([20.479, 19.415, -37.955])
# c = np.array([21.413, 19.090, -39.010])

CAposC = (CApos[0], CApos[1], CApos[2])
CposC = (Cpos[0], Cpos[1], Cpos[2])
NposC = (Npos[0], Npos[1], Npos[2])
OposC = (Opos[0], Opos[1], Opos[2])

CAposModel = (CApos[0], CApos[1], CApos[2])
CposModel = (CApos[0] + 0.002, CApos[1] + 1.445, CApos[2] + -0.518)
NposModel = (CApos[0] + -0.001, CApos[1] + 0.001, CApos[2] + 1.496)
OposModel = (CApos[0] + 0.873, CApos[1] + 1.839, CApos[2] + -1.295)

P = np.array([CAposC, CposC, NposC, OposC])
Q = np.array([CAposModel, CposModel, NposModel, OposModel])
print("P", P)  # Visualização
print("Q", Q)  # Visualização


def centroid(X):
    """
    Centroid is the mean position of all the points in all of the coordinate
    directions, from a vectorset X.
    https://en.wikipedia.org/wiki/Centroid
    C = sum(X)/len(X)
    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    C : float
        centroid
    """
    C = X.mean(axis=0)
    return C


Ptrans = P - centroid(P)
Qtrans = Q - centroid(Q)
print(centroid(P))
#print(centroid(Q))
print(Ptrans)
#print((CApos[1]+Cpos[1]+Npos[1]+Opos[1])/4) #Teste de cálculo


def kabsch(P, Q):
    """
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.
    The algorithm works in three steps:
    - a centroid translation of P and Q (assumed done before this function
      call)
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm
    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    U : matrix
        Rotation matrix (D,D)
    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


U = kabsch(Ptrans, Qtrans)

teste = np.dot(Ptrans, U)
#print("U", U)
print(teste)
print(Qtrans)

'''
def rotate(a, b, c):  #A função requer 3 coordenadas atômicas, a pos. de um átomo da base, do carbono alfa e de um átomo espelhado do alvo. Ex: (Nalvo, CApos, Nbase)
    ba = a - b
    bc = c - b
    vectorA = np.cross(ba, bc)
    axisVector = vectorA / np.linalg.norm(vectorA)
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    sin_angle = (np.linalg.norm(vectorA)) / (np.linalg.norm(ba) * np.linalg.norm(bc))
#Construção da matriz de rotação específica para essas duas coordenadas com origem da esfera no CApos.
    Axx = cosine_angle + (axisVector[0] ** 2) * (1 - cosine_angle)
    Axy = axisVector[0] * axisVector[1] * (1 - cosine_angle) - axisVector[2] * sin_angle
    Axz = (axisVector[0] * axisVector[2] * (1 - cosine_angle)) + axisVector[1] * sin_angle

    Ayx = (axisVector[1] * axisVector[0] * (1 - cosine_angle)) + axisVector[2] * sin_angle
    Ayy = cosine_angle + (axisVector[1] ** 2) * (1 - cosine_angle)
    Ayz = (axisVector[1] * axisVector[2] * (1 - cosine_angle)) - axisVector[0] * sin_angle

    Azx = (axisVector[2] * axisVector[0] * (1 - cosine_angle)) - axisVector[1] * sin_angle
    Azy = (axisVector[2] * axisVector[1] * (1 - cosine_angle)) + axisVector[0] * sin_angle
    Azz = cosine_angle + (axisVector[2] ** 2) * (1 - cosine_angle)

    return np.array([[Axx, Axy, Axz], [Ayx, Ayy, Ayz], [Azx, Azy, Azz]])


Nalvo = rotate(a, b, c).dot((a - b)) + b #Produto entre a matriz de rotação e o vetor ba(Nalvo - CApos) para obter a posição novamente, comprovando que a matriz fez a rotação desejada.
print(Nalvo)
# Para que o resultado seja igual, é preciso que os vetores sejam de tamanhos iguais, por isso há uma pequena diferença, evidenciada pela pequena dif na norma.
print("Norma de (Nalvo - CApos)", np.linalg.norm(a-b))
print("Norma de (Nbase - CApos)", np.linalg.norm(c-b))
'''

'''
cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
angle = np.arccos(cosine_angle)

angleX = np.degrees(angle)
# print(angleX)
NposAntiga = a - b
NposNova = c - b
NposRelativa = ((a - b) - (c - b))
print(NposAntiga, NposNova, NposRelativa)
heading = math.atan2(NposRelativa[2], NposRelativa[0])
pitch = math.atan2(NposRelativa[1], math.sqrt(NposRelativa[0] * NposRelativa[0] + NposRelativa[2] * NposRelativa[2]))
magnitude = math.sqrt(NposRelativa[0] * NposRelativa[0] + NposRelativa[1] * NposRelativa[1] + NposRelativa[2] * NposRelativa[2])
'''

'''
def rotation_matrix(k, alfa):
    """
    Return the rotation matrix associated with counterclockwise rotation about
       the given axis by theta radians.
    """
    k = np.asarray(k)
    k = k / math.sqrt(np.dot(k, k))
    a = math.cos(alfa / 2.0)
    b, c, d = -k * math.sin(alfa / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


angulo = (-math.pi * 136.900 / 180.0)
'''
# Mrot = (rotation_matrix(CApos, angulo))


'''
def rotate(roll, pitch, heading):
    cosa = math.cos(heading)
    sina = math.sin(heading)

    cosb = math.cos(pitch)
    sinb = math.sin(pitch)

    cosc = math.cos(roll)
    sinc = math.sin(roll)

    Axx = cosa * cosb
    Axy = cosa * sinb * sinc - sina * cosc
    Axz = cosa * sinb * cosc + sina * sinc

    Ayx = sina * cosb
    Ayy = sina * sinb * sinc + cosa * cosc
    Ayz = sina * sinb * cosc - cosa * sinc

    Azx = -sinb
    Azy = cosb * sinc
    Azz = cosb * cosc

    return np.array([[Axx, Axy, Axz], [Ayx, Ayy, Ayz], [Azx, Azy, Azz]])


Mrot = (rotate(0, pitch, heading))
print(Mrot)
print(np.dot(Mrot, NposAntiga))
Px = magnitude * math.sin(heading) * math.cos(pitch)
Py = magnitude * math.sin(pitch)
Pz = magnitude * math.cos(heading) * math.cos(pitch)
print(Px, Py, Pz)
'''
'''
with open (pdbinput,'r+') as fr:
    print(fr)
    frl=fr.readlines()

dictpdb = {}
dictpdb_obj = {}
siz = 0
for l in frl:

    if l.startswith('ATOM'):
        ch = l[21]
        atomindex = l[6:11]
        atomtype = l[12:16].replace(' ','')
        resn = int(l[22:26])
        rest = l[17:20]
        occ = float(l[56:60])
        wbf = (float(l[61:66]) * float(l[56:60]))

        if ch not in dictpdb:
            dictpdb[ch] = {}
            dictpdb_obj[ch]={}
        if resn not in dictpdb[ch]:
            dictpdb[ch][resn]={}
            dictpdb_obj[ch][resn]={}
            dictpdb[ch][resn]['restype'] = rest
        dictpdb[ch][resn][atomtype]=( float(l[30:38]), float(l[38:46]), float(l[46:54]) )
        dictpdb_obj[ch][resn][atomtype] = InfoFile('ATOM', int(l[5:11]), l[12:15], rest, ch, resn, float(l[30:38]), float(l[38:46]), float(l[46:54]), float(l[56:60]), (float(l[61:66]) * float(l[56:60])), str(l[77:78]), siz)


    siz = siz + len(l)


'''


def subtractV(v1, v2):
    if len(v1) == len(v2):
        subtractvector = []
        for i in range(len(v1)):
            subtractvector.append(v1[i] - v2[i])
        return tuple(i for i in subtractvector)
    else:
        print('Vetores diferentes')
        exit()


'''
#[ATOM, 1920, CA, ALA, C, 1, 20.0f, 20.0f, 20.0f, 20.0f, 20.0f, B]



modeloCA = dictpdb['C'][1]['CA']
modeloC = dictpdb['C'][1]['C']
modeloN = dictpdb['C'][1]['N']
modeloO = dictpdb['C'][1]['O']
print(modeloCA,modeloC,modeloN,modeloO)
print(dictpdb['A'][1]['CA'])
for atomtype in dictpdb['A'][1]:
    if atomtype == 'CA':
        ResultadoCA = addV(dictpdb['A'][1]['CA'], modeloCA)
        print('CA', ResultadoCA)
    if atomtype == 'C':
        ResultadoC = addV(dictpdb['A'][1]['CA'], modeloC)
        print('C', ResultadoC)
    if atomtype == 'N':
        ResultadoN = addV(dictpdb['A'][1]['CA'], modeloN)
        print('N', ResultadoN)
        print('N', ResultadoN[0])
    if atomtype == 'O':
        ResultadoO = addV(dictpdb['A'][1]['CA'], modeloO)
        print('O', ResultadoO)

obj = dictpdb_obj['A'][1]['C']
print(obj.print())
obj.setValues(float(ResultadoC[0]), float(ResultadoC[1]), float(ResultadoC[2]))
obj.save()
'''

# for atomtype in dictpdb['A'][1]:
#  TargetResn = dictpdb['A'][1][atomtype]
# print(TargetResn)
# if


'''
#for ch in dictpdb:
    #for resn in dictpdb[ch]:
        #print('found', ch,resn)
        #print('atomCA', dictpdb[ch][resn]['CA'])
        if 'CB' not in dictpdb[ch][resn]:
            vectorN = dictpdb[ch][resn]['N']
            vectorC = dictpdb[ch][resn]['C']
            vectorCA = dictpdb[ch][resn]['CA']
            vectorNS = subtractV(vectorN, vectorCA)
            vectorCS = subtractV(vectorC, vectorCA)

            def rotation_matrix(v1, alfa):
             """
             Return the rotation matrix associated with counterclockwise rotation about
                the given axis by theta radians.
             """
             v1 = np.asarray(v1)
             v1 = v1 / math.sqrt(np.dot(v1, v1))
             a = math.cos(alfa / 2.0)
             b, c, d = -v1 * math.sin(alfa / 2.0)
             aa, bb, cc, dd = a * a, b * b, c * c, d * d
             bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
             return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                                 [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                                 [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


            angulo = (-math.pi*120.0/180.0)
            vectorCB = (np.dot(rotation_matrix(vectorCS, angulo), vectorNS))
            vectorRES = addV(vectorCB,vectorCA)
            print(vectorRES)
            ##print 'atomCB', dictpdb[ch][resn]['CB']






exit()

#'''

'''
exit()

#dfiles=defaultdict(list)

dic_res=RJB_lib.return_dic_resnumb_list (pdb)
TotResN=0
for ch,listres in dic_res.iteritems():
    TotResN+=len(listres)

dic_pdb=RJB_lib.return_dic_chain_resnumb_restype (pdb)
#delete empty keys (chains with ligs or waters)
deletekeys=[]
for ch,lres in dic_res.iteritems():
    #print ch,lres
    if len(lres)==0 or (RemoveChains!=False and ch in RemoveChains):
        deletekeys.append(ch)
for delet in deletekeys:
    del dic_res[delet]


dic_pos_aa={}
for ch in dic_res:
    dic_pos_aa[ch]=defaultdict(list)

alig=False

if 'ALIGN' in typeee:
    alig=True
    ali=sys.argv[6]
    dic_ali=RJB_lib.generate_dict_count_alignment (ali)

##cc=0
##for resn,lmut in dic_ali.iteritems():
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
for ch,lres in dic_res.iteritems():
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
'''
# if 'MASPEC' not in typeee and 'TRYALL' not in typeee and 'ALIGN' not in typeee:
#   print 'Failure. Wrong option for typeee.'
#  exit()

'''

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

'''
'''
#check if coot script mutate is going to work:
##dchkcoot={}
##lseq=[]
##for ch,d1 in dic_pdb.iteritems():
##    sseq=''
##    dchkcoot[ch]={}
##    for resn,d2 in d1.iteritems():
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
for ch,dires in dic_pos_aa.iteritems():
    for i,lmut in dires.iteritems():
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
if not os.path.isfile (pdb[:-4]+'-areaimol.pdb'): RJB_lib.runAREAIMOLccp4 (pdbfile=pdb,outpdb=pdb[:-4]+'-areaimol.pdb',symmetry=sym)


#added 20181031 to remove partial occupancy atoms from dictionary
if dic_disorder!=False:
    for ch in dic_disorder:
        for resn in dic_disorder[ch]:
            del dic_pos_aa[ch][resn]




################added 20181109 to test maximum RAM memory usage
if 'SKIPTEST' in typeee:
    print 'Test of RAM memory usage for coot and polder jobs was selected to be skipped.'
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

    print '\nTesting RAM Memory usage of coot & phenix.polder processes with chain',ch,'Residue number',resn,'Residue type',a

    counterr=0
    ifmem=RJB_lib.GetFreeMemory()
    for i in range(9):
        varmem=float( RJB_lib.GetFreeMemory() )
        if varmem<ifmem: ifmem=varmem
        time.sleep(0.1)

    print '\nFree Memory before runnning external programs',ifmem,'Mb'


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


    print '\nFree memory running 1 process of coot =',vcootfmem
    cootfmem=ifmem-vcootfmem
    print 'Maximum RAM memory spent on coot',cootfmem,'Mb'

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
    print 'Running phenix.polder to calculate CC of refine in chain(s)',ch,'and residue',stresn,'with side chain',aa,'from file:',pdbinput,'to save in:',outf+'.log\n'
    process.start()
    vpolderfmem=RJB_lib.GetFreeMemory()

    #wait phenix.polder jobs finish
    while len(multiprocessing.active_children()) != 0:
        varmem=float( RJB_lib.GetFreeMemory() )
        if varmem<vpolderfmem: vpolderfmem=varmem
        time.sleep(0.1)

    print '\nFree memory running 1 process of polder =',vpolderfmem
    polderfmem=ifmem-vpolderfmem
    print 'Maximum RAM memory spent on polder',polderfmem,'Mb'


    '\nIt was found',nproc,'processors (counting HyperThreading if available).'
    '600 Mb will be left free as tolerance and to be free for other programs.'
    NewNProcCoot=   int ( (ifmem-600)/cootfmem   )
    NewNProcPolder= int ( (ifmem-600)/polderfmem )
    print 'There is',ifmem-600,'available RAM memory'
    print 'To be divided to',cootfmem,'for coot jobs'
    print 'To be divided to',polderfmem,'for polder jobs'

    if NewNProcCoot==0 or NewNProcPolder==0:
        print 'Run requires more RAM memory than available.'
        exit()
    else:
        if NewNProcCoot   > nproc: NewNProcCoot   = int(nproc)
        if NewNProcPolder > nproc: NewNProcPolder = int(nproc)

    print '\nTherefore, it will be used:'
    print NewNProcCoot  ,'processors for coot jobs.'
    print NewNProcPolder,'processors for polder jobs.'

################added 20181109 to test maximum RAM memory usage END





for ch,dires in dic_pos_aa.iteritems():
    #print ch
    RJB_lib.mkdir(output_folder+'/'+ch)
    for resn,lmut in dires.iteritems():
        print '\nEvaluating chain',ch,'and residue',resn,' with: ',', '.join(lmut),' (total: ',str(len(lmut))+')'
        print 'Running coot mutate, rotamer search and sphere refine'
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
                        print "FATAL ERROR: I cannot load correctly information of CPUs."
                        exit()


for ch, dires in dic_pos_aa.iteritems():
    # print ch
    RJB_lib.mkdir(output_folder + '/' + ch)
    for resn, lmut in dires.iteritems():
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
for ch, dires in dic_pos_aa.iteritems():
    # print ch
    for resn, lmut in dires.iteritems():
        print '\nEvaluating chain', ch, 'and residue', resn, ' with: ', ', '.join(lmut), ' (total: ', str(len(lmut)) + ')'
        print 'Running phenix.polder to calculate CC of refine in chain(s)\n'  # print lmut
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
                        print "FATAL ERROR: I cannot load correctly information of CPUs."
                        exit()





####RUNNING PHENIX.POLDER MAIN CHAIN
for ch, dires in dic_pos_aa.iteritems():
    # print ch
    for resn, lmut in dires.iteritems():
        print 'Running phenix.polder to calculate RSCC of main chain atoms of residue',resn,'in chain(s)', ch
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
                        print "FATAL ERROR: I cannot load correctly information of CPUs."
                        exit()


while 1:
    time.sleep(0.1)
    if len(multiprocessing.active_children()) == 0:
        break

###SUMMARY RESULTS PHENIX.POLDER BY CHAIN AND BY RESIDUE NUMBER
for ch, dires in dic_pos_aa.iteritems():
    # print ch
    for resn, lmut in dires.iteritems():

            stresn=str(resn)
            outfold=output_folder+'/'+ch+'/'+stresn+'/'
            if not os.path.isfile(outfold+ch+'_'+stresn+'_polder.log'):
                ou2=open(outfold+ch+'_'+stresn+'_polder.log','w')
                ou2.write('Residue\tCC1,3\tR\tRfree\tRimp\tRfImp\tRimp')
                lisel=[]
                for a in lmut:
                    aa=amino_acid_list_3L[amino_acid_list.index(a)]
                    outlog=outfold+stresn+a+'.log'
                    print outlog
                    di=RJB_lib.extract_CC_R_Rfree_from_polder_log (outlog)
                    if di!=False:
                        lisel.append( (a,di['cc13'],di) )
                    #d={'cc13':cc13,'rmodel':rmodel,'rfreemodel':rfreemodel,'rexcl':rexcl,'rfreeexcl':rfreeexcl,'rimpr':rimpr,'rfreeimpr':rfreeimpr}
                    #'\t'+di['cc13']+'\t'+di[]+'\t'++'\t'+)
                lisel=sorted(lisel,key=(lambda item: item[1]), reverse=True)
                #for i in lisel[0]:
                for i in lisel:
                    di=i[2]
                    for key,value in di.iteritems():
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
for ch, dires in dic_pos_aa.iteritems():
    dicmainchain[ch] = {}
    n = 0
    for resn, lmut in dires.iteritems():
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
    for k,v in dall.iteritems():
        v[e]=0

##for ch,dires in dic_pos_aa.iteritems():
##    n=0
##    for resn,lmut in dires.iteritems():
##        print ch,resn,lmut
##exit()

for ch,dires in dic_pos_aa.iteritems():
    n=0
    for resn,lmut in dires.iteritems():
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
            print 'true answer not found in POLDER evaluation',resn

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
        print
        ou.write('\t'+str(int(dall[stati][e])*100/TotResN)+'%')
    ou.write('\n\n')
ou.close()

maps=open(output_folder+'_polder_coot_open_maps','w')
maps.write( '(set-default-initial-contour-level-for-difference-map 3.0)\n')
lfmaps=[]
pathorig=os.getcwd()
dires=dic_pos_aa[ch]
for resn,lmut in dires.iteritems():
#for ch,dires in dic_pos_aa.iteritems():
    #print ch
    #for resn,lmut in dires.iteritems():
    for ch,dires2 in dic_pos_aa.iteritems():
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
for t in tab[1:]:
    ou.write('\t'+t)


if 'SINGLE' not in typeee:

    print 'Generating Final Table'

    for ch, dires in dic_pos_aa.iteritems():
        for resn, lmut in dires.iteritems():
            restype=dic_pdb[ch][resn]
            print 'Chain',ch,'and residue',resn,'is',restype,'in model'
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
                    ou.write('\n'+ch+'\t'+stresn)
                    ou.write('\t'+l['Residue'])
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
                    if l['Residue'].startswith(restype):
                        verif=True

            ####write base line of main chain polder RSCC
            ou.write('\n'+ch+'\t'+stresn)
            ou.write('\t\t' + dicmainchain[ch][resn] + '\tBaselineMainChain')

            ou.write('\n')

    ou.close()

if 'ALIGN' in typeee:
    for i in dic_ali:
        print i,dic_ali[i]



now = datetime.datetime.now()
date = now.isoformat()[:10] + ' ' + now.isoformat()[11:16]
print 'Finishing '+sys.argv[0]+' '+date




####writting table main chain

ou = open(output_folder + '/mainchain/mainchain_summary.log', 'w')
tab = ['Chain', 'ResN', 'PCC']
ou.write(tab[0])
for t in tab[1:]:
    ou.write('\t' + t)

for ch,dicresn in dicmainchain.iteritems():
    for resn, cc13 in dicresn.iteritems():
        stresn=str(resn)
        ou.write('\n'+ch+'\t'+stresn+'\t')
        ou.write(cc13)
ou.close()
'''
