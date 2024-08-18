import os
import glob
import numpy as np
from pathlib import Path
from sklearn import neighbors
from Bio.PDB import PDBParser
from utils.parser import read_wrl2,parsefile

#hfile = open('hydro.csv', 'r')
#hdic = {}
#for line in hfile:
#    ll = line.split(':')
#    if len(ll) > 1:
#        hdic[ll[0].upper()] = float(ll[1])
# print hdic
# print len(hdic)
hdic={"ALA":  1.800,
      "ARG": -4.500,
      "ASN": -3.500,
      "ASP": -3.500,
      "CYS": 2.500,
      "GLN": -3.500,
      "GLU": -3.500,
      "GLY": -0.400,
      "HIS": -3.200,
      "ILE": 4.500,
      "LEU": 3.800,
      "LYS": -3.900,
      "MET": 1.900,
      "PHE": 2.800,
      "PRO": -1.600,
      "SER": -0.800,
      "THR": -0.700,
      "TRP": -0.900,
      "TYR": -1.300,
      "VAL": 4.200}
edic = {}
charge = ['LYS', 'ARG', 'HIS']
necharge = ['ASP', 'GLU']
for aa in hdic:
    if aa in charge:
        edic[aa] = 1
    elif aa in necharge:
        edic[aa] = -1
    else:
        edic[aa] = 0

def get_pqr_from_pdb(pdb_path,output_dir,save_name=None):
    pdb2pqr = 'pdb2pqr'  # '/path/to/pdb2pqr-linux-bin64-2.1.0/pdb2pqr'
    apbsflag = '--whitespace --ff=AMBER --apbs-input'
    apbs = '/projectnb2/docking/imhaoyu/.conda/pkgs/apbs-1.5-hd8ae8d1_0/bin/apbs'
    assert os.path.exists(pdb_path), sys.stderr
    assert pdb_path.endswith(".pdb"), f"Not a pdb file!"
    # necessary!
    # https://groups.google.com/g/apbs-users/c/IFNWd_G42gc
    input_dir=os.path.dirname(pdb_path)
    os.chdir(input_dir)
    if save_name:
        pdb_basename=save_name
    else:
        pdb_basename = Path(pdb_path).stem
    command_format=f"pdb2pqr {apbsflag} {os.path.join(output_dir, pdb_basename + '.in')} {pdb_path} {os.path.join(output_dir, pdb_basename + '.pqr')}"
    try:
        os.system(command_format)
    except:
        print('pqr:' + pdb_path)
    command_format = f"apbs {pdb_basename}.in"
    try:
        os.chdir(input_dir)
        os.system(command_format)
        print("Complete!")
    except:
        print('apbs: ' + f)

def findvalue(dotcloud, gridsize, origin, delta, value):
    ind3d = np.floor((dotcloud - origin) / delta)
    ind1d = ind3d[:, 2] + ind3d[:, 1] * gridsize[2] + ind3d[:, 0] * gridsize[1] * gridsize[2]
    # print value.shape
    #     print ind3d
    # print value[ind1d.astype(int)]
    temp = value.reshape(gridsize)
    #     print temp[ind3d[:,0].astype(int),ind3d[:,1].astype(int),ind3d[:,2].astype(int)]
    return value[ind1d.astype(int)]


def getcontactbyabag(pdb_path):
    parser = PDBParser()
    structure = parser.get_structure('C', pdb_path)
    allchain=[]
    for chain in structure.get_chains():
        allchain.append(chain.get_id())
    newdick = []
    labeldick = [[],[]]
    for chain in allchain:
        try:
            test = structure[0][chain]
        except Exception as e:
            print(e)
            continue
        for resi in structure[0][chain]:
            if resi.get_resname() not in hdic.keys():
                continue
            cen = [0, 0, 0]
            count = 0
            for atom in resi:
                    # print atom.get_coord()
                    # print list(atom.get_vector())
                cen[0] += atom.get_coord()[0]
                cen[1] += atom.get_coord()[1]
                cen[2] += atom.get_coord()[2]
                count += 1
            cen = [coor * 1.0 / count for coor in cen]
            newdick.append(cen)
            labeldick[0].append(hdic[resi.get_resname()])
            labeldick[1].append(edic[resi.get_resname()])
    return newdick, labeldick


def add_EH_to_pts(pdb_path,pts_path,apbs_path,save_dir,save_name):
    newdict, labeldict = getcontactbyabag(pdb_path)
    file=np.loadtxt(pts_path)[:,0:3]
    coord=np.asarray(newdict)
    label=np.transpose(np.asarray(labeldict))
    clf = neighbors.KNeighborsClassifier(3)
    clf.fit(coord, label * 10)
    dist, indices=clf.kneighbors(file)
    apbs=open(apbs_path,"r")
    gl, orl, dl, vl = parsefile(apbs)
    av = findvalue(file,gl,orl, dl,vl)
    pred=np.sum(label[indices] * np.expand_dims(dist, 2), 1) / np.expand_dims(np.sum(dist, 1), 1) / 10
    np.savetxt(os.path.join(save_dir, f'{save_name}.pts'),
               np.concatenate((file, np.expand_dims(av, 1), np.expand_dims(pred[:, 0], 1)), axis=1))