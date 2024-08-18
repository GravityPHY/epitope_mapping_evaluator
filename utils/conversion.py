import os
import sys
import glob
import numpy as np
from pathlib import Path

import pymol
from pymol import cmd, stored

from sklearn import neighbors
from Bio.PDB import *
from utils.parser import read_wrl2

def pdb_to_wrl(pdb_path:str,
               output_dir:str,
               save_name=None):
    """
    convert the pdb file to wrl, and save in the output_dir
    :param pdb_path: path to a pdb file
    :param output_dir: path to an output directory, usually the directory to save your files
    :param save_name: the name of the wrl file save in output_dir, by
    default None and it means the name will be the pdb file name
    :return: None
    """
    assert os.path.exists(pdb_path), sys.stderr
    assert pdb_path.endswith(".pdb"), f"Not a pdb file!"
    if not save_name:
        save_name=Path(pdb_path).stem
    cmd.reinitialize()
    cmd.load(pdb_path)
    cmd.set('surface_quality', '0')
    cmd.show_as('surface', 'all')
    cmd.set_view('1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,300,1')
    cmd.save(os.path.join(output_dir, f"{save_name}.wrl"))
    cmd.delete('all')
    print(f"Files are saved in {os.path.abspath(output_dir)}")

def wrl_to_pts(wrl_path:str,
               output_dir:str,
               save_name=None):
    """
    convert a wrl file to pts, and save in the output_dir
    :param wrl_path: path to a wrl file
    :param output_dir: path to an output directory, usually the directory to save your files
    :return: the name of the pts file save in output_dir, by
    default None it means the name will be the pdb file name
    """
    assert os.path.exists(wrl_path), sys.stderr
    assert wrl_path.endswith(".wrl"), f"Not a wrl file!"
    if not save_name:
        save_name=Path(wrl_path).stem
    vb, _, cb, nb = read_wrl2(wrl_path)
    vecb = np.unique(vb, axis=0)
    np.savetxt(os.path.join(output_dir,f"{save_name}.pts"), np.vstack((vecb)),delimiter=' ')