"""
Accept a pdb file and make point cloud with feature file for it
output row format [x,y,z,E,H]
"""

import os
import sys
import glob
import numpy as np
from sklearn import neighbors
sys.path.append('../')
from utils.parser import read_wrl2,parsefile
from utils.conversion import pdb_to_wrl,wrl_to_pts
from utils.features import get_pqr_from_pdb,findvalue,getcontactbyabag,add_EH_to_pts


def make_pointclouds(pdb_path,save_dir, save_name):
    """

    :param pdb_path:
    :param save_dir:
    :param save_name:
    :return:
    """
    if not save_name:
        save_name=Path(pdb_path).stem
    pdb_to_wrl(pdb_path,save_dir,save_name)
    wrl_path=os.path.join(save_dir, f"{save_name}.wrl")
    wrl_to_pts(wrl_path,save_dir,save_name)
    pts_path=os.path.join(save_dir, f"{save_name}.pts")
    get_pqr_from_pdb(pdb_path,save_dir,save_name)
    apbs_path=os.path.join(save_dir, f"{save_name}.pqr.dx")
    add_EH_to_pts(pdb_path,pts_path,apbs_path,save_dir,save_name)


if __name__=="__main__":
    make_pointclouds("/projectnb2/docking/imhaoyu/24_epitope_mapping/epitope_mapping_evaluator/tests/point_clouds/5ZUF_l_b.pdb",
                     "/projectnb2/docking/imhaoyu/24_epitope_mapping/epitope_mapping_evaluator/tests/point_clouds", "5ZUF_test")
