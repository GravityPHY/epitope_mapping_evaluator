"""
Accept a pdb file and make point cloud with feature file for it
output row format [x,y,z,E,H]
"""

import os
import glob
import numpy as np
from sklearn import neighbors

from utils.parser import read_wrl2,parsefile
from utils.conversion import pdb_to_wrl,wrl_to_pts
from utils.features import get_pqr_from_pdb,findvalue,getcontactbyabag