"""
Convert pkl PyCLES files to hdf5
    input:  case_name/Output.case_hashtag/Visualization 
    output: case_name/slices
"""

import os
import hickle as hkl
import pickle as pkl
import numpy as np

fdir = "/Users/ajaruga/clones/UWLCM/src/cases/input_data/dycoms/Output.DYCOMS_RF01.04b3a/Visualization/"

it = 1
for fn in os.listdir(fdir):
    
    # read from pkl
    with open(fdir+fn, 'rb') as f:
        data = pkl.load(f)

    # dump data to hdf
    hkl.dump(data, fdir+'../../slices/'+str(it)+'.hdf', mode='w')
    it+=1
