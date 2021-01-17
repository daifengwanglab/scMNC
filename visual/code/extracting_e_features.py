# -*- coding: utf-8 -*-
"""
Code for AllenSDK Dataset

@author: Jiawei Huang
"""
import os
import pandas as pd
from ipfx.dataset.create import create_ephys_data_set
from ipfx.data_set_features import extract_data_set_features
from ipfx.utilities import drop_failed_sweeps

import numpy as np
from ipfx.feature_extractor import (
    SpikeFeatureExtractor, SpikeTrainFeatureExtractor
)
import ipfx.stimulus_protocol_analysis as spa
from ipfx.epochs import get_stim_epoch
from ipfx.dataset.create import create_ephys_data_set
from ipfx.utilities import drop_failed_sweeps
import matplotlib.pyplot as plt
from ipfx.stimulus_protocol_analysis import RampAnalysis
import seaborn as sns

nwb_file = os.path.join(
    "F:\\wisc\\visual-motor\\motor",
    "sub-599387254_ses-601506492_icephys.nwb"
)
data_set= create_ephys_data_set(nwb_file)
drop_failed_sweeps(data_set)
short_square_table = data_set.filtered_sweep_table(
    stimuli=data_set.ontology.long_square_names
)

#all features
cell_features, sweep_features, cell_record, sweep_records, _, _ = \
    extract_data_set_features(data_set, subthresh_min_amp=-100.0)
print(cell_record)

ID = []
title = cell_record.keys()
res = pd.DataFrame(columns=title)
path_out = "F:\\wisc\\Daifengresearch\\000020" 
files_out = os.listdir(path_out)
for file_out in files_out:
    path_in = path_out + "\\" +file_out
    files_in = os.listdir(path_in)
    for file_in in files_in:      
        try:
            data_set = create_ephys_data_set(nwb_file=path_in + "\\" +file_in)
            drop_failed_sweeps(data_set)
            cell_features, sweep_features, cell_record, sweep_records, _, _ = \
                extract_data_set_features(data_set, subthresh_min_amp=-100.0)
        except:
            cell_record = dict.fromkeys(title,np.nan)
        res = res.append([cell_record], ignore_index=True)
        ID.append(file_in)
        print(path_in + "\\" +file_in)
res["ID"] = ID
res.to_csv("F:\wisc\Daifengresearch\ef.csv",index=False,sep=',')



from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.ephys.extract_cell_features import extract_cell_features
from collections import defaultdict
from allensdk.core.nwb_data_set import NwbDataSet


def get_time_voltage_current_currindex0(nwb):
    df = nwb.sweep_table.to_dataframe()
    voltage = np.zeros((len(df['series'][0][0].data[:]), int((df.shape[0]+1)/2)))
    time = np.arange(len(df['series'][0][0].data[:]))/df['series'][0][0].rate
    current_initial = df['series'][1][0].data[12000]*df['series'][1][0].conversion
    curr_index_0 = int(-current_initial/20) # index of zero current stimulation
    current = np.linspace(current_initial, (int((df.shape[0]+1)/2)-1)*20+current_initial, \
                         int((df.shape[0]+1)/2))
    for i in range(curr_index_0):   # Find all voltage traces from minimum to 0 current stimulation
        voltage[:, i] = df['series'][0::2][i*2][0].data[0:12500]
    for i in range(10, 30):   # Find all voltage traces from 0 to highest current stimulation
        voltage[:, i] = df['series'][0::2][i*2][0].data[0:12500]
        # Find voltage trace for 0 current stimulation
    return time, voltage, current, curr_index_0

# initialize the cache
data_set = NwbDataSet(nwb_file)
sweep_numbers = data_set.get_sweep_numbers()
sweep_number = sweep_numbers[0] 
sweep_data = data_set.get_sweep(sweep_number)


nwb_file = os.path.join(
    "F:\\wisc\\visual-motor\\motor",
    "sub-mouse-AAYYT_ses-20180420-sample-2_slice-20180420-slice-2_cell-20180420-sample-2_icephys.nwb"
)

from pynwb import NWBHDF5IO
io_ = NWBHDF5IO(nwb_file, 'r+', load_namespaces=True)
nwb = io_.read()

time, voltage, current, curr_index_0 = get_time_voltage_current_currindex0(nwb)
print('time: ', time, '\nvoltage: ', voltage, '\ncurrent: ', current, '\ncurr_index_0: ', curr_index_0)


import seaborn as sns
sns.set_style()
plt.plot(time, voltage);
plt.xlabel('Time (s)', fontsize = 14)
plt.ylabel('Membrane voltage (mV)', fontsize = 14)
plt.title('Experimental traces', fontsize = 17)
sns.despine()


nwb_file = os.path.join(
    "F:\\wisc\\visual-motor\\motor",
    "example_file_path.nwb"
)

io_.write(nwb, link_data=False)
io_.close()