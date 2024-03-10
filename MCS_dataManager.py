# https://github.com/multichannelsystems/McsPyDataTools?tab=readme-ov-file
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import decimate
import glob

from McsPy import McsData
from McsPy import ureg, Q_


os.getcwd()
home = "/Users/j/Documents/LabData/DATA/32-MEA/32MEA_dataAnalysis"
batchFolder = "Python-analysis"
os.chdir(home)
if not os.path.exists(batchFolder):
    os.mkdir(home+"/"+batchFolder)

# Identify .h5 files
h5_files = glob.glob("*.h5")
# diagnosis steps
#h5_files[0]
#type(h5_files)

# Loop through .h5 files in home folder
for f in h5_files:
    print(f)
    # load .h5 file
    data = McsData.RawData(home+"/"+f)
    
    # Get raw data for Stream0.
    recording_stream0 = data.recordings[0].analog_streams[0]
    ids = [c.channel_id for c in recording_stream0.channel_infos.values()]
    channel_info = recording_stream0.channel_infos[ids[0]]
    sampling_frequency = channel_info.sampling_frequency.magnitude
    # get the timestamps for each sample
    to_idx = recording_stream0.channel_data.shape[1]
    time = recording_stream0.get_channel_sample_timestamps(0, 0, to_idx)
    # scale time to seconds:
    scale_factor_for_second = Q_(1,time[1]).to(ureg.s).magnitude
    time_in_sec = time[0] * scale_factor_for_second
    
    # create dictionary to store all the data of the array
    signal_array = {}
    signal_array["time"] =decimate(time_in_sec, q=10, zero_phase=True) 
    
    # loop through electrodes
    for i in ids:
        print(i)
        channel_id = i
    # get the signal
        signal = recording_stream0.get_channel_in_range(channel_id)
    # scale to mV
        scale_factor_for_mV = Q_(1,signal[1]).to(ureg.mV).magnitude
        signal_in_mV = signal[0] * scale_factor_for_mV
        signal_in_mV = np.array(signal_in_mV)
        #plt.plot(time_in_sec, signal_in_mV)
    # add to a dataframe
        column_name = str(i+1)
        signal_array[column_name] = decimate(signal_in_mV, q=10, zero_phase=True) 
        
      #### STILL NEED TO LOOP THROUGH MULTIPLE .H5 FILES
        
    # get the genotype and add to dictionary
    data_genotype = h5_files[0]
    data_genotype = data_genotype.split("_")[0][19:]
    signal_array["genotype"] = [data_genotype] * signal_array["time"].shape[0]
    
    #signal_array.keys()
    # merge dictionary in a data frame
    signal_array = pd.DataFrame(signal_array)
    
    # save array as .csv
    savename =  f.split(".")[0]
    savename = savename[0:len(savename)-8]+".csv"
    os.chdir(home+"/"+batchFolder)
    signal_array.to_csv( savename, index=False )
    os.chdir(home)
