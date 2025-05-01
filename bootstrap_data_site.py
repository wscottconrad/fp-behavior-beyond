# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 15:45:00 2025

@author: sconrad
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from perm_test_array import perm_test_array
from bootstrap_data import bootstrap_data
from consec_idx import consec_idx

# Load combined data
tankfolder = r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\'

with open(f'{tankfolder}allDatComb.pkl', 'rb') as f:
    d = pickle.load(f)

trialData = d['trialData']  # Extract site-specific data

thres = 5  # Consecutive threshold length
pre = 5
post = 25

# Colors for plotting
gree = [0.47, 0.67, 0.19]
gris = [0.65, 0.65, 0.65]
red = [0.78, 0, 0]

for site, data in trialData.items():
   
    
    if len(data['NR']) > 0:
        signal = data['NR'][~np.isnan(data['NR']).any(axis=1)]
    

    print(f"Processing site: {site}")

    # Bootstrapping
    print('bootstrapping ...')
    btsrp_NR = bootstrap_data(signal, 5000, 0.0001)
    
    # print('bootstrapping avoid...')
    # btsrp_avoid = bootstrap_data(avoidSignal, 5000, 0.0001)

    ts = np.linspace(-pre, post, signal.shape[1])

    # Plot approach signal
    plt.figure()
    plt.plot(ts, np.mean(signal, axis=0), color=gris, label='NR')
    plt.fill_between(ts,
                     np.mean(signal, axis=0) + np.std(signal, axis=0) / np.sqrt(len(signal)),
                     np.mean(signal, axis=0) - np.std(signal, axis=0) / np.sqrt(len(signal)),
                     color=gris, alpha=0.3)

    # Plot avoid signal
    # plt.plot(ts, np.mean(avoidSignal, axis=0), color=gris, label='Avoid')
    # plt.fill_between(ts,
    #                  np.mean(avoidSignal, axis=0) + np.std(avoidSignal, axis=0) / np.sqrt(len(avoidSignal)),
    #                  np.mean(avoidSignal, axis=0) - np.std(avoidSignal, axis=0) / np.sqrt(len(avoidSignal)),
    #                  color=gris, alpha=0.3)

    # Bootstrap significance
    ymax = np.max(np.mean(signal, axis=0) + np.std(signal, axis=0) / np.sqrt(len(signal)))
    
    tmp = np.where(btsrp_NR[1, :] < 0)[0]
    if len(tmp) > 1:
        id = tmp[consec_idx(tmp, thres)]
        plt.plot(ts[id], ymax * np.ones((len(ts[id]), 2)) + 1, 's', 
                 markersize=7, markerfacecolor=gris, color=gris)
        
    tmp = np.where(btsrp_NR[0, :] > 0)[0]
    if len(tmp) > 1:
        id = tmp[consec_idx(tmp, thres)]
        plt.plot(ts[id], ymax * np.ones((len(ts[id]), 2)) + 1, 's', 
                 markersize=7, markerfacecolor=gris, color=gris)
        

    # tmp = np.where(btsrp_avoid[1, :] < 0)[0]
    # if len(tmp) > 1:
    #     id = tmp[consec_idx(tmp, thres)]
    #     plt.plot(ts[id], np.max(np.mean(avoidSignal, axis=1)) * np.ones((len(ts[id]), 2)) + 1, 's', 
    #              markersize=7, markerfacecolor=gris, color=gris)
    
    # tmp = np.where(btsrp_avoid[0, :] > 0)[0]
    # if len(tmp) > 1:
    #     id = tmp[consec_idx(tmp, thres)]
    #     plt.plot(ts[id], np.max(np.mean(avoidSignal, axis=1)) * np.ones((len(ts[id]), 2)) + 1, 's', 
    #              markersize=7, markerfacecolor=gris, color=gris)

    # plt.box(False)
    # Add vertical dashed line at x = 150
    plt.axvline(x=0, linestyle='--', color='black', linewidth=1.5)
    plt.axhline(y=0, linestyle='--', color='black', linewidth=1.5)
    
    # Hide the top and right spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    # Set tick parameters to remove right and top ticks
    plt.gca().tick_params(axis='x', which='both', direction='out', bottom=True, top=False)
    plt.gca().tick_params(axis='y', which='both', direction='out', left=True, right=False)
    plt.title(f'Site: {site} - NR trials')
    plt.ylabel('Z-Score')
    plt.xlabel('Prey laser onset (s)')
    # plt.xlabel('Movement initiation (s)')
    plt.legend()
    plt.show()
