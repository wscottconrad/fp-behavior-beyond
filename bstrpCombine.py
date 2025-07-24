# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 12:50:11 2025

@author: sconrad
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from perm_test_array import perm_test_array
from bootstrap_data import  bootstrap_data
from consec_idx import consec_idx

# %% Bootstrap analysis and combine plots for experiments
# Scott Conrad 20/12/2024, adapted from Isis Alonso-Lozares
# Takes data created from combineData.py

tankfolder = r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\'

with open(f'{tankfolder}allDatComb.pkl', 'rb') as f:
    d = pickle.load(f)

trialSignal = d['trialSignal']
trialSignalinit = d['ZdFoFinit']
ITIsignal = d['ITIsignal']
trialSpeed = d['trialSpeed']
speedInit = d['speedTrialsMov']
ITIspeed = d['ITIspeed']

# Define variables to save
btsrp = {}
thres = 5  # Consecutive threshold length

# Creating index for trial types: 0 = no movement, 1 = avoid prey laser, 2 = approach
trialIndex = np.array([0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 2, 0, 2, 0, 2, 1, 2, 2, 0, 2, 0, 2, 1])

# plt.figure()
# plt.plot(np.mean(trialSignal[[3, 4, 29], :], axis=0))
# plt.title('animal 1 left')

# plt.figure()
# plt.plot(np.mean(trialSignal[[10, 11, 35], :], axis=0))
# plt.title('animal 1 right')

# plt.figure()
# plt.plot(np.mean(trialSignal[[14, 18, 20, 40, 41, 43, 45, 46], :], axis=0))
# plt.title('animal 2 left')

# plt.figure()
# plt.plot(np.mean(trialSignal[[21, 25, 27, 47, 48, 50, 52, 53], :], axis=0))
# plt.title('animal 2 right')

trialIndexSingle = np.array([0, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2, 2, 0, 2, 0, 1, 0, 0, 0, 0, 2, 2, 0, 2, 0, 2, 1])

trialIndex[trialIndex == 2] = 1
trialIndexSingle[trialIndexSingle == 2] = 1

nothing = trialSignal[trialIndex == 0, :]
movTrial = trialSignal[trialIndex == 1, :]

# for movement
MovTrialSpeed = trialSpeed[trialIndexSingle == 1, :]

# Colors for plotting
gree = [0.47, 0.67, 0.19]
lgris = [0.8, 0.8, 0.8]
gris = [0.65, 0.65, 0.65]
dgris = [0.3, 0.3, 0.3]
red = [0.78, 0, 0]

pre = 5
post = 25

# permTest_array
perm_alignITI, _ = perm_test_array(ITIsignal, trialSignalinit, 1000)


# bootstrapping
btsrp_nothing = bootstrap_data(nothing, 5000, 0.0001)
btsrp_movTrial = bootstrap_data(movTrial, 5000, 0.0001)
btsrp_movITI = bootstrap_data(ITIsignal, 5000, 0.0001)
btsrp_movInit = bootstrap_data(trialSignalinit, 5000, 0.0001)

ts = np.linspace(-pre, post, trialSignalinit.shape[1])

# plot avg moving tiral signals
plt.plot(ts, np.mean(trialSignalinit, axis = 0),
         color = gree)
plt.fill_between(ts,
                 np.mean(trialSignalinit, axis = 0) + np.std(trialSignalinit, axis=0) / np.sqrt(len(trialSignalinit)),
                 np.mean(trialSignalinit, axis = 0) - np.std(trialSignalinit, axis=0) / np.sqrt(len(trialSignalinit)),
                 color= gree, alpha=0.3)

# plot average ITI signal
plt.plot(ts, np.mean(ITIsignal, axis = 0),
         color = gris)
plt.fill_between(ts,
                 np.mean(ITIsignal, axis = 0) + np.std(ITIsignal, axis=0) / np.sqrt(len(ITIsignal)),
                 np.mean(ITIsignal, axis = 0) - np.std(ITIsignal, axis=0) / np.sqrt(len(ITIsignal)),
                 color= gris, alpha=0.3)

# plot bootstraps
ylinemax = max(np.mean(trialSignalinit, axis =1))
ylinemin = min(np.mean(trialSignalinit, axis =1))

tmp = np.where(btsrp_movInit[1,:]<0)[0]
if len(tmp) >1:
    id = tmp[consec_idx(tmp, thres)]
    plt.plot(ts[id], ylinemax * np.ones((len(ts[id]), 2))+0.5, 's', markersize=7, markerfacecolor= gree, color=gree)
    
tmp = np.where(btsrp_movInit[0,:]>0)[0]
if len(tmp) >1:    
    id = tmp[consec_idx(tmp, thres)]
    plt.plot(ts[id], ylinemax * np.ones((len(ts[id]), 2))+0.5, 's', markersize=7, markerfacecolor= gree, color=gree)

tmp = np.where(btsrp_movITI[1,:]<0)[0]
if len(tmp) >1:
    id = tmp[consec_idx(tmp, thres)]
    plt.plot(ts[id], ylinemax * np.ones((len(ts[id]), 2)), 's', markersize=7, markerfacecolor= gris, color=gris)
    
tmp = np.where(btsrp_movITI[0,:]>0)[0]
if len(tmp) >1:
    id = tmp[consec_idx(tmp, thres)]
    plt.plot(ts[id], ylinemax * np.ones((len(ts[id]), 2)), 's', markersize=7, markerfacecolor= gris, color=gris)

# plot permutations
tmp = np.where(perm_alignITI<0.05)[0]
id = tmp[consec_idx(tmp, thres)]
plt.plot(ts[id], ylinemin * np.ones((len(ts[id]), 2))-0.5, 's', markersize=7, markerfacecolor= red, color=red)



plt.box(False)
plt.title('ZI signal during \nmovement and ITI')
plt.ylabel('Z-Score')
plt.xlabel('Time')
plt.show()

# id = tmp[consec_idx(tmp, thres)]
# plt.plot(ts[id], (130)*np.ones(np.size(ts[id],2)), 's', 'MarkerSize', 5, 'MarkerFaceColor',red,'Color', red)

# # Bootstrap data
# btsrp['nothing'] = bootstrap_data(nothing, 5000, 0.0001)
# btsrp['movTrial'] = bootstrap_data(movTrial, 5000, 0.0001)

