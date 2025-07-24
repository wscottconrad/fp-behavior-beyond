# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 12:52:25 2025

@author: sconrad
"""

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

# Define where the stuff is
tankfolder = r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\'

siteList = ['ZI-L', 'ZI-R', 'SC-L', 'SC-R', 'SC to ZI-L', 'SC to ZI-R']

allDat = [[] for _ in range(17)]
trialDat = {site: {'approach': [], 'avoid': []} for site in siteList}  # Dictionary to separate data by site
trialDat_laser_algn = {site: {'approach': [], 'avoid': [], 'NR': []} for site in siteList}  # Dictionary to separate data by site


# Iterate through experiment folder
files = [f for f in os.listdir(tankfolder) if f.endswith('.pkl') and f != 'allDatComb.pkl']

for file in files:
    with open(os.path.join(tankfolder, file), 'rb') as f:
        data = pickle.load(f)
        
        # Extract site information
        site = data.get('site', 'Unknown')  # Use 'Unknown' if site is missing
        
        # Z-scored data
        allDat[0].append(data['ZdFoF'])
        allDat[4].append(data['ZdFoFITI'])

        
        # Store trial-specific data
        if not np.isscalar(data['ZdFoFApproach']):
            allDat[12].append(data['ZdFoFApproach'])
            allDat[13].append(data['ZdFoFApproach_trialOnset'])
            if site in trialDat:
                trialDat[site]['approach'].append(data['ZdFoFApproach'])
                trialDat_laser_algn[site]['approach'].append(data['ZdFoFApproach_trialOnset'])

        
        if not np.isscalar(data['ZdFoFAvoid']):
            allDat[14].append(data['ZdFoFAvoid'])
            allDat[15].append(data['ZdFoFAvoid_trialOnset'])
            if site in trialDat:
                trialDat[site]['avoid'].append(data['ZdFoFAvoid'])
                trialDat_laser_algn[site]['avoid'].append(data['ZdFoFAvoid_trialOnset'])
        
        # print([data['mouse'], data['session']])    

        if not np.isscalar(data['ZdFoFNR']):
            allDat[16].append(data['ZdFoFNR'])
            if site in trialDat:
                trialDat_laser_algn[site]['NR'].append(data['ZdFoFNR'])

        # Speed data
        if not np.isscalar(data['speedTrialsMov']):
            allDat[3].append(data['speedTrials']) 
            allDat[5].append(data['speedITI'])
            allDat[7].append(data['speedTrialsMov'][data['speedTrialsMov'][:, 0] != 0])

        # Traces for green and isosbestic channels
        allDat[1].append(data['Gdata'])
        allDat[2].append(data['Idata'])
        # allDat[8].append(data['InitGdata'])
        # allDat[9].append(data['InitIdata'])
        allDat[10].append(data['ITIGdata'])
        allDat[11].append(data['ITIIdata'])

# Combine data
allDatComb = {
    'trialSignal': np.vstack(allDat[0]),
    'Gtrace': np.vstack(allDat[1]),
    'Itrace': np.vstack(allDat[2]),
    'trialSpeed': np.concatenate(allDat[3]),
    'ITIsignal': np.vstack(allDat[4]),
    'ITIspeed': np.concatenate(allDat[5]),
    # 'ZdFoFinit': np.vstack(allDat[6]) if [data['mouse'], data['session']] != ['109436', '2025_03_13_'] else np.nan,
    'speedTrialsMov': np.vstack(allDat[7]),
    # 'InitGdata': np.vstack(allDat[8]),
    # 'InitIdata': np.vstack(allDat[9]),
    'ITIGdata': np.vstack(allDat[10]),
    'ITIIdata': np.vstack(allDat[11]),
    'approachSignal': np.vstack(allDat[12]),
    'approachSignal_trialOnset': np.vstack(allDat[13]),
    'avoidSignal': np.vstack(allDat[14]),
    'avoidSignal_trialOnset': np.vstack(allDat[15]),
    'NRsignal': np.vstack(allDat[16]),

    # SITE SPECIFIC DATA
    # movement aligned
    'trialData': {site: {'approach': np.vstack(values['approach']) if values['approach'] else np.array([]),
                        'avoid': np.vstack(values['avoid']) if values['avoid'] else np.array([])}
                  for site, values in trialDat.items()},  # Store separated site data
    # prey laser onset aligned             
    'trialData_trialOnset': {site: {'approach': np.vstack(values['approach']) if values['approach'] else np.array([]),
                        'avoid': np.vstack(values['avoid']) if values['avoid'] else np.array([]),
                        'NR': np.vstack(values['NR']) if values['NR'] else np.array([])}
                 for site, values in trialDat_laser_algn.items()}  # Store separated site data
}

# Save combined data
with open(os.path.join(tankfolder, 'allDatComb.pkl'), 'wb') as f:
    pickle.dump(allDatComb, f)

# Plotting
# traceTiming = np.arange(allDatComb['InitGdata'].shape[1])
# plt.figure()

# # Green signal
# plt.plot(traceTiming, np.mean(allDatComb['InitGdata'], axis=0), color='#64a858')
# plt.fill_between(traceTiming,
#                  np.mean(allDatComb['InitGdata'], axis=0) + np.std(allDatComb['InitGdata'], axis=0),
#                  np.mean(allDatComb['InitGdata'], axis=0) - np.std(allDatComb['InitGdata'], axis=0),
#                  color='#64a858', alpha=0.3)

# # Isosbestic signal
# plt.plot(traceTiming, np.mean(allDatComb['InitIdata'], axis=0), color='#7b8be8')
# plt.fill_between(traceTiming,
#                  np.mean(allDatComb['InitIdata'], axis=0) + np.std(allDatComb['InitIdata'], axis=0),
#                  np.mean(allDatComb['InitIdata'], axis=0) - np.std(allDatComb['InitIdata'], axis=0),
#                  color='#7b8be8', alpha=0.3)

# plt.box(False)
# plt.title('Filtered & Baseline Corrected\nGreen and Isosbestic Signals\nFor Prey Approach')
# plt.ylabel('Z-Score')
# plt.xlabel('Time')
# plt.show()
