# -*- coding: utf-8 -*-
"""
Created on 24/11/25

This code plot summary graphs per animal per site, with the purpose of data cleaning. 

perhaps best if run on data without trial time exclsion?

@author: sconrad
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import os

def build_filename(date, ID, channel):
    return f"{date}{ID} Channel {channel}.pkl"

nt = False
initiate_aligned = True
trial_type = ['fm_approach', 'NR'] #choose one or more trial type here ('approach', 'avoid', 'NR' (nt only), 'IR', 'ITI' )
plot_speed = False

focused_analysis = True

sr = 30
whole_trial = False

debug = True

if nt:
    exp_type = 'nt'
else:
    exp_type = 'fm laser'
    
filePath = 'W:\\Conrad\\Innate_approach\\Data_collection\\24.35.01\\'

r_log = pd.read_csv(f"{filePath}\\recordinglog.csv", sep=None, engine="python", encoding='utf-8-sig')
r_log = r_log[r_log['Exp'] == exp_type]

r_log = r_log[r_log['notes'] != 'no ttl alignment']
r_log = r_log[r_log['added to db?'] != 'neurotar data corrupt']

animal_ids = r_log['ID'].unique()
sites = pd.unique(r_log[['1','2']].values.ravel('K'))

# sites = ['PAG-L', 'PAG-R']

excluded_animals = ['105647', '118401', '118402']

animal_behavior = {str(animal_id): {} for animal_id in animal_ids}
# Load combined data
if nt == True:
    tankfolder = r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\'
else: 
    tankfolder = r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\freelymoving\\'

files = [f for f in os.listdir(tankfolder) if f.endswith('.pkl') and f != 'allDatComb.pkl' and f != 'approach_times_since_trial_start.pkl' ]

thres = round(1/3*sr)  # Consecutive threshold length
pre = 5
post = 25

for ID in animal_ids:
    
    if str(ID) in excluded_animals:
        continue
    
    all_animal_data = r_log[r_log['ID'] == ID]
    
    prey_speeds = {
        'snout':[],
        'hrC': [],
        'tail': []}
    
    IR_speeds = {
        'snout':[],
        'hrC': [],
        'tail': []}
    
    for site in sites:
        
        if site == 'exclude' or site == 'switched':
            continue
        
        if focused_analysis:
            if 'to' in site:
                continue
            
            if 'MLR' in site:
                continue
            
        # if site != 'MLR-R' and site != 'MLR-L':
        #     continue
        
        data__single_animal_site = []
        approach_data = []
        control_data = []


        animal_by_site_data_ch1 = all_animal_data[all_animal_data['1'] == site]
        animal_by_site_data_ch2 = all_animal_data[all_animal_data['2'] == site]
        
        if len(animal_by_site_data_ch1) == 0 and len(animal_by_site_data_ch2) == 0:
            continue

        dates_ch1 = animal_by_site_data_ch1['Date'].unique()
        dates_ch2 = animal_by_site_data_ch2['Date'].unique()
        
        ID = str(int(ID))
        
        for date in dates_ch1:
            filename = build_filename(date, ID, 1)
            data_path = os.path.join(tankfolder, filename)

            if os.path.exists(data_path):
                with open(data_path, 'rb') as f:
                    temp_data = pickle.load(f)
                data__single_animal_site.append(temp_data)
            

        for date in dates_ch2:
            filename = build_filename(date, ID, 2)
            data_path = os.path.join(tankfolder, filename)

            if os.path.exists(data_path):
                with open(data_path, 'rb') as f:
                    temp_data = pickle.load(f)
                data__single_animal_site.append(temp_data)
    
    
        for item in data__single_animal_site:
            if initiate_aligned:
                approach_data.append(item['ZdFoFApproach']) if not np.any(np.isnan(item['ZdFoFApproach'])) else True
                
            else:
                approach_data.append(item['ZdFoFApproach_trialOnset']) if not np.any(np.isnan(item['ZdFoFApproach_trialOnset'])) else True
            
            if nt:
                if trial_type[1] == 'NR' and initiate_aligned:
                    control_data.append(item['ZdFoFNR_yoked']) if not np.any(np.isnan(item['ZdFoFNR_yoked'])) else True
                
                elif trial_type[1] == 'NR' and not initiate_aligned:
                     control_data.append(item['ZdFoFNR']) if not np.any(np.isnan(item['ZdFoFNR'])) else True
                    
                elif trial_type[1] == 'ITI':
                    control_data.append(item['ZdFoFITI']) if not np.any(np.isnan(item['ZdFoFITI'])) else True

            else:
                control_data.append(item['IR_ZdFoFApproach']) if not np.any(np.isnan(item['IR_ZdFoFApproach'])) else True
                
                if not pd.isna(item['approach_speeds']):
                    for key, speeds in item['approach_speeds'].items():
                        prey_speeds[key].append(speeds)
                        
                    for key, speeds in item['IR_speeds'].items():
                        IR_speeds[key].append(speeds)

        

        
        if len(approach_data) > 0:
            approach_data = np.vstack(approach_data)  
         
        if len(trial_type) > 1:
            if len(control_data) > 0:
                control_data = np.vstack(control_data)
            plot_data = [approach_data, control_data]
        else:
            plot_data = [approach_data]

        
        plt.figure()
        for index, signal in enumerate(plot_data):
            
            animal_behavior[ID][trial_type[index]] = len(signal)
            
            if len(signal) == 0:
                print(f"No {trial_type[index]} trials for {ID} {site}\n")
                continue
            
            print(f"{len(signal)} {trial_type[index]} trials\n")
            # # Bootstrapping
            # print('bootstrapping ...')
            # btsrp_app = bootstrap_data(signal, 10000, 0.0001)
            
            # Colors for plotting
            if trial_type[index] == 'approach':
                plt_color = [0.47, 0.67, 0.19] #green
            elif trial_type[index] == 'NR':
                plt_color = [0.65, 0.65, 0.65] #grey
            elif trial_type[index] == 'ITI' or trial_type[index] == 'IR':
                plt_color = [0.416, 0.741, 0.741] #blue
            else:
                plt_color = [0.78, 0, 0] # red, avoid
            
            ts = np.linspace(-pre, post, signal.shape[1])
            plt.plot(ts, np.mean(signal, axis=0), color=plt_color, label=trial_type[index])
            plt.fill_between(ts,
                             np.mean(signal, axis=0) + np.std(signal, axis=0) / np.sqrt(len(signal)),
                             np.mean(signal, axis=0) - np.std(signal, axis=0) / np.sqrt(len(signal)),
                             color=plt_color, alpha=0.3)
         

            
        plt.axvline(x=0, linestyle='--', color='black', linewidth=1.5)
        plt.axhline(y=0, linestyle='--', color='black', linewidth=1.5)
        
        # Hide the top and right spines
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        
        # Set tick parameters to remove right and top ticks
        plt.gca().tick_params(axis='x', which='both', direction='out', bottom=True, top=False)
        plt.gca().tick_params(axis='y', which='both', direction='out', left=True, right=False)
        plt.title(f'{ID}: {site}')
        plt.ylabel('Z-Score')
        if not whole_trial:
            plt.xlim(-5,5)

        if initiate_aligned == True and 'NR' not in trial_type:
            plt.xlabel('Movement onset (s)')
        elif initiate_aligned == True and 'NR' in trial_type:
            plt.xlabel('Shuffled movement onset (s)')
        else:
            plt.xlabel('Prey Laser onset (s)')
        # plt.legend()
        # plt.savefig(f'{tankfolder}{site}_{trial_type}.png', transparent = True)
        plt.show()
        
    if not nt and plot_speed:
        
        # if site != 'MLR-R' and site != 'MLR-L':
        #     continue
        
        ncols = 3
        nrows = 1

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(5 * ncols, 4 * nrows),
            sharex=True,
            constrained_layout=True
        )
        
        ts = np.linspace(-pre, post, (pre+post)*sr) 
        # Flatten axes so you can index them with ax_idx
        axes = np.atleast_1d(axes).ravel()
        ax_idx = 0
        
        for body_part, speed_data in prey_speeds.items():
            ax = axes[ax_idx]
            ax_idx += 1
                
            
            arrays = [a for a in prey_speeds[body_part] if len(a) > 0] # filters for when no approach was made
            if len(arrays)==0:
                continue
            prey_data = np.vstack(arrays)
            
            arrays = [a for a in IR_speeds[body_part] if len(a) > 0] # filters for when no approach was made
            if len(arrays)==0:
                continue
            IR_data = np.vstack(arrays)
            
            # Prey
            ax.plot(ts, np.mean(prey_data, axis=0), color=[0.47, 0.67, 0.19], label='Prey')
            ax.fill_between(ts,
                             np.mean(prey_data, axis=0) + np.std(prey_data, axis=0) / np.sqrt(len(prey_data)),
                             np.mean(prey_data, axis=0) - np.std(prey_data, axis=0) / np.sqrt(len(prey_data)),
                             color=[0.47, 0.67, 0.19], alpha=0.3)
            
            # IR
            ax.plot(ts, np.mean(IR_data, axis=0), color=[0.416, 0.741, 0.741], label='IR')
            ax.fill_between(ts,
                             np.mean(IR_data, axis=0) + np.std(IR_data, axis=0) / np.sqrt(len(IR_data)),
                             np.mean(IR_data, axis=0) - np.std(IR_data, axis=0) / np.sqrt(len(IR_data)),
                             color=[0.416, 0.741, 0.741], alpha=0.3)
            
            ax.set_title(f'{ID}: {body_part} speed')
            ax.tick_params(axis='x', which='both', direction='out', bottom=True, top=False)
            ax.tick_params(axis='y', which='both', direction='out', left=True, right=False)
            ax.set_xlim(-5,5)
            ax.axvline(x=0, linestyle='--', color='black', linewidth=1.5)
            ax.set_xlabel('Movement onset (s)')
            ax.set_ylabel('Speed (m/s)')
            
        plt.show()

