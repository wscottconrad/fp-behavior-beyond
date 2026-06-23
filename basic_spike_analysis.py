# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:32:16 2026

@author: conrad
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
# from scipy.stats import sem
# from scipy.spatial.distance import mahalanobis
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import pandas as pd
# from sklearn.model_selection import cross_val_score
# import pickle
from perm_test_array import perm_test_array
from consec_idx import consec_idx
from sklearn import preprocessing



# -----------------------------
# Parameters
# -----------------------------

# user input

plot_trial_heat_maps = True

plot_stim = True
plot_resp = True

analyze_modded = False # only analyze neurons that showed significant zeta results for stimulus onset or action onset (freezing during stim or ITI)
plot_fraction_mod = False

remove_units_low_variance = False

perm_testing = True # default true
sqrt_transform = True # default true


# summary_data = 'Baseline_corrected'
summary_data = 'Average'
# summary_data = 'Z-score'



# parameter setup
general_baseline = 50
analysis_baseline = 25 # 50 = 1 seconds, try 25 to match MD analysis
all_trial_SD = True 
labelsX = np.arange(-1, 4, 1)

alpha = 0.05 # threshold for modulated neurons from zeta test, should be 0.05
sigma = 1.0  # smoothing, in bins

var_thresh = 10**-10

spike_path = r'W:\Haak\Innate_defense\Data_analysis\22.35.02\Scotts_analysis\collected_spikes.pickle'

def plot_trial_summary(
    trial_matrices,
    trial_types_to_plot,
    alignment_name,
    summary_data,
    labelsX,
    plot_fraction_mod,
    perm_testing=True
):


    
    for region in trial_matrices:
        
        if plot_fraction_mod:
            mod_list = ['pos_mod', 'neg_mod', 'no_mod']
        else:
            mod_list = ['no_mod']
        
        for index, mod in enumerate(mod_list):
        
            plot_data = []
            plt.figure()
            
            max_val = [] # just for plotting
            
            for trial_type in trial_matrices[region]:
    
                if trial_type not in trial_types_to_plot:
                    continue
    
                plt_color = 'grey' if 'nonfreeze' in trial_type else 'red'
                
                if plot_fraction_mod:
                    mod_selection = trial_matrices[region][trial_type][mod]
                    
                    plt.title(
                        f'Average ± SEM spiking activity of {mod_list[index]}, {region}, {alignment_name}'
                    )
                    
                else:
                    mod_selection = list(range(0,len(trial_matrices[region][trial_type][summary_data])))
                    
                    plt.title(
                        f'Average ± SEM spiking activity of {region}, {alignment_name}'
                    )
                
                
                data = trial_matrices[region][trial_type][summary_data][mod_selection]
    
                plot_data.append(data)
    
                mean_trace = np.mean(data, axis=0)
                sem_trace = np.std(data, axis=0) / np.sqrt(data.shape[0])
    
                x = np.arange(data.shape[1])
    
                plt.plot(x, mean_trace, color=plt_color)
    
                plt.fill_between(
                    x,
                    mean_trace + sem_trace,
                    mean_trace - sem_trace,
                    color=plt_color,
                    alpha=0.3
                )
                
                max_val.append(np.max(mean_trace + sem_trace))
    
            if perm_testing and len(plot_data) == 2:
    
                print(f'Permuting {region}...')
    
                perm_test, _ = perm_test_array(
                    plot_data[0],
                    plot_data[1],
                    1000
                )
    
                perm_hits = np.where(perm_test < 0.05)[0]
                perm_x = perm_hits[consec_idx(perm_hits, 1)]
    
                if len(perm_x) > 0:
                    plt.plot(
                        perm_x,
                        np.full(len(perm_x), np.max(max_val) + 0.01),
                        's',
                        markersize=7,
                        markerfacecolor=[0.659, 0.427, 0.91],
                        color=[0.659, 0.427, 0.91]
                    )
    
            plt.axvline(general_baseline, linestyle='--', color='black')
    
            plt.xticks(np.arange(0, 250, general_baseline), labelsX)
            plt.xlabel('Time (s)')
           
    
            plt.show()
            
            

collected_spikes = pd.read_pickle(spike_path)

trial_matrices = {}

for region in collected_spikes:
    
    zeta_results = collected_spikes[region]['zeta_results']
    zeta_results = zeta_results[~np.isnan(zeta_results['stim_freeze_p'])] # filter out neurons with only 2 or less trials
    
    if analyze_modded:
        
        # find which neuron were modulated in any condition
        mod_index = np.where((zeta_results.iloc[:,2:7] < alpha).any(axis=1))[0]
        
    else:
        mod_index = list(range(0, len(zeta_results)))
    
    trial_matrices[region] = {}
    
    for trial_type in collected_spikes[region]:
        
        if 'ITI' in trial_type or trial_type == 'zeta_results':
            continue
        
        n_neurons = len(mod_index)
        n_timepoints = collected_spikes[region][trial_type][0].shape[1]

        Z_trial_matrix = np.zeros((n_neurons, n_timepoints))
        avg_trial_matrix = np.zeros((n_neurons, n_timepoints))
        blc_avg_trial_matrix = np.zeros((n_neurons, n_timepoints))
        
        pos_mod = []
        neg_mod = []
        cmplx_mod = []
        no_mod = []
        
        remove_mask = []
        
        
        filtered_spikes = [
            collected_spikes[region][trial_type][i]
            for i in mod_index
        ]
        
        # baseline_means_collected = np.full((n_neurons, 7), np.nan) # oops shouldnt hard code

        for idx, neuron in enumerate(filtered_spikes):
            
            
            if sqrt_transform:
                neuron = np.sqrt(neuron)
            
            neuron = gaussian_filter1d(
                                neuron,
                                sigma=sigma,
                                axis=1
                            )
            baseline_means = np.mean(neuron[:, :analysis_baseline],axis = 1)
            
            # baseline_means_collected[idx, :len(baseline_means)] = np.mean(neuron[:baseline],axis = 1)
            
            baseline_means = baseline_means.reshape(-1,1)
            
            if all_trial_SD:
                baseline_SD = np.std(neuron[:, :analysis_baseline]) 
            else:
                baseline_SD = np.std(neuron[:, :analysis_baseline], axis = 1).reshape(-1,1)
            
                
            if baseline_SD < var_thresh and remove_units_low_variance:
                print(f'Neuron with low baseline variance found! removing neuron {idx}')
                
                # remove_mask.append(idx)
                
                post_event_SD = np.std(neuron[:, analysis_baseline:])
                
                if post_event_SD < var_thresh:
                    print('This neuron also has low variance after event')
                    
                continue
                
                
                
            Z_trial_matrix[idx] = np.mean((neuron - baseline_means)/ baseline_SD, axis = 0)
            avg_trial_matrix[idx] =  np.mean(neuron, axis = 0)
            blc_avg_trial_matrix[idx] =  np.mean(neuron - baseline_means, axis = 0)
            
            
            # mod_threshold = baseline_SD
            # post_event_activity = np.mean(neuron[:,general_baseline:general_baseline*2] - baseline_means)
            
            mod_threshold = baseline_SD*2
            post_event_activity = np.mean(neuron[:,general_baseline:general_baseline*2] - baseline_means, axis = 0)
            
            
            if np.any(post_event_activity > mod_threshold):
                
                if np.any(post_event_activity < -mod_threshold):
                    cmplx_mod.append(idx)
                else:
                    pos_mod.append(idx)
                    
            elif np.any(post_event_activity < -mod_threshold):
                neg_mod.append(idx)
            else:
                no_mod.append(idx)
                
            
            
            
            if np.mean(np.mean(neuron - baseline_means, axis = 0)[:analysis_baseline]) < -0.1:
                print(f'smth weird happened with neuron {idx}') # i dont think this is needed anymore, was an issue with an earlier bug when i was indexing baseline incorrectly
    
        
        # data_sets = {
        #     "Z_score": Z_trial_matrix,
        #     "Average": avg_trial_matrix,
        #     "Baseline_corrected": blc_avg_trial_matrix,
        # }
        
        trial_matrices[region][trial_type] = {
           "Z_score": Z_trial_matrix,
           "Average": avg_trial_matrix,
           "Baseline_corrected": blc_avg_trial_matrix,
           'pos_mod': pos_mod,
           'neg_mod': neg_mod,
           'no_mod': no_mod,
           'cmplx_mod': cmplx_mod,
           'fraction' : [len(pos_mod), len(neg_mod), len(no_mod), len(cmplx_mod)]           
       }
        
        if plot_trial_heat_maps:
            for name, data in trial_matrices[region][trial_type].items():
                
                if 'mod' in name or name == 'fraction':
                    continue
                
                
                fig, ax = plt.subplots(figsize=(8, 6))
    
                
                
                # Sort rows (neurons/components) by descending mean activity post event
                if 'align' not in trial_type :
                    sort_idx = np.argsort(np.mean(data[:,general_baseline:general_baseline+100], axis=1))[::-1]
                    ax.set_ylabel("Neuron (sorted by mean activity, 2 seconds after event)")

                    
                else:
                    # trying normalize step
                    if name != 'Z_score':
                        data = preprocessing.normalize(data, norm = 'max') # scales to unit norm
                        
                  
                    
                    # Time of maximum activity for each neuron
                    peak_times = np.argmax(data, axis=1)
                    
                    # Sort neurons by when they reach their maximum activity
                    sort_idx = np.argsort(peak_times)
                    
                    
                    # time = np.arange(data.shape[1])
                    # center_of_mass = (data * time).sum(axis=1) / (data.sum(axis=1) + 1e-12)
                    # sort_idx = np.argsort(center_of_mass)  
                    
                    ax.set_ylabel("Neuron (sorted by max normalized activity)")

                    
                sorted_data = data[sort_idx]
                
                im = ax.imshow(
                    sorted_data,
                    aspect="auto",
                    cmap="viridis",
                    interpolation= None
                )
               
                
                ax.set_title(f"{region} | {trial_type} | {name}")
                ax.set_xlabel("Time (s)")
                
                plt.xticks(np.arange(0, 250, general_baseline), labelsX)
                plt.axvline(general_baseline, color = 'white', linestyle = '--')
                plt.colorbar(im, ax=ax, label=name)
                plt.tight_layout()
                plt.show()
 
if plot_fraction_mod:
    for region in trial_matrices:
        
        
        fig, axs = plt.subplots(2, 2)

        for ax, trial_type in zip(axs.flat, trial_matrices[region]):
        
            fractions = [x/sum(trial_matrices[region][trial_type]['fraction']) for x in trial_matrices[region][trial_type]['fraction']]
        
            ax.pie(
                fractions,
                labels=['Excited', 'Suppressed', 'Unchanged', 'Complex'],
                autopct='%1.1f%%'
            )
        
            ax.set_title(trial_type)
        
        plt.tight_layout()
            
 
# if plot_fraction_mod:
             
if plot_stim:
    plot_trial_summary(
        trial_matrices,
        [
            'all_neuron_trials_freeze',
            'all_neuron_trials_nonfreeze'
        ],
        'stimulus aligned',
        summary_data,
        labelsX,
        plot_fraction_mod,
        perm_testing,
        
    )

if plot_resp:
    plot_trial_summary(
        trial_matrices,
        [
            'all_neuron_trials_freeze_aligned',
            'all_neuron_trials_nonfreeze_shuff'
        ],
        'response aligned',
        summary_data,
        labelsX,
        plot_fraction_mod,
        perm_testing,
        
    )

# # stimulus aligned
# if plot_stim:
#     for region in trial_matrices:
        
    
        
#         plot_data = []
#         plt.figure()
#         for idx, trial_type in enumerate(trial_matrices[region]):
            
#             if trial_type != 'all_neuron_trials_freeze' and trial_type != 'all_neuron_trials_nonfreeze':
#                 continue
            
#             if 'nonfreeze' in trial_type:
#                 plt_color = 'grey'
#             else:
#                 plt_color = 'red'
                
                
#             data = trial_matrices[region][trial_type][summary_data]
            
#             plt.plot(range(0,data.shape[1]), np.mean(data, axis = 0), color = plt_color)
            
#             plt.fill_between(range(0,data.shape[1]),
#                                 ( np.mean(data, axis = 0) + np.std(data, axis=0) / np.sqrt(len(data))),
#                                  (np.mean(data, axis = 0) - np.std(data, axis=0) / np.sqrt(len(data))),
#                                  color = plt_color, alpha=0.3)
        
#         if perm_testing and len(plot_data) == 2:
#             print('Permuting ...')
#             # ymax = ymax + 0.5
#             ts = np.arange(data.shape[1])

#             perm_test, _ = perm_test_array(
#                 plot_data[0],
#                 plot_data[1],
#                 1000
#             )
        
#             perm_hits = np.where(perm_test < 0.05)[0]
#             perm_x = perm_hits[consec_idx(perm_hits, 1)]
        
#             plt.plot(
#                 ts[perm_x],
#                 0.03 * np.ones((len(ts[perm_x]), 2)),
#                 's',
#                 markersize=7,
#                 markerfacecolor=[0.659, 0.427, 0.91],
#                 color=[0.659, 0.427, 0.91]
#             )
        
#         plt.axvline(50, linestyle = '--', color = 'black')
    
#         plt.xticks(np.arange(0, 250, 50), labelsX)
#         plt.xlabel('Time (s)')
#         plt.title(f'Average ± SEM spiking activity of {region}, stimulus aligned')
#         plt.show()
   
    
# # response aligned
# if plot_resp:

#     for region in trial_matrices:
#         plot_data = []
#         plt.figure()
#         for idx, trial_type in enumerate(trial_matrices[region]):
       
#             if trial_type != 'all_neuron_trials_freeze_aligned' and trial_type != 'all_neuron_trials_nonfreeze_shuff':
#                 continue
            
#             if 'nonfreeze' in trial_type:
#                 plt_color = 'grey'
#             else:
#                 plt_color = 'red'
                
#             data = trial_matrices[region][trial_type][summary_data]
            
#             plot_data.append(data)

            
#             plt.plot(range(0,data.shape[1]), np.mean(data, axis = 0), color = plt_color)
            
#             plt.fill_between(range(0,data.shape[1]),
#                                 ( np.mean(data, axis = 0) + np.std(data, axis=0) / np.sqrt(len(data))),
#                                  (np.mean(data, axis = 0) - np.std(data, axis=0) / np.sqrt(len(data))),
#                                  color = plt_color, alpha=0.3)
        
#         if perm_testing and len(plot_data) == 2:
#             print('Permuting ...')
#             # ymax = ymax + 0.5
#             ts = np.arange(data.shape[1])

#             perm_test, _ = perm_test_array(
#                 plot_data[0],
#                 plot_data[1],
#                 1000
#             )
        
#             perm_hits = np.where(perm_test < 0.05)[0]
#             perm_x = perm_hits[consec_idx(perm_hits, 1)]
        
#             plt.plot(
#                 ts[perm_x],
#                 0.03 * np.ones((len(ts[perm_x]), 2)),
#                 's',
#                 markersize=7,
#                 markerfacecolor=[0.659, 0.427, 0.91],
#                 color=[0.659, 0.427, 0.91]
#             )
            
#         plt.axvline(50, linestyle = '--', color = 'black')
    
#         plt.xticks(np.arange(0, 250, 50), labelsX)
#         plt.xlabel('Time (s)')
#         plt.title(f'Average ± SEM spiking activity of {region}, response aligned')
#         plt.show()

