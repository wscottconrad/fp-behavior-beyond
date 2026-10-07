# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 15:45:00 2025

@author: sconrad
"""

import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
import pickle
from perm_test_array import perm_test_array
# from bootstrap_data import bootstrap_data
from bootstrap_data_bias_corrected import bootstrap_data
from consec_idx import consec_idx
from scipy.stats import ttest_rel
from scipy.stats import shapiro
import scipy.stats as stats
from boxoff import boxoff
import math
# import statsmodels.formula.api as smf
from scipy.ndimage import median_filter





nt = False
initiate_aligned = True
perm_testing = True # if u wanna test difference between two signals
# site_specific = 'ZI-L' 
trial_type = ['approach', 'IR'] #choose one or more trial type here ('approach', 'avoid', 'NR' (nt only), ITI, IR )
correlation_list = ['approach', 'IR']

combine_hemispheres = True


plot_signal = True
calculate_auc = False
plot_speed = True
plot_correlations = False

plot_site = False # plot speed by site

filter_speed = True 
auto_ymax = True

my_ylim = [-0.75, 1.75]


focused_analysis = True # just for SC, ZI, and PAG. no projections
whole_trial = False # scaling. either whole trial, or if false +-5 seconds

sr = 30 # sampling rate, fps

def safe_concatenate(left, right, trial_length):
    arrays = [
        x[:, :trial_length]
        for x in (left, right)
        if x is not None
    ]

    if len(arrays) == 0:
        return None
    elif len(arrays) == 1:
        return arrays[0]
    else:
        return np.concatenate(arrays)

# Load combined data
if nt == True:
    tankfolder = r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\'
else: 
    tankfolder = r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\freelymoving\\'

with open(f'{tankfolder}allDatComb.pkl', 'rb') as f:
    d = pickle.load(f)

if initiate_aligned == True:
    trialData = d['trialData']  # for movement aligned data
else: 
    trialData = d['trialData_trialOnset']  # for prey laser onset data
 


thres = round(1/3*sr)  # Consecutive threshold length
pre = 5

if whole_trial:
    post = 25
else:
    post = 5
    
trial_length = sr*(pre+post)


data_list = [trialData]

if combine_hemispheres:
    combined_list = [{}, {}]
    data_list = [trialData, d['ITI']]
    
    
    for x, data_set in enumerate(data_list):
        for key in data_set.keys():
            # Only process "-L" entries to avoid duplicates
            if key.endswith('-L'):
                base = key[:-2]  # e.g. "PAG" from "PAG-L"
                left = data_set.get(f"{base}-L")
                right = data_set.get(f"{base}-R")
    
                if left is not None and right is not None:
                    combined_list[x][f"{base}-both"] = {}
    
                    for k in left.keys():
                        left_val = left[k]
                        right_val = right.get(k)
    
                        if right_val is None:
                            combined_list[x][f"{base}-both"][k] = left_val
                        elif left_val.size == 0:
                            combined_list[x][f"{base}-both"][k] = right_val
                        elif right_val.size == 0:
                            combined_list[x][f"{base}-both"][k] = left_val
                        else:
                            combined_list[x][f"{base}-both"][k] = np.concatenate(
                                [left_val, right_val]
                            )
                else:
                    # Handle missing side
                    combined_list[x][f"{base}-both"] = left if left is not None else right
                    
        data_list[x] = combined_list[x]
        
    d['ITI'] = data_list[1]


trialData = data_list[0]

ts = np.linspace(-pre, post, (pre+post)*sr) 

n_sites = len(trialData.keys())
ncols = 3

if focused_analysis:
    nrows = 1
    
else:
    nrows = math.ceil(n_sites / ncols)

if plot_signal:    

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5 * ncols, 4 * nrows),
        sharex=True,
        constrained_layout=True
    )
    
    # Flatten axes so you can index them with ax_idx
    axes = np.atleast_1d(axes).ravel()
    
    # Handle case of only one subplot
    if n_sites == 1:
        axes = [axes]
    
    ax_idx = 0
    
    for site, data in trialData.items():
        
        if focused_analysis:
            if 'to' in site: #skip projections
                continue
            
            if site == 'MLR-both':
                continue
       
        ax = axes[ax_idx]
        ax_idx += 1
            
    # if site == site_specific:
        if 'ITI' not in trial_type and 'IR' not in trial_type:
            if len(trial_type) > 1:
                plot_data = list(range(len(trial_type)))
                for i, signal in enumerate(trial_type):
                    plot_data[i] = data[trial_type[i]][~np.isnan(data[trial_type[i]]).any(axis=1)]
                    print(f"Total {trial_type[i]} trials for {site}: {len(data[trial_type[i]])}")
            else:
                if len(data[trial_type[0]]) > 0:
                    plot_data = [data[trial_type[0]][~np.isnan(data[trial_type[0]]).any(axis=1)]]
                    print(f"Total trials for {site}: {len(data[trial_type[0]])}")
                
        if 'ITI' in trial_type or 'IR' in trial_type:
            if len(trial_type) > 1:
                if len(data[trial_type[0]]) > 0:
                    signal = data[trial_type[0]][~np.isnan(data[trial_type[0]]).any(axis=1)]
                    
                    print(f"Total trials for {site}: {len(data[trial_type[0]])}")
                    
                    if 'ITI' in trial_type:
                        signal_ITI = d['ITI'][site]['ITI']
                        if np.any(np.isnan(signal_ITI)):
                            print('NAN detected in ITI')
                            
                        mask = np.all(signal_ITI <= 200, axis=1)
                        signal_ITI = signal_ITI[mask] # remove weird noise condition
                            
                        print(f"Total ITI traces for {site}: {len(signal_ITI)}")
                        plot_data = [signal, signal_ITI]
                    else:
                        signal_IR = data[trial_type[1]][~np.isnan(data[trial_type[1]]).any(axis=1)]
                        print(f"Total IR traces for {site}: {len(signal_IR)}")
                        plot_data = [signal, signal_IR]
                    
                    
            else:
                plot_data = [d['ITI'][site]['ITI']]
                print(f"Total ITI traces for {site}: {len(plot_data[0])}")
    
        
        
        # finding values to plot bCI and perm test data later...
        
        if len(plot_data) > 1:
            
            ymax_total = np.max([np.mean(plot_data[0], axis = 0), np.mean(plot_data[1], axis = 0)]) + 1
        else:
            ymax_total = np.max([np.mean(plot_data[0], axis = 0)]) + 1 + 2
    
        
        for index, signal in enumerate(plot_data):
        
            # Bootstrapping
            print('bootstrapping ...')
            signal = signal[:, :(pre+post)*sr]
            
            btsrp_app = bootstrap_data(signal, 10000, 0.0001)
            
            # Colors and settings for plotting
            y_adj = 0.75
            
            if trial_type[index] == 'approach':
                plt_color = [0.47, 0.67, 0.19] #green
                y_adj = 1
            elif trial_type[index] == 'NR' or trial_type[index] == 'ITI':
                plt_color = [0.65, 0.65, 0.65] #grey
            elif trial_type[index] == 'IR':
                plt_color = [0.416, 0.741, 0.741] #blue
            else:
                plt_color = [0.78, 0, 0] # red, avoid
              
            
            # Plot signal
            n_trials = len(signal)
            ax.plot(ts, np.nanmean(signal, axis=0), color=plt_color, label=f'{trial_type[index]}, (n={n_trials})')
            ax.fill_between(ts,
                             np.mean(signal, axis=0) + np.std(signal, axis=0) / np.sqrt(len(signal)),
                             np.mean(signal, axis=0) - np.std(signal, axis=0) / np.sqrt(len(signal)),
                             color=plt_color, alpha=0.3)
        
            ymax = ymax_total - index/2
            
            if not auto_ymax:
                ymax = my_ylim[1] - 1
            
            # Bootstrap significance
            tmp = np.where(btsrp_app[1, :] < 0)[0]
            if len(tmp) > 1:
                id = tmp[consec_idx(tmp, thres)]
                ax.plot(ts[id], ymax * np.ones((len(ts[id]), 2)) + y_adj, 's', 
                         markersize=7, markerfacecolor=plt_color, color=plt_color)
                
            tmp = np.where(btsrp_app[0, :] > 0)[0]
            if len(tmp) > 1:
                id = tmp[consec_idx(tmp, thres)]
                ax.plot(ts[id], ymax * np.ones((len(ts[id]), 2)) + y_adj, 's', 
                         markersize=7, markerfacecolor=plt_color, color=plt_color)
                
        
            
        if perm_testing and len(plot_data) == 2:
            print('Permuting ...')
            ymax = ymax + 0.5
        
            perm_test, _ = perm_test_array(
                plot_data[0][:, :(pre+post)*sr],
                plot_data[1][:, :(pre+post)*sr],
                1000
            )
        
            perm_hits = np.where(perm_test < 0.05)[0]
            perm_x = perm_hits[consec_idx(perm_hits, thres)]
        
            ax.plot(
                ts[perm_x],
                ymax * np.ones((len(ts[perm_x]), 2)),
                's',
                markersize=7,
                markerfacecolor=[0.659, 0.427, 0.91],
                color=[0.659, 0.427, 0.91]
            )
            
        ax.axvline(x=0, linestyle='--', color='black', linewidth=1.5)
        ax.axhline(y=0, linestyle='--', color='black', linewidth=1.5)
       
        if not auto_ymax:
            ax.set_ylim(my_ylim)
        
        # Hide the top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Set tick parameters to remove right and top ticks
        ax.tick_params(axis='x', which='both', direction='out', bottom=True, top=False)
        ax.tick_params(axis='y', which='both', direction='out', left=True, right=False)
        ax.set_title(f'{site}')
        ax.set_ylabel('Z-Score')
        ax.legend()
        # ax.set_xlim(-5,5)
        
        if initiate_aligned  and 'NR' not in trial_type:
            x_label = 'Movement onset (s)'
        elif initiate_aligned and 'NR' in trial_type:
            x_label = 'Shuffled movement onset (s)'
        else:
            x_label = 'Prey Laser onset (s)'
        ax.set_xlabel(x_label)
        
    for ax in axes[ax_idx:]: # hide any graphs not filled with data
        ax.set_visible(False)  
    for ax in axes[-ncols:]:
        if ax.get_visible():
            ax.set_xlabel(x_label)
            
    for ax in axes.flat:
        ax.tick_params(labelbottom=True)
            
      
        
    plt.show()
    
  
# %% AUC
for site, data in trialData.items():

    if calculate_auc:
        # AUC comparison
        for index, signal in enumerate(plot_data):
            # Define baseline windows (in frames)
            baseline1 = (-5, -2.5)
            baseline2 = (-2.5, 0)
            
            # Get indices for each baseline
            b1_idx = np.where((ts >= baseline1[0]) & (ts < baseline1[1]))[0]
            b2_idx = np.where((ts >= baseline2[0]) & (ts < baseline2[1]))[0]
            
            # Time step (assumes uniform sampling)
            dt = ts[1] - ts[0]
            
            # Compute AUC per trial
            auc_b1 = np.trapz(signal[:, b1_idx], dx=dt, axis=1)
            auc_b2 = np.trapz(signal[:, b2_idx], dx=dt, axis=1)
            
            # Store if you want later
            # e.g. auc_results[site][trial_type[index]] = (auc_b1, auc_b2)
            
            plt.figure()
            
            # Colors for plotting
            if trial_type[index] == 'approach':
                plt_color = [0.47, 0.67, 0.19] #green
            elif trial_type[index] == 'NR':
                plt_color = [0.65, 0.65, 0.65] #grey
            elif trial_type[index] == 'ITI' or trial_type[index] == 'IR':
                plt_color = [0.416, 0.741, 0.741] #blue
            else:
                plt_color = [0.78, 0, 0] # red, avoid
                
            # Means and SEMs
            means = [np.mean(auc_b1), np.mean(auc_b2)]
            sems  = [
                np.std(auc_b1) / np.sqrt(len(auc_b1)),
                np.std(auc_b2) / np.sqrt(len(auc_b2))
            ]
            
            # X positions
            x = np.arange(2)
            
            # Bar plot
            plt.bar(
                x,
                means,
                yerr=sems,
                color=plt_color,
                alpha=0.8,
                edgecolor='black',
                capsize=5
            )
            for i in range(len(auc_b1)):
                plt.plot(x, [auc_b1[i], auc_b2[i]], color='k', alpha=0.2, linewidth=0.8)
    
            # Formatting
            plt.xticks(x, ['-5 to -2.5 s', '-2.5 to 0 s'])
            plt.ylabel('Baseline AUC')
            plt.title(f'{site} | {trial_type[index]}')
            
            # Match your axis style
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.gca().tick_params(axis='x', which='both', direction='out')
            plt.gca().tick_params(axis='y', which='both', direction='out')
            
            boxoff()
            plt.show()
            
            diff = auc_b2 - auc_b1
            plt.hist(diff, bins=20)
            plt.title('AUC difference distribution')
            plt.show()
            
            stats.probplot(diff, plot=plt)
            plt.show()
            print(f'Normality test: {shapiro(diff)}')
            
            tstat, pval = ttest_rel(auc_b1, auc_b2)
            print(f"  Paired t-test: t = {tstat:.3f}, p = {pval:.4e}")
            
            # df = pd.DataFrame({
            #     'diff_auc': auc_b2 - auc_b1,
            #     'animal': animal_ids   # same length as auc arrays
            # })
            
            # model = smf.mixedlm("diff_auc ~ 1", df, groups=df["animal"])
            # res = model.fit()
            
            # print(res.summary())
            
# %% SPEED
            
# PLOTTING SPEED BY SITE
if plot_speed:
    
    control_key = 'ITIspeed'
    
    
    if nt:
        approach_key = 'speedTrialsMov'
        correction_factor = control_correction_factor = 1000
    else:
        approach_key = 'approach_snout_speed'
        correction_factor = 1

        if 'IR' in trial_type:
            control_key = 'IR_snout_speed'
            control_correction_factor = 1

        if 'ITI' in trial_type:
            control_correction_factor = 1000


        
    speedData = d['speedData']
    
    datasets = []
    
    if not combine_hemispheres:
        for site, speed in speedData.items():
        
            if not nt:
                
                if len(speed[approach_key]) > 0:
                    approach_speed = speed[approach_key][:,:trial_length] / correction_factor
                    
                else: 
                    approach_speed = []
        
                if len(speed[control_key]) > 0:
                    control_speed = speed[control_key][:,:trial_length] / control_correction_factor
                else:
                    control_speed = []
    
        
            else:
                if len(speed[approach_key]) == 0 or len(speed[control_key]) == 0:
                    continue
        
                approach_speed = speed['speedTrialsMov'][:,:trial_length] / correction_factor
                control_speed = speed['ITIspeed'][:,:trial_length] / control_correction_factor
        
            datasets.append((site, approach_speed, control_speed))
            
    else:
      
        base_site = []
        
        for site, speed in speedData.items():
                # Only process "-L" entries to avoid duplicates
                if site.endswith('-L'):
                    base_site.append(site[:-2])  # e.g. "PAG" from "PAG-L"
            
        for base in base_site:
                left_approach = speedData[f"{base}-L"][approach_key]/correction_factor if len(speedData[f"{base}-L"][approach_key]) > 0 else None 
                right_approach = speedData[f"{base}-R"][approach_key]/correction_factor if len(speedData[f"{base}-R"][approach_key]) > 0 else None 
                
                left_control = speedData[f"{base}-L"][control_key]/control_correction_factor if len(speedData[f"{base}-L"][control_key]) > 0 else None 
                right_control = speedData[f"{base}-R"][control_key]/control_correction_factor if len(speedData[f"{base}-R"][control_key]) > 0 else None 
                
                
                approach_speed = safe_concatenate(
                    left_approach, right_approach, trial_length
                )
        
                control_speed = safe_concatenate(
                    left_control, right_control, trial_length
                )
                
                datasets.append((base, approach_speed, control_speed))


    
    
    # ----------------------------------------------------------
    # Combine all sites if requested
    # ----------------------------------------------------------
    if not plot_site and not combine_hemispheres:
    
        all_approach = np.concatenate([d[1] for d in datasets], axis=0)
        all_baseline = np.concatenate([d[2] for d in datasets if len(d[2]) > 0], axis=0)
    
        datasets = [("All sites", all_approach, all_baseline)]
    
    n_sites = len(datasets)
    
    ncols = 3 if plot_site or focused_analysis else 1
    nrows = 1 if focused_analysis else math.ceil(n_sites / ncols)
  
 
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5*ncols, 4*nrows),
        sharex=True,
        constrained_layout=True
    )
    
    axes = np.atleast_1d(axes).ravel()
    
    ax_idx = 0
    
    for ax, (site, approach_speed, control_speed) in zip(axes, datasets):
        
        if focused_analysis:
            if 'to' in site: #skip projections
                continue
            
            if site == 'MLR-both':
                continue
            
        ax = axes[ax_idx]
        ax_idx += 1
            
        if len(control_speed) > 0:
            
            control_speed = control_speed[:, :(pre+post)*sr]

            if filter_speed:
                control_speed = median_filter(control_speed, 5)
                
            n_trials = len(control_speed)
            if 'ITI' in trial_type:
        
                label = f'ITI (n = {n_trials})' 
            
            elif 'IR' in trial_type:
                
                label = f'IR (n = {n_trials})'
                
            
        
            ax.plot(
                ts,
                np.nanmean(control_speed, axis=0),
                color=[0.6, 0.6, 0.6],
                label=label
            )
    
            ax.fill_between(
                ts,
                np.nanmean(control_speed, axis=0)
                + np.nanstd(control_speed, axis=0)/np.sqrt(len(control_speed)),
                np.nanmean(control_speed, axis=0)
                - np.nanstd(control_speed, axis=0)/np.sqrt(len(control_speed)),
                color=[0.6,0.6,0.6],
                alpha=0.3
            )
    
    
        # approach signal 
        n_trials = len(approach_speed)
        ax.plot(
            ts,
            np.nanmean(approach_speed, axis=0),
            color=[0.78,0,0],
            label=(f'Approach (n = {n_trials})')
        )
    
        ax.fill_between(
            ts,
            np.nanmean(approach_speed, axis=0)
            + np.nanstd(approach_speed, axis=0)/np.sqrt(len(approach_speed)),
            np.nanmean(approach_speed, axis=0)
            - np.nanstd(approach_speed, axis=0)/np.sqrt(len(approach_speed)),
            color=[0.78,0,0],
            alpha=0.3
        )
    
        ymin, ymax = ax.get_ylim()
        ax.vlines(0, ymin, ymax, linestyle='--', color='black')
    
        ax.set_title(site)
        ax.set_ylabel("Speed (m/s)")
        ax.set_xlabel("Time, approach onset (s)")
        ax.legend()
        boxoff()
    
    for ax in axes[len(datasets):]:
        ax.set_visible(False)
    
    plt.show()
    

    
# %% Cross correlations

if plot_correlations:
    
    corrData = d['lag correlations']
    combined_corr = {}

        
        # --------------------------------------------------------
        # Combine L/R speed data into -both
        # --------------------------------------------------------

    for key in corrData.keys():

        if not key.endswith('-L'):
            continue

        base = key[:-2]

        left = corrData.get(f'{base}-L')
        right = corrData.get(f'{base}-R')

        if left is None and right is None:
            continue

        combined_corr[f'{base}-both'] = {}
        
        for correlation_type in correlation_list:
            
            corr_key = f'lag_correlations_{correlation_type}'
        

            left_val = (
                left.get(correlation_type, np.array([]))
                if left is not None else np.array([])
            )

            right_val = (
                right.get(correlation_type, np.array([]))
                if right is not None else np.array([])
            )

            if len(left_val) == 0:
                combined_corr[f'{base}-both'][correlation_type] = right_val

            elif len(right_val) == 0:
                combined_corr[f'{base}-both'][correlation_type] = left_val

            else:
                combined_corr[f'{base}-both'][correlation_type] = np.concatenate(
                    [left_val, right_val],
                    axis=0
                )


    # ------------------------------------------------------------
    # Plot combined hemispheres
    # ------------------------------------------------------------

    sites = combined_corr.keys()
   
    n_sites = len(sites)
    ncols = 3
    if focused_analysis:
        nrows = 1
        sites = ['SC-both', 'ZI-both', 'PAG-both']
    else:
        nrows = math.ceil(n_sites / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5 * ncols, 4 * nrows),
        sharex=True,
        constrained_layout=True
    )

    axes = np.atleast_1d(axes).ravel()

    for ax, site in zip(axes, sites):

        for correlation_type in correlation_list:
            
            if correlation_type == 'approach':
                plt_color = [0.47, 0.67, 0.19]  # green
               
                    
            elif correlation_type == 'IR':
                plt_color = [0.416, 0.741, 0.741] #blue
                
            else:
                plt_color = [0.78, 0, 0] # red
            
            # if site not in combined_corr.keys():
            #     continue

            data = combined_corr[site][correlation_type]

            # if correlation_type not in data:
            #     continue
            
            # correlations = data[correlation_type]

            # Restrict to analysis window
            data = data[
                :,
                int(len(data[0])/2)-pre*sr:int(len(data[0])/2)+pre*sr
            ]
    

            # Handle missing data
            if not isinstance(data, np.ndarray) or data.size == 0:
                continue

            # Remove trials containing NaNs
            data = data[
                ~np.isnan(data).any(axis=1)
            ]

            if len(data) == 0:
                continue

            # Correlation array is assumed to be:
            # trials x lags

            mean_corr = np.nanmean(data, axis=0)
            sem_corr = (
                np.nanstd(data, axis=0)
                / np.sqrt(len(data))
            )

            # Create lag vector
            lags = np.arange(
                -data.shape[1] // 2 + 1,
                data.shape[1] // 2 + 1
            ) / sr

            ax.plot(
                lags,
                mean_corr,
                color = plt_color,
                label=f'{correlation_type} (n={len(data)})'
            )

            ax.fill_between(
                lags,
                mean_corr - sem_corr,
                mean_corr + sem_corr,
                color = plt_color,
                alpha=0.2
            )

        ax.axvline(
            0,
            linestyle='--',
            color='black',
            linewidth=1
        )

        ax.axhline(
            0,
            linestyle='--',
            color='black',
            linewidth=1
        )

        ax.set_title(site)
        ax.set_ylabel('Cross-correlation')
        ax.set_xlabel('Lag (s)')
        ax.legend()
        boxoff()

    for ax in axes[len(sites):]:
        ax.set_visible(False)
        

    plt.show()

