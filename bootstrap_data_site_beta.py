# -*- coding: utf-8 -*-
"""
Plot and statistically compare approach responses between
NT and FM experiments.

For each site:
    - NT approach traces are plotted
    - FM approach traces are plotted
    - Mean +/- SEM is shown
    - Bootstrap significance is calculated separately for NT and FM
    - Permutation test compares NT vs FM
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from perm_test_array import perm_test_array
from bootstrap_data_bias_corrected import bootstrap_data
from consec_idx import consec_idx
from scipy.stats import ttest_rel
from scipy.stats import shapiro
import scipy.stats as stats
from boxoff import boxoff
import math


# ============================================================
# SETTINGS
# ============================================================

initiate_aligned = True

perm_testing = True

combine_hemispheres = True
plot_site = True

focused_analysis = True
whole_trial = False

plot_signal = True
plot_speed = False
plot_correlations = False

# nt = True # will be phased out

correlation_type = 'approach'  # 'approach', 'ITI', or 'IR'

calculate_auc = False

sr = 30

trial_type = 'approach'

if trial_type == 'ITI':
    initiate_aligned = True


pre = 5

if whole_trial:
    post = 25
else:
    post = 5

thres = round(1/3*sr)  # Consecutive threshold length


# ============================================================
# DATA LOCATIONS
# ============================================================

tankfolders = {
    'nt': r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\',

    'fm': r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\\freelymoving\\'
}


# ============================================================
# TIME VECTOR
# ============================================================

ts = np.linspace(
    -pre,
    post,
    (pre + post) * sr
)


# ============================================================
# LOAD BOTH DATASETS
# ============================================================

datasets = {}

for exp_type in ['nt', 'fm']:

    print('\n====================================')
    print(f'Loading {exp_type.upper()} data')
    print('====================================')

    tankfolder = tankfolders[exp_type]

    # Load combined data
    with open(f'{tankfolder}allDatComb.pkl', 'rb') as f:
        d = pickle.load(f)


    # --------------------------------------------------------
    # Select alignment
    # --------------------------------------------------------
    
    if trial_type == 'approach':
        if initiate_aligned:
            trialData = d['trialData']
        else:
            trialData = d['trialData_trialOnset']
            
    elif trial_type == 'ITI':
        trialData = d['ITI']


    # --------------------------------------------------------
    # Combine hemispheres
    # --------------------------------------------------------

    if combine_hemispheres:

        combined_data = {}

        for key in trialData.keys():

            # Only process left hemisphere entries
            if key.endswith('-L'):

                base = key[:-2]

                left = trialData.get(f'{base}-L')
                right = trialData.get(f'{base}-R')


                if left is not None and right is not None:

                    combined_data[f'{base}-both'] = {}

                    for k in left.keys():

                        left_val = left[k]
                        right_val = right.get(k)


                        if right_val is None:

                            combined_data[f'{base}-both'][k] = left_val

                        elif left_val.size == 0:

                            combined_data[f'{base}-both'][k] = right_val

                        elif right_val.size == 0:

                            combined_data[f'{base}-both'][k] = left_val

                        else:

                            combined_data[f'{base}-both'][k] = np.concatenate(
                                [left_val, right_val]
                            )


                else:

                    # Handle missing side
                    combined_data[f'{base}-both'] = (
                        left if left is not None else right
                    )


        trialData = combined_data


    datasets[exp_type] = trialData


# ============================================================
# FIND SITES PRESENT IN BOTH DATASETS
# ============================================================

nt_sites = datasets['nt'].keys()
fm_sites = datasets['fm'].keys()

sites = [site for site in nt_sites if site in fm_sites]


print('\nSites present in both datasets:')
print(sites)


# ============================================================
# FIGURE SETUP
# ============================================================
# %%
if plot_signal:
    n_sites = len(sites)
    ncols = 3
    
    if focused_analysis:
        nrows = 1
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
    
    
    # ============================================================
    # MAIN LOOP
    # ============================================================
    
    ax_idx = 0
    
    
    for site in sites:
    
    
        # --------------------------------------------------------
        # Skip projections
        # --------------------------------------------------------
    
        if focused_analysis:
    
            if 'to' in site:
                continue
    
            if site == 'MLR-both':
                continue
    
    
        ax = axes[ax_idx]
        ax_idx += 1
    
    
        # ========================================================
        # EXTRACT APPROACH DATA
        # ========================================================
    
        plot_data = {}
    
    
        for exp_type in ['nt', 'fm']:
    
            data = datasets[exp_type][site]
    
    
            if trial_type not in data:
    
                print(
                    f'{site}: no {trial_type} data '
                    f'for {exp_type.upper()}'
                )
    
                plot_data[exp_type] = np.empty(
                    (0, (pre + post) * sr)
                )
    
                continue
    
    
            signal = data[trial_type]
    
    
            # Remove trials containing NaNs
            signal = signal[
                ~np.isnan(signal).any(axis=1)
            ]
    
    
            # Restrict to analysis window
            signal = signal[
                :,
                :(pre + post) * sr
            ]
    
    
            plot_data[exp_type] = signal
    
    
            print(
                f'{site} - {exp_type.upper()} '
                f'approach trials: {len(signal)}'
            )
    
    
        # ========================================================
        # DETERMINE Y LIMIT FOR SIGNIFICANCE MARKERS
        # ========================================================
    
        mean_values = []
    
        for exp_type in ['nt', 'fm']:
    
            signal = plot_data[exp_type]
    
            if len(signal) > 0:
    
                mean_values.append(
                    np.mean(signal, axis=0)
                )
    
    
        if len(mean_values) > 0:
    
            ymax_total = np.max(mean_values) + 1
    
        else:
    
            ymax_total = 1
    
    
        # ========================================================
        # BOOTSTRAPPING
        # ========================================================
    
        bootstrap_results = {}
    
        for exp_type in ['nt', 'fm']:
    
            signal = plot_data[exp_type]
            
            if trial_type == 'approach':
                if exp_type == 'nt':
                    plt_color = [0.47, 0.67, 0.19]  # green
                else:
                    plt_color = [1, 0.804, 0.004] # orange
                    
            else:
                if exp_type == 'nt':
                   plt_color = [0.416, 0.741, 0.741] #blue
                else:
                   plt_color = [0.78, 0, 0] # red
    
    
    
            if len(signal) == 0:
                continue
    
    
            print(
                f'Bootstrapping {site} '
                f'{exp_type.upper()}...'
            )
    
    
            btsrp_app = bootstrap_data(
                signal,
                10000,
                0.0001
            )
    
    
            bootstrap_results[exp_type] = btsrp_app
    
    
            # ----------------------------------------------------
            # Plot mean
            # ----------------------------------------------------
    
            mean_signal = np.nanmean(
                signal,
                axis=0
            )
    
    
            # ----------------------------------------------------
            # SEM
            # ----------------------------------------------------
    
            sem_signal = (
                np.nanstd(signal, axis=0)
                / np.sqrt(len(signal))
            )
    
    
    
    
            linestyle = '-'
    
      
            # ----------------------------------------------------
            # Plot mean response
            # ----------------------------------------------------
    
            ax.plot(
                ts,
                mean_signal,
                color=plt_color,
                linestyle=linestyle,
                linewidth=2,
                label=(
                    f'{exp_type.upper()} {trial_type} '
                    f'(n={len(signal)})'
                )
            )
    
    
            # ----------------------------------------------------
            # Plot SEM
            # ----------------------------------------------------
    
            ax.fill_between(
                ts,
                mean_signal + sem_signal,
                mean_signal - sem_signal,
                color=plt_color,
                alpha=0.2
            )
    
    
            # ====================================================
            # BOOTSTRAP SIGNIFICANCE
            # ====================================================
    
            # Negative significant periods
            tmp = np.where(
                btsrp_app[1, :] < 0
            )[0]
    
    
            if len(tmp) > 1:
    
                id = tmp[
                    consec_idx(tmp, thres)
                ]
    
    
                # Put NT and FM bootstrap markers
                # at slightly different heights
    
                if exp_type == 'nt':
                    marker_y = ymax_total + 1.0
                else:
                    marker_y = ymax_total + 0.5
    
    
                ax.plot(
                    ts[id],
                    marker_y * np.ones(len(id)),
                    's',
                    markersize=7,
                    markerfacecolor=plt_color,
                    color=plt_color
                )
    
    
            # Positive significant periods
            tmp = np.where(
                btsrp_app[0, :] > 0
            )[0]
    
    
            if len(tmp) > 1:
    
                id = tmp[
                    consec_idx(tmp, thres)
                ]
    
    
                if exp_type == 'nt':
                    marker_y = ymax_total + 1.0
                else:
                    marker_y = ymax_total + 0.5
    
    
                ax.plot(
                    ts[id],
                    marker_y * np.ones(len(id)),
                    's',
                    markersize=7,
                    markerfacecolor=plt_color,
                    color=plt_color
                )
    
    
        # ========================================================
        # PERMUTATION TEST
        # NT APPROACH VS FM APPROACH
        # ========================================================
    
        if perm_testing:
    
            nt_signal = plot_data['nt']
            fm_signal = plot_data['fm']
    
    
            # Only perform test if both datasets contain trials
    
            if (
                len(nt_signal) > 0
                and len(fm_signal) > 0
            ):
    
                print(
                    f'Permuting NT vs FM '
                    f'for {site}...'
                )
    
    
                perm_test, _ = perm_test_array(
                    nt_signal,
                    fm_signal,
                    1000
                )
    
    
                # ------------------------------------------------
                # Find significant time points
                # ------------------------------------------------
    
                perm_hits = np.where(
                    perm_test < 0.05
                )[0]
    
    
                if len(perm_hits) > 0:
    
                    perm_x = perm_hits[
                        consec_idx(
                            perm_hits,
                            thres
                        )
                    ]
    
    
                    # ------------------------------------------------
                    # Purple markers = NT vs FM permutation
                    # ------------------------------------------------
    
                    perm_y = ymax_total + 1.7
    
    
                    ax.plot(
                        ts[perm_x],
                        perm_y * np.ones(len(perm_x)),
                        's',
                        markersize=7,
                        markerfacecolor=[0.659, 0.427, 0.91],
                        color=[0.659, 0.427, 0.91]
                    )
    
    
        # ========================================================
        # AXIS FORMATTING
        # ========================================================
    
        ax.axvline(
            x=0,
            linestyle='--',
            color='black',
            linewidth=1.5
        )
    
    
        ax.axhline(
            y=0,
            linestyle='--',
            color='black',
            linewidth=1.5
        )
    
    
        # Hide top/right spines
    
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    
        ax.tick_params(
            axis='x',
            which='both',
            direction='out',
            bottom=True,
            top=False
        )
    
    
        ax.tick_params(
            axis='y',
            which='both',
            direction='out',
            left=True,
            right=False
        )
    
    
        # ========================================================
        # LABELS
        # ========================================================
    
        ax.set_title(site)
    
        ax.set_ylabel('Z-Score')
    
    
        if initiate_aligned:
    
            ax.set_xlabel(
                'Movement onset (s)'
            )
    
        else:
    
            ax.set_xlabel(
                'Prey Laser onset (s)'
            )
    
    
        ax.legend()
    
    
    # ============================================================
    # HIDE UNUSED AXES
    # ============================================================
    
    for ax in axes[ax_idx:]:
    
        ax.set_visible(False)
    
    
    # Make sure visible axes have x labels
    
    for ax in axes.flat:
    
        if ax.get_visible():
    
            if initiate_aligned:
    
                ax.set_xlabel(
                    'Movement onset (s)'
                )
    
            else:
    
                ax.set_xlabel(
                    'Prey Laser onset (s)'
                )
    
            ax.tick_params(
                labelbottom=True
            )
    
    
    plt.show()


    
  
# %%
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
            
# %%
            
# ============================================================
# PLOTTING SPEED BY COMBINED HEMISPHERE
# ============================================================

if plot_speed:

    # Use speedData from both experiments
    speed_datasets = {}

    for exp_type in ['nt', 'fm']:

        tankfolder = tankfolders[exp_type]

        with open(f'{tankfolder}allDatComb.pkl', 'rb') as f:
            d = pickle.load(f)

        speedData = d['speedData']

        # --------------------------------------------------------
        # Combine L/R speed data into -both
        # --------------------------------------------------------

        combined_speed = {}

        for key in speedData.keys():

            if not key.endswith('-L'):
                continue

            base = key[:-2]

            left = speedData.get(f'{base}-L')
            right = speedData.get(f'{base}-R')

            if left is None and right is None:
                continue

            combined_speed[f'{base}-both'] = {}

            # Speed variables that can be combined
            for speed_key in [
                'approach_snout_speed',
                'IR_snout_speed',
                'ITIspeed',
                'speedTrialsMov'
            ]:

                left_val = (
                    left.get(speed_key, np.array([]))
                    if left is not None else np.array([])
                )

                right_val = (
                    right.get(speed_key, np.array([]))
                    if right is not None else np.array([])
                )

                if len(left_val) == 0:
                    combined_speed[f'{base}-both'][speed_key] = right_val

                elif len(right_val) == 0:
                    combined_speed[f'{base}-both'][speed_key] = left_val

                else:
                    combined_speed[f'{base}-both'][speed_key] = np.concatenate(
                        [left_val, right_val],
                        axis=0
                    )

        speed_datasets[exp_type] = combined_speed

    # ------------------------------------------------------------
    # Plot combined hemispheres
    # ------------------------------------------------------------

    # sites = sorted(
    #     set(speed_datasets['nt'].keys()) |
    #     set(speed_datasets['fm'].keys())
    # )
    nt_sites = speed_datasets['nt'].keys()
    fm_sites = speed_datasets['fm'].keys()

    sites = [site for site in nt_sites if site in fm_sites]


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

        for exp_type in ['nt', 'fm']:

            if site not in speed_datasets[exp_type]:
                continue
            
            if focused_analysis:
        
                if 'to' in site:
                    continue
        
                if site == 'MLR-both':
                    continue

            speed = speed_datasets[exp_type][site]

            if exp_type == 'fm':
                approach_speed = speed['approach_snout_speed']
                baseline_speed = speed['ITIspeed']
            else:
                approach_speed = speed['speedTrialsMov'] / 1000
                baseline_speed = speed['ITIspeed'] / 1000

            if len(approach_speed) == 0:
                continue

            approach_speed = approach_speed[:, :(pre + post) * sr]

            # Plot approach speed
            ax.plot(
                ts,
                np.nanmean(approach_speed, axis=0),
                label=f'{exp_type.upper()} approach'
            )

            ax.fill_between(
                ts,
                np.nanmean(approach_speed, axis=0)
                + np.nanstd(approach_speed, axis=0)
                / np.sqrt(len(approach_speed)),
                np.nanmean(approach_speed, axis=0)
                - np.nanstd(approach_speed, axis=0)
                / np.sqrt(len(approach_speed)),
                alpha=0.2
            )

        ax.axvline(
            0,
            linestyle='--',
            color='black',
            linewidth=1
        )

        ax.set_title(site)
        ax.set_ylabel('Speed (m/s)')
        ax.set_xlabel('Time, approach onset (s)')
        ax.legend()
        boxoff()

    for ax in axes[len(sites):]:
        ax.set_visible(False)

    plt.show()
    
# %%
# ============================================================
# PLOTTING CORRELATIONS BY SITE
# ============================================================

if plot_correlations:

    corr_key = f'lag_correlations_{correlation_type}'
    
    corr_datasets = {}

    for exp_type in ['nt', 'fm']:

        tankfolder = tankfolders[exp_type]

        with open(f'{tankfolder}allDatComb.pkl', 'rb') as f:
            d = pickle.load(f)

        corrData = d['lag correlations']

        # --------------------------------------------------------
        # Combine L/R speed data into -both
        # --------------------------------------------------------

        combined_corr = {}

        for key in corrData.keys():

            if not key.endswith('-L'):
                continue

            base = key[:-2]

            left = corrData.get(f'{base}-L')
            right = corrData.get(f'{base}-R')

            if left is None and right is None:
                continue

            combined_corr[f'{base}-both'] = {}

            # Speed variables that can be combined
            # for corr_key in [
            #     'approach',
            #     'IR',
            #     'ITI'
            # ]:

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

        corr_datasets[exp_type] = combined_corr

    # ------------------------------------------------------------
    # Plot combined hemispheres
    # ------------------------------------------------------------

    nt_sites = corr_datasets['nt'].keys()
    fm_sites = corr_datasets['fm'].keys()

    sites = [site for site in nt_sites if site in fm_sites]
    
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

        for exp_type in ['nt', 'fm']:
            
            if correlation_type == 'approach':
                if exp_type == 'nt':
                    plt_color = [0.47, 0.67, 0.19]  # green
                else:
                    plt_color = [1, 0.804, 0.004] # orange
                    
            else:
                if exp_type == 'nt':
                   plt_color = [0.416, 0.741, 0.741] #blue
                else:
                   plt_color = [0.78, 0, 0] # red
            
            if site not in corr_datasets[exp_type]:
                continue

            data = corr_datasets[exp_type][site]

            if correlation_type not in data:
                continue
            
            correlations = data[correlation_type]

            # Restrict to analysis window
            correlations = correlations[
                :,
                int(len(correlations[0])/2)-pre*sr:int(len(correlations[0])/2)+pre*sr
            ]
    

            # Handle missing data
            if not isinstance(correlations, np.ndarray) or correlations.size == 0:
                continue

            # Remove trials containing NaNs
            correlations = correlations[
                ~np.isnan(correlations).any(axis=1)
            ]

            if len(correlations) == 0:
                continue

            # Correlation array is assumed to be:
            # trials x lags

            mean_corr = np.nanmean(correlations, axis=0)
            sem_corr = (
                np.nanstd(correlations, axis=0)
                / np.sqrt(len(correlations))
            )

            # Create lag vector
            lags = np.arange(
                -correlations.shape[1] // 2 + 1,
                correlations.shape[1] // 2 + 1
            ) / sr

            ax.plot(
                lags,
                mean_corr,
                color = plt_color,
                label=f'{exp_type.upper()} {correlation_type} (n={len(correlations)})'
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
# %%