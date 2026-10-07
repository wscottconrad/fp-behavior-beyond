# -*- coding: utf-8 -*-
"""
Created on 17 Aug 2026

@author: conrad
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.stats import sem
# from scipy.spatial.distance import mahalanobis
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import pandas as pd
# from sklearn.model_selection import cross_val_score
# import pickle
import random


def angle_between(v1, v2):
    
    # Calculate the dot product
    dot_prod = np.dot(v1, v2)
    
    # Calculate L2 norms (magnitudes)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    
    if norm_v1 == 0 or norm_v2 == 0:
        return np.nan
    
    # Clip cosine value to avoid numerical precision errors outside [-1.0, 1.0]
    cosine_angle = np.clip(dot_prod / (norm_v1 * norm_v2), -1.0, 1.0)
    
    # Calculate angles in radians and degrees
    angle_radians = np.arccos(cosine_angle)
    angle_degrees = np.degrees(angle_radians)
    
    return angle_degrees

def shuffle_condition_labels(f_trials, nf_trials, rng):
    """
    Pool freeze and non-freeze trial responses and randomly reassign
    freeze/non-freeze condition labels.

    f_trials and nf_trials:
        list of arrays, one array per neuron.
        Each array has shape (n_trials,).

    The same trial permutation is applied to every neuron so that
    population-vector trial identity is preserved.
    """

    n_freeze = f_trials[0].shape[0]
    n_nonfreeze = nf_trials[0].shape[0]

    # Pool trials for every neuron
    pooled_trials = [
        np.concatenate([f, nf])
        for f, nf in zip(f_trials, nf_trials)
    ]

    # Original condition labels
    labels = np.array(
        [0] * n_freeze +
        [1] * n_nonfreeze
    )

    # Shuffle the condition labels
    shuffled_labels = rng.permutation(labels)

    freeze_idx = np.where(shuffled_labels == 0)[0]
    nonfreeze_idx = np.where(shuffled_labels == 1)[0]

    # Reconstruct shuffled conditions
    shuffled_f = [
        neuron_trials[freeze_idx]
        for neuron_trials in pooled_trials
    ]

    shuffled_nf = [
        neuron_trials[nonfreeze_idx]
        for neuron_trials in pooled_trials
    ]

    return shuffled_f, shuffled_nf

# -----------------------------
# Parameters
# -----------------------------

# user input

plot_trial_heat_maps = True

plot_stim = True
plot_resp = True

baseline_correct = True 
sqrt_transform = True # default true


# summary_data = 'Baseline_corrected'
summary_data = 'Average'
# summary_data = 'Z-score'

shuffle_trials = True
n_shuffles = 10000
random_seed = 42
rng = np.random.default_rng(random_seed)





# parameter setup
general_baseline = 50
analysis_baseline = 25 # 50 = 1 seconds, try 25 to match MD analysis
end = 225 # stimulus end
all_trial_SD = True 
labelsX = np.arange(-1, 4, 1)

alpha = 0.05 # threshold for modulated neurons from zeta test, should be 0.05
sigma = 1.0  # smoothing, in bins

var_thresh = 10**-10

spike_path = r'W:\Haak\Innate_defense\Data_analysis\22.35.02\Scotts_analysis\collected_spikes.pickle'
collected_spikes = pd.read_pickle(spike_path)


# shouldn't hard code this
neuron_idx_animal= {'superior colliculus' : [[0,138],[139,294],[295,359],[360,493]],
                    'periaqueductal gray': [[0,89],[90,164],[165,207],[208,239]]}

# #test
# random.seed(1)
# test1 = np.random.rand(150,250)
# random.seed(21)
# test2 = np.random.rand(150,250)
# mean_test = angle_between(np.mean(test1,axis = 1), np.mean(test2,axis = 1))
# ts_test = (np.array([
#     angle_between(test1[:,t], test2[:,t])
#     for t in range(250)
# ]))


trial_matrices = {}

for region in collected_spikes:
    FvL_freeze = []
    FvL_nonfreeze = []
    FvF = []
    LvF = []
    PUvL_nonfreeze = [] # second to last (PenUltimate) vs last nonfreeze, hoping this is similar
    
    ts_FvL_freeze = []
    ts_FvL_nonfreeze = []
    ts_FvF = []
    ts_LvF = []
    
    f_trials = collected_spikes[region]['all_neuron_trials_freeze']
    nf_trials = collected_spikes[region]['all_neuron_trials_nonfreeze']
    
    if sqrt_transform:
        f_trials = [np.sqrt(x) for x in f_trials]
        nf_trials = [np.sqrt(x) for x in nf_trials]
        
    if baseline_correct:
        f_trials = [x-np.mean(x[:,0:general_baseline+1],axis=1).reshape(-1,1) for x in f_trials]
        nf_trials = [x-np.mean(x[:,0:general_baseline+1],axis=1).reshape(-1,1) for x in nf_trials]
        
    # takes average of post stimulus response
    f_avg = [np.mean(x[:,general_baseline-1:end+1], axis = 1) for x in f_trials]    
    nf_avg = [np.mean(x[:,general_baseline-1:end+1], axis = 1) for x in nf_trials]
    
    i = 1
    for start,stop in neuron_idx_animal[region]:
        
        # averaged PVs (one value per neuron, per trial)
        pv_freeze_array = np.array(f_avg[start:stop+1])        
        pv_nonfreeze_array = np.array(nf_avg[start:stop+1])

        FvL_freeze.append(angle_between(pv_freeze_array[:,0],pv_freeze_array[:,-1]))
        FvL_nonfreeze.append(angle_between(pv_nonfreeze_array[:,0],pv_nonfreeze_array[:,-1]))
        
        FvF.append(angle_between(pv_freeze_array[:,0],pv_nonfreeze_array[:,0]))
        LvF.append(angle_between(pv_freeze_array[:,-1],pv_nonfreeze_array[:,0]))
        
        PUvL_nonfreeze.append(angle_between(pv_nonfreeze_array[:,-2],pv_nonfreeze_array[:,-1]))
        
        # ---------------------------------------------------------
        # OBSERVED DATA
        # ---------------------------------------------------------
        
        FvF_observed = angle_between(
            pv_freeze_array[:, 0],
            pv_nonfreeze_array[:, 0]
        )
        
        LvF_observed = angle_between(
            pv_freeze_array[:, -1],
            pv_nonfreeze_array[:, 0]
        )
        
        # Positive value = freeze PV became closer to non-freeze
        observed_delta = FvF_observed - LvF_observed
        
        # ---------------------------------------------------------
        # PERMUTATION NULL DISTRIBUTION
        # ---------------------------------------------------------
        
        shuffle_delta = []

        for shuffle in range(n_shuffles):
        
            shuffled_f, shuffled_nf = shuffle_condition_labels(
                f_avg[start:stop+1],
                nf_avg[start:stop+1],
                rng
            )
        
            shuffled_f = np.array(shuffled_f)
            shuffled_nf = np.array(shuffled_nf)
        
            # First shuffled freeze vs first shuffled non-freeze
            shuffled_FvF = angle_between(
                shuffled_f[:, 0],
                shuffled_nf[:, 0]
            )
        
            # Last shuffled freeze vs first shuffled non-freeze
            shuffled_LvF = angle_between(
                shuffled_f[:, -1],
                shuffled_nf[:, 0]
            )
        
            shuffle_delta.append(
                shuffled_FvF - shuffled_LvF
            )
        
        shuffle_delta = np.array(shuffle_delta)
        
        p_value = (
            np.sum(shuffle_delta >= observed_delta) + 1
        ) / (n_shuffles + 1)
        
        print(f'P-value for animal {i} {region}: {p_value}')

        plt.figure(figsize=(8, 5))

        plt.hist(
            shuffle_delta,
            bins=50,
            alpha=0.7
        )
        
        plt.axvline(
            observed_delta,
            linewidth=2,
            label=f'Observed Δ = {observed_delta:.2f}°'
        )
        
        plt.axvline(
            0,
            linestyle='--',
            linewidth=1,
            label='No habituation'
        )
        
        plt.xlabel('FvF − LvF angle difference (degrees)')
        plt.ylabel('Permutation count')
        plt.title(f'{region} - Animal {i}: permutation test')
        plt.legend()
        
        plt.show()
        
        i +=1
        
        # # ---------------------------------------------------------
        # # Shuffle condition labels
        # # ---------------------------------------------------------
        
        # shuffle_FvL_freeze = []
        # shuffle_FvL_nonfreeze = []
        # shuffle_FvF = []
        # shuffle_LvF = []
        
        # for shuffle in range(n_shuffles):
        
        #     shuffled_f_trials, shuffled_nf_trials = (
        #         shuffle_condition_labels(
        #             f_avg[start:stop+1],
        #             nf_avg[start:stop+1],
        #             rng
        #         )
        #     )
        
        #     shuffled_f = np.array(shuffled_f_trials)
        #     shuffled_nf = np.array(shuffled_nf_trials)
        
        #     shuffle_FvL_freeze.append(
        #         angle_between(
        #             shuffled_f[:, 0],
        #             shuffled_f[:, -1]
        #         )
        #     )
        
        #     shuffle_FvL_nonfreeze.append(
        #         angle_between(
        #             shuffled_nf[:, 0],
        #             shuffled_nf[:, -1]
        #         )
        #     )
        
        #     shuffle_FvF.append(
        #         angle_between(
        #             shuffled_f[:, 0],
        #             shuffled_nf[:, 0]
        #         )
        #     )
        
        #     shuffle_LvF.append(
        #         angle_between(
        #             shuffled_f[:, -1],
        #             shuffled_nf[:, 0]
        #         )
        #     )
            
        # time-series PVs
        ts_pv_freeze = f_trials[start:stop+1]
        num_freeze_trials = ts_pv_freeze[0].shape[0]
        
        ts_pv_nonfreeze = nf_trials[start:stop+1]
        num_nonfreeze_trials = ts_pv_nonfreeze[0].shape[0]
        
        # First and last trial for every neuron
        freeze_first = np.array([
            neuron_trials[0, :]
            for neuron_trials in ts_pv_freeze
        ])
    
        freeze_last = np.array([
            neuron_trials[-1, :]
            for neuron_trials in ts_pv_freeze
        ])
    
        nonfreeze_first = np.array([
            neuron_trials[0, :]
            for neuron_trials in ts_pv_nonfreeze
        ])
    
        nonfreeze_last = np.array([
            neuron_trials[-1, :]
            for neuron_trials in ts_pv_nonfreeze
        ])
        
        # calc angle
        ts_FvL_freeze.append(np.array([
            angle_between(freeze_first[:, t], freeze_last[:, t])
            for t in range(freeze_first.shape[1])
        ]))
        
        ts_FvL_nonfreeze.append(np.array([
            angle_between(nonfreeze_first[:, t], nonfreeze_last[:, t])
            for t in range(nonfreeze_first.shape[1])
        ]))
        
        ts_FvF.append(np.array([
            angle_between(freeze_first[:, t], nonfreeze_first[:, t])
            for t in range(freeze_first.shape[1])
        ]))
        
        ts_LvF.append(np.array([
            angle_between(freeze_last[:, t], nonfreeze_first[:, t])
            for t in range(freeze_last.shape[1])
        ]))
        
    # averaged
    plt.figure()

    # Plot individual animal trajectories
    for i in range(len(FvF)):
        values = [
            PUvL_nonfreeze[i],
            FvF[i],
            LvF[i],
            FvL_freeze[i],
            FvL_nonfreeze[i]
        ]
        
        plt.plot(
            [1, 2, 3, 4, 5],
            values,
            marker='o',
            alpha=0.5
        )
    
    plt.xticks(
        [1, 2, 3, 4, 5],
        ['PUvL non-freeze', 'FvF', 'LvF', 'FvL freeze', 'FvL non-freeze']
    )
    
    plt.ylim(0, 180)
    plt.xlim(0.5, 5.5)
    plt.ylabel('Angle')
    plt.title(f'{region}: Angle differences')
    
    plt.show()
    
    
    # time series
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)

    # Put your four datasets into the corresponding axes
    datasets = [
        (ts_FvF, f'{region}: FvF'),
        (ts_LvF, f'{region}: LvF'),
        (ts_FvL_freeze, f'{region}: FvL freeze'),
        (ts_FvL_nonfreeze, f'{region}: FvL non-freeze')
    ]
    
    for ax, (data, title) in zip(axes.flat, datasets):
    
        # Plot each animal
        for animal_ts in data:
            ax.plot(animal_ts, alpha=0.3)
    
        # Mean across animals
        # data_array = np.array(data)
        # mean_ts = np.nanmean(data_array, axis=0)

        # sem_ts = (
        #     np.nanstd(data_array, axis=0, ddof=1)
        #     / np.sqrt(np.sum(~np.isnan(data_array), axis=0))
        # )
    
        # ax.plot(mean_ts, linewidth=2)
        # ax.fill_between(
        #     np.arange(len(mean_ts)),
        #     mean_ts - sem_ts,
        #     mean_ts + sem_ts,
        #     alpha=0.2
        # )
    
        ax.set_title(title)
        ax.set_ylabel('Angle')
        ax.set_ylim(0, 180)
    
    axes[1, 0].set_xlabel('Time')
    axes[1, 1].set_xlabel('Time')
    
    plt.tight_layout()
    plt.show()
    
    
  
    

