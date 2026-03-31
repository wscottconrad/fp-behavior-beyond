# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:29:51 2026

Calculates neural trajectories and mahalonobis distances between respones types.

@author: sconrad
"""

import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from sklearn.decomposition import PCA
from scipy.spatial.distance import mahalanobis
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
# from matplotlib.collections import LineCollection
from kneed import KneeLocator
import pandas as pd
import csv
from find_field import find_field

# -----------------------------
# Parameters
# -----------------------------

# user input
plot_speed_heatmap = False
plot_speed_trials = False 
plot_trajectories = True

recalc_freeze_times = False

region = 'periaqueductal gray' # or
# region = 'superior colliculus' # or
# region = 'retrosplenial'   # or 
# region = 'midbrain'   


# general
np.random.seed(10)

# spike analysis initialize
t_pre = 1.0 # s
t_post = 4.0 # s
bin_size = 0.02  # 20 ms
smooth_sigma = 1.0  # in bins

n_pop_trials = 30
n_repeats = 50  # increase to 200 to match Jercog 2021
latent_dims = 5 

time_bins = np.arange(-t_pre, t_post + bin_size, bin_size)
time_centers = time_bins[:-1] + bin_size / 2

# speed analysis initialize
# params are probably buried in sAP but i dont know where yet
fs = 1000  # Hz
s_pre = 4    # seconds before event onset
s_post = 6  # seconds after event onset

n_pre = int(s_pre * fs)
n_post = int(s_post * fs)
win_len = n_pre + n_post

stim_length = 3.5 # seconds. rough estimate, i think stim duration may vary slightly trial to trial?

freeze_thresh = 0.015  # m/s (adjust as appropriate)
min_freeze_samples = int(0.5 * fs) #500 ms
exclude_samples = int(0.25 * fs) # exclude freezing onset before 250ms after stim onset
event_end_sample = int(stim_length * fs)
pre_freeze_window = 1.0 # seconds to look back before freeze
pre_freeze_movement_thresh = 0.1  # cm/s  0.07


data_path = r'W:\Haak\Innate_defense\Data_analysis\22.35.02'
session_list = [
    '98332_20240326_AP.mat',
    '98335_20240321_AP.mat',
    '100131_20240508_AP.mat',
    '100132_20240508_AP.mat',
    '100134_20240514_AP.mat'
]

freeze_file_path = r'W:\Conrad\Innate_approach\Materials_and_methods\Code'

# -----------------------------
# Helper functions
# -----------------------------
def binned_firing_rate(spike_times, event_time, time_bins, bin_size):
    aligned = spike_times - event_time
    counts, _ = np.histogram(aligned, bins=time_bins) # bins spikes within trial window
    counts = np.sqrt(counts)  # variance-stabilizing transform # anscombe transform possible, but i havent seen it in papers
    return counts / bin_size


def neuron_trial_matrix(spike_times, event_onset): # for a single neuron, creates trial-based binned spike matrix
    mat = np.zeros((len(event_onset), len(time_bins) - 1))
    for i, ev in enumerate(event_onset):
        rate = binned_firing_rate(spike_times, ev, time_bins, bin_size)
        rate = gaussian_filter1d(rate, smooth_sigma)
        mat[i] = rate
    return mat


def sample_population_trial(neuron_trials):
    pop = np.zeros((len(neuron_trials), neuron_trials[0].shape[1]))
    for i, trials in enumerate(neuron_trials):
        idx = np.random.randint(trials.shape[0])
        pop[i] = trials[idx]
    return pop


def generate_population_trials(all_neuron_trials, n_pop_trials, n_repeats):
    all_population_trials = []
    for _ in range(n_repeats):
        pop_trials = np.array([
            sample_population_trial(all_neuron_trials)
            for _ in range(n_pop_trials)
        ])
        all_population_trials.append(pop_trials)
    return np.concatenate(all_population_trials, axis=0)

def generate_population_trials_by_repeat(all_neuron_trials, n_pop_trials, n_repeats):
    """
    Returns:
        pop_trials: shape (n_repeats, n_pop_trials, n_neurons, n_time)
    """
    pop_trials = []

    for _ in range(n_repeats):
        trials = np.array([
            sample_population_trial(all_neuron_trials)
            for _ in range(n_pop_trials)
        ])
        pop_trials.append(trials)

    return np.array(pop_trials)

def flatten_for_pca(pop): # literally just flattens list of lists into array
    if len(pop.shape) == 3:
        return pop.transpose(0, 2, 1).reshape(-1, pop.shape[1])
    elif len(pop.shape) == 4:
        n_neurons = pop.shape[2] # not needed, but handy to know its reshaped to n_neurons (for above case too)
        return pop.transpose(0, 1, 3, 2).reshape(-1, n_neurons)


# def get_latent(pop, pca, latent_dims):
#     X = flatten_for_pca(pop)
#     X_latent = pca.transform(X) # applies dim reduction (projected onto PC's from traning set)
#     return X_latent.reshape(pop.shape[0], pop.shape[2], latent_dims)


def get_latent(pop, pca_space, latent_dims):

    X = flatten_for_pca(pop)
    X_latent = pca_space.transform(X) # applies dim reduction (projected onto PC's from traning set)

    if len(pop.shape) == 3:
        # (n_trials, n_neurons, n_time)
        n_trials, n_neurons, n_time = pop.shape
        return X_latent.reshape(n_trials, n_time, latent_dims)

    elif len(pop.shape) == 4:
        # (n_repeats, n_trials, n_neurons, n_time)
        n_repeats, n_trials, n_neurons, n_time = pop.shape
        return X_latent.reshape(n_repeats, n_trials, n_time, latent_dims)

def mahalanobis_distance_trials(latent_A, latent_B):
    """
    latent_A, latent_B:
        shape = (n_trials, n_timepoints, latent_dims)
    Returns:
        dist[t] = Mahalanobis distance between condition means at time t
        
        do i want distance between trial means or per trial pairs, then take mean?
    """

    n_time = latent_A.shape[1]
    dist = np.zeros(n_time)

    for t in range(n_time):

        A_t = latent_A[:, t, :]   # (n_trials, dims)
        B_t = latent_B[:, t, :]

        mean_diff = A_t.mean(axis=0) - B_t.mean(axis=0)

        # pooled covariance
        # pooled = np.vstack([A_t, B_t])
        # cov = np.cov(pooled, rowvar=False)
        
        #
        cov_A = np.cov(A_t, rowvar=False)
        cov_B = np.cov(B_t, rowvar=False)
        
        nA = A_t.shape[0]
        nB = B_t.shape[0]
        
        cov = ((nA - 1) * cov_A + (nB - 1) * cov_B) / (nA + nB - 2)

        # regularization for stability
        # cov += np.eye(cov.shape[0]) * 1e-6

        inv_cov = np.linalg.pinv(cov)

        dist[t] = np.sqrt(mean_diff.T @ inv_cov @ mean_diff)

    return dist
   

# def permutation_test_mahalanobis(latent_A, latent_B, n_perm=1000):

#     nA = latent_A.shape[0]
#     combined = np.concatenate([latent_A, latent_B], axis=0)

#     perm_dists = np.zeros((n_perm, latent_A.shape[1]))

#     for p in range(n_perm):

#         perm_idx = np.random.permutation(combined.shape[0])

#         perm_A = combined[perm_idx[:nA]]
#         perm_B = combined[perm_idx[nA:]]

#         perm_dists[p] = mahalanobis_distance_trials(perm_A, perm_B)

#     return perm_dists
    
def perm_testing(pop_trials_A, pop_trials_B, pca_space, my_title, stim_trials, 
                 baseline_means = [], baseline_stds = []):
    

    
    baseline_start = -1000 # samples in ms (1s) baselining makes sense before stim onset but not before response onset?
    baseline_end = 0
    baseline_mask = (time_centers >= baseline_start) & (time_centers <= baseline_end)
    
    all_distances = []
    latent_trials_A = []
    latent_trials_B = []
    
    # this block finds the mean of mah dist for each repeat, then averages them together
    for r in range(n_repeats):
    
        pop_A = pop_trials_A[r]
        pop_B = pop_trials_B[r]
    
        latent_A = get_latent(pop_A, pca_space, latent_dims) # for a single repeat, projects into space determined from all repeats
        latent_B = get_latent(pop_B, pca_space, latent_dims)
        
        latent_trials_A.append(latent_A)
        latent_trials_B.append(latent_B)
    
        dist = mahalanobis_distance_trials(latent_A, latent_B)
    
        # normalization
        if stim_trials:
            baseline = dist[baseline_mask]
            baseline_mean = baseline.mean()
            baseline_std  = baseline.std()
            
            baseline_means.append(baseline_mean)
            baseline_stds.append(baseline_std)
            
        else: 
            baseline_mean = baseline_means[r]
            baseline_std  = baseline_stds[r]
        
        dist = (dist - baseline_mean) / baseline_std
    
        all_distances.append(dist)

    # for plotting
    mean_dist = np.mean(all_distances, axis = 0) 
    sd_dist   = np.std(all_distances, axis = 0)
    
    # n_perm = 1000
    # perm_mean_dists = np.zeros((n_perm, len(time_centers)))
    
    # for p in range(n_perm): # not sure if this is correct
    
    #     perm_all_distances = []
    
    #     for r in range(n_repeats):
    
    #         pop_A = pop_trials_A[r]
    #         pop_B = pop_trials_B[r]
                
    #         # concatenate trials
    #         combined = np.concatenate([latent_trials_A[r], latent_trials_B[r]], axis=0)
    #         nA = latent_trials_A[r].shape[0]
    
    #         # shuffle trial labels
    #         perm_idx = np.random.permutation(combined.shape[0])
    #         perm_A = combined[perm_idx[:nA]]
    #         perm_B = combined[perm_idx[nA:]]
    
    #         dist = mahalanobis_distance_trials(perm_A, perm_B)
    
    #         # baseline = dist[baseline_mask] # which baseline to use?
            
    #         # baseline_mean = baseline.mean()
    #         # baseline_std  = baseline.std()
            
    #         dist = (dist - baseline_mean) / baseline_std
    
    #         perm_all_distances.append(dist)
    
    #     perm_all_distances = np.array(perm_all_distances)
    #     perm_mean_dists[p] = perm_all_distances.mean(axis=0)
    
        
    # p_vals = np.mean(perm_mean_dists >= mean_dist, axis=0)
    # sig_mask = p_vals < 0.05
    
    z_dists = np.array(all_distances)   # shape (200, time)

    threshold = 1.65
    
    p_vals = np.mean(z_dists <= threshold, axis=0)
    # p_vals = np.mean((z_dists >= 1.65) | (z_dists <= -1.65), axis=0)
    
    sig_mask = p_vals < 0.05
    
    min_cluster_bins = int(0.0 / bin_size) # abitrary,jercog does point wise (0.0s)

    sig_mask_filtered = np.zeros_like(sig_mask)
    current_cluster = []
    
    for i, val in enumerate(sig_mask):
        if val:
            current_cluster.append(i)
        else:
            if len(current_cluster) >= min_cluster_bins:
                sig_mask_filtered[current_cluster] = True
            current_cluster = []
    
    sig_mask = sig_mask_filtered
    
    
    plt.figure(figsize=(8,4))
    
    plt.plot(time_centers, mean_dist, linewidth=2)
    plt.fill_between(time_centers,
                     mean_dist - sd_dist,
                     mean_dist + sd_dist,
                     alpha=0.3)
    
    plt.axvline(0, linestyle='--')
    
    sig_indices = np.where(sig_mask)[0]

    if len(sig_indices) > 0:
    
        start = sig_indices[0]
    
        for i in range(1, len(sig_indices)):
    
            # detect gap (end of a significant cluster)
            if sig_indices[i] != sig_indices[i-1] + 1:
    
                end = sig_indices[i-1]
    
                plt.hlines(
                    y=mean_dist.max() * 1.05,
                    xmin=time_centers[start],
                    xmax=time_centers[end],
                    linewidth=3
                )
    
                start = sig_indices[i]
    
        # final segment
        end = sig_indices[-1]
    
        plt.hlines(
            y=mean_dist.max() * 1.05,
            xmin=time_centers[start],
            xmax=time_centers[end],
            linewidth=3
        )
        
    # # mark significant timepoints
    # plt.scatter(time_centers[sig_mask],
    #             mean_dist[sig_mask],
    #             s=15)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Normalized distance (d)')
    plt.title(f'Mean ± s.d. normalized trajectory distance, {my_title}')
    plt.tight_layout()
    plt.show()

    return np.array(baseline_means), np.array(baseline_stds)

def plot_trajectory_3d(ax, traj, color_pre, color_post, label):
    """Plot a 3D trajectory with different colors pre/post stimulus onset."""
    
    # Project to 2D for LineCollection, then manually draw in 3D
    # Split trajectory at onset
    pre  = traj[:onset_idx + 1]
    post = traj[onset_idx:]

    ax.plot(pre[:, 0],  pre[:, 1],  pre[:, 2],
            linewidth=2, color=color_pre,  linestyle='--', alpha=0.5)
    ax.plot(post[:, 0], post[:, 1], post[:, 2],
            linewidth=2, color=color_post, label=label)
    
    # Mark stimulus onset point
    ax.scatter(traj[onset_idx, 0], traj[onset_idx, 1], traj[onset_idx, 2],
               color=color_post, s=40, zorder=5)

    
# -----------------------------
# Load & pool neurons across sessions
# -----------------------------
all_neuron_trials = []
 
# stimulus aligned
all_neuron_trials_freeze = []
all_neuron_trials_nonfreeze = []

# response aligned
all_neuron_trials_freeze_aligned = []
all_neuron_trials_nonfreeze_shuff = []

if recalc_freeze_times:
    pooled_freeze = []
else:
   # pooled_freeze = pd.read_csv(f"{freeze_file_path}\\freeze_times.csv", header = None)
    with open(f"{freeze_file_path}\\freeze_times.csv") as f:
        reader = csv.reader(f)
        pooled_freeze = [float(row[0]) for row in reader]

for sess in session_list:
    mat = sio.loadmat(
        os.path.join(data_path, sess),
        struct_as_record=False,
        squeeze_me=True
    )

    sAP = mat['sAP']
    
    # find_field(sAP, 'vecIsActive')
    
    vecIsActive = sAP.cellBlock[0].vecIsActive


    clusters = sAP.sCluster

    is_good = np.array([c.KilosortGood for c in clusters]) # try no filtering
    low_contam = np.array([c.Contamination for c in clusters]) < 0.1
    areas = np.array([str(c.Area).strip().lower() for c in clusters])
    in_area = np.array([region in a for a in areas])

    idx_keep = (is_good | low_contam) & in_area

    stim_onset = sAP.cellBlock[0].vecStimOnTime
    stim_off = sAP.cellBlock[0].vecStimOffTime
    stim_duration_trial = stim_off - stim_onset

    # speed data
    run_speed = sAP.PP_GetRunSpeed.vecSpeed_mps
    event_idx = (stim_onset * fs).astype(int)
    peri_event_speed = []
    valid_stim_onset = []
    valid_event_mask = []

    for idx in event_idx:
        start = idx - n_pre
        end   = idx + n_post

        # Skip events too close to start or end
        if start < 0 or end > len(run_speed):
            valid_event_mask.append(False)
            continue
        valid_event_mask.append(True)
        peri_event_speed.append(run_speed[start:end])
        

    peri_event_speed = np.array(peri_event_speed)
    valid_stim_onset = stim_onset[valid_event_mask]
    
    freeze_onsets = []  # list of lists (one list per event) # delete???
    freeze_mask = []
    active_mask = []
    freeze_onsets_sec = []

    for trial_velocity in peri_event_speed:
    
        # find immobility within trial
        immobile = trial_velocity < freeze_thresh
    
        # Find contiguous immobile periods
        diff = np.diff(immobile.astype(int))
        starts = np.where(diff == 1)[0] + 1
        ends   = np.where(diff == -1)[0] + 1
    
        # Handle edge cases
        if immobile[0]:
            starts = np.insert(starts, 0, 0)
        if immobile[-1]:
            ends = np.append(ends, len(immobile))
    
        trial_freeze_onsets = []
    
        for s, e in zip(starts, ends):
            duration = e - s
        
            if duration >= min_freeze_samples:
        
                # Convert to time relative to event onset
                onset_rel = s - n_pre  # in samples
                onset_sec = onset_rel / fs
        
                # Apply exclusion criteria
                if onset_rel >= exclude_samples and onset_rel <= event_end_sample: #should change event_end_sample to value specific to trial
                    
                    # Check for sufficient movement before freeze onset
                    pre_freeze_samples = int(pre_freeze_window * fs)
                    pre_freeze_start = max(0, s - pre_freeze_samples) #should change to average perhaps, ok for now
                    pre_freeze_velocity = trial_velocity[pre_freeze_start:s]
                    
                    # Require that peak speed in pre-freeze window exceeds threshold
                    if len(pre_freeze_velocity) > 0 and np.max(pre_freeze_velocity) >= pre_freeze_movement_thresh:
                        trial_freeze_onsets.append(onset_sec)
                        
        # print(f"{trial_freeze_onsets}")    
        
        if len(trial_freeze_onsets) > 1:
            print('multiple freeze onsets detected')
        
            
        # find decelerration point, work in progress
        if len(trial_freeze_onsets) == 1:
            freeze_point = int(trial_freeze_onsets[0]*fs + n_pre)
            exclusion_window = 500 # ms

            x = np.arange(len(trial_velocity[:freeze_point]))
            # knee = detect_knees(x, trial_velocity[:freeze_point][::-1],
                                        # exclusion_window)
            
            # decel_sample = int(x[-1] - (knee + exclusion_window))
    
            # acceleration = np.diff(trial_velocity)
            
            # # smooth acceleration
            # acceleration = gaussian_filter1d(acceleration, 3)
                        
            # # decel_sample = None
            
            # for i in range(freeze_point - exclusion_window, 0, -1):
            #     if np.all(acceleration[i:i+exclusion_window] < 0):
            #         decel_sample = i
            #         break
                
            # # decel_sample

            if plot_speed_trials:    
                plt.figure()
                plt.plot(gaussian_filter1d(trial_velocity,10))
                plt.axvline(freeze_point, color='r', linestyle='--', label='Freeze onset')
                # plt.axvline(decel_sample, color='g', linestyle='--', label='Deceleration onset')
                # plt.xlim([decel_sample-1000,freeze_point+1000])
                plt.ylabel('Speed')
                plt.xlabel('Time')
                plt.title('Speed during trial with freeze and deceleration')
                plt.legend()
        
        freeze_onsets.append(trial_freeze_onsets)
       
        # freeze_mask.append(len(trial_freeze_onsets) > 0)
        
        has_freeze = len(trial_freeze_onsets) > 0
        freeze_mask.append(has_freeze)
        
        if has_freeze:
            # take first freeze only (you already enforce single freeze earlier)
            freeze_onsets_sec.append(trial_freeze_onsets[0])
        
        # ---------------------------------------
        # 2. Detect movement-only (no freeze)
        # ---------------------------------------
        has_active_movement = False
    
        if not has_freeze:
            movement = trial_velocity[:int((s_pre+stim_length)*1000)] >= pre_freeze_movement_thresh
            diff_move = np.diff(movement.astype(int))
            move_starts = np.where(diff_move == 1)[0] + 1
            move_ends = np.where(diff_move == -1)[0] + 1
    
            if movement[0]:
                move_starts = np.insert(move_starts, 0, 0)
            if movement[-1]:
                move_ends = np.append(move_ends, len(movement))
    
            min_move_samples = int(pre_freeze_window * fs)
    
            for s, e in zip(move_starts, move_ends):
                duration = e - s
                if duration >= min_move_samples:
    
                    has_active_movement = True
                    break
    
        active_mask.append(has_active_movement)
        
    # session_freeze_mask = [len(onsets) > 0 for onsets in freeze_onsets]
    if recalc_freeze_times:
        pooled_freeze.append(freeze_onsets_sec)
        
    freeze_mask = np.array(freeze_mask)
    active_mask = np.array(active_mask)

    freeze_event_times = valid_stim_onset[freeze_mask]
    active_event_times = valid_stim_onset[active_mask]
    nonresponse_event_times = valid_stim_onset[~freeze_mask & ~active_mask] 
    
    
    # need to make trials for stopping ('psuedo-freezing') during ITIs
    
    freeze_onsets_sec = np.array(freeze_onsets_sec)
    freeze_aligned_times = freeze_event_times + freeze_onsets_sec
    
    active_offsets = np.full(len(valid_stim_onset), np.nan)

    if len(freeze_onsets_sec) > 0 and len(active_event_times) > 0 and recalc_freeze_times:
        shuffled_offsets = np.random.choice(
            freeze_onsets_sec,
            size=len(active_event_times),
            replace=True
        )
        nonfreeze_shuffled_times = active_event_times + shuffled_offsets
        active_offsets[active_mask] = shuffled_offsets
        
    elif len(active_event_times) > 0 and not recalc_freeze_times:
        shuffled_offsets = np.random.choice(
            pooled_freeze,
            size=len(active_event_times),
            replace=True
        )
        nonfreeze_shuffled_times = active_event_times + shuffled_offsets
        active_offsets[active_mask] = shuffled_offsets

    else:
        shuffled_offsets = np.array([])
        nonfreeze_shuffled_times = np.array([]) 

    
    if plot_speed_heatmap:
        ts = np.linspace(-s_pre, s_post, win_len)
    
        plt.figure(figsize=(6, 8))
        plt.imshow(
            peri_event_speed,
            aspect='auto',
            cmap='viridis',
            interpolation='none',
            origin='lower',
            extent=[ts[0], ts[-1], 0, peri_event_speed.shape[0]]
        )
    
        plt.colorbar(label='Run speed')
        plt.axvline(0, color='white', linestyle='--', linewidth=1)
        plt.axvline(0+stim_duration_trial[0], color='white', linestyle='--', linewidth=1)
        
        # Plot freeze onset markers
        for trial_idx, trial_freezes in enumerate(freeze_onsets):
        
            for onset_sec in trial_freezes:
        
                # Draw vertical dashed red line limited to this trial row
                plt.vlines(
                    onset_sec,
                    trial_idx,
                    trial_idx + 1,
                    colors='red',
                    linestyles='dashed',
                    linewidth=1.5
                )
                
   # Plot shuffled offsets for active trials
        for trial_idx, offset in enumerate(active_offsets):
        
            if not np.isnan(offset):
        
                plt.vlines(
                    offset,
                    trial_idx,
                    trial_idx + 1,
                    colors='white',
                    linestyles='solid',
                    linewidth=1.5
                )
  
    
        plt.xlabel('Time from event (s)')
        plt.ylabel('Trial')
        plt.title('Peri-event run speed')
    
        plt.show()
           
    print(f'Number of recorded neurons this session: {sum(idx_keep)}')
    
    for i, c in enumerate(clusters):
    
        if not idx_keep[i]:
            continue
        
        
        # trials = neuron_trial_matrix(c.SpikeTimes, stim_onset)
        # if trials.shape[0] > 5:  # minimal trial count. do i need to filter this? why 5?
        #     all_neuron_trials.append(trials)
            
        # -------------------------------
        # Stimulus-aligned 
        # -------------------------------
        trials_freeze = neuron_trial_matrix(
            c.SpikeTimes,
            freeze_event_times
        )
    
        trials_nonfreeze = neuron_trial_matrix(
            c.SpikeTimes,
            active_event_times
        )
    
        # -------------------------------
        # Freeze-onset aligned 
        # -------------------------------
        trials_freeze_aligned = neuron_trial_matrix(
            c.SpikeTimes,
            freeze_aligned_times
        )
    
        trials_nonfreeze_shuff = neuron_trial_matrix(
            c.SpikeTimes,
            nonfreeze_shuffled_times
        )
    
        # --------------------------------
        # Keep neuron only if all types exist; is this needed? is it ok to have neurons that only were recorded during one respones type?
        # --------------------------------
        if (trials_freeze.shape[0] > 2 and
            trials_nonfreeze.shape[0] > 2):
    
            all_neuron_trials_freeze.append(trials_freeze)
            all_neuron_trials_nonfreeze.append(trials_nonfreeze)
    
            all_neuron_trials_freeze_aligned.append(trials_freeze_aligned)
            all_neuron_trials_nonfreeze_shuff.append(trials_nonfreeze_shuff)

    print(f"Neurons x freeze trials: {sum(idx_keep)} x {sum(freeze_mask)}")
    print(f"Neurons x nonfreeze trials: {sum(idx_keep)} x {sum(active_mask)}")

if recalc_freeze_times:
    flattened_pooled_freeze = [x 
                               for xs in pooled_freeze
                               for x in xs]
    df = pd.DataFrame(flattened_pooled_freeze)
    df.to_csv(f"{freeze_file_path}\\freeze_times.csv", header=False, index=False)

    
# -----------------------------
# PCA — fit on all trials combined, project each type separately # change to gfpa!!!
# -----------------------------
population_freeze    = generate_population_trials_by_repeat(all_neuron_trials_freeze,    n_pop_trials, n_repeats)
population_nonfreeze = generate_population_trials_by_repeat(all_neuron_trials_nonfreeze, n_pop_trials, n_repeats)

X_freeze    = flatten_for_pca(population_freeze)
X_nonfreeze = flatten_for_pca(population_nonfreeze)
X_combined  = np.vstack([X_freeze, X_nonfreeze])



pca = PCA(n_components=latent_dims)
pca.fit(X_combined)

latent_freeze    = get_latent(population_freeze,    pca, latent_dims)
latent_nonfreeze = get_latent(population_nonfreeze, pca, latent_dims)

trajectory_freeze    = latent_freeze.mean(axis=1) # mean per repeat
trajectory_nonfreeze = latent_nonfreeze.mean(axis=1)


# -----------------------------
# Visualization: 3D trajectories — freeze vs nonfreeze
# -----------------------------
# Find stimulus onset index in trajectory time
if plot_trajectories:
    population_freeze_gen    = generate_population_trials(all_neuron_trials_freeze,    n_pop_trials, n_repeats)
    population_nonfreeze_gen = generate_population_trials(all_neuron_trials_nonfreeze, n_pop_trials, n_repeats)
    
    X_freeze_gen    = flatten_for_pca(population_freeze_gen)
    X_nonfreeze_gen = flatten_for_pca(population_nonfreeze_gen)
    X_combined_gen  = np.vstack([X_freeze_gen, X_nonfreeze_gen])

    # Fit PCA on combined data so both trial types share the same latent space
    pca_gen = PCA(n_components=latent_dims)
    pca_gen.fit(X_combined_gen)

    latent_freeze_gen    = get_latent(population_freeze_gen,    pca_gen, latent_dims)
    latent_nonfreeze_gen = get_latent(population_nonfreeze_gen, pca_gen, latent_dims)

    trajectory_freeze_gen    = latent_freeze_gen.mean(axis=0) # mean per repeat
    trajectory_nonfreeze_gen = latent_nonfreeze_gen.mean(axis=0)
    
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    onset_idx = np.argmin(np.abs(time_centers - 0.0))
    
    ax.plot(trajectory_freeze_gen[:, 0],    trajectory_freeze_gen[:, 1],    trajectory_freeze_gen[:, 2],
            linewidth=2, color='crimson',   label='Freeze')
    ax.plot(trajectory_nonfreeze_gen[:, 0], trajectory_nonfreeze_gen[:, 1], trajectory_nonfreeze_gen[:, 2],
            linewidth=2, color='steelblue', label='Non-freeze')
    
    # Trial start dots
    ax.scatter(*trajectory_freeze_gen[0, :3],    color='black', s=40, zorder=5)
    ax.scatter(*trajectory_nonfreeze_gen[0, :3], color='black', s=40, zorder=5)
    
    # Stimulus onset dots
    ax.scatter(*trajectory_freeze_gen[onset_idx, :3],    color='black', s=40, marker='D', zorder=5)
    ax.scatter(*trajectory_nonfreeze_gen[onset_idx, :3], color='black', s=40, marker='D', zorder=5)
    
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    ax.set_title(f'{region} Population Trajectory — Stimulus aligned')
    ax.legend()
    plt.tight_layout()
    plt.show()


#
# response aligned
#
population_freeze_aligned = generate_population_trials_by_repeat(
    all_neuron_trials_freeze_aligned,
    n_pop_trials,
    n_repeats
)

population_nonfreeze_shuff = generate_population_trials_by_repeat(
    all_neuron_trials_nonfreeze_shuff,
    n_pop_trials,
    n_repeats
)

X_freeze_align    = flatten_for_pca(population_freeze_aligned)
X_nonfreeze_shuff = flatten_for_pca(population_nonfreeze_shuff)
X_combined_align  = np.vstack([X_freeze_align, X_nonfreeze_shuff])

# Fit PCA on combined data so both trial types share the same latent space

# # Fit PCA on combined data so both trial types share the same latent space
# pca_full = PCA()
# pca_full.fit(X_combined)

# explained_var = pca_full.explained_variance_ratio_
# cum_var = np.cumsum(explained_var)
# latent_dims = np.argmax(cum_var >= 0.9) + 1

# # Scree plot
# plt.figure()
# plt.plot(range(1, len(explained_var) + 1), explained_var, marker='o', label='Individual variance')
# # plt.plot(range(1, len(cum_var) + 1), cum_var, marker='s', label='Cumulative variance')

# # plt.axhline(0.9, linestyle='--', label='90% variance threshold')

# plt.xlabel("Principal Component")
# plt.ylabel("Explained Variance Ratio")
# plt.title("PCA Scree Plot")
# plt.legend()
# plt.tight_layout()
# plt.show()

# print(latent_dims)

pca_align = PCA(n_components=latent_dims)
pca_align.fit(X_combined_align)


latent_freeze_aligned = get_latent(
    population_freeze_aligned,
    pca_align,
    latent_dims
)

latent_nonfreeze_shuff = get_latent(
    population_nonfreeze_shuff,
    pca_align,
    latent_dims
)





# plt.figure(figsize=(8,4))

# plt.plot(time_centers, mean_dist, linewidth=2)
# plt.fill_between(time_centers,
#                  mean_dist - sd_dist,
#                  mean_dist + sd_dist,
#                  alpha=0.3)

# plt.axvline(0, linestyle='--')
# plt.xlabel('Time (s)')
# plt.ylabel('Normalized distance (d)')
# plt.title('Mean ± s.d. normalized trajectory distance')
# plt.tight_layout()
# plt.show()

trajectory_freeze_aligned = latent_freeze_aligned.mean(axis=0)
trajectory_nonfreeze_shuff = latent_nonfreeze_shuff.mean(axis=0)


# plotting trajectories
if plot_trajectories:

    population_freeze_align_gen    = generate_population_trials(all_neuron_trials_freeze_aligned,    n_pop_trials, n_repeats)
    population_nonfreeze_shuff_gen = generate_population_trials(all_neuron_trials_nonfreeze_shuff, n_pop_trials, n_repeats)
    
    latent_freeze_align_gen    = get_latent(population_freeze_align_gen,    pca_gen, latent_dims)
    latent_nonfreeze_shuff_gen = get_latent(population_nonfreeze_shuff_gen, pca_gen, latent_dims)

    trajectory_freeze_align_gen    = latent_freeze_align_gen.mean(axis=0) # mean per repeat
    trajectory_nonfreeze_shuff_gen = latent_nonfreeze_shuff_gen.mean(axis=0)
    
    freeze_onset_idx = np.argmin(np.abs(time_centers - 0.0))
    
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot(trajectory_freeze_align_gen[:,0],
            trajectory_freeze_align_gen[:,1],
            trajectory_freeze_align_gen[:,2],
            color='darkred',
            linewidth=2,
            label='Freeze (aligned to freeze)')
    
    ax.plot(trajectory_nonfreeze_shuff_gen[:,0],
            trajectory_nonfreeze_shuff_gen[:,1],
            trajectory_nonfreeze_shuff_gen[:,2],
            color='gray',
            linewidth=2,
            label='Non-freeze (shuffled freeze time)')
    
    
    freeze_onset_idx = np.argmin(np.abs(time_centers - 0.0))
    trial_start_idx = 0   # first time bin
    
    ax.scatter(*trajectory_freeze_align_gen[trial_start_idx, :3],
               color='black', s=60, marker='o', zorder=5)
    
    ax.scatter(*trajectory_nonfreeze_shuff_gen[trial_start_idx, :3],
               color='black', s=60, marker='o', zorder=5)
    
    
    ax.scatter(*trajectory_freeze_align_gen[freeze_onset_idx, :3],
               color='black', s=60, marker='D', zorder=6)
    
    ax.scatter(*trajectory_nonfreeze_shuff_gen[freeze_onset_idx, :3],
               color='black', s=60, marker='X', zorder=6)
    
    
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    ax.set_title(f'{region} Population Trajectories — Freeze Onset Aligned')
    ax.legend()
    
    plt.tight_layout()
    plt.show()
   
    
# perm testing
pre_stim_mean, pre_stim_std = perm_testing(population_freeze, population_nonfreeze, pca, 
             my_title = f'Stimulus aligned, {region}', stim_trials = True)

perm_testing(population_freeze_aligned, population_nonfreeze_shuff, pca_align, 
             my_title = f'Response aligned, {region}', stim_trials = False,
             baseline_means = pre_stim_mean, baseline_stds = pre_stim_std)


# distances 

# real_dist = mahalanobis_distance_trials(
#     latent_freeze_aligned,
#     latent_nonfreeze_shuff
# )

# perm_dists = permutation_test_mahalanobis(
#     latent_freeze_aligned,
#     latent_nonfreeze_shuff,
#     n_perm=1000
# )

# # get p values
# p_vals = np.mean(perm_dists >= real_dist, axis=0)

# # correct for spurious false positives via clustering 
# alpha = 0.05
# sig_mask = p_vals < alpha

# # compute max cluster size under permutation
# max_cluster_sizes = []

# for p in range(perm_dists.shape[0]):
#     perm_sig = perm_dists[p] >= np.percentile(perm_dists, 95, axis=0)
    
#     clusters = np.split(perm_sig, np.where(~perm_sig)[0])
#     cluster_sizes = [len(c) for c in clusters if np.all(c)]
    
#     max_cluster_sizes.append(max(cluster_sizes) if len(cluster_sizes) > 0 else 0)

# cluster_threshold = np.percentile(max_cluster_sizes, 95)

# # plot results
# plt.figure(figsize=(8,4))

# plt.plot(time_centers, real_dist, color='black', linewidth=2)
# plt.fill_between(time_centers,
#                  np.percentile(perm_dists, 2.5, axis=0),
#                  np.percentile(perm_dists, 97.5, axis=0),
#                  alpha=0.3)

# plt.axvline(0, linestyle='--')
# plt.xlabel('Time (s)')
# plt.ylabel('Mahalanobis distance')
# plt.title('Freeze-aligned trajectory separation')
# plt.tight_layout()
# plt.show()