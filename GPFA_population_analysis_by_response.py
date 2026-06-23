# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:29:51 2026

Calculates neural trajectories and mahalonobis distances between respones types.

@author: sconrad
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy.ndimage import gaussian_filter1d
from scipy.stats import sem
# from scipy.spatial.distance import mahalanobis
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import pandas as pd
import neo
from elephant.gpfa import GPFA
import quantities as pq
from sklearn.model_selection import cross_val_score
import pickle

# -----------------------------
# Parameters
# -----------------------------

# user input

# plot_loadings = False

analyze_stim = False
analyze_resp = True
include_ITI_analysis = False

# recalc_stim = False

shuffle_trial_labels = True
find_latent_dims = False

project_train_data = False # if you want to project original training data

# region = 'periaqueductal gray' # or ~7-10?
region = 'superior colliculus' # or stim = 10? response = 13
# region = 'retrosplenial'   # or 
# region = 'midbrain'   

print(f'\nRunning {region} analysis with settings:\nShuffle labels:'
      f' {shuffle_trial_labels}\nFinding Latent dimensions: {find_latent_dims}\n'
      f'Analyze Stimulus Aligned: {analyze_stim}\n'
      f'Analyze Response Aligned: {analyze_resp}\n')
if analyze_resp:
    print(f'Analyzing ITI in Response Space: {include_ITI_analysis}\n')

spike_path = r'W:\Haak\Innate_defense\Data_analysis\22.35.02\Scotts_analysis\collected_spikes.pickle'

collected_spikes = pd.read_pickle(spike_path)

all_neuron_trials_freeze = collected_spikes[region]['all_neuron_trials_freeze']
all_neuron_trials_nonfreeze = collected_spikes[region]['all_neuron_trials_nonfreeze']
all_neuron_trials_freeze_aligned = collected_spikes[region]['all_neuron_trials_freeze_aligned']
all_neuron_trials_nonfreeze_shuff = collected_spikes[region]['all_neuron_trials_nonfreeze_shuff']
all_neuron_trials_ITI_freeze = collected_spikes[region]['all_neuron_trials_ITI_freeze']
all_neuron_trials_ITI_nonfreeze_shuff = collected_spikes[region]['all_neuron_trials_ITI_nonfreeze_shuff']

# general
np.random.seed(1)

# spike analysis initialize
t_pre = 1.0 # s
t_post = 4.0 # s
bin_size = 0.02  # 20 ms
# smooth_sigma = 1.0  # in bins

n_pop_trials = 30
n_repeats = 200  # increase to 200 to match Jercog 2021

if region == 'periaqueductal gray':
    latent_dims = 8 #10
else:
    latent_dims = 8 # 14 or 5? ; #SC

time_bins = np.arange(-t_pre, t_post + bin_size, bin_size)
time_centers = time_bins[:-1] + bin_size / 2

# speed analysis initialize
# params are probably buried in sAP but i dont know where yet
# fs = 1000  # Hz
# s_pre = 4    # seconds before event onset
# s_post = 6  # seconds after event onset


# -----------------------------
# Helper functions
# -----------------------------


def split_trials_per_neuron(all_neuron_trials):
    """
    Parameters
    ----------
    all_neuron_trials : TYPE
        DESCRIPTION.

    Returns
    -------
    fit_list : TYPE
        DESCRIPTION.
    traj_list : TYPE
        DESCRIPTION.

    """
    fit_list = []
    traj_list = []

    for neuron_trials in all_neuron_trials:
        n_trials = neuron_trials.shape[0]

        # shuffle indices
        indices = np.random.permutation(n_trials)

        half = n_trials // 2

        fit_idx = indices[:half]
        traj_idx = indices[half:]

        fit_list.append(neuron_trials[fit_idx])
        traj_list.append(neuron_trials[traj_idx])

    return fit_list, traj_list


def sample_iPV_trajectory(neuron_trials): # samples random time trial per time point
    """
    neuron_trials: list of arrays (n_trials × n_time)
    
    Returns:
        trajectory: (n_time × n_neurons)
    """
    
    n_neurons = len(neuron_trials)
    n_time = neuron_trials[0].shape[1]
    
    traj = np.zeros((n_time, n_neurons))
    
    for n, trials in enumerate(neuron_trials):
        
        n_trials = trials.shape[0]
        
        # sample a trial independently at each time bin
        sampled_trials = np.random.randint(0, n_trials, size=n_time)
        
        traj[:, n] = trials[sampled_trials, np.arange(n_time)]
        
        # for t in range(n_time):
        #     traj[t, n] = trials[sampled_trials[t], t]
    
    return traj

# def sample_iPV_trajectory(neuron_trials): # samples only one random trial 

#     n_neurons = len(neuron_trials)
#     n_time = neuron_trials[0].shape[1]

#     traj = np.zeros((n_time, n_neurons))

#     for n, trials in enumerate(neuron_trials):

#         trial_idx = np.random.randint(0, trials.shape[0])

#         traj[:, n] = trials[trial_idx]

#     return traj

def generate_iPV_trajectories(neuron_trials,
                              n_pop_trials):
    """
    neuron_trials shape:
        list of n number of x by y trial matrices, where
        n = number of neurons
        x = number of trials
        y = time

    n_pop_trials: number of desired population vectors to generate
    """
    # pop_trials = []
    # for p in range(n_pop_trials):

    #     traj = sample_iPV_trajectory(neuron_trials)
    #     pop_trials.append(traj)


    # return np.array(pop_trials)

    n_time = neuron_trials[0].shape[1]
    n_neurons = len(neuron_trials)

    pop = np.zeros((n_pop_trials, n_time, n_neurons))

    for p in range(n_pop_trials):
        pop[p] = sample_iPV_trajectory(neuron_trials)

    return pop

def remove_silent_neurons(pop, var_threshold=1e-10):
    """
    pop shape:
        (n_trials, n_time, n_neurons)

    Removes neurons with near-zero variance.
    """

    # flatten across trials and time
    flat = pop.reshape(-1, pop.shape[2])

    neuron_var = np.var(flat, axis=0)

    keep = neuron_var > var_threshold

    print(f"Keeping {keep.sum()} / {len(keep)} neurons")

    return pop[:, :, keep], keep

def flatten_repeats(pop_repeats):
    """
    pop_repeats:
        list of repeats

    Returns:
        single list of neo spike train trials
    """

    all_trials = []

    for pop in pop_repeats:

        formatted = format_ipv_for_gpfa(pop)

        all_trials.extend(formatted)

    return all_trials

def counts_to_spiketrains(trial_array, bin_size_s, t_start_s):
    """
    Convert a binned count array to a list of neo.SpikeTrain objects.
    
    trial_array: (n_neurons, n_time) — integer-ish spike counts per bin
    bin_size_s:  float, bin size in seconds
    t_start_s:   float, time of first bin edge in seconds (e.g. -t_pre)
    
    Returns: list of n_neurons neo.SpikeTrain objects
    """
    n_neurons, n_time = trial_array.shape
    t_stop = t_start_s + n_time * bin_size_s
    spiketrains = []

    for neuron_counts in trial_array:
        spike_times = []

        for t_idx, count in enumerate(neuron_counts):
            # Place 'count' spikes in the center of each bin
            bin_center = t_start_s + (t_idx + 0.5) * bin_size_s
            # Round to nearest integer count (iPV values may be floats 
            # due to variance-stabilizing sqrt transform)
            n_spikes = int(np.round(count))
            spike_times.extend([bin_center] * n_spikes)

        st = neo.SpikeTrain(
            spike_times * pq.s,
            t_start=t_start_s * pq.s,
            t_stop=t_stop * pq.s
        )
        spiketrains.append(st)

    return spiketrains


def format_ipv_for_gpfa(pop):
    """
    pop: (n_pop_trials, n_time, n_neurons)
    Returns: list of length n_pop_trials, each element is a 
             list of n_neurons neo.SpikeTrain objects
    """
    trials = []
    for r in range(pop.shape[0]):
        trial_array = pop[r].T          # (n_neurons, n_time)
        sts = counts_to_spiketrains(trial_array, bin_size, -t_pre)
        trials.append(sts)
    return trials

# def find_latent_dimensions(fit_trials, alignment_type, min_dim=1, max_dim = 15):
#     """
#     Returns: number of latent dimensions that minimizes GPFA prediction error
#     """
#     x_dims = list(range(min_dim,max_dim+1))
#     log_likelihoods = []
#     for x_dim in x_dims:
#         gpfa_cv = GPFA(x_dim=x_dim)
#         # estimate the log-likelihood for the given dimensionality as the mean of the log-likelihoods from 3 cross-vailidation folds
#         cv_log_likelihoods = cross_val_score(gpfa_cv, fit_trials, cv=4, n_jobs=4, verbose=True)
#         log_likelihoods.append(np.mean(cv_log_likelihoods))
    
#     plt.figure()
#     plt.xlabel(f'Dimensionality of latent variables for {region}, {alignment_type}')
#     plt.ylabel('Log-likelihood')
#     plt.plot(x_dims, log_likelihoods, '.-')
#     plt.plot(x_dims[np.argmax(log_likelihoods)], np.max(log_likelihoods), 'x', markersize=10, color='r')
#     plt.tight_layout()
#     plt.show()
    
    
#     return x_dims[np.argmax(log_likelihoods)]

def find_latent_dimensions(
        train_trials,
        test_trials,
        alignment_type,
        min_dim=1,
        max_dim=20):

    x_dims = list(range(min_dim, max_dim + 1))

    prediction_errors = []

    for x_dim in x_dims:

        print(f'Fitting x_dim = {x_dim}')

        gpfa_model = GPFA(
            bin_size=bin_size * pq.s,
            x_dim=x_dim
        )

        gpfa_model.fit(train_trials)

        error = gpfa_prediction_error(
            gpfa_model,
            test_trials,
            latent_dims_used=x_dim
        )

        prediction_errors.append(error)

    prediction_errors = np.array(prediction_errors)

    best_dim = x_dims[np.argmin(prediction_errors)]

    plt.figure()

    plt.plot(
        x_dims,
        prediction_errors,
        '.-'
    )

    plt.plot(
        best_dim,
        prediction_errors[np.argmin(prediction_errors)],
        'rx',
        markersize=10
    )

    plt.xlabel(
        f'Latent dimensionality ({alignment_type})'
    )

    plt.ylabel('Prediction error')

    plt.tight_layout()
    plt.show()

    return best_dim


def trials_to_matrix(pop):
    """
    pop shape:
        (n_trials, n_time, n_neurons)

    Returns:
        Y shape:
        (n_trials, n_neurons, n_time)
    """
    return np.transpose(pop, (0, 2, 1))


def gpfa_prediction_error(
        gpfa_model,
        test_trials,
        latent_dims_used=None):
    """
    Computes GPFA leave-neuron-out prediction error
    following Yu et al. 2009.

    Parameters
    ----------
    gpfa_model : fitted elephant GPFA object

    test_pop :
        shape (n_trials, n_time, n_neurons)

    latent_dims_used :
        optional reduced dimensionality p~
    """

    Y = trials_to_matrix(test_trials)

    n_trials, n_neurons, n_time = Y.shape

    total_error = 0.0

    # loading matrix
    C = gpfa_model.params_estimated['C']
    d = gpfa_model.params_estimated['d'].flatten()

    # orthonormalization
    U, S, Vt = np.linalg.svd(C, full_matrices=False)

    if latent_dims_used is None:
        latent_dims_used = C.shape[1]

    for trial_idx in range(n_trials):

        trial = test_trials[trial_idx]

        #
        # infer latent trajectory using all neurons
        #
        trial_fmt = format_ipv_for_gpfa(
            np.expand_dims(trial, axis=0)
        )

        latent = gpfa_model.transform(trial_fmt)[0]

        # latent shape:
        # (x_dim, n_time)

        #
        # orthonormalized trajectory
        #
        x_orth = (
            np.diag(S[:latent_dims_used])
            @ Vt[:latent_dims_used]
            @ latent
        )

        for j in range(n_neurons):

            # reduced loading vector
            u_j = U[j, :latent_dims_used]

            # prediction
            yhat_j = (
                u_j @ x_orth
            ) + d[j]

            # observed activity
            y_j = Y[trial_idx, j]

            # squared error
            total_error += np.sum(
                (yhat_j - y_j) ** 2
            )

    return total_error


def mahalanobis_distance_trials(latent_A, latent_B, repeat, stim_trials = True, 
                                baseline_means = [], baseline_stds = []):
    """
    latent_A, latent_B:
        shape = (n_trials, n_timepoints, latent_dims)
    Returns:
        dist[t] = Mahalanobis distance between condition means at time t
        
    """
    latent_A = np.stack(latent_A, axis=0)
    latent_B = np.stack(latent_B, axis=0)

    n_trials_A, n_dims, n_time = latent_A.shape
    
    dist = np.zeros(n_time)

    for t in range(n_time):

        A_t = latent_A[:, :, t]   # (n_trials, dims)
        B_t = latent_B[:, :, t]

        mean_diff = A_t.mean(axis=0) - B_t.mean(axis=0)
        
        cov_A = np.cov(A_t, rowvar=False)
        cov_B = np.cov(B_t, rowvar=False)
        
        nA = A_t.shape[0]
        nB = B_t.shape[0]
        
        cov = ((nA - 1) * cov_A + (nB - 1) * cov_B) / (nA + nB - 2)

        # regularization for stability
        # cov += np.eye(cov.shape[0]) * 1e-6

        inv_cov = np.linalg.pinv(cov)

        dist[t] = np.sqrt(mean_diff.T @ inv_cov @ mean_diff)
    
    
    # # Z-score normalization to baseline
    # baseline_start = -1000 # samples in ms (1s) baselining makes sense before stim onset but not before response onset?
    # baseline_end = 0
    # baseline_mask = (time_centers >= baseline_start) & (time_centers <= baseline_end)
    # if stim_trials:
    #     baseline = dist[baseline_mask]
    #     baseline_mean = baseline.mean()
    #     baseline_std  = baseline.std()
        
    #     baseline_means.append(baseline_mean)
    #     baseline_stds.append(baseline_std)
        
    # else: 
    #     baseline_mean = baseline_means[repeat]
    #     baseline_std  = baseline_stds[repeat]
    
    
    # dist = (dist - baseline_mean) / baseline_std # z-scored


    return dist 

def perm_testing2(all_distances, my_title):
    
    
    # Z-score/extract max  +avg etc here
    baseline = 50 # 1 sec
    baseline_means = np.mean(all_distances[:, :baseline], axis = 1).reshape(-1,1)
    baseline_SD = np.std(all_distances[:, :baseline], axis = 1).reshape(-1,1) #sd per trial, not just one value

    z_dists = (all_distances - baseline_means) / baseline_SD 

    # z_dists = np.array(all_distances)
    
    Zmean_dist = np.mean(z_dists, axis = 0) 
    Zsd_dist   = np.std(z_dists, axis = 0)
    
    threshold = 1.65
    
    p_vals = np.mean(z_dists <= threshold, axis=0)
    
    sig_mask = p_vals < 0.05
    
    min_cluster_bins = int(0.0 / bin_size) # abitrary,jercog does point wise (0.0s). change if you want consecutive thresholding

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
    
    plt.plot(time_centers, Zmean_dist, linewidth=2)
    plt.fill_between(time_centers,
                     Zmean_dist - Zsd_dist,
                     Zmean_dist + Zsd_dist,
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
                    y=max(Zmean_dist + Zsd_dist) * 1.05,
                    xmin=time_centers[start],
                    xmax=time_centers[end],
                    linewidth=3
                )
    
                start = sig_indices[i]
    
        # final segment
        end = sig_indices[-1]
    
        plt.hlines(
            y=max(Zmean_dist + Zsd_dist) * 1.05,
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
    
def perm_testing(pop_trials_A, pop_trials_B, gpfa_space, my_title, stim_trials = True, 
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
    
        latent_A = get_latent_gpfa(pop_A, gpfa_space)
        latent_B = get_latent_gpfa(pop_B, gpfa_space)
        
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
        
        dist = (dist - baseline_mean) / baseline_std # z-scored
    
        all_distances.append(dist)

    # for plotting
    mean_dist = np.mean(all_distances, axis = 0) 
    sd_dist   = np.std(all_distances, axis = 0)
    
    
    z_dists = np.array(all_distances)   

    threshold = 1.65
    
    p_vals = np.mean(z_dists <= threshold, axis=0)
    
    sig_mask = p_vals < 0.05
    
    min_cluster_bins = int(0.0 / bin_size) # abitrary,jercog does point wise (0.0s). change if you want consecutive thresholding

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


def get_latent_gpfa(pop, gpfa_model):
    """
    pop: (n_pop_trials, n_time, n_neurons)
    returns: (n_pop_trials, n_time, latent_dims)
    """
    trials = format_ipv_for_gpfa(pop)          # list of (n_neurons, n_time)
    latent = gpfa_model.transform(trials)      # list of (latent_dims, n_time)
    return np.array([x.T for x in latent])     # (n_pop_trials, n_time, latent_dims)


def plot_GPFA(latent_groups, main_title, dims=(0, 1, 2),
              t_pre=1.0, bin_size=1):

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    linewidth_single_trial = 0.5
    alpha_single_trial = 0.4
    linewidth_trial_average = 2.5

    ax.set_title(main_title)
    ax.set_xlabel(f'Dim {dims[0] + 1}')
    ax.set_ylabel(f'Dim {dims[1] + 1}')
    ax.set_zlabel(f'Dim {dims[2] + 1}')

    # Use matplotlib default color cycle
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, (latent, title) in enumerate(latent_groups):

        color = colors[i % len(colors)]

        # Single trials
        for trial in latent:
            ax.plot(
                trial[dims[0], :],
                trial[dims[1], :],
                trial[dims[2], :],
                lw=linewidth_single_trial,
                c=color,
                alpha=alpha_single_trial
            )

        # Trial average
        avg = np.mean(latent, axis=0)

        ax.plot(
            avg[dims[0], :],
            avg[dims[1], :],
            avg[dims[2], :],
            lw=linewidth_trial_average,
            c=color,
            label=f'{title} (mean)'
        )

        # Start marker
        ax.scatter(
            avg[dims[0], 0],
            avg[dims[1], 0],
            avg[dims[2], 0],
            c='k',
            marker='o'
        )

        # Event onset marker
        onset_idx = int(t_pre / bin_size)

        ax.scatter(
            avg[dims[0], onset_idx],
            avg[dims[1], onset_idx],
            avg[dims[2], onset_idx],
            c='k',
            marker='D'
        )

    ax.legend()
    plt.tight_layout()
    plt.show()
    
    
# -----------------------------
# GPFA — fit on all trials combined, project each type separately 
# -----------------------------
if not shuffle_trial_labels:
    MD_stim_path = rf'W:\Haak\Innate_defense\Data_analysis\22.35.02\Scotts_analysis\{region}_stimulus_aligned.pickle'
    MD_resp_path = rf'W:\Haak\Innate_defense\Data_analysis\22.35.02\Scotts_analysis\{region}_response_aligned.pickle'

else:
    MD_stim_path = rf'W:\Haak\Innate_defense\Data_analysis\22.35.02\Scotts_analysis\{region}_stimulus_aligned_shuffled_labels.pickle'
    MD_resp_path = rf'W:\Haak\Innate_defense\Data_analysis\22.35.02\Scotts_analysis\{region}_response_aligned_shuffled_labels.pickle'

    
# stim aligned

if shuffle_trial_labels:

    # Original lists
    freeze = all_neuron_trials_freeze
    nonfreeze = all_neuron_trials_nonfreeze
    
    # Store original sizes
    freeze_sizes = [arr.shape[0] for arr in freeze]
    nonfreeze_sizes = [arr.shape[0] for arr in nonfreeze]
    
    # Combine all trials
    all_trials = np.concatenate(freeze + nonfreeze, axis=0)
    
    # Shuffle trials
    np.random.shuffle(all_trials)
    
    # Split back
    idx = 0
    
    all_neuron_trials_freeze = []
    for size in freeze_sizes:
        all_neuron_trials_freeze.append(all_trials[idx:idx+size])
        idx += size
    
    all_neuron_trials_nonfreeze = []
    for size in nonfreeze_sizes:
        all_neuron_trials_nonfreeze.append(all_trials[idx:idx+size])
        idx += size

if analyze_stim:
    distances = np.zeros((n_repeats, 250)) # shouldnt hard code
    # fit model ONCE
    gpfa_model = GPFA(bin_size=bin_size * pq.s, x_dim=latent_dims)
    
    for repeat in list(range(n_repeats)):
        
        print(f'Repeat {repeat +1} of {n_repeats}')
        # split data 
        fit_neuron_trials_freeze, test_neuron_trials_freeze = split_trials_per_neuron(all_neuron_trials_freeze)    
        fit_neuron_trials_nonfreeze, test_neuron_trials_nonfreeze = split_trials_per_neuron(all_neuron_trials_nonfreeze)    
        
        #
        # Generate pseudo-populations
        #
        
        train_pop_freeze = generate_iPV_trajectories(
            fit_neuron_trials_freeze,
            n_pop_trials
        )
        
        train_pop_nonfreeze = generate_iPV_trajectories(
            fit_neuron_trials_nonfreeze,
            n_pop_trials
        )
       
        train_pop_freeze, keep_mask = remove_silent_neurons(train_pop_freeze)
        train_pop_nonfreeze = train_pop_nonfreeze[:, :, keep_mask]
         
        train_trials_combined = (
            format_ipv_for_gpfa(train_pop_freeze)
            + format_ipv_for_gpfa(train_pop_nonfreeze)
        )
        
        
        
        test_pop_freeze = generate_iPV_trajectories(
            test_neuron_trials_freeze,
            n_pop_trials
        )
        
        test_pop_nonfreeze = generate_iPV_trajectories(
            test_neuron_trials_nonfreeze,
            n_pop_trials
        )
        
        test_pop_freeze = test_pop_freeze[:, :, keep_mask]
        test_pop_nonfreeze = test_pop_nonfreeze[:, :, keep_mask]
        
        freeze_fmt = format_ipv_for_gpfa(test_pop_freeze)
        nonfreeze_fmt = format_ipv_for_gpfa(test_pop_nonfreeze)
        
     
        if find_latent_dims:
            
            test_trials_combined = (
                format_ipv_for_gpfa(test_pop_freeze)
                + format_ipv_for_gpfa(test_pop_nonfreeze)
            )
            
            # latent_dims = find_latent_dimensions(train_trials_combined, 'Stimulus aligned')
            
            latent_dims = find_latent_dimensions(
                train_trials_combined,
                test_trials_combined,
                'Stimulus aligned'
            )
            
        
        
        gpfa_model.fit(train_trials_combined)
        
        # Transform for each condition separately # TO DO rename variables     
        freeze_latent = gpfa_model.transform(freeze_fmt)
        nonfreeze_latent = gpfa_model.transform(nonfreeze_fmt)
    
        latent_freeze = np.array(freeze_latent)
        latent_nonfreeze = np.array(nonfreeze_latent)
            
        
        stim_test_groups = [
            (latent_freeze,'Freeze'),
            (latent_nonfreeze, 'Nonfreeze')
            ]
        
        # if project_train_data:
        #     train_freeze    = gpfa_model.transform(trials_freeze_fmt)     # list of (latent_dims, n_time)
        #     train_nonfreeze = gpfa_model.transform(trials_nonfreeze_fmt)
        #     stim_test_groups.append((train_freeze,'Freeze (Train)'))
        #     stim_test_groups.append((train_nonfreeze, 'Nonfreeze (Train)'))
        
        if repeat == 0:
            plot_GPFA(stim_test_groups, f'{region}: stimulus aligned')
    
    
        #Compute m distance
        distances[repeat] = mahalanobis_distance_trials(latent_freeze, latent_nonfreeze, repeat)
    
    # # Plot
    # plt.figure()
    # plt.plot(time_centers, dist)
    # plt.axvline(0, linestyle='--')
    # plt.title("GPFA trajectory distance")
    # plt.show()
    
    # # Permutation test
    # perm_testing(test_pop_freeze, test_pop_nonfreeze, gpfa_model, "GPFA distance (perm test)")
    distances = np.array(distances)
    

    # save file
    with open(MD_stim_path, 'wb') as handle:
        pickle.dump(distances, handle, protocol=pickle.HIGHEST_PROTOCOL)


    # loaded_distances = pd.read_pickle(MD_path)
    
    perm_testing2(distances, f'{region} stimulus aligned')

# # Z-score/extract max  +avg etc here
# baseline = 50 # 1 sec
# baseline_means = np.mean(distances[:, :baseline], axis = 1).reshape(-1,1)
# baseline_SD = np.std(distances[:, :baseline], axis = 1).reshape(-1,1)

# zBaseline = (distances - baseline_means) / baseline_SD

# mean_dist = np.mean(zBaseline, axis=0)
# sem_dist = sem(zBaseline, axis=0)

# plt.figure(figsize=(8,5))

# plt.plot(time_centers,
#          mean_dist,
#          color='black',
#          linewidth=2)

# plt.fill_between(time_centers,
#                  mean_dist - sem_dist,
#                  mean_dist + sem_dist,
#                  alpha=0.3)

# plt.axvline(0, linestyle='--')
# plt.axhline(0, linestyle='--')

# plt.xlabel('Time (s)')
# plt.ylabel('Mahalanobis distance')
# plt.title('Mean ± SEM distance')

# plt.tight_layout()
# plt.show()

#
# response aligned
#

if shuffle_trial_labels: 

    # Original lists
    freeze = all_neuron_trials_freeze_aligned
    nonfreeze = all_neuron_trials_nonfreeze_shuff
    
    # Store original sizes
    freeze_sizes = [arr.shape[0] for arr in freeze]
    nonfreeze_sizes = [arr.shape[0] for arr in nonfreeze]
    
    # Combine all trials
    all_trials = np.concatenate(freeze + nonfreeze, axis=0)
    
    # Shuffle trials
    np.random.shuffle(all_trials)
    
    # Split back
    idx = 0
    
    all_neuron_trials_freeze_aligned = []
    for size in freeze_sizes:
        all_neuron_trials_freeze_aligned.append(all_trials[idx:idx+size])
        idx += size
    
    all_neuron_trials_nonfreeze_shuff = []
    for size in nonfreeze_sizes:
        all_neuron_trials_nonfreeze_shuff.append(all_trials[idx:idx+size])
        idx += size 
        

    #ITI
    # Original lists
    freeze = all_neuron_trials_ITI_freeze
    nonfreeze = all_neuron_trials_ITI_nonfreeze_shuff
    
    # Store original sizes
    freeze_sizes = [arr.shape[0] for arr in freeze]
    nonfreeze_sizes = [arr.shape[0] for arr in nonfreeze]
    
    # Combine all trials
    all_trials = np.concatenate(freeze + nonfreeze, axis=0)
    
    # Shuffle trials
    np.random.shuffle(all_trials)
    
    # Split back
    idx = 0
    
    all_neuron_trials_ITI_freeze = []
    for size in freeze_sizes:
        all_neuron_trials_ITI_freeze.append(all_trials[idx:idx+size])
        idx += size
    
    all_neuron_trials_ITI_nonfreeze_shuff = []
    for size in nonfreeze_sizes:
        all_neuron_trials_ITI_nonfreeze_shuff.append(all_trials[idx:idx+size])
        idx += size    

if analyze_resp:  
    
    response_distances = np.zeros((n_repeats, 250)) # shouldnt hard code
    
    if not analyze_stim:
        # fit model ONCE
        gpfa_model = GPFA(bin_size=bin_size * pq.s, x_dim=latent_dims)
    
    for repeat in list(range(n_repeats)):
        
        print(f'Repeat {repeat +1} of {n_repeats}')

        # split data
        fit_neuron_trials_freeze_aln, test_neuron_trials_freeze_aln = split_trials_per_neuron(all_neuron_trials_freeze_aligned)    
        fit_neuron_trials_nonfreeze_aln, test_neuron_trials_nonfreeze_aln = split_trials_per_neuron(all_neuron_trials_nonfreeze_shuff) 
              

        train_pop_freeze_aln = generate_iPV_trajectories(
            fit_neuron_trials_freeze_aln,
            n_pop_trials
        )
        
        train_pop_nonfreeze_aln = generate_iPV_trajectories(
            fit_neuron_trials_nonfreeze_aln,
            n_pop_trials
        )
       
        train_pop_freeze_aln, keep_mask = remove_silent_neurons(train_pop_freeze_aln)
        train_pop_nonfreeze_aln = train_pop_nonfreeze_aln[:, :, keep_mask]
        
         
        align_trials_combined = (
            format_ipv_for_gpfa(train_pop_freeze_aln)
            + format_ipv_for_gpfa(train_pop_nonfreeze_aln)
        )
        
        test_pop_freeze_aln = generate_iPV_trajectories(
            test_neuron_trials_freeze_aln,
            n_pop_trials
        )
        
        test_pop_nonfreeze_aln = generate_iPV_trajectories(
            test_neuron_trials_nonfreeze_aln,
            n_pop_trials
        )
        
        test_pop_freeze_aln = test_pop_freeze_aln[:, :, keep_mask]
        test_pop_nonfreeze_aln = test_pop_nonfreeze_aln[:, :, keep_mask]
        
        freeze_fmt_aln = format_ipv_for_gpfa(test_pop_freeze_aln)
        nonfreeze_fmt_aln = format_ipv_for_gpfa(test_pop_nonfreeze_aln)
        
        
        if not include_ITI_analysis: 
        
            if find_latent_dims:
                latent_dims = find_latent_dimensions(align_trials_combined, 'Response aligned')
                               
            gpfa_model.fit(align_trials_combined)
            
            # Transform for each condition separately      
            freeze_latent_aln = gpfa_model.transform(freeze_fmt_aln)
            nonfreeze_latent_aln = gpfa_model.transform(nonfreeze_fmt_aln)
        
            latent_freeze_aln = np.array(freeze_latent_aln)
            latent_nonfreeze_aln = np.array(nonfreeze_latent_aln)
                
            
            resp_test_groups = [
                (latent_freeze_aln,'Freeze'),
                (latent_nonfreeze_aln, 'Nonfreeze')
                ]
            
            # if project_train_data:
            #     train_freeze    = gpfa_model.transform(trials_freeze_fmt)     # list of (latent_dims, n_time)
            #     train_nonfreeze = gpfa_model.transform(trials_nonfreeze_fmt)
            #     stim_test_groups.append((train_freeze,'Freeze (Train)'))
            #     stim_test_groups.append((train_nonfreeze, 'Nonfreeze (Train)'))
            
            if repeat == 0:
                plot_GPFA(resp_test_groups, f'{region}: response aligned')
        
        
            #Compute m distance
            response_distances[repeat] = mahalanobis_distance_trials(latent_freeze_aln, latent_nonfreeze_aln, repeat)        
        
        else:    # CODE BLOCK INCOMPLETE DO NOT RUN YET 
            print('this code isnt complete yet! skipping')
            continue
            #
            # ITI (response aligned)
            #
            # split data 
            fit_neuron_ITI_freeze, test_neuron_ITI_freeze = split_trials_per_neuron(all_neuron_trials_ITI_freeze)    
            fit_neuron_ITI_nonfreeze_shuff, test_neuron_ITI_nonfreeze_shuff = split_trials_per_neuron(all_neuron_trials_ITI_nonfreeze_shuff)    
            
            
            # create iPPV
            fit_pop_ITI_freeze = generate_iPV_trajectories(fit_neuron_ITI_freeze, n_pop_trials)
            fit_pop_ITI_nonfreeze = generate_iPV_trajectories(fit_neuron_ITI_nonfreeze_shuff, n_pop_trials)
            test_pop_ITI_freeze = generate_iPV_trajectories(all_neuron_trials_ITI_freeze, n_pop_trials)
            test_pop_ITI_nonfreeze = generate_iPV_trajectories(all_neuron_trials_ITI_nonfreeze_shuff, n_pop_trials)
            
            # Format and build combined pool for fitting (freeze + nonfreeze iPVs)
            fit_ITI_freeze_fmt    = format_ipv_for_gpfa(fit_pop_ITI_freeze)     # list of (n_neurons, n_time)
            fit_ITI_nonfreeze_fmt = format_ipv_for_gpfa(fit_pop_ITI_nonfreeze)
            
            # and for test
            test_ITI_freeze_fmt    = format_ipv_for_gpfa(test_pop_ITI_freeze)     
            test_ITI_nonfreeze_fmt = format_ipv_for_gpfa(test_pop_ITI_nonfreeze)
            all_fit_trials_combined      = align_trials_combined + fit_ITI_freeze_fmt + fit_ITI_nonfreeze_fmt
            
            if find_latent_dims:
                
                latent_dims = find_latent_dimensions(all_fit_trials_combined, 'Response aligned')
            
            gpfa_model.fit(all_fit_trials_combined)
        
        
            # Transform for each condition separately
            # latent_freeze    = gpfa_model.transform(align_trials_freeze_fmt_test)    
            # latent_nonfreeze = gpfa_model.transform(align_trials_nonfreeze_fmt_test)
            # latent_ITI_freeze    = gpfa_model.transform(test_ITI_freeze_fmt)    
            # latent_ITI_nonfreeze = gpfa_model.transform(test_ITI_nonfreeze_fmt)
        
            # Stack into arrays: (n_pop_trials, n_time, latent_dims) needed later?
            # latent_freeze_arr    = np.array([x.T for x in latent_freeze])
            # latent_nonfreeze_arr = np.array([x.T for x in latent_nonfreeze])
            
            # resp_test_groups = [
            #     (latent_freeze,'Freeze'),
            #     (latent_nonfreeze, 'Nonfreeze (shuffled onsets)'),
            #     (latent_ITI_freeze, 'ITI Freeze'),
            #     (latent_ITI_nonfreeze, 'ITI Nonfreeze (shuffled onsets)')
            #     ]
        
            # if project_train_data:
            #     train_freeze    = gpfa_model.transform(align_trials_freeze_fmt_fit)    
            #     train_nonfreeze = gpfa_model.transform(align_trials_nonfreeze_fmt_fit)
            #     train_ITI_freeze    = gpfa_model.transform(fit_ITI_freeze_fmt)
            #     train_ITI_nonfreeze = gpfa_model.transform(fit_ITI_nonfreeze_fmt)
                
            #     resp_test_groups.append((train_freeze,'Freeze (Train)'))
            #     resp_test_groups.append((train_nonfreeze, 'Nonfreeze (Train)'))
            #     resp_test_groups.append((train_ITI_freeze,'ITI Freeze (Train)'))
            #     resp_test_groups.append((train_ITI_nonfreeze, 'ITI Nonfreeze (Train)'))
        
            if repeat == 0:
                plot_GPFA(resp_test_groups, 'response aligned')
            
        

    # # Permutation test
    response_distances = np.array(response_distances)
    
    
    # save file
    with open(MD_resp_path, 'wb') as handle:
        pickle.dump(response_distances, handle, protocol=pickle.HIGHEST_PROTOCOL)

    perm_testing2(response_distances, f'{region} response aligned')

   
# # perm testing on mahalonibis distance
# # pre_stim_mean, pre_stim_std = perm_testing(population_freeze, population_nonfreeze, 
# #              my_title = f'Stimulus aligned, {region}', stim_trials = True)

# # perm_testing(population_freeze_aligned, population_nonfreeze_shuff,  
# #              my_title = f'Response aligned, {region}', stim_trials = False,
# #              baseline_means = pre_stim_mean, baseline_stds = pre_stim_std)

# # # check what baseline is first
# # perm_testing(population_freeze_aligned, population_ITI_freeze_aligned,  
# #              my_title = f'Response aligned incl ITI, {region}', stim_trials = True)

