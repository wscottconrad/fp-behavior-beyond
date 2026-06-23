# -*- coding: utf-8 -*-
"""
Created on Fri May  8 10:03:56 2026

creates neuron spike matrices aligend to stimulus and response onset, along with ITI 

WORKING: calculate neuron responsiveness with ZETA

@author: sconrad
"""

# analysis parameters
recalc_freeze_times = False
calc_zeta = True
calc_latenzy = True

remove_silent_neurons = False


plot_speed_heatmap = False
plot_avg_trial_and_ITI_freeze_speed = False
plot_avg_trial_and_ITI_nonfreeze_speed = False

# regions = ['superior colliculus'] 

regions = ['superior colliculus', 'periaqueductal gray'] 


import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt
from scipy.stats import sem
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import pandas as pd
import csv
import pickle
if calc_zeta:
    from zetapy import zetatest
    zeta_repetitions = 100 # could do 1000
    
if calc_latenzy:
    from latenzy import latenzy


# spike analysis initialize
t_pre = 1.0 # s
t_post = 4.0 # s
bin_size = 0.02  # 20 ms
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

alpha_sig = 0.05


data_path = r'W:\Haak\Innate_defense\Data_analysis\22.35.02'
session_list = [
    '98332_20240326_AP.mat',
    '98335_20240321_AP.mat',
    '100131_20240508_AP.mat',
    '100132_20240508_AP.mat',
    '100134_20240514_AP.mat'
]


freeze_file_path = r'W:\Conrad\Innate_approach\Materials_and_methods\Code'

save_path = r'W:\Haak\Innate_defense\Data_analysis\22.35.02\Scotts_analysis'

collected_spikes = {}
zeta_analysis = {}

def extract_peri_speed(signal, event_times_sec):
    event_idx = [ int(x) for x in (event_times_sec * fs)]
    segments = []

    for idx in event_idx:
        start = idx - int(s_pre * fs)
        end   = idx + int(s_post * fs)

        if start < 0 or end > len(signal):
            continue

        segments.append(signal[start:end])

    return np.array(segments)

def binned_firing_counts(spike_times, event_time, time_bins):
    aligned = spike_times - event_time
    counts, _ = np.histogram(aligned, bins=time_bins) # bins spikes within trial window
    # counts = np.sqrt(counts)  # variance-stabilizing transform # anscombe transform possible, but i havent seen it in papers
    return counts 

def neuron_trial_matrix(spike_times, event_onset): # for a single neuron, creates trial-based binned spike matrix
    mat = np.zeros((len(event_onset), len(time_bins) - 1))
    for i, ev in enumerate(event_onset):
        binned_counts = binned_firing_counts(spike_times, ev, time_bins)
        # rate = gaussian_filter1d(rate, smooth_sigma)
        mat[i] = binned_counts
    return mat

def run_zeta(spike_times, event_times, window=1.0):
    """
    spike_times : 1D np.array of spike timestamps in seconds
    event_times : 1D np.array of event timestamps in seconds
    window      : analysis window after event onset

    returns dictionary with zeta statistics
    """

    # require enough trials
    if len(event_times) < 3:
        return None

    try:
        dblZetaP, dZETA , dRate = zetatest(
            spike_times,
            event_times,
            # dblUseMaxDur=window,
            # boolReturnRate=False,
            resampling_number=1000
        )

        return {
            'p': dblZetaP,
            'zeta': dZETA['dblZETA'],
            'latency': dZETA.get('dblLatencyPeak', np.nan)
        }

    except Exception as e:
        print(f'ZETA failed: {e}')
        return None
    
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        # Returns string: Percentage (Actual Number)
        return '{p:.1f}%\n({v:d})'.format(p=pct, v=val)
    return my_autopct
   
def plot_wedge(sizes, labels, c, title):
    fig, ax = plt.subplots(figsize=(5,5))

    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=labels,
        colors = c,
        autopct = make_autopct(sizes),
        wedgeprops=dict(width=0.4)
    )
    
    ax.set_title(f'{title}')
    
    plt.show()    
        


for region in regions:
    # -----------------------------
    # Load & pool neurons across sessions
    # -----------------------------
     
    # stimulus aligned
    all_neuron_trials_freeze = []
    all_neuron_trials_nonfreeze = []
    
    # response aligned
    all_neuron_trials_freeze_aligned = []
    all_neuron_trials_nonfreeze_shuff = []
    
    # 'freeze' aligned during ITI
    all_neuron_trials_ITI_freeze = []
    all_neuron_trials_ITI_nonfreeze_shuff = []
    
    all_zeta_results = []
    
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
        
        # print(f'{sum(vecIsActive)}')
    
        clusters = sAP.sCluster
    
        is_good = np.array([c.KilosortGood for c in clusters]) # try no filtering
        low_contam = np.array([c.Contamination for c in clusters]) < 0.1
        areas = np.array([str(c.Area).strip().lower() for c in clusters])
        in_area = np.array([region in a for a in areas])
        # if remove_silent_neurons:
        #     not_silent = sum(np.array([len(c.SpikeTimes) for c in clusters]) > 10)
            
        #     plt.figure()
        #     plt.hist(not_silent, bins = "fd")
            
    
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
                movement = trial_velocity[:int((s_pre+stim_length)*fs)] >= pre_freeze_movement_thresh
                
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
            
        if recalc_freeze_times:
            pooled_freeze.append(freeze_onsets_sec)
            
        freeze_mask = np.array(freeze_mask)
        active_mask = np.array(active_mask)
    
        freeze_event_times = valid_stim_onset[freeze_mask]
        nonfreeze_event_times = valid_stim_onset[active_mask]
        nonresponse_event_times = valid_stim_onset[~freeze_mask & ~active_mask] 
        
        
        freeze_onsets_sec = np.array(freeze_onsets_sec)
        freeze_aligned_times = freeze_event_times + freeze_onsets_sec
        
        active_offsets = np.full(len(valid_stim_onset), np.nan)
    
        if len(freeze_onsets_sec) > 0 and len(nonfreeze_event_times) > 0 and recalc_freeze_times:
            shuffled_offsets = np.random.choice(
                freeze_onsets_sec,
                size=len(nonfreeze_event_times),
                replace=True
            )
            nonfreeze_shuffled_times = nonfreeze_event_times + shuffled_offsets
            active_offsets[active_mask] = shuffled_offsets
            
        elif len(nonfreeze_event_times) > 0 and not recalc_freeze_times:
            shuffled_offsets = np.random.choice(
                pooled_freeze,
                size=len(nonfreeze_event_times),
                replace=True
            )
            nonfreeze_shuffled_times = nonfreeze_event_times + shuffled_offsets
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
          
        # find 'freezing' during ITI    
        immobile_full = run_speed < freeze_thresh
    
        diff_full = np.diff(immobile_full.astype(int))
        iti_starts = np.where(diff_full == 1)[0] + 1
        iti_ends   = np.where(diff_full == -1)[0] + 1
        
        # edge cases
        if immobile_full[0]:
            iti_starts = np.insert(iti_starts, 0, 0)
        if immobile_full[-1]:
            iti_ends = np.append(iti_ends, len(immobile_full))
        
        exclude_sec_too_close = 5 # too close to trial start or end
        exclude_samples_too_close = int(exclude_sec_too_close * fs)
        
        exclude_mask = np.zeros(len(run_speed), dtype=bool)
        
        stim_on_idx = (stim_onset * fs).astype(int)
        stim_off_idx = (stim_off * fs).astype(int)
        
        for on, off in zip(stim_on_idx, stim_off_idx):
            start = max(0, on - exclude_samples_too_close)
            end   = min(len(run_speed), off + exclude_samples_too_close)
            exclude_mask[start:end] = True
            
        iti_freeze_onsets = []
    
        for s, e in zip(iti_starts, iti_ends):
            duration = e - s
        
            if duration < min_freeze_samples:
                continue
        
            # Reject if ANY part overlaps excluded mask
            if np.any(exclude_mask[s:e]):
                continue
        
            # require movement before freeze
            pre_freeze_samples = int(pre_freeze_window * fs)
            pre_start = max(0, s - pre_freeze_samples)
            pre_vel = run_speed[pre_start:s]
        
            if len(pre_vel) == 0 or np.max(pre_vel) < pre_freeze_movement_thresh:
                continue
        
            iti_freeze_onsets.append(s / fs)
            
        iti_freeze_onsets = np.array(iti_freeze_onsets)
        
        if plot_avg_trial_and_ITI_freeze_speed:
            
            # Trial-aligned freeze (you already computed this)
            trial_freeze_speed = extract_peri_speed(run_speed, freeze_aligned_times)
            
            # ITI freeze (from new code)
            iti_freeze_speed = extract_peri_speed(run_speed, iti_freeze_onsets)
            
            trial_speed_mean = np.mean(trial_freeze_speed, axis = 0)
            trial_speed_sem = sem(trial_freeze_speed)
            
            iti_speed_mean = np.mean(iti_freeze_speed, axis=0)
            iti_speed_sem = sem(iti_freeze_speed)
        
            ts = np.linspace(-s_pre, s_post, win_len)
            
            plt.figure(figsize=(6, 4))
             
            # Trial freezes
            if trial_speed_mean is not None:
                plt.plot(ts, trial_speed_mean, label='Trial freeze')
                plt.fill_between(
                    ts,
                    trial_speed_mean - trial_speed_sem,
                    trial_speed_mean + trial_speed_sem,
                    alpha=0.3
                )
            
            # ITI freezes
            if iti_speed_mean is not None:
                plt.plot(ts, iti_speed_mean, linestyle='--', label='ITI freeze')
                plt.fill_between(
                    ts,
                    iti_speed_mean - iti_speed_sem,
                    iti_speed_mean + iti_speed_sem,
                    alpha=0.3
                )
            
            plt.axvline(0, linestyle='--', linewidth=1)
            plt.xlabel('Time from freeze onset (s)')
            plt.ylabel('Run speed (m/s)')
            plt.title(f'{sess[:6]}: Speed aligned to freeze onset')
            plt.legend()
            plt.tight_layout()
            plt.show()
            
        # -----------------------------------------------
        # Detect movement-only periods during ITI (no freeze)
        # -----------------------------------------------
        iti_active_onsets = []
        
        # Build a movement mask: above threshold AND not excluded AND not immobile
        moving_full = (run_speed >= pre_freeze_movement_thresh) & ~immobile_full & ~exclude_mask
        
        # running bouts
        diff_move_full = np.diff(moving_full.astype(int))
        iti_move_starts = np.where(diff_move_full == 1)[0] + 1
        iti_move_ends   = np.where(diff_move_full == -1)[0] + 1
        
        # edge cases
        if moving_full[0]:
            iti_move_starts = np.insert(iti_move_starts, 0, 0)
        if moving_full[-1]:
            iti_move_ends = np.append(iti_move_ends, len(moving_full))
        
        # min_move_samples = int(pre_freeze_window * fs) # 1 sec
        
        for s, e in zip(iti_move_starts, iti_move_ends):
            duration = e - s
        
            if duration < min_move_samples:
                continue
        
            # Confirm no freeze occurs after movement onset within a reasonable window
            # (mirrors the trial-based active_mask logic — movement but no subsequent freeze)
            post_window_end = min(len(run_speed), e + int(pre_freeze_window * fs))
            post_velocity = run_speed[e:post_window_end]
        
            if len(post_velocity) > 0 and np.any(post_velocity < freeze_thresh):
                # Animal went immobile shortly after — skip, could be an ITI freeze
                continue
        
            iti_active_onsets.append(s / fs)
        
        iti_active_onsets = np.array(iti_active_onsets)
        
        # Apply shuffled offsets to create pseudo-aligned times
        if len(iti_active_onsets) > 0 and len(freeze_onsets_sec) > 0 and recalc_freeze_times:
            iti_shuffled_offsets = np.random.choice(
                freeze_onsets_sec,
                size=len(iti_active_onsets),
                replace=True
            )
        elif len(iti_active_onsets) > 0 and not recalc_freeze_times:
            iti_shuffled_offsets = np.random.choice(
                pooled_freeze,
                size=len(iti_active_onsets),
                replace=True
            )
        else:
            iti_shuffled_offsets = np.array([])
        
        iti_nonfreeze_shuffled_times = (
            iti_active_onsets + iti_shuffled_offsets
            if len(iti_shuffled_offsets) > 0
            else np.array([])
        )
        
        # verify
        if plot_avg_trial_and_ITI_nonfreeze_speed:
            
            # Trial-aligned nonfreeze 
            trial_nonfreeze_speed = extract_peri_speed(run_speed, nonfreeze_shuffled_times)
            
            # ITI nonfreeze (from new code)
            iti_nonfreeze_speed = extract_peri_speed(run_speed, iti_nonfreeze_shuffled_times)
            
            trial_speed_mean = np.mean(trial_nonfreeze_speed, axis = 0)
            trial_speed_sem = sem(trial_nonfreeze_speed)
            
            iti_speed_mean = np.mean(iti_nonfreeze_speed, axis=0)
            iti_speed_sem = sem(iti_nonfreeze_speed)
        
            ts = np.linspace(-s_pre, s_post, win_len)
            
            plt.figure(figsize=(6, 4))
             
            # Trial nonfreezes
            if trial_speed_mean is not None:
                plt.plot(ts, trial_speed_mean, label='Trial nonfreeze')
                plt.fill_between(
                    ts,
                    trial_speed_mean - trial_speed_sem,
                    trial_speed_mean + trial_speed_sem,
                    alpha=0.3
                )
            
            # ITI freezes
            if iti_speed_mean is not None:
                plt.plot(ts, iti_speed_mean, linestyle='--', label='ITI nonfreeze')
                plt.fill_between(
                    ts,
                    iti_speed_mean - iti_speed_sem,
                    iti_speed_mean + iti_speed_sem,
                    alpha=0.3
                )
            
            plt.axvline(0, linestyle='--', linewidth=1)
            plt.ylim(0,0.3)
            plt.xlabel('Time from shuffled onset (s)')
            plt.ylabel('Run speed (m/s)')
            plt.title(f'{sess[:6]}: Speed aligned to shuffled nonfreeze onset')
            plt.legend()
            plt.tight_layout()
            plt.show()
    
        print(f'Number of recorded neurons this session: {sum(idx_keep)}')
        
        # zeta_results = []
        counter = 0 
        
        for i, c in enumerate(clusters):
        
            if not idx_keep[i]:
                continue
            
                
            # -------------------------------
            # Stimulus-aligned 
            # -------------------------------
            trials_freeze = neuron_trial_matrix(
                c.SpikeTimes,
                freeze_event_times
            )
        
            trials_nonfreeze = neuron_trial_matrix(
                c.SpikeTimes,
                nonfreeze_event_times
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
            
            # ITI freeze aligned
            trials_iti_freeze = neuron_trial_matrix(
                c.SpikeTimes,
                iti_freeze_onsets
            )
            
            # ITI non-freeze (movement + shuffled offset) aligned
            trials_iti_nonfreeze_shuff = neuron_trial_matrix(
                c.SpikeTimes,
                iti_nonfreeze_shuffled_times
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
                
                all_neuron_trials_ITI_freeze.append(trials_iti_freeze)
                all_neuron_trials_ITI_nonfreeze_shuff.append(trials_iti_nonfreeze_shuff)  
                
            
            # --------------------------------
            # ZETA tests
            # --------------------------------
            if calc_zeta:
                # all stimulus, stim aligned
                if len(valid_stim_onset) > 2:
                    zeta_stim_all, dZETA , dRate = zetatest(
                        c.SpikeTimes,
                        valid_stim_onset,
                        resampling_number= zeta_repetitions
                    )
                else:
                    zeta_stim_all = np.nan
                    
                # all freezing, stim aligned
                if len(freeze_event_times) > 2:
                    zeta_stim_freeze, dZETA , dRate = zetatest(
                        c.SpikeTimes,
                        freeze_event_times,
                        resampling_number= zeta_repetitions
                    )
                else: 
                    zeta_stim_freeze = np.nan
                    
                # all nonfreezing, stim aligned
                if len(nonfreeze_event_times) > 2:
                    zeta_stim_nonfreeze, dZETA , dRate = zetatest(
                        c.SpikeTimes,
                        nonfreeze_event_times,
                        resampling_number= zeta_repetitions
                    )
                else:
                    zeta_stim_nonfreeze = np.nan
                
                # freezing aligned, stim
                if len(freeze_aligned_times) > 2:
                    zeta_freeze_aligned_stim, dZETA , dRate = zetatest(
                        c.SpikeTimes,
                        freeze_aligned_times,
                        resampling_number= zeta_repetitions
                    )
                else:
                    zeta_freeze_aligned_stim = np.nan
                    
                # freezing aligned, no stim
                if len(iti_freeze_onsets) > 2:
                    zeta_freeze_aligned_ITI, dZETA , dRate = zetatest(
                        c.SpikeTimes,
                        iti_freeze_onsets,
                        resampling_number= zeta_repetitions
                    )
                else: 
                    zeta_freeze_aligned_ITI = np.nan
                
                all_zeta_results.append({
                    'session': sess,
                    'cluster_idx': i,
                    'region': region,
                    'depth' : c.Depth,
                    
                    'stim_all': zeta_stim_all,
                    'stim_freeze': zeta_stim_freeze,
                    'stim_nonfreeze': zeta_stim_nonfreeze,
                    'freeze_aligned_stim': zeta_freeze_aligned_stim,
                    'freeze_aligned_ITI': zeta_freeze_aligned_ITI
                })
                
            if calc_latenzy:
                latenzy_inputs = {
                    'all_stim_latency': valid_stim_onset,
                    'stim_freeze_latency': freeze_event_times,
                    'stim_nonfreeze_latency': nonfreeze_event_times,
                    'freeze_aligned_latency': freeze_aligned_times,
                    'iti_freeze_latency': iti_freeze_onsets,
                }
                
                
                for result_key, event_times in latenzy_inputs.items():
                    all_zeta_results[counter][result_key] = latenzy(
                        c.SpikeTimes,
                        event_times,
                        [-1, 4]
                    )
        
       
            
            # fraction analysis skip for now
            # if False:
            #     sig_stim_all = df_zeta['stim_all_p'] < alpha_sig
            #     sig_stim_freeze = df_zeta['stim_freeze_p'] < alpha_sig
            #     sig_stim_nonfreeze = df_zeta['stim_nonfreeze_p'] < alpha_sig
            #     sig_freeze_aligned_stim = df_zeta['freeze_aligned_stim_p'] < alpha_sig
            #     sig_freeze_aligned_ITI = df_zeta['freeze_aligned_ITI_p'] < alpha_sig
                
            #     # any modulation
            #     sig_any = (
            #         sig_stim_all |
            #         sig_stim_freeze |
            #         sig_stim_nonfreeze |
            #         sig_freeze_aligned_stim |
            #         sig_freeze_aligned_ITI
            #     )
                
            #     # no response
            #     sig_none = ~sig_any
                
            #     # more categories
            #     stim_freeze_only = sig_stim_freeze & ~sig_stim_nonfreeze
            #     stim_nonfreeze_only = sig_stim_nonfreeze & ~sig_stim_freeze
            #     stim_both = sig_stim_freeze & sig_stim_nonfreeze
                
            #     freeze_aligned_stim_only = sig_freeze_aligned_stim & ~sig_freeze_aligned_ITI
            #     freeze_aligned_ITI_only = sig_freeze_aligned_ITI & ~sig_freeze_aligned_stim
            #     freeze_both = sig_freeze_aligned_stim & sig_freeze_aligned_ITI
                
            #     true_mixed = sig_stim_nonfreeze + freeze_aligned_ITI_only
                
            #     # -----------------------------------
            #     # fractions
            #     # -----------------------------------
                
            #     # n_cells = len(df_zeta)
                
            #     fractions = {
            #         'Stim all': sig_stim_all.mean(),
            #         'Stim freeze': sig_stim_freeze.mean(),
            #         'Stim nonfreeze': sig_stim_nonfreeze.mean(),
            #         'Stim freeze aligned': sig_freeze_aligned_stim.mean(),
            #         'ITI freeze aligned': sig_freeze_aligned_ITI.mean(),
            #         'No response': sig_none.mean()
            #     }
                
            #     counts = {
            #         'Stim all': sig_stim_all.sum(),
            #         'Stim freeze': sig_stim_freeze.sum(),
            #         'Stim nonfreeze': sig_stim_nonfreeze.sum(),
            #         'Stim freeze aligned': sig_freeze_aligned_stim.sum(),
            #         'ITI freeze aligned': sig_freeze_aligned_ITI.sum(),
            #         'No response': sig_none.sum()
            #     }
            
    
            # plt.figure(figsize=(7,4))
            
            # labels = list(fractions.keys())
            # values = list(fractions.values())
            
            # bars = plt.bar(labels, values)
            
            # for bar, label in zip(bars, labels):
            
            #     count = counts[label]
            
            #     plt.text(
            #         bar.get_x() + bar.get_width()/2,
            #         bar.get_height() + 0.01,
            #         f'n={count}',
            #         ha='center'
            #     )
            
            # plt.ylabel('Fraction of neurons')
            # plt.ylim(0,1)
            # plt.title(f'Overview of modulated neurons in {region}')
            
            # plt.tight_layout()
            # plt.show()
            
            # ###
            # # pi charts
            # ###
            
            # colors = [
            #     '#2fd7d9',
            #     '#e6e21e',
            #     '#832fd9'
            #     ]
            # # 
            # sizes = [
            #     stim_freeze_only.sum(),
            #     stim_nonfreeze_only.sum(),
            #     stim_both.sum()
            # ]
            
            # labels = [
            #     'Stim freeze only',
            #     'Stim nonfreeze only',
            #     'Both'
            # ]
            
        
            # plot_wedge(sizes, labels, colors, 'Stimulus-responsive neurons')
            
            
            # sizes = [
            #     freeze_aligned_stim_only.sum(),
            #     freeze_aligned_ITI_only.sum(),
            #     freeze_both.sum()
            # ]
            
            # labels = [
            #     'Freeze (stim) only',
            #     'Freeze (ITI) only',
            #     'Both'
            # ]
            
        
            # plot_wedge(sizes, labels, colors, 'Freeze-responsive neurons')
            
            # sizes = [
            #     stim_nonfreeze_only.sum(),
            #     freeze_aligned_ITI_only.sum(),
            #     true_mixed.sum()
            # ]
            
            # labels = [
            #     'Stim nonfreeze only',
            #     'Freeze (ITI) only',
            #     'Both'
            # ]
            
            # plot_wedge(sizes, labels, colors, 'Mixed neurons')
            
        
            counter += 1
        
        print(f"Freeze trials: {sum(freeze_mask)}")
        print(f"Nonfreeze trials: {sum(active_mask)}")
        
    
        
    if recalc_freeze_times:
        flattened_pooled_freeze = [x 
                                   for xs in pooled_freeze
                                   for x in xs]
        df = pd.DataFrame(flattened_pooled_freeze)
        df.to_csv(f"{freeze_file_path}\\freeze_times.csv", header=False, index=False)
    
    if calc_zeta:
        # collected_spikes[region] = {
        #     'zeta_results': all_zeta_results
        # }
        
        rows = []

        for z in all_zeta_results:
        
            rows.append({
                'session': z['session'],
                'cluster': z['cluster_idx'],
        
                'stim_all_p':
                    z['stim_all']
                    if z['stim_all'] else np.nan,
        
                'stim_freeze_p':
                    z['stim_freeze']
                    if z['stim_freeze'] else np.nan,
        
                'stim_nonfreeze_p':
                    z['stim_nonfreeze']
                    if z['stim_nonfreeze'] else np.nan,
                    
                'freeze_aligned_stim_p':
                    z['freeze_aligned_stim']
                    if z['freeze_aligned_stim'] else np.nan,
        
                'freeze_aligned_ITI_p':
                    z['freeze_aligned_ITI']
                    if z['freeze_aligned_ITI'] else np.nan
            })
        
        df_zeta = pd.DataFrame(rows)
        
    collected_spikes[region] = {
            'all_neuron_trials_freeze': all_neuron_trials_freeze,
            'all_neuron_trials_nonfreeze': all_neuron_trials_nonfreeze,
            'all_neuron_trials_freeze_aligned': all_neuron_trials_freeze_aligned,
            'all_neuron_trials_nonfreeze_shuff': all_neuron_trials_nonfreeze_shuff,
            'all_neuron_trials_ITI_freeze': all_neuron_trials_ITI_freeze,
            'all_neuron_trials_ITI_nonfreeze_shuff': all_neuron_trials_ITI_nonfreeze_shuff,
            'zeta_results': df_zeta
            }

    
        

# if not remove_silent_neurons:
#     # Store data (serialize)
#     with open(save_path + '\collected_spikes.pickle', 'wb') as handle:
#         pickle.dump(collected_spikes, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
# else:
#     with open(save_path + '\collected_spikes_silent_neurons_filtered.pickle', 'wb') as handle:
#         pickle.dump(collected_spikes, handle, protocol=pickle.HIGHEST_PROTOCOL)
            