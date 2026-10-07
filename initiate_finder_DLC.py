# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 15:09:59 2025

find approach initiation times 

@author: conrad
"""
def initiate_finder_DLC(ttlFile, eventTS_RWD, r_log, fIdx, l, 
                    stim_dur, sr, filter_mvmnt, plotknee, 
                    cleaner_baseline, exclude_misaligned_trials):
    
    import pickle
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from kneed import KneeLocator
    from scipy.signal import argrelextrema
    from ITI_approach_finder import ITI_approach_finder
    import math
    from scipy.ndimage import median_filter

    
    def scotts_knee(x, y, initiation, title, line_color, ax=None):

        if ax is None:
            fig, ax = plt.subplots()
    
        y = y[::-1]  
        
        ax.plot(x, y, color = line_color)
        ax.axvline(initiation, color='r', linestyle='--')
        ax.set_ylim(min(0, y.min()), max(0, y.max()))
        ax.set_title(title)
        
    def detect_knees(x, y, trim, curve="concave", direction="increasing"):
        x = x[:-trim] # initiate will never be this fast
        y = y[trim:]
        
        knee_index = 0
        
        data = y.copy()
        
        shift_counter = 0
    
        while len(data) > 2:  # Continue until there's enough data to find a knee
            kl = KneeLocator(x[:len(data)], data, curve=curve, direction=direction)
            if kl.knee is not None:
                
                # Remove the data before the knee
                if kl.knee == 0:
                    knee_index = 1
                    shift_counter += 1
                else:
                    # kl.plot_knee()
                    # knees = kl.knee
                    knee_index = kl.knee
                    
                    break 
                
                data = data[knee_index:]  # Slice the data after the knee
                x = x[:-knee_index]  # Slice the x-values as well
            else:
                break  # No more knees detected
                
        return knee_index
      
    # def initiate_finder(trial_starts, trial_ends, snout_distance, initiation_array, plotknee, trim_array, trial_type):
    #     for i, (start, end) in enumerate(zip(trial_starts, trial_ends), 1):
    #         if end-start > 600: # in case of IR trials lasting longer than 20 seconds
    #             end = start + 600
    #         y = snout_distance[start:end][::-1]
    #         x = np.arange(len(y))  # trial-relative time
            
    #         trim = trim_array[i-1]
                        
    #         local_min = argrelextrema(y, np.less)
    #         local_min = local_min[0]
            
    #         if local_min.size == 0: # case where animal always approaches towards laser
    #             initiation = 0
                
    #             if plotknee == True:
                    
    #                 title = f"{trial_type} approach only"
    #                 scotts_knee(x, y, initiation, title)
                     
    #         else:    
    #             knee = detect_knees(x, y, trim)
    #             initiation = x[-1] - (knee + trim) # flip it back
                
    #             if plotknee == True:
    #                 title = f"{trial_type} Approach initiation" 
    #                 scotts_knee(x, y, initiation, title)
                
    #         initiation_array[i-1] = int(initiation)
            
    #     return initiation_array
    

    def initiate_finder(trial_starts, trial_ends, snout_distance,
                        initiation_array, plotknee, trim_array, trial_type,
                        synced_trials = []):
    
        if plotknee:
            n_trials = len(trial_starts)
            ncols = 4
            nrows = math.ceil(n_trials / ncols)
    
            fig, axs = plt.subplots(nrows, ncols,
                                    figsize=(4*ncols, 3*nrows),
                                    squeeze=False)
            axs = axs.ravel()
    
        for i, (start, end) in enumerate(zip(trial_starts, trial_ends), 1):

            # if end - start > 600:
            #     end = start + 600
    
            y = snout_distance[start:end][::-1]
            x = np.arange(len(y))
    
            trim = trim_array[i-1]
    
            local_min = argrelextrema(y, np.less)[0]
            
            if i-1 not in synced_trials and trial_type == 'prey':
                line_color = 'red'
            else:
                line_color = 'black'
            
            if local_min.size == 0:
                initiation = 0
    
                if plotknee:
                    title = f"{trial_type} approach only"
                    scotts_knee(x, y, initiation, title, line_color, ax=axs[i-1])
    
            elif len(x) < 599 and np.any(y < 75): # 7.5 cm is liberal cutoff for approach
                knee = detect_knees(x, y, trim)
                initiation = x[-1] - (knee + trim)
    
                if plotknee:
                    title = f"{trial_type} Approach initiation"
                    scotts_knee(x, y, initiation, title, line_color, ax=axs[i-1])
                    
            else:
                initiation = x[-1] # i think this is end of trial?
                if plotknee:
                    title = f"{trial_type} no initiation or chase trial"
                    scotts_knee(x, y, initiation, title, line_color, ax=axs[i-1])
    
            initiation_array[i-1] = int(initiation)
    
        if plotknee:
            # Remove unused axes
            for ax in axs[len(trial_starts):]:
                fig.delaxes(ax)
    
            fig.tight_layout()
            plt.show()
    
        return initiation_array
    
    # Trace window and settings
    pre = 5  # 5 seconds before ttl
    post = 25  # 25 seconds after
    before = pre * sr
    after = post * sr
    trial_length = before + after
    filter_value = 3 # 5 didnt work so well
    
    mismatch_threshold = int(90*sr) # seconds to frame. some drift occurs between pi camera and rwd. how much do you want to tolerate? set from 1 second to 100ms 17/9/26
    
    # apr_times = [np.nan]
    dlc_savePath = 'W:\\Conrad\\Innate_approach\\Data_analysis\\24.35.01\\DLC'
    idn = str(int(r_log['ID'][l])) 
    date = str(r_log['Date'][l]).replace('_', '')
         
    
    if idn == '118583' and date == '20260514':
        switch_case = True
    else:
        switch_case = False 
        
    # approach_time_path = f"{dlc_savePath}\\approach_times_since_trial_start.pkl"
    # pd.to_pickle(apr_times, approach_time_path)

    # with open('W:\\Conrad\\Innate_approach\\Data_analysis\\24.35.01\\DLC\\prey_snout_distance_combined_DLC.pkl', 'rb') as f:
        # prey_snout_distance = pickle.load(f)
        
    with open(f"{dlc_savePath}\\{r_log['Date'][l]}{idn} DLC.pkl", 'rb') as f:
        all_data = pickle.load(f)
        prey_snout_distance = all_data['data'][3]
        IR_snout_distance = all_data['data'][4]
        snout_speed = all_data['data'][5]
        hrC_speed = all_data['data'][6] # head ring, center
        tail_speed = all_data['data'][7]
        # misaligned_trials = all_data['data'][8] # could use this to filter out trials with poor alignment 
        out_of_frame = all_data['data'][10]

    
    if filter_mvmnt:
        prey_snout_distance = median_filter(prey_snout_distance, size = (filter_value))
        IR_snout_distance = median_filter(IR_snout_distance, size = (filter_value))
        snout_speed = median_filter(snout_speed, size = (filter_value))
        hrC_speed = median_filter(hrC_speed, size = (filter_value)) 
        tail_speed = median_filter(tail_speed, size = (filter_value))
    
    # Detect where prey trials start and end
    valid = ~np.isnan(prey_snout_distance) # boolean mask
    prey_trial_starts = np.where(np.diff(valid.astype(int)) == 1)[0] + 1
    if l == 14:
        prey_trial_starts = prey_trial_starts[1:] # error where trial occured but not ttl registered by rwd. should save this trial somehow? 
        
    prey_trial_ends = np.where(np.diff(valid.astype(int)) == -1)[0] + 1
    if l == 14:
        prey_trial_ends = prey_trial_ends[1:] # error where trial occured but not ttl registered by rwd. should save this trial somehow? 
    
    # Detect where IR trials start and end
    valid = ~np.isnan(IR_snout_distance) # boolean mask
    IR_trial_starts = np.where(np.diff(valid.astype(int)) == 1)[0] + 1
    IR_trial_ends   = np.where(np.diff(valid.astype(int)) == -1)[0] + 1
    
    
    
    # convert from pandas to array
    if not filter_mvmnt:
        prey_snout_distance = prey_snout_distance.to_numpy()
        IR_snout_distance = IR_snout_distance.to_numpy()
    
    
    ttl = pd.read_csv(ttlFile)
     
# %% checks for drift between ttl pulses and camera log
# this code doesnt affect anything downstream, its just a check

    if np.any(ttl['Event'] == 'sync'):
        ttl_syncs = ttl.loc[ttl['Event']=='sync'].reset_index(drop = True)
        ttl_syncs['Time'] = (pd.to_datetime(ttl_syncs['DateTime']) - pd.to_datetime(ttl_syncs['DateTime'][0])).dt.total_seconds() + (ttl_syncs['Milliseconds'] / 1000) - ttl_syncs['Milliseconds'][0] / 1000
        ttl_syncs = np.diff(ttl_syncs['Time'][::2])
        
        
        fp_filePath = 'W:\\Conrad\\Innate_approach\\Data_collection\\24.35.01\\'
       
        cam_ttl_file_path = f"{fp_filePath}{idn}\\{idn}_{date}_001\\{idn}_{date}_001_pioverhead_triggers.csv"
        cam_ttl = pd.read_csv(cam_ttl_file_path)
        if 'input type' in cam_ttl.columns:
            cam_ttl = cam_ttl[cam_ttl['input type'] == 'recieved']            
        cam_ttl_diff = np.diff(cam_ttl[cam_ttl['ttl source']=='sync_ttl']['time elapsed (s)'])
        sync_check = (ttl_syncs - cam_ttl_diff)*1000
        
        # print(f"The difference between start and end of habituation timing is {round(sync_check[0])} ms")
        # print(f"The difference between end of habituation and end session is timing is {round(sync_check[1])} ms")
        
        if abs(sync_check[0]) > 2000 or abs(sync_check[1]) > 2000:
            print('\n large (>2 second) drift detected! inspect session further\n')
        
    
# %%

    # if ttl['Event'][0] == 'Received trigger':
    #     #subtracts neurotar ttl input from first laser ttl output
        
    #     ####
    #     ### CHANGE THIS SO IT TAKES LAST POSSIBLE NEUROTAR SIGNAL
    #     ttl['Time'] = (pd.to_datetime(ttl['DateTime']) - pd.to_datetime(ttl['DateTime'][0])).dt.total_seconds() + (ttl['Milliseconds'] / 1000) - ttl['Milliseconds'][0] / 1000
    #     ###
        
    #     ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True) # added because received trigger sometimes added twice
    #     ttl = ttl.loc[ttl['Event'] !='sync'].reset_index(drop = True)
    #     ttl['Time'] = ttl['Time'] - stim_dur/sr
        
    # elif ttl['Event'][0] == 'sync':
    #     ttl['Time'] = (pd.to_datetime(ttl['DateTime']) - pd.to_datetime(ttl['DateTime'][0])).dt.total_seconds() + (ttl['Milliseconds'] / 1000) - ttl['Milliseconds'][0] / 1000
    #     ttl = ttl.loc[ttl['Event'] !='sync'].reset_index(drop = True)
    #     ttl['Time'] = ttl['Time'] - stim_dur/sr
    # else:
    #     print('No sync pulse detected at start, something has gone wrong')
        

    # if fIdx.size == 1:
    #     expanded_fIdx = [fIdx*3 +offset for offset in range(3)]
    # else: 
    #     expanded_fIdx = [i + offset for i in fIdx*3 for offset in range(3)]
    
    # ttl = ttl.drop(expanded_fIdx, axis = 0).reset_index(drop = True)    
        
    # ttl = ttl['Time'][::3]
    
    if '2025' not in date: # all recieved triggers after 2025 (building move) is noise
        ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True)
    
    if ttl['Event'][0] == 'Received trigger':
        
        #subtracts neurotar ttl input from first laser ttl output
        ttl['Time'] = (pd.to_datetime(ttl['DateTime']) - pd.to_datetime(ttl['DateTime'][0])).dt.total_seconds() + (ttl['Milliseconds'] / 1000) - ttl['Milliseconds'][0] / 1000
        
        ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True) # added because received trigger sometimes added twice
        ttl = ttl.loc[ttl['Event'] !='sync'].reset_index(drop = True)
        ttl['Time'] = ttl['Time'] - stim_dur/sr
    
    elif ttl['Event'][0] == 'sync':
        ttl['Time'] = (pd.to_datetime(ttl['DateTime']) - pd.to_datetime(ttl['DateTime'][0])).dt.total_seconds() + (ttl['Milliseconds'] / 1000) - ttl['Milliseconds'][0] / 1000
        # ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True) # added because received trigger sometimes added twice
        ttl = ttl.loc[ttl['Event'] !='sync'].reset_index(drop = True)
        ttl['Time'] = ttl['Time'] - stim_dur/sr
    else:
        print('No sync pulse detected at start, something has gone wrong')
        
    fIdx = np.array([])    # code for removing failed laser presentations
    if isinstance(r_log['Failed'][l], int):
        fIdx = np.array(int(r_log['Failed'][l]))
        
    elif isinstance(r_log['Failed'][l], float) and  np.isnan(r_log['Failed'][l]) == False:
        fIdx = np.array(list(map(int, str(r_log['Failed'][l]).split('.'))))
        if fIdx[1] == 0:
            fIdx = np.delete(fIdx, 1)
        
    elif isinstance(r_log['Failed'][l], str) and len(r_log['Failed'][l]) == 1:
        fIdx = np.array(int(r_log['Failed'][l]))
        
    elif isinstance(r_log['Failed'][l], str) and len(r_log['Failed'][l]) > 1:
        fIdx = np.array(list(map(int, str(r_log['Failed'][l]).split('.'))))
     
        
    if ((sum(ttl['Event'] == 'optogenetics') != sum(ttl['Event'] == 'prey')) and
        (sum(ttl['Event'] == 'optogenetics') > 0)):
        
        print('\nDifferent prey key triggers were pressend in session, you need to write custom code to fix\n')
        # no cases detected from 2025 - 2026
        
    if ('optogenetics' == ttl['Event']).any(): # does not cover all situations, for example, a session where both 'b' and 'p' were pressed 
        expand = 3
        
    else:
        expand = 2
        
        
    if fIdx.size == 1:
        expanded_fIdx = [fIdx*expand +offset for offset in range(expand)]
    else: 
        expanded_fIdx = [i + offset for i in fIdx*expand for offset in range(expand)]
    
    ttl = ttl.drop(expanded_fIdx, axis = 0).reset_index(drop = True)    
        
    ttl = ttl['Time'][::expand]  
    


# %% this section detects misalignments between control pi ttl outputs and either DLC or RWD recieved inputs
# it then applies a correction to try to align the two

    drift = (ttl.iloc[-1] - ttl.iloc[0])*30 - (eventTS_RWD[-1] - eventTS_RWD[0])
    drift_ratio = drift/(eventTS_RWD[-1] - eventTS_RWD[0])
    
    if abs(drift) > 3: # 100 ms difference 
        print(f"The frame difference between last and first ttl pulse sent vs recieved is {drift}\n")
       
    # apply drift correction 
    eventTS_RWD = np.array(eventTS_RWD)
    for x, value in enumerate(eventTS_RWD):
        eventTS_RWD[x] = int(value + value*drift_ratio)
        
    # syncing between what i see on camera + ttl file
    ttl_file_diff = np.array(np.diff(ttl)*30, dtype = int)
    DLC_diff = np.diff(prey_trial_starts)
    ttl_check_pis = ttl_file_diff - DLC_diff
    
    ttl_check_pi_rwd = ttl_file_diff - np.diff(eventTS_RWD)
    if any(abs(ttl_check_pi_rwd) > 3): #90 ms
        print(f"The difference in frames between ttl file and recieved RWD system after drift correction is {ttl_check_pi_rwd}")
    # ttl_pi_drift.append(ttl_check[-1]-ttl_check[0]) # no drift between sync pi and camera :)
    
    # this may be why later trials are removed, since anchoring is the earliest case of matches
    anchor_point = np.where(abs(ttl_check_pis) == np.min(abs(ttl_check_pis)))[0][0]
    
    if np.min(abs(ttl_check_pis)) > 2:
        print(f'total session offset of {np.min(abs(ttl_check_pis))} frames detected') 
    
    ttl = ttl.reset_index(drop = True)*sr
    ttl_anchored = ttl - ttl[anchor_point] # now clipped
    
    prey_trial_starts_anchored = prey_trial_starts -prey_trial_starts[anchor_point] 
    unsynced_trial_mask = np.where(prey_trial_starts_anchored - ttl_anchored > mismatch_threshold)[0] 
    if unsynced_trial_mask.size > 0:
        print(f"DLC and RWD trials unsynced [previously removed from analysis]: trials {unsynced_trial_mask}") # removed these but didnt notice effect on my data, so im keeping them in
        
    synced_trials = np.where(prey_trial_starts_anchored - ttl_anchored < mismatch_threshold)[0]     

    # check that indexing is ok, probably will get cleaner data this way
    offset_cam_pi = prey_trial_starts_anchored - ttl_anchored
# %%
    
    
    # knee finder...
    trim = 10 # ignore the last part of the trial (frames)
    prey_trim_array = np.full(len(ttl), trim)
    IR_trim_array = np.full(len(IR_trial_starts), trim)
    
    if isinstance(r_log['IR Trimmer'][l], str):
        trim_temp = r_log['IR Trimmer'][l].replace(" ", "")
        trim_temp = np.array(list(map(str, trim_temp.split(','))))
        for trim_value in trim_temp:
            IR_trim_trial_idx = int(trim_value[0])
            IR_trim_array[IR_trim_trial_idx] = (IR_trial_ends[IR_trim_trial_idx] - 
                                             IR_trial_starts[IR_trim_trial_idx] - 
                                             int(trim_value[2:]))
            
    if isinstance(r_log['Prey Trimmer'][l], str):
        trim_temp = r_log['Prey Trimmer'][l].replace(" ", "")
        trim_temp = np.array(list(map(str, trim_temp.split(','))))
        for trim_value in trim_temp:
            prey_trim_trial_idx = int(trim_value[0])
            prey_trim_array[prey_trim_trial_idx] = (prey_trial_ends[prey_trim_trial_idx] - 
                                             prey_trial_starts[prey_trim_trial_idx] - 
                                             int(trim_value[2:]))

        
    approach_initiation = np.zeros(len(prey_trial_starts), dtype = int)
    IR_initiation = np.zeros(len(IR_trial_starts), dtype = int)
    
    approach_initiation = initiate_finder(prey_trial_starts, prey_trial_ends, prey_snout_distance, approach_initiation, plotknee, prey_trim_array, 'prey', synced_trials)
    IR_initiation = initiate_finder(IR_trial_starts, IR_trial_ends, IR_snout_distance, IR_initiation, plotknee, IR_trim_array, 'IR')
    
    # visually check for drift, useful dont delete
    # plt.figure()
    # plt.plot(offset_cam_pi[synced_trials])
    # boxoff()
    # plt.title("off set between ttl pulse pi and obesrved laser from camera, per trial")
    
    short_prey_trials_idx = np.where(approach_initiation[:4] < (stim_dur-1.5*sr))[0] # select trials < 18 seconds to filter out chase trials and nonapproaches / pointless because i filter later in main script?
    IR_app_idx = np.where(IR_initiation[:4] < (stim_dur-1.5*sr))[0]
    
    # approach_initiation = approach_initiation[synced_trials]
    
    # Compute the union of the two arrays
    app_idx = np.intersect1d(short_prey_trials_idx, synced_trials)        
    
    # Wrap back into a tuple to keep same structure as np.where returns
    initTrace = np.full(len(ttl), np.nan)
    initTrace[app_idx] = round(offset_cam_pi[app_idx]) + approach_initiation[app_idx] # init start (frames) relative to trial start
    approachTrials = eventTS_RWD[app_idx] + initTrace[app_idx] # init start in context of entire 
    approachTrials = [int(x) for x in approachTrials]
        
    
    IR_idx = IR_trial_starts 
    IR_initTrace = np.full(len(IR_idx), np.nan)

    IR_onset_idx = []
    IR_initTrace[IR_app_idx] = IR_initiation[IR_app_idx] # moving IR trials initiation times, relative to IR trial start
    IR_onset_idx = IR_trial_starts[IR_app_idx] # IR laser onset, for approach trials
    IR_approachTrials = IR_onset_idx + IR_initTrace[IR_app_idx] # moving IR trials initiation times, relative to recording start
    IR_approachTrials = [int(x) for x in IR_approachTrials]
    
    # ITI 
    speedITI, ITIidx = ITI_approach_finder(snout_speed, eventTS_RWD, ttl, cleaner_baseline,
                            sr, switch_case, exp_type = 'fm')
    
    
    # Filter out ITIs where animal is out of frame, or ITI near an out of frame period
    if len(out_of_frame) > 0:
        iti_out_of_view = np.array([
            any(start - trial_length <= idx < stop + trial_length for start, stop in out_of_frame)
            for idx in ITIidx
        ])
    
        
        if len(iti_out_of_view) > 0:
            # Keep only ITIs that are NOT during an out-of-view period
            speedITI = speedITI[~iti_out_of_view]
            ITIidx = ITIidx[~iti_out_of_view]



    approach_snout_speed = np.ones((len(approachTrials), trial_length))
    approach_hrC_speed = np.ones((len(approachTrials), trial_length))
    approach_tail_speed = np.ones((len(approachTrials), trial_length))
    approach_movement_before = []
    for index, frame in enumerate(approachTrials):
        approach_snout_speed[index] = snout_speed[frame-before:frame+after]
        approach_hrC_speed[index] = hrC_speed[frame-before:frame+after]
        approach_tail_speed[index] = tail_speed[frame-before:frame+after]
    
    if cleaner_baseline:
       move_mask = np.mean(approach_snout_speed[:,:before],axis = 1) > 25/1000
       if len(move_mask) > 0:
           approach_snout_speed = approach_snout_speed[~move_mask]
           approach_hrC_speed = approach_hrC_speed[~move_mask]
           approach_tail_speed = approach_tail_speed[~move_mask]
       approach_movement_before = move_mask               
                
    approach_speeds = {
        'snout': approach_snout_speed,
        'hrC': approach_hrC_speed,
        'tail': approach_tail_speed
        }
    
    
    
    IR_snout_speed = np.ones((len(IR_approachTrials), trial_length))
    IR_hrC_speed =np.ones((len(IR_approachTrials), trial_length))
    IR_tail_speed = np.ones((len(IR_approachTrials), trial_length))
    IR_movement_before = []
    for index, frame in enumerate(IR_approachTrials):
        IR_snout_speed[index] = snout_speed[frame-before:frame+after]
        IR_hrC_speed[index] = hrC_speed[frame-before:frame+after]
        IR_tail_speed[index] = tail_speed[frame-before:frame+after]
    
    if cleaner_baseline:
       move_mask = np.mean(IR_snout_speed[:,:before],axis = 1) > 25/1000
       if len(move_mask) > 0:
           IR_snout_speed = IR_snout_speed[~move_mask]
           IR_hrC_speed = IR_hrC_speed[~move_mask]
           IR_tail_speed = IR_tail_speed[~move_mask]
       IR_movement_before = move_mask   
       
    IR_speeds = {
        'snout': IR_snout_speed,
        'hrC': IR_hrC_speed,
        'tail': IR_tail_speed}
    
    
    
    
    return (eventTS_RWD, approachTrials, IR_approachTrials, app_idx, initTrace, IR_initTrace,
            approach_speeds, IR_speeds, IR_onset_idx, IR_app_idx, IR_idx,
            drift, drift_ratio, speedITI, ITIidx, approach_movement_before, 
            IR_movement_before)
