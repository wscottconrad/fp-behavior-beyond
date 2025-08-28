# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 13:44:54 2025

@author: sconrad
"""

def nt_ITI_movement(nt_file, ttl_file, event_timestamps, sr, 
                    eventtimestampsBehind, idn, d, fIdx, behindLaserIndex, 
                    driftTable, l, trialClass, setUp):
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.signal import resample
    import scipy.io
    import datetime
    
    aprchColors = ["#2E8B57", "#228B22", "#6B8E23", "#8FBC8F", "#20B2AA", "#556B2F", "#70a7ff", "#fca103"]
    avdColors = ["#8B0000", "#B22222", "#DC143C", "#FF4500", "#CD5C5C", "#FA8072", "#fca103"]
    rainbow_colors = ["#FF0000", "#FF7F00", "#FFFF00", "#00FF00", "#0000FF", "#4B0082", "#9400D3"]


    # # to visualize the palette
    # fig, ax = plt.subplots(figsize=(6, 2))

    # # Plot color swatches
    # for i, color in enumerate(aprchColors):
    #     ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))
    
    # # Remove axes and labels
    # ax.set_xlim(0, 6)
    # ax.set_ylim(0, 1)
    # ax.set_xticks([])
    # ax.set_yticks([])
    # ax.set_frame_on(False)
    
    # # Show the palette
    # plt.show()



    sr_nt = 100  # sample rate of neurotar
    stim_dur = 20*sr  # correction factor in seconds, since ttl file records ts of when ttl pulse ends, not starts   

    thresh = 30  # speed threshold for locomotion start
    start_idx = 5 * sr
    end_idx = 25 * sr
    set_up = 120 * sr

    # data = scipy.io.loadmat(nt_file)
    data = scipy.io.loadmat(nt_file, struct_as_record=False, squeeze_me=True)
    neurotar_data = data['neurotar_data']  # Extract main struct
    timestamps = neurotar_data.SW_timestamp_converted  # Extract SW_timestamp field, 5 is usually when ttl output stops
    timestamps = np.array([datetime.datetime(1, 1, 1) + datetime.timedelta(days=d - 367) for d in timestamps]) #changed from 366

    # neurotar clock is 2 seconds early compared to optottl?
    # print(neurotar_data.TTL_outputs[0:15])
    
    # to get names use: print(my_data.dtype)    
    speed = neurotar_data.Speed
    speed = speed[~np.isnan(speed)]

    ttl = pd.read_csv(ttl_file)
    # ttl_file = 'Z:\\Conrad\\Innate_approach\\Data_collection\\Neurotar\\107819_20250312_01_ttl'

  
    if ttl['Event'][0] == 'Received trigger':
        #subtracts neurotar ttl input from first laser ttl output
        ttl['Time'] = (pd.to_datetime(ttl['DateTime']) - pd.to_datetime(ttl['DateTime'][0])).dt.total_seconds() + (ttl['Milliseconds'] / 1000) - ttl['Milliseconds'][0] / 1000
        
        ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True) # added because received trigger sometimes added twice
        ttl = ttl.loc[ttl['Event'] !='sync'].reset_index(drop = True)
        ttl['Time'] = ttl['Time'] - stim_dur/sr
        
        if len(fIdx) > 0: # removes failed laser trials
            expanded_fIdx = [i + offset for i in fIdx*3 for offset in range(3)]
            print(expanded_fIdx)
            ttl = ttl.drop(expanded_fIdx, axis = 0).reset_index(drop = True)
            
        # extracts behind laser index    
        expanded_behindLaserIndex = [i + offset for i in behindLaserIndex*3 for offset in range(3)]   
        ttl = ttl.drop(expanded_behindLaserIndex, axis = 0).reset_index(drop = True)

        clip_start = (ttl['Time'][0] - setUp/sr)*sr_nt # 120 seconds from first ttl, ttl is in time relative to nt start, in nt frames, added 15/8/25
        # ttl1 = int(ttl['Time'][0]*sr_nt)
        
        # plt.figure()
        # plt.plot(speed[ttl1-150:ttl1+600])
        # clip_start = ttl['Time'][0]*sr - set_up - cor_fac # 10616.01, in frames,~20 seconds ahead?   , old way of clipping 


        # for calculating drift between rwd and ttl pi    
        drift = event_timestamps[-1]-event_timestamps[0] - (ttl['Time'].iloc[-1]-ttl['Time'][0])*sr
        drift2 = event_timestamps[3]-event_timestamps[0] - (ttl['Time'].iloc[9]-ttl['Time'][0])*sr
        
        driftTable[l] = [event_timestamps[-1]-event_timestamps[0], drift, event_timestamps[3]-event_timestamps[0], drift2]
        
        #applies drift correction factor to to RWD eventTS scale with time 
        newSpace = np.linspace(event_timestamps[0], event_timestamps[-1], event_timestamps[-1]-event_timestamps[0]-int(drift))
        adjustedTTL = np.array([np.abs(newSpace - t).argmin() for t in event_timestamps])
        ttl = adjustedTTL + event_timestamps[0] # ttl is now aligned to speedclipped data
        # ttl = adjustedTTL  # 15/8/25


    # old way of clipping, 
    # speed = resample(speed, int(len(speed) * sr / sr_nt)) #same frames as rwd
    # speed = speed[int(clip_start):]
    # fwdSpeed = resample(neurotar_data.Forward_speed, int(len(neurotar_data.Forward_speed) * sr / sr_nt))
    # fwdSpeed = fwdSpeed[int(clip_start):]
    
    # new way of clipping 15/8/25 not sure if better or worse
    # clip_end = (max(ttl) + 30 * sr)/sr*sr_nt

    speed = speed[int(clip_start):]
    speed = resample(speed, int(len(speed) * sr / sr_nt))
    fwdSpeed = neurotar_data.Forward_speed
    fwdSpeed = fwdSpeed[~np.isnan(fwdSpeed)]
    fwdSpeed = fwdSpeed[int(clip_start):]
    fwdSpeed = resample(fwdSpeed, int(len(fwdSpeed) * sr / sr_nt))
    
    
    speed_c = speed

    plt.figure()
    plt.plot(speed_c)
    plt.vlines(ttl, ymin=min(speed_c), ymax=max(speed_c), color='k')
    plt.title('Speed with TTL')
    
    # # make trial loop here?
    

    speed_trials = np.zeros((len(ttl), start_idx + end_idx))
    speed_trials_mov = np.zeros((len(ttl), start_idx + end_idx))
    


    trial_idx_init = np.full(len(ttl), np.nan)
    approachTrials = np.full(len(ttl), np.nan)
    avoidTrials = np.full(len(ttl), np.nan)
    initTrace = np.full(len(ttl), np.nan)
    # fwd_speed_trials_mov = np.zeros((len(ttl), start_idx + end_idx))
    approach_fwdSpeed = np.zeros((len(ttl), start_idx + end_idx))
    avoid_fwdSpeed = np.zeros((len(ttl), start_idx + end_idx))
    
    
    # speed during trials and finds movement initiation
    for t, timestamp in enumerate(ttl):
        speed_trials[t, :] = speed[timestamp - start_idx:timestamp + end_idx] # makes speed trace for trial

        if trialClass[t] == 1 or trialClass[t] == 2:
            possible_peak = np.where(speed[timestamp + 15:timestamp + 600] >= thresh * 2 + 10)[0] + 15 # finds all locations where pp could be
    
            for pp in possible_peak:
                
                if np.mean(speed[timestamp + pp:timestamp + pp + 150]) >= thresh - 5: # if 150 units after pp is sustained
                    init = np.where(speed[timestamp + pp - 15:timestamp + pp] >= thresh)[0] # is there an initiate before?
                    if init.size > 0:
                        trial_idx_init[t] = timestamp + pp - 15 + init[0] # tll + peak - 1sec + 
                        speed_trials_mov[t, :] = speed[int(trial_idx_init[t]) - start_idx:int(trial_idx_init[t]) + end_idx]
                        initTrace[t] = pp - 15 + init[0] # this is point of movement initiation within a trace
                        
                        
                        # filter into appraoch or avoid
                        if trialClass[t] == 2:
                            approach_fwdSpeed[t, :] = fwdSpeed[int(trial_idx_init[t]) - start_idx:int(trial_idx_init[t]) + end_idx]
                            approachTrials[t] = trial_idx_init[t]
                        elif trialClass[t] == 1: 
                            avoid_fwdSpeed[t, :] = fwdSpeed[int(trial_idx_init[t]) - start_idx:int(trial_idx_init[t]) + end_idx]
                            avoidTrials[t] = trial_idx_init[t]
     
                        break
 


        speed_c[timestamp - start_idx:timestamp + int(end_idx * 1.5)] = 0
        
    # if idn == '107818' and d == '2025_03_13_':
        
    #     for x in range(len(speed_trials)):
    #         plt.figure()
    #         plt.plot(speed_trials[x], color = rainbow_colors[x], alpha = 0.7)
        
        
    mov_trial_idx = np.where(~np.isnan(trial_idx_init))[0]
    app_idx = np.where(~np.isnan(approachTrials))[0]
    avd_idx = np.where(~np.isnan(avoidTrials))[0]

    
    
    trial_idx_init = trial_idx_init[~np.isnan(trial_idx_init)] #i'm not sure i need trial_idx_init
    approachTrials = approachTrials[~np.isnan(approachTrials)]
    avoidTrials = avoidTrials[~np.isnan(avoidTrials)]


    # same but for ITI
    pos_speed_idx = np.where(speed_c[:event_timestamps[-1]] >= thresh * 2 + 10)[0]
    pos_speed_idx = pos_speed_idx[pos_speed_idx > 150] # not during beginning
    pos_speed_idx = np.array([
        idx for idx in pos_speed_idx
        if not np.any((ttl - 150 <= idx) & (idx <= ttl + 750))
    ])

    # sus_thresh = 30
    iti_idx = []
    speed_iti = []
    iti_idx = []
    valid_speed_idx = []

    for idx in range(len(pos_speed_idx)):
        if ((np.mean(speed_c[pos_speed_idx[idx]:pos_speed_idx[idx] + 90]) >= thresh - 5) and #sustained?
            (np.mean(speed_c[pos_speed_idx[idx] - 90:pos_speed_idx[idx]]) < 10) and # low activity before?
            
        ((idx == 0) or ((idx != 0) and not (np.any(abs(valid_speed_idx - pos_speed_idx[idx]) < 30))))): # at least one second after, this might be buggy, changed 19/8/25
        
            ITIinit = np.where(speed_c[pos_speed_idx[idx] - 15:pos_speed_idx[idx]] >= thresh)[0] - 15 + pos_speed_idx[idx] # is there an initiate before?
            
            if ((ITIinit.size > 0) and (ITIinit[0] not in iti_idx) and 
                (all(ITIinit[0] > item + 300 for item in iti_idx))):
                valid_speed_idx.append(pos_speed_idx[idx])
                iti_idx.append(ITIinit[0]) # tll + peak - 1sec + 
                speed_iti.append(speed_c[int(iti_idx[-1]) - start_idx:int(iti_idx[-1]) + end_idx])   
                

    iti_idx = np.array(iti_idx)
    too_close = np.any(np.diff(iti_idx) < 150)
    print(too_close)  # True if any pair is < 150 apart
    speed_iti = np.array(speed_iti)

    if iti_idx.size != 0:
        plt.figure()
        plt.plot(np.mean(speed_iti, axis=0), color=[0.6, 0.6, 0.6])
        plt.plot(np.mean(speed_trials_mov, axis=0), color=[0.78, 0, 0])
        plt.vlines(150, ymin= 0, ymax=speed_iti.max(), linestyle='--' ,color = 'k')
        ax = plt.gca()
        ax.set_xlim([50, 200])
    
    # for x in app_idx:
    #     # if idn == '107819' and d == '2025_03_12_':
    #     #     print(approachTrials)
    #     #     print(approach_fwdSpeed[x])
    #     plt.plot(approach_fwdSpeed[x], color = aprchColors[x], alpha = 0.7)
    # for x in avd_idx:
    #     plt.plot(avoid_fwdSpeed[x], color = avdColors[x], alpha = 0.7) 
        
    plt.title('Mean ITI Movement Above Threshold')

    return (iti_idx, speed_iti, trial_idx_init, speed_trials, ttl, 
            mov_trial_idx, app_idx, avd_idx, speed_trials_mov, approachTrials, 
            avoidTrials, initTrace)
