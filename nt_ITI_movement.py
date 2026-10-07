# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 13:44:54 2025

@author: sconrad
"""

def nt_ITI_movement(nt_file, ttl_file, event_timestamps, sr, 
                    eventtimestampsBehind, idn, d, fIdx, behindLaserIndex, 
                    driftTable, l, trialClass, setUp, stim_dur, thresh, fwd_speed,
                    automate_initiate_finder, p_latency_ext, plot_angular_velocity,
                    labelsX, traceTiming, cleaner_baseline):
    
    import numpy as np
    import pandas as pd
    # import matplotlib.pyplot as plt
    # from scipy.signal import resample
    from scipy.signal import resample_poly
    import scipy.io
    import datetime
    from ITI_approach_finder import ITI_approach_finder


    sr_nt = 100  # sample rate of neurotar

    # thresh = 30  # speed threshold for locomotion start [old threshold before 15/10/25]
    start_idx = 5 * sr
    end_idx = 25 * sr


    data = scipy.io.loadmat(nt_file, struct_as_record=False, squeeze_me=True)
    neurotar_data = data['neurotar_data']  # Extract main struct
    timestamps = neurotar_data.SW_timestamp_converted  # Extract SW_timestamp field, 5 is usually when ttl output stops
    timestamps = np.array([datetime.datetime(1, 1, 1) + datetime.timedelta(days=d - 367) for d in timestamps]) #changed from 366

    # neurotar clock is 2 seconds early compared to optottl?
    # print(neurotar_data.TTL_outputs[0:15])
    
    # to get names use: print(my_data.dtype)
    if fwd_speed:
        speed = neurotar_data.Forward_speed
    else:
        speed = neurotar_data.Speed
        
    speed = speed[~np.isnan(speed)]
    angular_velocity = neurotar_data.Angular_velocity
    angular_velocity = angular_velocity[~np.isnan(angular_velocity)]

    ttl = pd.read_csv(ttl_file)
  
    if ttl['Event'][0] == 'Received trigger':
        #subtracts neurotar ttl input from first laser ttl output 
        # code to ignore false recieved triggers from plugging things into same outlet or noise
        if np.any(ttl['Event'] == 'sync'):
            nt_inputs = data['neurotar_data'].TTL_inputs
            
            if np.any(nt_inputs !=0):
                nt_first_sync = np.where(nt_inputs !=0)[0][0]
                nt_outputs = data['neurotar_data'].TTL_outputs
                nt_start = np.where(nt_outputs != 0)[0][0]
                nt_reference = (nt_first_sync - nt_start)/sr_nt
                
                first_sync = ttl.iloc[np.where(ttl['Event'] == 'sync')[0][0]]
                
                all_recieved_triggers = ttl.iloc[np.where(ttl['Event'] == 'Received trigger')[0]]
                art_ms = (pd.to_datetime(first_sync['DateTime']) -
                          pd.to_datetime(all_recieved_triggers['DateTime'])).dt.total_seconds() + (first_sync['Milliseconds'] - 
                                                                                                   all_recieved_triggers['Milliseconds'])/1000
                real_start_trigger = np.where(abs(art_ms-nt_reference) == np.min(abs(art_ms-nt_reference)))[0][0]
                
                ttl = ttl.iloc[real_start_trigger:].reset_index(drop =True)
                
                if abs(int(nt_reference)-int(art_ms[real_start_trigger])) > 3:
                    print(f"Cross talk between control pi and nt is {abs(int(nt_reference)-int(art_ms[real_start_trigger]))}! Need 2 fix")
                # print(f"NT first recieved sync - nt start trigger: {int(nt_reference)} frames")
                # print(f"Sync computer sent - nt start trigger: {int(art_ms[real_start_trigger])} frames (should be similar to above)")

        else:
            number_of_recieved_triggers = len(ttl.loc[ttl['Event'] == 'Received trigger'])
            # print("No sync inputs recieved by nt system. Skipping sync step.") 
            if number_of_recieved_triggers > 1:
                print(f"Noise detected: {number_of_recieved_triggers} recieved triggers (should be 1)")
                
                
            
        
        ttl['Time'] = (pd.to_datetime(ttl['DateTime']) - pd.to_datetime(ttl['DateTime'][0])).dt.total_seconds() + (ttl['Milliseconds'] / 1000) - ttl['Milliseconds'][0] / 1000
               
        ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True) # added because received trigger sometimes added twice
        ttl = ttl.loc[ttl['Event'] !='sync'].reset_index(drop = True)
        ttl['Time'] = ttl['Time'] - stim_dur/sr
        
        if fIdx.size !=0: # removes failed laser trials
            if fIdx.size == 1:
                expanded_fIdx = [fIdx*3 +offset for offset in range(3)]
            else: 
                expanded_fIdx = [i + offset for i in fIdx*3 for offset in range(3)]
            # print(expanded_fIdx)
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
        
        #applies drift correction factor to RWD eventTS to scale with time 
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
    final_length = int(len(speed) * sr / sr_nt) # makes sure speed vector size matches fp signal
    # speed_old = resample(speed, int(len(speed) * sr / sr_nt))

    
    speed = resample_poly(speed, sr, sr_nt)
    speed = speed[0:final_length]


    angular_velocity = angular_velocity[int(clip_start):]
    angular_velocity = resample_poly(angular_velocity, sr, sr_nt)
    angular_velocity = angular_velocity[0:final_length]

    # angular_velocity = resample(angular_velocity, int(len(angular_velocity) * sr / sr_nt))

    fwdSpeed = neurotar_data.Forward_speed
    fwdSpeed = fwdSpeed[~np.isnan(fwdSpeed)]
    fwdSpeed = fwdSpeed[int(clip_start):]
    fwdSpeed = resample_poly(fwdSpeed, sr, sr_nt)
    fwdSpeed = fwdSpeed[0:final_length]
    
    
    speed_c = speed
    
    speed_trials = np.zeros((len(ttl), start_idx + end_idx)) # speed traces, prey laser aligned
    speed_trials_mov = np.full((len(ttl), start_idx + end_idx), np.nan) #  speed traces, initiation aligned
    # speed_approach_trials_init = np.full((len(ttl), start_idx + end_idx), np.nan) #

    approachTrials = np.full(len(ttl), np.nan)
    avoidTrials = np.full(len(ttl), np.nan)
    initTrace = np.full(len(ttl), np.nan)
    # fwd_speed_trials_mov = np.zeros((len(ttl), start_idx + end_idx))
    approach_fwdSpeed = np.full((len(ttl), start_idx + end_idx), np.nan)
    avoid_fwdSpeed = np.full((len(ttl), start_idx + end_idx), np.nan)
    
# speed during trials and finds movement initiation
    
    movement_before = []

    if automate_initiate_finder:
        for t, timestamp in enumerate(ttl):
            speed_trials[t, :] = speed[timestamp - start_idx:timestamp + end_idx] # makes speed trace for trial
            
            if np.mean(speed[timestamp - start_idx:timestamp-(3*sr)]) > 5 and cleaner_baseline: # checks if there is movement before trial start
                movement_before.append(t)
            
            
            if trialClass[t] == 1 or trialClass[t] == 2:
                possible_peak = np.where(speed[timestamp + 15:timestamp + 600] >= thresh * 2 + 10)[0] + 15 # finds all locations where pp could be
        
                for pp in possible_peak:
                    
                    if np.mean(speed[timestamp + pp:timestamp + pp + 150]) >= thresh - 5: # if 150 units after pp is sustained
                        init = np.where(speed[timestamp + pp - 15:timestamp + pp] >= thresh)[0] # is there an initiate before?
                        if init.size > 0:
                            temp_idx = timestamp + pp - 15 + init[0] # tll + peak - 1sec + 
                            speed_trials_mov[t, :] = speed[int(temp_idx) - start_idx:int(temp_idx) + end_idx]
                            initTrace[t] = pp - 15 + init[0] # this is point of movement initiation within a trace                            
                            
                            # filter into approach or avoid
                            if trialClass[t] == 2 and temp_idx - timestamp > 0.5*sr:
                                approach_fwdSpeed[t, :] = fwdSpeed[int(temp_idx) - start_idx:int(temp_idx) + end_idx]
                                approachTrials[t] = temp_idx
                            elif trialClass[t] == 1 and temp_idx - timestamp > 0.5*sr: 
                                avoid_fwdSpeed[t, :] = fwdSpeed[int(temp_idx) - start_idx:int(temp_idx) + end_idx]
                                avoidTrials[t] = temp_idx
         
                            break
         
        
                speed_c[timestamp - start_idx:timestamp + int(end_idx * 1.5)] = 0
     
            
    
    else:
        for t, timestamp in enumerate(ttl):
            speed_trials[t, :] = speed[timestamp - start_idx:timestamp + end_idx] # makes speed trace for trial
    
            if trialClass[t] == 1 or trialClass[t] == 2:
                temp_idx = timestamp + p_latency_ext[t] # tll + manually defined initiate start 
                speed_trials_mov[t, :] = speed[int(temp_idx) - start_idx:int(temp_idx) + end_idx]
                initTrace[t] = p_latency_ext[t] # this is point of movement initiation within a trace             
                
                # # filter into appraoch or avoid
                # if trialClass[t] == 2 and temp_idx - timestamp > 0.5*sr:
                #     approach_fwdSpeed[t, :] = fwdSpeed[int(temp_idx) - start_idx:int(temp_idx) + end_idx]
                #     approachTrials[t] = temp_idx
                # elif trialClass[t] == 1 and temp_idx - timestamp > 0.5*sr: 
                #     avoid_fwdSpeed[t, :] = fwdSpeed[int(temp_idx) - start_idx:int(temp_idx) + end_idx]
                #     avoidTrials[t] = temp_idx
 
                # break
          
         
        
        
                speed_c[timestamp - start_idx:timestamp + int(end_idx * 1.5)] = 0
      
    # print(f"Number of trials with movement before stimulus onset: {len(movement_before)}")    
          
      
    trialClass[np.isnan(approachTrials)] = 0 # reset trial class, deletes trials where approach started too soon (usually <1 sec after prey laser onset)        
    mask = ~np.isnan(approach_fwdSpeed).any(axis=1)
    approach_fwdSpeed = approach_fwdSpeed[mask]     
           
    if len(trialClass[trialClass == 2])>0: 
        angular_velocity_apr = np.full((len(trialClass[trialClass == 2]), start_idx + end_idx), np.nan)
        angular_initiation = np.zeros((len(trialClass[trialClass == 2])), dtype=int)
        # plt.figure()
        for index, approach_index in enumerate(ttl[trialClass == 2]):
            angular_velocity_apr[index] = angular_velocity[approach_index - start_idx:approach_index + end_idx]
            angular_initiation[index] = int(initTrace[~np.isnan(initTrace)][index]+start_idx)
            # line, = plt.plot(angular_velocity_apr[index], alpha = 0.3)
            # color = line.get_color()

            # plt.plot(angular_initiation[index], angular_velocity_apr[index][angular_initiation[index]], 'D', color = color )
            
        
            
    else:
        angular_velocity_apr = np.nan
        angular_initiation = np.nan

        
    app_idx = np.where(~np.isnan(approachTrials))[0]
    avd_idx = np.where(~np.isnan(avoidTrials))[0]

    
    approachTrials = approachTrials[~np.isnan(approachTrials)]
    avoidTrials = avoidTrials[~np.isnan(avoidTrials)]



    speed_iti, iti_idx = ITI_approach_finder(speed_c, event_timestamps, ttl, cleaner_baseline,
                            sr, switch_case = False, exp_type = 'nt')

    # replaced this block with above function
    if False:
        # same but for ITI
        pos_speed_idx = np.where(speed_c[:event_timestamps[-1]] >= thresh * 2 + 10)[0]
        pos_speed_idx = pos_speed_idx[pos_speed_idx > 150] # not during beginning
        pos_speed_idx = np.array([
            idx for idx in pos_speed_idx
            if not np.any((ttl - 150 <= idx) & (idx <= ttl + 750)) # filters out potential ITI starts too close to ttls
        ])
    
        # sus_thresh = 30
        iti_idx = []
        speed_iti = []
        valid_speed_idx = []
    
        for idx in range(len(pos_speed_idx)):
            if ((np.mean(speed_c[pos_speed_idx[idx]:pos_speed_idx[idx] + 90]) >= 60) and #sustained?
    
                (np.mean(speed_c[pos_speed_idx[idx] - 90:pos_speed_idx[idx]]) < 25) and # low activity before?
                            
                ((idx == 0) or ((idx != 0) and not (np.any(abs(valid_speed_idx - pos_speed_idx[idx]) < 30))))): # at least one second after, this might be buggy, changed 19/8/25
            
                ITIinit = np.where(speed_c[pos_speed_idx[idx] - 15:pos_speed_idx[idx]] >= thresh)[0] - 15 + pos_speed_idx[idx] # is there an initiate before?
                
                if ((ITIinit.size > 0) and (ITIinit[0] not in iti_idx) and 
                    (all(ITIinit[0] > item + 300 for item in iti_idx))):
                    valid_speed_idx.append(pos_speed_idx[idx])
                    
                    
                    if (ITIinit[0] > 150): #filters out iti traces that start too soon in the session 
                    
                        if cleaner_baseline: 
                            if np.mean(speed_c[ITIinit[0] - start_idx:ITIinit[0]-(3*sr)]) < 5: # and with movement before 
                            
                                iti_idx.append(ITIinit[0]) # tll + peak - 1sec + 
                                speed_iti.append(speed_c[int(iti_idx[-1]) - start_idx:int(iti_idx[-1]) + end_idx])   
                                
                        else:
                            iti_idx.append(ITIinit[0]) # tll + peak - 1sec + 
                            speed_iti.append(speed_c[int(iti_idx[-1]) - start_idx:int(iti_idx[-1]) + end_idx])   
                        
    
        iti_idx = np.array(iti_idx)
        too_close = np.any(np.diff(iti_idx) < 150)
        if too_close:
            print('Found ITI idx that are <150 frames (5 seconds) apart!')  # True if any pair is < 150 apart
        
        if iti_idx.size != 0:
            # speed_iti = [x for x in speed_iti if len(x) > 0]
            speed_iti = np.array(speed_iti)

 

    return (iti_idx, speed_iti, approachTrials, speed_trials, ttl, 
            app_idx, avd_idx, speed_trials_mov, 
            avoidTrials, initTrace, approach_fwdSpeed, angular_initiation,movement_before)
