# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 11:14:01 2026

@author: conrad
"""
import numpy as np


def ITI_approach_finder(speed, event_timestamps, ttl, cleaner_baseline,
                        sr, switch_case, exp_type):
    
    ttl = np.asarray(ttl, dtype=np.int64)
    start_idx = 5 * sr # 5 seconds before
    end_idx = 25 * sr # 25 after
    
    trial_length = start_idx + end_idx
    
    # set parameters 
    if exp_type == 'nt':
        thresh_detect = 30 # possible peak detection, units are mm/s
        thresh_sustain = 60 # is the movement sustained
        thresh_still = 25 # is the animal still before
        strict_thresh_still = 5 # for cleaner baseline
        ceiling = 4000
        
        sustain_period = 3*sr
        
        pre_window = 3

    
    elif exp_type == 'fm':
       
        speed = speed*1000 # convert to m/s
        
        thresh_detect = 80 # possible peak detection
        thresh_sustain = 125 # is the movement sustained
        thresh_still = 75 # is the animal still before
        strict_thresh_still = 5 # for cleaner baseline
        ceiling = 4000
        
        sustain_period = 2*sr 

        cleaner_baseline = False # too stringent for ITI experiments        
        
        pre_window = 8 # seconds


    else:
        raise ValueError(f"Unknown exp_type: {exp_type}")
    
    if switch_case:
        last_ttl = event_timestamps[-3] # case where fibers switched mid session
    else:
        last_ttl = event_timestamps[-1]
        
    pos_speed_idx = np.where(speed[:last_ttl] >= thresh_detect)[0]
    pos_speed_idx = pos_speed_idx[pos_speed_idx > 150] # not during beginning
    # pos_speed_idx = np.array([
    #     idx for idx in pos_speed_idx
    #     if not np.any((ttl - 150 <= idx) & (idx <= ttl + 900)) # filters out potential ITI starts too close to ttls... changed from 750 to 900 16/9/26
    # ])
    
    # --------------------------------------------------
    # Remove candidates close to TTLs
    # --------------------------------------------------
    starts = ttl - pre_window * sr
    ends = ttl + trial_length
    
    j = np.searchsorted(starts, pos_speed_idx, side='right') - 1
    
    keep = (j < 0) | (pos_speed_idx > ends[j])
    
    pos_speed_idx = pos_speed_idx[keep]
    
    if len(pos_speed_idx) == 0:
        return np.array([]), np.array([], dtype=int)
    
    pos_speed_idx = pos_speed_idx[pos_speed_idx>pre_window * sr]
    
    iti_idx = []
    speed_iti = []
    valid_speed_idx = []

    for idx in range(len(pos_speed_idx)):
        if ((np.mean(speed[pos_speed_idx[idx]:pos_speed_idx[idx] + sustain_period]) >= thresh_sustain) and #sustained?

            (np.mean(speed[pos_speed_idx[idx] - (pre_window * sr):pos_speed_idx[idx]]) < thresh_still) and # low activity before?
                        
            ((idx == 0) or ((idx != 0) and not (np.any(abs(valid_speed_idx - pos_speed_idx[idx]) < 30)))) and # at least one second after, this might be buggy, changed 19/8/25
            
            not (np.any(speed[pos_speed_idx[idx]:pos_speed_idx[idx] + 90] > ceiling))): 
        
            ITIinit = np.where(speed[pos_speed_idx[idx] - 15:pos_speed_idx[idx]] >= thresh_detect)[0] - 15 + pos_speed_idx[idx] # is there an initiate before?
            
            if ((ITIinit.size > 0) and (ITIinit[0] not in iti_idx) and 
                (all(ITIinit[0] > item + 300 for item in iti_idx))):
                valid_speed_idx.append(pos_speed_idx[idx])
                
                
                if (ITIinit[0] > 150) and len(speed[int(ITIinit[0]) - start_idx:int(ITIinit[0]) + end_idx]) == trial_length: #filters out iti traces that start too soon in the session 
                
                    if cleaner_baseline: 
                        if np.mean(speed[ITIinit[0] - start_idx:ITIinit[0]-((3 * sr))]) < strict_thresh_still: # and with movement before 
                        
                            iti_idx.append(ITIinit[0]) # tll + peak - 1sec + 
                            speed_iti.append(speed[int(iti_idx[-1]) - start_idx:int(iti_idx[-1]) + end_idx])   
                            
                    else:
                        iti_idx.append(ITIinit[0]) # tll + peak - 1sec + 
                        speed_iti.append(speed[int(iti_idx[-1]) - start_idx:int(iti_idx[-1]) + end_idx])   
                    

    iti_idx = np.array(iti_idx)
    too_close = np.any(np.diff(iti_idx) < 150)
    if too_close:
        print('Found ITI idx that are <150 frames (5 seconds) apart!')  # True if any pair is < 150 apart
    
    if iti_idx.size != 0:
        speed_iti = np.array(speed_iti)
        
    return (speed_iti, iti_idx)