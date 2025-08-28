# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 13:11:28 2025
# extracts data from RWD system, processes it, and collates traces

# default traces made for TTL inputs (trials), but additional scripts can be
# added here to make traces around other events of interest (e.g. movement
# initiation, movement during iter trial intervals, etc)

#scott conrad 02/12/2024, adapted from Isis Alonso-Lozares

"""


import sys
sys.path.append('/Users/sconrad/Documents/GitHub/fp-behavior-beyond')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from lpFilter import lpFilter
# from IRLS_dFF import IRLS_dFF
from nt_ITI_movement import nt_ITI_movement
from scipy import signal
import glob

#from scipy.signal import butter, filtfilt

# Define helper functions

# Low-pass filter function
def lpFilter(data, sr, lowpass_cutoff, filt_order, db_atten):
    
    # Design a low-pass filter using Butterworth filter design
    nyquist = 0.5 * sr
    normal_cutoff = lowpass_cutoff / nyquist
    
    # Design a Butterworth filter
    b, a = signal.butter(filt_order, normal_cutoff, btype='low', analog=False)
    
    # Apply the filter to the data
    lp_data = signal.filtfilt(b, a, data)
    return lp_data

# IRLS dF/F function
def IRLS_dFF(exp_signal, iso_signal, IRLS_constant):
    
    import statsmodels.api as sm
    from statsmodels.robust.robust_linear_model import RLM
    from statsmodels.robust.norms import HuberT
    
    # Reshape and apply IRLS regression
    # Perform robust regression with bisquare weight function (equivalent to 'bisquare' in MATLAB)
    iso_signal_with_intercept = sm.add_constant(iso_signal)  # Add intercept term to iso_signal
    model = RLM(exp_signal, iso_signal_with_intercept, M=HuberT(IRLS_constant))  # Robust regression model
    results = model.fit()

    # Extract the coefficients (intercept, slope)
    IRLS_coeffs = results.params  # This will give the intercept and the slope
    
    # Compute the fitted isosbestic signal
    ft_iso_signal = np.polyval(IRLS_coeffs[::-1], iso_signal)  # Coefficients are reversed for numpy polyval
    
    # Calculate dFF
    dFF = (exp_signal - ft_iso_signal) / ft_iso_signal
    
    return dFF, ft_iso_signal

filePath = 'W:\\Conrad\\Innate_approach\\Data_collection\\24.35.01\\'
savePath = 'W:\\Conrad\\Innate_approach\\Data_analysis\\24.35.01\\'
ntFilePth = 'W:\\Conrad\\Innate_approach\\Data_collection\\Neurotar\\'

r_log = pd.read_csv(f"{filePath}\\recordinglog.csv", sep=";", encoding='cp1252')

r_log = r_log[r_log['Exp'] == 'nt']
r_log = r_log[r_log['notes'] != 'no ttl alignment']
r_log = r_log.reset_index()

animalIDs = r_log['ID'].unique()
dates = r_log['Date'].unique()
numChannels = 2


# Preprocessing initialization
lowpass_cutoff = 3  # Low-pass cut-off in Hz
filt_order = 4 # was 0.95 as filter_steepness
db_atten = 90  # dB attenuation
sr = 30  # Sampling rate (Hz or FPS)
ma = 90 * sr  # Moving average window, default = 90 seconds
setUp = 120 * sr  # start data trim from 2 minutes before first trial

# Trace window
pre = 5  # 5 seconds before ttl
post = 25  # 25 seconds after
before = pre * sr
after = post * sr
traceTiming = np.arange(-pre, post + 1, sr / (before + after))
labelsX = np.arange(-pre, post + 1, 5)

driftTable = np.zeros((len(r_log), 4))


for l in range(len(r_log)):
      
    idn = str(r_log['ID'][l])
    d = str(r_log['Date'][l]) 
    
    # # for debugging
    # if ((idn == '112742') and (d == '2025_06_05_')) == False:
    #     continue
    
    # Load data
    ntFile = f"{ntFilePth}Track_[{d.replace('_', '-', 2)}*{idn}_session*\\*.mat"
    ntFile = next(iter(glob.glob(ntFile, recursive = True)), None)
        
    ttlFile = f"{ntFilePth}{idn}_{d.replace('_', '')}_01_ttl"
    rawData = pd.read_csv(f"{filePath}{d}{idn}\\Fluorescence.csv", skiprows=1)
    # eventTS = signal.iloc[:,1]
      
    # rawData = pd.read_csv('Z:\\Conrad\\Innate_approach\\Data_collection\\24.35.01\\2025_03_12_107819\\Fluorescence.csv', skiprows=1) 
    timestamps = rawData.iloc[:, 0] / 1000
        
    eventTS = rawData.iloc[:, 1].str.contains('Input1*2*0', regex = False, na = False)
    eventTS = eventTS[eventTS].index
    
    
    fIdx = []
    # code for removing failed laser presentations
    if isinstance(r_log['Failed'][l], int):
        fIdx = np.array(int(r_log['Failed'][l]))
        eventTS = eventTS.delete(np.array(int(r_log['Failed'][l])))
    elif isinstance(r_log['Failed'][l], float) and  np.isnan(r_log['Failed'][l]) == False:
        fIdx = np.array(list(map(int, str(r_log['Failed'][l]).split('.'))))
        if fIdx[1] == 0:
            fIdx = np.delete(fIdx, 1)
        eventTS = eventTS.delete(fIdx)
        
        
            
    trialClass = np.array([int(x) for x in str(r_log['trial class'][l]).split('.')]) # 0 is NR, 1 is avoid, 2 is approach

        
        
    # code for selecting normal prey laser trials vs behind
    behindLaserIndex = []
    if pd.isnull(r_log['Behind trials'][l]) == False:
        behindLaserIndex =  np.array([int(x) for x in str(r_log['Behind trials'][l]).split('.')])
        eventTSBehind = eventTS[behindLaserIndex]
        eventTS = eventTS.delete(behindLaserIndex)
        
    else:
        eventTSBehind = np.nan
               
        
    clipStart = min(eventTS - setUp)
    clipEnd = max(eventTS) + 30 * sr
    eventTS = eventTS - clipStart
        
    for ch in range(1, numChannels + 1):
        site = r_log[str(ch)][l]
        
        print(f"________________________________________________\nCurrent run: {d}{idn} {site}\n------------------------------------------------")
        
        chIsos = rawData.iloc[int(clipStart):int(clipEnd), 2 * ch]
        chGreen = rawData.iloc[int(clipStart):int(clipEnd), 2 * ch + 1]
        
        # Apply low-pass filter
        lp_normDatG = lpFilter(chGreen, sr, lowpass_cutoff, filt_order, db_atten)
        lp_normDatI = lpFilter(chIsos, sr, lowpass_cutoff, filt_order, db_atten)
        
        # fig = plt.figure()
        # plt.plot(lp_normDatG)
        # plt.plot(lp_normDatI)
        
        # Fit isosbestic signal to green channel and compute dF/F
        dFoF, ft_iso_signal = IRLS_dFF(lp_normDatG, lp_normDatI, 3)
        plt.title(f'{idn} {d} {site} Channel {ch} dF/F')

        plt.plot(dFoF)
        
        # Prepare to store traces
        traces = np.full((len(eventTS), before + after), np.nan)
        tracesGraw = np.full((len(eventTS), before + after), np.nan)
        tracesIraw = np.full((len(eventTS), before + after), np.nan)
        
        # Compute the indices for the window
        startIdx = eventTS - before
        endIdx = eventTS + after 
        
        # get nt data
        if ch == 1:
            (ITIidx, speedITI, trialIdxInit, speedTrials, ttl, 
             movTrialIdx, app_idx, avd_idx,speedTrialsMov, 
             approachTrials, avoidTrials, initTrace) = nt_ITI_movement(ntFile, ttlFile, 
                                                                       eventTS, sr,
                                                          eventTSBehind, 
                                                          idn, d, fIdx, behindLaserIndex,
                                                          driftTable, l, trialClass, setUp);
        
        # Extract trial traces
        for m in range(len(eventTS)):
            traces[m] = dFoF[startIdx[m]:endIdx[m]]
            tracesGraw[m] = lp_normDatG[startIdx[m]:endIdx[m]]
            tracesIraw[m] = lp_normDatI[startIdx[m]:endIdx[m]]
        
        # # Align traces to movement if necessary. i do this to show isosbestic doesnt change with movement 
        # tracesInit = np.full((len(trialIdxInit), before + after), np.nan)
        tracesInitGraw = np.full((len(trialIdxInit), before + after), np.nan)
        tracesInitIraw = np.full((len(trialIdxInit), before + after), np.nan)
        
        trialIdxInit = trialIdxInit.astype(np.int64)
        
        if len(trialIdxInit) > 0:
            for m in range(len(trialIdxInit)):
                if np.size(dFoF[trialIdxInit[m] - before:trialIdxInit[m] + after]) == 900:
                    # tracesInit[m] = dFoF[trialIdxInit[m] - before:trialIdxInit[m] + after]
                    tracesInitGraw[m] = lp_normDatG[trialIdxInit[m] - before:trialIdxInit[m] + after]
                    tracesInitIraw[m] = lp_normDatI[trialIdxInit[m] - before:trialIdxInit[m] + after]
                
        else:
            # tracesInit = np.nan
            tracesInitGraw = np.nan
            tracesInitIraw = np.nan
        
        
        # traces for approach
        appTracesInit = np.full((len(approachTrials), before + after), np.nan)
        
        approachTrials = approachTrials.astype(np.int64)
        
        if len(approachTrials) > 0:
            for m in range(len(approachTrials)):
                if np.size(dFoF[approachTrials[m] - before:approachTrials[m] + after]) == 900:
                    appTracesInit[m] = dFoF[approachTrials[m] - before:approachTrials[m] + after]
                
        else:
            appTracesInit = np.nan
            
            
            
        # traces for avoid
        avdTracesInit = np.full((len(avoidTrials), before + after), np.nan)
        
        avoidTrials = avoidTrials.astype(np.int64)
        
        if len(avoidTrials) > 0:
            for m in range(len(avoidTrials)):
                if np.size(dFoF[avoidTrials[m] - before:avoidTrials[m] + after]) == 900:
                    avdTracesInit[m] = dFoF[avoidTrials[m] - before:avoidTrials[m] + after]
                
        else:
            avdTracesInit = np.nan
            
        # traces for No Response (NR)
        NR_idx = np.where(trialClass == 0)[0] # index of trials
        NRtrials = ttl[NR_idx] #affected by change in ttl 15/8/25
        NRtraces = np.full((len(NRtrials), before + after), np.nan)
        
        NRtrials = NRtrials.astype(np.int64)
        
        if len(NRtrials) > 0:
            for m in range(len(NRtrials)):
                if np.size(dFoF[NRtrials[m] - before:NRtrials[m] + after]) == 900:
                    NRtraces[m] = dFoF[NRtrials[m] - before:NRtrials[m] + after]
                
        else:
            NRtraces = np.nan

        
        # for ITIs
        tracesITI = np.full((len(ITIidx), before + after), np.nan)
        tracesITIGraw = np.full((len(ITIidx), before + after), np.nan)
        tracesITIIraw = np.full((len(ITIidx), before + after), np.nan)

        for m, idx in enumerate(ITIidx):
            if idx is not None:
                tracesITI[m] = dFoF[idx - before:idx + after]
                tracesITIGraw[m] = lp_normDatG[idx - before:idx + after]
                tracesITIIraw[m] = lp_normDatI[idx - before:idx + after]
            else:
                tracesITI = np.nan
                tracesITIGraw = np.nan
                tracesITIGraw = np.nan
        
        # Collate the data # i dont see a differnce in shape
        traceData = np.vstack(traces)
        traceDataG = np.vstack(tracesGraw)
        traceDataI = np.vstack(tracesIraw)
        
        
        #################################
        # Baseline correction (z-scoring)
        #################################
        
        # not taking trial type into consideration:
        traceDataSD = np.std(traceData[:, :pre * sr], axis=1) # axis 1 is along row
        ZdFoF = (traceData - np.mean(traceData[:, :pre * sr], axis=1).reshape(-1, 1)) / traceDataSD.reshape(-1, 1)
        
        traceDataSDG = np.std(traceDataG[:, :pre * sr], axis=1)
        Gdata = (traceDataG - np.mean(traceDataG[:, :pre * sr], axis=1).reshape(-1, 1)) / traceDataSDG.reshape(-1, 1)
        
        traceDataSDI = np.std(traceDataI[:, :pre * sr], axis=1)
        Idata = (traceDataI - np.mean(traceDataI[:, :pre * sr], axis=1).reshape(-1, 1)) / traceDataSDI.reshape(-1, 1)
        
        # APPROACH INITIATION and PREYLASER ALIGNED
        if approachTrials is not np.nan and len(approachTrials) > 0:
            ZdFoFApproach = (appTracesInit - np.mean(traceData[app_idx,:pre*sr],axis=1).reshape(-1, 1)) / traceDataSD[app_idx].reshape(-1, 1)
            ZdFoFApproach_trialOnset = (traceData[app_idx] - np.mean(traceData[app_idx,:pre*sr],axis=1).reshape(-1, 1)) / traceDataSD[app_idx].reshape(-1, 1)

        else:
            ZdFoFApproach = np.nan
            ZdFoFApproach_trialOnset = np.nan
           
            
        # AVOID INITIATION and PREYLASER ALIGNED
        if avoidTrials is not np.nan and len(avoidTrials) > 0:
            ZdFoFAvoid = (avdTracesInit - np.mean(traceData[avd_idx,:pre*sr],axis=1).reshape(-1, 1)) / traceDataSD[avd_idx].reshape(-1, 1)
            ZdFoFAvoid_trialOnset = (traceData[avd_idx] - np.mean(traceData[avd_idx,:pre*sr],axis=1).reshape(-1, 1)) / traceDataSD[avd_idx].reshape(-1, 1)

        else:
            ZdFoFAvoid = np.nan
            ZdFoFAvoid_trialOnset = np.nan
            
        # and for NR (only PREYLASER ALIGNED)
        if NRtrials is not np.nan and len(NRtrials) > 0:
            ZdFoFNR = (NRtraces - np.mean(traceData[NR_idx,:pre*sr],axis=1).reshape(-1, 1)) / traceDataSD[NR_idx].reshape(-1, 1)

        else:
            ZdFoFNR = np.nan

            
            
        # and for ITI
        if tracesITI is not np.nan and len(tracesITI) > 0:
            tracesITIV = np.vstack(tracesITI);

            traceDataSDITI = np.std(tracesITIV[:,:pre*sr],axis = 1)
            ZdFoFITI = (tracesITIV - np.mean(tracesITIV[:,:pre*sr], axis=1).reshape(-1,1))/ traceDataSDITI.reshape(-1, 1)

            tracesITIGraw = np.vstack(tracesITIGraw)
            traceDataSDG = np.std(tracesITIGraw[:,:pre*sr],axis=1)
            ITIGdata = (tracesITIGraw - np.mean(tracesITIGraw[:,:pre*sr],axis=1).reshape(-1,1)) /traceDataSDG.reshape(-1, 1)

            tracesITIIraw = np.vstack(tracesITIIraw)
            traceDataSDI = np.std(tracesITIIraw[:,:pre*sr], axis=1)
            ITIIdata = (tracesITIIraw - np.mean(tracesITIIraw[:,:pre*sr],axis=1).reshape(-1,1)) /traceDataSDI.reshape(-1, 1)
            
            plt.figure()
            plt.plot(np.mean(tracesITIV, axis = 0))
            plt.fill_between(range(0,tracesITIV.shape[1]),
                             np.mean(tracesITIV, axis = 0) + np.std(tracesITIV, axis=0) / np.sqrt(len(tracesITIV)),
                             np.mean(tracesITIV, axis = 0) - np.std(tracesITIV, axis=0) / np.sqrt(len(tracesITIV)),
                             alpha=0.3)
            plt.vlines(150, ymin = -0.01, ymax = 0.05, linestyle = '--')
            plt.title(f'{idn} {d} {site} Channel {ch} ITI')

            plt.show()

        else:
            ZdFoFITI = np.nan
            ITIGdata = np.nan
            ITIIdata = np.nan

# % figure;
# % imagesc(ZdFoFITI);  
# % colorbar;  
# % xline(150, '--w', 'Linewidth', 1.5)
# % % yline(length(ZdFoF(:,1))/2+0.5, '-k', 'Linewidth', 1.5)
# % xlabel('Time');
# % ylabel('ITI Trial');
# % title([animalIDs{idn}, ' Channel ' num2str(ch) ' dF/F']);
# % ax = gca;
# % box off
# % xticks(1:150:900);
# % xticklabels(labelsX)

# % collated
# % baselineSD = std(traceData(:,1:pre*sr)); % sd of baseline, for z-scoring later % check this!!!!!!!!!!!!!!!!
# % traceDataCol = traceDataCol - mean(traceDataCol(1:pre*sr)); % baseline subtraction

# % plot first trial of session
# % figure;
# % plot(traceTiming, ZdFoF(1,:), 'Color', '#7b8be8')
# % hold on;
# % xline(0, '--k'); yline(0, '--k');
# % box off
# % xlabel('Time (s)');
# % ylabel('Normalized dF/F');
# % title([animalIDs{idn}, 'Channel ' num2str(ch) ' averaged dF/F over session']);
        
        # Plotting heatmaps
        plt.figure()
        plt.imshow(ZdFoF, aspect='auto', cmap='viridis', interpolation='none')
        plt.colorbar()
        plt.axvline(150, linestyle='--', color='white', linewidth=1.5)
        for wee in range(len(initTrace)):
            if pd.isnull(initTrace[wee]) == False:
                plt.vlines(x = initTrace[wee] + 150, linestyle='--', color='red', linewidth=1.5,
                           ymin = wee - 0.5, ymax = wee + 0.5 )
        plt.xlabel('Time')
        plt.ylabel('Trial')
        plt.title(f'{idn} {d} {site} Channel {ch} Z dF/F')
        
        # Set x-ticks
        plt.xticks(np.arange(0, len(traceTiming), 150), labelsX)
        
        # Get y-tick positions and labels
        yticks = np.arange(ZdFoF.shape[0])  # Tick positions
        ytick_labels = [str(y) for y in yticks]  # Default labels
        
        # Apply y-ticks to the plot
        plt.yticks(yticks, ytick_labels)  
        
        # Get the current axis and apply custom tick colors
        ax = plt.gca()  # Get current axis
        tick_colors = ['green' if y in app_idx else 'red' if y in avd_idx else 'black' for y in yticks]
        for y, color in zip(yticks, tick_colors):
            ax.get_yticklabels()[y].set_color(color)  # Apply color change
        
        # Show plot
        plt.show()

        # Save data
        sesdat = {
            'session': d,
            'mouse': idn,
            'conversion': sr,
            'lp_normDat': lp_normDatG,
            
            'ZdFoF': ZdFoF,
            # 'ZdFoFinit': ZdFoFinit,
            'ZdFoFITI': ZdFoFITI,
            
            # preylater onset locked
            'ZdFoFApproach_trialOnset': ZdFoFApproach_trialOnset if len(app_idx) > 0  else np.nan,
            'ZdFoFAvoid_trialOnset': ZdFoFAvoid_trialOnset if len(avd_idx) > 0 else np.nan,
            'ZdFoFNR': ZdFoFNR if len(NR_idx) > 0 else np.nan,
            
            # movement locked:
            'ZdFoFApproach': ZdFoFApproach if len(app_idx) > 0  else np.nan,
            'ZdFoFAvoid': ZdFoFAvoid if len(avd_idx) > 0 else np.nan,
            
            
            'Gdata': Gdata,
            'Idata': Idata,
            
            # 'InitGdata': InitGdata,
            # 'InitIdata': InitIdata,
            
            'ITIGdata': ITIGdata,
            'ITIIdata': ITIIdata,
            
            'channel': ch,
            
            'site': site,
            
            
            # can add  else np.nan
    
            'speedTrials': speedTrials if ch == 1 else np.nan,
            'speedITI': speedITI if ch == 1 else np.nan,
            'speedTrialsMov': speedTrialsMov if ch == 1 else np.nan,
                
                
                
                }


        save_file = f"{savePath}{d}{idn} Channel {ch}.pkl"
        pd.to_pickle(sesdat, save_file)



col1 = driftTable[:, 0].reshape(-1, 1)  # First column
col3 = driftTable[:, 2].reshape(-1, 1)  # Third column

col2 = driftTable[:, 1].reshape(-1, 1)  # Third column
col4 = driftTable[:, 3].reshape(-1, 1)  # Third column

# Stack them vertically
xx = np.vstack((col1, col3))
yy = np.vstack((col2, col4))


xx = np.delete(xx, np.argmax(xx))
yy = np.delete(yy, np.argmax(yy))


plt.figure()
plt.scatter(xx, yy, color = 'black')
plt.plot(np.unique(xx), np.poly1d(np.polyfit(xx, yy, 1))(np.unique(xx)))





# go to combineData. next


# to do later, stamp within laser trials based on action, direction, other behavior
