# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 13:27:10 2025

@author: sconrad
"""




def lpFilter(data, samp_rate, lowpass_cutoff, filt_steepness, db_atten):
    from scipy import signal
    
    # Design a low-pass filter using Butterworth filter design
    nyquist = 0.5 * samp_rate
    normal_cutoff = lowpass_cutoff / nyquist
    order = int(filt_steepness)  # The steepness will determine the filter order
    
    # Design a Butterworth filter
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    
    # Apply the filter to the data
    lp_data = signal.filtfilt(b, a, data)
    
    return lp_data
