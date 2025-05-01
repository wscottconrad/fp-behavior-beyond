# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 13:33:25 2025

@author: sconrad
"""



def IRLS_dFF(exp_signal, iso_signal, IRLS_constant):
    
    import numpy as np
    import statsmodels.api as sm
    from statsmodels.robust import RLM
    
    # Check if the signals have the same length
    if len(exp_signal) != len(iso_signal):
        raise ValueError("exp_signal and iso_signal must have the same length.")
    
    # Reshape and apply IRLS regression
    # Perform robust regression with bisquare weight function (equivalent to 'bisquare' in MATLAB)
    iso_signal_with_intercept = sm.add_constant(iso_signal)  # Add intercept term to iso_signal
    model = RLM(exp_signal, iso_signal_with_intercept, M=sm.robust.norms.Bisquare(IRLS_constant))  # Robust regression model
    results = model.fit()

    # Extract the coefficients (intercept, slope)
    IRLS_coeffs = results.params  # This will give the intercept and the slope
    
    # Compute the fitted isosbestic signal
    ft_iso_signal = np.polyval(IRLS_coeffs[::-1], iso_signal)  # Coefficients are reversed for numpy polyval
    
    # Calculate dFF
    dFF = (exp_signal - ft_iso_signal) / ft_iso_signal
    
    return dFF, ft_iso_signal
