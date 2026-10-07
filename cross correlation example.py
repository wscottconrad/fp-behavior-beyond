# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:55:35 2026

@author: conrad
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
rng = np.random.default_rng()

x = np.arange(200) / 200
sig = np.sin(2 * np.pi * x)

shifted_signal = np.sin(2 * np.pi * (x-.10))

corr = signal.correlate(shifted_signal, sig)
lags = signal.correlation_lags(len(sig), len(shifted_signal))
corr /= np.max(corr)
fig, (ax_orig, ax_noise, ax_corr) = plt.subplots(3, 1, figsize=(4.8, 4.8))
ax_orig.plot(sig)
ax_orig.set_title('Original signal')
ax_orig.set_xlabel('Sample Number')
ax_noise.plot(shifted_signal)
ax_noise.set_title('Shifted signal')
ax_noise.set_xlabel('Sample Number')
ax_corr.plot(lags, corr)
ax_corr.set_title('Cross-correlated signal')
ax_corr.set_xlabel('Lag')
ax_orig.margins(0, 0.1)
ax_noise.margins(0, 0.1)
ax_corr.margins(0, 0.1)
fig.tight_layout()
plt.show()