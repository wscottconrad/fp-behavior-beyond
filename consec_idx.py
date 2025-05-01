# -*- coding: utf-8 -*-
"""
Adapted on Thu Mar  6 15:24:06 2025

@author: sconrad
"""

import numpy as np

def consec_idx(indices, threshold):
    """
    Derives logical index of consecutive indices >= threshold.
    e.g. if indices = [1, 3, 5, 6, 7, 8, 10, 11, 12, 15], threshold = 3
           c_idx = [0, 0, 1, 1, 1, 1, 1, 1, 1, 0]
     i.e. indices[c_idx] = [5, 6, 7, 8, 10, 11, 12]
     
     Adapted from matlab script written by Philip Jean-Richard-dit-Bressel, UNSW Sydney (copyright 2019)
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    """
    # indices = tmp
    # threshold = thres
    # Create a mask where consecutive indices are identified
    k = np.concatenate(([True], np.diff(indices) != 1)) #finds start
    

    # Calculate cumulative sum of non-consecutive points
    s = np.cumsum(k)

    # Calculate histogram of the cumulative sum
    x = np.bincount(s)[1:] # added 1: because a 0 was addedin first index for some reason

    # Find indices where the consecutive counts meet the threshold
    idx = np.where(k)[0]
    consecutive = idx[x >= threshold]

    # Initialize result with False values
    c_idx = np.zeros_like(indices, dtype=bool)

    # Iterate over the consecutive indices to mark the valid ones
    for c in range(len(consecutive)):
        x_idx = consecutive[c] == idx
        c_idx[consecutive[c]:((consecutive[c]+x[x_idx]-1)[0])] = True # added [0] 

    return c_idx
