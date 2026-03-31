# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:04:35 2026

@author: conrad
"""

import pickle


probability_path = 'E:/Conrad/DLC projects/DLC 06_01_2025/fm_subset_cropped/probabilities.pkl'

with open(probability_path, 'rb') as f:
     probability = pickle.load(f)
     