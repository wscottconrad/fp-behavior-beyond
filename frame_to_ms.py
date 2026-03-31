# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 10:46:38 2025

@author: conrad
"""
import math
frame_to_ms = False
fps = 30

if frame_to_ms:
    t1 = '21.04'
    t2 = '28.00'
    
    # split into seconds and frames
    t1s, t1f = map(int, t1.split('.'))
    t2s, t2f = map(int, t2.split('.'))
    
    # convert to seconds
    t1_seconds = t1s + t1f / fps
    t2_seconds = t2s + t2f / fps
    
    latency = t2_seconds - t1_seconds
    print(latency)

else:
    # time frame to total frame
    frame = '06:49:22'
    
    # split into seconds and frames
    t1m, t1s, t1f = map(int, frame.split(':'))
    
    # convert to seconds
    total_frames = (t1s + t1m*60)*fps + t1f
    
    print(total_frames)
    
    
    
frame_count = 3555
conv = frame_count/fps
leftover_frames = str(frame_count%fps)
frames_remain = frame_count-frame_count%fps
seconds = str(int((frames_remain/fps)%fps))
minutes = str(int(math.floor((frames_remain/fps)/60)))
time_stamp = minutes + ':' + seconds + ':' + leftover_frames
print(time_stamp)
