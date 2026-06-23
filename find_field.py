# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:04:40 2026

Looks for fields in matlab structs, return path

@author: conrad
"""
import numpy as np

def find_field(obj, keyword, path='', visited=None):
    if visited is None:
        visited = set()

    results = []

    if id(obj) in visited:
        return results
    visited.add(id(obj))

    # --- Case 1: MATLAB struct ---
    if hasattr(obj, '_fieldnames'):
        for field in obj._fieldnames:
            try:
                value = getattr(obj, field)
            except Exception:
                continue

            current_path = f"{path}.{field}" if path else field

            if keyword.lower() in field.lower():
                results.append(current_path)

            results.extend(find_field(value, keyword, current_path, visited))

    # --- Case 2: NumPy array (MATLAB cell arrays land here) ---
    elif isinstance(obj, np.ndarray):
        for idx, item in np.ndenumerate(obj):
            current_path = f"{path}{list(idx)}"
            results.extend(find_field(item, keyword, current_path, visited))

    return results