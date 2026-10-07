# -*- coding: utf-8 -*-
"""
Created on Wed May 14 12:40:35 2025

@author: conrad

if needed, this code corrects fish-eye distortion from dlc data then converts it to real-
world coordinates (adapted from Alexander Heimel's matlab code)
                   
data is then filtered out if there dlc gives low probability, if a point jumps 
to physiological improbable speeds, or if point jumps out of arena. 

future goal of this code is to get measures from dlc data such as animal location, speed,
angle of head relative to body, tail capture angle?

and possibly automatic behavior detection such as rearing, grooming, jumping,
scratching, gnawing (at object or wall)
"""
import numpy as np
import pandas as pd
import scipy.io
import glob
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import PchipInterpolator as pchip
from scipy.ndimage import gaussian_filter1d
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from boxoff import boxoff
import os


debug = False
plot_need = False
reanalyze = True

plot_angles = False
plot_traj = False

sr = 30
stim_dur = 20
setUp = 120 # seconds

def cart2pol(x, y):
    phi = np.arctan2(y, x)
    rho = np.sqrt(x**2 + y**2)
    return(phi, rho)

def pol2cart(phi, rho):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


def nt_change_overhead_to_camera_coordinates(overhead_x,overhead_y,params):
    distort = params['overhead_camera_distortion']
    
    distort_type = params['overhead_camera_distortion_method']

    overhead_x = overhead_x - params['overhead_camera_width']/2 + params['overhead_camera_image_offset'][0]
    overhead_y = overhead_y - params['overhead_camera_height']/2 + params['overhead_camera_image_offset'][1]
    
    distance_neurotar_center_to_camera_mm = distort[0] 
    focal_distance_pxl = distort[1]

    theta,overhead_r = cart2pol(overhead_x, overhead_y)

    if np.any(overhead_r>focal_distance_pxl):
        print('Point outside camera view')
        overhead_r[overhead_r > focal_distance_pxl] = focal_distance_pxl
        camera_x = np.nan
        camera_y = np.nan
     
    if distort_type == 'fisheye_othographic':
        camera_r = distance_neurotar_center_to_camera_mm * np.tan(np.arcsin(overhead_r / focal_distance_pxl))
        
    else:
        camera_r = distance_neurotar_center_to_camera_mm * overhead_r / focal_distance_pxl
        
        
    camera_x, camera_y = pol2cart(theta, camera_r)
    
    return(camera_x, camera_y)

def nt_change_camera_to_arena_coordinates(camera_x,camera_y,params):
    # invert overhead_center_position
    camera_center_x, camera_center_y = nt_change_overhead_to_camera_coordinates(
        params['overhead_arena_center'][0],
        params['overhead_arena_center'][1],
        params
        )


    # move center of neurotar to center position in camera coordinates
    camera_x = camera_x - camera_center_x
    camera_y = camera_y - camera_center_y
    
    alpha = -params['overhead_camera_angle']
    rotation = np.array([
        [np.cos(alpha),  np.sin(alpha)],
        [-np.sin(alpha), np.cos(alpha)]
    ])
    p = rotation @ np.array([camera_x, camera_y])
    
    arena_x = p[0, :]
    arena_y = p[1, :]
    
    return arena_x, arena_y

def nt_change_overhead_to_arena_coordinates(overhead_x,overhead_y,params, centered_desired):
    camera_x, camera_y = nt_change_overhead_to_camera_coordinates(overhead_x,overhead_y,params)
    arena_x,arena_y = nt_change_camera_to_arena_coordinates(camera_x,camera_y,params)
    
    if not centered_desired:
       # Map back to original video coordinate system
       arena_x = arena_x + params['overhead_camera_width'] / 2 \
                           - params['overhead_camera_image_offset'][0]
       arena_y = arena_y + params['overhead_camera_height'] / 2 \
                           - params['overhead_camera_image_offset'][1]
    return arena_x, arena_y

def timestamp_to_frame(ts):
    minutes, seconds, frames = map(int, ts.split(':'))
    return (minutes * 60 + seconds) * 30 + frames

def find_distance(item1, item2):
    dx = item1.iloc[:,0]-item2.iloc[:,0]
    dy = item1.iloc[:,1]-item2.iloc[:,1]
    distance = np.sqrt(dx**2 + dy**2)
    distance[np.isnan(dx) | np.isnan(dy)] = np.nan
    distance = distance.reset_index(drop = True)
    return distance

def calculate_angle(df1, df2, df3):
    """
    Calculates the signed angle ABC for every timestamp in df3. 
    
    Parameters
    ----------
    df1, df2, df3 : pandas.DataFrame
        Must contain time, x, y columns.
        df3 may have fewer timestamps.
        
    df1 should be snout, df2 should be hrB, and df3 could be back1 or laser

    Returns
    -------
    DataFrame
        df3 with an added 'angle_deg' column.
    """
    merged = (
        df1.join(df2, how="inner", lsuffix="_1", rsuffix="_2")
           .join(df3, how="inner")
    )

    # Extract x,y coordinates from each dataframe
    A = df1.columns[0][0]
    B = df2.columns[0][0]
    C = df3.columns[0][0]

    Ax = merged[(A, "x")]
    Ay = merged[(A, "y")]

    Bx = merged[(B, "x")]
    By = merged[(B, "y")]

    Cx = merged[(C, "x")]
    Cy = merged[(C, "y")]

    BAx = Ax - Bx
    BAy = Ay - By

    BCx = Cx - Bx
    BCy = Cy - By

    cross = BAx * BCy - BAy * BCx
    dot = BAx * BCx + BAy * BCy

    merged["angle_deg"] = np.degrees(np.arctan2(cross, dot))

    return merged[["angle_deg"]]

def plot_angles_trials(angles, relative, plot_angles, bout_id = np.nan):
    if plot_angles:
        
        if not relative: 
            bout_id = (angles.index.to_series().diff() > 1).cumsum()
            
        angles["bout"] = bout_id
        plt.figure(figsize=(8,5))

        for _, bout in angles.groupby("bout"):
            x = np.arange(len(bout))      # 0,1,2,...
            # y = bout["angle_deg"]
            y = gaussian_filter1d(bout["angle_deg"],1) # light smoothing
            plt.plot(x, y, alpha=0.7)
        
        plt.xlabel("Frames since bout start")
        plt.ylabel("Signed angle (°)")
        plt.show()
        
        return bout_id

def speed_calculator(x,y):

    displacement = (np.array([x[:-1], 
                          y[:-1]]) - 
                 np.array([x[1:], 
                          y[1:]]))
    
    distance = np.sqrt(displacement[0]**2 + 
                             displacement[1]**2)
    frame_to_sec = 1/sr
    speed = distance/frame_to_sec/1000 # speed is m/s
    
    return speed

def my_heatplot(data, title, bins=100):
    colors = [(0, 0, 1), (0, 1, 1), (0, 1, 0.75), (0, 1, 0), (0.75, 1, 0),
              (1, 1, 0), (1, 0.8, 0), (1, 0.7, 0), (1, 0, 0)]
    cm = LinearSegmentedColormap.from_list('sample', colors)
    x = data.iloc[:, 0]
    y = data.iloc[:, 1]

    # Bin coordinates into a 2D histogram
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)

    plt.imshow(
        heatmap.T, origin='lower',
        cmap=cm,
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        aspect='auto'
    )
    plt.colorbar(label="Counts")
    plt.title(title)
    plt.xlabel("X position")
    plt.ylabel("Y position")
    boxoff()
    plt.show()
    
def my_trajectory(data, title):
    x = data.iloc[:, 0]
    y = data.iloc[:, 1]
    
    fig, ax = plt.subplots()

    my_circle = Circle(xy=(0, 0), radius=175, facecolor='none', edgecolor='black')
    
    # Add the circle to the axes
    ax.add_patch(my_circle)
    plt.plot(x, y, color="blue", alpha=0.6)
    # plt.xlim((0, 582))
    # plt.ylim((0,752))
    
    # Ensure the aspect ratio is equal so it looks like a circle
    ax.set_aspect('equal', adjustable='box')
    plt.title(title)
    plt.xlabel("X position")
    plt.ylabel("Y position")
    plt.axis("equal")  # keep aspect ratio    boxoff()
    boxoff()
    plt.show()   

# W for my pc, Z for surf cloud
#dlc_filePath = 'W:\\vs03.herseninstituut.knaw.nl\\VS03-CSF-1\\Conrad\\Innate_approach\\Data_analysis\\24.35.01\\DLC'
dlc_filePath = 'W:\\Conrad\\Innate_approach\\Data_analysis\\24.35.01\\DLC'
dlc_savePath = 'W:\\Conrad\\Innate_approach\\Data_analysis\\24.35.01\\DLC'
file_list = os.listdir(dlc_savePath)


fp_filePath = 'W:\\Conrad\\Innate_approach\\Data_collection\\24.35.01\\'
fp_savePath = 'W:\\Conrad\\Innate_approach\\Data_analysis\\24.35.01\\'

metadata_filepath = 'W:\\Conrad\\Innate_approach\\Data_analysis\\db_FP.mat'
metadata = scipy.io.loadmat(metadata_filepath, struct_as_record=False, squeeze_me=True)
db = metadata['db']
# filter for freelymovingLaser & temporarily just the days that i analyzed for now

ntFilePth = 'W:\\Conrad\\Innate_approach\\Data_collection\\Neurotar\\'


r_log = pd.read_csv(f"{fp_filePath}\\recordinglog.csv", sep=None, engine="python", encoding='utf-8-sig')

# or use this:
# r_log = pd.read_csv(f"{fp_filePath}\\recordinglog.csv", sep=None, engine="python", encoding='cp1252') 
# r_log.columns = r_log.columns.str.replace('\ufeff', '').str.strip()

if not debug:
    r_log = r_log[r_log['Exp']=='fm laser'].reset_index()
    r_log = r_log[:-1] #temp filter for bad recording session
else: 
    r_log = r_log[r_log['ID']==118580].reset_index()
   


animalIDs = r_log['ID'].unique()
dates = r_log['Date'].unique()

combined_snout_distance = []
ttl_pi_drift = []

# initialize paramters, need to change on a file by file basis
# maybe load param file for immutable?
params = {
    'overhead_camera_distortion' : [320, 340],
    'overhead_camera_distortion_method' : 'fisheye_orthographic',
    'overhead_camera_image_offset' : [-4, -4],
    'overhead_camera_width' : 752,
    'overhead_camera_height' : 582,
    'arena_radius_mm' : 175,
    'overhead_arena_center' : [],
    'overhead_camera_angle' : 0,

    #not sure if i need these
    'picamera_time_multiplier' : 1.0002,
    'laser_time_multiplier' : 1.0002,
    'arena_shape' : 'circular'
    
    }

# first we need to correct fish-eye distrotion in dlc data and convert to real world coordinates
for l in range(len(r_log)):
    
    if l != 47:
        continue
    
    animal_id = str(int(r_log['ID'][l]))
    date = str(r_log['Date'][l]).replace('_', '')
    
    
    # if l == 47:
    #     continue
 
    #skip over analyzed data
    if not reanalyze:
        if f"{r_log['Date'][l]}{animal_id} DLC.pkl" in file_list:
            continue
    
    if r_log['first prey'][l] == "don't analyze":
        continue
    
    if r_log['added to dlc'][l] != 'added' and not pd.isna(r_log['added to dlc'][l]):
        continue
 
    record = [entry for entry in db if str(int(animal_id)) in str(entry.subject)
               and str(r_log['Date'][l]).replace('_','-')[:-1] in str(entry.date)][0]
    
    params['overhead_arena_center'] = record.measures.overhead_arena_center
    
    if '2025' in date or '2024' in date:
        params['overhead_camera_distortion_method'] = 'fisheye_orthographic'
    else:
        params['overhead_camera_distortion_method'] = 'normal'

    
    #accidentally indexed wrong on second batch, this fixes it
    if int(date) >= 20250512:
        index_correction_factor = 0
    else:
        index_correction_factor = 1

    print(f"________________________________________________\nCurrent run: {date}{animal_id}\n------------------------------------------------")
    
    if '2025' in date:
        dlcAnimal = f"{dlc_filePath}\\{animal_id}_{str(r_log['Date'][l]).replace('_', '')}*fmLaserMouseFP*.csv" #old model with fisheye
    else:
        dlcAnimal = f"{dlc_filePath}\\{animal_id}_{str(r_log['Date'][l]).replace('_', '')}*fJul2*.csv" # updated model (no fisheye)


    if "filtered" in dlcAnimal:
        print('Using filtered data! You must reanalyze since youre already filtering in this script')
        
    dlcPrey = f"{dlc_filePath}\\{animal_id}_{str(r_log['Date'][l]).replace('_', '')}*prey*.csv"
    dlcIR = f"{dlc_filePath}\\{animal_id}_{str(r_log['Date'][l]).replace('_', '')}*IR*.csv"
    
    ttlFile = f"{ntFilePth}{animal_id}_{date}_01_ttl"
    ttl = pd.read_csv(ttlFile)
    
    if '2025' not in date: # all recieved triggers is noise
        ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True)
    
    if ttl['Event'][0] == 'Received trigger':
        
        #subtracts neurotar ttl input from first laser ttl output
        ttl['Time'] = (pd.to_datetime(ttl['DateTime']) - pd.to_datetime(ttl['DateTime'][0])).dt.total_seconds() + (ttl['Milliseconds'] / 1000) - ttl['Milliseconds'][0] / 1000
        
        ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True) # added because received trigger sometimes added twice
        ttl = ttl.loc[ttl['Event'] !='sync'].reset_index(drop = True)
        ttl['Time'] = ttl['Time'] - stim_dur/sr
    
    elif ttl['Event'][0] == 'sync':
        ttl['Time'] = (pd.to_datetime(ttl['DateTime']) - pd.to_datetime(ttl['DateTime'][0])).dt.total_seconds() + (ttl['Milliseconds'] / 1000) - ttl['Milliseconds'][0] / 1000
        # ttl = ttl.loc[ttl['Event'] !='Received trigger'].reset_index(drop = True) # added because received trigger sometimes added twice
        ttl = ttl.loc[ttl['Event'] !='sync'].reset_index(drop = True)
        ttl['Time'] = ttl['Time'] - stim_dur/sr
    else:
        print('No sync pulse detected at start, something has gone wrong')
        
    fIdx = np.array([])    # code for removing failed laser presentations
    if isinstance(r_log['Failed'][l], int):
        fIdx = np.array(int(r_log['Failed'][l]))
        # eventTS = eventTS.delete(np.array(int(r_log['Failed'][l])))
    elif isinstance(r_log['Failed'][l], float) and  np.isnan(r_log['Failed'][l]) == False:
        fIdx = np.array(list(map(int, str(r_log['Failed'][l]).split('.'))))
        if fIdx[1] == 0:
            fIdx = np.delete(fIdx, 1)
        # eventTS = eventTS.delete(fIdx)
        
    elif isinstance(r_log['Failed'][l], str) and len(r_log['Failed'][l]) == 1:
        fIdx = np.array(int(r_log['Failed'][l]))
        
    elif isinstance(r_log['Failed'][l], str) and len(r_log['Failed'][l]) > 1:
        fIdx = np.array(list(map(int, str(r_log['Failed'][l]).split('.'))))
    
        
    
    if ('optogenetics' == ttl['Event']).any(): 
        if sum('optogenetics' == ttl['Event']) < 6:
            print('Mix of "p" and "b" laser presses detected, you need to edit something.') # as of 17/7/26 this case not detected :) 
        expand = 3
        
    else:
        expand = 2
        
        
    if fIdx.size == 1:
        expanded_fIdx = [fIdx*expand +offset for offset in range(expand)]
    else: 
        expanded_fIdx = [i + offset for i in fIdx*expand for offset in range(expand)]
    
    ttl = ttl.drop(expanded_fIdx, axis = 0).reset_index(drop = True)    
        
    ttl = ttl['Time'][::expand]           
    

    cam_ttl_file = f"{fp_filePath}{animal_id}\\{animal_id}_{date}_001\\{animal_id}_{date}_001_pioverhead_triggers.csv"
    cam_ttl_file_df = pd.read_csv(cam_ttl_file)
    
    # camttl updated several times, this code catches the first sync pulse
    if 'ttl source' in cam_ttl_file_df.columns:
        
        if 'input type' in cam_ttl_file_df.columns:
            recieve_pulses = cam_ttl_file_df[cam_ttl_file_df['input type'] == 'recieved']

            if '2025' in date: # for when i still used nt program to send initial sync
                if (cam_ttl_file_df['ttl source'] == 'neurotar').any():
                    
                    start_frame = recieve_pulses.loc[
                        recieve_pulses['ttl source'] == 'neurotar',
                        'frame'
                    ].iloc[0]
                else:
                    start_frame = recieve_pulses['frame'][1]
                    
            else: 
                start_frame = recieve_pulses.loc[
                    recieve_pulses['ttl source'] == 'sync_ttl',
                    'frame'
                ].iloc[0]
                
        else: 
            print('input type not found')
            if '2025' in date: # for when i still used nt program to send initial sync
                if (cam_ttl_file_df['ttl source'] == 'neurotar').any():
                    start_frame = cam_ttl_file_df.loc[
                        cam_ttl_file_df['ttl source'] == 'neurotar',
                        'frame'
                    ].iloc[0]
                else:
                    start_frame = cam_ttl_file_df['frame'][1]
                    
            else: 
                start_frame = cam_ttl_file_df.loc[
                    cam_ttl_file_df['ttl source'] == 'sync_ttl',
                    'frame'
                ].iloc[0]
                    
    else:
        start_frame = pd.read_csv(cam_ttl_file).index[1]
    
    df_pathList = [dlcAnimal, dlcPrey, dlcIR]
    df_list = [0, 0, 0]  
    
    
    for i in range(len(df_list)):
        
        df_path = glob.glob(df_pathList[i])
        
        if len(df_path) == 0:
            print('No data found, skipping')
            continue
        
        df = pd.read_csv(df_path[0], header =None, low_memory=False)
        df = df.iloc[:,1:]     # removes pointless 1st column
        df_list[i] = df 
        
        # cropping out parts where animal isn't in frame
        if i == 0:
            dlc_df = pd.read_csv(
                df_path[0],
                header=[0,1,2],  
                index_col=0,
                low_memory=False)
            
            if not pd.isna(r_log['animal hidden frames'][l]): 
                animal_out_of_view_index = list(map(int, r_log['animal hidden frames'][l].split(',')))
                
            else:
                animal_out_of_view_index = []
            
            N = len(dlc_df)
            final_segments = []
            
            cuts = [0] + animal_out_of_view_index + [N]
            
            for i in range(0, len(cuts), 2):
                start = cuts[i]
                end   = cuts[i+1]
            
                if start < end:
                    seg = dlc_df.iloc[start:end]
                    final_segments.append(seg)
            
            # fish-eye correct raw, cropped dlc coordinates for keypoint-moseq
            # Initialize a new DataFrame to hold transformed values
            transformed_segments = []

            for segment in final_segments:

                cor_dlc = pd.DataFrame(index=segment.index)  # preserves frame numbers
            
                bodyparts = segment.columns.get_level_values(1).unique()
                
                # Get the scorer name once
                scorer_name = segment.columns.get_level_values(0)[0]
            
                for bodypart in bodyparts:
            
                    bp_cols = segment.loc[:, (slice(None), bodypart, ['x', 'y', 'likelihood'])]
                    
                    # Separate x and y as arrays
                    overhead_x = pd.to_numeric(bp_cols.xs('x', level=2, axis=1).iloc[:,0], errors='coerce').values
                    overhead_y = pd.to_numeric(bp_cols.xs('y', level=2, axis=1).iloc[:,0], errors='coerce').values
                    
                    likelihood = bp_cols.xs('likelihood', level=2, axis=1).iloc[:,0].values
                    
                    arena_x, arena_y = nt_change_overhead_to_arena_coordinates(
                        overhead_x, overhead_y, params, centered_desired = False, 
                    )
                    
                    # Add columns in x, y, likelihood order
                    cor_dlc[(scorer_name, bodypart, 'x')] = arena_x
                    cor_dlc[(scorer_name, bodypart, 'y')] = arena_y
                    cor_dlc[(scorer_name, bodypart, 'likelihood')] = likelihood
                    
                    
                cor_dlc.columns = pd.MultiIndex.from_tuples(
                cor_dlc.columns, 
                names=['scorer', 'bodyparts', 'coords'])
            
                # Sort columns to ensure proper order (scorer, bodypart, then x/y/likelihood)
                coord_order = ['x', 'y', 'likelihood']
                cor_dlc = cor_dlc.reindex(
                    columns=pd.MultiIndex.from_product(
                        [
                            cor_dlc.columns.get_level_values(0).unique(),
                            cor_dlc.columns.get_level_values(1).unique(), 
                            coord_order
                        ],
                        names=['scorer', 'bodyparts', 'coords']
                    )
                )
                
                transformed_segments.append(cor_dlc)
        
                suffix_list = ['_01', '_02', '_03']
                
            for i, cor_dlc in enumerate(transformed_segments):
                file_name = f'\\{animal_id}_{date}'
                if len(transformed_segments) > 1:
                    file_name = file_name + suffix_list[i]
                kpms_file = dlc_savePath + file_name + '.csv'
                
                cor_dlc.to_csv(kpms_file)              
                    
    
    preyTrial_idx_filled = []
    IRTrial_idx_filled = []
    
  
    if pd.isna(r_log['prey trial times'][l]) and '2025' in date: # in the earlier cases for when r_log got corrupted
        reference_file = f"{dlc_savePath}\\backup\\{r_log['Date'][l]}{animal_id} DLC.pkl" 
        with open(reference_file, 'rb') as f:
            processed_data = pickle.load(f)
        
        old_unclip = int(r_log['first prey'][l]) - 120*sr


        prey_snout_distance = processed_data['data'][3]
        IR_snout_distance = processed_data['data'][4]
    
        # Detect where prey trials start and end
        valid = ~np.isnan(prey_snout_distance) # boolean mask
        prey_trial_starts = np.where(np.diff(valid.astype(int)) == 1)[0] + 1
        prey_trial_ends = np.where(np.diff(valid.astype(int)) == -1)[0] + 1 - index_correction_factor
        
        if l == 14:
            prey_trial_starts = prey_trial_starts[1:] # error where trial occured but not ttl registered by rwd. should save this trial somehow? 
            prey_trial_ends = prey_trial_ends[1:] # error where trial occured but not ttl registered by rwd. should save this trial somehow? 
                     
        preyTrial_idx = np.sort(np.concatenate((prey_trial_starts, prey_trial_ends))) + old_unclip
           
        # Detect where IR trials start and end
        valid = ~np.isnan(IR_snout_distance) # boolean mask
        IR_trial_starts = np.where(np.diff(valid.astype(int)) == 1)[0] + 1
        IR_trial_ends   = np.where(np.diff(valid.astype(int)) == -1)[0] + 1 - index_correction_factor
        IRTrial_idx = np.sort(np.concatenate((IR_trial_starts, IR_trial_ends))) + old_unclip
        
        
    
    else:    
        preyTrial_times = (r_log['prey trial times'][l]).split(',') 
        preyTrial_idx = [timestamp_to_frame(ts) for ts in preyTrial_times]

        IRTrial_times = (r_log['IR trial times'][l]).split(',') 
        IRTrial_idx = [timestamp_to_frame(ts) for ts in IRTrial_times]
        
     
    for i in range(0, len(preyTrial_idx), 2):
        start, end = preyTrial_idx[i], preyTrial_idx[i+1]
        preyTrial_idx_filled.extend(range(start, end + index_correction_factor))    
    
    for i in range(0, len(IRTrial_idx), 2):
        start, end = IRTrial_idx[i], IRTrial_idx[i+1]
        IRTrial_idx_filled.extend(range(start, end + index_correction_factor))
    
            
    
    # check it
    if len(preyTrial_idx_filled) != len(df_list[1])-3:
        # raise ValueError(f"Index length {len(preyTrial_idx_filled)} does not match prey length {len(df_list[1])}")
        print(f"Index length {len(preyTrial_idx_filled)} does not match prey length {len(df_list[1])-3}")
        print('Either DLC dropped frames (rerun, unlikely) or you need to check video times in recording log')
        continue
    

    if len(IRTrial_idx_filled) != len(df_list[2])-3:
        # raise ValueError(f"Index length {len(IRTrial_idx_filled)} does not match IR length {len(df_list[2])}")
        print(f"Index length {len(IRTrial_idx_filled)} does not match IR length {len(df_list[2])-3}")
        print('Either DLC dropped frames (rerun) or you need to check video times in recording log')
        continue
        
    if len(preyTrial_idx)%2 == 1 or len(IRTrial_idx)%2 == 1:
        print('trial start or end missing')
        

    
    raw_dlc = np.zeros([len(df_list[0]), 
                       np.shape(df_list[0])[1] + 
                       np.shape(df_list[1])[1] +
                       np.shape(df_list[2])[1]])
    
    # fisheye correct and merge laser position with animal position
    for i in range(len(df_list)):
        raw_dlc = df_list[i]
        
        raw_dlc.columns = pd.MultiIndex.from_arrays([raw_dlc.iloc[1], raw_dlc.iloc[2]])
        raw_dlc = raw_dlc[3:].reset_index(drop = True) # clips nonnumeric info 

        # general data clipping:
        if i == 0: 
            raw_dlc = raw_dlc[start_frame:].reset_index(drop = True) # clips to syncing start pulse
            # if l == 47:
            #     trim_start = preyTrial_idx_filled[0] - start_frame - 20*sr
            # else:
            trim_start = preyTrial_idx_filled[0] - start_frame - setUp*sr
            trim_factor = preyTrial_idx_filled[0] - setUp*sr

            if trim_start < 0:   
                trim_start = preyTrial_idx_filled[0] - start_frame - 6*sr
                trim_factor = preyTrial_idx_filled[0] -6*sr
                 
            trim_end = len(raw_dlc)-5*sr
            raw_dlc = raw_dlc[trim_start:trim_end].reset_index(drop = True)

        if i == 1:
            raw_dlc.index = [x - trim_factor for x in preyTrial_idx_filled]
        elif i == 2:
            raw_dlc.index = [x - trim_factor for x in IRTrial_idx_filled]
        
                
        # for trial info, uncropped
        # Initialize a new DataFrame to hold transformed values
        cor_dlc = pd.DataFrame(index=raw_dlc.index) 
    
        for bodypart in raw_dlc.columns.get_level_values(1).unique():
            if bodypart in ['pawFL', 'pawBL', 'pawFR', 'pawBR']: # these body parts arent often visible overhead, not good data
                continue
            
            try:
                # Get overhead coordinates
                overhead_x = pd.to_numeric(raw_dlc[bodypart]['x'], errors='coerce').values
                overhead_y = pd.to_numeric(raw_dlc[bodypart]['y'], errors='coerce').values
                        
                # Transform coordinates
                arena_x, arena_y = nt_change_overhead_to_arena_coordinates(overhead_x, overhead_y, params, centered_desired = True)
               
                
                # probably should plot body part data   

                
                # filter out low likelihood
 
                ###
                # i need to add something that detects long periods of low likelihood and flags it, eg for periods where i take animal out to fix cables
                ###
                
                low_prob = np.array(raw_dlc[bodypart]['likelihood'].values, dtype=float) < 0.5
                arena_x[low_prob] = np.nan
                arena_y[low_prob] = np.nan
                
                # num_nans_added = np.sum(low_prob)
                # print(f"Number of low probs: {num_nans_added}")
                
                
                
                # filter out out of arena, with some leeway (6cm)
                dist_from_center = np.sqrt(arena_x**2 + 
                                           arena_y**2)
                
                outside_circle = dist_from_center > params['arena_radius_mm'] + 60 # 6cm, changed from 3cm 3/7/26 because sometimes arena got moved during session
                
                # frames = np.arange(len(outside_circle))

                # plt.figure(figsize=(12, 2))
                # plt.scatter(frames[outside_circle],
                #             np.ones(outside_circle.sum()),
                #             s=2)
                # plt.yticks([])
                # plt.xlabel("Frame")
                # plt.title("Frames outside the arena")
                
                # num_nans_added = np.sum(outside_circle)
                # print(f"Number of out of arena: {num_nans_added}")
                
                
                
                arena_x[outside_circle] = np.nan
                arena_y[outside_circle] = np.nan
                
                # filter out jumps 
                dx = np.diff(arena_x)
                dy = np.diff(arena_y)
                distances = np.sqrt(dx**2 + dy**2)
                bad_steps = distances > 70 # roughly equates to 2 meters per second
                mask = np.zeros_like(arena_x, dtype=bool)
                mask[1:] = bad_steps  # offset by one since diff reduces length by 1
                
                # num_nans_added = np.sum(mask)
                # print(f"Number of jumps: {num_nans_added}")
                
                # Replace both x and y values with NaN where the jump is too large
                arena_x[mask] = np.nan
                arena_y[mask] = np.nan
                
                # interpolate
            
                
                ts= np.arange(len(arena_x))
                mask = ~np.isnan(arena_x)
                arena_x_func = pchip(ts[mask], arena_x[mask], extrapolate = True)
                arena_x = arena_x_func(ts)
          
                arena_y_func = pchip(ts[mask], arena_y[mask], extrapolate = True)
                arena_y = arena_y_func(ts)


                # Add transformed and filtered coordinates
                cor_dlc[(bodypart, 'x')] = arena_x
                cor_dlc[(bodypart, 'y')] = arena_y
                cor_dlc[(bodypart, 'likelihood')] = raw_dlc[bodypart]['likelihood'].values
                
        
            except KeyError:
                print(f"Skipping {bodypart}: missing expected columns")
                
        df_list[i] = cor_dlc
    
    
    
    # flag periods where animal is likely out of frame (due to untangling of fibers, etc)
    likelihood_threshold = 0.8
    min_out_of_view_duration = 30  # seconds
    min_out_of_view_frames = min_out_of_view_duration * sr
    
    snout_prob = pd.to_numeric(df_list[0][('snout', 'likelihood')], errors='coerce').values
    hrC_prob   = pd.to_numeric(df_list[0][('hrC', 'likelihood')], errors='coerce').values
    tailBase_prob   = pd.to_numeric(df_list[0][('tailBase', 'likelihood')], errors='coerce').values
    
    # Frame is suspicious if all three keypoints have low likelihood
    # probably_out_of_view = (
    #     (snout_prob < likelihood_threshold) &
    #     (hrC_prob < likelihood_threshold) &
    #     (tailBase_prob < likelihood_threshold)
    # )
    
    # plt.figure()
    # plt.plot(list(range(0,len(snout_prob))), snout_prob)
    # plt.plot(list(range(0,len(snout_prob))), hrC_prob)
    # plt.plot(list(range(0,len(snout_prob))), tailBase_prob)
    # plt.show()
    
    probably_out_of_view = (
        (np.vstack([snout_prob, hrC_prob, tailBase_prob]) < likelihood_threshold)
        .sum(axis=0)
        >= 2
    )
    
    # Find starts and ends of continuous low-likelihood periods
    changes = np.diff(probably_out_of_view.astype(int))
    
    starts = np.where(changes == 1)[0] + 1
    ends   = np.where(changes == -1)[0] + 1
    
    # Handle a run that starts at frame 0
    if probably_out_of_view[0]:
        starts = np.r_[0, starts]
    
    # Handle a run that continues to the final frame
    if probably_out_of_view[-1]:
        ends = np.r_[ends, len(probably_out_of_view)]
    
    # Keep only periods >= 30 seconds
    out_of_view_periods = [
        (start, end)
        for start, end in zip(starts, ends)
        if end - start >= min_out_of_view_frames
    ]
    
    print(f'{len(out_of_view_periods)} out of view periods detected')
       


         
    # calculate snout to laser distance (s2l) for prey    
    snout_speed = speed_calculator(df_list[0][('snout','x')], df_list[0][('snout','y')])
    hrC_speed = speed_calculator(df_list[0][('hrC','x')], df_list[0][('hrC','y')])
    tail_speed = speed_calculator(df_list[0][('tailBase','x')], df_list[0][('tailBase','y')])
    
    
    # snout to prey laser distance
    snout_pos =  df_list[0][[('snout','x'), ('snout','y')]]
    preyLaser_pos = df_list[1][[('preyLaser','x'), ('preyLaser','y')]]
    d_prey = find_distance(snout_pos, preyLaser_pos)
    
    # and for IR...
    IRLaser_pos = df_list[2][[('IR','x'), ('IR','y')]]
    d_IR = find_distance(snout_pos, IRLaser_pos)

    distances_trials = [d_prey, d_IR]  

    # calculate angle from head to snout 
    # create triangle
    
    hrL_pos = df_list[0][[('hrL','x'), ('hrL','y')]]
    hrR_pos = df_list[0][[('hrR','x'), ('hrR','y')]]
    
    earL_pos = df_list[0][[('earL','x'), ('earL','y')]]
    earR_pos = df_list[0][[('earR','x'), ('earR','y')]]

    hr_midpoint = pd.DataFrame((hrL_pos.values + hrR_pos.values) / 2)

    # my_heatplot(snout_pos.iloc[30:, :], 'snout_pos')
    # my_heatplot(hrL_pos.iloc[30:, :], 'hrL position')
    # my_trajectory(hr_midpoint.iloc[30:, :], 'midpoint position')
    
    if plot_traj:
        my_trajectory(snout_pos.iloc[30:, :], f' {date} {animal_id} snout position')

        
    snout2left_distance = find_distance(hrL_pos, snout_pos)
    snout2right_distance = find_distance(hrR_pos, snout_pos)
    hr_dis = find_distance(hrL_pos, hrR_pos)
    ear_dis = find_distance(earL_pos, earR_pos)
    
    # angle calculations
    hrB_pos =  df_list[0][[('hrB','x'), ('hrB','y')]]
    back1_pos =  df_list[0][[('back1','x'), ('back1','y')]]

    prey_orientation = calculate_angle(snout_pos, hrB_pos, preyLaser_pos)
    prey_bouts = plot_angles_trials(prey_orientation, False, plot_angles)
    
    IR_orientation = calculate_angle(snout_pos, hrB_pos, IRLaser_pos)
    IR_bouts = plot_angles_trials(IR_orientation, False, plot_angles)

    relative_head_angle = calculate_angle(snout_pos, hrB_pos, back1_pos)
    relative_head_angle['angle_deg'] = (relative_head_angle["angle_deg"] % 360) - 180
    plot_angles_trials(relative_head_angle, True, plot_angles, prey_bouts)
    plot_angles_trials(relative_head_angle, True, plot_angles, IR_bouts)
    
    head_angles = {
        'rel_to_prey': prey_orientation,
        'rel_to_IR': IR_orientation,
        'rel_to_self': relative_head_angle}


    if (snout2left_distance > 30).any() or (snout2right_distance > 30).any():
        print("snout distance too far, you need to correct it post interpolation")
        
    # create vector that passes through snout and mid point between ears
    
    
    # calculate angle between snout vector and vector to laser
    
    for t in range(2):
        distance = distances_trials[t]
        
        distance = np.array(distance) # need to convert for boolean mask in next step
        valid = ~np.isnan(distance)
        
        # Detect where trials start and end
        trial_starts = np.where(np.diff(valid.astype(int)) == 1)[0] + 1
        trial_ends   = np.where(np.diff(valid.astype(int)) == -1)[0] + 1
        
        misaligned_trials = []
        
        if t == 0:
            # check how off manually elected trial times are off from recorded tll pulses
            # so far i can think of 2 reasons for drift. case 1: frame rate uneven/dropped, case 2: laser started out of frame or took time to 'warm up'. in this case, values in ttl_check will be positive
            predicted = ttl.reset_index(drop = True)*30
            
            ttl_file_diff = np.array(np.diff(ttl)*30, dtype = int)
            DLC_diff = np.diff(trial_starts)
            ttl_check = ttl_file_diff - DLC_diff
            
            ttl_pi_drift.append(ttl_check[-1]-ttl_check[0]) 
            
            # i dont use ahnchor point for anyhting other than chekcing?
            anchor_point = np.where(abs(ttl_check) == np.min(abs(ttl_check)))[0][0] # this is where prey laser has same start time for both cam and sync pi systems, and should be used later for syncing rwd data
            print(f'Anchor point is {np.min(abs(ttl_check))} frames off')
            for index, value in enumerate(ttl_check):
                # ttl_correction.append(value)
                if abs(value) > 3:
                    print(f"Mismatch of {value} frames detected! Trial {index}")
                    misaligned_trials.append([index, value])
                    # in the case of a mismatched trial due to case 1: 
                    # case 2: using an anchor point 
                    
            
            # collect snout2prey laser distances per trial
            for i, (start, end) in enumerate(zip(trial_starts, trial_ends), 1):
                if end-start < 500: # all approach trials, edit this later
                   combined_snout_distance.append(distance[start:end])
                  
        
        # Handle case where data starts or ends with a trial
        if valid[0]:
            trial_starts = np.r_[0, trial_starts]
        if valid[-1]:
            trial_ends = np.r_[trial_ends, len(distance)]
        
        # Plot
        if plot_need:
             plt.figure(figsize=(12, 6))
         
             n_trials = len(trial_starts)
         
             # Choose colormap based on trial type
             if t == 0:  # prey trials
                 cmap = plt.cm.Greens
             else:       # IR trials
                 cmap = plt.cm.YlGnBu  # teal/blue-green
         
             # Generate gradient colors (dark → light)
             colors = cmap(np.linspace(0.9, 0.3, n_trials))
         
             for i, ((start, end), color) in enumerate(zip(zip(trial_starts, trial_ends), colors), 1):
                 if end - start > 600:  # cap long IR trials
                     end = start + 600
         
                 trial_dist = distance[start:end]
                 trial_time = np.arange(len(trial_dist))
         
                 plt.plot(trial_time, trial_dist, color=color, linewidth=2,
                          label=f"Trial {i}")
                 plt.plot(len(trial_time), trial_dist[-1],
                          marker='D', color='black', markersize=6)
         
             plt.xlabel("Time (frames)")
             plt.ylabel("Distance (mm)")
         
             if t == 0:
                 plt.title("Snout to Prey Laser Distance per Trial\n"
                           f"{r_log['Date'][l]}{animal_id}")
             else:
                 plt.title("Snout to IR Laser Distance per Trial\n"
                           f"{r_log['Date'][l]}{animal_id}")
         
             plt.legend()
             boxoff()
             plt.show()

    
            # snout speed
            # plt.figure(figsize=(12, 6))
            # for i, (start, end) in enumerate(zip(trial_starts, trial_ends), 1):
            #     if end-start > 600: # in case of IR trials lasting longer than 20 seconds
            #         end = start + 600
            #     trial_speed = snout_speed[start:end] # may be off by 1 frame
            #     trial_time = np.arange(len(trial_speed))  # trial-relative time
            #     plt.plot(trial_time, trial_speed, label=f"Trial {i}")
            #     plt.plot(len(trial_time), trial_speed[-1], marker = 'D', color = 'k')
            
            # plt.xlabel("Time (frames)")
            # plt.ylabel("Speed (m/s)")
            # plt.legend()
            # if t == 0:
            #     plt.title("Snout Speed per Prey Trial\n  " +
            #               r_log['Date'][l] + animal_id)
            # else: 
            #     plt.title("Snout Speed per IR Trial\n  " + 
            #               r_log['Date'][l] + animal_id)
    
            
            # plt.show()
            
            
    df_list.append(d_prey)
    df_list.append(d_IR)
    df_list.append(snout_speed)
    df_list.append(hrC_speed)
    df_list.append(tail_speed)
    df_list.append(misaligned_trials)
    df_list.append(head_angles)
    df_list.append(out_of_view_periods)

    dlcDat = {
        'session': r_log['Date'][l],
        'mouse': animal_id,    
        
        'data': df_list
        
        # 'site': site,
            }


    save_file = f"{dlc_savePath}\\{r_log['Date'][l]}{animal_id} DLC.pkl"
    pd.to_pickle(dlcDat, save_file) 
    