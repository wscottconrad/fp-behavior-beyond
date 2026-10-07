# -*- coding: utf-8 -*-
"""
Created on 24/11/25

This code plots summary graphs per animal per site, with the purpose of
data cleaning.

The analysis can now handle multiple experiment types simultaneously.
Select which trial types to analyze using the `trial_type` list.

Example:
    trial_type = ['nt_approach', 'fm_approach']

or:
    trial_type = ['nt_approach']

or:
    trial_type = ['fm_approach']

@author: sconrad
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import os

# =============================================================================
# SETTINGS
# =============================================================================

# Choose one or more trial types to analyze.
#
# Available:
#   'nt_approach'
#   'fm_approach'
#
# You can add other trial types later if needed.
trial_type = ['nt_approach', 'fm_approach']

plot_speed = False

focused_analysis = True

sr = 30
whole_trial = False

initiate_aligned = True

debug = True

# =============================================================================
# FUNCTIONS
# =============================================================================

def build_filename(date, ID, channel):
    return f"{date}{ID} Channel {channel}.pkl"


def get_id_string(ID):
    """
    Convert animal ID to the format used in filenames/dictionaries.

    Handles IDs that are stored as integers or strings.
    """
    try:
        return str(int(ID))
    except (ValueError, TypeError):
        return str(ID)


# =============================================================================
# EXPERIMENT CONFIGURATION
# =============================================================================
#
# Instead of using a global `nt` boolean, each trial type specifies:
#
#   - which experiment it belongs to
#   - where the data are stored
#   - which signal should be extracted
#
# This allows NT and FM data to be loaded and plotted together.
# =============================================================================

filePath = r'W:\Conrad\Innate_approach\Data_collection\24.35.01'

experiments = {

    'nt_approach': {
        'exp': 'nt',

        'tankfolder':
            r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad'
            r'\Innate_approach\Data_analysis\24.35.01\\',

        'signal_key': 'ZdFoFApproach_trialOnset',

        # Signal used if initiate_aligned == True
        'aligned_signal_key': 'ZdFoFApproach',

    },

    'fm_approach': {
        'exp': 'fm laser',

        'tankfolder':
            r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad'
            r'\Innate_approach\Data_analysis\24.35.01\freelymoving\\',

        'signal_key': 'ZdFoFApproach_trialOnset',

        # FM data currently only use IR_ZdFoFApproach regardless of
        # initiate_aligned.
        'aligned_signal_key': 'ZdFoFApproach',
    },
}


# =============================================================================
# CHECK TRIAL TYPES
# =============================================================================

for trial in trial_type:
    if trial not in experiments:
        raise ValueError(
            f"Unknown trial type '{trial}'. "
            f"Available trial types are: {list(experiments.keys())}"
        )


# =============================================================================
# LOAD RECORDING LOG
# =============================================================================

r_log = pd.read_csv(
    f"{filePath}\\recordinglog.csv",
    sep=None,
    engine="python",
    encoding='utf-8-sig'
)

r_log = r_log[r_log['notes'] != 'no ttl alignment']
r_log = r_log[r_log['added to db?'] != 'neurotar data corrupt']

r_log['_ID_str'] = r_log['ID'].apply(get_id_string)

# =============================================================================
# ANIMAL IDS
# =============================================================================

selected_exps = [
    experiments[trial]['exp']
    for trial in trial_type
]

selected_logs = r_log[r_log['Exp'].isin(selected_exps)]

animal_ids = selected_logs['_ID_str'].unique()


# =============================================================================
# SITES
# =============================================================================

sites = pd.unique(
    selected_logs[['1', '2']].values.ravel('K')
)

# Remove NaN sites
sites = [
    site for site in sites
    if not pd.isna(site)
]


excluded_animals = [
    '105647',
    '118401',
    '118402'
]


# =============================================================================
# ANIMAL BEHAVIOR DICTIONARY
# =============================================================================
#
# Example result:
#
# {
#     '105123': {
#         'nt_approach': 25,
#         'fm_approach': 18
#     },
#     '105124': {
#         'nt_approach': 31,
#         'fm_approach': 22
#     }
# }
#
# This is why each animal gets a nested dictionary.
# =============================================================================

animal_behavior = {
    str(animal_id): {}
    for animal_id in animal_ids
}


# =============================================================================
# PARAMETERS
# =============================================================================

thres = 5       # Consecutive threshold length
pre = 5
post = 25


# =============================================================================
# MAIN ANALYSIS LOOP
# =============================================================================

for ID in animal_ids:

    ID = str(ID)

    if ID in excluded_animals:
        continue


    # -------------------------------------------------------------------------
    # Create figure for this animal.
    #
    # There will be one figure per animal/site, with all selected trial types
    # plotted on the same axes.
    # -------------------------------------------------------------------------

    for site in sites:

        if site == 'exclude' or site == 'switched':
            continue

        if focused_analysis:

            if 'to' in str(site):
                continue

            if 'MLR' in str(site):
                continue

        plot_data = []

        plot_labels = []


        # =====================================================================
        # LOOP OVER SELECTED EXPERIMENT/TRIAL TYPES
        # =====================================================================

        for current_trial in trial_type:

            config = experiments[current_trial]

            exp_type = config['exp']
            tankfolder = config['tankfolder']


            # -----------------------------------------------------------------
            # Select recording-log entries for this animal AND experiment
            # -----------------------------------------------------------------

            animal_by_exp = r_log[
                (r_log['_ID_str'] == ID) &
                (r_log['Exp'] == exp_type)
            ]


            if len(animal_by_exp) == 0:

                if debug:
                    print(
                        f"No recording log entries for "
                        f"{ID}, {current_trial}"
                    )

                # Store zero trials in animal_behavior
                animal_behavior[ID][current_trial] = 0

                continue


            # -----------------------------------------------------------------
            # Find recordings from this animal at this site
            # -----------------------------------------------------------------

            animal_by_site_data_ch1 = animal_by_exp[
                animal_by_exp['1'] == site
            ]

            animal_by_site_data_ch2 = animal_by_exp[
                animal_by_exp['2'] == site
            ]


            if (
                len(animal_by_site_data_ch1) == 0
                and len(animal_by_site_data_ch2) == 0
            ):
                if debug:
                    print(
                        f"No {current_trial} recordings for "
                        f"{ID} at {site}"
                    )

                animal_behavior[ID][current_trial] = 0

                continue


            # -----------------------------------------------------------------
            # Get recording dates
            # -----------------------------------------------------------------

            dates_ch1 = animal_by_site_data_ch1['Date'].unique()
            dates_ch2 = animal_by_site_data_ch2['Date'].unique()


            # -----------------------------------------------------------------
            # Load all recordings for this animal/site/experiment
            # -----------------------------------------------------------------

            data_single_animal_site = []


            # Channel 1
            for date in dates_ch1:

                filename = build_filename(date, ID, 1)

                data_path = os.path.join(
                    tankfolder,
                    filename
                )


                if os.path.exists(data_path):

                    if debug:
                        print(
                            f"Loading {current_trial}: "
                            f"{filename}"
                        )

                    with open(data_path, 'rb') as f:
                        temp_data = pickle.load(f)

                    data_single_animal_site.append(temp_data)


                elif debug:
                    print(
                        f"File not found: {data_path}"
                    )


            # Channel 2
            for date in dates_ch2:

                filename = build_filename(date, ID, 2)

                data_path = os.path.join(
                    tankfolder,
                    filename
                )


                if os.path.exists(data_path):

                    if debug:
                        print(
                            f"Loading {current_trial}: "
                            f"{filename}"
                        )

                    with open(data_path, 'rb') as f:
                        temp_data = pickle.load(f)

                    data_single_animal_site.append(temp_data)


                elif debug:
                    print(
                        f"File not found: {data_path}"
                    )


            # -----------------------------------------------------------------
            # Extract signal
            # -----------------------------------------------------------------

            current_data = []


            for item in data_single_animal_site:

               
                if initiate_aligned:
                    signal_key = config['aligned_signal_key']
                else:
                    signal_key = config['signal_key']


                # -------------------------------------------------------------
                # Make sure the signal exists in this recording
                # -------------------------------------------------------------

                if signal_key not in item:

                    if debug:
                        print(
                            f"Missing '{signal_key}' in "
                            f"{current_trial} data for {ID}"
                        )

                    continue


                signal = item[signal_key]


                # -------------------------------------------------------------
                # Skip trials containing NaNs
                # -------------------------------------------------------------

                if not np.any(np.isnan(signal)):
                    current_data.append(signal)


            # -----------------------------------------------------------------
            # Convert to numpy array
            # -----------------------------------------------------------------

            if len(current_data) > 0:

                try:
                    current_data = np.vstack(current_data)

                except ValueError as e:

                    print(
                        f"Could not stack {current_trial} data "
                        f"for {ID}, {site}: {e}"
                    )

                    current_data = np.array([])

            else:

                current_data = np.array([])


            # -----------------------------------------------------------------
            # Save number of trials
            # -----------------------------------------------------------------

            animal_behavior[ID][current_trial] = len(current_data)


            # -----------------------------------------------------------------
            # Store data for plotting
            # -----------------------------------------------------------------

            plot_data.append(current_data)
            plot_labels.append(current_trial)


        # =====================================================================
        # PLOT
        # =====================================================================

        # Only make a figure if there is something to plot.
        if len(plot_data) == 0:
            continue


        plt.figure()


        # ---------------------------------------------------------------------
        # Plot each selected trial type
        # ---------------------------------------------------------------------

        for index, signal in enumerate(plot_data):

            current_trial = plot_labels[index]


            # -------------------------------------------------------------
            # Number of trials
            # -------------------------------------------------------------

            animal_behavior[ID][current_trial] = len(signal)


            if len(signal) == 0:

                print(
                    f"No {current_trial} trials for "
                    f"{ID} {site}\n"
                )

                continue


            print(
                f"{len(signal)} {current_trial} trials\n"
            )


            # -------------------------------------------------------------
            # Colors for plotting
            #
            #   approach -> green
            #   NR       -> grey
            #   ITI/IR   -> blue
            #   avoid    -> red
            # -------------------------------------------------------------

            if current_trial == 'nt_approach':

                plt_color = [0.47, 0.67, 0.19]       # green

            elif current_trial == 'fm_approach':

                plt_color = [1, 0.804, 0.004]     # orange

            elif current_trial == 'approach':

                plt_color = [0.47, 0.67, 0.19]       # green

            elif current_trial == 'NR':

                plt_color = [0.65, 0.65, 0.65]       # grey

            elif current_trial == 'ITI' or current_trial == 'IR':

                plt_color = [0.416, 0.741, 0.741]    # blue

            else:

                plt_color = [0.78, 0, 0]             # red


            # -------------------------------------------------------------
            # Time axis
            # -------------------------------------------------------------

            ts = np.linspace(
                -pre,
                post,
                signal.shape[1]
            )


            # -------------------------------------------------------------
            # Calc mean and SEM
            # -------------------------------------------------------------

            mean_signal = np.mean(
                signal,
                axis=0
            )

            sem_signal = (
                np.std(signal, axis=0)
                / np.sqrt(len(signal))
            )


            # -------------------------------------------------------------
            # Plot mean and SEM
            # -------------------------------------------------------------
            
            n_trials = len(signal)
            
            plt.plot(
                ts,
                mean_signal,
                color=plt_color,
                label=f'{current_trial} (n={n_trials})'
            )


            plt.fill_between(
                ts,
                mean_signal + sem_signal,
                mean_signal - sem_signal,
                color=plt_color,
                alpha=0.3
            )


        # =====================================================================
        # FIGURE FORMATTING
        # =====================================================================

        plt.axvline(
            x=0,
            linestyle='--',
            color='black',
            linewidth=1.5
        )

        plt.axhline(
            y=0,
            linestyle='--',
            color='black',
            linewidth=1.5
        )


        # ---------------------------------------------------------------------
        # Hide top and right spines
        # ---------------------------------------------------------------------

        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)


        # ---------------------------------------------------------------------
        # Tick parameters
        # ---------------------------------------------------------------------

        plt.gca().tick_params(
            axis='x',
            which='both',
            direction='out',
            bottom=True,
            top=False
        )

        plt.gca().tick_params(
            axis='y',
            which='both',
            direction='out',
            left=True,
            right=False
        )


        # ---------------------------------------------------------------------
        # Title and labels
        # ---------------------------------------------------------------------

        plt.title(
            f'{ID}: {site}'
        )

        plt.ylabel(
            'Z-Score'
        )


        # ---------------------------------------------------------------------
        # X limits
        # ---------------------------------------------------------------------

        if not whole_trial:

            plt.xlim(
                -5,
                5
            )


        # ---------------------------------------------------------------------
        # X-axis label
        # ---------------------------------------------------------------------

        if initiate_aligned:

            # If NR is ever added back as a trial type, retain the original
            # shuffled-onset label.

            if 'NR' not in trial_type:

                plt.xlabel(
                    'Movement onset (s)'
                )

            else:

                plt.xlabel(
                    'Shuffled movement onset (s)'
                )

        else:

            plt.xlabel(
                'Prey Laser onset (s)'
            )


        # ---------------------------------------------------------------------
        # Legend
        # ---------------------------------------------------------------------

        plt.legend()

        # ---------------------------------------------------------------------
        # Save figure
        # ---------------------------------------------------------------------

        # plt.savefig(
        #     f'{tankfolder}{site}_{trial_type}.png',
        #     transparent=True
        # )


        plt.show()


    # =========================================================================
    # SPEED ANALYSIS
    # =========================================================================
    #
    # This section is retained for FM data.
    #
    # Because speed information is currently only collected from the FM
    # recordings, this section is run if 'fm_approach' is selected and
    # plot_speed is True.
    # =========================================================================

    if plot_speed and 'fm_approach' in trial_type:

        # ---------------------------------------------------------------------
        # FM recording log for this animal
        # ---------------------------------------------------------------------

        fm_log = r_log[
            (r_log['_ID_str'] == ID) &
            (r_log['Exp'] == 'fm laser')
        ]


        # ---------------------------------------------------------------------
        # Speed containers
        # ---------------------------------------------------------------------

        prey_speeds = {
            'snout': [],
            'hrC': [],
            'tail': []
        }

        IR_speeds = {
            'snout': [],
            'hrC': [],
            'tail': []
        }


        # ---------------------------------------------------------------------
        # Loop over sites
        # ---------------------------------------------------------------------

        for site in sites:

            if site == 'exclude' or site == 'switched':
                continue


            if focused_analysis:

                if 'to' in str(site):
                    continue

                if 'MLR' in str(site):
                    continue


            # -----------------------------------------------------------------
            # Get recordings for this animal/site
            # -----------------------------------------------------------------

            animal_by_site_data_ch1 = fm_log[
                fm_log['1'] == site
            ]

            animal_by_site_data_ch2 = fm_log[
                fm_log['2'] == site
            ]


            if (
                len(animal_by_site_data_ch1) == 0
                and len(animal_by_site_data_ch2) == 0
            ):
                continue


            dates_ch1 = animal_by_site_data_ch1['Date'].unique()
            dates_ch2 = animal_by_site_data_ch2['Date'].unique()


            data_single_animal_site = []


            # -----------------------------------------------------------------
            # Channel 1
            # -----------------------------------------------------------------

            for date in dates_ch1:

                filename = build_filename(
                    date,
                    ID,
                    1
                )

                data_path = os.path.join(
                    experiments['fm_approach']['tankfolder'],
                    filename
                )


                if os.path.exists(data_path):

                    with open(data_path, 'rb') as f:
                        temp_data = pickle.load(f)

                    data_single_animal_site.append(temp_data)


            # -----------------------------------------------------------------
            # Channel 2
            # -----------------------------------------------------------------

            for date in dates_ch2:

                filename = build_filename(
                    date,
                    ID,
                    2
                )

                data_path = os.path.join(
                    experiments['fm_approach']['tankfolder'],
                    filename
                )


                if os.path.exists(data_path):

                    with open(data_path, 'rb') as f:
                        temp_data = pickle.load(f)

                    data_single_animal_site.append(temp_data)


            # -----------------------------------------------------------------
            # Extract speed information
            # -----------------------------------------------------------------

            for item in data_single_animal_site:

                if 'approach_speeds' not in item:
                    continue

                if 'IR_speeds' not in item:
                    continue


                if not pd.isna(item['approach_speeds']):

                    for key, speeds in item['approach_speeds'].items():

                        if key in prey_speeds:
                            prey_speeds[key].append(speeds)


                    for key, speeds in item['IR_speeds'].items():

                        if key in IR_speeds:
                            IR_speeds[key].append(speeds)


        # =====================================================================
        # PLOT SPEED
        # =====================================================================

        ncols = 3
        nrows = 1


        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(5 * ncols, 4 * nrows),
            sharex=True,
            constrained_layout=True
        )


        ts = np.linspace(
            -pre,
            post,
            (pre + post) * sr
        )


        # Flatten axes so you can index them with ax_idx
        axes = np.atleast_1d(axes).ravel()

        ax_idx = 0


        for body_part, speed_data in prey_speeds.items():

            ax = axes[ax_idx]

            ax_idx += 1


            # -----------------------------------------------------------------
            # Prey
            # -----------------------------------------------------------------

            arrays = [
                a
                for a in prey_speeds[body_part]
                if len(a) > 0
            ]


            if len(arrays) == 0:
                continue


            prey_data = np.vstack(arrays)


            # -----------------------------------------------------------------
            # IR
            # -----------------------------------------------------------------

            arrays = [
                a
                for a in IR_speeds[body_part]
                if len(a) > 0
            ]


            if len(arrays) == 0:
                continue


            IR_data = np.vstack(arrays)


            # -----------------------------------------------------------------
            # Prey plot
            # -----------------------------------------------------------------

            ax.plot(
                ts,
                np.mean(prey_data, axis=0),
                color=[0.47, 0.67, 0.19],
                label='Prey'
            )


            ax.fill_between(
                ts,
                np.mean(prey_data, axis=0)
                + np.std(prey_data, axis=0)
                / np.sqrt(len(prey_data)),

                np.mean(prey_data, axis=0)
                - np.std(prey_data, axis=0)
                / np.sqrt(len(prey_data)),

                color=[0.47, 0.67, 0.19],
                alpha=0.3
            )


            # -----------------------------------------------------------------
            # IR plot
            # -----------------------------------------------------------------

            ax.plot(
                ts,
                np.mean(IR_data, axis=0),
                color=[0.416, 0.741, 0.741],
                label='IR'
            )


            ax.fill_between(
                ts,
                np.mean(IR_data, axis=0)
                + np.std(IR_data, axis=0)
                / np.sqrt(len(IR_data)),

                np.mean(IR_data, axis=0)
                - np.std(IR_data, axis=0)
                / np.sqrt(len(IR_data)),

                color=[0.416, 0.741, 0.741],
                alpha=0.3
            )


            # -----------------------------------------------------------------
            # Axis formatting
            # -----------------------------------------------------------------

            ax.set_title(
                f'{ID}: {body_part} speed'
            )


            ax.tick_params(
                axis='x',
                which='both',
                direction='out',
                bottom=True,
                top=False
            )


            ax.tick_params(
                axis='y',
                which='both',
                direction='out',
                left=True,
                right=False
            )


            ax.set_xlim(
                -5,
                5
            )


            ax.axvline(
                x=0,
                linestyle='--',
                color='black',
                linewidth=1.5
            )


            ax.set_xlabel(
                'Movement onset (s)'
            )

            ax.set_ylabel(
                'Speed (m/s)'
            )


        plt.show()

