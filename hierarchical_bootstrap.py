# -*- coding: utf-8 -*-
"""
Hierarchical bootstrap of time-series calcium data

Hierarchy:
    mouse -> session -> trials

Site handling:
    PAG-L + PAG-R -> PAG
    SC-L  + SC-R  -> SC
    etc.

Important:
    - Only sessions with valid trials for the selected signal are included.
    - Mice with no valid sessions for the selected signal are excluded.
    - Trials are sampled within sessions.
    - Sessions are sampled within mice.
    - Mice are sampled at the outermost level.
    - Confidence intervals are calculated independently at every time point.
    - Significance is defined as the bootstrap CI being entirely above
      or below zero for at least `consecutive_points` consecutive samples.

The bootstrap gives equal weight to mice. This prevents mice with more
sessions/trials from automatically contributing more to the population
estimate.
"""

import os
import glob
import pickle
import re

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# USER SETTINGS
# ============================================================

# ------------------------------------------------------------
# Data location
# ------------------------------------------------------------

data_folder = r'YOUR_FOLDER_HERE'

# Example:
data_folder = r'\\vs03.herseninstituut.knaw.nl\VS03-CSF-1\Conrad\Innate_approach\Data_analysis\24.35.01\freelymoving'


# ------------------------------------------------------------
# Signal to analyze
# ------------------------------------------------------------

signal_name = 'ZdFoFApproach'

# Examples:
#
# 'ZdFoFApproach_trialOnset'
# 'ZdFoFAvoid_trialOnset'
# 'ZdFoFNR'
# 'ZdFoFNR_yoked'
# 'ZdFoFITI'
# 'ZdFoFApproach'
# 'ZdFoFAvoid'
# 'IR_ZdFoFApproach'
# 'angular_ZdFoFApproach'


# ------------------------------------------------------------
# Site selection
# ------------------------------------------------------------

# None:
#     plot all anatomical site groups separately
#
# 'PAG':
#     only analyze PAG-L and PAG-R
#
# 'SC':
#     only analyze SC-L and SC-R

site_group = None


# ------------------------------------------------------------
# Hemisphere handling
# ------------------------------------------------------------

# True:
#     PAG-L and PAG-R are pooled into PAG
#
# False:
#     PAG-L and PAG-R remain separate sites

combine_hemispheres = True


# ------------------------------------------------------------
# Time axis
# ------------------------------------------------------------

sr = 30

pre = 5
post = 5

# Number of time points expected:
n_timepoints = int((pre + post) * sr)


# ------------------------------------------------------------
# Bootstrap settings
# ------------------------------------------------------------

n_boot = 10000

confidence = 95

# Minimum number of consecutive time points for significance
thres = round(1/3*sr)  # Consecutive threshold length


# ------------------------------------------------------------
# Plot settings
# ------------------------------------------------------------

plot_signal = True

fig_width = 5
fig_height = 4

ncols = 3

# Plot the raw trial-level SEM as well as the hierarchical
# bootstrap confidence interval.
plot_sem = True

# Show individual mouse mean traces?
plot_mouse_traces = False

# Show individual session mean traces?
plot_session_traces = False


# ------------------------------------------------------------
# Optional exclusions
# ------------------------------------------------------------

exclude_mice = []

exclude_sessions = []


# ============================================================
# HELPER FUNCTIONS
# ============================================================


def get_site_group(site, combine_hemispheres=True):
    """
    Convert a site label such as PAG-L / PAG-R into a canonical
    anatomical site group.

    Examples:
        PAG-L -> PAG
        PAG-R -> PAG
        SC-L  -> SC
        SC-R  -> SC

    If combine_hemispheres=False, the original site name is retained.
    """

    if site is None:
        return None

    site = str(site)

    if not combine_hemispheres:
        return site

    # Remove only a terminal -L or -R
    site_group = re.sub(r'-(L|R)$', '', site)

    return site_group


def get_hemisphere(site):
    """
    Extract hemisphere from a site label.

    Returns:
        'L'
        'R'
        None
    """

    if site is None:
        return None

    match = re.search(r'-(L|R)$', str(site))

    if match is None:
        return None

    return match.group(1)


def get_valid_trials(sesdat, signal_name):
    """
    Extract valid trial x time data from one session.

    A session is considered valid only if:

        - signal exists
        - signal is a numpy array
        - signal is 2-dimensional
        - there is at least one trial
        - at least one trial contains no NaNs

    Trials containing NaNs are removed.

    Returns:
        ndarray of shape (n_trials, n_timepoints)
        or None if no valid trials exist.
    """

    if signal_name not in sesdat:
        return None

    data = sesdat[signal_name]

    # Missing signals are stored as scalar np.nan
    if not isinstance(data, np.ndarray):
        return None

    if data.ndim != 2:
        return None

    if data.shape[0] == 0:
        return None

    # Remove trials containing NaNs
    valid_trials = ~np.isnan(data).any(axis=1)

    data = data[valid_trials]

    if len(data) == 0:
        return None

    return data


def load_session_files(data_folder):
    """
    Find all pickle files in the specified folder.

    Returns a list of file paths.
    """

    files = glob.glob(os.path.join(data_folder, '*.pkl'))

    # Avoid accidentally loading a previously generated result file
    files = [
        f for f in files
        if not os.path.basename(f).startswith('hierarchical_bootstrap_')
    ]

    return sorted(files)


def build_hierarchy_old(
        files,
        signal_name,
        combine_hemispheres=True,
        site_group_filter=None,
        exclude_mice=None,
        exclude_sessions=None):
    """
    Construct:

        site_group
            mouse
                session
                    trials

    Only sessions with valid trials for signal_name are included.

    Returns
    -------
    hierarchy : dict

        hierarchy[site_group][mouse][session_id] = trial x time array

    metadata : list
        One entry per valid session.
    """

    if exclude_mice is None:
        exclude_mice = []

    if exclude_sessions is None:
        exclude_sessions = []

    hierarchy = {}
    metadata = []

    for file_path in files:

        try:
            with open(file_path, 'rb') as f:
                sesdat = pickle.load(f)

        except Exception as e:
            print(f'Could not load: {file_path}')
            print(f'    {e}')
            continue

        # ----------------------------------------------------
        # Required metadata
        # ----------------------------------------------------

        mouse = sesdat.get('mouse', None)
        session = sesdat.get('session', None)
        site = sesdat.get('site', None)

        if mouse is None or session is None or site is None:
            print(f'Skipping file with missing metadata: {file_path}')
            continue

        mouse = str(mouse)
        session = str(session)
        site = str(site)

        # ----------------------------------------------------
        # Exclusions
        # ----------------------------------------------------

        if mouse in exclude_mice:
            continue

        if session in exclude_sessions:
            continue

        # ----------------------------------------------------
        # Canonical anatomical site
        # ----------------------------------------------------

        this_site_group = get_site_group(
            site,
            combine_hemispheres=combine_hemispheres
        )

        # ----------------------------------------------------
        # Site filter
        # ----------------------------------------------------

        if (
            site_group_filter is not None
            and this_site_group != site_group_filter
        ):
            continue

        # ----------------------------------------------------
        # Extract valid trials
        # ----------------------------------------------------

        data = get_valid_trials(
            sesdat,
            signal_name
        )

        # THIS IS THE IMPORTANT FILTER:
        #
        # A session with no valid trials is never added to the
        # hierarchy and therefore can never be resampled.

        if data is None:
            continue

        # ----------------------------------------------------
        # Check time dimension
        # ----------------------------------------------------
        
        data 
        
        if data.shape[1] != n_timepoints:
            raise ValueError(
                f'\nUnexpected number of time points in:\n'
                f'{file_path}\n'
                f'Site: {site}\n'
                f'Mouse: {mouse}\n'
                f'Session: {session}\n'
                f'Expected: {n_timepoints}\n'
                f'Found: {data.shape[1]}'
            )

        # ----------------------------------------------------
        # Add to hierarchy
        # ----------------------------------------------------

        if this_site_group not in hierarchy:
            hierarchy[this_site_group] = {}

        if mouse not in hierarchy[this_site_group]:
            hierarchy[this_site_group][mouse] = {}

        # Prevent accidental overwriting
        if session in hierarchy[this_site_group][mouse]:

            raise ValueError(
                f'Duplicate session detected:\n'
                f'site = {this_site_group}\n'
                f'mouse = {mouse}\n'
                f'session = {session}'
            )

        hierarchy[this_site_group][mouse][session] = data

        metadata.append({
            'file': file_path,
            'mouse': mouse,
            'session': session,
            'site': site,
            'site_group': this_site_group,
            'hemisphere': get_hemisphere(site),
            'n_trials': data.shape[0],
        })

    return hierarchy, metadata

def build_hierarchy(
        files,
        signal_name,
        combine_hemispheres=True,
        site_group_filter=None,
        exclude_mice=None,
        exclude_sessions=None):
    """
    Construct:

        site_group
            mouse
                session
                    trials

    When combine_hemispheres=True, recordings from different
    hemispheres belonging to the same mouse + session + anatomical
    site group are pooled into ONE session. not sure if i want this

    Example:

        PAG-L + PAG-R
        mouse 116637
        session 22/10/25

    becomes:

        PAG
            116637
                22/10/25
                    all valid PAG trials

    Thus the bilateral recordings do not create two independent
    session-level observations.

    Only sessions with valid trials for signal_name are included.

    Returns
    -------
    hierarchy : dict

        hierarchy[site_group][mouse][session_id] = trial x time array

    metadata : list
        One entry per resulting session.
    """

    if exclude_mice is None:
        exclude_mice = []

    if exclude_sessions is None:
        exclude_sessions = []

    hierarchy = {}

    # --------------------------------------------------------
    # Keep metadata temporarily at the file/recording level.
    # This allows us to merge multiple hemisphere recordings
    # belonging to the same session.
    # --------------------------------------------------------

    session_metadata = {}

    for file_path in files:

        try:
            with open(file_path, 'rb') as f:
                sesdat = pickle.load(f)

        except Exception as e:
            print(f'Could not load: {file_path}')
            print(f'    {e}')
            continue

        # ----------------------------------------------------
        # Required metadata
        # ----------------------------------------------------

        mouse = sesdat.get('mouse', None)
        session = sesdat.get('session', None)
        site = sesdat.get('site', None)

        if mouse is None or session is None or site is None:
            print(
                f'Skipping file with missing metadata: '
                f'{file_path}'
            )
            continue

        mouse = str(mouse)
        session = str(session)
        site = str(site)

        # ----------------------------------------------------
        # Exclusions
        # ----------------------------------------------------

        if mouse in exclude_mice:
            continue

        if session in exclude_sessions:
            continue

        # ----------------------------------------------------
        # Canonical anatomical site
        # ----------------------------------------------------

        this_site_group = get_site_group(
            site,
            combine_hemispheres=combine_hemispheres
        )

        # ----------------------------------------------------
        # Site filter
        # ----------------------------------------------------

        if (
            site_group_filter is not None
            and this_site_group != site_group_filter
        ):
            continue

        # ----------------------------------------------------
        # Extract valid trials
        # ----------------------------------------------------

        data = get_valid_trials(
            sesdat,
            signal_name
        )

        if data is None:
            continue

        # ----------------------------------------------------
        # Check time dimension
        # ----------------------------------------------------

        if data.shape[1] < n_timepoints:
            print(
                f'Skipping file with too few time points:\n'
                f'    {file_path}\n'
                f'    Site: {site}\n'
                f'    Found: {data.shape[1]}\n'
                f'    Expected at least: {n_timepoints}'
            )
            continue

        data = data[:, :n_timepoints]

        # ----------------------------------------------------
        # Unique session key
        #
        # This deliberately ignores hemisphere when hemispheres
        # are being combined.
        # ----------------------------------------------------

        session_key = (
            this_site_group,
            mouse,
            session
        )

        # ----------------------------------------------------
        # Add or merge session
        # ----------------------------------------------------

        if session_key not in session_metadata:

            session_metadata[session_key] = {
                'data': data,
                'files': [file_path],
                'sites': [site],
                'hemispheres': [get_hemisphere(site)],
                'n_trials_by_site': {
                    site: data.shape[0]
                }
            }

        else:

            existing = session_metadata[session_key]

            # ------------------------------------------------
            # Same site recorded twice
            # ------------------------------------------------
            #
            # This is different from PAG-L + PAG-R.
            # We flag it because silently pooling two files
            # with the same site label may indicate duplicate
            # data.
            # ------------------------------------------------

            if site in existing['sites']:

                raise ValueError(
                    f'Duplicate recording detected:\n'
                    f'  site group = {this_site_group}\n'
                    f'  mouse      = {mouse}\n'
                    f'  session    = {session}\n'
                    f'  site       = {site}\n'
                    f'  existing file = {existing["files"]}\n'
                    f'  new file      = {file_path}'
                )

            # ------------------------------------------------
            # Different hemispheres from same session
            # ------------------------------------------------

            existing['data'] = np.concatenate(
                [
                    existing['data'],
                    data
                ],
                axis=0
            )

            existing['files'].append(file_path)
            existing['sites'].append(site)
            existing['hemispheres'].append(
                get_hemisphere(site)
            )
            existing['n_trials_by_site'][site] = data.shape[0]

    # ========================================================
    # Convert merged sessions into hierarchy
    # ========================================================

    metadata = []

    for (
        site_group,
        mouse,
        session
    ), info in session_metadata.items():

        data = info['data']

        if site_group not in hierarchy:
            hierarchy[site_group] = {}

        if mouse not in hierarchy[site_group]:
            hierarchy[site_group][mouse] = {}

        hierarchy[site_group][mouse][session] = data

        metadata.append({
            'file': info['files'],
            'files': info['files'],
            'mouse': mouse,
            'session': session,
            'site': info['sites'],
            'sites': info['sites'],
            'site_group': site_group,
            'hemisphere': info['hemispheres'],
            'n_trials': data.shape[0],
            'n_trials_by_site': info['n_trials_by_site'],
        })

    return hierarchy, metadata
# ============================================================
# HIERARCHICAL BOOTSTRAP
# ============================================================


def hierarchical_bootstrap(
        hierarchy,
        n_boot=10000,
        confidence=95,
        random_seed=None):
    """
    Hierarchical bootstrap.

    Hierarchy:

        mouse
            -> sessions
                -> trials

    The number of mice is fixed to the number of observed mice.

    For each bootstrap replicate:

        1. Sample mice with replacement.
        2. For each selected mouse, sample the same number of
           sessions as that mouse originally contains.
        3. For each selected session, sample the same number of
           trials as that session originally contains.
        4. Calculate the mean trace for that mouse.
        5. Average the mouse means.

    Thus mice have equal weight in the final population estimate.

    Returns
    -------

    mean_trace
        Population mean trace.

    lower_ci
        Lower bootstrap confidence bound at every time point.

    upper_ci
        Upper bootstrap confidence bound at every time point.

    bootstrap_means
        Bootstrap population mean traces.
    """

    rng = np.random.default_rng(random_seed)

    mice = list(hierarchy.keys())

    if len(mice) == 0:
        raise ValueError('No valid mice found.')

    n_mice = len(mice)

    bootstrap_means = np.empty(
        (n_boot, n_timepoints),
        dtype=float
    )

    # --------------------------------------------------------
    # Bootstrap
    # --------------------------------------------------------

    for b in range(n_boot):

        # Sample mice with replacement
        sampled_mice = rng.choice(
            mice,
            size=n_mice,
            replace=True
        )

        mouse_traces = []

        for mouse in sampled_mice:

            sessions = hierarchy[mouse]

            session_ids = list(sessions.keys())

            n_sessions = len(session_ids)

            # Sample sessions with replacement
            sampled_sessions = rng.choice(
                session_ids,
                size=n_sessions,
                replace=True
            )

            session_traces = []

            for session_id in sampled_sessions:

                trials = sessions[session_id]

                n_trials = trials.shape[0]

                # Sample trials with replacement
                trial_indices = rng.integers(
                    0,
                    n_trials,
                    size=n_trials
                )

                sampled_trials = trials[trial_indices]

                # Session mean trace
                session_trace = np.mean(
                    sampled_trials,
                    axis=0
                )

                session_traces.append(session_trace)

            # Mouse mean = equal weighting of its sampled sessions
            mouse_trace = np.mean(
                session_traces,
                axis=0
            )

            mouse_traces.append(mouse_trace)

        # Population mean = equal weighting of mice
        bootstrap_means[b] = np.mean(
            mouse_traces,
            axis=0
        )

    # --------------------------------------------------------
    # Bootstrap confidence interval
    # --------------------------------------------------------

    alpha = 100 - confidence

    lower_ci = np.percentile(
        bootstrap_means,
        alpha / 2,
        axis=0
    )

    upper_ci = np.percentile(
        bootstrap_means,
        100 - alpha / 2,
        axis=0
    )

    mean_trace = np.mean(
        bootstrap_means,
        axis=0
    )

    return (
        mean_trace,
        lower_ci,
        upper_ci,
        bootstrap_means
    )


# ============================================================
# CONSECUTIVE SIGNIFICANCE
# ============================================================


def find_consecutive_true(mask, minimum_length):
    """
    Find runs of True values with at least minimum_length samples.

    Returns
    -------
    indices : ndarray
        All indices belonging to qualifying runs.
    """

    mask = np.asarray(mask, dtype=bool)

    if len(mask) == 0:
        return np.array([], dtype=int)

    padded = np.concatenate([
        [False],
        mask,
        [False]
    ])

    changes = np.diff(
        padded.astype(int)
    )

    starts = np.where(changes == 1)[0]
    ends = np.where(changes == -1)[0]

    indices = []

    for start, end in zip(starts, ends):

        if (end - start) >= minimum_length:

            indices.extend(
                range(start, end)
            )

    return np.asarray(
        indices,
        dtype=int
    )


def get_significance_indices(
        lower_ci,
        upper_ci,
        minimum_length=5):
    """
    Determine where the bootstrap CI excludes zero.

    Positive:
        lower CI > 0

    Negative:
        upper CI < 0

    Only runs of minimum_length or longer are retained.
    """

    positive = lower_ci > 0
    negative = upper_ci < 0

    positive_indices = find_consecutive_true(
        positive,
        minimum_length
    )

    negative_indices = find_consecutive_true(
        negative,
        minimum_length
    )

    return positive_indices, negative_indices


# ============================================================
# SUMMARY INFORMATION
# ============================================================


def print_hierarchy_summary(hierarchy, site_name):
    """
    Print the number of mice, sessions and trials contributing
    to a particular site.
    """

    n_mice = len(hierarchy)

    n_sessions = 0
    n_trials = 0

    print()
    print('=' * 70)
    print(f'SITE: {site_name}')
    print('=' * 70)

    for mouse, sessions in hierarchy.items():

        mouse_trials = 0

        for session, data in sessions.items():

            n_sessions += 1
            n_trials += data.shape[0]
            mouse_trials += data.shape[0]

        print(
            f'Mouse {mouse}: '
            f'{len(sessions)} valid sessions, '
            f'{mouse_trials} valid trials'
        )

    print()
    print(f'Mice:    {n_mice}')
    print(f'Sessions: {n_sessions}')
    print(f'Trials:   {n_trials}')
    print('=' * 70)


# ============================================================
# PLOTTING
# ============================================================


def plot_bootstrap_result(
        site_name,
        mean_trace,
        lower_ci,
        upper_ci,
        hierarchy,
        ts,
        signal_name,
        confidence=95,
        consecutive_points=5):

    fig, ax = plt.subplots(
        figsize=(fig_width, fig_height)
    )

    # --------------------------------------------------------
    # Population mean
    # --------------------------------------------------------

    ax.plot(
        ts,
        mean_trace,
        linewidth=2,
        label='Mean'
    )

    # --------------------------------------------------------
    # Hierarchical bootstrap CI
    # --------------------------------------------------------

    ax.fill_between(
        ts,
        lower_ci,
        upper_ci,
        alpha=0.25,
        label=f'{confidence}% hierarchical bCI'
    )

    # --------------------------------------------------------
    # Zero line
    # --------------------------------------------------------

    ax.axhline(
        0,
        linestyle='--',
        color='black',
        linewidth=1
    )

    ax.axvline(
        0,
        linestyle='--',
        color='black',
        linewidth=1.5
    )

    # --------------------------------------------------------
    # Significant positive / negative periods
    # --------------------------------------------------------

    positive_indices, negative_indices = get_significance_indices(
        lower_ci,
        upper_ci,
        minimum_length=consecutive_points
    )

    # Put significance markers slightly above/below the plot
    # using axes coordinates rather than depending on the
    # magnitude of the signal.

    if len(positive_indices) > 0:

        ax.plot(
            ts[positive_indices],
            np.ones(len(positive_indices)) * 0.96,
            's',
            transform=ax.get_xaxis_transform(),
            markersize=4,
            markerfacecolor='black',
            markeredgecolor='black'
        )

    if len(negative_indices) > 0:

        ax.plot(
            ts[negative_indices],
            np.ones(len(negative_indices)) * 0.04,
            's',
            transform=ax.get_xaxis_transform(),
            markersize=4,
            markerfacecolor='black',
            markeredgecolor='black'
        )

    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    ax.set_title(
        f'{site_name}: {signal_name}'
    )

    ax.set_xlabel(
        'Time (s)'
    )

    ax.set_ylabel(
        'Z-Score'
    )

    ax.legend(
        frameon=False
    )

    # --------------------------------------------------------
    # Clean axes
    # --------------------------------------------------------

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

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

    plt.tight_layout()

    plt.show()

    return fig, ax


# ============================================================
# MAIN
# ============================================================


# ------------------------------------------------------------
# Find session files
# ------------------------------------------------------------

files = load_session_files(
    data_folder
)

print(
    f'Found {len(files)} pickle files.'
)


# ------------------------------------------------------------
# Build hierarchy
# ------------------------------------------------------------

hierarchy_all, metadata = build_hierarchy(
    files=files,
    signal_name=signal_name,
    combine_hemispheres=combine_hemispheres,
    site_group_filter=site_group,
    exclude_mice=exclude_mice,
    exclude_sessions=exclude_sessions
)


print()
print(
    f'Valid sessions containing "{signal_name}": '
    f'{len(metadata)}'
)


# ------------------------------------------------------------
# Time axis
# ------------------------------------------------------------

ts = np.linspace(
    -pre,
    post,
    n_timepoints,
    endpoint=False
)


# ------------------------------------------------------------
# Analyze each site
# ------------------------------------------------------------

results = {}


for this_site, mouse_data in hierarchy_all.items():

    print_hierarchy_summary(
        mouse_data,
        this_site
    )

    # Need at least two mice for a meaningful population
    # bootstrap.
    if len(mouse_data) < 2:

        print(
            f'Skipping {this_site}: '
            f'fewer than 2 valid mice.'
        )

        continue

    print(
        f'Bootstrapping {this_site}...'
    )

    (
        mean_trace,
        lower_ci,
        upper_ci,
        bootstrap_means
    ) = hierarchical_bootstrap(
        hierarchy=mouse_data,
        n_boot=n_boot,
        confidence=confidence
    )

    # --------------------------------------------------------
    # Significance
    # --------------------------------------------------------

    positive_indices, negative_indices = (
        get_significance_indices(
            lower_ci,
            upper_ci,
            minimum_length=consecutive_points
        )
    )

    # --------------------------------------------------------
    # Save result in memory
    # --------------------------------------------------------

    results[this_site] = {

        'signal_name': signal_name,

        'mean': mean_trace,

        'lower_ci': lower_ci,

        'upper_ci': upper_ci,

        'bootstrap_means': bootstrap_means,

        'positive_indices': positive_indices,

        'negative_indices': negative_indices,

        'time': ts,

        'hierarchy': mouse_data,
    }

    # --------------------------------------------------------
    # Report significant periods
    # --------------------------------------------------------

    print()

    if len(positive_indices) > 0:

        print(
            f'{this_site}: CI > 0 at '
            f'{len(positive_indices)} time points.'
        )

        print(
            f'First positive time: '
            f'{ts[positive_indices[0]]:.3f} s'
        )

        print(
            f'Last positive time: '
            f'{ts[positive_indices[-1]]:.3f} s'
        )

    else:

        print(
            f'{this_site}: no positive significant period.'
        )

    if len(negative_indices) > 0:

        print(
            f'{this_site}: CI < 0 at '
            f'{len(negative_indices)} time points.'
        )

        print(
            f'First negative time: '
            f'{ts[negative_indices[0]]:.3f} s'
        )

        print(
            f'Last negative time: '
            f'{ts[negative_indices[-1]]:.3f} s'
        )

    else:

        print(
            f'{this_site}: no negative significant period.'
        )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    if plot_signal:

        plot_bootstrap_result(
            site_name=this_site,
            mean_trace=mean_trace,
            lower_ci=lower_ci,
            upper_ci=upper_ci,
            hierarchy=mouse_data,
            ts=ts,
            signal_name=signal_name,
            confidence=confidence,
            consecutive_points=consecutive_points
        )


print()
print('Done.')