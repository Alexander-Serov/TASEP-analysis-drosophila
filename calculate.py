

import copy
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import shapiro
from scipy.stats import t as student
from tqdm import tqdm, trange

from constants import (LV, AP_hist_folder, L, default_intensity_threshold,
                       dt_new, gene_labels, intensity_thresholds, k, kV, l, lV,
                       output_slopes_folder, slope_length_mins)
from support import J_over_k_MC, bayesian_linear_fit, mixed_estimator_2
from support import rho_MC as rho_MC_func
from support import set_figure_size, welchs_test


def alpha_over_k_J(JoK, JoKV):
    """
    Calculate mean and variance of alpha/k from J/k
    """
    alpha_over_k = (1
                    - JoK * (l - 1)
                    - np.sqrt(1 + JoK**2 * (l - 1) **
                              2 - 2 * JoK * (l + 1))
                    ) / 2
    partial_JoK = (1 / 2) * (1 - l - (2 * JoK * (-1 + l)**2 - 2 * (1 + l)) /
                             (2 * np.sqrt(1 + JoK ** 2 * (-1 + l)**2 - 2 * JoK * (1 + l))))
    partial_l = (1 / 2) * (-JoK - (-2 * JoK + 2 * JoK**2 * (-1 + l)) /
                           (2 * np.sqrt(1 + JoK**2 * (-1 + l)**2 - 2 * JoK * (1 + l))))
    alpha_over_kV = partial_l**2 * lV + partial_JoK**2 * JoKV
    return alpha_over_k, alpha_over_kV


def alpha_over_k_rho(rho, rhoV):
    """
    Calculate mean and variance of alpha/k from rho
    """
    alpha_over_k = rho / (l - rho * (l - 1))
    partial_l = -(((1 - rho) * rho) / (l - (-1 + l) * rho)**2)
    partial_rho = l / (l + rho - l * rho)**2

    alpha_over_kV = partial_l**2 * lV + partial_rho**2 * rhoV

    return alpha_over_k, alpha_over_kV


def calculate_alpha(analyses_in):
    """
    Calculates alpha/k based on experimental slopes and steady state N.
    For alpha estimates, we use rho and JoK estimates that are already averaged over differnt embryos of the same group.
    Each group corresponds to a different gene-construct-nc combination.
    After that we construct a mixed estimator based on the two following the method described in Lavancier and Rocher (2016).

    The variance of individual alpha estimators is calculated as a sum of squares of variances with partial derivatives:
    V = \sum_i (d alpha / d \theta_i)**2 * var(theta_i)
    """

    analyses = analyses_in.copy()
    significance_level = 0.05

    genes = set(analyses.gene)
    constructs = set(analyses.construct)
    ncs = set(analyses.index.get_level_values(1))

    # Load alphas obtained from J/k
    JoK = analyses.JoK
    JoKV = analyses.JoKV
    analyses['alpha_over_k_J'], analyses['alpha_over_k_JV'] = alpha_over_k_J(JoK, JoKV)
    analyses['alpha_J'] = analyses['alpha_over_k_J'] * k

    # Load alphas obtained from rho
    rho = analyses.rho
    rhoV = analyses.rhoV
    analyses['alpha_over_k_rho'], analyses['alpha_over_k_rhoV'] = alpha_over_k_rho(rho, rhoV)
    analyses['alpha_rho'] = analyses['alpha_over_k_rho'] * k

    # I use aoK to denote alpha/k
    def tau(aoK):
        # Invert alpha to get tau and recalculate into seconds
        tau = 1 / aoK / k * 60
        return tau

    refuted_shapiro_alpha = []
    refuted_shapiro_tau = []
    for gene in genes:
        for construct in constructs:
            indices = (analyses.gene == gene) & (
                analyses.construct == construct)
            for nc in ncs:
                # %% Mixed estimator for alpha/k
                # Drop nans
                T1 = drop_nan(analyses.loc[(indices, nc), 'alpha_over_k_rho'].values)
                T2 = drop_nan(analyses.loc[(indices, nc), 'alpha_over_k_J'].values)

                # Calculate the mixed estimator
                aoK_mixed, aoK_V, weights = mixed_estimator_2(T1=T1, T2=T2)

                # Save to analyses
                analyses.loc[(indices, nc), 'alpha_over_k_comb'] = aoK_mixed
                analyses.loc[(indices, nc), 'alpha_over_k_combV'] = aoK_V
                analyses.loc[(indices, nc),
                             'kappa'] = weights[0]
                analyses.loc[(indices, nc), 'alpha_comb'] = aoK_mixed * k
                analyses.loc[(indices, nc), 'alpha_combV'] = aoK_V * k**2

                # Perform Shapiro-Wilk's normality test to see whether alpha or tau are more appropriate for comparison with other conditions. 1st part
                # Save to analyses for later re-use
                data_mixed = weights[0] * T1 + weights[1] * T2
                if len(data_mixed) >= 3:
                    shap_W, shap_p = shapiro(data_mixed)
                    analyses.loc[(indices, nc), 'shap_alpha_W'] = shap_W
                    analyses.loc[(indices, nc), 'shap_alpha_p'] = shap_p
                    refuted_shapiro_alpha.append(shap_p < significance_level)

                # %% Mixed estimator for tau
                # Calculate, then save to analyses
                T1 = drop_nan(tau(analyses.loc[(indices, nc), 'alpha_over_k_rho'].values))
                T2 = drop_nan(tau(analyses.loc[(indices, nc), 'alpha_over_k_J'].values))
                tau_mixed, tau_V, weights = mixed_estimator_2(T1=T1, T2=T2)
                analyses.loc[(indices, nc), 'tau'] = tau_mixed
                analyses.loc[(indices, nc), 'tauV'] = tau_V
                analyses.loc[(indices, nc),
                             'kappa_tau'] = weights[0]

                # Perform Shapiro-Wilk's normality test for the mixed estimator
                data_mixed = weights[0] * T1 + weights[1] * T2
                if len(data_mixed) >= 3:
                    shap_W, shap_p = shapiro(data_mixed)
                    analyses.loc[(indices, nc), 'shap_tau_W'] = shap_W
                    analyses.loc[(indices, nc), 'shap_tau_p'] = shap_p
                    refuted_shapiro_tau.append(shap_p < significance_level)

                # Print the mixed estimator for the current gene-construct-nc combination
                print(
                    f'{gene}, {construct}, nc{nc}: tau = {tau_mixed:.1f} +- {np.sqrt(tau_V):.1f} s, alpha = {aoK_mixed*k:.1f} +- {np.sqrt(aoK_V*k**2):.1f} min^{-1}')

    print('\nShapiro-Wilk normality test results across all genes, constructs, ncs:')
    print(
        f'For tau, normality refuted in {np.sum(refuted_shapiro_tau)} cases out of {len(refuted_shapiro_tau)} at p = {significance_level}')
    print(
        f'For alpha, normality refuted in {np.sum(refuted_shapiro_alpha)} cases out of {len(refuted_shapiro_alpha)} at p = {significance_level}')

    return analyses


def drop_nan(T):
    return T[~np.isnan(T)]


def calculate_free_travel_time(analyses_in):
    """
    Estimate the free passage time by fitting the non-scaled current-density diagram across all ncs for each gene and construct individually.

    The fitting procedure takes into account uncertainties in data points along both axes.
    the fit line must pass through 0.

    A Gaussian distribution of the slope values is used a prior.
    The width of the distribution `TV_prior` is defined by the a priori error in the elongation rate k and effective gene length L.
    Those are set in `constants.py`
    """

    analyses = analyses_in.copy()
    genes = set(analyses.gene)

    # Add columns to the analyses table
    columns = ['T', 'TV', 'L', 'LV']
    for column in columns:
        analyses[column] = np.nan

    # Make a prior on the travel time
    T_prior = L / k
    TV_prior = T_prior**2 * (LV / L**2 + kV / k**2)

    def prior_func(m, sigmaV):
        """
        The prior must be a function of 2 parameters, so sigmaV in the definition is necessary.
        But this particular prior does not depend on sigmaV
        This is a conjugated prior in the form of just the Gaussian distribution for the slope.
        """
        return np.exp(-(m - T_prior)**2 / 2 / TV_prior) * (2 * np.pi * TV_prior)**(-1 / 2)

    # Fit the diagram individually for each gene, but all constructs and ncs mixed
    # This makes sense because all constructs and nc supposedly have the same length and polymerase elongation rate
    for gene in genes:
        analyses_filtered = analyses[(analyses.gene == gene)]
        indices = analyses_filtered.index

        Ns = analyses_filtered['max'].values
        NsV = analyses_filtered['maxV'].values
        slopes = analyses_filtered['slope'].values
        slopesV = analyses_filtered['slopeV'].values

        # Fit N/s passing through zero. Give result in minutes
        fit = bayesian_linear_fit(x=slopes, Vx=slopesV, y=Ns, Vy=NsV, c=False, prior=prior_func)
        T, TV = [fit[key] for key in ['m', 'mV']]

        # Estimate L from T and k
        L_est = k * T
        LVest = T**2 * kV + k**2 * TV

        # Round up and print
        L_rnd, Lstd_rnd = np.round([L_est, np.sqrt(LVest)], decimals=-1)
        print(f'{gene}: T = {T:.2f} +- {np.sqrt(TV):.2f}, L = {L_rnd:.0f} +- {Lstd_rnd:.0f}')

        # Save to analyses
        analyses.loc[indices, 'T'] = T
        analyses.loc[indices, 'TV'] = TV
        analyses.loc[indices, 'L'] = L_est
        analyses.loc[indices, 'LV'] = LVest

    return analyses


def calculate_rho_and_J(analyses_in, I_est, IV_est=0):
    """
    The function uses the calibration coefficient to convert a.u. into polymerase numbers.
    It calculate:
    - the site occupation density rho,
    - the site density normalized to the maximal current regime r,
    - polymerase flux J/k in pol/min,
    - polymerase flux normalized to the MC regime

    One may provide the error on the calibration coefficient as `IV_est`
    """

    analyses = analyses_in.copy()
    T = analyses['T']   # Cannot use a dot because it is going to transpose
    TV = analyses.TV

    rho_MC, rho_MCV = rho_MC_func()
    JoK_MC, JoK_MCV = J_over_k_MC()

    # rho
    max = analyses['max']
    maxV = analyses['maxV']
    analyses['rho'] = rho = max * l / I_est / k / T
    analyses['rhoV'] = rhoV = rho**2 * (maxV / max**2 + lV / l**2
                                        + IV_est / I_est**2 + kV / k**2 + TV / T**2)

    analyses['r'] = r = rho / rho_MC
    analyses['rV'] = r**2 * (rhoV / rho**2 + rho_MCV / rho_MC**2)

    # J/k
    slope = analyses['slope']
    slopeV = analyses['slopeV']
    analyses['JoK'] = JoK = slope / k / I_est
    analyses['JoKV'] = JoKV = JoK**2 * (slopeV / slope**2 + kV / k**2 + IV_est / I_est**2)

    analyses['j'] = j = JoK / JoK_MC
    analyses['jV'] = j**2 * (JoKV / JoK**2 + JoK_MCV / JoK_MC**2)

    return analyses


def calculate_slopes(data, analyses_in, save_figures, pdf=False):
    """
    Slope and max. polymerase number calculation procedure.
    The characteristics are not calculated if less than 3 frames are present, i.e. nc duration >= 2*dt_new.

    The function plots and outputs slope detection figures for each data set in the `output_slopes_folder`.

    Also produces a histogram of recorded AP positions.
    One AP per trace is used, and the AP is the mean observed AP value.
    """
    # Plot paraemters
    markersize = 2
    lw = 1  # line width
    colors = {'trace': '#68972F',
              'slope': '#ED1C24',
              'threshold': '#EDB120',
              'max': '#0072BD'}
    alpha = 0.2     # transparency
    height_factor = 0.5

    # Constants
    ncs = range(11, 15)

    analyses = analyses_in.copy()
    dataset_ids = set(data.dataset_id)

    for dataset_id in tqdm(dataset_ids, desc='Calculating slopes'):
        dataset_data = data[data.dataset_id == dataset_id]

        dataset_name = dataset_data.dataset.iloc[0]
        dataset_id = dataset_data.dataset_id.iloc[0]
        gene_id = dataset_data.gene_id.iloc[0]
        gene = dataset_data.gene.iloc[0]
        construct_id = dataset_data.construct_id.iloc[0]

        # Check if there are data points after filtering
        if not dataset_data.shape[0]:
            logging.warning(
                f"No data points available for slope fitting after positional fitering. Skipping dataset={dataset_id}")
            continue

        # %% Calculate slopes and max values
        set_figure_size(num=5, rows=1, page_width_frac=0.5,
                        clear=True, height_factor=height_factor)
        fig, ax = plt.subplots(1, 1, num=5)

        # Get an average trace on a regulare time mesh
        avg_data, std_data = get_avg_on_regular_time_mesh(dataset_data, dt=dt_new)

        # Plot the average trace
        x = avg_data.index
        y = avg_data.intensity
        ax.plot(x.values, y.values, '-o', markersize=markersize, lw=lw, c=colors['trace'])

        # Add standard deviation as a shaded region
        ylow = y - std_data.intensity
        yup = y + std_data.intensity
        plt.fill_between(x.values, ylow.values, yup.values,
                         facecolor=colors['trace'], alpha=alpha, edgecolor=None)

        # Detect the time limits of the plot
        t_span = avg_data[~pd.isna(avg_data.intensity)].index.values
        xlims = np.array(t_span[[0, -1]])
        # Extend slightly for better presentation
        xlims = xlims + (xlims[1] - xlims[0]) * 0.1 / 2 * np.array([-1, 1])
        plt.xlim(xlims)

        # %% Calculate slopes and max intensity
        for nc in ncs:
            nc_data = dataset_data[dataset_data.nc == nc]
            if not nc_data.empty:
                # Add columns to analyses
                cols = ['slope', 'slopeV', 'max', 'maxV']
                for col in cols:
                    if col not in analyses:
                        analyses[col] = np.nan

                # Fit slope on a fixed time interval after the start of the nc
                start_slope = nc_data.time.min()
                end_slope = start_slope + slope_length_mins

                # Calculate max and slope only if the nc includes at least 3 frames.
                # Otherwise, return nans
                nc_length = nc_data.time.max() - nc_data.time.min()
                if nc_length >= 2 * dt_new:
                    slope, slope_V, coefs = get_slope(nc_data, start=start_slope, end=end_slope)
                    max_intensity, max_intensity_std = get_max_intensity(nc_data)
                else:
                    slope, slope_V = [np.nan] * 2
                    coefs = [np.nan] * 2
                    max_intensity, max_intensity_std = [np.nan] * 2
                if slope < 0:
                    continue

                # Save to analyses
                analyses.loc[(dataset_id, nc), cols] = (
                    slope, slope_V, max_intensity, max_intensity_std**2)

                # Plot current nc slope and detected max inensity
                plot_slopes(dataset_data, coefs, start_slope,
                            end_slope, max_intensity, nc, ax, lw=lw, colors=colors)

        # Save figure
        figname = "slopes_{id:02d}_{dataset_name}".format(id=dataset_id, dataset_name=dataset_name)
        figpath = os.path.join(output_slopes_folder, figname)

        if save_figures:
            if pdf:
                fig.savefig(figpath + '.pdf', pad_inches=0, bbox_inches='tight')

            # Save an enlarged png version for easier analysis
            factor = 3
            figsize = fig.get_size_inches()
            fig.set_figwidth(figsize[0] * factor)
            fig.set_figheight(figsize[1] * factor)
            fig.savefig(figpath + '.png', pad_inches=0, bbox_inches='tight')

        # %% Plot a histogram of recorded AP positions
        fig = plt.figure(num=10, clear=True)
        str_title = '%s, id=%i' % (dataset_name.replace('_', '\_'), dataset_id)

        MeanAPs = data.groupby(data['trace_id']).ap_mean.first()
        plt.hist(MeanAPs)

        # Adjust the plot
        plt.xlim([0, 1])
        plt.xlabel('AP')
        plt.title(str_title)

        # Save figure
        figpath = f"AP_hist_{dataset_id}.png"
        figpath = os.path.join(AP_hist_folder, figpath)
        if save_figures:
            fig.savefig(figpath)

    return analyses


def get_slope(nc_data, start, end):
    """
    Fit the first several minutes of a nuclear cycle with a straight line.
    Note that all points are fitted and not the average trace to avoid wrong weight assignment.
    A nc must have at least 4 points to be fitted.
    """
    cur_data = nc_data[(nc_data.time >= start)
                       & (nc_data.time <= end)]

    if len(cur_data) <= 3:
        slope = np.nan
        coefs = [np.nan, np.nan]
        slope_V = np.nan
    else:
        coefs, V = np.polyfit(cur_data.time, cur_data.intensity, 1, cov=True)
        slope = coefs[0]
        slope_V = V[0, 0]
    return slope, slope_V, coefs


def get_max_intensity(nc_data):
    """
    Get the maximum polymerase number for a cycle. Average within +- 1 time step of the regular mesh
    """

    if nc_data.empty:
        return np.nan

    # print(nc_data)
    avg_nc_data, _ = get_avg_on_regular_time_mesh(nc_data, dt_new)
    # print(avg_nc_data.intensity.dtype)
    max_time = avg_nc_data.intensity.idxmax()
    # print(max_time)

    # interval = max_time + np.array([-1, 1]) * dt_new
    int = nc_data[(nc_data.time >= max_time - dt_new) &
                  (nc_data.time <= max_time + dt_new)].intensity
    max_intensity = int.mean()
    max_intensity_std = np.std(int, ddof=1)
    return max_intensity, max_intensity_std


def plot_slopes(dataset_data, coefs, start, end, max_intensity, nc, ax, lw, colors):
    """
    Plot the initial slopes and measured max intensity of the nuclear cycle.
    Also plots the noise threshold.
    """
    dataset_name = dataset_data.dataset.iloc[0]
    dataset_id = dataset_data.dataset_id.iloc[0]

    # Plot the slope
    x_fit = np.asarray([start, end])
    y_fit = x_fit * coefs[0] + coefs[1]
    ax.plot(x_fit, y_fit, 'r-', lw=lw, c=colors['slope'])

    # Plot the noise threshold
    xlims = plt.xlim()
    intensity_threshold = intensity_thresholds.get(dataset_id, default_intensity_threshold)
    ax.plot(xlims, [intensity_threshold] * 2, 'g', lw=lw, c=colors['threshold'])

    # Print nc numbers over each detected nc
    _, ymax = plt.ylim()
    plt.text(start, ymax * 0.9, f'nc{nc}', fontsize=7)

    # Plot the max. polymerase number
    x_max = [dataset_data[dataset_data.nc == nc].time.min(),
             dataset_data[dataset_data.nc == nc].time.max()]
    y_max = np.asarray([1, 1]) * max_intensity
    plt.plot(x_max, y_max, 'm--', lw=lw, c=colors['max'])

    # Add plot title and adjust plot
    str_title = dataset_name.replace('_', '\_')
    plt.title(str_title)
    plt.xlabel('$t$, min')
    plt.ylabel('Fluo, a.u.')
    plt.ylim(ymin=-50)

    return


def detect_timestep(data):
    """
    Get the time steps present in the loaded data
    """
    same_trace = data.shift(-1).trace_id == data.trace_id
    dts = data.shift(-1).time - data.time
    dts[~same_trace] = np.nan

    return dts.mean(), dts.min(), dts.max()


def get_regular_time_mesh(data, dt):
    """
    Provide a new regular time mesh with a given dt and a time point at t=0
    """
    t_min = data.time.min()
    t_max = data.time.max()

    if t_max > 0 and t_min < 0:
        time_mesh_positive = np.arange(dt, t_max + dt, dt)
        time_mesh_negative = (-np.arange(dt, -t_min + dt, dt))[::-1]
        time_mesh = np.concatenate((time_mesh_negative, [0], time_mesh_positive))
    elif t_max > 0 and t_min >= 0:
        time_mesh = np.arange(t_min, t_max + dt, dt)
    elif t_max <= 0 and t_min < 0:
        time_mesh = (-np.arange(-t_max, -t_min + dt, dt))[::-1]

    return time_mesh


def get_avg_on_regular_time_mesh(data, dt):
    """
    Calculate an average time trace on a regular time mesh with a given dt
    """
    time_mesh = get_regular_time_mesh(data, dt)
    avg_data = pd.DataFrame(index=time_mesh[:-1], columns=['intensity'], dtype=float)
    std_data = avg_data.copy()

    for i in range(len(time_mesh) - 1):
        int = data[(data.time >= time_mesh[i])
                   & (data.time < time_mesh[i + 1])].intensity
        avg_data.loc[time_mesh[i]] = int.mean()
        std_data.loc[time_mesh[i]] = np.std(int, ddof=1)

    # Interpolate if values are missing. Do not extrapolate
    avg_data_not_nan = avg_data[~np.isnan(avg_data.intensity)]
    nan_time = avg_data[np.isnan(avg_data.intensity)].index.values
    # print(data)
    # print(avg_data)
    # print(avg_data_not_nan)
    # print(nan_time)

    # Interpolate if more than 1 time point present. Note this is average interpolation
    if avg_data_not_nan.count().intensity > 1:
        interp_func = interp1d(avg_data_not_nan.index.values, avg_data_not_nan.intensity.values)
        avg_data.loc[nan_time, 'intensity'] = interp_func(nan_time)
    return avg_data, std_data


def filter_by_AP(data_in):
    """
    This function contains the code used for filtering fluorescence traces by their AP positions

    For hb, we find the AP corresponding to maximum expression APmax, and then take all traces in the interval [0.05; APmax - 0.05].

    For kn, the expressing AP region moves with time. So instead of doing filtering over the whole trace altogether, we do frame by frame filtering. We do not need the connectivity of the traces, because they are anyway averaged out in a single trace for each data set.
    The kn expressing region is rather small, so the filtering works by finding the median AP (mAP) of all expressing nuclei for each individual frame and then keeping only traces withing mAP +- 0.02.

    For sn, no filtering is performed.
    """

    # Constants
    knirps_half_width_AP = 0.02
    hb_AP_margin = 0.05

    datasets = set(data_in.dataset)
    data = pd.DataFrame()

    for dataset in tqdm(datasets, desc='Parsing data sets'):
        # Take one data set
        dataset_data = data_in[data_in.dataset == dataset]
        gene = dataset_data.gene.iloc[0]

        if gene == 'hb':
            # Localize the maximum and keep the AP region described above
            max_AP = dataset_data.ap_mean.max()
            filtered_data = dataset_data[
                (dataset_data.ap_mean <= max_AP - hb_AP_margin) &
                (dataset_data.ap_mean >= hb_AP_margin)]
            data = pd.concat([data, filtered_data])

        elif gene == 'kn':
            frames = set(dataset_data.frame)
            # Individually process each frame
            for frame in frames:
                frame_data = dataset_data[dataset_data.frame == frame]
                median_AP = frame_data.ap.median()
                filtered_data = frame_data[
                    (frame_data.ap >= median_AP - knirps_half_width_AP) &
                    (frame_data.ap <= median_AP + knirps_half_width_AP)]
                data = pd.concat([data, filtered_data])

        elif gene == 'sn':
            filtered_data = dataset_data
            data = pd.concat([data, filtered_data])
    print('Data filtered!')

    return data


def identify_ncs(data):
    """
    Identify ncs in data traces by first thresholding and then numbering the ncs above the threshold.
    See README form more details.
    """
    # The sequential number of the last observed nc.
    # By default the last expected nc is nc14.
    # Add elements to the following dictionary if manual adjustements are necessary
    last_ncs = {
        1: 15,
        6: 14,
        54: 13,
        55: 13,
        56: 13,
        57: 13,
        58: 14,
        59: 14,
        60: 12,
        67: 13,
        68: 13,
        70: 13,
    }

    datasets_len = data.dataset_id.max() + 1
    ncs = range(11, 15)

    # Prepare output data frames
    data_out = copy.deepcopy(data)
    data_out['nc'] = np.nan
    iterables = pd.MultiIndex.from_product(
        [range(datasets_len), ncs], names=['dataset_id', 'nc'])
    nc_limits = pd.DataFrame(columns=['Tstart', 'Tend'], index=iterables)

    for dataset_id in trange(datasets_len, desc='Processing data sets'):
        # dataset_id = 54
        dataset_data = data[(data.dataset_id == dataset_id)]

        # Calculate average trace
        avg_trace, _ = get_avg_on_regular_time_mesh(dataset_data, dt=dt_new)

        # Threshold the average trace
        intensity_threshold = intensity_thresholds.get(dataset_id, default_intensity_threshold)
        avg_trace['is_expressing'] = avg_trace >= intensity_threshold

        # Detect the start of an nc as a sequence of 2 points, where the first one is below the threshold and the second one is above
        # Assign 1 to all such nc starts
        avg_trace['nc'] = 1 * (avg_trace.is_expressing &
                               np.logical_not(avg_trace.is_expressing.shift(1)))
        # Assign 1 also to the first time point
        avg_trace.loc[avg_trace.index[0], 'nc'] = 1 * \
            avg_trace.loc[avg_trace.index[0], 'is_expressing']

        # Number the cycles by performing a cumsum
        avg_trace['nc'] = avg_trace['nc'].cumsum() * avg_trace.is_expressing

        # Shift the numbers based on the last nc number
        if dataset_id in last_ncs:
            last_nc = last_ncs[dataset_id]
        else:
            last_nc = 14
        nc14_offset = avg_trace.nc.max()
        avg_trace['nc'] += (last_nc - nc14_offset)

        # If a frame is not expressing, set to nan
        avg_trace.loc[~avg_trace.is_expressing, 'nc'] = np.nan

        # %% Collect the nc time intervals in a separate table
        for nc in ncs:
            nc_data = avg_trace[avg_trace.nc == nc]
            if not nc_data.empty:
                Tstart = avg_trace[avg_trace.nc == nc].index[0]
                Tend = avg_trace[avg_trace.nc == nc].index[-1]
                nc_limits.loc[(dataset_id, nc), ['Tstart', 'Tend']] = [Tstart, Tend]

        # %% Use the table to label original non-averaged data set
        for nc in ncs:
            indices = ((data_out.time >= nc_limits.loc[(dataset_id, nc), 'Tstart'])
                       & (data_out.time <= nc_limits.loc[(dataset_id, nc), 'Tend'])
                       & (data_out.dataset_id == dataset_id)
                       )
            # Add as an extra column for the input data frame
            data_out.loc[indices, 'nc'] = nc
    return data_out, nc_limits


def perform_welchs_test(analyses_in):
    """
    Perform Welch's unequal variances t-test to evaluate whether alpha is the same across genes, constructs or ncs.
    Each of the three parts of the code corresponds to one changing parameter.

    All of the tests are performed for alpha since we saw that it is close to being normally distributed, while tau isn't.
    Similarities between other parameters are tested for historical reasons and are not re-used later on.
    """
    analyses = analyses_in.copy()

    genes = set(analyses.gene)
    constructs = set(analyses.construct)
    gene_ids = set(analyses.gene_id)
    ncs = set(analyses.index.get_level_values(1))
    quantities = ['r', 'j', 'alpha_over_k_comb', 'tau', 'alpha_comb']

    grouped = analyses.groupby(by=['gene', 'nc', 'construct'])
    idx = pd.IndexSlice

    # Initialize a temporary data frame
    calc = pd.DataFrame()
    for col in ['t', 'nu']:
        calc[col] = np.nan

    # %% Test for similarities across constructs
    # Some tests are only performed if all three constructs are availables
    has_all_constructs = 'bac' in constructs and 'no_pr' in constructs and 'no_sh' in constructs
    if not has_all_constructs:
        print(">> Skipping Welch's test across constructs: some of the constructs were not detected <<")
    else:
        print('>> p-values across constructs <<')
        for quantity in quantities:
            # Load/calculate means, variances and number of samples
            if quantity in ['r', 'j']:
                calc['means'] = grouped[quantity].mean()
                calc['vars'] = grouped[quantity].apply(
                    lambda group: np.var(group[~np.isnan(group)], ddof=1))

            elif quantity in ['alpha_over_k_comb', 'tau', 'alpha_comb']:
                calc['means'] = grouped[quantity].first()
                calc['vars'] = grouped[quantity + 'V'].first()

            calc['n'] = grouped[quantity].count()

            # Calculate independently for each gene and nc
            for gene in genes:
                for nc in ncs:
                    bac = calc.loc[(gene, nc, 'bac'), :]
                    no_sh = calc.loc[(gene, nc, 'no_sh'), :]

                    p, t, nu = welchs_test(
                        x1Mean=no_sh.means, x1V=no_sh.vars,
                        n1=no_sh.n, x2Mean=bac.means, x2V=bac.vars, n2=bac.n)

                    calc.loc[(gene, nc, 'bac'), ['t']] = t
                    calc.loc[(gene, nc, 'bac'), ['nu']] = nu
                    label = quantity + '_p_value'
                    calc.loc[(gene, nc, 'bac'), label] = p

                    # Copy the p-value back into the analyses table
                    gene_filter = analyses.gene == gene
                    analyses.loc[idx[gene_filter, nc], label] = p

                    # Print out results for `alpha_comb`
                    if not np.isnan(t) and quantity is 'alpha_comb':
                        print(f'{gene}, nc{nc}, bac, no_sh:\tp for alpha_comb: {p:.2g}')

    # %% Test across ncs for only alpha similarities.
    # Basically I compare nc14 to all other ncs for each gene-construct combination
    print('\n>> p-values across ncss <<')
    quantity = 'alpha_comb'
    nc1 = 14

    # Load means, vars and sample numbers
    calc['means'] = grouped[quantity].first()
    calc['vars'] = grouped[quantity + 'V'].first()
    calc['n'] = grouped[quantity].count()

    for gene in genes:
        for construct in constructs:
            for nc2 in ncs:
                if nc2 >= nc1:
                    continue
                data1 = calc.loc[(gene, nc1, construct), :]
                data2 = calc.loc[(gene, nc2, construct), :]

                p, t, nu = welchs_test(
                    x1Mean=data1.means, x1V=data1.vars,
                    n1=data1.n, x2Mean=data2.means, x2V=data2.vars, n2=data2.n)

                print(f'{gene}, {construct}, nc{nc2}, nc{nc1}:\tp for alpha: {p:.2g}')

    # %% Similarity test for alpha across genes
    if len(genes) < 2:
        print("\n>> Skipping Welch's test across genes: not enough genes <<")
    else:
        print('\n>> p-values across genes <<')

        grouped = analyses.groupby(by=['gene', 'nc', 'construct'])
        index = pd.MultiIndex.from_product([ncs, genes, genes])

        quantities = ['alpha_comb', 'tau']
        print('p-values across genes: ')
        for construct in constructs:
            for nc in ncs:
                for gene_id1 in gene_ids:
                    for gene_id2 in gene_ids:
                        if gene_id2 <= gene_id1:
                            continue
                        p = {}
                        for quantity in quantities:

                            # Create a temporary data frame
                            calc = pd.DataFrame()
                            calc['means'] = grouped[quantity].first()
                            calc['vars'] = grouped[quantity + 'V'].first()
                            calc['n'] = grouped[quantity].count()

                            gene1, gene2 = [gene_labels[i] for i in [gene_id1, gene_id2]]
                            data1 = calc.loc[(gene1, nc, construct), :]
                            data2 = calc.loc[(gene2, nc, construct), :]

                            p[quantity], t, nu = welchs_test(
                                x1Mean=data1.means, x1V=data1.vars,
                                n1=data1.n, x2Mean=data2.means, x2V=data2.vars, n2=data2.n)
                        print(
                            f'{construct}, nc{nc}, {gene1}, {gene2}:\tp for {quantities[0]} - {p[quantities[0]]:.2g},\tp for {quantities[1]} - {p[quantities[1]]:.2g}')
    return analyses


if __name__ == '__main__':
    """Some testing code"""
    alpha_over_k_J(0.00949, 3.11e-8)
    alpha_over_k_rho(0.1, 0.0049)
