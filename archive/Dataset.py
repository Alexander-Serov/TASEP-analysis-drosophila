
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


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
