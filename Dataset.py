
import numpy as np
import pandas as pd


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

    if t_max > 0:
        time_mesh_positive = np.arange(dt, t_max + dt, dt)

    if t_min < 0:
        time_mesh_negative = -np.arange(dt, -t_min + dt, dt)
        time_mesh_negative = time_mesh_negative[::-1]

    time_mesh = np.concatenate((time_mesh_negative, [0], time_mesh_positive))

    return time_mesh


def get_avg_on_regular_time_mesh(data, dt):
    """
    Calculate an average time trace on a regular time mes with a given dt
    """
    time_mesh = get_regular_time_mesh(data, dt)
    avg_data = pd.DataFrame(index=time_mesh[:-1], columns=['intensity'])

    for i in range(len(time_mesh) - 1):
        avg_data.loc[time_mesh[i]] = data[(data.time >= time_mesh[i])
                                          & (data.time < time_mesh[i + 1])].intensity.mean()

    # print(avg_data)
    return avg_data
