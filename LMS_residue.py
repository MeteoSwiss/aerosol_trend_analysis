import numpy as np
import pandas as pd

def LMS_residue(data, param, distribution):
    """
    LMS fit is applied on the montly mean (median if the distribution in 'log') of the data
    The LMS consists of:
    a constant, a linear slope, seasonalit�
    if not specified, end_year==2017 and the trend are computed for 10, 15, 20 and 25 years
    
    IN
    data: timetable
    param: parameter to  be analysied
    distribution: 'lin' o 'log'


    OUT: data_residue= residue of the fit
    programmed selon Weatherheat, JGR 1998 et 2000
    Martine Collaud Coen, 3.2026
    """

    data = data.copy()

    # if necessary, take the logarithm of the data
    if distribution == 'log':
        data[param] = np.log(data[param])

    data[param].replace([-np.inf, np.inf], np.nan, inplace=True)
    data[param] = np.real(data[param])

    datax = data.dropna(subset=[param])

    # Vector temps (dies)
    t = (datax.index - datax.index[0]).days.values

    # Matriu E
    E = np.column_stack([np.ones(len(t)), t,
        np.sin((2*np.pi/365.25)*t), np.sin((4*np.pi/365.25)*t), np.sin((8*np.pi/365.25)*t),
        np.cos((2*np.pi/365.25)*t), np.cos((4*np.pi/365.25)*t), np.cos((8*np.pi/365.25)*t)])

    # Ajust LMS
    y = datax[param].values
    b, *_ = np.linalg.lstsq(E, y, rcond=None)

    # Ajust model
    Eplot = E @ b

    # Residu
    x_residue = y - Eplot

    # Reconstrucció amb NaNs originals
    data_residue = pd.Series(np.nan, index=data.index)
    data_residue.loc[datax.index] = x_residue

    return data_residue