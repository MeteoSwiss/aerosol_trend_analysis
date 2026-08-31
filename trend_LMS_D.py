import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
from scipy import stats
from statsmodels.regression.linear_model import yule_walker

def trend_LMS_D(data,param,inst,station,distribution,end_year=2025,
                path="",period=(10, 20, 30, 40, 50),fig=True,):
    
    """LMS fit is applied on the montly mean (median if the distribution in 'log') of the data
    The LMS consists of:
     a constant, a linear slope, seasonality
     if not specified, end_year==2025
     the trend are computed for the whole time series and complete decades

    IN
     data:data either as a structure with time in data.start_time
                         or a timetable
                         or a matrice with the time as the first 6 columns (datevec)
     param = cell: either sturcture or timetable name or number of matrix column
     instrument= specify the instrument
    distribution: 'lin' or 'log'
    
    varargin:
    averageM= average method {'mean',@median}
           default=@nanmedian
     period= period for trend analysis, default: 10, 20, 30
     end_year, default: 2025
    fig: if 1, plot and save figure
    OUT: Tresult= table of results with
                   station= name of the measuring site
                   end_time= end year used for the analyis
                   length_period= analysed period
                   granularity= used average of the data
                   parameter = name of the analysed parameter
                   instrument = instrument measuring the parameter
                   MK_seasonality = define if the Mann-Kendall is applied on
                                    the whole year, on 4 meteorological seasons or on 12 months
                   method = trend analysis method
                   ss= statistical significance: 95=95% if the significance (Weatherhead) >2
                                         90=90% if the significance(Weatherhead) >1.67
                   slope: slope in units/y with log if request
                   UCL: upper confidence level in units/y with log if request
                   LCL: lower confidence level in units/y with log if request
                   median: median of the data with log if request
                   slopeP: Sen's slope in %/y with log if request
                   UCLP: upper confidence level in %/y with log if request
                   LCLP: lower confidence level in %/y with log if request
                   slopeR: Sen's slope in % for the whole analysed period without log
                   UCLR: upper confidence level in % for the whole analysed period without log
    r              LCLR: lower confidence level in % for the whole analysed period without log
                   significance =significance= slope/variance
                   variance=variance from Weatherhead

     create a figure containing
     1. monthly average, fitted data and trend, 2. residuals,
     3. normplot of residuals, 4.cummulative summation of residues

     from Weatherheat, JGR 1998 et 2000
     Martine Collaud Coen, 7. 2026

     ATTENTION, PATH FOR SAVING THE FIGURE IN LINES 2005-2009. PLEASE ADAPT FOR YOUR USAGE"""

    fig_obj = None
    dataGa = None
    x_residue = None

    # Check input format
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")
    
    if not isinstance(data.index, pd.DatetimeIndex):
        raise TypeError("data must have a DatetimeIndex")
    
    if isinstance(param, str):
        name = param
    elif isinstance(param, (list, tuple)) and len(param) == 1:
        name = param[0]
    else:
        raise ValueError("param must be a string or a list containing one string")
    
    if name not in data.columns:
        raise KeyError(f"{name} not found in data")
    
    # Keep only the selected variable
    dataT = data[[name]].copy()
    dataT[name] = dataT[name].replace([np.inf, -np.inf], np.nan)

    # monthly average
    dataG = dataT.resample("MS").median()
    dataG[name] = np.real(dataG[name])
    # Year column
    dataG["y"] = dataG.index.year

    # if necessary, take the logarithm of the data
    if distribution == "log":
        dataG.loc[dataG[name] <= 0, name] = np.nan
        dataG[name] = np.log(dataG[name])

    mask = dataG[name].notna() & np.isfinite(dataG[name])
    dataG = dataG.loc[mask]

    # take the number of periods to be analysed
    n_years = dataG["y"].max() - dataG["y"].min() + 1
    if len(period) > 1:
        if period[-1] - period[-2] < 10:
            nb_period = int(n_years // 5) - 1
            if nb_period == 0:
                nb_period = 1
        else:
            nb_period = int(n_years // (period[-1] - period[-2]))
    else:
        nb_period = 1

    # compute the LMS trend for all periods
    result = []
    for i in range(nb_period-1, -1, -1):
        start = pd.Timestamp(end_year - period[i] + 1,1,1)
        end = pd.Timestamp(end_year,1,1)
        dataGa = dataG.loc[start:end].copy()
        # center the time vector
        t = (dataGa.index - dataGa.index[0]).days.to_numpy()

        # define the matrix X with all fit dependencies
        E = np.column_stack([np.ones(len(t)),t,
            np.sin((2*np.pi/365.25)*t),
            np.sin((4*np.pi/365.25)*t),
            np.sin((8*np.pi/365.25)*t),
            np.cos((2*np.pi/365.25)*t),
            np.cos((4*np.pi/365.25)*t),
            np.cos((8*np.pi/365.25)*t),])
        
        y = dataGa[name].to_numpy()

        b, _, _, _ = np.linalg.lstsq(E,y,rcond=None)
        # trend per year
        slope = b[1] * 365.25
        median = np.nanmedian(dataGa[name])
        slopeP = (b[1] * 365.25 * 100 /abs(median))

        # compute residuals
        Eplot = E @ b
        x_residue = dataGa[name].to_numpy() - Eplot

        # compute autocorrelation coef. phi and variance of white noise
        rho, sigma = yule_walker(x_residue,order=1)
        phi = rho[0]
        # compute the statistical significance:
        variance = (sigma**0.5/((1 - phi) * period[i]**1.5))

        # computation of the statistical significance
        significance = abs(b[1])*365.25 / variance
        if significance >= 2:
            ss = 95
        elif significance >= 1.67:
            ss = 90
        else:
            ss = 0

        # computation of the upper and lower confidence limits
        UCL = (b[1] + 2*variance/365.25)*365.25
        LCL = (b[1] - 2*variance/365.25)*365.25
        UCLP = UCL*100/abs(median)
        LCLP = LCL*100/abs(median)

        # compute the slope without log for all periods for log distribution
        # slope in % from the value at the beginning==1
        if distribution == "log":
            slopeR = np.exp(slope * period[i]) - 1
            UCLR = np.exp((slope + 2*variance/365.25)* period[i]) - 1
            LCLR = np.exp((slope - 2*variance/365.25)* period[i]) - 1
        elif distribution == "lin":
            slopeR = (slope * period[i]- 1)
            UCLR = ((slope + 2*variance/365.25)* period[i]- 1)
            LCLR = ((slope - 2*variance/365.25)* period[i]- 1)

        # create the table of results and do the figure only for lin distribution
        if i == nb_period-1:
            if fig:
                fig_obj = plt.figure(num=203, figsize=(9, 5.5))
                fig_LMS(dataGa,name,Eplot,x_residue,b)
                plt.subplot(2,2,1)
                plt.title(f"{station} {name} LMS {inst}")

            result.append({
                "station": station,
                "end_time": end_year,
                "length_period": period[i],
                "granularity": "month",
                "parameter": name,
                "instrument": inst,
                "distribution": distribution,
                "method": "LMS",
                "significance": significance,
                "ss": ss,
                "slope": slope,
                "UCL": UCL,
                "LCL": LCL,
                "slopeP": slopeP,
                "UCLP": UCLP,
                "LCLP": LCLP,
                "slopeR": slopeR,
                "UCLR": UCLR,
                "LCLR": LCLR,})
        else:
            if fig:
                plt.subplot(2, 2, 1)
                plt.title(f"{station} {name} LMS {inst}")

                plt.plot(dataGa.index,b[0] + b[1]*t,"-k",linewidth=2)

            result.append({
                "station": station,
                "end_time": end_year,
                "length_period": period[i],
                "granularity": "month",
                "parameter": name,
                "instrument": inst,
                "distribution": distribution,
                "method": "LMS",
                "significance": significance,
                "ss": ss,
                "slope": slope,
                "UCL": UCL,
                "LCL": LCL,
                "slopeP": slopeP,
                "UCLP": UCLP,
                "LCLP": LCLP,
                "slopeR": slopeR,
                "UCLR": UCLR,
                "LCLR": LCLR,})
    
    Tresult = pd.DataFrame(result)
    
    if fig and fig_obj is not None:
        if os.name == "nt":  # Windows
            save_path = (f"C:/github_trend/result/{station}/" f"{station}_{name}_LMS_D.png")
        else:
            save_path = (f"C:/prod/pay/Aerosol_actris_trend/result/" f"{station}/{station}_{name}_LMS_D.png")

        fig_obj.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig_obj)

    if dataGa is not None and x_residue is not None:
        data_residue = pd.DataFrame({"time": dataGa.index,"residue": x_residue})
    else:
        data_residue = pd.DataFrame({"time": pd.DatetimeIndex([]),"residue": []})

    return Tresult, fig_obj, data_residue

#_______________________________________________________________________
def fig_LMS(data, name, Eplot, x_residue, b):

    # Time vector in days from first point
    t = (data.index - data.index[0]).days.to_numpy()

    # 1) Data and fitted model
    plt.subplot(2, 2, 1)
    plt.grid(True)
    plt.plot(data.index,data[name],"bo-",linewidth=0.5,markerfacecolor="none")
    # complete LMS fit
    plt.plot(data.index,Eplot,"-r",linewidth=2)
    # linear trend only
    plt.plot(data.index,b[0] + b[1]*t,"-r",linewidth=2)
    plt.ylabel("monthly data and fit")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%y"))

    # 2) Residuals
    plt.subplot(2, 2, 2)
    plt.grid(True)
    plt.plot(data.index,x_residue)
    plt.ylabel("Residue")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%y"))

    # 3) Normal probability plot
    plt.subplot(2, 2, 3)
    normplot_matlab(x_residue)
    plt.title("normplot of residues")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%y"))

    # 4) Cumulative residuals
    plt.subplot(2, 2, 4)
    plt.grid(True)
    plt.plot(data.index,np.cumsum(x_residue))
    plt.ylabel("cumsum of residue")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%y"))
    plt.tight_layout()

#_______________________________________________________________________
def normplot_matlab(x_residue):

    (quantiles, values), (slope, intercept, r) = stats.probplot(x_residue, dist='norm')

    plt.plot(values, quantiles,'+')
    plt.plot(quantiles * slope + intercept, quantiles, '-.r')

    prob_ticks = np.array([0.01, 0.1, 1, 5, 10, 25,
        50, 75, 90, 95, 99, 99.9]) / 100

    plt.yticks(stats.norm.ppf(prob_ticks),
        ["0.001", "0.01", "0.1", "0.5", "0.10", "0.25",
         "0.50", "0.75", "0.90", "0.95", "0.99", "0.999"])

    plt.ylabel("Probability")
    plt.xlabel("Data")
    plt.grid(True)