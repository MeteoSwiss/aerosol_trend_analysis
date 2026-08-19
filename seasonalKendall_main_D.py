import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mannkendall as mk
import os

def seasonalKendall_main_D(data, param, inst, station, resolution, granularity="daily",
    end_year=2025 ,path="", period=(10, 20, 30, 40, 50), fig=False,):
    
    """ compute the s.s of the selected data with PW and TFPW method
    compute the slope with VCTFPW method
     if not specified, end_year==2025 and the trends are computed for 10, 20 and 30 years
     if not specified, the granularity is one day with median as averaging
    
    IN: data:data either as a structure with time in data.start_time
                         or a timetable
                         or a matrice with the time as the first 6 columns (datevec)
     param = cell: either sturcture or timetable name or number of matrix column
     instrument= specify the instrument
     resolution= resolution is taken into account to determine the number of ties.
     varargin:
     granularity= {'daily','monthly','qualerly','yearly'} or length in [day]
                   default= 'daily'
     averageM= average method {'mean',@median}
           default=@median
     period= period for trend analysis, default: 10, 20, 30 and 40 years before end_time
     end_year, default: 2025
    fig= if 1, figure are plotted and saved
    
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
                   
                   period:year or 12 monthes + year or 4 meteorological seasons+year
                   ss= statistical significance: 95=95 for both PW and TFPW tests
                                         90=90 for both PW and TFPW tests
                                         10= test is s.s at 95 for TFPW and not s.s. for PW (false positive)
                                         20= test is s.s. at 95 for PW and not s.s. for TFPW (flase negative)
                                         0= not s.s. for both tests
                   slope: Sen's slope in units/y
                   UCL: upper confidence level in units/y
                   LCL: lower confidence level in units/y
    
     plot + save of the monthly and metSeason results
    
     dependency: call MK_tempAggr from https://github.com/mannkendall
    
    example: Tresult=seasonalKendall_main_D(test2, {'conc'},'cpc', 'JFJ', 1,'end_year',2018,'period',10,'fig',1);
    
    Martine Collaud Coen, July 2026
    
     ATTENTION, PATH FOR SAVING THE FIGURE IN LINES 2013-2017. PLEASE ADAPT FOR YOUR USAGE"""
    
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

    freq = {"daily": "D",
        "monthly": "MS",
        "quarterly": "QS",
        "yearly": "YS",}
    
    #aggregate data to the desired granularity
    if isinstance(granularity, (int, float)):
        dataG = dataT.resample(f"{granularity}D").median()
    else:
        dataG = dataT.resample(freq[granularity]).median()

    # Year column
    dataG["y"] = dataG.index.year

    # take the number of periods to be analysed
    n_years = dataG["y"].max() - dataG["y"].min() + 1
    # period can be either a single integer or a list/tuple
    if isinstance(period, (int, float)):
        period = [int(period)]
    elif isinstance(period, (list, tuple, np.ndarray)):
        period = list(period)
    else:
        raise TypeError("period must be an int, float, list, tuple or numpy array")
    
    if len(period) > 1:
        nb_period = int(n_years // 10)
        if period[-1] - period[-2] < 10:
            nb_period -= 1
        nb_period = max(nb_period, 1)
    else:
        nb_period = 1

    # use different distance between markers in plot as a function of the
    # number of present day trends
    if nb_period == 1:
        delta = 0.4
    elif nb_period == 2:
        delta = 0.3
    else:
        delta = 0.8 / (nb_period - 1)

    #for each period
    results = []

    for i in range(nb_period):
        start = pd.Timestamp(end_year - period[i] + 1, 1, 1)
        end = pd.Timestamp(end_year + 1, 1, 1)
        dataGTT = dataG.loc[start:end].copy()

        # 1) compute yearly trend
        # deseasonalse data for yearly computing is statistically needed
        # NOT USE FOR THIS ANALYSIS
        # # # # compute trend on all deseasonalised data
        # # # # with LMS fit of sin and cose
        # # # dataGTT_deseason=LMS_residue_season(dataGTT, param);

        #deseasonliased by dividing by monthly median 
        # and multiply by the total median to keep similar values

        dataGTT["m"] = dataGTT.index.month
        dataGTT_deseason = pd.DataFrame(np.nan,index=dataGTT.index,columns=[name],dtype=float)

        total_median = np.nanmedian(dataGTT[name])
        month_med = {}

        for m in range(1, 13):
            values = dataGTT.loc[dataGTT["m"] == m,name]
            month_med[m] = np.nanmedian(values)

            dataGTT_deseason.loc[dataGTT["m"] == m,name] = (values* total_median/ month_med[m])

        # Remove NaN and inf
        mask = (dataGTT_deseason[name].notna() & np.isfinite(dataGTT_deseason[name]))

        # compute the MK trend with the MK code from github
        obs_dates = dataGTT_deseason.index[mask].to_pydatetime()
        obs_values = dataGTT_deseason.loc[mask, name].to_numpy(dtype=float)
        result_y = format_mk_result(mk.mk_temp_aggr([obs_dates], [obs_values], resolution))

        # 2) compute trend with four meteorological seasons
        dataGTT["m"] = dataGTT.index.month

        multi_obs_dts = []
        multi_obs = []

        for m in range(1, 5):
            if m == 4:
                obs = dataGTT.loc[dataGTT["m"].isin([12, 1, 2]),name,]
            else:
                obs = dataGTT.loc[dataGTT["m"].between(m*3, m*3+2),name,]

            mask = obs.notna() & np.isfinite(obs)

            multi_obs_dts.append(obs.index[mask].to_pydatetime())
            multi_obs.append(obs.loc[mask].to_numpy(dtype=float))

        result_MS = format_mk_result(mk.mk_temp_aggr(multi_obs_dts, multi_obs, resolution))

        # 3) compute trend with 12 months
        multi_obs_dts = []
        multi_obs = []

        for m in range(1, 13):
            obs = dataGTT.loc[dataGTT["m"] == m,name,]

            mask = obs.notna() & np.isfinite(obs)

            multi_obs_dts.append(obs.index[mask].to_pydatetime())
            multi_obs.append(obs.loc[mask].to_numpy(dtype=float))

        result_mo = format_mk_result(mk.mk_temp_aggr(multi_obs_dts, multi_obs, resolution))

        # positions for plots
        t_mo = np.arange(1,14) + i*delta
        t_MS = np.arange(1,6) + i*delta

        if fig:

            fig_obj, axes = plt.subplots(2, 1, num=100, figsize=(10, 8))
            ax1, ax2 = axes
            title = (f"{station} {name} {inst}" f"{end_year-period[i]+1}-{end_year}")

            fig_obj.suptitle(title)
            plt.sca(ax1)
            fig_seasonKendall_trend(t_mo,result_mo,result_y,'units',delta)

            plt.sca(ax2)
            fig_seasonKendall_trend(t_MS,result_MS,result_y,'units',delta)

            ax2.text(t_MS[0], result_MS["UCL"].iloc[0] + abs(result_MS["UCL"].iloc[0])*0.3, f"{period[i]}y", fontsize=12, ha="center")

        for seasonality, result in [("y", result_y), ("MetSea", result_MS), ("month", result_mo)]:

            results.append({"station": station,
                "end_time": end_year,
                "length_period": period[i],
                "granularity": granularity,
                "parameter": name,
                "instrument": inst,
                "MK_seasonality": seasonality,
                "method": "MK",
                "ss": result["ss"].iloc[0],
                "slope": result["slope"].iloc[0],
                "UCL": result["UCL"].iloc[0],
                "LCL": result["LCL"].iloc[0],})
            
    Tresult = pd.DataFrame(results)

    if fig:
        if os.name == "nt":  # Windows
            save_path = (f"C:/github_trend/result/{station}/" f"{station}_{name}_MK.png")
        else:
            save_path = (f"C:/prod/pay/Aerosol_actris_trend/result/" f"{station}/{station}_{name}_MK.png")

        fig_obj.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig_obj)

    return Tresult


#_______________________________________________________________________
def format_mk_result(result):
    result = pd.DataFrame(result).T.reset_index(drop=True)
    result = result.rename(columns={"ucl": "UCL","lcl": "LCL"})

    return result

#_______________________________________________________________________
def fig_seasonKendall_trend(t, result, result_y, type, delta):
    if type == "units":
        s = "slope"
        U = "UCL"
        L = "LCL"
        plt.ylabel("Slope [units/y]")

    elif type == "%":
        s = "slopeP"
        U = "UCLP"
        L = "LCLP"
        plt.ylabel("Slope [%/y]", fontsize=14)

    sizeM = [15,20,30]

    plt.plot(t, result[s], "b.", markersize=sizeM[0])
    ss = result["ss"]
    plt.plot(np.array(t)[ss.values == 90], result.loc[ss == 90, s], "b.", markersize=sizeM[1])
    plt.plot(np.array(t)[ss.values == 95], result.loc[ss == 95, s], "b.", markersize=sizeM[2])
    plt.plot(np.array(t)[ss.values == 10], result.loc[ss == 10, s], "c^", markersize=sizeM[0])
    plt.plot(np.array(t)[ss.values == 20], result.loc[ss == 20, s], "cv", markersize=sizeM[0])

    for x, u, l in zip(t, result[U], result[L]):
        plt.plot([x,x], [l,u], linewidth=1,color='b')

    if result_y["ss"].iloc[0] == 95:
        a = sizeM[2]
    elif result_y["ss"].iloc[0] == 90:
        a = sizeM[1]
    else:
        a = sizeM[0]

    x_year = t[-1] + delta/2

    plt.plot(x_year, result_y[s].iloc[0], "g.", markersize=a)
    plt.plot([x_year, x_year], [result_y[L].iloc[0], result_y[U].iloc[0]],color='g')

    if len(t) == 13:
        plt.xlabel("Month", fontsize=14)
        plt.xticks(range(1,14),
                ["Jan","Feb","Mar","Apr",
                "May","Jun","Jul","Aug",
                "Sep","Oct","Nov","Dec",
                "year"])
        plt.xlim(0.8,14)

    elif len(t)==5:
        plt.xlabel("Season", fontsize=14)
        plt.xticks(range(1,6),["Spring", "Summer", "Fall", "Winter", "year"])
        plt.xlim(0.8,6)

    plt.grid(True)
    plt.axhline(0, color='r', linewidth=1)
    plt.gca().tick_params(labelsize=14) 