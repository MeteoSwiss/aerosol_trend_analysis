import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from all3_trend import all3_trend
from seasonalKendall_main_D import seasonalKendall_main_D

def all_trend_STN(data_tr,data_st):
    """Do all necessary trends for one station
        including the present day trends with all potential decadal trends (10y, 20y, 30y, ..)
        and time series of all potential 10 y trends
        create a figure for all potential 10 y trends
        save results

        Dependency: call all3_trend & seasonalKendall_main_D
                all3_trend call seasonalKendall_main_D and trend_LMS_D"""
    data_tr = data_tr.copy()
    data_tr["y"] = data_tr.index.year
    exclude = ("y", "Variable", "Time", "Proper")

    names = [col for col in data_tr.columns
        if not col.startswith(exclude)]
    for i, name in enumerate(names):
        # give the resolution and instrument as a function of the name
        if name.startswith(("Bs", "Bbs", "expS")):
            inst = "neph"
            resolution = 0.005 if name.startswith("BbsF") else 0.01
        elif name.startswith(("Ba", "expA")):
            inst = "abs"
            resolution = 0.01
        elif name.startswith("SSA"):
            inst = "abs+neph"
            resolution = 0.005
        elif name.startswith("U"):
            inst = "RH"
            resolution = 0.1
        elif name.startswith("N"):
            inst = "cpc"
            resolution = 1
        elif name.startswith(("AOD", "AE")):
            inst = "pfr"
            resolution = 0.0002 if name.startswith("AOD") else 0.01
        
        # restrict the time series to period with data
        valid = data_tr[name].notna()
        if not valid.any():
            continue
        start = valid.idxmax()
        end = valid[::-1].idxmax()
        data_trok = data_tr.loc[start:end].copy()

        # compute the overall trend with end year = last year for all data!
        Tresult_MK_25, Tresult_LMSlogi, Tresult_LMSlini = all3_trend(data_trok,[name],inst,data_st["name"],resolution,
            end_year=data_trok["y"].max(),fig=True,)
        
        # compute all the 10y trend
        start_year = data_trok["y"].min()
        end_year = data_trok["y"].max()
        nb_trend = end_year - start_year - 8

        # check which years contain data
        years_with_data = (data_trok.groupby("y")[name].apply(lambda x: x.notna().any()))
        
        Tresult_MK_10y = []
        if nb_trend > 1:
            for j in range(2, nb_trend + 1): # 10y trend for last year is already computed
                trend_end_year = end_year - j + 1
                # trends are not computed if there is no data in the last year
                if not years_with_data.get(trend_end_year, False):
                    continue
                # or if there is less than 8 years with data
                window_years = range(trend_end_year - 9,trend_end_year + 1)
                nb_years_with_data = sum(years_with_data.get(y, False) for y in window_years)

                if nb_years_with_data > 7:
                    T = seasonalKendall_main_D(data_trok,[name],inst,data_st["name"],resolution,period=10,end_year=trend_end_year)
                    Tresult_MK_10y.append(T)

            # Join present-day trend + 10-year trends
            if len(Tresult_MK_25) > 0:
                Tresult_MK_10y = pd.concat(Tresult_MK_10y,ignore_index=True)
                Tresult_MKi = pd.concat([Tresult_MK_25,Tresult_MK_10y],ignore_index=True)
            else:
                Tresult_MKi = Tresult_MK_10y.copy()

        elif nb_trend==1:
            Tresult_MKi = Tresult_MK_25.copy()
        
        # Initialize result containers
        if i == 0:
            Tresult_MK = Tresult_MKi.copy()
            Tresult_LMSlog = Tresult_LMSlogi.copy()
            Tresult_LMSlin = Tresult_LMSlini.copy()
        else:
            Tresult_MK = pd.concat([Tresult_MK, Tresult_MKi],ignore_index=True)
            Tresult_LMSlog = pd.concat([Tresult_LMSlog, Tresult_LMSlogi],ignore_index=True)
            Tresult_LMSlin = pd.concat([Tresult_LMSlin, Tresult_LMSlini],ignore_index=True)

    return (Tresult_MK,Tresult_LMSlog,Tresult_LMSlin)


#_______________________________________________________________________
def fig_seasonKendall_trend10(data, data_st, namesP, type):
    sizeM = [15, 20, 30]
    if type == "units":
        s = "slope"
        U = "UCL"
        L = "LCL"
        plt.ylabel("Slope [units/y]")

    elif type == "%":
        s = "slopeP"
        U = "UCLP"
        L = "LCLP"
        plt.ylabel("Slope [%/y]")

    for i in range(len(data)):
        result = data["results"].iloc[i]
        slope = result[s][4]
        ss = result["ss"][4]

        if ss == 90:
            sizeR = sizeM[1]
            color = "r" if slope > 0 else "b"
        elif ss == 95:
            sizeR = sizeM[2]
            color = "r" if slope > 0 else "b"
        else:
            sizeR = sizeM[0]
            color = "k"

        year = data["end_time"].iloc[i]

        plt.plot(year,slope,marker=".",color=color,markersize=sizeR)

        plt.plot([year, year],[result[U][4],result[L][4]],color=color)
        if result["Xhomo"][4] == 1:
            plt.plot(year,slope,"rs",markersize=sizeR)

    plt.plot([data["end_time"].min()-0.2,data["end_time"].max()+0.2],[0,0],"r-")

    plt.ylabel(f"{data_st["name"]}_{namesP}")
    plt.grid(True)