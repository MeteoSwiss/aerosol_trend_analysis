import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from LMS_residue import LMS_residue
from multiple_breakpoints_while_one import multiple_breakpoints_while_one
from multiple_breakpoints_while_one_snht import multiple_breakpoints_while_one_snht
from multiple_breakpoints_2y_period import multiple_breakpoints_2y_period

"""
 apply the SNHT or Pettitt and monthly data
 the tests are applied only for at least 10 years on yearly data -->
 subsequent division of the time serie also follows this rule
 tests on monthly data should not be applied on less than 2 years if divided
 find a way to extract the best potentiel change points (and not all)
 save fig all in one
 inputs: data = TimeTable
           params = parameters (as fieldnames) to be analysed
           alpha = confidence level, set at 0.05 as default
           model = snht or Pettitt
 outputs:  result= table of Tresult
 use: ...
 example:
 fev-mars 2026, mco"""

def change_point_analysis_v4(data, params, alpha, station, inst, model):

    data = data.copy()

    # remove lines with NaN
    # data = data.dropna(subset=params)

    # decide model
    if model == 'snht':
        func = multiple_breakpoints_while_one_snht
    elif model == 'Pettitt':
        func = multiple_breakpoints_while_one
    else:
        raise TypeError("model must be snht or Pettitt")


    # produce yearly and monthly medians
    data_m = data.resample("ME").median()

    data_m["month"] = data_m.index.month

    start_time = data_m.index.year.min()
    end_time = data_m.index.year.max()
    period = end_time - start_time + 1

    # deseasonlise with the median of the month
    month_med = data_m.groupby("month").median()

    data_m_deseason = pd.DataFrame(index=data_m.index)

    for p in params:
        data_m_deseason[p] = (data_m[p] / data_m["month"].map(month_med[p]))

    # deseasonalise and detrend by LMS on data
    data_m_residue = pd.DataFrame(index=data_m.index)
    # deseasonalise and detrend by LMS on log(data)
    data_m_residueLog = pd.DataFrame(index=data_m.index)

    if len(params) == 1:
        nb_subplot = 3
    else:
        nb_subplot = 4

    cmap = plt.get_cmap("turbo", len(params))
    couleur = cmap(range(len(params)))

    Tresult = []

    # loop per paràmetres
    for i, p in enumerate(params):
        data_m_residue[p] = LMS_residue(data_m, p, "lin")
        data_m_residueLog[p] = LMS_residue(data_m, p, "log")

        granu = "month"
        # --- 1. deseason
        result_m = func(data_m_deseason, p, 12, alpha)

        Tresult.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": granu,
            "parameter": f"{p}_deseason",
            "instrument": inst,
            "method": model,
            "results": result_m
        })

        # --- 2. residue
        result_m = func(data_m_residue, p, 12, alpha)

        Tresult.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": granu,
            "parameter": f"{p}_residue",
            "instrument": inst,
            "method": model,
            "results": result_m
        })

        # --- 3. residue log
        result_m = func(data_m_residueLog, p, 12, alpha)

        Tresult.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": granu,
            "parameter": f"{p}_residueLog",
            "instrument": inst,
            "method": model,
            "results": result_m
        })

        # --- 4. backward (invertit)
        data_inv = data_m_residueLog[p].iloc[::-1].copy()

        df_inv = pd.DataFrame({p: data_inv})

        result_m = func(df_inv, p, 12, alpha)

        # corregir temps
        corrected_times = []
        for t in result_m["time"]:
            if pd.notna(t):
                pos = data_m.index.get_indexer([t])[0]

                if pos != -1:
                    N = len(data_m.index)
                    corrected_times.append(data_m.index[N - pos - 1])
                else:
                    corrected_times.append(pd.NaT)
            else:
                corrected_times.append(pd.NaT)

        result_m["time"] = corrected_times

        Tresult.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": granu,
            "parameter": f"{p}_residueLog_inv",
            "instrument": inst,
            "method": model,
            "results": result_m
        })

        # # --- 5. finestres 2 anys
        # result_m = multiple_breakpoints_2y_period(
        #     data_m_residueLog, p, 12, alpha)

        # Tresult.append({
        #     "station": station,
        #     "start_time": start_time,
        #     "end_time": end_time,
        #     "length_period": period,
        #     "granularity": granu,
        #     "parameter": f"{p}_residueLog_2yper",
        #     "instrument": inst,
        #     "method": model,
        #     "results": result_m
        # })

    # 🔹 convertir a DataFrame final
    Tresult_df = pd.DataFrame(Tresult)
    fig, axes = plt.subplots(nb_subplot, 1, figsize=(12, 10), sharex=True)

    if nb_subplot == 1:
        axes = [axes]
    ax1_right = axes[0].twinx()
    ax2_right = axes[1].twinx()

    if nb_subplot >= 4:
        ax4_right = axes[3].twinx()

    for i, p in enumerate(params):
        # 1. MAIN PLOT (p-values)
        ax = axes[0]

        ax.set_title(station)

        # left axis: original data
        ax.plot(data_m.index,data_m[p],'.',color=couleur[i],label=p)
        ax.set_ylabel(p)

        ax2 = ax1_right

        # helper to extract results
        def get_select(name):
            mask = (
                (Tresult_df["parameter"] == f"{name}") &
                (Tresult_df["granularity"] == "month") &
                (Tresult_df["method"] == model)
            )
            if mask.sum() == 0:
                return None
            return Tresult_df.loc[mask, "results"].iloc[0]
        
        def clean_plot(res, val, marker, label):
            time = pd.to_datetime(res["time"], errors="coerce")
            pval = np.asarray(val, dtype=float)

            mask = time.notna() & np.isfinite(pval)

            ax2.plot(time[mask].to_numpy() + pd.Timedelta(days=15),pval[mask],marker,
                    markerfacecolor='none',markeredgecolor=couleur[i], markersize=10,label=label)

        # deseason
        res = get_select(f"{p}_deseason")
        if res is not None:
            clean_plot(res, res["pvalue"], 'o', "deseason")

        # residue
        res = get_select(f"{p}_residue")
        if res is not None:
            clean_plot(res, res["pvalue"], 's', "residue")

        # log residue forward
        res = get_select(f"{p}_residueLog")
        if res is not None:
            clean_plot(res, res["pvalue"], '^', "resLog forward")

        # backward
        res = get_select(f"{p}_residueLog_inv")
        if res is not None:
            clean_plot(res, res["pvalue"], 'v', "resLog backward")

        ax.grid(True)

        # 2. SECOND PLOT (PrctDiff)
        ax = axes[1]

        ax.plot(data_m_deseason.index, data_m_deseason[p], '.', color=couleur[i])
        ax.plot(data_m_residue.index, data_m_residue[p], '.', color=couleur[i])
        ax.set_ylabel(p)

        ax2 = ax2_right

        # deseason
        res = get_select(f"{p}_deseason")
        if res is not None:
            clean_plot(res, res["PrctDiff"][:, 1], 'o', 'deseason')

        res = get_select(f"{p}_residue")
        if res is not None:
            clean_plot(res, res["PrctDiff"][:, 1], 's', 'residue')

        ax.grid(True)

        # 3. CUMSUM PLOT
        ax = axes[2]

        ax.plot(data_m_deseason.index,np.nancumsum(data_m_deseason[p]),'-',color=couleur[i],label="deseason")
        ax.plot(data_m_residue.index,np.nancumsum(data_m_residue[p]),'--',color=couleur[i],label="residue")
        ax.plot(data_m_residueLog.index,np.nancumsum(data_m_residueLog[p]),':',color=couleur[i],label="residueLog")

        ax.set_ylabel(f"CumSum {p}")
        ax.grid(True)
        ax.legend(loc='upper left',bbox_to_anchor=(1, 1))

    # only if multiple parameters
    if len(params) > 1:
        for i in range(len(params) - 1):
            for j in range(i + 1, len(params)):

                p1 = params[i]
                p2 = params[j]
                ratio_name = f"{p1}_{p2}"

                # create ratio series
                data_m[ratio_name] = data_m[p1] / data_m[p2]

                # SNHT / breakpoint analysis
                result_m = func(data_m,ratio_name,12,alpha)

                # store result (like MATLAB table row)
                Tresult.append({
                    "station": station,
                    "start_time": start_time,
                    "end_time": end_time,
                    "length_period": period,
                    "granularity": "month",
                    "parameter": ratio_name,
                    "instrument": inst,
                    "method": model,
                    "results": result_m
                })

                # plotting (subplot 4)
                if nb_subplot >= 4:

                    ax = axes[3]  # subplot 4

                    ax2 = ax4_right

                    # ratio line
                    ax.plot(data_m.index,data_m[ratio_name],'-',color=couleur[i],label="ratio")

                    ax.set_ylabel(p)

                    # p-values
                    clean_plot(result_m, result_m["pvalue"], 'v', 'all ratios')

                    ax2.set_ylim(0, alpha)
                    ax2.set_ylabel("p-value")

    #legends
    legend_elements = [Line2D([0], [0], marker='o', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='deseason'),
        Line2D([0], [0], marker='s', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='residue'),
        Line2D([0], [0], marker='^', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='resLog forward'),
        Line2D([0], [0], marker='v', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='resLog backward')]

    ax1_right.legend(handles=legend_elements,loc='upper left',bbox_to_anchor=(1.09, 1))
    ax1_right.set_ylabel("p-value")

    legend_elements2 = [Line2D([0], [0], marker='o', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='deseason'),
        Line2D([0], [0], marker='s', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='residue')]

    ax2_right.legend(handles=legend_elements2,loc='upper left',bbox_to_anchor=(1.08, 1))
    ax2_right.set_ylabel("Median Diff")

    legend_elements3 = [Line2D([0], [0], linestyle='-', color='k', label='deseason'),
                        Line2D([0], [0], linestyle='--', color='k', label='residue'),
                        Line2D([0], [0], linestyle=':', color='k', label='residueLog'),]
    axes[2].legend(handles=legend_elements3,loc='upper left',bbox_to_anchor=(1.01, 1))
    
    if nb_subplot >= 4:
        legend_elements4 = [Line2D([0], [0], marker='v', linestyle='None',markerfacecolor='none',
                markeredgecolor='k', label='all ratios')]

        ax4_right.legend(handles=legend_elements4,loc='upper left',bbox_to_anchor=(1.08, 1))
        ax4_right.set_ylabel("p-value")

    plt.tight_layout()
    plt.show()
    return Tresult