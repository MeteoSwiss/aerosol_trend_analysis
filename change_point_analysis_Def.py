import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from LMS_residue import LMS_residue
from multiple_breakpoints_while_one import multiple_breakpoints_while_one
from multiple_breakpoints_while_one_snht import multiple_breakpoints_while_one_snht
from multiple_breakpoints_2y_period import multiple_breakpoints_2y_period

"""
 apply the SNHT and Pettitt and monthly data
 the tests are applied only for at least 10 years on yearly data -->
 subsequent division of the time serie also follows this rule
 tests on monthly data should not be applied on less than 2 years if divided

 find a way to extract the best potentiel change points (and not all)
 save fig all in one

 inputs: data = TimeTable
           params = parameters (as fieldnames) to be analysed
           alpha = confidence level, set at 0.05 as default

 outputs:  result= table of results
 use: ...
 example:
 fev-mars 2026, mco"""

def change_point_analysis_Def(data, params, alpha, station, inst):
    
    data = data.copy()

    # remove lines with NaN
    # data = data.dropna(subset=params)

    # produce yearly and monthly medians
    data_m = data.resample("ME").median()

    granu = 'month'
    data_m[granu] = data_m.index.month

    start_time = data_m.index.year.min()
    end_time = data_m.index.year.max()
    period = end_time - start_time + 1

    # deseasonalise and detrend by LMS on data
    data_m_residue = pd.DataFrame(index=data_m.index)

    # prepare figure
    nb_subplot = 3

    cmap = plt.get_cmap("turbo", len(params))
    couleur = cmap(range(len(params)))

    Tresult = []

    
    # loop per paràmetres
    for i, p in enumerate(params):
        data_m_residue[p] = LMS_residue(data_m, p, "lin")

        # search Pettitt breakpoints for monthly data detrend and deseasolized with LMS
        result_mPT = multiple_breakpoints_while_one(
            data_m_residue, p, 12, alpha)

        Tresult.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": granu,
            "parameter": f"{p}_residue",
            "instrument": inst,
            "method": "Pettitt",
            "results": result_mPT
        })

        # search SNHT breakpoints for monthly data detrend and deseasolized with LMS
        result_mSNHT = multiple_breakpoints_while_one_snht(
            data_m_residue, p, 12, alpha)

        Tresult.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": granu,
            "parameter": f"{p}_residue",
            "instrument": inst,
            "method": "SNHT",
            "results": result_mSNHT
        })

    # plot of data with breakpoints and cusum
    Tresult_df = pd.DataFrame(Tresult)
    fig, axes = plt.subplots(nb_subplot, 1, figsize=(12, 10), sharex=True)

    if nb_subplot == 1:
        axes = [axes]
    ax1_right = axes[0].twinx()
    ax2_right = axes[1].twinx()

    for i, p in enumerate(params):
        # 1. MAIN PLOT (p-values)
        ax = axes[0]

        ax.set_title(station +' '+ inst +' ')

        # left axis: original data
        ax.plot(data_m.index,data_m[p],'-',color=couleur[i],label=p)
        ax.set_ylabel(p)

        ax2 = ax1_right

        # helper to extract results
        def get_select(name, method):
            mask = (
                (Tresult_df["parameter"] == f"{name}") &
                (Tresult_df["granularity"] == "month") &
                (Tresult_df["method"] == method)
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
        
        # Pettitt
        selectPT = get_select(f"{p}_residue", 'Pettitt')
        if selectPT is not None:
            clean_plot(selectPT, selectPT["pvalue"], 'o', "residue")

        # SNHT
        selectSNHT = get_select(f"{p}_residue", 'SNHT')
        if selectSNHT is not None:
            clean_plot(selectSNHT, selectSNHT["pvalue"], 's', "residue")
            
        ax.grid(True)
        
        # 2. SECOND PLOT (PrctDiff)
        ax = axes[1]

        ax.set_title('Residues from deseason + detrend')

        ax.plot(data_m_residue.index, data_m_residue[p], '-', color=couleur[i])
        ax.set_ylabel(p)

        ax2 = ax2_right

        # deseason
        # res = get_select(f"{p}_deseason")
        if selectPT is not None:
            clean_plot(selectPT, selectPT["PrctDiff"][:, 1], 'o', 'Pettitt')

        # res = get_select(f"{p}_residue")
        if selectSNHT is not None:
            clean_plot(selectSNHT, selectSNHT["PrctDiff"][:, 1], 's', 'SNHT')

        ax.grid(True)

        # 3. CUMSUM PLOT
        ax = axes[2]

        ax.set_title('Cumulative sum')

        ax.plot(data_m_residue.index,np.nancumsum(data_m_residue[p]),'--',color=couleur[i],label="residue")

        ax.set_ylabel(f"CumSum {p}")
        ax.grid(True)
        ax.legend(loc='upper left',bbox_to_anchor=(1, 1))
    
    #legends
    legend_elements = [
        Line2D([0], [0], marker='o', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='Pettitt'),
        Line2D([0], [0], marker='s', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='SNHT')]

    ax1_right.legend(handles=legend_elements,loc='upper left',bbox_to_anchor=(1.09, 1))
    ax1_right.set_ylabel("p-value")

    legend_elements2 = [Line2D([0], [0], marker='o', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='Pettitt'),
        Line2D([0], [0], marker='s', linestyle='None',markerfacecolor='none',
               markeredgecolor='k', label='SNHT')]

    ax2_right.legend(handles=legend_elements2,loc='upper left',bbox_to_anchor=(1.08, 1))
    ax2_right.set_ylabel("Median Diff")

    legend_elements3 = [Line2D([0], [0], linestyle='--', color='k', label='residue')]
    axes[2].legend(handles=legend_elements3,loc='upper left',bbox_to_anchor=(1.01, 1))
    
    plt.tight_layout()
    plt.show()
    return Tresult