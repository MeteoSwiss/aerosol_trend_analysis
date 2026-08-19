import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_10y_in_two(STN_result, STN_st, seasonality):
    """plot the 10y trends of all the measured parameters in one fig (RH, scat, backscat,abs)
    and the 10y trends of all the calculated parameters in another one (Backscat, expS, expA, SSA)
    
    
    ATTENTION, PATH FOR SAVING THE FIGURE IN LINES 83-87. PLEASE ADAPT FOR YOUR USAGE
    """
     # Station name
    if isinstance(STN_st, dict):
        station_name = STN_st["name"]
    else:
        station_name = STN_st.name

    # MEASURED PARAMETERS
    mask_measured = ((STN_result["length_period"] == 10) 
        & (STN_result["MK_seasonality"] == seasonality)
        & (STN_result["parameter"].astype(str).str.startswith(("U", "Bs", "Bbs", "Ba")))
        & ~STN_result["parameter"].astype(str).str.startswith("BbsF"))
    select1 = STN_result.loc[mask_measured].copy()
    names = sorted(select1["parameter"].unique())

    Tmin = select1["end_time"].min() - 0.5
    Tmax = select1["end_time"].max() + 0.5

    fig_measured, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    measured_groups = [("U", "RH"),
        ("Bs", "Bs"),
        ("Bbs", "Bbs"),
        ("Ba", "Ba"),]

    for ax, (prefix, ylabel) in zip(axes, measured_groups):

        names_group = [name for name in names if name.startswith(prefix)]

        if names_group:
            plt.sca(ax)
            mask = [name in names_group for name in names]
            plot_each_variable(select1, names, mask, Tmin, Tmax, "units", seasonality)
            ax.set_ylabel(ylabel, fontsize=16)

    axes[0].set_title(station_name, fontsize=16)
    axes[-1].set_xlim(Tmin, Tmax)

    fig_measured.tight_layout()

    # CALCULATED PARAMETERS
    mask_calculated = ((STN_result["length_period"] == 10)
        & (STN_result["MK_seasonality"] == seasonality)
        & (STN_result["parameter"].astype(str).str.startswith(("BbsF", "exp", "SSA"))))
    select1 = STN_result.loc[mask_calculated].copy()
    names = sorted(select1["parameter"].unique())

    fig_calculated, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    calculated_groups = [("BbsF", "Backscat. fraction"),
        ("expS", "Scat exponent"),
        ("expA", "Abs. exponent"),
        ("SSA", "SSA"),]

    for ax, (prefix, ylabel) in zip(axes, calculated_groups):

        names_group = [name for name in names if name.startswith(prefix)]

        if names_group:
            plt.sca(ax)
            mask = [name in names_group for name in names]
            plot_each_variable(select1, names, mask, Tmin, Tmax, "units", seasonality)
            ax.set_ylabel(ylabel, fontsize=16)

    axes[0].set_title(station_name, fontsize=16)
    axes[-1].set_xlim(Tmin, Tmax)

    fig_calculated.tight_layout()

    if os.name == "nt":  # Windows
        save_path_measured = (f"C:/github_trend/result/{station_name}/" f"{station_name}_MK_all10y_var.png")
        save_path_calculated = (f"C:/github_trend/result/{station_name}/" f"{station_name}_MK_all10y_cal.png")
    else:
        save_path_measured = (f"/prod/pay/Aerosol_actris_trend/result/{station_name}/" f"{station_name}_MK_all10y_var.png")
        save_path_calculated = (f"/prod/pay/Aerosol_actris_trend/result/{station_name}/" f"{station_name}_MK_all10y_cal.png")

    fig_measured.savefig(save_path_measured, dpi=300, bbox_inches="tight")
    fig_calculated.savefig(save_path_calculated, dpi=300, bbox_inches="tight")

    plt.close(fig_measured)
    plt.close(fig_calculated)

    return fig_measured, fig_calculated


#_______________________________________________________________________
def plot_each_variable(select1, names, mask, Tmin, Tmax, type, seasonality):
    N = sorted(np.asarray(names)[mask])

    for i, name in enumerate(N):
        select2 = select1.loc[select1["parameter"] == name].copy()

        fig_seasonKendall_trend10_m(select2, [name], type, i + 1, seasonality)

    plt.xlim(Tmin, Tmax)
    plt.plot([Tmin, Tmax], [0, 0], "k-", linewidth=2)

    handles, labels = plt.gca().get_legend_handles_labels()
    if handles:
        plt.gca().legend(handles, labels, fontsize=12)


#_______________________________________________________________________
def fig_seasonKendall_trend10_m(data, namesP, type, nb, seasonality):
    # Select result column according to seasonality
    if seasonality == "MetSea":
        colonne_res = 4
    elif seasonality == "y":
        colonne_res = 0
    elif seasonality == "month":
        colonne_res = 12
    else:
        raise ValueError("No plot because no seasonality was given")
    
    # Marker and size
    sizeM = [8, 8, 8]
    markerM = ["o", "s", "v", "d",
                "o", "s", "v", "d",
                "o", "s", "v", "d"]
    marker = markerM[nb - 1]

    # Select slope / confidence limits
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
    else:
        raise ValueError("type must be either 'units' or '%'")

    # Plot each trend
    for i in range(len(data)):
        ss = data["ss"].iloc[i]
        slope = data[s].iloc[i]

        # Determine size and colour
        if ss == 90:
            sizeR = sizeM[1]
            if slope > 0:
                colorT = "r"
            elif slope < 0:
                colorT = "b"
        elif ss == 95:
            sizeR = sizeM[2]
            if slope > 0:
                colorT = "r"
            elif slope < 0:
                colorT = "b"
            elif slope == 0:
                colorT = [0.7, 0.7, 0.7]
            elif np.isnan(slope):
                colorT = "k"
        else:
            sizeR = sizeM[0]
            colorT = "k"

        # Plot
        plt.plot(data["end_time"].iloc[i], slope, marker=marker, color=colorT,
            markersize=sizeR, markerfacecolor=colorT, label=namesP[0] if i == 0 else None)

    plt.grid(True)