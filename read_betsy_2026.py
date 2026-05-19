import numpy as np
import pandas as pd
import read_spo_actris

def read_betsy_2026(filename, STN):
    """
    Equivalent Python version of MATLAB read_betsy_2026
    STN= station name
    """

    # 1. Read data
    station_dat = read_spo_actris(filename)  # assumed DataFrame

    # 2. Convert to numeric (similar to timetable2num)
    station_dat = station_dat.apply(pd.to_numeric, errors='coerce')

    names_var = station_dat.columns.tolist()

    # 3. Station groups
    west = ['BBE','BLI','CRG','GWA','GLR','GBN','GRE','HGC','IBB',
            'MRN','MZW','PAZ','RMN','SCN','SIA','TMO','UBW']
    east = ['ACA','CRO','DSW','GSM','GGW','LBW','MCN','NCC','SHN']

    # helper flags (MATLAB contains bug fixed properly here)
    is_west = STN in west
    is_east = STN in east

    # 4. Scattering thresholds
    if is_west:
        bs_mask = [c for c in names_var if c.startswith(("Bs", "Bbs"))]

        for col in bs_mask:
            station_dat.loc[station_dat[col] > 100, col] = np.nan

        t_mask = [c for c in names_var if c.startswith(("T", "U"))]
        for col in t_mask:
            station_dat.loc[(station_dat[col] >= 9999) |
                (station_dat[col] >= 99999),col] = np.nan

    elif is_east:
        bs_mask = [c for c in names_var if c.startswith(("Bs", "Bbs"))]

        for col in bs_mask:
            station_dat.loc[station_dat[col] > 500, col] = np.nan

        t_mask = [c for c in names_var if c.startswith(("T", "U"))]
        for col in t_mask:
            station_dat.loc[station_dat[col] >= 9999, col] = np.nan

    # 5. Special cases
    if STN == "AMY":
        u_cols = [c for c in names_var if c.startswith("U")]
        for col in u_cols:
            station_dat.loc[station_dat[col] > 990, col] = np.nan

    if STN == "PAL":
        if "BbsG_S11" in station_dat:
            station_dat.loc[station_dat["BbsG_S11"] == 100, "BbsG_S11"] = np.nan
        if "BbsR_S11" in station_dat:
            station_dat.loc[station_dat["BbsR_S11"] == 100, "BbsR_S11"] = np.nan

    # 6. RH filtering and "dry" variables
    u_cols = [c for c in names_var if "U" in c]

    for rh_col in u_cols:

        rh_values = station_dat[rh_col]

        # threshold
        indRH = rh_values > 50

        if STN in ["MCN", "MRN"]:
            indRH = (rh_values > 50) | rh_values.isna()

        # associated scattering variables
        suffix = rh_col[1:]

        assoc_cols = [c for c in names_var if suffix in c and c.startswith("B")]

        for col in assoc_cols:
            dry_col = col + "_dry"
            station_dat[dry_col] = station_dat[col].copy()
            station_dat.loc[indRH, dry_col] = np.nan

    # 7. Daily aggregation
    station_datD = station_dat.resample("D").median()
    station_datN = station_dat.resample("D").count()

    # 8. Data coverage filtering
    bs_cols = [c for c in station_dat.columns if c.startswith("B")]

    a = []

    for col in bs_cols:
        if col in station_datN.columns:

            indN = station_datN[col] < 6

            valid = station_datD[col].notna().sum()
            total_bad = indN.sum()

            if valid > 0:
                a.append(total_bad * 100 / valid)
            else:
                a.append(np.nan)

            station_datD.loc[indN, col] = np.nan

    a = np.array(a)

    if np.any(a > 20):
        print("WARNING: there are many days with low data coverage")

    return station_datD