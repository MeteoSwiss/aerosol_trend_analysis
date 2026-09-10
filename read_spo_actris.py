import os
import numpy as np
import pandas as pd

def read_spo_actris(filepath):
    """
    Python equivalent of MATLAB read_spo_actris
    Returns a pandas DataFrame indexed by datetime (UTC-like)
    """

    # 0. Constants
    N_HEADER_LINES = 1
    VARNAMES_LINE  = 3
    MISSING_FLAGS = [99999.99, 99999.999]
    MISSING_TOL = 1e-4

    # 1. Verify file
    if not os.path.isfile(filepath):
        raise FileNotFoundError(filepath)
    
    # 2. Read header lines
    with open(filepath, "r", encoding="utf-8") as f:
        header_lines = [next(f).strip() for _ in range(N_HEADER_LINES)]

    # 3. Extract column names 
    var_line = header_lines[0].strip()
    all_names = [x.strip() for x in var_line.split(",")]

    if all_names[0].lower() != "datetimeutc":
        print(f"Warning: expected DateTimeUTC, got {all_names[0]}")

    data_var_names = all_names[1:]

    # make names Python-safe
    # data_var_names = pd.io.parsers.ParserBase({'names': data_var_names})._maybe_dedup_names(data_var_names)

    # 4. Read full data
    df = pd.read_csv(filepath,skiprows=N_HEADER_LINES,header=None,names=["DateTimeUTC"] + data_var_names,
        na_values=["", " "],keep_default_na=True,encoding="utf-8")

    # 4. Convert datetime
    df["DateTimeUTC"] = pd.to_datetime(df["DateTimeUTC"],format="mixed",errors="coerce")

    df = df.set_index("DateTimeUTC")

    # 5. Force numeric conversion
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

        # replace sentinel values
        for flag in MISSING_FLAGS:
            df.loc[np.isclose(df[col], flag, atol=MISSING_TOL), col] = np.nan

    # 6. Report
    n_rows, n_vars = df.shape
    n_nan = df.isna().sum().sum()
    pct_nan = 100 * n_nan / (n_rows * n_vars)

    print(f"Read completed: {n_rows} rows x {n_vars} variables")
    print(f"NaN values: {n_nan} ({pct_nan:.1f}%)")
    print(f"Period: {df.index.min()} -> {df.index.max()}")

    return df