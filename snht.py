import numpy as np
import pandas as pd


def snht(data, var_name, alpha=0.05, n_boot=1000):

    # --------------------------------------------------
    # 1. input
    # --------------------------------------------------

    if isinstance(data, pd.DataFrame):

        if var_name is None:
            # agafa primera columna numèrica
            num_cols = data.select_dtypes(include=[np.number]).columns
            if len(num_cols) == 0:
                raise ValueError("No numeric columns found")
            var_name = num_cols[0]

        x = data[var_name].values
        time_index = data.index

    elif isinstance(data, pd.Series):
        x = data.values
        time_index = data.index

    else:
        x = np.asarray(data, dtype=float)
        time_index = None

    # remove NaN
    mask = ~np.isnan(x)
    x = x[mask]

    if time_index is not None:
        time_index = np.asarray(time_index)[mask]

    n = len(x)

    if n < 10:
        raise ValueError("Series must contain at least 10 samples")

    # --------------------------------------------------
    # 2. standardization
    # --------------------------------------------------

    mu = np.mean(x)
    sigma = np.std(x)

    if sigma == 0:
        raise ValueError("Constant series")

    z = (x - mu) / sigma

    # --------------------------------------------------
    # 3. compute T(k)
    # --------------------------------------------------

    cs = np.cumsum(z)

    k = np.arange(1, n)

    z1 = cs[:-1] / k
    z2 = (cs[-1] - cs[:-1]) / (n - k)

    T = k * z1**2 + (n - k) * z2**2

    # --------------------------------------------------
    # 4. max statistic
    # --------------------------------------------------

    T_max = np.max(T)
    cp = np.argmax(T) + 1

    # --------------------------------------------------
    # 5. bootstrap p-value
    # --------------------------------------------------

    z_sim = np.random.randn(n_boot, n)

    cs_sim = np.cumsum(z_sim, axis=1)

    T_boot = np.zeros(n_boot)

    for b in range(n_boot):

        z1b = cs_sim[b, :-1] / k
        z2b = (cs_sim[b, -1] - cs_sim[b, :-1]) / (n - k)

        T_sim = k * z1b**2 + (n - k) * z2b**2

        T_boot[b] = np.max(T_sim)

    p_value = np.mean(T_boot >= T_max)
    significant = p_value <= alpha

    # --------------------------------------------------
    # 6. breakpoint
    # --------------------------------------------------

    if significant:
        breakpoint = time_index[cp] if time_index is not None else cp
    else:
        breakpoint = np.nan

    # --------------------------------------------------
    # 7. percentiles
    # --------------------------------------------------

    before = x[:cp]
    after = x[cp:]

    p_before = np.percentile(before, [10, 50, 90])
    p_after = np.percentile(after, [10, 50, 90])
    p_all = np.percentile(x, [10, 50, 90])

    prct_diff = 100 * (p_after - p_before) / p_all

    return {
        "T_max": T_max,
        "breakpoint": breakpoint,
        "p_value": p_value,
        "significant": significant,
        "cp_index": cp,
        "T_series": T,
        "percentile_diff": prct_diff,
    }