import numpy as np
import pandas as pd
from snht import snht

def multiple_breakpoints_while_one_snht(data, param, nb_data_min, alpha):

    # --- checks equivalent to MATLAB arguments block
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")

    if not isinstance(data.index, pd.DatetimeIndex):
        raise TypeError("data must have a DatetimeIndex")

    if not isinstance(param, (str, list)):
        raise TypeError("param must be str or list")

    if not isinstance(nb_data_min, (int, float)):
        raise TypeError("nb_data_min must be numeric")

    if not isinstance(alpha, (int, float)):
        raise TypeError("alpha must be numeric")

    # while: process until each segment < nb_data_min
    nbboot = 10000
    # data = data.dropna(subset=[param]).copy()

    # Preallocació (aproximada)
    max_depth = int(np.ceil(np.log2(len(data[param]) / nb_data_min))) + 2
    nb_preallocation = 2**max_depth - 1

    # Resultats
    result = {
        "level": np.full(nb_preallocation, np.nan),
        "pvalue": np.full(nb_preallocation, np.nan),
        "time": [pd.NaT] * nb_preallocation,
        "pvalue_boot": np.full(nb_preallocation, np.nan),
        "PrctDiff": np.full((nb_preallocation, 3), np.nan)
    }

    # Arbre
    tree = [None] * nb_preallocation
    tree[0] = data

    i = 0
    iteration = True

    while iteration:
        i += 1
        iteration = False

        for k in range(2**(i-1)-1, 2**i-1):
            
            if tree[k] is not None and not tree[k].empty:

                n_valid = tree[k][param].notna().sum()
                if n_valid >= nb_data_min:

                    x = tree[k]

                    # SNHT 
                    try:
                        res = snht(x, param, alpha)
                    except ValueError:
                        continue

                    ttt = res.get("breakpoint", np.nan)
                    ppp = res.get("p_value", np.nan)
                    PrctDiff = res.get("percentile_diff", np.full(3, np.nan))

                    result["level"][k] = i
                    result["pvalue"][k] = ppp
                    result["PrctDiff"][k, :] = PrctDiff

                    if ttt is not None and not pd.isna(ttt):

                        result["time"][k] = pd.Timestamp(ttt)

                        left = 2 * k + 1
                        right = 2 * k + 2

                        if left < len(tree):
                            tree[left] = x[x.index < ttt]

                        if right < len(tree):
                            tree[right] = x[x.index > ttt]

                        left_ok = (left < len(tree)and tree[left] is not None
                            and len(tree[left][param]) >= nb_data_min)

                        right_ok = (right < len(tree)and tree[right] is not None
                            and len(tree[right][param]) >= nb_data_min)

                        if left_ok or right_ok:
                            iteration = True

                else:
                    result["level"][k] = i

            else:
                result["level"][k] = i
    
    pvalue = np.asarray(result["pvalue"])
    level = np.asarray(result["level"])
    prct = np.asarray(result["PrctDiff"])
    time = np.asarray(result["time"], dtype="object")
    if np.isnan(pvalue).sum() == len(pvalue):
        result["pvalue"] = pvalue
        result["level"] = level
        result["time"] = time
        result["PrctDiff"] = prct

    else:
        ind = ~np.isnan(pvalue)
        result["pvalue"] = pvalue[ind]
        result["level"] = level[ind]
        result["time"] = time[ind]
        result["PrctDiff"] = prct[ind, :]

    return result