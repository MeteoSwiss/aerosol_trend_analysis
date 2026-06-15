import numpy as np
import pandas as pd
from pettitt import pettitt

def multiple_breakpoints_while_one(data, param, nb_data_min, alpha):
    
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
    
    nbboot = 10000

    data = data.dropna(subset=[param]).copy()

    # Preallocació
    n = max(len(data[param]), 1)
    max_depth = max(int(np.ceil(np.log2(n / nb_data_min))) + 2,1)
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

                    x = tree[k][param]

                    # Pettitt
                    try:
                        a, PrctDiff = pettitt(x, alpha)
                    except ValueError:
                        continue

                    result["level"][k] = i
                    result["pvalue"][k] = a[2]
                    result["PrctDiff"][k, :] = PrctDiff
                    # Bootstrap
                    boot_stats = []
                    for _ in range(nbboot):
                        sample = np.random.choice(x, size=len(x), replace=True)
                        a_boot, _ = pettitt(sample, alpha)
                        boot_stats.append(a_boot[1])
    
                    boot_stats = np.array(boot_stats)
                    pvalue_boot = (1 + np.sum(boot_stats >= a[1])) / (nbboot + 1)
                    result["pvalue_boot"][k] = pvalue_boot
    
                    time_bp = None
                    if not np.isnan(a[0]):
                        time_bp = x.index[int(a[0])]
    
                        # Divisió
                        left = tree[k][tree[k].index < time_bp]
                        right = tree[k][tree[k].index > time_bp]
    
                        if 2*k+1 < len(tree):
                            tree[2*k+1] = left
                        if 2*k+2 < len(tree):
                            tree[2*k+2] = right
    
                        if (len(left) >= nb_data_min) or (len(right) >= nb_data_min):
                            iteration = True
                    result["time"][k] = pd.Timestamp(time_bp)

                else:
                    result["level"][k] = i
            else:
                result["level"][k] = i

    pvalue = np.asarray(result["pvalue"])
    pvalue_boot = np.asarray(result["pvalue_boot"])
    level = np.asarray(result["level"])
    prct = np.asarray(result["PrctDiff"])
    time = np.asarray(result["time"], dtype="object")

    valid = ~np.isnan(pvalue)

    if np.sum(valid) == 0:
        for key in result:
            result[key] = result[key][:1]

    else:
        result["level"] = level[valid]
        result["pvalue"] = pvalue[valid]
        result["pvalue_boot"] = pvalue_boot[valid]
        result["time"] = time[valid]
        result["PrctDiff"] = prct[valid]

    valid = ~pd.isna(time)

    result = {
        "level": level[valid],
        "pvalue": pvalue[valid],
        "pvalue_boot": pvalue_boot[valid],
        "time": time[valid],
        "PrctDiff": prct[valid, :]
    }
    return result