import numpy as np
import pandas as pd
from snht import snht

def multiple_breakpoints_while_one_snht(data, param, nb_data_min, alpha):

    data = data.dropna(subset=[param]).copy()

    # Preallocació (aproximada)
    nb_preallocation = 2 * int(np.ceil(len(data[param]) / (nb_data_min / 2)))

    # Resultats
    result = {
        "level": np.full(nb_preallocation, np.nan),
        "pvalue": np.full(nb_preallocation, np.nan),
        "time": np.array([np.datetime64("NaT")] * nb_preallocation),
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

        for k in range(2**(i-1), 2**i):

            if k in tree and len(tree[k]) > 0:

                if len(tree[k][param]) >= nb_data_min:

                    x = tree[k]

                    # SNHT 
                    ttt, ppp, PrctDiff = snht(x, param, alpha)

                    result["level"][k] = i
                    result["pvalue"][k] = ppp
                    result["PrctDiff"][k, :] = PrctDiff

                    if ttt is not None and not pd.isna(ttt):

                        result["time"][k] = ttt

                        # split segment (arbre binari)
                        tree[2 * k] = x[x.index < ttt]
                        tree[2 * k + 1] = x[x.index > ttt]

                        if (len(tree[2 * k][param]) >= nb_data_min or
                            len(tree[2 * k + 1][param]) >= nb_data_min):
                            iteration = True

                else:
                    result["level"][k] = i

            else:
                result["level"][k] = i

    pvalue = result["pvalue"]

    if np.isnan(pvalue).sum() == len(pvalue):

        result["pvalue"] = pvalue[1:]
        result["level"] = result["level"][1:]
        result["time"] = result["time"][1:]
        result["PrctDiff"] = result["PrctDiff"][1:, :]

    else:

        ind = ~np.isnan(pvalue)

        result["pvalue"] = pvalue[ind]
        result["level"] = result["level"][ind]
        result["time"] = result["time"][ind]
        result["PrctDiff"] = result["PrctDiff"][ind, :]

    return result