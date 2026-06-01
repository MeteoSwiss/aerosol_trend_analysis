import numpy as np
import pandas as pd
from pettitt import pettitt

def multiple_breakpoints_while_one(data, param, nb_data_min, alpha):
    nbboot = 10000

    data = data.dropna(subset=[param]).copy()

    # Preallocació (aproximada)
    nb_preallocation = int(2 * np.ceil(len(data) / (nb_data_min / 2)))

    # Resultats
    results = []

    # Arbre
    tree = [None] * nb_preallocation
    tree[0] = data

    i = 0
    iteration = True

    while iteration:
        i += 1
        iteration = False

        start = 2**(i-1) - 1
        end = 2**i - 1

        for k in range(start, end):
            if k >= len(tree):
                break

            segment = tree[k]

            if segment is None or len(segment) == 0:
                continue

            if len(segment[param]) >= nb_data_min:

                x = segment[param].values

                # Pettitt
                a, PrctDiff = pettitt(x, alpha)

                # Bootstrap
                boot_stats = []
                for _ in range(nbboot):
                    sample = np.random.choice(x, size=len(x), replace=True)
                    a_boot, _ = pettitt(sample, alpha)
                    boot_stats.append(a_boot[1])

                boot_stats = np.array(boot_stats)
                pvalue_boot = (1 + np.sum(boot_stats >= a[1])) / (nbboot + 1)

                time_bp = None
                if not np.isnan(a[0]):
                    time_bp = segment.index[int(a[0])]

                    # Divisió
                    left = segment[segment.index < time_bp]
                    right = segment[segment.index > time_bp]

                    if 2*k+1 < len(tree):
                        tree[2*k+1] = left
                    if 2*k+2 < len(tree):
                        tree[2*k+2] = right

                    if (len(left) >= nb_data_min) or (len(right) >= nb_data_min):
                        iteration = True

                results.append({
                    "level": i,
                    "time": time_bp,
                    "pvalue": a[2],
                    "pvalue_boot": pvalue_boot,
                    "PrctDiff": PrctDiff
                })

            else:
                results.append({
                    "level": i,
                    "time": None,
                    "pvalue": np.nan,
                    "pvalue_boot": np.nan,
                    "PrctDiff": [np.nan, np.nan, np.nan]
                })

    # DataFrame final
    result_df = pd.DataFrame(results)

    # Neteja
    result_df = result_df.dropna(subset=["pvalue"], how="all")

    return result_df