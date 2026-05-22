import numpy as np
import pandas as pd
from pettitt import pettitt

def multiple_breakpoints_2y_period(data, param, nb_data_min, alpha):
    """ multiple_breakpoints iteratively search for breakpoints in a time series
     (data) with a defined method (e.g. pettitt). The minimal used period and/or the number of iterations can
     be given.
       Input: data
               method
             nb_data_min
             nb_iteration
    
     Outpur: results
    
     example:
     mco, March 2026
    # """    
    nbboot = 10000

    data = data.dropna(subset=[param]).copy()

    data["year"] = data.index.year

    nb_2yper = data["year"].iloc[-1] - data["year"].iloc[0]

    results = []

    for i in range(nb_2yper):
        start_year = data["year"].iloc[0] + i
        end_year = start_year + 2

        data_2yper = data[(data["year"] >= start_year) &
                          (data["year"] < end_year)]

        if len(data_2yper[param]) >= nb_data_min:

            x = data_2yper[param].values

            # Pettitt test
            a, PrctDiff = pettitt(x, alpha)

            # Bootstrap
            boot_stats = []
            for _ in range(nbboot):
                sample = np.random.choice(x, size=len(x), replace=True)
                a_boot, _ = pettitt(sample, alpha)
                boot_stats.append(a_boot[1])

            boot_stats = np.array(boot_stats)
            pvalue_boot = (1 + np.sum(boot_stats >= a[1])) / (nbboot + 1)

            # Temps del breakpoint
            time_bp = None
            if not np.isnan(a[0]):
                time_bp = data_2yper.index[int(a[0])]

            results.append({
                "time": time_bp,
                "pvalue": a[2],
                "pvalue_boot": pvalue_boot,
                "PrctDiff": PrctDiff,
                "level": i
            })

    # Convertir a DataFrame
    result_df = pd.DataFrame(results)

    # Eliminar NaNs
    result_df = result_df.dropna(subset=["pvalue"])

    return result_df