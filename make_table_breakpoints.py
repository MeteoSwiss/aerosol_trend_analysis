import pandas as pd
import numpy as np

def make_table_breakpoints(Break_result, param):
    rows = []

    for i in range(len(Break_result)):
        station = Break_result.loc[i, "station"]

        # Saltar si està buit
        if station is None or (isinstance(station, (list, str)) and len(station) == 0):
            continue

        results = Break_result.loc[i, "results"]

        times = results["time"]

        for j in range(len(times)):
            row = {
                "param": Break_result.loc[i, "parameter"],
                "Break_point": results["time"][j],
                "p_value": results["pvalue"][j],
                "p_value_boot": results["pvalue_boot"][j],
                "level": results["level"][j],
                "minDiff": results["PrctDiff"][j][0],
                "medDiff": results["PrctDiff"][j][1],
                "maxDiff": results["PrctDiff"][j][2]
            }

            rows.append(row)

    # Crear DataFrame
    T = pd.DataFrame(rows)

    # Ordenar per data
    T = T.sort_values(by="Break_point").reset_index(drop=True)

    return T