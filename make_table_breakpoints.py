import pandas as pd
import numpy as np

def make_table_breakpoints(Break_result,model):
    rows = []

    for item in Break_result:
        station = item["station"]

        # Saltar si està buit
        if station is None or (isinstance(station, (list, str)) and len(station) == 0):
            continue

        results = item["results"]

        for j in range(len(results["time"])):
            row = {"param": item["parameter"],
                "Break_point": results["time"][j],
                "p_value": results["pvalue"][j],}

            if model == "Pettitt":
                row["p_value_boot"] = results["pvalue_boot"][j]

            row.update({"level": results["level"][j],
                "minDiff": results["PrctDiff"][j][0],
                "medDiff": results["PrctDiff"][j][1],
                "maxDiff": results["PrctDiff"][j][2]})

            rows.append(row)

    # Crear DataFrame
    T = pd.DataFrame(rows)

    # Ordenar per data
    T = T.sort_values(by="Break_point").reset_index(drop=True)

    return T