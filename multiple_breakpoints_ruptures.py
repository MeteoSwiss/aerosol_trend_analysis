import numpy as np
import pandas as pd
import ruptures as rpt

def multiple_breakpoints_ruptures(data, param, nb_data_min, pen=10, model="l2"):
    """
    Equivalent modern del teu codi amb ruptures

    data: DataFrame amb index datetime
    param: nom columna
    nb_data_min: mida mínima de segment
    pen: penalització (controla nº breakpoints)
    model: "l2", "rbf", "linear"
    """

    data = data.copy()
    data = data.dropna(subset=[param])

    signal = data[param].values

    # 🔍 Detecció global
    algo = rpt.Pelt(model=model, min_size=nb_data_min).fit(signal)
    bkpts = algo.predict(pen=pen)

    results = []

    start = 0
    level = 1  # simulació del "nivell" del teu arbre

    for b in bkpts:
        segment = data.iloc[start:b]
        x = segment[param].values

        if len(x) >= nb_data_min:

            # 📊 percentils (equivalent al teu PrctDiff parcial)
            p10, p50, p90 = np.percentile(x, [10, 50, 90])

            # 📍 breakpoint (últim punt del segment)
            time_bp = segment.index[-1]

            # 📈 diferència respecte segment anterior
            if len(results) > 0:
                prev = results[-1]
                prev_vals = np.array([prev["p10"], prev["p50"], prev["p90"]])
                curr_vals = np.array([p10, p50, p90])

                PrctDiff = np.round((curr_vals - prev_vals) * 100 / prev_vals)
            else:
                PrctDiff = [np.nan, np.nan, np.nan]

            results.append({
                "level": level,
                "time": time_bp,
                "pvalue": np.nan,          # ruptures no dona p-value directe
                "pvalue_boot": np.nan,     # no necessari
                "PrctDiff": PrctDiff,
                "p10": p10,
                "p50": p50,
                "p90": p90,
                "size": len(x)
            })

        start = b
        level += 1

    result_df = pd.DataFrame(results)

    return result_df