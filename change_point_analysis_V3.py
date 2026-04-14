import numpy as np
import pandas as pd
import LMS_residue
import multiple_breakpoints_while_one
import multiple_breakpoints_2y_period

"""
 apply the Pettitt and monthly data
 the tests are applied only for at least 10 years on yearly data -->
 subsequent division of the time serie also follows this rule
 tests on monthly data should not be applied on less than 2 years if divided
 find a way to extract the best potentiel change points (and not all)
 save fig all in one
 inputs: data = TimeTable
           param = parameters (as fieldnames) to be analysed
           alpha = confidence level, set at 0.05 as default
 outputs:  result= table of results
 use: ...
 example:
 fev-mars 2026, mco"""

def change_point_analysis_v3(data, params, alpha, station, inst):

    data = data.copy()

    # remove lines with NaN
    data = data.dropna(subset=params)

    # produce yearly and monthly medians
    data_m = data.resample("M").median()

    data_m["month"] = data_m.index.month

    start_time = data_m.index.year.min()
    end_time = data_m.index.year.max()
    period = end_time - start_time + 1

    # deseasonlise with the median of the month
    month_med = data_m.groupby("month").median()

    data_m_deseason = pd.DataFrame(index=data_m.index)

    for p in params:
        data_m_deseason[p] = (
            data_m[p] /
            data_m["month"].map(month_med[p])
        )

    # deseasonalise and detrend by LMS on data
    data_m_residue = pd.DataFrame(index=data_m.index)
    # deseasonalise and detrend by LMS on log(data)
    data_m_residueLog = pd.DataFrame(index=data_m.index)

    for p in params:
        data_m_residue[p] = LMS_residue(data_m, p, "lin")
        data_m_residueLog[p] = LMS_residue(data_m, p, "log")

    results = []

    # 🔹 5. loop per paràmetres
    for p in params:

        # --- 1. deseason
        result_m = multiple_breakpoints_while_one(
            data_m_deseason, p, 12, alpha
        )

        results.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": "month",
            "parameter": f"{p}_deseason",
            "instrument": inst,
            "method": "Pettitt",
            "results": result_m
        })

        # --- 2. residue
        result_m = multiple_breakpoints_while_one(
            data_m_residue, p, 12, alpha
        )

        results.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": "month",
            "parameter": f"{p}_residue",
            "instrument": inst,
            "method": "Pettitt",
            "results": result_m
        })

        # --- 3. residue log
        result_m = multiple_breakpoints_while_one(
            data_m_residueLog, p, 12, alpha
        )

        results.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": "month",
            "parameter": f"{p}_residueLog",
            "instrument": inst,
            "method": "Pettitt",
            "results": result_m
        })

        # --- 4. backward (invertit)
        data_inv = data_m_residueLog[p].iloc[::-1].copy()

        df_inv = pd.DataFrame({p: data_inv})

        result_m = multiple_breakpoints_while_one(
            df_inv, p, 12, alpha
        )

        # corregir temps
        corrected_times = []
        for t in result_m["time"]:
            if pd.notna(t):
                pos = data_m_residueLog.index.get_loc(t)
                N = len(data_m_residueLog)
                corrected_times.append(data_m_residueLog.index[N - pos - 1])
            else:
                corrected_times.append(pd.NaT)

        result_m["time"] = corrected_times

        results.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": "month",
            "parameter": f"{p}_residueLog_inv",
            "instrument": inst,
            "method": "Pettitt",
            "results": result_m
        })

        # --- 5. finestres 2 anys
        result_m = multiple_breakpoints_2y_period(
            data_m_residueLog, p, 12, alpha
        )

        results.append({
            "station": station,
            "start_time": start_time,
            "end_time": end_time,
            "length_period": period,
            "granularity": "month",
            "parameter": f"{p}_residueLog_2yper",
            "instrument": inst,
            "method": "Pettitt",
            "results": result_m
        })

    # 🔹 convertir a DataFrame final
    Tresult = pd.DataFrame(results)

    return Tresult