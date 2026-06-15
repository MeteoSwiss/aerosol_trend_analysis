import numpy as np
import pandas as pd


def timetable_naming_trend(ds):

    data = {}

    data["Time"] = pd.to_datetime(ds["time"].values)

    wavelength = ds["d_Wavelength"].values

    for var_name in ds.data_vars:

        da = ds[var_name]

        if "d_Wavelength" in da.dims:
            da = da.transpose("time", "d_Wavelength")

        values = da.values

        # variables amb dimensió longitud d'ona
        if values.ndim < 2:
            continue

        for j in range(values.shape[1]):

            # select variable type
            if "_scattering" in var_name:
                prefix = "Bs"
                inst = "S"

            elif "backscattering" in var_name:
                prefix = "Bbs"
                inst = "S"

            elif "absorption" in var_name:
                prefix = "Ba"
                inst = "A"

            elif "humidity" in var_name:
                prefix = "U"

            else:
                continue

            # select size cut
            if "pm10" in var_name:
                pm = "0"

            elif "pm2.5" in var_name:
                pm = "2"

            elif "pm1" in var_name:
                pm = "1"

            else:
                pm = ""

            # select wavelength
            wl = wavelength[j]

            wl_map = {370: "1",470: "2",520: "3",
                590: "4",660: "5",880: "6",950: "7",}

            if wl in wl_map:
                wv = wl_map[wl]

            else:

                if wl < 500:
                    wv = "B"
                elif wl < 600:
                    wv = "G"
                elif wl < 720:
                    wv = "R"
                else:
                    wv = "Q"

            col = values[:, j]

            # keep the paramerter only if there is data for the applied wavelength
            if np.isnan(col).all():
                continue

            name = f"{prefix}{wv}{pm}_{inst}"

            data[name] = col

    df = pd.DataFrame(data)

    df = df.set_index("Time")

    return df