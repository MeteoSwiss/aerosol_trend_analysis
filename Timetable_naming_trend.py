import numpy as np
import pandas as pd


def timetable_naming_trend(ds):

    data = {}

    data["Time"] = pd.to_datetime(ds["time"].values)

    for var_name in ds.data_vars:

        da = ds[var_name]
        
        if da.dims == ("d_Wavelength", "time"):

            wavelength = ds["d_Wavelength"].values

            da = da.transpose("time", "d_Wavelength")

            values = da.values
    
            # keep variable wavelengths
            mask = ~np.isnan(values).all(axis=0)
    
            wavelength = wavelength[mask]
            values = values[:, mask]
    
            for j, wl in enumerate(wavelength):
            
                col = values[:, j]
    
                # keep the paramerter only if there is data for the applied wavelength
                if np.isnan(col).all():
                    continue
                
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
                wl_map = {370: "1",470: "2",520: "3",
                        590: "4",660: "5",880: "6",950: "7",}
                
                if len(wavelength) in [3, 7, 10]:
                
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
    
                elif set(wl_map).issubset(set(wavelength)):
                    if wl in wl_map:
                        wv = wl_map[wl]
                    else:
                        wv = int(wl)
                else:
                    wv = int(wl)
    
                name = f"{prefix}{wv}{pm}_{inst}"
    
                data[name] = col
            
            # save overlapped 660nm data in a new column
            if prefix =='Ba':
                for i in [0,1,2]:
                    if (f'Ba1{i}_A' in data and f'Ba5{i}_A' in data and
                        np.count_nonzero(~np.isnan(data[f'Ba1{i}_A'])) < 
                        np.count_nonzero(~np.isnan(data[f'Ba5{i}_A']))):
    
                        mask = np.isnan(data[f'Ba1{i}_A'])
    
                        data[f'Ba660{i}_A'] = np.where(mask, data[f'Ba5{i}_A'], np.nan)
                        data[f'Ba5{i}_A'] = np.where(mask, np.nan, data[f'Ba5{i}_A'])

    df = pd.DataFrame(data)

    df = df.set_index("Time")

    return df