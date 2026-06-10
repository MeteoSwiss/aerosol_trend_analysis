import xarray as xr
import pandas as pd
import numpy as np
from Timetable_naming_trend import timetable_naming_trend

def ebas_nc_to_timetable(filename):
    """
    Read an EBAS NetCDF file and convert it to a pandas DataFrame.

    Parameters

    ----------
    filename : str

    Returns
    -------
    stn_rd : pandas.DataFrame
        Daily aggregated data.
    stn_st : dict
        Station metadata.
    """

    # Open NetCDF
    ds = xr.open_dataset(filename)

    # Modify the parameter's name according to the naming convention
    stn_rd = timetable_naming_trend(ds)  

    # compute the daily data
    stn_rd = stn_rd.resample("D").median()

    # take information on the station
    def get_attr_contains(attrs, text):
        for key, value in attrs.items():
            if text.lower() in key.lower():
                return value
        return None

    stn_st = {}

    attrs = ds.attrs

    stn_st = {
        "name": (get_attr_contains(attrs, "ebas_station_gaw_id")
            or get_attr_contains(attrs, "ebas_station_name")),
        "lat": get_attr_contains(attrs, "latitude"),
        "lon": get_attr_contains(attrs, "longitude"),
        "alt": get_attr_contains(attrs, "altitude"),
        "env": get_attr_contains(attrs, "ebas_station_setting"),
        "footp":get_attr_contains(attrs, "ebas_station_land_use")}

    ds.close()

    return stn_rd, stn_st