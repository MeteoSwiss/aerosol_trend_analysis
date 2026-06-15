import pandas as pd
import numpy as np
import re

def timetable2num(df):
    """
    Converteix les columnes de text d'un DataFrame a valors numèrics
    quan sigui possible.
    """

    df_out = df.copy()

    # Expressió regular equivalent a la de MATLAB
    regexstr = (
        r'(?P<prefix>.*?)'
        r'(?P<numbers>('
        r'[-]*(\d+[\,]*)+[\.]?\d*[eEdD]?[-+]*\d*'
        r'|'
        r'[-]*(\d+[\,]*)*[\.]\d+[eEdD]?[-+]*\d*'
        r'))'
        r'(?P<suffix>.*)'
    )

    for col_name in df.columns:

        col = df[col_name]

        # Només tractem columnes de text
        if not (col.dtype == object or pd.api.types.is_string_dtype(col)):
            continue

        numeric_col = np.full(len(col), np.nan)

        for i, value in enumerate(col):

            if pd.isna(value):
                continue

            try:
                value = str(value)

                result = re.match(regexstr, value)

                if result is None:
                    continue

                numbers = result.group("numbers")

                # Tractament de separadors de milers
                if "," in numbers:
                    thousands_regex = r'^\d+?(,\d{3})*\.?\d*$'

                    if not re.match(thousands_regex, numbers):
                        continue

                    numbers = numbers.replace(",", "")

                numeric_col[i] = float(numbers)

            except Exception:
                continue

        df_out[col_name] = numeric_col

    return df_out