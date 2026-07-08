import pandas as pd
import numpy as np

def compute_SSA_D(data, shortnamesSC, shortnamesabs):
    """
    compute the SSA from the given shortnames

    if there is several nephelometers with different wavelengths, lambda is a
    matrix with 4 columns corresponding to the 4 possibility (0, 1, 2 or
    without number)
    takes names
    """
    
    data_cal = pd.DataFrame(index=data.index)

    namesOK = data.columns.tolist()
    
    namesB = [n for n in namesOK if any(s in n for s in shortnamesSC)]
    namesBs = [n for n in namesB if "Bs" in n]
    
    namesB = [n for n in namesOK if any(s in n for s in shortnamesabs)]
    namesBa = [n for n in namesB if "Ba" in n]

    # compute SSAB
    if len(namesBs) > 0 and len(namesBa) > 0:
        if len(namesBs) == 1 and len(namesBa) == 1:
            data_cal['SSA'] = data[namesBs[0]] / (data[namesBs[0]] + data[namesBa[0]])
        
        # trhee wavelength + 1 size cuts
        elif len(namesBs) == 2 and len(namesBa) == 2:
            
            bsG1 = [x for x in namesBs if x.startswith("BsG1")]
            bsG0 = [x for x in namesBs if x.startswith("BsG0")]
            baG1 = [x for x in namesBa if x.startswith(("Ba3", "Bac3", "Bax3", "BaG1"))]
            baG0 = [x for x in namesBa if x.startswith(("Ba3", "Bac3", "Bax3", "BaG0"))]

            data_cal["SSAG"] = data[bsG1[0]] / (data[bsG1[0]] + data[baG1[0]])
            data_cal["SSAR"] = data[bsG0[0]] / (data[bsG0[0]] + data[baG0[0]])

        # trhee wavelength + 1 size cuts
        elif len(namesBs) == 3 and len(namesBa) == 3: 

            bsB = [x for x in namesBs if x.startswith("BsB")]
            bsG = [x for x in namesBs if x.startswith("BsG")]
            bsR = [x for x in namesBs if x.startswith("BsR")]

            baB = [x for x in namesBa if x.startswith(("Ba2", "Bac2", "Bax2", "BaB"))]
            baG = [x for x in namesBa if x.startswith(("Ba3", "Bac3", "Bax3", "BaG"))]
            baR = [x for x in namesBa if x.startswith(("Ba5", "Bac5", "Bax5", "BaR"))]

            data_cal["SSAB"] = data[bsB[0]] / (data[bsB[0]] + data[baB[0]])
            data_cal["SSAG"] = data[bsG[0]] / (data[bsG[0]] + data[baG[0]])
            data_cal["SSAR"] = data[bsR[0]] / (data[bsR[0]] + data[baR[0]])

        # trhee wavelength + two size cuts
        elif len(namesBs) == 6 and len(namesBa) == 6:

            bsB0 = [x for x in namesBs if x.startswith(("BsB_", "BsB0"))]
            bsG0 = [x for x in namesBs if x.startswith(("BsG_", "BsG0"))]
            bsR0 = [x for x in namesBs if x.startswith(("BsR_", "BsR0"))]

            baB0 = [x for x in namesBa if x.startswith(("Ba20", "Bac20", "Bax20", "BaB0", "Ba2_", "Bac2_", "Bax2_", "BaB_"))]
            baG0 = [x for x in namesBa if x.startswith(("Ba30", "Bac30", "Bax30", "BaG0", "Ba3_", "Bac3_", "Bax3_", "BaG_"))]
            baR0 = [x for x in namesBa if x.startswith(("Ba50", "Bac50", "Bax50", "BaR0", "Ba5_", "Bac5_", "Bax5_", "BaR_"))]

            data_cal["SSAB0"] = data[bsB0[0]] / (data[bsB0[0]] + data[baB0[0]])
            data_cal["SSAG0"] = data[bsG0[0]] / (data[bsG0[0]] + data[baG0[0]])
            data_cal["SSAR0"] = data[bsR0[0]] / (data[bsR0[0]] + data[baR0[0]])

            bsB1 = [x for x in namesBs if x.startswith("BsB1")]
            bsG1 = [x for x in namesBs if x.startswith("BsG1")]
            bsR1 = [x for x in namesBs if x.startswith("BsR1")]

            baB1 = [x for x in namesBa if x.startswith(("Ba21", "Bac21", "Bax21", "BaB1"))]
            baG1 = [x for x in namesBa if x.startswith(("Ba31", "Bac31", "Bax31", "BaG1"))]
            baR1 = [x for x in namesBa if x.startswith(("Ba51", "Bac51", "Bax51", "BaR1"))]

            data_cal["SSAB1"] = data[bsB1[0]] / (data[bsB1[0]] + data[baB1[0]])
            data_cal["SSAG1"] = data[bsG1[0]] / (data[bsG1[0]] + data[baG1[0]])
            data_cal["SSAR1"] = data[bsR1[0]] / (data[bsR1[0]] + data[baR1[0]])

        else:
            raise TypeError("The number of shortnames does not allow the instrument type to be identified.")
        
    return data_cal