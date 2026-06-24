import pandas as pd
import numpy as np

def Exp_AE_nan_regr(data):

    lams = np.array([370, 470, 520, 590, 660, 880, 950], dtype=float)
    x = np.log(lams)

    data = np.asarray(data, dtype=float)

    # Evitar logs de valors <= 0
    data[data <= 0] = np.nan

    y = np.log(data)

    aae = np.full(len(y), np.nan)

    for i in range(len(y)):
        mask = np.isfinite(y[i])

        # mínim 2 punts per fer regressió
        if np.sum(mask) >= 2:
            slope, intercept = np.polyfit(x[mask], y[i,mask], 1)
            aae[i] = -slope

    return aae

def compute_exp_ssa(data, lambdaSC, lambdaAE):
    """
    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame amb columnes tipus:
        BsB, BsG, BsR, BbsB, BbsG, BbsR

    lambdaSC : list or array
        Longituds d'ona [blue, green, red]

    Returns
    -------
    data_cal : pandas.DataFrame
    """

    data_cal = pd.DataFrame(index=data.index)

    names = data.columns.tolist()

    # =========================================================
    # SCATTERING VARIABLES
    # =========================================================

    namesBs = [c for c in names if c.startswith("Bs")]

    if len(namesBs) > 0:

        # GREEN SCATTERING
        namesBsG = [c for c in namesBs if "G" in c]

        if len(namesBsG) == 0:
            namesBsG = [namesBs[0]]

        # BACKSCATTERING
        namesBbs = [c for c in names if c.startswith("Bbs")]

        if len(namesBbs) > 0:

            namesBbsG = [c for c in namesBbs if "G" in c]

            if len(namesBbsG) == 0:
                namesBbsG = [namesBbs[0]]

            # Backscattering fraction
            if len(namesBsG) == 1 and len(namesBbsG) == 1:

                data_cal["BbsFG"] = (data[namesBbsG[0]] / data[namesBsG[0]])
            
            else:

                Cx = {} # CLASSIFY BsG CHANNELS

                Cx[0] = [("0_" in c) and ("dry" not in c) and ("Q" not in c)
                    for c in namesBsG]

                Cx[1] = [("1_" in c) and ("dry" not in c) and ("Q" not in c)
                    for c in namesBsG]

                Cx[2] = [("dry" in c)
                    for c in namesBsG]

                Cx[3] = [("_" in c) and not Cx[0][i] and not Cx[1][i] and not Cx[2][i]
                    for i, c in enumerate(namesBsG)]

                Cy = {} # CLASSIFY BbsG CHANNELS

                Cy[0] = [("0_" in c) and ("dry" not in c) and ("Q" not in c)
                    for c in namesBbsG]

                Cy[1] = [("1_" in c) and ("dry" not in c) and ("Q" not in c)
                    for c in namesBbsG]

                Cy[2] = [("dry" in c)
                    for c in namesBbsG]

                Cy[3] = [("_" in c) and not Cx[0][i] and not Cx[1][i] and not Cx[2][i]
                    for i, c in enumerate(namesBbsG)]
                
                for i in range(4):

                    # how many channels found
                    e = sum(Cx[i])
                    r = sum(Cy[i])

                    # only if exactly one Bs and one Bbs
                    if e == 1 and r == 1:

                        # get matching names
                        namesBsx = [namesBsG[j] for j in range(len(namesBsG))
                            if Cx[i][j]]

                        namesBbsx = [namesBbsG[j] for j in range(len(namesBbsG))
                            if Cy[i][j]]

                        # output variable name
                        if i == 3:
                            namesFG = "BbsFG"
                        else:
                            namesFG = f"BbsFG{i}"

                        # compute backscattering fraction
                        data_cal[namesFG] = (data[namesBbsx[0]] / data[namesBsx[0]])
        
        namesBsG=namesBs

        # COMPUTE SCATTERING ANGSTROM EXPONENT

        # CLASSIFY Bs CHANNELS
        Cs = {} 

        Cs[0] = [("0_" in c) and ("dry" not in c) and ("Q" not in c)
            for c in namesBs]

        Cs[1] = [("1_" in c) and ("dry" not in c) and ("Q" not in c)
            for c in namesBs]

        Cs[2] = [("dry" in c) and ("Q" not in c)
            for c in namesBs]

        Cs[3] = [("_" in c) and ("Q" not in c) and not Cs[0][i] and not Cs[1][i] and not Cs[2][i]
            for i, c in enumerate(namesBs)]
        
        # LOOP OVER CONFIGURATIONS        
        for i in range(4): 

            namesBsx = [namesBs[j] for j in range(len(namesBs)) if Cs[i][j]]

            if len(namesBsx) >= 3: # THREE WAVELENGTHS

                b = [c for c in namesBsx if c.startswith("BsB")]
                g = [c for c in namesBsx if c.startswith("BsG")]
                r = [c for c in namesBsx if c.startswith("BsR")]
                
                if len(b) > 0 and len(g) > 0: # BLUE-GREEN

                    data_cal[f"expS_bg{i}"] = (-np.log(data[b[0]] / data[g[0]])
                        / np.log(lambdaSC[0] / lambdaSC[1]))

                if len(b) > 0 and len(r) > 0: # BLUE-RED

                    data_cal[f"expS_br{i}"] = (-np.log(data[b[0]] / data[r[0]])
                        / np.log(lambdaSC[0] / lambdaSC[2]))

                if len(g) > 0 and len(r) > 0: # GREEN-RED

                    data_cal[f"expS_gr{i}"] = (-np.log(data[g[0]] / data[r[0]])
                        / np.log(lambdaSC[1] / lambdaSC[2]))
            
            elif len(namesBsx) == 2: # ONLY TWO WAVELENGTHS

                b = [c for c in namesBsx if c.startswith("BsB")]
                g = [c for c in namesBsx if c.startswith("BsG")]

                if len(b) > 0 and len(g) > 0:

                    data_cal[f"expS_bg{i}"] = (-np.log(data[b[0]] / data[g[0]])
                        / np.log(lambdaSC[0] / lambdaSC[1]))
    
    else:
        namesBsG=namesBs

    # =========================================================
    # ABSORPTION VARIABLES
    # =========================================================

    namesBa = [c for c in names if c.startswith("Ba")]

    if len(namesBa) > 0:

        # FIND RED CHANNELS
        namesBaR = [c for c in namesBa if ("R" in c) or ("5" in c)]

        if len(namesBaR) == 0:
            namesBaR = [namesBa[0]]

        # AE33 / 7 WAVELENGTHS
        namesBa7 = [c for c in namesBa if "Ba7" in c]

        # CLASSIFY Bs CHANNELS
        Ca = {}
        Ca[0] = [("0_" in c) for c in namesBa]
        Ca[1] = [("1_" in c) for c in namesBa]
        Ca[2] = [("dry" in c) for c in namesBa]
        Ca[3] = [("_" in c) and not Ca[0][i] and not Ca[1][i] and not Ca[2][i]
            for i, c in enumerate(namesBa)]
        
        for i in range(4):
            namesBax = [namesBa[j] for j in range(len(namesBa))
                    if Ca[i][j]]

            if len(namesBa7) > 0 and len(namesBax) == 7:

                b = [c for c in namesBax if c.startswith("Ba2")]
                g = [c for c in namesBax if c.startswith("Ba3")]
                r = [c for c in namesBax if c.startswith("Ba5")]

                # AAE
                data_cal[f"expA_bg{i}"] = (-np.log(data[b[0]] / data[g[0]])
                    / np.log(470 / 520))

                data_cal[f"expA_br{i}"] = (-np.log(data[b[0]] / data[r[0]])
                    / np.log(470 / 660))

                data_cal[f"expA_gr{i}"] = (-np.log(data[g[0]] / data[r[0]])
                    / np.log(520 / 660))

                a1 = [c for c in namesBax if c.startswith("Ba1")]
                a4 = [c for c in namesBax if c.startswith("Ba4")]
                a6 = [c for c in namesBax if c.startswith("Ba6")]
                a7 = [c for c in namesBax if c.startswith("Ba7")]

                all_abs = np.column_stack([data[a1[0]],data[b[0]],data[g[0]],
                    data[a4[0]],data[r[0]],data[a6[0]],data[a7[0]]])

                data_cal[f"expA_fit{i}"] = Exp_AE_nan_regr(all_abs)

            # MIXED AE + 3W INSTRUMENTS
            elif len(namesBax) in [10, 13, 16]:

                print("Warning: one AE and one/several 3w instrument?")

                # AE33 PART
                b = [c for c in namesBax if c.startswith("Ba2")]
                g = [c for c in namesBax if c.startswith("Ba3")]
                r = [c for c in namesBax if c.startswith("Ba5")]
                a1 = [c for c in namesBax if c.startswith("Ba1")]
                a4 = [c for c in namesBax if c.startswith("Ba4")]
                a6 = [c for c in namesBax if c.startswith("Ba6")]
                a7 = [c for c in namesBax if c.startswith("Ba7")]

                data_cal[f"expA_AE_bg{i}"] = (-np.log(data[b[0]] / data[g[0]])
                    / np.log(470 / 520))

                data_cal[f"expA_AE_br{i}"] = (-np.log(data[b[0]] / data[r[0]])
                    / np.log(470 / 660))

                data_cal[f"expA_AE_gr{i}"] = (-np.log(data[g[0]] / data[r[0]])
                    / np.log(520 / 660))

                all_abs = np.column_stack([data[a1[0]],data[b[0]],data[g[0]],
                    data[a4[0]],data[r[0]],data[a6[0]],data[a7[0]]])
                data_cal[f"expA_fit{i}"] = Exp_AE_nan_regr(all_abs)

                # CLASSIFY 3-WAVELENGTH INSTRUMENTS
                ba = [c for c in namesBax if c.startswith("BaB")]
                ga = [c for c in namesBax if c.startswith("BaG")]
                ra = [c for c in namesBax if c.startswith("BaR")]

                data_cal[f"expA_bg{i}"] = (-np.log(data[ba[0]] / data[ga[0]])
                    / np.log(lambdaAE[0] / lambdaAE[1]))
                
                data_cal[f"expA_br{i}"] = (-np.log(data[ba[0]] / data[ra[0]])
                    / np.log(lambdaAE[0] / lambdaAE[2]))
                
                data_cal[f"expA_gr{i}"] = (-np.log(data[ga[0]] / data[ra[0]])
                    / np.log(lambdaAE[1] / lambdaAE[2]))

            # 2 instruments x 2 size cuts x 3 wavelengths
            elif len(namesBax) == 12:

                Ca = {}
                Ca[0] = [("0_A11" in c) for c in namesBax]
                Ca[1] = [("1_A11" in c) for c in namesBax]
                Ca[2] = [("0_A12" in c) for c in namesBax]
                Ca[3] = [("1_A12" in c) for c in namesBax]

                for i in range(4):
                
                    namesBax = [namesBa[j] for j in range(len(namesBa))
                        if Ca[i][j]]

                    if len(namesBax) == 3:
                    
                        ba = [c for c in namesBax if c.startswith("BaB")]
                        ga = [c for c in namesBax if c.startswith("BaG")]
                        ra = [c for c in namesBax if c.startswith("BaR")]

                        data_cal[f"expA_bg{i}"] = (-np.log(data[ba[0]] / data[ga[0]])
                            / np.log(lambdaAE[0] / lambdaAE[1]))

                        data_cal[f"expA_br{i}"] = (-np.log(data[ba[0]] / data[ra[0]])
                            / np.log(lambdaAE[0] / lambdaAE[2]))

                        data_cal[f"expA_gr{i}"] = (-np.log(data[ga[0]] / data[ra[0]])
                            / np.log(lambdaAE[1] / lambdaAE[2]))

            else:
                if len(namesBax) == 3:
                
                    ba = [c for c in namesBax if c.startswith("BaB")]
                    ga = [c for c in namesBax if c.startswith("BaG")]
                    ra = [c for c in namesBax if c.startswith("BaR")]
                    data_cal[f"expA_bg{i}"] = (-np.log(data[ba[0]] / data[ga[0]])
                        / np.log(lambdaAE[0] / lambdaAE[1]))
                    data_cal[f"expA_br{i}"] = (-np.log(data[ba[0]] / data[ra[0]])
                        / np.log(lambdaAE[0] / lambdaAE[2]))
                    data_cal[f"expA_gr{i}"] = (-np.log(data[ga[0]] / data[ra[0]])
                        / np.log(lambdaAE[1] / lambdaAE[2]))            

    # =========================================================
    # COMPUTE SSA
    # =========================================================

    if len(namesBsG) > 0 and len(namesBa) > 0:

        # CLASSIFY ABSORPTION
        CaR = {}
        CaR[0] = [c.startswith("BaR0") for c in namesBaR]
        CaR[1] = [c.startswith("BaR1") for c in namesBaR]
        CaR[2] = [(("BacR" in c) or ("Ba3" in c)) for c in namesBaR]
        CaR[3] = [(c.startswith("BaR_") or "Ba5" in c) for c in namesBaR]

        CaG = [c.startswith("BaG_") for c in namesBa]

        # CLASSIFY SCATTERING
        CsG = {}
        CsG[0] = [c.startswith("BsG0") for c in namesBsG]
        CsG[1] = [c.startswith("BsG1") for c in namesBsG]
        CsG[2] = [(c.startswith("BsG") and c.endswith("dry")) for c in namesBsG]
        CsG[3] = [c.startswith("BsG_") for c in namesBsG]
        
        # LOOP OVER GROUPS
        for i in range(4):

            # CASE 1: matching BaR + BsG
            if sum(CaR[i]) > 0 and sum(CsG[i]) > 0:

                namesBsy = [namesBsG[j] for j in range(len(namesBsG)) if CsG[i][j]]
                namesBay = [namesBaR[j] for j in range(len(namesBaR)) if CaR[i][j]]

                for j in namesBsy:
                    for k in namesBay:

                        N_SSA = f"SSA{i}{j}{k}"

                        data_cal[N_SSA] = (data[j] / (data[j] + data[k]))

            # CASE 2: fallback Ba
            elif sum(CaG) > 0 and sum(CsG[i]) > 0:

                namesBsy = [namesBsG[j] for j in range(len(namesBsG)) if CsG[i][j]]
                namesBay = [namesBa[j] for j in range(len(namesBa)) if CaG[j]]

                for j in namesBsy:
                    for k in namesBay:

                        name = f"SSA{i}{j}{k}"

                        data_cal[name] = (data[j] / (data[j] + data[k]))

    return data_cal