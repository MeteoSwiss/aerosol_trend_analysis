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

def compute_exp_D(data, shortnames, lam):
    """
    if scattering data are given, compute the pair-wise expSC and BF
    if absorption data are given, compute the expAE, pair-wise for 3w CLAP or
    PSAP & 2wAE
    and by a fit for  7w AE
    if there is several nephelometers with different lambdas, lambda is a
    matrix with 4 columns corresponding to the 4 possibility (0, 1, 2 or
    without number)
    """
    
    data_cal = pd.DataFrame(index=data.index)
    if "Time" in data.columns:
        data_cal["Time"] = data["Time"]

    namesOk = data.columns.tolist()
    
    namesB = [n for n in namesOk if any(s in n for s in shortnames)]
    namesBs = [n for n in namesB if "Bs" in n]
    namesBbs = [n for n in namesB if "Bbs" in n]
    namesBa = [n for n in namesB if "Ba" in n]

    # compute expS and BF
    if len(namesBs) > 0:
        if len(namesBs) == 1 and len(namesBbs) == 1:
            data_cal['BbsFG'] = data[namesBbs[0]]/ data[namesBs[0]]
            if len(namesBbs) > 0 and len(namesBbs) != 1:
                raise TypeError('there is more backscat coef than scat coeff')
        
        elif len(namesBs) == 3:
            b = [x for x in namesBs if x.startswith("BsB")]
            g = [x for x in namesBs if x.startswith("BsG")]
            r = [x for x in namesBs if x.startswith("BsR")]

            data_cal['expS_bg'] = np.real(-np.log(data[b[0]] / data[g[0]])
                        / np.log(lam[0] / lam[1]))
            data_cal['expS_br'] = np.real(-np.log(data[b[0]] / data[r[0]])
                        / np.log(lam[0] / lam[2]))
            data_cal['expS_gr'] = np.real(-np.log(data[g[0]] / data[r[0]])
                        / np.log(lam[1] / lam[2]))
            
            if len(namesBbs) > 0:
                Bb = [x for x in namesBbs if x.startswith("BsB")]
                Bg = [x for x in namesBbs if x.startswith("BsG")]
                Br = [x for x in namesBbs if x.startswith("BsR")]
                data_cal['BsBFb'] = data[Bb[0]] / data[b[0]]
                data_cal['BsBFg'] = data[Bg[0]] / data[g[0]]
                data_cal['BsBFr'] = data[Br[0]] / data[r[0]]
        
        elif len(namesBs) == 6:

            Cs = [[("0_" in n) and ("dry" not in n) and ("Q" not in n) for n in namesBs],
                [("1_" in n) and ("dry" not in n) and ("Q" not in n) for n in namesBs]]

            for i in range(2):

                namesBsx = [n for n, keep in zip(namesBs, Cs[i]) if keep]

                b = [x for x in namesBsx if x.startswith("BsB")]
                g = [x for x in namesBsx if x.startswith("BsG")]
                r = [x for x in namesBsx if x.startswith("BsR")]

                data_cal[f"expS_bg{i}"] = np.real(-np.log(data[b[0]] / data[g[0]])
                      / np.log(lam[0] / lam[1]))

                data_cal[f"expS_br{i}"] = np.real(-np.log(data[b[0]] / data[r[0]])
                      / np.log(lam[0] / lam[2]))

                data_cal[f"expS_gr{i}"] = np.real(-np.log(data[g[0]] / data[r[0]])
                      / np.log(lam[1] / lam[2]))

        if len(namesBbs) > 0:

            Cbs = [["0_" in n for n in namesBbs],
                ["1_" in n for n in namesBbs]]

            for i in range(2):

                namesBbsx = [n for n, keep in zip(namesBbs, Cbs[i]) if keep]

                Bb = [x for x in namesBbsx if x.startswith("BbsB")]
                Bg = [x for x in namesBbsx if x.startswith("BbsG")]
                Br = [x for x in namesBbsx if x.startswith("BbsR")]

                namesBsx = [n for n, keep in zip(namesBs, Cs[i]) if keep]

                b = [x for x in namesBsx if x.startswith("BsB")]
                g = [x for x in namesBsx if x.startswith("BsG")]
                r = [x for x in namesBsx if x.startswith("BsR")]

                data_cal[f"BbsFb{i}"] = data[Bb[0]] / data[b[0]]
                data_cal[f"BbsFbg{i}"] = data[Bg[0]] / data[g[0]]
                data_cal[f"BbsFbr{i}"] = data[Br[0]] / data[r[0]]
    
    # absorption
    if len(namesBa) > 0:

        # 2 wv
        if len(namesBa) == 2:

            if len(lam) != 2:
                raise ValueError("The number of wavelengths is different from the number of absorption coefficients.")

            data_cal["expA"] = np.real(-np.log(data[namesBa[0]] / data[namesBa[1]])
                / np.log(lam[0] / lam[1]))

        # 3-wv PSAP or CLAP
        elif len(namesBa) == 3:

            if len(lam) != 3:
                raise ValueError("The number of wavelengths is different from the number of absorption coefficients.")

            ba = [x for x in namesBa if x.startswith("BaB")]
            ga = [x for x in namesBa if x.startswith("BaG")]
            ra = [x for x in namesBa if x.startswith("BaR")]

            data_cal["expA_bg"] = np.real(-np.log(data[ba[0]] / data[ga[0]])
                / np.log(lam[0] / lam[1]))

            data_cal["expA_br"] = np.real(-np.log(data[ba[0]] / data[ra[0]])
                / np.log(lam[0] / lam[2]))

            data_cal["expA_gr"] = np.real(-np.log(data[ga[0]] / data[ra[0]])
                / np.log(lam[1] / lam[2]))

        # 2 × 3 wv instruments
        elif len(namesBa) == 6:

            if len(lam) != 3:
                raise ValueError("The number of wavelengths is different from the number of absorption coefficients.")

            Ca = [["0_" in n for n in namesBa],
                ["1_" in n for n in namesBa]]

            for i in range(2):

                namesBax = [n for n, keep in zip(namesBa, Ca[i]) if keep]

                ba = [x for x in namesBax if x.startswith("BaB")]
                ga = [x for x in namesBax if x.startswith("BaG")]
                ra = [x for x in namesBax if x.startswith("BaR")]

                data_cal[f"expA_bg{i}"] = np.real(-np.log(data[ba[0]] / data[ga[0]])
                    / np.log(lam[0] / lam[1]))

                data_cal[f"expA_br{i}"] = np.real(-np.log(data[ba[0]] / data[ra[0]])
                    / np.log(lam[0] / lam[2]))

                data_cal[f"expA_gr{i}"] = np.real(-np.log(data[ga[0]] / data[ra[0]])
                    / np.log(lam[1] / lam[2]))

        # 7-wv AE
        elif len(namesBa) == 7:

            if len(lam) != 7:
                raise ValueError("The number of wavelengths is different from the number of absorption coefficients.")

            b = [x for x in namesBa if x.startswith("Ba2")]
            g = [x for x in namesBa if x.startswith("Ba3")]
            r = [x for x in namesBa if x.startswith("Ba5")]

            data_cal["expA_bgAE"] = np.real(-np.log(data[b[0]] / data[g[0]])
                         / np.log(470 / 520))

            data_cal["expA_brAE"] = np.real(-np.log(data[b[0]] / data[r[0]])
                         / np.log(470 / 660))

            data_cal["expA_grAE"] = np.real(-np.log(data[g[0]] / data[r[0]])
                         / np.log(520 / 660))

            a1 = [x for x in namesBa if x.startswith("Ba1")]
            a4 = [x for x in namesBa if x.startswith("Ba4")]
            a6 = [x for x in namesBa if x.startswith("Ba6")]
            a7 = [x for x in namesBa if x.startswith("Ba7")]

            all_abs = np.column_stack([
                data[a1[0]],
                data[b[0]],
                data[g[0]],
                data[a4[0]],
                data[r[0]],
                data[a6[0]],
                data[a7[0]]])

            data_cal["expA_fit"] = Exp_AE_nan_regr(all_abs)

        else:
            raise ValueError("The number of absorption coefficients does not allow the instrument type to be identified.")
        
    return data_cal

