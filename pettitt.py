import numpy as np

"""This code is used to find the change point in a univariate continuous time series
    using Pettitt Test.


    The test here assumed is two-tailed test. The hypothesis are as follow:
    H (Null Hypothesis): There is no change point in the series
    H(Alternative Hypothesis): There is a change point in the series

    Input: univariate data series
    Output:
    The output of the answer in row wise respectively,
    loc: location of the change point in the series, index value in
    the data set
    K: Pettitt Test Statistic for two tail test
    pvalue: p-value of the test

    Reference: Pohlert, Thorsten. "Non-Parametric Trend Tests and Change-Point Detection." (2016).
    """
def pettitt(data, alpha):
    data = np.asarray(data)
    m = len(data)

    # Matriu de comparacions
    #with vectorisation
    t1 = np.tile(data, (m, 1))
    v = np.sign(t1.T - t1)

    # Suma per files
    V = np.sum(v, axis=1)

    # Suma acumulada
    U = np.cumsum(V)

    # Localització del canvi
    K = np.max(np.abs(U))
    locs = np.where(np.abs(U) == K)[0]

    if len(locs) > 1:
        d = m - locs
        loc = locs[np.argmin(d)]
    elif len(locs) == 1:
        loc = locs[0]
    else:
        loc = np.nan

    # p-value
    pvalue = 2 * np.exp((-6 * K**2) / (m**3 + m**2))

    if np.isnan(loc):
        return np.array([np.nan, K, np.nan]), [np.nan, np.nan, np.nan]

    if pvalue <= alpha:
        before = data[:loc]
        after = data[loc:]
        #compute difference after-before the break for min (10%), med (50%) and max (90%) percentiles
        p_before = np.percentile(before, [10, 50, 90])
        p_after = np.percentile(after, [10, 50, 90])
        p_total = np.percentile(data, [10, 50, 90])

        PrctDiff = np.round((p_after - p_before) * 100 / p_total)

        return np.array([loc, K, pvalue]), PrctDiff

    else:
        return np.array([np.nan, K, np.nan]), [np.nan, np.nan, np.nan]