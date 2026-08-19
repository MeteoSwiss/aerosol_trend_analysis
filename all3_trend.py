from trend_LMS_D import trend_LMS_D
from seasonalKendall_main_D import seasonalKendall_main_D

def all3_trend(data,param,inst,station,resolution,**kwargs):
    """ do the MK, LM and MK trend and save the results and figures
    data: should be as time table"""

    #Dependency: call seasonalKendall_main_D & trend_LMS_D
    # Mann-Kendall seasonal trend
    Tresult_MK = seasonalKendall_main_D(data,param,inst,station,resolution,**kwargs)

    # LMS logarithmic
    kwargs_log = kwargs.copy()
    kwargs_log["fig"] = False
    Tresult_LMSlog, _, _  = trend_LMS_D(data,param,inst,station,"log",**kwargs_log)

    # LMS linear
    Tresult_LMSlin, _, _  = trend_LMS_D(data,param,inst,station,"lin",**kwargs)

    return (Tresult_MK,Tresult_LMSlog,Tresult_LMSlin)