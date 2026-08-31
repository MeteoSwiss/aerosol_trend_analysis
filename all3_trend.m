function [Tresult_MK,Tresult_LMSlog,Tresult_LMSlin]=all3_trend(data,param, inst, station, resolution, varargin)

% do the MK, LM and MK trend and save the results and figures
%data: should be as time table

%Dependency: call seasonalKendall_main_D & trend_LMS_D

Tresult_MK=seasonalKendall_main_D(data,param, inst, station, resolution,varargin{:});
%Tresult_GSMm=bootstrapconfidence_seas_D(data, param,inst, station, 'lin','granu','monthly',varargin{:});
%Tresult_GSMd=bootstrapconfidence_seas_D(data, param,inst, station, 'lin','granu','daily',varargin{:});
Tresult_LMSlog=trend_LMS_D(data, param,inst, station, 'log',varargin{1:end-2},'fig',0);
%Tresult_LMSlog2=trend_LMS_deseasonFit(data, param,inst, station, 'log',varargin{:});

Tresult_LMSlin=trend_LMS_D(data, param,inst, station, 'lin',varargin{:});
%Tresult_LMSlin2=trend_LMS_deseasonFit(data, param,inst, station, 'lin',varargin{:});
%Tresult=[Tresult_MK ;Tresult_GSMd;Tresult_LMS;Tresult_LMS2];


% old version: [Tresult_MK,Tresult_GSMd,Tresult_LMSlog,Tresult_LMSlin,Tresult_LMSlog2,Tresult_LMSlin2]=all3_trend(data,param, inst, station, resolution, varargin)