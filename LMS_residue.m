function data_residue=LMS_residue(data, param, distribution)

%LMS fit is applied on the montly mean (median if the distribution in 'log') of the data
%The LMS consists of:
%a constant, a linear slope, seasonalit�
% if not specified, end_year==2017 and the trend are computed for 10, 15,
% 20 and 25 years


%IN
% data: timetable
% param: parameter to  be analysied
%distribution: 'lin' or 'log'
%



%OUT: data_residue= residue of the fit
% programmed selon Weatherheat, JGR 1998 et 2000
%Martine Collaud Coen, 3.2026


% Set values from user input, or use defaults
%avM = '@nanmedian';
%end_time = max(year(data.Time));

% if necessary, take the logarithm of the data
if strcmp(distribution,'log')== 1
    data.(param)=log(data.(param));
end
data.(param)(data.(param)==-inf)=NaN;
data.(param)=real(data.(param));
ind=~isnan(data.(param));
datax=data(ind,:);

 
    %defini le vecteur temps centered
    t=datenum(datax.Time)-datenum(datax.Time(1));
    
    %defini la matrice X contenant toutes les dependences en t du fit
    E= [ones(size(t)) t sin((2*pi/365.25).*t) sin((4*pi/365.25).*t) sin((8*pi/365.25).*t) cos((2*pi/365.25).*t) cos((4*pi/365.25).*t) cos((8*pi/365.25).*t)];
b=E\datax.(param);
    Eplot=E*b;
    x_residue=datax.(param)-Eplot;
    data_residue=NaN(size(data.(param)));
    data_residue(ind)=x_residue;
    
 



