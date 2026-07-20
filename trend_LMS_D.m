function [Tresult,f,data_residue]=trend_LMS_D(data, param,inst, station, distribution,varargin)

%LMS fit is applied on the montly mean (median if the distribution in 'log') of the data
%The LMS consists of:
% a constant, a linear slope, seasonality
% if not specified, end_year==2025
% the trend are computed for the whole time series and complete decades


%IN
% data:data either as a structure with time in data.start_time
%                     or a timetable
%                     or a matrice with the time as the first 6 columns (datevec)
% param = cell: either sturcture or timetable name or number of matrix column
% instrument= specify the instrument
%distribution: 'lin' or 'log'
%
%varargin:
%averageM= average method {'mean',@median}
%       default=@nanmedian
% period= period for trend analysis, default: 10, 20, 30
% end_year, default: 2025
%fig: if 1, plot and save figure


%OUT: Tresult= table of results with
%               station= name of the measuring site
%               end_time= end year used for the analyis
%               length_period= analysed period
%               granularity= used average of the data
%               parameter = name of the analysed parameter
%               instrument = instrument measuring the parameter
%               MK_seasonality = define if the Mann-Kendall is applied on
%                                the whole year, on 4 meteorological seasons or on 12 months
%               method = trend analysis method
%               ss= statistical significance: 95=95% if the significance (Weatherhead) >2
%                                     90=90% if the significance(Weatherhead) >1.67
%               slope: slope in units/y with log if request
%               UCL: upper confidence level in units/y with log if request
%               LCL: lower confidence level in units/y with log if request
%               median: median of the data with log if request
%               slopeP: Sen's slope in %/y with log if request
%               UCLP: upper confidence level in %/y with log if request
%               LCLP: lower confidence level in %/y with log if request
%               slopeR: Sen's slope in % for the whole analysed period without log
%               UCLR: upper confidence level in % for the whole analysed period without log
%r              LCLR: lower confidence level in % for the whole analysed period without log
%               significance =significance= slope/variance
%               variance=variance from Weatherhead

% create a figure containing
% 1. monthly average, fitted data and trend, 2. residuals,
% 3. normplot of residuals, 4.cummulative summation of residues

% from Weatherheat, JGR 1998 et 2000
% Martine Collaud Coen, 7. 2026


%% ATTENTION, PATH FOR SAVING THE FIGURE IN LINES 2005-2009. PLEASE ADAPT FOR YOUR USAGE


% FOR JORDI: NOT NECESSARY IF YOU IMPOSE A FORMAT
if ischar(param{1})
    name=param{1};
elseif isnumeric(param{1})
    colonne=param{1};
end

% check arguments
if ~varg_proof(varargin, {'aveM', 'end_year','path','period','fig'},true)
    return
end

% Set values from user input, or use defaults
avM = varg_val(varargin, 'averageM', @nanmedian);
end_time = varg_val(varargin, 'end_year', 2017);
path = varg_val(varargin, 'path', '');
period = varg_val(varargin, 'period', [10 20 30 40 50]);
granu='month';
fig=varg_val(varargin, 'fig', 1);

%FOR JORDI: IF YOU IMPOSE A FORMAT, THIS IS ALSO NOT NEEDED TO TRANSLATE
%put the data in timetable
if isstruct(data)
    dataT=timetable(datetime(datevec(data.start_time)),data.(name)','VariableNames',{name});
elseif istimetable (data)
    dataT=timetable(data.Time,data.(name),'VariableNames',{name});
else
    dataT=timetable(datetime(data(:,1:6)), data(:,colonne),'VariableNames',{'parametre'});
    name=dataT.Properties.VariableNames{1};
end

dataT.(name)(isinf(dataT.(name)))=NaN;
% monthly average
dataG=retime(dataT,'monthly',avM);
dataG.(name)=real(dataG.(name));
dataG.y=year(dataG.Time);

% if necessary, take the logarithm of the data
if strcmp(distribution,'log')== 1
    dataG.(name)=real(log(dataG.(name)));
end

ind=~isnan(dataG.(name)) & ~isinf(dataG.(name));
dataG=dataG(ind,:);

%take the number of periods to be analysed
if length(period)>1
    if period(end)-period(end-1)<10
        nb_period=floor((max(dataG.y)-min(dataG.y)+1)/5)-1;
        if nb_period==0
            nb_period=1;
        end
    else
        nb_period=floor((max(dataG.y)-min(dataG.y)+1)/(period(end)-period(end-1)));
    end
else nb_period=1;
end

% compute the LMS trend for all periods
for i=nb_period:-1:1
    TT=timerange(datetime(end_time-period(i)+1,1,1,0,0,0),datetime(end_time,1,1,0,0,0));
    dataGa=dataG(TT,:);
    % center the time vector
    t=datenum(dataGa.Time)-datenum(dataGa.Time(1));

    % define the matrix X with all fit dependencies
    E= [ones(size(t)) t sin((2*pi/365.25).*t) sin((4*pi/365.25).*t) sin((8*pi/365.25).*t) cos((2*pi/365.25).*t) cos((4*pi/365.25).*t) cos((8*pi/365.25).*t)];

    b=E\dataGa.(name);

    result.slopeP=b(2)*365.25*100/abs(nanmedian(dataGa.(name)));

    % trend per year
    result.slope=b(2)*365.25;

    % compute residuals
    Eplot=E*b;
    x_residue=dataGa.(name)-Eplot;

    % compute autocorrelation coef. phi and variance of white noise
    [a,e]=arburg(x_residue,1);
    delta_b1=std(x_residue)/(sum(ind)+1)^0.5;

    % compute the statistical significance:
    %  95% of confidence level if real_T > 2,
    % 90% if real_T > 1.67
    result.variance=e.^0.5 / ((1+a(2))*period(i).^(3/2));

    %computation of the statistical significance
    result.significance=abs(b(2))*365.25 /result.variance;
    if result.significance>=2
        result.ss=95;
    elseif result.significance>=1.67
        result.ss=90;
    else
        result.ss=0;
    end

    % computation of the upper and lower confidence limits
    result.UCLP=(b(2)+2*result.variance/365.25)*365.25*100/abs(nanmedian(dataGa.(name)));
    result.LCLP=(b(2)-2*result.variance/365.25)*365.25*100/abs(nanmedian(dataGa.(name)));
    result.UCL=(b(2)+2*result.variance/365.25)*365.25;
    result.LCL=(b(2)-2*result.variance/365.25)*365.25;

    % compute the slope without log for all periods for log distribution
    %slope in % from the value at the beginning==1
    if strcmp('log', distribution)==1
        result.slopeR=(exp(result.slope*period(i))-1);
        result.UCLR=(exp((result.slope+2*result.variance/365.25)*period(i))-1);
        result.LCLR=(exp((result.slope-2*result.variance/365.25)*period(i))-1);
    elseif strcmp('lin',distribution)==1
        result.slopeR=(result.slope*period(i)-1);
        result.UCLR=((result.slope+2*result.variance/365.25)*period(i)-1);
        result.LCLR=((result.slope-2*result.variance/365.25)*period(i)-1);
    end

    %create the table of results and do the figure only for lin
    %distribution
    if i==nb_period
        if fig 
            f=figure(203);
            fig_LMS(dataGa,name,Eplot,x_residue,b);
            hold on;
            subplot(2,2,1);
            title(join([(station) ,(name),'LMS',{inst}]));
        end
        Tresult(i,:)=table({station}, end_time, period(i), {granu}, {name} ,{inst},  {distribution}, {'LMS'}, {result.significance},{result.ss},{result.slope}, {result.UCL},{result.LCL},{result.slopeP}, {result.UCLP},{result.LCLP},{result.slopeR}, {result.UCLR},{result.LCLR},...
            'VariableNames',{'station','end_time','length_period','granularity','parameter','instrument','MK_seasonality','method','significance','ss','slope','UCL','LCL','slopeP','UCLP','LCLP','slopeR','UCLR','LCLR'});
    else
        %add the slope of the other periods on the figure with the longest period
        if fig
            subplot(2,2,1)
            hold on;
            title(join([(station) ,(name),'LMS',{inst}]));

            plot(dataGa.Time, b(1)+b(2).*(datenum(dataGa.Time)-datenum(dataGa.Time(1))),'-k','LineWidth',2);
        end
        Tresult(i,:)={{station}, end_time, period(i), {granu}, {name} ,{inst},  {distribution}, {'LMS'},  {result.significance},{result.ss},{result.slope}, {result.UCL},{result.LCL},{result.slopeP}, {result.UCLP},{result.LCLP},{result.slopeR}, {result.UCLR},{result.LCLR}};
    end
end

% save the figure
if fig
    if ispc
        savefig(f,strcat('C:\github_trend\result\',(station),'\',(station),'_',(name),'_LMS_D.fig'));
    else
        savefig(f,strcat('/prod/pay/Aerosol_actris_trend/result/',(station),'/',(station),'_',(name),'_LMS_D.fig'));
    end
end

%_______________________________________________________________________

function fig_LMS(data,name,Eplot,x_residue,b)

subplot(2,2,1);
hold on
plot(data.Time,data.(name),'bo-','LineWidth',0.5);
plot(data.Time, Eplot,'-r','LineWidth',2);
plot(data.Time, b(1)+b(2).*(datenum(data.Time)-datenum(data.Time(1))),'-r','LineWidth',2);
ylabel('monthly data and fit');
datetick('x','yy','keeplimits');
grid on;
subplot(2,2,2);
plot(data.Time,x_residue);
ylabel('Residue');
datetick('x','yy','keeplimits');
grid on;
subplot(2,2,3);
normplot(real(x_residue));
title('normplot of residues');
subplot(2,2,4);
plot(data.Time,cumsum(x_residue));
datetick('x','yy','keeplimits');
ylabel('cumsum of residue');
grid on;


