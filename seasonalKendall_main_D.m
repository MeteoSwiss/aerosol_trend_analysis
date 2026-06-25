function Tresult=seasonalKendall_main_D(data,param, inst, station, resolution, varargin)

% compute the s.s of the selected data with PW and TFPW method
%compute the slope with VCTFPW method
% if not specified, end_year==2025 and the trends are computed for 10, 20 and 30 years
% if not specified, the granularity is one day with median as averaging


%IN: data:data either as a structure with time in data.start_time
%                     or a timetable
%                     or a matrice with the time as the first 6 columns (datevec)
% param = cell: either sturcture or timetable name or number of matrix column
% instrument= specify the instrument
% resolution= resolution is taken into account to determine the number of ties.
% varargin:
% granularity= {'daily','monthly','qualerly','yearly'} or length in [day]
%               default= 'daily'
% averageM= average method {'mean',@median}
%       default=@median
% period= period for trend analysis, default: 10, 20, 30 and 40 years before end_time
% end_year, default: 2025
%fig= if 1, figure are plotted and saved

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
%               results = structure of results with
%result.period:year or 12 monthes + year or 4 meteorological seasons+year
%result.ss= statistical significance: 95=95% for both PW and TFPW tests
%                                     90=90% for both PW and TFPW tests
%                                     10= test is s.s at 95% for TFPW and not s.s. for PW (false positive)
%                                     20= test is s.s. at 95% for PW and not s.s. for TFPW (flase negative)
%                                     0= not s.s. for both tests
%result.slope: Sen's slope in units/y
%result.UCL: upper confidence level in units/y
%result.LCL: lower confidence level in units/y
%result.median: median of the data
%result.slopeP: Sen's slope in %/y
%result.UCLP: upper confidence level in %/y
%result.LCLP: lower confidence level in %/y
%result.ak: first lag autocorrelation after detrending
%result.Xhomo: homogeneity of the trend per seasons/month: 1=homogeneous at 90% confidence level
%                                                          0= not homogeneous: the yearly slope should not be computeated
%                                                          NaN= not assigned for yearly trend

%plot of the monthly and metSeason results

%example: Tresult=seasonalKendall_main_D(test2, {'conc'},'cpc', 'JFJ', 1,'end_year',2018,'period',10,'fig',1);

%Martine Collaud Coen, March 2019

% % if nargin < 4
% %     error('Not enough arguments, you have to specify the name or the column of the parameter')
% % end
if ischar(param{1})
    name=param{1};
elseif isnumeric(param{1})
    colonne=param{1};
end
% check arguments
if ~varg_proof(varargin, {'granularity', 'aveM', 'end_year','path','period','fig'},true)
    return
end

% Set values from user input, or use defaults
granu = varg_val(varargin, 'granularity', 'daily');
avM = varg_val(varargin, 'aveM', @nanmedian);
end_time = varg_val(varargin, 'end_year', 2025);
path = varg_val(varargin, 'path', '');
periodx = varg_val(varargin, 'period', [10 20 30 40 50]);
fig=varg_val(varargin, 'fig', 0);
%put the data in timetable
if isstruct(data)
    dataT=timetable(datetime(datevec(data.start_time)),data.(name));
    % %     data=rmfield(data,'start_time');
    % %      N=fieldnames(data);
    % %     for i=1:length(N)
    % %     dataT.(N{i})=data.(N{i});
    % %     end
elseif istimetable (data)
    dataT=timetable(data.Time,data.(name),'VariableNames',{name});
elseif data(:,1)<2100 % data over 6 columns
    dataT=timetable(datetime(data(:,1:6)), data(:,colonne),'VariableNames',{'parametre'});
    name=dataT.Properties.VariableNames{1};
else % time in datenum
    dataT=timetable(datetime(datevec(data(:,1))), data(:,colonne),'VariableNames',{'parametre'});
    name=dataT.Properties.VariableNames{1};
end
if granu>29 & granu<31
    dataG=retime(dataT,'monthly',avM);
else
    dataG=retime(dataT,granu,avM);
end
dataG.y=year(dataG.Time);
%take the number of periods to be analysed
if length(periodx)>1
    if periodx(end)-periodx(end-1)<10
        nb_period=floor((max(dataG.y)-min(dataG.y)+1)/10)-1;
        if nb_period==0
            nb_period=1;
        end
    else
        nb_period=floor((max(dataG.y)-min(dataG.y)+1)/10);
    end
else
    nb_period=1;
end


if nb_period==1
    delta=0.4;
elseif nb_period==2
    delta=0.3;
else
    delta=0.8/(nb_period-1);
end
% for each period

% % Tresult=cell(3*nb_period,9); Tresult=array2table(Tresult);
% % Tresult.Properties.VariableNames={'station','end_year','period','granularity','parameter','inst','season','method','result'};
%marche pas Tresult.Properties.VariableTypes={'string','double','double','string','string','string','string','string','string'};
for i=1:nb_period
    TT=timerange(datetime(end_time-periodx(i)+1,1,1,0,0,0),datetime(end_time+1,1,1,0,0,0));
    dataGTT=dataG(TT,:);
    result_y=struct2table(MK_tempAggr(dataGTT(:,1),resolution));
    dataGTT.m=month(dataGTT.Time);
    % result_y.slopeP=result_y.slope.*100./abs(nanmedian(dataGTT(:,1)));
    % result_y.UCLP=result_y.UCL.*100./abs(nanmedian(dataGTT(:,1)));
    % result_y.LCLP=result_y.LCL.*100./abs(nanmedian(dataGTT(:,1)));
    %MK for 4 meteorological seasons
    data_seas=struct;
    for m=1:4
        if m==4
            data_seas(m).obs=dataGTT(dataGTT.m==1 | dataGTT.m==2 |dataGTT.m==12,1);
        else
            data_seas(m).obs=dataGTT(dataGTT.m>m*3-1 & dataGTT.m<m*3+3,1);
        end
    end

    result_MS=struct2table(MK_tempAggr(data_seas, resolution));
    data_m=struct;
    %MK for 12 months
    for m=1:12
        data_m(m).obs=dataGTT(dataGTT.m==m,1);
    end

    result_mo=struct2table(MK_tempAggr(data_m, resolution));

    t_mo=[1:1:13]+(i-1)*delta.*ones(1,13);
    t_MS=[1:1:5]+(i-1)*delta.*ones(1,5);

    %figure with slope in units/year
    if fig
        f=figure(100); %in units/year
        subplot(2,1,1);
        hold on;
        if isstr(param{1})
            title(join([(station) ,(param),{inst},strcat(num2str(end_time-periodx(i)+1),'-',num2str(end_time))]));
        else
            title(join([(station) ,num2str(param{1}),{inst},strcat(num2str(end_time-periodx(i)+1),'-',num2str(end_time))]));
        end
        fig_seasonKendall_trend(t_mo,result_mo,result_y,'units',delta);
        subplot(2,1,2);
        hold on;
        fig_seasonKendall_trend(t_MS,result_MS,result_y,'units',delta);
        text(t_MS(1),result_MS.UCL(1)+abs(result_MS.UCL(1))*30/100,join([num2str(periodx (i)),'y']),'FontSize',12,'HorizontalAlignment','center');

        %figure with slope in %/year
        % % f(2)=figure(101); %in %/year
        % % subplot(2,1,1);
        % % hold on;
        % % if isstr(param{1})
        % %     title(join([(station) ,(param),{inst},strcat(num2str(end_time-periodx(i)+1),'-',num2str(end_time))]));
        % % else
        % %     title(join([(station) ,num2str(param{1}),{inst},strcat(num2str(end_time-periodx(i)+1),'-',num2str(end_time))]));
        % % end
        %title(join([(station) ,(param),{inst},strcat(num2str(end_time-periodx(i)+1),'-',num2str(end_time))]));
        % % fig_seasonKendall_trend(t_mo,result_mo,result_y,'%',delta);
        % % subplot(2,1,2);
        % % hold on;
        % % fig_seasonKendall_trend(t_MS,result_MS,result_y,'%',delta);
        % % text(t_MS(1),result_MS.UCLP(1)+abs(result_MS.UCLP(1))*30/100,join([num2str(periodx (i)),'y']),'FontSize',12,'HorizontalAlignment','center');
    end
    %create the table of results
    if i==1
        % % Tresult=table({station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'y'}, {'MK'}, {result_y.ss},{result_y.slope},{result_y.UCL},{result_y.LCL},{result_y.slopeP},{result_y.UCLP},{result_y.LCLP},{result_y.Xhomo},...
        % %     'VariableNames',{'station','end_time','length_period','granularity','parameter','instrument','MK_seasonality','method','ss','slope','UCL','LCL','slopeP','UCLP','LCLP','Xhomo'});
        % % Tresult(2,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'MetSea'}, {'MK'}, {result_MS.ss},{result_MS.slope},{result_MS.UCL},{result_MS.LCL},{result_MS.slopeP},{result_MS.UCLP},{result_MS.LCLP},{result_MS.Xhomo}};
        % % Tresult(3,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'month'}, {'MK'}, {result_mo.ss},{result_mo.slope},{result_mo.UCL},{result_mo.LCL},{result_mo.slopeP},{result_mo.UCLP},{result_mo.LCLP},{result_mo.Xhomo}};
    Tresult=table({station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'y'}, {'MK'}, {result_y.ss},{result_y.slope},{result_y.UCL},{result_y.LCL},...
            'VariableNames',{'station','end_time','length_period','granularity','parameter','instrument','MK_seasonality','method','ss','slope','UCL','LCL'});
        Tresult(2,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'MetSea'}, {'MK'}, {result_MS.ss},{result_MS.slope},{result_MS.UCL},{result_MS.LCL}};
        Tresult(3,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'month'}, {'MK'}, {result_mo.ss},{result_mo.slope},{result_mo.UCL},{result_mo.LCL}};

    else
        % % Tresult((i-1)*3+1,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'y'}, {'MK'}, {result_y.slope},{result_y.UCL},{result_y.LCL},{result_y.slopeP},{result_y.UCLP},{result_y.LCLP},{result_y.Xhomo}};
        % % Tresult((i-1)*3+2,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'MetSea'}, {'MK'}, {result_MS.UCL},{result_MS.LCL},{result_MS.slopeP},{result_MS.UCLP},{result_MS.LCLP},{result_MS.Xhomo}};
        % % Tresult((i-1)*3+3,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'month'}, {'MK'}, {result_mo.ss},{result_mo.slope},{result_mo.UCL},{result_mo.LCL},{result_mo.slopeP},{result_mo.UCLP},{result_mo.LCLP},{result_mo.Xhomo}};
    Tresult((i-1)*3+1,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'y'}, {'MK'},{result_y.ss} {result_y.slope},{result_y.UCL},{result_y.LCL}};
        Tresult((i-1)*3+2,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'MetSea'}, {'MK'}, {result_MS.ss},{result_MS.slope},{result_MS.UCL},{result_MS.LCL}};
        Tresult((i-1)*3+3,:)={{station}, end_time, periodx(i), {granu}, {name} ,{inst},  {'month'}, {'MK'}, {result_mo.ss},{result_mo.slope},{result_mo.UCL},{result_mo.LCL}};

    end
end
if fig %&& exist f

    if ispc
        savefig(f,cell2mat(strcat('C:\github_trend\result\',(station),'\',(station),'_',(param),'_MK.fig')));
    else
        savefig(f,cell2mat(strcat('/prod/pay/Aerosol_actris_trend/result/',(station),'/',(station),'_',(param),'_MK.fig')));
    end
end
%save(cell2mat(strcat(path,(station),'_',(param),'_MK.mat')),'Tresult');
close all;
%_____________________________________________________
function fig_seasonKendall_trend(t,result,result_y,type,delta)
sizeM=[15 20 30];
if strcmp( type,'units')==1
    s='slope';
    U='UCL';
    L='LCL';
    ylabel('Slope [units/y]');
elseif strcmp( type,'%')==1
    s='slopeP';
    U='UCLP';
    L='LCLP';
    ylabel('Slope [%/y]','FontSize',14);
end
plot(t,result.(s),'b.','MarkerSize',sizeM(1));
hold on;
plot(t(result.ss==90),result.(s)(result.ss==90),'b.','MarkerSize',sizeM(2));
plot(t(result.ss==95),result.(s)(result.ss==95),'b.','MarkerSize',sizeM(3));
plot(t(result.ss==10),result.(s)(result.ss==10),'c^','MarkerSize',sizeM(1));
plot(t(result.ss==20),result.(s)(result.ss==20),'cv','MarkerSize',sizeM(1));
line(repmat(t,2,1),[result.(U),result.(L)]','color','b');
%plot(t,result.(s),'rs','MarkerSize',sizeM(3));
if result_y.ss==95
    a=sizeM(3);
elseif result_y.ss==90
    a=sizeM(2);
else
    a=sizeM(1);
end
plot(t(end)+delta/2,result_y.(s),'g.','MarkerSize',a);
line([t(end)+delta/2 t(end)+delta/2 ],[result_y.(U);result_y.(L)],'color','g');

if length(t)==13
    xlabel('Month','FontSize',14);
    xticks([1:1:13]);
    xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','year'});
    xlim([0.8 14]);
elseif length(t)==5
    xlabel('Season','FontSize',14);
    xticks([1:1:5]);
    xticklabels({'Spring','Summer','Fall','Winter','year'});
    xlim([0.8 6]);
end
grid on;
line(xlim,[0 0],'color','r','LineWidth',1);
ax=gca;
ax.FontSize=14;







