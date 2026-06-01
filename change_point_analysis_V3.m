function Tresult= change_point_analysis_V3(data, param, alpha, station, inst )

% apply the Pettitt and monthly data
% the tests are applied only for at least 10 years on yearly data -->
% subsequent division of the time serie also follows this rule
% tests on monthly data should not be applied on less than 2 years if divided

% find a way to extract the best potentiel change points (and not all)
% save fig all in one


% inputs: data = TimeTable
%           param = parameters (as fieldnames) to be analysed
%           alpha = confidence level, set at 0.05 as default

% outputs:  result= table of results


% use: ...  ... ...

% example:


% fev-mars 2026, mco
%____________________________________________________
% remove lines with NaN
for i=1:length(param)
    ind_nan(:,i)=isnan(data.(param{i}));
end
data(sum(ind_nan,2)==length(param),:)=[];

% produce yearly and monthly medians
data_m=retime(data,'monthly',@nanmedian);
data_m.month=month(data_m.Time);

start_time=min(year(data_m.Time));
end_time=max(year(data_m.Time));
period=max(year(data_m.Time))-min(year(data_m.Time))+1;

%deseasonlise with the median of the month
month_med=groupsummary(data_m,'month',@nanmedian);
data_m_deseason=timetable(data_m.Time);
data_m_deseason.month=data_m.month;
for i=1:length(param)
    data_m_deseason.(param{i})=data_m.(param{i})./month_med.(strcat('fun1_',param{i}))(data_m_deseason.month);
end
%N_m=length(data_m.Time);

%deseasonalise and detrend by LMS on data
data_m_residue=timetable(data_m.Time);
%deseasonalise and detrend by LMS on log(data)
data_m_residueLog=timetable(data_m.Time);
% remove nan from monhtly time series since it is the lowest used time granularity

% Prepare figure
% nb of subplot for the figure, i.e. two lines (param + ratio) if more than one param
if length(param)==1
    nb_subplot =3;
else
    nb_subplot=4; % plot of the ratio are needed for the visual inspection
end
% color order
couleur=['b' 'c' 'g' 'r' 'm' 'y' 'k' 'k'];
grandeur=[12 11 10 9 8 7 6 5 4];

figure;
hold on;

for i=1:length(param)
    % % forward direction, multiple break points
    % compute the residue of LMS fit with data and log(data)
    data_m_residue.(param{i})=LMS_residue(data_m, param{i}, 'lin');
    data_m_residueLog.(param{i})=LMS_residue(data_m, param{i}, 'log');
    % search breakpoints for monthly data with season
    granu='month';
   
    % % result_m = multiple_breakpoints_while_one(data_m,{param{i}}, 12, 0.05);
    % % Tresult(i,:)=table({station},start_time, end_time,period, {granu}, {param{i}} ,{inst},  {'Pettitt'}, {result_m},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});
    % search breakpoints for monthly data deseasonalised with median for each month
    result_m = multiple_breakpoints_while_one(data_m_deseason,{param{i}}, 12, 0.05);
    Tresult(0*length(param)+i,:)=table({station},start_time, end_time,period, {granu}, {strcat(param{i},'_deseason')} ,{inst},  {'Pettitt'}, {result_m},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});
    % % % search breakpoints for monthly data detrend and deseasolized with LMS
    result_m = multiple_breakpoints_while_one(data_m_residue,{param{i}}, 12, 0.05);
    Tresult(1*length(param)+i,:)=table({station},start_time, end_time,period, {granu}, {strcat(param{i},'_residue')} ,{inst},  {'Pettitt'}, {result_m},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});
    % % % search breakpoints for log of monthly data detrend and deseasolized with LMS
    result_m = multiple_breakpoints_while_one(data_m_residueLog,{param{i}}, 12, 0.05);
    Tresult(2*length(param)+i,:)=table({station},start_time, end_time,period, {granu}, {strcat(param{i},'_residueLog')} ,{inst},  {'Pettitt'}, {result_m},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});

   

% % backward direction (inv), multiple break points
 % produce the backward data
        data_m.(strcat(param{i},'_inv'))=flipud(data_m.(param{i}));
        %data_m_deseason.(param{i})=flipud(data_m_deseason.(param{i}));
        %data_m_residue.(param{i})=flipud(data_m_residue.(param{i}));
        data_m_residueLog.(strcat(param{i},'_inv'))=flipud(data_m_residueLog.(param{i}));
 result_m = multiple_breakpoints_while_one(data_m_residueLog,{strcat(param{i},'_inv')}, 12, 0.05);

% adapt the time of the results
       % for l=1:height(result_m.time)
            if ~isnat( result_m.time)
                for j=1:height(result_m.time)
                    pos=find(data_m.Time==result_m.time(j));
                    N=height(data_m.Time);
                    result_m.time(j)= data_m.Time(N-pos) ;
                end
            end
       % end
    Tresult(3*length(param)+i,:)=table({station},start_time, end_time,period, {granu}, {strcat(param{i},'_residueLog','_inv')} ,{inst},  {'Pettitt'}, {result_m},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});

% % Pettitt from all 2 y periods in forward direction only

 % % % result_m =multiple_breakpoints_2y_period(data_m_residueLog,{param{i}}, 12, 0.05);
 % % %    Tresult(4*length(param)+i,:)=table({station},start_time, end_time,period, {granu}, {strcat(param{i},'_residueLog','_2yper')} ,{inst},  {'Pettitt'}, {result_m},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});

    % % plot of data with breakpoints and cusum
    subplot(nb_subplot,1,1);
    title(station);
    yyaxis left
    hold on;
    plot(data_m.Time, data_m.(param{i}),'.','Color',couleur( i))
    ylabel(inst);
    yyaxis right
    hold on;
    % selectx= Tresult.results{strcmp(Tresult.parameter,param{i})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    % p1=plot(selectx.time+days(30/2),selectx.pvalue,'d','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur(i),'LineWidth',1.5);
    selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_deseason')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    p2=plot(selectx.time+days(30/2),selectx.pvalue_boot,'o','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5); %'MarkerFaceColor',couleur(i),
    selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residue')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    p3=plot(selectx.time+days(30/2),selectx.pvalue_boot,'s','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);
    selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residueLog')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    p4=plot(selectx.time+days(30/2),selectx.pvalue_boot,'^','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);
    selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residueLog','_inv')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    p5=plot(selectx.time+days(30/2),selectx.pvalue_boot,'v','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);
%%  selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residueLog','_2yper')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
   %% p6=plot(selectx.time+days(30/2),selectx.pvalue_boot,'*','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);
  
    %legend([p1 p2 p3 p4],{'break data','break deseason','break LMS','break LMS of log'});
    legend([p2(1) p3(1) p4(1) p5(1)],{'deseason', 'residue','resLog forward', 'resLog backward'},'Location','northeastoutside');
%%        legend([p2(1) p3(1) p4(1) p5(1) p6(1)],{'deseason', 'residue','resLog forward', 'resLog backward','resLog 2yPer'},'Location','northeastoutside');

    %ylim([0, alpha]);
    ylabel('p-value bootstrap');
    grid on;
% plot of break points (only without log  with deseason/detrend data
 subplot(nb_subplot,1,2);
 yyaxis left
    hold on;
    plot(data_m_deseason.Time, data_m_deseason.(param{i}),'.','Color',couleur( i));
     plot(data_m_residue.Time, data_m_residue.(param{i}),'.','Color',couleur( i+1));
    ylabel(inst);
    yyaxis right
    hold on;
    % selectx= Tresult.results{strcmp(Tresult.parameter,param{i})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    % p1=plot(selectx.time+days(30/2),selectx.pvalue,'d','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur(i),'LineWidth',1.5);
    selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_deseason')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    p2=plot(selectx.time+days(30/2),selectx.PrctDiff(:,2),'o','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i)); %'MarkerFaceColor',couleur(i),
    selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residue')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    p3=plot(selectx.time+days(30/2),selectx.PrctDiff(:,2),'s','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i));
  % %   selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residueLog')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
  % %   p4=plot(selectx.time+days(30/2),selectx.PrctDiff(:,2),'^','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);
  % %   selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residueLog','_inv')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
  % %   p5=plot(selectx.time+days(30/2),selectx.PrctDiff(:,2),'v','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);
  % % %% selectx= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residueLog','_2yper')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
   %% p6=plot(selectx.time+days(30/2),selectx.PrctDiff(:,2),'*','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);
  
    %legend([p1 p2 p3 p4],{'break data','break deseason','break LMS','break LMS of log'});
    legend([p2(1) p3(1)],{'deseason', 'residue'},'Location','northeastoutside');
    %    legend([p2(1) p3(1) p4(1) p5(1) p6(1)],{'deseason', 'residue','resLog forward', 'resLog backward','resLog 2yPer'},'Location','northeastoutside');

    %ylim([0, alpha]);
    ylabel('Median Diff');
    grid on;

    %plot(selectx.time_inv+days(30/2),selectx.pvalue_inv,'sm','MarkerSize',10, 'MarkerFaceColor','m');
    subplot(nb_subplot,1,3);
    hold on;
    %plot(data_m_deseason.Time,cumsum(data_m_deseason.(param{i}),"omitnan"),'-','Color',couleur(i));
    %plot(data_m_residue.Time,cumsum(data_m_residue.(param{i}),"omitnan"),'--','Color',couleur(i));
    plot(data_m_residueLog.Time,cumsum(data_m_residueLog.(param{i}),"omitnan"),':','Color',couleur(i));
  legend('residue LMS of log','Location','northeastoutside'); % legend('deseason','residue LMS','residue LMS of log');
    ylabel(strcat('CumSum ',param{1}));
    grid on;
end

%% if nb_param >1 do all the same for the ratio of the param
if length(param)>1
    for i=1:length(param)-1
        for j=i+1:length(param)
            data_m.(strcat(param{i},'_',param{j}))=data_m.(param{i})./data_m.(param{j});
            result_m = multiple_breakpoints_while_one(data_m,{strcat(param{i},'_',param{j})}, 12, 0.05);
            Tresult((2+i-1)*length(param)+(j-i),:)=table({station},start_time, end_time,period, {granu}, {(strcat(param{i},'_',param{j}))} ,{inst},  {'Pettitt'}, {result_m},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});

            % % if strcmp(direction,'inv')
            % %     %for l=1:height( Tresult)
            % %         if ~isempty( Tresult.station{(2+i-1)*length(param)+(j-i)}) && startsWith(Tresult.parameter{(2+i-1)*length(param)+(j-i)},(strcat(param{i},'_',param{j})))
            % %             for k=1:height(Tresult.results{(2+i-1)*length(param)+(j-i)}.time)
            % %                 pos=find(data_m.Time==Tresult.results{(2+i-1)*length(param)+(j-i)}.time(k));
            % %                 N=height(data_m.Time);
            % %                 if ~isempty(pos)
            % %                 Tresult.results{(2+i-1)*length(param)+(j-i)}.time(k)= data_m.Time(N-pos) ;
            % %                 end
            % %             end
            % %         end
            % %    % end
            % % end
            % subplot(nb_subplot,1,1);
            % hold on;
            % yyaxis left;
            % hold on;
            % %plot(data_m.Time,data_m.(strcat(param{i},'_',param{j})),'-');
            % %ylabel(param{1});
            % yyaxis right
            % p9=plot(result_m.time+days(30/2),result_m.pvalue,'vk','MarkerSize',12,'LineWidth',1);
            % ylim([0, alpha]);
            % ylabel('p-value');
            % legend([p1 p4 p9],{'break data','break LMS of log' 'all ratios'},'Location','northeastoutside')
       subplot(nb_subplot,1,4);
            hold on;
            yyaxis left;
            hold on;
            p10=plot(data_m.Time,data_m.(strcat(param{i},'_',param{j})),'-','Color',couleur(i));
            ylabel(inst);
            yyaxis right
            p9=plot(result_m.time+days(30/2),result_m.pvalue_boot,'vk','MarkerSize',12,'LineWidth',1,'Color',couleur(i));
            ylim([0, alpha]);
            ylabel('p-value');
            legend([p10(1) p9(1)],{'ratio', 'all ratios'},'Location','northeastoutside')
        end
    end

end