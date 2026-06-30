function Tresult= change_point_analysis_Def(data, param, alpha, station, inst )

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

%deseasonalise and detrend by LMS on data
data_m_residue=timetable(data_m.Time);

% remove nan from monhtly time series since it is the lowest used time granularity

% Prepare figure
% nb of subplot for the figure
nb_subplot =3;

% color order
couleur=['b' 'c' 'g' 'r' 'm' 'y' 'k' 'k' 'b' 'c' 'g' 'r' 'm' 'y' 'k'];
grandeur=[12 11 10 9 8 7 6 5 4];

figure;
hold on;

for i=1:length(param)
    % % forward direction, multiple break points
    % compute the residue of LMS fit with data
    data_m_residue.(param{i})=LMS_residue(data_m, param{i}, 'lin');

    % search breakpoints for monthly data with season
    granu='month';

    % search Pettitt breakpoints for monthly data detrend and deseasolized with LMS
    result_mPT = multiple_breakpoints_while_one(data_m_residue,{param{i}}, 12, alpha);
    Tresult(0*length(param)+i,:)=table({station},start_time, end_time,period, {granu}, {strcat(param{i},'_residue')} ,{inst},  {'Pettitt'}, {result_mPT},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});

    % search SNHT breakpoints for monthly data detrend and deseasolized with LMS
    result_mSNHT = multiple_breakpoints_while_one_snht(data_m_residue,{param{i}}, 12, alpha);
    Tresult(1*length(param)+i,:)=table({station},start_time, end_time,period, {granu}, {strcat(param{i},'_residue')} ,{inst},  {'SNHT'}, {result_mSNHT},'VariableNames',{'station','start_time','end_time','length_period','granularity','parameter','instrument','method','results'});


    % % plot of data with breakpoints and cusum
    subplot(nb_subplot,1,1);
    title(station);
    yyaxis left
    hold on;
    plot(data_m.Time, data_m.(param{i}),'.','Color',couleur( i))
    ylabel(inst);
    yyaxis right
    hold on;
    selectPT= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residue')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'Pettitt')==1 ,:};
    % I don't know why, but too low values are not plotted
    selectPT.pvalue(selectPT.pvalue<0.001)=0.001;
    p3=plot(selectPT.time+days(30/2),selectPT.pvalue,'o','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);

    selectSNHT= Tresult.results{strcmp(Tresult.parameter,{strcat(param{i},'_residue')})==1 & strcmp(Tresult.granularity,'month')==1 & strcmp(Tresult.method,'SNHT')==1 ,:};
   % I don't know why, but too low values are not plotted
    selectSNHT.pvalue(selectSNHT.pvalue<0.001)=0.001;
    p4=plot(selectSNHT.time+days(30/2),selectSNHT.pvalue,'s','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i), 'LineWidth',1.5);


    %legend([p1 p2 p3 p4],{'break data','break deseason','break LMS','break LMS of log'});
    legend([p3(1) p4(1)],{'Pettitt', 'SNHT'},'Location','northeastoutside');
    %%        legend([p2(1) p3(1) p4(1) p5(1) p6(1)],{'deseason', 'residue','resLog forward', 'resLog backward','resLog 2yPer'},'Location','northeastoutside');

    %ylim([0, alpha]);
    ylabel('p-value bootstrap');
    grid on;
    % plot of break points (only without log  with deseason/detrend data
    subplot(nb_subplot,1,2);
    title('Residues from deseason + detrend');
    yyaxis left
    hold on;
    plot(data_m_residue.Time, data_m_residue.(param{i}),'.','Color',couleur( i+1));
    ylabel(inst);
    legend(param,'Location','northeastoutside');
    yyaxis right
    hold on;
    plot(selectPT.time+days(30/2),selectPT.PrctDiff(:,2),'o','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i),'DisplayName', 'Pettitt'); %'MarkerFaceColor',couleur(i),
    plot(selectSNHT.time+days(30/2),selectSNHT.PrctDiff(:,2),'s','MarkerSize',grandeur(i), 'MarkerEdgeColor',couleur( i),'DisplayName', 'SNHT');

    %ylim([0, alpha]);
    ylabel('Median Diff');
    grid on;

    %plot(selectx.time_inv+days(30/2),selectx.pvalue_inv,'sm','MarkerSize',10, 'MarkerFaceColor','m');
    subplot(nb_subplot,1,3);
    title('Cumulative sum' );
    hold on;
    plot(data_m_residue.Time,cumsum(data_m_residue.(param{i}),"omitnan"),'--','Color',couleur(i),'DisplayName','residue');
    legend('Location','northeastoutside'); % legend('deseason','residue LMS','residue LMS of log');
    ylabel(strcat('CumSum ',param{1}));
    grid on;
end


end