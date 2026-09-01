function station_datD = read_betsy_2026(filename, STN)
% read the files send by betsy
% filename contains also the whole path
% STN= station name

%exemple:  SPO_rd=read_betsy_2026('/prod/pay/Aerosol_actris_trend/data/spo/spo_1979_2025','SPO');
%% read the data
if strcmp(STN,'SUM2')
    station_dat = importfile_SUM2("C:\github_trend\raw_data\sum_neph_psap_clap2_clap10_AE16_AE33", [2, Inf]);
    station_dat.Properties.DimensionNames{1}='Time';
elseif strcmp(STN,"SPL") || strcmp(STN, "UGR") || strcmp(STN, "APP")   || strcmp(STN, "CPR")
    station_dat = read_SPL_actris(filename, [2, Inf]);
    station_dat.Properties.DimensionNames{1}='Time';
elseif  strcmp(STN, "MLO")
    station_dat = importfile_MLO(filename, [2, Inf]);
elseif  strcmp(STN, "ATTO")
    station_dat = importfile_ATTO(filename, [2, Inf]);
    station_dat.Properties.DimensionNames{1}='Time';
elseif  strcmp(STN, "LLN")
    station_dat = importfile_LLN(filename, [2, Inf]);
    station_dat.Properties.DimensionNames{1}='Time';
elseif  strcmp(STN, "MBO")
    station_dat = importfile_MBO(filename, [2, Inf]);
else
    station_dat= read_spo_actris(filename);
end

%% probleme numeric negative: tentative to solve it
station_dat = timetable2num(station_dat);


%%
names_var=fieldnames(station_dat);
%% remove the variables containing g (STD) and N (nb of data)
c=contains(names_var,{'g','N'});
names_delete=names_var(c);
station_dat=removevars(station_dat,names_delete);
names_var=fieldnames(station_dat);
% Improve station, limit of scattering to 500 Mm-1 for the eastern stations (ACA, GSM, MCN, SHN)
% and 100 Mm-1 for western (BBE HCG IBB MRN MZW PAZ SCN SIA)
%check the RH for the scattering
west=['BBE','BLI','CRG','GWA','GLR','GBN', 'GRE','HGC', 'IBB', 'MRN', 'MZW', 'PAZ', 'RMN','SCN', 'SIA','TMO','UBW'];
east=['ACA', 'CRO','DSW','GSM','GGW','LBW', 'MCN','NCC', 'SHN'];
if contains(west,STN)
    CW=startsWith(names_var,["Bs","Bbs"]);
    N=names_var(CW);
    for i=1:length(N)
        ind=station_dat.(N{i})>100;
        station_dat.(N{i})(ind)=NaN;
    end
    Cc=startsWith(names_var,["T","U"]);
    N=names_var(Cc);
    for i=1:length(N)
        ind=station_dat.(N{i})>=9999 | station_dat.(N{i})>=99999 | station_dat.(N{i})>=99999;
        station_dat.(N{i})(ind)=NaN;
    end
elseif contains(east,STN)
    CE=startsWith(names_var,["Bs","Bbs"]);
    N=names_var(CE);
    for i=1:length(N)
        ind=station_dat.(N{i})>500;
        station_dat.(N{i})(ind)=NaN;
    end
    Cc=startsWith(names_var,["T","U"]);
    N=names_var(Cc);
    for i=1:length(N)
        ind=station_dat.(N{i})>=9999;
        station_dat.(N{i})(ind)=NaN;
    end
elseif strcmp(STN,'SUM2')

    Cc=startsWith(names_var,["B","U","X"]);
    N=names_var(Cc);
    for i=1:length(N)
        ind=station_dat.(N{i})>=9999;
        station_dat.(N{i})(ind)=NaN;
    end
else
    Cc=startsWith(names_var,["B","T","U"]) & ~startsWith(names_var,'Time');
    N=names_var(Cc);
    for i=1:length(N)
        ind=station_dat.(N{i})>=9999 ;
        station_dat.(N{i})(ind)=NaN;
    end
end

%specail reading
%erreur special pour AMY
if strcmp(STN,'AMY')==1
    Cu=startsWith(names_var,'U');
    Nu=names_var(Cu);
    for i=1:length(Nu)
        ind=station_dat.(Nu{i})>990;
        station_dat.(Nu{i})(ind)=NaN;
    end

    %erreur special pour PAL: false missing code for G and R Bbs
elseif strcmp(STN,'PAL')==1
    station_dat.BbsG_S11(station_dat.BbsG_S11==100)=NaN;
    station_dat.BbsR_S11(station_dat.BbsR_S11==100)=NaN;

    % time treatment for SPL
elseif strcmp(STN,'SPL')==1 || strcmp(STN,'UGR')==1  || strcmp(STN,'ATTO')==1  || strcmp(STN,'APP')==1 || strcmp(STN,'CPR')==1 

    station_dat(:,1)=[];
    station_dat.Properties.VariableNames{1}='y';
    station_dat.DOY=[];
    station_dat=table2timetable(station_dat);
    station_dat.Properties.DimensionNames{1}='Time';
    Cu=startsWith(names_var,["U";"B"]);
    Nu=names_var(Cu);
    for i=1:length(Nu)
        ind=station_dat.(Nu{i})>880 | station_dat.(Nu{i})<=-99;
        station_dat.(Nu{i})(ind)=NaN;
    end
elseif strcmp(STN,'LLN')==1

    station_dat(:,1)=[];
    station_dat.Properties.VariableNames{1}='y';
%    station_dat=table2timetable(station_dat);

elseif strcmp(STN, 'MBO')
    station_dat.MBO=[];
    station_dat.Properties.VariableNames{1}='y';
    station_dat.DOY=[];
    station_dat=table2timetable(station_dat);
    station_dat.Properties.DimensionNames{1}='Time';
    Cu=startsWith(names_var,["U";"B"]);
    Nu=names_var(Cu);
    for i=1:length(Nu)
        ind=station_dat.(Nu{i})>880 | station_dat.(Nu{i})<=-99;
        station_dat.(Nu{i})(ind)=NaN;
    end
elseif strcmp(STN,'MLO')==1
    station_dat.Properties.DimensionNames{1}='Time';
end

names_var=fieldnames(station_dat);
CU=contains(names_var,'U');
RH=names_var(CU);
for i=1:length(RH)
    N=RH{i}(2:end);
    %select RH too high
    if strcmp(STN,'GSN')==1 && strcmp(RH{i},'U_S11')
        RH{i}='U0_S11';

    end
    indRH=station_dat.(RH{i})>50 ;
    if strcmp(STN,'MCN')==1 | strcmp(STN,'MRN')==1
        indRH=station_dat.(RH{i})>50  | isnan(station_dat.(RH{i}));
    end
    % select variables associated with RH measurements
    C1=contains(names_var,RH{i}(2:end));
    CB=startsWith(names_var,'B');
    Cok=C1 & CB;
    S1=names_var(Cok);
    %special case of fragmentation of GSN scat but not RH

    for j=1:length(S1)
        S=strcat((S1{j}),'_dry');
        station_dat.(S)=station_dat.(S1{j});
        station_dat.(S)(indRH)=NaN;
    end
end


%test RMN without negatives
% % if strcmp(STN,'RMN')
% %     SR=names_var(C1);
% %     for i=1:length(SR)
% %         indRMN=station_dat.(SR{i})<0;
% %         station_dat.(SR{i})(indRMN)=NaN;
% %     end
% % end
% daily mean, apply a data coverage threshold ?
station_datD=retime(station_dat,'daily',@nanmedian);
station_datN=retime(station_dat,'daily','count');
names_var=fieldnames(station_dat);
CB=startsWith(names_var,'B');
SB=names_var(CB);
for i=1:length(SB)
    indN=station_datN.(SB{i})<6 ;
    a(i)=sum(indN)*100/sum(~isnan(station_datD.(SB{i})));
    station_datD.(SB{i})(indN)=NaN;
end
ind=a>20;
if sum(ind)>0
    warning('there is a lot of days with less than 50% data coverage');

end


end

