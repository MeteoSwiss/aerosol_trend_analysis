function [STN_rd, STN_st]=read_ebas_merged_STN(filename)
% read the Netcdf file and adapt it to the naming convention for the trend
% analysis 2026.

% INPUT
% filename: path and name of the Netcdf file

% OUTPUTfunction [STN_rd, STN_st]=read_ebas_merged_STN(filename)
% read the Netcdf file and adapt it to the naming convention for the trend
% analysis 2026.

% INPUT
% filename: path and name of the Netcdf file

% OUTPUT
%   STN_rd: data of the station in daily time granularity
%   STN_st: information on the station (accronym, lat, lon, alt, env) if
% available in the original file. ATTENTION: I don't know if these
% info are right, please check thereafter.

% DEPENDENCE: Timetable_naming_trend.m

% Exemple: [SMR_rd, SMR_st]=EbasNC2Timetable("C:\Users\mco\Downloads\merged_Hyytiälä.nc") ;

% take info from Netcdf
info=ncinfo(filename);

% Create a structure to store imported netCDF data
merged_data = struct();

% read all the attributes
for i=1:length(info.Attributes)

    merged_data.Attributes(i).Name = info.Attributes(i).Name;
    merged_data.Attributes(i).Value = ncreadatt(filename, "/", info.Attributes(i).Name);
end

%read all the variables
for i=1:length(info.Variables)
    merged_data.Variables(i).Name = info.Variables(i).Name;
    merged_data.Variables(i).Value = ncread(filename, info.Variables(i).Name);
end
%% take information on the station
STN_st=struct();

% accronym
ind_gaw_id=contains({merged_data.Attributes.Name},'ebas_station_gaw_id');
if sum(ind_gaw_id)>0
STN_st.name=merged_data.Attributes(ind_gaw_id).Value;
end

%latitute, longitude and altitude
ind_gaw_id=contains({merged_data.Attributes.Name},'latitude');
if sum(ind_gaw_id)>0
STN_st.lat=merged_data.Attributes(ind_gaw_id).Value;
end

ind_gaw_id=contains({merged_data.Attributes.Name},'longitude');
if sum(ind_gaw_id)>0
STN_st.lon=merged_data.Attributes(ind_gaw_id).Value;
end

ind_gaw_id=contains({merged_data.Attributes.Name},'altitude');
if sum(ind_gaw_id)>0
STN_st.alt=merged_data.Attributes(ind_gaw_id).Value;
end
% environment and footprint
ind_gaw_id=contains({merged_data.Attributes.Name},'ebas_station_setting');
if sum(ind_gaw_id)>0
STN_st.env=merged_data.Attributes(ind_gaw_id).Value;
end

%% Transform the structure in Timetable (necessary for matlab)
% Modify the parameter's name according to the naming convention
STN_rd_h=Timetable_naming_trend(merged_data);
names_var=fieldnames(STN_rd_h);
% % Improve station, limit of scattering to 500 Mm-1 for the eastern stations (ACA, GSM, MCN, SHN)
% and 100 Mm-1 for western (BBE HCG IBB MRN MZW PAZ SCN SIA)
west=['BBE','BLI','CRG','GWA','GLR','GBN', 'GRE','HGC', 'IBB', 'MRN', 'MZW', 'PAZ', 'RMN','SCN', 'SIA','TMO','UBW'];
east=['ACA', 'CRO','DSW','GSM','GGW','LBW', 'MCN','NCC', 'SHN'];
if contains(west,STN_st.name)
    CW=startsWith(names_var,["Bs","Bbs"]);
    N=names_var(CW);
    for i=1:length(N)
        ind=STN_rd_h.(N{i})>100;
        STN_rd_h.(N{i})(ind)=NaN;
    end

elseif contains(east,STN_st.name)
    CE=startsWith(names_var,["Bs","Bbs"]);
    N=names_var(CE);
    for i=1:length(N)
        ind=STN_rd_h.(N{i})>500;
        STN_rd_h.(N{i})(ind)=NaN;
    end
end

% compute the daily data
STN_rd=retime(STN_rd_h,'daily', @nanmedian);
%  25% minimal data coverage requirement
STN_datN=retime(STN_rd_h,'daily','count');
names_var=fieldnames(STN_rd_h);
CB=startsWith(names_var,'B');
SB=names_var(CB);
for i=1:length(SB)
    indN=STN_datN.(SB{i})<6 ;
    a(i)=sum(indN)*100/sum(~isnan(STN_rd.(SB{i})));
    STN_rd.(SB{i})(indN)=NaN;
end
ind=a>20;
if sum(ind)>0
    warning('there is a lot of days with less than 50% data coverage');
end



%   STN_rd: data of the station in daily time granularity
%   STN_st: information on the station (accronym, lat, lon, alt, env) if
% available in the original file. ATTENTION: I don't know if these
% info are right, please check thereafter.

% DEPENDENCE: Timetable_naming_trend.m

% Exemple: [SMR_rd, SMR_st]=EbasNC2Timetable("C:\Users\mco\Downloads\merged_Hyytiälä.nc") ;

% take info from Netcdf
info=ncinfo(filename);

% Create a structure to store imported netCDF data
merged_data = struct();

% read all the attributes
for i=1:length(info.Attributes)

    merged_data.Attributes(i).Name = info.Attributes(i).Name;
    merged_data.Attributes(i).Value = ncreadatt(filename, "/", info.Attributes(i).Name);
end

%read all the variables
for i=1:length(info.Variables)
    merged_data.Variables(i).Name = info.Variables(i).Name;
    merged_data.Variables(i).Value = ncread(filename, info.Variables(i).Name);
end

% Transform the structure in Timetable (necessary for matlab)
% Modify the parameter's name according to the naming convention
STN_rd=Timetable_naming_trend(merged_data);

% compute the daily data
STN_rd=retime(STN_rd,'daily', @nanmedian);

%% take information on the station
STN_st=struct();

% accronym
ind_gaw_id=contains({merged_data.Attributes.Name},'ebas_station_gaw_id');
if sum(ind_gaw_id)>0
STN_st.name=merged_data.Attributes(ind_gaw_id).Value;
end

%latitute, longitude and altitude
ind_gaw_id=contains({merged_data.Attributes.Name},'latitude');
if sum(ind_gaw_id)>0
STN_st.lat=merged_data.Attributes(ind_gaw_id).Value;
end

ind_gaw_id=contains({merged_data.Attributes.Name},'longitude');
if sum(ind_gaw_id)>0
STN_st.lon=merged_data.Attributes(ind_gaw_id).Value;
end

ind_gaw_id=contains({merged_data.Attributes.Name},'altitude');
if sum(ind_gaw_id)>0
STN_st.alt=merged_data.Attributes(ind_gaw_id).Value;
end
% environment and footprint
ind_gaw_id=contains({merged_data.Attributes.Name},'ebas_station_setting');
if sum(ind_gaw_id)>0
STN_st.env=merged_data.Attributes(ind_gaw_id).Value;
end