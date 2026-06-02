function [STN_rd, STN_st]=EbasNC2Timetable(filename)
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