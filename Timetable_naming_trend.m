function STN_rd=Timetable_naming_trend(merged_data)

% change the parameter's names from ebas into the naming convention for the
% trend analysis

%INPUT: merged_data=structure with data read from Netcdf file from ebas with all data
%from one station

% OUTPUT: STN_rd: timetable with correct parameter's names

% Exemple: STN_rd=Timetable_naming_trend(merged_data);


% read the wavelength (second dimention
Ni=contains({merged_data.Variables.Name},'Wavelength');
Wavelength=merged_data.Variables(Ni).Value;


%% compute time from ebas format (days from 1900 1 1)
% transform the structure with all wavelengths to a table with the right
% naming convention
for i=1:length(merged_data.Variables)
    if ~strcmp(merged_data.Variables(i).Name,'d_Wavelength')

        if strcmp(merged_data.Variables(i).Name,'time')
            STN_rds.Time=datetime(1900,1,1) + days(merged_data.Variables(i).Value); % compute time as datetime
        else

            for j=1:size(merged_data.Variables(i).Value,2) % take all data from an optical parameter

                % select variable type
                if contains(merged_data.Variables(i).Name,'_scattering')
                    var='Bs';
                elseif contains(merged_data.Variables(i).Name,'backscattering')
                    var='Bbs';
                elseif contains(merged_data.Variables(i).Name,'absorption')
                    var='Ba';
                elseif contains(merged_data.Variables(i).Name,'humidity')
                    var='U';
                else
                    warning('none of the variable has been identified')
                end

                % select size cut
                if contains(merged_data.Variables(i).Name,'pm10')
                    pm='0';
                elseif contains(merged_data.Variables(i).Name,'pm1')
                    pm='1';
                elseif contains(merged_data.Variables(i).Name,'pm2.5')
                    pm='2';
                else
                    warning('size cut is identified as TSP')
                    pm='';
                end

                %select wavelength
                if Wavelength(j)==370
                    wv='1';
                elseif Wavelength(j)==470
                    wv='2';
                elseif Wavelength(j)==520
                    wv='3';
                elseif Wavelength(j)==590
                    wv='4';
                elseif Wavelength(j)==660
                    wv='5';
                elseif Wavelength(j)==880
                    wv='6';
                elseif Wavelength(j)==950
                    wv='7';
                else
                    if Wavelength(j)<500 
                        wv='B';
                    elseif 500<Wavelength(j) && Wavelength(j)<600
                        wv='G';
                    elseif 600<=Wavelength(j) && Wavelength(j)<720
                        wv='R';
                    elseif Wavelength(j)>720
                        wv='Q';
                    end
                
                end

                % keep the paramerter only if there is data for the applied wavelength
                if sum(isnan(merged_data.Variables(i).Value(:,j)),1)< size(merged_data.Variables(i).Value(:,j),1)
                    name_var=strcat(var,wv,pm,'_I'); % parameter name
                    STN_rds.(name_var)=merged_data.Variables(i).Value(:,j);
                end
            end
        end
    end

end

% Convert the table into timetable
STN_rd=table2timetable(struct2table(STN_rds));