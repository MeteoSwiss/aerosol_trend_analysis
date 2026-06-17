function [Tresult_MK,Tresult_GSMd,Tresult_LMSlog,Tresult_LMSlin,Tresult_LMSlog2,Tresult_LMSlin2]= all_trend_STN(data_tr,data_st)

%pour une station
%fait tous les trends necessaires

names=fieldnames(data_tr);
a=startsWith(names,["y";"Variable";"Time";"Proper"]);
names=names(~a);
for i=1:length(names)
    % give the resolution and instrument as a function of the name
    if startsWith((names{i}),["Bs";"Bbs";"expS"])
        inst='neph';

    elseif startsWith((names{i}),["Ba";"expA"])
        inst='abs';
    elseif startsWith((names{i}),'SSA')
        inst='abs+neph';
    elseif startsWith((names{i}),'U')
        inst='RH';
    elseif startsWith((names{i}),'N')
        inst='cpc';
        %pour stelios
    elseif startsWith((names{i}),'AOD')
        inst='pfr';
    elseif startsWith((names{i}),'AE')
        inst='pfr';
    end

    if startsWith((names{i}),["Bs";"Bbs";"Ba"]) && ~startsWith((names{i}),'BbsF')
        resolution=0.01;
        %distribution='log';
    elseif startsWith((names{i}),'SSA')
        resolution=0.005;
        %distribution='lin';
    elseif startsWith((names{i}),'BbsF')
        resolution=0.005;
        % distribution='lin';
    elseif startsWith((names{i}),'exp')
        resolution=0.01;
        %distribution='lin';
    elseif startsWith((names{i}),'U')
        resolution=0.1;
        % distribution='lin';
    elseif startsWith((names{i}),'N')
        resolution=1;
        % distribution='log';
    elseif startsWith((names{i}),'AOD')
        resolution=0.0002;
    elseif startsWith((names{i}),'AE')
        resolution=0.01;
    end
    % restrict the time series to period with data
    ind=~isnan(data_tr.(names{i}));
    s=find(ind,1);
    e=find(ind,1,'last');
    PP=timerange(data_tr.Time(s),data_tr.Time(e));
    data_trok=data_tr(PP,:);

    %compute the overall trend with end year = last year for all data!

    [Tresult_MK_25,Tresult_GSMdi,Tresult_LMSlog,Tresult_LMSlin,Tresult_LMSlog2,Tresult_LMSlin2]=all3_trend(data_trok,{names{i}}, inst, data_st.name, resolution,'end_year',max(data_trok.y), 'fig',1);
    %     end
    %compute all the 10y trend

    s=data_tr.y(s);
    e=data_tr.y(e);
    nb_trend= e-s-8;
    if nb_trend>1
        % % % for j=1:nb_trend
        % % %     if e==max(data_trok.y)
        % % %         if j==2
        % % %             Tresult_MK_10y=seasonalKendall_main_D(data_trok,{names{i}}, inst, data_st.name, resolution, 'period',10,'end_year',e-j+1);
        % % %         elseif e-j+1~=max(data_trok.y)
        % % %             T=seasonalKendall_main_D(data_trok,{names{i}}, inst, data_st.name, resolution,'period',10,'end_year',e-j+1);
        % % %
        % % %             Tresult_MK_10y=[Tresult_MK_10y;T];
        % % %         end
        % % %     else
        % % %         if j==1
        % % %             Tresult_MK_10y=seasonalKendall_main_D(data_trok,{names{i}}, inst, data_st.name, resolution, 'period',10,'end_year',e-j+1);
        % % %         elseif e-j+1~=max(data_trok.y)
        % % %             T=seasonalKendall_main_D(data_trok,{names{i}}, inst, data_st.name, resolution, 'period',10,'end_year',e-j+1);
        % % %             Tresult_MK_10y=[Tresult_MK_10y;T];
        % % %         end
        % % %     end
        % % % end
        for j=2:nb_trend
            if j==2
                Tresult_MK_10y=seasonalKendall_main_D(data_trok,{names{i}}, inst, data_st.name, resolution, 'period',10,'end_year',e-j+1);
            else
                T=seasonalKendall_main_D(data_trok,{names{i}}, inst, data_st.name, resolution, 'period',10,'end_year',e-j+1);
                Tresult_MK_10y=[Tresult_MK_10y;T];
            end
        end

        if ~isempty( Tresult_MK_25)
            Tresult_MKi=[Tresult_MK_25;Tresult_MK_10y];
        else
            Tresult_MKi=Tresult_10y;
        end
    elseif nb_trend==1
        Tresult_MKi=Tresult_MK_25;
        % % %I think this is already done
        % % elseif nb_trend<1 %too short time series not ending in 2025, but trend should anyhow be calculated
        % %     % Tresulti=all3_trend(data_trok,{names{i}}, inst, data_st.name, resolution, 'period',10,'end_year',max(data_trok.y));
        % %     [Tresult_MK,Tresult_GSMd,Tresult_LMSlog,Tresult_LMSlin]=all3_trend(data_trok,{names{i}}, inst, data_st.name, resolution,'period',10,'end_year',max(data_trok.y), 'fig',1);
    end

    if i==1
        Tresult_MK=Tresult_MKi;
        % Tresult_GSMd=Tresult_GSMdi;
        % Tresult_LMSlog=Tresult_LMSlogi;
        % Tresult_LMSlin=Tresult_LMSlini;
    else
        Tresult_MK=[Tresult_MK;Tresult_MKi];
        % Tresult_GSMd=[Tresult_GSMd;Tresult_GSMdi];
        % Tresult_LMSlog=[Tresult_LMSlog;Tresult_LMSlogi];
        % Tresult_LMSlin=[Tresult_LMSlin;Tresult_LMSlini];
    end
end

%names_res=strcat((data_st.name),'_result_MK');
eval([strcat(data_st.name,'_result_MK'),'=Tresult_MK']);
eval([strcat(data_st.name,'_result_GSMd'),'=Tresult_GSMd']);
eval([strcat(data_st.name,'_result_LMSlog'),'=Tresult_LMSlog']);
eval([strcat(data_st.name,'_resultLMSlin'),'=Tresult_LMSlin']);
eval([strcat(data_st.name,'_result_LMSlog2'),'=Tresult_LMSlog2']);
eval([strcat(data_st.name,'_resultLMSlin2'),'=Tresult_LMSlin2']);

%-----------------------------
function fig_seasonKendall_trend10(data,data_st,namesP,type)
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
for i=1:size(data,1)
    if data.results{i}.ss(5)==90
        sizeR=sizeM(2);
        if data.results{i}.slope(5) >0
            colorT='r';
        elseif data.results{i}.slope(5) <0
            colorT='b';
        end
    elseif data.results{i}.ss(5)==95
        sizeR=sizeM(3);
        if data.results{i}.slope(5) >0
            colorT='r';
        elseif data.results{i}.slope(5) <0
            colorT='b';
        end
    else
        sizeR=sizeM(1);
        colorT='k';
    end
    plot(data.end_time(i),data.results{i}.(s)(5),'Marker','.','Color',colorT,'MarkerSize',sizeR);
    hold on;

    line(repmat(data.end_time(i),2,1),[data.results{i}.(U)(5);data.results{i}.(L)(5)],'color',colorT);
    if data.results{i}.Xhomo(5)==1
        plot(data.end_time(i),data.results{i}.(s)(5),'rs','MarkerSize',sizeR);
    end
    plot([min(data.end_time)-0.2 max(data.end_time)+0.2],[0 0],'r-');
end
ylabel(strcat(data_st.name,'_', namesP));
grid on;
% if result_y.ss==95
%     a=sizeM(3);
% elseif result_y.ss==90
%     a=sizeM(2);
% else
%     a=sizeM(1);
% end
% plot(t(end)+delta/2,result_y.(s),'g.','MarkerSize',a);
% line([t(end)+delta/2 t(end)+delta/2 ],[result_y.(U);result_y.(L)],'color','g');