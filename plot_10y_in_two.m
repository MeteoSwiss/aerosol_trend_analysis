function plot_10y_in_two(STN_result, STN_st, seasonality)

%plot all the measured parameters in one fig (RH, scat, backscat,abs)
%plot all the calculated parameters in another one (Backscat, expS, expA,
%SSA
select1=STN_result(STN_result.length_period==10 & strcmp(STN_result.MK_seasonality,seasonality)==1 & startsWith(STN_result.parameter,["U";"Bs";"Bbs";"Ba";])==1 & startsWith(STN_result.parameter,'BbsF')==0,:);
names=unique(select1.parameter);
Tmin=min(select1.end_time)-0.5;
Tmax=max(select1.end_time)+0.5;
% measured parameters
f=figure;
hold on;

% %     select2=select1(strcmp(select1.parameter,(names{i}))==1,:);
% %     [s,ia]=unique(select2.end_time);
% %     select2=select2(ia,:);
subplot(4,1,1);
c=startsWith(names,'U');
if sum(c)>0
plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'units',seasonality);
ylabel('RH','FontSize',16);
end
title(STN_st.name,'FontSize',16);

subplot(4,1,2);
c=startsWith(names,'Bs');
if sum(c)>0
plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'units',seasonality);
ylabel('Bs','FontSize',16);
end

subplot(4,1,3);
c=startsWith(names,'Bbs');
if sum(c)>0
plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'units',seasonality);
ylabel('Bbs','FontSize',16);
end

subplot(4,1,4);
c=startsWith(names,'Ba');
if sum(c)>0
plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'units',seasonality);
ylabel('Ba','FontSize',16);
end
hold off;

%calculated parameters
select1=STN_result(STN_result.length_period==10 & strcmp(STN_result.MK_seasonality,seasonality)==1 & startsWith(STN_result.parameter,["BbsF";"exp";"SSA"])==1 ,:);
names=unique(select1.parameter);
g=figure;
hold on;
subplot(4,1,1);
c=startsWith(names,'BbsF');
if sum(c)>0
plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'units',seasonality);
ylabel('Backscat. fraction','FontSize',16);
title(STN_st.name,'FontSize',16);
end

subplot(4,1,2);
c=startsWith(names,'expS');
if sum(c)>0
plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'units',seasonality);
ylabel('Scat exponent','FontSize',16);
end

subplot(4,1,3);
c=startsWith(names,'expA');
if sum(c)>0
plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'units',seasonality);
ylabel('Abs. exponent','FontSize',16);
end

subplot(4,1,4);
c=startsWith(names,'SSA');
if sum(c)>0
plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'units',seasonality);
ylabel('SSA','FontSize',16);
end

% % % figure in % for scat, abs and SSA
% % select1=STN_result(STN_result.length_period==10 & strcmp(STN_result.MK_seasonality,'MetSea')==1 & startsWith(STN_result.parameter,["SSA";"Bs";"Ba";])==1 & startsWith(STN_result.parameter,'BbsF')==0,:);
% % names=unique(select1.parameter);
% % Tmin=min(select1.end_time)-0.5;
% % Tmax=max(select1.end_time)+0.5;
% % figure;
% % hold on;
% % c=startsWith(names,'Bs');
% % if sum(c)>0
% % plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'%');
% % ylabel('Bs','FontSize',16);
% % end
% % c=startsWith(names,'Ba');
% % if sum(c)>0
% % plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'%');
% % ylabel('Ba','FontSize',16);
% % end
% % c=startsWith(names,'SSA');
% % if sum(c)>0
% % plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,'%');
% % ylabel('SSA','FontSize',16);
% % end


if ispc
savefig(f,strcat('C:\github_trend\result\',(STN_st.name),'\',(STN_st.name),'_',(STN_st.name),'_MK_all10y_var.fig'));
savefig(g,strcat('C:\github_trend\result\',(STN_st.name),'\',(STN_st.name),'_',(STN_st.name),'_MK_all10y_cal.fig'));
    else
savefig(f,strcat('/prod/pay/Aerosol_actris_trend/trend_2026/result/',(STN_st.name),'/',(STN_st.name),'_',(STN_st.name),'_MK_all10y_var.fig'));
savefig(g,strcat('/prod/pay/Aerosol_actris_trend/trend_2026/result/',(STN_st.name),'/',(STN_st.name),'_',(STN_st.name),'_MK_all10y_cal.fig'));
   
end
%-------------------------------------
function plot_each_variable(select1, names,c,STN_st,Tmin,Tmax,type, seasonality)
N=names(c);
N=sortrows(N);
for i=1:length(N)
    select2=select1(strcmp(select1.parameter,(N{i}))==1,:);
    hold on;
    fig_seasonKendall_trend10_m(select2,STN_st,N(i),type,i,seasonality);
    hold on;
end
xlim([Tmin,Tmax]);
plot([Tmin Tmax],[0 0],'k-','LineWidth',2);

hl=findobj(gca,'type','line');
nb=(length(hl)-1)/i;
if i==2
    legend([hl(round(length(hl))) hl(round(length(hl)-nb))],[N(1) N(2)]);
elseif i==4
    legend([hl(round(length(hl))) hl(round(length(hl)-nb)) hl(round(length(hl)-2*nb)) hl(round(length(hl)-3*nb))],...
        [N(1) N(2) N(3) N(4)],'FontSize',12);
end
%-----------------------------
function fig_seasonKendall_trend10_m(data,data_st,namesP,type,nb, seasonality)
switch seasonality
    case 'MetSea'
colonne_res=5;
    case 'y'
        colonne_res=1;
    case 'month'
        colonne_res=13;
    otherwise
        error(' no plot because no seasonality was given')
end
sizeM=[8 8 8];
markerM=['o';'s';'v';'d';'o';'s';'v';'d';'o';'s';'v';'d'];
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
    if data.ss{i}(colonne_res)==90
        sizeR=sizeM(2);
        if data.slope{i}(colonne_res) >0
            colorT='r';
        elseif data.slope{i}(colonne_res) <0
            colorT='b';
        end
    elseif data.ss{i}(colonne_res)==95
        sizeR=sizeM(3);
        if data.slope{i}(colonne_res) >0
            colorT='r';
        elseif data.slope{i}(colonne_res) <0
            colorT='b';
            elseif data.slope{i}(colonne_res) ==0
            colorT=[0.7 0.7 0.7];
        elseif isnan(data.slope{i}(colonne_res))
            colorT='k';
        end
    else
        sizeR=sizeM(1);
        colorT='k';
    end
    hold on;
    plot(data.end_time(i),data.(s){i}(colonne_res),'Marker',markerM(nb),'Color',colorT,'MarkerSize',sizeR,'MarkerFaceColor',colorT,'DisplayName',(namesP{1}));
    hold on;
    
    %   line(repmat(data.end_time(i),2,1),[data.results{i}.(U)(5);data.results{i}.(L)(5)],'color',colorT);
    %     if data.results{i}.Xhomo(5)==1
    %         plot(data.end_time(i),data.results{i}.(s)(5),'rs','MarkerSize',sizeR);
    %     end
    
end

%legend(p(1),(namesP{1}));

grid on;

