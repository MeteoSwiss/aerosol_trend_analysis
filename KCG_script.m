KCG_st.name='KCG'; 
KCG_st.lat=-40.68;
KCG_st.lon=144.69;
KCG_st.alt=94;
KCG_st.env='coast';
KCG_st.footp='RB';
%%
% % MAAP=table2timetable(MAAPdata);
% % Neph=table2timetable(Nephall);
% % % % sum(MAAP.numflag>0)
% % % % figure
% % % % plot(MAAP.DateGMThhmm,MAAP.numflag,'.b');
% % % % hold on;
% % % % plot(Neph.DateGMThhmm,Neph.flag,'.r');
% % % % %flag==1
%%
% % % remove NAN
% % MAAP.abs670(MAAP.numflag>0.99)=NaN;
% % Neph.RH_int(Neph.flag>0.9)=NaN;
% % Neph.sc450(Neph.flag>0.9)=NaN;
% % Neph.sc525(Neph.flag>0.9)=NaN;
% % Neph.sc635(Neph.flag>0.9)=NaN;
% % %remove not used variable
% % namesM=fieldnames(MAAP);
% % c=strcmp(namesM,'abs670')==1 | strcmp(namesM,'Properties')==1 |startsWith(namesM,["D";"V"]) ;
% % NM=namesM(~c);
% % for i=1:length(NM)
% %     MAAP.(NM{i})=[];
% % end
% % 
% % namesN=fieldnames(Neph);
% % c=strcmp(namesN,'sc450')==1 |strcmp(namesN,'sc525')==1|strcmp(namesN,'sc635')==1 | strcmp(namesN,'Properties')==1 |startsWith(namesN,["D";"V";"R"]) ;
% % NN=namesN(~c);
% % for i=1:length(NN)
% %     Neph.(NN{i})=[];
% % end
%%
% % KCG_rd=timetable(MAAP.DateGMThhmm);
% % KCG_rd.Time=MAAP.DateGMThhmm;
% % KCG_rd.BaR=MAAP.abs670;
% % KCG_rd=synchronize(KCG_rd,Neph);
% % KCG_rd.Properties.VariableNames={'BaR_A';'U';'BsB_S';'BsG_S';'BsR_S'};
%%
KCG_rd=read_betsy_2026('C:\github_trend\raw_data\kcg_maap_eco_polar_aeth_new', KCG_st.name); %[output:760372bd]

names_rd=fieldnames(KCG_rd);
%%
%[text] ## HOURLY
%Do hourly homogenizse KPS data
KCG_rd_h=read_betsy_2026('C:\github_trend\raw_data\kcg_maap_eco_polar_aeth_new', KCG_st.name);

names_rd=fieldnames(KCG_rd);
%%
datevec(KCG_rd.Time(1)) % 2008
datevec(KCG_rd.Time(end)) % 2025
 prctile(KCG_rd.RHsamp_S12_po,[5 50 95]) % 44%
  prctile(KCG_rd.rhsamp_S11_eco,[5 50 95]) % 40%
plotFigControl(KCG_rd,KCG_st.name);
%%
% homogeneise Bs eco+polar NOT RIGHT THE ORDER OF THE VAR
% common period jan 2018-may 2023
% the period 20.2.2021-16.8.2021 has much higher eco value: do not use
% try both periods before and after
Pt=timerange('2021-09-01','2023-05-01');
ratio_med_h=nanmedian(KCG_rd.BsG_S11_eco(Pt)./KCG_rd.BsG_S12_po(Pt)) % 1) 1.2529 2) 1.2226
ratio_STN_h=nanstd(KCG_rd.BsG_S11_eco(Pt)./KCG_rd.BsG_S12_po(Pt)) % 1) 0.8131 2) 0.9772
ratio_mean_h=nanmean(KCG_rd.BsG_S11_eco(Pt)./KCG_rd.BsG_S12_po(Pt)) % 1) 1.3934 2)1.3803
ratio_mean=nanmean(KCG_rd.BsG_S11_eco(Pt))/nanmean(KCG_rd.BsG_S12_po(Pt)) % 1) 1.1304  2)1.2279
ratio_median=nanmedian(KCG_rd.BsG_S11_eco(Pt))/nanmedian(KCG_rd.BsG_S12_po(Pt)) % 1) 1.2103 2) 1.2934
figure;
plot(KCG_rd.BsG_S11_eco(Pt),KCG_rd.BsG_S12_po(Pt),'.');

% fit between both data: 1) slope 1.013x+0.9224 , 1/slope=0.9872 2)
% 0.8921x-0.2811 1/slope=1.1210
fit_test=robustfit(KCG_rd.BsG_S11_eco(Pt),KCG_rd.BsG_S12_po(Pt)); % 1) y coor=-0.4878, slope=0.8967, 1/slope=1.1152
% 2) 0.4965  0.7159, 1/slope=1.3968
% 
%%
c_eco_po=1/1.10;
%%
P1=timerange('2009-01-01','2019-01-01');
P2=timerange('2019-01-01','2026-01-01');

KCG_rd.BsG_Sh=NaN(size(KCG_rd.BsG_S12_po));
KCG_rd.BsG_Sh(P1)=KCG_rd.BsG_S11_eco(P1).*c_eco_po;
KCG_rd.BsG_Sh(P2)=KCG_rd.BsG_S12_po(P2);


KCG_rd.BsB_Sh=NaN(size(KCG_rd.BsB_S12_po));
KCG_rd.BsB_Sh(P1)=KCG_rd.BsB_S11_eco(P1).*c_eco_po;
KCG_rd.BsB_Sh(P2)=KCG_rd.BsB_S12_po(P2);

KCG_rd.BsR_Sh=NaN(size(KCG_rd.BsR_S12_po));
KCG_rd.BsR_Sh(P1)=KCG_rd.BsR_S11_eco(P1).*c_eco_po;
KCG_rd.BsR_Sh(P2)=KCG_rd.BsR_S12_po(P2);
%Bbs
% KCG_rd.BbsG_Sh=NaN(size(KCG_rd.BbsG_S12_po));
% KCG_rd.BbsG_Sh(P1)=KCG_rd.BbsG_S11_eco(P1).*c_eco_po;
% KCG_rd.BbsG_Sh(P2)=KCG_rd.BbsG_S12_po(P2);
% 
% 
% KCG_rd.BbsB_Sh=NaN(size(KCG_rd.BbsB_S12_po));
% KCG_rd.BbsB_Sh(P1)=KCG_rd.BbsB_S11_eco(P1).*c_eco_po;
% KCG_rd.BbsB_Sh(P2)=KCG_rd.BbsB_S12_po(P2);
% 
% KCG_rd.BbsR_Sh=NaN(size(KCG_rd.BbsR_S12_po));
% KCG_rd.BbsR_Sh(P1)=KCG_rd.BbsR_S11_eco(P1).*c_eco_po;
% KCG_rd.BbsR_Sh(P2)=KCG_rd.BbsR_S12_po(P2);
%%
%% HOURLY
KCG_rd_h.BsG_Sh=NaN(size(KCG_rd_h.BsG_S12_po));
KCG_rd_h.BsG_Sh(P1)=KCG_rd_h.BsG_S11_eco(P1).*c_eco_po;
KCG_rd_h.BsG_Sh(P2)=KCG_rd_h.BsG_S12_po(P2);


KCG_rd_h.BsB_Sh=NaN(size(KCG_rd_h.BsB_S12_po));
KCG_rd_h.BsB_Sh(P1)=KCG_rd_h.BsB_S11_eco(P1).*c_eco_po;
KCG_rd_h.BsB_Sh(P2)=KCG_rd_h.BsB_S12_po(P2);

KCG_rd_h.BsR_Sh=NaN(size(KCG_rd_h.BsR_S12_po));
KCG_rd_h.BsR_Sh(P1)=KCG_rd_h.BsR_S11_eco(P1).*c_eco_po;
KCG_rd_h.BsR_Sh(P2)=KCG_rd_h.BsR_S12_po(P2);
p1=timerange('2010-05-13','2012-02-21');
p4=timerange('2016-02-27','2016-08-16');

names_sc={'BsB_Sh', 'BsG_Sh', 'BsR_Sh'};
for i=1:length(names_sc)
  
    KCG_rd_h.(names_sc{i})(p1)=NaN;
   % KCG_rd.(names_sc{i})(p2)=NaN;
  % KCG_rd.(names_sc{i})(p3)=NaN;
end
%%
% following 2019 to invalidate and take over from 2019 data
% in 2026 Melita said to invalidate from mai 2010 until Feb 2012
p1=timerange('2010-05-13','2012-02-21');
p4=timerange('2016-02-27','2016-08-16');

names_sc={'BsB_Sh', 'BsG_Sh', 'BsR_Sh'};
for i=1:length(names_sc)
   KCG_rd.(names_sc{i})(p1)=NaN;
  
   % KCG_rd.(names_sc{i})(p2)=NaN;
  % KCG_rd.(names_sc{i})(p3)=NaN;
end

 % KCG_rd.BsB_Sh(p4)=KCG_rdd.BsB_S(p4);
 %  KCG_rd.BsG_Sh(p4)=KCG_rdd.BsG_S(p4);
 %   KCG_rd.BsR_Sh(p4)=KCG_rdd.BsR_S(p4);
%%
% BP detection for scat coef and backscat

names_sc={'BsB_Sh', 'BsG_Sh', 'BsR_Sh'};
break_KCG_scat_BP=change_point_analysis_Def(KCG_rd,names_sc,0.05,'KCG','scattering');
T_KCG_scat_BP=make_table_breakpoints_def(break_KCG_scat_BP);
%%

names_bsc={'BbsB_S12_po', 'BbsG_S12_po', 'BbsR_S12_po'};
break_KCG_bscat_BP=change_point_analysis_Def(KCG_rd,names_bsc,0.05,'KCG','Backscattering');
T_KCG_bscat_BP=make_table_breakpoints_def(break_KCG_bscat_BP);
%%
% homo MAAP-AE33
KCG_rd.BaR_A_maap=KCG_rd.BC_maap.*6.6.*1.05e-3;
%%
% compute abs from AE3 myself assuming M8060, C=1.39 applied
% should be: babs=bATN/C and BC=babs/sigma
% their R program: abs=BC*sigma*1e-3 and BC=ATN/(flow*(1-z) *C*sigma)*1e9
% % % sigma=[18.47 14.54 13.14 11.58 10.35 7.77 7.19]; %for the 7 wavelengths with an exponent of -1
% % % H=1.76;
% % % KCG_rd.Bax1_A82=KCG_rd.BC1_A12_ae33 .*sigma(1);
% % % KCG_rd.Bax2_A82=KCG_rd.BC2_A12_ae33 .*sigma(2);
% % % KCG_rd.Bax3_A82=KCG_rd.BC3_A12_ae33 .*sigma(3);
% % % KCG_rd.Bax4_A82=KCG_rd.BC4_A12_ae33 .*sigma(4);
% % % KCG_rd.Bax5_A82=KCG_rd.BC5_A12_ae33 .*sigma(5);
% % % KCG_rd.Bax6_A82=KCG_rd.BC6_A12_ae33 .*sigma(6);
% % % KCG_rd.Bax7_A82=KCG_rd.BC7_A12_ae33 .*sigma(7);
%%

%Babs harmonised=Basb/1.76 for M8060
% --> let divide first abs AE33 by 1.76
H=1.76;
KCG_rd.abs1_A12_ae33=KCG_rd.abs1_A12_ae33./H;
KCG_rd.abs2_A12_ae33=KCG_rd.abs2_A12_ae33./H;
KCG_rd.abs3_A12_ae33=KCG_rd.abs3_A12_ae33./H;
KCG_rd.abs4_A12_ae33=KCG_rd.abs4_A12_ae33./H;
KCG_rd.abs5_A12_ae33=KCG_rd.abs5_A12_ae33./H;
KCG_rd.abs6_A12_ae33=KCG_rd.abs6_A12_ae33./H;
KCG_rd.abs7_A12_ae33=KCG_rd.abs7_A12_ae33./H;
%%
%[text] ## HOURLY

KCG_rd_h.abs1_A12_ae33=KCG_rd_h.abs1_A12_ae33./H;
KCG_rd_h.abs2_A12_ae33=KCG_rd_h.abs2_A12_ae33./H;
KCG_rd_h.abs3_A12_ae33=KCG_rd_h.abs3_A12_ae33./H;
KCG_rd_h.abs4_A12_ae33=KCG_rd_h.abs4_A12_ae33./H;
KCG_rd_h.abs5_A12_ae33=KCG_rd_h.abs5_A12_ae33./H;
KCG_rd_h.abs6_A12_ae33=KCG_rd_h.abs6_A12_ae33./H;
KCG_rd_h.abs7_A12_ae33=KCG_rd_h.abs7_A12_ae33./H;

KCG_rd_h.BaR_A_maap=KCG_rd_h.BC_maap.*6.6.*1.05e-3;
%%
% common period jan 2022-June 2024
% pummp of MAAP is no more good after March (or August) 2022 --> factor
% quite similar for mean and median, slope robust fit =1.08
% if I end the period on 1.3.2022, the robust fit slope =1.07 
% try just to keep H factor without homogeneisation, stop MAAP in Jan 2022

Pt=timerange('2022-01-01','2022-03-01');
ratio_med_h=nanmedian(KCG_rd.abs5_A12_ae33(Pt)./KCG_rd.BaR_A_maap(Pt)) % 1.27 %[output:60e2e113]
ratio_STN_h=nanstd(KCG_rd.abs5_A12_ae33(Pt)./KCG_rd.BaR_A_maap(Pt)) % 1.17 %[output:9c069efc]
ratio_mean_h=nanmean(KCG_rd.abs5_A12_ae33(Pt)./KCG_rd.BaR_A_maap(Pt)) % 1.51 %[output:16df2019]
ratio_mean=nanmean(KCG_rd.abs5_A12_ae33(Pt))/nanmean(KCG_rd.BaR_A_maap(Pt)) %1.05 %[output:2aefa845]
ratio_median=nanmedian(KCG_rd.abs5_A12_ae33(Pt))/nanmedian(KCG_rd.BaR_A_maap(Pt)) % 1.344 %[output:86f8bdcb]
figure; %[output:25aa6aa7]
plot(KCG_rd.BaR_A_maap(Pt),KCG_rd.abs5_A12_ae33(Pt),'.'); %[output:25aa6aa7]

% fit between both data: 1) slope 1.524x+0.04078 , 
fit_test=robustfit(KCG_rd.BaR_A_maap(Pt),KCG_rd.abs5_A12_ae33(Pt)); % 1) y coor=-0.0192, slope=0.963 %[output:8bf6a4cd]
% 
%%
c_maap_ae=1;
%%
P1=timerange('2008-01-01','2022-01-01');
P2=timerange('2022-01-01','2026-01-01');

KCG_rd.BaR_Ah=NaN(size(KCG_rd.BaR_A_maap));
KCG_rd.BaR_Ah(P1)=KCG_rd.BaR_A_maap(P1).*c_maap_ae;
KCG_rd.BaR_Ah(P2)=KCG_rd.abs5_A12_ae33(P2);
%%
%[text] ## HOURLY
P1=timerange('2008-01-01','2022-01-01');
P2=timerange('2022-01-01','2026-01-01');

KCG_rd_h.BaR_Ah=NaN(size(KCG_rd_h.BaR_A_maap));
KCG_rd_h.BaR_Ah(P1)=KCG_rd_h.BaR_A_maap(P1).*c_maap_ae;
KCG_rd_h.BaR_Ah(P2)=KCG_rd_h.abs5_A12_ae33(P2);
%%
% BP detection for abs
names_abs={'BaR_A_maap','BaR_Ah'};
break_KCG_abs_BP=change_point_analysis_Def(KCG_rd,names_abs,0.05,'KCG','Absorption'); %[output:8e1419b2] %[output:9d15429d] %[output:23a2d95e] %[output:21e0df38]
T_KCG_abs_BP=make_table_breakpoints_def(break_KCG_abs_BP);
%%
lambdaSC=[450;525;365];

lambdaAE7=[370 470 520 590 660 880 950];
names_neph={'BsB_Sh', 'BsG_Sh', 'BsR_Sh'};
KCG_cal=compute_exp_D(KCG_rd,names_neph,lambdaSC);
% KCG_expA=compute_exp_D(KCG_rd,names_abs,lambdaAE); not long enough

KCG_cal.SSAR=KCG_rd.BsG_Sh./(KCG_rd.BsG_Sh+KCG_rd.BaR_Ah.*(525/670));
KCG_cal.SSAR_maap=KCG_rd.BsG_Sh./(KCG_rd.BsG_Sh+KCG_rd.BaR_A_maap.*(525/670));

plotFigControl_cal(KCG_cal, KCG_st.name); %[output:030633c6] %[output:165b28c0] %[output:59fd8b89]
%%
%[text] ## HOURLY
KCG_cal_h=compute_exp_D(KCG_rd_h,names_neph,lambdaSC);
% KCG_expA=compute_exp_D(KCG_rd_h,names_abs,lambdaAE); not long enough

KCG_cal_h.SSAG=KCG_rd_h.BsG_Sh./(KCG_rd_h.BsG_Sh+KCG_rd_h.BaR_Ah.*(525/670));
KCG_cal_h.SSAG_maap=KCG_rd_h.BsG_Sh./(KCG_rd_h.BsG_Sh+KCG_rd_h.BaR_A_maap.*(525/670));

KCG_all=outerjoin(KCG_rd_h,KCG_cal_h);

pp3=timerange('2008-01-01','2016-01-01');
KCG_all.expS_bg(pp3)=NaN;

%%
% BP detection for expS, Backscat fraction, expA and SSA

names_exps={'expS_bg'};
break_KCG_exps_BP=change_point_analysis_Def(KCG_cal,names_exps,0.05,'KCG','scattering exp');
T_KCG_exps_BP=make_table_breakpoints_def(break_KCG_exps_BP);

% names_exps={'expA_bg', 'expA_br', 'expA_gr'};
% break_KCG_exps_BP=change_point_analysis_Def(KCG_cal,names_exps,0.05,'KCG','scattering exp');
% T_KCG_exps_BP=make_table_breakpoints_def(break_KCG_exps_BP);
% 
% names_BbF={'BbsFg'};
% break_KCG_BbF_BP=change_point_analysis_Def(KCG_cal,names_BbF,0.05,'KCG','backscat fraction');
% T_KCG_BbF_BP=make_table_breakpoints_def(break_KCG_BbF_BP);

names_SSA={'SSAR'};
break_KCG_SSA_BP=change_point_analysis_Def(KCG_cal,names_SSA,0.05,'KCG','SSA');
T_KCG_SSA_BP=make_table_breakpoints_def(break_KCG_SSA_BP);


%%
KCG_tr=outerjoin(KCG_rd, KCG_cal);
KCG_tr.y=year(KCG_tr.Time);
% begin at the beginning of a year:abs:2008, scat:2009
%end: 2018
P=timerange('2008-01-01','2026-01-01');
KCG_tr=KCG_tr(P,:);

% variables to delete
names=fieldnames(KCG_tr);
c=startsWith(names,["rh","RH","w","am", "p","t","BC","abs","AAE","Bbs","Bax"]) | endsWith(names,["_eco","_po","_br","_gr","maap"]) ;
N=names(c);
for i=1:length(N)
    KCG_tr.(N{i})=[];
end


%no neph data in 2009 at all
%no need to analyse dry neph

%neph problems 2010-2012:
% following answer to questions:invalid: 24.5.2010-1.8.2010 ,
% invalid october 2010-7.3.2011 and invalid: 28.5.2011-21.2.2012
% all is already invalidated in KCG_rd
% Melita said: invalidate May 2010-Feb 2012 (i.e. 24.5.2010-21.2.2012
% New data --> no need to invalidate neph 14.11.2014 to 31.12.2014 TEST NEW DATA
% 2016 NEW DATA fine



%compute trend only on green
KCG_tr.BsB_Sh=[];
KCG_tr.BsR_Sh=[];


%compute exp only on br and only since 2016
% remove calibration problem
% % pp1=timerange('2018-11-15','2019-01-01');
% % pp2=timerange('2016-03-07','2016-11-11');

pp3=timerange('2008-01-01','2016-01-01');
KCG_tr.expS_bg(pp3)=NaN;

%%
%after trend analysis, decide if trend or only dlm should be used!!
[KCG_result_MK,KCG_result_LMSlog,KCG_result_LMSlin]=all_trend_STN(KCG_tr,KCG_st); %[output:8abc7834] %[output:2306dc04] %[output:4d175c96] %[output:849e5967] %[output:73560081] %[output:009e7c8a] %[output:630232f8] %[output:7ffd71de] %[output:00782a05] %[output:58f06463] %[output:1ea3c794] %[output:224e6560] %[output:3369a97d] %[output:4de04fb7] %[output:7f63dcb3] %[output:5662530a] %[output:4b00554e] %[output:0dc135ce] %[output:2fad7a48] %[output:648473e8] %[output:57c86ac5] %[output:9ba05856] %[output:35fd3337] %[output:78b1b4d4] %[output:092d462d] %[output:81b97be6] %[output:765a963b] %[output:99e01b07] %[output:0896794e] %[output:882b6f8a] %[output:0e5af49d] %[output:7b9730be] %[output:5d56c84c] %[output:2a538191] %[output:512f4f2b] %[output:9dd71502] %[output:4baeac86] %[output:3dbf4763] %[output:734812ac] %[output:984dcb31] %[output:34874f36] %[output:096320c5] %[output:1c201b1e] %[output:90daa2ad] %[output:193d03a4] %[output:904868d3] %[output:298201df] %[output:9cf454ab] %[output:2bb1bda0] %[output:17cc76da] %[output:14905faf] %[output:7685ffc9] %[output:8fa502e8] %[output:80fd9422] %[output:57b82619] %[output:21a5cb93] %[output:162af99d] %[output:36dc2db9] %[output:4c54fce6] %[output:02cdc6f7] %[output:0f1f6a01] %[output:7714f6f6] %[output:95a7c342] %[output:60e44764] %[output:2d56ee40] %[output:740ba8cf] %[output:0ceea8f4] %[output:9359e771] %[output:355f85f0] %[output:4cfffd88]

writetable(KCG_result_MK,'KCG_res_MK_present_10y.txt'); %, 'delimiter',',' )
writetable(KCG_result_LMSlog,'KCG_res_LMSlog.txt'); 
writetable(KCG_result_LMSlin,'KCG_res_LMSlin.txt'); 
plot_10y_in_two(KCG_result_MK, KCG_st,'y'); %[output:09faf133] %[output:01cbb4fd]
%%
 timetable_to_netcdf(KCG_tr,'KCG_tr_def.nc'); %[output:722a53f7]
  timetable_to_netcdf(KCG_all,'KCG_all_csiro.nc');

%%
%Do hourly homogenizse KPS data
KCG_rd_h=read_betsy_2026('C:\github_trend\raw_data\kcg_maap_eco_polar_aeth_new', KCG_st.name); %[output:737ce877]

names_rd=fieldnames(KCG_rd);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":20.7}
%---
%[output:760372bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: there is a lot of days with less than 50% data coverage"}}
%---
%[output:60e2e113]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_med_h","value":"1.2740"}}
%---
%[output:9c069efc]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_STN_h","value":"0.8997"}}
%---
%[output:16df2019]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean_h","value":"1.5365"}}
%---
%[output:2aefa845]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean","value":"1.2789"}}
%---
%[output:86f8bdcb]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_median","value":"1.4672"}}
%---
%[output:25aa6aa7]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAALYAAABuCAYAAACKnelGAAAAAXNSR0IArs4c6QAAChJJREFUeF7tnT1sFUcQx4c0wSiKwAkNNsGS5ZoEF3aclFE604Bkm4ZYFERKQQH+wpEQUmxsIxdGihSEEKLBUFA5VZSGIg4uopjaimTApopMFEXBaXA0p6yzPu7e7r6b+9jb\/0mW4L29u5n\/\/N7c3Ly3e\/t2dnZ2CBsUqJkC+wB2zSIKdyIFADZAqKUCALuWYYVTRrC3trbo3Llz9OTJE+rv76eZmRlqaWlJVO7Bgwe0vLzccAwkhwJFKGAEe3Z2ljo6OujkyZM0Pj5OQ0ND1Nvb+4Zta2trNDw8TN3d3QC7iMjhHA0VaAi2ytZjY2MRzJyR19fXif+vb69evaKpqSnq7Oyk1dVVgA3oSlfACPbFixfp8uXL1NXVFYGdVGrw67wdO3aMFhcX94C9sbFB\/IcNCtgo0N7eTvyXdcsMNpcgd+\/epcnJyagO18FmoEdGRmhlZSWrndg\/EAV6enro+vXrmeE2gs03jo1KEc7WExMTe2Tnm8yFhQV6\/PgxnTlzJjK0ra3Nm9Bsbm5GVxkW2bfNZ9s5ATI39+7dS7yPc4mF2M0jn5RB1jO2AlvCUBenMNZPBSR5MYKtt\/vOnz8fZW9+7erVq3TlyhVqbW3dVRFg+wlUVawuFOwsTksamsUO7OuHApK8GDN2FkkkDc1iB\/b1QwFJXgC2HzEPwkqAHUSYw3MSYIcX8yA8BthBhDk8JwF2eDGvrcf93\/5Kz19u09FD+2my+5\/oCz2J7z1w81hbZPxw7MNvfqZnW9v0Qet++u6zfQDbj7DBSpMCKmPzOIBtUgvve6kAamwvwwajTQoAbJNCeN9LBQC2l2GD0SYFALZJIbzvpQKFgm0zS50n\/N68eTMS89q1azQwMBD9W9JQLyMFo50UkOTF2Mc2zVLXf4PNs05GR0dpbm4umiMpaaiTQhjspQKSvGSeGqYryPMfp6enaX5+PpqAIGmol5GC0U4KSPJiBNtmljpbr8qRpFKE5zyeOnXKyUkMDk+BSoLNYeD1RXhRnb6+vqjOVoZeuHCB+A8bFGikAE\/kLWQyr+2CObqxnLl547mRkp9AIFF\/BR4+fBgt11HIj6Bsbh4fPXq0O8lXX64BYNcfRkkPJXkxdkVsZqmj3ScZ3nCPVSjYWWSWNDSLHdjXDwUkeTFm7CySSBqaxQ7s64cCkrwAbD9iHoSVADuIMIfnJMAOL+ZBeAywgwhzeE4C7PBiHoTHADuIMIfnJMAOL+ZBeAywgwhzeE4C7PBiHoTHADuIMIfnJMAOL+ZBeAywgwhzeE4WCrZplrqaNbO0tBRFQj2Aif8taWh4YQ7PY0lejD+CMk000B9DrT4Eg4ODe6aGScyICC\/M4XlcGNjNTg3r6OgA2OFxmdnjQsG2naXOXvEHQR+vDJWepa4vFr701UeZBcUBqqFAJcFWtfbQ0NDu44LzmqWuLxa++vXH1YgKrMisQOVmqccztfJQ8hOoq6YvFl4nsEO\/ElVqlnra46fRFXFPYKFfiSQTobErYpqlfuvWrd0FKVUo1WpQkoa6YyK7RxHZtK5XIttISPJiBNvWqKRxkoZmsUNi39CzqYSGpmNI8gKwTWr\/937o2dRSpkzDALajfEWUEY4mYXiCAgDbEQuUEY6ClTQcYDsIz9n6p9\/+iB6QyVud2oMOMngxFGA7hAnZ2kGskocCbIcA4KbPQayShwLskgOA0+ejAMDOR1cctWQFAHbJAXA5PVqN9moBbHutSh+Jm1f7EABse61KH4mbV\/sQBA02Lu32oPg2MmiwcWn3DVd7ewsF2zRLXZndaAaN5GReXNrtQfFtZKFgm2aps3g6\/DrEkob6FiTY666AJC\/GR07rz23Ul1rQM\/WNGzfo9OnTdOnSpeh5j729vdHbkoa6y5S8R7M1erP7SdkdwnEkeTGCbTtLPb5Ugw629Cz1ZoKswHy2tR3tzj+KcvlBFGr7ZlR328c7sKvwLHUdTCW3C9io7d0gbWZ05Wap63V2VUsRgNkMauZ9JEu0Ss1SV643KkUkuyJ8PhcxXcaaw4gRcQUkS7TCSpF4x0MtOJm05IIk2CYYXcR0GQts3RWQvBIWCra7q\/\/v0ayhJhhdxHQZm8VX7JtdgWZ5STpzJWepS8BoyvrZw4AjSCtQW7BNMKa9H39dzXNspq0nHSwcz16B2oFt22NOK1HU6wpkvVfNr7m09ezDgJHSCtQObB3MRlk2rUSJfzDUMQC1NHr5Hq9WYOvLI6hMy\/J90nmQ4mtf25YqSn7bTG06br7hxNGVArUCWy8vjh7aH60BomdtHbrnL7eJ4Vdfh0sBaerCAL1iFKgN2EmL2cTLjfhX4TrYUkBKdGGKCX29z1IbsPXamksP3jgrq40zePy1+I2ha9lRbzT89s5LsKd+eTuClmFVtXO8Lcdh0evseK1t2z3xO7zhWu8l2F\/+uJNYH+sZWK+xVXiTfl6K0qGe8HsH9sHBBVrffmfPwpDxFp\/KzvoNIsqMegKc5pV3YP\/5+Sy9PvD+nm6G6n6kOek6ESAsBOrprXdgc8ZmsHnj3nLaN4V6uAB2PeFt5FXlwFYGsdHqwUr8b\/W6AlvdOKbdNKo1rFGChAe1zovE7\/cz\/7pPf8YjGzc9PU3z8\/PU2tq6C3a8FOFxSR0OldH1sEp9CRMmKn55XamMzcbwEg23b9+mlpYWGh8fJ\/V0XmXoX5+O0usD70Uqf\/\/F0T1qq27JW3\/\/\/sZ7PND0fl6hW1lZofb2dmpra8vrFLkd11fbNzc3aWRkhCqRsRnexcVFmpmZiQLFYPf19dHAwABtbGxEhrLQaZsO\/bs\/jL0xzPR+bnTgwKUo0NPTQ7yqASeVLFvmUqQR2GwYw81\/2KCAjQIMdFao+TwiYKeVIjaOYAwUyEOBzGA3unnMw2AcEwrYKJAZbD6J3u6TKPxtDMcYKNBIARGwk06Q1tuuWjjUKrFLS0t0\/PjxqLvDrco0n9SNMneAyt5sbNfjcOTIEbpz5w51dXWVbTrZ2L62tkbDw8P04sULY2ziDuUCtk\/lib7QplpZljs68Y3HTUxMUH9\/f9QBqgLYJtv1ODDMPP7+\/fsNP7xFEW+yPb4sNceGN15pzGbLBexGvW0bo4oao8RT7Um9w6ODy68\/ffo0Mmt5ebkSYNvarmvJGVD\/Aq0onePnacb2RkknyY\/cwE7rbZclZtJ541lB\/0AmlSOcZaoGtv5lmOpOpZVSSctAlxEPF91VOXL48GGnKw3Ajn1T2ggOn8FOuxpVHWxlnynpFFJj160UUaJVEWxTGcW2VyVTKx2bKUU4c4+OjtLc3JzVzW8uGbuON48KkKqUInFg02pQhpq3pBviMrK1niTW19ejm8Ek2xn+qakpOnv2bASya1LJBWw23pfett520jseSUC4ips3OCbbT5w4sdsuU7aYWpp52xzP2txmTdO9cu2+osTBeaBAmgK5ZWxIDgXKVOBf5ZZ6a3IgolAAAAAASUVORK5CYII=","height":98,"width":163}}
%---
%[output:8bf6a4cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8e1419b2]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BaR_A_maap\nTaille de la serie    : 197\nStatistique T_max     : 38.9781\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2007-12-01 00:00:00  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaR_A_maap\nTaille de la serie    : 196\nStatistique T_max     : 2.8679\np-valeur (bootstrap)  : 0.8950\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:9d15429d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:23a2d95e]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BaR_Ah\nTaille de la serie    : 206\nStatistique T_max     : 4.0672\np-valeur (bootstrap)  : 0.7210\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:21e0df38]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAALYAAABuCAYAAACKnelGAAAAAXNSR0IArs4c6QAAHCtJREFUeF7tXQl4VEW2PiEBQhSHRTCSBIKg8uSpKCgoIotjfMMnDorjhiMg4PI9xAUEFDEgIoLiuL0PFR2MozLMTPApfj4dhATDCC4oPhl9onwk7E4kIKAQs7381TlNdeUudbtv973p3Pq+fOnuW+up\/546deqcUyl1dXV1FKSAAklGgZQA2Ek2o8FwBAUCYAdASEoKBMBOymkNBhUAO8BAUlIgAHZSTmswqADYAQaSkgJJA+yKigoaP348DRgwgKZPny4m69tvv6Vx48ZRp06d6KWXXqIOHToQ5\/viiy9Eni5dutDSpUvp1FNPDU\/whg0b6IYbbgh\/N8qTlGhIokElLbAZ1JgrBi7\/NmLEiDD4ly9fTvfddx+9\/vrr4qXA92eeeSYC7AsWLKCVK1c2egGSCAdJN5SkBPbEiRMF9y4vL28EUHBj5t7qbDLw77jjDrr22muTbrKb04CSDti9e\/emQ4cO0caNGyNAbSSqqBNtxK2bExiSaaxJB2yWnTFJEDkeffRRatOmTVi2lmVwO2Crsvb8+fMDTt5E0J+UwIa8XFZWJmRnBqMOx2Ygs7zNcxiIKE0EzVI3kw7YzJEZyLKcjU0gwJufv5QOHGhHpaVEZWVE3boR5eYSZWWFtCiqjB0AOwC2ZxQw4sjMgW+99VahBSkuLqWbblpDO3bcbNhPgLtfv820ceMVYS2JrB4MRBHPptdxw55y7OJiorVrHfeZBg8mGjIkspwRsI8cOUIzZswQqrpFiwrpnnvOEVwaKS1tJ6WnfyT+jxo1ilavzg4\/O\/nkSiIaK55z4pfDeW+DEl5QwFNgz5lDNHu282GjTH6+fjmAuXv3UH5w5aVLG78YeIZ848YR4YVDKioyzqffcpDTKwr4AthOXB1SUkIvgy6wAdahQ0OgHTs2BGq7hPwAN16CbdvscgfP\/UiBpAc2rwpOQcrgdvIS+XGCm2ufkhrYsgjiVKwAxwa4kZysKGZAOryhkP616GrxuP01c6j97x5srphLyLg9BfZNk36mP\/1XBrW\/u4ge+W1Pum1wTnjQ1TV1VFVb24gIGa1SKXc2UdYDdbTK4Llc4B8lLejSS1IoK6eOvt3auC47Cp9xWqoQYVatrqOBg5yVb9miBaWl1stNRFRTsYt2z\/0P6nz7i+L7vxZPoC6z3qXUDlmNurBz507CX1NL2dnZhD+\/JM+A\/f2hShp4zT7a+m4X+t9dh+nev31DBeP+nU5q21pwym+31lF2Ti299t+VNLpda9rZogVl\/lJFGzNa1QvZ9Ww0nyirpk78tmR3BWVkZFBaWpqga3V1NY1ul04bF4byTr63Svzh959\/\/rlRXrPfhs\/qSHuXtaT+L9TQayMrLctP7NKB9rZqSdm1tfTa1i3UZuMy6px3M6V1yiVw64qCeyhr3oeUkt6W9swbTr8acTcdP2BUBA4A6HvvvZc++uiYNsYvQLHrR\/\/+\/emxxx7zDbg9A\/b7X++jsXccoV2rs+lwZTXdXPBPmnhRNv363zpS9\/qNXoNWTnDbXQ2cTxAXTLAB2ExsgDs1NQRqTqJMvdYFebNm1YkXQXDPmupGec1+2zu\/ZaithhXCqjxAzQkvXEpFGfXslEtFRALYP678A5088x2RBcDO6JPXSBxhvTsAkpXVmJvbgcur53gRn3rqqbDu36t+yO16Cuy7plfRP1dmhoE99LQOQhyBZo6B3YhIBsA2JeTLRDROqKSJNLQhjepBWdSBsqjDYcolIihVnAJbPdJHs6AHtJBQ+6PeMQ3\/HXYpLtnNTBHi0phmpb4ENvZsZsAuTamfWIVjm421lIGdS5QbhdquFB2B2g9sVzkQ0qRvGNg6oogRQECHgtCi0Sip7ysf\/e\/evTuc1+pgCYdamzZtomHDhhGMxyDO7du3T9jZwGx3zZo11KdPH1EX5zMadwBsiSpWoogVaJzqsXEwgw2gU60Ia1ScqgmN+q67eTQCSIM0Jbg03q3B9YtQmQR0GdwA9nvvvUeTJk0S3cDJ67x582jMmDERHkLcR7THIH722WfpsssuC+dD2aeffppg275ly5ZwvgDYNizNavPoJrBxkvjyy6ETRIBbJ8knkLqHOnb1yuq+zlP+1mjjiPIqsMGtIZYB1Og6\/nPiZ\/iOxQjPrIANTszubuDi7IwBTj1r1iyaO3cunX322XT77bfT\/v37RTOwjszLy6M9e\/bQl19+aSpDBxxbmX1Z3feXiWeJjaNdcsqxUR9zbR1wA9QF9Ws\/Dmbc4NZ245GfqwBhbm22ReAtAEtmRqIIDLcuvfRSmjNnDuXn5wu\/TzhUIHXr1q0Rx2ZR5Iorrgg4tpPJk\/PGYisybVrIwOn6668Xvooyx8Nn2RJPPmzp1OmneuOm3wjjJ3AodhMDqK65ZhqVly+ko0cHCFA\/91wlvfPONGFEhSTLq7ITgltWfyqw7fau6t5Y5dhMazPZe\/DgwQGwowWvVTkAbsWKfVRSUkJ5eZdRRkYbWrdunXDtGjhwoNi8DBo0SFQh5zn33IP0xz+OERse1iBgIzRlyhS6\/\/77Rf5HHnmEFi1aJDgUkgxufAdwe\/YMtdW162AqKdlOe\/f2EnmPP\/4HKilJo2++eY9KS0uFyStbD1533XWCA1q1FS2t3ODYsozN\/UDfZY7NvxvJ2AHHjnb2bMqB2MuWLaOrrrqKnnzyScFR4drF3BlcFpuaq6++mqZOnSpAB46NcnAkUPMzN+dmrVaJrl1r6fe\/r6FDh6ZGrARcFvXn5uaKJVynLackUoEN9V7Dqb6QsWXFjCxj8zMzjq2uZryiYRyQu7Hi4AXGynT33XeLgyiIIqD59u3bhdx92223RayCViKU03HHI79n6j6zwcjgAcDhs4gEIl944YVhn0PmoDKwrfKr7YGDz5u3jtq1a0e9erWhr756h159dYJhW\/hRXhHA1Zy0pTtxRpswFkdQB2TpbhZaEd123M4XbB5tKIpNDS\/9zLnjBWwnbbHDAsvzdn2LFjhmAJHBLdetqc6Ptjva5QJgW5CKOTXH87ATLYw4tq544KQtmVNztCi7vmkjQsloBRA+qMF\/6LJZrx1tW26WC4BtQk0ADTt0WR622wyqwLbLL8vJum2hjNGmS7ctp+DxI0B0xuDHfnsuY8tqMyYixwOB1oMPFVT7CRXY6gbJyN7CaVt44Z5\/\/vmIuWXVnlyXUVs6gFDzmB2pN3iqGVaJDaV8cBNNu7GWCYAdKwWTvLwRQFhXbTZ01T5LfhnlYJr4HYkDdrIKcNq0afTqq6+Kk0hWjUK7UlBQQKeffro41JGTVRBPt15wN6bZc47txiCSpQ4rYBvdgAVDRxnYKL927dqIaLMA6MyZM4UadPXq1WH1qA6wUY6jaC1ZsoQmT54svuusNF7PSQBsr2dAaj9WYDOnZUDKQ4ORU79+\/YSakrkw9g9WHDsAto\/A0ZS7EiuwMXb5+FwWG9h6Dzp4cHWIHgzshQsXhs0G1H1O0nNsOfiMV4HQ\/bhJcfNFcgJszEfXGTMIVjJ\/aQi8qfaFxQ3I1ThRhNsZuDbk7TPPPFOYuLrJsS+44AI68cQTw4FA3aSN07q0RRGjoOlOG4s1fwDsYxQ0AjZb7fFZANtUHz16VIAYpggANmuUUBs4udnm0ako4jtgc1R\/DFS2eGMyqtdbQB0HWw4EcITFG1RibPkm78pZbYd6cCSO1LZtW2G4BK4PgyUYK8Hjw8zTQ+4btyUbPiHAO4cOlnflRv2QwwnblZHpoI5ftuYze8aMAEZT2LShPStvFoBwwoQJ9OGHH9ILL7wgDK2QZK1Ii4oKOmn8eEr\/4gv6CXHeiOj8E04Ic2yZVrhwOSUlRYgc2Pghde7cWXx+8cUXhQ3I1q1b6bjjjqOLL75YaEAefvhhkW\/IkCEiL7g5xJVffvlF5P\/6668j8CGPHe3B26Zr167+4Nhbtmyp4wijIKZ6jwuDW+XYzD3lyeJoptiBI3Fd2E0D2ByMHc\/ku2FWrVoVcV2GWZsMVgAYwOOXBcfub731VviKDTaqR76OHTuGI6iyYQ\/qNyrD9h\/wFkf\/+EoPeVzoK1\/lgbpleskrCrfbt2\/fRm3J993wWM2Azb6OyPf+ggW0eeVKGt0Qzuq1ceOoX9++9Myjj4oXBzp\/jPm0006L6JfcfzYokwPjy\/2Wy5rNG0ej5StM7rrrLiHmQAyBPQ\/HJI91hY6lfAoDG1zTyq7YDNgq95QDq\/MNAYsXLxZcggGFSTQDhHz6qN4wYAQcJjL37\/HHHxftGB3sqOGA5TLcrszp1VVIXc1U0Yj3IWgf3BdWcfxymIlRKsfnyVT3MWrQTbktAAliBtpg+3KMA3OqPmNg81zguwx86LK57IMPPkgPPfSQ4NYqs2IPHNANJ7mgt+9EEfWEzQjgdsA2iiEdT2DL4JXfbPSdOTM7COA5XkCrMngGVyjZbYpfUtVQnwEuc29wYSNg84tntz+QOfYtt9xCQzkMVcPgYDcO7cUZZ5whTGohHvCpKHwaoa9ev359BJPr0aOHECfefvtt+uqrr8Tn1q1bh8thrEiox6jsnXfeKdSDSMhbWVkZ7sPll18e\/gxRBhzbd8BmajDxjbQedsA2CuPrBrDlG72sOLbZssVyJ8D4wAMPECZLDeyucj91RZHFBQ5LjIlmTsWrVqzAxtKPMR4+fFjIx3KCDMuJn\/Fv+G70nPOrz+RyyKNTVm1D7huenXfeeXTCCScIT3dfiCLr16+vY9kMHEpdlngAdsBGPjsZ26koorapK2N\/9tlnYTkYYbdYFuellfshy+UrVqwIL+Us8wO8LGdyGZa\/8XKoexKrF8+OYzOdrUKcgXtic8nH4qAHZGL087vvvhObP7y8OTk5Edwd5cw4dqtWrcQzo7KIG86rgsqxsWrI\/cnMzBR05j2F0QmlGfOJx+9C3SfvpnnZVj1PdIDN4JaJgUmw4ojvvjudevTYSYWFhSIAe05OZPy3HTtCz5D69etL69alCW9z9E+WT6urs+msswrFZaVwyIVcWFHxBrVtW0jwjOG7HlWZlrmterCBeniHzxycNSlGlzapWhZVNNMFttUky2cJWC2QwCWZQ8rzKPdRNtgCHfAS84vKADQqizzyBl1dyez6Ew\/A6taprcfWrdBpPrhq6ca6Rt1qfoAYv8GHccyY0H9O7HGO8Atmwd6d9jfI3zQo0KSBzaDGiyEDWiU9xwlBPvWKj6YxTUEvnVKgyQKbQa1zQwETBYqGANxOIdI08zdZYCPCkx2nNpoSgFs3IlTTnNKg16CAp8BG\/L4RNx+kNheURgR+56CUmVW19D9VNXRlRh6VUillVuXQ4DnFNKT+1rDitUQz7q9xPIt\/fi1ViC0DBh4rm94yldJbtgjX5fT2gaYarN0x8TQK+CUAvGfARuy+MUs30yllZ9LtU36JCPx+LIzwUMqs2kp7W+44RtI5s+n6PfdQ0aj\/rBeYi+nEyhxasG2F0Ly0SU8nSkmJ+Dw9dyT9kL6LOhzJoof+b5nI84cH21H567+m71tup6xyor+uuolyhk8QQdp1A0hyh5pysHYNnDrO4pcA8J4BG9x65pvf0cUH+4o4eRGB32mo4NCGwYTnzCYq7Ua0FEEJQqlTZQ7V1NZQaotU8V3+XN762EvR4UgXkaf8tociymeXE\/XsNISKqEj79gFum1VpTS1Yu2PEahTwUwB4T4G9ZN1O6rn9TLp\/Vq0A9rHA790bgG1ATQNga9A8MsvQIqIijrEUepRLubSNtmkHaVeB7Sd\/P8f0cKmAG7p6l7rinYwNjm0G7KFhjt14mKVzxhIVD6bcomMc2zExDICNOrwCtlHQSKsXRT4Y4bF75fwh0z4ANkwwLUQRK6DiMAaqvmg0IqiXD23MDoV0L0Jyk2OrMfesfBfRrhyUnT3LVUdexy+7CwUCYBOR1ebRDtg4TcYd7E5OLLlOnEJCK2J2UON08+jGZBoBG78NHz48wocRY2BbdL5tIBZgmzlJGDme4FYDWAkWFxcLc1gjC1A3aOHC+yWq8EzGRuPg2r\/JS6W07P30u74nUY9OGbbjQjBJ6KGj0UfzCaSdHlvn9gFdjs2qS76VwGiAVnfHwAYFhk6w6WCuzIZZslmubCsO81LYj1xyySV00UUXmdIUHjTt27cX8UR4lYBvJDQ9qA+Jo3ThMztiYMUwC0ssG9TZTmYcM3gGbDN3NNlgR+UKcnDIzMwBhEOa+fM3hJ0KrBwl+KQS9iRuHqvbcSlWXfINYmbANrs7hq0T5eDzALbMsdUQcQDp5s2bhc02DMvg0vXJJ5+I8MC4Zg\/PR48eLQzGjFzo2IKRL2liW3a+rwZjUO+swW92tIgjjhtV7QmwQVg5MDt7bMAQHmaXRsHb5WWTN1YrVx6kW27ZTR98EFLzqcHeebQANTi9lQgSLdHtJlO+Ac3s4jJVFJFlaIyVLz2SObYMbPUSJZjjnnPOOcIVD36IuD8GPo0HDhyg1NRU4WxQVVUlLk1CjGxYSjLHhtcPnHs5OE7AsaNFRsObHm2w9\/z8pdSx49N0yiktGgVsl637IH5YGUpF2307YOvUa6UVUc1Jof2A2xlc7dSwZHixcVHSm2++KQz\/wQzApREDHADfsWOHcARAwm3EeCHYhQ5msODmALYcZwR22BzsXpaxrWIj+kH16QnHVifbjWDvZWUpNHZsEaWn96LuuE2pIbkteqh9dwPYOuDXzQMuDcBCtoasDV9F2JXDyVYFNsvROnVjnLIoYlTGT7TwHNhOArDHeouBzgQ6zeOnyXTadyf5A2A7oJaTAOyoNpZg7w665ShrcwG2DlH8RAvPOHYig73rTEq0efw0mdGOwa1yfqKFJ8CWVXpM1HgFe3dr0szq8dNkxnusdvX7iRaeANuOQE3puRuTGW9bEYhwVvGt3aK3G7Rwqy8BsGOkpBuTGW9bkQDYMU5ycyweL2C7aSsCYD\/xxBO0a9cuERtbDs3g5py5QQu3+hNw7Bgp6cZkxttWRL7lDEf08+bNI4RFMwqOGQs53KBFLO3LZQWwdWL3udWgWo9RzD+rjRofEphFrIpXP6PdPLJtORwZ4KFjlIyO1Bl8btiKHDx4MELGNrLzcINuvgL25MmT6ziaJlt5ySF03RiwVR26wNbNF+\/+qvXbTWZ3CnkDsYeOLrDZFsQNWxFEi5I3j80C2N27d68zO9u3Ck3GYWR79+5NiAQK6zMc0eL3SZMmiflT70TkduTwwMjH8bkRiV9dllEHx8jjMGIcJRSgwpExJo0\/4+WU61c5npkFoNou91UFrlw3YgQibPEPP\/wQDiIJW41NmzaF6bHvpX20o0PI7xIeOmbABg3Ymg55ZFohMCcSZONobEUQyL3ZAXvkyJF1HFNZJboOsFEG5WEpBoMa3piA4wDsiBUnB2KHJZkZsDn8L99ZzoHFUYf6AsiiCNqeOnWqaEsOQqnGhuY+qi+yVaxptYwKbICutraWXnnlFRFjEGOWA85zfOxErSRGtiKgJ+yu453sVq94tx8hY8cKbDWGNIfolQGgC2y+O0UNEGkHbAwIgeRxLca5554bXgHUaKh8JI983BbKyj6EqsZAl2PDOhEaB145uE9yIPxETqwXbQXAXr48HOZX5sQMxPLycsF95XDAVhxbjsIP80ysIPLLIC\/xqMfoLhj1NFQ9CTUSo2RRJAC2zxwNYpWx3eTYMHqXXYusZHFVKwJgQhwxCv\/rhGvKDg2Qx636FAA7cl3wFce204qosiyAZ3WdhZEoInNblqPVi5ZQzkiMYO5tx7FVQLKoYdR\/VcZWb2OQJ8jocibuewBsHwPbSI8tL9WyvIvlGapBcED5ch0Ed1fVcerFSKwrx64eu3R4O5vJzhw4ftCgQVRSUhL2zOZrMtAPXNf26aefhi8TAonlzSYfPqgxOHS1IuptaOiT2vcA2D4HthcbjWRp02z5RYC2AilIG74Pqb\/vJV8ZuNWpoxGNsLpAnThs2DBxBR5cvbA554OrNWvWCDcwJM6XKFr7ShQBx07UwJOxHaPJ5LsZyxQgF+NGBqKI80erU0ejI2\/Zk0U9aJGdgKGmtHPlcns+AmA3UFTWRFiFTnB7Atysz2gycZsuwi2AQ6sJz5DGNvy3AjarSVmTw+IfODUOgubOnSsOxeCAu3\/\/flEj9Op5eXm0Z88e4Z2eSMfaANgNbl5TpkwxDLXgJvDiXZfRZCKqYOhIyTjJz41EET5tlYPS8D3p0NIwJ2aOzaIINuZ8FB9wbI9EEQACmz3onPm2WD5xBBwQNgBLOv6npaVRdnW1KVDwHH9GKd5B2WEKiss7cX8kYkMjQdxQZWm5b\/JzlP\/444\/pyiuvjOg+fsfJKUDLCZtmcOjvv\/9eXHD6xhtv0Pnnn08\/\/vij+G3gwIHi9BOXi8IjnfO59XLbBXUPOLYUSwSThwSNB0IEQE0HMM566Ts6cm0mLU\/vRf0zj4aW79JS6nPggHgRIE\/W1NSIADB703vRyD7tGs1fEJTdLUiH6rEL6h4A2wbYy5dvoBuzsoWgevOfSujvf+9CO8f3EMRdWP5XEQAGhld33jmZaHa9ddEYom25IblWTkzoICh77ADXCeoeALsB2GaiCDjtHYWFtLGwkNJ27hSzcvRof0pL20VpaaHvSDs++EA8n\/CPXfTMlFGNZs9PhI4dWt7WoENLnTyJGoVnHjSyVwcGq8bd05GNq7OzBbDNZD8\/ETpRE2rm34jNJ8fpi6YvOrTUyRNN29GU8QzY6Kys7ouHWspPhI5mctwsEwDbTWp6XFdTBzb6v3jx4rBZAdR87HTAen\/ZrY9teNipAPpu2PawKQC0JahDjrCK8MUw+WVzBbMrP3RoqZMnUZDwlGPHe5AqoRF9FX86ySiGNkIR66ZoyiMarBwRFv2HjTfb4hQUFNDMmTOFVgiARqjgzz\/\/PBzyF31jUeTGG28UUVPz8\/Mj1KlGwIZFJBI0UigfbVD3Zg3sWAO+wwBLFWPMTi1VQuP+Gly9p5OMDA1SUnRKhvJEUx59k68fkY\/PZbGNe4FxIzEXlzk2rvkArflFYFHECNg\/\/fQTseEZ6jPi2jqg1cmjT8HYciaUY7sV8N1u48kkSQaOzaeMOhcuwbMdBz3vv\/8+yRybHTEQAxv1MXfmFQHxr\/k3MzjpgFYnT2xw1S+dUGCr3QIhog34bnVqaQZsfbL4I6fMsdEjebXD9+eee07I4LKTM+xJVBkbeWECDJsSti+HVxFOMXFPjSxjI698n40TWgbAbqCWGwHfUZV8ailD0k+E9serEn0vdGipkyf6Hjgr6RnHTkTAdz8R2tm0+C+3Di118iRqZJ4AO1EB3\/1E6ERNaLza0aGlTp549U+tN+HATmTAdz8ROlETGq92dGipkyde\/fMU2EYqq3gGfPcToRM1ofFqR4eWOnni1T9PgZ2oQTnZySe6T\/FuL7AVCVE44aJIvCc20IoYU9gtWxHZoUJtiZ0u4mH34xQ3zQrYOE3XPFE39Fd0cKIeVXkck8g25Vja\/WYrYgcwO2cEu\/JuPW9WwIZLluaJOhm57js4UY+qPPomu5T50VbEzmnDzn3MLeDa1dOsgN0UOTYfqTcFWxE7sCXyebMCdiIJ60Zb8pF6U7AVcWPMbtURANstSsahnqZmKxIHEkRdZQDsqEnXvAr6SUetQ\/n\/B+QooFv3HB5PAAAAAElFTkSuQmCC","height":98,"width":163}}
%---
%[output:030633c6]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAALYAAABuCAYAAACKnelGAAAAAXNSR0IArs4c6QAAF2RJREFUeF7tnQl0FVWax7+EAAkkNFsgkCg4IptOQ4sQYFA8OEqDCKIgggygMCDNqrKpjSKLsjjsyH4EHERHUUFUIjiAh7A0g4jQQIA+MEBATliCCSFIAn3+N9zHfZWqelW33lbv3XtOTpL36u6\/+uq73\/3uVzG3bt26RSqpEYiwEYhRYEfYjKrusBFQYCsQInIEFNgROa2qUwpsCQaOHTtGL774Ip09e5aeeuopmjp1KiUkJHhKwrLl888\/pzfffJOKioooLi6ORo8ezfLgb57OnDlDH3\/8MWVkZNCJEyfYxw0bNqQnnniCunXrRmlpaYatQ17U8f3339ORI0fYdUlJSdS0aVPq0qULtW\/fnipWrCjRu8jIosCWmEczsAH1l19+Sa+\/\/jrduHFDF2rAvnLlSpo+fTq7Ri8B0smTJ1OnTp0oJibGc4mVvLi4Xr16NH\/+fKpfv75ED92fRYEtMYdGYAPqH3\/8kYYPH055eXm6UGuluVn1gHvx4sXUsmVLdhmgnjdvHvuxkpBvzpw5lJycbOXyiLpGgS0xnUZgb9++nYYMGWIINao6fPgw9evXj3Jychj4gwcPpp49e1LNmjWpuLiY9u\/fzyQ1fiNBJZk0aRKVL1+e9u3bx\/Lym6Z79+40YMAAqlOnDsXGxtLly5dp9erVDHz+JHjvvfeoR48eEr10dxYFtsT86YG9d+9eeuWVV+jixYuGOjWqgnowc+ZMVuvbb79Nffr08VI18HlWVha9++671KFDB3rssceYxIW0njZtGi1fvpzlBey4IQC0mPBE+Oijj2jXrl3s+4ceeshL\/5foriuzKLAlpk0LNuAcM2aMZwHYu3dveuutt7wWiqjm6tWrNG7cOPrmm2\/onnvuYZDWrVvXUgsg4SGdDxw4wBaIS5YsoerVq1vKG40XKbAlZl0Eu0aNGkyFgKTmCQs3gKeF9tKlS9S\/f3+mZrRp04ZJ70qVKllqgVjnc889RxMnTqRy5cpZyhuNFymwJWZdhIxnh74M89qVK1fYR9CF33jjDS+pLYKtZyY0a4pY56BBg2js2LFel0NNwUJTm+zWIzEcYZlFgS0xLVqwATUWjYDo5ZdfpuPHjzPIIZHbtm3rqcFfErtr165MB8eCkicFtvdEKrAdgs2hBthlypShVatW0TvvvMNK1Zrb8vPz6dVXX6XNmzfb1rFh8Rg4cCBhkYpNHOjntWrVUmAbzJ8C2yHYTz75JNt55Lt8WOQNHTqU9uzZw0rGjiOkODZZYLGAFJ81a1ap78RmZGdn04QJE6h169bMMsJNgaJVZNiwYYQfcSeTl+FrZ1Siy67LosCWmDJf4Hz77bc0cuRIZqKDqQ7S9YEHHmA1ae3YgBNWlCpVqrDr8T1szzDXcam\/YMEC9r0vO\/a1a9fo0KFDtH79elqzZg0rT+nYEhMcrVl8gQ3AYNb7+uuv2RBBqgPWxMREJrVFPxKzMYQ0XrRoEbVr145dBlAXLlzINmDwt6+E\/Fhk9u3bV1ey+8rv5u+VxJaYPV9go0itdIW68swzz3gA9eUrUrZsWWYLx66hqG4A6E8\/\/ZSZ+4z8TFAJdiOhzjzyyCOlNoAkuuy6LApsiSmzAjYAhOVixYoVrAY92\/apU6eYMxQ89KBXI2HjBp55vXr1MvXugy6\/YcMG9lQ4ePAgk+DVqlVjC1bcDM2bN\/eymkh009VZFNiunj7VeKMR8Ak2FjGQHnpJdLCB7piZmcn0P1gEID1eeukl5gsRzX7BCr3QjIBPsD\/55BNmnoLzuuhMj+ZiWxiPPCTu2Zaenk7PP\/88W53jMQsdD3ZdLJxUUiMQrBEwBRtSGFIZJzTM\/BqweYANCkhpvvpHB7Zt28Z8k\/FZx44dg9UnVY8aAfPDvNwbDaqEmdMN1BWoHKJpCmPLd9qqVq2qnHYUbEEdAVOJzV0lH3\/8cbabZpSgrkCif\/jhh3Tfffd5LuMSHzo3NikAuEpqBIIxAqZgHz16lB1AhaslDo+uXbuWtQnbvNgxS01NZf9jq3fHjh20bNmyUseQ8B124mD2gilLJTUCwRgBU7CxIISKgc0CAN6qVSvm6I5FIfyIsQuGw6KAF+qInlTGZgIsJVppHozOqTqidwRMwYbxHztmABcWEJ5wdAlSHE46OKI0e\/ZsR2DjaYANhpSUlOidCdVzNgLYZdVz7LI7PD7NfXoFAkIADbs1JDFiYxhJbF+qCKCGBxxcMfF0EBPqKSgooAoVKkh11kl+J3nRByf5neR1c91oe+XKlZnDl9Pk09zHwRJjW3C9GhIdYMNHWHbxyDeA4MeMg6digjPR+fPnmSSPj4+33Vcn+Z3kRUOd5HeS1811B0ViX79+nR1twmkQ7aIQ340fP559h7N9+C1r7uNgQ+rz+Bmc4MLCQjp37hxzqJcB20l+J3nRfif5neR1c922JZdJBlOJDSmKGBdTpkxh8S241Oa7jDg1DTNgbm4u26ABgOIuo5UNGgV26dlRYDtH3BRsBGZBqK5NmzYxq0iLFi3YCWtI10aNGtH7779POKWN9N133zHneiwyEcgFUtzKlnq4g40T6BcuXLA90tCTkQ86o92njZO8XMcO97oRl9AsNqHtAddk8Ll4xO4jJDdghmslbNdwitI6N+k5Qb3wwgvMEcosxEA4gw2f6rlz59Lu3budjrPKrxkB+BTNmDEjYHD7BDvQMxJKsKdlnKBTlwrp7qrxNLa99+YR1IGNGzeyw7eYAL4ZFejxiIbyISgQU1BvXeWv\/kct2Gv2\/EpD1hz2jOPY9nW94BbBDuQE+Gsi3VSOmTDzVz+iFmxADbh56tk8hRb0bOT5H9L86E+ZtGXpWwGVLP6aSDeVo8AOoLlPK7EBNeBG4t\/FXciixO3TFdh+vmsU2AEEG3MFqbz9eC61qVfZSw3h0lyB7WeibxenwA4w2EbTFk1ga0NF8DEJ5LpCgR0CsEUVJdIlNo8liB1fMcglB08v+KU\/ZLgCOwRgi4vKSAcbLsU4JGLkbmz0nVO4FdhRCjaAw44vT6Lk5DFNmjVr5nlbGVcn4IwGpzQk7BTD12fp0qWe135oJbA\/feXFNteuXZvVjbgq2J3W+gApsAMM9rWs2XSzIJtiK6RSQoORDAgsKKdlnGR\/h0Jiw82Xe03imJ2eusDB4OEvOFRcLxYD+vDP9MoRr3Pyrhptm\/XqF6W8AjuAYF8\/\/Tld3TfaM94JDUYwuJ9asI8y\/5ErBbbZTqaVxzcHAsfuxBci6YHATy1BNx41ahQLPsn1ZLNycK14mkkvbkyTJk0sn1E1qkt7symwhREQvdzmbDtnuP1tBA3yT\/jqMF2+EUf\/UiPRy6QHqAE3T+Xv6kbri15ju5EDr6ykWsW\/0v68eNq67xdLdmxfO5lWwDZSDbi0RbwWDrwYRF4LIodNqwYYQcjbZqYC4UbCq0e0byDDjaG9WVCemVRWEvu2xN5ymmjUVyXqAZJ2+9sIGjPYju3\/b6r2\/+M9WX9OepPW\/vpvlLd1JU1NW0yxiXG09zjRsO+SLIHtayfTKtiibq3No1UXOIjazwEwXvaEF6SKUQP0bhC9dokmQK7KGIFtdDMa3Vy+oLcyTlau8duWumyIMytOUJO35tKXB0rUAyTt9rdRR81g+98ty6hp3hRP1pzkQfRTXD\/K3TiUejTZyz4H2IPmx1gC22wn08pE4BozK4W2DFGPxWJNVC+sSGyE1MC5VfEpINahle4AG\/FlfvjhB7YY5QvRiJfYsiHOrICtldji9rcZNGaw5X7fnW4W\/p8nO1dF8veNomdTMm2DjQxGO5lWwfYFJFctuERFtFYcyYNkRuLvdLeqY5tFF9C2RVwgIuIXbgro9PgbFhjtuiAidGwnIc6sgI2TOdCx9ba\/zaDhOvbfL5LXtvnl\/3mHrmfPp3L33oknGBv\/EP3HsTfoL0mvU3rlLCmwrQJsdp3WwsAhFsE1soJwuDiUyMMlOf9MXGRy1QTXibZsPQsK2oXEF6how8mTJ9n\/RlYRsX6xz67RsZ2EODPrJCTgoVMX2AkU7eLPCkR6R6zytq6gKxuHUYXW3i\/\/LDpfSDvP3EsPNyuJU21XFbHSHqvXmC3ijHYFRX37wQcf9AQ6ArBnz55lVRuZ9PTeOKbdUtfq2CLYXI3i6wNux4YNXavn41qjOXdqVRLH1y86tpMQZ+jMrKWrKT29JaWlpTKnfyQcABDdSvGZVd2ad1ALNiQ1wI7\/UyHF1fR96t2Ojm0V2mBcZ7Zwk61fK7G1\/2vLBbz8ZVDa0HYc7CZ936Xa9Zt4sorzbdVAYNQfv4AtG+JMqwMbNZKb386VSaG8fx\/r5TdtNlGiue8vBaspafM0pn5opbVRGQrsOyMj6uP41EzHNvJB4aVxsPPbjKGi6g10h9+uENMW4jew0Vi7Ic5gtUjOWUyp8eaHZZ88\/w0V\/lJiFVlSqS9Vee7tUke59EaH3zidrmbQhEslOmLiEymWpHUoVRFZqcrzBUpii1YRUa3R2+Qx28nk17\/UuwHVTKtWqrvZhdUJVirx4IfdMQkp2FqTm1njoQPfzC+iPUX304aUeTSn+52orkb5Rnx2jGpdXk6dY36g2kUlp2XEBaOvwXKrxPbVr1B\/z8FePPQWNaun3xqA3aDVOOmmBhxssxBn2t0\/q704kfACVa3X3+flvxxcRw\/fmuPzukhTRaQ7HKSMVsCG+bXin2ZIt8gvYMsuHmXBLlOrC5X\/16k+O339wDgqPrfO53UKbOkhksroGrBlzX0wuxX\/vsH24Fi9m2VvHN6gPT8X0eAVZS3tPNruRBRnsAJ2mXKd6A9\/nic9Sn6R2HyDxm6Is5wFL9KNS19RXEpp01tsxTjmr6GXggE2dPrMtbn02tFaCmxpvPQz6oGN8eap6NdCSnp4FiU92k+6Zr+AjdplQpzBpgy49ZKZWS4YYP\/+j3za+bfK9NruYgW2NF7mYM\/rkEfNm8YRQC46V5niatRlGRIaP8osX06S38CWDXF2YNZ\/sk2TlHLe7wY3A9vqY8qJKvI73U+7bg5gkaACebDVyeS5NS+X2P9V\/xw1SSqR1JDOyUNKTv\/4I\/kNbNnG8E6uG\/YYVbtynN2tSNC9YyuVtm\/jkRXz+58tDYKRDg+z4c2rRUxSIOEm0qo9ZRpOoi2HEhXYshNrkk8PbEDtRPXQVhc2YGulovaEC294wY4LBFXEyt3NdXgtuCgDqoaY4v9Y2aPrA\/iYOoNod6WWEQ92KMMvLBnYnhoXHvGL6uEasNHQS190paKLuzxtBnDYgbR6dwNsqDlIHFxehi9BFNeyB\/107\/MRDbYKv+CLAgffm3n3iWDyKuLvf5RqT9hiqUa9\/NqMKK\/w71tLlRcNYKvwC5YwkrvILthWpTVaY2R1iUuuy1bgfPXNvf6Kcu4cPwsl2G4Kv8BVGRxB46dr4PPNDz3oUeEaf2w5pEtyGXXSCMoq3SdYNgXBu+\/sqnEU9\/M6EqHVW4Fr60vqv5AyC1OCroq4LfyC3vlIXzxENdhGaoQdsxD3x07MyqC85YM9420k9SG5rx3ayiR5QuexnsDvwTL3uTH8gt4JHwW2hMS2o4qIBw2urZ\/mgdaK8V8m8DtuDDwZoOpYqUMLgBvDL3CwcWpGjP9nBndUS2wMTM7Hf6W8\/ZspPiGBjZPdHSknb9+yC7ZWlbGjMnEItLq1Fo5wDL+gwDa4ha0e5rX75i1UF0ywtaqTHZVJBNtqIMhwCb+gwI5wsLUS247KxIfGjeEXEIVq3LhxpFQRDeCRIrHRLXHxKaNjowy3hV8YPny4AtuuTdOJKhFsVcSXJcDO924Kv6DAjnBVxA64gb42UId50W5u+fAVfiHsrSJ6J5B5o8VVeiBj9+EAQ7gvHgMNq53yAwU2j0SAtvDwC9qg7lbaGRbmPqzSZ82aRV26dKGE22Y33ni8N7158+bs30DG7lNgW8HlzjWBAtso\/IK91pmHGLZbltH1pm6rkMKQykeOHGHBD43eiR6M2H1KYvtrykNfTsgl9tWrV9mKF3fqxIkTqVy5crqjInuYF4VFklUk9Mi4owUhBzsnJ4cGDBhAiKU8dOhQw1GTDb+gwHYHiP5uZcjBPnr0qCdq55kzZ2jt2rWsjx06dGDxkFNTU9n\/srH7FNj+RsYd5YUcbCwI+\/TpQ2XLlmWAt2rVig4cOEArV65k+vbChQupfv36DGyZ2H0i2CNGjKD09HSvmSkqKmI+vg0bNqQ6derYnjUn+ZF39erVlJGRoQ7z2h558wwcbL05T0tLI\/w4TaaLR7yWDQ7jABcWEJ6ysrKYuad169Y0adIkmj17tjTYeBKMHj2adu\/e7bQvAcsfLLfVgHUgzAo2MyEDdvw4TTE7d+681atXL69yfL1qGNIMQGdmZrKI+Zh4I4ltFruPVwq48RNuCTfbnDlzlMT288RwsGfMmOFRZ3kVfpPYZmDD3FdQUEAVKlSgmJgYr+6JPg14IyzMgeILfnAxNxfu2bPH8nsD\/TyGjooLhi7oqIEuzRyMcTVURa5fv85eG3z8+HFatmwZJScne4YR340fP559t2TJEvYbuviiRYuoXbt2nuvy8\/PZ0SpEtDczF4br\/ARjAsK174FsVzDG1VTHXrVqFU2ePJmmTJlC3bp180htvssIUyDMgLm5uTRkyBCyG7svkIPnj7KDMQH+aKeTMkIZVySQaxdTsPPy8tjL6jdt2sSsIi1atGDv+EODGjVqxF4AX6NGDTauMrH7nExIMPJGOthRHVcEu4+Q3IA5OzubKftYbEL1wI4kT7Kx+4IBqGwdkQ62iisiS4bL84UKbBVXxDk4IY\/d57wLgSshFGCruCL+mU8Ftsk4BhtsFVfEP1CjlICBrdXN8c7tp59+mgYOHOhlOpQ9oADLDGznH3zwATMnikm2TJShbTc+s7p6v5Y1m24WZFNshVRKaDDS9iyFc1wRo85E1Sl1eAViWxQSCAtNvAKZ+5jcfffdNHfuXM+Ok90DCoAW\/iOw1mCXSu\/dknbL5JOmbTc2pqAaWAFbG\/Y4ocEI23CHc1wRBTYRwccEmzuQqG3btvWMycGDB5mPSe\/evb3s35DmONCQmJjIrt22bRvhkCg+69ixoyc\/tt3hePXZZ58RtvVx9F8LtpNDD9p221FFtG9PsPo6EREYMyuFFqxgxxUR7d18TuF6ETXhF27evMkkMvxIFi9e7KUm8IMLsbGxzLkKNnGrO5bc5oo8Xbt2JdRz8uTJUmDLHnrQa7cdsLUSG+8oBNx2UjjHFcFNt2PHDs+84UmMJ1nUgG02kb\/99huT1FWqVGEDtG7dOss+JgAbTjM4e4mNIvyt53jl5NCDtu12wEZe6Ng3LuymstXTbashvO5wjCvCpXXPnj0Jh3fF\/xXYt1UMgA03VUhqJwcUjHzAnZTpFGw70tns2nCLK8Kflgi9oMDWzBxO4wwePJhZRKCqYCveyQEFM7BlDz2EC9hObxB\/n1I3k9iuDb\/gdJCR\/\/Tp08zLD3c+P3WDz8MdbCwmYd2xYhXxxzj5qwx\/g412GenYUQv24cOHGRzFxcVMn4bzlKhPyh5QkJHYVg498Lah3bC5wz9GgU0evRo3OxaOcJDjOrfdG9Lu2sVu+bg+YBs0sDfDngxJXbduXZo5cybdddddXm10stAzAttJmWic2G6YIaFCuQ1sGRDs5METYcyYMTR9+nTCu2fsJleDLW6SwJ9bPKjAB0LWNGemxjgpE+WK7X722WfZukCBXaI2wnzLk5MxcS3YfKHYuHFjdnIG5j29xDdTZA4oGElsJ2Vq241Dy9xeK6NL2pVk0XK9K8HGsTHAvGbNGsJr0aCGaBNUks6dO7PIUrIHFMwWnjJl6rUb+vUXX3yhJLaf7zhXgs2jR8E3xCiJ7wGUPaBgBrZMmWbtdvLY9TMTEVGcK8GOiJG\/3QmzwC6R1M9g9wVPQmzSBVJgBMwqEuzBCkR9bgjmE4h+B6NMRP2CW4Q\/oj7ptVeB7WMWwzWYTzDgC2Qd\/gqMY9RGBXYgZ0+VHbIR+Ce3QD89G3TquwAAAABJRU5ErkJggg==","height":98,"width":163}}
%---
%[output:165b28c0]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAALYAAABuCAYAAACKnelGAAAAAXNSR0IArs4c6QAADmFJREFUeF7tnQmMTfcXx88MqmNXY61lRI0lCLGLjIREtbZqqK1hxi6WWqr2XSwVbVpkENS+dYLYmjZICEIQEWIXsdaMvdQSZP75\/prf\/O977pt3372\/ufd33zu\/5IXMu7\/tnM8977eeE5ednZ1NnFgCUSaBOAY7yjTK3RESYLAZhKiUAIMdlWrlTjHYNhi4evUqpaWl0b1796hjx440f\/58SkhIyCkJ05aMjAyaPHkyvXv3jvLnz0\/jxo0TefB\/me7cuUObNm2iP\/\/8k27cuCH+XLNmTWrbti117dqVKlasGLJ1yIs6\/vrrL7p06ZJ4rmjRolS\/fn3q3Lkzff7551S4cGEbvYuOLAy2DT3mBjag3rFjB02cOJHevn1rCjVgX7t2Lf3444\/iGbMESOfMmUMdOnSguLi4nEes5MXDn332GS1ZsoSSk5Nt9ND\/WRhsGzoMBTagPnz4MI0cOZKeP39uCnWwNc+tesC9fPlyatasmXgMUC9evFh8rCTk++WXX6h06dJWHo+qZxhsG+oMBfaRI0do2LBhIaFGVRcvXqTU1FR68OCBAH\/o0KHUs2dPKlu2LL1\/\/57Onj0rLDX+RcKQZPbs2VSwYEE6c+aMyCtfmm7dutGAAQOoSpUqFB8fT0+ePKGNGzcK8OUvwbx586h79+42eunvLAy2Df2ZgX369GkaPXo0PXr0KOSYGlVhePDTTz+JWqdPn059+vQJGGrg75cvX6a5c+fSF198QW3atBEWF9Z6wYIFtGrVKpEXsOOFANDGhF+E9evX0\/Hjx8X3jRo1Chj\/2+iuL7Mw2DbUFgw24Pzhhx9yJoDffvstTZs2LWCiiGr+\/fdfmjBhAu3du5eqVq0qIE1KSrLUAlh4WOdz586JCeKKFSsoMTHRUt5YfIjBtqF1I9hlypQRQwhYapkwcQN4wdA+fvyY+vfvL4YZLVu2FNa7WLFillpgrPObb76hWbNm0UcffWQpbyw+xGDb0LoRMpkd42Usrz179kz8CWPhSZMmBVhtI9hmy4S5NcVY5+DBg2n8+PEBj2OYgolmcIq0Hhvi0DILg21DLcFgA2pMGgHRkCFD6Nq1awJyWORWrVrl1KDKYnfp0kWMwTGhlInBDlQkg+0QbAk1wM6XLx+tW7eOZs6cKUoNXm578eIFjRkzhvbv3x\/xGBsrHoMGDSJMUrGJg\/F5+fLlGewQ+mOwHYLdvn17sfMod\/kwyRs+fDidPHlSlIwdR1hxbLJgxQJW\/Oeff\/7gO2Mz7t69SzNmzKAWLVqIlRG5FGhcFRkxYgThY9zJlGWE2xm10WXfZWGwbagsHDj79u2jUaNGiSU6LNXButapU0fUFLyODTixilKyZEnxPL7H2jOW66TVX7p0qfg+3Dr2q1ev6MKFC7Rr1y7avHmzKI\/H2DYUHKtZwoENwLCst3v3biEiWHXAWqRIEWG1jedIcpMhrPGyZcuodevW4jGAmp6eLjZg8P9wCfkxyezbt6+pZQ+X38\/fs8W2ob1wYKPIYOuK4crXX3+dA2i4syIFChQQa+HYNTQONwD01q1bxXJfqHMmqAS7kRjOpKSkfLABZKPLvsvCYNtQmRWwASBWLtasWSNqMFvbvnXrljgMhRN6GFcjYeMGJ\/N69eqV6+k+jOX37NkjfhXOnz8vLHipUqXEhBUvQ+PGjQNWTWx009dZGGxfq48bH0oCysHGkhbOQDRp0iQmD98wanpIQCnYWGvFyTScR47VU2V6qJVboQTsN2\/e0KFDhwTMN2\/eFFJlsBkuLyWgBGzM0nFjBJMXrMlu27ZNHLaPxXPAXiqT6\/6\/BJSAvX37dnH\/DzN5nHLD3T5sPDDYjJpXElACtrHxcimMwfZKpVwvJMBgMwdRKQHPwYYbAWwulCtXLioFzJ2KTALYZTU72BVZKR5bbECN0284honrVbGa8GK\/fPmSChUqpESpfpUj5ICjBPCn4hRuTy02TrBhwokzzLh0GqsJh6YyMzPFr9bHH38cq2IgyOHhw4dUqVIlx3LQAmx4Q5K+M2JRq69fv6a\/\/\/5bXByIZbBVyoHB1uBNUqlQDbpjuwkq5cBgExHG+vh4lTC2xE9wiRIlYsZiYxwd7JuQwVZIoJzAnjhxQmGpXFQ4CTRt2pQWLlwYALfWYIfrkPF7OXn0cowt2wAhf\/rpp5E0n5+1KQEYEfgUDNY7g21ToGbZdHi5FHbHF0WFkjmDrVB9DLZCYVosisG2KCgnjzHYTqRnLy+DbU9uEeVisCMSl5KHGWwlYsy9EAbbBSEHVcFguyBzBtsFITPY+gjZ\/ZbETo1ssV3QNVtsF4TMFlsfIbvfktipkS22C7pmi+2CkNli6yNk91sSOzWyxXZB12yxXRCylxYbXpxWr14tbrPgJgMO\/o8dO5bq1asX1lunWUwW2RezeCnyOx2g0qEN7qPlbY2uWWzpbw9RZRFQs0KFCvT777\/TqVOnCE7HEeEqt4TAnciH2ChwUG5MtWrVonbt2plm1wEqHdrgLWbu1+4a2HBFiwhXxmBAgB3enRAFC0cMg4E1imPLli20cuVK+u2338R9NatJB6h0aINVeUXLc66A\/fLly2x438fwAxFj4XVfpoMHD9J3330nhicNGjQwlSs89MNP36VLlyKKW4jCdIBKhzZEC7BW++EK2FlZWdmI+IpAPsGxA69cuSLcleXmh09Gm0VwoUiDauoAlQ5tsApEtDznCthXrlzJBryIdNWjR48A2clJIcbOmEiaJRkKuXnz5lS8eHExAf3nn39EDHCUWbt27ZD60AEqHdoQLcBa7YerYJv52pMBN7FCEmzNZSekVQfgnTp1oi+\/\/FK4EsCYG\/lzm3zqAJUObbAKRLQ85wuwEUQIcQy\/\/\/576tq1a87S4P3794XHVXj2keHcghUjO4gxeufOnT3RG1Z+4IUKcwlcMOWU9xJAvB14AMOvu9GfDK6GwUAqcZgjhyJmFtvKUCQ3MaDhCNYZavIpwUa4Nq9cnOEF\/PXXX4lvqec90MYaYEQwdzP6bEQAAXxq1Kjh2A1FnNPJIxqLCWRCQgLFx8cHSEc6hA91C12C7bWLM8CNj1cJlurp06eUmJjo2GedV32ItF4AHeyIVKmLM6fLfYsWLSI4fscadnJyck7\/ZHhlwI3vqlev\/kHfeXz7n0hU3s6OFDCdnlcpB+EJChs0U6ZMCZjoWd2gQeyZgQMHip1HY2xvTCrxt4YNG9Ls2bNNYw4y2Ay28cVSDjbOiWAwf\/bsWRGeGME2zbbUZSjl06dP51hhjIkwjsZ4Wq6K3LhxgzZs2CAmjtjNNFpyY0cYbAY7T8FG4cGHoGBpx4wZI9z7xsXFifrNwMbfATesPuDGDiSCLH311Vc0aNAgKl26dMhfOwabwc5zsL0YazHYDDaD7cWb51KdKseWLjU5T6pRKQflboQj6TFbbLbYbLEjeWN89qxKS+Wzrgc0V6Uc2GJrQIJKhWrQHdtNUCkHBtu2GtRlVKlQda1yvySVcmCw3dffBzWqVKgG3bHdBJVyYLBtq0FdRpUKVdcq90tSKQcG2339scUOIXMGWwMYVTZBpUJVtsvtslTKgS2229ozqU+lQjXoju0mqJQDg21bDeoyqlSoula5X5JKOTDY7uuPx9h+GWPjUsHRo0dp8eLFdPLkSXG6r1+\/fuK6F9wyhEq8pf6fZFRaKg3eU9tNUCkHJRYbLs6GDRsmLsPChcOFCxdo7dq1lJKSQjNnzgxwwmPsNYPNYBt50ApsnOMG1LDSuG0uPUnhZg0ua+JvcMlglhhsBltbsAEnhhzLli2j1q1b57QTV8twUeGTTz4J6SGKwWawtQUbDilx\/Sv4wq706Ycx96pVqwTgwYnBZrC1BXvBggV07Ngx4fkp+BoYvtu3bx+tWbOGqlatymC7sBpge+amQUatxtiAF5bXzCrD9QJWSsK5X4h1L0zv3r2jAwcOUM2aNalKlSoaIOZNE6QccKE8KSnJUSMcr4o4AfvOnTvidjx7YXKkw6jKjJW1hQsXUsWKFR31K0\/BDjcUQcsBNz6cWAKQAIB2CjXKcQy2k8kjq5IlkFcScAy2k+W+vOoUl8sScAy23KApX758wC6jlQ0aFj9LIK8k4BhsNOyPP\/6gUaNGiehi3bp1o2vXrlnaUs+rTnG5LAElYJsdgurdu7c4CFWsWDGWMkvAdQkoAdv1VnOFLIEwEmCwGZGolACDHZVq5U55AjaC62CrfceOHUIDbdu2FTuQlStXjiqNeBGf3u8CxHwtIyNDnD+aO3euCAFjJ7kOdlZWljinnZmZKf6F721cSkBsyPT09JBO4u10zss8XsWn97LPTuvGWRGcL0Ig3Hbt2tH8+fP9AbYxLg0OTSE6FJIMnYfIB9OmTTMN6+FUaG7n9yo+vdv9VFEfuEBoF8Qz2r9\/vyiyY8eO\/gH74cOHIspB3bp1aerUqQERshANAdHFcBIQcf78nGTkBy\/i0\/tRbjLsIn7NcbUQQ1XcxPKNxc4tNjvuTeImTqjQeX5SmAzD7UV8ej\/JSbb1+vXr4jx\/amqqOAA1adIk8ZVvwDbGdcQupTHJ7\/BzhNjtfk7SAnkRn97PckPb5a+dL8E2s8oSBrMIwX5TVm59yev49H6TVXB7GWyNNegUbCfx6TUWi6WmRR3Y0TgU8SI+vSV6NH7Il2Dz5JHEslZaWppYw+\/evXtIxOzGp9eYWUtN8yXYvNxHdPDgQcLl5XXr1lGDBg1Mle0kPr0lejR+yJdgyw2anTt3BuwyRusGjRfx6TVm1lLTfAk2enb79m0aPHiwWNYZMGCAcFoZjVvqXsWnt0SPxg\/5FmzINPgQVJs2bWjs2LFUrVo1jUUeedO8iE8feSv1yuFrsPUSJbcmGiXg+um+aBQi90k\/CTDY+umEW6RAAgy2AiFyEfpJgMHWTyfcIgUS+B\/9unt6N2iJWQAAAABJRU5ErkJggg==","height":98,"width":163}}
%---
%[output:59fd8b89]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAALYAAABuCAYAAACKnelGAAAAAXNSR0IArs4c6QAAF+ZJREFUeF7tnQvQVdP7x1eky0QpSq4VMYkhya1JGpqaFF1+imhI+iPlEipESheUy4hICPkTiaYoEzFpyi2NP5Eiyq0idB+9iv7zWX7rtc56176efc4+57xrzTReZ6\/rs7772c961nOpsmfPnj3CFUeBEqNAFQfsEttRtxxJAQdsB4SSpIADdkluq1uUA3YMDHz99dfi8ssvF+vWrRPnnXeeuOeee0TNmjXLe+LYMnPmTDF8+HCxe\/duUbVqVTFkyBDZhr9V+fHHH8ULL7wg5s+fL9asWSN\/btasmejQoYO44IILxGGHHeY5O9oyxptvvilWrlwp6+23336iRYsWomvXrqJjx46iVq1aMVZXGk0csGPsox+wAfWsWbPErbfeKnbt2mUFNWB\/9tlnxfjx42UdWwGkY8aMEV26dBFVqlQprxKmLZWbNm0qHnnkEXHMMcfEWGHxN3HAjrGHXsAG1IsWLRLXXXed2LZtmxXUJjf3Gx5wP\/744+L000+X1QD1ww8\/LP+FKbR76KGHRP369cNUL6k6DtgxttML2IsXLxYDBw70BDVDffnll6Jv375i48aNEvgDBgwQvXv3FgcddJD466+\/xKeffio5Nf+lIJKMHj1aVK9eXXzyySeyrXppevbsKfr37y8aNWok9tprL7Fp0ybx\/PPPS+CrL8Hdd98tLrzwwhirLO4mDtgx9s8G7GXLlonBgweL3377zVOmZijEgwceeECOeuedd4pLL700Q9Tg91WrVolx48aJTp06iXPOOUdyXLj1vffeK5566inZFrDzQgBovfBFeO6558QHH3wgn7dq1SpD\/o+x3KJs4oAdY9tMYAPOoUOHlh8A+\/TpI0aMGJFxUGSYHTt2iFtuuUXMnTtXNGnSRIK0cePGoWYAh4c7L1++XB4Qp0yZIg488MBQbStjJQfsGLuuA7tBgwZShIBTq8LBDeCZoP3999\/FFVdcIcWMNm3aSO5du3btUDPQx+zVq5e46667RLVq1UK1rYyVHLBj7LoOMtUceRn12pYtW+RPyMK33XZbBtfWgW1TE\/pNRR\/zqquuEsOGDcuojpjCQdMsUceJQY6CbOKAHWNbTGADag6NgOjqq68Wq1evliCHI5911lnlIyTFsbt37y5lcA6UqjhgZ26kA3aWwFagBth77723mDZtmhg1apTs1VS3bd++Xdx4441iwYIFkWVsNB5XXnml4JDKJQ7y+cEHH+yA7bF\/DthZArtz587y5lHd8nHIGzRokFi6dKnsmRtHuDiXLGgs4OIPPvhghWf6NH766ScxcuRI0bp1a6kZUapAXSty7bXXCv7pN5mqj6Cb0RhLLromDtgxtiwIOPPmzRM33HCDVNGhqoO7Hn\/88XIkU48NONGi1K1bV9bnObpn1HWK60+aNEk+D9Jj\/\/HHH2LFihVizpw5Yvr06bI\/J2PH2ODK2iQI2AAMtd5rr70mSQRXB6z77ruv5Nq6HYkfDeHGkydPFmeffbasBlAfe+wxeQHD30GF9hwyL7vsMitnD2pfzM8dx46xe0HApkuTuyKu9OjRoxygQbYi++yzj9SFc2uoixsA+qWXXpLqPi87EwbhNhJxpm3bthUugGIsueiaOGDH2LIwwAaAaC6eeeYZOYJNt\/39999LYygs9JCrKVzcYJl38cUX+1r3Icu\/\/vrr8qvw+eefSw5+wAEHyAMrL8Mpp5ySoTWJscyibuKAXdTb5ybvRQEHbIeNkqSAA3ZJbqtbVGLA5rS\/ZMkSeWJHh4u8169fP2m9Vpk9ORzE0qFAYsBWtsinnXaauOiii6Q+lYMRp3Ju4lB1ueIokC8KJAJsrnu5UoZLK30tC3j33XelNwm\/nXvuuflakxvHUSCZ8AvckiFy6JcJ0FbZRtSrV8+ZWTqw5ZUCiXDsF198UdpAPP300+Loo48uXwByN9wamZtrZQBeqmXTjFFi98a1omr9xqJurztLdZlFs65EgI1xznvvvSeefPLJCo6jPMN2gosKLh\/MogARlmJRgRO1\/7DzMOttW\/jPRQylxnHtxD71w3nGxB2v2Nqtq3qQmFL7slDTPqJeDTGsY0WshGr830qJARtxxMaVuf5FU2Jyc8YHDBsnXR5lvrJu3Z4jQ3HFuP1HnpBrEIoCAHtKnXDgHtaxcVbgThXYgFrndKGoI4So2eYSUfeqJwOrb3q8v\/hj8fOB9VyF\/FDgtVodxah6mZ4\/XiP3PqWhmNT72NgTyzmw\/USRuBy1Sq8JQrT6T\/CiP35F7JkxJLieq5EXCoysN0y8XqtjqLEANeCOWxIBdjaHx+UP\/o9YPf9\/xWn\/DQqjLwTDnt27dokaWviwms3biZrnh3vrMR\/d9PIoUWPdPzE6Mvr+Za087KmCXEzZbfzOb+pZFCLb5h62vd4WWX2XNs8wfSQ1to0eXrQw6bbx4NPF4fVqlE93WfUTQ8vYbZrun5UYwqCJADsbdR9tsWQjhp2KeKSosXPnTrF+\/XrpAlWjxr9EIm4d\/4IKG\/zrr7+K\/fffP6O9agfod36xUAK3bs9\/NBnqN1WnaoPGov41TwcNVfHFCRjbr8OgeQdNJpv2trY2OplzWPLNZsFB\/eSy\/xPLqrcQ29oP8+W4xCX0i00YtMag54kAW13QAED9ljHMBU1UYANo3K0+\/PDDoLW55wVMAW6oJ0yYkDNwJwJs6PfGG29IdyjiZRB6C0\/tMFfqUYGt6kOUQw89tIC3zk3NiwIwJWIK2r7SSVEtMWDbjKAuueQSaQjlFxQmLrBzSZSkiOv6sVPAb8+TolkgsAnLRUgBgISXB\/Yg3bp1k6EA9Cieca37HLCT2sri6Sd1YON+dP311wtcoTjgtWzZUsaOQ8Q44ogjxMSJE8vFgbjWfQ7YxQPIpGaaOrDxpyNMlxnRCB87YtARNoAYGps3b45t3eeAnRRciqefVIH9999\/S46M8wAx4XQDJhU1lBC2eF8TZDGudV+hAdsMnaDg4iXTYzJA9gK92GLrqed6mDNbn7a4gHrfJ554YtEblKUKbL\/3f+vWrZJTE8QFYM+ePTu2dV8hAVuBDn26HvRRzdEELLeqpo2M6gP6+dnOcD5BnDPz1yhgE0jHDNiu+ra1Kx5+LSTNvO4uklpH4OHRNhD6aYCNPhlOnY11XyEBG+7LLaoXIPVnCoD33XdfhYslr2fqa3DIIYfIEMM24zA\/YLMX+QBFUuDy6icfa4gM7K+++kqml4DjIKoQH9rGudSi\/Kz79I3CbptsV3rh5nHDhg0ydp3KyvXxxx\/LlymMuu\/e+WvE97\/vFGHNIIPmqs8tCIC2TdU3lKRHnFNwo9M5c1C\/qo9iTsGh1jB16tQKTIHgQLZ4hFFftkjA\/uGHH2S0UD6JhNpSGamSADZhuACsXsrKymSuFl4iFTIXeZ45BAF7+tINYuD0L8u7C2MGqcu3QcDRZfGw8fFMOtnoFgRsJdMHrT8qEPJZXwGblCWcGfSC+QMibralyvvvv78HeSfo8EOwRFR\/RO9HS3Lssf+aFPoBO8jRQC0SXTn5UvQCeH7++WfRsGHDclsP6nPpE7SxgBpwqxLWDFLNJ8yBTT8I6vVtc1N1dQ5t+yT7AVs9O\/nkkyvI5tkCIZ\/t\/fY8MY4dBGwuXtBRwyWRC3nLDj\/88Aw6ZGPdlysZ2+TYccwgTY2Hn7bD1GaY2gubmGM7DAZpRfzmkE9wZjNWQcjY+sXL2LFjrTkDmWghqvuQsRev3iyyNYPUxY6gLwUbrl4IBUIvFaICB4dJ5WFk49hBmpZsQJZG29SBrQ6KzZs3l17mXrJPPq37ckUUm6jgd1hExCLltKmuo40CMn8rPb+XessEspcoUkrgztUe6vvleXjk4AaYCSDO4ciWtg2R5Pzzz5fZq\/Jt3ReGc0blRn5nBVOF57c5uloPfbhfv6oukVdRM5J9jJzrNj22lz496jrTrp8qsPW8gl6E0LUB3FSSAg5tCZljKWeccYa4\/\/775eHPq+RKxo6zeV5c0eviRiU0Ml8yHcjMA7Weeemjz0+nAUZmXsCmjRozSGsTZ\/35apMqsG2LVMHMyfFtep2XkhGULQOX1xfCpkXRD3hhdOP6IZIkpdwT2Di2LuZgx5OLr1Y+wF1QwEaORt0HuOvUqZMB7GxCnBUSx87Hprox8nN7GuqCBj84Evzg+YB4ARfSOXahakUciAqTAgXDsREzSBaEbQOaEtPGoRD12IW5pW5WUKAggP3LL7\/IiKlkmCVf4YwZMyoAu1SMoBzs8kOB1IGNCMJNI44FOF+ix7YdhpKwFcmFEVR+tsmNEpUCqRtBvfPOO+KOO+6QjgYqAWeugJ20EVRUYrv6+aNAqkZQmFLisIunOdflpEymRAV2WkZQ+dsmN1JUCqRqBMVNo+nyZC5AGfqQpzBufOxCVffZLPd0mw6dFlHdyWirbjL52xaJNsi+hHZOj+39Snmq+7C9\/u677yq0XLRokZg1a5a4+eabxZFHHilOOOEEqdsuRCOoqJxE1fey11DGTfqtX1R3MjWGOpdww8sNru6Kpl\/E8LfNHoX2XNLYXoq4685Xu9QPj7aF2kSRUjCC0tfqZ+BkPoviTqbG0A2u1q5dW8FvMgywvV6ofIEzm3GKBtgsshSMoExuGia9SJgrcxMEehueYRti+k6aFoLKNc58OeI69pq25jYxy2Ybrn+t1HNCdDzxxBMyWgElyGa8IICNo8Fnn30mjZmYECEXSGBPWg4CC6pSCkZQai26\/UeQHBvFnczGifmNyy9T5AgCthrXJsYEcVPbV8ZU2dockk1A6mtXdArzJSkIYL\/66qvyEKmCTWKDjAUfSZTQbaswZ6VkBAUwbPFCvHwbbYZQXvE\/bICxcX0\/YKtny5Ytiyxje9mdm\/OyiWOmOa7XWQR6cAbzkv9TBzax+lD5NWvWLCM8MPptbiHh4mx2oRpBJZXJy7T28zMZDXIn83Pg1blvkFYkqcA55tpsXyjzxVWihleYiSCH5NSBPXfuXOlsQMYv3XmXYOoAvnXr1tIX8qOPPio4rYiZBiRsQia\/z7jidmgygrQRNncyL+dfNaYOVi+O7WUDHiR+6M\/1l0\/J1qasb8rXvMw4lSA20QYtDnWGDh0qxo8fn5EGMcgbKVVgqxyNnNq5VvdLGV2IRlBm4qb92vUV9QcGZybwC4TD5uub4hUbRIHI5Fx+B01zs\/1EERu4GWvcuHECD3b2C\/BxIKUuhzqTy5r23vq6ecH85H4d2LaDb0FzbEVYTuMQhWt1NB8Esenevbs0hFeB1wvRCMrk2IAacAeVoEObKT+GdSfzAouXlsPrUMnvNk8fE0y6nps2irPiembzv9TjlagX1vT6MQ+sRSljIzcjbhBHZMuWLeKkk06S+dBXrFghwwhzaATs+D0WqhEUMvYfKxYKEjJFyZbrFb8jige5qR0I4mLqwKpMgsnPYuOapuZGHWhJYaKLBQCVry0iA3O56aabZORcJXLo8ryp3VAvofLDJCCpLlqZ3J8+TS97P21NqkZQujyIfNWrV69ye5FVq1ZJP75OnTpJwilVYNQEpvqnvdCMoGyHN68rddYR5E7m9\/KbogugwFTYD9i65gag9ejRQ4oi7AVA9AI22ixTfgbIiC1oMhSXtq2fgyVxGwGm7ngMFvh\/NGaUIH\/MVI2gBg4cKA+GcAImTeIkVf78808xYsQI8c0330gHXv6pxZr50p0RVJDwk8xzJWOHAXYyI\/5r72ILzOk3RqpGUCRK4tP17bffWqOP6hyoFI2gktr8fPVTjMAOuvzKhna+Po9Y7KHxMFVb6jPFf5UjQikZQWVD0LTaOmBnUt4X2ASi7Nu3rwx1ixZEhXdFfcTBEtmK\/5aaEVRa4Kws46aqx4bIuIahAUGG6ty5c4ZW5KijjiqPj03dUjKCqiwAS2udqQObhWPcxEkYNRS+j8QUgYMjegSlw3N5HtOCTmGPWxDAzjWJCtWDJtfrrsz9O2Dv3CnWr18vVY01atSQWMgHUSoz6PKx9nzsYahIULlcbFyOTbg13R48l3N0fSdLAaxGScyVmrov2eXYe4sKbC6MIArh1lwpXgrAlCZMmCAwHchFKRiObePAaGXefvttaQ\/eqFGj8vXzxgPwoEJ7HCCaNm0amYDZtFUapco0dp3VC0SdrxeUb8mWo9uLLU3be24RgFZGdHolfk8C7KkD23HgoNezOJ4PabxRdDxge\/lk5\/+2r5iwtn7kycPg+JdtSR3YLABwh+HA2S7Wtc8dBeqsfltUn317+QBlXceILU3PiTxgyXDsyCt3DQqWAnHNhHOxoILg2LlYmOuzclPAAbty73\/Jrt4Bu2S3tnIvzAG7cu9\/ya7eAbtkt7ZyL8wBu3Lvf8mu3gG7ZLe2ci8sZ8DesWOHmDZtmjR04QqcjLPdunWTHjdBdtz9+vWT9t61atXy3B2uq3Fde\/TRR6VXtl4I9rNkyRJpQ7506VI5dpg+6SPNeTM+IQ+wocCPlNKhQwdpG0NU1aCSi7kzHxxNsMkn4zLmDdAST\/rq1av\/eyFTVibjdU+dOlWsXLlSkOSWhFx4voeZe9Daoj7PCbAJAca1KH54BGdp2bKlWL58ufTGYZETJ04stxOIGswS0GI\/QqBMbqlsIR+i9qmIlva8MQjDBY91EY6CtUIrItwSCQCvJa+Si7kzH0I7EOCf1C0NGjQQ8+bNE3PmzJHgHjx4sAQ3vq\/kKgLY7PeZZ54pzY2JyLtr1y4xefLk8hxGUQEat35OgM0C8XCHo\/LWqoIHDn6Sffr0EYMGDRKbN28WhHmAoxKLQoVRgzsQV4PfCNKjCtfuRHp9+eWXpduaLTBjNgEy05w33JCgNlu3bs1wuTNppnIBmRue9NzLysokWNesWSO\/fA0bNpRD8rLxJYaLE9ORkGpEfcU39vbbb8+IP6OCmh533HFi9OjRGRw+LmDDtksc2LiSwWUQBYgUpYsJfCoJAgMHIv0ETsFhvdv1AD6EWGMcIh2ZHBsuE7ZPnUhpz1vlqYd2OjNQEQGg2ZgxY6ziWS7mTvQvxIjmzZtXSCNCElti9sF84OQzZ86UXxSyNxOQRxUV\/xFxMEwQ\/bCgDVMvcWD7DQo3glOTLxJgz549O3RSJoCN7Nm1a1dx6qmnyr9tQXqyCZDpNfd8zRtOaAYnCrOJuaK5V78q5xAcGmB7FRVYichhJQ1sRAyAzWEIrppNMEuvkGHZ9Om1Qfma9xdffCHjkBNsnxcUrtmuXTtB8CI4Z5yS5NwZXyW1feWVVyqElzbnZ4bCU+E74qwjapu8cWw+XwMGDJAaET63HESyCWbpB2yvcGtx8sXkY97Is5xJ0Cbw+cbHs2fPnoLIqBzAtm\/fLj\/1BAaNUpKeO3ND7Bg+fLjcSz3WjDkvlaqcQy3nIiK45rPkBdik1iMOIOKEvshCB3a+5q0HgCRTBLRS3G3Dhg0SQBywObD5xSnXgZP03AH1W2+9Jb+27du3z8hwYQKWAzyxHUmdiNxNmpd8l5wDm2hSqP4IR4yWRM+M4AfsoGCWcTh2UJ868fM5bz7vaA1QpRFOTqX3VvOBbmiC0EI0adIkECNJz53DKaIHMjWBk0aOHClq165tnQeaEFSxKAb4Mrdt27Y8Sm\/gxBOskDNg84ajT4b7kOWXGH\/E0tZLNgc9L2Bn0ydzS2veREklMawtBUhYESoXc0dM4ivLRRgqPXTXZmo+taf6C8XXJarolCCuRc6ArV+SjB07NuO2US0grmqO9l7AzqZP+k1r3uT7ATRoRohRrQpgRQO0cOHCQI1J0nPnS4IOG2AjgqDi8zoAKnkeTs5LyqVOmiUnwFaL5CRPcibUe7aSTTBLL2Bn02ea80aW5pYPj3o0I0qW5lKK37kIQRSoVq2alZZJz10dFJGVES24UfQCtTooIm7avsxpADxxYPPpAszTp0+X9gSIIWZBJCEDFZsUN5iln3wep89CmDdcG+1Iq1atMrQiXM74aRZyMXcuv1QkXfbKdmjlUMhcuTLnq8JVeosWLSrsN\/EeyYgR9uCbxIuQOLBR7\/Tv31\/ahngVPRGozWApTDBLP2DH6bNQ5k0WZA5diBWEdSMdyjXXXONrSJSLuTM+dw1+BZOHLl26yNtkXkqvklROyiiATxzYUQZ3dR0FckUBB+xcUdb1myoFHLBTJb8bPFcUcMDOFWVdv6lSwAE7VfK7wXNFgf8HalHNPapj8bgAAAAASUVORK5CYII=","height":98,"width":163}}
%---
%[output:8abc7834]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2306dc04]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d175c96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:849e5967]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:73560081]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:009e7c8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:630232f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ffd71de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:00782a05]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58f06463]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ea3c794]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:224e6560]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3369a97d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4de04fb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f63dcb3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5662530a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b00554e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0dc135ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2fad7a48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:648473e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57c86ac5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9ba05856]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35fd3337]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:78b1b4d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:092d462d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:81b97be6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:765a963b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99e01b07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0896794e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:882b6f8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0e5af49d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b9730be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d56c84c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2a538191]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:512f4f2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9dd71502]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4baeac86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3dbf4763]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:734812ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:984dcb31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:34874f36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:096320c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1c201b1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90daa2ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:193d03a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:904868d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:298201df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9cf454ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2bb1bda0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:17cc76da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:14905faf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7685ffc9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:8fa502e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80fd9422]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:57b82619]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21a5cb93]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:162af99d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:36dc2db9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c54fce6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:02cdc6f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0f1f6a01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7714f6f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95a7c342]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:60e44764]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d56ee40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:740ba8cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0ceea8f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9359e771]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"144×12 table","name":"KCG_result_MK","rows":144,"type":"table","value":[["'KCG'","2025","10","'daily'","'BsG_Sh'","'neph'","'y'","'MK'","95","-0.0580","-0.0333","-0.0827"],["'KCG'","2025","10","'daily'","'BsG_Sh'","'neph'","'MetSea'","'MK'","[95;-1;-1;95;95]","[-0.0924;-0.0199;-0.0433;-0.0905;NaN]","[-0.0286;0.0317;0.0031;-0.0432;NaN]","[-0.1561;-0.0704;-0.0889;-0.1368;NaN]"],["'KCG'","2025","10","'daily'","'BsG_Sh'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'KCG'","2024","10","'daily'","'BsG_Sh'","'neph'","'y'","'MK'","95","-0.0296","-0.0047","-0.0546"],["'KCG'","2024","10","'daily'","'BsG_Sh'","'neph'","'MetSea'","'MK'","[0;0;95;95;95]","[-0.0295;0.0498;-0.0706;-0.0734;NaN]","[0.0321;0.1028;-0.0228;-0.0292;NaN]","[-0.0907;-0.0024;-0.1177;-0.1174;NaN]"],["'KCG'","2024","10","'daily'","'BsG_Sh'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'KCG'","2023","10","'daily'","'BsG_Sh'","'neph'","'y'","'MK'","0","0.0042","0.0283","-0.0200"],["'KCG'","2023","10","'daily'","'BsG_Sh'","'neph'","'MetSea'","'MK'","[0;95;0;95;0]","[0.0363;0.0756;-0.0350;-0.0623;NaN]","[0.0954;0.1250;0.0096;-0.0173;NaN]","[-0.0237;0.0265;-0.0802;-0.1069;NaN]"],["'KCG'","2023","10","'daily'","'BsG_Sh'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'KCG'","2022","10","'daily'","'BsG_Sh'","'neph'","'y'","'MK'","95","0.0353","0.0595","0.0109"],["'KCG'","2022","10","'daily'","'BsG_Sh'","'neph'","'MetSea'","'MK'","[0;0;0;95;95]","[-0.0007;0.0436;0.0300;0.0526;NaN]","[0.0637;0.0915;0.0748;0.1004;NaN]","[-0.0671;-0.0048;-0.0136;0.0046;NaN]"],["'KCG'","2022","10","'daily'","'BsG_Sh'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'KCG'","2021","10","'daily'","'BsG_Sh'","'neph'","'y'","'MK'","95","0.0500","0.0750","0.0252"],["'KCG'","2021","10","'daily'","'BsG_Sh'","'neph'","'MetSea'","'MK'","[0;95;95;95;95]","[-0.0047;0.0704;0.0522;0.0796;NaN]","[0.0568;0.1186;0.0979;0.1327;NaN]","[-0.0671;0.0235;0.0067;0.0256;NaN]"]]}}
%---
%[output:355f85f0]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"4×19 table","name":"KCG_result_LMSlog","rows":4,"type":"table","value":[["'KCG'","2025","10","'month'","'BsG_Sh'","'neph'","'log'","'LMS'","2.7393","95","-0.0274","-0.0074","-0.0473","-1.9681","-0.5312","-3.4050","-0.2393","-0.2389","-0.2398"],["'KCG'","2025","10","'month'","'BaR_Ah'","'abs'","'log'","'LMS'","2.1728","95","0.0269","0.0518","0.0021","0.8721","1.6749","0.0694","0.3093","0.3102","0.3084"],["'KCG'","2025","10","'month'","'expS_bg'","'neph'","'log'","'LMS'","4.1318","95","0.0998","0.1481","0.0515","33.0912","49.1092","17.0732","1.7135","1.7171","1.7100"],["'KCG'","2025","10","'month'","'SSAR'","'abs+neph'","'log'","'LMS'","2.7853","95","-5.3485e-04","-1.5079e-04","-9.1891e-04","-5.3747","-1.5153","-9.2341","-0.0053","-0.0053","-0.0053"]]}}
%---
%[output:4cfffd88]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"4×19 table","name":"KCG_resultLMSlin","rows":4,"type":"table","value":[["'KCG'","2025","10","'month'","'BsG_Sh'","'neph'","'lin'","'LMS'","2.5559","95","-0.1054","-0.0229","-0.1878","-2.6243","-0.5708","-4.6778","-2.0537","-2.0514","-2.0559"],["'KCG'","2025","10","'month'","'BaR_Ah'","'abs'","'lin'","'LMS'","2.4167","95","0.0019","0.0036","3.3565e-04","4.2790","7.8202","0.7378","-0.9805","-0.9805","-0.9806"],["'KCG'","2025","10","'month'","'expS_bg'","'neph'","'lin'","'LMS'","3.8056","95","0.0863","0.1316","0.0409","6.3799","9.7329","3.0270","-0.1374","-0.1361","-0.1386"],["'KCG'","2025","10","'month'","'SSAR'","'abs+neph'","'lin'","'LMS'","2.7823","95","-5.2500e-04","-1.4761e-04","-9.0238e-04","-0.0530","-0.0149","-0.0911","-1.0052","-1.0052","-1.0053"]]}}
%---
%[output:09faf133]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAHYAAABHCAYAAADSgbMZAAAAAXNSR0IArs4c6QAADWZJREFUeF7tXXlsjUsbf05uY4m6sVwEdfksCRJENLG7JUGs6a3at9YSO4mtdmrfL3qtRVVyLQkifL18tpbY\/3AJf1zUFrXFdqlQyYn3u7+p0Tlv3\/edmdMXp6dnkqbn9H1m3mfmN8\/MPNvUYxiGQaESdCPgCQEbdJiyDoWADU5c1YD98OEDTZs2jfr27UvNmjUL0qEIrm5JJfbVq1c0dOhQunbtGu3atSsEbCHB3xFYSOq6desoNjaWJk+eTAkJCSFggwFY3gcutSFgCwmqODxFR0cbWGa7detGS5cupZIlS\/pwv2zZMtq0aRP7W\/ny5Wn37t1Up04dysrKov3791OPHj0oIiKi8PS4iHDq2bNnj9G9e3fLwxE\/NHXp0oU2btzosxRfvHiR+vXrF9p3A3SieC5cuGDgpLt37166f\/8+A09cgidNmkRjxoyhRYsWhYANUBCt2PLcunXLwNIKYM+fP++zHHOp5MYpLNc4TKGEJDawUXYEVmSdL8stWrSg3r17h4ANbFzJcSkWJROf27RpQ\/Xq1WNLckhiAxtZx8PTsWPH2N66detWevPmDTNUrFixgjp06BACNrBxzVN3RowYwSQROmtiYiLNnTuXbt26RRMmTKDnz5+zbtSqVYsWLFjAjBRcYvG8adOmAd7Nws8eVEodtdLR8gTwoLdCv0WBvZjvsdBjp0yZQpcuXSr8o1YIegABwo9q8WRnZxvh4eGW9E7AogLAxU+ofP0R0JbY1NRUY9CgQbbAwvK0bds2ZpEKWA\/P\/fu5\/Neo8fVHuJC8wTN+\/HjDypQI\/rHfwkAxY8YM1p3FixfTqlWrqFy5cgXrnioQMrqMDKLERCL85gXgpqQQRUVZ8yhrs2A9C5jazFYMibQDix+SwHGB3XaqQKjQAdB58+wHEuDGxeU9V2lTbK2QTwApsFYjJ4K9ZMkSZrCQFlUgVOgglW3bSl+ZFBNDv\/31Fy0uXpz6\/P23Lf2UChVof3g4Mc0AJ3zdVUDKScEJdIMdPBMnTjSwxBYvXlzp7X4tz5AWBSCYBDpJobjccoly4Pp6kybUYOVKtXenpxOdPq23ChAx+3ru9v719nd\/gh08aWlpRufOnZVABRGkVftABVDFfVD5bS4QYq9VeTeAUZgshAkQFUUZGRlM38fvvPlWg1JSUijKbn\/3ozv+Bjt4cnJyDFVp5cBCtx03bhylpaURXHpOinNYVhaF1anjR5cCs4q3f39a9G9\/5jmsLMnJyTRgwAC\/OuD1eundu3cEFTQsLOxLG\/jcs2dP5SgWacyTeT+tXr06M1rAwR4XF0erV6+mRo0a2XYCwFZr08avTgZqJY8CY+aDZlZWLkgREV7H2h8\/fmSWvgoVKvhsj2XKlKEhQ4a4A6zVfopIRTjdR48ezQ4bO3fupMjISEdmS5iiMhTGxR0S1WVY823\/wd4qqdO\/f39mY8dKDeG+dKnElxoANznZa6mRZWZ66enTpxQZ+ROVKJFXx1WJtdpPsfRCYrt27cqWYyVbcXy82v6lOcBScozojh2uv1t+Hs\/lbN68dMezINj75ZdcWmzvqan5VfK5c\/PsLmvXrnVHYu1MiliOERZz79496diGCNwbgejoaHVgMzMzDXhtrIqTrRjPmjdv7h7XoZakI6CTZuVJSkoyxo4dawusnWoDYBEEh6W4atWqUqYYwYMHuWuN2QSI9UYsqnSf6+Tk5LCTJA4Y4klSu0281+60C3UoLo5SMzJ8VBz7jsOkOVhtXBSocFxITxcsaZI6nrFjxzJbcalSpfKROhkjCuyPxSBWry7vkgIdVIQXL14wYMUDh23jTm1+NlI0y8n5Uj0rLIyyhg0jGjyYTp8+zZwhzgXGihQisrFXy3ttSaGTFyk1KdrZiouCPzbC6yWAai5YHXjwgfkZVgzooP\/889pP+Oyr4UijauCSAuvEXVH2x2ZmZjJ3ptnyhMiTunXr0vTpzZQMXjrouyqxOi8uqrSwF5ttxSrmcS59KpbM3D1WfYSVJFbXs6D++uCmhAoNFd6qfD6LMT1W5h+RuZit2pcC649nIbjh0uudnRsYigB3F6tMALPiIOPCEVh\/PQuyl+o8R4bC9OnTWRXYpHlQgJVPmK8shw8fZvQ88hKf\/fIh6zCqQPv77\/+l1avHW\/YlJmYtvX49gXJy8hLLf\/zx1b\/24iUUHr7fpy9QQTdv3szaqVKlCvMoIZtDLFKJBfH3SqO8ffu2TzgOOvT48WOaOnUqzZo1K1\/IzvHjx7\/kH3Ge+\/TpQ+3bt\/86IT4KYHISnb5MmPAb3bz5P8u+2CXQmVmRenfsgA2lUWqg6hIpVi9xkpqlNJ\/E+pMfG0rxcAktjWYg9S9fvmR2el7E7cYH2EDOj4WT4cyZMzR4sHumOY1xDHhSc6KcD7CBmh8L6w78mbDk2NmyA37kvwGDOHegiHnN+K6UH2sl9jpLsegl4lchiKddc7QBvBhHjx6lsmXL0s2bN0MSK0wQjCXs1TzPColyVneDSPNjRTWhYcOG7K4nnfxYDqB4x4V4QkTiFxz38+fPZ0lgCA0pVqwY20uQ7AXVJbQU5yKLscTYi+qO7R6rshTzTIBhw4ZRu3btWCYeBxz6FCQONAhqExO3QPMAnhQin2x5MXse+wSyDVCvYsWKBGmFlKIunuEz3lezZs1vsLAFzysc82O3b9\/OFOEjR44QfJ4Ih8F+h6Aq0W0HfyxmEaQM+mX9+vV9Rig9PZ1u3LjBlGxII74jrgfSn52dTcuXL2efzfXwDNcndOzYMXhG3M+eaCdlRUZGGvBljhw50jI\/dvbs2XTnzh3GDiItILVYDsxuO+7d93jyx\/B9+vSJ8Hf+TPwuRgVY1fVzHIKumnYapVPujmoaJSTw1KlTTKLv3r3LpE8sVhLLJRh76oYNGwhRfTqJvUGHnKRD2hL79u1bo3Tp0pbNqkT9w3QnHnCwf6KI+TzmG2msDk92GX9FDcAv2QiqHnWbAfJcuXLFaNy4seVjv\/J0LFqyumqIn5btjNhFDlDdbEDJAHnmzJljYB+1CwJzNY2yyKGV12HH5C2VDEMxJVRhHD2xsbHGli1bmDHAqYhSB9syt1cqp1EqMOMmiejCE9194juwJcTHxzOPEYqoE7qVRSdN3lIJtSCij0eP0tQ\/\/1S+M1rJbccHoEmTJuzkPHPmTPez3N1E9bMyz68QhCqG0BVzHq+VRUwKhAafyMZzSt6CHzXOHP5v0\/6pn3+mYT\/8oJx87hk+fLiBxCq7C0Yw83HXE1Sdq1evUkxMDK1Zs8byXgq3ZrnG2FmSWt0ix2+\/efbsGasDoM33RyoBobgkYoK0lcW8EBHCmFSDVH\/VyQRwulwEA8BPuTzLDsAeOHDAJ40SgAJ8c8Qe0glbtWqVb\/D5TTNuqTd2qYfii8+ePevII57D1ykrcOZb9clcr1OnTkqB5QgBRwSyShnRoQPFJyYqXQbuef\/+vWG+o5i\/BEtwamoqW3qxr2LWc2B5GiUAFwE1M4ib3ECLgqUPiUXi3VAAF5angvyvAZgtEefboEEDy8x8vBM\/dgU87tu3T+nOKvQF9E4FExfXE6oW1f+PoyWxTv+eRfTAcCaRr\/P+\/XuWRgnj\/JMnT6T8Q8+FRHAXk1WFpKQkZqQQy6NHj9hXuxSSc+fO0cKFC\/NNFKR5tmzZktUFDRwQbhaogbxY8Yi\/YZKpFqS2SS86iIqiX8uUUU\/KcgJW1GMhsZipsB+vXLmSWZmQtQ0bsqxAKlUu+uLuOxXJVpFCSBesYG7fHgfnP\/rjtProOC1kEuuNiCBvcjL13bzZHWBFy9P169fZkoy9FFcb9OrVix4+fCjDVOs59i5ImpNk41SOg8+oUaOkbUNyVeikDZkIwIOMR34wk7WNviS3bElhf\/xhTxoVRd6ZM2n9+vXqwPp75R6W1tatW8v4Dj13cQS00igPHTpk2O1BofxYF1FxoSktYBMSEgxEL8BPai5OTgA8g56msse60Ce\/m4CpFOqQrNSuXZvwg5Acq4J2cM3PiRMnZE19eQ5XKFRBaA3iOKEtuD+tCniFHg5NRTTz4pyCrVC1OJoUnZwAODzgDoqDBw86voszpzK4qkzr0FWuXFl6cgePuKUFubUA4PXr1\/mAgMkVz3XOFdWqVfMBB2Ngm5gt6ZSr\/li8y8kJAHCh5yJqwqrgYIArg2C1GjhwoCProOWWKx3gnGghYbBlnzx5UsqjVZSGFRA61zNcuHDBra4wX7WOQUfJVizjzs6+ilglAIuyY8cOZnB3mgBox8nYweuqTALQiLekqfAo6yeeY\/tR4RGTCgEG36u4AqzIvFWuKH8uG1wV+yoHDO+RTRRMLKvixKMMCB0e3bx6T8aX+bnrwKoyYDe4KpLNAZNNFFVedOl0eNRt2y367wasUwesAMOhA3HNMMLjbmVzWiRikCtVqsTCWL9FGqU\/k8qtlFDX0ijdmkW67cAJgXhlhMACzEBNo1RZ2nXSKHGLu11KqHIa5eXLlw1Igs4NqLoAuUXPDSZOPmHxXdzBDg+U9lW8bjFt046\/fVFOo3RyAnzlvmk3LwJl929jeKOiDo5Adhm9NjMFrFCQvqikUf4fvqa7s2aOZyoAAAAASUVORK5CYII=","height":98,"width":163}}
%---
%[output:01cbb4fd]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAHYAAABHCAYAAADSgbMZAAAAAXNSR0IArs4c6QAACtpJREFUeF7tXWloVFcbfkOjRklLVCKuTVT8EeHTH3WHJlYwUlxI4m5+GHfcGnCPC2kUd4xLiwtEUwQ30B+i1nwKpioYtb+kQkmsySjBKC6NKBoxOO1zzBnPXO\/y3smdmTvTe0BMMufeszznPee8y\/NOgt\/v95NX4m4GEjxg4w5TMSAP2PjElQfs27dvac2aNTR9+nQaNmxYnE5FfA3LUmJfvHhBc+bMoTt37tDx48c9YGMEf1NgIan79u2jSZMm0YoVK2j16tUesPEArByDlFoP2BhBFZennJwcP7bZ8ePH07Zt26h9+\/ZBvd++fTsdPHhQ\/K1z58504sQJ6tevH9XX19OZM2do4sSJ1LNnz9gZ8X+kpwknT570T5gwQfdyJC9NY8eOpQMHDgRtxTdv3qQZM2Z4565LF0pCVVWVHzfdU6dOkc\/nE+CpW\/Dy5ctp8eLFtHnzZg9Yl4Ko162EmpoaP7ZWAHvjxo2g7VhKpTROYbvGZQrFk1h3o2wKrNp1uS2PGDGCpk6d6gHrblzJdCtWJRM\/Z2ZmUkZGhtiSPYl1N7Kml6dLly6Js7WsrIxevnwpDBU7d+6k7OxsD1h34\/pJ3VmwYIGQROisJSUlVFxcTDU1NVRYWEhPnz4Vw+jbty9t2rRJGCmkxOLzoUOHunyYsd89qJR21EpTyxPAg94K\/RYF9mJ5xkKPXblyJd26dSv2Zy0GRgABwj9uCRlYNABw8c8r4Z8BxyUWlqfDhw8Li5Tn4Qk\/gOTzfWwkPb1VjZlKLM5bGCjWrl0rGtmyZQvt2rWLOnXq1KpGvYd1ZuC334hKSojwvywAt7ycaORI21Nm6baTlyS82XPb2Z5f3gMA9McfjesC3IIC3rtaagWAVf2uqkMAW\/GhQ4dE9e7du1N5eTk9f\/5c2IlRtm7dKgwWbivSoHLu3DkaOHCgOE60O829e\/do1qxZ9OjRI9F9qRlEdCyQ0O++s2zyXUUFrfr1V3awQwBYAJienk6qQwAToj1XY2V7Vm3fcmzaBaje+rVeLcuZdqoCQFW3X4P3Xvn6a5r7xRfsXVMAq\/W3ykmZN29e4IyFPRkFk+H2C5We+VOqbSqAeo4Pp\/BivQcXpd69WVVRKTcnhx3sEABWXpJUh0BeXp7YqmTBVpWVlSV026VLl9KFCxcILj07ijN7FDYqNjc30+vXryk5OZkSExNtPBndqon19ZTYIjCcnizIzqZZJSWsKBZTYFXHu5SCrl270uPHj4WDvaCggEpLS8UZFs3y7t07YR1LTU2ldu3aRbMrttvu3acP+5mQJBZ2YBn6YrRFYQvGRePhw4e0aNEicdk4evQoDRo0iN25cFTEonvy5Alh0SUlJYWjibC9M+n771lnLFSe3JQUe1sxeq13ecLfr169GrAhA\/yFCxeKrXjcuHHCpOgGFaipqYkaGhqoW7duMQesuDhZ3Ypb9NncvXv5wMoIClXdUa\/9qroj\/+423TamgYX0\/PILkXKXCdoeACp02OJiys3N5QNrRfFQQVR1Vjf5Y2MeWCBpZHkqLrZtnMDrQjYpesCG59jFIn36+++UOnhwq44VS++Okc4KYGHMwO042v5YqDuNjY1C3Ym1y5N2eZiNBVoIt1gCa+aP7dWrF7cdr54DM2CH8RoysOin5491AC3mKzx\/LHOi4r0a6\/K0bNky2r17t4iHOnLkiOePjYFVYemPRaQi9NeEhARav349zZ49OwaG5XXR9TRKmDeLiooEUqpfVU+\/Vn2wqK8aWoz08UguAafGoucjl943OR5LiUXFaNEo4QhXw3GkrXrVqlVi99CG7Fy+fDnAP5J9njZtGo0ePTrqIT5OjcWIQKddoCED69EoIynrH9uCxKuLVCulao9C5se6yfIU+SmOTouQejUsSXvcBAHrZn5sXV0dXbt2jWbOnBmdmXR5q9pIkSBg3cqPRUQEOEOIiFiyZInLpzh63cO9A0XlNeN3Fj9Wdlvvlsnxx+oFjak3RO07YDqrqKigjh07UnV1tSexyrrBXGp95Hq5QSz5saqaMGDAABH+aIcfKwFUQ1rVGyKIX7BHb9y4UZDAEObStm1bcZbAuYDwUW8r\/nR5wtzr+cg\/uxVztmKpVsydO5dGjRolmHgScMQcQ+JQB\/ZMlbiFOg8ePBBtqmx5lT2PcwKBdHiuS5cuBGmFlOJZfIaf0V4fG7FB0dsY3dOyKT8W5kMAd\/HiRYKfEOEwOO9gfVJplD169BCrCFIG\/bJ\/\/\/5BI6ysrKS7d+8KgwGkEb8jIA7S\/+rVK9qxY4f4WfscPsOCGDNmjHtmLEo9se0EkOmAjPixGzZsoPv374vhgB8LqcV2oKVRSpcSTI\/a8uHDB2GSlJ+pv6uuKL1nozSPrms2KjRKSOCVK1eERNfW1grpU4uexEoJxpm6f\/9+ys\/Pj3p8suvQVDpkW2LNYp44Uf8w3akXHJyfKCqdQpuRRu\/ypJc8jE0ptEM95Nbl1nPpamC57VpLo\/ws1ZDPR+fPn6cfSksDRK8g8xiXUsith8nn1uXWk4C6dAFY2oodDTXlTBqXUsitBwC4dbn17CyUKEm0pdsOagi2WhS9uGI2jZIzaYihtQqeRkfAJTXjk7ZM5k95eXSjqopONDRYTu\/elBQqbGy0rEeVlYiid5zPatWw3ZzRpsCqVA\/pBoOh4fr16\/ZY7pxod4wMwMqtzWqkjM\/\/+OYb+t+XX\/IoFIz3iSrcPmIBhMBE1+tGKDmjLbditSHoqlBVbt++LYjEbdq0EU7wKVOm0JAhQwynhs1P4U5uDNRrzs+n5rKyVvcUkgqtATmjsXtyUwuzgZWEZ9AmIbGgUSKtLbK1gQxlVEAV7JWZ2eoBxuIL6mprbXcbccVv3ryhDh06BFFCU1JShGHIMWDVRJkgQsOkCNuupFFasu18PkrKyLA9wHh4oOnPPwPZX7iXZyPmILxckydPdgZYSCkk8\/379zR\/\/nw6duyYMD4gd7GkUXK8O6RjjYoIcPKMY6QCCEt\/\/H62liXb\/+uvZsEcHDw49TNWg2OkLEQoQvRPnz4t9E3s8diKJY0SlAOk37NktINN5mThXmDsAMt9p51xFBQIIp1R0TI2INXabEDq\/Qs0Fke2YtgnpaojO4dVgy0BWWMQ4eCVyM1Ajt0cFEZdg0N33bp1gXQ5qhMAZ+\/w4cMjNyqvJeHS5BbTWzFCP8+ePRtw22EbhrTCxAhgoeYgzS2bbYe9Ri\/7GDigrdD5BPWwJQdFxNh2JmPxpY9k2Vm4IKEepqeykp+Gz5atGKkKcK5q\/bFsYNWRwAGflmZnbIZ1oSI8e\/aMoBJEDFiTsTQ396Rvv3X+m01sCKz1V6BpQ2MgofLrWby0tsbrsq7Ovg5rtcpxpeHmzmQZKIzslB6N0hiKoqJhnIRrVlgGfR6yxKq5BWXwmYxJwrlqFnluq4f\/gcoc87iUPo55\/OMZy584yzNWppH3UtnyJ1XW5CSDycqydmiFkt3WFFg1zFF2lu2msz8PcfkEJxkMZwFAcbBTWGesnRc6Xdcp6qEbaJQ\/\/3yeSkt\/EFOkpYTm5e2lv\/8upKamT9\/P+9VXL\/5NIbiVkpPPBFFCHaNROg0W931OUQ9jjUZZWLibqqv\/r0sJdZRGyQUi3PUkVQRZWffs2WP5XQUyjWBaWprrUvGGOhY2jdIqM1u4wbLzfhUoozRF8n1qwmwEslvVt9MPJ+q2Ziwyuzv6YZQV\/R8FJNOkWjOOqgAAAABJRU5ErkJggg==","height":98,"width":163}}
%---
%[output:722a53f7]
%   data: {"dataType":"text","outputData":{"text":"Fichier NetCDF créé : KCG_tr_def.nc\n","truncated":false}}
%---
%[output:737ce877]
%   data: {"dataType":"warning","outputData":{"text":"Warning: there is a lot of days with less than 50% data coverage"}}
%---
