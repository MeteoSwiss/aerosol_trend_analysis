BRW_st.name='BRW';
BRW_st.lat=71.32;
BRW_st.lon=-156.6;
BRW_st.alt=11;
BRW_st.env='P';
BRW_st.footp='P';
%%

BRW_rd=read_betsy_2026('M:/pay-proj/pay/aerosol/matlab/github_trend/raw_data/brw_1976_2025',BRW_st.name); %[output:20e02795] %[output:0d6482bf]
%%
BRW_rd.U1_S11(BRW_rd.U1_S11>999)=NaN;
BRW_rd.U0_S11(BRW_rd.U0_S11>999)=NaN;
BRW_rd.U_S11(BRW_rd.U_S11>999)=NaN;
%%
names=fieldnames(BRW_rd);

names_delete=names(contains(names,{'g','N'}));
BRW_rd_short=removevars(BRW_rd,names_delete);
%%
% X1_A82 are AE33 data --> compute AE33 abs data

sigma=[18.47 14.54 13.14 11.58 10.35 7.77 7.19]; %for the 7 wavelengths with an exponent of -1
BRW_rd_short.Bax1_A82=BRW_rd.X1_A82 .*sigma(1);
BRW_rd_short.Bax2_A82=BRW_rd.X2_A82 .*sigma(2);
BRW_rd_short.Bax3_A82=BRW_rd.X3_A82 .*sigma(3);
BRW_rd_short.Bax4_A82=BRW_rd.X4_A82 .*sigma(4);
BRW_rd_short.Bax5_A82=BRW_rd.X5_A82 .*sigma(5);
BRW_rd_short.Bax6_A82=BRW_rd.X6_A82 .*sigma(6);
BRW_rd_short.Bax7_A82=BRW_rd.X7_A82 .*sigma(7);
%%
names_short=fieldnames(BRW_rd_short);
names_delete=names_short(startsWith  (names_short,{'X'}));
BRW_rd_short=removevars(BRW_rd_short,names_delete);
%%
datevec(BRW_rd.Time(1))
datevec(BRW_rd.Time(end))

prctile(BRW_rd.U_S11,[5 50 95])
prctile(BRW_rd.U1_S11,[5 50 95])
prctile(BRW_rd.U0_S11,[5 50 95]) % max>25%, no RH problem
plotFigControl(BRW_rd_short,BRW_st.name);
% 2 size cut, ok
% sc: ok, max in spring, some negatives
% abs: 2 size cut, max in spring, some negatives
%Clear seasonal cycle in BbsF with max in summer
%clear seasonal cycle in SSA with max in autumn
%%
% break point detection: scattering
% it has to be done on the merge Bs2_S11 (TSP) and Bs20_S11 (PM10 in between)
%it was said in 2019 that there is no blue data Apr 93- Oct 97. Not right blue has missing
%period --> use G
P1=timerange('1976-01-01','1978-01-01');
P2=timerange('1977-01-01','1997-10-01');
P3=timerange('1997-10-01','2026-01-01');
% B

BRW_rd_short.BsB0_S1S10=NaN(size(BRW_rd_short.Bs1_S11));
BRW_rd_short.BsB0_S1S10(P2)=BRW_rd_short.Bs1_S11(P2);
BRW_rd_short.BsB0_S1S10(P3)=BRW_rd_short.Bs10_S11(P3);

%G
BRW_rd_short.BsG0_S2S20=NaN(size(BRW_rd_short.Bs2_S11));
BRW_rd_short.BsG0_S2S20(P2)=BRW_rd_short.Bs2_S11(P2);
BRW_rd_short.BsG0_S2S20(P3)=BRW_rd_short.Bs20_S11(P3);

%R
BRW_rd_short.BsR0_S3S30=NaN(size(BRW_rd_short.Bs3_S11));
BRW_rd_short.BsR0_S3S30(P2)=BRW_rd_short.Bs3_S11(P2);
BRW_rd_short.BsR0_S3S30(P3)=BRW_rd_short.Bs30_S11(P3);

%4
BRW_rd_short.BsQ0_S4S40=NaN(size(BRW_rd_short.Bs4_S11));
BRW_rd_short.BsQ0_S4S40(P2)=BRW_rd_short.Bs4_S11(P2);
BRW_rd_short.BsQ0_S4S40(P3)=BRW_rd_short.Bs40_S11(P3);

%%

names_sc_short={'BsB_S1S10', 'BsG_S2S20','BsR_S3S30','BsQ_S4S40'};
break_BRW_scat_pettitt=change_point_analysis_V3(BRW_rd_short,names_sc_short(1:3),0.05,'BRW - Pettitt','neph SC');
T_BRW_scat_pettitt=make_table_breakpoints(break_BRW_scat_pettitt,names_sc_short(1:3));

break_BRW_scat_snht=change_point_analysis_SNHT_V1(BRW_rd_short,names_sc_short(1:3),0.05,'BRW - SNHT','neph SC');
T_BRW_scat_snht=make_table_breakpoints(break_BRW_scat_snht,names_sc_short(1:3));
%%
% idem variable for PM1
BRW_rd_short.BsB1_S11=BRW_rd_short.Bs11_S11; 
BRW_rd_short.BsG1_S11=BRW_rd_short.Bs21_S11; 
BRW_rd_short.BsR1_S11=BRW_rd_short.Bs31_S11; 
BRW_rd_short.BsQ1_S11=BRW_rd_short.Bs41_S11; 
%%

%%
% Backscatter begin with TSI in October 1997
% break point detection: 
% use at it the S0 data for TSP

BRW_rd_short.BbsB0_S11=BRW_rd_short.Bbs10_S11;
BRW_rd_short.BbsG0_S11=BRW_rd_short.Bbs20_S11;
BRW_rd_short.BbsR0_S11=BRW_rd_short.Bbs30_S11;
BRW_rd_short.BbsQ0_S11=BRW_rd_short.Bbs40_S11;
%use the 1 data for PM1
BRW_rd_short.BbsB1_S11=BRW_rd_short.Bbs11_S11;
BRW_rd_short.BbsG1_S11=BRW_rd_short.Bbs21_S11;
BRW_rd_short.BbsR1_S11=BRW_rd_short.Bbs31_S11;
BRW_rd_short.BbsQ1_S11=BRW_rd_short.Bbs41_S11;
%%
names_short=fieldnames(BRW_rd_short);

names_delete=names_short(contains (names_short,{'Bs1_S11','2_S11','3_S11','4_S11','10_S11','20_S11','30_S11','40_S11'}));
BRW_rd_short=removevars(BRW_rd_short,names_delete);
%%

names_bsc_short={'BbsB_S10', 'BbsG_S20','BbsR_S30','Bbs4_S40'};
break_BRW_bscat_pettitt=change_point_analysis_V3(BRW_rd_short,names_bsc_short(1:3),0.05,'BRW - Pettitt','neph BSC');
T_BRW_bscat_pettitt=make_table_breakpoints(break_BRW_bscat_pettitt,names_bsc_short(1:3));

break_BRW_bscat_snht=change_point_analysis_SNHT_V1(BRW_rd_short,names_bsc_short(1:3),0.05,'BRW - SNHT','neph BSC');
T_BRW_bscat_snht=make_table_breakpoints(break_BRW_bscat_snht,names_bsc_short(1:3));
%%
% Break point Absorption PSAP-CLAP
% use BaG0_A11
names_abs_clap_short={'BaG0_A11'};
break_BRW_abs_clap_pettitt=change_point_analysis_V3(BRW_rd_short,names_abs_clap_short,0.05,'BRW- Pettitt','PSAP-CLAP');
T_BRW_abs_clap_pettitt=make_table_breakpoints(break_BRW_abs_clap_pettitt,names_abs_clap_short);

break_BRW_abs_clap_snht=change_point_analysis_SNHT_V1(BRW_rd_short,names_abs_clap_short,0.05,'BRW-SNHT','PSAP-CLAP');
T_BRW_abs_clap_snht=make_table_breakpoints(break_BRW_abs_clap_snht,names_abs_clap_short);
%%
% remove non necessary neph variables
names_short=fieldnames(BRW_rd_short);
c=contains(names_short,{'dry','11_S11','21_S11','31_S11' ,'41_S11','Q'}); 
names_delete=names_short(c);
BRW_rd_short=removevars(BRW_rd_short,names_delete);

%%
% AE data are not homogeneised and the lack of AE from 2001-2010 don't
% allow an homogeneisation
%AE16 (1988-2001) could be homogenised with PSAP-CLAP:
% common period= 1.8.1997 - 9.12.2001
Pt=timerange('1997-08-01','2001-12-11');
ratio_med_h=nanmedian(BRW_rd_short.BaG0_A11(Pt)./BRW_rd_short.Bac1_A81(Pt)); %0.71
ratio_STN_h=nanstd(BRW_rd_short.BaG0_A11(Pt)./BRW_rd_short.Bac1_A81(Pt)); % 1.68 namely two time the mean !
ratio_mean_h=nanmean(BRW_rd_short.BaG0_A11(Pt)./BRW_rd_short.Bac1_A81(Pt)); %0.82
ratio_mean=nanmean(BRW_rd_short.BaG0_A11(Pt))/nanmean(BRW_rd_short.Bac1_A81(Pt)); %0.76
ratio_median=nanmedian(BRW_rd_short.BaG0_A11(Pt))/nanmedian(BRW_rd_short.Bac1_A81(Pt)); %0.77
figure;
plot(BRW_rd_short.BaG0_A11(Pt),BRW_rd_short.Bac1_A81(Pt),'.');

% fit between both data: slope 0.9955x+0.1021
fit_test=robustfit(BRW_rd_short.BaG0_A11(Pt),BRW_rd_short.Bac1_A81(Pt)); % y coor=0.076, slope=1.081, 1/slope=0.9251
%since the fit pass over (0,0), forcing him through 0,0 will lead to lower 1/slope
%use of the mean could be ok
%%
% AE16 has value of zero absorption which are probably not right --> NaN
BRW_rd_short.Bac1_A81(BRW_rd_short.Bac1_A81==0)=NaN;
%%

% lets correct AE16 with either 0.76 and 0.82 and detect for BP
BRW_rd_short.BaG0_homo1=BRW_rd_short.Bac1_A81;
P1a=timerange('1978-01-01','1997-08-01');
BRW_rd_short.BaG0_homo1(P1a)=BRW_rd_short.Bac1_A81(P1a).*0.76; %prefered cte by Martine
P2a=timerange('1997-08-01','2026-01-01');
BRW_rd_short.BaG0_homo1(P2a)=BRW_rd_short.BaG0_A11(P2a);

% % % BRW_rd_short.BaG0_homo2=BRW_rd_short.Bac1_A81;
% % % P1a=timerange('1978-01-01','1997-08-01');
% % % BRW_rd_short.BaG0_homo2(P1a)=BRW_rd_short.Bac1_A81(P1a).*0.82;
% % % P2a=timerange('1997-08-01','2026-01-01');
% % % BRW_rd_short.BaG0_homo2(P2a)=BRW_rd_short.BaG0_A11(P2a);
%%

figure; %[output:1343ce23]
plot(BRW_rd_short.Time,BRW_rd_short.BaG_homo1,'.'); %[output:1343ce23]
hold on; %[output:1343ce23]
plot(BRW_rd_short.Time,BRW_rd_short.BaG_homo2,'.'); %[output:1343ce23]
ylabel('abs coef homogenised');
legend('cte=0.76','cte=0.82');
%%

names_homo1={'Bac1_homo1';'Bac1_homo2'};
break_BRW_AEhomo_pettitt=change_point_analysis_V3(BRW_rd_short,names_homo1,0.05,'BRW- Pettitt','AE16+CLAP homo3'); %[output:58a62d12]
T_BRW_AEhomo_pettitt=make_table_breakpoints(break_BRW_AEhomo_pettitt,names_homo1); %[output:7283f847] %[output:38cde80a] %[output:728007cf] %[output:05c31929] %[output:37eec7a7] %[output:0276ef3a] %[output:487a130a] %[output:2b0e1e85] %[output:01715020]

break_BRW_AEhomo_snht=change_point_analysis_SNHT_V1(BRW_rd_short,names_homo1,0.05,'BRW- SNHT','AE16+CLAP homo'); %[output:748b463c] %[output:31dedadb] %[output:08d8b851]
T_BRW_AEhomo_snht=make_table_breakpoints(break_BRW_AEhomo_pettitt,names_homo1);
%%
% homogenise AE13 and AE33. 
Pt=timerange('2014-08-15','2016-08-19');
ratio_STN=nanstd(BRW_rd_short.Bac3_A81(Pt)./BRW_rd_short.Bax3_A82(Pt)) % NaN
ratio_mean=nanmean(BRW_rd_short.Bac3_A81(Pt))/nanmean(BRW_rd_short.Bax3_A82(Pt)) % 1.62 --> 0.606
ratio_median=nanmedian(BRW_rd_short.Bac3_A81(Pt))/nanmedian(BRW_rd_short.Bax3_A82(Pt)) % 1.65 -->0.617
figure;
plot(BRW_rd_short.Bac3_A81(Pt),BRW_rd_short.Bax3_A82(Pt),'.'); %basic fit : y=0.661x+0.0279

% fit between both data: slope 0.9955x+0.1021
fit_test=robustfit(BRW_rd_short.Bac3_A81(Pt),BRW_rd_short.Bax3_A82(Pt)) % y coor=-0.0123, slope=0.64,
cte_ae31_ae33=0.64;

%%
% complete ae31 ae33 time series
P1=timerange('2010-01-01','2016-08-19');
BRW_rd_short.Bacx1_A81=BRW_rd_short.Bac1_A81;
BRW_rd_short.Bacx2_A81=BRW_rd_short.Bac2_A81;
BRW_rd_short.Bacx3_A81=BRW_rd_short.Bac3_A81;
BRW_rd_short.Bacx4_A81=BRW_rd_short.Bac4_A81;
BRW_rd_short.Bacx5_A81=BRW_rd_short.Bac5_A81;
BRW_rd_short.Bacx6_A81=BRW_rd_short.Bac6_A81;
BRW_rd_short.Bacx7_A81=BRW_rd_short.Bac7_A81;

BRW_rd_short.Bacx1_A81(P1)=BRW_rd_short.Bacx1_A81(P1).*cte_ae31_ae33;
BRW_rd_short.Bacx2_A81(P1)=BRW_rd_short.Bacx2_A81(P1).*cte_ae31_ae33;
BRW_rd_short.Bacx3_A81(P1)=BRW_rd_short.Bacx3_A81(P1).*cte_ae31_ae33;
BRW_rd_short.Bacx4_A81(P1)=BRW_rd_short.Bacx4_A81(P1).*cte_ae31_ae33;
BRW_rd_short.Bacx5_A81(P1)=BRW_rd_short.Bacx5_A81(P1).*cte_ae31_ae33;
BRW_rd_short.Bacx6_A81(P1)=BRW_rd_short.Bacx6_A81(P1).*cte_ae31_ae33;
BRW_rd_short.Bacx7_A81(P1)=BRW_rd_short.Bacx7_A81(P1).*cte_ae31_ae33;
%%
names_short=fieldnames(BRW_rd_short);
names_delete=names_short(contains (names_short,{'Bas1_','Bas2_','Bas3_','Bas4_','Bas5_','Bas6_','Bas7_','Bac1_','Bac2_','Bac3_','Bac4_','Bac5_','Bac6_','Bac7_','Bax1_','Bax2_','Bax3_','Bax4_','Bax5_','Bax6_','Bax7_','Ba1_','Ba2_','Ba3_','Ba4_','Ba5_','Ba6_','Ba7_'}));
BRW_rd_short=removevars(BRW_rd_short,names_delete);
%%

break_BRW_absAE=change_point_analysis_Def(BRW_rd_short,{'Bacx3_A81'},0.05,'BRW','AE31+AE33');
T_BRW_absAE=make_table_breakpoints_def(break_BRW_absAE);
%%
% remove other non necessary variables
BRW_rd_short=removevars(BRW_rd_short,{'Bbs1_S11'});
%%
% remove Basx and Bax since I don't know what they are
% % % names_short=fieldnames(BRW_rd_short);
% % % c2=contains(names_short,{'Bas','Ba1_','Ba2_','Ba3_','Ba4_','Ba5_','Ba6_','Ba7_','X'}); 
% % % names_delete=names_short(c2);
% % % BRW_rd_short=removevars(BRW_rd_short,names_delete);
%%
lambdaSC=[450;550;700]*ones(1,4);
lambdaAE=[467;530;660];
lambdaAE7=[370 470 520 590 660 880 950];
% compute exp and SSA separately, BRW has too many instruments
%compute exp scat and BF
names_short=fieldnames(BRW_rd_short);
names_neph=names_short(startsWith(names_short,{'Bs','Bbs'}));
BRW_expS_cal=compute_exp_D(BRW_rd_short,names_neph,lambdaSC);
%compute abs exp
Ca3=startsWith(names_short,{'Ba'}) & contains(names_short,{'_A'});
names_3w=names_short(Ca3);
BRW_expA_cal=compute_exp_D(BRW_rd_short,names_3w,lambdaAE); %[output:7fa5c406]

% % Ca7=startsWith(names_short,{'Ba'}) & contains(names_short,{'_A8'});
% % names_7w=names_short(Ca7);
% % BRW_expAE_cal=compute_exp(BRW_rd_short,names_7w,lambdaAE7);

% compute expA from AE31+33
Ca7=startsWith(names_short,{'Bacx'}) & contains(names_short,{'_A8'});
names_7w3133=names_short(Ca7);
BRW_expAE_cal2=compute_exp_D(BRW_rd_short,names_7w3133,lambdaAE7); %[output:927a2d76] %[output:152ac259] %[output:939222bb] %[output:87395408] %[output:1b978798] %[output:929a6dcb] %[output:9775063d] %[output:59336d8b] %[output:1fc14498] %[output:201c41fc] %[output:7e53a1c6] %[output:951efc0e] %[output:1bb09c71] %[output:3a59b96d] %[output:9df110ac] %[output:7ed81e2c] %[output:7b77b11c] %[output:94982b7e] %[output:0efde159] %[output:9c0bd84a] %[output:75f06b72] %[output:7ae89966] %[output:1eaa267e] %[output:6acfa036] %[output:6100d8ff] %[output:112bdc49] %[output:9c2ee56a] %[output:65d976e0] %[output:16e76ae8] %[output:07769fac] %[output:9f34e255] %[output:617c23c6] %[output:10fb6a65] %[output:8110d059] %[output:17ce46f0] %[output:498574bb] %[output:7c5b7162] %[output:0a5bdfc0] %[output:7c3b09c1] %[output:4ef3699f] %[output:72605d9e] %[output:79b6a064] %[output:61eae4cc] %[output:2600e63d] %[output:9482da4d] %[output:394e1214] %[output:67a68fd2] %[output:9f8ff993] %[output:849f1652] %[output:041fcd47] %[output:9597f9c0] %[output:511ca6fa] %[output:18d48aa6] %[output:5f47d5f5] %[output:8c64e6ff] %[output:1c4d4d42] %[output:3bea7d55]
BRW_exp=synchronize(BRW_expS_cal,BRW_expA_cal);
BRW_exp=synchronize(BRW_exp,BRW_expAE_cal2);

clear BRW_expAE_cal2 BRW_expS_cal BRW_expA_cal ;
%%
%compute SSA
%BRW_SSA7w_cal=compute_SSA(BRW_rd_short,names_neph,names_7w3133,lambdaSC,lambdaAE7);
BRW_SSA_cal=compute_SSA_D(BRW_rd_short,names_neph,names_3w);
%BRW_SSA_cal=synchronize(BRW_SSA3w_cal,BRW_SSA7w_cal);
% put all cal together

BRW_cal=synchronize(BRW_exp,BRW_SSA_cal);
BRW_cal.SSAG0_AE=BRW_rd_short.BsG0_S2S20./(BRW_rd_short.BsG0_S2S20+BRW_rd_short.Bacx3_A81);
BRW_cal.SSAG0_aePsapClap=BRW_rd_short.BsG0_S2S20./(BRW_rd_short.BsG0_S2S20+BRW_rd_short.BaG0_homo1);
%BRW_cal=compute_exp_SSA(BRW_rd_short,lambdaSC, lambdaAE); %compute expAE from 3 lamba !!
plotFigControl_cal(BRW_cal,BRW_st.name); %[output:2a5c3baa] %[output:881c6099] %[output:658f5bda]
%%
% Break point for exp: wait answer from Betsy on what to use
names_cal=fieldnames(BRW_cal);
c=startsWith(names_cal,'expS');
names_expS=names_cal(c);
break_BRW_expS_pettitt=change_point_analysis_V3(BRW_cal,names_expS,0.05,'BRW- Pettitt','expS'); %[output:09638136] %[output:3b7a8450]
T_BRW_expS_pettitt=make_table_breakpoints(break_BRW_expS_pettitt,names_expS);

break_BRW_expS_snht=change_point_analysis_SNHT_V1(BRW_cal,names_expS,0.05,'BRW-SNHT','expS');
T_BRW_expS_snht=make_table_breakpoints(break_BRW_expS_snht,names_expS);
%%
%Break point for SSA green
c=startsWith(names_cal,'SSAG');
names_SSA=names_cal(c);
names_SSA_short={'SSAG_3011','SSAG_3012'};
break_BRW_SSAG_pettitt=change_point_analysis_V3(BRW_cal,names_SSA_short,0.05,'BRW- Pettitt','SSAG'); %[output:928ac5a7]
T_BRW_SSAG_pettitt=make_table_breakpoints(break_BRW_SSAG_pettitt,names_SSA_short);

break_BRW_SSAG_snht=change_point_analysis_SNHT_V1(BRW_cal,names_SSA_short,0.05,'BRW-SNHT','SSAG');
T_BRW_SSAG_snht=make_table_breakpoints(break_BRW_SSAG_snht,names_SSA_short);
%%
%Break point for expA 
c=startsWith(names_cal,'expA');
names_expA=names_cal(c);

break_BRW_expA_pettitt=change_point_analysis_V3(BRW_cal,names_expA,0.05,'BRW- Pettitt','expA'); %[output:0df8ef64] %[output:21361001]
T_BRW_expA_pettitt=make_table_breakpoints(break_BRW_SSA_pettitt,names_expA); %[output:91818b42] %[output:8af7b3ff]

break_BRW_expA_snht=change_point_analysis_SNHT_V1(BRW_cal,names_expA,0.05,'BRW-SNHT','expA'); %[output:52adb4e8] %[output:42366ccb]
T_BRW_expA_snht=make_table_breakpoints(break_BRW_expA_snht,names_expA);
%%

%Questions:
%problems in 2019:
%scat: b:no data from 1993-1997, but data for g,r and Q= 720 nm 
%scat: too high values in spring 1995: what happened ?
%scat: greatest "spikes" for radiance than for TSI
% otherwise the green seems to be homogeneous.
%in 2026: 
% BP in June 2007 and December 2012 not considered
% no BP detected at the change between nephelometer 1997/1998

%backscatter: no data before 1998 because MRP has no backscat
%backscat: some issue with the red wavelength since 2017: ratio G/R and B/R are lower.
%Backscat :negatives only between 2005 and 2010 and more low values between 2006 and 2015. Also visible in the ratio between parameters
% otherwise Backscat ok

%scat exponent 2026:
% lower SAE GR in Aug-1993-Nov1997 detected by both Pettitt and SNHT and visual inspection + change of neph in 1997 + period without blue: 
% do the trend analysis with and without  Aug 1993-Nov 1997 and compare the results for the full period and if there is a clear change in the 10y trends.
% lower SAE in Dec 2007-Feb 2010 detected by both Pettitt and SNHT and visual inspection: → similar decrease also happen in 2018-2020 and potentially in 2025. This can  be a natural cause. BP on 2007-2010 should not be considered


%BbsF 2019: not similar seasonality between PM10 and PM1: PM1 seems to have a different seasonality in 2003-2006 and 2008- with longer season with high values.

%2019 only PSAP and CLAP
%abs: strange low values in 2010-2011: bad data ? also visible on SSA.
% abs: the ratio between the wavelength == 1, 1.5 and 2 are favored: this is very strange. Do you have an explanation?
%abs: similarly the exp=0 0.83 1.31 ... 2 as well as -2 are favored.
%abs: for the green, it seems to be homogeneous between PSAP 1w,PSAP 3w and CLAP
%abs: I don't have the AE 31 from 2010 to 2018 (not sure it is necessary).

%2026 PSAP and CLAP: data since 1.1.1998. Missing data from Jan2010 to May 2011
% potential BP April 2015 not considered

% 2026 AE data: 
%%
BRW_tr=BRW_rd_short;
BRW_tr.y=year(BRW_tr.Time);
% begin at the beginning of a year: scat 1978, abs (BaG0_A11) + backscat + SSA: 1998
% end 2025
P=timerange('1978-01-01','2026-01-01');
BRW_tr=BRW_tr(P,:);

% %U >50% for only few cases --> no trend in U and dry, not PM1
names=fieldnames(BRW_tr);
c= startsWith(names,["T1";"T0";"T_";"P_";"P1_";"P0_";"N";"Uu";"U";"Bax"]) | endsWith(names,["dry";"1_A11";"1_S11"]) | contains(names,'Q');
N=names(c);
for i=1:length(N)
    BRW_tr.(N{i})=[];
end
%BsB_S: problem from 5.4.1993 to october 1997: to invalidate
% potential probel for the same period in exp scat GR --> keep and test trend results ?
Px=timerange('1993-04-05','1997-10-01');
BRW_tr.BsB0_S1S10(Px)=NaN;
%one time series for TSP and PM10

%backscat and abs: 1998-2026 
P2=timerange('1997-10-01','1998-01-01');
BRW_tr.BbsG0_S11(P2)=NaN;
%BRW_tr.BbsG1_S11(P2)=NaN;
BRW_tr.BbsB0_S11=[];
BRW_tr.BbsR0_S11=[];
% BRW_tr.BbsB1_S11=[];
% BRW_tr.BbsR1_S11=[];
%Abs:
BRW_tr.BaG0_A11(P2)=NaN;
%BRW_tr.BaG1_A11(P2)=NaN;
BRW_tr.BaB0_A11(P2)=NaN;
%BRW_tr.BaB1_A11(P2)=NaN;
BRW_tr.BaR0_A11(P2)=NaN;
%BRW_tr.BaR1_A11(P2)=NaN;
BRW_tr.BaG0_homo1(P2)=NaN;
%abs: invalidate 25.12.2009 to jan 2011 (end story leak)
P3=timerange('2009-12-25','2011-01-01');
BRW_tr.BaG0_A11(P3)=NaN;
%BRW_tr.BaG1_A11(P3)=NaN;
BRW_tr.BaB0_A11(P3)=NaN;
%BRW_tr.BaB1_A11(P3)=NaN;
BRW_tr.BaR0_A11(P3)=NaN;
%BRW_tr.BaR1_A11(P3)=NaN;
BRW_tr.BaG0_homo1(P3)=NaN;
%%
%compute necessary exp, backscat fraction and SSA: use G and expS bg and
%expA bg and fit

%exp + BbsF
names_short=fieldnames(BRW_rd_short);
names_neph=names_short(startsWith(names_short,{'Bs','Bbs'}));
BRW_expS_cal=compute_exp_D(BRW_rd_short,names_neph,lambdaSC);

%compute abs exp
names_3w=names_short(startsWith(names_short,{'Ba'}) & contains(names_short,{'0_A11'}));
BRW_expA_cal=compute_exp_D(BRW_rd_short,names_3w,lambdaAE); %[output:98450c23]

% compute expA from AE31+33
Ca7=startsWith(names_short,{'Bacx'}) & contains(names_short,{'_A8'});
names_7w3133=names_short(Ca7);
BRW_expAE_cal2=compute_exp_D(BRW_rd_short,names_7w3133,lambdaAE7); %[output:09a2f37b] %[output:2376b229] %[output:0ca296e1] %[output:3965de54] %[output:9b922ae4] %[output:5fa14e70] %[output:6df78528] %[output:6660adae] %[output:41b1540a] %[output:4930ca5f] %[output:47c93f6e] %[output:0364c28a] %[output:9f147c28] %[output:58e6774d] %[output:3e79bf00] %[output:1ee48409] %[output:7c760164] %[output:029bc713] %[output:75f9e04c] %[output:978da8e0] %[output:934940d8] %[output:7f42177f] %[output:774e4d2a] %[output:43dd921f] %[output:624e0fab] %[output:8053a40f] %[output:41331dce] %[output:64fa7b72] %[output:47778218] %[output:74372a75] %[output:96ea315b] %[output:05556658] %[output:35b96170] %[output:5ad1ec71] %[output:52ddc51e] %[output:6bb80cb3] %[output:34ab9ec3] %[output:0a152510] %[output:0115a49d] %[output:91f96e64] %[output:3b22d618] %[output:4793662f] %[output:1f0c69b2] %[output:82ee6987] %[output:6a609989] %[output:45e15f5a] %[output:3e7b5d93] %[output:091562a9] %[output:76de35e1] %[output:76f85a79] %[output:522c7e1d] %[output:986a34ff] %[output:5d14b4b3] %[output:2dfc0f80] %[output:73d22197] %[output:3c760eda] %[output:1e18f3f8]
BRW_exp=synchronize(BRW_expS_cal,BRW_expA_cal);
BRW_exp=synchronize(BRW_exp,BRW_expAE_cal2);
BRW_tr=synchronize(BRW_tr,BRW_exp);
clear BRW_expAE_cal2 BRW_expS_cal BRW_expA_cal BRW_exp;
%compute SSA

BRW_tr.SSAG0_AE=BRW_rd_short.BsG0_S2S20./(BRW_rd_short.BsG0_S2S20+BRW_rd_short.Bacx3_A81);
BRW_tr.SSAG0_aePsapClap=BRW_rd_short.BsG0_S2S20./(BRW_rd_short.BsG0_S2S20+BRW_rd_short.BaG0_homo1);

%%
names_tr=fieldnames(BRW_tr);
% compute neph only on G
%BRW_tr.BsB0_S1S10=[];
%BRW_tr.BsR0_S3S30=[];
%BRW_tr.BbsB0_S1S10=[];
%BRW_tr.BbsR0_S11=[];
%BRW_tr.BbsB0_S11=[];



%compute abs on homo and Bacx3 and SSA trends only on G, the longest time series
BRW_tr.BaB0_A11=[];
BRW_tr.BaG0_A11=[];
BRW_tr.BaR0_A11=[];

BRW_tr.Bacx1_A81=[];
BRW_tr.Bacx2_A81=[];
BRW_tr.Bacx4_A81=[];
BRW_tr.Bacx5_A81=[];
BRW_tr.Bacx6_A81=[];
BRW_tr.Bacx7_A81=[];

%expS only on bg
BRW_tr.expS_br0=[];
BRW_tr.expS_br1=[];
BRW_tr.expS_gr0=[];
BRW_tr.expS_gr1=[];
BRW_tr.expS_bg1=[];

%compute expA exponent only on BG and only since 1.1.2007 (whole year)
Pexp=timerange('2006-01-01','2007-01-01');
BRW_tr.expA_bg0(Pexp)=NaN;
%BRW_tr.expA_bg1(Pexp)=NaN;
BRW_tr.expA_br0=[]; %[output:2cb9f218]
%BRW_tr.expA_br1=[];
BRW_tr.expA_gr0=[];
BRW_tr.expA_bgAE=[];
BRW_tr.expA_brAE=[];
BRW_tr.expA_grAE=[];

% exp S do the trend analysis with and without Aug 1993-Nov 1997 and compare the results for the full period and if there is a clear change in the 10y trends.
%%
% remove BBsF PM1 and other wavelenths
BRW_tr.BbsFb0=[];
BRW_tr.BbsFr0=[];
BRW_tr.BbsFb1=[];
BRW_tr.BbsFg1=[];
BRW_tr.BbsFr1=[];


%%
[BRW_result_MK,BRW_result_LMSlog,BRW_result_LMSlin]=all_trend_STN(BRW_tr,BRW_st); %[output:1a6d931d] %[output:90bf4018] %[output:74accf82] %[output:8f6a210b] %[output:02c66054] %[output:275fae0b] %[output:68a84860] %[output:0e5e2927] %[output:6b329859] %[output:21853660] %[output:0ffe934c] %[output:8a0916d6] %[output:9d73ce37] %[output:633ef791] %[output:3c5fe4a0] %[output:47f9d4f4] %[output:3eb703dd] %[output:5dcd7821] %[output:27d03b35] %[output:1c17ca5d] %[output:767ec57e] %[output:806deac2] %[output:5d1ca422] %[output:0b31e208] %[output:60d29154] %[output:0087bb11] %[output:4ef60291] %[output:83b3a347] %[output:10f020f8] %[output:2aaaaf09] %[output:6b601709] %[output:4002a337] %[output:0da292be] %[output:34c77553] %[output:3f5cfa6a] %[output:2250e661] %[output:71f1d760] %[output:82b16049] %[output:32ca3c35] %[output:90a4cce6] %[output:862eeea8] %[output:78254fee] %[output:5af536c2] %[output:93e2212e] %[output:8d19b7b3] %[output:2856bd3b] %[output:564b1669] %[output:2f3e2c66] %[output:29e927b8] %[output:3813a8dc] %[output:2e2d94fb] %[output:4fb1dd73] %[output:590060cf] %[output:413172ea] %[output:29a3859e] %[output:484a7bec] %[output:6ed9ef45] %[output:1a6e0101] %[output:846b7d01] %[output:000e549e] %[output:1da35d86] %[output:82f20aa6] %[output:3d3b876c] %[output:31a277a7] %[output:8fe07ccf] %[output:40b2d253] %[output:17d1c91c] %[output:6f86febd] %[output:8e0ed078] %[output:7b911553] %[output:298ec69e] %[output:4ee779b1] %[output:248e1962] %[output:83bf22a4] %[output:84060908] %[output:08155abe] %[output:73710e2c] %[output:018dd32d] %[output:44d899d4] %[output:67388040] %[output:22a309cc] %[output:00e838d8] %[output:85d54d6a] %[output:993ae445] %[output:9b31a73d] %[output:36662d89] %[output:01d4841b] %[output:15201ad6] %[output:95d59e56] %[output:5cb5ef92] %[output:4a70d0e1] %[output:86b42ed4] %[output:1fd14807] %[output:34595fc6] %[output:938ffa44] %[output:44402cdf] %[output:01d04c02] %[output:7134ed46] %[output:9039e061] %[output:5ed933b6] %[output:565c79b9] %[output:7e2689a8] %[output:06e6c066] %[output:4171e70b] %[output:494ede43] %[output:7b004c68] %[output:874d8ec0] %[output:2dda597b] %[output:3ada91c6] %[output:128525a5] %[output:024dd3c2] %[output:9becd61c] %[output:348171fb] %[output:340be5b7] %[output:4c73e2a6] %[output:5b500bbc] %[output:93b6cf28] %[output:6c7e4d1e] %[output:937e65fa] %[output:9f12749d] %[output:14fedef7] %[output:7f8f71b8] %[output:43df0a2b] %[output:32d5edea] %[output:1568a178] %[output:393188c0] %[output:0e6bd3c2] %[output:3c75499d] %[output:92f7e111] %[output:32dc71af] %[output:87da0cfe] %[output:46e9006d] %[output:0f549bfe] %[output:5d19fc14] %[output:44484d29] %[output:13089376] %[output:1008e3de] %[output:0dfad22a] %[output:9e8c34ca] %[output:30f45b9b] %[output:51d0b936] %[output:89715877] %[output:0dc7802e] %[output:10733847] %[output:2b710032] %[output:5d5d6928] %[output:3753ca58] %[output:18baf1f7] %[output:8c3bab3d] %[output:65ff9912] %[output:6e9923ff] %[output:48540a2b] %[output:8db3c968] %[output:38ea5b7c] %[output:32e97787] %[output:2ff26752] %[output:42dc50e6] %[output:41376bf1] %[output:907db972] %[output:7981fd14] %[output:0b36dba7] %[output:34d214c2] %[output:5c93efdd] %[output:36b55a29] %[output:5d5bbf8d] %[output:4b16883e] %[output:7cbe1379] %[output:4d89d53b] %[output:9397a265] %[output:1b83538f] %[output:9bd1108c] %[output:93686c5b] %[output:860dd88f] %[output:1c489307] %[output:367decd9] %[output:546b111d] %[output:51142b55] %[output:8213be0b] %[output:90dc90d3] %[output:65f88126] %[output:5c0d4dc4] %[output:66e55553] %[output:3ce0b691] %[output:12129d59] %[output:95063375] %[output:2b869b28] %[output:031c6712] %[output:7cc6f0d7] %[output:0cc5f808] %[output:04817980] %[output:3080e591] %[output:95524bbf] %[output:54710bbd] %[output:9d06bcf8] %[output:077691bc] %[output:29b4bdd7] %[output:32035ee5] %[output:01853600] %[output:658072e6] %[output:2801f0a5] %[output:372b2fea] %[output:5554bccb] %[output:506bae5c] %[output:9126c750] %[output:4b24bfa1] %[output:89a59519] %[output:3494b737] %[output:2d1feeca] %[output:2ff4a3ca] %[output:6762415e] %[output:3e35ec14] %[output:96cdf2f1] %[output:5adb1609] %[output:89bbea40] %[output:8e63cdb7] %[output:585421c9] %[output:8a948b09] %[output:200ea141] %[output:8c92ee02] %[output:982c46b4] %[output:8f8834a4] %[output:704156ac] %[output:6a4c0247] %[output:47f151da] %[output:9d21ca6a] %[output:7fc53de5] %[output:772f5af7] %[output:7b54119c] %[output:51559519] %[output:69e42c04] %[output:3186ce1b] %[output:47547439] %[output:8e6bc093] %[output:9d134474] %[output:47231f80] %[output:2c0a50bc] %[output:078c1650] %[output:9d65e458] %[output:2c761628] %[output:62f1bbf8] %[output:59148269] %[output:9dce30cb] %[output:78ed7ead] %[output:4c07a296] %[output:06dc8a2f] %[output:445548e3] %[output:19eeab58] %[output:40debe50] %[output:544a3c5c] %[output:75ffc104] %[output:37c3ae91] %[output:67d65a07] %[output:196d887a] %[output:9809bce7] %[output:79fc8e81] %[output:34c31025] %[output:36558d63] %[output:4c536d58] %[output:2f1bc288] %[output:5f4a0cfc] %[output:3164b319] %[output:37e44cc3] %[output:062bea41] %[output:8a2e83b3] %[output:596aaa7e] %[output:05ede065] %[output:8d42046a] %[output:4ec52aa6] %[output:8d511a73] %[output:997613bb] %[output:6688afbb] %[output:4135bf82] %[output:26c9cd1d] %[output:2895fcf6] %[output:62d2836f] %[output:0ee813cc] %[output:142b0a16] %[output:67f38a16] %[output:34ac7265] %[output:4e51d712] %[output:64fe1588] %[output:3c2330cf] %[output:2a4c696d] %[output:6b30af8a] %[output:4fbb3e67] %[output:04e78057] %[output:4b6c3002] %[output:3024b889] %[output:43a1badd] %[output:481b0c93] %[output:3f3232c5] %[output:4cb41de4] %[output:213e149e] %[output:99a27963] %[output:59bca9f1] %[output:1edc1436] %[output:92e45a27] %[output:89a9978e] %[output:43a98f96] %[output:711bbcda] %[output:79e15497] %[output:19adfb2f] %[output:5637d336] %[output:7e671f54] %[output:12cdf73d] %[output:3bbddfad] %[output:76b100fe] %[output:390ccf2a] %[output:4469a233] %[output:30d4e6d7] %[output:713edae2] %[output:6dd81bf8] %[output:141220f3] %[output:06416cbf] %[output:3b25ddc8] %[output:313fa6c8] %[output:912241b9] %[output:14d949d6] %[output:64b1e93d] %[output:931bd380] %[output:4ad24684] %[output:233226df] %[output:9df1d2c2] %[output:9ab07325] %[output:216a19f9] %[output:7f6044e1] %[output:4a1f7854] %[output:8e9ee41e] %[output:7992a999] %[output:7b8ce3f5] %[output:8d883bbc] %[output:84193af1] %[output:297ca220] %[output:17703273] %[output:6d2fb422] %[output:831976a5] %[output:113f7888] %[output:233d1d1b] %[output:59c6040d] %[output:5b3e15e8] %[output:7400dfd6] %[output:801873b6] %[output:5c3fb85d] %[output:0ee00466] %[output:661c10e9] %[output:278dc67f] %[output:612d2e14] %[output:8a31ff0f] %[output:137e2f57] %[output:96b7b375] %[output:175c0342] %[output:7bf2a46a] %[output:325839f9] %[output:03ac8eff] %[output:9c6667b5] %[output:795d3d1f] %[output:22b85d3a] %[output:918473fb] %[output:11001239] %[output:984e4182] %[output:639ee66e] %[output:8b67ce71] %[output:9e78f5c5] %[output:56059f9f] %[output:87e9be89] %[output:6cb8b6ee] %[output:1d6c85c3] %[output:0dfe9da0] %[output:742b8293] %[output:7141f281] %[output:4afc3452] %[output:2623bc39] %[output:7fc0485c] %[output:4fb6c006] %[output:8561cd64] %[output:20d39efb] %[output:11431d9f] %[output:048e97ef] %[output:1852498f] %[output:9143a6f4] %[output:3380886d] %[output:47ffb32b] %[output:8588c821] %[output:135e38ff] %[output:454c7f7d] %[output:879c56fe] %[output:111bca2d] %[output:751c5705] %[output:5208834f] %[output:76b6c351] %[output:8bb63f0f] %[output:94b3a58d] %[output:8ce8b5ab] %[output:4bf4381a] %[output:79d6160f] %[output:969ef301] %[output:3773d76c] %[output:6dd507d9] %[output:9b0c1f72] %[output:2efecd0c] %[output:7599f9c9] %[output:0b9cde16] %[output:1cf8bdfa] %[output:8a595b1b] %[output:7c96b513] %[output:42879679] %[output:18d08194] %[output:0ec79177] %[output:27737cda] %[output:3b8b3559] %[output:7d2f958a] %[output:93559b24] %[output:470c9215] %[output:59c56c76] %[output:5237cd0e] %[output:508a7a92] %[output:4255bcda] %[output:7ec81959] %[output:0ade710f] %[output:5443e0d6] %[output:6c9e2418] %[output:875fd8b2] %[output:8b235608] %[output:3c8a13af] %[output:3e48f01d] %[output:3051d908] %[output:6b4dbad1] %[output:5644bfec] %[output:6c66ba1a] %[output:153e65a5] %[output:6e50782e] %[output:95bd2060] %[output:357ff63e] %[output:00b81aeb] %[output:12b8c40f] %[output:01d1f83d] %[output:33d2d9e1] %[output:6db989f8] %[output:25e66804] %[output:457f4c1a] %[output:7d574454] %[output:0db7d784] %[output:454b51c7] %[output:5c73b87f] %[output:39652ade] %[output:74357a07] %[output:5b2412d9] %[output:0903c4a8] %[output:234929f7] %[output:070a795f] %[output:99479e13] %[output:103d094a] %[output:0404a91a] %[output:4fa9468e] %[output:0bde47ba] %[output:85260890] %[output:10db9db0] %[output:4cdc6a3c] %[output:81963189] %[output:4b4c6db2] %[output:64f5ee8e] %[output:13139e4b] %[output:16cbb4eb] %[output:21ba8c85] %[output:3d0fb04c] %[output:8f126aca] %[output:8a2da152] %[output:311a487c] %[output:2170f141] %[output:6bddefc7] %[output:0da5f4ac] %[output:0817b4ed] %[output:0ee7b693] %[output:7a52dac7] %[output:7c16d95e] %[output:5195e48f] %[output:15c2e6a7] %[output:5d9bd20f] %[output:198efaba] %[output:80284d00] %[output:72e2ba8b] %[output:1835e793] %[output:4531b12e] %[output:1c381965] %[output:224b96be] %[output:635c9b22] %[output:445df048] %[output:6c09b771] %[output:80cf8149] %[output:3daa6ac9] %[output:90bb9036] %[output:758fd6f2] %[output:753ef725] %[output:7094ccf5] %[output:1ac20262] %[output:2169675d] %[output:187dec53] %[output:38b46c11] %[output:5c1cb559] %[output:7350707f] %[output:66461452] %[output:7cc48811] %[output:5b14209d] %[output:5b1610ef] %[output:112c2587] %[output:4de3aadc] %[output:008a21ee] %[output:43ad4fa1] %[output:4edaac14] %[output:8b903865] %[output:4d7cd259] %[output:15800c7c] %[output:6e9948b2] %[output:6f14ceca] %[output:883a501b] %[output:3c12e34d] %[output:7cb08fa9] %[output:469bc765] %[output:34071e90] %[output:0fadc2a5] %[output:97f618c7] %[output:1e270ebf] %[output:15b7c972] %[output:1d29dfc8] %[output:6a45d3fe] %[output:7c48b797] %[output:7a26972b] %[output:98a0eeee] %[output:1e52ecfd] %[output:34858ae1] %[output:6a5ad6af] %[output:1ccf8032] %[output:29198aa4] %[output:453264d7] %[output:6652ff6b] %[output:663158e8] %[output:833d3a1e] %[output:52afc33b] %[output:9965b3a1] %[output:65f3b8c0] %[output:4b734d75] %[output:8aa73ebc] %[output:99527e6a] %[output:8c00a7e1] %[output:82bf6ca1] %[output:83e751ea] %[output:959448c9] %[output:6046acdb] %[output:6d234253] %[output:2ab8be3f] %[output:9cfffedf] %[output:91fbd516] %[output:1ede8235] %[output:697cf7a9] %[output:547b319e] %[output:18430d68] %[output:16b30a3e] %[output:941a5882] %[output:3380b4bb] %[output:788649f7] %[output:1f54f3f4] %[output:06b5f3b9] %[output:2b2ef57a] %[output:1d6f898f] %[output:213aae46] %[output:22932d90] %[output:8f82340c] %[output:2f29075a] %[output:6b48a3f5] %[output:9c4369be] %[output:11a0a64f] %[output:5157b02b] %[output:9d6b948c] %[output:25fe6413] %[output:8386dd71] %[output:413d058a] %[output:3550837e] %[output:75fa8c8d] %[output:80ebbb5f] %[output:1aa4d602] %[output:12ea44e3] %[output:9a332fd1] %[output:091957ad] %[output:04d88352] %[output:2c6fb989] %[output:9f5e9c1d] %[output:72da9bba] %[output:1a9ba36d] %[output:42f90fcb] %[output:1e8345eb] %[output:15513712] %[output:8e53c04b] %[output:4c92f5a4] %[output:30b4837a] %[output:597192b1] %[output:73cfa7ea] %[output:14a98fd7] %[output:19e6b3fd] %[output:70003f97] %[output:06c51ea5] %[output:6c8bd1a7] %[output:1a0bcace] %[output:90364ea7] %[output:0251823f] %[output:40bdd972] %[output:7a040431] %[output:6acadad0] %[output:63115f0c] %[output:0801b24e] %[output:6d34fbf6] %[output:472e405b] %[output:31966b1b] %[output:5acb435e] %[output:7dea8a02] %[output:62a11b41] %[output:8f4f19a8] %[output:0608e2c1] %[output:4a8413df] %[output:60e8ca5b] %[output:9dd9441c] %[output:9dbcc9ad] %[output:5dfab478] %[output:582e4a1e] %[output:92d880f0] %[output:7d6d374f] %[output:948cde65] %[output:81af6044] %[output:8634c476] %[output:2b693640] %[output:3bafc0c4] %[output:971d3b0c] %[output:9bbea782] %[output:0aa815a4] %[output:623202dc] %[output:9e71f5e3] %[output:7dafc6ec] %[output:456551ba] %[output:38a604db] %[output:960579c1] %[output:5f59913f] %[output:6cc418b9] %[output:0c4e9014] %[output:21cc769b] %[output:9f24076b] %[output:636f52e4] %[output:0d7ea755] %[output:5e89a2d3] %[output:2483172f] %[output:88df6649] %[output:6334f9b0] %[output:330ad32a] %[output:2dc7fc8e] %[output:21c0df7a] %[output:8bb5f7ea] %[output:3051a516] %[output:336c33ea] %[output:0a085ee4] %[output:3c466ec0] %[output:015cb391] %[output:4b57dd35] %[output:18ac3c72] %[output:0da4bf13] %[output:30acf795] %[output:8ac12db5] %[output:3e21333c] %[output:19130dd3] %[output:0ac06016] %[output:0025699a] %[output:4c969a4b] %[output:3205fcd6] %[output:721da8fd] %[output:8d366bc5] %[output:983aa819] %[output:3f0b16eb] %[output:88db5cc2] %[output:6a51fa74] %[output:75d9574e] %[output:7f438b82] %[output:0f331783] %[output:9dc5979e] %[output:235cad5d] %[output:58b620ec] %[output:7b41c513] %[output:2ddb25f2] %[output:5ab2fc80] %[output:2f71bef0] %[output:37cf23e7] %[output:295c7b54] %[output:61077f08] %[output:2e4725a3] %[output:50878f63] %[output:7abb99b8] %[output:4e23b6c7] %[output:4278e853] %[output:7c936736] %[output:26a8f37b] %[output:8afe35af] %[output:6f7c22ed] %[output:50dc6dab] %[output:4d6c376d] %[output:3d9c5eab] %[output:49c72875] %[output:84ee9072] %[output:218cf845] %[output:2d1a72fb] %[output:8ad5b377] %[output:48e441d5] %[output:365b2f9b] %[output:36d3706e] %[output:18b338aa] %[output:5cf97442] %[output:4863cf71] %[output:2c592f84] %[output:34dd8e79] %[output:586a6115] %[output:488fb4b4] %[output:87fce924] %[output:0b207c75] %[output:3675ed46] %[output:43170eda] %[output:3251e38d] %[output:66d74be0] %[output:68b65577] %[output:7384198f] %[output:5ce446de] %[output:7698e34c] %[output:4474b738] %[output:191152a0] %[output:122b0f3b] %[output:5b03e27b] %[output:3e4b2265] %[output:2e2ff7bd] %[output:3965d6bd] %[output:0ffb08e2] %[output:2c147e00] %[output:34fde0c9] %[output:21ae4975] %[output:10b9104d] %[output:6566c6a0] %[output:5f904590] %[output:79bf07a8] %[output:0186dc17] %[output:363c88f7] %[output:36b180fc] %[output:58d87de9] %[output:8387a974] %[output:52c9e2c7] %[output:1831ac8b] %[output:66b7fcea] %[output:61607279] %[output:682774aa] %[output:28ced20c] %[output:752a74f9] %[output:6b516c0a] %[output:4094f7b9] %[output:506ff6e7] %[output:9336c533] %[output:8a8a2ede] %[output:967d09b3] %[output:8d5180de] %[output:74344efb] %[output:21aaddf9] %[output:71e11d2c] %[output:2d10f7e3] %[output:54be5c90] %[output:660ef67e] %[output:3398924f] %[output:169784af] %[output:56970fba] %[output:06b0e977] %[output:8a77ba81] %[output:4722ddaf] %[output:4e6de565] %[output:6a74debd] %[output:842a63d5] %[output:1d183d89] %[output:40de2c29] %[output:475024df] %[output:3a5a4286] %[output:76076ca3] %[output:8912acf4] %[output:0aefa498] %[output:2532de9b] %[output:6fa36ea3] %[output:8ceeb71c] %[output:80b37eb0] %[output:9417dabe] %[output:5a8c7a5d] %[output:87b7195e] %[output:8b698fa3] %[output:0dfda2b9] %[output:6c5ac0d7] %[output:1c185aa0] %[output:053ef558] %[output:7f12dde5] %[output:2e288b62] %[output:2edbd576] %[output:714aef69] %[output:49946398] %[output:7b78cc28] %[output:803b4b5d] %[output:7a4a7201] %[output:9327ad9a] %[output:1dabecb7] %[output:941aa8b6] %[output:93e78f5d] %[output:1448401d] %[output:49aab42f] %[output:4664b6ae] %[output:387f045f] %[output:24371ede] %[output:6cfb6040] %[output:3537b91a] %[output:0dcc4ea3] %[output:06638abe] %[output:88c538af] %[output:2462298d] %[output:4ef77851] %[output:8d5ffc95] %[output:0d81dbe2] %[output:0b080ef5] %[output:51290129] %[output:58623409] %[output:47d2303a] %[output:7c9fb410] %[output:2e7ed239] %[output:248e63cc] %[output:8638c5a2] %[output:6b284375] %[output:64c0ff0e] %[output:594ca13c] %[output:624d246b] %[output:2b0f3cfa] %[output:89588526] %[output:2505a68f] %[output:6e72b5f4] %[output:2f950a37] %[output:49dc6a9d] %[output:52f25e5a] %[output:92e68211] %[output:4b448e58] %[output:9f7c6739] %[output:7443477a] %[output:8027f5df] %[output:0a2a7959] %[output:9c60dee2] %[output:91c3db57] %[output:4f71c328] %[output:20a37044] %[output:4d127a48] %[output:4a03cd28] %[output:5aae1920] %[output:77fb2d40] %[output:18e5ffad] %[output:39dfb265] %[output:5caeece6] %[output:2c8dc615] %[output:07fd2a8d] %[output:6d4e90cc] %[output:0a545acc] %[output:8f37d48b] %[output:9f6022a0] %[output:275a974c] %[output:1f36175b] %[output:5906e80b] %[output:887499bc] %[output:1c4aea40] %[output:2793d475] %[output:1edb9035] %[output:6c7ae58d] %[output:3eb66204] %[output:0e363b82] %[output:6e92393d] %[output:11ccbdfa] %[output:587529f3] %[output:570e88e7] %[output:3900e64a] %[output:3575c0f9] %[output:2ad760b2] %[output:38cf78ae] %[output:0225c240] %[output:6a0488ba] %[output:4698c34d] %[output:89f194b4] %[output:7b7682da] %[output:31754962] %[output:2fbe3fed] %[output:2e2db140] %[output:12f9ea2b] %[output:16bbe101] %[output:31342d67] %[output:8e099250] %[output:82f43d66] %[output:68ac73ea] %[output:34721f8a] %[output:5942e1fe] %[output:006fe226] %[output:22545626] %[output:48cd084e] %[output:585bfd8c] %[output:1f6aef5e] %[output:2e2d6192] %[output:7557f411] %[output:95b6eee4] %[output:334cb600] %[output:20d579e6] %[output:254a64c0] %[output:13e374fe] %[output:4dd2fc4d] %[output:0855a80a] %[output:89775020] %[output:4faf8252] %[output:792b81d4] %[output:7155133c] %[output:58181a86] %[output:3c5bc103] %[output:4e1cfb09] %[output:1924c3bb] %[output:57f5f2e6] %[output:95d5428b] %[output:8bf89b19] %[output:4ddc9a86] %[output:6175de95] %[output:39fe8ddb] %[output:27c2871d] %[output:1837d5e0] %[output:7f25e7ac] %[output:0c2d92df] %[output:554cc82c] %[output:127911b7] %[output:29ac48cd] %[output:574e224b] %[output:45a4b885] %[output:7de238d1] %[output:4c7678b3] %[output:63c95266] %[output:2e8866ad] %[output:259a3e07] %[output:99abc3bc] %[output:87eaa7ef] %[output:85a3b8df] %[output:292ea369] %[output:92156ff3] %[output:2706e187] %[output:9fc54de2] %[output:597c3bd7] %[output:5ce550ef] %[output:2b2bdb07] %[output:5fb95bb7] %[output:1c1f91d3] %[output:88ecc3d7] %[output:6ce147c2] %[output:5f2f7ac8] %[output:93fe0ded] %[output:3b9f216d] %[output:5f667f47] %[output:5e2fa66c] %[output:2c58805a] %[output:83b4a87b] %[output:17d0ba2a] %[output:3f1f6d78] %[output:337611aa] %[output:66b4883d] %[output:8721a920] %[output:1d1581b1] %[output:60705e91] %[output:1cf2a9b5] %[output:516360a2] %[output:938e7d01] %[output:59b8d9ba] %[output:80e8f246] %[output:1950cda8] %[output:9a480775] %[output:976ba4c2] %[output:571549bb] %[output:4d5ebce1] %[output:822f0fdf] %[output:7bfd20d4] %[output:82c8b830] %[output:6a6f7ea2] %[output:34e17c9e] %[output:4ae1df08] %[output:8612d969] %[output:507d9e89] %[output:27533850] %[output:42aff600] %[output:2808f36b] %[output:1c1fd27c] %[output:0c83caf3] %[output:5be3ca54] %[output:88fda62e] %[output:77fb56cf] %[output:25ae451a] %[output:142011d3] %[output:149100a6] %[output:37a85aad] %[output:84cfd9d8] %[output:09cee69c] %[output:8792c568] %[output:9933e756] %[output:419746f4] %[output:79f87bf1] %[output:31b9d295] %[output:73151e9b] %[output:3969725e] %[output:6334ea98] %[output:00fb1fa9] %[output:4b8fcaa5] %[output:8b3d7c14] %[output:558d621c] %[output:62c51fd5] %[output:149ad87d] %[output:183cdfc6] %[output:73e4db2b] %[output:01a55b95] %[output:035add23] %[output:693605fe] %[output:8cdc2a0a] %[output:0ca93415] %[output:9553fce4] %[output:614b130b] %[output:3a5413b9] %[output:16f7fe47] %[output:3820925b] %[output:85d90bd8] %[output:9b26cc5a] %[output:75b252d4] %[output:5a30f061] %[output:84ad7612] %[output:2ebcbeb4] %[output:5b236c5b] %[output:1903177f] %[output:8b258047] %[output:12909183] %[output:49782fa3] %[output:49dc514c] %[output:26974022] %[output:16e80540] %[output:9213149d] %[output:16882105] %[output:49702400] %[output:93232c1a] %[output:20ffa0a1] %[output:16f84d9f] %[output:2cbbba97] %[output:20624876] %[output:47528640] %[output:44322dd4] %[output:52a7cbae] %[output:63b0d367] %[output:49be8e2b] %[output:2af5650f] %[output:9451822d] %[output:6ae73ee2] %[output:4fcee8dd] %[output:2a580ca8] %[output:393c3231] %[output:643bc4b2] %[output:451b5474] %[output:1e391f2b] %[output:29a99fc7] %[output:1ba3fb7a] %[output:91c5b0cb] %[output:1444427a] %[output:47034e22] %[output:77cce9d3] %[output:367bb50a] %[output:2bcf67b7] %[output:3a60695f] %[output:2c697b5a] %[output:6f3bbe72] %[output:21f632a2] %[output:13c77413] %[output:2f377d74] %[output:4be6877c] %[output:2bc77e66] %[output:9325c816] %[output:4459f20c] %[output:19803d9b] %[output:3aea5ef7] %[output:9d17b0a0] %[output:32b06ea8] %[output:831bf892] %[output:2b31a37d] %[output:2cdb3459] %[output:136f79ed] %[output:5490e1f0] %[output:797550e0] %[output:64978c26] %[output:54c675cf] %[output:46662cc7] %[output:457a8a73] %[output:3384717c] %[output:7c60e759] %[output:4ac55cb8] %[output:35d599c5] %[output:43bf38ba] %[output:770bc4c8] %[output:1f6ca1e1] %[output:1b2cfd90] %[output:1d2abde9] %[output:93151d9e] %[output:297189c7] %[output:134a2def] %[output:306ba9c8] %[output:65d6e7fb] %[output:03818a92] %[output:43c3bd7e] %[output:20ac780d] %[output:64547591] %[output:0828ae10] %[output:5a559b7e] %[output:9956869a] %[output:8e6ed4f4] %[output:645093b6] %[output:8e2234fc] %[output:06b92d0b] %[output:9047c90e] %[output:1fe031c4] %[output:1d2cbf9d] %[output:10e22b02] %[output:8e8e61eb] %[output:47b6c75e] %[output:3b262e29] %[output:54bc881e] %[output:64ead7e5] %[output:3de75cba] %[output:339585c1] %[output:5d38f3af] %[output:6d6ffd01] %[output:80e4e314] %[output:32b644ae] %[output:1d97f9c7] %[output:2941dc63] %[output:9c937f8c] %[output:0d577d2a] %[output:81d5fd03] %[output:98c2f4af] %[output:1093bee7] %[output:5f2679e1] %[output:74ec2644] %[output:69d38380] %[output:5b941da0] %[output:4f788cb6] %[output:4f4e71b4] %[output:12f2f3c8] %[output:9fce8c26] %[output:8e8acdd1] %[output:59e15eb0] %[output:62b3010b] %[output:5eabad07] %[output:4d2ce568] %[output:75696798] %[output:6afb3625] %[output:08d722eb] %[output:8c894d46] %[output:8429d51f] %[output:64b8861a] %[output:762dd105] %[output:4370906a] %[output:65bdfde5] %[output:1284b270] %[output:049bce27] %[output:5ec7eacd] %[output:2bd5c28d] %[output:64cf91d0] %[output:77bd1596] %[output:95c4cf8e] %[output:6e10772d] %[output:36cb4271] %[output:4432ddb3] %[output:318b240f] %[output:7c05cced] %[output:72691f80] %[output:29f350e8] %[output:6cccb46a] %[output:87053d98] %[output:90f70cbb] %[output:016de8a4] %[output:9be0c838] %[output:180fe098] %[output:9dc2d412] %[output:2d7488e6] %[output:64fd04d8] %[output:6adcddaa] %[output:3a65b52d] %[output:415f2770] %[output:262d8efe] %[output:161b01c6] %[output:2ca79f1f] %[output:1ecb8272] %[output:3806cca5] %[output:837753d3] %[output:9caa9acd] %[output:4b4fd6d8] %[output:3d115eb1] %[output:3bcb2e2d] %[output:698906f8] %[output:20d967a2] %[output:6b64e99d] %[output:134338bf] %[output:66d8092b] %[output:3ddad206] %[output:92a6b2b1] %[output:1be50e99] %[output:15341366] %[output:41c7f270] %[output:1cc42e2d] %[output:53ff2ee0] %[output:40f8d912] %[output:133a7222] %[output:50845773] %[output:34201cd3] %[output:07626141] %[output:600aaae9] %[output:0b18bfcb] %[output:4b98efdd] %[output:2908d658] %[output:7fbedf51] %[output:95794977] %[output:665a6e71] %[output:9a1b1d2e] %[output:216cc764] %[output:292d2644] %[output:2b697c46] %[output:31b41ca5] %[output:0506da1d] %[output:376c811a] %[output:488acef2] %[output:5cc8013d] %[output:0ff8c225] %[output:0724480f] %[output:235421d5] %[output:74f91fec] %[output:1c08d450] %[output:09afa47c] %[output:27b24d5f] %[output:3706599c] %[output:63019134] %[output:4b482c6e] %[output:47388617] %[output:4e851a44] %[output:27e50fef] %[output:448b63c3] %[output:4d63c7e0] %[output:08aa47d4] %[output:6d77c0d0] %[output:5db06ee2] %[output:59c85068] %[output:1c0712a6] %[output:0d629e19] %[output:6a6f634c] %[output:4132c1bb] %[output:45f9817a] %[output:640f8f1a] %[output:8ca598ed] %[output:759f3d8a] %[output:82089c30] %[output:5a3ee51e] %[output:7874c0ef] %[output:6a3aaa28] %[output:41fc1ba2] %[output:5b4f39bb] %[output:04f7f0cc] %[output:1e4e55ab] %[output:973717c7] %[output:4dec8759] %[output:8e9efc3b] %[output:6904849e] %[output:4946c2e5] %[output:8c6051ca] %[output:0951f871] %[output:15164b74] %[output:5a076e2f] %[output:020d1bbe] %[output:8fa77a60] %[output:96472f79] %[output:999ea3ab] %[output:025f8052] %[output:0b841185] %[output:23431ee5] %[output:1f1f4fc8] %[output:698074dc] %[output:4bf87e74] %[output:1623a8bf] %[output:31664ea2] %[output:26dd6d91] %[output:4b82ecf0] %[output:7aa0429e] %[output:33414c29] %[output:9ea249f3] %[output:1e15d10c] %[output:3fdedfa2] %[output:4a1dd1d8] %[output:8b95b246] %[output:14c2dd07] %[output:0876fae4] %[output:08d1e8d6] %[output:415c56cd] %[output:5708f700] %[output:86d52a14] %[output:3675ed81] %[output:28ce86f7] %[output:409b82bb] %[output:59b5e8bb] %[output:7570c643] %[output:69a79924] %[output:9ffd4381] %[output:39cb4046] %[output:7e5b1f12] %[output:46078d26] %[output:32350b6c] %[output:1126984d] %[output:611eb2e6] %[output:4bbdc540] %[output:5207df31] %[output:30055fa1] %[output:2373c263] %[output:3eb2ea88] %[output:45cdad1f] %[output:0af0c331] %[output:311128d8] %[output:08004ff4] %[output:742655d8] %[output:9e22d5b4] %[output:6044389a] %[output:4a20a12b] %[output:4de3e458] %[output:0975e141] %[output:6780b77f] %[output:343a86fd] %[output:82787bfe] %[output:7d4941d1] %[output:9d6748ca] %[output:33d34aaa] %[output:614a4538] %[output:55a82622] %[output:3f0a7ced] %[output:6fad765a] %[output:77ff3e93] %[output:3b03efe3] %[output:3646c97c] %[output:80b8fc47] %[output:78ab6f4b] %[output:08e32479] %[output:1e66b98f] %[output:38190043] %[output:38bb8c11] %[output:7811609c] %[output:82c8d038] %[output:8071b92e] %[output:8861bee4] %[output:03fda433] %[output:0861bfb8] %[output:09063253] %[output:38721af6] %[output:3878f7e8] %[output:61ea03b5] %[output:815ed6c6] %[output:41e2b519] %[output:49867fe3] %[output:3ac276c6] %[output:3f7ed5f3] %[output:0f01b1d9] %[output:51ead3b8] %[output:6954b36c] %[output:5eacbc6b] %[output:160eb0ef] %[output:70e85d0b] %[output:11acf8d2] %[output:4af6e2ee] %[output:4e997625] %[output:57109d9f] %[output:9014a0eb] %[output:5bbfa668] %[output:57aa0885] %[output:4823c304] %[output:4ce2860f] %[output:112a2c68] %[output:8e622b17] %[output:8138c763] %[output:98bf3f39] %[output:43a0ec4d] %[output:74bf838c] %[output:905c740a] %[output:158387e2] %[output:536df220] %[output:300f1eed] %[output:0c4ef80f] %[output:9304604c] %[output:196a9d10] %[output:46792483] %[output:37dba0d6] %[output:5fba7732] %[output:931b9515] %[output:734fb0ff] %[output:430e41ef] %[output:4ae5c49d] %[output:05a948ca] %[output:2beff939] %[output:1494941a] %[output:6d1544b8] %[output:4e83da6a] %[output:36736ebc] %[output:0ff5d822] %[output:3001f108] %[output:6eb06652] %[output:06160112] %[output:95c57c56] %[output:8c336f53] %[output:84290eb6] %[output:961d2d14] %[output:4bf02bcf] %[output:130d1473] %[output:5a28f104] %[output:17f8ba08] %[output:72771ad5] %[output:04da8cdc] %[output:4792d246] %[output:1379829e] %[output:2e096717] %[output:06615686] %[output:96b5604a] %[output:3a00294e] %[output:2be47f69] %[output:419be82b] %[output:4941d267] %[output:25fad3b3] %[output:1fe1a0ff] %[output:6ecd3117] %[output:31278e7c] %[output:976d4914] %[output:4a0aa4e3] %[output:463941c3] %[output:465729e7] %[output:3a654664] %[output:8368355e] %[output:9f8e31b5] %[output:97c13263] %[output:20ed1d48] %[output:41b33526] %[output:098163de] %[output:35f2d269] %[output:092ad1d4] %[output:67de3cc1] %[output:67c33165] %[output:920781d3] %[output:751c01f8] %[output:767fba1a] %[output:6d5cf8e8] %[output:88b98ccb] %[output:51acdd25] %[output:51f52304] %[output:4df58246] %[output:667645a4] %[output:333395fb] %[output:3cde3f37] %[output:02c0e161] %[output:2d18ebc1] %[output:2b78fa34] %[output:5846c457] %[output:5d317eb2] %[output:3da0030d] %[output:8f342e70] %[output:5921627a] %[output:2bb36146] %[output:74a00d3a] %[output:010de11b] %[output:1ca71d03] %[output:2ff7198c] %[output:3dbdefee] %[output:0ba9db33] %[output:82a4b16c] %[output:031ad614] %[output:8f755101] %[output:1c0bd465] %[output:8ee9e4a4] %[output:1b2792a1] %[output:2755528e] %[output:10f13a7d] %[output:50cf9a75] %[output:15e7d571] %[output:381d2c15] %[output:3ced0cd2] %[output:3034305f] %[output:6c9d9f9d] %[output:426c7b1f] %[output:7c9dfbc0] %[output:190b393c] %[output:74db290e] %[output:51974500] %[output:714530bf] %[output:6e4e60ed] %[output:15899809] %[output:333aefcc] %[output:71443a11] %[output:683b1ab2] %[output:4d02dcc1] %[output:5351b481] %[output:3deb909c] %[output:7d103324] %[output:1ae545aa] %[output:866198af] %[output:70023555] %[output:6404b111] %[output:8a6f2990] %[output:7bd9e406] %[output:0468ee64] %[output:4ec439ed] %[output:97afa3a5] %[output:0f886cb8] %[output:2ab8675c] %[output:21cc4258] %[output:8d4b1693] %[output:9700c28a] %[output:37ed52b2] %[output:6cd0d63d] %[output:06801188] %[output:8a77a511] %[output:76bb44b7] %[output:398ac965] %[output:55574715] %[output:5dc3be4b] %[output:6cbed04f] %[output:1c2d0c14] %[output:0858d5bf] %[output:1957b074] %[output:00160ca7] %[output:2231c7a9] %[output:6349738e] %[output:26d81c93] %[output:25138e8c] %[output:7d7e0380] %[output:7a5c39ca] %[output:453199f1] %[output:51dfbb4a] %[output:74f8b9e3] %[output:862b215d] %[output:8e3e4e8b] %[output:06bc1bb1] %[output:794f529f] %[output:0f1a0b54] %[output:3b33ae72] %[output:2dac6c13] %[output:45f25bc9] %[output:37771f7c] %[output:8d42183d] %[output:5c1f4605] %[output:82375a87] %[output:9779c5d3] %[output:4901e2f8] %[output:99d4fe90]

writetable(BRW_result_MK,'BRW_res_MK.txt'); %, 'delimiter',',' )
writetable(BRW_result_LMSlog,'BRW_res_LMSlog.txt'); 
writetable(BRW_result_LMSlin,'BRW_res_LMSlin.txt'); 
plot_10y_in_two(BRW_result_MK, BRW_st,'y'); %[output:2b5e2d53] %[output:8dcac29b]
%%


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":30.3}
%---
%[output:20e02795]
%   data: {"dataType":"text","outputData":{"text":"Lecture du fichier en cours...\nLecture terminée : 438312 lignes x 192 variables.\nValeurs NaN      : 52833570 (62.8 % des cellules).\nPériode couverte : 01-Jan-1976  -->  31-Dec-2025 23:00:00\n","truncated":false}}
%---
%[output:0d6482bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: there is a lot of days with less than 50% data coverage"}}
%---
%[output:1343ce23]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAIoAAABTCAYAAABAtp+DAAAAAXNSR0IArs4c6QAACJJJREFUeF7tnTFoHEcUQL8TkUSFi0gkMdLFvhRWk0KksiwCcgIur5EN5ztCGheGBOJCtqwTB7IbWbY5kA0JuFCTwkKFU8RNwJDEJIVVOSKkUiBHLF+cGCuksnFkFP44\/\/i32tmdmdvd2139hcNiPTP\/z\/9v\/vyd037t2d7e3ga5xAIhFtgjoAgjJhbYAcrKygrUajXVd3R0FJaWlmBgYADu3r0L1WpV3b948SKUy2WT8aVNTizQAcr6+jrMz89Do9FQcFy6dAlarRZMT09DvV6H2dlZNW3eJid2kGl0s\/VgFFleXobJyUlYXFxU0aW\/vx9mZmagUqnA2NiYGHiXWCAwR8GIUiwW4cCBAwqYhYUFZRYEZXx8vL39bGxsAH7kyoYFCoUC4Mfm0oKCuUqz2YRz586p\/EQHCgJy9uxZWF1dtZErbXtogUOHDsGVK1esYPEFhSIJJawICt7z23ooyUXBw8PDavoPHjxQEQYVsrlc+\/VCZlZ1xQV99epVuHHjhlXqsAMUBGJiYqJjkM3NTZiamvJNZgkUW8E2AEnb6Czg6q8OUPgjMKlWKpVUbrK2ttZ+POZQuAqObuoyko0FXP3V9YGbq2CbyUnb6Czg6i8BJTofZGIkASUTbuq9kgJK732QCQ0ElEy4qfdKCii990EmNBBQDN3UOv8BbP3VhL43izB0\/jvDXvlpJqAY+vL3T96BrUdN6HujCPu\/+M2wV36aCSiGvqSIgs0FFPNv\/+UcxRCwvDSTiJIXT8Y8DwElZgPnZXgBJS+ejHkeAkrMBs7L8AJKXjwZ8zwElJgNnJfhBZS8eDLmeQgoMRs4L8MLKHnxZMzzEFBiNnBehhdQ8uLJmOchoMRs4LwML6DkxZMxzyNyUPgrpag7f+eHl71wFRyzPWR4jQVc\/aV9pfT69etw6tQp9e6xvCmYH+4iAwUjCVYvuHPnjrIOvaRu8u7xsWPH8mPRnM4kMlDIPggGB0VXzYAEnz59GvAjV7otgC+oR\/KSuiso8pJ6ugEh7W7evKnKlNj6S\/urkN6IErb12ArOhlnzp2WsW48ks\/kBJlZQvI\/HUvYiu+BEDoqpKVwFm44v7aK1gKu\/5HWNaP2Q+tEElNS7KB0KCijp8EPqtRBQUu+idCgooKTDD6nXQkBJvYvSoaCAkg4\/pF4LASX1LkqHgpkGZbdXQUoSoUyDsturIEUNStDCyzQoq58dhn8fNWH\/66\/tyipIUYMStPAyDUrp83tw\/++nyl4\/1Q9HbbddN15Q+bFMgKILiUnlKEnJSTOZmQBFFxKTylGSkiOg+FjAhlDuKKzzSvVe8d+gkp7dRALeN0xOmh0clW42\/uIyY\/81A+6o+5tP20krKkFw8J\/9SnqaRAIdTH5worw0lg7tZkGYgpQKUPZ\/U9tRFZo76tN3v2onrV+3qk6g8EjEK0\/rYMpSXVmTBWEKhK5dz0FZ+fBlGPznV6Ufrwodtt3wCemqSXNn66JPGCgImOrrU948bCWH\/b\/OKbb9koC6Z6Dg5Jr3foR9r2x12AudTqvfb4shh3tB0UUMameSEOvyHy6T9MN7T3\/5fgfgXC8d7H619Lmzcd7ehWMLj0sESeWBGxkxbELoGL+8pNW3D4a2HrYNapOvmABBenH5uns8oumSYFP9dJEyaHvxc7ALWNwntCAI6p5FFBdQwqAKWoX8ycUvKpmMHQYPj4Reh0cFCs2Rj+8XgXTRTLeN4n2vT3h0\/\/NZH5S\/fR7dC2CmBjcFxXQ87kSdw7iRybi24we194s+Xpm6fMrPHt5t2ERX062bt6M\/LeONKNxGD5\/1wUc\/vx0fKLqyF3GC0m3EMHGIaxsvKLRFxAGuqY4EDeVcfv1iBSXoTcG4QDE1Ti\/b8f0\/S3ZAWIrvvW\/1h62MDtwwmujePc6SgeKAiiJL1uxg+4etjEHRlb3ImoGihgVX50u1H2Doy4o6QMzKlTgodI6SpIH4mQ06ii7dfa4btaF+3vMfPp7pnDA5bIz8AW+92nmWZNo\/6XaF4YISafM1hnFE0W09KHBjY0N9krroqwKS1\/p4+cXE\/\/8KwXuf68XbYD9vH2xL4yU1n6TlFAoFwI\/NZQRKUDJrI0zaZtcCRqDg9PjjsRTNya7DXTU3BsVVgPTLhwWcQfHWocUcBkuO4sXr0Oru6w7wvGZFObVaTd0eHR2FpaUlGBgYCK1769XDVB72s5WJfZ48eQIzMzNQqVRgbOzFn4k1lWkjj+TcunVLyaASrzbyXNB1AoWcT0rqzllQIb8keGRkBKampmB2dlbpPD8\/D41GQwHAr\/X19Y7\/w7FarRZMT09DvV7f0R\/7+o2ru++Vh+1sZeIYmMOdPHkS1tbW2kfjpnmdrbzbt29Ds9ls1\/9FuSdOnICjR48a2dQFEuxjDYpfHVq\/6FIsvvi2mCZF0OB9rGMb9BSlmwwCiec5k5OTsLi4qKJLf39\/eyXrwNTdp5UfZLwwmRjlrl27BsePH4czZ84oB+K4QYeU3cjz6ox27MampuBYg0IDh1WNHBoagomJiR1A0H3dAV6Q4twofv0RQJv75XI51E5hMmkMiioclDjmyHXmUevx48e+czeZY6gRXCKKHyi0YjFHQRCOHDkCe\/fuVauLtil+HwGyNSKPWrTKFxYWlDqYG4yPj6tIFSUoJjKjBMVGnjcn0tkkdaBwKmkVepV0DZPe8WxzIpetx1QmbQV+EcVme7WRxyPJwYMHleldtzqTaOKUo+i2HlrJeEKLyebly5eBh0N+f3Bw0DjxQgNiBOJ7sy5RjCKZJahMZVJC7AXFNJm1lYftL1y4AHNzcx3Jv408Uzh4u0hyFL714M\/8QI4\/HvvVp\/W258rxx0u6XyqVALccfMKoVqvqtsm4pgeGLjJRBy8otMr9dOxmjtyeNA4dR5jOMVFQXIRJn+xawDmiZHfKormLBf4D4PcYrcS6LM4AAAAASUVORK5CYII=","height":101,"width":168}}
%---
%[output:58a62d12]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANkAAACDCAYAAAAAq77nAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl4VdW1\/jMAgYAyk5AI4QGK1bYO7cdcEhFxxNeCIlgkEdDSKrQiiIiGyFRUrKitVsAkIg51aIVKVaAJaCtYJisqMjwSSJgJswQI5OXf966bfU\/OnXLPTe5Nzvk+Pm7O2dNZe\/9n7b322v+KKi8vL4d92RKwJRAyCUTZIAuZbO2CbQkoCYQcZCUlJRg1ahS+\/PJLN5H\/+Mc\/xsKFC9GyZUu8\/fbbePTRR6t0yezZszFgwACVv0ePHnjkkUdUmjlz5uDPf\/4z7r\/\/frd7a9ascZVp1r\/btm1DRkYG9uzZ4\/b4jTfeUOV7u5j3448\/xgMPPKCS6X9LuQ8++CCGDh2K06dP4\/nnn8eYMWPU+xnz2mOvfkmgxkCmg0SAJ\/cIshdeeAHZ2dno2rWr6gG599JLL2HBggXq3u9\/\/3v1\/+TJk7F06VIIUBs3bqzuSRr+7Q1kOhiYb9euXV7BaWyv8W9jXfwICOD5zPiRqF9DzH7bWgEZv\/Q6KJYsWVIFZBykw4cPB7XMqlWrFKgIQl6ijdq3b+92T8DjqVuNGofp9HoIWgEwn1GTDho0yO3eLbfcoor\/8MMP1f+33XYbRo8ejbFjx4L18xKtnJiYiMsuuwz5+fmutPxQePoI2MOxbkqgVkAmg50DlFNAb5qMwDp8+LALcIWFhQqQU6ZMwaxZs9TA7tixo+u5t2mfL5ARzKKBtm7d6irz0ksvddNGRk1mLNfWZHUTLNV9qxoDWaBrMtFSnD7qg7igoECtqSZNmqTWRwRVSkpKFU0Y6HTx6aefxsSJE11rPwHSXXfdVWVdaIOsusOtfuarMZDpazKZLq5fv15N9zZs2OACiUwH27Rp41onyaDu0KGDWj9JWaIxeJ\/X2LFzsG9fHAoKgMJCoGNHICUFSE11dK6Z4cM45TQaRWhcoQFDX1fZIKufYKnuW9cKyNhYfYqog4yaS6yNRushLYq8xBoo6crKknHZZbOwYkUfUzkQaOnpwPDhDuui2drNmzHDl+HDni5Wd\/jVj3y1AjLRZGLVW758udt0Twb1wYMHXRZHAZRxGjl8+KPYtOk5EGi8RHPxf145OVCajVdS0jm0aPEQpk7to0ztxktfS9lrsvoBgJp4yxoDmXFNpoPFzPAhVj8aR2iRY35aG+VvWui2bCnF5ZfHuQD0+usNXFNDXXgEWUYG4DTy4Ykn8pCVlVZFvgJ+WjJ50booYJS9Odk2mD9\/vtqr499Tp07F+PHjXRpS2i5alwYVSSt7gzXRuXYd4SGBkIPM39ckEOSf2XrKWA7TpqU58nAq6LTue62O6Qk0armdO\/1tmZ3OlkBwErAUZKfOluHe3K+x\/NsS1aqMnomYe0c30xbqoFq1yjGtM7tkPZWZWfmUebOyKvMQYLqBQ1KaAZf5eH\/aNEAv01ses7KDE7uduz5JwFKQvbxqN7YdOKWAtf\/EGQz64yZk9GyPX\/W7xCVTDvDcXGAa92dp9euniZtrJ+daSt1dxUUVgAIHiAgKajnTvBVATckH0lOBkSO91KGVqbtGe22XVrYZMOvTgLHfNXAJWAoyY\/UT3tmCrm3jXSBTU7wMoIBayWlW96vJBJ8soej04S0v02YA8FUHgVvRFoKSpv6sXD\/a5cwz\/\/Fy9O57wa+mN4iORmxMlF9pvSUqKioC\/9lXzUkgOTkZ\/BfsFTKQUZONzN6Mp4dchh+2b6qmaJ04+POcTSYY+M8JmKTz5Sg2DkbRdsa3ZL4KbYiOQNKIcjx14gzWNojB8\/ENApcH68iqbBfbMbi0DNeWnsX3p06hSXw81sc1dC87Dch94BR+eNUJNGnSBLGxsVXqLSsrw9m929Dq2w\/Q9oZ7EdtGV9GBNZPg4kb52rVrA8topw5KAt27dwedFIIFWkhAJmuzMX2Scf3lrRTAlAYTgHEKyIHtND5MAzDp3Hn8oEGMwp1+Jd1TjuLrooB0590cIMW5rsI0ICUT+ObceVXHD+6LqQQxk+cDj\/UpQwMDCGj8yGdFUqazaGkH\/yw9XYp9+\/cjMSEBjeIaqRS\/bhCjZq+82K7lj+90e663W+X\/9wdotPBuJE7LQ+MrAlHd7jIQayU7PCkpKaiBY2f2TwL8oM2bN8+1J+tfLvNUloPsnQ37cP\/r36AcUZh9exc1VaSxYRrXXhxnHNydnJrDOe7kG68AVgAkXaJpNYdvbeUUkX87p44Jw84hZlEsusZEIZ+glTrkXfOB3v8kGBxTu5udIC6gRiVaCHpnG3rvKMfyDpVTQKbdXlaGjlExWH6h8v6A6Gj8ixo3H\/jN+0cwfUZD\/OKiePVa1IKS9kzpGRz8\/APEzB9uGcj8OZITzGCw81ZKQHcc93UMypfcLAXZ5r0nMOjFjXg1\/Uq0im+Iie9+hxn9f4ifPtCwUsMYpmemDTQCzZiIIGM5KUD3redRtCsaxaOj3LWY5EkDklY4Dn+7pqOSn6pLrJZaOj0tgZN84QKKoqNdrdDL+XrREdyQ2FyVzbSrSk6rdGfOnMHxdR+i7eIRVUB2cs17ODB3iErX4s4stLjjCa\/95K3DCW5+X3jxY2X8zvgaAPZzcwmELciGvrIJy7cccWt1\/4s6492oDpWDWUaEZpIPuKOphaiNeMkUlFZDszKZlqNPZmuiSZlXB5kxnT+NygGSOjins870BJpc58+XIfZYMbq0SXE183xJMfZMvxFtxzrOyB14aTTaP\/4RYlp6ngaadbiAy2zng68qYvHnNew0VSUQtiBb8e1hzP+sCK+OvEK1mntmp\/6dgg9HXlw5yM2mdZ56WTfp6791kHFEyVaA2bJHRiHXX2J5lCmoNl10LbYM6zSvA1Av20tCYlz2vqnFSnIfQtLMfyMqrhn2zrwZF9\/2OzTtMdhjCWYdrilztbTsWPHdqdjdUN8NXt6ARrc1eqyMGzfOPtvmQeqRDbJANAbTyqDXLY0Ei0w7xVKiA0YXnFgiRzotkhyFHPUsV59bSbpANKyfeYwgO7b0D0h8bJlqJUHW5KobvE4ZjR0uy1KWy9fW7Za6ovYkEhtkvnV3WIPssQ+2Y8lvrkLTRjFKk3U81BlP92yqPq0cDAUayHwZtZVjrzORyiuyoWuU83eB7J9pI8qYVm1qc2Q6XaqQ7XQa1qeRfGZicTTrDrd2+JlH12TBgoyTAX4r9Nmu3k5R9PyOOM6SQ\/GOyKlv+n\/Sd5SajKfS5SS3+GrqvpfiqykO0yxLTkcYjw7RMOPP6XKm40HYrKws8LgTjxfpPqm+IRD6FGELMn1vjGKg4ePmhj\/G+LhYpTmk0ztVeGYUpFb+7UlkMpj4nHmJCaWIcoCdTg1Hczx9EtVoct5LdU4H8zln0gqh10hensN7hG5calknmjKr4lm\/yrNn3rrRVWQAeaQ8K6aLYhzVXtmtuWYgoxM2Lzo8cwC9+eabijaB9x977DE1baQTdL9+\/bBu3ToMHDjQxbdCzUdAZGZmuoiPWFa3bt0UWJiXZdIR+he\/+IUb4RDTsQ4etuUpeCmLx41mzpypTrhzH4q\/R44c6aoz9DDyXkPYgozN5rrszvn\/VW\/wlzE\/Uvtk3IQucH5SaRbIqNjrynEObk+kj+Lkoe+biYZKzwGyDWsnvQ61viLAmCHNoaHoRGx0ierEYzDOcuj9sVM++17kr7fL3zx6cVYYPnyBzEzTCYBojpbpIn\/zfJ1+UZtdc801Lh4V\/ZSBfqCV2oyaUOdEEQ0nJxZYLsvjRYoI1keNSkD9\/Oc\/rzj\/t8K1LnzxxRfdgG2DLEAJsB9zuCZKrVyQd6LXPU8te7CEubSFYZ+MgNlpMs\/U61DzSqahL2OuZw\/9tAIgX6ajHso1vqreLm9lexORbsJvO+Fdr0YPluNpTcZnxnWX\/hHQn+majNO83NxcpdV0TWbWZj4\/deqUOmokmkzS6cAQTSa0fUwjWotA5dlAoyZ77733bJAFiCWvyTvxWIrTrkxgZdILRBvgagGfArz497\/j8R07cHT8eMciTKyRTo0Tt2YNcgoLXee89PXDyZ98hYN\/ia9sRw7Q\/Ng8XN+\/P6688kp1f\/PmzVixciXKkpNxcvBgRx3OdrDsxOHD3c6RMQ+\/vr9+6inEjBqFhfRlKwBaPHQMzTddrcrUz8ZZKTMpKxDrotOf2jU1lzL0NRkH\/RVXXKGmifqajGmFHUxOoZutyXQNJes5WefdfPPNirlLtJ4Z45esyXQLp63JLBg5au2k+QhyXDuVjaN0bkBvP4djVx3DydatHfecG84Oi0klGJoeOoRfx8ejz7lz+NNLL+Huu+\/GF40bI7u8vDKvP22WuajmU5l07hzafPghnrzuOjS76CKcOH4ck5Ytw\/Y+fRQw2Y6OmeW49uSLmDjxpz5JUf1phq80ntYHYsY35rf3yXxJ1PfzsF6TeWu+AhoHtC9PekMhCaWlaH\/3fGzY8KDvvLIX5ocXfuLwM+hw4Z+4\/PLLkcN5o692Ob3w\/\/qHo5g373dq0S5krL67rfopfHl80FdavhfOWXn1K7NzKglELMik\/5Qvo+yBGc6TpdDU75y+nfjqK9xy6BBuatxYWcNIQzBzZixmLm5Quc8lhRYAF71fghsa7cGECSfx298ewdr9N5mmS1kFpKcAcXFzFJ0cF+Usv1mzZzyW3eSdAxiWeLaCzTjZ1QFStU74E4oxamWHh6J9dbFMK2Vuqe9ioMKmZuOpaF5G+jbd7CsmZ52mu3nz\/0VxcVds3vwVRo3qjw4dyrFs2ST06tXLBRqmX706GpmZ\/8Qll1yCW2650nWC2lv5v\/zlAkRHp6l1HNvVqtVxvPrqSLVwNzqLylqH9ZqR8wQqE7P0Vna4Fe2pD2VYKfNaBZmnzqIJmBpGBi1fmPdIQiO898OGDVPZA7kvAPG3fEkv7FlmIJM28H\/dsmblQDTrcC5Xnd8n06qscBT25BnCD5SY5K18z3Aqq06DTN\/PEaGzsydMmKDWQLxIzz137lz1O5D7jLASSPlMz8sIMt1c7QuAVgwcsw6nDckDLYpan3nyBrGiPTbIApNiWGky3Rwvr2GkhON9\/VyVnsfX\/eqUbwYy0V5i5q6NNZmAzIx0i8f1uOOhu2Hy3Rkh59NPP1XyY0wBozuVvokszMliZhdKPm5XpKam4tZbb1VliEaT8FBkWpYN6lBvbQQ21ANLXac1WWCiqB+pvWmyQEBGtydOaWUzWnenuvrqq7Fx40Y3z3yZLv7yl7\/EU089pTaj9em6GchIr86LU32jO1Yk9ZYNskjqLQvaahXICArxXSRRrH6J+5NoN12TcYNZ9wyR6aIZyOgdIhqe5UeqNrNBZsHAjaQirAaZUZMZZWH0L9Q1mb6uJchEa8k6lQYruRdJMja2tc6BTMhKuWDnb3UCxcmfKHtmtdFh4dIuq0FGWRpDCL\/88stqzSZ06p7WZMzbt29f5TrVqlUrN0fi\/v37q+g3utOwHra4NvqwunXWGZApjvosh5Nu6kjHJrT4\/ubQgTjXcayFG9TVAZ8\/IJE07AwJTEGAk83KV7vIXCxhmfzpTH\/aY1aOJ5DJ+VNjnlBbF\/1510hPUydAJgDjERijN5MLfBrDlWIZqPAG0cHnaZD7C16WR2CTcVg0pt4uHgjNy6wKJOW5RRBmmD+vMuj9\/Jh4eh9f+2RyUFw\/\/WPFPlmkAyWY9kc8yBQPY4UGS3UCTBeGGfiUI6w24L0Ncm\/glXqopdK4m5vrDhJju1RZGY5zaGYay9gus071pz2+QOurw+WQpqdDnMEMtvqa15fMA5FLreyTCVW38VyYN\/AplgECTaNhqgI+L+AVociUjVrDmN9Tu3jy2hvQjO3S6\/L0MTHrJE+g9eUgzH0xuTwdgpXngdLRBTKY6lLaiAYZ\/RUzVgHZ1A6GXvE0yJlMkcfkAHlaeFon7lzg85ZfqhLNJNNDAS9B5Kldqh4DwF1a0UO7VB4n77\/ZIVNPA9LsY+Ktw+WUtJSn83oY66jOqey6BJxA3iWiQaZCHvUDdhoQ5g18Ihx+sVMz3E86C\/imFTrcjMzA6wKE0yFZpyGQ\/KnOgBMFdKUIhFGbBDx5jnZxfSnH3jK9fEy8dbbZx8RThwszlZBvCR2cnIguO+g4ACM8\/IHyi4gpn7wfzZo1C2SMVkl79uxZvPLKK\/juu++U9bGm6MaLi4vxxRdfKLqDQC7mY\/wBK1iba3S6uKukVH3de77fCPMvXMD5svM4efIk4pvGY87shshJBbZ5iZYyJjoaO0aU4ZOFF9yCPCgO\/SggaWW5W369fAaFmDk9WkVjMa6vVP5OwN3by7F4RhT2\/rpEtckYSGLMqGjMX+gezaXJrBgMmXgae38Vh+LXolx03en8mJi8j7FNZh3P9vS+p1zVxfRbvtuCvqXfI6FHDzcqcH4PhCmPLlT8e0XZefzs5ElkNY1H6QezcG7LarT9TbYCGkEWCFPW9u3bkZ6e7jrlHMggjfS0JPdhoI+wBJlOpDPLyYUvAldrjgwgaYGD+ppXYfkFdIyKVs52axfHqEHk6SIVdkLvc1j59kk0auQIAsFrUrOGeP+mWPziH2V46sRZ3N28kaLVTjh7DmTxPdi4sapPBa94LUo9U\/WnAkWvRaO4c5Saco77+zk8\/0wD5A\/aCroH6XWoesY1dMu\/+OgZdPlTE6TfdxjLh7ZE8soLKrqMupyHU7ufOw+mk4v03eS7kPKlrWyPpOP7FI+KwuLnzii6b6a\/64pO2NewgRsVuJmc3Mrf+A6iXrtf0YR\/eSIOeS9MxhWH16Ldo3\/H7t270fD1sVh\/ojF+MO4l5Qli9F3s2bOnAll9C3QhwSYo37ADmRklXG7GlWjXzAEITvcUT6IZh7ScmPb1CTRw1jO54qY347w3lEWQNXBqG\/VIozRg48irv3ZWDBIePaeAHxsd7RZlpvh6d759Fe5pRpTyxO39s3L1XioYhQYy\/jR+OAh8lr\/s7AVcG9\/AUUcB0P3u84hd7Qxoob1Pm9OncSiuIYqc0Wn08ghOaaf83nG+DDExjnBOyVtWKXETZNv\/9kf87NhqnB29GIvf\/QD34l9ocfsEZK\/YgIE7FuCLdjdgcOYraFR+TpGufrF+PSZsTbRkoPnq1nB67snpvLpttHS6SC1mJDeV8EkOkKWhIGMkkMkNKqJNOHDTgLRMh1VDqQCdCE57NS56cguATB64NyCVY7uckyed4sxFgeqoKy0PyGO9Orr4OxOISq8gc2QYUIb05OdA2qbVT7t\/npTvfM452kig+\/zzKHoyGsUx1znKT8t2NpHp9bIq6086\/09nEAznPdJ3cYHHy9PHyPVQ3sE36z2XmLPXrMHeLRvxky\/\/hAN9f4snnngCD3U8hCk7EjCm\/WFc1+oUtvX4Hcb+8W8YkXgUqT+4BDEjXlCazIqveXUHaG3kC3uQGbnw0y5t6Yq02QmdUEA3ilV5QKbO5s4NqWwn+HSDtC5iMpKSmTSjYtNKGHaczwm+ThV+WARJiqf8xB\/rqFgsKapgHbx5DmxysaZ2cSUypsY9RwDkplbgUcwLzM82OfzgkzpkosF1q1GgVklSFz8m0h4d8FK3h3uuj4lE5zAONSMbpfehmIIUvLlmtjqackvHWFdEmeYPvoGGDRvixMqFuPD9ccXJ3+RHA1A8\/UasPdIA5\/r9SjFamYGMLeCnTnqCfysWMi9NiRR68MgHmfrKU6PI59rZK6T7dQOfAWAceIpVSni5tYFK8BWyyzM18Jr0tqrDDChO8GZRm3KomISfJRsrwSn0xKr4FKBgJ5CRj9hFA9AuMQnFMQ6nWXh9H2+gMH5Mgv+WO0D2pgKZnDan7+KzmZPwWMoB\/GFXa7yYdjHe3huP9\/+7H7O67Mfxq4agw81jTDWZMPUZSXt8RcWyQRZ8Xyr2YG\/TxTSk4bOizxSradnj2yqmhzmILRrjqrlsxHIFvtiiAZX3FA1bOmJnjgL6iRZzPHZQtGUjdkCsAp\/kJxpjPcRXZhoCwiy\/atfKbUCmoV1lyYgdsMgN4I6604HcTMR2HqDAF1sUi7LkMq\/v403MxvcJtksSGpW5zPc7XXFlKkvd9+xdaNrrDqW9JLqM\/nszksAjMUZN5i0wj9ORJqI5+MNak\/kyfLB7JcD49u3JGLEwGbGLijCjqEjN0srKktH38WRcv6gI04uKVCig9EIg4XQPTI8tQpculYHJ+SynENiypQemd3Y8k\/wEKidib65ZU2WcMg03pGcsKqqSn4l7ftxDzYGkXSmFUAxV9G+U+vV2FRQCnw+sWo8CseF9vIHG7H0kvVi7ArXyEWQ0Rcs+GUH1\/ecOTvzW972Co+\/PQNmhXW7NavubRTi2\/GUVaea7S24wBZkvGxUnyHKYNBI5+MMaZOwtMy58T4PLzZE30xn1hfGcORGj5z098OlfmF3pwMuyxNevIB\/INniACIkq7QeeaBRz8p3RKSuISo35o7KAnSMrTwek0HmY6w3ODJ2e+vk8HUAbCJeCuUC5l4WIP+3x9j58Xyu9D8z6gusxPU7akXeexJnd36Cg12+rBTIdhJHIwR\/2IKvOFIdgy811Bkt3okgNaHLnO8EnJMJMqw9ys\/NmZuA1y0\/gML\/uASIgY3sY+YWBKlQQdy0OGLWabGgzvTeQqdfRvfCr8T41DTIBnScTPqeL3gwcOsgikYM\/4kCm8xtywBkjfshC3Ow+B+eMGUVY\/C\/HNLF7u3bo1KmTmrrt3\/+2iwhGPxioC+ihh57H4cO3KpAcP34cWz9Zrsp56KGWyMpyGFD++teNuPfe1Wpqd+2112BVWiekFwAdOuRh0aJRKo2n8nlk\/64tQ10gEw5G0tUJnZyn9rDc8p07sX79BjRr9i66dy9VlHc8eaznYR0kq+HaiJvDn3\/+uWrTTfffjzseeUT93rJmDZ5y0gmkz56Nnw0d6vatC8RLTM\/oCdze1mRiFDGLixYpHPwRBTIBjrA5sfGB8CQy9pUZ5dvhw4ddtHByHJ4BDiZNmoSpU6cGRR2XQytmIKMyHyjPq6SN48liMRR4orJjm0kBINR2vt6BABszZgxuv\/12HJ83D6+VlKDdqFE4ftddOD1gAFpPmIASJ11ey1mzcGjuXFxw0tl5cxj2NetYvXq1x30yaiujdVEY0n3v3PmquXafRwzIhGyFDEm8yJJkptWEE0KCxIm2E\/psM1AaWXwpFNJsMwDdc889ZxkJqnS1r\/L5hX7++ecxZMgQPPzwwy6mYU8fFWP7Zd3l6R1Y\/owZM1zg5QZDizlzcFFKChI6dsTROXPQduFCRDVujEOTJ+PUsGEo7NHDtQVuEmnKr1HsDWQsQGjp5JvEegKJBuxXI2ohUcSATGRDkAjIzAYd2Yzo5W0Ek9wXDnyWQe9tMzpsYQQWTnudzttI262X4ym9kW7bV\/mS3owE1Z\/2Gz8sxjwJCQkqkLpoyMklJcieMAGdp0zBs4cPu+IEsJwHJ0\/Gxl69sKFiyqgTnFL2gfIuiu+i7fHRIyioW+pWZdYSHWQymEgZJiSZPEJBLSdTS\/0+wedrkPrizA8WZP6UHyzIfNWhg0xiMn\/g1FaZa9bgiDMYB+V75+TJyOvVC22GDnXbGSPIAuVd9KXJghp5YZw5ojWZUa5GTnpd+\/kzXfSX017A7S+XvlWc+f5MF\/15B5ku8uNErUK68uKuXZV7Y6M1a3CTM07AgcaN0XfyZJyoMLx81KOH29KSbQmUd9HWZGHohe9Lk8nahtM5bkrTUEFmWhoyRGPp90k5Zmb4CJTTnu2qDc58b4YPAT61tb5GM8tDw8d9992npsovvPCCskDyIshWl5TgqgkT8NaUKWCkgNWzZuGyuXOx2plGX1cKyPzlXZRonDU9XaztQBcRr8l0U73eeZ7um72wfk8GUbhy5nvqsEDege\/IbYGoKHefSprrM4cOVdqsvdOEf+GNN5Bj0GJiWDH6LgpbMJ+b8S5602QEuJkRlo4EobIu1lSgi4gDWRhPvSOmad42owPh+Aj0hb3V62lD2ni\/OgYXblfUZqALG2SBjpQ6kN6Xx4d+MMcsAEV1RWAVyAI1uNR2oAsbZNUdMRGczxfIQsW7aBXIAjW4iCarrUAXNsgiGCzVbbp0+vjx49G9e3fTYnx5xlenbm+MTYFMFwM1uJBZasWKFaitQBc2yKozWiI8Dy2upCfjkZfauMysi9UBGdseCYEuBGT8oPF4EY8LBXOFfDM6mMbZeSslIOfwalImco7NDGS1YV2sqXf3NT0PtB02yAKVWD1Kb\/VgixTRWf3eNsgipedroZ1WD7ZaeIVqVWn1e9sgq1Y31I9M\/hhc6qIkrKTopnxskNXFUWLRO9W2wcWi16hWMVYZPUICMt38yVO9xmMjZWVl4D+rL\/LWG7nrq1NHbRgYqtPOmsrDrzplUt8uWhQZFIP\/h5V10ZdDLDtr8eLF+MlPfoIOHTpg3bqmGD4sEVu3bVN9+P3336NJkyau3+SiJ7e7fs\/TbwqkdevWQY2F+vzlDkpwdTizFRrN0umir6MdouWmT5+OvQMHIrUA+PN3Z\/BVq1ZYGxeHGbFFGNXFweG+f\/9+XHzxxTh27BjatWun7n1UWoo\/nmyNkflA06aHcOWVTREbG4OoqF2KB6N58+ZBdbe0L1DqtaAqtTOHrQS8bWEE0mjLQebtkOVnRUW4ae1aDPsuEfPb90J2agE6796NM2cSUVQUi\/grDuKhFlfjqg+O4Z6f\/R8ubdhQEeDweWyXWBXd5FSbNvjgb81RkBGFQwfj0frDU0i8q1QB8s24uEDevUpaq61KQTXGzlzrErBqPNQoyNjoO++chMETJqB3797ok1ymNJZoKvl9R+vWKN1Sihdan0Tj\/fuRs7Mj3uoWhxkNinFLfLwS\/pfHjuHHF1+suBC3lJaiR0ICutkgq\/WBWZcaELYg80Z8I43mOaa+ffuq\/pBYXfrvdYeaomRDCZbfnYR1B+MxO26fYhbtHLNbxfUyy9e1a1fwmH4wl1VCDaYNtZE3FIckeTClG18AAAAM5klEQVT0448\/xgMPPOD2ShwflLPQ3\/EhXa3eeustdY8HeDMyMlyBB4XpzJdchF+F6aTsYLn3rRoPlmoyX4YPaXT\/\/v3Rp08fZb1ZuXIlunXrZvo7rls3rNm3D3dVEMn4SnvjjTcqQp5gLquEGkwbwilvMIckPYHsxRdfVH1OXheeCJeQuSdOnFAn5Un2M3LkSPCjyUtnIPYmG7P66iTIKARvHsyhtN7RQ53\/grnqAsj4DoGyUoXikKQ3kMXHxyuLMbd3mO79999X3SbtMFIyGPtUH2M8EU9wUmOS81LXfDrIlixZ4iLDlTSi\/ZiPFA+cCema16rxYKkm82eAh2ofyor9DDOhksWY\/4SWm9z2vMLlb7ZNpyrnO4TDIUlvIOMWDts4btw4cPBzCfCf\/\/xHgUz4W0gaxMs4XdRnS9R2QgN+zTXXVJmeCsh4Lu3ZZ5\/F3LlzXeUTyGyDANrICaorjGA5TmocZP4AsbbSmIEsq+JMx7RpFXz3zlDWQrMRLn+zfTqXP98hHA5JegPZwIEDFSC4Lt+wYYNiRibdgIBM73+dS5\/3jVNA+aiQ2Na4BtRBRkCJlmKZZ8+eVes+qdOsvRGryWoLQP7UW1c0WTgckvQFMoLrq6++UgxcaWlpioF58ODBio0rMzPTxcilfzQEZMI8Zmsy56jWD+l5C9wg7lee3LJ8uWv5AyJfaaz6cvmqJ5TPjYOytg5JEmS6lZDvzGnXunXrQE3GSygBOdUnyKhVjNZFYSFr3LixS2zGNZlQDHrSZDItFXauOrUmCySoAufLvALlR7RywNYFkFkpj7pelm65NE5LI3ZNxkEcysAQwQ4KG2TBSjCy8uvaVp9lyVtYNR5q1PDhK3BDoAEjjB7+wXaxVUINth12\/vCQgFXjocZA5iuoQrCBIazoFhGqN1YoK+qxy4gMCVh1eLNGQOZPUAXSUPMKJEigWZyvYLovlJvlwbTLzlt7Egi7oy5mojBzi\/HkflXbhg\/WH6rN8tobJnbNwUjACieHkGoy3dQqL2plYIhghGfntSVQUxIIKchq6iXsemwJhLMEbJCFc+\/YbasTErBBVie60X6JcJaADbJw7h27bXVCAjbI6kQ32i8RzhKwQRbOvWO3rU5IICQgO3W2DPfmfo0xfZJx\/eWt6oSgauIlTq55DwfmDlFVtbgzCy3ueKJKtWZpLnx\/HHtn3owzW\/+l0se27oCkmf9GTMukmmi2XYcPCVgOsv0nzmDQHzdh24Hv8ZcxP7JB5ucQPF9SjD3Tb0TbsQtUjgMvjUb7xz9yA4qnNEwveRtdah4k0M9m2MlCIAFLQUYN9sQH23Fv7yTc+9o3mHl7FxtkfnYaNVRJ7kNKA0XFNVOa6eLbfoemPQa7SvCUpkHLZFNQ+lm1nSzEErAUZNJW0WY2yDz3XtlBMkYCsW1S1P8E0LGlf0DiY8vU3wRZk6tucJsymqU52rEX9h45iRN52a7KmqVloGm\/e0I8dOpH8WHrVuUJZLtKSvHmf\/Zi2E8T0aGlO9tvTg5QWGhdx3XsCKSnW1ee1SUd+UsWTn+Tj7a\/yVZAqw7INk4ZiN9vjsL6\/9tndfPs8pwSCJmDsL6uGnB5S7w68grEN3Rw1Mu14tvDuHP+f9WfepoJ72zBws\/3IBpRiG8Yg2XjrsEP2zdV6T7bfhSD\/rQRS359Nfp0ceetJwtUWpo781J1eprsTbyysx0gEzYnuU9mJ\/ltLF9Ynzw9r057POUp\/TofDXPSkDgtD42vSFUgC3S6uGx8f4xbeVjFNSYvpX1ZK4GQcuETKF3bxmNEz0RTK6FRU728aje2HTiFJ2\/votIPubYdnvmksMqazBvIKB6CjFdeXvWFZUUZ1a\/d\/5ynv87H3mmVIKuO4eOT6RkY\/69SxZ1h9bEf\/9+k7qYM2aFNTwCae0c3lzS\/2nMSE9\/9DrkZV6Jds0agVpv\/WRHm\/LwrfvXGt5h8Uyc88t62gEEm2ky0UKDdxyknKdKYX3gRAy2jptIbQcZ6D702EceXPqOaEN9rKNr97i0QfPvmDkXChLeVpVE34RfeNAejZrxsgyxEnRZSkI3M3oynh1ympnnUUnlbS9ymjGZAzP58Dybd0AljXv\/a9coZPROhg9OXJmPGjAyAYBFeQ3\/lxyketRjBRZCF+1UdTWZ8J6sGQbjLigy\/ZBxu1aqV4mck+5TOXhWq9lsl3yrWRQLIF8j4UvqabPBVbbH7aKlLs\/G5bEinXdoSv+p3iZKDgOyxa8+gd2fPscR69uyB0aOLMGqU\/xEeP\/64B2bMKMKnn\/qfJ1Sd40+5rY9tR\/SCEQGtyeoryMifTxo54cf3R75WpAkpyLiZLOZ3WW\/pGsn4AjJdNBpIuLbjJXlpXezUCYj+\/hDi4tYiNrYIJ08OVr9LSys3UcvKklU+PueVmDgcZWVJOHHC4Q3Rps1EtyYw\/e7dq9X9pk3fs0K+IS\/jnsQjuG\/AVQFZF+siyIzc\/WT6FYpuHvAlMzD5G8kmRQPPe++956L3NvIoWt1pIQMZG+qP4UO03f+0jlPGDmqsLm2b4B+bDypQeTLj+3u8\/9FHeygqgEWLfGumESOSQaC9+eYaq+UcsvISGpWpWMSB7JP5AzLh7g9Zwy0omFZcseTq3P2kpThy5Ag6d+7sivbCCC8kLaUmk+miJ257q40\/IQWZbsLX11UE301XtlFeHPp00Zgm+\/O9qiuMa7JA+keMILQ0ejNi+JsukLprI60\/Jnx\/QCbc\/bXxDv7WydgCwt9vZDzWKSvat2+P7OxsU5AZue1Zd7hSBIbE48NfYftK58scz682DSW8gjH7+2pHTTz3x4TvD8giUZPp3P1CtS1xy+qsJquJQeVPHRwwXMN5MulHksnen\/fVzfNtJ7zr5rdolt+q6Yw\/bQtVGl2T6fHCqMVSU1Nx6623qhBHS5cuxbx587BixYq6sSYLlUADKVcCJRw5Mh5Hj47H4cMlKtKHDCwxdtCrgyAMJlBFoEEx+B780k6ePBnki5S1AOnvZNEuU51QWsTqAsgCGRM1ndYq+YbldNEYqKJly2NISPgOy5cnYerUqZgyZQomT26HTz45i8LCaCX76gaqCDQoBoGuf3HF28IMdKEeFFYNglC3M1LLt0q+YQkyY6c8\/vgOzJjRGc899yXy86chMzMbV1\/dHAMGLMaTT3ZWya1iHvYVFIOmZIb5GTJkCB5++GFX7GNjBMiaGFhWDYKaaGsk1mGVfCMCZATQa69lIC4uDn36PI6NG\/+AXbuicPvtv1VB5KwMVOErKIZYsESbSYBxI5GrMQxrKAaZVYMgFG2rC2VaJd+wB5kEquje\/RHlNtW9+z+wdu1N+OijM1i2bJKlIPMnKIYnkOmDSqaO\/ABYbVbW67FqELBMGpECvcL5KFGg7xJKw1JYg8wYqOLaa49jw4aLMGzYGSxceMFleLBiuuhvUAwxchg1mbGTWB4varpQXVaCjB8wCTrvT3u5dxmqbROJ9Tx69GgsWLDANJa0tzYy\/6ZNm9CzZ09XBE+upQO9rJJv2ILMLFDFpk1H0bdvGZYuPYmkpHOYNWuWimgfjOFDANqvXz+34yKegmJIZ5lNF2lqJqh8ATDQzvaU3qpBwPJ9bZdIG2pi2yRYkBk3uKsrb6vkG5Yg8xaoIjt7J5555lYlN\/0clZ7Hn\/sieG910ft7+PDhVeriDTMg6Sb8SFyTibeIpxMQAkTZNqnu4NXz6ZZa3p89ezYGDBigvO3NNBnTZ2VlYf369aBvY0pKCsSHUfwbJ06cCPYdvUXYv4xFvXXrVldfStATuu1JXGuzLZc6DTIrOq8+lGHVINBlxc1\/T8eF6F3DKeXOndZJd8eOHWjRooXaA+V2Sm5uLsaOHYvXX3\/dK8gyMzPVcReCh+DSZyT8TS+SQYMGqeni4MGDMXPmTLX1w31LiQ\/NdDSaefJ5tEq+YanJrOvCul2SVYNAlxKngwST0WfU3+lkdSSuzwCoZSZNmuQVZPqZMj3us8xuzEAm3vsEpjglU8ONGjVKaT1bk1Wn5+pBnlCAjGIz8xnlPQLNSi3GuqhVRJv4q8kEZLTi6oCTdby\/mky3\/ArwdEOVVfK1NVkEg9GqQWAUgfFkQyhPOuhrYhqfSAiUnp6uzo15WpMJsNhuurbRr5EX18FcoxG0XEtz7bZnzx6Pa7IlS5a41nO2JotgIISy6aECmWgzai5OGzl95PmvSKB1sFLeVsnX1mRW9koNl2XVIDBrtqzBaAQRsMlByxp+zVqrzir52iCrtS4MvmKrBoGnlohJ30qTffBvXXMlWCVfG2Q112eW1ySDYPz48RXuZqEJNMGpYn2bJkpHFRcXg3tuwfJa2iCzfOjXXIHcTOUgINOtfYVGAiGj6Q5Nc+1SQyEBf4mJQlF3fSgzbANO1Afh2+9oS8BfCfw\/KxJVlMTbRI8AAAAASUVORK5CYII=","height":101,"width":168}}
%---
%[output:7283f847]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:38cde80a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:728007cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:05c31929]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:37eec7a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:0276ef3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:487a130a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:2b0e1e85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:01715020]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:748b463c]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 422\nStatistique T_max     : 39.2389\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Mar-1994  (indice 73)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 72\nStatistique T_max     : 15.0788\np-valeur (bootstrap)  : 0.0090\nPoint de rupture      : 01-Oct-1991  (indice 44)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 349\nStatistique T_max     : 15.9354\np-valeur (bootstrap)  : 0.0090\nPoint de rupture      : 01-Apr-2015  (indice 222)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 43\nStatistique T_max     : 5.8363\np-valeur (bootstrap)  : 0.3240\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 28\nStatistique T_max     : 4.2832\np-valeur (bootstrap)  : 0.4820\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 221\nStatistique T_max     : 10.4318\np-valeur (bootstrap)  : 0.0740\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 127\nStatistique T_max     : 3.0458\np-valeur (bootstrap)  : 0.8680\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 422\nStatistique T_max     : 11.1204\np-valeur (bootstrap)  : 0.0860\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1\nTaille de la serie    : 416\nStatistique T_max     : 5.0187\np-valeur (bootstrap)  : 0.6410\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1_inv\nTaille de la serie    : 416\nStatistique T_max     : 5.0187\np-valeur (bootstrap)  : 0.6380\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 422\nStatistique T_max     : 51.8020\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Mar-1994  (indice 73)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 72\nStatistique T_max     : 14.7942\np-valeur (bootstrap)  : 0.0060\nPoint de rupture      : 01-Oct-1991  (indice 44)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 349\nStatistique T_max     : 16.9675\np-valeur (bootstrap)  : 0.0060\nPoint de rupture      : 01-Apr-2015  (indice 222)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 43\nStatistique T_max     : 5.9481\np-valeur (bootstrap)  : 0.3070\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 28\nStatistique T_max     : 4.1833\np-valeur (bootstrap)  : 0.5240\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 221\nStatistique T_max     : 9.6504\np-valeur (bootstrap)  : 0.1120\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 127\nStatistique T_max     : 2.9865\np-valeur (bootstrap)  : 0.8700\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 422\nStatistique T_max     : 12.6938\np-valeur (bootstrap)  : 0.0370\nPoint de rupture      : 01-Mar-1988  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 421\nStatistique T_max     : 5.1900\np-valeur (bootstrap)  : 0.5970\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2\nTaille de la serie    : 416\nStatistique T_max     : 5.8642\np-valeur (bootstrap)  : 0.4950\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo2_inv\nTaille de la serie    : 416\nStatistique T_max     : 5.8642\np-valeur (bootstrap)  : 0.4630\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1_Bac1_homo2\nTaille de la serie    : 416\nStatistique T_max     : 415.0000\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Jul-1997  (indice 100)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac1_homo1_Bac1_homo2\nTaille de la serie    : 99\nStatistique T_max     : 79.5488\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-May-1995  (indice 84)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n","truncated":false}}
%---
%[output:31dedadb]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANkAAACDCAYAAAAAq77nAAAAAXNSR0IArs4c6QAAHoRJREFUeF7tXQl4VdW1\/hMSCViRoWBCIgQM4CxQWxCLECdarfgcAcVHoiB1gtI8BkEaUEEtRcpgBcFC9AmOfE+ofSBgwOmBAo5VUSwBgoAKtTIlkOHx73vXZd+Tc4dzuCPZ5\/vywTlnT+ff+79r7bX3Wjultra2FuYyCBgEooZAiiFZ1LA1BRsEFAIRI9nevXtxxx134KOPPvKD9oILLsDTTz+N5s2b44UXXsD9999fB\/pHHnkEV1xxhcrfvXt3jB49WqV57LHHMGfOHAwdOtTv2dq1a31l2vXjV199hcLCQnzzzTd+rxcuXKjKD3Yx7\/Lly3HvvfeqZPq9lHvfffehX79+OHToEGbMmIEhQ4ao77PmNWPMIBAVkukkEeLJM5Js5syZmD9\/Pjp06KB6QJ49+eSTmDdvnnr26KOPqn\/HjBmDpUuXQojaqFEj9UzS8D4YyXQyMN+2bduCktPaXuu9tS7+CAjh+c76I2GGmEEg6iTjL71OiiVLltQhGQfpLbfcAkqZNWvWKFKRhLxEGrVu3drvmZAnUBdaJQ7T6fWQtEJgvqMk7du3r9+zq6++WhX\/2muvqX+vueYaDB48GHfddRdYPy+RyllZWejUqRNWr17tS8sfikA\/Ambo1S8EIq4u6pJMBjsHKFXAYJKMxNqzZ4+PcFu3blWEHDt2LCZPnqwGdtu2bX3vg6l9oUhGMosE+vLLL31lduzY0U8aWSWZtVwjyeoXWdx+bcRJ5nROJlKK6qM+iMvKytScatSoUWp+RFLl5ubWkYRO1cUpU6Zg5MiRvrmfEKl\/\/\/515oWGZG6HlcmnIxBxkumSTNTFDRs2KHVv48aNPpKIOtiyZUvfPEkGdZs2bdT8ScoSicHnoeZjfG9n+LCqnFajCI0rNGDo8ypDMkOWSCAQVZKxgbqKqJOMkkusjVbrIS2KvMQaqFslmXbYsGF+8yer1dBOXRSwghkzQhk+jLoYiSFX\/8qIKslEkolVb8WKFX7qngzq7777zmdxFELZqZGUPjRS0Hwe7ApGMubT51JmTlb\/Bn2svzjiJLPOyXSykEBDO3bET849F+np6epbKysrQZI1btwYzZo1w+HDh\/3u26Wk4O9eK6WonWL+DwRWKJIJ+WnJ5KUTV9bmZNlg7ty5aq2O9w888ACGDx+ujDAkulgsRerSoCJpZW0w1h1q6ks8BCJGsnA\/rR2AXAC9LRnKygDf31Ygt603US6wJdzCbdLp5W7dCrRtC+TmAr2tDTiOOkxWg0AwBGJGsgOHq3B7yT\/w8g2dcNL7O3H3j5WYetOZilglJcAELjFx4PeyNPcoIyeUAMXFzjoyaLlrgNzVQEHv8Mrdv\/YVfDv1RtWAZjdPRLOb\/lCnMaHS7Hq8v8qT+fvnnX2ISZ30CMSMZLPXbMdX3x7A4pvOxI2VVXhz2gb0bdcG82ZnoYwECiZZyoDcwvBJQYLlH00fbrnziwNLtuq9O\/DNQ79Cq7s8u1G+fXIwWo9fhgbNs32dHyqNELDxRf0MyZKeMs4\/IGYkk6ZRXSwAUP6XrzHvpTOAUgeN9pJt7vhaXNyzxjbjjm2p6DA4xVdudnUtBtbU+AnINTR+pDc4lj8fWOEtMz01FWkNUnzvSJC9Jb9H9qR3kZJxCnZOugqnXjMCP+l+Q8g0P+R0w7bPP8TeFycirWUuaioPoLmNFHSAgEkaAoGcnBzwL5GuuJCs21fVWNK\/Goc2nOTBQlRFqlOHj2DXSR6jiFKvDh9BtyM1ePXkhn6k6Da2Gs91rkRVVRUOHjyoDCe7dqbj1t9lYMebHpIMO3AEww4e8UuTlpam3jHf7xoCy5o3UffZ\/1mLNWO\/QKMNi9DqytsVKXiRZP9eOg1Z4\/6u7kmyxp2v9FMZ7dL80LYHHi7djnXr1iVSf5\/wbenWrRu44SCRiBZzkrWtrcW2B49KGs69qCJqBMupqsLrW7bhynZtUO4lA599cugwrjulsUrquwqBi+fVogxAVmUFVlTX4k9TGmHSpSmq3Bv278e82lQ0zGiIikMV2LV7N7IyM9U9L3l2Z7s2eJt1rQaKZpbi7tRLkTWhFI3O8eivbkn2eZPzceecZarDs7OPqZYn\/CiP4wfyB2369Om+9dU4NsWv6piS7KWNu9A\/tylqPs7wEIwM4Z82HyOphGAiyT6rOIyMRhk4O72BSq4uMo56n9cgkkvrZOFRw0mpR\/q9sXU7Mk87TeW7Kr0BNldVoW1KA6yo8aiZlRWV2LlrF3wk48N84OufptQhmRt18csz+uKOh2cnXIcnysCLRjv0TeChXJqiUX+gMmNGsk937kPfWR\/g0IRfYld6qqc9ZIxHKwt6cV7Fa4c2V1IP8j2ksiNd5\/0H0aBBmlI99XxSFvO0PHRIvfeppwuAhQvy8R9PFPskGY0a5ff3QPWebaqa1CatcPqfNoY0fGy\/6hEMHDrMlmT87IneZvPzRaiHwsG8D45AvSdZv6c+xIov\/oV\/P3ARappnHCPG8axXLdAW3Thqj3e0LgByvitD3shcH3dJsu1je6BGSNaiDU6f\/K5q\/66p\/ZBZ9IIinG7Cb1X0Mj5Fdh2PASEXm229CIMTG5AhXF0E6j3JVn6+B3PfLsfawedjW8ox691xDRYZrTRXRoJkZEEJkFt8bAE8HOui3TfYdTgFL7VcSi82mevtW49mnuAtIBjRuAWNu0+4b9P4qdmPGkOySJJM1EwvKdS8TJdqbplLBnCZoMCfZKGsi+GQjEWTZCQYJZauJfMzuLTBi+\/shLshWehONST7fA\/GvboZO2\/+Bba395dkYUzLFMJ+Rg+vZZIGD4qFMk2qhVuedJuv3IlArlfllK1c4VgXwyEZBS0lFv\/sNq+w+bTbUMJ5\/MKhYoiIBzcdX7kPlJKMHubilS37LvV9lLLvUjY\/syzxdLC6AdGDIRxPcaajU+vEiRPBPaTcrM02JZIHeL0n2e59lRg0\/1NMubETfjW6Mb551mv8sAysQL9XJAIlgSKEmP0nAqW9ju3WaMddHt4RGm6cO71c7irZIiPc25BIqYskEInE4kkk62VHMm6o5iWbkRctWqRCIPD5uHHjlNrIDc29evXC+vXr0adPH1\/sFEo+EqK4uNgXxIhlnXnmmYoszMtByU3N119\/vV\/wIKZjHXScpUe7lMWN0ZMmTVLe6lyH4v8HDRrkqzO0rIluinpPMsLLednNcz\/GgeVn4fATmX56UShSiCTwibNcz1YrnRSFhcCCQR59K1xDgq9c7mUsAbzhRXyjIdSWqUDDxtrhoUhmJ+mEQDRHi7rI\/zP2iX5RmnXt2tUXE0X3GNCdU+188UTCifcBy2V5vBjugfVRopJQ1113HVauXOmbF86aNcuP2NGlUOjSDcksGLXLB8o0c1ogUsi0SxkH9HW1AKTQyw1GNGu5uRPrSjFpstVyqG+pCpdkMidjeuu8S5em+jtdklHNKykpUVJNl2R29fP9gQMHVHg+kWSSTieGSDIJwcc0IrVIVLogWSXZK6+8YkgWmu9+KWK2TmZtFwM75fPnWyNaxtq1aDZ9Os7Iy8OokSOxsUkTzKmsxM6G3i1VoiaWAVljK5Gx9ixVrO4P9pe\/fIZ7XjrbV65Y8vTN\/VzDpnqmVE8aOiYC8wdF1v3FiXVR2qLPx9g0fU7GQX\/OOecoNVGfkzGdRPoSj3K7OZkuoWQ+J\/O8q666SkXhEqlnF71L5mS6hdNIsvDYFjeSsXmKaNR8OFEJd72sDGjzhxqc\/d0UzJhxvfpKRrOaOnWq+n9RUREuu+wh3DY+J3S53g3HpfM9PmaRvAKpLmLGt9YVrnobyTaeaGUZdTFIj048KtEmiFUghMtL1i2VeGLkZ3jmmQdVAB4JeDpgwABVA+cW8jw\/fzXW7f61x9KgizKSaw1QkBueP5mbwRisw0VVFaumdxrpphqTR0PAkCyM4eAjmw0p2pbWInPTMowY8aOakNPSpkca7tGjR8Dnp53WD2uoI6rJfGw8oxO1w8PohqRNkqiYx1VdtOtNTtrfeKMGrVsPwPbt5di2bQ2Kim7EJZfUqDWjYGQKRL5QgXeiMaoStcOj8a2JUmaiYp5QJKOqxwCmQgqCpqt\/JJmdWhjqeTx2ZNt1uDgOBBqUx7v1kuUG2hnCHy8xyScKKSLdDkOyEIjqa0KSlAOGhgwufvKyGjjCfc4TV2J92XU4jal2m4PZNs7PAu0GiUTbDckigaK7MhJCkulbguQzZMsO13p4IAUvPYipNRybSKtAz93B4z5XMJLZRd+SsAz6liuWwdNu3nrrLfXtPB\/Aup1KX0SWKMhiZhfsuB2rd+\/e+M1vfqPKEIkmRz0xarJs39JD+Ln\/+vjkNJIsPrjHrdZIkYzbnrggLIvR+naqLl264IMPPvDbmS\/q4sCBA\/HHP\/5RLUbrFlg7kjFUOi+q6dbtWHED0EXFhmQuQEvmLJEiGUlhDaQquMj2J5FuuiTjArO+M0TURTuScXeILGSz7GSVZoZkycwYF22PNMmskszaJOv+Ql2ScU4qc16STKSWbKuisUmeufjUhMliSJYwXRGbhkSaZGy19Tjg2bNnqzmbhEYPNCdj3p49e6qtUy1atPDbSHzZZZepk2z0QxH1I4hjg1ZkajEkiwyOSVNKIJLRgmi3gyva1sWkAe44GmpI5gUv1FqRFeNIrB0dR7+5zhpqnYykIha6b1myfqtrkCKc0ZDMCyjXipyEtXeaPsL95rq4UB0uTpqBnDhdV1yPM4bCPF7QxHydLBzSMJY99zEq97FBHlcUzs3VPTf3ctN+b2AQ30V493y4HeHXRps2hdogLDE9WF8oh1Wm+ddLD+Lgh6+rSMapjZv4RcfSD8Fgun+96PkZY9Qs8X0L9Dzc702GdIZkYUgyDtzCicBqHm3ErelelcrXwSTU0fe9ex2N6UHSMbLUAo83c6yOQrK2kSQXni9ggFVvm37727WYM+cW27iL4iUt32X1I7MO6Mov12HnpD5IzzlXkay2Yp\/tIRhV35f7DsQ49OVaX3jxI+Wf2z4nWU+ky5AsBMlk8DJGB1UoBsgh4eReAg4XlgGrC4HS4qPHH\/UG9PtoEc0XW8Tr4Km3yTpIKW2lTVmbBmDx4uHKhV8uiUwlzqQSDi5QlKqagz9i97T+SD\/9bFRuWqtIdvDjFbaHYBzZ\/g\/s37gMf8UvcdEvuqJy8YNoeu1IVH79Pio2v48WAx9D7eGD2PPMKPU8Pcfj9Bro4oGMTz31FDZt2qSsj7EKN75jxw689957KtyBk4v5Ro4cmXBRmxNCXSTBzhhegQGLG2L24SPY\/FUlrpt2iop1z8E4zhtaWwC\/IjUV71ye4juJxXrPdNVV1di\/fz9O\/snJkEMm9A4L9V4v4+fNmqC8PE1J0Vvn1WKu1p5A5Uib5t22BVdfcypuYhlpaThSXYucmhrlo0qljurzyqpqXLJ\/PyZW7kHVuwtxSv4g34EXVPN4pZ9+jk8ykWR2YeqY7p+fbMDwN\/b4vJydDNJkT8vgPuXl5fWbZAykc822H9FoeRkmX5uH3\/Y6XfWrijJVDGSf7g3HPTgF2Q\/UIqdNDcpTj0W1kkHAWPfV1VVo2KchnvufSvX41qYNfQdYSB6maZviya+Xw0H+190\/qhgW3FLUUMIbeMthWj1N0Zm52HD3SYoVbCPfPfeDp145jpfl3H5aE1WPvGebqnpWY+qTZeh\/TjsVDpxhwtfsPeQ3nn1lvPcUMlZO9cXi37D4aXz\/QjGKNqThmT8MRYN3SnDrqgp0bXIIRT9vhvOmvYc\/\/\/nPyHvzEWzY1wgXXnghUlJTMezVL+vdQRdy2ASB1fe4JsIPR8wkmYSEa194Hu7adxgjX96EksJz8fmGhsjXDo6wHiQRFKQFQHabWuX1rOLdOwhwysEuJExLTfUFwbLGzWe8\/IbvpGPd2nQ\/s6geU1\/K2a0diMH30qbM1kew68pjx0HpeeX7WEbav3cg+9syFfaEp8p8Mu1OnPzuXD8I0vMuQrM+d6Ns7n349sYn8OHnm3HTvmVodm0RqC6+vWwJhpfuTbiBFu3Bnigbw+2+M2YkoxRjcNNLRvwME1KgjrYd8sscvPNiC0zodXSi1ZszH+pPBUAvHn0pMyFRquyaXwoU5nqjgTL9fM+EaL4eMk3MEr4QppaCJKavPJZ6tecqBLjWRpVU3uvp+VzCG2sRgiSOicpnU76lRRLv4\/2VS9Fo4T3oMOYlbPrfEjRYPRv3bMrGoZpUTD5jFyqu9hyrm\/Hagxj7dSaG3HQNWm5ZhRHvVhqSafPgaBM8VPkxJRlj4efdfh7GVtcokuV3bI4XFpyO1fNl1YhWg1JgvjIfeNsucXftPoWkIsmYlnnk3kvYUF\/ve28Nmi0vvM8LaY3R26S\/D0ZepivW2igktJRfh\/a52IJSbHzmETRe9Tg6TPsQ\/1y5CNUvjUPLO59Eiytur3PARYMLrlKxEc8\/8jXGvLSuDskkrojE4ue9zAuDwZQs4cGNJPMGNg1MMl+IUS\/JqD+K3S0YySYA+YOAUll1st6HzbLgCRXx9TY5Kdd5m3JBkm1RJOM8jHMvXp+M+AUWf\/otnt3ZVN3b7V286KKLUFBQ4EcyEooIWwP2kHCWqHx1PsyQzElf26eNqSSzUxcn3d0CZaWrUV5+hWph1W0rkPZseViSrCpnLtI65AGlHsnFe+T3RtqzHY4fGa2Eum0Kv3hrG0PllGNYSTI3cfgD7ZkMtGWLRCvR4u+zfckYg99IMgBi+Lj5lrOxvmEDvP7Z97iucyus\/780tfg8wauZMRajfn9GeTnyykm6utdtVTnoXQLccYfn\/YDM7p4RUwwsWrs21HgO+\/348TnYPDcHi3aFX6ZYu3Y+shADV+b42hiqUv1gcTdx+O1I5jcltGkA9QDdWzsZY\/Abknk7VmLh8\/bFIefj8rNaeAKcrvEsLquDWiz3gQalCm3N3R7cHdLbMyPT70MNZifvw22TXiY7\/eZRo9Dg9jd9bXRSJ9O6icPvhmRWEiZjDH5DshCjK997LK3Y46z31uyyq4K2kdJSD8FoVJR7p4M5nPSh2mQt4+3yclzxTjm6z+6u2uj2chqHP5C6GGxTtpVkyRiD35AsxAjjjg8uSBfQvuDdBKzf69ln\/e1vGL\/qa+zffwO6fzEHM5cW4dUfmmLOnEqgpAAZGev8YuPr4Osx8wM9l7p0B0k6MRYXz0eXEU1x9Yvf4bNu3VQyvTx9LtO7oAAfN+uChyeV4\/bUt7BixTiVvnnr1hgxfz6yO4SeM4YbtdwKrdM5mRhF9BOjkjEGvyFZGD\/jcgBFrpdoXJTmgRRyr1TJsjIM+KICJ33TESXtU7Hsi9n477zLkbK1Hc59+U8BY+M7DR1HV38JPyeu+zyMoWfP8bjpiZ8i+\/UjeHT3bqwaP17F4GcaWuHoYbxhzx60nPU3fLz0R7R8YiTy81Nx4L778KqDdZtQG4aDwfnmm2\/WsS4yPaWV1booixLHIWjD6NnYJDEkCxNnvx3u1G+4Gdjr2Ui7iFqRWnBsx\/3qEqDJ4r14auxXeP75R0PGxncbBJUdyLDgPCxv8uSFONByFt7pkI7sWzej86mnonPTplj2xRdIS8vDuhcbIWt5JaYN\/wBjxgxA586dVdSojl7pxe8I5BmtL2O79eAJRDJCJ3EfRUqyDie+fWF2Y1ySGZI5hJ1kKymh5PJm5O53rz8Zn8gAPOmkRbj00lTHsfGdhvOWyMZ6vq1bU1BQUIp9LX6GVq1aqTZxD+L27Q9jypSrVcMZ7ThFO4T+h6FDccHo0fppUSqdSBT+nugBTjlwnMZdtFsncwh\/UiY3JItCt+nHrYqkCfcACickC6ceCSsuKqMcqsegrLJZlfOcs8aMwfc9emBWv35+YQdkvwuJqpvS+V1O4y4Gk2RR6IaEKdKQLMJdEauY+eHWI\/5iwUhGCO557DEsopfz6NE+b2jxL+N7q08ZB47TuItGktXjXfiR4lmsYuY7qUdi7VtJRnWxb9++mD59us8wsnP0aKzr3l1JMlr0xEva7hBAnWThxl2U0zhj7e4R74MujCSLEMN0IKXIaMTMd1MP22Mnydi+pUuXquYyLmK30aPVPnxeQjKqiXaHjeokY\/pw4i4GkmSs025ZgPPAaFoXY3XQhSFZhEiWTMUEizfhNMaHk+8OVG+gAEZ2z90YXIYMGYJ4HnRhSOZklJwgaYORjJ+Y4v1Oq7HjeD8\/UiRzanCJ90EXhmTHO3KSMH8okkUr7mKkSObU4CKSLF4HXRiSJSFJjrfJ0unDhw9HN+82LGuZoXbHu2lDoIhNTtVFIVm4BhdGllq5ciXiddCFIZmb0ZLkeRg1ieHJ6PISj8tqXXRLsnANLnZzMuaN1UEXQjL+oE2ZMgXilxcP7K11xsxpMxE+NtZtINH4F8tL\/NisJIuXdTFW3x5KPY9VO+zqMSSLJ\/pRqDuRB1sUPtdXZCJ\/tyFZNHs+DmUn8mCLJhyJ\/N2GZNHs+TiUHY7BJQ7NinqViRqimx9uSBb17o9tBfE2uMT2a\/1rS0SjR8xJpptZrV7FVVVV4J\/1Yhx7u1j2ki4exoV4DqRw6uaveqwNLuG0K9ppaFEk0RLJshhTknFHQFFREaxeytxcywHx3HPPoUuXLqjMykLLAwdQUVGBjIwMFaueril2RKvPv9rRHrDJWn4iSrOYqYuUYtzZ\/vTTT6NRo0bqKB7uUqebCN\/9+oV1+GH4cGRmVqCgrAy5W7dic1UOBuZUoVOnTopw1kskI9dFYnWsT7IOvvrQ7kDLF\/H+9piSjC78VsdKOjySLD2rcnDW2Y2RtawCLVsdxHOX5ylscjZXqZBql6elBSRZrN064t1ppn57BBLVwpgwJKMX8fjx47GlVy881L49Zs9Zjg4d0lGRmYmOHTsiz5DMcCsEAoZkIdRFkoxztvbt2+O\/9l+IMxpsx+C0cmRlZeG8884Lqi7WJ0kWDedI7k9cvnw57r33Xr9hTPWeA5cqvjim0j\/s+eefV8\/27NmDwsJC34GD9JeT0AvB+CB+d0wjZUci5n69J1kww4eAQ9WRUmvd7gzkpe3AKaecgry8PFx88cWGZCF+xY\/HOTIQyWbNmoVVq1Yp4nDuzDglPDlm3759GDVqlAryM2jQIHTwRuKy8ya3a7ZdfYZkEVJUAu2UDmYl5K8jJZyddTFRf7nChStRnCODkezkk09G48aNwR9Aplu8eLH6PNkQ3KtXL78zsQMZp\/icXuIkJyXmRx99pDzFRfLpJFuyZAnuv\/9+VZSkEenHfNR6MjMz60jeRB0PMZuThRp4gda79AMYAnWgri4ynBz\/5JB2FbeRrvde3\/t437NtDG\/Hi4MiEZwjg5GMR+SyjcOGDQMHP5dU3n\/\/fUUyCfw6Z84cP0JIP+naC6WdhP\/u2rVrHfVUSEZ\/tMcff9wXNFakI9sghNYjiOljwpAsFMtcvLcDdeJRn44JE4Baz\/HTkLCHiXLP9hV7I4qy\/YngHBmMZH369FGEoMvKxo0bce2116owA0Iyvdv0GPp8blUB5UeFQWKtc0CdZCSUzA9Z5uHDh9W8T+oM1F5DMhckCpXFDtRkk2SJ4BwZimQk1yeffIIePXogPz8fM2bMwA033ICZM2eiuLjYZxTRfzSEZLIBwUiyUKM5Qd8n6i9XuHBZB2U40aii4RxJkulWQrafKvj69etBScaLcymGG6f6TpKxHVbrokQO42YDufR5uLzn1CCQJBO11MzJwh1FUU6X7CSLMjwnVPG65dKqlloJnWhLOglj+HAzIgzJ3KCWnHl0acujrPS1O0OyKPapIVkUwU3CohN1PJwQkixYRKgkHCumyS4RSFTHzaQmmXF1cTkaT+Bs9drVJVr9apw2o4VscpYbbPNCvL4oqSVZvEAz9RoEnCBgSOYELZPWIOACAUMyF6CZLAYBJwgYkjlBy6Q1CLhAwJDMBWgmi0HACQKGZE7QMmkNAi4QMCRzAZrJYhBwgoAhmRO0TFqDgAsEDMlcgGayGAScIGBI5gQtk9Yg4AIBQzIXoJksBgEnCBiSOUHLpDUIuEDAkMwFaCaLQcAJAoZkTtAyaQ0CLhAwJHMBmsliEHCCgCGZE7RMWoOACwQMyVyAZrIYBJwgYEjmBC2T1iDgAgFDMhegmSwGAScIGJI5QcukNQi4QMCQzAVoJotBwAkChmRO0DJpDQIuEDAkcwGayWIQcIKAIZkTtExag4ALBAzJXIBmshgEnCBgSOYELZPWIOACAUMyF6CZLAYBJwgYkjlBy6Q1CLhAwJDMBWgmi0HACQKGZE7QMmkNAi4QMCRzAZrJYhBwgoAhmRO0TFqDgAsEDMlcgGayGAScIGBI5gQtk9Yg4AIBQzIXoJksBgEnCBiSOUHLpDUIuEDAkMwFaCaLQcAJAoZkTtAyaQ0CLhAwJHMBmsliEHCCgCGZE7RMWoOACwQMyVyAZrIYBJwg8P\/1mP\/RqehBUwAAAABJRU5ErkJggg==","height":101,"width":168}}
%---
%[output:08d8b851]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('snht', 'C:\\github_trend\\aerosol_trend_analysis\\snht.m', 103)\" style=\"font-weight:bold\">snht<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\snht.m',103,0)\">line 103<\/a>)\nsnht: la serie est constante, le test ne peut pas etre applique.\n\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('multiple_breakpoints_while_one_snht', 'C:\\github_trend\\aerosol_trend_analysis\\multiple_breakpoints_while_one_snht.m', 52)\" style=\"font-weight:bold\">multiple_breakpoints_while_one_snht<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\multiple_breakpoints_while_one_snht.m',52,0)\">line 52<\/a>)\n                    [~,ttt, ppp, PrctDiff]=snht(x, param{1},alpha);\n                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('change_point_analysis_SNHT_V1', 'C:\\github_trend\\aerosol_trend_analysis\\change_point_analysis_SNHT_V1.m', 191)\" style=\"font-weight:bold\">change_point_analysis_SNHT_V1<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\change_point_analysis_SNHT_V1.m',191,0)\">line 191<\/a>)\n            result_m = multiple_breakpoints_while_one_snht(data_m,{strcat(param{i},'_',param{j})}, 12, 0.05);\n            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"}}
%---
%[output:7fa5c406]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the number of abs coef does not allow to understand the instrument type"}}
%---
%[output:927a2d76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:152ac259]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:939222bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:87395408]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1b978798]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:929a6dcb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9775063d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:59336d8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1fc14498]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:201c41fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7e53a1c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:951efc0e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1bb09c71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3a59b96d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9df110ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7ed81e2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7b77b11c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:94982b7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0efde159]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9c0bd84a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:75f06b72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7ae89966]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1eaa267e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6acfa036]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6100d8ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:112bdc49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9c2ee56a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:65d976e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:16e76ae8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:07769fac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9f34e255]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:617c23c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:10fb6a65]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8110d059]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:17ce46f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:498574bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7c5b7162]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0a5bdfc0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7c3b09c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4ef3699f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:72605d9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:79b6a064]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:61eae4cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2600e63d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9482da4d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:394e1214]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:67a68fd2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9f8ff993]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:849f1652]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:041fcd47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9597f9c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:511ca6fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:18d48aa6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5f47d5f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8c64e6ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1c4d4d42]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3bea7d55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2a5c3baa]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdQAAAEaCAYAAACoxaaoAAAAAXNSR0IArs4c6QAAIABJREFUeF7sfQl4VUXSdhFCSNiTsIQdZFEWDSYqiAiKLK6ooxAQBxRQYVhkkH1GVEYR4y4gMIgK8ynojDqCijqiIsuIEJWRxQUQwi77EiBAwv9XY5\/07Zylz3Lv7Ruqn4dHc293nzpv9T3vqerqqlJnz549C9QIAUKAECAECAFCwBcCpYhQfeFHgwkBQoAQIAQIAYYAESotBEKAECAECAFCIAAEiFADAJGmIAQIAUKAECAEiFBpDRAChAAhQAgQAgEgQIQaAIg0BSGggsBTTz0FM2fOtOxasWJFaNWqFQwYMACuuuoqiIuLY30PHDgA\/fv3hzVr1liOrV27Nlx55ZXQp08faNGiBeTl5cGIESPgs88+g5o1a8Jrr70GTZs2Ncbv27cP7r\/\/fvj+++\/ZZ506dYLnnnsOKlSoYPRZtmwZmw8byjR+\/HiV26Q+hMB5iwAR6nmrerrxSCPgRKiiPEOHDgX8Fx8fr0SofCyS8rRp06Bdu3bw7LPPsv\/HNmXKFLjpppuMS3z33XeMLJF4sdWqVYuRbpMmTYw+r776Kjz++OOm4yONHV2PEIgFBIhQY0FLJGOJQMANoSIxzp07F9LT010RKgLVpk0bRqRImmhZYhs8eDA89NBDpmTJP3zllVegY8eO7M9Tp07BhAkT4O2334bU1FRGti1btiwReqCbIATChQARariQpXkJAQkBkVDffPNNRnxi2717NzzxxBPw4Ycfso+ffPJJyMrKCiFUJNjZs2dDSkqKMbSwsBA2btwIf\/nLXyAnJwfKly\/PyBjdt\/feey\/s3LmTWaeTJ09m3504cQLGjh0LCxcuhLJly0JBQQGcOXMGBg4cCKNGjYJSpUqB6BLOzMyEv\/\/975CcnEw6JQQIARsEiFBpeRACEULAiVBRjLfeegvGjRvHJEKX7e233+5IqFx8ef7mzZvDkCFDAPdCL7roIkbEuJ+6a9cutif7448\/ws033wyHDx+GpUuXMjfx1KlToVKlSrB27VpGxvv374c\/\/vGP8PDDDzP3MzVCgBCwRoAIlVYHIRAhBOwIFROWoSWJViRaqOKephiUZGahonWJAUvoot2wYYMxtnHjxszKRVcukuEbb7wBl19+OXz99dds\/xTH4d7qr7\/+ygKSxGuiDLiHK1rKEYKJLkMIxCwCRKgxqzoSPNYQUN1DLVOmDCNCtE7R\/aoS5Sti0atXL0au6M4ViXHSpEnQs2dPZoWKBLpnzx4jmpcHL\/GAJu4+vvTSS2MNbpKXEIg4AkSoEYecLni+IqBKqHiEBY+oNGjQgEHlhlBvvPFGeOSRR6BatWpsrOy6xaM0uE+Kx2n4UZmDBw8y9+7mzZuZe1fsc\/HFFzMLl893vuqO7psQUEGACFUFJepDCASAgCqh4qXE4y8qhIrkOHLkSEA3Lz+\/ivOIwUW4R4pu3D\/\/+c\/MvYzEiXusYpAS9hk9ejQMHz6cEWyPHj1g4sSJkJCQEAACNAUhULIRIEIt2fqlu9MIAaegpGPHjjEXLbpmjx49Cl26dIFnnnmGHWHhiR34HmqVKlXYvimSHf4X90jxiAwSZlJSknHX4vEX3CO95ZZbWHIJ7I+RwDzSmJ85xSMy9913H7su7rH+9a9\/hX79+mmEIolCCOiLABGqvrohyUoYAk6EireLiRbwSAsSKydP\/FwmVH5s5ueff4ZBgwaxwCJsDz74IDtzKkbkzp8\/v1iWIzHqF8eJiR5w3xTloP3TErYA6XbCjgARatghpgsQAucQcCJUtCb\/85\/\/MPJDC1WFUHHejz76iLlo0aJEVzFaoOIZVzkrEo6RXblyKkLsI5Mu6ZEQIATsESBCpRVCCEQIATd7qNwqHTNmDBw5csTSQsV+SKToJn799dfZneA+6IsvvmgkYhDPnfJb5RG\/\/G+c429\/+xv84x\/\/MNDAKGPsh9HC1AgBQsAZASJUZ4yoByEQCAJuCBWPqbz00kuASe+dzqGicDt27GDJ7vEcKjaM5MXMR3jsJj8\/n1m97733HvvOypUrHrHBfrR\/GojaaZLzCAEi1PNI2XSr0UVAhVDRzYpHXzDlID+qokKoeGeff\/45I1G0NnEsZkbi+XfFRPdY0QZTCVatWjUEEPGIjRy0FF3k6OqEQGwgQIQaG3oiKQkBQoAQIAQ0R4AIVXMFkXiEACFACBACsYEAEWps6ImkJAQIAUKAENAcASJUzRVE4hEChAAhQAjEBgJEqLGhJ5KSECAECAFCQHMEiFA1VxCJRwgQAoQAIRAbCBChxoaeSEpCgBAgBGIagT0TV8CprYchoX5lqDGhbUzfi5XwRKglUq10U4QAIUAI6IPAwblrYVv\/jw2BajzctkSS6nlFqNu3b4d33nkH7rjjDqhTp04gqw0P0WPeVcyhKiYkD2RyH5PoKBfJpKZQwolwUkNArZcO62lb\/0VwcO46Q+DkPi2g5szOWj471VA173VeEerXX38Nd911F7z55pshycP9AHjy5EnAXKk1a9aExMREP1MFOlZHuUgmNRUTToSTGgJqvXRYT7KFWnf29ZDUo7GWz041VIlQgQjVz1LxP1aHH7Z8FySTml4JJ8JJDQHrXriHeuyrbVChfV3m7tVxTfm9R7JQfSKo66LQUS6SSW2xEU6EkxoCar10XE8oua5yqaFKFipZqH5WSgBjdfwBkUxqiiWcCCc1BNR76bim1KXXlFDz8vJg7ty5bF8TS1ClpqbCbbfdxkpR8WobKPrZs2dh+fLlMGXKFFi1ahXr169fP+jTpw8rR6XSyOWrglL4+uj4AyKZ1PRNOBFOagio99JxTalLryGh7t27Fx588EH45ZdfWLBQRkYG\/PDDDzBnzhyoV6+eUQ8SRV+2bBkMHjwYWrduDT179oT169ezfu3bt4fHHnsMKlSo4IgFEaojRGHtoOMPiGRSUznhRDipIaDeS8c1pS69hoS6cOFCVvh46tSp0KFDB0NCrMvYv39\/uPvuu2HIkCFw6NAhRqZolT755JMGeS5ZsgSGDRvGPsMakk6NCNUJofB+r+MPiGRS0znhRDipIaDeS8c1pS69ZoRaWFjILFB0486cORNSUlIMCdENPHbsWIiLi4PJkyfDmjVrmGt3xowZ0LFjR6PfsWPHYMSIEWzsxIkTISEhwRYPIlS\/y8XfeB1\/QCSTmk4JJ8JJDQH1XjquKXXpNSNUO8GPHDnCLNPk5GRGqO+\/\/z6zYl977TVo0qSJMRT3VdE6xT3V2bNnh5Cy2fxEqH6Xi7\/xOv6ASCY1neqMU0FBAezbt0\/tRsLcC5MooCxVqlTR5lw6yrT7ww1Qdn8hlGucAuXb1w0zCtbTY0IdnlRHxzXlFxgtj82gKxcJddSoUcwyfeqpp2DFihXwyiuvhAQq4c3jdx999BG8\/vrr0LBhQ7JQf0dAx8VKMqn9XAkndZy+++47eO655yAnJ0dtEPWKKgIYA\/P0008zUtVxnfsFRztC\/fnnn2HQoEGMONElXL16dUaaaF2aWaFvvfUWi\/yVrVc7CxUDoVCx2NLS0tg\/rw0Xxe7du6FGjRqQlJTkdZrAx+koF8mkpmbCSR2nzz77DEaPHs0e0rVr11YbSL2igsDKlSvhxRdfhGnX\/wWufa43I1S7ZyemctUpnasKaFoR6rZt29ie6IEDB2D69OnQtGlTdg9BE6oITN++fZkV7LXl5+cDRivjC0DZsmW9ThP4OB3lIpnU1Ew4qeO0dOlStu0TZDpRtatTL7cI8C23x3+7FTLu7wRwfyPbZye6zXHbL5aaNoS6YcMGdoQG90Nwv7RZs2YGjnaE6sXlK77N+rVQT5w4AXv27GFWrk65fHWUi2RSezQQTuo4ffrppzBu3DgiVDXIotpLJNTMyy6D6u\/daPvsJAvVg7owsAjPmKJl2qBBA7YfUrdu6Kb5\/PnzKSjJJbY67k+QTGpKJJzUcfr444\/Zs4MsVDXMotlLJNSW+bWh8qgMKLyvoXaFRfxgFHULVUzY8MQTTxQLOsKbQ0XQsRl3aqaHshpehFNs40SEqqY\/HXrJhJrUoxHEP96KCDUo5fAApObNm7NzpFb+8oMHD7LEDlgiTcyKRIkdrDVBRKG2Sgmn2MappBMquv\/xTD4mwRFbuCxy3ELDvADYatWqpRTsqbaCzhlGmBEP91DRQiVCVUVOoR8GXiCJzps3D2655Rbm7pUbun67devGEjYsWrQIhg8fDu3atYPu3bvDxo0bKfWgDc5EFAqLUNOKF6Q7dd2VZELF4EzMGNemTRsYM2aMAQonpgceeCDkczXUrHshme7cuZOd\/ccTC25OUKhcmwhVBSWPfTAydsCAASx3r1VDouXKNUuO37t3b5Ygv1KlSkpSUGIHJZjC1omIQg3aWMTpqU9+hdwDJ6FeSiKM6Wp\/HlwNBedeiFNJJlQkNIwfsTouaPWdM3LFe5g9G7l1jJaqSOhe5scxMqFWnXYN5HdOhbhZv0LcnnxIqF+Z1UmN5Rb1PdRIgkeEGkm0i18rFokiGojFGk7zVu2GwfM2GFCN6dogIqQaLUJFosPIYt5ESxELfdx7772QmZlpGAOcmDD5BJ6Xx4Z9MI\/5rFmzWGpVbLLFGaSFKMqMBInXnjRpEjzzzDPMArYibztSd\/vb4M\/f55sNgTat20Dl0RmQO+1rODl2tTFVjYfbxjSpEqG6XRVSfx0ffiiijnKRTGqLLdZwQjJFUuWt1+VpMK1X0bE3tbt23ysahIpuUdzP5IlkzNyynDjwfGxWVhYjK\/FoDydddK\/yvVCzecR+fC73KJ07xy\/KbHZ92d3LrxMkqcsGDerv174L4fS7W43bSu7TAurOvsHLbWoxhgjVpxp0fPgRoaorVUf9xZpMsoWKZIqkGu7mlVC9uqc5EQ0dOpQRJW9mni9+dh5dpSNHjmRxItxtajcP9hWzvvG5RSzT09OVcpfjGKtrySRvRah4fVkmr3rl94KZkm54uT976Zct1Lqzr4fkPi29XiLq44hQfapAx4cfEaq6UnXUXyzKhCS1bOMhaNe4SkTcvXyNu91D9eOetrLWuHWJdZo50fLP0J0rEyAnOe5u5avVivxEa9HK1YyEiIGdItHjOCtClF8CIkmoGOV73ejuzOW7a9cuSHhtG5xcsQsqtK8b0+5exJsIVf3Za9pTx4cfEaq6UnXUH8mkpj8vFqof97S8dypLKbtleX\/5cyROzD+cnZ0dUj3LjJjNkBCP0nCXsRWhWr0EyKQeSZcvEurVWV2g2rRrGaHicUidssyprT7zXkSoftDTdK+SCFVdqURealjpipNfC9WNe9pNgI64Tymf51SxUDt37syOzIhWr6gp2ZpFQixfvjwsXryYBTnxACdVCzWSQUlIqOjyTerRmAhV7eenZy+K8o2uXnR9KOv2lkw4qa1TLxYqzuzVPe1EhNyFyy3I3NxcljIVLVFs\/Aig6h6qXQ5zMwuTBx6lpqYyMsY9W\/x\/jCiW933lPVQz4g3XsZlZA59iLl8d17nayrPuRRaqTwR1XRQ6ykUyqS02wkkdJ7cWqtrM1r3kiFlOOiJhWkX1clLjZIhjeAAS\/0wMXuIuYOwnnkU1iwhGubDxwCeUYcuWLexvqyhf8friSwC\/VpARvngtsyhf3V5m\/a4PIlSfCOr48COXr7pSddQfyaSmP68Wqtrs1r3szqFaZTES91MzMjKY1YhWJJIXHp\/BZnU0RkwHyKWSUw\/Ke6gioeIYs3OoeAZW3seNROpBLruO69zv2iBC9YmgrotCR7lIJrXFRjip4xRpC1VNMvteVq5jP3PLFqr8tzw3Ej\/2McvC5EcOu7FkoYYL2SjNS3uoUQL+98sSUajhTzip40SEeg4rcb8V\/7bbQ7XKEayGuvdeRKjesdNyJBFqdNVCRKGGP+GkjhMRahGhilG+ovvYLDmEn8xLatop3osI1Stymo4jQo2uYogo1PAnnNRxikVCVbu7kteLCLWE6ZQINboKJaJQw59wUseJCFUNKx16qRLqnokr4NTWwzFZfYaCknyuNB0ffnhLOspFMqktNsJJHSciVDWsdOilQqgH566Fbf0\/NsSNteozRKg+V5qODz8iVHWl6qg\/kklNf9E6NqMmHfWSEVAh1G39F8HBueuMobFWfYYI1ee61\/HhR4SqrlQd9ReLMi19eg8c3nYKKtdNgKtH1VBXgI+eRKg+wIvCUBVClS3UWKs+Q4Tqc2Hp+PAjQlVXqo76izWZ\/vfWQfhg2DYD9KtH1ogIqRKhqq9zHXqqECrKiXuox77aFpPVZ4hQfa40HR9+RKjqStVRf7EmE5Ipkipvl2Qlw80v1VVXgseeRKgegYvSMFVCjZJ4gVyWCNUnjDo+\/IhQ1ZWqo\/5iTSbZQkUyRVINdzsfCFUs1SbiKace9Iq1fEY1nOdTiVC9aknTcXRsJrqKiTWiiBZasYgT7qFuXXEM6retEBF3L39pLMlRvlYZjazyBbtdr3KFGacC527nl\/tzuZ9vNgTatG5jFBineqh+kY3SeCLUKAH\/+2VjkSiigZjuOL2+C2DnyUKolRgHAxsmRgMids2SbqHa1V91U5vVTEFWZO13XrvFwJ+\/WA+1ZX5tqDwqAwrva0gFxqP2C\/J5YSJUnwD6HK47USQmRo8cRGh1xmkVpMKkLYWGuAMbJDJS5VG+Fav+BpfdMh\/iytWGpAuH+1wx9sOjRah21Wa4lZeZmWnUP+Vu25ycHFauDRtWmxk\/fjxgxRcsCI6NFwXndx1E+TR+7SZNmhjFx7FE3IABA2DQoEHAa7jya4YjcT+fWybUpB6NIP7xVkSoYf2VhHFyItQwgqswtc5EoZPbScRJJ2sQa1dOP1QZPj0UZ2i7W1oC3Pm\/\/JAo38tumQeXdZsPSRc+qESqWPA798BJqJeSCGO6NlRYSee6RINQ5dqiZpYef87w\/Uir+qhYto3vhZrNw8kN+3nd2xT3YMV9V7OC4ohpON2+VoQaN+tXiNuTH5OZkeTFSkFJyj9f8446kgR\/2OhWvFdHrHSWycoa9LlkPQ3nOMkyTWxWDuIm7QmJ8r2w7edw7b0vQtm6d0L5S5+2vd68Vbth8LwNRp8xXRsok6pXQj349mNwZu8WiK\/WAJJ7PKKMhxXZmL2o8+ovWOB75MiRIBYOt5sH+\/Ki4yiYWWL79PR05bJrZgXQ+bzytfBzTuw9e\/aErKwsZWxUOsqEWnXaNXDs6DE4OXa1MTzWMiMRod51l\/FWqLIInPro+EAmQnXSWtH3OuqPy2RmDSKBRaPJVvPqQ2fgsirxzN0rR\/kimSKpIpkiqdo1JFMkVd56XZ4G03o1U7pFL4R69MvXYe+0e435k7s\/qkyqVi5YMxLin6E7VyZAK7eqk3Vo52qWC4zzG+SEWqtWLUBy583KQo0UobZu0wZqf3Ar\/Np3IZx+d2uRPvq0gLqzb1DSv46dyEL1qRUdH8hEqOpK1VF\/dtYgulij0Zxw+lOPNXD2xzNQWG83PD3u\/6BM1dZK7l7ZQkUyRVJVaV4IFckUSZW3itfcA9UGn9vXdGoyocn9Zbcs7y9\/jsQ5evRoyM7OBtzb5E2VzMzcuEERqhOpO2Fk971ood7wcn9I6tEYcqd9HWKhxlpmJPl+iVD9rBBNk9AToaor1Yko1GcKrqedNRjcVdzNZIeTG7ftjF9PFosQxj3UZRsPQbvGVZTdvXyNuz02I1uoSKZIqirNTQSsuP+J1qHoxlWxUDt37syKhFu5XmXiQ0IV66HyACcrC9VJBjlYSQUfpz6cUGcNfAquG93dKOCR8No2OLliV0xmRiJCJZev07oP2\/e6k1e0o3xP\/PQCFB7fAYUJNWBfuTu0in4UdXd26wwmJ4\/kVXXbLth9CiZsOG6sLx4h7HXBebFQ8Vq4h3pi\/ZeQ1PwaZXcvjlMlIU5iubm5MHXqVGaJYps8eTIkJSVZBv7Ibli+Dzt79mxISUkJgUmWRQyWSk1NZWSMLl50N48dOxZkl6+VNezmpcGt3jih3jdgMowb30PLilhu74kIlQjV75rxPJ4I1Rq6\/G3\/grzvRhkdTlbvCymtxkK0SZ4LxHWXcno5FPz4sCEnRvKO\/LYrbJ19GCrllYYj5Qugfv\/KpvugSKZIqryh+9rPnrBXQvW8gAFAjvI1C\/qxiuodOnQoC\/ThZIhycMuVfyYGL3HSw34iqZpFBKNc2Pg+KcqwZcsWGDZsmCmhYl8nOf3gZDaWE2q71P7w9C3XUmKHoAGO9Hx0bCbSiIdejwjVGv8pq7+E7QdyIa1wJ\/Q9NRNOJXeFKpe\/oB2hVtz9LMTt\/dC4EQw6+uLrEfC\/R\/canyGh9p5Uv9jNyhYqkqmfPeFoEKpIRPwGxfOjVlmMxP3UjIwMdg4VrUgkSjwWg83qaAyS5cyZM0PwlFMPynuoKoSKE0Yj9WCFSx6EW5tfCn+ZcAXodhLB7xOS9lB9IqgjSeAt6SgXyWS+2GSi6Zs\/A3pUrwKpze\/RjlBlCxUjef\/zdGvl5Pi4hypGCPv5+UWLUP3IjGPDkTxBtlDlv\/3KHMR4Tt7H2o2GM1UvhAc7pEGfSxK12trwe59EqD4R1JEkiFDVlaqD\/gZ8d4yRDG\/XJeTCkFoVtXrQyHuop\/etNCJ5qdqM+noLJ6EiYaHFi43vobZp08adcGHsLRPq7RdXgb9eU0Wrde739olQfSKowwPZ7BZ0lOt8lskqG5BsnaIuxzeIg8thv1YPGjvdUbUZdw+RcFmoYpSv18xK7u7EXW+ZUB9LKQ03966r1Tp3d0fFexOh+kRQR5IgC1VdqZHQn92xEjlQB5MlTG0Wr93ekhNOvNpMpR3PwJV9T7rOQqSusaKesery9XKvJWGMSKgt82vB\/11XC\/I7pxKhxqpyKSgpuppzeihHQ7pIyGR3rMQsUKdLlUKDUN\/cNAd2Hd8FNcvVhH4tHogGROyaKjj5yUIk3hhPtF+5boJtKTgi1KgtB08XFgk1I60Z\/Htia+1eHD3dmDCILFSfCKo8aHxewtNwHeVSlUl0j+LNe0mcrgqaqkyq85n1c8oGJAfqcJm+y18Nz\/0w2ZiyX\/P7o0aqTjihzpp+8mdI3\/a+Ia9qFiKRQCvXSwhJtH\/1yBohpCr2vXxoZXCb2MGPHmmsPwRkl+\/QRpXgnhtTyEL1B2v0RpOFGj3sVa0cmXxEid0kTle9UyeiUJ3HqZ+bbEBcple3z4Qvdv\/HmPqGBrfAXy5\/1OlSgXwv7\/na4cR1dnPeJ3Dnd4ch72RtKJ+4A67Nbu2YhUjef5WFvyQrGW5+qS77WO57+ZDKkHfxdzBixIhA83MHAiBNUgwBmVCv37AfJjRLggaTroPD2d\/Cqa2HjTEJ9StDjQltYw5FslB9qixSD2S3Yuool4pMsntUvG83idNV8VKRSXWuoPpZWajjL38UbmxwS1CXsZzHbM\/3wQ41Ld1zXGctNpeFG\/5bwZhXti7xCzn1oBwhLAuFZIqkik3u2+TWJChz4zoi1LCviGAuIBPq2M+2wC3NKkDl6xrCvsFfFrtILFaeIUL1uVZ0fCDjLekol4pMdhaqm8TpTmoNOs0fn8+usLabvUF+4L3H1wth54lCqJUUBwuu6e50W4F8b7bn++ztDQ1ClWu08v7X\/7cCtNxc1pBBtC7xQ7PUg22\/ORHi4q3XtjzkrsgzJWXZQu30VFXYWeUbItRAtB7+STihdrjwGkg+ch3c880uKH1FNUhoUAlOvL2pmADJMVh5hgjV5zpSIQmfl\/A0XEe5VGUS3aN4814Sp4ugyWQXdJo\/eT6zwtoyGcjWm2i53VMTGHnJtUd7lY6DMe0reVoPbgaZ7fnimUEzmTDTUbNjR+CNcVug7m9loO6eMsalROsSP7RKPcgjhOu3rQCHt52yTRIxfMJk2J23C9LK14TJ44fTHqobxUa5LyfUyYm14aJfbjWkSahfCU5tPVJMulisPEOE6nORqZKEz8u4Hq6jXNGQyYzsMLE7fs6b3zR\/mINXnM+ssLZd8gPZcrsnDeDWxIMg10NtsCQPHiofbxv56nqhWAyQ93y57mSZcPgNO89AxRHn0udhw+hctE6vHlUjZHar1IN2QUkiKX+0ZSFMWlW0h3zXBX2g1rb6JdpCFUu1iWDKqQe96j0aqQezU8pB0zWhFX6QVMs0qGzcRoX2dWkP1atSIzWOgpIihbT5daJBqGZkh9IFSagyaZsV1pYJVSQKOVMSHpsZVOVwMQv1iun7We5bHqQTTm3K7mmuO9lqRhmQ6FE23mRXryinHNFsZrlj\/60rjgFarCIpP7HqUVi0ZaEx3bVpnSFzf5sSS6hmSfDx5q3yBbtdD3J1m3DWQhXlfrbpLmi242o49WMXQ+RYdO+a4U0WqttVKPWPBkmoiKyjXNGQyYzsCo9vhxM\/vWjAqGqh2u2T4ndiOj5RRzJp4D7h3e81Yl2cMiW9+M0p+OyH41B9fT60+NdhRqY8SEdlHXjpY0ZyeEQFXb54lOeJ3OpQUCbdmNqMUK8b\/E5IiTcrOfBF452jr8GJlP2QdCAV7qh4r+ULA1qoj67bDIVxaRBXuBvGNqoFcT+WKbGEaldKzW+ZNSuy9juv3XrjLwJIqM0PtoT8b7NY99evqAkHLq4GzW5v4qoerpe1He4xRKg+EY4GSaiIrKNc0ZJJJjvZahUJ1Yo0VfZJ7UgDScrMgpP3FcvvPQOPfv0TNB5ai53PWzXlMCx9Zg8bahY167QWrAKh7F4OzNzTnbKrMULFozyLjl0KZxK7WhJqs04\/QYes0bB6QU84ur8GVL30augwIZP1f3XdzJBEFfOWvgPPH9lukGTblSeg06HezG3s5DJG13iDn76MOKHyyjEcALHaDLfyMjMzjfqn3G2bk5PDyrVhw2oz48ePh1mzZsGaNWvYZ+I8+DdeZ8qUKSHFyZ30LX\/Pr92kSRNYvHgxuxaWiBswYAAMGjQI5ELi4UiLyGUSCfXE8buh5X8vhE871YVJzaobYg9OLA337zkXlBbViR5bAAAgAElEQVSLR2eIUN2uUKl\/tEjCSWwd5dJFJitCLbX3g5CapGJwkco+qZVO7HLdzty0BabnVgkZevPSf0C\/tE2wL2E0fDZmnykRO+kfv\/9g+q6QsmqckPHlYMljOYzsKqbuYcSVdOFwY0ozeZvemmRYqJM37YRT5fqGECq+CORViwf877AarwJ6Ab547UGjD177aPdv2B7oqXJ9GIHeVO9yiCvYA+\/tb2b0u+DzTXDZ388FNskvELJrvFUFgNt2RJZQ5XqoZpYeJw6eT9eq7iiWbeN7oWbzcHLDfl5z84p7sOK+q+zu5QoIp9tXJNR8uBGuzBgMExpWhIUbjhn6x7OpeJyGt1g7OkOEavJkUjkCwYdFiiTcyISyRUoulQe7W6wOvv0YnNm7JWz5YK0ItWDDw5bBRSr7pFZY4L7hx3P2M7IZ2CCx2L7guwevDXGhXvvjmzD0yOeQs+lJ+OX9EyHTqlqpGKm7cNg202Ms305bBB9PrG3M26JbHtw660pLQj36XC0ok1EOkk4ehj83S4Y7cnJhx6kiq6JWYhzsPFlojL+39AdQ+a9xsPOnliEvA2v+OBsW7DoF+RVHWy4b0X0s7sWaucZxkou2fg3fTPyTq8QObn9LTmRjFpuBxIufY8HvkSNHMquQF\/+2Ii0zkpODhlCW9PT0kILjdr9BswLo2N+KUDmx9+zZkxVDD7KJhJpe8SSUqdEP3trUJcRCRTJFUuUt1vZWiVClFePWtRckcQXpbgxSrqB+VCoyBZUP1k5mVQtVDi6y2ye1up7Z2cuBDRON7hho896+ZiEu1E4\/jIHBNZvAnn0dYOY3i4z9RbQK4zOS4MqOlUCcw+zaeDY06ZOvYVf5jobliEdukKTMkingvi4PAhK\/39KhPHwzKNW4BCbvF0vNmV07ffMpuHD87pCvcO\/3m0bPw4y9LULuVR6PwU1IqtjEvWbZNc7HlV22Eva+NkiZUN3+vkX5rFywZiTEP0MXq0yAVm5VJ+vQztVstf44odaqVcsgdF0I9XTuZWwftWffi2F3pQR2C0imaUfy2RlVbLF2dIYIVVqJKq49kfhK1R8YkkGGf8entTvoL17a7ocuyxRXrg6UrXtHiJtO\/kGpkFdQRKk6j4pMe6fdC6cP\/BviKsRD4bEzcGZXFUhqcU1g1ionxTP7vzbEPlFnjFHM+9BnVzN3JWJcpdNS1Vuz7Gd19pIPyP7yF3jzbLWQ8W1+XQ5jrm0Oq3ZMhZd2ncsgc7ps1xDLDi1dO1L9Zc3\/weqdq+CpxMeMufmYbya8B5\/NPBcUJTe0gMV8ukimSKp2rfTpNSEWdqPPlkPmK+fSBWLjxGj28jCxWTn4adGr8H2ZOtDq9Hb4dc+tzNpFax6DsLhFbmWh1v37Eli38s\/KhKry+7a6V5nQ5H6yW5b3lz9H4hw9ejRkZ2cD7m3ypmodWrlxzeR2S6hOpO7nByFbqF\/vGQDbD7WFyZ0aFJv2npU7GamSy9cP4mEeq3JsBh+4YgSofKZQJr5S9R+AQ5WyWACJvAfHb8fsoL98q\/IPPT61DVS6ah7rJl9TZV6RvM5unaEUcekFfjfuMxVCzct5BvJ3TDNEObnmEJz83yH2d3L3RyG5xyPAr7n9ZCq8tOVWqJeSaBsdKL7kiMdlkDTPVrvJUn+i3sQ5VF+SUGZ0987YctK4Hzz2giTCm5nlhRZav9NHYdX1r7NcvrjniIR6tnRasXl4Qgj+BbpfkWhxPf1pf1dYU\/qyYmPWvTsD\/vfoPvjg6j8a1iuSFza0YOtdVcHIXrT+9hOwNutCY45SBbtD5MAvZEKtsHsHVP2pCsQV7IaE43Phwovqw2P3jYDpPyyDmfuK3MA4LmvPQrhzwbNs\/i8uugumdCzSPcqEySz4MaEXv\/o\/eK3gZkOWyts+hGsfOghzKz6nTKh+XPduImDF\/U+0DjEYiZOnioXauXNnViTcyvUqE59Ish06dGAYYZATWsdjx44F2UJ1kkEOVvLybJDHyISau\/Nm+HulvvBxsyIPCB\/D91LJ5RsE8mGaQ4VQnd5g5e8Lq90ER9MeYoR6KudeEC0ffhuclO3Ix4w0+TgrQrWbl5NXyunlUPDjwwaiMrk7EaLTURHx5cPpxUGJUKUkCWf2nIRjn55zH2L1kkrdbgoJHEJCxX9WifOtsMP5EL\/Szf5meBis9lDN5nC6Vw74fTkbYdWRqgb+SHhIqty6NLNQa\/ywFe75PgHKDPkYnt6+13TPkZMyErLc0BJNO7sTHt2aEvIVv\/YPX\/8KFf\/3IXzcZYDxPZIXP5aTu\/wYzCsoZGR78IIv4UitTEaiSKZl8j8JCUhy+qnGn\/wEEo9lA1bK2RL3R3ZMiDecr+KeHdB+03ro98MwRqZIqrxVW38Snq1alpE86uDldRthTtmBxvfXfPhPuGLGevhr9feVCRUHe3Hd4zhVEuLklpubC1OnTmWWKLbJkydDUlKSMc\/QoUND9inlfU2+Dzt79mxISQnVpSwLkv2KFSvYNdDNfNdddzFMrAjVyhp289LgpHsrQp1yw1Fo9ts+QJfv+yf+ZGqhph05Bddv2AcPP5AOyX2KXsLcXjPS\/cnlKyFu5l4tk9oGuFViRagycYnTlk64GRIaNgqxfM2ITSQnPh77yZl9+He4x4cNZeINLdvS5epAYUIN2FfuDqi4+1mI2\/uh8T1aZfx+8P\/FsTml+8K\/j2QZFp9srcsWmyyvWYYgEQeZUM3I2oy8uJVabfBrEFd5XUjg0MpDF0Lv78eAVeJ8WV+iPIjf2Wo3G4Qqexj4HqrZHHb3yo+qbG+zHF5JLmW6Z8jdrw98sxxW5rUIWYVIQt3yNsN1Wy+CMbWawrG0ogAi3hHHo2tUJCj+HRI2NrPv8POWa4\/Bnrh42Nu8aC8XyWvC8VIsYOof3xyDZ\/POBPYsSjg+B26tfhDK5PUu5t7mF8laNRlgWw146w\/3hlyX4zRl9Zcw+2irkO\/ue2YJJC\/50jWh+rkxOcrXLOjHKqqXEygnQ5SDW678MzF4iZMe9hNJVY4I5jL06tUL2rRpA+LfVoSKczrJ6Qcns7HcoJk55Cy0yDsEX\/7nTnZ0ZvgfmsL3tSvC\/YfnQM2Ccy\/Pu0qnwd8r97V8UQ5atqDmOy8JdfYjreGKS+owksSGhMXdhwNLjYba5daa4ouEsmX7RqiRF0pQxytcB\/F530N83rnzZMVaYRpAXGiQhh2xieOxH+4lyuOxDz7UC45vN7WK8Xs8Xxl3are1XJKg209WBSSoHSerwlWNqkBmwZyQHlxm\/FB0nfJO4veh91Cb7feKhGp1RMWMUAuPVAV8KeHuXpnI39l9FVS49BlGqk55e0VZcY\/UjOTlBA1mMpkFLOE6+mlFR3h2fXNm3R1o8iOz7swaBvd0q5nA8tvKDV2hXeIXwMlf\/wxfXFBUvUXsh+OxmQUI4XdojVoRKj\/iIl+Xu6Pl4yl+HzZ4P\/ft2wV7WtYPOSIjztti5S8A26rDujuL0s\/h9x1OA1x3STlTnK7+9Ge4csq8iBIqymQXHGSVxUjcT83IyGDnUNGli0SJx2KwWR2NQRKfOXNmiBrEIzCcYDGK2A2h4oTRSD2IhHpxqWPw47\/aQMo3HZiFGl9vNTx64KmQe0R38NFOYwCLYsRKO+8I9cu5veD+683Vg2TSuspPlrpj5HZ8e6zoVis5Zfeq7B7ne8ZHlvcq9oIgvnwg2Zm51bllLhI9t6i5i6\/wxPYQ\/fExh\/PLQ0qrscD3mhE4eZ9UtNblgCWRcF\/bPB7+ke6\/KkyNUjthz9laEdehShSvF6F6nNoJm6vXs4wQtiJ5\/BwbvqDIDV3GbV7+Aj75\/mlXLl8v8gc5JujkCXYWKhKsLk20UNHle2JlDSjYdwHsScuD0lU3Q60zoUbHwvJd4eeuzxOh6qJAWQ67\/TRdZS5Jcu1K7g\/l6vSGKkfegrNbQ9+48T6R4MzI0isGIhGjhWzn\/jW7Nidc7skQLWP8rkzV1szyFufFqNpPynTzKnKJHZdWuBN2xwX\/gnDmp9VwPPv+85pQudVstoeqK6Giy5cHG1ot+kdTxkDFa\/oSoer6VDCzfnSVtSTKhR6AM+VbwVVl3or47SEByhZqEELgvNj4iwARahCoqs\/R5O2FkPPJI+c9oYpRvhiQdPToUeB7qupohrenaKFmNgbAgEPczsJ\/SK6Jl1Rhx+V4m1G6N9Q+1RvuG3l5eAULcPaYcfmePXsWli9fznJbrlq1ClJTU6Ffv37Qp08fKF\/e\/pwc4kXWaYCrhqayROCvedNhRZo+braSrKpo7aHqjqnVOddoyy0TqigPkmt8jaIgOf5d6cR7oXKXCdEWXfn6MUOoy5Ytg8GDB0Pr1q3Z2az169fDnDlzoH379vDYY49BhQrmwRscCTt3nzJa1JEQsEEAA5KeSpzomAiBQAwGAYzyTdv3JozauyqmLNRg7j50FjlwKah6qUHKakeoVtcpVdgZkm\/7e5BihHWumCDUgwcPMjJFqxQj4Th5LlmyBIYNG8Y+u\/HGG22BIkIN6zqiyTExwWsPwoxOPUOOoxAw4UOg9\/Y3odXep+GBqaXOe0INH8rBzeyFUPPX1Ya409kxU2w8JggVFYGu3RkzZkDHjh0NDR87doyVbsJDzxMnToSEhHPn78waEWpwPwyayRyBfx0eAC\/XGUzwRAiBMScfgaobFhChRghvv5fxQqh4zVhy+8YEoc6fP59lHBHTdyHQuK+K1inuqZplExEXABGq358DjXdCIJYCksxSCDrdH37vdZzK3G77EKG6RSy6\/b0SqlPCmOjeVejVY4JQcX8AQ8JfeeUVqFYtNJE4fvfRRx\/B66+\/Dg0bNrTE9vDHQ6Hg1AfG9xhZJkaU6aQUkiU2EQiaUK3OZrpBBwkwrvBcgXJsKb9cyHLvlr7oAOwoe4ObqQBz5+JcB+v3czUuXJ2JUMOFbHjm9Uqo31f8C3S8tihNZnikC2bWmCFUVIaZFapa1Z5XMYlPS4Qzu08yMk1oZB\/IJEOMkWhx5eNNiZiFf+edsfzebqyKKr2+AFhdl8trFlmH8uA4bGbf49hTm47B7ip1oEH9ouLA\/D784KSCha593BKqnbWHqQAbLskLKZ3m5b6ve3sDJNZ\/HvY3+RFSf7kImnx0G9Rq+gPLSCTmxeVzX3rqQ9gVf2mxM6NIwq2nroOKVX+DL+5IhoIyraAwrkaxRPleZPQ65qGtYyDtwKfk8vUKYITHeSFUfJbsKuwL6Q9MjrC03i533hDqx6+8AE0\/+bOB0qH2LU3JwApGsepJubZVi5HxxpXxsHVreahfPw8aty7Kg4oLAgl8+U\/1oWbBHmh8xZmQsSpEyc9rWb0AYMrAOon7QkQXryvLhB2REI+v2AcHMi5iOIjWOr\/XnfFpcNEtxV8gcNzUvXew641o80kx0t28LpHNZ0a2HKdGrQtMv7fDw+9LibefiPqoZ0r\/AT4qV1SIAEdaWZlYkSXl10lwuM6NcLjuTcUucuHcvdBi5c\/w7UM1YcsFFxjf8wT1Kb9cxEqmHa8WH5JFCIm4+vp8+K15Wfbfiqm\/QeV2GwATK\/Q9VZRMYxhMg7UV24Zct2\/+DPjjgWksMcXTVZ8I+a7L9+9Ap3Ub4ej+GvDS7dXYNfH6YgUcK6QwNzH280rAbXd\/zbIsyYkhsBB7p7WPwNBFFUtsUJJ4vlTEN9xRvGKyfUzoH0TzQqj4nMK0o5jHOxZazBOqqst33qrdsGbGWMjM\/x5yyraCRq3PwB1py011hErEBztas9g4IeL\/Iyliiiw8hCx+P3b7A\/BB+a5wc94nMLnOTPYdH4fXw0TP2DAB9FUXbWEkhARRO3G\/I7HLKRFx3HtHO0C9qodg5aGLAAk1+6LZxr38c2NrSPn2R3afeN0\/lV0Cwy6YGUJ8W7ZWgLe\/z4TWlQ9B+rb3jfsRZa7f4BhkpecY83LLFKu74LyPHHiKfS8TPebWRblEmTjBI04HL8mCGw88ETI3XgSJfPP6RMAXFpmM8Tt8KcF2deYO196FSPwYFxeUhcnlZ4bUBk1bs5JF\/SL5iK3Wtwug\/pwWsKPLKci9KTRXKZJi+oJdcMElM+FofAuY1+WPxtCyR7NZxZcmH93KrE25ADgW6E5efRSOVCiAA5mV4Kc+RVskSJhIqqi\/\/1QeBHn1ahjz8u+2\/toY6qYegofSngwp\/YbfN938Geyq8UeYmnKrMpxo2SYdHgGn85PZmILU+4oVC7BLL5hesBqeP34f\/Pn0\/bAmZVDIda9b\/zy0W5gND\/1cs0QSqpwEn9+8Vb5gZaU4dOTzY6J+XiEniLn5vNPvOQ2XtyqeStLsGvgsxrzZWGkqFlpMEGoQQUmojKGzl8FbS9ZCrabpMLj6NFtCRSuMN0zSvC7jQfbntoMn4fpfZzJi5k0kTE6anLg5kWKyeRybe6CoNib2\/UPa8hDiWZpTm1l33NXKM\/xcnrDYuN5\/tzeCP278i\/E3T2TfusqPjMjwgdm6fgU4c+YMdGqZxq7ZcuPDIQQ2+sf+0CS9NyNbdIeL98plntf2DRCvi4v7hbU3wccNH2DdL9n2PkxKnsYIUGwLTj8EeyrcBLdVegs2\/\/Ile3FAQkScoPOD8FCVlXB09qBiLyWYLQXTjWF75rI5IS8ASKYP5z8JdZMTodOqcdC90UrjhQaJfn+TdlA7cV\/Ucy0\/VPZx+C6hyOJE6wxJJb\/iaAMirL6CdUKT8rpBQdX74ED5ovqo2AnroWblHoFruz0Lv359FqYenQBb2x1jVWdwvrjC3XDpxm1wyT8GwOFtp5j7llukWILt6N3xsPaCfMirXyNk7pZHV0Df0t\/At\/FYxaMhYAUXTK7f6sxqRrT44nao2evQLf7ZYrVUq+z9GdoUloHCC+vDp4fijHtBMiy39wwUHD8FF3RKYfL8UrHoe\/4CcMWuVvBjtYqwr1LDYuXf8CWg7e6V8EZWx2JHjrqeXgCdTm+G7FNXwN4qoRb10M8HQ8qqd0ssodqVUgtXmTXxPGu4CPXZprug9VWJhtGBi4kbJ\/KW2qqTvaFrj8djgUuZjDFBqEEcm8GbNd68\/vwctCz7L0tCRQtr08p4w5rFKif4AMKGlu7geRuKKRiLXMtkiZ3wc6yCwsc\/9cmvsGzjuYLZ2No1rgLDGrwPmPR9V1w6jPy2a8h3OE7O8oRkWLbuHex6OB7\/i3LxhgT7z\/7NjLJk7\/1wiMmM10HSRZLmD1Uc89vMSxkRoaX7WoN3jXmHN88JKe+GZDp13x0st+ayjQfZNbnFfUGzZEZ2PL8tl0W838vqJEKfSxKh9Lvj4MzXRekH0bW8q3QNRrjpAyezuW88OBSurLPJuCfUCQYn4LXxXip+9pTxUoPjsCpF9kWvmlbBieSvEd2luJfKW6vND8PGSisgL\/kNo6Zo+YO92dc3NLgFmha0hydOZISIiATTO\/NbaNF+CSu9t3pBT3gvrS3897qrjX63Vj8Ad6xNNgqB4xe1LlwL9dqWh46Tz1mQWG1GrGYzoEYcDGleyZhDjnzfX\/9v0CT9bvjn158Xk+nuyqfgjuQ8mH6ocgihIvmjvCv6VoSMe2vAvvcPM1JFgkfLHF8Aqh\/NhTl718NrLcaHlG9DMq63ZSNU+mo\/XF\/5OMxp0rZYUgwk02e6ZLD7kKvnZL37Glz0y\/iIE6pdtRme+D4zM9Ow7rjbNicnh51UwIbVZsaPHw+zZs1i9UuxYUFwrBjDm2p8iMr6FmXGYuN47UmTJoFYSFwsTffuu++yKjjhsFCRUNMrhhoW4j1w79\/LJ7Jgb7UHKJevioLd9OGJHbCIt5gVyU1iB7yeWGD8+Im10Opo6D4RlwldDC+sz2TEh4TFyVAmCSTLzb8dg\/Qa8dC0VhUY+e8txW4LCQAJ1W\/jFVPQihDJEOeVSR6vefvFVQxCTUxMBE5s8v3IZC3WPJW\/Q8sXv0c8zK7pdJ+8VFqFnz5hFipvn13xJCws1zUE68+\/eCVEP\/wlAu\/N7KUGr\/1MxichNWf9Ys7Hu60yNCfhAfg+\/jJouWM5JHzxX3j7yizIa3q7IQ63UMdf\/ih0OPgVPPReD9hUoy77vlVBDnT+3zJo8af6LIEJLw4vBzzxUmv\/GL4cCr7dygKNLus2n710iBGRL371f\/DtqTTISNgND7a\/25BB1i2v9oMdZPLCCjRTm8Wz9fSvbYfhjYIGxjxYdq36+pPwxb0XG5+hy7r83oIQchxQfjf8VrFeCCkiGV+++UeI730B3HCoCmR\/8TOszbrQmAeJ+Oldb7KAFFkmvMa1E3+DxFqTYOrhbyLm8pXroZq5Zflzhpdjs6o7ioTF90LN5uHkjP2sSruprHFZZnFeq71YHBMtQuX3hJ5BfMF2eq6oYBCpPjFhoSIYixYtguHDh0O7du2ge\/fusHHjRlepB2VCZXUDf3oB8re9E+ImlAt\/OylCrKf54pJdjIS5tWpGxk7zef1eJky5zqfVvLKVIp75kr\/bU\/4maHbdVGMqK5K2upYo04kFT8GJ9V9CUvNrWJ1Ts\/bq\/42FmoXfG25s\/nKCFqpokeNY\/A5dlWZ1WvF7fEnKz33HtJoN1o1NKNwPcHi1qRxuCZVPgi7yTT+Nh6fTeoSQy9nDq+DBJpuhX4sH4PvXfoCPxp41rntF1lLI7H8UDlXKAnyB5GXuZMt3XOE3kHVdF4bDT\/+dzDwP6O6\/8MqxxgPI7mUJdfv3Paks0AcDlu6vsZ9hhO2dtV\/B3\/ZeYsj0cLX\/wU2Nr2CEKluoGBiE7YuL7rJduvgCgK7v9\/YX7Rc3+mw53Hv0Erj5pXMvE6PfHQ+LK3QN2W9GIh5yxUXFCJVbxkfqrIe3D46PCKFyIuKFwvkNiy\/qvLoLEhJ+jhbnyJEjQSwcbjcP9hXP28v1SvGaWDTc6dw9l83qWjLJy8rTgVCxfNvJXtOLGTRen5GRGBczhGqWHL93794sQX6lSkVuLDvQzBY+9ufWn+yuVFGAKnGpzBVkH1W55IeuWDjb7jsuq1zQ2+4eRJl47VG57igfL1uh6MZeOPhS9rX83ZiuDdiPTqxZKsrBX5KsknucqDMGyuavZe5Vs2ZXyN2u5By6yDcu7AsV6qaEHH+5Y8lEGNK8MnuR+GDYNvjfWweNy16SlQydsqsZ3oWCDQ8bxIekyqN4kWiuzW7Nxlm92Ni9LP1782Z4dGuKcd37a+yDPzVvzP42+67fBXWYTKsgFSZtKTTG4T4mtikdp9ku34nNyjH376PrNrMjN6VPfw+XzSwDA7o2hIzB587DfrRlISPOM4lF2x7cEn9k3WZ4\/7cieQ0LtdtPMHXxGFeEuvTpPcwtXbluAlw9qigwy+n3Z+WC5dYl5hjPyspi0\/DP0J0rE6BVPVQr8uNy2bmakQAbNGhgXJ+PwWeeTNL4ndWzkI\/TgVCpfJvTiozy906LyIt4qsTlZW4\/Y9zIZfdCYfednQVkJjuXKeX0csOVif3MvAKyFYpuH7RCjR\/873vRohfAijD5\/GZuzoIK6cwalGUS5ceXjDP7VoZYv2i14j421kM1KwuI+9HXfJ0NLTaXhRv+W8EIHOqY9wJkrXqKTZ\/c\/VHYmveHEAv1+gk7oHn\/aw1C\/fvazfDq0aI6ohh0hP+QeLllZ7VO7F6IZBcqJy6cy+y78Q3jmEyI06ubtzOXNgYz9d76FizKvwGev7go6IrLk7ftAHS4uDqgy3hgw3MR86+sngCrt\/6bnYntsLERdHgkE\/CFBRu+WDz2vxWwqdNVxi0hEaNsZnuo124+Bh1SV7P0o6rHSPAa+BLD29UjayiTqkxoMu6yW5b3lz+3qgZjRsxmuhWP0vD7tiJUq5cApyLnESHU5DpQUNAecnachNJVN50LWgQwYlcwODKoLTM\/z1I3Y2PGQnVzU1Z9iVCDQLFoDjsLyI5QK+5+NsQaNEst5mWP1qpEnzj\/oc+uZi5+JMQqnZaCbDVjcBh6KrDx\/0fSlK1fO9c4v3fcc0YXbJ0vS8GyD\/ZB8ul\/QZfajxrQ4FEAjJBePi0edv58MdsHvWrwGSjd7G8GoU76tbDYviMGASGZIqk6NasXIjlgiRMXzmf2XZcqhUwmWXeII+7vLjpdFFh15vgpyNu2n0UOZ7U4Aw3qNGYvHrxxV\/6uuFbQ7+7JMOPXk7DzZCELaKrz\/F7j5QOJ+Mne5yLIzQgVZWq79StXhGrmEXB6MeFyu4msFfcpMQhIdOOqWKidO3eG\/v37s8pa3OoVdS1bs0iAWMZy8eLFLMiJBzjpbKFW7D8djl3YFSa+sRPeP3DKuL20I6cg7Wg+pP6hGfwzq4nTEtfqeyJUn+pwYwn6vJSr4ZGQS8UlLAptZaGKbmaxv9s9Whxrti\/O5zezqEvVHxgSvGUFst3Lg5WrWdxzxheERa+8AI8eOGedYsPD6niGF+fGSF5MmpB8YUNo+2hXQyY8oiJG6tZ8dzc03JMPL78deq7V1eL4vTMSGR6bES1IM0Id2CAR7qkJhoXKA6X4NTEIS8y6xCzT3W+GHAfjXgL5RanrNY1gQ8WiLRtugePcHW5fDFfNOJeMRXZD42deCFW2UFVfTPB6TkTII2a5BZmbm8tykGdnZ7N74BGzqnuofB\/WbL9UlkUMPMKqXEjGuH+L\/48RxfK+bzT3UEddXhruuKYxVLr3PbamPnhjGzxyoMBYwhO2z4Ibtv4Cx2++Ci59+AUvSztqY4hQfUIfCeLyImKk5HKz\/2xlDYrWi5d7NRtjJpcZKYrWoN2+rt3Lg5WrWX5RwBeE7eveh\/ikBKiVGAej7h3IRH93\/BT4cXYH4zYuH1IZGt9dyIKSMEIbie\/dzcfgl1\/2wuGfzx2P4vvGQeElzmPn8kWZECe03gtPbGfWvhyBjIQ6+uSEkGNp3KKXXfmprepD+bpFe6MYNXztqz9A9SqroGX9GcwtjnvNiPEdB+4NyZY0quI3UH77cVcWKt4n7qFuXXEM6retoOzu5fjIEbOcPLTrHAEAACAASURBVEXCtIrq5aTGyRDHcMuVfyYGL3EXMPYTSdUsIhjlYuvi92M3KMOWLVvY31ZRvuL15XUUTpdvlzvbwviD1QFrnWLFsNPvboXXr6gJ39epCJfFfQEDzxSRKNd\/ONZ5OOYkQvWJaqSIy62YOsoVCZl4kBTiJQc8mZHi2Wo3m1peZm5oq5cH2UIV91dFvcmuVLT8Ku87BAuHbYOWm8saXZvcmgSXjYs3CBW\/cNpTdrs+7PrbuXw5yeN4jqdsoabkHYeHtoyGK+sVZSLD9HGVr58Cb4zfCh9+8RscKV8AKy45AbKFikFONyR8xJKbYMKOMim3MUv+pdXrQ\/aSMYPSS\/V2wxfrK7gmVL9Y2QUHWWUxEvdTMzIymNWIViQSJR5PwWZ1NEYuHo595T1jeQ9VJFTsb3YOFc\/AovXcpElxt2o4CXXUxVlw1Ueppmoom\/EWlKlXFG2P2yKxknYQb4gI1eevKxIk4UVEHeUKt0xme6hywJNMilwmeW8QMVc9QqW6l2xm+e3\/fisLxsHAJd46PVUVqnXIDyFUL3vKXtYNjjn65esw9ZsfYW2tdtBy5zIY2CAJkrqNMXWNI54T914SsoeKQUSj9rwNR5f+2ciGU\/Hq5+HXPbeGBASVur0sjJtxIbPAuev5zmXdIa7yWkN0JNSK7Z4vtod6c+Xf4K\/N68HHH38ccUL1iisf5xQQ5GV+2UKV\/5bnROLHPqrHb7zIZHbNu+66Cx7\/7VZomV\/bdMr4eqshMaMo6QuSaaykHSRCDWCVhJskvIqoo1zhlsnM9epUS9FqXxdxdxrLdaO6l2xm+W1YswOe+mQL3PfvZKicF8eOcvRf1sCUvLzsKXtZP5iKEkmVN3ygYQAJ7neJFir\/3irA6eDbj4WcNVYJCLJ6ObE6G0uEek4L4n4r\/m23h2qVI9jLWnEzhlvvdoSK8yVc9CmUrroZElt0gLpTi3KUu7lWtPqSheoT+XCThFfxdJQr3DKZWahmAU\/i2VkxKOnkss4hST6sgqXMdKK6lywHAqErV7ZQ5T1Ur2vA6zgkUzG\/M1oJZdr0tCRUPD8qni9FVzYmrZCbSkCQ1csJfv7yuo3GcZ0\/tWgM6K4nQi0iVDHKV3QfmyWH8JN5yeu6kgk1rnY5iP9DfSj17SHIX76r2LTJfVpA3dnuavZ6lS2ocUSoPpEMN0l4FU9HucItk106PStrslT9B0zPoYqp+LzqwIp4C4\/vMPZ30ZWrsocapAwqc\/H8zrgfXP2B70KOF2GglNieWPUoLPr\/pMob5ij+y+VFx4P453ivb4zbAnV\/KwPbqp9mAUFyWk\/sa7dXLR5lwvUUi4Sqgn9J7CMTatVp10B+51Qo+5\/9sG\/wl8Vuue7s6yG5T8uYgoII1ae6wk0SXsXTUa5wy6Sylyn3Kax2ExxNe6jY+UpVd68b\/Vglwpj+xGY4+FJRoXazPVQ31\/Hb1+3xIrRQJ60qIlDMUXxjg1uKieE3sIqfV8UIaUwUQYTqV9ORHS8S6nWju0Pl0RmG1+Nw9rdw7KttkFC\/EpzaegQqtK8LNSaEVheKrLTerkaE6g03Y1S4ScKreDrKFW6ZVPYy5T6lL\/obHChzVbFMSW7cvao6siN88SjH5UMrK52NVb2u235Ox4tkCxXnn\/5xF1iTtxvSy6fBoOs\/Nb2kn8AqswhpPBtLFqpb7Uavvx2hmq2p6Enq\/cpEqN6xYyPDTRJexdNRrnDLpJoKUXQpinuo\/Hyll5zOKnpSIXwd1pTd8SKzoCRV3PHenAKr8HssR4gFJkR3sNXZWCJUlZWnRx9OqHeWuRuGtrsSqk0rSrFJhKqHjlxJQakHXcEVeOdwE6qKy1e+KVkmN8n+vQCkErwUbpxU5LY6XmRGqF5wN5PBqugB9rU6G0uEqqJNPfrw52+71P7wxCUZRKh6qMW7FESo3rELYmS4iULVAhTvRZSp1N4PQgqqq55D5fMFRcbhxklFl2b7lVbHZoIiVKc9VjlCmvZQVTSpTx\/+\/B1ysgfcO\/luSOrROKpbG+FAhly+PlHV4eFndgs6yhUJmZwsQJn0RJmwVJpYT9VNYJIbt6fTkosETnYyWO1XWhGqlxcZFQvVqdIIEarTStLre06ozzcbArd+OELb7TI\/qBGh+kGP9lBdoRdtonCKXpUtVDeBSUFZaQhotHFyyuVrtt\/l9CKjulCc9lhl70JJd\/mKpdrEe1ctV6eKu1VJOdXxKv04oU67\/i9ww8v9Q9Y5Rvme2noYEupXjsnoXn7\/RKgqK8GmT7Qfflai6ShXtGVSiV51IgYrt25QVpoOhKqay9fnT8f38JJuoVplNLLKF+wVUH6dvXv3hpSZ8zqf1Tgu95wnZ8DVWV0MQpXPodZ4uG3MkioRqs9VE22SIEJVV6Db6FV5Zie3rhMZq0qqw5oy26+0cvmq3lfQ\/Uo6odrVX3VTm9UOdzGLkly3NWh9yTEsfJ2f+ev3cOLtTcblYjFDElmobdoEsl50ePiZ3YiOcukgk5voVRnXIN26dotPB5xk+XSVKRouX7tqMzzxfWZmplH\/lLttc3JymAWIDavNjB8\/HrDiCxYEx8aLgnPs8TpTpkzxZTXya2NFGV58HEvEYW1WvC4mq8c0hPXr14eRI0f6upbTA9WKUGULNRYzJBGhEqE6rf\/Av9f1oaxqeamkNvQKGiZ2OLztFEuOH+3EDrH0ghZpQpVri5q5ZTlx8Hy5VvVRsWwb3ws1m4eTM\/bzmntX3IO123dFmaNFqHgUi2dKitUMSUSoRKhen\/2ex8U6oeKNH1neC87s\/9rAwO3RGjPw5KTx0U6O74VQ5WM2nheJi4GRdvlyguOFwrmoZsfxePUXLPCNRCUWDrebRyY1s8T26enpymXXzAqgm0EcbUKlxA4uFr4uXekcanQ1URIINRxuX7msmVmB8ehqzj7y2OyYDebaDXfzSqivrpsJu47vgprlappWxbGS28oFy63Lnj17QlZWFhvOP0O3qkyAVvVQrYhWdAGPGzfOEE90EcsFxnknTqi4P4rkbtWIUINZrRSU5BNHHUkCb0lHuXSSiUfrFibUgH3l7jCt82m2NIKM5uXzyxZqtJPju7VQzY7ZTGxWzucvy3m4F0KVE\/n3a36\/MqnKe6eyhLJblveXP7c6omJGzGYomLlxY4lQzY7NkIXqvN6160EWanRVEglCVclWJJPiyep9IaXVWFD9UQcVzStqQ6fk+G4J1arAuFVe3qBWoRdCVS01Zyajm8hacf9Tjp5VsVA7d+7MioSLVq8ok2zNIqGK9VC59aqjhYoFxuVqM6q\/vaDWTrjmIQvVJ7KRIAkvIuooV7hlcjrWwnGU3bankrtClctfUCZUL\/pwMybcOLmRhfd1kkk+ZmOXl9fL9a1I3m1QkmqpObPrORHhM888A23atAFOYrm5uTB16lTIzs5m02FkbVJSEqjuofJ92NmzZ0NKSkqISLIsYrBUamoqI2N08aK7eezYsaCTyxcJFc+hUnL8oH4JUZqHLNQoAf\/7ZZ0eyn6lU93flIn3RJ0xkNr8HiJUGwW41Z1TXl6\/usbxXixUHId7qN\/tzYFLq2Uqu3u5vHKUr1nQj1VULw9m4mSIc+IxGjzSwj8Tg5e4Cxj7iaRqFhGMcmHj+6Qow5YtW2DYsGFaEipmSqJcvkH8CqI4BxFqFMGPwL6um\/1N7rYtqJAOhyplKe+hRgJBt+Slo0x+ap+quoq9EqpfvOzOoVplMRL3UzMyMtg5VLQikSjxWAw2q6MxSJYzZ84MEVs+AiPvoepMqLMGPsVcvjquc79rg1y+PhHUdVHoKFckZHK7vxmkTCr7tyrLLUiZVK6n0sdJJrN7d5OXl8vgxlUcLUJVwcuuj5Xr2M+8soUq\/+1n7qDGWiV2MCsJGNQ1Iz0PEapPxJ0eND6n9zxcR7lKskyq+7cqCo01nIK8dzeuYiLUotUk7rfip3wPFfd0dWlEqLpoIiA5yOUbEJAep4k1onBzm6r7typzxhpOQd67G1cxEWoooYpRvl4zK6msT699xPJtbVq3gcqjM6geqlcwdRhHhBpdLcQaUbhBy83+rdO8sYZTkPeO2Ki6imOVUJ30X1K\/589fjPJtmV8bKo\/KgML7GmoVv+AXe3L5+kRQx4cf3pKOcpV0mdzu31otvVjEKah7d\/NzJEJ1g1b0+8qEmtSjEcQ\/3ooINfqq8SYBWajecAtqVCSIggfHoMxx5WpD0oXDbcWPhExu8SOZ1BAjQlXDSZdeMqFWnXYN5HdOJULVRUFu5SBCdYtYsP3DTRSy6xGld0peH26ZvCAYizIFFeHsBi8iVDdoRb8v7aFGXweBSkCEGiicricLN1HIwTEoYNm6d0L5S5+2lDXcMrkGKQbd9UFG+brBiwjVDVrR70tRvtHXQaASEKEGCqfrycJNXmYWKpIpkqpVC7dMrkGKQUINMsrXDV5EqG7Qin5fItTo6yBQCYhQA4XT9WSRIC+39UojIZNboGJNpqCjfFXxIkJVRUqPfkSoeughMCmIUAOD0tNEkSAKt9ZSJGRyC1YsykRRvm61rNZfLNUmjpBTD6rNZt3LqqSc33nF8USoQaKpwVxEqNFVQiSIwq21FAmZ3KJOMqkhVtItVLMk+IiMVb5gNdSK9+LX2bt3r5Gs3+tcduOIUMOBahTnJEKNIvgR2htES+nETy+yG3WK8MU+RF5qa0JXnNyWb1O7Wz162dVfdVObVYXksI9ctzVoFIhQg0Y0yvMRoUZXAeF+KHuJNg23TF4QJ5nUUIuWhWpXbYYnvs\/MzDTqn3K3bU5ODrMAsWG1mfHjx8OsWbNgzZo17DNeFJzfPV5nypQpvqxGfm0sEbd48WJ2LSwRh7VZ8f\/vuusuVuWmfv36MHLkSF\/XctIaEaoTQjH2PRFqdBUWbqJwu39KFqr6egi37tQlKeoZDUKV66GauWX5c4bn07Wqj4pl2\/heqNk8nJyxn9fcvOIerN2+K8pMhOplFYaOodSDPjHU8UFzvhKF2\/3T8xUnL0tex3UeaULlBMcLhXMczV7UefUXLPiNRCUWDrebRyY1Preos\/T09JCC43b6NCuAbtafCNXLr6L4GCJUnzjq+KA5n4nCbbSpF\/2FOyuQF5l8LmPH4brK5GUPdcavJ2HnyUKolRgHAxsmOt67kwuWW5c9e\/aErKws1p1\/hm5VmQCt6qFaEa14\/XHjxhnyyi5isxvhhIr7o0juVo0IVXkZ2HYkQvWJo44PmvOZUN2q063+vOzThlsmt\/N76e8WJy\/XcDvGi4W6YPcpmLDhuHGpgQ0SlUlV3juV5ZXdsry\/\/LnVERUzYrYjyYULFxouYyvsiFDdrip\/\/YlQ\/eGnZZQoEaq6Ut0ShZd9WnVpzvV0K5Pb+b3011UmtxYqkimSKm\/d0hJgYrNySpC4iawV9z\/l6FkVC7Vz586sSLho9YpCytasuFfaoUMH1hUtWLSOx44dyyJ4yUJVUrOvTkSovuDT8+FHD2V1pbolCi\/7tOrSEKG6wSoICxXJFElVpTkR4TPPPANt2rQBTm65ubkwdepUyM7OZtNjZG1SUhKo7qHyfdjZs2dDSkpKiIiyLEj2K1asCInexSAkIlQVzQbXhwjVJ5ZuH8g+L6c8XEe5SopMbvdplZX2e8eSgpPb+3bb3wuh4jVwD3X1oTNwWZV4ZXcvl02O8jUL+rGK6uXBTJwMcU48RoNHWvhnYvASdwFjP5FU5YhgLkOvXr1CCB3\/JkJ1u6r89SdC9Yeflu45slDVlUrkpYaVrji5dfmq3a19L7tzqFZZjMT91IyMDHYOFV26SJR4LAab1dEYJPGZM2eGCCUegeEEiy5d0UImQg1C2+7mIEJ1h1ex3jo+aIhQ1ZWqo\/5IJjX9ebVQ1WYPXy8r17HXK9pZqEiwujRK7KCLJgKSgxI7BASkx2mIKNSAI5zUcYqGhaomnXWvoAkVr2S1h0qE6ldb7saTheoOL7JQfeBFRKEGHuGkjhMR6jmsxChfTCd49OhR4HuqamiGvxdZqOHHOKJXIAs1onDHxMsHkZfamtAVp1gkVDXEvfeKRCk2L9IRoXpBTeMxRKjRVY6uD+Vdu3ZBzZo1ITFRPWtOOJEknNTQjdU9VLW7c9dLDlwKul6qO2nMexOhBoGiRnMQoUZXGUQUavgTTuo4kYWqhpUOvYhQddBCgDIQoQYIpoepiCjUQCOc1HEiQlXDSodeRKgR0EJeXh7MnTuX5aTcsWMHpKamwm233Qb3338\/VKtWzZDg7NmzsHz5clYfcNWqVaxfv379oE+fPlC+fHklSYlQlWAKWyciCjVoCSd1nIhQ1bDSoRcRapi1sHfvXnjwwQdZlhCMTMMDzz\/88APMmTMH6tWrBy+99BLUrl2bSbFs2TIYPHgwtG7dmuW3XL9+PevXvn17eOyxx6BChQqO0hKhOkIU1g5EFGrwEk7qOBGhqmGlQy8i1DBrAaslYNV6zHfJEzrjJdeuXcuyiNx9990wZMgQOHToECNTtEoxmwgnzyVLlsCwYcPYZzfeeKOjtESojhCFtQMRhRq8hJM6TkSoaljp0IsINYxaKCwsZBYounExrZaY\/BndwFghIS4uzkj2jK7dGTNmQMeOHQ2pjh07BiNGjGBjJ06cCAkJ9kmuiVDDqFCFqYkoFECiajNqIP2OExGqMlxR70iEGiUVHDlyhFmmycnJjFDff\/99ZsXyRNJcLNxXResU91TNKjLI4hOhRkmhv1+WCFUNf8JJHafzhVB5LmBMns+r1qih5L6XWBzd7PiNmETCana5qDr2I0J1r4tARqArFwl11KhRLOgIz1hhaaJXXnklJFAJL4bfffTRR\/D6669Dw4YNba9PhBqIejxPQkShBh3hpI7T+UConOAqVaoES5cudSwqroaedS8kbwz+xKBQjGWRCdyswo7KNYlQVVAKuM\/PP\/8MgwYNYspEl3D16tUZaaIyzKxQrnzZejUTiysUA6EwuAlbWloa++e14cNv9+7dUKNGDVbrUJemo1wkk9rqIJzUcfrss89g9OjRYScZNYnC0wufWyNHjgSst4rPQszPa1cs3I8UnCyxIHmDBg0YscrPVr+E+uqrr7J7cFrn8fHxgP9iqWmVy3fbtm1sTxTfyKZPnw5NmzZlWAZNqKKC+vbty6xgry0\/Px8wWhlfAMqWLet1msDH6SgXyaSmZsJJHSe02HDbJ5KZgezKt\/HE95mZmYZlxwkoJyeHkRM2LN+GAZmzZs2CNWvWsM8eeOCBYkQpk9eCBQtg\/vz5SltcIoqizEiWeO1JkyYxkhYT6ItWJD5\/MTgUT1VkZWUZ0\/kl1Oeee47VaXVa51WqVGHbfrHUwk6oXEEiKGYLZ8OGDewITUFBAdsvbdasmTHEjlC9uHyffvpp4ziOXwsVF9eePXuYlatL6joETke5SCa1RwPhpI7Tp59+CuPGjYsYocoFxuVi3yg5f+bx+qZWBcexDip\/ETCbB+fiBC0XJ+d\/qyAly8znFK\/P55GftWbPXr+EinkHLrvsMsdnFFmoJtp1IlQMLMIzpmiZoosB317q1q0bMhO+kVFQkspPp6gP7cOp4UU4xTZOXvZQ90xcAae2HoaE+pWhxoS2agCYkBsfaBabwYkIXbPorsVgIu6mlUlSnAf7ii5WJGPRInVLZlbXkkkeZeCkLlqkZvemEpRkViyd9lCVl5r3jmLChieeeKJY0BF\/46NjM+4wJqJQw4twim2c3BLqwblrYVv\/j42brvFwW2VStYrXMCMiMVJWjni1qocqk5+V1WpGhkjgaJCIrln+7JRJWrSiRVe52f1xGcTgJLekbvXioeNvT+3XYN0r7C5fOwF5AFLz5s3ZOVIrf\/nBgwdZYgesCCJmRaLEDtbo6rhYSSa1nyvhpI6TW0Ld1n8RHJy7zrhAcp8WUHf2DUoXlPdO5UGyVcb7y59blVeTidnMuydeU9w6syJUq5cAmdSdrE7cd+WWMxGqhoSKG9JIovPmzWPuEHy7khu6frt168YSNixatAiGDx8O7dq1g+7du8PGjRsp9aDNY4AeykrPSBZpSOXbnLHSFSe3hCpbqHVnXw\/JfVo6AwAAsvvVbpC4TymSEY5RsVDxuYfJbbCZnTs12+vEnOaLFy9mQU6cbHmEsBypK7tf7Y4UypazX0KdNfApuG50dy1\/e0oLwaZT1CxUjIwdMGAAy91r1cRDzGbJ8Xv37s0S5OP5LJVG51BVUApfH10fykSozjrXVXduCRXvFPdQj321DSq0r6vs7lUhQh4xywknNzeXxX5kZ2eHEKPKHioOwEhgq+AjOfBJDDzCFK0YnYt7tvj\/ZvPIbmO7wE\/xfvDoIh4PtCN7q9XEZZ52\/V\/ghpf7E6E6\/+z07kGEGl396PpQJkJ1Xhe66s4LoTrfrXUPOWLWzFqziuqVI3XxKtxy5CTLg5ecztfL18Uz+9h44BOO37JlC\/vbKsqXX58TsN35VvHZiXvCfgh1zpMz4OqsLkSofhaiDmOJUKOrBV0fykSozutCV91FmlARKbtzqPwZIx8NFPdTsaoWWo1oRaLFh8dXsPG9VqtgJFlLInFjPIkYlCQSqiwzP4eKZ2DRev72229NEziI1xODkyZMmMC267C4iV2TzwZTlK\/z7yymehChRldduj6UiVCd14WuuosGoTqjZd\/Dag\/Vz7xohYoWqvy3PDc+C7GPSg50P3KJY4lQg0JSk3mIUKOrCF0fykSozutCV90RoZ7TnbgHin\/b7aGqWsDOq8JdDyJUd3hp35sINboq0vWhTITqvC501R0RahGhilG+4lEds+M3ZokXnFeBvx5EqP7w0240EWp0VaLrQ5kI1Xld6Kq7WCRUZ7RLZg8i1BKmVyLU6CpU14cyEarzutBVd0SozrrTpQcRqi6aCEgOItSAgPQ4ja4PZSJUZ4XqqjsiVGfd6dKDCFUXTQQkBxFqQEB6nEbXhzIRqrNCddUdEaqz7nTpQYSqiyYCkoMINSAgPU6j60OZCNVZobrqjgjVWXe69CBC1UUTAclBhBoQkB6n0fWhTITqrFBddUeE6qw7XXoQoeqiiYDkIEINCEiP0+j6UCZCdVaorrojQnXWnS49iFB10URAchChBgSkx2l0fSgToTorVFfdnS+EytMMigVDnLWm3gMTQ8ycOdOoUiOP5N9bzShX1DHrR4Sqro+Y6EmEGl016fpQJkJ1Xhe66u58IFSe2Qirai1duhTkHLnO2rPvwVMhNmrUCDZt2mQk7BdH2VWjUb0+EaoqUjHSjwg1uorS9aFMhOq8LnTV3flAqPjcGjlyJGB5OCQ2u6owzpos3oPXecXKNHgds5JxRKhqyEatHqqaeMH2IkINFk+3s+n6UCZCddakrrqLBqHaVZvh1l5mZqZRGJyXWsvJyWHWHzasNjN+\/HjAii9YEBybXKEGP5PLtC1YsADmz5+vnNSej2\/SpIlRfFx0G8t5fZE4sfqNXNScCNX5N4I9iFDVcLLspeODBoXVUS6SSW2xEU7qOEWaUOXaomaJ5vmLO8+Xa1UfFYmLu2+tEtbLxcitipNbIcYJFUutmbmKZSPDyuggQlVbk0SoajgRofrEiUheHUAiVDWsECcvhPrUJ79C7oGTUC8lEcZ0bah2MQCwIjMzEuIExN2ovHA4XsxuHnS58qLj2Je7Y3mZNbOC5nY34NRftkitiN0pKMnMupbloj1U5aUWGx3J5RtdPRFRqOFPOKnj5JZQ563aDYPnbTAuMKZrA2VSRXKbMmVKsaAdTkI9e\/aErKwsNjf\/DN256enpIS5aq3qoMtFakZts8eL1kPDEAuP8BjmhYhQukrvYrIhdtsL5\/Pj89FM\/lQhVbV3HTC8i1OiqiohCDX\/CSR0nt4SKZIqkyluvy9NgWq9mSheU907lQXJJNN5f\/hyJbPTo0ZCdnQ24t8mbTMxmZdfEa4pWoRdCdXM\/5PJVWiK0h6oGk3UvHR9+KK2OcpFMaquNcFLHyS2hyhYqkimSqkqT3a92Y7j1h\/uk8hlNFQu1W7duMHbsWHYJOUDIzGJEwhProXKytbJQ7YqMm7mJiVBVVggFJamhZNNLx4cfEaq6WnXUH8mkpj8\/e6jLNh6Cdo2rKLt7USInIsRjLXikhRNSbm4uTJ06lVmiIjGq7KFif4wENjvCgt\/JgU+imzY1NRX69+\/PXLzobkZill2+Tt462a1MhKq2JikoSQ0ny146PvyIUNWVqqP+SCY1\/XklVLXZzXvJ+4tm1pxVVC8nR06oeAUegMQ\/48FLVvu1XCr5ui+99BL7iu+T4vgtW7bAsGHDTAnViSDNjtPQHqrzyiFCdcbItoeODz8iVHWl6qg\/kklNf9EgVJTM7hwqt\/zkqFdxPzUjI4NZn2hFYpAPuoWx8b1WO3esiIxI3EuWLAkJSrIjVNWjN+LLw7vvvstSE9o1ea9Y7ktBSWrrOmZ6Obk5vNyIjg8\/IlR1TeqoP5JJTX\/RIlQ16ax7WbmO\/cyL5CdaqPLffuYOaiwRalBIajIPEWp0FUFEoYY\/4aSOk9ugJLWZw9srXITKXbIoPd9DxT1dXRoRqi6aCEgOItSAgPQ4DRGFGnCEkzpORKjnsJKjfJ3cr2oIB9uLCDVYPKM+GxFqdFVARKGGP+GkjlMsEqra3ZW8XkSoJUynRKjRVSgRhRr+hJM6TkSoaljp0IsIVQctBCgDEWqAYHqYiohCDTTCSR0nIlQ1rHToRYSqgxYClIEINUAwPUxFRKEGGuGkjhMRqhpWOvQiQtVBCwHKQIQaIJgepiKiUAONcFLHiQhVDSsdehGh6qCFAGUgQg0QTA9TEVGogUY4qeNEhKqGlQ69iFB10EKAMhChBgimh6mIKNRAI5zUcSJCVcNKh15EqDpoIUAZiFADBNPDVEQUaqARTuo4nS+EytMMYq5fs+ozaohZ9+IFxK0KhTsVGJcr6phdiQjVr5Y0G0+EGl2FEFGo4U84qeN0aKXgRAAAFBdJREFUPhAqz+1bqVIlWLp0Kbz55pusqk1QjWduatSoEWzatKlYAXW8jlMyfRVZiFBVUIqhPkSo0VUWEYUa\/oSTOk7nA6Hic2vkyJGA5eGQ2JBMeVUZNaTse\/E6rzgnXsesZBwRqhrSVG1GDSfLXjo+\/FBYHeUimdQWG+GkjlM0CNWu2gy39jIzMw3XLC+1lpOTw6w\/bFhtZvz48TBr1ixYs2YN+8zM3SqXaVuwYAHMnz+fValJSUlxBIqPb9KkCSxevJhdS3Qbm5Vpw+o3sluZCNURataBCFUNJyJUnzgRyasDSISqhlU0qs3I9VDNSq3Jxb+t6qMicXH3rVXJNrnUmmrpNY4gJ9SFCxeauoplr52VF48IVW1NEqGq4USE6hMnIlR1AIlQ1bDySqgH334MzuzdAvHVGkByj0fULgYAVmRmRkKcgLgblRcOx4vZzYMuV150HPtydyy3SM0KmtvdgFN\/lFO0SK2I3SkoySqYSZSN9lCVl1psdKQ91OjqiYhCDX\/CSR0nty7fo1++Dnun3WtcILn7o8qkiuQ2ZcqUYkE7nIR69uwJWVlZbG7+GbpY09PTQ1y0VuXbZKK1IjfZ4lUhVIzClfddrYhdtsJxfrJQ1dYkWahqOJGF6hMnslDVASRCVcPKi4WKZIqkylvFa+6BaoPP7Ws6NXnvVO4vl0zj\/eXPkchGjx4N2dnZgHubvMnEzA0AK7lUrEJuoZoRqpv7IUJ1Wh3nvidCVcOJCNUnTkSo6gASoaph5YVQZQsVyRRJVaXJ7le7Mdz6Q3eqfEZTxULt1q0bjB07ll3C7NypTHDiXmmHDh3YOCRctI5xHplQraxfHGfmJiZCVVkhRKhqKNn00vHhR+SlrlYd9UcyqenPC6HizLiHemL9l5DU\/Bpldy+OcyJCPNaCR1o4IeXm5sLUqVOZJSoSo8oeKvbHSGCzIyz4nVng04oVKxj5opv5rrvuYkFIVoTqtP0lu5WJUNXWJFmoajiRheoTJyJ5dQCJUNWw8kqoarOb95L3F82sOauoXk6OnFDxCjwAiX\/Gg5es9mu5VOJ1J0yYABMnToRevXqFEDr+bUWoTgRpdpwGSVj1uI4ZehSU5GflaTjW6a3Mi8g6PvyIvNQ1qaP+SCY1\/UWDUFEyu3Oo\/Bkj72+K+6kZGRnM+uzfvz8jKHQLY+N7rXbuWBEZPueMGTNg+vTpLOhItJCtCFX16I348vDuu+\/CzJkzbRUj7xXLnYlQ1dZ1zPQiQo2uqogo1PAnnNRxchvlqzZzeHtZuY69XpVbq2YWapApCr3Kx8cRofpFULPxRKjRVQgRhRr+hJM6TkSo57BCa9VsD5UIVW0tBdWL9lB9Iqnjw49cvupK1VF\/JJOa\/qLl8lWTzrpX0BYqXkmM8sWApKNHjxp7qn7lDWo8WahBIanJPGShRlcRRBRq+BNO6jjFooWqdnfee1mdc\/U+YzAjiVCDwVGbWYhQo6sKIgo1\/AkndZyIUM9hJacGDLrEm5pG7HsRoQaBokZzEKFGVxlEFGr4E07qOBGhqmGlQy8iVB20EKAMRKgBgulhKiIKNdAIJ3WciFDVsNKhFxGqDloIUAYi1ADB9DAVEYUaaISTOk5EqGpY6dCLCFUHLQQoAxFqgGB6mIqIQg00wkkdJyJUNax06EWEqoMWApSBCDVAMD1MRUShBhrhpI4TEaoaVjr0IkKNsBbOnj0Lc+fOhVmzZhWrOYjfLV++nNUjXLVqFaSmpkK\/fv2gT58+UL58eSVJiVCVYApbJyIKNWgJJ3WczhdC5WkGMdevWfUZETGz+qwYBYxVaKKZ6IEIVW1dB9bru+++g3vuuQcqVqxYjFCXLVsGgwcPhtatWwMW8l2\/fj3MmTMH2rdvD4899hhUqFDBUQ4iVEeIwtqBiEINXsJJHafzgVA5QVaqVAmWLl3Kqsi4IcZwJJJQ01BoLyJUL6h5HHPw4EF48MEHAUm1cuXKIYSK3yGZolWKCZg5eS5ZsgSGDRvGPrvxxhsdr0yE6ghRWDsQUajBSzip43Q+ECo+t0aOHAlYHg4tTSRTTISv2ohQVZHy30+L1INnzpyBadOmwcqVK+HKK69keSl5WSO8RVxQ6NrFqgodO3Y07vrYsWMwYsQISElJYeWLEhISbBEhQvW\/YPzMQEShhh7hpI5TNAjVrtoMJ6\/MzEzDNcvTAubk5LDnGjasNjN+\/Hi2vYX1S7HJFWrwM7k83IIFC2D+\/Pm2ZdREl2\/9+vVZbVTeVFzGaui770UWqnvMPI1Ady5Wlcf90Z9\/\/pn9VyRUXEBYqFf8DC+E+6poneKeqkqdvnAQ6pYtW5jruW\/fvtCgQQNP9x+OQTrKRTKpaZpwUscJLbZPPvnEtRtU7QrFe8n1UM1KrfHnDC9nZlUfFcu2cfetVck2udSaSuk1eQ9VVwtVx3XudV3wcVG3UH\/77TfmtsUN84EDB8Lbb79djFBxEWMlhVdeeQWqVasWcs\/43UcffQSvv\/46NGzYMOIWajhI2q9SuVWPb6Zu91uCuLbVHDpiRTKpaVxnnPAOIrHOrcjMDBtewBtds+iu5YXDUVa7ebCvaDggGYsWqVlBc1mDsUKoOq4ptV+Dda+oEiq6ep977jlYu3YtvPjii5CcnMzcvbKFaldd3qmyvXjrXIG4V4vBTUG0HTt2wKhRo9j+b1BzllS5dMSKZFJbbTrj5JZQT\/z0AhQe3wFx5WpD0oXD1QD4vUSa\/GzCwWZRtfwzdOemp6eHeNCsLEaZaK2sVtniFUkarV5ecByDN7OysgwCxz1YN8FMysAodpSfv05rqk6dOoD\/YqlFlVA\/\/\/xzePjhh1kl+JYtWzLcwkmo27dvZ+SHe7XUCAFCoOQgoGqh5m\/7F+R9N8q48aQLH1QmVXnvVEaPu3j557y\/\/LlVNRiZmDkBWWlJ3HPFvhikiRax7haq6qpDIwX\/xVILO6GaLQpcCPjmdP\/990Pv3r1ZwFGpUqU8Eaobly9eAEkV\/1EjBAiB2EcAX47Ru6VKqEimSKq8la17J5S\/9GklIGT3q90gbm2ixVirVq0QN66KhdqtWzcWV4LN7Nyp6LXDuBMeeISu5dGjR8OQIUPY8UIdLdSnn34aateu7Yg5WagmEFkRKgbwjBs3zhZU7ir59NNPAwlKctQgdSAECIGYQsDtPpxsoSKZIqmqNCci5C5Vvs+Zm5vLnlvZ2dkhxKiyh4oDMBJ46NChjBTlJgc+4d\/z5s1j5IvX79+\/v7aEqvryo6IT3fqE3UK1uuFt27bB1q1bi3391VdfwXvvvcc28i+44AK45JJL2NnUII7N6AY+yUMIEAL+EHBLqHg13EM9vW8llKnaWtndy6WUo3zNgoSsono5OXJCxTl5ABL\/jAcvOcWGyNfFvVorQjXb4\/WHurfRXnTl7UrRGxU1QrW6ZbOFxBM71KxZMyQrktvEDtGDma5MCBAC4UAgGg9pu3OoXB75TKm4n5qRkcGsT7Qi8bgfuoWx8b1Wq2AkGT+RuPE7K0LF73gBcjlAKhw6sZozGrqK5P3htWKCUFHQRYsWwfDhw6Fdu3bQvXt32Lhxo+vUg5EGl65HCBAC4UUgFh\/S4TgXKrp8k5KSwgu6x9ljUVdubzVmCNUsOT4GNGGCfMxxSY0QIATOPwRi8SFNhOouF3EsrWrtCDWWwCNZCQFCILoIEKGew58s1OiuQ351IlQ99EBSEAKEgAcEYpFQPdxmiRhyPuiKCLVELFW6CULg\/ETgfHhIlxTNng+6IkL9fbXynMIY+CSn58rPz2cJuKdPnw4\/\/fQTyxmM+7e9evUCOQBg7969rN8777zDzoNh1QksPXfVVVcZySvwkkEUTI\/GDy0onPCM3ssvv8yOSGG77rrr2GH05s2bh9xWLOGUl5cHc+fOZUkGMK0alhu87bbbWAITMQe1m3vCCPdXX32VzYvrCdfmQw89xI6T8WQosbiegsLqfHhIR+N3Ho5rng+6IkIFADwTiym7UOHyoWPMN4z5O\/Ffp06d2GFpni8Tc\/eKxc2RTDFVFgYdDBo0iGVI+ec\/\/8mKAr\/wwgshNVuDKJgejkVvN2dQOGFmF8TnyJEj7L+Y5OODDz6AxYsXszJ+GMnNW6zgJOoes9bg0YgffviBRaLXq1cPXnrpJSM7jOo9YXnCRx55BPBstrieVq9eHbM4oV6DxArzgOtWBCLSv8tYuR4RaqxoyqOc+JaMFhKmLtu\/fz+bRSZUrGF4zz33sGhiPJgdHx\/P+qGliufI0PrApBPY8FzY448\/zrKjYPUcbEePHmUWBda4xM8xIjmogukeb9v1sCBxQmsf8zcjrmjJN23alMnDCyUsX76c5XZOS0uLKZwWLlzI6luKusf7wgc+rpO7776bWeCHDh1iHgu0XvHcYYUKFdj9m52pNpsTSRYzjB0+fNgoKBFr6ylIrLAWMhGq6590VAYQoUYF9shdlB92bty4MWDuTKx8IxPqs88+y0hXrsWKBPC3v\/0Ndu3axcbhg9GqzJxczxUXVixlfgoSJ7RO8AXlpptuYoUKRLclkg8eeMfrYSH5WMGpsLCQWaD8ZQAf8rzhywjmZI2Li2Np4TCjjYruCwoK2Dh08\/L1xefEohLoCUE38KWXXhozOKH8QWN18803Mzx1q\/YUuadY7FyJV5eh1IOxozNXkqIbFx9+uM+F7jmzN1270nH43X\/+8x9GtnXr1mXWCZKnSL64X4bJoDEfMe8XVMF0Vzfro3OQOKGlbpWjlJ\/PQ7c6WnOxhpMZxOjWxnvB0oRIqO+\/\/75SXmok1AEDBkDbtm3ZdoTY0GWOGGIdYczzWhJwwvvzghW+1D7xxBNUQcrH7zuSQ3GbDJ+HsVaWTRUj2kP9HSkrdwRaqJiliZMhB5ZbqLjvxwl006ZNzAWMgTUTJkxgbj0kUnwgossYXX3oMg6qYLqqkoPs5xcnvH+0UDHbFRKN2LiFeueddzLMYhknfl\/oysX7RGscLSnVe8L1haTJq4aIOPEXj9tvv51tJ6jOicF0OjevWJUpU6ZYBSlMw7d+\/XpWeaVixYrGbePn69atY\/qoXLlyCBz4HVavQc8Abjno2tDLgdtU+OLZpUsXlnLQzz3hyxvGeuBePd477vn7ndMKu1isIONmHRChOhAq\/sjvu+8+5t7t0aOH4aLke6g4XLRI8XOMFMb\/8jZixAg2R9myZdlHQRVMd6PooPpaEaoqTpiPGQsf4J41Wr78wcX3ULE4Ms+DGss4Id48+Or\/tXf2utQEYRw\/lUajEHcgIZS4BR83oRWdmkLUEr1Ko1NJqBRqLoCO66B583uSOXneY2d3lj1n9+G\/ydvw2p35zez8Z56vJcIXk\/DS0lLx2PP3uZP8ZK3X6Jx+ymp5efnL9M4xic7qJ3OqihMWtJubm9Hx8bEFvqU4keiculrv2t5HgtogqPiwCKK5vb01P+ve3p75Ta+vry0VglNpEtTn52fbES8sLJgpjgWUE+rd3Z2Z7\/wJFWGiMLb3t9GUpq9MtB3grv9\/TlDbcOIeiCZ9hwsie39\/P\/5WbfriRuSXmohoNlKInw++Ku3TXxLUn7L6K4LaNSfEFJcVp3UyGHzGQuk8rWLf9ZoT6X4S1AZB5deIBSYQBJCgmpQL+Pr6aqcsBHVxcdECIz4+PsanEf427QAxAV9cXIx2dnZqTyltP5g+68lWF6lXwim9gJjjSCV6fHw00xsmcYSUHTJ8k8k3t\/EYMqeXlxebC5jS8Kuvrq6Oh6nJJ8\/G4urqyqKem3zN3uQbkRNQumBVZcr+zgn1t88pz4ngMHLlT05OLEDw9PT0v5ropfN06G6EWa+PEtQCQc0NCv5VcgoRWkyYOZ8XIuwDTCIHkXwn9N1zmjyRe7Zvb2\/mXz08PAwZbMPmifnAyZTcWqJzCVbzV+nY\/\/agpC5ZVc2pnCCU8q+bp7NcpKfBidQ1rCYUVuF9Ozo6+lKgJhqnWY5J3bMkqA2CioBwWqLgwPr6+phlyv3DgX92dmYfS0dQOWnxz1+TgholHaRq4uQEtZTT5+enFSuYn583U3rK6+VZKT8xYjoI7fcFG4g89dWREsvSsf+taTOJQ5es5ubmvkzVnKCW8q+6Zx+LdtecUqEaBBVTL2uWfwfbztOhcOpjbKqeKUFtENQUUUmpwYODAwtK8mZcgmgo4pAENuXZ4T\/1Jl+c\/pzUMGtG\/mB6TlBLOaU0ooeHh9Hl5aWd5LhSSUP8z+fn55bXG4lTChYhwpsNFmkyVVebPrHBwCTnq0fVFXbAF+39YFXFIoaw8EyD1WS\/coLahn\/frLrm5NctioOQJlglpvQ7Eqe+x8k\/X4LaIKjs6FjQMI+wm9va2ho9PT2Z3xSB9dWT0m6SkoNUxmGBo+wgicwkoPvFLuoH03OC2oZTSi+iatT+\/r6NACX6yEP0ATz8PAInTGiIKH52Nkxpk+BfNEy\/BLWxoy\/tE4sapwjSJOBEARLSG6pKD5bes+\/FZ1qsSgX1L8+p9\/d3S+tjXjEXU5Uuz46yn5ubm\/ajKHOq7zktQa0YgTrfIIsACxmLPRG+a2trlgazvb09ToVJt2RXSSAKpwPKDq6srJi4EkDiC+lH\/WB6V5wQVU4RBCXBZXd313ynmND9FYFTMulTHCR3IbQUdqCvbfo0WRyfjy3go93Y2Ggsjo9VBfcDG5ehXNNk5ftYF1TThn9f3KbBiQ1\/KpOa6xflMMlQ4IrAqa\/xyT1XJ9ShjYjaIwIiIAIiEJKABDXksKnRIiACIiACQyMgQR3aiKg9IiACIiACIQlIUEMOmxotAiIgAiIwNAIS1KGNiNojAiIgAiIQkoAENeSwqdEiIAIiIAJDIyBBHdqIqD0iIAIiIAIhCUhQQw6bGi0CIiACIjA0Av8AaJmF2fbiwOAAAAAASUVORK5CYII=","height":282,"width":468}}
%---
%[output:881c6099]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdQAAAEaCAYAAACoxaaoAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQl4VdW1\/yKEkBDGhECAgNCCCCoIikFE9KGtLQrFVmSyIKCgDwmIDOIMVpksLQIVWhCwLYN1eC0tPEvbVxSoFBQoylDzBAFJIIwyI4T\/f23evjn35Az77LvPOfdkr\/19fK25e1q\/vfb6nbWHtStdvnz5MlAiBAgBQoAQIAQIgYQQqESEmhB+VJgQIAQIAUKAEGAIEKGSIhAChAAhQAgQAgoQIEJVACJVQQgQAoQAIUAIEKGSDhAChAAhQAgQAgoQIEJVACJVQQiIIDB16lSYN2+ebdYaNWrADTfcAA8\/\/DDceuutkJKSwvIePXoUhgwZAlu3brUt26hRI7jllltgwIABcO2118Lp06dh9OjR8Je\/\/AUaNGgACxcuhKuvvjpW\/vDhwzB06FDYsmUL+9tdd90FM2bMgOrVq8fyrF27ltWHCfv09NNPi4hJeQgBbREgQtV26EnwoBFwI1Rjf0aMGAH4LzU1VYhQeVkk5Tlz5kDnzp3hpz\/9Kfv\/mGbNmgX33HNPrInNmzczskTixdSwYUNGui1atIjleeONN+AnP\/mJZfmgsaP2CIEoIECEGoVRoj5WCAS8ECoS45tvvglt27b1RKgIVMeOHRmRImmiZ4lp+PDh8OSTT1qSJf\/j\/PnzoWvXruw\/L1y4AM8\/\/zy89dZbkJ2dzcj2uuuuqxDjQEIQAn4hQITqF7JULyFgQsBIqEuWLGHEZ0zFxcXw8ssvw5\/+9Cf258mTJ0Pv3r3jCBUJdsGCBZCVlRUrWlpaCoWFhfDMM8\/Axx9\/DJmZmYyMcfl20KBBcODAAeadTpkyhf129uxZeOqpp2DFihVQtWpVuHTpEly8eBEeffRRGDt2LFSqVAmMS8I33ngj\/PKXv4Q6derQmBIChIADAkSopB6EQEAIuBEqdmP58uUwYcIE1iNcsr3vvvtcCZV331x\/69at4fHHHwfcC73mmmsYEeN+alFREduT3blzJ9x7771w4sQJ+PDDD9ky8ezZs6FmzZrw6aefMjI+cuQI\/PjHP4bnnnuOLT9TIgQIAXsEiFBJOwiBgBBwIlQMWIaeJHqR6KEa9zSNh5KsPFT0LvHAEi7R7tixI1a2efPmzMvFpVwkw9\/+9rfQoUMH+Oijj9j+KZbDvdXdu3ezA0nGNrEPuIdr9JQDgomaIQQiiwARamSHjjoeNQRE91CrVKnCiBC9U1x+FTnla8Sib9++jFxxOddIjK+88gr06dOHeaFGAj148GDsNC8\/vMQPNPHl43bt2kUNbuovIRA4AkSogUNODeqKgCih4hUWvKLStGlTBpUXQu3WrRu88MILkJOTw8qal27xKg3uk+J1Gn5V5tixY2x594svvmDLu8Y8119\/PfNweX26jh3JTQiIIECEKoIS5SEEFCAgSqjYlPH6iwihIjmOGTMGcJmX31\/FeoyHi3CPFJdxn3jiCba8jMSJe6zGQ0qYZ9y4cTBq1ChGsA888ABMmjQJ0tLSFCBAVRACFRsBItSKPb4kXRIh4HYo6dSpU2yJFpdmT548Cd\/97nfh1VdfZVdYeGAHvodau3Zttm+KZIf\/i3ukeEUGCTMjIyMmtfH6C+6Rdu\/enQWXwPx4EpifNOZ3TvGKzCOPPMLaxT3WZ599FgYPHpxEKFJXCIHkRYAINXnHhnpWwRBwI1QUFwMt4JUWJFZOnvh3M6HyazP\/\/ve\/4bHHHmMHizCNHDmS3Tk1nshdtmxZuShHxlO\/WM4Y6AH3TbEftH9awRSQxPEdASJU3yGmBgiBKwi4ESp6k6tXr2bkhx6qCKFivStXrmRLtOhR4lIxeqDGO67mqEhYxryUaw5FiHnMpEvjSAgQAs4IEKGShhACASHgZQ+Ve6Xjx4+Hr7\/+2tZDxXxIpLhMvGjRIiYJ7oPOnDkzFojBeO+Ui8pP\/PL\/xjpeeukl+PWvfx1DA08ZYz48LUyJECAE3BEgQnXHiHIQAkoQ8EKoeE3ltddeAwx673YPFTv31VdfsWD3eA8VE57kxchHeO3m\/PnzzOt977332G92S7nGKzaYj\/ZPlQw7VaIRAkSoGg02iRouAiKEisusePUFQw7yqyoihIqS\/e1vf2Mkit4mlsXISDz+rjHQPb5og6EE69atGweI8YqN+dBSuMhR64RANBAInVDx8AOeNsTYpviVjacMe\/bsyb62jXffMJLMunXrWGSXjRs3snx4+hAjvuAXNyVCgBAgBAgBQiBMBEIl1JKSEnYq8fPPP4d+\/fpB+\/btYdu2bbB48WJo0qRJbMkLAcJ4pHh6MT8\/n0V72b59O8vXpUsXmDhxYtw7jmECSm0TAoQAIUAI6IlAqISKr13g3g6GQrv99ttjI4BLT3hN4MEHH2QXz48fP87IFL1SDMnGH0Fes2YNFBQUsL\/hMhklQoAQIAQIAUIgLARCI1R8cgoPXeAyLh7zNz5Hxe\/iYcQXDBaOF9dxaXfu3Lmx9xoRMLwIj9FesCxFcwlLhahdQoAQIAQIAUQgNEJ1gh+vCaBniu8vIqH+\/ve\/Z14sPnLcokWLWFHcV0XvFPdUzW9E0vASAoQAIUAIEAJBIpCUhIpLuUioePQfPVM8Hbl+\/XrLIN34G15sxzt4zZo1CxI7aosQIAQIAUKAEIghkHSEykOp4QlfXBKuV68eI1R8w9HKC8UHmfHkr9l7pTEmBAgBQoAQIASCRCCpCHXfvn1sTxTv3b3++utw9dVXMyxUEur+\/fsB\/\/GUm5sL+I8SIUAIEAKEQPIggHehjTGpk6dn9j1JGkLFCC94hebSpUtsv7RVq1axXjsRqpclXyRSXEbesGFDrO6BAwfGHlf2Y8Dwkv2ZM2egWrVqkVMOFXjoLj9iqDsGJL\/eNkB2DuCLSniOJkopdELFg0V4xxQ9U3xQecaMGdC4ceM4DPG1DBWHknDZGO+7Tp8+nYV0w+S3h4pvTR48eJC1k56eHiXdUNJX3eVHEHXHgOTX2wbIzgHyUCVMsDFgw8svvxwXHYlXh0So4toMJ1SMymR8jUOi28JFzp07BxicvEGDBloSqu7yo6LojgHJr7cN0GkOhOqh8gNIrVu3ZvdI7dz7Y8eOscAOSErGqEheAzsQoQp\/ByjLqLsx1cmY2CmN7jqgu\/w6zYHQCBVfwEASXbp0KXTv3p0t95oTLv326NED0tLSYNWqVezNR3yaqlevXlBYWOg59CARqjKeFK6IjAl5qLrrgO7yE6EKm0v5jBjH9+GHH2axe+0SEi0GdsjIyACr4Pj9+\/dnAfJr1qwp1BEiVCGYlGYiY0KEqrsO6C4\/EapSk5o8lRGhBj8WZEyIUHXXAd3lJ0IN3u4G0iIRaiAwxzVCxoQIVXcd0F1+ItTg7W4gLRKhBgIzEaoJZt0NKslPp3x10YHQDiUFb9qBhS\/Ee6h0bSY49HWZSE6I6o5B1OU3R1fzOnswsMXhw4cBAxXoeBcd8eIYXHfdddC8eXOvEEYmPxGqz0MVdWOSKDy6y6\/TcpedrkRZB6yiqyU6J3Qun5+fzwLr5OXlVUgYiFB9HtYoGxMV0OguPxFqtPeQraKrqZgXOtaBIV9nzpwZ6Aph0DgTofqMuO6Eorv8RKgVg1CD3Cby2SSFVn0YW25BC0uE6jPiuhOK7vIToRKh+mxiIlM9EWpkhkqso2EMqO6Eorv8RKhEqGLWqeLnCsP+Bo0qeag+I647oeguPxEqEarPJiYy1ROhRmaoxDoaxoDqTii6y0+ESoTqZJ3wPed58+ZZZpk8eTL07t079hvmPXDgQCwcq5jVK5\/LqU3MzfeLvbRnrLNhw4awcOFCaNGiRVzjYdhfWYxky5GHKoucYDndCUV3+YlQiVDdCBWJZsGCBZCVlRXL+vnnn8OgQYNgxIgRMVL1QnAybZrLiLZnzrd8+XKYNWtWOVIlQhUkjahkC2NAdScU3eUnQiVClSU3JCoj2YoSnJs9Ntdrl1+kPSubig\/KP\/XUU4Ce6vjx42PVh2F\/3bBQ\/Tt5qKoRNdWnO6HoLj8RKhFqIoRqXOLlBFejRg22LIvJannVvKRrfLULy3glVKf20BtdtmxZOQ\/b6u9EqD6TTdDVhzGguhOK7vIToepDqFPf3w17j56DJlnpMP7uZkLmzY7ccMl33LhxMG3atNheJCdKI0Gal1fN9Vl5i14IFfd33dqz2te1WvYNw\/4KDYLCTOShKgTTqirdCUV3+YlQ9SDUpRuLYfjSHTETMP7upkKk6nRAyOx9WhGhkTALCgosl1qt9kbtDkIZA1i4tYfLuXbLwkieY8aMidtHJUL1mWyCrj6MAdWdUHSXnwhVD0JFMkVS5alvh1yY07eVq4mz8xaPHj0KQ4YMYeX5gSU78jIur65evRomTJhguRTMO+PFQ7XzPvky769+9SvLk8dEqK5DH\/0MRKjBjyERarQJRYXGRFkHRG2G2UNFMkVSdUtO5IZEieRovMaC9RkP+uB\/m\/cr+QlhJEOrfVYvhOrWnh2h0pKv28hXgN9FJ4dKUaNsTFTgoLv85KFG+4PCi83APdS1hcehc\/PaQsu9qBtO5GZuW8RDNV69wfq5p1tSUhJbfvVCqG4eKnrEdCipzFLSHqoK1nCoQ3dCkZX\/w+kH4cS+C1CrcRrcNra+z6Pkb\/WyGPjbq+Bqj7L8XghVBlE3D9V4n1NkT9OqD1bEbHX31VxWpD2rpV26NiOjCREs4\/fksIIkysZExRDLyP+v5cfgjwX7Ys3fNqZ+pElVBgMV2CdLHVGW32+b4XTKFwM74AlbvsTLDzANGzYs7m8rVqxg3ie+MYr3PzFNmTIFMjIy2P+3us8qSqh4eMmuPYyExMlz7969sb1eCuyQLDPP5374PTmIUMsjIGNMkUyRVHlq07sO3PtaY5+1w7\/qZTDwrzfB1xxl+f22GU6nfO1CDxrvhbZt2zbuDignOCRZnsx5vC75OrXH2zDKQaEHg59jobTo9+QgQlVDqGYPFckUSTWqKcqEogLzKMsfhs1QgXky1qEDlrSH6rPmRdmYqIBGVn7cQ\/1y\/Sm4qlP1SC\/3IoayGKjAPxnqiLL8OpBAUDqiA5ZEqD5rU5SNiQpodJefCDXaHxQ6kICKeS5Shw5YEqGKaEICeXQnFN3lJ0IlQk3AfFSookSoFWo4gb3c0K9fv9hF6SDE051QdJefCJUINQg7E4U2wrC\/QeNCHqrPiOtOKLrLT4RKhOqziYlM9USokRkqsY6GMaC6E4ru8hOhEqGKWaeKnysM+xs0quSh+oy47oSiu\/xEqESoTiZG5h6qMWiDjPlyahPrM744I1o\/3UO9ghQRqqjGSObTnVB0l58IlQjVjVCtohbxAPcjRoyA3r17syrsYvl6NU2igR1E6zX3iyIliSIX8XxhLDnoTii6y0+ESoQqQ6icQI1km4yEamVTKZZvxIlStPtEqKJIqctHhBptQlGhCVHWAb9thpO3aCZQ\/t\/GUIBWYf7MS7oYD9gptq\/TGJufgnvuuedgy5Yt7L1VjDFsfjqO12X1d7+xVKGridZBS76JIuhSPsrGRAU0ustPHmq0Pyi8kMCxtybCxZI9kJrTFOo88ILQ9HEKjj9u3DiYNm0aYBB67rFisHojQZqXV831WXmLoku+XHYeU9gYJ5gHzHd6Us74Ug723wuWQuAlYSYiVJ8HRXdC0V1+IlQ9CPXk3xdByZxBMWtSp9eLQqTqdEDI7H26PadWUFDAXpvh3qOdaXNqkx9I4uSJdRi9W+6x8ldw7AjV6lk3IlSfySbo6sMYUN0JRXf5iVD1IFQkUyRVnmrc8RDkDF\/oauLsvEX+MDhWsGDBAsCHw0UeGMcHvydMmMBIFZ90496tsSMiHipvv0+fPrFDUVgH\/3vHjh3Zki8RavwQk4fqqvKJZdCdUHSXnwhVD0I1e6hIpkiqbsntgXEkR+41Yl5M\/H1UXrd5v9K87yni6Zr7aXXKGPOYl5BpyZcIVeqeldvEsPtdd0LRXX4iVD0IFccZ91DPbv87ZLS+Q2i5F8s4Eap5RU3EQ0VP1pi4R1lSUhLzWFV6qHQoiQiVCFX260CiHBFqtAlFYsjLFYmyDvi9TeTmoRoP9rjtoZo9Vz4QVsRsdffVOHCie6hWe6V0bUbFrIlAHX5PDisIomxMVAyp7vKThxrtDwq\/bYbTKd9BgwaxE72cKPlhIn7Clnu4K1asYN5nXl4eO5SEyemajIiHinU4nfI1n\/zdu3dvbK+XAjuosJwRqMPvyUGEWh4BItRoE4qKaR1lHfDbZjiduOWkxcfA6h5q27ZtY0Rm3ONEkuXJnEeUULG83T3UTp06xR1WMsphdyDKbyxV6GqiddChpEQRdCkfZWOiAhrd5ScPNdofFDqQgJd5zvdk0WvGk75ekg5YJhWhHjp0CPAu1ahRo8oN1uXLl2HdunWAewobN26E7OxsGDx4MAwYMAAyMzOFxjWMAdWdUHSXnwiVCFXIOCVZJrs9VC\/erVmkMOxv0LAmDaHu27eP7RUg6FavHaxduxaGDx8O+fn5gHejtm\/fDosXL4YuXbrAxIkToXr16q7YhTGguhOK7vIToRKhuhqmJM3AvdGtW7fGemgOY+il62HYXy\/9U5E3dEI9ffo0vPfeezBz5kw4cuQIk8lMqMeOHWNkil4p7itw8lyzZg3zaPFv3bp1c8UjjAHVnVB0l58IlQjV1TBpkiEM+xs0tKETKt\/Mbt68OfTo0QNmzJhRjlBxIHBpd+7cudC1a9cYRqdOnYLRo0ezKCKTJk2CtLQ0R\/zCGFDdCUV3+YlQiVCDNurJ2l4Y9jdoLEInVNwTRULs2bMnbNu2Dfr161eOUJctWwazZ88uF0oL91XRO8U9VR6eywnAMAZUd0LRXX4iVCLUoI16srYXhv0NGovQCdUosB3g6MWuX78e5s+fDzk5OXEY4W8rV66ERYsWQbNmzchDDVqDXNojQo02oahQpyjrgA4koGKMRerQAcvIEKpdZA+7S8RWA8wHdOTIkexwE6bc3Fz2z6+ExqS4uBjq168PGRkZfjWTtPXqLj\/3UEkHojkHNm3axLabrA5Kqph0MvdQjUEbZPrg1CbWJyIrt6W8ffOdWSf7+8YbbwhduUlNTQX8F6WkJaEaB2jgwIFswviVzp8\/DxhHEz3rqlWr+tVM0taru\/w4MLpjEGX58YQrntMQIRmZSegWKWnEiBGxAAp2sXy9tpvI1Rdsyxxu0C6QvrlfnITxnAwGm3BLtWvXhjp16rhlS6rfI0+oMku+06dPh0aNGgXioeJ9roMHDzIvOD09PakGP4jO6C4\/Yqw7BlGWH0kA77sHTaioN2biSwZCNT\/fxm2IXZB8o43hhPrmm2\/CTTfd5Gp+yEN1hcg5g90aOx1KShDYEItHef9MFWy6YxBl+f3e93PyFs0EahV60CrMn3lJ13x3VNRDRZLEsys1atRgHxTYFt6meP755+HVV1+NW7blXqr571aE6tfHiar5mkg9kfBQ6dpMIkMcbtkoG1NVyOmOQZTl90KoZ3f9HErPfAUp1RpBRstRQurjtOQ7btw4mDZtWuyRcE6URoI0nyEx12f18osXQsX3WI3B+M3LvVxIkWVfL1gKgZeEmSJBqDywQ4MGDeKiIlFghyTUKFOXomxMVaGrOwZRll+UBM7vextObx4bU5mMliOFSNXpgJDIw+BGwsQgN\/jaDJaze8qNLyXPmzfPUr2N3qPVgU87QuVLwRjFrnfv3pZ1i2Kpat6FUU8kCBWBWbVqFYvx27lzZ+jVqxcUFhZS6MEwNMZjm1E2ph5Ftc2uOwZRll+UBJBMkVR5qtr4fshsN91Vhey8RU5QWAG\/Y2+3h2rcv1y9ejWgV2n34gsnVLtbE8YOW+2LEqE6D2lkCNUqOH7\/\/v3ZgYGaNWu6Ki5mEJ0cQpUJZoqyMREU0TGb7vIjOLpjEGX5RW2G2UNFMkVSdUtOy69IaEiO3GvEvJjM3qeZ+MxProl4ulb99EKotOR7BcGkIlQ35Uv0d9HJkWg7xvJRNiYqcNBdfiLUaH9QeLEZuIf6zeENUKVuvtByr5u3aG5bxEPFqHPGxD1dvLqHj5C3aNGi3Olhu3luRah2h4\/oUBIRqgq+cK1Dd0LRXX4iVH0I1dUYWGRw81AxNKsTEVodOjI3Y0XMsku+dnulXq7N0ClfGU1JwjJevjZVdV93QtFdfiJUIlQnW+IW2AFP9PIlXn6AyXjqFv+2YsUKRrp5eXnsUBImYzQlq\/ussoSKdZuXokWWe7FcGPZXlR0XrYeWfEWRksynO6HoLj8RKhGqG6Hanbg1h\/OzuoeKEYeMD4NwjxVJlidzHi\/XZjAGgNXDI4mEHiQPVZJMkq1YGF9IuhOK7vIToRKhJpsdDKs\/YdjfoGUlD9VnxHUnFN3lJ0IlQvXZxESmeiLUyAyVWEfDGFDdCUV3+YlQiVDFrFPFzxWG\/Q0aVfJQfUZcd0LRXX4iVCJUn01MZKonQo3MUIl1NIwB1Z1QdJefCJUIVcw6VfxcYdjfoFElD9VnxHUnFN3lJ0IlQvXZxESmeiLUyAyVWEfDGFDdCUV3+YlQiVDFrFPFzxWG\/Q0aVfJQfUZcd0LRXX4iVCJUn01MZKonQo3MUIl1NIwB1Z1QdJefCJUI1ck6OT3fZhfYwRgFSczyxedyahNzGoPxHzhwIC7qklt7\/FFyqz6GYX\/d+qv6d\/JQVSNqqk93QtFdfiJUIlQ3QrUKA2gVzs8uOL5XEyYaKclre5wwjQ+gG\/tGhOp1pJI8fxgDqjuh6C4\/ESoRqgyhYhmrGLxePUartv0gVKPXS4Sa5ESoqntEqKqQFK+HCDXahCI+0vY5o6wDXmzGh9MPwol9F6BW4zS4bWx9IeicyM3sIVrF8rV6SNy8pGsmOK+EWqNGDbYMjMmuPR6g\/9133wU70veCpRB4SZiJlnx9HpQoGxMV0OguP3mo0f6gECWBfy0\/Bn8s2BebMreNqS9Eqk6vzYwbNw6mTZvG3jDlHisG0jcSJO5ZOj3xZvW8mxdCdWvPbCOclolFsVRhd8KqgwjVZ+R1JxRd5J+7+xwcOFcKDdNT4NFm6XFapQsGdlMpivLz8by4axMsG\/1Q7KCOnYxIpkiqPLXpXQfufa2xq3VxOiBk9gatiNBImAUFBez5NizHn3yzW\/K1e+HG+BKMW3tWbRChXr582XXUK0iGML6QomhMVA63DvL\/ofgCPL\/jTAy2R5umx5GqDhg46UzU5DeOJxLqmWlDXQnV7KEimSKpuiU7b5E\/5I3l+fNpdmRlfNx79erVMGHCBMulWd4XLx6q1fKt02PiRKhEqG46n9DvUTMmCQlrUVgH+ZFM0Qjz1CM3DSa1qhb7bx0wqEiEahxPUUJF+XEP9cv1p+CqTtWFlnv5Mq7dY9\/mh7yRrDCZPUMzwfETwkiGmEQ8XTtPVqQ9Y1kiVCJU1RwSVx8Z03NQVFQEDRo0gPT0+KVQX4EPsHKzh4pkiqTKE+lAtHRAxkOVVTcnb9G8oibioWZlZcV1hXu6JSUlsHDhQrYfSx6q7Gi5l6M9VHeMEspBxjRaxlR2sHHPbdPxi3BT7VTaQzWBGMU5wMcTPdQPnh7suuQrqzdO5OZ24AjbtDp0ZO6LFTHbecVmb9Ocz6098lDJQ5WdC0LlomhMhAQTzKS7\/AiT7hhEWX6\/z104nfIdNGgQO9HLl3j5AaZhw4bF\/Y1fWcnLy2OHkjAZIxVZ3WcVJVQ8vGTXHj99TEu+ZQiQhypIDG7Z7E55RtmYuMks8rvu8hOhRvuDIghCtTtxaxd60HgvtG3btrFDS0aPFUmWJ3Mer0u+Tu2ZbQB5qOShivCCYx6nU566E4ru8hOhEqEmbGAqSAV+f5wkA0zkoSoYBadTnroTiu7yE6ESoSowMRWiCiLUCjGMZUL4NaBOpzx1JxTd5SdCJUKtYGZUWhy\/7K90h3woSB6qIlDtTnnaEYpTZB1FXUqKamQJVSYuarn9nPd3w96j56BJVjqMv7uZJR5nd\/0cSs98BSnVGkFGy1G+YCaLgS+dCaHSKMuvAwkEpRI6YEmE6rM2WRkTt8g6Pncp0OpljKlsXFSjYEs3FsPwpTtifxp\/d9NypHp+39twevPYWJ6MliN9IVUZDAIdJJ8bi7L8OpCAz8Mfq14HLIlQfdYmK2PiFlnH5y4FWr2MMZWNi2oUDMkUSZWnvh1yYU7fVnGyI5kiqfJUtfH9kNluunJ8ZDBQ3okQK4yy\/DqQQFCqoQOWRKg+a5OIh2qOrONzlwKtXsaYoof6s\/+eBWezjkDG0Wx4uGoNuHXuE+X67bRszj3UWlfnQmq1qpDfsBr85q7cuDrMHiqSKZKq6iSDgeo+hFlfFOX3Ghw\/THyj0jYRalRGSrCfYQyo0x6qXWQdQXGSJpvTPqSMMV25ZwW8svHFmHz3fVYEg1sPgzoPvBD7m8iyed8\/7IUdNWrGypiD1p\/8+yI4+eETkJqbDheLz7HlXmMbqgCWwUBV28lQT9TkDzr0oNd7qMagDTLj6\/TCDdZnfHHGa\/0Y3Wn9+vVxgSV4HWHYX6\/9TzQ\/eaiJIuhSPmrGxCscbvuQMvK\/vPFFWLWn7GL6NYdOwmS4FXKGL4x1T2TZ3C1PyZxBgKTKU407HoprwysWdvllMFDVdjLUEzX5ZYPjy2DtFilpxIgR0Lt3b1a1U9AEL22LBnbwUifm5YRpftCcCNUrkhHJH8YXUtSMidehdNuHlJHf7KFinwbk3A1D73jF1kO1WjZ3C1qPZIqkyhMSNpKq6iSDgV0fpgqcXFbd\/0TrUyl\/on0RKR+0h2oXBtAqZKDVc2oiMhnz+EGoRq+XCNXriEQ0PxGq+oFz24eUNaYj\/j4UNpd8HOvw95t2h2c6lC0D4w9OAenxd7x6M3fPOTjUuir8qFvtckHrMc+xtybC2e1\/h4zWd\/iy3IttyGJgHi2Rk8vqRzjxGlXJn3hPxGuQCY7\/xmfzoOhMETSo1gAGXztMqDEncjPyjuiaAAAgAElEQVR7pPy\/jaEAzU+zYaPmJV0zwXkhVPNTcM899xxs2bIl7hFzrI\/HE3733XfBjvTDsL9Cg6AwEy35KgTTqqooGhOvkMz84Ddw4FwpNExPgZFdHowrLiu\/2Ut9usOL0K1pd+Guqbh6I9yYS0ZZDMzVipxcVtVnlfWokl9ln9zq4ucCPi4EeGj8O677imZ9Hdx6qBCpOi35jhs3DqZNm8aeXDMSpZEg3V6ksXodRpRQOQHymMK8LiRPY8B8s\/dLhOqmXRXk9zC+kKJoTLwMt9vhIFn50XN8Y\/s82JFTHTrd8SRcqDYwRtqPNnN\/V1XF1RsvODjllcXAykN9+uNj7NTyxTPn4ZUb6wBeB0r2pEr+oOQ0rrogoQ6bXcmVUM37\/lYrKlb9dzogJPIwuJEwCwoK2GszWM78CLmZ9OwOQvEDSbxeLGc8BMU9VuMrOESoZQiQh+rzLI2aMfEKh9vBHxn5zXub\/\/jhPJie+0Csa+bTulZ9Nnuo977WGNr0ruNVPCX5ZTCwavjX\/zwFPz190RMOSgT4v0pkljSxqCr5VcriVJfxXIAoocquqNh5i\/xhcOznggULAB8OF3lgfPXq1TBhwgRGqvxBcbOsIh4qb79Pnz6xQ1FYD\/97x44dLUmbXpuh12YSmqduIQSjZky8huJzO\/hjJ79TO\/z07fIO4+FQjSbwWaPO7H956pGbBngIyS395wNb4fLOi1DpmlT4xVtt3bL79rsKHcAPBPx42XN7pmccVAgmu6QZRUJFD\/UXnxVCcUpDuLTrY3j3ZytcPVSUEz84cN+\/Xc6NQsu9WMbtgXEkR+41Yl5MZu8Tl32XLVsWI17zvqeIp2vWEV6H8ZQx5qEHxp1nE3moJnxwqfFiyR5IzWnqekjFbbkzasbE7QqMlSq5YWBFJm7toIc6+587YXmHK48lm5NIIAzRAzx2H0QqT9OqIFRcwkas\/\/lYdgwOERxUySG7pBm1OYD9TZZTvuYtKhEPFT1ZY+IeZUlJScxjJQ9VxSemdR1EqAZczEuNdXq96Eiq5uXOtl9cgIePXobbxtaP1arCmMoOv1dv0+0KjFU\/ZJZ8530wGA4c\/ifUg1Lom3qWRScyhvwzkzS2iwee8N9NtVMtT+ua+2Y+wHPrt2vDiuHt4rLZfQyIkrGxMqeVChU6wJewP7u\/luOpZWOfZOSw0zXZJU07QlXx+IHsvHArlwz3ULGPbgeORDxGzGNFzHZXdTg2tIfqpiURJ1RcmV63bh3MmjULNm7cCNnZ2TB48GAYMGAAZGaWLYM5weB2KMnrRX+zQc4suQhN15wG3OPjpKrCmBplEjVEbl6gFU5uV2Csykz81z\/gvSNlMXLvy94BL7S5xfaDwmyY+1Y+Cw93mBQX8s9M0liZiDdmRSY89CAe4sHDTMZXZ+w+BrrP2Qzr\/vd4rDqrOMDGtmS8dK\/TFcnxtxP2QKdtV5a6rxpSC\/q\/cpVjNapPBcssaVoRqnF\/+\/Nu\/wU53c\/BNS2vEl4m9Yodzpn1+RlwOieVfZS5HWr7ry++gBe\/vOLpXdy1Cc5MGyq05Ou1X5jfLbCD8fAPP8BkPGFrvLKSl5fHDiVhMh4ksrrP6kaoRiK2OuXL\/2aWmfZQI7KHunbtWhg+fDjk5+cDbpRv374dFi9eDF26dIGJEydC9erVXfXZjVBlLvqjZ4IGFa+N8NS3cgqM73Il5J0ToXo95OHlKoiMt4n9Ra8WDdDJI\/WhXqe+cd62FcBzlyyDN2odhEtVboDK32yBwSfqw6P9+tgSKl86vFBtAJSm5EK7DIBf3dGb3SnlV28qlRbD63trxzWHHymYeB43o4h5b128C043KVst6FYL4JX2ZfUu\/+ufYXLKzbF2JpT+E0prtol7pQZ\/dCNUq5WKXlvPs3prNU6DDiNqQVFRETRo0ADS06\/I4XUpFskRx\/\/7\/yjT89vG1HccH5F4xq6TxkMGPoa5lw\/AoMp\/jL3cY54DuHy99FIpI7j9+WtZCymlxYDj82LnkUItGvXFSS8Qs9nvHI1bKnc71IZz58n3HoD\/rd+Y7aHu+eMIXwnVa+hB4z3Utm3bxvZO2fw9e5aRKl5t4cmcR2TJl5e1u4faqVOnuMNKPD8RagQI9dixY4xM0SvFLyNOnmvWrAE8Ko5\/69atm+tEREIdtGwN3HLvfXDzt5vEfanyCcpD0eGX7NhBj7rWiRkeWn8UtpxPieXFZcn57a4YPjtC\/eP212HS7guMVNCY4CR3uwxuNEToDSNx4+lVq2Q8WIG\/52U1gRE33eEqD5L8n9\/YzvJhYHqMoWtcwjZX8MmszfBEre3MOGL6ztE7YXL\/urFsZvl\/se5nsOh4CrsGwxPihXGNeUo7s5gR9KUqZQeJcDyMHy1Yxmn5F8lk4pcXIT27jIDMh5n4Xu2nDTvDdQfWwuM3XwNPFeXDX4rehCrpR+Gbc1lQsv+7ngkVxyZn+xVCxf\/\/UC5A8wdLY4TKia6g6e+hUfph6HL9ddDyFuv9Yo7Jrn9MgRc\/vAPOnmkMKZeKIe3Mm9Dymqtg4iOjHcf0+3MKoah12Xi4kQmvzGoJ22l1xOjVYR3fP\/8raFi3A5Rk3gz4gVTz+O9g+I0j2QeF+bQyb7PKZzvh\/hb1ISfvz\/AlfAcup+RaepRWWwK8DrN8OGdmtqwCJa3LrlqhLqEu8I8ys6y\/2F4If34rg1V5adcm+GzDE74RquuETLIMfE8WD0bhSV8vyezQ2G2TuB309NJm0HkjsYeKA4FLu3PnzoWuXbvGMDp16hSMHj2aHSmfNGkSpKWlOeJntzRnN0HdjDY2hgT08y9ugZTqZadQ8e8\/qHcUJl77LUtCRbJ77NMv4ZP0IbH+3pOyCp6u+YntQ9eoZO9vOA6765cZhhYnS+HO66tZLmFhfowSZExoRDgxWS198eVY7j1i2ev3t4Y2jdfAmezDlhFgzF4eluFGjU+MjHMn4IlWdZgxfefTD+ClkjaO45R67n32+8X0u13ng90YjVt7CP78Tbw+dKh5GH51Y\/O4Os2Rkp7+5Dj8+eCV9jHVPJEBE2p8At+5Jz5KE\/\/dybjzPPd9cwEGtDgdI1T0Ns\/vewemXbOAZVmcNgyO1Psh5NWJ\/8jj5fnH0eKqZR94lb\/ZylYE+lauDE9873FLnLYs3AZjz+TEkcnt3wDM\/G68928ubNYd1Jv7\/3UekJx4MnrH2D\/00t+v0iP2e6VLxXC5ctkdWRzT\/2xynH00Wi3pY8H6276EnErZsLP5urix5x9CONeQaDd+XTfu48qs48YT4HbkjWVQd3Ae4BjyZP6483vJ11XBQ8pgt4fqxbs1d91IqIeatmd6wJNxBco4HqIfgCHBVK7ZSBAqHgmfPXt2uXtV6Fyjd4p7qvyulhOw5onMDazdBOd1Yb5b0v9azoNEAnrxsy\/ivC1eBg1e92pb4KmbhpVb7rMildzSA9D20iZW\/Kq8TvCfrcsMv5vRNntebvmNxsPo5Q375zr4+Mj+ckSGxjD91DRWzBwB5uHNp+K8S8zDDw+ZDRW2hV6m8e9W42X2RkUmi3niWY0pGvmb6qRCg2oNy1XJDj1lpMRNcp4Jx6VxVhNLwrOS31x548rb4I7K65mHtuR\/F8M\/9++Gb458BPWy7oHSlPpxRGT8QOAfJJh31TftLWFAEmqf3QiqtK8W580hyT37Xgasvv7WuHLoMWMbrVqkl\/sQ4+1ZjU\/O9nOQWXKJedzXvn2C3end\/0QOG0\/sH\/7v1so3OQ4V\/8i0w8xMwsbKsM9biv\/k+pHFdY+XxZUP48qGiC4Z8+hKqIgB90a3bt0ag8QuTq8IrpxQF7yQD79u\/gJ8cqHsg8tuzotekRNpP4g8kSBU\/CrCJ4Hmz58POTk5cbjgbytXroRFixZBs2bNHDF75ONC9nXLEzewuLRkXHK0qgRJ5cbsPKhfo+ykKN4523fpets2S0\/thfbZeZBV5SLzzlJTryyLikxy\/vXsV36zsTJ\/qZuFwg8ExKtBZkN2z44nN3I01+NkNBNReNWG1K4vxnHBPCLyo8yIX7t6N8LmQ2XxiZ08cLOnJIoN758TCfO6ZGTBskiuiPfWbzmvCJn7bPWhJSpXWPl0JlTVmHNC7TR6GHwqGOvYuH2muj9+1BcZQrU7lWY+Wm4HkqjX5gfIVCchQAhEEwEiVHXjxgm12rhfQmpL59UMY6tRWvbVhlDdlnXVqQ3VRAgQAhUFASJUdSMpS6hRWvaNPKGKLvkSoaqbGFQTIaALAkSo6kZallC93kFX12PvNUWCUFUcSrJa8j135BRcn1YNju+7csrvTE5q7PoH\/nfd03vYtQZMN1zcBFWb3A+fnL+ykY4HHUQOOzQ7dgGaNQb4OrUs9qxdWdxn48m4T+l9WK+UMLeDB0qqlVyEqzq539nlbZrrwD0NnnBPzGo\/mEc1wnwXL15k\/\/j9S1HcjDIjLhlH6sLlyvXjxseYx9gv\/PuJfRfg8xplV5nMGLJ9wJafwoFd18V+KmqaAqnVrPcF7cbGaSyzPm8Zq7t++1Q4dOZKv7ANfs0IM2BfzOls1hE4m31F9zBhXea\/Gcs46YtVH814GetyGnPMh\/dMiyuVHeyyw6B9WjGkVMsrpwNu+owY8TmJeY366vSbVb3GsxFsH9uw\/2\/uN+otnze8Lr8J1em1GXPwBKc7nm6YGn93ahPz8djBInXiPVXzM3N25TihdnnlDcslXyu9E42MJtLXIPJEglBVXZt5dm0h\/Nene6BjfkdgZHrpLIucg\/frvlx\/ik1cDO2GkxAH8sS\/i6D9xcWQX3sXVKmbH7uozgeGP0JcemZ\/3Ik1\/juehhzZoQrk3H4+7lK\/mdzxFOmFT6rDuUqLYE\/WURjR8U7XO6kiyiESGlGkHreHvM3tGL8ozfdQ8aL91A++ZuHzMBnvB\/K+NC39NfsY4MEi8M5l3kedoc1vhrAL+sYA8XbXZma9cxQWZJURatMvvoAWJy\/DsdtbMgL7zrZ1kFntE\/jrnB\/FIFh\/\/RkofDQv7u4qHiiqcv59uOP8QTibddgy8PmE3x6GD1MOwrnan7KrInidpc2eOdBj8Z1w6HgHqNe+Mlw\/uQvTgZlriqD0vQvw9W3ZTHbsC+oJptqN0xiRcAL5S+3fwtEWOyE79Vqo8svroNGeg4CRhY602An1jz7EyuDvIsHY+RhiGREj5TbmRr3hY4rjYvxQGFr\/MDutLhMtzDgnzfegnX4z9ss8z\/IzP4N5N5ederaLlIX6+ZdtV650BBHYwep8iFVwepWEKhIpyc02WMUJdirjFlgHy3rRO7f+hfF7JAiVB3ZAg2SMiiQT2KFfv36evsBEBoXfE9ySemWjvUp2R2Youx24ZBklB\/PM2vR3Rtzo+Q68MA82XrgTZh8cDp2b144LjyfSvl0et5dgEqnbXNZuIlgZU6sPGH5snp08PfoS\/KS4XizoBRJqm18PgbwNnVmzeF3jwt01HInBbCz\/Y+cSFr6wxh1XiAiT1ZupHzxQNe7+asNP\/gB4bcopgALWg1eo\/vXjK\/dKMdU9fQFalZxk\/zug4\/NwquXdcR9VMthjFKtvDm+w\/LiTqU9lGRzT+VmV4k7+8r0vGUJV0Te3ONNO84PraFHlT+H1T\/wL7OB0r9MqZKDd491e8ErkLilvh5Mj\/rfTU3HGfokQqhc5kjFvJAgVgVu1ahWMGjUKOnfuDL169YLCwkLloQcTGSA7Y2dnTGTi5sr0L+wvPhljao66k7t1A3SZ3CAmvsjbpmZj+Xz2Ifhhm6vjILR7MxU9TvzYqbf9fOzOpV1EKqyQ13PFe9wFaTd8DcWnD8Ta6vetAfD9rB4JE6rM+AdZxo6gZHRARb9FPijd5ocXEpCJ8ONGqEYC5R6qMfSgFZmZl3TNd0dFCRVvUOB1Rd4eb+vIkSOAjgkuSV911VUwZswY27dXiVBVaLIPdVgFx+\/fvz8LkF+z5pW4uW7Jy+Rwq8v8u91kcjImyexxeJXfLr+MMTV7FvjgwHf+eBJqNa7ClkOdQiEa++FmLDHvqOenQPHpIsjNbAA\/n3Ql\/J\/M4+RGr\/uDO16HVf8\/6AdP\/5H7HRicN6zCEyrKa4W5jA6o0j8RHXBqS9RmuD2QYNeGU3B8894kJ0ojQbq9SGP1fqkXQsX3WI3B+M1yID5EqGWoRMZDVTHBRCeH17acJlOYxsSrHH7kl5HfilAfPlpqG7dYtt9Oj2aL7tNZtW2ud\/T1T0G7qjdpQahWeMjogOyYqi4najPclpedCNUuOL7Iw+BGwsS45hgYH8uZHyE3tu90KMl4IEnkjj8RavzIEqEqmIFOkynKxkQBNFIHUswHim5+\/Qi0a1sJHnvmWyq6FKsjkUez3TpifOqs37cHlgs\/6Va+Iv0e5TkgSqgiy8tWY2rnLfIDP1iGh1UVeWB89erVgF6l076mFw8Vb1g4hXUlQiVCVX4oyWkyRdmYqDDqMvIbTwLjPuaJXUXs\/c85fcveXVXRt0QezfbSvgwGXupP9rxRll+UUHEMZJaXncgNPUQkR+41Yl5MZu8T8xmJz\/zkmoina6VD5nqt8hChEqEqJ1SnyRRlY6LCUMvK\/5eCObB\/Sx0obVIEY+q1Z2SK75KqTq\/\/93dh6+liaJuZC49978+qq2f1yWKgojO4T1965ivbV4xUtOFWR5jyu\/XN7XcvhOpWlxcPFfOa2xbxUPHlLWOyutpCHqrMSImVoSVfMZykc0XZmEgLbSgoI7\/5BPTHlQfaPp+WSB\/N7WS0HFnurnEi9fOyMhioaDco+dz6Gpb8bv0S+T1MQnU7cIT9tzp0ZJbLiphF7qGShyqiIfF5iFC9Y+apRJSNiSdBbTLLyH9681hAMuCpauP7IbPddBXdiasjqHZkMFAhbFDyufU1LPnd+iXye1iEypdt8UQvX+Llh4mMp27xbytWrGDXVvLy8tihJExTpkyBjIwrj6Rb3WclQhUZfe95iFC9Y+apRJSNiSdBFRLq3\/5nPtxw8uVYjVtqPAN76l2CojNFlo+cy\/YzqLvAYelAUPK54R+W\/G79Evk9CEK1O+VrF3rQeA+1bdu2cYeGuMeKJMuTOQ8t+YqMvFweIlQ53IRLRdmYCAvpkFFG\/qnv74azu2ZCfu2dsOH4NfA\/Vb8Nl7IWxVoxP3KeSD+DuAssg0EiMhnLBiGfW1\/DlN+tb26\/+02obu1XpN91wJII1WeNjbIxUQGNjPzDl+6ApRuLY803bL4Maudsiv3395t2h2c6vKiie4HUIYNBIB0LqJEoy68DCQSkBuUOWQXVbpDtEKH6jHaUjYkKaGTkRzJFUuXpnlsKYTfMjf330x1ehG5Nu6voXiB1yGAQSMcCaiTK8hOhqlMSHbAkQlWnL5Y1RdmYqIBGVn5c9l1beDz2WIAxUMLga4ep6FpgdchiEFgHfW4oyvLrQAI+D3+seh2wJEL1WZuibExUQCMrP4b+w7cvazVOE47dq6K\/ftQhi4EffQmjzijLz0lg5MiRkJ+fHwZ8FabNr776CsaOHetLHIBkAYkI1eeRiLIxUQGNjPzm4PS3jakfaVKVwUAF9slSR5Tl379\/PyOBDRs2JAucke4HfpRMnz6dXfGpiIkI1edRjbIxUQGNjPxW75Q6PZ+mop9+1iGDgZ\/9CbruqMuPpIr\/ZBPKf\/z4cahbty6kpqaWqwafCsR\/+BYw\/quIiWNwww03QNOmTSuiiEwmIlSfhzbqxiRReGTkl3k+TbafMm9Yem1LBgOvbSRzfpL\/nNaPI6Bu6qIDRKg+WyJdFMkORiv5RUgskefTRIdU9g1L0fp5PtIBvQlF9\/EnQvVqMSKSP4xTZrpPJrP8QZGYiErKvmEpUrcxD+kAEWpRUZG27+ESoXq1GBHJT4Qa\/ECZySQoEhORVPYNS5G6iVDLEKAPCr0\/KIhQvVqMiOQnQg1+oNw81EmtqkGP3LTgO\/Z\/Lcq8Yem1s0QoehOK7uNPhOrVYkQkPxFq8ANlt4fKTzU+2iw9+E4F3KLuBpXk1\/uDggg1YIMTVHNEqEEhTct9tORLOsAR0P2Dggg1eLsbSItEqIHAHNcIGRN9rgzYaZfuOqC7\/ESowdvdQFokQg0EZiJUE8y6G1SSn5Z8ddEBuofqM8fookjkndgrEumA3oSi+\/iTh+ozyYRVPXmowSNPxoSWfHXXAd3lJ0IN3u4G0iIRaiAwK1nyFYmmFLw0ci3qblBJfr09dCJUObuR9KWIUIMfIhljmkzRlFQgJoOBinaTpQ6SnwhVFx2gPVSfrY4uiqRyDzWZoimpUA\/SAb0JRffxJw9VhRVJwjrIQw1+UGSMSVAhAYNCQwaDoPoWRDskv94fFESoQcyyENogQg0edFljGkRIwKDQkMUgqP753Q7JT4Sqiw7Qkq\/P1kQXRVK55OvzkARePemA3oSi+\/iThxq4yQmmQfJQg8HZ2IqsMTm76+dQeuYrSKnWCDJajgq+4wpblMVAYRdCrYrk1\/uDggg11OnnX+NEqP5hq9JDPb\/vbTi9eWysyoyWIyNNqkQoehOK7uNPhBq83Q2kRSLUQGCOa0TGmCCZIqnyVLXx\/ZDZbnrwnVfUogwGippOimpIfr0\/KIhQk2Iaqu8EEap6TN1qlDGmZg8VyRRJNapJBoOoymrVb5KfCFUXHaBDST5bLl0USeWS79KNxbDrH1Mgv\/ZO2HD8GsAl3\/F3N\/N5pPyrnnRAb0LRffzJQ\/XPtoRaM3mowcMvY0yGL90BSKo89e2QC3P6tgq+84palMFAUdNJUQ3Jr\/cHBRFqUkxD9Z0gQlWPqVuNMsYUyRRJlSckUyTVqCYZDKIqKy35lkdA9\/EnQq1Is9kgCxFq8AMra0ymvr8b1hYeh87Na0d6uVcnY6Jy2T94TfWvRdk54F+Pgq9ZFwxoD9Vn3dJFkciY2isS6YDeS566j79OH5VEqESoviJAxoTeQ9VdB3SXnwjVVxMbXuW05Bs89mRMiFB11wHd5SdCDd7ushYPHToEBQUFMGrUKOjYsWNcLy5fvgzr1q2DWbNmwcaNGyE7OxsGDx4MAwYMgMzMTKEeE6EKwaQ0k6wx+XD6QTix7wLUapwGt42tr7RPQVcmi0HQ\/fSrPZJf7yVvIlS\/ZpZDvfv27YPx48cDkt6SJUvKEeratWth+PDhkJ+fD3369IHt27fD4sWLoUuXLjBx4kSoXr26a6+JUF0hUp5Bxpj+a\/kx+GPBvlhfbhtTP9KkKoOB8oEIsUKSnwhVFx0IfQ\/19OnT8N5778HMmTPhyJEjbNqbCfXYsWOMTNErnTx5cow816xZwzxa\/Fu3bt1cTQYRqitEyjPITCQkUyRVntr0rgP3vtZYed+CqlAGg6D6FkQ7JD8Rqi46EDqhTp06FebNmwfNmzeHHj16wIwZM8oRKhIhLu3OnTsXunbtGrMBp06dgtGjR0NWVhZMmjQJ0tLSHO0DEWoQ5jO+DZmJZPZQkUyRVKOaZDCIqqxW\/Sb5iVB10YHQCRX3RJEQe\/bsCdu2bYN+\/fqVI9Rly5bB7NmzYeHChdCiRYvYnMV9VfROcU91wYIFrB6nRIQavJmWnUi4h\/rl+lNwVafqkV7u1Wn\/yE67ZHUgeG31p0Xd5ddpDoROqEYVtiM89GLXr18P8+fPh5ycnDitx99WrlwJixYtgmbNnOO9EqH6YzCcaiVjQqd8ddcB3eUnQg3e7rIWnQgVf7PyQpcvX85O\/pq9VysRiFCDH1hZY3LsrYlwsWQPpOY0hToPvBB8xxW2KIuBwi6EWhXJT0u+uuhAZDxUlYQ6cuRIdloYU25uLvvnV0JFKi4uhvr160NGRoZfzSRtvTLyn137Gzg275GYTDXuewZq\/PDZpJXRrWMyGLjVGaXfSX69bQD3UL3awdTUVMB\/UUq+Eyr3Co2gDBs2jF2RMScZD1VmydfY7sCBA9mBJ7\/S+fPnoaSkhC1VV61a1a9mkrZeGfkvvzUWYNM7ZTLd9COo9EB0HxiXwSBpB1SiYyS\/3jYAVUZGB2rXrg116kTrMGIkCFX1oaTp06dDo0aNAvFQz549CwcPHmRecHp6uoQ5inYRGfnNHmqdYb+CjM4PRhYIGQwiK6xFx0l+vW0AqoSMDpCHmqAVsPNQ6dpMgsCGWFx27wT3UM9u\/ztktL6D9lBDHD8VTcvqgIq2k6EO3eXnS75FRUXQoEGDCu1Y+O6helFoO0LlgR1wMIxRkSiwgxd0w8lLxoRO+equA7rLT4Qaju21PeWL3Vm1ahWL8du5c2fo1asXFBYWUujBkMbJS7NkTIhQddcB3eUnQvViMRXmdbrWYhUcv3\/\/\/ixAfs2aNYV6QddmhGBSmomMCRGq7jqgu\/xEqEpNavJURoQa\/FiQMSFC1V0HdJefCDV4uxtIi0SogcAc1wgZEyJU3XVAd\/mJUIO3u4G0SIQaCMxEqCaYdTeoJD9FStJFB5LqlK\/f5p4I1W+Ey9evaiKd3fVzKD3zFaRUawQZLUcFL0gCLarCIIEuhFqU5CdC1UUHiFB9NjW6KJIdjCrkP7\/vbTi9eWysiYyWIyNFqiow8FlNfa2e5CdC1UUHiFB9NSW0f6ZiIiGZIqnyVLXx\/ZDZLjqhCFVg4LOa+lo9yU+EqosOEKH6akqIUFVMJLOHimSKpBqVpAKDqMhq1U+SnwhVFx0gQvXZUumiSH4u+WLduIf6zeENUKVufqSWe3U64ei3Dvg8VX2rXncboNMcIEL1bRpdqVj3yaS7\/KQDNAdoDuijA0SoRKi+IkDGRB9jQh6qNQI0B\/SZA0SovtKJPopExtRekXQ3qCQ\/7aHqogNEqESoviKgy0RyAlF3DEh+IlRddIAI1Vc6IQ9Vl4lEhEoeOq3SkA4QoRKh+ooAESp9VOmuA7rLr9PBPCJUX+mEjCkZE9KBKOrA1Pd3w96j56BJVjqMv7tZQlYiivInJLBFYV0wIEJVrTmm+nRRJFruouWuiqIDSzcWw\/ClO2LijL+7aUKkqrsNIMY9UwcAAA1TSURBVA\/VZ5IJq3oKjh888mRMyEONmg4gmSKp8tS3Qy7M6dtKevJETX5pQR0K6oIBeah+aI+hTl0UqaJ4J36oA+lAtE65mj1UJFMkVdmk+\/iThyqrOUlejjzU4AeIjAl5qFHUAdxDXVt4HDo3r53Qcq9OZOJkXaKoAzLWkjxUGdQ8lNFFkchDpT1U0gFrBHS3ATp9VBCheiBHmay6Tybd5dfJmBChEqHqrgNEqDIs6aGM7oSiu\/xEqLTkTXNAHx0gQvVAjjJZdZ9MustPhKqPMdXdO6M9VAAiVBmW9FBGd0LRXX4iVCJUmgP66AARqgdylMmq+2TSXX4iVH2MKXmodDCPCFWGJT2U0Z1QdJefCJUIleaAPjpAhOqBHGWy6j6ZdJefCFUfY0oeKnmoRKgyLOmhjO6Eorv8RKhEqDQH9NEBIlQP5CiTVffJpLv8RKj6GFPyUMlDJUKVYUkPZXQnFN3lJ0IlQqU5oI8OEKF6IEeZrLpPJt3lJ0LVx5iSh0oeKhGqDEt6KKM7oeguPxEqESrNAX10gAjVAznKZNV9MukuPxGqPsaUPFTyUIlQZVjSQxndCUV3+YlQiVBpDuijA0SoHshRJqvuk0l3+YlQ9TGm5KGSh0qEKsOSHsroTii6y0+ESoRKc0AfHSBC9UCOMll1n0y6y0+Eqo8xJQ+VPFQiVBmW9FBGd0LRXX4iVCJUmgP66AARqgdylMmq+2TSXX4iVH2MKXmo5KESocqwpIcyuhOK7vIToRKh0hzQRweIUD2Qo0xW3SeT7vIToepjTMlDJQ+VCFWGJT2U0Z1QdJefCJUIleaAPjoQOqGePn0a3nzzTViyZAl89dVXkJ2dDT179oShQ4dCTk5OjLouX74M69atg1mzZsHGjRtZvsGDB8OAAQMgMzNTiOI++ugj6NevH2urY8eOQmUSzaT7ZJKVf+r7u2Hv0XPQJCsdxt\/dLNFhCLW8LAahdlph4yT\/OSgqKoIGDRpAenq6QmSjU5UuOhAqoZaUlMDIkSPh888\/Z0TXvn172LZtGyxevBiaNGkCr732GjRq1Ihpzdq1a2H48OGQn58Pffr0ge3bt7N8Xbp0gYkTJ0L16tVdtYsI1RUi5RlkJtLSjcUwfOmOWF\/G39000qQqg4HygQixQpKfCFUXHQiVUFesWAFPP\/00zJ49G26\/\/fbYlP\/0009hyJAh8OCDD8Ljjz8Ox48fZ2SKXunkyZNj5LlmzRooKChgf+vWrZurySBCdYVIeQaZiYRkiqTKU98OuTCnbyvlfQuqQhkMgupbEO2Q\/ESouuhAaIRaWlrKPFBcxp03bx5kZWXF5jYuAz\/11FOQkpICU6ZMga1bt7Kl3blz50LXrl1j+U6dOgWjR49mZSdNmgRpaWmO9oEINQjzGd+GzEQye6hIpkiqUU0yGERVVqt+k\/xEqLroQGiE6mQwvv76a+aZ1qlThxHq73\/\/e+bFLly4EFq0aBErivuq6J3inuqCBQviSNmq\/jAIdc+ePWxpeuDAgdC0adOKZCeFZJGVH\/dQ1xYeh87Na0d6uRdBksVACOAIZCL59bYBOs2BpCRUXMpFQh07dizzTKdOnQrr16+H+fPnxx1UwoHC31auXAmLFi2CZs2cD6+EQahhtJlMNlZ3+XEsdMeA5A\/+MGQy2QCd5kDSEeq\/\/\/1veOyxxxhx4pJwvXr1GGnipLTyQpcvX85O\/pq9VycPFQ9C4eGmIBKeXMYPgyDbDEIu0TZ0lx9x0h0Dkl9vGyA7B\/Ly8gD\/RSklFaHu27eP7YkePXoUXn\/9dbj66qsZlqoIdf\/+\/YzcNmzYEKUxor4SAoQAIaAdAuiE4L8oJd8JlS\/3GEEZNmwYjB8\/Pg6nHTt2MPAuXbrE9ktbtSo71elEqF6WfLFBJFX8R4kQIAQIAUIgeREgD9VibNwIFQ8W4R1T9Ezx0M6MGTOgcePGcTUtW7ZMyaGk5FUd6hkhQAgQAoRA1BHw3UN1A8gYsOHll18ud+gIyyMpq7g249YX+p0QIAQIAUKAEJBFIFRC5QeQWrduze6R4jUZq3Ts2DEW2AFDdxmjInkN7CALEpUjBAgBQoAQIATcEAiNUM+fP89IdOnSpdC9e3fLO5q49NujRw8WsGHVqlUwatQo6Ny5M\/Tq1QsKCws9hx50A4N+JwQIAUKAECAEZBEIjVAxju\/DDz\/MYvfaJSRaDOyQkZEBVsHx+\/fvzwLk16xZU1Z+KkcIEAKEACFACChBIDRCVdJ7qoQQIAQIAUKAEEgSBIhQk2QgqBuEACFACBAC0UaACDXa40e9JwQIAUKAEEgSBIhQbQbCj4fP8bTyG2+8wR5UP3v2LHvk\/Mknn4Q2bdpApUqVYj1R8Zi6Cv0KEwN8I3fQoEFw4MCBcqJYBQZRIa+5Dj\/kN7aBV8YwiMkvfvGLcg87VGQdEMWgourA3r174dVXXwW8pXDy5Em45ppr2FkQPDNStWrVpLIDfswBUfmTYfy92hUiVAvE\/Hj4HJ+ae+GFF+CDDz5gsYobNmwIv\/vd72DTpk0wZ84cdnqZJxWPqXtVBHP+ZMAAcbrvvvvKXafCKFrf+973EhXRsbwf8vMGkSz\/+te\/woQJE1isUqsY1RVVB7xiUNF0AO\/U4wfht771LejduzeLVY6Pe\/zhD39gpPrEE0\/ESDVsHfBjDniVP8zxlzEwRKgWqPnx8LlVnUiyaFRPnDgBM2fOZMTB79wm+pi6jDIYy4SJAfYDo2Ph60L46IE5claisomU90N+bBfDXmKcavyYunjxIrRt27YcoVZkHRDFoCLqAF4VfO6552D37t3sQY\/c3Ctv\/OIHFq5aodeKr2bdeOONSWEHVM8BL\/Inw\/iL2AlzHiJUEyJ+PHyO8YnxwXRc5sXQitWrV4+1+re\/\/Y3FMMYJ1a5du6SIChU2Bvyd2507d7Il0aCvRfkhP96lxkcfhgwZAlu3bmWeN7aDb4WaPdRkiAwWNgYVUQfww3nMmDGAgWzMscwxyA1ucRQUFDDPNWwd8GP8vcgf9vjLkCmWIUL1gJzsw+dIqHjntlOnTq4TKdnjFgeBAe7b4AdIZmYmC\/6BZJQsSVb+rKwsRqjTp0+HH\/zgB3DzzTez\/2\/1LGFF1QEvGFRUHbDT482bN7Pwqs8++ywj1GTWgUTmgKj8yTz+TraICNWDpZZ9+ByX9vDrEx9N79OnT1yLfOMdPRY8oKTqMXUPYnnKGgQGPOjHLbfcArVq1WIHuXAS33nnnQxD\/MIPK8nK36xZs3JdtntFqaLqgBcMdNEBVAq0D7hy9c4777AlXzwjkMw6oHIO2MmfzONPhKrA+iby8Dk2j4Q6YsQI9vVpTHwZEE\/84jKQqrdfFYhcroqgMODLXzipMPRkt27doKioiO2pIl7mQ1x+yGpVZyLyt2jRwhOhWnmuWMHy5cvZ\/hvuLVvV6TcWQWGgiw7g0ubbb78NzzzzDDusiDYiNTU1ae2A6vG3kz9Zx99tfpGH6oYQACT68HlFINQgMcDlr0cffZTtN91\/\/\/2xK0XFxcXM4FSpUoWRqt1jCgJD6jlLovJXBEINEgMddADJZPXq1TB27Fi466674h7+SMYPa9Xj7yR\/Mo6\/iNEgQnVBScXD53zJ18pDtVrytfNOvD6mLqIAInmCxsCpT7j8+7Of\/Sx2iEuk\/4nmUSG\/l+VOJ2MaZR3wgkFF1wE89INLvLhnes8998CLL74Yd\/gu2XRA9Rxwkz\/Zxl\/UhhCh2iCl8uHzqB5KCgsDHBI8lICPIqSkpMSNEC554lWjJUuWsMAYfiaV8uOBHHOyM5rJdCAlLAwqsg7g9RG8OoUBPR566CF29xR13ZiSRQf8GH8R+ZNl\/L3aFyJUG8RUPnwetWszHJKwMPjpT38K7777LtsnvPrqq2MjhJMbr9EgqQaxh6hSfquTynaEGvaVCeOUCAuDiqoDuFqFe+BIqLjUi2crcM\/UnJJFB1SPv6j8yTL+RKheEbDI78fD53hJGpd3jAdqnAI7hP2YepgY4CnCRx55JO6QBg4T7xNefH\/ppZfiwrQpGPa4KvyQX9RD5YEdKqIOiGJQEXWAH8B5\/vnn2SpLv379LMkUMUoGHVA9B7zInwzjL2NTyEM1oebXw+c4QfCLFC\/1Dxw4EJo3b24bejDsx9TDxgDbx31S3C\/lp3wxusxvfvMbdiAJvVSj5yqj+E5l\/JJflEwwX0XVAVEMKqIOYBCPoUOHMrJEvTYGeOG4YAjSDh06sP8MUwf8mANe5A97\/GVtChGqCTk\/Hz43B8dHT2v06NFw0003uQbHD\/Ix9WTAACcUevVIqhgxCUMx9uzZkxmknJwcWX0XKuen\/MYOOB08sQqOX1F0QBSDiqYDuHyKwRuc0uTJk2NX68LUAT\/mgFf5wxx\/IUNhkYkIVRY5KkcIEAKEACFACBgQIEIldSAECAFCgBAgBBQgQISqAESqghAgBAgBQoAQIEIlHSAECAFCgBAgBBQgQISqAESqghAgBAgBQoAQIEIlHSAECAFCgBAgBBQgQISqAESqghAgBAgBQoAQIEIlHSAECAFCgBAgBBQg8P8Acv4128TXLJ0AAAAASUVORK5CYII=","height":282,"width":468}}
%---
%[output:658f5bda]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdQAAAEaCAYAAACoxaaoAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQm0VsWV7iZehRsBuehVBkUwQBTsoEEDGiI82m6CCSaubhAwzSCKIgoOGMCBURtowESBIDSTdsQx2pFo2o5to40uDfIMaRUk2PJQJi+CCgRRDK\/30fqte+4ZqurUGeq\/31mLFXP\/Xbt2fXuf+s6uscGRI0eOEB4gAASAABAAAkAgEQINQKiJ8ENhIAAEgAAQAAIeAiBUBAIQAAJAAAgAAQsIgFAtgAgVQAAIAAEgAARAqIgBIAAEgAAQAAIWEAChWgARKoCACgKzZs2iRYsWhYo2adKEzjrrLLriiivou9\/9Ln3ta1\/zZPfs2UMjRoyg9evXh5Zt3bo1nXfeeTRkyBDq3LkzHThwgG688UZ69tlnqWXLlrR8+XLq2LFjqfzu3btp5MiR9Ic\/\/MH724UXXkh33XUXNW7cuCSzZs0aTx8\/bNMtt9yi0kzIAIF6iwAItd66Hg3PGoE4QpXtue6664j\/VVRUKBGqKMukvGDBAurRowfNnTvX+29+5s2bRz\/4wQ9KVbz22mseWTLx8tOqVSuPdDt06FCSWbZsGd1xxx2B5bPGDvUBARcQAKG64CXYWBYI6BAqE+P9999PXbp00SJUBqp79+4ekTJpcmbJz+jRo+mmm24KJEvxxyVLllDv3r29\/\/vpp5\/SpEmT6JFHHqHjjz\/eI9szzzyzLPyARgCBtBAAoaaFLPQCAR8CMqGuXLnSIz752blzJ91555301FNPeX+eMWMGXXrppbUIlQl26dKl1Lx581LRv\/zlL7R582a69dZbad26dXTsscd6ZMzDt8OHD6ft27d72enMmTO93w4ePEgTJkygVatWUcOGDenzzz+nw4cP09VXX00333wzNWjQgOQh4a5du9LixYupqqoKPgUCQCACARAqwgMIZIRAHKGyGQ8\/\/DBNnDjRs4iHbC+55JJYQhXm+\/V36tSJrr32WuK50NNPP90jYp5P3bFjhzcnu3HjRvrhD39IH330Ef3Xf\/2XN0w8f\/58atq0Kb3++useGX\/wwQf0D\/\/wD3T77bd7w894gAAQCEcAhIroAAIZIRBFqHxgGWeSnEVyhirPacqLkoIyVM4uecESD9Fu2LChVLZ9+\/ZelstDuUyGDzzwAJ177rn08ssve\/OnXI7nVt955x1vQZJcJ9vAc7hyppwRTKgGCDiLAAjVWdfBcNcQUJ1DPfrooz0i5OyUh19VVvnKWAwaNMgjVx7OlYnxH\/\/xH2ngwIFeFioT6K5du0qrecXiJbGgSQwfn3322a7BDXuBQOYIgFAzhxwV1lcEVAmVt7DwFpW2bdt6UOkQ6kUXXUSTJ0+m6upqr6x\/6Ja30vA8KW+nEVtl9u7d6w3v\/s\/\/\/I83vCvL\/NVf\/ZWX4Qp99dV3aDcQUEEAhKqCEmSAgAUEVAmVq5K3v6gQKpPjuHHjiId5xf5V1iMvLuI5Uh7GveGGG7zhZSZOnmOVFymxzE9\/+lO6\/vrrPYIdMGAATZs2jY455hgLCEAFEChvBECo5e1ftK5ACMQtStq\/f783RMtDs\/v27aO\/\/du\/pTlz5nhbWMTBDmIOtVmzZt68KZMd\/y\/PkfIWGSbMysrKUqvl7S88R9qvXz\/vcAmW55XAYqWx2HPKW2SuvPJKr16eY73tttvo8ssvLxCKMAUIFBcBEGpxfQPLygyBOELl5vJBC7ylhYlVkCf\/3U+oYtvMpk2baNSoUd7CIn7Gjh3r7TmVV+Q+9NBDdU45klf9cjn5oAeeN2U7MH9aZgGI5qSOAAg1dYhRARD4AoE4QuVs8ne\/+51HfpyhqhAq63366ae9IVrOKHmomDNQeY+r\/1QkLuMfyvUfRcgyftKFH4EAEIhGAISKCAECGSGgM4cqstLx48fTxx9\/HJqhshwTKQ8Tr1ixwmsJz4PefffdpYMY5H2noqlixa\/4\/6xj+vTp9C\/\/8i8lNHiVMcvxamE8QAAIxCMAQo3HCBJAwAoCOoTK21Tuuece4kPv4\/ahsnHbtm3zDrvnfaj88EpePvmIt90cOnTIy3qfeOIJ77ewoVx5iw3LYf7UituhpB4hAEKtR85GU\/NFQIVQeZiVt77wkYNiq4oKoXLLnnvuOY9EOdvksnwykjh\/Vz7onm+04aMETzjhhFqAyFts\/IuW8kUOtQMBNxBwilC3bt1Ks2fPpn\/\/93\/30OVVkPwl3qZNGzfQhpVAAAgAASBQtgg4Q6h8XBpvCTj55JNp6NChxEe18ZAY77njr+1vfOMbZeskNAwIAAEgAASKj4AThMorHvnqKV6cwSR64oknesjyEBVvJ\/jJT37ibVDn+SI8QAAIAAEgAATyQMAJQuVl\/8OGDfPItGfPniWcxAkvnKXyRci82AIPEAACQAAIAIE8EHCCUHljOp\/qIq6fygMo1AkEgAAQAAJAIAoBJwiVV0e+8cYbNHXqVHr88ceJCZbvcOzVq5e3oZ3vfcQDBIAAEAACQCBPBApPqGIPHV+GzAuR+ILk\/v37excf8y0YfP4pL0rC9VJ5hhHqBgJAAAgAgcITqnwTBu+x4xsyxDmlO3fu9Fb+8oHefJh348aNYz363nvvEf8TT4sWLYj\/4QECQAAIAIHiIMD9vHwmdXEsC7ek8IQqjkTj80qXL19e2qgumsSXJT\/66KPesWvt2rWLxJyJlPetvvLKKyU53oIzZMgQZV+xPX\/+85\/p61\/\/unPO5kbCfmVXpyYIH6QGrZJi4K8EU2pCqvjzjUpVVVWp2ZGG4sITKjd67ty53rFpTKgdOnSohcPDDz9M8+bNC\/zNDxjvZR08eLB3OAQf6caPbobKGfOuXbu8co0aNUrDJ6nqhP2pwqukHD5Qgik1IeCfGrRKilXxR4aqBKe+EJ8xypciy\/c3shaeU2VyXL16tdIKYEGoK1eurHUbh45Fn3zyCfFh4zyX6yKhwn4db6cjCx+kg6uqVuCvilQ6cq7jH4WKExkqz5VeddVV1L59e2+lr5gr5SFc\/nvXrl29g7yPOeaYyAgAoRK5Hsyu288B6nobYH86RKOqFfirIpW9nBOEyrBwlso3Zpxzzjm1VvnyoQ4LFy6kjh07xqIHQkVnHhskGQigQ8wA5IgqgD\/wTwsBZwiVh3f\/+Mc\/eqclrVmzxhtu7du3L11zzTXKh+ODUEGoab1IOnrRoeugZV8W+NvHVEej6\/g7P+Sr46woWRAqCNVWLCXR43qHAvuTeD95WeCfHMO0NDiTodoAAIQKQrURR0l1oENMimCy8sA\/GX5JS7uOPzLULxEAoYJQk3YGNsq73qHAfhtRYK4D+Jtjl3ZJZKiaCCOYNQGzLO46\/gyHjTb4T\/yyDHOkOt6Yv3v3buKN9y5uHYP9WUZL7br4PusTTjjB6a2HyFCRoZZiwEZnnt\/raIeM8rTfBqEGnfiVd5tQPxBQQaBbt240ffp0Ouqoo5zdyw9CBaGCUFXe9oxkkn7UBJ34lZHpqAYIGCPAR77efffdtGzZMjr11FNBqMZIFqQg5lDdz\/CSklERQjFpG2zEcRFwgA31CwERtyDUMvG7jY4oaWeYN5SwP28PJP+osRHH+aMAC+obAiDUMvO4jY4IhJRvULiOv405VBtxnK8XUXt9RACEWmZet9ERud6hw\/78gzqpD2zEcf4owIL6hgAItcw8bqMjStoZ5g0p7M\/bA24O+e7Zs4dGjBhB69evLwHYqlWrwGsT+XquCRMm0KpVq2qBHXXL05\/+9CcaPny4Jx90TaOsKEz\/jBkz6NJLL63jYKF7+\/bt3m98ocb48eOVA0H0G1EFwupWrqQeCIJQy8zJINTknXneIeH6B4GLQ76CkK677rpahMV3EU+cOJFkMhHE271791qkJd69MDKbNWsWsUxNTQ3169cvlPBEnX49wsbq6mrvKsfmzZt7oSr+PmfOHO\/KxjD7ouLaRr+R93tThPpBqEXwgkUbbLwYrnfosN9iQBmqSuoDG3GsYzqTHWd3M2fOpMrKylpF\/b8x4T300EO1SE0UCPtNkNzAgQNpy5YtHrHKpCjKhxF72O8ik+XfZdtZ\/7hx42IzYaE3a7x1fOOSLAjVJW8p2GrjxUjaGSqYmaoI7E8VXiXlSX1gI46VDP1SSGSPQSTn18OkOW\/ePGWy4vJyGf7\/PPQrMkpZv4od8+fPp06dOlHv3r1L2SgTtTwULBN40BCxv01Z463jG5dkQagueUvBVhsvRtLOUMHMVEVgf6rwKilP6gMbcaxk6JdC8hxi1Dwoi8vzlSrziv4sknXw\/Ks\/qzQZqvUP94o26+rKGm8d37gkC0J1yVsKttp4MZJ2hgpmpioC+1OFV0l5Uh\/oxvGsZ96hrXs+oTbNG9H4Pu2UbAzKPHm+VH54rjNoGDhoEU+XLl0ih3HljDQoy40b7g1qVBihChLnRVUqi5PiFiWFtc0I6DIuBEItM+fqdkRBzU\/aGeYNKezP2wPJF4bpxPGDa3fS6Ac3lBo9vk9bY1IVSnjoddGiRSWdUZmoWEQkhP2LiYKGcQURyouTwgjVv4KX6xEE98EHHwQOH4cRahgB6+Cdf3QV1wIQanF9Y2SZjRcDhGQEvbVCruPPQCRtg04cM5kyqYpn0LktaMGgM6z4Qwyd8spcna0uYtg4aCuObJic+akO08oLn8IINUwXk7t45MxVB28rwJapEhBqmTnWxouRtDPMG1LYn7cHsiVUf4bKZMqkqvqEZW2ivPxOdezY0dur6l8EJGT9WWbUAqagd1VlUZJMqFxvkD1Bi5L4b1OnTqUxY8bQHXfcQbfccgt16NDBM91Gv6GKdznLgVDLzLs2XgwQUr5BEYT\/3kem0uGaLVRR3ZaqBkzO10CF2pPGkG4c8xzqms0fUo\/2zbSHe8O2nsiEKm9BiSI9mZw5+wxafCT0CtJr06ZNaZ5WZR5VJlTe4sN1+OdKg7bN8N+ef\/55b06VdfAjVgDr4q0QAvVSBIRaZm638WIk7QzzhrTc7N+3egXVLPjihB1+qvpPKTypJvWBjTjWiUNBZF27dq21CCmI4AQRsn55m41\/mFWVHP1bcETb\/Qui5OFjeZ5WyIt53qjh3p49e3qHP7Bt9913H916663evtus8dbxjUuyIFSXvKVgq40XI2lnqGBmqiLlZj+TKZOqeJr0GkbVo5enimFS5Ul9YCOOddsQdNxf2NGDrNu\/cIn\/Jm+5URm+DVqcJOwO0h92ClPc0YNxxyrGrfJlm8JWPOviXM7yINQy866NjihpZ5g3pOVmvz9DZTJlUi3yk9QHNuK4yPhkbRsP8b700ku1sm\/+G5\/apLKtJmt7Xa0PhOqq50LsttERJe0M84a0HO3nOdSDb66myk69Cj\/cy\/5P6gMbcZx3HBalfpF5n3\/++bVOU2KMOQtWOR2qKG0puh0g1KJ7SNM+Gx1R0s5Q02Tr4rDfOqTaCpP6wEYcaxtdpgV4OPinP\/0p\/dM\/\/VNpVS83Vfd4wjKFx2qzQKhW4cxfmY2OKGlnmDcKsD9vDyBDzd8DsCAPBECoeaCeYp0g1OSdeYruUVLt+gcBhnyV3AyhMkQAhFpmTgWhglCLENJJPwpsxHERcIAN9QsBEGqZ+dtGR5S0M8wbUtiftweSf9TYiOP8UYAF9Q0BEGqZedxGRwRCyjcoXMcfQ775xg9qzw8BEGp+2KdSMwg1eXaUimM0lIJQcbasRrhAtEAIgFAL5Ay\/Kc899xxdffXVdP\/993vHhak8IFQQqkqcpC2T9KPARhyn3UboBwJ+BECoBY0JPsFk5MiRtHnz5lrHmcWZa6MjStoZxtmY9u+wP22E4\/Un9YGNOI63srZE3PF8snTQMYX8u3z0oL9+cTwg\/13nKjhZT9i9rHFHD+piEWb79u3bA1VF3RebtG6d8v67aYN8wrElX3agoz9OFoQah1AOvx86dIimTZtGTz75JB04cACEqumDpJ25ZnXWxV23nwFJ2oasCTXsIHvRQcuEEXb4vLA57LxdcbYv360qXyzuDyBRp1+PsLG6urrW6Ub+6+dU71XVCdyog\/7j2q1Tj6ls2IUFQbaBUE1R\/qJcgyNHjhxJpiLb0o8\/\/jgtWLCAhgwZ4t1fGPXV67fMRkeUtDPMFq26tcH+vD3gHqEy2XH2NXPmTO\/2Ffnx\/yZfn9a8efNasmG\/yacS8egTv6dBR\/7F3VDj\/z3s6jnbpBFnl8pFAGlFpcr1e4MHDy71o7axkduFDDUtLxvq3bRpE40aNYpuuukm4pdVDgQVlSDU5J25Cs5pyrj+QeBihqpDCFGXhofFhVyGZYYPH05z5sypszZCxY758+dTp06dqHfv3qHHB+oeKyj6DWG\/\/2aZOEL1Y+LXF3Rrj3+YmuuWRwLE73wR+j\/\/8z\/T+vXrPfP8mXtcn8dY7N27lxo1akStW7f2Pmb8Q762bVm2bBmdeuqp1LJlS6\/ecnqcyVD3799PkydP9hwwadIkeu2110CoBpHoOiG5br+LhCoTQNyIkNz5qswd+jMoxifo4nGToVr\/cK94XXR0iSFm0e6gjC+OUOUPAU4K\/ImAP3MPsttPjDLOwragdul+4PgJ1aYtyFANOuw0ivCo9COPPOINAy1evJjatm1rdOmvcOjYsWOpW7dunqktWrTw\/qk+3KHv3LmTTjrppDrDX6o68pSD\/Xmi\/0XdSX3w6quvelMeceQmWsq38Ryu2UIV1W2Nb+MJWtASdgeoPwNjO7p06RI5jCtnpEEkEEdaQV4NI1RBipwZRl3PFka8fltU5lDFx0XUkLhoQ9AQu9\/msDr9hKiS1cvYBZX3D\/eb2iLiYuHChV6GGteHVlRUEP9z6XEiQ3399de9oYzp06d7Qzn8xA1lBDkh6EUfOnSo1zmpPrwoihdO8OKHhg0bqhYrjBzsz98VSX3Aw3s33nijEqH674ut6j\/FmFTlDn\/RokUlIKMyUT8R+4ckgzr8oIvFwwgkaDhSkPcHH3wQOHysSqhhfUwYoYSt8pU\/fHQzeH+fJfAL+1jw45SUUP1ky9m1eHRtEW3heOGPmbg+tFmzZlRVVZX\/C6thQeEJlcf3OaNs37498XyB+GJJQqizZ8\/25gtMMlR+mXbt2uVltS6O\/8N+jbcjJdGkPuDYv\/zyy5UItWbBcGJSFQ9fvs6XsNt4RAbHH5g6W138Q5Ri\/s9vk5zVqg7TyhlgGKGq6gr6AJdt9BPKddddV+tO1TCMg7YVBc2PCoLm3y6++GJvKFxk1arXztka8rVhi8CTRxlPPvnk2D4UGaqNt9SnIy6oWTxoUj8qQ1UdKgvS4focHuxPIUg1VSb1gc7HpD9DZTJlUlV9wjIhUV62pWPHjjRixAgaOHBgILH4s6eozj6ojSrZlkyobGOQPf5FSX6CE0PZTPQqezJNhqNl\/LldnPHLxMm\/y6uqw7Ji\/+KtsOHooEVeXId\/TlhuM5Ne0Hy2qS2YQ1V961KU2717N23cuLFODW+99RbdeeeddOutt1Lnzp29lX1NmzaNtESnIwpTlLQzTBEqJdWwXwmmVIWS+kA3jnkO9eCbq6myUy\/t4V6VbRcy6USRnkzOnH0GddYCeEF6bdq0KRGLCnHJhMpbfOSsTv4I8NvMv\/F8qqiX\/zvsA8FPyCp2RQWUnDFfeeWV3kcAn\/4mz+\/6h8FV51BV\/MfDuCJDludQjz\/+eKu2gFBT7VaSKdftVLg2kzJ+K5N2hslanbw07E+OYVINSX1gI4512iA6765du9bKmoI69bBDBPzDrCokFJTBirb7F0TJJznJ87TyvN2ll15aIkw\/YQk8BAENGjTIIzX\/Kl+W8380qLRF6I9qE4+ciQ+NrVu3lhZxyRm0f5iZ9Yrh9qC5Z\/497mAHGUuZUEWGassWEKrOW5exrEmnYlIGhJqxY2OqS0pGRWhN0jbYiGNdHILm\/aKmWsQwplyPPNWiMnwbRhCC1OSFUfy3sFOY\/AuXwk5ZkhcVybb6p51096H6sY47AjAIa7bn+eefLx16IeaHOZvl3Q\/yHCd\/OAQ9cT4RSYecvdu0BYSq+9ZlKG\/SqZiUAaFm6FSFqpKSkUIVqYskbYONOE69kY5UIAjj\/PPP9+Z9\/RlqUZsRN7edpd2qtoBQs\/RKBnXZ6IiSdoYZNDOyCtiftweSn1ZlI47zR6EYFvgJ1RVsVUksC5RVbQGhZuGNDOuw8bKAkDJ0WEBVruPPTUraBhtxnK8Xi1W7PKTLw8H88OExYUOnRbBelcSysFXVFhBqFt7IsA4bHVHSzjDD5gZWBfvz9gAINX8PwII8EACh5oF6inWCUJN35im6R0m16x8EyFCV3AyhMkQAhFpmTgWhglCLENJJPwpsxHERcIAN9QsBEGqZ+dtGR5S0M8wbUtiftweSf9TYiOP8UYAF9Q0BEGqZedxGRwRCyjcoXMff5pCvfGtSvl5B7UAgHoFt27bRzTffTLgPNR4rJyRAqMmzo7wdDUIleu+997yO6ZVXXsnbHagfCGghwNdm8q1hRx11FC4Y10KugMIgVBBqEcLSxkcBkyr\/y+Nh+z\/88EM64YQTnLuvUowQlLv9s57ZQu\/u\/YROqWpE4\/u0zSNMAuvk4ww5bnbs2AFCLYxXDA0BoYJQDUPHajEbhGrVIE1lsF8TMMvicfg\/uHYnjX5wQ6lWJtTxfdpZtsJcXZz95przL1n4+1BtQgRCBaHajCdTXa53KLDf1PN2ysXhz2TKpCqeQee2oAWDzrBTuQUtcfZbqCI3FSBUTehdDwbYr+nwFMThgxRA1VBZ7vj7M1QmUybVojyu4x+FIwhVM8pcDwbYr+nwFMThgxRA1VBZH\/Cf9cw7tGbzh9SjfbNCDfeym1zHH4T6JQIY8nU\/mMvhZXS9DbBfg71TEAX+KYBqSSUyVE0gEcyagFkWdx3\/cvhCd90HsN\/yS6mpznX8kaEiQy3FgOvB7Lr9IFTN3jcFcddjCPanEBSWVCJD1QQSwawJmGVx1\/EHoVoOCAN1rscQ7DdwekZFQKiaQCOYNQGzLO46\/iBUywFhoM71GIL9Bk7PqAgIVRNoBLMmYJbFXccfhGo5IAzUuR5DsN\/A6RkVAaFqAo1g1gTMsrjr+INQLQeEgTrXYwj2Gzg9oyIgVE2gEcyagFkWdx1\/EKrlgDBQ53oMwX4Dp2dUBISqCTSCWRMwy+Ku4w9CtRwQBupcjyHYb+D0jIqAUDWBRjBrAmZZ3HX8QaiWA8JAnesxBPsNnJ5RERCqJtAIZk3ALIu7jj8I1XJAGKhzPYZgv4HTMyoCQtUEGsGsCZhlcdfxB6FaDggDda7HEOw3cHpGRUComkAjmDUBsyzuOv4gVMsBYaDO9RiC\/QZOz6gICFUTaASzJmCWxV3HH4RqOSAM1LkeQ7DfwOkZFQGhagKNYNYEzLK46\/iDUC0HhIE612MI9hs4PaMiIFRNoBHM8YDtfWQqHa7ZQhXVbalqwOT4AhoSruMPQtVwdkqirscQ7E8pMCyoBaFqgohgjgZs3+oVVLNgeEmoqv8Uq6TqOv4gVM0XLgVx12MI9qcQFJZUOkOohw4dolWrVtGyZcto48aN1KRJE+rZsyeNGzeO2rRpowQHLhhP\/4JxJlMmVfE06TWMqkcvV\/KPipDrnQkIVcXL6cq4HkOwP934SKLdCUI9ePAg3X777R6hDh48mL73ve\/Rjh07aMmSJfTZZ5\/RvffeS2eeeWYsDiDU9AnVn6EymTKp2npc70xAqLYiwVyP6zEE+819n3ZJJwh13bp1NGzYMLrttttowIAB1KBBAw+Xbdu20ciRI6lz5840ffp0atiwYSReINT0CZUdwHOoB99cTZWdelkd7i0HMiqHNqBDT7tbjtYP\/PPFP6p2Jwj1scceo8WLF9OCBQuoQ4cOpfYcOXKEZsyYQWvXrqWlS5dS8+bNQagxsYaXMf+XET7I1wfAH\/inhYAThBrW+E8\/\/ZQmTZpEb731FghVMULQmSgClaIYfJAiuAqqgb8CSCmKuI6\/8xlqWAOYSEeMGEF9+\/al8ePHU0VFBTJUZKgpdgV2VLveocB+O3FgqgX4myKXfjlnM9T333+fxowZQzU1NbRw4ULq2LFjLFqYQ81mDjXWEQkEXO9MuOmutwH2JwhgC0WBvwUQU1LhJKHu3bvXG+p94YUXvHnVHj16KMEjCHXs2LHUrVs3r0yLFi28f6oPB\/POnTvppJNOosrKStVihZGD\/fm7Aj7I1wfA3w38ecQxbtQx35bUrd05QuWVvRMnTqT169fTPffcQxdccEFp1W8cuIJQZbmhQ4fSkCFD4oqWfuf9sJwVV1dXx64qVlaaoSDszxDskKrgg3x9APzdwL9Zs2ZUVVWVr7GatTtFqBs2bCDOLj\/\/\/HOaM2cOnX322VrNFYQ6e\/Zsat26tVGGyntid+3a5WW1jRo10qq\/CMKwP38vwAf5+gD4u4E\/MtQU\/bRp0yYaNWoUNW3alObOnUunnXaadm2YQ8X8nXbQpFAAc2ApgKqhEvhrgJWCqOv4R0HiRIYqFiBxZnrXXXfRKaecYuRmECoI1ShwLBdyvUOB\/ZYDQlMd8NcELEPxwhMqH97ARwvyMC0fOXjWWWfVgee4447zTlBq3LhxJHQgVBBqhu9WaFXoEPP1AvAH\/mkhUHhCPXDgAE2YMIGeeuqpUAy6dOmCgx0UIwSdiSJQKYrBBymCq6Aa+CuAlKKI6\/g7P+Rry7fIUJGh2oqlJHpc71BgfxLvJy8L\/JNjmJaGwmeoNhsOQgWh2ownU13oEE2Rs1MO+NvB0VSL6\/gjQ\/0SARAqCNW0E7BZzvUOBfbbjAZ9XcBfH7OsSiBD1UQawawJmGVx1\/FnOFxvA+y3HNSa6oC\/JmAZioNQNcFGMGsCZlk8C\/z5PtfDNVuoorqt9ftcQaiWA8JAXRYxZGCWchHYrwxV5oIgVE3IEcyagFkWTxv\/fatXUM2C4SWrq\/pPsU6qabfBMuR11MH+tBGO1g\/888U\/qnYQqqZvEMyagFkWTxt\/JlMmVfE06TWMqkcvt9oguG6wAAAgAElEQVQKfxvSzoitGo8ha9twautL+x3QNkizgOv2g1C\/RACLkjB\/F\/fu+zNUJlMmVZuP3KF89vJDqWfENm1nXa53iLDfdkTo6XMdfxAqCLUUA2kHs81sK0iXDfvjbOTfD765mio79bI+3OsnpH1LR6WeEet1d\/HSNnwQX0t6ErA\/PWxVNLuOPwgVhJoJodqcfwzTlfRltGmjSucRJBOVoaaREZva6S8nPkSoqjXtP28EtWzZ0skbl5LGkC08TfWY2j\/rmXdo655PqE3zRjS+TzvT6hOXM7U\/ccUZKMAcqibIrgdDmvb75x8bde5Frab8pybCX4iHzWUmtT+LOdK4BvvbsPWadqVVxW1+8U5c8Vx+93+I0N+MpVZDZoJQDb2RhNxM3oEH1+6k0Q9uKFk7vk\/b3EjVxH5DmDMvBkLVhNz1YEjTfj9ZMbRRq2Sjhl75t72PTil5RywOCrI\/bghXdnEWc6SivjC7XJxDrePbc\/6OWo39JQhVs\/9g8aTkZvIOM5lyveIZdG4LWjDoDAPrkxcxsT95rdloAKFq4ux6MKRpfxChBq2SFUQjr6b1E69fF+8Jrezci\/zDjX6C5Kz46Oq2kXtI4+ZIS0ObRMZ7UaOGlmUf+OdQRTvT2gOrGe4lcX97GgyYTS37XesUoYqssFXTChrcqSK3Ieuk5Bb0DsdlvEym6++dQC0\/30k7jmpBXa6eSUyqeTxp9kF5tEeus94Tqk52w8C5FAxpLeoJC9o6w4JfEhIToSAIf+YpdPmHh4PIuVSvNNwYJecnc5koWVcQaQXZF\/VRIGMh64sa\/pZj6OCTs2pl4rK+JHtgVeJaRYbtkT+A2E8NTutGn3zvqtwISafTFEQjMkNRduR3mtGUH5+RyweBINSRH93nEdxp3zyd+k38eWCzgojS3wepZLw6awe4zu9vWUStDu\/y3pHFxw2pM\/caR+ByY\/yyUX2ojl6dOMhKtt4R6r\/d1Jd+dME5dPLJJ3sYy1mSyG7iwN+3fx9VNmpEFRVHx4nm+ntY28LsD8uK\/EQU1Si5Tr8c4\/vJG6tDi8v4f1azJVy26mSq7NzTw\/\/gG6u9+cewR9YZZJvf52H268oF2V\/Lli9j6LO3Xgy1v5SVG0RRmO9ZJz+MmWrsB8ltOngs3X\/SVVRRURFqnX\/xi04cGTQ5sMjKtTsC\/96qSQX1aN8sl3d4656DtObtD6nfgWdKtoX1PbL9Pb7RjNo0r\/yi35L6oDVv7\/UITzyMe49vVNVqt\/89iYotrlO2bV3DLrS94otsVtgQZlcQ2EGyYX2QX\/a0b56Rykp7W\/Hl11OvCDUsO0oLXBf1+rOxoKzTxXbBZvsILG46lBYfNzRS8Xc9EmhEIz++j5o8O8u+EdBY9ggkGanJGpx6Q6ggBvXQkgM4cuhVXSUkyxCBVcf2oanNx8e27If\/m4lN2QMyjQUKAoEIpHFaWVpQ1xtCBTGoh5AcwMBNHbf6JqlKqJP3zKo1hFjfcEJ7kyEAQk2GXyqlkaGqwyofLgBCVcetvklOaT6efnNsn9hmg1BjIYJABAIg1IKGh21y4In9qAUxBYUh0iz\/atugDxGW4efw+1ustF\/oi1qw5CKWss3cRlt45YUFx3vFiW3p3T2f0BMHvhk7fyrsdIVQxcKbVoe\/2q8ZhDUv0uGn66H1ka5guZaf76I4far+ZPtUdLHcjqNOCrWPRxYu+fpbke+uaKOo00Y74uzn+Kpp2Y2qd7xSy7ZNfX5G37\/ielWYcpWrN0O+jDIfjs+rfEf+zVnGnVvFN3vQ4c8+o6qufanixFNrHWwue5K\/qjhA+ExYfnTIIoiobZMYz5P6V3qynUFH34l9m\/y7\/3xb+TedNgqsZAIXumwSTxiRyatdw97AqA8mQS6MBz\/s5yi7xbw0t3HfH5\/1YqhRZWWdMoz\/4ff\/X6A+vz0s+8kbz9daresncd6TG7XyWsjL9YZ9LInY4G0P\/Ra8Ruu2fbWyVOjhRUgvvv1hLUiTEmpcR8yVMQHEEVxcT8tEs7z9ZLp146g6uvg3JrN1Dc\/yPiS4nZ3\/793U9dAfAkmzpmV3uuObv6B+f36GLvz9xMiqOdNn3T\/88zO1CFOQ2hftO4u6HfchdXn313HNoGXtJ3sjB99\/Z5Fnn4wf63nj22O9hWK8SCzIfq73qhN\/5smIvaq8NWf45qm1hu5lTL74wPhDnQ8NWYbtiJpLF8Qp9syyPrY3zz2zsWD7BOodoQ4ePJhWrlxJnT\/ZWIcMozpZ7pi586y8eDzt2LGjtAcvjFD8xCSOl5PxD+qUBWHJ8vIiobiha6FTkF9YZy938IL0kx4G75HF6hW1vi79HbWfdIOGc1RGEridXJZtDyNyGbftU\/5PLTku26hzz9APIj8+frIMG4YK+yiQ5f378KIOmvD\/Jq9Ur0XSPtzZ76JOWYeIB\/lDIMzvQR9SYuvLg5uPoikNLiuFs+h8+YxY3ks465mvtjIFEWrcFipWzD7e2+VSGlrzoxI58N\/92RLruqr6Lo\/g\/KQU1iFyR88ELGd8TGx7v3UpVf3x4Vodv381M5PpqtFne0f5idOH\/G2U3\/9VM673Yk8QG5OEICD+730XjvfIi3UJEhTELeznowKZTDs+c0OpSUx8si5BvExA\/MhHDco4MEn2aF9V+t1vu5gb55OU5MMfWJ8gYfnDgnXLH1G8v1aQobwKnPVx+\/iD8l8PfJO2f3y4RMICB3F6E8fQms0felub8jx3GIQagYB8fVsQoUZ98YvOK2pTcljHGEaCYQehR23CjiObMJ1ZnmEbRRB+O4LsDTsgQh5eF2QRNiTtJ4mwIwfDyKZqwORakaSLX9QRh0kOBwmLjSAcbB+0769DJhr5KDv\/QQP+DluMSsiEzR1s9Y6Xa2WGwsdx+lju2e\/M8AgibkUxk9pvvt6nNFwtd\/4ic+P6wkiBg0ImGs7UBZnIZTh+BBH47fd3UUJfkByTtyAVf+YWtGVJPqNXtk2uU9Qnfvdjxljyh0UQkQmiY30y2ckEyL8JMpT\/W+gT8T9zzX569LXdJdP8BK5LZkWQr7cZ6jfWLaqzsZ0Pcg8jv6izZOMcGXSMHuvzd9pCT1TnHXTUHpNH3HVjKkQW1w5bv8cd\/cf1+LMj\/\/C6IIsgPMIO5FepN6yNUQQZViasviSEGhYbNi8mCGuPvw55la\/cGfqP1vN32GFE\/29Lfl4rA5Pl5I788ren1hr6FO8my\/izPM6U5KHg9af8iBafeiud+HWi9bsO1zoQwZ+5CRxE9h2UMYWRpf+sXNl+1icOYvBnYGGEJWzh31dv3E2HDx\/2Tnni8kGkxX\/z2yaTs8hgRYYd9jFg652X9cjxf\/fzO5zMRMNwqbeE6s9Q5Zc3bOiSX1yTzlC3M46SDyL8uI3Pdcqc1o1aTV2dy7FrSV7Q0vxjm66lm050M8ek9cd9uKjoN4khoTcsNnRjTMVOv4y\/Dp7zWnCoZ51huaCO\/I6Kh73FJnHTCiofPXFtlUnJP4fJ7\/nR3Qd60zb\/+S7RuH\/9amhafBQIYpZPH4q6nSVIPs1sSyd+ooZOVT8GTGIlqoyO\/bbrTltfvSXU7t27l7KguDkk+XfTYFDpKGRnh8mrHkAv6yrnm0LiOte0XyAT\/aYxJOoKi40sroEL+qgJwiDtOTCd98kvq5IhmRxgn3abBc468RN3Nm7WHwPcBh37Td6vPMvUa0I1AT7vYDCZK\/OXcfGmkKjORKdzNfG57TJpxFDUvLsL9tu2MWmG5M\/e0sw4dduuGj8qh+aLurP6GACh6nq7wPLyoiTOUE0e1WA20a1aRl4BzCshVS6lFqTj0k0hQXgUAX9VP4XJpdGGLIe+07A\/KaY65VXsl1cq53kZd5J3wCTL1sHRVFYFf1PdeZdDhqrpgbyDwWQOVW5i3vZrwl1HvAj2i60jpneWptGGLIe+07A\/aVzolI+zXyez06nXlmyc\/aKeombZqvbbwitLPSBUTbTzDgbdOVTR+ZeaWdWa9p83wom7LJN8nWu6VVncxtBqWjGU1dB3WvYrOyGhYJz9Rc3sRLPj7Gc5MXcqVvHy34qSaavYn9DFuRUHoWpCn3cw6Myhhh4CIV3Qrdn83MXzxt\/G0GrebUjqxHK3v6iZnSqh5rV6VzWuXI+fqHY6Q6hHjhyhF198kebNm0dr166l448\/ni6\/\/HIaMmQIHXvssUq+LJc51KATbIIACD0E4py\/o1Zjf+ncthluY94vo42h1bzboPSyRAiZ2B+32jSpTTrlVezPcpGOju0q74A\/wxb6i7KwSgV\/XUyKIu8Moa5Zs4ZGjx5N3bp1o4EDB9Kbb75J9913H11wwQU0depUaty4cSym5UKosQ39UiA0QwWhqkIYKBd0\/J+OQtc7FF37izYnqWu\/jm+zkI2zP+5AhyxsjKojzv687UtSvxOEunfvXo9MOSudMWNGiTyff\/55GjNmjPe3iy66KBaH+kaoDEjQIRUEQo2NlSiBpPOoqh1K0sVPiRppMUMt2pykKv5p4ZdUr4r9LmfYSfHJs7wThMpEyEO79957L\/Xu3buE1\/79++nGG2+k5s2b07Rp0+iYY46JxLI+EioDUs4HO+Tx8iSdR1XpEJOSdpq4qNgv11+0OUld+9PE0kQ37DdBLZsyThDqQw89RPPnz6fly5dThw4dSsjwvCpnpzynunTpUo9Yo576SqjlfrBDNq\/KV7UknUdV6RCTknaamKjY76+\/SBmTif1p4qmrG\/brIpadvBOEOmvWLHrppZdoyZIlVF1dXQsd\/u3pp5+mFStWULt27UCoIQjgYAe7L1WSLSoqHaI8T8uWh10XZ7dVatpU7FfTlI8U7M8Hd1Gr6\/hHoecMoXJ2GZSFPvzww97KX3\/2GtRokaGOHTvWW9zET4sWLbx\/qg8Hw86dO+mkk06iyspK1WKFkYP9+btCxQd7F11BB9c8UDK2ssdlVHXVkvyN\/3KlNd6B\/FyhEj\/5WRdfs6r9FRUVxP9ceuolocoOGjp0qDc\/q\/ocOnSIampqvEy5YcOGqsUKIwf783eFig+OPHIz0au\/+srYc\/6O+AzmIjwq9hfBzjAbYH++3lHFv1mzZlRVVZWvsZq1O0+oJkO+s2fPptatWxtlqAcPHqRdu3Z5WS3fR+jaA\/vz95iKD4qcoarYnz\/K4RbA\/ny9o4o\/MtSU\/IRFSfaAdX3+wnX72ZMqbUi68MlexNTVpGJ\/mvUn1Q37kyKYrLzr+Ee13okMFdtmkgWwXNr1YHbdflVCZbkkC5\/sRQwINU0sTXS7\/g64br\/zhCoOdmjZsmWtU5FwsIP+6+h6MLtuvw6h6ns3mxKu+wD2ZxMnYbW4jr\/zhMoN+O1vf0vXX3899ejRg\/r370+bN2\/G0YMG74Xrwey6\/SBUg6C1XMT1GIL9lgPCojonhny5vUGH41922WXeAflNmzZVgqS+HuyAIV+l8MhMCB1iZlAHVgT8gX9aCDhDqDYAAKGqLYixgXVaOlzvDJGhphUZ6npdjyHYr+7rrCVBqJqII5g1AbMs7jr+IFTLAWGgzvUYgv0GTs+oCAhVE2gEsyZglsVdxx+EajkgDNS5HkOw38DpGRUBoWoCjWDWBMyyuOv4g1AtB4SBOtdjCPYbOD2jIiBUTaARzJqAhYib3vXpOv4gVDvxk0SL6zEE+5N4P92yIFRNfBHMmoAFiCe569N1\/EGoyeMnqQbXYwj2J42A9MqDUDWxRTBrAhYgnuSuT9fxB6Emj5+kGlyPIdifNALSKw9C1cQWwawJmEKGWj16uXffp8rjOv4gVBUvpyvjegzB\/nTjI4l2EKomeghmTcBCxE3PqXUdfxCqnfhJosX1GIL9SbyfblkQqia+CGZNwCyLu44\/CNVyQBiocz2GYL+B0zMqAkLVBBrBrAmYZXHX8QehWg4IA3WuxxDsN3B6RkVAqJpAI5g1AbMs7jr+IFTLAWGgzvUYgv0GTs+oCAhVE2gEsyZglsVdxx+EajkgDNS5HkOw38DpGRUBoWoCjWDWBMyyuOv4g1AtB4SBOtdjCPYbOD2jIiBUTaARzJqAWRZ3HX8QquWAMFDnegzZtn\/WM+\/Q1j2fUJvmjWh8n3YGiOoVsW2\/Xu3pSoNQNfF1PRhgv6bDUxCHD1IAVUMl8P8KrAfX7qTRD24o\/WF8n7apk6rr+EeFGghV40VEdqEJVgri5fAyut4G2J9CYGuotIk\/kymTqngGnduCFgw6Q8MafVGb9uvXnm4JEKomvq4HA+zXdHgK4vBBCqBqqAT+4RkqkymTapqP6\/gjQ\/0SgZdffpkGDx5MK1eupO7duxvFjOvBAPuN3G61EHxgFU5tZcC\/NmQ8h7pm84fUo32z1Id7y2GUD4QKQi3FADoT7f7XegH4wDqkWgqBvxZc1oVdxx+ECkIFoVrvFswVut6hwH5z39soCfxtoJiODsyhauKKYNYEzLK46\/iXw5CX6z6A\/ZZfSk11ruOPDBUZaqYZKt8kc7hmC1VUt6WqAZM1X7do8XJ4GV1vA+y3GtLayoC\/NmSZFUCGqgk1gjkasH2rVxBfIC6eqv5TrJKq6\/gjQ9V84VIQdz2GYH8KQWFJJQhVE0gEczRgTKZMquLhi8P5AnFbj+v4g1BtRYK5HtdjCPab+z7tkiBUTYQRzHoZKpMpk6qtx3X8Qai2IsFcj+sxBPvNfZ92SRCqJsII5njAeA714JurqbJTL6vDveVARuXQBrwD8e9AmhLAP010k+kGoWrih2DWBMyyuOv4g1AtB4SBOtdjCPYbOD2jIiBUTaARzJqAWRZ3HX8QquWAMFDnegzBfgOnZ1QEhKoJNIJZEzDL4q7jD0K1HBAG6lyPIdhv4PSMijhDqIcOHaJVq1bRsmXLaOPGjdSkSRPq2bMnjRs3jtq0aaMEF87yJcLLqBQqqQrBB6nCG6sc+MdClKqA6\/hHgeMEoR48eJBuv\/12j1D5cPvvfe97tGPHDlqyZAl99tlndO+999KZZ54ZGwQgVBBqbJBkIOB6hwL7MwiSiCqAf774O0+o69ato2HDhtFtt91GAwYMoAYNGnht2rZtG40cOZI6d+5M06dPp4YNG0YiDUIFoRbhVUSHmK8XgD\/wTwsBJzLUxx57jBYvXkwLFiygDh06lLA4cuQIzZgxg9auXUtLly6l5s2bg1BjIgWdSVqvkrpe+EAdqzQkgX8aqKrrdB1\/5zPUsAZ8+umnNGnSJHrrrbdAqIrx7Howu24\/u8n1NsB+xZctJTHgnxKwFtQ6kaGGtZOJdMSIEdS3b18aP348VVRUIENFhmrhtUhXBTrEdPGN0w784xBK93fX8S\/LDPX999+nMWPGUE1NDS1cuJA6duwYGwWYQ0V2FBskGQi43qHA\/gyCJKIK4J8v\/mVHqHv37vWGel944QVvXrVHjx5KCAtCHTt2LHXr1s0r06JFC++f6sPBvHPnTjrppJOosrJStVhh5GB\/\/q6AD\/L1AfB3A38ecYwbdcy3JXVrL8yQL2+NmTBhgrc1RjytWrWi5cuX11qIxCt7J06cSOvXr6d77rmHLrjggtKq3zhwBaHKckOHDqUhQ4bEFS39zvthOSuurq6OXVWsrDRDQdifIdghVcEH+foA+LuBf7NmzaiqqipfYzVrd4pQN2zYQJxdfv755zRnzhw6++yztZorCHX27NnUunVrowyViX\/Xrl1eVtuoUSOt+osgDPvz9wJ8kK8PgL8b+CNDTdFPmzZtolGjRlHTpk1p7ty5dNppp2nXhjlUzKFqB00KBTAHlgKoGiqBvwZYKYi6jn8UJIXJUKOMFAuQODO966676JRTTjFyMwgVhGoUOJYLud6hwH7LAaGpDvhrApaheOEJlQ9v4KMFeZiWjxw866yz6sBz3HHHeScoNW7cOBI6ECoINcN3K7QqdIj5egH4A\/+0ECg8oR44cMBbrPTUU0+FYtClSxcc7KAYIehMFIFKUQw+SBFcBdXAXwGkFEVcx9\/5IV9bvkWGigzVViwl0eN6hwL7k3g\/eVngnxzDtDQUPkO12XAQKgjVZjyZ6kKHaIqcnXLA3w6Oplpcxx8Z6pcIgFBBqKadgM1yrncosN9mNOjrAv76mGVVAhmqJtIIZk3ALIu7jj\/D4XobYL\/loNZUB\/w1ActQHISqCTaCWRMwy+Ku4w9CtRwQBupcjyHYb+D0jIqAUDWBRjBrAmZZ3HX8QaiWA8JAnesxBPsNnJ5RERCqJtAIZk3ALIu7jj8I1XJAGKhzPYZgv4HTMyoCQtUEGsGsCZhlcdfxB6FaDggDda7HEOw3cHpGRUComkAjmDUBsyzuOv4gVMsBYaDO9RiC\/QZOz6gICFUTaASzJmCWxV3HH4RqOSAM1LkeQ7DfwOkZFQGhagKNYNYEzLK46\/iDUC0HhIE612MI9hs4PaMiIFRNoBHMmoBZFncdfxCq5YAwUOd6DMF+A6dnVASEqgk0glkTMMviruMPQrUcEAbqXI8h2G\/g9IyKgFA1gUYwawJmWdx1\/EGolgPCQJ3rMQT7DZyeUREQqibQCGZNwCyLu44\/CNVyQBiocz2GYL+B0zMqAkLVBBrBrAmYZXHX8QehWg4IA3WuxxDsN3B6RkVAqJpAI5g1AbMs7jr+IFTLAWGgzvUYgv0GTs+oCAhVE2gEsyZglsVdxx+EajkgDNS5HkOw38DpGRUBoWoCjWDWBMyyuOv4g1AtB4SBOtdjCPYbOD2jIiBUTaARzJqAWRZ3HX8QquWAMFDnegzBfgOnZ1QEhKoJNIJZEzDL4q7jD0K1HBAG6lyPIdhv4PSMioBQNYFGMGsCZlncdfxBqJYDwkCd6zEE+w2cnlEREKom0AhmTcAsi7uOPwjVckAYqHM9hmC\/gdMzKgJC1QQawawJmGVx1\/EHoVoOCAN1rscQ7DdwekZFQKiaQCOYNQGzLO46\/iBUywFhoM71GIL9Bk7PqAgIVRNoBLMmYJbFXccfhGo5IAzUuR5DsN\/A6RkVAaFqAo1g1gTMsrjr+INQLQeEgTrXYwj2Gzg9oyIgVE2gEcyagFkWdx1\/EKrlgDBQ53oMwX4Dp2dUBISqCTSCWRMwy+Ku4w9CtRwQBupcjyHYb+D0jIqAUDWBRjBrAmZZ3HX8QaiWA8JAnesxBPsNnJ5RERCqJtAIZk3ALIu7jj8I1XJAGKhzPYZgv4HTMyriLKE+99xzdPXVV9P9999P3bt3V4Lr5ZdfpsGDB9PKlSuVy\/gVI5iVoE5NyHX8TQl17yNT6XDNFqqobktVAyanhq+KYtd9APtVvJyejOv4RyHjJKFu2bKFRo4cSZs3b9YiRxAqkevB7Lr9JoS6b\/UKqlkwvPQeV\/Wfkiupuu4D2J8eWapodh3\/siLUQ4cO0bRp0+jJJ5+kAwcOgFBVIliScT2YXbffhFCZTJlUxdOk1zCqHr1c0\/P2xF33Aey3FwsmmlzHv6wI9fHHH6cFCxbQkCFDaOrUqSBUzYh2PZhdt9+EUP0ZKpMpk2pej+s+gP15Rc4X9bqOf9kQ6qZNm2jUqFF00003UfPmzbXnQ20M+fJw83333UdDhw6ltm3b5huZBrXDfgPQLBcx8QHPoR58czVVduqV63AvQ2Fiv2UIE6mD\/YngS1zYdfzLglD3799PkydPpkaNGtGkSZPotddey4VQbZBy4ohMoAD2JwDPUlH4wBKQhmqAvyFwloq5jr\/zhHrkyBF65JFHaOnSpbR48WIvMzRxiigzduxY6tatm1F4bNu2jW6++WZKosOoYkuFYL8lIBOogQ8SgGehKPC3AGICFar4n3zyycT\/XHqcWOX7+uuv01VXXUXTp0+n3r17e\/iaEOp7773nkeErr7ziko9gKxAAAkCg3iHASQv\/c+kpDKEePHiQJkyYQKtWrSrh16pVK7r77ru9f+3bt6dbbrmFKioqjAmVCzKp8j88QAAIAAEgUFwEkKEm8E0Yod5www1eVhn1MPEuX76cOnTokMACFAUCQAAIAAEgYI5AYTLUsCbs3r2bNm7cWOfnt956i+6880669dZbqXPnztSpUydq2rSpORIoCQSAABAAAkAgAQKFJ9SwtpnMoSbACUWBABAAAkAACEQiAEJFgAABIAAEgAAQsIAACNUCiFABBIAAEAACQMBZQoXrgAAQAAJAAAgUCQEQapG8AVuAABAAAkDAWQRAqM66DoYDASAABIBAkRCot4T6\/vvv05gxY+j666+vc9k4XxH3zDPP0MKFC4m357Rr144uu+wyGjRoEFVWVtby39atW+kXv\/gFPfHEE97f\/\/qv\/5quvfZabxuP\/PDxiS+++CLNmzeP1q5dS8cffzxdfvnl3q05xx57rHZMZG3\/n\/70Jxo+fDht3769jq18itX48eNj28DX7fGF8HzBOx8\/xhj8+Mc\/9u62ra6uLpXXwWrv3r20bNkyTy\/vZebL5vnyhG9961vUoEEDI51hDcnT\/qLiL2O1Zs0amj9\/vvc+8OUVtuM\/DfxV7S8q\/tz\/zJkzh55\/\/nnat28fnX766V6\/0q9fP2rYsGHh41\/Vfhv4x3ZQFgTqJaG+++67HgHw1hvu3LkTFs\/hw4c90uN\/F154IQ0cONAjkXvvvdc7\/5evjGvcuLEnLm6\/+fjjj71bcPiM4d\/85jf0H\/\/xH94Vcz169Cjp5c5m9OjRng7W+eabb3q31lxwwQW1dKr4NC\/7uY2XXHIJVVVV1TLzjDPOoO9\/\/\/uRptfU1HjHiPGLMXjwYPr2t79N\/\/3f\/+1h0KZNG7rnnnuodevWng5VrMSFCS+88IKHPx\/w8eijj9Krr75qHf8i2F80\/IXD+QOIY37ixIne2at85rafUFV9GhZEaeCva3\/R8Of+iz9mTzvtNLr00kvpxBNPpKefftq7K5pJlQ\/FEaRaRPx17U+Cv0q\/akOmXhEqf+FyJslHGX7wwQcefn5CXbduHQ0bNswLyOuuu6501CFnqiNGjPCyKc4qOYu9\/fbbieU5k+3YsaOnjwn5rrvu8sHSM6MAAAhaSURBVLLRRYsWUYsWLYizKCZTzshmzJhRImT+quQsmf920UUXxfozL\/vZsIceeoiWLFninUh1yimnxNrqF+AjJfnoSM5gevbsWfqZz2lmXH\/yk594mf2HH36ojFWQTiZZ7tg\/+ugjz89M\/jbwz9P+ouLPdvExnhz\/\/CHDsd+lS5c6hFpU\/FXtLyL+ov955513vI9\/7mf44Y8bHq3hrHXFihXUtWvXQsa\/jv028NfusAwL1CtCnTVrlkdyfC7wxRdf7BGfn1Dnzp3rka7\/KEPuLPhw\/h07dnjl+IuZifcHP\/iBdzSiPLzIJMHDo1wfH+bPX2JMwpzlisP92V\/c+d94443e1\/y0adPomGOOiXRjXvbzS8qkzydWMSHqnkj1l7\/8xctAxUeGnL3wRwKf4fy1r32NZs6cSevXr1fC6vPPP\/fK8TAv+0OMGjCAzz33nJcNc8dy9tlnJ8Y\/b\/uLiD\/H6p49e7yPIfYZj1wwTnzXpT9DTRr\/aeCvY38R8ecPxnHjxnlTS\/7pFh454\/6HP9Y5cy0i\/jr2J8XfkBuNitUrQuUvOe7Med6Ohxt56NFPqExaHIBBw1b82+9+9zuPbPnWeQ5azmI5aOVHjPfz0C5nXZzdMRH5SVoECs+pBtXn92he9gvS47leFeLXiUQeLmeMOJNkQv31r3+thBUT6hVXXEHnn39+bIdiC\/+gdmVhfxHx5\/eICXX27Nn0ox\/9iL7zne94\/x307hQRfx37i4p\/2HvGd0XzB\/xtt93m9U1FxV\/V\/jTx1+mrVGTrFaHKgIQdXcgZ6m9\/+9s6Q5siQ+W5IiZGvvWGM9T+\/ft7hCA\/IkP9+7\/\/e6+zZyJ+6aWXvCFTefENl+HfeN6Dh2d48ZPqk6X9nI0zeZ133nl03HHHeYuAmEjCFmCptoHleNib8eMsnzsBVazYH\/xBw2X5wyXog4azJl6gpKpTB39RXxb2FxH\/IKzCPkaLiL+O\/a7gzzEpppx+9atfeX0Kr29wBf8w+9PEX6evUpEFofoWJXEHeeWVV3rDuwMGDCgN5Yo5VAaVCbVly5bekAvPxcpzGCKgeXhXrH6Nynoffvhhr7zubTlhhJqG\/WIIiQObh8p5vpeHvvkDgbMU\/wIslcBjGbGoiz8yeEiYF1WoYsXlw0YIxFAkLzYTHzRhow6m+GdpfxHxD7rZKYpQi4a\/jv2u4M8jXo899ph3YQgv4BFrQFTfKd3bupK8v0F1hdmfFv6q\/ZSOHAjVR6g8J8eLjXilnEweDzzwgJddvv322yXyE6vUePiIMzgmWc42xZ2rvHQ9rQ49jFDTsJ+HkK6++mrvA4KzbjFfvHPnTu+lPfrooz1S9a\/+jQpEXqnM88dMfvKiLtWXP29CzdL+IuKvQ0iqPtXp0JPir2O\/C\/gzGfF0FI\/08O4EeTeCC\/hH2Z8G\/jokqSMLQvURKoPHpPTggw9685qclYm9jbwox59N8vaXn\/\/857R69WpvOFTsAWOiUcmQbA\/5pmF\/VEDx8O\/Pfvaz0gIgleDbsGGDt2iI50F5bpmHpcQTN4cthsfFkG\/UHLY85Bs1L6475J61\/UXEX2fIVNWnqkPuNvDXsb\/o+POiLR7i5TlTXiQ5ZcqUWgsHi45\/nP228Vfpo0xlQKgBhBoGJs+v8n6uuAVEvJSd51evueaa1BYFmFxfZ2o\/48ELA\/hQC16NKz88ZMrbVPyLu4Iw5K9Qxo8zU96zy6tz\/VtwVBdQ5LEoKS\/7i4i\/f58p2xjWcav6NEinHEc28dexv8j48\/YTHuHhwzS4z+G9p\/7DZ4qMv4r9tvA3JUmdciBUH6EyUfEwLQ9hnnnmmSUsxV46PoSAV7p+9tlnNHnyZO+UIx4i5kVK4hF7Fm1t2whyaBihpmE\/E\/Hjjz\/uDXWL\/bZsE3dwnGEyqarMAcuby\/lyeP8CLdapusQ\/y20zAv+87C8i\/kFbvMIIVdWncdvGbOKvY39R8ReH0DCh8lAvrymQ+yERt0XFX9V+W\/jrEKOpLAjVR6hiywsfNcjzhjxfKCbLJ02a5O0l5YMJ+G+8TeDZZ5+lxYsXexkXP+JIwGbNmnmbq3l\/pCBjnmOV5zZ0D3aQnRxGqGnYLxY6yQsd2BaxKIE3j\/MiLvmoM39AClneN8cfJGHzrTpY8YcLD3PJi6KiDnZIgn+e9hcVf7+PwwhVx6dhHVka+KvaX0T85T6JR4h4C2AQmXIbi4i\/jv028DclSN1yINSAowe5g+YhFP7i4\/11v\/\/9770MjAlWPj2JFyjxyUl80MHQoUM97PkoPd5SIi+04b\/zVhw+N5iPI+StNps3bzY+elBkckH7aPmrz7b9PCzD86Q8XyoWavGw9i9\/+UtvQRJnqXLm6g9CLs8kyvPSvFBLfHzIcjz0y7o5c1DFijsK\/jLngwUYfz6wI+zoQVWdQS9Q3vYXFX9VQkoa\/2nhr2p\/EfHnAzS47+F3gN8b+WAT0S7ua84991zv\/xYt\/nXsT4q\/LikmkQehBsyhsgO5Y2ZS5O0hnTt39rbS9OnTp04WxqTKX+a8KInnLvr27evNnfLQsPwEHfjOWTAvYtI9eSiKUPm3NOxnnZwRMqny4qywg+2DglHsI+PDNMIeJlo+2IEx1MHKfzg+Z8s8R3vOOefEHo6vin8R7C8q\/rI\/oxa\/6PjUHyNp4q9qf9Hw5+Fv3rcd9fDpZuLQmaLhr2t\/EvyTEKRu2XpLqLpAQR4IAAEgAASAQBQCIFTEBxAAAkAACAABCwiAUC2ACBVAAAgAASAABECoiAEgAASAABAAAhYQAKFaABEqgAAQAAJAAAiAUBEDQAAIAAEgAAQsIABCtQAiVAABIAAEgAAQAKEiBoAAEAACQAAIWEDg\/wPJxFsIMjQmjwAAAABJRU5ErkJggg==","height":282,"width":468}}
%---
%[output:09638136]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ18VMW5\/x8kNQGhCAgkECRRAlhaQW0vEVCSCqIoICqCXFuIQDEtwlXKa6gh9kZALgrU3sAfaEjr5aXWXpX6goAJahGuWsTWokBLkHcjoCIkKVH+PBNmnd3sy9ndOXPO7v5OP9TkZM68fGd2zm+feeaZRufOnTtHuEAABEAABEAABEAghgg0goCJod5CVUEABEAABEAABAQBCBgMBBAAARAAARAAgZgjAAETc12GCoNAaAJfn\/mCjhQPotrdf\/YkTu7Sh9IKXqKLmn47dAZIAQIgAAIuJwAB4\/IOQvVAIFwCUrw0bp1OqQ+vFY\/Le\/wzREy4RJEeBEDAjQQgYNzYK6gTCERB4KsTh+jowhGUOmUdNW7VwZNToPtRFIVHQQAEQMAxAhAwjqFHwSBgH4Evtz1Ln69\/0svacvSJkZTc8TvUcvgj9hWMnEEABEDAEAEIGEOgUQwImCLAlpZDBb2p7tOPAxaZdNnl1KF4q5eFxlT9UA4IgAAI6CAAAaODIvIAARAAARAAARAwSgACxihuFAYCIAACIAACIKCDAASMDorIAwRcRuDkM4\/Syd8XNqgVlo5c1lGoDgiAQMQEIGAiRocHQcCdBNgH5vAvb6G2+Svo5J+epGa9h1Oz7LsITrzu7C\/UCgRAIDICEDCRccNTIOBaAup26S82r6TaA38X8WCwjdq1XYaKgQAIREAAAiYCaHgEBNxMQAata9rzZmraYyAdWTic0qY8Q2d2bqBTr63E7iM3dx7qBgIgYJkABIxlVEgIArFDwNcKI\/1h2k75g1hOwgUCIAACsU4AAsaFPXjixAkaO3Ys7dy506t2PXr0oJUrV1KrVq1o3bp1NHPmzAa1nzt3Lg0YMEA8n52dTdOnTxdp5s+fT8uWLaMJEyZ43du2bZsnT38o9uzZQ3l5eXT48GGvP69evVrkH+ziZzds2EATJ04UydTfZb4PPvggjRgxgqqrq2nJkiU0fvx40T7fZ13YTagSCIAACICAgwQgYByEH6hoKWBUAeJ7jwXMr371KyotLaWsrCyRlbxXUlJCK1asEPfmzZsn\/jtjxgxav349SRHUpEkTcU+m4d+DCRhVaPBzH3\/8cVDh41tff21Sy2OBJcUU3\/cVYC7sJtdWCb4uru0aVAwEQEAjAQgYjTB1ZeXvZc8WClVwvPDCCw0EDAuAUaNGEVtHtmzZIgQLCxy+pBWlffv2XvekMAlUd19LCadTy2FBJMUR\/40tQEOGDPG6d9ttt4nsX3zxRfHfwYMH07hx4yg\/P5+4fL6kNSktLY26du1KFRUVnrQswgIJLF3M4y0f7DiKtx5Fe0AABHwJQMC4cEz4EzBSSPDLn5eFgllgWLQcP37cI2b2798vxM6sWbPoscceE6KhU6dOnr8HWwoKJWBYKEnLye7duz15dunSxcuK4tsm33xhgdE3EIMdJYA4MPo4IycQAAFnCUDAOMvfb+mR+sBI6wovKakCobKyUviwTJs2TfijsGDJyMhoYMHxVxl\/vipyCWnBggU0depUj6+NrPfIkSMb+OFAwLhwoKFKIAACIBDDBCBgXNh5wZaQ3n33XbEE9Je\/\/MUjQLgJvETUpk0bj1+KzOPyyy8X\/irSn0ZaOvg+X6GWZ\/w58fouQ\/k6+LKjMDvjqn4sEDDmBlogHxj4xpjrA5QEAiBgPwEIGPsZh11CIIdXddlIFTBscZG7knx3GfHOI77kriF19xKnnTRpkpe\/iu\/uIn9LSLJBwRxzQznxYgkp7GER8gErp1And+lDaQUv0UVNvx0yPyQAARAAATcTgIBxYe8Es8DI3T8bN270WgKSz1RVVXl2Jkmx4m9pia0m7HDLW5iDXcEEDD+n+q7AB8YdgwmWFnf0A2oBAiBgLwEIGHv5RpR7IB8YVYj4c+KVu4PY0ZeXhjiODO9Kkr\/zTh65m0kuRckt2IEqGkrAyPx4xxNfqiiSsWfk1u3ly5eLWDT8++zZs2ny5MnCoZhFlKy7tBaxc7BMK2PfRAQTD4EACIAACMQlAQiYuOxWNCrRCfBp1Gfee1UsF515fyN9svBuuqhpC0or2EDJXXolOh60HwRAIA4IQMDEQSeiCSCgElCXkPi+PJmaf\/6kZBy1\/8Ur1LhVB0ADARAAgZgmAAET092HyoNAQwIsYKRoOXviIH2+\/klhiTl7cBcEDAYMCIBA3BCAgImbrkRDQOAbAl9ue9Zr2YiFDC8j4TBHjBIQAIF4IQABEy89iXZoJSAFAGfa8p4iajn8kQb5h0rD4fz5Sn14rda6ITMQAAEQAAEiCBiMAhDwIaAuwfCf\/PmNhEojxU3T60cYFzAIZIchDQIgkAgEIGB8evmvh7+kO0reo5Onz9KAq1rRb0Z3p0suTkqEsYA2XiDA4uNE2cPUoXgrNUppTkeKB1GLwQ9Rs+y7PIyCpZEC4lsdv0Nfn\/nCmICJlUB2Bw8eJP6HCwTilUB6ejrxP1z2EoCAUfgeO1VLo0v\/Rgvu7krfa9+MpjzzId363TbU\/6rW9PGJGlrz9hG69wdpdHmrFK9eqauro1OnTlHz5s0pKSm02LE7PVfO7jLcln80bW5SfZyq3\/wfap47mpLaZBCLE+n4yvmygGna82avZaRgaXjpqFnv4XT2wAdUe+DvxgSMHJRuDmTHwoXPz9q+fbu9MxtyBwEHCfTq1Yv4rDiIGHs7AQJG4btp13F6+W9VtHB4twbU39z7GQ357x206J6uNPy6dl5\/r62ppSNHj1Jqu3aU0sRb3PjrPrvTc5l2l+G2\/KNp82UnP6JTCwZR2pxyatI9JyoB862O3enLrc8I0cKxWJwQMPZOGdHlLgMW8uTeoQO2ckdHE0+7kQCL88WLF3uOb3FjHeOlThAwPgJm+ZsHqfJ4De355IzXEpIUMOP7dqCf3OBtGqytrSUO4c+HKSYnJ4ccG3anFy9zl9XJ7vpE1eaTH1JKyVAvARPpEhKLlzNvrfMaA074wfDSFVuOanf\/2VMXN5yDJAWM75lbIT80SAACMUIAY9xcR0HAKKyXbjlApW8dphd+1pOaJTem+8s+oNwureiBfh1JCpg\/PtCDemW28OqhmuoaOnrsGKWlplJySmgBY3d6rpzdZbgt\/2jafOnxD6l64W0eARPKQZfLspLGKQuMFC+NW6d7LV\/x0tZXxw86epgjJndzkztKcoYAxrg57hAwfiww0nGXfWD44iUlKWBe+Ok11Lfzpd4CpqaGjhw5QmlpaZSSEnoJqcbm9OJlbnMZbss\/mjZfeuIjOjnvVo+A4bzULdIydoqvb4m\/NOrAcErAuHkXEiZ3c5M7SnKGAMa4Oe4xKWDkTqE7rr5MiAt151BW26bCgtKueWhLiC9m1Yn3istS\/FpgIGDqqcW7gDH3EbSnJLa21Ox8xXP2Ue3u7XSkeCCl9LjFuFOx2kJM7vb0N3J1DwGMcXN9EZMChi0jWW0vEUs7p\/9VJ4TG+L7pYrcQLwPt+eS0X0dcK1jZkfee5e+LpHnXp3nygQXGmx4EjJXR5GwatgCd\/H2hpxKBAvKZrCUmd5O0UZYTBDDGzVGPOQHDVpJJaz+kJSO7CSsLW1+m\/uEjKsv7rvjd9++RomQhVL77hCcODAQMBEykYwnPfUMAkztGQ7wTwBg318MxL2B8hYYOASOXpL5\/eXMImABjERYYcx\/SSErCLqRIqOEZEIieAARM9Ayt5hBzAoYbJpeQhl3blob8+j3Ku769WE6Sf+P\/+ovlYgUKL0n9x9qPqFtaM\/q\/ys8gYCBgrAwbV6XBLiRXdQcqk2AEIGDMdXhMChi2srBw4Vgt0k9F3ut7ZYuIxQtjZ4sOX53bNiWOCSN3JCGQnfegjOdAduY+fvaUhF1I9nBFriBghQAEjBVKetLEpIDR0\/SGufDS0ZJN+2nRyK701j8+9ytgEMiunpvdgenCzT+qOvkEsrNrfJnMN153Ia1atYrKysoEyoyMDBo9ejTl5OSYRIuyQCAoAQgYcwMkpgWM3IG0cdcJD7FoDmBk68us5\/d60ef81o3vgUB2PmMyngPZmfv42VtSvO1CKioqojlz5jSANmbMGCotLQ0b5okTJ2j58uU0adIkatKkSdjP4wEQ8EcAAsbcuIhZASPFS\/qlyZ4lI3mP8UV7ijRvp\/a3hIQ4MPWDM96deNUgdYG2H\/tL43sitBu2LpubTkKXFOnkXlFRQbm5ucLqwoKlX79+tGXLFo+gKS8vD9sSAwETur+QInwCkY7x8EvCEzErYALtNtKxC4mHBQRM8A9HPAsYK8cEBEpTtWoKJXf8jji5WgaPuyx\/JTXLvguzDRFFOrmzeGERwxaYwsJvYtvwklJeXp4QNVasMNXV1TRjxgxav349DR48mNq3by8sMC+88ALNnDlT9NHcuXNpxIgRnrryvR49etDKlStp9+7dNGrUKJFuwoQJNH36dNqzZ4+ow+HDh8V9PueJ08tyZJ5DhgzxusfpunTpQmxZevfdd8XzXKd58+bBIhTDn5ZIx3gMN9mxqsesgPEnMvieGuQuXKq+S1IIZBeYYDwLGLasRHqYoypU5G6gFoMfMi5g4m0btRQwvpaWyspKyszMFJaZffv2hfzIr1tXf9CmFChr1qyhcePGEd8vKCgQwmH+\/PnCwvPOO+\/QwIEDKSsrSzzDFhsWGyygWrVqJZ7hq1u3bkKI8LP88mLL0J133kkbNmygiRMneurE6bm+LHpkXg8++CAVFxfTrFmzKD09XfzMfj2yzJANQgLXEYCAMdclMStg1J1IgXCFe6yAGsVX5i+3aCOQnTfleBcwn69\/Uhx6yBef6ty0583CqiIvFjmh0rAF5pOScdT+F69Q41YdjH2q43EbNVs42NrCVha2tshLWmDYkZfFTahLipPs7GwhItgHhn\/m\/NWLrTDXXnutx7LCFpXZs2fT5MmTPZYWTs9WGLbgqNYWaZnhspYtWyay5fz46tSpkyiPLUEsVoYNG0abNm3y+OE89dRTXqIpVHvwd\/cRgIAx1ycxK2BMIFKtORAwEDDhCBh1iSm5Sy8Tw9VTRjxuo5Y+MNxIXkaK1AdGtcDw0g\/vaGJrjGqB8ddZ\/PfTp0\/Tzp07PRYYmU4VHdICw1YWeUlrC4ugqqqqBhaYZ599FgLG6CfE3sIgYOzlq+YOAROAtXqw4\/faN\/PsQlp0T1cafl07r6fiKSZKart2lNIk9Ina8dTmy05+RKcWDPKcRh3tEpJTlhd1UPIOpFOvraQOxVuF9Uc6Fzf\/4VgvS5K5qaa+pGgmd2mF8a2zr19MsDapPjAsKLp37y6WjlQfGH6e\/VN4KUhaUPz5wKiWFek\/I\/1qBg0aRPn5+R5rDVtgAvnAqDuhYIExPSL1lxfNGNdfm\/jOMWYFDAuMx1\/ZR48O7UyXXJzk6SX2Y3nk+b007ZbMiE6k5ox8D4jke9ICgzgw9ajDjdNid\/qo6uQTByYaJ966Tw\/SiWeKqN1Da+mipt82Onv47oDyV3jSZZd7RI3Ryl0oLNrJnS0xbDVhXxL2e2FLjLqk5ESbUCYIqASiHeOgaZ1AzAoYbqI8s+iOqy8TW6l5yee59z+l5\/J7EltNIrl8LS8yDylg\/vhAD+qV2cIr63iKiZKWmkrJKckh0cVTmy89\/iFVL7zNY4HhxqtbpNtO+YNwwvVdmvFN0\/TqAcJfpnb3n734yedDQk2ABJjcE6CTE7yJGOPmBkBMCxiJibc837P8ffr9+Kup\/1WtI6YXbAs2fGC8scazE2\/EAwgPhiSAyT0kIiSIcQIY4+Y6MKYFjG4LDFtwSt864kX\/saGdxUGREDAQMOY+lvFbEib3+O1btKyeAMa4uZEQswLGLh8Yac3hLpDihX+GgIGAMfexjN+SMLnHb9+iZRAwpsdAzAoYBqUGnuMzi+YP60IPrP47Lbi7a0Q+MKr\/C+c\/9Q8fUVned4UzMAQMBIzpD2c05amB7JK79KG2+cvp6MJ7qG3+CjK9rVttBwRMNL2KZ2OBAMa4uV6KWQGjnoU0pk+65xTp3711hMp3n4joLCS2vhQ8v5de+FlPapbcmO4v+4DG900XfjUQMBAw5j6W0ZWkBrJreftDnl1Rn7+4iM6896oI0Gd6h5RsESb36PoWT7ufAMa4uT6KWQGjOtx+cuqsR8B8WfsVTVr7IS0Z2S3sbdTq+UfcBSxgcru08vKBQRyY+sEZz3FgzH387ClJ3S2lbus+V3OKji4cQalT1hmNDAwLjD39jFzdSQACxly\/xKyAYUTsdHvws1p6qH8nWvnGISoYdAWNWPE+9b2yheeE6nBQWhEwiANzQcDU1oqoom3atKHk5NDbrmMpDkw4Y8ataY8+MZK+On6QWg0vpM9efJIuG72QjhTfQhdnXU+pD691rNqY3B1Dj4INEcAYNwSaiGJawDAm1emWf1cdb8PFaGUJCXFg6qkmUhyYlvcUORq9NtxxLNOrcWr4nhvagck90t7Ec7FCAGPcXE\/FvIDRiQpOvEcoLS2NUlJCHyUQz3FgrETi1TnuEikvTO6J1NuJ2VaMcXP9DgHjw1q16KiB8eDE6w0qngWMlbOQzH1E46skTO7x1Z9oTUMCGOPmRgUEjEXWUsAUXFdLfa681Oupuro62rx5M3Xr1o06deoUMke703MF7C7DbflH0+Y+aUl0xcbZXoc5fr7+SbFbhy8+HqBpz5tjZhkp2JlIsX4WUsgPFxKAgMMEIGDMdQAEjEXWH5+oocm\/eo52PPMEXXTmU4tPIVksEGh3cR09clt3GjCjhJLaZIhzkGJZwARizo69zXoPF+c6OXVhcneKPMo1RQBj3BTpOHDiNYeK6ODBg+IfrvgjkJ6eTvyPr3hdQvI9jNKJXsTk7gR1lGmSAMa4OdqwwJhjjZJihEC8OvGqwqxxqw6O9AYmd0ewo1CDBDDGzcGGgDHHGiXFEAF1C3LbKX9wdNklXGzBfGCcbks0k3teXh5VVlY2wJGTk0OFhYXhYrKc\/sSJE7R8+XKaNGkSNWnSxPPcunXrhM9bdna25byQMP4JRDPG45+O3hZCwOjlidxAAASCEIhmcs\/NzRUChgWLvCoqKigjI4PKy8uNc4eAMY48JgqMZozHRANdVEkIGBd1BqoCAvFOIJrJnQUMi5XS0lIPJmmVsSpguPySkhJ64403aPXq1bR\/\/36aOXOmyG\/u3Lk0YsQImj9\/Pi1btkzcmzBhAo0fP95jgdm5cyeNGjWK2rdvL4TU7bffLvKQlpg9e\/bQhg0baOzYsTRjxgxav369SMt1zsrKivfuRfuIKJoxDoDhEYCACY8XUoOAawkEWzqSlY7lbdS6BMyWLVto+vTpxGKjrKyMCgoKxNIQC5drrrmGduzY4bVcJJeQ7rvvPnr88cfFchWnZ4Fy7733+hUwfMQGXyyI+PmioiLxXKtWrVw7flAxPQQgYPRwtJILBIwVSkTE5mI2X\/PEx5f6LU1+c+P7nE5+o+NvbzK9HNTqNz2LRUeVTK1Pjx49aOXKlWISDVSfcO9HVTmLD4fThurqas83X\/kN2uk+sNjMhEgWzeSuS8CwxYSFhTrWJXz+LPOlfoalBWbQoEHi8y0Fj1xC8meBOX36tMeKw\/nBCpMQw1s0MpoxnjiU9LQUAsYCRylWpCDhAcr3WAyo38Rat25N06ZNE9\/S+Gc2I\/PLs0uXLjRlyhSaNWuWKO2xxx6jhQsX2v5tjL9hqmVxnQ8fPizqOHv27Ab14br5q2eg+ya+TYbbho0bN3qEJn\/z5T4YOXIkDRgwwJE+sDC8bEly8plH6eTvGzq2xroFJlofGP7sSgHja4Hx7QgWw8XFxTRs2DDatGkTqRYYHvv8eerXr5\/Ijy8pitjCw0td8p4tHYxMXUsAAsZc10DAhGAtv2XxpMQXCxJ\/1hiesNSXJAsYKQaOHz\/uV\/CY3r3AH6w1a9bQnXfeSYsWLWogwLh9\/oRZoPum6y+\/3QRrg2+duD3cN+yj4K9tTrTB7o+3ug385J+e9ASv40B2yR2\/42hE4Wgmd16GYadd38vXLyYYX1XAcDrVuse\/L126VPjIsK8LX4F8YPhvN9xwA+Xn54svK+yLw18O2Mp50003efnAcFrV+ml3\/yN\/ZwlEM8adrXnslQ4BY7HP+OUnBYw\/CwybiFncyG\/9nFZdruGX7rx580QevHbeu3dv8Y3N5KW+zP3Vh1\/y4dw3XX8ppKQgCcWU+0IVkaHSm+wLO8tSA9Z9sXkl1R74O6U+vJYQyM5O6sgbBOoJQMCYGwkQMBZZqwJGvkh5p4LcjdC8eXNhTpbf8jmNXELin51+eapWI2mJ8RVUbhcwVtogRZX0hWEnS7ayBGqzEyLM4pCLONnXZ77wnN\/UtMdAOrJwOKVNeYbO7NxAp15bSR2KtxIC2UWMFw+CQFACEDDmBggEjEXWvgJGfUxaNvier6Mv31OFjeozY2r5QtZPvqwD+fBIYebr2xPovqn6q5aXUG3gOqmWF7l1NVCbTbbB4lDTkszXCiP9YWI5kJ0WMMgEBGwmAAFjM2AlewgYi6x9l5CkRYXPRpKOu\/58XXipyEkHUuloqL6o1Rc8N186+vLPbnPileKFRaDVNvjbshqozSYckS0OsYRIhsk9Ibo5oRuJMW6u+yFgLLIOtITEj3NALPlytbKNWk1vsfiIkvnbJjp48GDhiyMDcvnWX31GrWeg+xFVLIyHwm2Dur1dFiO3uTvVhjCaG\/dJMbnHfRcnfAMxxs0NAQgYc6xREggYISAD2jVu1ZHSCl6ii5p+20i5VgrB5G6FEtLEMgGMcXO9BwFjjjVKAgFjBKQjb+3uP3vKbHr9CLEbyckLk7uT9FG2CQIY4yYo15cBAWOONUoCAUcJ8AnbJ8oexi4kR3sBhcc7AQgYcz0MAWOONUoCAWME\/FlgLmragtIKNlByl17G6uFbECZ3x9CjYEMEMMYNgYYFxhxolAQCpghIHxguz8mYL\/7ai8nd1ChAOU4RwBg3Rx4WGHOsURIIJDwBTO4JPwTiHgDGuLkuhoAxxxolgYAjBPhgR75aDn\/EkfLVQjG5O94FqIDNBDDGbQasZA8BY441SgIBRwhAwDiCHYUmKAEIGHMdDwFjjjVKAgFHCMSjgJHRlh0BikJBIAgBCBhzwwMCxhxrlAQCxgj424WU3KWP44HtdEzuq1atory8PCovL6ecnBxjTFEQCFghoGOMWykHaRAHBmMABOKOgBQvjVunewLXyXvcWCej80Y7ufNhqbm5ueLQ1IyMDNq3b1\/c9R8aFNsEoh3jsd16s7WHBcYsb5QGArYTUE+ibtyqg6e8QPdtr5BSQLSTu7S+yCznzJlDhYWFtjahurqalixZQuPHj6doD\/\/kvGbMmEHvvvsulZaWkjwt3dYGENGePXtow4YNNHHiRLuLSvj8ox3jCQ8wDAAQMGHAQlIQiBUCHHX38\/VPellbjj4xkpI7fsfR3UjRTO5sdcnMzBRdwMKF\/\/HFVhi2xth16RQwfCr68uXLadKkSdSkSRO7qtwgXwgYY6gpmjFurpbxURIETHz0I1oBAl6WlkMFvanu048DUkm67HJHgtxFM7mz3wtbYMaMGSOsF7yUVFFRQQ\/+eDhNTz9MqVPWkWpx8td4Lr+kpITeeOMNcYr8\/v37aebMmSKpPLVcPdF8woQJQmz4s8CwGBk7dqw42b19+\/aiTunp6cLCsn79euKT3\/m+FCvS+iL\/Nm3aNGERUZ8\/fvy4p34\/+clPhGAbMmSIyJPzmj59OvGJ96dPnxZl8LOy7pyuuLiYPvjgA7r88stJzX\/UqFGUmpoKC4yBeSKaMW6genFVBARMXHUnGgMC7iYQ6eQurS9saWHnXf6vvJfahGjz\/d+lK\/9jJX1a+nBQIcPlb9myRQgBtkqUlZVRQUGBsIawcLnmmmtox44dXhaSQBYYTs\/1GDFihCevrl270sUXXyzucVlr1qyhefPmeawtqgWGRZHv8wMGDBDPqfXLz8+nxx9\/nL744gv6z\/\/8TyFwfvjDH1LPnj3FkpZsx0MPPURPPvkkjR49WixNcf369etH2dnZQvQwL84Xl70EIh3j9tYqPnOHgInPfkWrQMCVBCKd3KW1hXcdqT4vLEDYKrP69rbUK+WTkBant3cfEFYXKTDYMqFebIXhS1plgllgnnrqKRo4cKAQCyxMeGt3x44daejQoZ57vstFqoBZuXJlg+dvu+02OnnypKgfCye2qPTu3Zv++c9\/inpdccUVtHXrViG6WAAtW7ZM3GdrzyOPPEIrVqwQvjosyPhZKWawhGTu4xDpGDdXw\/gpCQImfvoSLQEB1xOIZHL3ddz1bSRbYBb\/G1HjYXMo+5ONIS0wUsD4WmB885UCgsXESy+91MCJ1y4LzNGjR4WA4YstJ1z2oEGDxO9y9xX\/3KlTJ2FdUS0wUsCwZUa1wKiWJ9cPkhivYCRjPMab7Fj1IWAcQ4+CQSDxCEQyubPvC7+4A10tL6qlCc3+TuvOXEm\/\/LcmlgWMFAjS2sK\/L126VCzRSN8SaYGRfi2yDnyfLR3h+MDws6oFhgWS7\/PsAyMFFqdnccLt\/6\/\/+i9R9M9\/\/nPha8PppPWIl4k6dOhAvNT09NNPe4SW6qPTo0cPuummm+ADY+AjF8kYN1CtuCwCAiYuuxWNAgF3ErBrcucdVmfeWhew0aacltliwpdcopL+Nu7sDdTKDgJ2jXE76hrreULAxHoPov4gEEME7Jzc3RDnxt\/OJFOxXmJoGMR1Ve0c43ENLoLGQcBEAA2PgAAIREbAjsldCpfL8p6gYwuH+90+bsoCExkVPBVPBOwY4\/HER2dbIGB00kReIAACQQlgcscAiXcCGOPmehgCxhxrlAQCCU8Ak3vCD4G4B4Axbq6LIWDMsUZJIGA7Ad9TqJteP8JzoCMfL\/DJwrtFHVreU+Q5UuDkM4\/Syd\/XnyfUdsofqFn2XeLnQPejaYSdk3vt7u10pHggfX3mc09oxsKKAAAgAElEQVQVsXQUTW\/h2UgI2DnGI6lPPD8DARPPvYu2JRwBFh21B\/4uRAv7hvCRAs1\/OJa+fdNYOvzLW6ht\/grB5JOScdT+F69Q3acHPT9X797mOT\/p7MFdfu9f1PTbUTG1a3KX4kwVZqoIU4VZVA3AwyAQgoBdYxzgGxKAgMGoAIE4JiAPcPxWx+50ouxhcf5Ro5TmdKR4ELUY\/BCdPfABnXnvVXHo47maUx6Rc2bnBr\/3k7v0ioqWHZO7tDpxe6T1SK2kv4Mto2oEHgaBIATsGOMA7p8ABAxGBgjEOIG6qvogb0ltvE9kZguMtLqcPXHQY13htCxgmva8WTzna7FpNfoJIWz83fcnEMLBZ8fkHmr7dKi\/h1P\/aNIePHiQ+B+u+CZw6NAhmjp1qjgslCMl47KPAASMD9tNu47TPcvfF3cfG9qZHujXUfz88YkaWvP2Ebr3B2l0eauUgD1SV1dHp06doubNm1NSUlJEPYc8vLHFEw9ume72nPpjMVX\/vYLa\/qzUI2J8rRKqFQICJqKPZVQPsXDhl9r27dujygcPxw4BCBj7+woCRmF87FQtjS79Gy24u6u4O\/UPH1FZ3nepXfNkenPvZzTkv3fQCz+9hvp2vjRgz9TU1NCRI0coLS2NUlICC51gXYs8vOnEEw9ume72nN22lqp+nUdpc8qpSfcc4fsiLS9yyYcFTLwuIYWysIT6u\/3TLIkTpjn0\/4IFC0TYf1zxS4BF6uLFi2GBMdDFEDAKZLa+FDy\/l174WU9qltyY7i\/7gMb3Taf+V7WGgIEY0yI87BAw5\/6xjY7MyRUC5uK0LDq6cESD84BUUcN1iCcnXumsXPfpxwGnTKd3I9mxdGbg\/YAiIiCAvo4AWoSPQMD4CJjlbx6k34zuLu6ygMnt0kosI0kLzKJ7utLw69oFxF1bU0tHjh6l1HbtKKVJZBYY5OGNN554cMt0tydp\/\/9R1aM\/FALm8w1LG5wJJHfmqNuo\/W2XvqhpC0or2EDSaiO3Ufvej3CuEY8l6uSeqO2OZqzE6rPoa3M9F5GA+evhL+mOkvfojqsvo4XDu5H8\/eTps5TVtqmwYPCyS6xdbIEJJWDG9+1AP7khPbCAqa2lqqoqatOmDSUnR8agFnl48Y0nHkLA6O7f\/f9HjRbd4llCcvPnLlEn90Rtt5vHol11Q1\/bRbZhvhEJmCnPfEhZbS8RlonT\/6rzWmpZuuUA7fnktBA2sXZZWUL64wM9qFdmi4BNq6muoaPHjlFaaiolp0QmYJCHN9544sEt090eXkI6WXwTBIyLJxy81FzcOZqrhr7WDDRIdmELGHZ0nbT2Q1oyspuwsrD1RXV29f27uaZEXxKceH2EAxySvYDocL4VAkYzV9UHhp143Xwl6uSeqO1281i0q27oa7vIarDA+AoUtriU7z4h\/EYuuTiJYlnAMB51G\/Xvx18tHHj5wi6kyPx5dL+sY31nFwRM\/W6ccLaY5uXlUWVlfaybYFdGRgaVlpaGShb236urq2nJkiU0fvx4atWqVdjP8wN4qUWELSYfQl+b67awLTBcNbmENOzatjTk1+9R3vXtPfFS+G98xeISUjDsEDAQMDrEGARM+AJm1apVxCIm1MXi5e5uzT3nPanpo9mFBAETijz+rhKAgDE3HiISMGxlYeGy55MzlHd9mhAr8l7fK1vEnXiBBQYxbXQJD135qGIqEZaQcnNzqaKigsaMGUOjR4\/2zJBFRUXiPltf9r77pjj7iSMJB4oYzC+XkpISeuONN4QVaP\/+\/TRz5kyR39y5c2nEiBE0f\/58WrZsmbg3YcIEmjRpkl8LzIkTJ2js2LG0c+dOat++vbD+pKen04wZM2j9+vU0ePBgcZ+f5zThWp7MvQZQkk4CEDA6aQbPKyIBY6569pWkLhVxKTLqbqBIvLDAwAIDC0z0n8dIJ3cWKSxi+Nq3b58QLLyslJmZKe6Vl5fTDVf7j4Hj++14y5YtNH36dNqzZw+VlZVRQUEBNWnSRAiXa665hnbs2CFEB9\/jK5AFhtNzPVj0yLy6du1KF198sbjHbV2zZg3NmzcvKgHDFiius7xYwOXkuNvXKfqRErs5RDrGY7fFztU8KgEjdyBt3HXC04IBV7Xy+MM416zQJfvbLQUnXm9uOl7YyKPhWNTNJBEsMExRtcKwtYOXlfjlzlYZ6fuinsbtbxbglwtbXaTAYKuIerEVhi9plQlmgXnqqado4MCBlJWVRWyNYWtQx44daejQoZ57y5cvj8oCw+3j+hYWFnqqyeWweIvU34frKuslRVroGRMprBKAgLFKKvp0EQsYKV7SL032LBnJe1wt6dQbfRXtyUHdCi5LsLKNGkcJhNcful\/WcOKt5x9sCUkNWCeD2IXXa\/aljmZyVy0uc+bMIf7HFhC2vvB\/g0XklT4wb+8+4BEwvhYY31az5aW4uFiInZdeeqmBE6\/dFhi2OrHlRRUvso4sYvr16xeRJQYCxr7xzTlHM8btrVn85R6xgAm02ygWdiGpPjzcpTL43l8PfhkykB0i8Yb3IdAddTbWoxszPd1M1Ei8fJSAPAeJy5JHBjRu5Y7zd6Kd3PnFzcJFXvyzvxd8oFGqWmA4zbp16zzWFv596dKlwkeGfVb4khYY6dci8+X7vCvJTh8YtjAFs7KE+rvKgMWYP9+cF154oYEPkOwjfr5Hjx60cuVK2r17t\/DhkUzkEhzX4fDhw+I++xRxepUVW7SGDBnidY\/TdenSRVis3n33XfE8+wvxUls8WISiHePhzbCJnTpiAcPY1Mi1vIWaL3+WDbcjllvBf5zdgX677VDQowQQiTe83tQedTbGoxsLAWNjJN6vTh33e2hjIKfW8Hoz+tQ6Jnf2e2FrDFtd2B9GveTxB741jWYXUjitZkHEl1yikv42kbQ7lEDhJTW2Plm5fOvFvjnjxo0TAk71AWKrzjvvvONZGuO85fIYC0XeRi7z6tatmxAiLDq4fdzWO++8kzZs2EATJ070VIvTc3+x6JF5Pfjgg8K6NWvWLOH4zD+zbw8vx8X6FUlfx3qbnap\/xALG14rhrwFuOVZAPerAn4+OXDqadnMmPf7qvqCHOSISb3hDVXfU2ViPbsz0dDM5++J8OvXsoyISLwuYz9c\/SWkFL4mOOlI8iJr2vJlaDn8kvI6zKbWOyV1uq2brBPu\/yEs9sPLkn56kZr2Hi91IR58YSckdv2OEgb+dSfxSjqTdwQQMCwL+u1UBw8tdLE6ys7OFiGAfGP7Zd3s6W0yuvfZacZ8tI2xRmT17Nk2ePNljaWHe\/ixTfI9FirqLS\/oUderUSZQnl+WGDRtGmzZt8jhLq\/5ENg09Y9lG0tfGKhdnBUUsYGKZg+8yl4xdM+2WTBpd+jdacHdX0Tw1wjB2IWEXkg5\/HiFgNEfiPbttLVX9Oi9hBAwz5OUH36UjFjDyJO4vNq+k2gN\/p9SH1wrfGH8ndJucwyJ5qen0gVEtMNL3h61EqgXGHw\/+++nTp8WSmrTAyHSq6JAWGBYw8pLWFhZBfD6crwXm2WefhYAxOQjjsKyEFDDcj+p2adVShEi834xy3S\/aWHfA1cHDDgGj7kJKhCWkQPPw12e+8FicmvYYSEcWDqe0Kc\/QmZ0b6NRrK6lD8VZyyhcoEgHD7WRLCC+VSbHGlhfe+s1XOL4\/qg8MC4ru3buLpSPVB4bzZP8UXgqScXD8+cBwOt\/dWjLmzaBBgyg\/P99jrQnmA6PuhIIFJg7VhYEmRSxg2Irx+Cv76NGhncURAvLinUiPPL+X2JoRiydSB2IOCwwsMLEgYOLdiTfUnOhrhTn5+\/rtx22n\/CFgcLtQeer4e6QChsuWlhh5nAIvBalLZzrqhzz0EYimr\/XVIjFyiljAMB7pW3LH1ZeJrdS8FPPc+5\/Sc\/k96Xvtm8UVQQgYCJhYEDB8mKO6jdrpF7fvJJCok3uitjuuXgIWG4O+tghKQ7KoBIwsXy67qIcfaqibq7KAgIGAiRUB46oPjk9lEnVyT9R2u3ks2lU39LVdZBvmG5WAiSULjG\/kXXUXlbozCUcJfDNIdLywkUfDD51uJokSiTfQtCiXjS7Le4KOLRxOdZ9+3CCpqW3UgeqIl5q5l5rTJaGvzfVAxAImlnxgeGmr9K0jnoMnGa+MV\/Oj69Po\/rIPaHzfdPpeejPsQlLGnu4XLZx46+Hq5proAsbcdBl5SXipRc4u1p5EX5vrsYgFDFdRPQuJrRjzh3WhB1b\/XWxDdosPDFteOrdtSi\/\/rUpQVU\/OLh7amfpf1ZqkdebW77ahguf3Bo0Dg0i84Q1O3VFnEYm3nr\/KVY3Eyz4wbr7smNyDHSEgWcAC4+ZREV91s2OMxxchfa2JWMCoZyGN6ZNOSzbtp0Uju9Lv3jpC5btPuO4sJBnrRQoYGe+FhRYi8fofUNojxsZ4FF0dPIT4sDESb6wIGA6M1qtXL30z2YWcmq\/\/BX19WSadvv5+T94X7y6n5I9eo1ODf6m9PKsZHjp0iKZOnSq2KXNAN1zxSwACxlzfRixg1GBwn5w66xEwX9Z+RZPWfkhLRnZz1TZqXQIGkXjDG5y6o84iEm89f5UrLyGdLL5JBLJzu4A5ePCgeJFv3749vIFkIXWbi7+igoxPqLiyLVX9q7HniUD3LWSpPQkEjHakrssQAsZcl0QsYLiKLAoOflZLD\/XvRCvfOEQFg66gESvep75XtvCcUG2uKd5LWi0v+ZbXdm5fATPk1+9RJEtIOI06vB7V7e8R6340QnxojsQbSz4w3H4WMfzPjostME12b6LTvccLK4z8vbpLf0ctMCzYFi9eDAuMHZ3usjwhYMx1SFQChqup7trh3x8b2pke6NfRXAsslqQKGCm+stpeQnDiDQxQ94s21sWHDh4QMBY\/sFEkU+PgXNS0BaUVbKDkLvqXq8KpIl5q4dCK7bToa3P9F7WAMVfV6EryFTDqNuq869M8FiMcJfANZx0vbOTRcNzqZhLMAhPoZS7D7tfu\/jOpDq6B7quOssld+ojDIi9q+u3oPpQJ9DReaonT2ehrc32dMAImWqQIZIdAdjqEh0kLjHo6M1sgTj7zqOdMIPWwQ\/W0Zk4jD0FU78ufW9z2H+K8oRaDH3I0NH+0n2fTz+OlZpq4c+Whr82xh4CxyFoKmILraqnPlZcGfKquro42b95M3bp1Iz5CPpILeXhTiyce3DLd7bnu0rN00YofhXTird29nT4pGUftpvyBqkrGUtOeN1PL4Y+Iowc+X\/8ktf2PdfTJohEN7rfJX0lHi2+mVqOfEKJFFTmRjG87n1EtSGwpapu\/nI4uvIfa5q9wdBkJLzU7e91deaOvzfUHBIxF1h+fqKHJv3qOdjzzBF105lOLTyEZCNhP4MdpJ+knA3pS25+VUlKbjIAFSuHR9oH\/52VFYQFzouxhSi14VQgbaV3x3J\/6Rzr2qzEeEcD5nHnvVdctI0nx0rh1OrW8\/SE68UwRtXtoLX3+4iLH64uXmv2fA7eUgL421xMQMGGwtnP3RBjVQFIQ8CKQmlxH6enpQcWLtLKw7wpf6jJQvAgY9STquk8PegTMuZpTdHThCEqdso4at+rgyOjBS80R7I4Uir42hx0CxhxrlAQCthNgX5Uzb60T5ciTqH2XfKSlIh6XkLj9Xx0\/SK2GF9JnLz5Jl41eSEeKb6GLs66n1IfX2s4\/UAF4qTmG3njB6GtzyCFgzLFGSSBgnACLF77Y10W9AjnrxoMTr7rzSrT9nqIG7TfdEZG+1PLy8igjo+GyYGVlJZWWltrajBMnTtDy5ctp0qRJ1KRJE09Z69atE\/59iCjsH3+kfW1rZ8Zp5hAwcdqxaBYIsNPukeKB9PWZzz0w5JbpRinNxTISb6NWt0X7OsHK7dLqNuqm149w1JoRiz0b6UutqKiICgsLGzQ50H0TbCBgglOOtK9N9F28lQEBE289ivaAAAi4jkCkLzVdAobLLykpoTfeeENEA96\/fz\/NnDlTcJo7dy6NGDGC5s+fT8uWLRP3JkyYQOPHj\/dYYHbu3EmjRo2i9u3bU05ODt1+++0iD2mJ2bNnD23YsIHGjh1LM2bMoPXr14u0bCXKyspyXX\/YWaFI+9rOOsVr3hAw8dqzaBcIgIBrCET6UtMpYLZs2ULTp08nFhtlZWVUUFAgloZYuFxzzTW0Y8cOr+UiuYR033330eOPPy4sQZyeBcq9997rV8C0adNGMGdBxM\/L+rdq1co1fWF3RSLta7vrFY\/5Q8A41KtshuV1bJ5Q+FK\/\/chvRHyf08lvSvytyF96p77pqHXr0aMHrVy5kniikh9g9dsd\/xzufZNdE25buG7V1dWeyVz6A6j96FS\/mOTmhrLU5a1A9VGjDTtR50hfajoFDFtMWFion0PJguccvtS5RlpgBg0aJOYhKXjkEpI\/C8zp06c9VhzOLxE\/A5H2tRPjMtbLhIBxoAflS04KEh7wfI8FgPoNp3Xr1jRt2jTx7Yd\/ZvMsCxgWC\/JbkFOOdPwt7rHHHqOFCxcK0cL1P3z4sKjv7NmzadasWYKsTMM\/T5kyxfJ9k9\/Ywm0L142\/XXJ\/sGldnjDsT9A4MLwSrkhfAeNGH51IX2p2CBhfC4zvgOFxXFxcTMOGDaNNmzaRaoGRn\/V+\/foJC4y0tnD72MIjHY5ZKCXqFWlfJyqvaNoNARMNvQield9e+MPOFwsSf9YYnggGDBjgeemzgJECQP3ZLevL\/KFds2YN3XnnnbRo0aIGYozb6k+kBbrvlDCTlqJgbWEBuWTJErr77rvp5z\/\/uehDri+LGtlHbumXCIZoTD+iOiHLhjhtfZFjin1IpNi1CjnQLqSKigoqLy+3mo2wukgLDD+kWhz596VLlwofGRbkfAXygeG\/3XDDDZSfny++VHH9+IsLfyZuuukmLx8YTqtaZi1XNsYTQsCY60AIGHOsvUril7kUMP4sMGx65Rej\/KbPaf0t0cjJRi4tOdQcIU5YdLFTH7\/8582bJ6rClqLevXuHfd\/Jb3Ch2iLrJvtGChhf07y65OdUv6BcIo4N8689b1GH4q0IZIcBYTsBCBjbEXsKgIAxxzqggOE\/yGUl6eXfvHlzYjOttFpwGrmEpFon5LIFiwSnXvqqBUlaYmJVwFhpSyABo3awG\/rFoaHtaLH+\/GHccHo2XmqODgujhaOvzeG2XcCoTo1qsxL926lqgfHtbmkB4Pu+jr7SaqM+Eywvu4eSrKt8qQfy55EizdfPJ9B9J5aQrLZF1s3XAuOvH\/31l919kmj5+4oWGYHYTRzwUnNTb9hbF\/S1vXzV3G0TMPIbKBfG38bVSI7B\/mau6c6W5LuEJJdd+Lwl6bh7\/PjxBn4jcjlGbokM9RK1s5XcBrYSqWJD9QPhsmPBiVeKKKttkQ7G\/paQ3NAvdva5G\/OO511IbuSNOgUnAAFjboTYJmBCxQAI9XdzCJwpyddqolqqVEc\/K9uonbBm+duKOXjwYCFWZdArJqu2RX3Gyn1TPRNJW7hu\/sSjbzAwp32TTDFEOXipYQzUE4CAMTcSbBMw8lstR3ZU45rwfflSduLFaw4tSgIBEAABvNQSbQxAwJjrcVsFDDeDYw7IrXayWYkY3Mhcl6IkEAABtxHAS81tPWJffdDX9rH1zdl2AWOuKSgJBEAABNxJAC81d\/aLHbVCX9tB1X+eEDDmWKMkEACBBCWAl1ridDz62lxf2y5g1JDr\/pqF5SRznY2SQAAEnCGAl5oz3J0oFX1tjrrtAoab4hsqX73HW1flFmJ1q7U5BCgJBEAABOwlgJeavXzdlDv62lxv2C5gAm2XlvcnTZokzpXho9pNHuBnDjFKAgEQSHQCeKklzghAX5vra9sFjAxa9+6771JpaSnxIXdyZ9J1111HQ4cOpeeff75BsDtzCFASCIAACNhLAC81e\/m6KXf0tbnesF3AyKb4BgvjQGZdunTB6b3m+holgQAIOEQALzWHwDtQLPraHHRjAsZck1ASCIAACLiLAF5q7uoPO2uDvraTrnfejgiYuro64di7YcMGuvjii+lHP\/oR3XjjjdSoUSNzLUdJIAACIGCIAF5qhkC7oBj0tblOcETAsA\/MsWPHqG\/fvvT1118LH5jrr7+eUlNTzbUcJYEACICAIQJueqnxBorly5cTb6Dwt\/OT\/\/7ee+\/RD3\/4w6B0uE3y8FJDGGOiGDf1dUwAi6KSjggYPnH5\/fffp5tvvpm++uoreu6552jQoEHUvHnzKJqCR0EABEDAnQTc9FILJWC4rvv376cRI0YEhInz7AKPMzf1tTs\/DfpqZVzAvPLKK7Rr1y6vFrRo0YLuueceatasmb6WIScQAAEQcAmBaF9qubm5VFlZSRkZGVReXh52q+Ru0PXr1xOfGs8BRGUICz5wly++P23aNJo4caI4UZ43WrCFRf07nzbPebGFpmnTprDA+OmJaPs67M5N4AeMC5gEZo2mgwAIJCiBaF9qmZmZHgGzb9++sCmyxYQvtqpwXTh4KIuVmpoauvLKK4UoKS4uptGjR9Px48eFBWbAgAF08uTJBn\/nUBh8YQnJfzdE29dhd24CP2CrgJHxXljZT58+3etk6h49etDKlSsRvC6BBx+aDgKJQiDal5q0wDCvSATM\/PnziaOeZ2dnk7qExJaWUaNGiW6Qx7pIASPFju\/fIWCCj9po+zpRPhM62mmrgOEPDZs8+YMgTZj33nuv+BD5O15AR4OQBwiAAAi4jYDTLzXVAsNfLMvKysS8zEtEvGTkzwJz7bXXip2ivn+HgIGAccvnyzYB43uEAH9oHnvsMVq4cKGwugQ6YsAtYFAPEAABENBFwGkBo\/rAsPW7e\/fulJ+f7\/F3YetLTk4O3X777aLJbHWZPXs2sc8MW2nUv\/MXUCwhBR4ZTve1rjEbC\/kYEzD8DWDr1q2eIwMgYGJheKCOIAACOgjgpaaDYmzkgb4210+2CRhuglxCYmewsWPH0siRIz1b8\/hvfLFvDC4QAAEQiGcCeKnFc+96tw19ba6vbRUwbGVh4cImyAkTJgixIu+xGRLixVxHoyQQAAHnCOCl5hx70yWjr80Rt1XAmGsGSgIBEAAB9xLAS829faO7Zuhr3UQD5wcBY441SgIBEEhQAnipJU7Ho6\/N9bURAaN6wMumcWwYjuro7ywOc81HSSAAAiBgPwG81Oxn7JYS0NfmesJ2ASPFC2\/Dkz4v8h43EyLGXGejJBAAAWcI4KXmDHcnSkVfm6Nuu4AJtF3a6W3U8jAyRq1GBZaDj+\/PnTs36IFm5roJJYEACMQyAbzUYrn3wqs7+jo8XtGktl3AcOXk2RuqtUWN0htNAyJ51jeoHtfl8OHD4mwQDt40a9Yska0aeC+ScnQ\/U0mVlEEZurNFfjYQ4BPX+R8uEGAChw4doqlTp9LkyZOpV69egBLHBIL1dXp6OvE\/XHoI2C5g1K3Ugars9LlIUmDdeeedtGjRInFGE\/vmzJgxg+TRB\/wyevbZZ+muu+7yGoB1dXV06tQpat68OSUlJQXsFavpOINAaYuoiDpRJxpDY0Q5VvO0mi5W8gynPXo+JuHlwmOFX1bbt28P70GkBgEQiGsCLF4XLFgAEaOpl20XMJrqaWs20hrUqVMncUorW4r4YgHTu3dvzwmuHF77N7\/5jTjLSV58muvRo0epXbt2QR2SrabjfP2lPZh0kLKSsqigroByKIf61vX1m84fqGjLdlueantCCUdbB06AzKUJmSeqDh06OFEFlAkCIOAyAvyFZvHixbR69Wqvd4jLqhlT1Ul4AaMeKqkudQUSME888YTwmZFXbW0tVVVVUZs2bSg5OTlg51tNxxkES7stZRuNShtFC6oW0O3Hbzdatto4Xe2JJE+1bBaOLVu2dNWHDmvgruoOVAYEXEEA84L+brBdwKhHt6tbpnkn0pIlS2j8+PHicEcnLl8\/HB5gfM\/fEpIcfL\/97W\/p+9\/\/vqe63I5jx45RamoqpaSkBGyG1XScgZr2aMpRejrpaUqvS\/csHclCrOZpNZ1v2Xa0R0eeanuaNWsWdOnOiXGFicoJ6igTBNxNAPOC\/v6xXcBwldlpNi8vjzj2C2+lZpHAp5yWlpaSPJpdf9OC58h16Nevn5cpj8XWlClT\/DrxBhp8vJxx5MgRSktLCypgrKRjHxd21C2pKaH\/Pfm\/9Graq1RKpZRLueK\/vg68VvJkClbThZPWyTzDKdv0uOLyMFE5QR1lgoC7CWBe0N8\/RgSMrLZbOlDdKi3rJgPr8blN7OvCl7pWaULAZFImFVIhjawZSW9XvS2WpbqldBOixt\/uI6svcqvpIGD0fMDcMs71tAa5gAAI6CCAeUEHRe88jAgYN1pgwkVpQsDIOoUSHGypqaAKernmZW3WH6tlh5vODlEUik+4fas7PSYq3USRHwjEPgHMC\/r70HYB42YfGMapOvGq5n\/+WQ1kZ6eAYTHCoqScykUPW31B605nR9l25Gm13fo\/LtZyxERljRNSgUAiEcC8oL+3bRcw+qusL0f2g1m2bBlNmDBB+OY46QPDIoa3Rzv9wrcqDqymc7o9+kaL9ZwwUVlnhZQgkCgEMC\/o72kjAkY9zJF9TdSIt0458bLlheO+bNmyRVBlAWNlF5LvHn6rL3J\/6VTRonatlTzzKI\/O1p2lnid70sTmE6N2IJblq2WvXZtC+\/cTFRYSVVYSZShBgK3UUea5d+9eEVPngQceoP79+wccxVbztJpO\/8fFWo6YqKxxQioQSCQCmBf097btAkY9zJEj3ZaVlVFBQQG98MILtHXrVscPc2QrjCpgTAWy48B045PGix1H\/nYXhQqOx1ur+fnvf\/l9evXsq9qD6NXVpdP3BjcjKiTqW1dHb3ZOooL0OpqdXuexEoWqoyqKBgwYQK++Wl\/PN998U\/ypb9++XiPaasC9WAlkh4BV+ics5AgCsUoAAkZ\/z9kuYNRDG48fP+4RMCxsioqKqLCw0IPTzx8AAB1\/SURBVLE4MIwzXAGjK5CdDEj3+oHXRYwX9bIaJE53Oq6DzHPduquoT5+zlJ1dQweTkujGH3eku1p\/SQsWVImqWi2bj2L46quviLn94x\/\/oD\/96U\/ieQ63z2Ix2nYjkJ3+SSHRcgzkpyettGrkbatseOPChg0baOLEiV6P8HzDLzKONSXjX3E5a9euFfd4juSQE3w2G19yeTtUufLIFk4n8w7UrlB54e\/2EICA0c\/VdgEjRQJ\/IO+\/\/3565plnKD8\/X3yweWLgpRsnL18B47ZAdqECv13d\/GrqR\/1oRd2KgBitBrLLOphEfTbVUcs7PqOlNam0sa6O+l6wuHDmleednq+6NYV2ldRQu3ahA\/hVVlaKQIVscWGHaF5Ceu6552jkyJF+62q1nghk5+QnJnHKtkPAPPXUU7R582Yx7\/H8x2O5uLhYnKfGS+slJSU0evRoT3wsf\/Gq\/PWAP8EEAeOusQoBo78\/jAgYrrZv7BV1h4\/+ZlnPURUwppx416asFctG0mnXt7ZWfTw43YzqGbS45WI6R+cCNjpYfixKisT\/EfUrJdo7u47Gjj1AWzIzLxwZWZ8tB9PLodE05\/zdORVEHfbW0U03HQgawI8tbBUVFeLQyZ\/+9Kc0bNgw4avD91etWkX79u3zqnM47bYSPND6KNCbEhOVXp7h5MbsWQS88cYbIo7T\/v37aebMmSILOedI532+xxYOFtnLly+nSZMmkYwD1b59e8rJyaHbb79d5MH+ciw4pFAYO3as8OvigJyc1jcoZyALDAuYSy65hJo2bSrOWON0f\/zjH0X9ZD18A2z6tl+dS6VPIX8h5LqrFhtVwPCSveSgblrgdsg2czRxX4tROOyRNjgBzAv6R4gxAaO\/6npyVAWMr9CyK5DdVSlXiWMBOGidv0v3izxYfo3OW1b4bOvK8wKmvND\/Fm4OpMdOw3Kb96pVRL\/4RR2tXn2EfvCDNkEdiLl93o7Ba4WJnK9z57xFl+526xkh4ediZaJip2j+l1O\/8YwqKur\/i9+t8QjUK8yeHfPZwsHiQPrcse8Vf9avueYa2rFjhxAr8mgT+aK\/77776PHHHxfL2upp9P4EDAeZ5ItFiLpMLpeFggkYPoqE68h1YGHBeb399tueY1V8BZZqpVa\/ZPEGCLYS8XXttdc2WLKS7Ro0aJBYwl24cKFYtpJWHa6DFEu+4STCH\/V4IhQBK\/NCqDzwd28CRgSMv8i3XA0+FFFdC3Zz59gRB6aIrR79vnlp8e8Vo4mojKhiDp+JVKNld5H\/HVB0YdO2N3U17dGj26isIo86jSkUgoutJvxCKC8vF6Lk7ber6LXXOopdSr6XamXxVz4vL2Wo25psiH\/j1HiyMlFxX8+ZwyKuvpaNWEkKUYffrfAIJmBYcLCw8DfvsBWGL9UaIS0f\/KLnFzlvMmABI5eQ\/AmY06dPixAM8vK1wgQTMAMHDhRi44YbbqC\/\/OUvNHToUGEB8ncunBQo3B6+fJeFpGDjDRK+PjeqgGGxIq0rnOe\/\/vUv4WcjywxUX6c+Q\/FYrpV5IR7bbWebbBcw0rlMrvna2RgdeauTnu5Adnw8wOk2p+mWlFtEVVdVEJWx5aNcrODQqlyinEKidtk1VDOjmhbvbEn7Sr\/ZvryKiPaff8lJvaAKg\/yUfHFekr9LpqtNS6P5KSni+czz5WVwhkVE6krON3nWUsq2blRZlEkZpeWiErztO7dRLtGYMbSrpISq3q6iDRvSqH\/ng+LvbD1gYcKXtLJIsaMu+XCaTF6iGjNGmN7lBQsMLDBWLVBWBIyvBcb3Gel\/wkubmzZtItUCo1oqWMDwJUURiwEpvqWw8M07lIBh4fLXv\/6VevfuTbm5ueJg27vuuot+9atfeW1s4PlICjIpYOR5bbDA6JjxzeUBAaOftREB44bdRlbQ2ekDs7duL2UlZQlLRmUuqxKijNILYqSSKLOSaF+OuC2sG0+dOkVTT7ehfZzuggMtt6HsgoBhAcLXa\/v20UdpH9GtKbdSBp0j6VXCgogvFit76+rEkQMfnG5D61JSqPRCnhdWckjREKLskzv+l9I++n9E+0cTFWYQVfI6x\/n1jooKyjs\/ea+qqBACJrlrV0r76CPadut8qswZTTmlY6jsvCKbc960wP\/YFM8X53nq1FPUsuVn50+O\/s8LrWGhwy37JrhMIgkYK+MRacIn4PvCZ2uDtLZwbkuXLhU+Muz3wVcgHxj+G1tIeMNB69atPTuD2Gp80003nfcT+8YHhtP6WpPl8SlyNxGn4SXpd955h9gCwxc77fKSVXp6uhAwbA3x3YUkz2iTy138nK8PzLx588SuvkAWGLlUBR+Y8MeTzicgYHTSrM\/LdgHDhcTK+qqpQHZHj6ZQ2fn395wcEoKDX+HsAyG\/fcoX+UdpaXTLeadXXl1gCVF\/0ED9xQLlbN15h9sD3zjS8j1pnWHjCnua8IrEpro6ml1XRzlzuC9ShMWHV2+8gtNVriKq2CL+sK\/Tj6nRlkcpYwxbR3hdq4JycyupfM4Yqhw9mspyc6lPXR11Xb2V2pxuQ0e7pYj8KnKLqKiikipolfiGyk66vFTS52wd5YzJpaTOb1Jl3jnK6MeaKJfKynJo9Oh+lJHBsqw07CMUNqRtEIddBnKG1v9xsZYjJiprnJDKWQLqDiffpSpnaxafpWNe0N+vxgSM+i1INsNtPjA8wIIFsuN1Zl6vfvLJJ0UTOHYD77DhgGwcj+T5558X9+U2Yf47X\/vH7KfLvryM5l2VTP37JxGVjqGSmhrxvL\/0d9xxB3GQON7BM38+i4Mx9PLLDdM\/\/fTT9Nlnn4lvg\/wNLVB+6en30YABSfTv\/76CsrKSzi\/zjKFdu2rogQfqy0+\/L51yMvPo6V8eJNqSQ9nTSzzt6duXY9RUUFZWsfjm+sADx6gisxP1LbiP6lbUUdK3kojKKykp6U2qyMygvMo8GjOnkCrm5NDojEoasy+XZbIQQiOnT6d5646KuDerCnnXxdO0b04GpResOM9x9oVhMYbq6lZQUlLngKOdBd5nn71HH9a8QremzKdzGYF3YOn\/yITOERNVaEZI4TwB1UrktrnYeTr6a4B5QT9T2wWMr9e8\/iboyzGUgLn11luFYFiw4BO6664v6YorrhCFv\/7662InwVVXXSV+\/+c\/\/yn+K\/\/Ov\/\/jq3\/QgKwB9YLh7Fl6\/cABr7+r6Xft2kVVVVV04403ivSrVx+m9PQ6uvHGy8XvvuWXlu6jK286QDde\/k16TjdqVHuR\/ne\/O0DpP\/rR+Y3Q9VtdPn79daqkDE9+Z\/f0paq3f0rtR41q0B42b7Ppm03U\/fsfoudbLqb\/ebM\/VeSUUtHr++n+jxtT8sC1lJR0UFhbmOG2bdupJrua6u7bRCeHtqDNP\/4NbRndjz67oyd1OteJ5t16K121dh5RZoawQL3x+gFhRm+0fz\/9q\/2r9PXXP6bk5GRRF45YPKt2FnVIP0sl1dm0ty6daj\/KpvS+P6LczErKqMyhwvJCseXVLRcmKrf0BOoBAu4hgHlBf1\/YLmDYUU6u78othvqboSfHUEtIN2yuo7o+dVSwqS\/9z\/+wRaOOpkw5RceOHSOOoRAs6NyAN5OEK8lfv\/wy5LZjfwHd3nwzSZTJSzW8e4V3KHG6Vq1aisbz7xzdNzeXDSY5dO7ceWfbykqqKCulHpN60NeP\/JNa\/+k6opw84T1cU10iFqZSbm1yfvVmDtWk1gcU9C2brTy8Ns+iSt01lJ\/fhHg79cmTJ4W1iMUbP+t7HTyYRFmdk+hnn35KnZOepsWXLqZdNbuoScpVVN4og3IyKqlRZSUt\/\/flNO7P4+jDlz9UWFZSZcUcOnrrA5S3L1f47uTklrNBiBqd43aUEuXlCkdgdgh2y4WJyi09gXqAgHsIYF7Q3xe2CxiushqbQX8T9OVo1Yk3NTVb+I\/wxY6wP\/tZFU2c2LyBMFlFq4RVorJiDFW\/XCP8WNgGEkzocJ6hnFml74rcyszWn27dUkR9Gp33mJlTWUqFGfxCr\/eK4XRJSVnnl3n4IMVSospMorxSovKc+gaEvROokjZtmk0Dxv+Zqnftom1Hj1JuXh6VFhbSGMUSciE+HnWo20tjDzSmzR03U\/+k\/iKIn9jRxDCE4iIqLS+l+36xkk5NuY5a3rGYKHcfVZRXUl5mrmDNrRm9r5xyOH5KWRld0r07\/duLP6fC0RU0Zox7rC9yvI8aNUpYriIJQ69vRCMnEAABtxCAgNHfE7YLGLmNWnr9q01w47qr6uHvL5DdgdcPUN\/0jVR+IYoKv39\/85t9fiPSZlayFaRQeOfuqqkRO4HS0tKiFjCSoT+hk0mZIuCcekAkp6uqepvatPmBKLuichUVUREVZpQ2cIANJZ44oB1ffAjlfXv30ve++13hOFxUWUmFGRlie7YIjHehkqxRciiPuld19zo1m4PjcV3H0BwRCnjVnFV0dvlyOtK1K51u8yJ1m7+dKnkXU0aOEEesYnj7NQukfy8ooAFZWSIejZuWjmS\/YKLSP1EhRxCIdQKYF\/T3oO0CRn+VnclRDr59\/2SvjTnnQ\/d\/s0VYFSZyNxHvxllVlEH7yut3+6SmmhEw\/uh4RcKN8hgDYVWiIrEM9HZVFf2gTRs6mpLi2QxdL1jqLU1ykzRbhVYfWU3DWtYfJSAvwaiskgr7ZVBZRp54kLeb76nbQ0e31TsYZ+SMIbbRyDyZel4FUWVuIwgYZz4KKBUEQCACAhAwEUAL8YjtAsbtFhjfLd6hAtnNXb2aqrOzxQuVX6zSsnLFFZn1AekqV1FFThFlZLCvRoaI7cJXKMuG2k9W0wZKx4KBrTBye7GajuPFcMA71UKjs2w1L7aycDlBz2LiGDgcl4YdasqK6LHHHqNh19QLHTXyr3A\/rqik3NxMyskpp+rqbLH1PNRynP6PTOgcMVGFZoQUIJBoBDAv6O9x2wVMoCpbPWVVf5O\/yVGeOaIebiajXHIqfpnK80P8DT5eJkmtqaH\/3XFSxEPp3DlJvIznKJX2FzU31Es3WgFTvzTzzVlLwlem6m36wYUlpGBMoy1bzZutNbzk9HLNy9T1SNeAy2eZuWKliCpWZdJd111H48aNo2633EJ5uez\/Uuk59JH9iVZxuGKqpNdee83Scpyd4ydQ3pionKCOMkHA3QQwL+jvH8cEjL8D0PQ3L3CO8pwTDgvOFx91EGoXkj\/HTA4SNyAp6by1o0IceVhK+3ziy5q3wMhWSwuIcJilXLE00zlIfJVwLEWc36MHHg0qirh8Xm5if5lg\/j+ZmeJ0AurXr4JmLF1K29eto9Jz50TUYXHQwAWHGnl0Eod879mzpzjbJZQYNDmmZFmYqJygjjJBwN0EMC\/o7x\/HBIwqFpzcXq2eRh0qDgwLmMa\/bUw12TVCDMgX\/juffkpLU2fQsZRjtLFuo99eYssGbzfmgHdqWHB\/ia2mDZWOT71m60dqTSr9+stf008v+amWsjk2Cx+LUFpZSkNbDPXkyff5YsFUnFQs\/GQko1BtF9uts4opNftDSvrdL6l\/cTG1KykR1qxx45Kob986IXJ882vevPn54wmS9H8yosgRE1UU8PAoCMQpAcwL+jvWdgETzAfGDdtMwxUwfCw9756S17om62h37W6aVj3NE3zNXzfV1taK4HS85VkGaQvUnVbThkonBQWnSz6SbHvZ96bdS9nV2XTXl3eJpnHEXb5C1ZPTsID5xz++oiuvPCDqyWfW9OrVy7MNmf\/Owfx882NB2LJlfSwct1xyopo8ebJoAy4QAAEQOHToEE2dOhXhFTQOBdsFjMa6RpyVGjLb93A0XwHDv69cuVJYFWbMmEH33nuveInKl9Jvf\/tb+vL7XwrnVP5fk5QmNPnkZJrTaE7Q5Qx\/wekCNchqWt3puD7R5MmCSYoWtW2R5sl+vRy8b+PGeuEiLzW\/Zs2auc4Cw1GLeaLavn17xGMWD4IACMQfAf5Cs2DBAhF5HFf0BBJCwATDpAoYK4Hs5q6eS0uzlxL7d\/BunuyabEvxXaw6x3JdrabVnS6cstk5d1PdJrFkxktVzIIdh\/1dVuvJvjDXXfclPf10khCD\/DvHxVNPyw6njtF\/PCLPgUUM\/8MFAiAAApIACxeIF33jAQJm\/nxBk514+QoVyM532cvqy9lqunBe0FbztJounLI\/rPmQdpzcIWK78BEGwU6Etlr+7Nl19OGHNULA8Ind0mnXd7hbzU\/fxwQ5gQAIgAAIuI2ArQJGLt3wsg0LhFg+\/TSQXwOfRr1582bq1q0bderUKWD\/Wk3HGVhNqzudHWWHm2dZGe8K60eVlRn0Sz4d28+ltrtPnz74RuO2WQX1AQEQAAEDBGwVMLw8wwcAjhgxQvhWqD4lvgHkDLQ1qiLg1xAVvrAerqnpRUlJh8QJ16EudpTlf7hAAARAAAQSi4BtAsY3zgtbX9TAcE7HgYmkm+HXEAk1e5\/BmrK9fJE7CIAACLiVgDEBwxaXrVu30rx588QOn1gUMG7tRNQLBEAABEAABBKNgG0ChkHKJaQBAwbQ2LFjaeTIkWI5Sf6N\/yudZxMNPNoLAiAAAiAAAiAQOQFbBYwaxE49b4jFDMdWgXiJvOPwJAiAAAiAAAgkMgFbBYwEq25NVmFzRFsOGufkUQKJ3PloOwiAAAiAAAjEKgHbBYy0wrC1ha0uuPQSYN+imTNnikxVQaiKxrlz53qW7sK9r7e2oXMLpz1yZ9v69etFxtLKxz\/Lk8b55\/bt21NpaSllZWWFrgBSgAAIgAAIxAQBIwKmqKiICgsLYWnRPCR8d3bxS\/vw4cM0bdo0mj17Ns2aNUuUKHd\/8c9TpkyxfN+0ZSzc9mzcuJEqKyvFUqQUyuxnNWTIEK8t+5qxIzsQAAEQAAEXELBdwHAbYy3miwv6JaIqyNO077zzTlq0aFGDM52kZcL3rKdA9522mIVqj2\/9VKdxKdRgdYloKOEhEAABEHA9AWMCRi5zqETgA6N3fMgXOEcEXrNmjdiyzhcHEOzdu7eIFBzOfbljTG8trecWqj1q\/dRzrI4fP06jRo3yFKQuLVkvHSlBAARAAATcTMB2AaO+WPBt2L6hoFq5pOUilgWMlfZIAeMb5VmlLP\/GAs5pQWZf7yNnEAABEEg8ArYLGH6BLFmyhMaPHw8fGJvGl3pkAxfBAobvWV0qctsSktX28BKSFYGsnjhuUxcgWxAAARAAAcMEbBcw8oW6ZcsWxH2xoXP55dyvXz+vHV7qS52LjBUnXimkwmmPPwdxFnByvGEXnA2DDlmCAAiAgAsI2C5g1GB2vu2FD0x0I8BffB0++ZuXjnbu3OnxA1m9erVH4KjPWLkfXQ3Dezrc9qhbpWVJcsu4+jf4wITXD0gNAiAAArFAwHYBEwsQUEcQAAEQAAEQAIHYIgABE1v9hdqCAAiAAAiAAAgQke0CBktIGGcgAAIgAAIgAAK6CdguYAJV2J\/zqe7GIT8QAAEQAAEQAIH4JOCYgGHLDI4YiM9BhVaBAAiAAAiAgN0EHBMwaqwS02fu2A0V+YMACIAACIAACNhLwHYBE8wHRt3Ga28zkTsIgAAIgAAIgEA8EbBdwMQTLLQFBEAABEAABEDAHQQgYNzRD6gFCIAACIAACIBAGARsFzC+59DIJaWRI0ficL0wOgpJQQAEQAAEQAAEviFgq4AJdIgeTgjGEAQBEAABEAABEIiGgG0CJtQ26T179ngOGcQupGi6EM+CAAiAAAiAQOIRcEzAhBI4idcVaDEIgAAIgAAIgIBVArYJGK4ALyFlZGT49XVZt24dbd26VZyc3KRJE6v1RToQAAEQAAEQAAEQsPcsJF4mysvLE5hLS0spKyuL\/N1DP4AACIAACIAACIBAOARstcDIirC1ZebMmZ56zZ07FzuQwuklpAUBEAABEAABEPAiYETAgDkIgAAIgAAIgAAI6CQAAaOTJvICARAAARAAARAwQgACxghmFAICIAACIAACIKCTAASMTprICwRAAARAAARAwAgBCBgjmFEICIAACIAACICATgIQMDppIi8QAAEQAAEQAAEjBCBgjGBGISAAAiAAAiAAAjoJQMDopIm8QAAEQAAEQAAEjBCAgDGCGYWAAAiAAAiAAAjoJAABo5Mm8gIBEAABEAABEDBCAALGCGYUAgIgAAIgAAIgoJMABIxOmsgLBEAABEAABEDACAEIGCOYUQgIgAAIgAAIgIBOAhAwOmkiLxAAARAAARAAASMEIGCMYEYhIAACIAACIAACOglAwOikibxAAARAAARAAASMEICAMYIZhYAACIAACIAACOgkAAGjkybyAgEQAAEQAAEQMEIAAsYIZhQCAiAAAiAAAiCgkwAEjE6ayAsEQAAEQAAEQMAIAQgYI5hRCAiAAAiAAAiAgE4CEDA6aSIvEAABEAABEAABIwQgYIxgRiEgAAIgAAIgAAI6CUDA6KSJvEAABEAABEAABIwQgIAxghmFgAAIgAAIgAAI6CQAAaOTJvICARAAARAAARAwQgACxghmFAICIAACIAACIKCTAASMTprICwRAAARAAARAwAgBCBgjmFEICIAACIAACICATgIQMDppIi8QAAEQAAEQAAEjBCBgjGBGISAAAiAAAiAAAjoJQMDopIm8QAAEQAAEQAAEjBCAgDGCGYWAAAiAAAiAAAjoJAABo5Mm8gIBEAABEAABEDBCAALGCGYUAgIgAAIgAAIgoJMABIxOmsgLBEAABEAABEDACAEIGCOYUQgIgAAIgAAIgIBOAhAwOmkiLxAAARAAARAAASMEIGCMYEYhIAACIAACIAACOglAwOikibxAAARAAARAAASMEICAMYIZhYAACIAACIAACOgkAAGjkybyAgEQAAEQAAEQMEIAAsYIZhQCAiAAAiAAAiCgkwAEjE6ayAsEQAAEQAAEQMAIAQgYI5hRCAiAAAiAAAiAgE4CEDA6aSIvEAABEAABEAABIwQgYIxgRiEgAAIgAAIgAAI6CUDA6KSJvEAABEAABEAABIwQgIAxghmFgAAIgAAIgAAI6CQAAaOTJvICARAAARAAARAwQgACxghmFAICIAACIAACIKCTAASMTprICwRAAARAAARAwAgBCBgjmFEICIAACIAACICATgIQMDppIi8QAAEQAAEQAAEjBCBgjGBGISAAAiAAAiAAAjoJQMDopIm8QAAEQAAEQAAEjBCAgDGCGYWAAAiAAAiAAAjoJAABo5Mm8gIBEAABEAABEDBCAALGCGYUAgIgAAIgAAIgoJMABIxOmsgLBEAABEAABEDACAEIGCOYUQgIgAAIgAAIgIBOAhAwOmkiLxAAARAAARAAASMEIGCMYEYhIAACIAACIAACOglAwOikibxAAARAAARAAASMEICAMYIZhYAACIAACIAACOgkAAGjkybyAgEQAAEQAAEQMEIAAsYIZhQCAiAAAiAAAiCgkwAEjE6ayAsEQAAEQAAEQMAIAQgYI5hRCAiAAAiAAAiAgE4CEDA6aSIvEAABEAABEAABIwQgYIxgRiEgAAIgAAIgAAI6CUDA6KSJvEAABEAABEAABIwQgIAxghmFgAAIgAAIgAAI6CQAAaOTJvICARAAARAAARAwQgACxghmFAICIAACIAACIKCTAASMTprICwRAAARAAARAwAgBCBgjmFEICIAACIAACICATgIQMDppIi8QAAEQAAEQAAEjBCBgjGBGISAAAiAAAiAAAjoJQMDopIm8QAAEQAAEQAAEjBCAgDGCGYWAAAiAAAiAAAjoJAABo5Mm8gIBEAABEAABEDBCAALGCGYUAgIgAAIgAAIgoJMABIxOmsgLBEAABEAABEDACAEIGCOYUQgIgAAIgAAIgIBOAhAwOmkiLxAAARAAARAAASMEIGCMYEYhIAACIAACIAACOglAwOikibxAAARAAARAAASMEICAMYIZhYAACIAACIAACOgkAAGjkybyAgEQAAEQAAEQMEIAAsYIZhQCAiAAAiAAAiCgkwAEjE6ayAsEQAAEQAAEQMAIAQgYI5hRCAiAAAiAAAiAgE4CEDA6aSIvEAABEAABEAABIwQgYIxgRiEgAAIgAAIgAAI6Cfx\/YML+vhnilVsAAAAASUVORK5CYII=","height":250,"width":415}}
%---
%[output:3b7a8450]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Index exceeds the number of array elements. Index must not exceed 8.\n\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('change_point_analysis_V3', 'C:\\github_trend\\aerosol_trend_analysis\\change_point_analysis_V3.m', 150)\" style=\"font-weight:bold\">change_point_analysis_V3<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\change_point_analysis_V3.m',150,0)\">line 150<\/a>)\n     plot(data_m_residue.Time, data_m_residue.(param{i}),'.','Color',couleur( i+1));\n     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"}}
%---
%[output:928ac5a7]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error: <a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\change_point_analysis_V3.m',227,102)\">File: change_point_analysis_V3.m Line: 227 Column: 102<\/a>\nUnexpected ')'. Remove trailing comma or add an expression after the comma."}}
%---
%[output:0df8ef64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:21361001]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:91818b42]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:8af7b3ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:52adb4e8]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 167\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Nov-2006  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 166\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Dec-2006  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 165\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Jan-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 164\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Feb-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 163\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Mar-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 162\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Apr-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 161\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-May-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 160\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Jul-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 159\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Aug-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 158\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Nov-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 157\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Dec-2007  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 156\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Jan-2008  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 155\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Feb-2008  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 154\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Mar-2008  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 153\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Apr-2008  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 152\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-May-2008  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 151\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Jul-2008  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 150\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Nov-2008  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 149\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Dec-2008  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 148\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Jan-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 147\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Feb-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 146\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Mar-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 145\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Apr-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 144\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-May-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 143\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Jun-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 142\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Jul-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 141\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Aug-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg1\nTaille de la serie    : 140\nStatistique T_max     : NaN\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 01-Sep-2009  (indice 1)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n","truncated":false}}
%---
%[output:42366ccb]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Out of memory.\n\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('multiple_breakpoints_while_one_snht', 'C:\\github_trend\\aerosol_trend_analysis\\multiple_breakpoints_while_one_snht.m', 57)\" style=\"font-weight:bold\">multiple_breakpoints_while_one_snht<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\multiple_breakpoints_while_one_snht.m',57,0)\">line 57<\/a>)\n                    result.PrctDiff(k,:)=PrctDiff;\n                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('change_point_analysis_SNHT_V1', 'C:\\github_trend\\aerosol_trend_analysis\\change_point_analysis_SNHT_V1.m', 80)\" style=\"font-weight:bold\">change_point_analysis_SNHT_V1<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\change_point_analysis_SNHT_V1.m',80,0)\">line 80<\/a>)\n    result_m = multiple_breakpoints_while_one_snht(data_m_deseason,{param{i}}, 12, 0.05);\n    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n<a href=\"matlab:helpview('matlab','error_nomem')\" style=\"font-weight:bold\">Related documentation<\/a>"}}
%---
%[output:98450c23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the number of abs coef does not allow to understand the instrument type"}}
%---
%[output:09a2f37b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2376b229]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0ca296e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3965de54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9b922ae4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5fa14e70]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6df78528]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6660adae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:41b1540a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4930ca5f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:47c93f6e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0364c28a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9f147c28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:58e6774d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3e79bf00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1ee48409]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7c760164]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:029bc713]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:75f9e04c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:978da8e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:934940d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7f42177f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:774e4d2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:43dd921f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:624e0fab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8053a40f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:41331dce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:64fa7b72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:47778218]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:74372a75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:96ea315b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:05556658]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:35b96170]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5ad1ec71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:52ddc51e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6bb80cb3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:34ab9ec3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0a152510]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0115a49d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:91f96e64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3b22d618]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4793662f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1f0c69b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:82ee6987]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6a609989]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:45e15f5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3e7b5d93]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:091562a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:76de35e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:76f85a79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:522c7e1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:986a34ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5d14b4b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2dfc0f80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:73d22197]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3c760eda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1e18f3f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2cb9f218]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('tabular\/dotAssign', 'C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\datatypes\\tabular\\@tabular\\dotAssign.m', 508)\" style=\"font-weight:bold\"> . <\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\datatypes\\tabular\\@tabular\\dotAssign.m',508,0)\">line 508<\/a>)\nCannot delete 'expA_br0' from the table because it does not exist. Assigning the literal value [] to a variable in a table deletes the variable. To create a new variable in the table with the value [], use T.expA_br0 = zeros(0) or assign first to a temporary workspace variable."}}
%---
%[output:1a6d931d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90bf4018]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74accf82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f6a210b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:02c66054]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:275fae0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:68a84860]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e5e2927]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b329859]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21853660]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ffe934c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a0916d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9d73ce37]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:633ef791]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c5fe4a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47f9d4f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3eb703dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5dcd7821]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27d03b35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c17ca5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:767ec57e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:806deac2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d1ca422]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b31e208]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:60d29154]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0087bb11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ef60291]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83b3a347]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:10f020f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2aaaaf09]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b601709]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4002a337]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0da292be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34c77553]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f5cfa6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2250e661]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71f1d760]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82b16049]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32ca3c35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90a4cce6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:862eeea8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:78254fee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5af536c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93e2212e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8d19b7b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2856bd3b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:564b1669]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f3e2c66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:29e927b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3813a8dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e2d94fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4fb1dd73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:590060cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:413172ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29a3859e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:484a7bec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6ed9ef45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1a6e0101]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:846b7d01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:000e549e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1da35d86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82f20aa6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d3b876c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31a277a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8fe07ccf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40b2d253]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17d1c91c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f86febd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e0ed078]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b911553]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:298ec69e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ee779b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:248e1962]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83bf22a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84060908]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:08155abe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:73710e2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:018dd32d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:44d899d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:67388040]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22a309cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00e838d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85d54d6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:993ae445]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b31a73d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36662d89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01d4841b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15201ad6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95d59e56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5cb5ef92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a70d0e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:86b42ed4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1fd14807]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34595fc6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:938ffa44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44402cdf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01d04c02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7134ed46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9039e061]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ed933b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:565c79b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e2689a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06e6c066]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4171e70b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:494ede43]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b004c68]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:874d8ec0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2dda597b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ada91c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:128525a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:024dd3c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9becd61c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:348171fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:340be5b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c73e2a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b500bbc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93b6cf28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c7e4d1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:937e65fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9f12749d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14fedef7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f8f71b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43df0a2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32d5edea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1568a178]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:393188c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0e6bd3c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c75499d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92f7e111]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32dc71af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87da0cfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46e9006d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f549bfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d19fc14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44484d29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13089376]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1008e3de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0dfad22a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9e8c34ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30f45b9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51d0b936]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89715877]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0dc7802e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10733847]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b710032]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d5d6928]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3753ca58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18baf1f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c3bab3d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65ff9912]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e9923ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48540a2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8db3c968]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38ea5b7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32e97787]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ff26752]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42dc50e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41376bf1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:907db972]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7981fd14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b36dba7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34d214c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5c93efdd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36b55a29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d5bbf8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b16883e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7cbe1379]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d89d53b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9397a265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b83538f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9bd1108c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93686c5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:860dd88f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c489307]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:367decd9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:546b111d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51142b55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8213be0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90dc90d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65f88126]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c0d4dc4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66e55553]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3ce0b691]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12129d59]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95063375]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b869b28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:031c6712]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7cc6f0d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0cc5f808]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:04817980]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3080e591]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95524bbf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54710bbd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d06bcf8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:077691bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29b4bdd7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32035ee5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01853600]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:658072e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2801f0a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:372b2fea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5554bccb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:506bae5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9126c750]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b24bfa1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89a59519]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3494b737]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d1feeca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ff4a3ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6762415e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e35ec14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:96cdf2f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5adb1609]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89bbea40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8e63cdb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:585421c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a948b09]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:200ea141]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c92ee02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:982c46b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f8834a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:704156ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a4c0247]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47f151da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d21ca6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7fc53de5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:772f5af7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b54119c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51559519]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:69e42c04]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3186ce1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47547439]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e6bc093]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d134474]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47231f80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c0a50bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:078c1650]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d65e458]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c761628]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62f1bbf8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59148269]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9dce30cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:78ed7ead]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c07a296]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06dc8a2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:445548e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19eeab58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40debe50]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:544a3c5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75ffc104]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37c3ae91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:67d65a07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:196d887a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9809bce7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79fc8e81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34c31025]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36558d63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c536d58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f1bc288]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f4a0cfc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3164b319]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37e44cc3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:062bea41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a2e83b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:596aaa7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05ede065]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d42046a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ec52aa6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d511a73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:997613bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6688afbb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4135bf82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26c9cd1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2895fcf6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62d2836f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ee813cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:142b0a16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:67f38a16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34ac7265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e51d712]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64fe1588]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c2330cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2a4c696d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b30af8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4fbb3e67]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:04e78057]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b6c3002]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3024b889]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43a1badd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:481b0c93]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f3232c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4cb41de4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:213e149e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99a27963]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59bca9f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1edc1436]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:92e45a27]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:89a9978e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43a98f96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:711bbcda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79e15497]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19adfb2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5637d336]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e671f54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12cdf73d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3bbddfad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76b100fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:390ccf2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4469a233]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30d4e6d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:713edae2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6dd81bf8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:141220f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06416cbf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b25ddc8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:313fa6c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:912241b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:14d949d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64b1e93d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:931bd380]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ad24684]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:233226df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9df1d2c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ab07325]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:216a19f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f6044e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a1f7854]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e9ee41e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7992a999]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b8ce3f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d883bbc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84193af1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:297ca220]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:17703273]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d2fb422]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:831976a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:113f7888]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:233d1d1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59c6040d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b3e15e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7400dfd6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:801873b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c3fb85d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ee00466]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:661c10e9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:278dc67f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:612d2e14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a31ff0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:137e2f57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:96b7b375]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:175c0342]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7bf2a46a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:325839f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03ac8eff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9c6667b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:795d3d1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22b85d3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:918473fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11001239]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:984e4182]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:639ee66e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b67ce71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e78f5c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:56059f9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87e9be89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cb8b6ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d6c85c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0dfe9da0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:742b8293]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7141f281]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4afc3452]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2623bc39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7fc0485c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4fb6c006]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8561cd64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:20d39efb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11431d9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:048e97ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1852498f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9143a6f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3380886d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47ffb32b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8588c821]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:135e38ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:454c7f7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:879c56fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:111bca2d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:751c5705]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5208834f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76b6c351]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8bb63f0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:94b3a58d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ce8b5ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4bf4381a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:79d6160f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:969ef301]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3773d76c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6dd507d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9b0c1f72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2efecd0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7599f9c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b9cde16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1cf8bdfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a595b1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c96b513]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42879679]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18d08194]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ec79177]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27737cda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b8b3559]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d2f958a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93559b24]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:470c9215]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59c56c76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5237cd0e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:508a7a92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4255bcda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7ec81959]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ade710f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5443e0d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c9e2418]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:875fd8b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b235608]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c8a13af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e48f01d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3051d908]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b4dbad1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5644bfec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c66ba1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:153e65a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e50782e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95bd2060]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:357ff63e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00b81aeb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12b8c40f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01d1f83d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33d2d9e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6db989f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:25e66804]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:457f4c1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7d574454]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0db7d784]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:454b51c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c73b87f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39652ade]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74357a07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b2412d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0903c4a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:234929f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:070a795f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99479e13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:103d094a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0404a91a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4fa9468e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0bde47ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85260890]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10db9db0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4cdc6a3c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:81963189]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b4c6db2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64f5ee8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13139e4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16cbb4eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21ba8c85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d0fb04c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f126aca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a2da152]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:311a487c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2170f141]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6bddefc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0da5f4ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0817b4ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ee7b693]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a52dac7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c16d95e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5195e48f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15c2e6a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d9bd20f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:198efaba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80284d00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72e2ba8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1835e793]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4531b12e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c381965]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:224b96be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:635c9b22]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:445df048]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c09b771]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80cf8149]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3daa6ac9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90bb9036]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:758fd6f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:753ef725]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7094ccf5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ac20262]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2169675d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:187dec53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:38b46c11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c1cb559]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7350707f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66461452]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7cc48811]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b14209d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b1610ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:112c2587]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4de3aadc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:008a21ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43ad4fa1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4edaac14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b903865]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d7cd259]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15800c7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e9948b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f14ceca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:883a501b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c12e34d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7cb08fa9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:469bc765]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34071e90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0fadc2a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97f618c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1e270ebf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15b7c972]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d29dfc8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a45d3fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c48b797]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a26972b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98a0eeee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e52ecfd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34858ae1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a5ad6af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ccf8032]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29198aa4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:453264d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6652ff6b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:663158e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:833d3a1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:52afc33b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9965b3a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65f3b8c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b734d75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8aa73ebc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99527e6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c00a7e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82bf6ca1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:83e751ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:959448c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6046acdb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d234253]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ab8be3f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9cfffedf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:91fbd516]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ede8235]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:697cf7a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:547b319e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18430d68]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16b30a3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:941a5882]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3380b4bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:788649f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f54f3f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06b5f3b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b2ef57a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d6f898f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:213aae46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22932d90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f82340c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f29075a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b48a3f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9c4369be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11a0a64f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5157b02b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d6b948c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:25fe6413]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8386dd71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:413d058a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3550837e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75fa8c8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80ebbb5f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1aa4d602]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12ea44e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a332fd1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:091957ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:04d88352]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c6fb989]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9f5e9c1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:72da9bba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1a9ba36d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42f90fcb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e8345eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15513712]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:8e53c04b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c92f5a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30b4837a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:597192b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:73cfa7ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14a98fd7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:19e6b3fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:70003f97]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06c51ea5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6c8bd1a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1a0bcace]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90364ea7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0251823f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:40bdd972]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a040431]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6acadad0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63115f0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0801b24e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6d34fbf6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:472e405b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31966b1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5acb435e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7dea8a02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62a11b41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f4f19a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0608e2c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4a8413df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60e8ca5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9dd9441c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9dbcc9ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5dfab478]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:582e4a1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92d880f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7d6d374f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:948cde65]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:81af6044]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8634c476]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b693640]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3bafc0c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:971d3b0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9bbea782]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0aa815a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:623202dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9e71f5e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7dafc6ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:456551ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38a604db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:960579c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f59913f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cc418b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c4e9014]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:21cc769b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9f24076b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:636f52e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d7ea755]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e89a2d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2483172f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88df6649]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6334f9b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:330ad32a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2dc7fc8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21c0df7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8bb5f7ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3051a516]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:336c33ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a085ee4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c466ec0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:015cb391]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b57dd35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18ac3c72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0da4bf13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30acf795]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8ac12db5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e21333c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:19130dd3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ac06016]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0025699a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c969a4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3205fcd6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:721da8fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d366bc5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:983aa819]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f0b16eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:88db5cc2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a51fa74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75d9574e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f438b82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0f331783]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9dc5979e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:235cad5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58b620ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b41c513]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ddb25f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ab2fc80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f71bef0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37cf23e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:295c7b54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61077f08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e4725a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:50878f63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7abb99b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e23b6c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4278e853]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c936736]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26a8f37b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8afe35af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f7c22ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:50dc6dab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d6c376d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d9c5eab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49c72875]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:84ee9072]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:218cf845]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d1a72fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ad5b377]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:48e441d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:365b2f9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36d3706e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18b338aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5cf97442]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4863cf71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c592f84]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34dd8e79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:586a6115]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:488fb4b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87fce924]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b207c75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3675ed46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43170eda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3251e38d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66d74be0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:68b65577]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7384198f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ce446de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7698e34c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4474b738]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:191152a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:122b0f3b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b03e27b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e4b2265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e2ff7bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3965d6bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ffb08e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c147e00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34fde0c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21ae4975]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10b9104d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6566c6a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f904590]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79bf07a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0186dc17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:363c88f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36b180fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58d87de9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8387a974]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:52c9e2c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1831ac8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66b7fcea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61607279]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:682774aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:28ced20c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:752a74f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b516c0a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4094f7b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:506ff6e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9336c533]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a8a2ede]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:967d09b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d5180de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74344efb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21aaddf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71e11d2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d10f7e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54be5c90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:660ef67e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3398924f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:169784af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:56970fba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06b0e977]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8a77ba81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4722ddaf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e6de565]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a74debd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:842a63d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d183d89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40de2c29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:475024df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3a5a4286]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76076ca3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8912acf4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0aefa498]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2532de9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6fa36ea3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ceeb71c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80b37eb0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9417dabe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a8c7a5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87b7195e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b698fa3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0dfda2b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c5ac0d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c185aa0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:053ef558]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f12dde5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e288b62]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2edbd576]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:714aef69]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49946398]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b78cc28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:803b4b5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a4a7201]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9327ad9a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1dabecb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:941aa8b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93e78f5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1448401d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49aab42f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4664b6ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:387f045f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:24371ede]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cfb6040]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3537b91a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0dcc4ea3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06638abe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88c538af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2462298d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ef77851]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8d5ffc95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d81dbe2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b080ef5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51290129]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58623409]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47d2303a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c9fb410]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e7ed239]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:248e63cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8638c5a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b284375]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:64c0ff0e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:594ca13c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:624d246b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b0f3cfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:89588526]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2505a68f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e72b5f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f950a37]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49dc6a9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:52f25e5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92e68211]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b448e58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f7c6739]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7443477a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8027f5df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a2a7959]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9c60dee2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:91c3db57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f71c328]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20a37044]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d127a48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a03cd28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5aae1920]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77fb2d40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18e5ffad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39dfb265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5caeece6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c8dc615]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:07fd2a8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d4e90cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a545acc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f37d48b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f6022a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:275a974c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f36175b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5906e80b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:887499bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c4aea40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2793d475]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1edb9035]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6c7ae58d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3eb66204]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e363b82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e92393d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:11ccbdfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:587529f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:570e88e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3900e64a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3575c0f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ad760b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38cf78ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0225c240]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a0488ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4698c34d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89f194b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b7682da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:31754962]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2fbe3fed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e2db140]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12f9ea2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16bbe101]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31342d67]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e099250]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82f43d66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:68ac73ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34721f8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5942e1fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:006fe226]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22545626]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48cd084e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:585bfd8c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f6aef5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e2d6192]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7557f411]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95b6eee4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:334cb600]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:20d579e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:254a64c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13e374fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4dd2fc4d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0855a80a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89775020]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4faf8252]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:792b81d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7155133c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58181a86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c5bc103]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e1cfb09]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1924c3bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57f5f2e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95d5428b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8bf89b19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ddc9a86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6175de95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39fe8ddb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27c2871d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1837d5e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f25e7ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c2d92df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:554cc82c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:127911b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29ac48cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:574e224b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:45a4b885]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7de238d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c7678b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63c95266]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e8866ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:259a3e07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99abc3bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87eaa7ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85a3b8df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:292ea369]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92156ff3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2706e187]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9fc54de2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:597c3bd7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ce550ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b2bdb07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5fb95bb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c1f91d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88ecc3d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ce147c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f2f7ac8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93fe0ded]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b9f216d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f667f47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e2fa66c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c58805a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83b4a87b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17d0ba2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f1f6d78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:337611aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66b4883d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8721a920]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d1581b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:60705e91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1cf2a9b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:516360a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:938e7d01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59b8d9ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80e8f246]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1950cda8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a480775]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:976ba4c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:571549bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d5ebce1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:822f0fdf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7bfd20d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82c8b830]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a6f7ea2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34e17c9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ae1df08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8612d969]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:507d9e89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27533850]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:42aff600]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2808f36b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c1fd27c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c83caf3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5be3ca54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:88fda62e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77fb56cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:25ae451a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:142011d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:149100a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37a85aad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84cfd9d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:09cee69c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8792c568]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9933e756]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:419746f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79f87bf1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:31b9d295]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:73151e9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3969725e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6334ea98]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:00fb1fa9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b8fcaa5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b3d7c14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:558d621c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:62c51fd5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:149ad87d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:183cdfc6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:73e4db2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:01a55b95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:035add23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:693605fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8cdc2a0a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ca93415]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9553fce4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:614b130b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3a5413b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16f7fe47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3820925b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85d90bd8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b26cc5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75b252d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a30f061]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84ad7612]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ebcbeb4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b236c5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1903177f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b258047]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12909183]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49782fa3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49dc514c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26974022]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16e80540]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9213149d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16882105]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49702400]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93232c1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:20ffa0a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16f84d9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2cbbba97]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20624876]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47528640]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:44322dd4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:52a7cbae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63b0d367]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49be8e2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2af5650f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9451822d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ae73ee2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4fcee8dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2a580ca8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:393c3231]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:643bc4b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:451b5474]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e391f2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29a99fc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ba3fb7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:91c5b0cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1444427a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47034e22]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77cce9d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:367bb50a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:2bcf67b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3a60695f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c697b5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f3bbe72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:21f632a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:13c77413]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f377d74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4be6877c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bc77e66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9325c816]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4459f20c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:19803d9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3aea5ef7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d17b0a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32b06ea8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:831bf892]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b31a37d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2cdb3459]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:136f79ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5490e1f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:797550e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64978c26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:54c675cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46662cc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:457a8a73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:3384717c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c60e759]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ac55cb8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:35d599c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43bf38ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:770bc4c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f6ca1e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1b2cfd90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d2abde9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93151d9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:297189c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:134a2def]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:306ba9c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65d6e7fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03818a92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:43c3bd7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20ac780d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:64547591]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0828ae10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5a559b7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9956869a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e6ed4f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:645093b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e2234fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06b92d0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9047c90e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1fe031c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d2cbf9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:10e22b02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e8e61eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47b6c75e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b262e29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:54bc881e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64ead7e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3de75cba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:339585c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d38f3af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d6ffd01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80e4e314]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32b644ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1d97f9c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2941dc63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9c937f8c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d577d2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:81d5fd03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98c2f4af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1093bee7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f2679e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:74ec2644]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:69d38380]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b941da0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f788cb6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f4e71b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12f2f3c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9fce8c26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e8acdd1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59e15eb0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62b3010b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5eabad07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d2ce568]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75696798]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6afb3625]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08d722eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c894d46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8429d51f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64b8861a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:762dd105]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4370906a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65bdfde5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1284b270]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:049bce27]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ec7eacd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2bd5c28d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64cf91d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77bd1596]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95c4cf8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e10772d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36cb4271]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4432ddb3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:318b240f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c05cced]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72691f80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29f350e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cccb46a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87053d98]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90f70cbb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:016de8a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9be0c838]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:180fe098]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9dc2d412]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d7488e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64fd04d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6adcddaa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3a65b52d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:415f2770]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:262d8efe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:161b01c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ca79f1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ecb8272]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3806cca5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:837753d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9caa9acd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b4fd6d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d115eb1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3bcb2e2d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:698906f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20d967a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b64e99d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:134338bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66d8092b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ddad206]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92a6b2b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1be50e99]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15341366]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41c7f270]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1cc42e2d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:53ff2ee0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40f8d912]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:133a7222]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:50845773]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34201cd3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:07626141]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:600aaae9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b18bfcb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b98efdd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2908d658]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7fbedf51]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95794977]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:665a6e71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a1b1d2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:216cc764]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:292d2644]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b697c46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31b41ca5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0506da1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:376c811a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:488acef2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5cc8013d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ff8c225]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0724480f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:235421d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74f91fec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c08d450]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:09afa47c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:27b24d5f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3706599c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63019134]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b482c6e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47388617]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e851a44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27e50fef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:448b63c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d63c7e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08aa47d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d77c0d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5db06ee2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59c85068]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c0712a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d629e19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a6f634c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4132c1bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:45f9817a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:640f8f1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ca598ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:759f3d8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82089c30]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a3ee51e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7874c0ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a3aaa28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41fc1ba2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b4f39bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:04f7f0cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1e4e55ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:973717c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4dec8759]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e9efc3b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6904849e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4946c2e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c6051ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0951f871]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15164b74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a076e2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:020d1bbe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8fa77a60]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96472f79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:999ea3ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:025f8052]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b841185]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:23431ee5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f1f4fc8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:698074dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4bf87e74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1623a8bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31664ea2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26dd6d91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b82ecf0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7aa0429e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33414c29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ea249f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e15d10c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3fdedfa2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a1dd1d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b95b246]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14c2dd07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0876fae4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08d1e8d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:415c56cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5708f700]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:86d52a14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3675ed81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:28ce86f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:409b82bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59b5e8bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7570c643]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:69a79924]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ffd4381]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:39cb4046]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e5b1f12]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46078d26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32350b6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1126984d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:611eb2e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4bbdc540]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5207df31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30055fa1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2373c263]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3eb2ea88]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:45cdad1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0af0c331]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:311128d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08004ff4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:742655d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9e22d5b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6044389a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a20a12b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4de3e458]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0975e141]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6780b77f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:343a86fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82787bfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d4941d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d6748ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33d34aaa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:614a4538]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:55a82622]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f0a7ced]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6fad765a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77ff3e93]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3b03efe3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3646c97c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80b8fc47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:78ab6f4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:08e32479]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e66b98f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38190043]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38bb8c11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7811609c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82c8d038]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8071b92e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8861bee4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03fda433]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0861bfb8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:09063253]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38721af6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3878f7e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61ea03b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:815ed6c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41e2b519]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49867fe3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ac276c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f7ed5f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f01b1d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51ead3b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6954b36c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5eacbc6b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:160eb0ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:70e85d0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11acf8d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4af6e2ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e997625]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57109d9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9014a0eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5bbfa668]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57aa0885]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4823c304]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ce2860f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:112a2c68]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e622b17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8138c763]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98bf3f39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43a0ec4d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74bf838c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:905c740a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:158387e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:536df220]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:300f1eed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0c4ef80f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9304604c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:196a9d10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46792483]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37dba0d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5fba7732]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:931b9515]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:734fb0ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:430e41ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ae5c49d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05a948ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2beff939]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1494941a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d1544b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e83da6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36736ebc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ff5d822]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:3001f108]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6eb06652]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:06160112]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95c57c56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c336f53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:84290eb6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:961d2d14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4bf02bcf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:130d1473]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:5a28f104]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:17f8ba08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72771ad5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:04da8cdc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4792d246]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1379829e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e096717]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:06615686]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96b5604a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3a00294e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:2be47f69]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:419be82b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4941d267]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:25fad3b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1fe1a0ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6ecd3117]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31278e7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:976d4914]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4a0aa4e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:463941c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:465729e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3a654664]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8368355e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9f8e31b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97c13263]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20ed1d48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:41b33526]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:098163de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:35f2d269]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:092ad1d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:67de3cc1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:67c33165]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:920781d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:751c01f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:767fba1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6d5cf8e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:88b98ccb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:51acdd25]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51f52304]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4df58246]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:667645a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:333395fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cde3f37]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:02c0e161]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d18ebc1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b78fa34]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5846c457]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d317eb2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3da0030d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f342e70]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5921627a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bb36146]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74a00d3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:010de11b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1ca71d03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ff7198c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3dbdefee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ba9db33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82a4b16c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:031ad614]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f755101]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c0bd465]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8ee9e4a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b2792a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2755528e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10f13a7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:50cf9a75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15e7d571]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:381d2c15]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ced0cd2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3034305f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c9d9f9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:426c7b1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c9dfbc0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:190b393c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74db290e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51974500]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:714530bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e4e60ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15899809]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:333aefcc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:71443a11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:683b1ab2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d02dcc1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5351b481]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3deb909c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d103324]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ae545aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:866198af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:70023555]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6404b111]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a6f2990]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7bd9e406]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0468ee64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ec439ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97afa3a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f886cb8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ab8675c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:21cc4258]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d4b1693]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9700c28a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37ed52b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6cd0d63d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06801188]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a77a511]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76bb44b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:398ac965]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55574715]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5dc3be4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cbed04f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c2d0c14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0858d5bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1957b074]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00160ca7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2231c7a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6349738e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26d81c93]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:25138e8c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d7e0380]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a5c39ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:453199f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51dfbb4a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:74f8b9e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:862b215d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e3e4e8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06bc1bb1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:794f529f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0f1a0b54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b33ae72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2dac6c13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:45f25bc9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37771f7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d42183d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c1f4605]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82375a87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9779c5d3]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"3594×12 table","name":"BRW_result_MK","rows":3594,"type":"table","value":[["'BRW'","2025","10","'daily'","'BsG0_S2S20'","'neph'","'y'","'MK'","-1","-0.0349","0.0198","-0.0901"],["'BRW'","2025","10","'daily'","'BsG0_S2S20'","'neph'","'MetSea'","'MK'","[-2;-1;0;0;-1]","[0.0989;-0.0297;-0.1580;-0.1420;NaN]","[0.1971;0.0263;0.0508;0.0375;NaN]","[-0.0007;-0.0890;-0.3692;-0.3233;NaN]"],["'BRW'","2025","10","'daily'","'BsG0_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'BRW'","2025","20","'daily'","'BsG0_S2S20'","'neph'","'y'","'MK'","95","-0.0566","-0.0352","-0.0781"],["'BRW'","2025","20","'daily'","'BsG0_S2S20'","'neph'","'MetSea'","'MK'","[-1;-1;-1;95;95]","[-0.0471;-0.0054;-0.0454;-0.1967;NaN]","[-0.0028;0.0152;0.0230;-0.1217;NaN]","[-0.0932;-0.0254;-0.1136;-0.2730;NaN]"],["'BRW'","2025","20","'daily'","'BsG0_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'BRW'","2025","30","'daily'","'BsG0_S2S20'","'neph'","'y'","'MK'","95","-0.0326","-0.0206","-0.0446"],["'BRW'","2025","30","'daily'","'BsG0_S2S20'","'neph'","'MetSea'","'MK'","[95;-1;0;95;95]","[-0.0408;0.0038;0.0187;-0.1261;NaN]","[-0.0164;0.0144;0.0532;-0.0876;NaN]","[-0.0651;-0.0067;-0.0150;-0.1645;NaN]"],["'BRW'","2025","30","'daily'","'BsG0_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'BRW'","2025","40","'daily'","'BsG0_S2S20'","'neph'","'y'","'MK'","95","-0.0333","-0.0262","-0.0404"],["'BRW'","2025","40","'daily'","'BsG0_S2S20'","'neph'","'MetSea'","'MK'","[95;-1;0;95;95]","[-0.0836;0.0027;0.0240;-0.0903;NaN]","[-0.0652;0.0104;0.0431;-0.0652;NaN]","[-0.1022;-0.0049;0.0050;-0.1158;NaN]"],["'BRW'","2025","40","'daily'","'BsG0_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'BRW'","2024","10","'daily'","'BsG0_S2S20'","'neph'","'y'","'MK'","0","-0.0273","0.0255","-0.0803"],["'BRW'","2024","10","'daily'","'BsG0_S2S20'","'neph'","'MetSea'","'MK'","[-2;0;0;95;0]","[0.0798;-0.0169;-0.0898;-0.2017;NaN]","[0.1851;0.0387;0.1031;-0.0336;NaN]","[-0.0230;-0.0741;-0.2816;-0.3691;NaN]"]]}}
%---
%[output:4901e2f8]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"23×19 table","name":"BRW_result_LMSlog","rows":23,"type":"table","value":[["'BRW'","2025","10","'month'","'BsG0_S2S20'","'neph'","'log'","'LMS'","0.3227","0","0.0050","0.0360","-0.0260","0.2938","2.1148","-1.5272","0.0513","0.0522","0.0504"],["'BRW'","2025","20","'month'","'BsG0_S2S20'","'neph'","'log'","'LMS'","0.8371","0","-0.0067","0.0094","-0.0228","-0.3724","0.5173","-1.2620","-0.1260","-0.1253","-0.1268"],["'BRW'","2025","30","'month'","'BsG0_S2S20'","'neph'","'log'","'LMS'","1.0222","0","-0.0042","0.0040","-0.0124","-0.2286","0.2186","-0.6757","-0.1187","-0.1181","-0.1193"],["'BRW'","2025","40","'month'","'BsG0_S2S20'","'neph'","'log'","'LMS'","2.2116","95","-0.0059","-5.6283e-04","-0.0112","-0.3200","-0.0306","-0.6095","-0.2097","-0.2092","-0.2101"],["'BRW'","2025","10","'month'","'BbsG0_S11'","'neph'","'log'","'LMS'","0.1625","0","0.0020","0.0272","-0.0231","0.4992","6.6443","-5.6459","0.0207","0.0214","0.0200"],["'BRW'","2025","20","'month'","'BbsG0_S11'","'neph'","'log'","'LMS'","0.4264","0","-0.0026","0.0094","-0.0146","-0.7171","2.6465","-4.0808","-0.0499","-0.0492","-0.0505"],["'BRW'","2025","10","'month'","'BaG0_homo1'","'abs'","'log'","'LMS'","0.6000","0","-0.0143","0.0333","-0.0619","-0.6335","1.4780","-2.7450","-0.1331","-0.1319","-0.1342"],["'BRW'","2025","20","'month'","'BaG0_homo1'","'abs'","'log'","'LMS'","3.2100","95","-0.0289","-0.0109","-0.0469","-1.4419","-0.5435","-2.3403","-0.4387","-0.4381","-0.4392"],["'BRW'","2025","30","'month'","'BaG0_homo1'","'abs'","'log'","'LMS'","2.6550","95","-0.0136","-0.0034","-0.0239","-0.7066","-0.1743","-1.2390","-0.3359","-0.3354","-0.3365"],["'BRW'","2025","10","'month'","'Bacx3_A81'","'abs'","'log'","'LMS'","1.5044","0","-0.0315","0.0104","-0.0734","-2.0191","0.6652","-4.7034","-0.2703","-0.2695","-0.2711"],["'BRW'","2025","10","'month'","'expS_bg0'","'neph'","'log'","'LMS'","3.5852","95","0.0766","0.1193","0.0338","38.9603","60.6946","17.2261","1.1501","1.1526","1.1476"],["'BRW'","2025","20","'month'","'expS_bg0'","'neph'","'log'","'LMS'","0.9290","0","0.0066","0.0209","-0.0077","3.9262","12.3786","-4.5263","0.1420","0.1429","0.1412"],["'BRW'","2025","30","'month'","'expS_bg0'","'neph'","'log'","'LMS'","2.2481","95","-0.0085","-9.3623e-04","-0.0160","-7.8118","-0.8620","-14.7615","-0.2247","-0.2242","-0.2252"],["'BRW'","2025","40","'month'","'expS_bg0'","'neph'","'log'","'LMS'","2.1700","95","-0.0054","-4.2629e-04","-0.0105","-5.5398","-0.4339","-10.6457","-0.1956","-0.1952","-0.1961"]]}}
%---
%[output:99d4fe90]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"23×19 table","name":"BRW_resultLMSlin","rows":23,"type":"table","value":[["'BRW'","2025","10","'month'","'BsG0_S2S20'","'neph'","'lin'","'LMS'","0.3199","0","-0.0300","0.1578","-0.2179","-0.5467","2.8716","-3.9650","-1.3004","-1.2953","-1.3056"],["'BRW'","2025","20","'month'","'BsG0_S2S20'","'neph'","'lin'","'LMS'","2.6771","95","-0.1338","-0.0338","-0.2337","-2.1912","-0.5542","-3.8283","-3.6755","-3.6700","-3.6810"],["'BRW'","2025","30","'month'","'BsG0_S2S20'","'neph'","'lin'","'LMS'","2.6046","95","-0.0653","-0.0151","-0.1154","-1.0342","-0.2401","-1.8283","-2.9577","-2.9536","-2.9618"],["'BRW'","2025","40","'month'","'BsG0_S2S20'","'neph'","'lin'","'LMS'","3.8894","95","-0.0706","-0.0343","-0.1069","-1.1228","-0.5454","-1.7002","-3.8228","-3.8188","-3.8268"],["'BRW'","2025","10","'month'","'BbsG0_S11'","'neph'","'lin'","'LMS'","0.4768","0","-0.0042","0.0134","-0.0218","-0.6316","2.0177","-3.2809","-1.0419","-1.0414","-1.0424"],["'BRW'","2025","20","'month'","'BbsG0_S11'","'neph'","'lin'","'LMS'","1.9336","90","-0.0088","3.0256e-04","-0.0179","-1.2592","0.0432","-2.5615","-1.1763","-1.1758","-1.1768"],["'BRW'","2025","10","'month'","'BaG0_homo1'","'abs'","'lin'","'LMS'","0.7423","0","-0.0044","0.0075","-0.0163","-4.1933","7.1049","-15.4915","-1.0440","-1.0437","-1.0444"],["'BRW'","2025","20","'month'","'BaG0_homo1'","'abs'","'lin'","'LMS'","3.3363","95","-0.0068","-0.0027","-0.0108","-5.1184","-2.0501","-8.1866","-1.1356","-1.1354","-1.1359"],["'BRW'","2025","30","'month'","'BaG0_homo1'","'abs'","'lin'","'LMS'","3.5566","95","-0.0045","-0.0019","-0.0070","-3.1234","-1.3670","-4.8798","-1.1335","-1.1333","-1.1337"],["'BRW'","2025","10","'month'","'Bacx3_A81'","'abs'","'lin'","'LMS'","1.0020","0","-0.0104","0.0103","-0.0310","-4.9325","4.9124","-14.7774","-1.1036","-1.1030","-1.1041"],["'BRW'","2025","10","'month'","'expS_bg0'","'neph'","'lin'","'LMS'","3.4999","95","0.0371","0.0584","0.0159","4.5212","7.1047","1.9376","-0.6285","-0.6280","-0.6291"],["'BRW'","2025","20","'month'","'expS_bg0'","'neph'","'lin'","'LMS'","0.3027","0","-0.0016","0.0089","-0.0120","-0.1875","1.0515","-1.4265","-1.0317","-1.0311","-1.0322"],["'BRW'","2025","30","'month'","'expS_bg0'","'neph'","'lin'","'LMS'","1.8405","90","-0.0047","4.1101e-04","-0.0099","-0.5286","0.0458","-1.1030","-1.1423","-1.1418","-1.1427"],["'BRW'","2025","40","'month'","'expS_bg0'","'neph'","'lin'","'LMS'","1.4507","0","-0.0023","8.8559e-04","-0.0056","-0.2580","0.0977","-0.6137","-1.0935","-1.0932","-1.0939"]]}}
%---
%[output:2b5e2d53]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdQAAAEaCAYAAACoxaaoAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQuwF8WV\/tElCgoqsDfyUi4qlGZ1iTEKQQyYglC6Bf8I7AISeSxGY3QXH6C89IIiKIQssFqKL0DLIFbQ2iXRWJgAFqLXxwbXbHSRFVQEDLK6ohH1Kn+\/hv7Rd27PdM9Mz6Pv73QVlfi7PT2nvz7TX5\/Tp08ftn\/\/\/v3EhRFgBBgBRoARYARSIXAYE2oq\/PhhRoARYAQYAUZAIMCEyorACDACjAAjwAg4QIAJ1QGI3AQjwAgwAowAI8CEyjrACDACjAAjwAg4QIAJ1QGI3AQjwAgwAowAI8CEyjrACDACjAAjwAg4QIAJ1QGI3AQjYELg9ttvpyVLloRWa9OmDX3729+mSy+9lM4991w6\/PDDK3XfeOMNGj9+PO3YsSP0+W7dulH\/\/v1p7NixdOKJJ4p6DzzwAM2ePVv8\/zlz5tDIkSMrz+O03Pz58+nuu+8Wv7Vv356WLl1Kp59+eqXO7t27hTyvvvoqnXHGGXTfffdRTU2Nqav8d0agahFgQq3aoeeO54mAiVBVWf7pn\/6J8K9FixbiZxtClc+DWO+66y7q0aMHvfjiizR69GhqaGigf\/iHf6Cbb76ZjjjiCFH1o48+oquuuoo2bNhQefXcuXNpxIgRlf\/+4x\/\/KIh8z549TZ7PEzt+FyPgCwJMqL6MFMvpNQJxCBXW6oMPPkg9e\/aMTah4YPjw4XTLLbfQ\/\/7v\/9KECRPo9ddfp7POOovuueceatu2rWhTJUsJ7CWXXEI33nhjhch\/85vfCGJHmTFjBv3jP\/6j12PAwjMCWSPAhJo1wtw+I0BEKqH+8pe\/pN69ezfCZdeuXXTrrbcSSAxFtRZVC3Xw4MF02223UatWrSrPwwJ9+eWXadq0abR161Y66aSThPv2r\/\/6r2nKlCm0evVq6tixo\/gNlivKypUraerUqYI8\/+qv\/oo+++wz4XIG6eI51SWMOg8\/\/DCdffbZPJaMACMQgQATKqsHI5ADAiZCVUkO\/3\/BggV00UUXNbFQdYSKSp9++mmFPDt16iTIs3v37nTHHXfQL37xC9EO9kB\/8IMfCBcwLNiHHnqITj31VDrnnHOERXz00UeL\/z3zzDPp448\/pmuvvZaefvppUef+++8XpMyFEWAEwhFgQmXtYARyQCCKUGENIuAIlicsVJUQIZrJQv3888\/pmWeeoenTpxMCiVRLE3ukY8aMET288sor6brrrqP333+fLrvsMtq0aZPYGx06dKhwDX\/yyScVy\/idd94R+6dvvvkmhZF4DrDxKxgBrxBgQvVquFhYXxGw3UP9xje+IUgN1ulhhx3WxEK16f\/kyZPppz\/9qXheJca\/+7u\/E6S9efNmQbKSQGG1ymheGbz0yiuvVAKaYKkigIkLI8AIRCPAhMoawgjkgIAtoQ4YMEDshdbW1lakso3yxV4nAouuueYaat26tXhe57p98sknxXEaeVQG+6o33XQTPfrooxX3rqyDNuAG7tu3bw4o8SsYAb8RYEL1e\/xYek8QsCVUdAdRvnfeeWeFxGwIFRYnrMzOnTtXLFu0pQYXYY8UZ2ERkIRAJZAk9liPOeaYSpAS6mC\/9LHHHhMEKwOcTjjhBE+QZjEZgeIQYEItDnt+cxUhYApKgiWJ\/VMkYNi7dy\/98Ic\/pJ\/\/\/OfC0tTtocI1jH1TWJbYf8V\/4\/\/jHKk8vyrhVY+\/gHR\/\/\/vfi71RuIXhHoZrWD1Gc\/nll9NLL70kIodhMSOoSVq8VTRk3FVGIDYCTKixIeMHGIH4CJgIFS1iTxPHXECAOIMKS7Fdu3aRQUkIOkKwEUgYRIo9UgQZqQV7pggw2rlzZ6PfZdQvflQTPRx55JH05Zdfimhg3j+NP9b8RPUiwIRavWPPPc8RAROhIlJ3zZo1Yv8U5GhLqHDpIn0g0giiIFMSzpKefPLJld7psiIFXbnBVITyYZV0c4SLX8UIeIkAE6qXw8ZC+4ZAnD1U9A3HWG644QZhdZqOzYCAcRwGZ0ZRRo0aJdy\/sDRR1HOnEjcZ8Ys9U1ngCoZLWBbsmy5btkyQNBdGgBEwI8CEasaIazACqRGIQ6hIrLB48WIRYIRiIlTUwR4oSBjnUEHCCxcupAsvvLAi9+OPPy5IVxadK1c9YoN6vH+aeti5gSpDgAm1ygacu1sMAjaEioxEIEEEFqm3utgQKly2ON4ya9Ys0cHTTjtNuH4lKatBRyBc1A2mP1SP2KANNWipGNT4rYyAXwgwofo1XiwtI8AIMAKMQEkRYEIt6cCwWIwAI8AIMAJ+IcCE6td4sbSMACPACDACJUWgWRKqvHkD0Y7BfaKSjgOLxQgwAowAI+A5As2OUOWlykjurbt30vPxYvEZAUaAEWAESopAsyJUWKY4bjB8+HCaNGmSOMfHFmpJNY\/FYgQYAUagmSHQrAhVjo20UnWEun37dsI\/Wbp06UL4x4URYAQYAUaAEUiDQFURKogUycDr6+srmE2cOJHwjwsjwAgwAowAI5AGgaoi1Oeff54uvvhikfdUHnhnCzWN+vCzjAAjwAgwAhKBqiRUDlbiD4ARYAQYAUbANQJMqK4R5fYYAUaAEWAEqhKBZkmoYSMpXb5soValrnOnGQFGgBHIFAEm1Ezh5cYZAUaAEWAEqgUBJtRqGWnuJyPACDACjECmCDChZgovN84IMAKMACNQLQgwoVbLSHM\/GQFGgBFgBDJFgAk1U3i5cUaAEWAEGIFqQYAJtVpGmvvJCDACjAAjkCkCTKiZwsuNMwKMACPACFQLAkyo1TLS3E9GgBFgBBiBTBHwilDVu04HDx5Mt912G7Vq1aoRQLfffjstWbJE\/NapUydaunQpde\/eXfw3J3bIVJe4cUaAEWAEqhoBrwgVZFlbW0tDhgyhKVOm0KhRoxrdd4r7UHW\/yxFmQq1qXefOMwKMACOQKQLeEGrwjtOVK1fStm3bxCXisqDOddddR9OmTatYpSp6TKiZ6hI3zggwAoxAVSPgFaGqZAlC3bhxYyO3ryRMOaKXX355I8KVf8f9p7169RLV+Pq2qtZ\/7jwjwAgwAs4QaFaEqqIi3b99+vShESNGiD8FCRe\/8QXjznSJG2IEGAFGoKoR8IpQJ0yYICzO3r17k87lGxxJ7LmiSLcwXzBe1brOnWcEGAFGIFMEvCFUoGAKSgJhrl+\/XhBocM9VtVD5+rZMdYobZwQYAUagKhHwilDVYzNyfxS\/zZo1i+rq6qhdu3aCdOWxmbA9VCbUqtR17jQjwAgwApki4BWhpkWCo3zTIsjPMwKMACPACIQhwITKusEIMAKMACPACDhAgAnVAYjcBCPACDACjAAjwITKOsAIMAKMACPACDhAgAnVAYjcBCPACDACjAAjwITKOsAIMAKMACPACDhAgAnVAYjcBCPACDACjAAjwITKOsAIMAKMACPACDhAgAnVAYjcBCPACDACjAAj0OwIVU2AP3fu3EpifAx12RI7bN++nVatWkXDhg0Tt96UoZRNJpYnWivKhg+kLZtMLI95ZmGMzBjZ1GhWhKreh4rOz5kzhxYsWCBSEpaRUMtG8IyR+ZMp25iVTR7WIf90iMfMPGa2NZoVoWJyQS7f+++\/n1q1akVTpkyhUaNGidtpWGnsVKJsEzTLEz1uZcOHvzPzd8Zj5idGZqmJmh2hrlixQlw6jgJC1d2Hql4wbgNSVnXeffddmjx5sriTVV54ntW7bNstm0wsT\/TIlQ0fSFs2mVge89dfVox8u8ikqggV+wQgsPr6erOGcQ1GgBFgBBiBwhCAkTF\/\/vzSxJfYANHsCDXK5QtAQKr4x4URYAQYAUagvAggULMswZq2KDUrQjUFJdmCwvUYAUaAEWAEGIG4CDQrQkXn1WMzvvnf4w4e12cEGAFGgBEoDwLNjlDLAy1LwggwAowAI1BNCDChVtNoc18ZAUaAEWAEMkOgqgg1KotSZggfbPjTTz8Vx3hWr15NPXv2FGdlZcKJ4Lshpzz+g\/O0WRQbeVS8OnXqREuXLqXu3btnIQ7ZyPPGG2\/Q+PHjaceOHUYMXQhpI5N8D\/bvJ0yYQDfccEPl3LMLGdQ2bORRMcKzl19+uZApi2IjD96LQMElS5YIEbLehjHJJMfplVdeqUCSpW6b5IEQqkxZyoJ32ciT93em000pp5pHIAsddt1m1RBq0QFLK1eupG3btonJDRNMbW1to7SIcmBRb+rUqTR48GBxnjYrQjXJo+IFEkX9Rx55JHIhkEY5TfIEPzBgiJIVWaBtk0xqfyVpZEkYNvLksRhTddWk06rMmKiXL19O06dPL0yvgzoK+VBGjBiRRn1Dn7UZM1WXUX\/jxo2ZffsmeeR3Js\/v56lPwcUpFj1Zfk9ZDHjVEKopi1IW4Mo2bZUUMr711lvisSw\/Klt5VEwwGQZTObrCLIk8UYsSF3LFkQnj9utf\/5r27t3bKDOXCzni6pA6Ybp8f7AtG3xQ59Zbb6WxY8dm5tlQ5bKRKS+dVq1BEzmpY5bl+NngE1xIZ\/ndh1mmixcvpuHDh9OkSZMy9fhk8X1UFaFGZVHKAtzgZCjdFyq569y+Wa9Sg9aeSZ6gteYaqzjySHdUTU1NZtayOhmaxgwT0KxZs+j666+nefPmZU6oUfKo7jz0IUv3oc2YSUKFLLA05P\/KVKBF6hHejUVZv379MnfRm3RIfl\/wTAUv9HCJke2YqRnmMDeA2LLc7tH1MY8tFJfYyraYUDNy9ehWzjYflvy48rBQbeXJ2u1j86EHld9mEZDmg7GVSU7K2BcP5o5O8\/4wi9B2zPB8lhjZ4CMnxZEjRwqXapbyxFkEoa5cCNXV1YXGMqQdPxuMglZjlotpG3nQZ3UPFVtPH330EU2bNi0XL4PEnAk1rfZl\/LwPLl8JQZYflTrxmFxRWVumQQveRh75DD56aRVmEShl6x5DIJIa4JKVFWYjT\/ATyhIjG3mCE3iW8sTVa8wH69evz3QP3gYjnYs1K722kUenQ1lt9URN+UyoGRNi2uZ9CUrKw0INEmXYfmTWARvqmNoES6j7cVkvOmwxCi4KsoxKLBtGJnmCGGZtocYZs6zdveriOCpwS0dyMn1q2CmANHOhacyKCP7T9YcJNc0o5\/RskVmU1P0tNYJXR1p5kIVJnu985zuVIypyeEzHfdIMo0keuAzzDue3kSlPQrWRR8Uoyz1U1SLEUbAwnVZlzlqeODLlFSxlM2ZFHZsJG7OgyzfL0wZhcwYTaprZlJ9lBBgBRoARYAQ8R6BqgpI8HycWnxFgBBgBRqDkCDChlnyAWDxGgBFgBBgBPxBgQvVjnFhKRoARYAQYgZIjwIRa8gFi8RgBRoARYAT8QIAJ1Y9xYikZAUaAEWAESo4AE2rJB4jFYwQYAUaAEfADASZUP8aJpWQEGAFGgBEoOQJMqCUfIBaPEfASgW3bDohdW+ul+Cw0I5AEASbUJKjxM4wAI6BHYN06olmziPC\/soBUly4l6t+fUWMEmjUCTKjNeni5c4xAjgiASGfODH8hSHXcuBwF4lcxAvki4BWhqjkv1TyUQciCCZ4jIWXXVL4ax2\/zG4Gw7wUW6fnnm\/u2du0hS5W\/PTNeXMMrBLwiVHkrypAhQ0LvnlRJF5cah15mzK4prxSVhS0YAdP3AjJV3bxh4sJCHTvW2i2Mm1pQankvtmAF4NfbIOANoQZvH1CvIZIdhWW6ePFiGj58uLhl\/oYbbtATakzXFH\/UNqrEdZotAqbvBW7eKFdvHGAOuoXXrVtHs2bNIvyvLCDVpUuXUn\/ei42DKNfNEQGvCPW6666r3BwfdcVZ2NU\/27dvp\/rbb6dhd9xhhLhhzRrCp4xrnoIf9b333kt9+\/Zt1AbaRunSpUtk2w0NDfTxxx9T69atqUWLFkY5ylLBV7mBn6+yl0HuFhs2UIuBA3NVw2XjxtH4ZctC34nv78c\/\/nGi78\/UkTJgbpIx7O++yh4mN+ZHn+ZIjEtVESruQ8U+T+99+4w6+3SXLjTwIEnqKs+fP5+GDRtGaHPRokVUX19fqQZSnTdvntY6\/uyzz2j37t1UU1NDRx55pFGOslTwVW7g56vsZZC746hR1FLR7Tz0EVQ63vAiuZ0T9\/szyV8GzE0yhv3dV9nD5D7uuOOobdu2SeEo5DmvCHXChAkVN67O5SsRDLVQN2ygLuedZw30YYaauPRaXhCuq3rXXXfRuEBUI9zS7733HnXo0IFatmxpLUvRFX2VG7j5Knvhcm\/bRi1PO60Q1TN9e6NHjxbeIMRVhBXd92fqTOGYmwSM+LuvsofJzRZqCmWwedQmKAntuCLUbkR08Hi6jXjaOmvXrm2057Nv3z7auXMndezY0StC9VVuDIqvshcuNwKCuuEryL+4+PYgtfr92cRCFI55Cqh9ld1XuXVD5Y2FqhLlK6+8QpdffrmwVkGeCF6oq6ujdu3aiT6GESrcQ72\/9z1rlTWtkm0agoWKQApZfFWeMsptM0EyodpoaUSdw2J8BdDz8SHOWkTpymMyFiLFeGtka\/j+xo4dax3gVEY9t4BLVPFVdl\/lLgWh7t+\/n3bs2EFbt24l+MhPPfXU3Dae4+yhIiDJ4lSdla6jz0yoVlBZVYobAerrB5u33NoFiu1xGETe4oxp2PGaujqi5cutjta4\/PZMCoXFrrotkzfmJvni\/N1X2X2VO3dChW\/8\/vvvF+Q5e\/ZsatWqFW3YsIGuvPJK2rt3r5DnoosuEqtHRL1mXUQAw9ChtGLnzshXwc2LdbaSPC2VaOi\/PEfnq\/K4lNvWstSBDl2ZGXFEIzhBoo0tW7YIN\/vZZ5\/NbnYNqJELFNQ3JWzQpBaEvux+8UWqUTG3SP7g+tuz+XBVt7BLPbd5t8s6vsruq9y5E+oTTzxBV199NZ177rkiEvbwww8X50Ofe+45mjhxonDN3nPPPTRjxgwaM2aMS93StgVCvfjii+mXgwZR77vv1r+vtpZElGEM95RJcLZQDyAUx7LUkS6eP980uSv7ZnHeZxrDIv6ex0LAaoGCzke5chF4BwtUKaGTJI7DRLS1rX9\/6hZxZAavwOJU6oeLcVG3ZXye3H2V3Ve5cyVUhEJPmzZNWAZ33nmnCH9+7bXXhHsFk+LNN99MIJobb7xRHGu47bbbhAWbtkjSRDtz584lROLKUiFUZFDC0RldEu+6OlpXW2s1cdvIikPoWAHL4qvypJXbauIeNy6SdIMH\/cPwh45h0o1rydqMZx518loIxFqgoOMh34suP2+kvkS5hb8eu2XLltH4ENLFuGJ8Ibt6PjztuMhFbx6LmLSyhj2f9hvNSi5Tu77KnSuhysAgpP5D8BDK448\/TkjOoBIdIndBdHANy6Ai0wCE\/R3vlMkfUGfOnDm0YMGCSruNCLV370PNwBoNpDYzfdQgStSJKrrMLr4qTxq5bSduTJQmTJPqhu65YAS2y7aTtmW78EjavvocFrY2pBQMrBPBRYZUgNb6EtJW2KICwYeSUG28FbY4Yatg+fLlXmdmssbcFhQH9Wy2d8ood9KuZxblGyRUZMO45ZZb6LHHHqMHH3yQzjzzTCGzS0IFYaI9kDOs3SlTptCoUaMqCRZCCTUEPdNHbSJdfPiYANTiq\/Kkkdt24k6qxEmfa0IUSRty9JztwiPJQiA4seG\/u8U4EqNuW9h0N42+BNuHrLpcvqbvz5VbWLcvb4NB3nVcYp5W9jheFp+9AkGcMiPUTz75RBAaiA1E+sEHH4ijLkcffXTFBYzfEKCEM5mwJtNmDgJhrlixQriPUfD+Pn36VNy++DsS6yNTERdGgBFgBBiB8iLgy0JGRTAzQsVLYIkiunfo0KH00Ucf0W9\/+1u65pprBLG+8MILIiAJJLdw4UK68MILU48sE2pqCLkBRoARYARKgQByNl966aWlkMVWiEwJFUngsYeJvQmUAQMGiGTzRx11lLAen3rqKfrZz35GV1xxRWrrFO3buHxB7oMGDbLFx0k9074g3Fllv0Hjiy++EAfH27RpE4mJ7opLU\/+dgHywkbgRoMHUkC5liWpL54K12c+UbYb1E3ok\/xanvShZk+inrb7kgbcr\/StKV2wxSos59EXnJpc6BTlcYWnTJ7ZQNShh7wVnTr\/88ks69thjxdEZ7Kc+++yzdMopp1CnTp3osDjZWCJGInFQks3oJqyT5b5YQpESPWban4m6LnPWLLvgl0SCKQ\/FjQANRmCb3u\/iPuyovSVXQTa6ABtT38L+nvTKNJO+JJUnyXNhmIMgoyLBg+9Sz5MnkSPrZ9JgbhMMB11wpaM2WJQtxsFG5kwtVBsBXNdRj80ELxiPG5TkQjbbgJyyK0\/Ux2q6LnPcuHW0bJmbvFMgQQR66e7KjBMBGocoTHdrB\/UkjHhNk5YLfXPZhlygBAPrbN6RZnK3aT9pnWCAU5zFfNzArKQyJn3OBnOdbtou+uN6f5L2Q32u7JgH+5groUprFRYqbhKA+zCOQqcdoLwJNetIyrR4yOdtLK8tWxoOZhuqaZRtyCL5jXjNuHHLaNkyfZ5X6VI0uZN0JJg0AlQXga3D07RYOHgftng0iniRdyvP1X0c3Qg7MiIXKHHaknVtJvck7bp+xnbBG9eb4VpOm\/bCvlGTbublQbLpQ7BO2b0CmRPqV199RUheL1265513nojsRWQtVujIniQLzqhOnTqVzjjjjCRYx36m7IQaR3lsSNAEkI3lZapjm+oVyXTGjl0XaVmajkHYkqDst+nYk4qPfuVuzrqHNpC3Y\/16+tp9GI54be2sr49vRlQ4+GjRVkDYAsWkS8G\/+0KoNtZZHG9GXJxs6pu+ddM3aloUHki0Gn2m3kbOLOpUtYWKICRE9T766KMVbM855xxxjOWOO+4QZ1Dbt28v7jHcvn077dmzR5yFQ7TvySefnMV4NGozb0LFy+NY4DbKY\/p4bEE0fWSwvN56K5okQCBRJBKURbkjQAQ\/6M4WxiFB275icn\/xxRe1uXyj936tcrmLHAd2mSrh9jZniMY5U51LO+6enw0+WVlevhAqMHK9kLPB3aaOzbdu+o6xkDXknzkoip1u2sjtqk5WuulKPl07Tl2+q1evFjl6kYQcHz8+qsWLF4sgpA8\/\/FD8dtVVV4mIXqQbxB4nzqhee+214ndTkckiYAEPHjw4NF0hkvIHkzqg7SIINYlLKXwPLprAVPdjFJa2blrTeMT9+9atxgQ7jZrM2loyTUZx+2euL7JEG6tJT4VuIRBngWZ6UZaWl0+ECpyyWMiZ8I\/6u0k38a1jIWeR2tpSDDvdtGksysuSZnvH5t1F13FGqDJ376ZNm0SmIml9vPzyy4JITzjhBPE7kjjI8v7779Nll11GnTt3FuQI13BUsblgXCXdMgQlxXEpEfXXpkuV16nafDxwP+ImLVl05GzrpnWtnKqF6rrtqPZ0k3tRiwoi802f0lOhkzvOAg33gJpy4iYJOLIZO98IVfYpypth028XdWx1094zYitVtG7abEfIRRreGBU4WFavgC1SYfWcEaokMrhzVXLE3ikO5wLoYAJ8aUm+\/fbbxly+wUvDV65cKdyGMk8wOoj2YBEPHz5c3GqDv2GfVpYiLFS820Z5iOoMe3B2bkW4eEDAUfnHwy77SKtMUc\/L6zKzfEdY23pisnPnupe329chWrikTF9UN5d+IWAOblItz6IsL58JFRd6YOHfsmVL7SCZ9jTj6EyZFrxE4bopdQpzbtxFWp7bO3Gwz6Kuc0I98cQTGxGnJMLg75IA4Zq1JVSZ+L579+4EQt24caPW7Rsk3yChwi3dq1cv8XOHDh3Ev6wLJjaceauvr6+8CouP6dOn0ymnXEoDB7ZwJkLcvU1nLw5pqEuXBrr33oZGlnPW71Tbx+S+a9cuOv7440UqzO3bW1D37u7wjteXcCsA+oDsMDLJR1Bu+R4s0JAMRVfQxujRo7XnK1250G36Gya7zbNF1omSe8OGFnTrrS3EYlUWWInQ7b59G2KJHdYWUn8XseA9ILxeN4M6FTWXJUl+AcxfeuklOuussxrdOIaTIPjnU6lKQlUHCG6xPO5ile+Ea\/zVV18Vkc0yd\/GoUR2pvl6\/GvZBmebP302TJ9eETPANNGzYXpo48cPCugLM4SmpqakRmINQv\/\/9E5zIg8UCCto0lS5dttD27d1DiXDYsGEiBkHVFVVu9UF4W3DHcHCBhufRTtEliHnR8ti+P0zuRYuO+xrvtqHN4BsYNuxjq9eY2rJqxHGlXr320cSJj8fWKQSXgnDTlDDMjzvuOHHtp0+ltIQq3cEIdEI2JUweiCCWblydy1cCb7JQ58+fL\/ZtUfKyUKVs6Nd7770n3guXElw+p53mL5nClfvkk\/soasU9cuS+Qr+JIOYQplWr9JjDOoGL\/dxzG4weBmnJIMoX6TfVtIBwp2Evc+TIkY1w0smtA9LFpOZ6gGxld\/3etO3p5IZFesEFZn3Bd6DGL+hksW0rbT\/iPK+zsvPUqTBdYQt1wgTKyuULBbEJSkI9E6EGg5XiKF\/ausG9JRBqjFu00r7e+nmbgAfUwX5tcBKxuC7TWg4XFdPsoaJvcMOZ7tbG0YQwV50k3sBNfqFHh2Sffd2HhPy+yp5GV2T8ghy\/rPdHbb9R6HDY0Zkw3XTx3dm24auu6Prn3ELFkZa4pWfPnsagJJUo8Q7cWANrFeSJaDKs8OUF5T4RKvrlKJVxXNhD68uPrGvX+CThTAiHDSWN8tUtGKIWC1GBYJhs4xafJxpfZU+74EUkex4BgVI3oY82CznXuhlXl6Pq+6or3hNq2kEsKspXlTvNCjhqRSpJME6iBRyxMVleZf4QbfUh7INNYlXavtOFle7zROOr7GkJ1T6Rgq0mNa0XtCrjfqMudDO59E2f9FVXMiVUlwBn1VZZCdXm3JlckQKbKBK0PWMaPMZi+sgOnM\/bTWef3TiXb1Zj5bLdqA827mTkUi5TWz5PNL7KrpO7KA+SzYI3uFj38Rv1VVeYUJ9\/ni5ODW5aAAAgAElEQVS++GKRoUk9n2qa2Fz+3aW1pCPBOORsCqAwWdYuccmyLdsP1rSoyFJGXdu2cuctl837fJU9jQfJBhfbOnEXvGi3OWFui1PZ6jnbQ82jY6bUg2pkMOSR+6xStjJYqDgLiAvXcVwnmMvWlbWUhSszSu48xj7NO3yV3Ve5MVa+yq6T22aRmkY\/g8+GBfuZ3tGcMDf1tax\/94pQTVG+6lEaSb44ijBixAiBfxkI1VaGtNaSK3Iu02Ik6Udki3nS9rN6zle5y\/KtJRmXMMxNi9SoSNo4cqSJuvVVX3yVWzeu3hCqTerBYAclAftIqHE+QlPdtOTs8wTps+w+TzS+yh4lt2mRGmevNe7+qOkbZz23QSj7Ol4Rqm3qQcAGAlbrqwqnph7MHuLGb3j33Xdp8uTJIiOOTH+YtwxJ3uer3Oirr7L7Knc1YI6rDXGsTC2XXNKFtmwxZw065ZTt9NBD2yuP6tqqpm80TM+RgSltFqYkOKZ5plkSatj1bcj+ATJT07WlAY+fZQQYAUZAIrBvXy\/auXNFJCAtWmynmprJ1LLloZzejKAeARgdaipOH3AqLaEmTT2os0zVgQCp4h8XRoARYARcI\/C733WhGTP0VqrcHx006HnXr22W7bGFmvGwmoKSdFmTMhaJm2cEGAFGoBECpr1Whqv5IlBaC1UHuXpsRpd6EFdfLVmypNGjc+fOrUT5Nt9h5J4xAoxAGRFwERBYxn6xTHoEvCJUHkRGgBFgBBgBRqCsCDChlnVkWC5GgBFgBBgBrxBgQnU8XMF7WrHvK93Qqvs57Hd5Dg5i5e2uTiu72ifcYbt06VLq3l1\/obYL2CHv1KlTRVPqjUVhGMb93YWMYW24kr3smKP\/uqj7IvQ8LuZhspcZ86hscWXHPEr2vDFP+u0zoSZFTvOcHHS5vwsFxm\/3338\/tWrViqZMmUKjRo0ST+p+79GjR+XsLOrMmTOHFixYULmWzqGoTZpKKzsITfYvjzzJb7zxRiN8IP+OHTvo+uuvpxkzZtC0adNEHyWG+P\/yXLLN7\/IqwCwwdyW7qlNlxBwYqnEPMoe2Gomfl57HxTxM9rAjeVnoCdqMK\/eaNWtE2kd5teWECRPExfUDBw7MfW5xJfuQIUNynVvSjCUTahr0lGex+u3atSutX79e\/AqF1ll8Mn+vVHpJrvgdz+uINuvJ0oXs6gebpVUaNlxYvKxYsYKGDh1KCxcutF7EhC1ussZc7UdS2dUFWBkxxyJr8eLFNHz4cJo0aZL4JoBr2EKzTJiHyW46ludoOgltxqQrQQzlyYii5pY4eh4me9FzS5wxZUKNg5ZFXSiwJFTdxAFXaL9+\/ZoQp\/wdpHDbbbeJNmDx9enTJ7co5bSy4yYfWYIXE1hAl6qKOnHoMMSEEud3ma4ylVCWD6eRvcyYSwyDaUMlKRSl58FFbNQ3p5PdB8zRR5X89+zZo9X\/Mum5KktQ9iIxt\/yMRTUm1DhoWdRVSUl+uNhDBWH279+f2rRpI1bq0sWq\/g6iLQuhxpUdfZJFusXyWgyonoCwybqshJpGdnUCKiPmZSVUG8zDZFengDJjHnRNF72IiYN5lFs9b8wtpvxGVZhQ4yJmqB8kVLV6MFm\/\/FtZ3DJpZA+udKPacgm57gKEOPvWcuEQ3OfOw\/2YVnadi0x6R1xiHGzLVm4pn87KK2JrQ7VM1QszomQJyq7DomyYh+Ux9wFzG5d6XnNLkm+ICTUJahHPBN2m0uJEukMEzMybN49U94v6e\/v27XMPHAgSvpwc1BWtjezoE\/aP1WAIuWfmGOJKc8AaVr1KLGEBL3ioLEFJcmJPK\/vmzZtLjbkM7AqSUhFBSXExD5Md30XZ9XzWrFlUV1fXKJjRB8wxRjrZi8A86ZzFhJoUuZDnwly+qC6jHOXHLY\/TqL+roe3q747F1DaXVnY1tD3rPVQVJ9mZwYMHi\/3nV155heSeiw22eWPuUvayY46x0Vl5PmAeJnuZMVdlk9+FPH5XdsyjZM8T8zRzLRNqGvQCz3711Vf09NNPixXsBRdcQOeeey4dFueSRIeycFOMACPACDAC+SLAhOoQ7zfffJM2bdpEP\/rRj+juu+8mWEwnnHCCwzdwU4wAI8AIMAJlRYAJNWJk1H1EHKJHUbOt6FyyDQ0N9PLLL9MLL7xA48ePp9atW5d17FkuRoARYAQYAYcIMKGGgCmJU+7LgVDVzB8ICEHA0c0330z4\/5999hmdeuqpIhAAATp4Hm7fk08+2eFwcVOMACPACDACZUWgaglVjXrDRbZqEgVYpm+99ZYYs40bN4pAFxAqSFL+N85DIWoUz33zm9+k\/fv30+7du+mLL76g008\/nTZs2EB\/+ctf6Ic\/\/GFZx57lYgQYAUaAEXCIQNUSKjCExQm3LHLA6ty3KoFKQg3myVSPhnz88cdi77Rt27aCkK+++upc8vA61AduihFgBBgBRiAhAlVNqMBMJlWXVqiKY1xClc9+\/vnndMQRRyQcEn6MEWAEGAFGwEcEqppQQZiPPPII\/eQnP6FXX31VJCUwEWrQ5YtbTYpITO6jsrHMjEDZEYAHCkVeYlF2eVm+ciFQtYSKPdTVq1fT2LFjxYiAXFHUFHpBC1UXlKSzbMs1xCwNI+AOAVvCsa3nTrJ0La1bt05k6cH\/ygJSxZ2+yMHNhRGwQaBqCdUGnCChSuLFpdZ5XKBtIyPXYQTyQMCWcGzruZY5DYGDSGfOnBkqEkh13LhxrkXm9pohAkyozXBQuUuMgEsEbAnHtp5L2eIQuI508fz5559vFGnt2rVsqRpR4gpMqKwDjAAjEIqALeHAwouy8uQL4hJTlOVpS+BRpBt084YBAQsVlioXRiAKgWZLqGoyZXbP8kfACCRDANabuq+YrJVDT9kSk8nytCV6vG\/ZsmVpxRbP46w5F0aglIQK5cT5z61bt9Jxxx0nsgy1aNHCyWiFXVCLa8hWrVpFw4YNIyRzKKIgNeHevXvFReOu+ptHP3yVG9j4KnvRcsM67Natm3P1MhGTjeW5fPlyp0Rv00nMVWWO\/i1aX2ww1NXxVW5dX3KxUEFwuMAZCjl79myRdQiZhK688kpBLigXXXSRiLJzkfs27JJaeX3RxIkTqVevXuK9HTp0EP\/yKvv27aNdu3bR8ccfL3DwpfgqN\/D1Vfai5cYCNIsjYa+99looMWFeGDhwYCk\/C8xjZS5F60tSbMLkhsHhk9GB\/udCqE888YTIGoTrzBYtWkSHH344TZo0iZ577jkCuYEA77nnHpoxYwaNGTMm6bhUngveNynv5tTdQ4ljMy7eaSs0cv4iRWFNTQ0deeSRto8VXs9XuQGcr7KXQe6TTjrJue7hViYUEDaK6i0aNWoU1dfXO39n2gaxAEfu7jKXMuhLEnzC5IbnElnnfCqZEyrAQvKDnTt30p133ikAwgoVexvYn0FyebiAbrzxRjHxuT7XKd2\/ffr0oa5du4qLp+fPn0+dO3cuxEKFPO+9956wilu2bOmNrvgqNwD2VfYyyI0LHlzuoeJMJxKoIIBJJU6Q6vTp0+mKK64o3TcBN++9995Lffv2LZ1sqkBl0JckAIXJzRaqBk1YnxMmTKDevXtXMhE9\/vjjIrG8vEkejyGICBYkXMO4scVlQdso\/fr1E4Sqy9vr8n1RbcG9gcVFx44dvSJUX+XGWPgqexnktgn+AeFggWyK8kU9\/HNJ0Fl\/t7JvdXV1Wb8qdftl0JcknfBVbl1fM7dQg4SKDehbbrmFHnvsMXrwwQfpzDPPFHK5JFQQ8\/r16wWBy\/fLtIJMqElU3l9SYkJNNt7qU4iSxSUSuqISjqleWckUVjMIU5cpCb\/7ktTBV2LyVe5CCPWTTz4RV5whAAdE+sEHHxD2NI8++uiKCxi\/IUAJVtucOXOc7C2qx2aCe6hsocafZH1Wel9lL5PcYcdYgoQTVc9lZC5IEPEPUUSPOqYjM8H0gsD8xRdfpLPPPtsrDxIvHOPPaVk8kbmFCqFhiSK6d+jQofTRRx\/Rb3\/7W7rmmmsEsb7wwgsiIAlW5cKFC+nCCy\/Mop+iTRmUxIQaH+IyTe5xpfdV9rzltk3fh3qm4yNBYnJ5BEclQRPRm6xmWJ+qOzdvzOPqMm8nuUTMfVu5ECruCV2wYAFhhYoyYMAAuvXWW+moo44S1utTTz1FP\/vZz0RAQpaRr0yoyRWIJ5rk2CV9Mi\/MTUkUksgflN0VoUbtaYYRvYl01f7lhXkSTE3P+Cq7r3IX4vKVL0UkL86cfvnll3TssceKozPYT3322WfplFNOEcnmDzvsMJPOWP9dPSIjg5+YUK3ha1LRZ6X3VfY85LZJopBkD1Ene5zvGykKs9jTNFnXeWCe\/CuMftJX2X2Vu1BCzUqJdO2qiR3wd+zLwkLevHkzR\/kmHAifld5X2bOW2yaCF+oSN\/8untHJbpvGEHufeKcsJhJMqNLax7LG3KWswbZ8ld1XuQsn1M8\/\/1xc5I2EDl988YXIinTOOefQ3\/zN3zjNiAFLFEFJOIKDYCi4lXFgHAVRvoMGDapkSspSwXVtwyr\/8MMPRd99Oofqq9wYA19lz1pu2yAhGQAU51vRyY7Ie1OQEN4BixhH3IooWWOeZZ98lT1M7h\/\/+MdOeSFL7GXbueyhwt37zDPPiHNqb731VpN+weWLZAs9e\/Z00mcQKrKaIEkECghVJnYYMmSIyFTEhRFgBBgBRqC8CCCZxqWXXlpeATWS5UKof\/jDH8Sq84gjjhCh7uedd57YM0XqMeyhYqV8zDHH0F133UU9evRIDSATamoIuQFGgBFgBApFwMeL3TMnVLh5b7rpJnG2K4ww\/\/SnP9Fll11GSHOGBAxpEyKbXL5Id1jUbTPYL8C5W+TyTdvPPLXdV7mBka+yu5AbblZdBiPb+0uljql7mjZ650J2m\/e4ruOr3M1Rz7HV4FvJnFBlpiIclJ46dao2kldmTwKx4kxq2oTIHJTkXg19DhzwVfa0cpsieG21JBgkZPNcWtlt3pFFHV\/lloTKaU2z0Ar7NjMnVFhjsD6\/+93vVnL56sRzmXoQ7avHZmQiBz42Y68YwZo80STHLumTaTC3jeA1yRbMJGSqL\/+eRnbbd2RRz0Zu2wQYWcgX1aaN7HnLZPO+LVu2iPzmPmanCvYvc0KF9Qmy\/K\/\/+q9KqsGgEDifimT5uA3Ghcs3bBCZUG3UW1\/H14+1WlfutkdUQJiSIIIjnyYxvK\/6EiV3Fgkwkn+RTZ\/0DfOy45lkbDInVAj15z\/\/Wdx\/ioQOuKJJvdAbEbdIOfg\/\/\/M\/tHjxYvrmN7+ZpB9WzzChWsGkreTbx6p2wlfZk8odNytRFkkUksqeXEPdPBkmt8l9XoYAGp8w9wHPJBrlnFDlnukrr7wSKk\/79u1FUBCifPfs2SPqdevWjb71rW+J+1HT7qGiPTU5PiKKofB4F982k0RN+LaZZKileyrpBBmXULdu3VrJzesqiUJS2dMhlv5pndy27vMkCTDSS3yohTJirnOP+4JnkrFxTqjYM0VUL8gybgHJuiBUeak4kjngHlZZ2EKNOyLl\/lhte1PGicZG9jRyx0nzh3Pirksa2V3LEqc9ndy27nMcDcTCvahSJsyj3LnBlJJheBWNZ5JxdE6oSYRw\/Ywa5du9e3cmVAcAl+ljjdsdX2VPI7ctCSSJ4LXBP43sNu1nVScod1xrP4vFiW0QVFkwN7lz44xdFnjGeX\/cus2SUNUIXwASvA914sSJldSD2M9V93TjAhi3PpR+165ddPzxx4u0iL4UX+UGvr7KnkZuWAg41x1V4BFCNposzvulkb3IbyIoNzxt6qLcJNtrr71mvNrO1Ib8+4YNG8StXBhLWRAohjHr27dvk2bKgDlkHjhwoG0XjfXU7Qhj5RJUKIRQ5c0ziABGcoM2bdo4vWlGxVW6f2XqQeyhqgWZm8aMGZPbUHz22Wci9SESO2R5VZ3rDvkqN3DwVfa0cq9atYomT56sVQWQ6bBhwwiLyyxKWtmzkMmmTZ3cJ510ks2jos6bb75pXTeq4qJFiwj\/wgpStWL81OIK8+3bW4hmu3RpiOyLrh622err651ggEbYQj0I5VdffUUITJLXsyHd4NFHHy3IBC6BJ554ogI69jmR9OGMM85INBBvvPEGjR8\/nnbs2EGDBw8WOXxV6w8BSihIuA1ChTJ27txZ\/Ja3hQqCf++998R7fUqO76vcGGNfZXchd5iVg0u1R44cmeh7s3nIhew273FdRyc3LH3VSgx7Jyz9J598MrVINt4FvATvUr0LaTGHITxzJlF9fctKH0Cq997bQGrSorB606dvoSuuOC11\/2UDWW1HOBNQ01AmFiouFJ89ezY9+uijlVfiVhkQ3R133EGPPfYYBSN9EeWLLEknn3xy6v7C5YuUazjTKqOO8f9ROMo3Gbxl2Z9JIr2vsruW21UEr80YuJbd5p0u6iSN8k2aAGPbtgNS19Yekt52\/zsYtJMG81mzDpBpWEGs1bhxRKZ6ROO\/7s2y1EORFM\/UL07ZQCaEunr1auFKQuYLDDoGGmdM5TU9+O2qq64SLk+4KZDJ6JZbbqFrr71W\/O6iqMdmgnuoMnOSi\/fEbSON0sd9l8v6vsoNDHyV3Ve5myPmuHYOXjBdSZIAA1YeyEnZHhWkWle3jcaP72b96aou0aT6AhnOP9\/8ShBuFOkqSwIiOrTva265cY0keMZ9R1b1nRMqCHLatGm0adMmcR8pwEF5+eWXBbmecMIJ4veOHTtW+vT++++L9IRww8KKhWs4i8LHZpKjmvRjTf5Gd0\/6KruvcjdHQkWfwo6CwH2Ouc22uLTy1KCdLVsaDqbwq4m1nQQyVYndth\/h9WCh6hcf8hm4c4Fb8AgN+CIununlddeCc0KVLlYEPajkiL1T3G0HwIJ7nDJw6O233xZk265dO3c9VFpiQk0OK0\/uybFL+iRjnhS55M\/ZYm7jPte5c22tQSKYjGYrDxZqmLULN23wwpagTPjvbvYGcQxgDwutG3TnAnPcRsa5fDWQSUI98cQTGxFn2O9oggk1hp4WVNV2oilIvMjX+iq7r3I3VwvVVrejCC7o5k1r5fXvv9Zq7zPcxUwU4sm27W5IPbD0wQ1ipYbOneuzngc7n5mFmjehrly5UiT5lsFH6lnUuXPn0ogRIyo30PAeavxvxWel91V2X+WuZkI1u3PjfHvRVt64cUtp5kzznaHwRi9LHycUR3Bau3adtTvXZz1vloQqA5Bk8FGZ70MF6S9fvpxw\/lXuL8fS1IIq+yo34PJVdl\/lrlbM7d25th9xtJW3bl2d471PW7mi68HNvHbtoTom97jPet7sCBWWKa59wzEZFFiosE5BstiPxXnUKVOmEA4co+DYTJGZkl566SWRSOKBBx5olGfYjSpn14qvcgMRX2X3Ve5qxfyCC1o6Jbgnn\/zt1+7cmY0SJSA2BTd2DRhwKXXvfiABQx4F51FHj8ZVnIfOqOreqzu3apIvTM+R9Af\/fCqZuXyjbpsJA6hnz56Jg5Jk8gZJqCtWrBB7uCgg1LBMSYMGDSL8y6sg7eC8efMEqX7729\/O67Wp3+Or3Oi4r7L7Knc1Yt7Q0Jn+\/u\/PTv2dyQZqa7fRv\/zLpkp7SIEIMpXF9fuiBG\/RYjv17buFJk78kH71q9b08MMDtNXVenGACNPzU089lfDPp1JVhHruueeKVGwuU2P5NNgsKyPACGSDQENDF3rnnWecNA5iqqmZTC1bRqfw27rVTZpDKTTeuXfvMNq379ANXZClbdtF1Lr1qkrf9u3rRR98MNFYLy0Y8CRmlRozrWxhzzsn1KwERbtRKQaDFqrO5YsUh1jpJblaLst+cduMACPgPwLf+94hIkraGxzbRxDRoEHPG5uYOrW3Uxfzc88deicWCCBTU7GtZ2pH93dY5KpVnqSNvJ\/xilCjwFEJNSwoKavzrXkPGr+PEWAEyoeAbYIEBO3U1YVlSjpAqDbFZRBUMJDI5v1cpykCzZJQ0U312EyRx2RY6RgBRqA6ELAhOFigwYQLSK6g5vKNgxaOw4SdI0WbIErTkRmdTHFk4LqHEGg2hMqDyggwAoxA0QiYCA7WJ6xTlyUqaYM8gxpFulnI5LJ\/PrXFhOrTaLGsjAAjUHoETASXZQfCrN0iZcqyv2VrmwnV8YgEMzapt97IjE14ZdjvugxPjkUMbS6t7GqfOnXqREuXLqXu3btnJj7kxT26KOqRqzAM4\/6emeBE5Er2smMODGVqUZwFR2AgShF6HhfzMNnjYJ7GnSv1L47cEmvc+IUik92omCOQaP78K0X2uKyLK9njYJ51n6LaZ0J1iH4wY1NUggldFHKPHj3ouuuuE7f1oMyZM4cWLFiQ2WUBatfTyg5Ckwk05KTpENomTSHiW8UH8uOC+euvv55mzJjRBEM0oMM27PcsA9hcya4mLSkj5sBQ5vDGuXQZy1BE0GBczMNk1y0OyqTna9asqaRgldjjIvmBAwfmPrfExTxM9iFDhuQ6t6QZTybUNOgpz+oyNuksPpluUM07DDLA78j4FHbcx5GY2mZcyK5+sFlapWE4YPGCZB5Dhw6lhQsXarNkhWXPKgJztR9JZVcXYGXEHIss3IM8fPhwmjRpkshiBuIPW2jmsSiQuJswD5NdXQyUEfMghkXPLXH0PEz2oueWOHMvE2octCzqms7DwhXar1+\/JsQpf9dleMrDNYOupZUdaR1lUV1NFrClrqJOHGFZsuL8nhfmEne5oIojIxZgZcZcYigtJZVQi9RzG8yjZPcBc\/RRJf89e\/aIBWcwe1yZ9FyVJSh7kZjHmZyYUOOgZVFXJSX54S5ZsoRAmLhUt02bNmKlLl2s6u8g2iKVPo3s8pYf9Fm6xZDuMY8PVvUESMtDl3YyDlnlITewSiO7KmMZMS8rodpgHia7OgWUGfOgazrsuyiTnktZotzqeWNuMeU3qsKEGhcxQ\/0gKanVpRUVVOKyuGXSyK7rE\/quEq1jqEVzQUzj7lvLNoIXKeThfkwru85FVibMpXw6C7UoN7st5mGyB3U46ptxqe9x5Na5pYt0s6eVvSjMk4wfE2oS1CKeCbpNpVWEdIcImEFifNX9ov7evn373AMHgoQvJ2R1RWsjO\/qEG39AoMEJ1DHEleaANax6lVjCAl7wUFmCkiSJp5V98+bNpcZcBnYF9aGIoKS4mIfJju+i7Ho+a9YsqquraxTM6APmGCOd7EVgnnTOYkJNilzIc2FuU1RXMzapYeDq70VmeEoru9qnrPdQVZzkUAwePFjsESGiVO652GCbN+YuZS875hgb3QLLB8zDZC8z5qps8ruQx\/XKjnmU7HlinoYSmFDToMfPMgKMACPACDACBxFgQmVVYAQYAUaAEWAEHCDAhOoARG6CEWAEGAFGgBFgQmUdYAQYAUaAEWAEHCDAhOoARG6CEWAEGAFGgBFgQmUdYAQYAUaAEWAEHCDAhOoARG6CEWAEGAFGgBFgQmUdYAQYAUaAEWAEHCDAhOoARG6CEWAEMkIAF4qi1NZm9AJulhFwhwATqjssuSVGgBFwhcC6dchDR4T\/lQWkunQpUf\/+rt7C7TACThFgQnUKJzfGCDACqREAkc6cGd4MSHXcuNSv4QYYAdcIMKG6RpTbYwQYATsEdO5cWKTnn29+fu1atlTNKHGNnBFgQs0ZcH4dI1D1CES5c4Nu3jCwYKHCUuXCCJQIgVISqrydAreGyBtEWrVq1Qi2qDpRF9SWCHsWhRGoPgRM7tw4iOzfH6c212UEMkeglIQqL6QdMmQITZkyhUaNGtXozkugElZHJVr16i48g3s9V61aRcOGDaMuXbpkDq7uBQ0NDbR3715q06YNtWjRohAZkrzUV7nRV19l91XuUMxt3bm2Crp1q\/Po32aHuS2WBdbzGfMgbKUj1ODdiStXrqRt27aJi6tlCavzz\/\/8z7R48WIaPnw4TZo0STyjXj4t7wOcOHEi9erVSzTXoUMH8S+vsm\/fPtq1axcdf\/zxFLS685IhyXt8lRt99VV2X+UOw7zlBRc0jtpNoojKM\/s+\/TRlC00fb26YOwcogwbDMIfB4ZPRAWhKSajXXXcdTZs2jbp3704g1I0bN4qLoyUBqbfP6+roLjRGZ3UXO48dO5bGjBmTgZrom\/zss89o9+7dVFNTQ0ceeWRu7037Il\/lRr99ld1XuXWYt9i+nU74\/vfTqmHl+X29etHOFSsq\/432URpSep6aE+bOwM64oTDMjzvuOGrbtm3Gb3fbfFUS6vz586lz586FWKjY333vvfeEVdyyZUu3o5lha77KDUh8ld1XubWYb9tGLU87zY2G1tZSw733UkPfvgcs3pkzqWV9faVtkCr+nuS8arPC3A3ambcShjlbqA6gT+PylW5hk4Ua3Ft1ILZ1E3Bv7Ny5kzp27OgVofoqNwbGV9l9lTsU88MOs\/5OQisiuQMifOvqDiR+cHxetdlhnh7xzFvwGfMgOKWzUCFgmqAkPM+E6v4b8FnpfZXdV7lDCRXnS9XMR2FqikxIkjCDmZLwOwg1JMBpe4sWhH+VsmwZUdeu1h8EAmTef\/99grvRJw+ScHd7IDuCQYMBoT7reWaE+vnnn9Nrr70mIlh1BVbZq6++Stdcc43RL65G6l5++eWVgCQQbb9+\/USgUVgdJlTruSNWRZ+V3lfZfZU7lFBtonx16QWRACKYy1dDziDSyTU1VO\/RVkqsj7AZVEYwKLbcVFL1Wc8zIdStW7fStddeSzg3GlV69uxJ999\/P7Vr164Q1ZBBSezyjQ+\/z0rvq+y+yh3pZofFOH68XgFVd26UioJgu3VrUuP5li3p4o4dxYQtYyTiazo\/kRUC9fX1tGjRIgrOvz7ruXNC3b9\/v1Dg++67j3Bu9Ac\/+AHde++9YgVy0UUX0YYNG8TZz4EDB1JdXZ04f1lUYUJNjrzPSu+r7KWU2\/L2l0jZwzIlSXeuScN3EGAAABPcSURBVE0NhFrkgtkkejX\/PWz+LaWeJxyo1HuoH3\/8sbBODz\/8cPr5z39ORx99NM2dO1ecHf3FL35BrVu3FqSKBA3\/+q\/\/SmeeeWZCUdM\/xoSaHEOfld5X2Usld8zbX6xl17lzbdRUE+AkLdQsCBXz2fLly8W8VltbSzhuh\/\/Norzxxhs0fvx42rFjR6V5mz7hiOHUqVMrz2AeHjFihPhvbJctWbKEOnXqREuXLhVHEvMo6Atwmz59uvBgXnzxxWyhRgEv9zKxrymjbDGwABHuXUSz4pzRjTfeKM6R4n+LOqzLhJr8E7KeIJO\/IrMnfZW9NHIniKbNXHbNHmoWhCqJdKYmmhi\/wevmuoCEnnrqKbrqqqtE0yophSWDCSbAUdOvoo3169eL+dmmLVf9kfOtTB\/LhGqBrI5QASSSM2Al1KNHj8oK6aWXXqJ77rnHGJRk8dpEVZhQE8EmHsp8gkwumvFJX2XPW26QB0ojy8smkAgPBW5\/yVx2jVxZEOqsrxcTOjKVSjdu3Dgxz7ksOkLFbxdeeKEgRNV6heWKORZygtx18Skg265du4pgThAtssn95Cc\/aVJXTXyjkuCvf\/1rWncw2hpeSFi7IEdpNeueO6ASa4UlzBZqDO2A9YmsRp988olw+cLFu3nzZjHot9xyi9hTlS4HAM9BSXwONYZ6Oama+eTuRMqmjeQlNyZLTMhy0pSkCqLon\/D2l1xkDwQ4uSZULDC6aYKfgiMF4ujv8NJznctXnnYAkYFAYaliPoXlOXTo0App4XfV9Qu3L4pKqLfeequIb\/m\/\/\/u\/yvysyz4nn1uxYoXIVId3g0wxh2OOx7tBzMHMdnhOuprZ5Ztgcli9erUAFUnsL7vsMjriiCMICgBlxIeKgbv66quFZSr3VRO8JvUjbKEmhzCXCTK5eJFP+ip7HnKbLDDYXtZXeSu3v+Qhuxh0ZW9XR6ggRWl5x1UvWFbLQNqGAisVe6pJCrwBwb3YoIUKqxIkiHcg2BPxKJhzUTDPgtR0Firmu7feeqsJocJCxQUhv\/\/970UaVBxl+c53viOCSZEPXSVrHFNEGyBISeDSdQy3NPZEdc\/J7T8bQt2yZYtIdnP22Wd7d\/Y3OOapg5LQIAKTZs+eTY8++ij17dtXhEb\/+7\/\/uxhktcB1knXeXNX9oG7KQw4m1CSf\/IFncpsgk4sY+qSvsruUW+fOhUV6vsVl3mu\/5q3+NuOi3P7iUnabV4vve+VKunjq1EZBL6YFg23bWdXT7cPqCFW6aeFmHTRokHClqgSnu0REJsiBdRrcQ73iiivo6aefptGjR1fmbbiM1Rzq6DOeNRFqUgs10jPi0OLPaux07TohVDT81Vdf0YsvvihWgz\/60Y\/Ef2N1d\/fdd4v3\/vSnPyWs5LK8YUV1W+Cdc+bMoQULFlT2CphQk6tWERNkcmkbP1lG2bX7lYEOu5A7atIKunnD8IaFarVLWISFqgit+759tVDDonxVdy72ORG1K63KOFG+7du3p02bNgmXLyxxtAU3LixOFHUPNYpQETil20OV83yYhQrrNmpvGtsN4AvfijNCLUPHMbDSx48BDd6lir9jvwFHe4q6DxU4wYXzjW98o7Bo56Rj5avcZcQ8uF8ZdgQjLebqe5KOO56zslADVkVa2ePKi1SB\/\/3f\/93kWEbcdmR9awve8R5qUnnjPAfjI0ioWSfckcSLW7Zef\/11o7iu96aNL3RQodkRqtxABzYg1D59+lQ2yDGgSD6BfQMujAAj0LwQAKHimJ7NmU3bnsMlHrUoySLK11Y23+pJQsV+KbwvpuIjtqkIFW5dBBzB4kMgEgp+e\/bZZ+mBBx6gDz\/8UCSZhp\/+vPPOy\/z+TwwYE6pJTfnvjEDzRCALQgVSYfuwWZ1DbZ6jcyiGxZZQgQMy8flUEhEqOvnMM88IHzj863BfXn\/99XTJJZfQv\/3bv4msGLj5QC0I04Zi4lhNVsXG5avL1JGVPLp2XeyL5SmvfJcLuW2S7djUkTJZZsHLJaDKlPugtnYWbds20zh06qrcBnMdBrbHPYzCHKyAIyFrTbe\/BBqzkd32\/bb1soyRCGZKQvSry6Mytn30ud4hC\/WXtG9f76+XKubr\/JAnPquMVFlgmYhQ\/\/CHP4gNY+yRnHzyyYSbZt555x3x2+9+9zs65phjxDGab33rW\/SnP\/1JpBxEwBKif7HZnVXhoKSskE0f5WsiHJyNR4S\/zfWWcUgXiGQ9udvmPiA6H+c8jIMkV+VRckdhUFtrd37SKMjBJA\/iPKq6P2qRLjBrzHWyZ0moNlhxnWgEmhLqgUQi+oK\/nd\/8LVRYngj8eeSRR0TGDXxomADgakVawZqaGnHw9\/TTT6\/gJA8q40jNzTffXHEPZ6GAasRZcC+lDB9cERONC5zTyG1POGZJEfgXdTQQxBwMDkwju1kiIttrPolwpjHkphXlRXJVvmVLw8HzeTWNzufZLE7Gjzev\/k19g2WARXKS9HpZY86Eahq98v29KaFGybiN+vcfL7It+VRiW6gfffRRJcfkHXfcIaxRFFioCPU+6aSTmiRvQBYlBAjt2rWLUw\/uq75MSfaE4+bTCWTBozBicvG2kItPIpo2E93atfspmKBIXhOKhi2OjorVvY01jAUxCFOXKQm\/Jz260NwIFeO8fDl97bY\/cDUr8jhklBu\/SXpBjLmrQKtgFiWZ0SiosGq2JvVOatRLmg84DqEe0Pd13rnVYxOqzN174okninRU8rxR2O8YAJmo+e233+bUg1VGqPEJJz3NwUKFpRrXNWx6s36\/Uns1Z0RTuMczytWFAyrhq3JMNFKOKHlra9fRtm0g1fACCzTozpW3qZiwMP29uRCqJFLdVgR+yyA3viCsuMnxTeOBv2OOllmNMCeH5f9VMzMhgQQ8ktgzRi7gYML7OHkFDh2b+SW9\/jr2UKN0kwh5QnwrTKg5j1gRE42LLiaVuwhCRX8x2dnsx9pgYyJmzU1iCS1UXAeGFApWpz6Noi9dukx4jXQljTvX+OIc9q11MmSxpWNyr8vFmw0mtnXiJscH0QUtT2Q4kicetm\/f3iTJDcg1KqE+vI8yI5P8\/zi7H0x4n5RQEZQUtTDEwpEJdcIEClquUCK2UA99SkmJyfZjzKpeGrnjEU5WPWjabtA1rHuzaUKFJQxX4MHLOAyrbuSVhYXatIDgtm1DZhh314FhQtq2TZ\/4Po0712aE0uiLTft5EKrtYtBGj+L0KUlyfJkVDu+RRInMR7jhC+2pN9FI8lVTs\/7mN78RyRYmTpwoEs7IORs5g4Pu5rQu3yOPPBDly4R60G0wQUOc7PK1+2SKmGjsJIuulUbuvPdQbftrsi5sg6lM1jDkObQHqie4cePqaOZMt6nWgkf4XLlzbfBNoy827dsSKiZtGxe5rj0slCxy44sguIS58YVeBPdi4ybHR\/Y31UWs9kV11wb7KK9169mzJz322GMiWBT5AvC\/an5emRNYd4MMW6iNUU3s8kWg0YABA8QZVGmF4gzqscce2+h3\/O2LL74QiZiRAIKvb6u+oCQbYrLdG0w62YY9F3Vu3HYhgAm1Xz+iEO+qmDBRJ7jfphKcrTVk23+ccikyQLIshGryMNjimVU93T5s3OT46vVtkFPeTPMf\/\/Ef9Je\/\/IX+\/Oc\/ixtp9uzZ0+iaN0m2mLORbAFZphBc+t3vfreROxjEi8KEataCxISKu\/HiFqyEmFCrj1ChJ4GrKxupjiScrl2jSQkkYWMxxNFLeUFKMOAoLsGBmMP2WkGkNnm+XbnGpTVc5IUdZSFUXy3UuMnxcbvX1KlTherLO1A3btwoAkfVPdQ1a9ZU6iF6Fwl5HnroIXHc8aijjqL\/\/M\/\/FPvuf\/zjH5skytclvE9qoXbocMjlK5MAwdWslqrYQ8WF4gAd\/xu3ICny3\/7t3xpTEEr3MUhb3noQHLioOtL\/j\/tZsWEvSxZBC3ExKGKiiSujrr4LuW0Ix1THhpijgpGCfYMVpzuiAhIMszh1+Cg3l4k\/W+Q+aNKMrUUcZc2HWcMudCBOGy70Jc77UNf1923jWcF7Xe+hxu13mvoffPCBuFsVGeRAaNhL7dChA5111llpmtU+q0b5vvbaoXm5CF1x3rmDDca2ULMSRG1X+uyRyD54Y4ysF1ZHJVpO7OButFwrvQ3hhNUxka4tMblDBzlH07dmM4Gr51HDFgI21nB6aaNbcK0vNvK6JlS806RLpn14G7mrpQ4TagEjLQkRN77LcHDsNckb4CFSWB3cC4jsTcOHD6dJkyaJZ3QWKiLZcEs9ClZj+JdXgfvl4YcfFvsRPuWoLKvcOtIFMV1wQcu8hpT69m2gNWsa565O+nJY4FdcoZe9S5cGGj26oclxIJvFSVJ5kj5XhL4gonXMmDHOkiDIvoftw2Z1DjUp5mV\/ThIqLk5R5+UwXYHFHHQDl72PpbNQ1Xy8OFSMDXG5F6AmkQjeEq\/WCRKuHAQ1LaH8bezYseIjzKvAjX3ttdeKbFLYU\/al+Cb3qlWtafLkGi28IKbt2xvv1yQdB7Q1b95u6t3bfB2V7Tuef77l13mv21J9\/SFixXsmTvyAhg372LaZQusVoS\/yna6yCqkABjMlIQityD3qQgc34cvl\/Buc+8J0BTeVtW3bNuHbinmsKgl12LBhlQvGYZ0ef\/zxuaGP1dhNN91EV155ZSb7FFl1xEe5n322Bc2e3eLgzRYHkGnRYrsIELrvvi6poUJbAwZspxtucEemqlDAfNq0e2jixP\/nla6gD0Xoi3xnFoSaWlm4gcoed3DuC9MVePB88uJhiAsnVPUAcadOncSNNLNnz664a2Gh2rp8pVs4zELFwE2ePJnq6+tZvasMgYaGLoJMUfD\/33nnGWsEamom0969w5oQc9u2i6h161XW7XDFfBBQt3TyeSO\/xQaBd999V8y\/tgXjiH8+lcIJVQdWmqAktBdGqHLlDGLlUt0IfO970blEVXSee+75yn+qxFzdCJav95iwf\/WrX\/GCuXxDU5EIsSu2JIlUh\/jnUykloaqRuupNB2rWj7A6JkL1aXBY1uwQMEVvyjcXnSAhOwSaZ8tYLFfTghlHu0wp\/DDSpjpIoZlH8ZEk4+BSSkKN0wGuywgkQSDOERUOPkmCMD+TBwLdupnJ0iYLmY9JFPLAN+47mFDjIsb1mw0CNkkisrieq9kAyB0pHAF4Wky5ipks8xsmJtT8sOY3lRABU5KIEorMIjECjEBJEWBCdTwwwahk7PsuWbJEvEW9Linsd\/WsrFrfsZja5tLKrvYJEdu4vBpnibMq6h2Qap7oMAxNvyPgaP78KytJwLOSG+26kr3smKOvulSgReh5XMzDZC8z5uqpCcivxqCUHfMo2fPGPOm3z4SaFDnNc3LQpRJDgfEbLgRAUgqZRhGP6n7v0aNH5dok1JF3HLZr186hlPqm0soOQgtLE5mF8LiRQ8UH8u\/YsUMk+54xYwZNmzZNvFa9J1ImA7H5PUvMXcmu6pSaeSYLvNFmXLmBoS4VqJq8JS89dyV7WJ7wsmCO5PfymKHEfuTIkTRw4MDc55a4mIfJHpWCNivck7bLhJoUucBz8m7B9evXi7\/gTKzO4pMHldWztfKYUNeuXbVEm\/Vk6UJ29YPN0ioNGy4sXlasWEG4ymrhwoXWi5iwxU3WmKv9SCq7ugArI+ZYZOlSgYYtNMuEeZjswUxujqYP62ZMuhLEsOi5JY6eh8le9NxiPThlSOwQR1gf6kKBJaHqJg64Qvv169eEOOXvIAVcuYQCi69Pnz65uCAlsaSRHTdWyKK6mvIYN3Xi0GGIxUqc3+Xdj2WXvcyYSwyD58IlKRSl51LXsbgN04so2X3AHH1UyR93oRY5t8TBXCd7kZjHmQPYQo2DlkVdlVClEmEPFYTZv39\/atOmjbBepYtV\/R1EW6TSp5FdvbxAusXyWgyonoCwybqshJpGdpX0y4h5WQnVBvMw2dUpoMyYB13TRS9i4mAe5VbPG3OLKb9RFSbUuIgZ6gdJSa0uraig9VMWt0wa2XV9ktauY4gbNRfENO6+tVz0BPe583A\/ppVd5yIrE+ZSPp2FqoshKBPmYbIHdTnqm3Gp97a6Arl1buki3expZS8K8yTjx4SaBLWIZ4IuX2lxInsLAmbmzZtHqvtF\/b19+\/a5Bw4ECV91+caRHX3C\/jEs1ajUjy7hVjNnyXbDAl7w97IEJUkSh0dCJZG4sm\/evLnUmMvArqA+FBGUFBfzMNlBTGXX81mzZlFdXR2pgXU+YI4x0sleBOZJ5ykm1KTIhTwX5jZFdfUWDDUMXP1dDW3P+9aMtLKrfcp6D1V3Fd\/gwYPF\/jOug5J7LjbY5o25S9nLjjn0XrfA8gHzMNnLjLkqm5yi5PG7smMeJXuemKehBCbUNOjxs4wAI8AIMAKMwEEEmFBZFRgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcIMKE6AJGbYAQYAUaAEWAEmFBZBxgBRoARYAQYAQcI\/H+pigpFM\/0\/OAAAAABJRU5ErkJggg==","height":282,"width":468}}
%---
%[output:8dcac29b]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdQAAAEaCAYAAACoxaaoAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQe4FsX1\/0cliiV2bGDBhkRJNIoYVMSKLUYFFWPFEsWG\/WcHC8ZILIgtggI2bKjRWFEpwQrGWLEQgYAV\/SfWEMX49zMwr3OX3Xdn67tz75nnuQ967+zume\/M7nfOmVMW+v77779XDWj\/+te\/1O9+9zu1ySabqDPPPFMttNBCC0iBaL\/\/\/e\/VSy+9pG644Qa13HLLNUBSeaQgIAgIAoKAIBCPwEKNItT\/9\/\/+nzr88MPVFltsof7v\/\/4vUtI\/\/OEP6rnnnlM33nijWn755eNHJD0EAUFAEBAEBIEGINAwQv3Pf\/6jzjjjDMW\/l112mfrpT3+6wPC\/\/PJLdeqpp6qf\/OQn6pJLLlFLLrlkAyCSRwoCgoAgIAgIAvEINIxQEe3mm29WF110kTrppJPUoYceqhZffPGaxBDtiBEj1BVXXKHOOeccdfDBB8ePRnoIAoKAICAICAINQqChhIoG2r9\/f3XfffdpDXWttdZSCy+8sPrf\/\/6npk+frr744gu11157qfPPP18ttdRSDYJIHisICAKCgCAgCMQj0FBCRbyvvvpK3XHHHeqmm25SH3zwQU3iVVddVR122GGqd+\/eYuqNn0fpIQgIAoKAINBgBBpOqPb40Vi\/+eYbteiii4pG2uCFIY8XBAQBQUAQSIZApQg1mejSWxDwBwG81f\/0pz9FCsyRx8Ybb6yOOOIIteWWW+qjD9Peeecd1adPH\/X+++9HXt++fXvVvXt3dcghh6g11lhD98Pqg48C7eKLL9bWHtMISRs0aJC6\/vrr9a9WWGEFNXz4cLXRRhvV+syePVvL8+qrr6pOnTqpYcOGqTZt2vgDukgqCJSMQMMJFY10ypQp+rw0rGEG5oXGcUniUEteHfK43BCII1T7Qccff7zip1WrVvrXLoRqrodYr7vuOrX++uurSZMmqQMOOEDNnTtX7bvvvuqCCy7Q1h\/a559\/ro477jg1ceLE2qOJ+d5vv\/1q\/\/\/aa69pIv\/0008XuD43YORGgkAzQqChhDpt2jR18sknq5dffrkupL\/4xS8kDrUZLbqWOJQkhIq2igc86z4podK\/V69e6sILL1Qm1vvNN99Um266aZPkKDZZmvk46KCD1Lnnnlsj8oceekgTOw1Pe3wapAkCgkA0Ag0jVGNywoy0xx57qO22204NHTpUtWvXTnv2snMePXq02nHHHbUncFicqkysIOALAjah3n777Tqhid0+\/PBDNXDgQAWJ0Wxt0dZQf\/3rX+uYbDvEDA30xRdfVGeddZZik7r22mtr8+2KK66oY70ffPBBhZMfv0Nzpd155506Qxla8CKLLKL++9\/\/apMzGcm4zjYJ0+e2225TnTt39gVukVMQaAgCDSNUHJDQTjkr+uMf\/6g9efmIEC5z+eWXa6ckSJUPwpAhQ3SKQmmCgK8IxBGqTXL8N8lO2FgGNdQwQqWPSZQCea622mqaPNdbbz119dVX6\/eJxuaVjSsEjAZ7yy23qA022EBtvvnmWiPmHeRf3jXzfj7xxBO6D5nKIGVpgoAgUEENNSz1ILvmkSNH1l5eds2YoNiN26YomVBBwDcE6hEq2iAOR2ieaKg2IboQKn4IEyZMUGeffbbCkcjWNNmUmqQoxx57rDrllFPUJ598ovNo\/\/3vf9dno3vvvbdOA0oIm9GMZ86cqc9P3333XRVF4r7NgcgrCBSNQMM01DBCJWcvL7xtmuJDNHnyZEmOX\/RKkPsXioDrGSppNiE1tFNTMCKJUxKDOO2009TRRx+tr7eJcbfddtOk\/fbbb2uSNQSK1mq8eY3zEn4NxqEJSxIOTNIEAUGgPgINI1S0T858eKkx+WLi5UVnV4w5ipecJsnxZQk3BwRcCXWHHXbQ7wVZw0xzJVTOOnEswiPeZBYLM90+8sgjOpzGhMpwrnreeeepu+66q2beNX2QATPwVltt1RymQcYgCBSKQMMIlVFx3oNGuv\/++2sTFC79Rx11lML1n3SDn332mTrxxBN1uIw5Vy0UjYw3J2UiZ07jx49Xu+yyi44nDCtLl\/ExcrmHCLgSKkPDAe+aa66pkZgLoaJxomW2bdu2yZqznYs4IyUWlqMV3j1IkjPWpZdeuuakRB\/OS++9915NsMbBafXVV\/cQdRFZECgXgYYSKrtndsq8uLzcgwcPVg888IAmU7sNGDDAi+T4nDdxLrXnnnvqgHnOnuRDVO6CrurT4pySeBc4PyUBAzHZO+20U81yE+bli2mYc1M0S85f+X\/+mzhSE79qsLDDXyDdp556Sp+NYhbGPMymzw6jYVPLMQuew2jMPmxmqzrvIlfLQqChhArUaHUEoOPdCxHx\/1SZMRlceOmDlWjKmiLOdEeNGtUkTMGEGyBDWPiDCWF44YUXtPlakvqXNVvVfk4coSI9xx94tUOAdux1vbAZnI5wNoKEIVLOSHEysps5SrFzZfN34\/XLf9uJHhZbbDH13XffaW9gOT+t9roS6aqFQEMJ9b333lPLLrtsJZPfG+K0PRz5sKFBENLARwqyJfsM\/82ZMOEFFEEnswzXY\/ZdZ511qjXjIk1DEIgjVDx1x4wZo89PIUdXQsWky+aTNII0jkuIJbXXXVhWpKApN5iK0IBkk25DgJOHCgIeIdAwQjXFwzmzgaTYFZfZ8DLm\/JYPGMkk0Ay6du2qTWZopjNmzNDiPPPMMzUNFZI0\/0\/cH9dz3UorraQD4QlZ+Pbbb3U+VDSHr7\/+WpvupAkCSc5QQYswlv\/7v\/\/TWmdcYgcImLXI+T0NnwTMv+adsuNOzUwYj1\/eP9MwBWMSNo3jCqxFkLQ0QUAQiEegYYQaFjYTL26+PewPVZj51iZQYmH5f0zTfOiM\/Py3yXrDJgFtAScqCBmHKjRWaYJAEkIlscJVV12lHYxocYRKH85AIWE2dZDwlVdeqXbdddca8NQchnRNCzPl2iE29JPzU1m3gkAyBBpGqOya+chAPJirlllmmWSS59QbGUxQvZ3OjdsnJVQjkilBl5OIcptmgIALoXJkAAliJbGrurgQKhYSwluMQ1\/Hjh216deQsu10BOHSN5j+0A6xAXLbaakZTIEMQRAoHIGGESrORxAZH5qpU6eqHj16qJ\/\/\/Oe1ahj2yDFd8be8zcIQJsXNjzzySF3RBm3TbmGEGjT5YjImxZs0QUAQEAQEgZaNQMMI1ZhM4yrNMD1FVJvh+cTiUT\/SaKP8a5evChJqmFNSMFF5y15OMnpBwHMEpk+fNwArscYCI3Lp4zkMIn46BBpGqHjFvvLKK9o7Nq4VpaHGPTdIqIZ4qdIRzLcady\/5uyAgCFQYgXHjlCL+nX9Ng1SHD1eqe\/d5v3HpU+EhimjFI1AaoeIVe+utt+pEBzvvvHPxI5MnCAKCgCDgggBEOmBAdE9IFa\/\/uD6HHtr0HqLJuqDfrPqURqhBr94qePk2q5mUwQgCgkByBNA6t902+XVhV4wdO0+bFU02Hzw9vEvphNqpUycdI0ewOW7+eBoGnYE8xFFEFgQEAR8RgExtM2+WMaChYiZOqslmeaZcWykESiNUQknI24vZF1d+EnLj3cu\/JFao1\/g7GYmI73RpdrFl+pObVEjbBTnpIwi0IAQwyTYiaYXRZFsQ1C1lqKURKoB++OGHauDAgeqxxx7TeUJdW1Iv37AEDL1799bVX0aPHq169uwZS+KusiXtx7jJbENFkWAS86T3KrO\/r3KDka+y+yq3N5g3ilDRZDmXpeV4zurrevFV7rDvb6mEagtQ5hkqsa7Ul1xzzTXVb3\/7W9WvXz\/VpUsXLc4qq6yif8pqc+bM0RuLlVdeWQUTSZQlQ5rn+Co3Y\/VVdl\/l9gnz1osvnuZ1yHzN3DFjVKuBAxfwKp47dKiam7L2rK\/rJUpuFA6flA4WRcMIlawslG0j\/KRIr187Zy9J6yFUuxGHSi3JshphQqSHIxNO3okqihyDr3KDia+y+yq3T5ivuv\/+qvXzzxf56iS+9+xBg9SXPXsmvs7X9RIlN4VTXI\/5EoNV0AUNI9SCxtPktuYslWThOD+R9B5CJdWhSclWtoaKTB999JHWilu3bl0GDLk8w1e5Gbyvsvsqt0+Yt5o4UbXaccf67wiORsY0G9XTpU+CN3HOI4\/8GP\/qeJ2v6yVKbtFQHSe+jG62ZmpSAxpCDUuEX4ZMxhRGXcpVV13VK0LFLOOj3IJ5WSu76XO8Wi8jRijVp084UBAlZ55rrhnfB2\/hPD2GzTmr4xR6hbk1Jl\/lDpuWZqmhQqYkCe\/fv3+Tai9CqI5vZkg3nxe9r7L7KreXm5io2NH+\/ecRKi2uT54xrTzv++8TvbC+rhdf5W4xhBpW2eP3v\/99zSlJNNRE76nu7POi91V2X+X2eb2A+exJk1Sbzp3rW5AwAYfl+62n7SZ97aZNm\/cMR09gX9eLr3K3GEKNWreioSZ9o3\/s7\/Oi91V2J7kdP7bpZz7dlU6yp7t1oVflInc9TTbKtBw2KuJV4\/ILNwPTaS6YF7oq3G\/eLE2+ZvhRTkmiobovENPT50Xvq+x15a54ertKYu6w+chd7qAmm1dmJs5XA7mDc5c9+Wci1RW+yt2iNFS7PJwhUNFQU613Mfmmhy3TlZEfGpdk7sFE7ZkkSX5xpT6SCTYfhcud5zlrIONS4bInXwZOV\/gqt3eEShHyzz77TMu9zDLLqIUXXthpgtBMr7rqKtWrVy916qmn6rSDdtiMJHZwgrFJJxa9jwkpGISvsofJ7RTmQXaoMWNSJwhIvjoWvKIqmLf+wx\/q5tbViRQOPLA2gFLkHjFCte7bNxTmue3aqVazZjlNwdwDDlBzhw0rV3YnyZJ1isK8RYbNuCZouO+++9S4H3ZnAwYMcA7WNVom03PjjTc28dh1mTJzfZBQ7WslsYMLkv4mR2B0zSng3TURAYkBSBBgmvlI88Euo1UB89bPPadWDSRyCRv7B7ffruZssYX+U1lyI9tygwc3SSrB3HzRs6f+vWub9u67ta4LzZihk8b8tFOnZpE0pkUmdnBJIfj9998rvGwnTZqUiBj\/9a9\/6co0tCTJ8c0KiyJUSezg+rr+2M\/XoHFG4KvsC8g9fbpq3bGj8+TN+c9\/5oV6DBiwwIcbzaxWONv5ju4dq4A5CRvQ6OOarek1Qm42O7WNTtI5njJlnhdwA+Y4DlfXv1MkBevXZptt1sSzukVoqOzgBg8erJ599lmNF1Vk4qrGmD7bbLONuvzyy9VSSy3linWmflGEKk5JyWH1+ZzDV9kXkDtpMnfOUQnjiGohji3JV0b4FQ3HPClW82M+Gy43cC60kPs0UCrO03JxWCzJF8C\/ppFzffjw4ao7dWU9bKm8fN966y1dy\/T99993HjKZgdAMu3bt6nxN1o5CqFkR\/PH6SnxoUg7HV9lD5U7ysXXBq6BSYmVjPn2+By8fZN2SEur8mM+y5Q6dorw8gc3NC5pjl+UV1Qci5fgveq83XB3aYMe6NONLRaj2g1xMvmGCGXPuhhtuqI4++uhI2YcMGaLGjh2rhg0bltsZqmioyZdKJT40ycXWV\/gqe6jceX9s7VJiKfENu6wszOtqOWDl2qqkobp4AifJHWzN8QIbD1d8cuzHnG3rMDd8933TVDMTKibgV155RRcK79ChgzPsLkRMnTzOTl999dVEZ69RQkjYjPP0LNCxrA9kegmjr\/RV9lC5XT62SUFMmOLO5fZlYB6r5fygqM5PGlhfZMyLaHFV2nzF5RdG5nrm\/MCIx40dWxnzKmRqm3mjJgcNFfOvTy0zoboOFnLk\/PT66693vaTW79e\/\/rW65JJLEtUPNeTJTXCI2m+\/\/WrVZkRDTTwF3mp5lfpIJoQ9kpRy\/tgqk+IuoXz1uhdNqM5azg8peOuexqHp8dGef2ZXtNyJIK6XcQl527d3vh09f3BdCm2QVlnmVTTk9gnkxqHVp5YboaJxsshnzJgROf6FFlpIa7OYe12cmbhRu3bt1Mknn5xoEuxKM9zj4osvVpdddpl6++23dfk2IdTkS7RSH5qE4vsqe1254xK1Jzlrnf\/Rckgk5Ix80Zg7azk4uUSVXjOVZEiAP78VLbczgMGOYbmDE8xxnJuTbV51XQeu\/eyhJCXUadOmqdq5eGrwyrswF0KFqPr27asYfL32i1\/8oma6dTH5poUB7ZQE+cSuLr744uqMM85Q1ESlQahHHXWU6tKlS9rbZ7oOTZ1YMQrn+lQP1Ve5mSxfZXeWO+xjS87YuBqegLPWWmp6\/+Fq5MimlcfgGngmLP+7ywvgLLvLzUL6uJzBmcvGRnnCcrZ4yCFN7l603CmHG36Z4xyjmUYUp6vdl7PKQw7p77QOWFZZ1kuSuWtxGioLEPIaOXKk+t3vfqezE3GeGtbQUE3GI5MF6Sc\/+UnuYTQQ6qhRo7SZmAah4l285pprqj322EMTmjRBQBAQBASB6iIAyaM5+9Qya6iYbyFS7OIXXnhhIRk6IF\/ImB+XJoTqgpL0EQQEAUGguggMHTpUHXHEEdUVMESyzISa1XT73nvvqSeffFJ9+umnocDNmjVLffTRRzo37\/LLL+8EbpzJlxjatm3bOt0r706c0ZCukbRaZALxpfkqN\/j6KnsucuPTwHmrFTyv7bj9++vKYC5WYfxfApbR2GWbi+w\/mCqNS8aaazZ9JBYxF0\/ReabMpmbdesLnJXcsQHl2qDPHfZKUi1PledR27z7jhyUZHofKmSlOUmeffbZX30imNDOh4lxEekC0yKQa6ksvvaS12ygyNWtuq6220tmZOHd0aeKU5IJSsj6VddZwGIavsucut3XWmjTvQVJny6yyxxWIcfHyTZN1J6vcDsuxsC7IHiyO7uq8Nc8Xujzz6rzj6\/BMSf379y\/N6zjvychMqAiERkhVl4EDB6pu3bo5mWY5e4WA77\/\/fnXaaaepLbfcUl100UXaqxcHpxdffFFdc801ar311kscMmNkwgGJJuXbsi8b3z80H3zwgSJbl0+OYEVinpRQTWSNq2dnFtldq9ONGDFCRWlgRsvh45ykZZE7yXOK6Bsmu8vGQymyS6Gdlpvuz2zSkJs87507d\/bq\/Qybw8yESmKHv\/\/979r8glftRhttpD1oF1100QWeh0PSvvvuq52QPv\/8c3XcccepVVZZpabZEtry8ssvq6uvvlo7NpHiEA2Wl2K77bbLvAYlsUN6CJvbhyY9EuVdWTTmji4JesD4hkB0QcuxFcLZBJi0srvmrTDZ9KIyJaXVctLKXd6qiH5SlOxxG4\/p00l\/kWzjkcd4zSbNZ8yDOGQmVLuQdxzIcWEzlHhDKyXQePXVV9fhDmixkDYZk8JIut4zqRxhQmbseqgShxo3Uwv+3edF76vsRcudVxbDsBz7aWV3lSksYyIxjlljFtPKnfyNyv+KerLX23iMHHlok41S\/pKF39HWUH20IBWioZrwF5d4ITtsxngHU7KHeqXGTHv88cdrQkXTpRGSg2aZtB6qTfRi8s3+ijTXD012ZIq7Q9GYu2qDLiMM5l93kT1oPk5qhk56rusyDhe5Xe7TiD6usgc3Hi7rwOTAqFfYJsmYrWyPXmdhy11DTQKi3dc4M\/3zn\/\/UWikOR++8844+EznzzDPVbrvtptLWUUUzxSuYmFjOdoMFxkVDTT5rri9r8jsXf4Wvspchd1wWQxcvYGYwqDFOnTpXoXV07txmgXOxekmekjilFpAx0euPe5b1ErcOmF+Oo+P6uaQYDmR79BrzwgkVEvziiy\/Ud999p89B0WBJ3hDWxo8fr4488ki1ySab6PPUjh076jSDkC35d0nAALmirSbN5cvzosq39evXr5YpiTNcfspqLHqK6a688sqJchOXJV\/Uc3yVm\/H4KntZcs+vQa6ef751bfrbtZurQ2UGDnQP7frPf+aoiRNb6WuCZ61Dh85VW21FEpjWdct3JlnnPC\/vVhbmecudxzqPWgdnnz1Xb5hMi+sH6fbt++NassfKujrggLlN1kAU5i2iwHjUQkArJDsRpll2ppyX4lxEQnyckTDlBsNeOCMlnuzSSy9Vpvj4hAkT1IknnqjPT2mAeuWVV6pdd9018RqMIlT7RsSoHXzwwYnvnfYC4m7vvfdetffeezcsFjaN7L7KzVh9lb0Rcs+a1Urx0aPx3926re68XHr2\/FKNHr1UZP+4vzs\/SCnVpcscNWrUB0kucerbCMydBHPolKfs9jqo9+iofs891\/qHUMflFtik9ev3L8U6sFuU3MTqu4ZKOsBTSpfMTklISaICvOpwKiJhAjsOwl8g03PPPVc988wzavPNN9dm2JVWWmmBgaGJksCB2qgQ6NNPP61uuukm3e+www7TITULL7xwJCDGVEzBc7syTRShUujcJHYoW0PlPJgx3XzzzYrzY1+ar3KDr6+yV0HuxRcP1zQauW4xGRqNN285qoB52jFVVfY4co6Su8VqqJADMaR45O6zzz4KwjKOREsuuaS688479d85zyQMJq6hnWI2Nnl\/4\/pH\/T2KUBt5hlqF0J00ePoqtyHURlcZ8hVzV6\/bNONLc01IgZg0t4m8RtZ5rnA63cxnzIMDzKyhfvXVVzo0haouJlNS0DOXsBc01X\/\/+99aayUO1TSckjgfXWyxxXSZNe5j7jl16lQdLkPAb5rmcoaa5r5ZrsG8QSIL+xw3y\/3KutZXucHHV9mrIPf48RSX2KKsZVZ7ztFHP6cdYObM+fHZrVrN+iG366zEaRCTCF8FzJPIa\/f1VfYoubFy8uNTy0yoYbl8w0Jdwn5nyr7NnDlTe\/Wef\/752pEJEzKFyG+55Ra14oorquuuu06tv\/76mXHFrAyZPf\/885nvJTcQBFoKAl9+2VPNnj0odLiQXOvWzyv65Nnat1+7dru5c9spniOtZSGA0sGPTy0zoZp4Ujx18cgl1jRInib8hdy9N9xwgz5oNmXfHnjgAZ2nl+xKwWoyhnDJkkToSx7J5CFVfqQJAoKAOwJTp7ZTN97YLizHvtpgg+fUr36VnxZL6MXvf\/+cu3DSs1ki0CI1VJPN6Nlnn9VkSaaSIKEaYiTJPaZfiDGMiIOrIoyIm+XKkUEJAh4hEFbPPK+z1mCMokewiKiCQPZqM2CI5km5nbXXXludcsopuhwb+X05LyUfLwSLaZickmiyNNeyb2kzJcncCgKCQHkIuGbbqRf4X7TDUXloyJNaKgKZTb4Ahyb5l7\/8RZ1zzjnaOzfYSOxAooa99tqrZtY1eXb5l6T4P\/3pTxe4jrNUPIO5HsclPIalCQKCQDURiMuiY7Lt1MuUZCcQqOYoRSpBIBqBXAjV3J54Uoj18ccf17GolMraeuut1R577BHqrWXCbU466SSt4eLhaxpEi0Z7xRVXaKIuM\/mCLBhBQBBIh0BSsgwzH6d7slwlCDQegVwJNelw7IQQaKicv5LAgXSFJHBG20WrxfvXDrVJ+pwy+xNzi+wm4T8m6z\/96U9aBLT0\/fbbT\/931O9NTFawfxljyCq7PabVVltNFzmgnm1RDXlxhKPZlYyiMEz6+6Lk5r55yV5lzA1ZBqs+Mf5GrPOkmCNnmOxVxtzI++CDD+rle9RRRzUpPmJqRNvfoqqs83qyl415WkwaSqgITczpHXfcoTMjkbLQNIpBk1God+\/e3ph6zaSbRcxHg9+RjhHt25SSM4Qa\/D2hQZxBn3XWWRoG4nIxhy+\/\/PJp59f5uqyyQ2h2qTznB6fsSHYsGx\/kJ1PW6aefri0aQQx5TBi2Ub8vEvO8ZLfXFOUJi25J5QbDsKpP\/K7sdZ6X7GEEWyTuSeUeM2ZMbUNvsOcbuuOOO1Ye8yjZsXCW+W3JMp+5ECqevoS\/UDVmGiUgIpqtRYR1QWMlMT51T33RSM042P2uueaaioT\/NDTUMI3P1GsMarH8nuvDCLjoj2UestsvbJFaadTaYvNCLmlyJJP72XUTE7W5KRpzexxpZbc3YFXEnPc9qupTI9Z5EsyjZLc3A1XEPLhuwbmR35YkmEfJ3uhvSxKCzYVQH3744VpC+\/bt2+vkDGGNuCIyH\/mW8DgJoCxgQ6hhGiqmUAoBBD8o5veQAg5YNHZlXbt2rZmJk8iRpm9W2Y05iWfbpqY0siS9xv5whGHIZiXJ741pPqkcafpnkb3KmBsMwzKWNXKdm02UIZp6ssQV2KjaOrfXrU3+n376aej6r9I6ryd7I9d5knc6M6EaE8hrr72mq8tQgi1Jw+RL+kF2e8SnoqVy9vbUU09pR6a+ffuqn\/3sZ0lu2dC+NimZF5czVAize\/fu2psZ7dWYWO3fQ7SN\/NBkkd2cGTNmsybK2gzYlgCj7QU3JVUl1Cyy2x+gKmJeVUJ1wTxKdvvjUmXMg6bpqPeiLEJNgnk9s3rZmCclk8yEanZw5Ns1mZJchSDlIDu8ddZZR2tlSyyxhE45SHJ909B480o96CpXln5BUrLvZTSR4CKuilkmi+xhYzKaehY8464NYpr03NpseoIm4jJMvlllDzORVQlzI1+Yltcok68r5lGyB9djvXcmbu0m+XsSucPM0lHvRZXWObK4mNTLwjzJ\/Ji+mQn1888\/18XBKb1maylxwhC7ikY7bNgw7UjSq1cv7dVLNZpFFllEn4OZe\/fo0UPn4A2mJox7RiP+HjSbGo2TdIeMk9qvtvnF\/v0KK6xQuuNAkPDNB9ne0brIzpg4P2YNBD+gRc0DWKPV2x+FKIcXZKiKU5Ih8ayyk4Gsypgbx67gemiEU1JSzKNk572oMuaMk6gIymnajnU+YB4leyMwT\/vNykyohhhJPYhTkuv5KKbdk08+WVFE1lSpATjiTY844ghNoDTcu998801NvlFns2kHX8R1UWZTnmWXjbPdwO3f2+EEZZeZyyq7Paaiz5ZsnMw8mlq4L7\/8sjJnLi7Ylo15nrJXHXPmJmyD5QPmUbJXGXNbNvNemBCZqmNeT\/YyMc\/CC5kJlYd\/\/PHHOqPR6quvro455hhFyEu9guCPjgtvAAAgAElEQVT2QkW7MJotoTPUTUVrJSG+2VWa2qpFhjJkAVGuFQQEAUFAEBAEMhMqSe7PO+889cYbb9QNmQFqO2zGJMffbLPNNKGaw+YXX3yxlhDAJN4nFMcXDVWWlCAgCAgCgkDLRCA3QnUpiWaHzRBvChH\/4x\/\/UEOGDFGffPKJOvzww\/V5GA5KBK2TdP\/4449Xu+yyS27l21rmNMuoBQFBQBAQBIpGIDOhZhEQUy5nbSahPmEzOCNR5o1D9YceekittNJK2vN3o402yvIouVYQEAQEAUFAECgUgYYSKg5NkydPVtdee61Opo9DEh69\/\/3vf3VSA5LtE4rTqVOnQkGQmwsCgoAgIAgIAlkRyJVQ0TSffvpphZcleVXJREKihy233DK0PFs94U0KwqwDlOsFAUFAEBAEBIEyEMiNUCdOnKgow0Y8YrARX4nr9vbbb1+LJcXh6Ntvv21Ssi1swGirkDSJI8JqppYBkjxDEBAEBAFBQBCIQyAXQiXAnBSBJGLgX0JeiBnl\/6mPOnLkSLXYYovphOUdOnTQMhGbRnjM0UcfHRlfSkpCwmgw\/XKthM3ETaf8XRAQBAQBQaBRCGQmVDRNgm4feeSRJoRpD+itt97SHry2t64J9l577bUXqHeKVkr1GpNV6De\/+Y1O\/uBbBZpGTao8VxAQBAQBbxGgkC5trbW8G0JmQjXpAcm5e+655+oE98EWFk+KQ9J9992nnY523333GqkarfSJJ57Q5cwIrSFFW1yiCO+QF4EFAUFAEBAEfkRg3DjyJirFv6ZBqsOHK9W9uxdIZSZUo2naGY\/CRo4WG8x4ZJPqTjvtpEjyQAKHzz77TO2\/\/\/46R3CbNm28AFKEFAQEAUFAEEiJAEQ6YED0xZDqoYemvHl5l2Um1LCcvEHxMeGivf773\/9Wl19+eRPTLaTKOes555yj41HXXXddbULeeOONvUiGX95UyZMEAUFAEGiGCKCRbrtt\/MDGjq28ppqZUO2qMVdddZXq1q1bEyLk7xMmTFAnnHCCTnqP1hmsGkOfMWPG6IT4O+ywwwJnqvFISw9BQBAQBAQBLxGATG0zb9Qg0FDRVCvcMhMqY3vvvfd02bV33nlHJ2bAfEvVGfL1Pv744+qxxx5Ta6yxhv5b2Bkr9\/jf\/\/6nCL0hhnXnnXfWBcdNW2aZZdS+++4rTkkVXkgimiAgCAgCiRHAAal9e\/fLvv\/evW8DeuZCqMj97rvvagejSZMmLTAMYkjJfHTBBRdowkza7KT6Sa+1+5NvePTo0apnz56KvMJFNBywMF0TMxu1eSjiuWXdU8ZXFtLFPUfmsDhsy7izl\/MX5bmblFCnTau0929uhGq0TLTVGTNm1NYVnrpt27bV\/4+zEebduIZme+KJJ6rXX39da7V4+q688spNPH2NMxQEbepgcl+I+8EHH9SPCNbkNPUAi6wzSgrFDz74QJewa926ddxQvfu7jM+7KVtAYJlDv+fQq\/lz8dxdaCH3CXHgD\/eb5d8zV0LNSzyyKq2yyirazHvWWWdpj1+8iO2G4xKpDffYYw9NovSByKdPn64r0xjC7d27t9pvv\/30pUKo2WfIq5c5xXCb+\/iApLmPUcaXYuEXcYmr567rGSqhMzgmVbjlRqhonjNnztQeu3\/961\/1S4uGxnkqcaZR4S\/EnVKujUxKF198sa6L2qdPH7XEEktor2BiUHFiMkXIwdKQJb+DaO+8884akYaRbpBQ+\/Xrp7p06aK7Qtz85NUY94cffqg1akrQNbcm4\/N\/RmUO\/Z5DH+av1cSJqtWOO8YCPXfMGIUJu\/Uuu9Tv60k8ai6ECiDDhw9XgwYN0uAEG+eJpBCEWG0PX5OyECLebbfdtHcv12Pu5YyT8m2LLrqorjZzzTXX1AgKQj3llFO09orzEoT6zDPP1OqoGtK1+9gaqi3fIYccoqvc5NUIESJVIhsINgnNrcn4\/J\/R5j6HC82Yod\/Bn3bqFPsOtpo1S0\/o3IJ8KopYLT7M36r7769aP\/987PC\/7NlTzR40SC01erRqc9pp4f0hUzx8+\/ePvV+jO+RCqA8\/\/LAmwU022UQXBN900001+RGj+ve\/\/11ddtll2mnpT3\/6U810a1IWkmJw8ODBWmOEbG2yROs94IADFF6+ZFUyyfHjCBUt15iBbVOxMflC\/OZcN28NlWd\/9NFHWuttjmeoMr5Gv7LZn99s55DzugEDmnzIIcq5Q4cuGL+YpG92yHO9Q+Xnb\/p01bpjR+cxz\/nPf3RfrdUOHLhgpiSI1IOkDowhM6Ea8qJc25AhQ0LNpyasBk0TT1+0ThyPCLWBhPEONpqrbc6FZNEgCceBtAnFMdonuYHDTL5BsrVnVc5Qndd4ZEc5n8qOYaPv0Czn0PW8DvCT9G30ZIU8v\/Lzl9Fzl\/HNnjRJtenc2TulJDOhpk09WO862+EIpyM01bvuuqtJtZkwp6T1119fm4379+8fWplGCDX716HyL3PGITb38QFPsxtjkkw7AOB5Vh4v5i+D564X44v4zmQmVJMcf4MNNmiiaQafBwFOnjxZ3XDDDVrTNJot\/2IStmud2iExePJuuOGG+nyUPL84KWHGtfuY8BiegVnZbngMi5dvRpaxLvd5sbug0NzH5y2h1qtA4uolitmQ+3ielceLNeo6JyGeu16MryhCRXu8\/vrr1T333KOuu+46hZYYbHi9crZKnVTqnxrz7s0336ydlShMfuihhzbxioVoR4wYoa644gqd5zcPxyHRUF0opX4fnxe7y+ib+\/i8I9S4OMak5kWXRWD6VDTm0Ys16mI1iPDc9WJ8RREq98X5COIbN26czsfL2eeyyy6rvvrqK+2UxNkqBccxxRIOY9rXX3+tk+Xff\/\/9WkNFG6VM2zfffKOooUo6QpyH7r333gXCbsISO5gwlTinJEnskOSr0rSvz4vdZdTNfXxeEarLWScaTpLUdS6LwPSpaFYeb9boiBFK9ekTjngdz11vxhcysswmX0JeSKqA41Ew6X3weSZLkulHSkHIlly\/N910k84wRKMfnr1os6+++qo66KCDnBI7BE3BQeIUDTXJ1yS8r8+L3WX0zX183hCqi4bDYAj0dzkTdZn8BT9Yaa4q\/Bqv1miUhaGO565X4wvMdm6ESpjIlltuGbqYIMinn35a4QlMmkCjSQaT3qPpEnKCCZgY06ikDVGJHahoQ8WbXr16qVNPPbXmBWyEEkLN\/q77vNhdRt\/cx+cNobqewSU5FzVFql3OUEPO9sjCRsOS1shWlTWaGA\/wc8CuKuNLM8eJCZWgYuJGn332Wf08zLNTp07VJt2ohPOmDw5FwXqoQaHjYkzpH9cnSLhBQpVMSWmWyrxrfMjSkn50zX98PswhyRZaWdWm4uZzziOPxGbaqcWjKuXedz4BUwVr4MCB+kjLNEh16NChaquttooTL\/e\/N\/odLBoPMz74ZKmllsodvyJvmJhQEYbzTeJA0ThdG8nicS7CMSmYQQjNFNKlEZ\/KWSuJItZee2116623qqeeekovXpO+MCuh2jJLpiTXGZzXz4csLclG1LR3cx+fD3MIoa7erZvzNM6cMEEnc4jKtAOZftGzp\/p3v376nvWy8gT7ojzwE9VIEkP1qjJbI9doGXiY8ZEFL8+0sGXMUSpCtQX7xz\/+ocNSyNlLLt6wRlYkMh1ROu3aa6+txYhCnnj58jfT7Go0nLXimERe3EceeaR2XVwu3zgNVTIlpV9alc\/Skn5o+srmPj5fxtg6QR5sl0w7c3r3brIy6mXlMX3RSHeJyzGrlP42dTfm5Izrz+XyRq3RsvAw42sxGqo96TgSEdLyySefaM2SRAx46ppGyAykSaajzTffXIfWmIxHxKQSX9qhQwf1q1\/9Shcjhzwh0Y8\/\/rhWBg1tlV2gfd+wxA4mzWAcoYqXr8trG97Hx\/ONJGc9Po4v6Wx6MUbXM9SIOMZEmXYizva23XbbJmbeKJwJ+SOXeVmtUfNXFh6NGl8e85dZQ0WjJDSGNIB8uA488EBtrsX2\/cQTT6izzz5b10HFE\/i4446rmW0JqSHf7qeffqoT3xM2c+GFF+oKM9zr2GOP1cXI+T0xrqj\/kGi9xA4GECHUPJaG\/4TKjprMWcGzLz5+URqFzy+z66x7MUYXL98C4xj5lrVPEI7jUufZdX7i+jVi\/srEoxHji8Pc9e+ZCdU8iOT3aKhUfVl++eV1HCq\/w1lp++231\/GktldvWOrBO+64Q40aNUoNGzZMEy9nq3jrcpZKfGtcWE7coMXLNw6h+L\/7stgh0gEDBkQOCFJFswg2HOywunT2MI9o\/OzN6+HLHKoGxjEmJZBp06aV5v3biPkrE49GjM\/13YnrlwuhmjJsLKp6jbjTG2+8URNuGKHiPYbWyscOjZSGVvr666\/rtIOQc1wzpEk\/O+0g\/y+EGode\/N99WOxopJin4trYsWNrmmoabZb7JzEnx8lT1t99mMMaFg2MY0yygW\/uGirzURYeXq3PwEubmVBNGbaRI0eqX\/7yl+rFF1\/UFQLw6qVKjKkms8466+gJQUvlLNSYfIlJxdSL5+9rr72mi4uTdcm4o0OoEKEh4nofHdv7l344SZEnGAIXQs3nc+3DYk961pNGm01LwPnMQra7+DCHoSMsOY7RdR1xfMDmLNiK2mw1av6y4uG6ahs1Plf56vXLTKiQGOemFPTFMQnHIwhyzTXX1CkFITUAIpcvoTYmqQNCkcsXwiOJw5577lnTWlmgJ598sva4pEg417toqBAvBAz58pxgTVSjofbo0UOnRyyiscHgHJgz5OZYD9WH8bEpc22YheuZhs196MOaps2YMaPuNXZfVznK7OfDHGbBI6\/xjR8\/XucTj2scHeDbYTcUjODZPf3MGoq7Z72\/5zW+pDJkwSPJs8z4NtpoI7XzzjsnubThfTMTKmEzePYSO4QzUpA0MQefe+65atKkSWqHHXZoUlkGT16yG73wwguaiCFNFiIVY9BoOUMlpSGZk3BoijM5QJicweI5TINQu3btukC1mTjTdMNnRQQQBAQBQaCFI\/DFF1+0jMQO9jwTFkNIC+Xb8NYN08rQNNEa2eFAlsYEy30w\/eKMRJIIHJAgUUxwhNnQdt11V+3sZJI61FtjLoS67777KvIPSxMEBAFBQBCoLgI4B7a4xA7EjJKGkFAX0nN169YtVJPk0J4zBQCyzb5h00lfdic0wmbiNFNzD1eTL1o0XsdFNMzTbArwcm7VqlURj2joPX0YX9DcFgUY6ePMOZcLqJjsXEyA5l5lxia6yG\/6NGIOyzIXMsZGjI\/nogi4rCeOtMjQlrblOb4y58V1vGZ8FEdpEakHbWAw9RKHynkBWih274033lj95Cc\/0Y5GkKNJyBBMhu8KsGs\/cUpyRSp9Px8cBly8fCFTCM\/FG9igxVFBktjEMkMpksxoI+bQ1aEljyQJjRhf0rCSLF7BScZXry47a6bMeXFdo0nG53rPsvplPkO165LGCW2HzYT1nTVrlk7jhYbHzgTnodVXX91ZQ+WedtiMlG+Lm5Hkf\/dlsaNJRjknQaZ8uDlKcP2gGE9OV2sJyGb5aCafGfcryp7DMsnGaKiYC4k0KMsxMOkYs2y2XOYvri47OCWVuaz17DI+99Vebs\/MhIrJl\/NPHIEeeughvYjRWldbbTXt7UssKE5HhNTQB80Vz17SDK6\/\/vr63BQzMeXdyI5kTL3AgMmU3\/GTh\/lU4lCzLy6fFntUaAtEapI6JNFmIdWkBGwjHqctZJ8dtzuUPYdJP9xZyKZRhMpzy9psxc2fS132eVXvkmWDyjovbqvTo8QjIQPKTKjck\/AX8vVCoPvss48i+byJHV1yySXVnXfeqf+OxkA6QtvLljNSwmbQJiFmzhbw8OVeODzRrr\/+ep183262ZkyNVcjanM3iBBUMmbG1V8nl67q0F+wX9zKnv3OxV\/LxiKpj6arNImFSAp53Dedr8\/41LSJrXrEgzL97I+awLLJpJKFm2Wwlmfh68+eSsZFnETZLGuQy58V1jI1Yn66yxfXLTKjBBA2LLrqoJkO0Us6ocM5BiyV05m9\/+5uO4eO\/qVBDvBHlgEg1CPESNkM8l0l8T8zWQQcdpEu+4UFsa6lRyfFtohWTb9z0J\/+7z4u93mhdtFlzfRICdtUWks9E+isaMYdlkU0jCTXNZivNLNabP9eaAmio5PMvc15cx9qI9ekqW1y\/zIRqCGzTTTfVZl4ckzi\/QFu866671E033aSzI2EGxiS8xhpr6N+tuOKKWjayKaG5kq8XQoVk8cLFaxgHJ8JcMANzrUk9GFW+jZjWq666SvXq1UubkrmHqUBTloaKJoSXKZp2lEYUNylV\/ntLGB\/e6hR1qDd\/LgScVFsoa94bMYdlkQ0YNmJ8aTZbaec7anwcKSTI56++\/z6dxSWt3K7XNXL+XGWM6peZUKlpethhh+l4UhKLE45CBRlMCbfddpv64x\/\/qBPmr7TSStrZiJSCpBY05llDjhAtZluTHQkzMPXwqFKDVguhQti0rAXG+\/XrV8uURBhPnrFOkydP1uXs2DTYZJ51oqpyvYxvwZmIMifvskvrJmbeqDk84IC5atiwuaVNcaPmEM2+b9++oePkXT\/ggAOcslbFAdWo8Rm52DyQLev555+vicr42KSFFWSIG0\/w71HjmzWrlVpvPfdQvSlT5iiOHsqaF9dxmvFxlGhS0Lpe2+h+mQkVjZIXgUxIpAnkhSHRAxmSHnjgAU2KQ4YM0YXFl1hiCW1iuPTSSxtKqDboaJIQYF6NknOkTbz88ssVXs3Nrcn43GaUj1u3bqu7dVZKvftu\/cISzjdy6NjIOcS3gmOeINmwyeW7kUdr5PiC8hO5wDcwz1ZvfGuv3d75UfaaK2NeXAUz40MpKbNwu6t89fplJlQ006OOOkqfmXbs2FGXWYMwSUmIhopjEZ6+OBwRl4pHL16+RkPl97vttpuumYqHMCZbHJgw10JIaKiYkAmnMRmWoky+XEOLq4fKi2sWOdopRc3zarxA5513nvZMxgze3JqMz21G585tp3bccV23zorzrHFaWyijVWUOiyAb8KvK+Iqay3rjO\/LIddXUqfEEvu66s9TQoVNDRSxqXlzxMOPDuTWvTZbrs7P2y0yohrw4E4VUOe+045Uw\/RIqQzJ6KtFsvfXWdU2+EK1xONpxxx11nmBiyTiPtVMWRjkl1SNUJgrCt3fHWQGU6wWBKASmTXvXGZz27dd27isdBYEoBObM6aI++GBUXYBatZql2rQ5TbVu\/aNJuoqIYrXgx6eWmVA5Q\/3d736ny7QdccQR+qwTr10qrqCx4qFLWMvQoUPVfffdpzp06KCOOeYYhTcwDQLmnJUz1iOPPFL\/nt+hxVKgHK2Vv+HsxD3wAuZs0vbmRUM22mk9QjW7V4hVmiBQNAJnnrmF0xkq4Qu\/\/\/1zRYsj928hCDz5ZDt1zjnhWipWEDx8e\/So\/nrDipi3ubzoJZCZUDlDJUyGfL433HCD9owM1jDlPBXSxSnJOCwlGVhchqUk95K+gkBZCLh4+TYyHrUsHOQ55SOQoi57+UI2wydmJlQweemll7T3GqEvOCY9+eSTOr8vjjlvvfWWJliIlKQPaZIdF50DuBnOqwypIghQTjOqPKvRFvr3r4iwIkazRMCxLnuzHHvZg8qFUDkz\/ctf\/qLOOeecJqkDzWA4Q8Uxaa+99kqUmaNsMOR5gkARCIi2UASqck9BoHoI5EKoZlizZ8\/WxEqeXrJd4EyEExKORb7Zwqs3VSJRc0BAtIXmMIsyBkEgHIFcCbUlgExeYgL5jRMU5myKptPQwkmpSHP5PYkqSM+43nrrVQI6xnbmmWdqWexza7uCjz3GpL9v9CDzGp89t1Waw6TjYz7C8l5Hzatv8xc1vqrOH\/ImmUMzdw8++KCeGts5sznMYb3xVXUOhVATfCXMJJqFG1XQ3BAqnsmEAZlE\/ZBUWNL+BCIU1pUUkBdffLH2riY8ibESI3z66adrUz6Zq2imD\/\/Nebnr7+2Qp8IGUefGeY3Pns8qZcJKOj7mIyzvdVxN4UbMHc\/Ma3xRhTMaNS77uUnHOGbMmNrm3sxl7969FeGGYe+mb+9g1PiweFb1O5oLoeLpS1YkEtjblWSCi9Rnb112jiTrp8I9DQ01TFs1+V+DWiy\/txd6VbTSqA8Jm4VRo0apvffeW1155ZU6bMkmk6hNQ9Tvq0Q+yJh2fJQcNB+rKs9h3Ph4F8PyXkdtEn2bv6jxBdOWVoFI076DwTmxi4rw38F31rc5jBpflb+jmQmVOFTiR\/HqpapM+\/bta0ns7YXy0Ucf6XPVe+65R3sD+9pYqIZQwz4+mABNxRx7QZvf\/\/a3v60NPRg\/WyVM7JcTYiXPMo2dYdeuXfXmIsnvjSm8KmPMMj4f5jBufGY+glnFDBEH59u3+as3Ph\/mz2xM2YhHvWv2nNgbBSIqwt5NX+cQLILjq+ocZiZUkzqQQVMujWQOwYYGe8EFF6hXX31V75oabXrI8lG3CdUses5QIUzyTlLfFe3VmIeDvzfPNqYnyKlqC93WvKM+sD4Tapbx2XNV1Tl0GZ\/PhJplfPa7X9X5Q8YkYwyasX3YFGUZX5XnMDGhQo7El1L0mxZMM1iPrIKFwLMQW6OuDRKqLYfRCoIEWe\/3Rttt1HiCzw3KmvSc2Gwyqmpuyjq+MDNUlebQdXxmHGEaapXNhVnHF7beqzR\/tmZqviP1zPBhJuyqm+2TzKGLib7eN7ns72piQkVAktX3799fESbD7ojDdBLfk1YwqhE2QxUWTMI+t6DJ15hWSGeIAw+FAWyTS\/D3nMGiwUYl8G8kNozNpHY0ckQ5qfB3n5ySzIcq6\/jI+lXVOUwyf1GFJqrqlJR0\/qLGB9lUdf6SjpH+559\/vv4W21a\/5jKHUeOr8hymIlT7ow95UC1mkUUW0TGoyy23XCM5ofBnR5l8eTAmb7Pzt926o35fpTNU283egGgsCpRTMmcW9ljsa1x+X\/jk1HlAnuOz57Yqc5hmfMAVtrGLmlff5i9qfFWcP2RNOof2OMzcmLC25jCH9cZX1TnMTKg4IlGqjYK6HKCT+J6E9gsvvHAj3z95tiAgCAgCgoAgUCoCmQkVL1\/qf77xxht1Q2YYFSEHFI3FUUeaICAICAKCgCDQnBDIhVAJ\/Mcs+PHHHyuclqIaSe5Hjx7tddhMc5p8GYsgIAgIAoJAfghkJlREMemySIJPjCl1TJdeemldx3Tq1Knq22+\/1fVSqUiz0047aQcmaYKAICAICAKCQHNCIDOh\/ve\/\/9Xp5958802d0xYSxQRMI\/b0u+++U9dee6166qmndL3Utm3bNif8ZCyCgCAgCAgCgoBGIDOhGi\/Bzp07q+OOO04988wz6pZbbtFnqtQ\/\/eUvf6k23HBDHSpC4gMKjUsTBAQBQUAQEASaGwK5ESqa5\/PPP69jMINthRVWUJtttpnWVkkKkabIeHMDXsYjCAgCgoAg0LwQyEyoX375pdY60Ug5Q+3bt69aZZVVtEY6cOBA9c9\/\/lONHDlSYRrGu5d4RZ9TDzav6ZfRCAKCgCAgCOSFQGZCxeFo\/\/3318nxCbalOgnZgfr06aN\/f9hhh6nXXntNUVYIRyWq0rRp0yYv+eU+goAgIAgIAoJAJRDITKiff\/65Ovzww9Vbb72l0xDuvPPO2imJkl9kTtpzzz3VJ598oh555BHtkPTnP\/9ZNNRKTL0IIQgIAt4jsO22Sk2fHj2MtdZSauxY74fpywAyE6pxSvrZz36mSPJASAxFqGfOnKmdlAiboS277LJq5ZVXVrfddpsQqi+rQ+QUBASBaiNAbvQ4QoVU4\/owyrg+QsyxayEzoXKGStJ7CPPCCy\/U56gLLbSQ\/kFjJZzmm2++UXfddZf64osvxCkpdkqkgyAgCLR4BFw1TxdCdSFLlz7TprX4aYkDIDOhUr7t6quvVsOGDVNXXXWV6tatmyZTGlmTMAm\/8sorql+\/fuqII47QWqv5e5xw8ndBQBAQBFokAi5ECcG59HMhS5c+PM+V6FvkpOUQhwpu7733ns6CRIYkEuOfdNJJ2rz74YcfqiuuuEKXe+P8lPhUEujHNWNGJp1hVA3VsD7c94wzzlAPPvigfkSwEgjOUqQ+7Nmzp6KcXCMamww0dQqRt2rVqhEipHqmr3KbjZ1gnmraU1\/k63qpjNwuRNkIQnWVK8HKqQzmCWSO6ppZQ+XG1IhE+4RY7YLj5qGUdoNo0WJJkB\/XTAHaPfbYQxMk3sJhhZ0hZ7vPjBkz1PTp05vUG8W72C7USxkyu9RYnCx5\/33OnDl6gwEerVu3zvv2hd3PV7kBxFfZfZVbMM\/hNXQlLpd+LtqnSx9XAk9oGvZ5nQdnOjOhsruAAAmHQRtF8yP21LQ111xTffXVV+rYY49V2223nSa7eppZsD4jeYINSZp7uvShb1hleAgV83OXLl307YiZ5aesxuJBc0eDX3zxxct6bObn+Cq3+bgL5pmXQKIb+LpeqiJ3644d6zoJzW3XTs195x3l0k9\/b+s4HHEvlz7Oz1t33djnzR0zpraeojBHJp+seAwoM6Hi2UtiB5Lfn3nmmaHno2itFL596aWXdD7fekXI7Wrz6623nk68TzrDSy65pEZAafow2LACvocccog6+OCDE30ssnQmwcXs2bN1LK5PRQJ8lZu58lV2X+UWzLN8IeZdu3q3bqrVrFmRN4IEZ06Y4NSPm8Tdy6VPns\/jXqZFrXMcXetxRXaU879DZkI12iIaH9ofsad\/\/etftZkNkybVZXbffXddBxVCu\/HGG+uGzaQhyyDp4l0cZio2hDpo0KBakv6yNVRk++ijj7RW7JPJ11e5eWV8ld1XuQXzeR9qLGs0F7+R4Kc9TvNUa62l5kyZolrvsktNGzSlM0O1uriQmHkCRzOMeV6M5oxc8wcfey\/TIWqdt0gNFTBOP\/10fX9oHEkAACAASURBVI46bdq00Hqo5O5t3769WmONNbSmueSSS0aC7WLOrdcnSMj2gwyhyhlq8p2Zz+ccvsruq9ysLl9lhwTJ+oZ2lHbDiy\/HiBEj1Lhx42ovGqQ6YMAAxRGYU+vTpz7BcZNAXCiEShKdBWR3uVdcrKp5XsZ7tZs7d55DqHXO6utaCZvHzBoqN+Vc9O6779ZA9e\/fX3Xt2lWbZ4lRJWE+mZNwWCIdIX+Pa2mdknB4Ov\/88\/UzwvIFC6HGIR\/9d58Xva+y+yq3r4RKFMBpp52mv1nSikGgy5w5atCSS84jVUsjZjPQRLP2NMNTZkJFQz311FNrlWbI18tubOGFF1b\/+9\/\/tNnjs88+07um7bffXp+lLrroonVnyw6JsUNfINptttlGe\/yG9eHv1GS1G88TL9\/sL4d83LNjmPQOgnlSxLL1DzsSynZHudpGgI3K4MGD1e2LLaa2mDMn1sRsa7G+IJmZUA2xbbrppjoUhLNSwkJM43ckyOd3L774YuwZapHAiYaaHl35uKfHLu2Vgnla5NJdV4XvQzrJ\/biqhq8QavSEkQmJ7EcbbLBBzcsXUy\/pBtFETe1TtMfJkyfHevlmXRq2J6+tnXLfKrwwvn4kfZXbV\/Ojz3L7KnuR3wesmyNHzlPKsGYecsi8f6vcwGPUqFFNIixc5Q37DguhOqBHSMz111+v7rnnHnXdddeFJm4gBvD444\/XcahHH310YakHbYckRCdJ\/2WXXVY7Ty3yhXGASnfxlZh8lVswd12Z+fbzbb3gP\/SHPzyn3nor38QvhkgHDFgQX37n4FKS78QkuBspZfl2E6URTKxT7zZ8h40vC7401MUmPPHTTz\/VkSBi8o1AjzPSr7\/+WmdHolwbXm0c6nfq1Em98MILatKkSeq7777Tpd2Iu8RZaIkllqjdjZy+yyyzjD5vzaNBmGjChOYwkcHQGf5OvdYNN9ywYakHGSdp8PAgpJCAT81XuQXzxqwyn9bLiBFKv5OrrvpYrpnUzj9f\/eDdG43\/oYcqNXx4vvNjwgZJwbraaqup4fMfYBQMnnbKKafooiYPP\/ywev\/993W61l\/84he1I7l33nlHPfbYY6pHjx76X6yQYc32ZYEsCQcMXkNYo\/Fupk+HDh3UKo8+GjvoQ8aOVd27d4\/tV6UOqc9QX331VXXRRRdpMiQTEkXEkzZ7ApNeG9bfNlHwdwgVj2PbKYlUhSRWkCYICAKCgI3APEJdNTdCRTslM2BcI\/olT96AwGh892xtEQ3x3nvv1Rt6NEY8bW2lw77OkCCaKZoqJEmiHQh42223rSXZsR1FTVY7lBabhM19IVUIFX8aLBhxbejQoTqlrU8tFaHiuUt2JP7dZZdd9PkpAKGx4smFZrrxxhvr7ElkUpowYYLOVvOb3\/xGrbjiijV80E733Xff2jlrVuCEULMiKNcLAi0XgTBChRTr5TuohxZnpmi+cQ0tlTPVNI1z2OBZbDDawWipJvOcIVs0WSqEHXnkkfpYzGilkN7hhx+uKE5imvFHgRxRSMhgxLcf0gteH6ehuhIqmjVFV3xqqQj15ptv1toptvFevXrVzkQBiolg58PZJRVVaBMnTlR9+\/bVJFxk+TYXk68kx0+3PH07E7NH6avsvsoN9j7IHtQgW7d+Tq26atMz1DiTbbq3Kb+rws5hbU3TfhLfx\/Hjx+vEOpBelIZqrjGWPVvLffTRRzWRUlkMBeqpp56qhTIaDRWCrXuGevvtaoszz4zP8JQwyX5+qKa\/U2JCRdM866yzdFYk8vKusMIK+gwV0y8TRl5czAiEyphmPIE5M7388sv1NUU0cUoqAtV59\/ThAxk1el9l91VuX9aLC6H6qKHaZ6jMBUdrKDQ4jeJfQgs7Q6VUJkl4KHKCSRiN1jRD0vjA4FxKoh4SYWCJNNqsOUNFaarr5QuhbrFF7d4+r\/PgNycxoZpDaJNGkMm79NJLde3Rhx56SCdWQIPF3GsaNnuIdurUqer+++9X66yzTmFffnsigykGxcs3Pew+L3pfZfdVbl8IFTkXWujHdyJMQ03\/xiiF5zD1uONa3meocc8zfw+afF2uI6IDh08Uqz333FM7PJlkO1Gasblv1PfX53VeCKGyQ5k5c6bO10sMKskdONyn4XaN+s8hNaaGP\/\/5z2rttdeuO3dFFRgXQnV5ZcL7+LzofZXdV7l9IlQIz6TczZtQwcG+f9ibVYSXr+tbnoZQg\/fm3JWUsngKxzmZCqGGzIwx33KAj\/kWkiSZNCYEnJQwE9x22236kPuJJ55QZ599di31YMeOHdU111wT64SUNpevFBh3fZWS95OPe3LMsl4hmGdFMP56W4ssglCRIOoctupxqPHoJeshhBqCF0mML7zwQsXhNBUVIEka5l4Kd9NIUk8+Xwi1bdu26sADD9RxqgcddJCOU+UsNaplrTZj3zeqwLhUm0n2IvikcYSNzFdi8lVu39YLnrgUUSmKUMEjmClpm23yDZVJ\/kaXf4UQagTmxmuXBA5UmoFk8folCX7Pnj21hy9thx12UCeeeKK64IIL1Ouvv64Pxbfaaqu6M5lHPVQeEFbGzUwoxE\/9VlrZ9VD5SGIGX3nllWuxXOUv7eRP9FVu83EXzJPPeZYrfFsvaKrHHz9ZfffdwbnFoWbBrzlea76\/HAkGnZLC3s8WUw8VAh0yZIj+MQ0vXwh07Nix8zOOrKprpC6yyCKacE844QSdfjC0+K21evIg1LgC4\/ZixZsNh6myWlR1+rKen\/Y5vsrNeH2V3Ve5fcWcuEuyBzXSgpX2\/fThOkOoHBVy3mpa1DqnQhnxrj61xF6+ZnDEIJHAgUQOH330kSJzEtmSOnfurE3ChMZg4p0yZYomVVyp8QTGQ8xuwTRZlPdB20XzZRdjYpv4f9OyFhgfNGiQNkU3QkONqk5f9UXjq9zg6qvsvsrtK+Z88An3E0Kd9zWqlxzfZFKCGFGidt1119hPmCFUokA222yzWv+odd5iNFQbOUNub7zxhtZA8fiFNIlNJVwG719MvpAvGqyd8CFqBtI6JUmB8dg1nbqDnOelhi71hYJ5auhSXVhIFAC25PHjo+Wp8GGqS3J8k10pKtevPXA5Q3VYlhDqmWeeqavIEOQb5nDEDoSAYrJ0EKeKB3C9JgXGHYAvuYt83EsGXJJplA54IYQal2qpAFffMpLjo6Hif0JuYBLrB0tlhk2eEKrDksb0y0\/c2SgaK2E1TELQ7OvwmFy6FPLCJJTMV2LyVW6mx1fZfZXbV8wL+T40gFCLTo7PuuSHfAIc7dWrRiMaakKCCHY3pl7Kti299NKabJOWKUub2MEmahbVM88806Q4biEvTEK8fP1I+iq3rx93n+X2VfbQ70OW3IMAEZchP0tmfO4fkh2\/qOT422+\/vU43u9FGG+kEPnkR6tSpc3UFms6d22iHVp9baqek4KAxM1DdHdMu4HBYjQ0ejy6qynC+6uqxlfYM1bhim+wdm266qRBqTqtTCDUnIBPcRjBPAFYOXUMJNU7DzOG5mW4RYjIuMjk+SXw42hszZoyucZpFQ+V4GXhNpiqzP6B8a57l7DLhm\/DiXAiV3QrFw++77z7tPcuHgEoGkOm5556rNcXNN99clwpaaaWV6oqYNbEDxE4VHPIFM\/mXXHJJzcQsGmrC1WF1l497euzSXimYp0Uu3XXNRUMtMjl+mzZtFiBUUg+SWN+OxAibARvfxx7bom7hdUjVs8ptesi5EKop50a4zD777KMISwE8tFVSE7JjIhTm1FNP1SXc6rWscagPPPCAvj27JzTmMEKVxA7JPzi+BerbI\/RVdl\/lNiZf35JpTJ48Wcek5xo2E6fhFuCU5Pp2p8nlS7k2NNQ4x9J6hNqhw+0\/ZNr7sdpMlLyNKhrgil9Yv8yE+tVXX+lybZxfQqiLLbaYwmRrCBXgCdxFU\/33v\/+ttdallloqUuYshEp1d8ib\/MEEaUcRqv1wSezgtnwkyYAbTnn2EszzRDP+XoUkdmhmhBqPYnQPo6F+8MHtas6ceEJtZOGAtOPMTKjGRMv5pVH5g4SKcGG\/4\/d5JnYgfzAhPHbDFEGyCJqZ0EYmdsDTeeTIkbreIPL60nyVG3x9ld1XuX3FvJDEDh7Hoeb9bUpKqDz\/++\/zlqLY+2Um1H\/961\/ajEv9U8iMONQgeeL5S5zSSy+9pL3E4pyTsjolGfKM0lBzNekknJ8qnOMmFLnJZqSR2KWR295I+Sa7r2vFV8x9xjvtu1HmdWkIddq0eY7MvrTMhGqqzzz77LOaLNG6goRKTl8qxpMYH9NvXMxq2sQONuhhabPCkuOXPVFUuqfijn2OW7YMaZ7nq9yM1VfZfZXbV8wN3r5tvNK8z424Jg2htjgNlYlB8zz00EN1XBJ1UZ988kntYct56VtvvaUJFpKk3BuabKParFmzNJmRf1iaICAICAJhCPi22fVlFs2GxfUMldAZHJN8apk1VAaLSfcvf\/mLOuecc3T+3mAjsQMm37322qtuLdQygINU+ZEmCAgCgoCNAB\/8e+65RzbcBS6LOXO6qNmzB6m5c9vVfQpmXh\/jUXMhVIPM7NmzNbE+\/vjjOhaVrBdbb7212mOPPXRcqjRBQBAQBKqMgGy4s83OjBn140chUrgA7ZPC7mENMsXDt3\/\/bLI04upcCbURA5BnCgKCgCAgCFQHAYiyT594sozKlASR+pjUgRELoVZnHYokgoAgIAg0CwSSkiUpk33y5o2aJCHUZrF8ZRCCgCAgCFQTgeZCli7oCqG6oJSgD5maCMi3k1xQA5Zm1wy0K0LYvzeu5cH+CURI3TWr7PaYVlttNTV8+HC13nrrpZYn7kLkNYk8KMZAqksyc0VhmPT3cc\/P8ve8ZK865mBkkrfsv\/\/+yhSwaMQ6T4p5lOxVxjyYx\/eoo46qfYuqjnk92cvGPO27LYSaFrmQ68ykm0XMAuZ3fOhJzUiKRj4qtLDfr7\/++jrs6KyzztJ9Lr74YnXZZZelypuZdFhZZYfQzPjMRzOpDEn6U1HIxgf533\/\/fXX66adrb\/Mghtw7DNuo36fJVeoqf16y22uqipiDoR1TbuI77fSiZa1z8gpTCJv3j9SnOE+SCnW33XbTnr1kVKPxe\/qQg5y0qsOGDVMzZ87UcfQU3Pjmm2\/U3XffrTcG\/H\/RDbltmR566CEt96677qpGjx6tdt999yZyv\/baa+qTTz7R4zLyIysl10h0E+zPOItqWWU\/6aSTFCXjcGot89uSBQ8h1CzoWdey+yUh\/\/jx4\/Vv0VDDND6TbjCoxfJ7rg8j2qI\/lnnIvuOOO9YIq0itNGq6TCKPvffeW1155ZXOm5iozU3RmNvjSCu7vQGrIuZssqgw1atXL10Yg3cCXKM2mmViTqgfNZsXXnhhXbN5kUUW0SF9\/Df\/mv\/m76Yfv6PZfXL6fDjfBrn5MfIhdz2ZjKymT1x\/Z0FSdEwqO5sHKpk1+tuSZKhCqEnQcugLIRpCDftwYArdZpttFiBO83uTLpF7sCvr2rWr2m+\/\/RyenL1LVtl\/+9vf1oSwTU3ZJYu\/g0lXaVcZsjFM+vuyMDeEbjZUYfNfT\/YqY24wDJZktLOYNWKd80xI0hAn\/22yt4WRpU2okAJ9TOMekG5ZzSZIZDHPjiL5pP2LHEdSWbAO8K389NNPVSPXeRJMhFCToOXQ1yYl87HkDBXC7N69u\/rpT3+qd+rGxGr\/nsVTFUJNKrtdC9GchZS1GbAtAVEf66oSahbZbdKvIuZVIVQ0HdaiySHOhx0yQlszmmoYoVJXmfe1Q4cOmriMhmo0RKNx8a7wd+5BmtV6JSqZJ+pDf\/7552rZZZdVU6dO1TH7mGd\/+ctfqkUXXVQtscQSavXVV1errrpqky+OeR7P4r+pjsPxASZgYv5XXnllRSrYjz\/+WH399ddqgw02UD\/72c\/UKqusovvHEfArr7yiM9txP2SgdvW3335bux\/36tixo35O0sazyZ6HzP\/4xz9U+\/btFXVUP\/jgA\/Xqq69qubfbbrta7WrmaNq0aQuY1cte50nHWQqhUoaKyaL9\/Oc\/1yXemmsLEqo9TqNFBbUfW7tqhMnXyJhF9rAxGU29yLkOYpr03NpsHILn3GWYH7PKHpSx3vzlOQeuchv5wjTUMtc51abatm2rTc9GMzWanSHUMJPvgAEDVJcuXdROO+3UhFBtLI22au43btw4vXGOapzzX3vttdqUCcFzJtujRw99TssPRAiZkfd8yy23rJF4UAPlOjbqFCdhs8CGcd1119X9\/\/a3v6k33nhDE9Lmm2+u5acZ07b5b6Oh8\/8Q59ChQ\/X5K\/eD8MwZMfd788031bvvvqs6d+6s72c2Fy7riudCmOR6h6z5\/kOMhlDvuusu\/ftLL71U+4vEmdTLWucuYwv2KYVQzQvFw40nZhphfbgmaDY1GicZWHCYYdFgwgj7\/QorrNAwpyRDLIYEbU3PRXbGxPkxmmrwA1rUvIE1Wr1NLFEOL8hQFackg3VW2dGGqoy5cewKroeynJJIg4r2d9NNN2mC4P8hBDb1ptmEyu9sLe78889Xm266qdpll11qhBrU9ILnq0888YTaYYcdIpf8Cy+8oChkTpHuMWPGaCK74IILtGaLsxNa6qRJk7TjFDKDIevEJj+0uokTJ2ripb\/RQP\/6179qrRZiheDJsU7ecgiS36HF2lq22QSY+\/FMcq1zv7AWvB\/Pj2uQI5opP2j7nKtD1DNmzNA+JmRNQvOFtHv27Km19qAJnW9R2es8blxRfy+FUNmNsQuh7bvvvnULjKcdSFWuizL5Ip9dxcJ2A7d\/b7u2l131Iqvs9piKPkO1cTJzj6fmJZdcos1K5szFBduyMc9T9qpjztyEbbDKwPyRRx7RWhXmVMywbGL5Fh199NFaSzJmWkgSYuP\/jaaK3GGEarQ7rqEFz1AhSZxoohrEgYmzW7du2ksXzey8887TplZIEq\/W+++\/X5tCyYEOCRmPdXNPTNj\/\/Oc\/tTcyBERDdvptttlmas8996w9fuzYseqjjz7SmiqkamvUhlzxIsZJkncmzrud+yEbG5M4RzgwwjqJZgruhx9+uGrTpo2Wlb\/dcccden74Pdo0xAqpm8bmnOdg\/SpznWfhklIINYuAcq0gIAgIAkkQePTRR7XJE40HjQ1ChygwO3J+Sc5xtE5DRlH37t+\/v9ZQCdtwbVEaqjl3xZJz3XXXaa0ToodIDjvsME2omFU7deqkiRayg5gZB5W87IazIv1QTiBd09B0ISDGZhpjx9MaLdWYfoNjibpf2Jhd7meuQ6NlA4G5Hc00eCZMP8h55MiR+m\/Ml02orphXqZ8Qao6zwe6PFwrzBIvaPgPJ8TFyK0FAEIhAAE2TsCkciTjLxOnHkCsaIASGSRanH5uMgreDfC+66CJNUMR0cubHT1iDLDF\/EhccpqGipXE2CUliIoZQKRrCEQ8kigYGoXK0wlkoplDMvXw\/IKTevXs3eSzl5bbddlutiWLS5twVjXvQoEFqww031OPjW8QmAtnYGHAWbGJQg2MIu1+wD\/dxvR\/XMg9oppjZjz32WK1pRzWsl2jIOCVhDjbOYz4u8twIlQlkYgkUZuGwILHpYzvHWwyCiTMn+AigLTPmC8bLQr\/++ut1sDgvrjRBQBAoHgE++pxFYl6EmG655RZFOBJepZgZMSliEl1xxRXrCoM5esqUKfq7hWbFOR+Ex73DnHEgD+4NiYdpqFdffbUm76WXXlp79+KZS6IFiBptGVPsZ599psmR7+djjz2mPXT5Gx67Bx54YBN5TzjhBLXzzjvr5A7IhGx4GD\/99NNaZupSIxN9IHlMymwMTPKK4ODD7hccZ5L7gR8bBbCH4DHzxjXGeuutt+qNCcdFvrZcCBUzAJPG+QNnWEzivffeq+MoIVYaNnxMD0y4L812zGFMNDt9WdgZJ+N98cUXFc4HeLH5bsLwZa5ETkGA80gcetDsJkyYoLVBzKbXXHON9qY9+eST62qlBsE\/\/vGPmvg4ayWsLUkL01BxRIRcMPMS7UAfSO+II46ohe7Yjjhoa5D+ww8\/rM3CfEfsxjhwxMNEahI18Pcwky\/X803+1a9+pQk4rEXdL6yvy\/047yT8B69qTLkmJKkejswP5ng2MmSl8rXlQqg333yzNo9gmsCTFW31uOOOUxTsJU8tuw8OzI888kj9ex+aIU7j6AKh2inj2Bni5GA89DDrsKtEC+echOvRystIT+YDniKjIFA0AphzMetiFcKL9JBDDtGPhASSpNgzGhabfzyCs2qoaJ8QPF62eNESWoOJlnAZiCeqQcSYhXHasRuaHGPCVGyfA5PXGsVlr732qnXHo\/j111\/XmjnaYliLul9YX5f7oSkTr8rGBM0YecwZcvCe5vf33Xef3vywQTDKS9HrpYj7ZyZUJhZNFJMAAKKRcRh98MEHK9LAnXvuuVruCy+8UB9AX3755ZXQ2mzXfV4YOysRLxQvJI1dk9G6IUnz\/5ynEIbBdbx4LAw2DuCAOYcXGxN4lCNAEZMp9xQEWiICOCGxoSUmE62UjzLvN6ZU\/h9zK\/8aD1dIjDjPZZZZJhIuvmF406JVct4JCeLhGpYVCVMt5MiZbZSXL57nDzzwgL4fWhv3QWu1NcygMHgZsyEPmnzxEsYKhsMUYzQExNkxzkp4CtPw7n3qqae0NsxmP2q8UfcLyuN6P3Od8ZI2Z7dsLJgH0+wQGTyNcZzy\/VgwM6Eal3hMECZbDu7QaKRDhgzR5wQ0zAAQVZXiUNE4MacQbB1mvrUJlEVrZ7UJCwXg7IazUw7VIeQTTzzR+wXSEj\/QMma\/EMA6dswxx+hQEs4wIVfMm2iDYQ3nHUJE8D41BBzWD2ceNLuos8ewa+rFoVJ9ie8GssZpYWwOUD6iTLUknscyiOaNExJEZM4qcW5CeaEYAKZrCJkz3npxo8H7oRmbs+Y09wObIKGiXCEvMbj4mhDKZAhUCHX+agoSKrZwzlOffPJJXb4LbQ3tDdMvuz48v6rkxWWqlBgt1H5JkhKquRYM6ply\/PpcibSCQLURgFA5j8RCBInw7kEyfLjjCJWNLw5MYQ1CIHFCEkKtF4eKloqFC20W7bReDmC0QTIXobXh\/RrWqJKDFzMewaQD5FuLVYyKM4wLxyK+RZxhYnrlp17L+35BQjWZpzgKI0YYz2Mh1MCMsBtEG8U9mlJjnB9i88dUYcy77HDw3MIcwZljVcgGwkSb5mwXs4edj5ZhhhFq0OTL2OMCnKv9ORLpBAG\/EUCTRDvju2JnArKT2JsR8ndS\/PGuQ1gQEOX+oggVs2pUuElSDZX+eL9yXgiZRBEl\/Uy6Pwgdc21U44gJokajxfRMw7TLMZS5Dl8PHJ3YZPC7eiEsed4vSKj2\/wc1UtFQrRnGGw3zJs4A7IhwRmI3ss8+++igXdzXMT0Qf8XZRRUamjW7WeO4AHnS7Jy0QUINc0oK02yrMD6RQRBoKQhgPsS7PoxAwzCAcNAAORPFUSi4kTbXBAnBBc+41IP4nOBsxLnmQQcdFKmloulyvoojUZpk9LasmMLx6UDZQeNlk5Clud4vjFCNxs+3Fy4w1kohVGtGWMwE52I+ZZeEbZxzAhpOO2h1mIHZ6SVJqpxl0vO4Nkio3NN4\/+JOj0lbtNM8kC7gHttuq9T06dE3XmuteX+L6zN2bAHCyS0biQB+HDTjEBTMRGRkM9VmOH+0m11thu+Z+X+XajPcB4cijr9wTuJMEa3Q3AczMNomZ6CkFoyLmXXFkW80Jl1CiwgHytpc7mfjBzY4j+G4RXws\/43yRWwujXAiznpbvFNSvYlhoeHwY8wxSScRT1qSQzz++OParOwKtinxwy6IVnRe2aTjkv4lINC+fTxZuhDqtGklCCuPqCICKAhoYy5xlMjP9w6PYo6+6jWcF9Fk8cClBZUM7kOCeqx+LgnoXbEjwQXhO8G4Vtfrg\/3i7mfjZzYdnEdjemb8eD2bhnMYZ+FV8q9Jg0tmL980D427BucCNEHOODibJVtJEu\/gMG9cYmTLLBodN0b5e8EICKEWDLDcvsUjML9AgFpooRYPhQEgV0Ll7BRVnt0Hpl\/ijMg5iVdXXPopTAgEDVNqiZgvk2EJzzSchjAVuO4Ug7Nr12\/EQYof09gJRpUrklXiLwKtO3asq6HObddu3nqqY\/Klz9x33vEXBJFcECgAgVYTJ6pWM2eikje5+9ytt1ZzA6lWW83\/1vIuJW28n2m\/+UmflVf\/XAgV8sP5iMN2zgOCDa8y1PmwM1S0UbzeuJ5zAxr9iRPjDAFHpyznrnYCB2K\/TjvtNF0j0DSckkhCUVYDK1NE2KfF4pvcq3frpszLHDa35gWP6zNzwoSylsYCz\/ENc3sAvsruq9xgn5fs9Uhw2cGD1XKDB0e+E7MHDVJf9uypWj\/3nO7X2vrW8s7NvvRSNWeLLZoSccQ3EYXMNxNwLoTKWQAH3QRSH3\/88Tr+i8N2zLUcNnOmQOJ4KsyTAAJvPBI6Y8bFjs5CgFy4HvduMixFed4l+bqZs1SImeeaGowEdhPUTStbQ0UmPAx5Lh6GvjTf5I7TUJWjU9KcKVMaNkW+YW4D5avsvsqtCXXqVEWChhU32yzdt2XcOPVDeMYCJDh36FClundXatw41doqDRf1YszZbz\/Ven7URFifOdddp5RVki4K8xapoRrSwvmIzEhh5lNMwXj+EgNF5hEyfECoNCo0oCWSLguvO2JY7axLab9mtmZqPHENoZZduNseA6ZwNHE8\/HwiVO\/kbgZnqN5hbi10X2X3Um6I8PzzNeHVGhvG4cPnEaFL4\/oBA6J7cq+RI5s+w+W+UX3wnp8vm5eYR4wrs4Yalnow+CxTk4\/CuRAwaa2IUeWHLCUmY4jLvVzmkPsQA0XAt+0ZLITqgl54n8osepdwGF5WIdT0k53DlZVZLwnH4p3cLkRoFyg3PgPGQgM+EDHvVZkNrqJT9gAAFXNJREFUmSDp+Zo1Skabzp29UjLC4MpMqCRsoIIM2ifVDqLOO0n0QBwUcVVUmsdZCc3RPkfMi1BxQsK8bDdSH0LenM2Khpr8zanMh8aFKAl1sYjXOLgtcGbNR6VeHCowNTBspjKYJ18u2ilRLDEpgEtyiSsRmljqKC02+PskMmTpi1xZNesszy\/g2syEivZJQngyXUBiVG4JNuz6kO4aa6yhsyhRs5CG2ZNg3j322EPXHSSgOS+TbxhWoqGmX0GV+UC6Eqo11Eyyu2rE6aGNvDKT3AXIk+SWyE7yAlKT+nS0webrk08+0REKlZfblQhdNo5JJjdD33Zz5yp+Yhvaq61Zx15QjQ6JCRWHIsrwQKSm4bWKow9B0JQ0wymJMm4EMHNWivcvCZzJmYmzErk0KUB+99131zx7SbFFhQi0WEqe5eGUFIRYCDX9oqvMx71sQk3xvPQoN72yMpinGBDOhXjUkxVImiBgEOgyZ44aNHu2G6la56y+IJiYUI1ZloTMSVswQUNU7ClFcskSQt7HPBPpC6EmnbEf+1fm456C4DLJnuJ56VFuPoQa5lGfFy5yHz8RIFxx8ODB6vYPPlBbzE\/kX3ck1jmrLyNOTKhoneTtRUtN2khKzfkp2muwEY9KyTfiUY0HsIlH5Rrbecn1uVFhM3KG6opgjoTqYjrlcXG5dfl7XJ\/AuacQavL5znpFFTavWccg1+eLQG1NuBIqjw8kj8hXovzvlphQ8xeh6R1NjCrESsIHEkUkTT3IHW1N2hBoFV7yTB\/3osGvc\/\/Mcrtoei6E6tJHCLWBK2Xeowt919hQEcLBv5wPHnLIvH8LaLaDo10QI5gvPPiNCgvbQ7zgdVH3rPfNo+oV+XgJVTQtzTeyALjq3jIVofIuFzS3RYy\/FELlvJVzU1IS7rXXXs7ZL7Ikx7\/qqqtUr1691KmnnqrPY+3EDv369dP1CGllJ3aAmHDSoiQTmZt8aVnljku04JoKME26wCyyu8hdVHrCLHI3el2RRpQMZLlagwyRhsVL8rv+\/XMdNgQwfvz4mj8HRMZG\/+yzz64ldjf5wYN9+X8cNfm+4JBpk6lJNMPv7JKQlGyjcc\/g\/eyBcc1jjz1Wu2+ugy7wZqkIVTTUH2eEM9Knn35aXXvttdqzt+xdlNFSg4Rqr5myUw\/i9Yh5m9zGiy22WIHLN99bZ5U7r1SAc9u2Va1+qLdbrwXTBWaR3UXuotITZpE739lPfjd8LE4++eR8CTUu5jLnMzebQIOb33qEhyLAhp4IhltvvVXnIicevt41IHz11VerHj166HBCvl0UBznhhBMW2Hj7Tqg3L7KI2mrq1PhFReIHz8onFqKhRuXnPeKIIxQ\/EC0J6oNxqHgJE4LD7paQGvpuueWWkUV442YkilAl9WAccgv+PWtKtjhNr8hUgFlkd5G7qPSEWeROPsP5XgF5HHbYYfkRKtopxwZxLWfPUNu8GqyBbDQuRLKVBZvw6DNjxgytdVIFy2igtunX3Bet0yZUk5yG2P6LL75Y49mxY0et1QZNvj6UqKxpqEcdpbY444z6M5k001Pcuijp77kRar1qMcSWbr311jquCzPIBRdcoNZaay11ww031My\/b7\/9turbt6+aFjj\/IjcwP2kSyUcRaq5mqIQTlfksMuHz8uqeWe68zlBTJFrIJHucMxUAp5DJZV4yye3ygAL7hJ6hxjmU1ZOHM9MRI+IlRkvlTDVN4yNe57wuKgMbj7K1WbRTO7EMNUAvueQShdZum5CNiEYzDRKq0VBfe+01fUxERrnddttNP8trk+\/tt6st3nxTqT59wmeJOWAeczbhp1kSSa\/JTKgkwCfBPWZdu1oMHwPMqXb2pClTpqhDDz1UJ8eHPA844ABt9sS0BckSl8q\/nH2SQJ4d2rPPPqtGjBihNt1006RjqzkmBU2+QqiJocye+cZXQk0OVW5XNDtCjTPZ5oZcyhsFzmFtjZI7GlMuJlzIbpttttG+GTRjoiWWnus4ZzVmYhyb6IsWe8YZZ+gqWvZ1xO6fddZZ6m9\/+1tNgzXmYf52xx136Bh9lBHImcQTXhNqhw5qi0cemZfyMCxTEkTqYVIHJi8VoYZ54lKpnsmGDNdZZx2dDD+Y5J6dGDs3\/mWBmYbDEiYM4k4vvPDC2tkiOzO0WzIpQYpJtVTRUFN+WEIuy\/xxF0JNPBmZMU\/8xPwuaC4aajCNqdmMR3nr2qRo0LS1V34HqT744IP6z3FevqR2hUhJloPpmPBBHBqDJt+gOTq\/mczvTrU1sdhiagurihPrfPakSS0zly\/pAfGSnThxokaaAuB483Xr1k3XMTW7tWAKwa+++kovJOJXIdSll166NlOUf+O8lHy7xmuOPxLzimMDgAevSTPNhbryOwrk60cys9wuplOXFGkpzKuZZXec27y7+So3OOT+riXJW+taYSXvCSvgfo8\/\/riuCc3ZKd\/eUaNGaSXDJ4dGA4tNqGeusoou1WmaKeFp\/p8jwbGeOSSl0lCN1oeZlmxGXbt2XSBRQ1iSe\/M78vlynmB7zVEvFbPG8OHD1UYbbVQD2ewCcVaidqpdOSbN2s39JU8hhK8fSV\/lZop8kp2iEVX70ARlCi77sI9fIe8amzK7RFlQkJy9fFO83nJJHQRsQt1\/zpwm6zxsTQX9aXwAN7HJF5IbNmyY\/uH8FDNs9+7dtQeaSRWYhFCN5krSfJyUqEZjWhZCtT3wjOZbyEuecJZ9+rjbQ\/NVbt8IlZzXNqFW4UOTRqbC3rWoc9gC4lATvtrSPQYBIdQ6AJmkC2iVJMBGZcfk27NnT12ajYT4W221VS0oOsp8a85PCY8599xzm5yT4uSEeYMzWbRaEuu7NDtLCf1xOUcLxpNYyre5ILhgHyHUdLglvSoNeSV9RtL+aWQqjFARPpgpCX+MZmTmTTo\/vvQXQnWcqbC4Uy79zW9+o73XSGJAtiTiP2+55ZYmXrscznMmO2TIEO0SbreHH35Ym5UJbjbZRlxEYuJwJsBMjGnZeNZxrRCqC4JCqOlQyn5VGvLK\/tT6d0gjU6GEWvSA5f6FICCEmhBWE4tKCTaIkjy8NByX0DRXWGEFdcwxxyg8go8++mj1zTffqOuuu057uuH9S5ouGuSLt9zpp5+uHZMgRvtsNU4sJo7De7RaGoTKWS8ecngMswGQJggIAsUhQMw5yVkaGaJW3OjkzmkQMISK5RGLV1zD+kmYpU8t8Rmq6+A4X3300UfVTTfdVKseQxzW3nvvrS666KIa2UKyV1xxhTYP02wvYs5nBw4cqENxyBbi2oRQXZGSfoJAMQg0F0KtYnL8sBmzfUbM34NRE1lmmthacgrQ7BAdfo+iYuJq6z0jKaES60v0h0+tMEI1IBCzSvwUJd9eeuklnQDiu+++09k+yPxBYXETbsM1xknp448\/1hmSOFtNQqbco57JF0Lv0KFDw+YILR7N3afE+IDlq9xlyT6unvepIgHPWvonrrncp57TEvc3z4oRSScFqlcJb9696DOurthhYyOzD5EAeWqoYEOmoahGbDsOknm1YO7dqiTHjyJUk+KQv9fL6pQUH0iTNUcuAFowoX9SQuX7i7JVr7Gm0FDznM+k407Tv3BCtYXiowyJ8lNkE6ek\/NEVp6T6mKY5Zwy7o8t9uM7FEzivXBpKJfc8LuIMlcxpA8IqzcwHkr\/1zzFdXVWT49vaoh3BEEWoaHomFaLJ+Wvfw+QhxmmTJBJsXCgNx705IiOVoknwH1yzRkNdf\/319bEe6RVpRi5DxnYqRjZZaLtmDcMLtGDinhYRNpP\/p7qYO9omEKmHmh1jIdRyCNWO+Yz60JRNqGut1TQ2NgyJ4MevORCq0cZMVqJGJce34++DJG\/SGiIrDpd247sH0QWr1hChAXlCpDT7Hsb3hD5sYLASjh49OrTqDdcaQsXxdLnlltO5AoKavNFuuc9pp52mnVJJBGSaz9+W4LtQqoZqHo4Z+Ouvv9ZhMLY5F3MvNQExDePQsPvuu6t27dplZ4P5dyjiJU8qnK+Lx1e5mZ8yZHfRLJPuuKPkdn1WXhpqiuRUoZmS+LDGmavrvU\/UIuVjHNVwYCF\/eNoWZ5ZvRHJ8NEO+iRAg2iJEFyROtEHMrraGamNga6NGQw1WrIF8afY9TIpYoi2CGiqaKGRNHndj8rXPm01BAPvvEyZM0E5GaKs77rijEGrahWpf9+qrr2qnJFJn2ekEOTMlPOaFF16odccbmMVCUfKk56hhsgqhpp\/BMkgpvXT1ryxDdleSSzLG5kaocSbbJNgU0TdoMq5CcnxyonMWjRLy\/PPP62gHU+Tc9sPg2xZFqAYryBgnT76nTzzxRE3rrKehYkKH0INnqERgkGuAaAxD5oZYRUMtYnWG3JNJYYHwL7stFvBSSy2lHV4GDx6srrnmGr2DwRb\/7rvv6kT55P5lV0ouy6xNCDU9gmWQUnrp\/CDUtdaaVtcJCAcgO31pcyNUHzXURifHJ4Mcuc67dOmiq3lBhrbGycoP0y7NG2Gy1pmzTTRUFBc7QT+\/QzuHEO0zVNuhLM7Ll+cZzRnnsLZt2+o87IaM0UqJ6Pj0008XcFTz+dsS\/PKUavK9+eabtXYaDIUxGZEw75LRyHj9koCfMm+QMIkdkmqpJnWhKZckhJqesnxe9GXIHpfvdh7y8YRqm1ebG6GmX33zrozTcPN2Ssoqbx7XY9GjhCWKB0cGp556amGJ8V203KRjssNqor6\/ZbyfSeVO2780QsVsQdYkFkUwZy9AU7GGXRM5gU2jdBFESkzb5ZdfvkAS\/nqDtndm4pSUdnn8eJ3Pi74qsrucaboQqgt58565FPhhhuuFzqA153WGmnUVtkRCzYpZkuuFUJOgFd63NEKtV23G1ElFg91kk01qkqZNjm8KAZMQgh1dlQqMY\/biDATnCZe4xOxTnM8dfJV7HmFUA\/OkhFqG3C6kWxVCLTsONZ83R+5iEIjSUMtY52XNQsMJ1WiuU6ZM0SkG8e41LS2hmuujCoyTO5gzCRrpDk3KwzJAnzx5stbGySDlkl2kDJlcnuGr3IytKrJ37Ni6rjbYrt1c9c4782LyqiS3y\/oI9jGY55nYIY0cck11EDCEGvz2Rb2fxKUGY1OrM5oGa6hR5ltzfvrzn\/9cXXDBBWrRRRetSfrJJ5\/o81NKuiU1+XKTKEK1oejRo4fip6xGBhm84yDVjTfeuKzHZn6Or3Iz8KrIfvzxu6tZs1pFzkWrVrPU3XdPrv29KnKnWTxGdiHUNOg1z2sMoQa\/fVHrfIMNNlD8+NRK01Dx5MVrl5RTttcuqcSIceJv++23XxPsjFPSQQcdpAOCo5yS7JgqE\/+ES3mQUGfNmqXvg\/u5NEGgbARmzpyg5s6NjquGUFdf\/ceA97LlK+J5tjWoiPvLPf1BAI9lvr+ujbXDj0+tNEIFFEOQnTp10ueakCxev2ivmHvtM0XAP+WUU9Trr7+uK9KY5PlJwA0SKtdCqvxIEwTKRmD\/\/beIdQAaNeq5ssUq5Hm8v1Sdks1rIfB6e1OO2lxJkqiPPBP7lAFaqYQKgVL3lB\/TghVlKNeGeZc0VVSsIfUVP2ls6WGEWgao8gxBIAyBohyAqop20ZvXPn3iPZSHD68qOi1TLh9JMslMlUqoCEbGj0mTJumMSJhwSUFFCitjzoVQCQhmZ0sA8oEHHlhY3FUSoKSvICAIVAuBpF7T1ZJepGmOCJROqHEgUlwczZTYU9tBKe46+bsgIAi0LARamsbfsmbXz9FWjlD9hFGkFgQEAUFAEGjpCAih5rwCgsV47VygpkYgj4z6vV12zu6fs5iht8squz2mYKmrIuQPq+lI+agoDJP+vgiZzT3zkr3qmDPeYApQfteIdZ4U8yjZq4y5wfrBBx\/US81Ul\/EB83qyl4152ndfCDUtciHXmUk3i5iPBr\/Dg5kwHlIrklfYEGrw95wl49lMikbaxRdfrHMb2\/UQcxS3ya2yyk5tRTO+MhJWECpl44P8FEWmCsY555yzAIYMNgzbqN8XiXlesttrqoqYg2FYClB+V\/Y6T4p5lOxhm4Oi3knum1RuuzKMwb53797aV6XqmEfJTiGVMr8tWeZTCDULeta1Jgk0cbU0woLCND4TGmSXQ4IMTLWHMAIu+mOZh+z2C7veeuvlhKr7bdi8UBx57733VldeeaXzJiZqc1M05vbI0spub8CqiDmbrKuuukqFpQBtxDpPgnmU7PZmoIqYB9dto78tSTCPkr3R3xb3r5BSQqhJ0HLoywI2hBqmoWIKpbxR8INifg8pXHLJJfoe7Mq6du26QMILBzFSdckqu1342DY1pRIm4UX2hyMMQ0pTJfl9MMlIQnESdc8ie5UxNxiGZSxr5Do3myizia0nS1y2taqtc3vd2uRP2TRfMGd+grI3cp0neZmFUJOg5dDXJiXz4lILEMLs3r27Lk2H9mpMrPbvIdpGLvossjMm04xZrKzNgG0JMNpecFNSVULNIrv98awi5lUlVBfMo2S3PwFVxjysdGUjvy1JMK9nVi8bc4dPfpMuQqhJEYvpHyQlu7vRRILaT1XMMllkDxuT0dRzhrjJ7YKYJj23Npue4Hl2GSbfrLKHmciqhLmRL0zLa5TJ1xXzKNmDa7neO5Pnuk8id5hZOuq9qNI6RxYXk3pZmKeZPyHUNKjVuSZoNjW7QrLG4DBDYnzb\/GL\/nor2ZTsOBAnffJBtTc9FdsbE+TGaalkZqsAard7+KEQ5vDCuqjglGRLPKvvbb79dacyNY1dwPTTCKSkp5lGy815UfZ1TN7Z\/\/\/5NnBl9wJw5CpO9EZinpQUh1LTIRVwXZTalu115w3YDt39vhxOUXakjq+z2mIo+W7JxMlNhCiO8\/PLLypy5uGBbNuZ5yl51zJmbsA2WD5hHyV5lzG3ZzHthwu+qjnk92cvEPAslCKFmQU+uFQQEAUFAEBAE5iMghCpLQRAQBAQBQUAQyAEBIdQcQJRbCAKCgCAgCAgCQqiyBgQBQUAQEAQEgRwQEELNAUS5hSAgCAgCgoAgIIQqa0AQEAQEAUFAEMgBASHUHECUWwgCgoAgIAgIAkKosgYEAUFAEBAEBIEcEBBCzQFEuYUgIAgIAoKAICCEKmtAEBAEBAFBQBDIAQEh1BxAlFsIAoKAICAICAJCqLIGBAFBQBAQBASBHBAQQs0BRLmFICAICAKCgCAghCprQBAQBAQBQUAQyAEBIdQcQJRbCAKCgCAgCAgCQqiyBgQBQUAQEAQEgRwQEELNAUS5hSAgCAgCgoAgIIQqa0AQEAQEAUFAEMgBASHUHECUWwgCgoAgIAgIAkKosgYEAUFAEBAEBIEcEPj\/IDtby0NK\/FwAAAAASUVORK5CYII=","height":282,"width":468}}
%---
