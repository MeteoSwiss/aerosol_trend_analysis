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
BRW_rd_short.BsB1_S11=BRW_rd.Bs11_S11; 
BRW_rd_short.BsG1_S11=BRW_rd.Bs21_S11; 
BRW_rd_short.BsR1_S11=BRW_rd.Bs31_S11; 
BRW_rd_short.BsQ1_S11=BRW_rd.Bs41_S11; 
%%

%%
% Backscatter begin with TSI in October 1997
% break point detection: 
% use at it the S0 data for TSP

BRW_rd_short.BbsB0_S11=BRW_rd.Bbs10_S11;
BRW_rd_short.BbsG0_S11=BRW_rd.Bbs20_S11;
BRW_rd_short.BbsR0_S11=BRW_rd.Bbs30_S11;
BRW_rd_short.BbsQ0_S11=BRW_rd.Bbs40_S11;
%use the 1 data for PM1
BRW_rd_short.BbsB1_S11=BRW_rd.Bbs11_S11;
BRW_rd_short.BbsG1_S11=BRW_rd.Bbs21_S11;
BRW_rd_short.BbsR1_S11=BRW_rd.Bbs31_S11;
BRW_rd_short.BbsQ1_S11=BRW_rd.Bbs41_S11;
%%
names_short=fieldnames(BRW_rd_short);

names_delete=names_short(contains (names_short,{'Bs1_S11','2_S11','3_S11','4_S11','10_S11','20_S11','30_S11','40_S11','11_S11','21_S11','31_S11','41_S11' }));
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
c=contains(names_short,{'dry','41_S11','Q'}); 
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

% PM1 absorption has not to be homogeneised since AE16 sampled only PM10/TSP
%%

%%

figure; %[output:1343ce23]
plot(BRW_rd_short.Time,BRW_rd_short.BaG0_homo1,'.'); %[output:1343ce23]
hold on; %[output:1343ce23]
plot(BRW_rd_short.Time,BRW_rd_short.BaG0_homo2,'.'); %[output:1343ce23]
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
c= startsWith(names,["T1";"T0";"T_";"P_";"P1_";"P0_";"N";"Uu";"U";"Bax"]) | endsWith(names,["dry"]) | contains(names,'Q');
N=names(c);
for i=1:length(N)
    BRW_tr.(N{i})=[];
end
%BsB_S: problem from 5.4.1993 to october 1997: to invalidate
% potential problem for the same period in exp scat GR --> keep and test
% trend results ? Results fine, keep these data 
% % Px=timerange('1993-04-05','1997-10-01');
% % BRW_tr.BsB0_S1S10(Px)=NaN;
%one time series for TSP and PM10

%backscat and abs: 1998-2026 
P2=timerange('1997-10-01','1998-01-01');
BRW_tr.BbsG0_S11(P2)=NaN;
BRW_tr.BbsG1_S11(P2)=NaN;

%Abs:
BRW_tr.BaG0_A11(P2)=NaN;
BRW_tr.BaG1_A11(P2)=NaN;
BRW_tr.BaB0_A11(P2)=NaN;
BRW_tr.BaB1_A11(P2)=NaN;
BRW_tr.BaR0_A11(P2)=NaN;
BRW_tr.BaR1_A11(P2)=NaN;
BRW_tr.BaG0_homo1(P2)=NaN;
%abs: invalidate 25.12.2009 to jan 2011 (end story leak)
P3=timerange('2009-12-25','2011-01-01');
BRW_tr.BaG0_A11(P3)=NaN;
BRW_tr.BaG1_A11(P3)=NaN;
BRW_tr.BaB0_A11(P3)=NaN;
BRW_tr.BaB1_A11(P3)=NaN;
BRW_tr.BaR0_A11(P3)=NaN;
BRW_tr.BaR1_A11(P3)=NaN;
BRW_tr.BaG0_homo1(P3)=NaN;
%%
%compute necessary exp, backscat fraction and SSA: use G and expS bg and
%expA bg and fit

%exp + BbsF
names_short=fieldnames(BRW_rd_short);
names_neph=names_short(startsWith(names_short,{'Bs','Bbs'}));
BRW_expS_cal=compute_exp_D(BRW_tr,names_neph,lambdaSC);

%compute abs exp
names_3w=names_short(startsWith(names_short,{'Ba'}) & contains(names_short,{'_A11'}));
BRW_expA_cal=compute_exp_D(BRW_tr,names_3w,lambdaAE);

% compute expA from AE31+33
Ca7=startsWith(names_short,{'Bacx'}) & contains(names_short,{'_A8'});
names_7w3133=names_short(Ca7);
%BRW_expAE_cal2=compute_exp_D(BRW_tr,names_7w3133,lambdaAE7);
BRW_exp=synchronize(BRW_expS_cal,BRW_expA_cal);
%BRW_exp=synchronize(BRW_exp,BRW_expAE_cal2);
BRW_tr=synchronize(BRW_tr,BRW_exp);
clear BRW_expAE_cal2 BRW_expS_cal BRW_expA_cal BRW_exp;
%compute SSA

BRW_tr.SSAG0_AE=BRW_tr.BsG0_S2S20./(BRW_tr.BsG0_S2S20+BRW_tr.Bacx3_A81);
BRW_tr.SSAG0_aePsapClap=BRW_tr.BsG0_S2S20./(BRW_tr.BsG0_S2S20+BRW_tr.BaG0_homo1);
BRW_tr.SSAG1_PsapClap=BRW_tr.BsG1_S11./(BRW_tr.BsG1_S11+BRW_tr.BaG1_A11);
%%

% trends in abs are very similar between AE31-AE33 and PSAP-CLAP --> use
% only the AE8+AE16+PSAP+CLAP results 
BRW_tr.SSAG0_AE=[];
BRW_tr.Bacx3_A81=[];
%BRW_tr.expA_fit=[]; %
BRW_tr.BbsB0_S11=[];
BRW_tr.BbsR0_S11=[];
 BRW_tr.BbsB1_S11=[];
 BRW_tr.BbsR1_S11=[];
 BWR_tr.BbsFb1_S11=[];
  BWR_tr.BbsFb0_S11=[];
   BWR_tr.BbsFr1_S11=[];
  BWR_tr.BbsFr0_S11=[];
%%
names_tr=fieldnames(BRW_tr);
% compute neph only on G
%BRW_tr.BsB0_S1S10=[];
%BRW_tr.BsR0_S3S30=[];
%BRW_tr.BbsB0_S1S10=[];
%BRW_tr.BbsR0_S11=[];
%BRW_tr.BbsB0_S11=[];
% compute neph only on G
%BRW_tr.BsB1_S1S10=[];
%BRW_tr.BsR1_S3S30=[];
%BRW_tr.BbsB1_S1S10=[];
%BRW_tr.BbsR1_S11=[];
%BRW_tr.BbsB1_S11=[];


%compute abs on homo and Bacx3 and SSA trends only on G, the longest time series
BRW_tr.BaB0_A11=[];
%%
BRW_tr.BaG0_A11=[];
BRW_tr.BaR0_A11=[];

BRW_tr.BaB1_A11=[];
BRW_tr.BaR1_A11=[];

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
BRW_tr.expA_bg1(Pexp)=NaN;
BRW_tr.expA_br0=[];
BRW_tr.expA_br1=[];
BRW_tr.expA_gr0=[];
BRW_tr.expA_gr1=[];
% BRW_tr.expA_bgAE=[]; %[output:7b24b1ba]
% BRW_tr.expA_brAE=[];
% BRW_tr.expA_grAE=[];

% exp S do the trend analysis with and without Aug 1993-Nov 1997 and compare the results for the full period and if there is a clear change in the 10y trends.
%%
% remove BBsF PM1 and other wavelenths
BRW_tr.BbsFb0=[];
BRW_tr.BbsFr0=[];
BRW_tr.BbsFb1=[];
%BRW_tr.BbsFg1=[];
BRW_tr.BbsFr1=[];

%%
% 27.8. too  many variables I don't know why
BRW_tr.BsB0_S1S10=[];
BRW_tr.BsG0_S3S30=[]; %[output:96dde337]
BRW_tr.BsB1_S11=[];
BRW_tr.BsR1_S11=[];
%%
[BRW_result_MK,BRW_result_LMSlog,BRW_result_LMSlin]=all_trend_STN(BRW_tr,BRW_st); %[output:1a6d931d] %[output:90bf4018] %[output:74accf82] %[output:8f6a210b] %[output:02c66054] %[output:275fae0b] %[output:68a84860] %[output:0e5e2927] %[output:6b329859] %[output:21853660] %[output:0ffe934c] %[output:8a0916d6] %[output:9d73ce37] %[output:633ef791] %[output:3c5fe4a0] %[output:47f9d4f4] %[output:3eb703dd] %[output:5dcd7821] %[output:27d03b35] %[output:1c17ca5d] %[output:767ec57e] %[output:806deac2] %[output:5d1ca422] %[output:0b31e208] %[output:60d29154] %[output:0087bb11] %[output:4ef60291] %[output:83b3a347] %[output:10f020f8] %[output:2aaaaf09] %[output:6b601709] %[output:4002a337] %[output:0da292be] %[output:34c77553] %[output:3f5cfa6a] %[output:2250e661] %[output:71f1d760] %[output:82b16049] %[output:32ca3c35] %[output:90a4cce6] %[output:862eeea8] %[output:78254fee] %[output:5af536c2] %[output:93e2212e] %[output:8d19b7b3] %[output:2856bd3b] %[output:564b1669] %[output:2f3e2c66] %[output:29e927b8] %[output:3813a8dc] %[output:2e2d94fb] %[output:4fb1dd73] %[output:590060cf] %[output:413172ea] %[output:29a3859e] %[output:484a7bec] %[output:6ed9ef45] %[output:1a6e0101] %[output:846b7d01] %[output:000e549e] %[output:1da35d86] %[output:82f20aa6] %[output:3d3b876c] %[output:31a277a7] %[output:8fe07ccf] %[output:40b2d253] %[output:17d1c91c] %[output:6f86febd] %[output:8e0ed078] %[output:7b911553] %[output:298ec69e] %[output:4ee779b1] %[output:248e1962] %[output:83bf22a4] %[output:84060908] %[output:08155abe] %[output:73710e2c] %[output:018dd32d] %[output:44d899d4] %[output:67388040] %[output:22a309cc] %[output:00e838d8] %[output:85d54d6a] %[output:993ae445] %[output:9b31a73d] %[output:36662d89] %[output:01d4841b] %[output:15201ad6] %[output:95d59e56] %[output:5cb5ef92] %[output:4a70d0e1] %[output:86b42ed4] %[output:1fd14807] %[output:34595fc6] %[output:938ffa44] %[output:44402cdf] %[output:01d04c02] %[output:7134ed46] %[output:9039e061] %[output:5ed933b6] %[output:565c79b9] %[output:7e2689a8] %[output:06e6c066] %[output:4171e70b] %[output:494ede43] %[output:7b004c68] %[output:874d8ec0] %[output:2dda597b] %[output:3ada91c6] %[output:128525a5] %[output:024dd3c2] %[output:9becd61c] %[output:348171fb] %[output:340be5b7] %[output:4c73e2a6] %[output:5b500bbc] %[output:93b6cf28] %[output:6c7e4d1e] %[output:937e65fa] %[output:9f12749d] %[output:14fedef7] %[output:7f8f71b8] %[output:43df0a2b] %[output:32d5edea] %[output:1568a178] %[output:393188c0] %[output:0e6bd3c2] %[output:3c75499d] %[output:92f7e111] %[output:32dc71af] %[output:87da0cfe] %[output:46e9006d] %[output:0f549bfe] %[output:5d19fc14] %[output:44484d29] %[output:13089376] %[output:1008e3de] %[output:0dfad22a] %[output:9e8c34ca] %[output:30f45b9b] %[output:51d0b936] %[output:89715877] %[output:0dc7802e] %[output:10733847] %[output:2b710032] %[output:5d5d6928] %[output:3753ca58] %[output:18baf1f7] %[output:8c3bab3d] %[output:65ff9912] %[output:6e9923ff] %[output:48540a2b] %[output:8db3c968] %[output:38ea5b7c] %[output:32e97787] %[output:2ff26752] %[output:42dc50e6] %[output:41376bf1] %[output:907db972] %[output:7981fd14] %[output:0b36dba7] %[output:34d214c2] %[output:5c93efdd] %[output:36b55a29] %[output:5d5bbf8d] %[output:4b16883e] %[output:7cbe1379] %[output:4d89d53b] %[output:9397a265] %[output:1b83538f] %[output:9bd1108c] %[output:93686c5b] %[output:860dd88f] %[output:1c489307] %[output:367decd9] %[output:546b111d] %[output:51142b55] %[output:8213be0b] %[output:90dc90d3] %[output:65f88126] %[output:5c0d4dc4] %[output:66e55553] %[output:3ce0b691] %[output:12129d59] %[output:95063375] %[output:2b869b28] %[output:031c6712] %[output:7cc6f0d7] %[output:0cc5f808] %[output:04817980] %[output:3080e591] %[output:95524bbf] %[output:54710bbd] %[output:9d06bcf8] %[output:077691bc] %[output:29b4bdd7] %[output:32035ee5] %[output:01853600] %[output:658072e6] %[output:2801f0a5] %[output:372b2fea] %[output:5554bccb] %[output:506bae5c] %[output:9126c750] %[output:4b24bfa1] %[output:89a59519] %[output:3494b737] %[output:2d1feeca] %[output:2ff4a3ca] %[output:6762415e] %[output:3e35ec14] %[output:96cdf2f1] %[output:5adb1609] %[output:89bbea40] %[output:8e63cdb7] %[output:585421c9] %[output:8a948b09] %[output:200ea141] %[output:8c92ee02] %[output:982c46b4] %[output:8f8834a4] %[output:704156ac] %[output:6a4c0247] %[output:47f151da] %[output:9d21ca6a] %[output:7fc53de5] %[output:772f5af7] %[output:7b54119c] %[output:51559519] %[output:69e42c04] %[output:3186ce1b] %[output:47547439] %[output:8e6bc093] %[output:9d134474] %[output:47231f80] %[output:2c0a50bc] %[output:078c1650] %[output:9d65e458] %[output:2c761628] %[output:62f1bbf8] %[output:59148269] %[output:9dce30cb] %[output:78ed7ead] %[output:4c07a296] %[output:06dc8a2f] %[output:445548e3] %[output:19eeab58] %[output:40debe50] %[output:544a3c5c] %[output:75ffc104] %[output:37c3ae91] %[output:67d65a07] %[output:196d887a] %[output:9809bce7] %[output:79fc8e81] %[output:34c31025] %[output:36558d63] %[output:4c536d58] %[output:2f1bc288] %[output:5f4a0cfc] %[output:3164b319] %[output:37e44cc3] %[output:062bea41] %[output:8a2e83b3] %[output:596aaa7e] %[output:05ede065] %[output:8d42046a] %[output:4ec52aa6] %[output:8d511a73] %[output:997613bb] %[output:6688afbb] %[output:4135bf82] %[output:26c9cd1d] %[output:2895fcf6] %[output:62d2836f] %[output:0ee813cc] %[output:142b0a16] %[output:67f38a16] %[output:34ac7265] %[output:4e51d712] %[output:64fe1588] %[output:3c2330cf] %[output:2a4c696d] %[output:6b30af8a] %[output:4fbb3e67] %[output:04e78057] %[output:4b6c3002] %[output:3024b889] %[output:43a1badd] %[output:481b0c93] %[output:3f3232c5] %[output:4cb41de4] %[output:213e149e] %[output:99a27963] %[output:59bca9f1] %[output:1edc1436] %[output:92e45a27] %[output:89a9978e] %[output:43a98f96] %[output:711bbcda] %[output:79e15497] %[output:19adfb2f] %[output:5637d336] %[output:7e671f54] %[output:12cdf73d] %[output:3bbddfad] %[output:76b100fe] %[output:390ccf2a] %[output:4469a233] %[output:30d4e6d7] %[output:713edae2] %[output:6dd81bf8] %[output:141220f3] %[output:06416cbf] %[output:3b25ddc8] %[output:313fa6c8] %[output:912241b9] %[output:14d949d6] %[output:64b1e93d] %[output:931bd380] %[output:4ad24684] %[output:233226df] %[output:9df1d2c2] %[output:9ab07325] %[output:216a19f9] %[output:7f6044e1] %[output:4a1f7854] %[output:8e9ee41e] %[output:7992a999] %[output:7b8ce3f5] %[output:8d883bbc] %[output:84193af1] %[output:297ca220] %[output:17703273] %[output:6d2fb422] %[output:831976a5] %[output:113f7888] %[output:233d1d1b] %[output:59c6040d] %[output:5b3e15e8] %[output:7400dfd6] %[output:801873b6] %[output:5c3fb85d] %[output:0ee00466] %[output:661c10e9] %[output:278dc67f] %[output:612d2e14] %[output:8a31ff0f] %[output:137e2f57] %[output:96b7b375] %[output:175c0342] %[output:7bf2a46a] %[output:325839f9] %[output:03ac8eff] %[output:9c6667b5] %[output:795d3d1f] %[output:22b85d3a] %[output:918473fb] %[output:11001239] %[output:984e4182] %[output:639ee66e] %[output:8b67ce71] %[output:9e78f5c5] %[output:56059f9f] %[output:87e9be89] %[output:6cb8b6ee] %[output:1d6c85c3] %[output:0dfe9da0] %[output:742b8293] %[output:7141f281] %[output:4afc3452] %[output:2623bc39] %[output:7fc0485c] %[output:4fb6c006] %[output:8561cd64] %[output:20d39efb] %[output:11431d9f] %[output:048e97ef] %[output:1852498f] %[output:9143a6f4] %[output:3380886d] %[output:47ffb32b] %[output:8588c821] %[output:135e38ff] %[output:454c7f7d] %[output:879c56fe] %[output:111bca2d] %[output:751c5705] %[output:5208834f] %[output:76b6c351] %[output:8bb63f0f] %[output:94b3a58d] %[output:8ce8b5ab] %[output:4bf4381a] %[output:79d6160f] %[output:969ef301] %[output:3773d76c] %[output:6dd507d9] %[output:9b0c1f72] %[output:2efecd0c] %[output:7599f9c9] %[output:0b9cde16] %[output:1cf8bdfa] %[output:8a595b1b] %[output:7c96b513] %[output:42879679] %[output:18d08194] %[output:0ec79177] %[output:27737cda] %[output:3b8b3559] %[output:7d2f958a] %[output:93559b24] %[output:470c9215] %[output:59c56c76] %[output:5237cd0e] %[output:508a7a92] %[output:4255bcda] %[output:7ec81959] %[output:0ade710f] %[output:5443e0d6] %[output:6c9e2418] %[output:875fd8b2] %[output:8b235608] %[output:3c8a13af] %[output:3e48f01d] %[output:3051d908] %[output:6b4dbad1] %[output:5644bfec] %[output:6c66ba1a] %[output:153e65a5] %[output:6e50782e] %[output:95bd2060] %[output:357ff63e] %[output:00b81aeb] %[output:12b8c40f] %[output:01d1f83d] %[output:33d2d9e1] %[output:6db989f8] %[output:25e66804] %[output:457f4c1a] %[output:7d574454] %[output:0db7d784] %[output:454b51c7] %[output:5c73b87f] %[output:39652ade] %[output:74357a07] %[output:5b2412d9] %[output:0903c4a8] %[output:234929f7] %[output:070a795f] %[output:99479e13] %[output:103d094a] %[output:0404a91a] %[output:4fa9468e] %[output:0bde47ba] %[output:85260890] %[output:10db9db0] %[output:4cdc6a3c] %[output:81963189] %[output:4b4c6db2] %[output:64f5ee8e] %[output:13139e4b] %[output:16cbb4eb] %[output:21ba8c85] %[output:3d0fb04c] %[output:8f126aca] %[output:8a2da152] %[output:311a487c] %[output:2170f141] %[output:6bddefc7] %[output:0da5f4ac] %[output:0817b4ed] %[output:0ee7b693] %[output:7a52dac7] %[output:7c16d95e] %[output:5195e48f] %[output:15c2e6a7] %[output:5d9bd20f] %[output:198efaba] %[output:80284d00] %[output:72e2ba8b] %[output:1835e793] %[output:4531b12e] %[output:1c381965] %[output:224b96be] %[output:635c9b22] %[output:445df048] %[output:6c09b771] %[output:80cf8149] %[output:3daa6ac9] %[output:90bb9036] %[output:758fd6f2] %[output:753ef725] %[output:7094ccf5] %[output:1ac20262] %[output:2169675d] %[output:187dec53] %[output:38b46c11] %[output:5c1cb559] %[output:7350707f] %[output:66461452] %[output:7cc48811] %[output:5b14209d] %[output:5b1610ef] %[output:112c2587] %[output:4de3aadc] %[output:008a21ee] %[output:43ad4fa1] %[output:4edaac14] %[output:8b903865] %[output:4d7cd259] %[output:15800c7c] %[output:6e9948b2] %[output:6f14ceca] %[output:883a501b] %[output:3c12e34d] %[output:7cb08fa9] %[output:469bc765] %[output:34071e90] %[output:0fadc2a5] %[output:97f618c7] %[output:1e270ebf] %[output:15b7c972] %[output:1d29dfc8] %[output:6a45d3fe] %[output:7c48b797] %[output:7a26972b] %[output:98a0eeee] %[output:1e52ecfd] %[output:34858ae1] %[output:6a5ad6af] %[output:1ccf8032] %[output:29198aa4] %[output:453264d7] %[output:6652ff6b] %[output:663158e8] %[output:833d3a1e] %[output:52afc33b] %[output:9965b3a1] %[output:65f3b8c0] %[output:4b734d75] %[output:8aa73ebc] %[output:99527e6a] %[output:8c00a7e1] %[output:82bf6ca1] %[output:83e751ea] %[output:959448c9] %[output:6046acdb] %[output:6d234253] %[output:2ab8be3f] %[output:9cfffedf] %[output:91fbd516] %[output:1ede8235] %[output:697cf7a9] %[output:547b319e] %[output:18430d68] %[output:16b30a3e] %[output:941a5882] %[output:3380b4bb] %[output:788649f7] %[output:1f54f3f4] %[output:06b5f3b9] %[output:2b2ef57a] %[output:1d6f898f] %[output:213aae46] %[output:22932d90] %[output:8f82340c] %[output:2f29075a] %[output:6b48a3f5] %[output:9c4369be] %[output:11a0a64f] %[output:5157b02b] %[output:9d6b948c] %[output:25fe6413] %[output:8386dd71] %[output:413d058a] %[output:3550837e] %[output:75fa8c8d] %[output:80ebbb5f] %[output:1aa4d602] %[output:12ea44e3] %[output:9a332fd1] %[output:091957ad] %[output:04d88352] %[output:2c6fb989] %[output:9f5e9c1d] %[output:72da9bba] %[output:1a9ba36d] %[output:42f90fcb] %[output:1e8345eb] %[output:15513712] %[output:8e53c04b] %[output:4c92f5a4] %[output:30b4837a] %[output:597192b1] %[output:73cfa7ea] %[output:14a98fd7] %[output:19e6b3fd] %[output:70003f97] %[output:06c51ea5] %[output:6c8bd1a7] %[output:1a0bcace] %[output:90364ea7] %[output:0251823f] %[output:40bdd972] %[output:7a040431] %[output:6acadad0] %[output:63115f0c] %[output:0801b24e] %[output:6d34fbf6] %[output:472e405b] %[output:31966b1b] %[output:5acb435e] %[output:7dea8a02] %[output:62a11b41] %[output:8f4f19a8] %[output:0608e2c1] %[output:4a8413df] %[output:60e8ca5b] %[output:9dd9441c] %[output:9dbcc9ad] %[output:5dfab478] %[output:582e4a1e] %[output:92d880f0] %[output:7d6d374f] %[output:948cde65] %[output:81af6044] %[output:8634c476] %[output:2b693640] %[output:3bafc0c4] %[output:971d3b0c] %[output:9bbea782] %[output:0aa815a4] %[output:623202dc] %[output:9e71f5e3] %[output:7dafc6ec] %[output:456551ba] %[output:38a604db] %[output:960579c1] %[output:5f59913f] %[output:6cc418b9] %[output:0c4e9014] %[output:21cc769b] %[output:9f24076b] %[output:636f52e4] %[output:0d7ea755] %[output:5e89a2d3] %[output:2483172f] %[output:88df6649] %[output:6334f9b0] %[output:330ad32a] %[output:2dc7fc8e] %[output:21c0df7a] %[output:8bb5f7ea] %[output:3051a516] %[output:336c33ea] %[output:0a085ee4] %[output:3c466ec0] %[output:015cb391] %[output:4b57dd35] %[output:18ac3c72] %[output:0da4bf13] %[output:30acf795] %[output:8ac12db5] %[output:3e21333c] %[output:19130dd3] %[output:0ac06016] %[output:0025699a] %[output:4c969a4b] %[output:3205fcd6] %[output:721da8fd] %[output:8d366bc5] %[output:983aa819] %[output:3f0b16eb] %[output:88db5cc2] %[output:6a51fa74] %[output:75d9574e] %[output:7f438b82] %[output:0f331783] %[output:9dc5979e] %[output:235cad5d] %[output:58b620ec] %[output:7b41c513] %[output:2ddb25f2] %[output:5ab2fc80] %[output:2f71bef0] %[output:37cf23e7] %[output:295c7b54] %[output:61077f08] %[output:2e4725a3] %[output:50878f63] %[output:7abb99b8] %[output:4e23b6c7] %[output:4278e853] %[output:7c936736] %[output:26a8f37b] %[output:8afe35af] %[output:6f7c22ed] %[output:50dc6dab] %[output:4d6c376d] %[output:3d9c5eab] %[output:49c72875] %[output:84ee9072] %[output:218cf845] %[output:2d1a72fb] %[output:8ad5b377] %[output:48e441d5] %[output:365b2f9b] %[output:36d3706e] %[output:18b338aa] %[output:5cf97442] %[output:4863cf71] %[output:2c592f84] %[output:34dd8e79] %[output:586a6115] %[output:488fb4b4] %[output:87fce924] %[output:0b207c75] %[output:3675ed46] %[output:43170eda] %[output:3251e38d] %[output:66d74be0] %[output:68b65577] %[output:7384198f] %[output:5ce446de] %[output:7698e34c] %[output:4474b738] %[output:191152a0] %[output:122b0f3b] %[output:5b03e27b] %[output:3e4b2265] %[output:2e2ff7bd] %[output:3965d6bd] %[output:0ffb08e2] %[output:2c147e00] %[output:34fde0c9] %[output:21ae4975] %[output:10b9104d] %[output:6566c6a0] %[output:5f904590] %[output:79bf07a8] %[output:0186dc17] %[output:363c88f7] %[output:36b180fc] %[output:58d87de9] %[output:8387a974] %[output:52c9e2c7] %[output:1831ac8b] %[output:66b7fcea] %[output:61607279] %[output:682774aa] %[output:28ced20c] %[output:752a74f9] %[output:6b516c0a] %[output:4094f7b9] %[output:506ff6e7] %[output:9336c533] %[output:8a8a2ede] %[output:967d09b3] %[output:8d5180de] %[output:74344efb] %[output:21aaddf9] %[output:71e11d2c] %[output:2d10f7e3] %[output:54be5c90] %[output:660ef67e] %[output:3398924f] %[output:169784af] %[output:56970fba] %[output:06b0e977] %[output:8a77ba81] %[output:4722ddaf] %[output:4e6de565] %[output:6a74debd] %[output:842a63d5] %[output:1d183d89] %[output:40de2c29] %[output:475024df] %[output:3a5a4286] %[output:76076ca3] %[output:8912acf4] %[output:0aefa498] %[output:2532de9b] %[output:6fa36ea3] %[output:8ceeb71c] %[output:80b37eb0] %[output:9417dabe] %[output:5a8c7a5d] %[output:87b7195e] %[output:8b698fa3] %[output:0dfda2b9] %[output:6c5ac0d7] %[output:1c185aa0] %[output:053ef558] %[output:7f12dde5] %[output:2e288b62] %[output:2edbd576] %[output:714aef69] %[output:49946398] %[output:7b78cc28] %[output:803b4b5d] %[output:7a4a7201] %[output:9327ad9a] %[output:1dabecb7] %[output:941aa8b6] %[output:93e78f5d] %[output:1448401d] %[output:49aab42f] %[output:4664b6ae] %[output:387f045f] %[output:24371ede] %[output:6cfb6040] %[output:3537b91a] %[output:0dcc4ea3] %[output:06638abe] %[output:88c538af] %[output:2462298d] %[output:4ef77851] %[output:8d5ffc95] %[output:0d81dbe2] %[output:0b080ef5] %[output:51290129] %[output:58623409] %[output:47d2303a] %[output:7c9fb410] %[output:2e7ed239] %[output:248e63cc] %[output:8638c5a2] %[output:6b284375] %[output:64c0ff0e] %[output:594ca13c] %[output:624d246b] %[output:2b0f3cfa] %[output:89588526] %[output:2505a68f] %[output:6e72b5f4] %[output:2f950a37] %[output:49dc6a9d] %[output:52f25e5a] %[output:92e68211] %[output:4b448e58] %[output:9f7c6739] %[output:7443477a] %[output:8027f5df] %[output:0a2a7959] %[output:9c60dee2] %[output:91c3db57] %[output:4f71c328] %[output:20a37044] %[output:4d127a48] %[output:4a03cd28] %[output:5aae1920] %[output:77fb2d40] %[output:18e5ffad] %[output:39dfb265] %[output:5caeece6] %[output:2c8dc615] %[output:07fd2a8d] %[output:6d4e90cc] %[output:0a545acc] %[output:8f37d48b] %[output:9f6022a0] %[output:275a974c] %[output:1f36175b] %[output:5906e80b] %[output:887499bc] %[output:1c4aea40] %[output:2793d475] %[output:1edb9035] %[output:6c7ae58d] %[output:3eb66204] %[output:0e363b82] %[output:6e92393d] %[output:11ccbdfa] %[output:587529f3] %[output:570e88e7] %[output:3900e64a] %[output:3575c0f9] %[output:2ad760b2] %[output:38cf78ae] %[output:0225c240] %[output:6a0488ba] %[output:4698c34d] %[output:89f194b4] %[output:7b7682da] %[output:31754962] %[output:2fbe3fed] %[output:2e2db140] %[output:12f9ea2b] %[output:16bbe101] %[output:31342d67] %[output:8e099250] %[output:82f43d66] %[output:68ac73ea] %[output:34721f8a] %[output:5942e1fe] %[output:006fe226] %[output:22545626] %[output:48cd084e] %[output:585bfd8c] %[output:1f6aef5e] %[output:2e2d6192] %[output:7557f411] %[output:95b6eee4] %[output:334cb600] %[output:20d579e6] %[output:254a64c0] %[output:13e374fe] %[output:4dd2fc4d] %[output:0855a80a] %[output:89775020] %[output:4faf8252] %[output:792b81d4] %[output:7155133c] %[output:58181a86] %[output:3c5bc103] %[output:4e1cfb09] %[output:1924c3bb] %[output:57f5f2e6] %[output:95d5428b] %[output:8bf89b19] %[output:4ddc9a86] %[output:6175de95] %[output:39fe8ddb] %[output:27c2871d] %[output:1837d5e0] %[output:7f25e7ac] %[output:0c2d92df] %[output:554cc82c] %[output:127911b7] %[output:29ac48cd] %[output:574e224b] %[output:45a4b885] %[output:7de238d1] %[output:4c7678b3] %[output:63c95266] %[output:2e8866ad] %[output:259a3e07] %[output:99abc3bc] %[output:87eaa7ef] %[output:85a3b8df] %[output:292ea369] %[output:92156ff3] %[output:2706e187] %[output:9fc54de2] %[output:597c3bd7] %[output:5ce550ef] %[output:2b2bdb07] %[output:5fb95bb7] %[output:1c1f91d3] %[output:88ecc3d7] %[output:6ce147c2] %[output:5f2f7ac8] %[output:93fe0ded] %[output:3b9f216d] %[output:5f667f47] %[output:5e2fa66c] %[output:2c58805a] %[output:83b4a87b] %[output:17d0ba2a] %[output:3f1f6d78] %[output:337611aa] %[output:66b4883d] %[output:8721a920] %[output:1d1581b1] %[output:60705e91] %[output:1cf2a9b5] %[output:516360a2] %[output:938e7d01] %[output:59b8d9ba] %[output:80e8f246] %[output:1950cda8] %[output:9a480775] %[output:976ba4c2] %[output:571549bb] %[output:4d5ebce1] %[output:822f0fdf] %[output:7bfd20d4] %[output:82c8b830] %[output:6a6f7ea2] %[output:34e17c9e] %[output:4ae1df08] %[output:8612d969] %[output:507d9e89] %[output:27533850] %[output:42aff600] %[output:2808f36b] %[output:1c1fd27c] %[output:0c83caf3] %[output:5be3ca54] %[output:88fda62e] %[output:77fb56cf] %[output:25ae451a] %[output:142011d3] %[output:149100a6] %[output:37a85aad] %[output:84cfd9d8] %[output:09cee69c] %[output:8792c568] %[output:9933e756] %[output:419746f4] %[output:79f87bf1] %[output:31b9d295] %[output:73151e9b] %[output:3969725e] %[output:6334ea98] %[output:00fb1fa9] %[output:4b8fcaa5] %[output:8b3d7c14] %[output:558d621c] %[output:62c51fd5] %[output:149ad87d] %[output:183cdfc6] %[output:73e4db2b] %[output:01a55b95] %[output:035add23] %[output:693605fe] %[output:8cdc2a0a] %[output:0ca93415] %[output:9553fce4] %[output:614b130b] %[output:3a5413b9] %[output:16f7fe47] %[output:3820925b] %[output:85d90bd8] %[output:9b26cc5a] %[output:75b252d4] %[output:5a30f061] %[output:84ad7612] %[output:2ebcbeb4] %[output:5b236c5b] %[output:1903177f] %[output:8b258047] %[output:12909183] %[output:49782fa3] %[output:49dc514c] %[output:26974022] %[output:16e80540] %[output:9213149d] %[output:16882105] %[output:49702400] %[output:93232c1a] %[output:20ffa0a1] %[output:16f84d9f] %[output:2cbbba97] %[output:20624876] %[output:47528640] %[output:44322dd4] %[output:52a7cbae] %[output:63b0d367] %[output:49be8e2b] %[output:2af5650f] %[output:9451822d] %[output:6ae73ee2] %[output:4fcee8dd] %[output:2a580ca8] %[output:393c3231] %[output:643bc4b2] %[output:451b5474] %[output:1e391f2b] %[output:29a99fc7] %[output:1ba3fb7a] %[output:91c5b0cb] %[output:1444427a] %[output:47034e22] %[output:77cce9d3] %[output:367bb50a] %[output:2bcf67b7] %[output:3a60695f] %[output:2c697b5a] %[output:6f3bbe72] %[output:21f632a2] %[output:13c77413] %[output:2f377d74] %[output:4be6877c] %[output:2bc77e66] %[output:9325c816] %[output:4459f20c] %[output:19803d9b] %[output:3aea5ef7] %[output:9d17b0a0] %[output:32b06ea8] %[output:831bf892] %[output:2b31a37d] %[output:2cdb3459] %[output:136f79ed] %[output:5490e1f0] %[output:797550e0] %[output:64978c26] %[output:54c675cf] %[output:46662cc7] %[output:457a8a73] %[output:3384717c] %[output:7c60e759] %[output:4ac55cb8] %[output:35d599c5] %[output:43bf38ba] %[output:770bc4c8] %[output:1f6ca1e1] %[output:1b2cfd90] %[output:1d2abde9] %[output:93151d9e] %[output:297189c7] %[output:134a2def] %[output:306ba9c8] %[output:65d6e7fb] %[output:03818a92] %[output:43c3bd7e] %[output:20ac780d] %[output:64547591] %[output:0828ae10] %[output:5a559b7e] %[output:9956869a] %[output:8e6ed4f4] %[output:645093b6] %[output:8e2234fc] %[output:06b92d0b] %[output:9047c90e] %[output:1fe031c4] %[output:1d2cbf9d] %[output:10e22b02] %[output:8e8e61eb] %[output:47b6c75e] %[output:3b262e29] %[output:54bc881e] %[output:64ead7e5] %[output:3de75cba] %[output:339585c1] %[output:5d38f3af] %[output:6d6ffd01] %[output:80e4e314] %[output:32b644ae] %[output:1d97f9c7] %[output:2941dc63] %[output:9c937f8c] %[output:0d577d2a] %[output:81d5fd03] %[output:98c2f4af] %[output:1093bee7] %[output:5f2679e1] %[output:74ec2644] %[output:69d38380] %[output:5b941da0] %[output:4f788cb6] %[output:4f4e71b4] %[output:12f2f3c8] %[output:67d4649f] %[output:9e0cdda0] %[output:0b6d3039] %[output:665d457a] %[output:76148732] %[output:9fbaba72] %[output:9de1d054] %[output:99598631] %[output:99a888c1] %[output:1935039b] %[output:2764cea3] %[output:8d9e61d0] %[output:761d0d20] %[output:3f61426e] %[output:43474f7f] %[output:0677bebb] %[output:7246af85] %[output:137ff240] %[output:07303485] %[output:63c4e12f] %[output:3ceac80d] %[output:1623d026] %[output:5862ddd2] %[output:46976e3e] %[output:95cf5821] %[output:1ebeed9f] %[output:374a7a19] %[output:30c12796] %[output:2e41aadb] %[output:2972a2dc] %[output:98fa4d2b] %[output:72571873] %[output:50a01290] %[output:326b8715] %[output:72fea0ab] %[output:2ffeb845] %[output:43e83f4b] %[output:75ee40ca] %[output:5c99078b] %[output:9097e384] %[output:812dc5eb] %[output:1a62274f] %[output:43eed418] %[output:204680d3] %[output:9775e68f] %[output:0c3855df] %[output:25e33875] %[output:5e806063] %[output:9ea7c506] %[output:30d15316] %[output:495a2236] %[output:10d7acc6] %[output:6071fc59] %[output:5d614f7c] %[output:9e7cb58a] %[output:0c48a3cc] %[output:37ce7784] %[output:583e5487] %[output:343ab644] %[output:497623aa] %[output:3e48d93f] %[output:5cd7547e] %[output:5d0b0ff4] %[output:9e54eab1] %[output:1b723b65] %[output:2047c511] %[output:87b50185] %[output:6ec126cb] %[output:01edad82] %[output:01f13939] %[output:2e91b473] %[output:90426e8c] %[output:2db378ef] %[output:6eb51b03] %[output:0b7cac90] %[output:27b08062] %[output:53448ab6] %[output:2f5eae9e] %[output:6d4851d3] %[output:68260f97] %[output:00433012] %[output:03a2d693] %[output:48128d6a] %[output:3042d64c] %[output:8a5e9051] %[output:23c19e0c] %[output:420497b8] %[output:7452df56] %[output:7a736a82] %[output:76fe587f] %[output:7dfde792] %[output:59cbd61f] %[output:40f048b9] %[output:18aadaba] %[output:413b3f92] %[output:5d10677a] %[output:32c834c1] %[output:3e9b45ce] %[output:5ac1acb4] %[output:925368cb] %[output:5b7b5d88] %[output:7ad62ba6] %[output:2117c05d] %[output:79165849] %[output:05ab930e] %[output:226a1551] %[output:4f13862e] %[output:505359ba] %[output:30a55861] %[output:9df1eb23] %[output:05db6a0b] %[output:50a68c5b] %[output:58bf8f3c] %[output:494a2544] %[output:5c7fa512] %[output:88439206] %[output:402cd7fa] %[output:6d727988] %[output:8e52e816] %[output:6c3c05c2] %[output:225e283c] %[output:4d075d6c] %[output:6a852b34] %[output:9dd8819e] %[output:9416bc8e] %[output:9a21e657] %[output:17bf66eb] %[output:5c9f6557] %[output:7deb8f7a] %[output:0b3568db] %[output:0bd65661] %[output:85f53734] %[output:387e1835] %[output:50f99c17] %[output:5c2c9cce] %[output:79f2d8a5] %[output:1640effb] %[output:30bad384] %[output:622f2875] %[output:2a6db916] %[output:0bc7e39b] %[output:1d9aa334] %[output:6822ab16] %[output:1d78aa78] %[output:513b5870] %[output:16d85161] %[output:3b8a08dd] %[output:6d447256] %[output:01b7c5cf] %[output:529e6bf9] %[output:0a2bdb55] %[output:0e8aaaf5] %[output:4db64809] %[output:23883739] %[output:57c19e61] %[output:761647bd] %[output:13e0d75a] %[output:1307ef03] %[output:6960a139] %[output:642d79d3] %[output:4817e881] %[output:6b411a3d] %[output:32d574fb] %[output:5d5d6f3e] %[output:80b83034] %[output:6cc388f7] %[output:9b58a23c] %[output:351821d8] %[output:799097b0] %[output:799a2222] %[output:751187aa] %[output:93c84f38] %[output:4ab507e0] %[output:27fc40c2] %[output:0cd7f1a6] %[output:0ce166ee] %[output:7af87bf1] %[output:7a0a17ca] %[output:92ab2c20] %[output:6288e5e7] %[output:0a6e18f7] %[output:0a176729] %[output:6e628886] %[output:659c02bd] %[output:33ff1abe] %[output:15b24b94] %[output:4dbe7346] %[output:26d28852] %[output:14486678] %[output:1c8dbbdd] %[output:25446015] %[output:8a26ff36] %[output:1079d9d5] %[output:532b75d8] %[output:594b51c6] %[output:1e030c7b] %[output:24f94620] %[output:6f22a0e4] %[output:9f875c39] %[output:23bc46ed] %[output:89e084da] %[output:60680960] %[output:95bbfda0] %[output:61c2c8ee] %[output:7f0eb7d1] %[output:618c5168] %[output:14be3dfa] %[output:38627e4c] %[output:0d726a71] %[output:2432710f] %[output:29d789ca] %[output:1a81a99d] %[output:0e5ba1ac] %[output:3fe4140f] %[output:22bc5f8e] %[output:6a2be139] %[output:98b7854f] %[output:97daa11c] %[output:3de97246] %[output:90d3f02f] %[output:3169c5e4] %[output:8cab4c98] %[output:6f8cacdb] %[output:669c3baa] %[output:75abd782] %[output:79558ee9] %[output:878b3c2f] %[output:39376864] %[output:494dc585] %[output:7bf9cbe3] %[output:5c76acaf] %[output:5ee6f087] %[output:35dffbc2] %[output:2f7a0691] %[output:1771c5f5] %[output:9a52888a] %[output:588a5069] %[output:1094d469] %[output:0bddb1a1] %[output:069a3b96] %[output:96e85c40] %[output:42e7542b] %[output:0dfbdc8b] %[output:6a4f5f2b] %[output:489cd351] %[output:9e31f403] %[output:1341bf5d] %[output:804f7253] %[output:40035a89] %[output:0670583c] %[output:4d2b7df1] %[output:92114223] %[output:8ece0d9f] %[output:1fc1adf7] %[output:1285432b] %[output:33a705b2] %[output:5d8dd318] %[output:1a765f3e] %[output:7f30baa3] %[output:254a4d00] %[output:35056cc3] %[output:12083e87] %[output:323b064b] %[output:855404f3] %[output:4655327e] %[output:8b0d5cbf] %[output:45eb5feb] %[output:8e7b94db] %[output:56a167fc] %[output:48fc1d9b] %[output:002181ae] %[output:53c4d703] %[output:96c83a29] %[output:1cc98398] %[output:2f373d5a] %[output:8dfa3675] %[output:3f7de47a] %[output:8ab75b05] %[output:78d01d21] %[output:13ccb897] %[output:6390c0fd] %[output:28c1feba] %[output:079374d7] %[output:5da4659e] %[output:1d882d06] %[output:6df06c49] %[output:64f922dc] %[output:9207c7f1] %[output:7ffa277b] %[output:6a9daaf7] %[output:6339c644] %[output:78bbd259] %[output:02d2b3ca] %[output:5845d840] %[output:350b65bf] %[output:2684b348] %[output:9ecb2ca2] %[output:89f08391] %[output:30ad95e5] %[output:2cc87424] %[output:4f28d897] %[output:0cb25264] %[output:2482bc12] %[output:8b90400e] %[output:6742217e] %[output:02656923] %[output:6e4ee37f] %[output:09f121b4] %[output:922375a8] %[output:5f085ef2] %[output:453aef2c] %[output:7b8de5bf] %[output:5afdf3ba] %[output:8a3f2782] %[output:1da5c1c7] %[output:3f6183e4] %[output:8283c67a] %[output:29ae0c19] %[output:8ca58c23] %[output:9fbe16fe] %[output:5f25d80d] %[output:8cf93ee7] %[output:6e8f14b8] %[output:76a0f9bf] %[output:9fd6209f] %[output:61d6acd8] %[output:3a8cce0c] %[output:197dc401] %[output:44bf8509] %[output:89e0773f] %[output:348f40ee] %[output:64fa1a60] %[output:81b04c90] %[output:51da59dc] %[output:25b6a4ef] %[output:3732c93f] %[output:142feeb8] %[output:60292b2d] %[output:646cf396] %[output:72a495af] %[output:5de2c74f] %[output:1781cb6c] %[output:01691c52] %[output:1be70811] %[output:847b0d16] %[output:2f97abd9] %[output:47dae02b] %[output:8677fd68] %[output:9cc0d729] %[output:5d61aaf4] %[output:58cee9be] %[output:01fcf444] %[output:194e2fc7] %[output:52a95ed6] %[output:4c5328bb] %[output:993c1df4] %[output:790ddd1c] %[output:706628e6] %[output:55853ada] %[output:5cd73f0e] %[output:1a1662c4] %[output:471c26ce] %[output:27f30a75] %[output:13eabec7] %[output:7c8f710c] %[output:99a6eb62] %[output:7e32c918] %[output:15c9d771] %[output:7752e6ba] %[output:9721c460] %[output:7ec915a9] %[output:9d631b47] %[output:581aec3f] %[output:6974ba62] %[output:71bd74e2] %[output:7c9b23aa] %[output:2ecec322] %[output:4d9f796f] %[output:571b729e] %[output:97e718cf] %[output:72485a3d] %[output:9a327883] %[output:419dedf1] %[output:88fbbca5] %[output:523e9971] %[output:86795c12] %[output:389fc8e6] %[output:0eff3809] %[output:35bcab0b] %[output:839c4d90] %[output:485f12a3] %[output:74a60b48] %[output:472bca99] %[output:8d8b7e28] %[output:23757618] %[output:4a24609d] %[output:6eba130f] %[output:008b2ee0] %[output:93b99ca4] %[output:4b81ff2e] %[output:7722cfc0] %[output:4a7c87ff] %[output:0df96489] %[output:166f89de] %[output:0ec7b580] %[output:82e38465] %[output:05c35e82] %[output:3e24637d] %[output:649e55f1] %[output:5985f94c] %[output:03d9874f] %[output:00b802f5] %[output:515739af] %[output:3e37efd3] %[output:06e48b48] %[output:1307e453] %[output:7a52dc21] %[output:78c32ff0] %[output:21bbdd1c] %[output:4709fb01] %[output:4c6945db] %[output:0dc8ca2c] %[output:7b0608da] %[output:4c5ec070] %[output:851ed762] %[output:59b338a9] %[output:6c484356] %[output:0c56a874] %[output:6aba0c4e] %[output:4381e90c] %[output:744b45c2] %[output:65d4d8a6] %[output:40e6b0aa] %[output:5ee4043f] %[output:1227cc3c] %[output:75aad05f] %[output:4df85a43] %[output:444fdf86] %[output:1d906467] %[output:8a15888e] %[output:9adaf37a] %[output:6b8efe90] %[output:821b0184] %[output:6471534b] %[output:95fac716] %[output:522e6250] %[output:1e36e552] %[output:13dc14f4] %[output:51ed18bf] %[output:45d7858c] %[output:6410adcd] %[output:6c817e89] %[output:1d3f73c0] %[output:071ee91f] %[output:93760aaa] %[output:07e0d602] %[output:5947e9a7] %[output:34d35ca0] %[output:1bc72563] %[output:31357896] %[output:60fbdd18] %[output:874a956f] %[output:2ac32125] %[output:7c425f8d] %[output:7828c641] %[output:7c925f22] %[output:62be5cc2] %[output:579a0b45] %[output:5675d56b] %[output:297119c8] %[output:06f94e52] %[output:239bcabe] %[output:54f79a2e] %[output:7b603af7] %[output:42791e4b] %[output:992c3f08] %[output:64c653d0] %[output:8b530051] %[output:20bcca1a] %[output:11bcec20] %[output:05d75d5e] %[output:6ced4381] %[output:6c3a5ff1] %[output:5e56df8b] %[output:43ee64f2] %[output:97434ce9] %[output:49e68588] %[output:59786baa] %[output:5b1e622c] %[output:07042ed2] %[output:8658390a] %[output:7a7a99f8] %[output:0a37aa50] %[output:975f28a7] %[output:2d2d6a85] %[output:278b81a0] %[output:6dc65d00] %[output:0b493abe] %[output:23f79076] %[output:89184c2e] %[output:407c20c9] %[output:177860ac] %[output:3cf3f165] %[output:78fa46bd] %[output:93fb6416] %[output:94863d74] %[output:11ba9b4c] %[output:4f97591e] %[output:66ba48f6] %[output:5d453ff0] %[output:5910b602] %[output:71fcc425] %[output:92906f6b] %[output:40c6b934] %[output:27de264b] %[output:145b9841] %[output:40b9fc49] %[output:976ecd41] %[output:009b9d83] %[output:62f64890] %[output:6b538123] %[output:8efabc44] %[output:4fb948ee] %[output:53c2b438] %[output:92a51c96] %[output:6f8aae92] %[output:26de4ece] %[output:9814d4b4] %[output:670bb8ab] %[output:244b0e9b] %[output:5a665ac9] %[output:86aa6b35] %[output:91ed1a99] %[output:395d78c8] %[output:799cd3a3] %[output:5630bb82] %[output:13743b4b] %[output:7a4199a0] %[output:34119526] %[output:4479e3b8] %[output:0beac469] %[output:4a10c6bd] %[output:5db1a1f1] %[output:614e2e40] %[output:8e06d8e3] %[output:22c53bcf] %[output:1c96f001] %[output:62ace964] %[output:8e1c68bf] %[output:8a8c186b] %[output:14e74ae1] %[output:69854e3e] %[output:55704137] %[output:63dca84c] %[output:39180e8b] %[output:7856fa49] %[output:81b8313d] %[output:2586ba6b] %[output:15196897] %[output:641aa428] %[output:07828f93] %[output:2814c7a3] %[output:0f5652ad] %[output:9cba3311] %[output:1df52f65] %[output:48719809] %[output:700678fd] %[output:293250da] %[output:20789b6f] %[output:9fd94820] %[output:03d062d5] %[output:164d12c1] %[output:2131ddf9] %[output:2be64074] %[output:578fdf12] %[output:68e7d81a] %[output:93bea9a8] %[output:7774f180] %[output:495fb57a] %[output:7c57f034] %[output:0b24147d] %[output:24369605] %[output:7db1dba9] %[output:80e5033a] %[output:1043d283] %[output:51176c6a] %[output:48b94a28] %[output:0f044448] %[output:3b44c548] %[output:8a7c4199] %[output:31a4b364] %[output:0b0787b0] %[output:0a83c526] %[output:773d97aa] %[output:47c2d859] %[output:60da6422] %[output:0c9f3a0b] %[output:4cd27df9] %[output:379436fe] %[output:68d7579f] %[output:4a11aa6f] %[output:765817dd] %[output:88ba3736] %[output:9b539c0c] %[output:88cc8ac2] %[output:55155896] %[output:3cc8d5b2] %[output:1759c647] %[output:1ecce73a] %[output:57f9c687] %[output:2f4789b9] %[output:4f682b8a] %[output:4d5d1de9] %[output:0b416b7d] %[output:270240f2] %[output:723c5045] %[output:86f73dfc] %[output:33e84adc] %[output:60d2ced1] %[output:0bad9679] %[output:26d9778e] %[output:699a18c9] %[output:89d21e33] %[output:2dbb107d] %[output:0f729257] %[output:1316a092] %[output:9ffd965d] %[output:89a5fbd4] %[output:58acac1f] %[output:6ca9239f] %[output:7988adf9] %[output:45e4edce] %[output:4a1f824b] %[output:959819b0] %[output:92120768] %[output:49f70b33] %[output:67432884] %[output:895c42a2] %[output:812df8d7] %[output:16e488e3] %[output:41f2923d] %[output:7f334eeb] %[output:5eeff1a4] %[output:2ab71e87] %[output:07e28390] %[output:22038812] %[output:736acd05] %[output:70f830f5] %[output:51b12777] %[output:35bf91a2] %[output:88828cfe] %[output:1db79aaf] %[output:68583dfb] %[output:23f92763] %[output:0f6e0587] %[output:066e3c5a] %[output:99a44186] %[output:572f7343] %[output:34176f95] %[output:0d287aa0] %[output:0c7129e5] %[output:149d493c] %[output:0e4f977d] %[output:1bd17cfd] %[output:2856ad6d] %[output:592dc527] %[output:6f471478] %[output:95738865] %[output:0c5e54ca] %[output:4b92eefb] %[output:17d8bf04] %[output:2f99ea5e] %[output:6ff8ae3f] %[output:9ff4fcc2] %[output:5fe4665b] %[output:9b77d868] %[output:6ce3a7a6] %[output:5b8b01e5] %[output:200e7cb7] %[output:7b7369f4] %[output:0dc3c9d3] %[output:39dec4a8] %[output:7953c2e6] %[output:49693945] %[output:16422f80] %[output:70430658] %[output:3cae4c1e] %[output:9edd77cf] %[output:6127f0f1] %[output:8ef2719e] %[output:279542f7] %[output:50a8d3bb] %[output:0bcd6429] %[output:4beeb40e] %[output:77648f37] %[output:59e0d2c5] %[output:580b4dc1] %[output:5fa5857f] %[output:1c977cae] %[output:8eafcaf9] %[output:6ca64294] %[output:72c65301] %[output:1d7f94b5] %[output:04130bff] %[output:0e30a1ab] %[output:18e9e11f] %[output:27674809] %[output:9d7002d4] %[output:418f95ce] %[output:2c811d35] %[output:8ee28afc] %[output:3f6f4d48] %[output:96757d87] %[output:04b1952f] %[output:9226e66e] %[output:928c2ff6] %[output:8f277d4a] %[output:940f7643] %[output:63fde13e] %[output:1b09f981] %[output:738a7a29] %[output:8a1aa569] %[output:9ef87133] %[output:3279a919] %[output:0ae29569] %[output:03c972c5] %[output:9f1c3495] %[output:5f1d3331] %[output:3ea129a5] %[output:30cf2c9a] %[output:1d783524] %[output:663e1ee1] %[output:52a01cad] %[output:1b0da752] %[output:32c370c1] %[output:0f352396] %[output:0be73d84] %[output:2825f824] %[output:455165c8] %[output:0014d4bd] %[output:03504aae] %[output:7669f20d] %[output:3fa3440e] %[output:6f6d0f3e] %[output:51f8cd54] %[output:75ce3df8] %[output:786c44c0] %[output:91e365c6] %[output:80b86e44] %[output:2eb70ff2] %[output:014ecd58] %[output:3369c2b9] %[output:815b35d4] %[output:9452b033] %[output:150c3479] %[output:6b4347c6] %[output:35a8ac56] %[output:8b890b43] %[output:97331f43] %[output:7ce9c94d] %[output:481f71e1] %[output:3371e760] %[output:9b850a5a] %[output:2e5c4886] %[output:25ba1ac9] %[output:938ee6b8] %[output:064e3e42] %[output:664dfb89] %[output:6217d75d] %[output:51b66828] %[output:222b49be] %[output:86d9d598] %[output:55cbaaa3] %[output:2feac792] %[output:8340fa75] %[output:086fc12c] %[output:1fce6c76] %[output:3c7c2096] %[output:59de9faf] %[output:6ac47599] %[output:539f56f8] %[output:6252c30e] %[output:368df944] %[output:8c333f62] %[output:679b3d76] %[output:6a5837b6] %[output:345c63db] %[output:5405f067] %[output:4da1c54a] %[output:7378e65c] %[output:39a03f01] %[output:5265e380] %[output:76aee1eb] %[output:55de75f6] %[output:864337b4] %[output:5e399c31] %[output:0e2c4766] %[output:9b3dcc7c] %[output:47176382] %[output:929e2be5] %[output:3b5456fa] %[output:699f0ac0] %[output:1d6e9e86] %[output:8ea78f57] %[output:6da35736] %[output:88d0fa8f] %[output:8f4c0949] %[output:1e55e069] %[output:4eccc0d3] %[output:22c2f0ab] %[output:77ab0f69] %[output:61b578e3] %[output:4418a538] %[output:577b3ad2] %[output:06185cef] %[output:8f508b40] %[output:586897f1] %[output:77f5e924] %[output:39e8a051] %[output:6cd885ef] %[output:907e2d0f] %[output:206a6d7c] %[output:928d7028] %[output:3ecd7a8f] %[output:4d583c4c] %[output:930a3f1e] %[output:285110be] %[output:120ecb89] %[output:0bfe3b70] %[output:4f2147ad] %[output:6c45117f] %[output:82dbbd37] %[output:544a9231] %[output:8be28180] %[output:0bc8f7c1] %[output:1fee6824] %[output:4072a0fc] %[output:505367bc] %[output:9a9bc25f] %[output:542fe636] %[output:934ab92d] %[output:917aab98] %[output:65a6f6c6] %[output:16fe68b3] %[output:2fe359fe] %[output:7067095d] %[output:97db2197] %[output:8573531d] %[output:361314f1] %[output:55f218c8] %[output:486ba7ca] %[output:70b0a8ae] %[output:8296999d] %[output:453bfa31] %[output:7c656203] %[output:0f271821] %[output:0be699fb] %[output:4d23f8e6] %[output:201d73d4] %[output:9b67e889] %[output:9c87c3f2] %[output:3bf84064] %[output:7f92fba3] %[output:17540ef7] %[output:42504ecd] %[output:8e31910e] %[output:3e024ea5] %[output:0abbff9f] %[output:2919cae3] %[output:83f77888] %[output:478ec03e] %[output:0f37217c] %[output:022ae550] %[output:12032b4a] %[output:05ebef81] %[output:655bcabd] %[output:05820577] %[output:1c6e7405] %[output:415a0e41] %[output:715bed63] %[output:2bf5cf70] %[output:5b31e1e5] %[output:60c616ff] %[output:78163ffb] %[output:79209768] %[output:297a5b07] %[output:806cb1a1] %[output:7a42be28] %[output:94612da5] %[output:6c64e8f1] %[output:98f9dcbe] %[output:6eb97f16] %[output:5de2cc45] %[output:44cb88ca] %[output:5d8e6622] %[output:5aa2a68b] %[output:8cae9f00] %[output:4c9f7998] %[output:8d7ee1ea] %[output:221d3d10] %[output:343ac351] %[output:308a2c72] %[output:552567d1] %[output:689da17a] %[output:24cc3bc6] %[output:36943925] %[output:5c0ca563] %[output:976ef4f7] %[output:2f5060a5] %[output:3371b0fc] %[output:47aa55af] %[output:016220bf] %[output:7f616a67] %[output:88404826] %[output:2220331e] %[output:849a735f] %[output:506c223f] %[output:0afabfdc] %[output:7d9d64a7] %[output:13c0459b] %[output:5ae7143e] %[output:30734195] %[output:80562a9d] %[output:8b3dbc4d] %[output:61e2410b] %[output:41311589] %[output:42c1bc90] %[output:0de4d03f] %[output:6b39fca8] %[output:4c95a6e3] %[output:0b12b34d] %[output:5fc18053] %[output:8b5adb59] %[output:020647e0] %[output:66698d0d] %[output:80c380d3] %[output:2e23abdd] %[output:806fad47] %[output:623ed63b] %[output:03bf74cd] %[output:84dcdca4] %[output:8760e01c] %[output:1145aa2c] %[output:3f7c03ed] %[output:99b1b93f] %[output:2b68b300] %[output:7c12e1d7] %[output:14aae969] %[output:4b601021] %[output:084d9598] %[output:19a8ed18] %[output:65a49127] %[output:611152ed] %[output:4a961c66] %[output:64d5fe13] %[output:22ab7464] %[output:1b5c0dd0] %[output:99c54a1a] %[output:91bfb35c] %[output:79611708] %[output:006b7de0] %[output:61a15a78] %[output:01d57ae2] %[output:1205d11e] %[output:8d5c0c10] %[output:35206373] %[output:5d7defcd] %[output:4c211be5] %[output:3be6e255] %[output:74414eb2] %[output:1e431866] %[output:4dfcf25e] %[output:0458d482] %[output:359c584c] %[output:5667cf60] %[output:9687c67e] %[output:33528839] %[output:1f0c8b52] %[output:7da91a73] %[output:179bf83b] %[output:02895f06] %[output:9677e3c6] %[output:90a77b41] %[output:1ed3464c] %[output:66601631] %[output:176a0ac5] %[output:31ca3fd8] %[output:0b6c612e] %[output:286834a9] %[output:9189f8ad] %[output:71d55744] %[output:51463443] %[output:9efa67fc] %[output:83f712c1] %[output:7f8609ea] %[output:991f70d2] %[output:3637cba7] %[output:77d350c2] %[output:9d172f85] %[output:5daf6efe] %[output:90fa2ae6] %[output:16ebe48a] %[output:76295b78] %[output:6406f22b] %[output:3cf51fb2] %[output:5d740cd9] %[output:878dcc2f] %[output:8203641d] %[output:4e78cad9] %[output:570cc0a8] %[output:3c7d5b4c] %[output:9a184af6] %[output:3f8bde14] %[output:746c199f] %[output:316f8a74] %[output:7fed47cc] %[output:3cbb0e2d] %[output:9a51ba0c] %[output:24903316] %[output:69491921] %[output:0994489c] %[output:5a8c7f1a] %[output:46a0cf4f] %[output:0922354f] %[output:4e17a5f3] %[output:16f186a4] %[output:3d27272b] %[output:2cdad14d] %[output:733d9ed6] %[output:6a6bcdd4] %[output:1ff773a3] %[output:563f32a9] %[output:4acf3a67] %[output:9c708e0e] %[output:92d63e72] %[output:6070f6f7] %[output:711b6604] %[output:10ee5620] %[output:77769871] %[output:4be26cd8] %[output:73f0e7a9] %[output:598219a6] %[output:0b58b3a8] %[output:7d5cc8dd] %[output:6d794174] %[output:106e1883] %[output:84f5e974] %[output:5f6cb2b1] %[output:908f0d98] %[output:4f994b9e] %[output:36f30acc] %[output:0b43684e] %[output:1f9adec4] %[output:2b8d9e1e] %[output:8383f4cb] %[output:5705358f] %[output:04412afc] %[output:21884795] %[output:3ed89835] %[output:263274ad] %[output:3e5a20c2] %[output:1b836e91] %[output:7e8f1cdb] %[output:8cf78d7e] %[output:46ee256a] %[output:0d9a0261] %[output:1eafd698] %[output:739af566] %[output:47d0aff7] %[output:41beb7f4] %[output:0c1e3a4d] %[output:564cb46b] %[output:40b12c44] %[output:0ff329d8] %[output:502bfe75] %[output:604335c8] %[output:75988d4b] %[output:3a99c051] %[output:4395a1b5] %[output:0b32b2f2] %[output:885f9a17] %[output:03c23071] %[output:67d8966c] %[output:02cc1c29] %[output:2a006c80] %[output:0343c1d2] %[output:4feaec59] %[output:7eaac0e3] %[output:3a7a6647] %[output:7e176b3e] %[output:9d30311a] %[output:9938efa8] %[output:590e3702] %[output:6b346380] %[output:2d40d8d6] %[output:64b2a5c8] %[output:7255f7d9] %[output:251aca26] %[output:699921ae] %[output:13a5a259] %[output:6769b710] %[output:9e57fca6] %[output:9f2f05ac] %[output:9f15c02d] %[output:214ac371] %[output:0582a2e6] %[output:60bb2e9c] %[output:2f12491e] %[output:102dac56] %[output:474ead07] %[output:6377e7b4] %[output:79f94a8f] %[output:49238d14] %[output:3b8657ad] %[output:83efc422] %[output:3c02e4fe] %[output:82fc39d2] %[output:018719c8] %[output:49d53b42] %[output:61f32cee] %[output:676c34f2] %[output:9b92e358] %[output:4300a1b9] %[output:8852d368] %[output:58d48393] %[output:92888563] %[output:6d270f07] %[output:6ede0ae1] %[output:9657f593] %[output:0db46d5d] %[output:6ace1df8] %[output:59c165cf] %[output:65954e35] %[output:92f3f4a7] %[output:2bc8cedd] %[output:68005111] %[output:0fb71617] %[output:838bae79] %[output:29ecd218] %[output:13bc4764] %[output:89542deb] %[output:870df706] %[output:2e0c0dd5] %[output:7cd00b02] %[output:41a19d04] %[output:476d05c5] %[output:592bd887] %[output:9b1c18e2] %[output:1f4720cd] %[output:865ad4db] %[output:4a29de84] %[output:3f59608f] %[output:60fcd37e] %[output:54c9233b] %[output:402f46e8] %[output:2b3ed03c] %[output:75952393] %[output:600d18ed] %[output:386244e9] %[output:04f8b022] %[output:356d799b] %[output:7887a45b] %[output:2bce111a] %[output:3c76a8a0] %[output:7a1aa8d3] %[output:5554110c] %[output:41a6ff2c] %[output:2167e7a6] %[output:44b3ad23] %[output:08c52778] %[output:1c076e3a] %[output:80fa801e] %[output:3ece19b5] %[output:72cd91b1] %[output:3dd42aa4] %[output:48b75621] %[output:5668f291] %[output:99da434b] %[output:22dbe781] %[output:4920a59e] %[output:0815c34a] %[output:4b4c0a1a] %[output:93eaa66f] %[output:9a111cee] %[output:30277a71] %[output:715b88c7] %[output:39cc4efc] %[output:5f6b0367] %[output:7d8f13c9] %[output:2bc82a96] %[output:58a9f715] %[output:22145044] %[output:25ec02e3] %[output:5a6959d1] %[output:2c0f3800] %[output:39559fbc] %[output:54d3baf9] %[output:8751b890] %[output:61c740af] %[output:53c0d092] %[output:0ddec3ae] %[output:49a6e21c] %[output:4f3a1d87] %[output:43f9d36c] %[output:0421ff82] %[output:6e1f1c09] %[output:4d8efceb] %[output:5e00dcdb] %[output:3055157d] %[output:00702609] %[output:24f86491] %[output:7d793125] %[output:4184e867] %[output:00f6d864] %[output:2de978aa] %[output:3d11c38e] %[output:5a4228c4] %[output:4c98fc52] %[output:046d2c79] %[output:073913c5] %[output:30b61fcc] %[output:2ee3a0e2] %[output:7c1a143c] %[output:5ae0f6ab] %[output:66bd86ee] %[output:79e2321c] %[output:709570a1] %[output:0ace55ed] %[output:087ae413] %[output:424eb910] %[output:821a2de2] %[output:75ac041a] %[output:00c007be] %[output:590b7eb5] %[output:7a205f10] %[output:3cad810f] %[output:35368708] %[output:02a95d27] %[output:8a5c3afa] %[output:5c54bba8] %[output:7e87085b] %[output:3c47212c] %[output:658aaf53] %[output:57a3eb5d] %[output:79d2af63] %[output:55efc2a8] %[output:5bd0bf4e] %[output:213ad7c3] %[output:02ff7c33] %[output:093a1bd3] %[output:1dffe339] %[output:24130a6a] %[output:9e2bebc2] %[output:169cdf94] %[output:721fa7a9] %[output:51b7d3d7] %[output:6aba5b27] %[output:7938cb75] %[output:2d52858a] %[output:5fcb7ac3] %[output:024cb691] %[output:31200853] %[output:0468cd59] %[output:7e464482] %[output:8832dcf1] %[output:0ad9b5a8] %[output:0df9c6a0] %[output:7bd85080] %[output:633bf38f] %[output:9b5a7284] %[output:1b4831cc] %[output:51f06e47] %[output:85d711a5] %[output:0ed7535c] %[output:91c278ed] %[output:89361c39] %[output:4008e51b] %[output:57b2f2b9] %[output:6c05fecd] %[output:62c107fd] %[output:111cf61e] %[output:9ea26089] %[output:1dee75a6] %[output:26ec1fc9] %[output:1babf915] %[output:706ee142] %[output:61250b9f] %[output:37a490b8] %[output:1afb9e40] %[output:28078cb0] %[output:17f31a0a] %[output:6a9028d9] %[output:68cd74c6] %[output:8b149626] %[output:893c0444] %[output:1fbe7c79] %[output:4d6ccfe5] %[output:38f56e36] %[output:09f191d2] %[output:1db13420] %[output:4ff1894c] %[output:9bf7b9d2] %[output:31c9bc36] %[output:4606d532] %[output:298dfcff] %[output:8d6a919c] %[output:535f6997] %[output:08ce0c22] %[output:5703d082] %[output:4b32fe83] %[output:66999e9a] %[output:22dca6cb] %[output:24fa7434] %[output:2c4b0100] %[output:4e2529f2] %[output:0588bc36] %[output:16b71ef6] %[output:49897c8f] %[output:47843b1e] %[output:4842cb19] %[output:4158372b] %[output:3aa935a3] %[output:23cd4afe] %[output:2add4f6a] %[output:865dac48] %[output:07096dca] %[output:3e320d26] %[output:82571362] %[output:1a94c68c] %[output:4ca6c68b] %[output:28d29e9c] %[output:37edbc0d] %[output:5381b18a] %[output:3b1e27a2] %[output:734e04a4] %[output:5b2986ff] %[output:79cf1fa8] %[output:093a6d00] %[output:1c1c4cf9] %[output:46057be4] %[output:6cdce246] %[output:63b387dd] %[output:1722da1e] %[output:388f03a4] %[output:8bfce059] %[output:28fb07cb] %[output:0b80b429] %[output:45487beb] %[output:3a16e380] %[output:38eaf689] %[output:33ea1b6c] %[output:589f6e2b] %[output:48c93fba] %[output:843fb569] %[output:65750b31] %[output:08ba0b54] %[output:909675fa] %[output:722d7810] %[output:94c3ed31] %[output:64c9ddee] %[output:97462cef] %[output:9398f049] %[output:7aea6f6f] %[output:16dfb84a] %[output:2f0a648d] %[output:0622d5c9] %[output:78548862] %[output:8e912464] %[output:750ed7b5] %[output:16bf5d11] %[output:90b5d124] %[output:5d16e62d] %[output:8f6d345a] %[output:707e2dcb] %[output:4192cb1b] %[output:063abb5b] %[output:91fdbe84] %[output:1577f30a] %[output:05fc40c4] %[output:0df41ef1] %[output:06ffde5e] %[output:86fed95a] %[output:57344ad5] %[output:689a5670] %[output:812f0c88] %[output:904bb949] %[output:51eb1117] %[output:646c7645] %[output:5ac928ac] %[output:13803640] %[output:5ac8edae] %[output:1b6f8c64] %[output:4d8f4e88] %[output:245c330f] %[output:721ac259] %[output:13bd1f48] %[output:5e00aa14] %[output:5bc22594] %[output:0f0b47c7] %[output:484c5697] %[output:24ee5ada] %[output:6cad27ea] %[output:22f346a8] %[output:920c0952] %[output:4652c92d] %[output:453377dc] %[output:86f6365b] %[output:30da2edf] %[output:558314cb] %[output:57e02254] %[output:33302f5f] %[output:461fa8ff] %[output:619d5838] %[output:7e83ce73] %[output:690b0601] %[output:1f5e6ae3] %[output:375af545] %[output:5acd97dc] %[output:8f775984] %[output:3f3fd046] %[output:6ab4b36e] %[output:620176b8] %[output:9bf0bc5e] %[output:6e2e667f] %[output:69e8bca5] %[output:64699c2f] %[output:8e185787] %[output:5efca3c3] %[output:3e09cbcd] %[output:6efe8e03] %[output:55e4428f] %[output:14ea59e2] %[output:4721f82f] %[output:14c3e4a5] %[output:7eba562a] %[output:8c52dd73] %[output:8f4d9b5d] %[output:85aed559] %[output:8c247199] %[output:3d096e62] %[output:3783947b] %[output:073e552b] %[output:1f368d7d] %[output:84c20b0b] %[output:90a67ff0] %[output:70b921c8] %[output:08a998df] %[output:6aa2d867] %[output:7f1cabb8] %[output:70a8dc83] %[output:0be9e313] %[output:6ce5d7b6] %[output:09bc691c] %[output:8addb485] %[output:31c7c518] %[output:783f8641] %[output:4b0b5def] %[output:1b5d39d0] %[output:1e576b5e] %[output:3d8f2e07] %[output:7cbcd797] %[output:7589a3fa] %[output:33763879] %[output:46289344] %[output:95ff2c8f] %[output:4402cfcf] %[output:1ed32b5c] %[output:33feed3c] %[output:8388b5ec] %[output:08f0e7f3] %[output:57d5ec0f] %[output:7fc8e9ac] %[output:6269ef0c] %[output:97e4734f] %[output:2786407f] %[output:40006de6] %[output:8c98918c] %[output:5dd675f5] %[output:00d7c778] %[output:74ef0e45] %[output:49759a53] %[output:250dadcb] %[output:455d35ea] %[output:42055a3e] %[output:5068a527] %[output:349c7b4b] %[output:2815bbe9] %[output:070a9fde] %[output:2b2cd0e6] %[output:75939c1f] %[output:023f4171] %[output:5b81c394] %[output:7a080eec] %[output:22a89b3a] %[output:1114074f] %[output:7a2805fe] %[output:34df9843] %[output:7ed08f86] %[output:18f1bef0] %[output:3ce26393] %[output:6703a486] %[output:62a987c5] %[output:87e70052] %[output:9dab620e] %[output:726fccaa] %[output:78eb80fe] %[output:83bc0459] %[output:5b385a18] %[output:0ee466de] %[output:3175a2c9] %[output:2b495cd3] %[output:77305892] %[output:3ad1cd0d] %[output:23b2d10e] %[output:030445d4] %[output:3337ae6c] %[output:020b80f4] %[output:51f83c4f] %[output:39ad28ad] %[output:11e3dcdb] %[output:6f952a58] %[output:7c47ac64] %[output:85c1e2eb] %[output:2d052ed3] %[output:0760d073] %[output:7a788e0d] %[output:65b706d0] %[output:66cabff3] %[output:8352370b] %[output:1f4ceda0] %[output:107a2530] %[output:682a480b] %[output:05956480] %[output:4d7dc8d5] %[output:6224b1c7] %[output:5016fdfb] %[output:99dce67a] %[output:8a0bbe20] %[output:8129ec73] %[output:1afad042] %[output:583d47d5] %[output:368c5043] %[output:49407185] %[output:9d6b6953] %[output:4c507994] %[output:24e9d3ee] %[output:55f08081] %[output:048ef0cf] %[output:395feee0] %[output:8b15310d]

writetable(BRW_result_MK,'BRW_res_MK.txt'); %, 'delimiter',',' )
writetable(BRW_result_LMSlog,'BRW_res_LMSlog.txt'); 
writetable(BRW_result_LMSlin,'BRW_res_LMSlin.txt'); 
plot_10y_in_two(BRW_result_MK, BRW_st,'y'); %[output:6dac349e] %[output:8db59851]
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
%[output:7b24b1ba]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('tabular\/dotAssign', 'C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\datatypes\\tabular\\@tabular\\dotAssign.m', 508)\" style=\"font-weight:bold\"> . <\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\datatypes\\tabular\\@tabular\\dotAssign.m',508,0)\">line 508<\/a>)\nCannot delete 'expA_bgAE' from the table because it does not exist. Assigning the literal value [] to a variable in a table deletes the variable. To create a new variable in the table with the value [], use T.expA_bgAE = zeros(0) or assign first to a temporary workspace variable."}}
%---
%[output:96dde337]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('tabular\/dotAssign', 'C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\datatypes\\tabular\\@tabular\\dotAssign.m', 508)\" style=\"font-weight:bold\"> . <\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\datatypes\\tabular\\@tabular\\dotAssign.m',508,0)\">line 508<\/a>)\nCannot delete 'BsG0_S3S30' from the table because it does not exist. Assigning the literal value [] to a variable in a table deletes the variable. To create a new variable in the table with the value [], use T.BsG0_S3S30 = zeros(0) or assign first to a temporary workspace variable."}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8fe07ccf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40b2d253]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17d1c91c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f86febd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8e0ed078]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b911553]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:298ec69e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ee779b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:248e1962]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83bf22a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84060908]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08155abe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:73710e2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:018dd32d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:44d899d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:67388040]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22a309cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00e838d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85d54d6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:993ae445]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9b31a73d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36662d89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01d4841b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15201ad6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95d59e56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5cb5ef92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a70d0e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:86b42ed4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1fd14807]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34595fc6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:938ffa44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2dda597b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ada91c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:128525a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:024dd3c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9becd61c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:348171fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:340be5b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c73e2a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b500bbc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93b6cf28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c7e4d1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:937e65fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f12749d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14fedef7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f8f71b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43df0a2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32d5edea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1568a178]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:393188c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e6bd3c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c75499d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92f7e111]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32dc71af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87da0cfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:46e9006d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f549bfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d19fc14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:44484d29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13089376]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1008e3de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0dfad22a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e8c34ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30f45b9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51d0b936]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89715877]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0dc7802e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:10733847]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b710032]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d5d6928]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3753ca58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18baf1f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c3bab3d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65ff9912]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e9923ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:48540a2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8db3c968]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38ea5b7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32e97787]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ff26752]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42dc50e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41376bf1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:907db972]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7981fd14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b36dba7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34d214c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c93efdd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:36b55a29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d5bbf8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b16883e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7cbe1379]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d89d53b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9397a265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b83538f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9bd1108c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93686c5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:860dd88f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c489307]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:367decd9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:546b111d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51142b55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8213be0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90dc90d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65f88126]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c0d4dc4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66e55553]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ce0b691]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12129d59]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95063375]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b869b28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:031c6712]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7cc6f0d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0cc5f808]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:04817980]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3080e591]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95524bbf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54710bbd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d06bcf8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:077691bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:29b4bdd7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32035ee5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01853600]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:658072e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2801f0a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:372b2fea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5554bccb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:506bae5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9126c750]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b24bfa1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89a59519]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3494b737]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d1feeca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ff4a3ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6762415e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e35ec14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96cdf2f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5adb1609]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89bbea40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37c3ae91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:67d65a07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:196d887a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9809bce7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:79fc8e81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34c31025]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36558d63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6688afbb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4135bf82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26c9cd1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2895fcf6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62d2836f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ee813cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:142b0a16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:67f38a16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34ac7265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e51d712]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64fe1588]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c2330cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2a4c696d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b30af8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4fbb3e67]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04e78057]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b6c3002]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3024b889]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43a1badd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:481b0c93]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f3232c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4cb41de4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:213e149e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99a27963]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59bca9f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1edc1436]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:278dc67f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:612d2e14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a31ff0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:137e2f57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96b7b375]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:175c0342]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7bf2a46a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:325839f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03ac8eff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9c6667b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:795d3d1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22b85d3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:918473fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11001239]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:984e4182]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:639ee66e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b67ce71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e78f5c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:56059f9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87e9be89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7141f281]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4afc3452]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2623bc39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7fc0485c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4fb6c006]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8561cd64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20d39efb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11431d9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:048e97ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1852498f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9143a6f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3380886d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47ffb32b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8588c821]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:135e38ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:454c7f7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:76b6c351]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8bb63f0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:94b3a58d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ce8b5ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4bf4381a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79d6160f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:969ef301]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3773d76c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6dd507d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1cf8bdfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a595b1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c96b513]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42879679]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18d08194]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ec79177]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:27737cda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b8b3559]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7d2f958a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93559b24]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:470c9215]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59c56c76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5237cd0e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:508a7a92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4255bcda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7ec81959]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ade710f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5443e0d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6c9e2418]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:457f4c1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d574454]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0db7d784]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:454b51c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5c73b87f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:39652ade]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:74357a07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b2412d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0903c4a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:234929f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:070a795f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99479e13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:103d094a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0404a91a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4fa9468e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0bde47ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85260890]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:64f5ee8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13139e4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16cbb4eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21ba8c85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3d0fb04c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f126aca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a2da152]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:311a487c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2170f141]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6bddefc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0da5f4ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0817b4ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ee7b693]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a52dac7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c16d95e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5195e48f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15c2e6a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d9bd20f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:198efaba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80284d00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:72e2ba8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1835e793]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4531b12e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c381965]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:224b96be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:635c9b22]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:445df048]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c09b771]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80cf8149]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3daa6ac9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90bb9036]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:758fd6f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:753ef725]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7094ccf5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ac20262]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2169675d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:187dec53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38b46c11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c1cb559]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7350707f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:66461452]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4de3aadc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:008a21ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43ad4fa1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4edaac14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f14ceca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:883a501b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c12e34d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7cb08fa9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:469bc765]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34071e90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0fadc2a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97f618c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e270ebf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15b7c972]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d29dfc8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a45d3fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c48b797]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a26972b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98a0eeee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e52ecfd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34858ae1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a5ad6af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ccf8032]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29198aa4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:453264d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6652ff6b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:663158e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:833d3a1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:52afc33b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9965b3a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65f3b8c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b734d75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8aa73ebc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99527e6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c00a7e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82bf6ca1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83e751ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:959448c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6046acdb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d234253]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ab8be3f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9cfffedf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:91fbd516]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ede8235]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:091957ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:04d88352]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c6fb989]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f5e9c1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72da9bba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1a9ba36d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42f90fcb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1e8345eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15513712]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e53c04b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c92f5a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30b4837a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:597192b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:73cfa7ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14a98fd7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19e6b3fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:70003f97]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06c51ea5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c8bd1a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1a0bcace]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90364ea7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0251823f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40bdd972]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a040431]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6acadad0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63115f0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0801b24e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7dea8a02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62a11b41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f4f19a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0608e2c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e71f5e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7dafc6ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:456551ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:38a604db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:960579c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f59913f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cc418b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0c4e9014]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21cc769b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9f24076b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:636f52e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0d7ea755]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e89a2d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2483172f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88df6649]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6334f9b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:330ad32a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2dc7fc8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21c0df7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8bb5f7ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3051a516]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:336c33ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a085ee4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c466ec0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:015cb391]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b57dd35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18ac3c72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0da4bf13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30acf795]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ac12db5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e21333c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3205fcd6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:721da8fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d366bc5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:983aa819]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f0b16eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88db5cc2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a51fa74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75d9574e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f438b82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f331783]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9dc5979e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:235cad5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58b620ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b41c513]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ddb25f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ab2fc80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2f71bef0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37cf23e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:295c7b54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61077f08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e4725a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:50878f63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7abb99b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e23b6c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4278e853]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c936736]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26a8f37b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8afe35af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f7c22ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:50dc6dab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d6c376d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d9c5eab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49c72875]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84ee9072]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:218cf845]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d1a72fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8ad5b377]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48e441d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:365b2f9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36d3706e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18b338aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5cf97442]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4863cf71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c592f84]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34dd8e79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:586a6115]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:488fb4b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87fce924]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b207c75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3675ed46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43170eda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3251e38d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66d74be0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:68b65577]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7384198f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ce446de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7698e34c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4474b738]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:191152a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:122b0f3b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b03e27b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e4b2265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e2ff7bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3965d6bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ffb08e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c147e00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34fde0c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:21ae4975]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:10b9104d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6566c6a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f904590]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:79bf07a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0186dc17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:363c88f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:36b180fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58d87de9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8387a974]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:52c9e2c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1831ac8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66b7fcea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61607279]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:682774aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:28ced20c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:752a74f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b516c0a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4094f7b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:506ff6e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9336c533]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8a8a2ede]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:967d09b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d5180de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74344efb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:21aaddf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:71e11d2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d10f7e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54be5c90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:660ef67e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3398924f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:169784af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:56970fba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06b0e977]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a77ba81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4722ddaf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e6de565]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a74debd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:842a63d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d183d89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40de2c29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:475024df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3a5a4286]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76076ca3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8912acf4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0aefa498]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2532de9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6fa36ea3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ceeb71c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80b37eb0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9417dabe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a8c7a5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87b7195e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b698fa3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0dfda2b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c5ac0d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c185aa0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:053ef558]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f12dde5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e288b62]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2edbd576]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:714aef69]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49946398]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b78cc28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:803b4b5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a4a7201]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9327ad9a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1dabecb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:941aa8b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93e78f5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1448401d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49aab42f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4664b6ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:387f045f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24371ede]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cfb6040]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3537b91a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0dcc4ea3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d5ffc95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d81dbe2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b080ef5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51290129]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58623409]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47d2303a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b284375]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64c0ff0e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:594ca13c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:624d246b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b0f3cfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89588526]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2505a68f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e72b5f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2f950a37]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49dc6a9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:52f25e5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92e68211]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b448e58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9f7c6739]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7443477a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8027f5df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0a2a7959]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9c60dee2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:91c3db57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f71c328]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:20a37044]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d127a48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a03cd28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5aae1920]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:77fb2d40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18e5ffad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39dfb265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5caeece6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c8dc615]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:07fd2a8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d4e90cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a545acc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f37d48b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9f6022a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:275a974c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f36175b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5906e80b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:887499bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c4aea40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2793d475]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1edb9035]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c7ae58d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3eb66204]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e363b82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e92393d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11ccbdfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:587529f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:570e88e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3900e64a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3575c0f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ad760b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38cf78ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0225c240]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a0488ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4698c34d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89f194b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b7682da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31754962]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2fbe3fed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e2db140]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12f9ea2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16bbe101]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31342d67]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e099250]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82f43d66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:68ac73ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34721f8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5942e1fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:006fe226]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22545626]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48cd084e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:585bfd8c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f6aef5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e2d6192]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7557f411]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95b6eee4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:334cb600]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20d579e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:254a64c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13e374fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4dd2fc4d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0855a80a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89775020]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4faf8252]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:792b81d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7155133c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58181a86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c5bc103]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e1cfb09]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1924c3bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57f5f2e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95d5428b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8bf89b19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ddc9a86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6175de95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39fe8ddb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:27c2871d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1837d5e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f25e7ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c2d92df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:554cc82c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:127911b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29ac48cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:574e224b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45a4b885]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7de238d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c7678b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63c95266]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e8866ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:259a3e07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99abc3bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87eaa7ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:85a3b8df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:292ea369]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92156ff3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2706e187]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9fc54de2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:597c3bd7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ce550ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b2bdb07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5fb95bb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c1f91d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88ecc3d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ce147c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f2f7ac8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93fe0ded]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b9f216d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f667f47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5e2fa66c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:337611aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66b4883d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8721a920]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d1581b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60705e91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1cf2a9b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:516360a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:938e7d01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59b8d9ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80e8f246]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1950cda8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a480775]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:976ba4c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:571549bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d5ebce1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:822f0fdf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7bfd20d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82c8b830]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a6f7ea2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34e17c9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42aff600]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2808f36b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c1fd27c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c83caf3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5be3ca54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:88fda62e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77fb56cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bcf67b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3a60695f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c697b5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f3bbe72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21f632a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9325c816]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4459f20c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19803d9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3aea5ef7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d17b0a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32b06ea8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:831bf892]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b31a37d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2cdb3459]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:136f79ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5490e1f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:797550e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64978c26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54c675cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:46662cc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:457a8a73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43bf38ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:770bc4c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f6ca1e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b2cfd90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d2abde9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93151d9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:297189c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:134a2def]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:306ba9c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65d6e7fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03818a92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43c3bd7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20ac780d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64547591]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0828ae10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a559b7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06b92d0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9047c90e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1fe031c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d2cbf9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10e22b02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e8e61eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47b6c75e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b262e29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d38f3af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6d6ffd01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80e4e314]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32b644ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d97f9c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2941dc63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9c937f8c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d577d2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:81d5fd03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:98c2f4af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1093bee7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f2679e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74ec2644]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69d38380]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b941da0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f788cb6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f4e71b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12f2f3c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:67d4649f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e0cdda0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b6d3039]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:665d457a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76148732]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9fbaba72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9de1d054]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99598631]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99a888c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1935039b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2764cea3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d9e61d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:761d0d20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f61426e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43474f7f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0677bebb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7246af85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:137ff240]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:07303485]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63c4e12f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ceac80d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1623d026]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5862ddd2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46976e3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95cf5821]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ebeed9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:374a7a19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30c12796]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e41aadb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2972a2dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98fa4d2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72571873]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:50a01290]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:326b8715]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72fea0ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ffeb845]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43e83f4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75ee40ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c99078b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9097e384]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:812dc5eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1a62274f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43eed418]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:204680d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9775e68f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0c3855df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:25e33875]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e806063]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ea7c506]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30d15316]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:495a2236]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10d7acc6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6071fc59]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d614f7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e7cb58a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c48a3cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37ce7784]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:583e5487]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:343ab644]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:497623aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e48d93f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5cd7547e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d0b0ff4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e54eab1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b723b65]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2047c511]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87b50185]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ec126cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01edad82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:01f13939]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e91b473]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90426e8c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2db378ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6eb51b03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b7cac90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27b08062]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:53448ab6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2f5eae9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d4851d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:68260f97]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00433012]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03a2d693]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48128d6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3042d64c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a5e9051]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:23c19e0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:420497b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7452df56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a736a82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:76fe587f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7dfde792]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59cbd61f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40f048b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18aadaba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:413b3f92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d10677a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32c834c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e9b45ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ac1acb4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:925368cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b7b5d88]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ad62ba6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2117c05d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79165849]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05ab930e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:226a1551]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f13862e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:505359ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30a55861]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9df1eb23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05db6a0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:50a68c5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58bf8f3c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:494a2544]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c7fa512]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88439206]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:402cd7fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6d727988]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e52e816]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c3c05c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:225e283c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d075d6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a852b34]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9dd8819e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9416bc8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a21e657]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17bf66eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c9f6557]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7deb8f7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b3568db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0bd65661]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85f53734]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:387e1835]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:50f99c17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c2c9cce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79f2d8a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1640effb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30bad384]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:622f2875]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2a6db916]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0bc7e39b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1d9aa334]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6822ab16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d78aa78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:513b5870]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16d85161]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b8a08dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d447256]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:01b7c5cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:529e6bf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a2bdb55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e8aaaf5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4db64809]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:23883739]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57c19e61]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:761647bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13e0d75a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1307ef03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6960a139]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:642d79d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4817e881]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b411a3d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32d574fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d5d6f3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80b83034]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cc388f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b58a23c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:351821d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:799097b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:799a2222]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:751187aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93c84f38]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ab507e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27fc40c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0cd7f1a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ce166ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7af87bf1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a0a17ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92ab2c20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6288e5e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a6e18f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a176729]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e628886]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:659c02bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33ff1abe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15b24b94]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4dbe7346]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26d28852]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14486678]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c8dbbdd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:25446015]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a26ff36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1079d9d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:532b75d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:594b51c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e030c7b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24f94620]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f22a0e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f875c39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:23bc46ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89e084da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60680960]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95bbfda0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61c2c8ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f0eb7d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:618c5168]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:14be3dfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38627e4c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d726a71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2432710f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:29d789ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1a81a99d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e5ba1ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3fe4140f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22bc5f8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a2be139]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98b7854f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97daa11c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3de97246]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90d3f02f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3169c5e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8cab4c98]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f8cacdb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:669c3baa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75abd782]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79558ee9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:878b3c2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39376864]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:494dc585]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7bf9cbe3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5c76acaf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ee6f087]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:35dffbc2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f7a0691]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1771c5f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a52888a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:588a5069]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1094d469]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0bddb1a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:069a3b96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:96e85c40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42e7542b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0dfbdc8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a4f5f2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:489cd351]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e31f403]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1341bf5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:804f7253]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40035a89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0670583c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d2b7df1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92114223]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ece0d9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1fc1adf7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1285432b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33a705b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d8dd318]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1a765f3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f30baa3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:254a4d00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:35056cc3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12083e87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:323b064b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:855404f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4655327e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b0d5cbf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45eb5feb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e7b94db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:56a167fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48fc1d9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:002181ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:53c4d703]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:96c83a29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1cc98398]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2f373d5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8dfa3675]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f7de47a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ab75b05]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:78d01d21]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13ccb897]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6390c0fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:28c1feba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:079374d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5da4659e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d882d06]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6df06c49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:64f922dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9207c7f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7ffa277b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a9daaf7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6339c644]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:78bbd259]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:02d2b3ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5845d840]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:350b65bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2684b348]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ecb2ca2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89f08391]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30ad95e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2cc87424]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f28d897]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0cb25264]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2482bc12]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b90400e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6742217e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:02656923]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e4ee37f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:09f121b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:922375a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f085ef2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:453aef2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b8de5bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5afdf3ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8a3f2782]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1da5c1c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f6183e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8283c67a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:29ae0c19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ca58c23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9fbe16fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f25d80d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8cf93ee7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e8f14b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76a0f9bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9fd6209f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61d6acd8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3a8cce0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:197dc401]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44bf8509]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89e0773f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:348f40ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64fa1a60]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:81b04c90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51da59dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:25b6a4ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3732c93f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:142feeb8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60292b2d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:646cf396]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72a495af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5de2c74f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1781cb6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01691c52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1be70811]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:847b0d16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f97abd9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47dae02b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8677fd68]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9cc0d729]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d61aaf4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58cee9be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01fcf444]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:194e2fc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:52a95ed6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c5328bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:993c1df4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:790ddd1c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:706628e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55853ada]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5cd73f0e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1a1662c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:471c26ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27f30a75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13eabec7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c8f710c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99a6eb62]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e32c918]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15c9d771]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7752e6ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9721c460]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7ec915a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d631b47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:581aec3f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6974ba62]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:71bd74e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c9b23aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ecec322]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d9f796f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:571b729e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97e718cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:72485a3d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a327883]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:419dedf1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88fbbca5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:523e9971]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:86795c12]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:389fc8e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0eff3809]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35bcab0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:839c4d90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:485f12a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74a60b48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:472bca99]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d8b7e28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:23757618]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a24609d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6eba130f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:008b2ee0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93b99ca4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b81ff2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7722cfc0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a7c87ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0df96489]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:166f89de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ec7b580]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82e38465]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05c35e82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e24637d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:649e55f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5985f94c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03d9874f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00b802f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:515739af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e37efd3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06e48b48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1307e453]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a52dc21]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:78c32ff0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21bbdd1c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4709fb01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c6945db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0dc8ca2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b0608da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c5ec070]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:851ed762]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59b338a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c484356]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0c56a874]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6aba0c4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4381e90c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:744b45c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65d4d8a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40e6b0aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ee4043f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1227cc3c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75aad05f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4df85a43]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:444fdf86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d906467]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8a15888e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9adaf37a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b8efe90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:821b0184]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6471534b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:95fac716]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:522e6250]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e36e552]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13dc14f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51ed18bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:45d7858c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6410adcd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6c817e89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d3f73c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:071ee91f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93760aaa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:07e0d602]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5947e9a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34d35ca0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1bc72563]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:31357896]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60fbdd18]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:874a956f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ac32125]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c425f8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7828c641]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c925f22]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62be5cc2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:579a0b45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5675d56b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:297119c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:06f94e52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:239bcabe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54f79a2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b603af7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42791e4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:992c3f08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64c653d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b530051]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20bcca1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:11bcec20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05d75d5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ced4381]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c3a5ff1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5e56df8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43ee64f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97434ce9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49e68588]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59786baa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b1e622c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:07042ed2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8658390a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a7a99f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a37aa50]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:975f28a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d2d6a85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:278b81a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6dc65d00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b493abe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:23f79076]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:89184c2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:407c20c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:177860ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cf3f165]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:78fa46bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93fb6416]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:94863d74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11ba9b4c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f97591e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66ba48f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d453ff0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5910b602]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71fcc425]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92906f6b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40c6b934]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27de264b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:145b9841]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40b9fc49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:976ecd41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:009b9d83]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:62f64890]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b538123]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8efabc44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4fb948ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:53c2b438]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92a51c96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f8aae92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26de4ece]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9814d4b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:670bb8ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:244b0e9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a665ac9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:86aa6b35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:91ed1a99]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:395d78c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:799cd3a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5630bb82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13743b4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a4199a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34119526]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4479e3b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0beac469]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a10c6bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5db1a1f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:614e2e40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8e06d8e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22c53bcf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c96f001]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62ace964]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8e1c68bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:8a8c186b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14e74ae1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:69854e3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55704137]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:63dca84c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39180e8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7856fa49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:81b8313d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2586ba6b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:15196897]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:641aa428]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:07828f93]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2814c7a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0f5652ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9cba3311]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1df52f65]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48719809]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:700678fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:293250da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:20789b6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9fd94820]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03d062d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:164d12c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2131ddf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2be64074]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:578fdf12]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:68e7d81a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93bea9a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7774f180]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:495fb57a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c57f034]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b24147d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24369605]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7db1dba9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:80e5033a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1043d283]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51176c6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48b94a28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0f044448]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b44c548]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a7c4199]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31a4b364]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b0787b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0a83c526]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:773d97aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47c2d859]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60da6422]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0c9f3a0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4cd27df9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:379436fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:68d7579f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4a11aa6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:765817dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88ba3736]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b539c0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88cc8ac2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:55155896]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cc8d5b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1759c647]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ecce73a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57f9c687]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f4789b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f682b8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d5d1de9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b416b7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:270240f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:723c5045]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:86f73dfc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:33e84adc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60d2ced1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0bad9679]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26d9778e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:699a18c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89d21e33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2dbb107d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f729257]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1316a092]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ffd965d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89a5fbd4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58acac1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ca9239f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7988adf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45e4edce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a1f824b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:959819b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92120768]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49f70b33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:67432884]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:895c42a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:812df8d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16e488e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41f2923d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f334eeb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5eeff1a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ab71e87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:07e28390]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22038812]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:736acd05]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:70f830f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51b12777]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:35bf91a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88828cfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1db79aaf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:68583dfb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:23f92763]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f6e0587]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:066e3c5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99a44186]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:572f7343]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34176f95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0d287aa0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c7129e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:149d493c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e4f977d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1bd17cfd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2856ad6d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:592dc527]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f471478]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95738865]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c5e54ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b92eefb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17d8bf04]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2f99ea5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ff8ae3f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ff4fcc2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5fe4665b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9b77d868]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ce3a7a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b8b01e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:200e7cb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b7369f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0dc3c9d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39dec4a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7953c2e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49693945]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16422f80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:70430658]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cae4c1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9edd77cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6127f0f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ef2719e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:279542f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:50a8d3bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0bcd6429]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4beeb40e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77648f37]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59e0d2c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:580b4dc1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5fa5857f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c977cae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8eafcaf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ca64294]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72c65301]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d7f94b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04130bff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e30a1ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18e9e11f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27674809]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9d7002d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:418f95ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c811d35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ee28afc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f6f4d48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:96757d87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:04b1952f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9226e66e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:928c2ff6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f277d4a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:940f7643]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63fde13e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1b09f981]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:738a7a29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a1aa569]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ef87133]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3279a919]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ae29569]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03c972c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9f1c3495]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f1d3331]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ea129a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30cf2c9a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1d783524]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:663e1ee1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:52a01cad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b0da752]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32c370c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f352396]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0be73d84]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2825f824]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:455165c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0014d4bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03504aae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7669f20d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3fa3440e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f6d0f3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51f8cd54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75ce3df8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:786c44c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:91e365c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80b86e44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2eb70ff2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:014ecd58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3369c2b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:815b35d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9452b033]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:150c3479]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b4347c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:35a8ac56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b890b43]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97331f43]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7ce9c94d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:481f71e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3371e760]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9b850a5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e5c4886]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:25ba1ac9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:938ee6b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:064e3e42]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:664dfb89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6217d75d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51b66828]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:222b49be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:86d9d598]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55cbaaa3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2feac792]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8340fa75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:086fc12c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1fce6c76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c7c2096]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59de9faf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ac47599]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:539f56f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6252c30e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:368df944]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c333f62]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:679b3d76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a5837b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:345c63db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5405f067]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4da1c54a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7378e65c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:39a03f01]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5265e380]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76aee1eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55de75f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:864337b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e399c31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0e2c4766]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b3dcc7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47176382]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:929e2be5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b5456fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:699f0ac0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1d6e9e86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ea78f57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6da35736]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:88d0fa8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f4c0949]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e55e069]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4eccc0d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22c2f0ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:77ab0f69]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61b578e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4418a538]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:577b3ad2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06185cef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f508b40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:586897f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77f5e924]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:39e8a051]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cd885ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:907e2d0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:206a6d7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:928d7028]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ecd7a8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d583c4c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:930a3f1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:285110be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:120ecb89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0bfe3b70]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f2147ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6c45117f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82dbbd37]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:544a9231]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8be28180]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0bc8f7c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1fee6824]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4072a0fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:505367bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a9bc25f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:542fe636]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:934ab92d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:917aab98]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65a6f6c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16fe68b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2fe359fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7067095d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97db2197]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8573531d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:361314f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55f218c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:486ba7ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:70b0a8ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8296999d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:453bfa31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c656203]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f271821]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0be699fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d23f8e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:201d73d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b67e889]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9c87c3f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3bf84064]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f92fba3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17540ef7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42504ecd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e31910e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e024ea5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0abbff9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2919cae3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83f77888]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:478ec03e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0f37217c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:022ae550]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12032b4a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05ebef81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:655bcabd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05820577]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c6e7405]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:415a0e41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:715bed63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2bf5cf70]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b31e1e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60c616ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:78163ffb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:79209768]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:297a5b07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:806cb1a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a42be28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:94612da5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6c64e8f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98f9dcbe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6eb97f16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5de2cc45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44cb88ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d8e6622]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5aa2a68b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8cae9f00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c9f7998]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d7ee1ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:221d3d10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:343ac351]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:308a2c72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:552567d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:689da17a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24cc3bc6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:36943925]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c0ca563]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:976ef4f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f5060a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3371b0fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47aa55af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:016220bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f616a67]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:88404826]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2220331e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:849a735f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:506c223f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0afabfdc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7d9d64a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13c0459b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ae7143e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30734195]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80562a9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b3dbc4d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61e2410b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:41311589]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:42c1bc90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0de4d03f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b39fca8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c95a6e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b12b34d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5fc18053]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b5adb59]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:020647e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66698d0d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80c380d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2e23abdd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:806fad47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:623ed63b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03bf74cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84dcdca4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8760e01c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1145aa2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f7c03ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99b1b93f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b68b300]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c12e1d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14aae969]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b601021]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:084d9598]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:19a8ed18]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65a49127]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:611152ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4a961c66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64d5fe13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22ab7464]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b5c0dd0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99c54a1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:91bfb35c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79611708]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:006b7de0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61a15a78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01d57ae2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1205d11e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d5c0c10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35206373]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:5d7defcd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c211be5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3be6e255]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74414eb2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1e431866]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4dfcf25e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0458d482]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:359c584c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5667cf60]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9687c67e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:33528839]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f0c8b52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7da91a73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:179bf83b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:02895f06]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9677e3c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90a77b41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ed3464c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66601631]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:176a0ac5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:31ca3fd8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b6c612e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:286834a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9189f8ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71d55744]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:51463443]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9efa67fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:83f712c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f8609ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:991f70d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:3637cba7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77d350c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9d172f85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5daf6efe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90fa2ae6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:16ebe48a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76295b78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6406f22b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cf51fb2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d740cd9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:878dcc2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8203641d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e78cad9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:570cc0a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c7d5b4c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9a184af6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f8bde14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:746c199f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:316f8a74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7fed47cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cbb0e2d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a51ba0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24903316]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69491921]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0994489c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a8c7f1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46a0cf4f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0922354f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e17a5f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:16f186a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d27272b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2cdad14d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:733d9ed6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a6bcdd4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ff773a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:563f32a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4acf3a67]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9c708e0e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92d63e72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6070f6f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:711b6604]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10ee5620]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77769871]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4be26cd8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:73f0e7a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:598219a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b58b3a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7d5cc8dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6d794174]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:106e1883]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84f5e974]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f6cb2b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:908f0d98]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f994b9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36f30acc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b43684e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f9adec4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b8d9e1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8383f4cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5705358f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04412afc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21884795]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ed89835]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:263274ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e5a20c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b836e91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e8f1cdb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8cf78d7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:46ee256a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d9a0261]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1eafd698]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:739af566]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47d0aff7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41beb7f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c1e3a4d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:564cb46b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:40b12c44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ff329d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:502bfe75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:604335c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75988d4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3a99c051]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4395a1b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0b32b2f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:885f9a17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03c23071]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:67d8966c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:02cc1c29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2a006c80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0343c1d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4feaec59]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7eaac0e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3a7a6647]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e176b3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d30311a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9938efa8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:590e3702]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b346380]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d40d8d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64b2a5c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7255f7d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:251aca26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:699921ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13a5a259]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6769b710]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e57fca6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9f2f05ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9f15c02d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:214ac371]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0582a2e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60bb2e9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f12491e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:102dac56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:474ead07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6377e7b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79f94a8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49238d14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b8657ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83efc422]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c02e4fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82fc39d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:018719c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49d53b42]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61f32cee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:676c34f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b92e358]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4300a1b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8852d368]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58d48393]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92888563]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d270f07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ede0ae1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9657f593]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0db46d5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ace1df8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59c165cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65954e35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:92f3f4a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bc8cedd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:68005111]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0fb71617]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:838bae79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29ecd218]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13bc4764]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89542deb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:870df706]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:2e0c0dd5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7cd00b02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41a19d04]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:476d05c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:592bd887]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b1c18e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f4720cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:865ad4db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4a29de84]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:3f59608f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60fcd37e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54c9233b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:402f46e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b3ed03c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75952393]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:600d18ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:386244e9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04f8b022]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:356d799b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7887a45b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bce111a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c76a8a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a1aa8d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5554110c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41a6ff2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2167e7a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44b3ad23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08c52778]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c076e3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80fa801e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3ece19b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:72cd91b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3dd42aa4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:48b75621]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5668f291]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99da434b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22dbe781]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4920a59e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0815c34a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b4c0a1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:93eaa66f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a111cee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30277a71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:715b88c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39cc4efc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f6b0367]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d8f13c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bc82a96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58a9f715]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22145044]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:25ec02e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a6959d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c0f3800]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39559fbc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:54d3baf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8751b890]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61c740af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:53c0d092]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ddec3ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49a6e21c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f3a1d87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43f9d36c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0421ff82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e1f1c09]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d8efceb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e00dcdb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3055157d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00702609]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24f86491]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7d793125]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4184e867]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00f6d864]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2de978aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d11c38e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5a4228c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c98fc52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:046d2c79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:073913c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30b61fcc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ee3a0e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c1a143c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ae0f6ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:66bd86ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:79e2321c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:709570a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ace55ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:087ae413]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:424eb910]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:821a2de2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75ac041a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:00c007be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:590b7eb5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a205f10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cad810f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35368708]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:02a95d27]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a5c3afa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c54bba8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7e87085b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c47212c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:658aaf53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57a3eb5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:79d2af63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55efc2a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5bd0bf4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:213ad7c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:02ff7c33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:093a1bd3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1dffe339]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24130a6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9e2bebc2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:169cdf94]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:721fa7a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51b7d3d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6aba5b27]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7938cb75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d52858a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5fcb7ac3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:024cb691]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31200853]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0468cd59]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e464482]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8832dcf1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ad9b5a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0df9c6a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7bd85080]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:633bf38f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b5a7284]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b4831cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51f06e47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:85d711a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ed7535c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:91c278ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:89361c39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4008e51b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57b2f2b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c05fecd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62c107fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:111cf61e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ea26089]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1dee75a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26ec1fc9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1babf915]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:706ee142]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61250b9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37a490b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1afb9e40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:28078cb0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17f31a0a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a9028d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:68cd74c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8b149626]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:893c0444]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1fbe7c79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d6ccfe5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38f56e36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:09f191d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1db13420]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ff1894c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9bf7b9d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31c9bc36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4606d532]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:298dfcff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8d6a919c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:535f6997]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08ce0c22]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5703d082]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b32fe83]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:66999e9a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22dca6cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:24fa7434]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c4b0100]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e2529f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0588bc36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16b71ef6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49897c8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47843b1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4842cb19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4158372b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3aa935a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:23cd4afe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2add4f6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:865dac48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:07096dca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e320d26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:82571362]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1a94c68c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ca6c68b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:28d29e9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37edbc0d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5381b18a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b1e27a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:734e04a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5b2986ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:79cf1fa8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:093a6d00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c1c4cf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46057be4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6cdce246]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:63b387dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1722da1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:388f03a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8bfce059]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:28fb07cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b80b429]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:45487beb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3a16e380]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38eaf689]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:33ea1b6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:589f6e2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:48c93fba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:843fb569]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65750b31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:08ba0b54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:909675fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:722d7810]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:94c3ed31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:64c9ddee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:97462cef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9398f049]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7aea6f6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16dfb84a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f0a648d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0622d5c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:78548862]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:8e912464]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:750ed7b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:16bf5d11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90b5d124]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d16e62d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:8f6d345a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:707e2dcb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4192cb1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:063abb5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:91fdbe84]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1577f30a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05fc40c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0df41ef1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06ffde5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:86fed95a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:57344ad5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:689a5670]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:812f0c88]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:904bb949]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51eb1117]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:646c7645]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:5ac928ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13803640]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5ac8edae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1b6f8c64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d8f4e88]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:245c330f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:721ac259]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13bd1f48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5e00aa14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5bc22594]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f0b47c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:484c5697]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:24ee5ada]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6cad27ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22f346a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:920c0952]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4652c92d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:453377dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:86f6365b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30da2edf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:558314cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57e02254]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33302f5f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:461fa8ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:619d5838]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e83ce73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:690b0601]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f5e6ae3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:375af545]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5acd97dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8f775984]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f3fd046]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6ab4b36e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:620176b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9bf0bc5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6e2e667f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69e8bca5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64699c2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e185787]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5efca3c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e09cbcd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6efe8e03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55e4428f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14ea59e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4721f82f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14c3e4a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7eba562a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c52dd73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f4d9b5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85aed559]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c247199]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d096e62]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3783947b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:073e552b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f368d7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84c20b0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90a67ff0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:70b921c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08a998df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6aa2d867]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f1cabb8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:70a8dc83]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0be9e313]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ce5d7b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:09bc691c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8addb485]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31c7c518]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:783f8641]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b0b5def]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b5d39d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e576b5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3d8f2e07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7cbcd797]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7589a3fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33763879]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46289344]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95ff2c8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4402cfcf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ed32b5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33feed3c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8388b5ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08f0e7f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57d5ec0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7fc8e9ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6269ef0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97e4734f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2786407f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:40006de6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c98918c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5dd675f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:00d7c778]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74ef0e45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49759a53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:250dadcb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:455d35ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:42055a3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5068a527]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:349c7b4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2815bbe9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:070a9fde]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b2cd0e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75939c1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:023f4171]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b81c394]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a080eec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22a89b3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1114074f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a2805fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:34df9843]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7ed08f86]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18f1bef0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3ce26393]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6703a486]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62a987c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87e70052]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9dab620e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:726fccaa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:78eb80fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83bc0459]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b385a18]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ee466de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3175a2c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b495cd3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:77305892]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ad1cd0d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:23b2d10e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:030445d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3337ae6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:020b80f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:51f83c4f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39ad28ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:11e3dcdb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f952a58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c47ac64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85c1e2eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d052ed3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0760d073]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7a788e0d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:65b706d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:66cabff3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8352370b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f4ceda0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:107a2530]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:682a480b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05956480]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d7dc8d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6224b1c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5016fdfb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99dce67a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a0bbe20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8129ec73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1afad042]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:583d47d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:368c5043]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:49407185]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9d6b6953]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4c507994]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24e9d3ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55f08081]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:048ef0cf]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"8493×12 table","name":"BRW_result_MK","rows":8493,"type":"table","value":[["'BRW'","2025","10","'daily'","'BaG1_A11'","'abs'","'y'","'MK'","0","-4.9905e-04","5.6204e-04","-0.0017"],["'BRW'","2025","10","'daily'","'BaG1_A11'","'abs'","'MetSea'","'MK'","[-1;-1;-1;0;-1]","[-0.0029;0.0005;0.0010;-0.0004;NaN]","[0.0014;0.0022;0.0029;0.0047;NaN]","[-0.0073;-0.0008;-0.0008;-0.0054;NaN]"],["'BRW'","2025","10","'daily'","'BaG1_A11'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'BRW'","2025","20","'daily'","'BaG1_A11'","'abs'","'y'","'MK'","95","-0.0020","-0.0016","-0.0024"],["'BRW'","2025","20","'daily'","'BaG1_A11'","'abs'","'MetSea'","'MK'","[95;95;-1;95;95]","[-0.0034;-0.0009;0.0004;-0.0077;NaN]","[-0.0010;-0.0001;0.0011;-0.0045;NaN]","[-0.0059;-0.0018;-0.0003;-0.0110;NaN]"],["'BRW'","2025","20","'daily'","'BaG1_A11'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'BRW'","2024","10","'daily'","'BaG1_A11'","'abs'","'y'","'MK'","95","-0.0015","-2.5079e-04","-0.0027"],["'BRW'","2024","10","'daily'","'BaG1_A11'","'abs'","'MetSea'","'MK'","[-1;-1;-1;95;95]","[-0.0066;0.0006;0.0007;-0.0057;NaN]","[-0.0015;0.0021;0.0029;-0.0005;NaN]","[-0.0116;-0.0009;-0.0008;-0.0112;NaN]"],["'BRW'","2024","10","'daily'","'BaG1_A11'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'BRW'","2023","10","'daily'","'BaG1_A11'","'abs'","'y'","'MK'","95","-0.0025","-0.0014","-0.0037"],["'BRW'","2023","10","'daily'","'BaG1_A11'","'abs'","'MetSea'","'MK'","[-1;-1;-1;95;95]","[-0.0027;0.0004;0.0004;-0.0071;NaN]","[0.0032;0.0017;0.0021;-0.0015;NaN]","[-0.0089;-0.0014;-0.0016;-0.0131;NaN]"],["'BRW'","2023","10","'daily'","'BaG1_A11'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'BRW'","2022","10","'daily'","'BaG1_A11'","'abs'","'y'","'MK'","95","-0.0031","-0.0021","-0.0041"],["'BRW'","2022","10","'daily'","'BaG1_A11'","'abs'","'MetSea'","'MK'","[-1;-1;-1;95;95]","[-0.0001;0.0004;-0.0001;-0.0173;NaN]","[0.0067;0.0024;0.0019;-0.0077;NaN]","[-0.0072;-0.0019;-0.0022;-0.0263;NaN]"]]}}
%---
%[output:395feee0]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"36×19 table","name":"BRW_result_LMSlog","rows":36,"type":"table","value":[["'BRW'","2025","10","'month'","'BaG1_A11'","'abs'","'log'","'LMS'","0.8948","0","-0.0199","0.0246","-0.0645","-0.7898","0.9754","-2.5549","-0.1808","-0.1798","-0.1818"],["'BRW'","2025","20","'month'","'BaG1_A11'","'abs'","'log'","'LMS'","5.0875","95","-0.0429","-0.0261","-0.0598","-1.8548","-1.1256","-2.5839","-0.5764","-0.5760","-0.5767"],["'BRW'","2025","10","'month'","'BsG0_S2S20'","'neph'","'log'","'LMS'","0.3227","0","0.0050","0.0360","-0.0260","0.2938","2.1148","-1.5272","0.0513","0.0522","0.0504"],["'BRW'","2025","20","'month'","'BsG0_S2S20'","'neph'","'log'","'LMS'","0.8371","0","-0.0067","0.0094","-0.0228","-0.3724","0.5173","-1.2620","-0.1260","-0.1253","-0.1268"],["'BRW'","2025","30","'month'","'BsG0_S2S20'","'neph'","'log'","'LMS'","1.0222","0","-0.0042","0.0040","-0.0124","-0.2286","0.2186","-0.6757","-0.1187","-0.1181","-0.1193"],["'BRW'","2025","40","'month'","'BsG0_S2S20'","'neph'","'log'","'LMS'","2.2116","95","-0.0059","-5.6283e-04","-0.0112","-0.3200","-0.0306","-0.6095","-0.2097","-0.2092","-0.2101"],["'BRW'","2025","10","'month'","'BsR0_S3S30'","'neph'","'log'","'LMS'","0.5618","0","0.0082","0.0372","-0.0209","0.5592","2.5502","-1.4318","0.0850","0.0859","0.0842"],["'BRW'","2025","20","'month'","'BsR0_S3S30'","'neph'","'log'","'LMS'","0.9253","0","-0.0069","0.0080","-0.0219","-0.4352","0.5055","-1.3759","-0.1291","-0.1284","-0.1298"],["'BRW'","2025","30","'month'","'BsR0_S3S30'","'neph'","'log'","'LMS'","0.5633","0","-0.0022","0.0056","-0.0101","-0.1374","0.3505","-0.6254","-0.0642","-0.0636","-0.0648"],["'BRW'","2025","40","'month'","'BsR0_S3S30'","'neph'","'log'","'LMS'","0.7666","0","-0.0020","0.0032","-0.0071","-0.1228","0.1976","-0.4433","-0.0755","-0.0750","-0.0760"],["'BRW'","2025","10","'month'","'BsB1_S11'","'neph'","'log'","'LMS'","1.0863","0","0.0155","0.0440","-0.0130","1.1514","3.2711","-0.9683","0.1674","0.1683","0.1665"],["'BRW'","2025","20","'month'","'BsB1_S11'","'neph'","'log'","'LMS'","1.4642","0","-0.0109","0.0040","-0.0259","-0.7598","0.2781","-1.7977","-0.1967","-0.1960","-0.1973"],["'BRW'","2025","10","'month'","'BsG1_S11'","'neph'","'log'","'LMS'","0.8172","0","0.0116","0.0400","-0.0168","1.1056","3.8115","-1.6004","0.1230","0.1238","0.1221"],["'BRW'","2025","20","'month'","'BsG1_S11'","'neph'","'log'","'LMS'","1.6766","90","-0.0127","0.0025","-0.0279","-1.1274","0.2175","-2.4723","-0.2243","-0.2237","-0.2250"]]}}
%---
%[output:8b15310d]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"36×19 table","name":"BRW_resultLMSlin","rows":36,"type":"table","value":[["'BRW'","2025","10","'month'","'BaG1_A11'","'abs'","'lin'","'LMS'","0.6296","0","-0.0033","0.0071","-0.0136","-4.0783","8.8766","-17.0331","-1.0326","-1.0323","-1.0329"],["'BRW'","2025","20","'month'","'BaG1_A11'","'abs'","'lin'","'LMS'","3.8136","95","-0.0069","-0.0033","-0.0105","-7.2761","-3.4602","-11.0920","-1.1382","-1.1380","-1.1384"],["'BRW'","2025","10","'month'","'BsG0_S2S20'","'neph'","'lin'","'LMS'","0.3199","0","-0.0300","0.1578","-0.2179","-0.5467","2.8716","-3.9650","-1.3004","-1.2953","-1.3056"],["'BRW'","2025","20","'month'","'BsG0_S2S20'","'neph'","'lin'","'LMS'","2.6771","95","-0.1338","-0.0338","-0.2337","-2.1912","-0.5542","-3.8283","-3.6755","-3.6700","-3.6810"],["'BRW'","2025","30","'month'","'BsG0_S2S20'","'neph'","'lin'","'LMS'","2.6046","95","-0.0653","-0.0151","-0.1154","-1.0342","-0.2401","-1.8283","-2.9577","-2.9536","-2.9618"],["'BRW'","2025","40","'month'","'BsG0_S2S20'","'neph'","'lin'","'LMS'","3.8894","95","-0.0706","-0.0343","-0.1069","-1.1228","-0.5454","-1.7002","-3.8228","-3.8188","-3.8268"],["'BRW'","2025","10","'month'","'BsR0_S3S30'","'neph'","'lin'","'LMS'","0.1089","0","0.0079","0.1532","-0.1374","0.1839","3.5613","-3.1936","-0.9209","-0.9169","-0.9249"],["'BRW'","2025","20","'month'","'BsR0_S3S30'","'neph'","'lin'","'LMS'","2.3821","95","-0.1047","-0.0168","-0.1927","-2.1398","-0.3432","-3.9363","-3.0948","-3.0900","-3.0997"],["'BRW'","2025","30","'month'","'BsR0_S3S30'","'neph'","'lin'","'LMS'","1.6636","0","-0.0370","0.0075","-0.0814","-0.7400","0.1497","-1.6297","-2.1095","-2.1058","-2.1131"],["'BRW'","2025","40","'month'","'BsR0_S3S30'","'neph'","'lin'","'LMS'","1.7158","90","-0.0257","0.0043","-0.0557","-0.5197","0.0861","-1.1255","-2.0280","-2.0247","-2.0313"],["'BRW'","2025","10","'month'","'BsB1_S11'","'neph'","'lin'","'LMS'","0.2615","0","0.0192","0.1658","-0.1275","0.4983","4.3095","-3.3129","-0.8083","-0.8043","-0.8123"],["'BRW'","2025","20","'month'","'BsB1_S11'","'neph'","'lin'","'LMS'","3.2046","95","-0.1023","-0.0384","-0.1661","-2.4208","-0.9100","-3.9316","-3.0456","-3.0421","-3.0491"],["'BRW'","2025","10","'month'","'BsG1_S11'","'neph'","'lin'","'LMS'","0.1021","0","0.0056","0.1147","-0.1035","0.1950","4.0170","-3.6269","-0.9443","-0.9413","-0.9473"],["'BRW'","2025","20","'month'","'BsG1_S11'","'neph'","'lin'","'LMS'","3.4309","95","-0.0857","-0.0357","-0.1357","-2.7778","-1.1585","-4.3971","-2.7139","-2.7112","-2.7167"]]}}
%---
%[output:6dac349e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXUAAADhCAYAAAA3fSUpAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQuYFcWV\/9GgYAQFjRuRUQefkGgIMajxESCRGN1PEpANIKuA4Cu6iw9QEQwgLwXxE6KrCGRAoxi\/D91IEjUaeUSNRNngmo3544OJoKOLrg9wdzCj\/PnVcO7U7VvdXd1dfW\/3vae+b76Z6a6qPvWr6l+dPnXq1B47d+7cSZIEAUFAEBAEqgKBPYTUq6IfpRGCgCAgCCgEhNRlIAgCgoAgUEUICKlXUWdKUwQBQUAQEFKXMSAICAKCQBUhIKReRZ1Z7U255ZZbaOHChb7N7NSpE33961+nsWPH0qmnnkp77rlnIe+rr75Ko0ePprffftu3fPfu3alfv340cuRIOuyww1S+n\/3sZzRjxgz196xZs2jYsGGF8vAxmDt3Lt19993q2oEHHkgNDQ103HHHFfJs3bpVyfPyyy\/T8ccfT4sXL6aDDjqo2rtK2ldBBITUKwi+PDoaAmGkrtf2L\/\/yL4Sfdu3aqcs2pM7lQe533XUXHXPMMfTCCy\/QiBEjqKWlhX70ox\/RTTfdRHvvvbfK+vHHH9MVV1xBzzzzTOHRs2fPpqFDhxb+\/\/Of\/6wmk\/fff7+kfLTWS25BwA4BIXU7nCRXBhCIQurQ2u+9917q1atXZFJHgSFDhtD06dPpf\/7nf2jMmDH017\/+lU444QS65557qEuXLqpOnbAZnvPPP59uvPHGwmTy61\/\/Wk0uSJMnT6YLL7wwA0iKCNWMgJB6NfdulbVNJ\/UHHniATj755KIWvvPOOzRz5kwCkSLpWrOuqZ9zzjl088030z777FMoD018\/fr1dMMNN9CmTZvoiCOOUKaUL33pS3T99dfTypUrqWvXruoaNHikX\/ziFzRx4kRF4F\/4whdox44dyvwD4kc53TyDPPfffz\/16dOnynpFmpM1BITUs9YjIo8vAmGkrhMt\/p43bx4NGjSoRFM3kToy\/d\/\/\/V+BwA855BBF4EcffTTdcccddNttt6l6YBP\/zne+o8wx0OTvu+8+6tGjB5144onqy2DfffdVv3v37k3bt2+nq6++mp566imVZ8mSJWpikCQIpImAkHqa6ErdThEIInVoxVgEhQYOTV0nZQgRpql\/+umntHbtWpo0aRJhcVPXuGEzv+CCC1RbLr\/8crrmmmvovffeo4svvpg2bNigbOWDBw9WZppPPvmk8IWwefNmZU9\/4403yG8icQqQVCYIyOYjGQN5QsDWpr7XXnspYoWWvscee5Ro6jZtnjBhAl166aWqvE7O\/\/iP\/6gmjo0bNyqiZxKH9s5eLryg+tJLLxUWWaGxY1FVkiCQNgKiqaeNsNTvDAFbUj\/jjDOUbby+vr7wbFvvF9i+sdh51VVXUceOHVV5kxnlscceU66O7MYIO\/tPfvITeuihhwqmFs6DOmCSOe2005xhIRUJAn4ICKnL2MgNArakjgbB++XOO+8sEKkNqUPzhrbdrVu3goaPuvQFT9jM4SuPRVIsnoKoYXPfb7\/9CgunyAP7+cMPP6xInhddDz300NxgLYLmFwEh9fz2Xc1JHrZQCo0a9nRsEtq2bRt973vfo1tvvVVp3CabOsw0sKNDw4Y9Hv\/jb\/iZs387g6y7JoL4n376aWUrh4kGphqYaXQXx0suuYRefPFF5VGDLwcstLLmX3MdJw0uKwJC6mWFWx6WBIEwUkfdsHHDBREkDB91aMwHHHBA4EIpFkKxAIqJAGQOmzkWPvUEGzoWPZuamoquszcMLuqbkdq3b0+fffaZ8pIRe3qSXpeyUREQUo+KmOSvGAJhpA4PlieffFLZ00HQtqQO8wq2+mPLPxJ2lMLX\/Mgjjyy01bR71GtW8YYN4MI68VcMPHlwzSAgpF4zXZ3\/hkaxqaO1cDG87rrrlPYd5tKISQCuivApRxo+fLgyxUDjRtL90hlJ9oSBDZ0TzDIwz3CCHX3p0qVqopAkCJQDASH1cqAsz3CCQBRSx+afBQsWqEVPpDBSRx7YxDERwE8dE8Htt99OZ599dkH2Rx55RBE\/J5NZRXd\/RD6xpzvpeqkkAgJC6hHAkqyVRcCG1LFzE0SMxU49GqINqcN8AtfDadOmqYb27NlTmWF4YtAXQkH6yOsNVaC7P6IOfSG1sujJ02sFASH1WulpaacgIAjUBAJC6jXRzdJIQUAQqBUEhNRrpaelnYKAIFATCAip10Q3SyMFAUGgVhAQUq+VnpZ2CgKCQE0gkFtS59jX8Cf2eiDURM9JIwUBQUAQMCCQS1LnI8YQ2tR0Ao70tCAgCAgCtYpA7kgdGjo2leAMyfHjx6sdg6Kp1+rwlXYLAoKAF4HckTo3gLV1L6lv2bKF8FNXV6d+JAkCgoAgUEsIVBWpg8wRBnXdunU0btw49SNJEBAEBIFaQqCqSP3555+n8847T0XbO+mkk0RTr6WRLG0VBAQBhUBVkrosnsroFgQEgVpFQEi9Vnte2i0ICAJViUBuSd3UG2x+EU29KseqNEoQEAQsEMgcqes+6Oecc446WmyfffYpaooegvWQQw6hhoYGOvroo0lI3aLHJYsgIAhUNQKZI3UQdn19PQ0cOFCdNendMRq0k1RIvarHqjROEBAELBDIFKl7fc9\/8YtfUGNjo9pgxAl5cPoMzqGEdq4nIXWLHpcsgoAgUNUIZI7UdcIGqT\/33HNFJhgmbu6VSy65pED6fA\/+6eeee664NFb10JXGCQKCgAmB3JG63gg2xZxyyinq+DKd8GXzkQx4QUAQqEUEMkfqfAI84rmYzC\/eToINHgkmGtl8VItDWNosCAgCOgKZInUIFrZQCuJes2aNInGvDV5s6jK4BQFBoNYRyByp6y6NbC\/HNZzwPmXKFDrggAMU8S9cuFD1ncmmLn7qtT6spf2CQO0ikDlST9IVoqknQU\/KCgKCQDUgkEtS1xdEZ8+erRZJkYTUq2FIShsEAUEgCQK5I3XdTx0NnzVrFs2bN0+ZZWqB1BFeeMWKFVXtsiltTPJKZ6es9GNl+iJ3pA7ihk19yZIlKnyAvutU91NH6N1qTG+99ZaKGQ+XTWljfntY+jG\/fadLzv2YpXW8XJL68uXL1YYkJJA6+6nrh2RUx5CRVggCgkDWEYByhTMcsnLSWlWROjqfj7PL+kAQ+QQBQaA6EMja0Zm5JHU\/80t1DBFphSAgCAgC8RHIHakHLZTGh0FKCgKCgCBQHQjkjtQBu+7SmKUFiuoYEtIKQUAQyDMCuSR1E+B+vut57RwOVrZy5Urq1auX8vaB26ae9Dbrh4Xkpc02beS2eENC5KWNkNO2nfpO6bwpKzZt1HeL53G8+o25oDMeKjFOq4LUq9Ekowcz43g4vMkKA8UbVx75H3zwQSP5V2Jg2TwzrI16HUx4eSM7tMGmnXqeV199lZYtW0aTJk0qOfXLBtdK5LFpox58zxRWuxJyJ32mPlFlZWxWBakH+a4n7bRKlDeFFGY3Tu\/RfiwfiEDfiFUJuaM8M0ob0b+\/+tWvaNu2bSUnYUV5ZiXy2rQTeWbOnEkjR44sOfilEjJHfaZNG72Tm00E1qhylDs\/2r1gwQIaMmQIjR8\/XgUZRHTZSqeqIXU\/3\/VKAxzn+d7POX3S8ppguP68vSS2beRgbtdeey3NmTMnt6TOxzKa+pJJHX0JbY9\/Z4EgbMavbV8ysU+cOJH08B42z8hynqyZBoXUMzhaorwkEB9EEabJZ62Ztm3EJ3vfvn3VuoLpzNqstcsrj007mRSGDRtWOOyF3Xb9JvEstdumjV5tvlrML+gHIfUURmMtm1\/ypqFz99t8suv2Sn3YZMV2aTOUbdrpJUWY0vjLxHsOr80zy53Hpo3eNaC8tTEIUyH1FEZcLS6U8qcsfusLqCnAm1qVNotr3kmAzRipCZVCxTbt1PPYmNtSEDNRlWFtNBF\/nr5GhNQTDY94havNd113ETvnnHMKh2\/j5UH6xje+QaNHj6a33367AJif62M8RNMvFdZGfbLKmttYFHRs2qnnyaO7n00bq9WlUTT1KG+D5BUEBAFBQBCIhEBVLJRGarFkFgQEAUGgihEQUq\/izpWmCQKCQO0hIKRee30uLRYEBIEqRkBIvYo7V5omCAgCtYeAkHrt9bm0WBAQBKoYASH1Ku5caZogIAjUHgJC6rXX59JiQUAQqGIEhNSruHOlaYKAIFB7CAip116fS4sFAUGgihEQUq\/izpWmCQKCQO0hIKRee30uLRYEBIEqRkBIvYo7V5qWAwSmTSOaMiUHgoqIeUFASD0vPSVy5hcBP+LG9alTW3+E2PPbvxmTXEg9Yx1SdnFEU0wXcj\/ibmwk6t697dmrVhH161csS1DfSL+l2285rl1I3a\/zauGlsdAUp02bRlN8tMigezl+J9yJ7kfcuN6\/PxF+c6qvJ9q0qe3\/oL6x6Dd3jZCa8oaAU1LfuXOnOrRh06ZN1LlzZ+rRowe1a9eucpgEEHMgIeXwpYlMvhaaIuqcuss0gB8vsQfdq1yHV+bJRuz9iBvkvXq1WVDcw09Dg78Wb9FvqFwm3MqMhSw8NTap46STJUuWKAKfMWMG7bPPPvTMM8\/Q5ZdfTtu2bVNtGzRokBpcHTt2tG6rfjqKfuKPt4LQk3ACiDmQkCxfGusGOczo96JGJl8LTbGxsZG6a+aBVatWUb\/d5oGgew6bm4uqfLGHJu5H3mEtA7F7tXgQPa6Hafi7Cd1vMhbCDwM\/\/\/djk\/pvfvMbuvLKK+nUU0+l+fPn05577knjx4+nP\/zhDzRu3Dh1wvY999xDkydPpgsuuMAaKZxbWF9fTwMHDvQ9PV4n\/qJDiPUB730pIEF9PTU2NPiSlXqRLF8aP5OEdUMjZvQjj1jk60c4IO1Vqwh19u\/fX\/1usw7Uqwk86J63SdWuLfpin4TQI44LZG+sr6fR9fWEiTdswpUvrBgA56xILFLfsWMH3XDDDdTU1ER33nkndenShV555RUaNWqUIoObbrqJYIq58cYbCXlvvvlmpcmHJe9Zf\/phtlwWGvqCBQtoyJAhahK57rrr6OSTT269Dc1SJ3PPA0FR2tKUuosJhFNDYyN5lqrULf2lSfOl8CNBvxc1FvmOGtXqbeGXpk6l\/qtX02qDlsmaut89kAqnNHEKG0dp3Pf2jS\/2wEBfAE1DGEOdS4loTX09LfWMfx7f+N0QpNBodVbTZBzZLLkbhzAMwu6XqduNj4lF6ky+IFOQKtIjjzxC11xzDc2ePbtwuj20bhwIDTPNAQccENpO1Is6MGEcffTRBFJ\/7rnnjJOC8bDXEFLvT0Q+1kwCzQV5Cy+uq6N+q1YpuTg9+eSTdNppp5W065ZbOtB11zUb2wtMGDM9A67jk3nSpEnq64bTli1baMCAASVaMybRs846Kzr5ElEb9Zq7JAinoE5k2SGzDU6hA8JhhpaWFjWOrr\/++sjrPKa+8cXeAl+HzSqqavQuJQXkbpuY8Ovq6ghjGclvHNrW6Tof+m379u3KhOu3Phf1nQprZxgG3vuQq6Jrhx7QnZA6gJ8+fTo9\/PDDdO+991Lv3r0LAyRLpL7LK1iRtymFkTrK\/HNdHd2\/ZUuhOF6GtWvXFlU3f37nXeaoLrtMUB\/QuHEfeu7NV6YqmKfwoxP3t7\/97cL\/uklp+PDhtG7duhKR8WyQZ5wU1lbTF43tc+bOnavaqMtmwsm2Plf5IJMJe71+vq9fQzu8fYP+QF5Twpde2KTpqk3eeuJOxqgH\/YY0YcIE4zhMS2auF++N933BPXzpz5kzh6699lpq3759iRjz5\/9yd7\/+wPqdMvUpf+0H3QOXod\/19xFje9GiRfT9738\/bYis649F6p988onSeGBSAZl\/8MEHdMkll9C+++5bMMfgGhZNu3btSrNmzTJ2iFdKG\/MLygBckB3SRRddpDR7ZXax+OwNGvh4GU3mF5bTRHbQ1FnLgQg9e3YoNOuxx5oLrsf4XO\/Zs6d27zG18Ijr0MS9JPjqq68qzR1aQRoprK3Q+KD5uUoY\/NAMGSvec+Oqfq7HVK8f9vqzp05dqrCeNGlE4UuprW\/Yg6uRbCbThl0mu11GrrKnNPoMhMVmt7gNYouQZuksqgp9hq\/bSZNadmHf4rlX2i\/IAOsgyq1b1\/a+1dW10KJFKL9a8YL3nUJbmDOi3oPJMehdhGkL5ucspFikDsGhkcPrZfDgwfTxxx\/T448\/TldddZUi9z\/+8Y9qkRTke\/vtt9PZZ59t3dawhVI20WDCgM0es\/fSqVOp86BBgfb0IGLme7Cua57CRplNLw7Id+TIKb5rrEG2b6xB+NmoTdetgQzJGNRWTF6Y\/DQvaiePBaljsTXMYzTunhtvvbzYqy\/4oiG6nXnKlFVKHn0Jgb0KYTddvbqvZ6oHKpju\/Ax5rVCFTZpOADVUskcKFQOvVseAekXwJnI2ETcw9cMWa\/JB99EMv7Jr1gQvC7X2TxRDlBvQwAOHH354xck9NqnDzjVv3jxatmyZQuSMM86gmTNn0he\/+EWlxT\/xxBP04x\/\/mC677DIrLZ1h1T1bMEHA\/oxrvDCxceNGNWPeeuutdPXVVytb2+0ffUQH\/frX1j0TpNFgroWmZUpBZFdfD8+DUj2\/1QvQTNx4WbyEY90IBxlNbbWjrPgPb2hYRaNHt+Hk3UgZd8+NyRN19OjuIfiGGaKC2hlMHEHjKD564SXhCOB6Mm59ajFWPPH5kS\/G\/dIAXoVSG3Q\/vKVhOZIYo8Lq9r+PSU93GIhfU\/ySsUkdj4SHC3zSP\/vsM9p\/\/\/2VWyPs688++ywdddRRdMghh9Aee7jVHaD9L1++XC16IWECOeWUU2jo9dcrTd3r3eJHzn6QtfnClOYIflnil4zffclK+kmcDimwrMVP9Wp9Jk9ULml7D\/lbXb3DWhLUZzbYxh8RNrXHyRPW4jh1tpZJilX8J8crmR4SQfLo+zniyZ28VCJST\/746DX4kvqXv6x8zN1OIdHlkxKCgCBQuwjAvKi7SVcCCWekzlo7NHW493Tq1Mm5lg6AQOowv8BNEgu10NSxaKpWr6dNo\/4BPthBVlDVEY2NRn0Ec374vA9NJkybCbbDlm8A1FM\/T4vs2uhCQhucXDynoN8bes+lDP59GrTo7rKFqKs8\/VfOFrlAqPzvGxZLsWhayRSJ1D\/\/\/HN66aWXCuaV008\/XXm8bN26Vdm8scuUE0h24sSJdPzxx0dqX1iYAF4ohT39tttuow8\/\/LDgBw\/Ch7Z+crPBR7xfP7X4F7pxxuvrXl9P3RsbQ0gdJIGODBv0rv0T8FzTdMNymAY17jVQP+pecL1L247eNgBscYo0ZCwye+2rNkviFtWqLGYrdmzXRhibYTbyrtyGbBoLcte1bUl4vp3hWTKVo7zf7VAMYX7JjaaOhVF4uzz00EOFbjvxxBOVbfuOO+5QPuoHHnhgweXr\/fffV9vx4QVz5JFHWnd1mPcLKvrtb39Ll156qaoTG3UuvPBC9TdI\/bzzzqO1mzdTXYvmGgUtHCtyIOjuxYtnJR2BlwnbvJF2rwZ1Hz3awYJbHJ8S1vz9iBs+Fl5SQRlcx2+\/e61kxEtfUQiBByzbq\/UBHG7DTrIwaT2EDBlNzqiuCMpMHKOwg9MbvAv\/w4sEv02uHbjHbnG8Yoz\/WfPzcxeZMkX5eowe7dIJ1Qtj7GkqScclKAuFZvd7nKAW26K8Yzep+6ft84LyWWvqK1euVM79ffr0US47zc3Nars+zC3QlnHtiiuuUJ4u2DCADTTwYYdGjes2ycZPPShMAJP6uF69aNy\/\/7t6ZEtdHbXAP3V3MCpo6tgNiMQbB7wd0W7GDGo3cyY1Y7fs1KnKv9yfsGy0vri6MJNgEHHrg9erCZvvwZ+XAybAUQ2kzilsQxMwgw89Jnh4O8E7Ca5cSPC3R1A3c6o0KXi\/klw4HZqJQx9XwOn+mTPp1BEjaPHixUXQ8DhrGTGCWjz3kLHdgAHUsnunpxdT0z2MbfSFd3MMdvqCdHjcc12Qc8SIERZ7ISr1hWXDGqY8cd+3uM9rdZOFPT0LyYrUOdbLhg0blKmDtbP169crMj\/00EPVdWw04vTee+\/RxRdfTN26dVPaPMw0YSlpmAAmdTxndvv2NPSvf6UPxo2jD7Xdm7hns7uw8\/z5hXLYURhMdmEaaJAuHGRC4b2JQcSNFnH9JjlK74HU167dbOwOrE+89dZbgV3Fu2hNOzBRsBQv\/eshbBSkcd\/0leRikjH3q3fHMDCFx5Y34b3qOnw4NS1fHsntNwwhHt\/nnntuYacoyujvBwgdOzVhJvXbtdz2nLDxHSaR6\/voz6D1qyjfnm5kyx2pswaNgaATNGzpY8eOVSTvDdrFoXHffPPNssV+4UGLLc8nnHAC1S9b1qptG5JfvAhTXmg4rKmbzQ4oFTRL+\/nMMrEEaeIsURBxM7H7Ra9B2bZ7+Pp\/5RVzbBoXQxzaOrR2JIyNUaMaaOrUsPUGF0821eGntUEe\/MQ9H3Qp1ddjR3CbEy1v0hk2bJhVY\/COvPvuu3TwwQdThw5tOyOtCodkwph97LHHSnLxFxY0ej0Gkd5nbYVAnMAnGzslixuDMW3aGMbrSNFNUaZ9I9yn2I+Dd99vI5tigDxp6kzqhx12WBF5+11HA21InfPAtAOfdmgYGGy4\/sYbb9DXvvY1Ze7B5yMnLvPoo4+qS4i5DjOQrokUheN18YZodZiiD4bEEdtdOq4JRW9AMTkHNc1vSzaXSXv8eXGyw8hxZ6lJBWaXaSUB0bBZbPVqfvFBXPqkg4kAm+rw2+zJwJtvYNXDztU4G05gwkSkU3zhuib1ICT9IgzCfIO28GTc2Giz+O++z1prDCJtnmT4neKv3dYvwlZLgr52Vl+yC7aNnIvvwXNFxwD\/6+ZZL0be+2mhEaVeK\/NLWqRuEhQHaxxxxBGKyBFTHYugvBCK\/ByOF7EdRo4cqSYA\/D106NDC52WapI7BwFoQm6HsCCuJCSVKl7bm9Z6OFr2GZCW8ONlhFP2ZfFiQ6TyK3eHhlceT90UFifN6eNuCcjE5tH3ig2AwOYBMWjV7F2dFAyNogBjHlfaYYOT1yXj16imxz\/mI3pOmEv6k3bqvzPxO8djntowaNZUaGoq\/yGzumU78UtNNwIlgbtqdrJZMkTomDwxwaDCvv\/66ss1ikRKkzdoFh\/DliQb2fBzUoZM6bJonnXRSMmR8SsPejEh2+jOGDz85KIw7tWvXGk2xpWUJ6GDX\/2MVObS01HmeUqyJc7m4Dfn97+NFcYz7PL2cF6fzz69TGHnb3IaNF4vW2o46qrUNr71Weh\/37ruv9b63D1Dvz3++hQ4\/vLUerPlgoRLR9PiAE0S4uPtuxOI3k0OPHs8rmZubkae1b1Dv2LFbaOTI5CiZxlLyWpPVAJlgUgWhNTefRNdfv\/usgmTVlpQGjqXj35vN3C8oO2PGFjWxNjc\/od6p1vdqTEnfnH\/+ErrvvjFG6YPuBcVL1zECz8AsjZ+spMyRum08dQDoXVjFYiYI1xSqNk3Am5oeCB2ghx7aGloX0StxqAjS5s1rA8th8HK5NOUvZ90giqam1kVDtO+ggyZQhw6toYW9eOB+167nqXxB91A2qF5un449X0O5Dz4YV0IOXbrMp44dVxSgQR9DllpLYWM0Lh7c963Yt00cuA7s8dumX\/D8pqZjqWvX\/xdXlMTlvAvjiStMWEEkUsfGo6ipV69eqSyU+p1RCmIPc8uL2oa08k+cGKzh47nLlz+f1uMrVu+SJXW7tOY6pWmdeWZb+6CNQ6NHwif09OlbijR1v3vcEL96bRs6ceITNHv2mbbZayIfj1FvGB3vMaomMEx59LUefWzjObNnm8c6ggOeeWZ2+yXXmrprUjctlGJ1no+oMx1nh8Hj1dBr4u2qskb6hdeNG6WR4QkK21tlEJa1OYY9eYS9TkEx03gx2bOXr3DGQFkbUEMPs9LUy4lH2I5SPQyvzRF55ZRdnuUGgbjx1N08XWrxQ6A0Xn0rqZu0eK4DXlZh8fMFcbcIZI7Uw+Kp4\/SShQsXFqGgn4vqFh6pTRAQBHQETBOuSYvfvYG7UFS+oMo3jjJH6uVruv+TvGYffD3wROI9WNt0vWhnq3YQd9y2JZVHlx\/7AeBbqx8MHUUuyIJAbUj6eolfm6NejyIL53UlkyucosqDdpjWiFyOI1cy+WEURxuPIpNuqgVefIAO\/naFkyt5XI2jOO8Cygipe5DjDuFB4xfqF8VMIYCPOeYYYg8e5MH5rDghKq6pKKk8IN6i8MRxRwqRivmitweyvf322+pQYARWU2fFam3G3yYs\/K7HwciVTCVhnGPiFFUetFn\/OuU9Fvq6UdJx5EomP+cEhiqKNh5VJpxtC79+PgltzJgxhJ272Lns4n1zJQ\/21rh632IOQSF1HTjM1DhjcA0OQdylPWEAmbRkPVIhb7XmtQCU9433HrGXXMijD\/q42rmf2HxgCc6pxVm03hj3fhOf33U+0T0iTEXZ48qkT8YucQqTB5MudkQPGTKExo8fX3ASCDw3IAlAuzVbxKLx6zc\/mdJ0UAjDyTs20njfdFjjypPm+2bb7aKpG5DiU8NB2KaXCyaMvn37lpA3Xzcetzd0qG2flORLKg\/CEXPSP1tjC7S7oP5imdqMCS7KdWwgS5qSyJQGTmHycJu9EUp9T\/gqA0ZBMqWBEU\/yUJb8xow+NvTJBSG+Xb9vSeVJCyPbd0NIPYTUuYNhOwdpIw4ETnUC4bNpRL8Osnc9yHRSjyqPHrSJP5\/Vma4JyUH\/gvEjoHKTehKZdDxc4WQjT7lJPYlM+qviCiPUGUUmrwkojckviTxpYWRL6MgnpG5B6noW1ry8pJjm56CX1KPIY5IT5XWyjzJgdC2G64667pCG+cXbL1FlMn3eJ8HJVh5+rklTd2XG4\/5NKpN3nASNS9sxFUUmk\/nHtZkqqTxpYGSLJecTUg8hdV0TwE5VLAoiDrX+2adfx+lPLhZuvMTNBBNVHsiJNQJ9gYk3d0UdLEw5DS9RAAAgAElEQVTG+BrRSdBvUQ\/5014odSXTxo0bneEEYrDFyBvLiPvG5UJpVIz8ZMLYq+RY8sZ\/Qrtc4hSl3\/BskzyuMYrzjgqpW2jquouSHgHS77ruYuUiYqSf+QWi28ijy5nEpq63i2FD6GPE0sduY7Yl6jL5YeEKI5cyucApjjxMTvDo0CfcSmLkJ5MLjFB3VJz05\/LYY\/diFzi5lMcVRnEIHWWE1OMiJ+UEAUFAEMggAkLqGewUEUkQEAQEgbgICKnHRU7KCQKCgCCQQQSE1DPYKSKSICAICAJxERBSj4uclBMEBAFBIIMICKlnsFNEJEFAEBAE4iIgpG6J3Oeff05PPfWU8tPFwdM4F3WPPfawLC3ZBAFBQBAoDwJC6pY4v\/HGG7Rhwwb64Q9\/SHfffTfBPxuHXksSBAQBQSBLCNQsqes7MxF2FUmPp2zaNNTS0kLr16+nP\/7xjzR69Gjq2LFjlvpSZBEEBAFBoDY3HzF5825IkLoeTxlbxhGU66abbiL8vWPHDurRo4eKiY5t9ygPE8yRRx4pQ0gQEAQEgUwhUJWauh4PAid9I2g9RyaEhv63v\/1NdcJzzz2ntriD1EHU\/D8iwSFmCcr9wz\/8A+3cuZO2bt1Kf\/\/73+m4446jZ555hv73f\/+Xvve972WqM0UYQUAQEARSI3UQIU7F2bRpE3Xu3Flpuu3atSsb4tC8YSKBDCZTik7iTOrek1X0OBzbt29XtvQuXbqoSeHKK6+MfZpR2UCQBwkCgkDNIeCE1KHZ4tQbEPiMGTOU5gtt9vLLL6dt27YpUAcNGqSimpXTDs3HrbE2rvduVFLnsp9++intvffeNTdQpMGCgCCQDwSckPpvfvMbpbnCzW\/+\/Pm05557qqO5\/vCHP9C4ceNUeMx77rlHnWN5wQUXlAUZkPaDDz5IF110Eb388ssl8cNNpO41v+DMTZdHm5Wl4fIQQUAQqGkEEpM6FhFBfk1NTXTnnXcq88Qrr7xCo0aNov79+6vFRphibrzxRrXgaNKaXfcAJpGVK1fSyJEjVdUgcCT9wAgvqZsWSsshq+u2S33ZQgAmPSQ+11aXLuheklYkqTdJ2SQyS1l3CCQmdT6xBYcm8Gk6jzzyiFpo5HjHEBemECxSwkwT59R4d01urclL6nxt4sSJ6ti6hoYG0dJdg15D9a1evVqZG\/GbE4gd4wrJ7x6OSwxLfsQb9EyuN0nZMLnkfjYQcE7q8OWePn06Pfzww3TvvfdS7969VUuzRurZgF+kqEYEQNhTp06N1TSQPr5yTSmItLHTOeiZuIfypkmm37Rp1KhNPt5n12OiWbUqVnukUPkRSEzqn3zyiXL9w+IoyPyDDz4gnK6z7777FswxuIZF065du9KsWbOoffv2iVqqnywiWnUiKKWwYwRAmjA7JkmrVq1SB5zrKclEESbLJpiHAjLBgNRokCmsXrlfGQQSkzrEhkYOr5fBgwfTxx9\/TI8\/\/jhdddVVityx+xKLpDC93H777XT22Wcnaqn3NPFElUlhQcAxAiB0XRuOUz009SlTpqiiMNkknSigYweRdtA9yABSnzZqVMF0FKdNUqZ8CDghdfhwz5s3j5YtW6YkP+OMM2jmzJn0xS9+UWnxTzzxBP34xz+myy67LLGWbjpRXIcLh0Dj5+CDD1Y\/kgSBciGAcZdFb6kwTTwMH5B6dyKCQiWpFAHsvynnHpywPnBC6ngIPFzgk\/7ZZ5\/R\/vvvr9waYV9\/9tln6aijjlKLjy6iGnoPiNUPUsZLNWHCBFq3bp3yfEnDfRJtwm5STFhZ6UiRKWyYkxqLafQbxhwSdi7j729\/+9vhwpQ5hytSX7t2rWpnOVNa\/ZakDV6ZsLkSXn9ZSc5IvRINYlOMHgIAJ9rPnTuXTjjhhFQ0dTzz3XffVXV36NChEs0ueabIFN4NrjGCSQSLj1AgODGxh0tT3hwwr+An7lJnJTV11\/3mAnmvTFWrqQMs7LbERh9sOkKcFOwePfHEE+mrX\/2qU61W19ahGfXs2VO5U\/J1U1gAF52JOpqbm5VPPhZ9s0LqIlN477rEKM1Fy\/CWxMsBUofGHicxqZt87XENC7tpJZf95krGLMqkt82Jpg7TCz7NoLlwsCz9ITC\/QHvu1atXYlx\/+9vfKnv94sWL6aOPPqIxY8aouhFcS0hdJhq\/AZbkRdR9u5MuWiZ+AWJW4ILUTY8GqSM8SFopSb\/VkkzOSf1Pf\/qT8q1FTBTYsk8\/\/XRlQ4eNETZ1LKDut99+dNddd9ExxxyTCGsQN0IPIGoiEsLfwpUSm5+E1IXUXZK6yS880eC1LJzW53wSfTrISTNNTR326\/fee08FBSz3lzHMaaY1hCxONE5JHSaXn\/zkJ\/TCCy\/4kvZf\/vIXuvjii1UMcphJkiww6odboCHesLqYUNJezMFAS9IGy3c7UjaRKRyuqBixhh5es7scGFdDhgwpstW7q11qioLASSedRC+++KJaZPcm71hK2wwVRe7E5hcOE9CnTx\/CFnuThwvvMgW5w2c9yUpxGKl\/61vfitJ+ySsIZAoBaKNYr4FJsVu3bpmSrZaEwQI4ghNi\/QyaeVhK2wwV9nynmjp2i0IL\/+Y3v1kSCVF\/kKswASB11IUYMtjFCk19+PDhBfNL2po6z9pZ0tRFpvAhHwWjVg3dZktO+HOj5mBST3OxP6pMtZifTblNTe13k3prYDa\/VFWkjpcFJPtf\/\/VfhbAA3obDfx0Bvg4\/\/PDE5hd8GeDw582bN6vHIOTAz3\/+czWpiE1dbOp+L52tHRSE3r07ttqEeXezT4j5iSZPEeQMM+mkSep4Nta38BvyYf3LT84kRK4fUMP18CSFQHrgAayB2ST93GDk9wYJXLhwoVUAPrghwsECbeYNYpATeEyaNEkpiHpqI\/UHqLkZsgb3N8qawjvYtNF1nsTmFwj03\/\/93yp+OjYdASB9JycWNBEe4PXXX6cFCxao4+GSJHTO2LFjldskkq7RCKkLqWeB1DkaoyniIsIIMLGbvh7eeecdZX5xqakzmZsCfuEahyRI8l7qZUGW2EV+xRVXqMs6eT766KPWpA5C59PIUI8eIgT\/I4gZ1uiCyBn52EQMLuLoq8wV+jnFSUkdziIchdMVlnHqiUzqDNBLL73k+7wDDzywsMMOBzUjQfv5yle+ouKrJ7GpB4UJ4I6Cd8wPfvCDVDYfQePDi\/flL3+5ZHaP0wEuyohM4ShGwahVa4uvqWOhHkQTlkwyYWEOO6FdknqYX71rMjKROq4h7hMTNTRsVsrg6gwzKs5AQII2PmDAABWeGBOOKVS3rvGD7KEw4kAcU96nn36ajj32WOXIAU0d\/QOtGhq7S00dssO9u9IpMqnDhg5vF94eHaUBADMpqQeFCdDvpRUmAAd9YMY\/6KCDEsexiYJdUF6RKRzJKBhhjWbdugdDYxfW1Z1e9B7oXlfYtxGWTDJBWbr66qudkXqbOSlYGpemA5P5hcN56No3a9g4vAZaN2v2Xu2ezxCGIwaTPn6zGYdNKyNGjFBf8F\/4whforbfeou9+97sEDxakcphf8Bz47Kdh0gobS\/r9yKQepfK080qYgFaE87CVOu2xEFZ\/FIxwvu6AAUcFknp9PdGiRU\/Raaedpog9jhutSSYoJhdeeGERqYOYw2zxfu2HJrp06dIweNQ+Ez4pLDSzJwNITCcyr6auE+p\/\/Md\/FJExa9iLFi0i1t6DNHXgwxscdVJHPf\/0T\/9EOFoTGxIRewqHzreuj5SP1HOpqUft8KT59VnfZP\/CIi2ShAmQ0AVBY812oZTrABfsPonOWC1IPelGSpNMpnWhMPNJ0ncsaXmvXd5E6kzeTz75pHoctHOTLRzmVTa7IK9uU0c5vO+YQEDoJps6vqLhOAFXRHjlHX\/88WUjdayhpLkRy7afUtHUOWIjO+h36tQpUYRG\/ZON\/dRhf8NzDj30UJozZ47sKJV4NIFjPs+knkdNffTo0UpT5qR7v0Cb\/v3vf1\/wWkEePb\/u4RLF+wU8AzLHaWs41wEHz+PrAx5yaZtfghbHbcnYVb7YpP75558T7H8cWhf+4TjtCDMlZlp0HCe4L8EexrNmFOH5lCO2yWEmHzRoUMGlEaSOM1GxQFIO7xf2JEjLHSwKNpxXZApHLSpGOLwoSFPHE5Nq6iaZXI9h21g1Lm3q4b2RTg4okfA+AUdgvQIbIrGLPU5qc2mcTc3NrXb51qjy5pSGF1EcuVEmFqnjUAycdPTQQw8VnotojDfffDPdcccd6nxSrwcMbFvYTYpYLbaJV7jxmYXEJpagzUcIvQvvF14gsX2WbT4swCBme5rPsJWF84lM4YjlBSOW06X3S9hpTK69X8J7I\/s52ki9dUcpNhvi8B9MkvoOU1wHocOVOyspFqnD9AFSw0yIAYFGwmaGmfLDDz9U17CSjc8ezJgYoAi6hVV9fYXbFgSv3Xz58uVqAkHSY7\/oh2TY1i35BIEsIuCS1NE+P7t8ljTMLPWDl9S9B2HosV\/AhfjJSopM6iDpG264gTZs2KC26vOq9\/r16xWZwxyC69hAwQlR1rBogVgWIGOYaaIkW1JHnXycXZT6Ja8gkBUEOOaIa1JH+7w7Svv27VtywHVWcKi0HEzqUExhIejRo4evSH7RHCvVhsikzpuP0BCdoGFLx05PkDyu69tu2fXwzTffVIRv2iAQ5OXiJXU\/80ulQJTnCgKuEHBtU3clV63Vk+d+iE3qhx12WBF5M9l7r2Mw2JB60KDRSV3fUYoys2bNUodemyaKWhuI0t78I5AmmZQr9kv+e4HK4nSRFk65I3UAoe8cTeMzNS2wpV5BIAyBNEi9ErFfglwaXQX0ApZhcV8Y77gujXnkl1yQetiL4Pq+N5AQu1XiOaYocd7r+qSj548rZ1J5dPlxIhUHNYojj+43jJgdbE7za3PU65WUyRVOUTHSv2aBKTzLXJJJ2OYl194v5QjopZ905heUi8dSkoBe6AckeNXp7zlbHzheDbtcI4+rcRTnXUAZIXUPcl6\/eL\/47dx53rjuOK4PYYaxmOzCPJRUHg6WxDHn4w4U1op0cxdkwwaTa6+9liZPnlzSZpQxYeF3PY4JDQTiQiZvbP64OEWVB23Wg+QBR5ekXqnYL94oja4DeoFUw4JycR8mCegFr74VK1aUjG19tyv337Bhw2jgwIFFZzzEHUdJysUm9aAojX4C6ZpdEqHTKmvyizdpyezxo29h1rcvu1rIdSEPot0xsXIcaVf48e7ewYMHq\/DK3gnOb+Lzu24bYztI\/rgy6ZOxS5zC5ME7AeLAEXYIXw2t00vq1bCj1HVAL4QZYEWDIy3utddedO+99zoN6AXPl8cff9x4KI8+Dvn9T\/N9s31vhdQNSIV528CEAXcwL3nzdT8\/ettO8eZLKg9\/OqJe\/TMxrjxcTp\/ITG2G7TTKdX5Rk8iVRKY0cAqTh9vM2p6J1MPMJ0nwclHWJvYLH1DhKqCXidSxX+b+++93GtALX05QYE37Yhg73XkDocbTGEdR+ikyqUepPK95dRJlrZJPWEHQHsSYwAzOphGQOV8H2adJ6lHlgZycvFEtk\/SPKR6Pd+CXm9STyKRPKK5wspHHhtTzqKnr5hc93rmrgF4mUocJDW7TLgN6hZG6fnCH90vT1TiK+p4KqYdo6iatGeYXr1aZhvlF14jxt07QXm3ZTx7Tdb+6bAcPt5XrjrrukIb5JalM3hfSO7HbYuPXL0Fn66JMkKYe9dl6\/krEfgk7zs5VQC+0U\/d++eijj5wH9AKpYzHUdCZy0IE9Nu9ukn4NKiukHkLqbBOFFordqlgURFRIfGaxRq5fR8wblwulTIBMxFHlgZwcopSJA5NDXPs1yA5fI3p5v70DkDnthVLGJ6lMGzdudIZTFIx4cTgtUgc+tRL7JY2AXvhChweM1\/EBuJpOZsL76fJ9i0P8QuoWmrruoqS7mvldd+1H72cOgug28uhyJrGpe0+dwvPZnQx2R7Ylms6N9crqCiOXMrnAKY48aWrqPLwl9ks0euR+xAEje++9d8nY1scK18zuyy7GUTRpi3MLqSdBT8oKAo4RSGPzEYsosV\/sOyvNfrCXIl5OIfV4uEkpQSAVBPJMJqkAUqFK89wPQuoVGjTyWEHAhACTSZbi9ddiT6UR175cOAqplwtpeY4gYIGAnAlgAVKZsuCgnblz58Y6VLxMIhofI6ReSfTl2YKAAQE5EyAbwyJrcdJtURFSt0VK8gkCgoAgkAMEhNRz0EkioiAgCAgCtggIqdsiJfkEAUFAEMgBAkLqOegkEVEQEAQEAVsEhNRtkZJ8goAgIAjkAAEh9Rx0kogoCAgCgoAtAkLqtkhJPkFAEBAEcoCAkHoOOklEFAQEAUHAFgEhdVukJJ8gkBUEGhtbJamvz4pEIkeGEBBSz1BniCiCQAEBE3GvXo0g3kT4zQnE3tBA1K+fgCcIKASE1GUgCAJZQsCPuEHaS5f6SwpiHzUqSy0RWSqEgBWpf\/rpp\/TKK6\/Qtm3bjGI2NTXRyy+\/TFdddRV16dLFqil80gsOV+CDFnDGoJ6C8gSdDWglgGQSBLKGALTwqVPjS7VqlWjs8dGrmpKhpL5p0ya6+uqr1YnaQalXr17qHD8+nisMIT5TcuDAgXT99dfT8OHDS45Y88ujk71+yk7YM+W+IJBZBKCh9++fTDxo6tDYxeaeDMeclw4k9Z07d6rQk4sXLyaQ73e+8x1atGiRCkU5aNAgeuaZZ2jFihU0YMAAmjJlCnXq1MkKDu9Zmfqp61yBX55\/\/dd\/pQULFtCQIUNo\/Pjx6jDmuOdtWgkrmQSBciAAQtdt5XGfCTON2NzjolcV5QJJffv27UpL33PPPenWW2+lfffdl3AOH47Fuu2226hjx46K2KFp\/\/SnP6XevXtbgeI9hRuk\/txzzxEOd2YTTFge0yHKCFmKSQaTzgknnEAHH3ywlTxRMuFwW+CCtrdr1y5K0dTyxpWp3ZYtSqaWujrnssWVybkguyvMmjwK991jqfOHH1KHnj3TanprHy9aRC3\/\/M8U1udZwylr8uj9xhwAHsgKF0C+QFJn4oQmDI0YCQS8bNkyZWrp2rUr7dixg2688UZFxvht07gwwsZzwvKYSF0\/9HfkyJF0wQUXOH9R0N6tW7fSQQcdRO3bt3def5wKo8rU4fnnqcv8+dRh3brC40DsW+fMoeaTTy5cwySJhEkyaooqU9T6o+avlDxBJMoydd2xg44aMCBqkyLnbz7ppNA+rxROfo3JmjyQ0ytT586drdcSI3dajAKRSR3Eec0111BDQwMdc8wx6pGwfb\/44ot0zz33WDUuifmFJ5cgUofJKC1NHQu07777rvoK6NChQwzI3ReJJNPUqdThllt8hWi+6y5avctNburUqbROI30QO0xv\/Sxd5yLJ5B6SkhrLLg9MIMDaM3FCY2b3Q12mzpYOBmlAhT5nz5nXXnuN3nnnHfrmN7+ZifFd9n6zANgrU640dcxIN9xwA33yySfK\/ILPjY0bN9Lo0aNp+vTpysbOpA6yL9dCKWvyY8aMKbKpl+Ow2ObmZoK3D75SskLq1jJZLsZhuU7zhC4a5pjMR40apUxwSPU+G2CsZbJ4iVxkcSoP7N+8GOknXNB9TIyrVlGRTGed5camHhOs1VOn0rTVq2m1Zo9H36K\/bSfymI8OLOa03xwJmEWZ9KaFer+sXLlSaebwTrn44otp7733pksuuYS6d+9O06ZNo48++oiuvPJKpaGznd0GO92DBfWxBg6tv2\/fvmrx0y+PkHoxwtaDzHIxDt7QowM6ES952MtvLZPNYHGQJ6o8\/fv3L0xcJY9vbCTs5VwVQ64t7doRfuBz3tKtG7333nuEz3el0SdxZ4whi15kza5\/\/Jwp8dWGd7ISCTb1AkZl\/jL2O84OXzNQ7Pr06ZMZxS4SqWNRcMaMGfTQQw\/RaaedRvPnz6dHH31UEbqe0PFp2LBNA0m3nWPhdujQoSqbaOoBXw\/QHLt3t34v97DO2ZaRtXhciUqiNo8L89QL+nqIKg+UFq7PJBtIfZON0FoekPmEIUOKzFoRq5DsZUQAB0+\/+OJcWrSojqZNa5vk9XHBX6r4vQr7BDKQQjV1yPj555\/TCy+8oAb5D3\/4Q\/X\/0qVL6e6771ZNuPTSS9UnuXfzUBrt0xdQUf+sWbNo3rx5yj9eSN0dqYP+d0cYidSNGNjQ5KOSaNBDwnbH46sBSkbQ10OYPN4JIw1Sf75DBzqva1flJtytW7dIuErm8iKA9SQosE1ND1BzM5wHYJD030cAUseeniwkK1LPgqAsA4gbJhrY7zGJ6BuXcO9HP\/qRemkwy6aR8DmITy94v2TFpm4tUwRN3V6nL0YZhA6N3VqmkE5atizYKjFqVCMtXer\/suELEp5QfvKAzE3hVBobUWfwtBb1FX5nN6nLhrk03ky3dbKC2NQ0m5qbmUuC3wpWaNxKEr22XJL68uXLlU87Ekj9lFNOUSYYdMS3vvWt6ChICUGgDAhACcACexqkjq9ouBrjN7RGTGR+i9hJmvrqq68qR4m33367UI1Ne+AKPXHixEIZ3WwKJW3hwoV0yCGHKIXg6KOPTiKidVm0BZhNmjSpxMrQRupN6qvTJsFaAfkrnXxJHSYWLIJiwxEWR5Fw7dlnn6Wf\/exn9OGHH6oFnhEjRtDpp59eNp9tgC2kXulhI8+Pg0AapM5kji8Sb8I17PR2mUCETzzxBF1xxRWq2iBi5Od6d4zrcZuQZ82aNcpRwqYuV21h0vaLOxWH1CEbduFXOpWQOoRau3at8lP+29\/+RnvttRdde+21dP7559Mvf\/lLNavhU1ZPCBkAmyZcHtNOYeaXs846S4UQiLNhxkb2v\/\/974UdpcAmCymyTAFeFvB8iWNLL8bBPygVvPlMru6wmyfbJR8WCCvsvrcng5HAQmnUmIgwvzzhWFPHe2cidG6Na+3RROq4dvbZZytS1rV4aPDYywIZMbmY4kKB8A8\/\/HDl7QayRwiQiy66qCSv7hzBRIx4VL\/61a8Kaylwu4bWj+v89WAqB2xgKsEXQbimfiY1N2NnOt6KgCiZuwGHXT2NL6QoPFNC6n\/605\/UoicAPvLIIwkRGjdv3qyu\/e53v6P99ttPuTh+5Stfob\/85S8qPAAWUbGoALDTTrJQWopw2CJgSQm\/lccpU9SwxYsZPyGut78XAGJOsfLILu5JgxO2yhq2tAsLeJRDJcLqC\/Z+MT2JF0ptzBU2+ENLx4JuWHJp6zWZX9glGWQKEsdaF8gUGvjgwYOLiFM3w8AEg6ST+syZM1VcKVgJeB+MaXc5l+Ovdjyb19qwlwbPxuQArsJeGxA4no3E3nJ25hdeKAWph2OdOU0dGjiAefDBB9WMiUUvCAngEAIAi4NYoDzuuOMK44g7Ge6ON910U8FUEzbQktzXZ1\/9BRHvlxgborBS6NlA5OdNEuTi19qfoDLYFO0PbMCjw\/bx2I2VMCfMqJ\/F+kJptG8XtWln1Cjq5\/n8MJE6MA3H1YwAtEx4oYUlKGSwscdJaIuueXo1dSh\/IGLUj69jrHFhbwsSyB7EatLU8a7CEuAldfDOueeeS08\/\/bQKxwGHh2984xtqNzOC+ekTBnznUQevp+lmHJiIzjvvPGM53hPjmtTBl1lwayzS1D\/++OOCreyOO+5QWjkSNHVob0cccUTJBiPsNkVHYmuxbZgA1JlGPHUh9RikHvKm88IbsoFAgrV4mDjc2nDtiCjY3ax1kkniQxw2YRRLWWLL3j1xmsZnmPnErv3p5fK2xUTqbDKBgnXmmWcqrZg1dRCoKQorh9WGlu61qV922WX01FNPqfU63h8D841X40bZMFIvl6aehd23PAqKSJ2J9rDDDiuJmIgt+d7rqIQXPd58882yhQnwi6cupO6e1L104a\/Fg8yjWpldkBG0aJiL\/AIbRP96KJbKPGFwqASvbzzIB\/dMyTQ+86ip+3m\/6KYVmGLhzcLadRTvlwMPPJA2bNigzC\/4GkFdMKlA80bSbepBpI7FXJNNnffTRNHUW78o\/Sf3NBal474dFSH1JAG9guKpl4PU2dsgLZexOB1ZKZl0Ld4yAkGc5gWWafVTD7J1xv96wItcX49NTW2LrNDIvMSt4xAkrOvxiQkF4QzCkkubetizXNwHP3hJ3fbwHRfPRx3cV+3bP0DXXXeyinfmp9AETeSu5IlST8VI3ftZ5CKeOnfEuHHjUtt89NZbb9GECRMozWdE6UDkzYJMa9Zgz0Bb2N6obYiav127LTR27BaCqRif79CUiv2JoaEn+3r4\/vefLyzqQiPE536SxP3kaqEUsoDU9a8Fr3yuvV+StD9PZYO4BP2IozvBYbD5+8WIqVR7q4rUEf8bhKuHjK0UsLX43O3bz6WtW+em0vQOHZ6nrl3Po5aWOgKh+yUs9iMUKnYBNjUtD5QF9XTsuELlbd0K3ppwvUuX+epeGsklqUM+P7t8lkwCaeCYZp262SbsOVDw8JOVZCR1LH6eccYZykcdCXZz+Kjvv\/\/+RddxDz7SWNTAJiW\/0Ltsd8eqOOxscH9EkDA+ii7KcXZB8dQhD4idD3fICsi1JMdrr9XRkiV1CX3OixGDGWT69C101FH+ZG7C+He\/q6PJk\/0P+IC7\/plnPl8oGjZhJO1HjifimtQhl3dHKTxDKhkyNylWlS7PpG4TpycXmnrYIdMmwMt58DSebzoko9IDQZ5fjAB7Sya1t4N8426MDHDJ53MhytZtrm3qZRO8Bh+U574q0tRxKMZ\/\/ud\/quOaoiYc7fa1r33NOlxAGvHUo8os+cuDgM3ZHK2LkqVnJoPMfZxJIgtvcMmPXEeSAnkmiiTtzmPZPPdV7gJ65XGAiMzqTAgK2qiqa+OVJt+0+isVosCMiVVqv4TDLSyPILRpdxoBvfBc27gv+mYnDvxlW9amfZwnlb6KIkCCvELqBvC8Nn6OIoespgpr\/6sAAAUsSURBVOhy3uv6IoueP24\/JZVHlz9pJDzd31g3ufm1Wb9+\/vlLaO1aeGu0IYFFyYsuepv+7d9OjAuP2tzCEQCjyuTXn0lw8pMHgfCwluTUph4WYyGJ7crQI64DeiHmS1hwLRaDv+6x05SjOdqWjTq4dFJHWfaR5\/GirxPivvf0NkSdREoyjqLKzPmF1D3IMQFyJ\/kFEEMxU1x3xL5gd03k0Q\/xiNNJSeUByekx5+PIwGXwQuvtgWwIwYqAb5MnT1Y7\/vQ2428TFq0vwWyaNeviovxxfJFdyeSNzR8XpyB5sBlGDzYV9xlF5TJA6kkCeoEcw4JrcXsROuDYY4+lu+66qxCWwLZsVKyZ1LFbdsWKFSVj+8knn1SL03Dc4Mlm2LBhNHDgQGfvW1SZhdQNiHHEOPg9I5m2OPP2ZtznTmWCx8YU+DL7HeIRtZNcyDNgwICioEZRZQjKj4GPuEAI2nT77beXHFziN\/H5XYfWljTFlUmfjF3G89blgZb++uuvF2vqsDUlCX6DU0SC4r9gQSJm3BfVF7zYsbtjXAf08guuhXU9nJnw3e9+Vx2jyanc5hdwwOOPP248lEcfq8wLab5vtu+GaOoGpNBBTOomTR2fVHAZ85I3X\/eL927bKd58SeXhT0fUq38mxpWHy+nxO0xtxgQX5Tq\/4EnkSiJTGjjp8iCeEjbZFZlfwjTtJGC4KOsx37gO6OVH6njOn\/\/8ZxUlFgG+sPcAqdykji9QfF2ZDuVhePUoku+\/\/37BVOP6fbPtTiH1EFJnrZJPZoHvb6dOnZQWr5\/YwtdB9mmSelR52K+fXwj9pCjbQWLKp9v5WRv1Dvxyk3oSmfQJhe2lfKJWXJy88hhJPYeaun5Ihh4DPU5ALz9Sx+H2vXv3pvXr19PXv\/71QqTIrJG6fuCH90vT1TiKOv6E1C1I3fSZ5dUqdY3MlflF14j5y8GkxcPs4yeP6bpfXbaDh9vKdUddd+CJyXTOrK0MYThElcn7QupfR3FkMmGE8NUl5pc4letlwjT9FBZKXQf0QnN0DxaEe0BE0D59+ih7NXYJw8RXKU0dmyZNY9Ub593UlUnHUZzhIaQeQuq6FoqdqlgUnDNnDuEzizVy\/ToizLlcKGUCZCKOKg\/k5NCmLjZsYZDia0QnQb+DSyCz30KpS4xcyMQHK+gLX7zjOeqL5ScPIp3mfaE0KhZx8kM7h5YLWzps6yDU4cOH7wrd0CVOdbHK8EIpvtDxBWJyAvCLE+\/yfYsjvJC6haauuwTq9lC\/67obnwv3Ne9sH1UePX8Sm7opHoYeBpVt0qaDSwCzzfWog9ilTC5wCpIHSoBzl8Yy+6lH7Z+85ud+xBcDzmj2jm19rHAb2d3RxThKgpuQehL0pKwgEAGBPG9oidDMqsia574SUq+KISiNyAMCeSaKPODrUsY895WQusuRIHUJAgEIMFFkKRa\/dJgZgTRi35cLayH1ciEtz6l5BCTef76GAA7AQOhdhNbNUxJSz1Nviay5R0Di\/eenC7MWJ90WOSF1W6QknyAgCAgCOUBASD0HnSQiCgKCgCBgi4CQui1Skk8QEAQEgRwgIKSeg04SEQUBQUAQsEVASN0WKcknCAgCgkAOEBBSz0EniYiCgCAgCNgiIKRui5TkEwQEAUEgBwgIqeegk0REQUAQEARsERBSt0VK8gkCgoAgkAMEhNRz0EkioiAgCAgCtggIqdsiJfkEAUFAEMgBAkLqOegkEVEQEAQEAVsEhNRtkZJ8goAgIAjkAIH\/D\/juq9Pb\/5yfAAAAAElFTkSuQmCC","height":225,"width":373}}
%---
%[output:8db59851]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXUAAADhCAYAAAA3fSUpAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQm8l1P+xw9CwwwZkRTdMKkZTFkqSgoxsq8TWcqUsSW7MTJFykioQSFDWWKQMYw1VIiyDH9ZmhEyomxjmWHujGv8e5\/b9+fc5z7Lebbf7\/nde87rdV\/33t\/vWb7nc87zeb7ne77LKt9+++23Ksf23\/\/+V73++uvqn\/\/8p+9dli1bphYuXKhOO+00td566+Uoibu0Q8Ah4BBo+giskiepv\/322+r0009X\/\/d\/\/xeK5E9\/+lP1+9\/\/Xv3whz9s+oi7HjoEHAIOgRwRyI3UWQBceuml6vrrr1f77bef2nXXXdXUqVNV+\/bt1YEHHqieeuopNXPmTNW\/f381atQo9YMf\/CDHbrpLOwQcAg6B5oFAbqT+r3\/9S2vpq666qpowYYJae+211cUXX6yWLFmiLr\/8cvX9739fE\/uvfvUrdeWVV6pu3bo1D8RdLx0CDgGHQI4I5Ebq\/\/jHP9QvfvEL1bNnT3XOOefoLvzhD39Q06dP16aWtm3bqv\/85z\/q\/PPPV9\/73vf07xYtWuTYVXdph4BDwCHQ9BEoK6nPnz9fnXHGGerGG29UnTp10uhecskl6vnnn1fXXXed2yht+vMtVQ+ZK9dee23gNTDhde3aVQ0dOlT16tVLrxKlvfHGG2rIkCHq\/fffDzy\/Y8eOqm\/fvuqYY45Rm266qT7uhhtuUBdddJH+e9y4cWrgwIGl88XEeM011+jP1l9\/fT23t9pqq9IxH330kZYHZ4Ctt95amyM32GCDVDi4kx0CYQjkRupo4b\/+9a\/Vl19+qc0vmFv+9re\/6QdrzJgx2sYupA7Zu41SN1GjEIgidfP84cOHK35k9WdD6nI+5D5lyhSteDz33HNq0KBBqq6uTh122GHqwgsvVGussYY+9IsvvlAnn3yyNiNKw8T485\/\/vPT\/K6+8ouf8J5980uj8qP667x0CSRDIjdQR5r777tOa+eGHH66OO+44\/TD88pe\/VDw0F1xwgfr888\/VqaeeqjV0sbMn6YQ7p3kgEIfU0dpvuukmhWcVLQ6pc\/whhxyilQ8xIy5atEhtt912DVaUJmHLCBx11FENTIn333+\/frnQRo4cqY499tjmMViulxVDIFdSZ7OUpesdd9yhevfurSZNmqTuvfdeTehmGz16tDr66KMrBoK7cXUgYJL6jBkz9H6N2ZYvX67Gjh2rIFKaqTWbpL7vvvuq3\/72t3ovRxqa+AsvvKBXl7jibrbZZtqU0rp1a72Zj4LCPpBpOmSP6Nxzz9WrgdVWW03vEWH+wZTIeaZ5hmNuvfVWtcMOO1QH2E7KqkUgV1IHlf\/97396CYvXywEHHKD\/nzZtmhI75PHHH68GDx7c4AGrWjSd4LkiEEXq3FyIlr8vu+wy7T7r1dT9SJ1j\/v3vf5cIfOONN9YE\/qMf\/UhdddVVeiVJwyaO6ZCXAJr8zTffrDp37qy6d++uVwZ4efEbby7xAHv00Uf1MeIgkCtI7uLNHoFcSf29995TrVq10hPdNYdAWgTCSB2tmE1QNHA0dZOUbUidyOcnnnhCnXfeeYrNTVPjxmYuK8mTTjpJmxQ\/\/vhjbVJ86aWXtK38oIMO0t5e7CHJCuHdd9\/V9vS33npLBb1I0mLizncIeBHIjdTRUs4880xN6HgNrLnmmg59h0AqBGxt6quvvromVrT0VVZZpZGmbiPEWWedpVhFcr5Jznvvvbd+cbDpD9ELiaO9i5eLbKgSSS2brMRssKnqmkMgbwRyI3U\/P\/W8O+Ou37QRsCX13XffXdvGa2pqSoDYbpRi+2azk1xEeGzR\/MwoDz74oN4vEjdGPGV+85vf6P0jMbXIMVwDkwz7Sq45BPJGIDdSx+bIQ\/jOO+\/odAHrrrtu3n1x12\/iCNiSOjDg\/XL11VeXiNSG1NG80bbbtWtX0vC5lrnhycoTX3ls92yeQtTY3NdZZ52SPZ9jsJ\/ffffdmuRl03WTTTZp4iPkulcEBHIjdTZEsXHyIC5evFjtueeeaptttin5+JqdxzTDd85EU4QpUVwZojZK0aixp2PuIyvoHnvsUYqR8PN+wUyDHR0Nm7nK\/\/yNn7k3utl0TYT4H3\/8cW0rx0SDqQYzjeniiOsuQXV41LBycC67xZ1XTU2y3EhdzC9RGRoBtBqyNPKSwoth7ty5aq+99tIRi2KvbWqToqj9iSJ15MbGjQsiJGzOqzCXRjZC2QDlRQCZYzNn49NsEjhHqmiziTcMn5nBSCgo33zzjfaScfb0os6opilXbqSOz+7LL7+sfXejWjVo6mhleDrglok7Jt4MbjkdNbLZfh9F6niwzJo1S9vTIWhbUse8wphiJqQRHIev+eabb17qgF\/0qNes4k0bICebxJ8tIu5qDoHGCGRG6vj43nLLLZrofvaznxUea1IT3HbbbQ2CUEwfZ7\/gFglQefbZZ7WrmmykFb6zTUTAODZ1uoyLIcnk0L6jgo94CeCqyGqMRhQ0phgxCZp+6QKneMKYLruYZTDPSON5IC6DF4VrDoFyIJAZqXu9XYrs\/SLkbfoO89BjiyVghaU2hE+eD\/5mtYFHA0U8yOHB+ZhgTE2uHIPV3O8Rh9QJ\/vnd736nNz1pUaTOMdjEeRHgp86LYOLEiWrAgAEl2P\/4xz9q4pfmZ1Yx3R85ztnTm\/usLX\/\/Q0mdjSd27wnkCNO+meyPPPKI9ufddttttYbDctWberdc3eOFwsPHMpyiHNhYd9ppJ70BhoaORw7t6aefLmnqELX8z6qD8zlvww031N4PPOhff\/21zsCHDfarr77SG3GulQ8BG1Ln5QsRM9ZmNkQbUmeccT2UNBZdunTRZhh5MZgboZA+x3pTFZjujyBjbqSWDyl3p+aMQCip22jbPAgEemCSwFf3rrvu0g8BLl54vfAbYg1rfI9WnGWNUvMh9jOlmCRODhD+J5UBy3XpN3\/LQ8vDit0VGXkpkIjMld9rzo+O67tDoJgINCB1zAwk3XrmmWe0tGw8RRGzHLPLLrtoQmTJ+vDDD+tdf9uWl\/cLmp2EjpvJm5ArLqlLX+ivpF617Z87ziHgEHAIlAuBRpr6X\/\/6V202CSsm4BWO7HV4DmDikGaj5efZSUj79ttvV8OGDdMFCqT6ktzTj9S95hfMNyR0cs0h4BBwCFQLAqnNL0EdtbXH5wEULxSi\/ahgI1o5v83iBV5S99so9aZnzUNWd02HQKEQWLKkXhwjxUJJvrDvCtWJ5i1MKKmLrzl28S233LJJIeUldSF\/8mN7M\/w1qY67zjgE\/BCYM0cp6hzwWxrEfuON9f8Ffde3b\/33jvALM68yc2ksTI+cIA4Bh0A8BCDs0aPjnSNHcx4vAr+XgRB+siu7sxIiYEXqmDPmzJlTcgX0uxcJu0g5mndAjlnIADnIseG1lyfEwp3mEGh+CEDG\/frl02+0\/MGDnRafD7qBV40kdYJvTjjhBF3iK6zl5cHivaef6yEV3k17eZkxdLdzCFQvAhC6qWVn3RO0dafFZ41q6PVCSV3S506fPl1XeaEYL\/Z1v0ZyK7T1VVddtawdwG2RvNmO1MsKu7tZU0AAO3il0heIFt8UcCxYH0JJ\/dNPP9VkTt4K6jEWLTWuGTmK6+HSpUvVzJkzdbATld832mijzOHmRYdnD2Ymb3rWzG9meUEnUzRQDqPGGLVYulS1qKDLbt2sWaqud2+FHLQ6nyDFahg3eKAoXACOubk0Rj9m6Y4Q2zqJlyTqkxQARxxxhL4w7oxSVzLdnRqejUcQKQMIQS\/KS87JFD3CDqN6jLwE2nGzzaLBy+mI2h499JVbLlhQugPE\/tH48aq2Z0\/9WTWMG3WYs4yGTwt3KKkTPUkeF3KJ56Gpc93PP\/9c9yGO6QYNnTD91VZbTW+UekmdQKi8NHVeJh988IFeBbRs2TIt\/pmc72SKhrHZY4Rde\/ToRgQqJB+NYHmPqJ0yRW+yVsO4VZWmzjCi\/VJAeuzYsapPnz6ZFoaQqFPuQ\/kvm1wqnEPFd0wtr776qjLzuoim7pfrJaspWVtbqyiUQBRtUUjdyRQ9uqkwyskHO5VM0V3+7og0Lotx7pP1sbNna5MMz9sGO+xQmOeN1CnItEOBZDKhjww+ojAE7oyQLhkKe\/To4Zv7JIlLIzZ7VgI024RevFyQxWwkFJMMjJhfHKln\/XTFv17ZCMtStETyhAXkZOCDnUgmi\/6SmI6mC2\/n6bJoIUuqQ5BfXqj1HaoPhsoA+yRywYNk8OS3NDC+8cYbVd8KyeTXDyubehFL0vllUhRNfcSIEWr\/\/ffPZaOUB3H58uWqTZs2ypskLMlEyeIcJ1M0inExannJJaEBOXVTp6q6I4+MvnHIEXFliroZKaFRehqQjlKKmNCVcZ+Bl0AjrmvXrpF5pu688+rJ1Md0UzdokNI4lbmZ2LNip0Vlgk0rIl52o0MCtCD2wfjkF6CFkrrYvEmvG9XK7dIYRurI6jZKo0Ys3++LtsEVR56W8+ertis33MNQWjZjhlrcvvdKUrHPSirXjCNT1GiRXZWfoAaxR1HO22+9pTcm\/7lwofrB1lv7OgJggze9VMBqvUmTym6rv2fECPXb+fPVAmOTFWIfP358aY9t6dIWgWMT9p0XQ9MBI2wcZs+eXQiNPTL4KGoyBX0vppWf\/OQnulBAULvyyisVYFDH0camznUAGa8XGlkYyaYon2N+cRulSUctu\/OKtsEVR54W\/furFk89FQnGo+0Hqf5Lbykd1759nZo6tc7aOhBHpjBh0MypxBXVZkdo7LWvv67+3aZNYkcAk\/Bb3HKLajFsWJRIib+fppQaEnD2OedMUXPmDFYLFnznyCBjwyko3H7fBaWx6d+\/vy6ME9XQ1NHYK91ikToaO7UcqZJOEBKa\/Oqrr+7bB5vUu\/igYksnNW6cjVKqElH9\/fzzz9faxA033KBfCG6j1G3eBj1Q1vbrmAE5q6jGq1jbuBprmSJYol+\/fg1MLkGHo6mHUs6336qsZNIyBO1JiK0\/JfutEnp+1Cus8clBaWyWLCGNgpHoLOS+NlaNlN2OPN2K1NEoqNkJ8bLrS0qAq666Sl1++eXaFXH48OHqBz\/4gf6f6kBxm1krNOpciBv71oQJExQ1IgkC4v64NfLdfvvtpw4++GC9oZtH40X02Wef6fsWxfvFyRQ90rEwGhKkAza+zwVqtFqiOjT6AoLo0PjjBsfFkimkixRBt22BpA7Rjhqli9vkMr8pISmAzJ2r1DR07XRtRRoytTJRsM+F2EWoT72dTQu\/m9yDdCp6g7qCLZLUiZ4cNWqUog4pZep4k2O7gsDRlCks0b17d13kFw2eYwnOsamaRL+5FuRsW20d4uYFQ65zmrf+6I477lhBON2tHQIOgeaMQFVo6hTXveiii3Tw0aGHHqrt1RArWvvaa6+ty8LxPb7spBSQZmN+STL4jtSToObOcQg4BPJGALdG9gcr3UI19S+\/\/FJrwrjuSUQppg8hdezY7JajsbNkQ3uX1LviOYPNPct0vGJ+ufrqq9W4ceMUG7Je8wsujVIBPmuAWamweiE0uCj5HpxM0aMcCyNMBRH5xd9RNWqaGuxreomSxut+LceLR1xyy8T0UNtv35oadYzX73vUqAbixsIpqqNJvrfCXimMN8GmF7lxHpuW\/qauIvmrW\/mpY6+WnOVeUgc+v8\/ijCcvAFwi+YlqrABOOeUU9fHHHyvqqXbt2rW0yeo2St1GadD8ib0BCLOG2KqxpY9WDQkxau7m\/32wT0gD0oHYA+y+sXHKo1OR2K9450beF5t6HlpzR9\/XCT7smJ6L0KyyNHbr1k1R5g3S9RI4NiQiOl988UV13XXXNUps895776nHHntMffLJJ779JXiAXCrY5G1cGtm05QXz5z\/\/WV9v5MiR6thjj9V\/u+AjFxAVRuqxg8Z88qUsUTXqAjVKa+nFbP6KEc9MWPCM9CXrgKjEGAXkqlk8aJDqEhnwxEalTchVfOl69OjZyDceMh86dGj8i+V0RmQ+dcwuzzzzjCZs3vZeUpciGr1799ZmGNMkAdFjZw8idOkT5xI4ESfTWaWCj3hJ3X333eqggw7KzcQTd6ydTNGIpcVIfLDnz2+pjjiibfQNK3REjx4D1YIFfyjdHUcEzJF4hNm0tDjZ3CPOMcjzwOTJasCJJ5aeN9Jrn3XWWSGXQY\/PXmvu0aNW3XbbMoVMkydPVieulKmqsjSCGsSMU\/1mm22m8A9H6yYfDPZzzB+QPAQ7bdo0hUYvDdcoXgj33HOPHoBevXrpDVUmGZWUXnjhBYVdnDzoeLLEDbkPI3Um8Pbbb69D+bNurCzIV4OfPJkgi9CcTNGjkCVGt97aQl1\/fX0kadHa7Nn1\/tT0N0nofJY4ZYFNkDzz5s3TfMLKQhoK5e67D1UPPTQli1s3uEaLFkvV6NFLVK9edRpbkwNQdivtxmgKG+nSiHkFUwdmDgKPvI2NUMwvBx54YAOb+BdffKFOPvlknX9FNlkvu+wyRR4ZfNwJXuKlgCbP8mXXXXf1HYg33nhD4Yf7\/vvvK9Of3Y\/UAZsXiBk6nPnougs6BFa40tbW9lixST9C1dbW5\/2m8eDX1bWvGD4tW85XbdvW1xNoTg0F0rQQBI3NeutN0mPkN27f\/\/5M9dlnIwJha9UKS4J\/GgZWQvwUpUWSugiK7znk\/sgjj+i3I4E3O++8sw728dMI\/Fwa8XVHOyeUdpNNNtGBDhA+HjRElq6xxhrWuPiRumgokuTH+mLuQIdACgQgcsiCdu65PXMt+RkkJvueY8YsVVtsUS+Ha\/UImGPjxcT73eLF7dXvf9++UUlV9j87d54fCCn8l2RVlNcYWZN6XAGkFB5mEPGcYSMT90NInTS+tKSeM0GkHldOd7xDIEsEKpXpFg\/MgjhfZAlnxa4V4iBUMZlsb5wbqUvVpL\/\/\/e9aO2cTVEwpeNLsvffeSjxnnnvuOevcL7Ydc8c5BCqFQIRHXmKx0Ma9aVNWRvdTJMg1h4BGIJLUMZHce++9mpjJaxDUyAfjTco1d+5cnUWRDVTs6126dNEpASB87PCYdCB4tPYkm6VuDB0CRUUgKJcV5BsR1xTYJVMbr2ZNsqhj1lTkiiT1Bx54QNcDhdzJz8IGp1\/DpuStXsQ506dP1zmOd9llF+0x88QTT5Sux3XY4Jg4caIaMGBAYTAl9QHVY8yAq2uvvVbLJ1WW+BvTkd\/nZv5l8\/ikHUwrjynnxhtvrM1feB0lacjCi5hmvsiD+hz380rKlBVOXoxGjbpRde3aqhRHwUYem26vvbZhqbsdOnyrOnWaoQ44oJu6884fN7DrYq8fNux9NXly9yTw6HPijhvn+BV3zwqjuDKJLPfdd5\/uD7WJTbOuFJxP87zFwShMniwxSjLgoaQugr\/yyivaYwVNO0lDI2fzktzqkDjuSKTLpRE4hLvjqquumuTSmZ8jAyKTRtISsArB7ZK0CZLLnWO9n3fq1Em7fkqOd1IZ4PVjE1jl15m08kC8IrMU6E4KGuYzsz\/IhlfS2Wefrb2jvH3mPn5YBH2eBKOsZDLHNg1OceX54osfqnXW+Yf6xS9+oT3DpBQje0a\/\/OXFaty4+nxKaeZRXJkYB9mzMmXyI\/lyzaVZs2aVFC2RbeDAgYpc51k8b3ExCpIHx5Gsnrek2FqlCaDAqkSUJr2ReR4aPO6RpO0tCpmL5tChQweF2YiGJuCnJYtPqleb53PO9yP7JETBvdPKY076pNp50JhLcjUCsVht2b74uF5WGHllSyqT+TLOEqcoeXjpEk19yCGH6KR4zDlJI10pjIJkgkyFQLPEiDGMwsn7\/IBN1s+bOZeSypPn82bLvaGkLr7maNiy1LG9sBzHRin2copZoG2gEUmiMKpyY7LhpVGkxoQRUvfT1DFhYE7yPnTyeVBq4KR9TCuPLE25v7lsTSqPnGc+WH595oUU53OKh6dtaWTKA6coeaTPXm+usGykeWMUJlMeGMlLXkg66vkxXy5Eq0cdnwQv23Hj2l558sLIth+hpI53CmYX0gSIB4vthTlOUgi8++672tuFStzY5MlySDGNm2++WbVu3VpNmTJlhT2xU5xL53qsSaIy4bCdQ9qk16QgCC85MY2Yn0P2WU+yNPKYL2NZPu+0004qLYGaK5ggAio3qaeRycQjK5xs5Ck3qaeRyXzossJIVsiy6o16mXlNQFHHJyGKOBiFmaSyxChOPyI3Sj\/88EO9LCRYiFwHbdu2tTKZYGKBjPCcIa8LlYi8WRiF9IkmhXyKksrWS6ImoPIG95JinsvBNPL4ySmrkDgTJQyDuPsO8qL0mmuSmKi8qwbpb1yZ\/Jb3aXDyzpMgeeS+fpp61uaXtDJ550vYvLSdW3Fk8jP\/ROFqK0fcecS42ZijssAobh8iszSS4+C1114LdWfkpl6XRgk+MjM8eoWLyvAYtzNZHe81d4jmzWYvm4J485jLPvPz9ddfP5ONGy+JCsGYmomNPMjJHgEvzSwCtsCG1YhJgubkRk7Z1OPvvDdK5QWRViYUjKxwioORbA57xyYI0ySbyXExCpKJuZcVRnFl4nhW+qQUMTHIEqc44xYkT9YYJeE0K1K3Cbv3ujTaVj4CSLPoRpJOZH1OkLmD+4h3gkxKcWk0Pzfd+MzPk8qZVh4xE3H\/NDZ1s1\/SF8nHg5eE2BJtsMgKoyxlygKnJPKAZVTW0TTzKEuZssCI\/saVybyvzD1xX8xiLmUpT1YYJeWLSPNL0guLPYnfuPRhh\/Y2bOuYdkgKxmYq5fFccwg4BBwCDoHkCORG6ogk9U1PO+00nb7XTK8L2ZOu94orrtA+zkcffXTyXrgzHQIOAYeAQ0AjYE3q+JUTNMQym4AT3I8IRiJwyE8L5+Jo4tjAyM7IMZyDXzrl69jt5pqk7MVWlmUdUze2DgGHgEOguSJgRepPPfWUQtv2q2DExiC2rd122823xig+6bfffruOIF22bFkJZ7xoiCYlKsyZXZrr9HP9dgg4BLJGIJLUxe2QQCQqFuF+iK85\/5NfndwuBBbhnrbllluGyofmTjIv8qY7zTzroXTXcwg4BBwCEeYX8TV\/8MEHA0mb6kXkrdhrr70K5WvuBtch4BBwCDRHBKzSBJCd0VtUWsCS6kWk5ZUydSaQmF9IFUCuCIKL0NbJEvj444+X6pX++Mc\/bo7Yuz47BBwCDoHMEbBK6EWgSVjulyBfc9ID4Be9+eaba5fFtdZaS6cHuPTSS0sd4YVRtDQBmaPsLugQcAg4BMqEQCipo1VT1KJVq1al4tFeuagvihb\/2Wef6XzpYiuXvDHXX3+9jsIkCx3eLhSaXm211XRWP0kYtueee+qC0d40AmXCwN3GIeAQcAg0GQSsEnpBzKQH7dOnTwPihbgpenHKKaeooUOH6upGQsx+LwSitvBH51hInIbnzKJFi3xNN00GZdcRh4BDwCFQJgQiNfWpU6eqmTNnKhJ7oVHvscceut4ouV0eeeQR9fDDDyvcGknPi1mF72h+aQJwa7zooosULwm8aGhFTBNQJuzdbRwCDgGHQOYIWNnUO3furN566y1FgWhvIxc6QUW4Ppo1SiWh1\/bbb6\/t8ZI24IUXXiiVU4vaZM28t+6CDgGHgEOgiSPQgNSxj5Mml\/zpNHzKKWSBXzoJu\/ifY6Thn07jGKlBKjZ1jiXD45tvvqmuvPJK9fHHH2vXRzZdpcj0iy++qIYPH+7cIZv4JHPdcwg4BMqHQCNNXfzOSQVg24gOxfRC8QWzYUPH+4UNUpoUme7du7dOH3D\/\/ferDTfcUHvEbLXVVra3c8c5BBwCDgGHQAACVuYXtGuKLRNB+uSTT6ra2lrVsmVLbV\/fZ5991AYbbOB7eTZSn3\/+eTV58mR9Dpuk2OXR9inOSkFqap9uvfXWboAcAg4Bh0B1ILBkSb2cNTWFlDeU1CFfTCSPPfaYTgeADdzbSNTF5ifkHsclUdIFFBIVJ5RDwCHgEPAiMGcOlTqU4rc0iP3GG5Xq27cweEXmfnnggQfUqaeeqqhghP17u+220yl0cVl86aWXdK50NlEpFmFWw+EF8PXXXzdIt+vXa14cZH9kwzUo22Nh0HKCOAQcAk0fAT9NHDIfPTq47xD74MGFwCaU1MVj5Z133lH777+\/DjDyNuzlmGXYSMVV0SyFxf\/HH3+83mj1a6QPQMvHDGN6zhQCGSeEQ8Ah0BCBgpsdUg9XkCYOWYcRutx49uxCaOyRNvVBgwap5cuXq88\/\/zwUM7R3ikyTEoAmfuqbbbZZo3zpaOfU\/Rw7dqz65ptvVLt27dTdd9\/dyDYv1yCHu5RN49rY4++77z59n0Tl2Zr65Ew9u90Fmi0Cfs9GlZgdUo1ZlCZuc3HIH429wi2S1NHQyYOOxk2ov5\/WTbQpphi0bQKRaGySUhyDjVDs7VIIQ7TzWbNmaa2eIsV4wVDf0q+iOz7w++23nyZyNmtZNVBgwyykTE52qSIfimdzmJwVnlDu9gVHIEihCXo2sBVPmxbcqQKZHRIjT9\/79Ut8eoMTv\/02m+ukuIoVqWN2QbP2cztEi8fWToQoxG9ulprEjqcMgUhkcuR6kk\/mZz\/7mfrDH\/5QImrpi7cQr98xHEtEKsQPqVMgm+hXTEE92rRRG220kapr315fsuUll4QuoeqmTlV1Rx6pj22xdKn+Leea+LJXwH4C\/vi4aBahOZmiR6G5Y9TiqadUi7FjG23yMe9bzJtnZ14IgLlu1ixV17t39CAkOKIc49Zyr70a4pJAztIpb79dca+YSFKntihEDZGRr6VHjx6akEmpi3ZOYBHaO37nZGGUxjnrrruuPhebO3VIsb9vscUW6rzzztP1SX\/961\/rlLwQ9tNPP10KSuIakPoZZ5wR6xj84icddJAa8emnqmdtbUkWyFmIOmy8Ph0xQn1v\/nzVcsGCBud+NH68qu3ZU3+G6Yg9ANw4JfgqzRzI4lwnUzSKzRmjVpMmqfUmTYoGKeFcqXvzAAAgAElEQVQR\/zr4YPXiiBH6bBSqLFte44YCSMMpcZM+fbITueiaOlr43nvvrQmWFsdl8ac\/\/Wlp8xONHXMLL4Xdd99de9MQbZqG1GUTF5NMyWyThV0sYHhrp0zRu9tEz4ILqw589YvQwOKDDz7QK5PmLBNmOf2g+vgPN1uM5sxRWhPNua1iXB9iJ2dU3wzc\/LIetzkrTC2jR49WCwzFDWLHEp7aKZH+slla4RaZ0Gv8+PHqnnvuUf369VMdOnQIJHb8zl977TUdZESDXCh+Qek6GsWmqXXKpiemmldffVVHoDIBcIkkspQUvtLCzC9eLV6fk6VdLGBQLlgxaKMNH1XIg4IfWUzeNPMAzNn3ILK3KKReTpl4UNmz4bc079jwMgYjXGebFUbYik2\/6jQTLeTcjkqplSE5paN4Nljpp2lZziPmCIQe1CD2xNIWyF89MvgIEwsPC5ug2NQxvwhRm+CsvvrqOmvjK6+8YjWGaO\/8rLrqqtpMg20cu7vZxF5ubpR26tRJP8CYe8R9Up9ThsnLdtEQn95lMXmtQAs4KMuJn0YO89xyyRT1oPIQM3\/DCB+5w7T8rDDxXicvjEp94YYdodv8m6mpm3ebPXt2aeXkt4KKkiwNRuaYMv4oplENPTuRxs7LYtSoqMuX5ftQUv\/kk0907nO066iGVj5hwgStcds0sjhihkFjP+aYY7SWDsFD5CQHw6RiujSK6yLfE+hktt+dfrraZ\/hwm9umPiZs8qKxV8JbMs3ETw1IBV80tg9qUB95GUM0UVp+tWDku2LJyrQQAgLroyC6BF8hVy4Rd3WbZG774WA7hmjqvk6JaOKsOlj1eCNKIfOUKxJb+WyOCyV1MiviSw7pbrrppr7XY3eaqFP8xiHb1q1bR96Xc9g0RUOHyOM00\/TCebhEXjFihGrVrVucyyQ+1m+ZycX69q2fDpWIIE4y8RMDYHliOWRC8zI1cEvRrA\/LewWWJUZRK5ZUpoUQxDC5sHo1Auet8LXFNgojrxIVhYONcL5OiYYmjkwfPfec2qBApjyzX5HeL6TLZacYjxXMICYJs2FIRCik3r17d11rVIpkyE3wSyfVLp4iEDBBSnjO4HeOnfPCCy\/Udk7bhocL2jrmIK4l\/us9d9xRXyLEo9b2FqHH+ZlfvjvBP\/CA+dChQya3D3yx4iaKm2VR7MW8uPOWaciQ8NHIAnHMN+wl5dGywojYjTBbsciONTnrnryz4uIhwfOhsNlgG4bR9OmNleYlS1aE8zey7scbvUY4sQF6zDGli\/jJlHbvIJ6E4UdHlrPDpk6gD0sooksxx0AcFMzATZE6o3igUMrOm62RwhknnHCCogA1XjS8RXF\/xD2SdLs333yz1uzjFJ6G1PGZ50VBg9TZcP35NdfoEQ4yj2QJmruWQ8Ah4BAwEWB\/sCgtMqEXguI1MGLECJ1Gl4bZBGKG3HfeeWetrTfYtCRwp65Oa9SkDqDwBhusXpdIIX28YXhx2ATzBJJ6mzZ6s9SRelGmlpPDIdB8EKgqUoe88TQh5J+gI2zakHOvXr00cT\/77LPa9EKqAHOTVMrZkd2RVAF+Pu4AQeFp0vted911jUw3flMi0PzSs6daOnasWjxyZPBMwg6yyy7+38+dGxpVx6Iu3G7Innm2u9\/svRirvsB+YeMDb1ZKNi\/Gcjxq5ZDJu8GZX78wq9nnzmbc2Ffz854zp+DixfUmqq5dWzUaN78N9+ApGj07BZsGPUFINvkMv37bcWPlno35yx9bLB7pPDFXPM8xDUNs4sJ1UV46fhhV2q3ZnPuRmvpNN92kNXFzU2q11VbTNUcJwEEDh5DPPPNMddxxx5Wu7Vd42u+hi1t42m+jlPS\/rBQkovS2LbdMtkMdkP9iyJIlEbb6zMIXGkFkk\/gtajMpP7ILvnI5ZErr\/WKPS7ZLa4jdz4lCckH5peyOSsFSv10ZvaOkNUreGAEFHuKMG1Hh6Yk9W2wbjmm8dTs2fkg9qsXBKOpaeXwfGXyETXzRokW6xigaOW8x\/ifRF2YVyB3f8TZt2uj\/pUapRHzyG9L15kqHnNlcwK+doJk\/\/elPjTxngrI0YteXOqokC2OVQIPUSQw2Y8aM+ijTkMkbCaZxbvTkZWslejJE3tPnAJvEb0WcZOWSKXpskqBunsMKrPJRgna9wLEweD2JNonfeFiLO25p3AfrPcLzxNbfV83PzRIyt93sjIuR3dhld1QoqVM0Go8X8i\/gU44nDAQKeeJ9AtHjX44JBjs79nPT9VG0\/NNOO00DhreKNFwl2UB9\/PHHVefOnXXKgLRZGhuRenY4adc5P3\/mmppRas6cxHFoVhJG7cEUcZKVU6agsWHO2XiFBA9Cfiswq4GPfVBQeJy9f3iaccMkA2HaraDKgW1jTd30kxd548KcBqO490pyfCipE75Put0+ffqoyy+\/vBQcJKSOyeOrr77Sfuxo72jb5E+XZtrj0dQBFJdIUgpQ4Jpc6gcddJAiTwybsWyWSkuSpTFPUhe5mAi8kHDxrJ\/A2WXtDBrAqMRvyES5QcYhyh6YZJIkOadSMnkf1CDCl2V2uPkgvxVYEkztzvE3OdiaFrIat+gVVN7Y+odE2eIQhnVWGNmNZ\/yjQkkd8h0zZoxOfoPdHMLw2sDFg2WbbbbRPudeMws+6bfffru64YYbNHHTsOvht44mRdpesjhmlaUR8wueOtj682jvvfeeTkxm3gN\/2Wuuqc\/imEd75pn5oZf1kykPOeJcs2gyIQ8rRjJ\/mnNj7ordRx50yVlU30e0SMxp+a7A4uBpe2zLll1Ube2i0uFsnGOu5IVv07IcNz9skaeubmqu2LZosVS1aDFM1dY+lBiHMKy8GJG\/KuvslDZjFXRM5EYpnimQLxr46aefrh566CH18ssvqyuuuEJr26QGQKvmzYynS1hDc+dFAdHnkaWRICkI18zAlgacOOfW1vZY4YEyYsVE+o7cmVx1delSkbZsOV+1bXtEHFHcsQkRwJsL0mEsly27LeFVKntax471K2XpS2Wl+e7upjzLls1o8JxkLWOrVitSDa9Xn2q4HDig4PFTlBZJ6mjVlJrD3AAhexuuimg6Rx55ZGAGR9k0JZXAxhtvrDdU8ajB3IIdPW6RDN8sjSsFg9glV3KlQIbIIXTaY4+1VyNHJiN2HBTGjFmqttii\/lqulQ+BNONWPikb3gkPmYsvDl\/VVUo2876LF7dXRx2V7JmQ6\/Bs8OOXhqVz5\/JiUHWaumkXJxoULRvTCcUt8I0mDYCfn7o5iKKhm59RXGOTTTbRbpFo7aeccopOyWu2WFkaizBbA2RIWs+2QInfCoxufqIlHbf8JAq+coEyv1p1nwp5SbM8ePuaxsnNStgqOyhSUxcPFmzrhx56qLr00ktL3i9rr7221rLRur1+6uAA6fMdgUveJhFYaPoEzTz44IPa1zxJlkYCmKxqlBZgcLwTMIg4Cpb4rQDIVVaEOOOGpEkJK2kvq1EBSPrSrMa+Jh3XJOeFkrok3sIVEVJHQ\/dulOLuiFsjkXF4yIifOsKwuUqOli233FLtuOOOOt86\/uxo9hTMIO3ugAEDtJZOmbu4GRuTdLjI5ziNo8ijEyyb37glJawwBIhZ4F4Fz\/yaaBDjvDQLlOU2UV\/zPikyS+Oxxx6r1l9\/fZ3Q62186wIaxI+f+uabb66PkBcCOdmvvvpq7RXDiwHyh+jJG3PXXXfpQtTXX3+9NucUpXlt\/GYOd3NVEPS5uFbSnyxWEWnlMeVkT4O0p0nxRhbSPtDMkoVBfY77eZI5kJVMWeHklWfUqBt1KgDBgo1YNvJee+272gMdOny7IohvhjrggG7qzjt\/3IC42Z8ZNux9NXly9xI8cRWAuBhxI7+SkVlhxPVtZerR4+eqTZt\/6+R97MvRpL4Cf2f1vNnKg1XA3Cf0ypMlRkmeh1BSx3yCyQVCpxxdx44dNVljI+dvyTPy\/vvv611mSJrPaX5pAnBtJMMiJI7JhetgtsGzBq+VODVQk3TW5hwZEJk0QblmuJZfCmCia6VgNseQbljSGNjc33tMWnkg3lJ64pXFs5PIwTlvvPFGg\/4gG2N\/9tln64ydBJCZfeZvPyyCPvcmhbORMyuZGqRxToFTXHm++OKHap11\/qED+yhGI9HQPD+\/\/OXFaty4+tQbaeZRXJkYBzOaW2TyrQtsM0g+x8SViRrH8BDOFSLbwIEDVf\/+\/TN53rKSx6zS5g2mTAhV7NNCSZ1CFgjJgzt58mS12267Bfqp9+7dW5thhOj9SB2TCwRjaooQA2YYNHY2YivZeFOTOxv\/WhoTyE9LlgAfmWRC8HzO+b753hMQRRbymJM+qXYeNCaSMZMAsokTJzbKcR\/04gv6PIuHIKlM5ss4S5yi5OGlS5Q2QX4oOOIRFpa4Lu0zklSmMK+zvGUKijbP8nkz+xCFUZA8eT5vthhHml94G5KEn3J1aF2PPfaYIsc69nP81CEwPz91P3s8eV6I4MPHnZeAPNxmhKqt4HkeR5+E1P0eLkwYlNzzkrd87pvv\/ec\/TyxyWnkIyJJmLlsTC7TyRPFO4sHy63Pcz7PY7E4jUx44RckjffZGUAemmE4xj2S808iUB0ZepSjq+TFfLph3o45PMs9tMeLaXnnywsi2H6GkTgEMil+wQUp+FzR3b6PgNHbjAw88sJH5BM8ZTA8syw844IDSsonEQgQysZzjRUE0XxE0dXPSC6nLhKNUH6SN7OwPoFGJacT8HLLPepKZpB5XHjP1giyfdVGRlORgrmCCCKjcpJ5GJhOPrHCykafcpJ5GJvPZzwojrhlHJq8JKI+XXxp58sLIltA5LrLyEWRLRkSSWWE+wYMFEpYCGZhngkJkP\/zwQ+3ZIjnXuRY5SiBINlSxqUv4Ni+PItjUhTRNUjcBlTe4lxTNN3tW5pegl0wcefzkDOqb7cTxYhB330EwblSSMIGJyqt9Sn\/jyuS3nE6Dky1Gcl8\/TT2PeYSJMAqjIJm888OrbNjOn7D5G2Z28jP\/ZG2mijNuNuaoLDCKi2uknzrEjK2PQKETTzxRp8mN43oouV+wy3MdiJwXBHVNabg0klzJWwovbkeyPN5r7hDNm0hVNgXHjx+vzGWf+TmeQllulHpfMqZmYiMPcrJHYG4wid02CWZgw2rEJMGgHPdcP++NUsEnrUzkMMoKpzgYyeawl9TD6gbkPW5BMjH3ssIoybjBG3CFuaGeJU5xxg35\/eTJGqMkYx3p\/UK052uvvRbqzsiNTfe2KEEIPBJTDqaMomjoQZqx6aJUytW+0vuFVQfN\/Nx0sTI\/j8Il6Psg84v3vjZyprGpm\/0SWffdd1\/toornhtgSbbDICqMsZTLxS4pTEnnA0kvqfFZJjIJkygIjb99s5pJ5Xzle3IWzwCnuuIXJkxVGSfnCitRtcqlggiFL43rrrZdUFneeQ8Ah4BBwCKREINL8kvL6DU7n5UA6AEwwRJ7uueee2qxTNE09yz67azkEHAIOgXIikDmpk+CL\/DBsqOL7ix2dIhvz5s1TJ510UgMPGnza+YyfohRMLif47l4OAYeAQyBrBKxJHRs4xIztlE1PdtG7dOmievXqVSqMQfk7ik+b6QSwmePSiJ2V80jW37VrV+3rjifM559\/rqZOnao331xzCDgEHAIOgXQIWJE6roxUjcGTwtvw9mDDgrS55HghpS4h47hNkTqA3OmkBSCjIy6NJnlLAQ582M1o1HRdcmc7BBwCDoHmi0DkRilkTiQodu8TTjhBkzfh\/AQmUYYObZviGRtttJF2d+RzsjO2bt1ao0pOBaJIye\/iDTDCtk4QUtGCj5rvdHA9dwg4BKodgVBSx0cdlzWiuO68806dQtfbcHc87LDDFAmRSKtLUBEubvxPE1etTTfdtMHnfCfRYdjhCURJktCp2gfAye8QcAg4BLJEoAGpYy4hp8s111zT4B5kaERT9\/NSweecH4KHHKlnOTTuWg4Bh4BDID4CjTT1ZcuW6aitjz76SH311VfafEJKAD8tXW7HJupqq62m1lhjDW1mcZp6\/IFwZzgEHAIOgSwQCDW\/vPvuu9r8gr38T3\/6k07s5W1S+QhTDaYW8qk7Us9iaNw1HAIOAYdAfARCSf2bb77RWRZvueUW7dWCv7lpgsHs8sQTT+ikXYMGDdKJv9Dqhw0bprV2Glr8hAkT1IYbbtjgc75jgxV3RjZKnU09\/uC5MxwCDgGHgBeBSJdGsijie44ZhgjQPfbYQ6cCoCoSAUYPP\/ywYhOUZF94wixcuDA2ynHyxsS+uDvBIeAQcAg0IwQiSR3zCpVZ0KjZSA1qW221lSb8r7\/+OjZ86667rvagMYtWx76IO8Eh4BBwCDgEwvOpg48UY6UYBnlaqIK01lprqVatWik2VSH6bt26qcGDB2tS97O7O5wdAg4Bh4BDoDwIhGrqaOmE+C9atEgXtsAuTipeGhkZsblTu\/Txxx\/XAUft2rUrj9TuLg4Bh0ChEKBeL03q9xZKuGYmTGSNUqqcb7\/99urII49U999\/v7rjjjvUxx9\/rDp37qz23ntvtfPOO2uip8wbtnfXHAIOgeaDwJw5c3SxCH5Lg9gpLg8n0Bzhl3c+RJL6scceq80tzz\/\/vK9NnYRd3bt3114xBC45u3h5B9DdzSFQMQT69VNLDDJvJEdNjRpSUxNK+BWTvQnfOJTUyc2CqyK5X7bbbjvtukiKXKrC4OJIWTtcHkkVgL0dLd6F+jfh2eK61uwQ6AdxrzStNOr8kiWqRik1OwAVDDIdA75Di589O+jMZgdzph0OJXUiSg8++GCdSnfixInapZFCFyToOvzwwxVaPAN+4IEHarLHPEOqgLybWXpKSlrlfU93fYdAc0SAYMJAUseGrpR6OwGpcwqkLiaa5ohtXn2ONL+gqS9fvlxRQBpSx34OwZOhkZS5RJFSRJoUvHzGZmqeLctCs3nK6a7tEGgKCORJ6njMYXt3LVsEQkmdNLonn3yy2njjjXVKXdwVx40bp0gfwOeLFy\/W0uDeyDE333xz7jVK0dIp7EoEKpkgf\/WrX+lVA9XtWUUcddRR2SLkruYQaMYImBugQTDUb4f6t++2T\/2\/byqaepFMSaGkThoAMjbeddddasqUKWqLLbYoZWskbS6ujhTO4Lvdd99dHX\/88bnXG4XUb7vtNp1fhgap77TTTrooB9\/tuOOOzfgRdF13CDgEKoEAXFmUFhlRymbpFVdcoXewzzrrLNWjRw+tmVOGjlwv06ZN00UzyOyIl4w0vGGIFGUzNcsWReq4WLZv3z7LWza4lkTVFqmmqpMpergdRtEYcYTgtHRpi5Un1PufBzVs6kEt\/Mz6s3B\/lH3YGp+LVcu4mSU87ZDO76hImzp+6tQljdvyyucSZn6RDVTqoWKOyaORfIxIWqo8kbysCM3JFD0KDqNojDhCcPrPf9qqLl2Y3\/ivBNNzmo3Smpq+qqZmtjK9IiF2zOwrXdxL8rjnzW78OCrSpZG6o6QKQDP\/0Y9+pHersa2TYfGFF15QL7\/8stbQBwwYoDdLpUXlc5GKSLwwSO9rpuuVa\/gdQyIxvG1I9UvDlfKee+7RrpSO1N2LJmjqNyVSZ++InzwamjHBhazG99oLUh8SSurIEOaY2C9UyNErzvYvOD96tFLUojflKbcSxYrfb9VfxLlkwhxK6gDKpuS9996rC0hjevFWP\/rb3\/5Wql16zjnnaNdGm8Z1WXrtt99+DTY7zXP9jiH3zLx587SnDY2iHKT6FZv6EUccoZymbjMC+R5TtIlfNHlMrTiOFgqZYwZdsGBBvgPorq757tJLL21E7EWcS9akjqZM6P+2226rzj333MBydviKM8nID0NBjagmGjgvAcwkrATwheV\/r5YedgzHCvE7UncmobB5V8QHMYlMshqFbFyupSimSf49fIYi66cgJhm35JLEPzNUU3\/zzTe1Bty\/f38Fcfs1tHk0ZdIIYAah8HRUM33NMelA6k8\/\/XSjiklnnHGGTihme4xM+BEjRqj999\/f6gUTJav3ewYUv32CrKS4dtxrZH28kyka0aaCEc\/Z0UcfnetqNBrNpn+EcMkNN9zQaH\/OO5ewTthaKMqBXCipv\/\/++2qfffbRgUck7SHnuenNArmNHj1aF8vApk3JOxvtIQtSx6XS9FEHLDPS9JhjjtGTP+tG5krqt1Jouyhphp1M0aPcVDBiD+r000\/PhdTxQpk+nQRceKUodcwx9b+L2kxPuLgKVlRUunxPPiucPszmnUvsP1A4qCgt0k8dDZ0CGRSWxl595pln6qRdjz76qNaiP\/vsM92X0047TZ100km+\/RICvu+++3SQEsuaiy66SJtbkphfvC8Fuam5NCVXjY0pKO5A0JcPPvhAX7vcGzdBsjqZokexqWDEHCc9R5b7RkLmbE56G5+NGhWNbyWOuOqqq\/SqGcUzjrcb\/IGSihs2L4OxY8cqlEAsAl4uuemmm3SWWrN551JVaep0BD917OnkdaGhnUJmRJjScLrHGwWQ4mRoTLpR2qlTp9KAeJOHOe8X5\/0SRC5FtIMmkSmPOX7BBWrFijuYlgcPrnczzKp5lTxJFUC0OgkCaZheWZGQggSLAQqh6SZNeU1KaZK6hN9EuPs104MOpRRlzHsO5t8OHTo0eDGE4bx4cZ12a95hhw0Ko9iZfY8MPuJgzC9EcfJmNMkcuzL29IEDB6q111471phjr8fzhQHGDo+HDW9NyH6XXXbRAPsdQ2k9NmTNxguF9AB5THhvp5I8iLGASXCwkykatKaCUdZzHC29Y1AqRQNWEiqK73g02uFHQKI09utMrZno9LvvvlsXq0dzxp3QNLGa55lEDC9B1GjavATILCnmGJNPxCHjoIMOavAiMK\/r1dRlRdSvX71ZimYmrRTzFL+LknTSmtR5m5JuF00dFyzelDvssIMaM2aMQnuO25Jq6rLM4v5ki8TMIj7uWU94vz41FXKIO15xjy8aTkWTBzyTyOQ3xyGZoOy4UeOGDX3atKijlEJbx8aepEF4pm2eZ99UzDDJwi\/iECGEj8KHEofiyKpctHM0bm9QpGRrhaDZ88LGvddee2nTsff8OJr6mmvOULW1PSPxpX9vB6WrTAJainMibeoEF7Eseu6550okzlIFTxc+Z2IOHz5cg2y7WZHWpZHBxg6Ghv\/SSy81InXn\/ZJiRmR0atG8TYomj5B6XE8qP++XKPNJRkOa+DJeu7yfZszFeWHNnTtXr\/oh3iBNXQRB06eZ2v5DDz2kyfytt95S\/\/vf\/3SpTVn5i6YOydva1Jctqyd1m5blasbmfkHHBJI6gLz33nu6MAZFMIS4+XzWrFnqxRdfVBSjZpJB\/CT0wh5GJaSoltb7BVMNjZeLJPfiheK8X5xHTtDca8reL9WmqZs2dcYLW\/kJJ5ygEwOSfZXmZ1Mn8pzU3+Si8m5syouC6PZdd91VcxeBWl27di1p9WJTx\/5u6\/0Sh9Sz3nuI4tFYpL5w4ULtnYL7IiABAuCwrIHkn3322dL12PndbbfddGgxSyqbykdpSH3o0KHar\/28887TOWn8SJ3AjLy8XwiSmj59up5URSmy62SKnv5NBaOsvV\/Iu4K9OKpVQgv1ml+iZOR7ssqi4PESp94DZh1TU+cY0fDDriekH4fUuV4RkjU20tSZ\/ESR8ps34\/nnn6\/zQBBkhCsiuWBIbo+5hSUONnXywmAO6dWrl68JJkuXRogUbxyzISeymcFHhPjm0dAACNPGxJPXPeLK7WSKRqypYCT9yNKlEVIPKzVaKQ00Cal7Z4LsveFBEyfJYFJSx65ead\/+RqSOXyZaOiR9yCGHlFID4MIDkeMShD2KzQ1MMYTTsnTiRYBG780NE\/S4pd0o5bre4AOXFyOa3NwRTQOBLEkdRILs8kX2U89zJJOSeuE0dZYsBBSRG\/i6665TrVu31rjxxsTDhMpGYgPDjYhoK5Y6uB5iS+d\/W19103+UQtaS98V0QQo6RgbTL6Iszwx2eU4id22HgA0CYTlJbM4PO8YbUUqWxKzcGNPKVu7zk5A6WBXBrbGBpi4kuummm5Y8Sgg+IvLqj3\/8o8aVghjULYTAMc2QswU7Oho8Wr6NTb3cA+Tu5xBoKgiUw223qWCVph9xSd2bBz7NvdOeG0nqYo7Zeuut1VdffaVDZl955RW9S43r0S233KJdGyF1or4cqacdEne+QyAYgTxJvdpyv+Q5T7ykDmmHxQIUyUzVgNSl0DRpANDEsY8T0UUqAIpM4wGDpo7POqQOgbMBwUYluWEgdaJMXXMIOATyQSAPUq\/W3C9pEnqJWdkv7wvfeUldinaw9+Ct1ERuHDaTi9IakDoeLniz4MBP7VEIms1RMi+SiZHv8Iqhw0LqTz31lD6GCC7yw5jVj4rSybhyePO7mxFwErnGNYM+j\/KBLbc8ppxm9F5cOTgebMT7yPQmCOpz3M8rKVNWOMXFSAgmLOuozLs8SD0qeKlS3i9RcyFpQi+uK6Zm3LQlmtW8n+BMROkhh2ygbr21v\/5axmHRolo1Zco5WpGlefcFJWI27fMWhYHf9428XyBpvFkwt\/Ab0wo5GdgIhciZsELqbKASJPCXv\/xF29lvv\/32QqWgTAKIPNgySPSVz+g7m8Ly4Ampez8nZYLkgecYSVKU1CyVVh6I10sWSXDhHNzDzP4gGyu1s88+W40cOVJvspt95m8\/LII+T4JRVjKZYxsn458Xy7jy0GfTIUC8WsxYDhNTKo1lWd2rErlfKp3QCy89Ik233HJLHfDkDWQyNXXSFMycObPR3CYAUwr7yPiRAyusklvS5y7ueY1IHW2dHC\/8YHbhBzMMD+7xxx+vxo8fr8vJbbXVVlozx2SDlr733nvrDdUiJYuPC4YkCSJUmYZHjp\/WLkFHZrUmcdEkytXvJZCEKLKQhwInZrGRuJiEHS\/LXxIkTZw40frFF\/RCTIKRV76kMpkvY0a\/hXQAABuGSURBVDMFa1q8ouThpQtx4D5MWmtJRx2mTHhJvdoiSiud0Es87STdSBipcyyWC6\/y5p2r8vzn+bzZzkXfNAH4n2M3J3KU3BQUv2CiQw6PPfaYLjiN5o6GTsQpLwLMNd26dbO9b6GPY4CE1P0eLpZURKl5yVs+lyhXroGWvNNOO1lFsQWBklYeSECauUxMOwjmi8yvz2YaBxOLoM9tIv2iZE4jUx44Rclj5i\/BjGmSehCmXlKPMp9EYZb3995NxEon9JI0vTakzgqUyHVcuoOeZ3NVhVUjj3kUZ4wiszSiqVPkmc6REtPbyP+CnYmc6raBR3EErMSxJomKVomNDNLu27evfqHx8MnkND+H7PMk9bjymHVfZdmb9iWDDOYKxtywsiHvvEg9jUzmCyUrnGzkyYLUq1lTN5\/vciX0iqOpR5G6XwU26VNW8yguB0aSulyQDQXInQ1TMt7hIbPzzjtrGxLZ1JpS85K62Tez0LXf51maX+T6aeTxar9h17IdQy8Gcfcd5MUUtaS1lUeuh1lM+htXJr\/lNNc1X4p5yCP39WYujWN+iSOX99hK5H4pQkIvcLDV1NkM9ZurQRXYvLyQZh4lGVtrUk9y8Wo9x2vuEM2baFX2FthXYJnl9zneP1lulAphycQwtWIbeZCTPQLIyUscScbHjPiV84M29fg+741SwUeSNiWViQ3IrHCKg5FsDnvHplwbpeDVlHK\/ZJ3QixU6m9d+TgCSvtfc4JfVRlbPW5Jn1JG6D2pB5hcONXNumLZB83PTjS+LHB1p5THlTGNTN\/slsBGjgL0Ru6PYEm2wyAqjLGXKAqck8oCl3wvXD6M8XBq5fxFzvxQhoRd7hWussUajue3dFwBDcXfMYh4lIXM5x5F6GvTcuQ6BMiOQF6nTDZf75bvBzBPnvKeMI\/W8EXbXdwhkiEA1k02GMOR+qWrG2ZF67tPD3cAhkB0CQjZFyuefXe+Kc6U88taXq3eO1MuFtLuPQyADBFzNgAxAtLwERXCoolZt3n2O1C0H2B3mECgKAq5mQHlGAjKvNkIHGUfqlvODKNtHH31Uu71RrZzSfU0l2MoSAneYQ8AhUAUIpCZ1qiW9\/PLLuqvbbLONWnPNNaug2\/FFpB7rSy+9pIvZ4guLK98mm2wS\/0LuDIeAQ8AhkCMCqUld\/GuRUdLx2sqLH+qTTz6po1Rx7k+Spc\/2Xt7jzCAeMvTRzJSpfv7l5Lgh7w05cYYMGWJdui+pjO48h4BDwCEQF4HUpE65uzvuuEPf97DDDrMiOlIOQKBTp07V+WTiVPmO20G\/44W8JXAGUjdTphJdSLTohRdeqPib1Ujnzp31S4cITc7HBLP55ptnIY67hkPAIeAQyAyB1KRuKwla7vPPP69uuOEGNWfOHJ3ZkbbDDjuoYcOG6URZWaXtNUOs2egwMyWiob\/zzjv63tRXJRoSUoeo5X\/JE895G264oU4\/zIvo66+\/1imHyTlPab899tjDtvvuOIeAQ8AhUBYEcid1yJCi1dOnT1fLli3TnSLLISHlhx9+uLZL57HhiOaNiYQiDn6mFJPEhdS9Se8lDSoysyLBlk7ueF4Kp556alnNRWWZDe4mDgGHQNUjYEXqeH6gmVJoGgJGy8YWzsYhmiymCNMezvF\/\/etftY393nvv1cejhZNvHeKkqELS7HdxEJfKPKKNm+fGJXU597\/\/\/a\/OBeGaQ8Ah4BAoIgKRpI4p4ze\/+Y0uLC3kePfdd2uThphQunfvrqu3oPFSJQQTC6ROo2QUlUV22203fQ0KAZBuNG9Sh7Qpr4dpZ+HChY3u50fqXvMLm7dZVsEp4gRwMjkEHAJNC4FIUr\/pppvURRddpKi\/R9pZtHAqhxBGS1YyzCuQ36BBg9Qzzzyjs\/WRfvbQQw\/VP+QXX3XVVTVq4imTN6lzH3Ig8zKhmeWzZPi8pO63Ueqn4Tet4a9Qb8gcRaupqZAA7rYOgaaLQCipf\/nll1ojZ4NwwoQJ2rPlxRdfVEcffbQ2oZx\/\/vkamTFjxmg7M2TKOQcffLDq16+f1nLNzc9ykbrNcHlJXcj\/3HPP1RWO\/CqM21zXHbMSAT\/ipiIDOV75LQ1iv\/FGpfr2ddA5BBwCGSAQSup+JIxJA82cwtQUm6Zhu0ZLRzPme+qb0tq2bauOPPJIXR0Jovz0009Tm1+8VVPS5AfPAD93CS8CQcQNaU+bFowXxD54sMPTIeAQSIlALFJnkxD7OsWn0WRx78PdDzMMGvx1112nWrVqpd59912F3f3OO+8sebz85Cc\/UQMGDFB33XWXdgVMalM36z7KSwfTUBZFi1Ni6U5PWwF59mynsTeXWUS5JVnN+fWZFRzzwbXYCISSOkE3aOUECF122WU68IaNToJuLr\/8cm2OIbkQ2vLWW2+tg3VMz5Ag33Q2VnEJ3HbbbVN7kgTVDI2NhDshHQK2xS7D7oKmPmpU\/RHO3p5uPIp+dseO4aQeNgcc4YeObuRG6QMPPKAJGH9yNHU2SEePHq03QfE9v\/nmm9UXX3yhpkyZonr37h14MzZU0fA5RzxjxF+dSFRzQ9V2PvoVfpUMdhtttJHix7XyINByr70a2sqzuG1NjaqbOlXVhcyrLG7jrlF+BFp26RJN6gFi1bVvr+reeKP8QgfckX3DrAIns+hUJKmjbZMGAI24trZWHXfccerEE0\/U92YTFTdATDL77LOPVRCR+LBD7gQlsQmbJE2A2NYJYJKK7Gauaez7bOhm3cADn\/211lqrMANZaZlaLF2qNunTJ2uoS9f76NJL1b8OPjjV9SuNkZ\/wTVmmtocfrlqsUACDGnMmTYPYg1pdu3Zq2W23pbl86LneccPkTFBiUVokqQcJii2daM0NNtggsQklaUIvNHRWD\/i9Y\/oRUpeqMCS232677XLR1JH5gw8+0Ndu2bJlIcax4jItWaK05pVjq33wwVT29opj5INNtcvUon\/\/QG07LWmnnUqBpM\/qb9asVJf3jlvVaeqpep\/DyRD6eeedp235r776aoMUAOWoK8hqhXQHePYUhdQrIpPXZXGVVXIYbeOS2NvxkEnYKoJRhKxVL5ONXTzheOV6WtB+jdjqI+IoijhuJl7Wmjq2dKJFKRRBp1hy4ItOigC0dbOxsQrpev3U\/\/73v+v8KST2ghSHDh2qi01IcJLNQI4dO1anHzAb3jd4vzhSL8OLJshlMcyTwWZgbY759lubo3yPKeKDWBUyhXmplGPME494whMhdrNf8gIwNmcXL16sFTuSERZFsYtF6tiPsH+PHz9e27+9jc1OIk6xqdP+\/Oc\/ay+Ympoa7eIotiZS2J5wwgnq7bffbnCJ4cOHK37ibjSIO6OZdMssyrv\/\/vvnYn7hQVy+fLlq06aNTotQhFYumVpecolasUtesS7Xvv56Yq8YwQjhP\/vss4r1wbwxz9bHH3+sFaSikEMjmYYMSbyhWQiQUwjRvq5O8VNqKwme\/FXS4DndCuSRE6mpP\/744+r444\/Xybgg365du+rEXmjjVDzC1ZGqQNdee61ad9111eDBg3UqAQic1AFUQsI1EqLHb53fhxxyiLZLX3DBBTpoadq0adoGHqeFkTrXyWujlL7gycPqpChVnvKSidUWjfTFLefPV22POCLOEDU6ls1ObK0tFyxIdJ2333or0XmcRMQzwW9s+i9IeP\/EN3cnViUCuvD0XXc1JPagnkDqHoW1Up0OJXXxMGFDlAhSPxdBzDJ4xOCnjknl+uuvV1dddZXaZZddSn0iGIk0uPilk1JAyPCVV17Rfu9EnKJxx9HWw0jdbZSmm07ku8dt1SQ\/9BEs2pHB\/ExufrypAEaNUrUDB5YEg9zZzLJ2hezbVy2aMmWlUmSfM8bsC9owc5T50a5du3QgubObNALM\/UmTJqkZy5apnrW1kX1Fd18ye7auC1HpFiui1E9YiSjF9MEDw4sAUl9nnXVKh6PtYz8X27d8QY7y008\/XdvoveeEAcO9cGWkkYWRACmas6nHs6n77QexeoLQgxrEHhnM\/+23ekw\/eu45tUGU3dEiaIkH5oKaGjXNs+wlqjnsIfL2RUjdL79+pR9Ed\/9iIVDikhikfsHgwTrSvtItlNQJKiIjI6XcSHQVVMxCcr988803arPNNiul6JXOYaIhJ4ykFpDPZSXABqptfVMJODrppJN0QjG0flL9ks\/dkbodqQftdQ4ejIbeL3JOErwdqo+8\/baq3Wgjey8hcsJguw1oF6z4POg1w5zC5OdtaOhs5JstT1LHzsreE7+xs2L+K9lbIxGNdwDPG+ZOmpl8zpsXyRv\/4ResxzW85wVdM248iTyPZu+8il28njc82qwpbMrM5wQziqtzknskIfWOSum0KZVuoaSOgHirsPnJJMK26m1sGmJrB0DKvG266aYNSF0yPWKmYeO0devWpUskIXXAZlLzEmCjkgAoCUDiOx5k5PSTNSuwkXv11VePZS7K6t5B14kjk2kZaXw9I4NiwM0wfoQaQFYuQePIpD0OArwpoiTy09YhV3NDi67kQepC5n6rGz4bJWkPMpoAzPG5c+eWcieRMpqXCW6+FKShSR4k77H8z7OMGRVlzSR0M4jPTEM9a6VPt3iXmfeO6pKUjhR5eKmwegKTtEXmzRxQyOGVOStSX3PZMrWRhfkF\/RxSxxEkr5d5FN7yfQNSZ4Pz888\/b\/C2IXoSGyTa9BlnnKE3Ssn5gumEcH+8Yjp27KhOO+00rTl7TSliT8d1ke9NuzluQZJLhtzlbMBGNSYKRaE5nuatP7rjjjtGXcJ930wRyIPUo8xVrCKyXJKbJO71vvKSuDnMvGApZEPW1FtuuUWbLWV1G0bUmEX33HNP7Z4MKVMs\/pRTTmnk+WVqzaaLMRvUfqTOdWS1IZlWzWvIqgCvOV5ErLzY2+Pa7MHRF+mDdzqLpt6pUyfNL9R4oIlc8kKQ+\/uZ40RTh6PgtKiGT18hNXXZfBQQojpifs8gQPhkYTS9WShWMWLEiAapeuU8ySvDJBHNIeqejtSjEHLfByGQNamjpaPQRLXZGW+gmfV3vbn\/TZOHaS7hnIcfflg\/Z6YGbRaQMc0wcl3OMUndT9P2vmhYSYujBLWIzQaBQrbelwP3hsCRmWZeQ5Q4juH+WAZmzpzp+3LhXCF1PNRwqebl5V3RSD3ioJdkElIf0revYqwr3Rpo6mjfuHyhrcdtuDOSihfAMU3gBkkCMBJ9MUF4K4r3DGadv\/zlL7qSEvfElMK5Nk3ML1dffbUaN26cdlMT84+zqQfb1LFsWPCPMQThEaK+lkM2WFeaG9IG1tgS5ncCoyv5GYWw1dfncfcjdT8zjc085BjMHigwUQ1tXapwRR3r\/Z6lfNhyPsykYRIWmq1optxj33331atdFDg\/TV00dC+pQ8ZoyJhlIGG0ZsjWS95oxZhATE3d7JuplYumbr6sOJYXAM28hnjWoRB6NXX6wgsDM5SYX8z9B+mz+b2sYLzX8trUv\/NM9x9BxmhOxMZ93LFPenykn3qcC0PWJOliU1UClShtd8UVV5QyOELCaO7Y3zHFECGK33rQJqz3\/kxiNHuCNjD\/sDqQTVZH6lmSOhqo\/1TuW1OjZnuj7iBzY8MyLakz7rZzon6OhG1QsWE6x5fUo8wnceZ\/Hsd67fLe0owmKUG4aMiyQSjmEuoYcB52dzHZiCaMZmzuS9EHc0MV5Ytm2tRxYSY2AnMtrn8oZGLXN01CXpu6Hz7IDwcceOCBOlpdTDthmjo2eV4qom1zXV4IKImYg5FZXihC7kk19Rlrrql6nnOOVgtwyw5qeeyfJJ1PqUkdIsduzoAwMCx3CM4BRML\/KY5B1Kk02Tj98MMPtYaNrT3Ow8skwKcdOxtt5MiR6thjj9V\/u4jS8CjX730vTgIyf02dDWjIQ29OQuwBeTSyiHLt37+\/fvlHN3xxwpa9+pFsMpq6qX2KRguRB3mxmMQsWJokx2cQO6ZSWpT3C04PuCkTnIPNmefe1LxFJn77aepeMy+aOmRuysBnaL+QsmlTN+3fUd4v3F9WELzsiE3AhVpeCKxcgkpXCpfgWScvSb\/4DZ4HXjK4bBelJSZ1wonnzZunJk+erMvXxXV3SgtApSJKmdBUdaJGa1ECWGxlOvzwtmrBAhtix9\/E37WRVRY\/Uc1WprDr+LnENT7eNixqlcy9X\/zcJv36k7VNPQr7vL9fuHChjgTHYQJvjzPPPDO36GobbT9uf21cHmXuUQxI7PxyH+Y2vEcKcjig6lPvooVjYmG5xVuahibOm4qfcuVDCSN1Cl9vv\/32Oj9L1o3QefLH4ycfN7VB1rLI9WxlmjevhRo5MriQCddr0WKp2n33W1cEhI5usOuPqYzxtS0baCtTFCa33nqrjlIObniwr6yWFHqxjqply+U6QC7L4CNcaCH3oJa190sUXk3t+0qTut9z7p3bUfse5R4TK009qCwdWcpwGdp5550bkDm2NlwhcVE0TSuYXlj6UM+Uh4skYFH+5ObmiWx08OLwI3WzSEa5gayW+\/3rXwerjz66NFDcVq0mrTChTSp9z9jHSd+QBw6YctiLMV3LkKmubqpNfKsWqWPHzZT0JUtS59pBdvki2VnzGJemfE27VWI9Arar13LhFUrqJO1ip5ilhqmV83Cxm+8XZcrSjKyNRHqaof\/Y0LGbPfvss6W+4SXDLjk2uTh2dS7gR+p8LuXsygVgNd5n8eL26ve\/b++XnkV17jy\/0F0yXzLnntuzQR+CBMf8f\/HF8\/Wmns7nMWNGqmhDv\/t4I0qx4RYhD0ihB7PAwgmp2+QJyjvYMS5MjUjdr9wc5IuWjJcKRafZ\/WbzgA1LszGx+Y7fuDqhqWB340HkYcINkeUo2j2ZHUnuhfskbmFdYlbOCSL1uAA09+ND9joLD41F2hi9j0s6Doi9HN5RhQfNCWiFQDXPlQakbrob0nPMK9T57NOnT8mDRcjUj9RvuukmraV73RQlcpQ3GnlgxBsGzwZS9PIiICgirrZuNTruoCaNQETaGJ3+XSL183xQy5n7pUkPaEE6l+dcybuLDUhdCBsfVGqA7rTTTlrTNlsQqXMO2RLZDffmeAEgXg64LIn7IdeUhGEEhbDL7L1X3p13128aCAQlKPO4zueiqZc79wsjVs6EXjJDwlIQZDmLzL6JRx17aGFpAbK8v9lf3CHzMNXlIa95zQakjp8rngb8YE9nMwq7IERMLvQ11lijZMv2aupC9t6EXtwM2zo+oWjyFNuQliShV96AcH1vsiBzoplZ5oI+NzdZsshKl1aeIBJIgqVffg4zQybXNPschEWWGJky\/ehH\/dVtt13cIGunyITPc9YPalTwUtbeL15yzTuhF2Mr+ErkZ5J5Y3MO8xRPEjMhGSkC8DaDk4JyvdhcO+4xpqbOueLvLnPbGxNgYpPl8xZXbo733ShF4CeffFInInrhhRe0TRyTCa6CuHARtdW7d+8GNvUgUhcN\/vXXX9eRn3i9FJnUZUBkkCQtgTcrpGhM3s8JUybxmeR4J5UBJqekWenSyuMXMZhkonCOmQmP\/iAbSZaI5CMIzNtnzvHDIujzJBjFkYkHk3HNSvuyTWWQpZ96uRN68Vy\/9NJLaq211vJNJyBzyS8fDdHk3oRakoyLeUOwk2jjXCcoWZgZNYv3HA4aNPNcVvr4j5PywPSSSzrXhdRZIZBnxju3zYhW4b6BAwfqvURvhG5SGZKeF+nS6OeXzs2oAUpHpeh0kClF7OnbbLONLmWHti+NUH\/s6aTjLYL5RYISmBg0NoL9tGTJxWGGKYuWgTYYlBo47iBlIQ9RmUKsZNrLsklyNQKxJk6c2CgdctCLL+jzNPmvpV9hMkEC7PmYpF6NuV\/KmdBL5oy5QuD+rLwJvKHeKwm\/iCr3prt98803GyXUIpsrPCCpfiXtAZYASTjmnaNmhkk4w5v0y1SkkFdSDKSZT0LqcMBDDz3km+rblFOe\/zyfN9tnN5LU5ULiq04WRt6wkttFfNVJeTthwgQNgOnNAkGybMLTxRu4IhulRx11lDrrrLMKs1HKAAmp+2nqhBbjsuYlb\/k8KDWw7aB4j0srj5lsKcsltPki8+szD3mcz20Dm8JwDJOJhHJogyapR5lPko5ZVudF+brnndAriNQhdnLKcH\/4AO1YtHIz9N6bUMtrSpHskcxRP02d5w8SF5v6J598onOwoOnT\/LI+ZhGwJKTOCpREYX6pvmWMzVw5yJfX82Y7p6xJ3bwg9nbIm7wIJNWSpRBaN0tx6pXyhuNFgGaEFo+Zwsw2x1IJDfLVV1\/VmRwx5xSlmSQqWqXkiWCPAVMU\/ZMJyySWzyH7PEk9rjym26nYAdkAT0ug5gpGtGPvxC83qUfJ5Efq1aaplzuhV5SmzrONycHcKxOtXnK3oDGL2ShIU2c++tnU+Yw9OXLAY\/41U+6KRp6nph5F6vJMmUVGhMeyfN7icGMiUpcb4NNOwh7S9RIlytuUsG6KVEvzZmIk1S6mFgaHlwNJvfipdNSidyklmrqf1mxu5sj3eZhfzGsnlcdL3t4XVpzJ4u2ruaHlZ3KSF5DtfkSa5XIQIZj39jO\/JOm\/nFOp3C\/lTOglffWaX0RTl+\/9kmuZWqsk1MKFGYIWm7rX\/h3m\/UKBD\/bzJPGYJP3C5GHa1LNYjZqaOvfzq7QWVB7QlkvSzL2wc1ORunlhtHKyMtJI8EXkKH7nAM6bVHzQpdg00X1EmDJQRJ8WqXnNHaJ5E60q6T2ZsH6fszmU5UapEKOQuqkV28iDnJIvO4uALT97pTm5kVM2h\/k7741SwcdMOctnfjJlvVHKfZpj7hcI3Uvqts9vUP5y2\/P9jguryJT0ukLqrNAx8fg5AfgVDDFfflk8b0nkz4zUbW9Oql40dHzTzU1T2\/PLcVyQ+YV7m\/ZYU6swPzfd9bLwtEgrjylnGi3GLx+GWWxBbIk2WGSFURyZGL+sXRq5psv9Yv9UVhupsz8IT3nntne1BALi7pjV82aPasMjy07qSQV15zkE0iKQZ5Sgy\/2SdnSKdX6ecyXvnjpSzxthd\/3CIFDND2phQGwmglTzXHGk3kwmqetmw8pYVO1xzSEQhADeebhZZ2E+LTfKjtTLjbi7X8UQcPn2KwZ9Vd6YFz+pd6NqPhStc47UizYiTp5cEXD59nOFt0ldvGh50m3BdaRui5Q7ziHgEHAIVAECjtSrYJCciA4Bh4BDwBYBR+q2SLnjHAIOAYdAFSDgSL0KBsmJ6BBwCDgEbBFwpG6LlDvOIeAQcAhUAQKO1KtgkJyIDgGHgEPAFgFH6rZIueMcAg4Bh0AVIOBIvQoGyYnoEHAIOARsEXCkbouUO84h4BBwCFQBAo7Uq2CQnIgOAYeAQ8AWAUfqtki54xwCDgGHQBUg4Ei9CgbJiegQcAg4BGwRcKRui5Q7ziHgEHAIVAEC\/w\/DWCrhZGAQeAAAAABJRU5ErkJggg==","height":225,"width":373}}
%---
