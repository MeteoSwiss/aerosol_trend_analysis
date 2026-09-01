MLO_st.name='MLO'; 
MLO_st.lat=19.5362;
MLO_st.lon=-155.576;
MLO_st.alt=3397;
MLO_st.env='Mt';
MLO_st.footp='Mx';
%%

if ispc
MLO_rd=read_betsy_2026('C:\github_trend\raw_data\mlo_time_mri_mse_tsi_psap_clap_AE31_AE33',MLO_st.name);

elseif isunix
%MLO_rd=read_betsy_2026('/proj/pay/aerosol/matlab/github_trend/raw_data/MLO_all/mlo_1974_2025',MLO_st.name);
end

%%
% compute the AE33 abs coef since only the eBC in microg/m3 are reported
% under X*_A82_ae33
% as given by NOAA documentation for AE33
% https://gml.noaa.gov/aero/docs/ae33_sop_20161227.pdf, abs=BC*sigma
sigma=[18.47 14.54 13.14 11.58 10.35 7.77 7.19]; %for the 7 wavelengths with an exponent of -1
MLO_rd.Ba1_A82_ae33=MLO_rd.X1_A82_ae33.*sigma(1);
MLO_rd.Ba2_A82_ae33=MLO_rd.X2_A82_ae33.*sigma(2);
MLO_rd.Ba3_A82_ae33=MLO_rd.X3_A82_ae33.*sigma(3);
MLO_rd.Ba4_A82_ae33=MLO_rd.X4_A82_ae33.*sigma(4);
MLO_rd.Ba5_A82_ae33=MLO_rd.X5_A82_ae33.*sigma(5);
MLO_rd.Ba6_A82_ae33=MLO_rd.X6_A82_ae33.*sigma(6);
MLO_rd.Ba7_A82_ae33=MLO_rd.X7_A82_ae33.*sigma(7);
%%
%remove the parameters with g and N correMLOnding to nb of data in hourly
%mean + STD
names=fieldnames(MLO_rd);
% % c=contains(names,{'g','N','X'}); 
% % names_delete=names(c);
% % MLO_rd=removevars(MLO_rd,names_delete);

%%
first_date=MLO_rd.Time(1) %[output:4600e4b1]
last_date=MLO_rd.Time(end) %[output:3ccf89c6]
names=fieldnames(MLO_rd);
c=startsWith(names,'Bs');
names_sc=names(c);
c=startsWith(names,'Bbs');
names_bsc=names(c);
c=contains(names,'A11');
names_clap=names(c);
c=contains(names,'A81');
names_AE31=names(c);
c=contains(names,'A82');
names_AE33=names(c);

%%

 prctile(MLO_rd.U_S11_tsi,[5 50 95]) %% < 32%
  prctile(MLO_rd.U0_S11_tsi,[5 50 95])  %% < 19%
   prctile(MLO_rd.U1_S11_tsi,[5 50 95]) % --> all RH < 19% OK
plotFigControl(MLO_rd,MLO_st.name);
% 2 size cut, since 2000
%Bs/N and Ba/N: clear seasonal cycle
%%
% homogeneisation MRI-MSE
P_mri_mse=timerange('1994-06-01','1998-06-01');
ratio_STN=nanstd(MLO_rd.Bs2_S11_mse(P_mri_mse)./MLO_rd.Bs2_S11_mri(P_mri_mse)); % 61.82 very high !
ratio_mean=nanmean(MLO_rd.Bs2_S11_mse(P_mri_mse))/nanmean(MLO_rd.Bs2_S11_mri(P_mri_mse)); %0.5601
ratio_median=nanmedian(MLO_rd.Bs2_S11_mse(P_mri_mse))/nanmedian(MLO_rd.Bs2_S11_mri(P_mri_mse)); %0.8407
figure;
plot(MLO_rd.Bs2_S11_mse(P_mri_mse),MLO_rd.Bs2_S11_mri(P_mri_mse),'.');

% fit between both data: slope 0.9955x+0.1021
fit_test=robustfit(MLO_rd.Bs2_S11_mse(P_mri_mse),MLO_rd.Bs2_S11_mri(P_mri_mse)); % 0.1368, slope=0.8725
% --> use a cte=0.85 (middle between median and robust fit)
cte_mri_mse=0.85;
%%
% homogeneisation MSE-TSI
P_mre_tsi=timerange('2000-05-01','2001-03-01');
ratio_STN=nanstd(MLO_rd.Bs20_S11_tsi(P_mre_tsi)./MLO_rd.Bs2_S11_mse(P_mre_tsi)) % 8.8 very high !
ratio_mean=nanmean(MLO_rd.Bs20_S11_tsi(P_mre_tsi))/nanmean(MLO_rd.Bs2_S11_mse(P_mre_tsi)) % 1.4932
ratio_median=nanmedian(MLO_rd.Bs20_S11_tsi(P_mre_tsi))/nanmedian(MLO_rd.Bs2_S11_mse(P_mre_tsi)) % 1.5152
% figure;
% plot(MLO_rd.Bs2_S11_mse(P_mre_tsi),MLO_rd.Bs2_S11_mri(P_mre_tsi),'.');

% fit between both data: slope 0.9955x+0.1021
fit_test=robustfit(MLO_rd.Bs2_S11_mse(P_mre_tsi),MLO_rd.Bs20_S11_tsi(P_mre_tsi)); % 0.1368, slope=1.4996
% --> use a cte=1.5 (middle between median and robust fit)
cte_mse_tsi=1.50;
%%
% break point detection: decision in 2019 to remove data before 1988
% it has to be done on the merge Bs2_S11 (TSP) and Bs20_S11 (PM10 in between)
% BP in March 2020--> try old neph 10% hibher --> old neph=old_neph/1.1 + BP analysis
P1=timerange('1988-01-01','1995-01-01');
P2=timerange('1995-01-01','2000-04-29');
P3=timerange('2000-04-29','2022-11-28');
P4=timerange('2023-04-07','2026-01-01');
% B
MLO_rd.BsB_S1S10=NaN(size(MLO_rd.Bs1_S11_mri));
MLO_rd.BsB_S1S10(P1)=MLO_rd.Bs1_S11_mri(P1).*cte_mri_mse.*cte_mse_tsi;
MLO_rd.BsB_S1S10(P2)=MLO_rd.Bs1_S11_mse(P2).*cte_mse_tsi;
MLO_rd.BsB_S1S10(P3)=MLO_rd.Bs10_S11_tsi(P3);
MLO_rd.BsB_S1S10(P4)=MLO_rd.Bs1_S11_tsi(P4);

%G
MLO_rd.BsG_S2S20=NaN(size(MLO_rd.Bs2_S11_mri));
MLO_rd.BsG_S2S20(P1)=MLO_rd.Bs2_S11_mri(P1).*cte_mri_mse.*cte_mse_tsi;
MLO_rd.BsG_S2S20(P2)=MLO_rd.Bs2_S11_mse(P2).*cte_mse_tsi;
MLO_rd.BsG_S2S20(P3)=MLO_rd.Bs20_S11_tsi(P3);
MLO_rd.BsG_S2S20(P4)=MLO_rd.Bs2_S11_tsi(P4);

%R
MLO_rd.BsR_S3S30=NaN(size(MLO_rd.Bs3_S11_mri));
MLO_rd.BsR_S3S30(P1)=MLO_rd.Bs3_S11_mri(P1).*cte_mri_mse.*cte_mse_tsi;
MLO_rd.BsR_S3S30(P2)=MLO_rd.Bs3_S11_mse(P2).*cte_mse_tsi;
MLO_rd.BsR_S3S30(P3)=MLO_rd.Bs30_S11_tsi(P3);
MLO_rd.BsR_S3S30(P4)=MLO_rd.Bs3_S11_tsi(P4);

%4
% % MLO_rd.Bs4_S4S40=NaN(size(MLO_rd.Bs4_S11_mri));
% % MLO_rd.Bs4_S4S40(P1)=MLO_rd.Bs4_S11_mri(P1).*cte_mri_mse.*cte_mse_tsi;
% % MLO_rd.Bs4_S4S40(P2)=MLO_rd.Bs4_S11_mse(P2).*cte_mse_tsi;
% % MLO_rd.Bs4_S4S40(P3)=MLO_rd.Bs40_S11_tsi(P3);
% % MLO_rd.Bs4_S4S40(P4)=MLO_rd.Bs4_S11_tsi(P4);

names_sc_short={'BsB_S1S10', 'BsG_S2S20','BsR_S3S30'};
break_MLO_scat=change_point_analysis_Def(MLO_rd,names_sc_short(1:3),0.05,'MLO','neph SC');
T_MLO_scat=make_table_breakpoints_def(break_MLO_scat);

% % break_MLO_scat_snht=change_point_analysis_SNHT_V1(MLO_rd,names_sc_short(1:3),0.05,'MLO','neph SC');
% % T_MLO_scat_snht=make_table_breakpoints(break_MLO_scat_snht,names_sc_short(1:3));

%%
% backscatter since April 2000
P1=timerange('1988-01-01','2022-11-28');
P3=timerange('2023-04-07','2026-01-01');
% B
MLO_rd.BbsB_S1S10=NaN(size(MLO_rd.Bbs1_S11_tsi));
MLO_rd.BbsB_S1S10(P1)=MLO_rd.Bbs10_S11_tsi(P1); 
MLO_rd.BbsB_S1S10(P3)=MLO_rd.Bbs1_S11_tsi(P3);

%G
MLO_rd.BbsG_S2S20=NaN(size(MLO_rd.Bbs2_S11_tsi));
MLO_rd.BbsG_S2S20(P1)=MLO_rd.Bbs20_S11_tsi(P1); 
MLO_rd.BbsG_S2S20(P3)=MLO_rd.Bbs2_S11_tsi(P3);

%R
MLO_rd.BbsR_S3S30=NaN(size(MLO_rd.Bbs3_S11_tsi));
MLO_rd.BbsR_S3S30(P1)=MLO_rd.Bbs30_S11_tsi(P1); 
MLO_rd.BbsR_S3S30(P3)=MLO_rd.Bbs3_S11_tsi(P3);


names_bsc_short={'BbsB_S1S10', 'BbsG_S2S20','BbsR_S3S30'};
break_MLO_bscat=change_point_analysis_Def(MLO_rd,names_bsc_short,0.05,'MLO','neph BSC'); %[output:35a82c92] %[output:774c0d08] %[output:2215ae78] %[output:151eb4ae]
T_MLO_bscat=make_table_breakpoints_def(break_MLO_bscat); %[output:6679b005] %[output:9b4abf9a] %[output:749b9ba6] %[output:7e40b023] %[output:349fc67e] %[output:8ba12aeb]

%%
% Absorbtion: two times series needed. in 2019 only clap/PSAP were used, no
% AE available
% 1) for PSAP+CLAP: 2000-2022
% take here only PM10
% consider B insteat of G

%%
% missing value still attributed to value>1e5
names_abs=names(startsWith(names,'Ba'));
for i=1:length(names_abs)
    MLO_rd.(names_abs{i})(MLO_rd.(names_abs{i})>1e5)=NaN;
end
%%
% homogeneisation AE16-PSAP
P_ae16_PSAP=timerange('2000-05-01','2000-08-04');
ratio_STN=nanstd(MLO_rd.BaG0_A11_psap(P_ae16_PSAP)./MLO_rd.Bac1_A81_ae31(P_ae16_PSAP)) % 1.07 %[output:6fbf83c4]
ratio_mean=nanmean(MLO_rd.BaG0_A11_psap(P_ae16_PSAP))/nanmean(MLO_rd.Bac1_A81_ae31(P_ae16_PSAP)) % 0.9062 %[output:63c6e05a]
ratio_median=nanmedian(MLO_rd.BaG0_A11_psap(P_ae16_PSAP))/nanmedian(MLO_rd.Bac1_A81_ae31(P_ae16_PSAP)) % 0.9167 %[output:33ae1221]
figure; %[output:784e460e]
plot(MLO_rd.BaG0_A11_psap(P_ae16_PSAP),MLO_rd.Bac1_A81_ae31(P_ae16_PSAP),'.'); %basic fit: y=0.807x-0.036 %[output:784e460e]
xlabel('BaG0_A11_psap'); %[output:784e460e]
ylabel('Bac1_A81_ae16'); %[output:784e460e]
grid on; %[output:784e460e]
legend('period 1.5.2000-4.8.2000','fit','period 1.5.2000 - 31.1.2001')
% fit between both data: slope 0.9955x+0.1021
fit_test=robustfit(MLO_rd.BaG0_A11_psap(P_ae16_PSAP),MLO_rd.Bac1_A81_ae31(P_ae16_PSAP)) % -0.0259, slope=0.8489 %[output:99f23bd4]
% --> use a cte=0.85 (middle between median and robust fit)
cte_ae16_PSAP=0.85;
%%
% homogeneisation PSAP-CLAP
P_PSAP_CLAP=timerange('2012-01-01','2013-08-07'); %[output:03ccf105]
ratio_STN=nanstd(MLO_rd.BaG0_A12_clap(P_PSAP_CLAP)./MLO_rd.BaG0_A11_psap(P_PSAP_CLAP)) % ^NaN
ratio_mean=nanmean(MLO_rd.BaG0_A12_clap(P_PSAP_CLAP))/nanmean(MLO_rd.BaG0_A11_psap(P_PSAP_CLAP)) % 0.7094 %[output:29764ddb]
ratio_median=nanmedian(MLO_rd.BaG0_A12_clap(P_PSAP_CLAP))/nanmedian(MLO_rd.BaG0_A11_psap(P_PSAP_CLAP)) % 0.7000 %[output:042c7f36]
figure; %[output:093f4a51]
plot(MLO_rd.BaG0_A12_clap(P_PSAP_CLAP),MLO_rd.BaG0_A11_psap(P_PSAP_CLAP),'.'); %basic fit: y=1.095x-0.089 %[output:093f4a51]
ylabel('BaG0_A11_psap'); %[output:093f4a51]
xlabel('Bac1_A12_clap'); %[output:093f4a51]
grid on; %[output:093f4a51]
legend('period 2012-01-01,2013-08-07','fit') %[output:45a6b617] %[output:093f4a51]
% fit between both data: slope 0.9955x+0.1021
fit_test=robustfit(MLO_rd.BaG0_A12_clap(P_PSAP_CLAP),MLO_rd.BaG0_A11_psap(P_PSAP_CLAP)) % 0.048, slope=0.9701 %[output:2f2715da]
% PSAP 3w-clap 3w was considered as homogeneous in 2020 and the slope seems
% around 1 --> consider it as homogeneous
%%
% homogeneisation CLAP-AE33
P_CLAP_AE33=timerange('2014-07-01','2022-11-29');
test_AE33=MLO_rd.Ba3_A82_ae33(P_CLAP_AE33);
test_clap=MLO_rd.BaG0_A12_clap(P_CLAP_AE33);
ratio_STN=nanstd(test_AE33(test_AE33>0.1 & test_clap>0.1 )./test_clap(test_AE33>0.1 & test_clap>0.1)) % 2.452.810.355
ratio_median=nanmedian(test_AE33(test_AE33>0.1 & test_clap>0.1 ))/nanmedian(test_clap(test_AE33>0.1 & test_clap>0.1)) % 2.70
ratio_medianV2=nanmedian(test_AE33(test_AE33>0.1 & test_clap>0.1 )./test_clap(test_AE33>0.1 & test_clap>0.1)) % 2.54

figure;
plot(test_AE33,test_clap,'.'); %basic fit: y=0.22x-0.038

hold on
plot(test_AE33(test_AE33>0.1 & test_clap>0.1 ),test_clap(test_AE33>0.1 & test_clap>0.1),'ro'); %basic fit: y = 0.2304*x + 0.1266
xlabel('test_AE33(test_AE33>0.1 & test_clap>0.1 )');
ylabel('test_clap(test_AE33>0.1 & test_clap>0.1)');
grid on;
legend('period 2014-07-01,2022-11-29','fit')
% fit between both data: slope 
fit_test=robustfit(test_AE33(test_AE33>0.1 & test_clap>0.1 ),test_clap(test_AE33>0.1 & test_clap>0.1)) % 0.1017, slope=0.2239


ratio_STN=nanstd(MLO_rd.Ba3_A82_ae33(P_CLAP_AE33 )./MLO_rd.BaG0_A12_clap(P_CLAP_AE33)) % ^NaN %[output:5909ecc4]
ratio_mean=nanmean(MLO_rd.Ba3_A82_ae33(P_CLAP_AE33))/nanmean(MLO_rd.BaG0_A12_clap(P_CLAP_AE33)) % 3.53 %[output:4717ed11]
ratio_median=nanmedian(MLO_rd.Ba3_A82_ae33(P_CLAP_AE33))/nanmedian(MLO_rd.BaG0_A12_clap(P_CLAP_AE33)) % 3.64 %[output:78ec53d4]
figure; %[output:8654cf04]
plot(MLO_rd.Ba3_A82_ae33(P_CLAP_AE33),MLO_rd.BaG0_A12_clap(P_CLAP_AE33),'.'); %basic fit: y=0.22x-0.038 %[output:8654cf04]
xlabel('Ba3_A82_ae33'); %[output:8654cf04]
ylabel('Bac1_A12_clap'); %[output:8654cf04]
grid on; %[output:8654cf04]
legend('period 2014-07-01,2022-11-29','fit') %[output:5e83376a] %[output:8654cf04]
% fit between both data: slope 
fit_test=robustfit(MLO_rd.Ba3_A82_ae33(P_CLAP_AE33),MLO_rd.BaG0_A12_clap(P_CLAP_AE33)) % 0.0055, slope=0.28 %[output:8766cbe6] %[output:272c8b13]
% clap 3w - AE33 was considered as homogeneous in 2020 and the slope seems
% around 1 --> consider it as homogeneous
figure;
normplot(MLO_rd.Ba3_A82_ae33(P_CLAP_AE33)); %[output:8d098593]
hold on;
normplot(MLO_rd.BaG0_A12_clap(P_CLAP_AE33))
%v --> it is not possible to homogeneise easily CLAP and AE33 
% --> I will do 1) trend with AE16+PSAP+Clap until 2022, 2) AE16+PSAP+Clap until 2022+Ae33, 
% and 3)E16+PSAP+Clap until 2015+Ae33

%%
   % time serie AE16+PSAP+CLAP
P1=timerange('1988-01-01','2000-04-29');
P2=timerange('2000-04-29','2013-01-01');
P3=timerange('2013-01-01','2022-11-29');
P4=timerange('2022-11-29',last_date);
% 
MLO_rd.BaG_ae_psap_clap=NaN(size(MLO_rd.BaG0_A12_clap));
MLO_rd.BaG_ae_psap_clap(P1)=MLO_rd.Bac1_A81_ae31(P1).*cte_ae16_PSAP; 
MLO_rd.BaG_ae_psap_clap(P2)=MLO_rd.BaG0_A11_psap(P2);
MLO_rd.BaG_ae_psap_clap(P3)=MLO_rd.BaG0_A12_clap(P3);

% do the same with B and R for abs exp
MLO_rd.BaB_ae_psap_clap=NaN(size(MLO_rd.BaB0_A12_clap));
MLO_rd.BaB_ae_psap_clap(P2)=MLO_rd.BaB0_A11_psap(P2);
MLO_rd.BaB_ae_psap_clap(P3)=MLO_rd.BaB0_A12_clap(P3);

MLO_rd.BaR_ae_psap_clap=NaN(size(MLO_rd.BaR0_A12_clap));
MLO_rd.BaR_ae_psap_clap(P2)=MLO_rd.BaR0_A11_psap(P2);
MLO_rd.BaR_ae_psap_clap(P3)=MLO_rd.BaR0_A12_clap(P3);

% time series AE16+PSAP+CLAP until 2022+AE33
P4=timerange('2022-11-29',last_date);
MLO_rd.BaG_ae_psap_clap_ae33=NaN(size(MLO_rd.BaG0_A12_clap));
MLO_rd.BaG_ae_psap_clap_ae33(P1)=MLO_rd.Bac1_A81_ae31(P1).*cte_ae16_PSAP; 
MLO_rd.BaG_ae_psap_clap_ae33(P2)=MLO_rd.BaG0_A11_psap(P2);
MLO_rd.BaG_ae_psap_clap_ae33(P3)=MLO_rd.BaG0_A12_clap(P3);
MLO_rd.BaG_ae_psap_clap_ae33(P4)=MLO_rd.Ba3_A82_ae33(P4);

% time series AE16+PSAP+CLAP until 2015+AE33
P3=timerange('2013-01-01','2015-01-01');
P4=timerange('2015-01-01',last_date);
MLO_rd.BaG_ae_psap_clapS_ae33=NaN(size(MLO_rd.BaG0_A12_clap));
MLO_rd.BaG_ae_psap_clapS_ae33(P1)=MLO_rd.Bac1_A81_ae31(P1).*cte_ae16_PSAP; 
MLO_rd.BaG_ae_psap_clapS_ae33(P2)=MLO_rd.BaG0_A11_psap(P2);
MLO_rd.BaG_ae_psap_clapS_ae33(P3)=MLO_rd.BaG0_A12_clap(P3);
MLO_rd.BaG_ae_psap_clapS_ae33(P4)=MLO_rd.Ba3_A82_ae33(P4);

names_abs_tr={'BaG_ae_psap_clap', 'BaG_ae_psap_clap_ae33','BaG_ae_psap_clapS_ae33'};
break_MLO_abs=change_point_analysis_Def(MLO_rd,names_abs_tr,0.05,'MLO','abs three homogeneisation'); %[output:83915ba8] %[output:1f3e4364] %[output:1a8192d8] %[output:470a0230]
T_MLO_abs=make_table_breakpoints_def(break_MLO_abs); %[output:8457fe63] %[output:2ff15fcb] %[output:3a6b3790] %[output:8baa350f] %[output:535ffb60] %[output:4571544a] %[output:1f5daac6] %[output:360bb907] %[output:84c34891] %[output:4252e3ba] %[output:8c6d783a] %[output:992c84b3] %[output:12b10a33] %[output:645e0844] %[output:8934ee95]

%%
% Absorbion: three times series needed. but much longer
% 1) for AE16+ AE31 + AE33: 1991-2025 AE31 is too noisy and no
% homogenisation with AE33 possible --> keep AE16+PSAP+CLAP and AE33 for
% 2015-2025

%%
lambdaSC=[450;550;700];
lambdaAE=[467;530;660];
lambdaAE7=[370 470 520 590 660 880 950];

names_neph={'BsB_S1S10', 'BsG_S2S20', 'BsR_S3S30','BbsB_S1S10', 'BbsG_S2S20', 'BbsR_S3S30'};
MLO_expSC=compute_exp_D(MLO_rd,names_neph,lambdaSC);
names_abs_tr={'BaG_ae_psap_clap','Ba3_A82_ae33'};
names_expA_clap={'BaB_ae_psap_clap','BaG_ae_psap_clap','BaR_ae_psap_clap'};

names_abs_ae={ 'Ba1_A82_ae33','Ba2_A82_ae33','Ba3_A82_ae33','Ba4_A82_ae33','Ba5_A82_ae33','Ba6_A82_ae33','Ba7_A82_ae33', };

MLO_expA=compute_exp_D(MLO_rd,names_abs_ae,lambdaAE7); %[output:25317e8d] %[output:6e2714c8] %[output:1dd95874] %[output:6cd8abec] %[output:77b11528] %[output:3f6f6768] %[output:5752878b] %[output:174e4178] %[output:1c613be9] %[output:21f652b9] %[output:3839e252] %[output:5f33d8e2] %[output:2b810821] %[output:02fb05d8] %[output:4c1a1aab] %[output:1eb51541] %[output:22a12065] %[output:512e6906] %[output:122d0140] %[output:3938b5c9] %[output:8da9cac7] %[output:502a80f4] %[output:1806be07] %[output:26d12cd8] %[output:00a31fde] %[output:7fd373ec] %[output:3f00a6d1] %[output:2c216b39] %[output:6369f23b] %[output:293bef6a] %[output:557de28f] %[output:8c8cfbfe] %[output:7c248ad9] %[output:44ac16be] %[output:1a829cad] %[output:6995da25] %[output:25778c6a] %[output:07eae6fc] %[output:000df71f] %[output:96144a9f] %[output:32b90bd2] %[output:7cc8af50] %[output:82f9b87d]
MLO_expA_clap=compute_exp_D(MLO_rd,names_expA_clap,lambdaAE);
MLO_expA_clap.Properties.VariableNames{1}='expAclap_bg';
MLO_expA_clap.Properties.VariableNames{2}='expAclap_br';
MLO_expA_clap.Properties.VariableNames{3}='expAclap_gr';

MLO_exp=synchronize(MLO_expSC,MLO_expA);
MLO_exp=synchronize(MLO_exp,MLO_expA_clap);
MLO_SSA1=compute_SSA_D(MLO_rd,{'BsG_S2S20'},{'BaG_ae_psap_clap'});
MLO_SSA2=compute_SSA_D(MLO_rd,{'BsG_S2S20'},{'Ba3_A82_ae33'});
MLO_SSA=synchronize(MLO_SSA1,MLO_SSA2);
MLO_cal=synchronize(MLO_exp,MLO_SSA);
clear MLO_expSC MLO_expA MLO_exp MLO_SSA MLO_SSA1 MLO_SSA2;
plotFigControl_cal(MLO_cal, MLO_st.name); %[output:12105f38] %[output:9e6249f2] %[output:8a90d031]
names_cal=fieldnames(MLO_cal);

% % MLO_cal=compute_exp_SSA(MLO_rd,lambdaSC, lambdaAE3);
% % plotFigControl_cal(MLO_cal, MLO_st.name);
%%
% statistical BP for computed value
%exp
names_exp={'expS_bg', 'expA_fit'};
break_MLO_exp=change_point_analysis_Def(MLO_cal,names_exp,0.05,'MLO','Exp S and A'); %[output:60440beb] %[output:1fbff167] %[output:83cb2646] %[output:1604147e]
T_MLO_exp=make_table_breakpoints_def(break_MLO_exp); %[output:40267d19] %[output:300021bf] %[output:76af1ee0] %[output:2806bc31] %[output:669feba5] %[output:368919af] %[output:7575622c] %[output:8e498ed9]

% BbsF
names_BbsF={'BbsFb', 'BbsFg'};
break_MLO_BbsF=change_point_analysis_Def(MLO_cal,names_BbsF,0.05,'MLO','BbsF');
T_MLO_exp=make_table_breakpoints_def(break_MLO_BbsF);

% SSA
names_SSA={'SSA_MLO_SSA1', 'SSA_MLO_SSA2'};
break_MLO_SSA=change_point_analysis_Def(MLO_cal,names_SSA,0.05,'MLO','SSA');
T_MLO_SSA=make_table_breakpoints_def(break_MLO_SSA);
%%
%Questions:
%problems:
 %scat: what is BsQ_S ?
 %scat minima are lower from 1985 to1995, reason ?
 %scat minima and maxima are lower before 1985
 %scat: strange high data in december 1993-jan 94
 %ratio of scattering is usually lower in 2000 than after
 % all the backscat ratio between several wavelength have less small values and negatives between 2001 and 2006, why ?
 % exp scat seems also to be lower in 2000 than after. More pronounced for PM10 than PM 1. Visible also on the BbsF
 
 %PSAP 1w to 3w change in 2006, the max seems to be lower with 1w.
 % exp =0 and exp= 3.16 seems to be favored
 
 
 %N: ruP_mri_mseure in seP_mri_mse 1991, not homogeneous before 1980, strange in 2011
%%
MLO_tr=MLO_rd;
MLO_tr.y=year(MLO_tr.Time);
% begin at the beginning of a year:1988 for scat, 2001 for backscat and
% 1992 for abs
%end: 2022/2025
P=timerange('1988-01-01','2026-01-01');
MLO_tr=MLO_tr(P,:);

%%% HAS TO BE CHANGE IF PM1 HAS TO BE USED
%no RH trend and no dry trend

%remove all not necessary variables
names=fieldnames(MLO_tr);
c= startsWith(names,["T1";"T0";"T_";"P_";"P1_";"P0_";"N";"U_";"U1";"U0";"X"]) |  contains(names,["Q","dry","mri","mse","tsi","A11_psap","A12_clap","ae31"]);
N=names(c);
for i=1:length(N)
    MLO_tr.(N{i})=[];
end

% begin Bbs in 2001  and R should not be used
names=fieldnames(MLO_tr);
P2=timerange('1988-01-01','2001-01-01');
c=startsWith(names,"Bbs");
N=names(c);
for i=1:length(N)
    MLO_tr.(N{i})(P2)=NaN;
end 

%begin AE33 in 1.1.2015
P5=timerange('1988-01-01','2015-01-01');
N=names(contains(names,"ae33"));
for i=1:length(N)
    MLO_tr.(N{i})(P5)=NaN;
end 


% expS broadens with instrument's change: OK but interpretation should be made with great caution

%old stuff 2020
% % %Exp Abs since 2007:
% % P3=timerange('1998-01-01','2007-01-01');
% % c=startsWith(names,["BaB";"BaR"]);
% % N=names(c);
% % for i=1:length(N)
% %     MLO_tr.(N{i})(P3)=NaN;
% % end 

% invalidate scat 17 October 1993 to 27 April 1994 (no Bbs at that time)
% yes, has to be done, stat BP
P4=timerange('1993-10-17','1994-04-27');
c=startsWith(names,"Bs");
N=names(c);
for i=1:length(N)
    MLO_tr.(N{i})(P4)=NaN;
end 

% Use the homogenized time series BsB_S1S10


%%MLO_cal_tr=compute_exp_SSA(MLO_tr,lambdaSC, lambdaAE);

%backscattering fraction and stop in Oct 2022 (higher thereafter)

%compute SSA with green from the AE16+PSAP+CLAP time series + from AE33 at
%520 nm
names_neph={'BsB_S1S10', 'BsG_S2S20', 'BsR_S3S30','BbsB_S1S10', 'BbsG_S2S20', 'BbsR_S3S30'};
MLO_expSC_tr=compute_exp_D(MLO_tr,names_neph,lambdaSC);

MLO_expA_tr=compute_exp_D(MLO_tr,names_abs_ae,lambdaAE7); %[output:7daedd72] %[output:23a2f8e6] %[output:6757de8d] %[output:297a0690] %[output:6ff3f15a] %[output:77a8cec1] %[output:803ba2bd] %[output:4225465a] %[output:43feaaec] %[output:841cf8f7] %[output:713e3948] %[output:6fce8c1e] %[output:015bbced] %[output:38dc8002] %[output:1251c42f] %[output:7fa8891f] %[output:30485f3c] %[output:65c6d840] %[output:155410f8] %[output:4570c4f6] %[output:919c996f] %[output:5bf15f4b] %[output:67838364] %[output:1f67c1c0] %[output:601b34ec] %[output:1757681e] %[output:6f129a6f] %[output:883f1fd6] %[output:78e033c7] %[output:59287da2] %[output:1ce6a2da] %[output:0b37b71a] %[output:4e5b1aed] %[output:83e132b3] %[output:816ce839] %[output:08af833b] %[output:45825452] %[output:71eb9b3e] %[output:979a5aa2] %[output:254a37f7] %[output:80e52097] %[output:5994e055] %[output:7f107dd8]
MLO_expA_clap_tr=compute_exp_D(MLO_tr,names_expA_clap,lambdaAE);
MLO_expA_clap_tr.Properties.VariableNames{1}='expAclap_bg';
MLO_expA_clap_tr.Properties.VariableNames{2}='expAclap_br';
MLO_expA_clap_tr.Properties.VariableNames{3}='expAclap_gr';

MLO_exp_tr=synchronize(MLO_expSC_tr,MLO_expA_tr);
MLO_cal_tr=synchronize(MLO_exp_tr,MLO_expA_clap_tr);

%MLO_cal_tr=synchronize(MLO_exp_tr,MLO_expA_tr);
MLO_cal_tr.SSA0G=MLO_tr.BsG_S2S20./(MLO_tr.BsG_S2S20+MLO_tr.BaG_ae_psap_clap);
%if SSAAE: MLO_cal2.SSA0B=MLO_tr.BsB_S1S10./(MLO_tr.BsB_S1S10+MLO_tr.BaG_ae_psap_clap);
MLO_cal_tr.SSA0AE=MLO_tr.BsG_S2S20./(MLO_tr.BsG_S2S20+MLO_tr.Ba3_A82_ae33); % if SSAAE take 470nm

%expS  trend on BG and expA fit from AE33 only:
MLO_cal_tr.expS_br=[]; %[output:3340cb3a]
MLO_cal_tr.expS_gr=[];
% MLO_cal2.expS_br1=[];
% MLO_cal2.expS_gr1=[];
%%

MLO_cal_tr.expA_bgAE=[];
MLO_cal_tr.expA_grAE=[];
MLO_cal_tr.expA_brAE=[];
MLO_cal_tr.expAclap_gr=[]; %[output:585aeada]
MLO_cal_tr.expAclap_br=[];
%backscattering fraction has to be invalidated after the erruption
P4=timerange('2022-11-29','2026-01-01');
MLO_cal_tr.BbsFb(P4)=NaN;
MLO_cal_tr.BbsFg(P4)=NaN;
MLO_cal_tr.BbsFr(P4)=NaN;


% MLO_cal2.expA_br1=[];
% MLO_cal2.expA_gr1=[];

MLO_tr=outerjoin(MLO_tr, MLO_cal_tr);

%trend calculated only on G for neph + abs (R should not be used)
MLO_tr.BsB_S1S10=[];
MLO_tr.BsR_S3S30=[];
MLO_tr.BbsB_S1S10=[];
MLO_tr.BbsR_S3S30=[];
MLO_tr.BbsFb=[];
MLO_tr.BbsFr=[];

MLO_tr.Ba1_A82_ae33=[];
MLO_tr.Ba2_A82_ae33=[];
MLO_tr.Ba4_A82_ae33=[];
MLO_tr.Ba5_A82_ae33=[];
MLO_tr.Ba6_A82_ae33=[];
MLO_tr.Ba7_A82_ae33=[];

MLO_tr.BaB_ae_psap_clap=[];
MLO_tr.BaR_ae_psap_clap=[];

MLO_tr.BaG_ae_psap_clapS_ae33=[];
MLO_tr.BaG_ae_psap_clap_ae33=[];
clear MLO_expA_tr MLO_expSC_tr MLO_cal2 MLO_rd_old;
%%
% IF TREND HAVE TO BE COMPUTE AGAIN REMOVE EXPS BEFORE 2001 + REMOVE
% CLAP==0 + COMPUTE SSA EXPA
[MLO_result_MK,MLO_result_LMSlog,MLO_result_LMSlin]=all_trend_STN(MLO_tr,MLO_st); %[output:7daca7d3] %[output:691d70ae] %[output:56b72ec6] %[output:4cbd50fa] %[output:8767859a] %[output:11b066ec] %[output:009c4dd0] %[output:29858109] %[output:8c609ceb] %[output:43878f87] %[output:92137970] %[output:0a50c7ad] %[output:9533c316] %[output:96569083] %[output:3f72c21c] %[output:8a64f282] %[output:819da75a] %[output:5ad1ea29] %[output:0ac26fcc] %[output:21f31491] %[output:47f7eba9] %[output:184e0113] %[output:8b5c357c] %[output:9dfc1fb3] %[output:72c69771] %[output:92bbcedd] %[output:40aa086e] %[output:18c8d889] %[output:53b4c1d4] %[output:77ee198a] %[output:6577f61e] %[output:21cf559b] %[output:6ba35594] %[output:914fd866] %[output:4a8a6645] %[output:94b07b84] %[output:16cee0e4] %[output:07b88436] %[output:94d83b3f] %[output:7261d18d] %[output:6c4e33ba] %[output:44f1c894] %[output:43e6e121] %[output:5ca5a6ff] %[output:32a14f56] %[output:71e76d47] %[output:7ab5e6f4] %[output:34001b96] %[output:5f3e610b] %[output:106c22b1] %[output:782857ca] %[output:41d65df9] %[output:1891fe05] %[output:38c80f41] %[output:2b079263] %[output:0874b400] %[output:881ddcf0] %[output:44031ae0] %[output:34396911] %[output:27c4bfa6] %[output:4a0ad7c1] %[output:65e8ab02] %[output:6d460f76] %[output:164e8471] %[output:50dc1af8] %[output:5f574c25] %[output:450a624e] %[output:74e66a45] %[output:6a5fabc0] %[output:7675f741] %[output:24aa0d05] %[output:57521f4e] %[output:435b04ff] %[output:5834e47b] %[output:51c0283c] %[output:03ea7eb2] %[output:8cb754ed] %[output:2521d800] %[output:96422d9e] %[output:6b93ccef] %[output:1ff1e6a5] %[output:5a750aa8] %[output:3a652dd5] %[output:6a45a07e] %[output:440cb725] %[output:03527de8] %[output:62b980ac] %[output:145cb312] %[output:1fb44726] %[output:9a72e9db] %[output:970d9f46] %[output:627b9ef2] %[output:7bb8a256] %[output:047d3508] %[output:53258313] %[output:00eb0146] %[output:5f73b0d6] %[output:61911cd4] %[output:9d41c5c9] %[output:0e4470d9] %[output:4e2b6092] %[output:6f127a96] %[output:4b83007b] %[output:0484a902] %[output:11101af8] %[output:37495a31] %[output:3c309555] %[output:74ccd857] %[output:14be6490] %[output:65bbdc35] %[output:196cb9c7] %[output:6e09b832] %[output:51081c51] %[output:5d738641] %[output:8f359d5a] %[output:4481adb1] %[output:40ddc0bb] %[output:547dbc7b] %[output:3725f127] %[output:75f627c0] %[output:4ccb3eb6] %[output:5df4caf1] %[output:03d4c46f] %[output:5d7ad6d2] %[output:27351042] %[output:4a092aea] %[output:1540ae5a] %[output:81dcf9e7] %[output:03814767] %[output:3c0856bc] %[output:644bdff9] %[output:49c60f87] %[output:08f093fe] %[output:6afaeda3] %[output:1537097e] %[output:50a0bede] %[output:3bd90208] %[output:6bb0d786] %[output:4925dd53] %[output:7f550899] %[output:77943ab5] %[output:5744ba15] %[output:8431292a] %[output:1c856963] %[output:5c9fe296] %[output:11153c03] %[output:2693967a] %[output:39a6b9c9] %[output:64f431c9] %[output:72d925c7] %[output:9e3df2f5] %[output:173985e8] %[output:2c0fdeaf] %[output:3ce2d07c] %[output:21497ac0] %[output:02bba6a4] %[output:5ee0f6eb] %[output:1f63a72f] %[output:96ec5789] %[output:9a056215] %[output:7cca7a2c] %[output:4549027b] %[output:4115420d] %[output:8d4861a9] %[output:0931055f] %[output:350777ed] %[output:1a4e8313] %[output:1ee5317e] %[output:84db3364] %[output:10d3ba1e] %[output:61806ca5] %[output:63e1ac1d] %[output:7aac0809] %[output:2fba3f96] %[output:2dd0c217] %[output:0a526f37] %[output:1e8ae956] %[output:1b27364e] %[output:6ba2f081] %[output:3218c0d2] %[output:42a6b02f] %[output:05cb2d68] %[output:29988949] %[output:0ce840fe] %[output:260a6be2] %[output:70604c23] %[output:5eda84ef] %[output:5135b724] %[output:69bd9765] %[output:87e0606b] %[output:8a8b83fc] %[output:9683eb8d] %[output:59ed2d67] %[output:074657ca] %[output:5763ca47] %[output:3285274e] %[output:01434088] %[output:0daecedc] %[output:4fdfe614] %[output:1bd8ff71] %[output:15ced2fa] %[output:990f2ee7] %[output:3937d464] %[output:029880fb] %[output:62b5eea5] %[output:2fcb0a8f] %[output:8712811b] %[output:1ceb2ef3] %[output:3f3ea777] %[output:175f3b96] %[output:93ea4fbc] %[output:8c1d9f74] %[output:2334239c] %[output:37456ea4] %[output:06cd1b08] %[output:8fa8128b] %[output:3241aa6e] %[output:36208924] %[output:97346451] %[output:63a05afe] %[output:9453c2ec] %[output:9362a015] %[output:18fb473c] %[output:1c95b9b9] %[output:0bcd305f] %[output:4b002ed1] %[output:7a09ab26] %[output:206afafe] %[output:1f8b5d75] %[output:5648a775] %[output:35eea1c6] %[output:4b8031c0] %[output:6863686a] %[output:81d0f071] %[output:92100cbc] %[output:0acc5e8e] %[output:31975bb4] %[output:230e981d] %[output:6f869c7e] %[output:7a886ec3] %[output:0ef3773d] %[output:2ab8c7fa] %[output:676164f7] %[output:455d64b0] %[output:2f3064ba] %[output:82c04625] %[output:1633a5b9] %[output:00329b5d] %[output:4e26f982] %[output:123b00ff] %[output:5c17d0a0] %[output:3742a577] %[output:93505a75] %[output:61120f9c] %[output:6a5714b3] %[output:33cd7afd] %[output:1884a694] %[output:1db9a7de] %[output:57f1e15f] %[output:9580b661] %[output:8e0e0ff8] %[output:52ac184d] %[output:93cb677e] %[output:40d86b5e] %[output:8f043f03] %[output:9ec8596d] %[output:17b5ad5e] %[output:876a8eb7] %[output:9d3b11fe] %[output:2551645f] %[output:3e462a53] %[output:4603b234] %[output:2cdbab0d] %[output:6ba43a7a] %[output:292005b6] %[output:44587408] %[output:2d2c25dc] %[output:9a4e322b] %[output:7209c994] %[output:4aaaea36] %[output:15765720] %[output:3cccc1c3] %[output:36b97662] %[output:47d109dd] %[output:6d6aed90] %[output:719ba86f] %[output:60041ab7] %[output:8288b7ce] %[output:47a2a2cf] %[output:2fe9d7c4] %[output:3e1bed74] %[output:6fd93c27] %[output:4e7457ff] %[output:458e336e] %[output:0cb08232] %[output:553ea746] %[output:1c44b927] %[output:1c30507b] %[output:759a6b80] %[output:5c9334e8] %[output:0587c776] %[output:3971c5ce] %[output:04efeab9] %[output:013402e7] %[output:967fbaed] %[output:3116fd48] %[output:906525eb] %[output:4b070030] %[output:89ad159d] %[output:1ad2fd08] %[output:750a3223] %[output:57d51fee] %[output:1b4f49ad] %[output:19cd64b9] %[output:6032e7cb] %[output:32c1fec6] %[output:887ba1af] %[output:472fba2e] %[output:7ed1c1f8] %[output:3e50198a] %[output:2d5d30b0] %[output:783e6de4] %[output:7a094f46] %[output:1ec23811] %[output:7050df14] %[output:7da68526] %[output:80d45bc7] %[output:66e98f44] %[output:691c68b4] %[output:1847470b] %[output:5b01f108] %[output:77ae7dc0] %[output:763e72c7] %[output:41e0e013] %[output:4767e172] %[output:091c8eab] %[output:0375b2ad] %[output:308a38e3] %[output:57970a44] %[output:0a9eaadd] %[output:83eefd76] %[output:243b56a7] %[output:7d25fdd6] %[output:62224dce] %[output:71005745] %[output:83d93634] %[output:9c7cc22d] %[output:64150d21] %[output:728da710] %[output:1190d817] %[output:381adbb8] %[output:85b092d7] %[output:1ae00d2a] %[output:9b025613] %[output:46dc88b7] %[output:8dec3a57] %[output:952f7d45] %[output:11939650] %[output:11838473] %[output:187cb860] %[output:4f7ac2c3] %[output:1d78b69b] %[output:0f71f93d] %[output:5e8f9139] %[output:98ce6d35] %[output:2ef1c20e] %[output:76107dbb] %[output:07c20a40] %[output:1b7752e3] %[output:730d9a8b] %[output:1f45c1a0] %[output:8a727372] %[output:62e057cd] %[output:378fafb1] %[output:45651540] %[output:18518a12] %[output:500ad4dc] %[output:26dd33b3] %[output:57de100a] %[output:4aaa3564] %[output:55bda12e] %[output:0bc04a08] %[output:90cfd8ca] %[output:5e1643bf] %[output:3ff75fa2] %[output:8036573e] %[output:385019d6] %[output:0c1aa960] %[output:4ba62846] %[output:6937e54a] %[output:8174113e] %[output:3fd50b0f] %[output:5413a336] %[output:3fbb82da] %[output:4f138eaf] %[output:26350180] %[output:5cf78d39] %[output:2297ef76] %[output:4db7e9e6] %[output:6216ff5b] %[output:3372c653] %[output:2684268e] %[output:58083934] %[output:2698c44d] %[output:1995cd3e] %[output:01fdb1b3] %[output:57ac2e2e] %[output:76244fe9] %[output:87ecc18b] %[output:35966a80] %[output:953d2ac2] %[output:31a1500c] %[output:99b1044f] %[output:019ba41b] %[output:6830611f] %[output:4cc68b26] %[output:0c6414aa] %[output:24a80196] %[output:67a5241c] %[output:37c42036] %[output:2d76c683] %[output:398012ae] %[output:1c2f3323] %[output:05353ff7]
% make only one time series with homo1 and AE33 data (abs, expA and SSA)
% not sure it is necessary

writetable(MLO_result_MK,'MLO_res_MK_10y.txt'); %, 'delimiter',',' )
writetable(MLO_result_LMSlog,'MLO_res_LMSlog.txt'); 
writetable(MLO_result_LMSlin,'MLO_res_LMSlin.txt'); 
plot_10y_in_two(MLO_result_MK, MLO_st,'y'); %[output:9a1762cd] %[output:24228ba8]

%%
% compute all trend only until 2022 for coherence for map
P10=timerange('2023-01-01','2026-01-01');
MLO_tr22=MLO_tr;
MLO_tr22(P10,:)=NaN; %[output:3d619543]
%%

[MLO_result_MK22,MLO_result_LMSlog22,MLO_result_LMSlin22]=all_trend_STN(MLO_tr,MLO_st); %[output:59e5c07c] %[output:98f1737e] %[output:6a259ee5] %[output:30107ed0] %[output:86913e46] %[output:86ae1794] %[output:0d65490a] %[output:1fb13ebd] %[output:17813802] %[output:62813a08] %[output:57d58dda] %[output:7dbfbb5b] %[output:281b8db5] %[output:9703e02a] %[output:75e1381e] %[output:433c385b] %[output:9012bd15] %[output:432e6906] %[output:98bdac5e] %[output:91e26909] %[output:3d068893] %[output:89132f54] %[output:20dd97a2] %[output:59bbc593] %[output:2c998b4e] %[output:61384b29] %[output:676b78ab] %[output:1de67dfe] %[output:1555129b] %[output:63316fac] %[output:6852646d] %[output:49ad3ea2] %[output:31881498] %[output:89c25b2e] %[output:431a590d] %[output:3d4f5acc] %[output:44b6fe82] %[output:8c1b98ab] %[output:0dceedde] %[output:9dafe738] %[output:4d5d3991] %[output:0b49fca7] %[output:57e28deb] %[output:229a40b4] %[output:9f66856c] %[output:9b56944a] %[output:98c72be3] %[output:8b70f32f] %[output:22372109] %[output:97dfb77a] %[output:93ed799a] %[output:9a30654e] %[output:5ad853e2] %[output:94afcfa6] %[output:3cb0572b] %[output:10869dee] %[output:511a9ba0] %[output:11c2cc34] %[output:28efcf19] %[output:8c9fc4d7] %[output:7d6adce1] %[output:92d1bc50] %[output:771ddec4] %[output:48187208] %[output:5b471001] %[output:55e80d49] %[output:2b9f215a] %[output:9adc4cad] %[output:355379ca] %[output:3a82d0c1] %[output:8373293f] %[output:3d0d2235] %[output:67ac51ca] %[output:9cfe8fa2] %[output:22d5e231] %[output:302840ca] %[output:550ea67e] %[output:372051f5] %[output:083af83f] %[output:860b4267] %[output:6b7bcdb1] %[output:2d88d7fc] %[output:515cd56f] %[output:37f4ed21] %[output:4b21e9b5] %[output:1b7fe61b] %[output:310232db] %[output:6487bf70] %[output:3f4d62cd] %[output:66532130] %[output:4c5523b3] %[output:06d28da8] %[output:00135a40] %[output:9a07706f] %[output:858bea75] %[output:628d544e] %[output:00aecbf3] %[output:9610ba5f] %[output:82fac96f] %[output:12419720] %[output:0e3b2499] %[output:87bf7e4b] %[output:3a16d8ca] %[output:77e60e6f] %[output:65915a7b] %[output:6bacbedd] %[output:93276de6] %[output:6db1c725] %[output:7ad67874] %[output:69b7c3e9] %[output:190530e9] %[output:46ed0118] %[output:370377e0] %[output:309bfdfb] %[output:3f1bb495] %[output:55542bd7] %[output:15a6f9f1] %[output:44b2e047] %[output:1359598b] %[output:8ce27118] %[output:9dbdb39d] %[output:31f421e2] %[output:35fc1edd] %[output:5558999d] %[output:4d1a220d] %[output:40ab9bdc] %[output:4f465fd1] %[output:394eb751] %[output:18a11ff6] %[output:5ae7b219] %[output:6dee4860] %[output:6f88328c] %[output:5a2be3a6] %[output:094299d1] %[output:68ff7093] %[output:94f290f7] %[output:892e4778] %[output:8ed20ad6] %[output:6ac1ba10] %[output:1f3cd523] %[output:4d1772be] %[output:583b9292] %[output:2ca426da] %[output:059a2990] %[output:2fd1e05c] %[output:4c1cef6b] %[output:45f57394] %[output:7ef74ba3] %[output:3358b62c] %[output:0ae26367] %[output:51c88e2c] %[output:04b6eb43] %[output:4d249b88] %[output:85986063] %[output:3f90aed5] %[output:346d758c] %[output:80b93dc6] %[output:26b2ce6c] %[output:67d1ee70] %[output:5055efa3] %[output:8633091b] %[output:28754ac1] %[output:66b79876] %[output:71e84bf7] %[output:463028d2] %[output:8a02e87f] %[output:71f84cf7] %[output:12724ccc] %[output:214e3706] %[output:3c89837d] %[output:3d516f1d] %[output:5f772d61] %[output:130f004f] %[output:6de7d553] %[output:400a7d65] %[output:4f965591] %[output:148c6383] %[output:76af6f19] %[output:77849cf2] %[output:757a5e2f] %[output:06978dce] %[output:3a4db530] %[output:2362cc63] %[output:41595d72] %[output:51b040a8] %[output:5e3c6c8f] %[output:327ea8ad] %[output:8e21ce04] %[output:3609f707] %[output:0a1dedf0] %[output:97c71e2b] %[output:31c0cb8e] %[output:69a4bc62] %[output:8c6a82ec] %[output:33cb0187] %[output:41263c8d] %[output:2e7dbdc0] %[output:47edae40] %[output:1f053cf7] %[output:95a971e1] %[output:334a2840] %[output:49204063] %[output:6804cc6d] %[output:75d5ffed] %[output:95a232f0] %[output:15857ebd] %[output:6be91978] %[output:8c694714] %[output:855e2d48] %[output:824a3d00] %[output:83d45c79] %[output:96a1fc1e] %[output:9db438f5] %[output:98fa302a] %[output:48031cf5] %[output:2261c2dd] %[output:95723cc6] %[output:295e44d9] %[output:42f76616] %[output:90e4dc8e] %[output:67858de1] %[output:32b5fb61] %[output:3eaa9cc1] %[output:635de8f2] %[output:0857736b] %[output:9293f623] %[output:12f37b2a] %[output:829ac9f7] %[output:71cf94ce] %[output:4703aac0] %[output:20edc517] %[output:36a2bd34] %[output:370831e2] %[output:19f3464a] %[output:2c29ba4b] %[output:1929e0bb] %[output:8c906ac7] %[output:43a60f55] %[output:39c35e9d] %[output:47cea7b7] %[output:7ca87854] %[output:2668f26a] %[output:7edeab23] %[output:0140b7a6] %[output:13fffbdc] %[output:01c2ea91] %[output:81f8a675] %[output:3fcd517c] %[output:41702c92] %[output:3292f211] %[output:1168d03b] %[output:72016153] %[output:530ae009] %[output:03355161] %[output:45adebb6] %[output:2104104e] %[output:5fa3211b] %[output:6b15c5ea] %[output:2d16e8a1] %[output:02e090a3] %[output:456c9889] %[output:96d5c9d8] %[output:0d93a46b] %[output:6aa4ba47] %[output:2dfe4951] %[output:2e506fda] %[output:033dfb7d] %[output:929df44a] %[output:0e38ca1a] %[output:85d6be29] %[output:12249d04] %[output:6b9ff3a1] %[output:2c9e9183] %[output:3885dcb7] %[output:73abecde] %[output:3f45b41b] %[output:69e60037] %[output:9e8543d6] %[output:6f6eb02f] %[output:25ce1c9c] %[output:8b5febd5] %[output:5033a840] %[output:3187f9a7] %[output:7708404f] %[output:710924d6] %[output:79e4fc6f] %[output:98dbb739] %[output:51447b7c] %[output:9655d290] %[output:135f7590] %[output:84d006e3] %[output:1c7575f9] %[output:04f69dfb] %[output:1f2ebffe] %[output:34eb889f] %[output:252917bf] %[output:9a7635e6] %[output:0721f5a4] %[output:345642f3] %[output:32a01465] %[output:0fc3f6c8] %[output:99312fe7] %[output:6c78418b] %[output:52d68318] %[output:0acd2952] %[output:7ee41ed6] %[output:62789c26] %[output:15501e2e] %[output:6949bef0] %[output:5f75a1eb] %[output:103734d1] %[output:147192cb] %[output:8dffcf54] %[output:56e3d87e] %[output:4b660e78] %[output:7df7f2ad] %[output:1da6e4a6] %[output:3bc478d1] %[output:8d00033c] %[output:884fc1cd] %[output:3627716b] %[output:8c5ffa0c] %[output:96e7b1b0] %[output:97fff6f7] %[output:4629151d] %[output:060eef26] %[output:87210ed5] %[output:4b014e8f] %[output:11524587] %[output:2cd24c98] %[output:0d156e82] %[output:328b7726] %[output:1da066bd] %[output:87f7d9b4] %[output:383cfe7c] %[output:6c2778f5] %[output:186116fe] %[output:24aa8ffb] %[output:9a0ecb32] %[output:5fda4f5e] %[output:22fbca55] %[output:60595b92] %[output:5921b86b] %[output:87e484a8] %[output:236ed6c2] %[output:8dafb288] %[output:8b934d7d] %[output:51e6f097] %[output:5d60b099] %[output:229611ec] %[output:54f1c4b7] %[output:4d07ab76] %[output:1fd01ed7] %[output:2d60b6cc] %[output:419ff267] %[output:2b07e8b6] %[output:24acbd94] %[output:2f07897e] %[output:104a27a1] %[output:08595ca8] %[output:0a7d86dc] %[output:4e7f48be] %[output:98e67a13] %[output:162edff5] %[output:7111247d] %[output:3f52226f] %[output:0cd432ee] %[output:3dd51346] %[output:792e1ccf] %[output:0291b4f0] %[output:665930fa] %[output:8e9e22bb] %[output:0aa4eb7c] %[output:10d457b6] %[output:9be7d1ca] %[output:2bd576fb] %[output:4bf0237d] %[output:266a4e2a] %[output:59768bd8] %[output:21d1bb56] %[output:3cc784db] %[output:9e7aa3d7] %[output:071e8ea0] %[output:8eb9d1ec] %[output:57b48c07] %[output:33fe7673] %[output:1a5f3c19] %[output:9ab16a12] %[output:17559201] %[output:97a9ea44] %[output:30ffd511] %[output:37ab7f95] %[output:38122745] %[output:131aff81] %[output:9b8c9542] %[output:2079614f] %[output:8cbf4444] %[output:682df635] %[output:2b62cd3e] %[output:74ea0433] %[output:450fa3e5] %[output:054a135e] %[output:8c7f84d9] %[output:05aa815e] %[output:1f8c873b] %[output:62cf130a] %[output:8dbac247] %[output:2c5c8e72] %[output:889e801b] %[output:73e99599] %[output:25befab1] %[output:70dab343] %[output:2076871c] %[output:94f630f5] %[output:3cfdb78e] %[output:59bcc8e3] %[output:236fa145] %[output:395f53f9] %[output:31ba1dda] %[output:22917fba] %[output:1528eb93] %[output:55a00260] %[output:3afcb33e] %[output:435725cb]
% make only one time series with homo1 and AE33 data (abs, expA and SSA)
% not sure it is necessary

writetable(MLO_result_MK22,'MLO_res_MK_present.txt'); %, 'delimiter',',' )
writetable(MLO_result_LMSlog22,'MLO_res_LMSlog_present.txt'); 
writetable(MLO_result_LMSlin22,'MLO_res_LMSlin_present.txt'); 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":31.4}
%---
%[output:4600e4b1]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"first_date","value":"   1974-01-01 00:00:00\n"}}
%---
%[output:3ccf89c6]
%   data: {"dataType":"textualVariable","outputData":{"header":"datetime","name":"last_date","value":"   2025-12-31 00:00:00\n"}}
%---
%[output:35a82c92]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsB_S1S10\nTaille de la serie    : 303\nStatistique T_max     : 17.1840\np-valeur (bootstrap)  : 0.0070\nPoint de rupture      : 2024-05-01 00:00:00  (indice 284)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsB_S1S10\nTaille de la serie    : 283\nStatistique T_max     : 16.2242\np-valeur (bootstrap)  : 0.0060\nPoint de rupture      : 2000-05-01 00:00:00  (indice 2)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsB_S1S10\nTaille de la serie    : 19\nStatistique T_max     : 1.4900\np-valeur (bootstrap)  : 0.9330\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsB_S1S10\nTaille de la serie    : 281\nStatistique T_max     : 16.4707\np-valeur (bootstrap)  : 0.0060\nPoint de rupture      : 2018-08-01 00:00:00  (indice 217)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsB_S1S10\nTaille de la serie    : 216\nStatistique T_max     : 7.0077\np-valeur (bootstrap)  : 0.3070\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsB_S1S10\nTaille de la serie    : 64\nStatistique T_max     : 14.0574\np-valeur (bootstrap)  : 0.0100\nPoint de rupture      : 2022-10-01 00:00:00  (indice 50)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsB_S1S10\nTaille de la serie    : 49\nStatistique T_max     : 5.5267\np-valeur (bootstrap)  : 0.3410\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsB_S1S10\nTaille de la serie    : 14\nStatistique T_max     : 7.6104\np-valeur (bootstrap)  : 0.0970\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:774c0d08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:2215ae78]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsG_S2S20\nTaille de la serie    : 303\nStatistique T_max     : 17.5789\np-valeur (bootstrap)  : 0.0030\nPoint de rupture      : 2000-05-01 00:00:00  (indice 2)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsG_S2S20\nTaille de la serie    : 301\nStatistique T_max     : 7.0647\np-valeur (bootstrap)  : 0.3180\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsR_S3S30\nTaille de la serie    : 303\nStatistique T_max     : 15.0520\np-valeur (bootstrap)  : 0.0140\nPoint de rupture      : 2000-05-01 00:00:00  (indice 2)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsR_S3S30\nTaille de la serie    : 301\nStatistique T_max     : 7.8805\np-valeur (bootstrap)  : 0.2320\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:151eb4ae]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl8FEX2\/gKDCciNIDGBTBAUUXe9QUBJkENcDxQRBDQJiALeIngAG7IaFREFcf\/oAgaUy1Vc13MVMOHwYFd3URQQ0EwgMRwKCgjBTMg\/X81Up6bTM90z6VyTKX+RZKaq+vXrqq9fvTOqtLS0FJEW4UCEAxEOhMCBqAiAhMC1yJBKcaDkQAH2zBqG9hNfQ8PWcZWaKzK4ZjkQAZCa5X+9vfqeZ4cjukM3tBr653rLg3C48QiAhMNTrGP3QAmkYEpPuH\/aVYFyxykdEZf5aUQyqSPPNAIgdeRBRciMcKA2ciACILXxqYQ5Tf50IBHdSN178BEAqXvPrM5SHOjoIm8q+oxeiJ3yPho0aV5n77M+Ee4XQCa+vk3wYdbQruLfvYeP49q\/bkLapadhXJ8O9YlHkXu1mQMRScNmhtbgdIYAogcPSd9vv7sxevG3SD6jdQREavChRS4d4UBt4UAFAKGkcc+KbXh+eFec2iy6Ap2bfzyCSW98h8Vp5xh+X1tuLEJH7ebAwdf\/gqObPhLHlaNfr8K+WTeiQZMWiJ3yIaLP6F67iY9Qp3EgaAAxA5hw4e3nn3+OESNGiNt58sknMWzYMPH7a6+9hkceeUT8vmzZMrRp0wZpaWm4++67tT56Hqhj+N0111yDp556Co0bNw4XdgV1H+oRhgN\/fOxKtBu\/QMyxb95tOG3avyJm3KA4WnOd\/R5hurQ72fCY8uLa3cjefgAvp5yNk09y1BzlVXxlFUDUDT9jxgy89NJLlgGE\/d955x1kZWWhS5cuOHbsGB5++GHs2rULCxcuROvWrav4Tmrf9AQQCRrFB\/Lx6zvPCUmkOH9rBEBq3+MKSJEhgPCYMnjeJjHwrfHn4dzTmsLoszp2r0GRKwHksssuw6FDh8RmZxszZgyaN2+O9evXm0ogO3bsMJRO\/H0eFIF1vPORz1f6HFsIJDzGtJv4Bpr2GFLH767+kB\/QjEtp49F\/7tS48cR1neuN8lQCCI8xOTk5eOaZZwQfHnzwQQEijz32mCmAyKMLjzo9evTQ+HjgwAExBz976KGH6s9qi9xp2HEg4gfi55FKAJk2bRo2bdqE0047DU6nEytWrMDYsWNx1113WQKQuXPnascXeSl5jOHf9VEXEnEkCx8ciQCICYBQgepyubBtm8cvpmvXrujTp49QsJopUatbAnG5APmTlwckJABOJ5CUVDsWbMSRrHY8BzupCKgDGfyHU4QjmdR\/HPytGF3aNcHbd54X9iZcKYEQQBISEnwsMvJvMwCpLh0IQWPxYmD6dOOlQRBJTQXS0+1cOqHPFXEkC513tW2kqRVGOo+N7R2Pfme1AfUiO\/b9pnmo1rYbsoseFUAuuOACoQxlozXl559\/tiSBsL+02kg9iNR\/8Ds7rDAEj+Rkj+TBJiUO\/su2aJHvd1lZtUcisetZReapOQ6Y+oHoHcfqmx8IJZBrr73Wx\/S6ffv2CgDy448\/+jzFO+64Q1OQVpUfCEEjMbEcOPyBA\/sR\/3JyPH2zs2seRE4cPYTCzKtwfPsnGt8icTA1BwShXtkUQPR+H\/UFQEJlaHWNUyUPHk8IHmaNkgpBhNJJbq5Z76r7XoJHwzbxaP\/ACu1CTDJU8nN+JJiu6lhv+8wBjzDXX9CuQgCdvzgZq5Tpj0RWx0X6+XIgI8Oj8wgWDCSIcGxN6UQiVpjwWc2GACIjb3fsO4q0S2OFvkN+1vv0FiHrP9R5\/z72D0KnEmnBc0A9ugR7HKEEQhBhq8l02pQ2ir76lxb7cnz7RhRmDkTMH6\/0kUqC505kRHVyoNrMuJQ8\/vzPnRjdKw6jX9mCzOs6+wUQqRTk25UbJNJ8OSBBwKr0oTfvSsVqsOBj93NgQN3Bv5ebhlrdlBF0jtT8\/Hzwpy61+Ph48CccWrUBiGSWlEL0ALLrQBGW\/6cQ3RNbIG1ga2FViOtQih3fn7DM5xJ3CY4cOYKTm54Mh6NycTp2zWXXPGSCnGv2c80w48lGmDKtFFOm+edPwa4GWLokytS8e89ZGWiWnAJHW6\/ppszjVrqa87r+NrZRH72vRyigYPWBEzgmTZqEjRs3Wh1SK\/p1794dM2fODAsQMT3C0O\/j5Vu7CamBR5pWJzfS4mNCeRr+AGTDzl9w7f\/9D2N7x2HtrNORv7uBmH7tl8csX+b48ePYv38\/2rZti+joiqkILE8EYOTgaOS5TiDB2QBL3zoezFCfvnbSJOdavug0vPR8U9wzqVj8GLX83VEYOTgGBbujxNcE4+69ShDfoRQbP2mAjZ821IbFN3EhaxHQb6gHQNRgN\/5tFCHrr8\/+RRO1bOvyWHLK+IUV4lvssMJIUzs3Y1xc3SgPQbCbM2eOcEJUwxtCXmA1PLACgOiTBlFp+tbXP2mgUdloXDMAeXPcH4UUEkorOlaEPXv3IrZ9e0THVA5Aup3RUEhBPCZs2V4SCjlijJ00ybnWrIrF3Xc2xshbSjF\/YUUJhHSTftnoiTpsRAkSvF6pY8c0wNJXPZIJj0N68y4liwOLHxDZ0aNimglza4tr7vcBASt9JEjox9plhZEAYrQZXXAhBx67dR7ykI7a4UUXiOaQF1kNDjQ149IP5PnVeZg9\/EwRvl9ZMy7HD1zyAQZcl49upzXVbn33gSLM+NCFmy9uj52dW8LrFyUeexKS4ES5eO1PR1JUVITCwkLExsYiJiamUmxNpFs4HbMA5JZfWptTT0N10CTv78MPYzF+fIyhBUbvWKZngnQwYz\/qQLq3zRFK1Y37k7T5CA4yxJ7jCSBNzhvgo5+w0ocSiD\/pxaiwVLAeqv42YwYyML3sP33jZzUNJBEAMclYZrZrCSA912Xihz89ZtZV+z4LWUhFqvY3naekdKD6M1gBkEAKWvU7V7YHQPi\/pLSKylyVBm7KHCKdsyLgVJYmlUnqXGedFaPxQCqbg\/FK5byrVrnRK3YDCqcn49Ytufhki1NIJRMHVh5A1COOUYYxO6wwRptxERYhDR6vYa4ZvngoiUhpRL+WZLiB6gioOgHqFyk9iRlc2bdvX3z11Vdo0qSJ8EzOy8sTCaU+\/vhjnHfeeWKY7KfOEQGQSgIImSkfcinKy\/JKHcjbE87HtM4tkSNFEGcUKgCIH+nAaLPqAcMf+JAu9TvkJoNiMFxOOJOzKzheqe7jyE6Gy+khmAs2F+WmI0s0mUg7cgGqc1ECoTVFNm58NRaGxxajIDrVtZ1jP3ghB2fkJGN7UjYG3eWJujv8WeWOMP4kD\/1mrKwVRr8Z+bwS4XHN1UsbBJBkeOzX6rojgHz44YciupqNkdKZmZlISUkRCaD0jdeUYPHCCy9g4MCBWj+Off7550W0Nr2VZb96ByDMvk6Fqb9W2YA6MwDp3bmltplRagAgHsHAc7woU7TKzRwf78Yrr+z2OcLoAUPd+HpvTB9QyE30AAghwZVreIyR\/OGi9fT1AIgzOVeTDj74oEhs5L17PUcOHhv0NHHJq\/fjj+96MJLz+H1QFr4gTWsuikLs9Gx0uzpJ0L3m7Z\/Q5ePkgGkG\/SlR3T\/l48DrGTj1\/hVVXppBDyASJHjkzVZAXLKBAMI+6gspEIDImCeOp1RCYGAeF0oeTPPAnDB\/\/OMfMX78eBw8eFBchukuBwwYII7SmzdvrqAsDXsJxMKaM+yiOon1P6u1YcpDKmSzPivE7xd\/iN9unoGvfzwssp2x+Uog1yNHOC9Q6+f7wNlX00\/QJVuRHAgg69b5AkggwNDfCBcYgYAg4Ha7scHh8S+gTKH3\/Fb7cvN7xnmbK1cDhK3ti9CliwP5+Q5Nx6Cnie9FKXCpgKj3g6kAIBmL4JpefrSTl5fShyqhiGOWNxZGf99DEhah65+SsPR9pwAQusXf2NWTMYxNZgnT6yhUMy77NPlD\/wrxLer4UNeWv3H6zSh1H\/50HUbfGx1hGP\/Uv39\/ZGRkID09XaSdZDwTGyOx9RKIPMIwZqreSyChPmSCA\/Oo3nJprCj9IKN35XzSunPJSYl4t\/RN5KSmYXpG+RGGfiDL\/l2IEZfE4uXWM8rJSJ+O1JwspCBVE8kTXeUbPdeZDfl3vDse6xq+4iOByM3Jzc2Dhbrx9W8pKUkQQLYWbUUXR3\/kO\/IrHEsEiHmlDiFxIFcIyIQBqntdyC4HkKJyCYTjKPWoNKXnAGlry88iPLerRxH19+Jitwjbp48LAcGVZAwgoT5DOU7qQSSABOMHIucw04HYacaVVhgp2ZKHlDKMXhCUQPjcKaWw6SUQOcafboS5YCIAUs5ZWxzJ9KZZo5B\/9knJ+gad8s7FvAPLgKw0IKocQPwu+tIoIC0LUzqkaE5T3Rp1Fm\/8uJIE7DjxPeTf7X+PQ8LlO8VxgY2+DwXr+nr7dsSOEx9rfeVY9br9G7Bvnvjom6NbcG7jMwWAmPVVaWDfRg09R5q4ko745ugHuLJRX+yN2Sv+XiVo6Kd9n5qRg8zHPL4awbS4W7LR6PQ8TQIR9+r1+dDPI78ToOMCJk46htduHyDuDcnUFvuamRa+cBD9d14e0hGmYes4SP8P0mFUpqEqzbhR8PBSBQn+LcHFA\/fl8qQ\/AKGyVJVAJE+NdCARCSSYlWvQV4LDzBvPFEcSI1+R1Vt\/xk3zv8axz5wo6rpWAIiRElU6km1c2g8QiknGomchbk0K4jucEE5dfVp3Q0HDPLGx1x7YgpEtr0R+g10oKXEDibnY82MjDUCQm+jTt0\/rQShouAtxJaWIv\/x74bAm5x3ZMhr5DRog\/sQJvLz3kJh3f+O9Yi7Pdcq\/5wvMd6yHBtkkfR8V\/g9XtDsXe04q0OhV6R\/yyA48P7MRSCc3+tKvtmLkL+NAL1LkOTGj9SOCvu69ToCOZGKu6DJQm8Y3bCoKvODDsXQqe+iekzQauvcs0Rzy+KEEmAu7\/46CdZ0ETeQtFvkeg97\/6zqcvfmWkPxAGp\/RA3Qma9F3DPbOG4PYia9XqPNiVzCdkT5BNeESLCiNEMwJIGz6440\/AGFfNTM\/\/1aTS8lMdcy4f\/\/994sjr5r2gXqRcePG+ZQEUecMW0cyFR\/kptdjhl6JagVA5BwiijRlEeBMQ5uiYjSLKXc55zEmtkW0Z6F7rRo8LCAnC0hO1Zy6pMTBBbKluDzpMx2t\/rAvGvkU8Qk91E5mJ8OZ5NEwxF3+Az55hXN7FLDCyOJVxm6JK0G3Rg21v\/976DcMcFyBwug9cCLBKzmUf0+yZBKfXpeVYtWacoeubgUeCYOt+4jtYrPKoxDplRKIoCHRo3A1awSI4Tf\/jtfG9kN+l\/VwTueGKNscXncHZyqw5W8lwoFMo8vdR0hUpL9Rv7XlOpDURR4JkFzIKAOP6eVOVpRSNi+vvBlXSCGzhhoCCO+VFpjDHy8UIEWpRbrAN+s7xnI8jD+FpGrKlXyVYBLxAzFbacF97\/cI489j1Gh6K0cYIwC59Nv9GHj2KeIroQP5TyFGXByLj1uv8BwlXEzqmVb2b5ZnoacASRmAi2ZT7+6nxUM2KlH5nRDNuTlysj0KBzqjJXnzhWYnC3CSjmnSVJuUlq3M6zHdurITtb68Tk6WMlaAj1fRm5wtwI1N\/puT7ulLvQz\/lQBSYZ7EbM+GF4jnvROCZ1KOd65UuIgVOpAR10kt14E4U11wZjmRQ1zwqlScpV7rkMuJpLRcwUYBqtMzgHQv8ugkEOYWmXtH1QCIlZyojlM6aqBitpSteKLSCzUBCRWcEc3mrqrv640VJliPUzMlKqWZD77Zj+ZbumoSyKqtP2kRudIK02rNpfjhn2cpEgg3Txaw2AMgfPM7c5Vju6tcJRnvXgU4E70Awj2ZBFeUxydDbDjhr+ExRegBhJIAjxGa5JPjUbJJ6UVICgqg8DspZYB9qUvwXoeKUqlkFQACeGiirkG9jvxb4gZpzEoTytFyCqmUdcHJvgQ1JXWhCiACFHPTPDQlevUaSTninsV1hUXLBUSVMU+TQKhVLktPJu\/Vm2jIipu6WR8zCcSODVoXN2NdpDnQswqoRJWbnvlAzJq\/HCJqnV2CzP89G42i5Ws9koXVJgCkzAMofbFHqqBYIY4l8rXt3VUubnoGsCiva6kg5L\/8nD\/qZ9xkAg2836s0eTdeUkY2XFnl4KON93Md0T\/b8\/YngMS547Axxhsxys1KsgUd3g3M66QugjOr7HynABMB0ANTXomJICad27xWGKG\/kPdEsOD9ZNCrzHss4XdeqUtMTsU1ASTdcy0BaN4mkwxVJpiOxxG2mgIQNQbGaHnpwyKsLkG7+oU9gKhA4I9plXEk8+hAXEh15qCPcgEqv6TjVpl8XbbSucC5w1IBV7on7bgQu70gUS7vq7K\/d0YjpYLYtWXjvYDhfeuKDWfWCExSoWva1ytlyM0MIH5nb48EogGO97yj6h9obZLApgKTvJ4mRZQBhLReSeAgIBBIpMShAogcz76pTN2eXg6kkgfeedQsZXofD1aLM\/MDUSvK1RSAGOk\/1Eem92rmd2q5Utb\/kWVI+TmbLP4lLTOTJ0\/GkiVLhGOZLE1KZezixYtx5plnCt8Rtalzhj2AmO2PynxPieaa0YdQtDwRf8KveNdZHnWrenP6KgSmewFEOberRAjJwfPertj4mva+7cWXqqKhMndiNFaZOxAQqICgvP3pcWu56c3fHMu5vI53SForgMKZmwwsTjF0NtOuRfDg0a1sTn8AUpU5PSzfs0HHQLEwNOOqAZh8OdEHSA8gnGPt2rUaSEggmDJlisiav2bNGvEdQ++tAAjHsWg6+86fPx\/33HOPTxH1CICE+MR9\/EDuPgmn7\/0Rn3Rro9WX8QUQ9SImAGJKT1WChunF\/XdQ9CaikwIgXPiafsVoBiMA0c9HuCz1xIWIuXKS4FzrCybCL4SKYUouCoBYOcJU4s7FUNWRjNnY242fjz2zbhK+J0bBd0bXCwQg9PXQAwjXmFEwHSUHufHV6zDW5aKLLsLy5cs1qYK+IYEkkAiAKByU3qOrth7QPvXnpm62oCh9TPnnTlx+6EI8c7sDZ\/+8B7ObNNKUqHYCiOkGVIjV9zUb6\/u9Hpz8g1WFcTrzqR5AtI1vBUDMmG\/wPS1WjvX9EXcy7V0u5Lcrl0DMFKQhXM5niOpI1urq+7XYmV\/fm42jmz6ynJXdDgAhYarXqXrckMFydBSjlMIjiwSQp59+GvQBUds111yjlSqt9xKIBI\/4ltFaEmX5GZn2csrZIj+I1UYAmb8hH513nYsnxjQQAHLXkWKtWPfpBy4Vb8oTrfdoUzryHXDHTxVHGEdmZrnp0dvDHe+ucHmOYYuP7y3mc\/jJlynHsj\/zU3r6esZWB4i4F5ZlElNSVnTe3cVH6iANO733x3tQ6XV3KM9C5ui8QViKpAct6XfvHAnkGCQxAd3g85GU5ELv3mtEmc4LWxajwYJbcPobvgBilg\/E6nM36qfqUtTgu9KiwzDKE+LvWnYBiDq\/6oFKZy8ZbUt9yLnnnisidyMSSDnHAvqBqBYUOSRY864cZwYg9AO5c\/lWfHTZtT4gwl3mdk8RANJhzUKftVS4rFDbWOoXHS7vgMJly+D2Jq7tcPnlFdagOtYDVB4w4tjqaEeODMH+\/TO1SyUmdjKgsfweHPldNBpVAGnbdhKaNl3pM3b3FWPgXjAl6NuQOhAryYKCnlw3QNaAaT00Hb+89xxOSZmFwswrcVKXSy1nZQ9GB8LLGx1hZJAcc3mwqSH5KoCoFQUpmfhTokaOMMqDlptelTakv8e4PsFtNP0RRphl\/dRy9Vlr7JMC3JaZjzFj6lb27UCbbM2aeEydWvnM3K++mo\/OnX35snNnPL7\/Pvi5+\/TxONxV9RFG8kW19PCzYJW1ZgmFjPhvZoXhGOlmrs\/3wetREokAiEUJxM68IFKJelNcN0wbHI0f1zcAyowFpq0PnbmArLLcnbWlyrwpzRY7yHwe\/pL\/+JuGzmQM1bda1sEiOVq36lCiBkuTUX8jAIn4gdjBWetz2BKNa\/VyMrbm8ITz0SahORo38mRet9JqsBKjFfJC6sN6tQQCAojV+jdqnVurJS1DIc7IDySUeapyTF00idZFmgM9w4AAolphaH2Zcf0ZGLdsC2TUbVUujqqaWy10zWxStPXTGUiNvGSkpTwTkw6eix9++GHcfPPNWir+QP2DoV1KITExnyM2doTIcOWPpu7dh2HBguKylHuNhDK0Q4fLRaYs6ehkF03B0B9K30AxMXbFwoRCV3WMqTcAolphUnvFa5nZX\/2sENnbDwRthamOh2N2DZrrnnjiCcyaNUuABs+zTKZLrfrUqVPx6KOPiinUPlJ5xjR28mzMzyZOnGjY34wG\/fdLluTjlls8+goeSYYNexHHjm2tQNOjj\/4NhYVP4JNPGqFVq1\/x5pst8Ic\/HBAp9oYPHy4yaNlFU7D3YFd\/Klab9hxaoYaMv\/nlZrz33nvBYk11oRUUFIhiWGEfzq9aW\/YdLtYA5MjxEhhZZ+rCw9PTyAVIJ6EbbrgBs2fPFm9+ehFKaYPSAFPU3XjjjXjwwQc1j0SpTNP3D7VQkFqvljQ2b34AN9xwSEhF559\/PfLzG2E9dUZekOFxR0b9EgSdTqdItcff7aKpJp5nsGUdIpXpauIp+V4z4BGGFpf8X47j\/n4JWLi+AFOu6oRhC75GZQps1\/wtl1Ogbj4CyVNPPSW+JID07NlTO8ZIKUS6NEvg8dc\/1HsUcUJ+LFMElQEDfsTrr5+jTa9KQnR2CnQPodJUneNU648MyjO7fqQ2rhmHqvZ7UyWqPqnQE9d11py\/qpa0qp2duhCXyyWkCjNAqC4AkTR17\/4QeLT55pvNGDPmCnTsWIr335\/sA2p6vYzZPVQtN4ObPZAORCZwDm7GSO+a4oApgNQUYVV5XSl5SEWp2ZHECEDsPi4EQ5MqecjaJWb3YMbPYI8PZvNFvq8fHKh3AMKNyszaqr7CTCmqBxCz\/sEunWBo4txGyX4rS5MEkFPSnsVPWQ+g\/cTXRKrBSItwIBAHTBMKMRGyvlUmH0hNPg7VzCnpkAFQtLKMGDFCfKzXkOsBhH3UuSqjUQ+WJjV3hbwHaXauLE20ghz9zFP\/xKgFY2I1Gm93SsOaXEuRa3s4YEtO1Agz6z4HpATSdswc7F94b0QCqfuPtFruQAMQWlzYZPrCbwoP48o5\/8O9fTtg0gBPXolIC18OVPcRRl8XV3K2slJO+D6h2nlnAkD04CFJnfOxCws3FGBCn45hYXmpnY+g9lBV1UcYeadqrM3Bd5\/TnMd4\/egO3SyXdag9nKu\/lETtOVRU6i9sn8F02\/cdRRRKwaonDbxVv+qqDqT+Pmbrd14d1hj1GofWLMTx3VtECH91XNs6JyI9rXDAL4DIwaHm\/5B+CszapGZ5skKUXX3CLe6gMnyRz4Nz0AGOHrdqU48weU8PRdbmw7jp1F\/RolF5wSy7jhfMSJabMQjPbI+Bo1Us7mywDomTX8fRrz70KTZVmfuNjK0eDmhHGBbGNsrxYVSm0gppMk0crRwy2MvKODv7RACknJtmAKLynZYe8k4G9dn5TORcR378HvePuAaNz+qFR3rH4thbj4mvIo5kVcHtqpsz6p577im95d6puGneFyjNmYUhV\/bBM49Nwxc\/7MfQMfciKu9zcXU1QtVfRKskUw1A42cEEcabpKWliejRl156SYsiVc2S0qTKMXQnZ2vWrJkwq1KKYbAbA90YAKdGoarsUWmT11KD4Bh8RpMtm2p+NaJDZte2MkaNotXfv8o7f99JwGVgHDOBk0Z\/90jaJSDs2rXL70ZXr0XfF7bmzZtrEojKK8n7t99+G4888ojoS57PmzcPCxYswKFDh8RnpJNlD\/773\/9q\/QI9N398MaKn6pZ5ZOaq4oAAEIq0XJDXDL0FuQ4nis6+ETHfvoFEtwvvvP4qtm\/fLnwkuOHatGkjgODuu+8WEaDcXHTK0ksZeglESgPqplDfdLxBORdT4RNAvvzyS7FY2XjNtm3bis2yatUqsXj1\/hf6a0pQYD8uZAlKvF9ulLlz54r5GUdidH+yWDKvbzRGxp4wHoP0SWlLvS\/SKq9D3qn8UiUkydcLL7ywwrWkt6m6CKwACOngEVLloZxf+r2QN2eccYZG133XJyHjtsHYdCAKmZ33IrphKZ4t7IhtiEXW4ldAWlS61bH+nhvXCr1+A9GjP1JV1YKPzGsvB3wARC5uZp9WF7p8k\/HteMEFF4jNQilAnzdDJc0fgOilARV8+EbkZpNvPblxuVn8bTzVo1SOl4WBjDaoXMySvmeeeUaQbeREJvsYjZHXVSUXCY7yM\/Xty2voj1TqseK2224DK7pLEPJ3\/NJLMJLnej2T7Cf5qz\/CMMpYPaaQ5l1fb8SdRW\/h9eZXY9MvDbXI3gdGXI3fd32Np559Hm373CyAQD+W6+HPf\/4z\/vKXv0BKRepLQV1TD44dhe8f7omZ209G47OTMOOZWRV0MvYu88hsVcUBQwCRxw0uCrXpN4j8zghIzABEvzk5V1UCiAoS6j2RdilpqGn6CXRyjNsdj4MH7xV\/u91x6NSpE7p124eNGyfgzjv\/JNL9qwAn723XrgZiXOvWrdChw2Vlvx\/E4cPv4pVX+iIpyakdQzivBBAJVmb6GzMJJBCAyI2uL0twVutGePGJKVj2v0INIGR6A\/dPu\/FA\/F60n\/wPPPqXzAolDQiWzBXKcgdGwK\/n0W8\/78N9Q5Jw4vB+THT+jJgGpWhy6TDLCZWrakNE5g2OAz46kFb\/noObrkpG72tGYFRKGo63PhMdkm\/F23eepxWAUqeXi9zIymIGIPoFbheAqEebQBKIPzZJvQA3xOjRf8Ftt+3G5Zd3xMyZ52g5OLh5J0x4Glu3XoI9e67E3LmH8fLLKTjzzIHYu3ecKIAtfxiG\/\/vv29G7txsTJpwMJr8ZMGAZmFT5xRePiyjbqgSQjh07+qQpkJtbL4GoJtRn5i+pACClvx\/Fvc23wPnqgA1mAAAgAElEQVTQ36F+L0s7cl69lKNKjhJA\/NFDoAolnD+45R7pbTcHogYPHlwad93DaH58L1bMuA+jbxuLzW2vQusd\/8CeHV\/hsjF\/waZvtmLD3x4WOgf5VpZHEX8aezMA4TxmOhCjN5kMv5c6C\/UIE6oOhApBqadgjRjqShYtSse55zaD271a3HO\/fv3w3XffYePGf2Pw4OuwdetWUPfRps0QIVls2dJOZAo7eLAF4uJ2Ij5+hxhD5eNbb\/0Tu3ePRq9exSKjWKdOu3DrrR3x9tuHUFh4B958814f3RL1BWYSiJWFEIwOZHTKLTjnxG5kvPqBD0BICUQFkP9s363pjIx0IP6em0rPiaIjSL15KLpF\/6pJIA2atEDslA8tV6azwoNIn6rlQFSPS3uW7t1TiP5X\/gn7CvNx+jkX4Ju2V+Glm0\/Hs4+na6LqI+mPYWzKSEGNqr3n30bBZFYARIIIrTJs8ogU6E0WCEA4hxpQFqoVpk2bC9Cx48tYubK5qHFqZIWR9\/fVV\/ejSZPGOHGiI6KjY5CcnIg5c371GUO9xvjxMzB+fAyKi4sFEDVqlIljx3qgU6c0zJrVHHFxOzTltF0AovrimFlhrh54hSZhSIAgP6lPevPNN6ECCKN0jSw47K8qqvW6K5WeS9p4\/EvaXdQfM559PqIDqdp9XmWzV3Ak0\/t9hOpIVmUUV8PEzJbOlIG6IusVrsxjCrOIsR+TI7Pl5panGzQiNTnZ05\/j2JiekJ9ZzcpeVbcfiZStKs6G97w+jmTXX9AOdF9Pu\/Q0zanMX5xMuLKF+UnXeuvVmAGIBBoCCRvH0uocqH4NgYP9CBhRUUBpqae0w6lH1+HMbI+vRrAFlur6s2BgnbjvoX+u67dS7+gXAEIpQxaRSrs0VkTkys+qIv8p37jBFImSG646ng43OH2uCCKBAEQCDf9V74Vg4nVdMSSX8yckeKQUXoPXuuwPBbj5si\/w6j\/aizH75t2G06b9q94k9IkASHWs7Kq5Ro1kJJNiv6GID4AvdJaGzvZ2CNTfbrZQqiAAmF1TAo08fvBvAokVAJkywYWxk50CRNh\/7h0rkZZaihWfXoqomGYozLwKLa6533J5A7t5YNd8VhMeH1n7irhk0z632nXpWj0PFfX8CYdWawCExbWX\/6cQj14cixOtYwSAyGp0+s3sdrtx+PBh4ebucDgCPodg+44a5caSJQ5kZjoCSiDsR\/PtBRe0EnoPKVX5k5YkHbNmNcMER39M\/k8Wulzi1ABkwLBO2LD1dHEvBJAm5w2o0yJ9XS25UB2bmjVsZs6cGRYgUuEIw1D9l2\/thtGvbMGOfUfR6uRGeGv8eTj3tKYBeaseg1jFTi3IrR94ddqvOLmny6fPhp2\/4Nr\/+x9aTb1UAAibPwApKioqM38WIjY2FjExnr7+WjKSsdO9E50dnZGtyTTGvTkvpQiCwIwZMaYA8thju9G3b6IlAJE0v\/RSLG7\/vjFu3ZKL\/jc5xbXeeXIlwg1ApDWMGyUuLpJbVa64jRs3Ys6cOeFTWOrI8eLS0Yu\/RfIZrYXilErTt77+SQMNq9G4HMeI3lsujQXnG9s7Hv3OalNhp27+8Qj63HwA\/Ucd9AMgj+JE6z1wwqlteL0EEgyAJCIRLLjM+XI1SPIPIJQsnnoqBkuXBpZAKHG8\/HIuOnVKFIpQHn3Y\/Flv9AAy5qdSDaz+es8HmP\/ENjzz4U214gij1sVVORVMOL8dfizVIQ1U9zXCjS8VzLjc4M+vzsPs4Wfi5JMcQplqVolOSh+Z13UWoEHQ2bHvNy09onxILJd534rvsGt1HFpelmcIIL9OHSEAJK4kATtOfC+GZj7WAFOmleelOF50HIV79qD9qacipnFgCaRbo34CQOJKSrX5JD2NB\/QTSgh3h44oXvMxOO9bb\/2CnTtiRRf1mvqFNnZMA6Sm7Ub\/vh2xZXsJlr7aAEtfjcL\/\/a3EUEF8VaOG2Ol2Iz\/TgbUbE5HdvbxceIdTDqJpzkRcl3mHuExNKlGlObd1yrOV0sGE20axC2jCjS+2AUhK1jda0W1\/Ugs\/Z3vxvfZwD\/sZ153XDic19JRspA5k2X8KUTTQo1Bju+c3T63ajU80wNL7jmufHz9+HPv37xfRudHR0QGfbemQQXDs3iVAImrlBz59Ey7sBsfuPLg7JCDvyy2Q804c70Svy6Nwz6Riv3N\/sq4U2auKsHVzS3Tv5QG352c2ws59Rw3H9GndGAUNo4AM4Ol303DxzCxhhaFFZvFizzFm36wbxdiazIlhV1YwuzaKdNhT47Kk4yKd2VjS06ykqOrAJh9OoFQJso++3IY+g74+Bozfr127NmCxc7v4YhegVXYeSwBy\/R2PYs\/n\/6hwLfkQKIGMWPA1fvrNjd1lQHBGu8aIbxWDxWnnCCmG7d95vyL15W9QePg4fmnRF+dGfYLPHu6ufS91IF989DA6HtoLOBOEVMDW5ImGiJtaKhSrq06cQNGxIuzZuxex7dsjOsYYQBpd0Rdw5QmAEM3pxLHtO33uofEZnT1mEO93ct6S4liMGxeDVWvKpR79zbPvoAEN0W9AND5ZH4XiklJ0OT0KO74vNRw3IKpUWJcogXT\/aw4WZQFnXZ2E1FQgJSU4s3ZlH7rZeJpVZZpBs77+vrdroxBAPvzwQ9x1113iUvx78eLFmDJlikjJYBVAGPvDeBzG7hBQMjMzkZKSItIT6JsKOGr0uFqPRz+H9MyVe4IezLI\/wwHk9dTUEWbAFyrvq3OcABDpA1LhwiW\/o\/F\/s9AsxoEPl72AhLYttS6qu\/kDUzOQNHszBv+xLR4f3BlXPPsFElrH4LXbz9P637JwM9779ifx98EWyWj1azbOPLUJPnvIU1VdAsimhSPQ8dAejzKB5g3WnqDXZnqiV4+RDSMdSLLe\/EvXUJcLLqKOMAuXzyf7Zicmwuntk5zrRLw7Hq\/sfkUoZxs3jvGYZZMWafeQilTt9+JiN7ZtK8LPP3uUy3p\/EL2fC\/sfOXIYK1c2Q\/eGS7BxfxIc7ZyGjmfHt29EYeZAnDj6q5hbRqkyFSAtNMe3fwJVH+Hvc9W7NPqMXoid8j4aNGkecH0F8kitCR2IEYDws6uuukq407M8qQyFUPO+yEhjGW2tB5AVK1aA6SmM8pAwV0qTJk2EC78sQhYIdAgWmzZtEmOkBKKnW0pLZL5RHFd1bno7rxXQjKuiqBp1KQmQ30+Y+AiufHE7Bp\/XTgBIv+e+RIdW0fi7AiAq0QQESiBjesfh3r6eHV5BAinTQRR\/79GBUAJBepSmF6Gu4spGDbDXa4GhfoPHA77hPbqOE9g7wCtdAEjM5eflOpW9Xt0HwYONIMM+7X+PQ8KJnWJezjP\/hxO46on1cDnXIm5qlo8OJe+HE3hx3hG88WZLdGmYh3RXGjKcWXh\/ewch1Ex4oqGYl5ITpaarCm7DztVxcDga4dqZwNtH04U3KnFS36jE\/PWd5ypseFUyUDOY+\/tc9mnxp\/uq3bfETAKpAPh+VrXREUa+5dX6xlIyYRwRN7GUWDit0RFG5mphICTnoV8Gj8UTJkzQQEV\/hNHPE+gIYwQgpIUSU70BEN6wTI6jZ5YqsqVOuA83L\/gaB35zC11G94RmaN6kEZ4a0hkPr9yJ54d39UkHQABJ\/H0dVtz+R1yS0MIQQIo7dMSPb\/1LfNf5\/5ogd1GU+D2+JAFFp8Why7ocLB84QEgQ1G+wqbqOkS2vxFMD1sPpdTOnniNq5RbRT+o+NOmE4fdOiHkpq8T8WACX04lRqz\/Gkn6d4MwBFnT4s7DlsHEu9k0esA\/9BzXFST0HI3O3xwxzaPgoxHcoxcon+6Fg9y4gIx1xa1JQsDsKHmKciH\/ahe3dXCKZjlEzOkJIKUP6h0iQaXffa9g3e5jmNyI\/bzt+IfZkDoBUhgZzLLGjZosZgDB0SDoMlquTK3JDvxFVSYBR1PIIw8+llDF\/\/nxNKjGSQHgVrl+2pk2bon379iJRFvUs9C1i2ke59vVlUFUKjXQkEQnEYEUbvQXUHCDUgVhRoq7e+jOGzv9K6EDuO2sbnrupq3Y1SiB8K+3of773WOFE750\/iO8LHi8rLDHdAyBs3Nz8SeK5wfs3\/xUShVefQR2IY\/268rtxeuZzIQlL+q3TgIWSR26id28r9+6RSpzITXR5971XVPBe4+BX32CA4woURu\/Bkn553vmSsDg1WxxnxP+ctP8kCVpTcxYhBYuRhBzkJCUhI9vYK0U9jpAcGeLeKP4sHylC5s5oP+Uj7J83RvNc1T6f9Cb2zk1Fu\/ELRHg8QeHopo9MjzF21WwxAxApgfAegwUQCRRMF8lGqUPVjchjiZSQJ0+ejCVLlmg6EAkgBJ\/zzz8fH330Ed59911cdNFFGDVqlKEEYjS\/XomrKlHrjQ6EsTCGr8EgPrRixmUfSinuklKs\/+0SND5\/BE75sj8aNlCAQYCAV+kpgCEBSMgD8tIFgKgSgxF5nu+dSM7NhdRvqP24kT3vPaV5\/cmlpOLb3wMsHCePO+J7pxNFzAfSvwsc+fk+dDmdSUjOBrISczSQyknyjFev4UpKgtNCCK4pUFQBgOyZNUyUtqxMzRYzALG6vMysMO+\/\/z7Wr1+vlQ7hvDLlJn\/3l3FOHoM2bNgg8uVS8li9ejVatGiBiy++WJBnZoXRW3LqpRVGBRBKCKEW0zZzJHv9v3vwwN+\/wwf3XoDl7U7GPzfdidlNHteczegxmp6cgySPUOGVLMpe5ElAnnM60qdP1z4z2uzqGKEylSGyVldqsP1kIIvBOAJZIBolCEklsVoRTm\/CFQrVWUPR\/q7F+HnFtApHFbuPMOpRqckfB4prx04MvmaLXQAS7GMJtj8Thr\/yyito1aqVMOOPGzdOWGqqqtUVvli9f02JqpcirE4g+\/mL6JVOaPcs34pV2w6K7scGUhLIwI35E\/HqmHM9wJCcCGeOTjooe9kvSvUFkGDpqrX9FSuTSiM38N7nhqP10PQKR49f35utmVerUomqrxx38O\/pgsRg\/FPCbaPYtY7CjS8+AKL3OPWXAVxlptRmFzc6WUsJYBQL8+7mfXjg9e+w\/0gxigYmonnM42U+VdPLi1l5za76B1UTAGJJgghxRbmbeAbShCslEP1UqhlXTfOn6kdUs6y\/z1WTbHUnLA63jRLi464wLNz44mPG5RHmg2\/2+7igS9NVz549hbLKXzM7wtz\/921466v9+PzhS\/CkIwpZn9yB0Qce0hSpLgaXCD2BrxQiASRl0fQKx4IKugkD4tzx8ch35GtjfcZQGUvbik7wEboSeQxRTL2mxxL99Q2OOQSQXVcBsdOz\/Vph7FqsNTlPuG0Uu3gZbnwJ7Ejm5drpzdzo8sMSZPx5qqHnnlUlqrTUzGt2Ev7+xThMLp6GyQM8uQAlgMSXeHS60oOUALLzspcxfXoG8htKa4eHMJpTN5we5QUdj5VEAhC\/E31K3OiSF4X1XVQQ8YynSTgpNwFLOnusMh4lKz8vhTM3V\/iTSG9Vzqe\/vlxUcpwKfvL6Drqvu1xeEzHg+KWgVgKIWht376yhcP+0q8KeqQlHMrs2bm2ZJ+wAxC4rjFUzLpW07tNbAte8ilGtnsKpzU7yPFvm9cvLA\/0\/Dg8bhcafrhcg4nLmIHURkDhrHmRci7SkEAB6f98R+Q3kYs\/Fhs6dvL4h5fEt9Al5\/aq96LzeY+Hhhvc4l5Vi7YEtOG3wlWjEeJkSNxJzc9GwoQNrDxwTfeV3nnEeYJOSiB5wlvTrqwEYQYjxNZLm9r8XY+O7H6PtwW2ImXedJQlE71BmlydqKB6qwW7AcNsowd6\/v\/7hxhdbEgpZ9QMhUxmRy3D\/6MHL8ErbGeV8ljoQXczK0gavYOQJT6YqRrTytOF2XybMsUwmtKXYN75FeqDSdZ2xL2rczC\/Xnu2JfQEBJFe4hW0pLtFoMIuxkZG92ckeK0u8240uuz0Bgpynm5c+6VzKz8ppduOj3F1o+fM2HJv1J1MAkc5cqq7DLk9UIw\/Vxmf0QMGUnoaSh2RQTUggZmZcK7Ewkn59NYFAlRVV\/Z\/q8+QvoE7vpSpjaFT6afalY1rYeqLKzb1q6wFtU\/U65TcceO8pXHvtNSLKkAxJTU1DYeGPcLdMhOOKB7Hoju6Y9MYOmIXzq0Dz5mmzkA6Pdl80BUBU5eIiLIIag8KugfKBJCLZm\/+DDkq5Pn0HxQwS33max31JdWIyyzOiOj8xwxjbDrdbS2oUyDlKzt3ywHc4+NSggABCyaN497difun8xd\/VTGWheqK2TZ0lgCKQh6pq4ZGPh9c78unrlivH2fWmtSOYjvegur3zb7nhb775ZsNoXtUHRDqRsfQoK++lp6dXCMqjVywb9YSqAxvLtMqgPc7JCPLHH388fBIKySOMBI\/4ltGaEpWf9U99CNGtYvHuc\/ehQUmxqPshma4+FDMlqj6vSAYyfACEZlxKB458B\/LXr9dwZTEWIwUp5UAjJBC3qGDftWtX4cqstlvi4zWIWJ+fH7Cvz0CTeSvTl2Mlzb1iHei0aqqpBMIxqveoBBCZK9XUwcyPJypd3PfOutGvh6q\/cP5gw\/yrEkCCDaZjEXg1klb\/LM3+JiDQPZ7lR5977jm\/UbxyHulWP2TIEKxcuRIsOk7vWPIkOztbzGVUS8mMjtr4vakZ95Gp6dh\/5k2YN\/oSHNqzC0888QRmzZolEFh111XNuGpmd2kafvpfucj6rFDjweEJ92Ogoy8udnpiYaaX\/SdbyznlUb\/ueDfaTmpbG3kXEk23xh7E7f3PQ7s7s+BoaxBJp8xaEwDCy0vHNlleQv4djCnYDEDoOCgzxQVKNWlHMJ3eDV09yvAYQxd2fwF1XOMTJ07Eo48+KgwIZgF1qmTDqn0ECxVAli9fLoq1hR2AcOHQjDt\/Q76WKYzMu27MRFw\/5gE8cPW5gsmffvqpqLVKRDWL1g1mh6ni\/\/LPPw9maJ3q2z7aLSI\/VfDw54lqBCCVDaazcoQhQ9W0hqGUnDQDEKupJu0IpvMngZDGvLw8vwF1vLb6wjRaaOpRRw82UnKpNxKIUV6QmG\/fgPvkdiiJPR9tvpiLSRNStRKXZB4bdSORZj8H9AFwValEbdpjiK03YAYgUgLhRQPlqjUCkGCD6aizYOAdc4eoa5Xr1+l0YvDgwRUC6pgThNKCfFmSzkABdW3atKkANvocImGtA\/G3elRttJptifVimVEpAh627jufyfQAYpcnanV4qJoBiFWumVlhrATTSQdIf1YYfUAdM\/2\/\/vrrWl1o0sr6xgQTAgutKLLJPSHTXsjPpYc2M5DJ4L56ZYWhO\/qM68\/AuGVbtFynVh96dfRTF4N8WNTLqGY2vZnOSPMeqL+V+wiGDv35WY3mVBegaja0QoPdffRA1W78fOyZdZOmeLVyPbsAxMq1KtunOgPq6hJfrPBVU6KqVpjUXvFaZvZXPytE9vYDAeu8WLmQnX30Z1NuPiaDYd6HqVOnCoUXm17hS6mJbxA1z6VUkOn7W6E3WDpUMVpKdkyrd+211\/pYt6xcu6r6SPBo2CYera6+Hwdez8Cp968AA\/ms5BORdIXbRrGL3+HGF0MrzL7DxRqAHDleYlrWwS7mhjoPHwrPqzfccANmz56NhQsXCiWvNDlTQuG5+cYbb8SDDz4ojl08fnEcwUffP9Rkt2Z06OeVZ3Aq+VRNf6h8sGOcaq51\/5SvAUhp0WHIPCENW5sXigq3jWIHbzlHuPHFxxOVvhz5vxzH\/f0SsHB9AaZc1QnDFnyNqiiwbdcD4TxyI9InRCq++DkBRA0ClG99FUAC9Q+WRqt0cF5VYy8zdevP1cFe367+tAqV\/JwvUgr88t5zOCVlFgozr8RJXS6tdkcyu+6ptswT1gBCJqtJhaKOH0bXbS+hMHerX\/6r+oeaeEiqM5uUAKjsqm4ACYaOQF6QVqOfq5rX+up00ifE6nXDbaNYvW+zfuHGF9NYmNqyoI0ejHzjSy272ZHESAKx4wgTDB16XwF\/98XP67KFK9w2ihkwWP0+3PhiCiB6cduoEI9V5tnZT5+vUk8n\/9Y7AekBRN3MRv2t0BsMHZzPyKWai0pm89bTaIWG2tjHro1iZsa1EkxnVNbBzNKl91blS0o\/jyp9G\/Xnc5HWNXm9sCssZVdO1OpcxPqISF7byE6vdxc22pzqXMG6FwdLh95XgHRLU7P6nZWyi9XJ71CuZSeA2F2Zjvejgrb+\/lQPUoKGBH01C7w6ByveSZd1tT\/Nw2qZB1bUGzBgAFJTU8PPlT3YnKjFxcUinJ5RpjSh8k3AvyOtbnIgUEU6eUc1Fc6vB5Bgg+loJlcr05kBiPoE1XANFRD8PWW1PwFHX7eG9Wf4ggj2ZVVbV1XAYDp\/RNNzb926dSIqkWUFzz33XPzwww+47777TItd11ZG1He69AASTOCcEe9MJRBvCktRmi9AeQs7gukCVaZjESlmZW\/YsCEKCgpwxRVXoHt3T7lVeSxRHRJVqVNvQND3V+vGSLd21qGZNGlS+AEIGfbi2t3Yse83n5yo+sXB1PdkOEP6+XAZX3D99deDtUaZLIV1ZSOtbnNAX9yKdxOM9CHf8AET5\/jJ\/6LnnB3BdHoJRNVbsbTl0qVLQSfDBg0aCGk6kbQpTV9ASn5lFBsjgYeSBwP16pUE4q\/Idpd2TfD2neeJ8pS\/\/\/678LVgrgMWH6b0cd555+Hll18WHpVMmBJp4ccB+ob8vuMzxGV+ClscyaQEQlZ5i6gbcc2OYDp9ZTp9kNuuXbtE5brCwkLcfvvtIkEU9RVTpkwRDokScKjPUEtdSl2JLPSt78\/7qTc6kGCWPEP6mRPknHPOER6UTJpCncjo0aMjepBgGFlL+xrpQ9TUilbINj3CWJnEGwGrVprjMKk\/oGRgJZjOSAciwxDobHj06FFR4vLQoUNCkubR\/K233sIjjzwiqJRKbf1RSLXkqFYYoxinemGFsfhMfbqxMiYfQJMmTRAVVV6mMpS55Bi9tSJQ7srKXMffG48L9u677w5YxkKK6RRTaeIjzdw0dIuvyspmdt+vnE8PGsEUkTKiyS4Aqar7lfMyU1xWVpbwDObxnGUtBw0aVGWXrSt8scoAS34gVierbD+J8BQp5UaUSjSaaavDsUpezwxArParLE+qa3xttcJU1\/1X13XqHYCo+UCodX7mmWdEQBqjWs2ccYJ9KGbMlQDDeemuzr9lXhLa4vn72WefjcOHD4tcDqSXn991112CFCnJ6K9D8XPu3LniTcSmSiB6KwDnYPCbjOyVIq6USAh89AlQpRF1fmYjo9hM+lSa9LzSX1eK7YFoZ2Jf3gdpe+yxx8SU06ZNw6ZNmzR+VJeEZPYsg10b4dI\/3PgSUALRu7FTTOfC50ajR6o+xWFlH7LZMcAKgJAGbhLa7GkFkA5m9AOQtOu9Af0BiD7MXr1\/PdCotPPaBFnySQKGBD3SIcFF0qj3CdDfZ6AxKu0EEJ7b9RXpOT8zZhEYq0uSC7eNUtm1LceHG18CAog+56nebGVnTlQy2A4AkVnS9EcMdaNZBRAZY6NKA1LqCgQg\/I5SAHN90HFISjRScpE0qjlB1LKhqrJOAiCtAVLKUc2jegCRkhQBXuWnpKm6ssiF20aJAIgxByIA0qOHkKSMjjByw+\/fv19IE\/KIYHTU0YOfTHJEEzclInUMfQ3UZuS6rneT17vqq5YISbtKX20GENb6WYu1GgtYtiMJSXbt0Vo9T7gBa60CEDPmWjnC2CmB0AnI39s+kAQiUyvyGEO\/mI4dO1bQ2VhVCKs6KB5NAtFUFwCE4JGHPJ+aQKwRxBIPWfDooGQL5IVqhBLkFfU9ffv2FTo6WgYpbUpL2ccffyx8lthkv+pGG7M1Xt30VPZ6pgCiKguNLmZnPhArVhi9roEbnG9wqUS1AiDq5pd6ji+\/\/LKCEtXo+CGlETMA0W98eUQxot9f0J+8F3XRSV0GrUR62ms7gOQgR0gePhUJvYuKINIHfXwkkUBeqEZR4eSTBIsXXngBAwcO1IrBy2JPXCfUPcl+ld1AwY6vVwASLHPs6q\/3A1FFfPWtRLGeJl9utGAARIr3jOOhTiMpKQk5OTmGVhiVlssuuwzr168XjkwETmlNIR3t2rXDF1984eMHolc6kz9mhYn8vX39OSeptNd2AElDWgUpQ10z+u8DAYiaxU19gVDyoOWJVig+I5ajpLs6GxXMjIalx+nmzZtrJB4lAiB2oURknrDmgNFGMQMQ1opRq9QZHWGkGV3Nq0IdFpuMP6G0JyUQeYRRvVEjEoh9S8+yI5lqspW1Mez2A7HvtiIz1TQHggUQ6kAIMHoAUUP5zXQjjFPRH2EiAFK1K8ESgKjmWpJTmVIIVXs7kdlrCweMAKSyOhB5b\/7cB4x0IBEAqdoVYRlAJGjwgchM5vn5+aa1Q\/XkB0ooXLW3GpndTg7w2fPHX2NuDea9uPfee7X8GuxLZSlNtlSYykbFKsFFr1zlHP\/+979Fugh927JlCx5\/\/HHtY+rATj31VPEZf9+zZ4+o4zx06FCUlJSgV69e+Nvf\/oa9e\/fiuuuuw7PPPiv6JTMquBqbEV\/obMifutgsAQhvTL5RKhNVqFom9JYHLkZG9TJNgBVmulwQuWhkY1AUXdibNWtWqYjgqp6HovpiLBZmS95DktOpbRxuIqkn4CaralpCXbB8VgSHjRs3hjpFZJzCASYwmjlzpqV1X9sYZxlAKku4NKPpizvJeYPRTnPj5eQAixcD6elAUhJEakVq15nQiPkcQm12z3M89jhmxMwQb1g2AgdczrIbSBIA6EpahFSUIsUFpC0GnHAhJz1R6AJ6FPWoVfekf1Zc9HFx5kWmQn0W9WEcQXjOnDk1YhGyg7+WAcQoLZzezdoKQf6yjksAoZadImagNn58DHbudGPkSDdmzIjB\/PluXHTRESG2UoyVbt9W6NH3IYDYNc8Hxz7A5MpHMVcAAA2WSURBVLaT4XDcJoAjCSnotXMnOm+4DdOmubFqlVtgyViHA\/0cDux0u4X0lINkjHTH4YUjL9hGi3pPvEZl8tcGA\/ahPIP6NKau89ISgEjw4PFF9aCU7tqMjLW6ac0AhMlcbr31Vr9rKD\/fgcsv71B2tj6Ie+\/9BZMmeTKgPf54PujkRc\/P6OjokNcgc0IEO8\/nMZ+jwFGA7kXMpenExpgYuEtK0GTfPszstg\/5jn6Cnu7bioDpQMEnDqxbt1uj8fOYGKxs1gxDDh9GgcOBSW2\/LctllYwn9zyJXjt62X5PLVu2RKtWrULmUV1f9KHcuJqpXV3r\/tIdWr1GXeelJQDxp\/UONpiOzGIuVTYqsGQRbP4tGcl8qxdddJEh\/3l0efxxBz75xCGkjt693Vo\/ghwVZO3bt9eOMIsWAfKHHVWdib8HbDRPoMVA1+xMRybyHVKhyJiObA9gFBVhUEwMHioqEkeusWMd6NfPgXnzyoAkQGv8eQxG9najwNFfHH0O\/nKwUscy\/T1FJBCr29u8XwRAmFLMQiOjmO5N5pPwF0nqbypZxOnOO+8UnoKUEphHVWbvMkNiggd1HtyI6VkeNcJiAOmACMNSdRc81rBNn17+Qz1JgOTfGtmc5z\/\/2Y\/33mtbFr9irkuh85NQ6KZli39JF9Udx455dDLHj8eCRy7SzXQjqanmzE5L8\/Qh0Lxw+AV82\/ZbASS5yNUG80jkhKJFDjCtXXodeQmzZ2V+hzXbg\/TPmzdP8yqm74hMX+ivTo9a+0X1g6In8NVXX+2TQFl60DIMRHor+\/OZquu8rLZYGDKKRx4mJHrggQfQtGlTkTZQVqw3YqQL5SDBDZiR4QGBDM8eRYowC3pAZLhXibpwYQcsXeoQGzUhwbNhKYUQQKQEIi13UgGrLme52das6YDMTM88KSnG0gs3cSISy04l05HiShcAwsZ\/hw8vV+oS0CQtVrYO6SWIbN1ahOhoj2KYiljVzBmFKGEOZQCaGZBUB4AI8PTev9k98lnoG5+v1WY0nmMDzcFnL58\/15o+2bGaEJn5Uf\/3v\/\/hnnvu0Y7m8ggzatQoPP3000hPTxffESAoVasZ2CWAyATj9IwN5LsSMHu9VabUUD9LRxg7aONDk\/4jnI+M79mzp5Z3VK9EHbRnD3KcXFVOHCvqoW1MLoKzYmIw0u3GVLdbAAl\/XouJwXUHDwowkWfUfG+hq3i3G+rvXOg8TnBTL1hQfgziRbjZNv3yC\/7Vvr1Q0I7c4RaAtGOHGxs2OBAfX3504rGli6MLPij6oEI4emWVsf37O9CjRxHGjMk3VAwTvMY6xiIBCVjgXhDwEam0\/NzsZwE4ditRCe6U+Kw0I5k3mHS6\/mTmQHOQNr4w2FSHM33aBH5PKYRNTaosJRB9BnZ5hDECkN9++03UTZLNSAoJawnEymKw2scqgFCJ2nd0X1ze4XKfqX\/I\/UH7++bYWMzcvx8EBtlmNWuGlU2bItvl0pSo7EelJFucty\/H9Dh2TIg2VMBSGTtkyBGsXNkU3bsX4fzzf8XTjRsjK2GRGLes8AyMOO1mT58vm6L7Q0XoMewY5rT6FfHuW4TuY1nhMvQoau9DbyjKWHWCzz+PKTPvtULfvrvw8Z87insVdAN4o1kzcV9x7gXYGDMe63avQ7zb44hEejbGbBQKXfmZpGVd4jq8c8o7QmK6MuZKq4+uQj9DabGOSSDS5d1fbRd507IEBJ3ZVq9eDVUC4fFblnzgfGyUNqSE4\/SKPGqyKD0zwx5A7MqJKo8wf\/3rX4X3KiMk9UeY\/gX9MffCuUIROihmEFa5d4hTPpWU84rmabzf4MjHavcCLHUsxRT3FPH2P\/XYqRWUqJQ6uPFmeCUWTkApm2ZTIZWMdcC1CAI49u716Dvmzj2CjKQM\/NT0DY\/PBlvGdKS60pGTRXnI0yiFM3qjjysFw9t7JCS1BauMNdrNco6i9u0xp2VLIWnx+vIEsECYfTMR716NXu5e2hSvxbyGrUVbtb+pP3nk+CPYEx2DcTHLMcUdb7sEEjIa1cBAVQLh5dWSDPz7xRdfFDoS6jrYZLSvrH8rdSD8jhHajPiVaRaYLIpRwKxwp+pA2Nco9UVYA4idOVEJRDxT\/vTTT\/juu+9EYhc1wS\/Pr8nOjLIdQjmYHlZAqbNcaaiuMyouPZs4CbSCyEbPTWa3etxR7uIsv9PAQHzgBPWUHluJR0ciz8g7491IQzI6OxzCmUvqOTxBXtbzZtmhdzCbg\/ofcisbOcJFnLR69CSpPkBD4LnF7Rbf3GZD\/eK6vuhrALP8XrKu89JUiaqGTVcmJyrBiD4k7777rmDm1KlTRSEq2cjI5Ayg\/bAiIM8J9+p47Fjlq58I9OCnuadh2\/FtuKLBFWjarClGuUdp3c+KOUtsLlXZONI9skyHMlX0oaRCqYRv90VIRvuibWIjpsJjMlniWCLcz1e5V1lee5XVgfBCRnPQD4Y6HGnC7u9wYJVylOO4xl5P3PllOiJKK4TYxYuAr4ccFvqhiBnX8mOs8o4RAMnIEBppq8WUzBzJTj\/9Lxg9+gqhrOSP1RZI57Cy6UrEueM0Zy\/+vbLZSqEnmLl\/pgCQEbGxwm8jD8no\/2N\/3Hri1mp3SNPfq9E90YmOIMJj18yZ+\/3yiM5pk71lRuNWu7FxUAy2bt0m7iniSGZ1VVV9vwiA2AwggRzJAj3OUHQOlEyGFQ0TSkWpLwllHiO67JhHPweBg0e9UaPcwqGOfjF0qPNn1tzgcIj72jnVLfpu3nxEOKRFJJCqBwarVwh7AAk1J6o+9SClFB6H3n77bcFbxtGw3ol6hKmMPdxMX2D0QKlcpD6lFOW+dKHMYzS3HfOYzcEi8gQPbz0swzXL4w79Xi688AiWLHFUyqPVrmdldXPVpn4RV3bjp1FtfiDUdLtcLuHCTlMt3678XZq4KovEZpvN6mKsTfOY0UJnMyqA9X4RUilMcOHvDNxj7E1lI5XrM4D4Wz8RV3aLruxGDjdkarBZ2aUOpEOHDiLJix5ArETj+nvj2xVFW1vmMVPE8jgzaFCMiOqVStUVK2KExy4lD+nK\/49\/\/IKuXcsjlev7ESbiym71dWrez5IE4k\/xaT69cQ8ZF8NgOpmeXwKUWTSuv2tW1nFLzlub5rFCi4xGpkKVbc6clkLBSgX0iBGxGDLkMMaN2+sTYVwVSlR6zFj0ZDcsIRWEJ7vfElSB5qCFTfrwRFzZQ925FcdZBhDVnFuZy\/tLaSgBhMF2F154YdCXoA\/ImjVr0LVrV5GdO9RWm+apLC1UuhJI9PPQQ1J6SYbCJ6PjpvRJsTKfUfRmlJWB3j7+oj8DzUF\/Ga8ne8SVPQhem3W1BCCcROowrFZU0ytRmTOE4CFzq+oLA0XS5Jk9Kvu+Z55S\/oTa\/AU+1iUJJOLKHurT9x0XFIDI4CJ1Cqs6ECu5Q8wS9dpzy5FZKpvEt7IK75p+AhFXdvuegCUAMdJZBEuCvtocx8vcC8HOFelfsxyo6wBSs9zzvXpd56UlAFHrilr1OK1NDylCi70cqOuL3l5uVG62us5LSwBCFqma68qxrOZGq1GX6tFLNVHrpSIjpW+g\/lbvLhha9Amt1Tq5dtBilWbZr64v+mDvtyr713VeWgIQNaRfz0yrOpCqfAhW5qZSl2kEZs2aJeJ2ZELoyZMni8A+mZ9V7WNUx0Y9zvG6an8rdLBPsLSsWrVKOOFRga2mkuzfv3+NVAms64ve6nOqjn51nZeWAKQ6GFnd1+CDY4a0G264AbNnzxapBdQUdQRGutrr69jIvCb6\/jI1Yyj3YUaLfm6CH82wNFfzdztpsUJ\/XV\/0Vu5R3yfiym7MtXoLIOomDJRqUe9EJzc7zdJs+tSMoSxOq7RwblUCUsuM2kWLFfrrI4D440vEld2CK3s4HGHUBaD6tJgBQlUDSDC06PUxZrRbAYNQ+tR1AIm4sofy1KtAApH5ICsjvtt3K9Zmkm97NQYn0DHACEDsOjYEQ4s\/93+7aLHGPU8vY0cyUe3X0jRGed1k6U8rE\/jLCxdoDiaTkgmlIq7sVrhsrU+ljjBWnMOskVE9vYwAz0wpqgcQs\/5W7yQYWjinUSiBXbRYpVn2M3ZlzxB5Vaw0NX2C7M8yFVab0XiODTQHaZNlMVRHMtWKJa8fycpu9UkAlQIQVaFY2\/1DjBaKrO2rJsldtmyZVquGbDQKJFTn0ve3wvpgaQnkhFdZWqzQq+8TDhJIxJU9lCdfcYwlAAmkAwllA9lDemSWmuJAOOhAJICQh5Gs7KGvJEsAEvr0kZHhyIG6DiC16ZnUdV5GAKQ2raY6QktdX\/S1ic11nZemACLD8qkvoCekGqZfV7xQa9OCCQda6vqir03PoK7z0hRAVFOj3g8h2BwhtenBRWgJnQNy0TOnSPfu3UOfKDISBQUFmDRpEuqqLjHowlL6WBG7MpVF1lLd4UAk+ZO9z4ogPHPmTDBPS11rQQEIJY5PP\/0UdONm3Ehd8wOpaw+nNtMbSf5k39OpbIIn+ygJfibLRxhGfrJGzPDhw7VM6jzesFlNcxg8eZEREQ5EOFCbOWAKIKoPiMxDIT+jC3sEPGrz443QFuFA1XLAFECq9vKR2SMciHCgLnPAMoAYuV\/zxiOm3Lr8+CO0RzhQOQ78Pyt\/l+TGYCZoAAAAAElFTkSuQmCC","height":164,"width":272}}
%---
%[output:6679b005]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:9b4abf9a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:749b9ba6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:7e40b023]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:349fc67e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:8ba12aeb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:6fbf83c4]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_STN","value":"1.0709"}}
%---
%[output:63c6e05a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean","value":"0.9062"}}
%---
%[output:33ae1221]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_median","value":"0.9167"}}
%---
%[output:784e460e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAGXdJREFUeF7tXX1sVkWXPyCG8opkKRhKQYTsYoLxFbW+UVmNxmyQNaFmIWvTmrRBdIlgNFmklEIEdSnfiV+bbEOaUv8AFWWRxsS8\/iN\/bIGsCA0xRnRdEEQrXwZQqxLYnNl3nnee2\/s8d+7M3HPn3ufchIh0Ps75nTO\/nvk6M+zq1atXgT9GgBFgBAwQGMYEYoAaV2EEGAGBABMIOwIjwAgYI8AEYgwdV2QEGAEmEPYBRoARMEaACcQYOq7ICDACTCDsA4wAI2CMABOIMXRckRFgBJhA2AcYAUbAGAEmEGPouCIjwAgwgbAPMAKMgDECXhLIyZMnAf\/wxwgwAvEQWHtwJJw4Pwg3jq2C3iV3xKtsUDo1Ajl37hwsXboU2tvbYfr06QXRkTiWLVsGBw4cMFCHqzAClY3Ahdkb4MofxsOU6io4vOrexMFIhUC+\/PJLWLBggVCuu7u7iED2798PTU1NsGnTJpg0aVIsAJB0Xn31VaO6sTqqgMKMpVsjU+G5eM9p+ObcIEyePDmfBIKRx9atW2HOnDmwZs0a2LhxYyiBbN++He65555YVjx27Bj09PRAS0sLTJ06NVZdLlyMAGPp1iOo8JS\/gE3Gj4nGqUQgKChGIa2trSUJ5LnnnoNHH30UampqtPUaHByE77\/\/HiZMmACjRo3SrscFhyLAWLr1Cgo80fcPHjwolgAqnkDQfBhJNDc3a1vy119\/hdOnT8MNN9wAI0eO1K7HBYciwFi69QoKPN98800RgeNX8QSCayB1dXWxIpBffvkFBgYGRJ2qqiq3HlBhrTGWbg1OgSdGIO+\/\/75YB6x4AjEBAMPE7777DiZOnMgEYun\/jKUlgIHqVHjyGshfdmGYQNw6cNzWqBw+rlxZLU+FZ8UQSClHsAGAykhZdeI4cjOWcdCKLkuFp834idZiaInUdmGYQEzMRVeHyuHpNEq3Jyo8mUB4CpOup\/+ldyqH90JZAiGo8GQCYQIhcOfoLqgcPlqSfJSgwpMJhAnEixFD5fBeKEsgBBWeTCBMIATuHN0FlcNHS5KPElR4MoEwgXgxYqgc3gtlCYSgwpMJhAmEwJ2ju6By+GhJ8lGCCk8mECYQL0YMlcN7oSyBEFR4MoEwgRC4c3QXVA4fLUk+SlDhyQTCBOLFiKFyeC+UJRCCCk8mECYQAneO7oLK4aMlyUcJKjyZQJhAvBgxVA7vhbIEQlDhyQTCBELgztFdUDl8tCT5KEGFZ6YJBPOdLly4EPr7+2Hu3Lmwfv360NSCGzZsgM7OTuEZ69atg4aGhoKX2ABAZaR8uHR5LRhLt1amwtNm\/Jho7PQ2LhIDJjOur6+HtrY2aGxsHJIYGRXEcl1dXXD27NkheVFtAKAykgnQWavDWLq1GBWeNuPHRGNnBCKjj+XLlwvSePvttwEzUeP\/qx8mU+7o6IAtW7YIApF\/r66uFsVsAKAykgnQWavDWLq1GBWeNuPHROMCgci3WjZv3gw333xzYSqCjc6cOVNEDHKQh3UUfCgKCaSvry90GiP7wpynwWmODQBURjIBOmt1GEu3FqPC02b8mGgsCEQd\/PggDU4\/amtrC9FDOTKQneoSiNoWvkIXfNpBAsDPOpiY010ddHh+IiNbeKb2rAMO\/hdffBFWr14tFj3Xrl0rnlSQT06qPy8VhehOYeQ6CS6cYqZqJKtZs2YVFlIlgaDp+FkHdw4ctyWKZwjiypTl8hR4pvqsA0YGK1asEOngx40bJ56exOkMfvjU5KJFi4asZwQNqrOIqkYgSCC4ayPXTdQ1EH7WId3hQvEMQboa0vZOgacXzzqoW6y66x\/qNEZu46qEg6SBn9yu5W1cWuc16Y1qzm4iWxbrUOGZyhqITwaxAYDKSD7hlZQsjKVbZKnwtBk\/Jho728Y16Tysjg0AVEZypavP7TCWbq1DhafN+DHRmAnEBLUKqEPl8BUApVCRCk8mEL4L48WYonJ4L5QlEIIKz1QIRL3DUgpLncNkLuxgAwCVkVzo6XsbjKVbC1HhaTN+TDQuTGHCzmSYNGhbxwYAKiPZ6piF+oylWytR4Wkzfkw0LloDCZ4mNWnQto4NAFRGstUxC\/UZS7dWosLTZvyYaMyLqCaoVUAdKoevACh5EZXSyDYMyk7vzlKMpTssK2oXRq6F9Pb2iqRAeNlt1apV0N7eXrgb4xba4taYQJJEV79tJhB9rHRKUuFpM3509AiWKZrCSPLAm7jz5s2Dnp4eWLlyJezZs6fk1XyTTsvVsQGAykiudfaxPcbSrVWo8LQZPyYaD1lElbdyMdmPJBAkFvnv5XKCmAgQrGMDAJWRXOjpexuMpVsLUeFpM35MNB6yiIoX3U6dOgVPPPEE7Ny5E55++ml45plnRJaxYHYxkw6j6tgAQGWkKB3y8HPG0q0VqfC0GT8mGofuwqg5ObDRYOJjk45069gAQGUkXV2yXI6xdGs9Kjxtxo+JxryNa4JaBdRx7fBz\/\/0QnDg\/CDeOrYLeJXdUAILFKrrGsxSA3hEIZlLavXu3mNb89ttvIhEy7sqMGTNmiA66zzqoEU4wUZENAFRGqgTvd43l7f+2D745NwhTqqvg8Kp7KwHCIh1d45kZAsEEyMePH4fff\/8dZsyYAZ999hk89NBDJd97iXrWQc3KHpY+kQnEj7Hl2uFlBILaMYFUJWZkm\/FjIlTkFAYH\/MDAAFx33XVw7bXXCjK55ZZbYNq0aUX96eZEDWYnCwptA4BrpzcBNC91GEu3lqTC02b8mGgcSSCXL18GPFSGCZa3bdsGEydOhMWLFw+JQOJkZT9y5Ah8\/PHHYlrEUxgTsyVfh8rhk9fEjx6o8PSOQHTh1yUQuU2M78GUS6rMzzroIp9MOXR4ftbBHbYUeKb2rIMLmOJMYeSLdfysgwvkk2mD4hmCZCT3s1UKPFN91sEF7DrPOqiLqNgnP+vgAnn3bVA8Q+Bean9bpMAztWcdXGUkU9sp96yDfIMGzR08pGYzh6OaZ\/rrpu4kYyzdYYktUeFpM35MNOaMZCaoVUAdKoevACiFilR4pkYgqCRnJKsUd47Wk8rhoyXJRwkqPFMlEB9MZQMAlZF8wClpGRhLtwhT4Wkzfkw0jjwHYtKoTR0bAKiMZKNfVuoylm4tRYWnzfgx0bhAILg7gg9q4+EuXNi86aabxKPa+FE96YB92QBAZSQToLNWh7F0azEqPG3Gj4nGgkDU8xj19fXQ1tYGBw8ehO7ubnECFXdN+vr6AA9\/4f2VJD8bAKiMlKT+vrTNWLq1BBWeNuPHRGNBILh4qmYcC95XCf7cpCPdOjYAUBlJV5csl2Ms3VqPCk+b8WOisVYEoh4\/5wjEBObs1aFy+OwhYyYxFZ6pEAhCwmsgZo6R11pUDp9X\/IJ6UeGZGoHoGvLKlSswbNgw8SeJzwYAKiMlobdvbTKWbi1ChafN+DHRWHsb9+LFi+J5h3379sHLL78MY8eONekvso4NAFRGilQiBwUYS7dGpMLTZvyYaFyWQDDa+OKLL2DHjh3w+eefi0ztDz\/8MAwfPtykL606NgBQGUlLkYwXYizdGpAKT5vxY6JxWQJ55ZVXYPz48TB\/\/nzo7+8X7ePzDkl+NgBQGSlJ\/dNoOyzhMWPp1hJUeNqMHxONyxIITls+\/PBD+Oijj2DcuHHwyCOPwP3332\/Sj3YdGwCojKStTEYKhiU8ZizdGo8KT5vxY6Kx1hrI1atX4cSJE\/DWW2\/Bjz\/+KN7JHT16tEl\/kXVsAKAyUqQSGSsQlvCYsXRrRCo8bcaPicZaBKI2jDlScQ0kbB1E91kH2R6eL8FPffHOBgAqI5kATVHH5dsrlY6la3tR4Wkzfkx0jk0g5TrRyUgm60tFOamyidnC67h8e0V1+H\/u+ryiH4VyYSEmkAgUdXOiYjPyaDzes\/npp584AnHhoQDg8u0V1eHv2Xyooh+FcmGeXBOIi5SGulnZ0RgYqTzwwAPijRmZYDkYmWzfvj32jg+VkVw4lO9thEUgKHMlPgrlwlZUvpnaFCYsQ3oc4HQJBBXcu3eviDrw0l4pAuFnHeKg774sOjw\/6+AOVwo8U3\/WwSaloe4UBqOPzs7OIsuo6yCSQbFAS0sLNDc3a1uRInW+tjAZL8hYujUgBZ4V8ayDapZyEcimTZugrq4OampqtC1JkTpfW5iMF2Qs3RqQAs\/UnnVQocJbua2trbBx48ZCMiE8\/9HV1QXV1dVlUdV91kE2Uo5AeA3ErQPHbY1qzh5XrqyWp8IztTUQNIxcB2lsbCxawOSMZFl1W3O5qRzeXMJs1aTCM1UCKZV5jDOSZctZXUhL5fAuZM1CG1R4pkogaAiMNl5\/\/fVCPlQ5LcFLdOqJ0aSMZgMAlZGS0t2ndhlLt9agwtNm\/JhoHHoSVc1Oho0Gn5806Ui3jg0AVEbS1SXL5RhLt9ajwtNm\/Jho7PQou4kAwTo2ACRlJJd3TFxgZNpGHD2SwtJU9qzXo8LTZvyYYMwEooGayzsmGt2VLRKHBIINxdGDyuFt8chKfSo8UycQXANZsWLFELtQPS5lA0BSRnJ5x8TW4eOQQLCvOHokhaWt\/lmtT4WnzfgxwbYoAlFPou7atUvcV8HFU3nLtqGhwaSPWHVsAKAyUiyFHBeOQwI2XVcCljb4xK1LhafN+ImrE5YfQiDygSnMQibvqfA2rh60NtMLvR7oSlE5PJ1G6fZEhWeqBKJeqLvzzjsLJ1I\/\/fRTkY1M5zSqrZlsAKAyUikdbaYXtri5rp82lq71Sbs9Kjxtxo8JRkMWUdVoA6MQuR5icrTcRCAbAKiMVEov2+mFTxFM2lia+I7PdajwtBk\/JvjxLowJagnV8SmCoXL4hKD0rlkqPJlA9u+HpqYmMIl4qIyUlHfaRDCuo5esY5mUjUzbpcIzNQJRb9Lilu3mzZvh+eefF+\/B1NbWFo62mwKoW88GACoj6eoSLOd6kKvtu45efMfS1AZp1aPC02b8mGAjpjDBbGS4bdvb21sgDb6N+1dobUjA9SBXDW4TvYQ5DpXDmzhtFutQ4ZkKgQS3afEuTE9PD6xcuRJGjRpVSIK8evXqsjlBdJ51COZfDd6zsQGAwkgmJKAObun8vucWpcAyi0RgKjMVnjbjx0Q3EYG4IhCdZx3UQ2ny0h5Ol+STmTYAUBjJ5De9CemYGNNlHQosXcrre1tUeNqMHxMMnRGIbk5UVciwBEY2AFAZKS7QYaRjMxWK279JeV+xNNHFhzpUeNqMHxOcCgSycOHCwgPaYQ1F3YXRzcquto0RSEdHB2zZsqUwNZIAUGVll48mTbx+BPznoj+aYGhUR761UjtmBHzS9iejNpKshA7PWdndIUyBZ+pZ2W3giksgpTLAU2dln9tzEk5duAw4kHtbJttAEKvuv+z6Hr67eFnUoexXV0iKLOK6suShHAWemc7KHmcKExZ5SCeRBEKVlV1GINj\/\/ufvyIOvOtGBIou4E0Ez0ggFnl5kZbexh84ianCHJ9ifzRyOap5pg1FW6jKWbi1FhafN+DHR2OlR9qhnHerr66GtrU2cMVE\/9dSpDQBURvJ9AdTEEYJ1qLB0IWsW2qDC02b8mODolEBMBMhiBJLFbdm4tqFy+LhyZbU8FZ5MIAnehQmLHEyiieC2rEkbvg8EKof3HQdX8lHhWfEEgoNx35Gv4N4\/\/h30Lom3qBllpLDIwUU04aINV47qqp0oLF31UyntUOFZ8QRiMxhVI8ndlRvHVhWIqNyBLnRk0+Pl2O5\/\/c+PYiz8\/d\/+TWzi83EQUTm8j7onIRMVnhVPIHe298I35wZh8uTJsQe0aiR5UGtKdVXsdkwcyIb4TPpLug6Vwyethy\/tU+FZ8QRiA0BYBKITWahrGFj+xPlBUCMX\/LeodQ6TOzK+OHeYHFQO7zMGLmWjwtNm\/Jjoy7swAKBGDwgiRkDB6UjeIowoZ6Fy+Cg58vJzKjyZQCx2YTAK+N\/Tl2DaDaOL1iHiRA8Yecj1DHX6k7cII2pgUjl8lBx5+TkVnkwgFgRSKkqIih5UckACwSkMRiFU6yc+DhIqh\/dR9yRkosKTCcSCQP7x1f8Wg3\/EiBFFC6dhBKKShpyyoOMgaUgSUddPoqKYJJwuzTapHD5NHSn7psKTCSSEQHQHb6lFVBlNSEJQt12DTlQq6oiKYiicURcHF7JQObwLWbPQBhWeTCAhBKI7eMO2caVzqZGFGnGoP5d\/DzsP4sMaiC4OLgYUlcO7kDULbVDhyQRSJgJRpxRhTlMqAlF3VsIijnLE4ZNzUpIYlcP7hG+SslDhyQQSYw0kGNKHGUmWkVEHRiLBz\/QEapIOl3bbVA6ftp5U\/VPhyQQSg0CCIb000uI9Z0S2L3VHRS6QRkUx0qEo1xuonDhOP1QOH0emLJelwrMiCATfmSn15m6cy3TBkF4aSaYp1F0gDXNMyvUGHwcGlcP7qHsSMlHhmXsCUdMZHj16FHbs2AHr168X78\/gFzVwy0UGagRy4Pgl0R5OWdTpi5yuREUYUXIk4WQ+tUnl8D7pnKQsVHjmnkDUV+4wT+TSpUuhvb0dpk+fLuwXdZmu3MAutYgadq4jiiAoFyyTdFzTtqkc3lS+rNWjwrMiCOTYsWOwfPly8eIdPieBfw8+LFXqWYd\/6jwiDovhF3wOAa\/w41F2JIyoJxrKtZM150xCXnR4ftbBHbIUeGb6WQddqDEC0SEQbK+lpQWam5t1m4a0nmjQFjBDBb\/99lvYtWsXzJs3DyZNmpQhyf0UlQLPTD\/roGu2qCmMDMHmz58Pd911F0yYMEG3afjXP18srHf8acZUEakM\/\/kMvPbIWO02uOD\/I3Dy5El44YUXYMmSJVBXV8ewWCJAgefAwAB88skn8N5774GaqNxS9LLVya\/zRy2iItDLli2DAwcOWOl9YfYGuPKH8YJAxvx5uVVbXJkRyBICd999N+C7SpiUK+mPnEBQIbmNW1tbC93d3YUFVKkskgj+sfnWHhxZWCv5j38YZtMU12UEMoUAEgcFeSAoqRBIpqzBwjICjEBJBJhA2DkYAUbAGAEmEGPouCIjwAjkhkDKHY9nM5dHQO58Yal169ZBQ0PDkAq4+L1gwQI4deqU+NmiRYvE+R3+9BDAM0\/BQ5N6Nf0ulQsCidrZ8dsE6UqnOjZK0tHRAVu2bIHq6uoiwZBkgtcO0pU8O71L8kWJwzYNsqPJUElzQSBRZ0uybKCkZUdi2LBhA3R1dYn7SPj4eWNjY+FksOxfPQCYtEx5ah8JeuvWrTBnzhxYs2YNbNy4cciuY5b1zQ2BlDvdmmUDJS27GllgX0ggs2bNKprG4J0l\/Pfe3l4hTqnt96RlzXL7GIW0trYygfhoxKjj8T7K7ItMOgQSlFWNWoJTHV\/08k0OJhDfLKLIw1MYc+PoTmHUHvI6GMxRjK6ZV8xyMYXhRdRoBy5VQmcRFacwa9euFZcbMe2CStgyj4u5BJVRkwnEcztHHY\/3XPxUxVO3cdVLWIgpfritq27j8hpIfHMxgcTHjGswAoxAzhHIxRQm5zZi9RgBbxFgAvHWNCwYI+A\/Akwg\/tuIJWQEvEWACcRb07BgjID\/CDCB+G8jIaF6WVCKHOdCW\/A0KbYRvDhneiGx1DF3PCI\/derU0Mt5GYGdxYxAgAkkIy4SPHsRltG+lCqyLGa+lzdog\/XVA2Vnz57VPnYtSSdIZkgenZ2dJW\/3ZgR2FpMJJB8+EHZ4S\/0NH4ww5s6dW3iwq9TBLzVyUP8u2wreiQkiif3j9f7bb78dfvjhB0FOsi6eFcEvKgKRRINlVRJSz6bgz+T5FHn0\/vrrrxf\/NnPmTHERkI\/Up+PnHIGkg3vsXsMiEDW\/BA5E\/ILv7eAAC7sgF0YG6mBX24sS1nQKE7yHI0+7Yn\/qxTNV9\/7+fmhqaipENnHkjNKDfx4fASaQ+JilUiNsDaRU8h8ZBeC1fEkg6hV99bc+Rir4fMNLL71UdAs3zsB0QSDljsSrRIMEItMPYNShXmPgKITeNZlA6DE36rHUFEZGHcGMYTLsLxeBqAPztddeK5puBAnk3XffhTFjxsDs2bOHyG9KINiQSowqIaokh+XklAwJRE1sxARi5E7OKjGBOIMy2YbCCEQSQDCCUCMQXDgttQaiEsiePXsKLwYG10DOnz8vkuJcc801Yp1i9OjRRcrqEMiVK1fg559\/Bvzv8OHDh7ShyoyNq1FGuQiEUwsk63dRrTOBRCHkyc\/LRSDPPvts0TpH8IHlcrswU6ZMEYut6tQguAuD7eEAP3PmDNx6660wY8aM2ASCMuBaDGbmOn78ODz44INw9OjRItKSayDYfzBLGnYo5cQ1ELmoGmeq5YkpcyUGE0hGzBm2BqLutKi7Fhgl4BfcAQlOC3TOgVy+fBneeOMN0RY+EH3hwgV46qmnYNiwvz7WpROBIIHs3r1bJGb+6quvoK+vDx577LGiTGdSHnVHCXdz2tvbYefOnSJXK5IObg\/jt3fv3sLUhtMKpOPITCDp4J6ZXjFVJKYyfPzxxwHJZNu2beLB85qamlg6IIG888478OSTT8LXX38Nhw4dMjpgxsmdY8GeeGEmkMQhznYHH3zwgSAL+cA2\/tbHKASjnzgfEghGErfddhucPn1akNC0adPiNCHKMoHEhizRCkwgicLLjUsEkEAwksGsZvzlBwEmkPzY0mtNLl26BIcPH4b77rvPazlZuHgIMIHEw4tLMwKMgILA\/wHlc0G5L1wEfAAAAABJRU5ErkJggg==","height":164,"width":272}}
%---
%[output:99f23bd4]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"fit_test","rows":2,"type":"double","value":[["0.0259"],["0.8489"]]}}
%---
%[output:03ccf105]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_STN","value":"NaN"}}
%---
%[output:29764ddb]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean","value":"0.7094"}}
%---
%[output:042c7f36]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_median","value":"0.7000"}}
%---
%[output:45a6b617]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:093f4a51]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAHopJREFUeF7tXX2sVbWWL4gRRiTKx8BFnkBm1GB8ovAmMHc0GE2ER4Jm9D0JaCCIPCIqJiJfFxXQAS4f1ygyYxi8AsaAijrI1URjYuSPAYmiEmMc0XFg5EO8oMQvUG9g5rd5PfT0dO\/u9nT3tJe1E6Kw2+61fl39ndXVdrXDyZMnTzJ6CAFCgBCwQKADEYgFalSFECAEEgSIQMgQCAFCwBoBIhBr6KgiIUAIEIGQDRAChIA1AkQg1tBRRUKAECACIRsgBAgBawS8EcixY8fYnDlzWEtLSyLs4MGDWXNzM+vevbu18FSRECAEaouANwL59ttv2YwZM1hDQwO7+OKLa6s1fZ0QIAScIOCNQD7\/\/HO2ePFi1tTURF6Hk66jRgiB2iPgjUBeeOEFNnfu3JLGS5YsYWPHjq09AiQBIUAIWCPgjUBECTGdmTx5Mps9ezYbPnx46dW+ffsY\/tBDCBAC9gj069eP4Y+PpyYEwgOq9fX1JS8ExDFz5ky2Y8cOH3rTNwiBdovAsGHD2PLly72QiDcCwRQGD6YtiIfMmjWLLVu2rBRQfffdd9n48eMTxS+88MJgOxcE98QTT5CcjnqI8HQE5F+b4Xhu2LChzLt3+5XTrXkjEHkZV46BcALxpbgtoHv27GHr169nEydOZAMGDLBtpvB6JKdbiGPB0\/c48kYguu70rbhOnrT3x48fZwcPHmR1dXWsc+fOts0UXo\/kdAtxLHj6HkdEIIZ2FoshkZyGHaspHgueRCCe5m625hWLIcUk55h\/\/ZC1HmPsdxd0Zi13X2XbNYXWiwVPIhAikLKBYLu03dbWxg4fPszOP\/\/8oKdakPOe5\/+LHfi+LdF7y7QrCyUC28ZDwjNrmZYIhAikZOO0tG073Nt3vaxlWiIQIpCS9ceytN2+h2tY2umWaYlAiEAqCCT0pe2whlj7lkZHELr3rtGhVRhDRH0G03wbgwgFNv71798\/12YkbAx888032T333FOBJvb\/LFq0KNk3I5\/C5kcadu3aldSbOnVqcrxBloOfoeLvxT1FYloIyDFp0iR24MCBsrY4jmhXdwYLemPPhyxH2jd1suL90qVL2erVq1nfvn3Z2rVrExzS2tPJqrMJ3XtDc9cWJwLRQlRe4EwhEBNY0giEE0Rra2tp4Ijtoh425c2bN4916dKl4pPyewzEESNGsL179yZlsasZA2br1q1s+vTpZUTFy15yySVs4cKFbP78+ck30sgM7fHBl0Zk8jdFkkmTFXUgH8qKZbZs2VKhw5QpU7Sy6ghC996kX\/OUJQLJg5JQJiYCgTG99tpr7J133kl+lfmvr\/jLz3\/BoSIG2s6dO9mYMWOSXbbwQPBeTATFp1NiGziC0KdPnwoP5O2332aXXnope+qpp5QeiPhrK\/46p3UJ94ref\/99NnLkyOSXHHKsWbOG3X777ey5555LiAREgbZBNEOGDCnzjtI8K7QD\/dHuxx9\/XOGB8PY4aaHtrNPk\/Dsoxz05eB0rV65kIArgKOswevTohGy4J6eSVUcQuveG5q4tTgSihShcDwT7J7767njq\/gkY08aNG1ljY2PiMvNfYgw4\/JLjJDT\/VbzrrruSs0n4pUaWOHEAcJeeDzKUEdtIc\/uBXNYURj4fleWNiN9YtWpV2eCDXvfeey97+eWXywgEg\/Hmm2+uIBDuSfCeFWU8cuRIyWOQTYOnpFB5KPJUhmMmkoD4HUz5RAKBDhMmTGDbt28vIxBZVh1B6N4bmru2OBGIFqJwCeTKf9nO\/vfb4+yi7p3ZRw\/+Y4Um4q8mXmLggTiefvrpUmpJ\/Du8DBxiFAdg2i8opgD45cV7HtewjYGIAovkJKe5xHQEHhH\/xZcJxMYDAbkgLgHdH3zwQXbfffclXhp\/5FiJKINIzPLUS5ZVJhDyQAwHnKo4QMYjziF9M6etGiFNYbgHAl3SCETngXAc+FSATwFMPBAeh5ADjzoPhMcp4Amp2uCBxnHjxpUFc0XPxWUMBPKm6SLKqordpMkqtkcxENtRJ9RLC1QRgVSCWy0mcgxEFb\/AVxHzQHoFMYbACSRPDARlrr\/++lyrMOKAUsVixEAnphOIr4gPdBBlcr0KkzbgQRBIgoUVIzFew4kFMqpkBTnSKowD4kAT3E1F8Ounn34iD0SDqwsC0QX7HHUtNeMJAZ1N6N67FtNrDERchpPX2n0rbgtkSFMYnQ5yDERXnt6Hj4BunOjeu9bQG4GIrqEqas8VRzDrpptuSpYFQ3xAIF9\/\/TXr3bu3cu+CS5mxXInIPO1EdYlq3G3xcfLMM89UbPKDXWIZHqlBfdmMNwLh80Cx+8TlMHFPAKL7GDghPr\/88gvDxqhevXqxc845p1ARMd++\/\/77kxUCHKCihxDYv39\/QhCPPfZYEgsSn2effTbZmIen3RGIqGiWB4LlxKFDhwbrgSCYdujQoUS+ojOS4TTuQw89RImmiTfKEMCPCZaZZS8dHsirr76a5Ow9YwnEl+K2NukzBgIZbfOBQM6jR4+ynj17sk6dOtmqW3g9ktMcYsoHosDMd\/DHvNtO1fBNICSnLQJu68XS777HkbcYiK47fSuukyftfSyGRHLa9rC6Xix4+h5HSgLhx6JXrFjBcJqRb6ABtOLGHZdd5FtxW9ljMSSS81QP684L5bWDWPD0PY4qCASbvWbMmMEaGhqSm61wEhM77\/g2ZQRAt23blhzQUh3BztshcjnfitvKGYshkZynelh3XiivHcSCp+9xpCSQrPwJWYee8naGqpxvxW1ljcWQSM5yDwR\/U50XymsHseDpexwppzD82DJWRHr06JFkecJ0Bg\/2++uOM+ftFLGcb8VtZESdWAyJ5LTtYYqBmCCXGUSVN38VFf+AwEQgJt2mL0sEosfIpEQsePoeR7QKY2JF5IEYoqUvHsvAjEVOIhDKyq4fdTlKxGLwJGeOzjQoQgRCBGJgLulFaWA6gbHUSCx4EoEQgTix\/FgMnuR00t2lRmpOIPJdHSr1igim+lbcttvI4G2RC391I2vTWSz97nscKYOoPL9jfX19Zup6E1MSL9JBPfnQnG\/FTWQXy8ZiSCSneQ9nbTqLBU\/f4yh1FUbckSrfKGbeNSzJ4s2zkEFJLBE3NzcnVwjQMq4Notl1YjH4kOTMSlIdkpxZPR8Mgbg36dMtQkmeLZxvh\/etuK1+sRgSyWnbw2FNtUzP8vgeR173gYjTGJrCuDVwuTUiEHN8Q4yBmJ7lCYpAxAHPU\/\/jEh4ctKtmWsMDtTigh7T34hSGcqKaG76qBgjEV+7WaiT2Keefmz9NbvKrO68T+4+pv68Qe\/iKD5OLuvp268Ten\/MPZe99yil++J9Xf5zIjEeWSVYgqJyonDxwEhfXA\/JrB3EpcLWncVVBWsqJWs0wrKzrM3drNZL7lHPM+n3swPdtCUG0TOxXIfZfXvmaHfyhLfl3+b1POW3xDConqnjqFhf8cAIR71iVryDMUly+BxUXGeEuVu7JcAKhnKi25lNez2fu1mok9ikn90ASj\/eBq4zE9imnkWBC4eByomKlBPeF3nHHHWzTpk0MFzDj5nBMO1TXGGYpTsu4tmZhV49iIHa4pdWKBc+gYiBibIIDK1867KqbfCtuK3cshpQmp2lU3xanvPVixzOvnr7K+R5HXldhskD0rbhth8Zi8CCK\/2n9kQ3s1ZW13H3aXRej+r+7oHMSoMN\/xTK22NjUKxJPl2RZpJw2uKXV8T2OjAgEgaTNmzcn05pff\/2VITaCVZlu3bpVjYFvxW0FjsWQ0pb\/xM1SwACrDhd171xVti5bLFGvSDxNl0Cz9ChSzmrwk+v6HkdGBIJky7is+bfffmODBg1in3zyCbvuuuuc5Eb1rbhtp8ViSH984r3Eu8CdMGmp\/LJ2XtriY1qvSDxd6leknKaYheTJGxMIbmU799xz2dlnn52QyWWXXcYGDhxYNQZEIFVDWNZALAZPcrrtd9\/jyIhA2traWEtLS7L0um7dOlZXV8emTZtGHohbG3DSGg1MJzCWGokFz6AJxG2XlLfmW3FbXWIxJJLTtofV9WLB0\/c4MvJA3HYJEUiReMZi8LZyulxhydMPtnLmadtlGSIQykjmxJ5iMXhbOeUVlqIJxVZOJ51p0EjNCYQykmX3ViyG1N7llFdYXC7ZqiwgFjxrTiAAr4iMZDoS9a24Tp6097EYUmhypnkIruR0uWRLBJJ\/dHjLSKYTiQhEh5DZewxMDKrWY6ymO0251GkegisCMUMnvXTRROdKzrR2fI8jb0FUeWokn6nxrbhtR4Zm8Fme0h8a30uOr9dypymXL81DCA1PE6IrOu5iY6O+x5E3AsHJ3gEDBiRJmrGjld+3KycUkjOV2YBYZJ3QDL5IAvExQELDMw\/R8bQAOAaAJwSC5nYQBIHwAY4zL\/AU+vfvn1yqjcfFlQ48xjJu3LiKjGREIG7oT5zCoEWbm+mLDkxCrtAIJE\/si2cuA3HwxwZfNz1d3krNCUQMoN54441szpw5bOfOnWzt2rXJDlQkBqo2IxkIavHixaypqYmyshdhRY4GZtGBSU4gIcVq8hCImJgoFOIIxgMRM5Eh45iYSQxCyu9N7T\/tugjOnJQT1RRRdXn8sseSE\/Xqxz8upRrU5f10g455KzHgGUROVJ0HwrOUNTY2Gp+BUXkeMnPi7xMnTmQTJkww72UPNWLIjQkYYpITuUi\/\/fWspPeycpUiGfK\/39zHQy9XfoLjOX87Y9\/8zJLEzLWSJQ2AYHKiFhEDQZs8ryq\/C0YEgnKiuh0XMeTw5BnH\/\/ZvGHvlL79nnTufjimIaIgxB9Ncprwd\/i0kT1JlZNehz\/H80\/OtpZUtW1l037J9H1xOVJ0iJ06cYB06dEj+ZD1yPlReVgyY+g7+6HTLMxdOM3jbtl3WiyE4yYO0\/BoFjqe8+uMiFlNtQJjjOW3L4VLm9tDiH7Af3+PIahn3hx9+YLjeYfv27ezRRx9lF1xwQdW271txW4FjGJjQLUQ5VcSAtIt4EPvgBFLtYFf1rUxCpkvUWXiatmVre3nq+R5HuQkE3sZnn32WXEn56aefJpnaR44cyTp27JhHL20Z34prBUopEOLAVIkaopwqYlDJ6cLj0PWvKUll4Wnalk62at77Hke5CeTxxx9nPXv2ZLfccgvbtWtXoiPfBFaNwryub8VtZQ5xYIZOIJwQeP5VyMvd\/1rhaUpSeTwQUS9b+6q2nu9xlJtAMG1544032FtvvcV69OjBRo8eza655ppq9S3V9624reC1MnhTeUOSk\/9CQwd512ZIckI+OgtjZmm5CYQ3e\/LkSfbVV1+x559\/nh09ejS5J7dr165mX1WUJgKpGsKyBkIamBiU\/\/nfRxP5qiEQH7EGk7MwbnvMTWu+x5ExgYhqIkcqYiAu4iC+FbftrpAGZpYOocmZdc2E6v4alW6uYw0qQspzFibk1Tff46gqArEdhKp6vhW3lT20gZmmRx45ffyic\/nSBqYJKZjGLXR9bPLtPHjqvufjve9xRARi2KuxGFIeOU0GkCFMuYvnub8md2OGBU0IKQ+ehp8vpHjNCYRSGmb3ayyGlEfOogjExLPJI2chI82w0VjkrDmBAFdKaZhuXbEYUh45TX6BTcabCTGZyNle7\/A1wVZXNggCgZBpp2Z1Cti+9624rZx5DN62bZf1aimnCTHJcoreC\/DA9Zxi4p5aXQheSzxN7ML3OKIYiEnvBLpFXKVCDAbPyaJXF8Za7r4q2couei\/Qi28+4yTCdfWdBSwGPIFNUASCE7SzZs1iy5YtKyUTwv6P5ubmUiIgw\/GXWty34rZyx2JIMcipOkwnEgj3NnhfiSRCBKK2YN\/jKNUDUaUdhMjVZiTTJRSilIa21FZerz0QiHjaVdzNCk2JQAInkLTMY9VkJON5RqA6T5HIYfDNnLbDNIaBCd1UcpqsjtjiY1JP5YGIu1b\/6e\/OT6Y2eOTt8Pg3n8fpY+l33+MoMwYCb+PJJ58sDXa+xItDdLNnzzaxlSQou2bNGjZq1Ci2YMGC0rSICMQIxtyFVQZvsjqS+0OWBTlRIBcIHt1x\/qzt8JYiGFUjAlHDpQ2iitnJ0IR8n4tRLzCWXOkgxlVkAqGcqKaIqsvD4OWcqDwrFx+wLr7EEwwjxZ+c6UtMPoxviWV4ljFkI9t8ex3r3bt3KUUm5Nyx91SekGH9u5baLUL+vBio8Mxb11e5IHKiFq2sjkDwfcqJWn0vuMqJinylB39oS80BOmb9vlJCZDmfKX\/HtYG3wcvwdnGeas3ITqxXr17s3te\/K30L38SlWGKd6lGxb8EVnvYS6GsGkxNVL6p9CR2BLF++nA0dOpT16VOb5Lk6zWLINQodZDltc4Lq8pGKXoaYI1T0IjimCHyijCjLhgl\/zw4dOpT097WrPi0t24orMCHkHo2h34PLiYoYyNy5cyvGVDWXS+kIhFZhdBSW7708Z7eNf5hsCoNkYvIglaQgEXE5FlOUf7uxJ6urq2Mh37dCMRDDGIi43PrKK6+wESNGJBnIxCsq85lyeSkiEBvUzOuk7fBESy5XL0SCQdvyhi\/8m0wa4t\/5FAUEEvIxeSIQCwJZuHAhmz9\/fpKFbM+ePcnKSzXLuFnDwPfyk\/mQPFUjFkPKI6eLZV15f4ZMGDJ58Pd8ioIYCOIiRCC2Fllez\/c40m4kq6+vZ0OGDCmtnHzwwQdJNjLXu1F9K27bXXkGpm3bLuvl2QeimtbkIRXR60jzODhBqLagi16QTzzz6JbWBz7lrMYOfI+jzGVc0duAF8LjIUXEKXwrbttJsRhSnn0gqviGLlYi7sdIwxBeRxaBiFMon3jqdMuyCZ9y2tom6vkeR9p9INUoY1LXt+ImsollYzGkLA8kKw4ikop48hV1xJOxKvzE2+p5PIT\/m3iitlYEYhoQjrHffY8jIhBDJomVQGzcd9XJWB1c3PvghMPJSt6KznN7xIqnDodavQ+CQMSsZFiyXbFiBXvggQeS+2D69u1bcY7FBVi+FbeVOUaDx\/Ioz4ouBjH5ORPZ61CdP8marvB3aV4G3svLu\/wwXIx4hrxa5HscVXggcjYyLNu2tLSUSKPa07hphuhb8TOJQMTNYPIgV8U05GkHx0r17+Kp2DxTBLkMEYitJarr+R5HFQQiL9Ni38b69evZvHnzkrMKtIx7nB08eDCqZUeRQMTYBGIRqmXYtPiGmNxHjHeo9pWkTZnkfSM4H8M3koX8yx4L0RGBbNjg9MpMt\/we5z4QTiAyFqo9GrqpisqDQXBVzleatuIhExZtJHNroUQgRCDWFiX+6m+aPCjxlKZtOVw62Wrd8F93k8orK6IHIyf4SZvOyB4IbSSrplcq6wZBIJMnTy5doK1Sr5qzMBQDcWswYmtyOkDc+IYTrS4eMf5hGvfI+n4sU4NY5Kw5gbgwNps2fCtus6wJvUI2pDznUkz7Jm1vh6vzNCHjKWIVi5y+x5HXfSDi6V55N6tvxW13JcZiSHmDo+IKS9q2dJAFX63hhJKXQHREHQuescjpexx5IxCs5ixevJg1NTWx3bt3s40bN7LGxsZSFirfiudZclT9WodsSOJg1e0azfJEZK8jK9ah82h0RB0ynuSB6HqXMW8EIu4fwV6TGTNmsIaGhuS6CDy+CUQPjbpELQ1etzSq8iBs9FRlPDchXFFOfB9khkfltdQSTxNsYpHT9zjySiBiSgAEapEeADlGRAKhnKjpZs2XY7H0iSTEeORdpiaDgpdFezzYKiY5zsp3yuuqyqjkTJMLA1PO3ZrnuzZ6VlNHJWc17RVRt13nRIUHkodAACzlRFWbF88jirc8t6iYk9R0xYWTBdrDhi7kIU1rW853yiVU5URVyZk2YFS5RrPyrBYx8PK0STlR1Sh59UC2bduWxD2ypjCUEzWPOZ8uw3+t06YvcjxDbj0r32havlOxjTxlsjRS5Rqttk0zBPOVppyoNSaQ0IKo+cymslRIc+G0xD4gDWzQ4h6JeEGTrd5F1QsJzywdY5Gz3cZA0Dl8GVd1ote34rYDIiRDSluq5dnP\/9D4XkIivq+BNME2JDyJQEx67lRZb1MYnWhEIDqEKt\/LG8fEEpiaEIGYY5pWIxai8z2OiEAMbSwWQ4KcRCCGnZtRPJZ+JwKhw3ROrB4GDw+l9dip5vLuHHXycYNGYhmYschJBEIEYjD80ovGYvAkp5PuLjVCBEIE4sSiaGA6gbHUSCx4EoEQgTix\/FgMnuR00t3kgfhmTttuI4O3RU5dj\/B0i6fvcUSrMIb9RwZvCJimOOHpFk8iEJrCOLEoGphOYKQYiAZG8kAM7YwGpiFg5IG4BUzTGnkg5IE4MTgiOicwkgdCHsiZaUhEIGdmv7d7DwQXU8nZyNDVvhW3NS8amLbI0SqMW+TUrfkeR15jIDjSP2nSpETztWvXltIZEoG4Ny0iOreYxoJnuyUQeB5r1qxho0aNYgsWLGDLli2LkkCQVQ1XfSJr2oABA9xaqcPWSE6HYDKWZNOLod\/bLYHw7oQXMmvWrFQCQU7UYcOGue19h63t37+fzZw5k5GcbkAlPN3gyFvheMrXprj9yunWvE5h8Nk0Atm3b18yMHfs2FGUrtQuIXBGIIAfYKQG7devX+H6FkYgS5cuZatXr04UENkwjUBQDiSCP\/QQAoSAPQIgDh\/kAQkLI5A09bMIxB4yqkkIEAK1QIAIpBao0zcJgXaCgHcCaSe4kRqEACFQiylMGupY5sVtdbt27UqKLFmyhI0dOzbYTkrbEBeCwFmXmIcgnyxDyFhC1lhsE3fXzJkzh7W0tFTEHovq92A8EARdsa8CpME3nK1YsaJ09WVRANi0m7UhzqY9l3V09++4\/JaLtkLGkusXi22Ktz9iPwjkbm5uZt27d3fRVco2giEQUTrOpOPGjQuOQHQb4grrqZwN6y4xz9mMl2KhY6kCIWTbFOUFgWzcuDG5CbJLly6F9WeQBCL+ihbJntWgGupqku4O4mp0LqpuqFiq9A3dNsVpjI\/NZMERSOjzYW5UoRo9EUhRNHcqFqI6CFrcF+1b5nGb2bNnF+rF14xAVBvNQmR3mw1x9t1efc2YpjChk7HYGyHaZpa1cE+kvr6+0MWImhGIrDw6CIeV5s2bV+icrfoheqqFUD2Q2IKoIWMpElwMtokfDzx8IUJ15syV\/fN2giAQefmJC+djDmcLaKgEAn2yLjG31bfIeiFjGZNtntHLuEUaKLVNCBACxSAQhAdSjGrUKiFACBSNABFI0QhT+4RAO0aACKQddy6pRggUjQARSNEIU\/uEQDtGgAgkgs4VD8dxcadOncqwSSjPo1pJkA8r2h7AEzeuibKI50fyyEhl4kSACCSCfhM3h+Fcg8kuQ152+PDhJcKR64sHr44cOaLMWauCiZOOTGZ8813oJ6oj6PrgRSQCCb6LTu3r2LZtW9nBKPEXXvYwxowZUyqrqguVRc9B\/P+8Oxjx\/QMHDrArr7ySffPNNwk58bp9+\/ZNUOWnq9MgFnf5iiTEM4vzenw\/ED8gdt555yVpMgcPHlz4adMIzKOmIhKB1BT+fB9XeSDimQwMRDwYxKJ3gQGG\/BC67czydENsTyeh7RRGPC2KbyxatCi5KgOPuINS1B25YsaPH1\/KFWMip04Pem+HABGIHW5ea6liIGnTA\/G4OScQMS2C+KsPT+Xhhx9mjzzySBnJmAxMFwSSddxcJBoQiJjjIrbzKV6NxtPHiEA8AV3NZ9KmMNzr4El5MKUQ3f4sD0QcmCtXriybbsgE8tJLL7Fu3bqxG264oUINWwLh06i5c+cmbYqEKJIc3vEpGQhEzHFBBFKNVbmpSwTiBsdCW1ERCCcA2YOQE96kxUBEAtmyZUty85oYx+DTnu+++y65UfCss85iiFN07dq1TNc8BHLixAn2888\/M\/y3Y8eOFW2IMqNx0cvI8kB8Zd0qtHMjb5wIJIIOzPJApk+fXhbnkK82zFqFueiii5Jgqzg1kFdh0B4G+OHDh9nll1\/OBg0aZEwgkAGxGFxrunfvXnbttdey3bt3l5EWj4Hg+5xAMLVBPTxcTsRAeFDVZKoVQTdHKSIRSATdpoqBiCst4qoFvAQ88gqIPC3Isw+kra2NrVq1KmkLl0t\/\/\/33bMqUKaxDhw4l1PJ4ICCQzZs3Jxerf\/HFF8mK0q233lqWAJjLI64oYTWnoaGBbdq0iTU1NSWkwy8r27p1a2lqU2TKvgjMo6YiEoHUFP6wP45pDTJ833bbbQxksm7dOjZhwgTWp08fI8FBIC+++CK788472Zdffsk+\/PBDqyQ3vvJ8Gil3hhcmAjnDDSBL\/ddffz0hi6FDhybF8KsPLwTej8kDAoEnccUVV7DW1taEhAYOHGjSRFKWCMQYssIrEIEUDjF9AAQCT4bv8yBE2g8CRCDtpy+D1eTHH39kH330Ebv66quDlZEEs0OACMQON6pFCBAC\/4\/A\/wG7CL4ixKuOngAAAABJRU5ErkJggg==","height":164,"width":272}}
%---
%[output:2f2715da]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"fit_test","rows":2,"type":"double","value":[["0.0481"],["0.9701"]]}}
%---
%[output:5909ecc4]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_STN","value":"NaN"}}
%---
%[output:4717ed11]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean","value":"3.2208"}}
%---
%[output:78ec53d4]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_median","value":"3.3945"}}
%---
%[output:5e83376a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:8654cf04]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAGzNJREFUeF7tXWuMVcWWXkw60hAHpYXwEIHJBP3jhAvMDciYEMPgK9NNruYGmh\/NIGGIo9JRBBowtm3kDYmtXh1EbCEoECMj05rAkCBoBjDXRggSIzoZkJeK4AMfoK1M1sY61tln711rv2rv0+fbibmXPlW1V3211rdXrapa1e3SpUuXCA8QAAJAIAIC3UAgEVBDFSAABBwEQCBQBCAABCIjAAKJDB0qAgEgAAKBDgABIBAZARBIZOhQEQgAARAIdAAIAIHICFgnkHPnztH06dPp4MGDjtBLliyhSZMmRe4AKgIBIJAdAtYJZNmyZTR06FCHND7++GOaNm0arVy5ksaMGZMdCngzEAACkRCwTiC6lD\/++CM1NTVRfX09CCTS8KESEMgWgUwJhD2QxYsX06pVq6impiZbJPB2IAAEQiOQGYFwLGT27Nm0YMECGjZsWJHgJ06cIP4PDxAAAuERGDRoEPF\/Np5MCCTI82DimDNnDr377rs2+o93AIEuh8Do0aNpxYoVVkjEOoEweaxbt44WLlxIPXr0KBm8ffv20ZQpUxwArr322i43uFE7xITa2toKXDQAgUmpNilMXnnlFStxRasEooKm7e3tRT3XO6sIxBYAUQ3adr2jR486xDt16lRnFQsPETAp1QLb9mOVQCRKbxsAiUx5KHPhwgU6ffo0DRgwgKqrq\/MgUuYyABMQiO8UBh5IMTQwllJjASYgEBCI8DsOYwGBSFTFtgePKYxkVAxlbCw7d3Z20pdffklXX301pjC\/jUelYhK0TAsC+W0VplymMFh2ToCB0UQoBIKWaUEgZUYgWHYOpfsoHBMB0zItCKRMCaRcPKaY+ovqGSNgIgjT70mLjxhITERtD5gu7ubNm2nIkCGiDUO8gW\/79u10\/\/33F\/WY25g\/f77zt5kzZ9K8efNI368zfPhwWrt2bdFZJT5RPW7cuJL3clu8N4PbcD8KJ\/67SuGgv5v\/PnDgQGprays52sDvW716tefv3O7u3bs93xmmb7p8tbW1tHTpUmejo1cbXirDmC1atMjZp6OOZnj9zT1+buz5d1N\/eaOl3wfLtj6CQMqYQMKI7kUg7l3BihiOHTvmNM0pF3QD1YnFrcBKcRUJ6bLxuaeWlhZqbm52jNJtaFyWZeH\/7rzzzqJu6e93y6uM2+udYfo2Y8aMgnx8qJPb5WfkyJFFu6b9iFPluDlz5kyBAL3+pnfMTz4uowjRa9e2iSBMv4fRGUlZEIgEpYAycQeM67\/xxhu0a9cuOnXqVOHrrCdeUl4Ai8GG2NHRQfyV5B2p7IHw75wWQe3wVcatt8Ffrf79+5d4IF4ezXvvvUe33Xab8yXlNtasWUOzZs2iI0eOUM+ePWnLli1FHogiCK5z6NChEm\/ATV5uz0l\/h\/t4g16WCeypp54iNnh+Dhw44Mjj54GE6Zv+Xh4TJlF3ois\/j2\/nzp10ww030HPPPVfwQLz+FqRqqm1+r\/IqvTwYk76Zfo+p7iXVQSAxETUNWO1f3qfjX12g63pXU\/t9Izxd+40bNzouMyuM+lKz0appgvoS3XvvvbR8+XLnS66+lKxsrHRq6qB\/7fU2gqYXLJT++zPPPFNEIEomlXJB\/xLrSn727FlPY\/YiEOXhqHf7TcXcBOL2XoKmMArsMH3j9tR46KRiwi\/KFMZLPv09XifWTfpm+j2muoNAkgbQNGB\/eGIvfXruAg2uqaYDj9zkSSD6146Nl4njhRdeKHgUXIm9DD5g+NprrznegJqfKwJxf7X468nKqObkfjEQblvPEsf\/dhOI8kCUQekEorLKsfekHo5xMKFx3ILlfuSRR2jv3r0F78fPq2CC0tvjdvjR+6Y8EEVmJgIJ0zc\/knC3oWIUenzIRCC8d0N5iXqMxd22Pk3kmNCoUaMKHwzGwqRvpt+T1n94IDERNQ2Y8kAcl9uHQEweiBLR7errbq\/JA\/EyNL+McCoG4I6BKDn8YgF+xhwUAwkiNmUwppiA1xQmbN\/chszvDpMxz0Qg7pw3fm3rGHpN7Uz6Zvo9prrDA0kawLgD5o6BeMUvWGb+as2dO5c2bNhQ4oFIYiBcZvz48UUxEH3lQeHC79fbC7MKE+QN6O\/SE2n7xRv0cTKtSigCcU+nOO6jP35946kX5+Z1e1Hs+Xi14ZW\/NyyBRMXepG+m35PWf3ggMRGNO2ASA4opYkVV55QHTLZdNUWmSd9MvyetDCCQmIjGHTAQSMwBqLDqJn0z\/Z40XCCQmIjaHrCY4qJ6mSNg0jfT70l3HwQSE1E1YI2NjcSHnPAAgTQROHnypJMzGDtRfVC2zaBxBxunceMiiPphEcBp3ADEyo1AuCs28oFwQqGvv\/6a+vTpQ1VVVWF1rkuWr1RMkA+kixGIDetERrJSlIFJKSa2P8CIgdiw\/gTeAWMBgUjUqEsQSNC5AffVDu6NSrYBkAxKHsqAQEAgEj20bT+JeyBq16DXEWsGIOhKS\/7dNgCSQclDGRAICESih7btJ1ECUWczeGsxP16JZUwXatsGQDIoeSgDAgGBSPTQtv2UEIjX7XH66UFJJ9gL8SMQdxYq\/VwEPBB\/dEEgIBCJ7WVKIIo8+Bix8h7U31h4lebN1JEgAtHrqoQ3\/C51QEnfmDVx4kQnCQ4eIiaQzz77jPr16+d5p3AlYgRMiked9YOTTQVtNEtaT4o8EP3YtX4Yye\/vfsJICUSR09ixYwvZn\/RTipzLoqGhIek+l2V7Fy9eJE6Z17dvX+revXtZ9iFpoYFJMaLr1693UjDyYyvJd8kUhg3YnZHJK1dCkDKYpjBcl3NNcDyEj6hzli2VL0ERCCfP4WQq8EAuI81k+\/nnnzt44G5cYOJlf+yBbN26lVpbW7MhED2Hph9BeOWHcJd1E4ieoMYdY0EMRPZdRgwEMRCJpmQaA5EImHYZ2wCk3Z+k2geBgEAkumTbfkqmMO5VEiW0xPOQdNBUxjYAJnny8jsIBAQi0UXb9lMSRJ09ezYtWLCgKHV\/2BiIpKN+ZWwDEEdWm3VBICAQib7Zth\/fVZgdO3Z4XhWQdqo42wBIBiUPZUAgIBCJHtq2nyIC0ZdV+VYutUKyf\/9+2rRpU8kVh5IOhS1jG4Cw8mVVHgQCApHonm37KYmB6Hs+2AtRd3faWle2DYBkUPJQBgQCApHooW37SfQsjKSDpjK2ATDJk5ffQSAgEIku2rYfEIhkVHJQBgQCApGoYSYEktQGMkkHTWVsA2CSJy+\/g0BAIBJdtG0\/8EAko5KDMiAQEIhEDXNHIJ2dnc6FzjfddBMdOXKEfv75Z7r99tupW7dukv6ELmMbgNACZlQBBAICkaiebfsxeiB84I0Pcd1888106dIlev31152j9wMGDJD0J3QZ2wCEFjCjCiAQEIhE9Wzbj5FAjh8\/Tm+\/\/Tbddddd9Ouvv9LLL79MkydPpl69ekn6E7qMbQBCC5hRBRAICESierbtx0gg7HUcPnyY3nrrLWfaUldXR4MHD5b0JVIZ2wBEEjKDSiAQEIhE7Wzbj5FAJEInWcY2AEnKnmZbIBAQiES\/bNtPgUA41jFt2jQ6deoUcUrDtra2QpKfsBnJJB31K2MbgDiy2qwLAgGBSPTNtv04BOJOLaiEUMl+QCCSoUu3DAgEBCLRsEwIxIsg1OYyDphOmDCBWlpaqLm5mXAaVzKMyZcBgYBAJFqVCYEoD6S+vr6QHZ2FVSQyfvx4J38pCEQyhOmUAYGAQCSalQmB6GTBHgcnPNYfTijEgq1duxYeiGQUUygDAgGBSNQqMwKRCGejjG0AbPQpiXeAQEAgEj2ybT9YxpWMSg7KgEBAIBI1rBgC8btk2zYAkkHJQxkQCAhEooe27cfxQGwf51d7ThgQfb8J\/9s2AJJByUMZEAgIRKKHtu2nMIXxumZSInDYMkxWa9ascU70PvbYY0W30oFA\/NEMSyC1f3mfjn91ga7rXU3t940IO0xlUT4sJmXRqZhCZkYgLLfftCJmnzyre11rCQJJjkD+8MRe+vTcBRpcU00HHrkpjSHMvE0QSOkQZEogNjXCRCCNjY00ceJE3I3726CEvYn+T6sPOR4IP+81\/dHm0Fp7V1hMrAmW0Yv4btyOjg6aM2dONnfjSvr97bffUs+ePamqqkpS3LeMiUC44tSpU6mhoSHWe7pKZdxEXzqSwKQYk\/Xr19O6deucP9q6RUG0jMt5QD766CNHOCaPhx56iK688spYtmkikBUrVtCoUaPggfyGMseoOLFT\/\/79qbq6Ohb2XaUyMCkeSfZAtm7dSq2trfkgkPPnz9O2bduchEKcQOiOO+5wMpMl8ZgIxBaDJtEXG21gvl+KMjDJeQzkySefpD59+tDdd99NBw8edKTldIZpPraDQGn2Jcm2YSwgEIk+2bafwCmMmrps2LCBPvnkEyeVYW1tbez4RxAQtgGQDEoeyoQhkEpYwuUxCYNJHsbQhgy27UcUA+GO83zznXfeoV27djlR3t69e6eCh20AUulECo2GMZZKWMIFgXgrmW37ERNICjbh2aRtAGz1K+57whCI8kD4nV11DwgIBAQCAgnBKmEIJESzZV0UmOQ8iJqFdsED8UYdxoIgqsQebdtPJofpEESVqEJxGRAICESiNZkQiAqSNjU10dixY0sykkkET6qMbQCSkjvtdkAgIBCJjtm2n6Igqs3DdH5g2AZAMih5KAMCAYFI9NC2\/WAVRjIqOSiTNIF0hb0iSWOSg2GOLULmBOLeYr5582batGmTlYTKjJ5tAGKPmKUGkjaWrrBXJGlMLA1lqq+xbT9FHojf9Q5MInv27KGlS5dSjx49uhQAqXYmwcaTNpausFckaUwSHK7MmsqUQPxuoMPNdJnpQ+HFMBbEQCRamCmBsIDsbTz99NOFXKUqXyofops3b56kD7HK2AYglrAWK4NAQCASdbNtP55BVP2ibRZa3ZEr6UDcMrYBiCuvrfogEBCIRNds2w9WYSSjopXJavUCBAICkahq5gTCU5j58+eXyDp8+HArKzG2AZAMil4mq9ULEAgIRKKrtu3HdyPZli1baNy4cU4CIb4bd+jQoVZ2qNoGQDIoepmsVi9AICAQia7atp8SAmlpaaHm5mbasWMHHT161AmcYhXm96ELM4UJU9akHCAQEIhJR\/j3TAlEv1xq5MiRNHfuXOfip\/3791vbTGYbAMmgRJ3CJDndAYGAQCS6att+SoKourfBXoiKh9hKcmwbAMmgRJ3CJDndAYGAQCS6att+sAojGZUYZZKaxoBAQCASNSxrAtEv6ebky15b39U0qb293cHDvbpjGwDJoMQpk9Q0BgQCApHooW37KXgguvGzUa9cuZIefvhh5zqHgQMHFnamBnVCrdbU1dUR5xapr68vuQbClDLANgCSQYlTJqlpDAgEBCLRQ9v24xCIHjydNGmSs2zLHkJbWxsNGzbM2d5uOkynCIhXbXjpl+uoVRy947zLdfHixbRq1SqqqakpwcQ2AJJBCVtGn7ZwXb6j9rre1dR+34iwTRXKg0BAIBLlsW0\/hZSGavmWjZqNnK+xXLhwoXP6VrKM6\/Ys\/EjHvVHNvU3eNgCSQQlbRp+2cN1Pz12gwTXVsTKk+xFIUjGWsH3MQ3mQauko2LYf6wSid9nttfBvCoDGxkaaOHFi6nfj\/nnth46HMOBvq+g\/Z\/5DJLtwt\/Gn1YecNvlhz0P9\/\/ea\/hipfa7ExsJ3n\/br168opcKYle87BDWwVxXFaT+yYBlW9MMkQ5EyfTXrR0dHh3Nvk61V00QJZPr06c7Gs6ApjI6we+qkEwj\/\/6lTp1JDQ0Mqg\/JvWz6j0+c76dS3nU77bIDtUwcVvUuVYXJ5\/q7+5P63Kly77kRRO6p8koL73USvZOJ3ueVP8v15bMsPkzzKakOm9evXOzMHfqwTCBu\/uv\/Wq7OSszCSICpPYfjhWIvXBdvKA1mxYgWNGjUqNQ9Efbl5aqEeNnwVr2BvRC+z7+ERJf9W9VQ59W9uk8sn+eAm+lI0gUkxJuyBbN26lVpbW+0SSFKKrq\/kzJw5s5A\/RCcN9zKu7RiIviqi+q1ub9NjFzz1+J\/\/\/dqJXfDj\/rd+45s75sHlk74RDvP9Ui0FJqWYZBIDSYpAkmgnbQCC9mW4yUUPfkrrJU0cClMYCwhEYl9p249bhorbiSrdl+EuJ60nGeQoZUAgIBCJ3oBA9u2jKVOmpD6HS3r5M257pvogEBAICESAgC0GTWqLuepS3PZM9UEgIBCB+WR7nF8iYNplbBFI0lMSv\/ZMnoXC0yQPCAQEIrE9W\/ajZKmIGIiXEQcZtm7Magt62Db8PBMpobiVhev935nv6O\/6XhlrS7xECculDEgVqzAlCKTBoO7lWd7rwSss\/HhtMVfl9d+9phimaQfXd3sWkjpeBhy1XrmQQRQ5K5FATB+gNOwnaGwqygNRpKGIQQHjXnp1L+fy715TDNO0wwv4KHW4nTta\/+pscuOds0x6cQ\/nRTHYvNWpRAIxfUhAICmuwuieBRvDP\/391YHTARPb2zQoZSxq23zcw3k2ZU\/rXZVIIKYPEAgkAQLxM3w3gZhIRC9vIpu0jES1y8bC\/eo4eflkLz9pbVpLuy9JtV+JBGLCDgSSAIH4uXn634NiIO4AqPp3llMHNpZ\/XPrXwhSm0smDxwQEgiBqKkFU05Iqv1Q\/Zu9njGEJx\/R1iPO78kDO\/Hi5FRAICMRLn+CBCDyQJGMTYZZzVV6PLIwXX9tScwEm8EBCeyBs8HxKlh9TIFEnBy7vlVrQFNWO4zUkWRfGAgKR6BM8EIMH4g6E6sft9ZyjOtHowLtJxxTVlgxa3DISjwoEAgKR6BkIxEAgytj0PR3KG1FxDfW\/7jJqANSKisRwJYMWt4zECwKBgEAkegYCEcRAGEivJD78d5XDQ\/1\/RS7uTWQcx5AYrmTQ4paReEEgEBCIRM9AICEJRBGEmzB0MtGXbNUg+O0uNQ1SXK8lan0QCAjEpJv8OwgkgED0L7XyNLymKSouYprmcMzEvW3dtEU8rtcStT4IBAQCAhEgoBj06smt9GvPPkU1\/GIaqpCeIFn3SNy\/e6Uq1MsELdNKphtB3YxaHwQCAhGYDzwQRSDf3rqshEBMAKoVlqAVGPcGMrcHwu\/IYp+HqW8gEBCISUcwhfltDlf37AHq7HODs8\/D5HW4QTXVCXumJWrMQjLYYcqAQEAgEn2piBiIfr2l+wIcBuBfXjoe2vsIAlcnFdPmM3c7UWMWXvLEISMQCAgEBELkXCalLtc+cuQIbdy4kZYuXVq4rtFv+iEBj8u4A6j6RrMoW9Gjxiy85I1DRiAQEIjEBrq8B6Jfus2XTM2ePZsWLFhAw4YNc\/DxOnIvAU4PoOZpj4cuexwyAoGAQCR2UBEEcvToUefWOq\/LtaMQiFd+jDhfe8lA2S7DmPG9p3xf8NChQ22\/PpfvAyalwwICeWJv6MApwzi0+jt6tq5vAdF\/\/68zhXbe+NfrcmkAYYQ6efKkc+t6Y2MjjR49OkzVLlsWmJQOrcLE6uXaNjUs7BTmb3740hHPvSfELTOX6\/Xf82x2Be8CArlEgD8wfDn9oEGDUpfPelJlSRBVBTtV7\/\/jn7vRoo7uznH8oIfL4QEClY4AE4cN8mCcrRMIv1Qt4w4cOJDa2toKAdRKH3j0HwiUGwKZEEi5gQR5gQAQ8EYABALNAAJAIDICuSOQoF2qkXtZxhXVUvfBgwedXixZsoQmTZpUxj1KXvRly5Y5jfLWgEp\/1DIu4zBz5szUMckVgZgCrJWoHGwcvO+DSYPxmTZtGq1cuZLGjBlTiXCU9FkZjA1jyTvguv306NGDFi1a5OwbUps005A\/VwRiWuJNA4ByapN37jY1NVF9fT0IhMjZiNjS0uIYyPfff5\/61zbvusL2w49NDzV3BBK0SzXvA5i2fPoXpqamJu3X5b599s7GjRtHx44dI6U3uRc6RQGZQA4dOkS7du2iU6dOVd4UhgEAgXhrGH9t3eeGUtTF3DfNU5fdu3c7XoeuN7kXPEUBmVCZOPhwKnur06dPd\/BJc7qbOw9kz549BQBgMJe1DZ5HqdWxsaxevbroh0qPg+hEqqa7Y8eOTXVKkysCQRC11FAYEz5Et3DhwkLKgxQ\/YmXZNDyQ0g8N\/6XiPBDuNHap\/m7D6ivS3t5eZNi2DkqVC5uAQH4fKX0bhI0l\/1x5IOWisJATCACBywiAQKAJQAAIREYABBIZOlQEAkAABAIdAAJAIDICIJDI0KEiEAACIJAy1QE92q66IN0Hoc7U8Kaj4cOH09q1a8m0s1Xfd6GvAumHt1iOqCtEQTLp77axslCmKpGJ2CCQTGCP\/1L93BAfnPJKUO31Fl4a1g9ZudvxqqOXOXHiBM2dO5eWL19O11xzjXMWpbm52SEgJhM2dgkh6e8JkolPIasdp9I+xkcXLUgRAIFIkcpZOS\/D10\/uuveQ1NbWFt2\/o7oj2eWqH5cPMmKJgUu8iSCZ9D7mbEgqUhwQSJkOu5cHom\/9lxq9xAPRPYuzZ88WPBD3MXE2fOWdeB0hd3syfqkJ\/GSSEFSZDmfZig0CKdOh84qB+MUHvNIA6ImKJHELFevwi5mYzl54\/e72JoJkUv3186TKdBjLXmwQSJkOod8UhrvDJzD1oKTqohdRuL\/qX3zxBb344ovOUfDevXs7VU3ejCIHTpLtlxXMb1u+V+A3yNOQeExlOqRlKTYIpCyH7fKZIXVymYOo\/LCXwHcNP\/roo\/T444+TOokZlIjI7Rm8+eab9MEHHzhHwDnXhpcx64Ti53lcunSJzp8\/T7\/88gtdddVVdPHiRXEypCBvJmqgtkyHOfdig0ByP0TeAgZ5ILNmzXKMVRGIft3h9ddfX5RXRI9b8F0iL730Et166620c+dOJ31iVVWVpwcyefJkqqurK3qPLinLxyTC9X\/44QdqaGgoIj2Vr4LbmTBhgq9M+\/fvL0oWhPyn+VJYEEi+xkMsjVcMRI8PuJPrcsPu3Kq8D4QfNbX58MMP6dVXX6VbbrmFtm3bRjNmzHDquKcfatrhNU1S7d14443E7fF\/33zzDT3wwAOF6ZDK46FPX9xt6dMtfeUGMRCxilgpCAKxAnP+X8LeAl\/yNWLECBoyZAgdPnyYjh8\/TlOmTAkt\/E8\/\/UTPP\/88jR8\/nnr16kXbt2+ne+65J3Q7qJB\/BEAg+R8jKxJ+9dVXTvyEk9B0796dvvvuO3r22WcdL0QFU6WCdHZ2OmTE5HH69Gkn\/vHggw\/SFVdcIW0C5coEARBImQwUxAQCeUQABJLHUYFMQKBMEACBlMlAQUwgkEcE\/h+hneO5cjY5EwAAAABJRU5ErkJggg==","height":164,"width":272}}
%---
%[output:8766cbe6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:272c8b13]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"fit_test","rows":2,"type":"double","value":[["0.0003"],["0.2787"]]}}
%---
%[output:8d098593]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('normplot', 'C:\\Program Files\\MATLAB\\R2025b\\toolbox\\stats\\stats\\normplot.m', 32)\" style=\"font-weight:bold\">normplot<\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025b\\toolbox\\stats\\stats\\normplot.m',32,0)\">line 32<\/a>)\nInvalid axes handle."}}
%---
%[output:83915ba8]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap\nTaille de la serie    : 352\nStatistique T_max     : 13.4869\np-valeur (bootstrap)  : 0.0260\nPoint de rupture      : 2017-12-01 00:00:00  (indice 293)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap\nTaille de la serie    : 292\nStatistique T_max     : 18.5135\np-valeur (bootstrap)  : 0.0020\nPoint de rupture      : 2008-10-01 00:00:00  (indice 183)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap\nTaille de la serie    : 59\nStatistique T_max     : 4.1240\np-valeur (bootstrap)  : 0.5900\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap\nTaille de la serie    : 182\nStatistique T_max     : 3.0544\np-valeur (bootstrap)  : 0.8930\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap\nTaille de la serie    : 109\nStatistique T_max     : 9.7080\np-valeur (bootstrap)  : 0.0990\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:1f3e4364]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:1a8192d8]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap_ae33\nTaille de la serie    : 385\nStatistique T_max     : 55.3086\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2022-11-01 00:00:00  (indice 352)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap_ae33\nTaille de la serie    : 351\nStatistique T_max     : 49.9353\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2015-03-01 00:00:00  (indice 260)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap_ae33\nTaille de la serie    : 33\nStatistique T_max     : 5.5399\np-valeur (bootstrap)  : 0.3290\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap_ae33\nTaille de la serie    : 259\nStatistique T_max     : 6.7403\np-valeur (bootstrap)  : 0.2980\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap_ae33\nTaille de la serie    : 91\nStatistique T_max     : 18.8748\np-valeur (bootstrap)  : 0.0010\nPoint de rupture      : 2017-12-01 00:00:00  (indice 33)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap_ae33\nTaille de la serie    : 32\nStatistique T_max     : 2.0951\np-valeur (bootstrap)  : 0.8960\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clap_ae33\nTaille de la serie    : 58\nStatistique T_max     : 7.5121\np-valeur (bootstrap)  : 0.1820\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 377\nStatistique T_max     : 13.0027\np-valeur (bootstrap)  : 0.0280\nPoint de rupture      : 2015-02-01 00:00:00  (indice 259)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 258\nStatistique T_max     : 22.9387\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2001-12-01 00:00:00  (indice 110)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 118\nStatistique T_max     : 33.2074\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2016-05-01 00:00:00  (indice 15)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 109\nStatistique T_max     : 12.0046\np-valeur (bootstrap)  : 0.0420\nPoint de rupture      : 1994-01-01 00:00:00  (indice 26)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 148\nStatistique T_max     : 5.0383\np-valeur (bootstrap)  : 0.5270\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 14\nStatistique T_max     : 6.7535\np-valeur (bootstrap)  : 0.1480\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 103\nStatistique T_max     : 11.3077\np-valeur (bootstrap)  : 0.0440\nPoint de rupture      : 2018-08-01 00:00:00  (indice 20)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 25\nStatistique T_max     : 3.1008\np-valeur (bootstrap)  : 0.6880\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 83\nStatistique T_max     : 4.6158\np-valeur (bootstrap)  : 0.5740\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 19\nStatistique T_max     : 1.1798\np-valeur (bootstrap)  : 0.9680\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaG_ae_psap_clapS_ae33\nTaille de la serie    : 83\nStatistique T_max     : 7.8861\np-valeur (bootstrap)  : 0.1760\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:470a0230]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAO8AAACQCAYAAADgBDytAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQlAVdXW\/kBQnAUnVFQsc9asLDUtsTRNm8ypzAJKG94rzcysNIEcKofUrKynKViZZnP5GjTFbND\/aWmjqSko5DxPmCg\/376s676Hc+49Fy56wbt7POGcffbZ03fW2msMysnJyYHNsuvISQxZsAEv3dEENSuWsflUoFpxmoHT+zOxc0p\/RA5fiFIRdZxdt7penMZW0voa5A14OfjXVmzHpt3HMKVvk5I2Fxf0eAjOzFFXI3vvNst5KNOoA2qN+i+Cy1W6oOfKXwbvFXhJea+bsgY7Dv+Tr\/+hpYLw9WNXoEWtiv4ytkA\/CjADAQpbgEk7T494Dd5bXlmH8bc2RJemVc9TlwOvDcxAYAY4A16DN3DmLfkb58CiZ3F83VeKRT7+8xLsntIHweUqo9aoL1GmUduSPwHFZIRegZdjemlZOmZ\/m4GMg\/+ga9MIvNCrER6c\/zsm9WmMlrUrFJNhB7ppNQM628w6f4\/tjhoPzVbVd88chNrPfOEiyArM5PmbAa\/Au3X\/MVwzcS2O\/3M6X4\/LlS6FlU9cgQYR5c\/faAJvLvQMELwC2FP7M3Do06mKAp\/K+CMA3kLPrm8b8Aq8uqpo95FTeGlpOqbd0RhHT54OqJB8uy7ntbWjq953YZUJYrLONYa\/hwrtep\/XvgVefnYGvAIvHxuy8Hd88ds+7D2a7WylfOlg9LqsOl7q38zv5nbVqlUYMGCA6tdzzz2H\/v37q98XLlyIp556Sv0+f\/58VK1aFfHx8XjkkUecdYyD0Z\/hvZtvvhnPP\/88ypYt63fjDnSo5M+AV+A99k827k35DcjJwZINB5yz07Z+RVQqF4o5sc1RvnSIX82aDl4dbC+88AJef\/112+Bl\/U8\/\/RRz587FJZdcghMnTuDJJ5\/Etm3b8MYbbyAiIsKvxl3QzgSMNAo6c+f+Oa\/Aa2Vh5c+WVwLea665BocPH1ZAY7nvvvtQqVIlrFy50iPl3bRpkylVtrp+7pex8G8MGGkUfg7PdQteg7fdc6txIOs0gqhnAkDbSv4E55xBeIUy+Oih1n4ldRbwknVOTU3F5MmT1Rw\/\/vjjCsBjx471CF5hl8let2vXzrlG+\/fvV23w2siRI8\/12hXJ+wJGGkUyrUXSqFfgZQ9oHvnRul34v\/Qjzg71bF4NNSqF4sYW1THr2wy\/Yp8FvM888wzWrVuH2rVrIzo6GgsWLMDgwYPx8MMP2wLvjBkznCyzDFxYZ\/4dOPsWyf4MNOpmBrwCrye2eUzPi\/Ds4i1+5bgg4KWwKi0tDRs2bFDT0aRJE3Tq1EkJszwJrC4kysu5OXP8MHaM74GTG79zbp2AXbP\/fUe8Aq8IrNZsO+Jkj3\/5+yhum7kObepVxIC2tTB\/9Q6\/pLwEb\/369V0kz\/K3J\/BeCGde2ZoC3FJVoxD52ALnjt354h04vS8j4JjgRxi2BV6huKSs9877HZt2H883hAbVyqJy2VKY3r+pX555Cd7LL79cCZ5YKDXet2+fLcrL+iKdlnOvnHd5T5c2pyFNtR+NaD9aZvtdCUib7c\/V+a5pC7znu5OFeb\/ONt9yyy0u6p2NGzfmA+\/ff\/\/t8roHHnjAKYyy0vPuKrsLKUhBat5\/0gABHIc4xCK2WAGaVDZr\/RdOW+aTG1djx\/huCLu0uws1Lsy6BJ4t\/Ax4BV5PZ94L0UmfgI1HPITieloSAXR91EcMYvyWQtM54cC7Cc7hhPdLQnjfMZ6GF7h\/DmfAFngJWroCmrHL0tfrGocjJb6Fi5GGnJGX\/LFfVYtvX0s58XOjc9PbLf66yTmGzuishsE+8r\/E3P\/slrmYqyizt0XMF\/mcFajM6hgFUSHV6qHO+O8DjgbeLoCf1LcFXukrHRO6TfsR1zeJwMwBzZ1DGL5oAzIOnswnqNKjbsgHIL59bYR1+lpRK7uloJvcbvsFqccPUAM0UI8SsAlIQDKSXcZFMPM6WWreY9mKraqOfLyWY7kCvd2iOw7wGTNPH6s6rC9OB+5c+wLSZrurcX7reQXewrLNBPklNco7wZujzDvclyAEwR\/BK0Al8AhAFh28AmgZHSk0AUtKy\/FwXFLszIPUJUXdn\/KYophBYRWVSqfyzcNcHAas6oRGRHn0DApImz3tSP+57xV42e3u09e4GGjIUCqXDcGqJ6+yDExH4MfO\/VX5\/a6t\/Z6iPqdOn8GpM2fczka50BDMOj0Hd52557zOWo\/QLs5zrVBTAWMndEJ6rpA5EUlAtEPanJg219Hf1BggPRqp9ZORGufgNhKTcpCYcBa8d52+B7POzLEcX2hwMEJKOeoTmOKmx78J3nKtb3A5j1rVYX1P51hfSpszMjLAn+JWoqKiwB9\/L16BVwB4VYPKeHlZOs4gCMFBQaBjwq2ta+LBTnVNxytn38Edo1T4HKFQCZ9tRpcmEahV+WwkyuzsbBw\/fhzlypVDSEgIGlYvjxeOvI7eWQPztW2s626yC1u3U0QzZJZK93494+cCyXnn2pw8wAblADlBqPPXtci8+BvV5uY9x1zalv6Wz9qLCuveRY0b7kVI9ehCgVcXOIktc0Tsi\/nc\/HwhbSZoR4wYgdWrV3s\/Z+f5ibZt22LSpEl+D2CvwSthcD78cbcziqQ7xwSd4kqkDQFv+GPL8MGDl6Jtg8rO5co6kYWdu3ahVmQkyoSVgTvKa6zrbs0LW1cory5VpuT4x8O\/4pnR\/+CVTw8BWx1nYF7\/\/dTmfH1vFtpQUe9fTvyJlmUbY9RfczD+4nvVM8dPnXWx5N\/S3yr7NuDElJ6olbgcZZvHKPAWlG3WfXGFPTZSbZnDwkqbRUVHENSpczaE7HnGpcfX82Mzffp0ZXWn27F7fPA8VPAKvKSgt7zyE36kXXNwEILO5CAnOAg4k4PSocH5okcKqJ\/v3RBPvr8ZRspL8H7yr8vQsWGVs+DNysKOHTtQq1YthIWFqbOh1Zk3y1DXLXh9VNcolDqRdQKXjHsLGW93cYJXzsHGvss5l8+UDSurWOvk6CQFaBFk8XeC\/\/Osz5Uga1fYLtT8LU2dqwnewgisDn\/t8KgiBVa62yl9UWv4oiKJSyXgNQOBaBtWYIUaK\/Xg\/mLU4q7f5wGfbl\/pFXgJxhtn\/IjwsiGYP6gVSH2f\/nizesG7g1vliyhJAdXcH3a4dGDCrQ2dAisdvELRIrMizyl409IA+UlPB+rXB6KisnHxxdudHxB9AAJebjb2edTSJRg\/PsRxthW2OC0ay9O2onNMEOJS56JTWhwQk4r46M6Iyo7CpuxNCrzRaTFIi3aozCi0ovRawPtH1h+4JOQSZIRkIGoPsHG3A7wsuhpIolsYz6pmdYxS5KLU3ZqBgGOj5N1MnSaCPJlrMUnVjWZ0gxnjrqbFGx1PrrvuOqxfv14du2hBl56eroIrLFu2DK1bt1aPST0zZJRo8HoTPZKUeszHm3FvhzrKrFJCxgoAKo+bjyvvWI3shutd9L7c4PwajwsZ5zXltTJPNFJpAjYlJVfIZKGWJYBjY4Fx484GFyDMHEqfeCAtMZc\/TnSANmlu7t8Alnd2CKzilwPRuWCmUIqCK4I3mvdSEbX0LoR0+c4h\/EqLzhNw8VOwXDP0WIE\/skYiDjHIPOI4M+rgLWoK4AtVkRkIknIFegQuR0uuQgn6kO4Esw5ggvfLL79UXl8s9OAaP348YmNjVTAEY+H7BKgvv\/wyunXr5qzHZ1966SXlRUarOql3QYGXgy1IxgTR8RrBG7w\/Emcidqo5rHO6vvpXFwoJdbOSNp\/MOokdO3ciu\/ZJDCn\/L7X5jWdSbognTo2G1I2sWRNr15TD4PuCFMX1VKKjgbgEIDWXeDpoJPW1lBpvzYWy44ybv5BCUn1EARWlzhR08SuR21jacsfjy\/Vnee+sNZO0F5WdjZCQrqj5W6qTbfbU38Le95WqyAhe0Ys7PlPLXdhkXWfO4wPruAOv2KRzrKTGBCX9qklx6fpJH+1LL70UDz30EA4ccER8YcijG264QXF1v\/zyi+WZtsRSXk7C3W\/8gsW\/7XE44efw6BsEKnvqR4Thq6GXm6qKrMDL9tqeugZDjo1S\/7KcPHkSvxxah0kXjcbast+ra1bS5vScLVhyYjGSa83AztKZlvuWH4Zbj92Ba7d3Q\/jhy9Dt6rNnbI+bXXDIigR7tICXjgc20O\/xBWyHQLd2ZIjak4qNu+Fkmz02WYgKvlIVGUEgVNfIHktXqTokTyP6cTO2mc4lXbt2RVJSEhISElToIdqbs9BDzEh5hW2mTfsFT3l\/3XEEt7z8E+bEtUDTyPJOvS0nb8R7fyrzSLMEZDp4GzY94rRM8maPGaWxfPbt4HkYXMohrbVT\/v3ji3j\/5oex8+9QVd0shG3X4OuQhvQ8DoCAIoXNJaJJuf9HArk1GYi2bx0m\/SKjSPUPz7BO0KdFI5rsNTqpd+Yv5EYcaqYle9LQpfq58VSipPnIsjecppOiVqp43X227ZuN4BVwWgkf5Sgl4DZSXpkbq7MwfbMD4HWDAppH3jj9J7w2sBmys3OcUTO27M2yDd6Mpp8qA42QzZeizP+6Y2S3aNSNCHO+lWfTo0eOokp4FaXnFW8d4xdbZ7X4sG7pZByCbBx1nTrWvBKT09kpIBIrKREaOaqQIsac5ZR5ns1TBzkbSUoEOsUCMbkHaDd2zXIEONs3stAJiE7hP+4+Bo4+UDxGhrKoip0YVt7YQvuC8upnXhk3BVM65ZXrZmfeAOU17JZJX27Fc1+mIbx8qHLI33XoJPrN+tlU2sxHl\/6xD31mrUcwghDbrhbW9XtICae4mbftz0I9DbisT2HVvO3znJJeXTWjmxGKuaHePbYpZyreF8ktgek0R6TBREou2DTpMJ\/hWYvlLHjz+GVyxnI8jUsG5hqA1jnHgXFVSJ475bHArGfmfEEAkprmnXHp17A8SPXb3DPp7Hk471PiIm0uSolxYT4URvDqThxGe279Qyz3rCgv+6RHBOXfeqAFiZjCSJ\/Dhg1T3I7uCspz8IMPPugSBlgfZ4k987rzLrqkRjl88u\/WLmyzbqDBCXr4q+X4Ju4Wp7DCyo2ubVZbBXDqeRWxzLMD1nWB7lzwjEDIBwxSy8QEp2pHQO8q8MoDTWoa0DmP4lGaHGPfGyr\/5jcRTPEITUm0ZbtnnyHkZ2sZDdi+v6YgMQOBzgHxbEu3SCtpc2E+HIV5tsSB104kDTPwkuqO+nizAnWFMqXQ\/4MlWHxHD4\/g5eTbBWqBFsoAXvM28uhcUi5YE3NZZyfLTCDzJw\/EVAU5uVmyzix0vl+R97vDm8hRhHZql9TtXNUTqbppcQXvDBsWVgWaEx8\/ZAUClyOM9k6jI4ePu2O7uRIHXn3kRh9d3mPCMbOA6wSvRJNkPQZsf3ewQ1FuzSbanueCVyw0eEkD8xTEqTnqWKwKuWZikBw42W1iOYE8N3W\/c4FYGmsYuq3qJeVy0VZ+wGfBy2\/ELzYcEwo+Mb570pOFFWUZ5HSo6xW9r+\/eXvCWSix4BbhRVcoop3opVv68ZuD9+u6bsC8sw7Y5nJXRhVGf62m5XNjsPPBG51jpaQk36mcprEoF4vMQl2dBFZ3rRJ8m\/sgaeGkspXTHtMlQllupQIzDWT86ORFpMbG5vxiETqmsm2RJedO08zF78ek5AK9upMGokTUemoWdU\/qpbIF2U3wWJxDoe6c49dtr88h\/zf8dp88A32w66Ezxed+bv6J0aCkkxzZ3OfMa2WZS3jK3zceb1Se6lQ7rkymCKSNbpUuF7fjDOgVWuVZNiSlblWUV9dRWheImxciKWpe\/U9Kci8pELEd9Bd94IDlHtDlOWTNppYOOpiIxL9IGVSQrEJfnkq+9NQlY3in3yGvhjy8aKj7BNocXMdusG2mE3zQM+xcloeawBTi0eJozZ29wuUqevpVOoZJu2+wpgoo\/REwpseAl5W3\/\/GrsOXoKr9zVDJ+v34PBnaJw26vrUK18CH54sq1LGByjwEp0wZEVHYIoK4W97AzdnM5olSP3WNcr8CYlIiY1Aamp7sHLE62DZgLRDfIoaqKDveUmo\/kmwRsdn6MsIB3jcZyElaWk4qIdDvgsBC\/D0UmbMkaqjLfmPW+GCB28bPOaIhZY6UYa2XsznODNyTqCnVP6I3L4Qlthc8xAYHTqMI7XqAPW80kxWL7kieJ1FslSIeqjJ554Am+99ZayuJLcUZRap6SkoHHjxsqwQy96m3K9xIJXBFdX1K+E575Iy\/MqAp7q3gBr0w+bBlsn9aUqiUWcF4yxnxwb+yw7yS+0hIrhdTq\/m8V6Empq9yOgbIk7L3fYFOdZiLkjIQSagl4eksnxpm11+OTKmT0uPgfJlE9RppUHXBLRmDw7XrFflo3pbDOv3bgUhqG17oVwALqhl5nTgUdS6EUFidEc0TcBBxdPRbXYKdgxvjtKX9LedvRId+A1+9gaPbD4\/IoVK5wAFRCOGjVKhdr9+uuv1T267dkBL59jNkfWnTVrFoYMGWKa3bHEgpfrL+fbYV3q442VmRjV4yL0n\/0zOl5c2eUc7Gmv6AAWMPBf\/WwakhGC8Onh6LC5gzM2sq7j676gO2ZeNdPTq5z3+\/zaB2tuXuS0aXbHNvOhhLnL8Wx8Hq0UUryVh1Qj\/fTchQrv34SQjC7IufZaHLrsMscDeYJsq6dJwZUmKc9C2n6kK8\/98VRD\/0Cwrrf65MKCVwerMYUqHQ\/atGmDd955x0lNabjhjvIGwJu34jo15SW6+VlF0XC3SaxcxOqcqoPqi6vj62u\/VuwP2SS6hnFxRo8ejaefflo1+\/R\/nsZFz16EiWUnut2LbO\/2w7djx\/M7cM+1b+CWWzyf2ZwN5rNt1kix27eKMYbD4cIhZs7jLpSNtCf4OO5LtTwjTXsP+UGtwoKXQ9BNIXUWV7yGaEFF6kw2WcA7ceJElYpVL3pq1wua8hqBK5Nkpuf1Zg9ZSZXZBjcCv7K33347pk2bpqgwv8bMj3vnnXdiR5kdePrPp7HhDkceIim1TtZC5Q8r47sbvnOpn5XVDp1tEk\/xKvrk9sP4sZIOetEFGVU8Bgsq4ySkOUyjySonRZvbYMkjnl0WvJnhc1vXF+DVe6ybRVIIJi5\/\/LC3bNlSuQ8GKK+bNeaZ97opa7Dj8D\/5ahUWvO62FheImf3oOUIQMyMfC8F79dVXu1xPD0pXlLp7k+6W9emczZKUKw2i1DmybCamX9UfQ\/9vIXaecIRsUaCNA0TGwTbLN2+Oda1b44usLGRmNHR0WWGYeiH+Ho169c5gW3CwQ3pFS0kpaUDI2xmIQRqWLOmoX1bqYN1ui6DNO0afW8QBKlJH5qirkb13W753F8a2mY2JwMrOmVe8hWStdJ9cHbx62hlSZCuB1QXPNtOrqPv0nzC1XyP0vTzynGwsLiKz+1E4IRTYHXitQK1flw0hA\/h2cSZKLeiPtbVfRIf9j+HYbQvRsefZuEtWffjmm2AkJCxD3bp10bNnCwX4Vq0cOXu7dh2PPXuaqVcwOgfvtW17wvnBMfbhnExmIV5CIVaFq\/vmC1Zn1aSvpc18j6idjM72fBc\/rgHweqC8A2b\/jEY1y7kEXS\/EnnD7qFBc2eiySEa2mY2wrt3rZoHFuDmP\/+DwDTUrQnX+t3G76bukTaEEIgk1tmVUcxTV3Pm6XW+TbpuBN6Dn9e2q2DLS0B0Scs7kIIhudUHBLj3xNdvMTU4fTWMm+uHDhzsFVhMmTMCUKVNUP7y5LjpAGYBszOr3TceeN4Y6dZne9EHaNIJXV3l4ArZvl9a3relRK0tFeI4GWZxULvpMFad+2wKvDM5br6KCbh+jyxfbEYkhQ50wIbbORvF3\/Rndqsfqut43AW+1+Bexd+5jCryksPIeqeupD6xnBlDd2MBdELWCzpcvn3N35pVgd3beJ\/M+dOhQMA5ycSmZmZkq3nSJC\/1aXBagIP20yzbboToFeX9JeyYQdL3oV9Qrylv03Tm\/b\/D2XHd+e+v\/bw+kOynaNfIavLpLIF0BX+jVCA\/O\/13lIJKMCEXbZd+3rrPNu6b0LbSaxPc9LPoWfR0Gp+h7HHiDV+AlcJksbOOuY9h+8B\/Uqlga7z3YCnfO+gV1wsOw6IFWLo4JgekNzEBgBopuBizBS6U4jSB0U7PO3W7EyvBeuL5FLXRtXh1JH23GkuGXI\/n7HZj9bQa+fPRyNIgoX3S9DbR8TmbAmKdIXuqNkcY56egF\/hJT8ApwaU8qble89ujwJ7Dsz31oeusQ9G9fHxP\/uxXv3N8Kd8\/9FaeyT2PJsDaWKT6Lgwhe5e8Z3w1njh9yboui2LAyv3wJDU6Mhvf6ntQti4wqLl\/tXb0\/40cOwYFJtynH+wOfTXUaZlCgV6ZuM9uhX33Vt0A71jNgCl6r8Jq8fut9wxHcZiB+2XdWz1utfChublXNrVeRv4NXvGiM3jNChbxRk3jacN6Al2omzh0NUM4VeA+9EqvUZUxMdnL778oNMCDM87Sq5\/5+UIMGDXJE96jrI9u3b4\/w8HA8++yzymODZfPmzfjtt99QrWZtbI\/uhbBfF6LUif24qsedWPDy+Hy9p1kh00yw8B2vv\/66U38m+lDqbVl0vZreD6NHiKS1cPcMU13IZje+h6FBxWLLeC9xcB\/c89REpzfLHXfcofxGlW65aQU8u+hbmEWREDBu27bNEmT6u2h8wlKpUiUn5dXnSsb8ySefOOePXNDMmTMxe\/ZsHD58WD1Prxs6qP\/444\/OevIs7\/PYw1KxYkU1v1bzovdnwphROPji7Sphd7lLuzkzCR5f\/6VLIPZzv1UDbzTOgJPy6l945ndhrFtuSqtyplRpHI0Zg+B9m1BhXQqOdBiBStEtVSxnSp3FnYubiay3AFI2kWwsso3cpDNmzFAbUfLQsF7VqlURHx+PRx55xBl7V1hN4zPisED1BJ\/R3yuUa8mSJc73sG1+CGjB9fjggfj8iVsxZNl+tcnlvVdccYUC10fvpOClFydjbnIymrQ561Qgc2MHvBw\/5QccIwv7KO2L4Qnf3ahRI2e\/ZN6k\/+JJtXbtWmdUCZ2j0Z+lsznnWOrKOzmX\/Hi560\/pE\/udUTNIfQ+864hA4UvuIwDFws+AAq9QBW5kOeOSEhBQ8rUX0Mz8ejPenjQC\/Xp0djoL0BKp4\/3P455br3dGi\/zsw\/edQGFWN32T6aDkRhKgT548WY3IzIJK6sjm058RE0ozSya5plMdvkPvz5WN6iLthX6YfrgZgkqXw6BBg1SSKvkAfP\/VZxj44BC89dpLuPqGm5yzbqTccsMYXsU4v0a2mXl0dNZY\/Jf54dDvCXhlLfi3ka2WZ8eMGaO4JuEG+Ix8rCQxl6y3N2x84bdcoAVfzYACrxEYbNwMvI+NTsKQt37G5kXjcKDiJbimz\/24s95BPPZQPOr2ScSrj\/bCs4u3qHA4qf\/9yBK8OkD1gZCllej2upSbFMndM7xH9lzPGCcb05jbRkCsU+GLqpYzBa98KKzAa5fyugOvgMzoQC79ZMgWI+UV8PJfo0aA1\/gsPWzomC51CdAAeH0FG\/9ox5Ly0tD\/ww8\/zN\/L8Po4fSYHXa7rhOBL++HqyjsxdthgBHd5EhMevBnvr92tYjiT8hJQcpZ1R3mtpkLOgdyMjKBBO1kBlBE4ZptUuAje01VfBDnPeaTw7B8p719PXo0X1h5XzfatcRCJWyIxsNYB9Kh+DOsOl8Hjm2rno7zegrdevXouvsjSZyPl1edDp6x2KK+deRHKa9YfSpv3j70OpSLqotao\/5qe8f1j6wZ6YXnm7dWrF2688UYlVPn9999RrVo1VKhSFSv\/DkJQ+mpsCqqPrOZ9UGrPBlT8bhIa9E8EajTFtP5NC3XmpfBFzr9RUVFOoYtQKNnw+pn3gw8+cFInUlShwnLuk2fkPMwPAFNFCiUSX2EBs5Gt94Wk3JdnXhkPwezpzGv1UXPXnzI5p7BjfA+c3PidEyHl2ve3HXwuAKtzMwMuqiI5H+bk5KBhw4aKdV68eDHmzZuHZs2aqTNu3weewKE9O3H9ddfieJPe6Flzr2Kbeead92Q\/l17rHj0ibQ4PX6diSB0\/fgKffPIxdu7cpZ7p3bs36taNUr9\/++23WLNmrfo9J+da\/O9\/FVSWcysJtTHWUfXq1SFURdhFkWrbkV4bjxG+AK9O+e1Km3Vwci4o7OKHSgcvfzeTVPO6LhTU2Wauo6f+GLefty6B52b7XthvcWukwRAzl19+uYoNxPMTqeKCBQtwz8hJuO\/d7c6Z8yYAHbVOhvC5blfA2\/oX9nL6ZvR6xgRpMbhcZdQa9aXtjAm+6UmgFXczYGkeqRtqCCvKhuakvImZGypgyR\/7ne1a5Soye7G3YPS2fmC5CzcD4qDAVuqM\/95WgPXCvTHwdEFnwCvHBCbX7jbtR3RpEoFXBzR3vtMqV1EAvAVdlsBzgRnwPAMu4LXSW+rNBFe7CHePnIyE3o5sfyySSYEqopoVy\/iUDT7XlNcYbFwGUxQ2zp6X5\/zXoHkoS3jfMee\/M4EeuMyAV5SXIG333GocyDqt0l3zh7m6+BMeVgqrnmrrFXg7p8UjLaUT6tTNQWbwNiDdEaA8LhaonxdtPCUZiIlLxzPHhykJaOWbhzkjGOpA022Sda8Y3SrI6rrMiLCMEbEv2o6SaLWfiqsjunE8R1fMU5cqdLrngoQONR788cfiFXhPnTqFoKAgvLJkE5K+zASCS6kxeRJYiQP\/4I5R+O7dqoiNBe4bmoUt\/x6Jp088jlaNclAhvi3K79qJ7N3RmBEdh0+OJ6jA5CtWAHUfn4ouj09TTvICRgLt77HdlfcLi2SIZ3Is+f3ExlU49OlUpa88lfGH6XXdVtlXxvfFOQSMP27S89knxt+aNGmSXwLYNnj37NmDOXPmoFGn3hg0exVKZ67ByUY3AqVKq7m1ih6pB61jorG3J1ZV9bsMPIhB2Y8icddIrFgUiRO\/rkca03e1UdzTAAAgAElEQVSlpaFVsxW48V+d8MJ7Mcg5nY1eozrgwYtngFEuhCrqqougsIpOqnxq+2\/OVJTMbCcAp2H98XVfKSDr1435ZkmdxZOmoJtGVEtc9Dp1PEdaLOh7As8V7QysXr0a06dP99tgdLbBS93ntt0HMGTJKZUxISRzDc6EN8CZcg4wmoGXFHfMx5txb4c6uHfe7xjSphGWvh2uKOo3P2ahc1oczvzSHVF3bUZGw4ZAaowjgx9Z5hVMSO1IrRm3PBmzW3VVkfx18ApV5fvJUtMThkV3Y5NnCGqz6xXa9XbugKLMFlC02yzQelHMgC\/0+0XRL2nTNnjpaTRmwiR8uL0S7m4biWpn9uHhhx9GmTLuBVR8kVDfKqtb4+rWZdApIQ2dU5McubdSmdgjN2lPDBN+OFwPHblCYh1AZi6Q+snIvqXoweurifbFohttstk3O+FIdYMNPqO7QPpqfAVth\/OSnp7udMksaDvn6jlfrGNR9tU2eNmJzAPHMGDmajSsXhaz7muLYObksVEI3pjhW3Bl6fpY+XJzpAl1JVllHhBHoh+tJT2FXl5iy85zkd6k\/jlhm8UFTh+aN9JmXyw6wcvkWfxAsrhLeSn91NOy8JpYUTEZm1mWCBtL59MqAfD6dDrhVQyrk3WuwvHWcQgKCfUqYwLBe0WvPRhw30FMuvNmLYe8pNdyTUvtGKKefS8OSJ2LhT\/1R4\/2\/ZQkuKgEVnJGLkwImKICLwHco0cPl9SXQpHpy8ugCcz+bifihm4eaZYhXraZmZulTt3F+2njxo347LPPkJqaqtKxmlF8HbxmARdor868VAzaYJfT8C0cXFvzxToWZf+8imFFW9lNu49jVY0+TkGVdM5KYMWUoH1mrcehIfepqtHRMUhL6wREk0VOA1LjEJ0S7UhWXz8t934q4jrHITUuDakJTGQtFHk52r4CfPzPa6g5bIFqyypDvKiEjCZ9VtdlDLq0uTAhYHyx6GZss0Q8oZ02wSq2z8xRy\/SnKSkpkGx4OsDMgKRn4XMX9siYRZ6gJ1AJWhZJCcPfJSAC65h9SAS8dArhMeziiy9W3MH48eMRGxurzG8lqZwdTqMogcG2fbGORdlHr2JY\/bltJ\/rePwIzX5yADs3qeuwXKS5DxcbdUgp3pg0DVsQ4clemxAIJDRQFPlHrC4Tt3Olsq\/NyICEJSIuJxoqEGCw\/sxTpwRl5yamXY2Hfe9H79bFFYrYnNr2FDQHjadGZGljyay+3mEUj26xvcvG2Eh9g8WN2B5iOHTsqpwaq+hgAgZoDoXDsAqnvq6++Cqb7oBcZwyDx591333UJZUSnBuOHRfyt9fOsZPLbvn27AkHPnj0VUKWO7rQilJ\/gZRpXsvjnIvCepw3saR09PV\/U9y3ZZnZcvqQS3fDppHH4YGMQ3ps8BC1qVXT2zcrCilR31Meb0WD4K1iYuRJImgt0ohg5Ndfygts3DjlBrsmpk+OAFZ2ATunRyoNhzuHJWFnpWN67lqPvj79jzj8nCm1EYTWxRupbkBAwnhadny0Br1XGezPw0u+Xvrh6flq+i5SXoDKeeYUyMrcx9ZVkp\/\/880+Ehoaqf1n0VKOkwLRj79u3LxYtWoTrr79euYKyyMeDrqJLly4F3S0lkod4SRkpL10uCV4mwuZ4SG35N51d5DwfoLwFh7gl5dUDvZk1n12lAY61H4qcMg4QmzknELzT1v8fLu37AyaUGgvELwc6pQBxPOtSopyKxCQ9rTSQ2ilPW5TLUcfEJOaefAlu8tSScjoO+xc969fmep7AK5SX8+YOvIxzxfOjFJE2G10ASbkETFbSZiadLleunPLPJiAbN27sEoWDbDD10mRbu3TpolxBmXeYZ2hxpySFFx9pnepL4nP9zCuxwGiw0rlzZxU0jxSdlFd8qdku+x4TE4ObbrpJ3ZMzr7tzeMG3u3dPelpH71rzfW2vpM2SXPv4P6fz9cTszKso7+5XkXT11Xgg5D5suzfRwRMrCbNDimxGeXk9JSEaCdFzMebww1hZqW9e7vgYIDUB646\/hEt7DPH9bPioRX9cdLpyhoWFKYpL32z6NRsLKS+NElq0aIHTp0+jT58+CAkJsTUrZpJkUtU333wT9K8+cuSIilzirj1+eIRttvXSIq7kj+uoD9kr8Ho7VwRv3L4nMePyflhTZiUmjCsFJJCSErjW4KXKt\/PWaMQhDl\/\/\/D+sbHXlWfAmJeC3Wz5Hs8tu9LY7lvV9navI3xfdaiII3nXr1uG6667zem59oQYKgNe7abcEr5H9kmZbtmqFUp2HIzU92\/kmK39enoW7rHoe\/a+MxN+l0rHw0FfYn73PI+WlgDleRUiNwbV\/9sW8bnty2WaaWyUo6XT2LZlFIrDybuqsaxdX8J48eRKHDh1CjRo1fDUVxbodf19HyzOvZJqnhFIy1I+b8By+zgzF5TE3uWRHcOfPS+p787bhahHHNbkdj2+dysNs7p+kwInIaZDsssAUWKXnWjknd4pGTNpcpEWnITUmXdWNmjQXwa36Ib1bOZ9uCl9nyPP3Rffp5JXgxvx9HT2qiih9FN0bVUUDHnoC82dORON6kc5l8+TP2xmdEYMYJCABQcqRUEoMEvPYZ7mSilSk5f0Xl5aItGjHFf4vOm0r4mJyOe8i3DBmOXmoTz76\/SLbAdj8fdGLcPpKVNP+vo5ex7BKmvYG9l75CD4bfi0OH89G31k\/IygIuLVlNcyNa2m6eIRfEpKwHA6tJv8moB2FUmcxh6QChdLnswIthdq0XNVS6lxEx8VYSmd9sWusXAK9dRX0xaIHbJu9W1E9kIQeYF\/XJ4uRi92WfbGOdt9VkHpex7BqPGAsrrv2ajz3xVlb5Idj6uL\/0g4hJb6FpTM+KS4FUNF5QE1MFV2vtKMDmCBfkXfOpVI0Tp13E+OKlupyAkl5j\/+wEOLcL397E\/rUF4sesG32bjvrwi6xHmMAf9GNUyctllyMRGqn+GId7bynoHW8kjaLxdSkPo2x69BJZ2qTLXuzMOK9P13Aa2b\/ejjiMHof7o0tc7eg7KqyuG3Yo5h51UiN8uoAJsXN9TZKnYuQLg0Rtno1HouIQFJnodgFHbLn53Szy4JETfTFopuBN2Db7Fg7b0w\/ZbV1Y5ALErycCAqg+s36GeHlQ1VSMYKYf9PRvktTh28vN9mECRMwZcoUZdUj+XMYQpaZD2ITYvFBxQ+wcuVKDOo6CPVO18OYZWOwp2dPZIaGAmkMhxONeluuRcThTzC1dWvUOXXKpU3PEDx\/NXwFXqORRsC2Ob93ldhW06RSWGeuvJ4SVdbjgmGb+XX7\/vvvVXoOyWJHq5f7n5mKR788mg8ZVo4JrCimljSenzZtmppYSd1BdzUWLgKv7y5XTv3+QLduLtf1+v7g3ubu0+AJvDzvUyjHI4TIAYztBWybTyiWNy4uTu2\/NWvWKJtnWolJIjqZM6PjBUFsdKiQfeiNP7GndTx\/5MHxZo\/SZlYStdG+Y\/\/gniGjUavbwy7JtT3FbSYYxYRO7F\/ZLr2UGNidVjXeXNftcX09gXrA8TKNOqDGQ7Owc0o\/FSvLGDLH6t2eFr0BGjjBu9VCBBewbd6kvKSYn4r6Z5p1fvXVV8qsk4RFvKdkDXQKLB5JzPTIZAHiJql7UtnZN57W0U4bRVnHo56X+XIFWKt+24T7hj6DznEj8cbgq539cqfn1Y3ldWcHfwSvALdU1SiE3zQM+xclKffDQ4unOeNfmSXXNi6Qp0UXysvn3IE3YNs8V1Fa7r9rrrlGqSwphNKTjnMO9bzGYi\/N1DB6alnW01Pd2AGVp3W000ZR1nHrVUT2RAzEJel1zdvG4LI2V7oEXbfS8wrFFUrJyRD22IpttnO9qNhmXSXEKJQCXgas2zmlPyKHL7Rl2eWPi15cbZvp2PDRRx+hQYMGWLZsmZKZhIeHFyUmnG374zrqA7clbdYjQFrNWvUKofhmxJVOVZHOxsgzPIsIC85rItTSWXM71+1Eiijo6lI1dHpfBiL6JuDg4qmoFjsFO8Z3R+lL2pd4I42AbbPrrilZ4N26DaX3bsCZkLI4VaM5gkLDnKPVBVa6YlwqCMsiwi9hd4SK6s\/ogdasrhcUnHaeM2ZN0AO623ne3xfdagwB2+YSCF4OiewLzxq33nor\/vnnH3z++ee45557UKFCBTv7+YKqU1zBe0Etko3B+vs62mKbOU5K8I4ePYrLLrtMDZvgbd26NWrVqmVjGi6sKv6+6BfWahR8tP6+jrbBy9hGzDh\/PKwmvt9yEKWO\/I1TVS7CmSr1kV3nCsuMCQWfuuL7pL8vevGd2XPbc39fR9vg5bRRcDVwzs8IDS6FVVsP4bK6FfDG3S3QfcZa9GsTiWdvbnhuZ9dHb\/NHl8CAY8LZxTXKULwNJG\/mnGAn9G2xBa8+OAqbaNr42MinsDaiB25o1xLVKpbBf9fvwYonrsSrqduVnfOXj16OBhHlnbOux+blRW8nvTDYlM3\/yCOPeIzQ\/\/1Xn2H9S\/9Ct7LbMSsjHL\/k1MG8\/66wFf\/YrI++WPSAY4JjZo1hab21UdajUOrOCVR96pEsJYifvp6+WMfC7GFPz3qM26zHA57zziJMefNzDBo+Bk3qVsWj8\/\/Ae\/9qjdQ\/92PGsm34fOjlKqqkAH\/btm1OG1MBEz8EjHRY1MUueI31dAsr6aM32RL4jC8WPeCY4DC+oHufXW8gT9TUCvhcrxIDXv1rxy+UBPPm4G+Kfwx7GvfD9NgrMenzNKzPPKoS9XZtEo6F9zsSbnvavDLJrEvbaf7NaJVUGzG0KX9v3ry5ClrGKIVcQF6X1B9CwY3voTXXjBkzQOsaFlooCeU1sqFsQ49iyPo0XJf+086a9rF8x3NtglEuYw2+u\/IpvDLnLdW+MXayzlXo\/apatarqh0SBFDWYu74zfvHUqVORlZWloi6y5OTkoGHDhvjrr7+U5RD7QAMa2XglOeg6x09TWolYacXBuQskL\/OtOyfIPly7dq2aT6O3kad9XNQEyFP7lmde8QS69957VQxf2okSPATYFbcMUp5EUoz5efksB657dugdsQNe1ufzjM5PSy\/RE9NYnYvIyRarLz0kqhl4aVLHxZecPeyftGEE+fNJo7Hyw7cwLno7tp4IxcS06pjcpTYue\/ZTPP3seDUMfnDYDxmj9NEIzOTkZLz\/\/vv5n9m3Dxu3bcOA4GBn8jD9w0PwUjjIjwsdOWTjMo4yrdVoq0uzQQJZPjr8uF0oQdfFAIihaT0FkjeCkmtmdE6wCvBebMGrU1ABnt0zqy\/Ay4+EHp1fKKi+ye2CV8wzdeorZp8KvLH34K6ITOeZ97cKLZG88EM1bHIBd9xxhwoULpRcKLb0UVzRWI\/vkkUX8PJD4WJX26ABVu3YgQG1a1uCl5S3e\/fuSExMVCal9PCimu7uu+9WP9dee636iOgsX0kMus75NDoiiON9ZGSkrUDy+kdNqLPsabZvle6lWIPXE9m2uu9v4BWwMUE4v8SkbC4U2gDe9UfCML7hLlQOPaMEWLv\/KYXmNcrhg0O1MTdlnpMl1wOiCwXkB0dfdF7XXdgUiNetw\/r0dLeUl+Alq8w+818Wpidp06YNbrvtNpV5gKlL2F5JD7pulDYL62snkLzOvXEO5SNqZMXN0qeWOPAyYXZ88i\/4+pe\/cX3L2ph4exM8OP93MLpGy9oOaytPg7bDNvuS8tLlkADyxF7zK2z88HAsjz\/+uAocXq9evXxndDPhm9n49RhL5GDc9Un\/uPAcpvdJuAGZH7sf2OLqmOBufHbGZHd+zOp52seFadsXz9rW8+7evRsTJ0\/GsvXpqBIegYY1yiH70rsw4+5WePOHHVi+cT\/mxDZH+dIhtqTN+oaUM6MkzBLhlSe2WT+vyrlWhA\/6PTOWV6iw8cxrBK8RdMKCm\/XfeOYlpWTyLgGamSCLxwFj34sCvHY2S8AxwXWWSgx4eV7ctvsAJn76G8YNvBrrf9+Eb\/aGY+agtjh68jSGLNiAl+5o4hKAzqjn1SV9+vmTrAzVSrq02Q7lFapEkEjOG+aHNZM2632hbyhD8IgqQgRC7AcDjjNqgy5s0wVcIpHU1RJccrvSZn0OpE\/Gvp8v8AYcE0ooeOmYQMHLd\/srY8fuvWhXKwfp5Vtj3J1Xof\/sn9Hx4sougdjtfOlLah1\/\/2KX1Hn39bj8fR1ts82cmL1796pkVW9+sQqjv9iJnEq1GUkHRlWRryexuLXn74te3ObzfPXX39fRK\/Cer0ksbu\/190UvbvN5vvrr7+toCl5dX6jHyG3WoiVKX\/84Zt57lWVw9fM10f70XqtFZ8TIFKSo4HMs\/FfSwOj9N3NKcBe2VBc0MdgBDTh0212Gj6H7JktBswD60\/yeq74UO\/ASrBTBi8BGNzmjgGXZhn3oFfcIHuxU91zN4Xl\/jyRhYXq01NwMhukpQP36QHq64984ZmzRysJVC\/HUgKecqikBKtO8pCNd5WySYkwFw+vuwr6aBQzXrYZefvlllYmeRviUUzBwAq28mLyauuFdu3apZNeB4nkG6AY7YsQIlXS8d+\/enh84xzVcKK8eY0o2iQ7eVet+Q+yQ0djbMh45ZSq6dNUqbrOuYNclstxY3FQ9e\/ZUdsIsC8IWqH9bH2yNhiENsTNsJ6KyHff0khGSgcisSGzO3qzqrQpbhe9CvkOd7DoYmD3QpW52drYKIsCIH+4SO5vVywgJwfiQECYXdS2pQNTSbGS8HaJSKZ04keVyv2tIV7CPM47OUJTVU4JqeXeVKlWUTMEdeMWqjC\/UVWukuM888wzGjh2Lpk2bqlCpq1evPsfbqWS+rm3btgrAsk\/9ZZT5wGtmHyudtTIjsxqMVcA5BpATUL\/44os40fYEMkMyMaL6CAd4D7VGSKkQPHbkMbTLaudsnvcJCtZl4aZvc6wN9pTfo64RvO\/seMelO1R\/UKdLI4syZcpYzrtZvelVqmD6oXBEzc7G\/Ad2gGBmvW+CgvBL1apYU6EC2mZl4Z0dvJfh8qG56NkViJz5FHr\/OhzDq7r\/asu7+cGkyZ8Z2yyOFGamfjT4EHtdUl5GNxGKUadOHX\/Za8WyH\/wATp8+3YWL8peBuIBXdJcMhG4W2JxUOHXtH1iY3Slf\/80oLwFqFuqVOlwB7\/XzrsfsjrNz20tENBLUOfC2gwexrsrHmJk10iUBKAEyPmQ87sq+C0uzl6LakWo4GH5QJeFemt0RmSEh6JCdrU6UzDXIwjGRVSQoSNWsilm9QRkh+C4kBH9EnqWser29FSogO3szkpGsKP+S7CXO5jdvzsbs7NmYXmU6kBqDtlltERMdjUSVm9i1SJv8spNDMFJeqW11Fmb+ZDPwmpn8+cvGKy798Odzbz6BlWwQTq54ZMi102dycLLjY5gYH+PMS+RuEThws0wIYrzP4OPZUVHIjpqVm8YzBtEpigtFVIds4lGxm59nZTkTgPJdvN80LExdI+UlNYzKzlag5bW3Q0JQJzvbCd6mYU1RM6tmbirvRKwO666eZz1+ftplnQXlWyEhSALw85Ej6rxoVeimt3PnTtSsWRP7Ku5THxumL9WBy2el3uaozXgr41t81\/BtZGSEICZ+q8otzrPy0oEDFcBZ97Ojn6FKlda4qUILS\/C6M6C\/0MFrleJEzx5YkA9GsQKvDFCXMvMa2baYHreZWlJZTYon8Lbf8CWQngAeKglAFm7w1q0PISSkFB577AjatXM9T7IO2dn0oCCEHzqEw1WroveRIy5A1Pvz1+m\/sOzMMjzX6Lm8y4QuBUZx6p18NiM0VFHYyJMnMfXAAdTPcwQwG5fOXv9U+ScMqDUA83fMN7D31fHgzp0os2OHC7tOziHj24YYMMARtG\/o\/kfxfsX38WbGm+gc3RlR2R2xMmQlTm3ahC+\/\/NLpv6z3w2ikr9tJ83dmFXjvvfeUtDlAeR0ZBXmsKGig\/mIJXitAMkvg57\/usWVNpbPNixcvVl4xjIpANk8m5T\/\/+Q+uuOIK5+sIjrVr92Hs2IuRy2XiySdPKCqll9WryyIy8gSCgrYpCrhmTQWMGxeCzMwQ3HWX4yNA6ta27Ql1RiXbfKrOKaypsEZRSkp4y4Y1RYfsu7AijzEfuSELSfFQrPUnnxyx\/EhLe3zvtCrT8HbI2\/j5yFnfZj64sGxZJeQanJmJnuXLuz1rv1D2BdWf1WG7EIetcIQRKFzx5w1nZ2Ts\/8yZM50mrOQq6N8sREQcSLifWERwJ8nF9MR4MTExStLONgTEciyhDb2YxoqLaHFyyLdlpGEnY4LZmVcEVnfddZdKDMXy4YcfOiNA0NOHAitGypBiFByREkdFOQAp5aKLGqhfSaHXraus7rdt66DQq1eHKerNsmXLVgVeEVjt2VMeI0ZUB6LTkDloHBCzAkN\/moo6mzvkCiXCVd2pUw+gfn2HC55Z0dt7reZrinJ+s\/2bfFW\/Cw3FiGrVMGnvXnQ4dcrtnp1e5WNMD2+Nz7Paobubc7mdjc86RvCmMdXx2VzobpsxOZKDedDtloI8z480f6Sw\/xKWRpKGSWIxylDo1\/zTTz9hyJAhziOOsM0DBw50JhfTU+eYgZdCTJbi6tNrCV4j2ywTS6BZRcgwW2AuAp3YWdq3b4\/Zs2erCZcNNm\/ePOWjKsWOgIngpEDowIEDOHasulr4jh1dAW7V3o03OoRWajPHJSOt\/grEpMzFyJFZiI5OQ0bDDKWCioNBeZvXoN6\/myvcnHcu\/zzf0Flv\/aFDiI+OxsysLBfBm7EysXFjWBiWK9Fb4YsRvEm5h\/nERHvtmp0YgoLsPctaBXmefUs4q\/pWe0MPDmeW0pPvEmqsU94ePXq4OO8L22wG3mPHjiluUIoZ9fVnLsZjlkCGGSGbyzODMXGYnSUVPbGexlMH75w5c1zOI7pAyK7gqDD1kpMdLHZkpEMQ9UXNLzAsfJgCbyd0UsIo6mpHZo1U+XT1\/l0efrmSfI\/OHp1vKqTehpo10Ss8XLHDd2gCMnlAdMlfZGVhU3a2TzJQlATKK+A1Ul7jREtQuV69emHp0qXQKa8kdhdpvFBZoexMOyvXrPZysQSv6BOXLFmihCD0rf1z204MeOgJzJ85EY3rRTrHa5UlUJ94nkNE8qyD1xPbbDWphdHfmrXpwl6X34Pp4dPxfoX3le6Wah4Km4YeGIrLDl2m2PCTtU6ia8Ou+YRV0rZre+UxIo9Fo05YyvsVKuC9ihWxOiwM03\/+Gf1q1FAqLRbdhVGnCLzOIkEARALN0LyMLMFYVsaYWnY+sv5UR6e87JeRC3zttdfUmZh7isXqzMt7dP9k\/DU9ECC5R8YD08+8rGvGVRY78Or6XjqyDxs+An83uRsZm39FWNo3ONZ+aD4LKybYHtWhDP51\/yAVKZG+sWSHGSxNL7xOpbc7gRUFTBQIeTKqKOp6maEOY5A6p+pgBVZgSsUpWLR3kRKADbp4EBqHNcaMIzNM970u2JJxPFKxIjZnZ+PzEyfUM5mhofg2JARtjh5FaGamEqjQyko\/87Ge\/hHkkeXrr79W4CU3VBLB628fEj0Kiz\/1zfLMq+sUSX3lfNF4wFgkj+zjtWOCrjbSKS+TlYnRPCeGulueZWne5w68JaWePuarrrpKff3dsYq0oOJHkVyMfBjJJZUkyutPACl2lLcoJs8IXto204QvYH97draHDh0K\/gi1lXjPOtssjgfU41IYKCFfCd6JEyc6YxuzjYCet\/A7OQBeizkkgPlT0orYw9KYXbctpjtgLGLzDVeu0zzSzPhd54IISHoNSQigli1bKoOOAOUtml1UrMCrB1wjCzd58mQVPZHCAStFdtFMW\/Ft1WrBKbnW3QFlhMbruicX6+hBwXXwylqxDilySRFY+dPKFxvwGh0TjIHXuKkY\/JvBvt2pZ\/xp8s9HXwoLXvbZGLxPWGBhm8USSKzYLnTwXvC2zUbDd6PghPcH\/nskwtrdgzcfao8f\/jqk0p6Elw\/FRw+1dsZtPh+A8ad3WoGXjhjUGRsLzSOXKxMN3xR\/pha+GaH9Vi4Y22ZP4NX1vFXCwxE791cVbJ1lxHt\/IiW+hddSaOMy6Do9Xe+mG+SbhVllO3auu1t2b98tLK2eB4nXKHSSpFjn46hR3MEbsG2293Fy64xvpLwEb9\/7R2DmixNwIqicysnLQOtb9mb5BLx834QJEzBlyhSVg0aSnVEYM3r0aDz99NNqVFKHvw8fPtz2dbZpVbx9N9vS5QPC1vLoMWjQIHW8oDqnoN4s9pbPvFY+C6s8V0o7bZqZZ3ph2mxq3unpedo5aabNLnrugG2z9arlAy+tTsRyxeyx6EbNsf7iwagSEaFY5V2HTirW+d3BrWz5+NrZQFJH1EtMX8lseTRQ0I3N5Wxo97o3QPL0bnIFzBTYp08fJdDTjSY4hwz09uijj4IhVM51kdhL8kGhn7JN02aYuWR4YdpcoOfZN92UJ2DbbG\/H2PIqsteU72uJLbVuF823kE1ltA9vr5tFB7Hqtad3S1tCfQW8OnvPtpkkjAnBzkcR8NIHw6ZTUYEopz62glBuM8obsG32vGMswWvchNJUy1atUKrzcKSmn\/XioWmk5Cny\/Ep7NfTws7qBx7kAr513W4FXRkf99ZYtW5TXSosWLc5pxEZ\/jrtkZ\/UDts12Zgmw9Coi6yfURJpihsB7U35DVJUyLs74wxdtQMbBkz4DsNF7SdQhdtljK3baDtts993SlpHyGqfd6Ehgb1kKV6u4C6wKN3rfPu3Pc2kJXrMoklbeQ1bXCzKNkvVcB5pVFEq27yuBlQBe3B+l7+4iYLKOGdssjuSegF2Q+bHzjD9vODv996c6\/jyXbp3xxRVQn8zXVmzH3B\/+xif\/bq3UQhJlI7597UIHYjdj1SUZsoQ2YV90m139GTvXrTZGQd5tBl75COghWsxy+BblBvXnDVeU4y6Ktv15Lr2OpJFdpYGpS6BV0PWimNBAm6OZekMAAAFjSURBVO5nwJ83XHFbO3+eS4+RNMzSaxS3BbjQ+uvPG66o1uKCN4+UidUN4d0ZNhTVQgTaLdwMXIjgtZqxC8Y8Up8AbgARvBjPvE9\/vDnfXAXY5sIBzpdPF3fwBswj7e0Gt6oiM0urMjUbYvYbb+CD347hxhbVlVUVVUWX1ChfaIGVvS4HanmagfzmkYxWbc9Mw8pxwtM75X5BnmdgP\/4nRSccAfNI65n3ysJKVwl9+ONubNp9TOl7fakqsrtJAvWsZ8AIXvoLM92LnZJjYuAYBPsWYgV5nn3T\/ZwD5pF2VsrCSMPqUTHS6NwoAh0uCcd9yb\/ijbgW+G7TARf1kb1XB2oV1QyUBMobMI\/0vDu8ZpubtWiJ0tc\/jpn3XgVSXzn\/FoVjgufuB2qYzUBJOPMKeDm+QOhX833uFdtsZv0UgI\/\/zUBxB68\/zag\/z6VX4PU2ubY\/LcKF1Bd\/3nDFbR38eS7\/H8g2mXtxs02LAAAAAElFTkSuQmCC","height":160,"width":266}}
%---
%[output:8457fe63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:2ff15fcb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:3a6b3790]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:8baa350f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:535ffb60]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:4571544a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:1f5daac6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:360bb907]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:84c34891]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:4252e3ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:8c6d783a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:992c84b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:12b10a33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:645e0844]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:8934ee95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:25317e8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6e2714c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1dd95874]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6cd8abec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:77b11528]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3f6f6768]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5752878b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:174e4178]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1c613be9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:21f652b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3839e252]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5f33d8e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2b810821]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:02fb05d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4c1a1aab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1eb51541]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:22a12065]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:512e6906]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:122d0140]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3938b5c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8da9cac7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:502a80f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1806be07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:26d12cd8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:00a31fde]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7fd373ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3f00a6d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2c216b39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6369f23b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:293bef6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:557de28f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8c8cfbfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7c248ad9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:44ac16be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1a829cad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6995da25]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:25778c6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:07eae6fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:000df71f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:96144a9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:32b90bd2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7cc8af50]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:82f9b87d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:12105f38]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQoAAACgCAYAAAD0OyL1AAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl8FEXWf0QMCQkQAoGEy0QB5ZBTBNkILl6fICAiIocgAgKLAovcKyooLMqhHAoip7oiivIJCuuu+IEcKwIiglzigpxBbggSIMTv9y+soaZTPdM9Pd3TPan6\/eYHmamu49Wrf7\/36tV7hX7\/\/fffSRVFAUUBRYEAFCikgELxh6KAokAwCiigCEYh9buigKIAKaBQTKAooCgQlAIKKIKSyFsVvvnmG+rYsaPfoGfNmkXNmjXz+27btm3UrVs3OnHihO\/7Xr160dChQ+nkyZPUvXt32rJlC7Vs2ZLGjRtH8fHxhggBk9fu3btp\/vz5tHr1ajp06BB7LiMjg+666y7q3LkzpaenU6FChQy1pyq5gwIKKNyxDmEbhQwoevfuTYMHD\/bbnAsXLqThw4f79WsVKLKzs+m1116jd999l3Jzc6VzKly4MPXo0YOeeeYZw+ATNuKohkKmgAKKkEnnzgdlQJGZmUnTpk2j4sWLs0FfunSJnn\/+efrwww\/DBhQAhqlTp7KPkQKgwAfAoYr7KaCAwv1rZGqEIlCkpaVRcnIyZWVl0dy5c6lmzZqsrSNHjjDV4tixYww89u7dy763IlGsWbOG+vTpQ+fPn2eb\/8knn6THH3+cypUrx9o+evQoLViwgKZPn86kjYSEBPZ\/gJgq7qeAAgr3r5GpEYpAUbt2bbr11lvpvffeo7\/\/\/e\/Uvn171hbqdOnShRo1asTUEdgSrACFVkJBX48++mg+OwTsF5BiuMqDOqNHj6bY2FhTc1SVnaeAAgrnaW5rj1qgeOyxx9jGFDcl1JBJkyYxW8HZs2d9KkioEgUkE7S1detWBkxvvfUWpaamSucJ6Qb98LowtKakpNhKE9W4dQoooLBOQ1e1oAUKGDHxKVGiBM2ePZuKFStGAwcOpC+\/\/JKwSTds2MA2thWJ4qeffmInKIcPH\/YDJBlhLl68SCNGjKDFixcztQQqUZUqVVxFQzWY\/BRQQBFlXCECRYsWLZg0MWzYMKZuvPPOO1SqVCm2qaEGYJNiw4YTKLhUEoisr7zyCutTAYV3mE8BhXfWytBIRaCAD8SYMWNo\/Pjx7MgSkkT16tWZmsBPQmBQDCdQtGnThsaOHUtFihSRjvfChQsMuJYuXaqAwtCKuqOSAgp3rEPYRqEFCjhLffrpp0zch9NVmTJl6IMPPmAnEiNHjqQZM2Ywe4UV1YOfouzcuTOojUJbV9kowrb0tjakgMJW8jrfuAwosIFxygH7wHXXXcf+xVsfhk7R8SpUYybaA+gsWrSITVjPR0Lra6FOPZznj1B7VEARKuVc+pwIFHzji29xDBs+DLBX1K1bNyxAgTZXrVpFPXv2ZD4S8KNo164dU3FuuOEGRimMYebMmcyXQvlRuJR5AgxLAYX31izgiGVAASco2AU+\/\/xz9uwtt9zCTkDgkBVMoghGHu6foTwzg1HK278roPD2+uUbvQwoUIn7TuD\/OA2B7QKSBcADqoKejSIYeURHLhgq4cINu4O66xGMct76XQGFt9Yr6Gj1gOKrr75iqgAKTj+efvpp9n9ZffH2aLAORaBAXX57FB6YK1eu9LmH89ujsEtUrVpV3R4NRliX\/a6AwmULooajKOBGCiigcOOqqDEpCriMAgooXLYgajiKAm6kgAIKN66KGpOigMsooIDCZQuihqMo4EYKKKBw46qoMSkKuIwCCihctiBqOIoCbqSAAgo3rooak6KAAQqsHn+Uzhy4RCUqxtKdg8saeCL0KgooQqedelJRIGIU+GHhKfqs3wFf\/3cOKmsrWLgGKA4ePEgff\/wxtW3blipUqKC7AHANPnfuHIvUZHcEZyf7woSd7E\/1ZX2PR5KGAAmABS+12pekB6dUtD4pnRZcAxTclfj9999nQV\/1Sk5ODruJiAtNcXFxthEGDTvZl9P9OTk31Zc5NpWpFFoaaiUKgATAwq6igCIAZZ1kcAUU4WFxJ9fMjr70VApZXwCUX9Zl0w2NE21VO7AyCigUUESVZGbH5nVSutVTKZycl2y+CigUUCigCFGYsWPz6qkUdvRlZtoKKBRQKKAws2OEunZtXplKYVdfRqeugCKKgAInR\/gYKbDYHz9+nJKSkmwHCif62r\/uPPMpSEiLoaJVL0TNvHACiI8Cij+4Wp16WDtlAUAg0c\/69euN4ISq4xEKNGzYkKVbKF26tGOnfcpGYZI5nEZxK\/1xoAVTlS9f3uRMVXU3UgCgP3nyZILLQJ06dbwBFAjQisjNGPShQ4dYxqmHHnqInnrqKb\/ckQiFtnbtWhY7EenqUA+ZrREuHjEa9YqSKKxJFEbp58YNocYkp4C4pp4ACiSh7d+\/PyHHZMeOHalevXosyez8+fOpUqVKNGXKFN9bbM2aNdS3b1+CyIS8Edu3b2f1mjRpQqNGjaLExEQpVYwyupW3rlmGdLIvq34URulnlgaqfuQo4DmgQPo3ZJpCJOemTZv6KLdt2zbq3r07de7cmQVrPX36NAMJSBEIuspBATkf+vXrx75r3ry5Agod3rMCTAooIreh7erZU0CRl5fHJAaoE8hRmZyc7KMLzxcRExPDwr9v2bKFqRhIU4f0dbxkZ2ezyM94dvTo0RQbG5uPtkYZ3cpmMrugTvalJAqzqxP99T0FFIGW4+zZs0ySKFmyJAMK5LiE1KFNZQ+7BaQJ2CyQeEYEG96+AgplowjEa2JyY7FesLtBXoaTqAEKqBQAChzLQZJAOvt169axBDApKSl+a4Tfli1bRvPmzSPkeNAWBRQKKPQ2Nc8zgsuCQ4cO9VXjPMNTJ3oZFGRjjwqg2L17N\/Xp04cBAlQTZMkGGGByMqkBqetwEqKVNpREcY1FrKg6RoHWi5sJvIMM7Hp8pfebF+cqjtnzQHHgwAFmcwDST58+nWV+QgkHUOB0pXXr1pSamipdZ2ymrKwsKlu2LMXHx9vKC072xW0Uoc5t48aNTKpzUhQX85Zi\/OKbHSdk3bp1o\/r16zO1FGvF1YdNmzaxFwYK6sBQ\/vbbbzMbl7Yd\/B3sJWOGCcQxlytXjvWNzO4TJkwIGN7ATB\/hqsuBYs6cOcyPQsYbiMlid1wWzMe0C\/eOHTvYUemVK1eYPaJatWo+ugQCCqOqBxrr2rUrY3pZuXjxIuG4FpJMkSJFwrUmEe8LA7AyN2wygLdTQIH1xGkYlxBl6gFndJ52kG9SPkYOJocPH\/aNW9aOWE+bwtAMA2jHLOvfTHt21+X0mzRpEkssLeN7uODDRmh3MQwUMEjCRwLMmJ6eThh8xYr+EXUgAlo1ZsKzEG8hPYkCb6WjR4+y3+0OXONkX1hoK\/2BqeDYZhQoXvliL+0\/mUOVkuNo6P35bUaBGI9vMCQ3bt++va+qTP3hLw\/YFgYNGkQtW7b02RkCtYO6opoq5kjlHdauXVvXOK4dv15fWvCye8OZaZ\/PGY6ONWrUkPK96yQK0ZFqzJgx+YyVIAAmpo5HzbCCf12nbBQLNmRR3wU7fJ0PvT\/dFFjoqQJcGoCjHQcQMeGxdmPzzasV+\/U2NR9wIJUHwIQXmQhgnDe14MO\/hxOhUYANfXXNP+k5GwU3XFavXp35QeiJOqdOnWIOVwhTJ3phKocrY0ziFFAAJAAWvHRokEpvdLimQgYbrXajautr1QNeX\/s9AGHIkCH06quvUpUqVXzNyABHNibxyJRvdD2g0AM3PbAKRgMnfvcUUEBvBjgsWLCAiY1Aa22BCtKqVSvmSLV8+XIaMGAAZWZmUrt27WjPnj3KhdsgVzkFFFqJAiABsDBaAp1C6In8sEPAeCiqE0YkinvvvZd5\/4pSitiHVvoAUOBO0YoVK5hxlBtYsemURGF0hfPXC2qjgAGlR48e7G6HXgGAcMu27FJYp06dmP5cvHhx3TaMHu9Z2UxmyeRkXxiblf6M0o\/TADaKNXtOU2blJFNqB54PtsG5KsHf+Pv372e2K0gOKJxXjNooAhnJtWMRDZa4SgCQgX0E\/8cJi9au4gUbRZ96r1G1cvUor9h5um9kuu22Odk+CQoUZjdXqPWNMrqVzWR2bE725TRQmKWFtr72BIGDgggEeqccfLPyTY5nuKTBvxONnlwVQT3Rl0J2QoJxoXDHLIxh37597G+9Uw+xf6t0CefzfE+0yHmZ0q7UZE1XapzgSDBd7TwUUARYWQUUgdk+kFFRz2tStFfgFjLe8njrAwCgnqDoHYFio+O+kVi0RkitjUIECjwn86OAD4fWThLODR9qWzKg4G3ZnfBHAYWJVVNAYYJYIVS1w5ColSi0f2uHic2IOnp3kEKYVtgeCQQUdif8UUBhYhkVUJggVghV7QIKfo0AQwpko9C7QxLCVGx5hAPF8NZv0KXVGSwmKC92J\/xRQGFiSRVQmCBWCFXtAgrx1ENUY2ROW1Y8PUOYsqlHtMej\/3ppH53dfr2yURhxfHFy8zrZFzjISn9GjcGmOFVVjigFPOVH4RSljDK6lc1kdi5O9qWAwuzqRH99vieajJ1DlWrfTvE5Z+iv1Ur6jkdlOUrtooo69QhAWQUUdrGdatcIBThQFB0ykwrffBt7pHd6HPXOiGOZzJF+kBe7T0EUUCigMMKzqk4EKCADiqaXiSbfl8RAAmDBi92nIAooFFBEYAuoLo1QQAYUt08\/QY83SGQnIAooGjXSpaOT6oCTfSkbhZGtU7DqiEBRIrkOpa86TzUWnSFIDyUqxtLqCUeVRKHHEk5uXif7UkBRsEDAyGxFoGj8VQYDChQABYqSKJREYejijxh8pmmJoyw5kxvjKhjZFKpOfgpwoGhw43iqvfXaNXwARaU\/JfoZM+12wFI2Co\/aKLRXxQEUW+aPiFqgKMjh+sVLYVyiADDgePSXddmOOGApoPgDKC7sep3yfjtEMUXLU\/zNA9i3blY9tMFnbsz7hU4uGR2VQFHQw\/WnPziViifXoYRjuT4bBYBixt4cOpyTR+XiYtiRqZ1FAQWC2h5YROc3D\/bROf7m\/gws3AwUWoni8YwztPQ154Lr2smU2rYLerh+0Y8CxszG6y9Q3vAyNDs5xkcq7l9h17oooCBiIAGw4KVIxUcooe54VwMFxioGn4mEjcJL4fq56oKQezz6lRhwya4NZqVd2fEoDJo4Iv22Tyna1zTB13yr1FgaXa2ole4CPquAQiJRACQAFpGWKMxEyjbqAh8uTvJauH5ZfM1w0cKudhRQSChrlNHt2rywUVw+vp6uL90w7DYKmf1Dxlzi3BZvPW0qUrZR+oWDqb0Yrl8WgSsctLCzDQUUEQQKo5s2XMZMPftHMKB4dvFeU5GyzQLFqQ9HUe6xfVQ4JZ1KPvqCKX73Yrh+DhQI8CvmLzU1cYcr63lmQv3Y3uYCbWt\/s29EbUrtoBdq3WHbCAuU6mFm04YLKPTsH8GAQitRBIuUbQYozq2cR8fe6OYbQsl2L5oCCy+G648WoIAxE5+vhx+hrNoNfWvYusxJGlXjRgUUnAJWVA8zmzZcQKEFJ9g\/8n47mO8oVtufnaoHQAJgwUuxu56glL5Xc4EaKV4M1x8tQAFposz2HGbMFIsyZmo41wpQyDYtjJZ6Ra8vM+oL2hbtHzFFK0iPYrVAoad68L5RX\/T5sCJRACQAFkaLF8P1I0vZsGHDWG4RL6seKdtzKOHYFb8TD6ybAoowAoV203LHKnwv09llQGFWfdG2HVslW3oUG0yigOrxcOoaP5DBM9znwwxQ8DFd2L6S4qvfZUrt4MvhtXD9\/fr1iwqgAP3heHXfP3+g0kfP0fGyxWjx4\/VI2SjCDBSyt6aezi4DCrPqC9o+t\/qvFJNYmHUdm1GZKOZaOj9+FKsFCiRgrvPyf3yJhL9\/7o58\/h54BhIRJIsNPxyiJ4Z+7KhnppfC9UcTUNz5r93Uc8LXPlYGUOQ9m65sFOLmtqJ66InWejq7EYlC3Oiy9s\/88xm6cukzXameSwTBJAokEh5QfVM+ieJITB1Ky\/ueNu0h6jWtkKNAYVRV0atnV3Bd9MfVi2Dh+q3Owc7nZace6K\/nhFV0579+8nW9+r4qdHTCvcrhym6g0EoUiffUoYTb72Hdnss+RwnJlalYjUG+Ych8LvQYRiuBaOvBZlGkYtt8LuOrV4yjfYf2+KpXSo6jzJtK+qktB3NK0\/rTN1Pb1LUKKP6glJh+EF\/xcP2NAtxItnOzW2mbA0WntN6UmnTtKLR0VjZV++GIr+kdtdLox\/faKKCwAhRGvRt\/fasuO41g6sFNifnWV3zzyxZfz8AZDCh4W4VLNaJCSfXpeNG2lHR2If3+i39GLD2GA1hUiDuugEIACr1w\/VY2bSSe5UDx8q+tqebF8gGHcGVwI6o7NtO2YUa1HwVA4pUv9vmIt\/DuVdS43Hm\/0wL8qDVQyqjN73+Ipw6oB\/uA9iQDm\/66ohVYM+IdEiOrmJtQmwqf32Kkql8dL6oepidZwB4wAxQlu9SgirMfsI1CYQcKWTZzZDLv0qULS0evV0Srfe2S3zA\/A5SDOaXof8+2p37pn\/q+E9UBLjHg9wpxJ1h9FPz\/\/Q1XxbNDOaXZv6jDy7GUXgRxHv1cPoH+rkoTeoUbDS\/smmzbYlhpWAGFFeq589moBoo1a9ZQ3759qWHDhvTYY4\/R9u3baf78+dSkSRMaNWoUJSbmF+uxTHhTf\/juZGrUsBGl\/P6N38pBD2+YtCvfauLNvfDH69j30NPNFC6yG30GUkMwMDHalh31FFDYQdXIthm1QHHq1CkGEqVKlWIZqTkorFq1inA8he+aN2+ej\/pGRP\/ILpn7e1dA4f41MjvCqAUKTAwqxowZM6hZs2Y+umRnZ9PAgQMpOTmZRo8eTbGxsX40M2rwM0voglRfAUX0rXbUAsUHH3xA06ZNo7lz5xIChPACuwWkiQ0bNkjTyyugsM7kCiis09BtLUQtUOAMe926dTRr1ixKSUnxozt+W7ZsGc2bN48yMjKURBFmrlRAEWaCuqC5qAYKTG727NlMzRCLXgwD1FEShXWuVEBhnYZua0EBhaCWFASgyMvO9d31sIsZFVDYRdnItVsggSKQ6qF36uH2I8lALIQjXZT1p2+hBoW30R0VfraV46IdKApyXg\/RM\/N42avuBaWPZvvx09GJ99C9\/erYxmNhdbgK1ZiJ2f20YjDt27KIGv\/5Ebry20Ff7Ep+r4JTIDc3l3JzL1NcXLwfUS7t3UN5F646TeENjlua8dXas78RC\/PAqRzWLtydeYEfBp6R+UeY9bMAMOAZeH4eLlSbpuxrzbrJrJxErQpPpLLnP7dtEdFwNANFQc\/rIQLF24Oa0Or7qlKbd7+jW7ZcdSjcWTuN4MLtmSjcoR6PYrJG4yno3R41ErUJXpz1cucz5y0xiO6u\/4yj\/\/60km10AMldDf+H\/nM4keqcG2N4cwMY8OnQIJUQO0IsZvxEcsp0pbiLPxKd2Wi4b1TcfPA26jlhk6dujxqdYEHP64FLYTVzyjNAwJVyWfFUPArucJWWlubnhRnM4SocQKG9AWo2apOYI2Po\/VdPZcwYWT\/O+hMN3dmdgQTAQltEySj3hL\/nKa9b6IZedLp4eyqWNZFijpmTQLacakTdR613FChUXg+jUBdaPb1r5rLWPBczc\/ny5TRgwADKzMykdu3a0Z49ewy5cFuVKEA8RKmyErVJuwCQXk5+P46KZq8I6r695PKzdDSxBXGQCcQaAI3\/\/rSK1uwuQuWuZFHhsnEENejOu4fRkSNHKH7LGIqNWZGvie+L\/Y1JQ7KLZpsPV6Cerx5yDChUXo\/QNr+Zp2RAgehWJSrGslSCYvFcKDzZpbBOnToRLoYVL15cl07hAAqjV8qNLhZXc7RveKgoh\/64fIa2oKqIYfWMtK\/NHfqnm5Loo+7VGFAk7vqCLm8YToVT45i9ZXOp+swoiotskFgANBcPfOwHXtsuPuJYhCuV18PICluvIwOKR7ZcpEv3F6MlWZf8OoB9AmBhVwmrMdPKIK0ChTYXJyJCGXm7BxozB4rky2vpys6RvqpDdnanT7L+xP4OtR8tUKCt\/k1TqUutOILqdmHJK7R7zT9p8fmbaWaJrqwvrVojBtDZfLgide4+2LBEYTZAsEgnN+f10FvPaIjCjcC6\/XddplOdS9LE87m+qcb+Np+Qe\/TJGr2sbMGAz0YNUGg3nsyoaJaKouH0919m+DKJ3fFRQ79YlmbbRX0ZULS5NYmeuyuJAQViZqLIbCey\/r7+aLBhiSKUAMFaoBg+fLjutOGu37791RMnFG7L0H4PyWTIkCH06quv+rn881MO3D4W29F2aCZNYDQABXKOcqlh2o3b6Vj1OLru8vcU+9s79EB6S\/pbgxdDYUVDz0QNUGglimAJc4xQR3bCEi7JRQYUEx5Kpz9XJD+gMDJO1Fn3fx8blijMBgjWjsHNeT1E8GjatCkbeq9evcjL4frTH5xKaVdqssQ\/sE+UqHg9rbvyBf3w+Gzf0oxo8CI1T29plF1M14saoDDz9jVKJRlQhEtykQFO\/6ZpzEYhShRGx2pUdUN7ZvObaMfg5rweADHcNxo3bhxt2bKFOnbsyNQxLwNFi5yXGVCIBYDxa69ldLLKTqqbUt9WtQP9RhVQGN1URusZkShClVwiCRSYv5kAwTJ6uTGvB5cmOnToQAimK\/4dbUBRq31JenBKRaOsbLmeAooAJNRz7jJqNwi0OjLJZGKbDEckCstc80cDbsvrwW0bCNUf7UABkDiz\/xKdOXCJqSN3Di4brmWVtqOAIgSgCMeKyGwqMGY6oXqEY\/zhaCPceT0CSRReDtevVT0qNU6gWo8l02f9DviW4c5BZW0FiwIHFGZ8LexINiRuMK1kou0v2BGm+Ds8M7k+7pVNEW6gAG31bBReoYnIH9zupAUKqB0oPyw85atutypSoIDC7ImF3UChfSuL\/RU69pluMmOZQTISLtxWpQo7gEI89QBwnjt3jrjNwup4nX5eDyi42rF6wlEFFHqLYmXzmj2xsNJXKEwl9ndlx0jdZMZoW3vE6bQLdyjzc\/oZPT8Np8cRan8yoIDa0XnxTUztUBJFgPRvVjavWV8LK32FwhyBJAptjlPtEaeTLtyhzM2pZ3Aa89Zb17Ks4WjUi2oH6CUDCq5iACREGwWkDK6S2EHrAqV6gIBmTizcBBSylIbiEacXbRR2MHQ0tRkIKDDP1eOP0i\/rsumGxom2GjLRV4EDCjOMFEmgCKZ6aOdhxuHKDA1U3chRQAYUdksOerNVQBGADyIJFFpjplb1UEARuQ3sVM9aoKjZI5buG5nuuwfk1DiURBGE0pEEClwKM+M9qSQKJ7eNM30poJDQ2SijO7l5newLJLHSn1H6OcPiqpdwUEALFBktrqc2b2YoicKIw5CVzWR28ZzsSwGF2dWJ\/vpaoGg0Mo4a96ikgEIBRU6BcuGO\/q1ubYYcKPrUe43q1K5DlTvnhXSz2Noorj6tjJkuNWbywDVGF7mgqB78IlrLli3ZVfL4eP+0DUbpZaQev2SG6+oyfwy9XCNi27i1KsucZ6R\/cU3r1KkT8kvESF\/B6iigUEARjEdc8zvfuIi9unr1asNh\/0KdAA\/5hzy6lSpVygdMHCjQvh2gpYBCGTOlvGvFJlIQJArMcdCgQTRhwgSCByY8LnGl3I4ihs5LT0+nqVOn0ty5c\/1C9imgsIPyQdo0yuhWNpPZaTnZlxeNmU7l9QBttJtyyZIlhMx0ZsV6cczlypWjESNG0NixYxn4iK7eIj9WrVqVunfvTtoYngoozO6oMNRXQOGt41En83qAvbQpAvRSBgRiRe2YeRuHDx\/Op8agLniSA5H2bxl4hdteEkz1mLE3h+X3KBcXQ70zrgZjtqsoG0UBtVHgnkCo0ZEikddDG9DX7Ntcb8xcwhCNlbIo4LIXmRFjpjbyuJmNHAgokNfj+R2\/+ZpDuH47wUIBRQEECu3NQ7PRkZzO66GXpFi2yfHmh01BG+af2ze0dgYZAMjmx8cgGjXNgpUZkEDdQEABkBCTAHkuU5hZYvD6SvVwTvWwGstAa5vQrnm483pw3tDjLYTj50ZNPaDQAzdt8JxgUgLsGhxsIgkUWolCZQrTcIeTBkYn+3LSmGk1loGTeT1atWpFw4YNY1wgO4KU2RISEhJoxYoVLFw\/BxGjEkWgF5ZWfYkkUIAesFFsPJ1LtyUVtlXtQF9K9SiAqgembCWWgZN5PTDWbt260TPPPCPNGsY3NpdiRINlqVKl2GkFpA38X9aOVn2RGS05i3Bg2L9\/PzNywngZCMRCla5lUrZyuPqDKkr1cE71sMrAeN6pvB56KoN283KJY8qUKewnrorg+X379rG\/9U49UB\/qBAeWQP4ZIp\/yXCF60o5VOgc79bDavpnnDUsU58+fp3feeYcdIx06dIgR9aGHHqKnnnqK4LnGiyybOTKZd+nShSAS6hUFFN4CCqyj3Xk99IyYWh4SpYJVq1b5GTNFoNCOmftRvP322yz\/6XfffSd1rBL7E42azz\/\/PI0ePZqWLl0acM+FGo7Pc0Bx7Ngx6t+\/PzvLxsWtevXq0datW2n+\/PnMtRUoXr58eUasNWvWUN++falhw4bMQWX79u2sXpMmTWjUqFGUmJgoJaoCCu8BhZk3kqyuHVG4ITWIEoX2b+04wHeoY9Zxy+rcjTzvOaAAYsKDbdq0acQTv2Ki27ZtYzpg586d6emnn6bTp08zkIC0AZ2RgwJQvl+\/fuy75s2bK6DQ4RIrxlOjQGuEQZ2qYxdQcEcpzCOQjcKoxOIUPWQgxm9Uu95GkZeXxySGtWvXsujGycnJvvlAHYExJyYmxpcUFirGjBkzqFmzZr562dnZNHDgQPYsRLXY2Nh8tDfK6FY2k9kFd7IvjM1Kf0bpZ5YGdta3CyjEUw\/xqFZ2zGrFIcpO2qBtz0kUegQ5e\/YskyRKlizJgOLTTz9lUofWqQV2CyzIhg0bdEU8o4xuZTOZXVgn+yqIQGF2PQpa\/ajlMzW6AAAI4klEQVQBCqgUAIrBgwczYyV0PaScnzVrlp+BEwuM35YtW0bz5s2jjIwMJVFIuN4KMBkF2oK22bw836gAit27d1OfPn0YIEA1KVOmDAMD8SKNuEjBjrk4UWA0bd26NaWmpkrXGJspKyuLypYta2vQEv6Gd6ovq\/1t3LiRgXWoFnYvb6hoHTvfE3PmzCHYKGS8WLhwYcLH7mL4eFQcyIEDB5jNAcag6dOnE67hcqnBKlCgna5duzKmlxUczX7yySf08MMP+05a7CKSk31hDlb6gxci1kQBhV3c4Hy7HCgmTZpEpUuXlvJ9UlISU\/3tLj6gkBl6RB96PpAdO3awo9IrV64we0S1atV8YwwkURhVPdq2bUu33XYbkxhk5eDBg4Tza5yu1K9f31b6ONkXJmKlP\/6sAgpbWcLRxvmeBK9jP8j4Hhfg8LG7GAYKGCThI4G3FgYGlKtYsaLf+BBIJFRjJhgdto7169fbPeeobl8BRfQsb7DLcJgpXtr42F0Mqx6iI9WYMWPyGSsxUEws1ONR\/kYFYKhingIA2MmTJyvVwzzpXPsEB4rx48frqtkVKlQgfOwuhoCCGy6rV6\/O\/CD0dKJTp04xlSAtLc3PC9OIw5XdE4329tWpR\/StsJvWNChQXLx4kYHDggULCCHSZfoQVBBcB4Yj1fLly2nAgAGUmZlJ7dq1oz179hhy4Y6+ZXZ2Rm5iKmdnHr29uWlNgwIF7nn06NGD3e3QK2KOBdmlsE6dOhEuhiHMuir2UMBNTGXPDK+2andeDxjd4YEsM+Sjf\/673hzFwDZW6eCmNQ0KFFYnq553hgJuYiq7Zmx3Xg\/uUn7TTTfRzz\/\/nM\/DmAOFngtAuOftpjVVQBHu1Y1Qe25iKrtIgDnamdeDR+5C7Ar0IwuWE8gFINzzdtOaKqAI9+pGqL1IMJWX8nrw6FRVqlTxhckTVWbtTVIAAsL4a8PvKaCIEIOrbsNDAaeBwmt5PcSguTJfEy399OipgCI8\/KpaiRAFzALFnB\/foiO\/HaG0omn0ZI1epkbtxbwewQLhaiUIvVgVwYyZekZQUwT+o7LZNQ2lD6PPREz1+PXXX1kwGxyliqncMHAcyX7xxRfsHsmuXbvYbVOcnHTo0CHfRTAEOn3zzTdp8eLFbM533303u9EKnw+xhBqizwgh7QgTCJ8UXAZC+EEwOWj07LPPUq1atahQoUK+YfF5IS3ezp07DTlcLdu3lMZueNHXxpPVnzIFFl7M6yHmEtXmK9UDPq3UpIyZRnZDGOvgUhkWC4ipFQNzc3NZ3EJ87rnnHhZOD7oiguEgvJ4YTo87giEuBm6ywsfjs88+YzroG2+8wXw5eAk1RF+wadsRJhCBfl544QX6+uuv2bxw5PbRRx8RbojqzQsX8zZt2mQIKMZseJGW77sW5\/GB9Jb0twbXgCPYnL2Y1yMQUJiZj1I9gnFHGH7HmxdvfrganzhxgrWoBQow+xNPPMH8LmB15ldoIVkgrBmC+cJNHFLHyJEj2eYQb7ACaHAPhUfkwnV17jEaSoi+YNO2I0ygrE2Ax\/Dhw+nMmTOMfvCOFeeF27Sgj5G7HlqJYkSDF6l5estgU\/X97sW8HnpAESgcnkxdUUBhmE1Cr8j1u8qVKzNPTmxoLWNPnDiRgYk2ShYA4KWXXqIjR46w5\/AmB6C0aNGCXSYTxXHE8kQOB\/SHkHxW76DozdiOMIG4lYvwgmBSzFMMRvzVV1+xC0BQR+rWres3r6JFi7LAx0aAAvOBjWLzsU1UN6W+KbUDz3oxrwcPrQ\/pTFQ9gtkBzOT9CH1nyJ8MNrZw9xeoPUdtFFAnEDcTYf7h6Slj7GBX1f\/9738zEEEAG73EMJyRobbAXmHlVmuoixFqmEAABTxhGzdu7MfQGAdULcwZth3k1hTnBQnNDFCEOi\/+nNfyeoBmAGAtUASTEGTHpsrhyir3mHheDy0hUeC+CMBAvMbOJQrYH\/AbVBJIFLhPAjAQC5coHnnkEV\/il1BD9JmYkl\/VUMMEYp4AA8wJQCcWDoBt2rRhhk0x9CA8CZ0ECozLS3k9ZEChZ8TUrrkIigiaBBfvQCVcAXsLrEQhElePCNhgPXv2ZGrGo48+6lMpuI0CbQAocEMV3nN4k0JS4aHzuI0Cxk9+VGUlRF8oQGElTCD605OUAr3d0KfTQBEKbWTAN2HChHwnX6G2bTavR6j9OPGcAgpNKHLxeBS6OYyUS5YsYXYM5AGBXeIf\/\/gHi4Eh+uCDkAADqDMQ1wEeCOCLmBb4wPOOp5LTExeDxfI0yxBWwwQqoDBLcf\/64ksBv\/C8HtojeGu9OPO0AooAQIElAFjgWjuyN8FoyX0I4CcA6UE0dCIT2euvv04rV66kEiVKsNMSAAROTHgOyWB2j0DRwc2wRDjCBHLVQ3bPQKZ6cABUEsXVlcJa6+X1MLOWbqirgCIIUOgtEuwX8IcIlv5t7969zH7xl7\/8JZ\/RD77+vBjJN2KEYcIZJtArxkwjdFF1rFFAAUUAoABxoC7AsahmzZo+SnOfAeQ6RSCdy5cvM6ckvD2gqoghy7kfguwY0WwGMyNLHc4wgU4djxqZl6oTWQoooAgAFFy8hst27969mTETb+xFixaxKMQwUiL\/Kb5DLMEvv\/ySZs6c6Yu8xV3DEcYcRjL4IdgZos+OMIEAuueee87PCzOQwxVsMw8++KBhh6vIsr\/q3SgFFFAEAAro6JAmcH8D1v\/bb7+dvv32W2aXAHCI3powbMJTE5GzkAsEBZnT4cMgemvieztC9NkVJhDABicy5OrAvOCgpufCzecF6Wvz5s2GHa6MMquqFzkKKKAIYqPABsTGwGbHiUeNGjXYken9999PRYoU8Vs5gAUMWDBmxsfH0wMPPMBsE1BRxGJHiD47wwRqL4UhhwlSJSDnidVLYZFjfdWzGQoooDBDLVXXEAXcxFSGBqwqBaWAm9bUURfuoJRRFUKmgJuYKuRJqAf9KOCmNVVAESXM6SamihKSRnwablpTBRQRZ4fwDMBNTBWeGclbUeH67aSuftv\/DyxnxdPZ82p8AAAAAElFTkSuQmCC","height":160,"width":266}}
%---
%[output:9e6249f2]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQoAAACgCAYAAAD0OyL1AAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQt0ldWV3kAMiQkRQsIzYFB52g4Cg1AbgUVZlYHhNRRBsSBgQQZ5lIpBFAusUqFYLEoFLAwPXTwU2kFGKAvtwAjULLSQQkGBGSyE1+INQRINYdZ36Ln898\/\/fuSee7P\/tbLEe885\/z7f2ee7++yzzz7Vbt26dYv4YQQYAUbAAoFqTBSsH4wAI2CHABOFHUL8PSPACBATBSsBI8AI2CLARGELUXwV+PTTT+nJJ5+MEnrp0qXUvXv3qM8OHDhAI0aMoAsXLkQ+HzNmDOXn59PFixdp1KhRVFhYSH369KE5c+ZQamqqIyDg8jp8+DCtXLmSPvnkEzp58qSo16xZM+rWrRs99dRTlJubS9WqVXPUHhdSAwEmCjXGITApjIji2WefpSlTpkRNznXr1tGLL74Y9V6\/RFFcXEyvv\/46vfPOO1RWVmbYp6SkJHrmmWdo\/PjxjsknMHC4Ic8IMFF4hk7NikZEkZeXRwsXLqSMjAwh9DfffEOvvPIKvffee4ERBYjhzTffFH9OHhAF\/kAc\/KiPABOF+mPkSkItUTRs2JAyMzPpzJkztHz5cvrOd74j2jp9+rRYWpw7d06Qx7Fjx8TnfiyKnTt30tixY+n69eti8o8cOZJ+\/OMfU6NGjUTbZ8+epTVr1tCiRYuEtZGWlib+DRLjR30EmCjUHyNXEmqJom3btvTd736X3n33XXr11Vdp8ODBoi2UGTZsGHXu3FksR+BL8EMUegsF73r88ccr+CHgv4AVI5c8KDNr1ixKTk521UcuXPkIMFFUPuahvlFPFEOGDBETUzspsQyZP3++8BVcvXo1sgTxalHAMkFb+\/fvF8S0ZMkSatCggWE\/Yd3gPbIsHK3Z2dmhYsKN+0eAicI\/hkq1oCcKODHxd88999CyZcuoVq1aNHnyZProo48Ik3TPnj1iYvuxKI4cOSJ2UE6dOhVFSEbAlJaW0rRp0+gPf\/iDWJZgSdS8eXOlMGRhKiLARJFgWqElit69ewtrYurUqWK5sWrVKqpbt66Y1FgGYJJiwgZJFNIqsYJ17ty54p1MFPGjfEwU8TNWjiTVEgViIGbPnk3z5s0TW5awJNq0aSOWCXInBA7FIIliwIAB9Mtf\/pJq1qxpKO+NGzcEcW3atImJwtGIqlGIiUKNcQhMCj1RIFhq48aNwtxH0FW9evVo7dq1Ykdi+vTptHjxYuGv8LP0kLsoX3zxha2PQl+WfRSBDX2oDTFRhApv5TduRBSYwNjlgH+gRo0a4r\/41YejUxt45dWZifZAOuvXrxcdNouR0Mda8K5H5euH1zcyUXhFTtF6WqKQE1\/7Kw6xEcMAf0W7du0CIQq0uWPHDvrJT34iYiQQRzFo0CCxxLn33nsFUpDh7bffFrEUHEehqPJYiMVEEX9jZimxEVEgCAp+gQ8\/\/FDUbdWqldgBQUCWnUVhB4+Mz+DITDuk4vt7Jor4Hr8K0hsRBQrJ2An8G7sh8F3AsgB5YKlg5qOwg0cbyAVHJUK44Xfgsx52yMXX90wU8TVettKaEcWf\/vQnsRTAg92P5557TvzbqLz29KjdC7VEgbLy9CgiMLdv3x4JD5enR+GXaNGiBZ8etQNWse+ZKBQbEBaHEVARASYKFUeFZWIEFEOAiUKxAWFxGAEVEWCiUHFUWCZGQDEEmCgUGxAWhxFQEQEmChVHhWViBBRDgIlCsQFhcRgBFRFgoqiEUbnx5W+o\/OuTVP3uxpTaclIlvJFfwQgEiwATRbB4Vmit9MR6ur53SuTz1JYTmSxCxpybDx4BZYiiqKiINmzYQAMHDqScnJzge+qhRYQhX7t2TWSF8potGiQBspBPzSY\/orR28zxIU7FKEPIFIohJIyyff3RVwVAZopChxKtXrxZJX1V4SkpKxKlHHJ5KSUnxJJLeogBJgCyCeIKQLwg5zNpg+fyjqwqGTBQWYxnUIMFH8e35Arorq1Ogy46g5POvzsYtsHz+kVUFQyaKSiAK\/+rCEzEMDFWZhFZ9U0VGx0SBnAZIdoKlAe6TRJLW\/v370+jRo6PSreP04K5du8RxY2R4RjlcBoMMSzjWbPYk6tIjDAWXbaqiRLz0CG+UVRljR0SBexsmTpxISMuOC3Dbt28v7mXARbRNmzalN954gxo3bizQwo1R48aNo06dOolUawcPHhTlunTpQjNnzqT09HRDVJko3CubKkrEROF+7JzWUGWMHREFMiYjOSuSn3Tt2jXSR9yIjavpcEM18htcvnxZkASsCOQpkKSANGkTJkwQn\/Xq1YuJwqmW2JRTRYm8EgV2uvAXqwc7CufPn6fatWt7dlaHLXtWVpbIc+rHoR6EjLZEUV5eLiwGLCeQ1h13WcpHplirXr26yJhUWFgolhjI7IyMz\/LBLddIloK6ZlfIsUXhfjjjmShAELiYqKCgwH3Hq1ANWOb4kUV+U687b0HAZUsUVi\/BdXSwJOrUqSOIAmnhYXXob3+C3wLWBHwWyNWoJRvZPhOF++GMZ6KQ4407R+Sy1T0CiV0DJLpgwQJxnULPnj3jlyiwpABR4JcBlgRugNq9e7fImai\/TxLfbd68mVasWEFIi6Z\/mCjcK31lEMXcrcfo+MUSapqZQvmPVRw3K6mt5FNxvN2PQLg1JEZxTRSHDx8W19yDELA0wcUyIAN0zshqQLZn7ISY3TUpQYHTtF+\/fqaX3IY7NNGtQ9FxqW79+vUpNTW1Ml\/t6F1hy\/f+3nM08f2jEVl+9oMc+tkPmjiSDYWs5Pvss8\/Ej4tKAXaOO1ZJBeWcgDX+wx\/+0FAHETHsNWrYTTc8LT1OnDghfA5Iwoor6ZAsFU8QRIF2hg8fLpQo1g8utsGOD8jQ7Iq8WMoYtnwzPjpPmw4VR7rYp3U6zeiR5bjLVvLBnwUdYqIwh1MSBe6PffTRRw11EI5YLP3DflwTxaFDh8RW6c2bN4U\/onXr1hEZrYjC6dIDa9YOHTooYVEg\/fzZs2eFLLF0JJkpQdjy6S2KBYMeoEHtsh3rpJV8mASIr2GisCeKX\/3qV9SjRw9DHVTOooBDEjES+BXIzc0VDpYmTaLNUNxpyc5Mx\/PId8HK8lHsPHqZ8h6oHTc+CnlbuhHA+usFUPbUqVPCGe9neWn1TsghCdHN++LSR6ENpMIN2XpnJcBAx3h71Pf8d9xAZRCFY2EMCsbKmWlm2SJgcMSIEeLCo8GDB0eWy0ERhZl\/TgtNQhOFdFy2adNGxEGYrYkuXbokAq4QHKKNwuSAKz\/TzbwuE4UxNnZLYO2EdjNxrUbR6p1VgijgkAI54HLZPn36iGWH\/sESpG\/fvpScnExbtmyhSZMmUV5enrio9ujRoxzCHQ5PiF0Fv8fgQxJNNBuUReF2i9aOKLQWhCQK5BzB8gBPo0aNKuzO6ZcWmAva5YpborB7n7TQcWQiLrZH4fXHVXQ422H2aEEzOhQ2dOhQ4bjKyMgwbUPFffV4nohhEoDTtoMgijV7ztC4NYcir8x\/LNfWV2K19HjhhRcIzsHmzZtHlh6IONbqsH4rX98enLS49BmEkp+fH2nH6dLD7n2ys3Hpo3CqHF7LMVG4Ry6eiczpeIMkQBbyeaJjA\/rtE3d22oxQs3Is6q0FI1LREgHCp\/Wk4Pad2p0du\/dJ4ok7i8K9+nqr4VRxvLXurVY8T0RvPQ62VhgWBUgCZGH1mFkU8vJl1JVBgWY+ClgV2MVDuW3bthFiGYyWJFIOt0sP\/S6L9n3yiANbFAajzEThfpLGM5G5GW\/4KNxs0doF\/mHSa7crgbz2lxz\/r5+4cscE\/g08TiwTM8vDyfvYojCZD24Ux\/2U8lYjnieitx4HWysIi8KLRFZEodczJxaF\/hCjtEzgv5NHEtii8DJSHuowUbgHLZ6JLMzxtrMotGeO3PgMtCNkRDhOnZn6ckbOUbYo2KJwzwgmNZgojIGxC7jCDod2twK7EGPGjIn6DMmaYC3g6gg4M\/FYbYe6sSis3id3Y5gomCiYKP4RxYsYgTDOengJ4dbGNbRt2zbqBLT8xQd5yEdfxg1RwM9h9T75DnZmsjMzELJgiyIQGJVthImCiSIQ5WSiCARGZRthomCiCEQ5mSgCgVHZRpgomCgCUU4mikBgVLYRJgomikCUk4kiEBiVbYSJgokiEOVkoggERmUbYaJgoghEOZkoAoFR2UaYKJgoAlFOJgpjGL3EUVRWKjw3A89EwUThRl9MyzJRmBOFUTi1Cqnw3Aw8EwUThRt9YaJwiZZdhqtYpsJz0xUmCiYKN\/rCROESLTuiiGUqPHRFf2R9+vTptG\/fvqiMWSjHRMFE4VL1jYtXlaXHpfdmUtm5rygpO5fqPP5zW+xUToUnJ7+8NkB7jkR7MI2JwmSYwzx2bKtZJgXieSJ67XOQ9YLIR3Ft+wo699sREbHqDJphSxaqpsKTpIDOaJ2n0sLQnmplomCiCGwuxjOROf1hAEmALORTq9vTlD1uuSWGqqbCkwlvhgwZErlXBB2Rn3fu3Dkq0xYvPXjpEQhZVAWi0FsUIAmQhdVjl7gmVqnwjHZd0A9OXONiOjj9hXHRpO+i8TwRfXc+gAaCWHpADPgobhzcTqltutkuO1Be1VR4bFEEoFRMFO5BjGciC3O87SyKWKXCYx+Fex2vUCNMxfEqXjxPRK99DrJeUBaFW5lUToVnteuhv0CZfRQK+ChufPkbKv\/6JFW\/uzGltpxkqIvxThRO+uh2EropH0uiQF5Ko8fsNvPKSoUHmcziKB555JEoJycTRYyJovTEerq+d0pEitSWEw3Jwi1RuN3vt3TIbT1Gxy+WUNPMFNMr9Kzkc9pHNxPfqKwVGcWKKPz2qbLrS98FEv5i50M+TBQxJgqQBCaSfGo2+RGltZtXQSo3ROFlv99MIZ3et2kln9M++pkUdmTERBGNrpmPwmypxEQRY6LQKzhIAmShf9wQhZf9frNJ6vS+TTcWhVkf\/RCFHRmZyQcrZM9fT9LT+RtCycLtp09h15XWQ2FhYeRV+pvR2aKwGAW\/zky3Zj+U9dvzBXRXVqdAfBRe9vudWhRm921aEQWski\/\/PIc61f6CCi63Iiyv8h9rFug8sCNcI\/lknc+PEo1ZWK3KEYWbAWCLImCLIkizXyuaG4sC9dzu99v5KOzu27SSz6lV4kZxzXwUZoRrJN\/VXU9Q2YVPiYnCHnkmioCJIkiz3w9R2A99sCXsLAqQhXyc3AIerHREevm0FggThT3aCU0Ut27dol27dhGCWvbs2UN169alkSNH0rBhwygtLc0UHT9LD79mv5nn3q1FYda5sLYp7eR76Bd\/juyc7Hv5e\/aaGXAJyHdx3xxKq3GZkjNyxXa0dCIzUdiDndBEsXPnTho3bhx16tSJcPjl4MGDtHLlSurSpQvNnDmT0tPTDRHSEkXL41tdHSv2Y\/abee7l5L5SmkaZD02llJQUy5E1IwO79q3iOOxUyYoo5m49RnO3fhVpIv+xXFMfBcrabcXayWL0vb7v1e\/OofKvi0TRMImCU+F5GS3rOtVuwQQI6Ll06ZIgCVgRCGyRpLBjxw6aMGGC+KxXr16WRPH26MeozbdrqXp6EpUXlwlHo5McBGjUyy+3kec+KauTozgL2RGrbUKj9kEON75cEMFBH8dh1w\/5fXlyfTp\/90Bq2LBhFJHJ7dUJuRupccp5OlmSReeyxxCWH\/rH6VasnkzsZET5h67Npo7JH0deiZ0l9B0+DRDFqJkFoTgz7SIzx48fHwlsQlltIhuvU8Hp3aNu2k9YiwIdwxJj8eLF1L179wgmxcXFNHnyZMrMzKRZs2ZRcnJyBbwkKIuG16WO7c5Hvi87XZvuqtM\/iiyMlPRI4btU9+\/TTSef\/oVS8aHIfe\/6deRrbCOWnS+wjbPQyqA1qdGQNi5DTxRJdTsLZ5722fPND2hfrZfEL76edFC+xt05ojgmGX6VtcFiJfWGV7B4+vx2L3W4uZJAFPIBUbT83u1bubWPmdNT9k+Wffe\/T1DJXy+L\/+3xb7nUkLZGmvm8xnD6S9LwiMUCbEGE2HHpVPvLSDnIflfd2wFFIIqwtkftMlxxKjw3dHW7bKAWxdq1a2nhwoXiunjt9e0wWmBNwGexbNkyQRj6RxLF715oTO0a3TZP8Xzzv8X09e7zJBOW6CcSJuUHZT+j2oeeNlTK8utlVFJ4OZIdSW6j\/r70hrBa8AxssOsfEzGHltxcRfWLP4wij2+OpFONlH+mU0n16Y+5Y2hSm8+jJqvWpL49oXPobFpv8UuOX9XIU95AfFde8llU9wsut6Sh+\/IJy4Ox6dMqEInVsJYnNxDv+q+vhwpiAGmt3nM60idZF++t2WSgsNAiFsn1Mtr+yS36PI0ilgciQQc\/eDOKKGUbZWdLhJWXfH\/08rGoJIt+f+b7BMJ4rf1W+p\/9Byq8X9+HsJceRsl1IYPegpD\/rw3hbtSoUQUd1i9n9LEPTi2KdevW0e7duyO3mRu9S2Il58Tz\/\/owFXZ8ie6rlx4hY6nHKOs085d7erhTI1CiAFgAYenSpZSdnR0lF77bvHkzrVixgpo1q7ifD+V9750F9PA\/5VD9lDtEAcWEggKM1Ae70bcXPo2sc+ULMNG0v1x6QKSCpzzYjUr+tl18rVd2bVtQfEkekqzk95\/XbEtNsq5QTsodq8fpAKAvIK6k+hX9HegDJrPWVHfaLsrZYSDb0ls0AlsDedy82+k4aNt0QxR2yxy9rCqnwgNRIB+GPu2dEd6SKH7d4jSVZbWkU0kNKO\/+2tQ0MzUqmQ\/qOsn85WVMZZ3AicKMyQGQ9nivVmi9leCnQ1w3PhBwShR2YeJGvVU1FR5ktZoH+r5oiaJtrRLLgXWS+cuPZihBFPp1vJ8Ocd34QMApUdiFiZsRhdEPlgyfRh25BDZzZmJCYymNctu2bRNWgNUywc3SQ7ZrtATX9qdKEoXV0oOJIj4md5BSOiUKuzBxN0Qhf9FjlQpPvr\/KE4VXZyYvPYKcgvHRllOiQG+cnMvR9lrVVHhhEoWTXKJ+NCPQpYef7dEjH98++t2o4oaI6B+ccGE8Jy6V0PGLNyJNw1GER\/9Zkzq3HZDlN4qinKmQ6+8niunm10XCwSl3UtzKKvunb1+2A0dn9dQcEYiGv+ppSVHvQn2juhEHanaucAjj0ff5+\/fXMewbHJ3ySWnTzbB9bT9lH2T7jVMuRDl9IYvExw1RuMVS1VR4noniX+6lTtF7A1GQOM0l6hZHbflAiUIGXCEASBuF6SbgauvvBkZ2PRA\/gAlodcLTT+dRVx9whKCknUcvic\/l80THBlHBSvpfONnG6Csr6futvqKujzS0nVR6ubVBV2i\/9MSGKELSHhOPHD7r+iB9W15EtXK6RU7A6pdx2FrGFrP+FwexDkaHzrR9Kym8UiGprdEyEQShHyO0f237Svr3uxdQUoMUKjtTQn9t+HPq+ijFPOBKe3+GdHxqdyHw2aZNm8QWaU5ODk2dejv+RHsXh56MwvRRzJ8\/n3r27GkbHex3LljVD5Qo8KItW7bQpEmTKC8vjwYNGkRHjx51HcKtzfITZudl2\/pJY0QeIAurB21gAje8eYaG9fqeCBDTTyp9vIUVWeA7O5PbLIQbx8v\/78h2On6+NmX+5Qvq0XeA4+hWJ3hLuVDWisSBY+HiqdShdB99XvMhavvsHOqfsVbEeoQdcKVqKjytk9SpMzMhicLoUNjQoUPFwbCMjAxTPfRzKMyJcrstg4m\/\/Yvz1LZ+Es3o39qWzY2Our99z7BIhGK3Tj2jDkUZyWOWactMdruzHnbH1N1i4rY85Jvxn4fobxeI8h6oHRWoFubSw62cqpZP2BBuP4CrRhToi93pTG1\/9UfdC5v0o1E0MVIEUZf6iE49Xm6zULmRz8\/Y2FlSZgfK9PJpLSwmCvsRYaIwwCjeiUJvUfzHAz+nt0q7Rnoq\/RxG\/gdteLW9+twpEWuisDtQxvko3IxmxbJMFAlIFOiSNsPVH5uNIbPEMV6CiIxULtZEYZdFy0g+6d8I8\/Sov+mpTm0migQlCn23zHYWvAQRqUgUdk5fJ1m4J06cKHKX8FMRgZMnT9KUKVMoIZ2ZXgc83pcebvttt6PhpL1YWxSQ0YwM7Xw8RUVFYhIUFBQ46WqVLQMSRS6Xdu3a2TrUwwQp8O1Rr8JWNaLwipO2ngpEYdUPO\/lAFvgL60FODAShIVANcSr6B\/JdvnyZsrKyKCnpdsoB1R7IVqNGjQrJiSpbTiYKC8TtFL2yB8tI0U+fPm2qRH5T3PmtH0v8nJw6jaV8TnVHFRmZKAIkCrd5E5wqi1k5N1m4rXJmGrVvt6PhRPZYKrkTh3Es5XOCn93yzWkbQZRjogiIKJz8ggUxYE6XHnY7Enay+K0fayV34jBmorDTgjvfM1EERBROfsGcD4uzkm4sCrf3etjtaDiRMNYT0c5hHGv54gFDKSMTRUBE4eQXzIliuCljp+hWOxJO3uO3vp18TmQIs4zq8sXaKtNiz0QREFGgGbtfsKCVXnVFZ\/n8j7gqGDJRBEgU\/tXCXQuqKJGZ1Cyfu\/E0Kq0KhkwUTBT+tdmkBVWUPF6JjJceBiPHAVfu5ytPRPeYaWuojh8TBROFPw3\/R23VFZ3l8z\/MqmDISw9eevjXZl56JDyGTBRMFAmv5Oyj8D\/ETBRMFP61iC2KhMeQiYKJwlLJ\/ZxfUWV9zRaFfx5jomCiMEXA7\/kVJgr\/E1QVDJkomChMEfB7fkUVJWeLwj9hMVEwUTi2KOIxS7jVFFGdyDiOguMo\/FO8y+sEvL7Qz\/kV1Sei6vIxUTBReJ23UfVUV3SWz\/8wq4IhLz146eFfm3l7NOExZKJgokh4JWdnpv8hZqJgovCvRWxRJDyGTBRMFAmv5GxR+B9iJgomCv9axBZFwmPIRMFEkfBKzhaF\/yFmomCi8K9FbFEkPIZMFEwUCa\/kbFH4H2LHRHH9+nVatWoVrV69mnDLct26dal\/\/\/40evRoys7Ojkhy69Yt2rVrF7355pu0Z88eUW7kyJE0bNgwSktLM5WYU+G5H0xVgnHidSKqjh9wVUVGR0Rx7tw5wvX0R44coSeffJLat29P+\/fvp5UrV1LTpk3pjTfeoMaNGwt92blzJ40bN05cZT9kyBA6ePCgKNelSxeaOXMmpaenG+oVEwUThXsE\/NVQZRJa9UIVGR0RxaZNm2jatGm0cOFC6tq1a6RfBw4coFGjRtFTTz1Fzz33nLgZGiQBK+LVV1+NkMKOHTvE1e34rFevXkwU\/vQ7UlsVJWKLIqABNWhGlTG2JYry8nJhMWA5sWTJEsrMzIx0B8uRqVOnUvXq1WnOnDlUWFgolhiLFy+m7t27R8oVFxfT5MmTRd1Zs2ZRcnJyBUjYonCvbKooEROF+7FzWkOVMbYlCqsOXb16VVgSderUEUSxceNGYXUsX76cmjdvHuW3gDUBn8WyZcuiyEYWYqJwqjp3yqmiREwU7sfOaQ1VxtgXUWBJAaKYMmWKsCTmzp1Lu3fvpqVLl0Y5OAEKvtu8eTOtWLGCmjVrFhcWxVdffSX8K8OHD6fc3FynY1tp5Vg+f1Crjh96p4qMnoni8OHDNHbsWEEIWJrUq1dPkAEsAyOrYd26dWInRG9t6C0KOE3hCFXhwe4OSFAlmbS4sHz+tER1\/NA7OxlzcnIIf2E\/nojixIkTwudw8eJFWrRoEbVo0ULI6YcoioqKxKQsKCgIu8\/cPiOQMAjgRwx\/YT8RopA+Au0Lx4wZQ\/n5+VEyHDp0SAh28+ZN4Y9o3bp15HsrorBbeqARkAX++GEEGAFnCFS6RWFHFAikQowELAms1+fPn09NmjSJ6s3atWs9OzOdwcKlGAFGIBYIOF56aAOpZs+eXcFZCeFBNl63R2PReX4nI8AIOEPAEVFIx2WbNm1EHAS2Q42eS5cuiYCrhg0bRkVhOgm4ciYul2IEGIFYIGBLFKWlpYIc1qxZQ3369DHcJsQSpG\/fviKQasuWLTRp0iTKy8ujQYMG0dGjRx2FcMei8\/xORoARcIaALVHgnMczzzwjznaYPSAQBFylpqaS0aGwoUOHioNhGRkZzqTiUowAI6AUArZEoZS0LAwjwAjEBAEmipjAzi9lBOILASaK+BovlpYRiAkCTBQxgZ0o7ERA2M5GQNxbb71V4RCe1+RCMYLK8LVh4Hf8+HF67bXXCLt0165do1atWgnfGnxwNWvWjMiRCPi5HUsmCreIBVA+zERAUOKPP\/6YXnzxRXEGwOjcjdfkQgF0PZAmwsAPMUCIRL7vvvto8ODB4uwSDjF+8MEHgix++tOfRsgi3vHzMghMFF5Q81knrERACH\/H2Zv333+fysrKqG3bthWIQsa6eEku5LPbgVUPGj+EAEyfPp2OHTsmDi42aNBAyArSRfpHWBk49dyhQwdKBPy8DAQThRfUfNQJKxEQDugh2xiSBw0YMIDwHhxR1lsU8R49GwZ+V65coeeff54QUKg\/24RgwxEjRogMbbA04h0\/r6rLROEVuRDq+UkEBKKYN28e9evXjx5++GHxb6Mj\/4l8HscPfmbDuXfvXnEs4eWXXxZEkcj4Wak0E0UIE95rk0EmAjI7yesnuZDXflVWvSDxg8xYvuHw44YNG8TSAyelExk\/JorK0lQf7wk6EZAVUXhNLuSje6FXDRo\/+CfWr19PL730kkjQNH78eEpKSvKVcyV0EEJ8AVsUIYLrtOkwEgFVJaIIGj+QxLZt20QipR49ekQdcPSTnMmpPqhYjokixqMSViIgL0ThJLlQjOGq8Pqg8YOzFEsN+CR69+5NM2bMiDqj5Dc5k2r4OZWHicIpUgGXCzsRkJlCJ4ozLgz8sE2K7WUEqT399NMidgIHHbVPouDnVp2ZKNwiFlD5sBMBmRFFomzvBY0fHJeIoQBRYMmBLVH4JPRPouDnVo2ZKNwiFkD5ykgEZEYUiZBcKGj8pOPylVdeERGtuDbTiCQw9ImjHWfiAAAAiElEQVSAnxcVZqLwgpqPOpWVCMhqLR3PyYXCwA+BabhsGySABExG9+MiEVPHjh3FyMczfl5Vl4nCK3Ie61VWIiAroojn5EJh4IdlDIKqrB7cdIeAKzzxjJ9HtSUmCq\/IcT1GoAohwERRhQabu8oIeEWAicIrclyPEahCCDBRVKHB5q4yAl4R+H\/kn6lMdaPv4AAAAABJRU5ErkJggg==","height":160,"width":266}}
%---
%[output:8a90d031]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQoAAACgCAYAAAD0OyL1AAAAAXNSR0IArs4c6QAAHKNJREFUeF7tXQtsVcW6\/sGmpVJraXmUp6DAAcQQKFqS0wOEeDVCyuOYCociyEMeFhBRUkBBWw8CAWsEFMqFy0MvL0GCBBpfBGIhcsm9poKA0BMMUCjyqlCEGh4335i1z9qr69m9ZnV2+VfSiHvNmvnnm\/m\/9c0\/s2bq3bt37x7xxQgwAoyADQL1mCi4fzACjIATAkwUTgjxfUaAESAmCu4EjAAj4IgAE4UjRNGV4Pvvv6fhw4eHGb1q1Srq169f2G9Hjhyh0aNH0+XLl0O\/T5gwgXJzc+nKlSs0duxYKikpoczMTFqwYAHFx8e7AgIhrxMnTtC6devou+++o7KyMvFcu3btqG\/fvjRixAhq27Yt1atXz1V+nEgNBJgo1GgH36wwI4qJEyfSjBkzwpxz8+bNNGvWrLByIyWKyspK+uCDD+iTTz6h27dvm9YpJiaGxo0bR1OmTHFNPr6BwxnVGAEmihpDp+aDZkSRkZFBy5Yto8TERGH0H3\/8QXPnzqUtW7b4RhQghqVLl4o\/NxeIAn8gDr7UR4CJQv028mShniiaN29OycnJVF5eTmvWrKGuXbuKvM6fPy+GFhcvXhTkcerUKfF7JIqiuLiYJk2aRDdu3BDOP2bMGHrxxRepRYsWIu8LFy7Qxo0bafny5UJtNGzYUPwbJMaX+ggwUajfRp4s1BNFt27d6IknnqBPP\/2U5s+fT0OHDhV5Ic3IkSOpV69eYjiCWEIkRGFUKCjrhRdeqBaHQPwCKkYb8iBNfn4+xcbGeqojJw4eASaK4DGXWqKRKIYNGyYcU++UGIYUFBSIWMG1a9dCQ5CaKgooE+R1+PBhQUyFhYWUmppqWk+oG5SjpUWgtUmTJlIx4cwjR4CJInIMlcrBSBQIYuLv4YcfptWrV9NDDz1E06dPp2+++YbgpIcOHRKOHYmiOHnypJhBOXfuXBghmQFTVVVFs2fPpu3bt4thCYZEHTp0UApDNqY6AkwUdaxX6IliwIABQk3MnDlTDDfWr19PKSkpwqkxDICTwmH9JApNldjBunDhQlEmE0X0dD4miuhpK1eW6okCayDmzZtHixYtElOWUBJdunQRwwRtJgQBRT+JYsiQIfTee+9RXFycqb03b94UxLVz504mClctqkYiJgo12sE3K4xEgcVSO3bsEHIfi66aNm1KmzZtEjMSc+bMoRUrVoh4RSRDD20W5fjx444xCmNajlH41vRSM2KikApv8JmbEQUcGLMciA888MAD4r946yPQqV94VdNgJvID6WzdulVU2GqNhHGtBc96BN8\/aloiE0VNkVP0OT1RaI6vf4vDbKxhQLyie\/fuvhAF8ty3bx+9\/PLLYo0E1lFkZWWJIc4jjzwikIINK1euFGspeB2Fop3HxiwmiuhrM1uLzYgCi6AQF9i1a5d4tlOnTmIGBAuynBSFEzza+gxememEVHTfZ6KI7varZr0ZUSCRtnYC\/8ZsCGIXUBYgDwwVrGIUTvDoF3IhUIkl3Ig78LceTshF130miuhqL0drrYhiz549YiiAC7MfkydPFv82S6\/\/etSpQD1RIK329ShWYO7duze0PFz7ehRxiY4dO\/LXo07AKnafiUKxBmFzGAEVEWCiULFV2CZGQDEEmCgUaxA2hxFQEQEmChVbhW1iBBRDgIlCsQZhcxgBFRFgolCxVdgmRkAxBJgoFGsQNocRUBEBJgoVW4VtYgRsEFj45Sk6feUWtUluQLnPtgsEKyaKQGDmQhgBfxDYeKiccjYeC2WW+2zbQMjCd6LAyrz9+\/eLpbzYPQkbpWCjVXy9iCXDVtfZs2dp27Zt9Pzzz1OrVq38QdUiFywvvn79utjtqbZ3gXZrS1BvEbf2SG0gXeZsTzjSIAmQhXb948lU+ugfnaU3h+9Egd2Yc3JyKD09XXzGfPToUXEYTO\/evSkvL48SEhJMK6UtJd6wYYPY9FWmY9y6dUt8zYiPoho0aCAdZLsC3NgS5FvEjT1BAlaX7alJHzf2BZAEyEL25StRXL16VZAEVAS+AdBIAZ8gT506VfzWv39\/0zoBtA\/+878pPb0XtWrVshprYjyGcZn+shqjaQ2AtMY02r2H6t2idwZ3diQKs7z0v5mVYaygXXo7R9Ce2\/+virC6o07oHPrxqb7ToXxtDKv\/N9Jb2aL\/HQf5gEChtswwNivLqaPa5eP0LOzp2CKJ3hxQfW9NY32s8nIazzs5rRk+jzZNsJT9xn5j1iZmysCsPprt2j39c1E59NC2gceuSfoj7NDQ+BAJZ0yYbc9uZEmnjqO\/\/9fHkkRn1l96IPG7Po3+ntmzxrLN8oLjGi+7vIx54FlNMloRhRtMrOplhR\/SW9lu9ruWj7FuZvVx02aR5uOmre3sMHMqM+fzUo5Vu9cEI6v20fqwWRuZvTTctIXXNL4qCmyxhs+ZjTsrI24BNYGYBfZBAGHoL+O4y2slojU9Ou6rfZqbDoPuV0xktqVxPA+SWPjlLzKLDCxv2crCV6LA7soHDhwQ+xEYz2rAvd27d9PatWvFgbVMFH+qiveHtGOiCMid9EThRrEFZJYvxcgOavpOFBh+mKkG7KSEmRCzcxzqWqO5bXkEooY8kWRKFCpjYieR3dY9knT1f79E9X\/\/9ynsbvOCM2mBPygJu+GW2zxlp7v7YArdfbCxYzH3BVEABZWkdovEGDp3zfw0bscWc5kg\/ZEE2j7hCUKMAqdnNWvWrNrp3pDGn\/1wSZotXuuZ1rIBZbRPEgE82FZcWkHnr98Osw9ptMt4T\/tdS2N13w5CkMR\/\/PY5HTx40CXS0Z0Ms4dfP\/x3R7KQPfsRmKKwG3qgKaFE\/p6\/iW437kQD0h6hZs1S6X\/L\/pzl0Dq0vmNrnU1Lo3UH\/f9ndk4Qndhth0Se+EM573xzKayHaXlZOYHeQfT2wh4z0nnn6caEPLGDNY7kw1DN6iyMwoMVISxQvlXdnfDRP4d\/T0hPIi1v7d7\/nK4U9v76e7iDjX8qSaQ3u\/R5GNNo9zTbtHKRz85jlWE4azbo69H0QaIzV6tC2Fw89RNd2\/2eOKukZcuW0c0CDtaDDD\/88EP6y9C5dCbusWqEDJzQTn07NZa+6MpXoqhpMFMjiuHDh4sOkJaWZnl2pVPP+HDfefGmw5sPgULt0n7H\/2MRD\/4wBahNu9qlN94z5ml136xs\/KZPj30mcdI3zuqs7TUdsE2zZ9u\/HqDvf6kUVXCqn1Ob2N23ai\/tGSM+eKFgAZ+23iaSslV\/VltbhB3Te\/bsSWZYYQo7iEWDvhJFTadH9UQRRAdQaRGPSragHVS3x7gwT3Vnj8Q+lerqK1FoC66w4lG\/CtPNgqsgQVHJGVSyhYkiErf2\/9kgfcLJel+JAoUVFRXRtGnTxNmWOASmtLS0Rku4nQyP5L5KzqmSLUwUfw69tLNR9X3MSunqz0XR0tsd1Kzf4dxJPddpojD7KCw7O1uMKxMTEy39O0hQVHJOlWy534lCc2J8a5Sbmxvqq1rfNBIAAvTG5QBaHnjYbpkAgtdt2rQR56vEx8eb+kWQPuH04vVdUTgVaHU\/SFBUck6VbLnfiQLqAAF5KwfX3zt58iSNHj2aFi9eLD5i1F9W9zS10qJFC2rbtq3luiItryB9wslvmSii4OtRp0b0877qxCXTeewWBRox1sgAp6wNHTrUVRPobcchSGPHjhVfWFs9L7OurgzWJWKiYKII6zN1iSicvgi1cv5z586Jb5PsCEAfy8jMzLQdQmjlGIcqZkMXvU1MFCZ0FiQoKjmDSrbUpaFHTffw0B+xqHXTbt26mQ5HrI5eNAtSamn1CsKpzzvd96oKIknPioIVRZ1UFH7sBGWc0bCbzdCGIlAjuIzkYjas0cjDKqjJRMGKIoQAKwr795wRH7fO4+dOUPphhtOUJmqjEYxGLFZTrlrNEdw0+1jSbV0jUQpun2VFwYqiTioKVEr7cE37kM3OKcyGBmYzGVrwEvEFqAez6U2NGPA87peUlBA+TzAjGbugKBMFKwpWFC5fZTVVFC6zD0tmF1w0TnnaObF+GhTrMezy1dKePn26WhyEiYKJgonCpScHSRRWi6WsFmKBAAoLC6spBT0xoJqYBjUu4nIzu8FEwUTBRKEgUWgmaQSgN9EqNmE2S6IPerpZm2EV1GSiYKJgolCYKFyaJj0ZEwUTBROFSzcLcujh0qTAkjFRMFEwUbh0NyYK89kSl\/D5loynR3l6NKwzqb6uQ6W3rG9eaJGRSnVlomCiiEqiePXVV8WxlXX5KisroxkzZiix7R8TBRNFVBEFDrOG89xPu3BjH1nZB3c7ES4TBRNFVBEFjAVZ4K82LgzNKioqqHHjxoFsaguCqG2SAM5MFEwUUUcUtUEQWpmqxXCCwoKJgomCicKDtzFReADr119\/palTp4pNdI3bgJntmYn9MkeOHEkNGza0LCXICK9Kja2SLWgctsfeEVTDx4PbRpTUs6I4c+aM2HgUjm22rLW4uJhycnJERBqbdBw9epR34bZpItU6HtvDRGGGgGuiuHHjBm3fvl0ccXb58p8HxBqJQjvXIyUlRWwllpCQINLxuR7WnY8dM7ocU7X2ikgmeHjYNVFoH8q0b9+eBg4cSAUFBdWIAioDQ4wVK1ZQv379QmZUVlbS9OnTKTk5mfLz8yk2NraaiTz0aK7EkYKqOQLb48GbJSZ1TRRLly4Vjj548GA6fPiw6UYcfpw96mYHoUjxUKnzqWQLxyice5Zq7eVssT8pXBOFvjirtz9Ux4EDB2jVqlXidG795eY0c6tdgPyp6r9zUamxVbKFicK5p6nWXs4W+5PCd6Iwnpykmen0XT4PPXjoYdalVXNM1ezxhwacc1GOKLCGf9CgQZSamupsfQ1ToLHLy8upWbNmlse51TBrz4+pZIumKFTBhu1x7k4xMTGBrBANEYXTTj1uhx5WisLt0APljBo1SgRFZV1VVVV08eJFMTyKi4uTVYyrfFWyBQazPfbNpho+SUlJ1KhRI1d9LZJEvhKFH8FMfACTlpYmVVFgQ9MLFy6IMhrU8spMlWxBR2J77N1JNXwCVxRe2MYqnsDTo15Q\/DOtamNetse+DVXDx3uPq9kTvsYotAVXzZs3p7y8PF5w5aJNVOt4bA8ThRkCvhIFCigqKhLfgGRkZFBWVhaVlpbyEm6bvseOGV2OqVp7uXgX+ZLEd6Iw+ygsOzub8GFYYmKipdE8PcrTo2adQzXHVM0eX1jARSY1IgoX+XpOwkTBRMFE4dltAnuAiaKWZz1Ue0OxPdE1FAqKKZgomCjC+hoTBROFb8FMGSzGQw8eevDQQ4Zn+ZMnKwpWFKwoPPiSaorLg+kRJWWiYKJgovDgQkwUHsCSkZSHHjz04KGHDM\/yJ09WFKwoWFF48CVWFB7AkpGUFQUrClYUMjzLnzxZUbCiYEXhwZdYUXgAS0ZSVhSsKFhRyPAsf\/JkRcGKghWFB19iReEBLBlJWVGwomBFIcOz\/MmTFQUrClYUHnyJFYUHsGQkZUXBioIVhQzP8idPVhSsKFhRePAlVhQewJKRlBUFKwpWFDI8y588WVGwomBF4cGXWFF4AEtGUlYUrChYUcjwLH\/ydK0obty4QevXrxcnmJeVlVFKSoo4sHj8+PFh54ya7ZmJ\/TJxoE\/Dhg0trWaiYKJgovDHqWXk4ooocKoWjvo7efKkOMW8R48e4kTzdevWUZs2bWjJkiXUsmVLYV9xcTHl5ORQeno6DRs2jI4ePcq7cNu0nGpSlu2xdzPV8JFBCmZ5uiKKnTt30uzZs2nZsmXUp0+fUD5HjhyhsWPH0ogRI2jy5MlUUVEhSAJqY\/78+Xyuh4tWVK3jsT1MFDUiirt37wrFsH\/\/fiosLKTk5ORQPhiOzJw5k+rXr08LFiygkpISMcRYsWIF9evXL5SusrKSpk+fLp7Nz8+n2NjYarbw0IOHHjz0cPFmqaUkrhSFlW3Xrl0TSgKHpIIoduzYIVTHmjVrqEOHDqHHELeAwjh06BCtXr06jGy0REwUTBRMFLXEAi6KjYgo9u3bJ4hixowZQkngxPIDBw7QqlWrwgKcsMPtaeYIlvbq1cuF6TVPopK8VskWIMr28NCjRkMPK9hOnDhBkyZNEoSAoUnTpk0FGUAZmKmGzZs309KlS6upDVYUt+j8+fOE81pr+2R1Jgrnl49qROpssT8paqQozpw5I2IOV65coeXLl1PHjh2FNX4QBWZXBg0aRKmpqf7U0CQXNHZ5eTk1a9aM4uPjpZXjJmOVbNGIQhVs2B7nHhQTE0P4k32FiEKLEegLnDBhAuXm5obZcOzYMTFVeufOHRGP6Ny5c+i+HVG4HXogs1GjRomhjKyrqqqKMOULNRQXFyerGFf5qmQLDGZ77JtNNXySkpJEjFD25ZooEJDEGgkoibZt21JBQQG1bt06zL5NmzZFHMxctGgRpaWlSVUUN2\/epAsXLogyalvuq2QLGpPtsXc51fAJXFE4MZJ+IdW8efOqBSvxPFQJT486IRl+X7UxL9tj336q4eOtt9U8tasYhRa47NKli1gHYSV1rl69KhZcITCXl5fHC65ctItqHY\/tYaIwQ8CRKDAmAzls3LiRMjMzxbDDeGEIMnDgQLGQqqioiKZNm0YZGRmUlZVFpaWlvITbpu+xY0aXY6rWXi7eRb4kcSQKBP3GjRsnvu2wukAgWHCFGQSzj8Kys7MJH4YlJiZa5sELrnh61KxzqOaYqtnjCwu4yMSRKFzk4UsSJgomCiYKX1xJSiZMFLxxTVjHUu2NyfZI8XvPmTJRMFEwUXhwG9WIy4PpESVlomCiYKLw4EJMFB7AkpGUYxQco+AYhQzP8idPVhSsKFhRePAlVhQewJKRlBUFKwpWFDI8y588WVGwomBF4cGXWFF4AEtGUlYUrChYUcjwLH\/yZEXBioIVhQdfYkXhASwZSVlRsKJgRSHDs\/zJkxUFKwpWFB58iRWFB7BkJGVFwYqCFYUMz\/InT1YUrChYUXjwJVYUHsCSkZQVBSsKVhQyPMufPFlRsKJgReHBl1hReABLRlJWFKwoWFHI8Cx\/8mRFwYqCFYUHX2JF4QEsGUlZUbCiYEUhw7P8yZMVBSsKVhQefIkVhQewZCRlRcGKghWFDM\/yJ09WFKwoWFF48CVWFB7AkpFUtqJY+OUpOn3lFrVJbkBD\/1JPnDWCM07NzimRUT+rPH\/55RdlbIGNbI9966uGT1B9VTlFgQOQ09PTfa3\/d+djaP5310N59nn4ApWsmy0OW\/a7LK+Gl5WV0YwZM5SwBbazPfYtqBo+rVq1IvzJvpQhirNnzwqHOXjwoO91\/r3HGPqjzV9D+cae3k8P\/t9\/+V4OZ8gIBI0AXnb4k30pQxSoKMgCf35fRkUx628P0d+a3\/a7GM6PEQgcgftOUchGGDGK4tIKymifRLnPtpNdHOfPCNQpBJRSFHUKWa4MI1CHEGCiqEONyVVhBGQhwEQhC1nOlxGoQwgwUdShxuSqMAKyEGCikIUs58sI1CEEmCjqUGNyVRgBWQgwUchClvNlBOoQAlFFFDdu3KD169fThg0bxFLjlJQUGjx4MI0fP56aNGkSapZ79+7R\/v37aenSpXTo0CGRbsyYMTRy5Ehq2LChZfMVFxfTsmXL6OOPP6bk5ORq6U6fPk2LFi2ir776Stx75plnxGrSNm3aKNElZOCDOi9evJj27dtH169fp06dOgksMzMzKS4uLmLMgwKuNrEx1nHPnj00ceJE0Zd79eoVFAQRlRM1RHHx4kWxVPXkyZM0fPhw6tGjBx0+fFh8UAVHXbJkCbVs2VKAAYfPyckR33EMGzaMjh49KtL17t2b8vLyKCEhIQw0EMu3335Ls2bNEuvmV69eXY0o8NHalClTxH18TIZnUGb9+vVp5cqV9Nhjj0XUEJE+LAMf1HnChAn06KOP0tChQ6lp06a0e\/du+uKLLwRZvPbaayGy8Ip5pPX18nxtY6O3FR+V4cVWWloqXnhMFF5a0kXanTt30uzZs8Ubv0+fPqEnjhw5QmPHjqURI0bQ5MmTqaKiQpAEVMT8+fNDpIA34tSpU8Vv\/fv3Dz2PJePLly+nzz77jG7fvk3dunWrRhR4k77++ut07do1QQ5wGFzGsuvVq+eiJnKS+I1PVVUVzZkzh06dOiWUWWpqqjAcBIk3IVTG2rVrKS0tja5eveoJczkIWOdam9jorQKm+fn5gmihcJgofO4Jd+\/eFQ6K4URhYWHY2x6Az5w5U7zZFyxYQCUlJWKIsWLFCurXr1\/IksrKSpo+fbp4Fo0VGxtLV65cESSDZ4YMGUIoB4xvVBQ\/\/PADvfTSS8IGPUndvHkzVPY\/\/\/lP22GNz5CEZScDn99++43eeOMN6tKlC+Xm5oaVd+LECRo9erQgXigNKA+3mMvEwSzv2sZGb9Pnn39OH330kcAKypaJIsDegLc8lESjRo0EUezYsUOojjVr1lCHDh1CluBNCDWBmIVGBCAKxBwGDRpETz31lPg3Or2RKDZt2iTeovi9efPmAdYu8qIiwceqdBAnOvtbb70liAL4uMU88hr5l0MQ2GjWglwnTZoklCleVhg+M1H415aOOWFIAaJAUBGdd+HChXTgwAFatWpVWIATGeEextiQzO3aVf8wDPfNiAK\/\/\/TTT+ItgLcCHANv3L59+9K0adPEW1fVy098UEcMzwoKCmjbtm0Cx86dO0eEeW3iFgQ2qB\/U7Ntvv00NGjSguXPnEoiWiSLAltdYGjMeWuzAytlh1ubNm8V426g2NJPNnsW4ErGR48ePi\/E5FEVWVhZdvnxZkBE6AYKZ3bt3D7Dm7oryGx\/Uf+vWrfTmm2+KtyOCuzExMYIozAjWDebuauJ\/qqCwAWZbtmwRahT9BDuqyd7NzX+0iKJm1sNY+TNnzoiYA4YPCEZ27NhRJImk05o9q8UhEBDDlBbKhHPgKi8vF86CwCmCe8bZFBkN5jZPv\/FBh\/\/666+Fcnv66afDZo8iwdxtffxMFyQ2CHhj5ujdd98NxcyYKPxsTZu8jh07JqZK79y5I8bGkL92qkB\/z+vQA1IbjYznoES6du0aZhnKx4yJ1XAmIEjCivEbHwQEMdRATGLAgAH0zjvvUGJioi+YB41PkNhgNgj9tH379kKVai8YJgrJrY63Gubr8VaHhMNYuXXr1mGlRhJYs3ozvv\/++7R9+3bTIYvTcEYyJGHZy8AHQy8oNixCw8wP1k7Ex8f7hnlQ+NQGNhoh2NWxRYsWlkPhoLBxU05UDT30i3rmzZtXLViJCkcyVWdFFLt27RIOYlxJh86HmZK9e\/cqMSPiNz5QU4jpgCgw5MCUqPZW1HeuSDB300n9SFMb2Fy6dEnEtozXzz\/\/TOi\/iPU8\/vjjIhiuV2h+1NfvPKKGKLTgE0DFOghMh5pd2uIfBB31qzCtFlzp87AiCsQiMM6EhNTnicVa+B2LjiDLsTajti6\/8dECl4jSY8UqovRmJIH6RoJ5EHjVJjZm9eOhh6RW11a0bdy4UXxjYHYWB4YgAwcOFM5aVFQkpi0zMjLEDAWWy9ot4dbHMKyi91AVGGf27NkzbNYDC730wVRJENhmKwMfbakxSAC4mgVqge+TTz4pbKsp5rLxUgEbYx2ZKCS1Otbqjxs3TnzbYXWBQLDgCuNns4\/CsrOzxfcJdhLPLnqPPH\/88UcxDQsZiznx5557jl555ZVa\/yhMBj6oI9al2F1YwIYFV7hqirmkLhPKVgVsmChktzLnzwgwAkogEDUxCiXQYiMYgfsUASaK+7ThudqMgBcEmCi8oMVpGYH7FIH\/B7bBKxBUhCboAAAAAElFTkSuQmCC","height":160,"width":266}}
%---
%[output:60440beb]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 443\nStatistique T_max     : 23.8281\np-valeur (bootstrap)  : 0.0010\nPoint de rupture      : 1992-02-01 00:00:00  (indice 50)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 49\nStatistique T_max     : 4.4865\np-valeur (bootstrap)  : 0.5060\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 393\nStatistique T_max     : 50.6578\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2000-01-01 00:00:00  (indice 88)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 87\nStatistique T_max     : 4.3450\np-valeur (bootstrap)  : 0.6290\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 305\nStatistique T_max     : 38.5016\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2007-04-01 00:00:00  (indice 85)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 84\nStatistique T_max     : 7.3516\np-valeur (bootstrap)  : 0.2280\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 220\nStatistique T_max     : 12.9215\np-valeur (bootstrap)  : 0.0210\nPoint de rupture      : 2025-03-01 00:00:00  (indice 211)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 210\nStatistique T_max     : 4.8588\np-valeur (bootstrap)  : 0.6050\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:1fbff167]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:83cb2646]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expA_fit\nTaille de la serie    : 125\nStatistique T_max     : 10.1641\np-valeur (bootstrap)  : 0.0770\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:1604147e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl4FFW2\/pM0EJawb5EAQUFRccTlPVBQQQV3dAQERQUEVNwQWRRRARFlQEYQZ1CRgYwbPHd944CgoILCc3QYF0RECRB2TBCQdKRJHv\/tnM7tSlVXVXdnaVL1fX6G6ruce6ruX+eeNamoqKgI3uVxwOOAx4EoOJAUbwD5ZvtBXDNnLRrXroZ37+yIZmk1oiDL61LVOHAkdxt2zuiH5qMWIaVhi6q2\/IRdb1wAREAj77fDihGPX90Wt1\/QMmGZ4hFeMRzY+ef+qNHyFDTo+0jFEODN6poDUQOIDho9Tm6Ivw08FY+8sxGXdWiCi09u5JoQr0PV5gAlkG3jz0Vg75ZSjPA1boUWUz7zJJNK+IpEBSC7DhTgnoXr8XT\/9mFHlFGvrfcApBI+ZI8kjwNlxYGoAITEEER6\/WUtftx9KHRs+XH3bx6AlNWTOsbHtdKBeLqRyv3gowYQ47Ke\/XgrHnxnY+i2pwep3A++slAX6egiNNY4sQvSx7+P5Fp1KwvZHh3FHIgbgOgcJZjM\/3y7Z4XxXjPHHPAkDcesqlQNywRAKtUKPWI8DngcKDMOeABSZqz1BnbLgbzXHsWhtR+o48qhr5di94w+SK5VD+njl6DGiZ3cDue1LwcOeABiw+TVq1fjhhtuUK2eeOIJ9OvXT\/29aNEijBs3Tv39yiuvoFGjRhg8eDDuvvvuUBvj0Hof\/nbVVVdh6tSpqFmzZjk86so9hX6EIaXbJ1+KpsNfUETvnjMUxz282DPjVsJH6AGICwDRN\/yf\/vQnPPfcc44BhO3fe+89zJ8\/H+3atUN+fj4eeOABbNmyBfPmzUPDhg0r4etRfiQRQAQ0Dufm4Nf3nlKSyOGc7z0AKb\/H4HqmqADEaMI1m7Vd01rHhBJVJJDzzjsP+\/fvV5ud15AhQ1C3bl18+umnthLIjz\/+aCqdWN13\/RSPkQ4HV78RdmwhkPAY03TU66jTufcxsspjaxlRAYjOAjqPtWtaO8x1fdn3v+Cf3+7BjL7tE55bAiA8xqxYsQJPPvmkWtPo0aMViEyePNkWQOTowqNO586dQzzJzc1VY\/De\/fffn\/C88hZQ9TgQE4BYeaRa3U9E9gqAPPzww1i7di2OO+44ZGZmYuHChRg2bBjuuusuRwAye\/bs0PFF+CDHGP67qutCPEeyRNwdQEwAwiVTApn\/+Y5QAJ38e\/A56ceUBEIFanZ2NtavX6+edPv27XHBBRcoBaudEtWTQKw3h+dIlpjAIVTHDCAciEeW6+Z+rcZsULsa3h7eEacdVyexOVNMvUggBJDWrVuHWWTk33YA4ulA7F8Fz5HMnkeVsUVcAKQyLixeNOkAcuaZZyplKC9aU3755RdHEgjbi9VG9CCi\/+BvVlaY7GxA\/tu8GWjdGsjMBLp1c7a6FQA+BsD\/y8WuAwFkOhvCa+VxICIHPACxeUF0AOnVq1eY6XXDhg2lAGT79u1hI952220hBalTPxCCRlYWMHGiOXEEkUGDgAkTzH\/PBkCY04HD2JIAsrySAUnhof3YMeVyFGxYFSLXi4Op3AgWM4Doxxd9qZHMuIyVWb4hV+UQqV3dV7k5VM7UETy6dw9KHrxE4uD\/eS1YEP7b\/PnhEglBo3sxzewy6ChQXKCtgRKJjksEEYcCjRpFTK38u8F1k0yT\/9i1YeIgXs3vWxiiTMAjpVFG2H22PfJLjhdMV87vodPpYgIQ8QeZcnVbx0mEJBHR2a3SPAAxPCWCRps2JcBhBAdpznY8Sa0oFjGWLw+CCDGnuLsCCQsBRQ1DCWVB8YBOQUR39mJXMw9RuzYCLrXO6RcGFJ4VxumWrVztYgYQs8RCVkv87fcA7l34A9qn18H\/Ze\/zAERjlC558HhC8LC7KKkQRCidbNoUlDyIKZQ6HHQPtacEQhCxu7j5c7PuU9nBklLT1HGj3lUjw5y8IrURkKjW8hRQ4tAlEM5NacP\/n8Wh2JeCDWuwY8olSD390lJt7Wj1fi8fDsQEICSRxxEmEnLiNMa2vNo2rYW5K3NKAciWXD8ChUVIr1e9fFZfiWaZMjkZUyYnoUXLIvz4U6FjynpclIxVnyah20RgxQSgxZEi\/FjovP8p1VKU5LL0SBG66P1+2YJqKUnwNSlRtxIcxMWcBBJAanXsGXaMidSGAFHn3L44vPU7FGxdZwoKDKjL+58S2cnqmBSJQTk5OeB\/iXZlZGSA\/yXSFROARHJpN+pAeHR5etlmzOx\/Ej7\/6ddSAELwuPPV73FKem3cel6QiYFAAIcOHUKtWrXg8znTlZRXn3jSl7M1Cd3OCgbUvfyWH526hANApDWtWZWMAX9MDSo2JgDX+gOYduB3x7wbm1Ydb6b6VL\/H8w6V8HvJn1Az+zM0vXN+CERiAZBqLU\/Fwc9eU6BBkLACkFg3D4FjzJgxWLNmTaxDlXv\/Tp06Yfr06QkFIjEBiBsOGzOWsS+TMS8adroaZuXGfej113\/jzdtPR6c29dQ9f74fO3ftQnrz5qiR6qw8RHn1iSd9PIZc3jNFHUXWbThSiq2R1tQjORmr2iYFzyzdSiQJp3x4OTkZw1KSMOBIEWb\/dijE76L\/W4S8Z29B+sTlqHlqUM0ayxGG4HHo80VhazPqQeJhhRGrGTdiixaJUx6CgDdr1izllKiHO7jZYxXRNmYAMQMGLiSSFYaWG+MRRgDk3TvOQNe29YMA4vdjx44dSE9PR2rq0a+sg6u8+sSTvkmTgiZb\/mdmmo20JipNs6kRpXNHt6Aug9vdjg8ErY9pkmkNTKTSJBsY\/0IABw8eQJ06aUBeDg6syEJat4Go1jR4jOEGP\/Xne3H1lNvUv6NRorKfmQQSLyuMAIjZRuRRTUzbm22UzA5etbg2iUR3XCeK82AxAQiPMAPnf4vpfU7CglU5oYTKZgF2Ot0egISDolMAubllS+T4fMoJTJSeSnF6FICUrbZbUBAhHtgBiMyp0IaD6fbfCC\/ZuCHrMPTXU1ULiZI1WlB0M65ZJK0ZgMTLCmO1EckiM7caO2tVnPeb5XBVFkDECvPWV7tDytRogumqsgRC3w6aZa2sLwIG5xcDCN9C3XLCvgsIIINKLDBOAWTi0cKm3ET8735N4iv6aTV2TOwedoRJSrKWkuKx0eJhhTHbiDRXB\/2Hg\/whABMvRRoR0JU1SOiB7hSoOwQa10qvYgZaXnjhhfjPf\/6jdHb0Ut68ebNKLvXRRx+hY8eOqpu0M45RJQGEZtlbsr5D9xMboku7Bhiy4FvMG9QBq37Mc51U2QmAOHHttts4Zi96NH3ieYQRABFzrJFG0kfAWKMd47gJNhU3VOZc\/l0slrDYsd2alASyAMjcFPQfYdfOFQwgcryJxQpj3IiRfGN0oUsvEE0AWbJkiYq05sWo6SlTpmDgwIEqGZTZ5heweOaZZ3DJJZeE2rHv008\/rSK36bks7TwAKeaALm1QCpHSDv8z7A9hzmUCNku\/z1U9jdG6kQCkoCAdixalOnLttts4lRFASBMdyAiQ4hSm08k1tfP51PFFLgEQcT4j+BAMuCkINv+00R8pAOFgE0qkGZ13VhKIfoRx44lqjLqNxjzrRMoxAoiAhJWvi\/jO6FJIJACR+CfSQqmEwMCcLpQ8mPKB+WFOP\/10DB8+HHl5eYpkpr7s2bOn0ud98803porSKimBOHmg0kb3FxHz7+BzjgslIrICkC++2IMbbkhHTk5w89i5ds+Z48dJJ5W94jWeEogCVB5DFgQ9SgkiRgDhBtiVmqqkBcUHSg2aR6qKjZlf4onaNRDA37duNVVAc4xJ2UfnI+gUO50ZFa9WADKy81\/xRNZZigY3StQ9C0aF6t6Kg1jj4fPinmnMuBFF92Gl6zD73ewIw2jsHj16YNKkSZgwYYJKQcnYJl6MyjZKIHKEYfyUJ4G4QQqHbY2KVjMA2bgxgHbtSoBDXLslVoRgws1mdO1+5ZUd+OMfG5Sp5SbeAKJLITqI8AvJDd\/M70c1nw8rRQo5CgCZWUd\/m1jiicoxdLE8IxBQxpmLNcmFhhfqBBQQUeJhdG\/xM3Migdx79lN4csl1UXmiyqshFhejF2s8zbhihRH9h5V3rkgguju\/UQIRuq10I8wL4wGIw40fj2a69UbyhgiANPjwHBTuLzHZEhwG3FSEufNKnKtOOTFFgYbRa1O8MjMyAvj6uwKk1nRm+i3wF2DHzp249ZZWStLhuPTClIA245r5+9IPC8F+l\/ZIxq5dqaoP79ldMlfzZs1K0Sf+ICJpLV1WhB4nJKnN3vz3w\/AlJwWPMVq4Leclb\/QQf\/48jP4hKUe1nlYXG3UHDml+J5dXS8HGQCDktNfsuxVKNyJ+IFSi3nfe3zB9cR81qltPVCGFEohReilLM65wwRjzI+Ci65NIoxWAUFmqSyCyHko9HoDYvflx+l10IcO6ZoTpSAgglCwCOUEfELmaH3dYbc6d26uF7m3bGnT55vXxl\/kYcE0N5GxNRkbL4AZe81kK7hlzWP1nduntX367AJs3FeKb\/+Rh6sTmYfNYLZlzcy7OSVp46feM\/UgjL5m3SdN8cN4aNWqE0c57Ic9SGYSfTnp26wk8KHdrkXB0Xw+F4AKKtr8t2o9vft2HT1o2x5c1S0ID1lRLCY5cLLu3CBQFAak7sO2TpLB5MvZkY8Pu7LgCiB5sp9d6KUszrm7ClQhlLllYaDzeWAEI2aaX+eC\/9URTkrWO2fdHjhypvIH1FBDUi9x+++1h5UF0EJLsdlXKkcwNrphJHtJfAsPcjMdNy2MMN7GZtHDo99JenSKlyDwcQ0DAzdxu2grY6VKN3NPnlvVw7FUEpmLvUkdziXKEjSmd\/VyE9MV+7HgiNQwUtolUIrtKzA\/sb8wylA10G1yik6EEEssRxkzy0NdG\/5ADH81TwXopDVtAFK9pFw4xTRtgxhcrZaRuypV+AiaRopYd8T4OjaqcElXC8q\/5Q2MVSCf\/zvvtsKkXqp1viFgh4vAsQkNIDg19TKtjSTznjXks3WQgwCAbPF6pxIwAYkF0ZjawqXhOAki0StTA3hzkvjYJzUYuDCuS7SQnqq9xqxCo2PHWiScqvVCPOuAq3U+82GlHl93vVQ5AdCWo8VhiFqEryZZ1Rj5+dduQFaYsAMTuoVXK3\/UDueZh6ppWM4lCH8QpgGj+JgSQaDxRa\/2hR6lMYySlLOq9JOpGTFS6o3JlN0oTlD7GvP4DsgZ3QLO0GrCTNsw2wwmn+9UxpFXDVNPjCKWJSNKDSBvxkjDs5nO9oZ12EE0fZW6CiZt0YU7nYDszADGADnUgDOcXh7Wy9kR1Q75VW7ONqMfAmPWrDJJIlQYQY4rCaABEN+Pe1KN+KbCwAwi7340vjl175ZhVDFpuQMk4rg5EtqAkRxeaaJnW0JC0lHiiaNHkbnUvmp1nFRxiM5ZVwF80JJRFHztXdrM5ja7sbKOXLmUtIClJyvu8pBCYWGbGjh2Ll156STmWSZlSKmOzsrJw0kknKd8R\/dLH5P0qBSBcsBxh\/nhmU\/T6y1roTmH8jZcxyZCeP1U\/vrCtDiAbV9ZXjlVyyaa02oBWG9yoA9HbWY2l92GWL1788tJMKikEjS8hzcY0\/0pflR1My2vK9pzbbs5ssbjQxXxSiWdpiA\/842jMTPZ8w9l9BZAdQVLRQSb0NzOZfQxkrihelwSJaACVsTtb+ZFINC6n\/+8mn+Ck5cEsq2XlTRoLsEQCEGMS6WLjUygAUeblGB9\/\/HEIJAQIxo8frzLof\/jhh+o3WkucAAj7sYA6286dOxf33HNPqYLqVQ5A9GRC4pYu97qeUK8UeOgWGD4o\/chjBBCG8+s6kfx8f0SnMCv9iTG2hADAja0Dk4AKQUCSFglwSDsR3SWhsYTey+9duwYBhJexL+\/pYGIEOx102rCMAyWf4pTqkqaQvhm8FH1MqjwouPE3FYNGm0lAtkNTQnHkfkhqIUYq+uYDgYzieXJ8CNQOyjW6Gdcu32ksG599dUcyZmNvOnwuds64Dk2HvwDd5BtpnkgAwrXqSlOJkzELpqPkIBtfn4+xLmeffTZeffXVkFRB35BIEogHILG+GcXFp8a\/s1EV3K5TI0UF4em+IEZPVHHt5tTRAIi+MYX8UAh7sUs877MdPV55\/fhjwBSoBEC4+QkiRknCjj6dfRxLBzAddMy8InVPVHqZXpaaqrxNdacoPWBMn4ueqHr8DH+LlC\/VzhPVSUKhaF8V3ZGswZUjQxabX\/8xE4fWfuA4K3s8ACT4LgQLojMiVz9uSLAcXdUppfDIIgAybdo00AdEv6666qpQ2VJPAon27SgGEEkixGEkivf2C1qqUQVAxp2XhvPSgxv6+us7K9Gfx4cnnlhtOftNN2Vg48aSXJI+Xw5eemmVilHgNW5ccBz96\/\/55yXj0eGHYinLVUoffbJzzumMoUODOTZfeKF0zsqFCxeZ9jMjmGPJxTGHDCnJ3Tmuc2cFDo\/l5OCi4pye13furKQFX04OXlq1Cs\/266faTF29OlSuYWNGBm7KyFBtjBfvBTIy1H+8+O9PLfKF6nw4q\/5hJL9wU6mMZHY5UaN9RXRHMt3kW+Q\/gJ0z+qH5qEXKN8TuiheA6PPoHqh0kZdoW+pDTjvtNBW560kgdk8mxt\/1JEJmACI5Ub\/8RxZS17+rZsvLG4F9+0YEX3xfDtLTb1D\/N16BQAb27JkGvz+4OevXn4UGDWaFmm3d+gnYRr\/atDne8Yo2bfpZjckxDh7sjdTU1aG5OIjbsWTiJk3GoE6dN0J07Jk+HQd790adN95AkzFj1P0dr7wS2vwtzz9f\/dvfubP6ne0Un0aMwL4RI1B\/1tF1zypZt75Af6dO2PHqq0hdvRrpN9xgu\/ab0\/Nwa4+OrnOi2g4coYHUgGnYdwL2\/eMpNB44AzumXIrq7c5xnJXdjQ5EPbviYEJKZnJJkBxzefDSQ\/J1ANGrC1IysVKiekeYWN6K4r4EkEhHGDYjiDApbvKhvarXvHkZpb74BJCLLw6CiEgVAhxytCj1ZbeRQOyWJxIIdSiUdBj5yqOMXLo042QsafPiizlo27YEEEWS4O+frzaXuB7OyMCyjAxcnJODyTk5oFPUjcUSxotHpZS2FtLFvKMZv1\/IyMDQnBwMcZCxvHmNgErua8zKblfWwW79dr\/r2czY1q2itiysMKRDgvOM+T44HyURD0DsnmyMv9spUc2GF51FURGg6y+MbQkcUurRzlchUt4NqyVyTAENUczqilvS5\/QSHYhV8iDRg1jlrzAmwRFrbKTaLmJtkMRB0bqWlLUS1SkPI7Xz\/EDiwUXnY0TlSKYPb0wUxN+Ybd2sbKVuxjUmHLIDEPk9lAxY5WEoXWzaDkAi5d2IBCDym6Qd1JW8TgFErzxnlb5QV4jaJcERk6ye08OUj0dvMmDMaQGpSK+PXb5T569e2bRMVHNootIdE4AIeGTUrxEy28o9vh5Oat\/yvJmdnR2yuesOPGed9S5ef70DuEH1+4x6lPOpHh3J+\/379wvLbq4XtGamqBNPfAXTpgVrsFCXQb1KpPE6deoXKjfJ41PLlueH2osUQiXv7bcvUpmneHEe+gvQoUjoo\/6kU6e\/qvWI9MGz9QMPPIDrr78+lMqf7a+dNUvpK3iZ1beduGkTPpYamADe3b8fV9WtGxYpyjV16tcPLxw+jCnVqinlKfUmIzp2DPHayDvhadls7cijRoqJiVcsTEWsy+mcVRJArDxOnXqiCihIwlo5T3Lz0fGme\/cVRwsEXeb0GYTaibckTXGPP\/44ZsyYoTYz53vrrT9ENWadOnvBbGedOuWHxpw9u6FlmsVIRJO+u+\/ODaXCk\/M1lXKjRo3Cgw8+iG3VquF6vx8727e3XX\/d3FwM9fmwbs4cDBgwQLWf8a9\/IbdXL2xJTkaDX3\/Fm\/Xq4Q+5wTn79++vsmvJXGyv88l2wnJsINXs6nTu7WhW2YgjRowACzUlyrVt2zZVEKvK1YUxK9FgV9aBD5WSAU2mtKXzomefURq57bZXkZt7mvp937596Nq1q\/p75cqVqF+\/PurVq49Vq1aiV6+rUa2aD8uWLUOHDqfhxhszwhLsyEvEl2v69C9Qt+5V6mu9fn2wSrzdlZ5egBNOeAzvvDNSAZtIDVTeknwCgpOrbt1c9O17CH37NsHatTPRp08fjB49OuTVaARQztNt0CB8fcYZoQzinIdSCT1E+f\/e+\/fjq7p1LaeXEhDiQEUQzczMVLzn3wLWRknIyXrKo41VnhCrub3KdOXxVErmiOkIQ0mj56yvsDnXj2R6EqIIyQjPghWpwBTJ0GMLzDYQnXiYMs74sst9egROnTpVrYib4Nxzzw0db4ys1DeP9JsyxXc043ZJsiK9Dzd8z57bMWrUQeV5GGmeSEre1q2LcOONARw4MDqMPjEDils01283j9maap96KtZ27IjFfj+OPz5ont77r3+h72+\/YZLmeqtLOHSEcjtX+b6awdl05zUnfiDs49XGLb8nFRWA6Lk\/SKoxrsUN+cbgJDnWECC6deuGtLQ09YU2u09gcboJdOnGbKPWr38Ntm1rh2+\/\/QZDhlyEVq2K8P77Y9WG59fa6TwTJizHZ59VU9JSIPAT1q17Hy+9NNQU4GIFELs16WBq1LdEA1ZunqvbtpF0IGUR9u+WPq+9OQccA4gOGmJleeSdjaFqdNEy2Agg+jgiMRiVe27FcOM4ZpIOFZkiERnFeqv7xtRzTueRfmYA4vRY4WYuXfKQuiZWPCjLdHpujyPRvlNev\/LjgCMAsVKK3vHKd9h78DDmD+qA2tVLapbQEkNwGXtpG5UfJNJlPMLIl55iKN2DGV+gi9v6\/UaNGtkqAjk+JRV9Y+gbSklQxYpW\/m2mWLS6L2HbAjBO55F+RgCxokufx+1cbG+WCNjpXPF8FQVAGg\/+M\/bOv8+xe3o8afDGii8HHAEIp9Sjb+XY8uPu33Byeh08sSQbktqQCtS3v96Lt4d3hGRcdwogsjmee+451UXXSOtmXP2+boo0arD134QGCW5iISAmsTXOYzVevOfhvEYA4b1I8xh\/d7ImnW\/SXszWdnPF91ULjkaryqHPg\/VUrC5q0XTfPDdm3LKg2RszwrMqKnLqBhU+CJMISRU6\/jKwczqyVu+AEwcx74FUTQ6IBNJkyCzsmTfCk0COgdfAsQRitVbqRi6Z9SVq+JKx+oH\/xrTFm1xJIMcAD70lOOSA3RGGWdn1urgyrCeBOGRwBTSLCUB4rCFgPHp126h1IBWwZm\/KCuSA3REmpX5z1Djhv1Gn282q7CXb12h5iuOyDhW4tCo5dUwAUiU55i06Zg6YWWP0e\/s\/nIeCretUCL9nuYmZ3WU6QBiASBamu+++u5Qzlij8aM0YdMe9uGfhejzdv72ysvAY8\/SyzZjZ\/yQliTh1ZbdamfgsMLuTMflsmXJDGzxRYxPKgj\/yPDg2nenojRvp2UkxaaP1SD\/C7JrRF4G9W0oNk1yrHprd\/RJ2zb0D6aNew6H\/LAkrNsUObugpC354Y5ZwwLEEUp4AIkBGi4lkvy7vh+YBSAnH3WxYWn3IOwkmdPPcjFKI6EOMjmRu6HEzv9fWPQeS2rRpUyRmPaMEoueFrH7CuTj0yw4EmrSHv30v1PyKqWiB\/DMHIylQgLQ1T+OGq7rh8QkPYfOefbh20D3I+36laqNHuxqjY40vmgAVzay8CCLXXnutyk\/JoDuaeCX4TjdRinmWfejSzoterDTtUophgBr9PZjjUvob2aXTJnPpgW4MRBO6rEzMxhyYTvro0bvG9eu8s\/pNnhOD5JiakTRarVH\/gm\/ZssVyo+tz0b+FV926dUMSiM4rWfO7774bikgmz+fMmXM0IdQL2L9\/v+pPOlke4auvvgq1i\/TcrPhiRo\/7V9\/rEQ8OJC1cuLBo9uzZ6sHy4kblEUYvCsxNvmHDBuU30e6CPujWezBen\/MEzm5dDzNnTMPcFT9h9oR7MfCai\/DQg+PwyOTH8dbij\/H6q1nI256t+nHD0fFLxmc0KDcXj0RGKcMogYg0oG8K\/UtHumUspswngHz55Zdha2rSpInaLEuXLlUvr9FnxDingBPb8UUWUKIIz40iPKOTm9n6hH+kzayP7jBHnoi0pa+LtMo85J3OL11CEr6eddZZpeYSz1P9ZZEveCQAIR08QurvhYwvPjTkzYknnhhGl06\/BB7KsyAtQvefT\/sVrZP2YfzGZjg9zY\/Bp6XhmZSe+Oqb70u9i\/REjkRPpCNVPDaJN4Y1B5J++eWXIr0QjmzwM888M7TZ+QD5Rbp50GDk1mmL2VMewK13jQoByLLvtuPe4beh35XdMPruO0Lt3nrucVQ7\/FsohFzGpBSgf1mN5FkBiFEa0MGHX0RuNvnqycblZrHaeLp3qvSXAkJmG1R0Q0Lfk08+qUg3c0izkubYR+bVJSijVKV\/fTmH8Uili\/FDhw4FK78LCFkdv4wSjPDdqGfSj6sEd+ORgToO\/ZjCdfCZEij13wRA5Fnw35NH3ILPP1yMZ0YPwvGDHlfAwL53nRTA1Hn\/g51FdZH13kdIrpkWem7MfK4\/Q+8IU3kgTR1h5GXl1zQSgPAhHqrfDl\/UvQj1v85CQaAQ+WcMRNKRw6j9+Sx1vCloeQ7qfD4TKfm5Yas0O3YYjzfSwQ5AzJS9ZQkgOkjoiyIIiqShp\/Mn0EXqw98oBZEnxs2hHxvZTp6NLo3wS85NdMe0aWqe+669VgGIAJyd\/sZOAokEII888ggeffTRUuULhE4WThJwMQJI4aEDuLd3N3y4xR\/2brAvc4pOGX0HCjZ+gVEtd6KgMElJJ2e3b4N7n3\/PA5DKgxlhlCTdc889RXz5+TKLroEvokgL8lXTX6rTLrkRw+4epQZK7ToMf7\/pRDw5\/h71ZTVuCKt1y0tuZmWxAxDjC8454gEg+tEmkgRitSbRC3BDPPTQQ2BSG6NFy\/j1NEpIMrZuieKz4bk16wE8AAAgAElEQVSfks6MN95Q+UE2Hjmijml1f\/lFPatXtm9Hj2rV8FqHDrZlEp0CSKtWrcJSGIgkYZRAdH6YHWGkX\/X8XEy46TJ8m9wSf8t6MVQCkr+Tpvvvuwf536\/CU6+8h99TGyjQOLNVA1x98EM8srsDWh\/f1pQe7whTceiStGHDhiIzqcNKB6JLEvKlEf1IpN\/0r7IcRaw09nYAQnbZ6UDkpXV6hIlWB0KFoOgpmMVcdCXypRY6dL3Jm2++GfpKi06GvBP9jfRh4KCuM7p+3Dh8MW2ayi7G8gxG3Ut\/6n6uvBL\/tW4dnrzyypiyW8VTB6IDyD\/HXo17PspVtOn6E6597D134NDaJRjdaif8R4ISSMeGRZjwwtuY+dYKS52MByAVCCA8wogmXH9hRe8hVgS2odJN9A66qG38Tf96Go8puvaev5mlcHMCIAIiEngn4GX3hY8k3uvBZdFaYXSJyqhzkLUaq55RwStfewE83QpF3cKamjXxfEEB8vv3L2UJMh7pztq\/Hztuuw1vjhgRFoXs5jXTn6FTKww3ss5D6pMIljqAsEjUJyfciPGPBY9f8u4dyduOe3tfiKLDfvz572\/hSPP2pscWSsueFcbNkyzbto79QIQMPbO6Tppd5rGyXcaxPToztbNcLstZOr1YHsJNe6fjxtIuUtIgGdeLe4mFw+Xf1xWASEj\/lKvb4uKTG5U\/tcfwjMFS1uHFn2W5rGHFnKZu6rmwfkxWcdU1K7bpxazZptY5\/UIV4PTyDXpxJz3gTXfwsrof7SPjeLwa9H0k2iG8fuXAAdcAoruwu6HPWD9m8DnpoVIQbsY5ltoSNLjJjeDBzU+wmFC8WEoT9NLRK8s74QPLNrIivdXFTarHnGwbfy7SLhyCuhcNwfbJl6Lp8BdU191zhuK4hxeDNWvl7\/wNqyF1cg\/nfG96P7mWdbJnO\/o9ALHjUOX43RWAkGTmAWEioRl97csN6EvU+4kkM\/ic4yDFtSsHO5xRwQJRTKLMSwp2szbMwIHBQldOLlVAKhvInBQEC+M43TN5aAlCSCYyQ9KHEVwizSUjkKTAniBM6aUqjX0l8rVay1NhVsLy8NbvcGjtB0gf\/z5Y9FpAhvEqZvdrnBh9WQWnAJKoCZTNnhuV8PwvkS5XAGLMSqYv1K0OxFj6gXVxX\/1iB67\/r3S0apgKVoo\/cOCAckfPyfHZbli9vc9Xkl7R6mG4aS9tf\/klDcziTpAwggUBICsrWDN37twAzjorSLsZLQSBwSuAQR8DAwYE0KhRSVsFThyHdRuyiS4XYCK6haQR0vLtgQN48ihfWAtGP9YwATuBTK4VFwAtlgXQ+OAB+LZ8icM\/rEJat4Go1jSIcqzyt7xYUaKXrTycmxOSLthux5TLUatjT9XHKLE0HPhnEFjM7jut5SL0Go9UvF\/jxC4KsMykmUQu4WD2XrKOzfTp0xMKRFwBSLyQUa+TK2kPV27ch15\/\/TdmXncS+p7VDAX+Avz7q71HY1\/SkdkmGQOPFrTWv+5qwy4AFr2SgrnzinD2f+Vjx86daN6sGVJrptqSyvGdtmfb006pjtZtUtD1vELbsVd+moxVnyYh99ffStGyLTkZ3bKBpYWFaj1GOlYlJ4M6j+zsQVg1dDBaPDwf2y7MwqHDR4IbWKP7qtq1VMnKLoVBmqZMTsb4h0voG5acjGcOBflSf8MS5C+4A+kTl6PmqUGUoRQ1YQIgG7feVSNVDg7qP+R4Ul4AIjSkNMoI6WHkHmkwAxGx+HDTtWjRwva5VOYGa9aswaxZs2IyvVfE+lwDiDGVoRDtVAIRXciwrhlhilgBkGFdW+DW8zJwwVk10aRpPs49P8n0K64za82qZKz5LAUrvtoAmkRr1IicyFltxIIC7Nmzx1H7zZsKMf2xQsz46xHHYw+4pgZGjitCl\/NL6uTkbE1Ct9k18fJTfnQ6HNzoOh17aqXi\/rTqeHlfQWh5HGfNyym4tnkA0w78XoruAfVrhNo\/Pb0a7hlzONR3bFp1TN57ILjOvPVInXN1KQAZ2WUu9sy5Fawp3KDfJKW01GuxJKWmKQmE4BLPI4yupNWfpa9hC7R4Yg2kBkykfCB2HrcVsaGinTNR1+IKQHTJYcGqnFBJByeV6MhY9r\/xb9\/Al5KMURe1NgWQN28\/Hc2S6mHSw4cxYfJ2pDdvjhqp9oBwec8UDLxlO665pr6j9v58P3bu2uVo\/MEDC9G48UE8NqWG47Effuj3o9Xn6mPph5pEMCQZ2x5NwvstgtIEL52OK2stQjZWIBPZWFr4UahNrcdTlDqEUoiRbpFY5hYWhkkg1Hi8nJyMUb8dUuus\/8t65M+4IgxAHh6xHYP3tEP6fYuQXLdJSBHKictSiSrmXB5\/9GOOLvnIkSVSRrJE3XRmIJOoa3ENIGKFeeur3SFlqpMEQmxz20vrsCXPj+y9\/lLJl0UCefeOMzBvSn00a+Y\/6l6\/A+np6UhNtT+SUBz\/5z\/96lzvpL3f78eOHfbjc7zFi\/3w+\/MxdWpNx2M\/8EA+2revifbtU5VeguNksYrnhKBFJQQgxXTc3PJmrPQNUR4fVJpu0uwn1KsMHhS0qDQ3oVusLXIkIXhQx8t5ZJ31c39A3tTLwgDkvvP\/hruacc6Sq\/a5\/dBs5EIlheye0Uf9YGauZeKf9PFLIIpSMeMa75ttFquMZLQCmSUZkjGMPiKVadPpgZHRJMGqTGtxI0W5AhA5fnQ\/sSG6tGuAIQu+xbxBHbDqxzzM\/3w73r2zo2UdmBGL1uHFNbvCaNMr2gmAjDr1DARy6iM\/3z2A1K69Bx07NsDFF9srUZ0CCDdlp05+LF7sHkAIOMOHp4KZEjhO6wuAzd1KzLNKAikGhPNb3owcCwBhu6RJwKAJwBwTABFrS9YkYOAEhOroDrIBkAeHrsOdzW5V+gVeoiwtD98L3YTs5oXV21aWTUc6WONZ0lLQMzgrKwvjx4+PmL2tMq7F7bNwBSByDNGlECnt4LScg5UzmgBIkx9Oxsx7m+DDpb9jyLAdYUrRy6ulKJ+JFkeKlBJSv6hAHHDTJkyZ3AZz59krOp0qUYcNScYzf8nHQ+MLMPmx6o4VtNL+7rtqKXo4jvp\/cjJ43JBL6OjZZiByfBcfvT0RLY60xo+FP4Wt75RVKfhrlyO4o1oKJm\/bgSvS6oZomZKcjAGFhZj2eAoOP1SEjMLCkNVGxq+96jkE\/jEtTAKpKACJ5JHqxhM1XgBilkKTzDcWHLvvvvvw\/vvvq\/QDdKnXo9idAoYe3sBYpubNm+Ouu+6yDYB0u7HLq71rAImVMDsAqbf6dCx7J1UpLa8ftD1MyXlBw5rYlpKkAOTj3PwwUqhAvPb6TXjq8dZK2Wl3OVWijr2nOiZPP6DoGfNQsmMlqrS\/pV9dvPx2ATjOtKd\/h670JI0D6l+KzUU\/IyXFh20ptKlQniBItkZGYSu8vG+x+veAfTXwcv0C\/BT4Ha\/iCL5vXF\/xIeNIEd5I9aETrTTZwD1ND6t7IYAqKFBz7M3\/P7TYTff2EisMAeT2apehxZTPoCtL3Zpf7XhdVr\/bAQgd8PjBodE6kls\/47N4SfyXVPJjegvG8tCdYODAgcq8ymBJlkFlTJjezxjfJHllCDjdu3cPSSJ6pURjbWMJjizL8qLxfhblCiCMo+kz9z9IRpIqRPXUdSXOaCKB\/PrTOUh\/swYyHg5g\/uCtYUpOkUDIhHXFZk1hiEggt97SCh98VGL5sGIYlZE9fcnYUcNcv8KXjlJOj4uS8d7\/HsIpU6sjZbJPvYzqfnJy6OU0SkMcm0pUKl3\/2KsW3v\/giBqHClVKCzS7irvGKdXaIjvki7ocmRgc+jd1IesOb1S\/3rEqBe93CVei7q0VTG7cI+VjHD4yCJg0AZhQ7OGG1vT0UD6t21I+Ue0y9gAbdpcAyIT792Po711CHqc7Z1yHwtwtqiqcrsvQ\/TN0CcHqvi5hWPlxxKMGjB2AUDckABLJI9dYvU\/XYeggQUlFTxhN0FiyZImSIPSLUoaAEKOtaQVr0KABLrvsMjBfClNeMImX3t9uLfHe+PEar9wARCw4D1zWBvcu+gE1fEl4984zQjoTfi0+yfWjsHdq8HMxCcgYHzA14Zp9UahjuPnmTbjlljYh5yhhktmXiLqHdj4fciI4nXGe7MFA1\/EBrHzZF1JeqPvaExDnU6GLY1OJ+uVdacgp7sdxNk0I+rJIoBv\/v0L9i6PxP0og3PT8L3jRBzU7ewIys0qUJ3Qma+vzKTZ1x2Blt1F9J008SiPHKOltpLQbNoW+xuTZqEtKlKXJjVoh\/b7\/UYpRbvADH81T0oleZkG3iuh6DP2+\/F3vintDJmBdqtGd1vL+9ynUObdvVDVg7DadcJbciAQgOkjoG0t0G7Vr18Yll1xiKYFIHykCrwPN4sWLFXD8\/PPPKCwsxEcffRSq1exJIC5gjNLH+Hc2Yt7AUzB4wXdIq5GCBy87PmTKla+FkuDpwc0gEQkGMZlH9xhXG3oSEBgQQNuXfVhe3E\/fmvqWUnv1qFJSgYAgAUUCen8aXdG5l2k94RWBnrDxi+lRQoBEwa0AMoutMPJV1EEoCCu8Q9sJ\/w7Ch3J0p1tqZhaUB1po+3dXv2YrooPHHkyaD0wo\/tvi2RCQBELEamPWtGDDGmXWbTbqdeyZM0R5ooqPCJ3Mmt67CLtn9it1v8nwedg5pSfERGumLI1XDRg7AHH6ehrTT1C3wQxvTI\/JPLq8WHTdqAPRE2gbpRhJ3UBdyYUXXoht27aBnrMdO3ZUaQokr6ynA3H4lAggc1fm4G8DT1U9bsn6DrTmSCxMCEBcbtjQ9JTcizdsZvH5IHyDFu\/PSUCLlkWodksSsnWwkMg2uoFyD4f5iBui25ysmfRwLfrh2zLGXv9WkiiCAgeYEMSUwZnAcgEVIYw\/yN8kml2MEkhpQnUTcSQAkY3f9Pbnw6QIcTJrPv4DBSy69yrjZ5qPeRO7Zg9SxyKRZCRORnw75OhDUKp1+iXYMaOvZQ2YSKyOF4A4eZxsYzzCOOn37LPPKv0HdW7XXHONShjNfCZGHUp5r8UJ7U7aOD7CsHjUNXPW4po\/NFaBdPLvvN8Ow4kXqh2AnJDrx88Ni\/UR3Ct6OKqTlZhtWOnH8eTc4SQmnhuf0oa+V90k5JAEHvoYpMU0sYcE6\/PHzSUoyMkpVk1qDQyk9FFyrAkuhoOzvXZk6b4c6CbikjXTJhaLUnosjN7azJW9FFDEACCcy2kNmEQHECP9urJVT5x9zAOI7m1qdEd3EqErR5iFt3bA6Nc2YN+hAMZdWnKEoRKV+zzQtn6pjWbUOZR6qaw2rDQcXHwymR+uEdDHCRNGBIx0AKBeloDiBIBID\/e72cGb9zm+kkyCWpDgJRQUy03ZBJD5xccWHVzYlsChImZKsYISRolSNvxnKlEZjStOanqdWnEWo+Th\/+lLHN68Vh1Dav2hB3ZMvgSBvO048ssWJKc1hq9JGzQb\/RZ2\/7lPqftN734ROydfhMLDBSg8sBcpjVqh+glnI33MG04+A67aJOqmM1tkoq7FkQRi9DSl9DHm9R+QNbiDUoI69US94YWvsfe3ABh5275ZLbx9R4njmVhhGjx0DgobpiI7XhuWW7MNkM0NaxJqr9+S\/S6ifVJ3IHM5VGSw2uLZwMopPjVOZrE+JKTP0I5AGUsDyGkb7swWBoLiKspY\/pDSs5gSNc7RyEHqPZYT+YKq0hJFq7x+VrDK+\/JbEGAyAhlAXo76W7fCGF9kCaFnZO2hzxeFPFC3T7wQv+\/4CZnPbUb2ba1RPf0EHDfxI1jdlzYEk60jTkFaz1vR+ObprsDBSeNE3XRVHkAocSzfkKv0GU5r4VJqGbzgWyxbn6f4N+6STIy5hJqP4KW7sndtW195aH7xxR4VjduunU9Fjco3W9kssoFuRyX7FfwQLwcyMgIhi0qYNGFQXuoPLyMQtPIYBQUBELqQ02rSuXOJ2\/vJqalqbuo0dQFA5fQYCFBwEMARC49updHljezs7gCtK9mtlft6ULxZYTiybEI2QhqisHcv2Cf8CkJGcEWZ6I5AYKP6++v\/zC3lyq73pNJ0x5RLUHjo19DtlLTGOG7yp9g1ezCYs\/RI7hb40tsjuVZ9pVzdNaNPqftN71qgJJCiwwU4cmAvfBkdUL3lqWGFshsP\/jOsauNWhCOZE7Aq6zaJCoZhEohV4SFhXlrDpnjuhXkYsfgg9GRAPN7wcpJkyM6RTA\/nl3D7nbtSVei+JN0RetYMS1Ebllfz3w+DgCB+HT9qIe78fUWxPlL6HjkSlCrW+X8v5V2qh8WLJ6rQwhB6OTjIHMaXi96fPVKSwnxMpK34jwQ3eDesKvbRMPM+lXHDfUVKZjPro4\/Plt8Wh\/M3zvsBB6ZfHuaJaqRb9BJNh8\/F7jnDlIJUHWGKo3El1F8pS2PUgcRjQybqpjNbe6KupdQRRrdNy0J577sNP+G9PenY\/e8l6Hv7OMwccLo6uvT6y1p0PaGeI\/DgeHYAIuH8TjxF6dWZk5ysyCQgLPxuU6nwfHqAZrQsCgtzZ\/tI4+th8eKJ6jT0325s48tzQcNTsC1ls\/I8\/Th3nem+kjbyY\/PfWyjPVV5WfaRtaJ3\/9zxSl82ICCDUidAnQweNeAIIs5i5DZiLBDSVbdORHilZaiw1wT3UunVr5f6+fPlyXH755WFLq2xrcQrwpSQQ8aCTcpccSDzrWLuDnngTJkwIKwpkNplupelxcsOw4w5Bx5iYWY4wDOfv1KZeqbB1uwVZhedTguAxRE+0w7EihfPrEgjTBLz5zm+OQ\/\/txjauQ6QL8To1W+fl1S4OKUblePRN\/g+OUwuocP4flyB\/fnhCIV2J2vjW57HvzcdKRcM2vfNF\/Lr02bj6gXCNZmH6tP4c\/Oy1UEIhu2de2TYdq+vt3LkTV155pWU5DSvv1cq2Fjvey+9hACIONXoxZDE7sbDy1VdfjXfeeSdUod3pJHo7Owlk\/FkF6HJCfaVHYKX59u3bK+S2u4ztx3XuHPTv7A4MbZuDIUOCSkS5Io0\/b15GqP1NN2Vg\/vzsmGiJRPtNGTeFwOHTnE\/tlhk1X7qk+3D80ociSiAyuTFDmZXHaSyeqKwP03zUolDiIM4dKXmQGWPiteliDabT3dLpsWrm3k4JhE5jjK2hc5mxNnS81mL7AsW5gakVRhYjc0kVMXrjPfjggzCr+O6ULisAoWVmxOy38e\/X\/ozkQ3udDmfZbusnnyDABLWDgfpvzzoaizDL8Zh5eSNU+0AgAwcO9HbV1\/Ek5dywWfUAHrniVPR4YE7ExMokywggesyLHttidV+PhdFLRehLFulHSkbIv63aRwMg3dFdgTOlOwYRWl3xCKaTIwodxCiJEEi4TxhMR4MA\/zv++OPRqFEjU4A5pgCknN\/t0HTxzLAtEgiNGm2n5GDy5HAJJNIaRQJhkuQuXTLQtq3zvhXFOyfzVras33rSIieJiIxrtNt0bdAmBCB6gibjOLEG0zGKVtzTZWxKGBdddBGef\/55dOjQAQcPHqzaAEIxLyUlBdWrV3fyrlaqNsz3WVQS4W5LG824zLo+eHBJ1nLbTl6DcueAHYCIBELCIgFIvIPp9GjctWvXqvgXRuXyKF7lJBAiJxWnRUVFKhag9fEnYOnh07HshxJ\/AV1JWu5vkYMJCQQ0AUsJA7suBBC2J4joZRLs+nm\/ly8H7ADEKTWxBNOxiPpTTz2l8oXoR3sBJSb5NgIIC6YzEE8ymJHOeK3F6Zrj1c7WE5VKVILIGWecATqDXTNuLtqedCrm3Hq+okHc2vm3OJbFi7hoxtGLd0usQe\/eDbF\/\/36sW7cOqalrcPHFFymxktfWrTl4442gm7XcZ4xI377r8OSTV6r7RoWXU7rMaKHCTdcxGceWl1mS1nCuWPNtOqW3PNsZ9Sf0PWE+EgnCc0JLeW+6aILpGL5PANGtmmZrK++1OOGvkzYRAYTAwfRtK1aswBVXXIHdB\/y4d\/pLmD\/lLhzfvKQ2rhNXdifExNqGYCdp6PjAuPFI\/9ixY\/HQQw9h4MAJmD27IX7++SN8+GEQAEUxvG1bNYwcuRa5ub3wzjv7MWvWSKUw5qWP6ZRGO1rMxtYd+SQc3AxQnNJQWdvpNWAaXDkSua9NUomcf\/3HzFCFOydlMct700UDIE6fQXmvxSlddu0iAsgPP\/ygsinp1zc7\/fhoTz0sGHKGcmNXm\/C19WjXtHalK1Mpjj3XXnstZs6cqXI70MHnjjumISVlCHburIFvvvkGXbp0UY5Zhw79FWPG\/JdaE8FH2utp7OwYavW7FS0yNqUlHhX79OmD0aNHK\/GWGn2CSjysX9HSXRb9dHMt6+0KgNDRzMy8G4mniZgGsMpIIMaFRiptKW2dhPaXxUtpNiZBIDMzUymvxEOQ7bhpzz33XNf3JeNUNPTb0SJjixQiAGI0qd92221hZ+doaKkMfWi2PfJLDhr2nYB9\/3gKjQfOwI4pl6J6u3MS1pEsFr4mvAQiDmOi3LHKWxALk8qzrzFdXEUCiBNarABE55kcZwh+sYBZeT6HSHMZq9OJT4hT+hJ10x2TEoh8IfliGs\/dZvExTh9yRbTT18L5+aKZHUmsjirxPMI4pUUycRslECP\/OB4vXYNfETyuDHNGAhAGaeuplZjjzUkql4paV6KCodKB6HZrccvVFYfy+x2jxuH5Nfvw6NVtQ\/oPMpyWmEfe2Yixl7axLCxVXg9GT5svc+p6BN4z1vswKjTZRvQOens7TbrZZpf0dXa0yNhmRxgpWmQHLuXF48oyj9WmI3gw5ZKewpY5nCTjrE6\/LmnL\/UjHRD4D+nYwzynzmtaqVQss\/7B582YlFYrVhWNJOyf8OqYAhBLHZ599Fop50QFmm796WGpDKlDf\/nov3h7eEacdV8cJr8qsjVFfwIkk8a0kseU9sXCIdEJFnNP7TomPhhaObQYSuhn3WNGBOOVjpHZmm05S6prlvyaIXGCQRIzBbZS+p0yZUsqvQ+jgnAIWuss6f9etNBs2bAi1c7LWhAYQLlBE7R49eii33P79+4fO2WZiM1MUXjf361I1bp0wy2tTNTkQqSKdcCTWhEJMXavXHjZy2vh7JAChZCEfFwI367mIy\/rDDz+MyZMnq\/B8ZnDPyytOlDVuHHr27KnqLtPCp3+s3IJhIrwlITOu7oMgXzm5x\/O5nLmNyZUrkwSSCAyvyjQaAcRN4JwZ38y+2nYAYkyMb3aEoXMfP6R6agvxLKVFzyiByBGmV69eocJTVU4CcfJi04w7bfGmSq0DcbIOr03Fc0D3RI1G+tCPn\/pXPhKAmCXFt8rPYaUboU7LA5CS98fWlb3iXzWPgqrCAfqG\/P7j56oaXkrDFrbLLgsdiExqNCxE0oF4Eoj2qIyBRfypx6VXILfDTXjmptOVlYXHmKeXbcbM\/ic5Tqps+zZ4DaoUB8z0IVZ1dK0YY6V4pBQilXPYVxLm82+jctVKAtElHJmfRxseYagX4d\/Z2dkqOdDIkSNVoiceYeikuGXLFqUXuf322x3HUSW8EpVMEvBgcWHRefDevaPG4l+bf8WbC55G6yb1yw1AjHkaog1qi2Znigh799132zpt6Zp50sx\/0w3erdk3GjoTqY8RNKQWTbRriLTpxBojCbBpfRkU7UTl0O+YABArse2HLTtxw\/CxeGXONJzUqnmZA4gAGZFcNqLRU7asn6lTAHHarqzpTYTxy8MKkwh8MKPxmAAQEduMmaUfmfw4Pt1ZHQun3l0uRxg7ZgrAkN6pU6cqyYnmNVqLxNR26qmn4sCBA0rEpKmN9++66y717ESSMc5DTfvs2bNV\/VJezNsgEohRqSaaej0TFa1XwkMC39y5c8OkEX18ZgejuEv6dJqML5dxXlEYRqL9q6++UusgbTQ18qLZkY5Nwo9jQUKye08SCUwSdS2u6sLwgXAzjp7yNBb8a3+Z6UDsjgFOAIS0cpPQnMYzqziUMeKVm4ggIXZ+2ZRWACJnW8nRQfpkDCPQ6LRzbkbWci4BDAE90iFHHaHR6DNgXGekPjrtBJBx48YpoBTaSS\/HZ0Yss4Q2ibTZhNZE3XTHtARitrjyjsKNB4CI74rxiKFvNKcAIoFrujRAPZGZpKLTTl6KU96ZZ54ZkmjEWU9oFH8b3XmPfXWFtgCg1BtxIoGQPmbJMqNJ9+1JRPAQSc8snF+UqMZ1mbmys41VwiajA6Uc8Zlf5qWXXlJSrei5+G5kZWXhpJNOUmVP9EvelUjJyBMVDCulGbeyAYhseOZG4aaUI4IdgOhJjZiVihKR3ofJjvTLzE3d6BJvdM03k550+qoigNBl3cqV3Xif\/JVYIz4LAYLx48er58XSInpuFjqXRQIQ9iPIE2x4hGUtJWORqWNWArFadFlmYoqGmU6OMPGUQMR0Z3fUoaRiBD++oDzGMDdmq1atSulsnEbV6p7CujnRA5AbSrmLuwEQHTCMG52xLmeffbbKJSNShQcg4Tu2lARitHbo5\/1Y6sG4EZOdWGGMugaKsnq8ghMA0fUXoiuQolr6b2bHD5FGIulAKIEYN74ch8zoN+pAjKEEupgrugwqeY20exKIcwlEpA7qhSgR6scNCZbjUZdSCo8sAiDTpk0LKcDl3daPmVVSAtE3eWU4kxn9QHQRX9dH8MHR5KtbYZwAiIj3zz33nHpxunXrpvK\/Wh1N2I7Xeeedh08\/\/VR9+ahUFmsK6WjatCn+9a9\/hfmBmIGw0WHPysfFaIXReSD8MdLuAYg7ANHfe92Vgc9XCkSR16eddpoqCuUdYUo4ViklEDfSite26nLA6iPn5ghjrAmjH9d1ABFpkNymZGKlRPV0ICaKn\/LWgVTdLeGt3A0H7FzZjWPRO9WswKVR0pWjpDHfB+djWw9AIkggRtOhlESItSaumxfDa+txwAkHKsMx2wmdTtok6lrCjjB6LAxLIdCuTZHs3XffDeomQ6IAAA1OSURBVMtQ5oQhXpuqxYF41jV2yrlt27ZhzJgxGDFiBDp16uS0W6Vsp6+ld+\/eyvEwEa5SnqiSRIWaZwEQAoueXCXahVk57HA8voCsEKczjxGOdEdPS0uDzxesQcOLpSpZPY4e53rpSav2VvS6bc9x3PaJZ3uJKqUozkvE8XjOEc2z5bPjRl6zZk003b0+Bg4QDKdPn54QIFJKicpNTnPWLbfcgtdee02FJTOGJFbPRbsKa2YinN\/vV6nh0tPTkZqaGsbmBQuAj4vTbtPxLzMTiNTe7C11255juO0Tz\/ZUDjJhMLOL8\/9SLzyec0Szm+XZ8aVv0cI+j0c0c1SVPgThWbNmOU6FWNF8MfVENXo\/xiOM3q7CmlsAEcaxEPbEicCmTUDz5taAU54AwsLcBLR4Aw7T8emSh6ypsgCI0\/yfFf3SV+b5E00XUm6u7EZQMrpty+8Eq6uvvlo9Y27EZ5\/Nx5dfpuHw4YA6tphdNWumYu7cAPr0OYidO3eiWbNmjtyHufHctBdAsOszfHiqov2f\/\/QricWuvb4mq\/Y5Ph+G+Xy43+8vVd8kljmMx8NoNleivfTRrNHYx8pRjGZhei5LnR+3cyUaL0sBiHGjCwOk0n08kuSYVViTeQcOHIibb74ZOTk+vPFGHWzY8DuOPz7laJ2WA5bP4vzzW6JTJz8eeyxH1fKl23iNGjVsn11BQYGr9hwwUp\/rr09HRkYAnTvnq7nHjGmCO+\/ci59\/PoKxY\/PRurUcOqxJMxt\/TJMmWJOaihF5eeh98GCpzm7Xobcn2DZo0MCWV5EaJNpLH9NibTpXaQAp68JFOjidf\/75OPnkk0OZz+S3v\/\/97zh48GxQx8HNOHRoDpo3b15KB6I\/x6FDfXj5ZR9yc\/Owa9cu2\/bSl0Dmpj37SR+\/vznatw\/Xy\/B3Ah\/plr\/Xr\/eDEgnvz5njxyCbtFhmNC1MTVUp+jr7\/aavr9t16O3r1KkTpqCOZnMlOoCQ\/jlz5oS8i5k0mekQeMnx3Vibh27tEiwnNYfEI\/jKK69UiZdFEpG0iYzMFq9lqwjdROOlpRUmHpKG\/jJ+8MEHqmDPCy+8gF9\/\/VWFuVPpxhoavIRxt932Ch54oLPSa9x\/vzOdBsGGlpnvv\/ejRg1zpavZxnCrO5AjDBW7S5akY82aVFxwASKCgszB9rffXgI4POJkZQGGyG\/XSlrS9Kzfj0sslM3xWrdbCYTr439OLt2SJu2tjqtm45n1Z7tIY1BHJXoqPSLXGFxH4DjjjDPw73\/\/OyyyVo4wN954IxgXw2A7BuMRIJg3xgxAKBnzYjxUpKTNZikKnPCxItqUOsKUVR1cPiTa63nE4HXCCSeobFlyVuTv1103Fjt2vIIbb8xQJlq3G7ys2wuAPPPMAfzlL01KmZHdbtakpOBLTAWwXG7XwNwX+X4\/nqhkACLKbScvdZHJyY68cXqZ9WffSGPwAyXgzXdPSjWYHeEphfASqUSCNimBXH755eCeERd2OcKYAchvv\/0GiafieGZSSEJLIFwUGSCM0h9grDoQMkZSJXJcIrVeZZ6\/d58EdO7XGQOLxXz6N+zLy0OdtLRSR5j+JuL8Sz6fZXvjy8j+RuUjjwqRrs3FX9SJk4DZZx0Ej076Zdbfag2d1\/sVePTo4VNHCH6tu82H0nG0bZvjWBHMNVO5unTjRpzVqJFr5XFZKVETTQIRAIkU3s9nLaUv\/\/jHP2LZsmXQJRDJ\/yK1Y0TaEAkns1jkkYhss3ctoQHEztTq9ItgxRg7AOmxbRsO9u7taJqf9c92cY\/j27Rx1JeN2N+ofHTT\/\/v160spat301+mfNau+onvWzAbI2BjAmPfXoWfP6qHxV69OxRtvpCnlbIsWgZCOhUeo119Pw7aVPtz088+4ubAQe\/bUVkcrKpVFF2NkCtf9zYEDSnk6aXAmXnzRFxLnHTNQa5hoL71xjboEYvYRffbZZ5WOhLoOXroEwoRBet1lRmrTd0rSLdCnih\/fiy66SB3b9Ty4Zh\/lRONlKVd25t3UU7VF80JZAQiPMLt370ZSUpKSKB599FH06dNHNdeVqEziImjvRsmpKwcXLkyNSmFpRvvJJwfNxPyANGp0wJXi1a2CsyaloMHAM70P4skr62BAIICHAgH86U+p6kzPLzulGiplSQ\/pWtAN+MHvx9JAQPGVZu2uXQNYudKn\/s82vNieY2R0DWBcQQH+90AjdMsClpul73Lx4BPtpXextHJvmmi8NDXj6ine4sVBIjHFPro8E3kff\/xxzJgxI5RTMlpHMp0+XX\/ATU+Lh1FJadXe6OlqtW63Ogq37QkWU3w+dKMSkvoRAyGimBQFIH9e7\/ejh8+nikpfXOzyLwpE6iLU30fBg0ekFZOAzOXAnXv24JoGDdBWCxGI9lkn2ksf7TrLo1+i8dJxVvZYdSA8HtE9\/uuvv1bPwei1GG8A2bkzVVlmCCBWWnq3m5t0u+0TTfsH8vMxq0EDVQgpUqV5eaE5x1t5eXg+PV3lAqWru1zEjqwVwIrMICBxzE4urTZ2GyfRXnq79VTk74nGy0rtiRqLhyVNai+95MOwYT5l3tW\/2PrGc+MlKgDipk+0a9iYkRGSJuxeaJnj2WbN8GVaGloDeCEQPLYM9fnwss+njkI3BgIKXHSaykqJakez97s5BzwAcfBmOPFE5TCxeFiKJyq9Q6l8HDFiX8i7lX9HGp9Ky86dzZ224kFTJBa5Hd9sHW\/UqWPqsSrzep6oDl7SCmqS0ACiJwA28s\/tEYbmsOuuu045jRlLEXDsSJ6o8VCiik6DikRKIWKRuPhin3JQ42Wm4KSOgfqTpUsDSgFpvNwqRcu6vdU6Ir3\/nidq7OjgxcIEeejoCENvPNq23QQIieuvBM059UTVdSPR6A\/Mwv8JCqJU1F3JzcaXXCMmVmLFsHjRZPUKux2\/PGiy226J9tW0W08sv1fpWBgrxlm53Vq1FybSmsOLtU+ceKLShVePxo1Wf+A0GpfmzhdfzMEVV9RWDlg0\/c6bF8D48YGIitfy0IE4XYMASLQ0lZUOhD53Dj3ZS0UXc00Wgdemr5yuNNYbRBqDsUXFGRfUu+nFwkQHm44kEDKYEoXbgsx6aUAnnqgEEKbR53+86O\/w2WefKbd3J4lq3LZ\/\/fU6+PjjNujUqZnSj+zatRP9+q231H+UB01u1xArTXRwat++fXRvT3EvMwmEyY8mOhzVLEbZhSd7KLGScbpIY5A2cX\/xYmEcPiiTZo7NuNEki3EDIBWVFs\/v74T8\/M6oWXM1UlOrXko+Ovfxv1guMwBJNAnEi4WJ7g1wJIFEGtpY4Gnq1KmheAwjgIgUo0ct6nqVikjMGx3bjp1eTN4bawLfRNeB6K7sXiyMu3c7ZgCJNJ0OIHqcDfsYPVHdke21rkwcOJYAhHw1BpR6sTDWb1sIQIw1cXXJwq0JV6bTAYT39FDpaI5ElWnTeLSUcCDRAaQyPctE42UIQLjZGW7MUGNjBvWyyhFiHFfP+qQncra6rwOSMfGz\/hXRAdCqj9v7bsfnS2qWmT5eazDW2tVzzkaaIx6bJ9Fe+nisuazGSDReKgAxmmkpfehHDLdmXCfMNfqJ6JYeXUfCscx0JyeeeOLRPKmjwIp5xiORkX4pVaFX2dP78G+zsazus2aOzh+78ZknQnfSE+kr0rHO7RqWLl2K7OxsZTKXufr3748ePXpY8snJc3LSJtFeeidrqqg2icZLUwDh15XmU1GIxhtAzPxEzKQRScAiG0PAhPeZb9JOKSsvgZiQWW1v5syZyhztBKSswMvoUGc3PiUgpklg6oLRo0erTc4xrEDTzGHPbg5jH5Eo3fAp2k0jL\/2xUCEuWh7Eq59UqEuUI36pIwy\/WEx8wq+XZE4y6jLixSw7Kw1TvtED1ggUcj9SgiKdRn0zmfXhJnNz35hRym58aW9MWm3nG+NmDTpNumRDackpn6J9rhVlgo+W3sreLyEr0+kitpyf5V6sVemsHpgRmORYI9mt6SXJr7XZfQKLk42hSzZWGzYWAHEyfqwA4mYOo57FDUjFsrE8E3ws3AvvGw\/TevyoiTxSmZpx7RYRSbLRlbpWX2K7I4xxDLd6FrsjjNPx5XhhJoHEaw2cwywlpZtjkt3z8n73OGDkQIUUlhIirDxV+TWjwpPp8nURXL\/PnJNWSlTZ+MYAQCulpZWy1Oq+njxX1z3Y+boYAcSuPfnjZg1mBdDt5vC2hMeBWDhg6souSr5YBnbS1+oIw766Ekk34+r3rfxK9PtCh1lKAbdjSftoxicdZoW74rUGnUeyZjFte\/43Tt5Gr000HCi3wlLREOf18TjgcaByc6DcCktVbjZ41Hkc8DgQDQdMAaQsCktFQ5zXx+OAx4HKzYFSRxhRTLZr165yU+5R53HA40CFc6DcCktV+Eo9AjwOeByIOwfKrbBU3Cn3BvQ44HGgwjngOCNZtCH9Fb5CjwCPAx4HyowDFeqJWmar8gb2OOBxoFw44AFIubDZm8TjwLHJAe8Ic2w+V29VHgfKhQOOJBCzmIxyoc6bxOOAx4FKzYH\/B+\/hlOR9DfWAAAAAAElFTkSuQmCC","height":164,"width":272}}
%---
%[output:40267d19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:300021bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:76af1ee0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:2806bc31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:669feba5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:368919af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:7575622c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:8e498ed9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:7daedd72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:23a2f8e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6757de8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:297a0690]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6ff3f15a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:77a8cec1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:803ba2bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4225465a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:43feaaec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:841cf8f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:713e3948]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6fce8c1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:015bbced]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:38dc8002]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1251c42f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7fa8891f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:30485f3c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:65c6d840]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:155410f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4570c4f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:919c996f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5bf15f4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:67838364]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1f67c1c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:601b34ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1757681e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6f129a6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:883f1fd6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:78e033c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:59287da2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1ce6a2da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0b37b71a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4e5b1aed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:83e132b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:816ce839]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:08af833b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:45825452]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:71eb9b3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:979a5aa2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:254a37f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:80e52097]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5994e055]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7f107dd8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3340cb3a]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using . (line 508)\nCannot delete 'expS_br' from the table because it does not exist. Assigning the literal value [] to a variable in a table deletes the variable. To create a new variable in the table with the value [], use T.expS_br = zeros(0) or assign first to a temporary workspace variable."}}
%---
%[output:585aeada]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using . (line 508)\nCannot delete 'expAclap_grAE' from the table because it does not exist. Assigning the literal value [] to a variable in a table deletes the variable. To create a new variable in the table with the value [], use T.expAclap_grAE = zeros(0) or assign first to a temporary workspace variable."}}
%---
%[output:7daca7d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:691d70ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:56b72ec6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4cbd50fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8767859a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11b066ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:009c4dd0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29858109]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c609ceb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:43878f87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:92137970]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0a50c7ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9533c316]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:96569083]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f72c21c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a64f282]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:819da75a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ad1ea29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ac26fcc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:21f31491]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47f7eba9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:184e0113]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b5c357c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9dfc1fb3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:72c69771]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:92bbcedd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:40aa086e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18c8d889]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:53b4c1d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:77ee198a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6577f61e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:21cf559b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6ba35594]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:914fd866]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4a8a6645]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:94b07b84]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16cee0e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:07b88436]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:94d83b3f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7261d18d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6c4e33ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44f1c894]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43e6e121]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ca5a6ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32a14f56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71e76d47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ab5e6f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34001b96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f3e610b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:106c22b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:782857ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:41d65df9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1891fe05]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:38c80f41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b079263]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0874b400]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:881ddcf0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44031ae0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34396911]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:27c4bfa6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4a0ad7c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65e8ab02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6d460f76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:164e8471]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:50dc1af8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:5f574c25]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:450a624e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:74e66a45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a5fabc0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7675f741]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:24aa0d05]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57521f4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:435b04ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5834e47b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:51c0283c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:03ea7eb2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8cb754ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:2521d800]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96422d9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6b93ccef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1ff1e6a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5a750aa8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:3a652dd5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a45a07e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:440cb725]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:03527de8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:62b980ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:145cb312]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1fb44726]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9a72e9db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:970d9f46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:627b9ef2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7bb8a256]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:047d3508]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:53258313]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:00eb0146]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f73b0d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:61911cd4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9d41c5c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0e4470d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4e2b6092]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f127a96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4b83007b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0484a902]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:11101af8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37495a31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c309555]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:74ccd857]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:14be6490]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65bbdc35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:196cb9c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e09b832]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51081c51]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d738641]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f359d5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4481adb1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:40ddc0bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:547dbc7b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3725f127]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75f627c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ccb3eb6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5df4caf1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03d4c46f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d7ad6d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:27351042]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4a092aea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1540ae5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:81dcf9e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03814767]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c0856bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:644bdff9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49c60f87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:08f093fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6afaeda3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1537097e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:50a0bede]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3bd90208]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6bb0d786]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4925dd53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f550899]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:77943ab5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5744ba15]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8431292a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c856963]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5c9fe296]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:11153c03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2693967a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:39a6b9c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:64f431c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:72d925c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9e3df2f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:173985e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c0fdeaf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3ce2d07c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:21497ac0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:02bba6a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ee0f6eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f63a72f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96ec5789]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a056215]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7cca7a2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4549027b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4115420d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8d4861a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0931055f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:350777ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1a4e8313]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ee5317e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:84db3364]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:10d3ba1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61806ca5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:63e1ac1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7aac0809]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2fba3f96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2dd0c217]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0a526f37]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1e8ae956]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b27364e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ba2f081]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3218c0d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:42a6b02f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05cb2d68]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:29988949]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ce840fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:260a6be2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:70604c23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5eda84ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5135b724]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69bd9765]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87e0606b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8a8b83fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9683eb8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59ed2d67]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:074657ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5763ca47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3285274e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:01434088]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0daecedc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4fdfe614]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1bd8ff71]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15ced2fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:990f2ee7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3937d464]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:029880fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:62b5eea5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2fcb0a8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8712811b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ceb2ef3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f3ea777]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:175f3b96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93ea4fbc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c1d9f74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2334239c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37456ea4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06cd1b08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8fa8128b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3241aa6e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:36208924]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97346451]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:63a05afe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9453c2ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9362a015]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18fb473c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c95b9b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0bcd305f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b002ed1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a09ab26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:206afafe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f8b5d75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5648a775]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35eea1c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b8031c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6863686a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:81d0f071]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:92100cbc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0acc5e8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:31975bb4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:230e981d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f869c7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a886ec3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ef3773d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ab8c7fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:676164f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:455d64b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2f3064ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82c04625]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1633a5b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:00329b5d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e26f982]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:123b00ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5c17d0a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3742a577]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93505a75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61120f9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a5714b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:33cd7afd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1884a694]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1db9a7de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57f1e15f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9580b661]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8e0e0ff8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:52ac184d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93cb677e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:40d86b5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f043f03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9ec8596d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:17b5ad5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:876a8eb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9d3b11fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2551645f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e462a53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4603b234]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2cdbab0d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ba43a7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:292005b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:44587408]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d2c25dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a4e322b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7209c994]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4aaaea36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15765720]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cccc1c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36b97662]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47d109dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6d6aed90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:719ba86f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:60041ab7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8288b7ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47a2a2cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2fe9d7c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e1bed74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6fd93c27]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e7457ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:458e336e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0cb08232]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:553ea746]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c44b927]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c30507b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:759a6b80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5c9334e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0587c776]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3971c5ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04efeab9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:013402e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:967fbaed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3116fd48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:906525eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b070030]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:89ad159d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ad2fd08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:750a3223]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57d51fee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1b4f49ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19cd64b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6032e7cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32c1fec6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:887ba1af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:472fba2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ed1c1f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e50198a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d5d30b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:783e6de4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a094f46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ec23811]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7050df14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7da68526]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80d45bc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:66e98f44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:691c68b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1847470b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:5b01f108]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:77ae7dc0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:763e72c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:41e0e013]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4767e172]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:091c8eab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0375b2ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:308a38e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57970a44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0a9eaadd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:83eefd76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:243b56a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7d25fdd6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:62224dce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:71005745]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:83d93634]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9c7cc22d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:64150d21]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:728da710]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1190d817]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:381adbb8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85b092d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ae00d2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b025613]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:46dc88b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8dec3a57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:952f7d45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11939650]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11838473]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:187cb860]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f7ac2c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d78b69b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0f71f93d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e8f9139]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:98ce6d35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ef1c20e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:76107dbb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:07c20a40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1b7752e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:730d9a8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f45c1a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a727372]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62e057cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:378fafb1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45651540]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18518a12]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:500ad4dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26dd33b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57de100a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4aaa3564]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55bda12e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0bc04a08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90cfd8ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e1643bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ff75fa2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8036573e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:385019d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c1aa960]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ba62846]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6937e54a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8174113e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3fd50b0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5413a336]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3fbb82da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f138eaf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26350180]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5cf78d39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2297ef76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4db7e9e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6216ff5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3372c653]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2684268e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58083934]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2698c44d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1995cd3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01fdb1b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57ac2e2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76244fe9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:87ecc18b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:35966a80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:953d2ac2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:31a1500c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:99b1044f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:019ba41b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6830611f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4cc68b26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c6414aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:24a80196]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:67a5241c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37c42036]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:2d76c683]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:398012ae]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"2562×12 table","name":"MLO_result_MK","rows":2562,"type":"table","value":[["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'y'","'MK'","0","-6.2725e-04","0.0028","-0.0042"],["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'MetSea'","'MK'","[0;-1;0;0;-1]","[0.0205;-0.0071;0.0009;-0.0010;NaN]","[0.0464;0.0032;0.0077;0.0060;NaN]","[-0.0059;-0.0182;-0.0063;-0.0084;NaN]"],["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'y'","'MK'","95","-0.0069","-0.0034","-0.0104"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'MetSea'","'MK'","[0;-1;-1;-1;95]","[-0.0032;-0.0075;-0.0008;-0.0053;NaN]","[0.0252;0.0028;0.0063;0.0014;NaN]","[-0.0278;-0.0186;-0.0083;-0.0122;NaN]"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0100","-0.0040","-0.0160"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[0;0;-1;-1;-1]","[-0.0253;-0.0047;-0.0078;-0.0041;NaN]","[0.0106;0.0095;0.0001;0.0039;NaN]","[-0.0629;-0.0192;-0.0159;-0.0122;NaN]"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0110","-0.0083","-0.0138"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[95;-1;95;95;95]","[-0.0429;-0.0034;-0.0042;-0.0055;NaN]","[-0.0249;0.0025;-0.0006;-0.0014;NaN]","[-0.0622;-0.0093;-0.0078;-0.0097;NaN]"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","30","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0058","-0.0048","-0.0068"],["'MLO'","2025","30","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[95;0;95;-1;95]","[-0.0218;-0.0030;-0.0041;-0.0017;NaN]","[-0.0132;-0.0004;-0.0024;-0.0003;NaN]","[-0.0307;-0.0056;-0.0059;-0.0032;NaN]"]]}}
%---
%[output:1c2f3323]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"20×19 table","name":"MLO_result_LMSlog","rows":20,"type":"table","value":[["'MLO'","2025","10","'month'","'Ba3_A82_ae33'","'abs'","'log'","'LMS'","3.1262","95","-0.0726","-0.0261","-0.1190","-5.1766","-1.8649","-8.4883","-0.5160","-0.5154","-0.5166"],["'MLO'","2025","10","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","4.8948","95","-0.0465","-0.0275","-0.0655","-34.7014","-20.5226","-48.8802","-0.3717","-0.3714","-0.3721"],["'MLO'","2025","20","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","5.1774","95","-0.0239","-0.0147","-0.0331","-55.0205","-33.7666","-76.2744","-0.3799","-0.3796","-0.3803"],["'MLO'","2025","30","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","1.7301","90","-0.0053","8.2207e-04","-0.0114","-14.7924","2.3074","-31.8922","-0.1462","-0.1458","-0.1467"],["'MLO'","2025","10","'month'","'BbsG_S2S20'","'neph'","'log'","'LMS'","0.1620","0","0.0026","0.0343","-0.0291","0.1419","1.8936","-1.6098","0.0260","0.0269","0.0251"],["'MLO'","2025","20","'month'","'BbsG_S2S20'","'neph'","'log'","'LMS'","1.0755","0","-0.0056","0.0048","-0.0161","-0.3119","0.2681","-0.8919","-0.1063","-0.1058","-0.1068"],["'MLO'","2022","10","'month'","'BaG_ae_psap_clap'","'abs'","'log'","'LMS'","3.8480","95","-0.0820","-0.0394","-0.1247","-3.1666","-1.5208","-4.8125","-0.5597","-0.5592","-0.5602"],["'MLO'","2022","20","'month'","'BaG_ae_psap_clap'","'abs'","'log'","'LMS'","0.1892","0","-0.0020","0.0192","-0.0232","-0.0794","0.7600","-0.9189","-0.0393","-0.0382","-0.0404"],["'MLO'","2022","30","'month'","'BaG_ae_psap_clap'","'abs'","'log'","'LMS'","4.0497","95","0.0256","0.0383","0.0130","0.9632","1.4389","0.4875","1.1564","1.1587","1.1542"],["'MLO'","2025","10","'month'","'expS_bg'","'neph'","'log'","'LMS'","0.2871","0","-0.0021","0.0123","-0.0164","-0.8308","4.9568","-6.6183","-0.0204","-0.0200","-0.0207"],["'MLO'","2025","20","'month'","'expS_bg'","'neph'","'log'","'LMS'","0.8500","0","-0.0022","0.0029","-0.0072","-0.8707","1.1780","-2.9194","-0.0422","-0.0419","-0.0425"],["'MLO'","2025","30","'month'","'expS_bg'","'neph'","'log'","'LMS'","1.5589","0","-0.0029","8.2071e-04","-0.0066","-1.0720","0.3033","-2.4473","-0.0833","-0.0831","-0.0836"],["'MLO'","2022","10","'month'","'BbsFg'","'neph'","'log'","'LMS'","0.9001","0","0.0048","0.0154","-0.0058","0.2499","0.8052","-0.3054","0.0488","0.0491","0.0485"],["'MLO'","2022","20","'month'","'BbsFg'","'neph'","'log'","'LMS'","3.0163","95","0.0060","0.0100","0.0020","0.3049","0.5071","0.1027","0.1274","0.1277","0.1272"]]}}
%---
%[output:05353ff7]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"20×19 table","name":"MLO_resultLMSlin","rows":20,"type":"table","value":[["'MLO'","2025","10","'month'","'Ba3_A82_ae33'","'abs'","'lin'","'LMS'","3.4046","95","-0.0351","-0.0145","-0.0557","-14.2518","-5.8797","-22.6240","-1.3511","-1.3506","-1.3517"],["'MLO'","2025","10","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","3.3361","95","-0.0672","-0.0269","-0.1076","-7.6846","-3.0776","-12.2916","-1.6724","-1.6713","-1.6735"],["'MLO'","2025","20","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","3.2030","95","-0.0386","-0.0145","-0.0627","-4.0329","-1.5147","-6.5511","-1.7723","-1.7710","-1.7736"],["'MLO'","2025","30","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","2.4870","95","-0.0170","-0.0033","-0.0306","-1.7593","-0.3445","-3.1740","-1.5093","-1.5082","-1.5104"],["'MLO'","2025","10","'month'","'BbsG_S2S20'","'neph'","'lin'","'LMS'","0.6096","0","-0.0018","0.0041","-0.0077","-1.1014","2.5122","-4.7150","-1.0180","-1.0179","-1.0182"],["'MLO'","2025","20","'month'","'BbsG_S2S20'","'neph'","'lin'","'LMS'","1.2605","0","-0.0018","0.0011","-0.0046","-1.0851","0.6366","-2.8067","-1.0358","-1.0357","-1.0360"],["'MLO'","2022","10","'month'","'BaG_ae_psap_clap'","'abs'","'lin'","'LMS'","3.7449","95","-0.0090","-0.0042","-0.0138","-12.0274","-5.6040","-18.4507","-1.0902","-1.0901","-1.0903"],["'MLO'","2022","20","'month'","'BaG_ae_psap_clap'","'abs'","'lin'","'LMS'","0.0689","0","-1.2006e-04","0.0034","-0.0036","-0.1501","4.2058","-4.5059","-1.0024","-1.0022","-1.0026"],["'MLO'","2022","30","'month'","'BaG_ae_psap_clap'","'abs'","'lin'","'LMS'","3.1722","95","0.0028","0.0046","0.0011","4.0663","6.6300","1.5026","-0.9146","-0.9145","-0.9148"],["'MLO'","2025","10","'month'","'expS_bg'","'neph'","'lin'","'LMS'","0.4663","0","-0.0042","0.0139","-0.0223","-0.3296","1.0841","-1.7434","-1.0422","-1.0417","-1.0427"],["'MLO'","2025","20","'month'","'expS_bg'","'neph'","'lin'","'LMS'","1.0456","0","-0.0034","0.0031","-0.0099","-0.2657","0.2426","-0.7740","-1.0681","-1.0677","-1.0684"],["'MLO'","2025","30","'month'","'expS_bg'","'neph'","'lin'","'LMS'","1.4352","0","-0.0036","0.0014","-0.0086","-0.2738","0.1078","-0.6554","-1.1076","-1.1072","-1.1080"],["'MLO'","2022","10","'month'","'BbsFg'","'neph'","'lin'","'LMS'","0.4407","0","3.6164e-04","0.0020","-0.0013","0.2434","1.3479","-0.8612","-0.9964","-0.9963","-0.9964"],["'MLO'","2022","20","'month'","'BbsFg'","'neph'","'lin'","'LMS'","2.9793","95","8.7158e-04","0.0015","2.8649e-04","0.6228","1.0409","0.2047","-0.9826","-0.9825","-0.9826"]]}}
%---
%[output:9a1762cd]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAYUAAADqCAYAAABN9llzAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQvYjkX6v9VXKeS0lFCkrI6ujkSFlvZSaRVb1FUonZUOSKhQKFaLtEVZh9olu2rbDptOaCtkbdraVBSVQ638KWrJF\/\/9DfM23\/M9h5nn8D7zPN891\/VeH+87z8zcv5lnfjP3fc89lXbt2rWLODECjAAjwAgwAkRUiUmBxwEjwAgwAoyARIBJgccCI8AIMAKMQAEBJgUeDIwAI8AIMAJMCjwGGAFGgBFgBMojwDsFHhWMACPACDACvFPgMcAIMAKMACPAOwUeAzlGYNGiRXTJJZeUkfCxxx6js846q8x377\/\/PvXq1Ys2btxY+P6aa66h22+\/XfxfLWfUqFF08cUXG6G2Zs0a+uMf\/0hz586lVatWiWfr169Pp512Gl166aV03HHH0V577WVUJmdmBIqFAKuPioU015M4Am6kcO2111L\/\/v2pUqVKhfqffPJJuuOOO8q0Jw5SKC0tpenTp9Po0aNpx44dnvJecMEFNGTIEKpZs2bimHAFjIApAkwKpohxfmsRcCOF008\/nSZOnEgHHnigaPcPP\/xAd911F82ePTtWUsAZ0D\/\/+c80ePBgAjkEJRDDsGHDqGrVqkFZ+XdGoKgIMCkUFW6uLEkEVFKoV68e1apVi7788kuaOnUqHXvssaLq9evX05VXXkkbNmwQRCHVO1F3CitWrBAqqXXr1ol6OnfuTH369KFGjRoJVdGmTZvo+eefpzFjxtCWLVtEnjCqqSTx47IZASDApMDjIDcIqKTQvHlzobt\/4oknyky+yHP55ZdTy5YthUrp73\/\/u5A\/KilgN\/LAAw+IsqCyuvXWW6mkpKQctgsWLBBk8d1335FzF5ObjmBBMo0Ak0Kmu48bryLgJIVu3boJ28FFF11Ew4cPp3333VeokjB59+7dm7799tuCGikKKWCCHzhwoNgJ1KlTh6ZNm0ZHHXWUa+ds3bqV+vXrRy+99BJhN4NdTNOmTbkjGQFrEGBSsKYruCFREXCSAgzM+FSvXp2mTJlC1apVEyv4V155heCVtGTJEpo0aVLkncL\/\/d\/\/CZXUu+++q7X6Hzt2LD300EOiXngpYdfCiRGwBQEmBVt6gtsRGQGVFM4991yxS8AKHt\/PmDGDateuLfT+MApjhf7000\/HTgqdOnWi++67j\/bff39PeVTvJyaFyN3OBcSMAJNCzIBycekhoJICJucRI0YIw+7jjz8udghHH320UBtJXf7DDz8cOymcdNJJNHnyZF930\/vvv79QL5NCeuOFa3ZHgEmBR0ZuEHCSAlbszzzzDA0aNEgcYKtbty7NmjWLLrvsMrrzzjvpkUceKRiHo9gUYCeQaikdm4LMyzaF3Ay9XAnCpJCr7qzYwriRwocffii8jbZv30577723+Dty5EiCEVpV40QhBaijYMD+7W9\/KzrA6wyC8ywDex9V7PFqq\/RMCrb2DLfLGAGVFOQkL88lgByQqlSpIuwLJ5xwQmykgHKXL19OPXv2FOcfkNq3b08333yz8CwCGSGkBnYpDz74YOG0M59TMO5ifqAICDApFAFkrqI4CLiRguouilY0a9ZMeCJBdaOzU\/Br+SGHHCIM1kceeaQwXvOJ5uL0M9eSLAJMCsniy6UXEQE3UkD16sEyeCXB1oAdA84V3HjjjaKFXuojXVJAPoS3ANHgTATHPipix3NVsSLApBArnFxYmgh4kcJrr70mvI6QYOTFiWIkr\/xuMZTc5FJ3CurviJKKXcOrr75KH330kSALGSW1S5cudMopp3CU1DQHCtftiwCTAg8QRoARYAQYgQICTAo8GBgBRoARYASYFHgMMAKMACPACJRHIJM7hf\/+978ifEH37t05bgyPakaAEWAEYkQgc6SgBh\/jEAExjgQuihFgBBiBrN2ngB3ChAkTqGvXriL8MO7U5QiTPI4ZAUaAEYgPgcztFCC63C04SQGugPggNWjQQHw4MQKMACPACOgjkBtSABkgdv7ixYuF9H379hUfTowAI8AIMAL6COSGFOSBI4RKxkEh3inoDwLOyQgwAoyARCB3pMDGZx7cjAAjwAiER4BJITx2\/CQjwAgwArlDIJOk4NYLUn3EO4XcjVEWiBFgBIqIAJNCEcHmqhgBRoARsB0BJgXbe4jbxwgwAoxAERGwihTU08q4eB1x7\/fff\/8ycKiXnquhi1l9VMRRw1UxAoxAbhGwihQw4Tdq1IjOP\/9819hGfjGPmBRyO0ZZMEaAESgiAtaQgvOUMm6wWr16tQhlIRPy3HbbbTRo0CBxBaKaJCngwFqLFi34nEIRBxFXxQgwAvlBwCpSUCd8kMJbb71VRoXkvBHL7wpFPtGcn0HKkjACjEDxEMgUKaiwSFVSq1at6OKLLy5crcgnmos3eLgmRoARyB8CVpHClVdeWYh86qY+csIPGwQSVExsU8jf4GSJGAFGoPgIWEMKED3I0IyJf8GCBYIEnDYIJoXiDx6ukRFgBPKHgFWkoLqkSnsBvhs2bBjdfffdVKtWLUEckyZNEj3hZlPgE835G6QsESPACBQPAatIIYrYvFOIgh4\/ywgwAozAbgSYFHgkMAKMACPACBQQyBwpqG6po0aNEp5HSLxT4FHNCDACjEB0BDJFCurhNYg+cuRIGjt2rLA1VDRSwE1zc+bMoS5dulSIa0dZ3ugvu80lcP\/a0zuZIgVM\/DA0T5kyRcREGjhwIHXv3p1atmxZ4UihopEgy2vPpJFES7h\/k0A1XJmZI4WZM2eKU85IIAXn4TUZ5iIcHNl5au3ateJOapY3O31m0lLuXxO0spdX9q+N3pK5IQVsPzFJLl68OHsjhFvMCDACFQ4BxGhDBAbcJ29TyhwpeKmPACqIAR9OjAAjwAjYjgDIwDZCAGaZIgU\/Q7PtA4DbxwgwAoxAFhDIFCkAUNUl1UZ9XBY6ndvICDACjIAXApkjBe5KRoARYAQYgeQQyA0peB1qSw664pUsw4Q\/++yz1Lx5c+GSi7MZalLlV68pLV4r46tJR15ZmzMwYnytKF5JuvKqcb+yvEvWkVeNg5b18ew1kvxukize6CtfUy5IIe+2BjWMuIwkK09yo0udN9Ih\/6xZs1zJI83Bplt3kLxqOXKizPIkqSOvmmfFihU0ffp0Gjx4cLk7zHUxTjOfjrxqWHy3C7fSbH8cdaukZ9vYzQUp+B1qi6MD0yzD7TIheVYDB\/jcEiYN9bR3mu03rdtEXvT7c889R1u2bCkcYjStL+38OvIiz4gRI6hHjx7lrqFNu\/2m9evIizJV4tC5W8W0HWnmBwYTJkygrl27Ur9+\/Qp3yKTZJrXu3JCC16E2W4AO2w7nFlMlQKcKSdaR5ZdIV14ZUn3AgAE0evTozJOCejJful3L\/pWkgP7FqlL+xUn+rCXd\/pXEcMcdd5Aa4yxr8vq111bVJ5OC5aPM5CWCKCCNoJ2EzSLryouJs02bNsLGooY7sVk2t7bpyCsnj27duhWunnUSR1bk1pHXuZvIo\/oI\/cWkkOCoZfXRbnCzvEOQw0NHvaDqY9VhZZtuVmfI68jrnEihHpQ7pCOPPFKnGmvy6MjrtJFlWV7eKaQ09Cq6oVkSAv6qBuiUuiNytTqGSCeJSPVL5MpTKEBHXjWPjgoxBTG0qwyS1404srozYlLQHhbxZ8zzoTbVha9Tp04iICCMzHi5kE488UTq1asXrVu3rgCsl+tq\/MjHX2KQvCrx2erWZ4KKjrxqnqy7aOrIWxFcUll9ZPKWcF5GgBFgBBiBVBDIhaE5FeS4UkaAEWAEcogAk0IOO5VFYgQYAUYgLAJMCmGR4+cYAUaAEcghAkwKOexUFokRYAQYgbAIMCmERY6fYwQYAUYghwgwKeSwU1kkRoARYATCIsCkEBY5fo4RYAQYgRwiwKSQw05lkRgBRoARCIsAk0JY5Pg5RoARYARyiACTQg47lUViBBgBRiAsAkwKYZGz\/bnVq3e3sFEj21vK7WMEGAGLEGBSsKgzYmnK\/PlEw4YR4a9MIIapU4nati1fhS556OaLRYgIhei2M+58Bk1evafuRkzYBqhx1mIhwKRQLKSLUQ\/IYOhQ75pADD177v5dlzx08ym1pjLp6bYz7nwGcs+fP5+GDRtG+PsTXzeiqVOnUls3wi7GmOE6GAEHAkwKeRkSmGjatQuWZt48ogUL9MjDhGQEz6Q06em2M+58e9DWkRtkMNSHsEEMPSVhB\/ci52AEEkOASSExaOMv2HcFDkJQVUZe1UNlIVUnfk3EBOa365DPgmTathUrYNNJT3dH4ZtPlwx15dHNZyA31ETtNAh73rx5ZXYMuvjEP9K4xIqMQChS2LVrl7jQZdWqVVSjRg1q1qwZlZSUVGQcE5U9cCWKSb5x40Tb4Fl4z540v0cPo0kvUB6DFbjYHemQYdzoGMgNUpATvF8zsFPAjkEXnzAiMdGEQa1iPRNICrglacqUKYIA7r33XnHj1xtvvEE33HADbdmyRaB1wQUXiJVi1apVU0NvzZo1NGfOHOrSpQs1aNAgtXaoFZeWlgqMqlWrFpo0tVbg0EenRQpE1K5t2zJ6ci\/wMelhgtTZUeRJbpPBCGx08HGWicl+69atngu0JInGRD6dvHG8Nzr1FCtP1uQJJIUXXniBbr75ZmrdujWNHz+e9tprL+rXrx8tXLiQ+vbtS7hSbvLkyTRkyBC6\/PLLE8dZvXZz1KhRhTuJ5fc2Xd6+bds2Wr9+PdWrV48qV67siY2XIwxeZG21g4Z6IqnOwR5ljwNsLFUETYwFzRUMtL16xVJnmELillu3DaqaSWey1yJYhz1D1zlLt80m+XTfG5My08y7cmWpmAdOOaWO7zyQZhvVun1JYfv27TRo0CAh0EMPPUQ1a9ak5cuXC4MYJqvhw4cTVEl33nknIa+8Ozgp4UBAt912m2gT0siRI2ns2LFUq1YtyiIpBDnCAGPVU8VvBT4Vb3EaahQiqpRUhweUK9Qt06alVHvKck+dqmXHMbVnBI3JYoBtQgpxk5dJeUF5bcAyTH\/5koK8WLply5Z0++23i\/KffvppMTGrq\/T7779fTMpQM2GCTiqhDtSFeqDGGjhwIHXv3p3QPkkKe++9tzXqI+AA9ds+++zjqj7ym8OhEdIhBIl1W10DcgKdo5yISKB0\/yJdTl4UrQ2pyq2pstO1ZyAfPkFjsljg+r03ahucx3GiHv0wKS8obxCWcEO+++67iwWpdj1GpADd2D333ENPPfUUzZgxg0444QRRUTFJYebMmWJHggRSaNWqlVAhSVKA7YMTI8AIMAK2IwBSgCrQtuRLCt99952YeLEqBxls2rSJrrnmGqpSpUpBnYTvYHSG3hzqnP322y8xGTHxB5EC1Fm2JJAoktMzS8cjdLcM+pr6cqdj1Urclk9ejXDm1cwX5F2ju2I17buC3Jrt9HTHTVBuKZMXRlKGIAxNsTHPHxwSRYVJF\/Kf5N\/9L6\/VvE55uu+OzsZZJ49su0leXdzbtm2UPVKAcNgRwOvowgsvpG+\/\/ZZefPFFuuWWWwQ5vP3228LIjMl63LhxdM4552jhIdVS7777LnXq1MnVFuGWB\/WhbhBRpUqVqEmTJoKsVPWRTYZmNwOTufdosMbec8WByrzeQC+FJ7azqtFRM9+0adOol4fRF5Me9P9Qh5moxIIGUxm5NdvpeZI7QbmlisDLKIzfpZ1OB5+kCJZIb0GFdZfuOUD0oY5uXbc8XQ9k3UlcNx9eiSTMdhatYQuvW6D3EdzcYMydPn26eKh9+\/Y0YsQIOuCAA8QuYu7cuXT99dfTddddp71LgLoJA\/v8888vYxdQJwG3PB988AFh8oFN4ZtvvqErrrhCeEbhb7EMzUHGpaCXAIPQzHvU38cFOEYOk+BHHns6Bca\/DUuWUJ1TTvH0oAia9HS8qSSB+Lll7l5t+sgdkzyF8RhQXpDcbuSGXYFzd6eLD\/pbxystiFTL\/g7rjJ4qA9FSdJy+dA\/P453QcZ7TrddM7nRzQ9sd1Q4StwSBpIAKoZKBv\/2PP\/5I1atXF26pUI28+eabdMQRR9AhhxwiVu46Se4AYLjGCv\/JJ58UB3ukIRtl+OVRXVLPO+884Sqr2hTgJtuiRQs6+OCDxSeu9MYbJTRiREm5OHOPPlpKp5++W02EdP\/9lX0PAiP\/VVfpH\/R7+OFHBOG6JZzHuPTSS3392uOSH6Tw5Zdf0kEHHSTUiUHJbdLDMyB1HXl08wW1w+t3U3l06\/GSW\/d5Xbk7dOggzgsFJYwRnOHxS7vzPP6\/Nb2e2X63E0RQzbsnOx11j24+3XqDW2ZPjkzuFOKGT3UrPfLIIwUpvPXWW2VUSGHyqGSBNvfo0SO2cxPjx9f43xmNmp5QjBmzgbp02UqLFlWmSy6pFxtkLVpso5kz14tdEM6ILF68uFA2XmQQIA7rFSPB5XjDhg1Up04d7R2hV7t05dHNF0b+OOUJU7\/fMzpyO8e716Jh9OjRtHbtWurfv7\/nwgJjaNGi+2jxYu+zNPLhBg1Kac0a\/UVN3NjkqTyQnIV2ZtLaKagdIXcN2CnAgIrTurq7BJSDCf+mm26ir7\/+mj7++GM6\/vjjhQvpmDFjCivQoDxwV4MrKgzejz32mHhOviQop379+rHtFLAi6tgx+GX529+2iZ0EdhRBSWdlhDzOXQjKxaovjRPbwPyrr74SuPodxAuS3fm7rjy6+XTrT0oe3fp18\/nJ\/cQTT9BVV13lWpRUw8kdOHYVUPs6I7TCntGtWzcxbjt08B+7GJODB5vtdHXl1M2n8+7olmWSTwaxjWuH5BfN3qRdSeR1JYWdO3cSjMBSPXTGGWeICRgrRZyOxClnmaACuuOOO+i4445zbR9ePtgenn32WaFmwooXNohf\/epXwmiMyf2www4T38sEUujYsaNrHvx27bXX0meffSbUT\/LAXFSbgpetQNewZbq19dOPYsDAsGWTC7PJgaIkBmrcZeZFnrjsGcAX5wC9bAXqmNTUFMfdZaI8rKyD7A9SRx+X6kpO4Khfp26826hbB8tEQIpYaDlSgGEZ3kazZ88uFH3qqaeKyXfixInijELt2rXFahWrmI0bN1Ljxo2FFxK8gYKSnPA7d+4sjMQghUMPPZQmTJhQjhSceWB8vOuuu8RK9ZhjjqFly5ZFJgU\/zwhzo3CQ9D\/9DgMTBo7bfThORxj9UpPLmZdJVCKUR3mWLFlCp\/g4AuiMDh0nLt2Fkk59JnmkukWHvNq0SWYC16lbLuZ0sDSRv1h5y5ECVvTQVWNwwU0OLw8mbKiLNm\/eLL7r06eP0CtDLwsXULiF3nrrreL7oARSuPHGG4Ua6aOPPqIzzzxTGK9BOqinTZs21LRpU9c8devWFZ5HasLO5bLLLiuoj0wMzUFGYd0oykEyu\/3+3\/9uK\/O1hsNMmGpieyYpw2xsDTQsiOUJBsxrTOqoVGF7QNKxP+jYKZAH6lRVjYP3U7WDIA\/UW9KjGhP4dde5q36R99JLSwtOIZApqDyJmElePIOx9o9\/fE0nnVS7jJMG1O82RpcuQwoy1hFW4Jh8pcvc0qVLBRk0bNhQfI+DajLBNnD11VcLPT4mdqiZ\/FIYI7I0Rvfu3VsYpgcPHizUW\/Igm2pTkHUHGZrjNgqjXhiGdYx10oAc\/Frak8Nmw2wYlFieMKj99MycOVWpf\/86roVgwu3SZYt4H4IcL5B39OgNtHZtSWB5fftudq0PxCNJyJkB7zmcRJzk0bfvJuEc4pb8ynPm18nrNdZw7QDiydmWypCCdAWFakid4GFLwIQMknAGvZM2g88\/\/9w19pGbTQFqIHz\/6aefCkMzdiWY6GVCOzCpO\/PA\/RX2CzXhwBzUWl6GZi9bAYxqOkZh3Q7DCgarFB1jnZsBWbeetPJlxTCriw\/Lo4uUdz4vN22oT7p1270TfuKJEk8XbGmnuP323Xl1yovSap0JPEr5Xs96jbVM7BQkKUDHr07+Xt8DhCBScAMK9y8cfvjhgghwgA2Hz\/BRk1se2BhklFTYMhCdFSG7oXJyGpqLaStQPQlMdI5JDMCkysyjDl4nrHlSeMZdbtr9E8fheYkJZFmyZENmQk0H9WXafRPUPufvrjuFJElB7gIA1CeffCJsCkcddZQ4vIZTzNKmgJ2CMw9+W7BggcirxkFS1Uewccyd29L3AJmprcAvv5unUFYNTH6DJ2sDO+hFYHmCEErmdx3bGfdNMtjrlpoKKcjVvunhNcRfgvcTVFXOuElyp3DuuWNo4sR4D3TBKGxqXEIHwDvr\/vufpNtu61IupIFuB9mSD7L84Q9\/EKfHywXfs6WRBu1geQzAKnLWitI3rD7aM7CiGJpPPvlkQQrqfQrO0Nnr1\/+Rtm1rGdswdjMK6+omQVzwynrggQeoefPmsbUpjYLyJAvwY3nSGEV6dVaUvsmUoRmdYpow6bldsuNmaMY5iDCxj6A+Uj2OoG5Ckuqkiy4aQF988bpp0z3zl5SsoaFDV1Pr1j\/FNjIpHCsenKtAaPGTTjrJ5FHr8uZJFrmLy0vf5E2eijLW5MVGtr3sruqjOEnBTeCwUVJxfkGqnuAhpV6yg4F0yy3j6c9\/HqONcY0a42nz5r6u+UEIVavO+Z\/L2E8nrbUL5oyMACPACAQggDNV+NiWjGMfxSGAelcC7mWQ8VmkoRnhK7zyqIHA1GflaqlhwwbaTVy4cBGtXNmApkxpUC76KdzqmjVbpF0WZ2QEGAFGwAQBLGzTiGMW1MZUSCGoUVF+1z2C7xahUMczIkrb+FlGgBFgBGxHIHekAC8h3aBV8si87Z3E7WMEGAFGoFgI5I4UAFxeD5AVa1BwPYwAI1BxEcglKaA783iArOIOU5acEWAEioVAbklBBTBuWwGC8skYTKorrmoEHzVqlDjohWT6fbE6X9YTlzxwFJg0aZIoFndn4C5hHFAsZjKVBW2TbtMI4w4nB78+K6YsqCsueWzoG1N5VHd2PKs6lni9Uzb3j588tvQP8KsQpBDnQFmxYgWNHDmSxo4dS7Vq1RKhOdatW0cDBgwQcZgGDRokqpN58G\/pRqvzPcosZopLHoQagYuwOrEWUw7UZSoLsFa93BAiRXq+ufWZ7X3jJY8b6RW7b8L0z8svv1y4v132E26Jw\/3UWewfL3kQ\/y3td0cdD0wKEd8OGYMJITjGjRtX5rQ1JkgkEId6Ctvve7lSjdis0I+HlUc9Q1Ls3YGXsEGyYJeHOzy6du1K\/fr1KxyoxHNufWZ733jJ44wiEHpwxPxgUP848Zbnm3BTYxb7x0seleRseHeYFCIOdHWgytPWKFIerMMANvleqpwiNiv041HkueSSSwr1Os+QhG5QhAeDZJFYy1WoPGWvBltU+9L2vvGTx7a+kYslnOr1ekdUvFViQ4Rkt3cqK\/0D2Z3y2NQ\/TAoRJh3oe1evXl0uamtWSSGKPOoLKdUVMi5VBIhDP6ojS5ZIIYo8Kog29A3aYyKPU\/1lI2lHkce2\/mFSCDntyFWoakw2URN5qZXSUlFElcdtawwZ5Wn1kDCHekxXFtlmt52CTeqJqPI4QVRjhoUCOOJDJvK4qb5sU+9Flce2\/mFSCDHA1XAc8nF18OK7rBiaJTkh2KA6sZvK8\/HHHxfuunBOsiEgDv2ISd9Iw7GzvV6yF9vQbNo3XvJgEpX3kKTZN6byID\/uYL\/77ruFU0fQu2Z7\/3jJY1P\/oI1MCobTj+oKJx\/t1KmTuKkOgQSlblB6siCP+ozO94ZNipQ9TnlUt7o0bAphZAF4bhOlV59FAtvw4TjlSbtvnO+Bzrujtlnml67eWewfP3ls6B+JMZOC4YvK2RkBRoARyDMCTAp57l2WjRFgBBgBQwSYFAwB4+yMACPACOQZASYFzd7duXMnvfLKK8Jg17FjR2rdujVVqlRJ82nOxggwAoxANhBgUtDsp08\/\/ZSWLVtGnTt3pkceeYRgXG7YsKHm05yNEWAEGIFsIFBhSUE9AIO4PUhq8DHVS0h2ZWlpKS1dupTefvtt6tWrF1WtWjUbvcytZAQYAUZAE4EKSQpy8peupCAFNZgafO5xjH748OGEf2\/fvp2aNWsmfKVxxB7PQ4XUpEkTTZg5GyPACDAC2UAgl6SgHj7CHagyDhFOH2OH8Nlnn4neeeutt8T5ApACJnr5fxyrRxRGPFe3bl3atWsXbdiwgXbs2EHHHnssvfHGG\/T999\/T2WefnY1e5lYyAowAI6CJQC5JAbJj5Q8VD8Jau6mCVBKQpCDjGLkdZtq6dauwJdSsWVOQys0331zmlKUm3pyNEWAEGAGrESgqKWDFjUl61apVVKNGDaGSKSkpSQwgedeB3A2oFZmSgnz2hx9+oH333TexNnPBjAAjwAikiUBipAAVDO4QAAHce++9QkUDtcsNN9xAW7ZsETJfcMEFIrZJEgZbTPqzZs2iq666it57771ygdncSMGpPsKFOTbEN09zgHDdjAAjULEQSIwUXnjhBaFigT\/\/+PHjaa+99hIXmSxcuJD69u0r4s1MnjxZ3FZ2+eWXx4o6yn722WepR48eolwQAJIa3tlJCm6GZrcdRqwN5cIYAUaAEbAMgURIAd46WGWvX7+eHnroIaGHX758OfXs2ZPatWsnvHqgSrrzzjuFZ08ak6+TFCR54O7ltO4XtmxscHMYAUYgDgRwSTxSo0ZxlJZ4GYmQgjTUIhSzjKf\/9NNPC48e9UJ76PzhDQQ1k0nYW5suuU68h7gCRoARyCYC8+cj9jcR\/soEYpg6lahtW2tlKgop4NDXPffcQ0899RTNmDGDTjjhBAFIGFKw5RJya3uUG8YIMALpIwAyGDrUux0ghp4902+nSwsSIYXvvvtO+PjDuAwy2LRpEyG+fpUqVQrqJHwHo3O9evXEhTT77befFkC2XkKu1XjOxAgwAvlHADuDdu2C5Zw3z8odQyKkADSwI4DX0YUXXkjffvstvfjii3TLLbcIckCYCBiZoToaN24cnXPOOcEY+ygbAAAf80lEQVQA7snhvHhEvcxlzZo1hA\/SwQcfLD6cGAFGgBEoJgKVO3YsqzLyqhw7BewYLEuJkQIOe40dO5amT58uRG7fvj2NGDGCDjjgALGLmDt3Ll1\/\/fV03XXXae8SnNipl5DDy6l\/\/\/60ePFikQ2eR3F7NZn2HdRmOPkMmZM8j2HarjD58yQL5K+o8pTsWTSVNmgQZhgU5Zks9w3wbXjmmfo47dqln7dIORMjBbQfHkY4k\/Djjz9S9erVhVsqOvzNN9+kI444Qnj5RA0\/LS8hxx3DuApzzJgxVL9+fSt2CiCtr776SrSlcuXKRerSZKrJkyxAKE15cHIeqVGM3iiB8kClMXQoVd6zaBLE2KABlT76qHUqjEBZkhni8ZS6ejVVPuoo\/bJWrbLOKylRUtBHRj+n1yXXKAGk4BbSQr\/0eHNu27ZNuOXCbpJ1UsiTLOjlNOSZP3++OKyJvzKBGKZOnUptI3qjrFy5Uoy1U045pfxYy5jRM42+ifXNN7lnpaLtFAA0wkLgRDEOrSGgHE4vn3rqqXTMMceEVqm4XXItbQ1MCrEO70JhmX9RHbAkIY\/fDgBkMNTHGwXEgHM8atLZUQQSTUijp07dyYy0dAg7jCyeGMHIrLqhehWOhQCMzZalxHYKUB29\/vrr4kWQUUlV2aE+gqqnefPmsUDCpBALjJ6FJDGJJtti\/9LjlCdoYsbvOLQZlObNmyd2DEHlyXK0iAY2PZ0Jao\/RU7fuIFmi\/B5n30Rph9ezgRjpELHF5xUSI4V33nlHrHwQPA5G3zPOOEPYEOAdBJsCDNAHHnggPfzww9S0adPIfcekEBlC3wJsf1FNpdeVJ2jFrDMxY6yrKiOvtuJ9gTpJZ0ehTTREpHtMatj\/djI6dZtibZpft29My40jv05\/ix3ftGlEvXq5VwlCQJ67746jSbGXkQgpQGV011130ZIlSzwn\/Q8++ICuvvpqcVkNTj2beOeobqnyhDSTQuxjo0yBNr+oYST31cH\/byINXA3uyaOzAwjTPr9nsKNw2iY8iYaIdJ0eGxPRnoAMntXL3UzcMqnl2TrWtIl4z45P7NDcTjSDDCw9uIZ+SIQUZJgLGL0QS8jNw0iecgY54MwC4iPpJPXwGvLj4BtcX3FDGhuadRAMl8fWF9VUGp3JXnc1CELQ2QGYtjEov1QxBeWTv+s6PVbSKBCrYNg\/1BS0m9Io1poFiF+YIt3+dmKEd2fDkiVUx80JwBScIuRPhBRwWhm7gJNPPrlcyGpVpjBhLrAjwHOIl4QT0zjz0L17d1Esk0JyIyYPpKAz2UN9o7P6x8SIS5yykFbB\/TWgofCHCrZ67C4E9kIkHYINg08aYy0oTBGIr3Fj7KX0ksQIudOQR6+V7rkSIQXsAjBx\/\/vf\/y6EtXBWj\/MLCJB32GGHGamPQAq4PxmRVZHkVZsoB6SA+xs4MQKMACNgOwLY8UEdZ1tKhBQg5H\/+8x9xfwIOrQ0ePLhMyAncd4zwFp988glNmDBB3IOsm5gUdJHifIwAI2AzArkmBWlDePfddz37oHbt2tSgQQPhfbRx40aRD9uxo48+WtyvoGtTCFIf4WSzLQnnMrB1rFatmi1NCt2OLMuCrX8Sun+omqQ+3Q9Y5MMEMA0eKS5J\/o42xlGerEKWW6gSCnNnGOc9h+ZM6tZpI+osV7\/m6DMZa1GvKvDoknIt3X343Kx\/ZCFe8gCfuy30QIplpwAbAryNZDA6zb4X2UAUJqTAhmYTdOPLmzW9qCq5qT5YFzVs\/YPsD+qJZS8dPCYGGCd1vFvU8kAyXnYN5EOZrpMOZlJHiA3dulGeiS1F1a3rTuA6Yy3IBuDsQ7e68Z2BmYDmzQs+b+J2Ql1HHt0xV4x8sZBCMRqq1qG6pMoTzOySmmwvZG1gO9GIGmPLWZ7c+oedmEFUbrGPTMsLIhqTUaFTN84cmRhcYeNbvbqR0V0zK1eW7gnZUcc1PIxJ1A4\/8gAvmpACzJXz55sTcdbenUySgttAZ1Iwef3N82ZtYDsl1HUn1FEJOVeDcU7MQknhESNJ7ijceg\/9g3NBrrGPDLpbp24Tgh06dJf2XTM6q3+dw8IQF\/bbBQuC77kxcSCTYYp0MFIhz9q7UzRSkBFT4ZmEg2rQs5sMrqBxzaQQhFC037M2sJ3S6qpH4GqKVXwotQwOf3nsAMKir1teEv3jVbcuwe4+Sx3sXaM7geO8l25YIewCpMooLPbqc15hinT6J4m+iUMmrzJiJYWdO3cSjM0yNDZCW+C2NXgbwUf8hRdeKLQD9zfjYNtxxx0Xi3xMCrHA6FlI1ga2myA66hGpgzddDSaLfnDpxewfXYJdvRqH3IKDbOhO4DgzZ7KyD0Ztt2kliDyihikqZt\/oyByUJzZSwKU6uGlt9uzZhToRDRXnCSZOnCjuZ3Z6IEE3idPMTZo0CWpnmd\/VKKmIp4TVHTyabDu8hlUE4t5ADxtn7HwjsGLKnBdZwkz2OqvBmGAOXUyx+yeIYNu27UnTpsUb2werdZ3YfqYg+pFNHGGKit03pvI788dGCs8++yz17dtX6DTh9QB2xBkEqIs2b94svuvTp4+4ZW379u3i3gPc33zrrbeK73WTvG0Np5ix25DJxp2CjW3SxdmZL0+yQDbIc9FFF4lFjDqOwuKT9nNp9I8fwYIUTIy4aeIHAzJ2C0mFKUqjb6LgGQspYJIfNGgQLVu2TISfkKvipUuXCjJo2LCh+B6Xzcj09ddfi1AYuCUNuwmomXSS6pJ65JFHliMFEFOLFi10iko8z9q1a8UVoTa1KazQeZIFGLA8YUeC+3MIj4+oAmo67bSfFm1x1VZSsoZKS+O9SnThwkVlmvfZZ0QOUSI132uswR0fH9tSLKQgD69BQHWChy2hd+\/egiTwPWIVySRX\/J9\/\/rkgjFq1amlho7qj4oFrrrlGhMnAGQn1jmatwjgTI8AIJIbA+vV\/pG3b4iWGevW60\/r1M33bDOJA0iGPypUXUb16lySGgV\/BWCziY1uKlRQOPfTQMpO\/JAvn9wAhLCmoAMoyWrVqRRdffLEghjAH6GzrFG4PI5AHBFaubECXXea\/EpZn6IKMvcADNoVRoxbRq682oCFD3MuVNoAmTdZo1X3PPWvoiCN2k0ixU4XYKSRBCitWrBDugevWraNOnTqV23HA6IyE3QInRoARsAsBnbtmEJkm6GI6pweQzlUFOnVbGGUi9Q7M3E4B6qMFCxYIEpA7Efw7D8bC1EcDN4ARSACBpCdwl6gdBSl06k5A5EwXmTlSANqqS6q0KWS6F7jxjEAFQSDNCdyv7goCv5aYsZKCX5RUr9Y0b97cyNCsJRVnYgQYgUwjwBN4et3HpBAC+yeffFKcxkZSSc3t7mjkMf0+RJMiPRKXPG6HClW34UiN1HzYVBYU63b2xavPNJsRW7a45LGhbwCKiTyyX3AGCknVCmSxf\/zksaV\/gHMspBDbG5CBgmD4lvdCw40WnQkj+IABA2jIkCHivAaSzIN\/44Y53e91XXPjgiouedSrUdOy75jKAqzVu0BkxF2v8Oy2942XPF4HPuMaQ7rlmPbPyy+\/LGJJqfbDbt26UYcOHVzfKdv7x0ue888\/v3CtcFrvjtqHTAq6I9ojn7wJ7sILLxS3ybndHe11p7Tb92kPirDyNG3atPCiFnt34NWFQbJgl4dT9127dhW3BEqHBa+LnGzvGy95vA58Rhz6kR8P6h8n3nhfcOYJh+Sy+O54yaOSnA3vDpNCxKGtDlSvu6NNvsd5izRTFHkQe0omGxwAgmSRWDu92LyufLW9b\/zksa1vME50+wd5VWJDnDO3dyor\/eMmj039w6QQYQaGflRub\/3ujs4KKUSRR30hnYcKI0Ac+lEdWbJEClHkUUG0oW+kbSHo3ZH941R\/2UjaJv3jp86zoX+YFEJOO3KVIweu393RWVAfRZXHbWsMaNM4VKgri2yz207BJvVEVHmcQzztA58m8ripvmxT70WVx7b+YVIIQQoYBG3atClzYM7LOGm7oVlu46PK8\/HHH1txqNCkb6Rh0kkKthiaTfvGSx6bDnya9A\/kxz0suONCNSJntX+85LGpf9BGJgVDUnAG5MPjMvwGzmlI3aD0ZMHvbndK+31v2KRI2eOUJ+1DhWFkkfrdK6+8smBoznLfeMmTdt84MZWD1u\/dUdss848aNUrEOfN6pyK9DIYPm443P3ls6B8pPpOC4UDg7IwAI8AI5BkBJoU89y7LxggwAoyAIQJMCoaAcXZGgBFgBPKMAJOCZu\/u3LmTXnnlFWFM7dixI7Vu3ZoqVaqk+TRnYwQYAUYgGwgwKWj206effiquG+3cuTM98sgjwriMa0Y5MQKMACOQJwQqLCmoB2DkNaFqsC7Ve0h2eGlpKeHe6bfffltc\/FO1atU8jQWWhRFgBBiBiumSKid\/9SY3NVgXfO5xCnn48OGEf2\/fvp2aNWsmfKVxxB7PQ4XUpEkTHkKMACPACOQKgVzuFNTDLbgHdeDAgSTvccYO4bPPPhOd+NZbbxWu98REL\/+Po+aIbIrn6tatS7t27aINGzbQjh076Nhjj6U33niDvv\/+ezr77LNzNRhYGEaAEWAEckkK6Fb1bmc3VZBKAlAfqbFL3K753Lp1q7Al1KxZU5DKzTffXOaUJQ8lRoARYATygEBuSQGdI+86uO+++0jaDWSnmZKCfO6HH36gfffdNw99zzIwAowAI1AOgdySAib9WbNm0VVXXUXvvfdeucBsbqTgVB\/hYhwb4pvzuGUEGAFGoFgI5JIUoP7BFX49evQQOIIAkNTwzk5ScDM0u+0witUxXA8jwAgwAmkgkEtS0AHSSQqSPHD38iGHHEJTp07lXYIOkJyHEWAEcoVAhSWFXPUiC5MeAu3aEa1e7V1\/o0ZE8+al1z6umREwRCCTpKCGmeVVvWGPc\/ZgBEwm+saNg0lh1SoikzKDW8g5GIHEEMgcKfhdZZcYSlxwPhDQnZh1J3qgoptXN18+kGYpMoxAbKQAV83ly5fTli1bXOFYv3698AK65ZZbhK9\/2OR2PV\/Ysvi5nCAQ92RvMoHr5tXNl5MuYTGyi0AspLBq1Sq69dZbCTeP+aXmzZvTlClTIh36ct52dM011wh30zVr1tCcOXMIJ5jr169PBx98sPikmRArCYfeECOppKQkzaZErttmWSofdZSvCqe0QQMqXbGC4s4HUOMus3LHjsGyvPxyuf60uX9MB1+eZIHsXvJgTrBxXohMCggBMWbMGHrsscfo\/PPPp7POOoseffRRMTlfcMEFIiQEJusOHTqIu1arVatmOkY880tVEkJYHHbYYYWrMPEA3FEvv\/zy2OoKUxBiJiE8Rp06dWi\/\/fYLU4TRM\/W6d6eStWs9nymtX5\/Wz5xpVKbMbLMsDc88k0rWrPGWu0ED+uL11ynufKgw7jJ1y3MKW+z+CTWINB\/KkywQ2UueGjVqRNKaaMJpnC0yKWAljF3CXnvtRb\/5zW+oSpUqhHtUV69eTQ888IBYJYMYEEfowQcfpBNOOMG4kX4PwOiMhIvncT8yCMqWnQJI66uvvhI7lsqVK4eSu6RDh0BDZumelWPQqpUaNaJty5eTbpnOfFjxlFnZNGpEsu5Qwvk8pCuLzmpdyq1bpm6+tOtW4YtjrMXdh2HLy5MswGDlypX05Zdf0sknn1xmHsjtTkHGCWrZsmXh1DDOAEyfPl2oiurVqyeY8s477xShJvA3ypYJ6iNcdAOVkRqjCOCDFNziHIUdnGGea9eunSBEmZwTaaNGjWgeXBTj1oOjQl29ddz5COKUlduJnZRbN5+2LCnLrd3OBDBXMd62bRvBbof3LewCJMx4T+KZvMgyf\/58GjZsGOGvTHgPcAaqbdu2SUAXS5mRdwpupICJG1FGIXzTpk1FQ7Gi\/8c\/\/kGTJ0+OvGXq27evOLGMdOaZZ9K0adNI2hrSJoXGjRuXIQW3yRE2mNgnk5QnR125dfNp45Oy3IHkjvZp9PeakhJa8\/e\/E\/Xq5b8zRHl7zj1gwvFbgCAr3sGsJSykvv76a4J6JYsEB9U5FsRDhw71hB790rNnTyu7JjIpYBeAGEHfffedUB9BXYQ7CHAJzT333CNsDJIUMHFHNTSr3kcod+TIkTR27FhRpw07hdgnPd0VZsqTo67cuvlSJYWgXZyc6E1eaZ8y4STRv04dWhxSxWjSDM6bPAJQFy9ZskQYmP0SNAY27hgikwKExqodO4Pu3bvT1VdfLaKIwisIEwBWM998840INQ1XVGlnCNs1IBbsOkAuUEfBVoF6kUAKv\/zlL6lFixZhi4\/8nHP15lagWL0NGxa8IjTJh4riLlO3PFF12VWrl9y6+bRlSULuyKPArADYnSZNmlSwh5k9zbltQmDx4sU0fvx4ocqDGswvYadg404uFlKAsfnee++l2bNn0+mnny5A+etf\/yomCjVhOxXVI0i9RhNlywt0pPeRUM1wYgQyhABUJLAFpK36zBBk1jZVqrF1SAFCwHvTthQLKUConTt3ii0TdJy43B7\/h64fF9MgXXvttUKH5rzXwBQQJgVTxDi\/7QgkSQp4H+H0gb8wcsJVG385JYOAKSlgEWtbf8RGCslAXL7UIPUR7lWGoSetBFuKavxztkN6HwQaFPHiQn0UZHiU+VCRrrpHt0zdfIRm6smtm09blpAdja39pk2bxBmSKN5wIasv8xhUDlAfxblTkGTgZuzEdzgzFGdSbzqU5QbJI+9KR364sauh7dW2qWXLw6ryjJJ0OAmqSy1PfVaNnaYejJX1mGJkSgq53imYghc2PwzN2Il88cUXoggcCnviiSeEUSeXhuYiGz0L\/eII4iaNZuUm0D3qOl0Dsm6+sOND9zmb3B6T8JyD6tbP+yVufTYm7rlz51KfPn1EF+D\/2KEMHjzYVTuA9xiHXG+66SbCJI32gqhq1apVpgvx24gRI8QOBxdewZ6IM0nynnUQCcryet5tPMh72vGsdHFHOyZMmCAu5YI2Q61Td0whnwkpwMgs3NMtS5F2ClARwYiMA2vyikp89+abb9Lvf\/972rx5s3Aru\/TSS+mMM86I5VQvBknv3r1p4cKFAkq5QkjixQrTV6ofvtdEKuweSUz2YRqs+UzQJBp0\/gDVQG7dfJrNCp0tSJ7QBYd4MO6xi10CyDcoxen94kYK+O6cc84pc1+6+s7K9gVN6hMnThQOJCAF9d\/yecwJckJ3kkrQ6l8996SWlzQp2HxeIRQpYMvz+uuvi5UIGHufffahAQMG0GWXXUbPPPOMWB043bEQ8gJsDpfVKMkrIF7cL1aUNspnbZp4osqTJ1mAhU3yuI1dTOx+aki\/\/sQKHfa8oITdgrydMCiv83dMaqou3E19JFUwiImG80pYgTsnYalCUtVHzz\/\/PH344YeE80jYmarqHqeayC9qspNs1BsY5XNLly4tc6GW7Iuo6iOQmLSnumEH7ONW4Zn2oVf+UKTwzjvvCKMxgG3SpAkhQirUOfju1VdfpQMPPFC4qB599NH0wQcfiPAWMELDK6lTp06R2u4VEE9+j4EEl1QbAuJh4sHx9oMOOiiygT0SaDE8nCdZJCnY0jc41AmvPHXCC1L\/xNClkYpw2iWcOwVV7QMbH7wEpf7fbcLFhA0PQgTNfOqpp4TmAdoF\/MVcIu9Lh\/oIZCTVRupvWIi+9tpr4pAsIiw0bNhQ2LrWrVtXkNVZt9cuQ1UxmQAl5yFoSvDOACfYjGQCFlg0Y67MTZgLAI+OmTVrltiyQS+GncPMmTNFCAsY7nCG4Nhjjy0AIVcRcFeFIViqmkzAdsvLAfGiIqj\/fEUJUqaPSHw5sZJG\/DCVFLK4U1BtCupkC7mk+kfuFC688MIyNgdpK6hevXohXAcWmogXpNoL5Gr\/xBNPLBxclSojuMZDZY1F6UsvvUSnnXYajR49upytQt0xyN3EjTfeKBav0q7hdq+7To9LUsB5LBAcEt4dXBtw3HHHlVGh5yYg3rffflswJkG\/hw5AQgeClQ8\/\/PByB9Rw2hkrBazMTMNcqNtS7DLuu+++MqvuvAfE0xmIxciTtyBlNsmDieSKK66IzfsIsXZguwlKcdsUnKtySXKqlxHeYXj8wLCLs0y4Ex0JK3iooB9\/\/HGxsDzggAPoX\/\/6l5hT3n\/\/\/UIEZDkHYEEKjy2ZMAHjrNTf\/vY3scMAqSJqMwgXDigyQU2F79Wdi9Muibxuc00QnvhdksKMGTMEoSF5jbXc7BRkrKNDDz20zATt9b0EBZ3w+eefRw5zYXtAPHXg2KS31hnQfnnyJAvktEmeJOxhIAU1EJuzb+P2Poo6vvA8XIShYsIkjgkTtgWogU866SSt4hEu5C9\/+YswskONNGTIkMhx1rQqVjK59aVNY01HHmObQhqkgJUGmB+RUZFsDojHpKAz7NLPY9OLmgQpAGEvu0QS5xTS71E7WsCkoKhyktopQD2EbaI0ENkeEE8dmvIAUR5OkeZJFvSRTfIkRQqqnPJEM3z8bQzCZseUHr0VXp5k8AjLyjxg9U5BeiTg\/gQk7BSCTjSbnGyMPgT8S0jyZU+67c7y8yQLZLNJHpvaUuxxlbf63Poya\/0bmhRgPG7fvr04o4AEYwrOKMB7QP0ev+3YsYNeeeUVccgtTOhsaUyWpABPJxickZwB8aRLqg2Dbe3atdS\/f3+h7kozcmscWORJFuBhkzyyLTYtaOIYMxWxDKdrvN9Yg3tqmiF5vPonNCnAqm+a4CGQFCm0bt1aTMCqT7Bp+zg\/I5AmAkwKaaIfT93Oc1R+pWKxiI9tyZgU4HMLVzH8NU2IU3T88ce7hrvwcz117hTc7lPAYRV4H+DDiRHIEgIyBn+cpIAbIPdoXV2haNOGKM4bIcMExEPDVHdV\/N8vMF6x+zTMATZJCvKueL8252anUOyOQX0qKXgZmp0xT9JoJ9fJCIRBIAmdMwLm+twGKX6LM1CqaUA8SQiqV6FfyIowuEZ9JgopxEnwUeUwfd54p2BaQRz5VVJAeeoWLcvgx4ENl5F9BPJKCn4B8RALySSyqVe4a2fvS29FfC89FtXdiFRh4\/re5557TpzlQBgMtx2KSgpqufJgGw7fgdTkITrMRUg2RGuO8lZkghSiCMjPMgK2I+Duxhh8W6ufXNOnE\/nFxMOd8T16hEcG13iod\/WYBsRzhrkIulvBLTSFM9S2Go5bXuYFMgEByJATMpwGJJcOK16huyUpdOjQQRysQ5w3NabTP\/\/5z8L5KRkq\/OyzzxZxjbK8WGVSCPFeuK08oL5SdzDqysP0+xBNivRIXPKoqyn18pJIjTN82FQWFO+mtvDqM8PmaGV3I4Ug9Y9WwQlmcqqfTAPi4d4Ct52CnIgRJw2B8SpVqiSC3yHAnBrWAuPrd7\/7nfAiQ9BNxDnCB1cCq6Ez4LHoJCy5okeEZ3mxjwzJjXA9aMO5554rxoXMo44HObZBCgjiB3umjPWEmEzYoeA0NcJuIMm5wHkxkBqcz4Z3Rw4XJgXDFwcDbOTIkTR27FhxIQg6E9tPxG3BQEA0RySZB\/+WkRx1vi+2bSQuebAyg3tw9+7dxUuSRjKVBVjLQ5fwppOru2LbrfKyUzAJiIfJ2hmpAGNGRkGFCzf656OPPhJu7\/iLpN7Ohn56+eWX6de\/\/jX96U9\/ol\/84hf0s5\/9TOSTK3qE7Ic7PGItYYz67RSuvPJKEcMNwfswlrA7wP8x0UvZdHcK2JkgGJ\/6zqOt0oYix123bt0KsZjSfHfU95VJIeLshRca21Bsh8eNGydcbtUJUg503e\/TmlAlDGHlgY5YDWMcEdZYHg+SBS8uAqt17dqV+vXrJw5HAn+vA5JJ9U1ebAqmAfHwnnh5H+E2RQTFQ8BNTPQ\/\/\/nPywSxQ9\/BwweTLM5FIU4SQmVDpSTd5bESBxk4w3Yj9DZW+KpNAQuC2rVrC+9FxI1C4E\/sQLBTgPoIhIFysUvAifDzzjtP\/CZtCnL3sHHjRmFTgKoJsvktliQBonyb3h0mhYjTj+xYDDKvQ3Um33vdUxuxmdqPR5FHjUYZ9pIS7YZqZAySRWItV20qKbj1WVJ9kwdS0OgOoywIzV+5cmWxQ8C9LG73sKDfcEcLwvT\/+OOPgtx179t28yzCLkBGad2yZUshMJ9Xw2XEBXWxIPuyVatW4rpRJHnA1rnLkUQgiUTWk\/a7w6RgNFTLZla3v3JV6nbSOiukEEUedcCr91wkNZEGdZuOLHkmhWKfUwjqjyR+ByksW7aMzjrrLOPiw7ibOisJSwp+rrc2vDtMCsbDafcD6g1Q+L+XyiEr6qOo8jhVK0434pAwh3pMVxbZZredgtcByVANCngoiZ1CEu20rUwcoMUd8XXr1rWmabIvvdRHXtcJqwKk+e6gHUwKIYaTNFapE6GXcRLF22xolqSF6JlR5IHbHwIXQgXjnGRDQBz6EZO+kUZ9Z3ttMDSHBoAfTBUBSQpuhmY0zMvjyoZ3RwLHpGA4hNxim8jDLOotT6qfstdhOxsO4cUpj9vBIUN4I2UPIwsqdCOxYvYN7xQidbtVD8u+VF1S5Vygvh+y0dJdNe13RwWRScGqIcWNqYgIyIkkD9F0K2L\/qTLnIeItk0JFH8Usf+oIwA2SI\/ym3g2xNQBnLOAua2NYbB0hmRR0UOI8jEDCCHCE34QBLmLxtkY\/1YWASUEXKc7HCDACjEAFQIBJoQJ0MovICDACjIAuAkwKukhxPkaAEWAEKgACTAoVoJNZREaAEWAEdBFgUtBFivMxAowAI1ABEGBSqACdzCIyAowAI6CLAJOCLlKcjxFgBBiBCoAAk0IF6GQWkRFgBBgBXQSYFHSR4nyMACPACFQABJgUKkAns4iMACPACOgiwKSgixTnYwQYAUagAiDApFABOplFZAQYAUZAFwEmBV2kOB8jwAgwAhUAASaFCtDJLCIjwAgwAroI\/D9uaVd8kaRUZgAAAABJRU5ErkJggg==","height":234,"width":389}}
%---
%[output:24228ba8]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAYUAAADqCAYAAABN9llzAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQm4FcWVx0tEZTTGFRSBAGpU3MYNjbiBEQnuuIsbGjRqRMVt4o4iZgJuqAQVI5uKccQFdxABd1Q0RhJUUDCgqOi4JBqijE5+BedZr+nuqu7bfW\/3fVXfd78H91Z31\/lXdZ06+wrff\/\/996oK7R\/\/+If65ptv1Morr6x+9KMfVeGJ\/hEeAY+AR8AjkBSBFarFFJIOzPf3CHgEPAIegeojkCtTQDKYNWuW+vvf\/x5K2cKFC9Ubb7yh+vfvr9Zaa63qU++f6BHwCHgEPAKNEMiNKcydO1edffbZ6vXXX4+F\/D\/\/8z\/VH\/7wB7X22mv7qfEIeAQ8Ah6BGiOQC1PATDFkyBB12223qQMOOEDtueeeasSIEapt27aqV69e6tlnn1Xjx49X3bt3V5dddplaffXVawyDf7xHwCPgEfAIgEAuTAGjMlJCs2bN1NVXX61WW2019dvf\/lbNmzdPXXvttdrQDGP4zW9+o2688Ua17bbb+tnwCHgEPAIegQIgkAtT+N\/\/\/V\/1y1\/+Uv3sZz9T\/\/Vf\/6XJ\/OMf\/6hGjx6tVUWtW7dW\/\/rXv9Qll1yi\/uM\/\/kP\/bd68eQHg8EPwCHgEPAJNG4GqMYUXX3xRnXPOOWrkyJFqk0020aj\/7ne\/U6+88oq69dZbvaG5aa\/DTKhnjfXu3bvRvVBhor4028yZM9UJJ5ygPv3004avf\/WrXzUcYMz7IOEeccQRica3YMECddddd6knnnhCYVujtWnTRu28887q6KOPVltttZWWon3zCBQRgVyYAlLAhRdeqL766iutPkJd9Pbbb+sXceDAgQ0vKUyBF9Abmou4NMo3pjCmcMopp6jzzjtPrbDCCg0EIbVecMEFjQjMgiksWbJES8ODBw9W3377bSSA2NUuvvhifxAq3xJrEiPOhSmA3EMPPaQlg6OOOkqdfPLJOmiNF69jx47q8ssvV1988YU666yz9IshdoYmgbgnMjcEwpjCrrvuqm666Sb14x\/\/WD8XN+lLL71U3XPPPZkyBZwr7r33XnXRRRcpmIOtwRh4D3wgpw0p\/3u1EciNKWBsvvLKK\/XLx4s5dOhQNWHCBP0imG3AgAHquOOOqzbd\/nl1iIDJFLBb4eb84YcfapXllltuqSkmNgZ716JFizSjEPVOpZLC7NmztST8wQcf6OccdNBB6vTTT1cdOnTQqqLPPvtMPfLII9orT+J20qim6nDaPEkFQyA3pgCd3333nXr55Ze11xEvCf8fNWqUuvnmmzUMiPZ9+vTRxmbfPAKVImAyBeJf0N3fcccd2vNN7AL04RCCEwQqpWeeeUY\/tlKmgDSCxCvrGu+7MOeJadOmaWaBajUoxVRKv7\/eI5AFArkxhffff1+tueaa2h3VN49ANRAIMoUjjzxS2w4OP\/xwdcUVV2gVpmzeffv2VV9++WWDGqkSpsAGj3s1kkDLli31wadTp06hJCNBn3vuuWrixInaC890vKgGRv4ZHgEbArkwBVn4MISrrrpKrbLKKrZx+N89AhUjEGQKGJj5rLHGGtqZgSBJTvBPPvmkDqxEir3lllsqlhTEBZvofZfT\/zXXXKOGDRumn4uXElKLbx6BoiCQC1MIi1MoCsF+HPWLgMkU9t13Xy0lcILn+zFjxqh11llH6\/0xCnNCv\/\/++zNnCvvvv7\/67\/\/+71iVqOn95JlC\/a7HslKWC1PA+wJ30\/fee08b1jiplb1hD+GEiU64Z8+eapdddmnk5lh2+uph\/CZTYHMeNGiQXn9jx47VEsLmm2+uUBvJaX748OGZM4Xtt9\/eGnfDuyESimcK9bDy6ouGXJgCGyheGCz+OXPmqB49eqitt95a63SDDdUSvxVdxfTuu++qP\/3pT9pgjqGcTaddu3b1tRpKTk2QKXBif\/DBB3XMDAFsrVq1Unfffbc69thjdRQ98yjG4UpsCpLWhUODi01BVFjeplDyBVenw8+FKZg6VhtutcqSygYybty4RqK+TaxHApoxY4Z66aWXtBrC+5jbZre6v4cxhTfffFN7GxFQueKKK+q\/2LkwQpvzXQlTQB2FAfu6667TBEfFIARjGVzsD9VF0D\/NI5BTQjxevD\/\/+c\/6BbS1WkgKshmY+l\/8zNksMAISfQ3DwGOFf0PHZpttpv3eSY3A9aiQNtpoIxt5\/vcqImAyBdnkJS4B5kDD+QH7AkkYs2IK3Je6IbhXE\/9A22uvvXRwJildYEasG6QUEkBKtLOPU6ji4vCPckYgE0nhn\/\/8p\/YHR53yi1\/8wvnheXVEUiGaGrUB6boxNnbp0kX7qrNxYOugPf\/88w2SAhuE\/B96uJ7rUDlwwuNl52UmCIoMr19\/\/bXae++98yLB3zcFAmFMwXQX5ZYwd0nK6MIU4oaxwQYbaIP1T3\/6U71GfERziknzlxQOgUyYQtDbqAjeR2aEaZgxz2QCBM\/xf4LsyOoq4+ff4i6I3hgdNGk5YCqcAn1hoGKt5zCmwAjNwDK8krA1IDEQV9CvXz9NRJT6yJUp0A\/1IusICdPnPirW2vCjcUcgU6ZABCl5ZQgKCqbOdh9Sdj0xdGPwDnMRTMoUZFTkzgkzmGc3an+ntAhEMYWnnnpKex3RMPISUUyL6h+WQylsTKakYP5OllSkhsmTJ6u33npLMwvJknrIIYeozp07+yypaSfZX5c7ApkwBTZK8hyhQmLxk1MGryP+or6Ja\/zOySrrGs1s+uhwTzrpJF0HWuo6yFjCmEJQfYT6CdWAbx4Bj4BHoKkgEMsUUJmQ0I4TUZytgCCgxx9\/XOtVp06d6pQlUgDOw\/sI9Q9ZWo8\/\/nj9GBgAzcyLH2QKYYZmWxBSU1kknk6PgEeg6SAQyxRcbAMwArwoSBkgdRFcrqs1xEGmIMyDKNgotUCtx+yf7xHwCHgE8kagEVPA9ZIU1y+88IJ+LmohmxpI+uyxxx4NdRFcJYy0xJkRoX4DT4uiv84j4BHwCCyPwHKSAoYxjMSSF94FNCIzSSeA22feDXdRXEUp3uMTieWNtr9\/k0Vg3rylpHfoYIcgSV\/73XyPGiNQsfqo2uM3YxC8Ebja6PvnlR4B2wY+dapSFMLirzQYw8iRSnXt2pj8JH1LD1zTISCWKUhkMl5Em266aSFQCboLmv7luALyoa2\/\/vr645tHwCOgVPNnn1XNBw1abrNfMmKEWrLrrhqiFr\/7nVIDBkTCpfsec0zivh7\/cAQowhRWiKnWeGXiklorIkSVhNqKrKXkzp8+fboeDp5HtS7ziX86kc+rrrpqISc\/ybzVEy3QXW\/0NF+wQK+1ZhtuuNxaW3PoULXW0KGR071oyBC1pE0b1bp3b+uSWHjXXbqPa9\/FKWpF1NvcRNFDEbKsXfGtE+jQwYkpoLLB1VTSQ4Tdl\/TYVLiqdpI4jM40DN29e\/fWtg1iJYogKcC0PvroIz2WFi1aOExHcbvUEy2gXDf0oMIZMEC1WHYY0gyvbVvFqV6re6ZOVS169rQuLK6BsdjakqOPVuq997TkYWv0XXLbbT90s6mulvWsm7mx0FNaSYGEcKeeempDgfOohZBHvEHYs1AfUdMgmI6CvjCFIuWnX7x4sS4UjyG+7EyhnmhhrdQFPej+Y9Q92g4wenRjlZFtJ8\/69++\/X\/p8VztFvcyNgWPZ1lqspCDFckaPHq1OPvlkdeihh+oo5bBGEXSkhWbNmmW9rJa7n+mSKjYFsTV4ppAP\/GVb2DYUSk8PG223bjYya\/97nz5KjRoVPQ4YF32MtmTOHH2Yatm5s\/UwRb4yWgcXL6kaoVG2tRbLFD777DPNDDp27KgGDhyYaSEcCvF88cUXepqyYCaeKeS74su2sG1olIGe2A0PhmB6CNkILvLvU6Y0qLpcJQrU2ZdffrlWa0uDMZC1tmvQS0opneyyVsyjDGvNXB41c0mVqGcGI5HQlaxbzxQqQc9+bdkWto2iItNj3fDY4Dp2tJGYz++y4WbJkJAUOOnbVGHLJAqYwYCYvjAGalvQrFjmg1KjuxZ5rYWRH8sUiFYm6ymn+qwlBaQQ7k3LIiGeZwr5ru6yLWwbGrWkJ+7U6rThsTFnzRTYlMUQHAWexCvwey1UV1OmKOSCbg7PnjJlirY9ujIP23qp5PdarrU047Z6H7HZnnvuuboI+u67716IYvVmrIJUr\/JMIc30u19TtoVto6wW9NhOrfzuuuF1ddgYbRgYepelwWkwhRNOCL8MhsDp+7LLlv6OnSCuL4wrzpbgPDijY58+qtu8eY1URtH8q0ODyijuUTCPMHVTmuFFXVOLtVbJ+K3BaxSrZ7Gi4qHq2E477RRaT6BaLqlmRDOEmyU0vfdRJUsh\/tqyLWwbEnnQE+dx6SIB4NBh6sijaEA1MpKHuahw2JzJFuy62Ud5CsEMAgbhSK8i6bvCCrZpSPx71nfUWMIQc2x5rLUch6ucbAqvv\/66dQymS6qohrbYYgt1yimnRF5LvVo49W233eZcxQyJAO8jmBQV0yQPEg+BKVA32VbDwUpMhh3+\/ve\/aw+KlVZaKcO71uZW9UQLCGZJT3CP5mBt2jtHZXxq7rMsBsG6EugnqqEoJhLc7OWmEOXq1RPWN2OatYBiJTh5B7E\/JL\/S\/YqwtUZsVTWe7T7KpT1jmYJ4CJEe29ZMl1SX1Nm4u2JLoABOEkMzTGHcuHG6mhpN6i+3b99eM4W5c+fahup\/9wh4BDwCNUcAtRWH4qI1q03BZcBs8Ndee62uYZy07b\/\/\/qHlMqPu45lCUoR9f4+AR6CICNQFU0BiQAz6v\/\/7Px3EhiQhahGCTS677DK1aNEipzoMTBJqHmrmEgfh2mzqIx+85opksn5l04vaqHOlJ85O4B4qgNIjwoBrG2jI72GbCfQsevllp4CvFI+s7BKbncIlEG+Z55OL95EEsomXV9zgq7Exu661ykDO7monSYFcJKhsUPOw+WM\/uOmmm7R0gIG5X79+jRI7uaiP0pLgDc1pkavsurItbBu1c+Ys0Wu5c+eWoVGztswMyUMFsjGRRgVolWZ+ouwUNm8mw\/MJ+8wJEYZz8EFPj77e5skVF+xmWz9Jfi\/N3CwjysoUqKKGBEAdZhLNQSAnfBjCJZdcoih2v+OOO6obbrhBtWrVSt9WbBFIEXkkyDNdUkUy8C6pSZZp8r5lW9hRFNo2e65zSSmUPFQAaXhZ4ZqIwXFqJbuvbcPjfQy2upgfm0RhEB3l3gs2Yrx1YR5hWCZ\/O+KvKNvcWJnCmDFj1JVXXqmD1w477DCdhZQNGKlhtdVWU9Q65ndiGUiJkbTBQDBS86mkeaZQCXr2a8u2sMMoctns0VK4hABgH3Tp98M44te3eWp12fDqkiksIyqpKgw1UVTuozRY2t+GZD3K9u7EMoWvvvpKe\/fg+ikRzbiDClNYe+21FYV4kBg+\/\/zzhhrNAtn777+vJk+erD799NNQFCmIQ2pppAzuVUnzTKES9OzXlm1hBylyUVtzjUtgL\/3QZiQLFbCrPMJOrXEbnklj2eenGrS4Yml\/G5L1KNvcOMUpUAuZVNW0IFOI+u61117TkkMUQxBYd911VzWUIiBrreWMtJkldYMNNtDBJzzHB685Q5i4Y9kWdpBAd6OwOzQu0oJZyTLPU2vZ56caTMF9ZrPtWba5ccqSuu2226oLLrhAq3iCTAGPJFJNwARuvfVWvbnjoopk8cADD+hqaFRFQ8WELYLaDDNmzFDDhg1T1Fgm3gBJxLVJtbWjjjpKwaykFVFS4GRClCp64iKn9nXBvsy0JDcKuyCiFCExSCCuwcLmXbM+tZZ5foJo1xMt0FY2eqz1FNjcX3jhBb3hs7EFmYIU4eHEjxqJakJffvmlOv3003XFMVE7XXPNNYrIaLyWcGd96623tCSByLznnnu6vYVKKdP7CKZSZKZQREblDHSgY5lpyYspSExnAvtoWvit15V5foLE1RMt0FY2eqyGZiQArPkbbrihOuecc7SNgHxIeB+xscMk2Kix9CNR0MJcUvFeQjpA1dOuXbsGaQKbBJHNK6+8snXhmwBL52CRnTPPPFPnZypCw6aCpFSkMaXFpey07LzzD1JlWgzM6zbeeIEaO3b58pXvvadU+\/ZZPCHZPco+Pya19UQLdEXRg+akSCl5ZA6sTAH10MMPP6wuvvhiHbgWbLidoj7q1atXgweRFOfZYYcdGmwRcEviGWAKJNajhdknkrwKokrq0qWLVlGxAU83atUmuZfvW98ILFx4l1q82M4YmjdfoJYsaRsLBn1atjxPtWgxvb5B89TligCHRT5Fa1amIAMmUhnmMHHiRB2rQJK33XbbTR1wwAHLcTupw\/C3v\/1NSwfYGWbPnq39r7FN7LvvvkpsES+\/\/HJs7iO57oMPPlBhKTFgLDQM4Xgz8fHNIxBEYM6cturYY+M3e4zCAwcuoC69uvji8L6SQbpHjxc9yB6BihAoraSQlmoKXJx00klapYR9oVOnTjqlBQwDyQImA4NAakhibEbi4N4wAVFT8W\/T6Jx2zP66+kYgQdCsNSt0fSPlqWvKCDhLCklBwgMJz5vBgwfrkHNsEE8\/\/bQ666yztD2BhlH6+uuvV\/vss0+i25suqWJTSHQD37nJIpDGKJwkg3STBdYTXjcIWJkCG\/iECRO0GiguLbVZT8FEB4kAlQ61FWACzz33nLr99tt1lxNPPFHbApo1a1Y3gHpCyoOA3+zLM1d+pNVDwMoUHn300YbTPdlMcScNa+jHXGotw2QwWJNIr6zMgNQeqL5oJjMMKxNKn6TfV2\/6lz4pK3rCggpNt+Fq0JWUFsYUFvsSNWfVoMF8Rlb0FGFukq41mZeHHnpIQ2JqBco4P3H0FGV+wDmWKQgRM2fO1PEF2AWSNAzN2AtWWWUVXTaTIDVJnTFnzhzNRDp37pzkljXvi+FbSoCSmoPJxAh+\/vnnaw+tCy+8UI9R+vBvXHldv6803UdSgLKix6yCVyv7TlJawFrsUsTQSHLFqEy8RZ+bKHqiAj6TrpVK+yedn0mTJunAL9N+eOSRR6ru3buHvlNFn58oenDWkQqStXp3zLl1SnPBxi0Rza4LQ4La5s+fr72NqFGLlEHWVYrxjB07Vq277rpq+PDhapNNNnG9beH6SdGfgw8+WNtHwsqERpUPDfu+1osiLT3MoTC\/aksHUYvCRgtSHnm3Dj30UJ3QURwWomp2FH1uouiJCvis9ctkm58g3rwvBNBSZbGM704UPSaTK8K7E8sUJDIZe4DkPnJZSKiImDRsEeQ1IpgsmAVVmAbRzNwbe0MZm7lQo8qEJvn+iCOOqCkMldBD7ilpRXAAsNEiWAe92KKq+xV9buLoKdrcsE5c54e+JmMjz1nYO1WW+Qmjp0jzE8sUiCVAbUSaC4k3cNmxJHjNzJkUvC4sZ5LLvYvUB32viLdxZULLwhQqocd8Ic2gwlq9qC60lIkpVEKP+c4UYW4YTxJ6guqvIjLtSugp2vxYDc0ff\/yxFq1JTXHaaaep1q1bWw3ErpXXKo1oriWDkFOObCxxZULLoD6qlJ4w0Zj5SSJhZjWfrrTImMMkhSKpJyqlJ4irGfCZFeZJ7pOEnjDVV9HUe5XSU7T5sWZJvfTSS9Vf\/\/rXWHdUiDK9cISz85dEeKuvvvpyawbbAsyGNBkYoynYU5bGIiD2wtwIo4yT0FRkQ7OI8ZXSgzqwCEGFSeZGDJNBplAUQ3PSuYmip0gBn0nmB\/qxRZI00zQil3V+ougp0vwwRiem4JI6IuiSKhXb+vfvrxPqmemxYRYk0Lvuuuu0x85xxx1XFn7QyL1UBi3pN\/BgEd2geLLQJ6x8aNz31QTDHFul9NQ6qDANLaLf\/eUvf9lgaC7z3ETRU+u5CWLqstbMMUt\/siEgnUe9U0V+d+LoKcL8CHZW9VFakM3azkgKeA0Ql0D5TfTwxCqQRI+TQB51nNOO21\/nEfAIeASaMgK5MQVAJSbh7rvv1hHMCxcubMAZuwTRzPgcl0lt1JQXiqfdI+ARaBoIODMFTvakqEBFQrAWJ3+C2UhTEWYzCMKH5EAyPOomeMmgaSwuT6VHwCNQPgScmMKzzz6rsA2E1VteZ511dNbTn\/\/858vFIpQPDj9ij4BHwCPQtBGwMgUJMiOQjfrKBJsRmcz\/qa9AJlTSWBDJu+mmmzZCE\/URqS6I0iM4DWmBIjtPPfVUQ73mzTffvGnPgKfeI+AR8AgUCAFrjWas4o899ljopg8dlOTEc6Nnz56NIpNJb0FU60YbbaRdTldddVWd3mLIkCEN5JNgr+xpLgo0l34oHgGPgEegYgSc0lyweV9yySWhqShIaTFw4EAdx0D0M1KERELfdtttOlEcuWWwSZx88slqxRVX1DmCJIVGjx49dBnNYBqMiinzN\/AIeAQ8Ah6BxAg4JcQjSCsuMjUYmYyaiCpra665pmYYqJfwKyYeoW\/fvpoJ0LBFvPnmmw3MJPHo\/QUeAY+AR8AjkCkCsUwhbHMPPv1f\/\/qXliI+\/\/xzXV0Nz6KwNBe4pV555ZUK6QG7BK3MaS4ynQV\/M4+AR8AjUBAEnBLisZGTYnj33XdvpOZBTUSJzTPOOENLANRiRg0kCfF22GEHLWFI2osZM2ZoQzOG5zC1U0Ew8cPwCHgEPAJNFgGr99H777+vbQEUyED\/v\/fee6u11lpLb\/wTJ05UTzzxhI5ZuPXWW\/VfGvEI5Ex655131I033qg++eQTbYxGDYXRmZQXr732murXr99yBuomOxOecI+AR8AjUAAErOqje+65RxuHKcv58ssvLzdkCvCgDiJp3oABAzTDoGFDwPsIAzMNl1QMzLvuuqtOcPXII4+oVq1aaY+kLbfcsgBQ+CF4BDwCHgGPgLOhGeMwUsN7773XgBoVkDbYYANtG4BhEKsg2QxRLb3yyivq97\/\/vVq8eLE2MiNpYIOg9NyiRYt0NbetttrKz4JHwCPgEfAIFASBRkyBDZtKaRTVETUQtZRxMyULalhDVUQfUi+LodlGm6S7sPXzv3sEPAIeAY9AdRFYTlKQYDTyG7k2EtwRlNalSxd9CUbkb7\/9tlG6bPNeZrH0HXfcUdsdWrZs2ehxZh9JTY0twkwxi5QihmvXsTaZfvPmLSV1mZ0nkm7Xfk0GOE+oR6BpI2BVH5HNlGhkVEFs9sFGMjxcTffbb78GzyQ2dDyWTjnlFC1lBBsbO99zz1dffVVv9BiwzSbVjA444ACtbjrqqKN0IR\/5d62LqBd22UydSmUSpfgrDcYwcqRSXbv+8J1rP7nCM4\/CTrkfmEcgSwRimQLqpFtuuUWf5LfbbjvtLbT99ttrCYAYhj\/96U+6stq7776r+wXLG2644YbL1Uv48MMPdZEM4howQrPRc29iHaQFK2FJ\/dOTTjqpoYoZbq3BRjEgKQi0\/vrrKz711JovWKDJWRKhymvxu9+pf1v7I0leMmKEWnLMMcq1Hzdq\/uyzqvmgQcsxGX2vXXcNfZZtnPU0J54Wj0BaBHC+4VO0FssUJL4AVRKMIWyTFZdVDMZXXHGFTo2Nkfn+++\/XhmQkCCmkQ3I8XFUp24jqBwnjo48+0jYMcVUFoGBdVpjC888\/rw4++GB1wgknNGCId5NEWsMMMIZPnz5d\/3788cfXvKIbktXXX3+tJa1KJr\/Fiy+qtYYOVS2W0SaMYdHgwWrxz36m6aVP6969revrszPP1PeytYV33aWfF9d30ZAh6h+HHNJwK5dx2p5brd+zmptqjdf2nHqip55o0e9qxD5Axgfx1rTNbzV\/t6qPiC\/gNL\/11ls38jwyB0mdBSSAcePGqXbt2umfTMaAaohANnIjISGst956WsLYaaedlGz4LkzB7CMMCzuGWZ4P20abNm00A6u1pMAYYXqMo0WLFunmdcCApSf7iLZ4+HCl+vRRzbt316f6rBrSiJz44+65+LHHlqqlHMfZ6F6uKinXfgmIz2RuEjwv7671RE890cK8R9FTSkmBpHWczFEPffHFF7HrmtPwhAkTFCojaTAG0mtThxlV0cYbb6wuuugiHa\/ACR91k6iGzNxKUeqjYP4l7A40vpearWZt5LxfRNv9ccWl4hyG+FimELXpoffv1s32mKX2AkOCsl+QYY8+fRDL3MY5ZcpSBuJqz3Dtl4Ic57lJce9aXFJTelyZtmO\/mtKSw+SVjZ5YSQEPoqOPPlobhCmdSWRz0HD88ccfa+Pvbrvtps466ywd6GY2GMOkSZO0amevvfbSqqRhw4bp6GfTiBw0HIcZmrkvqieYQJBxlJIp2DY9GIJpMI5asLLR5rCgnW7p+nwYCEbvGLuHZnD0w1ju0s9pgMt3KtuLaiOzJvTY1q8M2rXfsv64uHOYIjDWKmE7Mhobfnn+vmQZPS1d6MlzII73jmUKpLJAfYRRGQbBZozKB10YBXQwNGNrQErApsDCNBsTSi6k7777TlG9jVKev\/jFL3TMw0MPPaRVKzvvvLM2UpNID0ZAvAMMwnRJNW0Hpkuq+X3pmIJt02NDjNsUHSe4lN1caRfJIwWRNdlEU4zT9ZKs6Zm3bLOV1DXLjcO2flMw96lTp+pDI3+l8XzczruannP8mJDRWOlZ9kDXfk7zknCMTvesQqdYpkDuIk7zGEuT1jtAQqDZrsNeYUZCp6W5VEzBVS2UBAxO4HJqSnJdmfsiUbD5pGhZb6IphpDpJa702DY9p43Zdf0mYO6XT5um0+RENRhDH+ab5sqQNO9wYzSu\/czxxWKZYIyZLoQMbuZkaCYCGU4d5UEDA8B9VRgB40LNdOCBB6rVVlstdpgwjTXWWEM1a9asInJKxRRc1UJJEOHUbLM\/wDh4sWwSiAS8uTAZV9VRElqS9F12+EhyCX1dN9Gk961Vfxs9Lpsep3SnjTnj9Tu1QwfVzWGtTZkyRelIG9s6p8+UKcqV0TjTvWxyrVi6Ms0KJN0815k1Syqb7bnnnqsGDRq0XOrsNAPDPQujc1pGIJs\/z6ZID55HtCKSgX96AAAgAElEQVQyhTlzlizTjbb8QTfK4u\/YMQ100dewMbPARo2KNjgLQ7jsMrd+e+xhf\/kkKC4YLJctdfF3mzvXHrUdcgfbJlpNErJ4VuhaW3Zjl00PNU03h812CqqcHJwaVnAAAUlhJO+Pg53NldHABOMYoQxLM6SuXbV6y8o4R492GqM+oKWUdB3gSt0llinw4pDimpoJqHjIZopNgViEYEOK6N27t1pnnXUafiIuATdSKq9dddVVOugNWwSGaYxJxDVgTHJtZvwC13BPXFtJwlckphCrSuyQMVMIRitHPRxmIOI34Ln0c2UyLiejvNRbTVxSsKmtOdW6bPYwBVGHxL2Pfbp2VSMdNmXXd1r6cUxalpgl9tKlSmm35sJo3O7Eq9NHxz65YDnl36+XkTsg\/hEp16\/ruNP0s9oUOImzGdsaaqIHHnhAbbTRRrrr22+\/rU499VQ1f\/58te+++2oOi0oJozXpsseOHavWXXddNXz4cLXJJpvYbq9\/Z+PH0AyDgsGYKS+EKZCJFcZVq0YS2TjtDL+1H\/BDAJ51nDb1DJstG35YYzDt21sfoeL6\/VvXqyWQYGNcuKJKs52Okqik7CNe2iOOdss9kFiJmcHBwerh4jqeKvdzWWtTp45uZLjNYojprDjxT3Z9I7A6OKxo\/bB\/+685MRpXTJAUTCN41HUwBOPNiLw91We6ppR0Xcecpl8sU6A4DhwS97D+\/ftrzyFT9487KplRcTklVQWFdpAUeOHYvIlbIOsqm3TQ4CxMg1oMeDW5RPyy8RMgh\/RBgykEg9fmArJvHgGPgEeg4AjAPKaUTVLAcIzbKZs2ouUxxxyjYxE4XT355JM6EI2gNpLVUYpTMp1KOc5tt91Wp7oI80Di3tgEUE\/BTFzCvT1TKPgq98PzCHgEnBFA8sBWUbRmNTQzYDZ5opA5pXfs2FGnssDOQP6iwYMHa7uAufFLjAHxBsEoZBMApAk2eleXVBf1EXEVpLmQJjWBglqUKK1IcIJQw6NCjXOOQJOyLNOD8\/yOHBCjZxKjsIvqx\/mJlXXEvoTqjxgVF6ku9mlMCqAGM7kG1WAO\/ZLMI7ZzaRhmoWfLLX8US0\/U+jHpixoDa0eeiS3etoZEG2e7n01ttDz2rsoZ+xphI0O3rlWOUXpSc\/069KNsl83Yi80DrYVWHTl4z43GbObg0WSneGkPV9WR3M+qYluGUZ8o1a\/rwHLo58QUyN1B+ghyFyEZcMqnrCYbBNIB2VPNk77kJeIvhmDSawcbLyReTSuttJJWB9lcV7k+iaHZZoBz9apztY8mzTShVYnzItJcB43COUx80lvm5q3Di2ur+cBgI\/q5zqM4etjWheDi2s\/FHR3yHBx7tAMZDMEWyM0BJJkDm92M62JoXi6QzMVZAUAd+o0aNapRsktzfQpDoIyvbg4OEFP32MNqFJb7ujAk4iSCgXVR75CWAGCcUV5apidg0hexCv2tTIHNGzUREcgSh4CRl00CwzH5kSiUc8MNN+iay9LGjBmjs6Bii4DDc400mAWL4LrrrtN5kSjV6dpMl1Qzz5HpffTEEz+zGntthw3X8Ug\/Np0we2zUfZZTJbpujkkHllH\/3JhCBeNL6t1ri6VKGoTr4nQFea4HC9d+LiEpjWGN98ORzR4VsZmFOHZjDs6b6\/qN6Rfl\/w8zaAhcs3Ft40Dlymhc+7l4cjVinA7MsILln9ulsUwBJsApnoI5pKpATYT66M0339SBaaiUiHamsfn\/+te\/bhgozITJJIU2kgJgYaTmPiw+YhV69eq1XL2FtJQKU\/jVr+5Sv\/nN0nTS1W42RyEZj4QVVHt8lTyvHpiCC\/02xiH3YGOuVXgGBxBHd32t1jz+ePdTeKKN2QXQFH1Ya9R8d8p9xP0zYDSudLsyEJNs6Fn08suqLnIfUaOAUphIA6iIyDWERCB2ALyTkCJImIcaCW8jM101MQl33323uv3227UHkzSyhlLRjSR7Lmojl3UlTGHhwrvU4sW1YQouKqSwImgu9NW6TxGZApiskKUzegKQXQ8ACW6ZqKuLtGCuNddNzxwEh7fI3EeJRpusc15rzZUeW7+kWOZFTzJU3Xtb4xQOOeQQhdsotgGMyUHjMKogRLu\/\/OUv6sEHH2yIUwgOAcmBdBkEvvEXgzAJ8sz6y+Y1SWs0wxQOP\/x8NX\/+0+7UZ9wTlZCDujMyrCDj4WR6u6IubFebQqZgFOBm2KTQTqRRW9s2vVqTV9S1FoaLC5ZlokcftL43ExYFqP7000+1dEAMgriWBpmCuJZSPQ3XUqQAWwtLi+2SOjuuRnMeTCFJvJWpEiqpKjF22lj8o0eP1p4ntTg9Rg3OVa9vW5Nl+11sUn6tFX\/mivruRCEXyxQIQhs4cKAul8mG3759e\/1\/9H2ohPA+IliMyGW8kaizHHRXFAlBBoB7K7EOqKPYXEibzX0qrdEs6qO5c991XiVx6h5xEEiSAiiY3VcY1T33DG6oX+08uIJ1LFIakSA0LtJZ1o4FjMFVheRiQE57ABEs\/For2AtjDKfI704YalbvI4LLOB3iaYSqCLUO9ROILcDQPGvWLO2VxCmSYDVz88f7CEOz2UzBBHUUMQXkVCIyWjyU0tRoTmpT2HjjBWrs2AUNLoCmHaJ58wWqb98FDVkcyOBw883hdgr67rXXglCVEPWrKS505pln1jT1RhavS9FpEVfOqHk89ti2as6ctllAoe\/B+mGN2JwaWB9XXrlA20Jtawgm43o\/M+aC8RR9fpIAX0+0xM0NdWX4FK1ZmQJeQn379lXTp08PTW+NNxFlNu+8805de1kakgWeS5tuuqkupDNx4kT9O8V4xo8frwOHyIlE3iOS46FSsjGFuBrNGMXZgNkcFi4cF4szL2rLluepFi2mN+q3ZElbxW9hbfHinf4dxHdmIyM2fddaa6j60Y\/GF21em\/R4wuaR+XNZF8zl55+f6bx+\/vGPQ9SiRUNC+7M+uB9rhOayhpLcr0lPch0Qz2GRT9GalSlIvMH555+vmQISAgFsm222mU6lTSOojUA0ynXSJBMqNgkM1KSzoFYzicdwbx0yZIg+OT\/22GMKQzYSSKU1mnkujIHP5Mlt1cUXh3NgUQv16PFi6rmIYx6pb+ovzB0B13Xh2k8GjATyhz+0DQ3Q3myz8HUWt4bS3C938PwDMkeglJKCbO6c4LElkAI7aGimuA72ADJOogIiL1JYmgtcU0mTAVPAHoFB+plnntExC6TQRpowW9IazcEZq0cDXOarsgne0HVduPYLQugaw+UKfdb3c32u79d0EXCqvGbmMArLVxT8LowpUKOZrKaEi5NJVVxSsSk88sgj2mZRSY3muCn0L1bTXeBZrAu\/fvz6aUoIxDIFPIVOOukkHZFMwZy4tNQYnwle23DDDRvUR6aEMXPmTB1CT2oLPJVoSRPiNaWJ8bR6BDwCHoFaIGB1ScXITEZUGikuUClhJObfuJ8SmIEKiEjme+65pyF9NrYI7AkXXnihOuiggxpUSiSLOvvss7Ud4ZxzztHXY5NAUvDNI+AR8Ah4BGqLQCxTYOMmeO3555\/XHkYkr5s8ebKusYD94K233tKnfdRF5AQxXVIpwHPGGWeol156SSfMY+PHbZW4BKqzwVhwPSNnErUYwmou1Baa6Kf\/8Y9\/1MF8NALqJPW3rX40\/cPqSge\/rzbdWdHDWmB+aaRVR1X405\/+tKrkJKWFwUlWX+qCSBBl1FxWlRilVFb0FGFuwC4JPTIvJOOksReJQ0oZ5yeOnqLMDzg72RSogUx+I9xTg43U12x0JLcLbuyS++iDDz7Q3kkwAtLPPvroo\/o2++yzj06aJ8V5qv3CpXne7NmzG9WGZjKhD+8smCaSEU3qR\/NvJCLX78G6mi0reoLlUatJgzwrKS1gbaZTkay7cSnaq0lXVvSEMb1q0pF2fqjoiBYCRiDzRL607t27h75TRX93oug54IADGpUWrsXcmM+MZQokwuMUv8UWW+gEdriVEm+Aygf30t12201BUJIADILXhLlgqyiThBA2WVIN7uCDD9bBfMH60VwTVlc66vtguo9qL5C09BBvIsyv2tJBFEY2WpDySPB46KGH6kMLmw\/4xxVzqvZ8mM9LS08wGLSWNCShJyr1DZkVomq115I22\/xE0WMyuSK8O9bU2ah9SHMxbNgwp5KZtZyUWjxbXGdZqGH1o5N+f8QRR9SCjIZnVkJP7969G+5jivq1IshGi2Atp1CTKUTVAq8VLXKIIDVM1JqKo6doc5OEHvqajI34pzLPTxg9RZofa\/AatgFOUZTgPO2003R8AUFsSRtBZQSroUIilqFHjx76nmWWFNCPingbVT+6TEyhEnpMZibqii5duqhaMTkXWsrEFCqhx3xXizA3jCcJPUH1V1yt9qT7Ulb9K6GnaPNjdUm99NJL1V\/\/+tdYd1SIQhTHXjBixAitYkKdADMh6vm5557TBXhMmwSeS3zHp+Kav1nNbIL7yClUNpYolUNZ1EeV0hMmGkN7XI3uBHAn6upKi4w5TFIoknqiUnqC4HG\/Ws2NKSHY3h3mJ0z1VTT1XpL5cVHl1Xp+nJgCp3xbI9Ppu+++q+bPn9\/QFZsBBlYMeBhjSay3zTbbaO8lPJFIlwET2SOY3cv2sBr\/bgbZyVCijJP8XmRDs7ykzIG5sSel5+2331bTpk1rZBQUdUw1pyvJ3IhhMsgUimJoTjo3UfSwiRZhbpLSQ38OmjijmEbkss5PFD1Fmh\/GaFUfubzQGI+xPdx4443aA4cTAGm3hw4dqtNaUF2N383Nn+yrFOchhiEs5bbLc2vRx3SFk+dLoSCKBoluMKx+NP1dvq8mXVnSY7rV1cKmkIYWsA4yBb6LqgVe9LmJoqfWcxPE1OXdMccs\/cWlu4zzE0dPEeZHMM6EKWAnICCN8pxkR1133XX1\/XGpI4qZKOdggJpc44PXqrnN+Gd5BDwCHoF4BJyZAvYAbAOchlEF4QXRqVMntcsuu6hvv\/1W5zL6yU9+otNlmymww75nSGI8In2GBH\/5yfIIeAQ8Ah6B2iLgZFMgr9FTTz2lcAULNpLbde7cWessiVsgt5FnCrWdVP90j4BHwCOQFgFrRPPRRx+tPvzwQ+0hRNnNPffcU+cpIrCNYDYMxl9\/\/bWi2A42A88U0k6Fv84j4BHwCNQegUZMAeMwOY1uvvlmPTIMyMQRSAnNsJgC+Y1aCzAMkuB5SaH2E+tH4BHwCHgE0iCwnKSwcOFC7QK2aNEirffHWEx+IyqtRTEFvDf4YGNAmlh55ZX1WLBDXH311apVq1Y6Bbd8z2\/ffPONdkfF0OxtCmmmzl\/jEfAIeASyRyBWffTOO+\/o3EYEot1\/\/\/2RT8eXmLTZqJGSRiibWUazJ8\/f0SPgEfAIeASSIGBNiEdhaQLT8CpCRRRslOOkotqsWbNUt27dGkkDLgNZY4011OGHH65TX\/jmEfAIeAQ8ArVFwCkhHgFoZJMkZYUpCWBPoABPv379tERB9DKeSr55BDwCHgGPQDkRsMYpUIKT4iMkxhOjskmqGKKxOdx5552NwtHLCYkftUfAI+ARaLoIWJmCVEoyvZCCcGFgPuWUU9Tee+8dqmJquvB6yj0CHgGPQLkQiGUK2AtQCb355ptq+PDhOp\/RFVdcoSkkb1GbNm3Ugw8+qKZMmaLTW\/B\/3zwCHgGPgEegvAhYg9dIU0HEMjWJkRZuv\/12Rc1UcSP97LPPdPrrrl27qpNPPrm8SPiRewQ8Ah4Bj4BbjeaddtpJZ\/8kgnnChAm6tgIqowMPPFDtt99+mlHgvkrgm\/ci8qvKI1AlBLp1U2revOiHdeig1JQpSmXdjye63rNKUPjHZIdArKRAJtOzzjpL1z0gER7qo2CjZsKOO+6oM6T6ILTsJsbfySNgRaBjRztTmDtXqaz7MTDXe1qJ8B2KhoDVJfWMM85QjzzyiI5oRoW0\/vrr66hlCtVvtdVWOq0Fldk6duyo7r77bqc6zpK\/HkYjtQgkNYYAFNXHzDu+wQYbqJEjR6oiFLsu2sT68TQBBFw35qz7eaZQ14srlimQ5uL0009Xzz\/\/vE6P3bNnT0V5zuuvv16rkiiQQ1U2MqhK1LNLnIKUryO2gcA3XF7DyjmSntvsQ\/RzVP+6niVPnEcgDIGsN3vX+3mmUNfr0cnQvPnmmysMykQ0X3XVVbrkJsxizpw5Ghwintdbbz2nOIVglSuz4HVQSpByjtKH\/ElS2jIoHcCcxo8fr9q2bau9oJBo+NSyoW5DBYedpYx1qE3s6okW6KoHelp06hSrPlrStq1aMnu2yrof+Dnfs2dP+xgnTWr0mtbD3Li8O+wJRdwXrGku2PxFdSRZU\/FCQorAVZX26KOP6rrLuKWutdZasftwsHA1Gz6SSLA4j7n5Sx9UVlRyk2aWfAyWYqQe9HHHHVdLnqBw6SWxYMuWLUsfv1FPtLAo6oGedrvvrprH1E+HKcx\/+mmVdT\/wy+Oe8rLWw9yYG08UPRymbftlLTYwq02BNNr33nuvjlNARSQNbk4WVJgD+ZFIm00Amy0hXiVMwWQcUrmtS5cuuia0MIUhQ4YURlJgjB999JGWWFq0aFGL+c3smfVEC6DUAz2207rq0EEtnjXLeqpP2s9FUkhzT1ms9TA35ouHRoWaNDvssEOjfaCUkgKEof4gYO2xxx7TbqhDhw7VmxwnYIzQb7zxhkK9RHEd054AcyDZXbNmzRptTJWoj1AnmQ3bBI3vhSncddddy9knMtsZE96ItOCkIm\/dunXpmUI90cI01gU9rjaArPsBYB73XPZ+lWVuxBsYz9+wNnXqVEUGaf5Kw06KcwxxXUVtVpsCldfeeuutxOOPS4md1tDMIMaNG6eD51BltWvXTg0ePFgzAc8UEk9RogvK8qK6ElUX9CyLFcCetqB583DSiVNA5RoXz8CVSfrRP497LqMALQQu7qhXiihhv\/eeUqNGKWXs9QhGasAApdq3V9quiXv+AL6IaDAGskIUscUyhc8\/\/1y7nyIN9OjRQxMbbNRtxqYAB9xll10aJIO4lNimu6lpF4BZUNKTTT6sD9\/16tVLG7ppMAXqPKy99tqeKeS8uupiEzUwqhd6YAjnnXeemj59es4rwN\/eFQE0KS+\/\/HJoXJd5D9IDFVFiiGUKeByRumLbbbdtSHMRBIYT+29\/+1v12muvORmaXYEN64c0AOOACxPXYLqniqQA8yICuwiNEw+MFe+jIp54kmBUT7RAd73Qg83qlltuUWJLSzKnvm\/2CMCcUbGjNubgEdeQFJAYitZimQJSALmPdt55Z623j2ps1NRVoO4C+vO8Ghs\/6iMMzjSYQtDQTAoO3zwCTQUBDhu8c0WypTUV7KMOrqQEcmEKXB9WjqDW+MUyhQ8++EDnPMLvHy8kUloEGyfhY489VoPwwAMPhKqYsiLSM4WskPT3qRcE8mQK8+bNU6NHj1b8RT2Mmzd\/fYtGQDQWrkyBQ2zRMLUamkl6hw6fxHdIBKuttloDIiwWbA4Yotu3b6\/uu+++XP1uXdRHeEqF2T5qsZARH1HBEadQxCCVJJjUEy3QXS\/0oK5AfZSlpCDMIMxQyneXXXZZkqVTtb7moTGYNodYJ\/YoHGDQ5e+zzz65jCspUyidpMCAX3jhBR1\/gGvqSiutpDbeeGP9lyR5LB5cTzfccEOtOnLleGlzH+G\/TGoNMTQTYX3HHXdo\/1\/vfZTLGm+4ab0YZoWgeqEnj3WPG2Wc50xRdeE33XSTjgfgABtMmyPzPnv2bPXEE0\/ojAx5tCRMASMzDKpozVp5jQG\/\/\/77OlvqjBkz1HfffdfgYbTiiitq9RLeD0lSZqd1SYXL9+3bVzMqmnk6yuPlqHSy6mXjkZN1vcRc1BM9Wa97Dnokt7S1LD1nJBAVV3NJcsnzSalDwk0aGQ7OPvts7emIWpu+ptu7bPY4moRt+kgKeAWhzeBanGMIes26uTKFIscrODEFgIMZ3HPPPboSGw0ux4QlDdOuJHgtLvdR1i9HFovFM4UsUMznHvUyN2Hrno2dT5qGDWEUTviWhrSAjSFNY0M0tQps2DQ2afYHJBVUVDi6sImTOYFnoRY2PQ7N60Q9hISAxABzID8aTIS55oNGY5111qmKpMDzscOGNWgHv6Kq4ZyYwldffaVdp2688UbtWonhmbxHVGQbOHBgo\/QXtkVSSZoLl9xHpNzAJbUICfFYiIizJAsM6jhtOBXt93qiRSSFepibV155Ref4MqVmm\/qn1msraJcw0+EzNjMlvrnxI1HccMMNisMhsUkiHaCtwEuSVPzSkAR+\/vOfazf5LbfcUqu\/q8kUKDzGOwOtZgwJjO2iiy7STKF0aS6QDGAGVFQDYIIxYAJICm+\/\/bYOGnv11Ve1BEG6CybFZeOrhCm45D6SReET4mX76jeVJGXZopb\/3dgIUauYTKHMkoKJGFLQtGnTtHOLBM+GSQpyjaiDTGmDRJ3bbLONmjRpkjY0V0tSoAol6i0a7w4BwNSfwQ4qrVQJ8SDgyiuv1LYDGAPeRf369dO1Dc4\/\/3z10ksvabowRNOHv3vttZdWJwXdVoP6QgI7uHcwLbYZB+GiYuL5YbmPfEK8fDaiektSVi\/0sHGeeOKJmXkfkaenG+kzLC0vmwKPZTPFq5EknASq0sJsChToor4Ledc4BJrp9EXCwPMvyBTItMy1cbFXNvqjfhd13pgxY7QDDC1qrZVGUuCUQRQzf+HOX3\/9tbbUI4Ihug0bNkyLPkgG7777rtb\/YYDEOAUQiHW2ltbQzH05OTCZQcbhbQo21Cv7vV508IJCvdCTx7qHKZhJ3IIrp1beR0H1kcuKpgAYTMFlX3K5n61P2HyUba0tZ1NgY+ckP2jQIF1yEykAjsbGDyNAJ2ZKBM8++6zm6ocddpg699xzG2VKjQIwKvcRDOaZZ57RzGbvvffWz0RqCeZHwi+bBrdH8qDl8XLYFoDtd\/H3roegn3qihXmrF3ryWvdRdolaximkYQq2dzTr36MM\/xjwy7IPNGIK6L6wGRBlh4Fm3XXXbcAMYjFo4ZpK+cwf\/\/jH2p4gJTsxQKNHS+KaGpwQFwnCZChFd0nN64XNeiG73K+eaCnqIcJlHoJ98pyXYEQzySqLmMAtDW55XRM2H3nOUR50NGIKsuH+5Cc\/aVQJjY3\/tNNO0\/mNkBzQ+eH2BRPACERtBTw50P+lFdNc7AhyUjj00EO1VCJ2CfMlF++jPMBKek\/iO4jhKNKYktIg\/euJFmiqF3qEjiwjmtOuEX\/dDxoL852PWmtoXYqSfcGcOytTwJXr4osvVhMmTFAEqyEhwDRgCJdccokupQkjwAXU1aYQtnhcvJLkuiAD4XufQti\/kk0ZAc8UijH7IhW4jAbGwadorRFT+PLLL7VR2VQFiY0BVdJuu+2mq6nhnopUgJTAX2wMpL8YP358avVRpUxBGAPMwTePQFNBQFI1e6aw9JQuWZTj3OPNQLes14kwBZdU5qWQFMgxTzDa448\/rqMakQjwCybPEaklMECjZ4RwURXhmoa3EgQiTYRlUg0DPo2rapykkPXk+vt5BMqAQB76aiqKTZsWTf0ee5DRoHjouOQ+YtTVYAplZtLLeR+JNxGBFngVwSRgChiQYQQAKkyBjR3\/YQrsUAWNNBhpbQpMlouhmX5h6qPiLVE\/Io9A\/gjkwRT+nQ9Pl5aMavyWZaLUPHMfmY4pkitJAtk22WSTRpHQBOkSi4XrfVh+JRMPc8x4QRKFTRAvwYREWNcVU0BaIJ0FH4zKZEHlA2C4nVIT+bnnntNxC4888ojOS0JlNvqMGDEicS4kE+ioEpySC0UYjmcK+W82\/gnlQKAemEKeuY\/koEm0szyHmSW6mcA2crdJygzcRvv3769Iv4+HJXmU5JrttttOx2xJlLI5ZlNtVZdMAcBwNcVuALEEsbAJo04iQA2QYQoAiYQAeBMnTtT2BgzPZa8bUI6twI\/SI7AUgXC\/eOIw0iM0evTSwvRRjXrzKXPh6VtSp8es1ZNX7iMYgZkcT+gx1UfmsyVCmjIAwfxKBPISQ\/XJJ5\/o+KjJkyc3JN1jf+RAXLeSQnAhsPETwUgyKVRFgEE+Ebjpe++9p43MwjSQGJpCY1FdcMEFmlQzfa\/peWCm5k36fbUxzIoe8wUzk5pVk56ktDA2UQXI6dDcbPl3XmmWXXBxoYf7NE6IF6\/+cXlunn2C6ifz1G0+N4vcR2zWxFdw6pf7kaUUSYH9i7\/8RnI9m6RAfw6\/9957rz44r7zyyjqza5ikgMcmNlhz\/ZgqJ76PCsqt1bsj2FuzpKJCevjhh7VbKqqiYKPgDi9Nr169tAqp3huLR\/K8Iy2xEaJ\/JCcUGElq8WAueNfvK7HJpME+K3rw9jCTlaUZS6XXJKUFrMOCIU1POMZkznelY0xyvQs9f\/7zn3UsTOOEeOWSFIKbZZa5j7p3795gNwjaFJgL9P80GAfZn7GjUrhLbApIBSTixKsR93xSZoA1GVjZ96jNEGZT4Fmo2s31gy0DRx0zTc+RRx6p7Ri1fnfMdWllCtKZADWYA6oicnngtgrXhKAiBmAkefkq6SunBFKCXH\/99doYb26Q3BvG4fp9VMWoSsaY5Nq09GC0Q5KE+ZmJyZI8O+u+Nlp4cTEqBoMhuS5szoo4NzNnztSSfJaGzWobml3nvRppLqKeAYPG3XXTTTfVsVq479NEUpCcbKIV2GijjbQ9Iu6wJPYOGFeR3h1npuA6cU2tn0wsYqj4SIMBnL9Lly5aPE3yfR7VoJLMSSX0yKmL55micZLnZ9nXRouZapm8XhIhb6oDzLks4txg+\/NMIbtV48p4TAnTVPcIU+DdR3UVtX5MaZRiQkV6dzxTqGA9cRIQcTBqIykTU6iEHnPDFHUAL0atNlIXWsrEFKLoyYMplDVOoYJXObNLXZhCmA1LBlCEd8czhZTLwXR14xZRKoeyqLLeLJwAABg4SURBVI8qpSeoWjFrXaSEOPVlrrTImMPSsBdJfRRHTx7qo9TA+wsbvMGi1EfBzA1hkNXy3WE8nimkWMhMmng0yOVRxkl+F30h\/7YZoKttaBamVSk9VOOLqnWRAuLUlySZm6i4l6IYml3mJszQnBo8f2HFCIikEGZo5ubBmCs5UBbh3RHiPVNIuAzCEl6Jh4IErnDLsLTert8nHFJF3bOkx3RJrYVNIQ0tgBcWDGneK0sjbpLJcqWHe9ZDJt4k2BS1r2RENV1SZf0E4zGgQdyda\/3umHh6plDU1eXH5RFwQMBnB3YAqcpddtppJ0VCvLJ6ZXqmUOUF4x\/nEcgaARiDzw6cNarp71fU7KeuFHmm4IqU7+cR8Ah4BJoAAhUzBUp4Yuyibb311mqVVVapS9gIa3\/yySe1MbVnz55ql112aRIR3HU5mZ4oj4BHIBKBipmCGOl4QiXlOIs+R++++67O+XTQQQepm2++WYe2ky7cN4+AR8AjUE8IVMwUyAdCHQXa4YcfnrryWrVBDavSZCYfC\/M4Ia34jBkz1EsvvaROOOGE0tBabWz98zwCHoHyIlAxUygj6bL5iysp+UnM5GP43JOagkyw\/BsV2WabbabThROSzvWokAhQ8c0j4BHwCNQTAjVjCoRzk5+cBHskUcsyaMsMPsITQPIQSfIqUuDSnn\/+eZ0BEabARi\/\/l4pyXNeqVStdbIiEgN9++60uLkR1Ogpu7L333vW0FjwtHgGPgEfALaIZIyub4GqrraaNq6hR2NDRsbNpcmp23dTZXNmASRZFKm6zHkGW88HJHxUPKXDDVEEmExCmEExrKwnSGBdqMmwJVGqCqZx11lnONGdJl7+XR8Aj4BHIEwGrpMCp+9JLL1Urrrhiw6n6vvvu06dvmANtxx131CmIYRBhjX6vvPKKuv3223UlN7muc+fOusJR165dc6nYJrUORBowx5aUKci133zzjS6u4ZtHwCPgEahHBKxMYcyYMbqCEMUgKCSD1HD66acrwrkJ0ebkj\/qHzZ3vzcZv999\/v65otHDhQv3T6quvrtPEUuUK7528CvOw6d999916XG+88YZOi2xjCkH1UZFqA9Tj4vM0eQQ8AsVDIJYpfPXVV1oiQJd+9dVXa28bynMed9xxiqIy1GSmDRw4UG\/61157rVp11VXVW2+9pd1TJ0yYoKUC6jZTqhP1DNcFN+isYUG6oSLS8csKyYaV+wtKCmGG5jAJI+ux+vt5BDwCJUVACmGbBadLSoo57FimIDEIpBiWjZzTNyfoG2+8Ue277776XqhpnnvuOb3hs9nCFGhUKWJjpnQd6icKmZj3qiV+QabAWMQrqdY1UmuJi3+2R8AjYEGAghOUp+OvNBjDyJFKde1aevgSMQX06dgXJk+erEaOHKk9cfDMueyyyxR2BozR66yzjjrssMP0hwIzzZo10yCFMZjSo+cJ8Ah4BLJBoFu3+MLSbLpTpix9VpK+2Yzuh7vYapXCGPr0yfqpVb1fLFPAPx+pAC+ha665Rvvoc9rHPx9VEeokEnHxHb9Ru5R6t926ddN1elEbSfNMoarz6h\/mEcgXAdeN2bVfx452pjB37lKakvTNEgUkA+ixNZhXiSUGq6H50Ucf1e6XGIWRFDAwDxgwQEsCGJDHjh2rvvjiC3XyySdrFRLlAWmtW7dWxxxzjDrggAMU6pjPPvssE\/WRlKvDZkCrRd5+25rwv3sE6h4B14056361ZAowBFNlFDXJSApIDCVtVqaAoZg0FtgNFi9erDf\/0047TZOLERqPHVRK++23n\/5u\/vz5WpX0P\/\/zPw0eR1tssYXaZ5991L333qsDvioxNJu1akX6wDOqVrWASzrvftgegeURcD3VJ9mYa8kUktBjWw8YlaHFtX3\/vWvPwvWzMoWoEWNLIDCsZcuWoX77UbEJxDQgeWy33XaZ+PtL\/Vqylo4fP14dcsghhSluAQao3nDDNVVphVsFDgOqJ1ogt0nR47o5um7gZWEKSeixvQNJmQKqrmVeSWVba6mZgg1D83fiFTBOo24SzySJVyCJnmmQTnJfM50FNg3iH6Qs4frrr6\/41LIhWX344YdqvfXW06k0ytzqiRbmoSnR06JTp1h9\/ZK2bdWS2bOVaz\/wc+2bdb88nu36XrZI8A4v\/uc\/G24btdY4KBbxsOjMFLAlPP7447qmAESuueaa2qBMigukBZdG4BtMAeZAUBvxD2nTXIhtgSA43FyD9WxxhSWeopYNQz0MEXzKXmeinmhhTTQletrtvrtqvmBB5KsAU5j\/9NPKtR83cu2bdb88nu26R7Q+6ijVYvp0a\/fFO+2kFo4b19Avaq2xh5I2p2jNyhQQfdjEBw8erDfxYOPET8QzNoUk0cmVJMQjiR4GZhoBZpLoDkmB2qht2rTRUkKtJQVo\/Oijj\/Q4WrRoUbS5TzQeGy3du3fXwYlRrUOHDmrSpEmqe\/fmKqablrgnTVqaPiXPZqMnz2fnce84emyndUBfPGuW9fQv\/VxO60nv6TrGPJ7tOh\/Nn31WNe\/ePb57hw5qyYgRasmuuzb0mzNnjtYY7LDDDo32gdJKCk899ZQ65ZRTdERyv3791DbbbKMT46Erp+IarqoUoLnlllv0iT3v9s4772g10bBhw9Qaa6yhrrrqKj0GUlzzfVjyu7zHFHV\/JCoivfHEKjtTsNHSsWNHK1OYO3duzbwJg3Nko6dWaybtc2PpMXTrC5o3V3yWa7hRnnBCvFsoF0msgGvfrPsxhjzu6Qr86NFKjRoV3btrV9V2xAht1yTP2+WXX67\/SuNwRIwX+d6K2mIlBVHRYFAmgjns5I1aCY+krbbaStcfQFQidiEYp\/C3v\/1NZxklMR6bZN++fXVJSwlucwUIm4G4o8o15DdClQVTIHKaCalVC56CJc2HOR5Ow3Gn5R8WUK2oWPpcV1roay78qFHzIrh49FXrfWF9r7TSSoXU67rMvPP8LAOdk+lXLVuqD0sutbpgU8s+7JMcnFG1x70LvA8E\/hatJYpoDhs8XkgkxiM+AaMxQW1ww1tvvbVBX8Yp\/tRTT1WcFM2G5MEnibHFrJjGvaRWAsZqmELwGUUD3I\/HI1ArBJBWOZCJirVW46jn506fPl0NHTpUawiQ3uIaTGGKSF4FAiWWKXz55Zc68ylVxy644IJImwFuoU8\/\/bT6+OOPNWkwgKOPPlobV5EckCCIW+AvEc\/o2RGrXnjhBTVq1Ci1\/fbbO0Mybdo0ddFFF2l3WBrR1Ugdnik4Q+g7NlEEhCkUScVab1MhDi91yxSQAlD5PPzww9pmEKaWwYDCaZ8FR9Gdm266Se2xxx4Nc00wG8VuiEsgm6p44cycOVNHOBPxTDCbq7SADeHBBx9Ujz32mObEJOUjunrPPfcshE3B1TU6634AnvU9Xe+39NnZ2hSSPds9Q4IszHqwKThjtCxO4cUPP1S9W7fOxe6GkwEOKfxFU4D3H3+bWkvCFMCGPbZorZGkgMsoKSvMgZLkDnETm8A555yjDc3kPKISGe6leCXBLOhHqm2YAjmQpGGo5iSPismMOub6s88+W2\/swWviQDJjE+gHY+nTp4\/afPPNPVNwTB\/jupm49vNMoTavdZL5YYSyYWUpKQgzIPVNsPFd1jpztBIcUGlmNmOz0iK\/BdPfmJkQZJymG7u5PwXd3ekfvD\/fhbnTJ2UKqLuLxjwbMQVJG\/H6668nWuUwEeoo7L777tqmYAZq4RlEum3Jqio3FuBhNtRecC3naS5u\/r311ltrplWU4LVOnVrEGpHbtl2iZs9eorLuBxZZ39P1fkuf3SnW+4iDAy+W6z1d+yWh21zURQ1e69nTvn7EZTcJRtCOkwexO1kyBdTAYQxBsObAxrufRWPDRX0saXJYT0gn\/fv3V9ddd52WTnBwocE80FjgEcleQ2VIPCalD3sdY4dpsV8NGjSo4Te5L+72Z5xxhv6d75544onlCokF6UrKFAovKXB6J88R0kKSBugkzoNzUoxHmIIU6cFDCcPzuuuu23DbSpgCN+F6Atew8t92220KRoahWVqtgtd2372dWrAgxOVv2cBgCk8\/PV9l3Y\/bZ31P1\/vxbOaCeY5r2J2OOqq1ev\/9aHy4Pgk+SegOPjvoGdamzRI1btzSCoG1akkwT9IXenhHkM6zYgpICagNbQ1jahYumLJZY1MMZggwmUDYRk1ddeyO\/EVjEdzkkST4HSaC5qJHjx6aCfAXRpMHUyilodk22fJ7lCpI7Am4nlKlzbQbYIiRNNwEoLG5hzVTbNt\/\/\/11sBoMgdgJJphJ5DthCrUOXrOd3lCzzpq12HpiTtoP7Gr5bNe14trPlZY86K5lgF0edAvmnGJPPPHERkyBjT0u6DBuvjil4yhia0gLUgXR1jf4O6oVU71i7gfBYlimaslUH8kmj7SKxIALO5oF8+Qv1RkJwhwxYoSWEHCtlz5h6qOwDM2ukkKR4xWsEc22SUT8QQVE3QU2ZhaKeBMRT0BcgVmlTe4nKbkBP1jbOe6ZiH1kZcWwTfZVjNsmU8jqFGSjO+p3M\/cYJ1Fa0IiOZ64tRxnX0S+J3ti1r2s\/V1rSYhV3nesYuYdr36z71ZrupPMTZlOwqX\/yoDHJPePsEqYKKKh+lk0eBxccXcRbkWezR1AMLMgURJLA01Ka2A2CTCSKBsEYCQMnnbAGQ4BRZm1vSYJrXN\/UTIENj\/oJv\/\/973WMAiIWk4R3Ead4ai8MHz5cG4MwDEngG0zk1VdfVeeff742VmNPoIKbazNPA3INC5vn12NCPNupFQywUdBsp0yxZyS5p2BcbR18kjG60p11P9c1Sz9XW4HrGIPPdpmfMJtCmSSFYK11sRWQGZmDJ5usMAc2ZzQJ0Ce2haXvylI7BG7zOMmYNoVevXqpMWPGNLqPyVyS2BRuv\/127UQDUyN2QRrSCuovmEJp01wEFx8J3khmB7CogGgYZFAFtWrVSgMqOZLgxhiAdl2WB4RCO0gOz5JDpHlzbdwhbiFJziRTr4hkMm7cuEaSgoy3VjYFE69qJ11z1dcn2cykb7VpSTJGV7pddfCu\/ZKM0fWerv2Cz3aZn6xtCkSxk0nA1rKyKfCc4KFQNANB9Q6qZg6esvGbkoTYH7if2CHxPkJSMA3Z\/C7SCAZ69i5T4gir5S6SAg43SBk05uaNN97QWR\/MxJilTYgHUVG1ETp37qyZwW677dZg+IFpMEGkr0C9A8OQJoZngtyIbcDWYGMIQZsCSaWCYhcFfKjy5hPi2V7P9L\/XQwI52yk8jR3HFdG8n+0yP2E2BdfxR\/WDKcSlOMnS+6jSsVbjemEKSBzsVbSouSmlpIAL14QJE7SKyJQKEIs4icdFOec1AcEaChiwL774Ys2BfUK8vFBfWn+g7Mn98rAp2GxDUm8+j2ebs+0yP3nEKTCGKLtEHnEK+a3wbO4chrHL3GTz9GzuspxNIazmAUnDEMdQ9ZBWggR4eP1UUlYz7fBNX2UzD5J4HxWpyA7eC3feead2gStagEpS\/OuBFttpPWkMCRhmfc8kthRzDl3mJ484BRlDMKIZPX4WbqhJ12mt+wtTwKYgWaOj5qYUkoKp8wdc1EPo0ghKEzWQBLhFMQWYCtHNuJiaqiFUR+TTf+2113RSLuovpMlmKowAjycYgbipFiVOwVyUosM19Yu1XrRpn18PtJi2B9Kz0IKZf\/OKkUhrK3CdL5f5ydqm4Dq2ptQvzKYQNTelsCnIho9hhDrKXbp00SktzBbHFDCmUHAHVY6ZugIbAq6nL730UsOtkD4w7mDxt9kVzOfjDUCENF5LBLAEs6RKjWZedspg1rJxQsB99te\/\/nWipH+1HHPUs+uJFmiMo+ekkza2pjafNGmOhuqEE7rG9m3efIGir2u\/tHPvMj\/SJ0u3bbJyT5sWPWrSoGWZCr0aaS6CVRzD4hGiKJZrzXc+am6CMRhp5z7r6xqpjzCIEB3MB3sC4g0iIAEv+PuuvPLK2hqPcTkoKSA+olbiL0nu0CfCUDBSk0qWojgYnbiWojwkxyNymuAXUiREtaCh+eCDD1b33Xef9jiCKbBIaGRlPe+88xq5f2UNlr+fRyCIwPz5T6slS6Lrd8AU2rXbXbn2qwbCWTKFf2e5+Pe7Hj1qfsuqZEA10lwQjyAejRI1zR7DBm7mbrMxBZd5RNXNp2gtNE5BSmWSs2TGjBl6Y0d9xCkcbwPcvHAzNW0KWNuREoJuphK5jKqIPEiihsItFV9hGAnBa67Sgmlo5p4iKTBhcGQ+vnkEqoXAUUf9zFpedNy4F9UFF8T3Y7z0y7NJrv+yMoW801ygXUDFbeZQSjofIim41Kxg\/0qjQk86pqT9rcFrYXEJPOTAAw\/UUcwUpUfdxL\/J+BfMcQRI2CXYvJE4pEmtBiKT0bkH1VRxhJjiXRLRLik4vr9HwIaAq1eR7T7V+D3MM4bqbS5VAKPGZ6tO2aePUscfn546vLfMDNx5prnA7on+X4Lg0uwzeXl4pUcw+ZVWpiC3lFiFe++9V5fDlAA1jNFUXENSQMQStY5ch22BiGZ+p86ztEoT4iUn1V\/hEcgegbIzBZv6J3vEkt0xTv2UdZqLKEnBNRkelDUppmBOJfaGxx9\/XOF2RU2FqNTZIkHMmjVLG4bxOvJMIdlL4XsXGwFbnAKjD1ShrRlBZZcU8k5zgdoozKYQVoshahKbLFMQQHA\/Jb\/I2LFjdeoLwrqRDEQVJPYEah5QihNDtbRPPvlE2xNIp51UfVSzt2rZg1kkkjTLLLQRVbQj6ffVpi8reqI8Q6pJT1JaGFtYUZWoOasmLTwrK3rMuSmrTQE88kxzIYbkLLyPCKjFxkqTAj6yztC00EzVdxHeHVnbzuqjuJcB1RJRjRMnTmzkTUQeEdLU4mkUtNyLofnYY4\/VXkOuhuZqv5TB5yFKUhIUozn5VJhM8qGQZ4WFgG2FJn34NxXrXL9PUmwoCyyyokfcg6mrIEE7WYwvyT2S0gLWZmEp2SyD1f3M+U4ynkr7ZkWPbEYcYNioyswUKsU07+uFoYA1DjnmXoARG+9MHHRk3R155JHaWxObay3fHROXTJgCN5RNnqRPEA2jYAFiUEZ1ZEb0UoyFjfIvf\/mLzqQqCfPynrA87i\/BdLjKXn\/99Y3iJ5hkOd2YcRVx39dqQxVs0tKzySabNDA\/qX6VB95J7mmjhReX\/PpE6p977rl63YI\/18Hsg3NW9LmJokeYHAGjHMCyZArVjlNIMv+16CtMgcwPSHlxhyVxdaWGgxwci\/DuZMYUYAKkr+UjLZgJlVTZqIrGjx+v4yBIiscnWG+gFpOZ9pkysWRYFP9m7mUG1SX53sUXOu1YXa6rhB6z8l0RvMJstAjWcmozmULYnBV9buLoMecmS6bgsqaaUh9hCgT+UqzH3AvM9RPM4VakdyczpgDx2BiorUDkMuogOCAnSFENSYU2\/KWJcCazqZlKtmyLxzRAmXmYysoUKqHHXPCiruDFqNVG6kJLmZhCJfSYUqBZc6Rs71sZxosWBGksjimE2bCEtiK8O5kyBduk4aWEhEBsgml0tl1XxN+DUY5RKoeyqI8qpSeoWpFI81okTXSlRcYcJikUSX1UKT3y\/hDY6aP+899NJEMDDjhB9ZEpIUSpimr57oBOVZlC\/tNRnSeEFQmPMk4yoiIbmoVpmdWp+C4pPW+\/\/XZDgZLgJludWVn6lCRzI0b94HiLYmjOih4zPcTMmTNV\/\/79tVFz8803r+bU6GehlsP2YT6bgyLp+dEc0O644w512mmn6X9TzAtXUbMuS1R\/s0+1CAujB+9LHFHCnEtwyDErxDFOc35q+e4IZp4pJFw9QXc1Lg\/L1Grqbc1rXL5POKSKumdJT1Th9IoGmODiNLQIAyQnl9gU5EUVPW+tdPBZ0lPruQliKtMa9+6Eld4V986odyrBcqm4a9L5iaOnCPPjmULFS8LfwCPgEfAI1B8CXlKovzn1FHkEPAIegdQIeKaQGjp\/oUfAI+ARqD8EPFOovzn1FHkEPAIegdQIeKaQGjp\/oUfAI+ARqD8EPFOovzn1FHkEPAIegdQIeKaQGjp\/oUfAI+ARqD8EPFOovzn1FHkEPAIegdQIeKaQGjp\/oUfAI+ARqD8EPFOovzn1FHkEPAIegdQIeKaQGjp\/oUfAI+ARqD8EPFOovzn1FHkEPAIegdQIeKaQGjp\/oUfAI+ARqD8EPFOovzn1FHkEPAIegdQIeKaQGjp\/oUfAI+ARqD8EPFOovzn1FHkEPAIegdQIeKaQGjp\/oUfAI+ARqD8E\/h8oL02oDauonwAAAABJRU5ErkJggg==","height":234,"width":389}}
%---
%[output:3d619543]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('tabular\/parenAssign', 'C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\datatypes\\tabular\\@tabular\\parenAssign.m', 115)\" style=\"font-weight:bold\"> () <\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\datatypes\\tabular\\@tabular\\parenAssign.m',115,0)\">line 115<\/a>)\nRight hand side of an assignment into a table must be another table or a cell array."}}
%---
%[output:59e5c07c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98f1737e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a259ee5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:30107ed0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:86913e46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:86ae1794]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0d65490a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1fb13ebd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:17813802]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62813a08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57d58dda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7dbfbb5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:281b8db5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9703e02a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75e1381e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:433c385b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9012bd15]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:432e6906]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:98bdac5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:91e26909]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3d068893]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:89132f54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:20dd97a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59bbc593]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c998b4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61384b29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:676b78ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1de67dfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1555129b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:63316fac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6852646d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49ad3ea2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:31881498]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:89c25b2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:431a590d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3d4f5acc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44b6fe82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c1b98ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0dceedde]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9dafe738]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d5d3991]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b49fca7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:57e28deb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:229a40b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f66856c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9b56944a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:98c72be3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b70f32f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22372109]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97dfb77a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93ed799a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a30654e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ad853e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:94afcfa6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3cb0572b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:10869dee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:511a9ba0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:11c2cc34]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:28efcf19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:8c9fc4d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7d6adce1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:92d1bc50]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:771ddec4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:48187208]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:5b471001]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:55e80d49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b9f215a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9adc4cad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:355379ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:3a82d0c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:8373293f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3d0d2235]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:67ac51ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9cfe8fa2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:22d5e231]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:302840ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:550ea67e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:372051f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:083af83f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:860b4267]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6b7bcdb1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d88d7fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:515cd56f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37f4ed21]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4b21e9b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1b7fe61b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:310232db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6487bf70]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f4d62cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:66532130]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4c5523b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06d28da8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:00135a40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a07706f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:858bea75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:628d544e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:00aecbf3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:9610ba5f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82fac96f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:12419720]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0e3b2499]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87bf7e4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:3a16d8ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:77e60e6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65915a7b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6bacbedd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93276de6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6db1c725]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ad67874]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69b7c3e9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:190530e9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:46ed0118]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:370377e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:309bfdfb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f1bb495]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:55542bd7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15a6f9f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44b2e047]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1359598b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8ce27118]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9dbdb39d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:31f421e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35fc1edd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5558999d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d1a220d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:40ab9bdc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f465fd1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:394eb751]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18a11ff6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ae7b219]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6dee4860]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f88328c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5a2be3a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:094299d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:68ff7093]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:94f290f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:892e4778]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8ed20ad6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6ac1ba10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f3cd523]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d1772be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:583b9292]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ca426da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:059a2990]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2fd1e05c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c1cef6b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45f57394]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ef74ba3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3358b62c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ae26367]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51c88e2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04b6eb43]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4d249b88]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:85986063]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f90aed5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:346d758c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80b93dc6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:26b2ce6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:67d1ee70]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5055efa3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8633091b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:28754ac1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:66b79876]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71e84bf7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:463028d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8a02e87f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71f84cf7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12724ccc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:214e3706]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c89837d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3d516f1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f772d61]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:130f004f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6de7d553]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:400a7d65]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f965591]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:148c6383]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76af6f19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77849cf2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:757a5e2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06978dce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3a4db530]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2362cc63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:41595d72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51b040a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5e3c6c8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:327ea8ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8e21ce04]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3609f707]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0a1dedf0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97c71e2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:31c0cb8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69a4bc62]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c6a82ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:33cb0187]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:41263c8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e7dbdc0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47edae40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f053cf7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95a971e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:334a2840]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49204063]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6804cc6d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75d5ffed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95a232f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15857ebd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6be91978]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c694714]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:855e2d48]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:824a3d00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:83d45c79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96a1fc1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9db438f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:98fa302a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:48031cf5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2261c2dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95723cc6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:295e44d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:42f76616]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90e4dc8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:67858de1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32b5fb61]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3eaa9cc1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:635de8f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0857736b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9293f623]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12f37b2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:829ac9f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71cf94ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4703aac0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:20edc517]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:36a2bd34]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:370831e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19f3464a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c29ba4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1929e0bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c906ac7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43a60f55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:39c35e9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47cea7b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ca87854]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2668f26a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7edeab23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0140b7a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13fffbdc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:01c2ea91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:81f8a675]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3fcd517c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:41702c92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3292f211]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1168d03b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:72016153]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:530ae009]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03355161]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45adebb6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2104104e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5fa3211b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b15c5ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d16e8a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:02e090a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:456c9889]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96d5c9d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0d93a46b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6aa4ba47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2dfe4951]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e506fda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:033dfb7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:929df44a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0e38ca1a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85d6be29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12249d04]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6b9ff3a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c9e9183]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3885dcb7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:73abecde]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f45b41b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69e60037]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e8543d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6f6eb02f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:25ce1c9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b5febd5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5033a840]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3187f9a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7708404f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:710924d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:79e4fc6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:98dbb739]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51447b7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9655d290]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:135f7590]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:84d006e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c7575f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04f69dfb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f2ebffe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34eb889f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:252917bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a7635e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0721f5a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:345642f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32a01465]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0fc3f6c8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99312fe7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6c78418b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:52d68318]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0acd2952]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ee41ed6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:62789c26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15501e2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6949bef0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f75a1eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:103734d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:147192cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8dffcf54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:56e3d87e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4b660e78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7df7f2ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1da6e4a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3bc478d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8d00033c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:884fc1cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3627716b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c5ffa0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96e7b1b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97fff6f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4629151d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:060eef26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87210ed5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4b014e8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:11524587]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:2cd24c98]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0d156e82]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:328b7726]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1da066bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87f7d9b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:383cfe7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6c2778f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:186116fe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:24aa8ffb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a0ecb32]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:5fda4f5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:22fbca55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:60595b92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:5921b86b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87e484a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:236ed6c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:8dafb288]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b934d7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:51e6f097]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d60b099]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:229611ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54f1c4b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d07ab76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1fd01ed7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d60b6cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:419ff267]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b07e8b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:24acbd94]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f07897e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:104a27a1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:08595ca8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0a7d86dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4e7f48be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:98e67a13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:162edff5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7111247d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f52226f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0cd432ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3dd51346]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:792e1ccf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0291b4f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:665930fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8e9e22bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0aa4eb7c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10d457b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9be7d1ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bd576fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4bf0237d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:266a4e2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59768bd8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:21d1bb56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3cc784db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9e7aa3d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:071e8ea0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8eb9d1ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57b48c07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:33fe7673]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1a5f3c19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ab16a12]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17559201]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97a9ea44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30ffd511]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:37ab7f95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:38122745]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:131aff81]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9b8c9542]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2079614f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8cbf4444]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:682df635]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b62cd3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74ea0433]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:450fa3e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:054a135e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c7f84d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05aa815e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f8c873b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62cf130a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8dbac247]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c5c8e72]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:889e801b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:73e99599]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:25befab1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:70dab343]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2076871c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:94f630f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3cfdb78e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:59bcc8e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:236fa145]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:395f53f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:31ba1dda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:22917fba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1528eb93]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:55a00260]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"2562×12 table","name":"MLO_result_MK","rows":2562,"type":"table","value":[["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'y'","'MK'","0","-6.2725e-04","0.0028","-0.0042"],["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'MetSea'","'MK'","[0;-1;0;0;-1]","[0.0205;-0.0071;0.0009;-0.0010;NaN]","[0.0464;0.0032;0.0077;0.0060;NaN]","[-0.0059;-0.0182;-0.0063;-0.0084;NaN]"],["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'y'","'MK'","95","-0.0069","-0.0034","-0.0104"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'MetSea'","'MK'","[0;-1;-1;-1;95]","[-0.0032;-0.0075;-0.0008;-0.0053;NaN]","[0.0252;0.0028;0.0063;0.0014;NaN]","[-0.0278;-0.0186;-0.0083;-0.0122;NaN]"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0100","-0.0040","-0.0160"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[0;0;-1;-1;-1]","[-0.0253;-0.0047;-0.0078;-0.0041;NaN]","[0.0106;0.0095;0.0001;0.0039;NaN]","[-0.0629;-0.0192;-0.0159;-0.0122;NaN]"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0110","-0.0083","-0.0138"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[95;-1;95;95;95]","[-0.0429;-0.0034;-0.0042;-0.0055;NaN]","[-0.0249;0.0025;-0.0006;-0.0014;NaN]","[-0.0622;-0.0093;-0.0078;-0.0097;NaN]"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","30","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0058","-0.0048","-0.0068"],["'MLO'","2025","30","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[95;0;95;-1;95]","[-0.0218;-0.0030;-0.0041;-0.0017;NaN]","[-0.0132;-0.0004;-0.0024;-0.0003;NaN]","[-0.0307;-0.0056;-0.0059;-0.0032;NaN]"]]}}
%---
%[output:3afcb33e]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"20×19 table","name":"MLO_result_LMSlog","rows":20,"type":"table","value":[["'MLO'","2025","10","'month'","'Ba3_A82_ae33'","'abs'","'log'","'LMS'","3.1262","95","-0.0726","-0.0261","-0.1190","-5.1766","-1.8649","-8.4883","-0.5160","-0.5154","-0.5166"],["'MLO'","2025","10","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","4.8948","95","-0.0465","-0.0275","-0.0655","-34.7014","-20.5226","-48.8802","-0.3717","-0.3714","-0.3721"],["'MLO'","2025","20","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","5.1774","95","-0.0239","-0.0147","-0.0331","-55.0205","-33.7666","-76.2744","-0.3799","-0.3796","-0.3803"],["'MLO'","2025","30","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","1.7301","90","-0.0053","8.2207e-04","-0.0114","-14.7924","2.3074","-31.8922","-0.1462","-0.1458","-0.1467"],["'MLO'","2025","10","'month'","'BbsG_S2S20'","'neph'","'log'","'LMS'","0.1620","0","0.0026","0.0343","-0.0291","0.1419","1.8936","-1.6098","0.0260","0.0269","0.0251"],["'MLO'","2025","20","'month'","'BbsG_S2S20'","'neph'","'log'","'LMS'","1.0755","0","-0.0056","0.0048","-0.0161","-0.3119","0.2681","-0.8919","-0.1063","-0.1058","-0.1068"],["'MLO'","2022","10","'month'","'BaG_ae_psap_clap'","'abs'","'log'","'LMS'","3.8480","95","-0.0820","-0.0394","-0.1247","-3.1666","-1.5208","-4.8125","-0.5597","-0.5592","-0.5602"],["'MLO'","2022","20","'month'","'BaG_ae_psap_clap'","'abs'","'log'","'LMS'","0.1892","0","-0.0020","0.0192","-0.0232","-0.0794","0.7600","-0.9189","-0.0393","-0.0382","-0.0404"],["'MLO'","2022","30","'month'","'BaG_ae_psap_clap'","'abs'","'log'","'LMS'","4.0497","95","0.0256","0.0383","0.0130","0.9632","1.4389","0.4875","1.1564","1.1587","1.1542"],["'MLO'","2025","10","'month'","'expS_bg'","'neph'","'log'","'LMS'","0.2871","0","-0.0021","0.0123","-0.0164","-0.8308","4.9568","-6.6183","-0.0204","-0.0200","-0.0207"],["'MLO'","2025","20","'month'","'expS_bg'","'neph'","'log'","'LMS'","0.8500","0","-0.0022","0.0029","-0.0072","-0.8707","1.1780","-2.9194","-0.0422","-0.0419","-0.0425"],["'MLO'","2025","30","'month'","'expS_bg'","'neph'","'log'","'LMS'","1.5589","0","-0.0029","8.2071e-04","-0.0066","-1.0720","0.3033","-2.4473","-0.0833","-0.0831","-0.0836"],["'MLO'","2022","10","'month'","'BbsFg'","'neph'","'log'","'LMS'","0.9001","0","0.0048","0.0154","-0.0058","0.2499","0.8052","-0.3054","0.0488","0.0491","0.0485"],["'MLO'","2022","20","'month'","'BbsFg'","'neph'","'log'","'LMS'","3.0163","95","0.0060","0.0100","0.0020","0.3049","0.5071","0.1027","0.1274","0.1277","0.1272"]]}}
%---
%[output:435725cb]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"20×19 table","name":"MLO_resultLMSlin","rows":20,"type":"table","value":[["'MLO'","2025","10","'month'","'Ba3_A82_ae33'","'abs'","'lin'","'LMS'","3.4046","95","-0.0351","-0.0145","-0.0557","-14.2518","-5.8797","-22.6240","-1.3511","-1.3506","-1.3517"],["'MLO'","2025","10","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","3.3361","95","-0.0672","-0.0269","-0.1076","-7.6846","-3.0776","-12.2916","-1.6724","-1.6713","-1.6735"],["'MLO'","2025","20","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","3.2030","95","-0.0386","-0.0145","-0.0627","-4.0329","-1.5147","-6.5511","-1.7723","-1.7710","-1.7736"],["'MLO'","2025","30","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","2.4870","95","-0.0170","-0.0033","-0.0306","-1.7593","-0.3445","-3.1740","-1.5093","-1.5082","-1.5104"],["'MLO'","2025","10","'month'","'BbsG_S2S20'","'neph'","'lin'","'LMS'","0.6096","0","-0.0018","0.0041","-0.0077","-1.1014","2.5122","-4.7150","-1.0180","-1.0179","-1.0182"],["'MLO'","2025","20","'month'","'BbsG_S2S20'","'neph'","'lin'","'LMS'","1.2605","0","-0.0018","0.0011","-0.0046","-1.0851","0.6366","-2.8067","-1.0358","-1.0357","-1.0360"],["'MLO'","2022","10","'month'","'BaG_ae_psap_clap'","'abs'","'lin'","'LMS'","3.7449","95","-0.0090","-0.0042","-0.0138","-12.0274","-5.6040","-18.4507","-1.0902","-1.0901","-1.0903"],["'MLO'","2022","20","'month'","'BaG_ae_psap_clap'","'abs'","'lin'","'LMS'","0.0689","0","-1.2006e-04","0.0034","-0.0036","-0.1501","4.2058","-4.5059","-1.0024","-1.0022","-1.0026"],["'MLO'","2022","30","'month'","'BaG_ae_psap_clap'","'abs'","'lin'","'LMS'","3.1722","95","0.0028","0.0046","0.0011","4.0663","6.6300","1.5026","-0.9146","-0.9145","-0.9148"],["'MLO'","2025","10","'month'","'expS_bg'","'neph'","'lin'","'LMS'","0.4663","0","-0.0042","0.0139","-0.0223","-0.3296","1.0841","-1.7434","-1.0422","-1.0417","-1.0427"],["'MLO'","2025","20","'month'","'expS_bg'","'neph'","'lin'","'LMS'","1.0456","0","-0.0034","0.0031","-0.0099","-0.2657","0.2426","-0.7740","-1.0681","-1.0677","-1.0684"],["'MLO'","2025","30","'month'","'expS_bg'","'neph'","'lin'","'LMS'","1.4352","0","-0.0036","0.0014","-0.0086","-0.2738","0.1078","-0.6554","-1.1076","-1.1072","-1.1080"],["'MLO'","2022","10","'month'","'BbsFg'","'neph'","'lin'","'LMS'","0.4407","0","3.6164e-04","0.0020","-0.0013","0.2434","1.3479","-0.8612","-0.9964","-0.9963","-0.9964"],["'MLO'","2022","20","'month'","'BbsFg'","'neph'","'lin'","'LMS'","2.9793","95","8.7158e-04","0.0015","2.8649e-04","0.6228","1.0409","0.2047","-0.9826","-0.9825","-0.9826"]]}}
%---
