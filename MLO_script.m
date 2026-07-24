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
% homogeneisation MRI-MSE
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
MLO_rd.Ba_ae_psap_clap=NaN(size(MLO_rd.BaG0_A12_clap));
MLO_rd.Ba_ae_psap_clap(P1)=MLO_rd.Bac1_A81_ae31(P1).*cte_ae16_PSAP; 
MLO_rd.Ba_ae_psap_clap(P2)=MLO_rd.BaG0_A11_psap(P2);
MLO_rd.Ba_ae_psap_clap(P3)=MLO_rd.BaG0_A12_clap(P3);

% time series AE16+PSAP+CLAP until 2022+AE33
P4=timerange('2022-11-29',last_date);
MLO_rd.Ba_ae_psap_clap_ae33=NaN(size(MLO_rd.BaG0_A12_clap));
MLO_rd.Ba_ae_psap_clap_ae33(P1)=MLO_rd.Bac1_A81_ae31(P1).*cte_ae16_PSAP; 
MLO_rd.Ba_ae_psap_clap_ae33(P2)=MLO_rd.BaG0_A11_psap(P2);
MLO_rd.Ba_ae_psap_clap_ae33(P3)=MLO_rd.BaG0_A12_clap(P3);
MLO_rd.Ba_ae_psap_clap_ae33(P4)=MLO_rd.Ba3_A82_ae33(P4);

% time series AE16+PSAP+CLAP until 2015+AE33
P3=timerange('2013-01-01','2015-01-01');
P4=timerange('2015-01-01',last_date);
MLO_rd.Ba_ae_psap_clapS_ae33=NaN(size(MLO_rd.BaG0_A12_clap));
MLO_rd.Ba_ae_psap_clapS_ae33(P1)=MLO_rd.Bac1_A81_ae31(P1).*cte_ae16_PSAP; 
MLO_rd.Ba_ae_psap_clapS_ae33(P2)=MLO_rd.BaG0_A11_psap(P2);
MLO_rd.Ba_ae_psap_clapS_ae33(P3)=MLO_rd.BaG0_A12_clap(P3);
MLO_rd.Ba_ae_psap_clapS_ae33(P4)=MLO_rd.Ba3_A82_ae33(P4);

names_abs_tr={'Ba_ae_psap_clap', 'Ba_ae_psap_clap_ae33','Ba_ae_psap_clapS_ae33'};
break_MLO_abs=change_point_analysis_Def(MLO_rd,names_abs_tr,0.05,'MLO','abs three homogeneisation'); %[output:83915ba8] %[output:1f3e4364] %[output:1a8192d8] %[output:3e9ce198]
T_MLO_abs=make_table_breakpoints_def(break_MLO_abs); %[output:1417b02e] %[output:51adbc5f] %[output:4e38aa20] %[output:07abaa44] %[output:278538d2] %[output:1956b86f] %[output:809982d3] %[output:10204c9c] %[output:17960a08] %[output:52d36b7e] %[output:65eca5af] %[output:3fd1ff8e] %[output:1a435c1c] %[output:51a3da26]

%%
% AbsorP_mri_mseion: three times series needed. but much longer
% 1) for AE16+ AE31 + AE33: 1991-2025 AE31 is too noisy and no
% homogenisation with AE33 possible --> keep AE16+PSAP+CLAP and AE33 for
% 2015-2025

%%
lambdaSC=[450;550;700]*ones(1,4);
lambdaAE=[467;530;660]*ones(1,4);
lambdaAE7=[370 470 520 590 660 880 950];

names_neph={'BsB_S1S10', 'BsG_S2S20', 'BsR_S3S30','BbsB_S1S10', 'BbsG_S2S20', 'BbsR_S3S30'};
MLO_expSC=compute_exp_D(MLO_rd,names_neph,lambdaSC);
names_abs_tr={'Ba_ae_psap_clap','Ba3_A82_ae33'};

names_abs_ae={ 'Ba1_A82_ae33','Ba2_A82_ae33','Ba3_A82_ae33','Ba4_A82_ae33','Ba5_A82_ae33','Ba6_A82_ae33','Ba7_A82_ae33', };
MLO_expA=compute_exp_D(MLO_rd,names_abs_ae,lambdaAE7); %[output:1abda206] %[output:1bd90c64] %[output:67262dce] %[output:8779b381] %[output:23cb9a2e] %[output:0ff9556a] %[output:51721cb1] %[output:0b5f4ec4] %[output:56c2e67b] %[output:9eb5aeb1] %[output:034c3c7d] %[output:50631bc8] %[output:65cb25dc] %[output:061681c4] %[output:12cf2d52] %[output:4c6eea10] %[output:42042d8a] %[output:9b831cf5] %[output:62a07a15] %[output:6a36301f] %[output:9892eba9] %[output:9dae0082] %[output:5aae826e] %[output:9c710b4a] %[output:8fd68342] %[output:171979cf] %[output:61206dc3] %[output:590fdeac] %[output:107586d7] %[output:885de158] %[output:7a329501] %[output:9b556a7d] %[output:004c3265] %[output:781eb837]
MLO_exp=synchronize(MLO_expSC,MLO_expA);
MLO_SSA1=compute_SSA_D(MLO_rd,{'BsG_S2S20'},{'Ba_ae_psap_clap'});
MLO_SSA2=compute_SSA_D(MLO_rd,{'BsG_S2S20'},{'Ba3_A82_ae33'});
MLO_SSA=synchronize(MLO_SSA1,MLO_SSA2);
MLO_cal=synchronize(MLO_exp,MLO_SSA);
clear MLO_expSC MLO_expA MLO_exp MLO_SSA MLO_SSA1 MLO_SSA2;
plotFigControl_cal(MLO_cal, MLO_st.name); %[output:8e104d42] %[output:48ab55ad] %[output:741b5124]
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
names_abs_ae={ 'Ba1_A82_ae33','Ba2_A82_ae33','Ba3_A82_ae33','Ba4_A82_ae33','Ba5_A82_ae33','Ba6_A82_ae33','Ba7_A82_ae33', };
MLO_expA_tr=compute_exp_D(MLO_tr,names_abs_ae,lambdaAE7); %[output:7daedd72] %[output:23a2f8e6] %[output:6757de8d] %[output:297a0690] %[output:6ff3f15a] %[output:77a8cec1] %[output:803ba2bd] %[output:4225465a] %[output:43feaaec] %[output:841cf8f7] %[output:713e3948] %[output:6fce8c1e] %[output:015bbced] %[output:38dc8002] %[output:1251c42f] %[output:7fa8891f] %[output:30485f3c] %[output:65c6d840] %[output:155410f8] %[output:4570c4f6] %[output:919c996f] %[output:5bf15f4b] %[output:67838364] %[output:1f67c1c0] %[output:601b34ec] %[output:1757681e] %[output:6f129a6f] %[output:883f1fd6] %[output:78e033c7] %[output:59287da2] %[output:1ce6a2da] %[output:0b37b71a]
MLO_cal_tr=synchronize(MLO_expSC_tr,MLO_expA_tr);
MLO_cal_tr.SSA0G=MLO_tr.BsG_S2S20./(MLO_tr.BsG_S2S20+MLO_tr.Ba_ae_psap_clap);
%if SSAAE: MLO_cal2.SSA0B=MLO_tr.BsB_S1S10./(MLO_tr.BsB_S1S10+MLO_tr.Ba_ae_psap_clap);
MLO_cal_tr.SSA0AE=MLO_tr.BsG_S2S20./(MLO_tr.BsG_S2S20+MLO_tr.Ba3_A82_ae33); % if SSAAE take 470nm

%expS  trend on BG and expA fit from AE33 only:
MLO_cal_tr.expS_br=[];
MLO_cal_tr.expS_gr=[];
% MLO_cal2.expS_br1=[];
% MLO_cal2.expS_gr1=[];

MLO_cal_tr.expA_bgAE=[];
MLO_cal_tr.expA_grAE=[];
MLO_cal_tr.expA_brAE=[];

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
clear MLO_expA_tr MLO_expSC_tr MLO_cal2 MLO_rd_old;
%%
[MLO_result_MK,MLO_result_LMSlog,MLO_result_LMSlin]=all_trend_STN(MLO_tr,MLO_st); %[output:7daca7d3] %[output:691d70ae] %[output:56b72ec6] %[output:4cbd50fa] %[output:8767859a] %[output:11b066ec] %[output:009c4dd0] %[output:29858109] %[output:8c609ceb] %[output:43878f87] %[output:92137970] %[output:0a50c7ad] %[output:9533c316] %[output:96569083] %[output:3f72c21c] %[output:8a64f282] %[output:819da75a] %[output:5ad1ea29] %[output:0ac26fcc] %[output:21f31491] %[output:47f7eba9] %[output:184e0113] %[output:8b5c357c] %[output:9dfc1fb3] %[output:72c69771] %[output:92bbcedd] %[output:40aa086e] %[output:18c8d889] %[output:53b4c1d4] %[output:77ee198a] %[output:6577f61e] %[output:21cf559b] %[output:6ba35594] %[output:914fd866] %[output:4a8a6645] %[output:94b07b84] %[output:16cee0e4] %[output:07b88436] %[output:94d83b3f] %[output:7261d18d] %[output:6c4e33ba] %[output:44f1c894] %[output:43e6e121] %[output:5ca5a6ff] %[output:32a14f56] %[output:71e76d47] %[output:7ab5e6f4] %[output:34001b96] %[output:5f3e610b] %[output:106c22b1] %[output:782857ca] %[output:41d65df9] %[output:1891fe05] %[output:38c80f41] %[output:2b079263] %[output:0874b400] %[output:881ddcf0] %[output:44031ae0] %[output:34396911] %[output:27c4bfa6] %[output:4a0ad7c1] %[output:65e8ab02] %[output:6d460f76] %[output:164e8471] %[output:50dc1af8] %[output:5f574c25] %[output:450a624e] %[output:74e66a45] %[output:6a5fabc0] %[output:7675f741] %[output:24aa0d05] %[output:57521f4e] %[output:435b04ff] %[output:5834e47b] %[output:51c0283c] %[output:03ea7eb2] %[output:8cb754ed] %[output:2521d800] %[output:96422d9e] %[output:6b93ccef] %[output:1ff1e6a5] %[output:5a750aa8] %[output:3a652dd5] %[output:6a45a07e] %[output:440cb725] %[output:03527de8] %[output:62b980ac] %[output:145cb312] %[output:1fb44726] %[output:9a72e9db] %[output:970d9f46] %[output:627b9ef2] %[output:7bb8a256] %[output:047d3508] %[output:53258313] %[output:00eb0146] %[output:5f73b0d6] %[output:61911cd4] %[output:9d41c5c9] %[output:0e4470d9] %[output:4e2b6092] %[output:6f127a96] %[output:4b83007b] %[output:0484a902] %[output:11101af8] %[output:37495a31] %[output:3c309555] %[output:74ccd857] %[output:14be6490] %[output:65bbdc35] %[output:196cb9c7] %[output:6e09b832] %[output:51081c51] %[output:5d738641] %[output:8f359d5a] %[output:4481adb1] %[output:40ddc0bb] %[output:547dbc7b] %[output:3725f127] %[output:75f627c0] %[output:4ccb3eb6] %[output:5df4caf1] %[output:03d4c46f] %[output:5d7ad6d2] %[output:27351042] %[output:4a092aea] %[output:1540ae5a] %[output:81dcf9e7] %[output:03814767] %[output:3c0856bc] %[output:644bdff9] %[output:49c60f87] %[output:08f093fe] %[output:6afaeda3] %[output:1537097e] %[output:50a0bede] %[output:3bd90208] %[output:6bb0d786] %[output:4925dd53] %[output:7f550899] %[output:77943ab5] %[output:5744ba15] %[output:8431292a] %[output:1c856963] %[output:5c9fe296] %[output:11153c03] %[output:2693967a] %[output:39a6b9c9] %[output:64f431c9] %[output:72d925c7] %[output:9e3df2f5] %[output:173985e8] %[output:2c0fdeaf] %[output:3ce2d07c] %[output:21497ac0] %[output:02bba6a4] %[output:5ee0f6eb] %[output:1f63a72f] %[output:96ec5789] %[output:9a056215] %[output:7cca7a2c] %[output:4549027b] %[output:4115420d] %[output:8d4861a9] %[output:0931055f] %[output:350777ed] %[output:1a4e8313] %[output:1ee5317e] %[output:84db3364] %[output:10d3ba1e] %[output:61806ca5] %[output:63e1ac1d] %[output:7aac0809] %[output:2fba3f96] %[output:2dd0c217] %[output:0a526f37] %[output:1e8ae956] %[output:1b27364e] %[output:6ba2f081] %[output:3218c0d2] %[output:42a6b02f] %[output:05cb2d68] %[output:29988949] %[output:0ce840fe] %[output:260a6be2] %[output:70604c23] %[output:5eda84ef] %[output:5135b724] %[output:69bd9765] %[output:87e0606b] %[output:8a8b83fc] %[output:9683eb8d] %[output:59ed2d67] %[output:074657ca] %[output:5763ca47] %[output:3285274e] %[output:01434088] %[output:0daecedc] %[output:4fdfe614] %[output:1bd8ff71] %[output:15ced2fa] %[output:990f2ee7] %[output:3937d464] %[output:029880fb] %[output:62b5eea5] %[output:2fcb0a8f] %[output:8712811b] %[output:1ceb2ef3] %[output:3f3ea777] %[output:175f3b96] %[output:93ea4fbc] %[output:8c1d9f74] %[output:2334239c] %[output:37456ea4] %[output:06cd1b08] %[output:8fa8128b] %[output:3241aa6e] %[output:36208924] %[output:97346451] %[output:63a05afe] %[output:9453c2ec] %[output:9362a015] %[output:18fb473c] %[output:1c95b9b9] %[output:0bcd305f] %[output:4b002ed1] %[output:7a09ab26] %[output:206afafe] %[output:1f8b5d75] %[output:5648a775] %[output:35eea1c6] %[output:4b8031c0] %[output:6863686a] %[output:81d0f071] %[output:92100cbc] %[output:0acc5e8e] %[output:31975bb4] %[output:230e981d] %[output:6f869c7e] %[output:7a886ec3] %[output:0ef3773d] %[output:2ab8c7fa] %[output:676164f7] %[output:455d64b0] %[output:2f3064ba] %[output:82c04625] %[output:1633a5b9] %[output:00329b5d] %[output:4e26f982] %[output:123b00ff] %[output:5c17d0a0] %[output:3742a577] %[output:93505a75] %[output:61120f9c] %[output:6a5714b3] %[output:33cd7afd] %[output:1884a694] %[output:1db9a7de] %[output:57f1e15f] %[output:9580b661] %[output:8e0e0ff8] %[output:52ac184d] %[output:93cb677e] %[output:40d86b5e] %[output:8f043f03] %[output:9ec8596d] %[output:17b5ad5e] %[output:876a8eb7] %[output:9d3b11fe] %[output:2551645f] %[output:3e462a53] %[output:4603b234] %[output:2cdbab0d] %[output:6ba43a7a] %[output:292005b6] %[output:44587408] %[output:2d2c25dc] %[output:9a4e322b] %[output:7209c994] %[output:4aaaea36] %[output:15765720] %[output:3cccc1c3] %[output:36b97662] %[output:47d109dd] %[output:6d6aed90] %[output:719ba86f] %[output:60041ab7] %[output:8288b7ce] %[output:47a2a2cf] %[output:2fe9d7c4] %[output:3e1bed74] %[output:6fd93c27] %[output:4e7457ff] %[output:458e336e] %[output:0cb08232] %[output:553ea746] %[output:1c44b927] %[output:1c30507b] %[output:759a6b80] %[output:5c9334e8] %[output:0587c776] %[output:3971c5ce] %[output:04efeab9] %[output:013402e7] %[output:967fbaed] %[output:3116fd48] %[output:906525eb] %[output:4b070030] %[output:89ad159d] %[output:1ad2fd08] %[output:750a3223] %[output:57d51fee] %[output:1b4f49ad] %[output:19cd64b9] %[output:6032e7cb] %[output:32c1fec6] %[output:887ba1af] %[output:472fba2e] %[output:7ed1c1f8] %[output:3e50198a] %[output:2d5d30b0] %[output:783e6de4] %[output:7a094f46] %[output:1ec23811] %[output:7050df14] %[output:7da68526] %[output:80d45bc7] %[output:66e98f44] %[output:691c68b4] %[output:1847470b] %[output:5b01f108] %[output:77ae7dc0] %[output:763e72c7] %[output:41e0e013] %[output:4767e172] %[output:091c8eab] %[output:0375b2ad] %[output:308a38e3] %[output:57970a44] %[output:0a9eaadd] %[output:83eefd76] %[output:243b56a7] %[output:7d25fdd6] %[output:62224dce] %[output:71005745] %[output:83d93634] %[output:9c7cc22d] %[output:64150d21] %[output:728da710] %[output:1190d817] %[output:381adbb8] %[output:85b092d7] %[output:1ae00d2a] %[output:9b025613] %[output:46dc88b7] %[output:8dec3a57] %[output:952f7d45] %[output:11939650] %[output:11838473] %[output:187cb860] %[output:4f7ac2c3] %[output:1d78b69b] %[output:0f71f93d] %[output:5e8f9139] %[output:98ce6d35] %[output:2ef1c20e] %[output:76107dbb] %[output:07c20a40] %[output:1b7752e3] %[output:730d9a8b] %[output:1f45c1a0] %[output:8a727372] %[output:62e057cd] %[output:378fafb1] %[output:45651540] %[output:18518a12] %[output:500ad4dc] %[output:26dd33b3] %[output:57de100a] %[output:4aaa3564] %[output:55bda12e] %[output:0bc04a08] %[output:90cfd8ca] %[output:5e1643bf] %[output:3ff75fa2] %[output:8036573e] %[output:385019d6] %[output:0c1aa960] %[output:4ba62846] %[output:6937e54a] %[output:8174113e] %[output:3fd50b0f] %[output:5413a336] %[output:3fbb82da] %[output:4f138eaf] %[output:26350180] %[output:5cf78d39] %[output:2297ef76] %[output:4db7e9e6] %[output:6216ff5b] %[output:3372c653] %[output:2684268e] %[output:58083934] %[output:2698c44d] %[output:1995cd3e] %[output:01fdb1b3] %[output:57ac2e2e] %[output:76244fe9] %[output:87ecc18b] %[output:35966a80] %[output:953d2ac2] %[output:31a1500c] %[output:99b1044f] %[output:019ba41b] %[output:6830611f] %[output:14818eb6] %[output:793e416e] %[output:57b48bc0]

writetable(MLO_result_MK,'MLO_res_MK.txt'); %, 'delimiter',',' ) %[output:6b7c2cfc]
writetable(MLO_result_LMSlog,'MLO_res_LMSlog.txt'); 
writetable(MLO_result_LMSlin,'MLO_res_LMSlin.txt'); 
plot_10y_in_two(MLO_result_MK, MLO_st,'y');

%%


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":20.7}
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
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap\nTaille de la serie    : 352\nStatistique T_max     : 13.4869\np-valeur (bootstrap)  : 0.0260\nPoint de rupture      : 2017-12-01 00:00:00  (indice 293)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap\nTaille de la serie    : 292\nStatistique T_max     : 18.5135\np-valeur (bootstrap)  : 0.0020\nPoint de rupture      : 2008-10-01 00:00:00  (indice 183)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap\nTaille de la serie    : 59\nStatistique T_max     : 4.1240\np-valeur (bootstrap)  : 0.5940\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap\nTaille de la serie    : 182\nStatistique T_max     : 3.0544\np-valeur (bootstrap)  : 0.8800\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap\nTaille de la serie    : 109\nStatistique T_max     : 9.7080\np-valeur (bootstrap)  : 0.0810\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:1f3e4364]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:1a8192d8]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap_ae33\nTaille de la serie    : 385\nStatistique T_max     : 55.3086\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2022-11-01 00:00:00  (indice 352)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap_ae33\nTaille de la serie    : 351\nStatistique T_max     : 49.9353\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2015-03-01 00:00:00  (indice 260)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap_ae33\nTaille de la serie    : 33\nStatistique T_max     : 5.5399\np-valeur (bootstrap)  : 0.3560\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap_ae33\nTaille de la serie    : 259\nStatistique T_max     : 6.7403\np-valeur (bootstrap)  : 0.3310\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap_ae33\nTaille de la serie    : 91\nStatistique T_max     : 18.8748\np-valeur (bootstrap)  : 0.0010\nPoint de rupture      : 2017-12-01 00:00:00  (indice 33)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap_ae33\nTaille de la serie    : 32\nStatistique T_max     : 2.0951\np-valeur (bootstrap)  : 0.8930\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clap_ae33\nTaille de la serie    : 58\nStatistique T_max     : 7.5121\np-valeur (bootstrap)  : 0.1690\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 377\nStatistique T_max     : 13.0027\np-valeur (bootstrap)  : 0.0290\nPoint de rupture      : 2015-02-01 00:00:00  (indice 259)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 258\nStatistique T_max     : 22.9387\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2001-12-01 00:00:00  (indice 110)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 118\nStatistique T_max     : 33.2074\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2016-05-01 00:00:00  (indice 15)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 109\nStatistique T_max     : 12.0046\np-valeur (bootstrap)  : 0.0410\nPoint de rupture      : 1994-01-01 00:00:00  (indice 26)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 148\nStatistique T_max     : 5.0383\np-valeur (bootstrap)  : 0.5430\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 14\nStatistique T_max     : 6.7535\np-valeur (bootstrap)  : 0.1440\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 103\nStatistique T_max     : 11.3077\np-valeur (bootstrap)  : 0.0540\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 25\nStatistique T_max     : 3.1008\np-valeur (bootstrap)  : 0.6870\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba_ae_psap_clapS_ae33\nTaille de la serie    : 83\nStatistique T_max     : 4.6158\np-valeur (bootstrap)  : 0.5330\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:3e9ce198]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQl4U1X+9lsabNj3TQpNEQQUhHEZQNC2CirKuDGKotBWREEZ+RgXBNS2KigiCuhfdEDbiiIuzKjjhoJNFRRmHGVcGaiSQqHsoGzFBvr1PelJTy735t4shaTN8eERkrMvb377L66ioqICIZblP+3G9Qu+9fYy46quGJfSKcReY83r6g4c3bMF22aPQPu7X0d8y451dRuiYt1x4QCQqFhpbJJRtQPbnroBCZ3OQIvrHoqqede1yYYMIFrqQ25gt7YN8e6dfdGuSUJd29PYekPcAVIgW6adD\/euTcf1ZGvdGR2nfxGjTELc43A1DwlAtu8\/giv\/by2mX9UVg3u2CtecYv3EdiC2A1GyAyEDyF1L1mHeDT1ilEaUHHg0TNNIBhKTjUTe6YUEIFzO84WbsWHHQcy+rkfkrS42o6jaAX+si1xIwukD0WHaB6jXsGlUra22TjYkAJEszIYdh47bn5gMpLZemZpfV4zSqPk9DtcIIQFIuCYR6ye2A7EdiM4diAFIdJ5brZ\/13jcfxqG1Hwt25dC3n2DH7D+jXsNm6DBtGRJO71fr1x8tCwwZQA7+7sYt+T\/gk5\/2YEjPlph5zekYt\/hHzPpzd\/Q+tXG07IPpPFevXo2RI0eKeo899hhGjBgh\/v76669jypQp4u+LFy9Gq1atkJmZib\/85S\/eOtrO1Tb87k9\/+hMef\/xxNGjQwHQedaGCysJwvVsfuQxtxy8US98x\/1ac+uBHMTVuhFyEkABEgkdi8wRkDEzEvOXFmHNDdyz6shQF6\/fgpfQz0egUW4QsNbRpqACiPviZM2fihRdesAwgrP\/Pf\/4Tubm56NatGw4fPoz7778fmzZtwosvvoiWLVuGNtFa0JoAIkGjfE8Jfv3n04ISKS\/5KQYgEXa+IQEIhahSjbtjf7kXQA4cOer9vLYYkkkAueCCC\/Dbb7+Jx84yZswYNG3aFJ9\/\/rkpBbJhwwZd6sTo8wi7Kyd0OgdWL\/VhWwgkZGPa3v0WGvcffkLnEhvMeAf8AogVK9O731yHkn1HMGlwEl78fAumXd4FIxZ+i0GnNatVql0JIGRjnE4nnnzySbGr99xzjwCRRx55xBRAJOtCVqd\/\/\/7eU9mzZ4\/og59Nnjw5dl9jOxA1O2AIIIFYmdYFZzoJIA8++CDWrl2LU089FQ6HA0uWLMHYsWMxYcIESwDyzDPPeNkXeUskG8N\/x2QhQMyQLGrwA34BJGZlWn2QEkAoQHW5XFi3bp34skePHkhJSRECVjMhaowC8f8wYoZk0QMccqZ+WRgzK1N\/hmQcoEWj+nh7fN9aoY1RASQpKclHIyP\/bQYgMRmItQcSMySztk+RUMuUhTGzMtUDGfnZ0F5tsGBlSa3QxqgAcvbZZwthKAu1Kbt377ZEgbC+1NpIOYiUf\/C7mBYmEp5EbA6B7EDYtDCqtkVqZx66ogsefv+XWuFspwLIlVde6aN6Xb9+\/XEAsnXrVp9zuP32270C0lDtQFxwib4dcARy1lFV99ih31A6\/XIcWb\/KO++YH0zkHaFfAFGNxOTUaSwm7Tvk919t2u9lVb7begBXz1+Lczs3wch+HbB4TWmtoEBO9tERNPKRD2fVf3I+BJEMZCALWeKj2gAuEjziWyWi\/V+XeLeeQYaO7i6JOdOd7MuojG8IIKqRmOppK9W2qpGYVgvzxtiz0DuxMdJzv\/exSNUCUuaADrVK1VtT50rQyETmceAgwUJSIwQT1jUCl3DOT9ppCFnX9Tm6kcP06mgpC70AQTEtTDhPqmb7ClgLoxqPBWokpspLpAA2c8Cp3vipfBDqAzBbeipSazUZz\/VzP9KQJraC6yWlwf\/Lkoc8AS4qaKiUiASXXOT6tDPbW3\/fq5airKdnXm5Uh\/Wllak\/nxZSG2X\/\/cjr+3Jk\/RqUTr8U9j6X+VAloawj1jb0HfDLwpDaePvbXcexJ1ef1ToslAP779a2kRdAtI\/BbHl8FCTfa2shoCYjWSwvG9leNkWuVwUXUh+SIqmAJ042\/01wkaBcgIKwgAgpiz35fxWhBePsTYSsotmfJvlYiBrVqd8y0bI5Oh3q9r7hYc1YjCid2nr+0bAuUyEqqYap7xR516KNuG4mJzHaBFIgWhZHAsihcrfp3jWsb8OCoy\/hpmOjTetGa4VX672MsfG3YODRC\/HJsU99lrGl3iZ0i+8iPptc\/qAAlyH1LsKq+M9w09HRWHDsJW\/9sfVuwavxL4t\/W9nb4\/Zr9ybUj4+DrY1HaEtwkP4p\/DcBpGHfS3zYGKM6rH+iQaGkpAT8E20lMTER\/BPJxRRA\/E0+EDmJ2o9sN3ZQojeW6qY9ZRj78xx8fN4UFO08qDus2+3GoUOH0LBhQ\/To0Awz97+A4WU3W9pfta3NFriDXyjt\/bW9qfllKKm3CYnHOuPVfR\/5rIXfran\/Of62dSkujBsMzlvWZ8Ut8cXoV36Bt91S+yuY3OR20Qf3UB03vfUw0de1ZTfjif0e5z+zIts3KpiDxiVr0PbOXAEioQCIGmVdGo61TH\/qOP+WcGlhCBz33nsv1qxZY7bciPu+X79+mDVrVkSDSEgAYiQP8Scn0aM8eHIri\/bhsn8\/joM3zjT8lSw7XIZt27ejQ\/v2aNGkUUAUiNo2wR54pPhQ2vtre0b9roLVIAvyQfkKFCpC0PH1x4hL\/UjpTDRu3Bh2ux0z6z\/iZVX43SdHP8XAYxd6Lz8pM0lpqOO+1fB1Qc3I76y8Ftm++YZlOJx7BzpkF6DBmakCQIJlYVRHOAkSWuolnFoYqX7nQ+zYMXpyzBDw5s6dK6ybVb8pK+d2IuscByCqDcctL\/8IM0OyQOQk\/oBFBRDJw2s3oqysDKWlpejQoQMa2BsgEBmI2pYPMdASSnt\/bSkglbILsiGqQNR0ji4HsvM3+lTLTk8GHC5k521EeVEiDhzYj2HDmmDwYBviECfqGu2v0X433\/M\/7H18qBdAQhGi\/rbC48VMSkQIRmdfhw53v+kTJCicWhgJIHoPUQrtC1EoADwd6REjlPc3b9N7cQIrhESByHmayUlkPYJN7pelPsuTMpW6CiDqZkgZkHzg8sH\/svEXAZozZ9qRnU0LMhewMRmoBBAk+wIIKjwggbjqhINsk5UFIZDlo9mIjZYeigQ+LYCwe1VFK13stQ9fr46WNTESjIZLC6P3EKVNDQXT2kKhPH+YTnaJegAJhj0JZdNjAAJoAUQ++ER3opB\/uHIqNU7ZWXBUeIBAS01IrQx\/TQkSBIAGDTygQwDxR4FISohtqa1h4WdF7iIk7bXh5TtcXgoklHMOpG04tDB6DzEHOUKrxbVSJZ6CFBSjWHzGooKI9F9SLYtVq2LteuiaQG\/tiy66CP\/973+FvI6uDsXFxSJC3aeffoq+ffuKZrKe3p5ELYCYOchxsao1qpWYIdoN0hOixgDkeACRF927fzmVFzw7C14qAxCPXdqF8METRKTKVwWQlCxfexIJErJvCVYSfPi5F8B2AoV34IQDSCBgY1RX+xClalwCpeoOoKrNJZVGAFm2bJkI18DC0AvTp09Henq6iCinLRxPgsWzzz6LSy+91FuPbefNmyfCP9D9QdarVQAiF2PFYCyQmCFqv8xmR9kKLVZlRrsYgBwPINwzSTWI\/asCEJUCIXgQDNRfVf6bD0MFEGeWB1xY+J20G5EPSZXF8PGwXOC+wEPp7C2pNQAi98mIVaEMipSgBGF\/ACKdKLlXpEoIDAwMRcqDcWMYZKpPnz4YP3489u7dK\/aU8XMvueQSIcv77rvvDIWkUUuBBIL6VkBG7Y+Ux0PvFOGWgR1BAa2aEjMGIPoAsty9HENsQzzb6EwFMnM98g+lqGBAQaw0rpMA4sjNgSvDQ7KrRWp\/JGBoz96fDCSQexJM3XCpcbUPUQKEkQBespESYPRYGMaEGTJkCHJycpCVlSXi2NJBkoWhHbQUiGRh6IRZZygQboZWOCovgpo0yixmiN7l0aNcVACxeuECMSQ7UnYEpdu2oX27drA3sOPy+oPFr2vHo52PM9LSG1\/b3uocWc9qW2k4php7sW3vBt1RYjM3hJp29CFMO\/aQd20c2+V0AKlOdDyaJIzBWLhm8V0l38+y4djPusvhHp0MGUhNqHGlFiZUCkRulJFshMGlYgACQLXXyFtVAsb2ILuhNT+3AjLa22kEIMOWfIjLrtmCsztXpy3cUs9zyTseS8KRI0dw6OAhEcQ43hYvjKgSjyZZestsu3PnTrRp0wYJCQlIaXmGMMTiwyrc86NpH9r22gYlm+OwZVMcSjbXw5bNcejYqQKJnY6h38BjwviruOIXJMV1Oc5YTO1HGoI5vz7k7at441H8vcd0bDu9yq09tdpZzoMCDo9WpqpwPVyXWviZLIGsWe5R4gmWgdSkGlc1\/9ea9qsyEPmdloWR+0hhqUqByM\/1ZCB1kgJR2ZN\/fL3Dm\/9W\/ZybpjVHN32JVeBEOYiWhbnyuW\/w93F90C55H\/KRJ3h2Wk+qj2PYnqvwuP0xBGoMpjXmkhQIL82P5T\/DAf9AZGQM5nIB+XnAzBnxuksnkGwZ8zCQ5WEhfiyvdguQDQgJhQBeqTRd30JjL2eF54M8D53gBQrZgCwMQSOH6sYMIN0FOPIrv\/VVS1J7UzL9VkFxOLLyUD++nld781x5kRhChSO66KULGYmnXBKXipW2leLv\/X7wCGxpSHYiCjUw+z99UfjbxLfsKOKkbpl2PppcNEbX89doTnqyBMnGsA1lHUlI8quFUYWo6jhqqg9+TtZGRqeToS+ZwmPSpEnCKliNI0O5yLhx43xyDOn1HXWGZHIRUlOSdnpLDOzWAmPyvseLGb2wasNe5H65Fe\/e2ReNE+KFTOO+y5IRiGeuEQVCALn7jl+xsOskU9f1QL1LVWOubfZtlmJrqAeqZwzmdAIMTEYQMS188Bn5wiUuJaWSq0j1QAP9aKsfMRGDn1CISUDgN8oTJ7WRnwHkZQEbyY5UAYjP4NWA4+Dc2Dw7RwAIZR3Hj3n8zAkgDmTiw7L56GbrJtgnUiHrd9QsgFiJiarn\/u9v742EkSqIqO31nBZNz7YGKtQKIaqWCpFOdar2hGrcD7\/fGZB3rhGAXL4sF\/vv\/Ks4Dj3X9XVl64SAcI3d49cQiHcpAeDfO\/+NH9r8gJn2mT7m4Nrzl0F6VMtELYAQNJJ9ZZmWrxHtMlKyaGdRVfiqiRckObKqXfPNO6QncopSjdSBJkqZE3AUst8cYR2iYYD8DuHx700W4MFS0wBivt7Aa5hZojJIE1dJWxBpFxL4KOFvUSsAxGxb\/NmMqIJWs374\/Wd71iGlZU9R1ehXQD7iBzs9iFdtr4q6gZhly19TK\/NhHVVSrwLI0G3j4UzL8sgf+GQrksXlU20rpA2F+J6u9mm5Hi0K\/53thCvL8\/cMF+BMq6JisnMEqxN8IaBIACFL5gl1QEqkIJd2HZ7CsEMeZ\/9iL+smBaqeGtVtPTDnPOEsTPB74NsyWh6idr3RMu+wmLKH47ClEZTWRoF9S3Xjh2Ufen1hhtqHChmJmemxat+ghvvT2kHINch5yH9LgFIBpOfMJXBlK3FIKuIESKjqUG0\/or+0Ag+I0NCzkoURzz2nEjBV83TxwAkunK2HHUl0l6PEx4PY86g9hWbp1f40voZRZHFSq6sKzMqDo9ISU91XzlsFPE+\/VZMU41QHNNIaoIXj7Guyj2h5iLUSQKhhkTluv\/z5V1y\/4NvjUjUEGg9EtVxVY4uolpDycms3dZB7EF7e\/LLwC1liX+J1POMjV02x1cfheYa+Qgr5yCQwqWCirS8BjX1SpWkr6QpXssfUm+pRAQhVAKL2w6+l4ZYQZlINy7o5BZ63SQ9cV6a3L0dBJlyplIF4mDeGAvLO28WXX1BNXPD9p6R5xndlAMm5HutU2ojkKcBGYKqOx1O1nWxcTeVI4DseQNTGHsBSQdJKSMOaBAarfRv5wviLfKeNdKfmP2ZCMZnXmJ+zyGyCUjNz33334ZVXXhGGZTLXMbU5+fn56N69u7AdUYvap\/w8WoDPUkhDLkpqW\/j3e9\/6H\/Izewkh6i35P4DJtc3iprKdqhpW+zncpNQbecvsYki\/EBVkeOB6QGHUlxq9SwsmeuCl1hEgQEqiSiiqNS03mz9yKqoedR7gzK\/uy2scRnShfkRhZTimq0ByJJVARKogxwMgLJxPQZoHPIoVGYiU1vpMyhdA+JW6H9VVVQDxCHclgFjxxjXdB5MKqiEZo7G3Hb8A22Zfj7bjF\/p47pqNo\/cQzSLfqawr2xcWFnpBQgLBtGnTRBqOFStWiO\/ocm8FQNiuQYMGou6CBQtw1113iX9rS60AEAka23894s3v8suuMi+AcNF62euMLFRJfUx7p8irwSH4MKhQ1577LQOI2YWp8e+lPwrlFSwa3xTT8amiFRxKZqX+t\/Kxs72P7ENqYBQA4ZiUuno1qBoQIHBkCJ2vpuiRIMcDiP6cjQHESjwQ033wU0E1JGsxbBL2vJmDdpOW4Nf35+DQ2o8DisruD0D05Gd0HVABRAUM7UOnr8u5556L1157zUtV0DbEHwVSZwCE5yvZDZlhjkBCNkbVwgQSD4T9yURT7J8AQjXxuJROvj4fody+mm4bLgBx5lQRGsEACMGCQFDFmuVlV1InQpCirJ7f6QGIVBV7aA9ZbFUh\/9yJMpxkddvEnS7Y2njkIJSXWIlIFsoxqIZk7l0lXgCpKNuPbbNHoP3drwvbECslVADhGKrVqcpuSGc5GoqRSiHLIgHkiSeeAG1A1PKnP\/3Jm\/u41lMgVg5H1rEaD8QfgGhlIGbjqwJRs7r8XpWDaH1C9Nob1XflefxRHBlOuByVrEaVe72VOYh5ODcKSoKCTFeeq7qvXKnUpVDUw6J450AKw1GpxamiQOSzd1UFCHLkVcCVQUGuahdCwSob+PLbDmGdRp2xlAtVxxN5bfVqTOk\/pYodpCzF05a9fCJAxHpMVKv7YVRP5oBpeV0W9r3\/NFqnz0bp9MtwSrcBAUVlDweAqHNULVBp5CW9bSkP6d27t\/DcrdMUSKARyQK5KEYsDE3kpY+C9C4161frum5WXxUSWlH9Si9YrXYlL89jPOZwVKpgGZ6jMtZGRXXsHrNpII7vlxwH5aBKX3SQ8wBGAbKFGLXauS43r8KTuKFKPkraIEWme3A5UODaiLRUX9Kb1fWYFV25qs6s1bbaNjXNwsjpqIJafhZMVPZQAUQ6yTGWB4vqkq8CiJqilJSJkRC1TrEwpq8hwApGQlRpxSofrZlqVs913WwqamyNQABEzyaFBmQ0JKM1KeWZVgFEGJ8RCaqUOMQdb18FOZXyWT7VVGSjwAdANqICycSWKtKDtRjdg5oEB7KxEVmCBdR6mIYLQKRCV+7xiRCimp2n1e9DBRCOo2ph+G9pXq6N98GxWDcGIFZPR5GTaJsYGZKpalxVlsL22iRKfBDagC8yzwk\/V13XrUw5XABFCoSUiCxWAITgwXZCgVL1IoXdaFVfAowKPN6yqchAQeXapbYguzLPi2peJtgf8YnHfI1Mjx6AEKvENKXJSOUYVQph0+2SbY3q64UrNO30JFQIVQtzEqYshox6LQwXodp4MArZzGtOx7jFP3rTVQYTUMjsQFQQYV0fFWpVYy14qMmqGcBFZrlXnZ3o3HR4xOGAghb784uQlAOnlJzcRQSOMRq3X78RyK8UOwiDMQA9h\/4LP33wxyqwAP536hGUliZ41LJUx1ZZg\/C5M76Y1cJonvElg7F06VK4ExNxYPhwIepIyqpAWpeHkVdlf9Bl0yaM7uxx6d9cUiLqswy++GL06tVLMFIEHhI8VaZoVqcQlnr+fGLC4QsjgykbTTYSMh5GPYCoOV8yBiZi3vJizLmhOxZ9Weo1Ljtw5KiuGjfUW2QU9Fb6qMhE0hyHEvIZM2Zg9uzZwmiHJCTjV1KQ9cADD2Dq1KliOqwz8emJeC3hNTzR4Am\/U+xY3hFt3m+D3NRc9G3uiV+pLa+8UoJRo4JL+tO8+VxMXzUEd55xhqdb6QtDSiGXziteLxmDefJpUyCi9SDW+MK4gM6djyEjpxANG\/4LH02ebNkXRgpqNSGbQz3akNpTsNr4\/OuOyyHjr9NoeYjaNUTLvC0Zku3YX+4FEC1oBONMF8gtClTbwo2nXv7aa6\/FnDlzBFVA\/f3999+PG2+8UQydlZuFvnP6HgckDXc0xIQmE4QvjqzvLycH2ZE0s7detdikpAp067YKZ5xxBr7++mnMmnUF2vfvj+G\/\/Yavm1bHP\/FFFA31oXrjmpEGLqDpVXuQl12Ma675g3fLiVUXV1Igv1RRIHpnQfAgl6Vxywvk2GqkrlGcECsAMnHiRDBRU7SULVu2iIRYUevOz42mjUfJviOYNDgJL36+BdMu74IRC7\/FOZ0b4z+bDujmjJEHFKgzXbgOlhSIw+EQcRkIJI8\/\/rjomoBw\/vnnH\/f5HU\/cgXPOOQfntvYYBGnrS+m7v\/nlSH8WnUpNm+7BJZdsxZtv9hLfSmm9tF4k4L2wbBk6Tp2Kz+Lj8csvv2BLSddqqkTGCnA5vJofciJizCqNjo9DLlmWggrc3MmN\/fvvEWtW18Dx7777bqRnZeEVmw0flZWhSxdPisxdX32F6w4eRI5VVAzXoVnsR9X+WLUDiWWms7i5QVYzdabTRl3X5sYNctwaaUZZiMvlEqbFkhIxAxAjYJGfWwEQOW6\/fpNB1ub777\/DmDEXo3PnCnzwwX0+j1gPQPSAq127ESikRbuIs+lRG1PQqi1ZWQX44ov6GDRoENzun\/Hjjx9UzuFWH9BUVZAqZWW0R1bWXCMHWNWpPxmIzEETyPix3LiB7FZgdU0BJLDuTl5tSXnIyy9VanosDOta\/dwsraDVcWU\/egCiNx+zcbnTgYwtKQ\/KhGQ6AqM9sjL2yTvp2MiRtAN+ASSYnC8nY3F8SAxmq1589cFwTlLQyr+ThFeFqxTAGn0uvSn11hXIuLIfLYAYzZP1\/fH8gYzNuevF7\/Q39sk4x9iY0bcDfoWo2rilkbg8bVxKzlH6HDA\/x8iRI8W0VWGU2sbK53rrDmZc9qMFEH5mNB8JIK0zn8Ku3L96fUACHVtrCMUxqdYmtWY09sk465oIaXgy1lGXxrSkhfEX79TMVqQubWZNrJWqy0NfenKO6JVA7SJqYo6xPuvuDvhlYcxyvlixFWl0iq3u7m6IK5cUSJsxc7HzxYkBeaGGOPRJb67NiysnFAPMk340PhMwZWGYglJbpIqWn8t4IP5sRSJrydEzGyMWJnpWENxMVV+bve897TUeIzWW0OmMgNI6BDeDWCurOxCyFsbIVmTQac0CitRudcJ1rV5dZGFU4fFvK17Ekc0\/Chf+YAzJ6tp9OdHrDRpA6NZMmwIGTWnZpj029r4DFU09QV5OlK1ItJj7hnqoVh6OPA+ORdsXvTB5ch6qS7o\/LVMo8w5kPtpxZESyhn0vQcM+l6J09nXocPebOPTfZT7JpkKZX6xteHbAL4D4E5DaD5YiMzNTaDxkUNnwTMl6L7UdQFQWZvvs6+Detem4zZEygd8btBSAbgVAqJXh3knnP+s7br1mKADCUbRUyN43PMGNgjEksz7rWM1AdyAuOTlZhMORaj3+XfVubdHpdMybvwB5X\/3mdab7eK0LRwtm47tvvxXjEUToe0JAuf322\/HCCy+I\/xNYVBWiVK+yjbzsTZo0ESpWhoqjbQbtNegMJ9trF6TOTY4lVbFSRUr1LYuqotWbhwxuO2bMGJi1Ub1tteOoe2f0nQyLd8MNN4hAvBzPaI2cu3yAmzZtMnzo6li0g2Fh3mBJgah7Jff+3XffxZQpU0Rd7vn8+fOxcOFC\/Pbbb+IzzpNRx7\/++mtvPX\/nZrQvevMJ9HLG6kf+DsTt3r27Yv369cJegg+uVatWAgjSx47Hh\/s64\/cVT6JHn3Pxa7drBYBIZ7pJ59lwz4TbcdnlV+C7Npfj08++QJNVs3DD6DH4uuVlwuX\/g0X\/5\/2l41bwodLYi5GoCSD\/+c9\/xGVl4ZhMfM1fxU8++URcXq0jkXyEkuqRoMB6vMjqLzAfyjPPPCP6Z8xK7fr+8pe\/eHOVyl9tbRtpYk5TaJXaUn\/BOVc5DvdOrlGa02vHpd8NH7g6lrQMVa+LFQDhPMhCqnso+5c2MNyb008\/3Wde6vylo6E8C85FpezUtkbnxr2kTYm\/+fhjqSL\/mcRmaLQDgoWRv2T8dTz77LPFYyEV0Hf4Xah\/etpxznQUkI7rG4+MjExUJPXDgGvG4pwGpXjgrjF4KX8R1h3rDEmlnD9ggJfF4S8iH5v81ZMPl4\/F6OGp1qWyvczLoV50CXzyMkuwefLJJ8Xa9QzKZB29NnJclXLRUlXqry\/H0LJUKhl\/6623ggmVJfgZsV9aCkYeXNtTjmLGaVvhaOgJenygaWc8tKMXBgy6UOyvlmWYN2+eD5siwxwQvNTvJIDIs+C\/tSyObPvQQw\/h4YcfhqSK1B8FBhRWzzAUFkYalMW37BRQBPbYMz\/xO+BlYTi0HtshLmvfdLgdF4rZSQEpH9\/o9AzEOfrjHy\/MwIrCL7wA0rPvubj1mU+w673H8P8m3uX1Bq1JAFFBQt1GshcyK7oaJZu\/zP7a8DtSQdwT7eNQo3SzngQSlRrhL7kegEiwMpPfyLYbf1yLnEYr4bgpy0d9+UveVEx4Mg8DLr4MD859yWcs+dC1UcHlPJmPRMpAtADC9UjhuLqPbMtQfYw2rgf84QQQ9q\/mhZHzaDhgREABlU\/8c6p7I\/oVospLrpc5Sz6iRt3OR\/Pzb8KlrXfi0btvw1Pzc\/HYN3ac1x7Y8s7jgmWRQtZwAIjK2vijQIyOUsoF+CAYcIhxIuSjlm20v55aCklbjw+VQEO+X7IsXHeoADL5nrvx85oVeH7GNHS6ZLTPkkipZIy4Bu3Kd2LOUifqNWziw8JpKRBh+IkFAAAgAElEQVS1sR4L448CsbIvEkA6d+7sExJB7TfU5xWMO3+oY8ba+9+B42Qg8lf5\/xbmI\/tfp2Bzwcuw7VyHgwMmoiKhibe35uU70Oarebj26itxzpW34sYZbwgZyP6B9+KRW4eJXC\/qRWVDrQxE75dMKztQWZhgZSAUCEo5RWJiovehyV9qOQ9VLvH3v\/\/d+ystZTIECSkHkG2kfIQgNGTIEFMZiFUKhP2bqW8fz3kA77y5BLl5eYhv1k6wnsHKQNSHbiYDMTq3cMpA9CiQeg2bocO0ZQFlposBQM3ugKkWhsOrwkxp3t6n6X5Mn3wnRl1\/DR6YOsUwCKyeDMHsF94fea86fwWrhVEpKiPNjTaZEAW88tdVUiRSc6NNGKSn1dHKW8xYGCsAcmDrz5g08k9YsalMUD8sZloYsizqHlKeRLDUUgp6GhzWUQXVWspMtQ0KRQsjZSAcr+P0LywnkarZpxLrXW8HAnamU\/PGME0lUyT87eYz4c\/hrqa3ngG09ILtGI3LUIQFVakVanpuofRvRoGYfR\/K2LG2sR2wsgOGACKNyL7atB9vj++L3qc2xndbD+Dq+WtxbucmGNmvAxavKUW7pvVxRvsmgmWRxZ8THk3fc78sFVVlykz2HUpheD9NwnO\/3QVaP5S5hdI25t7u2T061on7ct1DoWxnrG0N7ICpKbs2qBBzufRObAwm3r5\/aDImL93gGxu14hgQVw\/N7PFYPaWfD2UiQYkJtZmNLlwlUEAItH645hnrJ7gdiAFIcPt2IlqZAkggkyDl0bVtQ3z4\/U7RbPZ1PXyaq5npQqU61I4DBYRA6weyB7G64d+BGICEf0\/D1WPIAKIb9rDiGJo1qH8cBaKtmzmggw\/IuI9WoPzYsYDXNv2Repj2oP92l9cfLHLPdjzaGak5TtP6epOoX68ebPGe7HHRXKIhyHDFkUPY\/cpklG\/+3rvV9Tv1QqubZyIuoWE0b7\/luVNjyD+RXCw70+Ho72jwdS4StvzLu54hl12Bb069DjOGn+nDklDOoUeBqBsh2Zm001sK+cmmPWV46YstGNyjJdo0isehQ4fQsGFD2GzVAYncbrfu5089Fodb7\/zVb\/2L256FLfHF6Hg0CcOnbEDG7bst989xfy\/dgFY\/vYO2l9zizVLP9Uh1Y7M\/TfImPFLTPqoJodUgOapTmNHnNXFxojnNQU3sRyT3yTw2s2bNimgQMRWiJjZPwKPDkoT6jurPHV2uErli\/u\/6Lpg2ZSq+Kv4Vf8+bh6Q2zb1nYQYgEjyYlTqxhV1QISuL9uHK577B38f1QZ\/2Cdi2fTs6tG+PBHuCt9+yw2W6nz847QjGjN3mt\/6+K89EZpYLqS4HnMjCpUN3ovk7b8NWsgVwJAHpGSBQ7Nu3D\/ua70NWPY\/gjkWM+8U7SHjxJnTILkCDMz35FVQhpwQEo8TT7l0l2DH\/Vpz64Ec4vH41fv3n08JMu7zkJ93P6zXUJpsKz1WXKlxezI4dPeEXYiXydmDNmjWYO3du9CaWkuraeTf0QP3yg96o3uX1G3mjkPHz8XdPRZuUTDybUZ31ywhAyML8\/ZttIinV\/3YcRGJzO+Zc111QLxJA3r3jDzg30Y7S0lJ06NABdru9+iGXlel+fv\/9ZZWWoCb1e\/ZEToYLWXkO5GRsxOjRG9HpwgthKynxJF3ZuBFlVf0v67AM4+zjfMdd+Q9UvDDSCyCkPHY8fxtaXPlX0NW+ZfpTggJRrSXj7E1QOv1ykDop3\/wDDq39WIBGRdl+bH3kMrQdv1DEuND7POF03yxqKlWjXvdAQ\/xZsT+JvOdU92YULefkl4UhELz97S6hxj246XvMf2kRVra6Fpf2bi9AgFaq7kZtvX4y3mNWZCD8TIY9bJwQj8FPfYX\/7Tgsqg7p0QKv3+bJPRsOADlypANmzvQAjkzo1qfPXowb1wA9xg9FTroLWTmoBpDRoz0AwqIDIGlIE3KTRHciXnbe4QMgcq2SClEBRFIXrEMAYWAcFjWy1pZp5wvQIbDofU4wMhojlOcUjoup9QXifCI9BSMNBukDREviaPAMDsc5hXJPrLY1FaJSszLtjW\/Q6Mu5sO0zTrWs9Uz1NwFSN9qUERJAqCY+N7EBSrdtQ\/t27WBvUE2BHCk7ovv52MwKHDh0GD16NiQnIggK8WDLjuCrr3bho4\/aYekbCbipKAvTjj0ECl1vGrUR3W5JFwDi7tQZ5Ss+FfU57kftPsK4BrfjjPpdBYC0\/70jPnPOQvyCagrkRAPIttkjwhJUORwXkwCybNkyTJgwQWwD\/52fn49p06ZF7OOMAYhVSAisnimABNadtdr+AOT5m3qiV7sE7Ny5U8QHSUioloEcOXLkuM9LNsdh1qPHMGrMDvTu09yw\/s4ddqT+8wm82mcy1nwRj2tv3Ih+wy6CfesWMenD51+Aja+\/I\/r\/rNMyjKwYg5uaX4aSeptw9KhbgMiuw\/9C1zapKBCppz1FjwLZk\/9XYYIdThaGglZJqVjbZf1aNQUgBJHLL79cgIkMB+GPMlFN5dWATOqs9dwg1HbyR4ttnnrqKTAhdWFhoTcvkEppqACiBlWSfTAmznvvvQen0ylCWRjNKZS9D6RtOM4pkPGCrWs5M53N9Rkar80\/bpxGHbrhrdfy0b1ze+93qvxEz8TdH4AEIkQlm\/LEjHi0anPIVIhKYezYeg9jy605KP+5Ah+VDUDjf6+pXo\/Dgb3\/\/V4IaVd0XIYx9W9D\/YsvAlzFKE\/siNM\/24QSWwkccGAjPJQYAx4f\/vJ1iJBupzREx6xPYWudiJIp5+Pobk\/4wXpN26LTk1+DQtTSx69Exf4doMK5foceSHx8jSUhqj+L1JMhA9FjYWQoCPoHMQiR9Lnhg9aGvNRSLEZZ9rQsB31t+ND56FlkO44nsw0yjIK\/\/hiXhZoobR\/sTwaQ4jh6mfyCfWTBtIt6AFEf+dnt4rwHRMers8473+upq8pAjqEC9eCxk2jd2IaV952Hto2rWRB+TrnKi19uFfUanRKPD+46W5jJByIDGd1ptHjMrrxUZKemYPkDg3BO61Zo3LgJ6tevVvuWl7txoPV+7BrWBFtsNjhdOUChJ7amI8UFR6ELjmyX4Hkq\/yb+T5mJfdzreLzHOCA5WQhT3ImJGLQG2FK\/RKhwCSBSfdskLR37lj7qFaLysW+eej6OSQBp1RmdZnzhAZCZV+LYbzuEHUP9tl1x6qOfg9oWqcY9Ed6mZhczjfIj7g+g0Fm+T0DLwvDBTZ8+Henp6ULlqMYTkR7M\/NXn42dsFtaXAZ5kz\/zF79SpE3766SchOGeITJVSkAClR+EQQFSw4RqLi4vxxz\/+UTgKdunSBQMHDhSR2ygDkR7UpDQklcT\/s43MrczYJ5deeqk3j3AwIBBKG7NzCqXvcLa15EynamHo2u5yucSvyv82bcPI8fdh8fwnBAVCcOjWthFGDeiAW\/J\/gNZkXapv\/3xOOzz5cTGmX9XVaz9iBUCK3EW4v+x+9J69FHCmwoV0pBc4sNy9HPnIRy5yMdg22Ls\/RW43utGOxAWk5gCOrBxkIQt0pivMzIPDUQxnVgpchQ7A5cD8yWVo1GgnXlg2BatefQW5yEQqnEIOQhBxfF0C2x9ThSeeqq5VNSaBamG02hajww1XoiWzi5msAIiRxEsPQBh\/hDFBKEyVD49jkQK57bbbQPbz2LFjoHqyV69eujITth06dCh++eUXUfe8884T2yEB6pprrsHy5cu9glB\/FAgBg4Bw9dVXC7aUhezPzTffjFdeeeW4PmIUSHCwYikznbQDOf\/880XIw\/vuu09EpmKcjSVLlog4plTvqoJRPYc6f6bsZkJUCjTPOXw++u7rg\/b\/9wxe3TIYAzv9gsQH04X9xvcHfsC+prsEZfLJsU+xql49vHCkDKu3jkXJ9EcB2oA4CuHM8VAgoiRV\/dSm8jfXyd8gwJUk2JbsSg9fpyMFTqTCsTEZBckuOMRPs0flq1WrymhZ\/DwQLYyqbTE6wnAmWjIDEEmBcC7+AESVc8hfccZu0YYBoO3QRRddhC+++AI0jKJ6nkCg1pPt+dBZl7KMdevWCYpBhkxQY7HISGv8zOFwiDgsqgxE1n3nnXcEJVNeXo79+\/cL0CILw7CM2j6SkpJ8ZCAnW6tkdk7BPffwtzIFkKnvFIlR447sh\/3bxSg7ayROPfADDn6+wEv+8eJowYEAUrB+D15KPxMyvaU\/U3YJIHpC1JL4YsxrOAOP7Jonfk3OGZ+JcVsfRR4yUfyfH8WvmxS63tLuKlxb9gD+br8YM745jEX93oIjuxg3rcjF5L+NwpoZt+CqZ5fhuabbkeokIuSR+IDLkYuMvGLkOVKA1Hw4XPy8QFA6qYWV\/v\/8uLCw0o4kTwCIWtQ8JvU7nVkjACK1MKEmWjoZF\/O7777Dl19+icaNG2NjJfjec889PsJuuZfPP\/+8kJ3wPEk5tG3b1tKNN9KwfPjhhzhw4IAApJtuukkI5Y2KZHskC2Np4BqsdDLOKZjl+GVh6HHL6OpWHN+sAIg6Qa0pOwFk2JIPcfX129C6004039ccwxoPQ3d7d4ytd4sQXt598D6vJeqEOxtgwYse\/xfVQvWthq8jT7AdwJZb01DgzEHaxlRk5+QhO51yDn6jX\/iVy5EtiPiMPAfyMip5H2c6UguBVBTClZSEvIxUZMNBPIHHHtVTKFBlaXz+dQhEC2OFhQlnoqVIvpgffPAB0tLSAlYFh0NFGwOQYOAD8GvK\/tA7RbjvsmRLwYK0mhWjmCAqFdIvqQnOOLUx\/nKdHTnONCHUTHF63njRwEQ489Lxqu1VZCBDrG5y2WSvJerQoXZvUCBpQUrL1SX2JchHIZykLgpTUeDMRlqBA6k5NGFPB7KKK8GoUv6RBxRkesCEFEhmrgPO1NRK4MiDs7hSmpKdD4fDieSNqchIcyDd6fLIQ1IdyCzYiAElP+KO+bfh9\/WrhNg4rmEztB7\/Ihqc3t9rZcq+pfm6kSm7VZP1cCVaimQAIaXQrl07H\/+n4K519LeK5HNSd9cvC\/OP1UV4dPIE7N283vBEOiT3RLsr78ei8QNwS973+PjHPWjRpD56tG2Evw5O8nGye\/PrbXjw7Q14c9wfsPvg7xi58Du4HxmFK065GCk5zmrzUTmaw4H8lCqKweFAQRqEvIMOdpmuLGSlFsJRkOs1QR\/dqRNW2oYI2sCRCWQ5nQIE0goykF0pHcwm5+HMBRzJyE7OQBZIbRBAHEgrKIDDlSwsVdMKPKCR68xDXEEFHI5kpCZnCaGqlIEQNH4e1Rwo+9UzW3szJM39QYTfU+Ujek5zJ0LbYnRg0XIxox8CQltBtJyTXxaGQtENOw4dtxPd2jbEu3d6TNClmTr\/fuPCb7HnoFt41rZtbEPhvX\/0qfNdyQGMXfQDfi07Kj5v69iNpPSlWNbsNa\/KVDtYTjYEVZJWAAEgqZR1CpEnmYpUZJGUgAdYSmw25GcAzpQMuApd2Jidh8xKvMjLyEZGHsEiF3nFqUBSjqA+SFGw5GVkID+drEomnELWkQFnUp6oE5ddAUd6JXWUlo4CVzWAXPjbTrT7oRBLug\/wMRizIhQN7WqF1jpaLmZoq4z+1tFyTkFbovLB5r76OpZsaoGR3Y9h16+H8HXF6Xgpoxd+2VWGe9\/6H\/Ize\/mwP2RfFqwsEYJVll6u4TjzlNPwnuMpQwDJ83AvoPKE4JGbWX050oRNKPUG1YX1c7JygZxiwb4QULKzKu3b8xzIQBLyitmP5ztZMnNz4SguhiuJgOKEozgDeUkuFOQ4kZxRgYykNDgLU7GRIFSlhZla8iMWNWqB4oRGohvp8xLOsHuB5MYl5WOlRMvFtLKW2lwnWs7JejyQqtMa0rOlAICtxRuxfft2lLXqiesX\/BctdvwbC+8ZjqMJLXD9gm9BnxZt2EIVQIrjinH9t0\/C0aqBXwAR4FH5bvPSPQBCFkPKQWn8lYlcHxDxAEgGXJnpyHWkibqUgZAUyXC4KjU3HvCQlAzZl5ysLC+ApBQ6UVwJEtkuD8uS6ahAdqU9CMWyBdkU0HjUuASQoh0uPN7zQjgqjtUIgNTEAwnHxYw50\/mejBqNXo34r0a\/95cHWe+cw3FONXF\/tH1aigeihiaksRjjgTx8cTN89eUqoWen0c+rr74KpsZkWgGjQgBhJHeyP283eBVLPqpAp37r8XKbmYYUCPsi+0JBJ8GDJcPDtYhCVobmYZIS8YANAaMAqalpQvjJ9o60XDjSHXAWu5DryhR9sE\/aXOZneMgcglRBWjZyMiq1NXkZcGRkIxukQDLhKHQiK89jtSoB5PDaZZgw6EYkn9LQ67YfThamJoIqh+NixpzpfG+4qsGRxnO0eJXGdVRNS0tdvTzItRJAjPxZ5OdzR3THrk1FWLZ8Bf7xzQ6ss\/fBsUYePbukUqT9h9wgVdU7v8kMfPWxAynDNmNWwiPY1GeAEKLWO7RLVE90e3LACpBI9YDIxmQIlmTMox5zdVmHlAhBxEFDsOxsj8wkqwKO\/BxkIA\/ZlJtmZgDp6cLqNBtx3j6SNwKDij4R\/aXneyiWzNRsZDidSHI4BQWSmpKGrBynh2qpApCsw7\/hp+8KsKi5xwdIalusshJ6l8bfZ1QTJ3Q6wycyOYW1B754M6B0jzUFIHXFmc6MqpAAovr\/qKb+dQZAeJnVeCBqWoerz2otoohJWw5GLdOjUlQjMvk4pBr39\/OW4Y6URDQ6dbcwL6fgdeIzb+ObN58SILK4tNQHRDqVu8UD5gN3FANzmjfH\/9u3z+fNOZMcWDm4BDd3WiDUrc70VNi3D0XZ5P64dfqjKHJ\/AmdKaiXNEodFQ2y4sFMnYaJeungx3Ilu3Dq9mxjTiexKG49snLbShlE3bYAjNQ2fdyupDi23cSNIDG3YWYyH7\/DEDlC1LYGCg1l9o\/wvweSFMQMQGQOFdjeq17E6x7rqTEd7E9XJjta0LDQ+k2zMf\/7zH2FBK4FC7nedY2HkhaE9h7BGNYiJuqfXKDw7qo+PsNTMG1f2nSOeIQSAsPgL9rsicYXwSXHCia4lXSs1KKnCuEx+VpRYhAWJC2BfbUf\/KVMENUO7jrRcqoFZ24lHVi7CBYMSsajkFXDsR0sWCUVu+\/79xdhkhorcy+FwZlQCUB7GbXsNN7RfhoKcPAFeJZ9\/7gURynIJHYphvBkOhPQ9KZBDX74OGWNV\/jvQhNNmAJKMZLGnqtexduJ1xZnu4osvxubNm\/HVV1+JHM909tOa8GuBgUAiWZeWLVt6ty5QQzWzcwrpMoWxsSUtjERXCogkaSY\/o5p315mj8M+Jnhww0qAsc8CpPsmm9ObMXzs+bQkg\/tbFB0\/w4K8iLzid51iSKjUrvOzsR6+wTVoVUHm+l\/UEP1L1p6ol4wM40ysFJoVCY8OS5ypGRYYvTBCSqPsxDq8UxhNSulLtS4K1JTG7mJIC4bAybIEVAKmNznQMEVCvXj0h1\/v4448xYMAA4QOWlZUFFRxUSkRSKcyBzHzMsq5ax8rtMDsnK32ciDqWAERLukmQKNq83esfoybe5sSlrYi\/lJd83AQGI1JZ3QDWS0GKIVD42yyOk+nMh6MwHa6s\/EpwyPWYkDldcKRmVtI1qVU+7J5fXpcwVvdAW7pmRMXt7oRRH+G8COG4mP60MLXJmY4UB2OEXHDBBcIDnaEI6NynhiJgGAJ+roYwkI54qrxEzZ9s5TzDcU5Wxgm1juWAQvE71+EUVyEOn52Jbh2aC01K3nNzhDdkKA5IcYgTpup8uEZFsilGv4hWNoF9pFVqX1xOBxwZTrio76XjnKNQONRVD58lmCM9ikaCB2mYdI0vjJU5REKdk3Exo9WZjmzt22+\/jeTkZHz66ad44IEH0KJFixNyjCfjnIJZmKklKmN2MKCQmnFeb6BAYqJq2\/NxX7HjCmx9bSsarG6A3k17Y8KsCfi66dd497d3UXp7Kexr7D5h5lR0V8PPGX0ufxkpLG10RyNcedWVaNigITaXbMa7c9+FrcSGp69+Gv1GjMDC8nLkFRfj1\/\/+V8RMPe2ii9C7d2+BMWd98w3uHu4JeHwiw95Jh7oj61ch4fSBaDt+AbbNvl5EdrfikCf3PJIvZsyZrvplRPI5qe\/XUkAhf2yIV8iqQQUrLIxsQpJ4xowZmPj0RLzT\/B2sXLkSccVxuO3S2\/D+ve9j6tSpoirrzJ49W\/xdhrCz8vnu3bu9bcm7MhANo1Exrgl\/VbT9y6BJIyZPFnlicjIzhY0L407ojavyw8GguFkbCR7xrRLRYtgk7HkzB+0mLcGv78\/xpoSw6pQXyRcz5kxXiwCES6HKlXlutTlu5TJVu468VSUY2quNsD6VkcmYcS6YwktO3pNGanPmzBEBi2iMQz7zxhtvFF0SBKx+Tgm6Wsz619bnWGTVGHRGb1xt\/WDW7K+Nqq6lV68EEOaXCTRaeyQDSLj3LZr7i5ZzOo4CkQJSrROdWVDlgo3l2LDjoAAbq2pcowNWH6wMdMu6BBBGReNDDuRzrYzGrH+1PgXIkuogJaM3bigyIKuXnGrbo7tL0PK6LOx7\/2m0Tp+N0umX4ZRuA064IZnVOcfqBb8DUQsgektWHxGD1KakpAi9+KMzHsOKLfUx4voRGNitBcbkfY8XM3ph1Ya9yP1yqxC0+mN\/9MairELGXJWUwuOPPy6qhgNArPQvAUGqqkn1cL1G8zkRAML1a8Moqnl3rV7VcFzMuuQLo11roAZhfDtSfqjKCc3SWoTjnKzeiVDqBazGVYMqc3OmPJCFnd2vx\/xb\/oh\/fL3DY3QG6DrTmU1UUgbyQXITA2FVjFgbyWJY7Z\/1VdBUrQpPBgtjtm+BfB+Oi1lXfGH0TND1Ukb4238CBSlmGS+WdSlLk1HkjVJIhOOcArkXwda1BCDyl9hfUOVQBYlGuTwCEZZyE4yEnIH2r5cXRAUVjiWFuqGuPdjDC6ZdOC6mHoDUVl+YQADDjKrQMybT2ljJMw3HOQVzPwJtYwlA2Km6UFIhU2guHqacqKrqVS5AGt6ohjtqpGy1jdnngfavZkST85EqW6NxA934k1U\/HBezrvnC6GXI056fv2RZko1hGwr+5Q+OBBw9c4BwnNOJuGOW44HQw3bmNadj3OIfLQdaPhELqM1jnBR3\/rQ0T2hJeh0XVKfwVPe5rvjC6CWWklQErU+tJMvSCuT1EnyrbE6toUBUT9uMgYmYt7wYc27ojkVflh6XrqE2P+KTuTYtgATqOKc3d9NftqpsfDJsgV4fdSWxlNafhXvB\/WPCqksuucRSsiyVBZJUCkGFoCGTkeuxSabndDIvpjK2JUOyHfvLMXnmAiQf24j7HnoY4+a+j\/8tfhBq9KUIWU+tnYZqiSoXGWheXPkA6MthmDhJUiCsrMl\/I8etK74wFHxq1ypZ66KiImHebpYsi2k3pRZGfS+qvERPsxP1AMLLIqOPjTmvKR56KAcvPDUDty76AfFr\/oZFz8wQ9ykaBYm1BWVoG\/L7hi\/RcfoXIhq8lXIyLma0+sL420+ra7JyJkFRisF2HOZ2pkJUWqOOeHYVGvxnIQ73uh63n9cUB39cAdpm0NnoRACIVqh5In1Q5C8QyVkzew9J3rIe58x\/q0KzUM5OTx5Cn5gO0z4QCbqtlpMBIFbnFvOFqd6pSD4n9TxNAUQr1JFkGK0y\/ZLCVm+Nn3pSfbxp0ybvQ5QPmqSkGjYuDMPpdmEVQKzWC2SeWtAIR9SzSL6YMV+YWgwgRhefjzw+Ph6nnHJKIG\/DUl2zyy4Bhp2RIuK\/yW+Sd2WmeP79zDPPFImVmUyZloD8fMKECWJ8rWpWygXInzIYDEPTsTAKlaRAtDwx+6BhkOqtTJ6WRVIglLyr1Ijaf2Jiok8sCZW6UgHEdciGqUXtsaO8vuj7yW5b0bfpEXwX1xmTvor3yjTUvpn8nOvg3B555BHRbtSoUXj55ZcRFxcn9iNcFJKlA41VsrwDZnffckc1XDFoNW5yU4iM6CxHjx4VMSDpscqsceEqZmyAFQDhXPhI1q9fLygmKQRjBC2CCkFCS00ZAYgMHCNN2zk\/2YcWaNS5c2wmlOZYEjAk6HEeElzkHLUCTu06\/bXRAgjtdbRBbyoqKjBr1izhqHiiKLlw3Ym60k\/UA4iZGnfawAQcPXIYf\/jDH8SZMhP6WWedhY4drQnzrFyEcAAIqRGyOloWQ31oVgFEykBUKkSydP4AhN+RCiDAnn322V6KRlIuco7S4Ij1VHmLmndEG9lKe9H0KBAZ5Jf7uXz5clCD8Le\/\/Q3PPfecoNZOBCto5bxjdap3IOoBRPWopRpX2oEcOHJUpLOcltocqwuW4YorrhB5Yd577z0R2q1Zs2Zo3LhxWO5CpAGIfPA7d+4U1IRkEfRYHe3cZQySvn37CopIbcPYJGrxp9aT9bSWunrslzo\/Uoj+AETGmeX\/Wfh\/vXi1\/qxQ9Q6doLh27VpcdNFFIhxgw4YNBcVHWwqCJFWh3BMWWS8slyfKO4l6AOH+SzXupMFJePHzLZh2eReMWPgtBp3WDFmXnooff\/xR5KRVS5s2bdC9e\/ewHJ\/ZJlphYcJJgdApShUc+5OVaAGEayEbw\/3p3LnzcTIbq1SA6t1J1sTfnKwCyIjJI0TA6mIU+wS41otZ688KVS\/niaqZevbZZ6Fad6oRzMm+SVAJy+WJ8k7M7n6kLM9UCyPzuMgJz7iqq2m09XAtzooWRitr4APnL7gUoloBEJX9kHIOmd9D\/U6P\/ZDUiD8Whr4P2ocvWRS9+WtlILKtXIt6uVq1auVlibRztwog7Sa3M4xsn8d4scwGCE\/2Pn8AIllB1lPPgJTHgw8+KAS5FNyOHz8ee\/fuFf1RRkOrztLSUtC2wtDALVyXKkr6qTUAIvdby4er4QCtZtsK9uy0diAqia+S1CTrqfJVtTBWAESS9y+88IKwrk1NTYWzMjOdEWvCelBJvgQAABL4SURBVCyM1v3555+LS8+HISNzcx5t27YV+URULYcqdJV7pu4r+zSycfEXl0Luj3buVgFk3eR1yK2MVG9UMkUGYs\/3eiyM1ETpJV0ihSQpC0mBSBaGgCfTQcQoEN\/dj3oAUWUgTW3HxOPgBWWYwfz8fGHHT2eiL774QpDjDDkYK5G\/A3oXUwUIvRVoAWTZsmVeVbisbyQbYfCpGIAEfi+iHkC4ZBnb9Prejbwp\/fjrIQHEKBhK4NsVa3GidkDvYjIbnb+UGVYAxF9cixiABH66UQ8g2tio9h\/eQtyh3Thy2hC02bkaS56eiuzJf42pAQO\/Gye1hd7F9Je0y0wGoi5G9i0\/U4W8\/DtDVdJuZtKkSUL4LmU2ZDspFxk3btwJTZVxUg\/CZPCoBxC99eldEDP\/kEg+pLo4N6OLyZSW6Uj3Cku5N4FkDqyLe1mTa64VAKLVwMgNCyTnS01ucm3omzoOj37Dt+R48o4jKckT20eWVCUFMG01pL2G1b3wdzFJiUg7EPbHbIFW8hYbje0vWbrV+dbVevQLuvfeezFx4kQRMoAWzPwTacVSZjrmetErBw4cwO+\/\/y7IURpD9ezZEwkJCUGtkZdt6dKlGD58uO5GcQz6tDRp0iRgc3mjtjILvdGEa2JMjkVTLeID\/88\/apJuOebZZ3tSKDI4mFoIJvwsNcMFZ26y+Cob2ejo7ojh+4eHdX+COsiqRjxPPoA1a9aE0k2sbdUOEETofhBpIGIpoJBeaobNmzfjjTfeAP0qmIR4x44dQvcfrFOd\/GV86aWXhFxFlpW2lRhiGyKsIuevm4927dr51fgssS\/BDWU3oMRWIn5NB7kHoaysDNu2bRNtdzfZjbG2saJ7kujzy+ZjsG2w+Hei2xfh1XaBapmM2pbYbBhrszGvNyeAEevKkCcyfXsK2729dh9yUj3ZgqcVuZGoGOvlMY2vA1jZ9RW8ansVh\/t\/CMopZtpnIteViz7N+gSsEVPnGgxA671yeZ689OF0b6iLiEIQnjt3bkTayJgm1zbKTEe13aFDh4QHLI2ZKGmn52uHDh38nrEqR9HLaUsHPdpUsKy2r8bSJkuRWJ6IVfVX4alvnhKWnKRy+N3IDiPFoydYLC5djC22LZjbYi76lfXD0sZLMWvnLAw\/MFyEnqPBF9v27NFTfM92\/Q\/3F23ZhmXi3omYuG+id\/5qu0ApK6O2q+12jOzQAYk3uzG89X5MnLjPZ79ku8+Sk\/F8+\/Yg4Ezcuxf9ysrQv6zMcG9\/PvozRiWOEvMfcXhEQO9MnStBNhwJpKOFhw9oo05S5UjeS8uZ6dS9owwkd0Qi3n4tHwMHDhQ+IQSS2267Da1btzbcZn9pEeQm0dW89bmt0dPeEw6kYqB7gXjsJO0PHjiAx+x22O12MQapk1W2VeLv\/BUe6B6Ihe6FeNT2KOrb6mNy2WTxHee2fft2tG\/fHl81\/kpQJZQ9OCv\/UKSw3b4E5e5BWGWbjg\/L5nvnr7aTY1q9Q2rbbXa7GE+yLGvsdnxS5MagRF83AO1cOaaXYqlq\/VNZf0GZaAvHe\/Dog3i\/9ftwpeViRL\/2oLwk1eHwkaHozV+dK\/2Yptume6uRSiuAb3Bl1bBPDdPHz1lolh\/Jl97qGQZaj\/fbatDkQPqO5L00NWX3t1A+asaVoFkyL7uZEx03wigxk9yk2xffjuf7fwQX0sXzJiFPIePCKjJ+g9vtBRB1blKeIYWSVRyCAIgitxsJpaU4r00bb1vWy6\/qgECiFj7QLAA3lJXh2f3PYkKTCbpj+tsbsgU0zyZFxr3JrAIsR6X8wpkJbMz1FY7KvrTt+DnXloZMpCIdziqRK5+0CiRqu5kz7SCr40rNA7Jy4MjJRVauxzmOhaDAfZX\/1htTzoeCVVWQynMqLCz0evCq6QxodbtixQrvdzUdcCqQR3gy6+pFXQ9kPrUCQNQgsOriAwlKw40wyi3L7wb0H+ChCVxZgJPgATgKgYEdPb\/UBw7sR1YWfHj8lTYbhthswtCagCCpCi0o\/LRunaH8RLaRoMOxZtrtVX6pcZjmnoYH3A94lz3ebkdKFcBoLwIpBs7p3AMHYCsp8Rlz6FC7EIAOHOjGwoXHUx\/sSys7IYtFmc2IshFCxcq58g+pmGluNwZVAauRzIUCVspKSkpsghIhaJBa+7DsQx8AWbtvLT4qW4cGdjuyHNV6IS2AaPOfqOunqfq5554rznjYsGFCJhatvi28j\/Pnz\/e6KpBFl7mQJOutzRdD\/ytJgch8RtK9gPvBPmSWOulTxDAP0gXCKEh5VAKIaspev\/ygN+PboiVv4t9HkpA3+c\/Ie26OyFpv1RbEFEDWlQHFDiQuT0TJShsSq0j84cP34\/XXG2DbNjvuvHNX5Vz2+7xbyhWWNmkiPhu+fz+22Gzo6PYIH\/nQWh88iNM2b\/bKT6ygP4Hg56NH0Wjnv\/FEzycE60NWitqONfZxSHQPEjQA61E+QTmFHG9uixb4pKjIK3eh\/GTp0sZYvbpBpVpur3ddevPQyk7mNp8r5ECfbf7M77T9yWsIQmte7wH2VXLr9EoVUDZ+GT3a2x\/b3tjhRqz9prn4rCK1mmXRAgi\/14uHQt8e1deFIEKjMQkgBE6tRsloQaqqWtZxan8R\/OxGMO0Jrqq6XKW09JJGMQ7ON998g7vuusv7gyZZmJtvvhlPPPEEsrKyxHcECAah0gMQyuVY+IaiMUudJRlI3JH9sH+7GGVnjUR86Tfoav8Nn+TNxJEDv3lN3K2kd7TCwjz22GJcdZUnPoRawiWPUGUZMn+S0V2UY7oT3UIoy0Lyn\/9RjpKOLNiqNDhSLUuFyk1uN1rt3++VuwQiP9GukxooRnkjxeCvhLo\/l9S\/BJfaLxUqYTWqnB6AqPNQLz3BQrrr035BBRDatWQr2iZ\/a6moOP7buDgrsO+pE0x7zo3UrSy8q9IEX2tAyTqkQlgkVSK9j0mBXH755T55XyQLowcgBw8ehHTOZH96VEjUUyDSmS7U3LhWhKhGJK8\/Pp2Xqyp8KfLzgXSKT6oKf40oq4mP3+yVR8jv2E795WFdGm6xZGR42AlVjiHbyeA7BBIZdEcKSCll4K8l5\/HQQxuPG9PsGWjHpIWoXmAfbT\/+9sdsTK5jlHsUFtkWedka2UYLINr8rmpMDxVAmP5UZWGikQKRAOKPbeM+ySTc11xzjYj4plIg\/GGViaPYHwupDUnhkIKXnxmdU1QCiHYx4cqNq6K5Xk5brR2InIc\/mww+2PHj7QIoKF94\/XWPlkYtV1+9Dy+9VOEjP2E7\/pGkNf+\/Zo2nrZAXpAL79x\/AggXlPu0oT1i50oZihxPZjkxscG\/wsiUUYC5f7kZRkRuLFpXgnHNaBWSXoV0ntVE3uW\/ykcHoXTS2u7WSdVtQ7jtXM\/Dg98vdyzHe3kMIVsm8BKKFYXt5jmrAoEi+9Fb2RKVAWF8rA3z++eeFjISyDhaVAiFbo+Z0ZtgH+vrI2C00uqTs8OKLLxahLqUMhP3oyRQjeS9D0sJYOQijOjwQNRSA3CTVDkRta2aTwUctZSb8Owv\/vXq1XdiBrF\/\/O0aPPmbZUnbu3OZYtao+3O6jeO210uPadenisQJl6devzDs2AYhyjmHDduO7\/d+hd5PelsdkX+o6dzbaiQs7Xei1Z\/G332z3cr16+LRzZ7xWWhrQ0XhkIB2QcuyY0D7F7EAC2r4arxwDEM0WSyHcOeec440lotqBUJKvLaHy+NIOJBR5hDonCVISuIYM8Wg55s8vw5Cu3cSvOakCCl9pj0INise21H9R10ljuaH2ofjE\/YmQufgrst0rXbsK7dEDbreF0Tw9Us3du3FjfOJ2C0VvOCLrR\/KlNzuDSPs+kvfSlAJRQ\/GRvHryySdFbE+SaMHkxpX84mmnnSaC6MpgRHKTgmFhzA48WJP0YNtREPmz+2e0K2snbEBo8Lbdvh0L3At8gICfa4FBHZPaF6pwD5cdNluij\/p3ZoMG2G63C2GuVPP660Cqwj8sK8Ngm01EUgu1SGewaFXjhrr+cLaPWgB55uMNyJv7CCaOHIpTz74Eo\/\/fQ0jY8i88+8JCXDGwj+ALA41IJoVw1IdLmxCqukJlYczIe2nKHohJuhnbFMiY1OLc26bSu3LvRGFeT5Uw\/00V62ulr3m7CnZMbTuql6lO7n\/4MIYfOHDcVPk91c73tmmDpY0bY+j27Xhq715ht0LKMFwlBiCh72RUAgjtQMa\/9C+0+d8b+Ov90zDx7S34yzk2vPfWYvyceBUW3XYOaB+ixsE02ypVmk0KRg9AaMoeDSyM2Vr1WC6CBakOOvvJMt4+HuXuciEkJYsTLKumP54N0202lFfaxEhLXo5Lkwo69JFVIctCO5ZUl0uY+tOaOBwUSCQ7gJmdXaR9H7UActNzX6Djz2\/h6lsmYcl3B8FkUvMX5FoGEG3AYwIDjWvUwgDE9DQ026RQ1JTBtg22HdcXSFuqSqlKpcl4\/7L+uqpjs0vtb7y0qsbUsFA+wn\/T1pQaa\/4\/kLmazUN+b3aeVvuJpnoxXxjNaf1jdREm3PsAEvqn451JF+Cb737C\/TOfw8KnHsXlfTsaWs5ZOXTVKlVlYWqDDEQCiAwhYCUUAAGELvlP73sa9m1207AF2j32q+aucho0OpeadOePsTAeFbA0YbfyNrR1IhmM\/QpRVQGq0cID8YVR+zACkNGjR3szlan1aeNBp72mTZsGpBZlH8G2DbZdsGNS6+JMdQpV7r69+\/D02qct37dwzZXn2aNHD8vjGlWM5EtvZXExXxgruwRYTq4tuxvSsyVeSj8TjU4JXxJt9h2LYOXZ4QPDD2DnrJ1oPrc5Wsz1RCU7kYUm6PwTatECiLTUtdKvErXRWz0AV5gqn2PfkczaC8dNpUnMF8bKSfkBEDW59uzrenhlFNpug6VA9KYXi6Fp7dBqsla4Ym9qAYS+QhZdYaDjCoMAXGGCas+5qdK5mC+MtVtmKaQhtS00uWWgmOQz\/yCSa8+7oQf0Qh1aGzZWq7bvQG2gQGK+MOa31C8L83zhZuR+uRUvj0zCc7Mfwx13T8HoxcXIHHDqCcuPa76EWI1I3IHaIANRk33HfGH0b5kld342tbk+g+3gDpSd+WfE0jpE4pONrDlFO4BE0m5G8l6amrLLjQxHRLJwH4o6J1UWo3r86gVu5jysfK4330DHZB8ygTaDysiI80ZzDHVMbbJuNRF5IGOGelaRfOlDXduJbh\/Je2kJQNQ4HjKr\/IneRO14NFKbMWMGZs+eDRlzgW7S9913Hx544AFMnTpVNJF1+Pe7777b8ud6AZICHZN9qKpwaRPhLy5KqOtkHA6mkKS8So59ww03YMiQIbrrtxIIKpizjuRLH8x6TmabSN5LUxnI1HeKgKO\/w77uXRzpeikqEppEJAsj7UquvfZazJkzBwzwq4aT4wXQC+hs9Lmam8bo8piNSapo3rx5+POf\/ywcEPmo2S\/bGQWXNruoZmNq581xGLSGhkzBjmk2J73vI\/nSB7Oek9kmkvfSrxYmPfd7zPpzd+StKkGXY5twwPUNdnS5Ct3aNoo4Iar6UPQCN6vOe7wMDOLCCGtGn1uJ82o2puxDUgIqgBgFlza7qFbHZD8qpbN7927DgNZmYwbzfSRf+mDWczLbRPJeWlLjvlq4Hq\/Ouhd7N68\/bh\/DaQcS7CFRLiHJdtXC1QpQBAsgVsYMN4AEMqZW7mK0L1aAMphzieRLH8x6rLSJ+cIouyQNydJOb4mB3VpgTN73eDGjF1Zt2CtUu+\/e2Tci7EDkL7J8CEbsQThZGKtjSnZCjwIJlJ0IZEw9mVUobJOVx6OtUxcBxGif6qwvjJra4R9f74CQhwB4Y+xZMEq4HcxlC7aNDFar8v1GAkqOEaoQVYJQSkqKT\/5eM6GoFkDM6mv3I9B16oVYCHTMYM9EtjvekIyhqDWZwg0GkQmv1K\/pbGi1BNOeoRTUiHExXxhru21ZC0NLVBlAVu36ZLEwqkpSzoehARjhTA1oqxe4mfWtfG70q6p+bjYm62oBhJ+p8\/fnsRroOtVkR3KeUmVtdUxrV8d\/LS2AMGQBI7VZKRU6xuhxARizB9OeczPKwBfLC2N8apYAxKi53i+jlQsSq1P7d6A2UCCxvDDm9zQkADHKpGU+bKxGbd+BaJeBcP4xXxjzWxoSgKiCuZoySDJfQqxGJO5AbQIQ7m\/MF0b\/llkCEH+BhWIRpyLx+Z78OUU7gJz8HayeQSTvpSUAiaTNjM0lOnYgki99dOxgDECi7Zxi8w3jDkgAYXSzfv36hbHnutdVJOfYMaVAZGR1qitpiq1GWj9ZKty6d4Wib8WxEJXhPTOC8KxZs8CIcZFUTAFEtYDUmkerptWRtKjYXCJjB2IhKsN3DuEKNRm+GXl6Mo3Krlo1at3ZY2rccB9HrL\/YDkTXDgQEINpUljEAia7Djs02tgPh3gHLLAwD0tCcncFppOMa2RsWykZiJbYDsR2oezvw\/wGgu\/kg71ucsAAAAABJRU5ErkJggg==","height":164,"width":272}}
%---
%[output:1417b02e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:51adbc5f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:4e38aa20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:07abaa44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:278538d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:1956b86f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:809982d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:10204c9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:17960a08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:52d36b7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:65eca5af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:3fd1ff8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:1a435c1c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:51a3da26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:1abda206]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1bd90c64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:67262dce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8779b381]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:23cb9a2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0ff9556a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:51721cb1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0b5f4ec4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:56c2e67b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9eb5aeb1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:034c3c7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:50631bc8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:65cb25dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:061681c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:12cf2d52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4c6eea10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:42042d8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9b831cf5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:62a07a15]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6a36301f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9892eba9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9dae0082]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5aae826e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9c710b4a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8fd68342]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:171979cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:61206dc3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:590fdeac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:107586d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:885de158]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7a329501]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9b556a7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:004c3265]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:781eb837]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8e104d42]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl4VcXZ\/ogxJBAhrAmrSRVkK6sVqyiW6k\/Fsv2IrC6ICggCP7JXsdCCCtQWEEVlc6mgQimgUFu1UJaKEW0qJoC0gCwJsgUIEDDE\/3kH5zr35Cxzzj333HNuZp4nDyR31m9m3vtt830Vvvvuu+9IFUUBRQFFAQcUqKAAxAHVVBNFAUUBRgEFIOogKAooCjimgAIQx6RTDRUFFAUUgKgzoCigKOCYAgpAHJNONVQUUBRQABLHZ+Djjz+m\/v37h61w4cKF1KlTp7C\/7dixgwYNGkTHjx8P\/X3IkCE0YcIEOnHiBA0ePJhycnKoa9eu9Mwzz1BKSooU1WDg2717N7366qu0adMmOnToEGuXlZVFt912Gw0cOJAyMzOpQoUKUv2pSv6jgAIQ\/+2JazPSA5ChQ4fSuHHjwi7tW2+9RZMmTQobN1IAKSoqot\/\/\/vf0+uuvU0lJie6aEhMT6aGHHqLHHntMGpRcI47qyBUKKABxhYz+7EQPQDp06EDPP\/88ValShU364sWLNGXKFHr77bddAxAAxrx589iPTAGA4AeAokqwKKAAJFj7ZWu2IoDUqVOHqlevTgUFBbRkyRJq0aIF6ys\/P5+JKEePHmWgsnfvXvb3SDiQzZs307Bhw+js2bMMFB588EG69957qW7duqzvI0eO0LJly+jFF19k3EnlypXZ\/wFuqgSLAgpAgrVftmYrAkirVq3oxz\/+Mb3xxhv09NNPU58+fVhfqHPffffRjTfeyMQa6CoiARAtR4Ox7rnnnjJ6DuhHwPVw0Ql1pk2bRklJSbbWqCrHlgIKQGJL\/6iOrgWQvn37sgsrXlaIM8899xzTRZw+fTokyjjlQMDJoK8vvviCAdZLL71EGRkZuusEN4RxeF0oeGvVqhVVmqjO3aWAAhB36emr3rQAAuUpfqpWrUqLFi2iq666isaMGUMffPAB4fJmZ2ezCx8JB\/LVV18xi87hw4fDgEqPMBcuXKDJkyfTqlWrmHgD0apRo0a+oqGajDkFFIDE8QkRAeSuu+5i3MfEiROZ2PLaa69RjRo12GWHOIHLi4vsJoBwLsaMxM8++ywbUwFIMA+iApBg7pvUrEUAgQ\/H9OnTadasWcy0Cs6jWbNmTNzglhkoMt0EkJ49e9KMGTOoYsWKuvM9f\/48A7S1a9cqAJHaUf9VUgDivz1xbUZaAIET2OrVq5nYAGey2rVr0\/Lly5mF5Mknn6QFCxYwfUgkIgy36uzcudNSB6Ktq3Qgrm29Zx0pAPGM1N4PpAcguNiwukD\/cMUVV7B\/wSVAwSo6lDlVoqI\/gNGKFSvYgo18PLS+IsoK4\/35cGNEBSBuUNGnfYgAwgFB\/NbHtOGDAX1ImzZtXAEQ9Llx40Z6+OGHmY8H\/EB69+7NRKWrr76aUQpzePnll5kviPID8enhkZyWAhBJQgWxmh6AwLkLeof33nuPLalJkybMIgNHMysOxIoG3L9EeaJaUSp+PlcAEj97WWYlegCCStz3A\/+HdQa6EXAiABWIHEY6ECtSiQ5qUJDClR16DfUWxopywf1cAUhw985y5kYA8tFHHzGRAgXWmBEjRrD\/69UXX+NaDSgCCOry17jwON2wYUPITZ6\/xoXeo3Hjxuo1rhVhffy5AhAfb46amqKA3ymgAMTvO6TmpyjgYwooAPHx5qipKQr4nQIKQPy+Q2p+igI+poACEB9vjpqaooDfKaAAxO87pOanKCBJgU2zjtCpAxepaoMkumVcumSryKopAImMfqq1ooAvKPDvt07SuyMPhOZyy9h0T0DENwBy8OBBWrlyJfXq1Yvq169vuilwTDpz5gyLZ+FFHE0vx\/NyLBBZjefe\/Y8lLf8yJp8AIry07FONfjm3gXuLM+jJNwDCnZjefPNNFl7PrBQXF7P3FHC\/Tk5OjjqRvBzPy7FAODWee8cnlrTcvfp8GAcC8ACIRLsoAJGgsJcHw8uxFIBIbL6NKrHeO+hA9m8toqtvSvVEfAFpFIBIHBAvD4aXYykAkdh8G1Xife\/0SKEAROKAeHkwvBxLAYjE5tuoEu97pwDExmEQq3p5MLwcSwGIwwNh0Cze904BiMPz4uXB8HIsBSAOD4QCkBAFTEUYBJ9BtCpYRpAYGVG8e\/ToQY888khY\/g48296yZQuL\/4DUAKiHbGQInYc4EzJFWWEuUymWAHLs2DGCOT2aBaZOjJOWluaJBc3L8bwaC24O+PH6rNjiQJAgaNSoUYQ8H8jw3rZtW5YACJnWGzZsSHPnzqV69eqxPpHKcPjw4dS+fXsWWzM3N5fVu\/XWW2nq1KmUmppqeSYVgMQWQC5dusRimW7bts1yr1SF2FIA9wzR9WvWrOmpO4MtAEGofUTvRvSqjh07htru2LGD5VIdOHAgC0RTWFjIwANcBwLKcLBAXMyRI0eyv3Xp0sWS4gpAYgsg+\/fvZ1wjDib\/YrDcNFXBcwoA4OfMmcOkgtatW\/sTQEpLSxmHAbEEeUKQlJkXHlMzISGBhcLLyclhogpSAiBVAC9FRUUs2hXayuQ8VQDiDwCRceTz\/NaoAUMUEO+JbwHEbL+QPxWcR7Vq1RiAIM8IuBRtWkLoRcB9QCeCoL0iCOn1rwBEAYjCCWsKBB5AIJoAQJBjFZwHUhNu3bqVBc\/VJkbGZ+vWraOlS5cS4mCaFQUgCkCsr4+qEWgA2b17Nw0bNowBBUQcZDYDSGBRelwG0gTAMiOTNFkBiAKQaMKDmEZTHCdoIltgAeTAgQNMp4Eo3cihimjaKG4DCCw\/3bt3p4yMDMPzBPNVQUEBpaenU0pKSjTPHevby\/G8HEtcG8y3MM8H7ULJbD6PLI9HmhMmTCijT5BJAi4zjhd1OIAsXryYKVH17gFeqHvxSh3rlXJlz8vLYyZdmPqg72jatGmIVmYA4kSEQcf3338\/E4+MCtInwswMTsgocbObm+nleF6OBRrx8Q4fPkyTJk2KSwABJ4wcwEZcstFnbp4ht\/riAIIcxkgKpncP4GMDHaUXxRRAoAiFjwc4j8zMTJZ4uUGD8BgDIL6bSlSYEdu1a2fKgYAdPXLkCKvjxXN+L8fzciwcMD4ezLgQT73kQMRMeJiLyAnA\/2jQoEHsLEBZD06TiyHbt29nYjEK6sDd4JVXXmEWQW0\/+N2OKG116cQ5161bl42N3MKzZ8+2DENh1bfM5xxA4ODZvHlz3XvgGw5EdBCbPn16GSUpFowFKTOuzNbL1fHau5CPx\/1AvAIQcKfwNeL6MT0xg18WnrCKX14+Rw4y4J743\/T6Eetpk1\/J7crlWto5641vpz8ndQOjA+EK02bNmjE\/DiOW6OTJk8yRDMF9RK9T5Ujm5HjEzpXdLoA8+\/5e+vpEMTWsnkwTOptb2LSU4BcPaTT79OkT+lhPkc5FZOguxo4dS127dg3pMcz6QV1ReS9m3eMDtmrVSsrFAPWNxtKCmrNdl28VCACBXAzQQPZ0bBjEF22BKNOtWzdKSkqi9evX0+jRo6lDhw4sE\/uePXuUK7v8mQirGQQOZFl2AQ1flhea94TOmbZAxEik4NwDnkNwYBFTa2ovPL\/UWvHB6LLzCZuJTgAsnHcR2DinrQUl\/nc89fCKcwsEgEAxg9ypePtiVAAsXDbVe0w3YMAA5hpdpUoVqaukzLiXyRQEAAF4AER46feTDJrf7wfFutWGay+wtr5WzOD1tX8HUIwfP55mzpxJjRo1CnWjB0R6cxJNuxwAjADECPSMQMyKBk4\/DwSAOF1cJO0UgAQHQLQcCMADICJbzKwi2j5EPQOUlqJYIsOB3HHHHeztlsjViGNouRUACF6Qf\/jhh0wpyxW7OJ+KAym7w1JmXNmDEUk9BSDBARDMFDqQzXsKqcO1abbEF7S1uvhcJOEcwtdff80sfeA0UDjnK6sDMXM10M5FVJTigSjAB\/oX\/B8WH63eRulAIH\/4oCgACRaARHpktBYNDhYiQBhZXfgl5pcfbThnwv8mKlu5SIN6oi+InsUG80LhDmeYw759+9jvRlYYcfxI6WLVXokwBhRSAFK+AASrNVNm8vOg9RIV9SGIUQOuAFwCgAHmXBQjUy0AAK\/LxaJVfmp1ICKAaOfM\/UDgg6LVw1gBgdPPFYAoADE9O0FQojo9\/G63i4YCU8uBaH\/XrgEXGnVkXpy7sX4RQL7dlEUH805QRqMq9LNJl4N7eV2UDkSC4l5eai\/HwtJj5UgmQXbLKtECEP44FBMw04EYvbGxnHgEFTiA3DpjMRXubESVj5ZQ8xWnyKtUltqpKwCR2EwvL7WXYykAKbv5WiuMKA7pOaNF4tkqcfTKVOFzqDT+ZUq87nr2OQCk3xUJnqSyVADiYNe8vNRejhV0AHGwlYFvogcgmRvP0rSmlTxJZakAxMER8vJSezmWAhAHhyHGTfQA5PHKiXTvDdaBy6MxdSXCSFDVy0vt5VgKQCQ232dVRACp27wt\/Sz1Av1f02qevErXI4UCEIkD4uWldjJWJA\/bgqxEldi6uKsiAsiUzjfQT+g4e8jqRVgLBSAOj5OTS+1wKNtvYSJ92KYAxOlOxaYdB5AxL7xGrb5opMy4fBv84Eh2ftcfqPTcIUqoVI9SrhsdOiF+BhDZh21Wa7P7nD8210eNamXGXbC3mA4Xl1Ld5AQampUcdYIpEeZ7El84sILOfj4uRPCU60aFQMTPACLzsE1mbQpAon7XXBnAzIxbOjmdpuSdC40zNDM56iCiAOR7cgM8cNF4qdjgbqrcZhb7NZYAIqPfsHrYJrO2eAeQeIvKLvqBwIz70IlS+mRYDVpTcDF0hrtlJDHzbjSLAhADDgTgARBxE0CMxAhxg0WwWvVFYUSBe3i\/Wg5Eb23xDCDxGJVdCyB3vHuGPhlWnXJ+lKQABOH3zUqkHIHRRcbfvz22ja6s2d51HYiZGGEEII+v2hsWuOfma9Jo7fA2jr5QjNYG7iX362OUWrSf1s8b71lkLUeLcNgoHqOyiwByw4vHCVzIvo6VGRfCC7gPcCHRLOWOA5G9yEaX2qm5zEyMkOVAUM9u+ECzwyPqTxKP7aLUzTM9BZAgRWXnIhAin\/FgQ2JUvmheUrFvIx0I3Nm\/vLsq++FF6UAMdiUSDkT2IlsBiIw4IvahJ0aUnjtYxuqjXVvX+Z\/Tlv8UhroSwwfanYOWnKIFx2sACVpUdr3wh16BhhWA1MotpqyNZ8O4D7RROpAoAIiZPuDk21Op5Og+SqyVSdXueSo0uvZSO+Fi0Nk3L7UhSiggKs2gql1+pWv10Y5lZGVxOgeRpLHiQIIYlV0v4JFfAMRoHgpAogAg6FJPH3Bmw1I6On9QaMRqvX8dAhHtpXbCxfD+k1umUUJqIiWmJ7N\/eeFWHz3uSs\/KYncOnFvBeKKfy\/T3vqINOy\/rQHJenSwtwhiBrcylCmJUdg4gCCIkpseUWa+bdfREGN5\/z9c\/o5pHztCx9Kto1b1tqWeNPHqq5U\/dHL5MX+VOB2JETYAHLjkvV932ANUafjn7mRUHIlo1zPr\/9sSfKblVmm4VKwDR5mA5vaUflRz\/OAyAAAx6jnAAj\/O75oSNy\/1cnHiimoGtzGkNYlR2vwPILX\/dTQ\/P\/keI\/ACQ0sczaWrzH8lsieM6CkC+J532UqTe3poq33B7iLCnLlSm6q0nht4cGFk1jHbi7PbZdOHQfMONSqxxI1W5eVkZsMI4b2z8NNSuwzVp1LB6SpjPCvuwNOOyePR94YCkFXX45xgP1iYAzs5jFWnxjtbSVhgzsJU5iUGMyu53AHl49ka65a9fhcif17IOfflGT+UHoncg7SpRZZyxMA50FFBsoiRdU\/Z5tOidqjcvIzEBdbUcg157XOoKae3oWKVe7IHUd\/sXlOEcZC6oCCKXzh0M41SM2r\/waRNa\/MYuKRFGC7bg1MCxyZYgRmVHUquJEyeSX0SYAXWGUkbadSGS1ywooqb\/zg\/bgkvjbqQ2MzrIboujeq5wIHqJpZBUCjlzkWNDpjh9CzNnY75pikWAx7Pv7wtN4a2fb6Sb6p4t897F6JtanDv\/VhetH+zL\/9whVk30ZMXvAIQrKtXX\/cyMJiWVW1HFKpl0KX+1DOlcqbN9D9GQ5ytIAQgGhA7kfO4GSml2W5jCWXYyQYvKPnLkSF8ByG+\/6U4tLpjHQa12X3NqsOhO2S1xVM8VABGTcCOBT25uru3UlriUb78+h25sfyPVr3\/50qEcLK5Bc\/d1Z\/8fmbma6icfZ\/8\/U3SGUpKT6e3Pj4XqcvYebVBQ983sy6hcWlTC\/u197bZQ\/aO1htB1P53IlKoXDqwMcR9GlORsv1af4IjyPmtkF0DcmH6QorIrANHf8YgBhCfXRuIdxIdMTb3M+ttJrm317b+t8Do6WFyTemVscePcCuBUk66s0Z7Sz77nar9B7CwWABIpnaIVVBnz4pYWq6jska7BbnvOqccNB4IFQVRZsGABderUKUSPoqIiGjNmDFWvXp0l6kYSbqOiNUnaJaqqHzkFFIBcpqGYxQ6\/86jsVs8rIt8BuR7iDkCWL1\/O0g6KOUtBCuhFwJFkZ2db5sxQACJ3eKJZSwHIDwAi5sb1Ouq61R7HHYAAsbdu3UoLFy6kWrVqha0fn61bt46WLl1KWVlZigOxOh0x\/DyIABJDcsVs6LgEEJ6IB+KKWIw8DrXUVxxIzM5jaGAFILHfA5kZKADRoZKVElWGsKpOZBRQABIZ\/bxqXa4ARFaEAfHjlQvhFqT6yceoXvJxwr9+LApA\/LgrZedkB0CO\/O52umNk66guLGIzrhtKVKwQhNnwWj96dMDlYEJws0ZBgB9tKSkpoZKSbyk5OSXsI7ykxQ8veFWLH14OnCwmeGZqL3FCpfq6PiBGfzfaEZiaARhwVGtQLZn9Hz4sHa5NY85urc9Md90U7dbpUADiFiWj248egODdCx7Ria7s+Bs8UX0f0tANMy4HkP79+0t5Qhq5ssu4WMMztW3Jq9Q+bVcIqBCBfdc\/n6H\/frWB+ZsAYG5r\/wsCgIiBlq2OxsqCm2nCzsEkxuwQ23z094UMRGQK3OYBnuKDOdl2JUf30\/m85T8Aabp1dG4FIDLUjX0dDiBwZW9RXI92tqrDXt6i4DVuk5z80N8Gnyilx3qF6yXdXkHEHAh3JMPbjalTpzpyJHMLQNBPJC7WMs\/mzTZg\/M7B9KeCm2l+v6YMRPSKntdryZFi9ryfl+La97OHe5fyniz7aM7iBPD3OiIdrqi30xKIFIC4fbWi05\/Rc\/4rvs2hS1e2ChsUCbcn3FolOhP5vteIAQT9rF+\/nkaPHk0dOnSg3r170549e2y7sjt9C+M0xKAdqp75cjZzdU+4+MNrV217uMWP\/awzE1cmdDY2WfN2\/DXvir1X0+6\/F9DNTfYxECmp3JpatB3CHtNt\/ecbuhzLv676FeOgtG9v0Lfee51\/Hk615HziHUDiOSo7wKPCpQIqSe4cdiwDE1BI7zHdgAEDCA\/qqlSRQ0C3AET25a0dAIHIdDx3KaUcfDbUjL\/PaZ+2k4k7YiIqO31rE0P1\/HEaPXFbGgMQRGWHaIUxuGi1rbAJe78DDkePm+GpKETR68+n+zAdzP9mbDFU4sYzgMR7VPaKZ2ZSja+a0OG23cKOXiBEGDuXxayuGwASaZpHo\/lxnUvNcyvpu8LtTMH703fah14B\/+sJ51GfzAAE3BUXqxpWT2bj6XE42tgketHKnv\/mUdq8p5Bmt32f6pTmsKWK+pV4BpB4jsqekbON6uRso0brerCAyt80q8j2tnbuBSa+tOxTza0rqtuPKyKMGzN0A0Bk0zzana9VnNJIIqVbAYjduaK+WdxXsb8Tf\/5pKAhRLABERWW3v7v8nmT+ch7VudSCEI1drzS8qTINXHWN\/QFstogrAJFJ82iTPqy6FkDcBCrtnGf3yKSfNaCIM65bRUyDxerUXx6jSjfVZGv0GkBUVHYnJ\/GyuwOslXcV\/5YBiFEB5\/HLuQ2cDWKjVVwBCNZtlebRBm1CVa04EDOri9V4WgAZ1TGD7muZHDGAWI3LwxIiyHNiRjIDkOHvpEiZ0a36tvpcRWW3opDx5wpADGjjhgjjfFvMW8pGSncyfjREGJl5aH1m\/ttuKD3y8vvSABJJTho\/R2U3op3fYqJacSDgPqKt\/wCt4o4Dkbk8di01dmOwysyB17ESYWQuqkwdvTm9MG4o1cr\/mPJrXU8tuw9kVrM333yTrGJfRJqTxs9R2YMIINB3XH1TKu3fWsSmj\/\/fMi7dzjF0XLfcAYgTS000AUQrdo3qWIfy8\/OZCFPh6Lu6yafE3XZ6mbV06Fj1iHReGLs5abSn089R2UV\/kY4dO7KpDxkyhPwWVFnkQLzSd+ihTLkDECcK0GgDiLgx4lhaT1TuJCbWt3uZOfeFdJkwC\/Pyo9L9dGLNNEcciExeHHHOfo7KDnBDfJtnnnmGcnJymMISXJmfAcQrcUUBCBHLdg8Q4UVGARorANFyIDx3jBkHYnaZtWsX+7k36xSt\/f0YKQBBOysrjxVP7Meo7Jz76NevHxPjxN\/9DCC3jE33TGTR7mu540C0IoOM23msAASOZNp8Mnq5aWQvs5b7gnNanasSqVV6Iv0i45S0DsQKHGQ\/91tUdu6xioDKQQIQJcJ8\/5w\/0te4sgfXbr1YAohdEcVsbXrcF1znoXPZv3+\/5wBidx+09d2Oym7GgVgpliNdi2x7PTOuCCCbZh2hUwcuUtUGSZ5wJeWSA5HdLF4vlgAi61Uquyatn4yT3LiyY0W7ntsAgvka6UCCACD\/fuskvTvyQIjsXog2CkAkTnksAcQNfYPZEhWAhFNHtMKAIz5z5gxxnYjEUYl6FTMOBOABEOHFC9FGAYjElscaQCSm6LhKkAHE8aIlG4LDGT9+PM2cOZMaNWok2Sq61fQAhFthtByIF9YZBSAS+60ARIJIcVIF1qGXXnoptBoZxzovly6jA4FDmVfOZApAJHZfAYgEkVQVTyhg9BbGC32H3gIVgEhsuwIQCSKpKp5QwAhAvNB3KABxuMUKQBwSTjVznQJGAOKFvkMBiMPtVADikHCqmesU0AOQWIkvWJwSYSS2WAGIBJFUFU8oYAYgC\/YW0+HiUqqbnEBDs6xTebgxYQUgElRUACJBpIBU4e7zXbt2ZQ\/mUlLCk5O5uQzuGo9HeXrWHKNI8eIc8AZn0aJFxPNOG1lhSien05S8c6GmQzOTPQERBSASJ0YBiASRAlCFX2hkCti0aZP0w0GnS+OBk2rVqkUNGzYsA1gcQNC\/LJgZ+YGsaFmR1hRcDE3Vi5QOSoSRPBkKQCQJ5fNquHxjx46l2bNnE\/w94J6Oh3PRKGIEs8zMTJo3bx4tWbIkzCEtEgCZ1H0+XbnnR1S95Xf0P09m0l8LE8I4EKS0BIhEuygORILCCkAkiOSgildR2TE17WVds2YNIa+zKB7ILEGcc926dWny5Mk0Y8YMBkriexkxRGfjxo1p8ODB1LdvX+rTp09omEgABCJR69atQ8Gn8HIbOpBPC0vo+rRET8QXSw7k7Nmz9NprrzFW79ChQ1SjRg3q0aMHPfLIIwS2jBe9xFIIj3ffffdR5cqVZfYlFG1axvPPywuNyXs5npdjiWvz+jWul1HZsU5tIGejwM5mh1U7Z97H4cOHy4hDqAsQ4QCl\/V0P1GT0MSIwaQFE6qK5XMmQAzl69CiNGjWKER6Pitq2bUtffPEFS1kJeW7u3LlUr149Np3NmzfT8OHDqX379gxlc3Nz4yq1pZeX2suxYgUgsYjKrg2jaPfb32jOnCMRv\/i4rkXkOPSChssoUZ9++ukwriUwALJ27VrGnj3\/\/PPEY0PiwO3YsYOxYwMHDqQRI0ZQYWEhAw9wJ1hsamoqA5WNGzfSyJEj2d+6dOliiXtBi8puuSCHFYICIJHEnfA6KrtRaku9yw9OAToLUdTAVnL9iVaPoXdu9dbH5yAqU+2CGJ8Hj5vjWw6ktLSUcRhbtmxhD4u4CQkLgFgzceJESkhICMWNhKiyYMEC6tSpU+jaFBUV0ZgxY1jbadOmUVKSuUJHAchl0gUBQCKNO+F1VHZ+towwHUGTuTLVCECMQE8bk8SKq4DehINQ3AKI2Zfn6dOnGedRrVo1BiCrV69mXIoWmaEXAfeRnZ0tpahSABIcAIk07oSXUdm7devGvvBQ9EyleroK6O0+\/PBDFlSZg4ssB2J2jrVikNsAEglX6JBhtu+JCtEEADJu3DimJMUGIIr1woULwxSrmBA+W7duHS1dupSysrJM56gAJDgAEmncCS+jsoOqgwYNoscee6yMWCKKBFzXICpKIZZDXAd3gv\/r9aMVg\/SUpfzgc8D4+uuv2ZcqlKZm4KZ3YYx0ILtXn\/c8GhnmZ8uMu3v3bho2bBgDCog4tWvXZiAhapvFRRuxfWaEgeK2e\/fulJGRYQg4YPMLCgooPT09qp6EfAJejuflWFxkAi0PHjzIrGsyVjC0w7ddJHEnvIrKbnUGtVwAzjUKF2nQft++fex3IysM6oMD54Bj5l8iAgCP9G7EHZndk8WLFzMzLr8HG39VSLkrz4SaePU6VxpADhw4wHQaUAa9+OKLBNs25zLcBBD0ef\/99zPuxqhcuHCBYCUCkFWsWNEp9yXdzsvxvBwLBODjwRQ5adIkaQCRJp5JxWhHZTdSnmqnJHIR4LBFJaoIIGin5wfyyiuvsKhln332ma7DmDieqEydMmUK0w\/CYGFWRFDnAPTcc89RkyZNQvfg8AcJ9PFvfsjz49Xr3ApZWVnfiZMXFUr873l5ecyke+nSJabvaNq0aaiJGQfiRISZNWsWtWvXzpQDwbfGkSNHWB040ES7eDmel2OBbnw8+IGAu5TlQKJNc5n+oxFUGWdW5EC0v2vnhQuNOnYd0mTWZ8aBwD+refPmYfcge96piLhCJ3MyBRAoQuHjAc4DqAzUa9CgQdg48OZTSlQnpNdvEwQrjHurjazbjRU1AAAK9klEQVSnaAEI56gxOzMdiCyHE9kqw1sHxg8E0xYdxKZPn15GSYo6WJAy47p3RBSAyNMyWgAiWmFERy49c7DW0Ut+9s5qBgZAuMK0WbNmTE6D2VavnDx5kjmSIRn01KlTlSOZs3MRaqUAJEICxnnzQAAIFGsAjWXLlhHiJkB80RaIMrCxw0Fs\/fr1NHr0aOrQoQP17t2b9uzZo1zZHR5kBSAOCVdOmgUCQGDheOihh9jbF6MiBmTRe0w3YMAAlioRsRdkivIDuUwlBSAyp6X81gkEgMRiexSAKACJxbkL2pgKQAx2TAGIApCgXeZYzFcBiAIQ03OnRJhYXMvgjKkARAGIApDg3FffzVQBiAIQBSC+u5bBmZACEAUgCkBieF+jndaBJ+fWexKCZWuTd2tJIcYL0SOTAhAFIApAYgQg0U7rwD1jr7nmGvrPf\/5TJkYOBxCjx6cyZFEAogBEAYjMTYlCHR4UKFppHXigJDz9R\/oIvRgkZo9PZZasAEQBiAKQ7ykQpLQOPG5Io0aNQtHKRGdK7cM6AAVCJGijoCkAkYFJB3Xs+IEgwAuiwyNuiJ6bvYPhTZt4OZ6XY2HRfDzElvAyHkjQ0jqIsU71Qh5oz6\/ReXYTQBDOwst7oHdJpAMKuX0ptf1xgiPuCNJDmBXkqEFIRZm6bszby\/G8HAu04eP16tWLVq5cKR0PZPGXL1H+uXyqU6kOPdh8iC0yBzGtg1X8Ui3HYfTU30qJaqR85QQW70n9+vV17wH+jh8vim8ABCH1AArbtm3zYt1qDAMKyAQUWrdvLc3I\/nWohwebPWILRIKY1kFMValNh2kEiFouy00lqtkBxhcrfrwovgEQLBYggh9VvKcAgHvOnDlSHMj07F\/T+n0\/hOG7M7Mr\/eonPwCK1eyDmNbBDEDsrMctEQaR+3hiNy29yyUHYnXo1OfRpYAdHZSWA5n8k19Tl8yu0hMMYloHIwAxi0qmJ\/a4BSAynKL0hkRQ0VccSATrUE0jpIAdAMFQ0IF8fnQ7tanVzpb4grZBTOvAI6jD0UsUYazoZiftg8wWWo0n04ebdRSAuEnNAPfl9cEMWloHpGlFDhctgFhxFHrmXbccyZA+ItZFAUisd8An43sNIFh2kNI66AGIkfJUu6UiWP7pT39i6WLNilmc1Vjsk9lcFYD45ALHehp+O5gy9IhWUGWMLebKFX+XmVc06\/htnxSARHO3A9S33w6mDOmiBSBcxMAceFoHP4gLmI\/f9kkBiMxJLQd1\/HYwZUgeLQAxSusgM6do1\/HbPikAifaOB6R\/vx3MgJDN82n6bZ9iBiDffPMNQTGFdBBa9hBpJd5\/\/32Wg3fXrl2UlZVFiPLer1+\/Msm0ken8hRdeoFWrVrHN\/PnPf04jRowg5LMRi17keESNR1IsfOO4Uc6ePUtIOQgbPVzEkWy5R48eLGk18vjyYmcuyLuDRMroF34FoNXjjz9OLVu2pAoVKjjqU2+tfjuYbuxHPPbht32KCYAgUTeUVCCG1iGmpKSEJSjGz+233059+\/ZlrxoXLFjA3siIyat48qvTp0+zvK54WPfuu++y15Lz589neWp4EbPsoc\/c3FzbuWvMDiRSYcB9GGx1\/\/79qW3btiwtBh47NWzYkJD1nXsOys6lqKiInnrqKfrHP\/7B1gcT4jvvvEOffvqp6+vz28GMx8vvxpr8tk+eAgi+ocEpwGX6+PHjjJ5aANm+fTs98MADLKcM4ikkJiayeuBEoNDCtzm4BnApTz75JKE+OJXGjRuzegAg5PDdsmULM5fhxSLPngeOACay1NRUVheZ2MEF4W9dunSJaH+RYX3y5MksT3DHjh1Dfe3YsYPNe+DAgYwzKiwsZJn8ZOai1ydABa9mT506xeiIjIFurM9vBzOizYjjxn7bJ08BhL9EvPbaa1lWO1x0LYD87ne\/YyCzZMkSQuwFXgAMv\/nNbyg\/P5+1wzc+gOauu+5ij\/BEdh6XdtCgQSx8XKdOnVzL32t0LktLSxmHwUGrevXqoaoATTggJSQksNgQOTk5UrmEL126xNpBbMF6Oeih448++ohxOxBr2rRp48r6\/HYw4xgDIlqa3\/bJUwCBWILLBb0A2Huw+loAMfPsw2d\/+9vfGLgg\/QFAQi\/qE9fOQ1TBt\/7y5csZZ6AFJegiwH1kZ2fTokWL2NzcLhCvMAdwCgCQ1atXS80FAILsgDfddFOY6zTmB9ENawf31KdPH1fW57eD6fY+xEt\/ftsnTwFE3EQjQoADQa5dXHbk39VyINBv4DOINuBAkIsXF1QsnAO5++672eUD8GzdupUWLlwYpsxEG3y2bt06Wrp0KVPWul0gJmF+4JIgesnOBRwXQAJtAYRi4QDZs2dPplCV7dNsfX47mG7vQ7z057d98h2A4MI9\/PDDTFy55557QqIJ14HgIABA6tSpw+JOQpcCzga6DhSuA4HSlQdnMeNqjGJTuHHguJIXFhiIOLVr12aX3egthDgXjG\/EYdl5XyG7Pr8dTDfor9dHPEVl94Nzm+8ABDI\/lKNr1qxhehIoN6H3+OMf\/8i4BzHaNQ49QAKiB9h9gAq4CR5XBDErOQcic2lFnUukBxiWpjFjxhAuu6jkVQASKWWdt4+3qOwKQHR0IDgeAJFly5YxvQSUpdz3YefOnYzbEHUZMMf+4Q9\/oA0bNlDVqlWZ9QbAAd0I2lkBSDREmLy8PKbkhB4DupemTZuGTr2VjoeLU1yEMdPxiCKMEUDKrq88cCBYI7jWeInKrgDEAECMvmOgH4EPhZXCc+\/evUw\/8uijj7qmZJT53oNSFvMD5wGfFFhPRD0O+pBV6JYHJaqKym5fae83oPedCAMCgWuAI1iLFi1C95b7OsApa9q0afTtt98yJyt4kULk4f4iaMD9J9w0c8oAiOggNn369DIKW\/SB9UGZCh0NTMy8wL8DwANxDOuLdzOuisruzOqnAOT7G2NECG5hgOv60KFDmRIV3+wrVqygKVOmsIsHRy38DXEhP\/jgA3r55ZdD6R24i3xaWhpjVeE\/wcEHOhLRk9VNRzKuMIULPQAAZlu9YmcuAMInnngizOvUzJEskvXZPZgL9hbT4eJSqpucQEOzkmXwNVRHRWU3jgciG5W93Ic0NDqwkP3BfeB9C6wQN9xwA33yySdM7wFAEb1ToVCFZ2qVKlVYjhgUuI7D90JUXOLvMA3j3Q3c22H63bNnj2uu7PCKBWhAbwP9i16uGogyUAonJSVJzwVgA\/MvnM+wPjjgGbmyR7o+OwCypuAiTck7FwKEoZnJtkDEz1HZjZDQb1HZFYB8\/LGuIxk2EBcSFwUgAAtM8+bNmWm3c+fOVLFixbA9BoiAHYYSNSUlhe68806m+4CoIxa9B2zgcqB0BQBFUqDohRUIznFGRcxiZmcu2sd07dq1Y6LO9ddfb\/mYzs767AAIwAMgwku3jCSa1rSSNAntRDFHp7y+NlIXOJnx48fTzJkzw7yWubUF\/jNwtONrM5qg1bc+2vktKnu5BxDp06YqekIBOwCi5UAAHgAR2eLnqOxiBjr+pgkAYxRUWUVlx9ehKuWeAnYABMSCDuTTwhK6Pi3RlviCtn6Oyg5wg9cyf7fEn1uoqOz6VyRmVphyf2N9RgC7ABLp9P0YlR1KeuiyEHcGPhacG8HvRgCiorIrDiTSuxAX7b0GEFG3wQko6iL4fLT6CVEfgpgrULQjXAJ8gxA3BoXrSszEC3HTeJ+w8EHvBjcCGQBRUdmJ\/h\/buIohNReBVQAAAABJRU5ErkJggg==","height":164,"width":272}}
%---
%[output:48ab55ad]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQuQFcW5\/lmWfbiwwsLKa8HFBHncWCugLol7kUusMgXhFYMgGBAwEAoFQoQFFQPURUB8RCARDIZHLB4K1yARyusjcAPELSRCQFDYezGwBKjlsfKQBZfl1tfUHOfMzkz3zOk5Z\/bsP1Vb4jk93X9\/\/fd3\/v7777\/rXb9+\/TrxwwgwAoyADwTqMYH4QI1fYQQYAYEAEwgrAiPACPhGgAnEN3T8IiPACDCBsA4wAoyAbwSYQHxDxy8yAowAE0gS68DHH39MQ4cOjerhsmXLqFevXlGf7d+\/n0aOHElnzpyJfD527FgqLi6ms2fP0ujRo2nv3r3Ut29fmjdvHmVmZiqhhg2+Q4cO0cqVK+mvf\/0rHT9+XLzXrl076tmzJz3yyCOUn59P9erVU6qPC4UPASaQ8I2JNonsCOQXv\/gFTZkyJWrSrlu3jqZPnx7VbqwEcvHiRXr55Zfpj3\/8I1VVVdn2KTU1lR577DF64oknlElJGzhckRYEmEC0wBjOSuwIpKioiBYvXkzZ2dlC6KtXr9Kzzz5Lb775pjYCAWEsWrRI\/Kk8IBD8gVD4qV0IMIHUrvHyJK2ZQFq2bEk5OTl08uRJWr58OX3ve98TdZ04cUIsUcrLywWpHDlyRHweiwWyfft2GjduHF26dEmQwqhRo+hnP\/sZtWrVStR96tQpWrNmDb366qvCOsnKyhL\/BrnxU7sQYAKpXePlSVozgRQUFNAdd9xBb7zxBs2dO5cGDx4s6kKZ4cOHU\/fu3cWyBr6KWAjEatGgrYceeqiGnwP+EVg9xtIJZWbPnk1paWme+siFE4sAE0hi8Q+0dSuBDBkyRExY82TFcuall14Svojz589HljJ+LRBYMqhr3759grCWLl1KLVq0sO0nrCG0Y5SFgzc3NzdQTLhyvQgwgejFM1S1WQkEzlP83XzzzfT6669To0aNaPLkyfTBBx8QJu+uXbvEhI\/FAjl8+LDY0fnXv\/4VRVR2wFy5coWeeuopevvtt8XyBkur9u3bhwpDFsYdASaQJNYQM4H06dNHWB\/Tpk0Ty5ZVq1ZR06ZNxWTHcgKTFxNZJ4EYVowbxPPnzxdtMoHUTkVkAqmd46YktZlAEMMxZ84cWrBggdhaheXRuXNnsdwwdmbgyNRJIAMHDqTnnnuO0tPTbeW9fPmyILRNmzYxgSiNaPgKMYGEb0y0SWQlEASBbdy4USwbEEx2yy230Nq1a8UOyYwZM2jJkiXCHxLLEsbY1fn888+lPhBrWfaBaBv6uFXEBBI3qOPfkB2BYGJj1wX+h\/r164v\/wkqAg9UcUObXiYr6QEbr168XHXaK8bDGivAuTPz1Q0eLTCA6UAxpHWYCMQjB\/KsPsRGDAX9Ily5dtBAI6ty2bRv9\/Oc\/FzEeiAMZNGiQWCrdeuutAinI8Nprr4lYEI4DCanyKIrFBKIIVG0sZkcgCO6C3+Hdd98VXerYsaPYkUGgmcwCkWFgxJdwJKoMqeT5ngkkecayRk\/sCASFjNgP\/Bu7M\/CNwBIBqWDJ4eQDkUFlDlCDgxSh7PBr8FkYGXK193smkNo7dlLJnQjko48+EksKPNiNefzxx8W\/7cqbT+PKGjQTCMoap3ERcbp169ZImLxxGhd+j9tvv51P48qADfH3TCAhHhwWjREIOwJMIGEfIZaPEQgxAkwgIR4cFo0RCDsCTCBhHyGWjxEIMQJMICEeHBaNEQg7AkwgcRqhy1\/8hqq\/Pk4pN7WmzA6T4tQqN8MIBIsAE0iw+IrarxxbT5c+nRJpKbPDRCaROODOTQSPQGgIpKysjDZs2EAPPvgg5eXlBd9zxRYQBHXhwgWRO8Nvzk6QB0jEeNLb\/JSyuixQlEBeTIeM8lZiK8Eyxoaf8XbYcAwNgRhBTKtXrxbp9cLyVFZWirMbCPXOyMjwJZbVAgF5gER0PTpk1CWLUz0sox6Ew4YjE4hkXHUNGHwg35wuoQbNCrUvX3TJqEfF7WthGfWgGzYcmUDiRCB61IcnJ+MYu0WsE0MmECYQnfrkWFfYfjntBGUZvauCK4Hg6DdyRcAvgWsJkUNzwIABNGbMmKjs2Tg0tWPHDnH6Eol5UQ53gSBxDU55qjzJ7ANR6X8sZVjxY0Hv23cZR+84OhII0vNPnDiRkGUb96t27dpVpN\/HPadt27alhQsXUuvWrUWLuEho\/PjxVFhYKDJbHThwQJTr0aMHzZo1ixo2bCiVjAlEChH\/uvuHSOlNJhAlmKIKORIIEt0idyZyR9x3332Rl3ARM24yw8XIOAZeUVEhyANWB45zG2SBrFQTJkwQn\/Xu3VsqGROIFKKkJxBs5eMvUQ+2SE+fPk2NGzf2veMWtOyQEakokUHO766gThltCaS6ulpYGFiWIEs3rkQ0HiOjVUpKikhEg1vbsVRBQl7zre+4XBm5JvCuyo1jTCD+hzUZfjlBHLizpqSkxD8QdeRNWPr4Yc7Pz094jz07UXF7GSyPJk2aCAJBlm9YKdZLgeAXQSfhE0HKPDMJ2fWaCcS\/LiQDgRjjj2snjKWxf0SS900Q7CuvvCL8kmGIl\/JMIFiagEDwawHLAxcD7dy5U6Sus15LiO82b95MK1asIGShcnuYQPwrfbwIZP57R+jo2Upqm5NBxQ+4j6e1NzIZwzr+\/kclmDfDhpMnAjl06JC4dR1EgSUO7hUBSaBTdlYGkvRiZ0blykIDGDhu+\/fv73ifajDD4lwrFB93uDZv3pwyMzPj3bxSe\/GQ8a1Py2niW6UReX71wzz61Q\/bKMmHQjIZP\/nkE\/GDFJZfVuWOxbmgMU\/+8Ic\/OFogOHLh99iF1+4oE8ixY8eETwM5MnGDGXJZ4tFNIKhzxIgRQpnC8OCeE+xIgTSdblhLtJzxkHHmB6dp08GLka727dSQZt7fTLnrMhnhS4N+MYG4Q2oQCC4AKygosC0MJzBcDPF4lAjk4MGDYkv32rVrwt\/RqVOniGxuBOJnCYM1cLdu3UJjgSC7+KlTp4Q8YfB62ylFPGS0WiCvDPouDeqSq6yjMhkxMRA7FBSBGHfw2glsTQaNsrgcHD6+WKxOtzYhh9FXL+0ZBIL4rLvuussW\/9BYIHCEIsYDvwzw+IL12rSJNltxNSI7UZXnUSAFZf4FXY3CB7K9tIKKvtu41vlAnH7oEOeEC8ZxncXgwYMjVrUuAnFa3pvHxA+BBEW0XnXF1QIxB4jhYmarkxSNASDexvUKu97y8SKQWKSWyRi0c1BmKZsnupcJ7YaJW5tJTyCGwxQ3uCOOw2lNde7cORFIhuPu5qhTDiSLZbp5e1c2Ob3VFkxpmYyJJhCzxWEQCHLA4JceT6tWrWpsBliXKH379o1a9nglEFl7xg82IsNDbYHA4QXSwN2lAMUuYAVLmX79+lFaWhpt2bKFJk2aREVFReIe1NLSUg5lD2Ye2tYqm5xxFMWxKZmMXgnE65ay2xJm6tSp9Pzzz1P79u0jSxgEUJoJwbqjaK0PPh5cGQqiKS4ujtSjuoSRtWcA6xWnoMfedgmDXQfcXIazL06PGVy7w3TDhg0TTrHs7GylPoQNGENomeIrdS7gQskgo5fxX7PrJI1fczCCavED+VKfjJtD02pd2JGNmSBwRMNKFnZD7Nam2YKQtWcQUq2xQALWd9vqvShQPOVLhskZT7yc2pLh6GX8QR4gEeN5+O4W9NuHv90ZdJrMdtaAcXUn3jFimZx8ILBCsGmAcu+\/\/z5Nnz7ddmljtO91CWPd9TG3Z0Rye8EpHuOutI0bD0HCBgxbIHpHXSeBWC0QkAdIxO2RxSuBDMzbqqjL\/MuP\/7dOaGMHB\/4TPCqWjBO5qbTHFojLCDOB+J+wssnpv2Z9b8pk9Dr+XreU3QjE2raKBWI922VYMlj+G5HXbIHo0x9pTV4VSFqhpgIyxdfUTEzVJIOMQY+\/zAIxH7nw4pMwD5wdEak6Ua3l7JyybIGwBRITUTi9zAQih1UWSIZNAfPuCXZFxo4dG\/UZcuTAusC1I3Ci4jH7LaxteLFA3NozdoeYQJhA5JruowQTiBw0P6Hs5rgMnDsxHxg1LASQivFYy3ghEPhR3Noz2gjaUpMjGV2CnagSxJJhcnpViiDKy3AM28QIAgMddYYNJyYQJhAdei2tgwlECpFSASYQB5jCBowhpkzxlUY94ELJIGNYxz\/gofNcfdhwYguELRDPSuznBRnJhW1i+OljPN4JG05MIEwg8dB7kZHM7Y7hsE2MuIDio5Gw4cQEwgTiQ429v8IE4h0zuzeYQNgHokeTTLXIJqf2Bn1UKJMxbBPDRxfj8krYcGILhC2QuCh+ognETxxIvFIaehkAJhC2QLzoi1JZ2eRUqiTgQjIZg54YskjURKY09AJ90Dh5kQVl2QJhC8SrzvgqH1YCQWfsQtDjmRPVC6BMIGyBeNEXpbKyyalUScCFZDIGPTFkOVETmdIQ0FtTA8yYMYP27NkTleEM5YLGyasasAXCFohXnfFVXjeBnHtzFlWVf0mpufnU5KFfS2UKc0pDgxSM6yXM52zMB\/qYQFyGOWzMaogqU3yp5sahQDLI6GX8L2xdQeW\/HRlBtsmgmVISCWtKQ4Ms0Bmz09awSMynhJlAmEACoZO6RiAgD5CI8TTq+Sjljl\/uiq2TBZLolIZG+0OGDIncS4OOGJ\/jAm3Oiaowbbz8AilUp61IMkxObWDEUJEMRy\/jb7VAQB4gEbdHllAoUSkN7S62Qj84oZBHZfOiQB6rjqm4TPFjqlzTy8kgo9fxhw\/k8oGtlNm5p3T5ApjDmtKQLRBNk8CrAmlqVlpNMkxOaSfjUECGY9DjL7NAEpXSkH0gmpQvaAXyK6ZM8f3Wq\/O9ZJAx6PGXBZIlMqWh2y6M9eLvoHHyqpe8jStBLBkmp1elCKK8DMegJ4afUPZ4pTQE3k5xID\/4wQ+inKtB4+R17JlAmEC86oyv8okmEF9CJ\/AlwzeCHRjsxBgPE4jDoCQCmMtf\/Iaqvz5OKTe1pswOk2wlkyl+AnUs0rRMRpV+Bt0PmYyJGP+g+6xSv5MPxGnJFTac6qwFcuXYerr06ZTIGGd2mGhLIjLFt1MSr1GSboqmcom0m4yq\/VRRdlkZN6KS4Ri2iSHrq87vDWtj7969kWrNd0+b2wobTnWWQEAemFzGk97mp5TVZUENvZApvvUFP1GSTsqoeom0m4yq\/Yx1QsiISoZj2CZGrHgE9X7YcKqzBGJVeJAHSMT6yBTfWt5PlKSTsqleIu3FAnHqZ6wKLyMqJxkNq2V3KdGjxRsi99PGKk+yvs8EEpAPxM+yAcr7zekSatCsUJsPxE+UpKoF4nSJtBuBwIr54m\/zqLDx51RS0ZGwVCt+oJ32+SUjZDsZze+AQMYurscEIhkZJpAACETnsiFWCwTve42SlPlAtpdWUNF3GztOfDcCUbVidDCKGyHbyXh+x8NUdeZj0TQTiNoIMIEEQCA6lw06CERNFfSVklkgIBHjcbJi9EljX5NVRqvFEjSB+IkD4ZSGcq1ICh+IjmWD0w6CVx+IE+RBbqXKZLzzP\/9GR89WUtucDNrzzPflWhFACch4ds88yqpfQWnZ+WL73OzEjgeB4NfbfL8tuml3mA1kwxnJ1JRAC4Fcv36dduzYQThPsGvXLmratCmNGjWKhg8fTllZWUqSmE2zDkff85QsJtZlg9MOgjHpv7qSRTl3TqOMjAzXvjiRhFP92KIdkL2O8jLOuMaiyAD0YoEUP5Dv6gNR2TaWyWP3vRWDlJvyqPrrskjR\/Vd+GqgTVZaRzEwuTCDqI6yFQLZv307jx4+nwsJCQl6DAwcO0MqVK6lHjx40a9YsatiwoVQig0BeG\/MA3bZ7SaS8SrIYc+V+funtdhBSmxUqxYkYbbttY9rV\/07Vr4Rzc0L+xoj41lgU1b7IfCC55UupdcZpOl7ZjMpzxxKWMXaP6raxmWQgvywYD+XvvDCH7k77MNIsdrwQwGc4sfee605Dhw4NzIkqI5BEpjRct24d7dy5k4zQ+VatWtHy5cupffv2NYYp6Xwg586dE+QBqwMHfwyy2LZtG02YMEF81rt3b2UCeaUonbq1P0kpDVOp+mIVNcgZYJssxjq5oKRQUCiq8TgFhxnfGxMB7\/Vr8GLkPWx1Vp0ukcaJmGWwmuTmuBIrgaD+HaUVUbKi8V1Xf0h7Gj0tLAQrIaU27V5jt8hovzqtOZ2+6UFq2bJlDStpxYbFUX0DgXT4\/jTb8bBzuL7Q9T1BEObnjb8co8p\/VFDadxrS4ILdka921x9Bfzo\/WCyVjJ0eg5RANGayhAXSoOm3Idpet3FVyTUy1vPni3yidkuYqVOn0vPPPx+ZsIa\/xBzMhUnudmLXLn+HG2mZ8UTdyEdiTV9oN0hJRyDoEJYqS5YsoV69ekX6fPHiRZo8eTLl5OTQ7NmzKS0tzZVEDGB+\/2Q36pL3SaRs1YnG1KDJgKicD9bJdSqrD637LJV+0mIH5WWcjrxrKGn1pSqq3FsRyZ9pbPm++PfrlFHQWPw6P9hiB5knqbWNq4cbUv2MuyJ1yExykMRvDnSrQWpU3YJyBvyNKv57EFVXfttPCD3189H0XyfvJSwz+qW+SM0vvVsDM0NG9M0cSVuVVUDlaXfRn78eJiYqJn1ZZVP6v8NbqbDxF1G\/\/OaAObM1geXU\/+zbLyyVhV\/2p5fu\/SyKfMzCAE+QPEjEeMoqmwn58e6932lMIJ8vj5cK\/wtkMI+NtWNefCCyoDU7RQtrSkPIaiUnt4lizJOCEc9R9+6FUctRQ6\/xvmquWOkvu6RAzEuYtWvX0uLFi2uYXPCLwPqAT8TK+k7MunXVw\/Tje4ha5XxbAlZI1alKyvi3ntQgN1988c2Zj6PWzyogoI6qkzfqufrP7ZFXzBMAk7P+TXmR78xOvqv\/ezHyOeqglJOuMqRk3EVvfZleY+KgP6nNugsZUptH+1RKKjoQJiF+weEfME98lT6iDOpwe8\/6y79614lI1SBR40E93291ybGP6AcekIj1Mfphrk8mvxcCkQWtORGInQWS6JSGBoFgHqnOEyz1ht7RnHKb3ExF32lMbXMyRZfNaR7x\/16X\/7Ixsvs+ZgIBs2P9tmzZMsrNzY1qA99t3ryZVqxYQe3auQcvwSS9\/MUrfvrA7yQBAl4IRBa05oVAjAmcqJSGfgnkxdtPUEGjSteRV8kVG6vqaCEQO2Y3gDGvG52EtSpErJ3i92sfAl4IBL1TiSI2oxDWlIZMIA7OKS8EYjVJa5\/6s8SxIuCVQLy2F9aUhkwgLgSiuoRhC8TrdEi+8okiELv7VwyHq3lXBJ9t2rRJ+Pry8vJo2rQbO1nmaFW7KzKdrHPzCMKJ6tUHorKEUclWH6smxbyE0eVEPfzhjeP1ZgequXNwcFqf6stlUU4+uzIq7xw7V0lHz16OFIVTqmXVSRHMZjzlLbtTmyY3nJ7Wdo0yG07eK91tiHXArAFY1vqAgZN8RtkTKXdG+ts640zU7oiBoVMdcEanZKXaOk\/NssCRisdav1EGTlizA3bXnioat6JBoHEgS5cutYXfmnfUCCSLV0pDPwTy24EdXH0gqtnqY9XHmAlE9zbue79\/kHKvR++yuMVzeF0LAzDrO3bxDyjX6IP51O3KHtqdfidduL84EoBlF6Ox\/sit9OTfHxDjMTdvKQ1stE38W2Wy2ZEAdoOsu02Y3Nn3rhHyXz64Lorg6LZCapTXM3KqGGWuHNsQRbBi0mbcRTk\/eZuwfWsc0sO2r92pZAMnyHLt6zKxG4b4DyhnWvuLUXEy6AO200EYON2MrVzUjx0lbA93++Y9unqkVJAGdsP+0fLXdN+\/k2gX9a7fVEHztv0zMAKJdaKE5f2kiwMxAskQxGSOOvUbSLZ69WqRA9IPMfgdZGsEphGp6XYIzSqfuY4xX62kMedXCnGy+\/2YUm7eHxENMRhXjm6InEK1k9kgTNlug3Hqt95thVT572NrBJJ99JdlUcFqX+88TTf\/aJH0EiZVHA0MUN4tJQK+B5Z2hMz5QFTRvlEu6QgEndqyZQtNmjSJioqKaNCgQVRaWuo7lN0gEG+wxl7a\/ItsRFHis62fn6aC5qk0c0An6VkYlL+wdSWNKp0VEeizrhNpdxaJfBw9C38kLASZ09gcxapCpE6h7OZcIPiVL7r1P5QuYYodzZo12JH0T1psjwTDBe0DCaJPiagzKQnE7jDdsGHDxIG67OxsJZzDBowhtOykq7Vz1tQCm7IeoFk5xaKYcZBNRiBes4a5yWhHjEoDorkQZJz5p4P02RmK5DYx48AEogZ42OZJzD4QtW7LS4UNGL8EYk0tMDOnmP6cdcM38vDdLYQfxWnXCQ7S9DYPOmZHc0LRK8nJR8NfCbeTvKoZySZOnCgOZfJjj8Dx48dpypQpofEVMYFINNXP5DR8E7vTC2hEef9IC+ZkPnZOTqfEzrLJ5EdGWZ1ev5ed5HXLiQpH6qnKPHp2WRmVlJR4bbrOlS8oKKCXX36Z8vNvHO1I5MMEEgCBmKt0W0LInKSqihEGApGlTlSRsaysjPDn9ak6U0LfnP6YGjTrTqlN\/VsvOEphpFiEDNj1gkM7TA9wTE9Ppy5dukh9cvGQmwkkYAKRDaKKk1RWh8rklNUR6\/d2TlIs2YwnDDLK+qiL0GXtxPJ92HBkAkkwgcSiTGGbnG7WVqIVXyV\/iJF2seG1L6Tb0jrGzU8dicbRKjMTCBOIHz32\/E4iFV81f0giZVQFNGwyMoHUAQLRkec01joSqfiq+UMSKSMTiCoCDuWSZRsX3VMxl2OEK+p1nUmV7eSS7bCo9CWRk1PVt5FIGVUwRJmwycgWiGYLRNVcVlUYlXKypMogAOMxYlFU6jXKyHZYVOpKtOKrOKsTLWNtwJF9ICqjZCrjValUzWWPYrgW92KB+LlYSrbDotIXrziq1Km7DMvoHVG2QAK2QLyGpXsfQrlZqyOcPdY6eHL6Gdma74QNRyYQzQRi+EBkl3brUacbtYRNqez6xjLqGfGw4cgEEgCB6FEV9VrCplRMIOpj57Vk2MaaCYQJxKsO+yofNsVnkvM1jDVeYgJhAtGjSYxjncSRCYQVv04qPlsgeoadCYQJRI8mMY51EkcmEFb8Oqn4bIHoGXYmECYQPZrEONZJHJlAWPGVFD\/W8z28C6MEs7RQ2HBkAmECkSqtjvM9YVN8XsJIh12pABMIE4hUUXSc72ECkcKsVCBsODKBMIFIFVf1OLxbRWFTfLZApMOuVIAJhAlESVFUjsMzgShBGVOhsBExEwgTSEwKrfpy2BSfLRDVkXMvxwTCBKJHkxjHOokjEwgrfp1UfLZA9Aw7EwgTiB5NYhzrJI5MIKz4dVLx2QLRM+xMIEwgejSJcayTODKBsOLXScVnC0TPsDOBMIHo0STGsU7iyATCil8nFZ8tED3DzgTCBKJHkxjHOomjK4FcunSJVq1aRatXr6bjx49T06ZNacCAATRmzBjKzc2NAHb9+nXasWMHLVq0iHbt2iXKjRo1ioYPH05ZWVlKwCbT1ZZKHdZYiKM89YDJOHrH0ZFAysvLaeLEiXT48GEaOnQode3alfbt20crV66ktm3b0sKFC6l169aixe3bt9P48eOpsLCQhgwZQgcOHBDlevToQbNmzaKGDRtKJWMCkULkWIAV3z925jcZR+84OhLIpk2b6KmnnqLFixfTfffdF6l5\/\/79NHr0aHrkkUfo8ccfp4qKCkEesDrmzp0bIYtt27bRhAkTxGe9e\/eWSsYEIoWICcQ\/REpvMoEowRRVyJZAqqurhYWBZcnSpUspJycn8hKWNdOmTaOUlBSaN28e7d27VyxVlixZQr169YqUu3jxIk2ePFm8O3v2bEpLS3OVjgnE++AZb7Di+8eOLZDYsPPsRD1\/\/rywPJo0aSIIZOPGjcJKWb58ObVv3z7KLwLrAz6R119\/PYqE7ERmAvE\/kEwg\/rFjAokNO88EgqUJCGTKlCnC8pg\/fz7t3LmTli1bFuVYhVj4bvPmzbRixQpq165drbRAvvzyS+HPGTFiBOXn58eGdkBvs4x6gGUcvePoiUAOHTpE48aNE0SBJc4tt9wiSALWg52VsW7dOrEzY7VO3CwQOG7hjA3Lg90nkGXY5DLjwzLq0ZZkwTEvL4\/wF49HmUCOHTsmfBpnz56lV199lW6\/\/XYhny4CKSsrExO1pKQkHv3mNhiBpEUAP3b4i8dTr127dtfNDY0dO5aKi4uj2j548KAQ6Nq1a8Lf0alTp8j3bgTiZQmDCkEi+OOHEWAE\/CMQVwvEjUAQIIYYD1geWP+\/9NJL1KZNm6ierV27VosT1T9c\/CYjwAgkCgHXJYw5QGzOnDk1nKQQGv4PHdu4iQKA22UEGAH\/CDgSiOEw7dy5s4jjwLat3XPu3DkRSNayZcuoqFOvgWT+u8BvMgKMQKIQsCWQK1euCNJYs2YN9e3b13b7EkuZfv36iQCxLVu20KRJk6ioqIgGDRpEpaWlnkPZEwUAt8sIMAL+EbAlEJyDeeyxx8TZF6cHxIJAsszMTLI7TDds2DBxoC47O9u\/dPwmI8AIhBoB5W3cUPeChWMEGIGEIMAEkhDYuVFGIDkQYAJJjnHkXjACCUGACSQhsH\/baNBJm7AVj+C\/3\/3udzUONOpIBJVg+CLNB4Hj0aNH6YUXXiDsKF64cIE6duwo\/Hrw\/6Wnp0faTiYcvY4nE4hXxDSWDzJpE5T6ww8\/pOnTp4tzEXZnlXQkgtIIh++qgsAR8U2Iyr7tttto8ODB4twXDoa+8847gkR++ctfRkgkWXD0MwBMIH5Q0\/ROUEmbcBwA55XeeustqqqqooKCghoEYsTvxJoIShMUMVWjG0eEMcyYMYOOHDkiDoO2aNFCyAdSRopPWCU4Yd6tWzdKJhz9DAITiB\/UNLwTVNImHHZExjgkeho4cCChHRxTt1ogyRJBHASOX331FT355JOEIErruTA\/g4pfAAAB8klEQVQEWI4cOVJk24Nlkiw4+lVpJhC\/yAX4XixJm0AgCxYsoP79+9M999wj\/m2XbqEunGGKBUen4f3000\/F0Y1nnnlGEEhdwNFN1ZlAAiQCv1XrTNrkdFpaVyIov32Mx3s6cYS8WA7iQOmGDRvEEgan0usCjkwg8dBWTW3oTtrkRiA6EkFp6rb2anTjCP\/H+vXr6emnnxZJtZ544glKTU3Vlg9HOwBxqpAtkDgBrdJMEEmb6iKB6MYR5PH++++LhFf3339\/1KFRXQm1VPQjjGWYQEIyKkElbfJDIF4TQYUEQiGGbhzhpMWSBT6PPn360MyZM6POd+lMqBUmHFVlYQJRRSqgckEnbXJS8GRz\/gWBI7ZzsR2OILxHH31UxH7g8Kj5STYcvao5E4hXxDSXDzppkxOBJNv2o24c4TBFDAgIBEsXbN3C52F9kg1Hr+rNBOIVMY3l45G0yYlAkikRlG4cDYfps88+KyJ5cbWrHXlAFZIJRz+qzQTiBzUN78QraZPbGj0ZEkEFgSMC73CBPMgBSbPs7nZG8qy7775baEIy4OhXpZlA\/CIX43vxStrkRiDJkAgqCByxHEKwmNuDWxcRSIYnGXD0q85MIH6R4\/cYAUaAmEBYCRgBRsA3AkwgvqHjFxkBRoAJhHWAEWAEfCPw\/4krfU8RdkACAAAAAElFTkSuQmCC","height":164,"width":272}}
%---
%[output:741b5124]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARAAAACkCAYAAABfJnHCAAAAAXNSR0IArs4c6QAAHsBJREFUeF7tXQtsVtWWXkjT0lArD4vlIZY7wBVGRUCFYRohXBMTuIBcU0HgqlDUS3jIRQkvFVsHlKliEBTq8L4OBdHrYCOM8RHQQiRoDIEpDvQGBuRlgSKUYE2Rybfl\/Jz\/9Dz2f84+j\/\/\/10kaSs8++\/Httb6z9trr7NXs6tWrV4kvRoARYARcINCMCcQFavwII8AICASYQFgQGAFGwDUCTCCuoeMHGQFGgAmEZYARYARcI8AE4ho6fpARYASYQFJYBr7++msaM2ZM3AhXrlxJgwcPjvvb\/v37afz48XT27NnY359++mmaNWsWnTt3joqLi2nv3r00bNgwevXVVyk7O1sKNWzwHTx4kNatW0dfffUVHT9+XDzXpUsXGjRoEI0bN44KCgqoWbNmUvVxoeghwAQSvTlR1iMzAvnLX\/5CM2fOjFPaTZs20Zw5c+La9Uog9fX19MYbb9Df\/vY3amxsNB1TRkYGTZw4kaZOnSpNSsrA4YqUIMAEogTGaFZiRiCFhYW0bNkyys3NFZ3+5Zdf6MUXX6T33ntPGYGAMJYuXSp+ZC4QCH5AKHwlFwJMIMk1Xwn1Vk8g7du3pzZt2tCpU6dozZo1dMcdd4i6Tp48KZYotbW1glQOHz4s\/u7FAqmqqqJJkybRpUuXBClMmDCB\/vznP1OHDh1E3adPn6aKigpavny5sE5atmwpfge58ZVcCDCBJNd8JdRbPYH06tWL7rzzTnr33XfplVdeoVGjRom6UOaxxx6j\/v37i2UNfBVeCMRo0aCtRx55pImfA\/4RWD3a0gllSktLKTMzM6ExcuFwEWACCRd\/X1s3Esjo0aOFwuqVFcuZxYsXC1\/EhQsXYksZtxYILBnUtW\/fPkFY5eXllJ+fbzpOWENoRysLB29eXp6vmHDlahFgAlGLZ6RqMxIInKf4uemmm2jVqlV044030owZM+izzz4jKO+ePXuEwnuxQA4dOiR2dE6cOBFHVGbANDQ00Ny5c+nDDz8Uyxssrbp16xYpDLkz9ggwgaSwhOgJZOjQocL6mD17tli2rF+\/ntq2bSuUHcsJKC8UWSWBaFaMHcSLFi0SbTKBJKcgMoEk57xJ9VpPIIjhWLBgAZWVlYmtVVgePXv2FMsNbWcGjkyVBDJy5EhauHAhZWVlmfb38uXLgtAqKyuZQKRmNHqFmECiNyfKemQkEASBbdmyRSwbEEzWrl072rhxo9gheeGFF2jFihXCH+JlCaPt6nz\/\/feOPhBjWfaBKJv6wCpiAgkM6uAbMiMQKDZ2XeB\/aN68ufgXVgIcrPqAMrdOVNQHMnr\/\/ffFgK1iPIyxIrwLE7x8qGiRCUQFihGtQ08gGiHo3\/roNmIw4A\/p3bu3EgJBnTt27KAnn3xSxHggDqSoqEgslW677TaBFPrwzjvviFgQjgOJqPBIdosJRBKoZCxmRiAI7oLf4eOPPxZDuv3228WODALNnCwQJwy0+BKORHVCKnXuM4Gkzlw2GYkZgaCQFvuB37E7A98ILBGQCpYcVj4QJ6j0AWpwkCKUHX4N\/hbGCbnkvc8Ekrxz59hzKwL54osvxJICF3ZjpkyZIn43K6\/\/GtepQT2BoKz2NS4iTrdv3x4Lk9e+xoXfo3v37vw1rhOwEb7PBBLhyeGuMQJRR4AJJOozxP1jBCKMABNIhCeHu8YIRB0BJpCozxD3jxGIMAJMIBGeHO4aIxB1BJhAoj5DKdq\/RZ8cpqPnfqbObVrQrAe7pOgoU39YTCCpP8eRG2HFnlM0ueJArF+zHixgEoncLMl1KDIE8sMPP9AHH3xADz\/8MHXq1Emu9xalELh08eJFcd5FkOdshtUuYEimtkEeIBHtevTefHrr0R6u5jyZxu1qgBGTcWN3IkMgWhDThg0bxPF6Xq6ff\/5ZfG+B8OwWLVp4qSqhZ8NqF51MpraNFgjIAyTi5kqmcbsZn9UzYY5b3ycmEIWzGuakJlvb8IFU1Zynwq6tPC1fkm3cqsQtzHEzgaiaRUM9YU6qVdtBOCujOG6fpjiu2nQdNxOIT9IVNYEKylkZtXH7NL1Nqk3XcTOB+CRhURMolc5KO8iiNm6fppcJxARY9oEolLaoKZJKZ6UbAuHlk0LhitBymS0Qn+ZVJYFgWxs\/she2M8+cOUOtWrWK23la9MkR2vmP8\/Sv\/wRnZYFsdQmVM2sbbaJt7cIui9udFrvOWI07oQG4LBxW2whzuPnmm0PZaTRCxRaIS+Exe0wVgYA4kL9l9+7dCnvHVaUKAv369aOXX35ZnGkbdKgCE4iPUqSKQLSYGKRg6Nixo4895qqTDQG8VJYsWUKrV68WZ8wygVybQQ4kuy7KKrFINgXh\/tojoMkGE4gBJ5VKo8oSSFSYVbWrEotEx8Dlo42AE4HUvVdCjbVHKCOvgFo\/Mt\/3wbAPRCHETCAKweSqTBGwI5CL29dS7VvjY8+1LnrJdxJRQiA4PHfnzp3iFG4kaEbO1QkTJogERjjtW+ZS+dZVpcgy\/daXUdWuSixkx2B2eLJVvlp9Skp9\/XbfMWlJt1HeSxJtrW3Ug9Pks7OzbYdo1tdevXqJVBZt2rSRhSeunD79hXbDLg+wHluZb7208rNmzWryXZgdgYA8QCLadeOgJyhv8hpXY5R9SAmBVFVV0eTJkwneYWQ4q66upnXr1tH9999PJSUllJOT49gflUqD+IPqo2eoe4dWNG9ocNnek5VANOVGSodRo0bF5kpTFP1p65pw44NHCLh2afNnpUhIoo0ytbW1hDy9+mcdhUNXIBECseqT2bhk+6CNQ09AGiaow4yY0B5ernl5edS5c2db4tMTnhnZJGKBgDxAIn5engmkrq5OkAesDgiaRhbITjZt2jTxtyFDhjiOAUr\/xn\/8J\/Xr1586dXK384DDaTq3yY47awLxD\/i73aXdxwE3Xq\/6+noRh+HlGIGMM\/9LlW\/MIJm3ldf+4nkoxYkTJ0wF23gPyoB8ulaKYnZPUzC8XI4cOSKIxK0FIEsgGilakZWbF5ZW52uvvdbEMrC6p\/UX1lxBQYEgEisLTKsDc4HLiUA6H\/wvyr58hlq07yaWKvB\/1G1+KSYSSbGEwURgqYLEzEjYrF1QJOQcgZlYWlpKmZmZlrJujJhUoRTJXAcIJKfq3wMlEFml1t6miSxD9M9gXsaPH09mSigzZ7IEYkd0aEe2Hn2frCw1u37riQo5cIqLi4WVrrf08Lye8P70pz9ZYqTV99bI39Pvj\/53rOkW\/zyIfv6f7XFdSYolDN44yHRmFCj4RWB9wCfi9LYxfrMhI0ipXCZoAtEnlHKyevRvSWMiKbM5MSoqyiC1Ji4ZH4axThnFlymDehMlQ\/3yApaNTP+NSx6zJZBxjHaWjjZXr3c\/Sb1utLeYk4JAAMiuXbtECkOs8fQX7m3dupXWrl1LyEZmdTGBxCOTKIGo+ObEzDFopSR6wtF6buWYNFOGRBVXj44MOVj5aYzy56YfVpn6zIhXv3TTLA6ZpVPaEYiV+Ss7QSqXMB1yM+jEhcakNlASIRA\/PtkH8ZeXl8cwtLM0jMRjdKKavXHN\/BNWDk+jgxK7Lk4WjBcCMY5dZmdJ81kYSdRM\/rW+2TlTmUCuiZ4sgaD41FVVtGnHfvFkfv4tdMstzkfcnbzYGEcWfTu2IBBI5YH6mPDj\/+1vzLAlFDyH69vj8SahGRnZEdSVK43U2HiFzv3S3JLEtLb0HTKOo\/Mv\/6ALWxdK+UD8\/GRfE3bsnDj5PMx2D5zy6uoVDgQC5cXHgK+\/\/npsi9X4dxkCkbFSzJYwRoexnSKbWTNz5swhjUCttrq156y2yHFfFYEkxS6M3ZpOdgkD0PTff\/Tt25fy850JZPN3tfTM5prYXC4p6ip+N\/6tqHf80sov8wRCc\/r0adH38t114uBg\/c4OdoQ2F5sfHrxkx8nYEX\/9sv5PxNE4+SMwDq+f7DspiawTUC\/42naw3QvEaMrj\/xUVFeIg7D\/+8Y+xXQ7IULt27ejLL78UxCJDIBo5WO0W4b6Vb2bAgAFxDk60jwvbzna7Vcb69u7dS2PGjDGdQydnrAoCCWIHBrh43sZV4UTVE4iM0ugJwOxszQUfH6Lt35+hB+7I93TeZqJEY4wDcavcMutkJwxk++70tkZfnnvuuZgFYvfC0As+rAu7pYbRlIfCgUDuvvtu+vHHH4XCogziiB599FGxpEqEQJyWMWYYa8sxq1gWu3nRb9dqZGO1tNfKHj161HSDQYZA3p7wB+p+cR+1uBZIl91zkJjyy9XbCb8HEcauhEBUbON6IRAzRVEV0CWrhFo5s3bdHB6cKIEk2k9jeU1gYfnpdxbM3pRWQVNGhXV6yxqXEGfPnhUEAstr8+bNNG\/ePAKpIJ4I25oLFy5MiED0MmUkBLtAMuOyy8ynAzIzvuj0xIq2sV1rDLbT4243xzIEkjIf02mBZPisWB91mmggmUqliRKBuFFulVjItm+2ZrdbpxudjWhHr1SJbFdit2fgwIGCQObPn0+LFy+mxx9\/nP7+97+L4Ks+ffo0IZDKykrToRkV3ksou\/ascVljtgulb1fG92fnTE0rAsEsbtu2jaZPn06FhYVUVFRENTU1oYayM4HI0kZ0ymk+EFhAH330kQh5\/\/bbb2nu3Lmik5oF4vb7FaeRasslEJi+DZABomfdht47tZvofbtQ9kTrUlHesw8EnTD7mG7s2LHCHM3NzZXqp8q3LhOIFOSRKqQnEJzIhmhVbUmF\/\/tNIABD7zDF\/83iOMIGLSUJRAWoTCDXUVSJhYq5CaIOPYFg6aAP+YZJHwSBaCQiGwMTBC7GNphALFBXqTRsgYQh2txmEAgwgTCBOMqZSjJ1bIwLJBUCTCBMII4CywTiCFHaFmACYQJxFH4mEEeI0rYAEwgTiKPwM4E4QpS2BZhAmEAchZ8JxBGitC3ABMIE4ij8TCCOEKVtASYQJhBH4Q+DQPhUdsdpiRXw61R24+cBZh\/1MYEwgThKatAEwqeyO05JrIAfp7KbfZ1r9a0MEwgTiKO0Bk0gfCp7f8c5QQG\/TmW3qhdyoD9KAX1gAmECcRTWMAgEbTodfo2Oy3xpahwgn8r+28FCdqeymwmFmRwwgTCBRI5ANKFEx5wOdOJT2WcTjhLw81R2TUDMyJoJhAlEOYGoSKjMp7LLZTAM4lR2\/XLJmBiLCYQJRCmB+JFQmU9lv34ifRinstulymQCYQJRSiB+JlTmU9kPSWfRM56n6vZUdqc8u0lLIDjoFrlucfIYznrUX2YHCuEwIaS8bNmypaPC6L3LTmtwmcrS6XN+owWS6FH+fCr7bxnycOEQ6DBPZdfmAgnarBzaSUkgx44dE0e6ofNmCl5VVSUSbPfr10\/k\/ayuruYjDVvYJ\/S2I8JEd2HgA3F7Gjefyr4h9kIM81R2q4OtjXKSVARy6dIl+vDDD2nJkiWEU7PNvPTaocpt27YVuXBzcnJEOT5UOTgCkbHK7MrwqezX0QnjVHZte9cuW53Ww6QiEM2Z1rVrVxo+fLg4LdtogWBAWKqsWLGCBg8eHJuJ+vp6mjFjhjigtrS0lDIzM23lPNG3rl1l6bSE8Uoe2vN8KntTJIM6lf2ee+4Rp9FbXfrUoklFIEuXLhUE8NBDD9G+fftMM22FnVjKDHQmEFW0Elw9UAykdeBT2e0xTyoC0Q\/FykKAlbJr1y5auXIlwfmjv9yktmQn6vVwZRVYBEcB3lrSEwifym6NZUoSiFUYdCJhzxowzzzzDI0YMUIqN64VzLBATp06RbfccovIpRrUpardb775RiwL05VA+FR2ZwJZvnw53XbbbaYynpGRQfgJ4pLOC2NngagkEAwaWcmgQG6vhoYGkZgIFlFWVpbbahJ+TlW7SOkI\/1E6EUjCYKfpA5oewi+CzIFmMt6qVStq3bp1IAg169Kly1V9S4kmFrZLYehmCVNWViYSCiHDvdsLb7DTp0+LOlp42E5NtH1V7UJIEEfDBJLoDKR+eY1A3nnnHerUqZOpjAdqgXglEHaiXhdaVc5blTtSqa9S6TXClPOB8DYuE0h6qXC4o005AtECydq3b08lJSUcSHbyJAELL0snvUMZ0b18MQIaAsePH6eZM2fS6tWrhRPVq6x5RdazExUd2LZtm\/hGprCwkIqKiqimpoZD2T34XrCNCSHZvXu31\/nl51MQAbxUXn75ZWrevHlqEIjZx3Rjx44VjsDc3FypKVS57lfli5DquK6QynZBIviRvdD2+fPn6eabbw5sC0\/rG7cdLOZwnmKeTyqwdmXly6qctAXitSGn55lAnBCyv6+SvBLtCbftbcmaKN4oHybm+v4ygbiZPYtnwpxUbju9lDjM+WYCUUga+qrCnFRumwnEJ7G2rZYtEIWosxKzEisUJ9uqwpQ1tkB8muUwJ5XbZvLySazZAgkKWFZiVuJ0kDW2QHyaZSYQJhCfRKtJtWHKGhOIT7Mc5qRy20xePok1L2GCApaVmJU4HWSNLRCfZpkJhAnEJ9HiJYwTsByJ6oSQ\/X0mLyYvbxLk7mmOA3GHm+lTrMSsxArFieNAEgGTLZBE0GpalsmLycubBLl7mi0Qd7ixBaJDgMkrvciLnagKSUNfFStSeilSus43EwgTiFIE0lWR0nXc0gSC3Ljr168Xp4PjKDXkv0WWuqeeeiouiZTZgUI4TAipGVq2bCklrOwDkYLJslC6CjOPO3irT4pAkFcFSZ6QeHnMmDHUp08fkd5y3bp1hCTAb775JnXs2FHUVVVVRZMnTyYctTZ69Giqrq7mIw09HGnohkpYkYJXpHTFXIpAKisrae7cubRs2TIaOHBg7Jn9+\/dTcXExjRs3jqZMmSKO0QN5wDpBspucnBxRdseOHTRt2jTxtyFDhjjqBFsgjhBFdlsvXRUpXcftSCC\/\/vqrsDB27txJ5eXlIsG2dmFZM3v2bLrhhhtEImRkUcNSZcWKFTR48OBYufr6epFdDc+WlpZSZmamrQIwgTCBuEEgXZU4zHE7EojdRF64cEFYHkidBwLZsmWLsFLWrFlD3bp1iz0Kvwisjz179tCqVaviSMisfiYQN+pz\/ZkwBYrbTq\/lkycCwdIEBIK0A7A8kL5y165dtHLlyjjHKhpxk9pSRTrHsAQ6rHaBNbedXkoc5ny7JpCDBw\/SpEmTBFFgidOuXTtBEiqTa8NxO2LECE+5cQHuqVOnTDOXe3vP2z8dVrsagYQxZm47eDlzwjzQ3LhXsdaQuI4dOyZ8GufOnaPly5dT9+7dxVOqCQR1Pv7448K6cXs1NDQQdpHMMpe7rVPmubDaRd+47TzKysqSmSZlZaKKeatWrYSLIYirmUxy7QMHDogt3StXrgh\/R48ePWJ9syMQN0uYsrIy6tu3rycL5PLly3T69GnTzOV+ghpWuxgTt53vKZ2oG7mIKuaBWiB2BALjBDEesDwKCgpo8eLFdOutt8ZhvXHjRnaiXkMkzHUpt80+EDck6PUZ24\/p9AFiCxYsaOIkRePwf\/A27m\/TwErMSuxVIWWfD1PW9H20JBDNYdqzZ08Rx2G1pqqrqxOBZMgSXlJSwoFkIeUrDVOguO30Ik5HAoFzCKRRUVFBw4YNE8sX44WlzPDhw0WA2LZt22j69OlUWFhIRUVFVFNTw6HsHMou+zL1VI7JK3jyciQQ7GBMnDhRfPtidYFYEEiWnZ1NZh\/TjR07lvBBXW5urpSAcCCZFEyWhViRglekdMXckUC8ibK7p5lA3OGmPZWuwszjDp44mUC86SpbAQYEWImDV+IwMWcCYQJRikCYwsxtB09eTCBK1ed6ZSzMwQszYx485kwgTCBKEWAlDl6Jw8ScCUSp+rAFEqYwc9vBkxcTCBOIUgRYiYNX4jAxZwJRqj5sgYQpzNx28OTFBMIEohQBVuLglThMzJlAlKoPWyBhCjO3HTx5MYEwgShFgJU4eCUOE3MmEKXqwxZImMLMbQdPXkwgTCBKEWAlDl6Jw8ScCUSp+rAFEqYwc9vBkxcTCBOIUgRYiYNX4jAxZwJRqj5sgYQpzNx28OTFBMIEohQBVuLglThMzJlAlKoPWyBhCjO3HTx5MYEwgShFgJU4eCUOE\/OUJ5AjR46IQ52R4c7sQGil2qOrLKx20QVuO36u694rocbaI5SRV0CtH5nvy5SnK+aRJhBkwOvXr5+nCT9+\/LhI\/q2irkQ6Ela76CO3fV1ubqr5nLK2PB+buksDJtKlfylOZCqlykYV806dOhF+grhsE0sF0QGtjR9++EEo\/e7du4NslttKQQRmFtTSg23rYyP75GwOlR3JS8GRmg8JL078BHFFhkAwWJAIfvhiBLwgYLRAGkb8G\/3U9Q9eqkyqZ9PSAkmqGeLORh4B+EAuV2+n7J6DfPOBRB6EADoYKQskgPFyE4wAI6AQASYQhWByVYxAuiHABJJuM87jZQQUIsAEohBMrooRSDcEmEDSbcZ5vIyAQgSYQBSCyVUxAumGABNIus04j5cRUIgAE8g1MC9dukTr16+nDRs2iLDwtm3b0kMPPURPPfUU5eVdj2K8evUq7dy5k5YuXUp79uwR5SZMmECPPfYYtWzZMm5q6urqaPXq1aLey5cvU\/\/+\/enZZ5+lu+66i5o1axYrm0idCudeVOXHuPV9rKqqomXLltHbb79Nbdq0iet+qo376NGj9Nprr9GOHTvo4sWLdPvttwvZGDZsGGVlZUVivlXLDxMIEdXW1orQ30OHDtGYMWOoT58+tG\/fPvFBXufOnenNN9+kjh07CuyhEJMnTxbf64wePZqqq6tFufvvv59KSkooJydHlKuvr6f58+fTl19+SZMmTaIOHTrQ5s2b6ZtvvqG33nqLCgsLY3MpW6fqyfdj3FofQQ6ff\/45zZkzR3yXsWrVqiYEkkrj\/vrrr+npp5+m3\/3udzRq1Chq164dbd26lT766CNBIn\/9619jJBLWuFXLD+pjAiGiyspKmjt3rnhTDhw4MIbz\/v37qbi4mMaNG0dTpkyh8+fPC\/KA1fHKK6\/EyAJvnGnTpom\/DRkyRDxvVidIBQr1008\/0ZIlS6h169YEK0W2TtUC4Me40Ud8jrB8+XJBmI2NjdSrV68mBJJK425oaKAXXniBDh8+LCzT\/Px8MVUgUVifsErWrl1Lffv2DXW+VcsPEwgR\/frrr8LCwLKkvLw87i0J83727Nl0ww030Kuvvkp79+4VS5UVK1bQ4MGDY\/MBYpgxY4Z4trS0lK5cuSKew7Jl8eLFMaLBA1988YWwdiBYvXv3Jry5ZOrMzMxUOv9+jBt9PHfunCBdYDVy5EiBLz57N1ogqTRuvBCee+456tmzJ82aNStung4ePEjjx48XLxhYJmGNW6nw6CpjC8QG2QsXLgjLA5YCCGTLli3CSlmzZg1169Ytbk0L6wM+ESgKCGTixIk0YMAAR4HauHGjVJ1G\/4FfAoF63Y4bfQSBlJWV0YgRI+i+++4Tv0NpjASSauO2mo\/vvvtOvCCef\/55QSBRHLcXWWICsUEPSxMQCI4ZgBAsWrSIdu3aRStXroxzrKIK3MOaF6YqzHa8dfAs\/CT6C34W3MPbGQ5V2Tq7dOniZZ4TetbtuM36iPGZEUiqjxuAQw5ggX7wwQdCLnr06BHJ+U5IOAyFmUAs0IPpCecndmCwxIFTzEoZUMWmTZvE+hfWCS6QxNSpU8VbR39pJj52ZGDuytapt3i8TLjTs17GbdZHOwIxIxYjlsk6bvg\/3n\/\/fZo3b56QI8hCRkZG5ObbSR6c7jOBmCB07Ngx4dOAssMZ2L17d1FKVtmTlUC8jjtZCUT1uEEen376qbBcH3jggbjdOVkZCoo4nQjC6T4TiAGhAwcOCCcn\/Bjwd8Ds1C67yTdbwphZIGZLGKs3sb5Ov5cwKsad6BImFccNpzGWLPB5DB06lF566SXKzc1NWIb8nm8nYpC9zwRyDSm8NbA\/D8sDBzFj7XrrrbfG4SjrAEsmJ6rKcZs5eq1IVxZLv5zHfowb27mwWBE098QTT4jYj+zsbFcy5Ne4ZYlBthwTyDWk9ME9CxYsaOIkRTHZLbhk2MbVBETluM22mq0IRBZL1dvXfo0bDlP4wEAgWLrABwafh\/EKe9yyxCBbjgmEiDTHIfbxEceBbVuzSwt+at++fdy61iqQDGasPurULpBMpk7ZSZUt58e4jW1bEUgiWMqOR7ac6nFrDtMXX3xRBAoimtmMPNC\/MMcti08i5dKeQGB2gjQqKirENwtmeWSwlBk+fDjhbbht2zaaPn26CEUvKiqimpoa01B2CAreRAioQn6arl27Woayy9aZyMQ6lfVr3LIEgnKpMm4EyuGbKcw55ET7nEGPBeTl3nvvFX8KY9xO8uD2ftoTCL4HQdAXvn2xukAsCCTDetbsA7CxY8eK7x30zjLtbaP\/mA6hzPCx3HPPPY4f01nV6Xaijc\/5OW59W3aO50SwjPK4sQxEnJDdhUBDbUs\/jHGrws9YT9oTiF\/Acr2MQDogwASSDrPMY2QEfEKACcQnYLlaRiAdEGACSYdZ5jEyAj4h8P\/AtG58sQvJdgAAAABJRU5ErkJggg==","height":164,"width":272}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:455d64b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:2f3064ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82c04625]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2551645f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e462a53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4603b234]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2cdbab0d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6ba43a7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:292005b6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44587408]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d2c25dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a4e322b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7209c994]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4aaaea36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15765720]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3cccc1c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:36b97662]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47d109dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d6aed90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:719ba86f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60041ab7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8288b7ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47a2a2cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2fe9d7c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3e1bed74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6fd93c27]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7ed1c1f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:3e50198a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d5d30b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:783e6de4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a094f46]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1ec23811]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:7050df14]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7da68526]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:80d45bc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:66e98f44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:691c68b4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:1847470b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b01f108]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:77ae7dc0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:763e72c7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:41e0e013]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:4767e172]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:091c8eab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0375b2ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:308a38e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:57970a44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:0a9eaadd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:83eefd76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:243b56a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d25fdd6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62224dce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:71005745]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:83d93634]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9c7cc22d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:64150d21]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:728da710]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1190d817]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:381adbb8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85b092d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1ae00d2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9b025613]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:46dc88b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8dec3a57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:952f7d45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:11939650]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:11838473]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
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
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2ef1c20e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76107dbb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:07c20a40]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1b7752e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:730d9a8b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1f45c1a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a727372]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:62e057cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:378fafb1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:45651540]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:18518a12]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:500ad4dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26dd33b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57de100a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4aaa3564]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:55bda12e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0bc04a08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:90cfd8ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5e1643bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3ff75fa2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8036573e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:385019d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c1aa960]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ba62846]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6937e54a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8174113e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3fd50b0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5413a336]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3fbb82da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f138eaf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:26350180]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5cf78d39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2297ef76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4db7e9e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6216ff5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3372c653]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2684268e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:58083934]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2698c44d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1995cd3e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01fdb1b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:57ac2e2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:76244fe9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87ecc18b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:35966a80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:953d2ac2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:31a1500c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99b1044f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:019ba41b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:6830611f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:14818eb6]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"2211×12 table","name":"MLO_result_MK","rows":2211,"type":"table","value":[["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'y'","'MK'","0","-6.2725e-04","0.0028","-0.0042"],["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'MetSea'","'MK'","[0;-1;0;0;-1]","[0.0205;-0.0071;0.0009;-0.0010;NaN]","[0.0464;0.0032;0.0077;0.0060;NaN]","[-0.0059;-0.0182;-0.0063;-0.0084;NaN]"],["'MLO'","2025","10","'daily'","'Ba3_A82_ae33'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'y'","'MK'","95","-0.0069","-0.0034","-0.0104"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'MetSea'","'MK'","[0;-1;-1;-1;95]","[-0.0032;-0.0075;-0.0008;-0.0053;NaN]","[0.0252;0.0028;0.0063;0.0014;NaN]","[-0.0278;-0.0186;-0.0083;-0.0122;NaN]"],["'MLO'","2024","10","'daily'","'Ba3_A82_ae33'","'abs'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0100","-0.0040","-0.0160"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[0;0;-1;-1;-1]","[-0.0253;-0.0047;-0.0078;-0.0041;NaN]","[0.0106;0.0095;0.0001;0.0039;NaN]","[-0.0629;-0.0192;-0.0159;-0.0122;NaN]"],["'MLO'","2025","10","'daily'","'BsG_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0110","-0.0083","-0.0138"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[95;-1;95;95;95]","[-0.0429;-0.0034;-0.0042;-0.0055;NaN]","[-0.0249;0.0025;-0.0006;-0.0014;NaN]","[-0.0622;-0.0093;-0.0078;-0.0097;NaN]"],["'MLO'","2025","20","'daily'","'BsG_S2S20'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MLO'","2025","30","'daily'","'BsG_S2S20'","'neph'","'y'","'MK'","95","-0.0058","-0.0048","-0.0068"],["'MLO'","2025","30","'daily'","'BsG_S2S20'","'neph'","'MetSea'","'MK'","[95;0;95;-1;95]","[-0.0218;-0.0030;-0.0041;-0.0017;NaN]","[-0.0132;-0.0004;-0.0024;-0.0003;NaN]","[-0.0307;-0.0056;-0.0059;-0.0032;NaN]"]]}}
%---
%[output:793e416e]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"19×19 table","name":"MLO_result_LMSlog","rows":19,"type":"table","value":[["'MLO'","2025","10","'month'","'Ba3_A82_ae33'","'abs'","'log'","'LMS'","3.1262","95","-0.0726","-0.0261","-0.1190","-5.1766","-1.8649","-8.4883","-0.5160","-0.5154","-0.5166"],["'MLO'","2025","10","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","4.8948","95","-0.0465","-0.0275","-0.0655","-34.7014","-20.5226","-48.8802","-0.3717","-0.3714","-0.3721"],["'MLO'","2025","20","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","5.1774","95","-0.0239","-0.0147","-0.0331","-55.0205","-33.7666","-76.2744","-0.3799","-0.3796","-0.3803"],["'MLO'","2025","30","'month'","'BsG_S2S20'","'neph'","'log'","'LMS'","1.7301","90","-0.0053","8.2207e-04","-0.0114","-14.7924","2.3074","-31.8922","-0.1462","-0.1458","-0.1467"],["'MLO'","2025","10","'month'","'BbsG_S2S20'","'neph'","'log'","'LMS'","0.1620","0","0.0026","0.0343","-0.0291","0.1419","1.8936","-1.6098","0.0260","0.0269","0.0251"],["'MLO'","2025","20","'month'","'BbsG_S2S20'","'neph'","'log'","'LMS'","1.0755","0","-0.0056","0.0048","-0.0161","-0.3119","0.2681","-0.8919","-0.1063","-0.1058","-0.1068"],["'MLO'","2022","10","'month'","'Ba_ae_psap_clap'","'abs'","'log'","'LMS'","3.8480","95","-0.0820","-0.0394","-0.1247","-3.1666","-1.5208","-4.8125","-0.5597","-0.5592","-0.5602"],["'MLO'","2022","20","'month'","'Ba_ae_psap_clap'","'abs'","'log'","'LMS'","0.1892","0","-0.0020","0.0192","-0.0232","-0.0794","0.7600","-0.9189","-0.0393","-0.0382","-0.0404"],["'MLO'","2022","30","'month'","'Ba_ae_psap_clap'","'abs'","'log'","'LMS'","4.0497","95","0.0256","0.0383","0.0130","0.9632","1.4389","0.4875","1.1564","1.1587","1.1542"],["'MLO'","2025","10","'month'","'expS_bg'","'neph'","'log'","'LMS'","0.2871","0","-0.0021","0.0123","-0.0164","-0.8308","4.9568","-6.6183","-0.0204","-0.0200","-0.0207"],["'MLO'","2025","20","'month'","'expS_bg'","'neph'","'log'","'LMS'","0.8500","0","-0.0022","0.0029","-0.0072","-0.8707","1.1780","-2.9194","-0.0422","-0.0419","-0.0425"],["'MLO'","2025","30","'month'","'expS_bg'","'neph'","'log'","'LMS'","1.5589","0","-0.0029","8.2071e-04","-0.0066","-1.0720","0.3033","-2.4473","-0.0833","-0.0831","-0.0836"],["'MLO'","2025","10","'month'","'BbsFg'","'neph'","'log'","'LMS'","3.7701","95","0.0499","0.0763","0.0234","2.7091","4.1463","1.2720","0.6465","0.6477","0.6453"],["'MLO'","2025","20","'month'","'BbsFg'","'neph'","'log'","'LMS'","7.8729","95","0.0244","0.0306","0.0182","1.2507","1.5684","0.9330","0.6289","0.6294","0.6283"]]}}
%---
%[output:57b48bc0]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"19×19 table","name":"MLO_resultLMSlin","rows":19,"type":"table","value":[["'MLO'","2025","10","'month'","'Ba3_A82_ae33'","'abs'","'lin'","'LMS'","3.4046","95","-0.0351","-0.0145","-0.0557","-14.2518","-5.8797","-22.6240","-1.3511","-1.3506","-1.3517"],["'MLO'","2025","10","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","3.3361","95","-0.0672","-0.0269","-0.1076","-7.6846","-3.0776","-12.2916","-1.6724","-1.6713","-1.6735"],["'MLO'","2025","20","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","3.2030","95","-0.0386","-0.0145","-0.0627","-4.0329","-1.5147","-6.5511","-1.7723","-1.7710","-1.7736"],["'MLO'","2025","30","'month'","'BsG_S2S20'","'neph'","'lin'","'LMS'","2.4870","95","-0.0170","-0.0033","-0.0306","-1.7593","-0.3445","-3.1740","-1.5093","-1.5082","-1.5104"],["'MLO'","2025","10","'month'","'BbsG_S2S20'","'neph'","'lin'","'LMS'","0.6096","0","-0.0018","0.0041","-0.0077","-1.1014","2.5122","-4.7150","-1.0180","-1.0179","-1.0182"],["'MLO'","2025","20","'month'","'BbsG_S2S20'","'neph'","'lin'","'LMS'","1.2605","0","-0.0018","0.0011","-0.0046","-1.0851","0.6366","-2.8067","-1.0358","-1.0357","-1.0360"],["'MLO'","2022","10","'month'","'Ba_ae_psap_clap'","'abs'","'lin'","'LMS'","3.7449","95","-0.0090","-0.0042","-0.0138","-12.0274","-5.6040","-18.4507","-1.0902","-1.0901","-1.0903"],["'MLO'","2022","20","'month'","'Ba_ae_psap_clap'","'abs'","'lin'","'LMS'","0.0689","0","-1.2006e-04","0.0034","-0.0036","-0.1501","4.2058","-4.5059","-1.0024","-1.0022","-1.0026"],["'MLO'","2022","30","'month'","'Ba_ae_psap_clap'","'abs'","'lin'","'LMS'","3.1722","95","0.0028","0.0046","0.0011","4.0663","6.6300","1.5026","-0.9146","-0.9145","-0.9148"],["'MLO'","2025","10","'month'","'expS_bg'","'neph'","'lin'","'LMS'","0.4663","0","-0.0042","0.0139","-0.0223","-0.3296","1.0841","-1.7434","-1.0422","-1.0417","-1.0427"],["'MLO'","2025","20","'month'","'expS_bg'","'neph'","'lin'","'LMS'","1.0456","0","-0.0034","0.0031","-0.0099","-0.2657","0.2426","-0.7740","-1.0681","-1.0677","-1.0684"],["'MLO'","2025","30","'month'","'expS_bg'","'neph'","'lin'","'LMS'","1.4352","0","-0.0036","0.0014","-0.0086","-0.2738","0.1078","-0.6554","-1.1076","-1.1072","-1.1080"],["'MLO'","2025","10","'month'","'BbsFg'","'neph'","'lin'","'LMS'","3.5173","95","0.0098","0.0154","0.0042","6.1681","9.6754","2.6609","-0.9021","-0.9020","-0.9023"],["'MLO'","2025","20","'month'","'BbsFg'","'neph'","'lin'","'LMS'","6.4245","95","0.0041","0.0054","0.0029","2.9154","3.8230","2.0078","-0.9171","-0.9170","-0.9172"]]}}
%---
%[output:6b7c2cfc]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'MLO_result_GSMd'."}}
%---
