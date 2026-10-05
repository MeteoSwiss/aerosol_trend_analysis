APP_st.name='APP';
APP_st.lat=36.2130;
APP_st.lon=-81.6920;
APP_st.alt=1076;
APP_st.env='con';
APP_st.footp='RB';
%%

APP_rd=read_betsy_2026('app_tsi_psap_clap_eco',APP_st.name); % not full %[output:6af139ed]
names_rd=fieldnames(APP_rd);

% % APP_rd_nc=read_ebas_merged_STN('C:\github_trend\raw_data\merged_Appalachian State University.nc', APP_st.name);% not full
% % names_rd_nc=fieldnames(APP_rd_nc);
% % plotFigControl(APP_rd_nc,APP_st.name);
%%
datevec(APP_rd.Time(1)) %data since 2009 %[output:514e2ecf]
datevec(APP_rd.Time(end)) %2025 %[output:5c62887c]
prctile(APP_rd.U0_S11_TSI,[5 50 95]) % < 38% ok no need to compute dry trend %[output:9b43f96e]
prctile(APP_rd.U1_S11_TSI,[5 50 95]) %[output:342d5901]
plotFigControl(APP_rd,APP_st.name); %[output:807fa570] %[output:4564360e] %[output:9dc72aa1] %[output:04266fdc] %[output:448bf181] %[output:9e3bfa85] %[output:139a5ea0] %[output:3075d01f] %[output:80619bc2]
% 2 size cut, ok
% sc: ok, max in summer, no negatives
% abs: 2 size cut, max in summer, no negatives
%very high seasonal cycle of abs ratio and abs exponent
%N:july/August 2011: too low: problem and rupture, often lower in summer

%problems:
%ratio B/R and G/R with higher max and min values since mid 2016: change PSAP --> CLAP
%should be due to the red wavelength change (660 to 653)

% ratio BbsG0/BsG1 and BbsG1/BsG1 and BbsG0/BsG0 rupture in the second part of september 2015?
% %and BbsG1/BsG0 very low in summer 2015
% SSA also show a rupture in september 2015 as well and the backscattering fraction and the exp BG
%%
% change from TSI to ecotech
% common period: 30.9.2023 - 2.1.2024
Pt=timerange('2023-09-30','2024-01-02');
ratio_med_h=nanmedian(APP_rd.BsG0_S11_TSI(Pt)./APP_rd.BsG0_S13_eco(Pt)) % 1.0772 %[output:054c6af0]
ratio_STN_h=nanstd(APP_rd.BsG0_S11_TSI(Pt)./APP_rd.BsG0_S13_eco(Pt)) % 0.2640 %[output:24547367]
ratio_mean_h=nanmean(APP_rd.BsG0_S11_TSI(Pt)./APP_rd.BsG0_S13_eco(Pt)) % 1.0958 %[output:25e48712]
ratio_mean=nanmean(APP_rd.BsG0_S11_TSI(Pt))/nanmean(APP_rd.BsG0_S13_eco(Pt)) % 0.8914 %[output:5075afc3]
ratio_median=nanmedian(APP_rd.BsG0_S11_TSI(Pt))/nanmedian(APP_rd.BsG0_S13_eco(Pt)) % 0.8763 %[output:9ef7d6fd]
figure; %[output:82396b59]
plot(APP_rd.BsG0_S11_TSI(Pt),APP_rd.BsG0_S13_eco(Pt),'.'); %[output:82396b59]

% fit between both data: slope 1.028x-0.567 , 1/1.028=0.9728
fit_test=robustfit(APP_rd.BsG0_S11_TSI(Pt),APP_rd.BsG0_S13_eco(Pt)) % y coor=-1.6968, slope=1.0405, 1/slope=0.9611 %[output:5fdc60cc]

cte_tsi_eco=0.95;
%%
% homogeneisation scat
P1=timerange('2009-01-01','2023-10-01');
P2=timerange('2023-10-01','2026-01-01');

% B
APP_rd.BsB0_homo1=NaN(size(APP_rd.BsB0_S11_TSI));
APP_rd.BsB0_homo1(P1)=APP_rd.BsB0_S11_TSI(P1).*cte_tsi_eco;
APP_rd.BsB0_homo1(P2)=APP_rd.BsB0_S13_eco(P2);

% G
APP_rd.BsG0_homo1=NaN(size(APP_rd.BsG0_S11_TSI));
APP_rd.BsG0_homo1(P1)=APP_rd.BsG0_S11_TSI(P1).*cte_tsi_eco;
APP_rd.BsG0_homo1(P2)=APP_rd.BsG0_S13_eco(P2);

% R
APP_rd.BsR0_homo1=NaN(size(APP_rd.BsR0_S11_TSI));
APP_rd.BsR0_homo1(P1)=APP_rd.BsR0_S11_TSI(P1).*cte_tsi_eco;
APP_rd.BsR0_homo1(P2)=APP_rd.BsR0_S13_eco(P2);

% Bbs
% B
APP_rd.BbsB0_homo1=NaN(size(APP_rd.BbsB0_S11_TSI));
APP_rd.BbsB0_homo1(P1)=APP_rd.BbsB0_S11_TSI(P1).*cte_tsi_eco;
APP_rd.BbsB0_homo1(P2)=APP_rd.BbsB0_S13_eco(P2);

% G
APP_rd.BbsG0_homo1=NaN(size(APP_rd.BbsG0_S11_TSI));
APP_rd.BbsG0_homo1(P1)=APP_rd.BbsG0_S11_TSI(P1).*cte_tsi_eco;
APP_rd.BbsG0_homo1(P2)=APP_rd.BbsG0_S13_eco(P2);

% R
APP_rd.BbsR0_homo1=NaN(size(APP_rd.BbsR0_S11_TSI));
APP_rd.BbsR0_homo1(P1)=APP_rd.BbsR0_S11_TSI(P1).*cte_tsi_eco;
APP_rd.BbsR0_homo1(P2)=APP_rd.BbsR0_S13_eco(P2);
%%
% BP detection for scat coef and backscat

names_sc={'BsB0_homo1', 'BsG0_homo1', 'BsR0_homo1'};
break_APP_scat_BP=change_point_analysis_Def(APP_rd,names_sc,0.05,'APP','scattering'); %[output:0c16ee21] %[output:23593bf4] %[output:50842089] %[output:9777a6f5]
T_APP_scat_BP=make_table_breakpoints_def(break_APP_scat_BP); %[output:554728ec] %[output:402ae825] %[output:9e326ede] %[output:5810a729]

names_bsc={'BbsB0_homo1', 'BbsG0_homo1', 'BbsR0_homo1'};
break_APP_bscat_BP=change_point_analysis_Def(APP_rd,names_bsc,0.05,'APP','Backscattering'); %[output:9ec141c0] %[output:17dd0a36] %[output:7233dc42] %[output:9c097411]
T_APP_bscat_BP=make_table_breakpoints_def(break_APP_bscat_BP); %[output:15e5ec7a]
%%
%Abs homogeneisation between PSAP and CLAP
% common period= 4.2.2012 - 15.1.2016
Pt=timerange('2012-02-01','2016-03-01');
ratio_med_h=nanmedian(APP_rd.BaG0_A11_psap(Pt)./APP_rd.BaG0_A12_clap(Pt)) % 1.0421 %[output:3af32782]
ratio_STN_h=nanstd(APP_rd.BaG0_A11_psap(Pt)./APP_rd.BaG0_A12_clap(Pt)) % 0.2887 %[output:1bf2e1f7]
ratio_mean_h=nanmean(APP_rd.BaG0_A11_psap(Pt)./APP_rd.BaG0_A12_clap(Pt)) % 1.0902 %[output:4534507f]
ratio_mean=nanmean(APP_rd.BaG0_A11_psap(Pt))/nanmean(APP_rd.BaG0_A12_clap(Pt)) % 1.0597 %[output:783a1690]
ratio_median=nanmedian(APP_rd.BaG0_A11_psap(Pt))/nanmedian(APP_rd.BaG0_A12_clap(Pt)) % 1.0893 %[output:6101957b]
figure; %[output:0788e32b]
plot(APP_rd.BaG0_A11_psap(Pt),APP_rd.BaG0_A12_clap(Pt),'.'); %[output:0788e32b]

% fit between both data: slope 0.8657x+0.1804, 1/0.8657=1.1551
fit_test=robustfit(APP_rd.BaG0_A11_psap(Pt),APP_rd.BaG0_A12_clap(Pt)) % y coor=0.1130, slope=0.9029, 1/slope=1.1075 %[output:54d2849a]
%since the fit pass over (0,0), forcing him through 0,0 will lead to lower 1/slope

cte_psap_clap=0.95;
%%
% homogeneisation abs
%P1=timerange('2009-01-01','2011-01-01'); data remove due to BP at end or
%PSAP-1w
P1=timerange('2011-01-01','2013-01-01');
P2=timerange('2013-01-01','2026-01-01');


% B
APP_rd.BaB0_S11S12=NaN(size(APP_rd.BaB0_A11_psap));
APP_rd.BaB0_S11S12(P1)=APP_rd.BaB0_A11_psap(P1).*cte_psap_clap;
APP_rd.BaB0_S11S12(P2)=APP_rd.BaB0_A12_clap(P2);

% G
APP_rd.BaG0_S11S12=NaN(size(APP_rd.BaG0_A11_psap));
APP_rd.BaG0_S11S12(P1)=APP_rd.BaG0_A11_psap(P1).*cte_psap_clap;
APP_rd.BaG0_S11S12(P2)=APP_rd.BaG0_A12_clap(P2);

% R
APP_rd.BaR0_S11S12=NaN(size(APP_rd.BaR0_A11_psap));
APP_rd.BaR0_S11S12(P1)=APP_rd.BaR0_A11_psap(P1).*cte_psap_clap;
APP_rd.BaR0_S11S12(P2)=APP_rd.BaR0_A12_clap(P2);
%%

names_abs={'BaB0_S11S12', 'BaG0_S11S12', 'BaR0_S11S12'};
break_APP_abs_BP=change_point_analysis_Def(APP_rd,names_abs,0.05,'APP','Absorption'); %[output:06f55808] %[output:94d2eb33] %[output:627fa8c7] %[output:1c2a28bc]
T_APP_abs_BP=make_table_breakpoints_def(break_APP_abs_BP); %[output:8eb7bf5c] %[output:298d53cc] %[output:2b5eedfc]

%%
lambdaSCtsi=[450;550;700];
lambdaSCeco=[450;528;652];
lambdaAE=[467;530;660];
names_sc1={'BsB0_S11_TSI', 'BsG0_S11_TSI', 'BsR0_S11_TSI'};
APP_expS1=compute_exp_D(APP_rd,names_sc1, lambdaSCtsi);
names_sc2={'BsB0_S13_eco', 'BsG0_S13_eco', 'BsR0_S13_eco'};
APP_expS2=compute_exp_D(APP_rd,names_sc2, lambdaSCeco);
P1=timerange('2009-01-01','2023-10-01');
P2=timerange('2023-10-01','2026-01-01');
APP_expS1.expS_bg(P2)=APP_expS2.expS_bg(P2);
APP_expS1.expS_br(P2)=APP_expS2.expS_br(P2);
APP_expS1.expS_gr(P2)=APP_expS2.expS_gr(P2);

APP_expA=compute_exp_D(APP_rd,names_abs, lambdaAE);
APP_cal=outerjoin(APP_expS1,APP_expA);
%names_SSA={'BsB0_homo1', 'BsG0_homo1', 'BsR0_homo1','BaB0_S11S12', 'BaB0_S11S12', 'BaB0_S11S12' };
APP_SSA=compute_SSA_D(APP_rd, names_sc, names_abs);
APP_cal=outerjoin(APP_cal,APP_SSA);

APP_cal.BbsFG0_homo=APP_rd.BbsG0_homo1./APP_rd.BsG0_homo1;
APP_cal.BbsFB0_homo=APP_rd.BbsB0_homo1./APP_rd.BsB0_homo1;
APP_cal.BbsFR0_homo=APP_rd.BbsR0_homo1./APP_rd.BsR0_homo1;

names_cal=fieldnames(APP_cal);
plotFigControl_cal(APP_cal,APP_st.name); %[output:2fe376f1] %[output:8ddf09ed] %[output:747a1fce]
%%
% BP detection for expS, Backscat fraction, expA and SSA

names_exps={'expS_bg', 'expS_br', 'expS_gr'};
break_APP_exps_BP=change_point_analysis_Def(APP_cal,names_exps,0.05,'APP','scattering exp'); %[output:77017bdb] %[output:5d36622a] %[output:69cbfed6] %[output:2dac2ad6]
T_APP_exps_BP=make_table_breakpoints_def(break_APP_exps_BP); %[output:7674b0e4] %[output:713075df] %[output:599eed85] %[output:97a709a7] %[output:9f15906d] %[output:3b28743a] %[output:0984f883] %[output:2f6cc226] %[output:9e6ab244] %[output:5bfaff3a] %[output:5fd7fc21]

names_exps={'expA_bg', 'expA_br', 'expA_gr'};
break_APP_expa_BP=change_point_analysis_Def(APP_cal,names_exps,0.05,'APP','Absorption exp'); %[output:269b3aec] %[output:8c096be6] %[output:2d823996] %[output:5389054d]
T_APP_expa_BP=make_table_breakpoints_def(break_APP_expa_BP); %[output:7d33eba5] %[output:2c2d1b64] %[output:4ccbbbdf] %[output:71b1e216] %[output:0b10475a] %[output:81f74aba] %[output:7ab9e0b1]

names_BbF={'BbsFG0_homo'};
break_APP_BbF_BP=change_point_analysis_Def(APP_cal,names_BbF,0.05,'APP','backscat fraction'); %[output:66b3b500] %[output:67188cfd]
T_APP_BbF_BP=make_table_breakpoints_def(break_APP_BbF_BP); %[output:2412e4f9] %[output:20242383] %[output:17ab63e2] %[output:4e4b1496]

names_SSA={'SSAB','SSAG','SSAR'};
break_APP_SSA_BP=change_point_analysis_Def(APP_cal,names_SSA,0.05,'APP','SSA'); %[output:8fa6af9b] %[output:476babb8] %[output:5cc08166] %[output:73358e43]
T_APP_SSA_BP=make_table_breakpoints_def(break_APP_SSA_BP); %[output:66e85348] %[output:4c901919]
%%
%Questions:

%%

APP_tr=outerjoin(APP_rd,APP_cal);
APP_tr.y=year(APP_tr.Time);
% begin at the beginning of a year: 2010 for neph and 2011 for abs
%end: 2025
P=timerange('2010-01-01','2026-01-01');
APP_tr=APP_tr(P,:);
% large RH change in october 2013: trend on U to do, but RH<50%, no need to trend on dry
names=fieldnames(APP_tr);
c= startsWith(names,["DOY";"U";"P0_";"N";"Uu"]) | endsWith(names,["dry";"_S11_TSI";"S13_eco";"A11_psap";"A12_clap"]) | contains(names,'Q');
N=names(c);
for i=1:length(N)
    APP_tr.(N{i})=[];
end

%scat: 2-17.1.2015: scat + back red to invalidate
P2=timerange('2015-01-02','2015-01-17');
APP_tr.BsR0_homo1(P2)=NaN;
%APP_tr.BsR1_S(P2)=NaN;
APP_tr.BbsR0_homo1(P2)=NaN;
APP_tr.expS_br(P2)=NaN;
APP_tr.expS_gr(P2)=NaN;
APP_tr.SSAR(P2)=NaN;

%APP_tr.BbsR1_S(P2)=NaN;
%scat: 15-21.9.2015: scat + backscat to invalidate, all wavelengths
P3=timerange('2015-09-15','2015-09-21');
APP_tr.BsB0_homo1(P3)=NaN;
%APP_tr.BsB1_homo1(P3)=NaN;
APP_tr.BbsB0_homo1(P3)=NaN;
%APP_tr.BbsB1_homo1(P3)=NaN;
APP_tr.BsG0_homo1(P3)=NaN;
%APP_tr.BsG1_homo1(P3)=NaN;
APP_tr.BbsG0_homo1(P3)=NaN;
%APP_tr.BbsG1_homo1(P3)=NaN;
APP_tr.BsR0_homo1(P3)=NaN;
%APP_tr.BsR1_homo1(P3)=NaN;
APP_tr.BbsR0_homo1(P3)=NaN;
%APP_tr.BbsR1_homo1(P3)=NaN;
APP_tr.expS_bg(P3)=NaN;
APP_tr.expS_br(P3)=NaN;
APP_tr.expS_gr(P3)=NaN;
APP_tr.SSAB(P3)=NaN;
APP_tr.SSAG(P3)=NaN;
APP_tr.SSAR(P3)=NaN;

% BbsF not ok 
P4=timerange('2023-11-11','2023-12-31');
APP_tr.BbsFB0_homo(P4)=NaN;
APP_tr.BbsFG0_homo(P4)=NaN;
APP_tr.BbsFR0_homo(P4)=NaN;
APP_tr.expS_bg(P4)=NaN;
APP_tr.expS_br(P4)=NaN;
APP_tr.expS_gr(P4)=NaN;

% PROBABLY RESTRIC BbsF AND Exps TO BEFORE 2022 AND ExpA FROM 2013

% trend only on scat G
APP_tr.BsB0_homo1=[];
%APP_tr.BsB1_S=[];
APP_tr.BsR0_homo1=[];
%APP_tr.BsR1_S=[];
APP_tr.BbsB0_homo1=[];
%APP_tr.BbsB1_S=[];
APP_tr.BbsR0_homo1=[];
%APP_tr.BbsR1_S=[];
APP_tr.BbsFB0_homo=[];
APP_tr.BbsFR0_homo=[];


% trend on BG expS
APP_tr.expS_br=[];
APP_tr.expS_gr=[];
%APP_tr.expS_br1=[];
%APP_tr.expS_gr1=[];
%abs: 7.3.2016 PSAP to CLAP: ok for coef, BUT break point for ExpA at least with red: use BG et caution if positive trend
% trend only on abs G 
APP_tr.BaB0_S11S12=[];
%APP_tr.BaB1_A=[];
APP_tr.BaR0_S11S12=[];
%APP_tr.BaR1_A=[];

%trend only on ExpS bg
APP_tr.expA_br=[];
APP_tr.expA_gr=[];
%APP_tr.expA_br1=[];
%APP_tr.expA_gr1=[];

%%
% restricted trend with 
APP_tr_V2=APP_tr;
P5=timerange('2023-01-01','2025-01-01');
%APP_tr_V2.BbsFB0_homo(P5)=NaN;
APP_tr_V2.BbsFG0_homo(P5)=NaN;
%APP_tr_V2.BbsFR0_homo(P5)=NaN;
APP_tr_V2.expS_bg(P5)=NaN;
%APP_tr_V2.expS_br(P5)=NaN;
%APP_tr_V2.expS_gr(P5)=NaN;

P6=timerange('2009-01-01','2013-01-01');

APP_tr_V2.expA_bg(P6)=NaN;

%%
% trend computing

[APP_result_MK,APP_result_LMSlog,APP_result_LMSlin]=all_trend_STN(APP_tr,APP_st); %[output:39c29f9f] %[output:51561f91] %[output:55b61ecf] %[output:153dd0c9] %[output:8cd5c8a0] %[output:5b5d7fac] %[output:0c72d7e0] %[output:0e92f4b7] %[output:13bf07f7] %[output:3b526546] %[output:09d2e8cf] %[output:61d5c65d] %[output:854e323f] %[output:93c842ab] %[output:4292f193] %[output:7fd07625] %[output:6311627c] %[output:4333cf51] %[output:412c8e4c] %[output:2fdfec66] %[output:107dc2fd] %[output:8c4db4ca] %[output:2079213c] %[output:6eb425fd] %[output:8673b863] %[output:019c3e9c] %[output:1b97fc30] %[output:05a9a300] %[output:8ad946a0] %[output:30f38b1f] %[output:3d6d6bc5] %[output:25020dcc] %[output:7a90955c] %[output:2b6c7968] %[output:6afba108] %[output:27211ad6] %[output:5643eb4a] %[output:84dfb973] %[output:5511e390] %[output:9f4fd056] %[output:33daa2f1] %[output:5c87ca2a] %[output:70356c85] %[output:37765ff1] %[output:83f69a0c] %[output:30b6f92a] %[output:645802aa] %[output:3023b4ce] %[output:332f588f] %[output:18ad7d9b] %[output:37f52a9d] %[output:62c2f4cd] %[output:3a51b7a4] %[output:847d426f] %[output:53643fab] %[output:514002c3] %[output:69ca112e] %[output:84fa34af] %[output:4c9eca56] %[output:51ebd9dd] %[output:35697a8f] %[output:304accda] %[output:42676b03] %[output:8cab75b1] %[output:823c3cef] %[output:6aa20f99] %[output:7f7ad05e] %[output:95dd3faf] %[output:51e6399a] %[output:6984b449] %[output:4a0f24bc] %[output:90f9bde8] %[output:7988a91f] %[output:34d47b0f] %[output:9fc1ddc7] %[output:29d9ebf3] %[output:51946b13] %[output:73a7e3d4] %[output:7228ae3b] %[output:10ec561a] %[output:7d06daef] %[output:2ed4dcee] %[output:8fce9b6b] %[output:38a36ee6] %[output:05038c20] %[output:7be4cc95] %[output:96a47bbf] %[output:78564df0] %[output:824d1e64] %[output:40846abf] %[output:88e4cfa1] %[output:1c653a11] %[output:3bef9a26] %[output:7c70dc28] %[output:53c32b96] %[output:54b107ba] %[output:198be025] %[output:33e1fcda] %[output:888dab7f] %[output:49db6cb5] %[output:649ffd32] %[output:6ce2d67e] %[output:14a83052]

writetable(APP_result_MK,'APP_res_MK.txt'); %, 'delimiter',',' )
writetable(APP_result_LMSlog,'APP_res_LMSlog.txt'); 
writetable(APP_result_LMSlin,'APP_res_LMSlin.txt'); 
plot_10y_in_two(APP_result_MK, APP_st,'y'); %[output:929ca2c8] %[output:5e50f5b2]
%%
% with restricted expS, expA and Bbsf
[APP_result_MK2,APP_result_LMSlog2,APP_result_LMSlin2]=all_trend_STN(APP_tr_V2,APP_st); %[output:83b20349] %[output:20702d7a] %[output:5d1cd89b] %[output:4bdda39b] %[output:4375a328] %[output:9fc16f2b] %[output:3559a0a8] %[output:13a73f18] %[output:0a26820c] %[output:057e9057] %[output:4713ce00] %[output:5739c69a] %[output:64b973cb] %[output:975d8141] %[output:22c5c274] %[output:6179cbfe] %[output:72c89fc3] %[output:974c0935] %[output:715e5e51] %[output:4e757775] %[output:4f78cb31] %[output:0217a16d] %[output:8209d2f6] %[output:489bd53b] %[output:03700d07] %[output:3eec10cd] %[output:79ce1524] %[output:639e87fc] %[output:1bf85273] %[output:539c218e] %[output:6d130db8] %[output:2f74aa87] %[output:78edf7d2] %[output:6438a825] %[output:41f5bdd7] %[output:7724eb22] %[output:85b5cbf9] %[output:143e9611] %[output:548ffb1f] %[output:5a251bca] %[output:5fb69ed2] %[output:2d32bbd9] %[output:6477dd42] %[output:0c081830] %[output:00564cc3] %[output:2483e4a0] %[output:964403f4] %[output:218196bb] %[output:836d1211] %[output:74c64c8e] %[output:6350f5ab] %[output:87e6bf5b] %[output:80126126] %[output:70179f68] %[output:359fd545] %[output:3832f91a] %[output:67ec98b2] %[output:1c680544] %[output:662cec20] %[output:48c661c3] %[output:734adad5] %[output:61faa441] %[output:5ce14ffa] %[output:58c7b668] %[output:356bc868] %[output:6108f4fb] %[output:7cd0802e] %[output:5e13cf95] %[output:49f09336] %[output:9ec7c045] %[output:7bbc295d] %[output:06635f4c] %[output:1c5545a4] %[output:82bf21cc] %[output:7c474475] %[output:412e81f2] %[output:22b26ff3] %[output:1489397c] %[output:4e7d2b49] %[output:1cd9348e] %[output:22d9ec64] %[output:07c308a6] %[output:78a4c413] %[output:3559a220] %[output:18d2a38f] %[output:2049bf38] %[output:05fc8b66] %[output:09c01672] %[output:251df1f1] %[output:5eeaac3f] %[output:15e66d17] %[output:19f9922f] %[output:558013dd] %[output:39a87371] %[output:6d120068] %[output:2bc3125a] %[output:6b1231d8] %[output:156164f3] %[output:389777bd] %[output:60cc417f] %[output:1db86f63] %[output:09b60491] %[output:409e85e2] %[output:4c30d1b1] %[output:6e70fdf7] %[output:44910ae2]

writetable(APP_result_MK2,'APP_res_MK2.txt'); %, 'delimiter',',' )
writetable(APP_result_LMSlog2,'APP_res_LMSlog2.txt'); 
writetable(APP_result_LMSlin2,'APP_res_LMSlin2.txt'); 
plot_10y_in_two(APP_result_MK2, APP_st,'y'); %[output:2120334f] %[output:41032e27]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":12}
%---
%[output:6af139ed]
%   data: {"dataType":"warning","outputData":{"text":"Warning: there is a lot of days with less than 50% data coverage"}}
%---
%[output:514e2ecf]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"ans","rows":1,"type":"double","value":[["2009","1","1","0","0","0"]]}}
%---
%[output:5c62887c]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"ans","rows":1,"type":"double","value":[["2025","12","31","0","0","0"]]}}
%---
%[output:9b43f96e]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":1,"type":"double","value":[["-3.1500","15.8000","38.9000"]]}}
%---
%[output:342d5901]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":1,"type":"double","value":[["-3.2000","15.5500","38.5950"]]}}
%---
%[output:807fa570]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAABoNJREFUeF7tnHtIVFkcx38jarkqbIqsjwRH0BRyk3XBIFQQ8Q+hJSEqBPNVGWSFEqiQYKblkxAVDK1MQiaI9g+pXSEQ13wlguKaC0K4iuKiKAuKtIgu38Oe2ev13pk717k6ztwD\/uHMeX7O7\/5e59wx7Ozs7JBeDpSAQYd+oLzZYDr0g2d+NKBvbW1RTU0NPXv2jCHy9\/enFy9e0OnTpyWRoe7Tp08lv3N3d6dTp07RhQsX6MqVK+Tt7c3qra6uUl5eHk1MTEi28\/X1pdjYWLp27RqdO3eO3NzcVG\/XkZD0+fl5ysnJoS9fvpgXWlhYSAUFBWQwGPYs3hJ0YeW0tDR6\/PgxAag16LwdNu3+\/fuUmZkpObaSnTgS0Lu7u+nu3bu71hMTE8OkOTAwUDV0NCwvL6erV68qho42AQEB1NHRQdHR0UoY76nj8NC\/fv1KZWVl9ObNGwoPD6ezZ89SV1cXW0h7ezslJydbhI66aMPL2toavXr1ipqamghqKz09nR49ekQbGxtm9XLmzBmmyvz8\/Fiz7e1tmp2dpaqqKurt7WWfNTQ0sLZqisNDn56epuzsbFpeXqZLly4xqbxx4wYtLi7SxYsX6eHDh3Ts2LFdaxeqFzF0VERfubm5NDU1RefPn6fq6mra3NyUhc477+\/vp6ysLPYv1NLly5fVMHd8Q9rZ2clUAAqkMyUlxSz5co+5Jei6pFuRk\/X1dSoqKqIPHz4w1QKPJTQ0lN69e0e3b99mrWHUILXCotSQok1jYyOTdqWGFG2MRiNTP2FhYc4n6WNjY0y1QN\/CW4Buh\/ewtLRE+fn5NDk5yfR1S0sLnThxwgxAKXSoK2yaj4+PYugeHh5MtUCfS3lOSnbBYXU6shPNzc305MmTPUZT6LdjE9ra2igpKUkRdO5vYxMTEhLM9sCapIeEhLD60OmRkZGqgWOSDgt9ZWWFGczx8XGrwoOgpri4mD0FKNYMqVSHQuhi78XqBGys4LDQ+\/r66Pr168yts1bEOlaHbo2YxPfisF9JFzzI0SVdCS2JOsKwPyoqinkKQUFBe2oKDW1qairV19czo6hLugrwwrBf6LWIu4LPfevWLRoeHmaJK4TmcXFxOnRbmQvDfrSVC\/V5v8+fP6fKykr2LzeoCNF5llEqItUNqa274gT1HdZ7cQK2skvQFDoCnIGBAZYzGR0dZYcPCNmRtOKHB84MV25tmkL\/+PEjM3Tx8fHslObz58\/08uVLSkxMpAcPHjBPwxWLZtC5ZwHpRq6CA0bQc+fOHfYZTm5csWgGHW4c1Ehra+uugwaeOcQBQUVFBXl6erocd82gm0wmlrBCOjYiIsIMFnoeUg4dLzydcSXymkFHVDg4OMj8bBw2CAu+e\/\/+PQtmkDcRl78qBumfP\/8+9H34OfMHWvnOV3Iewcfd6KbxuKo5agodKkZKml+\/fs08GvFTgBWsdf5O83m\/qlqMPRv1p0ZQ273\/08VSfd8MO64KvMNBn8\/7hdY6p+zJT1VfbfcSqT810mLbnwI9qSL6G5v7PxToltTLUZJ0AAd4W4tm0PdjSKHT13+bt3Utdq8Pnf7H93uzmxjox2\/dVakWtNUMuu4yysuAZtB5cIQ8uDD61IMjIsPQ0JBm99M\/ffrEfHVcgcPB8cLCAvX09BDOIHHK7+XlZXeVYO8OT548SfizZzEYjUbNoNtzoofVF\/JGdXV1dgWvqaQfFih7jTsyMsIuIyk9CFE6rmY6XekEHLkenIGMjAwd+kFukg79IGn\/N5YO3RWgI1RHwXU1qYKkFbKIuNttq+uHtqWlpbu6xYVQ4Vj8mhsuiMrNQet9OFBJ51DEIPgi+WT4hXpboGMzxdlHDhj986yky0DH2wglJSWEyz4oUtCFt6dshT4zM8Ne2MJNLOErKRhL\/J1LQOfA5+bmWBRZW1tLwcHBux5tAMeGIA\/+9u1b9gqKLeqFg8WFfmuvjrgEdKF+5Bsghi6sgw2wFbrwSbL2lOjQJSyWGujoRu7yvTjq06ErhC732oncJXuubvDEoAjrOTV0ZBnFRm2\/6oUDg6sn7lvOzRN7TE4Nvbq6ekfsB2sF3ZJK4mNiU4TvdbqMn64VdEvBhnhMp5Z0qd970Qo6pJjrfbHhFAdNTg3dZDLtiH1mLaEDPJd4oX53qTSAlE7fb05DjSHd75hatNcs9yIl6ftdgA7dMkFNTo506IcAfb9PiqO010y96L9WJ7\/FWkH\/F34Ek6RdYqivAAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:4564360e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAACD1JREFUeF7tnGtsTVkUx9dtqmqoDMJUEb2ipR7FGEFQiUg\/EA0iiKSKelS8gnq\/6\/2OeKtniSAeCakZIcSUekUQLVUiJg0xUZqJek3EnfxXZt\/Z99x97rnn9t726Jyd+NB79j5nn99eZ+211l6Lw+VyuchulUrAYUOvVN78MBt65TO3HvSvX7\/S2rVrad++fYyjQYMGdODAAWrXrp0SD\/ru3r1beS08PJxatWpFAwcOpOHDh1Pt2rW537t37yg9PZ0ePHigHBcVFUUdO3aksWPHUo8ePSgsLCyoS2M5SS8pKaHRo0fT8+fP3S86ffp0mjx5MjkcDq+X9wVd7tyvXz9avXo1AagRdDEOi7Zw4UJKTU1VPjvQlbAc9HPnztG0adM83qd9+\/YszdHR0QFDx8ClS5fSyJEj\/YaOMQ0bNqSDBw9SQkJCoIy9xlkK+pcvX2jRokV08uRJatGiBXXr1o2OHj3Kk967dy\/16dPHJ3T0xRjRysrK6MiRI7R161aC2ho0aBCtWrWKPnz44FYvHTp0YFVWv359Hvbt2zd68eIFrVy5kq5cucK\/bdy4kccGq1kK+uPHj2nUqFH05s0bGjp0KEvl+PHj6dWrVzRkyBBavnw51axZ0+PdZfWihY6OuNeYMWOosLCQBgwYQGvWrKFPnz7pQhc3z8vLo7S0NP4TamnYsGHBYm6tjTQnJ4dVABqks2\/fvm7J1\/vMfUG3Jd1ATsrLy2nGjBl06dIlVi2wWJo1a0a5ubk0ZcoUHo1NDVIrN383UozZsmULS7u\/GynGOJ1OVj+xsbHVT9Lv3r3LqgX6FtYCdDush9evX9OECRPo4cOHrK+3b99O9erVcwPwFzrUFRatTp06fkOvUaMGqxboc5XlFOgqWEKnIxKxbds22rx5s9emKdvtWITs7Gzq3bu3X9CFvY1F7NWrl3s\/MJL0Jk2acH\/o9Pj4+KACx8QtAb20tJQ3zPv37xsKD5yaOXPm8FeAZrSRqm4oQ9daL4YTCEIHS0C\/evUqjRs3js06o6bVsTZ0I2KK61q3359bCCfHlnR\/aCn6yG5\/69at2VJo3LixV095o01OTqYNGzbwpmhLegDgZbdftlq0t4LNPWnSJLp58yYHruCad+7c2YZulrns9mOsnqsv7rt\/\/35asWIF\/yk2VLjoIsqo8kjtjdTsqlTT\/pawXqopW93XMg0djsz169c5NnLnzh0+ZIBrjuCUOCT4v0E0+76moV+7do03tK5du\/JpzKNHj+jQoUOUlJREy5YtY4vCbr4JmIIuLAhIN2ISAjCcm6lTp\/JvOKGxWxChw1yDGtm1a5fHgYKIEOIgICsriyIiImzuPgiYkvRjx45xYAph17i4OPdtoech5dDx8imMTV5NwBR0eH\/5+flsT+NQQW64dv78eXZaEB\/Rtj+z8unvP\/76btbhTOrPVPpTlO58YyLDKMMZGdD7mIYOFaOS5uPHj7NFo\/0KMKuynAIqSf8toAlWxaC85DjKzvwvfKw3h4zYyIDAVwr0kvRfqSynsCr4BfTM7MwkykuONxybEh1BWQk\/GPbTdggadF\/qpbpKOoADvNlmCnpFNlLo9PLfS8zOr8r6Q6cXJXpHO8WEfvkxPCDVgvGmoNsmY3BkwBR04Rwh3i17n7ZzZG4xHDdu3DCVn3779m221ZHqhgPily9f0oULFwhnjTjNr1WrlrkZWLx306ZNCf+C2RxOp9MU9GA+\/Hu4F2JM69evDyp405L+PYAK1hxv3brFCUr+Ho74+1xTOt3fm1aXfjAcRowYYUOvzAW1oVcm7X+fJUOPiYnhDARYau\/fvydkLuDwBrmR2kxieaqXL1+mjIwMQnKsSOO21YuPxRTQkQMJ3Y7EVqRMN2rUiIN7Z8+eZfCoFFGBR547MteePXvmoaKU0J8+fUqzZ8+mdevWeYRwxfxEWhrS2+QkfH+EUZVHCCnSBsp8BdD8eU4w+gjoMI1xZoCAnqgGQTgb0gvpF+kg8jOR6YCzBSwMkmLlzdgLuoCCZHpVxBAJ9XPnziXkq5jd1bGYqCdC6rOcZA\/A8+bN80i+txL0xMRE6t69O+dQyq24uJjfB6dm2qKB06dPc4YxDn3gSOpCFyuLG6ukT0BDZQSaWegIimEsqiG0TpT2mpWg673nvXv3GCrUjwwdizFx4kSaOXMml9VoLSC3pAvgOAFq3rw5ZWZmeki6AI6NY\/DgwbzC+LTMqBeA1YvHa9WB1aEjB3PTpk106tQpj0IwqKElS5ZQZGQkLV68mLAwutDllwYYLXT5ulgAs9DlL8noK7EydOhzFKMtWLCAJRrqEqnb+P3EiRN8yLNnzx6u3lCZncqNNFTQsXBCf8uLKAqwZJVjVegAe\/HiRZo1axbXRMmBv4KCAq4aQUGaqAQMKXStvpeh6lWnaUtX5H5WhI5yR6gT6PD+\/ftzUVrdunX5VRGBRf1ry5Ytaf78+e6iBSX0jx8\/urSbWkUlHcBwgK3aMFWmnMpishr0Tp060c6dO2nHjh0cTYVtLnOTVaeeuSqME0dxcbFLTqfAgFBAN9oHtBJhJeiHDx8mhLQBHWoFRoQovxGAUcJTVFTkxfvJkydcCAz937ZtW2rTpo365CgU0IV9j1mpvgDtM60EHXoaDhB8CVgiWuC+HLGQ6nT5wXrqRUg7kvll8CqnyUrQUW33+fNnSklJUeZq9uzZk7p06aJkr4T+9u1bl6iLF6NCIeni3rJHK36zehjAKKTgq4xdCV2l040eYnTd7EZqdL+quh6y0K5K0iv6kjZ03wRDEtq1oVcB9Ip+KVYZHzL1Yv9vdfpLHCro\/wDCVn0aqgzSpQAAAABJRU5ErkJggg==","height":56,"width":93}}
%---
%[output:9dc72aa1]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAACVtJREFUeF7tnAlIVVkYxz+jxWhhUgrbBm2yfaURjaygwsAoCiotxvaVymi1fd9tpWy3ncioYGoKoihCs41oX5VosqRB0yKdciZy+H8z53m8nnvfvc+b87R7QOi9d+655\/7Od\/\/nO9\/5Tj6FhYWF5JQyJeDjQC9T3nwzB3rZM\/c+6F++fKG1a9dSYmIi4\/D396f9+\/dTmzZtlHhQd9euXcrfKleuTM2bN6f+\/ftTdHQ01ahRg+vl5OTQ6NGj6d69e8rratWqRR06dKAxY8ZQly5dqFKlSrYOjddZekZGBo0cOZJevHjhetBp06bR5MmTycfHp8TDG0GXK0dGRtLq1asJQN1BF9dh0BYsWEAxMTHKe3s6El4H\/cyZMzR16tRiz9O2bVu25oCAAI+h48IlS5bQsGHDTEPHNXXr1qUDBw5Qy5YtPWVc4jqvgl5QUEALFy6kEydOUJMmTSgsLIyOHj3Knd67dy\/16NHDEDrq4hpRcnNz6ciRI7R161aCbA0YMIBWrVpF+fn5Lnlp3749S5mfnx9f9vXrV3r58iWtXLmSLl++zN9t2LCBr7WreBX0J0+e0IgRIygrK4sGDx7MVjlu3DjKzMykgQMH0vLly6latWrFnl2WFy10VERbo0aNokePHlHfvn1pzZo19OnTJ13oovHk5GQaPnw4f4QsRUVF2cXcuybSQ4cOsQSgwDp79erlsny919wIumPpbuwkLy+Ppk+fThcvXmRpgcfSuHFjOnv2LE2ZMoWvxqQGq5WL2YkU12zZsoWt3exEimuCgoJYfgIDAyuepd++fZulBXoLbwHaDu\/h7du3NH78eHrw4AHrdUJCAtWpU8cFwCx0yBUGrWbNmqahV6lShaUFeq7ynDwdBa\/QdEQitm3bRps2bSoxacp+OwZhz5491L17d1PQhb+NQezatatrPnBn6Q0bNuT60PRmzZrZChwd9wro2dnZPGHevXvXrfFgURMXF8dvAYq7iVTVoAxd67247YANFbwC+pUrV2js2LHs1rkrWo11oLsjpvhdu+w304RY5DiWboaWoo687G\/RogV7CvXr1y9RU55oIyIiaP369TwpOpbuAXh52S97Ldqm4HNPmjSJrl+\/zoErLM07derkQLfKXF7241q9pb5od9++fbRixQr+KCZULNFFlFG1InUmUqujUkHre4X3UkHZ6j6WZehYyFy9epVjI7du3eJNBizNEZwSmwTfG0Srz2sZekpKCk9ooaGhvBvz+PFjOnjwIHXr1o2WLl3KHoVTjAlYgi48CFg3YhICMBY3sbGx\/B12aJxiI3S4a5CRnTt3FttQEBFCbAQsW7aMqlat6nA3IGDJ0o8dO8aBKYRdg4ODXc1C52Hl0Hh5F8YhryZgCTpWf6mpqexPY1NBLvjt3LlzvGhBfERbkuP\/oA8Zf7m+Tov8lf70z3Z9zsqIoL8LikK23j5gP\/r5Ulzvks9ppt+WoUNiVNaclJTEHo32LUAn7ifl0m+xGa7+vA5Nofsx\/6ZYoLzP+pky06PN9Ner6sT1DvQIfJlAB3CAF+X+L4n0OizF9TkzPYreZ4V4FVAznRkSEkAJQ6xnCdgG3UheKqqlAzjAWy2WoJdmIoWm\/56aV0zTc4KfFmn66wjK\/\/CT1f7\/b\/XDm\/7gkbSgw5agOy6jPWNsCbpYHCHeLa8+ncWRtcHwuXbtmqX89Js3b7KvjlQ3bBC\/efOGzp8\/T9hrxG5+9erVrfXAy2s3atSI8Gdn8QkKCrIE3c6bl4e2EGOKj4+3FbxlSy8PoOzq440bNzhByezmiNn7WtJ0s41WlHpwHIYOHepAL8sBdaCXJe3\/7iVDb9CgAWcgwFP7+PEjIXMBmzfIjdRmEstdvXTpEk2YMIGQHCvSuB15MRhMAR05kNB2JLYiZbpevXoc3Dt9+jSDx0kRFXjkuSNzLT09vZhEubyX0uZgIwyAgpQ3VUFADBFK5IeXF7dSQIdrjD0DBPTEaRCEs2G9sH6RDiI\/NzIdsLeAgUFSrDwZs6WnpaXxOR+kJHuS\/A6gc+fO5exaFXTReZGUX96gt2vXjjp37lzi2Z4\/f87csGum5Xbq1CnOMMamDxaSJaBjhAAOsRUrmxA40TBnzhxCwhCKCrqcgVVeoeu5jHfu3GGokB8ZOgZj4sSJNGPGDD5Wo\/WAXJourB2vi3xuR0\/yBPBXr17xCnXdunWEyUa2dADHgCDGjpHHMRYr8qLNPVcNmnjLRD\/RBzmmL4wJCao4oSeKGd\/byHtBDubGjRvp5MmTxQ6CQYYWL15Mvr6+tGjRIsLAuIXuicSIAdBClwcMAK1AR315w0R1D9Xbqb1ODIqcEm3WFdSrBz3HYbT58+ezRYMZUrfx\/fHjx1ktdu\/ezac3VG24LF3kbCOtwqqu2w3dTHt6\/dW+sXo7WmaMQAUMYC9cuECzZs3iM1Fy4O\/hw4cssTiQJk4CegQdF82cOVO5DScs2QwkMw8pvxnCQrVyoSd3WikS8qE3VxltL4p7aIHhuCPkBBrep08fPpRWu3Ztro4ILM6\/Nm3alObNm+c6tGAIXeXBACbOU6Igqqj3BnwL6Lin6BNkCUWl1\/Ca5N\/wb3gUYm7Sc1XNGJMMrGPHjrRjxw7avn07R1Phm8temKhrtIYT\/TecSPHQyN4aNGgQbd68mQ+xikOucuOeQkf7OEyLlGdMSkZWLeQE50IxUQq42jlIJS8qr8yKpR8+fJgQ0gZ0yAoGVRy\/ERxwhOfp06KdMPH9s2fP2HCh\/61bt6ZWrVoV7RypXkN8hwILx+uLRYLKsykNdDwAvBJ4Pe4WULLloV+qYJSQJVleVFkKZuRO3A86jQUQ3ircUwvcyLp15cVIWnDCDIlFuBhxB9XipzTQZ8+eze4m7gFrhm+L86TIrUGRXUzZM3n37h1bnBgwrRzJ0AFLdje1A6MHTQDDabvPnz9Tv379lLma4eHhFBKizmZQQhebGNowgFZP0TG9k2ilgQ55EbIloGMiwm6NvPBS3V\/bR\/QPWo6JH2+keHsgLz179mQJU80N7qAbWTJ+MwqhGE6k2oZVsRR8B9\/Tqkup12lAky1d+9ndw5r53ZOVtmjXrD9vph9yHWWUUUxasBRZw91prtWba2XN7vbRn3IDHSMs\/nch2VuxGipwNwjCe4Gvi3DBtzhIW26gu4Nl1+8Cup4ratd9PG2nTOXF005ave57hf4P2dKjKeEHfgIAAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:04266fdc]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAABwBJREFUeF7tnH9IVUkUx88TrQ0VNiXwRy0aayUYBbFr4GogUhAEBqEmaD\/U\/KM0CkGFArPyByoSJSj9MgkzCP8J3QJB3Mo0CYxWA4VoEcVFSRaMcIne8h12ZLzO\/fGu7+p9zzsg+N78uDOfOffMOWdmnsvtdrvJSatKwOVAX1Xe7GEO9NVn7hvQv337RrW1tXT37l2GKDw8nO7fv08JCQlSZCjb0tIizQsMDKSdO3dSeno6ZWVlUXBwMCv3+fNnysvLo3fv3knrhYaG0t69eyk\/P5+SkpIoICDA9HT5hKRPTEzQqVOn6OPHj4sDvXDhAp07d45cLteywWtBFwsfPnyYqqurCUD1oPN6mLRLly5RTk6O9NlGZsInoD99+pTOnz+\/ZDy7d+9m0hwREWEaOipWVFRQbm6uYeios2XLFmptbaX4+HgjjJeVsT30hYUFunz5Mj158oS2b99O+\/fvp\/b2djaQO3fuUGpqqiZ0lEUdnubm5ujhw4d08+ZNgto6evQoVVVV0ZcvXxbVy549e5gqCwsLY9W+f\/9Onz59ouvXr1Nvby\/7rqGhgdU1k2wP\/cOHD3Ty5EmamZmhjIwMJpVnzpyhqakpOnbsGF29epU2bty4ZOyielFCR0G0dfr0aRoZGaEjR45QTU0Nff36VRU6b\/zFixd04sQJ9hFqKTMz0wxz+y+kbW1tTAUgQTrT0tIWJV\/tNdeC7ki6jpzMz8\/TxYsXqaenh6kWWCzbtm2jrq4uKioqYrWxqEFqxWR0IUWdGzduMGk3upCiTmxsLFM\/MTEx\/ifpb9++ZaoF+hbWAnQ7rIfp6WkqLCyk9+\/fM33d1NREmzdvXgRgFDrUFSYtJCTEMPSgoCCmWqDPZZaTkVmwrU5HdOLWrVvU2Ni4bNEU7XZMwu3bt+nAgQOGoHN7G5OYnJy8uB7oSXp0dDQrD52+Y8cO08DRSdtCn52dZQvm8PCwrvDAqSktLWVvAZLeQiprUISutF50O+BhAdtC7+vro4KCAmbW6SWljnWg6xGT5CvdfiNNcCfHkXQjtCRlRLd\/165dzFKIjIxcVlJcaA8ePEj19fVsUXQk3QR40e0XrRZlU7C5z549SwMDAyxwBdd83759DnRPmYtuP+qqufq83Xv37tG1a9fYR76gwkXnUUaZR+ospJ7Oih+Ut6314gdsVYdgKXQ4OK9evWIxk6GhIbb5AJcdQSu+eeDPcNXGZin0ly9fsoUuMTGR7dKMjo7SgwcPKCUlha5cucIsjfWYLIPOLQtIN2IVHDCcnuLiYvYddm7WY7IMOsw4qJHm5uYlGw08cogNgsrKStqwYcO6424Z9I6ODhawQjg2Li5uESz0PKQcOl7cnVlP5C2DDq+wv7+f2dnYbBAT8rq7u5kzg7iJMv1d2U\/\/\/vXPms9Da2IUTYfK38Sfwn6g0kPL+26k05ZCh4qRSfPjx4+ZRaN8C9DhubY\/aSLvmZG+W1rmWXw41aRpb1KUHooxBd520Cfyfqe5thFLgRppHMABXisd\/yWCmo57fiJgTaBrqRdfknQAB3hPk2XQV7KQQqfP\/zHh6Vi8Xr7110gajg6Vtvvbzz+aUi1ozDLojsmoLgOWQefOEeLgovfpOEdErtevX1t2Pv3NmzfMVscROGwcT05O0vPnzwl7kNjl37Rpk9dVgrcb3Lp1K+HPm8kVGxtrGXRvdnSt2kLcqK6uzqvgLZX0tQLlrecODg6yw0hGN0KMPtcynW60A3YuB2MgOzvbgb6ak+RAX03a\/z\/Lge5v0MfHx9n1Epz5RsLhTBxTM5sQ0EKEEee+ZWahXr7Wc1G3vLx8SRFZfxFqUAu4GR2XZZI+NjbmBnAc1MEJWH6mD\/+bAc87yg\/bK6Hr5WsBkYHk\/UU9MaJpa+jFxcXMThelEmBKSkqkoVc9KPy8iQy6ePJKbVLU2udvIxcOsZwsz9bQ09PT3dg0Fq9ycOlRfq8HHCezECPv7OxkqkqcSEDQytd75TlYXAYwcu3E1tCTkpLcSulZqYrBgJXQRah6+bIJwJ2gsrIyNnFG3hKfg84HGBUVZUqv60HVy1eTerWD+zKP0S+gK60GrdtlelBl+WpXVtQO6CstLmU5n4OuVC8Ajk0Jbh3o6Xwz0Ll087ZhOYn3P42YkaLpaGvoRhZSDABJNCG1bG0roGu1ydWhaIXZGjpMRqXuVpqMntrWVkDXclRka5CtoWMTA5E0rqPVLBfReoBEaVkQVkDHM7neVy6cMsC2ho4f2TETBpCpHK53rYKO9rnEizre58IAer9sBAnHDxHg\/qR4PA6Df\/TokWp8Rc\/Z0TMLPVlIzT5Lr55lsRc96OiYctFcqR2vNVgz1osePLP5awqdgxejeyuNRDqSblYU\/Lzemku6n\/OVDs8q6P8BpqwyszC8IKwAAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:448bf181]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAABYBJREFUeF7tm1kodG8YwJ\/5siRcIEUoo2xFlAtKKElSihJS9u3CFiUSZV9CEorskih3opSS7KRIKBdSloS4ISnx9Tz\/78z\/GGfMmTMzZuh9ay6Ydznv733mWd8je39\/fwfWvpWAjEH\/Vt60GIP+\/cx\/BvTX11dobW2FoaEhQmRnZwcjIyPg4+MjiAz79vf3C35nYmICnp6eEBsbC0lJSWBpaUn97u\/vISsrC\/b39wXHWVtbg7+\/P2RnZ0NwcDD8+fNH8nH9CEk\/Pz+HjIwMOD09VWy0pKQECgoKQCaTfdr8V9D5naOjo6G5uRkQqDro3Dg8tKqqKkhJSRFcW8xJ\/Ajos7OzUFxc\/GE\/vr6+JM0ODg6SoePAmpoaSE1NFQ0dx9jb28Po6Ch4e3uLYfypj9FDf3l5gerqapiZmQE3NzcICgqCyclJ2sjg4CCEh4d\/CR374hiuPTw8wMTEBHR3dwOqrbi4OGhqaoKnpyeFevHz8yNVZmtrS8Pe3t7g7OwMGhsbYWlpif7X0dFBY6U0o4d+fHwM6enpcHt7CwkJCSSVubm5cHV1BfHx8VBfXw\/m5uYf9s5XL8rQsSPOlZmZCYeHhxATEwMtLS3w\/PysEjo3+crKCqSlpdGfqJYSExOlMDd+Qzo+Pk4qABtKZ0REhELyVf3Mv4LOJF2NnDw+PkJpaSksLi6SakGPxcXFBebm5qCwsJBGo1FDqeU3sYYUx3R1dZG0izWkOEYul5P6cXV1\/X2Svru7S6oF9S16C6jb0Xu4vr6GvLw8ODg4IH3d29sLNjY2CgBioaO6wkOzsrISDd3U1JRUC+pzIc9JzCkYrU7H7ERPTw90dnZ+Mpp8vx0PYWBgAMLCwkRB5\/xtPMSQkBCFPVAn6U5OTtQfdbqHh4dk4PiQRgv97u6ODObe3p5a4cGgpry8nH4F2NQZUqEJ+dCVvRe1D6BhB6OFvry8DDk5OeTWqWvKOpZBV0dM4HvlsF\/MFFyQwyRdDC2BPvyw38vLizwFR0fHTz35hjYyMhLa29vJKDJJlwCeH\/bzvRblqdDnzs\/Ph83NTUpcYWgeEBDAoGvKnB\/241hVoT437\/DwMDQ0NNCfnEHFEJ3LMgpFpMyQanoqv6C\/0Xovv4Ctyi3oBDoGMmtra5Qb2dnZoSIDhuaYnOKKBL8ZoqZ70wn01dVVMmiBgYFUjTk6OoKxsTEIDQ2F2tpa8ihY+5+A1tA5DwKlG3MSHGAMboqKiuh\/WKFhTYfQ0V1DNdLX1\/ehoMBlCLEQUFdXB2ZmZoz7PwJaS\/rU1BQlpjDt6u7urgCLeh6lHHU8vwrDyOsg4YXR3\/r6OvnTWFTgN\/xufn6eghbMj7D2HwGtJR3BoooRkubp6WnyaJR\/BXz4FxcXgJ\/f0JydnQE\/6prBoWPlBj+\/oeGNBeVbC0L70it0MeqFSboEcWOGVHNoWks6cxkNAJ0LjjDfzY8+WXCk+jBkGxsbWt9P397eJl8dr7phgfjy8hIWFhYAa41YzbewsNBcHIx8hFhPRdCQyuVyraEbOR+9PB7mmdra2kS5iMoPoBNJ18uujHjSra0tcnPFFkg+QWdvYmh+uug8JCcnM+iao5M+gkGXzk7ySFXQMcmHV7oxF4XXr1U5EFr76ZKf\/AcPFIKOd3Uw14Rp7KioKLp+zaDr8JD50NGLOTk5oZcE8HYxNu7OO4OuJ+hYMcP3oW5ubqhUiTEKVs+YpOsQOE7Fl3SsIWAtAYNADJgqKytpNQZdj9D57zPhKzQVFRUMuo5503SqvBcGXR+0\/83JoOsRrqqpGXQG3QAEDLAkk3QG3QAEDLAkk3QG3QAEDLAkS+0y6AYgYIAltZZ0XdwGMMC+DbokZhLLysqkl+vYbQBp56fNbYC\/NHFOlcjA84sAAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:9e3bfa85]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAABYBJREFUeF7tm1kodG8YwJ\/5siRcIEUoo2xFlAtKKElSihJS9u3CFiUSZV9CEorskih3opSS7KRIKBdSloS4ISnx9Tz\/78z\/GGfMmTMzZuh9ay6Ydznv733mWd8je39\/fwfWvpWAjEH\/Vt60GIP+\/cx\/BvTX11dobW2FoaEhQmRnZwcjIyPg4+MjiAz79vf3C35nYmICnp6eEBsbC0lJSWBpaUn97u\/vISsrC\/b39wXHWVtbg7+\/P2RnZ0NwcDD8+fNH8nH9CEk\/Pz+HjIwMOD09VWy0pKQECgoKQCaTfdr8V9D5naOjo6G5uRkQqDro3Dg8tKqqKkhJSRFcW8xJ\/Ajos7OzUFxc\/GE\/vr6+JM0ODg6SoePAmpoaSE1NFQ0dx9jb28Po6Ch4e3uLYfypj9FDf3l5gerqapiZmQE3NzcICgqCyclJ2sjg4CCEh4d\/CR374hiuPTw8wMTEBHR3dwOqrbi4OGhqaoKnpyeFevHz8yNVZmtrS8Pe3t7g7OwMGhsbYWlpif7X0dFBY6U0o4d+fHwM6enpcHt7CwkJCSSVubm5cHV1BfHx8VBfXw\/m5uYf9s5XL8rQsSPOlZmZCYeHhxATEwMtLS3w\/PysEjo3+crKCqSlpdGfqJYSExOlMDd+Qzo+Pk4qABtKZ0REhELyVf3Mv4LOJF2NnDw+PkJpaSksLi6SakGPxcXFBebm5qCwsJBGo1FDqeU3sYYUx3R1dZG0izWkOEYul5P6cXV1\/X2Svru7S6oF9S16C6jb0Xu4vr6GvLw8ODg4IH3d29sLNjY2CgBioaO6wkOzsrISDd3U1JRUC+pzIc9JzCkYrU7H7ERPTw90dnZ+Mpp8vx0PYWBgAMLCwkRB5\/xtPMSQkBCFPVAn6U5OTtQfdbqHh4dk4PiQRgv97u6ODObe3p5a4cGgpry8nH4F2NQZUqEJ+dCVvRe1D6BhB6OFvry8DDk5OeTWqWvKOpZBV0dM4HvlsF\/MFFyQwyRdDC2BPvyw38vLizwFR0fHTz35hjYyMhLa29vJKDJJlwCeH\/bzvRblqdDnzs\/Ph83NTUpcYWgeEBDAoGvKnB\/241hVoT437\/DwMDQ0NNCfnEHFEJ3LMgpFpMyQanoqv6C\/0Xovv4Ctyi3oBDoGMmtra5Qb2dnZoSIDhuaYnOKKBL8ZoqZ70wn01dVVMmiBgYFUjTk6OoKxsTEIDQ2F2tpa8ihY+5+A1tA5DwKlG3MSHGAMboqKiuh\/WKFhTYfQ0V1DNdLX1\/ehoMBlCLEQUFdXB2ZmZoz7PwJaS\/rU1BQlpjDt6u7urgCLeh6lHHU8vwrDyOsg4YXR3\/r6OvnTWFTgN\/xufn6eghbMj7D2HwGtJR3BoooRkubp6WnyaJR\/BXz4FxcXgJ\/f0JydnQE\/6prBoWPlBj+\/oeGNBeVbC0L70it0MeqFSboEcWOGVHNoWks6cxkNAJ0LjjDfzY8+WXCk+jBkGxsbWt9P397eJl8dr7phgfjy8hIWFhYAa41YzbewsNBcHIx8hFhPRdCQyuVyraEbOR+9PB7mmdra2kS5iMoPoBNJ18uujHjSra0tcnPFFkg+QWdvYmh+uug8JCcnM+iao5M+gkGXzk7ySFXQMcmHV7oxF4XXr1U5EFr76ZKf\/AcPFIKOd3Uw14Rp7KioKLp+zaDr8JD50NGLOTk5oZcE8HYxNu7OO4OuJ+hYMcP3oW5ubqhUiTEKVs+YpOsQOE7Fl3SsIWAtAYNADJgqKytpNQZdj9D57zPhKzQVFRUMuo5503SqvBcGXR+0\/83JoOsRrqqpGXQG3QAEDLAkk3QG3QAEDLAkk3QG3QAEDLAkS+0y6AYgYIAltZZ0XdwGMMC+DbokZhLLysqkl+vYbQBp56fNbYC\/NHFOlcjA84sAAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:139a5ea0]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAABYBJREFUeF7tm1kodG8YwJ\/5siRcIEUoo2xFlAtKKElSihJS9u3CFiUSZV9CEorskih3opSS7KRIKBdSloS4ISnx9Tz\/78z\/GGfMmTMzZuh9ay6Ydznv733mWd8je39\/fwfWvpWAjEH\/Vt60GIP+\/cx\/BvTX11dobW2FoaEhQmRnZwcjIyPg4+MjiAz79vf3C35nYmICnp6eEBsbC0lJSWBpaUn97u\/vISsrC\/b39wXHWVtbg7+\/P2RnZ0NwcDD8+fNH8nH9CEk\/Pz+HjIwMOD09VWy0pKQECgoKQCaTfdr8V9D5naOjo6G5uRkQqDro3Dg8tKqqKkhJSRFcW8xJ\/Ajos7OzUFxc\/GE\/vr6+JM0ODg6SoePAmpoaSE1NFQ0dx9jb28Po6Ch4e3uLYfypj9FDf3l5gerqapiZmQE3NzcICgqCyclJ2sjg4CCEh4d\/CR374hiuPTw8wMTEBHR3dwOqrbi4OGhqaoKnpyeFevHz8yNVZmtrS8Pe3t7g7OwMGhsbYWlpif7X0dFBY6U0o4d+fHwM6enpcHt7CwkJCSSVubm5cHV1BfHx8VBfXw\/m5uYf9s5XL8rQsSPOlZmZCYeHhxATEwMtLS3w\/PysEjo3+crKCqSlpdGfqJYSExOlMDd+Qzo+Pk4qABtKZ0REhELyVf3Mv4LOJF2NnDw+PkJpaSksLi6SakGPxcXFBebm5qCwsJBGo1FDqeU3sYYUx3R1dZG0izWkOEYul5P6cXV1\/X2Svru7S6oF9S16C6jb0Xu4vr6GvLw8ODg4IH3d29sLNjY2CgBioaO6wkOzsrISDd3U1JRUC+pzIc9JzCkYrU7H7ERPTw90dnZ+Mpp8vx0PYWBgAMLCwkRB5\/xtPMSQkBCFPVAn6U5OTtQfdbqHh4dk4PiQRgv97u6ODObe3p5a4cGgpry8nH4F2NQZUqEJ+dCVvRe1D6BhB6OFvry8DDk5OeTWqWvKOpZBV0dM4HvlsF\/MFFyQwyRdDC2BPvyw38vLizwFR0fHTz35hjYyMhLa29vJKDJJlwCeH\/bzvRblqdDnzs\/Ph83NTUpcYWgeEBDAoGvKnB\/241hVoT437\/DwMDQ0NNCfnEHFEJ3LMgpFpMyQanoqv6C\/0Xovv4Ctyi3oBDoGMmtra5Qb2dnZoSIDhuaYnOKKBL8ZoqZ70wn01dVVMmiBgYFUjTk6OoKxsTEIDQ2F2tpa8ihY+5+A1tA5DwKlG3MSHGAMboqKiuh\/WKFhTYfQ0V1DNdLX1\/ehoMBlCLEQUFdXB2ZmZoz7PwJaS\/rU1BQlpjDt6u7urgCLeh6lHHU8vwrDyOsg4YXR3\/r6OvnTWFTgN\/xufn6eghbMj7D2HwGtJR3BoooRkubp6WnyaJR\/BXz4FxcXgJ\/f0JydnQE\/6prBoWPlBj+\/oeGNBeVbC0L70it0MeqFSboEcWOGVHNoWks6cxkNAJ0LjjDfzY8+WXCk+jBkGxsbWt9P397eJl8dr7phgfjy8hIWFhYAa41YzbewsNBcHIx8hFhPRdCQyuVyraEbOR+9PB7mmdra2kS5iMoPoBNJ18uujHjSra0tcnPFFkg+QWdvYmh+uug8JCcnM+iao5M+gkGXzk7ySFXQMcmHV7oxF4XXr1U5EFr76ZKf\/AcPFIKOd3Uw14Rp7KioKLp+zaDr8JD50NGLOTk5oZcE8HYxNu7OO4OuJ+hYMcP3oW5ubqhUiTEKVs+YpOsQOE7Fl3SsIWAtAYNADJgqKytpNQZdj9D57zPhKzQVFRUMuo5503SqvBcGXR+0\/83JoOsRrqqpGXQG3QAEDLAkk3QG3QAEDLAkk3QG3QAEDLAkS+0y6AYgYIAltZZ0XdwGMMC+DbokZhLLysqkl+vYbQBp56fNbYC\/NHFOlcjA84sAAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:3075d01f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Limiting legend entries to 50. Specify a vector of graphics objects to display more than 50 entries."}}
%---
%[output:80619bc2]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAABrFJREFUeF7tm3tIVFkcx3\/j+thS\/0hXaNJyJrBsK3STtiAqqAhWECqkJOihpf5RKkbbY7fASntLGyko2hORDGH\/kAJ3AwnzVckWYe0SRCBKrKKwaCIrunzPdmbPXO\/Mfey944xzDwjjzO+ce8\/n\/O7vdc61TU1NTZHVfErAZkH3KW92MQu675kHBvSJiQm6dOkS3bx5kyGKjY2l27dv04oVK2SRQba6ulr2t9DQUFq6dClt27aNsrKyKDIykskNDQ3RgQMH6NWrV7L9oqOjKTU1lQ4ePEjr1q2jkJAQ3csVEJre29tL2dnZ9P79e9dEi4uL6fDhw2Sz2aZN3ht0UTg9PZ0uXLhAAKoEnffDop06dYr27Nkje201KxEQ0JuamqioqMhtPitXrmTaPH\/+fN3Q0bGkpIT27t2rGjr6xMXF0Z07d2jZsmVqGE+T8Xvo4+PjdPr0aWpsbKTFixfT2rVrqb6+nk2ktraWNm3a5BU6ZNGHt+HhYaqrq6MbN24QzNb27dvp\/PnzNDo66jIvKSkpzJTFxMSwbpOTk\/ThwwcqKyujlpYW9l15eTnrq6f5PfS3b9\/S\/v37aWBggHbu3Mm0Mi8vj\/r7+ykzM5POnTtHERERbnMXzYsUOgQxVk5ODvX09FBGRgZdvHiRxsbGPELng7e2ttK+ffvYvzBLu3bt0sPc\/x3pvXv3mAlAg3Zu2bLFpfmeHnNv0C1NV9CTkZEROnLkCD1+\/JiZFkQsCxcupIcPH1JBQQHrDacGrRWbWkeKPtevX2fartaRoo\/T6WTmx+FwzD5N7+7uZqYF9hbRAmw7ooePHz9Sfn4+vX79mtnryspKmjdvnguAWugwV1i0qKgo1dDDwsKYaYE9l4uc1KyC39p0VCcqKiro2rVr05ymGLdjEWpqamjjxo2qoPN4G4u4fv16lz9Q0vT4+HgmD5u+ZMkS3cBxk34LfXBwkDnMly9fKioPkprjx4+zpwBNyZHKDShCl0YvijegUcBvoT958oRyc3NZWKfUpDbWgq5ETOZ3adqvZgie5FiaroaWjIyY9icnJ7NIwW63T5MUHe3WrVvp6tWrzClamq4DvJj2i1GLdCjE3IcOHaLOzk5WuEJqnpaWZkHXylxM+9HXU6rPx7116xaVlpayf7lDRYrOq4xyGanlSLWuyiyQ99voZRaw9TgFU6EjwWlra2M1k+fPn7PNB6TsKFrxzYPZDNfT3EyF\/vTpU+bo1qxZw3Zp3rx5Q3fv3qUNGzbQmTNnWKQRjM006DyygHajVsEBI+kpLCxk32HnJhibadARxsGMVFVVuW008MohNgjOnj1L4eHhQcfdNOj3799nBSuUY5OSklxgYeeh5bDx4u5MMJE3DTqywvb2dhZnY7NBbPjt0aNHLJlB3UTaxv74iSY\/9bGv\/3qwmib+jGafw5N\/pZC5wz5bn5+X22kgUv5JtM+1U87yfF33Yip0mBg5bW5oaGARjfQpwAzGextp9Lfv2WRGW5JpuGIz+xy66AV9uapB1yT1dGp1xFLNt4leu+Z8nacLvN9BB3CARxuq2EyfWpLZ54hVDRS26IUefrr6ADjAe2vfOTLox9X\/biVqaTMC3Zt5CSRN\/2F1CaU7MrTwZrKmQf8\/jhQ2\/e\/BLpdNH++J\/2zTf6EvvvrvwJHm2WrsAJv+e5x8LvFNXJou02IqdCtk9LzCpmk6T45QBxezTys5IrJ1dHSYdj792bNnLFbHEThsHPf19VFzczNhDxK7\/HPmzNH4wPtePCEhgfBnZLM5nU7ToBt5ozM1FupGV65cMRS8qZo+U6CMum5XVxc7jKR2I0TtdU2z6WpvwJ\/lEAzs3r3bgu7LRbKg+5L252tZ0Gcz9Hfv3tGxY8fo8uXLbiVZPmd+\/AzH2MTD9mqYyJ0XXLBgwbTCl7eCmJrrGCXjE03nUHBoXq4CiIPzJ06cIJxL0erRsZh4bwhHnMXD9AB88uRJt0P2QQOdryq0RE77ODS8AYGmFTqKXOiLtx6kSZH0t6CAzoFjRycxMZGOHj3qpukcOA7P79ixg2ksjrBpMS8A66m+LjUHQQFdnDTASKGLv\/MF0ApdfJKUnhILukQN9ULHMNx+i0PyF61Ek2NBVwFdau\/FLp7eQpO+oiLKWdBVQBdFAAwb0nIOUy6Uk4uYLOgGQFcySdJ42IJuAHQe32MouSdA6rwt6AZAxxBc23FoXwQvlzRZ0A2CjmHEjJYPG\/RlACNqFlodqRHXNGsMn9RejLh5C7oyRcN3jizoMwBd+ZKBIxEw5iVwkCrfqVnQ\/wGWa7ekK9XRjwAAAABJRU5ErkJggg==","height":56,"width":93}}
%---
%[output:054c6af0]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_med_h","value":"1.0772"}}
%---
%[output:24547367]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_STN_h","value":"0.2640"}}
%---
%[output:25e48712]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean_h","value":"1.0958"}}
%---
%[output:5075afc3]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean","value":"0.8914"}}
%---
%[output:9ef7d6fd]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_median","value":"0.8763"}}
%---
%[output:82396b59]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAABIdJREFUeF7tnD1MFEEUx9\/ZCBRGURs+hMRAizaCaEJj7LSgIVCZEELQgoIghJgQEgUuhFKRglwHoaFQKzsKEQqjtCQmRKFSTqNGqDjzRt85bHbmZvZmP2ZvJiEhy\/Bm5rdv\/+\/N3rzLFAqFArgWKYGMgx4pbzaYgx498+DQ9\/b2AH9sbE\/enYbP347g1O+v8OzuRSNLIJuN56rg5YOrUpuBPB1hj46OwtbWlpEJR23kx+0sHNdcYNDPvB4zMjzZvFRbBR8eXS8NfXNzE\/r6+oodZ2ZmoKenB\/jrdA070fW5uTmor683MukojZBX4pjPb2WMDH3\/xRf4lD+ChoYGNeirq6uwu7sLY2P\/73o+n4eRkRGYmJhgk5qenob5+Xmora0tQl9eXoaOjg4jk7bdCDmiChMmL9lsFpqbm5l3U0MjeH1paQmqq6thfHwcent7GWSdAWyHqTp\/HSaZg4ODQn9\/P2xvbzP7bW1tDPTOzg6srKzA7Owsu47QOzs7T8iOyl1VnbTt\/bSge\/N0lJqNjQ3o7u6GtbU1B13RG8qCTrIyNDQECwsLTl7CgI7yMjU1BZOTkyxIoo5jGxgYcIH0H\/A7T9+zvF6Wg2t7Op8akqbzWQqOzeu3zgCKjpLoblcev2XpILYbl8\/6bn50mATaHOkMkGSaKh6Mfd58\/F5chmjzo8OkoqGTByNIlA4\/CeH7EHm\/HaeDLnm8eO\/GbgiaGkqI15OpP\/aRbe8ddAl03nN5iKpwRaYddI4MDxMvU0BUeTGlE48cdI4Wn3ngZYQt02Yd0HxfB93j6ZR9EPBSr16DgHfQPdREOh4EbiI1XSX3NblYFVvlBkmVMWL19Ci8SgVC1H1ihR6FV4mAxvmUxQo9LA9TARrnU5ZK6CpA43zKUgNdtGUPI+Ur9wlNDXR+Y8NvalTOlpQLUff\/UwOdPJ227gTC9BZeF7Bff2uhi4Kl9\/2J942fSpA1AVZmw0rofh8W0CJLyYlKkHXQAUD09s8PTik5iTNrofla4emyt39+4JOYsST+LaNXd\/2gJx2sdZruTf9syEh04kAi5cUbKL0LKqXbOgDi6BsLdK98+KV5\/EdlNup24uTFe5yBPyuSNq8OdXPEn\/zyKwrgT33Jdo9p8+rQoKsWBfBVDTgZPwmxOStRjQ9GNF2lKGB4eBhyP68VQeMEsY6H2qt7japztr7f\/v4+q8NSObMvPFaH0EVFAXyh16+bD+G45vwJaKaKp2y7E+3t7YB1WFh3JGuBoKNBm0saw7qZCLsUcBxbCl1UcxTWpCvFrhC6LJBWCpyw1ik9Ks2njCoBIqxJps2u9vl0Ue6eNjB+60G5XVxcZH+qq6uDXC4HLS0twiJnERMt6JUsOYeHhydqaQloECZa0GW5e9o9nYeL3k0tCBNt6KLcPe3QeVnFtQ4ODrKyftl+xoi8BBkgjTeDpAYryJuamoSbSGPQXe7+FyXV23Z1dQm\/Q8EI9CBBIy2ejk\/5+vo6kxTkgN+ngL+3trYKi5yNQEcjlZy78ykjaXoQJlqBNC1eG\/c6HPQY7sAf+bRMpJ1JME0AAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:5fdc60cc]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"fit_test","rows":2,"type":"double","value":[["-1.6968"],["1.0405"]]}}
%---
%[output:0c16ee21]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BsB0_homo1\nTaille de la serie    : 198\nStatistique T_max     : 11.2430\np-valeur (bootstrap)  : 0.0560\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:23593bf4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:50842089]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BsG0_homo1\nTaille de la serie    : 198\nStatistique T_max     : 11.2132\np-valeur (bootstrap)  : 0.0560\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BsR0_homo1\nTaille de la serie    : 198\nStatistique T_max     : 14.4693\np-valeur (bootstrap)  : 0.0090\nPoint de rupture      : 2023-05-01 00:00:00  (indice 167)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BsR0_homo1\nTaille de la serie    : 166\nStatistique T_max     : 21.5133\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2011-08-01 00:00:00  (indice 27)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BsR0_homo1\nTaille de la serie    : 31\nStatistique T_max     : 8.5848\np-valeur (bootstrap)  : 0.1070\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BsR0_homo1\nTaille de la serie    : 26\nStatistique T_max     : 3.8376\np-valeur (bootstrap)  : 0.5430\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BsR0_homo1\nTaille de la serie    : 139\nStatistique T_max     : 2.7607\np-valeur (bootstrap)  : 0.9220\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:9777a6f5]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAGTVJREFUeF7tXAdYVEfXfpciTRFY6SBgAVsEPxuWKMSCnxUrsSAg9o4FsAIKAipWFJUPAYMoIeazJrYEjA0VC1bUqFgoSrXQy\/6eWe56WRbEFvM\/X+Z5eNi9d+bO3HfOnDnnPWdWIBKJROCVlBTA1hag\/87OQHg4\/67sz1Q\/Ph4wNQUePXp\/\/U+pUZ6TiowgB2i7bkRm2BzozY+BvJbhpzzyL28rkAbdxwfw9q4KYDyAUwDoP1dsADgBMK28wAFPbb28vtx7cKBrjV6JjAB7oPgVOKlRMu8G\/SW\/QE5V\/csN4DM8WdDJP0F0\/0UB+rTUgpdNa7Q2V2CPbTDyKuTaFkF35r9wT0Opxq4I9DgAKfHiFUKl6tr5DKOUekT6muEovPhzjQ9WEBrBcFXC33YFCELin4gcu+hjQuQtdKpnhgXO6lDRKkHw3XK4NlJhL2ZcIYKrnAA9ea9Jku\/N+07Au5iJ1VJcHGBDS+ELFSbtgUOZhOt7\/FcCLrcK\/u4qR3Didpaod0shtp16il1bVHA+thEaDctE1j5tBlnryxnYoKoIqsPpeg5LQ2MRTOOA3fICdsnUBUiJEO8DtB98icIBK1BSRdHN3wBRRbVuFBo1hqHfuVol\/dmzZ6C\/r1EE11Nfi74xqM9ADw5SxK1DemJ9YQOMrahAadgN2JprYWpPY5hVSjJ\/oLMXlmKPSzkyWyqLlf5bFUN6fb7dPrwIGsGqao7ygebI5dXe701C9ToEauqSrijLelJrW7qZG7sCuT++20Bq6ke6YwJ74cKFuHDhwtfAHFVAD48ALqYYM9Abi0S4XVrO1I4EdNsUpJhWbqekzPn6hlMnPsB46wJY5qzCQOtgmJTq40XIRBgsO1pF8gjctJX9oDPtP+zFuTqZEfOhZNyKTVLxvQtI97NDo2lhqG89\/LMBlJCQgDFjxmDNmjUwNPzrLR9Bw7m\/ieQEAogqRDDONMMNczPAC6h3KR2qe+6gha4aVg5uhmYqQphFRgBeLnV++bD8dXAWuCLdrz8aDnKrAhxJeU7kPKYGBMoNZNapKHgl8zoNgJuQioKXkvHURa1QZQ706OhoWFtbV3mfFKQgHvE4hVMwhSmc4MT+f84isNtwSbRvmhWT6EndjeBUTxNpTeUw8\/QrGD85i23LZsBtix9mX7wHeEWQncLrnwYj\/Z2UOS15AXaWhMCpbAwDTtWqbxUVQ6C\/PLSemXhUZNUhYGWtEk4tyQsbozxbrIaocN915v9U68qQBTqBHYnIt8YB3zwQP9cZzthauBWenp7se0BAAFRUxEaGrFJYWIhNmzZh0qRJ0NLSqlZFYL7stCjzTSlcuugjaGQLWBWUIklVEQrLrkFNeSO0wk6haFcZ0rs\/\/YDJFk+GqddOPPAY\/lGg89WPknlnSd9M+lfaoTw\/D\/VMLaE3b6\/kXsa671GSkgR5NQ3oLztWo70uC3Qf+DDASaptYIOe6InHeCyZBMdCR2h6arK+Jk6ciGnTpiEtLU3S95QpU+Dh4cG+BwYGstUUFhbGQM\/JycG1a9fw3XffISkpCQJNt9+Zb9FcRxUHZ1jBZv8+FHstZY3zpuRBY7sG0qLTUGxdDLkcOei66kI5SZnd567TZ81ATVaXyutBr5EdkA3RhZ1wP2yHm\/E2uJNTyu5xS5qk1X\/JfOx5JMeut9RSxLZVS2DcdzwyLp2A65SpuJMndnv8\/f3h4ODAPmc+uAXnUfa4k1te5d79+\/fh4jQe\/YWvcDathN3nA8GXGA70Ll26oKCgAN5h3rDSsmKAxyFOok4IrHGu45CclIz8nvmwhS101XUZ6KdOnYK2tjYWLVrEHt24cWNs3rwZt2\/fllwzMDBASEgIWxnZ2dnQ19fHlStXIPFIyXo5cjYJyYf9kT1kMHI9dKEZOK8S9L4osZwLoee3rAMCVPWgKjTXayNDLhzy6\/rBYIwBm4Spwhs45jwQT2enomBwAYQdZsFe6xesjzmKI6cT2cDCw8OR+eguwle6Yc2WMKQ+z4Tr3IWwHzUGbqMHwdd9Bm6KDBAWHoETJ05I2giFQkxwckSbiqeY59AXpw\/Hwi1Rnk2khkIZnMc7ol0zYyxoWY6LFo7YsjOK9dW8efMqq1QadIswCwRqBTI1Eo53vAdJ7KFDh6Afro\/92A8LFwv0at+LgR4VFYXY2FjWNwE+ZMgQ9OrVC\/b29hg9ejTrz8nJCZMnT8aAAQOQl5fHVse2bdvegX7yTjbcAnYAyYeQHBWJsubRUErwrwQzGuXC09Bz2YTcWbnIPxoMhbuK0HsagszMtW+JFxtWj0rb6LbYb70fAgigcF8BeiObYYxxewQeikJy4hm4ODsjaN0GdO07ECTtfm6TEPNcvEImT3CGY+FxBP9+m11roVoEv2bP0dzzR6ajzx0\/jHFTZ2NtszS0a1iMwgoBglKErO1InTx4P9RDL0M5TOuoiye9lrO60pslSS+tmgcPHlSZiFKDUswOn42FzRey61TP1dWVbbS6HrqYUDgBHT07whrWDPSxY8fi9evXrK5AIECPHj2YnpeXl2crgNTJoEGDoKioiKNHjyI\/Px979uxhKkeQ8apIpNtACfNjk\/Ho4lGknvkRN1cFo8B2PpQSzktAB+bCYMyLKgOlLy+srVEQthdCz5ZocKhQcp+kngo3GfyGpC6o0NIkFUAbDveCpBeZqnBxkehMS0tLNlgm9Rs3YFWTVBgVP0WRDNDH6edicEshng3dIhN06rcuks4HPdkjGZGFkbDxtME3+Abjxo3D8uXLcffu3Sp4WFhYMInOyMiQXKcJIpXDgc4kna\/TJ+r\/CX+fZagYPwuPvee9lXSlStDXoVw4EnouemJJLxwJmKYANmS5kLVCO774Ty1GDTqLdFBkWYTspdnQmaODjvlDsCvEVmJR0O7OtwToOx90bsRcPVriNDk9e\/Zk9jUnvfzncJvbrFmzmCTXZhZKgz4hbAIGaw1m3ZJOp42UA13YWIidATvZvVGeo6ALsU739fVFSUmJZLPkVoePjw\/bPEnSSVBI+ql+FUnns4wkYaPHOUHQ8hskRuyQbI5p0etQYjkMQk\/xUs4OaAXVg4bQ3PwQGeE\/QOlKO2hu1kRGeAbKjcoxfnI\/HNG6hKzlWWi0ohH+ldgMEdP74A+F1hL9\/PPPP0t2eJJgTupnz55dZULIeySpJzD79OlTZXL4wJK+5+rVFXS+6nGBCyJAJjGJjjdMYIJ9gftw7dA19l72sEe6Szrat29fTaebm5tLxkWCwel0etaKFSuYejl37hz7TBZPNWp3S8xxBC2ayjrnrJfUoH0oOZgElJpB\/moTKD5MhbLyBRTs7AGDoe2QAgFMA00h2C7mYPSFGrgYfQNlzcuYxdO9Z1s8zRc7MdyL8lUI7fKkB2lDop2ek3wyr6iQbuRsY04CuXvSz\/sQSZfW93zgqV9BoYAJWoNDDdgqo6Kuri4ZS0xMjMRS4Y+REwaqTxs5CRgV7h2qgU4b6pIDfyLd3QFPnwiAeC8xaW5DvCKpk0pnKIVdxFt\/AugZ\/5ZW5LPt71Sdqagxziw1ruaRVlGGH\/GF81aL750F8eg600KRETSK0Qp8u17Wo9\/nkZKTRM4S2eqc3f4RQ6yxSTXQn78uhlP4TTSxvYyQBnsA03AxkWVDZHkKG4QX8zghdh6IVoy3qdTx8czsiqvojsdyvqy+4\/UZ8As\/VY17+ZSX4ACXFxpBc6AbcmJ9oOu2Fy+PbEDBtePvDWTUBvqnjKuubauBTg13n8zF+M0qqDjQUizZJNWmFIcTqw++PUsSYQYzwDYO8IqEs40XIpj4e8MoE9jtBfxrfO1ueV0Hy9Xj6N0Gtk7I2j65WnMFLUMY+l+okdr9W4DO101knmWq\/YizvqmAqdmH4eFDwL+N9yEe2vnaWK22mk3QlygZAUNQmHQMDYe4o\/hBIho5BSHdrx\/qNe9ShRqoq3rhiK6axkoWDZ\/4Isdp+\/btrDrtSZwTRtepcJQA7UFk0bi7uzOHisxjwb1790TugSF4YjIcEa6WCPY5hCilf0Pe4hHyR\/fHdwkeONnf7d1YyuUB+ZPi6Gj5MwgqerN7IrnuwGljtjIMu92H2S4zHBh6QCbh8zkmgfHu3r1Rnp4siZHWlU+XRe0SqyiL7OLGSvdIx1Mhu5s2cs5KSU1NZRaKo6Mjjhw5wlx9uteqVSvmQEVGRuL7779nfsbAgQMheFNcKiKGMfHJa+yfZoWAxbmIHayFVQ0PY2vOajR5uhU+4+tj\/iVXXOr8GMvOncFKA0NUbCb+JQWGxnGAUwTGVThht8ZKmGwMxYkl3RA6ZwQG5ooZxLqC8aGTQUGM4qe33yvZ0s\/lQO\/cufNXCWSwGGnGqxIkPMrDmhEWsBv7CgVBwATVrTj96DeYPN6CCV30cWxnALQa98Oqne1QMmo74EQ5FynwfeaLspwzMP95G1Z+3xZ3lF\/h1\/p+0Jo5G5k93NC7d2+Z9OyHAixdXzrCxL\/\/Pl6db9J9jUCGYGLkTdGS\/k0wNfo2ejdIwbpdXaC1Vg4dtJfg1q3bDPSezTRwZfcPiIAX5IL0UaHVj72jdmER5AU2yFBOgFGmKZ5pAyaP4mD0ShXFisex+OJmDBk1tsZAxMcC\/ymAU5980GUFMj52XHVtx8J1Og0U0dv3GJxalCJicyvkK7qiw9jGuK57A62zI2HZUAcxfhnIuzoX6UpKCHchszweR3VNsMBNiHyHK2jsZQqVMxq4v24S8XFYrtgfk0v9MaX5RJl8el0HWFs96RgpV\/dDJP1TQJfmiKj\/2p7HTbaEe2G6V00RwnPqyHg9DYVGffBmyiO0uugK1TPfAJ1+xjPXKAY2gX5fQQGORka4GSBCvsNjZh4avrBBp7uzcKFHZ2jdnIjX+hr4w3j7FwGdH+TIPbwe9buOhIq5NZ7Mt4K6rROE49fWOF+fS9IJ9GPHjmHmzJmsL\/pOm+aSJUtkRpao38ePH0Nw\/vx50ShHZ8gpKKNcQQXF+eZQUn2A7hOdsXfuJKhu+h7C\/1yGwpsy8abo\/AJ9Cgrx22+6yM4WBzOMu+ahzU3gvwY6EF1tg2YaCXhu0ggdIppg0ktfJLsPQWRyEVY1y4CBphq21BuII\/HnWFv3juroW5EEvnSS83NuYR8sOJmGF6WKrB4nQRxgUds2weRUAP5oOg6bg7cg6N9muF+\/Bbbu+hEOLepjc1IBa7ds2TJGPhFpxrGV9+7dY8TZ+yTzfatRFuh0rX\/\/\/lVUGL9fBjoRXqGRu7Fm3Uas3xqK7DclWDp1DvLmTQIqsqB0PgItOmxDswUPoGo1D6fK1PDfuxlINQQ8BXpwy8xF57f8sX23hii6UR\/n4g4iOfEcpnrMxJ0fZ2BkYj5UgrZBpWV3rNnxA1ZPd8DFq9cRse8IrsUGY8qOYwxQ03MbJFkA6buXw3vXL1Bu0Q1uBulswq7lyTPWjgNsfftStKl4hl+z1BCVoQm\/Zhm480YJ657qYGGrMjjuiMOSwE0MbHo+R4gRR8KxlZ8DdD4FTc\/jolXBwcGws7OrEkCRqBcCnZwjLqJz\/mEefBfORJriPNhfGomzSw\/gmUs0LFzuoCCtAKqqpQx0Ki56ehgtn4ssW21EZDVGk5t5iL1MSTPAwLG98HT4UOikHoXy4SLMXuDBWEKXsQ5oU\/4EXjsP4v56Z3gmlmPc9PkYYKLAAtU6c2PwZO0orL2njGNJj\/DvLlaYZ\/QcZl6\/spgn35vsaG6M0GkDEPVUBRvHWOP8T6GIStfEdj9PWI6YUSVWSWPi6OO6gC4mPcS0k\/iNqhdpSSeizs\/Pj0WMGD6VMQGZki4Nup\/HDIwc4gb\/ODvo1neFam4aynY9h8XmLDz5QxMheRlINQD8X+nhzSh5CETe6L9xA9IavkTY5cvMIRoXOA7308rwPF0TRmk3mMdGxclhGDLy8qu8AUnHjJ7NWEqG3pLjyAxxxb2mg+Hqu01Sb6BdLwSu28ScEj6nzh87OSWc8FCIjh8g\/lDQyRfnQK8pEVkW6LKyAGiMVExMTMQ6XZakr1gwHR4L5kJdzRSL3B2R5h6NLlOtcSnGrpJDHwDSL3ozYyD\/xgXxj9ciRvMlEho2lIDuneCNnxb8hLISXdh1acNozazTsZi62A\/f2o+F+7y5VUxJLg+GA53Lk3l6fBemui\/HnVdy8FvqDtWj\/iwuurZ5GqzUi\/FLphqTbtovbucrIeqFNiJ2\/YAWHbrXCDqF2rgMr5qsDU7SCazaQJdWL9zz+DSB9F4iKCgoEFEU5\/Lly0was\/NLMN7JBVOmz0SLNpaY4TQNGnNmo19qX8Rcao+yl2owzhJHwAsE3pgkSscAgTZ2ALjzdlMlSS+JmA4fzVj8uq8x5G4LWDS\/r9xd5j1GlrZjKiI0eAMSvIZj9u85TOe2QapEvdwLGMHUTvf+w9kKiP+Pv6QeBaApzjpzykR0vBGCwMsFuPVG+R3o6ZpY3VGAbpvOY21olCRQwpd04kX27dvHgP8Uk\/F9G21N9wXdunUT2djYIJ5MwUoVQNGjJ4Z2KDPtAdvM07h6VhxRUfq2NYpP30JFdDQElpYQeXaA3KFCVAyqAHQmQZCYiMuVuR6UR7Jx8Qs0OZ6E7c4d0URXkyUb8UNw9ExPW1NMDvsdlLPCpdOR\/Z18LRELT+VI4qTSeSVENukolqJTwwJcfKn6QZJOoH9NplEmtSs9Q7bxgJcNC1mwKD8xh7WlmhFjF\/E2Oc0m\/hGOGojzEWWlv3HpdFwggkvolw5Q8BP9mTe6qDPKclJrFLT3OUfU8KuCbm9vL6LNqaZwE6mGzp0d4GL6bhcnUCm0xWVC0QTQNWLqKA\/QCV5vP9kweuZL5KlzfLpcAyFERa8hp67DSC+OeVQ0MIe+54FaV\/9XBX3v3r2iwYMHs2Aw0ZEUZJ0\/fz4WL17MBr1q1SoEBQXhupYWfCjhn5dLyeX\/sZ0ZJkz6KSmNC9x9GSYdYCvBpzfK8p5DZ9IWPA+dDv35sShIOobXJ3dATsOgymEBWeh\/VdDJIyUzhvQcJc00bdqUubDkzt64cQPOldn9ZH3oOjiAcnad4wH3zuI0CpooSsihExi+Uc8QaQrUz96HbXp6klS4j91wamtHUk0uv+hNFpRb26LoVhwqRICcAFDt4vBeuvdzgP6h3Au9DwU1BCdOnBDt3r0b3377LbOBy8rK2H8CnSIeCgoKsLKywtWrV5nE56irI\/j1a4TeLkDpnUIWIVKyUEK5nBy0Cn7BeisrGJaWSlYI5Qp+yRMPhbf\/QF7lwQCBshqE49ZA0YjCjLUXCjyQ9TJnzhwQr\/4xhZ5x8eJFDB06lDXnBzPq1atX7ZEU0Ni6das4ckSOBBnwlJvRsmVLFlaiDKYtW7awUFTXrl2RmJjIuPErVxrA0LAMbdrUZylkXCGHiNTSzp072UqhVUD1yST7WicePgbIv6INI7xIPRDoKSkpbLYIpBcvXjB106ZNGwYm5W7o6elhxQpbWFvrwdc3leXpcUVNTY1lqFJQgJmCnp6sfmho6Fc78fBXAPgxfQj4G6muri5evXrFYnkTJkxg6cB\/\/vknk3yyixs16oAVZrYoKzPCxdaXwF9CSkpKIJJHFuhfwwH5GDD+qjYCzmQkV5WCp8QPLF26FA8fPmRj6NChA2bMmMEk3cJuCmb0J0lOQbt5czG8vRAD2zRg9Qh02ozJweLUCz2Tcv7+Ab3qdDLniHgCYt5IzdCuTim9BCCpCP5GOmzjMAxWHwzt\/I7IVLvE8mEOXD8APb0iqKqqMimnCaNCpiYRU+RJ\/i+BTtYJqVTKyeQfkSH1TQJNGDOdzhH6BBb9VACdODt48KDEZKRre\/fuxW9Gf8I75RS84k0g36sX1re9irY5gzD+URycHsZhgFor3PYIYaAT0FT4jOBftXz\/jv1UAV36twFI0knyuTRfvi1O90a5u2PBiBEse9UvtTl+HaUDn8e74GXiJEkz5l5a2hYmW57+6lJkebL0ow91LR\/Tnn5Qgv744yfj4PTp00yIyJ\/hjrtwR3L4bCKXa89JOkdDkwVI\/BblvDBqV5oG4DJmKcRFhRLdqUMyCTkQyZrp1KkT20itLC2JrodR2UncL7sPZWVxCI+KNOjcjz3UBThZvy8gEGf11al8THvpH5Og8ZMxwR1U4Mc\/Cex27dox\/4WvSjj1QgcHVq9eDS8vL8keR44kA10WDUCWy\/Xr19nLtW3bltneNYG+X0MD3nEizLFqiIAWKrWC\/v9R0gkkfr47f8b5J0roOl\/SKU5KKoULUnPqRRKY5tvptKHWpl5IRxM5Rl4q32SUl28Ke3urWkGvk4j+jSpx0XsC\/X2Rfi5UR97pyZMn2REZTtJJYDljhYEu7ZEOGzaMmYfEtVAhnU4eKXVc25l6mmVSO2Tt1KRe\/kZ41mkofNCpAT\/Rlr7T+SHS+dwBhZp0OtUlmoVO11GA\/INAp8Y1\/XqEkZER6I9fPgepVCd0\/p9VqkYD1KZePvTd\/gFdNmIyN1JZfLqsM+7vm4R\/QK8BdM5k5McgP1fa2T+g1wA6l4LBGf1cugDfRKSm\/PP59J0LMHNBDLrGnyyqT27v\/5pHWicagKwXLiTHmTZ01pGOaxCPIh22437VgbKlaNfmeBXqTFot\/S9yLzWpXD4N8H8JFSc2DJRxaAAAAABJRU5ErkJggg==","height":56,"width":93}}
%---
%[output:554728ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:402ae825]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:9e326ede]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:5810a729]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:9ec141c0]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsB0_homo1\nTaille de la serie    : 198\nStatistique T_max     : 5.9089\np-valeur (bootstrap)  : 0.4720\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:17dd0a36]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:7233dc42]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsG0_homo1\nTaille de la serie    : 198\nStatistique T_max     : 4.8267\np-valeur (bootstrap)  : 0.6080\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsR0_homo1\nTaille de la serie    : 198\nStatistique T_max     : 4.9913\np-valeur (bootstrap)  : 0.5690\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:9c097411]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAGeFJREFUeF7tXAdYFFfXfheWLiAgvYhGsaARo0bUqCCRJPbEgn3BhiXqhw0UUVBRUIm9BWmKRkUjkcSI0QhBxYYlRuwGpa0gRRAWpP2eu979lrICppj\/+XJ5eNiduTNz550z7z3nPeciEBeUVBlrq0FRS0kBHB0B+uvqCoSFKewq20H94+IAa2vg99\/r79+YHhW56RAHucBw8iZkh8yFyfyDUNY3b8wp3nlfQVVVVdWbRuHnB\/j6Nh5ADjwdu3z5n3efHHT9MSshDhgGlBaA34CaTS+Yeh+HkqbOn3fBv+BMArL0o1ezcD+rCEEj27JLnLqdg1HBv6KyQB3PQ3qwbWfOAA4ODR8BWToBT+3Nj7Xh5+Q9M9cNh+TStwoPFBpYwHz1hX\/sGyCYd+h2VVhiJtx6mDLQnxaWQhT2G9aNaIPLicqYPEIDllZVePJY0Ch0Xlw4gnZOXZBWbI1jvhEYvFxU63jqkxU0gm3XG+UHvZHLqvURfzWafTeZd6DadmbtgZ8zCzf1PCoDl78F\/3TKEfyU\/Kzqx9+y2U0R6GTl3t89wLFZdti8VgWrVylhwqxi7NmqyfqUV1Sh\/8cCxvHmllX46XQl266ipAShsvTB0M1nrPwU3o9\/xt6jhghyXoq538yoZnm8j9GM3eyYrB1TYOZzQtaHPxDNHi7VQOfACtQ0UfLbaaBKen35JmxmBXP\/82+09LS0NNDvu2iM0+dH3akGevDZNISKbLF6pRIDfaDrc3wfpsvo4v7DKqTvEQDWgKGkBDHIgtoLMYzufQ8j50kQGlqDAMuNmIcQnZtYsVYH83qHwm+tLprYD5fdI+9D4AjUtZHpPwC6gz1YHw6simV7VBYX1LJ0fpK8qBXIO\/TfCaOut6UuUAnshQsX4uLFi+8CczBOX3tC6mJwS+egR+4VYPoUZfT4tAjnf9RCixZSLwbcg3kM2B17AbFYDc3LzjHe17B1YKA\/j9kA7ydnEBGpgq8G+cN1YkU1+uB9aOKjRqBr2jmzPkQrTXqORFnqLZSmJisE\/W0Ru3DhAsaOHYt169bB3Pzv93zqtPSpkcnIk5ShLNkUxbHtYGBSjiuJQjj6ASnhim\/Ve2YKVm2TWjqBPi4hDvEJQkSNmw2noYYNAl3F0hYvzkcxoMmSFYFeeu8iMv0\/QWXxc9mAGkIr1JmDvn\/\/ftjb21e7oRSkIA5xiEc8rGENEUTs75\/ZaoH+e24RnIKSsGTAe9DMu4elcyYjI2M\/TE3tpVZOzTrllSsT9\/qLK0AfX+8j33z7wl+ge34xeu07B+vmlTg7ro+MOvjgFdELAV6ceLDaPdbkdc73ygZWqMh5IuvLvxvNP1yNymoCVhfoBHYEIuD76qdmc4Urtku2w8vLi+0KCAiAhoaGwucgkUiwefNmTJ06Ffr6+rX61QKdvJch267hfpYEytl3oH1uHQO9tNReCnaYmxzgZAH0+ypySrEG3F75la+fhYV2GtIKLTC8eTjW9\/OE5fqrtSbStMU9ZaAp6RixPtTSvXui\/NkTkHsitOwAqw03ZQMnjs9c+QkqivKhat2pGvUQLb1MuQFlraYw9YlV6K\/XBbof\/BjgZNUOcEBf9MVjPJY9hAmSCdDz0mPjmDJlCmbMmIGMjAzZuNzd3eHp6cm+BwYGsrcpJCSEgZ6bm4vr16+jX79+uHHjBgQvSsuqtFSFsoOXbwzB3s3+7LvEajA0nsTIQFfqHQNjA3eo31Bn+zP2Z6DUvpR91gvUQ9NdTdnnwsLByMkJgLnGUxzr\/x68Uyxxp0CJ7dvcTx+frf2OffYZ0QsHUqTb2+pUIuLID6j4biWK9Vti\/ne32QCp+U4dgYmL17LP2Q9vwXXUMNzOq2Df16xZAxcXF9y\/fx9uookYYFCAcxkv2X55IOTNjYPeo0cPFBcXwzfEF3b6dgzwMzgjoxMCa\/zk8bhz4w6K+hbBEY4w1jFmoMfHx8PQ0BCLFy9mp7ayssKWLVuQnJws22ZmZoYdO3awNyMnJwempqa4evWqdCLlMgAbuJsbBg8ezJ7awEEBuJ38NV6U7MUz8QcwMPACtIqQc+I0NG9kQ2+LHsRhYijnKMNsrBkyvqrEB2t\/QmqqG\/LyZmPudEcUSnwhuRWH1X4++FmsxAYWFhaGtISj2Ld7OzYePo30nAK4jhmJAU594LMpVGYpwVs34vC8odj7SAnh+w\/CwMAAk0QT0KEyFfNcnJHwfRQ8riiDuLmpsByuEyegcytLLGhXgUttJmBbaCS7VuvWreUxl3E6B71NSBsE6geCaCRM5iVILTYmJgamYaaIRjTauLWBUxcnBnpkZCSioqLYtQnwoUOHwsnJCcOGDcOYMWPY9UQiEaZNm4aBAwciPz+fvR07d+5ENRng4MGDDJRpPhswJ6ocJZECmJmNhY3zBkzo3R5eXlIwi4pcIGweB5MmXyL7o48AsQnMft0jtf6M\/TDWfx9PMjUhFN5H9w9FGG9RDJcRX+CZ3Rj2UNevX48OSGeT7V4NZ3wdKp2dRX07YHnYMXazu3btQsc272FVi0y0XxXLqOn8ye8xfvocrG+Vgc66pZBUChCUYsCOHWmUD99HJnAyV8KMbsZ44rSM9a05WZL10pvx8OHDag+izKwMc8LmYGHrhWw79Zs8eTKbaI09jTFJMgndvLrBHvYM9HHjxqGwsJD1FQgE6NOnD+N5ZWVl9gYQnZDxqqio4MSJEygqKsI333zDKKcW6IFBG\/HCfi68+g7FzHE3GOjTfLfCwUafuVk1W1bWGhQXD2FvgbZ2jGx32fMtyC4yYMfXbEQJpY+S4Bt8GNMmuWKKmwiioc74sGMbBjqnigzxU3Zop06d2GB\/+uknbNm0EatbpsOiNBUldYA+3jQPQ9oZIO3zbXWCTuerSS91Wbo86Hc87yBCEgEHLwd0REeMHz8ey5Ytw927d6vdWps2bZhFi8Vi2XZ6QEQ5HPQ6LZ04qtLRE0+TnVF17SkDbf6anXD+oAWc3ZyRPTsbRS5Fr2REByD+tRjTPIVNstYOKQg4GMA4TROWuJ\/+Fdq1m4PpbYow3sNH5lHQ7L5gyjiU\/n6V0UtJlZCBbt+rF6MXcgcpQtVbGA3vwM3sFSd+7tu3L3vw3HrpPNyj4JPb7NmzmSW\/yS2sCfqkkEkYoj+EAUWcThMpB93AygChAaFs3yivUTCGlNNXrVqFly9fyiZL\/nb4+fmxyZMsnQyFrJ\/6V7P0mpw+ZrwIqu\/1gON767Bj2zo0bboLoz03wmdif3T16ors4mzke+ZD\/ZQ6dCN0kb0uG8J0IeN35TBl3La4Dc95c3AtLh1n7+yFQ08vtFKKw4aDJ\/BDwhUZpx\/ZH4GE775ByK6dOJ14BT4bgzFNNBbTP3fGonmzod6mFwLXB7FQnSiJwOzfv7\/slac5Rx5Y4nver6Ggy1OPG9wQDinNkRfTHM1xJPAIrsdcZ\/PWMAxDplsmunTpUovTbWxsZOMiw+CcTudasWIFo5fz58+zz+Tx1PJeth08iaDF09nF8\/PdZaCvdh\/Cnv7w4fNx5cogtl9Pzw7Z2dko8YqFleRnKB9KZdtp1jZU90f06b7w+2Acruin4beH0n3ce3mUUwzXsS7IzMmHkUoZjKxawtrGFnOb3EDB3QvwfmCMO8VSL+mzHnZYv3sfsxpugdyz4cBxJ6Axll6T7+WBZ1wtEcDAywDaMdrsLaOmo6Mj89NpDuTeC\/E399+5MVB\/msi\/\/VaqiPL9tfR0LniNa9oFc2cIoWXwEtEJhWilYQDS1in4EYmkf3mjoCkiAggPlyY5SALmenpjJeFqJPmGL8xf9x+A0nvnQDq60YxgiINGgQQ0NZvubzxNfREpBUkULJGvzv32ho6rIf2quYx0gLy0+6m9JjJSlbBpRzmuXRSyZIQ82DUvQOC7uUkfyoGlR7C9h2LZtiGDU9SHA65sYAG9QR7IjfKDsccBPP9hI4qvn6w3kfEm0P\/IuBp6bC16oQN5EqMoth1eJps0Ku1GwDs5lGJasxkobOcOH5\/asm1DB6eoH1chtR1FeLZrWq1uQn1zmK+5qFDafeegk7Qrz03cPaMZ2CsiARuWv4eX6UZob5OLW7eMZDfIPQeaNLhodPToNUya9Au6a1kgS2KOy\/cta8m2fxRwfrw4YCgkN2KhO3QRSh9eQTNREDL9P4Vq6x71qpKKtBcSuhQ18mjkhS8eS1B\/msN4EEbbqXFJgOYg8mgWLVrEAirSYwT37t2rWr16NYKCgpirQwfRDDvEbS5G+D+DapoFspKMgDIhrKwqER+vhPHmZbiWl4fKR49wBmCgX7+ej88+E0MsbosxHQ\/BVPkq+np7YNAAjWqy7Z8FOll7uu\/HqMi8I8uRNlRPr0vaJVWxLrGLj5f2EcdTI7+bJnLupaSnpzMPZcKECfjhhx9YqE\/72rdvzwKoiIgIjB49msUZgwYNqi4D0AlpQHu+OYh0yyG4FClEa0cjzO9vgWEDtaCkU4LVhzLRY9FYWBWIUWJigqZHj0JYKIbyte\/wvpsISs2scPKrQ8g6PA+B+Uk4dkTrLwGdxvom6fdND5eD3r1793eSyKjlvZClP1JrBz09PRwI1EHXORXYOaY1Vi5bils5rrhqlI\/rIVLQi4xNcOuHH6GZcR2ae9ywqcQTAYeWQ3I5Gqnb3eGRfhmnYk1lWSECgkJ\/nrGX90DktXBF2+WBZJbO1cgaCNenq8u7dO8ikVHNeyFuv\/ggC4WW\/TC25XO4TdOqBvopb2\/cEergdJfPEJEiQoq1Ne6Y9WK3LLRLh2v+aizZvQmVeRm4Na87plxNwuXvxSy61Or8CfJjgph7x0GXt1SSZdUs27NEh6LtHNs\/Ajh\/m7mkUVci48+iQEXnkYFOFm5tbY37FWYI9veApM1g5F2eC+2R19DPRhfmtyMQNTMbmQHTYGiYDIjiUFB4CUtmARNzBdhZDHzbxBwlZdEIDtbG3SNfYH\/0Quzv6widwQugrKHNxsBdOvosn6Lj2Saj\/xxE1kYXWepOPq1Xs56lZo6U32RjLP2PgM4DMnldnZ+PDLh58+a1MlM0RkYvBDhFXDQh0olICnhi\/gme350Gm09TsW+KLlYu240mh20hsvaDQbNipFUAfnmq8MjOg0tREco1gVXDP4XDpB3w9DSGqup+9Ow59tXE\/N8MC4FUE3SejOaZJBPvk8jeMVmWaZLPMMlXcslXE+R9v4HlVDVs7PFkvh10HEUwmLheocHK08sfBT02NhZffvkluxZhR5Omt7c3jh07phj0xMTEKppZSZ6kVqGqBaFQFe4zZ+FmejdEn3iKZnmzoCQoZ\/vXZGWhv0SCycbGuKEuDdPd8\/PZ33MGBrhZLu0nFHaEh0c4Zs7Ux541i7A9\/BusbiWGiQawXWs4jidcYv0WdlTGJ6oPoKTdDEoqajDxOY3sba54mJmLRYnFyCpTYf0id25GT+dBMr2FvjePD8Av743Hlq3bEPRZC9xv0hbb9xyCS9sm2HKjmB3n4+PDxCcSzbg7fO\/ePZli+meDTsAPGDAAZOkpKSlMoqZG16Hrk0DHLJ3r6ORrJj7Kx4oFM+G54D\/4fNgQ2HX2QouWo7F3d28cdGyBmPJyhL2WLt1MTDA7T2rpgXp6uKCri5CkJBzz\/BKrTt6DrsQDP4Y\/h\/eaDVBv3R1LnFpia8geXM8FQvcdxAXfUfBIkiYhLGIWoOJZKsxWncejxT2x\/q4qND8YgNnqVxB2OR3J2u8jJCwcHLANXcrQoTINPz7TQqRYD\/6txLj9Qg1fpRphYftyTPj6jEyhpPNzQYw0Eq5WcjBqJqcbyul10QvPVnHQyV\/nbwAJcZRxYhHp90ePyNQ\/At3fcxa8F3pg+MiR6OeeidwrBTDWdWf+u1lZ2RtB941Owtq1gFg8Gba2o9G9+wfYsUOqEnaT3MDszYdgZ6CE5WExuL3MGUuuCzHRczV6\/boFFTmpMFt5Ho+W9MS6O6r4Ob0CTlbqWNC6CK033WQ5T\/nAppuNJYJnDERkqgY2jbVH4uFgRGbqYZe\/FzqNmFUtV0lA8qREQ0CnikDKtZPERLFIXY3AlKcXChj9\/f1Zxoh8dc7p8onq4OBgqZ8edzy6TtD7OTmhb79JUK7MRTP9CMzN6I8tFRUKQY9X18WgpUlo314fFy9Kg6xbt+xQUhLCIraC2O2YtS0aWaXVS\/RcjPNZtgelL2DicwpZW12R9CQHHhekVEVt0CdOCPxqMwtK5DV1+beUbpSnAylFJ58gbizoLeRAV1R4XBfovAqAAiFqPH\/LuZ6U0jfSS+uW1iwnSUmMgd2d4e6wAukVh3BULM2YcHrpXmSI2Xod8JuuENHRIejUSZ9Z5IIFC1BYaAgHBysmaz5LiML0Jf7oPWwcFs37T7WqLkUTaerJPZi+aBluFyjBf+kiaJ5Yw\/Ki61tnwE6nFMeztZh103yRXKSGyCxDhO\/Zi7ZdP1IIOqXaeIWXIk7nlk73+SbQScNX5L0cP34cCQkJMpmAYebmBkFxcXEVkXtSUhKzxpyil5gocoP7zC\/RtkMnzBLNgLfPLAwf\/BnG9+yJrNJStKgainR0Q17FOljl9cFLww9gjjAU4Qnj9JfhM5H5S5RME6f0nLPSXVY4FFHWmT0QSjpfWD4cc37OZZzOc6bkMt4LGAGvKxX4aMBwzOrbCnG718j6UQLazdUVX7pPQbebOxCYVIxbL9T\/C3qmHtZ2E6DX5kSsD46UlULIWzrx7JEjRxjwf2QibSj31+wn6NWrV5WDgwPi4uIY6NS4y1hu3QdGsRfwUhLMtqv17o3ShARU7q+EoJMATl018bNEgsGVlfhuKiC40hlJr2s9qD+93rt3x+DA3G5oaazHAh8ulJE3Qc3L0RrTQn5mpXTywdGd61ewMD5XZkU160rIK6Dkx4e6xbj0XLNRll4z6\/S2E+lbg17fogDSx+vT0RVdnGTekKA0TC7oUGf5Gy8c5YmIuuSBmoX+LBpd3B3luekK77m+4IgOfJfybr0rMai4Pz7+7VZTUCaJkh6NWUzQEOvherqStgGqSgpB1WFU+8iVRxUzG5h6SQuaFLV3Cjr3XsiR5xqwfMRGfLxzpwuryG1M41mkxh7XkGswQczvY5TnP4XR1G14GjwTpvOjUHwjFoWnvoZSU7NqiwXqOuc7Bd1nhX9VZNhu0NIjikptbW2hq6uL5cuXM\/eMJpvycjPo6BytN4nBb6S83AItW74qxvRt8adbOQeQrJpC\/qoXz6Bu64iSW2dQWQUoCYCaBad\/Fehvrb1ERUVVkW9LxY0E9Jw5c3D27Fns3buXTYRCoZCFtSEhD1FYOA9xcQJkZWWxErGbN28iPDwcH374IQoKCtibIhItf+WjayAuLhz79k1ltYJ\/5YoHSfIvyH+9MECgrgWD8eugYtGu3heGEg9kUHPnzgXp6m\/T6ByXLl3C559\/zg6XT2acO3cOxsbGLJFRy3shegnfvpFtJ9BCQ0PZL1Uwbdu2jfmYPXv2fFV2cQW2tgOxbl172NiowcfnISsh442yTgMHXoZEMgqxsS+xY4cnPv74Y+aSvasVD28D5N9xjExl5KCTpT958gQeHh6sQJIaWQKBbmJigmWPJuCjXuXo9aAc\/v6tcPGiVLjS0tLCsmWxiIycwb6T70\/9Kex9F4mCvwO8t71GNdBJT6fX4vbt23jw4AHj+A4dOmD+\/PmsYIaB3jsUfZs3x+nQ3YAoBSfumsLEpARqamrYunUrA7gm6O8iAHlbQP6O42Sg0yRIEm\/Xrl0xatQoTJo0CYcPH2YFkfv27WMSpe1AW7jbrJEuDgDVMcZh1b4EDOqgzUAneqIAi\/QFsnSSMqnm71\/Qqz9KJnhNnzia1Vhv2rSJBQ1Lly7Fo0ePWE8CjrQKynJ\/sekLDNGRFlryNjV9KqakT0GykRFitmxBwHRpSR5VGJAwRZHk\/xLoVHJBlEo0Lb9ERj6TJKAkhnwJNLmOtOKMMh\/knbi6ujJ38sCBA4i1j2VlCm0DE6HS2gbiTz7F+7k6mPj7RIj6TMTA7Gwkv\/YECGhq8org3\/Hq\/lOvUQ10nsTgwRFZOrmQiYmJbPzvv\/8+q+EgidL\/7Fn4pgBBuVdZ9aq\/igp+\/DAZHXO\/xU39mziDMFZmzFvNAIQCJtlisXrQqSuKpei4oe1tjqfoWb5skMZPy1dIKSQjevz4saxglC+7kS86oreaiom4pXMZmjxA0reo5oXOIVjpv7oqdHewbH3OyZMnmf9KtdXE01RbvXLlSjg7O8v0CppYyTcPDbXEzh2UjvKFxYOPcN\/iJ6i\/TuERODVB5\/\/coSHA1fX\/BF5nFBtyeJ3\/j6C+42v+8wgaP2V65LM\/lP8k2iCwO3fujGvXrlWjEk4vtHBg7dq1LPbhcxwZLwP9+Olfqn69fJ7dCJ2cQCcdnHz1+kBXVVWFXadoICUOw7duQeSqVm8E\/f+jpRNI8vXu8k+crJ0aL5eWt3SeJ+UPidMLA51XA3DQa9ILLYYaMmRItQtTntHOzg4EuleJPVsxPStbCwEBGm8EvUEm+g\/qRFhw0OUz\/XWtIeWpOopOT506xZbIcEvn5YqUJlQIOi1IokwPNXL9KCKlp13XmvrMTOmE6esbz\/x5kg0Ucfo\/CM8GDUUedDpAvtCWvtP6IeJ8vkBBEadT3969ezPphBLkdVo68RVfL0Ogy1fmKvrvERYWFqBf+fYulbwGofqOOtUCnSYCstglS5awIclX9DZ2jP+CXjdiLDiSF7zkvQ76\/EcCm39BVwB6fem6xlr3v\/RSP2Kyuhfu9sivxKiZQaLJlLe6VmLU7E\/FNv9rEWnDZIBrv1Xt2rS+1koMWq5BGkxd3C6\/rJDTT11zwf+i9qLIzuVlgP8D8s8VGjhM84MAAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:15e5ec7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:3af32782]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_med_h","value":"1.0421"}}
%---
%[output:1bf2e1f7]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_STN_h","value":"0.2887"}}
%---
%[output:4534507f]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean_h","value":"1.0902"}}
%---
%[output:783a1690]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean","value":"1.0597"}}
%---
%[output:6101957b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_median","value":"1.0893"}}
%---
%[output:0788e32b]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAABSdJREFUeF7tnD1MFEEUgB9WYGH0gg0\/0kisTNSGn9pSqUgQKhMKEi2BQC6Whp8QiiPRSEGwgtCKlR2FASohtpoQgUpzWgmF8cwbM+e7cefn7c3sLne7Ccne7dy+N9+8fe\/N7BtaKpVKBfIjUQItOfREeQthNdC3trbg6OgIZmZmxMW9vT0YGxsT5\/Pz8zAyMpK8hg0osQp9cXERVldXYWJiQkAvl8swOTkJxWJRdHtubg6Wl5ehUCg0IIZkuySgo4X39PTAzs6OkI7Q0cpxINbW1qCtrQ1mZ2dhdHQU+vv7k9WwAaXVuBeETKFvbm7CwsKC+A6hDw4OVl3MyckJ4F9+\/CPQ1dUF+Gc7YkFH2NPT07C\/v2+7f1Nd7+vrg6WlJSt4I3Sde5EBFgV0dnY2HNgnb77Cl\/I53Ci0wsuh6079QwMslUqwsbFhdcFa6KZAKqG7CHDSOGON7jzfrUI\/eDbgpB2HiRa6mjJSwBwBThpnrNHDFx\/g+Pu50IpCl993X2uF7ad3a7TmMIk1OeIIyBjPutQxPQEcJjl0xjDongC8BV7b\/fgJBm7f\/O8pUEXk0BnQTU05cSCH7gn6veK2CL6Yp9uCb1NCp24CmUcFRu5YNKVPN2UWKkDpCuT3mI\/brNM2CE0J3eRT1QGRn9Ed4JFDt5mU5npUZqHClYClOzFlI1w1mtLSVUgI9P3nH1VLxhOflq3Ku\/DQXQMddRsIAWeRaMV4SODUffi07GDQ5YsNFNDR0QHr6+vQ29tbfaMUYu2FWqgp0EW1k4ClVaOvpuc+shSd2\/Fi6WdnZ9oXFxwBnKwCO0QDouygXPGj0Gg76jbkObVwzsSF68tlew4TbZ5OVxnRuunBEeDa4aiMAiGjy1BBRrkQnZVjKhjSrXiFTl9K443lu1M850B37bCaO9ssDl0HHQzaXrqVenNvmw5xDdFpRipdjXxdx4FuUlwNhDT4cTqstvWRd3Plc5g4QUcF1PenWJpRbyClrgddhg\/oaVg59+nXQseRw+oAWY4xPj4uzrEagDOqLpaucxNca0sLuDfo0rqxFqYen66CU3NwmdJxAUe1T8OteA2kJghxLV2XW\/sATlNEX\/fj3IfDxNmnx43U8nehgefQFRPxCVxNFeln1b1wJ2Ycy1bbZsrSfQKX1kwzHQQtJ1F4nebmrhOzemCn5tN1FuUTuGrRLiuGrhOzCwldtaiodWwfHctaxhI3znkJpNSifE1yTINk8uOhBtd231R8ekjrlpMe2nGdH7fBCXU9FejcBStO59Oc9LjqmTh0nwFTZihqZ5NcMXQFHcSn6\/Ycye+vPirB78vt2iVWm\/Kqb6btsw45SJ5uK5UeenkAv9pv2biK8gb1JUSWfbO1Q5oGXtyLac8RXnvw+lhYOT0u\/fwWqZIK\/u3j7rh9y+zvTk9Pxe4Ul+Vu49Kubs8Rbn9B6OqSLEK\/8u7vdshmPGJtf1EDg2mjV9QSLf7+1f2WZuQt+hxro5cKPd\/SGMZ+nKoBUHS+edffABiXAWjK6BIg\/KnV2Hdir73ocvekMemqz5LQI6oQi8OFBT0r\/y\/AVH0WGjoywJf0h4eH1fSQy4UF3ZS7h+4svb+p+iykHjjYKysrMDw8DFNTUzXVEZykgw3dlEaG7LCaWcl\/iYLf0+qzJHSQ1k5LUjhcLiR0ClatPmtI6JzHKAkAKINWnyUhM8rSOVxYls4NGKEAmKrPQslUYwqteONyYUFHwVnJ3WnKmLZP53JhQ0\/CkhpdRg49hRH+AyL\/9aRJ0sDYAAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:54d2849a]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"fit_test","rows":2,"type":"double","value":[["0.1130"],["0.9029"]]}}
%---
%[output:06f55808]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BaB0_S11S12\nTaille de la serie    : 179\nStatistique T_max     : 8.3169\np-valeur (bootstrap)  : 0.1920\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:94d2eb33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:627fa8c7]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BaG0_S11S12\nTaille de la serie    : 179\nStatistique T_max     : 7.9204\np-valeur (bootstrap)  : 0.1940\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaR0_S11S12\nTaille de la serie    : 179\nStatistique T_max     : 9.0756\np-valeur (bootstrap)  : 0.1260\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:1c2a28bc]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAGd1JREFUeF7tXAdYVMcW\/hdpUqVJEWVVxBpL1NgFNMaEh5jEgmIBbKhRbEFBNICKHbuoUVokMUaJRtGIUUBNbLFhImpsIF2kSJW6zzPrrJfdpcT4Ur6X8dvP3Xtn5s7958x\/zpxzBpFEIpGgjpKUBNjbA\/S\/mxsQFlZXbek9qh8fD4jFwKNH9dd\/EzWqctOQGeQMk8mbkR0yB2YL9qORYbM30fUb70NUF+g7z6RgW5Aabh01UwpgEpLYgMQQKwxsYPdsnLtmAs8OAfD3F8Fg1GcKdYouRuFJ0Eh23WB0AKtTXVKAjEAHlP32E7uuatwCzQLPKwAoX09ElTX0YOZ9GLn7lv4zQf8lvQiOq27j8bae7OXj4gCxnRTkCEQg\/uU\/jiQB7wY3+MEPJHVRHkvhfDCU3X489y1YLD1RAziqk778fTSdsYfVebJjCqtDhV\/XsOmlVMo44I2MLGE2\/2tZnSQPK1TnPq5VMmubwDcuyvV0KCoqq5Boq6vWqFZcXom5X9+F2hNj7NheCrhFwM5PCrN8IbC5xHOpD050R8etIRh4\/CGSHzfCgXGz8f4sO+j0HiFrTlKeGzGfSbFIU5dJt\/6weVAztJRNQG30wKlEnkLY9TUfgfjSfNGhvy+9ZBY8l5jqatTAkmiFyvnsX7G3x0xALJVwZYVRi3sYkCSGOM5dNjEH9nRC9JlriIhUwwbHQLhNrKpBMQT6s6MbYe57nHVLoGt1fY99z\/vGT\/YoTjvyz87cMAbPE07A3DcGtCJKr0QjY+2H0Ow0CJUZd1H5VFHihZKempoK+vwVRYHTiVa2nErGiLGZ+EDtXTamFhIrPBYlvxyfPwAChSbCHeDSHx8HuwA7Bnw4wlndz3yfY9lKDcwfEIols1MbBLqQ+0ly03z7wtB1Q41VwoHKO7CsQRMkDyyB7eXlhUuXLv0VmENkMC9W0qapFo580hUk8STl3ud+xLMlLtIBBfhDu9Acxes9XgBsByTFAUkBUsyTrQDbM4CdFGRxkgSPxMDAgk44p3cLCHdjq4AkfapvuwbRi5CCOHfTClCmiF8XsYsXL8LFxQXr1q1Ds2Z\/voWjIOnE54MurMJl289e2H12gH0cEPfCBrSLBwJe2Iv+bq\/oJuml1cLvJ7lBIg5DaPFGTNaeL8VEJMH+kZMwYtfyBinSgtMhrBmBXPbbJWQEjYL5ggOMQhpSWJvAoagueSarLq9AOehfffUVevfuXaNb0k+ku87gDLPKXOGq1DpryFhqq6MAOkl6YIeRSDW5DI8Da3By0U6kf5WOst5lTOphGy+dAF4IePq8vCZuKWG2vL0dM+Jgqf0It09fZVIurwCFJmPTBQdZHXlTsDZOl69Hz1I1b4fKjDsy81OehvgzlIFOYJNV5v\/in3whqyy4NBje3t7s1pQpUzBjxgykp6fLqnp4eGDRokXsd2lpKbZs2YKpU6fC0NAQubm5uHHjBgYNGoSEhASI5BUpPbwlWrLZ3XdxH1uGMtAbMr2il3utRy3ZinALeIQwP0U7viFd1VZHmclI15I\/aY1GusawXH0JKlp6NZoLFfflm4nsvYSSHoAX+wn4s\/e2gx1sYYtkJMsmYULpBBh4G8hAP3PmDGbNmiUDOTAwEK6urmjTpg3WrFkDmtiQkBAGOn1PTk6Gs7Mztm3bBgVJ375\/O4J8glhnTiOm4EjUHhnoKkcHwDT0PjQTNNl94WQYrDFAk11N2PXCwmHIebwVksYS5OTmYPLkyWyGqQhflAa3a9cudr1Lly6yQZJkCNusWrWKDZgK3ZvkOgE3b91mv\/m9O1d+hLurK4a10sRliNl9ofQJV9nPv6Uw0Pv06YOSkhL4h\/ijq2FXBngc4mR0Qs8aP3k87iTcQbFtMexhD1M9UybpBLqJiQl8fHzYOFq0aIGtW7ciOjoau3fvZtd0dHQQHByM2bNnIz8\/H9bW1rh\/\/35N0O\/duwd3d3ckDktE3qI8tBniicoHRxm45V3KYdRjNFCsg5wTp6GVkA2DrQbIDMtEo5xGsHCxYPWqjKpg5tAFeSvuwMTJBE7eTmwAq1evxpEjR9jAwsLCkJOTg3379rHrZE3Qc4cNG8aWqFBSfvjhB1kbIyMjNhnEw9NaleHsgRDMu9KITaReSSYmzfRENxsxNuyPqfEskr66QG8b0hZrDNewzV0YXvk5aBxHjx6FeZg5DuMw2rq3xeDugxnoNF6aFO5Fad68OQYPHsyukXQTndD7qKmp4cSJE0hMTMQ333yDmJiYmqDv37+fveC9z5PxpEMqNPyHw+KLm6\/AHCtG3p0VKC52hqpVPMx0ZiG7f3\/ALJPVk5d+1blz0P37KHg4DYT7ui\/BJ3X9+vUyBSaUdi6Z\/JpQ+qnv8yejMX66J9a3SUdXvTI8rxYhKMmIPXdU03z4PzTD4GYqWHbwPLg085XFQVeftAsu0+fiwYMHNeinwqICnmGe8GrjJVtRfIJNF5liUukk9PTuid7ozUBfsWIFysvLZauTxnz37l32EXI91SXAi4uLmZAR5dTYkRLoGzZtxv2ZlsidcAAaFzWYBM8InY+lmp+y7\/LlyZNVKPk2GEbHEqF7VFd2O31DNTB3EyxUXloxgoZECVRoaRLQpHD4C5Kk88nhg+fgC6WepJcUFlduk0Y5wWPaFIxvmg0Hk2LcKNDAp\/csZBNEz+NWjDy9KJN0TnG0qu4suoOI0gjYedvhLbyF8ePH47PPPmMAC0vHjh2ZCWpubs4knQBu3LgxmyQO+s6dO2sqUgKdgCjs54Wne2fKQO87dTU+GFkBb3dv5I7RRcnVhVJzkgpZLv4BgGs4U5za+7XR1KcpqvWNkZG4A+3be2J622KMn7dUZqcLwSJ6od9C0PmL8Hq0xGlybG1tayhAYT\/coiD+JP6vyyzk9zinTwqZBCdDKQ0Sp5Mi5aAbtTBC6GqpD2m092iYQsrp8pLOx0x9jx07VjYXy5YtY\/Ry\/vx50HcSpBqKlCRs7HhXqLfuA\/vlHbHj4KdMOap83gTBlvsw9eM0ZIkvo2zLLqg9T4POdzooGFcAjWgbqK6dgoLqERCFuMHo+K9ocaYnriXshZ2tH6wRi2VTR+KsakcZP3\/77bcyDU8SzKXe09NTJr1CvicwhwwZUmNyhMAS3xPP\/h7QhUrdHa920mTFWMEKUWuicOPoDaa3PsSHyHDPQPfu3ZkeIsOAW0A2NjY1xiXUSSTpfDVSO\/qtxHo5iSCf6Wym8j3yGehpQWooP+KDFr5hyL1yDbrf6UI1TRWZOzNR0aGC1bXcZwnRopbIzl4LXesN6FU9FNGJk7FtYzqOR9rjdq60Hn9RIYVYWFgwS4AsAKHkc4uHFBIfsLxlI9\/f64JOYxMCT79FpSIYeRsx2qRVRkVPT082Fs4MdF04Ri4MdJ2MBhIwKrWCfup2Dny\/u49xTbpjzvttpMskwA8IIz\/Lq8JdubRjY1UQwEwttrmwj4PlVWukFloi9kQxbM4OZR5E4Ra\/Rmev8UO4OdKw6YemM3YjM2g0cxXXt3utb0dKmyTar5Ctzu321xhirU0UJP1RbjEGB11FfmkVa\/TsbhWaHZ2Cxyo0jJfexiQxAgKkfVJEiYqdHRDuao8kcTwm3vXCF+3WwqpZGe7G3qjXVft7X0i4OTJwnIfcAwEwnfc1nh3bhJIbJ5nnUn5zJHxGXaD\/3rG8Tn2lboC4uzn44U4eik\/boHLZXKg2qYKLrT78ksIY2BSGc3WV\/s+Le1IAwiNeaCJyE7iHodOzp1jc2Qu9TOKh3deZgfKmitDmrnyaKgNd8ryQhex07V3xdNc0hcfJWy9CTuc+l9rGSMpVWYTsdd5JwQ1AnXCuytMfjuK1yRCpA+VW12G+eCfzMGpmZmL+lvnQcdTB\/Yr7CE8OR5p1mtSSAdBo0CmsMw1EbGMfVD+xw3yrpygYtQkfuc14nTEqbUP+9KqcVBiO8kP+sY0wdg1CRuD7ULfqivLkG7W6g6kzZZJO7mji9NoKbZpo88SLcH9BOom4m7sAqA73w5AOCggIwMKFCxEZGcnMYwXQScGtXLkSQUFBGHboOzwqykHW8Z6otngk43XLMiukqWSy50vUyl7yjFjqfQx3Q98kVZTvmIQMFX3EnSyDKGwKorL1MXXzQbZbe1NF6DCjPsk5pjd4MpP2ugLTdYEugYTtkIUBDtr+k64ijqdCmx1S8tw0TEtLY2bhhAkTcOzYMVy7do3d69ChAwoLCxEREYExY8aArDRHR0dF64UDQi7ellE+GGtphR883sHt0qZAXARFodGsKpac6TApNUPku8nQzspkQm4f9sLbGG+H8TZ34Pt1CdpvTcDT1GSUbBuLc1kStJq2UcGV+qYmQNgPBTfKUhJrxE+F9+sCPSU15X8e4FCaDXDgWiY8IhNRqH0TzQyfouf5MXBZfA8fd\/uYxf7dwpPgSvzNAhev3OsPB7TAJ\/\/pgBOLTiBiSyu4znmA77\/LQ8cL\/2GgNxo6V+a4+l+ATX3yaFNd4Tp5FwG14\/Ry4eKF\/3mAQwH0XzMK4bTtOkLdOsHLvRBXpvij74BqPLL8AdtGFePtXwqYVJ9h1oqUcbh7PVdPDwtvjcKIghbYk7ca2fP2oUXr97DSrBdKUm9DHZUMa60+zkwKhWafMNBQ23U+UQ0Btq6cl7oknYOuLMDxpgRFIRtgWOBR3PwqAKVth6FwqATPh0Zg0q09EHXZh5Vt96HpC1coL4etNLDDSRU+vxTjUoIBzuhrouiqNeKwD6PutodK8CCo5QdjlXlvpBeWodH0r9DLdrCsvZAGSDFqNO\/AIka1Xa\/vpakPnb6jUJFyq0bslLery3qRl\/Q\/Arq874ieL+xPQZF+EfMzAuZ7YPjEd\/HEKwuHR4+DZuRkhJUHo08zFxgXFbF3uKeqCnczM8zOy8PO6GL4BQBnXMWwdQuDOMkOLQPCsbaHO6IC4hA5ax0SE2+jf\/AVmSKVj3\/yIEPTufvxZJMzywygCRAGH+qyvTm1UApGVWkBTGdFIC96I5sE2pQJJ7UuSde8qAlzF\/MaINU32fL3CXRy4fIgB\/0mZerr6\/vKDcBdumT2fHnqGsI3r0Th8HToHTJCZsEeGBtPh9pLCV\/15AmGkIPK1BQJmtJgxhhJPo47AJKfGyM3QwxDQ1voGt9Aq83G+PadmZjd+QPE5elgpXUmLAy0sV3dEcfiz7O2\/lNHYqLPWgYu5cGY+Z5E9o7JeNptLD7ZtE\/mJo3cuQV933NUMPeEY79wcDeCv\/gGzu10sDVBuiLn9G+OhwZdmV+8vaEawvd\/iwfZhQqRI26nX7p4CTEuMW8cdALewcGBeVCZpMcfPyxzRF14mI+AhR7outAGuxx3o0fvhei8qAl+gSb0f96MyjPa8NHLRLN04FM1M9j0zEPGtGLkHTCA0XF9jPC6ChWV37Anchb6r3CGzqlZKCyUxhbJ97De+xNcPBOLkMhvcGXTTHjG5rIX7IQ0GegpWydhQ6opVI2bY8nHvbFpmQ9uaXVEaMRe\/PbbbwywjT2q8JbkMY5nayMyw4BNaGKxBjY8boqFPfUw3MkJK08\/ZJMbsWYptLJu4pPgw\/jQeRzshjoqgM6l9U3sVpXRC48VsBgp5TIKpSX1WSncpo1F0Tt9UK01Dro5xnhfWw3nro1E2ZMSaEi0MHS41I8cH2OGj7XyMPGaCVwH6uNxiQ4O7zkCjeipcIm5AN3Gb8G1nytCT0pDVn0NyjDdexm66Ffgsy3heBDmDa\/YLDg0eQrHVo2hatISpp8ewuO1I7Dqx0zEZqlgkEk5vPubo9Xy06AdZ\/T0vixatGWQIRy2xOLA0e9lAkP2MQVh9mxaDd1D3tirboefYk9imcFVNl6\/4v7o966DzEU8Z84c9OpVM8uAbG7KiamN0+1fZvzQNjCuFt6RpxdyQfMYqqWlpSLo+3\/OROimaXj2dj+gkROaXN6C6txSmJ+QIOPcTOivCkCbiv8gDc2QLfoSrbIHwqL4bVSIN+Gnrmq4GhKDAsMC2C\/+EFrf2WKBZzvs3rYRm0f3QNHFA1j8yBJPnkszBXghKXB+GAx189aw8I9Fuv8gXL71APOvNpLVIS\/ePIsM3JMYwePzGAa63RQfHEuuVABduDvkAWLqiPvsx40bV68tXhvoLQWg15aQrAx0nh1A\/SrQCw1ulPsoDJw9EC5WLmwZOkzdiu0ntSGe4InKz4uRvT6NgWHiZYKp5XmYmFuMdfoGSGjZAUdCjmBBwVbEzTWFcWokmjc3gZmoCJ7qP+H58zL4PjBDv6HD4bNsNdKXDUVVXjqqch+z9AkVrSYwXXAQWUEjZdeLDNvC52o1budWwqt9BczVyzE\/QRvrrdPRWa8MJ57WpJfIJyYI\/2Iv2vXoXyPWKgSdtui1pdVR1tfmzZvrlXTqry7QybcvDNvxSSRWEZWUlEjIyX716lXmP6BCDVrMbgHPtr6Y6OoOj5mzYGPhCA+vUVDVS4GoUBsqz\/NZ3fxp+SiaUATdL3XR4WIH9Pd3wf7r2dhiNBuhodIsAB6xJ1MwKPwgbqmKsTt4Ky76jVDK6fc2usHnpjqjgk9srXF6w3xGKZv6aqC16wp4BGzE5D5WGKiShM1FXXAt4SajlPh1s7A3WQ2rrLNgWZaC3akGSCjURKB1FhvrkqTmGPDReHj7raiFGJT7ZWqt\/Jo3RP369ZPY2dkhPj6+BujFs0phUxgKw3uHEHPwC9Z9jx4DcOXKOWT5iNC\/rTWqd1XjwYUHaN2nNa43vw7dm8aQGKxDfOAw5oHk0XS+3An03CsxDKjo76Vp0Q21XqY4D8eE56eYv3zToXiWutFUrQLv6JXgcoGWTJE2VNJrw+tNKNL65kJ0M61Q8uGOG8grrsCQ9oYIde2IorIqDDsZAYljJH7WOKvQB\/nQI+KTEC+W+gLIh26XROEMN+ZXJ5u45MJ+do9nVdF3Ap37u+k3z9QV2uO12ekmM0KQGfge8x42tumNFK+3od6iM1R0jZgtTtcee9pAUvZq8yY\/8Ibkp\/8poHPrJSkpibkjswrL4LbmIO5+tZSNuc2qNohxjlE6edy2pUwoOgxQXxGCThud37sj5Rscnb6jkbqoF4zc1qMqL0Pm3OJ+dpoIzdbdFbKEi84fqNUJxsf+p4C+evVqCS1VbkfOCr+EwzsDUdppNBuH5t290P4iH+NaOsjCVwQ2JVhSoiWBTQ7+hhR50OVDbjziU9t1oc9FpGOMJg6zodVlqCzJtCQhBoWnPodKEwuFQwG1HSSQH\/efAvrx02clN3+W7g5J0umh5HBPSUmBSCSCvr4+42ab92ywInUFoqKioJaqhhkOMxBgHyDLPaEdHxVhKpswQCtMjWvIBNVXRwgiZfrygwREZyTRRG88+ZTTHXe01dX3mwC9Xt8L0QuBykE\/ePAgy9s+dOgQDAwM4OTkBG1tbWzfvh1LlixhE1NZWcnqEJCxsbEsOZKuFxQUMFuYHPaULrFgwQL2Oysri0VNZs6cCV3dVwlJ9QH7R+6XJp5F\/ssTHSJNbRiNXwc1y\/b1dsk3R8o2TvU2flmB+rh8+TI++ugjdkUY5FBXV5dujoSg00zzHENqwLOT5s2bh7Vr12LgwIHo1q0bmxBhoYiQlZUViziJxWL2ffny5Sxl4a868dBQkP7sekpBX7p0qTS7VCRi+SjDhw9nudWU5nv27GN0774VQ4aooUcPqceRCq0GmgjypC1evJgliFJ9ymz6q048\/NlgNvR5CqBTYuXo0aMZj587d45lmlL6b3Z2NgPxu27PUDTCGJb3feEYpYPp06UbDw0NDXaf8vwo\/49WDAf9j\/imG\/oi\/6R6SiWduPvhw4fsPSjfj5xVVHx2+eDH5T+i0pIiQGJA9AgnTtyFmdlzBjpxfadOnVhdAp1WDE3iv6DXFAmFcB0Hi4B\/5513WB5e3759mWIcfHYwfu3xK3zKd6FZRTNMbzsUAcnJeD8rC1paWizpXfOlj12Y\/vb\/BDq9Nx0KoJxMylvkhXwupOeIBRRADw0NZQqQ+JwKJb1TqjIpV\/p4394BP7EbWrVqBYcPPsCCX3\/F4JQUfN2uHb4PDoZOVJTM10J9UXbr\/xPotdFcnaCTpBOf89xqknTK4eA8TV5HmhRKrBneuTMcc65iWE4Ohncejjg8YnlQvMjbvOQ+4Gl49XEwuRPkC\/2Rh4aW12lP\/iJh1hqNf8eOHUy3keCQacyPu\/B9h\/yhBkom4pLOM3spGYn8W5TzQn3IghjcDcBpgZLaqbRt25Y9kB9YItDpge3atcMSyyW4oX8DqaqpLLvrXuU9WFu\/OvIuDzql5PkrHl5TiqOyv83xcvE1CPfXaU9j8xN4M2j8dLaIH1QQxjkJbDKdr1+\/XoNKOL2QQUEmtp+fH6MZLrwMdHk3ADtINWkSbt6UHmfp3LkziCaUgR5tFM2AN3s+BqpD9uLeD5UyTufKVHiK7Z8o6fxUnHB3zWddeKKErvFTJSTpFA8lSuHBaE4vDPQLFy5IaDapcDdAffRCkReS9IuZ7bCteyEsTzfH3sEpjIK4IlUGeoNE9G9USXgUUT6iLz9MHpKjXeipU6eY6cwlnQSWMKUcdxm91LUj5dYLHSlRdqY+O3sddHUPYs6crmzbr6paO738jfBs0FCEoFMD4SEA+k3nh4jz+eEFoaST9cI5neoOGDCAHfilEyP1ugGEoFPj2sJcFHClj7C8CedRg9D5h1VSCnpt9PJ73+1f0JUjpgA6KVKiCfKfUOFp06+T4vwv6A0EXagA6fsf2dj8C3odoMvfEioM4allZUEJ4VlPbjbxUwj\/jzvS13IDCE9icFOH8jcomkT+GHnaodMF8hsrHsTgf1Thj6yW36tH\/q71hW6A\/wIFb86ONiXNhQAAAABJRU5ErkJggg==","height":56,"width":93}}
%---
%[output:8eb7bf5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:298d53cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:2b5eedfc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:2fe376f1]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAACQtJREFUeF7tXHlIlVkUP89xea2MRWFWoNFmC0VlBW0QLVA0tG\/Qvmi0WpntpW3aapTRXkZEhRUU1QyVEa0SYdJqhUQyIST1RzkiEzr8DnPlvuu3vU2f9l0QfN+72\/e75\/7Oueee8xwVFRUVZJdqRcBhg16tePNgNujVj3nggP7z509KS0ujEydOMAxNmzalU6dOUZcuXTRhQd0jR45ofhccHEwdOnSg0aNH0+TJk6lBgwZc7+vXrzRnzhzKy8vTbNeoUSPq3r07zZ07l\/r160dBQUF+WZKAkfTCwkKaNWsWFRQUVL5oQkICLVq0iBwOR5WXNwJdrjxixAjasWMHAVAz0EU7LNr69etp2rRpmmN7uxIBA\/q1a9do6dKlLu\/TtWtXluaIiAiPQUfDzZs30\/Tp0y2DjjbNmjWj06dPU0xMjLcYV2kfEKCXlZXRhg0bKCsri9q0aUN9+\/alc+fO8WSPHz9OgwcPNgQdddFGlG\/fvtHZs2fpwIEDBNoaM2YMbd++nUpKSirppVu3bkxlTZo04Wbl5eX08eNH2rZtG929e5ef7dmzh9v6ugQE6G\/evKGZM2fSly9faOLEiSyV8+fPp8+fP9P48eNpy5YtFBYW5vLuMr2ooKMi+po9eza9evWKRo0aRampqVRaWqoLuuj8\/v37NGPGDP4IWpo0aZKvMQ8MRXrmzBmmABRI55AhQyolX2+bG4FuS7qJnPz48YOWL19Ot2\/fZmqBxdK6dWu6fv06LV68mFtDqUFq5WJVkaLN\/v37WdqtKlK0iY6OZvqJioqqe5L+7NkzphbwLawFcDush6KiIoqLi6MXL14wX2dkZFB4eHglAFZBB11h0Ro2bGgZ9JCQEKYW8LmW5eTtKtQop8MDcfDgQdq3b18VpSnb7ViEY8eO0aBBgyyBLuxtLOKAAQMq9YGZpLds2ZLrg9Pbt2\/vF8DxAjUKenFxMSvM58+fmwoPDjVJSUm8C1DMFKlWhzLoqvViOgEfVqhR0O\/du0fz5s1js86sqBxrg26GmMb36rHfShfikGNLuhW0NOrIx\/6OHTuypdCiRYsqNWVFO2zYMNq9ezcrRVvSPQBePvbLVovaFWzuhQsX0pMnT9hxhaN5z549bdDdxVw+9qOt3lFf9Hvy5EnaunUrfxQKFUd04WXUOpHaitTdVanj9U2tF9jSDx8+5OP506dP2c+N0yH8I8JPXccx8vnrmYL+4MED5tQ+ffrwhcDr168pMzOTBg4cSMnJyazU7OIeAoagCyUG6caxWAAM+3rJkiX8DJcEdnEPAUPQYTGARg4fPuzi0xZOKviiU1JSKDQ01L1Rf\/HahqCfP3+efSPw\/LVr164SKvA8pBwcL18E\/OJYWn59Q9BxAHn06BGbdPBrywXf3bhxg+1mHNHtYh0BU9BBMVrSfOHCBbZo1F0ghi7NT6fyf\/52mUlmWBwVOSKtzy6Aa0Y6gyg+2unRDP0CellhFpXkJrpM6K+QPyjNmezRJAO1UXyU0yPgPQbdiF4AOIB3oSNnMgH4ulT+iAillJj6br+SXxTpryLpABzAu1v8ZjKC0\/8tzqnC6Xm\/9XJ3jgFZv9fvwR5RC17G0uEILlf59GkfjryTA1M3wM2bN2nZsmXUv39\/mjBhAn348MF2A3iHOTkeP35sGp+OG\/krV67Q27dvqXHjxhyXguN\/\/fruKxEv51vtzVu1akX482VxREdHm4LuywFrW19w9O3atcunwFuS9NoGlK\/mm5OTw4FKVi9JrI5ryulWO6qL9XAanzp1qg16dS6uDXp1ov3\/WAL09PR0jrWEqfz9+3dC9AJuzxAfqUYTy9PMzs6m+Ph4QoCsHMpt04vBYgrQYaXBtY2w6ebNm7N39erVqww8skW0gEesO6LXYGKrOqHWgI7Y8tWrVxNCN+TiayUn9y1AR1wjXNgiIwT3CZBexOCIkBC5HaIdcLmDhUFgbK0EXcQgYosinlEUAQqie+XnvmIi0T9oBFaMXN69e8c5Uri2VBMHLl++zFHGuHXDSb5Wgg7fPW6x9Pz6et95C76RIs3NzWVQEYYtg47FWLBgAa1YsYJTa7SsH6\/oBWCsWbOm8t1kiXv\/\/j1LAqKxkHpSr149Tj8BRSBUDpcfKKizdu1aDoUWqYaq5JpdmLgDrjznyMhIHhv5SKAKWdmhTz3QEYe5d+9eunTpkksyGO6ON23aRE6nkzZu3EhYGJ+CDn86+FXcHGlRgJi0yN0RLyy2m1gY5BaJZ1r9yPW8yQNS56w1vhany\/QAPkdC2rp161iikS2C8G08v3jxIu\/Go0ePcgaH3qI5Uv8sqEga7t4dp5gsBpS3ltYgeFE8B+euXLmSzSzBv0b9oK58FSj6lkFxJ8ZcbyxVEIxAB7C3bt2ixMRE9j\/JnteXL19y5giS0kQ2oC7o4QnZFUnDo8gd4PW2u5BSBCWJxTAKxBdAqFtbDyABiBGtYZEhZapyAwDqQhpRiPpd7969mU7A4SNHjuTENDj\/UBAfhBzYtm3bMl2JxAVD0KfERlDGFOtJqupLq5yqUoCorz4HuKtWraKdO3e6hHhoLZ4Wb8tmpKAAPdD1BEVv4WXQEdEGPXTo0CHOj4JtDh2lWlFGugX6Q+xcByQdgAN4q8XImlD7kHlTHhj1rEj60KFDOVJX3j3yGOquAOiIsbxz5w4rZqGUvZH0sWPHss0NWoHiF5Is5oE0Hri91ZKfn8\/JwOD\/zp07U6dOnXh3eMXperQgngtJ\/PTpEwctQaJRhDVjldOFXtAyGdWFk5UlwgFFaDX+B2CqHrLC6ci2A3CwRFTAjQRVl148\/b0X1RIQAMug6lkr4sUFYGgjtp54JitcQTeoJwOvZelgXihCWWMOOJLjs571Io+vpUghnePGjdMMlsWNWmxsrCb2PgcdoxgpNL3ToszvPXr0YOmDNAJMmI4oemahVu6oetpTOV0GXZ2zsNNxRlD1iszpRtJsZML6BXSrOkCvnpES87RvVdLVz2q\/AEb8zoz4cQZVQfrav+PVidRTYEQ7f4EuQgExjhGn6\/l0bNDdXFnVepG3v9YByxN6cHNKVarXqKR7O3l\/t7dvjvyNsEb\/\/gL9P60GPcImSVFLAAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:8ddf09ed]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAByFJREFUeF7tnFlIVVsYx78jmqgplFcoJzIoE4y69KCgFUj0EAQJNiCUmUMPDQ54LSEhbTQbrphQlJUSaiD1INlLENEoIRWidSkkshuB00tekRt6+S\/uOiyXa5+9z869PSf3gh46e42\/9e1vWPtbuqampqbIKbYScDnQbeXNBnOg28\/cf6D\/+PGDampqqLGxkWGKjIykmzdvUnJyshIb6l69elX5LDAwkBITE2nbtm20a9cuCgsLY\/VGRkYoLy+P3r59q2wXHh5Oa9eupfz8fEpLS6OAgABTW+Y3kj4wMEC5ubnU39\/vXmhJSQkdPHiQXC7XjMV7gi5W3rJlC505c4YAVA86b4dNO3bsGO3evVs5tt5O+A30jo4OKioqmrae1atXM2lesmSJaehoePz4cdqzZ49h6GgTFRVFt27doqSkJD3GM577BfSJiQmqrKyk9vZ2Wr58OaWmplJLSwtbzPXr1ykjI8MjdNRFG15GR0fp9u3bVF9fT1BbmZmZdPr0aRobG3OrlzVr1jBVtnjxYtZscnKSPn36RKdOnaJHjx6x3y5cuMDaelv8Avq7d+9o7969NDg4SDt27GBSWVhYSF+\/fqWsrCw6ceIEBQcHT1u7qF5k6KiIvvbt20e9vb20detWOnv2LI2Pj2tC550\/efKEcnJy2H+hlnbu3Oktc\/8wpM3NzUwFoEA6N23a5JZ8rdfcE3RH0nXk5Pv371RaWkoPHz5kqgUeS1xcHN2\/f58OHTrEWsOoQWrFYtSQok1dXR2TdqOGFG0SEhKY+lm2bNmvJ+nd3d1MtUDfwluAbof38O3bN9q\/fz\/19PQwfd3Q0ECLFi1yAzAKHeoKm7Zw4ULD0IOCgphqgT5XeU56u+DTOh0nFJcvX6ZLly7NMJqi345NuHbtGm3cuNEQdO5vYxPXr1\/vtgd6kh4TE8PqQ6evXLnSFHBM0KehDw0NMYP55s0bPeFhBvDIkSPsLUDRM6SqDkXosveiOwEvKvg09MePH1NBQQFz6\/SKrGMd6HrEFM\/lsN9IFzzIcSTdCC1FHTHsX7VqFfMUli5dOqOmaGg3b95M58+fZ0bRkXQT4MWwX\/Ra5K7gcx84cIBevnzJDq4Qmq9bt86B7i1zMexHW61Qn\/d748YNOnnyJPsvN6gI0fkpoyoidQypt7vi5\/Ut917gaz979oyF769evWLn4IgecX7Cz7H9nKHX07cc+tOnT5nOTUlJYR8M+vr6qKmpiTZs2EBVVVXM6M23Yil0buQg3QibOWD434cPH2a\/4SPCfCuWQodHATVy5cqVaWfe\/BALZ9XV1dW0YMGCecXdUuhtbW3s7AQngytWrHCDhZ6HlEPHix8K5gt5S6EjQHn+\/Dlz+XDuLRY86+zsZH41Qvj5VCyHDhWjkuY7d+4wj0Z+Czj88b\/+pMl\/\/vbZvQgIjaGQxGJT8\/NJ6BMD7TT2+g9TC7KzUUhikSnwcwbdk3oBcID39RIcl0Vhv9d6PU1LoZs1pP4i6QAO8N4WS6H\/jMsInf7vUJe367GtftBvKaZUCyZoKXQeHOFIVow+neDI4lTpBw8eUHFxMaWnp9P27dvp48ePzjHAixcvLM9Pxxf7e\/fu0fv37ykiIoLlrSD8Dw0NtU0dmB0oNjaW8G82iyshIcFy6LM5Ybv7wkFdbW3trIJ32SHpdoOarfG6urpYIpLRjyBGx7XUkBqdhK\/Wg\/eVnZ3tQLdzgxzodtL+fywH+q8EXfReVAaD77a4Zk8pZ2K+SXR0tOYp4hww9HpI2yWd5\/UheV4+fuVg5aR4\/I5EfSTYh4SEkN7xrdcUbG5gO3QA1DoLx9rl56oJ4mbD0aNHCRKP5E5\/K7ZC\/\/DhA7vJhqR7resdqHP37l2W4IlvnZBqnCrKHyy0ftfaANSvqKhwP0YOurhh\/G3C1UJce0Hhamx4eJi5eLzI6lJOhdZTf7ZCN6MWZNXCF+5NX+gD6XRcnXFI8fHxbpXFVRvfDP42oR2\/O6RSbVyQUIdvIt9greDHMujckIoD66kWlZRqQcfEy8rKdA0qh4IEUPEmnLxw1dxUGyv3p5of3zCsh9shcW2WQVf9bQAt6KorJdyYGoGOBeHqIPIM+VVBvkitzeHSjkQlqDrVOHrQ4W1p2RZP6s9W6EZUghEYACr25Qm61phGxtGDjqsqSCzFGyQbdJ+BbsSQyjCMGFIYOkg6UpkvXrzoNoDIiZn3kq5yCWU9LkNXQZNdRnkzRVXx5csXpcek0uliLCC\/TTypSaXTZRfYp3Q6FiJ6BbJ151ZfdLl4\/c+fP7vdRvm1B4jy8nI6d+4cy\/gChNbW1mmeiRHvxQx0n\/Je9P7IDp8sFsqLJ\/\/W0zEA+hINqQydS63op+tFvUYlHfVkP13vBp2thtSqyNEIdKvGNtOvA90MtZ9s40D\/SYBmmv8S0M0sfC7bONDngL5V0P8DJ0LHszwhVNQAAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:747a1fce]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAB15JREFUeF7tnHtIVUkcx39XNNfVYrMNLFNUtkwoUmIzcDWSMAiCEullD3spYuUWiAWFaVmJhhvlH2IvRUxD9I8oCISIFioktoiyJXHLmyFktn\/kurKhy3dgLuNxzst7vPcoZ0BIz8zvzHzmd36POb+Ta2xsbIyc5lMCLge6T3mzmznQfc\/c3tC\/fftGFRUVdO3aNYZm3rx5dOPGDVq2bJkUFfrW1tZKrwUGBlJ8fDxt2rSJtm3bRqGhoazf4OAg7d+\/n168eCEdN3v2bEpMTKQDBw5QSkoKBQQEeL1NttZ0t9tNe\/fupZ6eHs9Cjx49SocOHSKXyzVh8VrQxc4bNmyg8+fPE4DqQefjsGknT56kXbt2Se9tZidsDf3OnTtUWFg4bj3Lly9n2hwRETFp6Bh4+vRp2r17t2HoGDN\/\/ny6efMmJSQkmGE8oa9toY+MjNCpU6eotbWV4uLiaPXq1dTU1MQWcPXqVUpPT9eEjr4Yw9uXL1+osbGRLl++TDBbmzdvpnPnztHQ0JDHvKxYsYKZsvDwcDZsdHSU3r17R+Xl5fTgwQP2t4sXL7Kx3jTbQu\/q6qKcnBz69OkTbdmyhWllbm4uffz4kbKysujMmTMUHBw8bu2ieVFCR0fI2rdvH7169Yo2btxIFy5coOHhYVXoXPijR49oz5497FeYpa1bt3rD3L6OtKGhgZkANGjnunXrPJqv9phrQXc0XUdPvn79SseOHaOOjg5mWhCxREVF0d27d+nw4cNsNJwatFZsRh0pxly6dIlpu1FHijGxsbHM\/MTExMw8TX\/27BkzLbC3iBZg2xE99Pf3U15eHr18+ZLZ65qaGpo7d64HgFHoMFfYtLCwMMPQg4KCmGmBPZdFTmZ2wXY2HacSV65coerq6glOU4zbsQl1dXW0Zs0aQ9B5vI1NTE1N9fgDPU2PjIxk\/WHTlyxZ4jVwTNZ20AcGBpjDfP78ua7yIKkpLi5mTwGaniOVCRShK6MX3QlMsoPtoD98+JAOHjzIwjq9prSxDnQ9YpLryrTfiAie5DiaboSWpI+Y9i9dupRFCgsWLJjQU3S0GRkZVFVVxZyio+mTAC+m\/WLUohSFmLugoICePHnCDq6Qmq9cudKBbpa5mPZjrFqqz+Vev36dzp49y37lDhUpOj9llGWkjiM1uyszqL\/l0QuyyZKSElq1apXXZxQziPO4pVgKHbYWj3x7e7slB0MOdA0CsMeIr5Emv3\/\/nvW04jTOga5BoKWlhU6cOMFep+3cuZNu375NR44cccyLCjNLzEtbWxs7596xYwd9\/vyZvWLDaaC3586Ophsk8PbtWwe6DitLNF28hxnow3\/+RqP\/9BncTvt1C\/g+kkLifzU9Mb9BH3G30tAfRaYnbLcBIfGFpsH7DTqAA\/x0b8FRWRSaVGlqGX6DPlM0HcAB3kzzG3RMEjb9v4GnZuZrq75BPyabNi1YgF+h24qgDyfjQPchbH4rBv3Dhw+EHytaX18fK+LJzMyktWvXWiHS7zIWLVpE+LGqudxu91hRURE9fTp9batVMNTkJCcnU2VlpWXgXY8fPx5D+g6hKDdw2ngCUEYUJhl9KWKEnwe6lUKN3Hi69MErQSillXwc6Dq770D3w+PBoaM64c2bN+z4GjWUqCLmX3OYndaUazpKkY8fP0540y82tceVn82LfVG\/iEouWRMrtGQy+QEcjp6VzUhFF4eOaoP8\/Hx6\/fo11dfXU1paGpWWlrLSD7NtSqFzICj2FKHxhShhom4F18TCfC4DCxP\/zheKTUIpNcqno6OjWbgaEhLi4aB16im7nxIgnyvuzUNgvCXDSxq8HcOnNGbblEIHkObmZlVY4jUOB4VD4hcUWJDaNf4ULVy4kJUvAz7KqhcvXmwIupFjaJlN56Xc+GKjrKyMZs2aZYr7lEOXgZDN0AgANS2EWUFFLepf8OWc+MZKS67WRvN7yaCjshha3tnZKVUovR0wBb3i\/l\/UO\/gvRYd\/R8XrY\/VkezQU9lTvRbVo+\/mnKaKZkN1MaR5k5kLPvECumr\/ANbXoBfe6d+8eqy5DIauZZhj6rc5+KrjV5ZFdvD7GEHg+aXFSag5MrVZc5iB5X1GzZYC0HCnmpOWktaBzX6I0Z0bgG4YO4ADP2\/afI6hmu7lP+5SRidaClbCUGyVbNN8I0aFqaTqfj9ZTqKbpPoGu1HQAB\/jJNNGUGMn0OBy+SWphKJ8LHCvXQC3osk3S8huig\/eJecFkYNN\/7\/6bfvnpB13TInv8xQUpYWARsP3KkA9jOGT8G9fxSblaaq6UqwVdKVfmQ\/zuSM1qtVYcrIwctNJtMTSE09OSy\/v29vayyEKrDkctjxDXOe1CRrXERm2xvKhfaXJEyACC0FCZcKmBQtquVvwEuciUtZzhtEuOOAjZp4ZqtlwW7YgO14gDE201\/vcKpO\/eHgMkJSWxb6G6u7vtfQxg1hzZsb\/swCs7O5sdes2ZM2dSUzYcMk5K+gwY5Bzt+mETHegOdD8Q8MMtHU2fIdD\/BzllXrPqswh5AAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:77017bdb]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 198\nStatistique T_max     : 16.1830\np-valeur (bootstrap)  : 0.0020\nPoint de rupture      : 2023-12-01 00:00:00  (indice 174)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 173\nStatistique T_max     : 13.6541\np-valeur (bootstrap)  : 0.0330\nPoint de rupture      : 2013-01-01 00:00:00  (indice 44)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 24\nStatistique T_max     : 3.3092\np-valeur (bootstrap)  : 0.6280\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 43\nStatistique T_max     : 2.3770\np-valeur (bootstrap)  : 0.8690\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg\nTaille de la serie    : 129\nStatistique T_max     : 3.6938\np-valeur (bootstrap)  : 0.7570\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:5d36622a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:69cbfed6]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expS_br\nTaille de la serie    : 198\nStatistique T_max     : 54.4894\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2023-12-01 00:00:00  (indice 174)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_br\nTaille de la serie    : 173\nStatistique T_max     : 55.7647\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2013-01-01 00:00:00  (indice 44)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_br\nTaille de la serie    : 24\nStatistique T_max     : 6.0197\np-valeur (bootstrap)  : 0.2520\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_br\nTaille de la serie    : 43\nStatistique T_max     : 7.9682\np-valeur (bootstrap)  : 0.1460\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_br\nTaille de la serie    : 129\nStatistique T_max     : 6.0524\np-valeur (bootstrap)  : 0.3830\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_gr\nTaille de la serie    : 198\nStatistique T_max     : 64.4408\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2023-12-01 00:00:00  (indice 174)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_gr\nTaille de la serie    : 173\nStatistique T_max     : 68.7796\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2013-01-01 00:00:00  (indice 44)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_gr\nTaille de la serie    : 24\nStatistique T_max     : 9.5710\np-valeur (bootstrap)  : 0.0600\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_gr\nTaille de la serie    : 43\nStatistique T_max     : 11.5585\np-valeur (bootstrap)  : 0.0280\nPoint de rupture      : 2010-08-01 00:00:00  (indice 15)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_gr\nTaille de la serie    : 129\nStatistique T_max     : 8.5233\np-valeur (bootstrap)  : 0.1710\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_gr\nTaille de la serie    : 14\nStatistique T_max     : 6.5172\np-valeur (bootstrap)  : 0.1860\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_gr\nTaille de la serie    : 28\nStatistique T_max     : 3.0664\np-valeur (bootstrap)  : 0.7240\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:2dac2ad6]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAGbJJREFUeF7tXAdYVMfafldqFIgUYUGQNdhiA2NDE8OiURML0XgBRQys3aioeBVsAbwW7DUaYwGMHY3GXuOiRjFWjEHEtigCioJRAQFhf75Z5nhYdpEYNbn\/zTwPD7tz5syZ884373xtVtJ50Vn1nH\/VRxMHM1SmKAF48oYqQKYEZCmAzBmQyQDIgXBRR0c1VW+8FGXdRcY8X9QYsAiZq0dBOnYzDKxqvvHnvsoDJGq1Wq3rRqUSiIsrveIBqGT0pwKBDtCXUvTpv6gQ8IGBwLWvUrC+hjO7slMZg+7ygHKPeRq\/Dffn\/YvVW\/pEwNL7axTnPkb69C7IT\/6Z1Rva1ELN6SfLAajdTkKNTSwgDd2BrI1T\/jtBj4gAwqMBRKkAeUyJvBLcGsg1hcQ6EFCFaURfC3xHizSk\/ujAxNz9+hmcsHIoAxxJZtp\/PoXtsFWst\/vLB8Jhyn72mdeb1GutU5A44AbWjpAGbxLaqIY4ozjrtl7h0zeBryKtf+YeScbjZ2o7c5NyfUQoVQiXE8gl6GsjWla2GfhHEQBPiYxJ+Y1LmTh+voam1S3N\/OxNPIbPGn4s3ElSnhUTzKRYYmrOpPvd7mNgZOUoTIA+euBUok0hrH5WT9DStQ\/Z\/t9FLyS0tUFirvgDExoFeAYiXA6M7bwNnXzfw6nbzYAozYLoo0rABplrGdB\/37UA9pP2sjoCvapbJ\/Y5e0uY0I7TjvZAMub3xrOE\/bCfdAC0IvLO7kb67B4wbdwez9Ov4vmD8hIvlvTU1FTQ319RynF6TsFz1H6YhEz7Jmw8bX4bjs5blyLcg++IxCVimhENW6lGeJwGdAK073El4lwMGfCVBZ14nReS3LuT2sIqYD7M3HuVwyc7dmqlJkj7RgJ73LhxOH369F+BOcqBHpxwHwtcSfo0Uu4acgEJs90AroZERwMpKYAzqS7E9+IJKJF0ZRTGmmqo43qXc+gQZMMoplZ+LlJMqr6UXsTgcu6mFSCejD+LVHx8PPz8\/DBnzhzUrPn2NZxynE4MHs52Rg2Y3hEpiA2vVQq6SkPvbNMs1WAIdJmgRCIQaqwSbZL1p7ohZZ0JPArzoTR6sXfo20gfH1nNnksg5yefRvo8b9iPjWUUUpnC7pneGcW5vwvNtTdQDvqGDRvg7u5eplsVSENTIg5xkEGGAASw\/6+zlJN0T\/bQ2uwZJvEmcPBzQNr4TsgfelUHrdBgtNQWhEOGMPSZloicjbH4ITAAqeNkTHefXKpL8w1QrDLajt3KKERbFdTH6drtaLyG9g3wPD1JUD85UJyG+DN0gU5gxyCmZJxiK0PTQyACsSxvGUJDQ9n3gQMHYtiwYUhLSxPmYsiQIQgJCWHf8\/LysHjxYgwaNAhWVlbIysrCxYsX0b59eyQkJKCcpOsEfUMa8t3z\/8Bkc3WSLCYyjWQIVwJhr8lK0qUyUl3KcBcYmNvAMfI0qlS1KDNemmC+cf9yKZHRi1jSI0q0NAKcpFoOOTzggRSkCJPQL68fLEMtBdDj4uIwYsQIAeTp06cjICAAdevWxaxZs0ATu3r1agY6fU5JSYGvry+WLl1antObb96M7AkTWGePhjxC9RXVkVYKepWsKrAbYAfTBFN2ndczOphlydpSedL9CR5GPoT6HTWqZH0DO+9DML2VwK6JX5QGt2LFClbv6uoqDJIkY8CAAUwqqMycOZMNmApd6x\/QD5d+u1LmWtLZE1AEBKD7e6b4BTJ2XSx9YjXzTPIdBnqbNm2Qm5uL8NXhcLNyY4AfxVGBTuhZ\/gP8kZSQhByPHHjCE3YWdkzSCfQaNWpgQilWtWrVwpIlS7B7926sXLmSjc3MzAzLli3DyJEj8ejRI9SpUwfXr18vC\/q1a9fQLVCBTK+7yA5RCUCmdeqEgi9Gw3pPO9YZAVp1Zzgsl6xERlQGDB4aaGhoQxqKbvwK6RwvZE+4g1yvXFiHdgUOv4f1CyKRnb2TDSwqKgoPHz7Exo0bERkZyVQ3hUKB7t27syUqlpRDhw4J91hbW7PJIB4e\/F4+jsWuxpizBmwiLXIz0P+rIDSrJ8P8zQewc+eLZ5H0VQR6\/dX1MctqFqORKKbjagqNY9euXbCPsscO7EB9RX10aN6BgU7jpUnhBr2TkxM6dOjA6ki6iU7ofYyMjLB\/\/34kJiZiy5YtOHDgQFnQN2\/ejCmzFuLWZkc8r7tV4PSnDb5H9qKTkCoWI3tkNnJ8c2B4LRlSRT9kznUroY9oBrpG+psiX7qDba6G105AqqgDieF43FYqQJNKg507d66wgYmlnUsmrxNLP\/V98uBu+A8Nwty6aXCzyMezYgnmqazZc71tHyH8phQdalbB1K0nwaWZrywOunH\/FfAbOho3btwQwKUPhQ6FCIoKwri644QVxSfYLsQO\/fP6o2VoS7jDnYE+bdo0FBQUCKuTxnz16lX2J+Z6akuA5+TkMCEjyinD6QT6nDlLkDC3GM\/lpwTQm87fgP3S0XDwu19moPTl\/kx35HptgnWoNcx3mQvX0zZ8CSBSmAzxjUQXVGhpEtC04fAXJEnnk8MHz8EXSz1JL21YfHPr7+2FIYMHwt82E11q5ODiYxP8+5qDMEH0PK7FaNOLLknnFEerKikkCTF5MZCHytEETeDv74+vv\/6aASwujRo1Yiqovb09k3QC+J133mGTxEH\/9ttvy0s6AcG5WtBeiDasiyBVSAVJL4d+aUW1zdVgO8EWhkUuSPnhGGxH2cLFKwC7P24mGDhisIhe6LsYdN43b0dLnCbHw8OjzAYo7odrFMSfxP8VqYX8Guf0\/qv7w8vKiz2WOJ02Ug66dS1rrIlcw675hPrADhpO15Z0Pmbqu0+fPgI8U6dOZfRy8uRJ0GcSpDIqI0lYT+9+SOuTjOyQ7BecviENBa4FTJqpaDi9KiyXWCIjygAm5+vDckki4\/cixyJYj+wP\/FgND5IWw2aqDaoWdEL8Ry44ZthI4OcffvhB2OFJgrnUBwUFCdIr5nsCs2PHjmUmRwws8T1R1x8BXbypK6BANHN9kJ0SDmc4Y9usbbi46yJ7rx7ogXRFOpo3b872IdrkuQZUr169MuMS70kk6Xw10n30vZyeHh4ej7Vr\/djDX2gvnZDvvgJVHlWBxfcWQDFglGqE3\/3aoMDyRyAmGpamwYL2UmjlANycjdSsD0Eaj4t\/BzxPus765C8qphAHBwemCZAGIJZ8rr3QhsQHrK3ZaPf3qqDT2MTA03dJnkSgTVplVCwsLISxEB1z7UU8Ri4M1J6UBhIwKgLoT\/ML1dWMDYXlQFa+Ql5bY+KXlkDcQnSpwaQJSZRYodFHgZQ4ICAaiAkEwrhRIQMibmH00x5YOOdHpqMHH5mCiGrv6vSf6KOpl9WLjSOTeh\/CdthKZMzzYa7il1mvL7NIyUgiY4l0da63v2w8f+S6Ttdu7egIqAJfWGaBJUZDtMpZ4+WlotIYEFRkchWUYWRS8UkKR8N+QZhra4P6t4tZm6vVmqPT\/B+Zq1VX4ILaiJ1X3HLUVy82jiy7jUFWbATsxmzC73sWIvfiQea51DaOxKBUBPofAe9V20q0JZ0tMwUQHSWB0enOKC4uQM3lYaji4gFVgGfpClBBDTVUKsAzRgmVUgaEKQB5iQNMeRQjtyjgUkeJ0cGANGE+DkyfAyfvobDoMEBn4OL5g1TBh56XHC9YjoWpV3TWq589YaE5cifQvRx0Xm\/uGYAHKwaXw0RbexFzOve56AOShOx1+WAkhxIfqB9eOixwE1fPXMdfQqp\/DBAXBVDUTfAoppQwlTOgkmvcLszLWFoHiu+RL5wiTXxVBODQ1FA0\/6whjJwa6QxcFN75TZBQAo5HjnITDuisJ\/ogf3rRw1RYeYfh0Z4FsAmYh\/Tpn8LY2Q0FKRf1uoNpVLoknTZR4nR9hYwmMp54EdsXtCcRd3MXALXhfhjagyIiIjB+\/HisW7eOqceSASvj1UU\/f4t58+YxS4o6S0nLwK8un+BE0FkNf6tTAAmFgAh7OVQyGWSqaBY3ZRxfZKBpZ0AWK80EWXXkhVTCrOALKCM+ROOwVSjMShWkmO4UBy7y7ySy0JvYh06Toaueu3\/FVEX9kXOMVhNfBfoiTxWBTitYO8BB5j9pNMTxVMjYoU2eq4Z3795lamG\/fv2wZ88enD9\/nl1r2LAhnjx5gpiYGPTu3RukpXXr1g2Sz5edV68f0AR8M6UBfb3xNJLs1bgTuggGl21QlNkF8JxaGheVQ6ZSITBC408nqo8OlANKuSYlICYFCCPgKRAtZ\/tAvOM2tB7Qi\/G5vmjRq4CuTyppf+D96WpTEeh3Uu+88QCHpNe3F9QxisYC6MHTl2BfQWskKfYxUGseN0C7g99gk4uJhrNlEQiP0PiXlR4ekMfFQekRDqU8EIgOg7M6CiscJiJIJkVy\/e0M+KzYqcw\/ri8u+ir0og9wvlIqCtdpuwioL04vp+JPvfEAhyR4yxX1PO8G7B1I71z\/az4OObXDsyH9GWDGZ9Jh2GAxDA89xeMvxgLKKASqJPAoTc9IkclKQNcIeonKA3+n61CtaImuxWOw8MsmuJd1Ed73f8TGOi3YM3RlALzKRkphvJfFQfVNTEWSzkHXFeB4VW1F+z7JlnPpau8PpIzLZTIZ8qu\/hylhw1GnnxV+Gt4G5uczMO7sVHxdLweQRwARUQiHhKnl8SYmSDE0RJ53DhZWleF+mgXmLgzB9u7VMat+V4wfDmyT32IqZsb23kyt0xW4oEFxlbFK1XeFYHNF9boAoM3VrK03aOWIg9u8bUXai7akv1HQSWWMnLcIvxi2wNqh7lj\/y2UsCBuATgO7YE3gTFgeCUCT\/TNwbM5BBEbHITpFA7rfNEMopFKMzM5moKeUME54mAxHVUcR0s0Isxp+gH2r\/DDbwhVQBuJi7mK4dgl6XcKisx+eglGU9xh2I2KQvXsBmwTaeGlCTJwaMpqrSNJN401h72dfxu\/\/ugct+en4z2qFvx8kEpYjhWf1WsMw5zLqftULR\/2OwfTwPdQYYg5DdSa7ft99JkbU8sPVtXZIMNUEM9o6PcIdJ+Bajgue9hiBd0POoZHDXrju88Jsu3dRLagBGv48A+s3bYWjoyPzRZATi4o4QCF+OW1PI5c8bcCIErmP\/tTWlVi2dgt8G5hhSUIu627UR064aenGnve+lRGiN\/+AG5lPykWOuJ5+Ov40Dvgd+FOgix11XJ2kscyYMYNpiYKeLgz85iNM+s8gFH\/VDuk934XtiC145NQQ+eE98LH\/Qty8aIZen5zG8G\/AJL2PQTYeKHIQnW+J\/CsNYbx6JzoMOoTfroSi1eYOWFqnNZz8z6Cb6QPMXbWexQ55KCs5Obncy9PgtL2Quu5Z0KIITdS3sTezGtalW2JGnQwk5phg\/m1bjG9pgc+9vDDjyE3sUZ5EzKwpqHrvEoYv24Eevn0h79xN53Pp2a\/DWiVBoELeTq6nh4WFscAN+WEkw9ZfVnsYXH4R0ckpgM\/gwSjsIUXGxKaofkSF95b2w\/3MITBKK0KhhRXc+yagx3Zg5mMpTHyy4aXOwdmzljhX1Aye62Lhk\/gQ3wV9CvPh9bDNxR81g5dinMWvaFs9D5NuOuCjz77AxFmL8OBuCgI+74TPzO7C631rIWcx5+F9jO4lx5Hbz9BeWoyFW4\/AzMGF6fC7h7Zl0aLF7a3QZfFPiN21Txg76cckPKsWRsJ8eyi+N5bj558OYqrlOQZCWM5H+PCTLoKLeNSoUWjdumyWAenclBOjj9N5BiHpb5SVoquIDSe6LjaeaEIkPt9dVHetmojvln\/DrCoqn\/l\/iYzgAuR17IhaPvthkmqGJt9GYd+93rBckg+nYVfRcwew6VJNjHpSiFtD0hCdL0P+ITfc67AaNyKt8E1wb2zJiUd2lTA4\/rIAsZFhMJY1RYDvF8h4lFNmrOQrDzC6IHAubaqnTsVjxI4koR158cY4pOOa2hpDvjvAQJcPnIA9Kc\/LgS62Dvmqoo64z75v374v1cX1gU55EmSFEOgac7F8EUu6+CqNhWKrzDj6wiJJAD35fi7+PWIw7oR2xXPnvSxaZDOuIxp95oDY0zthueQxnIap0OoMcPiwFHZfNUVxw71IOOgK00O1Yfz1aqi6W6FH\/FCc+3ccii0aoJ1FLURFaYIVir6+aFx0G2FrdiJz+QCWSsd1eDKcbEdvxv2FvkL9nYNrMXTidFzJKsS49wthb1yA4IRqmFsnDU0t8rH\/QVl6WXe\/BqLXfo8GLT4qE2sVg04mur60Osr6WrRo0UslnfrTB7qY06kduVYoZWP58uWacN3g6HNq0\/MxOHfuHJP0zWcysHL6GGSNHY1n7RQsxpk9shnyOtaC3YAjMMh8hIwoSg5tBaniMLJH1kRO62RYTh4K06vxMD4QDonVBdzOSoKd\/wmYJqWV2SwjFF44c+UWorftRnxYLwT9lMVesDHuMr+MdNJBXFsQiAmXjBkVDPeogyPzgxmlLGxrApeAaRgSsQAD2jjj4yoqLHrqivMJlxilKOeMwPcpRphZ5x4c8+9gZaolEp6YYnqde0zgJquc0K6nP0LDpumR0dfD6Xo7L70gad6qjfrTju2hVCoFeqEITK2RI2HYvjt2r5mG6is06cgm9o2Qn\/4b0jbYotjmPKrPqg6zI2Z40qg9Ot4+jGvGBTh\/oA6KraYyh5jsO38Yxqqwd6\/GGUQ6+r3t8xlQu\/dp0qLDB\/0LX06YLVirBDqtgAfN+mD4wo1CkHeg7+fo9+ww85cv3K5kqRu2RoVoZZGLXx5XFTbSykq6PmBex0b6UtDpUMC9J\/kI2pSExb0bQJw27ZZbCMed1+HqXR8bDG6zHChNobAWsRr9aUJcL77L4fT4PhyOJ+JuLlA7cguOnavBjB+xP0Q7T5H7ZbTphdfXGLYaGdM7Me\/hO\/XccWfcBzCu1RRVzK2ZLk51t4PqQZ2vURV1lcrkp7810PUNkpy2BDPPBOFhiogSfzuduAgLe5HayLwC1CAGoOjT0aOlx2FKrU26rJ0EKp4IsfGir563MWvrg9SQ1rAOnIui7PQynkjyMNJEmLo0L\/M8mrynJ2PLHCLQ9d5vBXRdQQzxYAh48jIT8OKsOApgxMSUGJsENLnQyRcWAwQElLjYRQ0rSujkhwHoqAuF3HjERzsUx+vFziyJmQ2qdxmJqq6dhSRT8r8\/OfwdqlR3KHcoQN9BAm3g3wrodBJDuXdHmSCGm5sb8wFToQym53Z2cF+1Cg8aNxZI5XpREdt82zo4YLKjI6sXB2T1WZov47vKXheDSJm+3NdCoT6S6NxTm4VEUloh9L1qG9+\/h6QnJyeruXlKQQxKhNy3bx8z1bdt24aJEycyHHibx1ZWLC8vQqFgjnyuz5LlNXbs2HLtKVfwrzjxkJd4DI9KT3RITKvB2n8OjBzff+mccuNIl+H00psr2aBcCsbBgwcxZcoU0ENjY2OFLCWahE8++QQ3b95EixYtmJVlYEARI8DX1hfxpvH4vMfnWLNmDcvt4O1pUv6qEw+VxOCtNysHOiX7UOYpL2TFBQcHMxClUinLSO3fvz+aNm3KkiMzq2UipFUIMkwzMHzEcHa6gYq4\/V914uGto1nJB5YBncxXlUol+Caoj8aNG8Pb2xtnz54VQB89ejTLdHpo\/hCTHSfjV6tf2eMCwgOwPGR5OdDfpG+6ku\/5t2omgM6DGJS6RtzctWtXrFq1Ck+fPoWnpyceP37MzFnK4aPrrVq1grGxMT6t\/yl6Pd2Fb2xa4oMlH+BQ30MCvfD2\/4Beds5ZslH0soVMuilDlTiZrD0CavLkySzyTTl7SUlJ7D85p8SgZ2SYQip9BjfXHvhotSPWfKxJtqSNl7f\/XwKdFAqiYKJp2tt4IRZxdnZmGEv2HjmmHj7wxRFyUhHNzc2ZhHOVkeo2bdK4Anx8xmPMmJ5o2bIlRjZpgruDDNHN7CH2LbNF8ydPkFjqKiWgqWgfM\/lbrfO3OJgyoGv\/NgDp2qS98KR5U1NTDB06lM2cQhWN6JgUuG73wr133kfGRlOMvXwZHQzuwMbGBk2aNAG150Xb0CCDiv4qU8QGFm9Pv1dQ2fIq95OVzX5UorTQ+MkzePz4cbby6dwQTxjldoj2oQZKJuKSzjN7SdOTy+Us54X6KKe9kL93zJgxLMO0Xbt2TAuhxBlaFkIWVLQaUvdnkNrnY0VpYny1atXg4uJSIejs9wbKH17TiaOun4kojShWCvdXuZ\/GRq4NMeiEBz+oQAbjpEmTGG0Q2M2aNcOFCxfKUAmnFzo4MHv2bFDEiKvQhKNO0EnCfXx8WKfEQWLD6UTqCbRzpCyuQEif7Uez35thZobmVIWJiQnLJvj\/Jun8VJzY2uaTIj5RQnX8VAlJepcuXVhKC58kTi86QecnCXigmpYExSip0DUC\/bnjczjecERkfKQgFXTsg1ZDRaBXSkT\/Ro3ERxEpUC6WdO1hUuCCrPmePXvi8OHD7IgMl3SerkjKik7QxZ3xCEjbtm1ZkJXM+WGnh2F3r92ovqg6LBdpzlXyWSatxtDwRa7723Aevck5EoNOzxEfAqDvdH6IOJ8fXhBLOu2BnNOpLVE1RY\/oxIiEfC8UtKCzMOLTBPxliGao8CxUfWEuSq2gP3H5bwf9TU1ouY2UB08JZH7UhD5rn6GvzID+AV03Sjp\/TkrX2c7KgKzd5h\/Q\/wDorwKwrnv+Ab0C0MUbhPiUckVBCb7Jch2eazdkgVIhdYpUzv81i7RSbgDtIAZRC22qdFyDfC\/aQQz+Uxr8BxMqCmL8L\/pe9LGE2A3wf+QBYbP9EkZ4AAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:7674b0e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:713075df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:599eed85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:97a709a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:9f15906d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:3b28743a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:0984f883]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:2f6cc226]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:9e6ab244]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:5bfaff3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:5fd7fc21]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:269b3aec]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg\nTaille de la serie    : 179\nStatistique T_max     : 5.9016\np-valeur (bootstrap)  : 0.4350\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:8c096be6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:2d823996]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expA_br\nTaille de la serie    : 179\nStatistique T_max     : 20.0706\np-valeur (bootstrap)  : 0.0020\nPoint de rupture      : 2012-10-01 00:00:00  (indice 22)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_br\nTaille de la serie    : 21\nStatistique T_max     : 1.7454\np-valeur (bootstrap)  : 0.9220\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_br\nTaille de la serie    : 157\nStatistique T_max     : 22.1125\np-valeur (bootstrap)  : 0.0010\nPoint de rupture      : 2015-03-01 00:00:00  (indice 29)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_br\nTaille de la serie    : 28\nStatistique T_max     : 1.8032\np-valeur (bootstrap)  : 0.9240\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_br\nTaille de la serie    : 128\nStatistique T_max     : 5.9643\np-valeur (bootstrap)  : 0.3940\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_gr\nTaille de la serie    : 179\nStatistique T_max     : 34.5244\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2012-12-01 00:00:00  (indice 24)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_gr\nTaille de la serie    : 23\nStatistique T_max     : 2.9102\np-valeur (bootstrap)  : 0.7140\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_gr\nTaille de la serie    : 155\nStatistique T_max     : 37.5932\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2016-07-01 00:00:00  (indice 43)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_gr\nTaille de la serie    : 42\nStatistique T_max     : 5.1002\np-valeur (bootstrap)  : 0.4030\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_gr\nTaille de la serie    : 112\nStatistique T_max     : 4.6342\np-valeur (bootstrap)  : 0.5830\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:5389054d]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAGZtJREFUeF7tXAdYVMfaflcQEASllwBiVNSoF70aCzYUS\/RaExELhmaJMVhjwAoYsaOiiSWKgLFEidFYohgLGgsmRMUoGrBQFFCksywguL\/fLHM8u+wi+qvJ\/9\/M88Cyc+bMmfPON+98bZDI5XI5XqKkAPAGEFvDPQ4ATgOgz7dVKnMfICvUHea+YcgOnwarWXugZfLO23r8Sz1HUlz2RG6goy3cJC2vgE\/UDfx8M5fVeXexxmC3FjhTBbQYbK9n4PYETQMB7MDaBIkevyuoFzrdiIXxyGAYuy2sNrDiuH14FDqC1fM2T0sKkRkyEGVJ51m9tpk93gm5UA1A1XYSaqxrBKuAA8jdveDvDXpWYanc0lBXAGTTmXQkP5Ii1K0F4qXl6J1TiiJ7IxXACOioatNAwHvBC3fLZmO7rj6753j6TbRYPRI2C44pAUeSmfHlB7CYvJW1e7RxPGtDhdfrOnZSK0EccC1TW1jN\/E5okzKpEZ7mpmmUOk0T+FJi+hoaSzTRC0l0r6oHmJZWIEdPsRocEIsURjAKCX9OIvw7YF9mjjTdvc8mxQU9Kyuwc6ELGgyegfqdPxKGTFKeGzWTSbFEz5BJN7Wpa2IrTIAmeuBUokohrH75cBBfWvvv\/\/vSizrQCb7GVfA0uZSBuU2M4WtWr0qyFVPRtbIH5j0NxM46vbBTiy1uRCAFUfBGrMD4xOwuiN7qB9cG5koUQ6AXHFoD63k\/sXsJdP22\/djfeXsDhcnRRE1Zq0ehNOEYrOfFgFaELP4wMlcMg17r3qjI\/BMVj6tLvFjS79+\/D\/r5K4paSSdYSdIdknKxpVKO8y1NEcQku2oqgoPwTvhCmFuWYtP2bAx2KEF2vRbPyCUF9+AAj+xe2GnO2V8Ot8s\/YfOd+FqBLuZ+ktwH85xh4rlaaZVwoPKiF9VqglSBJbBnz56NS5cu\/RWYoxronFb0Csvgsv0GZrk2QkhLU8QiskpvIWZ5rpdY2TxBls4DAPcY1ZBsd4jbh74mYxHnWPassRfcLrthW7msVvQipiDO3bQC1G3Er4pYXFwcxowZg5UrV+Kdd96+hiNR3Uhnl1Vgla42TK5lQx55HXsn\/AshLYsQq1FRdGEUQnoMfZJGszX3AcJiOmHWaJoMwDb3AVIgr9VGWngynN1DIJclXUJmqBusZ0UzCqlNYfeE9MfTkgKhueoGykHftWsXOnfurNRtClIYPZ7BGaaRecKTfb7OUk3SbR9J8cDCAIZfX4FuXBwMz69EceQTZPeoif9oUAR3IPsdAYA427CzQh3sXvgIZ43MoboBilVGi1nfs5Wgqgpq4nTVdvQcbesWqMi8VU1F5TTEn6EOdAI76plGFqSk9CqgJo1sg2wDAgIC2Pfx48dj8uTJyMjIEOZi0qRJ8Pf3Z99lMhnWrVuHCRMmwMTEBLm5ubh69Sp69+6NhIQEVJN0zud2T+8h7\/RimE84jYxdGSjrTFTxokLgR+A0XJjsN0Zj0MsE4R4CX6O0qFMZqS51ShNoGZrBdtkl1NFXVnPFG\/ev1xIZvYglPRjBDHCSahp9T\/REKlKFSRgnGwfjAGMB9DNnzuCzzz4TQA4JCYGnpyeaNWuG5cuXgyY2PDycgU5\/p6amwt3dHV999RWqGUft9+xB3pw5rLP8SflouLmhAHqd3Dqw9LWEXoIeuy6eDOPlxqwtlaLBHbFkWQQm1tMH3ePiOwT3EhLYNfGL0uA2b97M6p2cnIRBkmT4+voyqaCydOlSNmAqdM3Hcxyu3bipdO1W\/Dl4e3pi8Lt6+BUO7LpY+sSr7LekdAZ6ly5dUFJSgqDwILQ1acsAP43TAp3Qszx8PXAr4RakPaXohV6wNLJkkk6gm5ubY04VVvb29li\/fj0OHz6MLVu2sLHVr18fGzZsgJ+fH\/Lz89G0aVPcvn1beSNNTk7GEG9vPBzcA3n+S8GBJHDLncphGmAKoA1ylsVC\/6A+jNcbIysiC1o5WrAZY8MmodK0ElbeVtDym4C0IUGwD7DHUAzFsmXLcPDgQTawiIgI5OTkYPfu3ayetAlvb28MHjyYLVGxpPz888\/CPaampmwyiIcnvluGs9HhmBGvxSbSqCQLPp9ORTtHB6zeE6P0LJK+mkBvHt4cy02WMxqJYOSoKDSOQ4cOwTrCGgdwAM29m8O1vSsDncZLk8K9KHZ2dnB1dWV1JN1EJ\/Q+devWxbFjx5CYmIi9e\/ciJiZGGfQ9e\/Zg5fr1SIgYiopmX0A3TrcamHl+FpC6h0I7uS8DN3tVNhsgga4q\/drJ2mg20gqz+\/SA98qdoEmlwa5atUrYwMTSziWT14mln\/q+cPwwPD6ZilXNMtDWqAylTyUITSFBANws8hF01wqu79TBou8vgEszX1kcdB2fzRjzyXTcuXNHAJf+eGLzBFMjpmJ2s9nCiuITbOlvCR+ZD94PeB+d0ZmBvnjxYpSXlwurk8b8559\/sh8x11NbAlwqlTIhI8pR2kgJ9LVh6\/DH1+WQtotXAn0t1mLFmBVKA6Uvj5Yao2TIFbYKDA8ZCtdJ6i0rTKH18XMXA79IdEGFliYBTRsOf0GSdD45fPAcfLHUk\/TShsU3Nx+3IZg0cTw8LLIx0FyKq4W6+DzZRpggeh7XYlTpRZ2kc4qjVXXL\/xaiZFFwCXBBG7SBh4cHFi5cyAAWl1atWjEV1Nramkk6AVyvXj02SRz0TZs2KXM6gU5AcK4WS\/ofpn+g9xwf5HrloGQQ8SltlbQUafNUWKRUDPYYwGKOBUqdStHR3xuPZ52Gh20JPGYsEPR0MVhEL\/RdDDrvi7ejJU6T07NnT6UNUNwP1yiIP4n\/a1IL+TXO6T7hPhhiMoQ9ljidNlIOuqm9KbYt28aujQwYCUsoOF1V0vmYqe\/Ro0cLeCxatIjRy4ULF0B\/kyApaS8kYSPHjUXq8NvI889jnF7\/SGtk7OqHrtZzcOn3M9DKSUdJnxRo5fWGweHfUNJ9HXQvPxb4vdK2sor7gXoLYjAscCtkN2KxaMIInNVuJfDzDz\/8IOzwJMFc6qdOnSpIr5jvCcy+ffsqTY4YWOJ7oq6XAV28qXvDG5HMACRPaRAaoRH2Ld+Hq4eusn1rGIYh0zsT7du3Z\/sQbfJcA3J0dFQal3hPIknnq5Huo+\/V9PSv9xxH6JxPUNqpE\/KmtEVpN3\/UKTyLp0Y9IJHJYPDjj9B++JANrmiMCSrNuzOPo+HueJjNu8HqiR9poF3fleH7rAfwHNoPN3OfsGv8RcUUYmNjwzQB0gDEks+1F9qQ+IBVNRvV\/l4VdBqbGHj6LpFJBNqkVUbFyMhIGAtnBqoXj5ELA9WT0kACRkUA3XjGKXkzC30cnNIW5OI9cTMHg2zaoKzBpGe3NBJMfxdm4BOlKAqRiif7TGHaOKq86bQ0FQ4vByyUJSBQrnBmqXoZhY5e8Q+xcaTr2BUWk7cgK3QkcxW\/yHp9kUVKRhLZF6Src739FYep9rZqevrDojI4fX8Q2UMHo8Hii6g0KEDx5M\/Rw8aO8Z26ImGcTg7V59xOi\/R4+kj0kBW+0FX7si8kNo6MB81AbnQwLGd8h4Ija1Fy9TjzXKoaR+Jn1AT6y47lVdpLfk58LO\/zzIvIC0WO3j+djlvdcyHX76DUp6oeS9JAS1Ih2eTwok9v1C+xR\/QiYzjeURg3Oo5dYRtyjv2tLlpE9WKPITfXNdWLde6Kx\/cF0OWlRSxkZ9jLE483T6yGh6r2IuZ07nPRBCKt4Nflg5HM3HtTTlEiXkjSW91+jALDp6hoal+LieTOIPK6EOH0woLodhiRroMS5+mY5uWOsOYZeNJ2OJynrVEbLSLgeORIlhQn+Nmf3L+ptp6kmPzplTn3YeIWiPwja2DmGYrMkA+g06gtylOvanQH0wupk3TaREmANBUymkjoeBHbF7QnEXdzFwC14X4Y2oOCg4PxxRdfYMeOHUw9lgzdcFm+07cNeJw09NwtzNVtiDq5MtSLSUFZlx9R4hamZizcyUW8n1rF\/8T5KVh4Wo5pTk6YNWsW5s6dC+yajnNJmRjgPQ3l+xZWixY9Sb8h0AJJKw\/XlSTEqK3nnC1eNTRAco4Zufoyaa8pMF0T6HLImYUsDnCQ+U8aDXE8FTJ2aJPnquGDBw+YWjhu3DgcOXIEly9fZtfee+89FBUVISoqCqNGjQJpaYMGDYJk5DdX5ds8Wwmg8+gndf55eQX+s\/4U7trtRPrI7WxzfPeyP+7+m1YGD+ZRywi29Cqyo3DfPBbm0q8QlNoYP8+bh62hi1EcNhK7ylqjv3NbWP55RG20qCw9kcU7xYELmgx19WKfuzrJJKri96m7XhPo6ffT33iAoxq90CD7F5ThbFwG6t7Oh1anQyhvJsXoIj+E21YwznZIcaj6Aa6PaI\/4+vuq+DwC3cJPwsItET8Y+cF93RWszIyGwQA\/BO+Lw0B7LbTO\/fWNgs4nraZwnaqLgN6Z08vFuItvPMAhOXTtoXxQGwtBIKIvZ8HLuB48zeohMekkLj1di3dbmCHRUIdRh0tsLDwpEeCZGV9mZYUKW1voZmZhwuKmVcArbLrZGY0QbxOFOyN6of6H8xFyqZCB3vL2vv81vWib2bIw3ovioJr4uSZJ56CrC3DUYoOrVRNmkVLLqd\/dwrpRLRBz4RoWzPNDnuckyAb8iCY\/zcefU78Wov9BwbFolAK4xAJZmbpI1daGu1SKnmPHIll6Hk3D7fCLySIkjxqADut6YkCdbVj2ZV98U9kZw4YNg9WRudXSLl5lI1X3drS51nd2A9GSOLjN29akvahK+hsFXTUbYHvMbwieOQkPF\/pC6v4YHc\/0xp1OY5Cj1wdekZEIDFaESA84aePHOCt4P8nDx7nm6HX6NFIcglm4onveJ+iUMIpR0T6Xe7jj1hg3tBuj69fx0Ek6LSQYqVMN6+g3ECL8YpVRtV4d6DwFo1JWCMvPopB3eA2bBNoDaEJ07d5jYcCaJF0vTg\/WY6yV\/P61Et+XaMTcAGTOcj\/3zhOXERm2BCU9KpH9jRTaydaw8ciGVnYl63ZQs0doYSbDruuWyChSBDOMBuTDLq0rftfOQoVdKmx1GyLrdlN4zGmGlR3b4dN+ixGfVRc79h2Era0t80WQE4uKOEAhHreqp5FLnipg4rFf\/H4LNmzfC\/cW9bE+oYR1N62bHe4at2XPa2lSF5F7fsCd7KJqkSOup1+Ku4SYMTH\/K9DFjjquTtJYlixZgtDQUIXvRWngd\/Ox6PNP4RTQBPs\/1IFpQB7qulkjs9s76OFxCjlXc9CnTx\/2QjGxsZCR3vmfVWjt1wC\/a9fFx+GuqBObim9X38TSzUvx8bvnMXroURgZNMTqXYdZ7JCHspKSkqq9PPWr6oVUd8+aDpVoI0\/DT9kG2JFpjCVNs5Ao1cXqNAt88b4Rhg4ZgiUn7+JI7AVELV8A\/YfXMGXDAQxzHwuX\/oPUPpee\/TqsVcKTCnk7uZ4eGBjIAjfkh2FugMP79z2P6EjLMcbDC48nFaPIjygCaP1HNvInn0PdjELo6+ujpScZQcDNqChMLi6GQa9URDwxRrK0CeaHe6PdtVJ8FLYbHw3vjU1dzNB+6EaMtXiEbg1lmHfXBt0GfIi5y8Pw+EEqc4YNqP8AQ1qaCjmL0pxHmP6RC06mlaK31VOs\/f4k6ts0Yerk4U+cWbRoXW8TDFx3CtGHjgpjJ\/2YVuzWtctguD8A3+q44Pyp41hk\/Dsbb6C0G7r2GSi4iKdNm4ZOnZSzDEjnppwYTZxOijKp1TxJVh2riA0nui42nmhCmBsg59oJYeCyJ5X4cMxY5Aw2ReGMf8PS9xxMsvNhFLEfhuv7ouCQBMuqEn0DJBKMk8vx2FuCU8kFSJbaYXP4Yhhol2HExgNwuncbd9t3h3XoHmxZugB6rV3g6f4hsvKlSmMlX7ln3SsC55KeffFiHD47cEtoR168GTaZSJabYtI3MQx0l\/FzcCS1ohroYuuQryrqiPvsx44d+0JdXBPolG7FQSfHh7oilnTxdRoLxVYlk3del\/fUuq4k6WPHj8LAyUOxpdUqFobrt6sfrDtvwso97dH280dYrKWFD3bsgPns2cibPh1SdxmMl59Ag7jz+D38d6w3WY+lcXJYzYiCiY45utjaYtXWnYw2vMe6o3VlGgK3HUT2Rl+WSkebG4\/WW0zfg0dr3YX69OPb8cncEOYant3yCax1yjEzwQCrmmbgX0ZlOPZYmV52PDJH5PZv0aJDN6VYqxh0MtE1pdVR1ldYWNgLJZ360wS6mNOpHUW+KGVj48aNinDd4LCLcrvkXUi4cpn5Dy7ezcf8ueMgfb8bCmclwso7G3l+mUjrm4bWvm4wSE9CxPVs9DlxAtbjxiF\/4kQUj3sKw53H0SAuDQfCD2CYyXSk5XaCve9vkCTcU9osg72H4Leb9xC57zDiAj\/C1FO57AVb4wFLKLWadxzJa7ww55oOo4IpPZvi5OqZjFLWOuuiiediTApeA98ujdCjTgrCip1wOeEao5TYlZ\/h29S6WNr0IWzL0rHlvjESivQQ0lTh\/5+fYofuwz0QELhYg4y+Hk7X2HnVBYmj0\/vyYQP74PzZswx0Km7ebujh1wPr3U3hNDEJhScusvpW3bvjxi+\/wDrMGksbLGXpExcvXoSFqwVkDjIUxpfBObw9dpnoMLdA1+V6SDl0SHAGkTQ\/3L+aAXX4qCItOmjCCHw8ZwWTdA46rYDH7UZjytrdQpB3vPtQjCs9wXT8tftj2bMt6j5BR6MS\/FqoL2yktZV0TcC8jo30haCrehn5DYrkG2BmtDtWuw1gTEZetpSqZGlKmX4esFAcCqAkaip2hSsw\/FgipkUHgeviqv4Q1TxFTfTC680nhyMrpB\/zHtZz7Iz02f+Gjv2\/UMfQlOniVJc21RHyMoWqqK7UJj\/9rYCu6k+nyNHR69m47DYZ5ZU9kHpxAno3N0ETc\/2qzF3F\/k1ZixRboWQ0nqvuwNyfLlURpeevTYBTUU0CFU+E2HjRVM\/b1Hceifv+nWDqtQqVeZlKTjHyMNJE6DVpXy1LuPhCtNIhAnUT81ZAJ5WxuKxScANQyG5W9C1seXQUsv5RmHrlO3Yqgwr3QEayYFaVAwaN2MEXL5ZqSY5+5VJTQic\/DEBHXSjkxiM+qqE4Xi92Zknqm6HhQD\/oO\/UXkkzJFVx04hvUaWhT7VCApoMEqsC\/FdBVs3ZpEDzgeu\/uPdifsseoLqOgX08f5Pbct488ioDTkKFo3LgxKipv42DBQZh2UESfxKls4gCtJsvzRfyn6boYRMr05b4WojOS6JKLe4REUloh9F2\/i\/vfQ9JVfS9kfnNzldLDAjYFYEefHXC3cEfixkT4u\/vjvvZ9hMaHwuA\/BnA95wr9X\/VZpIRH6slhT+kSFMQg3fjhw4csavLpp5\/C0PB5QtKrAl6b+2SJZ5FfdaJDomcAU4+VqGvb8oW3cuNIneH0wptr2aBa1i7dJ5ZQCi9lZWXBeYwz5ibNxdTcqejevTumvDcFcXpxOFBwAO0K2kFHRwdmZmbMt+Dg4IBGjRrhyy+\/ZCkLf9WJh1pi8NabVct7IWkVwmxVVpyXlxcLPVGaL2UqfbDsA2wcuZEN1jXdFbOuz2KAU0rZvHnzWIiO\/Ay8\/V914uGto1nLB1YDnaScfAdkOVG8b\/v27UxiqRCIZ8+mId3XFVb+m5D13QewGnUMG25ugKOOI7O4KM+P8v+oHw76m\/RN1\/I9\/1bNqtELgUXZpRYWFizjKj4+Hs7Ozowupk3bj8TEjig2aw\/vqO3IyuyPo+4yLL5Px1xSYF1mDZ8mPuwFqZ8FCxaw7Nh\/QFeec7WSTmDNnz8fHTt2ZL5vAp02Rtdvz+Jqfj4Q64WfvjiK0NDWOLnVDnAIZlHGIEQgsCoLTJz+9t8EOr03HQqgnEzKW+SFNEISXGKBaqCTeU05d1QkEgmePn2Krl27YuvWrXArdsMRc8pXvIevExNhkVgCN7cOsAooRVZzPVjPHQ09vUuCr2Xbtm0su\/W\/CXRNPFYj6ARUZGQk4\/SCggL4+Phg+vTp7JOfy6EUup+OHhX69\/fvhD\/+MIFMVgo9PUU0SawFcdBTUgD6qU1xeZ42KTSPrem\/QKh0+ir3OzgA9CMeP+1Tv\/zyCxMcOjfEj7twu0P1UANpe1zSeWYv+dNdXFxYzgv1ofbwrlhlpIYk6RQFEWdBRcdHw6rUShgg6d\/NmzevEfTgZywUJP6PDTWgr+5\/c0jEqZIvmLlXuZ\/GFvj8sDbbl8j\/zQ8qUNIQaWdEGwR2u3btcOXKFSUq4fRCCsWKFStAESOeLk0JSBpB5+8jVh8pZYxyFinbyVxqjqiz3A2gaE0qY5s2bf7fSTo\/FScWRI6P+EQJ1fFTJSTpAwcOZJY9nyROLwx0dW4A6oA74ml2+AFXenCXzl1ge8cWi84qnFi8kKZDK6ImeqkNrfyd2oiPIpKlLpZ01XESXnSscfjw4Thx4gRTnbmkk2VPK4Ny3BnoSUlJcjrBQMcyeGI7dcANJJJwXija0t22O7LisliagrjQLNM92trP\/3fM23AevclJEoNOzxEfAqDvdH6IOJ8fXhBLOmkvnNOpLVnxFD2iEyNqLVLKMiUuohlSLZrCXJRaQT\/i8n8d9Dc1odVAV41k04Nf1UP4D+jqp03jP9l5HbP8D+j\/gP465Oi19KFW0sUbhvjUsrqghGq6gTiI8d9okb6SG0A1iEEcT5oNHd8gfww7WSHKy6PTBSkpKWqDGPyfKvzjBlBoPtz38j97wnOzYjYe2wAAAABJRU5ErkJggg==","height":56,"width":93}}
%---
%[output:7d33eba5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:2c2d1b64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:4ccbbbdf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:71b1e216]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:0b10475a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:81f74aba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:7ab9e0b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:66b3b500]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG0_homo\nTaille de la serie    : 198\nStatistique T_max     : 32.9093\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2023-02-01 00:00:00  (indice 164)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG0_homo\nTaille de la serie    : 163\nStatistique T_max     : 53.9837\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2012-05-01 00:00:00  (indice 36)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG0_homo\nTaille de la serie    : 34\nStatistique T_max     : 4.5367\np-valeur (bootstrap)  : 0.4560\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG0_homo\nTaille de la serie    : 35\nStatistique T_max     : 6.0858\np-valeur (bootstrap)  : 0.2620\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG0_homo\nTaille de la serie    : 127\nStatistique T_max     : 14.8477\np-valeur (bootstrap)  : 0.0150\nPoint de rupture      : 2015-09-01 00:00:00  (indice 40)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG0_homo\nTaille de la serie    : 39\nStatistique T_max     : 5.1927\np-valeur (bootstrap)  : 0.4250\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG0_homo\nTaille de la serie    : 87\nStatistique T_max     : 7.2846\np-valeur (bootstrap)  : 0.2090\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:67188cfd]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAF41JREFUeF7tXAlYVVXXfi+zIrNMyqA54ZBDOWCGQ2qYiZoDmEFwxYkylAxFnAcUVFCcTRHpNwMtNIHUtMwhhX40jBRFTRAEAUFRBlHg\/q593ddzD\/cCKn1m37+fh4d7z9l7n33evfZaa79r7SuRyWQyqCmLFwOLFgHqa9RsKJHI2yxcqLrTqqKbuBXqBnPvcBRETIfVzBhomjZXN4R\/5XVJScUjmb6OluLlSh9WYkLUBRxJK0L5mRZ4kNjybwG9qTQMtyM\/f27Qq8vuITdoKCrSf1WMXbdtH1jP\/QEajQ3\/0ZMlOXLxtmxQezPFILccz8KV\/FKEjnWAf2AlVq\/QwuZfsjC1n229XqQuSadOboWNQ9mZGLX9aTW1Q\/Og02pXAAdc08wGVp9HK\/qhfqsKs\/\/xwEs+35MmI4BVFa5eGhp0ehZXM8+jXtS1fZE+6yVRDVRJ4vplimyHZ0cIVQzvm0v6ocRS6JbrIyMDyMwEfvkF7DMvzW1l+OUYoKUpQV2SzoEh9ZIXOhaVt2\/UeJW6JJ2vlgfnD8F67mHotu2FivQk5AY5Q6\/LECXpV4dTdnY26O9lFLWgk27vNeY2LsRZ1WtcLVoA9EcTQoZ0pvN3yA8dw9qauC6GydgFNfopSVRfhyYnZ+kQWPhsZ6CqKnf2LsGdPU8ttrrniNsS2P7+\/khKSqrXuzV0JYlv9AVZuFsHpX7z7ldg9OrLOLfBAeVFOs\/8zAWz7mHiwz4MMCr5myei2fxDSjpaCKq4Dpdaus4l+ZkHUUuDxMREjB8\/HqtWrULz5v95z0my8tBfMn\/nloohEuAe4VdxZH5HxbX+\/eUSXN\/i9EYBdvXszoyhRM+AeRlGLn5o4jha0cXd2GAURc8FZNWKa5oGTWE97xCK9q+C0TveyNvsDeuZe9VKen3HI67HQd+9ezccHR2VbpPapHc9fly+cj095f8bsjBDOmtIS\/hGX8K6cQ4IiLqBr+baovqeHnR1E9Gs2Xjk5OxGRYXy4OoahI1RLiZ8aor+fSvR9oQzGnd9V6FiuFrRNLODbVgq64omRqKtgwcXjsFi5rfQNrVBbuhYtaC\/iMuoCnQCOypKrhrFxcsL2LSpHAEBAezWxIkT4ePjg5ycHEXVKVOmYPbs2ex7eXk51q1bh0mTJsHU1BRFRUVISUnBO++8g\/Pnz6OGTuceCzV+EdD5aBYEVmAyBipA52DptuuNistnmHvHQaeJ0bbtiOK4NTB1W4K8jV4qQX9Rl1EV6Py9SappZffrJ3ca+CR4eJTDxOQp6MePH8e0adMUIAcFBcHT0xNt2rRBSEgI6BkREREMdPqcmZkJNzc3bNiwQRl0mu1OnWJgYTGHdXb37hQYG29FWNhu+Pk54saNIlhaekNP7zy7X1a+E2bmnXEx1RgmJiGsLhU7i754kCVFUk4\/2NvdRkfzPkgresTu7dqyDvbHg2EwwBMhK4LwzXUNdr29qTa2LJ+LZt0H4lLQGASm6iH16uO3fjz5K1asYAOmQlIzwdMDf1xIU7p35coVSD0\/xlDzUvxWbcvuC6VPKL0c9N69e6OsrAyLFkWga1dTpkaOHXuqTuhZ7u7euHTpPEpL+2HAAMDS0pBJOoFubm6OOXPkWNnZ2WH9+vWIj4\/Htm3b2LUmTZpg06ZN+Oyzz3D37l20bt0aV69ehZKf\/vXXVxAQIEVpqQvu3JmtAJLUy8OHXWBmJp\/pwsJgNG58ACYm63HrViQ0NQsVaqiqygytW3nBt\/1tfPp9Cmvj0PQKDv60Hgknk7E+fC2Cu1VDNvBTfLVmKVZtjMDNvAJ4z\/DHSNfx8P\/MBws93kPqIzMsaX4Nf3bzweY9CYiMjISZmRm8vb2ZHp78WgVO7I2AX7ImSDcba1XC62MPdGv3GsJ2x+PAgQMMBGpH0lcb6O3aRSAkxBSkRiIjn9YkiY2Li4O1dST27wfatZNi4MA3GehSqZQJAGdRbG1tMXDgQHaNpJvUiYuLC7S1tXHo0CFcvHgRe\/bsweHDh6G0Ix0\/PgYnT8qBrKxso6ReCEwrKylKTD+EY9v3YKKdi5SUL9Cp03xcTi7BfcgnhOv\/gih\/WHtPZm3u3PkMJSVuuJR8ClJPTwT2NMDQ0ASUpyciyG8SYvKMWVuSzC8muTPQd6eVwMGwGju2bYNVj8Hs\/ukf4+E+1Rer2+Sgq2EFHlRLEJoh302PtbiLRX9ZYaCNJpbs\/RX\/m57FPBSxsSRQaNVcu3ZNaSIePWoGX99I+PvLJ4jq8Qm2tJyNCRPK0aNHAMjuEujLli3Dw4cPFSqEJujy5cvsT6jrqS4BXlpaim+++YbVlxDhRR4LGdLKIxfx25mNKkGngZBRFZdPPvkE0dFeKLgZAM1GPytuT2rTBYuPz1bZxr9DJapKihB2wwJulnfhalmMuVct0cXgASbZ3EG2ri0C\/2qOnFt5rL8uXbqwwR45ckRJeslgiY0bLWUCtTYPRaxeVEm6EPRLl2YjKqoc\/fsH4PXXAXd3dyxYsIABLCwdO3ZkLqi1tTWTdBpzo0aN2CRx0Lds2SIHnTf09Y1BfPwchbQKDSmX9IoegzG2mwW62+hBR0fuw3\/xxWjcuqXHPuvry22CpY4xzl3fBguL6UzSQ\/1aYcrS7qyOEKzg4GD2nUsV9wCE9WiJ0yro16+fkvQ2FOgTJkRg+HBTNjbS6WRIOehmZnbYsSOY3XN1DYClpWpJ5xjShH744YeKuViyZAlTL6dPnwZ9plWgBPqqVVewbl3tOt3m9Uzcyb8JlBRAQnt+ZnAnw9g4AdnZkaiqspHrcfNrOJUWhaZNl7A6nn3fRw\/XOwpJjY2NVVh4kmAySASsr6+vQnppQmj3SPqTJHjw4MFKkyOUZtL3vF59JV2oeqRSYOdOOVbksdjbA999F4KUlDi28keOBHJzpXjzzTdB4yLXj6uvtm3bKo1L6L2QpPPVSO3ouxLoTIXYnIaejjt7eFbWL2jSZB\/u3x\/NwKTSuHE+02WVlfLvZHzI6icknEJqqhasrceD9KOngyOCjq6Ghoayx8NflHkbUimb+WbNmjFPgDwAoeTTi1Ehg8QHzCWQ3xP39yzqRazvhcDTcyWSciZABgZxbJVRMTQ0VIwlJiZG4b0Ix8iFgeqTIScBo6IWdOcPKvDjfl3F8ujplYXfdj6ldYUBDb6h4BIiJMFihrvB7cBT+vaHjSfw3id9Ff2+6Afh5oh4dAufbbgV6lorVyNUAaqMLN3n70T\/CWfut7\/oeIXta0g6GdWWLSQKzqV5z0Jo3DJF1g25KhFHkWhwtLHgwFMde7sqnBzVFXZr5btNG8McZFyXNViESLg5Mhnmh6K9i2HpF43ihLUoS\/mxTj69NiPbkOCq66sG6FRRKO2WzavRvo2GgnsRSzotSTEvQyrHf\/hBdBz1HnvuR+9nYle8fY0xCFlG4c26qF0hb155O1sBuuzBfRYKrIujV0cD1MYvkXFtKA6GgS7UTeSefTxnDcZJ9aBhWI7qci08PP9UvfCHX7tWjZMn5btJKgQ0l3Zb277o1WsTvv22Exvo9es155yAuzn3LZh6hikRYfWVNB4lMh27EHcT1qCpZyhyg4ZAp03vOvl0VaDT2EmA1BXaNNE78kLGcutW+Q6cbBLfhNF1KtwLIxu0ePFizJo1C7t27WJ8jCQ9PV22fPlyhIaGsp3U\/BWh2HmuO27+0Js1buKSivtxnYCyJ3r+GABVrBsFNQbUHLK6IHVDRHlopfyxYhzyKuQxXoMBUjTp93Gd83bz5k3Gp0+fPh29esm5emIVaazkMooLbf\/p3hNbyjY7ZMi5a0j9kVvo4eGBhIQEnDt3jt3r0KED7t+\/j6ioKIwbN47tM4YNG1bTe9mw71fMW9kUQ9pp4+jBMtytboHqcm3ISnWfgs2jRjTzwgkg+vcJBWxoWISxY8vg7m7D\/F5VhYIQFVkX65RMdSi+7GBEnbOrpkINnT558WYc2D8c36zPhM9kLRRVtcdrbauQFG8MkJrgIBO4JNn0XQS8k1M1xji4YNgdOYOoKqLD1cvzhuuo35cdjGgQ0Em3\/3SpCEd+dMNyv+tYvKAKGq\/ZwmdGOa7GheOb5ctRUKAPPM6HYUB7Pnks6UKKmh0HnH6S6\/rNA0bAdkxPDBo0SGXk6HkHTO34hJ29no+Z6dY1+JUX6fs\/0VYh6WQAWrRogQrj1zB\/6nRYWkpRfM8IGi3fRuhyCQ7vCMa6A+FysJ+oC+vME7j3IBel7dxgsi8EevmJsLffipwMfVhUn0NyvhOqSuX5KeLIEb2cOMbJX7gu74XXawjXT7hJ4\/3yTRMJob29fY3oknhiOB1BdAUvnC8iOyl0VOg+SzbasCaU7biIMk1MuQC3DyZh4DAp9v\/WG7r6BohdU4awZTsQ+3ApEAkY61RD6\/pf0Js6BY2IzBkmt6DZFi0eb9cLgfxrGFXtgnnZmVi6sKpG5IhLKw8834lfgyZvjWVeDHklurYdVAayxS\/bUKAT3coDEjQJZPjmzp3L6OH6gi6MFAnHSYBnZGQoRZUk62NPyVb7faTgUSoqOkNbJw\/3jDzwoPFQaOf9DDPZMsX9kq6e8BrUGUnfbkJWhpweHdqvNwwuHcShu51xOXMuHKznwcKoEF3eOYmk38wxzro7vrpShaiYWNjY2DAugktF0LxZGGKQrzCo3Ku5\/0EwJs4IUNCkXPqEQNOzySsgCoHuk9dAHDqRZ0uXLmVjmz9\/PmP86HlC6ePAEMhi0Ona0KFDmYQSYNw1pGdQH8LxU4Bl+PDhSuE53jd3FxcuXMg8Q15YEKNn9XkFEbV2XTWO\/uIJHY3puHF\/ICx0lqKjw0coLnZEZk4QKqvi0LLtcpgb6ODCH\/6wbynF+a5S6CQtxKObaci7EoGj27ZjzurN6Nk3Fp1eb4GDO4ehWVtLrN7+NRscAbdtw1okLhwN35+LEBUyH7ZHF7LQXNn5wyg4EoENmu9CotOY8RW8DVGl6enpjGha070K1ffy4XfJUgl0Is44EBwcAosTYsSRCJlMVeqFR5yEUspXAJFpwlAdAalKvfBn8lVDRBdXMyxG+n7ji\/hy80bm4AcGytk0l1E+CNrogpHvacFU6yZOnaJISQ6qGpli\/CRfdLVpgtUrgzF6jCvCfaehYkk4jBIScePPCNxKz4XX2OEwsPXH\/GnmGDUlCIvmTYbzSFclNq7g2gV4uY2C+yd+TNp5DovBtN1YFpvIpFNIJAk9FnopdZLONypdQ0JQ\/DhWaRsRgf2PCTxV9LFY0glAHu+klcPVizDYTOE4Lv3PJenqQPec5IPQY61Q+WswtLVu48CBSHz0ySnk5WzECC9fdLDWx5drgzFkxFhsH+OLRnOXQKM4DbJrEcjIMMXS6RNwOCEDI5slY\/utDohL+B8GEmcWhXpPVSxTyNRRXQ6+kFKtC\/SWISGoSEyETkQEzj0D6Fw\/02aGCkm3UNeT1FJRtdsUqhGqI9bpzJCK1QvleWhoSEGgbzr7EMXxS1FdvRtSqSMSTkUy0J3eHQEHC118G7MHVs0m43DaZNiNCESjzDRcORCBqip5BHzixC\/Qpo050tPtkJysPlghNpDC70Iql\/R\/40MrWFyUQnbE54vVizAuyiVdVgfoYkEQei8\/\/PADTp48qdjqiwWnNknn71HDe0n4PVt2OHIlzp49q1AvWVlSePv4YPMFPZTFh0G\/kS927BgMNzdvPHyUi8LbHqw\/IyMC2Bj5SV\/B4OuvoXc8EXmH5aATWM7O3igsPA8TkxU4e1YezQ9ePA8n4\/YiYutGJIV4KwCkmCcvJYZ2WJDfCb3f7sv0r7pgBYXFKO2hadOmSoZUyIPwVAjqW5V6qW3C\/657kj59+sh69XFC\/OGfELUzEtcu6GDFCimm+\/tge25r6Bw9gpt\/7WDPt7Z2Qm7uScybN4+lE5BeO3PmDBq5uKDUwgIPDyajtUYETpyQW2oXlxAUFMTB0TESa9fKA75ioyNMrxC+pNjACVUQJ5uMNStxp1Lz1QNd1UkMIng48UOfqRBbSNy5Oh6F6hBTx3NH6Dvx7BT2EgYCatv+c9Dr2hzxPMiMPp\/Dc9bSV29HWlZWJuOGgb90YGAIoqOVacvJk9swwFUda+GBDLIHwkkh0ImjVsXc0bNUbYSIOSw5vbdWEoz78pn9A+A+xffVA\/3MmTMyYRIlX\/7t23+ImBhHBWAELPHNYgB5wqU4rMXr0ySpWh0vktjPI0dphp0xeeuhVw\/04OBgGXEuPD2McjeMjIxAu6g\/\/yzBxx\/\/BS2tbGzZYoVx49wYgPT36NEjxMb+jrS0norgdFZWNr777ju2YLp0GY7MzJYqAxh8RfFjMJyF5N8b93ark+6lSTs4eyR8fyp89UCXSqWyyspKrF27lm1VKQWC8vF4cXV1xZAhU+Hnl4L27d9DVpYmKisfIClJi2VvRUf7oGfPnrh37x4+\/XQVrK0DEBtriN69A7F16xyWK1jbiYfyiydw90liv0RPH2buq6Bt075ejoOqYARvSEEJ4YkRivrwIASvQ+1px1tYWKh4Hu0HhHkrwoFQQIJyEbt168YyxPT09FBcXIy8vDwMGDAAv\/\/+O3MwqPB6ql5E4uvrK+OpAVSBXCxK\/xo1ahTCwsKYf\/rWW28hOTmZ0bQUNXFwcIClpSU0NTUZ1fv++x3w0UeVOHXqZ6xf74hBg7QZP0H1yed9WSce6jVzL6GShNSLkIsgd6yqqgrbt8tPUbz99tto164d8vPzYWVlxTJSZ8yYwZJuKDly4sTWSE01Ra9eD9C9+xfsdAMVAp3Xf1knHl4CnvV6pERoSGnndOzYMaSlpbGMqn379iE1NRUUQ6UtMQdx5syZTKXwtLp582yQkmKEESNmqARd1YmHeo3uX1qJBTF4JisPplLAtqCggL1yy5YtWbIkZSkRrUnZqmLQk5ObYNEie7zxxmi2q+WpZLz+fxPohCVpA7KNQldcGBBhoJNKUQQxnuh0UjmkIrS0tNC1a1dmJIhSpZ0hrYIePXooJJ0mp1evnnB2Dsb69fJzRbQ6eP3\/JtDVLU4l0Em9EDi8UBIvpftS1IRUi5eXF0t8j46Wn0ymuhMmTEDnzp2ZTudl377OeP31IqxZM4Jd4tSruvS1V0FzkFOxefNmRnjR+9ARFu5ac\/pCmP9CAkl5LVzSOSNKzkj\/\/v1Z+gX1USMbgB5EHfHcapJ2Uju0geIJOeHh6xAXNwp\/\/tkU33\/\/B8NPX18frVq1Ym4UL+JwGm2YhPmOtQGvakP1LCf8nqc9PwsrHD8FLGjVi6ldwohcR9IAQlXC1QvlsK9cuZLtd7i6JRwZ6GIagBqNHDkSWVlZ7Nm6urosM6l79+7YtSsbHh42mD370OMAw5vo27cKU6fKE\/epHm2yagNdeIisLklX9csbTzKz62rK7j9Pe3FilPCAlpjfp2eQtFPh0i+UdB7uo1grjxpRQISBrooGoJMDxB5yNcFpAnqwk5MNY\/XMzcuwalUBrKwe1Bv0V1HS+ak4VUEM4ezziNMHH3yAo0ePstMaXNJp08ntJgNd7KeTpJN3EhgYqPKAlKvrLLRo4YlhwwzQvLn8xBwVsgM0ObVJer1E9B9USSjpNCxxMIKOspDO57ny6nQ6tXVycmJnTylWq3JHKjSsQh67tjQ2qkeTRd6OOp3+D8LzpQ6lBujiJUOGlGgAfo5T3a9HUGoF\/QlLQ+SlvFR0\/qaH11Av4ueIU3+fZRz\/D7pqtGoYUgKKu0k8KEwuk\/iHC+oDPgddmJJcn3b\/9joKGoAS13kmktjhFxJizwLIq5rK\/Czv+Dx1VZ7EEP6QADeq4gAyjzDxjRM9XOjL8vpiG1B9\/zZub5+GqmK5f2\/suhiNOjzfATCijMPDw2v8bgv9mgWFDsVFfJ3IPSrEhVOhU4MUhKGdI7l9ROqRV0Y8Orl\/VGjV8uR+AwMDdk14KIBIQKpPMYnRo0crUSV8PDVOYpCU0zFDOq5BUX9yHakIT2sIc1E4ryJ0NcX1hS\/fEIFp3p86m0GbMFWxXFXXhaua+uXvQ79W4ezsrHCb+U6drvNjLDyxSOzDqyO9+Lj\/D7KSRgwBlyfkAAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:2412e4f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:20242383]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:17ab63e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:4e4b1496]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:8fa6af9b]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : SSAB\nTaille de la serie    : 179\nStatistique T_max     : 3.6531\np-valeur (bootstrap)  : 0.7820\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:476babb8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:5cc08166]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : SSAG\nTaille de la serie    : 179\nStatistique T_max     : 4.9030\np-valeur (bootstrap)  : 0.5680\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : SSAR\nTaille de la serie    : 179\nStatistique T_max     : 15.0768\np-valeur (bootstrap)  : 0.0070\nPoint de rupture      : 2023-06-01 00:00:00  (indice 149)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : SSAR\nTaille de la serie    : 148\nStatistique T_max     : 19.9408\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2019-11-01 00:00:00  (indice 107)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : SSAR\nTaille de la serie    : 30\nStatistique T_max     : 3.3676\np-valeur (bootstrap)  : 0.6530\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : SSAR\nTaille de la serie    : 106\nStatistique T_max     : 4.0345\np-valeur (bootstrap)  : 0.6740\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : SSAR\nTaille de la serie    : 41\nStatistique T_max     : 2.9111\np-valeur (bootstrap)  : 0.7830\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:73358e43]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAAGRxJREFUeF7tXAdYVEfXfhdRio0FqSJgwxIFo0axxF7QWGJBUFFAJIJdUcQGohSNgIhiA0QMoqgxhZiIsWBJwNjQ\/DawgCACCtjoZf\/vzObe3F2WZkh5vi\/j48PuvTNz575z5pwz7zmzIolEIkEdSzziMRRDWe0hGAJPeLK\/XElFKhzgAKpH5TzOy9yv42NqrFaR9wxZAdbQdtyOF+FLoOcag0aarRui67+sD1F1oBeUlmNO5B04DTTEiC5aIEDboi0b2AZsYIBXVwj4AzjAbj+yAsTTvCC28qhS\/V3il8gJmMquy9chcDM3WULHJQwqpn35thzorRwC8TJieRXQKwvf4LnPWJQk\/8S3UTEdAP2130NJvcVfBmxND1IIevbbEkwISUJKTiGOOpkx0AlEAtMEJjUCTg+rLMiHZ8kqZGiWYdabT+DjmQaD9adkJFIIKrXJ2T2Xr1OSfAXPfUazceuvjZMBna5lBdqgMCGm2vcSqajDKDiZfx7Vr8jN+McAXwV0knCPbx5izoDWmHPwLnwmdmCgk1rh1EZ9xSVj7QC0HL8MzSym8E1JyvMil6O1z88QqTZn0kl11Ewt8OKAK1oOc0T2bkfoux6rAjp1wkm8UL0oulZd3fq+Q0PWr1a9cNJOoHfo8pZXLfTwwrLyKmNorKQE5UYidp0AbW4hVRtUCHT1HqNkVAzVeR27jUkfFQJdWIdJe4BVFdCF6iU7wArlL59WxUOkBH23b6DWexy4VaNqbgm95UdqxS4jIwP0\/88s7wV6j48KsefgCzYu5zbaeKGmCpO8NGY+lbVNqoCe4N8bnYzHNwjodQEj\/9hG5B\/93eZUZ1Pk+yKwV65ciStXrtTlMe9dR1RYWChRU1Or0oFQ0km9iCCVYlZSJdAzKKMPyGryGfuLVBMMcTiP8+erSnrqPCNoOQRWUS+5+5xR+fYlyH9qpPV7HQZajCcqAbSa4QeNye6QN5DKrYyYampIzyUxMREzZszA1q1b0br1n+cRidKfZ0sM9XRkQHc9dh\/hCZlQggh2FvrYNq0zUy\/kwUiLCYe+7GSlmkBi8gSlj69DpV1vds+oqBXOO7yDoe9lNGnXi69PdTLWDITO3BB2LSdsAasjKS9nRrXVTD9khTigsdgABt6XICl+K+PNcGqjsvA13+cfnQgO9OjoaFhYWMi8G7072bQLuMCcCTvYsb\/vU0QnElIkkyw68G3P3MvF2m8eItyuK2b6HsPbk95w9duDl9YJzFWsrdBApmf0gZ\/hUVZ18VHA9XJnaM7wZkaSfGwyfkXJiXj5m6RTPWV9aZ2y9DsoTDoNTeuNyN4xC0qN1aC79JB0Yn7zcKgtuZryaoNTK9RX+fP7\/FDr6jIqAp3AjkQkvIq8oOWuxfrM3ZwLiZoE9rBHBCL456SkpMDBwQGZmZngtj\/Ozs5YtWqVDGyiBRGJkp32v\/vBey6kIyWnAAFWncENwmZVEEzn5TGXsT7FRGKEs8UnoeLtXGdDSv2XpN9lRo+M5rO1\/aFpF8gmg+lpCZg6ajbQGrrLZA0jqaAM976oePsSxiGPeL+8ri6jItC94MWEzaTIBDruOqB\/5pvN4aPmw6AQAk+gx8XFoaCggGEXEhKCXbt2wc7ODh07duShE00Kuij5wqUfmjZRZhc9g8LxRbC0wzlznbA\/LBQfLl6BE0sXQClPCbqOulC9pcruZ0ZnosSihH0WbxFDY68G+\/x2\/FsmDeuVV2FRjjPsJo7CvTyyAQC3dMl72bzJC9H33rHrXTQbY8fCKRA3VUV28i2svlaJW7dusXsbnKZi9urP2ecXj+7AftqnuJdfwb77+fnB2toaTMrsZmOs1hv8lFnK7s+bN49JmSJXsqioCO7u7nj69CnCw8OhqanJCxk3xlt5tzDGcQx73z6D+0ANamjRogU2b96MbLVsdIvpBp3VUtU8fvx4zJ07Fzt37sSPP\/7Irunr66N79+4oLS1Ffn4+kpKS2HXRtH1Jkv12HzDQueXRtGN\/XBWPgUbyCeDuSRh6fILL1iEyy0v9W3WId4iRFZGFRrmNYDDDgE1ChVYF9Bz0kL8oH5MsLSFeo4OiO\/Hw9VqPc1lK2LFjByIiIpBx6SscCtuFoONn8Sz3DeynW2Hs8EFYONAEAQeO446yCfYF+uH4ismIymyBA1HR0NLSwhy7WehWkQaX7uq4K+6Fefvi2ETSPVra3TWVsEj1Km4O8sLuI98iLGgz1MNnovkwRxnvqS6gj98yHkmxSbCIsIAvfFn\/vXr1YqCTQJDRpXdearoUDxwfoFOnTrh8+TKePXvGS\/XGjRvx888\/48KFCzh48CB69+4N0eIjdyTbrbuySjExMQwUy\/l+CLxWjkYv7qP5T1uhtE8Daca\/8mAWWBdAOUWZfX\/hL3UdCXQqo6JHYa\/FXvaZ6vQap4dZbSvgHH0Vj3ML2cD9\/f3xkWkbZhi\/aDIMYcdiWf3P7GZg2fTx2LhwNg4\/LEO39m3gZ1aKLpvi8PpkEK49es5ADhqgBjPVV9BYdBie4V9CUlqI8a9OYf0NYJhmARwN8nH7nQpcUwzg3zETPVqUgDOyr6EGR0dHfhVx6BgYGGDp0qVwc3Njk2hqaoqhjkORbpGOwFWBsC6yZiuDCoEeHByM2MRYJIQnYLbmbHTe0hkPHjzgJZvUC60g8gwXLFjAVlFxcbFU0r2O35R4TOnBg+7nH4SXvRfhhPsnuHH1Cj53c2azKQSWn0Yybn45KJxQyFZB89jm\/K2a2visc4NSUzFWr14Na91XmKb7Gp4FAzFgxFimDq7vXIYFO44jp6wx68\/c3ByhO4NwfPlEHEjKhW+HLJiol6O4UoSAVKlxs9J5hQ2P9WCrn48JXbSQMSkEts6LeXUmHDN9rk3SCfTRjqPxyOIRJq6aiF1Fu3jQPTw8QBIcGysVFq6QpA8dOpQZUSHoy5cvx6+\/\/oqjR49CLBZDtPzoPYmbZVssPnIfI5r8H\/y81qPRCHfkNuvIS7q82igZINXjXNFTke5QaWM0PmY8vlv9HYrNi5G7LhcfLjaHvWEpbJet5\/107oU5qaHvJH3kpgktPVePXo708+DBg9mS5nSusB\/Spy4uLli0aBHT8TW5f3UF3crRCreNbjP7dAqncNL9JC\/pm4I34UDiAWSHZ+Os5lnGpsobUk7SaYW8fPkSCQkJUklXpNOHjR6LuMbD8e5KNIpuxTJJLzUvZdIsKhTh1apX0LHXQUm\/EpS3LodmlgjNbmrgyMb9MDHrCyN3I9b5S4+XWOThwnT6RqepuKj8Aa\/TT5w4wUsDGR6SegJ28eLFMsuYdomkkgjMkSNHykyOEFhOp9cVdHnJp+\/yE7VlyxYciD2AhxEPWfVODp3Qtldb5r343\/JnKtUs2gz7TffLjIvaCSVdqJZI3ciATh2HxJxGwGpn9hBbh7mIighD8cFyPB+YzryXQa6DcG3cNShnKMNnpA\/MWpih3LAcPpd9EG8SD+2V2miU3ogZWMMOBriaeVGh9yL0aUmfamtrw8jIiOlLTvI574U8A7pOA87Ly5PRyZzUc\/01JOjClVYwuIBhUtmikvfTp8dMR+LqRHZdOEZu8ug6OQ0kYNyqZqCTeiGfnCtCP52u0e40oetOJHwg3TmSXzoYgxUJCtKQxnatHJfuUbQcnhJPnkEUsowKO2iAi0K6gDZFOi6hyAqYVoWXr4ukC+twmyT6S+9Pm0BhAKc+Qxf9ePelhLgVrnA70m8X9EAzlUYskDHUVBMug6UqQ34XJv8wbjNh+YsyQrZKdb1SCx202XqjQXkSRS\/JAd5IyxDiccuQd8yLbaDI86Fdbm2BjNrsQH2AramuKDc3V1LWuCkzpME2naHbXAUk7Wu+keoxh376bHdaHz69ZQEw5FcVBPpLDS6RWYZ+DUtOUb\/CyJPwJZU1W0PbaQ9exe1ioBNvw9EPNRFk1dEANcURSNrry8GI+vglSihCNLKLJrhNEvnrZNiK2w6Hesd+iPLugb2tduD2m9u40eJGjRMufi3GqJcGKM15CPuKEKyx3wTfDtlQGemC0a6BDSUsMhSBUG1xW35NK0+8OrkNrewC8NzHEk069quVT1cEOhcxq27gxL3Q6ucKGdG9e6X7FLJVpNOJAqDrVMg7E+2OfyqZ1U+fj4caK+fB19cXjm7esIt+AnFZDlp23wuX0U7MWG41C4Tx7FloLNh1sd5MTLD7+2lInpGM9vOKoP8mBd+90cW+viEQRczFpWwJJu44w7bbDVGqixIpWgHV8enyAQvaSRKfvmTJEvTtK+WjiFUk7oWC7PKFVj\/d42zc3bt32aZr+vTprCr1d+rUKcyaNQsnT57EjRs32D1epwsNKIXslh55gM76zXA6KRV53dxha2CDoIfRsFOJxpyZA2H0JqsqdkOGYIulJV6Yn8eYxDQs7VyJ+I9Po3DndAZ6u8+2VaFM\/8gEEKvIkWP17eevClgoGpfo8p2nkgFd2zA9fj45j6mYLxKes7oddNThHnkBPT89A11tLbz+JhgrvNR5wFNNpHyySaqUZ680MoLjsGHovLAcrc+cRIxBNxwZewx5\/pPwVa4YusNns41LQxSOgVQUrqsLr\/5XBSzqBPrK0Sbw+vYRLj9+xShUbeUiJPmMxpIcZ7QJD4engFIv0NVD2ydP8EvXLgz4vBYtkBwXh7vdE3Dz3CYUp5gjYKQX8sNcsK\/CAmYfj2ow0BW9TH0m4q\/yVBSCzrmMnHppUvAc0VtXoqjTeJSbDGJtBnXUQPrMuZiw7YYM6F8bqyCurzb0umYg5ZwYuVktcCjhBnZo7kBsViSOLnrC2mssiobXl4lMn8lHZP6I1MvHQrm+5CWdjGuz\/lYy4cKGAL06oync+NGYOIqZG18VQ5r59DE8lztjhvVk9JvkhOXHHiDQqhMkxlcQeGMCzg9l4VCsnK+Me5F6mKaSj9l5BezaBTsTHLCXZn9NfjMZJ5acwJo1a9izyDgHBAQ0qCHlkpHyv9vGg0oAq7TpKkPjKjK6fxR0ak90LccVEdCRkZFYtmwZtm3bJhO4oMkh3ogEjna5zGV88sspaKd+j69jonDozA0c2O4LlBWgTLc7SjpYQnxtFyqK30r1tmsO0qcXwchWF6L70mCGxOaVdBKTVPE0KhstQ1tCJ1EHwS7BoHCVnnIpKpWU4Nv+OfRUK7D9XU\/8mJLHmrh91AKjKm\/x1Cv50bTJ+XnlSKw4k8kzjdx2nwMrak8wjC9sxsX2ttixMwQBY9riqdlM7NgehKn6RQi++Yb1v379elw7E4sfEpJg9kEX7I\/8QmHAor4rjgN57dq1jJ4QFiHI8v3S+FkKBsejk0+Z8PgVfFYtwNqVyzB23DiMnrEAoydMhe+8CVjToy2+QjnEK7PQ5ypw+aQeHMry8UHbAgQ0EyOloiW2RgVDO1kbs1bMwoaIDZhsOBkr5s5E6bP7CPoyHjvDIngy6Pq+tXwQwuTnIF5Cnx\/ywIaD30O18wAsM3iOyPvFSHrViPHTycnJjGmM2h8K4x89cbrUGGEXk+HTIQvpZrPgufcIlhvlYLhWIfxTtXA+vznj1Ns7BWBB0GHGkZB01ibplLFJ7gG5ClWdRSmU8vwR55PTPaHqEaoXiizxki5+fBLHog8iKf0tNq6Yj1UrlqKwzQCExKdjZd8miPRzZQFXdfUyfPUgC89aAwub6cHMLB+jVApw9aYYma9bIvz6dTagCY4T0NSmKXb23Ak768mw\/1AbU7ce59m4lYtccG\/9aLhfq4DtfFd8YqzMko90lsbgqf80+CerIu7WE4zp1wPLDbPR1vMHFvMUgkWBkFCXTxCVrobtMyyQcDwUUc\/F2OvjDvOpC9iLc2wfjUlIH9cGOmVtcqBLLVPNhYg4Ly8veHp6VlGhJNRUevbsyRhTUUDsbUmr1\/8Hn63bsHZLCEx11THbzgGTZjnhyPPWULoQCN0mRdgc6YN9nuNw+7aYgU7FQU8Pi\/LzYV1QgC1iMRJbSkGnDZDtFlukZKbA6oOJ+DJ4P3yM01ibNY8MkFPaSOYNSBIWDO7A0uz01p7Gi92OSG4\/AY7ee\/h640YPx5bAYD5Mxqkb4SqlzQcXDuR2ge8LOifpNABFoHNAci4w6WqKJk2ZMoWNQQg+jSEtLQ2pqalMt7MUjNK06\/Ddug197L3gNrodbGxnI721JSRNdVi4jvRg6scpeNprIY6W1A30DYkb8PWKr9Fcuzl6GvXkKVt760noXpkGz8g4Bi6X48jlNnKgc9fTTx+Es5sH7r1RYkFoY2NjmUAGA317EHzbPcPtjDwm6VxkKTRDjFtvVeHTIRtEgnnkdEO\/gYPqpF5qk2x5FULf5Wlm0gxUSKVRGPDzzz9nkyG6+ihH8sV2b8T\/dAVGk9fB59MOmDXbHv0n2KFf3968fr875hauDo1BZmU5dr+VSvq6CqmkDygpwSFDQxCzTJJeemA+9lUcw8VoXTx+rcpH7KmN76oluBz3DcJDgnE9YiMWn8tjg+2GZ7x6Sd48lamdgWOnsBUQH+bH19NQLoeDvT0WzpuLj37djS3XC3HnnSoD+m6BCgP9849EGBCcAP\/QqPdWL3UB\/X3riPpY9JeMHD4UsafOMNDnDDTAmiUuWL1iKQx6joLHJl9kJX7F+m\/3cTs8vvQYOtE6yDbPxvDeEpwrUsL4ykqo2diwwCyXzkBE0dUtV\/HTt5exq58Wen3mx7JvyUh6hB7H2afSIK37UBN8Fn6OpT9zrh753\/eTrmHlhTxmR6gIjRFnpIhQGjJkCOLj4xmxROpl+xZveOvfZzHUP1vS3xt0+c1RR52mPK3LdUoMZIyTOYgrp1LTgQCuDYFOQQ0ig\/olq7F8c0qBU1JvyeecywccOL67uuvUd027Tu7ZjTT0oNK+D5oNmc02RIp899oM6fsCWpd2zJDOG9VF5tQF15ACGqGXM3jKl1g17thLbZ3TBBEPrYidq61tXe\/Lg8kmZMMINDYwhWr7XjKZYvJ8+t8KOsenc8EK4QvLg04gEph1AZLqkZS\/b0irNuAV7TLZCvEagfJX2dBxCkF26HyW3154Kw5vz4XLZPn+raArOnPEBTHoxSnnRJh21s+iH5p92QyWnSzRrVs3lJWV4cyZM3iQ\/IDh1LtXbxgONGRSfjjxMPM0qHDpb7WBWZ\/73DEYji\/nvouatYJRQBLenA3n89R1XI\/\/c7gXedBplyXkSchokTEjl2fdunXM3UoTpcGu1A4Ohg5QSVBB0cMijJ0\/Ft+XfY\/9j\/Zj\/ov5WGK+BK6urmxDkp2djaioKMyfPx\/Nm\/+ekFQfgKurW3T3Il79dgBApNoUWrZb0diwS61dKwpY1NpIUKG6gIWVlRWOHTsGS0tLPsf98OHDTHi7du3K8hpZjFQ+miNMIXByckJWVhYmT57M\/MxBgwbhww8\/hOMwR7ZJds52xqhMU7YhemL8BNcDrsPExIT505s2bWIJl3\/2yYb6gPVPqFvl+AttZ0lCOXaQJNXe3p7NEvEGlAxJYB4aeQg\/6Pzwn3fYgCu\/jEHTpk1ZyhgRQNQ2NzeXr\/9nn2z4JwBZnzFUAZ2knFQK6XGK61GmKYFMZWn+UlwZcwUW9y3QQ6MH9ujRNn0DTsXZIDZWC5mZ0XBz68MoTOqHmyRFJxvqM8j\/troKQScdpKOjwzKurl27hv79+zN14ZDqgPs20hMOesWnkKWqwlTMp0ka+FpDA4vychHcU5pDQ6ATrfro0aNqkzj\/28Cs6\/tU0ekcWGQ0+\/Tpw\/IKCXQuj\/CUjR6Kre9De9d3eHG1KWYbPcbBtu1gaLsOI5Sf8cdBhOlv\/0uSTu8dGhrKcjKFPDt5hCS4pAWqSPr+\/ft5dSISiVjaL4EeFhYGWgFLWrbEuylTEHb2Ebo3z0WWahbWGQ6EaZMmuGFmJuMeUl\/e3t7\/Svpvuf\/Vgi5vSIXuI62CaW5uGGO5GY4jpMdlJppNRJaqJZbkB2KzmhpUVaXRJE7FCFObKWngt8SBWlfikN9\/54GvGy\/9nYc6lfdpT8kNvyU48OPfvXs3Ll26xASH6FlKwqLC7TvkgxXk7XGSzp3W4DiicePGsT5qdRmFqoFzJcm7adeuHXv4oQ6HcKjjIZy8ewfD2rWrEXSv\/1A3G2o\/oMf6VfTbHCLBUdbakH+f9jQ2T8HvTAjjoPLhOQKbXOebN2\/KqBJOvdja2vJULqkZUtMUmGeg1+enRzjQyZvhToud0j0FT2NP9C3uy3ah\/22STiAJDxkIJ5uknQon\/cSEcpI+duxYdpyIi6FyOr1aSa9OijjQZ86cic6dpenVGR0y4N7XHVPeTUGUclSNoNcmnf+0+1zEhzu9R9F+RYFoGjdFjnx8fDBp0iRGiwglnTaOXLC63qBXl4r25PETrCpeBW9lbygrS3W9Ip3+TwO1tvEIQae6Qk6Kvu\/Zswek87nDC0JJJ++F0+lU9+OPP2bHc+jESL3UC5Psan4lwtDQEPRfWP5OJq82QP\/O+\/UGvT6D\/Rd0xWj9C3p9pKiB6v4LegMBWZ9uFPrpNQUx5IMSwhNo9GBhAPl\/cUdaJxpAnk+vLYghnxBKZ0ApiYaCGxzfYmNjw3M1ZMH\/l7iX6iReyL38P+7pnz3t12QWAAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:66e85348]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:4c901919]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:39c29f9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51561f91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:55b61ecf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:153dd0c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8cd5c8a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b5d7fac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0c72d7e0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0e92f4b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13bf07f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3b526546]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:09d2e8cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61d5c65d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:854e323f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93c842ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4292f193]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7fd07625]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6311627c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4333cf51]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:412c8e4c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2fdfec66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:107dc2fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c4db4ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2079213c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6eb425fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8673b863]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:019c3e9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1b97fc30]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05a9a300]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8ad946a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30f38b1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3d6d6bc5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:25020dcc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a90955c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b6c7968]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6afba108]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:27211ad6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5643eb4a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:84dfb973]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5511e390]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f4fd056]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:33daa2f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5c87ca2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:70356c85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37765ff1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:83f69a0c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30b6f92a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:645802aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3023b4ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:332f588f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18ad7d9b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37f52a9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:62c2f4cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3a51b7a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:847d426f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:53643fab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:514002c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69ca112e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:84fa34af]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c9eca56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51ebd9dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35697a8f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:304accda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:42676b03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8cab75b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:823c3cef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6aa20f99]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7f7ad05e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95dd3faf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51e6399a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6984b449]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4a0f24bc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90f9bde8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7988a91f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34d47b0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9fc1ddc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:29d9ebf3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:51946b13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:73a7e3d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7228ae3b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:10ec561a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d06daef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ed4dcee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8fce9b6b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:38a36ee6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05038c20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7be4cc95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96a47bbf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:78564df0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:824d1e64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:40846abf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:88e4cfa1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c653a11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3bef9a26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c70dc28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:53c32b96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:54b107ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:198be025]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:33e1fcda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:888dab7f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49db6cb5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:649ffd32]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"774×12 table","name":"APP_result_MK","rows":774,"type":"table","value":[["'APP'","2025","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","0","-0.0227","0.1021","-0.1486"],["'APP'","2025","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[0;0;0;0;0]","[-0.0478;0.1573;-0.2912;0.1578;NaN]","[0.2739;0.5934;0.0330;0.4384;NaN]","[-0.3664;-0.2713;-0.6057;-0.1178;NaN]"],["'APP'","2025","10","'daily'","'BsG0_homo1'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'APP'","2024","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","-1","-0.0617","0.0677","-0.1920"],["'APP'","2024","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[0;0;0;0;0]","[0.0136;0.1902;-0.2479;0.0084;NaN]","[0.3422;0.6510;0.0964;0.2981;NaN]","[-0.3176;-0.2668;-0.5913;-0.2760;NaN]"],["'APP'","2024","10","'daily'","'BsG0_homo1'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'APP'","2023","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","95","-0.2088","-0.0746","-0.3414"],["'APP'","2023","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[-1;-1;-1;-1;-1]","[-0.2107;-0.1642;-0.1439;-0.0516;NaN]","[0.1329;0.2960;0.2053;0.2303;NaN]","[-0.5443;-0.6233;-0.4986;-0.3388;NaN]"],["'APP'","2023","10","'daily'","'BsG0_homo1'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'APP'","2022","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","95","-0.3428","-0.2134","-0.4718"],["'APP'","2022","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[-1;95;95;-1;95]","[-0.2341;-0.4979;-0.5698;-0.0247;NaN]","[0.0920;-0.0528;-0.2139;0.2376;NaN]","[-0.5515;-0.9432;-0.9234;-0.2860;NaN]"],["'APP'","2022","10","'daily'","'BsG0_homo1'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'APP'","2021","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","95","-0.4428","-0.3117","-0.5738"],["'APP'","2021","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[95;95;95;-1;95]","[-0.3975;-0.4450;-0.5352;-0.1785;NaN]","[-0.0596;0.0055;-0.1652;0.0844;NaN]","[-0.7388;-0.8990;-0.8976;-0.4402;NaN]"]]}}
%---
%[output:6ce2d67e]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"9×19 table","name":"APP_result_LMSlog","rows":9,"type":"table","value":[["'APP'","2025","10","'month'","'BsG0_homo1'","'neph'","'log'","'LMS'","1.3906","0","-0.0147","0.0064","-0.0358","-0.5349","0.2344","-1.3043","-0.1367","-0.1362","-0.1372"],["'APP'","2025","10","'month'","'BbsG0_homo1'","'neph'","'log'","'LMS'","2.9420","95","-0.0334","-0.0107","-0.0561","-3.5651","-1.1415","-5.9888","-0.2841","-0.2836","-0.2845"],["'APP'","2025","10","'month'","'BaG0_S11S12'","'abs'","'log'","'LMS'","0.8603","0","-0.0092","0.0122","-0.0306","-1.6976","2.2489","-5.6441","-0.0880","-0.0874","-0.0885"],["'APP'","2025","10","'month'","'expS_bg'","'neph'","'log'","'LMS'","6.8667","95","-0.0134","-0.0095","-0.0174","-2.0831","-1.4764","-2.6898","-0.1258","-0.1257","-0.1259"],["'APP'","2025","10","'month'","'expA_bg'","'abs'","'log'","'LMS'","3.7175","95","0.0147","0.0227","0.0068","6.2984","9.6869","2.9099","0.1587","0.1589","0.1584"],["'APP'","2025","10","'month'","'SSAB'","'abs+neph'","'log'","'LMS'","1.5648","0","-8.3190e-04","2.3138e-04","-0.0019","-0.9483","0.2638","-2.1604","-0.0083","-0.0083","-0.0083"],["'APP'","2025","10","'month'","'SSAG'","'abs+neph'","'log'","'LMS'","0.1203","0","8.3079e-05","0.0015","-0.0013","0.0819","1.4444","-1.2805","8.3113e-04","8.6899e-04","7.9327e-04"],["'APP'","2025","10","'month'","'SSAR'","'abs+neph'","'log'","'LMS'","3.4559","95","0.0028","0.0045","0.0012","2.3100","3.6469","0.9731","0.0289","0.0289","0.0288"],["'APP'","2025","10","'month'","'BbsFG0_homo'","'neph'","'log'","'LMS'","3.7526","95","-0.0117","-0.0054","-0.0179","-0.6458","-0.3016","-0.9899","-0.1101","-0.1099","-0.1102"]]}}
%---
%[output:14a83052]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"9×19 table","name":"APP_resultLMSlin","rows":9,"type":"table","value":[["'APP'","2025","10","'month'","'BsG0_homo1'","'neph'","'lin'","'LMS'","0.7031","0","-0.1546","0.2851","-0.5942","-0.9900","1.8262","-3.8062","-2.5455","-2.5335","-2.5575"],["'APP'","2025","10","'month'","'BbsG0_homo1'","'neph'","'lin'","'LMS'","2.6396","95","-0.0645","-0.0156","-0.1133","-2.5256","-0.6119","-4.4392","-1.6448","-1.6435","-1.6462"],["'APP'","2025","10","'month'","'BaG0_S11S12'","'abs'","'lin'","'LMS'","0.7021","0","-0.0117","0.0216","-0.0450","-0.6790","1.2553","-2.6134","-1.1168","-1.1159","-1.1177"],["'APP'","2025","10","'month'","'expS_bg'","'neph'","'lin'","'LMS'","6.8952","95","-0.0249","-0.0177","-0.0322","-1.3082","-0.9287","-1.6876","-1.2494","-1.2492","-1.2496"],["'APP'","2025","10","'month'","'expA_bg'","'abs'","'lin'","'LMS'","4.0235","95","0.0153","0.0228","0.0077","1.2072","1.8073","0.6071","-0.8475","-0.8473","-0.8477"],["'APP'","2025","10","'month'","'SSAB'","'abs+neph'","'lin'","'LMS'","1.5402","0","-7.5543e-04","2.2554e-04","-0.0017","-0.0825","0.0246","-0.1896","-1.0076","-1.0075","-1.0076"],["'APP'","2025","10","'month'","'SSAG'","'abs+neph'","'lin'","'LMS'","0.1349","0","8.4493e-05","0.0013","-0.0012","0.0094","0.1480","-0.1292","-0.9992","-0.9991","-0.9992"],["'APP'","2025","10","'month'","'SSAR'","'abs+neph'","'lin'","'LMS'","3.4718","95","0.0025","0.0040","0.0011","0.2880","0.4539","0.1221","-0.9745","-0.9745","-0.9746"],["'APP'","2025","10","'month'","'BbsFG0_homo'","'neph'","'lin'","'LMS'","3.7174","95","-0.0018","-8.1473e-04","-0.0027","-1.0728","-0.4956","-1.6499","-1.0176","-1.0176","-1.0177"]]}}
%---
%[output:929ca2c8]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAACt9JREFUeF7tXG9IVNsWXwMDGRh46\/nopqHRDQ1u\/6CH2tWsyLqXSszSTC9p\/7CstLLS1FBRS+8ryOSRlZZBmYF+Kijyg\/nopvUt6fnBJvLi0H0U\/fngLYOBef323D33zJlzzt6jR8des2Bomlmzzzq\/vc5aa\/\/22lqcTqeTAjKhCFj0QP\/48SMVFxfT1q1bKTY2dkKN+n+\/mCbob9++pZ07d9KTJ0+otbU1ALrJXuAFOjz83LlztHnzZjpy5AgVFRUFQB9v0Pn43NsDoJuMOBFZPnz44Jw6darHyHV1ddTY2Mg+mzFjBt24cYPmzZtHdrudOjo6aNOmTRQeHm6+NV\/JiJaenh6nMlHyBLpu3To6f\/68R3jp7e2lzMzMQJwfo3NYamtrnQghyrBSWFhI+\/bto5qamgDoYwRY6+eW\/Px8Z21tLfEQw72Zl+8bNmxgiRUS8HRzZsALdCW4eL9w4UJWq2\/ZsiUAujmYk254KSkpYZfYtWsXrVq1isrLywOgmwW6OpFevnyZLly4QHfu3KGRkRFav3497d+\/n3bs2OEGvaCggGJiYkwyYfIM43A4mDFWq9XQKD09VHQyVR1bHKEmr6ysZN48MDBAJ06coOfPn7MLz507l3k7wgtKxqNHj9KjR48mD1ImWALnevfuHXMyLgA+NDSUgoKC3J+J9OCMeInEa0WKZNnQ0EA9PT3st8qYjv8DeLy+JPntN5e1ERHeVsPZ7t69q3s7e\/bsoezsbOaUIr3S0lJ5T1de8d69e2z5jzAzZcoU5uVVVVW0Zs2aSYXz4KDLnMhIfbPu3yeqrCTCv1ygf+UK0YoV+Pw+rVy5UnhfFRUVhJdIurq6aAUGFoiXpwN0hJDm5ma\/gC4CUwQkv1+AbYQTgL96dSUD3izJycmhKxjYV9DV4SUuLo6Sk5NZTIcgifBEIhrcl+8fPLBSTY3VyysvXXJQfLwrwdXVBRkCCd2ff3YQxkpKMk6GLtvg5eaBjhFltic0Yzq4FiyYIODUly1bRhEREYwCePHihS9YfnW68Pbc3FxDZtaL8IKnIwO\/evWKLBYLy96gCfAYonzcu3cv+z4sLIwBapSklIjbbA4aHh6m778P9ijJuruJWlrMmxvEbB6i5EbdLqcmqTVnzhwhN+VFeL18+ZI2btzI4vqiRYvo5MmTbEUK8isvL4\/NIjY2RkZiDZMUt1EUg5HHTAyrktAo1eYQ0Z9ZeRS\/Vv4ESTQkJES4B6G5IsVCqK+vj40HgCEIOaB08fhERbV+Lp\/0t\/CQS3JyXJWDUTLDdxJFwRihEP3cYqgQGRnJ7llUvUAPSbS+vl4Muh7hxS2BZycmJrpBz8g4T7\/\/fkN0JwxMfwGKxItkKhLoZWU1sZCpJVhdZmVlMcBbWlqEetXV1ZSWluY76MqLc2595syZ9PjxYxZekpOnsdAyWSUmZuRzznlHmZnfGpoYHu6gX355TbGxI4zegIcqV9oAHLkLTzcXkR5CC6KEaLfNK7yoLcUuEpa\/g4OD9OOPuZSR4R\/ArVY7ORzGu1XQqagYpB9+cND161ZqaorXBB56q1fbqajor2U\/V8RqW8SfoGR++vQpzZ8\/n61luCDEHDp0SAy6mvDCbHZ3d7MfKvdJMXB6+jEaGvr3hDs5QAoNPcpAf\/36n7pABgd30Dff1Lu\/HxmJ+cypFHg8mRgLOtA1W\/BkoMoTerqa8Jo+fTrl5+fT7du3mU3Lly9n8YyLxTjvjOo+qqvtVFam7cUoAZGU167tZWPbbOHU3BzutYgqLyeKjnbpaAkmDIAbiWz5q6fnE8uoNATeje06zqejZDxz5gxhMiBml3igKrq6XGWjFk8CMAG6WlCL4yVBdbjrdj2eRlTW8mvL6ok8T3NFijgO7gVbeOouL1xYxBFx7xRVL0ryiRsqAtOXG5fRFZW1suUv1xMBju+laQAl99LU5KC8PO3hsVD9XGWxchFRSUYPI5nJvWA8GZ4Gky7D0ciWv3hiZZ48S0pKihPtc9iABt+C98pNDNAA4JQR55HZZ8+eLTOZX60OJ7wQLbADB5k1axZbOKF3iHl6W1ubEywiDyOfPn1iZQ8mICEhwSO8oLKBLrK0zRamyZkg\/iYmemM+UdwLPI2HKH\/MfFvbVIaRUfOtm3u5efMmq8XRw5ienk6YKTCLykTKWzBkaQDctCiump2Y\/QG08poNDbcpM3OZuxjh3q3UsQwMDDjxBUB\/+PAhpaamMq4BDCMEzKKy7yU1tV6KBujsdNCvv1q\/AO7F3GkaGHhGb968YTQ4F1ApyoYuS19fn3PBggVu0HnjESgAMI3YsOWbGPB0eKYMDYBltt0u5j\/MvWX\/jsbLX6UVnErBngQvRizNzc1O8AU8vChXojzBIhHg8wcP7JSQMLkbR3n1IEMXy3DvYyl\/OfgI1RDu7Zbs7GwnuHJlIj179iydPn2aDh48yHgXbEyjydSfoMtyL1jdQoqLjTkijAddJN3GRm1dztFggXb1qljv0iVX34selcIbdS1Llixxvn\/\/nm1O8JlQljsoJcHAQVwl48R7+mi4l+HhTdI8zV8czcw\/nTOSUQZqjkakh6oOjCzCiBJDr5iu5tP5I6F3KGA8qo3x4l5keRoQfODM1c1G4Mexl8BFpCfNvahbpXEBXhri\/e7du+n48ePuC8vSAPiBzF7laLkXPr5R3ws32ohaQBOR0a4QFjWo5mT1ZFK51x4pJ7zQn15WVsb4F\/AwnPDCoFjeb9fZz+WJBw4iw9Hwxh+lsQBJBKaoP0ZmjeDXZiN1LyPiERLp4cOHKTg4mA4cOODVUiDDCspMDpKULyJabPGxZIisSdVspNX3otds1N7eTiUlF1niQb2PSghPBbb38Fgi+QQHnyabbYYbW6v1vxQS8i8qKAhhKzcI1wfvjc9yc9d6zYUMiTVZmo0MuRf14V3Ecz3QA4SX+JlE8y3O4BpyL1qg6\/HpSsKLNxuJzXBp6BFe+A6hwCjp+kJiySx4PG02t9mora2NkpKSjLkXNbWLZWtKSgoNDQ0x27Dxeu3aNVq6dOmoTmKIYrBMNSQ7saPTM7fZ6OLFi2LuRU3toqsL7dG8P115TN3Xg15yyczfHV7mNhupW6U1uRfeDcC5F9TlfI9UTUv6Aro\/PdjVey72e+hlZ7fQdp36l3d34YQKNudl9LSu6sW9aFG7ysGVS1gOOg4MiP4yxk8\/BUnduBga3zTQtVVa6iBc30jAgqK1mh8OwAJJ3WyEkxWowLigrjfSw5EZvITcixp05ZlS9aPBQd+2bRstXrxY954cjjBKS\/uHb2iZoI2ytbDwKevaam8PpuvXV2uOCr34eBsVFLz3+l622ai\/v5+ioqI8mo2io6PZdqeQe1GHFyXZDouUj4bsQS\/U2mY3JQUF9dK0aR3SJBZs90ezkdRBL3UihbFaHV48nMge9IqLM6\/9DmVgVZWdvvvOPq7NRmN92KQJL14y6lG7alpS1jBZNtKVzMRcjhZdIMPRyNo7kXq6f05qrEbIVC+eJ9186\/Aaq33+\/P24gY6bGi3h9aV6sOxEjhl01Pecb8fCitPAvNJBMsMOTH\/\/39028V2ZwsK\/uXerlBz+qVOn3Ju4sjdilp7ofnAdbh+v7m7dusUurxeivZqNxvInAp89e+bRF4NKB2eWjh07xrh4dRNqa+tD+uOP\/3i0YWdkZHhwFTBe3bRqFqCicXy9n87OTtYrpNzMx\/2Imo3+B5f\/hMLpEiugAAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:5e50f5b2]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAADrNJREFUeF7tXA1sVFUW\/oYM0KZFKQE1OIYSQNEIkkBC0SpDKKs1CxogAaRIsYJBqk3Eyk\/RUmihtSkKaAoWLKiUViBmI1GWrqFgra0xBGF1BRo6sBMwlmmJYmjDbLr97vRO38y8n\/umRXDdmxDaznn35zvnnXvOd88dR0dHRweiaTU1QF4ewP9lS0wEyssBtxtVVVWoq6vDm08+if6FhSFyl\/r3x3927IArLS1yZI8n8Df29T\/aHFGBTrDXrjWGhMCnpweUoiLHniyUGDHYH1g59kEnOFOmWNsgwTYDXPZw5Ahw9OiNU471TH93CfugE3CtS+nplOlGpNWa9WVXOT2d1w183h7oBGf48Bs4HZOu7SjH7Q7pyNOl1MRbZJ\/QBb2+vh7PPPOMmPjGjRsxZ86cwCJuJuiKqvbPnw\/\/jh1Cura2FgUFBajRvJkEvqysDMnJyRE9er1O8TeXy286mt\/vx9WrVxEfHw+nM\/AMG3\/W\/m7USQToLS0tWL58OVavXi2e2bBhA0pKSjBo0CBQGUmTJiku\/+aJNZ07h82bN4t\/Rq24uBizZs0SH9fXx3TKJqChISYoTuDffLMZSUltEV2cP+9Ac3MzxowZgP79+wc\/HzhwIBISEiwXHgE6gS0qKsLOnTsRGxuLlStXYt68eUhKShKgcxNNaouciOVIv6PAodJSpC5dajni559\/3vkWuFFU1A12+EOlpW0iEGPjC8PYIFw5ZWV+RsnRWzqB3bt3LwoZWwMC9Icffli4GH427w9g6V2RviXodDUqe7jcCsxkpUxTU5PluLqWbgb6pD8A6JarvoECKrmmbfcyhe4lKal72uEZqVGEYCVnFIayP9mnVahKyxX7vaqtcxmhkU5P9MFpNjUdsezC9kbKqCYrKwsTJ07s7vz8eWDYMNPBuONfOXEC8Q8+iJgYHR\/alSBp9wuv0wnv888DCxcG+t69G0nbtumOI2RTUoDcXNh7G6NjQYwWq0KqWIaMFRUVQcv2er3Izs5GQ0ODpTZ7IuDy+0EQ9drEtjZktbaGbOaU3ZyQgAPx8eKRS5cuoU1ps6eVW1umnbXQpVulA\/aSIwAEnv9uhWaknMbGRixYsMB0ioFEqRweT++5Fw7YK5YekhzpLEOVd1KV6y1l7tq1C4sWLdLtjoCnp6dj8uRcSxpJJXKRgzBsJFth1Wz5dG1nqqSgqpzsW1U5KnLMRPPy8iIy0tzcXAE6265dgIFuhJug2OTJ1hyfhtW2why2ohfZmypjqyrHflWVoyqnXTkVxIjG7dbn6I36zM0NgK6qHMqrNFtxugRHhdnNyPBi506X5Ryqq\/346iunKQvMjC8tzS8yRzO2WMrJQWtrnSgocEacs1AuOTmSXyH3smtXDSoqloguxowZg9LSUkGBfPPNN5g7txStrVlobx8XXFd8\/BXExJQgPv4AFi9ejFWrVonPmNVv375d\/Dx06FCUl5dj1KhR4ndboHMDveeeeyyB\/DMLVFZWYsaMGSH0STgejqeffrrju+++w\/Tp00Xqz5+ppXfffVeQXa2trXjppZdE2Fhb68Wjj\/7jlsaUFsyUYc8e\/ZBTO\/n58wNvmVkuxf4ox8Y3R69ffk45MowDBgzAtGnTgqShtG7tuI7KysoOrWbuvfdevPzyy7h8+TJOnz6NcePGCfKLr1gAdGuXcbO1MnFiWwgpZTQfMomSzjWbc0XFJdEfmUijVlzcjIwMp8DN5\/MFqXHKv\/DCC1ixYkXwUcfXX3\/dQSvmQTI3GwJOgYMHDwqhNWvW4LnnnhM\/k\/CaNElDAdxsdH+n8VWVw\/0pJSX0Dbt27VoIaSh8+pkzZzr4CsjTe7oYUrrk1TMyMoQCJNdC0LmJtrX9+YBX0S8jHZ7Jhze6azZp7bZBnzlzMy5d2ms6B6fTi5QULw4dMlcO5dj8\/t5zWSNHetHY2Hv9qYCtlWFGSuM8evSoAFnPeCPci9SGnrDkXr74woXm5mLd+RBIhk8JCZtx9eosS7nY2HolJbLPK1eyLJU9ZEi2COtu5NsY4DIBDyLj\/q1bDyIz868hIWOETw\/fSKUr0QOdA0nuhdbEODycsWWCMHp0fRAcFTkqcc0afeuUWeHjj9dDVY5jLlhgbu120nu5GDdqkIs88H\/ZCPwilKOmiyJW4l5kyBiuDSPQ9UyNIZcVsyYsw0ROJStkH6py0aT3ehYs15XuyUMujAusCLzHna7GvYTH6dxE2eSuK89H7fq2nsj3hhJVFUTl7F6kb8F5ieVITHdj4bAaJC6yLrDylB8R8lYtIk6ne5FWzkRJy6dbdWb6uQpDBZ68NOHYsWNYKA8uejRo18MeT+At0+NeVAii3bvVCqyMwpewNejG6Vu2bMHs2bPx6quvhoSMgUjDL\/7Jtn\/\/frzxxhvi13CugkweDx14wHDHDz8En+GhQ\/aQIXhw2TKRubGR23j22WfhcDgwdepUvPPOOz2G21lbC2dBQUSRq7+sDP7kZIjPp03r8TghHXQ5dVPuxU6czk10z549eOyxx3DXXXeZTvbEiRM4XFKCTXV1iO2aSFVcHFbdcYd4ruLiRYzIz8fVrtoTHuh++eWXaG9vF5u1ERcuB5UHKS6X\/oY5sLPuJcGk7qW5uBgD9u9HTC+fgh3cuhVTMzLMuRc7oMvKrw8++AATJkwwBF1a\/\/Rff0WhzydAP+t0YsOgQSi5fBln+vbF3gEDsM7nwz+3bEH7Aw\/g+++\/F3TDqVOnROXUL7\/8AgLKaiy3pkyOHPnatWtDjgwj5GpqEJOa2rsWrNhbVWWlNfcSTgOYxekSdFZO7du3D6+99hruvvtuvP766yJrpUuim7hw4QL6vvUW6i5cCIJOK6+LjRW\/X3M4sHzwYKz0+TBo1iyUDB+Ot99+W2TCdC\/9+vUT+8r169fFUgl8WlqaiH0JuFGTcgJwq8oBRRBtibndOPvee9bci16cLoqK5s0T45EjliV2EvRNmzbh9ttvF9zxzz\/\/LMruSIzJ5vR6UZeaGgSZlk7QPf36YUVrK1r69EHGnXdiRUsLktrb4dCsrE+fPkFL1y6YFQhmZXJSdl9xMWZnZ9vCqjeE\/S4XuFfEPPFESHe63Et4nC5rGZctWyYsmLV677\/\/frCWkSUY0r0QbJ68r1+\/XlhpsHk8+NuECcqgsw7YTqWKGUjzk5PxUW1tb+DY3YfbDRamOhcv1u+3K4Pz5+ToFpBGcC\/hNzGsahlJA7Pwsm\/fvjh+\/Lgoubt48SJSWHOiacczM9EUG4sZPh\/6dnTgeFxc8PfrDgc+GTwYf2lpwRC\/H\/rHx9HjpsM5Rd8Zn5QHpazNYWAf3rjndIW4jNgsuRc90I3K6uTJUbgLiIuLE9P47bffgtPh32j95JYZmZDgJydPvpkKI9kvP+sZIrfW0zw5Yt2nNmSM4F7sgM7lqda9nKiqwr8++ii4kfJZGTIOvX4d5T\/9hFEuF6Z0Hhhrror1CoL\/\/vBDuCzqXshbeNevh3fkSLgaG+HauTPy0lpuLupHj1aeE6MooxBW24lj+PDhoq7soYceEiHbmTNnxMHFuXPnhNzYsWNFzBlSv6g6DQUCZEpNTUiJhGrXRnIML4+w+ERhbJbgRbRwDkLPnVhNUpYQGMg5fD5fB197vg70zUuWLBGZodb5f\/zxxxgxYoTVUPqfWzBUjLtZlGrWZHGQWbjI5ynHU\/dgXK\/KjpkNbveOlULFUbAaQNalz5w5E\/n5+WBpGmNmFmPKg+noUNc8ZcBkqVRjsUBIVU53ntFaMItjVW4JBrQevEdrhlUQdFo2LWXYsGGGlwLYkV3uhc+sW7dOJE6MWRmGfvbZZ2JOPA6U3Mu2bdtE3E8aQDYmXjk5OcFqLP5dZqRMwszkrAyEnAu5l95sp8aPx5hvvxVdWta9yENpZqNmNzFUuZcDBw7orkUSZXqXobqKoS0xkFyNkaDR2OHyTODIz7B1FXEZjs3Eh\/JmjTK\/njyJw4cPW9e9FBYWdtDC5Q06qzhdmxwZTSI1NdXW5qhasOxJTMQhTSmD3vi0MDuXApTGTk+H\/5FHjJMjeoCyMhHPU+mWdS+Se5ELkFw6WUK2++67T3Dq8nYdQWeGahbNcNH333+\/pdUKNygKlq3vQzBjZRLVm+Elw7sPvV7TuxgivT97FjSk3M5IS++IgnNbkJwsoiaeB1jWvciQkQCwyosu5sUXX8TJkycFaAwZw2kAci8MMc0aLw+ovOp8tVUySOaBKpmrjJNVauiZWY9uaEC+ietY43Ih49gxQWeP9Hp1rxAw9mp0uQRm4VcadbmXaGiAp556CuPHjzcFnbw7L8+qNBXQVQAXb45wlStQVGR9pXHOnFJUVXUyl\/DA6PIOx5VyZIjIcWplO2Ob4MlpTk4Z8vOfj1iyEvfy\/9t1KqaiL0PKo1e5l0BZ3a1\/Yzp6yHr+ZNTci9GNafrJpUuXKl304jVuO40bql5Nfec3xmhoX0qZkcDaXkKfjJyLlO2Wk26OLoN7SGCkSDkZZobyjQG5pqZ0y3IU29dfVAkvJjraL0JQUQCv62ijA0Yq3QQBPxG3RE264sLljWUrll7KdsvRX1PxoZt2pJz+BAJySsVGet9spP0WjGhLMKLhVLSQaUPE7ptw3NbMQXe7A6DX1NgHnc9R8aFMkDroHDuqi14qFqkqQyLLyNrld68wrtXKhVublCsvD8h1gy7B764nDFAfAdBVZNPTm7B2rbZPHSckLnupyXHssK+a0YXK9j1SVcClnMoNt4Bldt+Ek9ZGwI1vwtEXU0WBSF\/WPEq2NpLZJbiB+4ZaWVUGWFVOBZ8eg07eRl5ukpy8NnvlJORd1B9\/\/FFcnPr000\/F3LQnKlqX9sorryAzM1N3\/qFsbQBIgqi9CdetcP0v1AuX1WOAWX08duwnqK6eH8zGWSYeXhF8220tnefIG0WlsnY9loSXinb0ZM6ePRvyJTySk2dpBg9Cwr+op7q6WvAi2rrtuXPnhnAVHEf7xT5mc1OteWQfKrJcz+rV72H79lUCaLP1ZGW9hdOn\/667HquLXv8FYb4y4DqXcmAAAAAASUVORK5CYII=","height":56,"width":93}}
%---
%[output:83b20349]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:20702d7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d1cd89b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4bdda39b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4375a328]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9fc16f2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3559a0a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13a73f18]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0a26820c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:057e9057]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4713ce00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5739c69a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:64b973cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:975d8141]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22c5c274]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6179cbfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:72c89fc3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:974c0935]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:715e5e51]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e757775]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f78cb31]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0217a16d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8209d2f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:489bd53b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03700d07]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3eec10cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:79ce1524]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:639e87fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1bf85273]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:539c218e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6d130db8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2f74aa87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:78edf7d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6438a825]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:41f5bdd7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7724eb22]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:85b5cbf9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:143e9611]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:548ffb1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a251bca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5fb69ed2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d32bbd9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6477dd42]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c081830]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:00564cc3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2483e4a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:964403f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:218196bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:836d1211]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:74c64c8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6350f5ab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87e6bf5b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80126126]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:70179f68]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:359fd545]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3832f91a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:67ec98b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c680544]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:662cec20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:48c661c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:734adad5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61faa441]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ce14ffa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58c7b668]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:356bc868]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6108f4fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7cd0802e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5e13cf95]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49f09336]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9ec7c045]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7bbc295d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06635f4c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c5545a4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82bf21cc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c474475]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:412e81f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22b26ff3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1489397c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e7d2b49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1cd9348e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22d9ec64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:07c308a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:78a4c413]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3559a220]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18d2a38f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2049bf38]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05fc8b66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:09c01672]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:251df1f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5eeaac3f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:15e66d17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19f9922f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:558013dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39a87371]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6d120068]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2bc3125a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b1231d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:156164f3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:389777bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:60cc417f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1db86f63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:09b60491]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:409e85e2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c30d1b1]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"702×12 table","name":"APP_result_MK","rows":702,"type":"table","value":[["'APP'","2025","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","0","-0.0227","0.1021","-0.1486"],["'APP'","2025","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[0;0;0;0;0]","[-0.0478;0.1573;-0.2912;0.1578;NaN]","[0.2739;0.5934;0.0330;0.4384;NaN]","[-0.3664;-0.2713;-0.6057;-0.1178;NaN]"],["'APP'","2025","10","'daily'","'BsG0_homo1'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'APP'","2024","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","-1","-0.0617","0.0677","-0.1920"],["'APP'","2024","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[0;0;0;0;0]","[0.0136;0.1902;-0.2479;0.0084;NaN]","[0.3422;0.6510;0.0964;0.2981;NaN]","[-0.3176;-0.2668;-0.5913;-0.2760;NaN]"],["'APP'","2024","10","'daily'","'BsG0_homo1'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'APP'","2023","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","95","-0.2088","-0.0746","-0.3414"],["'APP'","2023","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[-1;-1;-1;-1;-1]","[-0.2107;-0.1642;-0.1439;-0.0516;NaN]","[0.1329;0.2960;0.2053;0.2303;NaN]","[-0.5443;-0.6233;-0.4986;-0.3388;NaN]"],["'APP'","2023","10","'daily'","'BsG0_homo1'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'APP'","2022","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","95","-0.3428","-0.2134","-0.4718"],["'APP'","2022","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[-1;95;95;-1;95]","[-0.2341;-0.4979;-0.5698;-0.0247;NaN]","[0.0920;-0.0528;-0.2139;0.2376;NaN]","[-0.5515;-0.9432;-0.9234;-0.2860;NaN]"],["'APP'","2022","10","'daily'","'BsG0_homo1'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'APP'","2021","10","'daily'","'BsG0_homo1'","'neph'","'y'","'MK'","95","-0.4428","-0.3117","-0.5738"],["'APP'","2021","10","'daily'","'BsG0_homo1'","'neph'","'MetSea'","'MK'","[95;95;95;-1;95]","[-0.3975;-0.4450;-0.5352;-0.1785;NaN]","[-0.0596;0.0055;-0.1652;0.0844;NaN]","[-0.7388;-0.8990;-0.8976;-0.4402;NaN]"]]}}
%---
%[output:6e70fdf7]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"9×19 table","name":"APP_result_LMSlog","rows":9,"type":"table","value":[["'APP'","2025","10","'month'","'BsG0_homo1'","'neph'","'log'","'LMS'","1.3906","0","-0.0147","0.0064","-0.0358","-0.5349","0.2344","-1.3043","-0.1367","-0.1362","-0.1372"],["'APP'","2025","10","'month'","'BbsG0_homo1'","'neph'","'log'","'LMS'","2.9420","95","-0.0334","-0.0107","-0.0561","-3.5651","-1.1415","-5.9888","-0.2841","-0.2836","-0.2845"],["'APP'","2025","10","'month'","'BaG0_S11S12'","'abs'","'log'","'LMS'","0.8603","0","-0.0092","0.0122","-0.0306","-1.6976","2.2489","-5.6441","-0.0880","-0.0874","-0.0885"],["'APP'","2025","10","'month'","'expS_bg'","'neph'","'log'","'LMS'","4.5050","95","-0.0075","-0.0042","-0.0108","-1.1471","-0.6378","-1.6563","-0.0723","-0.0722","-0.0724"],["'APP'","2025","10","'month'","'expA_bg'","'abs'","'log'","'LMS'","3.7175","95","0.0147","0.0227","0.0068","6.2984","9.6869","2.9099","0.1587","0.1589","0.1584"],["'APP'","2025","10","'month'","'SSAB'","'abs+neph'","'log'","'LMS'","1.5648","0","-8.3190e-04","2.3138e-04","-0.0019","-0.9483","0.2638","-2.1604","-0.0083","-0.0083","-0.0083"],["'APP'","2025","10","'month'","'SSAG'","'abs+neph'","'log'","'LMS'","0.1203","0","8.3079e-05","0.0015","-0.0013","0.0819","1.4444","-1.2805","8.3113e-04","8.6899e-04","7.9327e-04"],["'APP'","2025","10","'month'","'SSAR'","'abs+neph'","'log'","'LMS'","3.4559","95","0.0028","0.0045","0.0012","2.3100","3.6469","0.9731","0.0289","0.0289","0.0288"],["'APP'","2025","10","'month'","'BbsFG0_homo'","'neph'","'log'","'LMS'","0.0433","0","9.5328e-05","0.0045","-0.0043","0.0053","0.2514","-0.2408","9.5373e-04","0.0011","8.3310e-04"]]}}
%---
%[output:44910ae2]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"9×19 table","name":"APP_resultLMSlin","rows":9,"type":"table","value":[["'APP'","2025","10","'month'","'BsG0_homo1'","'neph'","'lin'","'LMS'","0.7031","0","-0.1546","0.2851","-0.5942","-0.9900","1.8262","-3.8062","-2.5455","-2.5335","-2.5575"],["'APP'","2025","10","'month'","'BbsG0_homo1'","'neph'","'lin'","'LMS'","2.6396","95","-0.0645","-0.0156","-0.1133","-2.5256","-0.6119","-4.4392","-1.6448","-1.6435","-1.6462"],["'APP'","2025","10","'month'","'BaG0_S11S12'","'abs'","'lin'","'LMS'","0.7021","0","-0.0117","0.0216","-0.0450","-0.6790","1.2553","-2.6134","-1.1168","-1.1159","-1.1177"],["'APP'","2025","10","'month'","'expS_bg'","'neph'","'lin'","'LMS'","4.5234","95","-0.0142","-0.0079","-0.0205","-0.7384","-0.4119","-1.0649","-1.1420","-1.1418","-1.1422"],["'APP'","2025","10","'month'","'expA_bg'","'abs'","'lin'","'LMS'","4.0235","95","0.0153","0.0228","0.0077","1.2072","1.8073","0.6071","-0.8475","-0.8473","-0.8477"],["'APP'","2025","10","'month'","'SSAB'","'abs+neph'","'lin'","'LMS'","1.5402","0","-7.5543e-04","2.2554e-04","-0.0017","-0.0825","0.0246","-0.1896","-1.0076","-1.0075","-1.0076"],["'APP'","2025","10","'month'","'SSAG'","'abs+neph'","'lin'","'LMS'","0.1349","0","8.4493e-05","0.0013","-0.0012","0.0094","0.1480","-0.1292","-0.9992","-0.9991","-0.9992"],["'APP'","2025","10","'month'","'SSAR'","'abs+neph'","'lin'","'LMS'","3.4718","95","0.0025","0.0040","0.0011","0.2880","0.4539","0.1221","-0.9745","-0.9745","-0.9746"],["'APP'","2025","10","'month'","'BbsFG0_homo'","'neph'","'lin'","'LMS'","0.1615","0","5.4998e-05","7.3607e-04","-6.2608e-04","0.0329","0.4403","-0.3745","-0.9995","-0.9994","-0.9995"]]}}
%---
%[output:2120334f]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAACt9JREFUeF7tXG9IVNsWXwMDGRh46\/nopqHRDQ1u\/6CH2tWsyLqXSszSTC9p\/7CstLLS1FBRS+8ryOSRlZZBmYF+Kijyg\/nopvUt6fnBJvLi0H0U\/fngLYOBef323D33zJlzzt6jR8des2Bomlmzzzq\/vc5aa\/\/22lqcTqeTAjKhCFj0QP\/48SMVFxfT1q1bKTY2dkKN+n+\/mCbob9++pZ07d9KTJ0+otbU1ALrJXuAFOjz83LlztHnzZjpy5AgVFRUFQB9v0Pn43NsDoJuMOBFZPnz44Jw6darHyHV1ddTY2Mg+mzFjBt24cYPmzZtHdrudOjo6aNOmTRQeHm6+NV\/JiJaenh6nMlHyBLpu3To6f\/68R3jp7e2lzMzMQJwfo3NYamtrnQghyrBSWFhI+\/bto5qamgDoYwRY6+eW\/Px8Z21tLfEQw72Zl+8bNmxgiRUS8HRzZsALdCW4eL9w4UJWq2\/ZsiUAujmYk254KSkpYZfYtWsXrVq1isrLywOgmwW6OpFevnyZLly4QHfu3KGRkRFav3497d+\/n3bs2OEGvaCggGJiYkwyYfIM43A4mDFWq9XQKD09VHQyVR1bHKEmr6ysZN48MDBAJ06coOfPn7MLz507l3k7wgtKxqNHj9KjR48mD1ImWALnevfuHXMyLgA+NDSUgoKC3J+J9OCMeInEa0WKZNnQ0EA9PT3st8qYjv8DeLy+JPntN5e1ERHeVsPZ7t69q3s7e\/bsoezsbOaUIr3S0lJ5T1de8d69e2z5jzAzZcoU5uVVVVW0Zs2aSYXz4KDLnMhIfbPu3yeqrCTCv1ygf+UK0YoV+Pw+rVy5UnhfFRUVhJdIurq6aAUGFoiXpwN0hJDm5ma\/gC4CUwQkv1+AbYQTgL96dSUD3izJycmhKxjYV9DV4SUuLo6Sk5NZTIcgifBEIhrcl+8fPLBSTY3VyysvXXJQfLwrwdXVBRkCCd2ff3YQxkpKMk6GLtvg5eaBjhFltic0Yzq4FiyYIODUly1bRhEREYwCePHihS9YfnW68Pbc3FxDZtaL8IKnIwO\/evWKLBYLy96gCfAYonzcu3cv+z4sLIwBapSklIjbbA4aHh6m778P9ijJuruJWlrMmxvEbB6i5EbdLqcmqTVnzhwhN+VFeL18+ZI2btzI4vqiRYvo5MmTbEUK8isvL4\/NIjY2RkZiDZMUt1EUg5HHTAyrktAo1eYQ0Z9ZeRS\/Vv4ESTQkJES4B6G5IsVCqK+vj40HgCEIOaB08fhERbV+Lp\/0t\/CQS3JyXJWDUTLDdxJFwRihEP3cYqgQGRnJ7llUvUAPSbS+vl4Muh7hxS2BZycmJrpBz8g4T7\/\/fkN0JwxMfwGKxItkKhLoZWU1sZCpJVhdZmVlMcBbWlqEetXV1ZSWluY76MqLc2595syZ9PjxYxZekpOnsdAyWSUmZuRzznlHmZnfGpoYHu6gX355TbGxI4zegIcqV9oAHLkLTzcXkR5CC6KEaLfNK7yoLcUuEpa\/g4OD9OOPuZSR4R\/ArVY7ORzGu1XQqagYpB9+cND161ZqaorXBB56q1fbqajor2U\/V8RqW8SfoGR++vQpzZ8\/n61luCDEHDp0SAy6mvDCbHZ3d7MfKvdJMXB6+jEaGvr3hDs5QAoNPcpAf\/36n7pABgd30Dff1Lu\/HxmJ+cypFHg8mRgLOtA1W\/BkoMoTerqa8Jo+fTrl5+fT7du3mU3Lly9n8YyLxTjvjOo+qqvtVFam7cUoAZGU167tZWPbbOHU3BzutYgqLyeKjnbpaAkmDIAbiWz5q6fnE8uoNATeje06zqejZDxz5gxhMiBml3igKrq6XGWjFk8CMAG6WlCL4yVBdbjrdj2eRlTW8mvL6ok8T3NFijgO7gVbeOouL1xYxBFx7xRVL0ryiRsqAtOXG5fRFZW1suUv1xMBju+laQAl99LU5KC8PO3hsVD9XGWxchFRSUYPI5nJvWA8GZ4Gky7D0ciWv3hiZZ48S0pKihPtc9iABt+C98pNDNAA4JQR55HZZ8+eLTOZX60OJ7wQLbADB5k1axZbOKF3iHl6W1ubEywiDyOfPn1iZQ8mICEhwSO8oLKBLrK0zRamyZkg\/iYmemM+UdwLPI2HKH\/MfFvbVIaRUfOtm3u5efMmq8XRw5ienk6YKTCLykTKWzBkaQDctCiump2Y\/QG08poNDbcpM3OZuxjh3q3UsQwMDDjxBUB\/+PAhpaamMq4BDCMEzKKy7yU1tV6KBujsdNCvv1q\/AO7F3GkaGHhGb968YTQ4F1ApyoYuS19fn3PBggVu0HnjESgAMI3YsOWbGPB0eKYMDYBltt0u5j\/MvWX\/jsbLX6UVnErBngQvRizNzc1O8AU8vChXojzBIhHg8wcP7JSQMLkbR3n1IEMXy3DvYyl\/OfgI1RDu7Zbs7GwnuHJlIj179iydPn2aDh48yHgXbEyjydSfoMtyL1jdQoqLjTkijAddJN3GRm1dztFggXb1qljv0iVX34selcIbdS1Llixxvn\/\/nm1O8JlQljsoJcHAQVwl48R7+mi4l+HhTdI8zV8czcw\/nTOSUQZqjkakh6oOjCzCiBJDr5iu5tP5I6F3KGA8qo3x4l5keRoQfODM1c1G4Mexl8BFpCfNvahbpXEBXhri\/e7du+n48ePuC8vSAPiBzF7laLkXPr5R3ws32ohaQBOR0a4QFjWo5mT1ZFK51x4pJ7zQn15WVsb4F\/AwnPDCoFjeb9fZz+WJBw4iw9Hwxh+lsQBJBKaoP0ZmjeDXZiN1LyPiERLp4cOHKTg4mA4cOODVUiDDCspMDpKULyJabPGxZIisSdVspNX3otds1N7eTiUlF1niQb2PSghPBbb38Fgi+QQHnyabbYYbW6v1vxQS8i8qKAhhKzcI1wfvjc9yc9d6zYUMiTVZmo0MuRf14V3Ecz3QA4SX+JlE8y3O4BpyL1qg6\/HpSsKLNxuJzXBp6BFe+A6hwCjp+kJiySx4PG02t9mora2NkpKSjLkXNbWLZWtKSgoNDQ0x27Dxeu3aNVq6dOmoTmKIYrBMNSQ7saPTM7fZ6OLFi2LuRU3toqsL7dG8P115TN3Xg15yyczfHV7mNhupW6U1uRfeDcC5F9TlfI9UTUv6Aro\/PdjVey72e+hlZ7fQdp36l3d34YQKNudl9LSu6sW9aFG7ysGVS1gOOg4MiP4yxk8\/BUnduBga3zTQtVVa6iBc30jAgqK1mh8OwAJJ3WyEkxWowLigrjfSw5EZvITcixp05ZlS9aPBQd+2bRstXrxY954cjjBKS\/uHb2iZoI2ytbDwKevaam8PpuvXV2uOCr34eBsVFLz3+l622ai\/v5+ioqI8mo2io6PZdqeQe1GHFyXZDouUj4bsQS\/U2mY3JQUF9dK0aR3SJBZs90ezkdRBL3UihbFaHV48nMge9IqLM6\/9DmVgVZWdvvvOPq7NRmN92KQJL14y6lG7alpS1jBZNtKVzMRcjhZdIMPRyNo7kXq6f05qrEbIVC+eJ9186\/Aaq33+\/P24gY6bGi3h9aV6sOxEjhl01Pecb8fCitPAvNJBMsMOTH\/\/39028V2ZwsK\/uXerlBz+qVOn3Ju4sjdilp7ofnAdbh+v7m7dusUurxeivZqNxvInAp89e+bRF4NKB2eWjh07xrh4dRNqa+tD+uOP\/3i0YWdkZHhwFTBe3bRqFqCicXy9n87OTtYrpNzMx\/2Imo3+B5f\/hMLpEiugAAAAAElFTkSuQmCC","height":56,"width":93}}
%---
%[output:41032e27]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAF0AAAA4CAYAAACfSPQqAAAAAXNSR0IArs4c6QAADfpJREFUeF7tXA1sVFUWPoODlFB+iiERGMM0gApJ0QRhqRQdglV+FoGWhJ8CFgoiAhJFQCjSlhahNLgZamiwYKugtEADSKO7NkAplbaGuMhqVIp0sBOaWAvoApbsbGb57uudvnnz3rt3pjN23eUmTTt9592f75577jnfPXcsXq\/XS6GUykqirCwi\/ObFbicqKiJyOKi0tJTOnj1L2ydPpm7btvnJNXXrRv\/es4ds8+aF0vIf\/h1LSKAD7MxM48ED+NRUZVJk5NQ1uVzKJ0zg\/2hhoEMrXS4XrVu3jg2ztraW5s6dy\/7eunUrzZo1q3340Ozx48VwAGwzwHkNp06xlcFWjMnKCWjwDzw5lm3btnl3795NS5cuZaBfu3aNVq9eTRs2bGDjfOutt2jHjh3Ut29fZdwAXG1SxPCbS2BFQKtlV0Swk9PR\/kXgfUtNTY339OnTrGqADi3Pzc2lvXv3Uvfu3emNN96gOXPm0JgxY4igXbGxEeiGRJVYEein7ORIVNlZIsy8AGQ16AcOHKBt2PyIGOhPPvmkYmI6E3SsBm5STNDyVFSQJyEhonh6PB66efMmRUdHk9Vq9bWFv9WfjTqhC3p+fj7V1NSwd0aMGME0HaBjFYyJj4\/ogDpa+c3kZGrOy+toNabvX7lioebmZoqL60ndunXzyfbp04diYmKEbQeA\/tlnn9GaNWuYeUGFixcvpuzsbHr22WcZ6LDpY1pbhRV3psD1a9ci0vznn3elnBwr1dVF+eq32TxUUNBKY8f+i5njqKj2Z9KaDtBff\/11eu+993RBdyYl0YGmJtNBeWw2+mdyMsU4nUI5q9sddoAaq6oIfQhncTr7kNNprMV5ec2UlmYNTdOhzWrzEh8fT88\/\/7zPvIwbN45sHo9wPBh0JAAVNkwUdsDRptvdbruN+oBtp6GhQdjFgOAIoBttpHiGSbhXjBGQCfB1QTdyGQE6ND0SxWj1uFXegWiFqWXD1UePJzgz5fWGoOlmwZHb7WabbF1dXbjG5KvnT62ttOr6db9NGiA6Y2KoLDqaySXfvEl5zc26bUMWcpAPZwHojY1V0lXCuogYDF3uRU0DfPTRR0pg1FYAPH4iWaDRRlpru3SJbHv3BhJtGRlU++ijEelWfHz7+EUNyNCHQtADuBdRq7\/ncwRLIrUKQ39kmQ9QSAicRSUAdCH3oqlRlneSlRN1uDOey3B8KlZb2MWgNlJ1bbK8k6wcr1t2cmTlhAhIChQXEy1cqC8MwMHbZWTIVRaUy4gqwTts2WI15Z0KCz308MNnadq0v9ONG6sMe5KcXE779z\/DnhcW1tPatbeotbXdfmIwqCshoT0uqK62sva1ZydaObnhB0odPnyYNm3axB7ExcVRQUEBY1i\/+OILmj27gK5fX0V37jzuezE6+gZFRe2g6OgyWrJkCa1fv549gwcI9hZlwIABVFRUREOHDmWfgwIdG+hDDz0U6nj+L94rKSlhwaQfO6sZeVDmRfHTT5DHM\/B3B1CSZGTL\/OmnI9M9MMswM9qibhP8S2Jiou9Mgmu3+h3L9OnTvV999RVNnTqV0bm\/\/fYbe+G1116jt99+m27cuMHILyyx6mo3jRsXXLAQmeEb15qS4qE9e8Q0RbD9ys2NEppUgA8aoKWlxXfyhnb4ARFv01JSUuLVLgeQXi+99BKT2bhxIy1atIj9\/UcAHf28fFkcFQYDem1tFM2d21\/4yqefttLEif4sI5TY70wCNh0nRwh++DnpK6+8Qjt37qSZM2cythGnSTw4UrgX+UBB2MsICXz7bWtY3ffERCthAxcVaDrO5LVFfUjENtKLFy96YXd4ygRMDOwS\/PW0tLQA0BEoqD0MUUc643lNTW3YmgUNEIxJRUQK5cQRKD9z1uIYFOjwXlauLKOjR43dQIzWanUzF8rMXeRy+B0sqWSGaFRULfXvr2QyhKPocS92UtJEXBSYJpKfX04rVvzZz2UMsOla88LTMPQ0HQ0B+Pffh63X31B5oPDcc7V04oRNKDd4sJvmzzffnHmdoowOyGVnu2nIkPByQ9ykOqiSMiiL8JsXAL+QiqiSHOxfUtyL3kaKl41A540ZRZqIymDbgpGTjfZk5ULWcIMwFybVUZlFGWScYAXgXY5UOe5l+PDhXuywgwcPpo8\/\/pjZc5SrV6\/SlClTaOXKlT7vxWgwsryTmVw4JzFo0AVchau4kuwLxQlWrqJTZE9VNN6s6LqMXMvhv2upXVGFHX0OP\/fo0Sp69dUXhFXJTrYvdUOPkZRJ\/YM9lUmwMnJfNCMJymXEu+Be8MOLGVeR2mZnNm\/ezFxQrKg333yTPvnkE\/Y6dnUEYijgNhYsWEAWi4UmTJhA77zzjhB0kYC1upqsW7YEcO+ewkKWG8OeJyaKqgnueZtRN+VegnEZsYl++OGH9NRTT9GDDz5o2pnz589TeXk58\/WRloCJOnbsGO3atYu9h1S9xx9vJ45wtnjmzBnq1asXi+pmzJhhWj8\/SLEZnPr3cTpNsxGQG9Pz8GGKCvMpWHl+Pk1ISzPnXoIBnZ8offDBB\/TEE08YgsK1f\/LkySxnBvvEN998w2iG7du3U2NjIx08eJBp\/aVLl+jOnTt0\/\/33s\/B51KhRbCXMM0ijrry7zDMzM\/2ODAF8YWEhOXCKgFJZSVGTJgWnoWGSLi0pEXMvwbiMHHSn00mHDh2itWvX0sCBAxl4iFphQmAmfvzxRzYEyHPQkWGAgAEnUdBqhMZYBf369WOf6+vr2bswQfgbCayxmrxJLFkAblQAPCaLAS5jg8MEtK8ah4Pq3303NO4FYCGVDgUcMc\/g5aCDCOvduzfjjn\/66acAU4H3oK1ffvmlz7zAtPzwww+0fPlyam1tZXWibrWJwXu\/\/PILnTx5MsC8qM9tzbA6lJdHM9esCTecwvqQ54O9ImriRD9ZXe6Fs4zaVGmAAw1Gah2yvcAyas0L7HJTU5NPm9WtwcRoNf27775jE3X79m1atmwZ20RHjx4tHBAEQJdWV1cLZVMSEmi\/hJywIrWAw0GelBSyLlmi\/1pb9OZJT9dNIA3gXrQ3McxSpfEMjGRycjJ17dqVaTIyeuHTP\/OMcgLEC55hQ4Q8ZM+dO0cXL15kFPJ9991HR44cYfmRMC8yZaHRWZnOyzqck0wTxjKcMDci1LGXvKC4uPDYhNyLHuhGGV785KhLly4sTfjXX39lDfXo0YP9vnXrlq\/j+B82UGyOsNlIIcZq+fnnn9kk9OzZ0\/esY4j8d72NkyNkOKtdxgDuJRjQMTzZvJdTp07R119\/zQh8eCYo+B82uwceeIAxcEbunh6MwaTzNe7bR7b5881nw24nd3Y2uYcMoXDl0mA8MmOyxMbGstt1jz32GDshggnAwcXly5dZp5GfDk9DnXDUGbo1fvx4grsoKnAbMbnsXC2U43ttmKt3PifqhJp80pG1tLS0eLHssRxgm1988UUWGaqNP3xqcDOdWQA4gDcrdrudnbqr\/XXdy2NaVs6sUtlMI16HRMaR72CaZ+smJSXdTXzPYUELQnIsa5Bena3pGFNxcTEZbagAHJtYhlHySagafOWK3C1BdFAy48gHOjQbHR80aJBhqjTqjST3ouVq9BSQR6QIpHhBgJaens5Aly3gXMC9hLP8Y+RIijt3jlUpzHtR3yM1y0+X5V7Kysp0x8KJMr3LUGJOUakSd4rMilHb2ndwYQH8DIpoqmQuOLDbJxcuEA71hXkvuEcKDecXdEV+Oi71iriXSZMmSW16PjMIz0ZC5Vx2O\/217YKxkTg0DBeRZQvYGmHbqankGTvWODiCBSgsZKc3mHRh3gvnXngnOZcOlhDlkUceYZy6OiJFhGpm4zHoYcOGSY0bp4wIZkTUP2BEKqHYf5FqlgnBvdvndpu2zcL7+nqCImVU8kM5\/zbQt\/kJCcxrksp74S4jqkG0CP\/55ZdfpgsXLrCa4TJqaQBwL3AxzQouD8gsdSxtmQgSiVUG+Zt+3eB+skwOPSLrR+vqKMck336jzUZpVVWMzh7iduuuCvhUl2w2hpn2SqMu9xIKDTBt2jQaOXKkKejg3WW4ElQiA7oM4KhLMZXrKDd3mVDlZ80qoNLSu8wluWiQgTTa5XI4\/wfHqZa969v4Tk7T0wspJ2dxQE1S3Mu9i17C+TIUAOURVu7l3u068WSEzL0Y3a6DnQQlK3PRC9e4gynYUPVy6u9+Y0xbag8zHn6fAutX1+L\/prFsuxw3czAZ2EMUH4jX6V8f9iL\/BF5FrqEhVZjSF\/T1F1nCC7y5DFeiBgOum9qLgafSHvjjCWAwcwcxcJ48im\/rkJFtl4O9xsT7b9q8Trn6pJKN9L7ZyOx2naz2ynIliCL5EZwaMrWLiM0R263LhW3NHEiHQwG9slIOJCJ\/OUy8P8MjDzraDumilyyoMnJmzKACpJLPrZbTahuXKypS5NpB5+C35xMq1IcCuoxsamoDZWaq69QxQuw+kZwc2uZn42b4hPYdXjKIt8lA47OysvxMDYAEMaXmStRyXNu0cv5sLWwspkjx9LWXrQKZXYCrxJ5qWVkGWFZOBpoOgw7ehl9u4py8OnpFJ\/hdVJyR4uLU8ePHWd\/UJypqk4ZbICtWrNDtv38GnAIkQNRjazuSqofM4xEjjlBFRYovGk9KcrKLXupU8V69rt09R97KspTV4xESXjKzoyeDVAn1d3xxTh6pGTgI0X4PWEVFhe8L2jjdMHv2bD+uAu0EfG+YQQel0+pglCTu+SqpH+\/S7t3rGdBm41m16i\/0\/fd\/0x2PiPD6D\/TY1dGUjIcnAAAAAElFTkSuQmCC","height":56,"width":93}}
%---
