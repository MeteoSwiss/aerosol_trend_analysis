LLN_st.name='LLN'; 
LLN_st.lat=23.4686;
LLN_st.lon=120.8736;
LLN_st.alt=2862;
LLN_st.env='Mt';
LLN_st.footp='F';
%%

LLN_rd=read_ebas_merged_STN('C:\github_trend\raw_data\merged_Lulin.nc',LLN_st.name); %[output:9fa424f5]
% files have only similar data, no more
% LLN_rd2=read_ebas_merged_STN('C:\github_trend\raw_data\LLN_DJ7T-6ZFV.nc',LLN_st.name);
% LLN_rd3=read_ebas_merged_STN('C:\github_trend\raw_data\LLN_V3VX-EYMA.nc',LLN_st.name);
% LLN_rd4=read_ebas_merged_STN('C:\github_trend\raw_data\LLN_Z5CU-4RW4.nc',LLN_st.name);

names_rd=fieldnames(LLN_rd);
%%
% new data from Betsy since they were not yet submitted to ebas 
% keep the old way and multiply AE31 and AE33 with the right constant
% ATTENTION, redo the importfile_LLN since it is now for Ae31_AE33
%%% LLN_rd2=read_betsy_2026('C:\github_trend\raw_data\lln_tsi_psap_clap_AE31_AE33',LLN_st.name); %[output:838b327f]
names_rd2=fieldnames(LLN_rd2);
%%
% adapt AE31 for C=3.5 and AE33 for C= 2.45 (new filter M8060 confirmed by Jene)
P10=timerange('2007-01-01','2023-11-01');
P11=timerange(2023-11-01','2026-01-01');
for i=1:length(names_ae)
LLN_rd2.(names_ae{i})(P10)=LLN_rd2.(names_ae{i})(P10)./3.5;
LLN_rd2.(names_ae{i})(P11)=LLN_rd2.(names_ae{i})(P11)./2.45;
end
%%
datevec(LLN_rd2.Time(1)) % 2007 %[output:7200f62b]
datevec(LLN_rd2.Time(end)) % 2025 %[output:5edb9f45]
prctile(LLN_rd2.U1_S11,[5 50 95]) % mostly < 40% since 2012, but compute RH trend anyhow. %[output:0a79c560]
prctile(LLN_rd2.U0_S11,[5 50 95]) %< 40% %[output:01e0b367]
plotFigControl(LLN_rd2,LLN_st.name); %[output:24552bb3] %[output:5f489c33] %[output:5d98ae53] %[output:3874273a] %[output:2b90f89b] %[output:43492d8f] %[output:2a1c4d19] %[output:065be3b7]

% 2 size cut, ok
% sc: ok, a double seasonality with max in spring and second max, lower, in autumn, no negatives
% abs: 2 size cut, a double seasonality with max in spring and autumn,, no negatives
%SSA: seasonal cycle with min in winter. low values : max=0.95
%%
% statistical BP detection for scat and backscat
names_sc={'BsB0_S', 'BsGtsi0_S', 'BsRtsi0_S'};
break_LLN_scat_BP=change_point_analysis_Def(LLN_rd,names_sc,0.05,'LLN','scattering'); %[output:085f9d5e] %[output:12369ae2] %[output:17485e5f] %[output:15045ae0]
T_LLN_scat_BP=make_table_breakpoints_def(break_LLN_scat_BP);

names_bsc={'BbsB0_S', 'BbsGtsi0_S', 'BbsRtsi0_S'};
break_LLN_bscat_BP=change_point_analysis_Def(LLN_rd,names_bsc,0.05,'LLN','Backscattering'); %[output:0207f450] %[output:69a175b3] %[output:7eb82644]
T_LLN_bscat_BP=make_table_breakpoints_def(break_LLN_bscat_BP);
%%
% homogeneisation PSAP CLAP
% common period 1.8.2011 to 31.12.2023
P_PSAP_CLAP=timerange('2011-08-01','2023-12-31');
ratio_STN=nanstd(LLN_rd2.BaG0_A12clap(P_PSAP_CLAP)./LLN_rd2.BaG0_A11psap(P_PSAP_CLAP)) % ^NaN %[output:5ec3f60e]
ratio_mean=nanmean(LLN_rd2.BaG0_A12clap(P_PSAP_CLAP))/nanmean(LLN_rd2.BaG0_A11psap(P_PSAP_CLAP)) % 0.8380 %[output:2a69195c]
ratio_median=nanmedian(LLN_rd2.BaG0_A12clap(P_PSAP_CLAP))/nanmedian(LLN_rd2.BaG0_A11psap(P_PSAP_CLAP)) % 0.9620 %[output:6851978f]
figure; %[output:85ce0528]
plot(LLN_rd2.BaG0_A12clap(P_PSAP_CLAP),LLN_rd2.BaG0_A11psap(P_PSAP_CLAP),'.'); %basic fit: y=1.078x-0.05362 %[output:85ce0528]
ylabel('BaG0_A11_psap'); %[output:85ce0528]
xlabel('Bac1_A12_clap'); %[output:85ce0528]
grid on; %[output:85ce0528]
legend('period 2012-01-01,2013-08-07','fit') %[output:3036bb08] %[output:85ce0528]
% fit between both data: slope 0.9955x+0.1021
fit_test=robustfit(LLN_rd2.BaG0_A12clap(P_PSAP_CLAP),LLN_rd2.BaG0_A11psap(P_PSAP_CLAP)) % 0.0295, slope=1.0639 %[output:9ee9f073] %[output:5a5b97b3]
% 1/slope fit and robust fit =0.93: use it as constant
cte_PSAP_CLAP=0.93;
%%
   % time serie PSAP+CLAP
P1=timerange('2008-10-11','2011-08-15');
P2=timerange('2011-08-15','2026-01-01');

% 
LLN_rd2.BaG0_psap_clap=NaN(size(LLN_rd2.BaG0_A12clap));
LLN_rd2.BaG0_psap_clap(P1)=LLN_rd2.BaG0_A11psap(P1).*cte_PSAP_CLAP; 
LLN_rd2.BaG0_psap_clap(P2)=LLN_rd2.BaG0_A12clap(P2);

LLN_rd2.BaG1_psap_clap=NaN(size(LLN_rd2.BaG1_A12clap));
LLN_rd2.BaG1_psap_clap(P1)=LLN_rd2.BaG1_A11psap(P1).*cte_PSAP_CLAP; 
LLN_rd2.BaG1_psap_clap(P2)=LLN_rd2.BaG1_A12clap(P2);

% do the same with B and R for abs exp
LLN_rd2.BaB0_psap_clap=NaN(size(LLN_rd2.BaB0_A12clap));
LLN_rd2.BaB0_psap_clap(P1)=LLN_rd2.BaB0_A11psap(P1).*cte_PSAP_CLAP; 
LLN_rd2.BaB0_psap_clap(P2)=LLN_rd2.BaB0_A12clap(P2);

LLN_rd2.BaR0_psap_clap=NaN(size(LLN_rd2.BaR0_A12clap));
LLN_rd2.BaR0_psap_clap(P1)=LLN_rd2.BaR0_A11psap(P1).*cte_PSAP_CLAP; 
LLN_rd2.BaR0_psap_clap(P2)=LLN_rd2.BaR0_A12clap(P2);
%%
% Bp detection of homogenise abs + AE
names_abs_tr={'BaB0_psap_clap', 'BaG0_psap_clap','BaR0_psap_clap', 'Bac3_A81ae31'};
break_LLN_abs=change_point_analysis_Def(LLN_rd2,names_abs_tr,0.05,'LLN','abs three homogeneisation'); %[output:3bd33ba9] %[output:0a0bb36c] %[output:413735f7] %[output:42f88105]
T_LLN_abs=make_table_breakpoints_def(break_LLN_abs);


%%
lambdaSC=[450;550;700];
lambdaAE=[467;530;660];
lambdaAE7=[370 470 520 590 660 880 950];
LLN_cal=compute_exp_SSA(LLN_rd2,lambdaSC, lambdaAE);
names_clap={'BaB0_psap_clap', 'BaG0_psap_clap','BaR0_psap_clap'};
LLN_expAclap=compute_exp_D(LLN_rd2,names_clap, lambdaAE);
names_ae={ 'Bac1_A81ae31','Bac2_A81ae31','Bac3_A81ae31','Bac4_A81ae31','Bac5_A81ae31','Bac6_A81ae31','Bac7_A81ae31', };
LLN_expAAE=compute_exp_D(LLN_rd2,names_ae, lambdaAE7); %[output:07104b66] %[output:5993295b] %[output:8ff1526f] %[output:796f32f2] %[output:610d3fbe] %[output:82feb1d5] %[output:04097e2b] %[output:03ece419] %[output:2c35047f] %[output:776256bf] %[output:11816353] %[output:977393f4] %[output:7a7ed79e] %[output:94698f8e] %[output:8c6b983d] %[output:01aa338a] %[output:68863217] %[output:45d1a085] %[output:32252d6f] %[output:00ce1993] %[output:4963f04c] %[output:64249436] %[output:05163949] %[output:931afc20] %[output:194e8c1f] %[output:18144483] %[output:91a029a5] %[output:4788142a] %[output:2d9bfcc4] %[output:34cfd756] %[output:8941d291] %[output:137fdab4] %[output:75eac030] %[output:945a7157] %[output:6d90b171] %[output:6e156428] %[output:45de0648] %[output:8acdac76] %[output:78df5589] %[output:48a509d2] %[output:79ab9148] %[output:0575dff8] %[output:9e6e61dd] %[output:44679edb] %[output:32971274] %[output:5b462c57] %[output:4439e1a3] %[output:5c2dab60] %[output:80b89144] %[output:9a6b778b] %[output:92f26cd4] %[output:8bd617f0] %[output:6c6b7f0f] %[output:6ab5fc69] %[output:0415e9d0] %[output:585f473b] %[output:2015c8d7] %[output:5aa2591a] %[output:14b9eb2f] %[output:0e8cf68b] %[output:93b17453] %[output:07e48f11] %[output:19232560] %[output:60f6926a] %[output:37418d3a] %[output:9251ad8a] %[output:227c8651] %[output:2a325c28] %[output:6d87e43b] %[output:9ff13350] %[output:63268dba] %[output:84a2d3f0] %[output:59723cec] %[output:95b9da20] %[output:80b5adde] %[output:233439d3] %[output:22bb79b7] %[output:16261198] %[output:35eaf959] %[output:3e7a3b38] %[output:2e6b76c5] %[output:9285ce9c] %[output:18f33a24] %[output:4c26c83e]
LLN_expA=outerjoin(LLN_expAclap,LLN_expAAE);
LLN_cal=outerjoin(LLN_cal,LLN_expA);
LLN_cal.SSAG0_homo=LLN_rd2.BsG0_S11./(LLN_rd2.BsG0_S11+LLN_rd2.BaG0_psap_clap);
LLN_cal.SSA30_ae=LLN_rd2.BsG0_S11./(LLN_rd2.BsG0_S11+LLN_rd2.Bac3_A81ae31);
LLN_cal.SSAG1_homo=LLN_rd2.BsG1_S11./(LLN_rd2.BsG1_S11+LLN_rd2.BaG1_psap_clap);
plotFigControl_cal(LLN_cal, LLN_st.name); %[output:296de27f] %[output:6d0fa995] %[output:47d3fbe7]
%%
% statistical BP for computed value
%exp
names_exp={'expS_bg0','expS_bg1', 'expA_bg','expA_fit'};
break_LLN_exp=change_point_analysis_Def(LLN_cal,names_exp,0.05,'LLN','Exp S and A'); %[output:7a5e52d5] %[output:3d3e90ea] %[output:6a568048] %[output:2b042988]
T_LLN_exp=make_table_breakpoints_def(break_LLN_exp); %[output:12dbb617] %[output:6baa8017] %[output:58c0f89d] %[output:4d9e7c12] %[output:83f879fb] %[output:0c55df08] %[output:9ce16302] %[output:1c656e10] %[output:1594056c] %[output:6a0f39a7] %[output:57b0d064] %[output:333aac54]

% BbsF
names_BbsF={'BbsFG0', 'BbsFG1' };
break_LLN_BbsF=change_point_analysis_Def(LLN_cal,names_BbsF,0.05,'LLN','BbsF'); %[output:671e0000] %[output:139bcf99] %[output:7e199d6d] %[output:880f0b10]
T_LLN_BbsF=make_table_breakpoints_def(break_LLN_BbsF);

% SSA
names_SSA={'SSAG0_homo','SSAG1_homo', 'SSA30_ae'};
break_LLN_SSA=change_point_analysis_Def(LLN_cal,names_SSA,0.05,'LLN','SSA'); %[output:7ca75747] %[output:9365cd62] %[output:08489f7f] %[output:7ac8d054]
T_LLN_SSA=make_table_breakpoints_def(break_LLN_SSA); %[output:88f33f76] %[output:02ee8dfe] %[output:6b77bd44] %[output:243631ac] %[output:4d8b2271] %[output:60759f32] %[output:854eaa9f]
%%
%Questions:

%%
LLN_tr=outerjoin(LLN_rd2,LLN_cal);
LLN_tr.y=year(LLN_tr.Time);
% begin at the beginning of a year: 2009
%end: 2025
P=timerange('2009-01-01','2026-01-01');
LLN_tr=LLN_tr(P,:);

names=fieldnames(LLN_tr);
c=startsWith(names,["SSA0";"SSA1";"expS_bg2";"expS_br2";"expS_gr2"])==1 | endsWith(names,["A11psap"; "A12clap";"1_S11_dry"]);
N=names(c);
for i=1:length(N)
    LLN_tr.(N{i})=[];
end

%AE31 should stop end of 2024 since 2025 is not complete THIS DOES NOT
%CHANGE MUCH TREND. SEE IF FULL 2024 HAS TO BE REMOVED. Yes for expA since
%AE33 leads obviously to high expA.
P1=timerange('2025-01-01','2026-01-01');
c=contains(names,["AE";"81";"ae"])==1  %[output:59d3a130]
N=names(c);
for i=1:length(N)
    LLN_tr.(N{i})(P1)=NaN;
end


names=fieldnames(LLN_tr);
%problem with RH that decreases at the end of 2013 (rupture)
%compute trend in RH, but not in dry (no sense)
%make one time serie with PSAP and CLAP: change in 2012

%LLN_cal2=compute_exp_SSA(LLN_tr,lambdaSC,lambdaAE);

% plotFigControl_cal(LLN_cal2, LLN_st.name);
% LLN_tr=outerjoin(LLN_tr, LLN_cal2);
% LLN_tr.expA_bg1=real(-log(LLN_tr.BaB1_A11./LLN_tr.BaG1_A11)/log(lambdaAE(1)/lambdaAE(2)));
% use  green and expA bg
%compute abs and SSA trends only on G, the longest time series
LLN_tr.BaB0_psap_clap=[];
LLN_tr.BaR0_psap_clap=[];

LLN_tr.Bac1_A81ae31=[];
LLN_tr.Bac2_A81ae31=[];
LLN_tr.Bac4_A81ae31=[];
LLN_tr.Bac5_A81ae31=[];
LLN_tr.Bac6_A81ae31=[];
LLN_tr.Bac7_A81ae31=[];

LLN_tr.BsB0_S11=[];
LLN_tr.BsR0_S11=[];
LLN_tr.BsB1_S11=[];
LLN_tr.BsR1_S11=[];
LLN_tr.BbsB0_S11=[];
LLN_tr.BbsR0_S11=[];
LLN_tr.BbsB1_S11=[];
LLN_tr.BbsR1_S11=[];

LLN_tr.BsB0_S11_dry=[];
LLN_tr.BsR0_S11_dry=[];
LLN_tr.BbsB0_S11_dry=[];
LLN_tr.BbsR0_S11_dry=[];



LLN_tr.expA_br=[];
LLN_tr.expA_gr=[];
LLN_tr.expA_bgAE=[];
LLN_tr.expA_brAE=[];
LLN_tr.expA_grAE=[];

LLN_tr.expS_br0=[];
LLN_tr.expS_br1=[];
LLN_tr.expS_gr0=[];
LLN_tr.expS_gr1=[];

%%
[LLN_result_MK,LLN_result_LMSlog,LLN_result_LMSlin]=all_trend_STN(LLN_tr,LLN_st); %[output:31ea7669] %[output:6f852e64] %[output:1c9e8752] %[output:293586e5] %[output:123f3ad8] %[output:42132fa3] %[output:6353cdaa] %[output:96b54ad6] %[output:8801b32f] %[output:4cd076ee] %[output:65fa3869] %[output:013471c1] %[output:08148c64] %[output:8e48a980] %[output:3583ea30] %[output:553c6d0f] %[output:4f2e2e2f] %[output:49423185] %[output:0fd34b35] %[output:5de37e85] %[output:14e8daa9] %[output:7b7116bb] %[output:2e7fb7e3] %[output:656fb5bb] %[output:99849bd9] %[output:6ead520c] %[output:4e74eefe] %[output:0fd7a76d] %[output:37fad56a] %[output:55c4a5f0] %[output:94ac4a83] %[output:2134c1d8] %[output:5462d411] %[output:02656171] %[output:6e6cfc9d] %[output:0ecbb07c] %[output:65e26383] %[output:4e2a3163] %[output:0afaa83c] %[output:2294c197] %[output:0b5938fb] %[output:7cf75b5f] %[output:04b104f6] %[output:34002a4f] %[output:773da4a2] %[output:7b31cf5a] %[output:670bb4d4] %[output:399774f9] %[output:7579b753] %[output:25e94361] %[output:5b07f81b] %[output:9711eebc] %[output:789c8619] %[output:319c2495] %[output:5b755ac4] %[output:3396abf6] %[output:59e5c958] %[output:160e7865] %[output:1f821d55] %[output:43f484e8] %[output:0384aac9] %[output:763fd066] %[output:2d2f8cd3] %[output:38746add] %[output:3aa46ccc] %[output:9a403478] %[output:034ad521] %[output:7832eb2a] %[output:740d9a63] %[output:635fccc1] %[output:8493fed0] %[output:48cb69f5] %[output:2cac285e] %[output:61e85ed8] %[output:37f0b234] %[output:5df85d49] %[output:5f05ee8d] %[output:4e5e17fb] %[output:35a75de3] %[output:2cc0df15] %[output:7e0e5485] %[output:45ac046a] %[output:546934d8] %[output:8041e268] %[output:9f745a68] %[output:16a5f179] %[output:87fb0582] %[output:7df2602c] %[output:58d22566] %[output:7cd57c45] %[output:696de59c] %[output:35847a0a] %[output:8cec2449] %[output:828859a2] %[output:47ef7e5e] %[output:5e26f736] %[output:9fcbbc89] %[output:24b0785b] %[output:2c0c6292] %[output:371cb65b] %[output:049c0398] %[output:2be58b7d] %[output:4f17db34] %[output:03423f8a] %[output:17947987] %[output:4368b427] %[output:1d061903] %[output:8a0e0d1c] %[output:181ad1ff] %[output:15f1e605] %[output:61cf842e] %[output:8130ed1b] %[output:80f2e0d5] %[output:7aa61f2a] %[output:6f3579db] %[output:0f197060] %[output:3738de88] %[output:5ed8dd4f] %[output:0982242e] %[output:3c40ded9] %[output:75d565b9] %[output:897947ca] %[output:44bacd6c] %[output:8dfd0d99] %[output:416b9254] %[output:402e0fa7] %[output:126097c0] %[output:344eb300] %[output:9a49d798] %[output:86c24713] %[output:4fca1f5a] %[output:458428a9] %[output:7bcf3233] %[output:1d374ea9] %[output:0c55ea45] %[output:4a4c51da] %[output:4c03b7c6] %[output:82fb2d74] %[output:0a175418] %[output:19d9720a] %[output:9f5ebeb3] %[output:28b367ef] %[output:600c9891] %[output:5a0b3f92] %[output:9bdfd3f5] %[output:9a4087b7] %[output:8c2ec768] %[output:8fc2e19e] %[output:1d6f9b4b] %[output:60001037] %[output:520f45a3] %[output:1c87610f] %[output:1a31aab9] %[output:191751aa] %[output:6577cca2] %[output:606aea90] %[output:0b03daeb] %[output:97c009fd] %[output:5e7d2bc7] %[output:4f8483cd] %[output:59d715d8] %[output:3b8bd2b0] %[output:1c9303e7] %[output:3125e871] %[output:62fd9624] %[output:8ee1cd2c] %[output:8ffb7531] %[output:318d36f2] %[output:52da2a22] %[output:3a5ceb78] %[output:04b307a2] %[output:48ed9cf6] %[output:6d83a0a0] %[output:30f3aac4] %[output:0d62ea9e] %[output:71597e64] %[output:68921c4e] %[output:8f2e2fba] %[output:415ca299] %[output:2cc8c10d] %[output:5876c111] %[output:2fb6da6a] %[output:29b37f41] %[output:13f3f9d9] %[output:799394b5] %[output:219049db] %[output:596658d4] %[output:22d84503] %[output:05c96b02] %[output:2e914c57] %[output:32e40d39] %[output:20845ffb] %[output:8ecfe50b] %[output:9a028c52] %[output:3ea31299] %[output:2392a08e] %[output:34088c29] %[output:413775fa] %[output:5f994791] %[output:668f2837] %[output:1322b157] %[output:8647a23c] %[output:960b51ad] %[output:2a6ccdeb] %[output:309d3545] %[output:6edcbdfe] %[output:14ab9601] %[output:67565d66] %[output:1629796f] %[output:2b8595ea] %[output:4ca27c84] %[output:2bd8b0d8] %[output:214b4aad] %[output:9a64ade7] %[output:159b4bfd] %[output:42d2544a] %[output:8f2bbc20] %[output:5ff0a4f9] %[output:3f50eee2] %[output:6b50a1dd] %[output:26e76127] %[output:2ef9f7ee] %[output:6b96a88f] %[output:77e967f1] %[output:09afd9ee] %[output:93487815] %[output:95ed6b35] %[output:095f36e5] %[output:0f87b71d] %[output:3d60336b] %[output:7d2c88f0] %[output:97e5bf60] %[output:69feb6dd] %[output:6585dd1f] %[output:5a3e6287] %[output:45daeb50] %[output:669ba1db] %[output:6c2b64f8] %[output:7a008c00] %[output:829089e7] %[output:6aa36e28] %[output:3aaf4034] %[output:95246d47] %[output:29ae58c2] %[output:7083da54] %[output:2d1994a8] %[output:05d8ad38] %[output:1a88a226] %[output:0ba54330] %[output:9fe470db] %[output:2bb41207] %[output:0cea6083] %[output:16dda99f] %[output:12d92b20] %[output:45313845] %[output:18a2c6e4] %[output:12b7761b] %[output:5da8725d] %[output:2e2ec90f] %[output:3a4fe584] %[output:3e25e0b3] %[output:10ef21ce] %[output:7548a969] %[output:9d3d04e9] %[output:947f7fff] %[output:17290c99] %[output:537a0183] %[output:7e2856f6] %[output:87bd088d] %[output:4ce9bc13] %[output:84525409] %[output:39f91975] %[output:8f91a60c] %[output:70a5c384] %[output:6e1d5c24] %[output:2ddc2551] %[output:2e8f6f4e] %[output:5a469311] %[output:620319d9] %[output:157801df] %[output:5d3706a3] %[output:5e51df3a] %[output:924a27ac] %[output:6f768d52] %[output:443f3dcc] %[output:588a78d3] %[output:99d35ce7] %[output:8ca35004] %[output:0fcadb06] %[output:432f908a] %[output:2a74e4e4] %[output:2cdb7819] %[output:13e84aea] %[output:02ee6ae4] %[output:34f3f410] %[output:54b2f813] %[output:4e709e51] %[output:2debc682] %[output:3de94c27] %[output:22413a52] %[output:96023bd8] %[output:10142ea3] %[output:25860191] %[output:53db19dc] %[output:8fd15a49] %[output:01f3075f] %[output:4b1488d9] %[output:4fdb4704] %[output:46f2845b] %[output:6ec4d61f] %[output:2a14a5cd] %[output:5c98f522] %[output:617fc710] %[output:2cbc7bd5] %[output:123ae74d] %[output:54bfb3c4] %[output:3c656e4e] %[output:7ba46003] %[output:7bda5753] %[output:5bfd5bff] %[output:5c5c5e16] %[output:4a47b4e6] %[output:158222fd] %[output:03bbad0a] %[output:9c3360ef] %[output:18b43b1e] %[output:13d441cb] %[output:0fa89510] %[output:70630c73] %[output:6a998272] %[output:7718e7eb] %[output:02baf68d] %[output:91fb0d2e] %[output:03a2474f] %[output:607eb337] %[output:37fd6946] %[output:8628d639] %[output:596e4d87] %[output:41293588] %[output:46854725] %[output:4a33397f] %[output:3c30dfe3] %[output:663bd1a5] %[output:852dae32]

writetable(LLN_result_MK,'LLN_res_MK.txt'); %, 'delimiter',',' )
writetable(LLN_result_LMSlog,'LLN_res_LMSlog.txt'); 
writetable(LLN_result_LMSlin,'LLN_res_LMSlin.txt'); 
plot_10y_in_two(LLN_result_MK, LLN_st,'y');
%%
SS=strcmp(LLN_result_MK.parameter,'Bac3_A81ae31')==1 | strcmp(LLN_result_MK.parameter,'expA_fit')==1 |strcmp(LLN_result_MK.parameter,'SSA30_ae')==1 ;
LLN_result_MK_paper=LLN_result_MK(~SS,:);

SS=strcmp(LLN_result_LMSlin.parameter,'Bac3_A81ae31')==1 | strcmp(LLN_result_LMSlin.parameter,'expA_fit')==1 |strcmp(LLN_result_LMSlin.parameter,'SSA30_ae')==1 ;
LLN_result_LMSlin_paper=LLN_result_LMSlin(~SS,:);

SS=strcmp(LLN_result_LMSlog.parameter,'Bac3_A81ae31')==1 | strcmp(LLN_result_LMSlog.parameter,'expA_fit')==1 |strcmp(LLN_result_LMSlog.parameter,'SSA30_ae')==1 ;
LLN_result_LMSlog_paper=LLN_result_LMSlog(~SS,:);

writetable(LLN_result_MK_paper,'LLN_res_MK_paper.txt'); %, 'delimiter',',' )
writetable(LLN_result_LMSlog_paper,'LLN_res_LMSlog_paper.txt'); 
writetable(LLN_result_LMSlin_paper,'LLN_res_LMSlin_paper.txt'); 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":15}
%---
%[output:9fa424f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: there is a lot of days with less than 50% data coverage"}}
%---
%[output:838b327f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: there is a lot of days with less than 50% data coverage"}}
%---
%[output:7200f62b]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"ans","rows":1,"type":"double","value":[["2007","1","1","0","0","0"]]}}
%---
%[output:5edb9f45]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"ans","rows":1,"type":"double","value":[["2025","12","31","0","0","0"]]}}
%---
%[output:0a79c560]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":1,"type":"double","value":[["5.5000","20.8000","41.6500"]]}}
%---
%[output:01e0b367]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":1,"type":"double","value":[["5.6000","21.0500","41.9900"]]}}
%---
%[output:24552bb3]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJ0AAABeCAYAAAA5WhbMAAAAAXNSR0IArs4c6QAADotJREFUeF7tnWlsFMkVgN8QwMyCF4xhicEGG8ESEzCyjGPEEkAc4gaJGySz2Chc5oplcwQEBnGYG8IN4lS4ERIiBsGvRQbEJa2QEGBxLGLNkXBFMmA7OJ7kFanZntrurqrpnp72UC0h4enqqupXX7+qV\/Velcfn8\/lAXUoCDkrAo6BzUNqqKCIBBZ0CwXEJKOgcF7kqUEGnGHBcAgo6A5Ffu3YNxo8fT+6uWrUKxowZY9o45eXlMH\/+fDh79ix06tQJ9u7dC40bNzZ95sGDB5CVlQXPnz+H6OhoOHDgAKSmpgY8o00zZMgQKCwsBK\/X6zgodhaooHMJdFiNQYMGEcAbNGjgr5WCzk7cXZ6X05oOxVG7dm2iyYYPH66gczkfIaleOKDDF0lOTobdu3dDixYtyHspTReS5nVnpuGCDqWRk5MDs2fPJppPQedOPkJSK6eha926NaAx8uLFiwCjQkEXkuZ1Z6ZOQ4eWabdu3WDhwoVQVVUFffr0gfXr18PLly\/9Fq6yXt3Jim21Cgd0ixYtgmXLlkFRUZHfqOjYsaOCzrZWdXlG4YAOLddHjx7BpEmT4NWrV8SomDt3LtF+OJenNJ3LobFavXBBV6dOHdi2bRts3ryZvEJ6ejo8fvwY3rx5o6Cz2qhufz5c0OFqw7Nnz2Dy5Mlw7969ADEpTed2aizWL5zQYdXPnTsHc+bMIUYFvRR0FhvV7Y9roTOra\/PmzWH\/\/v0QHx\/vX3vlvduUKVNg3rx5pnNw2rVcBR1PohFyP9zQoRjv3LnjNyrwb6XpIgQuo9dwA3QYSbBz505Yu3YtqaaCLsKhU68XOgko16bQyVblbCABBZ1Cw3EJKOgcF7kqUEGnGHBcAgo6x0WuCnQtdOUlm6D64zOo9VUL8Labo1oqgiTgSugqfz4FH37M94vZ2262Ak9BFzoJlJaWQsnfx0Fq81J\/IVEJI6F+6lqyDllWVkY8a9GVm3eFOj2WH+oy3JY\/T+Yi9\/2ajno2dO\/enawLai988fPnz8P27duhpKQEkpKSYNq0aWSGPCoqyp8UZ9CvXLkCW7ZsgZs3b0JsbCxkZ2fDhAkToH79+iL1AVwJOLJpPBSM\/2VfHwQOwauoqCDu3HFxcVCvXj1ufqFOjxUIdRluy58rdIEEBLrKykrisXr06FGgi9H0WQTp1KlTsHjxYhIa17t3b7hx4wZZ5J46dSrMnDnTr3UuX75MgkoyMjJg7NixcPfuXTh48CAgyEuXLg2I5zSqG11+urBnBDSrVwp1mmT4u1Y3NoDb6hTq+ggwxU1CoEM3GtRuCB96rWo13ZMnT8hv\/fr1g9zcXAIYgnjixAlYuXKlPyr93bt3BDjUbtqA4UuXLsGsWbPIbwMHDuRWiEJ35MgR6NKlS0D6UAtUNn+l6bjNqZvA8\/DhQx86DI4cORIuXrxIGloL3fHjx2HdunVEs3Xo0MGfCXZzCCOCNGPGDNItYjeKC9S9evXyp3v\/\/j2BFbdYQG1at25d05oq6ALFI\/shhDp9cJgFPuXJzc314fgItRT+00KHGg01FI7P2L05Pnz4QPzHatWqRaLScR+Offv2kd9ooDAt6vDhw3D\/\/n3i749GgNmFY8v8\/HwS94ndNDu2fP36NTRq1EhoTIdjUbP0VW+uQ\/XHUqj1VTzUjs0gRoFM\/tSQkHlGtoxwpkcfQfxn9+UZMGCADw2Er7\/+mmguLXTUkbC6upqApTUG6D0MGCkoKCBjQfy\/uiJHAvjRo1uV3eB5ioqKfNhFvn371hA6FKPebkGrV68m3er06dMJdFhBVstFThN8WW9y\/fp1EhykN7a2KgnPp0+ffGgc2AFdKCpo9QXV88FJgI6tccjEGnTIi8g8qVHJ\/nk6Wei03StqOTRGFHTBNbAbn6LQbdiwgey3p71wTB0TExN0tT1JSUm6u6ujATF69GghQ2LEiBEwceJEBV3QzeC+Byl0hw4dgs6dOwdU0LKmKy4uJtDh8hJOjbRv357sOtmqVStISEgAkSkTrBTuWumEptOLktJKRLYO2g1q2KY32lGTfYadUJdFCMfGeLErQTQfbIOrV686ugun2dSV7Pux6U27V0xMJ4cHDx7sX31gJ4dxUtlp6PSMm2AERQHClRV2i1dqKGmni2h6\/EBxrKM3LJFpFARqwYIFv1oJonnQd3I6KCcYWYq+Nxc6BGzPnj2wceNGsj1p\/\/794datWwHLYPi3G6CjWtDI2tYTihl07D2j\/LGB8vLyiEzatm0rJHtWY+tpS4R+165dJL8vCjp8Yb0Ff9ygedSoUWTBP5RfBduCZmDRexgATbsqve5T2wWLQMdqNVxX1mpFqu3Y343oo\/V8+vQpbN26FdasWQPaOuNzCBxumo0gnz59msyBOrnJdSjb1BZ\/ulBWUAY6rAc2Fu0O9YBitRKve9WOtdiuldbNSher96Gw74zvFFHQZWdn+9A7BLtRXFvFsU2PHj3I8ha9eK5NstCtvvATPH1bAS0b14N5\/ZKEuiOaiGdIaAf\/It2emSGBZWq7PSPoRMDhaT1W02nTRxx0mZmZPuwmcYuqkydPQnFxMeDmfJmZmeDxeAiMPNcmmTHd0ZsvIefoL7sRzeuXKAWeWfdKoUhLSyNdEabFpb3bt28bDtTNNB0d5NNzJGSg047HECAjq1oE2IiDrqyszEfPLUCNhu5K6I6EXVRiYqLfejVzbZKxXhE4BI9e49J\/C9vGJQtrO56xoB0L0UE9C4Ce9tKzXmm32bJlSwIxejXjWJaO8Yy6VxYSI1jx+S8SOvYUROxq0UWJfpl2z9Oxmg6BQ\/BELx50WF\/0XDayJNkpCjNNx5ZFNaeZITF06FDiadO1a9cAY8NoLk5B97+dH3G\/W\/zqETr0MhBxbZJdkcAx3eWH\/4JubRpJda1azWA0LaI3t6YFmgWJai8zTUc9b4wAYceOvLk3vfp8UWM6raZDh0ucqERTHueHGjZsSL5anmuTk2uvZpqOGjR0DKZn4LCajWe90mkL2lWzZRhZrvR3HE+yBomC7v8nW+N4DrulHTt2wIoVK4gnMXqhInQ8rUJdm2SXoES7VL1GQhj0LrYOelt+aQ+Yc2IZjH4obJf7RY\/pEDhsrOXLl5MoLxpswxs\/sf50TkAXDKhOPoMaDoOQlixZEnAKIna5uKRotL7qZB1FyqIf694lGfCHlHgS9O5pNZU8annBv6KiwofaDb2HsZvEeAcaVigy+4+Tlk52ryICC3ca1miQXbEId\/2xfArdrhk+SGvzuUYV33wPlc2+J+ECllyb8vLyfGfOnCFzc7h+qnXOE42RkDUk3CDUUNeBnaYROTM21HWSyV8Put\/EDYOojoXWNV1KSooPu1X0IsHJYPaye8pE5sVV2vBJgEL3p5Zj4duYbyA69h\/wXU4V2WnB6uW5cOGCr2\/fvrrAYeZuc22y+sLqeTEJUOgGVSyHuP98Dj39LucT9FicJpaBSSpPYWGhTy8WddiwYYAn87nNtcnyG6sMhCSgB13KmBgY\/NcEoefNEhm6q2stUTe5Nll+Y5WBkARCCh27DCZUIyaRrJdJMGWoZ5yVAG3Tb3t1hkYN48D7NhZGRGfZo+kUdM42Zk0pjUIXOyUaolp\/3pZt6L+TIad9nOWNKmucE2dNabSaXk896JIftoX+N7oQS\/aP+c2C3qiyxkHHc+K0siqityRGz\/7Si31gl9iszsXxosJ49+0EXQ+6+GvdIOVvk0gxXbJ+gl6Fw4IqssZCh2\/LxgwEO7bUgsxCa+QxIuP2LtIyPM8U3n2RMmTS6EHnfdMEYh\/8jozvejxqA2P2\/DMobRdR0PHWio2ELuKDp\/XRM\/IswXyOHTv2qx2uzBqbFxXGuy8DkkxaPei0z6PWm9ZwQlDzdhEJnUw0mEhQDabB7RXQGQI3CDLyBDbzENZrcF5UGO++DESyaXnQYX5ozf65\/wzZrCGioAsmGiyYLtko4MfMN4\/XMjwPYt59Xv6y90WgG5A4BBamF8hmHR7o3p1YClWvnkDtpokQM3qJVKV5hoRsNBiva9WrnBF0VrxJeFDx7ksJUSBxREFX9sMBeLUty\/\/aMaMKpMCzOxrMCDo6cNe2Dw3okYHOrqgwBR3zpch0UQgcgkev6J4ToWnOfoFv73MSnrEgGw0mUne2wUW7VzujwtwI3V\/SC2Bg4hDhtqMJHR\/TsZoOgUPwRC8edLzukp16EDEk2AYXMSSwm7czKixc0PVZ1NO\/IqFto9SmaZD9+ymizRaQznHosHQc05Xf\/QG87XtKda2img41EbsxN31rPWh5oLINbjR2Y6dMZObWeFDx7gfV+iYPifQAwZYZFuiCrSwPOjZSS09wRhYmHXuxqwraVQftxDEFiv5mlK9dUWEKOgtjOivAaaGzKxpMWx8WEHrPaGlNdhnMSlSYgi6M0FmF1qnna3pUmOpenSLF5nJqclSYgs5mGJzMrqZGhSnonKRElUUkoKBTIDguAQWd4yJXBSroFAOOS0BB57jIVYEKOsWA4xJQ0DkuclWggk7DAM+JUzYajLelv9H5YLJLYDyMeZFeTp8PpqDTgQ5\/siMazMzFnN1d3ev1+uewtMcyWXFTxwx53igUACePalLQCULH87XT0zY8YFjhOxkJhvUN1\/lgCjpJ6GSiwUSh4x1gYnckGAUuXOeDKegEoUNB2Xk2mLbro2NFUVd13hhOe1\/EbcnpU3MUdBKGhGw0mMgpiPTYJxzTyQTliIKnoBOVlCad7FdRXrIJqj8+C2r3H7ujwXjWK3sOqwx0dkWC0a7WyZMQZdtUBhvH3dUrfz4FH37M99fR22621H4YPGNBNhpMT9NRYwErycZaiHavdkaCKeh0kJb5KhA4BI9eUQkjpTZP5kHHC7IRPRvMCLxwRIIp6CxCx2o63K0bwRO9eNDZeTYY\/Zi0XWw4IsEUdBahw8dxTPfp9XWo0yRDqmslz5aXGx4dFUw0GG\/KRC9KzOlIsEiD7r\/6ih71UAe7pgAAAABJRU5ErkJggg==","height":76,"width":126}}
%---
%[output:5f489c33]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJ0AAABeCAYAAAA5WhbMAAAAAXNSR0IArs4c6QAADphJREFUeF7tnWlsFEcWgN84NsaAxWEuAwab5XLCEbC9mA0JZxYWBEsguyaImwSDbAwBzK1dYLnPcN+nxCkQPxBIIJC4ZYIAISMgECACzLFg2OUyDsazesXW0FN0d1XP9DHjqZKiCHdVV9Wrb169qnqv2uV2u90gk5SAjRJwSehslLasikhAQidBsF0CEjrbRS4rlNBJBmyXgIQOAHJycqBPnz5E+LNnz4a0tDTdgSgoKIAJEybA\/v37oVmzZrBhwwaoVKmSbpkbN27AoEGD4P79+xAdHQ2bN2+G5s2be5VR5unWrRvMmTMHoqKibIfC6goldA5Ah4PatWtXAni5cuU8Yyyhsxr3AHq\/3ZoOux4eHk40Wc+ePSV0AcSCbU1xAjrsXGJiIqxduxZq1qxJ+io1nW1D7nxFTkGHPc\/IyICRI0cSzSehc54F21pgN3R169YFXIw8ePDAa1EhobNtyJ2vyG7ocGXaunVrmDx5MhQVFUHHjh1h4cKF8PDhQ88KV65enefC0hY4Ad2UKVNg+vTpcODAAc+iokmTJhI6S0c6gF7uBHS4cr158yYMGTIEHj9+TBYV48aNI9oP9\/KkpgsgQKxoilPQRUREwIoVK2DJkiWkWykpKXDr1i3Iz8+X0Fkx0IH0Tqegw9OGvLw8GDp0KFy9etVLJFLTBRIhFrTFSeiwOwcPHoRRo0aRRQVNEjoLBjqQXqmETq9dNWrUgE2bNkGtWrU8Z6+8fqSnp8P48eN19+CUZ7kSOp5ES8hzp6FDMV6+fNmzqMB\/S01XQuDS6kYgQIdRA6tXr4b58+eTZkroSjh0snv2SkC6Ntkrb1mbjJGQDDghAanpnJB6iNfpgQ4N2a1bt8K6devItkD9+vU9osFnp0+fhmXLlsG5c+cgJiYGBg8eDP3794eyZcuGuAhl941KwAPdxYsXYeDAgcTVhoXu1KlTxO+rZcuW0Lt3b7hy5Qps2bIFvvrqK5g2bZqXy7XRBsj8oScBAt2zZ8+IIyGCV758eS\/o8BkCh9pN6dN\/\/PhxyMrKIn\/r0qWLIckV\/PITFL\/Og7AyNSGq4ShDZWXm4JeA6+3bt248dD579iy0atUKdu3a5QUd7mHhNIp7SO3bt\/f0+OXLlzB69GgSBYUuOqVKlRKSRuHdPfDqYrYnb1TDkRI8IcmVnEyukydPujGcDu2169evk\/8rp9edO3fC8uXLVe081HJo44mE4N27dw\/27t0LnevmQBV3jkeCkXHfQtnm88m544sXL8j0jq7bIsmuMtgWu+rypR472ycyLrw8rrS0NHebNm1g2LBhsHv37o+gmzt3Lpw5cwbWr18PVapU8XofPsPDaozhTEhI0K2L7vpvntsLGkfu8eRF4BC8N2\/eEPft2NhYKF26NK\/d5LldZeysy5c+2dk+oYHhZHL169fPjf5cFStWJFMrq+kQLARGTZup5deqj0KH2vHPjW7D2yc5EFE51TO1orDRXbtatWrCAcZ2laGDKtv3fnRxJhKdjdR4cOXm5robN25MntkB3YABA6D572nw6oEbysa6oMn3kaTuwsJC4kGL2jQy8v3feMmuMrJ93iNRoUIFoqR8Ta7i4mK3y+XyCTpfptfsrkvg55cXoaBSPkQ9jYH+DYZAqx9jSHTUo0ePoHr16sLTq11lUDh21UXrKf\/fXRBR\/Bg+KVNLaKFlV\/tM0XQJCQmqlyLSOzoOHz5M3KnHjh3rCQqmhG\/btg2uXbtGfPtxAaCX0EM2Ozsb6tT5HF52vunJ+kVYJ+j9ZU9iqD958gTwVyRq04mUKco\/C8Wv70FYmVoQHtNSqB62DDZWpC62\/76WefbbESj97y2e16HNWyqul658fa1LT+boN4j\/mZ3I6pW+9MSJE7Bv3z4CGMZmNm3aFBA6dLfB1adMoSUBPAzAsTcbPK+zVzWb7siRI8SHHyun1x+EluhDs7e4b4sz3Pbt2yE1NdVUIXCho6tOKyo3tSfyZaZKwMpxl9CZOlQl52W2QacmMisrLzlDVPJ6YuW4c\/3prKxcOVRqEVHK58rp\/enTpySIBT1eeLdm8nDAbZ81a9aoZlO7lZNtp+hNnFrtwEtzcPU\/b948L3cymp\/2FSPKzLat9GRj5bgHHHQoCPbaU1YAZkOnduJCb1AaMWKEB2wKHIYiIgSY9E5seMDTfuCmOOtOhmWVgNttU1sK3cqVK93YIdxHQ\/elHj16kNUqPWe1snI1TacGHRU+fYb\/NlPTaR3zsUCpre59\/QEoI9BoPK3ScVZ5bRj2u0RBl5yc7MZLnlu0aAG5ubnEObN27dqwdOlSskUSSNBRDUMHukOHDnD06FG4dOkS4VctbI+dPtk8epoKn+FlNlTzsv+mPxqtv2tpOuU5dJ06dci+qFLTUeCwrXg9LF6QvWDBgpIzvR47dox4mdBEg3779u0LmZmZZFBR81n9S2O1mfJWcRwkHFjqdEChQ9io3aVWngXKyPTI2lpqZZXQaWlL3hSL5VjolGUogCUKuvz8fLfyOvpXr16RKxPCwsLIL3zHjh0wY8YMYejmHroNd56+gdqVSsP4TvruTkYWEkqDnUKHhjW1rfBdygHC\/NgPpf2lBoDeQkI57elBh9Mu+h2K+BWybQhJ6Nhvgz1\/\/pxoOPQiQOgWL15MfOlENN2Ocw8hY8eH24fGd4oXBk9P01GYkpKSSJu0bDoWDIRh4sSJBDw1Q11vIUDBxjwIE2peLYhZ6Fh7TC9aX0IHABj7gNDh4Ty6qWPsBH6kQwQ6BA7Bo+m7lOqw4rtE3gxDnutBR+HAdiA8uODBhQS7jaCmjVgAWAB5\/oIILfZdT3Mq34F3y7E2mJ7NF\/LQobv68OHDycoVFxJVq1Y1BB2r6RA4BE8k8aBTrhwpdOw+nd4UiG1Q26LQg45dRIksJNDOVNqfdNrX2osLaeju3r1LAm1wYFatWgUNGjQgrBjRdEQjHboNp379D7SuV0F4ahXVdNRYx\/yo6fRsOq2NVDWQRD2j1Ww3dsuEt\/cWbDbdxo0bP1o1++05jDYd3gKJcL17944E4eD9tzQZhU5Eq6nl4a1ecVuHrlSVq1c67dO\/4XYP2n2Y0AbDpNxsZjWblqZTblvQxQpbB9p5WuWpPYn169mUga7pFi1aREwLZfLbc\/jEiRNu1HDx8fGAFcTFxXlVMGvWLOGFhK\/AKTUd2m1qSe0YjN2noxcQ0vJqR2vssZUdx2Dslo+yf4EO3ZrJSZCcnOTlvey3pmvatKkbnfVmzpz5UbQXCgfVq5EtE3\/AKwllUcNhUp4Jo4YcM2YMTJo0SfV8NRD7Tc2QNZluSKoHYGZ8siszM9ONwdJagRbUiVNk9RqIwrO7TWoH+P7s49ndflofhe6H2r0hqZ4bvsgoIvHJZiRXVlaWG6dWNuE02717d7hw4QL5FqqETlzc7DaNv54o4jWbl5NC16B9MsS+awwDE7pCm38kmVKBSyswh25o4haAhM4UWQfVSyh0MenREFk3HHqEfw1jv3m\/QPM3BYxrk78dkeXNlQALXYdPCmHCZ98LhUPyWiKh40koRJ+z0P3xWXn4V0qSKXadhC5EoeJ1m4Wu0o1GMKZotCl2nYSOJ\/0Qfc5Ch2Lo\/nsiZHwa6\/e9ghK6EIWK12016BJ\/rQ+df06F6JhH8GV2NZ\/tOwkdT\/oh+lwNuqj8yhBzoxG5g6ZftWbQfs5ffZIOFzq7NoedigZTSk3tSIw9WqP57Y4K40WN+TT6OoXUoFNmb3PzDzDlmz\/5pO240Nl1DCZy4M8e7psRgoiCpAJmnS21PEaMuL2LwMDzTOE9F6nDaB4edP4sLLjQ2X3gj8JhQxCtjAbTcn2ng6T23M6oMF7UmFGYRPPzoMP39IoeBD92zhR9pScfF7pAcG1iNYuZ0WAit4liHrwKrV27dsRtXcSZUxlYpDUqvKgw3nPDo22ggAh0f4nvBpNTphp46\/usQQGdVdFgPG9lNWk6FRXGc4EyPPKcAkEF3bPd06Do8W8QXiUeKv79n8Ky4C0krIgG402tRqHzx5uEBxXvubCgBTMGDXQvjm2GxysGebpV8W9ThcFzIhpMCzqlZzLtDPtVa7XQRiujwgIRukkpU6FLfDdBjD9k406vRhYSCByCR1N024FQJWOTUKN4Ux3aUWZHg\/HqpA1XDjj9lLoadFZGhTkFXdv0eAgr9+EW1vOlPydiGZHaAQZ\/li40tmwmLnRGtkxYTYfAIXgiiQeAVdFgIgsJdsBFFhJmR4U5BR1uUzW8cwgKrhyDvMJEmF7rB8NBV4ahM7o5jDYdNjDq07bCUys2igedUotgfrOiwUTsOnbAnYgKcxI6s68o42o6alBa7TkssjlsRTQYAkw9fTHeV3k1hHJxo9w4diIqTEInMl8azMNbvVoVDaZspjJskP5d60oIX47B2K0fZd08qHjPDYqbm91KZRMwmo4rhSDLEOxRYRK6IAOOTtnsVRL+7OPZLQIJnd0SN6m+YI4Kk9CZBIF8jbgEJHTispI5TZKAhM4kQcrXiEtAQicuK5nTJAlI6EwSpHyNuAQkdOKykjlNkoCEziRByteISyAkoPPlGMyMwBy9SxFxiALh+2B2R4Jhv0MKOuxwIHwbDNtBz2KV576hEAkmoVNxe9L6joT45PEhp97t6mqeL6EQCSahU0Bn97fBKHR37tzxuDyJOHAGeySYhO7\/9oXd3wZDwWvd2M5zVVd+9kpE+\/LclnjPReowmiekbDqt29WtiAbDgeB9BRE\/QExtOr3wQ388SHhQ8Z4bBUokf1BBV\/DLT1D8Os\/wdVJORINR6LS+bM1+\/8EIdGZ+H0xCp\/MzKby7B15dzPbkMHINvEiMhNnRYFqajq5a2e0S0UBrs78PJqHTgQ6BQ\/Boioz7Vvi6UB50VkWD8b54w4InspAI9kiwoFpIsJoOvzuA4IkkHnRWRYNpQUfbc\/78ea\/PdoZCJFhQQYeNRZvu7ZOzEFG5paG7y5yKBtNbSKhFiYVCJJjV0P0PvFyc9Sc9FAkAAAAASUVORK5CYII=","height":76,"width":126}}
%---
%[output:5d98ae53]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJ0AAABeCAYAAAA5WhbMAAAAAXNSR0IArs4c6QAADpFJREFUeF7tnWlsFEcWx58RYAyxMLBouWUjjnUIIgQQ+cASBEhoieADx0KQDIuJMIjEIMQNYgOIcCYhcojIwfmBw1xiubQfogQwCYfIBikCkSyEG0uAQTKHlyWe1b9WNdQU1V3VM9M9PZ4qyZLtrq6uevXrV6+q3qvOikQiEbLJSiBACWRZ6AKUtn0Uk4CFzoIQuAQsdIGL3D7QQmcZCFwCFjoHkZ86dYrGjRvHrq5YsYLGjBnj2jlPnz6lefPm0cGDB6lHjx60ceNGat68ues9v\/76K02cOJFu375Nubm5tGXLFurZs2fMPWKeYcOG0cqVKyknJydwUJL5QAtdSKBDNd5++20G+CuvvBKtlYUumbiHvKygNR3EUb9+fabJRowYYaELOR++VC8V0KEhhYWF9OWXX1Lbtm1Zu6ym86V7w1loqqCDNKZNm0bTp09nms9CF04+fKlV0NB17NiRMBm5c+dOzKTCQudL94az0KChw8y0X79+tHDhQnr+\/DkNHjyYPvroI6qsrIzOcO3sNZysJK1WqYBu0aJFtHTpUjp8+HB0UtG9e3cLXdJ6NeQFpQI6zFwvX75MkyZNort377JJxZw5c5j2w1qe1XQhhybR6qUKugYNGtD69evp008\/ZU3o06cPXblyhe7fv2+hS7RTw35\/qqDDbsOtW7do8uTJdPHixRgxWU0XdmoSrF8qoUPVjxw5QjNmzGCTCp4sdAl2athvF6Fzq2ubNm1o8+bN1K5du+jeq65tJSUlNHfuXNc1OHEv10Knk2gduZ5q6CDGn3\/+OTqpwN9W09URuJyaEQboEEmwYcMGWrNmDaumha6OQ2eb558ErGuTf7K1JTtIwEJn0QhcAha6wEVuH2ihswwELgELXeAitw9MO+g2\/FZDt2tqqU2jejSloJHtwTSUQFZxcXGkoqKCsCb02muv0fvvv09vvfUW1atXL9ocbMUcPXqUPv\/8c7p06RIVFBTQ1KlT2bpRdnZ2NB\/KOHnyJJWVldHZs2epRYsWVFxcTOPHj6cmTZokLJ5\/VD6jxRefRMuZkt\/IgpewVIMvIKuoqCgyevRognfD7t276cSJEwS\/rqKiIsrKymIw7tmzhxYvXswCRgYNGkRnzpxhWz9TpkxhkMKtGgnwwtW6b9++NHbsWLpw4QJt3bqV+vfvT0uWLImJcnJq6s2bN2nv3r00cuRItrUkJgAH8Hga3qohLS1s7Cg1vCzV1dXME5fX0U3EXvOjLK\/3pHv+ZCCaVV1dHeEhbxDIhx9+SMeOHWNxm\/n5+XT16lW2FTNkyBCaOXMm6zyAWF5ezvLyWM0HDx4w4KDdxDA6lFVaWsr+N3ToUG2d+U7A9u3b6c0334zJL2s6AAfwnFJNTQ1z\/27dujU1aqQfir3mx3O93pPu+bUdaJDhJZsO2grDIe\/0Xbt20dq1a5lmw\/DLEzoTMAKk9957jwAL7sO2zcCBA6P5Hj16xGBF4DG8Yhs2dIYEN7lBhw775OIDuky51DuvvnZo9buDLXQGhCmyvAQdXKUxZAI6DJPQULDP5Ij1x48fM68K2H7weD1w4AB99tlnDM7OnTvH2HlOZaiqrIPOT83lFVILXRKgg1aaP38+Xb9+nb744gtq2rQpA6u2tpaBJU4GuOsN3KgRp4mf77\/\/nr7++mtq2bJlTG1WrVrF\/MMwFGMS4pb8gO5ATTOqoobaGa+FTm+CxIdZ7F1M08F4v3btGu3fv59prHfffZfNYJ89e8ZgQoInqzw07tixg00W4Md\/6NCh6O8w3MX07bff0r59+xjAPIjYqfLwmp09ezaL+4SmFRNszsrDFyn7fi017tScmvRv7yoD5P\/u5kPa+fCFMGEDOtmByH\/v3j3Ky8szsgHxcK\/3pFP+Vq1aEX7kBLveZGLm1DlZN27ciMyaNYvNSG2yEhAl0KtXL2aPy+DhpWzWrFncwso6fvx4ZMKECcxnS6eF4n6KvTHtJHD69GkWHLRt2zbq3bt3TP0T1nQTJkyIHD9+PDpbTTvp2Ar7IgE32zrRB2Z169Yt8uTJEwtdopKsY\/f7Cl1ZWVnk448\/ttDVMWgSbY6v0JWWlkZweqRqByDRisv3qyKcxDzx1EE8YEZ+ntOJmHL8g8lJm26ywJIQEiK8VAkL7FhOSqdTNH2FrqCggH1HIp4O9wolhw73yR0QbyM5dFjQlo9oBQwoV1zYxt+YrfNFbLf7TdoHoLC2ycMK5Xt4u9ItqCbe\/jCRWdYPP\/wQwdm6qYbODUi3hrhBI1+rqqpiW3fY0xW1EsDZuXOn0TnBvC6y1lZBB+ixyI5koXvRi6GDDsHLMhDQJDzx4Ga+1WYCHfaOARrPy\/\/mZTr93wl2Dhx2brD1t3r1apLrDeBgtkCjYmEcOzeiducvwPDhw+mnn35ieVVwqkwSHeAoRzYtUB\/U4fXXX6dly5axZ8myFNubEZoOjYRgxKFQpYHkIVM3vIq2ljy0ytCphmjdcMGhkKET7+MdroLu\/Pnz0dPbZW2PMrCLI5ojHFa4jnFzQi5fNWpwrStqXMgXvo\/yfjme5yt0FRUVEfjOeRleV\/3zN7peVUMdmjeiuUPc91JF4esmEuLbqRIuypK1kttEAvlFreAEndOzdMDheqLQyUO92L4uXbowc0AETK6Tk5aW26qyb93q7it0o0aNipw7d84Yuh1nK2najhenCc0dkm8MnpvdxoWHrRfVLE+0jyB4\/pK4aTpu5PPZqVfonPKrXqR4NZ0MlGx38nabfJtCfqnF4VOlbdEOFYy+a7pOnTpFfv\/9d2PoABzA4+mdPq1o\/TuFJkohqhXE4UIehrgdBJuNQyPaH\/gdH\/yQ7TTV0Mg7sEOHDgxkDGXizNVteEUHLl++nGXBaZhOHy9JVNPBfhWdVVXlyUs8KnuNT1j4tV9++SWmrU7LOk6TKF81XXl5eQQNNx1eZU0H4ACeSdLNUEUbg8Mlw+Q0vKqgk58HbxoRWBk6cYKB58DVHq7869atY+f\/qr6Akyh0Ok3nNJzCfQz2L+BSrT6ohld5MpMyTRfPkglsuop\/P6R+nfKMh1bR\/nHTdHxdzUmYXPuZDK\/yUOVku6nedvwPCRoOWgKuXrL7vN82nep5eKb4cv7444\/KyYA4e8aoESqbLh7oTLSaKo\/J4jC3v7hGw2yLL6GIkwYT6GTB8w7DEozb\/XxohfcNOgwvAmI9VDsOiWo6DPm8LrI5gLLldUVZhrgfmk7cVeFDo2zTYQgWJ1Yq+fB+42Vs2rTppZctYS+TVEDH16RkMOUhXp6Zwl7BEAi7jM\/6\/NgGU5XpZMgnCh2i67755htmbyLJa3AcRH4deeSFZtnmw3WEfiJQipsdfCIBB1vIGcltcsLLxL488okpYX+6IKGLV0MGfZ\/K6Mb\/EB2n+xqiaV0TWaYxfYaYz2n26lQWh84XfzoLXazYOQzyrDLZm\/bpAp3pBNPLixDoNpiXiqUqr2pnBHXxulWmq7+FLqANf11H2OvhkYCv63R2eA1PR4epJha6MPVGhtTFQpchHR2mZlrowtQbGVIXC12GdHSYmmmhC1NvZEhd6gx0OifORBciVdtXJi7ZnCO\/o8J0UWNh4rnOQQfhJisaDGWJMMvgOkVrBR0VposaCxNwqAuHbuobn1BhmzeoafuG9OfZf0xKNQPdkXDzMtH52rm11s3XH\/fJ14OMCjOJGktKTya5EA5dl4G9Ka9pa8qpakHFr5YkBbzQQec1GswJILEPkAfeEjicG4cEBRUVZhI1xuuyYMEC+uqrrxy9TYKOCuPQtSjJpeyO\/z9TumNNd+r6p3z2e+vGram4W0lcqIcGunijweKxPVIRFebkAqXyEZRfJNUo4HdUmAo6mbDiVyfHBV5KoHPyp4snGkw3tKpeRa8BOiavs86vTged7G4v1hGHhwcdFWYC3V\/yh9HCPh+YiCcmT1zQPShfQs\/vXqX6LfOp2V\/\/bvxQP6LBnKATg3p4BbmDpFfokhEVpoPOKQBcdMKE529QUWGhgq76uy10d\/3EKGjNRn9gDJ5usiC7T5tEg5kMr3KHexlekxUV5gYdjs\/FKQHiAeGq4VP2EPYzKswEugV9PqCh+cOMlQ7P6FnTATiAx1PugL9Ry2mbjR6sgy6eaDCTiYTc4V4mEsmKCktU08kC5m3wKyqMQzd40YDoREKsQ8+WveKy51CGZ+hkTQfgAJ5J0kEnRiyZRoOplkTkusgdnoqoMB10bjadqAHFtvkZFWYygpj0uSqPZ+hQCGy6pxe+o5xXBxgPrbjPj2gw3igeCS\/vKohDkrhwbBLKmMyoMB10aId8fBmPhFNpc7+jwkIHXbyEe90GM4kGE+uiipzCdaftNdlGkoFNZlSYDjrMThE8jYBoJLkuQUeF1Rno4oU1VfcFERWW7NgLnaxMo8IsdDpJ+nA9qKgwC5301UEf+jJtilTtkKDyyYYk2eXpBGw1nU5C9nrKJGCH15SJPnMfbKHL3L5PWcstdCkTfeY+2EKXuX2fspZb6FIm+sx9sIUuc\/s+ZS2vM9B53QYzkbjuSP+wfB8snSLBIPc6Bx0alaxoMLcj\/eXjVHNychjHsj+dWxkm4OsivXTXTZ4RdJ6MgE7n9uQkdB0wsvBsJJgZvhkFnddoMFPo5AO0\/f4+WDpHgmXU8Jrsb4NBeLLfnBdXdTOdoP9Uk86tCe5M8gnr\/CDvVESC1UnokhkN5qbpVJ9+8hqUYwJeXYsECyV0Ty+to9ont6he47aU03WGSb+wPH5Eg+lmr\/IR+V6hy8RIsNBB958be+jxv2ZHQcvpOt0YPN1kIZ5oMJWm45MFVFL8lKdq5soboionUyPBQgcdgAN4PGW3H0VNeq4x0nY66OKJBnMaXp3As5FgsR87RsepPuEUqtmrrOkAHMAzSTro4okGc7PpuODEIdZGgp2K0f5O9miooGO22aV19N97p6nBH\/oaD606m443MpnfBuNvMSLjxUAXGwlWEv3OmdP3wfyE7n9uJ9sTeh9JuAAAAABJRU5ErkJggg==","height":76,"width":126}}
%---
%[output:3874273a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJ0AAABeCAYAAAA5WhbMAAAAAXNSR0IArs4c6QAAD79JREFUeF7tnQtsVUUax7+rtbSRIsWAFAWBaBdwDasLBQwIa1QWDAhdWXyhCypIihZEaRUID5cAQR6KRSCUV1jAEmAJKppo0IguhGWTbhEEawXU8lIgVqSElrv7HzKXuaf3nJlzz7PtTEKAe+bM45vf+Wbmm5lvItFoNEo6aAn4KIGIhs5HaeusmAQ0dBoE3yWgofNd5DpDDZ1mwHcJaOiIaPfu3fT4448z4c+ePZuGDx9u2RAXLlygwsJC2r59O3Xt2pWKi4upRYsWlu988803NHLkSKqsrKSMjAxavXo13XXXXXHviHEGDRpEc+bMofT0dN+h8DpDDV0A0KFRH3roIQZ406ZNY22sofMa9xCl77emQ9VTUlKYJsvNzdXQhYgF34oSBHSoXOfOnWn58uV08803s7pqTedbkwefUVDQoeZ5eXmUn5\/PNJ+GLngWfCuB39B17NiRMBk5fvx43KRCQ+dbkwefkd\/QYWbau3dvmjx5MtXU1ND9999P8+fPpxMnTsRmuHr2GjwXnpYgCOimTJlCM2fOpPfffz82qbjzzjs1dJ62dIgSDwI6zFy\/\/fZbeuaZZ+j06dNsUjFp0iSm\/WDL05ouRIB4UZSgoLvuuuuoqKiI3nzzTVat7t27U0VFBf38888aOi8aOkxpBgUdVht+\/PFHGj16NB08eDBOJFrThYkQD8oSJHSozgcffEDjx49nkwoeNHQeNHSYkhShsypXmzZtaNWqVXTLLbfE1l5l9RgzZgwVFBRY2uDEtVwNnUyiDeR50NBBjPv3749NKvB\/rekaCFxm1QgDdDg1sHTpUpo3bx4rpoaugUOnq+evBPTWJn\/lrXPTZyQ0A0FIQGu6IKTeyPPU0DVyAIKovoYuCKk38jx9g+5syQyqOX2EUlq2p8y\/TmvkYm\/c1fcFuqpPV9PpopExSWcOm67Ba8Tc+QIdgAN4PGT0+xu1zFvF\/ov1xqqqKraDFlu2ZUHHt5aQ1\/KRtY\/Kc1+gM2o6AAfwEKqrq9m27aysLEpLS5OWWce3FpHX8pE2kEIEX6BDOcoWPkflH62jfqMK47pWr4Wk03cXUgWmpFF8g46vb65fv5569uwZK5iGwl0ovJanlCiFCJEOHTok9E8nulfAYvQXX3xBixcvpr1799KNN95Io0aNoqeeeoquv\/56hWyuum7Q0NXv4YRSY0siRXJycqIPP\/xwHZ8ZOK2E7dMIu3btYucze\/ToQY8++igdOHCA1qxZQ\/feey\/NmDEjzjWCWX5ONJ1obkkfXJD0GPCnn36iH374wVIkGIgjXvPmzZXGmA05PvYN4o\/bITJixIjo22+\/Tc2aNUuY9tmzZxlw0G6i743PPvuMXnzxRfbbwIEDpeVKFjrjJCR9UAFV9xkTm3jI7H+8u6mtraWpU6fSnj17pGXVEa5IAEoGW63cBi9SUFAQxVG41NTUhLIGLOhGsdfrvvvui8X59ddf6aWXXmLeiqze5y8kC53R3JLSczjV5s5m0F3avTGh\/W\/uR9\/RsTPV1K5FGuX3zWKa8ejRo2xIACFyNw4aLnMJ4OPEgSHjcMgNmUUWL14cHTdunGlaGzduJGhCbNO+\/fbbY\/EwzoOWwxhPxVWWHehEaEaVz4yz8YnQnZkzgKq\/+jRWJphhPs6ZTaVLCymr9gT7vXf331OLgXkx6LwQohsNEbY0zNrLjXJGiouLoxjnbN68maU3YMAAeuGFF2LaYO7cufTll1\/SihUrqGXLlnF54hkOlcDXWocOHSzLowrdhr0nYtAcv7Y101Z\/PrIsljaH7sy\/tlHm1glxeQK6fU26UvZH8b\/TA\/l0LHsI03QaOjVsPIUuOzs7Cmd9vXr1orKyMjZBwPjunXfeoezsbAJYKEAibfbuu++yGa1RCyaqlip022ePpzv+c+UcaKIA6P7ZeTJlrR9Bf7xYWgc6\/CCufrAI3f5Cx3qM09Cp8cZieQrd559\/HsVMlYdDhw6xAyL33HMPvf7667Ro0SJXoUOXjNkyDxjow4fHTTfdRLRvM51d9pylaFbeNo0tnY0+MqtOvLQ7+rHfxC6X\/dCxB1U88Hd2vtRrTcdPdiHbRJ40zZ6LznPwLj9FZoOTuKhQCOihzLx5yp5z6FauXBlnV0UmWK5UWbI0K3sd4zAaFLDBLgcNhkYy03TJdK9PP\/00m5jwcPHiReZWAV136rYpRP++0s3zUJnSmtrUXBmfIaD7zKo9GfebSsOUtupHE3ccDSV0HLg33niDNfCZM2fYh49\/4\/ii3cCBMTvcI3uO\/HicBQsWMBe3YoA5KTMz026xYvEjly9fjkYikbgEABP86QK6ffv2uTqRWLt2LXXr1i2WH778kydPUuvWrenCmnF0Ydc\/kq6M1YulKbfRxD21oYMOZYb\/YgRRK6HRX375ZaWhi1hvtN2yZVfGwImgkz3naXHojO2F54413alTp6LiBAGaB\/as8vJy5iUSf\/tlMqkqHlt3POYSgl+ld6L8XRdDBx0+Omg1GN1FB9tc2xl\/txKHqCy2bNnCHPGIIMuei2kbx3SiRaGgv\/WkUdZkkZKSkugjjzxCXNvx1Ydnn32WYEo5d+4cMw7DLiauPnhhHDba3WSFt\/O8tCqNJh7OCh10sBxgIse7Vl4np10sADNCZ9SIVs9F6L67tn2cGapfrz9Q1zFz7Ig\/Lm6Ez15zcnKotLSUNQrcVkEIrVq1YpF37NjBfG1gwjFs2DCm\/bxaBquc\/qe6E4Gkq3f1RbvQJftl251ImEHH04Eri2TGdW5Cd3rnxjoWBScbcSNLliyJAjR4D4KlHvcpGBfyEy34P\/HEE8wEYbZ8ZuRE1WRiXIFwgTeWhB3oYCvM23DVi1JB\/\/ak2qV4CR1mnK+++mpMJFZ3XrgF3fLR\/anjvqV1mkHciGu3jUK3tSkM0AE4gMfDY91bU9FjnZVk6xZ0xu4VwGF1iNtLZWM+t6Cbn32cumZUa+iUWt8QyYmmA3AATzVYNbgRJjN4jL8jTQSxq7WytWnoLCzcxk2HpcsKKePjKwJ2M9iBDvliTLer\/Bz1vq25ctfKy2u1UmO0yZmN3YwmExXbmp2JggxKnp\/WdA4otAudg6yYu37Y3o4dOxa3fMi1V7t27erY5DCW5mM0s5mr0Y+dlWcnGVSy5zLosPrTZvrOpMSkx3RJiU3tJdEQy98wG\/wnswyWqMvl+cigkj2XQYd8kp3BaujU+Ak0FjTcrFmzCEuI4vYygLFhwwZPbktUgS7ZGWzooPtqy1JK3zDW9Ub2s3t1vfD\/v0fMOGlwaseTlbFRQffhikV198PJJKTwvL5Dhyoa7XROd6JYia1RQRcGO50Cww0+iobOhSZuCJrOBTEoJ6GhUxaVeUQNnT0haujsySthbA2dPSFq6OzJKxYbu4z5OQoNnT0hBg6dn24lrCYScKgIx4oqAcCVth1Co8pnsOgaOhWpXY0TOHR+upUwgw7WbwBX56SXiSyntyigQec\/1JrOHmux2IFC57dbiUTQ8XU+uJA4u2m6khjFrtVPTWd3axNuQhQDX3fFbhLRu5VSpYVIstNesucq0Il+Bu2UT7oigcz9OiMBp4hG3yWoDF9ucWLD86t7dQKduKDv5KikbEeK7Dlkngg64\/DGs7XXINxKGOHiX1QiIM2+MGhH8fxr2KEzLvgnC53stJfsOZeniqbzbO3Vb7cSvNLoSi8c+JTSu\/SL89xpBNIIF97nXyDSqPrvx1Rz6RIdqsmkvK2HQncwB90rBw5blXJzcxMe1FHpvmSnvWTPxTwCh87ssLUXbiVkwsWmz8q1hZRWWRoDkgOKd42Q8k2i3GtTshpEVi7+3En3ijSMmzxV8zXGk21dkj1vUNDl5+czv2c8wKPAJ598Qp06daJbb71VKuNk4zdt2pSKiopCqenESmvosG3bwoGOHbcSOGr3yiuvBO6UUFXTyZwtmn0dXmq6IE6DmW1XFyd4Uk1hiCCdvbo1kUC+AE\/mftVuBVTj23Hy5+SyFa+gC8NpMOP4OfQmE1U4vIpndu42UX5Wl62olM\/OaTBjembda1CnwYqG\/o4dQeRjZbMJnopceByppuPGYaduJewUyou4dqCzumxFpWx2ToOpQqdiWxPTkk0UZM\/tyEtFJmIcKXSI7IZbCbsFczu+XSE6+aLtngZTnUgEcRpMdQxsp72UoHPDrYSdQnkR1y50bpTBzmkwnp+d2asfp8ECg86NBgg6jSCgc6vOQZ4G09A5aMX6DB2qHdRpMA1dI4aOgyd6bfLjNJiGrpFD56D6tl\/1smdQmkjYLnGCF86fP0\/wX8t94eHapyFDhjCP56L72WR3KWOjKS5ZWbJkCbvFRwxIE\/dgwKWDF1+uG\/IJWxr1Hjp4T8eaK2ZmcBRz9913x+6sgDOZt956K3ZZit1dygAKa7fodnCHVaL7LpDm888\/T7\/99puGTpHueg8dPLW\/9tprTBP17ds3Vu39+\/czJ89PPvlknH9j1cvvsKSGS1Y2bdrE7paA63kjdNy4jUwhSK3p1Kir19BdvnyZaTLcSwFX82LXhy4XLrWuueYa5gQGPo9Vdynzbd14Z+jQoYR8jhw5Ugc6vvN5woQJ7DI6DV0jgM6qir\/88gvTcLgIA9Bt27ZN+c4KQAeIcPsOnHTj34n2\/fENC4AOu1yMW6vUmqDxxYIPasjLi4\/Ut4lEombDtQCADpWDhnOyS9lsCxZPE7cAYSKh73tV\/4A8u+81ipF4AOHw4cM0duxYNnNF94vrA6z27sl2KVtBxzUgJhLi1qqdO3cSLvlAF5\/MHbBIa+HChXTDDTewC5dxfRH8xeHm70mTJlFGRkacZGX5Wb176dIleu+999j4lYe2bduyjxZ\/exE8u9k6COi+\/\/57dkExukh+2yKE5jV0RlOKDGSrhvSzDrjQGXCXlJQQLprp06cPuzh53bp1hBu7MUHDjZX1JfjevR48eJCNq7iwcFEKD052KatoOiN0dnY+iw3qdx1wp+7EiRPZDeLi7UanTp1iGrZJkyYMPKNmDSuEvkGHXhz2Mmi49u3bE27XM3YLTnYpm0HnJE1jowVVB6uxrln9wgocyuUbdKLRF\/5zjbdkczuaqsnEKFQz6Nw6LI78gqoD6oZJF2yQ2EwrBg2dyefFJw1dunRhXYTZXaFOdimbQeckTbE6QdYBxnV0r\/hYE3WvsHPipJuTO1j91Iyeazpc5QnQMDPDYWJ0rcaAbnbw4MGUmpqa9C5lq\/Gg053PQdcBE4lp06ax2SsM4Q8++CBhPIeJUEVFBQNOvJ3cT4CSyctz6LDuims8y8rKTMsnXsKR7C5lK+iSTZMXOAx1wEbOrVu3shnr119\/zSYNWFKEyaQ+zVx9HdMl80XodxqmBDzXdA1TbLpWTiSgoXMiPf1uUhL4H\/5wLwTNp5UTAAAAAElFTkSuQmCC","height":76,"width":126}}
%---
%[output:2b90f89b]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJ0AAABeCAYAAAA5WhbMAAAAAXNSR0IArs4c6QAADnpJREFUeF7tXVmIVjkWPr\/tvoBabVuuY1W3itPjiK2iLaJQOA8qKpa4PijqjKWIlpT7wriA4y6IigvuD1o+6IuUPjQIgrsMo+CUiFiIqCXt9lAy6rTT\/3AC5ycVk5vcLffWX7kglP9dknz5cpKTnC\/JZLPZLLjLIWARgYwjnUW0XVIMAUc6RwTrCDjSWYfcJehI5zhgHQFHOgC4desWzJw5k4G\/detWmDZtmmdFfPz4EVavXg0XL16EAQMGwLFjx6Bjx46e7zx+\/BjmzJkDL1++hHbt2sHJkydh4MCB9d7hnxk\/fjxs27YNWrVqZZ0UcSfoSJcA6bBSx40bxwjetm3bXB070sVN9xR937alw6I3bdqUWbLS0lJHuhRxwVpWkiAdFq5fv35w5MgR6NatGyurs3TWqjz5hJIiHZZ80aJFUF5eziyfI13yXLCWA9ukKy4uBnRGamtr6zkVjnTWqjz5hGyTDj3TESNGwLp16+DLly8wevRo2L17N7x69Srn4TrvNXlexJqDJEi3fv162Lx5M1RVVeWciv79+zvSxVrTKfp4EqRDz\/XJkycwb948eP36NXMqVq5cyawfzuU5S5cigsSRlaRI16xZMzhw4ADs3buXFWvIkCFQU1MDb9++daSLo6LT9M2kSIerDS9evID58+fDw4cP60HiLF2aGBJDXpIkHRbn0qVLsHTpUuZU0OVIF0NFp+mTPOm88tW1a1c4ceIEdO\/ePbf2qitHWVkZrFq1ynMOjl\/LdaTTIZon95MmHcL44MGDnFOB\/3eWLk\/IpSpGGkiHqoFDhw7Bzp07WTYd6fKcdK54dhFwoU128XapOY2E40ASCDhLlwTqjTxNR7pGToAkip+5d+9eFiMc0IPDGfLJkyfDwoULoVOnTrn8oGd1\/fp12LdvH9y9excKCgpg7ty5MGvWLGjTpk0S+XZpNmAEMn369MlimM2UKVPYQvPBgwehd+\/ebD2QiHft2jUWbDh06FCYPn06VFdXw6lTp2DkyJGwadOmenH+DRiLBpX14\/8+DLX\/qYUurbvA3B\/LGlTeMxUVFVmeOFeuXIEFCxaw+C6cK3r\/\/j0jHFo3Xkhy9epVWLJkCftt7NixDarQDT2zl55ehH\/c3Zgrxtw\/zm9QxMtUV1dnMayGrjdv3rAF6OHDh0NFRQXcuXOHdaM4cVlSUpJ77sOHD+w+Su8wLqx58+aedfn8+XM4f\/48675xGUm8cN2xrq6ORdJi6LbJZesdzIuttEzS2XJ3I1x+ejEH0Zhe42HVwPVW8DOpF90zmbq6uiwvgxNfqKyshP3797M1R+x26cJxHlo5HOOZ6D5p1v\/MmTMwbNiwr\/L16dMnFr7dpUsXaNmypS7f7L6td2ymZVIm0dKtHbIRSgr\/YgU\/o4rRPJSpqanJHj58GC5fvswqcdKkSbB48eKcQmn79u1w48YNOHr0aD3nAr+L9zBCAoXDRUVFnkk50kXbkHBM96\/X\/4SBnQaxrtWErGIFBXknEtKVlJRkUWmO4zJyENCBQCL26NGDEQsJI7Nm586dYx6taAVlGSPSoXWcOHGi1NKhRqBz587GqnYEzcY7ZOlspBWkTDbzh2nh8Md0CCTjQqaysjI7depUyGQy7P6jR49YtMOYMWNYSA5Np0RFutmzZ7Mxonh9\/vyZhW0j4Vu0aGHUoGy9g5mxlVaQdGzmD9Nq3749dOjQwaiOpKSrqqrK8vtw4ED2+PHjbPpk+fLlbL8OtIAYv4+DfP46e\/Ys3L59m8WWFRYWemYCI2RXrFjBPOFBgwZ99Symi2HaWCBT0tl6BzNrK60g6cSVP3T4ZE5faEtXVFTkduIM3Gbz+0Wcl8VQKxnxwpQ8c\/PmTUe6MAjm6bvYg+ECgWq2IUyx3dprGPTy+F3dbEOYojvShUEvj991pMvjyk1r0WIlHTkSfvtufrMXETjV7pTiO6SUCgo8ziHihVM7sgvnEXFiW7WjJb7fq1cvmDBhgqe6yy82fF5kOJGqjF\/hoXdEvYbJzqBe+OkwUt2PlXS4pX+QBAhMXL0Qt0uVTSjT87t27WLLYO\/evWPzgfi3ijReYCKh1qxZAyriUplUAheU\/W3ZsgVw3pAkhZieSNAg2OB3eFmhSFpV3jEtnKaiyXYvjE0aqg4jr\/tBy22SLzamI4BkoKs+4gWIeE\/1fRFkkwyLGlEZ6ZD0uKKCl4p0mEcMz8K9Q\/DCuUZZ+YNgg9\/RrdaI91WNEJ\/D9W+T9W3CT4eR7j5+h0iHc7biWnnoeTqedGj2yerIugW+xZqQTrRqGIvHW0UCWvxdRT4C69mzZywIYceOHcDnGd9DwuGENlqLCxcusEluWfeKlYkX5seLWHTPDzYmVhyf2bNnDwuYxZ04xZ6AMFD9HhQjEwx50mEecbjEX6FXJKh7xcqi1iQjlB\/TL44TVMCZVI4OXJF0\/POYDxnpCPgZM2awVuxFOiy3X2yCdE0qqx+mi5U1GB4fr\/tUhtOnT8PgwYPrVUNoS0eOBD\/4N+n2vBwJzCHf7alIpwPFq7s1eVdFOiQ7Bq5u2LCBxQPKtnXg0\/aLja5rlZVLhbnf3sCUVPicCenCOFGq+mNjOiIFroliV4SZwUH+\/fv3lQN1rxZIA1TyvPyQjh+PYaZVhQ5DOqxgjHymoYSXpQuCjYp0hAtfGdQ4\/ZAuKowSJZ04FiJXXiyczHrJvFdqnT179mQkxqhhPLiDxngEuti9ipbJazwThnSYzqhRo3IDZJ2zwI8TTbAx6V7F\/Jt2r1FilDjpdF2C6F57WTqxEslyejkSNFeGYfK8s6GaRwpKOrFr5bsZmfdq6ony0zcmY1Ux\/yaOBHbz6GVHhVHipPMK1pRVDFkvL0tHc3CqwomtWzev5Ge8QtZbdCQwTQzJ4j1aE0unCmRVkVbXiEVMVGM3ccokSowSJR0f1YtWRtY9iJZNNzlM0xbUHYlpqKwB\/Y7jSdEhiYJ0\/FQJfU\/nveK5YTQ+NcGGvkvDE3FVgV914MesRCj6TYVxVBglRjrVMpi4HINA8uDZWAYjUMTuROd58ZXOWzr8Hq1C8EtQOu9VdGZ02PCNQyQI3VM5SH6XwcJgZEK62CaHvaYmbN2TjbVoPPX06dNAS2W28m4rHZsYEflxcvhJp0fw5vOv8G2L76C069Tw4eppOk5ddBrCzFHZIoLtdGxhRKSbsXkK\/NKhK\/zepBCa\/P4K\/tq9KVOfhRLmpIl0NPindVOxS7ddwWlNT5zKChuJIisnka7vmsVQ+8Oc3CMTv3sHm34sDgWNC+IMBV\/+vkyk++HvR+HXP\/yUK+iEwuawuV\/rUAV3pAsFX\/6+rLJ0kwoewoY\/\/xyq4I50oeDL35eJdNOXFcPlvvsh+00hZP73Cr7\/7wq2q0CY3aIc6fKXN6FKRqQbXlEGD7ityL757T4jHzoVC3q1DLRblCNdqKrJ35eJdK1XHoGmfeuHNlGpgzoVjnT5y5tQJTMhXVCnwpEuVNXk78uxks6pwaJXg+kCXFVqOb9LYDrK65RgXmo5E9Lh1AlaO78Xs3Qm8V\/ih3UL\/mJURmNSg3lhI8Ya4ubieIkRN2HC1PF7umgUqnOVcInuF8\/\/G7z\/aXiu+v\/U5Df298\/fj4AFRWZ77onccWqwGNRgOsKIjdymEgwJYKKWozweX18GNUPbsw0YB7QthKmfH0Ozb4dCq75L\/Rq43PP1SOdH8eTUYHOY6IcuU6UcWTU+XMokgFO2Za5Y6xQ1EoVaLkjvZ8rCXPfqV\/Gk614xA6KcUReubpppfC5o5DD\/blxqMJ2lE+PmTEPVbeHDN4xYhDlODRa9GsykFyARFI7pVKQLE2UTplHypCsvL2fnh\/CXarNE00bh1GCc1UTQRGF2EDWYznsVdyXwQ7qolGA0tlOJ0WWBqkQqJCL+C3rl5un8Kp78eGiNTQ0mw4asFlaUuEWEafcapRLMlHS4EyfuQBC5pSMX22undKcGU+8ib4qNingmjkTUSjBT0sUypqMgTqcGk2+gQ5UjzjvyLV8U9ngp5ajb4rvYJJRgiZOOgPCjeNJ5r41ZDabzXmUqMdtKsERJ59Rg9feSkw2O\/arBdKTj1Wf8t\/0ug4VRgiVKurRoJGwqnYJ6XUm+Zxuf2CeHkwSTT9uW0ikt5fWbD5v4NBrSkcl3ajA1HW0owWJfkUhL9+q31bvn40WgUVm6eKF0XzdFwJHOFCn3XGQIONJFBqX7kCkCjnSmSLnnIkPAkS4yKN2HTBFwpDNFyj0XGQKxki6oGowvnSx+rCGcf4VloKUkjCIuKChgG3LzYehUTpWCy6uWvWLS8D2dKIa+HWZXJgomwChuVci7TDUWK+mCqsH4CsPF\/YZ4\/hWWgV9ewuPckXSm553pzIpXxVFDReLxZ6NFqQpTrfHy+fY6pwx1HLGFNuk2elaBq9vMOe3nX2G5sJLpTAmTMHMZIVX46KxFnPiIvY9IHnHLW1k0sxXS+VGDmWxbn\/bzr5As\/JkSJqTjxUVi9ykOKUxIh0cBiBtbiwImr\/M0ZITnrWhpaelXZ3iYqMZ0eddZea\/7gdVgQTJlGpbtp0A6AYrXfWwUy5Ytg7Vr1wJufK2LEcR8UVcoK7+4\/b4OI1EiEAc+OsJ6HbcQm6ULqgbTda0y4vgRoJgSLwzpsELo+E1UZfkR1JgcielFOlkEcRz4pJJ0Qc8Gy4fzr8QzJbwsneq8M\/R0Vd6lznuVnS\/BHzJMDU8Wzm6qCkst6WhsI4aYe50Npus6eO+Wxoqm3UeUqieVJZSdKeFH4YaWUXb+BE8kGUaipJH2MiGnRkY6MV9+8Ek16XTdpehamzgSaT7\/SuxasdK9SGfi4YvaB1XDVBEvDlVYqknnVw2GLVRH1DSff4V5Fw9FMbF0dN6ZbMwpNkSv3kDsrmnOEI88FQ\/uC3M+WGpJR+D4UYMR6A31\/Cvx+E2dpRM9TVmDE0nmRTrVxG3UqrA0ku7\/1ZeKMcGxEl8AAAAASUVORK5CYII=","height":76,"width":126}}
%---
%[output:43492d8f]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJ0AAABeCAYAAAA5WhbMAAAAAXNSR0IArs4c6QAAD5FJREFUeF7tXWlsVdcRHjdgthgbkMVmFrsikRtSSm2zNQmCmJJAHFRcAbEUVpWllCXsq9jCvlRsEvv2g02CPwj6oxI0FVuEKkTLooJACAFGUEAVFKhK42ouPZfzxme7y7v3+fncP+B3zzrz3Zk558zMyaiqqqoC+1gKREiBDAu6CKltu3IoYEFngRA5BSzoIie57dCCzmIgcgpY0BGSnz9\/HioqKpxfly9fDoMHD1Yy5eXLlzBr1iw4duwYdOrUCXbu3AlNmzZV1rlx4waMGDEC7t+\/D1lZWbBnzx7o3LlzQh2+TFlZGaxYsQIaNGgQOUCS0aEFXcygw+779+\/vAPzdd991R2NBlwy4p2ibUUs6JEOdOnUcSTZw4EALuhTFRVKHFQfocEKFhYWwbds2aN26tTM\/K+mSyubUajwu0CEVxo8fD5MmTXIknwVdauEiqaOJGnQFBQWAi5HKysqERYUFXVLZnFqNRw06XJl+9NFHMHfuXHj9+jWUlpbC2rVr4cGDB+4K165eUwsjoY8mDtDNmzcPFi9eDMePH3cXFR9++KEFXejcTdEG4wAdrlxv3rwJo0aNgkePHjmLihkzZjjSD\/fyrKRLUbCENay4QFe3bl3YvHkzrF+\/3plKSUkJ3Lp1Cx4\/fmxBFxZzU7WduECHpw337t2D0aNHw7Vr1xLIYyVdqqIlpHHFCTqcwokTJ2Dy5MnOooI9FnQhMTdVm+FBpxpjq1atYPfu3ZCXl+eevermNGbMGJg5c6ZyD44\/y7Wg01E0Td7HDTok4+XLl91FBf5tJV2agEs2jVQAHUYQbNmyBVavXu0M04IuzUFnp5d8CljXpuTT2PZAKGBBZyEROQUs6CInue3Qgs5iIHIKaEGHK6l9+\/bB9u3bnX2pDh06uIPEd2fOnIGNGzfChQsXoFmzZjBy5EgYOnQoNGrUKPLJ2A5rBgW0oLt48SIMHz7c8fWioDt9+rTjeNi1a1cYMmQIXL16Ffbu3QuffPIJLFq0KMHnv2aQw44yCgooQff06VPHkxWBl52dnQA6fIeAQ+nGB5V89913MHHiROe3fv36RTEH20cNo4AUdHj2h14P33\/\/PXTv3h0OHTqUADrcREU1ipuYvXv3dqf9\/PlzmDJlihOGhz5imZmZSpLcvXsXjhw5AuXl5c6REv\/gGJ49e+ZIWXThNnmiqJMufSA9\/czFhA+qMlLQoerEeE60165fv+78y6vXgwcPwqZNm4R2Hko5tPFMYkDZCcD+\/fuhW7duCWN99eqV48bdsmVLqF+\/vtFco6iTLn0gQf3MxYgRikJC0D18+NBRkT179oSxY8fC4cOHq4Fu5cqVcPbsWdixYwfk5uYmdIHv0FsCg4jz8\/OVY4wTdLuubIXKF5WQm5kLfRr3Mwa3H0Z5reO1vF8A+ekndNChuF23bp1z6IwOhU2aNHFUK5V0CCwEjEiaicrLBhoX6E7cPgbLLix0h\/WrVoNgfNEkI4nqh1Fe63gtX6NBd\/LkSZg\/fz5s3boVOnbs6DAlCtChSh4wYEA19YoBKs2bNzdOqYDMMqmz8uK38Me7f3D765HzMcwuXmDUj2kf\/GS81vFanoHOZO5BxoV10b42tbFFwibj3LlzblJEHDBGImFEUt++fd3yp06dgqNHjzo2HgsGPnDggLNFgr78aOjzD77DBQiWb9GihVIao7fs9OnTnf74PpmRi6vkxo0bQ7169YykOkpqkzrXnl2Bo\/cPu20WNegCvZqXGvVj2gc\/YK91vJZPBr3wYxfxLycnx9GAfp+M\/Px8m4nTL\/XSvB7uvy5ZsqTarkKoki7NaWin54ECqKnQphftKnhoRlhUeyIRtANbv2ZSQLXACzojC7qgFEzT+hZ0acrYVJ5W5KDDPTjcMhE9ouyUNILJNCOljOiYPAZXxatWrUrwamHlnzx54gSuYGQVPcXAMlh\/2bJlzkocvWNM54JbQ3jSYnKSYgoYPhEOq8MiyXiPHfaOxmiYZAOVjUVHJ6yHvMYHack\/sYBOtPHLCDhhwgQ3LSoDHBKSDVy1caxjFiMUplegXi1Ylwe4zMjFsaPjAY5HNhbRXMIEnWqc2M\/s2bOBhSTygJs2bZo7b9EYdfRj703oJBsHtpEyoGNfBg9I0cYxAw66O+ly9oq+LvxNJA2o1JCBDoGGR3goBVUfAH0XJuh0pzL0PaMZjpmXOn7GpKMT1UwU\/Dzovl5eAYsGf2uKdaNyzj4dDXHTMQqTurDEy1iW\/5v1KvtdNir2ZaE6adeuHfBfPFOZmBwax4ppUvH\/a9asqaZekXnoy7dgwQLH08XLXBiDP\/30U+cokD0U3CKVyZeRAYifO5bBPsaNG+dsuLM26Zxkv6tMExWdGODu3LnjOGygCcNrKarmm43Jgt98NhqGvj\/K7TLwPt2LFy+q8OTARD1SW0ukWnnQyc5mdZ8D1qOg4+uoGIF18USE\/yhkpgK1G5m64W1SqmZEKo+O149qks05iIrVAVbFPzYHBF1pSSmMzv+dy4LAJxKia5pUCwle7akG7UctyGwbClIVMbFffJhaN50L1pGpRF5qX7p0SflBqNpRfWwy0Pk1VXjtINIIvH2sk3TzyxfCL\/M+D0\/SyUAnkg6MANg7rvAw0xCVkmxkFHRUJami1v1KOvwIli5dCsOGDXNXvTL1SueCqlj2ofDjQU9pXDkj+ES2kAp0TJLy4GNtmIKOjRv7x0e1UxCGpJu+fiqMKxuvU06e3js2HTXadW5LuPJCGwYnLAMd3wbmWKM2mMrm8ws6JDLGaGAyQXbRh+lc0IBHUKCPIL0oRDQeKkF5AJqoV6olTNRrnz59HMDzCzSVRgkDdN\/MK4cJI9+ktwjryXj8+HEVywDJtihUjKIENVlI4FeJ5fj9L9VenF\/QUdWKRPIyFxNJJ9pbo1sPJgsJCjqThQRKWWqHqvbiwgBdWf9CWPhxF2gyaEFYmHtzC6IISLJFALV7RIyidohu743Oxg\/oGBO\/+uqrhBWtTtLxzqkmNp3oqiTWN86DSUndlgkFncx24+nLzBm8EsrkID4M0H3xZS\/47MZ78MOPKyC7TSZ8PL15YPA5oKOM0W2ooj3G9pIYsdq2batdMfI2jW5X3uvqlW6VMMp4mQsbH29vst8Yk0WqU7bCZCqYnirwpw48eGhfsnbpqYUMgGGAblBeOTT++9dw5dfZ8K\/cOtC+23P4aZs\/wcgPxvgGn2PTUWM0imMwJBxVuUFWr\/wpBE8NL3MR7dOJPg7KdOxPdlxFDX82NtVpCrsQT9Uua0dlpoQBuvd6F0Oj5sPgZukvXLK+859L8EXrNrDogwJfwIvEy0RkayEzpk6dCnPmzBGer\/qaTZpXEn1YMrMiKCnYh5X1zXj4UWF\/qHqnugf42Pb1YWy+WZQeP55IQCf6GoPs4wUlaE2tL1o06Oxfv3NloGs4YxvUeb9Y2MyXLTJhcWFDz11EAjocFd2nC+qJ4nmmaVKBqmuVbRxkymkBuiAEsHWjp4AFXfQ0r\/U9moAOVSuqWK9PZOrV68Bs+XgpoAOd30UEzsqCLl7epmzvFnQpy5r0HZgOdH5XrlbSpS9mAs+MB11Odkt43uLNNe\/s8WvPWdAFZk36NsBAV1KwGjI6\/Rlu9S6A\/9b9mTPh4pw6sLXL2xMKr1QQ2nRejo6wQxsNlkh2kTs7X0K2R0mP1+KKBMOxsrGM+\/nv4YeSf8KxgjdXguIzp2Qh9Gtf5hVrbnkp6EReJjYazIzOKhdzkYMEY3KqRILxoGNnxJjL7+Kjv0Dn3KJAh\/1S9aoLZuEBaaPBqgNRF9dAPVVSLRKMB92uXbuqBT8FDsyx0WBvI9vCigYzBR1TnyYOnKKgcgp31o4sYs40EowHHUasoTmAD6ZW+8e\/H0Kb7Lbw26JJZmJfUCrDRoO9zSIQVjSYDnQiHz2R\/6CuHRXXg7g18aDDO0SKi4sBM5eu+9sKt8uRPxntW816XkjYaDB9NJgKLOxdUVGR6\/Qq8xSJKxKMBx1e6YB56vZV7oC\/vr7ogu7z9mUwt+Rt+lwvYs\/TQsJGg+127s3QRYPpVq+6dBKMgSLQ0Z2FZHgN86BjY2lYlAk5g97eghRkBWujwf4f2RZmNJhI0olCHhlDTSLBMI6XBkGpVGhY6hUvOmYpf69mXQplBWujwSoq3CCXsKLBZOpVBjyThQQL9+zRo0dCfhgEIj4061JYoDMJAPKiWt0tExsN9uaivbCiwfAWIIzz5bNb8VINYyB4FWsSCcaCwUXZnkRMT3nQ2WiwN7nZwooG0606RVFippFg1HNYlmUgpUFno8EGu4IirGgwHehkueO8HoOxdqjKxQmlNOhEuUy86mhdeRsNpqOQ\/r0srhdpe\/v27Wo2nb5FdQlRfG\/QNln9SJw4bTRYOOyii4Yg+3i6EdV40PHiHhMo4mOjwXRsF7+n+3RBPFFUI0gL0Pkjsa0VFwUs6OKifC3u14KuFjM\/rqlb0MVF+VrcrwVdLWZ+XFO3oIuL8rW4Xwu6Wsz8uKYeK+hEkWGy876oo8JUCQGRWUg4dqfEhg0bjO8IUyXh9gICeqxF68oyzPM0D5qVSZWTmI1H5KkSC+hYp5QwsvzBUd8RZpLHmD9+06WB5T1CwgadyD2Ij2fg3ZJo37rcxaqPIMj9YJGDTpcdXPQ+yqgwXoLIJAG9U8JLhFsUoEOwUJqJGK26IEYFOOq9TIGvux8sqaATRYOZfF1YJisrC3r16uXc2WCS2l+UmZwSjk1WdkeY7j1rD4nO3ymhA53ovjOcHzILH5Pcw7SMjnHUnUnmROo1a2kYUWG6sXsxM2jZatFgovT0ug7iuiNM5ubNpAjvfaFSr\/RuBmZT0SzrfPp\/EVMoOHSMw34wPT9\/f4focj8TISDjkV8XJ93YdZhQva\/mZaJTraLG4rojTAU6ZCi7fhPHrEqVYXJjEJ2jifRRMY694xdlMm2hmqeO+TUadKLU9IxZeXl5sdwRJmOGyPdMJulEcQsy5vNAY1dPoWSSeXroVq+0ninoorgfLKmSjt4NZqpeeYarQMczO+w7wmSgw9\/Zzdb8tgD+LroqndpWsmAXKt2oMY598UASMU4U98qPUadeWQhksu8HSyroRHeDmdgQlOEmC4mw7wiTgY6qVqZeZaCjBDaRdBgoQx8a+yBjnAx4JgsJ\/HCjuB8sqaAT3Q1mYtdRhosIluw7wkSgk12Kolq90o9MVNZk64LSTcU4Jl2pZKTpJWi\/Jntv\/MeQsjadiMhssLm5uQkqiZ80v7qL444wEejoVolOvYo2aZnEokY+v9IUaQMKMhXoRODhE9wwM0Cmdai9GHakfzIl3f8AKqvzMQjK6PoAAAAASUVORK5CYII=","height":76,"width":126}}
%---
%[output:2a1c4d19]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJ0AAABeCAYAAAA5WhbMAAAAAXNSR0IArs4c6QAAEq1JREFUeF7tXWmsVUUSrjfKKsqiJmyyzYghSBCFoPCEUckQEBg1A4qiKEQYgz6IyiKLuyMIKKgvIKK4PyWoGQnG5ZcEVELcghFRcIQBQWTxB\/syb\/J1Use+\/Xo755670ich4d3TS3X1d6q6qqu6K2pra2spPIEDeeRARQBdHrkduhIcCKALQMg7BwLo8s7y0GEAXcBA3jkQQKew\/PPPP6cbb7xR\/Pr444\/T9ddfb52Uw4cP07Rp02jlypXUvXt3euGFF6hFixbWOj\/++CPddttt9Msvv9CZZ55JL730EvXo0SOjjlxm6NChNHv2bGrUqFHeAZKLDgPoCgw6dH\/11VcLgDdp0iSiJoAuF3Av0jbzLenAhtNPP11Isuuuuy6ArkhxkVOyCgE6DKhLly60ZMkSatOmjRhfkHQ5nebiarxQoAMXJkyYQBMnThSSL4CuuHCRU2ryDbpOnToRjJGdO3dmGBUBdDmd5uJqPN+gg2VaWVlJM2bMoBMnTtCAAQNo\/vz5tGvXrsjCDdZrcWEkdWoKAbqZM2fSww8\/TKtWrYqMim7dugXQpT67RdpgIUAHy3XLli00duxY+u2334RRMWXKFCH94MsLkq5IwZIWWYUCXb169ai6upoWLlwohtKrVy\/66aefaO\/evQF0aU1usbZTKNBht2HHjh00btw42rhxYwZ7gqQrVrSkRFchQYchvP\/++zRp0iRhVPATQJfS5BZrMzLobDS2bt2ali1bRm3bto32Xl1jGj9+PE2dOtXqg5P3cgPoXBwtk\/eFBh3Y+O2330ZGBf4Okq5MwGUaRjGADhkEixcvprlz5woyA+jKHHRheLnnQAhtyj2PQw8KB7Sgw2K2pqZGBCRiT7Br167ClB80aJDwmPMDNbB27Vp65plnaP369XT22WfTmDFj6JZbbqEzzjgjMDtwQMuBOqAD4GbNmiUiYQG0nj170po1a+jll18WURD4x8DD7\/i7d+\/edMMNN9B3330nyvXr148eeuihjKDEwP\/AAeZAHdB98skndPvtt9MjjzxCI0aMoIqKCuEzgqf83XffFW6C888\/n\/bv3y8AB+kmR72iflVVlfht8ODBgdOBA3U4UAd0UKkfffQRLViwgFq1ahVVYKvujTfeoEsvvZTwN9QorKwrr7wyKnfgwAG6++67RZ4ANrHr168f2B44kMEBb0PirbfeokcffZReeeUVkUTy5ptv0rPPPhtJPnmdBymHNZ5PkkqYj+LjwOFNC+h\/h3bQnxq3oUYXTEqdQCfojh07RqtXr6YHH3xQSDioXewTzpkzhz799FNaunQpnXvuuRmE4R22c5Dl1LFjx9SJDg3mjgNH\/7uCDn41Oeqg0QUTUweeFXRYw91zzz2CgD59+tBTTz0VAQzAgorVSTNIRVi0vP6zsWj79u309ttv0yWXXCKMFvnBWhLqGllSstXsYnnSemg3ad2k9bLpM5u6JnqPbphGJ3f+O2Lxaa3+Tg26zc5gOeYiznyo82UFHaIddu\/eTZs3b6bXXnuNEH4Dldq5c2ch6dIAHa8VR48eLdaI8nP06FERXwZJ2qBBAxfWovdJ66GBpHWT1sumz2zqmuitv\/9DarR9TsTLw22n0rHmAzN436xZM2revLn3fMQCnVx406ZNYj8QEg8qFoaGCXRx1CuDDmtFVdLBffPrr79Sy5YtqWHDht6DTFoPHSStm7ReNn1mU9dGb+3WxXR8zzqqd05vqmj\/zzp8z4mkO3jwoDAYYKkixgtuEez\/7du3j7Zu3SpUKixcSL277rqL3nnnncg5jMx1SKevv\/7ay5BQrWJ5hEeOHBHOaVjRcUCXtB76Tlo3ab1s+symbjb0en\/9hoIZ6hUi94EHHqAPPvhA6OybbrqJLr74YtqwYYMwCrADgVAe\/B+S7+abbxYqF4kl7BzGOg6q8LLLLhO+OpfLJIAuOdBLFnSfffZZxqGI8LvBWh0+fDhde+21EVY\/\/vjjyDC4\/\/77hTRDQknjxo3F7gPreBgfy5cvF0cljBo1yvlRQJJOnjxZ5HtiZ0N+sNjds2cPYQ0RR9IlrccL81LpM9f0QsDgX9pPRceOHcNJnGlztUzagxBAeFXawKtQJV2Z8CsMI0sOrFu3Tmx98g5Uls1lVHc6h9PsLLRVOhywrbWzHUUAXbYcLNP6DLoXX3xR7ETJT05cJmXKxzCsGBxg0D355JPisEf5yZtzOAa9zqKc8TRy5EjhA+RTKdWKvidb2jrU5Tz4tIutPOwtZ3MCJvyacKgjA0yVFipdPjSp4wSNCLzIRWCFzWlfkpIOkwE3C3yCyGAH6OBkVo9atW21uZDNEw7XjroHjHafe+454\/GuzPBsEmLkVEJ1MQ6w3HfffRmLdNCEwFmf\/Woeez5AVzaGBCYVwZ5yDqgOdHxclu6dC3QuwJreMyDRflLQycd8oR154hiMyJvF+PnhjwROdtc5xyUPukOHDtXyQc08GE4KlidWngz8rqoDvMdhLxdddJHYm8XDCcmINFbb6t+\/v1A5NmDxu3nz5mWoJxctPmBFGWzfIUqaD6aWpQ3eYTxx1Sv3DcDiOFdIcZV+3QfDoANPZDCyVOQ68uHbOkmnllfngOtg3HfeeWdEiirRcmq9VlVVCecwM1f3xTGguAx\/rXI9BoIsHXQhTmgf4VLTp08XYe82gKBNPPIk+NASJ7TKJDHVflySVffe9NHYyspSXVW5Kq9U0OlAqEp0BqUsNHQAyynorrnmmlqbSDcxDkTde++90RpEp650qgTtIXkHx2AhGFRVReqEyFI3G1rigiafoGM+bdu2LTIKTGOVDZz33nsvMiQwPhgt6lyq7Zg+SHW8kctk5nj665g\/lgFx+agrXzF79uxaLKp9rCf1nA1ZdJsmSfelgRBet9gkHX+VuvscXLTowq5UtQw6THdFuEDneo+2fSUd0yWrOPWj1k2ezZBQx8ptm+qoYGTQLRk3kAZMW5QG1qI2hHPYZb7LA2Bw\/vDDD3UknaoK8bc8SEi2xx57jBCwyes8G+hY1bdr1y5S\/z60+KhX18LdBirQjLAuPFgXqWtW5q4P6HSAY765oq916hVWsbyexv\/ldaXJFaSCnDExaFQVzZ+eGVOXusuEGYVoXfh\/AC7cIKNbaKrqVbfwliUd3COyapWlgc5CVdeO33zzjRctPoZENqDDxPHz888\/Z6w5ZZFgA53NpZIEdCbXk0696nx7Jkl3oHIKjRncm8b3bhYNLSfOYZmAL7\/8UpvvoC5yfdZ0aFedJB9JxxadbT2i+rhcLpOkoANYWFrDsc3+Rt3VTCbQMeC++OILo1\/OZ\/2KuWEAmYSD6hOMu6brOmQYte41lhYO\/8MDkbWkgyEhm+km6SKvfVj0qms6rA3lhb8OmOwqUVWQyTksg4n7ddGCtn2csyaXDuqb1KtqCKGcOiaXelX5YlowqeVUt4qsXlnSwXvA1r5spMlrOqhg1cugOquZ13\/7Rx\/qP3webdt3mLbtO0LtWjSkqQOzy\/Cr2Lt3by2sHqguflSnqLrmw3ucWYIMfwYLTxIuWMMA8cjGibwLIUsFm\/WqM258aNGpOKh+fmxg4zIm0Kl+MJTX+TXlpYPsp3NZ6yrv4\/jp1LbBP\/SNZZCsLSAdr7rqKsK+qrz+k9emMuhO6zGDatbvivg3dWCHrICXWpSJjzVn+qJL5XfdXqrqdyz2sfhuncmgW35ibMawIO1G9mqZGHgBdDFQgonAxyVvsLMaR5ac7\/ZVjC5TL5oG6JiopBIvgC7GtOp2SNjSzFW0RwzyvIqmCTpIu+qRXbz6lQulBrrYPYcKRc0Bm3plwgPoinoKS4+4ALrSm7OSpziAruSnsPQGwKBrOmAMbW3SVzsArOegYuM+YU0Xl2OnSHkGXffR\/6IT51xAa7f8njHyvn9uRisn9EjEjQC6RGwr\/0pqPN3Q6q8ygJfUiADnAujKHz+JRqiCDjsSE2r+uCgvqWotGOjymQ3mw3FsH+F+1SeeeMIYpmRqh3cpeBvRFZeYTV8+Y+Eytu0zXTtqyJMucnjOh\/+hNZt\/p8q\/NEu8G1Ew0OUjG8x3gmxZY642dHkNtuiWbPpy0SK\/VyNJ5LwNOfSf6+iy33Iaro4gzjgDSqNsPrLBfOiUgwd8ggDUNnVhQqZQrWz78hkPyphCtkw7EXJQrBxskFPQlWM2mMx8k9qTw6Tat2+fEQXtO8G+oPPtS1WJug9BdzWnT26qKWmHQ8fU7LcoRyIXx0qUYzaYTp3Y1J5PPoIOiCb1akuaNvVlAgXKc4CBLnzfRyLp4hDV8ZgSc3JyrES5ZoOpofOuIx7k0HtfScflTCpK144OdCaVGCeTSwYn9yvH17kSx02gQ8wkTr6Xnw4dOhD+JX3KLhuMb5pWM+htDHJJOtN7nQRRw8PVfl19obwukwtWMZLi8agJ4L5t2iSwCXQ6vuHUVPxL+pRdNhjyFhAJrWbKJwUd50Sgfrdu3aKYOV3CuSz5TKcD2NSrLZOLQQfg6B6XIWQ6QcBEM39QOImzTZs2GV1meyxsHedwqWeDpQ06zonAGcy4xmD+\/PniGArTmSSYHVvMmg50JotXVq8MujgSXEZKUtD5GClxJZ52R6LUs8F0ofM2kNjUE6cbIipYTsJJU9KZjAFdJpcpfZB\/5wQd9fwUV1qmSb3mBHTlmA2ms15tCdgm0Mnphkhakf2L+LrTWtPp6NVlcumsVx2YVEud64Fm01l2eQVduWaDxdmeMoFOl7mlbnOpZVxrK9++dJlcAI2vn071+Zky1lxrupxIurR2JMoxG0yXE4Hf4C4ohSScuGstubyP\/y9p+6lFmZQb6Ex+PdNZIEknoFjrBdAVYGbAdDXdEGT4HIpTAHJT77IkQJf6qEODBeVAAF1B2X9qds6gwylbl19+ecQEXFoIBzXul\/j+++8Jx4jgLBdsI+JIN58ntTWdT2ehTOlwgEE3ZMgQcT8YbraE5Txr1iwBOhwfBzDiatSlS5fS8ePHCZcZXnjhhc5BBtA5WXRqFmDQ4fbJV199VWz642izW2+9Vdx+OWLECKqoqBDMwU2W48aNo65du4pDzl23kAfQnZqYco6aQYdtRahOuIhWrFhBS5Ysoerq6oywfsQB4\/i29evXe12kEkDnZP+pWUCWdLhkxuaXPHbsGOEOYFw87XN7TwDdqYkp56gZdDjq9fXXX6cuXcwH5QBsiOwZNGiQOJARJ3Xangh0PjFZTkodBdBHTU2NiAd7+umnxVVJusd04nm2\/av1Tacw+fQTJ9vKFkDq05dvmThbf2jTFdgKYwEXsGC+TEDavXs3VVVViZvOFy1aRJ07d3aSm1fQqREbtmjXJFczOUcrFWDQuPYkdW2asq1sh3XD4svFPibTFzd03nXQNks6qMsrrrhCy9r9+\/cLtbp69WqxzqusrPSagryBTo3YsOUsuA6p9hqZoZC6YR4XdKa4NB3NajBALkGni+FzhcHzkbg6ulzOYVisCDpF4hO0Vr9+\/SJr1jU\/VtDpoixk1cfvceXS888\/H51brJtI3QHROkkHgnX7uC5aeKA2tSffTIN7IJBcHTco0hd0criS7Y4wVSViHDr+MQh4nK6cB1l9yrfo+NJlA93GjRtFuPrJkyfFmcY9esQ700TkSGDxp67pdHuMKiG2mC\/dxWrycf4maabLgPehhcEq5wHYAhdtQZ2uLzVuMrNpv1YniXQ0q8GctgBSpl139ZM6Lts+sgl0uDrgjjvuoLPOOktEUXfq1MnFrjrvKxDEyZeUuC4jUSfKllisZlepR9+ryScyZWpMmk8ksImBpqiQbEAHWuOoThttukhgebygU5fzYfugZN7aVHpc0LHRAAmH1MTzzjsvNuBQwWtNp4p2Fv8molWG6I7zN0k6V5SriZa41rcLdKb3Ogmii+iVZ8MnMkUFMatP0y1BLvrlD0O+W8KXLuYzlgbw0fXs2VNsc2FLDNtfuGJVfZo2bSp2Kpo0aWIFoxF0KhOwlhs2bJhIg+N1kE4V6tYSapg3q0LTmk5VJz602MLRdRxwTRq7d1AXzlG++8KkdmzAcqlXPoWAP2ZZ0jHoTLPoMoRsfPGRdOgXc489WMz9qlWrjIByHR7EFbWg49xRFJJzLE3q1ZUEoqpWF+jkiTXle6q0pA06phk3cMuRwiaJarIUZYmjuyxZ9+GpoMsmETxb0OXC4q7o27dvLZiBhwdnSuOTLR\/bVejyxKAt+VJhRrvNLSIzKi4t6sSaQGKTdOqpUvIlemlJOpMxoKpqXtOp97iqIDdFbtv47CPpcgI6nGUCaQYxzqBjSSdfeiv7t9Q1HYC0bNmyjJuqeR2hukpcoFOBrVtD6WjRWa+2XE9XSiIkHD4s1b+Y5poOgNBdpIedGtfdXSqYbFavaYenUKD7P8rLWUBq5Iv5AAAAAElFTkSuQmCC","height":76,"width":126}}
%---
%[output:065be3b7]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJ0AAABeCAYAAAA5WhbMAAAAAXNSR0IArs4c6QAADeJJREFUeF7tnQeIlcsVx8+KrhrX2J6ubdEVDLFFjIoilmfBGjvWoMaCBSvWtTyNxt4eVuwV7ArGqBAiYmL3QfDZcEWf+NTVrOUZFdfVuOE\/Mte543zlznfLd+\/OB8Lbd2fmmznz+87MOTNnJikvLy+PzGMkEEUJJBnooiht8yomAQOdASHqEjDQRV3k5oUGOsNA1CVgoPMg8osXL1L\/\/v1ZCYsWLaI+ffrYlvbu3TvKyMigY8eOUd26dWnr1q1UunRp2zx37tyhwYMH0+PHj6l48eK0Y8cOqlevXlAeMU3nzp1p8eLFVLRoUQ8ti2xWA50H+UYbOlS1U6dODPCUlJRAzQ10Hjox3rLGArqCBQsyTdajRw8DXbwBE476xgI61LtGjRq0adMmqlSpEmuG0XTh6M04KSNW0EE8o0ePpvHjxxM0n4EuToAJRzWjDV21atUIxkhWVlaQUWGgC0dvxkkZ0YYOlmnTpk1p5syZ9PHjR2rTpg2tWLGCnjx5ErBwjfUaJ\/DoVjMW0M2aNYvmzZtHx48fZ0MrjIo6deoY6HQ7Md7yxQI6QHb37l0aOnQoZWdnM6Ni6tSpTPvBl2c0XbxRFGJ9YwVdoUKFaN26dbRq1SpW44YNG9K9e\/fo+fPnBroQ+zDukscKOqw2PHr0iIYPH063bt0KkpvRdHGHUWgVjiV0qOmJEydowoQJzKjgj4EutD6Mu9QidHaVr1ixIm3fvp0qV64cWHt1auyIESNo2rRptj44cS3XQOck0QT5PdbQQYzXr18PGBX422i6BIHLqhl+gA7RBhs2bKBly5axahroEhw60zw9CZitTXpyM7k8SMBA50F4JqueBLShU1lOvAryrliY9CdPnqT169fT7du3KT09nUaNGsXmH4ULFw7UHPOTc+fO0Zo1a+jKlStUpkwZGjJkCA0cOJCKFSum10KTy3cS0IYOSzDDhg2jcuXKUc2aNYMaVqJECerduzfb3QqQDh06RLNnz2YbD1u3bk2XL19mLoSRI0fS2LFj2RoinrNnz7ItO40aNaK+ffvSzZs3aefOndS8eXOaO3du0G5Z30nSVMi1BLShy8zMZIvMM2bMYFuorZ779+8zk75du3Y0ceJEBhhAPHDgAC1cuDCw5\/\/ly5cMOGg3cTv2mTNnaNy4cez\/dezY0XXDTEL\/SkAbOmglDJG7du36KlBEbO7+\/ftp+fLlTLPVrl078BP2hAFGgDRmzBiC+wHDKMz\/Vq1aBdK9efOGwYoAFuyuSE5O9q80o1izfy17Sq9+zqUSacnUbEpqFN\/s\/VXa0O3bt4+2bNnCYEpLS1PWBBoNGgrzMzny6e3bt8w7X6BAAbY95+jRo7R27VpWXvXq1YPmeVZleG9+fJbw4\/6X9LdxPwcq32xyalyBpwUdh+natWvUq1cvBsqNGzeoVq1abBG6Q4cObBjlxsanT58YWKIxwH\/Ddhzs98e\/8+fPM5DLli0bRMOSJUvYOiPC72CE2D0PHz6kw4cPU8+ePdmyk\/jAoHn9+jXbdcvnkaqy\/J4OwAE8\/vyuTyn6w2r1hy+25cL3zy21o9s2h+Mz1YKOaylsJISlOmjQICpSpAgdPHiQTp8+zYwD\/Pvw4QPTZnhUsZiACcMqtODmzZsD\/y3HgmKIhkUra0GVAPgqwZ49e6hx48ZBSXJycthW7woVKrD6Wj1+TydrOgAH8FQPb0v2mcL0j2nPLLWj2zbHDDpM+idPnkwlS5YMsirxtaxcuZJ2797NtBKs2lhBhyG5a9euX0GHrd2pqam2wcjoAD+nQ6NOL35E9\/75in7Tsozt0Mrb8u+lefTT8Q8BedTsWZzar6wQ+Nttm5EBo4TdSOEEppamsysUC9CwaqH9YChYQacaXrnWkzVdKMMr13R4PwwT8Xn\/\/j3bbYvhW\/QPyu3xezrUN9Q6vrlckgAefxp\/V4TSOxUK\/O22PGSAsilVSq1ZnYDD79rQ5ebmsvJla5JHJsHxi23U0TYkOHSwqhs0aBAkA4D+9OlTKl++vO3w6vd0aNTLA3+mX+7foJS031LZfn+x7GuxLef\/uJEe3ylLFatnU6vD47VkExZNd+HChZAPRcSqwvz585lGE90bqBD8d5i\/DRgwgFq2bMnmeJiTYW+YaAS8ePGCli5dyhzB3bt3Z45gAAr3iHhWB4SGlQxM\/rE64aTWsaN2ypQpLCYUZYsPhv9nz56xL9VuTueHdP+7nE259\/9LyVV\/TcWaBxsJOTfPMOj4U\/zbP1HKt4OU4PG2FHl0ld4dW2yZR9VmGGKyMeZGkzmlSUpPTw8ZOqdCze+JIQF8tNgyFW7wknQ0XWKI1LTCTgKXLl1igT8qL4BXyWnP6by+2OT3twTsXE9ea26g8yrBBM2vgm7DTznUaPMP9M3T15RcpQSlzm6i1XoDnZbYEj8Th+77GmPo9+VrsAb\/9UkuNfv7nUDjU79rogVeUmZmZh5cG7AkxTVPXjJcIHa\/o3J79+5lFuvq1atp48aNyh5xc1KlVVda7d1zKhMWMnyFsJzl1YlQsIH1jSU6qxMunWQUyrv8kpZDN\/8\/Xan2+89HkslPqYG1KG1rh5CrnNStW7c8OExVS0y806x+x9vQIXhw9Km4rCU6eLnvDktjTkekyi1A+dOnTycekid+DHDZwNGrOkZVBNXLZJgL3yrgxY2MQu4VH2SIKHRwmfC4TFHTiZFOqt8hF3TsggUL2OoD8lpBh7R2v1nJ2AlWq9\/Fo7NQti50qDPX3Cro3MjIB\/xoVSGi0MFlgnVUUdPxF2L4qlKlCltnVWlCdC529uLwFhx14AQddpTwIUo1ZMrazA2o2A6FNV7upBZXRLBTGdoQ+\/lCHV7xbhxIjXYfOXKEHU4jDq9uZaTV4z7IFHXoxDbj5VbQYejDzmDMmey0mTzn4cAhD+9IPkxhmzqGYP43YOHlh9oXHEAd6MR3AUAZOrcyCrXOfknvBrq0re2p1MAvG3Pd1p05h62gQiF20KEzWrRoEdAi4nAkV0AcomXAVJV1GlrdNNBA50ZK6jSi9Vr6VQX6scSv6Jsnb+hZ+RRq\/z6HUpqnaVmueJs2dAAHwTJz5swJ3IVgNRxyyPBCPunngFrdp+B2voYyncpQaTpuoHCR21nC+VnT8fkwfHQ\/\/PKRGpQsSCPTrfciusFcGzp8CQiaEYc+uzkY72RxUi8fyyDC43Z4RbnYOq+yYK00nZzHSfMa6II3w7oByy6NNnTy0OpkoTotq3BARBeIG0NCBzqUi0f8YOx8cQY6H0AHzTBp0iQWfii6WZw0ndOWc3lbupt5nQ50Tr43+Ss10PkAOtlVwjvJyTkMXxe0i2roVFm03JDB\/Vuyn4yXcfXq1a8cx7w+doaE7LKxO+3IQOcD6MRVCNm14HYZTISGl+HU8XLZsl9P1lChWK+qIVf8mPKry0TXsW47p8tvN1vLqyhcOOIasp9vEPQ6iXebn09Btm3b9pVj3XeBOW4bFct0stHAh1r4EnUd0bFsTyTezaFDdB+8CuITs8CcSDQ0mmXKfjqnoTqadfPDu+wCnIym80MPJWAdnFxcXppsNnF6kV4C5zXQJXDn+rVpBjq\/9kwC18tAl8Cd69emGej82jMJXC8DXQJ3rl+bFlHoTDSYc7fn52iwiCyDmWgwe+icdqQkejRYRJbBTDSYNXQmGqw\/O+Qy7MtgJhpMDZ2JBrtI2FKmOufP8zKYl8AcEw32GVi74CXnGaM\/U0TUkPACnYkGS3zoVIdLej0sUTtGwkSDfdFQiazpVHoYIOKf7qMNnYkGyx\/Q4STOSpWCD9CJmaYz0WD5A7qI+Ol05nQmGix4YEnk4dU30JloMAOd7nwO+bTmdCYazEDnCToTDfbF7cFPFDXRYJ99j3AOR2R4zW\/QATETDeaspwx0zjIKOYWJBnO30cFoupDRMhl0JWA0na7kTD5tCRjotEVnMupKwECnKzmTT1sCBjpt0ZmMuhIw0OlKzuTTloCBTlt0JqOuBCIKnddoMN1G2eXDDpaqVatSly5dKCMjg10ionq8+JB07xsL1z1gdocwqhzYkZCzXZkRhc5rNFi4hSEeWoh9W4AOj3wZnBeh6N43Fq7IL6v3c1k6RaCFW+aq8rzI16l+7Dp13bvBnArX+V3cwYL8VtBZnVHs9E6nw7OtfuedgPKt7kpzeresXVVn4jlFoDm9I1y\/RxQ6L9FgqgZyweLU9VOnThEOom7bti1lZWVRkyZNgk66VN3fIO5gsQNLdXpmpO4b4x3gdFeaXYfzuj148IBwnxmuOpVP\/nSKQAsXUG7KiTp0YqVC3aBodeWl6uR1uWyet1+\/fuycWzvokBdl8ktLonXfWKjysPsw7Y6bdTrR3Q04XtLEJXTyPEzVCFmwcrCP1WSfC1N1ww6\/0E4lcKeh1U0nGejcSMk+jdYmTrlI8QZpgIB5mPwVy8OhamiVg33sNB0HqH79+gEjIxr3jdlBJ1\/IZ2VduzlY22g6i\/teOXxuoENa8XabzMzMr670lIN9nIwFcQ7Eb+4RJ\/x4Z7jvG7OCTobEzcUpZnjdvj3oyiXRfFddzSkOfbh0BBNkHIdvpelQntgRuMxOvBBEFUfrBJ18rZOsgSNx35gKOl5PGEridfFWvjij6TTvexW\/bFiofGuzHXQiRBxSfuM0OlPeLu4EndsL68Q7ydzM6+zuG7PSdE6+N\/GDMNBpQGdlZcLqtIMOgufDn3xHqyrYx8l6Bej8ntZo3TdmN6eTr56yupsiv0P3f8PsuuZ\/fPVJAAAAAElFTkSuQmCC","height":76,"width":126}}
%---
%[output:085f9d5e]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BsB0_S\nTaille de la serie    : 193\nStatistique T_max     : 2.0589\np-valeur (bootstrap)  : 0.9790\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:12369ae2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:17485e5f]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BsGtsi0_S\nTaille de la serie    : 193\nStatistique T_max     : 2.2282\np-valeur (bootstrap)  : 0.9590\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BsRtsi0_S\nTaille de la serie    : 193\nStatistique T_max     : 4.1322\np-valeur (bootstrap)  : 0.7080\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:15045ae0]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAX0AAADmCAYAAAA5r1VlAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQmcVcWV9wGe0iAuQFgaGnhEwBB11KgfuNJtbLeJGBeCSwJ0CInEbQzBDSO0DkZFxi0Z4w9JAzGCUfIlMYlDULtRXPgyJhATMaL2AxsbREAFsYmtfPzrvnNfveq71F3e69uPujP8Yr9b67\/q\/uvUqXNOddqzZ88eMo9BwCBgEDAI7BMIdDKkv0+Ms+mkQaCoCHy+6yNqnn0O7X7jBbveriNOovIZf6TO3Q8qaltMZfkIGNI3M8IgYBCIFQEm\/C69K6j\/D5aIsvk3\/Lch\/ljhDlyYIf3AkJkMBgGDgBcCn23bSJvmjqf+0x6jLr0G2kndfjdoFhcBQ\/rFxdvUZhDYJxDY+fJS+vDJe\/Kk+k3\/dTF1HfRl6jnuln0Cg6R20pB+UkfGtMsg0EERgES\/ccaJ1Pr+BtcepL4wmAbOfjFvJ9BBu9vhmm1Iv8MNmWmwQcAgYBAIj4Ah\/fDYmZwGAYOAQaDDIWBIv8MNmWmwQSD5CGx\/\/Fba\/quZbRpq1DrtP3aG9Nt\/DEwLDAIlhQB0+u\/edhb1nfowbf\/9PdTjxHHUY\/SFZA5ykzHMhvSTMQ6mFQaBkkFANs386Jn5tPud14S9vjHZTMYQG9JPxjiYVhgESgYBdsTqfvQZ1P2oM6l57jgqn\/Y47VqzjHY8O99Y7bTzSBvSb+cBMNUbBEoRAVXaZ\/1+32lPCFWPedoPAUP67Yd9h675scceowceeIDq6upo+PDheX1Zt24d1dTU0FVXXUXjx4\/Pe7dt2zaaPHmy+G3+\/PnUq1cv4t9Gjx5N119\/fYfGxTTeIJB0BAzpJ32EEtq+qKS\/Zs0a+vGPfywWBUP6CR3kkM0yuvuQwBUpmyH9IgFdatXEQfoDBgwQO4XevXsL6d9I+qUzS4ylTnLH0pB+cscm0S2LSvro3JYtW+jcc8+lKVOmGNJP9GgHa5xXGAZjpx8My0KkNqRfCFT3gTKjkj6kejxPPvkk3X333XTnnXcaSX8fmDemi+2PgCH99h+DDtmCOEifJXxZ6jcHuR1yOuQ12oRWTvYYGtJP9vgktnVxkD4IHuXceOONop\/f+973jPVOYkfcv2E60TXN7Vn+OBY6hSH9QiNcouXHRfpsuQNrHkP6pTFZjPVOsscxFtLfvGM3jf3palr33i7H3vY8YD\/6zdSj6cgBPZKNhmmdNgKyhC5nevTRR4U1Duz033333bzyQOpOh7Yvv\/wyXXrppYb0tdE3CQ0C4RGIhfRR\/c9WvEPr3vuY5o77kt0a\/u3sI\/rQvJVN9POJh9MB+6fCt9bkNAgYBDoEAoiyuWv1n8TNWbv+tpzem3sRde5+MJXPWEZdR4zqEH0o1UbGQvqQ9K9e8jrdf\/GXqN+BXW2s+Pdb\/v2LdOsf3rbfYzG46bdvinTD+3an311xtMgn7xiqR\/ayF4mn126lb8z7m0h\/+3nD6PIxg0p1PEy\/DALtjgDHzjn43GtDhUyQ1TvoDEfcxH+\/9+B3aMCP\/sfcmNWOoxwL6X\/8r1b69sJ\/0P9u2GGrcV59dyd9\/cHVdNzgA+nSUeX06KpmQeJvv99C05\/4Jy2sOUIQ\/bTHX6emD3aLd7f89k0a3vcA+tYJ5aK8KSdX0JEVPWhi3d9pzkWHCZjkvO2Im6naIFCSCMiHsWHj5MihlT\/d1mTflftp01pD+gmYNbGQPvdDlsjx26+m\/FseaTvp9JEHqp87zx9O4x9+lWafN4xOH9nbVhdBNTTjt2+K3UCPrl3sxQBp+NmwrYXwb8AhXan84P0TAKtpgkEgwQhs3UC0dT3t1zdNqT5pu6GQ8N\/72Xep59gf0Oa546jXxP8KJemjQFyMLqt0QP74O+xC0tTURPjXkZ6KigrCv6Q9sZJ+mM5B0od0f\/5X+toSPRYHqIDq39hGE0YPpEUvbxQ7ATzYAVSN6GWreED2VyxeSy+89QFNOXkgffeUYCC3trbSrl27qHv37pRKFea8oRh1AJti1FOMOkxfwnxJAcb\/D7Op0x9up7LDK6nvFXV5xI+aWdqPQvrheuCcC2Q\/ffp0WrVqVZzFFrysUaNG0Zw5cxJH\/O1K+vLhL\/T5rMYJQvor3\/yAxv73X+mnl4ykUUMPDizpt3zSQps2b6by\/v2pa1nuPCLOGVGMOtDeYtRTjDo6Ql+qv9qZMhmidJpo+TOfe06XxGG2dQO1PLeItj8+i8pn1VO3wyvz2h+V9ON2zmLrLhDowIED4\/w0C1YWFqj77ruPYM3G3ucFqyxgwbGRvqra4XbIB7Vy21jC50NZPsQNqt5h0v\/d94+hk4cdErD7RC0tLdTc3Ezl5eVUVlYWOL9OhmLUIYjS9EVnOOw0UfAaOpRs0m9s9K42Sj1BOhSknk\/+0UDNs6piJf1COWcx6SeRQN3GJ8ltjoX0VcL2m6ggfOjqZb088vBCEOQg15B+Du0gH73fGLm9L0YdHWEBq6qySB+PIf382RK3c5YXgWYoQw3UIBqwntbTTGp7GXvYuR4l3z5B+k4mm06gOe0I2Dxz5+7PbCevmhPKbZt\/OQ8Oh+XFwpC+If2wH2cpLWBBF8pCSPphx8EvnxuB1lItzdr7f+qD32Ty50t9ZGdBL+9veImvXr2aTjvtNIKnOM77tm7dSuvXrxf3Pzz77LN09NFHi2o5ndqGkid9dNjJOctvMON4b0jfkH7YeWRI31m9ExZPOR\/b+u9+4wX757Bxd5wIdAEtoBqqEWVPokmUprSQ+Fnqr6M68TsekP6yZcvoyiuvFH9\/8sknNHv2bJo4cWKbW9\/wHvUxwf\/kJz+hM888006HvPfff7\/wLH\/jjTfsdPsc6XuFYXDT6ccxsVCGIX1D+mHnkiH9wpA+E36X3hXU\/wdL7OHBxSqfbW0SXrqdux+kPWwq6UOlM5SGivyqVA\/Sr6Iq8W4P7fElfUjwCAGCRw4TAgn\/Rz\/6Ed1222101FFH0dSpU2n79u0iHQIEnnHGGeIs8NVXX3U8rN0nJH3tEYw5oSF9Q\/php5Qh\/cKQfqGsd\/ggl4m9kiqpnurbDD9IH2lY2ndS7+CqzurqaqqtraWZM2eKu5oRTwrPkCFD2kj6rN4ZO3bsvi3py2EWvr3oNceAa0bSL45VDSZrMUisGHWYvoRbxoKMjZdOP1zt+bkg1bes+R871s7uN1ZR8+wzqeyos\/Kkf526VKmZdfmqlM9lqe9V9Q6nc9P1jxkzxpC+zsC0Vxoj6RtJP+zcC0KSYeso1gIWtJ5Ckz7ag6Br23+Vs6bp+Y1a6jnulsBQqqTP+nzo7CHNqw9L+tgFYDfgRvo4sJUlfS7HSadvJH0FZbeAa4FHN0QGQ\/qG9ENMG5FFJv2zzy6z7e7r22oMwlbRpp5C+YIkkfQjgSZldtKPd6JOIgUTOyfnBQEHu41kOVC4kT7ecdmcH2ofqHeg58d\/ZzIZcaXntddeKzzeod654YYbaMOGDULPf\/nll4t0sOqRn31Cp686W8U14H7lGNJvS\/oTJgyipqaU8BbdFwjMb464vZdJf+TIHOn72d2jPBwVwkwfkWv81oioOwr2CfAbzyD1FFrSL7T1jmyuCYKH1I8DXpA+HjfVT9i5EjRfyZN+KVjvlBJRnnpqjvR1CCzIhA5CLEHKVdMWox4nSR\/t0MEMtiNM+j4OuZHPWnS9f4NgVkjSL7T1jirVy3OHF4D2dtIqedKP8nFHzRuXpJ9UotSV8uTtPRawlSut4HGVlXrSvm49QYglytgWo54odbCkLxYJn45GqUfsKjS9f4PUU0jSL7T1jgw3e+TCG3cIDRE6fBB\/ez+G9As4AnGRPkv6upJekC4F+RjVcnWlPJn0EUcoqLpCt54ofSkWZrr1lFJf1PH3OzsoJOmjLYW03tEd3\/ZMt0+QPl+ksnztNkJYhTvPH0GXP\/qauPykkHfjxkX6SQ24pivlqR89H0zqLmK69ZQSUZZSX5JG+mhPoax3ULYcc8eJ3FWJ\/84776SHHnpIJB0wYADV1dUJL1v8juf6668X\/8sWPddddx098sgjwvMWNvx4cCC8cOFCOuyww4Rtv\/zIZeL3kid9JvyKQ7rSpJMq6P6n19O9Fx9Gv3ipWcTEd7sbVw3dEOa6RCfS11VVBP1QwkoOpUQuUfpixiV4JNdCYeYm6fPlJ5jrbiaWkOJ3vWQ5MhXj3lu\/MAxO36UchgH5V6xYYRM7k\/eMGTNo\/vz59Mwzz4h3CIGsQ\/rI161bN5F23rx5dPXVV4u\/5afkSV822Xxvx6c26SOAmlsgNlj71L3UTHJgtbiibOqqKgzpB1\/GopB+McdFhyyj9CUIclHqEZjBFD1NVAlrLI+Kg9TjRPryNYeoxuk+26j35wbBjdN6kT5MNmUdPnYAsNNXY+9AQmeyltuA2DrHHXccLV682JbeYbvvJekb0s8iyHfdXnv6EJr\/\/Eaacc4XafzDf6OTDz3YjpbJYEPCH9a3Oz319y3ip7njvmRfih5HPH1dVYVM+hMGDaKmVErLBC\/oxA3yMQYtW05fjHqi1BFmXMKq3XQWmDzrnbIybRPMoGMUFbOGLOnjeNLr0DhIPU6kDyl\/28If0MDZL1KnsgOpefY5pF6OLi8MXUeMcoWi0CabbI8PW3yV9BGXRyZ9Vs3U1NQQIm3KqhgOqAbnK+wGoM5h0r\/rrruEjb78nHvuuXTHHXcYSZ9BUcMm337eMPtaQ6fZgYVCJv0oN2ddf2aaxh3bP\/DNWbtbdlPzpk10xtDBgvQHfraH1n3ufROS3JdzzuginHoGDtrjeoMS19G\/Xz8q61aYi1rQJq7nu9m+4Df0Z3mA\/vgRWrH7EhYzvtlKfPBvOY+n3JevHHSAIP2g4++HlzwuofvSubNom+iLx1hqj83WDdT6+gra8tOavEtUQPofPnmPCIiGB6Tf\/egz8rxoZfUP0nQ\/YXybsArFMNkMSvryOMmeuIjnw1E0od8\/8sgjRUROI+nrzOwQaeIkfVQf5o7c3bt305YtW2jal9L0ymbrUvVRl31Gv\/zNbq0ejTm2G218p5Mg\/RWvfOKYh+vo06cPde0a7ErGyw7pSk2dO1PF55\/TLz\/wbhPXc\/HhQ2nT\/vuJtoDEVmxzbpfcWN161L5c9vWu1PROZ6oY9Lk2ZjrARsFMp3xBxtmxx7h8u99BAmc8enidRU2dN1DF54Pplx\/8j2eVxeiL2h\/PeZa9Ixd55OsSdUhf7iiTu7owFMNkMwjpcyA19pqVwyPLpI\/FYPLkyaKL2AG4HeQa9Q6RUM046e79wjOopD\/2p6sprHpHviO3OisZYTvsJ+XK95ce82\/dbVf81974TIs7WNJHYrc8Ue5I\/fJ+XWy1w2ufereJ65ElfdEun3xIo1uP2pcvj7B2OvAW1cVMB9gomOmUjzRR6vjyfsOEBQlUC699+qZnlVHq0e2L3J+aiZZzHh7HO3xd7sjVUe+o7cGhLh45hDL+huXOjmfnC1VRl14D7cvWDzxtcuD4O0F0+qjbSb0jW+8gDUfsVOPloy6kNaTvMvO8PHE5C9+KdcD+1iSUH5n08XtsB7khvSWDmjm6fZDqIWIQXataZrGcgHTrUfsSRE8fiMCy9\/0m1VNaBPXK3pdYWdPoGe4ijvEPEu6BHQ2Z9J08jMMe5GJh2Pni44Lk1QvUde7ITX1hsL0Q6M6HqNY7uvXEmW6fst7pd6C++kIlfXkRiXJdoi6BCemoAJeJq4eIhajDaYIWo54odYSJVxPWU1qnrih9ERKl5uXoUeqxwj1YvYGvqVP8eJ4LXI\/saIh3uqSPtLLOvu+0J6jH6AsFuW+aO576T3tMSO2yyaaTTj9O8kRZTgQa1E4\/7jb5lbfPkr6fescPOJ33bs5ZOiZ7UUlf1CFZVrA5nSr9RvnodTBQP\/qwFi86dUXpS5h4NWE8pYG\/jsVLlL4AK91dTpR6QPcN4pYoS5XEkSOjLvqF9sjVmUu6aZJMoG59SHKbO+3Zs8e6UyzEE1W9E6LKNlncSD+oFDZh0ARqSjWJD8tLmpIbIOoA06etaItu5nRRPnqZXOKMshgW+zbqnQDRJou1A2Pb9nSl1ctCjYsuhpHHn6rE+YHVF3ejzSD1dETSv+aaa2jUKHczUd3xKEa6jRs30vTp0x2vUixG\/V51RCJ9LrgYEr1bJxw9crP6VmG\/W1XvGTWRP5RTB50qSB\/0XUmNvuFybTLOSvoyuaiqhSAfo1M\/gy5gQSV93V0R2qb2JYj0HmSy26qKEP4TYSTws8vOtg9mdRd93f5EHf+g9bDPCfK5nQV0JNJvamoSBLpq1SpdKBKRDgvUnDlzqKKiIhHt4UbEQvrt2SMn0sfpPVtWeElGMolB0l8pSB8fil8u7x6rRKjqWv0kdrX0oCTGuxarL\/47F3lRQdvYGscpFr+bpO8lUYeZH7nFOOc05xfNMmg9cl9Glo20pOlMmvwOZrkenXMDp4UySDtxYMxz2W8xUjFj0nfCrdCkLztndR1xEvWdOo82zf0G9Z36MHk5dblhA+LHv470gOyTRvjALzbSh5ftTb9ta7rWHnfkWtelWdthP6ndyStTl8DcPnpVjWF\/jCHi3If76HnXYpG+3xImdOAzq4jSFulhdwTydzoAjLKAhemLLLXqkn4YMoakz9Y46arGWGPqR5H0wwkw1kLJT7FJX3bO6vm1a2nb47XU79ol9OEf7qVdq\/8knL86dz+oI\/F3SbU1FtKHeoe9aRe80ERnH9GHTh\/Z2zbBvHzMoIKB5qrTz+qawf2VNe4x5fFBIub85s1lgW6ayh0W+ltWOFlV6FzUAdD4o9dRO+VUItb5RO6j96ZL+bAQpA\/SEwufQ7YoC5jdFw1pOhpR6lm8tNm1aMatVyV9PyEhSl94obTq8B7HIPUUUtKXrX1a32+ySX9Py448K6CCkYIp2BOB2EifnbP+71\/eo3XvfWzH03ELuBbXuLha78DqIeu7nq7KJzBZh\/3UUy00fHjKvl5Qm4w5EFalv2VFJP20tGvxk9mDfPQy\/mwWmObjwkyjHdzLzecglFWNZNvuJ03n7cAC3l+ra\/ES5VA6yPwNOy5B6pDVSDrqvUKSPtoCs87PtjZRr3Ez6YM\/3ENfmDiXmmefRfsPP6GNI1fQfpr00RCIhfQ5tHLViF500vCeNHnB32n+pCPohXXbqe6ld+l3VxxNQez3g3RJJX1b\/12f04Wqh7myDnvt2pykn6mronSlpeJI19R7OtzI9Qi1iCSJqWoMXf2024Gqm9WLEyH77Vqc6uDyOb4L+sLWSG4+Bzm1i5Vb5+wA5YY5n1gJXW52TPwWZa9xkeuH+goLfnNzM\/ECpmOJJc9NXXVVFJVYkG+hrVGCu3qv0KSPdqtxetxCNQfpo0kbHYFYSB\/NkC14IO2zfv9XU\/5NqHoK9aikzyRFjUMtHbWDNc5Q8HqW2Nb2byFY9m0WURazR7BZFYcXwVhU50x4qh5WlfTd1AEqwfpZ1Tilb2iwkHbSyQuVlMd7+QCa2yhj1ZjOWe+cmrWqgQJKoOmjssmRcc7EVSdMsGxVlc40Etrg9Vjmmtb5RGU6\/xBb7T8WfJC+7cHaaOWzzkH8rjvPqd78zk1UlRjar3uNZZDvhlWVq+442xJeso+TWqgYpB+k7SZt8RCIjfSL1+T8mvAhP\/eXFjr1K1b0ygZcaAO7+XTOtllVizC5gaug+slJeFIeH4KxynBW7agu+rZE6eMLoBKsvYDVWzsQVZp2TC+H4VUIUo7NLiT57Puc3jhn4MfaYxkrnI3k+sKHhVWUyYYjgMrGzfontxj7+zVgHOXziSBWVQKTtPO4qP3Hgi9I\/7NsX1zyqXM7yG4CFlA26XM9mHc17ipHJ6spne8L9eiqKg3p6yBammliI\/1X391pq3VwPSKseQqt2sGQyNIu\/pa36CSZu8mSW546g\/X+EkHyLsDr2GyoKLvBjuWtSlNOKqRVS0baFjJ8WMrmm05epEzq5EZiSowhi\/Cykr7klMQkJbyHsw5LsjNZ7rDYUgeInULW1NPSj2fLzBAxUcrObJlMvb1zwuaHTT7lBQB1C6ixytrepTlp2ktVBaGb5VavQ3kxH7I7MJxPWLuV3Ciq+KzNhuDgXUsaY6oRT0ddbFVqcFOJ5XZH1jSQdy26vhheNCRIvyll3wvhtSsqJOl7xeAJE3unNKm3\/XoVC+mzTn\/KyRV5qhwQv9d1iXF0W\/4ABelLpGLdpNnWfZ0lW7yx7iKymArqAFBFLo\/7Ft\/PlE5VIYmPcVg12xQRZSwyYonbIqR8QsypW3JafSYxa5HIV0dYC5FFdrKKIrd45BYEixCtJ4dHTiXAKouc+au1DIAooZrgBcxaJBpzpC+pzsR4ZPuIZFbpzuogJ6JkAhPjk1Uh4awFajf5bED2e2AMnFQ0chwbvF\/b8pQl6Us7MBpaby9acj35deT65SQYuKvEclZVavucFhJ5IRRYZqOZYifg1H91TqSrcn1Rdw+FJH237xqHuz1OHCfi+Zin\/RCIhfTDhlbW7bZ8OYt6MYtNdLaIa8ml+Kjwf06REGXJ1qLdHFFaRJVbKFzj68gEW9UoPkjxZMMyyPp+kMvwVLXt8StfeofawfUWgWdlapwpgEilRUD9gJ3UGKrEzguEraJh8hQCv6XvziOW7OG3TTDoVz2fjQCYShLtEDuGKmJpGhhaiw1vl\/L\/m3cPFq45qZ119JafQDY7TGxrLTXS8JQltdoLBXI3NArVCMc8khdOa0nJOeYJDCWvVN6diXHKpOnkb62jRYveIT43EO3M5PKgHrVdyCovtkDR7m12l5Opz\/k8VNbUi77IC2VuUZZ2IQ670jYCjSokZKdcvuCQU22pfZEXjSh35Op+t2o6tzj7Ycsz+cIhEAvpo2pI9XOeXk+\/mXo0Qb3DcXmcrksM0lTZBwD5pj\/xT1pYc4RtDZRP4FbJTDIgPVvKBU9l7xi1pFdWWuTnEQuFvQiAzLMfcJYo8BGrUpNVR1YatdVEuYUE5MNBs7i2HAZZEZgJ3yZOFJRbwEBSTKwW0ViLBNprkY7V1lRqpe0FyKRnt81WkmRRAsFl+yeWmbS1sxEqHV7ExMJjK3isHUqbHYlFtO5PbjGwxsbqi1iWq+otxzDp4DGdQQ\/rxELZWtGaXShsBY+lw8u2gQkc5UAaJ2nhkhcjCzPlqDrTSBUVrZRKSTswXsAgRKCeSsl2qiHbj0ppXKS+iM2IOCiSsGiopHQ6f7wYJ6j4ePGS28YOhfKuRUy\/7Ejb48+Yi0P0esJiI6sc7f5zX1j1l115Ww\/P0LD\/k7bns84duUG+XTWtHK8f0TrN0z4IxEb6aD70+l9\/cDVt\/\/hT0Ru\/6xJ1ugwpf8Zv3xRmnz26dqFvL\/wHyWqkQ7edQG\/3etkuSiZ8\/CirPPA3S7i8WKSyl03AXTonGefUJDm5NFsFpHBbb25Jy3kSZDaZRcY2c9p6bC9yzOVpu\/vIffQWYefvTnJ1pVam6PnUL2j06NFS3\/PVNmq7eKF07j\/Xlt8XOY\/TIuo0tjw2eX0BIckk6YCfvFPILZrKIpPd6WFs+DzCaUeXjzGXlrPWcdsFqv1x7IvrhM4fL26D2B1KZzDyuFg52u5APfufxZFVc237ojQQ65O02Ie5RMWpy146fQ7XrPPtmzSFQSBW0i9EE0H681Y20c8nHi6KB+nDH4C9fGGyWUWV9HmvzdR51\/sizfNNz9tNuWT0JXlSKEge779V8S3790GnDiI5gp9THkic6sNlOaVHWiwkMskjvbwYyGXKiw\/n4d+4vW9W5MJc8Du1Hvx+99K7RTRCuV1u6ZFfrgd\/O\/XHqS8y1jKe+N0LL6SV+yK3Qc3HbQNuQfI4Ycb1OPWF54xb25zGSp5HTu1264tav1u7ZFzc+u\/ULh4Xp744LV682Ae9LrEQ37sps\/AIdHjS37Ctha5YvJZeeOsDKnv9t1T2+u\/yUGt+tLkNAYHkvR6nPCAfJyJDWW51+NWtvncqS26rU3r0I+rvKMOrHiYftf9eODr1Xa4nDGZx5XEay6D9V\/uuO2fc5osbxn6YqX3RaVf\/rpYAk+pjbTMM6ReeaJNUQ2ykL1vqvPTWh\/SNeX+jngfsZ+v4w3baT72DckH8+AdJn6X9sPWZfAaBUkcApN9v\/1bar2\/aJn70Oap6p1DXJZb6eBS7f7GQvmy9gw5w8DX8t3rwGrSDfge5Qcsz6Q0CBgFnBAp9kGtwTwYCsZE+E\/3mD3fbOvi332+JTPqASTbZLHRYh2QMi2mFQaB9EHC6IzdMS7Y\/fitt\/xXscPMf45wVBs1488RC+jIxs0oH5A8VT3uR9GOPPUY33nijQOuoo46i+fPnU69evexLlvH7j3\/8Yxo\/fryN6CeffEI33HADXXLJJcL6BQ\/fdemUHr\/FUQ\/X++STT4o6v\/e979H1119vtyuOOuRpg\/JefPFFuuOOO6hbt26x17Nu3Tqqqamhd999l84999y8euLqi1yOilfQcfHCP87xd6snzvH3K4uxcRr\/uKhF3jFs\/\/09tkMWnLO6Dvoy9Rx3S1xVmXJCIBAb6Yeou2BZQDq33347zZ07VxD9nXfeKQjouuuuo5tvvpluuukmUbecZtu2bTR58mRas2aNfa8lfps2bZpjeuSPqx4QGNz\/QfTcjosvvlgsSHHVwWAzIR977LF5ZBxXPTJmw4cPF9iPGTNGLKJx1YFyMJZ33XUX9e7dW4wbsOOFOmg9y5cvd8S\/uro61vF3qwdjE9f4u9XBwo3b+Mf5McpOWB89M592v\/OaCKdsnLPiRDl8WSVJ+iockNYWL15MF1xwAd17771C6oeEy1I9dgL3338\/XXTRRfTDH\/7QJhDkA2mp6Zlc4qpHLQd1ptPpvF0IpwnbF+SHFDh79mw69NBDafXq1W0kfbkdYetBvhUrVuTtVNymZ9g65IUFpM8LMxYZp8evHnU8Gf8hQ4bEOv5u9ci7TbQ\/yvh71RFk\/MNTChHfnNX96DOo+1FnUvPccVQ+7XHatWacVLtQAAAgAElEQVQZ7Xh2Pg2c\/SIZ56woCEfLu0+QvvwRg\/yh1sAD0j\/xxBNtcmUpm6VGJgu39G5kDbIIUo9cjiopx1kHdhR45PbJ6h25rqiYbdiwQeyaVPVOHHWgDB4r\/Der7tw+Bd2+cLm8iGzdutVzHMOOjVqPvFjFNf5OdQQZ\/2i0QnlSPaR91u8b56yoyEbPX\/KkL6tO\/Eg8CulHqYeH0elMQR7iKHVgW79w4UKaMWOGIGNelJxIP0o9yLtkyZK83ZG8sHJ\/otQh78BQnqreCYuZir\/ffImrnkKMv9qXIOMfnVZMCUlGoKRJX90m+6lrnEhfR70TtR4vyY8nT9Q65INPLtNJCo9aj0yUWFBQHh75YDpqHfKCgbKd6uDfZTWZ1\/g7Sdh+8yXM2LhJ8roSPquCgvZFd\/yTTFambfEgULKkLx8gMlR+B7Mq6fulZ2Lhg8oo9dTW1tLMmTPFwbOT2iBqHXKZKjHL5BW1HhkzhBtQVWhxjIsT4am7iSD1oP9O+Mc9\/l71xDX+bnXojH88lGKpdjbOOJG69BpE5TP+SJ27HxRX0aacGBAoSdKXzexUqRaqjUsvvVT8\/Oijj9oWHyxtq6oCuSw1fVz1gKAeeuihvOFkc9K46vD76OOsRy5LNqeMsw4vk82g9ejiH3X83eqB5U5c4+\/VF54Dbot+DHxiF8GHubvfeMH+rfsJ482l6HGCHLKskiT9kFiYbAYBg0ABETChlQsIboCiDekHAMskNQgYBPQQcJL0O3c\/mMpnLKOuI0bpFWJSFQQBQ\/oFgdUUahDYdxFgnT4QMDb5yZsHhvSTNyamRQYBg4BBoGAIlAzpO5mkATU1vk4QJNll\/aqrrmrjHcuWPvCAlM0Rg5Svm1aOpzJgwACqq6sjN+9T3TKDpuPDUfUwM2g5pZCexwN9UeMXqf1DWnh7T5kyxdEyKw48grQnjvqCloHga3hMzJ2gyBUmfUmR\/gMPPJBHiFgI1N\/igrGYpM+LD+zqC73AuOFjSD+HTBCShTUNsPPzGo4yL4O0J0o9YfMa0g+LXGHylTTpq0Ql7wZkxyQ52Jq8O1AlfTV6JEINQNK\/+uqrhT06S374CGH6ybsANfKhvPtwizrJw622De1GDCFEsYQ5JEz92CxSNtfj\/qEcbtuBBx4ozFSxW0DQOQScQyA6pyiVyKeaRaIulvTVdsk7AKd2wFFLN48cFdVtbNA+v3FD0LpnnnlGeCC79RHl8PhgPN3IWa4Lvgx4DjroIFvSd5pbv\/vd7+xIr8D8wQcfpIcffpg++ugjkR\/zCbu2v\/zlL3Y6r3Fzw8WpPYWhi3ClGtIPh1uhcpU06cuSPuKowD4f5DRixIg8UpalMUQp5N0BQAe5Qr0zduxYQZ5MDG+88YYoD2TiR\/py+ZwP7UCwMC4fER3lhUIecFXS58VMJjK5DuTlsrhtr7zyiiAY7lOfPn0EwaG\/CEGtqm3UOpnIkQ7kIy9yIDfGTMZZ7h\/jxwujmofDQjQ1NQlMeFfjNjYcXZMXVnmB53o5kqhcl5NaTIf00Q6EvpYx5PLZ98NvbnGQPx4LtEVutzwv3caNVY1e7XGLpxSURKJExXSy3uk64iTjrBV0EAqQvqRIn+PnM06y\/lvdZuNvSLnQyULnCilWlqRYEmNS\/spXvmITNFzhZfWOF+lDlyuTuRw6mctEO7zOHtxIX5W65fMFXvBYunTbhbipbVTVmBOpMgFx++6++24BvZPzm9uuCXk4MqS8Q1B3L+rYqO2WVRzf+c53aOrUqfbC4dZHdafgNG\/kHYW6c2NMMX9kFY46t\/gdkz7n4zAVTnlvueUWuvXWW20hQ17I1TlVKPUOk\/4Xav6L3q\/7AfWf9phWdEwm\/C69K2xnLP4N\/TBeugVg8gBFlhTpqxI6S7P8sfElJYwPEwmkU770A+\/U30FuYUmfVTEgdvlxUsngvRP5+5G+04FzIUlfJna5T2g7S\/Qy1licvPLgHRZsYKISmqxSk8dG3pFBYnYifV6U\/M4j\/CR99fxGrovJ2W1uzZs3z14QVNJHf7Bjcsr7k5\/8RNwX4LRYF4v0UTcuPtn1khWd1elxugnLbYcQZecQgNNMUh8ESpL0QQKsY5XJ1e9ATda9Ix8Ttkz6rHYIK+m7jQcTk5N1jh\/pOx0qx0H6strHS9J36xOPARZRXF5zzTXXCFWZ021lTuQmH1qrYwM9NqvrIIEXg\/QHDx6cF5ab26xK+jIe8g5TR9LnvKr0Lp8TMem7tSdu9U6fyffRlvnXaEv66AO8bz988p48qd7cnJWM9ahkSZ+JcMuWLUIP66bTV1UzrFNWpXsvnT7IyUlv77TgyDp9ln5ZTeNm6eFH+ijHT6fvRqpuUnBYnT4OJXnHxQHXUDdLxNwOWc\/+61\/\/2paG+YzB6axEHhv1DMRrUfKT9HU+xTh1+owByNlPp+82bsXU6QdV77BzVuv7GwLtEHTGwaSJjkAspL95x24a+9PVtO69XY4t4ntzjxzQI3qLXUpwMs\/kD4otImRrCi\/rHX4nk4ysx+fLQdh6B6TvZtmDd7rWO+iakx28Dukz8XPgLl5wvCRGtM2LENXAaWGsd+Sdi5v1jowd0kMtx1IsS7jAHI\/XuDF2qrorDtKXx1DXekcmdbQdwgcWOJn08d9Olj\/4XdcirJDWO2HUOwX7yE3BsSAQC+mjJT9b8Q6te+9jmjvuS3bD+Lezj+hD81Y20c8nHk4H7J+KpeGmEIOAQUAPAT5EPfjca6nH6Av1MkmpjC4+MGSJzhAL6UPSv3rJ63T\/xV+ifgd2tTvMv9\/y71+kW\/\/wdpv3iUbGNM4gUAIIyKqWoFcVytY7m+eOIyd1jdNBbgnAVtJdiIX0P\/5XK3174T\/ofzfsoN9MPZqgxnn13Z309QdX03GDD6RLR5XTo6ua20j66u5AVhNVj+xlp3967Vb6xry\/iYG4\/bxhdPmYQSU9KKZzBoE4EICE\/97Pvks9x\/6AQNq9Jv5XKEk\/jraYMpKDQCykz92RyRm\/\/WrKv9GRFT1oYt3fac5Fh4nFgJ9pj79OdS81U80J5bZKCL8N73sAfeuEcrGITDm5Ii8\/8k5\/4p+0sOaIvB1FcuA0LSkFBKqIKENEaSKq7wAdat2SoU\/fy9B+fdOU6oNW5z8s7SeJ9HFehn+l\/sCYAf+S9MRK+rodg4Q\/rG93eurvW0QWnAOwlD\/7vGF0+sje9hkBzgNm\/PZN+t0VR1OPrl3sxQBp8GzY1kKL\/9xMg3t1o1FDD6byg\/fXbQZ91voZ7dy5kw7ocQClUvGfNZjyvYciqfh8eb8ugvQHfraH1n3+uWsnktL+nUtvpR1Lb6Wywyup7xV1bYg\/KunvfmMVNc8+kz7f9aGNRRS1Dsh++vTptGrVKu1vtaMmHDVqFM2ZMydRxN8upC9L+zLpyzsCLAz1b2yjCaMH0qKXNwpVDx7sAKpG9LJVPCvf\/IDG\/vdfxbspJw+k756iv6ru3r2bYNIJa5GuXXNnEXFNMFO+N5JJxeeyQ7pSU+fOovErtn3i2onEtH\/reqKXH6FOf7idymfVU7fDK\/PaHIX0YW\/\/3tyLqOc3avOiZCKezvZfzaSg5wRoGFtTgQwHDhwY1+eWuHKwqN13332OFnnt2djYSF9V7XCnhvftLqR0+YC3EKT\/00tGBpb0Wz5poU2bN1N5\/\/7UtSx+0jfle09tg098+LT8YwVtn\/3VWEnfz+rHyQFLh8ziMKHVqae90yS1n7GQvqqa0QUbOnxZ0oetf1D1Dkv6v\/v+MXTysEN0qxbpWlpaqLm5mcrLy6msrCxQXp3EpnwfUjP4ewIUZP588o8Gap5VFSvp+5lq+r1361xSyVDnmw6SJqn9jI30nUw2\/QCSSR9pwxzkGtI3i5bbPAtCmn5z1el9ksr3Iv0wfUMeP1L3e19I0ldjMqEurwt+ZEdDpJWdBuV3US5dUvtb0qSPzjo5Z\/lNNpX0ZZNN2apHVh3BIogPcVG+IX1D+ob0ifZF0l+2bBldeeWVYvixCCxcuJBmzJhBTrGHQMDr16+34z7h7xUrVogAf7W1tTRz5kyRb\/bs2TRx4sRYbqYradL3CsPgpdP3WxR03hvSN6QfN+lXURVlKENpSlO9h9HmviDpb5xxoqNTFmMexoonDjIEyaukj9\/OOecc+6AYbeSIuYh5pZI+3uPuBbkchMQYMmSIHe5bh4MKuaOJUr9b3ljUO4VomG6ZhvQN6cdN+kNpqE36jdToOhVLnfR1v8Gg6fxIX8dPwkm9w\/GmEJb6zDPPzJPWVfUO2gxVDkKmq6SPd8cdd5yIiQQbe1j4ff\/733fcQXj13a+fQXGLK30k0pfDLHx70WuOAdeMpG9IOW5S1p38YUmZJX3UUyjS1yG2IO0vhHpHF+eg6fzIcKjkHOe25KqSPgLisWoG7eH7MdwkfQT\/g1pnwoQJ9NJLL9lqIpb0EZW3f\/\/+YlHAXRi4ahTXYwZ5\/PoZpKw400Yi\/TgbErYsI+mbRaWjLSporw6x7aukzwuiteg6P06kj3sNoKPv1auXnQkkjgcqG1m9A9LHBTff\/OY3xWU1qk4fUV7\/9Kc\/0e9\/\/3sh9SNd0HsKSpr03QKuBSVy1daf4+x4xd4xpG9IvyOSvg6x7aukr8MbXtY78rWbsqTP13iifD\/rnZUrV4rzAEj4Tz\/9NB188MF0\/PHH6zTNTlPSpI9esrlllGBoThZAWFDYUxf1qLF3DOkb0vcj\/QmDBlFTKqUfS6eqiiiTIUqnierdo+8EIWW1jTqHxUHKLyX1TiBmLVBiHPwuWrSIevbsSfC8vvzyy\/N2EDrVljTpx2W947RwQMr3ir3DpH\/9mWkad2z\/QLF3drfspuZNm6h\/v35U1i1+56xiln\/BQQfY8WKWe8SL0ZmsnKaY7S8k\/mcMHSxI3y+WDve724hhgvRbBw2hT996yxWyKPh8eb9h4rB44GdDaN3nznUEKb\/19edoy62nOTpnBRlzt7QcjkF9317WO3H0qdBllDTpxwGeunDwAfCrTTvtC1hQj4m9k0ObY79M+1KaXulmBZoDsXnFiwkyVomJLROk0VJaGZ9N++8n3uhgM+DrZ9F+71hX\/a1\/5TV30o8Qu+myQ86ips5WHSu2OdcRCP83nqdO955VENKPErvHCbykkmHIaeaaLan9TOxBbtCAa1Fi79QMzm3\/o0jJ1Z072yF5UU6xYsuwJMtBdV\/79LNY5m+x2r8vxj7a76unEWXWE6WH0KfPPOs4XkHw94q9E3UyhPW8das3qWQYFSc1f1L7GRvp80Uqy9duI1yAcuf5I+jyR19rE0dfF1hW61x3xlC660+NrqGV49DpnyrpfN2tsv1brlpkBNHJyqXrmPMhPZcPnfXKVLUwdKv0cShCvqDld\/TYRBMGTaCmVJOvs5X\/COenCDu+opShQ3PnBo3Osy5I+YXW6SOq5u53XqP+P1gSFKY26ZNKhpE7phSQ1H7GQvpM+BWHdKVJJ1XQ\/U+vp3svPox+8VKzCI+sczeuagHEIRquO2towQ9yT5VIoTHCtRmqRUaQj1aeLzrmfPmkP4FWplaKIuBF6mVbLvhGww5aLj+JpK+zcDH+pw461SZ9P2yCfPhhx1fUwYfF+O8Ckr6si1fDI3Nf5cvPO3c\/mMpnLKOuI0bZUMhXLqr4tJdOP2jsHbQbZpqTJ0+mNWvWiG6wZQ9MPPFb9+7dHcMvsHnn1VdfLS5+YR8AdgbraDuaWEhfJuz3dnxqk\/7O3Z853p3rBpJsmik7dRU69k5cpKBaZIQlBR1zPpmU5fYLDvHwIrVIPxdmwGuRC9t+XeKMUr7OwpXbCVmSvg42um2PsijqLFhBy3eS9EHW7952FvWd+rDo1nsPfocG\/Oh\/qEuvXAx7v\/DJQfDQTRuHBOwUhsEr9g6Ie9q0aXTTTTfZxI52LF68mO644w6aP39+Gy9etT+yA9jw4cMJpqFjxoxxDdkQRz91MQ2SLhbSR4WQzJs+2E3Xnj6E5j+\/kWac80Ua\/\/Df6ORDD7avQwzSMN20cah3PptwKqVwdVs6Ten68Aoe1X0\/LKnpmPPJpMDqC11S43byvsDtWkDd9uuSmDqmuuU7zQWdhTFK+TrzL2z5OgtWHKQPKX\/bwh\/QwNkvUqeyA6l59jl08LnX5t2TKy8MsnSv9p8vTUmSpB8k9s7y5ctF08ePH99maJmcIfkjhAOctZ588kmRDqEaqqurbUeuRx55hCDxw1FLDeKmFlzypI8OuzlX6XxAYdPEQfqDTs2RvttWW6d9qvt+eFIIEfvl7LNtHXFVfb3nHa\/cTlwJiH0BSN9pqfNqv0z0KIfvlEU5uotAWHzQap2FkcsfNGGCvah72d3rjLGcJmz7495p4Y5cSPpbflqTZ70jX3KCdoP0ux99Rt4NWKopZvcTxrfR28sLw\/bf30M9ThwnFg6ohboO+nJeeToYxkGGQWPvyFK5rOZhFQ\/CMiNeDx55McHfsveuSvqI1HnRRRc5xumJo586eAZNE5ukH7TiuNI7kr6Pc40sJfZraaGV1dUWKQgO1JP0HUlHqTcsKYSK\/TJypE36Qxsb80hYxjqfrL3VPHL7zy4ry1tIZGkVFMxRKaEuKoQkq84XnaBoNunHtKirbQg7vpmq3CGu187SC3+5Ldt\/VUvbH58lfpKvS9QhfbkcVvWoC4NsvfPRM\/PtA92wVj1+ZKizoAeNveMk6UNdw6EbmPRZbfPQQw9pS\/o9evRwjNPj18+4ODBoObGQvlsYhrjCM3h1ypH0fSwjWNJCuRWtFbSudTmppOYHpCPpKPWGJQX7kM\/HIxQE\/mZrKz1fXU3phgbLgxQ6+8Z0HgnLfUHf66oylM4QZdIZqqpve\/jLC0NFaysteucdcbPYyLKz7TIRcLghK8+L8MNVKMzyYAWJxS3JOo2FF3Gq7bcl\/QCLuttCCf9cJ3wC3bzmMz\/98MfCKi\/ey7dkaEf9QkH8Kun7qXdUbCG945GtdOTFoPtRZ1Lz3HFUPu1x2rVmGe14dr5QH8nnBH7fjh8Z6izoQWPvQEWj6vRRxu23305z584VF7CokTk5KNtVV11FS5cupe985zt0zz332PH2efdwzDHHOMbp8eunH06Feh+J9L08cbnBMN\/Usd7x6qBf7J0qqqSKYR\/mTPKyEjekz6pG0FIuLjqkiJlVDYL08LRWVFDF8nU0UpJkdWR9R9JRLDJ0SF\/eddi6dQdSkD9ytNtSqVi\/Ng61SFyQPqT8rCSZSaepNntGweqXDA2lek6P3Gn8S1NNvaXmAakxaWNBfO6dRVnSH2mTPpaJOsYQRC8akw1bINXPi4BjHz2uq1QlPUd1kYRRlbSzASU3CHSsBZ3bL5Oy2h55p8Lx8712RE74OJXPeLZpvzJP1Pd++GM+cx+t2Y1ZkKF+\/8iI+P98MbrOQS52AztffFyQvJcTlirt41J0PFEuRne76Upnpxs09g4sdLysd0Dg0OPfeuutdMstt4jImnhknb6b9Y5bnJ6SJH0m6kJK9Dqxdz6\/cCQN\/mhTm4NYJ2LGb+mGLONnSb91HUg\/J8mqFi2O2802pGN9uvIC40b6cnloCUveFknWtzHnQ8myZG2RfkZI2CBc67PPPlBPcduQr9J6Zy8KaCP6jwUCRJ19QPxcP+OGBbHLoucE6Z+dxQfJ66soh2E2n12WXH92EZJ3F6ijpt5qP0h5UXZRkUlTHbeMpDYT+OCRiFO0PZMRC5gYA2lB5\/bnl5+NrSOAsShT3qmIUmRVnfLeCR+n8tGe2noskLlF0W6\/JOXYdYnmWLs07s+3llfQ89Wp3Fhld3M23tLfVn2VeRe\/yDp7JmhVLSObbDrp9OOWOJNKhmH76RanJ6n9jCTp+5F+HIuBTuydwcdmST\/7EbcOGiy8HDmGisV6acLvKbjXZ+OqpN5ZL161nnIqDX9ug02+MokiT1OXXJ6Bz1hxUmyPSiJqOhQ62garrEFDaFgXS2fbMmAgNS5Y1Ca2z+YzrNgu\/DBJcVtUD83NZ5yeV36qSyeR3yY3JUYM2pZ6\/rk2cxht47z4b3iDwiuUcWCc5LbtPH4Uffr0M9TzvK9licgi2HQWQ+RR49PI2OAd+ssLLeptOtRqP\/DZs+C5PHzO2a8L1Q0dKspnPOQxc4qFkzfOUq8F\/g355eO1nF7GxHWeZPFW24O\/GR85dpBc\/punDKFh71j9BSl\/881Tafnn+R64Xu1Bnrz5IY0hzzceU7FYVUazQAtLfEHyJZUMg\/RBJ21S+xmJ9Iuh3gHpz1vZJFREeNTYO0Ligz5becTHnCV1+ZUgO+nhNOrHxUk4vb1AZPN\/Omgwvfub\/yHEaen24vOuc+CDY\/+PSNe1a1c7zZBjv+zYNk7wyYmniDx4\/MoHsaW6pNrEiOH4MV4YIK6MX\/loA9qDWDRyWYwLx6bh+hgXGRDur4w9lwV8UqkUlb27kZB3Y5dONExasETQs0GDbYy5DC\/85bFnfNQxl9sil89Eit\/UPsvv5DH3Gy+5PZhn\/QZY4yuPkTrP1AmFfBWf5ea0a\/uzuysdUvJLwzuCL9T8F22eO87x2sT2cs7ya3sS3pck6TOwcUj0boPkS\/pZVYbY4gs1Rk6CdiwznaZP3nhTvNrvtEpKrbQ8Wds86jZaTZAtx5bSPNJjJwHpnSVgQXicXi5XansbEpDVMYpqhqV2uQ4hxRNZEr+SHn8zBq7tV+vD3\/LiKpWBelRpleu3dxLZ\/ordlsMuxAt\/xwWcd25cltpeFOg1F7LtzxsTqREYMxs\/MbFy6jCn3YEjnm71q3MnW34eNuhfa2vOqiydpoZ0hipZvlHaT1utXWuqr3Wuk+QnqWQYN2ZJ7WckST9ukJzK81PvsH4VKhmYKtZBZ+\/xsQvVTfajgM6dKiupbPPmXNWclyUmaVFpU66kS7bLlQ6R5W056hEkJJfPRJK10oFuV1jhKI\/YhaDdCrG1vvlmHimI95xfXoSyevdMtm6xPHLbs4eeavsFZwp1t3VOITCTDyAVIhRNlnF30j2r7VPTy9hw+XJ\/XMaV8clrL3TjWZWb07yy8XTCi+vO9snGX+4jynfCp6FBHIyLVxgyPj\/K6urtOVFZaSEv8JWerBoIJbSsXUtN1cPFGNsqR\/E++5d0ftLanWjDOfkmm1G\/T6\/wC1y2kfTdUS550kdUzJt+a0nQ8hP1jly\/g9xDt50gDr4G9yoT4Qdkc0znjz0XpsDpoJUP1WxiFB+xRboNQr+a06Xb5aex4OSHP7CshCCZKYelMskyiWQXGJz3zbQXLVjFZPM6lI+saH9T1sfAJg6JFPiAU7xrbLTN\/MSf0kE02i5+U1y0Tmk9hW6rXmkdFqcr214oIh0Yi0VEIlomN+twVTolkUgtjxQZTF4s1ANheRHKwmKVbZWOPuB\/xSFoti9oP\/6G+oh\/l89rGoVDgaVr54VNFiJkUkf59VU5aR9Vn7KuIocPqDtvvCw3NbZyQvmwJKsbmrUcS+Nvy5\/Cai8Oj62zIO7P2pa1VEmVtLnMEkoapfoh4DRmVZtC9fNe\/KQvfz9Ojliy5U+QBSapZBikDzppk9rPWCR9mZgXvNBEZx\/Rh04f2TuW27QAbpDYO\/lOQ85Dw5tfJ9LPc++3pXwSHyyoQbaasXlK2j3INZ7S2kq3VQ\/PEb8g76y1SNYKJk+SzouLA6ue7CLjUj7aX0lEm8usC2Dq80jDMsGUSUVum0wwXgG\/UsMtSVM82K1IN0nZVjZM641YdC1vYhAqWwxZi0qj3T7gWCV+y5l+2hKwpJNW28gLurAeynbMInvneEPy+Oass3Jpc+232iMLDdZo5y8i8txA69e2tFAePizdZ8cr3+rKKt8L95x\/g5XW1+RX2oWe8v8ytOL78Ur6PF\/cnLAK5ZylQ6hBA64xActlc8A0OGnV19fTOeec41g1X5Y+evRoEW8HjlvydYtu7S150r96yet0\/8Vfov\/7l\/do3Xsfi3g7hdT1M9Cqc5ZOTBbO6\/dRqaSA+JWyfTxTqltkS5Q\/PDWcflHdlDOtlPSt6keOduUTC3u6+pMaTAadynObkDoOMNZOYjgNW5klfeWQMD+Gj9VGNkdlM1RZnSa3D6nh\/LUyNVwsTUzkbumtsnOLmNovvzAS7HzHUjX+V22PHEJCTueGoYoPnyvxAqKOJ8r3wl2du7rzE5gN+n\/1tOL7QwtyiQr6wWadHKmT\/w5j4hkHGQYNuKbGyZGDp6F\/augFpzFHGQi7cP311xPq9wrwhvxx9FNnAQyaJhZJn0MrV43oRScN70mTF\/yd5k86gl5Yt53qXnpXxMLvd2DOeiVoI73SxxF7xy10sPyByrTLHyc7R3lJmsNT1dSUWpkjNUlN4heuwG8BU0nBrzwZR50Fwlq0qukXrOJRdhxOpMl1OJGbE6lZ+GBRsUR3eQH167\/fPPIjzTjKlxd1sZA0WjcauB2l6uCuK5TIGL\/2jzpqnlVVMNJHm2Sbf6cQzH7jwe\/jIEMn0sdvkNZlqV6+GH39+vV2wDUmfQRge\/jhh4VjFhyxhgwZQnyBOudFCAf8jvz4X0j8cggHOH65LRIoy80JTRevuNPFQvpolCzVQ9pn\/f6vpvybUPUU6iks6XtL2n6kYZFmippS1dQo6XJZneKX3w8zldSClKeTNtd+SxpXdzRei4wOucnqKb8F1A8Lp\/d+pB+mTDkPL4pY1Pmxzgzcg1vr4K5P+rn5WQzSj4pX3KTPce25XFbXIFKmGlLBSb0DkgfpywuIU15W76ikP3v2bLrsssvopZdeoi5dutDGjRvpq1\/9Ko0aZd1FEMfiFhfmcjmxkX4hGqdTZmFJPz96pE57VFKAzn1VWRnVV1VRJR\/0xWRSVwxSs9p\/dpb01cNq9zMDHXLL17nnArrFZXBYDHysRT23KGL8\/e4z0J1Hfu2XMX7tHw0Fl\/R12+2XzpcMfQImovygAdfgNcuSvqqakcuSzwp0JP1x48bRH\/\/4R3E5S+fOnbQVvecAACAASURBVEX4hqE4a9kXSP\/Vd3faap0jB\/QgWPMEVe24hWb2i70z9r\/\/Sr\/7\/jF08rBD\/OZb3vsgH1UYIuLy47qOUe2cX\/sDgeGQuJjtd9K5x9X+Qt78ZR2kWyE88FTiwDpqw7P5g4xvoa9L5KBru994gbqOOIn6Tp1Hm+Z+Q1zQ4hWH3wkKX9LXuEoyaMA1Vs9wPH1VP++k04eEzw\/y43HS6W\/ZsoUQcrm5uZm++93v0pFHHln6pM86\/SknV+Spcvhyc92Aa0jPh8AMtp\/JZiEl\/ajfLn+0Qe+w1a03CCnolqnuVDCR0f6mVEq8CrP4udVdrPYXkvQZn5VZfNzuJoiCv077C0n6TPhdeldQz69dS9ser6V+1y6hD\/9wL+1a\/Scqn\/FH6tz9IO0u+pK+xlWSQQOuyZI+NxSWOIh19JWvfEVcgXjuueeKVxxWWZX0nax3DjzwQEH2iLT50Ucf0ZIlS0QUTnjg+\/ZTG7F4E8ai3okrtDJu3xre9wC6fMwgu5d+zlkdgfSt0MS5KJXF2v5HnSqlQso6pBkGK6edUJwLYxD83Uhf545cv77Lppmt7zfZpL+nZQdtmjue+k97LNbQyn7tSdJ7eE3X1dWJCJ67d++m448\/ns7GhUb7gnoHUvqcp9fTb6YeTVDvcFwe3esS1Tg+7NT1atNOz9g7TPrXn5mmccf2p\/KD99eeExs3bKBfPrKIvnHptyidDVugnVkjoVz+lC\/WUIYsV\/l1n1tB26I+xWx\/ofHpyOUvqPk2fTqwIju2n0cdVjt\/kPFtff052nLraXnWOzqhlXUbCxPNz7Y2Ua9xM+mDP9xDX5g4l5pnn0X7Dz+hzU1bfmUmVQL2a3fQ90ntZyySPoMBvf7XH1xN2z\/+VPx0+3nD8qT2IKCxamjC6IG06OWNrgHXmPRR9thDO9HYYZ20q9ncvJHm3Hozfes736ejjj1eO59uQrn8ZdOvpabOnWm\/jU00d+3bukV4pitm+wuNjym\/7VDrjm\/\/\/Vup3863qNOi70W+RMVrwqlXK7LNftDJnFQyDNoPv\/RJ7WespO8HgvxeXiCcLlphtc51Zwylu\/7UKGz9e3TtIqJsymcHG7a10DUP\/IZe+cNC6rxrK3Xe9X6QZhQt7TvPPScubIF3K+7kNY9BIC4EJpRvpwkDPqCywyup7xV1lOpjxacIel1iXO3xKyepZOjX7qDvk9rPdiN9FUD1XAD6fTzXnTWUJtb9neZcdJj4e\/oT\/6SFNUfkOXs1NTUR\/iX5uXH0aDskwuKXX05yU03bOhgC\/bu2Ur\/9W2m\/vmmb8A3pt\/8gljzpy5Y6L731IX1j3t+o5wH72Tp+nSGQTTPlQG1esXd0yjVpDAL7IgKQ9IPekVsMnJJKhnH3Pan9jEXSl6V0AOcnmccNrinPIGAQaItA1IPc9gitvIAW0ApaYXdmIk0UkUblx8lck71xneYBLGtWr15Np512Gq1Zs4a6d+9OW7dutZ21nn32WTr66KNFVk4Xx3wqedJnot\/84W7b2ubt91sc1TFxAGrKMAgYBPwRcLoj1z+XlUIl\/TDB1ZzqciNDEP56Wk8zybp0HU8t1QrHtzqqs3\/z8sYdPhze0fmPHGxNDbMgx9BxsuXXxSpIP6OUGUfeWCR9NIRVMKzSAflDxVPo2DtOIMCT7sYbbxSv2MECQZHk+Bscd4PzY\/BvuOEGuuSSS0RAJTxu6eMon+tDoCc8sqQSR\/kyLijvxRdfpDvuuIO6desmXsVVhyx1wbmF64irfLmcsBh5YR3HGLuVH9cYe5XD4+w0xnEQBJche+Tyb2EuUJG\/KzkQWQM1CAlfJnyZ+MfQGFvi9yJ9SPAcMA3zZcqUKSJEAiT8H\/3oR3TbbbcJTpg6dSpt375dVAGuOOOMM4ST1auvvhpbgLSSlvTjnFxRy8KEuP3222nu3LkEoofXHeJhXHfddXTzzTfTTTfdJKqQ02D7xxODJyJ+mzZtWpv0mFRxlI+PFDdZIUwr13\/xxRcL78A4ymccmZSPPfZYm5ALgREkLGA9ZswY6t27dyx9QDsxbnfddZcoE2MEvIKWjyiJTlhXV1fHMsZu5WMM4hhjt\/I5pIDTGEf9jnTyw3b\/X+teooGzX4zsnFVDNXnSvFq\/\/N5JvQMhDuNZW1tLM2fOFN8+h1GQQzCwpM\/qnbFjx9L9998vFgcj6euMegdIg9V28eLFdMEFF9C9995L8+fPF9IuS\/VY9THoF110Ef3whz8UpAJJH\/lAYmp63gVw18OWr0LHLuH8IcdRPoePPfTQQ4WuUpb05frD9kGOX+I1FcKWLy+8IHpehNUtvF\/56pgx1iCDOMbYrXx1LMOOsVf5umMc9VN10u8jBk\/QEAxukr4f6eOehvpsVCNV0lcFHAh5\/EDahyDCwdYM6RNFVu98+umn4jo6uIwDbHxI+Dspj\/yBg\/xBfHhA+ieeeKIdX5ulbZn0vdJz\/8KWL+Mjk5tKaFHKlyUd7gurd+T6w9bBZLthwwaxfZbVO3GUjzJ4XPDfWICdYpfrtp\/L48UD0l6cY6yWL49lHGPsVL7uGIf5HlWi7zvtCeox+sIwReXlcVJ7eJE+dPp470f6wFiW9GXByZB+bggikf7KlSvpueeeEwGGEKQI0eXefvtt+o\/\/+A8RcKi9H1mFwgQVJ+lHKZ+xcTpL4HdRypfDx4KQ3Ug\/Sh3IiwBT8m5IXkjRjyjly7stlMXqHVnyDVK+irXfnAjafrexjGuM1XJ0xzjsd1hM650oOn25f\/IZDX6XL0bBf0PdhnO0a6+9lhA3B+odCIAQXKDnv\/zyy0UedZcWBsOS0+kjuNCiRYvEwScmH8A8\/\/zzBQlgO4UgV+35qFtpP3WNk6TvtfWPWr6XVIh3UcuXD0B5HFRJPGodMmliB4Hy8GC3FFcfWCfO5YUt30nS9psTQfBxk+R1JHw53K\/bnHMqR2eM2\/MbdKvbjQwhzeMSGj7MhYS\/kBaKYpwOeJPYN6cFqGRuzvrXv\/4lpMcLL7yQfv3rXwspH7auP\/\/5z8Xq2adPn3YbEz5QlCVCt4NZVheopO+VPq7ynbaiTG5YOKO0X518qqQfVx9YVVJRUZGnMoujfCdS5p1EkPKBhRPWcY2xV\/lxjLFb+X5j3G4foE\/FXhIwS\/x8PwGsdibRpKR2xbNdJSfpo7cwA4SVzBFHHCEO2ZYuXUrQ8X\/7299uN72+ur1DO1nChZqDzbnU1VclfeSTy+L0cZUP0uK43Txz1Ds6VQk9aPs5vyqRx9UHFSM2qYyzfCeTzaDlu2ENCTuOMXYrH7uUOMbYq\/1uY5xklkwqGcaNWVL7GUmnL4O0Z88e2rVrl\/B269RJP9Jl3ECb8gwCBoFkI5BUMowbtaT2MxDpO0kcAMrLBTpuIE15BgGDQMdGIKlkGDeqSe2nFumz1QBAUW29vd7FDaIpzyBgEOj4CPgd5Ko9VMMw4L0sgA4YMEDcXsUOgnjPxgRsxgknP9xjCycsPsdj66fDDjtMOHTJj1xmWMQ7NOm72b8yGH7vw4Jm8hkEDAKlh4AbGSLOjlsYBvl35OcLyoGObLoK8+FnnnnGdrLUIf0ZM2YIh02knTdvHl199dV2uJIo6Hdo0pdXVtWGlQ\/ajIonyvTQy6uq1+KyJ9apnV3fr7rqKl8bZkx2doZBm\/G3m2OVTt0mTWkhEJX0ZZJXnQ3hcXvccccJy0KW3mFB5SXpG9L3mF9OMS\/i2AaV1pSOvzesQoMDCZMnjwUsk3grG3\/NuRJ1SV83XSHbaspONgJRSZ+l+5qaGhEFQOYgOcwCdgNQ5zDpI4YTBzhkhGTfFSPpJ3ve7FOt89smqucq+Bveq7Dz5yiDhx9+OO3YsUNMesQbwu9XXnmlwJF3DGo92MU98MADQl+KBx8ZS\/qqAMABrzhwHdJj94eHJX1snWWpXy6f7fz5o3Tbxaj1qqa0\/Ldc9l\/+8hfRD7QNURbxIOIi4hExHmYnUrxPKg7Sl1srq5cx\/meeeaat34f\/0LJly4ykLwGmdZBbvOlganJCwE9FokP6KBfEhkiC8FVgCQfB5kB8IHYOS+tEnDLps+s6h6FG+7gMdXGQ2466EdQOdTHJIz2MA9AOXhC4jaovhdpPrzwq6SN8LhYSbjvai\/IRyA2LWbF2TGaG5\/xf1PHV1elzvCH2ZJZj4sukL8dtwg7A7SDXqHdcZqUcfhiS4t133y0+YDgMGRVPYT\/lOEgfUj\/UQKr6RSZHXdJXQ\/rKW2wv0sc7SNscQpp3DgiJyzsTNdS0HANFjiuvhpTw2qWwpC9bePACw21ifAo7kqZ03vlB8FBJn8MwqCjBS5eDrfE79XyLy1IvScE4I60h\/RyqWpI+f2yyCzxLdjCTKvQFDvv6p5I00meS3rJli5DaZVL1In35fgOE7MDOQ1YdySFxWT2knleo3riqt7WXeseQfjK+JD91ZTJaGb0VSe2nFumrJpnq6bkx2Yw+QbxK8Js8OuqdOCV9hM+WJTUv3b+6YKEv2CEiNtPgwYOFakc+g9A9lJZ3nnL4CkP6hZ2LcZTuN5\/jqCMJZSS1n4b0kzA7fNqgY72j6s5ByvJ1cTqk76S3f+WVV9oc5OJ2L1U1w1K\/n6SvkjWrb5za7xYfifsif1Ssm8dBM+vtue1GvZOsSZ5UMowbpaT205B+3CNdwPJUPabsG6HeVQvzTtl6R4f02aMRQcJwTlNZWUkNDQ2O1jtyW0455RR6\/vnnhY4W5z2ITw71H1Qvffv2pf\/93\/\/Ns9OXD375ohH1Hlhd6x0ZA26T2nZD+gWclCGKZjK85ppraNSoUSFK6BhZNm7cSNOnT4\/tzt24eq1N+rIpnlPl8gXkcTXOlGMQMAiUHgJNTU2CDFetWlV6nVN6hEVtzpw5wlotKY8W6SelsaYdBgGDQGkgAOLHv1J\/QPZJInzgbUi\/1Ged6Z9BwCBgEJAQCEz6snkmX+oR1k5fvt8UbZIP+WQ7bNlMr5jxZsxMMQjsiwjsK1J4Mca2w0v6smkmAOOr8vDft99+u7hFi8OW+gHKh27qQVw6nbatL+DxOWLECM96MEFxYxeubQy6jcLFyAhNcOCBB4qbvjIZhBogqq\/3a73+e7UO\/Zz6KQtdx4IFROvXE1199XYbK\/3W6acsdD\/QElOH93jsS\/p2\/ZkZPmWH1+nLd4rCe5PvXcVECUL6kPBh642ASHhkL0z8NyxNeBeAu2LdLotG3ihmUS0tLdTc3Cwucd+0qUwQfjpNlA01E36kpZxyHWVlZbGUqRZSjDqGDiXq16+FGhqIOnI\/ioFVR66DvyccPg4cOLAg83VfKRQH1ffdd1\/HtN6RB4knBat0VNf9IAMKMpdJn3cOspfvBRdcIC5ehxMPHpgDsmewTPow\/zrvvPOof\/\/+2k3Ax7lp0ybq168fLV16INXWIjZ3q3Z+nYRyHWoYWJ38OmmKVUfPnofQvHmt9M1vxosR97FY\/eAx7+jjEaQf2Mnin98TRYjyK3tfe59ULAPr9OMcuDhJH+2aOHEiTZgwQbuJu3fvJjgVwTt00qQ0VVS00pw5W+z8S5f2oOnT+9Bzz70j3oV55Dq6du0apgjfPMWqo6oqTePHf0LXXPNBmzY1NaVCY8SFFasfPOYdfTyC9OOQQw6hnj17+s6lpBKVb8MTmCCpWCaK9OELEFa9g+3oscceG0jSh0PQ5s2bRR5IsZfNaKWTb26lSdIE6tatLJJ0K9dRKLVIHHXgPAP\/oL554YUULV9uLXIZIvrPVIqgiGv9Vit9tUuKHn44fwFEvpEjy0Sek08Otziirjj64fft76t1GEnfe2a4xdJnVTRUzkGfDk\/6hYiyKUv6ABR\/Bz3IjQIs61537y6nkWeXEeEANy3+n3Bj5vq9jZo1lGjWJCLlCk3t8e8o+t2qKovw+UwjXWkRfm22p5XAogFB9Yn2NOZ3H\/lwHtKo\/K4NUjZhR8HKr18duR9Rvic\/XDri+32W9AsVZVMlfXlhka16ZJNNNR5LlEnKH+c\/\/1lOZ68qo0kziXBdCIhuFs\/QvRLwpNrwh7sdhQBwUFtZaS1uIP5OEtljARzd0kLDUylqmp2ixolWGn6wYOBhqydI\/njkNDoffEfByq8vHbkfTt8T7wL9+o33mEPqA6FA93HKj7xeZWCeyXMNfXjwwQft0CC4uhP3KeBhk285jAjHqOL7cWVTdIQi+drXviau\/4TxCSR+hDzBxSzQTHDIESez9SjcpItXmHRa6p0kR9mMAmwe6Z9dRnV7KE+1M5SqKINlYGi6jXSrC3ZHIAB81CB9kDZ\/dFj4xlgbH\/EP\/fhmaystfaKHWBgnSTowsWDMJJo5yUqLA\/FZsyzJHx8j\/h4yJD+PE34dASudce\/I\/XD6nng8dfq+Z0\/bVJ1YgtAowCk\/snmVgbkm78TRB744XY0IDLI\/5phj6K9\/\/WveBeis3vnmN79JuFYR9+vioB+kDtNxJ9LHWSAeBA10ijQchZs0oAqdxJD+XpPNh5btpjsvHykuaqiknKhSRVWECxyoppHqxqQFIepIryBRTpdkAkA7a\/dKSQsgiVXlSNqNkH+yYwdN3zvRZy0gmpjFAhJYFVaIemvBnJldQPAh4h8WEt4JuH3QXF+SsdL9woBp\/\/45U+BCneMUCqtSkfRB0iBj9f4Flvbxvyz9y5L+OeecI8zF+TYtVu84kf7HH39MCE7IjyrtG9LX\/WoCposCLH84Ny77Jz12+dltSB9NGUpDKVM7idILZoqDTpZe3ZrJDl5M\/A8+2EJdu1q+AEkjACHhQxWTXedwpOGwOxddtbG6sZw2by4TixqctsQjRH8iCHlCIktbZL8wqxZDOtb5ey2ahSIyeayi1ME7Iiz+Tg58eM8L3KJF79Dxx\/dJ3Jj7fV5Rvie\/sov1Hn1g0lclfbUNUF3Pnj2bzj\/\/fHr66adJlvT50h\/4CqE8lup5J4HzR\/7NqW9JxVJb0k9qlM0owL7Z+ib9ecufaU3PNfRY2WPUSG1PIiHpQ+JPL5hlHXJWTqS0UGI4Pw2ZDFU11O7Vk+BUNCPy3d71u3T++T0TRwANIKlsV\/A\/XuewMlm+\/LLlyAYntoUriBZkz0IoS+68COBcBOcBeLAYIL2sFlIRjELIuoQSpQ4QOi92sioMdfNiz4vajh076ZFHUokbcz+conxPfmUX671M+qgT0jpL9fj7Zz\/7mdD5Q3ePx02nj3cIGz516lT7LmXc7oaIwl\/96lfzdPpIq0YaTiqWWqRfrMEKU09YYDOUoYW0cO+B7SxB4vg\/9R5Obs8CWiDSio9baPnr8tRAcrtrMwtoVrqW0lQploaGhjSNuuH6RHqy1u5V68yqtAR1tNVNykf\/sUBOfmey5b1cVkZDa4nSM6080A7xgtFpr7UT\/kCZQAzvsZlYULuX8DPeB+JRCFl37kSpgxcu1MXe22yxhHcgfCwGbM0ER79hw\/wdonTbLqeL0g+v+sJ+T2H6UOp5koplYNJ30pE5rXLFGtCwwNbutdGZRQuoonU5pVK\/pEZbJnVuOcg+Q+m9Mj90GjVisZjpkAe7ApBco7D\/tCyBJjQ2JlK9UwO7\/HRbCR89xQ5nPa0XvRbnGkTUr6Wf+O8lZWWiX9Dh4w0vGtxfIIRFwELKeir3bn4a9h64een1C0VkcZAlq3Z4tyKf2whhIOvnANUP+gH\/Dr+dTZRvpFBYhf2eovSlVPMmFctApM8mlexAlYTBCgMsiAyWORWt8+jR5sPo+D76uldYJ4IG6yhNmYY01S4kqs+aOXK5lVQnjoTxFOrjjEpkTGKTZmWobmZOXWUthtbuB88kmiQWN\/QDh9wjy0aKnY7uY+FlLQJ+Kp6kYsWkrlo4uWGAfgwfnqLTT0\/FGscp6pjrjFmY70mn3H0xTVKxDEz6tbW1wpxJN5pmoQc7KLAg5qpMDWXSlTS+5Xr6cTbgmu4hKwgMRCbUGpDu0pB1oc6xjNMze6X\/WYIorSdpRMaqB0ikC9CB+iqalbaIHQ8Orp3UV+jHT3b8hK488EpHPbW1D3I\/60DZMP3Doa6bE1fSsJLnrlgk9zqwVdYTTUxbuxzsZFS1GC9yzyZ0d+f3PQb9nvzK25ffJxXLQKSPAVRj4Lf3oLoBy3p4Jz19VaaWGtIz6amWFjosIOlzf6G6Bs2nodvOZKiybiGNoYk48hWKHdaPq0SmQ45BMW1Th2QyirKwLDEdw1oHzmbC+3ZSA2Vm1oizDD\/C9iJkPuzGTsDtXIT7BF04njF1e81FFxClsVPKhrJOMulbC6I1riB77HXQFXms8R6LAebFjp076ZGUOcgNOpfbK70Jw+CBvHoSzkkLeUeufI6gXqLiRPogfKgpQEJOqoisA6kgfQ6trCvpc39BpDioBOnDHp111SACWfkhExlUI5CooTKJ85Hr2FS2iapqG6y+z0xT7V4KmiVaOiZ7uAxtvPUeZI3\/dTqbUNvnRshYxGqoRiwawF31dUA5+B0P+s27JHG4ix9xf8FESycO+\/bDDkumeSuaCmvUrHWqvQA4Kbr4LGNdaysN04hsGWYueC2QHFLDz7zYqd6kSqdhMIojzz4bhoHBk+PpI\/xxMR65TtSnxu13mqSWUxWkMEuzbruRw6EIEhusSIjo4hicaJy8WWVcbPv28nJ6rGyqoPsgOnEdjPPMKcteFiam3Hv8V1pQFCRQUHPGJvwgC5BKMp2okyiNFw+Yu\/JiK5u+YmHBooDFAb8vFIuDpePHAgnHsDSf9orxbU6keSsWKUvSr6WZ2QUUiEK4mEgT86y5hE4\/laLb9l7ME+\/ynpsNeQv9Jstvgh+x0HdKhzpILgXSRx9MGAZ35gik3oEjw\/33309Tpkwpmk4fA6hzicqPH\/0xvT76dUE8lgpllrDIgRoDumRhf4K9OFQfDUSNk+LTt6MOt5ghNun\/s5weO2uJoLo9wo0pvkcmgDvL7qRZe81GKdNo20vWCwesBmF2CqL3U+U4tcxJTYXy2ISV8wB\/Xgyshcb6PzZ5xYLHpx+8HE1qgP+D5djUu3dy7ds7UQOls\/ZI2LVgQeP+yWNqHXpnaGTZlwIceQebD\/J4nH3ny9QwZCE1IniUWJyqxDyf1FBHdZPcDHGd63MifYxXNpySbyOdagsQesfVbNirDA4Vwo1j5ykYnJgwDG2HLBDpI7sMqO8MiCEB6uMbulCc0yUqY+deRdvvWU2tFRW2Nv2yN2fQw+JvS9LH85+PWCGC553cKkgaH06QiyjCdIfruPnpClp65ROC9Ne2rA1FvG71y\/24oOcF1No6kFbu9wil64han26ldUoo5Cj9wIUzYS8fqU5V07DUMBpDY4SE\/FTLU2KJ+GUqRVCFwJlpypQUvfrqzoLat4cd825l3Whe6zyqaK2gKakpNKN1hqOq7pHUI+I99jPzHq4oyKUzPOatFa10ZI9zRZuaUivF0F7fcj29kHqBhtAQerj1YfFblNDKeQEIfSaPkzgTIPSOqzjkVYbsBMgcZcIwuA9UINKXo2CqRRZKp+9H+itXNtEpJ\/+IDrnvA\/r6wefRGWeeSdf9ow\/ddfgWERlSfnDRB\/6NHm39XsxLO557bijdfFnF3kOATjRnyxy6cOeFYbjXMQ\/3Y3f5bqoeVk0X7nySKm4+ae9VbT1p1KgWWry4OXJdcWDVlGqi6X2m08bURrpwx4V0zQfX5LXrrc\/eourhp4uLbC68cGfkNjsVELYf9\/32Q7rvP46hR5sfpdEtowl9AdG61VGVrqJNZWfRhT98gOZckbuYR6dTKJvLX9pjKd3X8z4xXzBv+EE\/LoGjXNezKJX6BT33The675APxXnNNR+cR0vff4WmHz+O3m58m5b26EGTU6nQl6h0REnfhGGIifR1JmzcaXTVO1EvUQl6kKvbT7604\/XX+9P5Rx9ClekMvdl6ipAY8cAmnp2fuMzhqeG0vHW5UFX9MvVLkVYOBKfWjTpe2foKLeu\/jO7u8QRVNjTSU6Nb6M47y+ikk6JdbMJ1xXn5iBthYicwuWkyff2QrxcsfAH3Y\/78CnriiR60dm2+YOA2rv+5ciXNPr1ajMvJrSd7Dj\/qeGj3Q3Rj\/xvpsqc\/pYcDXCwDdREO\/HMPlBfQS9ZSRevTtK51nQhlgYtumn6REqE+ZjWk6fqL2\/YDDmLj17bQqi+V8VURvtO2VHT6TProsAnDkD\/sgSR93xlTgARhDnJ1m1FsE8FuS8qoEheycDyf7EGoerAL1QefTWBRwP\/hENRNF4+DUraQETYyC6DL1UVBL10xsHq69WkC8aezvdVrWbBU3I+HHiqnxx4ry\/MZkL1q1VL9fAzk9FzHF4d+kSpr66l+pr5enc1fbRuwzCSq2xvJdMUkyPC1VFljnVNhsfrzn7cQLgCCE5jTU7PXJBaH5sJ4QHM+lALpB5sRhUudVCwDkX57qHcwJLLJZiEuUSlWBMyRS8ook41GCYkOjxuR4+NfQSvEwStbv8BKxMncE+qGxc2Laf6gt2h26nQR+hix7eN8ikH61gFoJa0q60azqN4nMEa43uUdrG8qEwf6bPgCHwI4jzmFT8A7LApO0TXVlnAdpw46lZoeOZ32TNL3YLZ9HjJ7BLlXSn4MTguL19zl6y+9gtypbU8qUYUb7fbNlVQsA5G+G4SwrkH40TD3SEYdlijAFovI2Bdg5J1llJlpmSt6+662RYWlf6dIoLb0Wl5Od24qE5JhkA9dZwyKhRWint5UfhOtTJ1OdcKnId6H+3HJ3tvSVp1VJqxFoDwR9wIstCy7cMgP4pefKuGJZ6X1e3JmujfSY0u6UWNlndY9DChXtKNqr+9CI1HGYxwLNR5Rvic\/XPa190nFMhbSd7o1plgDHAXYQn04bhLZ1CVltGBSONL3wjPPF+BlS3\/rZkIadlyKiRUOpKHXniQ8GuKlfVsK\/2yQsBDC4gs+bxDRQREKO0OZmso29yawQ18Q0l9Wvoymlk2l+kwjVercvpMdHMQo8nOsKtR4za0HIwAADGdJREFURPmews6tUs2XVCxjIX35sLXYMXmiAFuoD8eV9KeW0YK6nFenzmQH2cD62ksrbEuv5Xul1xIgfagsQJbQR+\/JWrrDTltfM+6OrK1vX7SI0jMXiNMDjpdjxdMZSpkFlid3epLldY1rI63A2sEkffSj26aRNCnt7BmuM\/5uaQo1d6N8T1H60955TRgGlxHw0umruvZiDWKUSVqoD8eN9GFNM2si0ax07nIRP5xgI80erH4E8EWEgazdGxKCo735FR7gfdGxevllmlVZI1Q8aZokQjkjPLN8F2qA5oukwBFmvFv+vIVO\/WyC2A1xrCBeVKBThyqtYcHe2BAIUiS5JU2iWaI9fo+M1cg7l4j4RoV0yIvT8izK9+SHS0d8v8+HYUjioEWZpMUmMp2PExEJ5EtJgLkc8wV\/q1Kvl0VKXGNWbKyWLCmjmhVw1bUUK+kFeyiTva0L5xV+B+FqvxEIYoGIEVRJi945jXDI6nYnAh+m4j0cnKDMr8msoFnpMXstr\/z3GzJWU8vKaAFV0R7lNACH89hNsFdv0HhMhRqPKN9TXHMtajkmDIM3grGod6IOUpT8USZpoT4cN0nfifQ5OBcOdyFX8gM9M+uPrRBpuduo4IGIhw8hbZ3+jeW0alW+GWIUbIP0I4561PEQgcNE3IzM3qAHlfbhqgiPLZQy2DH5B7Djew5QBki2tbVVOD85BYfL4e8fKtqtz3mk\/3IZLai07g+WH+uyHev\/8ID0YZ2F373Mc7mMQs1d5zAMuXb6jbOTP4nqh+JVhps\/ilcZfPMdl2vCMMRE+ohhUVNTQ+eeey5xTAv8zXdGzp8\/v2jxeOQudXTSR1+yd4kLEmejEStQWu5vEZxsb2J2OZdjuds6\/UvKxY1NOmaFfh+v+r5QJOO3sPDl67I1kggLkEFApTFEkxaKo16QP9QyeNjvAX+NEUHRaijTUEmNldZlMDd8cgOtWdOT6iv9VTVBcUJ6FSt5p8YKI7RrPS0URI9HqJTEEoeoUf7mqoUaD6fviS\/W0cHCSY2F4Hy6j5sazKsMdceGPpgwDO6Ia0v6MMvE7e\/jx48neBwiBs4ll1wizDTbM8Z+KZC+FXjYcqLhB\/If4s0jWBCbD7KjjRchT51qXdMX91MokvEjfad+CNkYl7ojrhh4e1JV9i7e3P3FYgclNgmIOYlFYqbwXWiPfvACjgNhtKkSl+9U5q6ptIKJNVCVCOY2yWqrj01vofpRKpK+CcMQkfRVk0xI\/XKIY2Oy6Q5wlI8TpIbHT3KPUofu4pDEOtimPd2YoUwa19ekKdOpnibV73WySlshmyvTlu0739bVHv3g29ZE4L+9DmB7IyfYl6hbUr51XIz1C5ajsxqI\/Jx4C9WPKEKU7lwqdDpZ0kddJgxDPuJakr5K6gDxxRdfpDvuuENEXDSkr0\/6QTw7QfqwMHGzWGGPS4QIePbZZF6+HuQDD0Nk6DtIEgcgUI\/BwgceyWBRLJaqv0KYOoL0AWmd6hB8D6\/e7AIEW3y5fdaNB5bXr45zXaH6UQqkH3S8CpU+qVhqkT5AYfVOdXU1TZ48mS6++GKh6uF3+F\/o+oM8qlpINgnF2QEvKnIYBp2bs3TbUKgPR67fqQ5YVoLIQUgLFxJNnGhJfvzgw8ffIAbcyiWTPt6B6Dg9dN79+rWIyJTHH69\/wbsuRpyuvbDSaad14GudgTRCys\/qy5y8kturH+plO0I1lcntQPjMBh3Rca4rVD+SSlQ68yBpaZKKpTbpy4T8ve99TxA8\/wa9flDCxyLy0EMPEZclLyxjx461zwxGjBhB06ZNo5tuukmMqc7NWbqDX6gPx4\/0QUqQ+qB75ztj+bpF5AUhgNSRTpYGQRKcHioLkBoWBFwzGPbax6RjpdM+VvMIu3sfl9n2GnOW4nmx5gvqMX4YR7bigq6nURECnDAoVD+SSlQ68yBpaZKKpTbpxwkoOzysWIErTawdAi8g+G\/5cBgxfXRuzrrmmmvovPPOo\/79+2s3FR9O2As1dCtxq4OlUERIxMUhiJooS\/uwwsHz1FMtYkdQXW3dBXDZZa10882tQkrk9O3ZD10cdNJF6Qd2RMDJL\/xElDp0+oA0YevAJT8La4nWrWv1rSpoHUEvUcH3NGrUKN92mATuCGzcuJGmT59O7eW46taywKQvq1rkQsNcogIyl0mfJXrcv8vnBhdccAH9+te\/FqoePE43Z1166aXi3cSJE2nChAna8xCDgrJRx8CBA7XzBUmoU8eppw6ia67ZnndxyPTpfWjp0h709tuw4CdB+PIFMHIbdOoI0mantKYOfQTDYHXffYcEuvQmaB2HHHKI1iUqTU1NgqhWrVql32GT0hUBLJy466Mie4tfEqAKRPqqNB61A3GS\/oUXXkjHHXcc4To\/3QcT\/JZbbqErrriCjj32WN1sgdLp1PHLX1oSfNhHp46wZXM+U4c+gmGwWrFir33+rDRNmpShiRP9b6QNWgfMrfFP50HZ+Gee6AiA7JNE+OhRYNKvra2lmTNnajtisU3\/k08+SQMGDKC6ujqCJI9HJX0cEAdV7xjJJPrENCUkA4GdOy+ksrJVlErFT7hQ1+CfeQwCgUgfcMXpiCWTPi8CkEaCHOQin5FMzEQ2CHgjkESJ04xZ+yAQivRvvPHGNq2NqtNHgU4WQvhdPkdI2qFI+wybqdUgYBAwCIRDIBDpy\/fVsoomXLUml0HAIGAQMAi0BwKBSB\/6+fvvv5+mTJmirdNvj06xGop3JPIuxMvRS40ppO4yVMcw2b07bB3ymQfqk\/0W4uwHj4PqTR1HH1A2B+RDAD7ZsS7OPshtjYKTF+Ze8yMIVm51xDnefmUx9rL3fHt9j6be5CAQiPSZBGFfH9QZq5hdVmMD4ewAZHTdddfttXG\/2dHRS1YtsQpJ3tmg\/bJjWFx1yGck3Ab2do6rDsaeiRmWSjCBxVmI3Kc4cMIOUL4zOa4+oByM31133UW9e\/cWXuF86B+0juXLl1Mmk8nzDwHm8DZ3cwSMqw6MhVPd8G6Pqw72lFfHGyFTzGMQCET6XjdnhdHpFwt+SG+LFy8W9vj33nsvIQw0PgCOFIq2Ywdz0UUX0Q9\/+EObTORrIOX0ThfAh61DxUCOZqq+i1IHpMLZs2fToYceSqtXr7ZDXMh1hC1fjl\/uN6Zh65AXYJC+7NMRFCd1\/BjzIUOGeDoCBsHKrQ4mZC4rynh71aEz3n5jZd6XJgKBSL+jQiB\/1CB\/N0cv1Q+BCcotvYxH2DrkMvzOTKLUgR0FHhAbY6BKfmHLZ5w2bNhAa9asaaPeiQsnHh+U53V\/g24\/UI6M+datW21s8E51BAzTD7UO+SwsrvF2qkNnvDvq92zaHQ2Bkid9WX3iR+JhST9KHTx8TucJ8tBGqQPb\/IULF9KMGTMEKTuRfpTykXfJkiV5O6gTTzzRDsjH\/YhSh7zrQnmyeicsTirmfvMjTD\/cxjXO8VbL0hnvaLRhcndkBAKRfkdT76hbZz91jRPpe8X9wcBHrcNLEnRTAQTthxpPHOXKh61R+yCTJXYPqv9FHDip\/iFR63CSsv1wDdoPN0leV8JnVZBXu5zK8hvvjkxYpu3REQhE+m7VyQd30ZsUTwlObfI6mGXylSVIv\/Rx1eHl5RxHHTKiTgSNoHayftiv3+riKKeHE5CqFomjD07EJ+8mgtQBPJww9+t3XHXENd5u\/fAa73i+LlNKR0YgFtJvz0tUnMB3CgrHki3UGxygTXX0coot5OYYFlcdHGJa7gebhsZVhxsJyFhwmjhwks0p4+yDm8lm0Dp0MZfnR1x1wHIHIcXjGG+vfnD56iLfkcnKtD0eBGIhfVkK69WrVzwtM6UYBAwCBgGDQOwIBCJ9L52+CY8Q+9iYAg0CBgGDQOwIBCL92Gs3BRoEDAIGAYNAUREIRPqqxYTqQVrUlpvKDAIGAYOAQSAwAtqk72Qih9rYRtjJLjtwa0wGg4BBwCBgECgoAlqk72edo8YMKWiLTeEGAYOAQcAgEBqBWEjfb1EI3TqT0SBgEDAIGARiRUCL9FGjV2AoNVxvrC00hRkEDAIGAYNAbAhokz6HaUXNfM+t02+xtcwUZBAwCBgEDAKxI6BN+lyzGtdDvVgk9haaAg0CBgGDgEEgNgQCk35sNZuCDAIGAYOAQaDoCBjSLzrkpkKDgEHAINB+CBjSbz\/sTc0GAYOAQaDoCBjSLzrkpkKDgEHAINB+CBjSbz\/sTc0GAYOAQaDoCPx\/GLdc0LVu+wUAAAAASUVORK5CYII=","height":140,"width":232}}
%---
%[output:0207f450]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsB0_S\nTaille de la serie    : 193\nStatistique T_max     : 2.3029\np-valeur (bootstrap)  : 0.9610\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:69a175b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:7eb82644]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsGtsi0_S\nTaille de la serie    : 193\nStatistique T_max     : 2.3533\np-valeur (bootstrap)  : 0.9520\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BbsRtsi0_S\nTaille de la serie    : 193\nStatistique T_max     : 3.0378\np-valeur (bootstrap)  : 0.8890\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:5ec3f60e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_STN","value":"NaN"}}
%---
%[output:2a69195c]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_mean","value":"0.8380"}}
%---
%[output:6851978f]
%   data: {"dataType":"textualVariable","outputData":{"name":"ratio_median","value":"0.9620"}}
%---
%[output:3036bb08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:85ce0528]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAPoAAACXCAYAAAA8hka5AAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQu0V0XVHxCRFFBAiacvVqKlCVKJr0XLfOBSXKJhal+gEWBq5kpCEFdqiaBivpBAJSRfFfmkVBALC1T6BCVXaZYJgmBiaoKoRPCt3\/Dt675zZ87MOWfOY+6dWYsF\/M\/MnD37zG\/2Y\/bsabVt27ZtIpbIgciBZs2BVhHozfr7xsFFDkgORKDHiRA50AI4EIHejD7yli1bxPz588WJJ56YalTvvPOO+MlPfiKOOOIIgT54gWXXqlUrof7dpk0bMW\/ePHH99dcb3\/X666+Lv\/3tb2Ljxo2iY8eO4oknnhDdu3cXGzZsEF\/4whdkv7vuuqv40pe+lEjvf\/7zH\/HSSy+J3\/\/+92KnnXYSe+yxh9hll13Eu+++KzZt2iTw\/DOf+Yz48pe\/7DRu3\/2pL0X\/y5cvF7Nnz5aPzjnnHDnGX\/\/615K\/oP3II48UO++8sxO9PipFoPvgYk362Lx5swTs0KFDxZ577ulMFYB+8803iyuuuMK5DSpeeOGFsp2p\/OY3vxHPPPOM2Lp1qwTojjvuKCf4Y489JgGO93bu3FlcfPHFAgsHgK8rr776qlzA\/vvf\/4oDDjhAHHbYYRIs\/\/znP8Ubb7zR8Oy73\/2uaN++vbEf6lvtb+DAgbJd1v5Uml9++WXxy1\/+Ui5CKOj761\/\/upg1a5bA4odFD7zefffdU\/E7T+UI9Dzcq0lbAACT6l\/\/+pf4xS9+Ibp06SKBc+aZZ0qA2QoAd8stt4jLL7\/cVrXRcwDrpptuMra54YYb5PuPPvpoSdeBBx4oF6Fbb71VvPfee+Lkk08WDz74oOjXr58EMaTy5z\/\/+Ub9rVu3Ttx+++3ii1\/8ojjhhBPks48++kgsXbpUDBo0SP4fUvLFF18Ujz\/+uNRm1D54h7r+1AGk6Y+3hSR\/7rnnxNy5c8XnPvc5KcmhCQH0+H3YsGHy\/7\/61a\/EhAkTItBTzbZYWUqJP\/\/5z2LFihWidevWckL99re\/FQcddJA4+OCDrWAvCugzZ84U7dq1E1\/72tfEtGnTpJZxyimniNtuu03ssMMO4pvf\/KZcKAB6gKt3797iO9\/5TqMvCnX3tddekxIR0h8FC9rUqVPF5MmTG9X90Y9+JPr06SPOOuss46zQ9Weq7NIfbwtJDpB369ZNLmKf\/vSn5eN\/\/\/vf8ndoDAD\/Aw88IL7xjW9IraasEiV6WZwu6D2QYpDmAwYMED169JAgwETDbwADpCQke1KBrXvjjTeKK6+8MhWVJokOE6Jt27Zi0aJFYvHixWLixIlSvcZ7Tj31VKnCQm2FnX7vvfeK\/fffXy4Cd999dxNTYPz48WLUqFHyOVR\/lA8++ED6Bn7wgx80ohfvWrhwYaIJouvPNGiX\/tCWJPn9998vYAZ89atf1XYJSf7ss8\/K74HvAhW+rBKBXhanC3jPX\/\/6V\/HII4+IESNGSMn55JNPShAMGTJE\/PGPf5SSEgvAvvvum\/j2V155RQKObHQsFiZ7GaonTAMUHdA\/\/vhj8Y9\/\/ENKK9jCsMd\/+MMfSucU3gNt484775TSGwsTSTf0e8cdd4hrr722Ea14x9VXXy1tcqIL4\/rZz34mxo4dK0FGTr2\/\/\/3v8ndIYlPR9afWhVYE+l36Q1uS5HA0YnwmSU2Sfc2aNeIrX\/mKOOqoowqYFfouawN0DB5\/YtFzQPVMAxjXXXednDBQCQF6qMWQKvDmAtyDBw82spOADAkDNR8qPnnrYTfrnGN4J1Rs2NkmoENF79Chg5ToADpA+a1vfUsuDlgA+vfvL4HRqVMnaZfDWQc1d8mSJXIcUGl54Q4\/ogs2OswU7BKgDyxy+DdMEJgLWFhMRdefuqi59Edtfve738nFDP4C192ORx99VDz99NPSsQiTpIxSC6AD4N\/\/\/velgyUWPQcADl4AoJ\/+9Kfi3HPPFTNmzBB9+\/YVxx9\/vHSqYdIDUKTq6nqEjYwC9R6TFuokJGyaopPokMigA+YDQIetNdjoWEgg7T\/1qU9JPwIKtI\/169eLrl27yq0n0HTaaac1IuF73\/ueuOSSS6SqTzSjzx\/\/+MdNVPcXXnhBLnRJEl3Xn2nMSf0RLZi3UMOPPfZYsc8++zixb\/Xq1eJPf\/qTdERCi0EBBnr16uXUPkulWgAdUgUOFHxkfChIqp49e2YZT2VtsEjBsVQU7apEh4MKfwAa8AqTjcADNZcX2gPnv3GJDqkETzW82ugvj0SHegqvOrzk2AaDDwEg32uvvaSU1xW0waKABQJ1eYHdDlohLXfbbTf5SHXGYXz4gwUPY4cpYyq6\/kx1p0+fbuyP+IfYAIwVW4TwPbgUaDhYqBBjANMKBRiAfV9UqRXQ4UWF5MKH2nvvvYsacyH9rly5UsyZM6dw2gEcSEY4tgB+2OKQLrALocZCisNupQJw2VRKqNvQEC699FIJdGz9oB1JXuoL++FwtJG32+SMg1ZB\/UFdh9oO+kgKqh\/grbfekvv\/0E7IU011sG2GwBz4HbCLgIK+YedfdNFF8v8wJ6AZYBsOi8Whhx5q\/Ma6\/tTK6A+LCS0+Sf2hLd4N\/wYWOBTwCYsXBR\/BDML4iZ\/nn3++rI\/AHyotCuhQRbHyw6kBOymkArsRe7RF037VVVeJyy67TECtBKAXLFggvdqmAg0D6iXK2rVrteo5JirA7svrjv6gikO6AqBQ15MK7FXYuSaVGx572OTYb8f8gIqMhQhgQpQdvON\/+ctf5Jw57rjjpDc\/qaj9Qahgv\/\/999+Xf9L2pwId\/8c3gjqPAqmP70bAjkCPQLeubQR0ChCBCom9Z11gDKQvzAkCOgJidGD2vY+O9yJkFZFxZ5xxhtXrDyDA9AE4dAVgBjih6sLJ+NnPflYuqHD2QfIuW7ZMggjAhyNw9OjR0ploKr770wEdvIbqj3LeeedJvjcroGPSjBw5UjpQyOZAVBTUQRRs46i2CNnoUaJbcS6j1+A0A8DJToTKriuwWxEiCmcYwAQAYh9ZLb6BTv3fdddd0smEd5q261AXKjIAmhTognro680335TeeajHADO0BvgXYL\/j+R\/+8AepWWHrCrZzUvHVnw7oXJVXJXizkOjXXHON\/HAEaKzC2AdFgAP2Ue+77z4xZcqURupcWUAvcgsP9tjbb78tJ1wRZgc546CyA7AuBeBCUAaCSEAbpD4BnQOvKKA\/\/\/zz4qGHHhJjxoyRXnjV5qcxTJo0STruECqbpwD8+MaIvIOtr0bZpe07qT\/OPxXo0DIQzEM2O4CNLT+KPwge6AAsbDOoRtiugeSGNIcNBnB\/+OGH0jsJpw95icH8MoAe+haeur3mMmmxIECqQvJRfpGzzz5bNuWOsaJOr8GHgG0kaBNwpGGh4afg8G\/YxFiEEE3m65DHqlWr5P48VHksMnmLrj\/OP9APkwo2ORZ8bP9hceWn14455hh5uAXfAvyAmcLHG4wzDpMFdsi4ceNkdBMHOjzSUOV1an1ZQKfFpKjtr7yTydbedpTT1l59ziUSPPiQglDzk1Rs3gdACrVZjTdX62CvHAs7iq5v9AMHGqS5r2Ob6BPefhyqgZMxb9H1x8cCCQ5BhnnOj\/TiMAvK6aef3miBg\/MPmhVJeMJAENtrUNkRHYR9UAwiC9Ch0mArAyd\/bJ7atB8Pp4eGDx+u9RGk7SvWjxzwzYEgJDpJamyB8AI7HWqPq+qOtgjjhNfU5XhlGmaDNkRF6ZyBafqJdSMH0nCA\/ELYGtQVOBfxByUIiU6DgB3OJXoaZxwcFlBrcG7Zt0MLKyaORVYFdPgqsAfs8jHBM2wnXXDBBY3mBt+9gO0Jc4j4jaASaFOIzqLjnGhMmpZuB4T6UCcgmTn4Hao5vPZUQBsCg7B3r9O6yBmLnQHYqKov5qmnnpJ0q8X32FS66X2msSXRTW11deg3qqN+A3ofDtvo9vcffvhhuSNV9Lz0HhmnAh0MoI+o+\/hkn2BrpcjtNWJ40QxNs9qb6uqArgKMAAyNCQVgxBgJSHwB4DsgHKTqIoB+yNeCbTwAGR5xRCoCsMRDBMGoOyf0Hen9Kr00B3SLi4+xgW5ExuHQCsYOfxHGwBc909hgY5vo5guESx1VMtvmne25j\/mEPrwDPQthZXjdfTCUdhVwzhpRZiQ1uOlCKzr4gMmGYA4AA5oKJDr5MCCBUQiEvA8setiOUiU65y1pCPA9IOwTQOSTHVuZcG7hGCh8JzpNQqdlqIsM1QHdiInHe0wSnfcHsCHNFM6So2BbEPSYJLrPsXFAc6CbxganGGlQSdoUaWR8bNS\/qZ1t3tmeZ8GTrk0EOuPKkFufF6vf\/Uj07tROzDu\/fxN+4aNQHACXGpAkBCaSTt\/+9rfl7gNJFQIBJDDfhSDJw\/tAXaqj+2j8ObzKHOiqJNNJbdKydO\/QgYE0BvydpLqrQOfagCrxTRM479hIc9Cp7kljs2krtrHx78DHZgOy7XkEekoOuDC031XPiNff+Ujs2bmdeOGyw7RAB1DJZsXHBcBx1JAkNBpB+mEbDyfxoEpCBeZA55IBYEB\/eE4qskk6oG8AF9oBp0En0cl+1gFd7YPsTNCN\/V3sQZM2oUp9DnQ4mpAaibQb0Jck9bhpoft8PsZG\/RLd0CAQwJU0Nr74kk8J++JYqOlADexraDPQjFSJzjUp1W9hm3e25ymnubF6lOgaiS5VTQPQbRKdulM\/fhqJrgOEzvdBkpkkrq4dB7qpDz47kmx0m0Tn79dJfhPQfYzNZOdzk8U0NmS\/IR6a1P6kseEZFwBRohvWm1BtdJ19jSHCJkfgEM4+qxLdxUZHHWSO4TY69xYTG\/F+3p\/N627qw3T2AO8pyusOcJNqD2eYGueedmy06NGZijQ7CtxxiX5MDluTZz5pR8UmsW3PC5Xo6sBp8uo8rT4ICQnoppXbBx9aWh9w6mFR5A6z5sYDG5Btz33xo4nqTiDHVhjtd9JveGkRYI9A9\/U5Yz9144ANyLbnvsbTBOgmG8X0uw9CQgG6j7HGPloWB2xAtj33xS2tM45vI3HvLff2+iIA\/ZQJdKQ\/sqUG8jm22FfL5gAOCiHxh8nurwzoprh1\/rl0Tp88n7MMoId+TDUPf2PbajkAwYLtVl2W18qAXgVLygA6xlVk4gnkjEPSQpwxTkpjVAV\/Xd4ZMv0m2q+Zv1IGQKFMO2N\/FzYUUgcAN6VyjkCPySELmXSmTstKblnEoFTabRGORdCQtc9Kgc5PEhWpslPfZUn0rB\/DpV3IQMH4QqZfpd0W4ejyPcuqUxnQYaNTuid+IEINT\/TJiAh0n9zM1ldzAjpJdHBCF+GYjUPFtKoU6HQwArmtdAcwfAc4RKAXM4nS9NqcgJ5m3FXXrQzoFBxz+OGHi0MOOaQhBxzigX\/+8583SWzgg1ER6D64mK+PCPR8\/MvaujKgg2AeHAOpnpSTnQaIQwV0kkndfot53bNOg\/LaRaCXx2v+pkqBnnbI6ukj2PMoCKFNk0qqyAwzaceUtn7IQMnrjKvayx0y74MCugoK7rirS173tMBNWz\/kyZYX6FV7uUPmfaVA5+BE9lQcITTle+OAIPUdd0xRkkKeMaTKvO5pgZu2fsiTLS\/QuZcb2Xl0WXqKlPoh874yoHP7HBOAttrwb7payeZ1B\/GQ6gC76rlX72VDvzTYIvO6pwVu2vqYbEjbi6uAfOekT0tLlvq+6B849XmZpadHxzbiufGf3H1m+j0LrWobX7T7oCVtH2XdN6A9vUbgRkIAyqiC8FFXoEOy040t8NbXIa972g+Qtj7uLMcdXNBmfOekT0tLlvq+6B\/9wJti3YYtkoR5I3o1kGL6nSrQ8+4d2ojbTu2WOAS1ri\/as\/Atbxukgca5\/KKzExtPr3F1nTKAmIjh2UKQV4ur\/nyBsF2yWGRe97wfxNYePMAtn8je6jsnve3dPp7noX\/ozBcb1PUHxxzUiJykZ7wiSXzk63t2bNPEnEl189Dug3d5+ggur7vL9lqVed3zfAyXtiHbiRhfHvpVZxy3x2GvQ5VHOaLPbtrsuniWJpoNdZe8+l5Dn3NHHiCvSsad6aEtspXZ6KZJjZs5ESUH+7Nnz57GK3BdQKHWiQEzWbjmt00eoKsg5cCHc45Aacqum2Uk\/B3QACLQk7nolAUWtuf06dNlwgZMCKjjSMzvyxaNQM8y1f22yQN0lRIV+CZpnccTz\/uMQLfPBSegQy3ftGmTzDiKsmDBAnnjKSS7jxKB7oOL+frwCXRXSnztv1dBu+sYbfVqpbpj2+iWW24RRx11lLznefHixQIXu8Pm7tu3r20s1ucR6FYWFV6hCrCkscuTGFAF7b4+SK2AvmXLFrFhwwYJcl7atm0r2rdvn3vMEei5WZi7Ax1Y8qjWuQlK0UEEup1ZTqq7vZt8NSLQ8\/HPR2sdWHyp1ib6+EKCOkn33kWJnu8rR6Dn419D65ClCgaRJNHxHAkckiR8FunPFxK8I+neuwj0fBO1AehVZH8l0purRM8y+fN9zuytXRYqAibeou6J82fYRjPdSMsp1MXI06KSZiQutKfpr8y6ldjoPOkE3dZZxqCbK9CLVn19fhsbWHiQCt6r7okTaCk4xueeuW2cNtpt7at8XgnQMWCeMw7hrGWU5gp0X15ll2+gag9ptYkksOhAbpK8ZY6Z+BKBbp8h0Ua388ipRtWTTdUe0moTSfRztVwnzZ0YVGClqnmfZ2iVSfQ8RGdt21wlelZ+ZGnnGo1m6tsm0XnMeplquQsvItDtXIoS3c4jpxohTTadWo\/fXlu\/UeyzR3vjwRPV4QbGuDjdnBiYo1JIvFeHWbpEj173HDMt5+mvfG9O31rnISeJnSStVYdbXdT4CHT7HGgk0aPX3c6wLKpv9l7NLdM623hPOg85Ej4A7Lg3znbpATcTTE65IsZcF977HFvpEp2Ij173bJ+xbKmS1tmmG1WWo555FphsnLW3Kpv3dorca1QGdHcS\/dWMzjg3XupCRrNIVJ1Edz3q6WOBcRute60IdDuvUjvjcMDl5ZdfFn369GmUBJHAileqmWTiBQ72D4EaJmmZNRhFjTybd\/72FE1q6Cl+A9BRf\/2HyQ62KvbJbdyLQLdxSAhnoOP02uOPPy4effRRMWTIEHHiiSc2JJ5Q1X0Am65vQr45SippyxnX0i9w0ElLHqwCRxkVmx3NAY1\/cyebLvSUbPS1729JrEuLhX1qlVcjAt3OayvQAc67775bdOzYUXTt2lXst99+YuDAgYk989tZkO7ZNQtsSwe6Tlqq0tcF4PRxdG1VrYHqID0zwA6Jzs0B1UOf5v326eenRgS6nY9WoCMVLVJJnXbaaTIvV+vWra1A55c2pLnAIeZ1b\/rBKIsqnvA86VRz2KyXpLccIDVlYOVt1fzq6H\/pqo3bVfpubcT9ow5sZJLx95tosE+zYmsA6KHm1K8sr7vuk0Bth2TGfen9+\/cXw4cPl\/nLdQW2OuWCRyLJNEBHf6eccooYPXq0t3x0xU6xT3qvKrf4kDlrBNRtSGSeR900bqqP5yTFkYcdfXTdWYiH\/qd75H1Zk0YIUWled9M4kWFm9erVYuHChWLo0KGiU6dOjapyUNODNHevxbzu6WeYmg\/dlked6tOb6EgpouJQFl90UHApk2Ned\/u8saru9i621wCgUdTjrfE2VVcOZqun2tC6Y6LqthxPv4y3wu4O2c4NmfbK9tFxZ9rMmTObzLoxY8bIa5B1hV\/eQM\/5Hem0vRYvcMgG5qRWuq03qk+OM5eTbS6x7v6p99NjBLqdjw0SncJf0WTKlCmNHDJJz+yvsNcoKmCmzCiuLJPNJ30c8GqGF5eTbXUMhLHPnO01svDete+i65Uu0fktqrrbUm3P8zCkKKCXOXmzTLYi6NOle3JZUE646X+dY93zfOsi2mbhfRF0ZOmzdKCDSFLbJ0+e3MjWJtU7SX3PMkhqUxTQy4ziosl23iNvy9tEXY5v+gQ6H6tqp7u8J2SwhEx7JUAH8HT2tsm2zgNu3rYooPuiz6Ufmmy0feWSnMEFgHi3i0TWZVRFW9BB2VXxf1PAS8hgCZn2yoDuMql912lOQCeJngQq4p8L0NUQWBNQeV\/8YkN6l23hCRksIdMegd7uk7hu3wtLEf1lmWxJQDd50wFiRMKhcPNA7Us9M25beLLQXwQfs\/QZMu2VAV0NX50wYYLkPd8uy\/IxktqEJNFNanSWyZbkQ9ClbYI01yVqBG9V9dxF3effJAv9vudB1v5Cpr0SoPNTZ\/C88yAYOOpQTHvpWT8S2oUEdJMUdp1srgA0qew6p5tOPdeli0o6eeZKv+k7u44rzzwxtc1LexE0ufZZOtB1mWU40HmEm277zXVgunohAd0khW2TzeVMuQoW265Bknru8j6fEt3F35BnjiS1tfG+qPf66LcSoF955ZXi8ssvF81lH93Hh3DtQzfZOHB5umQugZPsbNd3J9WzLRbUNi9YXN\/jY0xqH3lpL4Im1z5LB7otMaTuwIrrYGz16izRXVVS3WRT7WkTH+geM9VznvV2URu\/dc9DBkvItJcOdHx82kPHv2fPni1wJZPutywTKalNnYHuqpLqAmYwZpMkp6AW1KGtL93CYNsW8\/UtQgZLyLRXAnSaNDzHG35TI+V8TS7qp85Ad1VJ6VAIznUTeGkrjINa5R2lh+LXEvPINtQvI6tLyGAJmfZKge4byLb+6gx0F\/uXA5rndUPbJJBT3+oVxK6Li42vaZ6HDJaQaY9Ar1HAjMlO52q97gIEvkVGoOMLAZfcPBAGdV1i5dMA2VY3ZLCETHulQFdVd5okRQXN1F2im\/akueSlvOjqoRad3U0SXN0C42Asyzand4YMlpBprwzofD8dOeIGDRokk0EiYGbvvfdukkFGlRRqzjg8Dz2vu2lPmkv6uSMPkIdPlr2xPTw1jYNNl2e9LNvcJ9Bddyhs2kXa5xHodo41SSXFz50jIeTKlStlNJzLeXQCNPK+U\/KK5pJKSo1U4wdHKNGDmqJJPVxCajuFq5J6XoVNrk4NH2Bx3aGwT8t0NXzQnu6N\/mpXJtH5fvohhxwixo0bJ6699lqxfPnyhksZdAE1IHjVqlWSA5THnbLA1j2vu04Sqb+57olzB5u6taZK+bLV86Tp6QMsVS1YPmj3B910PVUGdJCpSnU62HLvvfc65XRXgc61gpEjR0oNgV8CUaWNboopN+VZS\/MZ6Sw4qfKkjquAcFF5Xeqkoa0IiZ7n\/XnaRqDbuectCyy9iqd3TpvXvcwLHOjiA\/J8I8c5Cl2SoF50oKZJtrN2e950XEOMd+Df1Ld66YL6Ll3fLnVcaDLVAVhCvQQhZNprdYFDmgmkA7qr6o73qBc4jH7gTZmaCdtXt53aLQ0piXX5xQdUkV+AoF50QIEwpk5x+cFbmz55SgsH6KYLEtQLE+jSBRojWpsuYXCpk4c5VV1AkYdmahsy7ZVc4ACVHar1ihUr5PnzqVOnirFjx8r\/u6aTUoGexhmnu8BBvaDAx8RAHyRV8W9sjalXD9mATeq4qhHwyDj0S4WPg++Z8zq6sdkuZPDFj5AvQQiZ9ocffljANHYxi\/N86ybpng8\/\/HC5hYbttHnz5jXEvKsANr1UVy9PXvc8Dp4kuzavs00dv3pJIZ7z0NWs4yjLk53Gzi3aX5B2QqehPW3fRdcv3Rmnbp9BEuOCxYkTJ8oc7y7ba1mZUpQzLgkkphtOso4B7Qb0bCeQ3KGdx6i+rAtE2nGkAUtZi4\/rGNLQ7tpnWfUi0B3B4iK1VemK\/+vCU9N+XL4vjrbkcPMJ9LQ0Za2fBixlLT6uY0lDu2ufZdWLQHcEelbpkhbo6mEVmgh06owuKYRnvbkDvSwQuL4nAt3OqQYbnTviTM3qGOueJF2SAmFcTpWpfOBgVw+dhDzZMM6Q6Q+Z9tIlun1NKK6GDxtdB2pdCmQKU807GjWqLeTJFoGedzZkbx+B7qi6E4t1KjyX9rpLDUyfx3SUFPV11xKbgFI373TSdAx5oQqZ9gj0lEBXQU2XHNDkTqOqc2mtmgZp0j1n9R9klw\/ZW4YMlpBpj0BPCXTypOtytGWZ\/mrWF+rD5BPQTba6eaejRM8yE4ptE4FuAbqqFqf1ouPzqdtj6iGUNLnaQpYq0UYvFsxJvUegM6BzUINpqtSG9M0qydGW+sTfPDw1L9CjjV4OgEJeZCPQGdBNVwJnmUaqFEcfJjU9Tf+6yRZt9DQczF43At3OO+\/HVO2vbFojaXuNVHKdJzztu9TUyjwjTBrprXtvtNHTfg1\/9SPQ7bysNdBVu5vb0PahNa2h7n37dJaFPNmijZ5lNvlp02JV90nLdhLrP9zOxDRbYklsN+Vq423y2tMR6H4mfpZeQuZ9iwX6BYt2EC5nwW0TQlX1bR71vPZ0yJMtSnTbbCrueYsFeh6JrjradJcX0ifjVyDB045CQTZZ7PUI9OLAYOs5ZN43K6C75nXf7YybxOa2u8k8a1lUd1OWVV3UnHqNUhZw8wkY8mSLEt22lBT3vNkAPU0qqfePu0Zs3Xn3zFwloJuAjaQQuquJ8cII9I\/EunXrRPfu3YM7ZhvyIttsgM5TSyG318UXXywuvfRSeSUzFRqsK9BN3nddnnTdCbY8KrppFQp5skWJnlm25G7YrIDuktf95OkviC2793VmnC58VSeZfW6hJREXge786bxXDJn3LQ7oJ925OrXaDrBv2bL9PnIUypvufSY5dLhmzRpxzz33yMSauKMutBIy\/SHTHmxed3WCu6jupgMprTe9Lbsju53+r\/52x+AdK8cVLj\/A1VXDhw8X\/fr1q5yetASETH\/ItM+fP1\/gT2npntNODNf6Ls44sqMBZBXUHRdcIjYeOU5s3bmLfCX+j6L7zZWmWC9yoE4cOPTQQ8V1110nevXqVRhZpYTA2vK6m5JGzDimVWEDjx1HDtSFAwB4kSDHOEsBel0YGumIHGipHIhAb6lfPo67RXEgAr1Ffe442JbKgQj0lvrJL1XGAAAFWUlEQVTl47hbFAdqA3RbPHxdvwouvlCj\/SgIAjRPnjxZ7q3XqaiXdXAa6047+IgIy\/Hjx8tLQFH41lQI9NNcwEWmKJdcsn0nqUjaawF0ly24OgGFaAHd55xzjvzv7NmzZVgvBz5+v\/rqq8X1118vOnfuXJshYIIhqAcLEI0BV2Tvt99+DYtWXWkHXRAKFG0JcGA8s2bNkvylRbfO9HNQjxkzRgK96HlTC6C7BNXUBiX\/Twg+zO233y4GDx4srrjiChksA6DziYdbaCF5zjzzTDFw4MC6DUHSQ9IRNKIQaEKgnQBz3333iSlTpogVK1YEQT\/dTIz58sEHH0igFz1vagN0Wzx8LVEihJSI48aNawR0mnigGUCnO+frOAaTNlV32rn6Tqo7wBIC77GYDho0SKxataqRZlIk7RHoOdEXMtBV\/0IoQOGfjPwNZOcWCZacU6XBDn\/qqaekFFdNkCJprw3Qn376aal+mY6y+mByEX3ogB6C+sslOfkPilYfi+A\/SXZoTXvttVftVXfMjZkzZzZiBex0SPgi500tgB6qMw5fSwV60U4VH2ABzXPmzBETJ04UsMWphEA7OePwNzkTyXTq0qVLMM441alYNO9rAXQa9IQJE0SPHj0aPNg+JnXRfahA5x5V\/LvoU0lpx6duTVF7bueeddZZ8ue60U60um6v1ZV+GgdX3YueN7UBetoJG+tHDkQOuHMgAt2dV7Fm5ECwHIhAD\/bTRcIjB9w5EIHuzqtYM3IgWA5EoNfg01EY6tq1axuoOfjgg2VYZ97QWV0svm3IPEQWdVX6TPH7PMqurpGAtrE31+cR6DX4sjrPvXrgIQuZFEyyfv16550M2uclMAO8kyZNEiNGjJAhvjpaVW94nUN+s\/CxObSJQK\/BVzRt0VGkFPa6+ek+dQuSPyOAkhRGDv25c+c2yaWvDpukMfpGoUMvpnq6sF5VovPTWOiHtruwmOyyyy7iySeflPHpdLCjBp+i2ZIQgV6DT2uS6AQ2HrEGVZ5Lex62qosq9KG6cxbxkFNVPedAR\/AKPwPADy7dfPPN8ogpTvyh3siRI2VIaFT3i5uMEejF8da5Z52NPmTIEBkSzCPXqEMeaKHa0+pLfQM9yaRIstH5ggSgo1B8uho44sy4WNGZAxHozqwqrqJJdVdjnynRAiiBunvhhRdaT8fpgL5x40YxY8YMGULau3fvJgMzLR74HQ5D0wKkAl2N66bFC0DnpkEEenFzi3qOQC+ex9Y36IDOAbp8+XJBh37IXqdjvVkk+rJly8TChQvFvvvuK4YNG+YEdJMkB7hxJdKuu+4qPv7444bz9+iUFiqYG0kS3Yfj0crkFl4hAr0GE8Am0Z944okGoANYsGlhz9JRR37yT7V3VYm+bds2aRvjtNRjjz0mzj77bNG+fftGXFAXDxMQAd6lS5fKBQOe\/dNPP10elDElscBLoA1AoqMtZYWJNnrxkzACvXgeW9+gs9G5Z53neMP++qhRo+Q1PqRCcxVZ3eNWgY7ri6ZNmyaOPfZYsWTJEnH88ceLAQMGGIGu5pejinjPpk2bZIadbt26ic2bN4utW7c2SHTQSXndMBby\/iOtFjLzbNiwQSxatEiaAnXMq2f9aIFViEAP7IPlJRd2Pra2cD8cQAawnXvuuaJNmzapu4ZmAKDDcw5TApLdpR+buZGakNjAyoEIdCuLmk8F2NBQl6Fad+rUSd5EC6fcSSedlOkG2FdeeUXu73fo0EH06dNH9tOqlf0arQj08udUBHr5PI9vjBwonQMR6KWzPL4wcqB8DkSgl8\/z+MbIgdI5EIFeOsvjCyMHyufA\/wEzHWCU6EsGDQAAAABJRU5ErkJggg==","height":140,"width":232}}
%---
%[output:9ee9f073]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5a5b97b3]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"fit_test","rows":2,"type":"double","value":[["0.0295"],["1.0639"]]}}
%---
%[output:3bd33ba9]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BaB0_psap_clap\nTaille de la serie    : 203\nStatistique T_max     : 1.5349\np-valeur (bootstrap)  : 0.9970\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:0a0bb36c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:413735f7]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BaG0_psap_clap\nTaille de la serie    : 203\nStatistique T_max     : 1.4949\np-valeur (bootstrap)  : 0.9960\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : BaR0_psap_clap\nTaille de la serie    : 203\nStatistique T_max     : 1.8244\np-valeur (bootstrap)  : 0.9890\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Bac3_A81ae31\nTaille de la serie    : 164\nStatistique T_max     : 11.1344\np-valeur (bootstrap)  : 0.0570\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:42f88105]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAM8AAAB9CAYAAAAFgzchAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl8Tdf2\/iIJSYiaMiBIDEEpWvpM9Yiihla1NfZp8fAvRVGqRc3VUVR1oEVRqkrbR6ujuS3VPlVUzZ7QkAQ1RSYS+ffbsa99T86959zcm+TeuPtXv96cs88e1tnfWWuvvQaf7OzsbGhKytVM\/HvpH4goUwKxPepY7o5ZfRDxFzPwfr96KFncT\/uY928vBW4pCvjogScpOQNPrTyIub3rICy4hIUgtq7fUhQrhMleWD0dqbu\/Q8WJXyF173qcie2OYkG3oeLEb1EiumkhjMjbJSmgCx7eGPnxfny08wyuZWWjQ91y+HfLCPR9\/3fcX788Fve\/w0u9AqJA1vlTSIzthfAxH4seT8\/oiNChC8XvM\/MGodKkb+BbrnIBjcbbjUoBXfBIse1S6jX8ciLZUr9LvQq4ev26V2wrwDVE8EjAXDsfj0tfvC440LX4A17wFOB70OvKUGw7k3wNczecwJzetXElI0tXnCvkORT57q\/s+NRKVCOIKLqFjvkEpZo9UuTn764TtCm2SeXA6HbVsOiHU5jYuTp6LdyLe2rcZqVEcNeJecflpUB+U8AKPFIhMLlLdfz7g\/04ciY1V\/+1QoPw+bBGVoqE\/B7krdy+uudR9za2rt\/KtCroudvkPAU9EG9\/1hQgOE5NbIHMcydtkqZEdEux\/ykWVNpLvkKggC549iUko9Mbu5B8NQvF4GMZ1nVkoxiA0kF++OLJu3BHpVKFMGTbXX788cd48803sXjxYtSqVcuq4pEjRzBgwACMGDECvXr1srp3\/vx5DBw4UFxbtGgRypUrB3mtWbNmePbZZwttnl4OU2ikN+w4l9jW9e3duuKabCmiTHG8+FAtLNuR4HZaN2fBs2fPHrz00ksCXO4CHsM36K1QaBSwyXkenrcbnw1thPoVgy2DU\/dE07\/8X65D1EKbxY2OXQGeSpUqCc5Vvnx5wY0Km\/NwatdTLyNhZmdkHN5mIbFXZCvs1WbjkFSe8+w8mYw1QxsJ8ez301fQbd5uNKkajEebVsSKn4se5+HrOHv2LB544AEMHjzYLcAjgeNbPgLhT6+0rJjE2b2R9Ve8d89TiBiyqTCYv\/VPTFh7NNfQoioE4rZAX7zRq26R2\/OQy7B88cUXmDVrFl555ZVC5zxebVshosOga5uHpP0W78Nr3WtjybZ4dKofgnZ1y4NnP7VCS2JI6ypuOSNnxTaCR3IclQsVpsKA4yCXSd\/zjcWWLePwz0iYeR8CGna04kZu+VKK8KAMLQz+s+sMjpxJEQej7m4Y6grwEChsZ\/z48eK1P\/HEE4WqbZNrj8ahF1ZNsSzFsj2noWyPyUV4abr\/1GzatpHzxF9Ix6EzafC5no3sYj6oHRqIrGxg3fA73fKQ1FXgkZo2at\/cBTzuv5RuvRHaBM+ji\/bi+Nk0rB\/dBOQ+cv\/TKKIUvhh+p1v686gcQ32VK1asENoznvOcPn3a6i0THHrKgR07duDRRx91C\/B4tW3uCUxDsc3rz1O4L86rbStc+tvr3a5h6Jq953Kpqrs1qOA1DC3A9+nVthUgsR3syq6qeumO0ziYlAI6ahfz8UFgcV9885T7meU4OGf96nFxwJYtwIkTQLVqQJs2QGSkS5p2thGvts1ZCubP84ZiG7uVamv+fuaTQ1g6oL6VwmDDgb\/Qc8FeywhffLCmUGer1+U1u9OIiQG4iLloN2\/OnxlrW2V\/AwbkAEdbOI7Fi3OAVMjFGW1bfHw8+M\/TSkREBPjPXYtN8HR7ZzcOJqWC8UGK+QA5UUJ8EOjvg+9GNbYy2+GBqlRny4lSrW0PdCfPp+Oj\/yagz90VUbVcADIzM4GoKPjxJXPRHj9uoRnvJScnIzg4GH5+uQOP2Ltv99klS3KAw8I+VW4zderNd0YgKwDKa3+ZZ+OQvHkpgmP6wS+kYLgaQfPMM8\/g559\/dtc1aHNcTZs2xWuvvea2ALIJnraxO5Fw+SrKlvQX+56kSxmCu6wa3EAcmKpF7\/CUXGfi2qPC96dUCV8RjWfwPRGWZ388ehFd3\/lNtNeyZhlkpGfA7962CEhMFE1fO3bM0gXvJSQmIjwsDAGBAbmIbe++rXv+J0\/Cr2Z10VbahEnA1JtnKLzG59LHj0f423NFncyNm3Ct1T\/F77z0x+eufDodyZ9OR8WpmxFY7yY3k56irOPq8xupNeQirFzZc2IdEOxvvPEGqCmVlh\/uhn5DsU3VtukNnhxGtcSWznK\/x1\/Bgh\/jheU1C8ETE13OYp0gwTP\/X3XRuGppZGRkCLuykJAQlChxM2KPWKx27hndt\/VspW4dEbj9B8Q9MRypz0+32Wf0kvdQcW4s0lq0wuk13+SAJ69j3bEcPh88YQUeNUYB29YL6uGMqlqCR3cRyn0eO+Zeb4r1B6QwF6vdcRfmwJS+bXqSPjx\/D85euZZrmEaepBThNh8+j8ebVcYHO04ZguezIQ3RNOo2pKelw7\/dvTmcJ7Iarm3cZOmb9xKTklAxPBwlAnKA5X9vWyDuhKgrOMOxY\/CrUcPqOV7Xe5bPBUbXQGZEBI5s2mrVruxUfa5Mg3piL5ZZpZroL3ndV7nGo\/ecHKvl3h9bcWHmvVbgIdc5v\/RpVJ65HT4BwcJ6+rYHRltiEzirqra5CKdNA1TRVA6S19wARB4HHhUptgxDjcAjxbVxHaIwZd1RxF9MF+50dcJKYkbXmlZiG2La4B6\/S2K\/kf711\/CrVUt3z5Oeno6EhARUrFgRAQEBABULcoMvNWJS0cC9kqJ4YLtWz3KSN\/Y66U2bIuGjj262qxIgJgaZR4\/Cr2bNnP0Qn2HhWA8cyN3mjWdzjVVpM+2PLUiYGpMLPDIiDqsSPEGNOlhMb5xVVesuQnWv179\/zvxIT0lTKkl4\/UaRjoTqAbM9ywtaaOzevRtt27YFrTSCgoLw119\/4cSJE8JXatOmTWjUqJFoXdbTfqU9FjzqZt\/IW1Rr78b9D8uAeyqj89xdmNO7rvh71MoD+OqpuyyKBoptVRvXRdXL5DQ5C5Kb8oCkpBw6KgqDXAsyKipHK3djMYv\/S\/DIhSDvRUbmgICKBgk0BXjc\/AuAUCmgavvYnmyTC0l+pU2AxzIPjdbwwqppuLB6qkPg4TSoaUvetEhwJ8YxkC7awW0HGtq35VqEnBPpx6LlMqQLacCiBJIleL799lsMHz5c3EpLS8PMmTPRr1+\/XB67vM8+JVDeeust3HfffZZ6fHbu3LnCquPw4cOWekUGPPTnGf+fwzh5PgPfH7mI9nXL4ZWHojFkxX5haa0FlKqSVvc8Iz8+KJQOLDTrmdCpuhXnsYDnxkIXC\/nGIo\/ZDMThBkAQifWtMlFTqlslcLQUN\/pb5VLautSmScDYa0cL9MhIyLFGIhKbufZUrsg+b4w3LSUOCa1hSmwLjG5mGMPAr0JVC6BsDTkXeCRAOF+94wDJ1RXuYw885Cg0Y2JRTZ3IcSZNmoQZM2agYcOGGDp0KC5cuCDq0ei2Q4cOgnv\/\/vvvukoBj+U8BE\/Xt39DakYW3nu8viVu27KfEsR+xkysagLKnsIgBjFYHLUFkRIfmrcfxzUXCbEwWY5HwWZdI8y4\/L4CCMH0lLHaGyfrtfoFqBnSBpuRMzEzCgNnxp9rEcq9jq29jc59PbGN7urt27fHtGnTMGXKFBH3gbaFLNWqVcvFeaTY1rVr16LNedQAID7wuWFVnXPWE1zcF1+PvCl+2XqxRuCJQhQWx8TlAoQKJi62qBvHPfyit9E5x3RmYbnyWTlWW+DhfRbOhxzqOG6eY6mqalcHMswFHrnfoShK7qItkvMoZ1taziMfsbUXat269a0LHhLH2aCHRuc85Dw3xbKbb5AgUQEkwcMa2ntOL35V\/NOa4ti7Z6NjjtXeGNW5qOBxeh52GtAVf3xuRETSHP5KRYr2kNoWeKgYUDmPHIbenueW4TwkAkW3B976Dbvjr1heTXRoIDY+3cSUO4KRhUF+LhjTbasbZNXCgMBRtGtC8+QG6lvT81Iq6oJHVVNz3pyfOmeNSGcLPOxGti+7pDhHsY37IP6Oi4sTbu2jR48WViQU25577jmcPHlS7IOGDBliiVikzs+j9zyuyM+jKhL0LBPyshhc\/gwXjdSyaRt3E9s29ZCUUXNChy5AYmxPkS3BKMWIzUWoqqvlvCWQ3OBD4bHgIdfos2CvFdeR9A2iZfXIO61s21y+oAujQYJo6VLRc1ZYaSTuex\/h0ws\/fYd6SFr2\/tE4v3oawkavxKUv51hy9tiLGGrKwsANLck9GjzNXv4FmVnZ4mxGDT2VF38e1cOTaksZlVNl+TLYoMQNzwPI3vv06WOxbbJXX8Wbs\/3dHV1F5MT5LaITrq99ATPjQvHwgCcLJZaBekiaeS7eAp7s9GRL3h57+Xk8YRHqfSs9Ydx2bdsaVyuNl7+5qUs25VagoQTl5RdffBGxsbFCnclwTjypHjduHJ5\/\/nlMmDBBPKHWUWMISJssXhszZoxufbVLV\/V3ZFoX+O77yiYPNHPG4ioGKmO0lesxBRe\/fB0V+sUiYWZHFK\/V3DB6jicswiIFHk5Gz83AFYuBL\/Ojjz7Cww8\/jDlz5gguFBgYaOEy5Ew8ge7evTvGjh0rvva0quVzBJ62vpHFbV76I+eJf+URvJcYigHlT6Hs8A8wdsYsAVxtDGxX0MRMG6o6m\/XNWl\/rgkc1CNXrXMcRkLR\/9913RW0ZVZW04HUWGZ5LauD4cVy+fLmwJOBHk4UftqVLl6J27dribEgtapu87gmgt8l5bMWsNrJtM1oIJHZkZKTQyBBEL7\/8sniEIlqLFi0sQdgl91HBY6++rX7z0p8U24o\/MhV7Z\/bG8uy7EbvwQ8siMJqjO903tG3TG6zGto1tbN261QIQCYKJEyeKj9nGjRstHzkz4OFz\/GCy7oIFC\/DUU0+Jv9XiseDhJLTeoXJizoCHexGqLgkIyRHyEzzO9PfHxPYIPLzBLcQ2Z8BoFzw851HPt6TmUccwlBxDLnp1PLRda9KkifgQSm7Csx97nKdIg0f66Mx88KYVtDMvkM9KDiDTexiJYXqcxxGxzdn+CLwzR\/bi4bStIpnurAXLrcQTZ+lh9nl7eXrM7LvsgofGt1rw0GhUAx4pcsnQXaqIJQ0\/eQhK7kQxTYLn1VdfFWc8amEccH4wiyTnoZZr1JhxWP\/Nl5Y5l617Dz58dw5GfXZc1zDUaCFwIdNkQ92fGCkAtOAxqq+Owdn+KLYdmdoZ8+LLYnTIMVw\/n9v\/38zCNaKLM\/epRCjVoodhTlJXgUcdq2pZQIWOtJom3e+44w5hgX3LcR6pHg4JC8fmwHbCMPTNnjUwfNQzoOatzgMj8OOJVFOGoZLY2hNoXpdfH1reSotcraejFjzqJpK\/bbnnurI\/dwy7K+lqNumVK8AjDT6l1KC6FajgUZOEkSPZUhgUSbHt0MlEPDhwLJKie8CnRDCui5AfPiiWdhkB+z5CaoM+KBVcxpRhqDNfVe+zxhRQPVAdPueR1gXaPQ+7tSG2qdo29eOl9deR4vgtBx75defmz+cf\/0ZCKsBs2OMnz0T3Ng3wwdnoWyobtsxGcD31kmU1F7S4Zm\/PY8YC2xXaNmMou76Gx2nb1MNJW+RQLQRcTzL3aVGeq2jPU2T8NDML1x1m46pznoKei8eBpyAI5AlEkfZkaiAOlTYElow74EwmarnHZNtSA2XrHaj7DHnoaOZ9OUJvR8Zjpm9n6jgybmf6ceZZhw9J2ZmM5WYU30BvYJ5AFKPNuNF9sy\/E3mLVimsL4stiT3IAZtZMwm3+3I0CZkRISe+RI0eCQQTtlatXr1qsCOhSXbx4cbNTcXm9U6dOiWCNbh23LSoqKpuEoukF\/8\/Tf5nYSaVIdN16OJoSiPLBgahbNQTff\/0ZKlashIkTJwi7NNqr2YqootVasS\/VZo2Jc6l5Uzei\/K1uUrXnA2aeUUVMrUiqGqJq7818fhz+eWw5kh96GYNGPYfevXuLU3SZr2fs4L66RpkSDPRVkcav2lWl9kX1PUvp0qUtnEellZzz559\/bnknPGOZN28eFi5ciMuXL4vneeLPJMS7du2y1JPPchFy\/DyHYfHx8RFRYO39lvdcjggHG3T7iKESPDz154vt1q0b\/EoE4IRfJBo2\/gdSfvscSZevYsTAPlj51fc4un8vZr\/9Lt7fdhrn1r2E8LBQsVDWr18vXpz2SyFddfky2YcEBOtxcdMsh4ViCxfJm2++KRaCDCyh5tYZMWKExZnK1jPShIdhZnmop\/bLr7Acq+xHm\/Vafqnnti2HmkPmCPA0btzYanwL57yM4P88Jw5PVU2XGfBw\/jw45BxZOEbZvlTdc87R0dGWhML\/Vz0Dr86eY+E8JXyzERtXHn+klcKSD5ahTpN7rGzB1Gdp+sJ\/DPEk7c84BtoWxsTECMuA7du3W93j84XNeUgbt49VTfDIBX\/s2DFhkNm2YxesiSuB65fPoGdMA6xd9h7q1I7G7oQMFPPzx4RpL+Cr304ja3MsWjRvbjG34ZmNFjzabG2q2CYTThEUPEOQQGMyXRa9MyBZR+8ZeQircizJDeU1rcJDK0YSAM8+\/RRS96xH99tOYer\/whFTLgWDIy5g9+USGHukEmbVOo0m1cMskWtsKVq0xo6ynkxPrxXbaBArAS4t0E\/FHcUI\/5+wOvBebN+6CfOmP4MKrXpgVLdWKBZ0G2avWCdO6zk\/7bOUBiZPnozp06cLz01+OFjItdX8q7bGo7U3c5BxFPnqPip4GAaIi7hTp064HNUBP34Yi\/NlbkeJY+uBgLLwL1sRYWVKotQ\/\/w9TOlXBrIlPWbJF29rL2AOPChCV0hSppLuuat5BYNp7hvfI\/fSyvWmDVUgQkWNKLkQrYXVBDxo0SLgKS6Aa7deMOI898MhFrjVnaVCvLmZGn8cXpe\/Ftk3f4bUOVVDt6eUC4GkHtuH1FV\/At2wlwcG1z3KOPG+hmYzk1ByjFzyuwbUVePhyGciOQQa56N\/44D9YNnem6OlKg75oEpgAGoZSxFJfgjT0tMV5VHHOHuexNSW5D+BioA8QN79yQctntF9x7fi09bjQCDLuO9RxFwR4qlatamVNLhe2lvPwutT8LTxdDrsS0zA9\/BCqjV6O8VOm4VrSMcz5ZCNKVaqRi\/OYoYtMJ2lrPF7OYx9kVuBhVZX9qyyecrO6PzELnrzuebj5lRyBsq\/sW36h9fY8n332mUV0kXswAkQ7drkfIgAZe0x+ibUfAa1YacR5zHzP8rLnoYJiymOdsK9YFcx+pCHSv3hZ7HkC6rTCrIUfCrFNHZt2z2Pvvdkbjxc8JsAzf\/58ocGRGi+9RzLLRKJ8SBhaRIdg1quvmOY8bEu1N5OaPUe1ber+QbvHkG2pohnrM+OC\/KpKsMs5mtHeafdXrgCP5GzkfGa1bSo4SE8qG\/ihkB8Qucj1NHWsYw88RuMx80G4VevYTSV\/OCkF8RevWsLtdpy7E9VDSmL1Ew1MhZ+6VYnqnfetQQFd8Bw\/n4L75uzCvXXK4ck21TDn22OY3as25n53BB\/s\/AvfPX03osqVdCsKqW7K1EJVnPitCMukHjaq123V105K6\/4s75s5oHQFgeT4fctVQcWJX8EZiwZXjMfbxk0K2LQwaP\/6TiRfvY4p7Spg6aL3MXHsMAxasAMlTv0X3703CdVCyuSi4++nr2Dgkn1Y1L++iLjjcE5SGf5JDbhOf3qDOGJqvGcCRs0qcHbJGDHO8KdXiuupu79DyNBFSIrtbol7ps1CICcmF265frMN\/Wbyc1HpJbcKat7LMPhHfo7J2zZgN5X8tmMXceToUfhkpuN6uRqoUNIfd\/kcROwTnXKl6GOEUQZK3HkyWaRhDA32N0wEDGYUYHAexnFuMy0nuwCBoibQJZAY4ZKn8UrOGHsvj9bQMsPa5Y2LkPHnfgt45G\/1ebW+eujpCjOc\/Eqmm7b\/eyR\/+w4qDHoLxYIrFMm17PaHpNnSVkOH\/IJzvLUdgYfW4lqZ6ujVOAyRxS+KPC3a1IeMtpN4+Sp2HL8ovE2Zw9ReTlIm9EUUUPUyc4wuwbmRv8NvxgykpKSIgPKlSpYUvzOzsuDr64vyT48CWkbDv6dxMlxyEhUk0hLaVsQZbX2VFPbuGa1YT06mazS3grjv9uY5tsBjFQDkehZwLU3Qq1ZEBXw+\/C6rVPIU1+ZuOIGJnatbcvgQPEY5Scl57om\/iJSwPti+dgEqVKiAc+fOiX7kbxor0kCxcsr\/ENX7Qfh\/uRl+d9tO7a5aPLMdNdOaFNvUvYM9C2lXxQ\/wtGS6BQEMoz48MqEvVZdDRo7BDxtykteyqGpd7aQpro1aeQhPtatmJaqZAY\/Ihv3nx4js1Awp\/2iK0LAwJCUliTSM8vfVjAzB5UKvHIffiC4IGZaT8o8+\/Kk\/5eSDkb41Wi6h3QsJ57bYHqg4ZrVQJjjDVYxePu+7QrVtpp+iWMcTaGe155E6\/9sqhGFNdmv0bVYJI1qFW50TaA\/OyHW6zduNCyk3k\/\/SZWFipyjM3nDCICcpOc9oZEYsEO8\/s2o2Dr17VIht9YfWEvuhlNB0HF\/6JyqnxiH5tc6oXko\/r0zS672Rsu1j+o0Lmy9q2\/wj6uL0tHbITDwM4Q1aPAi+QaUR8cpOcC\/EUvregSIDmy2lgBT5tAvUjLbNExaAuwLPE2hnBR4ZFeXJMePxwIKjOHclA\/ApBp+MZATsXYH0Bo+iZpUwfD6skZXYJl+AmlbExyfbbk5SoSwQyaoG8NjP8g7\/apSMXfNPoNXDtRFw2h8pYenYvvaoENsC5j2Iqvsj4ffHzcRQfDBt5zokvNYNoHgpSzE\/RLy0A+dXTUfavg3IzkiFj38A\/MOjUb7\/60iKfTgHUDeKb3AFRMzabWUlrXKuC+tet0SrIdcrUeV2x\/OB5nGl6iWRsufnoh58yi5VdxFpJKs1XM3j8AwfIxBkjlLDyjcqeBx4pKhBM\/WaXYZjyjc5IZcC\/vgEmSVDkRn5T2HbZgY8FNvs5SSlsiAn5ag1eHglvdI1ARxZ+Lefrx\/w11b4hQ6wSvar9zLMatu430neuEjYj+l5jaraNlVrZ1YL56oFoM2Po0bs1DOh0Xqdqgl4ZXw1miIZtWN2oRvVuyXAYyaGgV9IDXyxehlqVw23SzOjtIo3wTPt75xvdAqzrQS42dEWINIYPGa0bRIAzHVzZt5gXfDI8xWmdg9qeJ9lv5S651ur7NS2CJGf4OHC79y5s8WNQ6Z5l35SMuM0XRsInpUrVwqnOPpMMdQx3RDsuXbruXXoZZ9gRut169Zhy5YtwiFSm+1CfpAl59FzcOSYGElWxsJWrefd2pNUT9tGC4O2s37Fq92j0eMu+yCxtXDMg4ctUIa7kbnXDiQzg2LgFxpnl\/OY1badmf9\/QgwLatBeaORsxSvQcp8Lq3IClJsJAOJK8MhonVoxjLZ6NASV9m+M2ikNYVUXBdUFQwWPXkp4vRjSBBqBwnZYZGBJ\/pYOiKxjL80ijXCZEbtGjRpW6ehpBCzDMEtuyGzZ\/fv3d283bD3wyORWdcKD8M6j9Yy4su59o5ykxEv8j\/GIyIy4sflZini\/GVZt5dwD4v1YbxEQWQ2Zv7SBX0ikU9q2kP97B38tHIbMcyet+jMDCEeIYQQe9YzY3qfDXip3aXEugSKtyFXOwzHLwIX8reU8zIH0xx9\/YP\/+\/WjevLn4t2rVKotLt9wvafdekkOo+xkZx+3PP\/8U2sYuXboIoMg6qpGw3HMRPNox3XXXXcJlxKM4j8qa1YWSWSYKKc1H2lUYqPXN5CTlIWLEYxE39j5LEO83E\/EvvMC3K5pqNr4ZsOXv\/KCR04D+kYgfOFC45uoVimosZXtMttxWxS5e1wsUaBQpxxGwaOsagUdu+2hgYa0CsW5JDzwSHGrETvYnOY8eeLhAWWTGA\/mVp8Pfjz\/+iB49emD16tW49957xTkbi9wvPfTQQ9iwYYMlo4E9zkMXD4KHYXjZB7kN\/yYgGIqXh+zqPqxIcB5GDO3adyjORT+CEn\/+hGuh9ZEVUse0wkC7eBzOSWrLtq1fP+uA5JqO7AUnZFWqosllVMNQ2YQeeKSoVmHAbCTF9sjFofisK1TVkvOwPSPwaMU2+UXWuiHwa04PWLpeq2KbPW0b\/Za4gNu1a4cvv\/wSVapUERkPpPuG6hOlcjiZKkbd88iYE\/wwMkYCg5SQo5HzSN8ptstxtmnTBvfff7+4J\/c8khupMSyMcjA584Fz5lldVTUJt+7rbzH\/y504GH6\/RVXdtPswLBva3OuOYJLiRpzHZDP5Xo17nDfeeAP169dHVlaWiGPh5+dnql89TRq5yrJly4Q\/VXJysvDUtdcePwBSbJOdegLtdA9JG\/+jGVafCMbpb9\/GyoVzxVeJ2prKDz6Hs9cCHAr0buoNuFkle2Y5cqiu4DzuMm2Ch9F12rZt6\/CQ8qKG1nZSJMDDSZGQo597HtvL3I+hNc\/irZepSs7JSlCqWj08PG83PhvayG42bGcT6ko2rW4u9VSgDr\/pPDygdyDKvdOV7asNXQI84etJkmRkZODSpUsIDQ3NA4Xy5xFPoJ0V55FJrY6cSQWyr+M6rQtEpgQIkxn+q1qmOH56rqlN0c1VCXUdyceTP68PwpGOWbG18dkK+pA0v+bnzu16HHgkMQmip1YeRJNqpfGSk9mwSYSCTODr6gUhDVClO4P824wzmicsAFfTy1XteQLtcjnDUeRiBEmGl5IRLJ2xgcpLQl2Z\/TovCXxd9fLUdsy6bGv79oQFkB\/0ckWbnkA7m9o2EmDMmDEifToL41HHxsY6lBHamYS6kmPZSvjriheU3224agE4ahgq56U9s1P3je5uHOoq2uXnO84FHgmY73Ydx4tvvY+0uwbAJ+UsAvetQlrjQaYPSZ1NqEviOZLANz+JlNe2XbUAHDUM5XjVDxf\/lpbWtCZg0R6U6mVHApaXAAARAElEQVS6zuu89Z5zVCvnKtq5cg7atnKJbXLQxUqWx+RZ7yD6tiyriJpmBuNsQl2Kbe6gMOBc1eAbJaJbgoakibE9LcFD7NHDVQtADzz2DENp66ZnY6ZyJHczDuUWQT1sZTTTBQsWeJZ5DglMhUG3d3YjtLQ\/VgxsgJ+OXULPBXtN5eVR1cvyZeVnAl8zYM5rHQkc3\/IRKHv\/aJxfPQ1ho1fi0pdzRBQeo1BQRuCJQQziEIdIRGKzHcNYPbFNWgzoGYYyA8LSpUshuYkqvlF0Y3E341BaFKgGpjRu5d7bo2zbJHj6LNgLGoZOeaAm+ryxGQ+EnUNQyVL45HQFLH\/ibl1nuLwuUnd9TlVJZ56Lt4AnOz1ZV4WtnYcReKIQZQHPcTsGOo4ahvKrbc+6WQse2sEVtnEowaMamNLfiHZ2HgMeq3OeGyuhWOo5+P\/5M65WaQaf65kITPwVn84eg3\/UzJurgrsCxda4qJrO+ise5XpMwcUvX0eFfrFImNkRxWs1d\/qQVHIe9u0oeOwZhnLhafc87EPVfLqbcWiR4TwktDTqLHPtDF7tWgVlImoL0e3Zumfx+P2tULFiRU\/DQZ7Hq40aaiuElaOcx+yA7Gnb9AxDKfLQv8cRbVthG4eSFuqeh9kwXnjhBc\/hPHovk2n5mCKkZs2aIh0fX+Sdd96JevXqoWPHjmbf\/y1Zz0hscyeiFLZxqFYb5wm0sxkxVEYAXX\/gvOUdt69bDnec\/QrvL8yJdsNSkDZn2qwF9haf+jK0WdMKatF6wgKQtChs41CPB4\/qX64uMOnP0eHRJ3HxzGl8\/ekyRISHWvznZd7P\/F6UZsFjtl5+j9eTwONuxqGeQDvBebRZ1S5nFsOTK\/Yj6zrw\/cEziDr6ISJLXsNvO39GRI+pWP18r1zaNnuZ2WQGMop69O+gQxV94XmdXoUqB9MSTU3LyHp0CtPmI5XBL8gFVYcr1ifwWdguc3Ly7ED+ZnAMtX2tS\/PTVc+gc0hKLozFpfphwvEInEmnyWyOxbk0KVKzzHE87733HhgxlH4yTNZFL8sZM3LczSdNmiRcASQ9bGXQzm+Qu2P7HgMeGa+NTnBcUBTZ2r3+KxIvZWDxgPpYvOEgTnwyCScvXkfIAxOwaUJMLqtqM+DhS+ICYSAJLjJ5\/kPNkcwQrfUgtAUembOUKlYuXL0MZxJkqtjGvseOHSsSRKkZ52gGJNMavjfrBWwf3xmjd\/qK5L0t7utm0azZS8Ir5yXBpAUP947afKtqtu+C4uDuCBbtmDwGPBy4FNlKthqMP8s3RTGerv9t01Y87geU2r0U18pG5fgktB6Nr8e1zuXPYwY8MuuyVqxSAWIWPEw8zKLNBqemaNcDD5\/h159hmOhTLzmZNr2iDMPV6+EH0ergQmQc3ib6S7\/ug9kJVbEpsViuMMTaF67HeTg+Jg62lb5Spnv3hAWen2P0KPBoF6IkTPnQcNTsPglRV3bhk6+3oMPQF\/BW\/6a56FbQ4JGL\/ezZs4KLqDlMteKdrTTrjRo1EpxQBZwUAeUEVd9\/ee2rpztj+JqDFhpoLShscR4167Yj4HHUMFQvYqhqGa+2pzc\/V4PCXnw4W315HHi0E1EDePiePYjgba\/Z1LsXNHhoXqLuL+ztjbTg4Yuh6KaXs1RyR0kLPZds2rjRNOdieqbgYjSRIZfRjslVnMdRw1C9xco5y8g6aqw2PTtEL3jMUcCmqpr5diasPXqzlayrCNy1GCWv\/mWJGKqX6VpuxqX8zy+bVBiYEdtUriH3Nb\/++qsVd6A4pidySS5kxHnUyKhaM32OX7vn6TDtQ0tmOPmsnIv6hVSzZ3PsjGJD+yxVYZAXsc1Rw1AqY7ShpyR4pMGlPEi1Ze2sJw5zfvJjQRrL\/Rv7Ise2pfhQwUyxXI0ERC5NQ9bZs2eDZ4oEODn5I4884jlBD9UFdXv9O5B610D4\/3chjhzYh3Ih4fho+VLM352FU1uX4ZevPrKASmX7KsFJgJMnT4rNvCPgkfsBhl6V4YkYytXWXkaGaG3VqhV++OEHwRm5eJgBmi+T46Bv\/s6dO4WIRoWI3ONJJQX7ZNGKO2a1bbbCOtWtW1fEOrMJnsWLsSg9HahaFQPLlxe00tvzOGoYai9iKOdJjaMKHmmqo35v9YJyHDt2DGXLlhU0lDHfRo8eLcJcScWNDK4o96SSrhLMjN+mjXDK9SFdYeT7p1TgERYGctG0aNECnPTkGS\/i40\/WILL78\/jvhWCEJm3HPcEJOFa9NwKCgrCkX71bwjDUHPO2XctQbo+KAhirLjLSbghhRw1DCQyV86iimTaUruQ8t99+u9g3MpkYQ0\/RmkRyGXW\/pI01PXnyZCxcuFB8ICWoGNiQ64hACgoKEh+w999\/X9Sh2ZD8sJFyUjLRAnrz5s2e4ZKgVVX\/vGcf+jzzOmp1+DeeaF8Loxb\/gn9e+QqZDf+Fbacy8f24JlbZsNUvtkpo1T2hIC0RnF30rnreEDwxMTngYTluO+yhsxFD1eicVM\/r7XkY1ZOckuXAgQPo1KmThQxS5JNBDskhjTgPH2Yoq8zMTAFIhi4jePQinHo059GCZ8fuPzB82puoGNMPYzrVxuSVv+HYdwuQ2qAP7qpZCSsGNbDiPKrXqCpbO+vG7apFXFjtGILH5MAcNQzlXkvlBuxGjWqk7jukyLlp0ybccccdYkT8zcJzKRZtFE9ea926tUjqzFC9y5cvt+x5pOZx7dq1IouDv7+\/CG1F8ZogWb9+vaVd1mXbffv2xbx58yx7Ho6J7asKIZOkKtBqQmGgBQ8J\/cTkN1C+1eMY0yka877ZD\/89K7An7EH8s34VxPaoYxmk9ll5gwvH092onX0TWvDkV2ZsZ8fJ5xmrmu+SkT2pwueiNlPIVT799FMRNjc4ONjyyKFDh4R\/Di1KWrZsifBw2y4srEOLaioJmH+WhcqDZ555BiNHjgQT++qVws6WbQGPqkXRGygDvY+aPgdP35\/zdZJFgod\/cwOufqXsRb\/hQiLRGUWfhM08ehQBAQFWManJ8q9cuYJSpUpZwrWavcbx5FddLjAz4WhV8PBFczEwUa23uIYChZ0t26aq2uz0pEgxa9YsoS2SoavoCky7rp9++kk01aBBA6GRkVoYubC4F6pTpw6CZ89G9Q8+wMEDByxdk90zOzYj9svU9WavsZH8qsvxyCwC9uikgof1KIZ4M2ObXVn267lDtmxd8Migh3N71zHUqqmBOqhmlOLaY489JlyBqR7mwqccToNIJi1ikQuL1\/hcxQkTUPbSJRxbuNAKPEGrViG1Z08r8FBm55mDCihtPQkevbraa47WDQsLsyuGyAnogced3Ypds6wLphVX7SedGa1u0EO5UVQblpnF5DmJek9VGEjOw3QVDEBBFaU98Lz66qsCPCHjxqHcpUuIW7LE0nTAjh2IHDAAh775Br6+vij16afIysxEmbVr4Vu9uqgb7+eHbdv8MaV\/JE4NHoxrLVogvVkz0YY7ch4veJxZrjefdTvwSC7yxMixGDjtPZwpVcdUfh71gFWCjBYGtOWSYhuzjfHUXRXbPu3aFV2GdcG6MeUxqP2HSAxIxGPr70HLzMoWKnUc+jH6Hr0nJ+HV1q1IDw+HX3w8Xpjhhw19K+PnxCQgMo7\/YfEA8ROZEREYsLg1KkdEoOzFS5iCKVZcas+lS0gLO5Czx0IkpmKqeK69f3t0KdlF1D3lfwqZ1zLhc9IH5DRruq3BIP8PkdmyJfxq1oTfoEGGq8DLeXISB6hnOJJoeoewhgRVKrgleGTUFaoUZZ5IWxo1I3nfnsKAk0dMDMLD00UzISkhWPP6g5h1\/04kluCijkNWVhWcLflfJAwdgjJr1uBit25IGjIEPyz7E4smTRIpGSOuRSAkNQSpqam4VPaSyJo9bMROATCC0X\/QcpSqX98KPI0aNkR4ekVxPzw9HB1XJmLKNGBLG2DAjaz2bPvdw++iREIJ1PzxR6TN\/waXRnVDwzlzckBn51zGK7YZw6DIgUe1NMgoUx2TJjyHK3cOgO\/5YwiI+96htIpGqmr55XinZ080uXIFp6Ojkd26tWWROyJ26YlnBBtL+PjxYh91rWXLG3utABy+ehWtWx\/H+fOlUXfVq6jzU4IQ+WoMGoT5z\/bCwWeboWl6GG4\/c3suhQUbcUZh4MliG98Zz2OkGRRV0VLEl4fg2izaqi2dGvtczQonY8jJw2BqfqUVgq046W7HebgwyGXGPz8Fu0K6onPZk1g2d6ZYdI6+dKOIn3Lyw4YNE7ZoGzduFFo3mTeTfVLVrL1u9pp8\/udXXhHBS8Ju+P\/w+XseewxxkZEode6cOOHP7NsXR\/\/1L5xYuhRldu9G+TffFHPW64vXedLOf0bFSGyjcYE0MDBqq02b3DW2bDF66ub9vDzPKarTlAfgtL+TFgYysCJBw8Awv\/32m8VuTq4nim08COX+lg6X0kSH2lcCUAse2rWxUMS3d45Y2IeoTmvbjEQ3TlAPfDznKaxzj0euXBFjorLhZ7Hvcazw4I7\/jIoReKb9nTds6lSjVnLuZ9MRUVN8crzATZW8PM+xTZlys3nVAls1vZI1ZDRSyY1UuzVaG1BUU6OYEjR64ElJSREWCbLocR+35DxywPTl+XrfWStrAlNvyYFK7nzibm8aZk+2jcDjiZxHRvXUch4tvbRZtFXOQ42tNFZleyzkMqoNnbxm6z24FXjIHh\/vPwD79\/0OWhOkNh6IoF8Xwe\/icWQFlsOV5qOQXboyaoUG4fNhjQzPfxzAUJGtagQeT5u4ynk4dm1Qxfnz54s9kZpFW2\/Pw2fpQkK7OOkDRX8giu9MY6\/ueVhX75jEbcCjdUlQg2nwDEZNeEV51VvMUaCogcfcrAumltuAR88wVI2ynxdVtSRhQSb3daQvreOb6tCmyvPOuFJ4wZN\/QHJb8NAl4anpb2H1uy+hWkgZi8bjyTHjMfnbJJgx2yHZVDN4KeeSPY8bNw6MRayXdU49cJUaPiPNXV760p5jyYg6DCziKlcKL3huIfAQHI+vOIGjR4+gxP82Iv2O3oBvcfhkJCNg7wqkN3gU7RpVw\/v96tnMhm2PXFxMBZXc16gvma5ejlfNIOAqVwoveG4R8Bi5JNizbTNLooJM7mu2L45d5Wxqqgve42GddE83O09ZzwueW8g8x9HF4Uj9gkzu60hfaq5OV2fh9oLH9gopcuY5nKre4ZctdaFZ8BRkcl9H+tK6U8j5e8U2\/TfrNc+xpkuubNgU32h+od0TmAWKtl5BJvd1pC+OUy\/1oBnlhFlaGHEehv64Ef7DsEkd6xw4YJ2DvDxPAyTVCMlrnmMAHntZlA3fsKaCHhfLr+S+jvall05FqqXVthy16VNJYASev61z\/naHMFd0rHPggHWOCDOuLUbPc2yKdY6QSqSFgdc8B9B1hqMrwm3NHsXmw+eFZs2RbNjmlsKtUcsIPJ7IebzmOTfXrmlP0lp162Pya29jxoYzWDqgvtc8xwT+jcBjogm3quI1zzEQ23hASE\/Sqd+n47XutZF0KQMLfowXHOh\/59LxzCeHvOAxuaSLGnhMTrtAqrmNhYGcrRqQe1dStsh+XbakP9YMbSRAxL9XDW6AdnXLFwiBPL0TL3jy7w26HXg4VVWjkn9TvzVa9oIn\/96z24FHtSvTTptuCk36TcWyoc3zZJqTf2R035a94Mm\/d+N24JFTlWnkI8qUsHKGG7P6IOIvZuTZti3\/SOmeLXvBk3\/vxW3BYyvooSPBEPOPbJ7Tshc8t5Btmz2xjQmvku4cjoFt62BI6yqes4ILcaRe8NgmflGwbft\/LfbHLScyqh8AAAAASUVORK5CYII=","height":140,"width":232}}
%---
%[output:07104b66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5993295b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8ff1526f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:796f32f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:610d3fbe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:82feb1d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:04097e2b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:03ece419]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2c35047f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:776256bf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:11816353]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:977393f4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7a7ed79e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:94698f8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8c6b983d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:01aa338a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:68863217]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:45d1a085]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:32252d6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:00ce1993]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4963f04c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:64249436]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:05163949]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:931afc20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:194e8c1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:18144483]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:91a029a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4788142a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2d9bfcc4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:34cfd756]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8941d291]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:137fdab4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:75eac030]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:945a7157]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6d90b171]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6e156428]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:45de0648]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8acdac76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:78df5589]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:48a509d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:79ab9148]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0575dff8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9e6e61dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:44679edb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:32971274]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5b462c57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4439e1a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5c2dab60]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:80b89144]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9a6b778b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:92f26cd4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8bd617f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6c6b7f0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6ab5fc69]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0415e9d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:585f473b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2015c8d7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5aa2591a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:14b9eb2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0e8cf68b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:93b17453]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:07e48f11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:19232560]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:60f6926a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:37418d3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9251ad8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:227c8651]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2a325c28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6d87e43b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9ff13350]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:63268dba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:84a2d3f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:59723cec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:95b9da20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:80b5adde]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:233439d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:22bb79b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:16261198]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:35eaf959]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3e7a3b38]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2e6b76c5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9285ce9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:18f33a24]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4c26c83e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:296de27f]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAM8AAAB9CAYAAAAFgzchAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQd4FFXXfgkhJBAMoSX0ovTeESkaPkGRIlgQRPwoH0VEpSN+KKAgIIqAIKKooILKJyogRbAgRboiHURCCyUQAgmQQAL\/\/948dzOZzOzO7M5uNsnc58kDu3vnzC3nvefcc885N8+dO3fuwC72CNgjYHoE8tjgMT1m9gP2CIgRsMFjM4I9Am6OgA0eNwfOfsweARs8Ng\/YI+DmCNjgcXPg3Hls69at6NGjh3j0zTffRLdu3ZySuXHjBsaMGYMVK1agbt26WLBgAYoUKeL0maNHj6J3796IiYlBoUKF8Omnn6J+\/foZnlHW6dixI6ZMmYKQkBB3upSrn7HB48Pp9zV42LVHHnlEADU0NNTRUxs81ky6DR5rxtEQlaSkJFy4cEHUpQRRMrQWAZ4iXLp0CdevX0e+fPkQERGBgIAAp++6efMmYmNjkZqa6qhXtGhRFCxY0PFZWadAgQLg73ny5DHUB7tS+gjY4PEhNxAMt2\/fFm8kCIwwLOvLo7i8efO6bK3yHcrKymeVddgGV4B0+dJcWsEGj48m\/vLly\/jggw9w6tQpARoj4CGTcy8yefJk0Upn4OH+iNLpv\/\/9L\/j\/qlWrolOnTihZsqR4VgkSZ+BZuXIlDh8+jJMnTzreeddddwn1r3HjxqZHy2p66gZwT8jCfZuviw0ekyNOlef9999Hly5dUK5cOcNPx8XFYebMmXj11VdNSZ7nn39ePKcHnh9++AG\/\/\/67kGj58+dHYGAgmjdvjtWrVyMsLAx8L1XE4cOHIygoSIBWCzzHjh3D2rVrhbpXvXp13HvvvULVO3\/+PM6cOeP47cUXXxTqpiupqabXrFkz8Zy79NQDfejQIXz99ddCpWUh7aeffloYVQh8An78+PEoVqyY4TkyW9EGj8kRI+O98847Ajh16tQRK7yRQiaeNWsWxo0bZwo8Q4YMwbvvvqsLHgKZ+6IqVapg37594l9KiQ8\/\/BAXL17E3XffjX\/++QcjRoxAcHCwsMBpgWfu3Lm4cuUKBgwYIKSdtL6lpKQIYJ44cUIAKDIyEh06dBC0nBU1PdYn4AhO\/pmlp37XhAkTBI3Ro0eL\/syePVu07d\/\/\/rcAUHx8vJDCNniMcKeX63A1279\/P\/bs2SNW7yeeeAI\/\/\/wzateuLczIXPWdFW+Bh6ogGZNm7\/feew9ly5YV6tr8+fOFFOrXrx9mzJgBqo0EQvny5UFAyr0XGXrVqlU4fvy4WLmlKZyAnD59urDUKcvrr78uAClN7lp9pqqmpqc3NkboKZ+lxFm6dKkACvtJIwoLgc\/vKdloql+2bBmeeeYZIX29VWzJY2BkubJRdSIztmnTBqVKlRKTQsvZwoULUa1aNfC8JCvA87\/\/\/Q9nz54VEuOrr74SbXz00UeF5KHa1bVrV8ybNw+09JHp9+7di2nTpmUAD8HFBeC+++5zSBypZnKFV5YlS5YINY6STK+8\/fbbmejp1TVCTw1eLgRUfylF5T6Qc3T16lWwvbQeDh06FLQkerPY4HExutw8L1++HM8++6xY4X\/66Sdcu3ZNgGX79u1CPWjYsCEqVarklNKRI0fwxRdfOPY8ZAA9KxcZgQzAoqW2JScnC1WMAObegvubiRMnYvfu3WKz\/9hjj4nDUaqWZcqUEaswpQqlzUcffSTAIwslDxmNRgmCjRKH37FfixYtEiC5deuW+I7v+\/vvv8X3lBh6hfsiNT11XfbdKD0+KyUODSCU+noSRUqg06dPi4WuZcuWXsOPDR4XQ8uVkSsaJ+yPP\/4QzMqVvmnTpuL\/VB2oKrkq3AOQMceOHSuqvvLKK0LVUwOIDE6jhFSXtMDz22+\/Yd26dUJKVKhQAX\/99RdeeukloZrxIPbJJ5\/E559\/jsqVK6Nw4cJCnaHBgMDasWOHAKQaPNIo8fLLLwvDAsFCkPIdBDOZlfsJ7leoHr722mtOwaOmp+6nGXp8EcHK\/vG9lDh6Cw\/HLyEhQSxSbDdB7Oo8zdXc6f1ug8fJyHH1\/fjjjzFw4ECh+tA40K5dO7E5pYpDtxeah10V7h3IjE899ZTQ1VmMmKr1JA8lB9tBWtzzJCYmij0PVS+qZxKUfCfN1jw0LVGihHDz4XeUTErwULpw483NtVSDSJOGEWkdlPX\/\/PNPfPPNN04lz7BhwzLR0xsjI\/RGjhwpFoIHH3wQFStWdDXc4nceCXBR+fbbb4W09UaxweNkVDkBO3fuRJMmTYTUadGihdCjP\/vsM7EpNXq4qAUeV5MpwaWUPPKZTZs2CTWGIKCrzZo1a4SFrV69emLPIw9VCRT+0VBAEFGK0mBAJlSChyofJWitWrWE1GHhKk\/roFrC0Liwa9cuh9VQqx80Eavp6fXXCD2qgfyjFW3SpEmuhk78Tqsmzfxc\/ChNvVFs8GiMKpmRKgsZiMDh3oYrMlUNqi2UNtwHyEKGI\/M6K1R1SE+qbVI90lPbpk6dKhhfCzx8D6UfJSPpcf9DqREeHi4ArbSkSfo0btCsTcNC8eLFM4CHVkRKJe7jaD1kIW3um6gOshCAVDsp6Sj1qLbqFRol1PTUdc3Q47OUngTlnDlzBCn2kfsb0mHhosH5kf0dPHiwqK\/sq9UAssGjMaJUVcgctCrxTOPLL7\/UVRcIJqoHXL2dFdahZUyu5HJjrvUMJQcn3Rl4uHf5\/vvvhTmWey4aMwhqPe8Bvp+HimQoLgRqgL3xxhtif\/Tcc885DlLJnPyOhSZgvo\/GEu59CFRnRU1PXdcsPTV4uBBxMaBhgIWGkUGDBjnaZYPH6qXCID1OPA\/YuIJSCkkm1jrL4eaem2Pq5SwEh9q8y+9p+qW0cLbR1moe1RW5+Vb+zvcePHgQ9DDgXsqVtY8Ghm3btol+aRVusnkISjWQtGrUqCFce2gUIdCpqhHQBBQ37P379xervV6xmp4aPPzMsaQhhoWg57hLSWODxyCzW11Ngod0pXeynl8Zf+fGmuA5cOCAMGtLfytlu6wGD2nTOEDjAc+ZaGFzVmiupiRVWtqU9SmxSO+XX35BdHS0iAciONg\/urrQ5EvvBYJo48aNaNu2Lehyowcgq+lpgUepxqnBYoPHalQYpEcLEFd2I4VSicy1fv164Q5D6eQr8LB9NGR89913Yi9D65ueEYMbbfqrRUVFGemWbh2qe1SV6MHA0309MBp9iVF6avBQGlK9lnsggoWqszwfs8FjdAayuB7PT+hpIM9E1OZdqbbResXV0kx54YUXhNVLr1Ci0RxLdaV9+\/aakoDnVLT40aWGUsrTQqlCCUtDCmnWrFnTI5JG6KnVNKqFNLrQEMNCyxrPxqhSaqlxHjVQ52HbYOCNUdWgyQ0uV2uqTq48kuXjZCqe+qv9y9TqFjfx0oqnRZt0qGZR6ljlskKatPLRtUcysCdD6YoeJQ3DDqhSso+sz39pBGGh2iq\/4788PKYGICWRJ23Te9YGjzdG1aaZK0bABk+umGa7kxyBG4ffRUjVtHMrK4oNHitG0abhtyOQfOp\/uHE4LZjw9vXTCKn6omUAynLw0HJDXym6mvCgy1nhaTI3itwUOjtjcDWTOZFOTuwT59HTfl37YyQIIFnyl30cBeu\/5YpFDP2e5eCR6ZgWL14szg2cFZ5D0KOZh3euIhlzGx17bLRnnMAhgGQhcAggK4oNHhfhxNkFhDZ49GeKe51bF7chX7GmlqlsfJsNHhs8Gbgup4LQCkmjpmGDxwaPDR43kaUJHh66MdSW+xAe6vGgiXHxdAZUunjzMGrz5s3C4ZFevqzXp08f9OrVK0OGSmdtyw57Hho1pPeusi\/czNIlh57HnuzBrKBjBQ25QfeXPrnTHhqdXBme3MRKpscygYduEPTkpXctXS8aNGggvIvpfsKYeLqKlC5dWhBiUBZ9iOi+T89eOkayXqtWrYSHq5HwV38HD0FDp096JNvF\/0eAvPjWW2\/5BECZwMMgJrp60OWidevWjtFiTrC+ffuiZ8+ewo+IwVIEDqWNMpH4hg0bQH8sfkdfK1fF38Ej28cJkYuGqz7Zv2fNCHCBY\/iGEcutFS3MAB56uFKyUBVjxKDyOguqcvQVotcufYyYv4zqGWP7lZ66jH1nDDufpZerDOvVa2x2AY+vJsSKSc2tNMzwkhVjZNhgQM9cShxGEBI8jCqkdPrkk09EcgZZuA+i1OEeyMh9MmY6nBWWIDPts2JCbBruj4Cv58oweKiOETzU\/ylxGGO\/ZcsWkZlEHSfO35jYgTHwrrKdmOmwDR73GSs3PGmGl6wYD0PgYcI+xocTJFTrmMaIAGFjtaQLM1fSAqeWSloNlh2mkaJz586O1ExadQmec+fOiSAsT24yM0OH2XO4WORUtU15+5xyzL3VX\/INtwQszLxqhEeMMrrkJWbM0fJWoUuXJ25d6na4BA\/TL3EPw6ArJlxgKC6L1eAhTWblJKPqFWa0oTWQIHaVG9rZgJuhw70d++8tZjLKGN6oxzmlEYiMxrxtskgmZHSq8ntP20CeYXi3vMbRzCJr5N2y3QyLZw47deGRgqvEJUbeI+s4BQ8TTFAiyCyRvHpCFmfgcUdtozWLaWtlUkCtTnCVZNYV1vHkXMUMHU4Iz65yInjIvMwMpKc96P1mhsHUgFSOo5R6lEBWgFSCh2eUjRo1ytRMn0gebvp5hsMVlxF5RLI6pSwH1jYYuMNGrp8hUysT9SklgLxPlAuNXMElEzI5B9UgFiZl5JEDE75TerKoJYmVK7+yzQQD381Utwz\/pmTTA6ozALseqYw1\/GLPozz8ZOIIrcRxbKhtqjY7va7rU2rzrE3uBbRUK8kk8kZtybhyVVde2Cu\/06KjrGfkdm691qvbrPV+tcomaVkJ4CwHjzQOMG8Xz2n0dETG5POQlOEBSm+CnHpIalRtm7r2OE7GJaFckWCMbmcsr7JkJMl0zEijvGZeiymk2kx1h7mmme1Tqj7O6LCucpOuvKFbtsPotfWsr\/cuNaD1wMP3q9vkeonRrqEcp7rhW3H7+hkEFChtqSe18s0Z9jzcSBMwzGnMyaDKpi7y8iQefvJqC6ZjZQ5n3iLA6ydymnuOmdVsyY5zGLzkoGPIRrerYApAequwlBp0gZKgkt9RJVMzu2RoqTK5AqdSCuipi2R+8oMS1HxOj\/nV4+ZL8Hw69THUyp8eAGdl9KgueGjJ4k1i9GXTKwSV1LW1HEN5Dww32MxlZqSYYU5\/P+chcAggWbo3jsSc7ulGFlfjod7rqOurVStZX\/09wTNq1CiREFF5gK0FQq02Kc3XUuLqgUcP8GoA+1Jt+3BUadQvlZaGl8XK6FFd8LiaXG\/8npPAo5Y8BA4BZLSY2Twr9xXq8xIjkoc3JdBMrZRmynaq1TEyP9Ng8XIvSjtpfDAqeXxpMFBLHiujR02Bh9n16ehJ9Ux98JQbQhLMgJsDyz3Ppr\/j0eKewqZUNj7riumlGiYlA+9JpcVT3vQmNQKjex5nxw1akkMaMugMTOBxj8X\/07Kn3qep9zxaIPOWqZrSknseb0SPGgYPD0g5QOy41oY5N4QkmAWPUSmjV09tuZIMxvrKw0XuTdTWNcnAkvH5jDQOyO+UhgWpxrGe8qxHyzLHdrFIowTBwQSE\/KxnbVO+Xwl4+S4rLW18l6\/nSjcYjilc6d7NTI0savBIa1tuCUkwam3zFDx83tk5j97pv3L\/wxgsSgNKBzIqT\/VZ9MzRSpcZ2X51f9V7HiV41G2W5zw8Y1Lvu3zhnuOrudIEj+zgPffcI2744iGpukH2OY8VMPEODT31z5O3qSWP+rOaNvmDdYx41nvSLuWzfiF56NTJeByGXtPyxohSNXhsDwOrptx6Ot4Cj3QEZoud7Xn0fOas72lGikrwVD25Fimx0QgsXgHBdcO8cubj0jFUD82+CEk4P3ELbp64gqDyYYh4tbm4P8bXedt8vZpZwWDeAo\/S2qZUAbUOWj3xWHB3DGQ75vdvh0q75gkyQXeHokDzYg6SVp75eAQeuRIpI06l\/ms2JCHTnmrRPpzqu8bR6YhxzRE2qoENHnc5Kxc8J8Ezs0V+1LxxSPSYwCGAZLHyzMdvwKOO5zk7YB0SFh92dLpQj6oIn9nSjufJBSBwt4sSPHP7tEGVPWnXx6slj5VnPl4BjzshCeyoMp7n1rITSBqz0zGOwVMa4fYjkU7jeT7YFo+zCSkoWSgQA5qmXUSrVex4HnfZ07+fk+DpM3kWelxdj+DzCxEQGphB8viF2ma1wUArnmfeE6tQZG8s4moXx8Cl7eEsDmfx9rMY8V20Y3aHtymD4W3KZpptuY+6XvgOyk+KchkXxAnJqfE8\/g0F862T4Hl82kd4pU4MUg+Ny0TEL9Q2b5uqtZwsX2xdUnfPY8Sv7LLGPoqGCGclOxoMzLNdznhCztX0jz9Hu7BvM9yOIHvoF2qbt0MStMDwdpeKuuAx4ld2qu9qXF6038Ep4b1qYvGT1ZyGENjgyT7AknM1p0tVNK4XiICwfY7GBxZt5vtE786Yx5shCfQRm7o2XQ2jg6Uz8HCUWn8djZik2ygVHIANT2YOp1BLno0D6+HU9c9RMvUczuaNRPiTr2XyR7PBk\/3A83aVs\/inWTiu1imCEriNyqs642aB5xBWNggtR0ZY1iG3DQZsgTdDEih5Cq2f6mDshH+Ndgqe5edu4tWD1x0DM7BCMAZWDM4wULxq4uK8tY7v\/qxZFDXXfuH4zHfUHTAlwzM5HTzezp6jPgPy5vmPfNejbW9ie5t0kBA8lVc9Kua15YgIywDkEjyWwVSHkB5z7l82DyFLBjmeOl11GJr8d1IGte3y1xMcp8gzao8EASRL+7wxmNyqhuOz+pIj+UPSnngk\/RUvPha6\/98oPjgtB4AsORk83s6eo\/ak1vP2torH5FxV6xmEq7ULOsiW2doCdT7vKz7X6RaODrMyG5LcaYPfgid2Tm8k\/Pqpo0\/Hz3fC7Qdn4Z6et0Xo962tX4J1ZNnY8mXMqD3K8XnIz4PRtXZamiy6aeQJ24eAuy5mGqObxxJxfUva98E170ep8b9AgpLf7UkIRv\/5a+3sOSa5Sw+YZmKWTL7S4VWtBk+Ro9XQbGZaai2\/AY9V8Ty\/LuqOLg83Q7kq6Td37flgjFDbZLkQ3winKi5Go5cDBXgSFgwS4Pqq8WhcKFQOgSUqYF2R5qibuhP1Unai246pSDlbWACHhSfN35QvjgsIEHpw98AbWBg0ADHXS6DYuX\/wzKU5QgIRQAsLNBU0SyScRLWfZ2L4kZI+BU92yp4j1T5GrMpAOYY9MCKZiTL1QsHV35sFilZ9KXlKdW0NNP0rQxUCqOjRqnisUG\/\/kDxWxPNc2zUdyWfmODr6RalPEFuwCQI2rEeDj1ejVvk0HyWWRW2+x9WGzVDu1gU8+eMQrE4thdlR6c+yzrPJ8\/DszbSMlFIl4ynz701LY2ZKuiivEtoDu4PTRLl8rvvBt\/FDXNsMNAme7V8s9Bl4slv2HK2QbY6nWmWT4+xN1U2CJ7HFKITXO4HiZX\/MhLGOcV0xesArVmDV\/WsVrYjnoeTgZatk7iUpITiT\/2GsKjAOVw6fRcqNmyh97CyaRExE+dRoxJ5qiy9qDUVgSBDCqpZEzTObEHtXOTxwaDFKJJzChUJl8VXjMYi8HYPm57ah2S9H0CDqc6ScT0JgRDBmB4fhp9T8YtBu5W+H9t\/9C6XClyGgwGXx7IaGvdDu1nKcCyiFr\/9pLN7PdxXEGVyf1t8n4MmO2XO0gvWcgcdoHgV3uFsJntDqV1Hqnq8ykWm+szCGF3xdOBp7Wtze81hxSMo9y62477CmZG18HHYNnY6UQNGkIkg8lRaAR\/PxtzVLIl9wHDrsShaWN5bgoqEILJAGBIJHln2lWghACTAEPYAaBaaiULEL4ufzdwKw73Y+8f\/Gv3fBgHW\/IajaukzPRtyOwbEzabRZDgcVx5r1mw2DR2nEoOnbTMmO2XP0Qqn1JI8vwNOjdgSC7g7Ej\/emW1\/lPFS7kIBJ1zsgYuQ7ZqZGs67b4LHCPYfgWRywASvzVEajhMP4z\/YTmRr59T01EJM3Ei8d\/tlUZ7+pXR\/fV8+j+UyfVXlxf+IOQ\/RoMDC656EkVRoxwp8YL86OjJbsmD3HLHh8obbxnKduoSSx8PJPXXoVb4f+9082Oi269dwGjxXxPN9tGI3pF9aLxhE4LaPTJI6ybKxQVHzU+s1Z7\/nch03Ka1bRe5dWZTPgUVsItUzfztpsxhLlL9lz9MBjNJmJxxysICDVNgkePR54uEJHvNJ4vMev9gg8nsbzTNoxHqujV\/gcPA32FDYsycyARy15eGZEABktRhnOn7LnMOEibwxUJ2vXU8\/MLBBGx03WU4OHi6dcfJW0cgR4VkWvwOQdaSsAJYuW2vZu1Sg0TN5jqeSJj20k3tXxWrq3gd5EmQEPaXDPc+PArwipcb8plU2+P7tlz2FaMi3wsD96ObTVaarMgkSvvgQP9zzh9U9oqmx8dmzj8WhfoaPHr\/WK5DETzzNq2VhsSU1j4n67ktDq2AFHp+bf9SxWNsyPqJi\/M0iKXfnT7l4pmXoepVLSM3QqR4Og2103HiGXiuFG0cyHo7TecR\/VLm674zHSVdNcfaMC3j6Qx7DBwOMZyWbZc5yBh2ORFe45deq2xcUe2nva+sUbYvb9862YJvdN1VYYDNiD3XNWY\/bOHbhU+TDiy\/fBw2fWoeaZzdgdXBeLyj+OthEv4MqFxii\/MxUNk\/\/Ervz1sLDSYGFxK3Y+EU\/vn4HSQQfE9yyyzpGQcgip\/Tuq\/10Za5roXwPfZf9Z8I\/AGV9kNMoVi0f78z8g8sppQfPjW80QummaT8Hj6cx6K4cB2yXztrnKnuNpH9x5XgK1Q7U62N37lCYJq1Q2Endb8lhhqmYDDiyfi9PLTyHmSG1Ety6I76LSr6\/n78\/EvYfHlv2FV5N6Yl2DkgI0POd5fE8yNuS7gxa7DmBrMHAmMA9qHb2Ec3cFod7pBAQ0+BZNo\/ahTd5kjLkcggMFQhyDye\/oZbAnAah6PF6Ah2VFwXbIG1YXPdp8jVnRnbEtviqwfw8O7\/jRBo\/iGk2Olcyeo3V9oTuMb8UzEjxVukYisWmy\/4LHqngextikJn2C\/LVikLyvFMY0noBDddPNi602bUO\/iXuxsW1lRI8Mwdn8DdAiMhS9A1di6KZQbIls5hikQQtX4v4\/d+FIjQuIfngHuiUlYO6NbtgfF4uzNdMSQrC8GHhNgIo+bfRtkyX672dQuUF3FBt1WKRq5QHrts1JGPztYRs8qlzV3vSOdhdIDvBENUJiu2OaZEZVegidGk5y9xUZnnNb8pCKFfE86hibP96IwowmlRyN\/M\/0DWj541H8HZWIVuMvIiFyuPBtSz04TnhRTw2e4Kgbcmg3aiZsQZngi+gXs0AAY02FAVhdYQDqV16MCwlrUTtPivBro+sO49vvXH0Aib\/uRurFSrh5qC3KLngI4b1qOWjmZK9qSzjIj4g4Ikm7F8ePDW44DsXZRO59qcK\/+XQIGE1qRfEIPFbF8yx\/ZAb61u+K0FZlhdvEvONJ2LLnEu5efhhdPtst+rlo4FaUf7IiOpR5ToAnT+xK4dpD584\/AxshOu4e7N163DEm\/a8sRP+rCzOEGcQt64KUS1uRci5JSJXQqIdQpOu3YF6DxN9OOd6vHFgbPFawmW9oyLka168oQqsnZPBlrPNZX5TZ1gL3Db6F1q82tKRBHoHHihY4Y841Q2YjccMpHK1xHj88tg8PRD6IPmUGCPDwQl8Gt8lM+CN2t8twN06H\/zdDj4+bKuJz5FmLO+cwNnismGXf0JBz1fbx5pjcYrPwl9wSWxonisSBMT0hcUV961WdlVeMKA9ROfxq8CinRJ3DYErJbXiifHKmsxaewyT8tR4p5RqiVK8phrLnaKUb9g072G8xMwISPB2HvoNqlefi3PUYJFwsge3hVxxkGJIw9KHnzZDVretU8mT1FSPLfhiOd67\/6mj8cwHNcG+DEQ7Jo+6V0btxWO\/AyYuoUqowXnmkstOBtCWPJXzmEyJyrlqP7IKjxdL5Rvlyr5uqr127Bn+4YoRq1vxfx+Jg8VBUj01Er2avIrFqO13wyEFydqmu2XtDbfD4hO8teYmcq1pDqyAuMlaTptetbf50xYjS3SWk02inuaq5B+INyG9tvCXOaVjUl+quePMl\/HP4kAh3mB\/2rLj20Nm9oTZ4LOFrnxCRc9XloZbY9kB62ikeSzAkhZbWfo0nijtKrSiaaps\/XTGy8a3zuHLqpkgb1HhImC541Ak+CB7+KcGhNhjQ\/afuwClO7w3N6eDJidlz3rjQGSc6xOFI9QuoX6Ih\/hNVzCtXLLq0tukxjxUhCUS\/M+b866vLWPlCuptF4+fDHAlAaG1TFpqtCSBZvjl3H0Yf6iukirxUVx0ycLlONzQc92Wu3fPk1Ow5BE+t5NJiXpnYct2wkzh7\/SxKFiiJPjUHWCF0BA2PwONpSIIr8BA4BJAslTuHOBKAqMGjljzLbw3H+dBHMiQxXPPRu6iydqiD3pF2M\/BQv5dyLXichQd4GjqQldlzlODZOy4F71df6pjjPjX6WwYgvwGP+ooR9nbNsLM48E2Co+OVOuRDvZF5EBERgZCQdF81WSHt3Gcr8hVrhpCqmUGx9I9Y7P9onMN5tMqzr6FHk8yRhko07dy5E7169fKpe46dPcc94eCI57l\/OKocDkVo67L4tOPPWHd6tYOgpda2ihUr3iHlAQMGODxmlU13prbpSR4zIQlKl3XlFSNsw\/EfbmHr60mO5tQflQehTeJRvHhx5M+fnmfAzFDzGpLtJxNRrfBtvNAqwiWdPXv2YNiwYT4Dj509x8xsZqwreYl36DJIj2XjxV8wPzo9w5JVsTxCbXMXPFaFJMgOa10xwgbumH0FJ7YkonzzUNTqF4Tz588jMjLS5eGmsylwdlWJ+jm2z1dXjNjZc9wggx9WAAAGxklEQVQHjnILsGjRIjRq1MhBbPGxhfgjdhcYy+MXex6rQhLMWLOyw52k0lweUKC0purojD38OXuOXrvNJgDxRfacLL1K3ojaZlVIQk4Cj9poYfYWMn\/OnmMVeHyRPcfvwcPBtCIkISeBR20uN3sLmRkLl6+z5yjPg1q3TgtY5D5ZLwGI0WQmnilq2nuebAEeq0ISjDpe+rvappY8Zm8hM8pwWZE9h8DesmULpkyZAhpR5Jz5Y\/YcvwGPlSuDFq2cJHnYP2WYhJa53NV4+mP2HAnW7t27g2HXys964GE\/syp7jg0eDS7zd8njChhGf3d2ziMXG\/XRgnyG4dENGjRA7969RZ6BBQsWICYmRrxaL3Ra+jIq26dkQLnJZ\/IPM+BRWsAkbW+Gb5tZiI3OhbN6Lg9JrXiJMxpmOpxbwOPpmFudPceZ5PHHBCC25NHgoOjoaCxcuBA8TK1QIfOdo0aZzgwdM+A2+n5v17MaPFIF09rz2ODx9mw6oS+Zk+45TZs2ddqSM2fOYOTIkTBS1xkhM3RkXV+tZlZMhTfAo7S20ViQkJAAuQeyos1W0HDFS2XKlAH\/rCpZrradPn1aAGLbNv3EhFZ11hM62Qk8nvTTyLME56hRozBt2jTwRjh\/KUpXL602cdHln1Uly8HDjhBA\/PPHQlDPnDnTZ75t\/jgGbJPaqOCPi4nS1at06bSQBGXJcZLHX5lFtis77nn8fUy91T5fz5VfSB5vDaYVdH09IVa0ObfS8PVc2eBxwWm+npDcyvhW9NvXc2WDxwaPFXzrFzRs8PjFNKQ3wtcT4mfdz1bN8fVc2ZLHljzZCiDOGmuDx8+m0tcT4uvuezv1lF5\/lF7aWvko3BkHX8+VLXlyseTxduopvaGVTN6xY0cR4mCDxwUTMoUvY8t5uEaXl6JFi+LRRx9F\/\/79RUIPWbRihJhDgBlsChYsqPuWTZs24b333sPcuXNRpEiRDPXcpUkivl7N3Flx3X3Gm6mn9NqkPGy1wWNg5mJjY4VbBN066BdFl\/m9e\/cKJ89y5cph1qxZkCfCBMHgwYOFn9tTTz2FAwcOiHqtWrXChAkTEBoamgkYP\/30E15++WXht0QXfDV4zNJUviArwJOdUk\/JsVK2mdfKjx07FpMnT8b06dNFGAOLMlZp2bJlIlTCljwuALRixQoxmJQMMoSXj+zbt0\/EnPTs2RPPP\/884uPjBXAolRj3IYGyYcMG8NZlfte+fXvH2+jS8\/7772Pp0qVISUkRIcFq8MhcC0Zpqrvia\/Bkt9RTalDQ100ZIq7nxsN+2uBxAZzbt28LybJ582Z88MEHGaQCVbkxY8YgICDAEd5L9WzevHmIiopyUE5MTBS50yhRJk6ciKCgIEh9nSHBXbp0Ad\/DUAM1eDzN8mMWPMrc2i1HRhiQy+lVsmPqKb02q6NI1QNhg8cUa2SufPXqVSFxwsPDBXi+\/\/57IZ0++eSTDN663LNQ6uzYscMBDoKHud46d+6MJk2aiP+T0dXg8TS\/nBnwqHNrtxwRATMA8ufUU2R2xlB169Ytw0RyfEaMGJFpzlyNmw0eD8FDdYzgYUgCJY4nyeP5rBZ4PKFp1mCgzq1dp1s4Oswqa3iU\/Dn1lB549ADvKqbIBo9htshc8ciRIxg0aJCwtFGtK1GihACPFgD4tN4kScrOwOMuTbPgUUseAocAMlr8OfUUx5fWThpnqCrLHAq25DE6uxbV43WN3MNQ9eKGv0qVKg5rjLuM7g\/gYSe455Hpgc2obHzWn1NPKQ0ZNL7Q0MOEIPw\/k40MGTIkg0pn73ksAouSzMGDB4XZOjU1Vexvqlev7vjZmeRxlTzeHfC4omlW8lgxXP6Yekpa1PgvASM1ARpo+FndZrkIsJ56\/6rUFGxrm0GO4aaf5y2UONx4MpN92bIZ9wOebO71wOMJzawAj2RMnlvJokwzJTfivkw9JcGjNBhQskjwqNssz3k+\/PBD3VBte89jEDispjyonDRpUgavAknGE7OyHng8oZlV4DExrJpVXW3W3aHP8VVKHvVnNU2OO+toHVq7834jz7iy8BmhYaaOT3zbpHGgRo0a4pyGpmmt4knyeD3weELTBk\/6LCnHl9862\/Po+cyZYUx36uY48CQnJwvALFmyBPRl0sq3RvWtU6dO4vDT3eTxzvZL7tK0wZMRPEprmzLzp2RaJcN7MzOoHrByHHjo19avXz\/hy6ZXlA6C7iaPdwYed2lmV\/C4s2rnhGdyHHiy+6T4ekKy+3hlZft9PVc+2fNk5YB6+m5fT4in7c3Nz\/t6rmzwuOA2X09IbmZ+T\/vu67n6P5PmXcSnnzujAAAAAElFTkSuQmCC","height":140,"width":232}}
%---
%[output:6d0fa995]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAM8AAAB9CAYAAAAFgzchAAAAAXNSR0IArs4c6QAAGrVJREFUeF7tXXlwVsWy71gsSZEAEZGdijEgmxcwKLx3QSjkiYiAYrFcsVhEWYwsAhp2kaUAWe4DVJYCQf4QEMESChBRH3qBKw8pXpBNiBoRBWRTEwQkyKvf4Hw538nZ5nz9Lfkyp4oqzTfT093TvzMzfXq6E27evHmT9KM1oDWgrIEEDR5lnekOWgNCAxo82hC0BnxqQIPHp+J0N60BDR5tA1oDPjWgweNTcV67ffHFF\/TUU0+J5jNnzqRevXo5dr1y5QqNHTuWNm\/eTE2bNqUVK1bQ7bff7tjnxIkTNGDAAPrpp58oJSWFVq1aRc2bNw\/qY2zTpUsXmjVrFiUlJXkVQ7ez0IAGT5jNItLggTidO3cWQE1OTg5Ip8HDP9EaPPw6DaIYDfCUKVNGrCzdu3fX4Anj\/GrwhFG5IB0N8GDchg0b0rJly6hWrVpCQr3y8E+0Bg+\/TqO+8kgGsrKyaMSIEYSVSIOHf6I1ePh1GlXwpKenE5wOp0+fDnIeaPDwT7QGD79OowoeeNJat25NEyZMoMLCQurQoQPNmzePzpw5E\/DIaW8bz6Rr8PDo0ZZKpM88AMbEiRNp6tSptGXLFrFlg\/Pg3nvv1eBhnmsNHmaFmslFAzwAyzfffEMDBw6kc+fOCefByy+\/LFYjfAvSKw\/PpGvw8OgxplYegKds2bL0xhtv0IIFCwRv999\/P3377bd04cIFDR6mOdfgYVKkHZlorTyIHvjxxx9p0KBBdPTo0SD29MrDM+kaPDx6jLmVR4bebN26lUaOHCmcB\/LR4OGZdA0eHj16Ao\/TUDVr1qSVK1dS7dq1A7FtbqwNHjyYsrOzHb\/hGGPlNHjcNKr2uwaPmr6UWxu3bdEAD8Y8dOhQwHmA\/9crj\/I0WnbQ4OHRY8yuPGAMN+2XLFlCc+bMEXxq8PBMugYPjx41lVKoAQ2eUjjpWmQeDWjw8OhRUymFGtDgKYWTrkXm0YAGD48eNZVSqAENnlI46VpkHg1o8PDoUVMphRrQ4CmFk65F5tFARMBz6tQp2rBhAz355JMi\/MTtQRxWfn6+uAmJ+yh+Hy46GJ+LFhcdTp44aZnlu\/Tuq1R4Lo\/KVE2j1J6vKE0lp66UBvbYWBk8n376KQ0ZMoRWr15NrVq18jSMDFF55513PPW5evWquEZco0YNSkxM9DSGVSMuOqDNRYuLDidPnLSM8l3\/Yi2de2NAYGpSe0xRAhCnrnwbkUNHJfDk5eWJEPfc3FzyCgSMrcFTNAOcBhGLtIw85a8YSvk7VwWET2nXn6pmrfRsx5zyeR5UoaFn8Fy7dk1c7d20aRNdvnxZg0dBycamnAYRi7ScVh4ABwDy+nDK53VMlXaewbNx40ZxM7Fv37706quvavCoaNnQltMgYpGWmSecea4c2UlJjdopbdk4t5I+p8q1myfwHD9+nIYOHUqjR48WeZORe1lv21x1a9kAxnXgwAHCSh7KeU4e8s+fP0+VK1eOGVo45EeTJzikvDil\/M1ecC9X8BQUFNArr7wiJmfy5Mli4jV4\/DsxcG4cN24c7d27l2P+NA2TBlq2bCmuXkQCQI7gwT2Qd999V2TqR+rWtLQ05cO\/dhgEz650nmCCZSpcjQAeDeCFhIQnKruiUEZ2BA9uIOKq77Rp06h9+\/ZiHFXPmbEPUr9269aNqlev7sgztjZI0letWrWQymBw0ZH7bw6evvzyS3FujNQEh2IcJa2vtM233nrL8pMIvhmG8t3QrA9b8Fy6dEnkOc7IyKDx48cHBg0FPBi8X79+wnicHpwHkG+satWqVL58ed9zyEUHDHDRysnJoVGjRmnw+J5V+47SNufPny9qG5kfnA1TU1PZRrYFj5e79zJpRb169RwZMm5VMjMzXVceJK04e\/asaBfKoZqLDoTjogVdPPPMMxo8bCZcREjaGT7gt2jRotgIEVt54DE5duxYMQa+\/vprmjFjhsg+2bhxY2rUqBFVrFjRE3i8blW4XLBcdOS2jSPqwc\/K7dXOZs+eTUuXLrVsbq5Kh7bIHhpqhTinMcGI1ZxbvZjtquCZs\/84VcsLp26tlOrqbTN38sOgah8uo+eiU5LAA12bSzHKCgnDhg0LlHXkBI\/VmFbGdvHixUAWH6TZMu5YJAiNIJfAwQ4HKbbwoJ3deKp25vWlZNdOg8ejBrmAGM4JdjIs82\/RAI8Tf1bAWLduHS1atEjks5NAkwDs3bt3sfqu4dStXnk8AsWqWTyAx7hNk+BB5Dq2VniszrDmbZk5bZUbIKQurVY\/s57RBpEszz33nPgYbwdwu7\/HPHj82J+qUFyGykUnWtu22du\/o5MXr1Ld2xMpu+Ndrqq3M2QYJaokvPbaa4E3uASFEQzmN72Znuo2ysiw1SriJJDVWLK9nZyqduaqUJcGyts2PwOqCsVl9Fx0ogGeNfvOUNaaogTt2R3TXAHkdHg3rypWBmg02OHDh4u0v8bzhtXcO41pdBZ4XaHkGE7gARDXrl1b7Gynamd+bNnYR4PHowa5gOh1ggEcAEg+\/7i\/Or3xj4aO3NoZqDwnoLN0JthtfYyGuWPHDhFK5PRJwiso7NpZgQ9Og65du9qCV4PHwQy4DJWLTiysPAAOAOT0OBkyDA5AkKsB2uKRXixJ12yYxlqmVucir+Dxsm0zOgOcwKO3bRo8rrdqcebZlfsLtc6o7LplgzpV3LheVh4c2o2PlavZK3i8OAzMnrQS7zBAWMrmzZsJ8UL4aArvTNu2bWnMmDFUt25dT5ser1sVSYxrxeCiE42Vx5NiTY3cVh6j29ftzGNekeRQ5rn0Ch43cON3M3istmclxlWNQ9ukSZMEeHANoU2bNiK3wPLly+n69esi836TJk1c51mDp0hFqrpwVa6hgZO3bcCAAaI6gvFjI6IRZI0fadyYa2OdIPzdGIVgHkMFPMZoAXPkgdxWGs9XEih4SUseVFZXFd35aevoMNi\/fz\/1799fVFfu2bMnJSQkiDFkuT6E5yDi2i14U9VguFYMLjolaeVRDc8xfucxh75YFcYyt1EBjzRQ8znK6iwl25bY8Jz33ntP3OPB9WtjKAXu+cAjsm\/fvmLuQisEa\/BEZuXx8\/aMpz6qdhaq7L5c1X\/88Ye4VYogUXMclQaP85REeoJDNZCS1D\/SuvUFHoBm4MCB1KlTJ7GHdrtgpCoU13aLi05J2baVJEMPB6+qdhYqD8rg+fnnnwlfn3FZbfHixVS\/fn1XHqRQ+iYpkb5J6mouvhtIO4v6TVIrCXC7FNu1zz\/\/XJyDWrdu7UlQ4\/2NeLpJenPHAqJLp4hSa1PCf43wpAt9k9STmnw1knYW9ZukZu7hYcMXakz+woUL6cEHHwx439wklUIh6UW83CRFJkxkxJRPyhMTKKX7RDdViLso+iapq5p8NZB2FvWbpEbujx49KvIZ3Lhxg+bOnUvNmzdXEk51L8p1VuGiY3XmQQ5mP6lkVXWhpOhS3jjSunU988iEh7hqPW\/ePEpPT1eeIlWhuIyei44VeAAcYxJzr6lkVXWhrOxS3CHSunUEj3QOYMXBPrJOnTq+pkZVKC6jR6rXS999RSl1GlDVp6Yr824sj5HUNbtY5QY\/qWRVdaHCtNP1gJKaw8Aov10wq2wTTt1azYMtePAhFOE3OKcgLKdZs2bF+leqVElEHiQnJzvOsapQHOAxrwyq5S3M\/ZO6ZNPVNoNDLnuiqgtV8IB+vOQwMMouw3eM4URm3YRTt0rgQSUEXIbasmWL7fw5ZTIxdlIVigM8fs8kkm9z\/zKtetGN7jOVwWMu7qSqCw7wgIZVTBpX9hwrwFrx7RbK4xSsipg7PCUCPCqT5tbWyWCsKodxgMfvmUTKYu6fMnAxFdzTUQk8VqvfxyltKTvLW9421apqTsZpDu+X\/x\/LOQxkXNvJkyfp9ddfF9fInW62hvPFpLTyuAFC5Xc7oey2VhzgAX\/n3plI+TkfU2pmJ1HeQsUYcYszZ8lYyrz2f7S\/fDPCmafv3xKVwGO1+k3\/sxtt\/qd7xlA\/206nqOqSlsPAbF9O17Jj7syjAg63tnbgsdta+QWPGRxGOqol\/nANOuXj2VTjxq2r0EjC0bRBGtXsO8tzFlOr1e\/DK2meVh4\/204nh0FJy2EQKnjuObnddy1UN3uWv7u6qr0ScmpnB57DG5dQ0pqiD43HO\/6THnl2pGv9T6sVxOpNbfSQqZb4y1k6VoDH\/Kg6Ht58aQhVPf0FnavRip6fs4R2bV5HfUcUXYe205ufbWc85TDo1atXkGpUVp5lgzpS+v4lgf6qc+bV5l3BA6\/b7t27RfI5XEGoUqWK+EKOZO0VKlTwNI4deOTbXW6N8jtkiyQXTiuP3XbG6k2Nc4pMkasCHvMYRiFV6mpaZcDJPPAaDVq23VOualVXuNtN0pKUwyAU8CxoXZ4aXylKFa0yZ54M+q9GruDZtWsXZWVlEYoGIUvjkSNH6O233xbhOSiv6Oamxjh24JHGNejXt8X2KP2eBtRl3H9bgkeuNlcO7xTLsXwSG7ejmlP+R3ztN360ROnypMbtKL8gX3znoUs\/eo4IMAPRCjxezk9WGXCe+H6hZ\/CoTCTaqtyyLAk5DIzyq6w8ZvDAFgAg1VL2bvp3BA8CQQEcrDb4yCaB8tlnn4nIavzt0UcfdRvDFjwwwNVb\/01dLm8PAkPZqmmUX7aSOF9c2TRbgMUYCmO3lQI9tDOCy445KPOtjMn0SN5Sqll4VjSDkpdV6ksd\/nccpR5cZ9kVYE1p1y8IqIfvGyHORJKO7Hjy4hVan3erihxeDmhDR7bT6OM1PK08roo1NXBbeUpaDgO\/4JlX\/zQ1TblaTH3cK5AjeLBiYHuGj6WyuBU4QqlF1JhBdhVUyC5XrpzjPIPOh6M7UbcHWwSVu3MCBAjCUK8e3ulqQ3KVQUPzymTXGX3e\/\/2eIOCi7f7yTSnzWo7rmMYGP5WpTjULi3KsOXXOyU+MOHjkteeSlsOAGzygxwkgR\/AgKyP86+aM9qrXsEWYzPopSgYZr43DDZ54y2Eg7UBl22a38khaXA4ER\/BgG7Bnzx6RLQdV2owPftu6dSutWrWK7rrLPo+y0+E7XgESrZWnNOrTKLM8W7uBh2v1cQWPXeiFlwyQEMzp8F0aJzucK09p1Gdcg0evPMEmrcETPoh7XXm8Xh9x49T3yuN12wYG8MHx+L8+LMZLjRtnHQ\/acBjgKfw5z9KDhkO\/F8+amxLk73I8Y3u7sY1trPiA4+HvGanFeNfg8Tob6u3cwIP5TWrUjs1lHRGHAYTqPnUtZXbuV0wjE79+XnyBlw8MkW6vTYV1M4NCYazc0Dj44blyZKfwksEVnJqzLhhQ6S0pMSlJKE22hSv8uuFbEf7upFQx9sGPqfD6dREnJ+kY+8kPmvgbeDly3wiRX9oYI4ffLpVNpVmffS9u5uLbmX74NIBUAS+99BJN7NGWHknKK\/bi4nIUSI4j5qpGul6r4q5WYShlW\/UudvEMDHuN95KGnJDekuUODsb2G2+HvsaE7WmJl8UElzn\/NZ\/VaEoBDeCFhDtotWvX9mwvftXn6SNpjRo1gqIJuD6SSqbNYSh2hqoa7xWKwZsVykULdKa+8y86+sMFUf3g73dX9jt3VHDwk7+ixh+hxEZtfdORL4dffvmF7rjjDtr7fUGgOgP4K9j5drHojOR2\/ejqkc\/o0rtFnyBSe04hSm9Jl\/Z\/SITMQqYHv6vwCV1JnuxyAx5cOk7sIrCb+NvgmQI0+IdH1V5UFeganrNt2zYaOXKkSDPVo0cPys3NZQvPsWPWyVBV4r24DD7UlccoZyzy5Caf04pv9+JL\/vcKun58tzBqGLefs4abrrxc21CxF3bwWAWG9unTRwSHIimIl0f1kpKb0ryM6WYQXmnIdvHMk5uu7KLfrXTIpSc3nlS28apz7bW968rjlZBTOw2eIu1E0rhU5s6JL7vo92iDJ9zbMjf9afC4aeiv37mMnouOlzezR9FEMye+zFcrnEo8Rlq+cG7L3PSnweOmIQ0eoQGvJR4jDR6P0xeWZho8HtXKZRRcdNxWC49iBZpx8cVFh1s+VX14ae8KHl2T9JYauYyCiw4nT5y0YlU+L2BQbeMInkjUJMV24OTFqyI6AF\/k430iY9W4uPjiosNpB6qg8No+qjVJre74A0BcE8BFh3MiY5Gn0iCfV0CotHMET7hrklrd8XdLAKIiXCwaaizypMGjYlVFbV3PPFZkuWqS2rlAuQyMi05pMC4uXXHR4dS5P2i49\/IFHs6apFYuUK4J4KLDOZGxyFNpkM8dCuotlMETSk1Sq6hqK5a5DIyLTrwbl3TapCRcpSmPN\/ScETWcc8epc3VYeOuhBJ5Qa5Lqgr63XN5nzpyhatWqUVJSkrdZsmnFQWv9gXM0Yn1uYITRD9Wm0Q\/5q8MkDT6W5DOqDpHZbpXbVSYk4ffff7+JUiKyhAM6m\/Ma428cNUlBJ54K+qooWrbFdzNUEkdClfLly\/shEejDQWvKx+dp89GCAM0uDZNpSoc7fPPFwVM4dAWalStXptTUVN+ymTt6Ag9XTdJ4Kujrdwbw7ezs2bNUvXr1kLZHGJ+DlnnlWdAjg3o0D86UpCIrB09yPE5aoMm+8tzEnQOHR9ckvaUcrvMTFx1OnnDm2XnsPDWtVkafeRTeFHFdkzQWDTUWeeIEYqzKp4AJz01jsiZpXl6euK2K81FaWppnYYwNZSJ25Ly+7eGRvulImhw8gRbonF49ltIqlaHEGvVCyuTCxZPkK1Sdc9LhpuXLiFw62YInHDVJvWaMkVlQvLY3y1gp9xMq\/8HEwJ8v\/+ezdPk\/Boakv1B5koNz8sbFE3jjosVFh5MnqXtjfoOQjOGvzkquar8Dnjp1SmSM2bt3r18SSv1eSjtHHasUeZC2X0imOXn+D8FKg7s0jmXeOOWMRVp4GeMf1xMR8IBZAAj\/IvGY3+7Xuk2nXzMeisTQrmPEMm+uzJfwBiVy5YmGzqN5PddN3ljmzY13\/XuRBiK28milaw3EmwY0eOJtRrU8EdOABk\/EVK0HijcNaPDE24xqeSKmAQ2eiKlaDxRvGtDgibcZ1fJETAMaPBFTtR4o3jTACh6E9KxevVrU4UGYRpUqVejxxx+nQYMGBRUEtkoej8TxKFtfoUIFWx3v2rVLVOd+8803RRl74+OXpsqERlM+8Hny5ElRe+ajjz4SbD\/88MMicqNu3boqYti2DYd84Hnu3LmEsjT5+fnUoEEDUSQAZe2N95k48gOyKEGBCBt4cMELoQ8nTpwgFLK677776KuvvhIBnpjchQsXUq1atQRrAEFWVpaojNa7d286cuSIY9kSAOOTTz6hcePGidorK1asKAYeVZoKOhJNoy0fkuUPGzZMyI+AWegEOr3tttto2bJldPfdd6uKFNQ+HPKB58GDB1N6ejr16tWL7rzzTlFBfdOmTQJAL774ogAQV37AkBTgozMbeHATdfz48WJlaNu2qNDSoUOHaODAgfT000\/TCy+8IIoVAThYlWbOnEnJycmCbbuCWQjpWbx4Ma1fv54KCwupadOmxcCD6+EqNH3oSdy0jZZ8eGOPHj2afvvtNwEYGCEes24TEhL8iCb6cMuHlWTSpEn03Xff0aJFi8TlPzwAPXYnWI1WrVpFmZmZtH\/\/furfvz9NnDiRevbsSVIO7F6wa2ncuDFNmzYt5Ju3vpVj05EFPH\/++aeY1N27d9PSpUuDVgUZnY035KxZsygnJ0dsz5YsWULt27cPsFVQUECjRo0SfadOnUrlypWjixcvCuChzxNPPEEYB2H45pUHbzivNP0oMNryHThwQBgXdGx8MeGNjSv00O306dMdt7xOcodDvl9\/\/ZXGjBlDjRo1ouzs7KDhccFywIABNHz4cLEiceUH9DO3ofRhAY8TA3hbYsXB3XGA54MPPhCr08qVK6levXqBrngjYSXat29fABwAD\/b43bp1owceeED8N4BiBs\/atWs90wxFWVZ9IyUf3taQGyUuI\/mEIp8dn3gZ4GWHlQbgsXtU8wNGUi8YK+zgwXYM4MHBFgpDCfo9e\/bQ8uXLg5wIYMatPD1+twJPKDRDVXik5Dt8+LCoC7tx40bCywJv9nbt2omSl3i7h+vhlA88Yus9f\/582rBhg9i2NWzY0JZ11fyA4dKBHd2wgkfmP0CmGLlXtwMAGFy3bp3YH5tXJcm8E3isQOWFZigKj4R8ODvgrHXs2DFxXsDKg9qwFy5cEC8gbHfhMGjevHkoolj25ZYP\/GOLNmHCBBo6dKhwgNilgvKTH5BdAS4EwwaeH374QZxhsPXCgb9+\/fqClXgBT6Tkk+caHOiHDBkidCoNDvnRYIBwvuAALp0vHEbELR+As2PHDrED6dChQ1B1dTO\/fvMDcsitQiMs4DGmqsL5xrg0O4EnlG2b3crjRlNFWbJtJOXDNgeeJrh4sSI3adIkiGXoF55IbIHuuutWiZZQH2754JDANg1nnM6dO9OUKVNsi0GHkh8wVLlV+7OCB28XfG\/B2xGJO7C3rVMnOPtkKId7O+CFQlNFYdGSb968efT+++9bbmfdtrrRlg\/bTuw88GEbHkN827HLlBpqfkAVWTnasoLH+KFyxowZxRwCYDgUt7IdeEKhqaLEaMm3ZcsWYXTwuLVq1SrAMsAMD+TOnTtZPHHc8mHVxBkW4MF2De5puzMOR35AlbnkaMsGHik8PD\/4TmOX1lR+0MTBF94jt4+kRiHtwBMKTa9KjKZ8ONvgS31GRkaQzvABGX\/Hh0ZsifBtzO\/DLZ90DkyePFlEhiDqxM05cOPGDcvdil+Zwt2PBTxYmgGYNWvWiJglq1xr2L517dpVTPC2bduEi7V169bCc5Sbm+sYniOV4HRe8kvTi4JjQT6sPvC6tWjRIsjbhg+kRoeMF3nMbcIhHz5mIzoALzbMu5UzA\/MPefDBHCtomzZtqFmzZsVEqFSpkog84HSI+NGTuQ8LeBAX9eyzz4pYNrsHoMJHUux3rYI4+\/TpI+KdKlasaEvDCTx+aXpRYqzId\/DgQeHyx\/YqMTGROnXqRM8\/\/3zIgaHhkA884rue04OP4o899piIksDLwe6xCsnyMm\/hbsMCnnAzqelrDcSiBjR4YnFWNE8lQgMaPCVimjSTsaiB\/wdJ39BbbptibAAAAABJRU5ErkJggg==","height":140,"width":232}}
%---
%[output:47d3fbe7]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAM8AAAB9CAYAAAAFgzchAAAAAXNSR0IArs4c6QAAF0NJREFUeF7tXQtsVcXWXuWHPmiRZ4WKhbb3AgH0GkQsQW4hqKAgPggoVkUERbFKGwUBNWJR8IEvXgIKKuaPqOCDVEuuUQPhEWs1BEUISJXLQ2gQCgrysJY\/3\/jPYc7u7NfZj7PPPjMJudfumTkza+abNbNmrW9Szp49e5ZUUhJQErAtgRQFHtsyUwWUBJgEFHjURFASiFECCjwxCk4VUxJQ4FFzQEkgRgko8MQoODeLffXVV1RcXMyqfOaZZ+iWW24xrP7kyZM0bdo0qqiooEsuuYSWLVtGbdq0MSzz448\/0l133UW\/\/PILtWjRgt566y3q1atXVBkxz\/Dhw+nZZ5+ljIwMN7saqroUeAIwnH6DB10eNmwYA2pWVlZEAgo89iaDAo89eXmSOx7gadq0KdMsI0aMUOCJcVQVeGIUnJvF4gEetL979+702muvUceOHVl3lOaxN6oKPPbk5UnueIEHnSkpKaHS0lKCJlLgsTe8Cjz25OVJbr\/BU1BQQDA6HDhwIMp4oMBjb3gVeOzJy5PcfoMHlrT+\/fvTY489RvX19XTVVVfRiy++SAcPHoxY5JS1zXyoFXjMZeR5jniA5\/HHH6eZM2fSp59+yrZsMB5cfPHFCjw2RluBx4awvMoaD\/AALDU1NTR+\/Hg6dOgQMx488sgjTBvhLkhpHvPRVuAxl5HnOeIFnmbNmtHChQtp7ty5rI99+vShn376iQ4fPqzAY2HUFXgsCMnrLPECD7wH9u\/fTxMmTKDt27dHdVNpHvNRV+Axl5HnOeIJHnSusrKSysrKmPGAJwUe82FX4DGXkec5RPAY\/dgFF1xAb775Jl144YUR3zazxt177700depUwzsc0VdOgcdMoue+K\/BYl5VnOeMNHnRs69atEeMB\/ltpHvPhVuAxl5HnOYIAHkTjL168mObMmcP6q8BjPuwKPOYyUjmUBKQS8A08WNk2btxI8+fPp+rqamrbti2NGzeOxowZQ5mZmWp4lAQSTgK+gWfDhg3MCbGwsJBGjx5N27Zto+XLl1NRURGVl5dHxZUknBRVg5NSAr6Ap66ujgEH2kYMwFq3bh1NmjSJ\/W3o0KFJOQCq04krAV\/AgwMxtmc4kA4aNCgirePHj9NDDz3EQojhZ5Wampq4klQtTzoJ+AKed999lxYsWMDuKLp06RIRMs5B0Do4A1mJw0+60VEdDrQEfAHPc889R5s2baKlS5dSdnZ2lEDwDTfcIKTIz8+PWVi1MzfRmf8eo9TOLan9E\/1iridsBeveL6f6Q7upaXYetb55Rty7d3LHK9Twx35qOFFPp7Yc1W0Xz9ekeUfK6FbWqN38O\/+glw\/fzeqKVSi+gQdbN5l2ee+995gFTquVZB3iE0H27cjbWyN\/zhqQy0BkJ6Vf0oqaZDa1UyQh8v6+9q1IO9N7DqRm2Xlxbffpvasiv19fe4oajteTrF1ivqZt+9L\/NL8wqt3id\/5Blg\/fxLxpuSMps9ffd1lOU8KAB5Pg0MK7nPZXWj71H1nUvF87T+pWlQZPAhndSqXazG5L4w4eq9s2AEdcRe121Cg\/gAMAqZQcEnBL+\/gCHjcMBkrzJMfE9qOX2LYBQE6TL+Bxy1T9\/cv30K7\/\/C8V9u0b1e+9R07RniMnI38rzGtB9X\/+Sek22C7T\/9WKmrZPP7cfr6+n+vo\/KT09NsZMuPfrlcdh2cr5yqgOKwPvtDx+w2kdTsu70QZRVs3aFbqyZUOdvoCHX5Lm5OREeRPYvSTlDpTvvPMO9dUA6Ln\/\/Ewbdh2l\/v9sRaUDchgzDH4vPf0cIKxMOJ7n1KlTjupwWh7tcFqH0\/JhaYOdcbeT1xfwoEFr1qxhAVdgbRk1ahTt2rXLtnuOEXjETgd90uzbt4\/wzyxh1f7111+pVatWMS0CTsvzVT+R2oBYJ\/zzI\/kGHplj6G233cacQ8877zxLfQ0DeACaKVOmUFVVlaU+q0z2JADfSYRV+AEg38BjTwTy3GEAD+8DBpjT3LohG1UHsQUJZCaybb0X8lHg0ZGq062fXnmrC4AXgx32Ov2WrQKPAk9oMKXAYzCUVoXjVGt4aWWy2ofQzGgfO+K3bJXmCZnmkTHhoIt65wD4Fk6fPj1KCpxxRyaaI0eOMKKQLVu2GJ4t4DmyZMkSVgVn\/RE96sW6eZ1g+dFeQSCf2XdelwKP0jwxr9V8kmECYiJqJ5UWFJjgWoddXgfKGjnywju+U6dO0qcXUS8oe\/mzjEbOvyLYZQA3+y4KS4FHgSdm8GCSwhVKb9KL3\/hzIi+88EKj1V7vG5\/I0CR5eXlSb3jZBBbLiaAWnzSRaUez71pBBQI8uJP57rvvaN68eQTuAfz3RRddRA8++CANGDCAmjRpEmk3LuJwAfrqq6\/Sjh07WEzOxIkTGXVRWlpaJJ8bBCBWhZOsZx474R18YmJMzR4Qlm2LunbtyrZv4KMQy+sBWPt3\/vuYJ3jaEY8Ni0A2+y5bYazOj5hXJ01B6Znnww8\/ZIyUUP\/wBgAh+MqVK2n9+vWEpynuuOMOSklJYaBatWoVPfHEE0wAV155JX399dcsNue+++5jYMPzFUhuEIBYFU5YwAOXoz1HTlGnNuk0dYh5oKC4Upu9qi1uh6xytGm3ebJtn3bLxuebEbCNtCDKm30PzJlHj6wDGmb27NkEfzRsC6C2d+\/ezVafIUOGMC4CAAWAev\/991le\/ly5WwQgyQSeFdUHqWTFOfL1qUPyLAFIRqCo99y8ePgXF1XZ2YPnFTWNbDz0wIO8kydPlgY9moHD7HtgwIOG4o3KO++8s5E6h\/YAkQcXLlYTqFpoGmzreIJTJkAFRpwHHniAHUrdIABJJvAAOAAQT7f26UALb+1ua8ehtaQZWdG05wst4GSagwNKNBwkNXiMRgeviGErBvDAh0iPvOPEiRNs24ezESwuq1evdoUAJJnAo9U8AA4AFEuyY7FC\/Rx0HGx65m\/eFtEUndTbNr3BAU0U7gP27NnD7PctW7ZkAGloaGAAEVk\/ubBhrsRT5fjnBgFIMoEH4yCGWZideWTbKnEstQYCvUmOMnz88P8xtrjTKS4ult7raOu1ajCQtU1m+UuoM48MPDjvgKRj0aJFNGvWLBo5ciSLNQF4uIDxUJKYxMPk66+\/3ug+wcpBUtuWZAOPXS0jO8DzOrTnBiNZak3LRvXyvFhUcRbeuXNno7ONnqlar23afifMmUfbcAAH27Snn36amaC5BU27OvkJHpzJbrjhBurQQb6NAbDxsnP79u1J2y6rE9JpHXrlv\/nmm6hzo9X2WMmnd8Gpd3nKvQC0BgIRLPhdnF+1F69ie0Qg4qyERZWDCYSWZiZ0M3CYfedt4e144403pJ4KMGhx668VeZrlScnPzz\/LM2nNm6dPn2baBnc42P\/i8M\/vbozAI9u2aW+y+W9aJQBBftGSBIMGjBCyhHbjkVrcgot3TWbCEL87rUOvPLZAsEx66TYvusbwPun9nsw6JxoWzCY+6tczHFh1zzEDh9l3LXheeuklAoi1CUGFrVu3tjMNDPPqgue3336jp556ih32cbeDPa+IWiO2T68NBoiF6d27t67mAXhra2vZ91jDsJ3WoVcekxUBgF6Cx7XZkWAV8YXg7bffpssuu6xR613XPGeBAk2CcWDGjBn0+eefs+3addddxy5FtUmZqvVnl4rn8R95Vs\/EbrWskYcBNw7gghMr\/NVXXy0FDhrAL0kBLn4WMrok9ZIARBRIWDwM3BrkZKkn7uDZvHkzjR07ltq1a0eDBw+WvlyAw3pBQQHzJoAl7eWXX6Zhw4bRNddcQzgQy9xzFAHI31PY7wFOFuDEQ7aNNA88bx999FFDmYv7dZljKJz84BNn5hiajAQgCjzewdlv2apgOJ2xdLr1U2ce70CiV7MCj4HMrQrH6cRHE5zWocCjwOO\/BBR4AiXzRGqM1cXVrT6pbVvItm1h4zDQXuLqhVcEwmDgFiq9qMfqyuJ0y5Wo27awcRhw727RQAXviYqKCmlckNX54dbcVJonRJonTBwGes6kRt7jCjzqzBORgN33RK34ofHKg85hoDcN9LSr2raZ6FOrK0sYtm3ax7xaj3rS9EHesHIYiNPCCPRW54ev27b9+\/fThAkTqKioKIoPDI1Q7DnyoXBqqtY+I9li4FjKLnnTdNzDyGHAO62NHUK4g5gCBx641s+cOZNWrFjBwhJE3i3FnqM\/l52CR6t5ABwAyE4KC4cB77Ne\/BH\/HjjwVFZWMsAARAiKEsGj2HO8Aw9qxpnn5La1lNFjoOmWzQhUic5hgL6ZASdwZ56amhq2XUPY9WeffdYomlCFJHgLHjtaJqwcBnaAHxjNw2N6EExWUlLC\/omhuPEMhjMLJAuDwcAOcMRtjV7EbiJyGHDgfPvtt9J7Ha2MAgEefpZBaAFCsPHsoTaOnXdMsed4YzCIBTxh4zAwuhCVyScQ4AHnNLZrON+AuFBmW1cEIMbTOx4EINqDtdjCROMw0JIwaqUtowj2nQBEG4b9+++\/08MPP8weREVcD+K+gwYeCFIRgMSim8JdhoMnbgQgAMonn3zCiAr5g7N2weM1e44iAAk3CGLtne8EICL1FFhy4HQHeiSjBIqqm2++WdHtGgjJ6T1PrBMomcvF9cyDuxy8y4P\/FRO2cqBC7dGjByN\/79y5M+Xm5jIyO0X0HhyDQTIDJ3D3PHww9JzxFHtOcO55kh04CQcexZ6jwBMk0MZ126YnCCM3cMWeo7ZtQQFQIMGTaMJJVg8DjFPYwrC19z1GT0Aq8Bgg1apwkhU8YQvD1sbuJFxIQlC0jp0DYbKCJ0xh2Bhv2QNcRi8mWF1c3ZrTisNAR5JOAejGPc\/JHa9Qwx\/7qUnzjpTRrcx0zJMhDFuBx3QayDNYXVmcTnz8utM6nILn9N5VdGLzlIggMrqVmgIo7GHYsneAxJlidX7EOP0aFVOaJ6CaB8ABgHhKyx1Jmb3mmI57GMOwRSOI+HiwVhgKPMpgwCSg1TwADgBkJ4UtDFs892pfMbRzJrYjQ6O8upoHr7vhhS2873j48GHq2bMnew5wwIAB7Il4nhQBiFy8TrdtzOy84xX689cqatau0HTLZjTIdqIxUU\/QnpLX9k3vgeFAaB685zl58mTatm0bc\/2HTxueHlm7di17WnzEiBGsP4oARH\/KugEeOytkWMOwZTIILHgACGgcxEQsXLiQ+vfvz9rP43wQno2\/42FURQASHPCgJWF6Sl6PMdQoCDPumgcP+eLV627durFIUvERX6hzhGXjNbiuXbsqr2oD1eC35kFTwhaGLQODEYtO3MGzc+dOwstueAl70KBButNDEYAYb6riAR7eojA9Ja91zwn0KwlA7z333MMMBQcPHmQRpT\/88APTRPfffz9de+21TBspApDggsfOWSlMeeOueT766COaPXs2Awtodu+++25q27YtrVy5ktavX88MCQAXVtZp06Yx2cOIkJGRETUO4v4b2zw9SiQ7t+JcOKWlpYRHhTt06CAdez3yDTsTxWkd8SQAsdPPMOWNOwEIN1P27t2b5s+fH5mgMEnDiPDBBx8QnpnPy8uLG3gw4IoAJEzT3p2+xJ0AJDs7m2kb8BmMGzcuqld4Zn7MmDGMuxrPxutpHkUA8ndoQG1tLVt8QBzJEwYYcjUjbnRnOiVXLXElAMGtbUFBAeGJdxgMwFcgJtFFXBGAqDNP0KAZ9zPP3r17mbVt8ODBNGXKFEpJSYlaNaF5Fi9ezCxxigAkWPc8QZvMfrcn7uDB2QaH\/S+\/\/JIWLVrE7nOQwKgDQ0J1dTWzwIEUMUgEIKL7fm3ajbR8+XJ2LsLZLJaEvjmpQ6+83wMcS98TtYzfspX6tv3888\/Mj+3o0aONrG2zZs1iryZAI8WLAATWtsLCwsgYZ5+totaHXon895a6vjS+vIq0+exMClgaoXljrUOvPP+7V2eeRArD1t5Had9\/0o4X8iOJz9yIeURrrDg\/eB4s+PjnVtJ1DIV\/GzTMxx9\/TMeOHSNY33DPc8UVV5g6hmLbN2rUKEpLS4u0E0DbuHEjs+BBe8H8jbMVDs8gkreS9u3bxyZ0VVVVVPYZxWdp+OXn\/lTxNVH5O+e2m1bq9juPF+BJlDBsWTi1WayO1llVNl6ycAwxHxZC\/HMrJVQ8DzoNAOGfmLSapy67jA6lnNNMbgnLjXoA\/Llz53pibUuUMGy9aFBMftwj4nWOLl26MHFrNamRduLgAR0zp4oWx8w3zePGRPGzDrfc971us5f7cjsXzvF+DdtIc3CtLGqoBQsW0PPPP08IhjPbtnmh1WXtTTjN4\/Xk9rp+O+BZP6eWju09Qy1zU+nfU9qbNi1Rw7B5x4zAr+dlLTvzKPCYTpXEzGAVPN+9V0efTNob6eS\/J7e3BKBEDMNGJznwwcsm0ywKPIk5311ttVXwADgAEE\/\/uqU1XTcv11ZbEiUMWy+UQuysAo+toQ9nZqvg0WoeAAcAiiUFOQzbCnBEw4E688QyA0JSxip40F2cef676Th17pdlumVLxDBsvlWDP+WyZcuoTZs2uqOsNE9IAOCkG3bAY\/d3EikMmwMH94eykBZt3xV47M6GEOb3EjyJEoZtdiEqG3YFHgtg4JRXMDfClQWeCDfeeCN7nRvqnSeZxwK8FeC4mpmZqftLGzZsINwZgItBu02ItU7tjxn1oaamhoqLiz25JOXtCHoYttaQoZWfjJNNgccEPHAJgvsEVDom2KWXXkrff\/89c9Ds1KkTzZs3L3JzDBCUlJQwH7fRo0czmizkKyoqovLycsrKyor6NQDjiy++oOnTpzP\/Jtke226dsu6Y9QGuS5MmTfIUPBbWqFBm8VKrywQWqEtSPCaM5+uhGUCuyNPWrVtp\/PjxdPvttzNmHzisAjjQSlilOFDWrVvHJib+NnTo0Eh5uPPAQxyh5PAal5FI1NXV2apTb\/aZ9WHgwIGsHX5d5IUSJTqdSlrwNDQ0MM0C59ElS5ZEbamwDULUKphKcbjEa91iXBGXJTjl4A2O7RiiXVNTUyN0TChz0003EX4H4QJazQPBW61Tb0Ja6QNAumnTJgUeD1CdtOAxkiXnkgPRIsCzevVqpp1EB0KUl9Fh4XAKR0EQhlx++eXs\/8vISMCIarXOWMad9wEAU+CJRYLmZRR4JDLCdgzbNYQjQDvgQIwJuHTp0igjAoriW2VlJSMpyc\/Pb1SbnjnXSZ3mw0rE+4BQDbRNbdusSM1eHgUejbxAwjhx4kQGEmzrzj\/\/fENaWTPPYiPwuEGPJRtusQ9jx45lcVEKPPaAYSW3Ao8gJfAp4AyDrZcYEm50GRg08Gj7gL54baq2MtHCmEeB5\/9Hdfv27cxs\/ddff7GzSPfu3SPjbQQeJ9s2Pc1jVqfeRJT1wesBDlsYttZLPNB0u\/FekXDox30LNA7IO0C0mJsb7U3s5HCvBzwndWplZtQHL8ETtjBsfpkqbnExfrgO0BqLMAZeylaGi0Dd86CB4kUlyEZErwLeASdmZaO3XZyaqnn7jPrg5QCHKQwbF9m4ntB6URs5wHop28CDhx+s8ZgW7mlgmpYlfqGZk5MT5U2gd0kq1qEHHid1ivWb9cHLATY774ntTIQwbNnY62nXpNY84IUDYFasWEGIJpTxrWH7dv3117PLzzVr1lBZWRl7fAvm3127dhm65\/CBMDovxVonr9tKH86cOcMuga1Y2974YQkd+OMA5TTPoXE97zXdUYc5DJt33gj0Xi5MgdY88AkDRzZ82fQSQMXd12OlsjICT6x18vZa6UO\/fv0sXZJW7q6g2dVPRkQxrscESwAKaxg2BCGjrBLnStKCx3RZDUkGqwM8q\/pJWrO7ItLra\/OG02N9zoHJijjCFIaN\/hq9CpfU2zYrkyEMeayCR6t5Hu3zJA3NGx6TCMIQhm0GHAWemKZGYhWyCh70CmeezYe+pV7ZvU23bGENw7YDfDuydWPW\/B+pi3+1mhEjMwAAAABJRU5ErkJggg==","height":140,"width":232}}
%---
%[output:7a5e52d5]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg0\nTaille de la serie    : 205\nStatistique T_max     : 11.1194\np-valeur (bootstrap)  : 0.0680\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:3d3e90ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:6a568048]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 205\nStatistique T_max     : 17.8181\np-valeur (bootstrap)  : 0.0050\nPoint de rupture      : 2022-10-01 00:00:00  (indice 167)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 166\nStatistique T_max     : 38.1434\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2015-10-01 00:00:00  (indice 83)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 38\nStatistique T_max     : 27.1135\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2023-06-01 00:00:00  (indice 8)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 82\nStatistique T_max     : 33.4485\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2009-05-01 00:00:00  (indice 8)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 83\nStatistique T_max     : 17.5536\np-valeur (bootstrap)  : 0.0050\nPoint de rupture      : 2019-02-01 00:00:00  (indice 40)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 30\nStatistique T_max     : 7.7091\np-valeur (bootstrap)  : 0.1270\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 74\nStatistique T_max     : 10.2116\np-valeur (bootstrap)  : 0.0740\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 39\nStatistique T_max     : 6.5241\np-valeur (bootstrap)  : 0.2140\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 43\nStatistique T_max     : 20.1763\np-valeur (bootstrap)  : 0.0010\nPoint de rupture      : 2021-04-01 00:00:00  (indice 26)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 25\nStatistique T_max     : 6.0542\np-valeur (bootstrap)  : 0.2550\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expS_bg1\nTaille de la serie    : 17\nStatistique T_max     : 6.1249\np-valeur (bootstrap)  : 0.2070\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_bg\nTaille de la serie    : 203\nStatistique T_max     : 7.5150\np-valeur (bootstrap)  : 0.2510\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_fit\nTaille de la serie    : 164\nStatistique T_max     : 47.4248\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2022-12-01 00:00:00  (indice 148)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_fit\nTaille de la serie    : 147\nStatistique T_max     : 31.8607\np-valeur (bootstrap)  : 0.0000\nPoint de rupture      : 2015-11-01 00:00:00  (indice 63)\nDecision (alpha=0.05)  : RUPTURE SIGNIFICATIVE detectee\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_fit\nTaille de la serie    : 16\nStatistique T_max     : 3.5202\np-valeur (bootstrap)  : 0.5610\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_fit\nTaille de la serie    : 62\nStatistique T_max     : 8.4592\np-valeur (bootstrap)  : 0.1320\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : expA_fit\nTaille de la serie    : 84\nStatistique T_max     : 9.9363\np-valeur (bootstrap)  : 0.0810\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:2b042988]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAM8AAAB9CAYAAAAFgzchAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl4VEW2\/oFAQiCybwKmM8qioqLjjKymgyzqc0cRFU0i8JQR4SHgqCxJ1KCijOAyioDpKKuO+4IOYjojCoyO4jLKoCMdDRJ2kF0Cefmr+zTVl3vTtxegQ7q+Lx\/0vVV1q86t\/55Tp85Sq6KiogLxEqdAnAIhU6BWHDwh0+yoNdj28v3Ys+rvaDPhXez5agk2TrsWtZMboc2E95HY8YKj9tx4x+FRIA6e8OgW9VYHt65D2bTr0XrsItX3Lw9cjJYjZqv\/b3xmGE6e9B7qNG0b9efGOwyfAnHwhE+7qLYkeAQwB7aWYsdbjysOdKD0uzh4okrp6HUWFng27NyPK55ehe837rEcSYeWyXjzjq5olZIYvdGe4D3tWvFKgKhGEFF0azn2b2jYbeAJPvvqN72wwKNPc+zLq9GhZQPcnt7ef\/mD77Zg8TebMO26ztWPIvERxylgkwIRgYccaNTC1XhicOcADmN13eaYamQ1fc+j722srtdIIsXYpCMCD+dCzlOwfD2mXHma4j7yO7t7mzjnsfGyCY51E3qgfPNPlrUTO\/ZU+5\/aySfZ6DFe5VhRIGLwcKAU0wbN+kqNuUmDunh9RFecdXLDqM5h0aJFePLJJ1FQUIAOHToE9P39998jOzsbd955J66\/\/vqAe1u3bsXQoUPVtTlz5qBp06aQa926dcOf\/\/znqI4z3M7iHCZcyh2\/dlEBz7EYfqTg+fLLL\/HQQw8pcMUieI4FDePPiC4FahR4Tj75ZMW5mjVrprhRLHEevtZDe37F+vxLsX\/Nx\/63HBfZorvgo9lbxODRRTZ9YNFWVUfKeTi2TZs24fLLL8fw4cNjDjwCnDrN2qH1XQv9pCz7y2Ac3FIa3\/NEc9VHqa+IwCPnPflXnoa+pzeL0pDMu4kUPOQyLG+99RYee+wxPPLIIzHFeeLatqO6fI5K5xGDx0xVHepI3QCKAfBfKU4AmQAcvgtVgcft9uDOO6ehd++b0arVH5GaCjgcgNOJgP2NcBydC8WKwoBjIpfZ9+V7flu2\/WtWYn3+ACSdc3EANwqVvvH6R4cCEYGHQ3q2+Gd8v3F3WGppD4BsA2iM0yR4igCsNNG2eTxAYSGQm2tOHAJo0KC9WLPmRj+XIQjvvfde1eC2226LGW2bzIDGodteyvFPqMmgPDS5bvLRefvxXiOiQETgqcpMJ9ieh1wmwzd0AiSrEiTp2lTIiXRMTC4qwtJJk\/yqagInIwPgvyyJievRrds+ZGSkqd8u1+F7yckbMWTIUsyceYOfE1H7FovgiehtxhsfUwpEBJ5wR8r17l3iXoAc\/s4e2SM5k8t3uc0NNyBp5UqUl7fDzz\/\/Q11t02Y\/pk7dhGefHYxffvkloINBg+7GsmVD8fHHddX1oiKvKLdixQrceOONMQeeuLYt3BV1fNpFDB6Kbfe98cMRo6+K85DjkPOQ2xTYmLfU5z6oQOM4WVlAgY0OyKHcbu8+aO1aGw+0qCKGm7xtJU6Z1TGCIqH5KWib\/0mAi0Fc2xb+ezleLSMCD8W2zIJv8Oi1neD6uBSXdGmhtG5mxqIyQRHXKKqFso7JqcixnHmVwMsNHQgCIO6PcqpidRZvQncZYBUzHxurOqwv7gZWTm1xbdvxgkD4z40YPKJte+3zjX7FgZlh6O7fynFr4b\/x5iknYd+ANJzp2YFvHI1sj9wvvlGGyz4sgll1wL2Q\/JWUeP\/PfRBLOI7n5ChbC+9SHKNWUoo6zGx0+ZgAVwGrOnWbtrPlkxPXttleDjFRMSLwCCAyOjZFzw5NMNT1DeZkdcHH329DwfJfAvx5RCv3+XWdlcjWef53GN22YYArgyVFPB643W4Ul5TA40mFx+1E0VqHHxwERl5JpUyX6tNtF1dyKQsNHJ\/hCLbR8g1ENH38SWCIgxp\/EzzJXfsHaMKs6rC+XQ1aJNq20tJS8K+6lXbt2oF\/1a1EBB5OVucy5D6y\/3lp+NmmB6eyf7nO\/RP61KrlB89PW\/cp7duhpkl+Gjo8HhRkZ8PBDYuxOBzIQAHcWc4j1XTkMFz53CQZVXi8R\/mP97lfYp0qCm9TVR4JeHRVs1hRN838S1Qd3Aia8ePHY+XKldVtDeKCCy7Ao48+Wu0AFDF4Qn1TAp4uSzyYf2Zzv\/W1XJf+slwuBRwWj8MBt9OJEu74PUCO6zBbyfZpDBzkSqmpqh7rV1l86m2C0+lxo6q2sjeLRGzTvUBFMWDkWqHS0VhfNIhchG3bVp9YBwT7jBkzMH\/+fHUWV53KMQfPkEOHMK92bbQ6cAj169b208q3ntVvLuq1aV5ldl7lDj9X3+H7uAaBRYCZFYKHoCKQzIpwNKcJRzO2FfBEojD4dekcNQxyIGU1MO06tBn78hERcSJRVQt4TBch5VqZKzeA4WhMjtKqrnLcR+mZ0er2mIKHIt7lfy\/BpwM7qvFzacvydsEFD7hxAYoyXHC6PcjLdSA3JwsOTzo8jsNAKMrIgNnCNxKF4Mng4Q6lPicXkANF2RmK2wQrbocTGWuLlHQnWkFdDS1xBYxaMrM6RlCYqbkjVVVbLsK8PHMTjHDVjsEIF+L9OHhsEIzAuXH2V\/jXzzux+0\/nofy0xgo43E+weGgZ4DMXcPjYkIdSmgPIyyEX4alQIFfSH+vKylJiXbrbHQAsgsfNJzkAR6EHa3O9HI3AKhbOlAukZhFQjgBullGUC48zxxQ80TabiVRVbboIyZl9oi94KEZxlhxIuBBFXl73FXEq1A+bq7LCoF\/UqlWr0KdPH9BiIzk5GVu2bEFJSYnym\/rwww\/RtWtX1bvUMy6VOHhsgOdP8\/+NhZ9tVDWpFNgxsbuf+ygAUUwTWxsBlMOBwsqXq4ttwnUICKdmSkoxjQDKzctDjmbs5uc+iqMd5lgB4mDlIarT6YYn0+EXFwVg2UVFCjx2xDYbZKiyCjVtOz+c4z9AFeVCSp+hQe3bjliEpKVP9FXGf7qoRvDw4ItF09sTPO+\/\/z5Gjhypbu3duxf5+fnIzMw8wnuX9\/lMAcpTTz2FAQMG+Oux7RNPPKHcP9asWeOvFwcPgK9\/2YWrnlmFq85uroxC5fe23QdgZl1Av58Jb\/yg1NcNE+ug1f5D2J1Sz8tNaG2Ql6e4Bvcj\/CMQ8nJyAjb\/AXsh5MKNdBT5LOQEDATPhNmzkaCpbNN8ZgWyjyI3K8zS9lJq61R5MJQFBR7hfLyasXYtihwOpW0Lds4TDnjsxDAws0gIuggFIOSuFF2NRU6NNe5TFXjIUWjSxEJuJBbq5DiTJk3CAw88gHPOOQcjRozAtm3bVD0a4Pbv3x\/r16\/H119\/baoUqJGcR7cikPOe4b3aKfW0maU1wTNrWSmezzxTEZbgOXv3ASS2bhDgiiBaNoKH3IQgoHj3Qbt2SqQSDVxepVVcLnIUeMiBhPMMy89XnGfXwIHovMgbfVM0ctI2f1gCJs6aeNiqjuApARKG8owkDRPnlFf24V1teQUFyMnKsqWqDgc80WpzxCKUvY7V3sbkvpnYRtf1fv36IY8cPSdHxYCgZTpLamrqEZxHxLYrrrgiznnMXq7RgoBcZ\/zf\/oPC7C4qBJWZhYEZeJxb9uJtRyN17FLpWYC5vkO+vd29It3JgwcrQ9Bto0dj++jRASLZI0kX49l9z2Ctz8RU7W2cTiQ9sgIJD67DPfvew4Ryb7haciUWEedmNG5cKcnkqj5V8QAJ80qVA1H7Cy\/E6G3bMHr7dnXLlZuLrJyc6gce2e9YGQAK5xFrWQBGziPv3movlJ6eHgdPqF8\/IzjIaYrWbFVcpUG9BEvwGMW2a\/aV44UWyf7Hywl5N\/rbuN0oHTYMpUOHYk7lCfRsA+dR90pL0e2995RoJ6IZKpVLNNXOReXex+fUYOQ8bEuRkH2qUrmHJnfzOIEFK1ag3Zw5aDfbBzyN8xwNsS1U2lvVNxV\/atXyVtcA4v0i+BQJBktZK\/BQMaBzHhmD2Z4nznlsvFER2zqf3AADZ36F5skJ+Mf4PyjOw3ssesRQ3YiU9y5oWA93JNbGo4kJRz5N39Dm5qqFLseiFbIg1IFppcEaLQ2ojaOlATXSbiCzUivtLM6Dw3eYKlzJ39bQp6PS94cfaEpvSi2tqXepMChwOo+JwsAG2S2rmIJHV1OTXpykbuRnEOmswMOHSv8yAIpzFNu4D+L\/PR6PcnEfM2YMysvLQbHtnnvuwU8\/\/aT2Qbfffrs\/epE+iRq356GYdsVfv8COPeWKDn9MTcF7o89XHIcxrHud2sjUs1QPFlLvkXSMqFvb2pdHA5BuYaBr0pTqlXuj1By4fBpXv9kaVbQuVwBX8mvbqIyo3CgLILNclSpwHTxaWyoMqjrniWTBG9vq50GMmtNyxCyUTRuksiUESzFiuQh1dbU8UIAUA4elNQI8ujat3+lNlYg2+Y0f\/G4IoS4inrZwvVfpHWB0FzU+hJok31kNpTUWmrI5fV9Xv7UAD1gZ08DtVupqFnduroqboMx9ePzhi6OQI19mh0MBzOV0mrpORDtIoX5I2uSyMdj6ch5ajVmIHe9M9+fsqSpiqC0LA1oXMMADaRbMhCnUFxpm\/RMePFaxpymeiQ+PHdrp8QM2rFiBc5s0Qb+6dfFzaSleeeUV1UXfiy5Cly5d\/N0dOHAA7uJi3Lx0KRIXJqE4JxWZhYVBLQzEvs2\/F\/L1aMc6QdqevmABytu1w4Pz5h0R60DA0zz7L9hccJfKqxNJ\/hwdjOWbS\/3gqdi305+3p6r+q+sirK7j5nKybZ5jjFfA2NQM\/GEXPJSnp0yZgmnTpil156k\/\/YTNDRvipEaNsH3bNjRq5PXt2bFjh7pfu3ZtHDp0CJs3b8Zvv\/2GFi1aYP7yRO+pDj+cHo\/iIvzXWOzYtlm13d2ihVJtf9CrlxoTx3DewIEYPHjwEaF86X+zZ7lXbWtW7JzP6O0kRlvT63Kw\/Z3H0TxzGtbnX4x6HboHjZ5TXRdhdR13SOAxLg6j+7UEerfDgViHRFuwYAGuueYaTJ8+XcWRrl+\/vtpk3nDDDerAjSfU1157LcaNG6e+\/LS6\/etfv8Udj57hdycgeLIYQocSWGqq3wK7qnGctHUrXCUluLpJE2\/4HZbUVCxYvhy1+\/Q5AiSM8eZwOAKuC6doMXQGNs0ZHTHnkfHqtnG8ZtcMyHQR6gahZgQxEd8415kzZ6raEmGVscF5nUVCdYkG7u6778bcuXPVoSk\/eiz8UBYWFqJTp07qbEgvep+yDqh0qLFW1QSS0fktGIhkQVJjQxA9\/PDDqgnB06NHD\/9ClbjSAh4BXUrKY8ifV\/fIsDtUwFU6w2WlA8tOW4ZfOnZEUlISNmzciIy0NAwpL8e7d98d8Aw+l88ZO3Ys7rvvvgBTFKvr0RbbgtEr2P2gtm1mHRhs29hHcXGxHyACggkTJqiP29KlS\/0fMTvgYTt+EFl31qxZGDVqlPqtlxrJeYQAYl2w5LutfpqIQoFnPmaFex+qNgkIAUOo4JH6Q4bMRu3aGWqfpAc7DOUZtMMSjqf7lFhdlzlFW2wLBpCq7lcJHp7z6AoCUcSYGIaSY8ii159H27Xzzz9ffeiEm\/DspyrOEwdPFW9MgNOucaJfNS3X2EwOTfUujCIQXzqvGcU2WcRmnKeq+nxWKM8IleMYyREtrVtVNm529k5Vgoe2fUbw0GjUAB4RuZiuhZbVuoglhp88BCV3opgm4Jk6dao649ELY4LzAxfnPBYACjUzHBc1TTr0r7u+ePkYXakg4hQzGojYFqx+KM9g\/2Yn51Yn6joZdLFtw7TrTJNT2Vn0wbgNuVvDHtcFddmOFnj08eh04J5ErKZJ47POOktZYMc5T7A3WMV93WZNxDSz0FPGE2p2KV8nWuaKxa5x42jkPGyn96XXD\/UZ+uZYpiin5bJp1q8bE2dFQDbbTe1ytmiARww+ZZ66W4EOHj1hGDmSlcIgLrb5gnwEy37N1XAIFagNnz2Vb3kEC7trexXV0Iq6K0TI5zxiXWDc85CWFmKb8YMiHyejv46I23HwhLAw49mvA4klmQwO7dnhvxGOuFbVnsdOKvloaNtCWAZRq1pjtG2h7nGiRuEY7UjOZIxnMRJ7zc6ij9bUonXOE63x2O2nRoBH1LbUqhys3xTDJzyO4Zf\/EVc9\/QX+s3Ev9OzX1LjR7u3ui9MCUszbJahVvVgitNiiGaOGytj1OG9HI4u1vA8+j1ot2TeaHTbqexc5yIz0XRjbG8djPM+pDu80VJrYNs8Rhyhu8n9\/xTB\/9uuU+t6znGu7tlDqaop1r3+1+ahkxI4l8ATbyAe7Ly8qXDdsK\/CMHj1aBRHUC89mvv32W6UZS0lJCXWN2KpPEypRstBNu149r4t9sLJu3ToVrLHaWhjoBps0i5GU6zJx0a7w68ZCANGshucBJBSJNuDaW7CoPB3XJxTj\/b+94K+nH37yIl8eCcUzBJ7mUzXNMwWrKC3GZFR8lhDaOC79BegbX+OZA1XfMherNjodjM+hRu7afr2UwWa9W2fif8dN9PcnGbdXf7YM2VlZuCHzVhR98qm6X1UkGgED\/V+M9Dd7D1T5s5x00kmK8\/D85ZZbblGHzywVFRWoVauWss2jnaBck+vGRc3rej22lTZSl9dYj8Xq\/3IvGGj0+9U2YuiaNWsqCII777xT+apbZYnWOY9YBlC9fM1N2XgNGSoQSMsf38BL7xZh2x9H48XsLnhswih1pkOzDJ7g\/+tf\/1LJqVj4TBp7cqEsWbJEBYswfn2MzxRAsB4XN\/tk4eJ588038eSTT6r+JVgF6zHztcxPHLSs2oiZED1U2Yag41z5XHI9GSufM3v6wzi4YDwmrGmK7r0u9FtLiJ1W44RyBZ4\/9M7AI395ImB8tBUzFoLnrhsvQ8n33yH\/tA1oVPeQv4ooIB6bNVcdRuo0\/P3vfx8gtk2cOBHt27cHDy7POOMMZSeoc57ExET1sWNEG86NMaLJlR588EEY2w4cODCgLgdEWvDDmZGRofr95JNP\/OY8vNexY0f1kbDLedhntY1VLeDh11++mmZfDVnIZ\/XoqzhMnU2rkfLxo9jZczwOtuiMU1PK0eqLp9C7Zw+MHDNW7Xk67\/kXCmc9g2eeeQazfW7NXOhcKDpIrcQxYx5SvZ4OCp5LyPiYrJfF7NxI6vBDYWwjB7c6xxJOIdd0bkRx6+3be2DMZ3XwWIdf0PWk\/dh3qBamebyJja9ruR25Je1w1fU34Z6cB\/1nU1bnWMIJhfZtmjWGa\/4ipLw1CYntz0DFRSMDaGYU22hEKwDnvoZj5jslvfV7YnwrHxD+1j8OetvJkyfj\/vvvV96g\/HCwyHszZhQPd88TCoeKtbpqz2PU7ZuBiAvvlsws\/NrqfMyYMhkNd3yvFqgsBuPC5ERl8R8N8OgA0YnKsQuH0U1GOM6q2vAeuZ8eVomA4tfZGABDQESOKdyO3ERfQMOGDVPuxwLUYPu1Xb\/8F2NuvBybUtLwl4HnIHnrj8oNQfZOIh7KmPRnySI3msjIOGmUKcAygofzJgc3a8szHHIwAZr+0YuDx+DPIy\/YaDZO4gl4ajm64blpuRj\/9Jv4z\/xJCjwNTumCx974Ctvfn4oe3bsj60\/\/B+bt6Vvvm4g5jy7OGTnPjffei8ysLIy87DLLj5LsmbiQKJZwQy0LWhoZv5pGzmisx4VGkHHfoX9AIgHP7i0b8X8DnSiraIhZj+Vj3+yhKp71ni\/fV4EQ6497V+2tTjnllAALdFnYRs6jE0TnLHY4jx26CHisxmNX2xZr3CSU8dRavnx5hb4AjCxcOpOvb7uu6Wja40as\/XZVAHie+KAEqaVv4l+f\/hMPz3gG9772PX5b+pgS42TPY\/YF0\/dPdvY8T7\/zDvrPn4\/9bdqoPdRJW7agb9++eG3bNux8+20s7d0bn3\/+uZ8jUJ6WvZF8oWUc+j7p1VdfDdjXCBcyjl32Q2Z7xKrEymCch2PSNXQMDi85feS8iO8m2J6HNOS+QxeLg4FHH5veNth7q2o8NQI8FNt0jRZfopnaUBfb1qdddcSeJyWxDtzjf4+XZj3lV1kOyR6G+yfdFyDOhLLn4Vh0e7VBlarWaZddhhm\/\/oo\/nXFGQGZr1p32yit489xzGX4NKy2cuqw0dLpoRs5LZYZ8VYUT6dpG3WLYTHtnFGPtgCfYV08\/azNq27hY9feoaxh1GlLZwA+FvucREZsfDBZpKyKd1Ucv2HiCzae637d9zqNPdPIDUzC3wBvXTC+3DhuOiffdaxq3zQ6hxKHaLLsOA3TQ59NG\/l7lqs1wUozwIV7adJrMzIyZuBd2yBGvE+MUCAk8\/NIMvHUUvizdib3nZQN1Dh+EpdQpR99dr6NhvQTcNTEPd79egicGdw5qYSDRQo3AUVlBfNF1eI9prvSIyyqVYt7hvKNaGDcvYIqAtRoKeY0e17SVpLbcInVPWK\/L6DotnYRi4yaHpXWatkebCe\/iaFglhDW5eCNLCoQEnv\/8VIbr\/nc88h\/MxahX12Hb3nLU5mFaxSFUVNRC7QO7kPTVfOw7+0ac1r5VQE5S4wgk1BNjp+l5evR6xICAiNcFCwzppsf9EFd8AsPjTfEDB88QHd7U8zpQ2I5h2ehaHw0ARTNNollyq+Tu1wcN\/hFf38eHAiGBh4ahV992H8pWvIYBWeORfHqGMsnh9cH3PImSxX+t8hRdD3rYY9R5aOFoZEsMkzTya33chuJXYZ4HOY7Cyow6HnjcHm9y31wPSnLWwsX\/a1xHZWEw0DeXOZ98sSlUrLcw6W\/XDCfM7lWM7G+fvQvlNz2N2inNw+0mpttV20NSKgzsUlZcrM9q8CveevIebNqw3t+0doNmmPuiC926erMgGIsebveXOrVwdf26WFurwi\/WMaEvc8ac8muZipXG+MoJCQlAhgd5zkKkMn+o26O8iV1uB7LgUmlIstLTMXHOD1jX92M44FRpSpi0ypWZWgkrJnTKY3YsJjKtzAjnRGZhtkohotKI8F8we1yOCtPL0vNTD5y5QEKLIHlNtQnSinr\/z99GzCHMOM\/G2k0wfX8P\/POrb+2+pmpXr9qa54QCHr6VULNfy5vU8\/O8Xj8Bkzbtxfzt+\/0ZszOY2zqN+Xm4s2EaNy7ePHgq87q5nTnILnCCydsyM\/K8QdnV6k9XAMnIApzpLuTkERhsX6g4T7HTgdzUHByY1B75CX2R43H5MsI5kMlQvJ5cryjo4HPZa2ZlRgUgx51tntPGZFlGGntAupR++Ltt\/if+AIqiKatuiXrtIjie0NcGpXR37dsSamPVFxsxctcBfyp5T4aHa1jb2fiS5iATaQUOeAgQNwMdUogTZpmNPEcmcoucqEg7rHJQAMitheyCtXCUOJCem4diym2VdVlyK+Pr5uZ5MwOrR6bDx5HSkZeThZzCbDh9Nni8rSsE7MZRs0ESW1WioeK29aDjVKk6zy+kPY\/Q13hWYqS7mYWCgCcn60xMqXMkeFSaHX\/wT12\/5kFaEdN\/MOVhni9C6GGFtQcZKqZ0QfZhMYvgyXRloDArE6keJ7JcXsBRY8eWGYzE5AEKsj3qmiPX4UtmVQuurCKUpHtUQisWO+kUxfnNSIdQtG3GtuyT5T\/t+wdYMRynNX7UHlvjwENK6nHRhLJyjQd4YqEsJ80itj1957l4LSkBb6zaiOnJdf1iG8g5\/Lt8hm1n9E+nWuQZuYA7h+DhUvdu7x1eVUGlGOZGiYOcxrvYy9uVg5nf0t3MqM1Q7g64HB5keXJ84HEo0cyVyRSKhB73OAKebLidQHFmph88wdIp6uDa9vbj\/kg3jHpDg06mjw+nxBp4rCKJ2p2bflCr207WOPBYhWaS6zTroK2VpOEjgUVhcP+Vp+GaPeWotWUvPjunxWGFwTndK+PlFqBlw1PVqU5CnWUorVOi8DRzzd+xcF0\/FGSnoTjdg\/dubIs6dRJQt04aPK4iFBRnY8jS573gOViOcV\/tQuM3XsewB+hX5EAh0jHhYCZ6H\/xBnRWNyK+rrufkeZBNRURuMSbMKQC2ePXcHXZXBKQV2fHW4+rshWV9\/qVI7trfDwqjSY0oDiLVwsUSeLjArSKJ2jHD0dcL6+tJgmsceMQsQ\/xzaFEs5ij0L7nyyivxxhtv+IPeydfp5c\/LcNvc77DzT+ciAxX4+2lN\/B+uZT9sVzl\/Zlz7O\/R7eAy2\/nUObm11JUoqfkTZgLlAkVPlI730vO74rpdTmc9sSk6C88n6KHJloOVHr4K+KiUVh\/BowiFMWbgIKQvnYtmpv0OxOx3jP7oSN12ViOmzSnHT2R+j34s\/ITPPjfd\/uwl1HqiD3Nw8bL2qN35dugzz\/r3Wn\/pEd6c2A49oyAio5HMGYP206wIMOvXNf7CvtJm27dvG52PU0i0ReVrqZjQiUnMsEiOP\/2eo4bvuugvvvvuucmWgDZ1YZdM\/yiqSqHFOukhPm8nWrVureG96lm1KKAyzTAvxGgEes+zXV+YsAJY+5KefGCWaxXyWuAa39myLAT\/uQPPzWmFVMjmAtxA8NKtpldYIDw67Fa1vuhlJiUkqy9ikiW2w4uMOKicpVdZld9yJrtu3q\/v9Sk7Bkjk3Y9+ESUg49TSVsGpl3dr47r4JXlW3pwR5a4dg4F97YOGCurjsii8w7uodWFenFuAqRkJCKW668HlklfyM378+EwkvvoCpGysCwBMsnWIwg85goOF9PT8PXRHk2rujL4oYPHo8Np0LEBS0c9u5c6dKFy9GtHSg48LW2xlt\/2gjx48mwUbHOOFAetBJEePpPGcED+dHn6oaAZ5Qs19bLRhJVfLPey7w68xYl+c86U2TlM6AIFk9cyZWJiWpbnY1G4hNnz6KNiffgF\/WL0ST6dOx74ILkLBuHc5v1gzO4mI81qWLypyd8sor2DlwIL4eNw43tmmDgTt3olnz83HXjq4KKElJK7H2x0zAk46Efglof\/AWlTC4yYwZeHTTJlwziLHeAAATGUlEQVR6yfmY+lKRHzx2FAZ2wBGsjpWY98mStzHktlGWnIcfHK\/eMNB8SX+e0V9LV+joADEGCrGTo5Q+TZs2bUKTJk1wySWXqIDukjFB2tdozmMn+\/Wd8\/+NjISvsbx4qXLBvfnmm3HhhRcqX3e9CHg6DT0LdVokB1gYyEJg\/Y86dEDpiy8yw7tSJz\/brRuWr1iB1feuRt+EeSh94AHvYWp2NibNmoXSfgl48UWmggfK58xRi+njvkMxeHU3VOTWQoazAg89tELd795tNRzMw+jIxoKk21TC4KFMDnzDDXAVFaHEUalA0Aatq6qPZjgps4g77469FCNfW20JHlFScr6S\/tEIVGMkULkve5kGDRoo0cqK80h9s0ii7733ngLNjz\/+qOIlfPjhh\/6QysJ5JK4198A1bs9jJ\/v1sCeX4H\/Pb4DL+2coIj7lWoTp3zTAdqRAz5og4Bl0w+nYlXqSErPMbNtoKeApLFRqaN7n4qBJTUGmByVphch0uNQ1D\/OLpqaqxV6YV4Icn+WBIz0ThYUOeCq1apmFrsqj1GJlEcqvNEFawAzYPpsc6vZyaGXKY1lfPpkq0z1qvjfMDBeNWNVWh61f7kzC2DVtgnIejt0KPPqeh\/W4l6GXKz18xb3abM+juzVYRRLl3qhPnz5gFBz6OnXt2lX5EkkoLO55Ro4cGeBaUuO0bSK2XX1eS5W0N7v7yf4DTt7bs20DLm+3G\/3798fBgwfx+uuv49JLLz0i1JGAJ\/\/K05SaWqyqyV0IBgUI36eOVgBOmkLzwNLhYHZ5ZUmdw3TX2dlwZnnPdtw+K1Gn04k8TybctE7wBpdR2jomwCpgP5Xqa7czC0hPPWzLVlICBy1Kea6TmYk8XzaBYOAJJoZF6\/6x3BOEE9\/t2WefVdxk\/\/79uOqqq1RwEgnmb8XxdNocy\/lF651IP7YPSfW0ihLgUK6l7fkGFzT5NWBsTEk4aNAgNGzY0FRsE\/DoN039ecSXQFl+emu4PE4U01rA4VBWPOk5ynXHm93aZ67tSPX2rFwVCCAvlLwg8vWlgOp0wp2ZqUQ1FmKOfQUDT7jx1kJ9gcdycYUDHuN8dMWCWRgzY\/1jOb9QaR+svm3wBOvIzn2r+HBWB2js0yzBFOsPGnQ3du4ciJNPvgEtW7b0goQNqAtwALv\/vRsHfpit9kWHevfGyJQUjKhfH2UrVvgj64j4YBRrTn3hBSzp1StgShw7wyyJB6l+0+xAlPuXXZ+8HLGxaHVeXHbWRHWeX8jgMcZSs\/ulMSb0ldBIjGLJwBwMgMii5+fRzwz0QIeiCjfWl5fFZ9333HMqw1nXxo1x24IF+LJxY\/Ts2VNtaHv37q2qfvTRR+g\/YAB+\/O9\/sX37dvTq1Qt79u7FJ1OmYGSXLv7UjvoZlhE8VlqySA9JZS7VeXHFwWOggJ51zcgVzEx2rAjIRRFOQl+2C5YZzkw0qOpZerIttjXOkSfip556KlatWmXJeZgVW4xGJd1iNBzZ4uCxA8HjUyckzmM0yzFyEzsZ1WSaskDDTehrlcPUjIx2n8W2xjSLsunVx2lmkqKrs2snN0KbCe8jsWNgzOhwXnEcPOFQ7di0iQg8xn2AXfCEkmzXKhu2XfCE8iwjJ9WzQVP1ajR2PRavKFbAE6lhqOxfdbs2XouV+YXzLkMCDx\/w4JSHsHRdXXyFjmiwfAb2OS5ExiVXqeS9Tz0+TY2BsdisSijJdoUT6DlJQxHbQnmWWWJfY0gujkc\/+wiH4KG2iYXFFalhqP4eaY0gpj01Cjy0T7v5meUoeTUf235eo+IVMC41r33\/0oO4\/n+cKvRUVcA52gl9dbHQ7rPYxiyxrz4P2aOZadt0g87Ejj3RcsQslE0bhJYjZkcsukUDPMfbMJSHpFTUdOrUSR3M0o5Ogt1HY36hfpCiVT8kzhNJZjhdHS2Dj3ZCX+k31GdZJfbVE\/hagUc36Gxy2RhsfTkPrcYsxI53pmPPqr9HHEYq2OKi+zr9mujfVBQQnOvwEjnehqEiiRDENVpss5v9OlrojvV+dJV0+eZSP3gq9u1U+Xtaj13kj0cQzlyCgScNaX7wrLUw0DnehqHkPDV+z6NbGeDgb6j\/eQES1\/3Tvyb2t\/0j2g74E94afUHQYIfhLKRYbUPV9MEtpWh6XQ62v\/M4mmdOw\/r8i1GvQ\/ejfkgqnIe0sQLP8TYMjXMebeWKDE3Tdp0werIpOx6GsQqGcMZljBoarUAhwTiPnbHGgmFojec88qKCuWHr7td2Xm68jjUFogEeu\/QNx7Ytbhhqk7rGAO9nn322Ovs4VKeuygTXYPWrOLPjqX6zFpvd2qpmljzLqiEXXElJiRqHVcoUWw+NgUqxDh4jieyaaxmVO9U2oW+wNaInmx2X\/wRuLvgGKcum4tcyr5Wz0eGN1+xY1AZ7rn7fLnjs1gvl2cez7rEEz\/GYZ3Weny1VtdkEJe70S8PPRk9HckByXbP0e2eeeabylZfAEvQuFC2MWDcbn6PnJOWLlcS8ej5RBqtgYR+SkFjy6PAcioX90ulLTy\/I3Jt6\/+JFKekFrfKzGlMsSj7SVb8mYtz3J\/vzk767qQHmbmwB1wsv4sv\/rlPJtnjY+8ADD6gxTZo0SdnK6YE2OCZjqc6Lyw4Yq\/P8bIFHF33W7auHq55ZpbJfM8g7HeFe\/3wdem15FY3rJyjDSTPwkJBcwMzCzKgqcsbDEFWS7UzPYi0BKCTnpw4eyTkqgSrMMpRJ6kR97Hz2uHHj1Am3njHOmPRWxmgUJfS0ifl\/HoWHbr4IX2z4TWWvXv+7ARj5+mF3aR2YzFTHxFHGfKl6tm7JvB0Hjx3IxUadkMDDdInPrdyOMb2bYl5hgUqXSK3a5h07cXXWKHRqlYS09u0wePBgtUgJAGPiV6NYpS8yu+CRw0uziC5GDqWDh\/f49ef4zjvvPD8nE44lyXLFno719INSXWulm+nwoPTdUX0w6sOtfs7z\/t72mLe1LQoKXwhI88iTdbMxybPj4IkNYNgZRUjgIecQ0UL36yGIyA3ox04wtG3b1jLleDTAI4td7KT0HKRVgUdPk05fe86HXEjaiAgohJNU8johjdYLRisJ4VbPDe2D55evDRDbxKaruoKnKhMlO4vNrM4JL7Yt\/vAj3DEsE0\/PLsQlfXqDMdyYwHf64E74+vPPlBjW8JSz0P2MUzB92lRTsc2Yll3EqnA4D90D9CTEVe2NjNo2vixyRbOco1Zff7OXrjvq5U+8G8nvPYQxn9Xxc56\/1z4Hc0uTUeAqPGE4D1PLl5WV4bLLLlNSRTTKCQ+ekk3bcU3WKLSpswMvuArAfc\/Ep17Cv1+agrPOOgvPvrAQfS+5Avs3ejB\/\/jwVqJCLm1\/uUMQ2s32NRCXV75mJXMKFgnEefdHrSgHj3kgHpywSacuFM274ELx9ew8\/YE4dPg13TF+g0tTLnkzGrnPG4yG2RWoYSo5tFYPNCCCziKGiGKqRnEcMQjuVvROQyPeWrGzkTp6k4lAz9NSmt6Zg44YyFQOMX6hQ9zyysGbOnAlaLzAajtvtDhCtdEUA67HQrZou1RSZqCKnpQOVEBSpGN\/gs88+UyKaiJy6gkGse42n8Ha1bde32o7h7bapccwqbYJFGxqjZd0D+GOjPfjnrpOOvtiWkeENaMIgDkV61tbDSzUahqF6iFxyIL5juxFDq3JROeE5j1nQQxHbGtRLCDv7dTTY\/oneR9DFlZZ2GDxrzSO3RWoYSi4ssdiE3vJxIaiCRQyNc56Fq\/3ZrbnneXzxGtRf5cLid95Gv4v\/Bxt+dyWSv3wReZMn+n01TvSFfSzmFxQ8wnk4GAvwRNswVDfPshMxtMZzHgY6\/H7jHu968VlUVyQ3w\/723ZH441K065ON4e1\/xjef\/9M0SMaxWGgn4jOCgsfGpCMxDJ08eTIef\/zxAAc2PlIAmZycbCtiqNUwozE\/GyQ4KlVsqarNNoXieamnn+BLCuaRWdUsqGjgX6wUZllQmRbCKAw9y79IC9X\/48ePx+jRo8HEt0ez\/Pbbb3jllVeUNi0lJcXWo958802V2uXAgQMqdNfixYvVvvOMM85AkW8PxiwKVsXO\/GI1W3ZY4CEhJO7arbfeipdfflnFPqZsG4q6VycoF9q8efNUbDW+OH7RZOESUHv27PFfC\/ab\/QarI\/e3NWqE1IqKI+qzj8aNG6vo\/6EWzoULnslq4yVyCsRqtuywwUOSGA8MrTRUQj69vrGu3HvuueeU6QwDhCf5UoyQo23YsMF\/LdhvPk+vU1aWhFatAvv4obwcT+\/ahbebN0fP8nI8Y3gG+wiX88hcTtQM1pHDwX4PsZwtOyLw2CdBYDw0ttMjg+pAvL8y1G3Ptm3Rpk0bP3j27duH9evXo83+\/eravtatvb99dXh\/06efosUf\/hDQ5tNNm\/BOixZ4pCwJ7XqX4x\/\/+DmgTYeEBJTOrUyAleXB4tWt0SnxcJ+hzM1YtzrL8ZHM+2i0jWVahg0esyAbJJ6VKwLrVxXpk\/dHLFqE7aNHo115OUNOo++yZSpAO0WszZs3o\/W93sg827dtU+dIrVq1UrZ1DIDRuf7p3nfH846sLBW3mqlHWDzOYlWnVeu9WLx3sWpDwLGPfv1OQ+lH\/bwBrosZRiNdBdPg852Z3kDyoZZYfuGhzuV4149lWoYFHmMgQjsE1u2iWJ8HmT169PAbXvI+MjKUaMaUIMaycsQlWN8xEXCkonPn29XZAk1suFlleosRZafDyfQGjPXu8qaK1wuxwb8uZ96BnWPHqjbSx9vNdmJGkzdQmjDP186DtcwaReRYqH+rmnMsv3A77yrcOlwXdPsQg2HpRz9gDbXvWKZl2OAJVatmBzx3TrsCCzr1R\/68DorGrrUqacgRZeu2rQF7oPqnl6F1630o4x6pyA2HuwR\/7tYKu3buwu0bLkZSWRmY3Kf8gw+QRK6zeHHAnoiioMftVvWaP\/00NjfcrDJfMUGWPwNWCG89ll94CNOIWtU4eAykDCWoO5vaEdt4kv3888\/js8\/ORJcuKUjqXIZiN+AqBMp7\/oDSYacxEyK+Wxy450mrzE3lKQGYsod5dZjzyr9P0vZFxn2S\/hu1aqk8p3C5vPumcOQ1H41OBPBwDgxQKGZPdGunTxKLKHt0ywWxYxTOI5nhxMyK6m\/2IVmwxVaOlgtiTqXnSpXlFsu0DIvzcGJmoWh53WrPo4ezZT0rhcEdd9yBXbt2oXPnzorQLNzzLF261H\/N+Ls3c5Nyj1Ra6k9KFazNEfd\/+AFL\/\/vfgOeGe75g9sK13FxBv+iS7lGvyKx4dks47fmt0L8XnENxcbGKjqTH7OZ+kaA599xz8cUXXwSIaCK2DRkyBFOnToXkICU46LhoBh6K3iz0mzILLHPCgccsrrOdF6srGYxemrF4NsKDSf6FWsxeONNB5ppLoUd0X1Fx5BMNeZGrHFI47Tk2XzpW1TfnIEFUzJRD5D4swo10zsN0mvy4Mj8SwSZimxl4du\/eDTHwZX9G7nPCgSecEEV2FmC0TuXtPMtOnZrOeQQ8Rs5jpJ2E0b366qvxwQcfQOc84oDIuOHsT7iMcDaHj93pHrt6\/ycceOTLJGzdzkKsSXVi+YXbfQ8652Ebo5jOeG3cE+nBVui7ZdzzsC1dRmiB0qxZM+X6To9divcXXXSRstaWPQ\/rGsX+WKZlWHse3eHJ+DKiHXLK7suOpXqx\/MJjiU52xhLLtAwLPHYmbVYnmgl9qZljMTMJCuU5RotjiVtgdd3O3GP5hdsZfyzViWVahgWecDjPsUroG+pzlixZAo\/Ho7RKetQcLiCz61ayeXWR02MJGHbGcsKBx2rSVGEaE0pZ1ZVD02uuuQbTp09XbtLUzIhak+IfY7pde+21KmAHFzcttoOdFxmfF+w5VSXz1fsyZpmr6sXH8gu3s2BjqU4s0zIszmNFXLs5SdnebpLdSHOS2n0Ox2Slgg9VNR\/LL\/xoAiNunhMBdXWuYBY6VroOJcluJOAJ5TnGZL4yVqvrcc5jf6HEzXM0WlW15wkW7T6UJLvCDcJJ6BvKc6LFcYREJwLniZvnBP84RFVsC\/Y4sz1RMLMdI+cJVl9EwkiT+YYigprts4xx32jkfaStuDnFnCaXQ7DOQTjt6Xmhe1\/EzXOCrWbgmIHHzMTjaCT0DfU5Vsl8qWnTzUZIymCeslVxnkrrHNi0zoGJdQ5qBX+X\/hrhtOfYaFSrzyFunlM10Y8ZeEJ499W+qqlhaDXkPHHznCiDRw\/wLha3usmFHpmz2qMgzAmcKHseAQ\/JEDfPOXIxhMx59I24URMVqo9PmGsz5pudCOCJFSLHMi1DAo9xE208zY9kkx0rLysa44jlFx6N+R3LPmKZlhGBh5zmk08+8UcIjYPHu6xi+YUfy4UfjWfFMi1DAg+JIWKbJJjSs6fxHktVsYmjQdBY7yOWX3is0844vlimZcjg0Q9IxQJZz1tT04ET5zzRhecJBR79HEDcAnRyxf15aq7YFrdts\/HhCCdum41uT5gqsfy1PB5EPlFt2\/4fL0bPD99gOucAAAAASUVORK5CYII=","height":140,"width":232}}
%---
%[output:12dbb617]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:6baa8017]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:58c0f89d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:4d9e7c12]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:83f879fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:0c55df08]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:9ce16302]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:1c656e10]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:1594056c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:6a0f39a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:57b0d064]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:333aac54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:671e0000]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG0\nTaille de la serie    : 205\nStatistique T_max     : 4.6931\np-valeur (bootstrap)  : 0.6430\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:139bcf99]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:7e199d6d]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : BbsFG1\nTaille de la serie    : 205\nStatistique T_max     : 3.0202\np-valeur (bootstrap)  : 0.8770\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:880f0b10]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAM8AAAB9CAYAAAAFgzchAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl4VdW1XiEhCTMJQxISpoqIWBWtfUyioIiKcxHFYkUULBSVp4BTREREEUURfUUfWMBPRcVnW63UARQqVahVsDhQoMqQkEBCGAIxgRvy+HeyDvvu7DPde+7Nuck938enOXfPe\/1n7b322utPqKqqqqL4Ex+B+Ai4HoGEOHhcj5mvMuxb9giVbfiQsnKXU9m\/PqI9c66lRk1bUVbuB5TSvbev2lrfGhMHTwzPaGVJPhXOuZ4yJ70herFrxiXUfvxC8f975o+hDlPfp8T07Bjuob+bHgePv+fHsnUADwPmaEkeHXj3GaGBjuZ9HwdPFOY1Dp4oDHIkqzi09v+ClmoAEZZu7Se9Rc37DItk1Q2+7Dh4GrwIxAcg1BGIgyfUkfNBPnnPI+9tzN77oMn1qglx8MTgdAIc+bn9KFC8w7T1Kd37i\/1Po6YtY7CHsdHkegmeN954g5577jlatGgRnXzyyUEzsWXLFho9ejTdcccddP311wf9VlJSQrfeeqt499JLL1F6ejrxuz59+tC9997rq1mNa5i6nY44eKTxZ6B8\/fXX9Pjjjwtw+Rk8dSs68drj4DEBT4cOHYTmatOmjdBGftQ8aPqxsoNUMHMoVWz+u9GT+JItOsCOg0cDHrwqKiqiK664gsaOHetb8DBwEtvkUObdrxs9KXx6BFXuzYvveSKMoTh4NOCBlsHz7rvv0lNPPUVPPPGELzVP3NoWYXTYFG8Knt2lFXTn65to3ogelNEipW5b6bL2cA0GAA9rHFkL+c1ggLZBy5R\/\/b7hy1axeR0VzLyYUs+8JEgbuRzCeHIHI+AYPBt3HaJ5K7bT3BGnULPkJAdF110SL8ADoKCc+++\/X3Tkt7\/9re+sbTzCcA7d9+Y0Y8DTrptOacMfqrsJaCA1x8FjsmwDeGTrm5\/B00Bk1XfdrLfgYY0hj\/hrr70mrGc459m1a1fQZAAcOuPA2rVr6de\/\/rVvNU\/c2lZ3mLIEz5X\/s4G27Ckzbd3J7ZvSOxN61doTYb\/EeS86NZ3+MOq0oKXe4SMBumXJt\/TR9yWi7NF9s2jO8B51Nwp1VDMLfqsr7grJiTNubaujiaupNiLWtknLNtHJ7ZvRb\/pmCZCMPTeHBp\/axujpC6t30pY9hwVgGGij+3agced3rNvRiGLtsotNqB7QcWtbFCdMU5Uj8EDYP9lcIjTI5\/85QNct+BelNWtMfxrfi07v0DyoWAbDzKu6CcDIQDHrKoOtoYAHGmPPC7dR2pV30+45wyl91NMhaZ64tc3n4JFN1mjqqEXf0JPXniJaPeWtf9OS0T8PWrYhPacBsGTg6ax0anqUu6OknALHqiirVXJERmfokETato0ou2MVfbTymFHH0MaJtI2Isiur6KNjJ96H24jGjRpRUmJCrWJY+4QDHhQajrUtLy+P8C\/WnpycHMK\/unxsNY8s3LsPVNCCNXlCA\/1QXB42eHjvIy\/rAJwJS7+nnlnN6LYBORQIBKisrIyaNm1KSUknTOThvL+wdwvK35kgwLNyXakof2yHdPqySTVYAZ6Ve6rfq\/Xid9Q9snUqFSY3ppwakOU1akSZR47Sgl0loiz5t92NE6kLEX0izTQusRXNuZYQfaVZv+sp464THgKcDGc4ZZ9XX7GORFwCgGbKlCm0bt26upTBkOru3bs3Pfnkk3UKIFvwoGcrvt8btFQDiLB0e3PsGUF7GaR1umzTaRzkX7N1P135+\/X09rgzqXfXVlT+UzkV7t5NWZmZlJJ64rA21PejO3WkdQWNjQkDgLIqyoWw59WAE4L+1cHDol6kx3u8wwPNlBMIUH5iQnV6vMBTkwAAQlnqg59\/rHnJ16fbjnyc9rw4lqoSk6njY58FxRsI15jgRCLZkgghzM6OnVgHAPuzzz5LsJ6yN4iT\/nqdxhF43FZqZzCw8l5g8Lzzu7NoarfWQjYzystpFRGlpqYaTRlERFsDAeqWlBT0RS8vL6eCggLKysoKSs\/vz+tYDQadcDMOGAvQMLq0QXmRCf8GBpcIsKjlMXigdUqW3E2ZU96mwtlXEyU3p\/RfPxq075HjE1hFwQnHVM3g0Qkh2o4xx7OdiE4cwbqVBu\/TW7Xb+9rMS4wIeGRTNZuhZcDMfv9HWvR5QVCrHruqm7C2MXjmXncKTevToVoAt9V82Luc2I\/0rNmfsKBjqdU4MUGkh9DLSzyuyBEYQhl9VMpoYRWllCNrHoAHwTra\/24hFcy4SICnxYAbgrwCODYBF9O07\/W13G3CNVWbCeH045U+rBkHvPMDiGIKPOq5DMZVd37D4213zsPpdHsegAdaJaNrK7E00j0Aitlvoch+NPKo4Cl57UEK7Nokqk5s04laDr7V1KWGQdK015CgNOGaqnVCuBjnbjUDcnPNNwEaiLXQouPfMrznhy8XyofOVt4Y8NrYsGEDXXDBBYR7U9hT7t27l7Zv3y7uT3388cfUq1cvUTynU+cnZsDDAp7TOiXoIBNLs7z9FbUOQNFRu2Ub0sgAk\/dOAA5PVMSEmjVFxCqoXbAMnoMrFtDehRMoc\/LbRMlNqWDWFdT2lrnUcvBY0xbBeIBHvnqAv2FpK\/34JcqeWb1nYgteiwvMwciVqEKIYela86OqZTAnmBs8cohZgOeDDz6g22+\/Xfz2008\/0cyZM2nUqFG1bvHid9TJQHn++efp4osvNtIh77x584Snx+bNm410MQses\/2J1Xt4F1id8wCQD\/15K93SP5tuefk7Iy0GCZMn7xUiId9Cc+1MMDb5rupQG6dbpkngzCmqzpDUrothMAB4ihdMoKwp5uDBsu3QZ8sEWFSTtpMYBkltOxmAMuufCh4GCLZvsmWQ8\/OHTdY+VuCBRoFrEx7Z\/QkaZ+rUqTRjxgw688wzafz48bRv3z6RDm5VQ4YMEfvWjRs3ao0CMaN50CGc1Sz6fJfhimPlFeDmnEe1zAnwbJPAY7J\/sBR2SXCT8vIoIJ0F4G9+co5VUV4j\/bIQadR8\/Lcogz\/PSLclR4BQLlvUUdOOzJQArf5dV8p6+BNqclq1VUEYDLBsK6hetiVl9RB7npYX3mpEAIUWkU3Vuj2PK9BrEqtCyHsds72N7nfdsg1X2C+66CKaPn06TZs2TcSCgIc6ns6dO9fSPLxsu\/LKK+uH5pGXVWaTpPNtCxc8gwYRrVlTI+T4\/GkAFCTIusZZAKjjeec5krmC114zBZBchpxOV3ZGcoAeuuw0uui++UL7MHg4uif+xjVqdT\/jqJFhJlLBw\/sd7GmgXdSHNQ+mhY2LqubhPGZ7ofPPP7\/+gyfUeXF6zoPydZoH7+1Ove\/v08dY2i1duzbUpkYtn3oazqZq7FMSUlsI8ITqHBpOJ3TLH9bFMkBQBwNL3rvhvRl4YBiQNQ+3U7fnqXeaJ5xJcWIwsAJPOHXHQl75DAftraug7DrwyGZqAAVaCIoc4MGjLunMwIO0XD7PCZZzWLZhH4T\/37Ztm7jqftddd4njBSzb7rvvPtqxY4fYB40bN86IYiTPa0zseeSlG5Zof7ipp9jg45qCmWOoDAqk053z8LVuM80TCwAIt43yOU6oXtVog3xIiqg57ccvoMI51wm2BDuKETMhlM3V3E8GUvyc58TMmx6Ssol6UPd0cXgJbfKnfxUbntR2Dp\/hCldDze8mkKF8SJp2+V1Usmy68JE78N5cg7PHKmKoEw8DeBd0rtnjhGK\/icQ8+l7zqKZoNYaBFwFC5BgBMFlylE5Z3XPwQZ4EnAVAtd9www2GX5NZeq\/Kj2YsAwZP29FPU\/GiuwX3jhnHjgy0QHGeAZ6q8tIgq52ZAPtFCN0CzC\/tDvkadrjgwVr5scceozlz5ghTJsI74ZT6nnvuoQcffJAeeOABMaZyGjmmAPtj4d2kSZNqpccm1Ivy0U60afbs2UYARMQ3iKRDomyi1gmWfIbDMdrSh0+j\/e89Q21HzaGCmZdQ8sl9baPn+EUI4+BxOwJKekzk0qVL6Ve\/+hXNnTtXaKEmTZoYWgaaCafP1157LU2ePFlEsoEAIx+Ap6ZXhTvU8mVwIv4BA1WNgR1m943srE3a3fosFb000VLzcCbVD85p9BwdeGSHUF2fYKJWl28Y\/xdffFEk50irGB+8x8Mhu9gCh4\/RK6+8IjwJ8OHEg4\/UkiVL6JRTThFnQ\/Ijl4n3fgF9RGIYhCJIGOguXboIawxANGvWLFEMlmj9+vUzgrKz9pHBY5We2xJq+civCwAfSh+d5HGzbHNSnlUaO982XV7Vtw1lrF692gAIgyA3N1d80FauXGl86JyAB\/nw0UTaBQsW0J133in+lh\/fgyfciXGTH3sKmC0BCNYQXoIn3PJZs6FPiFvtp2Wbm3FW01qBRz2fhkbCIanOMRQag4VergO+a+ecc474GLI2wdmPleaJg8fFjLJGYLoPu2WYTvNYLdvCLV8GHrqlLkVcdNVVUidWNysft1B829BANlPj7pG8PGOnURU8vOTicF7yEosdP7H\/hHbCMo3Bgz0kznjkB7HB8dGMax4HogJBhLuGvD8xMwDw2lgFj1V6L8rXgVleRjropqsk8rINwUF0BFZ2wIARoXm\/4bZBRaw0jxvwyB2UPQtg1GGvaczF6aefLjyw45rHlUjUTqyePiMFf3ngdcveuOotRxU8yCeXxem9LD+apuowh1Vkd6K15HGTx9it5mGHT145yNcKZPDI+0ZoJDODQXzZ5oUExMsIeQRk3zmzMyI78Oh8cuFMrlu2ydY2lMtgVO\/rsBaPgyfkqY1ndDICzHZwrOyAkVxdrlnteZy4\/HhhbXPSF6\/TxK1tXo9oPSqPz23U8xqOz+YEGE6Gw6tzHid1eZmm3oGH3WZgQVEPtbwcOKuy\/DKo4fTXLuQUBw+xY7rm+UBb2IKltovHa\/HixfT5558HHVqG0wddXiftcVqnX+bZs+g5fPkJm\/66IoHyy6A6FQJdOrvNPn4veOIaOrZ\/FwVK8k2rCqR1oucTh1BCclNb8PziF7+g\/fv3CytYixYtwmm+ad4jR44YXgi4kp2cHHo02Pz8fBGssc7jts2aNasKwi5rDoyA7JBp5mDJIyX7nOEdAAQ3G9j+MVBw3eCIKvLmkq1ryANPAjyYPAwKtBf82+CfBp83s4gsqiUMdcl+bzjUhPVO3sji\/3Xt4PMFJ3lkR1a1\/\/LYmf3GH5sRI0aIU3i0EX2cPPZGrVMnz8\/2H7bSzO4l1CP3rVoOo3JdOALA07JlSwM88lhh7O+44w667bbbxAE1nqqqKvEvIaH29XS8b9SokZGO0+A9P3jHf5v9P9LqyneLWF9EDGXwQJjw5YZLBSKXwFQs89lgoHEv3YwZWtU8rAVkoZfrwGBxWXDBAHi+\/PJLwUCNB8Br166daM9HH30kAkOoXxq1TgYE0kG4GZBYtrzzzjv03HPPifI5MIXaP76MxUsdNQ+7AeGmK9rHWlbuF9rK9ahM2rJmZJ4gfPXl9i2cO4ta\/PG+Wj5tZuCR41QvyEujT0qa0WPdCikpLYse+CGbfvHL\/xLls\/kffe7evbsx9iNHjhRha7\/77juheVJSUsTHDjKAjypuweK3Rx99VDjsduzYUTjJ9uzZk4YNGxaUlj9K+HAOGjRIeBZ89tlnxkoE44S6w9U8qMcXsaoBHpXUib9g+CqeffbZBhmUej1A\/lqYgUfVAjIlO9Mfzp8\/nxYuXCiKw0RDUGSQmi3HVPpEnXAC9DiD4PaBoBeP7hyJ0+jy8EGurLFUbSprI9Shtlte948ZM0bclmQActpXXphHOX+cYByOHjjaiHK3ZtCmshPRUlF2h8wMWjjvSWr26nhx8e3HZbPpnlXF1H\/wUBrVeD1VZXSnWf+ottRhTOFUyx9H2Ytd\/Y2dcTkf\/pY\/DnLehx56iB555BFx8xMfOfmDqMqUl3set1oqUukTAB5eYunY0vD1Ue34OhDZgUcVTHQokuCRASIPHtrOGkZ2DwHIrfLgN2g\/HYOcGuyCQSRrIXgZ68DDQLXbr3FeFlQIsbw\/2v7O\/9CdL7xLA665USz9tj1xHT17sKfY87CQq+4w3E44YDKwVPCg39Dgurw4r4EW0n30GgR4pj4ys+qblgPoyMqnaED\/fpabfZ5gnTXNDjyszbzWPPJyzkrzmH19eB8AQcKyZOLEiWIvwCfmyKd+NVXNyGXL+0aADPsOXv6i316DR44kWtGlL90ydiyddPa59MCFP6PS1UuCDAaq5pHHQ9YsTjSP2l8r8HTq1CnIQ17WaJHSCNEqNwHg+XfmZXRK4Xu08asvaPLMefT0sr\/T+iUPiD0Gf415+aWqcG6oHXiQzm7Po5sE9rSWhdCsTqd7nq+++srYl2DtzHsj\/kJzO+Q9z9tvv218nXkPBoDwfo3z8H5It0e0Ared5jETCFn7PDppPP31s\/Viz9Nm9FyaMHcp8Z7KbM\/DKwszzQMwyW2T90tq39WPCuYDGkvex3J71GsG0RJ4L+sxwDP76s709KPTDPX88IyZdNPIG0RdspUGf+tMhE7AwwDii1O8Z7D7slsJluzDxpY9t9Y2WZOq1jEuS16aIT2MGfxVZaFhq57qHayz3qnL2FDBIwuDrPmcWNv4\/EceQwg6PhSqhlAtdfKVEbOPnl17vBTkuigrofBgedWdr2+ieSN6CIY3NVZBXTQqXmd8BGJhBAR4IsV6jQEwY0wASK+ev4H2HT5qybgQC4Ooa6Odp4CTfqnXqzkP+7jh7\/zcfpSY3pHsPA6c1BdP424EwvIw4OCGxX9\/hf6wcEGtmrGM2vOzq2oxY5+e0zyItxTlXPrzdrVY5tx1xT+pvWK6BjDs+Ep15FaRiGvtn9H1T0uCwOOW9fqKZ9dR1\/+8RpktU6jbZbfTtgOVBg0J1rv\/PekeWr1lHz395BM0tFe2wYwNoPz1m6IgypJIDgmuD3P4al30\/3DqDkgsCCjHK6ZrOzcdqzbzlQSa8DbtOnQ0nO75Nq8vDkmravwpQmG9Hvn7zyj7P2\/RM7Nm0JsbDxt088x6\/e8dhTT8tik0\/+nHqH\/PjgYz9k19sunltfm0bW+5iD4qE2WB0HfpFwWCjxT\/KgOVdOjQIWrWvJnB9ub23S\/TWgp6RFCLbDpyVJQ3PO0qykvKo+zKTvT+kY9CruPQ\/z1CRzf9jdpPWGQEcofEecF0De+Bip3fWYaQ0mkekP8eG\/MK5f7+1Zgk63WCWF+458jgYQp4N6zXPy9aTm+\/uoguvnkKlWT1M8iu2DqTds5VtHjuDJJp5fuf1Jpe\/UehoCxpnpJItyz5ljgyKdMq3tKvA43ul01HjlRQUVExtW3bVriO4HH77t6eXQyS3eX5++iWtkNpQ\/N\/iLIyj3ak5fnrQ66j8T9epaRXx9WiEAmX6dpJfALe8+C\/TG7FgscWtFgj63UCHF8S+obKet05qYSuveE3VFqyx+g7zLlPPf8iTXz\/UC2iKyzbHn\/\/B\/pqZyklUAL17tyCenZoLpZxMnhWXT6G8hO3U7ufsmhJyXJKTq4GT0VFBRUXF1O7dm1t393cZigVNSmg7MrOtKTkryLvZdlnUmEyPJK7CG30zo4vbcsblX6p0ZbExEQqTD7B9dN58076hKr5d\/zGdF3XnsdOwOA2jRdmfbd16tJ7YjD4Td8soT3GnptTa9OvY0zIbJ1MQ+d9RXNHnEqgaxzxv1\/T+PM70pQhXQV4BtFAyui6Xwir0A5HsumrQ99TSkoqXZk6mLbRdqqszKHOCatp+ZFjxrt2P2XQ8iMraHiry0Ua8WWuKUN8nSsRdRnv4DVc7UmcXXkefXHgL3Rl6kVUkFJojBHqTKiqrHnXmfITd5hy1uUUEW3eUw0eL5muZadPefLMAoAgPZ604Q\/5JjCgF0KqllEvwKNjvUZHZeuZLg003MQ3NlHBwSNiXHJaJ9PTx7XO4FPbCPCcn34qHUs\/IcjBgl8NiuqnWnOoAJH\/tps8ACrnWCda1\/hTi6QchElP+KiCxwuma5mGZN9fnjGi4SAyTkrHnlryX6\/BoyOoYo915g7lqEZmg6f6RcoH0vLhLPKrh+bwTlCdbZHOl+Bxy3rNAybHNrO7+4M8AM+CNXlif4RH3vMMOq53VkWe0tcOUy5\/70ID6UfB42lQJnrIdH1w5UuG4cDKChcJ8MhkvRwNFHw68IKXw+WaDZgaBITTAQB8xYNddTjGHqeBfyHSyRFJfQmeUFiv0RGA5fXXXzcYDuRQRGYBAq3A05W60raIU\/q6xIaD5Awfr5iuZafPpmdeTAVzhlPWpGVU9vUHQezXOmsbeHp2DHmERo4eE9ZtS5W4Cn\/jH+7qQPPAC1\/WDnwPDMPFGkMOP8XDaMeYLQNOF3bXd5rHLes1BkIOOMiBz2XwqEwIPHgAT+6ftwZZ23i\/dFJJX\/rBWLKFz4udlJdUS\/QDOQFTOCA9\/179\/zmO+Lm7ENiuf6Rwma7lhslaBtpn35vVAdA5AIjMz8MU8\/xu\/Z4KunPlXlPwODn7MuMVZYdQpnnhOS8qKgqihkdbdddZVLJf1T+xWbNmRjk6ekbfgQcddcN6zeDR8U6yAJjxUsqkv0g75a1\/05LRPxe+dTjnubCsWOieQM4AIbgMAFnodaAwQ0TH8zoG\/VTwWoEBEBkonCgYQJ9SeR84yFa3QyUTBts1E\/UCPNFkujZbwuH9X++92hI8iMHGB8eIDqp7zGjisZxavny5sWzjdIgOymF3Q9E8XA7awpFGfa15wvFtY69Zs\/CzapxneYLYLI53b449I8hKx4S+9\/e531jCLV27lOS\/kQ\/vwn3UOnTlyWmgX1gfgkwY4MFpd10xXeui6cCo8M3RdLrtxfdtNQ\/66wY8WK7deOONhNu\/o0aNIqw41KihvJzHf3XaCO91ex6WF7Bl4ImZPU+oQsiqHfnhzo7B1L0LtfxYyxcu07WbWNUYG\/i\/6eJZf12aSpM2Z4W952FNwvPA1raZM2fSt99+K2Ij8BUMAIuvm1hpHi6rXlnbwhFU9c6PVbyDcOrxe14\/M137feycts9Xe55QWa+ddrahpQuH6drKLcfYk7XtVMsdRx1jvwhYJObeL31LOFRxtEo+Z4kW67VXZLuy6le1nZs61Lh1csgsqzoiIRxcpu5AFMA89Nmyess36mQ8fQMeu5ukVsS9ZtY0uwGINJkvTr3d1oG4BMxOJ4fegllVRxhsd7JuNwZ2v1tZ0grnXG\/LVeoXAbPrZyi\/+6Vvttew3YBH5qN0E+CBLS+RIvOVrTtmdagEwDKHqRXzXCiT7zQPs2JzwHf+28llNy8ErK7cc2SDgupdwHOpCwjjdFy9SucL8IRKtiubOzEgKvmvPEhO60Ae+fAX0UWdEAZ7NSFqOfL+Cfd0snI\/oJTuvW2r8wo8deGeAxM1L7l1YZa96JvtADpIEFYMA3XZFormCZds14lgu6mD9z58eu4GoA7GO2pJvBCwunLPgVxt2LCBmjZtWsuvzVeahy\/DhTKrapgmXRk6r1hOFy7ZLgTEbknlpg6du5GTOkIZu0jnsQMPHHDhQwi3ItxF0j115Z7DLBs6p9B6A55wBMALsl078l83daAvOlcjuzrCGQO7vLLTJ5w9249fQIVzrhOxqe2WbnbgYQdc9skzA4+8bGOHzki759x+++2iOXHwaGZFPVlGEi\/JfOWvk1y9VR2qAyPyselbbm+0bmbKTp9pl99FJcumU8Zdr9OB9+ZS2YYPbUNN2YGHNQ\/6CZ88p+CJhntOg9M82FeATkKms6grhji7L3os\/C6bqgPFeQZ4qspLtfw9ap\/swONkDKysbZF0z+G2NQjNIxsO0PG6OBdxIgyxlgam6cq9eZQ+fBrtf+8ZajtqDhXMvISST+4bPySt4Y9SjxiiOcdhxTDghpqZduEdDWa3OXPmUKQPFaM5aNGsS40aqpL8mrXFC80TzX66qcsvffMEPPIeg5dqMvtaXX4d3ExKfUrrFwGLxJj6pW+egccJAVYkBhJl6oizrL7I27dvF3dFzOhSItXOaJbrFwGLRJ\/90rewwcOHivALy8rKog8\/\/NAgde3WrZvgAnXjqhPKYDsFj9N0obTBb3msBGwxEa2WGjyKiAZqOmB2zmPGds6HmxdccIG454NDTqxA+GP18ccfU69evURNOARFulCeegMe7kjfvn3pjDPOCGKWw5d9586dYnwQKFDlG2XqvdNOO41KS0uNYBJ4z7Z+1VzMpmKZjxTl49KWyiXKNJEog8mImUMHbh+83ASfpkwtiP2ZXD4TYDG1oNldJVXYuK3qZMtlM9EWOHxmzJgh2nT7KZX07Z4y+mRfC+rRtJxmdttNrRofM+TMLG6bLIhmAgbgIHhXdTSE6md6zZXsairlE4\/ZNWy+QaoKPupkoKhRc\/CR5XBVCBTC6Ro0eACQNWvWUHZ2tjAOyIYBfInAtlZZWUnJycmm4MEAyizc8s1EZhZT91Bm4GG+UXav0bGTMchUBvDJkyeL27AyW5xKeCszhct7OZWgS6YxVPOo4AE1JHOl3nv3nfSXD1bSUyfvotZJlZS78yS6euRoS7pLnQDqwLOqRuPIwJEBBPJ5WQNZgYfnA\/llnlZ8nKZOnSo+BPAuAWnxvn37RDXo55AhQ6igoIA2btwY8i3XeqN5WAAh8BB0gAAAYtcdeDH\/85\/\/FINnpnmYp1RdVslC5hQ8zCWqMrnJ1H468KB9+PrLDOA6akT5uoKOtxRjIDPDsXaTvYB1moevsMuAPvZTKY26agid0biIxuZUC6ATraOrE+9G46q8xade\/V23bGMtLntjcAyDzp0719I8vGzDRy2ueaTBhyBdeumltGdPdYzqhASEsQ1+Tj\/9dIIFzkrzeAkeXp4h8AQEUuYfFQIkLe\/MKNKxLsdHQAacGVO43FvVc0L1aNAtOeX2If6DGW8r7zNw9nNky+ch3SS1Aw\/CUclebqrm4b6a7YVA5agu2+Lgsfha2alQK75RlW7cC82Dr5\/ZV94OPOgK94JpAAAIq0lEQVQLlm46vlGZxdtunS47zOJLbdUmK\/BU7iugm4ddRmc0OWhoHvi4OWGB082LFXgQEQi\/OwGP2SVI3Z4nDh4LaWFw7Nixw1iyMQjw5VXZlnn9L6+TnWgeWfB5X\/Pll18GaQcss84++2xDu6hayA48qtDzskzdG+kuYnFe7ossvG3atDHapLZdBc+s6Q\/Sp398hR7tUm1omXb4XOo\/eKgv9jxWWlb+SOD\/YX3FEhaheQOBAHG\/ISfYB40bN87wHbT7GKm\/232w3ZYXavqwTdUsNDA9yss2cGbBiKAufQAoDCCEzI3m4SUNQhthGThw4EBatWpVLfDw+Q2HQBowYAB9+umnYnOKDSwuzPG+pH379mI\/xvs0DKJsYOAoqGp8A6fWNvkiF5+DqW23A0\/u1gw6s0W5oXmc7nvMBAzaBVHR2GgAjbOkRnp0hoRQBSuS+eoFeNSgh6rgyc6ikT7rieRkxWLZVgLGVjcO3ggr280x1Ml6AR67m6ShBgiJoXn0bVP9ImCRGCC\/9C2sZVscPJEQDW\/K9IuAedOb4FL80jdfggfe2NhkZmZmRmLsI1pmUlKSQTzspiKOz+0mj1Xa\/Px8mjJlCk2cOJFAflufHru+RYspO2zw4GCRXV50E2QVw0CXHkKESe\/RowfddNNNAkRlZWXCTwqCiUf3zux9uGnd5EcbWrduTWlpaa5klfsMotr4E\/4IRIspOyzwhNNN+UBRtl7x+5dffpnOOeccglFi9+7dQgulpqaKKnXvzN47STsrNZVeXZNEzx0qpR49zOsqTE2lVYuJbr7ZvA2haB7uc31krg5HRkLJG02mbM\/Ao56uc8d1mscqqIa6nt1UXk4pBQXCY5vBU15eLvyj5HeoT\/feLO3WQICeLy2l1HVN6IlLUqlLDVnNY0sL6Kyz0qiwMJUGDjxRZkVWFl1amCroRQaOJpo\/v5xSUmq3IZQJ98saPpS2+y1PNMfSE\/DwWQ8ORJ1cfEMHzUJGcecnPPkkfTNsmBDW2woK6Nt27YQ5NScQECApLCykjIyMoOsOuvdmaS+FFkF5WwOUmxMQZf8mEKC3iptTl+lEGd+X06pV1eB5obyc\/tS6NXWrWTa+Cra5JUQ5RwO04NyAWE6ee24121w4midawUX8JvBeticmwWPFEKcODjpoFqyQO\/\/jxz9QzooAXd6imNZe0oQ2tGpFvcvLaV1qKj1eWEj9t2wRbjQpKSlG8RUVFYJMSX6ve4cM\/6msFBpNTfvyjkZ0MD2dJrUpFeUi\/91paXRBo0Y07NAh8W5Ku3YiP57iL5oJv5Zhw0pp4sT9Ie15ojnhXgpqOGXpGN9QHs4G4c7k5COsqz+aY+mJ5uFOc6B0u0F1Ap6THvkbbZ2aQzk5AWp7ziHa8KfWtG3b8dPxh4kyepfTuMJCuvrq1mIpB+306JokWr2KqNeBA1T0RTPKz0+i\/iMDQlONGZNH5ZmZ1KNmz+Rmf2S3Z0K5zY9rq8WLiR58sFoLsWHDbhz492hOuNM21VW6Bgse3NdQH92ex8myjZcwmzaVC6G85JJUGjy42tqGd82aFdEvf9nOAA8udGEfAyCBO1SmAsbSLK\/b8aVWV6KqquoWOt0fudlfhSpwsQ4etB80i+wGBc9qlgU2BsnX9NmvERcQQQ4May38Bdl16fLLLxfe2ayB2Lsbll12rzILaxbNsfRE8+jC1FoJkhODAZ9PwFS8cuVKYbrGYOLRveP3S1avJlz\/XtW5s\/DhAohW5ORQj\/Hr6MZzs+l8+KKYlOG2LrVdKDeUMwbdhEPL4p+TB4YN9cF+zemjy4+8VmWANrSGOjQosqcarxygOeuss2j9+vUCKOymxcs28JvOnj2bpk2bJn4DOHCRUQceLLHxwH\/RyrM7WgwKnoBHvmLrNMSUbJ2TN8qxfuYB0OOfm0cHnunHVenDDzsrhbWpnFpztcq0MF1+JLYqA22bVuNJivbzPR6d1RXaBw9rI1nzDB06VOxzcnNzBXh42aYDz+HDhw3OU5Sn0z4xp3nQETRax6XibPqDU3l92h5KG0LN01A1D4PHjimD411fc801tGLFCsGszZoHH16OL47yWMuwbMks2WbzE3PgsWJLcOthEKrQxnK+aE54JMZJ1jwoXyV3fuGFF8SeSA6+guso6p4HeXGFBPd9+A4UbvBChi688EJxTZ73PEhrtp+OqWVbJCakIZUZ6+Dx01xFcyw92fOEonnckO3y5KjEU7xcxJcGj9klNTd11QWxbzQn3E+CHom2RHMsPQGP2SDo+HGQ1i3ZrhyNB6qfDQxOuHPc1lUXxL7RnPBICKyfyozmWEYUPE4vw\/GhaSQJfXmC7eqqC2LfaE64nwQ9Em2J5lhGFDzyYaiVCdsp2a7qQ2flqWClDWG1wZmRHZ+pGfsDyrYiD3YrFNGccLdti7X00RxLT8Bjteexc3Z0Q7YbLnjc1BVNYt9oTrhfwBD3bQtzJtyQ7aIqHXjsCH25iW7q0nlMWLkUhTkM4owsWubVcNsa6fwN0rfN7aC6IdvlJZ8KHicGA7TLTV1IH21iX617To1rkZNx1TEcuPDO0TIkoF6rMuD6hH9s8Yz7tjmZKSWNGuBQDsVqdkCqc+HwmtBXNhCwKZvf+Y3YV+uec7yxDr1zqMbXNWhmagc+Np9cXX6ktioDbeM4b7J3Sdy3zQWI5OWQuk+Q9xguimxwSeuD5on7trkUW13oKZmD1Kmp2mW19S55rO95ZPecuG+bQ\/FUwaFGCI2Dx9lA1ifwoMdx3zZn8y424zg34aDq4LeRA6SjGDMaPodV1PtksQ4eP01QNMcy7HMe+YyHA5urjAF+Glw\/tiWaE+7H\/nvZpmiOZdjg8bLjDbUsnvD6GN0z2nPK0UTtDue9aJdn4NGZn9HA+H0e+2mK9duz9j2MboqYihjqNm5bdIcyNmqL5duzfhvhUG7zhtIHTzRP3KoWytDH88T6CPw\/NVvTAL5h5VQAAAAASUVORK5CYII=","height":140,"width":232}}
%---
%[output:7ca75747]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : SSAG0_homo\nTaille de la serie    : 203\nStatistique T_max     : 2.8508\np-valeur (bootstrap)  : 0.8980\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:9365cd62]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:08489f7f]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : SSAG1_homo\nTaille de la serie    : 203\nStatistique T_max     : 7.0442\np-valeur (bootstrap)  : 0.2930\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : SSA30_ae\nTaille de la serie    : 163\nStatistique T_max     : 8.5234\np-valeur (bootstrap)  : 0.1610\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:7ac8d054]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAM8AAAB9CAYAAAAFgzchAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl8FMXS\/ycEEkjCGQIhAYKCF6CIvI9TTcQDFUFERHxoEkARuURAQZQQAQ8ggqAPFXghylN8CO8AfRxqAiiHoqIgIiAECIQzXAECJOSjelOb3s7szuzubLJLtv3xMzvT3dNTXf+p6urqqoCioqIi+IufAn4KOE2BAD94nKaZxxqcWPQazm1eiahxX+Lcr6twJPVRBFargahxKxB8XVuPPdffsWsU8IPHNbqZ3qow9wAOpfZG\/ZGfib4PTuyCyEFzxd9HZg9Ag1eXo1LtaNOf6+\/QdQr4weM67UxtSeBhwFzKzcappdOFBLqU\/bsfPKZS2rzO\/OAxj5Zu95S3YbGNqkYgItUtcuTnCGvX0+3+\/R2YSwE\/eMylp7+3CkQBP3i8ZLLlNY+8trF33UuGXaGH4QdPOU8\/gePAuA4oOLbP7kiCr+so1j+B1aqX82j9j5cp4HPg+eyzzzBr1iykpaWhWbNmNrO5c+dOJCUlYejQoejdu7fNvdzcXPTv319cmzdvHmrXrg2+1q5dO7z00kvlyhl+CVOu5Hfp4RUOPL\/88gveeOMNAS5vAo9Ls+dvVK4UqJDgadCggZBcderUEdLIGyQPccHlc6eRM\/kBXNjxnZUp\/CpbueLD4cMrHHiIGkePHsVDDz2Ep59+2mvAw8CpVCcG9V9YaJ20Q28\/jsLj2f41jxdiqMKBh6QMlaVLl2LatGl46623vELy+K1tXogOnSH5NHgqV26GrCyIf3v3AiEhh7FgwSS88konuwYDAg9LHFkKlbfBgMZCUib\/l+VWX7YLOzYiZ\/J9CLmli4008j02uzpH7BJ4thzMw8OzNyMitDL+O7gV6oUHlxl1yNo2bdoitG\/\/wRWrWx3N59atexaDB4ciObnktmocoH7Gjh0rKgwcOLDcrW08UnIOPfHPkoHXeiwFtXqNLzP6+h9knAKGwcOAOXH2kuj99e5N8eydDY0\/yaSa7767DKNG3YYLF6JEj7GxQFyc5f\/Hjx\/Hhx9etLmXlma5r4KHf5P1zZvAYxKZ\/N2UAQUcgkcGzD031sbfE5pj\/H924f4WdXH3jdpffU+OmdSzJk1KnhASsgHh4YsRFJSNJUuGIz+\/PoYOTUV2dlOcPDncWvGZZz7FG2\/cV8o4sGHDBjzxxBNeAx6\/tc2T3GN+33bBc\/jMBQxbuB0zH7\/BRi0buWi7aeA5e7EA\/dJ\/w9OdYnTBSMCJj7esb7QKSZ4sFN\/MirWpQvcSE2GjxrlCSnbcpLb21CmtOioogiIaIXryOpsjBn5rmyszUr5tHEoeAlC39zZj55FzVlVt55GzpoBH7vufT9\/sEDwqcGRVjQAzH+nAnZlX9LPMEmoSgNITETs\/2Qo4asdqnLNkl48MUFutMzb26lB9Pm5g71Cb39rm7IyUf33Dax4a6vur9+Pl\/+yyjtrVdQ9JHFL\/+nWMRr+PtmFy96alwMNWNPp\/UlIJoYj5SYpQyUQm4rOSgFhF4vBvqkQgSrqy8MmMs3aSkWFZBzlTSKLkpr8gJEZASLjYzKzx0AibowL26lSuHWPoTI7f2ubMjJR\/XafAIw+XgJS2\/qBb1jaWPjJ44sn0nH5FBZtgXz0TdzIykRUbX1zpitUgMw2xq+NK2hGA9hDqiqURAygrVhgXsEef+FQto7gaAYMPqNElAk+1VvfaWMLs1aH6Ri1o7ljbsrOzQf98rcTExID++VpxGTyOXlRWydjQEFolyKYJrZ3mrT+IQAQgtEolfDmsNcIvhqEJ4cHOusbagQCGZDngG8WqGuYnWK4Q96fFl0gmuR7dlCVUqReyQIcEFAHIHfDIpmb2oq6d8LapB9wINKNHj8bGjRt9jQfRtm1bTJ061ecAZMja9vDNEUjtdQNk61uzyGp2pQ4Bo1lkKJ5sH6VpEGBDwaO31cO0lXuF2ta0ap0SSxoxdXKKRd2iQmsZUrsmpAAJ84uB4UBVM4V9LOKJhZQ7apt8CpQNA6rUcnfIbDkkJoyO9p1YBwT2d955B5988onw9PCl4hA8DALaz1EtY6S2kfGAQCUXVRXTqkd1EtK2Ysz9TfDS4p0Y1uY6vPJ0LbGwb9rpJHb\/50Zc\/vUGID7DIh0YSInzgTRpAUTgmp8IZMTbGgv0ZoD6pT6tBgaCiK24izlKnexBUN1YoeG5YzA4\/fU8MSKSQMJrILUXokYuKhURxx1TNYNHiwnJqELrQyp7sRfJkHaP9Wjl4fuOxu3hR7vdvWFTNUmd0Z\/\/gfSkFsJ0bc+UzcCY+uj1aNkgTBgZMnbkij0iVt2++v04HpvzKy6jSKhtUTtvxLZl9cVa5ONVJxGPOBT8owcw4cokk7RpnAWkJNuqaqSiNdljAQCB54pBo\/7FaByqcsCGKNF\/H48Dd31kVdGiCxvjwIAJpUBX\/2JHAaBKlSzqZUDuXqx+DoiakIGqzS3WBdkMzXEFVCuZVh0VFFpmbndN1faYMAUpmHDlP7XQNW8AUYUAjwoCd8BDE8kAKjwdjNPzOoi5JStYUMxJdPvbzzhV2B6XZ4ag\/ti5qDRxkgDPgdf6WaxnJI2yYhH91Dc4sOYa0XZg7kj02fsM1jRcgfERw6y80nbnk9h43UcW4BWvcYa9eBEzp1ZG4sDj2DK9G36sug4PnuqF0bsnoW7duggODgZ2rEXAjC52wWO224y7pmotJpyP+UiCRVInIhGxiBUSiKVQGtLEdS58mPDgwYPWa468L8hLY\/PmzbjrrrtAnhrVqlUTXh579+4VvoXffPMNWrVqJfrieiqIr0rw0Euy2tajdaTY70lq38DqkkP3qLijtpF0+mF9JfR\/tCrC617C6SOV8e0uC3h6ht6KuZNqonuPM5ifXgldq92P7yqtsahyxZImevzfLYDKjMOSM8twww2H8VT0E9hUbSNiMxORFVe8PgooAooCSuaNVL2kNLzx1lGEDfkcQ6s9J+5t2\/knourXR3BIMPJ\/W40TkztbwWNEbXNXDyBL25lv5lk3UNm4EH5Xf13\/NpUJSVVrAotRRZUywsQPi7QuQknAWALPihUrMGTIEHHv\/PnzmDx5MhISEkqd2qX79EwGyrvvvov77rvPWo\/azpw5Uzjh7tixw1qvwoBHtpoltY8SQOFrna6tUQo4TBg9g8HUFXvwt9X7sWFsW7wxKQDvTKmCe\/qewMqPa5UCT0xMAVbt3IUbQ260SB1S1ViKsBSKz8Ds3u3Q\/MEfcEfDO8QXNjZ+DzIzigEjAc46eQFF+OSTHPToUQtVQ6qKy7v37EZUVBRCQkJw\/rdM5EyIt4LHiMHAFfAYiWGg5ZGgx4QMkDjEIcNqcC9pReChOrL0cQQekijkykSFpBF7ppPEefXVVzFx4kTccsstGDRoEE6cOCHqkePtvffei5ycHGzZskXTKHDVSh5XmIHa2AMdu\/tsyc7DC4u2I\/vkRZxfH4v8DU3wYOIpLEurYQXP2NvD8bcX2ggjwpxVu\/D03c1KwEProORiPb4YUImJWXjwuU3o9Zde9kGWMsHaLqjyfixatAk1a9ZEUlyScO2ZtWwWWoS1EOCJOLULgXOftAGP3j6Pq\/Qyo53KhLzWsbe20bqvpbbRkfV77rkHKSkpSE5OFrEfyCOdSuPGjUtJHlbbunXrVrEljxmTqtUHrXfmfJstjAgLPg7AswMqoX2Xs1j3v1Dsy83H4E9\/x49fpOPM2tuRl9cT5ACaf769pStSwYS+IalhAUUIC1ssLuedebSknlQnKDsIIfe\/ibyl74q1T1jLaah79j1Rd89uy45pk2tK9o6eijqBZ+5phcjBacLiZmSfx1P0MtKvCh5e79CahqSLWljykFQi6URFlTzcxt5a6M477\/SDx97k8L6Os\/s8epNN4Bn3n11in+jf\/wwS4KkXfRmHsgNFUwIQbfoFnjuGPn3aWXzTWFWLz8DdQU3RadwuTIgr9jCgPSAyW5MU4nqSlCE17u5vByAzKQG7dlqOUazfsF78fzVWY0y7MULV+3TDp9ah1w8uEJt2BBwBSgPuOXrv7cn7WupPACwfGBkg9JuBRe+8R3K1sAceMgzIkoffQ2vN45c8xdRxZZ\/HCIPI5myq37ZlFZzPrSKsbarPWWamxZtamKxJVcuMQ2xSBvbsgVgQWz2pqSNeA0mDIAZJzEpGZlKixcqUES+AwkzDX2A9021ZGAyM0M5eHS3wyGZqQYcr\/1kcaS2GFPWd7YGH6nL\/\/HxS50hto3UQ\/Z2VlSWOto8YMQIFBQUgtW3MmDHYt2+fWAc9++yz1qhF8jtclWse1RRtdJ\/HKAOwqZrq\/2VPB6z8d7AADgFILVYAsRp2xVo2YW8a5jdOQVaiHSe44k7IcJCVWXxEoXgzlZmGmYsYi77O9H9HRWsPx+j7Gqkn7wdR1JzIQXNwKPUxkS1BL8WIPSaUzdU8BgaSf5\/HyKzYr2N4k9ToPo+rw6FDbqSe2QMQ9ZtyBUWlVDX2cSPzc3qCxQNB9lljlS4xXUguYhxaAxBwSBLRb2Iieb\/D0Tt4KjihvElaq+sI5C5KQb0RC3HqixnWnD2OIoYa8TAg74LGaCzWOHofClfn0dl2V6XkISK4ss9jhHhy\/AAyb1IEz19\/rW1RzwBxMnT69M0YMqSrtTvaNyA1oHq3mzGp8wCHjyHGiMtKxPxYx1Ip6kIUCl4pQNjiMKtKwc8hFYSKuklI4Nk6tjO+rNQKfWvsQ4PRi0zJmyODsuBYthU8RflnrHl7HOXn8VUm9NVxE28YPgznzD6PI84mvfr1119HamqqMHtS6Cfa0X7xxRcxYMACLF9eEva2UaPL6NcvEOfOncfChetx5Eg15Oe3s0iWYkmi9Sz+qtqsh6SK1XOrY0DQAOQOz8XLL78s7vCYVq1aJfR3iqbDcQ4ef\/xxazQeGv\/3z3dA+2rH7L6mkX0ZrcYco612r2Sc\/GI6IhJSkTO5C6o0a68bPcdXmdBXx60LHiNSxN06RLxPP\/0UjzzyCGbMmCGk0Ntvh1\/Z2a6s2XVU1AUMHBiMunW3YfDgm2xVNC1A0fmdzESkxSUI37lDGw6Vel7VqlWFVOvTp08pz14Cd2xsrAAPSaUZKWPQ5fj\/8L\/Kt+HJGntNkzz8svK6iq4ZdQPSYkLZIVSLmFrqG73vBx98IKpzZFWKCU7XqXCILrbA0UdvwYIFYtOUPoZU6AOTnp6O66+\/XuwNyUXuk677weMGgpg5yXJDIHrzzTdFb8TMNWs+jMBAiy4XEXEWCxe+ialTHxQMzqALD59WGmjFax6xMFZiFzh6XocOHWzivRGDjBw5UkgnYiBSN6ucz8WNP72PfxXcaKra5gYJRVM93zat\/lXfNupj9erVVoAwCMaNGyc+al9\/\/bW4R\/Q3Ah5qRx8mqjtnzhwMGzZM\/JaLHzwuzjwxI6tIDAYZPDIzswrFk6fW79t3rgBaixYt0LhxSUgqeWjOPI\/XPiyNZEbaPakbqu74ynS1zUUy6oJHtSSSRCITvZZjKEkMZnp5POS71qZNG\/GBY2lCez+OJI8fPO7MqIO2sjrEX066Rl84LTVKCzyO6quPduZ5qsShvmQjB\/2uW6UQ0269hNapaxHa4Fq3qeTIx83IGsqR5KE9Ldm6xk6jKnhY5aI0LbQOlVUsdvykTVCSTqSmMXimTJki9njkQrHA6UPolzxus4ZtB8TI5NohnxyUGZZqy0YF+q2CR6++\/ERnnkfttHbTuT9i8p0THhBq2xOX16Hw+P5S1DHC7EZJSkaEsA69dI9smwUeeVyyZwEdsmOvaaJny5YthQe2X\/IYnUkT6qk71dQlf6XIQ5c9d9UTkSp4WFrZq89DdfZ58oKZ++B8PnKfvD5TdXgTSGTThdF9JTPAww6fnBhMPlYgg0dOFEYSyZ7BwK+2mc0N\/v6cooDsU+fsPg97F2h5T5Brk5bapn48+COmntchsFJdP3ik6VQ3CWXJ4OmvrFNcVc6VOYPB5XOnrCNxR11ztOYxkkreDGtbeZD0qrG2MXBooSin3OCNTF4AlgeRvemZvBej7sFwzDUjzG72+5i1z2P2uPT682nwDBs2rIgsJQQY2qT8+OOPrYee+MXtuaTrEcbIfV8jHvugqdFC+V3lcz\/OZq\/mjxf1pfehktcjtDnpCTo6Mx4jc61VxxPjdnUszrYL6NixYxEt2FnSkJQhMHG2aT4IJddx9iGO6vsa8fQW8Hr3VVrI6lr+5QCkZlmyT4yMPY6QQMvBPy11kOaJaMeZvZmOw4cPF0EEzSgXL160ehuQj1+VKlXM6NamjwMHDohgjT4Zt61JkyZihggc5CJDNv477rgDa9assb6k\/JstY3STvACohIeHi5cn6UW78WRmpn0Ce5FX1MRS5A7CxJPz5lDfMlHlxay6j0CJecla56gNO6HKaeS5jWxRU8fA98SHJOEpPFD3LL6\/3BC\/\/va7puPoodTeqD\/yM1ysWtt6poWZXOU++VlkvqdSvXp1q+SRacXv\/N\/\/\/teamItoPnv2bHHkedOmTTh58qT1EUVFRQgMtBwwpL8DAiyH4+hvLnSNf9v7m+pyW9PRQ+e5fDViqCx5+OslM738hSPCcfZocrUg8Pz4449CSlEh4FHoJmIUcrCkABDqF0WVZAwIqkfMzYAktYWYZNasWaJ\/DkBB9SiLNT1r6NCh1kNX9HytNmxSppOp1IYlqPxeNFZ+jpohW5aM4rkJT+HGKqcw\/ZOl+GLtJms7ct+hIkseI+CRJf3pFX\/D4Fmfo3lYvpA82\/OqYNTOBpjephAd3vgSz4x6xZo\/VR4\/byp\/\/\/33wp2ITsBu27YNkyZRislX0LBhQ9BG5k033YSePXsKaUIRbWRtgz6c8fHxwoNg3bp1Nveuu+468ZHwhOQhmvlsrGoVPOSOovWV4a82Mxp97ebOnSsYhpiW9GM5Lbs9dYy+pMyoxHClmLMYFLTXwECjxLtUtPZ0uA4BSW3Dm7CyxOIPA1+TpRE9Qx23rPcPGDAAgwY+gzur5qB\/7b3YfDpYMPe0ZgfRqvoF60c5r3ojjD\/SQkgmuahOkXK2ulFP98XuCfdhZl4rXDywHZNeHom\/b9iHtf9agJnPPoRrEl+3eqATvUnSsNrG4OG5oN+qWsdGn\/Hjx+O1114TJzzpIyd\/ENXs4GWx5vGEJCurPm3WPCtXrhRfGAol1K9fPyvzMmPSoJj5PQEeGSAyAUht4mO9shsISSFHbegeST85VBIBir64alAL9ePAa75S4Bk0SEg8Aqreeo3bMqOy1zG\/mwqerLcewzunb8KlQ7sxrHUNvH+sYSm3Fx4nOVraAw\/1TxJcdZmhtrQvQ1JI66PnB49zsNMED4l1irel5vE0CzyyOqdKnieeGIvQ0KkICLBEsqHTpZzMKiHB8jevA4gZSC2hRbIMcGqnfjVVychkkve1CGS07iAJx+pmWYCnUaNGeH38OBya0gPTdoQgqFYUBgeuwZI6PUSgkslND6P5VNtMclpqmyPJo76vI\/DQeGQHXblf59jr6q5tAx5iZDr41b17dxG5ngpP0qjJM\/Fk2lYErXkbj3eNx9hRI2zWJ0bVNntrntTUxZg+\/Rb88stJREfvQkSbCBxvuRh\/Fv6JC3+0s2RJoBC70Zcwb95lfPnli2J8rIbwBMvrpCVLlli\/zrwGI4Dweo3b8HqIAEgxyuypn\/Jay4jkMcI68pqn8NRhJCUm4i+3x2Nspyis\/2iKUAvnvfIsbnl4gM249MAjf5RozaKuVbXAQxJZHg\/VoXXibbfdpms6N\/KuV1sdTfCwhUZ+WbbIUBQUIrLel92RSiP7mxEzv\/feF6hbd5mwMu0u3I1jQ4aiSs46VM6ujG7Du+HPmD8tkW+kLG8tWgzFkiXDxDkbexY6WTWj9QYZM\/irymBna5sR6526vtJT24wwiyz5jFrbaE0j05BUTPpQqBJCy1LHKp098OiNx8g7VZQ6NsewiQmfGfw8fop4ECcrR6IsMmDLaRMzY+cjMXkvEpCAlFhLFE8KzHEnLCZcirEmIv4XgygjOc7p9IgVZWL97+l5CtiAZ+\/Rk+gzfCLS3hyN6xvVtz7dzAzY6itRvlEKLRUbl4XYtBQBFIpsQ6DRCo1EgKLI\/5kCdWkoTwCpR6b53VzxcePN0kq1GyJq3Jdw1jvB86zif4JKgVIBQP73zVqMmvU5cpo8LOpS0l6zMmCrDyfQpKQUGwUykgRYSOJweChH08UBD+NSMrD0vuM4kmoJs2v0zL+7rOCJ9Ihaya2qte+tG\/zD3Xfxt3eNAqXUNnmnXu6yoGYTnG0\/HJMfu9WaZsS1R1paUdDDh\/qdRkF2LVzT5jCe+NdCIU1E7Lb0NGRlWjK1xSETCRlXUsIr8QhJ+synOlmx+OyLAHSfPFDU10rx7s447bV11g1Hqx8jCXjPb1uDMyv+hogB7yIwPMITr+L1fXrrJqpTCX3NyIBNM8XhdkO+uxkH9gdi94OzEXoxFNcmrcE11eshOX41Gv26QUxqFmLxQf0xWBgyELNe246OwZ8hPD4BC+pmWhI3Zcbh0TcfxMJP+thN8e4p7iAv6gv7t7kkGXw5Aa+n6GmvX2913ykFnrI4z8OB3ptsbYOD+wOx9q7XcO\/z29GgRoiwqmU0ySqVqJpAlIQ0jKuWgrv\/D8jKSEOTrHghedpn3IxvX5ooaK+V4t0Tk21WzAFfS8DrCVo66tObE\/7agIeA8+zwUVh7MBA9EgZbM2B3f2YsLpw6iuj7nsPS4W1FTlJ3CqcYyVlwM9auCUDd19\/D5RFvomVIUwEeKQNiqcfEIwPJsenAngQkzV8tLHIxJ4Ox9xnLrnlZgced96e2Zpi53R2DL7T3ZjqVWvN07z8SPfq\/gBe6trRmwH68ZSgyP52FxvcMwIHzle1mhOPJkJNbsbmbk\/nyeofy8zB4KtfMR3hOM9wc0hRrcrfjz9sOWSXP2Xr1EXr4kM08E4DiilYLKx2tk9rtDMa3I18tF7XtxD9LZ5Y2Ym3zZqbwJlB5M51swENM32PgyyjK2oCP0ucjPzQKwz74CvmrpuLh7t2Q+Nzz4OxujqSPXlpFkjxDpx9ByC9Nsfl7S2TQyD\/uwHPXdca3u07g1Z4b0ej0IdTJO4fASoE4G1kf585sFvVijgALIuOQIrLMkzkhEyNWROKNpBfE\/bI0GByc2EVkMDixbLo1wg1FuwlueJPTOUTdYVh7UT5VFVx1gtUKsUXjUNvJDq3yPbU\/d97BXlufAg+B4+4qW\/FGyqvW95kwcTKe+msfu+nj5RdnqTO5e1PcfWMdkUqeTN1y4l+q0+Gx42h403nsXtEA+3+rirbjvsbGSXeLCP5rcvNxufahK5nb6mLJlvdQpXJlDGjyFE4UbhGP6vh9HDZuSkBWAtm5Y\/HlthRcn1F8FuahUTi91OKFHVitBqLGrRDpOWQzsCwZ7F2X1zSU7kPde5Gtbae\/nmc1HBi1wpnFFNSPvSif5KpEhaPhqHXp97Jly1C\/fn2bJL5q6GE5vji5OXGfan8VGjz08u5mRpATV7VsECbAk7EjV6RQlFW35nGnkNP6J1w+HYJT8ywpEzOKSrI00+\/aeZH4z5al4hzJxUuXsPf4any\/7T3MPFec\/vzOTMStThYbpVTkBFQEGDm7tMzgsnSQLWbydf67xoPPi3WUeuyaQVet1b2odst9yEnthaiRi3DulxU2Ga09\/UWVI5lqhbKVgSWPhY9x9+3b1yZ0lB4g5Cg6jsLoqgceyaWLI7aS4zFJLbo2atQoMSz2YlfbqY66ngCoq32Wsra5mgGbB2AEPJQoq2fvAvz8dRjW\/3kKk+edRGbqtcLVJmnxciTUvh9xmUBCOtBuRzT+PTsZW2pswvIzP6PZR81R9YFYkCuPCOAem2zXRYci3JAaV2\/k5zg6uz+I0Wv1Gm\/NLxr5\/Gc4MqN3qet1B83Docn3onbC2yLYoD2TtCp9eP1jJACIWZKH6K768DEj0j3ZB05Ws+QscFSHU8LLsdtkFY1VNwp0yMEPObbFU089haCgIAEIHo+aCYPzl3J4ZXoOHbrjA4zUrnXr1qUyaJA\/Ih3q88Zj2rr7PHl5eZg+fTq2bt0KOtNOG1Z0mKpWrVqCUGcvFqBf+m9Y9XsuaoVWRlrCTRj9+U44UttIGg3tG47wXj9bQX9LYCMBINoMPXjTBjR4r71wYaNC+XMKdnWypBdJsKQEjE1JQ1qCY982ZvrIZz+0kR4cC63+uJUCVCxVrNdHL8HhWYnWjGzUz7nNK011m9EDD4W3p21iIoFGsjy7H0tHwVpkKUWH6TgbAnUmB57UklYscVTw0LmiBx54AD\/88AP2798vzl2FhYWJaKJDhgwR42RA0t8U0J\/OVMkgZfASeNR2u3fvFkHifRI8RHA6Ai2HxtUTc3oGA2pPPm2UfWLX+ZLM2EcPBonrZEWjEsfip9i9IGt1Y0uKxCvHEyZQVvjShi7r0LSyV5cCiReDh04zMXgsubq1i6Mon8R0clhjVrOI2amdHNCdQxLzUXg53YpsWPjpp5\/EQPhIBgGNjnATn1CsC1pn8REGOQcTSx498KgSy6ckj5wBO+X+GHy3bR8GTf0Ep+q0wjX1auAfA1qiaXSkNbCE1pTaU\/1kSx0BZPVqoONjJeDhNRG56KSnFzuMFksfPhRHB+LkA3K0Njm3\/jMxDFaXVDVLXp+4q7axqhaR9DYOp\/ZCwbF9pUhghqmaJQ917gg8dN9elE97VjMZADx4WSrRNfkkqp61jT6u586dA50bWrhwoTgTJgcp4aPvBFhH4CFA+vSaR86AvXnrNoz5ZBMeaBGB1o2qY9mvR5GTV4gPhz+E4GD3NkppgiiNolYCX0eSjY8w2GtHwKFCIJGLPcOAqwYDPemrd19PbdNr7033CXgqeMyKLuvNdCq1zyNLB7MzYKsTrgcErfrkhU3SR005T3UdhcANCAkX654LO76DbHpWM1CzSVo2VcuezY7ccni8ZkgebwJHeY7FZ8DgorD3AAAPAklEQVQzfuLrWJBmiYgjl34DnsYrL481tM\/jiND2EvmOHHkaJ060Ek3VjASsesTFJeLUqVuFQSEkZIM1ko5avywnWmtDlNZaeesW6TqLejNTlCUN9Z7lzXQSkkc+Uh3wf\/1wQ0xt9GgdiW7vbUbf1rWx64t3xTvSvUpVgnXdc7QI4iiR7\/PPz8Dly6+K+ARnzoyyJvulhWqfPmPx\/fc3oWrVZ7FuXTCqVy9JdUjPUfP46E2GWfftbYaW9SapWe\/jrf14PXhk8+alyqECNDuPnANnwP5j3yF07z8Kd\/91KN5NNCeUKxFFTeRL2a4TEzMQEdEGtWrVxu+\/b0PjxrFXsmBPsclF6kxGOE8yBRsr+AAe\/zZygM2bmcKTNHO2b2+mk3XNwxYbVQ1iVcte6FxnicH1Oc2hViJfZ3KRklVITcTr6phcaScfxZbdgfT68mam0Bt7Wd73ZjrZGAzUQIBEJDXKpRmEcyaxrl4i3\/IGj6v0MJMpytIxlN9Xz43HVbqo7cykk1lj4n50PQzMfqAziXXp2e4m8jV7\/Gb1ZxZTlLVjKEU99ZQ2okVbs+hk1rzJ\/ZQpeJxJrMuhad1J5OsJgsl9qmbuyEFzcCj1Matbj6Pnm8UUZe0YSvOxefNmVKtWzcabW33XCucYKqtTWmZlNdayM8wpOyhyO08m8nVmbK7UZeBUqhODWl1HIHdRCuqNWIhTX8ww5AdnFnho7GXpGMrZIGSJx5ukFdYxlMBCrhWcR0b2mSKJQUVOtegKw11NbWSTdMGxbCt4ivLPgPPzOJuAV6ZPPOLFEXMKw0UJeY0WTzuG2gPP2rVrK6ZjqNaJQhk86h6N0Ym82uuRabrweDZq90rGyS+mIyIhFTmTu6BKs\/Zub5JyXDoCzx4H3m1l7RhqDzwV1jFUL+eo3v2rHSSO3k+NGmo06KKe2saSh57tCDx0vywdQ5kWWmpbhXQMZQ8De\/sl8lqoIgPFzHfXA4+Zz\/J0XxXeMZQXnURoNZmvfE39ypWlb5mapUDPmsWnI9UsaZ5mJiP9X03gMfK+rtbxZjqVMlXLVjZ6YQaHVpYzT2fKVgluFDxG67k6oWa082amMOP9zOrDm+lkeJ9H7yUc5evhdH3NmzfHmTNnRLo\/MmvSdT6qyyBVnyPnMKUJ4US+cv5RyrzNQOfkVJx3h9yKqFC\/ZEmU0xHyhh\/nSKUj5vIhMC2pSla274a1x4s\/FOHIJUvYLM5JyjlKKQFv1\/fX4fNV31rzr9IBNHoOxQKfONES3ZRyHpGqvH79ekEPexmzzWJEX+xHj+\/K850Mg0dP9TECHnpRYhDKxEypC3mfh87TE6DUrNd81p2ZWwYP5yjl48JaGc041aI8dno2RWyhZzFYqF81SS6PUT07n3fwTzz\/aGcUXTgrMlYvDGiPLYWRNu\/FbWTgE3gonSTnV6XkwAQaSgt5++2322TqLk+G8LZn+zR4WAW6dOkSIiMjkZ6eDq3NUiPgUZPpMnPLTCanjHcEHo5FprVBqEooGTx0j77+dFSYAk6wJFPTKbJnA9XjZ1Fb+Whz1\/u7YHjYLwjYvU7wHEueBe\/PRId7u1qTHxNQWfLwepLyqNJR5Q8\/\/FCc+ZdTOXobA5fneHwOPKqDKH0tKWQQqz7eAB5m9qNHjwopIjOnI\/DQ2DmteqtWrYTEoPbchlVAZhgtb3LVW0L1lHCktnkCPGY6hmqFm6K9HZknzPawdwROnwCPTBxZneLoK3ovUdaSh44yyFmrHa2NVJWT3oVUN60cpSwdjXxtj\/75GxJ7PojfTwfihUZH0CC4QCTgXfD3OegQ17lMJI\/ZjqHUH1spuW9KgDx58mQkJCSIHLBaPopG6OVKHT2+c6VPs9qINY+9TVCZSEasberagpibvlJsMDCitmmta3788Ucb6UDqnpbKxVJIT\/LICYBlo4DW+NU1DwPm5qqn8XTMCey9\/y30n\/S+iCsmZ8vmNRmP3VNqm9mOoTJjMXho\/sjQQiDiZMIMMLm+PanlDrN6PXjo5dSM0qyqyXG\/qJ66zyOLcFV67du3T8R7cwY8\/GWjgHx0liguLg6ZmZmlwEPrEHkstOgm\/ypiYo49RkYIkqK0Vtu0aZONNUs2MLDLiRqqyai1rXe9kwJIVOZk18Jnh2uiQf16iLurs3XsToOHQgtxjC2dEENmOobK6zoGPn0UVPBoBUaUXYXM8krxCfCoXwd7+z3ufEX8bUsooMsUTZqUgGePXuS2kn7ddQyVQ0ZxPOuePXti8eLFpSQPSVcyepBlkv5evny5TRRSMw5S6tKpHJnKkKla9bgux\/FeNY\/WZQqWPPTGDsBjtmMoqWNUSLIzEF988UXMnj271Jrn+uuvx4ULF3D58mVQBjcKx8xtzZooXTqZ9SAX+jEEHhf69TfRoYCZTGGmY+iIESNELHJSeanwmk\/L2rZlyxZ88803oJyhOTk56NKli80msxkbv2bSyWymtIJHdbWRiWUGEcweuC\/05yjb9YEDBzB69GjQfg8xny+WPXv24LfffhNGBAIPSavKlS1eF2YVppNXB3qXYwvwwpl37z3hVS3vlcgLc2I40q9p45D2YShLA0XdpxQWBQUFbv3+Nuhb5OfniwRaRvqjOvTPleLPdu0K1ey38VrwqItM9fCbWZYTJo18+I6uyYELGVSd+3fGo48+iqDsIERERIjY2KRfHzt2zOnfiXcmIuJsBAoLCxF5NhKnap3CprBNSNmbgo67Otrtj55L\/1wp\/B7+bNeuUK+kjddnw1bBISceIpFsNniIsewFLmSma5zcGNt7b0fzo83xes7rVvCQ+w6ZThlMWr8vRUdjzs5obAoLFYA5+gClWUxAq1trYcSIX1Cz1UnkVs\/F29XfRuf9nREePhQ9z5wR4KT+LkVfQuOixqhXr55IOehK8WZd3ZX3Ka823kxHm6CHsbGxYLcX2afL7BgGRBCKFkrOmFTk2Gt077EXH8P+NfuBrDQgcy+QOAGTsiehoLAAZ5behidvb2gjiS5ERWFcdLToq\/mRI1jZoAGuqVQJhQUFKCgsxLWVKiGG0jKmB2BIyx9tJE1ibCIOBa9HzFcFyK4chMJOhTgaOg31ty9HyKEbsCfO4rrjbPHmSXf2XcqzvjfT0QoeeZOUNz75mjMuK0RouS929ZH3D4gg5BD5xx9\/iIUm7QfQM5988klrGsBKa17Drpi+wHwgpu9XyA66xzKHWXFAfBpyf6ohJMXhw4eFhHi4Zk1sXBgCxM1HTMxX+PXMLJv7LKm4vvx7YPD1OJAdhDZt8gTYQkOPitRSQUELsCBogUu8482T7tILGWhkL0epnJfHQDc2VbyZjh4xVbPxQT02wFRZuXKlsDSRUyad8aF6ZJzo16+fFTx\/\/etfsfyll1DjxAk0CQzEsjpnhAp2c251VJ1yBFMHh2qugUbWHon4gHj0zOvp1Brp5xo1cOupUzZrKjPWPN640HWWgd2t7wePQQqqQQq1LHUEHkrSSuChomUwIPB06tsXQdnZNmrW1rytGNN2DNrktRFGgPPnzosAfD+E\/SD6uv\/w\/eh9vrdbBgY2UFRk8NAXnzZG2eWJNk\/pPBIVto7K+0vsw8huPHQYkXwb2cWqa9euwuGUM8NxQmE6isEHELU8EnxC8pi1z6OGsVKND0R8IggdcKPDYFTat2+PuXPnWp0OieiDBw8WPmpff\/01brjhBkF0MlXT761ttqJTdCe8EvMKQjaGIL9tPppmN0Xfgr7I+yKvVH21vdHftAakf64UrUknVzX6Z6RoJe\/iXK167bXaUhtH7ek15VeVvbVV51MCza233oqff\/7Z6rLD6jqBh9LTT5kyBcnJyWJOCRykWWiBhzzbqcgeDdSOj734BHjM2ucxAh7Z3YeJy5F7nNkfKYgpEKZsTxXawKR\/rhStSaesdpSI2EgpKipdKyDASEs63q1dz1F7NUEyjV8+mkAfNLmQ9KHC0kiWPGrCYFbbtMBz9uxZh\/5wXg8ed\/Z5VDf0d955R6hkFF2UT4LSQTo52qhsbSPwqNY8RzvzxtjHnFp0TJv+uVKuBsnD4HF07IFoQzxA53169OiBr776ykby8OFDzobNUoYlG0t2+bSuTG+fA4+7+zx6BgNZOnEcgfLMseMKOPTaePOk642d7suSh36rXvbvv\/++WBPJgVbkcz+85qG2dFxk0KBB1vNOdFqXVPLOnTuL4+dy0BXVFcyb6eiRfR57Zu+UlBShB9PXSHbPKctjvUYYx4w63jzpZrxfWfXhzXT0yD6PHmHtZWCQAWUvsS\/72\/HXkXVxR8EXnXmeeiBOC9haRhD1nb150vXmx5vuezMdPbLP44j4jhL7Uhiml19+WTSXzdeyJON9E0f+cfLznX3eqlWrRLATWqNpRdBhqyQ5rpKHhLz56yu6ujeBQ28sPgUe+esvv5injiWw8eCRRx7BjBkzxN6PbN6k51JcN3ISpaAdbIigdvb84xxNiN7zyMghF9UKSQvja6+9ViR48oNHj\/Xdv+8z4FE3ON1\/df0emDnLKrGv0efRyLXM7nRdHqtf8pTMcYV2zzHbe1oPOrL3gWy+pnZqol4V2Hr1tZ7tzPPUM02yuZYsSezY6geP3ixbLHXsWaBf27aGz0geGrYnDr5pEUxWh3jx70gN0wKPM2qbM89TJQ7ThTcE+X20nF75njdPuhEG9rvn6FNJN0sCd2HmmkcraJ6eAUAFj159dd2ihtCy157aySZ1LRKqm7z26shBGakOeeYY9M5BnEanmfrzKWpotaXrjtqTE5LsiOR3z9Entg14tL64+l04V0PLIOHJxL7OPk8NpkFvp5rBXQXPFe8cGPTOgZaHjUHvHM229B6O2tO4kqWp9Lvn6PO1DXg4ThftFLuT+Vr\/sVd\/DU33HB+TPH73HMd8Wkptk8X11c\/innvDq2HNI4fU9bvnlOaVUmob+Rqxv5Jc3cw1j+dY1nt69nXweAslvZmOZe5h4C2T4ulxePOke\/rdzezfm+noB4+ZM60suFVrm4cedVV36zPgkX3I1Bnxq23O8ag3T7pzb1K+tb2ZjoYkT1kmMyrfqTLv6d486ea9pW1PFdo9xx5Ry9ptx1OTW5b9VkTw2KNvhXHP0SKA7MHs3\/8xBkFfB4\/fPUd\/nv8fV7LLDyKSADMAAAAASUVORK5CYII=","height":140,"width":232}}
%---
%[output:88f33f76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:02ee8dfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:6b77bd44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:243631ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:4d8b2271]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:60759f32]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:854eaa9f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:59d3a130]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"84×1 logical array","name":"c","rows":84,"type":"logical","value":[["0"],["0"],["0"],["0"],["0"],["0"],["0"],["0"],["0"],["0"],["0"],["0"],["0"],["0"],["0"]]}}
%---
%[output:31ea7669]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f852e64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c9e8752]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:293586e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:123f3ad8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:42132fa3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6353cdaa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96b54ad6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8801b32f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4cd076ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65fa3869]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:013471c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:08148c64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8e48a980]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3583ea30]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:553c6d0f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f2e2e2f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49423185]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0fd34b35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5de37e85]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:14e8daa9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b7116bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e7fb7e3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:656fb5bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99849bd9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6ead520c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e74eefe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0fd7a76d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37fad56a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:55c4a5f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:94ac4a83]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2134c1d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5462d411]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:02656171]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e6cfc9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ecbb07c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:65e26383]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e2a3163]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0afaa83c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2294c197]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b5938fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7cf75b5f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04b104f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34002a4f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:773da4a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b31cf5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:670bb4d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:399774f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7579b753]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:25e94361]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b07f81b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9711eebc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:789c8619]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:319c2495]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b755ac4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3396abf6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:59e5c958]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:160e7865]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1f821d55]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:43f484e8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0384aac9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:763fd066]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2d2f8cd3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:38746add]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3aa46ccc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a403478]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:034ad521]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7832eb2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:740d9a63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:635fccc1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8493fed0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:48cb69f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2cac285e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61e85ed8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37f0b234]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5df85d49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f05ee8d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e5e17fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35a75de3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2cc0df15]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7e0e5485]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45ac046a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:546934d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8041e268]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f745a68]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16a5f179]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87fb0582]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7df2602c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58d22566]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7cd57c45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:696de59c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35847a0a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8cec2449]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:828859a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47ef7e5e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5e26f736]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9fcbbc89]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:24b0785b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2c0c6292]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:371cb65b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:049c0398]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2be58b7d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f17db34]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03423f8a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:17947987]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4368b427]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1d061903]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8a0e0d1c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:181ad1ff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:15f1e605]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:61cf842e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8130ed1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:80f2e0d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7aa61f2a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f3579db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0f197060]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3738de88]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ed8dd4f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0982242e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c40ded9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:75d565b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:897947ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:44bacd6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8dfd0d99]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:416b9254]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:402e0fa7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:126097c0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:344eb300]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a49d798]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:86c24713]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4fca1f5a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:458428a9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7bcf3233]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1d374ea9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0c55ea45]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4a4c51da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c03b7c6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82fb2d74]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0a175418]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19d9720a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9f5ebeb3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:28b367ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:600c9891]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5a0b3f92]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9bdfd3f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a4087b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8c2ec768]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8fc2e19e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d6f9b4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:60001037]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:520f45a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1c87610f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1a31aab9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:191751aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6577cca2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:606aea90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0b03daeb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97c009fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5e7d2bc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f8483cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:59d715d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3b8bd2b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1c9303e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3125e871]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:62fd9624]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ee1cd2c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ffb7531]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:318d36f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:52da2a22]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3a5ceb78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04b307a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:48ed9cf6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6d83a0a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:30f3aac4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0d62ea9e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:71597e64]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:68921c4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f2e2fba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:415ca299]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2cc8c10d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5876c111]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2fb6da6a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:29b37f41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13f3f9d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:799394b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:219049db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:596658d4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22d84503]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05c96b02]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e914c57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32e40d39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:20845ffb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8ecfe50b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a028c52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3ea31299]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2392a08e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34088c29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:413775fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f994791]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:668f2837]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1322b157]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8647a23c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:960b51ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2a6ccdeb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:309d3545]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6edcbdfe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:14ab9601]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:67565d66]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1629796f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b8595ea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4ca27c84]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2bd8b0d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:214b4aad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9a64ade7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:159b4bfd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:42d2544a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f2bbc20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5ff0a4f9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3f50eee2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b50a1dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:26e76127]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ef9f7ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6b96a88f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:77e967f1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:09afd9ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:93487815]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95ed6b35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:095f36e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0f87b71d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3d60336b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7d2c88f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:97e5bf60]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:69feb6dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6585dd1f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5a3e6287]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45daeb50]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:669ba1db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6c2b64f8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7a008c00]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:829089e7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6aa36e28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3aaf4034]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95246d47]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:29ae58c2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7083da54]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2d1994a8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05d8ad38]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1a88a226]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ba54330]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9fe470db]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bb41207]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0cea6083]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16dda99f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12d92b20]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:45313845]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18a2c6e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:12b7761b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5da8725d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e2ec90f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3a4fe584]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e25e0b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:10ef21ce]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7548a969]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9d3d04e9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:947f7fff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:17290c99]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:537a0183]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e2856f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:87bd088d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4ce9bc13]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84525409]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39f91975]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8f91a60c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:70a5c384]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e1d5c24]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2ddc2551]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2e8f6f4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5a469311]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:620319d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:157801df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5d3706a3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5e51df3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:924a27ac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6f768d52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:443f3dcc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:588a78d3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:99d35ce7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8ca35004]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0fcadb06]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:432f908a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2a74e4e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2cdb7819]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:13e84aea]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:02ee6ae4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34f3f410]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:54b2f813]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4e709e51]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2debc682]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3de94c27]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:22413a52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:96023bd8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:10142ea3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:25860191]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:53db19dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8fd15a49]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:01f3075f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b1488d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4fdb4704]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:46f2845b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ec4d61f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2a14a5cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c98f522]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:617fc710]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2cbc7bd5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:123ae74d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:54bfb3c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c656e4e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7ba46003]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7bda5753]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5bfd5bff]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5c5c5e16]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a47b4e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:158222fd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:03bbad0a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9c3360ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:18b43b1e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:13d441cb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0fa89510]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:70630c73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a998272]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7718e7eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:02baf68d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:91fb0d2e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03a2474f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:607eb337]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:37fd6946]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8628d639]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:596e4d87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:41293588]
%   data: {"dataType":"warning","outputData":{"text":"Warning: no s.s. autocorrelation"}}
%---
%[output:46854725]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4a33397f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c30dfe3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:663bd1a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:852dae32]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Output argument \"ak_ss\" (and possibly others) not assigned a value in the execution with \"nanprewhite_AR\" function.\n\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('prewhite', 'C:\\github_trend\\Matlab\\prewhite.m', 66)\" style=\"font-weight:bold\">prewhite<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\Matlab\\prewhite.m',66,0)\">line 66<\/a>)\n[c.PW, dataARremoved, c.ss]= nanprewhite_AR(data,'alpha_ak',alpha_ak);\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('MK_tempAggr', 'C:\\github_trend\\Matlab\\MK_tempAggr.m', 132)\" style=\"font-weight:bold\">MK_tempAggr<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\Matlab\\MK_tempAggr.m',132,0)\">line 132<\/a>)\n    dataPW=prewhite(data, (obs{1}), resolution,'alpha_ak',alpha_ak);\n    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('seasonalKendall_main_D', 'C:\\github_trend\\aerosol_trend_analysis\\seasonalKendall_main_D.m', 149)\" style=\"font-weight:bold\">seasonalKendall_main_D<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\seasonalKendall_main_D.m',149,0)\">line 149<\/a>)\n    result_y=struct2table(MK_tempAggr(dataGTT_deseason(:,1),resolution));\n    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('all3_trend', 'C:\\github_trend\\aerosol_trend_analysis\\all3_trend.m', 8)\" style=\"font-weight:bold\">all3_trend<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\all3_trend.m',8,0)\">line 8<\/a>)\nTresult_MK=seasonalKendall_main_D(data,param, inst, station, resolution,varargin{:});\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('all_trend_STN', 'C:\\github_trend\\aerosol_trend_analysis\\all_trend_STN.m', 68)\" style=\"font-weight:bold\">all_trend_STN<\/a> (<a href=\"matlab: opentoline('C:\\github_trend\\aerosol_trend_analysis\\all_trend_STN.m',68,0)\">line 68<\/a>)\n        [Tresult_MK_25,Tresult_LMSlogi,Tresult_LMSlini]=all3_trend(data_trok,{names{i}}, inst, data_st.name, resolution,'end_year',max(data_trok.y), 'fig',1);\n        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"}}
%---
