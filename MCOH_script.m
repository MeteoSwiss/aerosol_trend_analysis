MCOH_st.name='MCOH'; 
MCOH_st.lat=7.3;
MCOH_st.lon=73.30;
MCOH_st.alt=16;
MCOH_st.env='Coast';
MCOH_st.footp='RB';
%%
% import neph:
Ysheet={'2004', '2005','2006', '2007', '2008', '2011', '2012', '2013','2014', '2015', '2020', '2021'};
for i=1:length(Ysheet)
    if i==1
datax = importfile_MCOH_neph("C:\github_trend\raw_data\MCOH\MCOH_Nephelometer_Hourly_STD_Daily_Monthly_Aug_2026.xlsx", strcat('Hourly_', Ysheet{i}), [2, Inf]);
    else
        datay = importfile_MCOH_neph("C:\github_trend\raw_data\MCOH\MCOH_Nephelometer_Hourly_STD_Daily_Monthly_Aug_2026.xlsx", strcat('Hourly_', Ysheet{i}), [2, Inf]);
datax=[datax;datay];
    end
end
%%

datax.Properties.DimensionNames{1}='Time';
names_neph=fieldnames(datax);
names_del=names_neph(contains(names_neph,["std","N_records","alpha","CF"])) %[output:1287d862]
for i=1:length(names_del)
    datax.(names_del{i})=[];
end
names_neph=fieldnames(datax);
%%
MCOH_rd_neph_h=timetable(datax.Time);
MCOH_rd_neph_h.BsB_S=datax.TSC_450_mean;
MCOH_rd_neph_h.BsG_S=datax.TSC_550_mean;
MCOH_rd_neph_h.BsR_S=datax.TSC_700_mean;
MCOH_rd_neph_h.BbsB_S=datax.BSC_450_mean;
MCOH_rd_neph_h.BbsG_S=datax.BSC_550_mean;
MCOH_rd_neph_h.BbsR_S=datax.BSC_700_mean;
MCOH_rd_neph_h.U_S=datax.RH_mean;


MCOH_rd_neph=retime(MCOH_rd_neph_h,'daily',@nanmean);
%%
% read Aethalometer
names_abs=fieldnames(x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4);
x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.Properties.DimensionNames{1}='Time';
names_del=names_abs(contains(names_abs,["STD","BC"])) %[output:0e3ab0e6]
for i=1:length(names_del)
    x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.(names_del{i})=[];
end
names_abs=fieldnames(x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4);
%%
MCOH_rd_abs_h=timetable(x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.Time);
MCOH_rd_abs_h.Ba1_A=x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.BABS_370;
MCOH_rd_abs_h.Ba2_A=x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.BABS_470;
MCOH_rd_abs_h.Ba3_A=x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.BABS_520;
MCOH_rd_abs_h.Ba4_A=x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.BABS_590;
MCOH_rd_abs_h.Ba5_A=x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.BABS_660;
MCOH_rd_abs_h.Ba6_A=x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.BABS_880;
MCOH_rd_abs_h.Ba7_A=x1_Final_Hourly_Daily_Monthly_Reference_444_487_No_MSCFS4.BABS_950;
MCOH_rd_abs=retime(MCOH_rd_abs_h,'daily',@nanmedian);
%%
% do the _rd
MCOH_rd=synchronize(MCOH_rd_neph,MCOH_rd_abs);
names_rd=fieldnames(MCOH_rd);
%%
datevec(MCOH_rd.Time(1)) % begin in 2004 %[output:3ae1e748]
datevec(MCOH_rd.Time(end)) % end in 2025 %[output:2a073352]
prctile(MCOH_rd.U_S,[5 50 95]) % in 50-80% compute RH but nothing is dry.  %[output:9572b09a]

plotFigControl(MCOH_rd,MCOH_st.name); %--> absorption cannot be used %[output:0fd6e28d] %[output:61cf64fb] %[output:282b8119] %[output:50ffdf32] %[output:51acd8aa] %[output:73b57dff] %[output:7596822a] %[output:96ba5ace]
%%
% BP detection for abs
names_abs=names_rd(startsWith(names_rd,'Ba'));
break_MCOH_abs_BP=change_point_analysis_Def(MCOH_rd,names_abs,0.05,'MCOH','Absorption'); %[output:69c2359d] %[output:7178928f] %[output:73fb85a8] %[output:86992b43] %[output:3da873bd] %[output:937c0a96]
T_MCOH_abs_BP=make_table_breakpoints_def(break_MCOH_abs_BP);
%%
lambdaSC=[450;550;700];
lambdaAE7=[370 470 520 590 660 880 950];
% % lambdaAE7=[370 470 520 590 660 880 950];
MCOH_cal=compute_exp_SSA(MCOH_rd,lambdaSC, lambdaAE7); %[output:8e8d30b0] %[output:57028429] %[output:11874e63] %[output:9b14cb52] %[output:00570387] %[output:9c491c5c] %[output:7c87fdab] %[output:53da11d8] %[output:60fab80a] %[output:133f7c9a] %[output:9d1ea025] %[output:076d2b23] %[output:2f137536] %[output:6908aff7] %[output:0d9b3f26] %[output:2cad31eb] %[output:2fb9451e] %[output:02dcab76] %[output:20bfb149] %[output:692db228] %[output:1801bd33] %[output:5f4d4999] %[output:249f9816] %[output:863778f2] %[output:9ee48031] %[output:57d3e313] %[output:15451798] %[output:979bdba9] %[output:7d96775e] %[output:3c0673da] %[output:837ca3f7] %[output:01d970f2] %[output:044fa5ba] %[output:1ed5920c] %[output:11b46dbe] %[output:674f9d41] %[output:383f4db8] %[output:563da939] %[output:5a5a0a8c] %[output:03fdd8b2] %[output:89247c3f] %[output:2093555c] %[output:92fd0be8] %[output:96788f96] %[output:50ad103c] %[output:04143249] %[output:92dd5f37] %[output:467cd496] %[output:6657d4be] %[output:10d57844] %[output:6da3ae6c] %[output:93851e33] %[output:220dff2d] %[output:41602578] %[output:6af10b35] %[output:2991f550] %[output:333b04cf] %[output:4983c5b7] %[output:28f1614a] %[output:6ddd829a] %[output:3f7884f7] %[output:41d4fa27] %[output:0d9fefb2] %[output:1fea820b] %[output:9fb6bbf0] %[output:3d998d9c] %[output:29d1b1b7] %[output:04d483b2] %[output:076161a5] %[output:77adc073] %[output:16996c90] %[output:5e1c8e91] %[output:854d5457] %[output:3840bef9] %[output:0e108e19] %[output:712df5c4] %[output:0fea583f] %[output:0df4f7b1] %[output:9f1faa57] %[output:91c96e39] %[output:17bd0fda] %[output:329b88b1] %[output:504b269e] %[output:51d1cdc4] %[output:1617207d] %[output:2c9fbdab] %[output:2fa1a655] %[output:90647205] %[output:5bedff53] %[output:26a0dc79] %[output:206dfb78] %[output:677d5102] %[output:2f7b30d9] %[output:0ae922b2] %[output:016035ee] %[output:817310d6] %[output:79004a0b] %[output:5c48ecde] %[output:9916eabd] %[output:8150bc17] %[output:1a440d44] %[output:94d4fefc] %[output:4804dee6] %[output:62fe1f9d] %[output:44b6d2f2] %[output:1ba547e1] %[output:0c4761c1] %[output:0be1ba44] %[output:111d9300] %[output:085d80fc] %[output:3d93dd73] %[output:78fb2c1d] %[output:2cedc9b5] %[output:0827a9fa] %[output:1b6875ae] %[output:88007c1c] %[output:931eadd5] %[output:61c5ad7b] %[output:31339437] %[output:3905f5f6] %[output:6bacc0e9] %[output:2cad2eaf] %[output:3bb615ef] %[output:5fe998d6] %[output:6c0c3242] %[output:97cc1070] %[output:92811090] %[output:558f4e50] %[output:8cfc6482] %[output:6a9ae863] %[output:763f99d1] %[output:913a087c] %[output:861546cd] %[output:4196beb6] %[output:9a3f74d5] %[output:5df0932f] %[output:2d3c4887] %[output:77b5f251] %[output:3e5144bd] %[output:42c6d1fc]
plotFigControl_cal(MCOH_cal, MCOH_st.name); %[output:8a71b0ce] %[output:7a046379] %[output:7b704838]
%%
% BP detection for abs exp

break_MCOH_expA_BP=change_point_analysis_Def(MCOH_cal,{'expA_fit'},0.05,'MCOH','Absorption exp'); %[output:01ec9e76] %[output:957844cb]
T_MCOH_rcpA_BP=make_table_breakpoints_def(break_MCOH_expA_BP); %[output:6c0abbfa] %[output:5e11db53] %[output:1f8ba7d8]
%%
MCOH_tr=outerjoin(MCOH_rd,MCOH_cal);
MCOH_tr.y=year(MCOH_tr.Time);
% begin at the beginning of a year: 2005 
%end: 2021 for scat and 2025 for abs
P=timerange('2005-01-01','2026-01-01');
MCOH_tr=MCOH_tr(P,:);

%do not compute trend on B, R and exp with red
MCOH_tr.BsB_S=[];
MCOH_tr.BsR_S=[];
MCOH_tr.BbsB_S=[];
MCOH_tr.BbsR_S=[];

MCOH_tr.expS_br3=[];
MCOH_tr.expS_gr3=[];

% compute abs trend only for 520 nm and fitted expA
MCOH_tr.Ba1_A=[];
MCOH_tr.Ba2_A=[];
MCOH_tr.Ba4_A=[];
MCOH_tr.Ba5_A=[];
MCOH_tr.Ba6_A=[];
MCOH_tr.Ba7_A=[];


MCOH_tr.expA_bg=[];
MCOH_tr.expA_br=[];
MCOH_tr.expA_gr=[];
%%
[MCOH_result_MK,MCOH_result_LMSlog,MCOH_result_LMSlin]=all_trend_STN(MCOH_tr,MCOH_st); %[output:8d83e427] %[output:2b313d91] %[output:557e7fe9] %[output:4b5e15b2] %[output:5773662c] %[output:1d59e938] %[output:9ade79aa] %[output:6bb7324c] %[output:4f4fe050] %[output:7f604946] %[output:6ea0556c] %[output:0eb45221] %[output:04afdb29] %[output:458de9dc] %[output:795f3b3a] %[output:77b0daf7] %[output:821b9aec] %[output:02fbbbdf] %[output:898fb1a0] %[output:892e4de9] %[output:5562f0fb] %[output:19210897] %[output:8a1cf6f2] %[output:9d95dc6f] %[output:7b903a60] %[output:8bb61f75] %[output:55511abb] %[output:45ba36de] %[output:2532efe0] %[output:6c02a654] %[output:139ec4a6] %[output:251f6ee6] %[output:2bfac671] %[output:3f70d40b] %[output:408c013c] %[output:6a88a66c] %[output:32712235] %[output:0ae0c948] %[output:8a2501cf] %[output:4d40bc11] %[output:5c3d6a53] %[output:5d69597c] %[output:8855113b] %[output:946682dd] %[output:8048bf8e] %[output:8ff0b9a2] %[output:6c1dd886] %[output:68537892] %[output:60606612] %[output:4098bad7] %[output:86521e56] %[output:206ff973] %[output:29894092] %[output:6e5d61ba] %[output:5bfa5163] %[output:53b6b0d5] %[output:3654cece] %[output:47973b33] %[output:674b4041] %[output:9c708581] %[output:97b7ffc7] %[output:27196540] %[output:7d0809ec] %[output:590930b9] %[output:2349716f] %[output:995f949d] %[output:0deb8e24] %[output:41231d7e] %[output:33fb96a5] %[output:5f29d1a0] %[output:88e9a44c] %[output:3b78f7b1] %[output:84ec0601] %[output:34b9deb0] %[output:7c6038e6] %[output:73af79b9] %[output:16c46dc7] %[output:82b652f5] %[output:879494cd] %[output:022d8ae5] %[output:020ce108] %[output:32dbbaae] %[output:95daa554] %[output:567b5e80] %[output:7834d797] %[output:05d81d15] %[output:8d3c3037] %[output:9b8f104d] %[output:279a5baf] %[output:58aae59b] %[output:35440d17] %[output:3fbf2f6c] %[output:1307ab25] %[output:5540db3a] %[output:7b3871fc] %[output:2f54b455] %[output:24eec615] %[output:6645ab9d] %[output:39069aa2] %[output:0954c4dc] %[output:483a055b] %[output:3516ade9] %[output:75f2c6bb] %[output:0c619943] %[output:111a2789] %[output:7c4e216c] %[output:2275369c] %[output:9a7a5115] %[output:06b549ad] %[output:05a760f6] %[output:4f5c3002] %[output:8674442c] %[output:19970e91] %[output:3994f613] %[output:35d371dd] %[output:7947507c] %[output:7ee5a05c] %[output:80435562] %[output:5f7b889d] %[output:8625f8f7] %[output:27e94020] %[output:47edea7b] %[output:083d465c] %[output:97d8ea97] %[output:2c8d35a7] %[output:14c7b5b8] %[output:586fed4b] %[output:7e23854f] %[output:90395e18] %[output:69efb019] %[output:14bc5711] %[output:130df383] %[output:8c8be2ca] %[output:5f3166b0] %[output:3592293e] %[output:2368f304] %[output:705a29e1] %[output:829fd44e] %[output:196b0863] %[output:96bda9e4] %[output:14114271] %[output:2c4909e5] %[output:7b28ec87] %[output:1486087a] %[output:34417615] %[output:73c3f070] %[output:8b51281f] %[output:0ecc1f28] %[output:1ed99c03] %[output:2b415f7a] %[output:95137503] %[output:23a2b67d] %[output:49b372d0] %[output:3c4fcc60] %[output:38fd57b3] %[output:3ba6c0cd] %[output:03023ca5] %[output:0505da1b] %[output:9159c2b9] %[output:5b1a0de6] %[output:54292d7b] %[output:0a6fc3df] %[output:01532a58] %[output:47469dac] %[output:03139ea1] %[output:4c1bd746] %[output:85df4659] %[output:74ebdf3b] %[output:7e1c809b] %[output:012edad7] %[output:6fa4dace] %[output:61b7ef84] %[output:44d0f9c9] %[output:911b8834] %[output:468e70d8] %[output:3ac308da] %[output:5fa7c808] %[output:84d765d2] %[output:3552ebe9] %[output:12b850f0] %[output:860be683] %[output:74e1c2d8] %[output:5cdba6c3] %[output:116f9b59] %[output:167f3238] %[output:3c8c1774] %[output:6a2785e4] %[output:805f4cd3] %[output:3e0c14c0] %[output:4aa06cfe] %[output:51e91622]

writetable(MCOH_result_MK,'MCOH_res_MK.txt'); %, 'delimiter',',' )
writetable(MCOH_result_LMSlog,'MCOH_res_LMSlog.txt'); 
writetable(MCOH_result_LMSlin,'MCOH_res_LMSlin.txt'); 
plot_10y_in_two(MCOH_result_MK, MCOH_st,'y'); %[output:0f4d311e] %[output:1bdfd500]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":35.2}
%---
%[output:1287d862]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"23×1 cell array","name":"names_del","rows":23,"type":"cell","value":[["'CF_450_mean'"],["'CF_550_mean'"],["'CF_700_mean'"],["'TSC_450_std'"],["'TSC_550_std'"],["'TSC_700_std'"],["'BSC_450_std'"],["'BSC_550_std'"],["'BSC_700_std'"],["'CF_450_std'"]]}}
%---
%[output:0e3ab0e6]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"21×1 cell array","name":"names_del","rows":21,"type":"cell","value":[["'BC_370'"],["'STD_BC_370'"],["'STD_BABS_370'"],["'BC_470'"],["'STD_BC_470'"],["'STD_BABS_470'"],["'BC_520'"],["'STD_BC_520'"],["'STD_BABS_520'"],["'BC_590'"]]}}
%---
%[output:3ae1e748]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"ans","rows":1,"type":"double","value":[["2004","10","1","0","0","0"]]}}
%---
%[output:2a073352]
%   data: {"dataType":"matrix","outputData":{"columns":6,"name":"ans","rows":1,"type":"double","value":[["2025","12","31","0","0","0"]]}}
%---
%[output:9572b09a]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":1,"type":"double","value":[["53.8710","70.3200","86.6378"]]}}
%---
%[output:0fd6e28d]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANsAAACECAYAAAAZW15iAAAAAXNSR0IArs4c6QAAHn9JREFUeF7tXQuYTlX3X8PQyHVo5H4bl5KIKXzli24kSbooUilCUUhySR9FF7mkcu+qni\/d8ReqR6FQkzyMMUSSGDU1cpdxnf\/z23372O+Zc9n7nDPvnHdm7+eZ5+E9e6+z99rrd9baa++9Vlxubm4u6aI5oDmQ7xyI02DLdx7rF2gOMA5osGlB0ByIEgc02KLEaP0azQENtgKQgWPHjtGIESNo0aJFxtuvvfZamjJlCpUpU8b47dSpUzRu3Dh65513jN+qVatGb775JjVo0MD4Dcvubdu20dy5c+mbb76hPXv2UHx8PDVp0oQ6d+5MXbp0oYoVK1qO9Pjx46wN3rFhwwY6fPgwlShRgi655BLq1KkTde3alcqWLRvRdt++fdS7d29KS0ujZs2a0euvv56HvjhGuzoFwPoCfaUGWwGw3wpsViDau3cv9e3bl4GAF3O9I0eO0IsvvsjAAnBaFbSZPHkytWrVKuLxjh07aOTIkbR27VpbLlSqVImee+45uuaaayguLo7V02DzJjQabN745quVFdhA8JVXXmHahJf169fTPffcQ0ePHrUEG8CFNvhzK3Xr1qU5c+ZQcnIyq\/rnn3\/SI488Qt9\/\/71bU6bZpk+fTm3atNFgc+WWfQUNNh\/M89rUDLbLL7+c1qxZQ3fffTc9+eSTzAREeeONN2j8+PF08cUXU3p6OvtN1GwrV66kBx54gGk0AGLYsGEMrImJiXTixAn6+uuvaezYsfTbb7+xtkOGDKGBAwfS6dOn6aWXXmIA4jSHDh1K1113HTNjQQ9mKbTh8uXLWZ3LLruMpk2bRklJSVqzeZx4DTaPjPPTzAy2Pn360FdffUXlypVj2ue8884jsQ7WR1gXiWCrXbs2\/ec\/\/6EPPviAgXPWrFl09dVX5+kW1mOvvfYadevWja688koGyt27d9N9991HMCMBXtBu1KhRnrYwUWFmLl68mD3jmlebkd5mX4PNG998tTKDja+bvv32W3r77bepefPmBiBOnjxJDz30ENNwMCe5ZqtQoQIBpNB4rVu3ZloKGk2mrFq1ipmnKADy8OHDDW1qbi9qT7wPfd2\/f792kMgw2lRHg80D0\/w2MYMNDoiff\/6ZaaDRo0fT\/fffzzQdhBteygcffJAefvhhZg5ysKEP0E74DVrr6aefppIlS0p17f3332egEbWVXcNffvmFevXqxcAPz+bzzz\/PtC73Rsq8UHsj\/+GSBpuMtARcxww2aLODBw8yQEGgn332WaapYBoCcLfddhsDoB3Y+vXrx7STbBHB9u677zLNaFdEk1GDTZbD1vU02Pzxz1NrM9gg8HBMQFOVL1+erY3gzl+2bBn7d8OGDQ0tFrRmgxMEe2l25aeffjLercHmabqNRhps\/vjnqbUV2OCS52sweA2hfbDhjA1sAFB0aOA3gBOm3I8\/\/qi8ZuMmqsyaTayr12yepluDzR\/b\/LW2AluLFi0M72Lp0qWZMwQuf6zjihcvbqyRuGarVasW2yb46KOPmHNj6tSpdMMNN+Tp2KZNm5iWxJbAVVddxZwoO3fuZPSwHjPvoYkEzHtx2hvpb961ZvPHP0+tRbCJ+2Z8X40T5Y4PuOC5Q8Jpn23UqFEMVNB62GdLTU1lx722b9\/OSGLth\/8DvBMmTDC2E3BKZNCgQRF7dPByvvDCC8bpEr3P5mmqIxppsPnnoTIFO7CJLnkQ5Z5J0Ukhgk3lBAk2o7GfhvOSKPoEifK0+W6gweabheoE7MAmbjaDKryUOCJ16NAhdvIDYPRyNtLqfCPo433wYn733Xe2g9BnI9Xn166FBltwvJSmZAc2mIuPPvoo80LWq1ePOUdq1qwZcZrE7tR\/RkYGAyeOV\/31119sHYdTIfAgwnx0OvW\/evVqdhJl3bp1rC0\/9X\/99ddb3hjQJ0ikp1qbkd5YpVtpDvjjgNZs\/vinW2sOSHNAg02aVbqi5oA\/Dmiw+eOfbq05IM0BDTZpVumKmgP+OKDB5o9\/urXmgDQHNNikWaUrag7444AGmz\/+hab1hM9\/oV37cqhWxQQa3qFuaPqlO3KWAxpshUAa5q3NogHzthgjGd6hjgZcCOc15sGWmZlJH3\/8Md16661Uo0YNWxbjHCFiIuKUOw+oozofYaGBfot9GfThTwTA8dL9sio0vfuFrsML63gKcn5cmeajQsyDDef6evToQW43jnNycuj333+nqlWrUkJCgieWhYUGOi\/2ZX76gQjNBqABcG4lrOMpyPlx45mf5xpsCtwLs3BizbZq+wFqU7+CtAkZ5vEoTItRNYjxeHmvbBsNNllOmbRJQX99gxCssNAwa+qC5q2CSChV1WBTYJcWzrzMCoInGmwKQliQVfWazXoNCscR\/pwKHCTIJ4AYlF61SRA0uMMnmn2BM83JoZYfMq01mwJXg\/iKB0HDTRMAZAhFjrAIulhzAElGJk6cGFXAabApSGMQQAmChhvYuLaHMFWvXl1hhEWjKj5CyHXg5sEOmhsabAocDQIoQdCQBVu0hUmBlQVaVXbpEXQnNdgUOBoEUIKgocGmMGkWVTXYPPJPlnFBCHlYaGiweRSW\/zWTlRl\/b8nb2lazIfjMmDFjqGXLlnTHHXfkaYlUsgjaOX\/+fPasffv2bFGO4KFigbdq6dKlNGPGDNq6dSshKR8SRSAQzTnnnGNURapaBJ4BTWTCRFQnxLdHthUELbUrsowLC1CC6IcGmz8YyMqMv7dIgg0pgZCiCEBChhUz2HjMwT\/++INlr0T6V+RzRsi1mTNnstj0KAAQIvYij9gtt9zCUsUi0yWiRvXv358lkuDn4BCmbcCAASwV7Z133kmbN29mNJFT7KmnnorINS0OQ5ZxQQh5WGiEDWx2mVT5PKmuHcX8AmaRDSIjjqzM5CvYEFse+bgAsF9\/\/ZW9yww2AAgZKBGLXkyil5WVxcADoAFc0Fo8zHWHDh1YiDYAC+0RNg2ZWt566y2WiwzgBtB4jEKexB19AZjRB6vQ2uifLOPCApQg+hFWsKFfSClVqlQpQ05l50cUbA42yJP5Q49IzqAJ2bMLz+cGEi99cqMp8zzCjOSphCD0PXv2ZKCAsIsD5knVEYdeTEmLlyF8Nr5iPN4h6E2aNIn9n0fiRT0cCEY4bQAIwUcxeJiL5uyZPI4imGqXf0yWcUEIeVhoxBLYuNazAqKdgDqBzemZjMCrfKBl6cnWiwDbJ598wnKA4RQ9gnUic4r564Jcy\/jdDEK8kIfPBuBgDkIjYf1l\/gohacSIESOoWLFi7Eu4cOFCpi0BygYNGhh9hxa0o8ErabDlPUEiyxNZIXGr5wQo\/gzBZXkOOSszUTQ1ZcCGj7hTXjmnPkebP7wvtg4SuwHzjvLQ2FbrJ+T8QjRdAOrMmTMMUKKTg08AgI0c0vhDAndkbEFMerHAbFiyZAkzOeFcMRdZxoVFKwXRj1jSbJgfnsQDFoqVXKHOY489Znxs3cxIjF8l+aNXmXH7yKg+9ww2q0WvyKSbbrqJgc3OfBBt71dffdXWDocpCg+lWetpzWZ\/P0\/2A8R56DekgpuDRHRqmIFlJbBODhLUV820qsEmLHSDABvMzS5duth+bKBR4MQ5\/\/zzIxbwKl+nsNDgmk0cz7H\/m0CnsndSfFIdyqh2DVsDy3gBgwip4GRGcuCkpKTkycdtBxonzcb9CqLjTsw9AN64eSz5xwg+BitTFI48r7fFneQpUM3GByGakVaazcqMtPMwyZqR9957LxMwuwJPa3Z2NjNTxf09FbCFhQb6HNGX9E8p94NhxlDSKrejoUt\/lQIbYpd4Cakg8s3NCYI5XLRoUYR1gt9mz55tkBGB5wQ2Dizs52J5gnfD2YbtIu7IAyDfe+89W48ll9MpU6YwYJoLbkEgaWTQRRlsYXWQYA156aWX2vIHk4J9wSpVqni+ThIWGhik2JdjcwfSsVX\/PQu2+Po0NPW0FNjMmk02pIIK2NyWAlxbccA5gc0MbNxwePzxx1niRu5c44DEus5Kc4l+ByuZCY1m067\/go9jws1IHlPl5HfvUfb0+wz535HSn\/rO+VwKbGjkJaSCCtjc9sasAGTlCcc7OZAAIoBJXC\/KmM2gobqmDUrDKWs2vqm9YMGCiNMiTpvaN954o3FaxGlTG8F4xNMielPbfprNXs39HzxFxzavoFKN29HWWh2kgiAFJUROZiQXbL7GshJ0syZz80aaTVIRQHxMTsCLGbBhMMhYCZUPJvfp04e59e2Oa8H5wROoYzvghx9+sDyuhfOTgwcPZpk2b7\/9dpYHWh\/XkgebWDPawuTmjTQLPu+f2GfR4eH3uBbam03LguQPf7eyZuMNzQeRce5x6NChlJycHCEhVgeRYSIAUG4Hke+66y52GLlcuXK2UicrWEHsb4WFhtmMNIc0kOVJUJqtIOlgrLCAxH03Dv7u3bs7rtlkzc6gxqfvsylwUoNNgVlRqmrlDHHbyyuoj5EGm4JQaLApMCuKVc37bFZ5x0NtRkaRV75eJfuVCgtQguiHNiN9iUz4vJH+hhO91hpsBX8QOXqzHcybZGUmmLedpaLNSAWOBqGVgqChNZvCpFlU1WDzyD9ZxgUh5GGhocHmUVj+10xWZvy9JW9rrdkUOKrBpsCsEFfVYPM4ObKMCwtQguiH1mwehUVrtugwLgghDwsNDbboyIy\/t2gzslAmQyzKJ0i8AELWGvJC26mNXrMpcFRrNmtmqZ6NVGC5EUYBITR4cdu0dqOvwebGIZvnsowLC1CC6EfYzEiZU\/+q5xCdrs6Y77+pio6szKjSdauvNZsbh4TnQQAlCBqxBDa3W9x27Je5cOoUm8ZpWjXYFIRerCrLuCCEPCw0YhFsKqHszBdErUQDdRDWAKHsVdNiycqMR5G0baY1mwJHNdic12x4ahURWTWUXX6DIb\/p24mUJ7A5LYjNkY10Yo1I1gcB2LBqNtygtiqqoezcTEiF76Nl1ZgCG6JU4YZ25cqVqXHjxhEDKl++PHXr1o0lwtCJNfLOdRjBhpAKPAxeYrcxyrIcdCg7O7Bxx4jYQauoXNxziUxJZk2LtjEFNh5ha9SoUdSpUyfbydGJNcIPtsMr3ooIFpR4+1hSBZybE0Q1lJ0MGMxhzXkoBTEsOd4L4AWR7EP5C2TRwJMZiZj+WJgifByy0Dh5lHRijXCbkYjKBcDxUrZdL0oa8KaSbLmBzc0sNLvyZRwkZrABoOLaEAOwi0UiA2YlBkhW9gQ2BMBEXH6ercbqXU5JMXRiDX\/h8IJcs5k1G4AGwKkUN7CphrJDyik3gJrBxgGKJY5dqHo+ppgBGwdReno6C9qDgWVkZNBFF11Effv2pY4dO7LQzZwZOrFGuDUbeieGwVM1IdFeZlNbJZQd5xiPmmzOEShG5zJvlovrOruTJjEDNq6VFi9ezEI3I+w3zuZ9+OGHtHz5chYfEn8nT57UiTUs1EMYHSQqWsyqrupxLbdQduI7zPFF+DOZEylWpmVMOUiQJRTpfRAPXQyoChc\/Nhnfeecdlt4JXspoZrEp6ok1RAFFbE7ZxBp+gRaW9tBoKGLiTgAV4RXhyBPz\/nGwhyaxhhcmbtq0iSVKhLZDsgM7sOnEGv4SfGBunJJ8pKWlsbTKMl9\/L\/McxjZWzhC7BBscbKFJrOHE0BMnTrDHJUuWjKjG3a\/Y30BE2mhmHi3KiTVgxh\/bOpVO\/51Jxc+tQWn7W7PgtkUJbNz7iI8932ezSx3FwWYnM1FPrGEHNqTtRaRi5LhGmh6xrFu3jnr16kWjR49m6lzn1M7LxfxYs8Vlf0pH1wspo\/a3pt5PpRY5sMlq3JhxkCCBBnbt4Z59+eWX2SkSFCSbHzNmDG3YsIHlxapTpw7xTW2dWOOsGOQH2E5veZKO7\/7IeMn632rQAy\/s0WCzQV\/MgA39X7ZsGVt41qtXj63PihcvzryRqampzHTs2rUrxcXFseNaOrFG5IznB9jMmm3T8duo1\/CPNdgKA9gAoo0bNzLNhtMkKMg+88gjj1DTpk0Z0HjRiTXyH2x8zXZybyqVOK8VW7P16NFDg60wgE3WNo5GPVmTIAiNEhYa4KtTX2R5Eo35CeM7Coo\/no5rhYmBsowLC1CC6IcGmz8J5DLTecgUeunhm\/0RU2itwabArCCAEgQNDTaFSbOoysF2pM3jNPSuDjS8Q11\/BCVba7BJMspNwGXJaLDJcir\/6olgu71jW5re\/cL8e5lAWYNNYAYSue\/al0O1KiZYfu3wfPOuvdSwWgV6olMDTxOkweaJbYE24mBr17wplbrtuQiwucmAn44UGbC5AQXPJ3y+0+Dl8A51IgA3b20WDZi3xfa57CQECbaxC7bQ\/pPxVK9ymYi+yq5jZfvsVk\/1ILIbPaec2mhrdzLEjS5\/zvkzueHvVK\/rIGrW73n2KKg5tutHzIMNIHnx1f9Sq1at6fKmyZYayY2J5udgVvfLqkR88QA01OPF\/Fx2ooMCm9OYCgps4EEQt6I52HB7RDxYDPr8FkCtWrUsQx7IzIMItjad7zQuywY1x4USbFYgMWskDBxMLLtsAlU9nUW\/F69Ch68d7ggktLkiuQItGnD2FrobDZlJRp2gwOYkGB8uXUHDB0TvbKTTfTa3i6VWfHMCG+r7\/Zg4gc1JTmTnuFCCjQtc34NzbYGEgX\/22lRq+PkQgwcZLQZR55FTjf+LQMKPAOSc8veSCNw5E5+na78fabQBYLn5oTIJQYENGh2XPq0+IF2feo\/S5o6K2qa2DNhU4kbKgs18qVR2HjjYXul4mFo2vtHQbGY58TrHhRZs+BL1PTTXEQTmOBuLSnegnO4zDZPTzGQQQ51tHV5kGhAa9M\/pvajz0c+N93iJ1RGkZrMDP\/o66JUFVGbVC6EAm\/kCpxWQUAd3JHk4Azew8dvYXm81cLDNHphLrVveTGXbvMjmNYh4LE6Aj+k1mxVIrECwbtydlLjxn8uFKOvOaUYrOs0zTEkzk3k9\/mWD5oNmDAvYnMCPvu5Y8BJtXftF1MEWVNxIJ7DxZykpKb7XbABbjWad6MJrphlgO7lvARUrE09njpyiEhVvVg5+VGjBBpBw5mCQYNDv1IFSnnzPGDMXzNuTUw0m5mw8wNZt3Aw00zHTggC3W9ydrmj0q0EjdX9bpvlUN0SDMCM7T19P93x7X0R\/6PQlVG3scmYyH\/toJA3dVlUabLgLd+bvPVTs3OpUqtFgWWvMqBd03Eg3b6QYK5J3wtzGLmYk6ouabWeZK2hf7XFsHrcsG0\/n\/\/26MS6AjWs9ZaZYNIhpzXbws4fp9IlPI4b14fZWdOCi5wwQQDD7\/Pkw\/Ttlj1EvJ+1AxFfLig4qc1rQijfmTKCEZhUMGng2MrNfxLpOZkL8gg0fj6WvTaVnE6fTuZefZ7zyxM9H6JzqA2nH1i20a\/VCabDhao54F65Uo0HKgHNzgqjGjbTSbNwLiQHjClfFihWNsavEjLQC24ayTzArZ8uXA+n8o4sNuufUvI1KN58oM61SdWIabBAS8R4XRvxtZjLNyp3IPIlcMF+44HUqmVzGYMipP3Jod1ZbQwNa0eG0ntg3jvrlTqY7mq2LYOjOX8vQBxtS8ng23bjuF2xcq7Vt+wfFn58Q8brUfW1p5TdEV5xaSIO+SZTSbOaxexEwN7C5haUzx420MyPtAKcSM9IMthMJybQk8RVqUz+Ryq3tSf+q8bMGm5UQW4EEX\/j301KYAwSnQWD+QauJYBO1FswHO7D9vWYvTcu+lf6dkhkxCWgPwB75IovMns38BJudVuPv\/DjrCqq8ZyedW34P9ZsWJwU2s2bDlxyAUyluYFONG5mZmcli2Vjts3ETUDQlVWJGmsF2cdwR+uPc3swDfcnhZyLkxMuHx4lvodFsuCO3evVqQs4thF6oVKkSi6OBKFGlS5e2HIMVSLBuA+DWxd9LJZLqMMdGmfZV8mgBaEB80WA+2IENdLZ\/H0\/1GudEmJC8MwBc6oFGtLnpm9JrNz+aja8drbQa+gRtW+1UFqXnlpEGG9phzcbvwgW9ZuPgUIkb6eaNtIsnKRMz0gy2ZmVzKO7E9WxKc0t+FiEnhRZsuIQ6YMAAatWqFYttsnnzZpo7dy5deeWVESHzRNTZgQR1Vm+tTTsyEgiOEXGtJQLlo21dmBk4sPKMPOYorwfAwTtlNtnEfqw9cQ3B7pdxlvgBG5wfTbY\/GbFWE\/sB8KOkHU5QApuKFrOqq3pcyy1upBvYnLKS8v7ZxYw0gy2lPrH1e2ZOpQjnCOoVSrAhFiWABm2GLyAy4KCsXLmS3f7GbzfccEOeeXYCGzQchM9sPopEdu2twDTTZSUyqE7tI75kDiYc92o5EfIDNnhNzV9fq3et205RBZsvxgXQWCVmpBXYip1bg\/XizN+ZEb0plGDDVwjm4qxZs+jqq682BowgQoh\/CM8TonmZQ+c5gS2AOVQmkZ3Ujxr9a4RjOz9gO7puEh3fM921X0UNbCoxI63AZsfQQgk2JOqYNm1anoQITsk5wKCwgU1mcnyBzcL7qjXbPxww77M53QwQ99lgRhYpsGHBu2bNGpYZJykpKWLseLZkyRIW0rxu3cgbtUUNbIdXDWGb+G6lqGk2N36Ynxd5sIEB5s1KMMlpjyZsYEMIuZI1b3Wce0Qb+\/LLL+mCCy6g2rVrK8lJYvZLlJT7nWsbDTZnFmmwxTjYuIC7IsFHhTE9cqlzS3cCGmwabLYccNr0dDIjzRuy5hfAy2T2MFn9JraLr9SazhzLzNOO17Frvz9pMGXHtXJHgo8aSbmplJh99mqQHSkNtmDA5mWD3+nNodjU9uogwcDEDVn8H5uzKAhWig1a\/tzqNyShQDIKXngbka7IPCuaIl0fOJJuKo5H7L\/477XpmTr8uANHuRn5+phW1LLpWRkQeSjKgvTkuFQMBdi8uv6DYkJho8OFadCgQeyQgC6RHNizZw8NGzZM6jhbkLwLBdj4pnbVqlUjTou4bWoHyYjCRAtnCyFMyL2gizUH8BGaOHEi1ajxz4Z2NEoowIaBLl26lAYPHsxyBiBX9\/bt212Pa0WDQbH6DgAOf7pYcwAgiybQ0IvQgM3qIDLywOEwcrly5bTMaA7EPAdCA7aY56QegOZALDhI9CxpDhQFDmjNVhRmWY8xFBzQYFOYhqNHjxKSniOEGtzHuBJ08803U9++fSPOdHq5CItu4E4fDmTPmDEjIsYGnnml6TS8ghwP+rVr1y7mEfziiy9YN9u3b8+8qIh2XBiLBpvkrGZnZxP2rXC6HFk9W7RoQenp6cxjCuFAFtbq1aszaqoXYQEknJccOXIk85BZnRFVpek2rIIeD\/YCEfYA40WqaPAAPCxWrBjNmTOHkpOT3YYQc8812CSnDDERR40axTRP27ZtjVabNm2i3r17U8+ePWngwIF04MABpYuwcM\/PnDmT5STHIWWrqyFeL9c6Da0gx3P48GGWk\/3QoUMMYJUrV2ZdNfNSTBctOU2hrqbBJjE9Z86cYUKBGCmzZ8+OMPFgio0YMYJ9kZFUIi0tTfoiLA9UgzZdu3YlvGfnzp15NFvQJ2wKejzr16+nXr16MZ6KHy4e7gC8HD9+vG3sGYkpC2UVDTaf04KvMzRaYmIiA9vChQulL8ICbFizdOnShVq2bMn+bXXVyM\/ZUdXhRWs8WPvCXMapoaJSNNh8zjSOlAFsWNgjtIPXi7Doht3tBz80VYcXrfFkZGSwo3mffPIJ4WNy8OBBateuHTtF1LhxY9Vux0R9DTYf07Rt2zZ68MEHmSeSrz2crgu5BSt1ApuXy7WqQ4vGeI4fP87Wvj\/++CNzikCz4XjeX3\/9xW7qI+4MHCTNm59N16U6jrDW12DzODO7d+9mwYhgCsLB0bBhQ0YpVsEWrfGIYej69+\/PeBgfH894l5WVxTyU2FKZNGmSEWXN4xSFrpkGm4cp2bJlC9sGOH36NFufXXjh2QToXi\/CupmRdprN6XKt7NCiOR54XMeNG8fiyiBFVJMmTSK6CX7CM2sVc0Z2PGGtp8GmMDMwe7Dfha9xnTp1aMqUKVSzZs0ICn6cGXZA9UPTaXgFNZ7JkyfT\/Pnz80RTQ1\/dTG2F6QpdVQ02hSkRN5afeeaZPJHAQMqPm94ObH5oOg2voMazePFiGjJkCDuN07p1a6OLAD88sitWrCiUnkoNNkmwcecBPGUIGAtXv1XxcxHWDmx+aNoNryDHg7UZEmPUr18\/4rIwNvjxOxIdjh49Ok9QXsmpCm01DTaJqYEHDQCbN28eIckeTEhzgTl50003MQHxehHWab3nlabV8MIwHmg3eCUvvfTSCG8kNrRFh5PE9MRMFQ02ianCOcI+ffqws5B2Rcx06fUirBPYvNK06m9YxrNx40a2ZQJzNiEhgTp27EgPPfSQPogsIZO6iuaA5oADB7Rm0+KhORAlDmiwRYnR+jWaAxpsWgY0B6LEgf8H9badgoBLrmYAAAAASUVORK5CYII=","height":58,"width":97}}
%---
%[output:61cf64fb]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANsAAACECAYAAAAZW15iAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl0FkW2vsEQgmwJ+76IoDJoANlEBtlFEYFBxBEGUBiQfRMJPJFdQAFFQJYXBHRGcB8eIuoBlRFQRA5GliCiZoQIypIIQYgm4b2veNXp7r+36u5\/Sf6qcziH\/F3rrfrqLnXrVszVq1evkkySApICQadAjARb0GksG5AUYBSQYJMLQVIgRBSQYAsRoWUzkgISbD6tgcuXL1NycjJt2bJFqbFz5860ZMkSKl26tPJbbm4uzZkzh1555RXlt+rVq9O6deuoQYMGym9QpY8dO0YbNmygTz\/9lDIyMig2NpYaN25MPXr0oJ49e1L58uUNe5+Tk8PKoI2vvvqKLl68SMWLF6cmTZpQ9+7dqXfv3lSmTBlN2fPnz9OQIUMoNTWVkpKSaO3atQH1q8dolscnchbJaiTYfJpWI7AZgejs2bM0bNgwBgKe9Pmys7PpueeeY2ABOI0SyixevJhatWql+fz999\/T1KlTad++faYjq1ChAs2fP586depEMTExLJ8Em08LwaIaCTafaGwENlS9bNkyxk14OnDgAA0cOJAuXbpkCDaAC2Xwzy7Vq1eP1qxZQ\/Xr12dZf\/nlFxo7dix98cUXdkUZZ1uxYgW1bdtWgs2WWv5kkGDzh46kB1ubNm1oz5499Le\/\/Y2mT5\/ORECkl156iebOnUu33norHTx4kP2m5mw7d+6kv\/\/974yjARCTJ09mYE1MTKTff\/+d\/v3vf9PMmTPpp59+YmUnTJhAo0ePpry8PFq6dCkDEK9z0qRJ1KVLFybGoj6IpeCGH3\/8McvTokULWr58OVWqVElyNp\/WgVU1Emw+EVkPtqFDh9JHH31EZcuWZdynYsWKGkBCP4JepAZbnTp16KmnnqLXX3+dgXPVqlXUsWPHgB5CH0tJSaEHH3yQ2rVrx0B54sQJeuSRRwhiJMCLum+66aaAshBRIWZu3bqVfeOcV4qRPi0Ei2ok2HyisR5sXG\/67LPP6OWXX6amTZsqgPjjjz9o5MiRjMNBnOScLSEhgQBScLzWrVszLgWO5iTt2rWLiadIAPKUKVMUbqovr+aeaA99zczMlAYSJ4T2kEeCzQPx1EX1YIMB4rvvvmMc6Mknn6RHH32UcTosblgpR4wYQWPGjGHiIAcb6gN3wm\/gWrNnz6a4uDhHPXzttdcYaNTcyqzgDz\/8QIMHD2bgh2VzwYIFjOtya6STBqU10gmVtHkk2MRpZlhCDzZws19\/\/ZUBCgv66aefZpwKoiEA98ADDzAAmoFt+PDhjDs5TWqwvfrqq4wzmiW1yCjB5pTC3vNJsHmnIatBDzYseBgmwKnKlSvHdCOY87dv387+37BhQ4WL+c3ZYATBWZpZ+vbbb5W2Jdh8WgAOqpFgc0AkJ1mMwAaTPNfBYDUE98GBMw6wAUC1QQO\/AZwQ5Y4ePSqss3ER1YnOps4rdTYns+tPHgk2f+hoyNmaNWumWBdLlSrFjCEw+UOPu+666xQdiXO22rVrs2OCN998kxk3nn\/+ebr33nsDenjo0CHGJXEk0KFDB2ZESU9PZ\/VBH9Ofoakr0J\/FSWukTwvAQTUSbA6I5CSLmrOpz834uRqvgxs+YILnBgmrc7Zp06YxUIHr4Zxt7969zN3r+PHjrErofvgb4F24cKFynAAvkXHjxmnO6GDlfOaZZxTvEnnO5mRm\/csjweYTLc3ApjbJoylumVQbKdRgE\/EgwWE0ztPgL4kkPUh8mswgVSPB5hNhzcCmPmxGU7BSwkXqwoULzPMDYHTjG2nk34j60R6smJ9\/\/rnpyKRvpE+TLliNBJsgwcyym4EN4uLEiROZFfKGG25gxpFatWppdDwzr\/\/Dhw8zcMK96ty5c0yPg1cILIgQH628\/nfv3s08Ufbv38\/Kcq\/\/bt26Gd4YkB4kPi0Ei2ok2IJPY9mCpACjgASbXAiSAiGigARbiAgtm5EUkGCTa0BSIEQUkGALEaFlM5ICEmxyDUgKhIgCEmwhIrRsRlJAgi1K18DCD36gH89fodrl42nK3fWilAqhHbYEW2jpHRGtbdx3mkZtTFP6MuXuuhJwIZiZkILt5MmT9NZbb1GfPn2oZs2alsODjyDiHcKDnQfL8UKPSK4v1H0D0AA4nv7aoiqt+OstpuQNdf9E59nv\/om27zR\/SMEGf72HH36Y7G4So\/NXrlyhU6dOUbVq1Sg+Pt7peEzzRXJ9oe6bnrMBaACcWQp1\/0Qn2+\/+ibbvNL8Em1NK6fL5OcF+1uV0o4LOtut4FrW9McFWhAxH\/0Smxe\/+ibQtkleCTYRaqrx+TrCfdTkFm8iwo61\/IrQRySvBJkItCbYiL9K7XA6OikmwOSJTYCY\/d3vRumBowj+zBIMB3hRAHEo\/9N2iWB8MdHZGOpdLw7RY1IANOsqRH89Sw+oJ9F\/dC16LcUtQUYBYtSNSF0CGkOQIjyCTewrgQZJnn302pICLCrAF41xJBCB2S0KkLm7RxUKpUaOGXdXyuwEFsFHhXQQnVnE\/CRgVYBM9V3JCYBGA2NUnUpfI8Yldu9H6PVw0jAqwQYRc+EG6srbsDnGdLEIRgNjVJ1JXuBaK3RgK0\/dw0TAqwCY5W2GCQvD7KsGmo7HIbm\/nVCs5W\/AXcGFqQYLNJdj0QDJyqpWczR0UzF5T5bWpDQw8OtdDDz1E\/fr1c9fg\/5dCsNnVq1cb1oHXgbzWz8E2as6LNKl\/N099FSlcqMVIvZURA7+zfgJtGdVUQwNRX0AnBBThvHb1idQVyl2Zgw39x7NSJUuWVIai74ffYEP9CECrDtfHHwTBy0BeAMf73mPCElo6ppfd9Pj2PWLB5uRcTM+xOFX03K0omv5DYba2Apv+G3\/fzS\/OZgQ2zC+4ntk3p6jgYFu44iXqe097TTE7lcRpG0b5IhJsTsFhBja9tRH5ymxfSNXyTtOp66rSxc5TLK+UOCGoCDeyq0+krkjhbBxsCDCLCMycs3Xq1Il27NhBqampbNj8SSo1V9SLifo8VoDCN7xpp+e0djRWfzcTI52uO5G21HkjEmxOdSwzsOmvjLyf8jw1\/GCCMm6ALWn4Arc0Y+VEAGLXkEhdkQI29IM\/5AFRTx1RmetVRpxRDyQ9aK24F8TIJ554gj0O0qCBey8gMzHS6bqzm0+z7xEJNv0OY6SHYUBOdbYzKx6hi5+sV2hQpv1gqjRqnVuaBQVsM\/+VRpl\/xNINlUtbXnkRBZsXscjOQKJ+6peDDS+eql9M5XrWokWLCPmTk5PZ2wZWr6paGUiMQrWLTqSZGBmVYAPxmsz9jMXIMNPDzMBmdGCdujqZKv2ymoqVjqX87FwqXr5XRIFNRHwRAZtIvUYL1kpn4yC6\/fbbNW9y63U2PefizxFbgcZMjOSARl\/VxhM9OO30WTMaBsOQFvFipBHHMgIRdqJKZ1ZTjfizbEwZVyrSTXckB9w6\/vX9MZT3+7vKuAG2Mm2fE90QNflFRD+7hkR2VBGwidQrCjYu7m3ZsoU9FoKXcfDeHDiW+j1vIzFR\/cww6tEDz0pn42DlgNLrcGpOavauuJXpX+RSrd286r9rxEi8jIlXUzCQjIwMRsBevXrRsGHDCG+B8XT16lXCKyl4tXLfvn0sHx5jHzhwIOGFTbPkdKEY6WJGV\/e\/Tf0HVfjPdE1zJW8aRyVvGq\/57dKByZRz4k3ltxK1HqBSTZ8VpVXQwCZy6O6Uhkac3y78gZ4gVpwNebHwsQbUYLPjbPo2OLc6c+YMqwe6mBXY1OPnYmmbNm00RwEoj2QmqlqZ\/r2I3XYLSgEbBouXKrEzIE4InqjFS5UbNmwgPD\/7wgsvKF7meFNs1KhRhGsKIO6RI0dYvnbt2tGsWbPYK5lGyelCAWdLXZWsWA\/XlBtERofVaTtGU5VLWzVNGQEp0sEmwoGc0pATxctObQc2NSjQHjiblc5mx2nU3MrMvK8GOIDJOd3w4cMt9UD1IjHT2byK3Y7BBnEAT8ouX76c7rrrLqUc3m8GEQcMGMAe78vKymJA4w\/qcWDt3LmTxo4dS7BEGb0DjQqdLhQYM2DU4GlN2UEB5noQJvvA49Sn6m7NGK+Lu4\/KdVtWqDgbxvLNZwuYOGwmCvMBOaWh3cQ7+e7kUJtbHtXWSA4a\/hs2a5jqkWAgQVKb7vWczIyzcRERRwVqrqVuG3XbAc9MjBTZ9JzQT5+Hcbb8\/HzGuSAawk1GfWoP0RIEKlasGCMQzk8gLq5atYo6duyo1Mcf\/UPZ2bNnU1xcXEB\/nC4UAO2P8\/9SDBq7v6lDn3TfqDkb4\/ra2LqbNe0Y6WORztkg4qKPPBmJwuEEGzZio2TkrqU\/Z9MvfCMLp9qqiXa8uGvx+vWipRFnmzlyEA18fIbyKSRgs0Ipf442MTGRgW3z5s2M+3H5mpeFHoddDjqc3s1GdKFc2r+IcjJWKN36\/bts+vn6IZqzMYhH5f8zPYCzOREjYyu0prJ3bnSzOSll\/DSQiGwGTjcsT4MrJIXB0aC2zJgxQ8MgIFqmp6fb6mxrht1NnZNXKqMNuzUS4iHER1zFB0fDrrNnzx5KSUnRGE34jvTee+\/R+vXrqV69wJDWTheKfvGhbj3Hwi7U+Ph06pe0X7M0YPgA4NTp4KtDqUbpHZrfrLiHk7UmweaESsHPozeGOPHR5OtQDza2hgVC\/ImOzvJQ+9ixYzRixAgGKoiZlStXtrQU6ZVXfWecgu3irglMjFQnPccCUeI3jqD+XQ8zcRMJ4uaR29ZpDoWxW2W9PzoAlF4tkr6CTcfJjfROUelAdCEU5vx6sdPuZgBfhy891ZPaD9YeAV3+5nnK\/y2Dil1fI8Cq7ZVGpmA7ceIEe3gdO8XKlSupYcOGrC27MxBuCjZyp+GDBDF69uxp2vecg8mUd0qri10tdjeV6vy8UmbcG99Sk\/RZGhB9fzieTtHd1GVCgStW37VpdH\/xxQHipr4+UUICbKdPn6YqVapovOFF60H+zNVD6eqFrRRbNZ5yT1+hmLLdKXF4imFVX375JZMw7A5u3fQjWsrwdZgyrye16T1fGXbeqX9RzsGpyt9epR89PQ3BlpaWxo4B8vLymH52yy0FceDtnESdiJGDBg1iC8YsnT++lupd\/qfyGV4fF060opgHC87Gfvr4NWqQt5RiqxSEJs\/9+Qplf3iaqMs4iukyjrakZdPM7Wdp4c1rNWBDvksJk4ia93G9vnJycgjHJeD6JUqUcF0PCn65+Z9058XFBQahMpOoec\/+hnXCQIVNUILNPck52JZN70qN2k8pANfJhRSX+YHyt1fpxxJsMHLgDA2TWbduXVqyZAnVqlVLU2bTpk2eDSQ4OG\/evLkQZ4sp051K3rFIKXNp+3iKyS8gDP9wJTVL4Qzgfm8cOEv\/bLKQWiV8o5QtFt+c4tu94n62\/u\/AFFavn3\/+mapWreo5NuMPhzdRlYxZSn9OJQ6h+i0eN+wfFgocCCTY3E8fB9uqsXdTu8cKxMjv9i2iaplrQwM29WH1vHnzAgwg6AU6GmzTv94aiXb1u4yREQX5YLn8qsx\/Ubeh45mj8o6PU+iZmwsIiDw\/l+pOt3Ra7n62fPb6l9ZIT1MhXJiD7eWl86ltj4Jb5SLzINwoESliJDeGNGrUiJ2TwdRvlDIzM9mhNl6XUXuL+HmorfdlNAKbkREF+T47WZ\/eS1ymnMkZeZm8dfpOOl9nju2DElYE9dNAkrZ9LlX5rWBDkAYSN0vZeRkFbMvGUdvu45SCH32cQk0uzlP+xqbdscNQ5xXb5GRgg\/4BgG3cuJFd9oMIqU8QJ++\/\/352WL1t2zYaP348tW3blvr27UvHjx\/31V3LiGs55WxPHB1CJWr1UcBmBFyADYS0epPMjsJ+gm3\/nIeoGn2gGEisbiU4teja9T+av3Marl\/Yh9r1LbADHH57FVW5Mk\/RnTMu9vd871FNZwY2KPpDhw5lvpBmSX2b1sgRuX\/\/\/kyXKFu2rGkdTheKF7C9kN6TzlQabgk2cL8LLf5h+SaZ3WL0C2wQdbelPE8zz19znkXa3nI+DZt8za1Jn5zS0K7\/0fzdDGzBuB0SALZQEZ4PMmnQ01S9YZLpe85ewIaxAGy4aoOU9WFfyr\/ypWaIsG7GVmxNJerc6fosxS+wcRehDyuNpJrxZ+nklYrU9cyLho7XGIQEm\/fVaga2kOls3odgXwMG+UDyCvqt2aNKZiNvfjjl4p6aOjkVI1GG54VD88VPJ1B8UoJp59yepfgFNiMjjp47qzsfSrDZ3dQOVig79XiN\/CTtHI3tVmLUgK3Xou30e+07FXoYXQp1wtnMDCSomCu2cGjGpdG4+sZXftTAtJsg\/Xe\/wGbk42mlU4YDbBh7KEPZqTm4PhiQ0f030bnjNJw9pTsNGF5wQyQkBhLRzrrNj0E6AZsRkPTmejPTP0RErthC4cXlUiuwGflSOhmfH2Djjq\/6c8DCALZghrIzi2fC58Xuu938cbB1faAN1etUYJXGekk8PFUxVJ2pPNx\/A4ld5\/z6zsE2uNxxdjE0\/rYEqlu7NLW7tbFGdzICkt5cbwY29JU7Let3KvU4AMoD5YZQl+4zXQ3PD7BBX8s58VbAOeDerJvo3eLPGVpLI4WzBTOUnZ2PLSYMecqUKUMdOnQQdpdTg+26pgVWaf3VLj9i1agXV8ija6XMmkzTSu9i3Ob6NhWVvnDdCWIVzjpaxGm99PW7vZUYyXU2K0DCZWvv2X7UY2qBv6UI6vwAm9k1IfRDbeQJt85mJkYGI5Sd3e1wkTkyy8vB1vHPzahEtzmmlms\/YtWEFWwHpnWhnh1iGdi4tz7Xnf4ndxJht9f7MuI7wFa66SLFXG8FJA42oxgl6sHjCOBQ\/f92dbjtB9jMOJuVLhkOzmZ2eTQYoey8iohOwMhp2LlbSyr+5xkMbEabfFB9I5101EseDPLw2t7Uu0eB8zCvD7rTxF1\/Yi5W82uupr43ap+x1QPDCWfDYr7vjwkav0h9\/804iN04\/QIbxvtR44lUu2KWpkmziRYFW+brsyj3TDrFVqpLiQ8W3Eq2Gx++hyOUnRnY9KEP0D+jqFzqBznM\/Ef1YmTbGxPZJv+Xqrs1In1QPEicEN6PPBjkN+\/+lXq01NaGs6Vt5V6n2uVLskA\/DzbZT3XrZGsyQew7lT9IUVidcDYs5LL7BtAdNb8z7b7b3csPsGE3vfzNUhpd8S0Nl\/eLs+ljuST2nSkEODuRDmKk36Hs7NrkE4m19Pjjj2sicqnDkluFtONgG3BHDTp33xpWJdYKEs47sfH9eDaB1tV925OXkX7RhVxnMwIb11Gq51RlgX6gyxlZENUytBOwYbEhn5U10u3u5QfYrAw4fnA2r5Gg7RZ+sELZOTGQqMFWs2ZNFifHaUg7DrZl91ykmrcNoSPNxjPOdt+lD2hBzYJgvoXeGmkGNiyu3\/acZTe0cYlSfU+N7xDqBehEjMRiuxr3vmFdqFN9TCDKuf0Am5MNQ98vETFSz9kQch2h150mO7AFK5SdE71Nz9lEQtpxGq4efZVat7wWsBdSRuNvp9MdtQuitRV6A4kZ2KCz5f6Srgn0o18UarBZLVR+Jqf3pjdaZG4J6gfYjDxljDYWdb9FwIZy0NkuH\/mESjZqLyRC2ulsvB\/BCGWHtrkYiMu56gBSaq8WswNv\/oKOmaeJGmxtOlwL2ItwCDiGyf\/tpEJutyqG2WYWcjEybVNf6tnxWswQnni0KysA6fUYffg3dX38TG505Rc1kZCNiOCWoH6ALdiczSkHM8vnxl3Lj1B26v5wjqX+zegZKv0YrELa6cEWW7GVJpSg3Ybnlq4hB5sRZ8PN6YSubzBgqOMnWnE2iEg5P00PMCygDD8UfjrxRU2M\/8IENrM7baKcze3CKEzlREPa6cGGsarD0xcZsBlxNlgai5d7jIk559\/uTRT7leFcq122jO6pqQvBpF8zLtZSLNVzS5EFFmzOZibeSrAZz5JISLuo5mwgH19cRjEeOXnVLltORE6zHctIhBUBGvL6AbZgWyNFx1TY8zsNaRc1OpuVNRIeJEbx+\/kiUIc8gEKLMyqzVOz6muyTWuE1y+vmmo0EW+GFphHYMJoid5\/NDGw47wLnujdzjOkhtPpg246ziSwFN+HI\/QBbpBtIRGhYmPKqwday0bWHWGD6v+3ULGpVfqcyFLfGMzNaRISBBJ07VawJUX4VqvTzZtNzMbW4aRSBy8uEi3I3CTYv1A5vWf0527slphj65BZqsGH3qJg2knrfcNQ1tSEeIqBP9vZU24uhIo2IEtYPsDk5mNePQRpIRGbVOK9ejHzi6KMsMkCrhKMaP1qrKGduehEyzsYvShp59LvpeO6pBMre\/hWV7lrV0kNEfbPAqp1wgM1KjDRzI5Ngc7NatGX0YPupWFLAC7ZqKcp7i9dqCBnY+Htq8KxGYBuvCfpd03P7mdOokWuXaP1hAZvuQQ11nyPhPpsoDQtLfjXYWtxmbkgTXRN24w8Z2Ozultl1VP8dXtn6aymidajzi4ZH8EOMtDor9MMR2Qs9inJZNdhuv9F8pIUWbH5aD0EeOBE7FRGdLBxpIHFCpaKRR4ItzPMouov5wdmk6T88ky7BFh66K62GA2yRLka6cUR+6KGHqF+\/gscq3Eyr1ZvaqM\/usUO7NiXY7CgU5O+RBjazl3ZCaY20us+m74eT53WdTqHdg5tTp0719GRWkQebneOw04kIVr5wgM1KjMTNha+rrAwIRhQpYAtm3EgrsNldaHWyPoo82Jxc5HRCqGDlCQfYrByRMU4j83+kgQ1Bd6ZMmcKegx4yZAjp77MZ3T3Ti4n6PE7A9uOPP2oulYqsiyIPNr+tkSLEdZI3HGCzOw4x6lOkgA39CEbcSMyVFdg4qGvXrh0QEt3JPCOPBJtTSgUpXzjAZrcBRQrYQhk30g5snCt6eea4yIPN7kpMkDDkuNpIBJuRb54oZwPd83\/LoGLX1xB+HisccSM52NTxH9WTqI8Vyb\/pRVMrMBZ5sNnt4o5REaSM4QCb0RPEfHhXUrOozJ+fC4iGJQI2fZgJ0YN7O2NEMOJGmnE2HovEyOyPfjiNGSkiRhZaR2QJNu0uYRtANrcJlf\/LOwFbiwjYvF6GtANbsOJGmulsRoAzC+yjD5OgJqRTzuY28poZPwiZb6QEm3YK7CyRfvhG6jmbqP+nHdiCFTfSDGy8P\/v371ciIYOqIjEjRTibqLRjJ3S5ApvRm9p4T3vgwIFUqlQpwzaDDTYQJu+3k5R77nO7MRt+FyWsV3ctO3qYAUOEs2Gg0Nn+OLuXilds5avOxvsRjLiRVtZIs3iS+rcArF4ndcrZ3EbL9pWz7dq1i0aNGkWtWrUiuOccOXKENmzYQO3ataNZs2ZR6dKBL33aLS5XCFEVAliQjEKSOak7ksCGC7IJnT817LYo2JyM3SyPG3ctP+JGWoFNzZnM4kdaxYwU4WxuH13xDWyZmZkMaBUqVGA+ahxYO3fupLFjx7Lf7r333oD23ILNKcfyCjZRZTjYnM3MmBFKsHkBaqjKisaMFAGb6AZsN2ZhMRKTDXFx1apV1LFjR6X+7OxsmjhxIpUvX55mz55NcXFxmrbdgg0gwOVQO47lBWwIJFSywXyhOPjBBpsfOpvd5BeV7yIxIwsV2DZt2kTLly\/XKKgYAPQ4cLV9+\/YZutHolXVEtNLrVxCf1KHnuPm7WLnDvoDtcs0pVLZENuWd2qxph0dkFll8wQabjIgsMhvXvE7UZ3NWNwOc6mxh52wY1J49eyglJYXw6IE64dt7771H69evp3r16gVQ69sdkyk99U1q230cU9a58o6MXIHnv13JyKbsrPpUfeACykubbgs29jDH2b2afABRiTp3MgNBXukkyirbj6pVqxZQnxuiBgNsuBCbfymXck9fMTxjU+\/KXjwoxJZx0cvNwbZ2RitKSjQ3qLlZF1bUEhYj7a4\/LFu2LIDr8Q7wQY4bN44ZV6xSbm4u7dixg26++WZqXvs0JZ4pePv6p\/NE1csXlE7NbE3QcSpd3avJdyjnAYqr1YdltKpPnc\/p0lLXV6dOHafFlHz6vuLDj2lX6IePsyj+T+0psa\/xK6EZGRk0efJkT1dMhDtbxAqo1+FtCZ8zCQtvvFWJL3jBBkMWPSqxI1NIwXby5Em2UPbu1T7ha9dJfB\/W7SohXsT+40Rr3o8J+JvXoc9nVrfTfE765jbPjIeval5hnflqDL37hbPaJGdzRiejXBxs+m98TQB4dZP6CB+V2PXIV7DZiZHoDACHfzJdowDEGH4OBg5tl7BRLV26VHI2O0JZfOdge\/bZZ6lGjRoBOQE2\/PM7CYPNrYHE745Ha33S9O995sNFQ2GwuTX9eyeRrAEUENF7JcWMKRAuvVcYbPxQG1Y9tbeI3aG2nHh\/KOBF7\/WnB0WjFhjoIEYGQ1w0o5Aw2FDRtm3baPz48dS2bVvq27cvHT9+3NZdq2hMUWSMQuq93uchWHqZVc9cgc3IEbl\/\/\/4EZ+SyZct6p4SsQVKgCFLAFdiKIB3kkCQFgk4BCbagk1g2IClwjQISbHIlSAqEiAJRB7ZLly7Ryy+\/zA6FYQLGVaFevXrRsGHDNL6ebi7IYs5w1w+O2i+++CK7AaFObuv0shbCOV70G\/EdYfX78MMP2TC6du3KvIgQii7aUlSB7cyZMwS\/TNz2ffjhh6lZs2Z08OBBZknF5L\/wwguKR4HoBVkACb6cCI0NS9fatWsDwCZap9fFGO7x4kxwzJgxjB6DBg1iN0NA42LFitGaNWuofv36XodYqMpHFdgQ\/3DatGmM89x1113KRB06dIhF8x0wYACNHj2asrKyhC7IwhS\/cuVKeuONN5jDc1JSUgDY3F669bKawjneixcv0qRJk+jChQsMYJUrV2ZD0dM6JibGyxALVdmoAVt+fj6b9N27d7N7T2oRD6JWcnIy23EXLFhAqampji\/IqmNf9O7dm9BOenp6ANhC7XkT7vELtAXNAAACmklEQVQeOHCABg8ezGiu3th4yALQeu7cuaYxawoVihx2NmrAZkUP7L7gaImJiQxsmzdvdnxBFmCDTtKzZ09q2bIl+z+ApRcjI8mnNFTjhW4MOsDbSCZpjWRrAK5mABsUd4R88HJB1uy+n5c6\/V6ooRrv4cOHmUvf22+\/Tdhsfv31V2rfvj3zPmrUqJHfw4r4+qKesx07doxGjBjBLJFct\/ByQdYKbEYcDytEHey0QYMGQV00oRhvTk4O042PHj3KjCLgbHDrO3fuHLvhj3g1MJA0bdo0qGONtMqjGmwnTpxgQYogCsLA0bBhQzY\/RRVsoRqvOgTeY489xmgcGxvLaHv69GlmocSRy6JFiwzDHkYaSPzqT9SCLS0tjR0D5OXlMf3slltuUWhqBTa7C7JuOJtdnX5MdijHC4vsnDlzWDyadevWUePGjTVDAL1huTWLVePHeCOxjqgDG8QanHdht61bty4tWbKEatWqpZkbL8YMM7B5qdPLwgnXeBcvXkzvvPOOYTyaUIrNXmjnd9moA5v6YHnevHkBEcJAYC9mejOweanTy6SHa7xbt26lCRMmMG+d1q0Lwj0A\/LDYfvLJJ1FnqYwqsHHjACxhCCQLU79R8nJB1gxsXup0C7Zwjhe6GeLt33jjjZpLxnAAwO+33347PfnkkwHBfN2OtTCUixqwwUIGgG3cuJEQIx4ipD5BnLz\/\/vvZAnB7QdZK33Nbp5uFFAnjBXeDVbJ58+YaayQOtNUGKTfjK4xlogZs8BMcOnQo84U0S+qHGtxekLUCm9s63SysSBnv119\/zY5UIM7Gx8fTPffcQyNHjpSOyG4mVZaRFJAUcEaBqOFszsghc0kKBI8CEmzBo62sWVJAQwEJNrkgJAVCRIH\/BZX\/hXPAo3Q6AAAAAElFTkSuQmCC","height":58,"width":97}}
%---
%[output:282b8119]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANsAAACECAYAAAAZW15iAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQd4FUXXfoEQCEmAVAgk9N5baKJ0QZHipxSpKl1AkCJNFFA+6R0BBSkioIh8NClSpROQFggl1CRAEkJCSIBAAv9zhn\/27u7de+\/u3ZubBHaex0dyd+bszJl595Q5cybbixcvXsAoBgcMDqQ7B7IZYEt3HhsvMDjAOGCAzVgIBgecxAEDbE5itPEagwMG2By0Bh4\/foxRo0Zh8+bNAsVmzZph5syZ8PDwEH5LTU3Ft99+i19++UX4rVChQli2bBlKly4t\/Eam9OXLl7FixQocOHAAUVFRcHFxQaVKldC6dWu0bdsW3t7eir1PSUlhbegdp0+fxsOHD5EzZ05Uq1YNrVq1wvvvvw9PT09J2\/v376Nnz544c+YMqlatiqVLl5rRF4\/RUh0HsfOVJGOAzUHTqgQ2JRDdu3cPffr0YSDgRV4vKSkJs2bNYmAhcCoVajNjxgzUqVNH8vjatWsYPXo0QkJCLI7Mx8cH33\/\/PZo2bYps2bKxegbYHLQQrJAxwOYgHiuBjUjPmzePSRNeTp06he7duyM5OVkRbAQuakP\/2SrFixfHjz\/+iJIlS7KqMTEx+Pzzz3H8+HFbTZlkW7BgARo0aGCAzSa3HFPBAJtj+Ag52OrXr4\/Dhw+jW7duGDduHFMBqfz888\/47rvvULlyZZw7d479JpZs+\/fvR+\/evZlEI0CMGDGCgdXLywtPnz7FP\/\/8g\/Hjx+P27dus7RdffIGBAwciLS0Nc+bMYQDiNIcNG4bmzZszNZbokVpK0nDv3r2sTnBwMObPnw8\/Pz9DsjloHVgjY4DNQUyWg61Xr17Ys2cP8ubNy6SPr6+vBJBkH5FdJAZb0aJF8fXXX+P3339n4Fy0aBGaNGli1kOyx5YsWYIOHTrgrbfeYqCMiIjAJ598AlIjCbxEu2zZsmZtSUUlNXPr1q3sGZe8hhrpoIVghYwBNgfxWA42bjcdOXIEK1euRPXq1QVAPHv2DJ999hmTcKROcsmWP39+EEhJ4tWtW5dJKZJoasrBgweZekqFgDxy5EhBmsrbi6UnvY\/6Gh8fbzhI1DBaRx0DbDqYJ24qBxs5IK5evcok0FdffYVPP\/2USTpa3OSl7N+\/PwYNGsTUQQ42okfSiX4jqTVx4kS4urqq6uFvv\/3GQCOWVpYaXr9+HR9\/\/DEDP3k2J0+ezKQu90aqeaHhjVTDJWkdA2zaeabYQg42kmYPHjxggKIF\/d\/\/\/pdJKlINCXAffvghA6AlsPXt25dJJ7VFDLbVq1czyWipiFVGA2xqOay\/ngE2\/TxkFORgowVPjgmSVPny5WO2Ebnzd+3axf5dpkwZQYo5WrKRE4T20iyVK1euCO82wOagBaCCjAE2FUxSU0UJbOSS5zYYeQ1J+tCGM21gEwDFDg36jcBJqtzFixc122xcRVVjs4nrGjabmtl1TB0DbI7ho6Jkq1GjhuBddHd3Z84QcvmTHZcjRw7BRuKSrUiRImyb4I8\/\/mDOjdmzZ+Pdd98162FoaCiTkrQl0LhxY+ZEuXHjBqNH9ph8D01MQL4XZ3gjHbQAVJAxwKaCSWqqiCWbeN+M76txGtzxQS547pCwts82ZswYBiqSerTPduzYMRbuFR4ezkiS7Ud\/E3inTJkibCdQlMjgwYMle3Tk5Zw6daoQXWLss6mZWcfVMcDmIF5aApvYJU+v4p5JsZNCDDYtESS0GU37aRQvScWIIHHQZKYTGQNsDmKsJbCJN5vpVeSlpBCpxMREFvlBYLQnNlIpvpHo0\/vIi3n06FGLIzNiIx006RrJGGDTyDBL1S2BjdTFoUOHMi9kiRIlmHMkKChIYuNZivo\/f\/48AyeFV8XFxTE7jqJCyINI6qO1qP9Dhw6xSJSTJ0+ytjzqv2XLloonBowIEgctBCtkDLClP4+NNxgcYBywCDY6T0Vf1Z9++knxrBV9OcmTRUc5SC2hDVoKFyKvm1EMDhgcMOeARbDRURAK6SE3svxgI9kZAwYMYGepOnXqhAsXLrBDjhQUO2HCBMlhSYPpBgcMDrzkgCLYKCiV3MYEONp8FYONnhHQuJHNTyFTcCudpaKYQKW9IYPhBgdedw6YgY1czxTDR\/s59erVY1EPYrCRl4vURfnxD+4IIKNdSwDt6z4BxvhfHw6YgY1URMqlQfYYHTak\/4vBtnbtWnbgUClnBkk1suGU8le8Piw1RmpwQJkDErDxTdGGDRuiX79+zHUsBxtFKdAJZAo5ok1VcaFnf\/31F5YvXw46sm8UgwMGB0wcEMBG6iNlgqK4OzpeT\/F2pEIqgY1USSXppVTfYLbBAYMDMgcJRYJTEOzixYuF8J\/0AFtkZCTLh1GwYEFjDgwOZEoOUPAAzxnjyA4yyUYR45RerUuXLsz5wdObaQWbLTWSgEYJbEqVKsUiIPLkyZMug9LLIPoYPHr0yOifnYzM6vyj9BRq01FoYREDm\/iUr6XG\/Bj8zp077XaQkPrZuXNnlgAnMDCQSbfcuXNr6a9T6lLoVXR0tNE\/O7md1fmXrpKNgldv3rxpxlpKm7ZhwwYMHz6cxfVVqVKF7b3Z6\/rnYKNjJ5RJKiAgIFOC7cmTJ7hz546u\/k3ZcR237j9BEe\/cGNnCsc4iR\/TPThypamb0T4U3Ul5FSY3km9oEFHG0iJpN7dcFbGtC7mLAmjCBnSNbFHMo4IzFrArzFitlFP+sBiJb8i5u27YNQ4YMYUdF2rdvzw4yqgnXel3ARkAjwPHyUXBBLPiovL4VImqdUYtF7QCM\/jlIshEZClKWByKTc4WCkSkpqaXyuoBNLtkIaAQ4RxWti5kcU\/Sfswo5SOhOA3I0ZEabnPpHJ9spl6cz++fUIzavC9hoUZPNdjA8AQ1K5XeoCkm0tYCNe4Ap\/M4oJg5QED1FPBUrVsxpbDHApsBqLYvZaTNlpxrJP3DTpk1D4cKFM6K7me6d9OGhwA1b+TUd3XEDbK8J2Jy9sBy9UB1Jj3+AnM0TA2wG2By5jrMELQNsmWiaSI0c\/78wxD9zQQl\/D4fbXHqHqkXNzaiFpXeM6dk+o3hiSDaFWXXEPllm2dTOqIWVnmDRSzujeGKATWHm9O6TOQKs1hZUZpVslm5f5WNxhI1Ee790xItu3nFzc7MLdwbY7GJb+jTSu0+mF6y2RpXZwUb9l4PBEQuc0+CXgRhgs7JSstI+G9ls5+Ng1z6ZXrC+imDjUk8JiLbGS8\/pRAkd\/6JigE0Fx7Ia2PQ4SDLLprYjJIqKqWVVrAGKP6OEtOJ75+QnTpQS1hLQNm\/ezFJx\/Pnnn+xOO0ONtDErWQVs6W1zqV28luplVTWSX\/zBMzkT0CinjfjUP9WxlAmASzgDbCpWUFYBmyNsrqzqjdTTb1sOEvHVwDzdOeUd7dixo7B6+EWN06dPV7w9lcBogO0VAptem8uaZNSzmDmL00uy6ZXo1tRIDqKaNWsqqoBim4zGaclzaQtstp4TbWeq1mJYGK5\/hY+E3kVnSTLqpZveYNMr0W05QcS2V+nSpVlCqdGjR7Nhya86tkeyEaApzSIVuiGI3qFUDLCpkIzOqqJ30RGoziwahYC0u7iToyCq9pvMjtjopZveYNMr0W2BTXw+ksZC1xwPGjTIYWok0eeF8uqIHTHitWOAzVlIUvEevYvu4b7liF3wifAmr\/bj4dXhG3agVHyC295zbumlRjIHhI6jQbbAJnZ+UAJgykcjVxe5tNOqRtK7J02ahB49erDU+JRF4JtvvlG8VssAmwoQOKuKXnWPgEaA48Wz0cfwG7CM\/alnMae3ZNPLX2tg4wuczpCRQ4TbcLRnxiUQ\/40cIFrBRm0pW8DYsWNZZAkBm5IN161b12xYBtj0zrQD2+tV9+SSjYBGgONg05sIKD0lmx422vJGygEkBhe9l7yVZKtRgikCiZIaaMkBopQhrm\/fvoo0DLDpmWUHt9ULNm6z1Uw5jZO5qjEVkjJs6ZWYmV2yOXgaVJPj2wgETi7J6Ldhw4ZhzJgxZo4SA2yqWZv+FUnVm7LjhvAirQl7LIFVL4gNsCnPPYFHvmHOpWz9+vUlDhiiYIAt\/TGk+g16QWEJrHrpGmBTnkICGhW52qkUoWKATTUMnFORQOG5a4rgun\/YbKSmVHSW2uuVmAbYHDP\/hmRzDB8dQmX7ktkos+MLgRaBrWrfycLftqJAzv+5CG5r+pu1NySbQ6ZHNxEDbLpZ6DgC1lz3apwc1P7Z\/f8hu4cLnielIqd3O+b6N8DmuDnSQylTgC05ORkrV65kexxRUVFsc7Bdu3bshhvxxYdKSVopQSvdAeDu7m6RD1klEFku2S63mIWWvYawcakBTPLJ6UiJWiDwgcDm2WAW22PT43gx1Eg9EDO1tQY2W1qLnh4IsZGxsbHs0nra+6Cd\/Ro1auDcuXNso7BIkSKYO3eukHeQrgKmS+wp0SVFbV+4cOGVSj9uSQ1UAzaarGoPJyHYdbcwL7mCPoR79WmqgKpmMjPrPpuavmeGOpbApkZr0dN\/AWx0OI\/2JCiQk3beeaGbSHv27ImuXbuy4M6EhAQGNJJ6FA3g4eHBqmamizX0fp2sqZHWJBufrP8UPISp5ZYKPDztORZNGvdiYNPjeDEkm56lbluyqdFa9PSAge358+dMclH+fjp6zg\/3EWFSLelC++zZs7OjEWfOnNF1ZdR\/Jq5FcOP3UD7IB+PblVfMtW4LLNaeO+LrdGbxKAYKXsThVtbiG8WT9XmxjaiT\/yKOJZRDrF9f5s205XhRO5GGZFPLKeV6aiWbvbGrlnpn84hNYmIik2h0EyOBbePGjXZfhqgGCLbq2HruiK+TJQcHZ2K1744Id6+d\/qqewFs52ArnvoeoJ74C2KxJTC3LxwCbFm6Z17Vls6XXHQ02wUbqIYGNruclBwhtIFIqsSVLlkicJjQkW9f8yoHwRsn82DyguoQbtsBi67kjnBBhu75DgUcmNZA7OKij8uMzPBSLP6P+kRrZ4XIyHsYVgKdPNIr2zIey9UZBbguKHS9alk9mBZvW2Ei1Y5bHUCrlKVFLi+plCm+kvMN0DKJ\/\/\/4MVKRm+vv7M0BZyhFh6T43Tlculeh3+UWBtsBiC2y2nquZlORTI5AS8YdQlTs46Ae5Kni+xmC0Hj1bqEv9d9txFEm\/VhR+e2PAMzT8uiaz2T7a1ASFUu\/itktBrGmzR9NmeWa32dRE\/WvNHSk\/LUA8sHUMx9YcWwNb\/O8TkBp7Ay5+xVhMqyOLRclGV\/8OHToUFNC5cOFClClThr1XD9hokG1+OI1U37LCGNpX98Wc9qYTtYPXXcG6U\/csPl93KhaD14ULz6e3K4bOtQOEv221p4pz9t8R1MDBDU1tOZHkXUOQ7fkOgeaL7C3g3uwloOIX98Ljg78Kzza7t0BS+3kY1jRI+G370Du4sP6h8HeFDzzRcmYAVs6bjWbHX55MZrQaDUPFnt9J5tNW36gySba7d++iQIECNhOVnjhxgmkkWhe5PYtMTXYtoqslM5ZSlL+lTF1q+yzeghIfwXl8cBXiF\/cWyPBziGrp2qqnCLawsDC2DZCWlsbss\/LlTbdmWgObLTWSBtlu+i48LfKG0K\/xzXzRuvxLjyaVxccS8OPxBOFvekZ1eNkcloTxu0xg7FM7P\/rWyW93ezl9IpTz5ATkcd0v0HySsylSyo1lf7\/4fQRwYr0EbP++8a3Qxxd\/z8HhHT64dbCtUKdO5zM437waa9s62QRi1PoA2TpMUz02XjElJQW0VUMaR65cuazOMTm06KOZWcBmTyo7+QBtHVK1teg52GbOnMmO9fAin1uxY8wWTTXPJWCjzWraQ6PJoUviqDNBQaYvNhGk1GMEQMrhJ87xQG1pKyAkJESSmkzcCSU1cljTQIlUsCWZHP2c+ifuA6mBtR5PQb08JrDlCGiLXJVfhmudX\/oVvPbNkIDtQtNpTDrTObaHS1+GaYXe7IeYhGD45w9BcJeHGO89EhV2S8Hm1qALvPouEWjZGhuvSIstOjoaBQsWtHlzJi0sCjjIaLApRebbk8qOeMCP1Fg682Zr4XOwVewyEbVq1RQuTtm1fD6bI17staktvV8CNvFmNR0xF0eNcALUUVJLFi1ahCZNmgh0k5KSGEhp22DixIlwdXU1e6fcnqIK8uMrem02W+3lz8V94B+DKeWWoujxIMHB8caAVLYpTUXuqVx3tQ5oUshNLPY2EtiSnxSGe+4o1OvxBLtqf8\/ykvRJXCHwRf7llH+MlBxIfNvDM9sTi1snYsY70xlgy0HiiFR2NDZ54iBb4JI\/5zxJavAlM2m434D4X3pXYwTmvofIJ7640myvQ69nFsDGnSEVKlRgYCFXv1KJj49nm9oBAQEsz4OWTW01YLO18WvLAdJ6wSkcumpSQ+VgVuoD30\/hzz6+FQ7fA3WE4XMHB\/2w4\/evEJzbZLOtC6+DvG\/OZpPCwXY9ug1WlayJx95xcLvvg+5+LdB4ah2c\/LYTArBDiJm8gxaoOW6thM3y\/osdSLa2PTgh8T5kw3zRirk+LC1QPQ4CZ6Sy484Rnl5BPg4tqew42PgaIacYOcd4cSs7GG5lX4bpOaIwsJENQABbs2YNy6OudM8wqZNt2rRhEmvbtm0YMmQIGjRogPbt2yM8PFxVuJYcSD\/m62Em2X6cNlniRJB7+6yBkS\/GPg9WWDweY23B8mctj3ig0jWTLVSloxfemxvE3P5Jp4bjg4KHBN6HPG2KFh++VAW5p3K9VwtsbB4r1PnA8xN80XIgzq3uhcIepjAuAmpCxe8l979Z+1jY+tDQC+XjI7CdWTFGlRppKVGR2oVmy5bSm8rOFtC0prLjYOMfW2teaLU8sFaPgY2M7V69erFYSEtFfJmBUiByly5dmG2QN29eizTkbnPy5HEVjC+UmAUfS5wIVOfJRwuFBWktCoMW48N9KzD+vin6Q348ht5jaVOaqSg7riP112vItsNPGAcHG9H3i10Mig7h5UhkSYSW\/In1j4DyTWJnbCz8AnvdTAHZ7xRrjRLoZwbU9XffAIVy0WQrAYV+E0tmNWqmHJAlnt\/E\/U0TVYFN76a7LbDpSWXHk7hakmjEK62p7ORq5IPtg5D2dIswtzlc30O+lvMcgTNGw+amtsPe9P\/2jjjrFNEWg4GDZXLgYkHVGhXZF56NeggLUq6K8eMrRIva0zm09iWPmR1v4eOQL1glmzF+bhKi\/TcIaiBJJrFkE9tzVSouQXSenrhQYwh2713CYiJ3p+XCnFQT2L4s0RI7LndlAcpiqSgHmy3Jb88+pRawWUtUpGYd2AKbvansONCsOXnsSWUnVyPpY3PyoK9ga9dscE\/IiqZm\/LbqOBVsJJWq3JkgAOHJ2QTEV+ko2C20mPiC5R2fe6Mti74gUNDzsL+\/w5BKW4VxiaM7aLFWCh+HjlVPKj7ngBQHA5MqKzaQyYnx1BU43sQk5bkaSO\/fNOM4Ku8tZpJ6FX9CnQ+SMN77SwFMcrB9XLwV3F98hsQDQ1D8QnHB8XKzdoREsskjTOSSX63NS\/3kRQvYqA3ZbI8v7INbhUaaN3XVbGprTWWndgPbnlR2nSsXYA4SfhL\/+NcbsGtxSYF39Xo8Zra2o4pTwbZn7xK2IHl5ejUJv52pKaiJSjZRtHsrlG\/6MqU0LTa5dBBHdyg9F7cnGtYWNFdRf6pdFAeK+Qj9JDVwbPB49n63meHwvlVUeFa8wCY063tVABtJvXUlryOs1BWz9nOHrpFEltwMjsDtdwoJIWtyNU5J8ouBRM9tOYC0gk3PwrLljdSayq53797sxAntFSoVsXfTnlR2M8rcQVXPJ4J2teXzCJz9Ld70If1\/W10PT8RtnQo2uQFKHSEnwdmAb5iaSIvZ5\/oCieSKSmqKyp1fOiBooV06MlliM9kCG6lq94t+K9h81hY0f0aexJ01X5iBhcCYuvoq1hf1E1TMlpdcke\/T3Iiv2gkhUzYzL2ZknYM4280UW\/lW9EBUrvMOTn4WIgHq\/SI38fObHoJklduj1AH5iQNxWnMlB5M8dvNxtXbYPMs5m9qOWpRa6dibyo6DjfP43wXbMO9EiDC3g2oFo8aAd7R2x2J954JNdoKZeiW2W7gaWPJmaUHVul7huuCxUyv5xDaVXFVjatKl2UyVpSJOW8D30BYWLSBxcFTK2xyLWkxmrv3lZ+5IPI3vJXji\/Uotmao7t8XfSDpdgNG98u7\/EFf6EnyulEXKxXfgUmU93jniAdoW4IWk4pngcEGNIedM03\/fZfs81C9Ss8Vgk9tUpGb6D1gu2QuS17lUpCUGbLikykHisFXlZEL2prKTg23TybGYem270HtS\/3vVmuiw0TgVbOTtOfaHBwMS+2r7REMMBgJbya3\/w5O\/WwkDfFYjBDHd6wuSj2w+ud3DJZdSe1LVXLqUEBwsco8TvYjbfTydATk3yO7iJSG2FrqVHAWvs7\/hpvsS7PfLJzxrmiMFXxZ\/B88fVMSWzyMlYOKVQkukoE7ZcSh5q7REKr59MjtiahwUPLJKajZ5xHjqciWpLI\/fk9c571YOgw+mvNJgszeV3YL2j1HV44mQI2ZSyHhsu7HZTKNxFNqcCjYS09snFpb0XQwGsqfOjr9n9vWvMCCJZbeiL3\/smn3w3lZToMHbNyjlhW1LZqPVufy4evVd4TktdK\/PPQSwKamyXBUN2z0QBZK3Mk+iHGwu9z9GLdfdKFxqreRZ7fh8+Da4Jh4dvofQLYlYWbAqU0N4oU1tkmwlGm1FvrijEqn49slsyHfXVwCbUt\/EDiB59Ipc8tE7X0ewaQUDjyBZPPAFapYCKLigYptFWLR6LVblNMWqdn02Av06d9JKPnOokXIDlHpFdsvjoaWEcKdj6z3MbKImXS4J2ankDgoOJqJFnkj\/fxvg0mFTGJkS2A4tcBHU1Fpt1oLAtil1GNsHU3JwkGS7Hd4JFMYV6XtaArbAow3QKq0DWpZegp2P\/wapoPKScrElrsY3Q9ty43DF67HwuGHsA1Q77IF9rdaw8YslP0l93jceKiZPJEQOJrHkI8J0ytwvxrR1cvJKwVdesmlFAwdb7yKdULPUC6ZdkXZUcasL5lTohcc+9+AW54vBF5awLR9HFadKtv0TT2L59a3Cl5+++qRKcclFDoLz53dhQ\/07wvjIJhpQvjrLTkXPH\/z8xExyPezqggal8sP75jh4b6tlFWwnpmzGzpkmb2LZ+nuYN3FM\/GfItf4gc3Cc7boUkXUPStRIAtskrwWIzJkd+0teFZ4R2LIf7sok14ugVRIVU6yG1k0KMpOKpIL2PR+FWP++THLLXc\/UtwoDwwVvLEk++YdCLPnofXLARvhvx6S1D15pNVIrGDjYWj35DgFplcC1o8ovzkgkG9\/y0UrfUn2ngk0upqlTpEr1D6rMJBfZLPvjZ0skB7eJ6OvOF9KJzR9J1MTt9ZLY3yR5CGxb8icKHiWf0OqC5KQ6StKV9lPuti+KzZ9HIFe5bYisc4h93cRgaRcRgMaHTuGn4GISIHpfKQePzUNQrud+piauK3HdjNckGQe7JJtJRVJBv0w+J0gnpb7RQigxqALzpm6aEYLQqSZbksDY8usoIUiaXixX1aOK\/41tMQsMsIlmRQ420q5q\/hCMa1j06thscgOUxk+q1JAXbzKwkc208VootuQ3HbwksA3N05iFzdCXnWw+OZjIfU7ll1KTcOp6GWyvfUxgbctLOdHar5GQ0VhpQVM41sNuObFo9Rqk1V+lCJbex28yFZViHsVSjyqTmtinwE3kXPUufu28ggG11fpK8In1QJxfEvbWf4C6SYFI8r2A414PJFKx++XsCG6ZwMb\/v283YOO1cxLJT\/YetzmV+h700UN0m206H0h11j9cJnxs0g6WwMG0rw2wKYCtTJNaTLKRdtXyqzT8UPzFqwO2Wdvns4UgLgSm4f7NsSXXSBZhcckjRrIga0W6Y6xrSbYYCWxLdl6XSA8CU9L1cnDxL4bm4WFYWKSgBAycPqmhSl9++o2kw1LvYrgSOxqPAiIUwVbjjBc+PPUAS9ol4n7pi5I6ee4Ewefwp\/B6FMn21+r+UxzdF5ku4UusuR8prtnNpCKpoO2vFcfbw26y8W86MRYRvqcltEkyB5Ybxmw6JTWcwFh1vK+wj\/jrjjlYmLhSoJF2sDiiN\/9rgE0BbD59PZGrhIugXf1Yu+irAzYlycbVxC8vfspyd8Rli8CVd02BvrQg++bozMJmCGzjQk5KwEjtC92oxPakCl5uhA317piBjVzz3Mlwetk5\/DXKtGHt6RODT3\/dgQMnCzMVlkeOyCVTTOJINHuwFEdrn8LTPN3xPHtBZH9+F66PVuLNG3Go8ksvrOkUxoBIQCPA8eJS5ATOuFcSnvPfSQX9z\/ZWaD\/nFGj89P78fickYCPanfJ2YpJZSQ0nqVq7VFfB2yrncVKYJxKX3zLAZgVspF296ToQYYEuWJOWJsxtW\/94FjnkqOJUm00JbNx1Pj\/mM5wZfw9JrWdLJActyD53JzOvENkjtOkoVuMIbIH3qiE1+glT81YEh0jaE\/3eJSsITgaxmkVMJCcNbV5W8JmP2dkOMAeHkmSKaxGPf6smYLtLG6R4finw3\/XRCjS\/OQ+Fw\/Jg9Vs52e\/y9jc\/2IkTpV2wu0q02byV\/qste\/\/tdwpj48pZiHxrlxnYSI0mya7EP5Kqpf2+F8Am1x7ynkvGxVVPDbBZARutoc4BXbHWvSE2xJlSgLT1v48JFUs4CmvOjfqnxbIxxksiFUhy9c\/XHb+XzoPrZ35GYs09ZoMLTv0Isz4azpwbPxYcZQampqnZWbaq6Ct5JfYaESL6+S69LTgZlKQDeZ0+xSUsjDjH7EUlyeQ9cDeWe1TEH4H90XqdH3yjH+JeAU9sbh+LVo+\/ZX2mvTmSem\/sfohP5uUWxvGkw0ls\/eCcxPEjlm4tw8ajW89DTHKJz8FRHbFNqwQ2seSj+v12jEJo4t9pnIkIAAAP2UlEQVTCuw2wmWOFO0i4Gsk\/+FNyT0COaUeFuU0bURcTy+fJmmCbcPaI5Mvh8mQHCp3axCTXtnpJgs0kVuG2fhAKHsGRL+6YxDXLwUROBrfUW1hcIZBJPXH7I1VSmbeQOxmUFiyFY83wTxNUVDnYXCucxp1Cjxn9gBd+EiBt6pQD99\/cjHp+tzE9x3+Y1Os9fT\/e3GkKRM7T+CL+6XAZP+dLlvSNxkaFVMGuV08ipna42dYBLYQJFaqyrQ+lDwWBjUs+Cmf7KWQcnvr9KyyQ\/NdTcGHRI6dINq2ByGpXMY995AHJ4gBktTTE9eRg4x98r7wxwIjzQlUCW\/X\/NrDnFYptnKpGfh32SPhyUG9IMtAX\/6McOfD2Ixfsip6KiAs1JM4Fev7rW7nR0GsICpX6jRmwcjD1CAlmgyMVssydFLP2f1YMFOwaJbARmFser4tIn9O41qQEPpnnhjd2JwoMI7CQzUX039udX2KPHXi7NHbXCES\/O3\/gxyZNUHlXDpQ7cwflz5r2CtMq38e1ov447XYR7deVEujS2PjHZEBcHO5WvIhtuXozyU+FbMISe66xTXPyOCqBTSz5KFxNbvdlBNio7\/J0dfbmQlFK7uOoHCRcsnFTpWpyKOJXmsDm1b0igpZm0UDkv+eeRoFhUpuE1DB84YaBcXGY+WgvCi9oKFnMR9+6jrkfpGJs8AS2DxK\/MtQMTEnVHqOk710mOeRSidov7ByL+nkmMLuGwJZ9xnXmlqdCrnkCMzk4Ej4MQ53zgXhz52X4Rr\/cu6NCYDn12UFGv+tPbVF\/r0m1CKsSgJ97VMJ7i2LxLO89NDsZafZVuxzgij3VXFAt+bxkbIcbP8Kq3huZ5CbJuuVRaRwI6itpX\/DMMXj\/5s08jqneK9jHRuygaRC1UZB8JNmWHemHhJKNBVXdbf8mRPxx3amSTQlstg6WWhIfShm4OADpBqWOHTtqljwc+J5fDGDeSHJyFYhph0b\/VEfDP02e5iwNtoie2yRfDuISgSF2sh+GuCQzNc5vTxBara8sMJCeH\/skDku7rmVA8RsVawbGLU0TkCfOl6mQSmBb2e8oeOT+nhkr4DvKlB+EXrS3fgL+rOGDGtmSJSoi70RErkL4pe8xXG\/ojqETKkqkFtWhD0Z4Lne8u68cykbeMpt8GgP1Qd43koorPruByvFHMOBBD4xOK4T6p+OYzcAlP9HOduIsk8wk2cnmfZqnh\/COEnuuYnBOsIzLFPFPqvquyqa0EG77NyJ65YRMAzZH5I10FNi6BPSDS9laTLt6dPMCm5+WYaa41iwNtkPvb4DnFlOokxhs7xyvg\/EPVuHdfWUlkuNy+Ric6vcnhvt2Y5uOZXs+QpkwfwkYaSFTLBttJg\/5tqnic2rwaYU+aD6ziBngySZb0fUgim9vJQEyfwmXTLlK3VcEI4FpTbeS+HBdSYmtxttzsMn7RlJx6veBePtMD5Aqc9etFHr8YDoFzoG8v+4l5nH0r5gAzH4sAePhaj4YFO7NvLUUhDw5bz3k3VtKMPLX17yDpBnjMhxsjsobySXkrVu3LOYntSXquGT7LqYtKqUUZh\/LuHzRaLUnP\/yuBQrNsyzYSMVJ7bEMVa+8dI\/zQqrU6WGPUC2uGpZ5nFWUHOTNK9+\/IEL\/agb3aWGS9nwh049ky4mlIv3GVTX6N5247r6ojhnYoj77B0cbXjNTYcVgOVKkABqd9UX1s48UJVdY1QBFIPIPyuUK0RL1lxPZ0K0Gkt46BJfYpYpgp\/Hd\/GAHOuT4EDEpnRTV8Hut\/Zi3lgKVN83PhYqrUoQ+Eth+idngVLDRXX9KRW\/eSJ6LhGjrSTorBxvx+G7wMXy4rhRSb9XK+mCjIOJi807g6YVqknmgr\/u4QcXB9qvCgxQXLNlMZZY+w\/35TcyAQpJv9rjdZntb\/CVEf8bEe0wvr+5XE1\/+8r4ZjdgmEfin727k3NATbdamKYKJJFe\/eQXMVEgO6Oxp0YpSkYON\/i\/e6Bb3b\/9wH\/jlGg7XP1qxbQX5x+hOv99ZlM3F+ZWRa\/19s48N2bR03q55+EVcm3IFZe48FersD7yMWdl2qV6cdLD2+aMoZM9TWHPORGfkjaSB8VttKOPbyJEjJfzQkjeSSzb6IEd1XY1ui+riSUilrA82CiK+NsYTOc55my3m3Q1yY2n7c\/h8vYtFNe5So4ooF30JpfeYSxbSucnhobSY6WWkJmx7\/6VjpvGGYIlHkH4jNZHA1GF7giKY6MtnTXIRoA81y4teMy8pftHpg3CpUjYzIPHKp4qUxYWGJ1E8MrtEhabnRPvAF6lo7PcL4qaWNRsjLZQ\/3\/uXqZkUCF1swwlJHS1g05uk1JYTRG\/eSDFzlW5M0po3koONePxXd2+0\/fU2Sp0y5SDJsmokRX\/Iv7qcedxmqvlLB4m9xZ\/TYl\/dJg6dN\/kog7F8DO77JVuVLGTXWQKbLTARWKzRJ7r\/61wIDXYnSryYishT+JHAvrNuKzS9elRZclbzwb63TqPJ6VSzMZKT5Ub7syh+PorZtPIP1r4iiZj9YpUqyaY3SaktsOnJGylnmxLYtOaN5GAj2rubvmB2btVQU97TLAs22iOK23JK8etOauKDinfgvdZ0p5mYuZa8eWoXs9iuIzWy2DZTdIctNY+\/I84vGT6xplyQ8nffK+BhF9D4+23ZfEreTC75SA2temc5\/nrggXeOFJFIx9BcUfjKf6MqsMklG8WT0sFatcUW2OzJG2lJLZTfpmRP3kgx2GiNUBFrR5kCbEoZkSkbMl244e6uvCDJbe8z1lxN4mqc231fBKXctqiGkV0md52rXQTcrvPJURHzt36kaPfZklx6wGSrn6QKWrP5uJNH7s3kdMkjmT3\/bZxqnopGv5WSSEctYCN6ZLM9u3cMOX3rONRmk19qqGR3iW8Y5Q4Q\/tugQYOEPTWlNOT25I3MEmAT33ZDG4sXLlywmet\/6r4rKDrsuEWbSP5VkS9QW3aZrQXNI1Hkeym22jnjOdkMZNO1W638saHnRxo\/wqdzTGfhxP0iMIaWvI2CkU3RZvtNSZe1gk3PeLWGa8mv7yVv5fTp0zF8+HCIr4OSh2spXfNrT95IMdhIc6Ei1l4yXLLxW2x8fHzYfWxqb7FRih7hE2vLgcBVLYr2kLv21S4OUhP+7BKNJrNbSTYu1bbPzPUIjD8MBgZNP2Nm8zoTbBnFI3vzRorBptT3DAebvfezKUWPaJkcpS+PlvYEtosdw9DhVGHkWWK6bVILDWt1bdl0jnqPJTokGcuGms7p8XqvA9jszRuZ6cFm782jesGmd7ESGMjrV67t2wj8JxzJ+81PZOt9R2Zs\/zqAzd68kZkebDSww4cPY8mSJWY3k1q7UzujwcaBkPvjWnB7\/tjMSZIZgeKIPr0OYNPKJ3kEiaX2Ga5GWrvAXmnvgw8ks4AtucnLy+7d95huJ9U6WVmpvgE289kywOakFbzH\/SJCc93G5\/dNiVwd8eoYl4fwT\/V0BCmH0jDA9oqCzZoaSefQInqaLi1wbxiUIXZTtgmVkdwkH4rvfM7Oxj29aTokqnWV0xioeLwVBALbi28s39yqlbaj6htgsx9sQUtbwqu7KVZS75xoPqltr4OEOho98TCS\/olgi7PA1\/WFv+0ZBNGgQvSsFdeieQVA8feK6\/M+yWnYom+NFn\/nxSoBKJQ7O\/Iej9I8RE5D3H8tRHi7U3fD8OWjlaoiSLTQz8p1uRo5q\/xAdsRGqSjNr94xawabva5\/vR012tvHAXvTEdj3tqzRKqN4ohlsfFM7ICAAEyZMUL2pnTWm4dXrZUYtrMzMyYziiWawERO3bduGIUOGoEGDBmjfvj3Cw8NthmtlZua\/yn3LqIWVmXmaUTyxC2xKgchdunQBBSPnzWs6opCZGf669C2jFlZm5m9G8cQusGVmRhp9k3LAmQtLayCyPXPF3yFPHKSFljN5Iu6XATYts5QF6zpzYVk7z+aofvDo\/r59+5qlRFA7PY7qi9r38XoG2LRyLIvVd+bCsgY2WwdL1bBVfCTHAJsajqVjneTkZKxc+XJPKSoqCnQMqF27dujTp48kjtOew6\/UbTrHN3\/+fPzwww\/w9pbmUrGXZjqyg5HObGCzN28kB2u1atVw+vRpGGpkeq8cK\/RjY2MxePBglnmpc+fOqFGjBs6dO8e8pEWKFMHcuXNRuPDLDUyth18JSLt378bo0aMRGBiomK9QK01nsSqzgE1v3kieHuGbb75hW04G2Jy1ghTeQ7kKx4wZwyRPw4YNhRqhoaHo2bMnunbtioEDByIhIQEDBgxgUk\/N4dfIyEgsXLgQ69atQ2pqKpQudbD3QK0z2JURYHN03kgaA53eXrZsGfvYjRo1ygCbMxaP0jueP3\/OJNehQ4ewePFiiYpHqiVNTvbs2dllD3QTCuVKWbRoEZo0MQUjJyUlYejQoaztxIkT4erqCvFx\/Pfffx\/0nhs3bphJtswcVaMVbAemReNBxFPkC3LFmyMKaJrS9MgbKU81bssbqSVvpJ5Er5oY8\/+VX3kHSWJiIpNoXl5eDGwbN25k0o++kqVLlxZ4RqoiSbqQkBABTDTR06ZNQ9u2bVG7dm32b1q8S5culQBaT7yoPZOmpY0WsJ39LZ7dgcfLm8MLaAKcLSeIPXkj5eCxBjateSMNsGlZSSrq7t+\/n4FtxIgRTKLZe\/iVXmXpLJ8emiqGoKuKFrAR0AhwvFTp6MXuEFBbbIFNa97IMmXKMBOA38sm74c88Y\/WvJEG2NTOrIp6ly9fRv\/+\/ZknktRMf39\/i4AhctYOv9oCm5LEU0NTxTB0VdECNrlkI6AR4NQWW2CzJ2+k\/N2WJJs9eSMNsKmdWRv1IiIimA1GqiA5OOgraQ0waoBhTbK9CmAjHpDNdvNwEorW99CkQlJbNZvapKrTnWpq80aqBZs9eSMNsDkAbGFhYWwbIC0tjdln5cubLiW3ltbB2uFXeyWbLZoOGK5VElokm96+aA3XUps3UtwvS5LNnryRBth0zDg5OWi\/iyRasWLFMHPmTAQFSW0OPc4MS0DVQ1PHcFU1dSbYVHUoHSrZmzfSAJuOyRBvLE+aNMks+xeR1uOmtwQ2PTR1DFdV09cBbEob5lwC1q9f3+wq4IziySvj+ufOkAoVKrB9MnL1KxU9h18tgU0PTVWI0VEpoxaWji5rbmpv3khDsmlmNZCSksIAtmbNGtAFeaRCygupk23atGGb1fYefrVm79lL047hamryOoBNE0OcHC8q7tsrIdkoLrJXr14sFtJSIRDSprabmxvsPfxqDWz20tS6ULTWN8BmzrGM4skrATatC\/B1qp9RCysz8zijeGKALTOvCgf0LaMWlgO6nm4kMoonBtjSbUozB+GMWliZY\/TKvcgonhhgy8yrwgF9y6iF5YCupxuJjOLJ\/wFkMa\/N2pj+pQAAAABJRU5ErkJggg==","height":58,"width":97}}
%---
%[output:50ffdf32]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANsAAACECAYAAAAZW15iAAAAAXNSR0IArs4c6QAAH0NJREFUeF7tXQvYTcXXX4Rc3vS6lXuEovjkVohIF1GSHroguRSVaySUbrqHlFuo\/klRKql89BT66gtfJQm5FJWokHveRLl8z29qjjnzzuw9e+855z3nfWeex4OzZ689s2b9Zq1Zs2ZNvuPHjx8nVxwHHAcSzoF8DmwJ57H7gOMA44ADmxMEx4EkccCBLUmMdp9xHHBgS7AM\/PnnnzR8+HCaN29e7EuXXnopPf3005SRkRH77ciRI\/Twww\/TK6+8EvutfPny9NJLL1GNGjViv2GJ\/d1339HLL79Mn376Kf3yyy9UoEABql27NrVr147at29PJUuWVPbq8OHD7B184+uvv6YDBw5QwYIF6bzzzqMrr7ySOnToQKecckrcu3v27KFevXrRqlWrqG7duvTiiy9moy\/2UVcnwWxOC\/IObAkeJhXYVCDatWsX9e7dm4GAF7leVlYWjRs3joEF4FQVvDN27Fi64IIL4h7\/8MMPNGLECFq+fLm2x6VKlaLHH3+cLrnkEsqXLx+r58BmT0Ac2OzxUklJBTZUnDBhAtMmvKxcuZK6detGf\/zxhxJsABfewR+\/UrVqVZo2bRpVq1aNVf3tt99owIAB9MUXX\/i9yjTbpEmTqFmzZg5svtwKVsGBLRi\/QtUGUI4dO8bezZ8\/P\/v3SSedxP7wAvOQ1xE\/gvrQMvJz\/rvufbyDOiigK+7wqN7Fc7GO2LajR4\/GmiT+LrZT\/IauTijm5aKXHNiSMJgi2CCIEEyAAf\/m5hoXVg4s3iwODFGYZbCIgANgQIPTxTMRLLp3ZVCK9RzY7AiJA5sdPiqpQPAPHjxIJ598ckxrAWBci8CxIWstrvngVOGaEHW4wIsaS\/XRJ554Iu5nUSPiXThJYNqiwCTFv0uUKEG33nor+5trV\/E7OrAtXLiQvvrqK9q7d29Me4JGo0aN6LLLLgvMWZEefxn06tevH4qe2AD08\/nnn2c\/oa9FihQJ3L6oLziwReWgx\/t\/\/\/03zZw5k1q2bElwXKAAYAAABJibkqJWg5Dj\/wAbgCOagnjfC2z8HbFJItgAjDlz5rDHEGKABM6QxYsXs\/b9+uuvdMYZZ9BNN91ERYsWjX1bBBve5YKL95s2bUrnnntu7P21a9fS559\/HgOwiVD70Vu2bFkgevKQiPR53wE4gA99RpEnqUSIhQNbIrgq0Bw9ejQVK1aMmjRpwoQSYEOBaQkg4f+i1kok2F544QXat28f0xQA2Jlnnsk8oHCmwFuJ3+GthCBi+6B06dKsrTLYxo8f76kh\/vrrr5gWgWPGr\/jRE7WSCT2VRsPEgH6hwJuLSQD\/nzx5Mu3cudOBzW+QUvX5\/v376csvv6TNmzezfawGDRrQRx99xPapLrzwQgYy7rrnYJMdGonQbAAbNNp1113HNKcINm5e3X333UwQDx06RNWrV2dgFMEGUw+ahtfjYyBqVWhsAOSpp55ims\/LpNTRU2knE3oqoOG3rl27xvYHZZMSgBsyZEjCxclpNsss\/v777wkCVK5cuZh5VahQIQY8mHANGzakVq1aMQHmXkn8DQBypwSeiWAzXbP5mZGLFi2ib7\/9lm0DQJuhcM2WmZlJV111FT3wwAMMbNByn332GZvxRbCNGTNGuYYC+LBHh8K9kXwNNmzYMC2Xn3zySeM1mQk9\/iG\/NZrfc8tiwcg5sFnkKsyR2bNnU+PGjVnUx8aNG5mg1qtXjz755BP2\/4svvphq1arFfscf7iARtwKwef3QQw\/F1mzctOSueZ1HEUINgZe9kVu3bmUaDaYi1pAQcAjujz\/+GAMb9uaqVKnCQNixY0dm+r766qvZwHbvvfcysPI1KGefCmxYD8FE9FoPYYJQ0VMNiwk9vGcKJLEe+iz3yaJoOLDZZubcuXOZVoDZ+MEHHzDyWCNVqFCBECECZwSEmmsqUWPAnASIIABTpkxh70II8RsiP7z22SDMeA6B54LNzdL33nuPlixZQoULF2ZmIxwYffr0YV7JpUuXMrDNmDGDzj77bAZIOA1AB5EsHIxiOwE2\/g3uMUVb5X02sY4f2FT05LExpWcKNFkD8jVdIgHnNJslxCFK45133qE2bdow50PZsmVZFAaiMSpVqsTWLTDVuLsfwimGXGFtB4BAEwCY\/fv3Z2sMUYvJm9Ny00Uzkr8HLQaQn3rqqcwxAGFE5Erz5s1ZKJYYR4ln0H516tRhWg1bFljfyWBLZc3m52xRDTf6jdjV9evXMzM6UcWBzRJnV69eTVivYa0Gc+faa69llKFVYFbyjWxxbw1g46ahCLbdu3dTv379soFN1iA6sIkAhaMG4AdAQBdmItaMAJuqiCDl+1vptGYD2OBl5Ga4yfDySeqNN95gk0uiigNbRM7C3IIwYj0GoP3000\/MlERcIsCEtRhMNg42BPvWrFmTfZVHlkCjAWx8vTF16lQGQtGMFJvJN8U5ULmJdc8997D1mFz4bA96cO2jnfgjFnkdCAGE0KIfItjgiYRnlWteTkNes0Fr8vhPOIV0BacQVPTk+qb0+HsqZ5GuDUHqRhEXB7Yo3CNiEfY4MrNixQrm2kewL8wwmGByxglscmPNhHWSV4Fgwc0ddKNVJzTcscDNR78uQ0Nj7cY9kzrwiu50sQ4HBn4zMcv4ZGCLHr4r80JcX+K5yNscBxsE5a233mJ7Ko899li28BbMyu+\/\/z7bFIQ7GQv\/22+\/nZ2pgqDxAjpYiGOWw6yKmb1nz54swh0er3QvON+F816YndetW8f61Llz5zgeiH2EMwJ9RwF\/4QVTlTAC4PWO6X4W2gLTc\/78+VqgyBEfcLxAAwKkWPNBZsAHnJfDRr6faWabng5sOidLGF6HkVulZgOQ4MIeNWoUXXHFFWwWEMNuOBDvv\/9+tjbBAhwzOg463nbbbczE4JESWLP07duXna+64YYbmEDi4ONFF13E7GrxAGWYDuT0O88++yxdf\/31hKgJbjbqXPOogy2AHj16EDyXa9asIfAwGWCDQEODYFL0E354MLdv3+6rgQFKAIuHPKEf8OZhIxumI7Q4nDJY+\/l9k4PcFj2VZkspsPFTwDCNsAGKAk0lgw0btJjNW7duTYMHD47F+2GBCS04ffp0trcEdyqAxg8lcmBB4LB+wJ5Q27Ztcxovkb6\/YMGCmFeRR9p7JSwDL8AXxClis1s2b3hjwsy2fu9w4ffb18K6D4DROVGCMAyAhAY3BZwfbVN6KQ82LPIx68KNDS2EI\/cAiAw2aD1EEkCT4Tg+L9u2bWMgBIDgTUMEAkwm7BvBA8YLNm0BUridoT0hdHmlYHsA6zvsewGUI0eOTIpm4x8RI0dUH+brNTg8dOkVgo4VB4jpmtGPvgm9lAcbGI34ue7du1PFihUJ3i0UEWwQEGgkrL\/kfBQ4ZYxOwozCO++++y5NnDhRmUdDR8OP0en8HLzbsWNHzPOIvkDT6czIMH31c6pwZ4kXbayzEPlvs3CA+LXP9Jt+9FQWg9fmua12ebVf643kx\/llsPHf4cpGA0UnB3+GAcUMij+wwQHgMmXKxLUDpgpMMJicWEe44jiQ2zkQGmwyCDmjACKYj9B6cB\/zf8smCUxReCjl7FG5neGuf3mXAykNtp9\/\/pk5HxD65IrjQLI4AE8696bb\/KZVsKnMSJ1m8zMjAbShQ4eyaAu+L6XqOMCI1AM4WRyVQbZo2aKD\/tqiZYtObm8T+ocYVgRl2y6BwZYsBwlAis1hbAJ7hfsA4HA6QPvBwxel2KJliw76YouWLTq5vU3oX8poNjQmGa5\/DrZZs2axQF5dwYlibDkg1i8q2GzRskUHfbZFyxaddG3Tkx\/8SFv2HKLKJQvTsNY545ALrNnAbL6pjZO9PFoEGk+3qQ0giNEiJpvaDmz\/TC+2QGKLTjq2CUB78oPNsfl6WOsqOQK4UGADsOBpRCpsbFQipAv7HqpwLcRPDho0iJ3t6tSpE23atMkoXMuBzYEtyJJAN5m8tnw79X1tfRwpaLcbG5VlgEumxgsFNr5IlgOREX0CQPkFInfp0oUFIxcvXlzLTwe21AMbnFbcQ4wDrnAkRDHd4bSxQYfLo4oWNNrS7\/dlk7NjRUtRkzrV454lWuOl7BEbB7bUAhv3DiMnZG4oCIxfeOq1dKzoP+n6UKDtJt1YK2Hdc2ATWGtrXWOLTtj1kco0itomPvkhDyZyqqRzwYSB0xpZze6mI6XPdmBzmi2cZpPXKNw0sgU2P+9wOoCQy5YMNqzlLqyWmTCPpdNsuUCziZoM7m0ATjaNHNhODLQObOJEkYj1mwNbmoNN5W0ThYavQxzYgoEtEes3B7Y0AxvXYmg2zB5Zk8lmHBb8EJx0ApvuxAnvm99zP1PWRLM5sCm4GFWIRJK2aNmiIztI5q7Zp9wzAuB44QDk\/0\/HNZsfmPyeO7D5cUB67hwk2R0kQ+b+GLceM2FpOpqRfmDye+7HF6fZHNiUMiJqSZVm8xMsB7bsHHJgiwC2B99ZT3v\/LkBnnpbBqIQNOLVl\/tmiY2JGOrDFZ37z4weeO7CFBJuXNy6o+9YWSGzRAUvgEFm3ZRedVT6Tfv39iDMj\/z12xHOMyMmoHNhMOBASbO0mrVTGvoFcUI+SLZDYouPn1jdhKw+6HdiiXKSjSKZraN6mKAG+fmsy1XNkYUZmt1WrVrEmIDu1nJCKt81Es3EvrgmPTeuktevfTxgBNu6dMznHJIMkrMDYAhui1cUNatNBVdUb2KIsdfuvwtpzf359DQI2XRSLafuDgg31ATSkX0TCXBScuXz99deVgDMBW9CJ2qRvaQ02P2HUucF1jPFyRgQxSVMRbB3qZNLIlplKsJmAIwjY5HEJI7hIm4EsbSozkWsxHCrGBZDId4o8l7gfAZdQovA6eC4fPjYBG2gEGfNcCzY+C8MRojo+oeu4n6bjIJmx+hDB8yfuXwURGFtg89PcJgPM68hgMwnxEukHAZvc7jAmmVf2NZ5MGImCASSuCXHHmknspinYgoy5yViklWbjAhLWtEKQqQhOeeYCSKZ9tJEeXLQrG++CMD4VwTbmmip0cSVimk3eQpAtABU4goCNO3aWbNpHzapnhjoVzQG0ZcuWOFOQa6zKlStn03q8jXzwdMBzYPNxkISZ5WUhkhEkAwgg6fPKGpq3PislwOZnJpvMprxOp3qlaXizDAY2eXOca3wvcAQFW5C2edWFOYn76sSCbNp8baZ7V2VaBnGQoG6QCdakv2mj2WwKHmeMPIMDbPBurvjlRPgTrxuE8bY0m+0+9z4\/kx68plY2zQY+bNnzp+f+ZE6BzUSI0TbktcH6jBeuGW+88cbQa7YgY27SzpQFG0zGcc\/PpE5tWtKEXs3YfpOYtMWkc351ZDPSS3sGYXyqgq1drQyaelMdlsoA\/OSarHLJInExlyrHQCqDTeUMQXvvuusuZcZtZ0YKyJCFHmutrXsPxTks\/IBk8hx05\/WtF6vqpUlyAmy2JxgRbCJ\/TLyHqQw29EXeZ8M9cbrU9g5swujbNp+8gCfO4l6aTQamF8100Gxi+2VQqyaWVAebyeSqWrOVP6uudhIPMsGafD8lzcgwzhCTzqrqiAz1+67pvostsPm1J2ifc7NmC8ILp9kkbp33yP9ZNxv9wOanUU1nOltgQ3tt8kEHNpN9sdyq2cSEP7J8mI63KdCdZhPSl\/mBzXRz1hbYkqXZICyiwwT\/l09O5FawdWrTgvVVFRyRJ8DmJ\/SmM4lJPVMzEnt2X49sYkIycgoC\/hHbfJA1myoeUhe6lZvBBn6rAiXyBNhsz+heCBHTl\/0zw5\/ICS+\/l+w1WyLBpst\/r\/NMcrANHDiQkOA0nQvuisd1ZEhlB83WrHqJbOkm0D9TS8aUFylpRtpeq5gyw6+e6Uxny4xMFNhUJ7553+QjS\/z33JoR+foWtRmoVLw2nVz95IY\/T0mwJVOzmTIqyExnC2yJ2mdT5TLhGl42p8QJhuf6R\/\/27dtHpUuXjnQBpS06GBvQglWy80+iSiVwLVSVbEMr5v1Hrn+kHuf9U8mc6eRqKkMpCTbbM7opM\/zqmc50tsBmmw8PXlqaereqkS020qvfKoGz1T9bdPiaS7ytRjVWKn7qNDpoOrD5ISKBz02Zb0uIbGt4XWykV8C2at1iq3+26GDITaJgVJYCxlS3ZjMdb1ORc5rNlFMBZjqbQqTbZ4MgBD3PJ3ojOV2d+Qi26Dywtvpni44p2HSaLU97I22bTwHw5FnVdKazJUR+gdE6IdF1goMtSEo8lTlmq3+26KC\/JiFnKrniJx5UXmjT8TaVL6fZTDmVA5rNLzBaZ\/74gS1Istd0WbOFMSN5vKuOzw5sAcBhu6rpvoutGdvLG8kFAXWgAcUUDrp+5yUHiThWfPMeUSKqVBc6C8KBzTaCDOmlWtS\/KExeGhDrLrjCG1YsHMuulRs1my7yxc8U5zeNqiYtBzZDcNiuZqrV+J7Ptm3btGnjTNumA5HsuPACm3yLzf9sJbrrHX2UjNy2dDcjTXjD13yyheDAZiqplurxJKfDWlc1pmjLjPQSFJNzeKI25pmVkfIBmZVNS7qDTdZs3IsrJiLSab8gE6wJP52DxIdLYWa3ZIBNbpeXGRRlvy5d9tm8jgmJJxpUk6ZuUjMNYjABGuo4sKUp2FQg0Lm\/w26l6NaptiYTW3S46Y4LVtbupsDp8\/K0g8R2TKDpzKOql1OazU8bmYQj8baHBZuu77ZAYouOjXWyc5BEQYmld8PY7TaEyA8gKiDo9pr8aOlYlZfABh6YnFiPIlZpaUbyXX\/T\/aUoDApjtycDbKpJQCcsDmzmEuC3vjOnlL1m2oFNXkc8On8jfbxhlzKxahTG8HdT0Yz02vNTCYufSeo0mw1J8aeRdmCThR9aRJef37\/7\/jVyCmxe69agpm1YsKWjg4TfQBtkq8ZfCuzUyBVg0+Xnt8GioIJtY7HuBY4wZm1YMxJ9SZdAZF0EiQ0ZsEUjJcFmGmLDBVvWbH4XagRhXhjhjrpmM93MNu1HFLCl+6a2KY+SUS802OTrecTGireMHD9+nJYuXUoTJkyg5cuXU6lSpahnz57UrVs3KlasmLaPuut7VWYkQqNmrTtCn23OYnssCDZVZUsKw9CcMCP9wBF0AghrRoJf6QK2RHsSw8iO\/E5osOEK1XHjxlH79u2pSJEicXSbNWtGjRo1Yr8tWbKE+vbtyzIy4RrWdevW0csvv0wXXXQRPfTQQ5SRkaHsh2kIjUqLyO9G0XQ5YUb6gS3MBBA22Wu6gA1yEHZT2waQTGiEAhu0FbTXhg0baOLEiVS8eHHlt\/bu3cuABm2G+hxYuN5nwIAB7Le2bdsq34XmvHry1yRmrFUt2HUmm+iVU2k6CBFAiJtcVJd24FthL\/JLpBkJZgWdAPKCZovKcxOwRK0TCmx\/\/PEHDR8+nJmBo0aNokKFCmkBA3NxypQp1KpVq1idrKwsGjx4MJUsWVL7PsB27ajX6VDN9rH3ws6yXoleQNz24jrqwPsdmTFNFssZ56cpvYQoLM9NBDMqn8Rv2KRl0vYwdUKBbefOnXTLLbfQZZddRv369dN+F6YmNJ98dQ\/XjFjDvfjiiwx0cgHYrhmziP6qfGFksJmkKeNmVvniBejL4Y3YHWZhS9SB9wOHW7NlH5moPA871kHeCwW27777jnr06EG9evUi5BKcM2cO+2abNm2of\/\/+VKFCBfZ\/XNG6bNkyeuGFF6hMmTJx7cKzBQsW0PTp06lq1ezHVwC2jsMn0cH6PWPvRYlAlwVYpJVOmg3MCLNm45v\/yKtocqqbM91ptiBw8q4bCmxwesA8LFiwIANdkyZNaM2aNczxgfXbc889R2eddRYDG0Cj0l6zZ89mHkq\/C+vO7TKKDmVW066fMKNt376dTj\/99GyOGrnrulCcgW9upDdXnri0HvdPP9upRmguB2mT6iNye+Q6uIy+8\/nlArWPt+l\/fz2Jhv\/3z8bvqngRtX\/847bogJ5NWgUKFIiUfFbH3FBgmzdvHj3xxBMMTPA88vLtt98ybde0aVN6+OGH6ZlnnokMtptvvpkBW1cOHz5MMGuhOU8++WRjIRIr4sL6BxedABtydSATVdgStU1oC9okF5i4V9XMoD4XZAZumtim\/vP3Goe3gQ\/gh1ii9o\/TskUH9GzSyszMpBIlSgTmsd8LocCmI3rkyBEGMuyrQWPNmjVLCzYTM7Jz5840Y8YMatiwobYfuKh8x44dVLZs2UjrrGc\/2cZiLOueXoCGX1E1Eq2obdJptigaV2zTsHlbs2nyM0\/LUHpmVd+M2j8+mLbogJ5NWiml2eDgOHjwIBUtWpTy5csXBwSACJoPYFuxYkUkBwnABsA2btxYCzabC2NbtKLS0TlIgrr8RaaJbZIT\/ojrMpOUcFH7J5qRNnK1cDPSFi0\/DRX2eWDNBnV9zz330KZNm7I5PvDsvvvuY8+mTZvG\/o7i+ndgOzGsQT2QskA4sIWFiL33AoMNn4Zp98gjj9Cjjz5KHTt2jGk3Hi2CbQFsCeCmE2xqlytXLi5axHRT24Htn4FuUKEwzetbL5Jpawo2k8zCTrOFA2AosB04cIBGjBhBCxcuZN7I888\/n1atWsVMvlq1atGYMWPotNNOYy16\/\/33adCgQcyR0qlTJ6btTMK1TG+6tDXwNk2RqG2STTndXdhBhlxsk5x+3CsHpXP9B+Gyd91QYANJRJFAwwFguMkRe2vQRHKAsSoQuUuXLiwYWRfmBfoObNtjI2cbbNiw122DmAT0Rp1M3JrNHoCtUHJgSyzYdINkssHvwBZOxENrtnCfM38rL4NNFnieoz8ZIWTOG2kuo0FrOrAJHEulGZubeWKO\/mSAzZmRQSFkXt+BLUXBZntdE2Qi8cswFYSWlyjaomPTuWUOneA1Hdgc2JRSw69Zwpk\/OXmOLZDYouPAFhz4cW\/k5TWbyAhbAhmEjp+TJAgtp9lOcMBpNqfZsuFBzv+iy\/uCYIVkrCNN5m1bE4DJt8LWcWBzYIuTHZODtrYE2xYdZ0aGhf+\/7zkz8h9G2BJIUzp+KSRyok0momTaPxNaiaqT8ppt4MCBLDOXruBYz+LFi6lmzZp0xhlnROKTLVq26KAztmiZ0pm59gjNXPt3HB+blztKI5qfON9nSstvMGzRsckn0KpYsSL7Y7ukLNiQbmHo0KH0+eef2+6zo+fBAaShEPO+oGqhLUup6Ff\/yTN8wwSPP7ZLyoINHQXg8MeV5HHg020F6PFPD8R9cETzU6h5OfOrgZPX2sR8Kc9ptsSw0VE14QDf1EbdsLkzTb6T1+qktGbLa4Ph+pu7OeDAlrvH1\/UuhTjgwJZCg+Gakrs54MCWu8fX9S6FOODAlkKD4ZqSuzngwJa7x9f1LoU44MCWQoPhmpK7OeDAJo2vnMgId8tdc8011Lt377jLQcLeqIp0f7jZZ\/Lkydlu7wlL01REc7JvaOOWLVto9OjR9OGHH7ImX3755SxKqHLlyqZdSOt6DmzC8OHOAITpbNy4kWUKq1+\/fuzCEAjE+PHjYzf0BL1RFUBCDCdSACJCQXXZSFCaQSQvp\/uGwHLccIS+4\/4G8AP8zJ8\/P0voW61atSDdScu6DmzCsCFtOrI9Q\/O0aNEi9uSbb75hF4Z07do1Lvms6Y2qCDnDzT5vvvkmCyyuW7duNrBFuaXVRPJysm\/IMzpkyBD6\/fffGcB4TlGZr3Iqe5N+pVMdB7Z\/R+vYsWNMEHApyNSpU+NMPH7TKmZh3N6DhLSmadX37NnDgIp3OnToQPjO5s2bs4ENM78pzaACltN9W7lyJXXv3p3xV5zEcBkGbrAFX5FhGzfZ5ubiwGYwupiRkU4d1wgBbO+++67xhSEAG9Yp7du3Z5mj8W\/VnXVRbmk16IK2SrL6hoS+MJ1xujuvFgc2g5HH3QQAGxbz0D5RblTVXRAZhaZBF7RVktW3tWvXsvse3n77bcLEsn\/\/fmrZsiVLTX\/OOedE6ULavOvA5jNUuNL49ttvZ55Ivt6IcqOqF9jC3tIaVtqS0Td+69GGDRuYUwSaDXc+7N69m92ClJWVxRwk9erVC9uNtHnPgc1jqLZu3UqDBw8mmIL86mJUzw1gS1bf+LoMDprbbruN8ROXDaLgemZ4KOFowmUsGRnhb3tNB8Q5sGlGaf369Wwb4OjRo2x9htt5ePECm9+NqmE0mx\/NoIKWzL7x22gXLFjALsisXbt2XHPBW3hpp0+fTlWrVg3albSq78AmDRdMHex3YQauUqUKPf3001SpUqW4WlGcGTqwRaFpKnE51bexY8fS3LlzGdhq1KgR19zZs2fThAkTlM9M+5Uu9RzYpJESN5Zx2SPWanKJ4qbXgS0KTVNhy6m+zZ8\/n+688052xZh4ZTPAD+\/sxx9\/nCc8lQ5sgqRyhwG8Y6NGjWKuflXhG9BhblTVgS0KTROw5WTfsDbr06cPVa9ePe4GWmz24\/cGDRrQyJEjqVChQiZdSds6Dmz\/Dh28ZgDYa6+9Ru3atWMmpFxgTl599dVMKMLeqOq13gtL00\/6UqFv0G6IzmnYsGGcNxIb2qLzya8v6fzcge3f0UPsIO4CX7NmjXY8AUJsahcpUoS5sRFtgvXG8uXLmUfN5EZVL7CFpekngKnSt9WrV7PtE5izSFvepk0buuOOO1wgst8AuueOA44DwTjgNFswfrnajgOhOeDAFpp17kXHgWAccGALxi9X23EgNAf+H147pnOlaGmqAAAAAElFTkSuQmCC","height":58,"width":97}}
%---
%[output:51acd8aa]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANsAAACECAYAAAAZW15iAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQmYVcWxrsEBh7DIAIMgoCAIAvpkFfQhsglBIstTlMUAKiICMiIiCERBMMGoKIthkQQwRhRcnvoQlxiXIBHR4AICCUY+EMSMOGyGkfV9f2Nf+\/btc7r73HO3ud3fx5c4t7u6u7r+U9XV3VU5J0+ePEmuOA44DiScAzkObAnnsevAcYBxwIHNCYLjQJI44MCWJEa7bhwHHNgSIAOHDx+miRMn0ssvvxyh3rVrV5o1axZVrFgx8rdjx47R9OnT6Y9\/\/GPkb2eddRYtWbKEzjvvvMjfsK3+xz\/+QcuWLaO\/\/vWvtGvXLsrNzaULLriArrrqKurduzdVrVpVOZMffviBtUEfH3\/8MR08eJDKli1LzZs3p549e1Lfvn2pUqVKUW2\/++47uummm+iTTz6hiy66iH7\/+9\/H0Bfn6FUnAazNaJIObAlYPhXYVCD69ttvafjw4QwEvMj1Dh06RI888ggDC8CpKmjz8MMPU9u2baN+\/te\/\/kV33303rV+\/3nOW1apVo9\/85jfUpUsXysnJYfUc2BIgFG7PlhimqsCGnubOncu0CS8bNmygwYMH0\/fff68EG8CFNvinK\/Xr16dFixZRgwYNWNV\/\/\/vfNGbMGPrggw90TZlme+yxx6h9+\/YObFpuBa\/gNFtw3nm2lMF26aWX0tq1a+mXv\/wl\/epXv2ImIMof\/vAHmjFjBl144YX02Wefsb+Jmu2dd96hm2++mWk0AGL8+PEMrPn5+XTkyBF69913aerUqbR7927WduzYsTR69Gg6fvw4zZ49mwGI0xw3bhxdccUVzIwFPZil0IZvvfUWq9OmTRuaN28eFRQUOM2WAJkASQe2BDBWBtuwYcPoL3\/5C1WuXJlpn+rVq5NYB\/sj7ItEsJ1zzjl0zz330IoVKxg4FyxYQJ07d44ZLfZjixcvpmuvvZY6dOjAQLlz50664YYbCGYkwAvajRs3jmkLExVm5qpVq9hvXPM6MzIBQuHAlhimymDj+6a\/\/e1v9MQTT1CLFi0igDh69CiNHDmSaTiYk1yzValShQBSaLx27doxLQWNZlLWrFnDzFMUAHnChAkRbSq3F7Un+sNYi4uLnYPEhNGWdZxms2SYSXUZbHBAfPHFF0wDTZkyhW688Uam6SDc8FLeeuutdNtttzFzkIMN\/UA74W\/QWvfddx+VK1fOpHt65plnGGhEbeXV8Msvv6ShQ4cy8MOzOXPmTKZ1uTfSpEPnjTThkjMjzbhkWUsGG7TZ\/v37GaAg0L\/+9a+ZpoJpCMBdc801DIBeYLvllluYdjItItieeuopphm9imgyOrCZcjhYPafZgvHNt5UMNgg8HBPQVGeccQbbG8Gd\/+c\/\/5n9\/0aNGkW0WNiaDU4QnKV5lX\/+85+Rvh3YEiAMAkkHtgTwVwU2uOT5HgxeQ2gfHDjjABsAFB0a+BvACVNuy5Yt1ns2bqKa7NnEum7PlgBhcGBLLFNVYGvZsmXEu1ihQgXmDIHLH\/u40047LbJH4prt7LPPZscEzz77LHNuPProo3TllVfGDHzjxo1MS+JIoFOnTsyJsn37dkYP+zH5DE0kIJ\/FOW9kYuXCabYE8FcEm3huxs\/VeJfc8QEXPHdI+J2zTZo0iYEKWg\/nbOvWrWPXvbZt28ZIYu+H\/wZ4H3jggchxAm6JFBYWRp3Rwcv529\/+NnK7xJ2zJUAQJJIObAngsRfYRJc8uuWeSdFJIYLN5gYJDqNxnob7kijuBkkCFjZOkg5scTJQ1dwLbOJhM9rBS4krUgcOHGA3PwDGIHcjVfcbQR\/9wYv5\/vvve87S3Y1MgAB4kHRgSwCvvcAGc\/GOO+5gXshzzz2XOUfq1q0bdZvE69b\/pk2bGDhxvWrv3r1sH4dbIfAgwnz0u\/X\/3nvvsZsoH330EWvLb\/3\/\/Oc\/V74YcDdIEiAU7gZJYpjqqDoOqDjgNJuTC8eBJHHAgS1JjHbdOA44sDkZcBxIEgdCAxvc1KtXr6bf\/e53tHXrVsJjxkGDBtGAAQOofPnySZqO68ZxIH05EArYADRcrMW\/fv36Ubdu3djTkMcff5w9WMRBqwNc+gqBG1lyOBAK2DZv3syeaVx99dXMtc1fIuP+H64cwWXtd\/M8OVN1vTgOpJYDoYDt1VdfpTlz5rDrP\/wGA6bFb5Tjacl1112X2plmcO\/FK6bRsaLtlFtQj\/KvvTeDZ5LdQw8FbF4sxM0FvBjmF2Wzm9XBZn\/w7aVU9NgNkcb5\/aY6wAVjZcpbJQRs2MMh5iCe+mOvhuAzuLunKl999RU999xzzAStU6eOL0NAF3EPcZOdm6pBOZgptAA0AI6XSh2HUsGoJcbTzpR5Gk9IUTHMOcYzDl3b0MGGGIXwQoIB559\/PnsagseRftpv4MCBLHYhgo36lZKSEtqzZw+deeaZcTtcMoXW4TVPUvHCm38C203zCYAzLZkyT9P5qOqFOUfQx4c83o+5apyhgw2XX7FXw61zOEiguWBG8piE8iBgagJsQ4YMiQSp8WI8HlsWFRUxLXn66afHsz7s4Wam0Dr5xmyiL94natCOcq4otJp3Js3TamJC5TDnCLIItmQaXMlmzKGDTewcWgjOEVx89YoOxcEGj2Xr1q19x44Lvt988w3VrFmT8vLybOYZUzcbaGHS2TDPMOeYUZpNlmoE\/nz66adj4tfzehxsusA0qA9z4euvv6ZatWrFDbZsoOV4Ftf3OPTGoWi2+fPn0yuvvMJuj+DJiFgQcOall16ipUuXslslXmZkOoMNpjD+eRXsTxG3H+ZHvBo3TFoYb5j0ShMtOON0Drmw0RYK2HigT9wUwVN\/nqABIa4RExEZWR566KGoDC6ZotkAMoT9RggCV0oPB5CE5MEHH0wq4EIBG2xmBBGFCx+ODjhDEBfjySefpDJlyhA0n5dHMt3NSD4+LEzt2rVLj7Rl8Uzw4cRxlIk1FSabQgEbBgSPEPKRAViI6oTn9n369GEpkbzO2NAuU8CW7IUJc5EdrWgO2MhcmLwLDWxBB2Uz8VQ4NWzGF5QHrl1yOZCqNXVg03g2U7UwyRW\/7OotVWvqwObAll1Is9y6hMkcB7Y0A5tX1lK+6LZ7RzGWvyw4QbPP8DHiYTD25jzbTlj0ZU+1SDfomEUaTrNpsq2AWdmwZ+OCjPkifZP46DaIkPg9c0LUZNBUJaj3+6Ij1N20adPo3nvvZaHxADbVM6p46CNCNK7TIdwfjo54Ac2FCxeyu7RBn20F4WMYGi5Ks+3YsYOdh+HcDLfrcZEYqYwQm1B3F9HvC+qX8shm4tkONj8gegmDH9iCvjfEmkFGEAA2EfR1INX9rgOGjczpaNn8HgEbBgBQIHgovhg1atRgt0Jw+wOAQ+YVP8Ahmi8OsJGeSL7E2aRJE0JAUFWxmbgD22GaOHEii5rM87WpPnKiqWkCBnxgbV7SQ9gvv\/xy1saWPl9vLgtyUFqTDwDqPP\/88yzfuFdwWj8Q2MicDZh0dRnYcEaG8AU4H0MmE1z0RTl58iQLaYDFwHWrVq1aedLD\/UdkZOFRfnUdy3a5yV4k28EGIeEJMyBkKsFEnTvvvDNifunMSKyDTaJFmJDjxo0jJPmAeWdDXyXkeBkC2eGmLP4bMiibj6byZFIvpWCDbYwFatq0aQzjceUKNvmYMWM8bWSAEjY0conh4jEStZsWm4lnAtgeeO1L2vFdCZ1dNY8mdI+9C6rji85BIjoIZGCpaPuZ96hvm9UU9JYtW0aTJ09m+0kb+jKwVOON10TU8Re\/28icCT3TOlpv5IYNG9g7M2Rc8dqQItcYzBvkHbPJ\/Ww78XQH2\/L1e2jU8s0R3k\/oXs8acH77Mi7YsDDk3NdeoPHTPDwdsOhsEOP8YyKy9w9tULgs2NAXgenl4PACG3eMiIIt0hBTG6OOnwMlLcGGW96zZs1idx5hRmLvpSrQjMhaeckll7AsmshDhswsXbp0YdlZoDG9is3E0x1sABoAx8uANjXpsQFqnnnxQ+cEgdDhWpxoZsmCKALPDwwcWEi8KIK3f\/\/+ETCJ2gia7P7772f3X7mH0IY+2qs0twwanRnJx83HKWtM+XeZ1zYyZ6q1TOp5ajaYhsh6CXMBjg+4dr2einNTE6Dr1asXy5CJd2fYw2HieDiqe6mdrmERPvzwQ6bZTfaUsmYD0AA4m6IDm25Pw7\/wHHB+YJD7wguHu+66i0VJ42Digot9Hc7URBMS87Khr4odKrvyTRwkMphAA0Xce4IPa9eujTk+Ea0pKAWVYyipYREAtDfeeIM9LenatSs7U0G2S68CU3PEiBFs34f0RfyJjc1L7XQNi4DARYiFaQI28Ad7tjXb9lH7hlWsTUi014FNt6dRAcjrHIwLLQQOgipqHdV8IcBIISwKtYlm4\/RV8iOPgfFQc\/4ng41rKhxRyWeTqj55fVhtMJPlkrSwCCdOnGBmI\/ZoSCk7depUK4eHPHB8PRCDBF7NFi1axEyMTzxdwyJgfDj6MAWbjRZT1TU51OZml8ockoVf5y2UTVLxy8\/Hx+cuuvz5bzb0VVpZNQc\/0HPNLR8ZyOapH\/B0MpcUzYYjADyRwYtrRDjG2Zpp2HA4SVAX79fEwpnjJaw29nO679niBZqo2QACVZH5KJ9boY24B4r3uhbaw7TExxcfRNwaEc+2bOnLjgyM10s2VLRVySJVfFKZlryejcyFsaacRmTPxvM3A2wwH2F6mIbzQugDHDJi0y4+EoU5iqMAMNjr3MRm4tkAtjAX15aWeDOEtxXvQdocfNv2HbQ+xic7bbh2Xr58ue+eLVnWShTYuDPknnvuobvvvpuFljMFGgjxsAiyI4WHRYCrGiETVDdQHNiCiln47URnCAeWyVle+COxoyg7Q\/gHQrxpI1K0kTm7kfjXZpoNm168qC4uLmbeRJUzBN7ENm3aRDbRyM\/MtRXMT+zLsD\/j3kjcRkFYBISxg3bL9LAIyf4KhrnINrTkczZTs82mj0TUlc3TsO7jhjlWBjbca4R726\/wfQD\/aohgQzseFgGAw00SFxYhzGVytMLkQEo1W5gTsaVlM3G3Z7Plrquv4oCNzIXJQe11rTA7i3fiDmyJXo3soO\/A5h6PZoekp8EsHdgc2NJADLNjCBkPNtlBghxqeGCIK1y46OpVbCbuzMjsAEOiZ2kjc2GOJZQ9GzyUeHyKWw84o7vssssiF5GPHj1KCxYsiEr\/K07AZuIObGEuffbSspG5MLkUCthwDIDrXbjSI8b637VrFzu\/a9asmTvUDnPVHK24OJDRYMNTnEWLFrGnNGIkJP6CG9lIvSI42Uw8GzSb7qW27eG67d1FEylOZCg71V1PcUymN\/v95mEjcyb8MK0Timbz6uzIkSOEK2Bbt241AlvjHa\/RsaLtlFtQT5mkPQywFa+YxvpAOVj2DDpr8EzfNE\/JXhiTW\/82gNPdysf80imUnR+\/+VwAOJu4KbJ8JntNef8JBRtAhvh\/PXr0YMxR3bfkE180vDud+9GCCF9UidrjBRsSwSMhvFjy+01VApvXSfbC+IFN99ZN9dEziX6livmo0wyJCmWn47fu8ayJlhFlruvE+SZNQqmTMLAhpzaCBOH1tknKqNntT6dmh7d4AoFrpIOHDlL5vDzKq3WeL0h4fVFLAmgAnFjwO4Cdf+29SobqFj+UVRCImICtNIey0\/Fb92TLZD14Hw83+po63jjRV45M6JnWSQjYcKEZ5uO7777rGxIBg+QTn9SmAnU+sTFq3OXbD6L8WxbT4ZceoOKVU2PmVKnvZKr0P1Ni\/n54zZNUvPDmn7Tkj\/Xkv4sNvWjZhEUwZbpfPZ0ZWdpD2enAporBYst3EWztelzNZCzqA5yba\/XqxbT\/0MEGDySe6SCcwJw5c6hDhw6RMAmqQfGJj69XRN2rHYqqknPtg+y\/T64Yr55P66uJ1xErsPofPvfTn4R6J9+YTYR\/csmvQ4R6VxRG\/WIbFsGU8V71dA6S0h7Kzg9s\/Dfb8Hsyr0WwXdSpZ4wMJS0sQjzCsnnzZiosLKTjx4+zwK6qMAheE5fNSJh3NWZtpuKFw+jwmj8ph8U1n\/yj3EasV7xiKh1++VSAGFWRNRwWxiYsgsp8teFptoey03kjVSHqZI+rzmMZZUbeOZfKt7+eLRFk42jRdipfqyEVDJxhs2xGdUPTbPyhKAK04uU2wpibFEx86W296K56RbHKpt8p01FlQuLvKicK\/i7vzXg9lYNE7lSmqTNrxPYyfZ3zRcUfnROktIeyU\/Fbjpcphurgv4kh1MGj3bt3ewb\/Ue3Zwlg7nbyHAjbuDIFGQ8SiunXr6vqN\/A5NsPlPM6hmuWMxbSD4KLJTg1f0AhtoigDNa9aRzpr6VgwIlZqt41AqGLUk8pMN2LxAbswMg+haOm9cpoey8+K3F+BQX9zHgtc8booYkk9cAxFs7a\/qz9ZbXjsuMzZrp6sbN9hwcI3rWEjwjmtazZs3j+kTgVtxs0R+Aa7TNGBCyaZ3rMGm8jpCy\/hpSS8A24BNng\/Gzz8YuoXgv5toNr+zsUwPZefHb78Izqr0Ul48NwEb2gaxTPzWOW6w8dDjq1at8uzHK4GdChQiEUw2t8Y5MWdjvI7q6wOtBqHnB9ciiPy0ZBhgY2BeMY0Of\/42lW\/aMZBLWeeNxN3T0hzKzg9sJiHusAa6UA6mYPOynEw\/nHK9uMEWtGMu+PIhs0gPk81rdrkn2OSvj5+mBDArdRziSwv04tmzxcMLWbNlayg7nSXBY6T4aTKVaWlrRqpkId71TSnYMPgdI+vHaCFRc5UtqOdpRsoM0WlKAK5k09u+PEs12OJd0HjaQ0j5zRAZ\/Ejpm46h7OREHxi3nNZK5onTbAGlRDQlZcdIEJLZDLZMDGWncoboUlNlJdh0msgULNz5odqrmdIIa89m21+61c\/EUHbyOZsuyX1Wgs1EE+FwW3Z2yAJqUsdUqLNZs5nyKNPrpR3YEEsSwVUR998kb7Hfuymv6zU6zQYQlW\/W0XfPFvbCO7CFzdH0o5c2YMO52ZtvvsnuN9apU8f4rVOQBPY6sKVimeTjBJ13LBVjdH3Gx4G0ABuS4eE5zMqVKwmJNnS2rzjlIAns0xFsmJN4mOnAFp9gp2PrlINN3Bj37duXkKcNOQBMXvEGTWCvu0GSqoUSTUkHtlStQuL6TQuw4cpV79696eKLL2bXr\/yuBYmsCJrAPiywhekgwbzEa1YObIkT+lRRVoFN5axL2g0SXapVkVFBE9jHXP78rypUpmIunTh0jEo+3Re1Fnk\/\/sb\/qKpjsngiHWU\/P15a5rQc2Ey4mll1+JrO7XGQLm76C\/Zx3T21U8yFh7QEW9AE9iLYyjWoSD+7tHpk1Uo+2RcBHACSd1GVmBUV68g\/clCJgJL7QJtj35TQsT0lP\/XVrCNVnbg6Qi7ZL7UzS2wzc7QcbAtHn6TWLXoSHW\/OXvbLMpOWYAuawJ5d2t36KNNmuWfmsf\/lBSA49PoeUgGE1znyxSH6z9pvY1bcC7gAM35TlSjQXVEYebGd7JfamSm+mTVqEWwXtunOZKhc8WtRH3vIVoUWD1q\/2vDjhOfdSBsz0q8DvwT2+1+9jY4f+T\/P5tBcAKAXQFRgw9cJ9VXA9dKQ4gDQZ07lnpG4FFgYm5fa8YqdLiyCTRg7eSzx5KhWzStT40eKYLv4st4MbDknXouRs\/KNC6l849vjXdJI+9DAFiSB\/fcbxtMPO5\/1nAzAhGIKNj8tqAOuqC3FL1qy92wmT2xsAWfyNCVIXA94sKdNm8aS2u\/du5flYVeFxQvy4fbje7zxI0WwtahRnQ68+CHTarKcnV73GqbdwiqhgC1oAnuYkIe3KoLv\/Dg7Hdi4qcmZoTMTsX\/zAi6nUSavNVXptjLC33QCm+5hqZdQmLzunjt3biRts6lwia8Ewo5PqeO7bk5+c4gyI3NObUUyBmxBE9jrNBvApAMINyVV5qPMcNASzUvVgpw4UJ1OK\/eLyMNP3aKbCqZpvbDjRvLzUzyP8YoijDoIZ3HrrbdS7dq1TYfKwhEgUxFom4BNjBPC+co7kx986vgeT\/zIjAGbKqd20AT2OrCZrjpACQdLGIVrS36LRLfoYfQp0tCZkbZxIxM1fvnNmC7MOebIwa4ak\/wsRjfueOJHimC7tNM1tPfxZ6nceYfSz4wMM4F9WGALW+CxvytbtQ87f9Eteth96xwktnEj4zG3\/OYGcC1btowmT55MiHbldxEddMQ9oe69Ger78Z3\/FmSfKdKG65\/v2Sp2qxnzwU7ani1sIVLRO7hmLB397n+T0ZVVH6Lb1xZs2Iee+M8uKvOz2oE8WWHHjfQCGzfDRMaIwquLxSi\/kPbTbHKgHpG2Kg6kCAivhVO1k+ekow2wXZhziJ21ime8vM9SBbZ01Wyn1x5FFVrdyXhuAzZ4VjEnXoK4jnVOENu4kSbj533yHAK6WIyY3\/33309DhgyJpAjzAxvfNyID7cyZM5kmVGlwERyqcfvFj5S1Je+zf\/\/+dN1110VhVt6z4UeV48yBzUpHBassgsREWHkv8scjyGLpwKYzC+W4kSYOEhlsmLNfLEbMVzQh8d9+YNPNCe3R38KFC30jh4n9tGrVKioIK9qL+0L8f\/Bi7dq1McFaHdiC4SIhrUT3vw3YZM2GMxoAzqboBFN3ZqVqrwOoDDZdBCvQw4sQ0btpotl0HlGkF+N1\/Piuih\/J6+tCj8fs2Wp6e6iDfCz91jql0bXS1YyEDV\/2jBHM\/W8DNjAae7aj366jstXbhr5n42OxiRvJF1\/WHPzvogtePiwX90Cia150+XM6Om8kQvMtWbKEmZ0q8Mt89uO71yG9bJp6AU\/UbK0aesPDgc1GTcRRF97ISu0fsQZbHF2ypjpvpAwI+bwKNLwcA3IwHz5Wkxsp3LTE0yscgOPWiBguwzadsMpBI45D95HTaV9umsqmpazZHNjildgQ2vOvmm7RQ+gq7UgEicWYykngAyU7bDioli9f7rtn8wPbx5UmU+dOw0KbmjMjPViZzWALEosxNIkMSEh2hsj7UJGsiRk5Z3tvKiq4hR4b0CTgiGKblVqwfVVy6m1cnbzYJzgm3MtmsIE\/trEYTXia6Dqyaep16K0D27p9jWnQxxMY0Aa0qRnasEst2J7b8990dtU8alPuzUDMynawBWJahjTSgQ2y890502lC9\/qhzqjUgu2jfc3p7Pw8Kjj5fiCGFeW0o+KCQkLa4vHjx5OJEyFQR65R0jmgA9s3FXpSky7zQh9XRoINJqLOPHz5g1O8uupiPc9U9NB+2lM5kcYObHo+ZkoNHdjCdoxwvqQUbPIhMNR32ypbPYEEW3rdvvOpXu2G1Kvsw75rW1xw6oVtftGjnvVyq7Uj0NxVUj2GHtoX5bSldevW0ezZs51myxQkGYzTD2xwjDS+ZGKoe7W0ABsGIR4CY6Jrtu2jh1q+RrVOfBLFtq\/LXER3\/r07tW9YhdnSvJ1YaVdJNaqdtzfqQFlVD23kQ2exnvhbNrr+DeQ1o6vwNV36wNXU\/KyvInPBh\/fvuUNC36ulDdjSfdX4whQWFlLbtm3TfbhufAYcSNU+PKVmpAFfUl4FIdnhIIE56Urp4QA+nLgNg3wWySoObAacBuDwz5XSwwGALJlAA+cc2EqP\/LiZpDkHHNjSfIHc8EoPBxzYSs9aupmkOQcc2NJ8gdzwSg8Hsh5siOT8xBNPsENruISrVatGffr0oeHDh1NBQUFkpZGD7r333mNvudavX8\/qISz54MGDqUKFCp4S4ZcuOShNW\/FL5Rwx1h07djDP3+uvv86G3q1bN+bhRUySbCpZDTakusL5GW64Dxw4kFq2bEmfffYZi68BQZgzZ04kaClAM2rUKHbWhiAyn3\/+OavXoUMHFoK7YsXohB0m6ZJtaQYRzFTPEeeUCEkOzx8CBIEv4GuZMmVo0aJF1KBBgyDTysg2WQ02PNWfNGkSzZs3j0X25WXjxo2EeBjXX389jR49mvbt28eABm2GV9AcWIgEPWbMGPa3K6+8MtLeJF1ycXGxFc2g0pXKOR48eJDGjRtHBw4cYACrUaMGm4bM35ycn+6gBp1nJrTLWrAhjTEEAKYhojqJT\/x5JlV8fRF6DWmjYC4uWLCAOnfuHFnXQ4cO0R133MHa3nfffVSuXDkyTZeML74pzaCClOo5IpXY0KFDGZ\/Fjxl\/2An+zpgxw9cMDzr3dGyXtWDzWwx8iaHR8vPzGdhefPFFpv14wBrelucSxx6O5x4H2EzSJT\/99NPGNBMhOMmaI\/bD4E2tWrUSMY2MounAplgumIcAGzbx0D6IJoX4g4sXL45ymqApfnvllVdo6dKlVL9+7GNDr9Bz8dAMQ8KSNcdNmzaxPe3zzz9P+MDs37+fOnbsSLfffjs1bdo0jKlkDA0HNmmpkLIY2VzgieT7DL9YjbqYjH5ggynJNaI4DB3NeKUrGXNEshXsh7ds2cKcItBs\/fr1Y3nc8NGCCQ4HSYsWLeKdTsa0d2ATlmrnzp1sDwZTcP78+dSoUSP2a2kCW7LmKIbkGzFiBONrbu6pNM579uxhHko4nJBGSvbkZgx6LAfqwPYjwzZv3syOAY4fP872Uk2a\/BRVyQ9s8ZiRXppNR9NyjSPVkznHY8eO0fTp05mJjb3uBRdcEDVs8HjlypWe5nfQOaZzu6wHG0wcnHfhy1uvXj2WFLBu3bpRaxaPM8MLqPHQtBWoVM0RGWlfeOEFZUbTRJvKtjxKRv2sB5t4sIxAn+KtEb4A8bjpvcAWD01bwUjVHFetWkVjx45lN3QQw58XgB8e27fffjurPJVZDTbuKIBXDOdkcPWrCj+AxiZfvC3idagt0vACWzw0bcCWyjlib4bYjQ0bNoziGw798XdkopkyZQo7n8yGkrVgg7cMAEN4aiRggAkpF5iTvXr1YsKwevVq5q5u374986pt27bN97oWp+W33wtK01Qw02GO0G75r9+yAAAAg0lEQVTwSrZu3TrKG4kDbdEJZTqnTK6XtWDDncFhw4axu5BeRcyCoro0PGjQIHYZuXLlyp40\/MAWlKapwKXLHD\/99FN2jAJzNi8vj3r06EEjR450F5FNF9LVcxxwHLDjQNZqNjs2udqOA\/FzwIEtfh46Co4DRhxwYDNik6vkOBA\/B\/4fGQdHghDTY\/QAAAAASUVORK5CYII=","height":58,"width":97}}
%---
%[output:73b57dff]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANsAAACECAYAAAAZW15iAAAAAXNSR0IArs4c6QAAHWZJREFUeF7tXQm0jtW7f+g4HRlylDIPCZGuMcqVsYluhpVKkzlzhHAMZW4wZpalQt2I4rZc9delKESyDiGSInOd6uAccXI49\/62uz\/729877Pf73vcb917rrMX37ncPz35+7\/PsZ+\/92\/ny8vLySCctAS0BzyWQT4PNcxnrCrQEmAQ02LQiaAmESQIabGEStK5GS0CDzQMduHDhAqWlpdGaNWt8pd9\/\/\/00ffp0Kly4sO+33NxcmjBhAr333nu+30qXLk3vvvsuValSxfcbptUHDx6kJUuW0Ndff00nTpygpKQkqlmzJj366KPUtm1bKl68uGFPcnJy2DuoY9euXZSVlUUFChSg2rVr0yOPPELt27enIkWK+L37119\/Uffu3Wn37t1Uq1YtevvttwPKF\/tolscD0cZ0kRpsHgyfEdiMQPTHH39Qz549GQh4kvNlZ2fTjBkzGFgATqOEd6ZNm0YNGzb0e\/zLL7\/QiBEjaMeOHaa9vOmmm+i1116jli1bUr58+Vg+DTYPlELP2bwRqhHYUNPs2bOZNeEpPT2dOnXqROfPnzcEG8CFd\/BnlypVqkQLFy6kypUrs6y\/\/\/47DRgwgL799lu7V5llmzt3LjVu3FiDzVZawWfQli142Zm+KYOtUaNGtHXrVnruuefo5ZdfZi4g0jvvvEMTJ06ku+66i\/bs2cN+Ey3bpk2b6Pnnn2cWDYAYOnQoA2tqair9888\/9NVXX9HYsWPp5MmT7N1BgwZR\/\/796fLlyzRz5kwGIF7mkCFD6IEHHmBuLMqDWwpr+OWXX7I8d999N82ZM4dKlCihLZsHOoEiNdg8EKwMth49etAXX3xBRYsWZdbn5ptvJjEP5keYF4lgq1ChAr3yyiu0YsUKBs4FCxZQixYtAlqL+diiRYvoiSeeoCZNmjBQHjt2jLp27UpwIwFelF2tWrWAd+Giws1cu3Yte8Ytr3YjPVAKDTZvhCqDjc+bvvnmG1q6dCnVqVPHB4hLly5R3759mYWDO8ktW7FixQgghcW75557mJWCRVNJmzdvZu4pEoA8fPhwnzWV3xetJ+pDWzMzM3WAREXQDvNoy+ZQYCrZZbAhAPHzzz8zCzR69Gjq1q0bs3RQbkQp+\/TpQy+88AJzBznYUA+sE36D1Ro\/fjwlJyerVE8ffvghA41orcxePHz4MHXp0oWBH5HN119\/nVldHo1UqVBHI1WkpN1INSk5zCWDDdbs7NmzDFBQ6FdffZVZKriGAFyHDh0YAM3A1qtXL2adVJMItg8++IBZRrMkuowabKoSDi6ftmzByc3yLRlsUHgEJmCpbrzxRjY3Qjh\/\/fr17N9Vq1b1WTG3LRuCIFhLM0s\/\/fSTr24NNg+UQShSg80D+RqBDSF5PgdD1BDWBwvOWMAGAMWABn4DOOHKHThwwPGcjbuoKnM2Ma+es3mgDBps3grVCGx169b1RRcLFSrEgiEI+WMed9111\/nmSNyylS9fni0TfPTRRyy48eabb1Lr1q0DGr53715mJbEk0Lx5cxZEOXLkCCsP8zF5DU0sQF6L09FIb\/VCWzYP5CuCTVw34+tqvEoe+EAIngckrNbZRo4cyUAFq4d1tu3bt7PtXocOHWJFYu6H\/wO8b7zxhm85AbtEBg4c6LdGhyjn5MmTfbtL9DqbB4ogFanB5oGMzcAmhuRRLY9MikEKEWxOdpBgMRrradgviaR3kHgwsCEWqcEWogCNXjcDm7jYjPcQpcQWqXPnzrGdHwBjMHsjjfY3onzUhyjmtm3bTHup90Z6oAAmRWqweSBrM7DBXRw8eDCLQt52220sOFKuXDm\/3SRmu\/737dvHwIntVX\/++Sebx2FXCCKIcB+tdv1v2bKF7UTZuXMne5fv+n\/44YcNTwzoHSQeKIXeQeKNUHWpWgJGEtCWTeuFlkCYJKDBFiZB62q0BDTYtA5oCYRJAhpsYRK0rkZLQINN64CWQJgkoMEWJkHrarQENNi0DjAJZK4YR7kZRyipREVKfWKMlooHEtBg80CosVZk1sbFlDG3q6\/ZqY+P1YDzYBAjDrbjx4\/Txx9\/TI899hiVLVvWsovYKwjeQ+xk56Q5wcokEcqCbFT6CaABcDwVadaFSvR7N0C0KmWpjke0lqXa\/mDyRRxs2Lf39NNPk92JYnTu4sWLdOrUKSpVqhSlpKQE01\/fO4lQlqrMZMsGoAFwckoUmYWkWBYvuwY2HIQEAzCOkeDAI6xP06ZN6aWXXiKczTJLGmz+knFToVXBtmzHadq9II3q5eyindfXZi7k8IcqabC5jDpXwIaNtzjoCLDBSt13333MAuFgJNijwLXBj37I7ddgizzY+i3bTwAcT0\/dXZLmPlVdgy0awYbd5GBowvksHIjkNNbgpAe99p133skONV5\/\/fUBzddgizzYADQAjicADYDTbqS7aHPFsuHoPshHwRglXwgBGjdwzRtdzoCueAU2BF7wZ5YwQQfXPvgZQ53\/RWtZ6Ltq295Yd4S2\/HyG\/r1yMRr+UEVDsamWpaKikS4LwTi7gJxKP5zkcQVsZhXi6D5YfX\/88cewgg0gA1U3aAN00hIwkgAuIZkyZUpYAecp2AAycGu0atXKlJXXC8vGy4Qwy5Qpo7VNS8BPAvgI4y4ElQi4m6LzDGycAyMjI4Pmz5\/PuBGNEgcG3E3cM2aVEKk7ffo03XrrrVSwYEHTrN999x2j3w63MN0cGF2WdxLgOofIuRGBLdZwQ13HNWq9J2ADVzzcR9yyIl5FZAW2zp07+\/jpzcSM5QWAF+Q2RsEW\/h4u8QP9gAabdwobyyVzsOFySlCnywnzeNV7FZzIwXWwIQIJnnko\/KxZs9jNKjw6aQU28GvUr1\/fsu1YYvjtt9+oZMmSlkENCBN03hpsTlQhcfJysJnpXExYtv379zN+QtwPNnXqVHZbi13ycs6mwWYn\/cR87kTn3JSQa5YNl+vhNhbcQQZ+ebBHqSQnHVfdXeGkTJU22uUxu2mUvyeCnjNXdezYkZ588km7oi2fg4j1rbfeMsyDOXAw5UN2y5YtY7fZwDNRLR9twcUgeM9qPu2kw3wcxXfcuDEn3PrB2+8K2HgwBBYNfjDo2VSTk45HO9jQZ1nZ5P65DTaUL69h8ssycGuOU8DhDgIkvAcAqZbvJti4jDA\/B92fuHbLPzDBfkzQNyc6p6rHKvlCBlteXh7bjoUwO7Zp1a5dO6BeXByBnSWgzZaTk47HIti41eNA5HefuWXZjMCAusyAYqUUaNukSZMIwSoouFUZ8jM3wWbXdrvndorvROfsynLyPGSw4YKItLQ031WxRpVbmX4nHY9lsIF8FezE\/KvdsmVL2rBhAwskIfHrmkQXTHYT5Tx2YHDq1sEiLlmyhEaNGsVcQSflc7BhAzrcZiQjwlnZNZTzqFhl5Fm1ahW7b9yMnNYKBE50zgmY7PKGDDa7CuyeO+l4LIIN\/eOXXEAxRLZh7grJ1s9I0XkeDlor6wVlHDZsGLs4Q3TB7MYCLiRuwOEXL1q5kXL5\/MMgfhBQHm7G4a6g0Vgjz\/Lly32usPyOXZuDee5E54Ip3+wdDTa4XOsO09G\/LlL54imGR0vsBG4XIBEtOwcbFlPF20T5Fx1RXOSHtyACy6gNVgESI6ti1w+Uh2NRfKHXSflGwJQ\/EDKwzPpk5hrbtV\/1uQabzXW0EKSqZdu85kPqNHCE0jqbvOMdm3CNznLZzXUADj4vE11BDqJ69er53Vctz9mMFBPrlVagMbM8HNBojxg8kcEjR0nHjRtHY8aM8blmTso3m7OJAMM9A\/yecLMAh1mdRsAXyxCvNka\/rQIoGmwugm39632o58J1SmBTPcsVLNi4u4ezfnCncGsM9ovCqolbhYzcRPEKXqM5kNWciisfB5QMBtGSoh1QwE2bNvlZW6flo43y3d+yNTPyAmTQiK6nkdzliK5ch13EV4PNRbCFYtnMznKFAjZxHsLBZmfZ5PqMwuFWYBAViruljRo18lsKwPscILILaTUnxDNZYVUsm1EwQw7lqwRIZDCJ\/eByg8y3bt1quO6nweYi2JwKE3O2zYfOUOPbcZYrkA7Abi5gFOAQ3xFBgd9h2azmbEabY80U3Gx+IwcauKXr1atXgPWB8g4ZMoRws6m8pqVavsqczcpKifKwC+3LYOPjbRTRNarTqX7Yjb\/qcx0gUZWURT4rsPGB5a6SGI3kLh7\/DVwtWBRHMpoDykpoFS3E3AjKJ7p1Yt2ogwNPDvnzrjopn1soEcz4jbvPALFRpNFI8UVXU95yxz8a8lxWdk+tgKfBFkHLFire7KKRRtu15HU22eIYlSmvV1pFC+12WPDy4VryJO82cVI+dyPFdTaj9VU5kIG6zfawynNWo3mr2dgZuZY8rwZbDIMtVLB6\/T4smhxpRJ3yuprX7fCifHnXiwgovsdT3qsZdWDD\/c5z5syhefPmKa3SG32FeMeN5gnBfGVUQ\/+REqYXyuRWmfKX3i5i51a94ShHDoYYRXbFdkRKPwLmbNjriG1EWOMBIYoZUY8sRIATu\/7bt28fcPCuevXqhPubQ52sarCFprqyW2jnaoZWW3jflt1Ttz7wbvbCD2wgygGFwcqVKxkrk5PjDNhyA55Ifim7aiOdfGU02FSlqvNZScCJzrkpSR\/YxEgVrNOVK1fYPjkVywZriK8kmJDheuJMm2py0nENNlWp6nxRDzYckwHpToMGDdiRGdU9anznf6FChWj8+PGUnJysPNoabMqi0hldkoATnXOpSlaM6Tqb3cKi2Agc8uvRowfde++9hLNrYC06d+4cIbzdv39\/qlGjhmmbnXRcWzY3hz5xy3Kic25KyRWwgRIBi6gAXZs2bah169Y+rn+4p1YMW046rsHm5tAnbllOdM5NKbkCtvT0dOrduze7saZDhw4+Ni1wPOJofoECBRjgjOjBxB0WmjfSzaHVZZlJgOtc1PBGOnEjrYYVHZoxYwaBNsyIbYt3XPNGanCESwJc56KGN9Ip2BAkwUp9\/vz5\/WQmH\/WQBco7rnkjw6Vquh47nQs7b6QTsIG6DpwQWGMTacaxJIClAABOZkniQ+7Ef47WOVsweyPdIvxRpZpThRjGI5JUdlwfzNqrurM\/qkP\/cuOcgA2HDkG+gh0kmKNxnnTOJYlTyuG8n80JgFWV0Cqfyq5\/eYe\/W2AzWp5RORNm1p9IU9lZjR3vl3yawekYhls\/ePscB0i4YuECRG6twMGPeRnmZzwaefjwYXr\/\/fdZcATWze5iDRX24mi3bBCqzBupqeycQcEOCG4QAtnV4azF6rldARuqk+\/Uxonkdu3asZtHcRGGXWTITbBt\/fJjerb7UCVaBHVRmee0smzyplhNZbeNXQXNk3wuzQ4IdjEAlfG0q0OljGDyxOXh0a9WDqUuwz+OCrBhYOGS821vmsruab9xkflD7IAgH0gNRumt6shcMY5yM45QUomKlPrEmGCKN30nLsEWKcuGU8lGSVPZXSWoVaGyswICf2a1o18FHbychT0fovvT5vteydq4mDLmdvX9P\/Xxsa4CLi7BtvKzjTS8n\/qVUaF+zazcSE1ld42EVYXKzi4aaXQsyAmNHZDE65hW9RQ165bmA9TJsc3p4r6NPrAVadaFSvR7VwW\/SnmiCmzVjq4zNeEcEFkFbqTSnV43vZ8N+Wb+zwFaumG7khvpxtdMhfBHU9m9zQ4h21HZGVk2+YMlnrx2SmMng63xox0ZoGQ9QL64Bdvc9tWo2tF\/GX5VZEGk3NmMCtZoFmDieb7dWSk05GApJbDBbcB7PAUjYDuwaSq7a\/TisgmQqezM3EgzwDmlsTMDm6wHcQu2xS+0oWEVMwJMMfeZjQSBzLJPzfM5AZsMZHzlADgnyQ5smsru6pzNKMl07FZzNu4uiq4kz+9ksVt0I7llSxiwpY98gB66KTtgLLiVkX1pMyvE8zkBG8qC63nhh42G1lIFdCqL2prKTo3KzgpsZhR3TmjsEt6y\/WtIK+pU+kyAXsNdLNKss1+ESMyE56XHfsl+Ei2UU7CpAMoqTzDbtTSV3TWJimusdqF\/q4sSeYlWNHZmYNs5oSOlfn\/1IshQphRWehLxAAmsyv7\/nEglk3OD0nlu\/UQ3INxgC6rhYXwpXqnsgqGxMwLb+gav0e9zu9Cj59fFL9iMIkDB6CDmbliI5IEODbZAKcYrlZ1TGjsjsI0tPoyarX2K6uVcvZgypixbdnY2u3oIXCZWdzqL1ijl34pR\/sJJvo5eyc6li99fcy35c\/l3vIDVfoCN59l5iOiFz67egGnGmx8MqGP9nXilsnNCY2cEtt9LNaSCy\/oEDG8wkemwupGZmZk0ceJEWr16teUdWWgUXMjMlWMpuXJhuqHRzQHt\/OfnbPp76x8Bz3N\/u0i5py\/6gVEsA2DrNSefBlusfx08aj+fF85ulUUNavwHXco44reYzas9+NAMerjHi661wrU5GzYi46gNom6\/\/vora6AdCSi3bAAawGKULu4+Q0klUyjp1pSAxxx0eID3uWXUYHNNP+KyIA62t\/rnUa1CxSh7\/S5fP0UPKuOWXlSr19WLTtxIroGNm3Ls9n\/22WdpxYoVNGDAAFs38tJf\/+UHFLlTcBtF91Kl0xpsKlJK3Dwi2O7Kd9V7QgLQUmoV8wmmQPF2VKTxDNcE5RrYcFL75MmT7PgE3wOHg6RWc7bzO6dSzom5rnWGF6TB5rpI46pAI7AZTWWuL9eBCtWZ4lrfXQOb2CLVk8Ln04dSzrGPXOuMBpvroozLAmWwwXsSpyG80xpsCsOvLZuCkBI4iwi2OiXNpykabApKosGmIKQEziKCrd7t5oLQYFNQEg02BSElcBYNNhcHX4PNRWHGYVEabC4OarjBFsxGZDeo7ESRGd1\/HQp9ABQyXvkjNdiIKP8NZZn+XPn7eEjQixTY0GiZyk7exe729bpm57lUdsdbCTme+SNVwYawP+ZtbqWIhv4R9kf4nyd07uK+TXT5n\/8OqX\/RBDYveSPlg5ey0OyemwlZ3k1vRdgrP8P\/sd4qf3ScDqjdUZtQ+CNVwFaw2kAqWM29rVrof0TBhgb8tGEoHdn9ETV+5GrnsHMfABS3b2FbFk9YE0m56w5L6xeNYAM\/Ik4ru8kbqaJwyFOkSBFq3rw5u4tBJWGddMmSJTRq1Cj2jh3YRHBxsKFObARHkrkh8ZtM7BNO\/kg7sHkBNM\/ApjKgPI\/RFwwnrum6XWxPJDYc554qRgXvbEYXTh2i3PL16IZbsik5\/wbTaqIJbOifF7yRdlQMTsZAzguA4opnTmVgBjaActiwYTR58mSqUqUKK4bPHUWaAvmjYDTm4eSP5PW\/PaYh1Urd5ut+0k33UIGb\/+8EgMsWjVfgiWVzMtBGgjfjBeH040VOT6P8GWt91eB0ACwewIkEsPVbWTBsu\/7tAiRe8EYG6yKqjA0A07RpU9\/xJKPgCy9HtkhGwJRZoSPNHynqHMB26Y\/tnoIsqsGGxhnxgnCwFb+0hS4fePka2H4qTH9v2+v7\/76Cd9DAzTnKYLvw45t05e8TlP+GMkF91SLBG2kGNpFx2QoQ4u034rk\/o1PdZpaN14V6OOOz2ZxNBFik+SPt5oMqH6Rg8kSlZTPriHixRt6vC3xfJLiZWV8PYqcDYOF+KNCRei5cpwQ2OUgTjL9u59KJlNk4FdG9e3fmookHW+WvP2TA95hiTiTPfezqFN103AjLL0GRwcDrmDp1KmsPFBFHpUQ2LKs5m8y9b8b\/IVuzSPJHarApnKo2u8VGBszuzHuo+zg1klZ5M3QwW3TsFN8r3kiVAAkUi4OtbNmylJaWRo0aNfI7jSECRHYh+TwM5XDrJX4MZcVVsWwga5VTOPkjNdhCAJsMmPSTZen5ySeCsmzBrK3Ygc0r3kiVeZsINgQxuCUyWvBGeUOGDKGRI0f6Ah52YJMBrzJnM\/JcwskfqcEWAthky7Y3p4OjW2wwZwtlkhwp3kjR1cS1XKLlEd00mcBUntdx4Mkhfw4Kq2hk165dSbyckFsoEcyiG80BP3v2bL\/baGUAWAEiVP5IDbYQwAalEAEDNxKHWMNF+GMXjZQDEJizucEbKVoImfQGz1RYgnnb4VryJB\/4tYpGytQXRutsYjSW12HU3nDxR2qwhQg2q3lEMJGjeHxH80deHVUNNg22sOBb80dqsCm5fNF6p3ZYUOJiJYnOH6ktm7ZsLsJJF2UlgYQH28CBA6lhw4aWWpKbm0sbNmygO+64gypUqGCa98SJEzR0aPgusNeqHVsS4GAz0zmsR+LP7RTxHSTHjx9nwNi+fbvbfVNyTV2vVBcY9RKQTxzIDQYI8ed2ijjY0CEADn9uJQB35syZGmxuCTTOyuFgmzJlCpUpUyagd3Fr2bwYx0j55F70RZfpvgQipR9RYdncFqedT+52fbq82JJApOb0cQk2L+eBsaVWurVmEkAwDm6kF4EQszrjEmxezAO12saXBLyal1lJKW7BFl+qoXsTDxLQYIuHUdR9iAkJaLDFxDDpRsaDBDTY4mEUdR9iQgIJD7bz58\/T0qVL2QI4QsLgCGnXrh317NmTcCCTp7y8PNqyZQvh0OOOHTtYvm7dulGnTp2oUKFCpoO9efNmmjNnDs2bN49kOoBgy3SqWZHsI9p69OhRFvn7\/PPPWdMffPBBtmuofPnyTrsS0\/kTGmwZGRlsWw5OKOOwad26dWnPnj2MoBSKMGvWLN8OA4CmX79+bP8mePp\/+OEHlq9JkyY0btw4KlzY\/05wAAn7OEeMGMHCy0b8HU7LDEbTIt1HrHniBlrIoHPnzgS5QK758+enhQsXUuXKlYPpVky+k9BgW7NmDePbgOUBTyJPe\/fuZQxYuBu8f\/\/+dObMGQY0WDOcTObAAgsV7g3Hb61bt\/a9j3W++fPn08qVKwmbp41OKmdmZjoqM1jtimQfs7KyGKfJuXPnGMBuueUW1g1Zvvny5Qu2ezH1XsKC7cqVK0wB4BqCQ1F08eB2gYUKX19w1u\/evZu5iwsWLKAWLVr4Bjg7O5sGDx7M3h0\/fjwlJyf76MXxTvv27Qn1gF1Ytmz44quWGaxGRbqP6enp1KVLFyZn8WPGqRgg34kTJ1q64cH2PRrfS1iwWQ0GvsSwaKmpqQxsn3zyCbN+nHtRnMfBqmEOx8EE6gHMT9q2bUsNGjRg\/zaigVu+fLlymV4oTrj6iPkwZFOqVCkvuhFTZWqwGQwX3EOADZN4WB+cbN66dSstWrTIL2iCV\/Hs008\/pcWLF1OlSpUCSjNjpgqlTDc0LFx93LdvH5vTrlq1ivCBOXv2LDVr1oxefPFFqlGjhhtdiZkyNNikoTp48CD16dOHgYrPM+wYgWVaNrFIK7CZEZ+qkK+GomHh6GNOTg6bDx84cIAFRWDZHn\/8cQL1OD5acMERIKlTp04oXYmpdzXYhOE6duwYm4PBFUSAo2rVquxpPIEtXH0U6f169+7N5JqUlMTkefr0aRahRMAJtOdyJDemEOSgsRps\/y+s\/fv3s2WAy5cvs7lU9erVfWK0AlsobqSZZbMr08H4+mUNZx8RhZ0wYQJzsTHXrVmzpl9bIGNEa83c72D7GM3vJTzY4OJgvQtf3ooVK9L06dOpXLlyfmMWSjDDDKihlOlUoSLVx2nTptHq1asDAktov9euslMZhSN\/woNNXFieNGlSQAAEgxBKmN4MbKGU6VQxItXHtWvX0qBBg9gOHfHGHoAfUdqNGzcmVKQyocHGAwWIimGdDKF+o8QXoDHJF3eLmC1qqwRIQinTCdgi2UfMzcD5f\/vtt\/vJDYv++L1evXo0evRotj6ZCClhwYZoGQC2bNkyxokPF1JOcCfbtGnDlOGzzz5j4erGjRuzqNqhQ4cst2vxsqzme8GWqaqY0dBHWDdEJevXr+8XjcSCthiEUu1TLOdLWLBhz2CPHj3YXkizJF5MYbRp+JlnnmGbkYsWLWpahhXYgi1TVeGipY\/ff\/89W0aBO5uSkkKtWrWivn376o3IqgOp82kJaAk4k0DCWjZnYtK5tQRCl4AGW+gy1CVoCShJQINNSUw6k5ZA6BL4XyJozHMHgjBxAAAAAElFTkSuQmCC","height":58,"width":97}}
%---
%[output:7596822a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANsAAACECAYAAAAZW15iAAAAAXNSR0IArs4c6QAAExNJREFUeF7tXX1MlVUYfzDETIJQWUp+ZCqoS4ukRMZ0STVm84ONovVhIl9SICU2PoZS0PzIATlluYaY5pqmwzZW\/WM0CS1TRqyPMdSoKSpqFCQZZbP9Tntv917ee3m9nHPvfe99zsase9\/znOf8zvO7zznPed5zAm7evHmTuDACjIByBAKYbMox5gYYAYEAk40NgRFwEwJMNjcBzc0wAkw2BTZw\/fp1KiwspPr6eov0xx57jCorKyk4ONjy2Y0bN6i8vJzef\/99y2cRERG0e\/dumj59uuUzLKvb29tpz5499MUXX1BnZycFBgbS\/fffT0uWLKFly5bR6NGjdXvS398v6qCNb775hn7\/\/XcaPnw4Pfjgg\/Tkk09SUlIS3XnnnTZ1u7u7KS0tjVpbW+mBBx6gXbt2DZBv3UdHzyiA1tQimWwKhk+PbHokunr1KmVmZgoSaMX+uWvXrlFVVZUgC8ipV1CnoqKC5s2bZ\/P1jz\/+SEVFRXTy5EmHvRwzZgxt2rSJEhISKCAgQDzHZFNgFLxmUwOqHtnQ0vbt24U30UpLSwutWLGC+vr6dMkGcqEO\/gYrU6ZMoXfffZemTp0qHr18+TKtWbOGvv7668GqCs9WXV1N8fHxTLZB0XL9AfZsrmPnsKY92eLi4uj48eP0wgsv0Pr168UUEKW2tpbefPNNmj17Nn377bfiM2vPdvToUcrIyBAeDYR47bXXBFnDwsLor7\/+osbGRnr99dfpwoULou6rr75KOTk59M8\/\/9C2bdsEgTSZ+fn59Pjjj4tpLORhWgpv+Pnnn4tnHn74YdqxYweFh4ezZ1NgExDJZFMArD3Z0tPTqaGhgUJCQoT3GTt2LFk\/g\/UR1kXWZJs8eTJt2LCBPvzwQ0HOnTt30qJFiwZoi\/VYTU0NPf3007RgwQJBynPnzlFqaiphGgnyQnZUVNSAupiiYpr58ccfi+80z8vTSAVGwWRTA6o92bR105dffkl79+6l6OhoCyH+\/vtveumll4SHw3RS82x33XUXgaTweLGxscJLwaMZKU1NTWJ6igIiFxQUWLypfX1r74n2oOuvv\/7KARIjQN\/iM+zZbhEwI4\/bkw0BiLNnzwoPVFJSQqtWrRKeDsaNKGV2djbl5uaK6aBGNrQD74TP4LXKysooKCjISPN04MABQRprb+WoYkdHB61cuVKQH5HNzZs3C6+rRSONNMjRSCMo8TTSGEq3+JQ92eDNenp6BKFg0Bs3bhSeClNDEC45OVkQ0BHZsrKyhHcyWqzJ9sEHHwjP6KhYTxmZbEYRdu059myu4ea0lj3ZYPAITMBThYaGirURwvlHjhwR\/x0ZGWnxYrI9G4Ig2EtzVE6fPm1pm8mmwBisRDLZFOCrRzaE5LU1GKKG8D7YcMYGNghoHdDAZyAnpnJtbW23vGbTpqhG1mzWz\/KaTYExMNnUgqpHtoceesgSXRw1apQIhiDkj3XcbbfdZlkjaZ5t0qRJYpvg0KFDIrjx9ttv0+LFiwco\/t133wkviS2BRx99VARRfvrpJyEP6zH7PTRrAfZ7cRyNVGsX7NkU4GtNNut9M21fTWtSC3wgBK8FJJztsxUXFwtSwethn+3EiRMi3evMmTNCJNZ++H+Qd8uWLZbtBGSJ5OXl2ezRIcr51ltvWbJLeJ9NgSHYiWSyKcDYEdmsQ\/JoVotMWgcprMl2Kxkk2IzGfhryJVE4g0TBwA5RJJNtiADqVXdENuvNZtRDlBIpUr29vSLzA2R0JTdSL78R8tEeophfffWVw15ybqQCA3AgksmmAGtHZMN0ce3atSIKed9994ngyMSJE22ySRxl\/X\/\/\/feCnEiv+uWXX8Q6DlkhiCBi+ugs6\/\/YsWMiE6W5uVnU1bL+ExMTdd8Y4AwSBUbBGSRqQGWpjIAeAko8G37BS0tL6ZFHHqGUlBRGnhFgBFR4NuTVIc\/v8OHD4j0pJhvbGSPwHwLSPBs2aJHUCoL9\/PPPQjiTjc2MEfgfAWlk0\/LxEN16\/vnnxYIcLy+yZ2NzYwQke7a6ujqRSPvss8+KiBfSj5B4y2RjU2MEJJPNGlAtuZXJxmbGCCiYRjLZ2KwYAecISFuzDYVs58+fF+dijBs3jseLEfA4AkgY0M6JkamMx8kGouEgmxkzZlhe5ZfZQX+RhR+rP\/74g+644w4lhuJPOCLDZsKECdJx9DjZkLeHoApSkWJiYvxlTKX3EyliXV1dYnZw++23S5fvLwKBI87zRBqdbBy9hmyDvb7vL4Ptaj\/\/\/PNPunjxIo0fP166kbiqkxnrqcSRyWZGi9DRWaWR+AhEhrqhEkcmm6Eh8P6HVBqJ9\/denoYqcWSyyRsnlyUhSIS\/oRQESLDWwHmTstcaQ9HLW+siAII\/+2I6st0KwFqAxF\/XbFo0FkcccHEfAriEZOvWrQMIx2Rz3xi4vSXtxwYDf88997i9fX9sED9suAtB7weeyebDFuHvnt0TQ+sMcyabJ0bETW0y2dwEtFUzTLZBjsl2\/5C4p0Umm3twtm6FycZkc3omv\/tN0ndbZLIx2ZhsbuI3k43JxmRjsqlFwN\/XLP7ef7XWpS+dPRt7NvZsbmIek43JxmRjsqlFwN+nUf7ef7XWxdNIGwT83dj8vf9MNjci4O\/G5u\/9d6OpWZriNRuv2XjN5ibmMdmYbEw2JptaBPx9GuXv\/VdrXRwg4QCJwQx0TxiiP7TJ00ieRvI00k1MZ7Ix2ZhsTDa1CPj7msXf+6\/WunjNxms2nTVbXl4e4RAaLuoR6OzsFEfe8xkk6rH2qhb4dC3PDAefrhUb6xnkPdyqjHMjcVDNb7\/9RmPHjpV+IYSH4VHSvMfPjcTl87W1teKSC1wwEBsbS\/n5+TRnzhwKCAhw2mntAkTcPmpfsrKyqKCgQLc+r1nk2JLKU6HkaGgOKSpxtJyIfO3aNSotLaXGxkbKzs6miIgIOnjwIJ06dYqqq6spPj7eKVpNTU2iXlJSEoWFhdk8O3PmTEpMTGSyKbQ3lUaiUG2vE60SRwvZ6uvrqbi4mHbs2EELFy4UIICARUVF1NPTIw61tCeRNVL79++nmpoa2r17t7hux2hhz2YUKefPqTQSORqaQ4pKHAXZMGUsLCwUU8fKykoKDg62INPQ0ECIlGFqGR0drYvYzZs3adOmTdTW1ibIGhISYhhZJpthqJw+qNJI5GhoDikqcRRku3LlCqWnp1NcXNyAtVV7ezulpqbSmjVrKCUlRRexvr4+QdZRo0ZRWVkZBQUFGUaWyWYYKiabHKg8hqMgmxbcyMnJoWeeecZGGe07rMUQLNErGlnnz59PoaGhIsjS29tLCQkJBJmzZs1y2EEmmxwLUvmLLEdDc0hRiaMN2XJzcwd4r+7ubkpLSxORSUcRRc37gXRLly6lxYsXi1swsYZDfWcBFiabHCNUaSRyNDSHFJU4SiFbS0sLrV69mtatW0fJycmWbYJLly4RCIwLwUE4vQCLRjas+ZYtW2aOEfFCLWEkwPvuu++mkSNHeqGG5lAJOMJpKLtTW5sq6nk2I9NIZzBiSllVVeUwwKKR7cUXX6QVK1aYY0S8UMv+\/n5hJOHh4TRixAgv1NAcKgFH\/EVFRUm\/VFJKgAQwIkiCX9Rhw4bZoHrgwAGxfeDoskONbIh2xsTEmGNEvFBLRJK7urpo3Lhx0o3EC7urTCXgiBtclXm2oYb+KyoqqK6uTuyxRUZGWoDAlgC2AkA4fDd9+vQBIPGaTY7dqFxryNHQHFJU4mizqV1SUmITzDC6qX306FHKyMgQGSSYigYGBgpkETjBZ3PnzqXy8nLd6Q2TTY4RqjQSORqaQ4pKHC1kQ14kXjtobW0lrJ+mTZumm66lecHm5maLt8IcF+syrM+0aGRHRwft27dPBEfg3aw9njXsTDY5RqjSSORoaA4pKnG0kA1Q2CciwyOtXbtWrKW0RGQ9sqEuCIeULxAOmSRjxoyh5cuXU2Zmpli0OypMNjlGqNJI5GhoDikqcbQhmyfgYLLJQV2lkcjR0BxSVOLIZDOHDQyqpUojGbRxH3pAJY5MNh8xFJVG4iMQGeqGShyZbIaGwPsfUmkk3t97eRqqxJHJJm+cPCpJpZF4tGNublwljkw2Nw+mquZUGokqnb1RrkocmWzeOOIu6KTSSFxQx7RVVOLIZDOtWdgqrtJIfAQiQ91QiSOTzdAQeP9DKo3E+3svT0OVODLZ5I2TRyWpNBKPdszNjavEkcnm5sFU1ZxKI1GlszfKVYkjk80bR9wFnVQaiQvqmLaKShyZbKY1Cw6QqBg6JpsKVH1Mpkoj8TGonHZHJY7s2XzEklQaiY9AZKgbKnFkshkaAu9\/SKWReH\/v5WmoEkcmm7xx8qgklUbi0Y65uXGVODLZ3DyYqppTaSSqdPZGuSpxZLJ544i7oJNKI3FBHdNWUYkjk820ZsGhfxVDx2RTgaqPyVRpJD4GFYf+HZ2Y7E8DPZS+MtmGgt7\/dVXiyNNIOWPkcSkqjcTjnXOjAipxZLK5cSBVNqXSSFTq7W2yVeLIZPO20XZRH5VG4qJKpqymEkdpZMMlGseOHaPt27fTyZMnxYnIq1atEtdA4fpfR4UPaZVjkyqNRI6G5pCiEkdpZGtqaqKXX36Z5s2bJ64K\/uGHH2jPnj20YMECeuONNyg4OFgXbSabHCNUaSRyNDSHFJU4SiEb7ggA0eDNcIOoRizcboOL7\/EZrv7VK0w2OUao0kjkaGgOKSpxlEI2EAbTxZ07d9KiRYssqOLKKVzMMXr0aCorK6OgoKABiDPZ5BihSiORo6E5pKjEUQrZ9u\/fL66Fsr\/wEOs4eDWs4Xbt2iVIZ1+YbHKMUKWRyNHQHFJU4iiFbFu2bKHjx49TTU3NgOuh8N0nn3xC7733Hk2ZMoXJpsjmVBqJIpW9UqxKHKWRDR5Kz3vhil9EKAe75jcvL08EV7i4hsCNGzfos88+oxkzZtDkyZNdE8K1SMMRF4Lee++9UhHxONnOnz8vbjw9ceKE1I6xMEbAVQTwo79161aaMGGCqyJ06ykn22DTSGgFwuGPCyPgDQiAZLKJhn5JIdtQAiTeAC7rwAi4AwEpZBtK6N8dneQ2GAFvQEAK2bRN7fHjx9tkixjZ1PYGEFgHRsAdCEghGxT99NNP6ZVXXqH4+Hh66qmn6MyZM4bStdzRSW6DEfAGBKSRTS8R+bnnnhPJyCEhId7QV9aBEfAoAtLI5tFecOOMgAkQYLKZYJBYRd9AgMnmG+PIvTABAh4jW2dnp0jjOnz4sIDpiSeeEJkkkyZNMgFs7lUR0d7a2lrau3cvXb9+nWJjYyk\/P5\/mzJlDAQEBTpU5ffo0paam0oULFwY8l5WVRQUFBe7tjAlaQ\/zh0KFDIt9348aNNHLkSClae4Rsly9fFu+5dXV1iX9hMHjRtLe3l9555x2KjIyU0jlfEILXlEpLS6mxsZGys7MpIiKCDh48SKdOnaLq6moR\/XVW8FIv6iUlJVFYWJjNozNnzqTExERfgElaH5AbiXxevBIGbDZv3mxesuFXA6\/joENIXI6KihJAXbp0iXJzcwXRNmzYQCNGjJAGoJkF1dfXU3FxscBs4cKFoisgYFFREfX09NC2bdsGkMi6v8juwdsYSASfOHGimaFQqjvssr29nSoqKujIkSOirSVLlpibbFevXqXMzEyaPXs2rV+\/ngIDAy0gYqqE8yPZMP6DBFPGwsJC8W9lZaXN0RINDQ2ENyUwtYyOjtY1RO19wra2NkFW3oJxzFdtuo1ZF471wDIHJw6Y2rPh1wNrCEwfU1JSbHqPKQ\/e+OYDW\/+D5cqVK5Senk5xcXED1lbOcNRA7evrE2TFgUuO3pRX6i5MJPzs2bNiBrBy5UqRhIzZBIqpyaa9mY1fZPv1hvYdXDnWGP5etF\/bnJwc8WtrXbTvgBOCJXpFI+v8+fMpNDRUBFmwLk5ISCDInDVrlr9DrNt\/bUbhM2TT816aAWHtZu\/1\/NEqnOHR3d1NaWlpIjLpKKKoeT+QbunSpeLQpYsXL4pfcNQ3EmDxR9yZbH446kMlW0tLC61evZrWrVtHycnJlm0CLRg1fPhwQTj7KKUfQm3TZb8gG08jbc3cGdmMTCOdkQZTyqqqKqcBFn8lnc+QjQMkxk14qAEStIQgCTZlhw0bZtMwtl6wfcDBqIHj4TNk49C\/cbINNfSPQFNdXZ3YSrFOFLDe63R0EJNxLX3vSZ8hmzbQH330kU22CG9q6xstNrVLSkpsghlGN7Xx8m5GRobIIEHQSdvTxOwCn82dO5fKy8s5gcAOep8hG\/p17tw5Ql4eOoV9JOwDcbqWPtmQF4mc0dbWVsLxatOmTdNN19IMpLm52XJsYH9\/v1iXYX2mRSM7Ojpo3759hOAINro5Nc6Hp5Fa1+wTkbH3g\/2iqVOn+t68ZIg9sk9EhkfCse4xMTGWCKMe2dAsCAfvCMIhkwT3MSxfvlxk8YSHhw9RM9+s7lOezTeHiHvFCDhHwCNZ\/zwojIA\/IsBk88dR5z57BAEmm0dg50b9EYF\/ATcMymV85P2XAAAAAElFTkSuQmCC","height":58,"width":97}}
%---
%[output:96ba5ace]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAANsAAACECAYAAAAZW15iAAAAAXNSR0IArs4c6QAAHY1JREFUeF7tXQu0jsXX3yQ5i+gUXQjnuFVouVb+vhOiCCGVEuUSud9VLknuqVxK7lHo+6KrLAurVgqhUEtyjVOsXCNObgu5fes3ffN+8847z\/PMc3nf9znnzKzVWjnvzH5m9uzf7D179szOc\/Xq1atkiuGA4UDcOZDHgC3uPDYfMBxgHDBgM4JgOJAgDhiwJYjR5jOGAwZsSZCBc+fO0eDBg2np0qWRrz\/44IM0adIkKlSoUORvly5dotGjR9MHH3wQ+Vvx4sXp\/fffp\/Lly0f+hm337t27af78+fTdd9\/RwYMHKV++fFS5cmVq1qwZtWjRgm688UblSC9cuMDa4Bs\/\/\/wznT59mq699lqqWrUqNW3alFq2bEnXX399VNsTJ05Qp06daMuWLVSlShWaO3duDH1xjFZ1ksD6pH7SgC0J7FeBTQWiv\/76i7p06cJAwItc78yZMzR58mQGFoBTVdBm4sSJdN9990X9\/Pvvv9OQIUNo06ZNlly46aab6LXXXqMGDRpQnjx5WD0DNm9CY8DmjW++WqnABoLvvPMO0ya8bN68mdq1a0dnz55Vgg3gQhv851TS09Np9uzZVLZsWVb16NGj1KdPH9q4caNTU6bZpk2bRhkZGQZsjtyyrmDA5oN5XpvKYKtduzatX7+enn32WXrllVeYCYjy3nvv0ZgxY+juu++mrVu3sr+Jmm316tX0\/PPPM40GQLz44osMrKmpqfTPP\/\/QmjVraMSIEXTo0CHWtn\/\/\/tSrVy+6fPkyvf322wxAnObAgQPpoYceYmYs6MEshTb89ttvWZ177rmHpk6dSsWKFTOazePEG7B5ZJyfZjLYOnfuTN988w0VLlyYaZ+iRYuSWAf7I+yLRLCVLl2ahg8fTh9\/\/DED58yZM6l+\/fox3cJ+bM6cOfTkk09SnTp1GCj3799PHTt2JJiRAC9o33HHHTFtYaLCzFy2bBn7jWteY0Z6m30DNm9889VKBhvfN33\/\/fe0YMECqlatWgQQFy9epB49ejANB3OSa7YbbriBAFJovFq1ajEtBY2mU9auXcvMUxQAedCgQRFtKrcXtSe+h75mZWUZB4kOo6U6BmwemOa3iQw2OCB+++03poGGDRtGzz33HNN0EG54Kbt37069e\/dm5iAHG\/oA7YS\/QWuNGjWK8ufPr9W1jz76iIFG1FZWDffu3UsdOnRg4Idnc\/z48Uzrcm+kzgeNN\/JfLhmw6UhLwHVksEGbnTx5kgEKAj1u3DimqWAaAnBPPPEEA6AV2Lp27cq0k24Rwfbhhx8yzWhVRJPRgE2Xw+p6Bmz++OeptQw2CDwcE9BURYoUYXsjuPO\/\/vpr9v8VKlSIaLGgNRucIDhLsyp79uyJfNuAzdN0RxoZsPnjn6fWKrDBJc\/3YPAaQvvgwBkH2ACg6NDA3wBOmHK7du1yvWfjJqrOnk2sa\/ZsnqbbgM0f2\/y1VoGtevXqEe9iwYIFmTMELn\/s46655prIHolrtlKlSrFjgk8\/\/ZQ5N9566y1q0qRJTMe2bdvGtCSOBB544AHmRNm3bx+jh\/2YfIYmEpDP4ow30t+8G83mj3+eWotgE8\/N+LkaJ8odH3DBc4eE3Tnb0KFDGaig9XDOtmHDBhbulZmZyUhi74d\/A7yvv\/565DgBUSJ9+\/aNOqODl\/ONN96IRJeYczZPUx3VyIDNPw9dU7ACm+iSB1HumRSdFCLY3ESQ4DAa52mIl0QxESSup813AwM23yx0T8AKbOJhM6jCS4kQqVOnTrHID4DRS2ykKr4R9PE9eDF\/+OEHy0GY2Ej382vVwoAtOF5qU7ICG8zFAQMGMC9kmTJlmHOkZMmSUdEkVlH\/27dvZ+BEeNXx48fZPg5RIfAgwny0i\/pft24di0T56aefWFse9f\/www8rbwyYCBLtqTZmpDdWmVaGA\/44YDSbP\/6Z1oYD2hwwYNNmlaloOOCPAwZs\/vhnWhsOaHPAgE2bVaai4YA\/Dhiw+eOfaW04oM0BAzZtVpmKhgP+OGDA5o9\/Oa511scj6dKxfZSvWBqlPvlqjhtfMgdkwJZM7ofs26dXzaNj0zpGepXaaoQBXIBzlO3BduDAAfrss8\/o8ccfp9tvv92SNYgjxJuIiHLnD+q45WNYaKDf8egLgAbA8VKgUj0qPuLfB3+sShD9iNd43M5vvOtne7Ahrq9NmzbkdOP4\/PnzdPjwYbrtttuoQIECnvgaFhrofDz6Ims2fMdJuwXRj3iNx9Mkx7GRAZsL5gYhWEHQiKdwHhrxAJ3fvirClevrdaBiPd+35FLYx+NieuNe1YDNBYuDEKwgaMQTbLJ2A9AAOKsS9vG4mN64VzVgc8HiIAQrCBrxBBtowyN5bscqSqlYz9FBkh3G42KK41o1V4MNzhX8p1vgDMD7+3iz0eu+Lwga3KEQhr5k1\/HAmWbnUNOVCTf1ci3YADI8142nA0zJfRxAkpE333wzoYDLtWDjXkwwvESJErlP2nLxiLHAIteBkwc7aBblerAlmuFBT6Ch554DusdF7inbtzBgc3gROGiGG3rJ54ABm8c50GWc7DXTbeexW6ZZiDmQrLk3ms1othDDIj5dM2DzyFddxiVTs1llGuVDFveN\/OWq1q1b01NPPeWRK7HN8CjrrFmzon5wm5BDbAy+L1y4kGW1mTJlSgxtXhcZeryOw4pvfmiiX7oyExjz\/4+Q0WwJ0GxcaMBzCGdKSkpkHuWJDxpsnD5PisG\/zb9z7Ngx9mRe+fLlXckWchGgAEgAMr4jJ7LnSTmQncct4HimHXlB4DT5o7NWT\/TZDcaAzdVU\/39lXcaFQbOpwCYDkec+C0KzcUAhJZQqpZTT71ZTgj6OHTuW2rdvz0BqBTa0t\/vNir4TSJ1+dxIlXZlxouP2d6PZkqzZONjw+CoAwQHQoEEDWrlyJW3ZsoXNqayZuCCLpqFcB9oByTDsNBfq4NoRkm6IGtdOkCDs8+fPp5dffpm1cQIb8spxja4yDWXtpQNQ5PeuWLGiMrWxEwgM2Jw4ZPG7LuPCqtnQf57kAiaR+Now35uozFBZIGXQ2pmuHlkdaQaAIhMO15Z2ZuRLL73EEnRAA6r6JJvNXrWtmzHpyowbmjp1jWbT1Gyvf7mX\/jhxnkrdWIAGNUrX4W2kjpODREyDayVs3HSaMGECof7gwYPZu\/9WGUfjKbQAV926dSMZS1XOFz548bl0nf2oXxNRZ2IM2HS4pKijyzg\/mm3hpiPUc+HOyNcHNUpzBTg7LcOFq0aNGlH5quU9m6y5uANB9fY\/OmoFNlFzqgChMk9lb+nIkSPp1VdfjeQPsNJs\/FugyZ0nHJhWebatwCZmQOX9Vi1S3Oy2y+OtKzMeRdKymdFsGpoNQAPgeHn6nltp2tN3ac+Fk0kHAVy6dCnbWyFrDHKxQWOJua5lsOHjsgCKwHP6Ju88BO+FF16I7OvQF3GPJWpU9Af1V69eHaVR7fZYfFEQAcuF3Q40Vk4d3gZ0Fy1aFMkxB56JC5T4u+yxDAXYkFYWk46kfEgfi40zzAVMBjJdiuXgwYNs87148WL254YNG7IoerkermCsWLGCpk+fTr\/++iulp6dT9+7d2Yb\/uuuui5C8evUqIZsKaG7atIkJHZK2t2vXjpCJ06roMi5IzQagAXC6xUnwRUcGB5uTZpO\/rXLl6zhIRLDhygnM09q1a0e56gEmFCwAsgnJNaHK9Y\/fnOZH5crXcZCIYELmHXFvKGp2edHS6ZPuvLqtF9FsEAikjQXY8KbH\/fffz97sQJrZixcv0syZM2MS6f3555\/Up08fypMnD\/NOIY\/YjBkzWMJ1FAAIaWiHDx9Ojz32GMHDtnHjRraKduvWjXD+wh\/fQe6xnj17Eq4+QNB27NjBaNapU4dgtiCbpqo4TSZv4wdsTKC+3EtrM\/+mjHI3uDIh0dYJbKJwoT5WaXlllzWMDi909m2yZrM63+ICPHDgQEKGU\/Fczkmz6XhExTo6+zYRbPCIYpGA7OoEluvKjFswOdWPgA25uTp06MCyXSK9LACEAg3WpUsXqlSpEksRmz9\/foLbFYOFHY4cYChHjhxh4AHQAC5oLZ67uVGjRizvGIAFACIX2Lhx42jevHlUrVo1ysrKYkDjifc4sGCuAMzwyqnyRbtZpfyCzYmRdr\/rHGpzz6O4p+KCw\/8GqwEudBQIF4p4SK4SeqtDYNFpY3Xgzfc\/3DUvu\/z5mJ0OtUFfPNYQFxIr3jgdxqNv8pGBbJ5aAS\/pYIMGmj17Nk2bNi1q1QI4IAgw7QCuK1euMPAhuTo0ofgsHMxPDJAn8QMg4T3Dv3l6WUwQNCZWbwAIGTUxeJiL0J7169ePyC1PDgibe9SoUQzoctFlXBjAhpVXVVThWvI5myxYKg+nnVOAayzx+6qzO7l\/\/DswLXmRo0HsvJFyaJXKQWPXD69hZlgYZNOS919XZvwssKq2jg4SJEKHpsJ+C2DDVfyOHTsyjSMzneeEhvDAHBRBKm5Sz549y1bmvHnzspV5yZIlTFvKh68y0FWhObqMSybYgp60eNEDEGRPI74ln6vF6\/te6aqcNnyRePrpp6McTW6sIa\/9sWrnCDaADFqocePGzBT48ccf2Z6O53sWCXPBnzhxIiFFLAAFTQhAiU4Ozgh4vaBN8d\/69evZ\/hAxb2LByrZ8+XJmcsK5kh01W9CTFk96ojME39E5G4tnf3Ro8z6KzhB5L6qSU539nc73devYgu3o0aNMgyFYlTs+7DSJuLFt3ry5cl\/BOyba+e+++64ykJWvqnYbbKPZdKdav55suvmNstf\/sveasnlqdf4YSs0GpwXMxzVr1rB9XEZGBuNEWMEGgWjRooXlbMGMhBPnlltuYfF80NDYJyZ6dfMuTqZlUBzgMgwfg3iWyenDD+H1iXq7Pio1GzyQQ4YMYUGwuKsE9zv3TtqBTWVG4uPytRKVGWl1TqNrRiICHeCxKjhDhIaGmQpPKcYGD6kBW1AinH3ocDmdNGkSC32TC54qTE1NDXxAMWDbuXMn9e3bly5fvsw8iXDNi2X37t2hdJBgD1mzZk1LBgHgOBe89dZb2ZuPYDgOzQ3YApep0BPkYLOSmYRoNgAJ0R2FCxcmODnKlCkTwzh4I43rP\/TyZDpowwHdfX7QTIxoNu4MgUaDei1ZsqTyW3DHw03\/xRdfREWL2B1qP\/LII5FoEbtDbWSYEaNFcsqhdtCTZuj540BSwQYA4EAZD5YiTKtq1aoxoylSpAiLLEF0x\/79+9npPUyzzp07M7e+VbgWPI2TJ0+mpk2bsuMAOCZU4VqIn+zXrx9zxLRq1YoyMzNDFa7lb3pN6zBxIKlg44fMy5Yts+SJHJ0gByIj4gFxc2XLlo2ioQpExqE4AOUUiNy2bVu2r4JZa1V0GWcOtcMk7snti67MBN1Lx0PtoD8YND1dxhmwBc357EtPV2aCHqEBm8Z9Nr9Md7qp7ccjqrpUaXeg6zQWMcwJgeGwQhDpIxe7OEynb1jFUvrhg9M3xd8N2NxwS6iry7hkajadqH+3giYCWG5rd03Gic1ifCTuiQFsqqfodO6cqb6lugTr1Kegf9eVmaC\/azRbAjUbJs\/qgF\/1m91kO10Mdfrdbg\/Mb2Lb3SvTuXOm+oYqjjFooXaiZ8DmxCGL33UZF1bNpvPcAYauuoZj93QAhBpHODg3dZMSS7yJrQM2BD6IIU+yiSibm3ZXXzyKgOtmujLjmrBDA6PZkqzZMPHiU3YqAZcj2OMlLACoeBPbDmzy7QDImfx+icp85hoXYXP8cqqfPaYXQMSLf059MWBLINisLo+Kq7\/d1RA+mV5NRCdhkG9iq5wvIg3xQqvVsw3yeMTHjfjTCqpHgZz66ud3AzaP3NNlnF8zEkndLx3bR\/mKpTkmdZeH4uUpO9W1fyewqW5jqwDBvYuqF5TxDX4p2E6z8W+prt\/I3lcnzWXFn3g5U3RlxqNIWjYzmk1Ds51eNY+OTesYYWJqqxGuAKfz4A9\/yo6v9vLeRwSNjrDIgqrSPKLZh8GJ7\/fj33Zgk99F4c+Q8+fQubZGvK34VJ6VJKq8mxgnMuWgiO9U+gWBDv\/8fkPV3oBNA2wAGgDHy\/X1OlCxnu9rz4cT2JzMQtmVr\/Nqlgw2eW\/IwcTf6cC\/xff7ncAmjwmaGDf45WMIHbOY7\/fka1bcWQPvaFpamutMOFYTZMCmLbrRFXUZ58eMlDUbgAbA6RYnsDmdWdk5GqySZshgc0oRpXpnREezcY+o1YIh7tH4u5Tgm3gEojIX5fM+eSHQ5b2qnq7M+PmG0WyHDxNuFvD7bKqV2IrB2LOd27GKUirWc2VCgp7OoTbf+6gEwUrouakp75s4DXzb6sAbv4l7KdXjq07eSNH05d8U+8L\/Jn5HVU\/lNEE9ft4H\/skmrh8gGLB55J4u4\/xoNo9dizRzG64lgoUTsXoHRPU0nApkVis8BB23PfDOi7wvsvNGqsK15H7DAYNAcrwJKkahyHRVtFQhXUG9haIrM37nXW5v9mwae7agmZ4MemKmUP59+VwtGf1SfVN18I2\/4WFfXGr2km1U\/I4Bm8eZ1mVcMjWbx6EF2kwlwHbJJwL9uEti6BeeNhT3dUGGeenKjMtuO1Y3mi2XaDZIgo755igxca4gvsCsenlZzLDjtSsGbB45p8u43K7ZPLI3RzbTlZmgB280Wy7SbEELT3alZ8DmceZ0GWc0m0cG58BmujIT9NCNZjOaLWiZCj09AzaPU6TLOKPZPDI4BzbTlZmgh57rNRtef0Z6K1NyDwfwMhxSUrt9isIvh3It2A4cOMAYvmHDBr88NO2zIQewwCJyBvGaiSq5FmxgMACH\/3QLTNG\/\/\/6bihYt6jnLSRA00N8g6Mg0zu9YTVkfj4hhR4FK9Si11atKNgXRj3iNx25eAbJEAg19ydVg0wUZryfv+9y250KFNMc8INoLjaDoqMbzR490dklWLHZXioLgSTzH45W\/8WhnwOaCq0EIVhA04iWc8lUizhq7K0VhHo+LqU1IVQM2F2wOQrCCoKECm5dnG+S+yJdkOWvsbqbHazwupiVS1WtfvPDOS\/8M2Fxwzetkip8IgoYMtos\/LPL0bIMu2HKyGen3yQsX4mP2bG6YFQRQgqAhg+303O6enm2Q+4IVPuuTWAdJTgab3ycv3MiP0WwuuBUEUIKg4aTZ+B7LyTwymo3YIiU+5uT2yQsX4mM0mxtmBQGUIGhY7dnEZxt0zCOj2YhEbY5nCktN3+tGJFzVzfWazWn1D3q\/FS+wybOuYx7lds2m8r66fabQDdpyNdh0Vv\/sCjYd8yi3g03lfXX7TKEBm4IDKo2is\/pnV7Ch304vgumCLaees6kcQmbPZrN86EZwq\/Yn2OOc374qQt2J0UGYgEHQUO3ZdFdY0WxOaT6IxGiWQyMeiOIHaCJUq\/iIby3JJ3s8fhZDebF1Gqsuj63qhcaMvHr1Kq1bt449qbZp0yZC1ks8g9auXTsqWLCg5ThVYFPtw0ShkM+lOHEnez0IwQqChhewcZ6ILzunNBtE5+\/vykLHrHjiJIDJGo9KIJz6IstFrgXb2rVr2fuCiMZu3bo17dixgz2HXadOHRo5ciQVKlRICTgwcNHMCZTRrDXd3f\/dGFcut8HFiTgxvnHMCg7iTva602TqrHxB0HALNqswrHx3ZNDlTgsY2OSzOnEs2T2CRLXQYEyyZYMxOy24OnMcas2WlZXFgAZthoc4ObDwIm6fPn3Y35o0aRIzBpWDA0G04urNTaFrK\/wXnflPJyr065dMsFQFrl8ALvVJ6wh30ezikwhautlt4g02lVa3CsNiPHioLxVvN94WbNn5UNtyoSmWFhNwrWM2Z3uwwRSEuThz5kyqX79+ZDxnzpyhAQMGsEc5R40aRfnz548aq8rBAaFXRUGwhmXuI8o6SJRlf60GqxuKHP1OqSUYYK3MLp1VMUiwHVowmFLO\/UUFbivPFgh5ww8T8NpiaXRu+yqlYLFB1nycivf9b7LS9nwhSalUT7mgBDkev7chdJxgOmBxsnB0aKjqhGLPtmjRIpo6dSrJSSKwj4NWwx5u7ty5MS\/hysIFoClB4pI7oBMDNE6jzH2UUrwcXfx1naurKLx5UMKpApbo7NEecurtlFKpLp1b+z9aTeQFJajxBEFHRcMqBM1psPEAXCjAhnfd8QLunDlzCOlfxYLfli9fTvPmzaP09HRbzebEwHj\/rjNBfoVKtf+I97hE+rLDxO94glyE0BdZ29ua0A6M07FU3PA+NGCTc3PxQdjlLvPDSDdM0q17styDdKHFaNvqly5dopUrV9Kdd95JpUuX1iXN6hXJXEnXLRnmqk08Kp+t3ZnO\/qfTv1aEj\/GIfQuCjswf9PPS0X1UJPNrT2zQWTzdEDZgc8Mth7pfHi9Eb+6L1swBkqcX045Ro5vOBEnSE614j9NTp4hi+IN+onjlWa4Dm50ZaeVp8jpZfttdaDGGTpZr4JeMZfuwaLZ4j9MrA2X+oJ8oXq0BpyAHt\/0MhWbz6iDBYHlIktuBB13fS5JEL30Iarw\/XVeFSt1YgIod3sA8lheld0cOX3ML3Xb5z5guJmqcXngjyoPYTzueJXKcoQCbV9e\/1wkx7QwHksGBUICNH2rj\/EqMFnE61E4Gw8w3DQe8ciAUYEPnV6xYQf369aOMjAxq1aoVZWZmaoVreR24aWc4kGgOhAZsqkDktm3bsmDkwoULJ5ov5nuGA4FzIDRgC3xkhqDhQMg4YMAWsgkx3cm5HDBgy7lza0YWMg4YsLmYkLNnz9KCBQtYqiGkHcKVoEcffZS6dOkSFdPp9SIs7vQhIHv69OkxQddeadoNL5njQb\/++OMPlknmq6++Yt1s2LAhyyxUqlQpF7OSfaoasGnO1bFjxwi53Pbs2UNt2rSh6tWr09atW5nHFMIxZcoUKlGiBKPm9iIsgIR4ySFDhrDMKqobDm5pOg0r2ePB2Wrv3r3ZeNu3b0\/gAXiYN29emj17NpUtW9ZpCNnudwM2zSlbunQpDR06lGmeunXrRlpt27aNOnXqRM888wz16tWLpZRycxEWKatmzJhBn3zyCQvqrVKlSgzYvF6utRtaMsdz+vRpGjhwIJ06dYoB7Oabb2ZdlXmZJ08ezdnJHtUM2DTm6cqVK0wo8EbKrFmzokw8mGKDBw9mK\/L48eNpy5Yt2hdhT5w4wYCKNi1btiR8Z9++fTFgCzrCJtnj2bx5M3Xo0IHxVFy4zp07F+HlmDFjbN+e0Zi20FUxYPM5JVidodFSU1MZ2JYsWaJ9ERZgw56lRYsWdO+997L\/V1018hM76nZ4iRoP9r4wlxE1lFuKAZvPmUZIGcCGjT2edvB6ERbdQFsV2PzQdDu8RI1n+\/btLDTv888\/JywmJ0+epHr16rEooooVK7rtdraob8DmY5p2795N3bt3Z55IvvewAgw+Y3cR1glsXi7Xuh1aIsZz4cIFtvfdtWsXc4pAsyE87\/jx4+ymPt6dgYOkWrVqbrsf+voGbB6naP\/+\/ewxIpiCcHBUqFCBUcquYEvUePi+DA6abt26MR7my5eP8e7IkSPMQ4kjlQkTJlg+X+hxypLezIDNwxTs3LmTHQNcvnyZ7c\/uuuuuCBU7sNldhPWq2Zxo6gwvkeOBx3X06NHsXRk88FS5cuWoLoKf8Myq3pzRGUuY6xiwuZgdmD0478JqnJaWRpMmTaKSJUtGUfDjzLACqh+adsNL1ngmTpxIixcvjnlNTcfUdjFdoatqwOZiSsSD5bFjx8a8BAZSftz0VmDzQ9NueMkaz7Jly6h\/\/\/4sGqdWrVqRLgL88MiuWrUqR3oqDdg0wcadB\/CU4cFYuPpVxc9FWCuw+aFpNbxkjgd7s65du1K5cuWiLgvjgB9\/r1GjBg0bNizmUV7NqQptNQM2jamBBw0AW7hwITVr1oyZkHKBOdm8eXMmIF4vwtrt97zSVA0vDOOBdoNXsmbNmlHeSAQHiA4njenJNlUM2DSmCnGEnTt3ZrGQVgUgxKF2SkoKc2nLGXl0LsLagc0rTVV\/wzKeX375hR2ZwJwtUKAANW7cmHr06GECkTVk0lQxHDAcsOGA0WxGPAwHEsQBA7YEMdp8xnDAgM3IgOFAgjjwv+uBIGTUiHhFAAAAAElFTkSuQmCC","height":58,"width":97}}
%---
%[output:69c2359d]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : Ba1_A\nTaille de la serie    : 240\nStatistique T_max     : 10.0881\np-valeur (bootstrap)  : 0.1040\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:7178928f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:73fb85a8]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : Ba2_A\nTaille de la serie    : 240\nStatistique T_max     : 9.5956\np-valeur (bootstrap)  : 0.1030\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:86992b43]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Ignoring extra legend entries."}}
%---
%[output:3da873bd]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : Ba3_A\nTaille de la serie    : 240\nStatistique T_max     : 9.8589\np-valeur (bootstrap)  : 0.1020\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba4_A\nTaille de la serie    : 240\nStatistique T_max     : 9.6456\np-valeur (bootstrap)  : 0.1240\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba5_A\nTaille de la serie    : 240\nStatistique T_max     : 10.6189\np-valeur (bootstrap)  : 0.0800\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba6_A\nTaille de la serie    : 240\nStatistique T_max     : 9.4663\np-valeur (bootstrap)  : 0.1060\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n\n=== Resultats du SNHT ===\nVariable analysee     : Ba7_A\nTaille de la serie    : 240\nStatistique T_max     : 9.7270\np-valeur (bootstrap)  : 0.1110\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:937c0a96]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAALIAAABrCAYAAADTu9iPAAAAAXNSR0IArs4c6QAAF7tJREFUeF7tXQlwVsWWPiFGlgHlgfLAEEjK9ck4xqVKxKGSKMYZXKO8517oKKXUIC64I0JE3FCfqOVS6CDlbqHjUr43PkTCYPmkygV1XEZ08iNGAQt86qjRBDLzdTh\/Tjp9l779\/\/mXdFel4L+3l9Onv3vuud2nvy7p6OjoIJ+8BgpcAyUeyAU+gl58pQEPZA+EotCAB3JRDKPvhAeyx0BRaMADuSiG0XfCA9ljoCg0UJBAfvPNN+mMM85QA3DzzTfTqaeeqv7\/9NNP0zXXXKP+\/8QTT9D48ePV\/2+99VZ68MEH0wMmy+DiunXr6Nxzz6WvvvpK5dljjz1oyZIltPfee6frveeeeyKvFQUiCrQTBQ\/k448\/nm655RYaOHBgN8ACyAceeCBdffXV9MUXX9DDDz9Mw4YNS4MW5a666irih8L0QPDDgAfEAzm\/EV7QQJ44cSJ9\/\/33CqRI5513Hu2yyy60evVqZZGRYLmldZbD8fPPPyugI\/HDgP\/r11988UUP5PzGcWHOI7MVBUibmpro9ttvV2q+\/PLLFZjnz5+vwLtq1Splcdka62PBLgVbZ3kf7giXXb58uQeyB3LmNcBAnjNnDq1du1b5tJWVlfTUU0\/RtGnTaMaMGVZAvuiii9J+Nksr3Yl33nkn7XvL3ui+dOZ76muMq4GCdi3g16ZSKfrkk09Uf\/fbbz+qqalJuxPeIseFQeHnK3ggjx07ttsMBv\/2PnLhg9OmBwUP5IMPPlhNnSFhymzLli1pi2wza3HBBReoWQwknq7zsxY2UMpt3oIH8gknnNBtiu3TTz\/tMVOhzyNL0EL9fh45tyDMROsFCeRMdNzXUVwa8EAurvHss73xQO6zQ19cHXcG8qsfb6E\/LH5faeWmE\/eiC2sqiktDedCb\/33zWdp8xxQlyW\/+0Ei\/+f31eSBVfongBORNP\/xCU5f8Fy2csq\/q1RXL\/puWnvv39Nsh\/amujiiVIiqv6KDlK7bT5PrS7r\/LJlGKUlS+bQwt3\/4aTdZ\/a\/kH1k9SFbZXjKGj\/7NfumzZpFWUWkJElUTl2zqorPTI7vVq9bioX++DTV02Zcv69aOdSktU9du2ttBX8\/+JRkx\/SP3efP\/5tMec\/6DSYeWBzX\/55ZeEv0JLo0ePJvwlSU5AhjWe\/cJn9OK\/VtPg\/qX0L0s\/pGn\/OJom\/W44VVXVqcUKJIC5ZUPnwKR\/l67v+r1tLLWI35Ule1Jbezm1bPginb9sQ1f+VKXoKn7I39TZpmoH9WrtJlESl3Gpy6ZsWWmJWqlcuXIlwRpvXXoZlS94g0oGDKGvF0ymXY+\/lAaPP8XYFQD4iiuuoDVr1rh0NSdlDzvsMFq4cGEiMDsDefHrX9K\/TR2nOg4g1+0zTLkXVVVVaSDbamWnnXZSRdrb222LFk1+ALm5uVkB+buX\/kijZv9J9Q1AHlRdH+he8PI9AFFeHmy1801RePAWLVoUGOAVJW\/WgFxX12WRo4Qw3cdAskVPUr4YyrgA2RTxB30iyApp\/fr1NHfu3LxREz+AQZGKUYI6AznItYhq2N+PpwFb1yIIEI2NjTRv3rwejeJaPgA6p0AO+9iLN0w+V5QGbD\/2TIB45JFH0sv455xzjvK\/YZnZOmNpH9eR9FVOXNNXQqXMW7duVRGIRx55JL333ns0aNAgFSYAi4+dO6+99hpVV1erIpzP1OecAhkCyem3Z6b9g\/rQ8ymzGpDTbyNmLQv80EOrOiDgTuB7BUm3vgAyXEAkJpwCkF955RUVCouETQYLFiygqVOnprd+yd6hPQbtvffeS8ccc0w6H8refffdKrQWoQOcLy+B7Dpkcp8dgnw4CD5oXx4rFzs7Tj\/99PS+vLD8tjLayMS7SV566SXVjLRemZQpbh90IDNYa2tr1SyIngBk5GGrHAZkDsjifgKg2MgAS4zYcGxowBhOnz6dvv32W9UU9lDW19fT119\/TR988EHgx1zOLXJcBZvyQWk33XQT3XHHHWo\/HYJ7sAH0yiuvpOuuu46uvfZaVUzmwauMlccfBrg2a9YsY35b+Wxlwu4RWD1EzrFsp512Gh199NEZk8mmDzog2DcO8oX1+ybXAnHf6A\/ywp\/GWOFhR0LYrG6R2bVAQFefscj6a+rJJ5+kk08+me666y5lnbGplK0vnnYoZsqUKWpbE8CDndIYPDwEen7eRW0DBD0v6g6TSW8DcsAHxQBnS6aw\/uhAZv8YPjCsbpBFhrWG1dYtMucP8p2xkcEDWdOqBAHAgw2hSADyhAkT0tuR2PJJIIfldwFyXJnQhnwzwCplSyYbICNvSUnnYhSDlcszyHnOGteDgIy+SYvMdZh85D5tkfGq4lc0W8FcA9lGJvaV2W+P6oPLw2ULZDn1BtDCOkPXADKSdDuCgIx80ufHb7gcvCOHt53hW+HSSy9Vi1kyVhx+84UXXtiNh0R\/G4fteI\/Sl9M8clTlce+z1WOilShXwWSRM\/0at5FJWmImdYnqQ1zd2OYL+miSU3BcJ4O6z88j2yrZlB+AgZ8lfc2ojzcdyFH5beW0kQl1m165mZYpbh\/Cvv55ZQ8+LSwpfGKAOR9SQc9a6K8qKJSZgzClw7RY+rKlDmT9tZd0mdP0+oySSd9Gxa9cvF1k\/1xksgGaKyBs2spkXle588K1yKRC+npdroDIlf5c5fZAztXIZaldEyBksJCpWd3FkG8ZSUKD60i825xnMjDv\/9hjj6kVPMwxI+GjcenSpbTvvvv2iOUwEdt4IGcJEIVabVSshalfMtYC5UFsw2BlQM6ePVvN069YsSI9fx8HyCiHtQDkXbx4Mc2cOVP91pMHcqEiLktyhwEZ88jy4w6WGkvUetAQLCkDUIqJWIpDDz1UzY\/zTAc+dMMsckEDWcYf5IofzfUJzxLOsl6tqd+whFjuP\/HEExUZI1tEDiiSQMbYAZywyps2berGFc1BQVjwwH24Egzk2267jTjehDspKX+lRcZ9nQXVdbyy4iOHsVxmfSR3NOCqmN6SM9PtuAJZspDC35Ureph54eg25DvggANUpJytRc4KkH\/66acOPIUIxMFcLnwjPaIriBVeRqvxgMigHlzDU4nYCdBaITIMzPEcISY\/KvjplZ0cMmSIipaCVYdFQfAQgoqC4mNl1Bq3JQOLONgIbQQx2utWJE4ZqQe9\/1J3Qff4wUewEXxQTD2GxQDz+EgCc9Y\/Axl6vu+++1Rde+65J73\/\/vtqLNgiS12NGzeOnnnmGQIPNDP+Q+f333+\/Gi+UHTNmjIpeQ4AXLDLn23XXXWnZsmWKCRX7BXfeeWficRs+fDi98MILavw+\/\/xzNZ36zTffqHUDJHBZMy+1q+EpMQFZPpWSggqCAZCgYUU0lAS\/tCy6RWYhdX41XJck3XiQ8DGA187bb7+dDnJBm7vvvrvKi2gzKFGfl9XblPxtzAEHGaE4SdwtueJk\/3h5NagMx1Fg8CAfcywH8SqjbqkvOXDc7iGHHNJDPl4plPqNA2TIwzp8\/vnnVbDVEUccoZajwV562WWX0Z133kmnnHIKHX744d3kf+6552jz5s0qRgNxyrvttpvSN6zvyJEjlQXGb5CsYzm6oaFBuSvQBcAOVwWLLpdccokyjGeffbZyQZ599ll6\/PHHacCAAUpn3F+UzTiQOcaUrbMMTWTCQFhF\/RyOOEDWrSO3gbLMRwwr8NBDnVvfATr9bRHUYf14BBNQmAeZQc8E4aaFF85jKmM6m0R\/y+hvK11uyYp\/\/vnnqxhe\/TiIoIUgWFmZ5HcI2kG8x1577UWTJk1SQGIfGWXw4A4dOlSBCT7w9u3b6cwzzyQ8TNA3AM8GhiMPeSz4eAu+L0Nvr7\/+errhhhvSx1ygDD+4OqZMJwVkHMjsBvDBMKwwkzuAeyZAB1lkHhgdJNkGMupnsEoAQHa2vPJDRVLS6j4syiDhrQCd6IOkhzsyoHXWexOQ+aGJGtQoi6wDGRYYFAFI8HFhWT\/88MNuXWM5AfogILPbp3\/UoSw+BPHBZzJAOQGy3qg+kLovZpqViAIyW\/lMW2TpcoRZ5KA+sd+IgUFg\/8UXX6zcKA5mQjndmuhvDK5bfmcA8JKAHP3uDSDvv\/\/+3UJhGWTS6vICBsstXaM4FlnvbxiQ4WfLiEZp6aMe3qAx4+slW7Zs6WA\/WFpdfiqlj8yWjS2r\/oXLlUYBGflkWfkaYh\/ZpBB5ClOmfGQcq8AnNoHlhqeF+FXJcki\/Gj4k64d9duhOl539Z9M3RdiD5jKoXBb+J4fFsr732Wcf9Rb57LPP6MYbb1QPa0VFhbKkAD0sOXz\/jz76SPnB\/fv3Tx\/rhnL4kMM9U1n42nwEHPL+8ssvPep94403unFQszyot6WlRb01ksaklFRVVXXAN8MXsO2sBRRkajgOkFm5svNyxsQWyMjPg4j\/J521kG8YfZZBd43gfiE\/PkTZ2rCFZj82ziyI7mq5AFkyDeFDjT\/Y0pZrR5A9fON+\/fqpy5xP\/z\/ycFA+\/4s8YWVxP6jNMHlQLmdMQ1Hm3t\/PjQY891tu9O5b9Rpw1oDzyp6nlXUeA19BBjTgBGTPNJSBEYhRhedHjlaSE5DDaGXrqE7xFI9sq6A\/t63oIcnQY+sV33FreQW1\/mUFNQwqU4SwI9u205\/btvXI33BsvfoKH1leQX8tK02X\/efXS2lJ3RKqTFVSa\/l2GvD6UUSUota2Cmo1tButkuzkGHpsGcRSMrb+pWf\/uNWGQfVKb5VUSStppedHjjkczkAOpJWt66KVBZj1NKBlQ\/oSwLwR4NyRRnf0o+G\/HUWbWjr5kZE2ivw96JElYXJlFz8ywJwvaUBLV\/9ay4OBvLGsUy9qi\/7KTlpZz48cPYrZA7LnR47WfkgOz49spz5nIAfRynp+ZLuBMOXOND+yu0TZq8Fl7hxSOQHZf+xlb2C5ZlvXwhUQaNdEjxVnxS2KuTNMW65yOwEZgnla2eyCORP8yLYS6mxDct+eab8d6udVUMQbI4zTFH6a10C2VZLPb68BF35k+9Z68r8ByPibPHlyD2vNlhqE3tgxjTDcIC7logayDRcxK0LnWcN1GWcRFisdZ2BtZMp3fuQ4\/dXzhLHWI4YEwT4cDC93XJtci5dfflkF8iOSkA85MsmUc9ciiaK4jC0XMe8h8\/zIwVqPAgT46TFBiSnMnrTfnfWGkX1zhKCJ2FwHMn4jUhDRbRMnTlS7S4JSlNxROHP2kaMasLmPznh+ZBuN9cwbBQgcwsBAbg5oygRkJuyWG1DRVphF\/vjjjxVT\/ahRo2jDhg2Kub5PADkuF7HnR3a3yKghDMjYUyd3CbEvLN0uhKgijJVJV6RFxlarRx99VIW44oAcbGBFnYMHDzYKH\/UARj3eeWORbbiIewvINjLlMz9yFAiycR9niMD9wBYz+MbwleFaIOC\/6Hxk7pANFzHKeH7k5BY5G6DNRJ0Fb5FtuIh5f5nnR\/ZA1jWQU9dCTpmxYJ4f2c2+uVo2t9aTl3aVO6dATt5tXzJbX\/+50qwHcq40n6ftugIiV91yldtb5FyNXJbadQUExEoSNKTz7jG\/ctxuusrtgRxX0wWSzxUQDGR5HnVU0JB+3\/QBH6U+V7kVkPUDXVxjFaKElvdN9Flh\/h+fshlEDmPTdjHmdQVEEJCjgoakLmGdcWqUzcmzrnKXzJw5s0PSk\/Y2t3FcIMfNV4zgtOlTFCB4LyXvCTTVnTRoCHXJRST87rWgITANBQVNh\/GcMUccuHV\/+OEHtYoDzjRcnzFjhtIPW3ZduZI5E\/mYqhYca7oS+UBvyVMMFiEkpq2SxHt84LdOg8VBLkFvG71d1kmY7Ey3Bdnmz5+vZJozZw6tXbs2rQ9Q4er8ajbAtM0bBeQqqkpvbm0OWKROGjSkL2z1atDQSSed1BGk7DhAhqJRnjnieB4YQSYAD4KsJQcxXjdBQGZmTD4GF4rhOnTA6xzOOGQdbUn+Np0mVfLYydee3k9J8qeXkbIDyCBO1Fk98RAw5zFTxdoCMmn+KCCzRUb9NkAOCxpizjseN5a9V4OGXIHMfHH6q18OeFwgM+ultI7MxRYGZNyDVQTjO3M4m4gDJdeziWETD43kamOrL89I1oHMlh87IoKIGW2\/4JOC2CRvkrrCZi1MQUPAAN6qMoHQG\/X0WtBQPgGZWfB5u4xkygwDsiScrq6uVm8IWGcuE8T1LBWvrzLqK4wy+ovBK+UrJiAnAb9epteDhlx95ExaZD5p3gSaKCADiHAvTMyYkoc5apAkAydchjCZPJCjtBn\/fpRLFFVT5KyF6TyRMKZ2JsaO8oP1M0JQzuQWsHWOArIOQHYdgs5D0Q9xjzrfA\/KxD8+yeyBHwSv+fWcgm+aR5aE10l+SPMpBRw6YgMyvXXAhw+fFkbFNTU3dXv9cTs5pY3vM6tWrFQczH2jDfuyIESPorbfeUm4EzwrIj0PexRt2QpVUs+4X6gf3mGT3QI4P1KicGQFyVCP+fuFowBUQueqpq9x+iTpXI5eldl0BAbFcYy2SrAy7yu0MZPm1n6QDWRrPPlutKyAYyDaxFvJ4Xj4CGGdV2ywEucrtBGR0YNasWepUUiScTIqTMW060GcRl6WOuwIiCMhxYy3kkb82OHCV2wnIaBwfWPjg4qOs5OpOa2urOt0nabrr398lWl9JNDZFVEt0aut4GjUqeX0mOXByEQ5PdEn5dGYHn44EQhQcLpMkoQ6simIhixM+9DG2OIoXJ0HxCU8gbMF1pJUrV6rDKTERABJLJGAEsTw49am0tItaV5eL5V64cKHKa5ucgcxH2aJhHO01YcIEdS4dBnfRokUK4DgCFqmsrIza2trSMiJwBaTW6CD\/4WZzcw399ckLOgkYRALV6ifjW+l3AweCQFj9Vf\/tb+nT7jnrd98dRH837htqW7cuVB+\/\/vqrkhd\/SZM8RSlpHb5clwaSnuyUVSDjyY2TqqpAG9I9gdhDp1hqb28n\/pO59fLNzf9DHaAapZI4zavjuZImfiXCkpSXlyetxpcjojVr1ijjF4f5U1eYM5DDXAtY1KZUqtOyrsc\/KaqsIVqaIqpc32lux44lqq2tpI0b7V\/vCtSj22mAVra9fbSihMJ9pNSqHW1PTVFlqnaHDjr5dlKVKaqt5Gv2eHL17exbLN4SLrp0ArL\/2OsiT0xiRQoVknKWQtLMJgmolzrIGZAhhJx+60uDyQPgovxCBXKQ3AUN5GIbDNv+FDqQIT84jTkUAFvJEGONxOsCMmyA42wwO4E4ZMxaIMyVQw+OO+44Qh281YmD9BHLgskAhBjIY5LzxiKHPZ2sEMRJcExE2AJKJniPM8VtLAcvSPHyjcRvI\/VJoM22BOmo1uCaNzXFf5SSlN8x2ZNuBOPBjJqmTaQHHXQQvfvuu2miQhRk1+Kss85SB6dj8UNOv5qAjKhEJMxoBc01uxgFJx85SOW54j22bXf58uXq7D4Evsuge32nShi0dOU3NhLNmxcPjKbJkh1nnseqIEl5yDZ3blf1kJ839Oox2WyV8S8bJWmRwWAPwzF79mwFZHYtTED+8ccfCYFXnEzGIe+ArI8CBMwF73FUu\/ouX95zhgB\/XrGMOgtDV34hWmQGctS2f6aNbWhooFdffZWkRebNDTU1NerBYOvLFh9rAHwt6CnNeyDnivc4brv8umTw8tYsVrgM6TQ9pHIrVCxTmkeZpEWGWNI1w+8HHnhA+dDwhZGCfGTcQ9jt9OnT0\/sVsTMHbuVRRx2ltqKxj4y80t1kdeQ1kG04hjPJe2zTrsk3Z+XyPV6xLDYg59EzlZ4BSzL7lRUfmZWTK95jm3blXHiQG4H6kEybSF2sSD6BKB9kcdFl1oCcK95jm3YxeI2NjeqrW0ZqyS95\/S3hLXL2IJ93QDZ9\/fYG77FtuzpVGIYoaO40aEu\/i\/KzB4nCrNlFl1mzyIWpSnupXZRv31pxl3DRpQeyIzZclO\/YdM6KF2WsRc60mScN90UgB6nex1rkCSiTiFHoQPaxFklGvQjL9FjZ23GyaJyumqKgLUItsPurR4oqj\/W1zjW2zuRjLeKMVB\/I0yPW4v\/7HDPUQu1i0VO8PS2dpZKUh2wi1EIB2cda9AGgRnWxGCyyj7WIGuU+cL8YfGQGMoarUGMt\/g9moPVkkNsVXQAAAABJRU5ErkJggg==","height":58,"width":97}}
%---
%[output:8e8d30b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:57028429]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:11874e63]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9b14cb52]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:00570387]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9c491c5c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7c87fdab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:53da11d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:60fab80a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:133f7c9a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9d1ea025]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:076d2b23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2f137536]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6908aff7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0d9b3f26]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2cad31eb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2fb9451e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:02dcab76]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:20bfb149]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:692db228]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1801bd33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5f4d4999]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:249f9816]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:863778f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9ee48031]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:57d3e313]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:15451798]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:979bdba9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:7d96775e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3c0673da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:837ca3f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:01d970f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:044fa5ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1ed5920c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:11b46dbe]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:674f9d41]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:383f4db8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:563da939]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5a5a0a8c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:03fdd8b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:89247c3f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2093555c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:92fd0be8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:96788f96]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:50ad103c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:04143249]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:92dd5f37]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:467cd496]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6657d4be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:10d57844]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6da3ae6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:93851e33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:220dff2d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:41602578]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6af10b35]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2991f550]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:333b04cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4983c5b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:28f1614a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6ddd829a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3f7884f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:41d4fa27]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0d9fefb2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1fea820b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9fb6bbf0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3d998d9c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:29d1b1b7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:04d483b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:076161a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:77adc073]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:16996c90]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5e1c8e91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:854d5457]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3840bef9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0e108e19]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:712df5c4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0fea583f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0df4f7b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9f1faa57]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:91c96e39]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:17bd0fda]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:329b88b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:504b269e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:51d1cdc4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1617207d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2c9fbdab]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2fa1a655]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:90647205]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5bedff53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:26a0dc79]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:206dfb78]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:677d5102]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2f7b30d9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0ae922b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:016035ee]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:817310d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:79004a0b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5c48ecde]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9916eabd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8150bc17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1a440d44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:94d4fefc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4804dee6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:62fe1f9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:44b6d2f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1ba547e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0c4761c1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0be1ba44]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:111d9300]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:085d80fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3d93dd73]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:78fb2c1d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2cedc9b5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:0827a9fa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:1b6875ae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:88007c1c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:931eadd5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:61c5ad7b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:31339437]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3905f5f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6bacc0e9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2cad2eaf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3bb615ef]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5fe998d6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6c0c3242]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:97cc1070]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:92811090]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:558f4e50]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8cfc6482]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:6a9ae863]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:763f99d1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:913a087c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:861546cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:4196beb6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:9a3f74d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:5df0932f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:2d3c4887]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:77b5f251]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:3e5144bd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:42c6d1fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Iteration limit reached."}}
%---
%[output:8a71b0ce]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAIcAAABRCAYAAAD4gLU4AAAAAXNSR0IArs4c6QAAD3VJREFUeF7tnQmMVEUax79RQFjOQXFQGUAEVBZWxUxEghyCgdVVJLqgRgmoqIBiUBBRmR0ERS5RFE9QVCJC1A24mLjegngQFw8EV8Ah4IGgDAiEWza\/wq+3pua9fq+nj+mGV8kk0911var\/+676vq\/yDh48eFCiEq2AxwrkReCIcOG3AhE4Imz4rkAEjggchxc4du3aJXfeeae89tprsQfr0aOHPPjgg1KnTp3Yd\/v375dx48bJCy+8EPvuxBNPlGeffVZatWoV+w6x69tvv5XnnntOFi9eLD\/88INUq1ZN2rZtKxdffLH07t1bGjZs6LmIe\/bsMW0Y4\/PPP5ft27dL9erV5cwzz5SLLrpI+vTpI3Xr1i3XdsuWLXLdddfJF198IWeccYbMmjWrQv\/2M\/rVSTeuc5JyeIHDa9N\/+eUXueGGG8ymaXHr7dixQ6ZNm2Y2FzB5FdpMnTpVzjnnnHI\/f\/fddzJ69GhZtmyZ7z4de+yxMmHCBOnevbvk5eWZehE40ghrL3Aw3COPPGLeVi3Lly+X\/v37y86dOz3BARhow19QOfnkk+Wpp56SU045xVTdtGmTDBs2TD799NOgpoZyzJgxQzp16hSBI3C1kqzggqNjx46ydOlSueaaa2TMmDGGJVCeeeYZGT9+vLRr106++uor851NOd5\/\/30ZNGiQoRhs4MiRIw248vPzZe\/evfLBBx9ISUmJ\/Pjjj6bt8OHD5eabb5YDBw7Iww8\/bDZc+7z99tvlggsuMGyN\/mBTUJt3333X1CkqKpJHH31UGjVqFFGOJPc\/bnMXHNdff7288847Uq9ePfN2H3fccWLXgb\/D121wNGvWTIqLi2X+\/PkGTE888YScf\/75FcZFnpg5c6b07dtXOnfubEC0YcMGGThwoMBWABt9n3rqqRXawrJgO4sWLTK\/KWWL2Eoa0eGCQ\/n+Rx99JM8\/\/7ycddZZsQ3ct2+fDBkyxFAQ2ItSjgYNGgiggqJ06NDBUAEoRpiyZMkSw64oAG\/UqFExauW2t6kT4zHXsrKySCANs9CVqeOCA4Fv7dq15g2\/55575NprrzWUhM1Aixk8eLDccssthj0oOBiXt5\/voAr33nuv1KhRI9R05s2bZzbZpgZ+DUtLS2XAgAEGrGg+DzzwgKFqqq2EGTDSVsKs0h91XHBALbZt22YAwAbcf\/\/9hhLAKgDI5ZdfbgDjB44bb7zRvP1hiw2OF1980VAev2KzkAgcYVc4iXouONggBEEoQf369Q1vRz196623zP+tW7eOUYlUUw6ETmwZfmX16tWxsSNwJLHpYZt6gQMVU2UItArebgxUGLwAjC1A8h1ggrR\/8803CcscyrLCyBx23UjmCLvDSdTzAkf79u1j2kft2rWN8IkKixxy9NFHx3i8Uo6mTZsatffll182wuRDDz0kF154YYVZrVixwlAhVNxu3boZoXXdunWmP+QJ14Zhd+DaQiJtJYlND9vUBodtt1C7hvajgiYqpQqA8ewcd911lwEBVAU7xyeffGLM72vWrDFdIrvwGbBNnDgxph5jBb311lvL2UjQgiZNmhSznkZ2jrC7m2Q9P3DYKiZDqOZiC4U2OBKxkGK8wp7BeQslspAmuYnpau4HDts4xdhoMZisf\/vtN2PZBDyVOVvxOh+hf8ZDy\/n44499HzU6W0kXCnz69QMH7OO2224zWkqLFi2MMFpYWFjOWup3Kvv1118bMGHu\/vXXX40cgtUTDQN2Eu9U9sMPPzSW1s8++8y01VPZXr16eZ7oRhbSDAMmGi71K5CTR\/apX4aoR68ViMAR4cJ3BSJwROCIwBFhIPEViChH4mt2xLQIBAdmaFQ8DrdwvEVvv\/TSS41vJoYhLTjpotJhIsanknqchOL3gDk7Krm3AnHBsXnzZmMW5mTxqquuEs4vMAvjpc3ZxPTp0+Wkk04yT42BaejQocYJ94orrpCVK1eaenhPjR07tpxXeO4t06EZL578s2zbsFfqF9aQ80YW5OpjhJ53XHDg+s95A76PXbp0iXXKYRRnFVdffbWxPG7dutUAQ62BGh6AFxROuDjjeB1qhZ5lFlT8cl6Z\/GvYhthMzhtRcNgDxBccv\/\/+u6EMsIonn3yynIUQVkPcyFFHHWU8m4i\/gH24fphqscS66Odp9f3338srr7wiPXv2lJYtW\/q62yWCD85MiB\/hxFSdjRNp79alv38OKZXVC3bFfvpLv3z52\/TCSnWbjvml8nn1oQJlDq+n17MKjq8Bx4IFCwx18QoWgmogg3gF7tA35xKwLAKSMDfXrFmzUgtuN9q9e7f89NNPcsIJJ6Ssv3+PWycrZu5NCTjSMb9UPm9S4IBdwE5w5YdicHxNaAC+E7aQyiD89vrrr8vs2bOF2A+35Ao4oByli\/ZF4Ij36hKPgcMuIIDtHH\/88QYAbLIXdcAjCw3GpSo6Rq6AY06ftbLpPwcicPiBgyNqTj05VXz88ceNbyYFvwkcY+64444KcaGccr766qtGRlHNxu4f9VgpUNeuXVPCBuDphEISfpAKNrV+6U5ZPOXncsuCzNGuX7hQBi8ZJlXza9y4sXnOjRs3SkFBgdSqVcvIWamQtULLHKtWrTJqLdFeyBenn366eWYESjYXcEQl8yuA6QDvebzToObHHHOMAUvYGJx4Mw4EB8YtbBhQjObNmxvBER8Jly1MnjzZkzJkfrmOnBF5IQnLJMqvSZMmAhWBUmaMctjGrfvuu6+CwKkyQ1D8xpGzZZl7Ul17fGcJ70yVdhZKW1Hhs02bNsZO4UWqInBkDgx+ml7GwUHMB4CYO3eucZWDpbgF9gIpQ52NKEfmQaIvJsoBlCNjAinnKhok5PfYgOayyy4zsaAROKoOHBga8Y3NuEAa9MiHO1vxSxSTjpfBdjxm3YMCqHXtq0wgPZLBoZtFoLQdaK2bkmgAdry19BoL4yKHn0EGxMHtp0mTeq2lcat60m30oVPyVJRAVTZokMOZcmDdfemll3wtv36\/Ba2Z1+9eYylgcIHo169fhWa69sPLBkqLfYfSUTXrWEfqF1aXGs3qS0Fxx8pMJdYmJ8Fhp0DgSew3WKPazz77bHMoiMVQWQNxJbyFFAKrcUd4+umnzamy2w+fg0z\/iay8PWfkA8YmVcSUKVN8UziEBcf4Tb2l7Z6KFKNgTMekAJJz4HBJrRc51jcKQY03TjdG5QQ7LYJ+59WPXU\/7SgQQWteds9f4br8K6PXr1weeaPuBI7\/\/n6Vw1l8rM2XTJuPgmPhGqazfsluaNqwpo3pWPKWN9yS6qCRpscmsF2vTw0BkhREjRhh1XOWGeP1Q1+bx2rc9ryBB0a7rN5YLWLsNc8eHhhJP8NW5+YGjcFYvye9\/KLa3MiWj4Ji7bKMMnbsqNs9RPZsnBBA\/Mu9FfuOFHOqGuSTdbyN1wvHYGRuKLciVDdhAF3D0F0ZW0\/nYwLY3WfuY9Kf+0vq\/\/0\/OW7tLodTpXJgUS8k45QAYAETLlUWNZcaVhw7wwhR3c9w2LunX+u73LDonyKRIsDMZB\/F4Hc9Wb\/XN9gOHH6D9AOo+Uzy5xzafN164U6p\/uSMloNA5VCnlABgAJGyJpz24fdh83Q2eDkM5yCmKn6yfpuBSGcCBl\/3bb79tBFwVkpOhHDxTWHBk\/GwlzKaFIY\/l+OkbpbJkzVbp1LJBQiyFPoI2VdmELczhXgCFoKj2ElbmiOfE5M7FFjpxtNYUlPyPZuTKSa7MQXsS2ukcdc3izaHKzOdhgBGWd4btK0w9V\/JXINib76ed6AbpxtJGhU8v\/q5shnq2l5ufwYp6KvQyB9JD8dlPW\/Ea3waRH1vUdbK1ssh8\/seqxBMM\/ayX9kITf8PbrJmNNX21n7pqaw+6Ma4W4cocNjiUPWjuUrVzYGOx5R7XfO6VS8R+gSLzeRhykmCdsMJgIt0CDptyuJ\/dvthYzSvmlxgmaPwq9ecImlxVsJUwcwqqky5wqJM148eTOfzObILm7QUwwjoIKCPVpvrM4hXGX7Ilo9pKspNNVft0gcPWVmwW5WVIS8bi6soc7rrg68tfsuWIBEeyi5Yt7RV0rv9uRDmyZYeqcB6JmhESnWpEORJdsSyqH4EjizYj26YSgSPbdiSL5hOBI4s2I9umEoEj23Yki+YTgSOLNiPbphKBo4p3JN2hCe65TSIhDxE4qhAc6Q5NcI\/jXd\/XoEf3Agc+p7gukGCHVFCnnXaayeqINxkR+Fpsfxd3HPVFiewccXYgnaEJfiZ8P78Or2m64OAzG8uNEbgrkliHrEoLFy40AOF6MwUIAfIk4eF+OjcGmvQapODKSXDkUmiCsiXcEdVLTC8CJGzCLZUFB3fpci0ZV4uRSYkYZgopNMgjCzUh9RYhGxRibkjTpdeOeIEv58CRa6EJXv6mfsQqnkd6PMrB5nM3DI7MZERwr0ElWwL+K6T9hKIAGA7+uPwQTzlu8s55cORiaIKXp5q7EfapbSIhltqOMAZ8Xv3K8uXLTSYE0nMBDk0VyilyvMuWM045yuaPlf2b10m1Rs0lv+8\/gmSucr\/nYmiCggOvrqCLjcMEMtkLouAg2xIpx70K+dH4nVyvsBXkCc2gcO6555prVcnvQfrQ7t27myyRUB9KRsGx\/b3ZsnnGwNgz5P+9JCGA5GJoQiLgYGGCYme8wMGRPakw3AL74GrUu+++O3ZlOymhlM0AkksuucRklyaPKTIIGhq3eXM3XkbBATAAiJa6XQdIo6GHYlfDlFwMTagqcACMN9980yTz69GjR7n887CZm266ycgo3F+Xl5dnlp+MhDg4c0cdAMkoOFzKATAASNiSi6EJhE6SZtNlK8nGs7BmfmyF1OSwEWQM7sktKSnxFTrdtYfFcMkyQm5GwcFEkDl2rXxParXpmhBL0YfItdAENAQvcHgJql4ZAuK9OF4CKem6SAP12GOPmYxL2Da8VGaEUr4nf71dbI0p4+AISyXi1cul0AQ\/cNhg16BpvkvEt1TBoffn6iXKgAN2gvrqlax26tSpJnEwNg5NNMzYsCJUW9aX33ISHMkCLF0OxsxLNZKg0IRkn8FmK5zHkKwW4bO4uFiIj8Er3S+LMab1QYMGlRNS6U+zR2Io41r2CBwdOqRin8rlf6dDDU0gZVS6im0+xyLK7VllZWVGA9E7b+yx0UCKiooE1oNcgXyh2gqW1Tlz5hhhFOoBRYnAkaLNcwOpE2EPlQWPDQ5YCoaueMWeEwAh3xgAwVLqdT3bEQmOym5GtrWLjuyzbUeyaD4ROLJoM7JtKhE4sm1Hsmg+6QbH\/wDgX3D9pOJo+wAAAABJRU5ErkJggg==","height":191,"width":318}}
%---
%[output:7a046379]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAIcAAABRCAYAAAD4gLU4AAAAAXNSR0IArs4c6QAAD4FJREFUeF7tnQuwVlMbx59DkhGOxCSku0K66KjpCx9q3KbbKCJEbkNConJPiKKIymWEMiMJM2pojNtEucWUESmlpty+JCEV4nzzW9+33llnnbX3Xvt993bet7PXTDOcd132etZ\/P7f1PM8uq6ysrJSsZRRwUKAsA0eGiyAKZODIsBFIgQwcGTh2PnBs27ZNRo8eLfPmzcttrkePHjJp0iSpX79+7m87duyQO+64Q55++unc3xo3bixPPvmktGrVKvc3VK+VK1fKjBkz5J133pFvvvlG6tSpI0ceeaT06tVL+vTpIw0aNHAS8vfff1djWGPp0qXy66+\/ym677SYdOnSQ008\/Xfr16yd77bVXlbGbNm2Siy66SD755BNp3769TJ8+vdr85h6D+qSJ7ZLlHC5wuA5948aNcumll6pD083ut2XLFrn\/\/vvV4QImV2PMxIkTpUuXLlV+\/uqrr+SGG26QxYsXB57TfvvtJ3fffbecdNJJUlZWpvpl4EgR1i5wsNxDDz2k3lbdlixZIueff7789ttvTnAABsbwL6o1a9ZMHnvsMWnRooXqumHDBrnqqqvkww8\/jBqqOMfUqVOle\/fuGTgiqVVgBxsc3bp1k3fffVfOO+88ueWWW5RIoD3xxBNy5513Srt27eTTTz9VfzM5x4IFC+SSSy5RHIMDvP766xW49t13X\/njjz\/k7bffljFjxsi3336rxg4fPlyuvPJK+euvv2Ty5MnqwPWcI0aMkJ49eyqxxnyIKbjNW2+9pfpUVFTIlClTZP\/99884R4HnHzrcBsfFF18sb775puy9997q7W7YsKGYfZDvyHUTHIceeqjceuut8txzzykwPfLII3LiiSdWWxd94vHHH5czzzxTjjvuOAWi9evXy4UXXiiIFcDG3Icddli1sYgsxM7LL7+sftOcLRMrKaLDBoeW+++9957MnDlTOnbsmDvAP\/\/8U6644grFQRAvmnOUl5cLoIKjdO3aVXEBOIZPW7hwoRJXNIA3atSoHLeyx5vcifV41p9++ilTSH0InU8fGxwofKtXr1Zv+M033yxDhgxRnITDwIq5\/PLLZdiwYUo8aHCwLm8\/f4MrjB07VurWrev1OLNnz1aHbHKDoIFr1qyRCy64QIEVy+eee+5RXE1bKz4LZtaKD5X+38cGB9zi559\/VgDgAMaNG6c4AaICgPTv318BJggcl112mXr7fZsJjmeeeUZxnqBmipAMHL4ULqCfDQ4OCEUQTrDPPvso2Y55+vrrr6v\/bt26dY5LJM05UDrxZQS1L7\/8Mrd2Bo4CDt13qAscmJhah8Cq4O3GQYXDC8CYCiR\/A0yw9i+++CK2zqFFlo\/OYfbNdA7fEy6gnwscnTp1ylkfe+65p1I+MWHRQ3bdddecjNeco0mTJsrsff7555Uy+cADD8hpp51W7amWLVumuBAm7gknnKCU1rVr16r50CdsH4Y5ge0LyayVAg7dd6gJDtNvof0aeh6taGJSagUwzM9x4403KhDAVfBzfPDBB8r9vmrVKjUlugv\/D9jGjx+fM4\/xgl599dVVfCRYQRMmTMh5TzM\/h+\/pFtgvCBymickS2nIxlUITHHE8pDiv8Gdw30LLPKQFHmJaw4PAYTqnWBsrBpf1L7\/8ojybgCefuxXX\/Qjzsx5Wzvvvvx+41exuJS0UBMwbBA7Ex7XXXquslObNmytl9JBDDqniLQ26lf3ss88UmHB3\/\/jjj0oPweuJhYE4CbuVXbRokfK0fvzxx2qsvpU95ZRTnDe6mYf0HwZMtlyyFCjZK\/tkyZDN5qJABo4MF4EUyMCRgSMDR4aB+BTIOEd8mtWaERk4as1Rx99oBo74NIscMf7VNbJu03Zp0qCejDq5WWT\/Yu2QgSPhk5m1+HsZOmt5btZRJzctWYDUKDi+\/vpreeGFF+SMM86Qgw8+OEdQ7jvI\/eC2UwcKu87Qtx9jffv69guaE2AAEN3Ormgkkwe08tpPWs+ZL\/5z4CCmkRtN3Me4polsIpr6qKOOyuVaBC1CzATJRYwnNoJDPf744+W6664TrsWDGvcR55xzjtiRVNu3b5fvvvtODjzwQKlXr17geN9+TODb17df0JyIlPGvrq0Cjon9mnntJ63nLAgc3EfcdtttKgyfWEvuHubMmSMfffRRlVwL1yIAiZgIwMFBH3vssYoQxFAQ2EuYnr7FtMfvjODoNXWJLFq9eecBBwdLHAM5FbzxNB1ST1wm+RlBUdlcNBE8y9U4sRM6o4t0QjLNjjjiCBX\/sPvuu1fD1s4GDlvfYMNTz24r\/dqVlybn2Lp1ayU5p3AAO8+U8DYCWHSov4tzEEVFngjBvHbuKRHhpAm68kCZKwwcZKohrsLECvoB6Y6kGIT107Lcp28hcyJOTK6BtTJlYBul7\/isncRzoruZ+lu+IoVxZRs2bKgkrpGMMTv6mowt4i5J+TvrrLNirUMUFQlDK1asiA0OQvAI+ycKK2vxKEAu77333psIQMpWrlxZCQAIhBk4cGCVJ9FR00RWo5zGaYCCsLxTTz01MOEniHPov7PJgw46KM6ytbovLxMqQFSqhC+RcuAg38PmDjogBcslTk6HDp\/74Ycf5OGHH1ZpAa6mQYCVY+Z9oAiTTZbUJn2JUer9guiJOyDMJRC078TBgUmMOMHyMbPKw8CBrkNGl27UrCCaKwNHPLhqcNj0RCfzTfM0VwwFR1yxgoWCrsDhPvjggyrpWFsvYeBA4e3cuXOuC5skOy0DR37gsOmZN+dISiFdvny5smwoTXDfffepROaoppHefvA4ady6fe4uIkgXiZrP53fSCR599FFnV6wrU7TSl\/RJclv32GMPn+mdfcLWZIDrJdA0MCeMypdNmm5lhZqyPDxWDc4zyh+QGkhgr09jM\/1HT5WtnYbkuuNu\/nf5f2TU0HQ4BwfFurZ5rbmkqXslCQ7Xmi4aaT0Pfc0uTaVBZoNYz5M4OCg1iRMMJ5apI\/g6wbTyCcdA1hHp7dvYTO9pS2VHw6p1LepsXCH1F05IRawEgYNntn+rCXCEPZ\/rGU1aa3D0Gj5JJg\/r63sMgf3U3QpKJBVt0BUGDx4sLVu2dLrPdToAXlFQTT\/c45icuM0pkGY3clTxnJpF3HQf+x5C\/70mwWGKEQ0O7opg\/TRXWoMtNnSytBZFUQeu9+3iXjY96fPiiy+qakR2qoQGx5buI1XtMpxwtHxDBwIv3o4++mhlMaAoaqXSBgeEwruqq9a4IBgkJ12u5nzAETd2IkysjBw5UqUvak+vPnTzsEnOJt9Vs3x7Pk0jaKPNf19w2HPHffVNcNjcGHGNKz9Oq7Ere\/tq23xoX86RT+xEmHJocwXXoZqHj+eYl8MEgov4YWuayqgviIIOOAwcjIkbW1LS4HDFTkS9HUEHoBVBiKiV1SCdgzf82WefVf1ee+01Zb67xE1OfAYowfYhBz2bC1wupTQKHHG5R42BIwmxYs8BMCBAWAt7O3W1Hv0201e9cVbFHxMcyH2zOItLL\/HlCD5iRYOYqw7bo12U4KD6L7miyGJuYUkcxomFC5w6Ga7mAgeKE7GXvmKFedE5Fq7aLN1blnuF44UdlG0K+nAOWyl0maK+4PBRSEsOHGSuDx06VGnIIPrzzz9XZaLxkN5+++1OS8UWCf9qUS7zhnZUhz1n\/gLZNHfsP27KRimbgNGlcNrgd4HM188RBaSSAgdmMMDQZQa0yUp5RRQ2ZKOrWk6YSEjamWMeXpQTDMvEtDLwpprF5BiPXwhrhbgJFFKa6UV1+Ut8wWFWD7A9p1rsBek3tljhhTPjS3zErkmrgnUOHgjxYRd41aUQYLtBJRyHTV8osxcskxGDTq4iEtIGR1z3uennsE1zV5ltu08UN3CJXVuPCfKxmGPN64iuXbsomsYVu4mCA62d8ELXVwjyjQRLExxx7PxS65s03QrmHLwV1BwnoJiySDYLf+WVV+Spp54SispHyWb9e9KbLLVDzvd5k6ZbIuAIkqdRppneDLe55qcquPrHnZ9d2ceDSRA9840rrVFwkNQECIJiRTNw5AcOexQvH\/\/itlTBgcgJEys8LADhn9l0LKTNUeJurrb11xzXjr2tMc5RiEIadHhRHKW2HXqc\/SYafV7op0MLMWXDNu3iKHGIVFv75sslXPQqWKxoJxh5raY3NMoJVlsPr5T2XTA42Oz8+fPlmmuuUcVgBwwYoEpBR7nPS4lItfVZEwGH6+Jt0KBB6vKNuNKslSYFEgFHaW49e+ooCmTgiKJQLf49A0ctPvyorRc9OPigDhlceEtx8hAa0LdvX1X7w7zLySfgCOIQi8LF4bRp06pFc+c7ZxTR09jTunXrVDIZViIls9q0aaN0PkIQzNoocfZU1OAgsQcvKdfXVA3iS0x84AZLiHJSpFzqLPy4AUcQ6Y033lDxn\/gGXDVE4s4ZBQp+T2NP+JqIOSGZjNDBAw44QHmm586dqwDCJ800QOLsqajB4ao4BIH5rBblHc4991xVOmLz5s2xAo5wsJH9T2krCqu40ifyDWKKAkjSe6LADWW3+JwYYZqNGv0vhhbww3HhJtyKk2oSd09FC46\/\/\/5bcQZiUwnOMWM1YctEYO2yyy4qAotkLN+AI\/M7J9QdYR2KxdicIw3Pbxp7oiwXhfkOP\/zwyOI7cfdUtOAIewP1V5coKwA4XnrpJe+AI8DBxVSfPn3kmGOOUf\/tCjlI484orT0FzUvpLF4aUl0RN3H3VJLgQOlCnHDdz+YLCTgKyxXJN4gpSrS4fk9yT8yPuCR3mTqviJW2bdvGplPJgUNn9GOpIHZQvqJyUcz0RftgwsCRbxBTXHAkvSf0DQr53XTTTbnPtFOjIy6dSgocfGyP\/F1Eg1lOKu6mzcOraXAkvSeAQRYeXLVHjx5VLkPj0qlkwGEWh8EvAZvULWzTUQFH+YAjak5fzpH0nlB4ESPoGHwbd8yYMVXutuLSqejBwZuAbQ7HaNq0qbMGSFxFy4dzFDJnFDjS2BMmLdwUZx5Fg\/Ft2NWI4u6p6MFhOm3uuuuuahHuHERcE80HHIXMGQWOpPekP5wMOBAnlA51VQ+Mu6eiBodW1LDhSYwKqohXSMBREKstZM4wcCS9J618UsERby+e5KCyknH3VLTggE0CiFmzZqn7AUSK3Sgx1bt3b6lbt27eAUdhcjjpIKY09oQDj3smDh5auCooEYRVUVGhyBdnT0ULDu4gKLvNXUpQMyvu5BtwFAaOfOcMet409oSIwtcT1sxaHnH2VLTgiJLb2e\/pUyADR\/o0LtkVMnCU7NGl\/+AZONKnccmu8F9ksKf+m7RttQAAAABJRU5ErkJggg==","height":191,"width":318}}
%---
%[output:7b704838]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAIcAAABRCAYAAAD4gLU4AAAAAXNSR0IArs4c6QAAEPFJREFUeF7tnQm0V9Mex383c6JUvNCoDFmmEJYXosyvUkvmeSZJysuQ65oLRQhFpqySDI8MyypDFCpWhshqlqEkKZKQ7nuf\/d7+v3333Wf6n3Nu92rvtVqr+z\/n7PG7f\/Nv75Ly8vJy8cXPgGMGSjw4PC6CZsCDw2MjcAY8ODw4\/nrg+PXXX+Wqq66ScePGFQbXsWNHGTx4sNSpU6fw25o1a+Smm26SkSNHFn7bbrvt5NFHH5Udd9yx8Bui16xZs+Txxx+Xd955R7755hvZcMMNZbfddpNOnTpJly5dpH79+s6J\/O2339Q3tPHRRx\/Jzz\/\/LBtttJHstddecuyxx0rXrl1liy22qPDtsmXL5Nxzz5WPP\/5Y9txzTxkxYkSl+s0xBr2TJ7ZrLOVwgcO16EuXLpULLrhALZou9nsrV66Uu+66Sy0uYHIVvhk0aJDsv\/\/+FR7PmzdPrr76apk2bVrgOjVo0EBuu+026dChg5SUlKj3PDhyhLULHDR37733qt2qy\/Tp0+WMM86QX375xQkOwMA3\/IsqLVq0kOHDh0vLli3Vq0uWLJHLLrtMpk6dGvWpohxDhw6Vdu3aeXBEzlbKF2xwHHjggfLuu+\/K6aefLtddd51iCZRHHnlEbr75Ztl9993l008\/Vb+ZlGPixIly\/vnnK4rBAl555ZUKXFtttZX8\/vvv8vbbb0tZWZl8++236tvevXvLpZdeKn\/++acMGTJELbius0+fPnL44YcrtkZ9sCmozZtvvqneadu2rdx3332y9dZbe8qRcv1DP7fBcd5558kbb7whW265pdrdDRs2FPMd+Dt83QRHs2bNpLS0VJ5++mkFpgcffFAOO+ywSu0iTzz88MNywgknyMEHH6xA9NVXX8nZZ58tsBXARt0777xzpW9hWbCdl19+WT3TlM2zlRzRYYND8\/333ntPnnjiCWnTpk1hAf\/44w+55JJLFAWBvWjKUa9ePQFUUJQDDjhAUQEoRpwyadIkxa4oAK9fv34FamV\/b1In2qOvP\/74oxdI40x0Me\/Y4EDgmzt3rtrh\/fv3l3POOUdREhYDLebiiy+Wnj17KvagwUG77H5+gyrceOONsvHGG8fqzpgxY9Qim9Qg6MP58+fLWWedpcCK5jNgwABF1bS2EqdBr63EmaX\/vWODA2qxYsUKBQAW4NZbb1WUAFYBQI4\/\/ngFmCBwXHjhhWr3xy0mOEaNGqUoT1AxWYgHR9wZTvGeDQ4WCEEQSlC3bl3F21FPJ0yYoP6\/0047FahE1pQDoRNbRlCZPXt2oW0PjhSLHvdTFzhQMbUMgVbB7sZAhcELwJgCJL8BJkj7F198kVjm0CwrjsxhvutljrgrnOI9Fzj23nvvgvax+eabK+ETFRY5ZIMNNijweE05mjZtqtTeZ555RgmTd999txxzzDGVejVjxgxFhVBxDz30UCW0LliwQNWHPGHbMMwKbFuI11ZSLHrcT01wmHYLbdfQ9WhBE5VSC4Bhdo5rrrlGgQCqgp1jypQpyvw+Z84cVSWyC38DtoEDBxbUY6ygvXr1qmAjQQu6\/fbbC9ZTb+eIu7op3wsCh6li0oTWXEyh0ARHEgspxivsGfhbKN5CmnIR8\/o8CBymcYq20WIwWf\/000\/Ksgl4ivGtuPwj1E97aDnvv\/9+4FC9byUvFATUGwQO2McVV1yhtJQddthBCaNNmjSpYC0N8sp+9tlnCkyYu3\/44Qclh2D1RMOAnYR5ZSdPnqwsrR9++KH6VntljzrqKKdH11tIqxgwvrlsZyBzlz079\/rrr5f99ttPTjzxxGx762ur0hnIFBz4C\/BfPP\/88yp+wYOjStcy88YyAQeGJpxLAOLLL79UnfTgyHytqrzCTMCh\/QxI5aeddpoSzAiC8ZSjytcz0wYzAcdzzz2nHFqnnHKKktQxU+MA8+DIdK2qvLJMwGH2WjuZPDiqfC0zb3CdgmPga\/Nl4bLV0rT+pmpg\/D9psb\/l735HtkhaTYX3db\/iVBLW9yz6EqcPvBPV52L6ss7AMXraYukxembcsSd6r0+HxtKnQ5NE3\/DykImLFEDpW1al2L7EaV\/3l3fj9Lnfkc0TbZx1Bg6AEWdAcSbJfqdT6zpS1rFhhZ\/xoaxatUpq167tDOcbN3OllE1YWkxzod\/Ql\/7t64W2bVYQ1U\/97rApy2X41OWJ+nty20Yy9OTWsb\/5S4Kje5uGMqT7\/xOWmA3M7d999500atRINt30v2zMLL3Gzpax07MHB30Z8I\/GoW2b\/YjqJ++OmrpI+v5rQexF1i96cIio3cFEmGX16tWyaNEi2XbbbZ3gyIuS0Y9BXVuEtp2kn7xbbF89OP5DJVy81YND1IapEWyl54hJMmbijMSkMc4Hf29ZTwHE5uWkRpKO4GIrA19bIJPnJuPhcfuCUBrWdpJ+8i79pL+usrZ2A1lbu6K8lZqt4BchigqXNXyPaGoyuPbYY49CfmfQZKxdu1aZz0liJqSuVq1acsghh6hsL1ceyNdff60yy4iy8iXbGSCXd3zdbk6AFEU5tCeV1D\/yO4h3GDt2rHzwwQcV8juDhoGFlIz39u3by0knnSSff\/65ylZv3bq1Shm0AUJgDNbUO+64Q7bffvtsZ2c9ro3NxnyvbPdPWdOwcvadSxYLmy6lrXCMAbGT5HGy4yk6jY9cENcC60p1FjuLjLNNH39AxPVFF11UCMw1O6HBEZXvsR6vc1FD1\/Nav3OpzKvVrFAHBjCoRlLjYMmqVavK2fWwEvtsCxaYoFmdXujqsTaXE4IH1dDl+++\/V2kCJDjbyUIeHEWtfeRHQeBIyk50QyVLliwpD1pEssRxooV5WHXM5hFHHKHkCH3+hAYNgAE4nnJErm3qFzQ4OvUeLCPn1y3Ul5SdFMAxa9ascgBg73xe0AtMNhfCqatg0SPtEGGWUH1AgmeWvFNkD\/M8C\/193pQj6OyOIDZmpjbqPoalR5rxn0F1mu9QZ1Suq34fKutKrYx6ThvmvE5c8TeZNGe5tGuF5lacr6lEg8PlRdUdorNheaQE+zz00EOKLelCYC5JQq5jCfIER1CfdZv2ogNonpnHLuk6GIvrOCbARHISqQokRpEYvdlmmxXG7uoD7SDb2cdN8ZEJZhfYop7ntelSg0MfmUSADxHaBx10kLIGPvnkk+qAE4Rc8lSriq2wcE899VTgoprPNGW88847K+3WoGd6odDomjdvrkBiL7irDxowsFkzzsXMo2WObHBEPXfNK1TcpD5E0evDbJLwrlBwxGErr7zyimI5sBHAoWUOnfCzySabKICYB6blSTn0rnbtUHtiiok9MfsO6MmisxfctQAucOj2SX3o1q2bku9MoEY9t9vRfYOCw8Z0wfAX99wRs87UAinkkuOWyEeFzJqFXQow7IVKAg4z5iMO7zR3WlQcq0mudfa7yR5ci2yzIRdbsr\/T7SxcuNBJ0Uz5zkXF4jw3ZQ60y3333bfQjaIpR1pVlsnBOgpvxqmVJTjsmI+48QgafGZfggRCW3DU37h4v2v3RwGd+Rk2bJiqNsyuE8bikoIjK\/tRwQhGTql52l1cIxhCFmzllltucbIVTOn2cUpRE6oXyPY+FqOv25pImBZi83cbUC6WpQHjEkxNcJoswiXcV1tw4FfBRsGBqWeeeaa0atXKaT7X5JGUP80qtOn9pZdeUgeYoMoibzCRHKZmAi6pVG1TjmL19Tgagc0KNKg0mILUY\/2dK8XSVadLgI1DGaLAY7KVTCkHFXPOBIjmsFWESk7l428ka3a\/OcEmOPTvsBXUWU7v5TRgJoujl7Q53pyouJSDb5A54urrQRqBbtsWQCH5RM3bqqg5Vv7PczYO\/iDXxMcVbMOE5ajFj3qeGzjSOt5YbOwkjRs3VpQHcNxzzz0KVFVtBAsTEO0JDgOpqbKyScLqtQVONokLdEn6ZlOddQaONI43KAUyB0ccAIhtttlGjQvXPWoeSU5YX7WKmwfCzYkMMmAFGce0wGhTBHMhqZ+xhBkDTaCR3GXn7mg2FaRBRS1+1PM85jW1443jozlGEWCYLETvJqgH+bMcw5RU5rB3T5K\/TS1BfxfEi13ajSm4xrGd2IKpfZRklEwStfhRz3MBR1rHG7YM9GqXKhu0mElkjiSAWN\/fzXpeCxbSYh1v7FAOPbnhhhuEoB\/AQgwIgT+XX3657LrrrpXWLOtBrO+gyIsip\/Kt4HAjSIijGhFCMYJ1795deWWxmCLoIpBy1HSx2opf+PgzoDcdMTjm1R8oCvxLWlKBw9T9ifriuCXt4Fm8eLHSYBDOMAmbF+R4ypF0meK975Kd+BKw8C9pSeV407cg4XzDKKZP2dOdwK9CLOpjjz0m3FVikz8b4Uk779+vOAPcLoUx047NLZpypBVIiTDnJB+XFzRIyvfR5\/nBGnYCOIphI3avUquy3CPCUdJoLGYMATIInXzrrbecmgwA4Z8v2c5AsVTC1YvUjjdkC2wC+GPQWLRswcLz+z777KMOio17VUW2U+VrSzMDChxpHG80DvVAayGGwNRWMIA98MADlSLB0nTYf1t1M1DIsrcz3tjxaB8suDZ9u7yydBUW8sknnygrKScEk2549NFHq9uRcGX7UjNnIPMjGGrmNPheB8ocfmr8DHhweAwkmoFqz1a4UAc1GY8qRh4srscdd5y6ZdoMaEbu4XB6Iq0IWOI97nTjBkfTI2zPDjISxrr777+\/0sH3xdYZtQJ5jIngZSzRxPMSRrHLLruo8RM4TQaALknGVK3BQb4tVlTc1URhcRMTF9yQwY+giwCss\/RZ5B49eiifgpnpzz2wpoptTtLrr7+ubnjENuBKXkpaZxQoeJ7HmDCbYzbglggi94ipwWr94osvKoBgh9IASTKmag0OVxASE2wHEi1fvlwBQ99rom0t7CLyfAmwMa\/nwgaDio1pHxeAKzId7S1JnXGAwTtZjwnnJ1eREeYJ1eTMM61BQnGhJrgv0D6TjqnagoMDYaAMsApC+827TiDLnAyAHUXHd8I+7Bul9d0rfKvvjDVTEQiIph3ua7MpB7sxbp1xgZHHmAiP6Nu3rwqNsKPa7UT4pGOqtuAIm3B96xJZXIDjhRdecCZPwV+hGsggevEBB2b9Ll26qGs\/+L+dK0vbQQlZrjrjgiOvMQXVS5QeAMdCDbtJOqYaCQ7YBcFJeCAZfFjWHc\/gv7ZnWE9oUNBvmjqLAUuWY6J92CVpkc8++6waO6csJR1TjQMHpJKjqdBUdEBzWFR3VPxnGDhcFIWJj6ozKTiyHhPUjetQr7322sI17cTZJJ2nGgUODorBpA9rMH02SQdtLt66BkfWYwIY48ePV1S1Y8eOFTS1pPNUY8Axc+ZMpdbqYx0gk1Gsgedp2EoQ5YiqMy7lyHpMCLywEWQM7sYtKytTyWnFzlO1Bwc7Ad0cisF5GPBRbns0S1JBKw7lSFNnFDjyGBMqLdQUYx6pItg27BMDko6p2oPDNNqQrG0f88BCJFXR4oAjTZ1R4Mh6TPriZMABOyGhynVYS9IxVWtwaEENHR47RdABJNq4Q\/S7aQ0NMoLFAUeaOsPAkfWYtPBZWlqqrL1YkoNO8Uk6pmoLDsgkgBg9erTyD8BS7AJ76dy5s4oye\/XVV1WeDLdPE3DE3fOY2YPM53H4cLF1BoEjjzFhwMPPxMIzF2aUv+4Hc9K2bVv1Z5IxVVtw6HNM8aUEFfM0HpdD6dRTT1W+BVMos+sKk+CLrTOov3mMCRaFrSesmPm5ScZUbcERxbf98\/xnwIMj\/zmusS14cNTYpcu\/4x4c+c9xjW3h3+gRag3GQa4IAAAAAElFTkSuQmCC","height":191,"width":318}}
%---
%[output:01ec9e76]
%   data: {"dataType":"text","outputData":{"text":"\n=== Resultats du SNHT ===\nVariable analysee     : expA_fit\nTaille de la serie    : 240\nStatistique T_max     : 8.6009\np-valeur (bootstrap)  : 0.1780\nDecision (alpha=0.05)  : Pas de rupture significative\n=========================\n\n","truncated":false}}
%---
%[output:957844cb]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAALIAAABrCAYAAADTu9iPAAAAAXNSR0IArs4c6QAAHihJREFUeF7tXQl0VMW2PQGEGJn5AhkEwhLFCUFRAoqIqDx9DkyfCA8fIiCDDA+RGYQICIJ5DGE9JgFRlqCIKPAcPiBEUMEBRf08FWTwB4LIIIMaBMxnVzzN6cq9fe\/tezvdnXStlZXuvjXXrnNPnTq1Ky4\/Pz+fYiHWA1HeA3ExIEf5CMaqr3ogBuQYEIpFD8SAXCyGMdaIGJBjGCgWPRADcrEYxlgjYkCOYaBY9EBUAnnLli3UuXNnNQCTJk2i9PR09fmVV16hESNGqM8vv\/wypaWlqc\/PPvsszZ071zdgMg1+3LlzJ3Xr1o0OHDig4iQlJdGiRYuoXr16vnyzsrIsfysWiIjSRkQ9kO+\/\/36aPHkyXXzxxX6ABZCvv\/56Gj58OP3www+0YMECqlq1qg+0SDds2DDiSWE0IXgyYILEgBzZCI9qIDdv3pxOnDihQIrQvXt3qlixIm3atElJZARIbimd5XD89ttvCugIPBnwWf991apVMSBHNo6j047MUhQg3bhxIz333HOqm5988kkF5vHjxyvwZmdnK4nL0lgfC1YpWDrL51BHOO3atWtjQI4B2fseYCCPGTOGvvjiC6XT1qlTh5YtW0Y9e\/akfv36OQJy\/\/79fXo211aqE9u2bfPp3rI1ui7tfUtjOdrtgYhULfbuJeK\/ffuIatcmqlOH6PbbC5ol9dq9e\/fSN998o36vX78+tWjRwqdOxCSyXRhEf7yIAjLAu3gx0bhxxh0LMD\/yCFHr1gVWCyzQateu7WfB4O9OdeQ+fZ6lgwfj1QTatesMffzxq1ShwlFasqQHGenIs2atoblz36UuXUZTXl6NQpNNb4HV5Ix+KIW3BREDZAx0y5YFkhiBJTD+I7zwwoVniYmniegRmjatDd1www3KdIYAk9mRI0d8EtmO1WLlys\/p0Uffp59\/Hmg4Eii\/ceOvKTe3t8r\/oovq2ZpsY8cWZGd3cnL88MIhekuPCCBjsFNTLwB40aILaoTsWsQDZjduLPj1qac20PDhaX4mtu+++66QpUK3I\/fq1UuZ3vTJU6ZMDsXHb6UqVY5TmzYP0urV1XwTq1y5XMrIKEVz5tTw++2vf72EGjSoWGiyYQIAnBkZ9iYn4pu1O3rhVXQ1DzuQJZigNmAwrQIkN8CMwd+zxyq28fNgJw+\/LexONqfxN2wwnsTBtbLkpAo7kCGxoBM7BSWDGWmdvpadTh4JekDD6iiC0\/jI04vJWXJgW7ilYQWyHHCnkggSGYNvB1h6s51OHo7P+VhNHqfxOV83k7Mkgxhtdw3kdf85Qh3nf6n68ZkHL6feLS6z3acMRqfSmAuAXo3J4GQSOJ08Mj4AzBYVM6nsNL7sLLPJeWrLCjqU2UFFrdIxg6r891O2+7ikRHQF5B9Pnqaui76mqR2uVP015LVvaXG3a6lGhXJ+FgizzpQWChmnZvIfFH9RKd+iyml6q8EzK9dOOXbSyjh24styOT5s5pig547upwPj\/0LV+zyvoh2a3YOSxrxDpasmmzYzJyeH8BdtISUlhfAXTHAFZEjjUW\/uolWPN6Ty5UrTo4v\/l3remkJ3XlVNWSF4UJxWDNaD\/6qRRAf3l3KatNjE57cUpPHRxU9Q8sQPKS6+AuVOvJcq3T+Iyqe1N2wrADxkyBDaunVr1PVFkyZNaOrUqUGB2TWQ52\/OoYVdr1GdBiC3vKKqUi8YyAClWTh7tmD2GcXBb\/w8mPRGaYzy08vW62T23ajeev6B2mDWL5yHBPLx1dMocdRbqkkAckLDu03VC971BCCSk82ldqShHBNvxowZpg5eVvUNGZDl5sbSpVsM67F+fQqNHp2iLBZmcQI14OGHU2jXrhR66aUcuvxy61fpiBFpPhs0b7Qgf1k24vCbBL9zHe+8M4fGj8+hTp0KnnN6PS3buPk56rVuXUEb9WDU5jFjUlR8mPdgjoREDgbIRh5\/qDfql51dUP+uXS+0wwoooX7OE9DMU9GqfNdANlMtrArm58Es2JCWF1VOF4rYIcSmit10eny2LJgtMJ3G1\/tJz9+pamEEiEC7i7rtXj9kgPrxBpLRmB49elQ5bt1xxx20fft2SkhIULur+\/btU45Y7733HjVs2FAl5XhG+YQVyIEWe3aBDFBh8HlxYyed3OGzu4ki83U6eTg+pKSdSeA0PtfNaHI6XewZAUKaG9HPLVoQwRmLLTCyDwHkd999V3kQIsA3e+LEidS1a1ffiRnZlyiPQTtr1ixq3bq1Lx7Szpw5U3kkYseV40UckFEhaX57tWcDtdBzGnjg7YBZShe7UlWvj9PJw\/E5H6vJ4zQ+v2F4+13PX5rfqg9+zXShh3x0IMvJgbeIVHGkqRA7pHgWCMjsx8JSGgCF\/zckMVxq4QcO\/5Y+ffrQsWPHVHfh6Nndd99Nubm59NVXX5nqwGGVyE4BaxQffr+DBr1BublL1eNatf6g7OxSdPBg4XN5UhJXrHiUJk06SH37Xq3SmZ3jM6tjoMkjz\/5hYMaOXURt21b26c6pqXV9ZwX5NMnq1atVUXgNp6cP8zM\/yvhG9fFicnK+OiBYGptNPp50vMljpFrAy\/Cuu+6ijIwMGjt2rDoyhj5CgLehLpFZtXjggQeiRyK7ATM67ZlnnqHMzEz68suqvp065Fm+\/GHq1i1OZf\/22\/+hGjWa0AcfXKS+JyQcomrVOtCrr05RB0yhpw0ePJhGjhypnnOe6HCzIDcfIIlYWsk6If3w4XNoyZI7af\/+y31ZJSefoaZNR9LcuSMIp0fgEw0nJNSjU6cRlJMznnbsqG4YX6+TnJxeOA7pQGag8uJR7w\/W6RnoukTm+Ga6M\/y\/SzyQ9U7FIPzjH8fO20DvMcRf7dr5VK\/eBzRnTqI61gTwAMhIBw83HGnCIVScw+vUqZPvFLUdMCMO+zvjMwDGq3x8r1kzj5KS5tO2bf192cGdtFevcuq7Hh\/PEhPn2Y4vJ5Mb4eCFRJY6MtcFk1RKZPkGiAFZGzGAEUeW8LqaOvUTuuqq3ioGnNybNk2ijIwC5wp0KnQzCeSlS5eqA6QIAHKzZs0KHV8yA4juGyHjwaWzVav\/o8GDTxGXMXFimfMLoIK3gx74TdK5805b8XnyOHV8MmuLDmT55tEtLUbb9WYSGeVJ9Q3f5cEGfMabCSrWoEGD6OzZswTVgk+xQ2\/u3bu3H32DbEPU68jcGOhc\/IpGowIB02sgcx3YxorvOF61Z88GSkj42EcboNepcuU2VKpUweRKSvqd3nlnLg0ZcpPvLREovn58y40UtgKEXHxCF0bZZlYLr+rhNJ9iAWSWxEy0YqUqGAE5GNUiUGc7qZPU0ZnUxaoNTgfabnwzQOiWFM7PypPPbrlu40U9kAEYLBiYFYhVh0CLNx3IwSz2rEBst07Ix0h39LpOdoESCBBsHcF\/2JLlgV67+YcqXlQDWde50EnMHATbJNNi6duWOpB1\/S3YbU4jPdCqTvoxKtYd8XaR7XNTJyfgcQsIJ2V5GddtvV1tUXvZkFhe3vSA2RY1+4AYlYKNKCNfEG9qZC+XGJDt9VOJiWUECLYVm3WCbmOWbxlJQoPfEWAtYhUQatXQoUNpyZIlaiua7eSwfixevJiuvPJKtYkigxGxTQzIJQai9hoaCMhGp1ri4goO\/GJDhFUrENswWBmQo0aNUnb69evX+8yebFsOBGSkg20fcefPn08DBgxQ3\/UQA7K98S0xsdwCWQJXBxycgho3bqxMoyxlrSRyVAJZ+h2EixfN7cyOdsTL9jNBze7dRP\/+92T644\/CkpAlcnr6BU81+EowX7QcR\/Zuw3NIbagSDOQpU6aozRAZJOUvJPLs2bOV81Dp0qX92E\/5TRCIOdVqXDxd7AVit7SqiFfPY0AucLaS\/NB2gPzjjxfYR6U\/iNyaRp7spgl9+brrrlMun3ZVi5ACuU2bNvmw4UIn0j25zNjgMdN1qlY2icFshoDZ2K5dOzWz4REGxnh20JaLCZ61SMNcxRUqVFADAWkARyA4AYFN3szBW3qrcVls7tLrZcZkr0sPdk9EvczSyH7Qy5F9Z\/aMJ\/5DDz2kdE\/0XSAndh4fSVyuCwA4MT322GMUFxen7PMIhw9X9Elk2Vdoc1bWZOrXbxWtWVPA9F+pUiXq0qWLIkcH9\/S5c+eUQ\/zChQuVjov8EQBoZvivXr262pWFA\/2KFSuUKydUijfeeIP69u1Ljz\/+uGobDpbCBaFKlSreS2QJZMkJLKmnqlWrpgAJ+lW482GQGfyyI3WJzNJRDo4sA2k5LywCAOTPPvtMcawhoMxLL71UTRp0IHxbdXusXiZPEimRkBf8MCQZoeSIk+1j\/wCzNLztjEOeqB9zK5vxKSNv2V\/yjcHl3njjjYXqxzuEsn\/tAHngwIGqnTizB0ChjlWq3KiAXLbsdkpK6kwHDrxMZ85cQTVqdKe8vDT65z+HkZTIkJzPP19wahsBbUaApIcKgdMgaFOtWrUUgAFknLkD2DERUGbHjh3p0KFDamH39ttv0\/Tp02n58uW0efNmpWdLYnVPVAsGMjtJM0BZikBaMFEgZqB+\/4YdIOvSUU4C5iGWnYdGYtDMACB3AfVrEYyAwvzHDHomBjfacOE4RmmM7iTR3zL620pXdSQbfo8ePZQTun4NhNkGEL\/tuM\/1dQjGDBsxu3btUuBjHfmXX4iaNJlM2dkzae9eqB4LKCGhKq1b9ywdP36Apk6dTG++OdNHbM4ehDyZ+VoLtIXfxJi4wMNTTz1FTz\/9tO96CymcdEyZ3RDgKZBZDeDXBXeWkTqAZ0aANpPIPDA6SJBPKIHMUkR\/\/aLuLHnlAkVS0RqlwW94K6BP9EHS\/XUZ0DrbvRGQedJY6fdWEtkMyKg3A05fkHE9oTYwUHUgIz3elkZpsQCElGbQSwEUFiDrheoDyd+5s42sElZAZikfSCKfOFH1vAvgaDp37iytXbtT3Qny97\/XogYNPqfBg9sXUi1Y5+PJEkgim7WJ88Cgjh49mvB61lnsdWmivzE4b7nO6NhxKH39dXvKzn6frrjiLtWWW245Q\/v3P0M44cISORRAhl568803+26zwuSD\/rpjxw61OMM6RAZIcH5Wrlw5v3Rly5ZVEt4o7e+\/\/+4X9\/Tp0wrYV199Nd13333qM1QPlI\/At2vhO\/LlsH\/\/fsXHEexWflxqamq+lLo8K6WOzJKNC5H6oFzhWgEZ+ZjpyA8+OPD84vBLqljxCL3wQkuqXv1Xn2qBo0MTJuQoXuIZM074jjchv2B1ZFyn8Nxzy8\/reQspIeES2rChgNbzb39LokOHpiqg6Xr1ggXraenSK+iWW26h7duP0+7du+nWW89SVlZjmj27YLcLabZu\/ZHS03fQbbfVolGjUigjo5tvTcF8zBdf3JumTPmJ5sx5yDdprCSy2USUv2MivvbaawRAyYDFX35+vloE4r\/RM\/5NxsFnBDtp9byNypJl6O1xRdASjNUCFTCaOXaAzGCWMxNAHTPmHJUpM8EHICMdGRLu0kvXUGZmRd81DMhPOufYsVrUrJmmeItzcjbTt9+OVAQx\/IYBkffs2b\/RrFknqVKlgYovOTNzBa1a1YgqVDhCH33Ui44c2abiYyFauXJDSk4eSS+9VJpq1RpL+\/Zlnz94OZA6dPia\/vWvob5dLd0KgjrAdvvrr0Np9OhblW7rBZDRb1g4wwqCNwwCjuizBNywYYOyPiA0bdrU9zuk7YQJE9TvsGC9\/\/776rOUnEZpdSkrJTJO6bDE\/uijjwzrI8EcNsosOxLCKg5cCgEqO7zInBe4H7B9z3eKWJUhn3N5SB\/IUYbP0oHEBIQmduPzuTu7dXLTFrtllIR4rjdE3LBxooPh8G0FEqOBAABwdMdJcDpp+CiQXbZPtAVpnE6yYNripN0lIa4rILslaOGjRcGcV4NHl1PHcKeThuPjv9WkkW1xCky9LTEaWedTzxWQA7Fx2qkKVApsPgWjIrA0t6uSOJ00Mr6dSSPbYie+3j+YLGiLU2YhO\/1cEuK4BrIZGyeOZP55QZNpP+4FB5uNO0PMMnCSfm8GUR3s2P55V5\/V4BaKj4ujAtS1UF0s4hcqvxvRnkUFpIVOaGSRT7TyIet9ELbFHiSyKa2sDSCT08HWW47ZYldPdlqWHt8qvdP4WlvqtCTas8E5+2Y08yHrw+nK\/JavGxWtRJV4Hki1AMY2W7Cmn304hcq8ZE0Ha1olm0A+Cxb0bkRlxueQ+gxOZq1uhX7XgIm6YtLo6bhuhdoi6mZWJqdVLO1BAplNdtHGh6yPaVj5ka0We1avPHAB9+iRo3gWnAbwMoBHuHt364kwIi2NNrY8r1qAxE+8KZZuucDbjDisCuF3cC+Dd5lDpxFpPukv0\/FzvS3gWZ40qSB\/PW+9rSD3Bg8zFr1OVQsvbM9O+z4U8d22w5WOjAa5YeN0ugCTHeh0QRWsxQKWET5GH8i6IttiJ75ZW5wu9twCwEtQoi7sHaifLoEbABiksEmDjZV7773Xr2i37XANZLcd4dRUhfJ4s8LKJCbr5nTSOLVaoCxui5NJZtQWNzSybsfDTXo4EB08eFD5WEgPRZmnGSVX1APZKSh5UwM7bk7Ndk4nDeIzubedScNtwcDZjQ+zXTBtYXC4BQDyMTqiht+Z1RSfQZjzxBNP0FtvvaXcN+EJJ11WGaBwuDciQYRErlmzJr3++usqre496bYdYZXI7HWWl9eEzpwZRd9+m6SOk0vfCdlgAGXtWn+ONXSyWXxdugSaNDonMvxumeo2MbGT8rnguhhxIvPd1v36\/UTr139PiYmdA\/puM9On000dvU1uAYD8mOsY\/h7yaBMOHwB4J0+eVIz1WJRKplOZjlUHSGLJXA\/g5+Xlqb+6desSDhMYAd1tO8IGZCMe4lWrbqB77qlPR48O9OM6HjhwGr35ZmVauPAPql+\/D+3cudbntOSUmgoqA6SglJh6XeCh9803edSgwQB68cVsKlu2NS1ZkuOTUDonMhyCWrXqTnFxLWnNmlk0bFg6ZWVVpd2736P162\/zcT2wWgQvPqgfOoN8MK91KwCwPR+LXDNLpc6WJF10JVjlVQoQOCyFcUBBOkWhHZj0rVq1onnz5tG1115Lp06dKp5ANpIsWCg0a9aZxo3bQ40ataXSpcuoo09166ZS\/\/4V6Ysvpp\/3KuvgmhuZ1RNIQ30Bh2dwGV2+PIHmzdtJy5ZNVoz1mZkVaOfO\/6Hp06v46X+I37PnZtq+vRHNnfsdvfji0z6e5r59p1Dp0t3p3Dn\/SxDdqBJOJXLqn1YaANns\/nkJVpk\/JglOS19yySXqjJ6ZROY0TEIppTrO++EsHyY\/FnvFTiLrAyK5kYuKUlaS+sn6lC27lO64o5TqeFkXAPPUqfZUvvy1vuhwWzx2bDrNnNlO3Wbkhqc5lBIZeZsBWVeV+B4QHD+DimWmI+OIFk6eTJs2rdBlOTw54OqqA1medfRK1w+baiEHLRK4kbk+TurCAGB2fGl+Qn5OCcdDAeRg8jRLo6sWdvLm68kCXYMh1zlBnxBxs7NnpyFWcZzwECOvUHIjO6lLpHAiW+nIVv3v5HkwQLabv9t2hFUiRxI3spO6YHAihRPZLQDsAi3U8dy2I2xAliYz7qRwcSM7rUskcSK7BUCoAWo3f7ftCBuQ7TYwFi9wD7gFQKT0r9t2xIAcKSMZZD3cAiDIYj1P5rYdMSB7PiRFm6EZAHTTIr5jS1+3m4frEnW9l2JALlrcRFxpRgDg7W+4ukrgmu1qhuMS9ZAAWV+8BOJ383okjWi0zMrAoPEtm2YkMV7XL9LzMwJyIO87PENghvpwXaLuOZAHDBiQL2lKi5rj2C6Q7caLdOB5XT8jIPNBVrOy5HMj1aIoLlH3HMigzDLbTQnEd8Zccddcc43yjmK3Pvzer18\/VU+W7HpnSwZNxGPKWuzV6x3LnSqdUphHDPliC1US8PHN9VlZWYqelv0DmIDP7G2jl2vEJQfPLll30G6hHNRt\/Pjxqs1jxoxRfMK6m6PXAOb8vACykTdaqC9R9xzIoMzSSbu5EDtARlykZ644tgXPnDlTDSbAJLmIdTBIIDNDJm\/5Qn3gPHTA61zOuGRdAhfxwcOGejDgJZ+ddPzW2xkojQ5ksHPq7J6YBMx9zJSxRQnk1FRcP2xeoi6Rw3GJesQBmZk19Ve\/HHC7QGbvKSkN2KUwEJDxDFJRcjkbkZJLzmcuC2ml04xkrsezQG8TlsiYQCDmNiNo5BuSQgFmI4kciC\/ESkeWddQ3iry8RL1YA5nZ8H\/66SclXSVQAgEZ6gQTT8PTCm8IyXpvxvkcaND0XUZWNYxUi0gDMtqF0y1wF+VFHX4zslqEYnIFk6dr85tbHdlLiQy3SXmzTyBdWrdaoCOgXsBtEFcCGLHe2+lged+HlEDRBmS0FZIZpjgORv7XdvqkKOK4BrKV1cLoXpFAjO1MWm0GQtaD9btCkI6veJBqAUtnK4msA5BVB7N7UfRL3K3u+UCd9LpHqmpRFMDzugzXQIYbp25HlpfXSH0Vr1uY6gACs6sHjIDM+iM4kaHz3n777bRx40a\/1z+nk3Vp3rw5bdq0ye+qLSz+UA+woH\/66ad+t0vJxSFfJhPopio5GPoqXb\/Ax6juMSB7B2dPgOxddWI5FXUPuAVAUdfXrDy37Yj5WkTKSAZZDzMAYLFnRJeABZ9OVRCOS9Q9t1q4PSEiTTRFubUd5LgXu2RmQMZCz4gZSf8d6cNxiXpEAdnpUfxih6IIaJBbIGNtsHjxYnVTaVFeoh5RQEYn4rUEuy3fzca7cqgoSDn024XsjH1ubjkV7eDBeGU+ghdXWloepaVduKlo5cpK1KjR8UL3gCDtli3xKk16eh4lJvrfbqSXj6u44uPj7VTLNI4VWaOrzC0S87VeuFINtKwcwJ0BO7IejH5HHjBXYuMKO5I4NJucnEwrV65UV5wdP35cXcELKiyAHhtPy5Ytow8\/\/NAv+2bNmqnrg3HtGNwW1qxZQ+3bt\/e7hsysOdwOsIoijdPgSkcGkM2Ov2NwZ8w4TkeOlKeGDY8rUNasmWf4Pzl5rboxHn8cHnusp\/qIG5cQwFiJfC67bBd98sk1dPjwYQXknj39OxNxOS0+o8yaNU9To0Y\/q\/+pqdm+MnDj0E033aQGK9hQnPiJg+0DL9MFy5EcUiBjwbFrlz85iVGjU1PrFvoZNFplyuz3ARkRcnNfVvFAX1Whwgq\/Z0b5AvwnT7Yn5IXP+IuP36L+X3bZbSoJJBn+gg38ao92fuJg2+9lOjccya6BHEi14FcugAPJKv\/LDmCp62WnmAEbqsr336eou0sAajd0\/yjDrdko1G2Opvzd9KUrIMcWeyUTyBh3uM7iYkq5QJREhsFMoLABWUokfA6WJSaYRkdKGjedHylt8KoeUQ1krzohWvOJdiCj\/uB4Y1cAHCWDjzUC7wvIDRP2s2GJDGsGHL3Y9QCWDeQBBzC4MvBRKviywBoCFwPJ9inH3U1fulItAoHPiG84EPcx8tK51HSJb2fDxUm5ZjzHKNdst0tvs9750MGlx1mgPjLbebM7qYNJDw84eXWx3BDRbcrog0aNGtHnn3\/up0awatGlSxeaMmUKjR071s\/8agRkeCUi6BzMkhMu4oBsxDcMn+ChQ3GB+Gg\/7uPMzEzFHyy911hFcaqDOy3XiOcYNlL9pEogYOmdj52zcePsQTE\/v3C8uDh7aRErmPSom9zxQ\/35QK\/uSM9SGf9ZSkuJjHtAIDh4M4VVCyMg\/\/LLLwTHKw5GUjnigGwktWBvbteuHU2fPr3QBgpoTHG8yAvuY\/1VFahc\/Z4LJjGEgz+uGhg5cqQ6+eEEyNEokRnIgXb5+I05ceJEatu2La1bt46kRObDDS1atFATg6UvS\/w6f74G5Mkcq7eb\/SlNFDLVQlYiHNzHrB6gA3WeYyO6Vyn9+WgWt0G6dHrZ+U4GKlRxpURGGVI1w\/c5c+YoHRq6MIKZjoxncLvt06eP77wi3sIQUmCulzoy4sr7R7htES2RnfANG1HGBkuc7aRcI92cO5efYfvVSJq46fxQgTNa83XTlyGVyE74htH5XnEfOynXiOdYBwLyQzA6ROqm86MVcKGqt5u+DBmQnfAN88pVB7LTxR6rE9DT9KNMrPMijrx2y4jnWK7k9ToVN9UiVKAMJt+IA7LR6rcouI+dlhuI51i3nZod6ecyde+zYAaypKdhD7hgNtZCJpFLyqDEvN+8HemweL9524TozS2c\/sjR22vGNQ\/WiSsmkYsbEkpoe2JALqEDX9yaHQNycRvREtqeGJBL6MAXt2bHgFzcRrSEticG5BI68MWt2f8PafToc3m9BakAAAAASUVORK5CYII=","height":58,"width":97}}
%---
%[output:6c0abbfa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:5e11db53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:1f8ba7d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: The assignment added rows to the table, but did not assign values to all of the table's existing variables. Those variables are extended with rows containing default values."}}
%---
%[output:8d83e427]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2b313d91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:557e7fe9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4b5e15b2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5773662c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:1d59e938]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9ade79aa]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6bb7324c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4f4fe050]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7f604946]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6ea0556c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0eb45221]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:04afdb29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:458de9dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:795f3b3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:77b0daf7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:821b9aec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:02fbbbdf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:898fb1a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:892e4de9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5562f0fb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:19210897]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a1cf6f2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9d95dc6f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b903a60]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8bb61f75]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:55511abb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:45ba36de]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2532efe0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c02a654]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:139ec4a6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:251f6ee6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2bfac671]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3f70d40b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:408c013c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6a88a66c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:32712235]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0ae0c948]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8a2501cf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4d40bc11]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5c3d6a53]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5d69597c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8855113b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:946682dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8048bf8e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8ff0b9a2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6c1dd886]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:68537892]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:60606612]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4098bad7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:86521e56]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:206ff973]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:29894092]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:6e5d61ba]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5bfa5163]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:53b6b0d5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3654cece]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47973b33]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:674b4041]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9c708581]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97b7ffc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:27196540]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7d0809ec]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:590930b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2349716f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:995f949d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0deb8e24]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:41231d7e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:33fb96a5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5f29d1a0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:88e9a44c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3b78f7b1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:84ec0601]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34b9deb0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7c6038e6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:73af79b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:16c46dc7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:82b652f5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:879494cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:022d8ae5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:020ce108]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:32dbbaae]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95daa554]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:567b5e80]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7834d797]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:05d81d15]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8d3c3037]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9b8f104d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:279a5baf]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:58aae59b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:35440d17]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3fbf2f6c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1307ab25]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5540db3a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7b3871fc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2f54b455]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:24eec615]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6645ab9d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:39069aa2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0954c4dc]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:483a055b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3516ade9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:75f2c6bb]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:0c619943]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:111a2789]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7c4e216c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2275369c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:9a7a5115]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:06b549ad]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:05a760f6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:4f5c3002]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8674442c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:19970e91]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3994f613]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:35d371dd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7947507c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7ee5a05c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:80435562]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f7b889d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:8625f8f7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:27e94020]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:47edea7b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:083d465c]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:97d8ea97]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c8d35a7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14c7b5b8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:586fed4b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e23854f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:90395e18]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:69efb019]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:14bc5711]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:130df383]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8c8be2ca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5f3166b0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3592293e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:2368f304]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:705a29e1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:829fd44e]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:196b0863]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:96bda9e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:14114271]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2c4909e5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:7b28ec87]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1486087a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:34417615]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:73c3f070]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:8b51281f]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0ecc1f28]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:1ed99c03]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:2b415f7a]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:95137503]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:23a2b67d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:49b372d0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3c4fcc60]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:38fd57b3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3ba6c0cd]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03023ca5]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0505da1b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:9159c2b9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5b1a0de6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:54292d7b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:0a6fc3df]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:01532a58]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:47469dac]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:03139ea1]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:4c1bd746]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:85df4659]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74ebdf3b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:7e1c809b]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:012edad7]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6fa4dace]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:61b7ef84]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:44d0f9c9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:911b8834]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:468e70d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3ac308da]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:5fa7c808]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:84d765d2]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3552ebe9]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:12b850f0]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:860be683]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:74e1c2d8]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:5cdba6c3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:116f9b59]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:167f3238]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:3c8c1774]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:6a2785e4]
%   data: {"dataType":"warning","outputData":{"text":"Warning: more than a third of data is NaN! autocorrelation is not reliable"}}
%---
%[output:805f4cd3]
%   data: {"dataType":"warning","outputData":{"text":"Warning: the trends for the temporal aggregation are not homogeneous"}}
%---
%[output:3e0c14c0]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","ss","slope","UCL","LCL"],"columns":12,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell"],"header":"435×12 table","name":"MCOH_result_MK","rows":435,"type":"table","value":[["'MCOH'","2021","10","'daily'","'BsG_S'","'neph'","'y'","'MK'","95","-0.4082","-0.1404","-0.6807"],["'MCOH'","2021","10","'daily'","'BsG_S'","'neph'","'MetSea'","'MK'","[-1;-1;-1;0;-1]","[-3.0642;0.2017;0.1492;-0.2732;NaN]","[-1.0151;0.6041;0.7941;1.4522;NaN]","[-5.3050;-0.2185;-0.5134;-2.0698;NaN]"],["'MCOH'","2021","10","'daily'","'BsG_S'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MCOH'","2015","10","'daily'","'BsG_S'","'neph'","'y'","'MK'","0","0.1113","0.3370","-0.1116"],["'MCOH'","2015","10","'daily'","'BsG_S'","'neph'","'MetSea'","'MK'","[0;-1;0;0;-1]","[-0.1337;0.1513;0.7799;-1.1294;NaN]","[0.6708;0.5573;1.5057;0.6869;NaN]","[-0.9465;-0.2447;0.0547;-2.8937;NaN]"],["'MCOH'","2015","10","'daily'","'BsG_S'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MCOH'","2014","10","'daily'","'BsG_S'","'neph'","'y'","'MK'","95","0.2614","0.4557","0.0666"],["'MCOH'","2014","10","'daily'","'BsG_S'","'neph'","'MetSea'","'MK'","[-2;-1;-1;0;-1]","[0.6454;-0.0393;0.3857;0.6676;NaN]","[1.3634;0.2950;0.9813;2.5342;NaN]","[-0.0723;-0.3854;-0.2094;-1.1371;NaN]"],["'MCOH'","2014","10","'daily'","'BsG_S'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MCOH'","2021","10","'daily'","'BbsG_S'","'neph'","'y'","'MK'","95","-0.0677","-0.0080","-0.1277"],["'MCOH'","2021","10","'daily'","'BbsG_S'","'neph'","'MetSea'","'MK'","[-1;-1;-1;-1;-1]","[-0.1242;-0.0436;-0.0370;-0.1044;NaN]","[0.2790;0.0705;0.1296;0.2113;NaN]","[-0.6109;-0.1620;-0.2052;-0.4192;NaN]"],["'MCOH'","2021","10","'daily'","'BbsG_S'","'neph'","'month'","'MK'","13×1 double","13×1 double","13×1 double","13×1 double"],["'MCOH'","2015","10","'daily'","'BsG_S'","'neph'","'y'","'MK'","0","0.1113","0.3370","-0.1116"],["'MCOH'","2015","10","'daily'","'BsG_S'","'neph'","'MetSea'","'MK'","[0;-1;0;0;-1]","[-0.1337;0.1513;0.7799;-1.1294;NaN]","[0.6708;0.5573;1.5057;0.6869;NaN]","[-0.9465;-0.2447;0.0547;-2.8937;NaN]"]]}}
%---
%[output:4aa06cfe]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"10×19 table","name":"MCOH_result_LMSlog","rows":10,"type":"table","value":[["'MCOH'","2021","10","'month'","'BsG_S'","'neph'","'log'","'LMS'","1.9755","90","-0.0189","2.3484e-04","-0.0381","-0.4698","0.0058","-0.9454","-0.1725","-0.1721","-0.1729"],["'MCOH'","2021","10","'month'","'BbsG_S'","'neph'","'log'","'LMS'","5.4097","95","-0.0381","-0.0240","-0.0522","-1.8807","-1.1854","-2.5759","-0.3169","-0.3166","-0.3172"],["'MCOH'","2021","10","'month'","'U_S'","'RH'","'log'","'LMS'","0.2621","0","0.0042","0.0359","-0.0276","0.0963","0.8314","-0.6387","0.0425","0.0434","0.0416"],["'MCOH'","2025","10","'month'","'Ba3_A'","'abs'","'log'","'LMS'","0.4513","0","-0.0077","0.0265","-0.0420","-1.8065","6.1984","-9.8114","-0.0744","-0.0736","-0.0753"],["'MCOH'","2025","20","'month'","'Ba3_A'","'abs'","'log'","'LMS'","1.5678","0","-0.0220","0.0061","-0.0500","-3.6642","1.0102","-8.3387","-0.3554","-0.3544","-0.3564"],["'MCOH'","2021","10","'month'","'BbsFG'","'neph'","'log'","'LMS'","3.5776","95","-0.0245","-0.0108","-0.0382","-1.1903","-0.5249","-1.8557","-0.2172","-0.2169","-0.2175"],["'MCOH'","2021","10","'month'","'expS_bg3'","'neph'","'log'","'LMS'","3.0143","95","0.0175","0.0291","0.0059","69.5646","115.7204","23.4088","0.1911","0.1915","0.1907"],["'MCOH'","2025","10","'month'","'expA_fit'","'abs'","'log'","'LMS'","2.2098","95","-0.0139","-0.0013","-0.0264","-9.7036","-0.9214","-18.4859","-0.1294","-0.1291","-0.1296"],["'MCOH'","2025","20","'month'","'expA_fit'","'abs'","'log'","'LMS'","1.4405","0","-0.0079","0.0031","-0.0188","-5.0200","1.9498","-11.9898","-0.1458","-0.1452","-0.1463"],["'MCOH'","2021","10","'month'","'SSA311'","'abs+neph'","'log'","'LMS'","3.2769","95","0.0028","0.0045","0.0011","8.3263","13.4082","3.2445","0.0283","0.0284","0.0283"]]}}
%---
%[output:51e91622]
%   data: {"dataType":"tabular","outputData":{"columnNames":["station","end_time","length_period","granularity","parameter","instrument","MK_seasonality","method","significance","ss","slope","UCL","LCL","slopeP","UCLP","LCLP","slopeR","UCLR","LCLR"],"columns":19,"dataTypes":["cellstr","double","double","cellstr","cellstr","cellstr","cellstr","cellstr","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell","cell"],"header":"10×19 table","name":"MCOH_resultLMSlin","rows":10,"type":"table","value":[["'MCOH'","2021","10","'month'","'BsG_S'","'neph'","'lin'","'LMS'","2.5380","95","-2.0121","-0.4265","-3.5976","-3.5291","-0.7481","-6.3100","-21.1209","-21.0775","-21.1643"],["'MCOH'","2021","10","'month'","'BbsG_S'","'neph'","'lin'","'LMS'","4.2797","95","-0.3384","-0.1803","-0.4965","-4.4594","-2.3754","-6.5434","-4.3839","-4.3796","-4.3882"],["'MCOH'","2021","10","'month'","'U_S'","'RH'","'lin'","'LMS'","0.3854","0","0.4573","2.8304","-1.9159","0.6066","3.7550","-2.5418","3.5726","3.6375","3.5076"],["'MCOH'","2025","10","'month'","'Ba3_A'","'abs'","'lin'","'LMS'","0.1089","0","0.0036","0.0694","-0.0622","0.2334","4.5219","-4.0551","-0.9642","-0.9624","-0.9660"],["'MCOH'","2025","20","'month'","'Ba3_A'","'abs'","'lin'","'LMS'","3.3148","95","-0.0829","-0.0329","-0.1329","-4.5533","-1.8061","-7.3005","-2.6579","-2.6552","-2.6607"],["'MCOH'","2021","10","'month'","'BbsFG'","'neph'","'lin'","'LMS'","3.2848","95","-0.0029","-0.0011","-0.0047","-2.2772","-0.8907","-3.6638","-1.0291","-1.0291","-1.0292"],["'MCOH'","2021","10","'month'","'expS_bg3'","'neph'","'lin'","'LMS'","4.5602","95","0.0195","0.0281","0.0110","2.0019","2.8799","1.1239","-0.8048","-0.8045","-0.8050"],["'MCOH'","2025","10","'month'","'expA_fit'","'abs'","'lin'","'LMS'","2.2473","95","-0.0157","-0.0017","-0.0298","-1.3650","-0.1502","-2.5798","-1.1574","-1.1571","-1.1578"],["'MCOH'","2025","20","'month'","'expA_fit'","'abs'","'lin'","'LMS'","3.0394","95","-0.0098","-0.0034","-0.0162","-0.8374","-0.2864","-1.3885","-1.1959","-1.1956","-1.1963"],["'MCOH'","2021","10","'month'","'SSA311'","'abs+neph'","'lin'","'LMS'","3.2517","95","0.0027","0.0043","0.0010","0.2761","0.4459","0.1063","-0.9733","-0.9733","-0.9733"]]}}
%---
%[output:0f4d311e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAHkAAABJCAYAAAAZgInkAAAAAXNSR0IArs4c6QAADptJREFUeF7tXWloFNsSrpCIccUVt0QTNV55ICouiSLP+CO890RFjBqNV40bwSsSee7GJYqJGy7RHy5RE8WV6wL6rgr+SHy4JCI8xR\/+UIyXO3FBoyIuEQfy8p3kjGd6uvucnu6emYQckMtNnz5dp+osVfVV1UTV1tbWUnNr0hyIahZyk5Yvm5ySkMvLyykzM5O9sG3bNsrIyGj6nBFm+OLFC\/Z\/CQkJjXLeUiG\/f\/+eli9fTuvWrWMTLCgooN27d1OnTp0a5YStEF1WVkabN28m\/Jc3CLq4uJhSU1OtDBXWvlIhYxfv2LGDjh07Rq1ataI1a9bQzJkzKSUlJayEu\/1xCDcvL8\/wMxB0VlaW22Q4Mr6SkM+ePUvbt29nH4SQR48e7TuyvV4v4V9Tardv36a0tDTplG7evEljxoyR9nO6Q0xMDOGfarMlZI\/HQ\/Hx8arfau7nIAesGEVKQjY6rnGUY8V\/\/vw5gHwcZWPHjvX9vaamhvXr0KGDpVXoIF8MhxJpq6qqMj2mtYPg2Ha76fHOylUhFbKZ4mXl3gKhr169oh49elBsbKzbfLE0vkjb69evKTExUfn9yspK17Vuu7yLqq6url2wYAE9evSIJk6cyO5eKFi8ffv2jRYuXEj37t1jf5owYQLt37+faZzjxo2TMqO0tJRponYJlX7IRgctbVFRUcqjWTk2lQfVdLTLu6jt27fXwiyYNGmSruZ8\/vx5gp24evVqwq7GgpgxYwYdOnTIz7QwmgCOFRxpz5552U4eMaJrRO9knDJYvKLZZDQ3LF4sYrebbSFPnjy5FgKESSQK1Ihw3M9t2rShpUuXKs8tNbWWBFOT4FPAVRYppqaWiSqnVCjtZdtCzsrKqoWjIykpiQn57t27AUc2lya\/n+fOnat0VP9cBZXwFwUsiqIiL\/36a\/jNLzARd3G3bt18V1VJSQktXrzYcCFjY5jZ0co7QKGjlj7LJpSqkHE3i44QK\/cWkTEGcubMK0pJqVGYqntdvn\/\/Tm\/fvqWuXbtSy5YtfR+C9VBYWEgVFRW+v8XFxVFOTg6lp6e7R5BmZC19sFA6duyo\/P0oleNa1LCx49FU7y0iuP+M761Zs7x09Gh4dzMW8Js3b6h79+6G+gJ8AhBwOJqWPss7WaZ4QcAwlTZt2uTnr1a5t+qPaNiR5n7ee\/fKw8E73zfhsXv37h2z4SPNvAORWvqw2KwsOD8TKjs726dFc8EWFRXR4cOH\/YTAkSjcW\/PmzWPP+I1bj9fwBt\/vJqkA4+P\/TjExHmm\/5g71HMB1gX+qTeoMkQ1UVufIr9NA\/PYqBD2PsqiM7WJ5C\/dOllMYWT0s72RbQQMNAjZiwTwqphLiSA3f4\/5aNsyoEJiakSWlEFMT\/E6G4avg8RpHeVRGwGN\/YrL8rk5ISI0oeznEvA\/Z54IXMgQsejgMSC5hR7d+S00tptLSxoHJhkwiLnwoOCEjHMaCE9\/MEww\/+KJFi9jUEHGC4AQ0ABlQ+Pr3729r2hcuXKCNGzeyMQYNGkQHDx5kVsL9+\/d9oP+WLVto6tSpBFNlw4YNdO3aNdYfLlxExaDp9bdFGBHZoY0ryaABXkiuHPfs2ZO5kbmpi+dKQg6I8UpOtiRkYDr+WvdP9sCpsGvXLrv8svU+TJSvX79S69atIw4GxcS09MHUmz9\/PsMQjDAHkSFSIRtCjZ07KzNWhulgB4WzqThDIok+OENw6sF\/Djyfx+CJu9eSkA1jvNauZXdynQFl2rCDcS+btVD5gI1o+PHjBwtoaNu2LbVo0SKc8tT9th59QAaBOVRXV\/siafGyeIzzwaQ7GULWjfHq1o1p17JdGnEcayIEAd\/XBlNyfEGMwVO6kw2FnJFBnvx8il+\/vomwrXFNwyhYAUoYGlAySzvZLMbr31Om0PHevSlWQGq8cXH0bNYs+lfDB43YB89Nbm5u2GOYAeXBd92lSxef73rHjlgqMblncnNDB5Nq6cOdzK0ObMJbt275BXXw+ACfkGXhP6Lidf36dTp9+jRdvnyZoKpzrXvPnj00ePBgivF4CALm7eLFi7Ry5UrDLWLVB+vWXtNCeeXlsZSZ2UP6uVDBpHpQ47Bhw3zZLKIJpXsny1AozFQ0oUaNGkVHjx5l4Dr\/+8mTJ2n48OG6TEEMc35+fkAWAlAtmACR0LTadVpaDN2+LY9rDhVMahtqVMGT8REIql+\/fvTw4UNf5AgX8vHjx6UZFThyHjx4QFiBYqBgJAhZjLyorm5HSUlyAXO6v31zP+AhJJEhCAtC69Onj0\/TFnfynDlzaMiQIabygkH\/4cMHat++vV\/0RSQIWaQtOrovTZs2QpmsAwf+Q3Fx7gY9aHk3cOBAwj\/VFiUL\/3n69CmdOHGCKUgI2+XmFISMaAncuQiP4akyRukbsueqBIeiX2Xlc+XPJCb2Ve7rVEeruow0\/Ae7eC0cH0JDfDZin9BOnTrFfM84UniDoLdu3coyKKD5wdlh9NypiTs5ztq1KSrYC4s23bYt9FEtlvFkFcWLM1C0mbGTZRkUiLlG9IhRi9TMQBUUNdLCis0WuTT8R8xDFoWMI1olg0K2w3iGhaxfqJ9jbTZENul+Glmtm+SRTaEmW\/d7UremEZXq0Zrm8+QZFhHBDQ0R2NEIftEmBkC4jSQ1mc0oKCHDOW4lKUwmQO6ik+GisnH0nos6BRw2wKtxOumVyOC+36tXr7KhRMfC5cv\/o+XL62OtnSqp4RRtMr4pCVnLkOTkZEeF\/OTJE9czA80WCDR\/jkJZSe4OZtEF846WPtA4bdo09\/HkzhbwZNnEnj9XN1lkYwXz3CiDIpix3HhHz6155MgR9\/FkmFUqmX+ySQP4Dnc1IeyUL1++MDw5OjpaRnLIn2vpw05GsrzreDKSw5zQrkPOsSbywaDxZJ5\/LCaki14uBJ4hYAwKyezZs5k\/e30znhyWZRM0ngwtTZuQjhng73Brwp2GbEYIF94sKGQIJOvbt69u5l+vXr3YsYLFoZcZiOeR0HDnffr0iWnd4nFdUdGKCgs7BJCYk\/ORkpNDF5empa9du3as4gOaEp7MMyj4LtYmpMNluWzZMgYs9O7dm4ARa\/FkDjWaZf7BPMGpACAj0pLK9AL54OxA4IBRO3iwJmS2sh7UOGDAAHU8WRSyGPUnJqTv3buXYchoeniyKtSoTfSOhF0MGrRQHrBkYMqydvOml8aMcReB0qPPcuqqTMgoCgOB66FQ2sgQM6ZEspmipW3mzB5UUSGvUJSe\/pl27XorWwu2n9tOQpcd17ijjVAoLuQlS5awYACzBjMA4aMIDBez+W1zwIEBRNqio\/tRWpp61kZpqZjj5QAxOkNoeQeZWCnm6ufx0lO8xLBPLQol4snuTC\/0o3q9cfTXX\/9V\/rBKbnVcQ1lKj4VSiWYEWMaTATVCU4YwRRNKLyBMK2QQAkHjX1Nqo0apF4dludV\/\/lk\/\/T59\/NlQh6UjdztFwNohaM\/WrURCtULD9w2YahlPtpWf3JQkK8xFMWGTshLKqDhBB6ZC\/aoGARuyCH0ASuvBXJr6V3brbSsBFE1UlobTCgwaCEygz6PNdYUyjEsl2+ZZXU3tsoQER+ptNwsZGZc6lenrgwagVCHbyz+BPpWy6uoZuSjghhWCopZGap2VqBolIZv9PAE0v3Pnzvkq2w8dOtQ0B5ivcJ4LPH36dBo5ciT7c6jzkxH3jQhSMf4MTh9UHAR9Hz9+ND5t6xaH2+nzZgn8IIxH1djGk82q5ELhAgiPzAqt23Lnzp1+sdjhzgH2eOqdGzx8FoGIPBgx2KMV1cnc\/pECs4RCKMzgvezXBKQ72eznCWSBfHCR8rrMEDI8XmK+UbDMtfIeTmK4KEXnBgTt8aC6gT0bFwJWq29khWL\/vjiyjRL40RO+B0fyk\/VSV5uhxuAF5+Sbly5d8pW8wLiO5ierlkJ2ckLNYwVyQIQcbeUna1NXUed5ypQpFngeqIkiMF0shdyg4DLTUa\/pPYcGzP9ugRhNV5VaCMajG+rYmIgN4kAVLhOzo1qv3rY0P1kvaACrg5\/5PHUVFXus\/WSQcSlk8CI\/PyYg7BVlkoHwABHSe46QWJNf9LEo8+DqJZiVhq1BEjhWcV5eQO62NzeX0RdrUmoZhpuZkcbrbQO2leYnix4vI9+1XuoqKuWot+B+DhIoz8WLbdU\/E3RPs\/pE+oPCtXggPZ3+WV4eIMQPOTn0WVMqWZu7zYRcV2q5Y2Gh7vsniExzuwH04MpEfJw0P1mGQkFD1ktdHT9+vGIgn3kp5KDl4uiLgTsZQjTzyYMvYsmGFlVV9MNGpIve+3fu3GE5ZVrzFLndgIBVm8+E0ta0FoMGrly5wsYTU1fV0mTUSiGrEutOP9x8gT+YAqQH8eWRUtT88ePHrNgcYFrLRc3tBA3IA\/nUSiG7Izz5qCi\/HBODjMwbvs6IusAu0YYJh7OouWN4slGMl1nQACaOHxy5ceNGQGqq11sUAsefXJBmPTp0KKwr81+fggtGRmL2hB79lvFkFcWLf0iGJ4uMUs3xtScm87ehuUMLB6qkl7Q2cGDoc4udmK9lPNlu0IAR0So5vk5M2GwMbXopzNZG+hPItlgl9V3bGV2W42tnbLwr26mNKb3ULi\/M3ndVyPiwWY4vnpslekNIZgXTmneq2tJwRMiqebbJyRm+45L7WX\/5JZt+\/\/1vfncmtN7ExBP0xx\/TqaoqyVYiuCpt0KbN8pPNMHU1Vgf2coo223iybAKoDlRQUMAAf6SZ4IMvX76kVatWsTwpZN6hiX1E9+mZM2eYixQM\/u23nZSd\/Q9Dl6nVO9UqbfjRa73fpRTLDWvnIuOP0XOnaHOk3rXVSXANHADGvn37WFABishwYBvZ\/vB9o1L8ihUrmNeIR4rKcFGrtGj7y2jT+uO5mxdOIKP6onZp0louRnwzos2RetdWJyEyRreEckYGG1Jrl4vHIZ7r4aJWadH2V6WN08cXHYB5s7nYpQvv26EtMzPTR4Ie3\/4P7t3TvEIo2kUAAAAASUVORK5CYII=","height":191,"width":318}}
%---
%[output:1bdfd500]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAHkAAABJCAYAAAAZgInkAAAAAXNSR0IArs4c6QAAEpVJREFUeF7tnQuUV9MXx\/ePiSKUpJipplatsCpvTQ9MJVap6a+HoZZHNPSXymtKROUvPdAySZOkCYUirKiFtJqsYsYzi6U10YMmIj0oGkzdf58znd\/c3537\/D1mppq91qypueecu8\/ZZ5\/H3t+9b8gwDEPiQAsXLpRXX31VsrKy5Ouvv5bRo0dHtMrzjz76SCZPnix16tQR8\/\/37dsn9957rzzwwAPSqlWrOHBT04R5BELxEPLOnTvl7bfflptuukm1jQChzMzM8LusQv7uu+\/ksccekyeffFLWr18vr7zySngC1IgoviMQFyErljZvLuMsNdWWw3dmzJDPP\/9cRs2cqTRZT4YxY8bImWeeKXl5eTVaHF\/ZhlsLLV++3OjWrZuEQiHXV7Ck3n\/\/\/XL99ddLWlpaedn8fJEJE0T4rQlB5+WJpKeX\/d3teYI6VtNs+QiEOnXqZDz\/\/PPSunVrx3FhOb711lvlq6++kpdffrlcyAhv\/Hjn8bz5ZpF585yfMxEoU0MJHQEl5B49eqiDUlJSUoWXocHTp0+X\/v37y3333afKKU1GQ7t0iZ25lSvLNL6GEjYCoYcfftjYtGmTzJgxQ04++WRPbQ4LGQGbl+goWSwdNEhK58yJsvbRWQ1ltFNIp9EITZ482SgoKBCW7FNPPdWfkBs3FmnePG4jvGnjxri1FU1DpaWl8tdff8kJJ5wQaPCieVc0daz81atXT+rXr++7qVB2draxe\/dumTZtmtStW9e2IpOAAxfEPfiBgQPjKuSSdescT+W+exJDQbakX375RRo3biy1a9eOoaXEVLXyF1iT27VrZwwZMkTuvPNO2xM2hy4MFcOGDZOHHnpIjj\/+eJk7d66c2qBB\/HpkY485cOCAfPDBB7Jq1SrhzNCpUyfPG0C0DJWUlMjPP\/8sZ5xxRrUUcqz8hbp3727Mnj1bUh3ut2jxlClT5IknnpB77rlHafvw4cMlbcwYtSe7nJ39jTnvHTeuQtkdO3ZIcXGxtGvXTtasWSNt2rQRlqlEEMshqxl9q46abMffzQFuJaENGzYYLVq0cBw7hKytURTirtyxY0fJbNRIna7db9flzbICnHTSSYLwtCX1xBNPlNNPP10V+umnn+Tvv\/+uwMdxxx2n9kqEUEPlIxDEGu1p8XIUcmamFE+cKE3GjvUcey3MPXv2hIXMvsJB77fffpNatWqFJwBCxTCDwOnIMcccozSMfenff\/\/1fNfRUiDuQma5fuaZZ5StedeuXWXLdVqaMAH+26OHZLVtK28WFcmVO3dKvf37ZUnTptI8I0MuuOAC2VxQIDs\/\/VSkqEg21akjGTt2SK1mzeSL3r2Fq1tGRoYS3ptvvilXXHGFEjbEhGBfZp\/csGGDKnPWWWclRIbseXv37lXbAZNv1Sp3G45mghXz8ssTwlJEo1b+eBhoufZyUHDwGjFihNK4oqIiOe+888LXLYQ8cOBAZQVr0KCBDL7hBvnp118jrWKH2J0\/f77kL10qj8+cqY7\/OCw2b96sjCvaoha+g4uoQZ81a5Yq+8MPP8hdd93lesWLZaitBxu\/JgCEjNEu0RTzwWvJkiVGr169HE+uLJMM\/jvvvKP6MnbsWLnlllvUv7WQR44cKe3bt1d7N5Ph9ttvF5ZdM3FS\/uKLL2To0KHK6LJy5UrZtm2bupqhtVOnTlX\/PueccyLqcegIcvGPZsB5B3yXaXJLufTSFN\/NfPxxge+y0RY088fBMCUlRf34pRBXqGeffTbS6WBT207bOP1mZ2dLYWGhWlrZS\/lhv7A6PMzPad5aXr\/Sy1Hit2PRlistTZEtWz70Xb1Jk8skKanYd\/l4FESp+PFLoZ49expt27aVRx55pIL2mRuxEzLPEfS6desUIOCqq65SVdBSqIvJts3fvvnmm7CWU2\/BggVyxx13yJYtW2TFihW2K4DfjsSzXIcOJi8b3lMpc6Nulopu1MrQZGvfAmuy1ayJM3\/w4MHqStO7d++wI99JyH4H1woaoB5\/q47+ZL0np0u+jJMJwm9NCHqw5En+wb\/iVzk0n\/0OQ5WUU7brzz77TDCIuNlDYxVylfQuypfid1nVZYKME2c36oTUPLk87+bDwoEW6t+\/v9G1a1d1IHLbD48mIft2ox4mbtLQkCFDjHHjximrkvnwc8oppyhDhB\/Syy5lzz333ApXLP4+adKkCMyXHdKE+ziHQChekKAgvIV5mj1b0kpKvLse4x0qGt7A0kHcYDRY0mvcQs2bN6+A1jQLSnfcrnFeZgbkYcHihezno0aNUtctEJgQvx988EF1z7ZDmjjCi7yH2rGEX940oHD58uWy68svZeiUKf7fGiXYNRrerHaF6667ThmTbGFZph6Epk2bVkHIaPG1116rzIl2RgsaB4nJCXnx4sVy9dVXO97bVq9eLRMnTpR8E8AAZ8hzzz0nnTt39j+YCSzJPRTji3JQbNsmSQFgwZXhJjXzp92MoFwZx+7du3vCmT1t19axRVNpHCFrYwi+aLTfSjk5OcKPEz3++OPSr1+\/BIrPX9PYybdv3y4NGzZUrtTmLg4ba4uVAXiw8ofRBo8gqyMOH6yOmszLeHj7xayJtQebMR10I+1b1iB4LWRO5hdeeGFEVdyDXMG8iG0AX7GVMJbgS+b+jU37kksuSZg\/mUEENNCoUSM1BidlZEjS6tVerEtp586yZ8kSz3KxFrDyB6R57dq1FQxYestTXkIT5j20aNEigyV35syZrrZhuz3TrMnYriFtbsNEiSXMi9BkNNpKGEgwsiBgIjMwrOCsSARZNaV2QYGcYdIOu3eWpqTI9qlTpcQMT04EcyLKI2deadzgP6y0kDmCRR280JLc3FzHe\/L777+vTnMQYS56liDkvn37qtnPoUATgma\/NpOXP5lT9UUXXRRRZ\/\/+\/WrGYvNmScJlmQiyagrvqLVggdQdPtzxdSWjR8s+SyhQInijTTtN1uAGZMCK5+Toob5CawL9YT+yI9x8DDCuRg5k+iTKSXrCQdz1eDfc9aEGvfzJTDC0Hnssk+Wff\/4RgAy8D7DAsmXL5LLLLpOmTZsmZBytmqJfgkbXz8mR2qYVCQ3eNXKk7K3Es4SdJrM96mup+QpluycfOHDAcDOCMPD6+qQ7D5ivQ4cOYds09+nTTjtNXY3QPpZulvc\/\/\/xTabn2IrGXaGQIgtf\/X7RokcybN0+uueYaJVjo999\/V20hWJZtOsr1KxFk9fLYvSOpuFgQcFWQlT\/Gk23ML3merp2QIfh6zdciXgzaEZSHHZTHLFS8VPwfdyQgBCYJhx4miRkChJkVQdMm5TiM1VDZCMQdGWLFeBFeCjrESgjl2GOPjcBx6TLRCLlGoM4jEEjITz\/9tEGck440tDar0ZqA7ymDdSU9PV0duMyEEAED\/PHHH2ENtT63LtcszRq3xbJuB+SrEXTFEeBGwjnGL6nTNes7FhSNrzJXNt+N+bs+eOkrE39jucX2jdUI0qdg9mSISfHJJ5+o\/ZmDVLNmzZQfmTuw9ie\/9dZbMmDAgJrwVR+SC+xPRpM5gnOFcQqT0fdh3q+jGrm3mvdkJ96wjvFjNWsSj4wB4u6771YhKhy0NFGe50wOTfqK5oQP9zE2R22REKdrBpBDk9OSbTc6CM2M\/IhmBEEccqp2IgSNULmq2U0S8ySI5v1HSx3P0zUDYdZkfTeLh5BjHWQmQRBoaqzvO1zrqwgKTsVO\/mOnPZnN389yneiBATuGRjtls4iXz9ZuosfaN7+8tW+fKcByc3NHh20Wgf3JwGBffPFF2z3Z7+k61g5HWz89nUwFebbZLJKTy5PPuPu6Zytfd1HRexWw4LhVze488+HTLdTXqz9+\/MmFhXUkO3uvbNtWHlRA6Nj06X\/I3Lk3iW9\/cmpqqoHhn1nVpEmTCrzZGUOc7sleHat5Hr8RGDBgqfTrt8efP5krFPdbTJd+hYwri1lUQ1U5Aoa89FKxtGxZ7O1P7tixo4F5EU3WEYZm1s3L9dKlS9VVC6QHjonK2JMxaeJxIcrCmRKTXKZz538lOblUtm4ty6XCv+3I67kT32vWJElxca0oZwqHzsgwHUd\/8qOPPmrgznOC5OqD16BBg9S+BRGcRoK1WK9QfnqH7xSbNZY0e+IuXQbmjze1b18G5issLM8+kJJSKlOnbpe0tBIpKKh9EPlS3\/G55qe4uGySUFdTTk49VTdWskLMbP3JXbt2NcBakUXAKeYIYwmAewjv05w5c9SdGo0GrOdE5nuwnT\/5yiuvVEBAyM6pgRUNBwV17YVMRAMRZ5WfPSgzs0QWLnROPZGbW6Ly1uGJtU6SQYNKZcqU+KSteOONL2X9+vfd\/cnEQmGQOP\/88x2FpdMoYo40p0fUoAHcgWYUCGY3tB4hs6SD3MTkqeOTk5OTFQiB1QC0IaZSkCC4Idk68G3feOONyhJGeikAdvZCBvxeMUtBrNpxONVHkz39yYsXLzbw4zr5lNG0F154QQmNZG1WIQMoII6KgxvtAOgj9QSxyT179lTIDrQURwQCZy\/nTg4QAMQHucFwTGBk6dOnTziqcePGjep9xDAzgfitnRlMouTkB6WwcOjhJI+48+o3TKeCxcsaCwUkB\/C9mQDogcI0A\/kwqOB+JBzVDsrz2muvCfBcJgRCxiHx448\/qthnnBb85scKAeLZu+++q07zW7duFVYBaM2aWtK7t322oriPZjVs0JzZ0os9X2ZN3Yj5zsyerIXM0gpqgxhmDmoIxBqfDHrz22+\/Vc\/wdlEXlCix0Szl7PPgt4OkQn799bqyYIF\/hITXYBxOz9nrbfLp2HYhJiFXh\/jkkpL2B1EjI6WkpDzclHjhunUXy+7d\/mN4q6OA6Qc\/1r4NGVIsubmR4bVu\/AcSsl1D1Sk+mQByc0D4ihUpMnasMy6LeyY2b2uCX5bCiL+5xCcnanLAw\/\/+V2bsgMx9C+xPJnSVkBhOzn4D3KLpWFXFJztlYmapM6fCQtg6lZnODesVnxzNOPitE2Q59mpTIUPMAW5eFQ7n52ZBevUjv8sESc93j08evzl6SxsTiknGhLKuJNYJ6MWr1\/PQe++9ZwCeJ\/QlFq+K14u8Mtt71q\/MAj7TPOePXykT8tNthQS7h+xHtpxbNTXIBAw6FDHvybzQzS+a07evjNy1KzLeNzVV\/p41S0YtWxaRCd8rzjZo57x404FiGggRS3yynZC8tgq\/\/mQiVtxCiL3GLZSTk2NgpOAKE40mu\/lFP+3VS\/6zdq2jbLIbNpR+S5aowK0jPT7ZOgn8+JMZOHPsdNTxyRkZGSqCQmW+dch37TRTdHwyKBFQ\/lAYaOdzyVOZVao4cz28Y3Ll\/k6kRKBc3ps2JTyNcwR\/h74uoEOIfcUnc\/AidBQvlB2Qz03DtO0aCO73338f1liFtsRtYAqCc1TnGFMyRLOEW+tUyHjn8dGViPpRZhoIwrddRj4cRr7jk1u3bm1oJ4GdF8oak2xmzivgDUH7OX+W7NsXpM9xL8sgYo4lVIeJntS9u+\/45NLly+POj90kjOAvKUnI2BTxVZ+D4aqO\/uSsrCyDTjp9g8IMYOPlGkDmF62Jp9fLEbjlww+rLJiMPh3x8ckjRowwcBR4fYOCwTDPFGvAm9N0LoPZuVNVa7Ld5wmS5s+XpKwsZ8bHjxdilCuD3D5P4Cs+ediwYQb7KQcv0ip5ERs+EYYageBVnueuH4P06y\/z86IoyzhmofW6A0X5vqDV7Pgjftt3fHKzZs0U7hp\/sF3mVe154jmRdID9CG2xBry5Mb5J5ae0oSD+sqAjE6C8r1TDibRWeK10MX4jI9SiRQv18S8C3uwSw3DwAgwAcgNCyGC8zAFvXuPpqMnxNNB6MeHy3JeQY2g\/1qqx8hdq2bKlSoVM6kU7sgPXk\/SFxKd+0JposC1Ah2X60FdaYx2EWOsf8R8a6datm0HIqE78Yh0wp0wDXDcqA60ZqwCP1PqBgtAnTZpkgKdy+viX24dGAABqFKfdYILn0jmwj9TBrop+ReVPRpBOVyin5VpfxFmy7UJLwYXVRBxWxRSo+E71GT+QkF7gep0I1ZziydocBvSaIPHqIVgzF6FWrVoZfKbe6dO6FLbLNFD9ulLDkdMIhNLS0hS43u0j2V7DF8Qvqtuqzvmu7dI+V2V8csz+5IsvvtgAPH\/22Wd7ydL2eVC\/KO7M6pzv2s5nW1XxyQx4XPzJbdq0MW677TbHr64Glbw+jWMRe+qpp9SBTqeG4n4Nnszuy+pu3q6gPDiV9+LN6tXRPltAjvzb2hdr+Vj4jJY3X\/7kPn36GISsun0\/OQjz5oGx\/ZDnoRS91m9aOHm7grzbq6xf3mjHPOmI0XLri9d7\/TyPhTevfNf\/B08p4tojlaWxAAAAAElFTkSuQmCC","height":191,"width":318}}
%---
