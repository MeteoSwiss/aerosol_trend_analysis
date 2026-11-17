FKL_st.name='FKL'; 
FKL_st.lat=35.34;
FKL_st.lon=25.67;
FKL_st.alt=150;
FKL_st.env='coast';
FKL_st.footp='RB';
%%
%2026: read data from import file

names_AE19=fieldnames(AE_19);
names_BC=names_AE19(startsWith(names_AE19,'equivalent'));
names_del=names_AE19(~startsWith(names_AE19,'equivalent'));
for i=1:length(names_del)-3
    AE_19.(names_del{i})=[];
end

names_AE_20_S00=fieldnames(AE_20_S00);
names_BC=names_AE_20_S00(startsWith(names_AE_20_S00,'equivalent'));
names_del=names_AE_20_S00(~startsWith(names_AE_20_S00,'equivalent'));
for i=1:length(names_del)-3
    AE_20_S00.(names_del{i})=[];
end

AE33_S00=[AE_19;AE_20_S00];

%%


names_AE_20_S08=fieldnames(AE_20_S08);
names_BC=names_AE_20_S08(startsWith(names_AE_20_S08,'equivalent'));
names_del=names_AE_20_S08(~startsWith(names_AE_20_S08,'equivalent'));
for i=1:length(names_del)-3
    AE_20_S08.(names_del{i})=[];
end

names_AE_21=fieldnames(AE_21);
names_BC=names_AE_21(startsWith(names_AE_21,'equivalent'));
names_del=names_AE_21(~startsWith(names_AE_21,'equivalent'));
for i=1:length(names_del)-3
    AE_21.(names_del{i})=[];
end

names_AE_22=fieldnames(AE_22);
names_BC=names_AE_22(startsWith(names_AE_22,'equivalent'));
names_del=names_AE_22(~startsWith(names_AE_22,'equivalent'));
for i=1:length(names_del)-3
    AE_22.(names_del{i})=[];
end

names_AE_23=fieldnames(AE_23);
names_BC=names_AE_23(startsWith(names_AE_23,'equivalent'));
names_del=names_AE_23(~startsWith(names_AE_23,'equivalent'));
for i=1:length(names_del)-3
    AE_23.(names_del{i})=[];
end

names_AE_24=fieldnames(AE_24);
names_BC=names_AE_24(startsWith(names_AE_24,'equivalent'));
names_del=names_AE_24(~startsWith(names_AE_24,'equivalent'));
for i=1:length(names_del)-3
    AE_24.(names_del{i})=[];
end

names_AE_25=fieldnames(AE_25);
names_BC=names_AE_25(startsWith(names_AE_25,'equivalent'));
names_del=names_AE_25(~startsWith(names_AE_25,'equivalent'));
for i=1:length(names_del)-3
    AE_25.(names_del{i})=[];
end

AE33_S08=[AE_20_S08;AE_21;AE_22;AE_23;AE_24;AE_25];
%%
AE19_25=outerjoin(AE33_S00,AE33_S08);
AE19_25.Properties.DimensionNames{1}='Time';
AE19_25.Properties.VariableNames{1}='Ba3_A00';
AE19_25.Properties.VariableNames{2}='Ba1_A00';
AE19_25.Properties.VariableNames{3}='Ba2_A00';
AE19_25.Properties.VariableNames{4}='Ba4_A00';
AE19_25.Properties.VariableNames{5}='Ba5_A00';
AE19_25.Properties.VariableNames{6}='Ba6_A00';
AE19_25.Properties.VariableNames{7}='Ba7_A00';


AE19_25.Properties.VariableNames{8}='Ba3_A08';
AE19_25.Properties.VariableNames{9}='Ba1_A08';
AE19_25.Properties.VariableNames{10}='Ba2_A08';
AE19_25.Properties.VariableNames{11}='Ba4_A08';
AE19_25.Properties.VariableNames{12}='Ba5_A08';
AE19_25.Properties.VariableNames{13}='Ba6_A08';
AE19_25.Properties.VariableNames{14}='Ba7_A08';


FKL_rd2=outerjoin(FKL_rd,AE19_25);
FKL_rd2=retime(FKL_rd2,'daily',@nanmedian);
plotFigControl(FKL_rd2,FKL_st.name); %[output:6025f3a3] %[output:1c89a639]
%%
cte = 14625e-9/3.5;

FKL_rd2.Ba1_A00= cte * FKL_rd2.Ba1_A00.*1e6./ 370;
FKL_rd2.Ba2_A00= cte * FKL_rd2.Ba2_A00.*1e6./ 470;
FKL_rd2.Ba3_A00= cte * FKL_rd2.Ba3_A00.*1e6./ 520;
FKL_rd2.Ba4_A00= cte * FKL_rd2.Ba4_A00.*1e6./ 590;
FKL_rd2.Ba5_A00= cte * FKL_rd2.Ba5_A00.*1e6./ 660;
FKL_rd2.Ba6_A00= cte * FKL_rd2.Ba6_A00.*1e6./ 880;
FKL_rd2.Ba7_A00= cte * FKL_rd2.Ba7_A00.*1e6./ 950;

FKL_rd2.Ba1_A08= cte * FKL_rd2.Ba1_A08.*1e6./ 370;
FKL_rd2.Ba2_A08= cte * FKL_rd2.Ba2_A08.*1e6./ 470;
FKL_rd2.Ba3_A08= cte * FKL_rd2.Ba3_A08.*1e6./ 520;
FKL_rd2.Ba4_A08= cte * FKL_rd2.Ba4_A08.*1e6./ 590;
FKL_rd2.Ba5_A08= cte * FKL_rd2.Ba5_A08.*1e6./ 660;
FKL_rd2.Ba6_A08= cte * FKL_rd2.Ba6_A08.*1e6./ 880;
FKL_rd2.Ba7_A08= cte * FKL_rd2.Ba7_A08.*1e6./ 950;
%%
neph_eco=[neph_19;neph_20];
neph_eco.Properties.DimensionNames{1}='Time';
neph_eco.Var1=[];
neph_eco.Var2=[];
neph_eco.TemperatureInlet=[];
neph_eco.TemperatureOutlet=[];
neph_eco.Pressure=[];
neph_20_M903.Var1=[];
neph_20_M903.Var2=[];
neph_20_M903.Temperature=[];
neph_20_M903.Properties.DimensionNames{1}='Time';
 neph_20_M=retime(neph_20_M903,'daily',@nanmedian);
 neph_eco_d=retime(neph_eco,'daily',@nanmedian);
%%
figure
%%

% % % FKL_rd=read_betsy('FKL_2004_2017',FKL_st.name);
%%
% % % FKL_AE3_2011 =FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2011 AE31',4,8763);
% % % FKL_AE3_2012 = FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2012 AE31',4,8763);
% % % FKL_AE3_2013 = FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2013 AE31',4,8763);
% % % FKL_AE3_2014 = FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2014 AE31',4,8763);
% % % FKL_AE3_2014V2 = FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2014 AE33',4,8763);
% % % 
% % % FKL_AE3_2015 = FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2015 AE33',4,8763);
% % % FKL_AE3_2016 = FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2016 AE33',4,8763);
% % % FKL_AE3_2017 = FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2017 AE33',4,8763);
% % % FKL_AE3_2018 = FKL_AE3_importfile('Finokalia aethalometer dry 2011 - 2018.xlsx','2018 AE33',4,8763);
% % % 
% % % FKL_AE22_2003=FKL_AE22_importfile('Finokalia aethalometer wet and dry 2003 - 2010.xlsx','2003 AE22',5,8763);
% % % FKL_AE22_2004=FKL_AE22_importfile('Finokalia aethalometer wet and dry 2003 - 2010.xlsx','2004 AE22',5,8763);
% % % FKL_AE22_2005=FKL_AE22_importfile('Finokalia aethalometer wet and dry 2003 - 2010.xlsx','2005 AE22',5,8763);
% % % FKL_AE22_2006=FKL_AE22_importfile('Finokalia aethalometer wet and dry 2003 - 2010.xlsx','2006 AE22',5,8763);
% % % FKL_AE22_2007=FKL_AE22_importfile('Finokalia aethalometer wet and dry 2003 - 2010.xlsx','2007 AE22',5,8763);
% % % FKL_AE22_2008=FKL_AE22_importfile('Finokalia aethalometer wet and dry 2003 - 2010.xlsx','2008 AE22',5,8763);
% % % FKL_AE22_2009=FKL_AE22_importfile('Finokalia aethalometer wet and dry 2003 - 2010.xlsx','2009 AE22',5,8763);
% % % FKL_AE22_2010=FKL_AE22_importfile('Finokalia aethalometer wet and dry 2003 - 2010.xlsx','2010 AE22',5,8763);
% % % 
% % % FKL_PSAP=FKL_PSAP_importfile('Finokalia psap wet 2000 - 2002.xlsx','Φύλλο1',5,18075);
% % % 
% % % FKL_AE31=[FKL_AE3_2011;FKL_AE3_2012;FKL_AE3_2013;FKL_AE3_2014];
% % % FKL_AE33=[FKL_AE3_2014V2;FKL_AE3_2015;FKL_AE3_2016;FKL_AE3_2017;FKL_AE3_2018];
% % % FKL_AE22=[FKL_AE22_2003;FKL_AE22_2004;FKL_AE22_2005;FKL_AE22_2006;FKL_AE22_2007;FKL_AE22_2008;FKL_AE22_2009;FKL_AE22_2010];
% % % 
% % % AE22=timetable(FKL_AE22.VarName4, FKL_AE22.black_carbonngm3Wavelength370nm,FKL_AE22.black_carbonngm3Wavelength880nm);
% % % AE31=timetable(FKL_AE31.VarName4, FKL_AE31.equivalent_black_carbonngm1, FKL_AE31.equivalent_black_carbonngm2,FKL_AE31.equivalent_black_carbonngm3,FKL_AE31.equivalent_black_carbonngm4,FKL_AE31.equivalent_black_carbonngm5,FKL_AE31.equivalent_black_carbonngm6,FKL_AE31.equivalent_black_carbonngm7);
% % % ind=isnat(AE31.Time);
% % % AE31(ind,:)=[];
% % % AE33=timetable(FKL_AE33.VarName4, FKL_AE33.bypass_fraction3, FKL_AE33.sensing_zero_signal4,FKL_AE33.reference_beam_signal4,FKL_AE33.attenuation_coefficient4,FKL_AE33.reference_zero_signal5,FKL_AE33.equivalent_black_carbonngm6,FKL_AE33.sensing_beam_signal6);
% % % ind=isnat(AE33.Time);
% % % AE33(ind,:)=[];
% % % AE3=[AE31;AE33];
% % % PSAP=timetable(FKL_PSAP.TimeUTC2,FKL_PSAP.Mm);
% % % 
% % % % nommer les variables
% % % AE3.Properties.VariableNames={'Ba1_A'; 'Ba2_A'; 'Ba3_A' ;'Ba4_A' ;'Ba5_A' ;'Ba6_A' ;'Ba7_A'};
% % % AE22.Properties.VariableNames={'Ba1_A'; 'Ba6_A'};
% % % PSAP.Properties.VariableNames={'BaG_A'};
% % % %enlever missing > 1e4
% % % names3=fieldnames(AE3);
% % % c= startsWith(names3,'Ba');
% % % names3=names3(c);
% % % for i=1:length(names3)
% % %     ind=AE3.(names3{i})>9999;
% % %     AE3.(names3{i})(ind)=NaN;
% % % end
% % % 
% % % names22=fieldnames(AE22);
% % % c= startsWith(names22,'Ba');
% % % names22=names22(c);
% % % for i=1:length(names22)
% % %     ind=AE22.(names22{i})>9999;
% % %     AE22.(names22{i})(ind)=NaN;
% % % end
% % % 
% % % %PSAP abs sinon BC
% % % cte = 14625e-9/3.5;
% % % 
% % % AE3.Ba1_A= cte * AE3.Ba1_A.*1e6./ 370;
% % % AE3.Ba2_A= cte * AE3.Ba2_A.*1e6./ 470;
% % % AE3.Ba3_A= cte * AE3.Ba3_A.*1e6./ 520;
% % % AE3.Ba4_A= cte * AE3.Ba4_A.*1e6./ 590;
% % % AE3.Ba5_A= cte * AE3.Ba5_A.*1e6./ 660;
% % % AE3.Ba6_A= cte * AE3.Ba6_A.*1e6./ 880;
% % % AE3.Ba7_A= cte * AE3.Ba7_A.*1e6./ 950;
% % % 
% % % AE22.Ba1_A= cte * AE22.Ba1_A.*1e6./ 370;
% % % AE22.Ba6_A= cte * AE22.Ba6_A.*1e6./ 880;
% % % %median diurnale
% % % %delete à la mains les NaT dans AE22 car j'ai fait lire plus que le nb de
% % % %lignes nécessaires
% % % ind=isnat(AE22.Time);
% % % AE22(ind,:)=[];
% % % 
% % % % % ind=isnat(AE3.Time);
% % % % % AE3(ind,:)=[];
% % % TT=synchronize(AE3,AE22);
% % % FKL_rd_mco=synchronize(TT,PSAP);
% % % 
% % % % % %IN Mm-1
% % % % % names=fieldnames(FKL_rd_mco);
% % % % % c=startsWith(names,'Ba')==1 & contains(names,'G')==0;
% % % % % namesc=names(c);
% % % % % for i=1:length(namesc)
% % % % %     FKL_rd_mco.(namesc{i})=FKL_rd_mco.(namesc{i}).*1e6;
% % % % % end
% % % 
% % % % daily mean, apply a data coverage threshold ?
% % %  station_datD=retime(FKL_rd_mco,'daily',@nanmedian);
% % %  station_datN=retime(FKL_rd_mco,'daily','count');
% % %  namesOK=fieldnames(FKL_rd_mco);
% % % CB=startsWith(namesOK,'B');
% % %  SB=namesOK(CB);
% % %  for i=1:length(SB)
% % %  indN=station_datN.(SB{i})<6 ;
% % %  a(i)=sum(indN)*100/sum(~isnan(station_datD.(SB{i})));
% % %  station_datD.(SB{i})(indN)=NaN;
% % %  end
% % %  FKL_rd_mco2=station_datD;
%%
% read data prepared by Betsy
FKL_rd_B=read_betsy_2026('C:\github_trend\raw_data\FKL\fkl_neph_aethBC_aethbap_fixed',FKL_st.name); %[output:05eacf68]
%%
FKL_rd_B.Properties.VariableNames{45}='Ba1_A13';
FKL_rd_B.Properties.VariableNames{45}='Ba2_A13';
FKL_rd_B.Properties.VariableNames{45}='Ba3_A13';
FKL_rd_B.Properties.VariableNames{45}='Ba4_A13';
FKL_rd_B.Properties.VariableNames{45}='Ba5_A13';
FKL_rd_B.Properties.VariableNames{45}='Ba6_A13';
FKL_rd_B.Properties.VariableNames{45}='Ba7_A13';
%%
names_B=fieldnames(FKL_rd_B);
names_X=names_B(startsWith(names_B,'X'));
for i=1:length(names_X)
    FKL_rd_B.(names_X{i})(FKL_rd_B.(names_X{i})>9999)=NaN;
end
%%
% compute abs from A13 (with AE33 correction) just see about H
sigma=[18.47 14.54 13.14 11.58 10.35 7.77 7.19]; %for the 7 wavelengths with an exponent of -1
FKL_rd_B.Ba1_A13x=FKL_rd_B.X370_A13.*sigma(1).*1e-3./1.57;
FKL_rd_B.Ba2_A13x=FKL_rd_B.X470_A13.*sigma(2).*1e-3./1.57;
FKL_rd_B.Ba3_A13x=FKL_rd_B.X520_A13.*sigma(3).*1e-3./1.57;
FKL_rd_B.Ba4_A13x=FKL_rd_B.X590_A13.*sigma(4).*1e-3./1.57;
FKL_rd_B.Ba5_A13x=FKL_rd_B.X660_A13.*sigma(5).*1e-3./1.57;
FKL_rd_B.Ba6_A13x=FKL_rd_B.X880_A13.*sigma(6).*1e-3./1.57;
FKL_rd_B.Ba7_A13x=FKL_rd_B.X950_A13.*sigma(7).*1e-3./1.57;
%%
%AE31 correspond to _A12

cte = 14625e-9/3.5;

FKL_rd_B.Ba1_A12= cte * FKL_rd_B.X370_A12 .*1e6./ 370;
FKL_rd_B.Ba2_A12= cte * FKL_rd_B.X470_A12.*1e6./ 470;
FKL_rd_B.Ba3_A12= cte * FKL_rd_B.X520_A12.*1e6./ 520;
FKL_rd_B.Ba4_A12= cte * FKL_rd_B.X590_A12.*1e6./ 590;
FKL_rd_B.Ba5_A12= cte * FKL_rd_B.X660_A12.*1e6./ 660;
FKL_rd_B.Ba6_A12= cte * FKL_rd_B.X880_A12.*1e6./ 880;
FKL_rd_B.Ba7_A12= cte * FKL_rd_B.X950_A12.*1e6./ 950;
%%

%AE22 correspond to _A11
FKL_rd_B.Ba1_A11= cte * FKL_rd_B.X370_A11.*16.6.*1e6./(3.5* 370);
FKL_rd_B.Ba6_A11= cte * FKL_rd_B.X880_A11.*1e6./(3.5* 880);

%%
% TO DO abs
%1) compute cte between _A13 from Betsy and A_AE3 from _mco2 and adjust
%mco2 read from 11.2014 to 2018 and rd2 _A00
% use 880 nm
% FIK_rd is not usefull for abs

% neph: try to use 2012-2025 or eventually from 15.9.2008 with invalidating
% 9.2011 until 2012 since it is once again PM1
%%
datevec(FKL_rd2.Time(1))
datevec(FKL_rd2.Time(end))
%prctile(FKL_rd.U_S,[5 50 95])
plotFigControl(FKL_rd2,FKL_st.name);
plotFigControl(FKL_rd_B,FKL_st.name);

%%
datevec(FKL_rd.Time(1))
datevec(FKL_rd.Time(end))
prctile(FKL_rd.U_S,[5 50 95])
plotFigControl(FKL_rd,FKL_st.name);
% several size cut, NTO OK
% sc: so many instrumental and size cut changes that no homogeneous time
% series is possible
% abs: 2 size cut, max in summer, no negatives

%%
lambdaSC=[450;550;700]*ones(1,4);
lambdaAE=[467;530;660]*ones(1,4);
lambdaAE7=[370 470 520 590 660 880 950];
FKL_cal=compute_exp_SSA(FKL_rd,lambdaSC, lambdaAE);
plotFigControl_cal(FKL_cal,FKL_st.name);
%%
%Questions:
%problems:

%%
FKL_tr=FKL_rd_mco2;
FKL_tr.y=year(FKL_tr.Time);
% begin at the beginning of a year: 2004
%end 2018

P=timerange('2004-01-01','2019-01-01');
FKL_tr=FKL_tr(P,:);
%put AE22 and AE3 together
P1=timerange('2004-01-01','2011-01-01');
FKL_tr.Ba1_A=FKL_tr.Ba1_A_AE3;
FKL_tr.Ba1_A(P1)=FKL_tr.Ba1_A_AE22(P1);
FKL_tr.Ba1_A_AE3=[];
FKL_tr.Ba1_A_AE22=[];

FKL_tr.Ba6_A=FKL_tr.Ba6_A_AE3;
FKL_tr.Ba6_A(P1)=FKL_tr.Ba6_A_AE22(P1);
FKL_tr.Ba6_A_AE3=[];
FKL_tr.Ba6_A_AE22=[];
% 880 is better for difference between AE31 and AE33
%use it for the trend analysis
FKL_tr.Ba1_A=[];
FKL_tr.Ba2_A=[];
FKL_tr.Ba3_A=[];
FKL_tr.Ba4_A=[];
FKL_tr.Ba5_A=[];
FKL_tr.Ba7_A=[];
%PSAP: not possible to homogeneise
FKL_tr.BaG_A=[];
FKL_tr=FKL_tr; FKL_st=FKL_st;
%%
FKL_result= all_trend_STN(FKL_tr,FKL_st); 
%%
FKL_result_D=FKL_result;

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":33.9}
%---
%[output:6025f3a3]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7t3X\/sVfV9x\/E3jglMtANGRRCjFWRViyWQQQiKcYuZa\/oHSSvMGZBCEFp+VCjhR6lE0VGD0qIQDeWHkoaoNLiU6pJmIZW4pk6IofEHgutcixMEYVQ6RtWyfI69331\/3O\/9\/Drve7nv+7wJaSqfz\/ve+zjnfe6Lzzn33B7nzp07JzwQQAABBBBAAIEmEuhBgGmircVLRQABBBBAAIFCgADDjoAAAggggAACTSdAgGm6TcYLRgABBBBAAAECDPsAAggkCTzzzDOybNky79wvf\/nL8t3vflfOnDkjM2bMkP3798sNN9wgmzdvlv79+3eY\/9JLL8k3vvEN+fDDD4v\/Pm3aNFm0aJH8yZ\/8iSxdulR27drV7VzvC2EAAgiYEiDAmNqcvBkE6idQdoBx4eWee+6RDz74oEN46du3bxF+CDD127Y8EwLNIECAaYatxGtE4DwUKDPAvPnmmzJ37lz5j\/\/4jy7hxf0HAsx5uAPwkhBosAABpsEbgKdHoFkF2geY1atXy+TJk2u+lRMnTlQ9hXTw4EGZM2dOt+GFANOsewivGwFdAQKMri\/VETArUEaAOX78eLHy8vbbb1ddeangsQJjdjfijSGQLECASaZjIgKtLZAbYFatWiUPPPCA\/Nu\/\/VsBecstt8jatWvlkksu6QJLgGntfY13j0A1AQIM+wUCCCQJhFwDU\/kGUp8+faT9KaSLL75Y3J\/\/+q\/\/antu9\/+ffPJJGTVqFAEmaYswCYHWEiDAtNb25t0iUJpAToBp\/yIGDBjQ9s2jv\/mbv5FHHnmkCDftH6zAlLbZKISAGQECjJlNyRtBoL4CZQSY22+\/XebNm1fcT8Z9jdo93AXB7r\/36NGj7Q0RYOq7bXk2BJpBgADTDFuJ14jAeSiQcw2MezsupKxYsULcfV52794ts2fPlo8\/\/liuuuoq2bhxo1x99dUEmPNwu\/OSEDhfBAgw58uW4HUg0GQCOQHmC1\/4gmzatEkGDhxYvOv2Kyzu\/\/\/93\/+93HvvvdKrV68uf9\/dXXybjI+XiwACmQIEmExApiPQqgI5AaZaCHn11VflrrvuKn5GoGfPnvLEE08U30zqHHAIMK26x\/G+EegoQIBhj0AAgSSBsgOMO320bt062bBhQ\/F6Ro8eLY899pgMGjSoywpNrRd89913y5IlS5LeE5MQQKB5BAgwzbOteKUInFcCZQcY9+beeecdmTVrVtuN7RYsWFD8uONHH33U9ltIPgQCjE+Iv0fAhgABxsZ25F0gUHcBjQBz7tw52bZtm9x3333F+6ncG+Yv\/\/IvCTB138I8IQLntwAB5vzePrw6BBBAAAEEEKgi0JIB5vTp07Jy5Ur5q7\/6K+8P0LHXIIAAAggggMD5J9ByAebkyZPF768899xzxQ2zfL+ge\/5tMl4RAggggAACCLRMgDl79qy8+OKLRWj5z\/\/8z2LLE2BoAAQQQAABBJpToGUCTOWCQ\/e7K3feeac8++yzMn\/+fFZgmnO\/5VUjgAACCLS4QMsEmJ07dxa\/fHvHHXcUPxw3ffr04jdYOIXU4h3A20cAAQQQaEqBlgkw7bfOoUOHCDBNubvyohFAAAEEEPhUgAAzeTL7AgIIIIAAAgg0mQABJiDAHD58WNyfysPd2tz94YEAAggggMD5LuB+W8z9sfYgwHgCjAsuixcvlpdffrlt20+bNk2mTp1qbV+o6\/txv3vzP\/\/zP\/Jnf\/ZnJhurrpidngxbPX1ssdUT0Kns9tk\/\/dM\/lcsvv9zcsZYA4wkwv\/jFL4oLf9esWSNDhgwp9jBWYPIb7cyZM3L06NHCsnfv3vkFqdAmgK3ezoAttnoCOpXdPnv8+HEZOnSouWMtASYwwGzfvl3GjRuns4e1YNX\/\/d\/\/lffee08uu+wyc03V6M2Jrd4WwBZbPQGdypb3WQIMAUanazxVLTdVQ0DbPSm2elsAW2z1BHQqW95nCTAEGJ2uIcA0xNU9qeUDVsNQ\/\/jE2OptAWx1bC27EmAIMDpdQ4BpiCsBRpfd8oeBrpy\/OrZ+o5QRll1bMsDE7ASVi3i5BiZGzT\/WclP5373uCGz1fLHFVk9Ap7LlfZYA49lnCDA0lY6AXlXLByw9tbDK2IY5pYzCNkXNP8eyKwGGAOPvAIURlptKgSuqJLZRXFGDsY3iihqMbRRX8GDLrgQYAkxwI5Q50HJTlemUUgvbFLWwOdiGOaWMwjZFzT\/HsisBhgDj7wCFEZabSoErqiS2UVxRg7GtztX551aiUP842N0x1t1w7c\/\/\/M+5N1QAoLuzrvvje1jeZwkwBBjf\/q\/y95abSgUsoii2EViRQ7HtClbt51YiWRmeIDB27NjiDvG+EGN5nyXAEGASWid\/iuWmytfJq4Btnl+t2dh21an2cyt6W4DKTsD9Nt+6desk5NuxlvdZAgwBpiFHBMtN1RDQdk+Krd4WwLb7ABPyYaq3ZVqrcsy3Yy3vswQYAkxDOt9yUzUElABTF3b2WwJMXXa0Ej+XLO+zBJgSd5TzYcdultdguakavQ2w1dsC2BJg9Pau8MqswHxqRYAhwIR3TYkj+SAoEbNTKWyx1RMgwNTTtrvnIsAQYIL2w5gdJagggwoBPmT1dgRssdUTIMDU05YAU1ubFRhWYBrSj3zI6rFji62eAAGmnrYEGAJM1v7GCkwWX7eT+ZDVcWV1S88V2+q2HCN197lq1WPMLR9rWYFhBab+3ccpJFVzywcsVbiA4tiyAhOwm6gPIcB8SkyAIcCoN1u1J+CDQI8dW2z1BAgw9bTlFBKnkLL2t5ikm\/VELTaZD1m9DY4ttnoCBJh62hJgCDBZ+xsBJouv28l8yOq4uqrYYqsnQICppy0BhgCTtb8RYLL4CDA6fDWrEmD00LElwOjtXeGVYz6XLO+zXAPj2WdidpTw3Y+Rlpuq0VsXW70tgC0BRm\/vCq8c87lkeZ8lwBBgwrumxJGWm6pEpqRS2CaxBU3ClgATtKMoDyLAfApMgCHAKLda9fJ8EOixY4utngABpp623T0XAYYAE7QfxuwoQQUZVAjwIau3I2CLrZ4AAaaetgSY2tqswLAC05B+5ENWjx1bbPUECDD1tCXAEGCy9jdWYLL4up3Mh6yOK6tbeq7YVrflGKm7z1WrHmNu+VjLCgwrMPXvPk4hqZpbPmCpwgUUx5YVmIDdRH0IAeZTYgIMAUa92ao9AR8EeuzYYqsn0H2AWbBggYwdO7aeT92yz\/Xuu+\/K4sWLZfv27TJu3LiaDpaPBwQYAkxDDgKWm6ohoO2eFFu9LYBtV9vDhw8XH6Yvv\/yyHjyVuwi4sLhmzRq5\/PLLCTDsH9UFYpbqMAwX4IMg3Cp2JLaxYuHjsa1u5UKM+5PzcLb\/\/d\/\/LX\/xF38hPXv2zCnVEnNdcPGFFwdheZ9lBYYVmIY0u+WmaggoKzB1YWe\/1WPGVsfWsisBhgCj0zWeqpabqiGgBJi6sLPf6jFjq2Nr2ZUAQ4DR6RoCTENcrS8ZNwz1j09s+cMA20YL6Dy\/5X2WAEOA0ekaAkxDXAkwuuyWPwx05fzVsfUbpYyw7EqAIcCk9ET2HMtNlY2TWQDbTMAa07HFVk9Ap7LlfZYAQ4DR6RpWYBriygqMLrvlDwNdOX91bP1GKSMsuxJgCDApPZE9x3JTZeNkFsA2E5AVGD1AbOtua\/l4QIAhwNS9oVgl0CW3fMDSlfNXx9ZvlDoC21S52vMsuxJgCDA6XcMppIa4Eg512S1\/GOjK+atj6zdKGWHZlQBDgEnpiew5lpsqGyezALaZgJzm0APEtu62lo8HBBgCTN0bilUCXXLLByxdOX91bP1GqSOwTZXjFJKOXMlVT548KVu2bJFt27bJmTNnil\/hXLRokYwcOVJ69OhR89n+8Ic\/yIsvvihr166V119\/XQYMGCBTpkyRr33ta9KvX79u5\/JbSCVvxD+W42Cl40o41HPFFltdAZ3qlo+1TbMCc\/r0aVm5cqXs2bNH5syZI4MHD5YdO3bI3r17ZcOGDTJhwoSaW3\/nzp2ydOlSufnmm4vg8sYbb8hTTz0ln\/\/852XdunXdhhgCDE2lI6BX1fIBS08trDK2YU4po7BNUfPPsezaNAFm165dsnz5clm\/fr1MnDix2Gou1CxbtkxOnTpVM4QcP35cZs2aJUOGDJHVq1dL3759i\/m7d++W2bNny\/e+9z350pe+VHVPIMD4GyRlhOWmSvEocw62ZWp2rIUttnoCOpUt77NNEWDc6SK3euL+150CqgSQSghZsGBBcVpp1KhRVfeAQ4cOyfTp02Xu3LnF6kvlcezYMZk5c6aMHz9elixZQoDR6Z+qVS03VR0Zsa0zNvutHji2OraWXZsiwNQKGgcPHizCyfz582Xy5MlV94Df\/OY3xZhbb71VFi9e3Ha9TCXYuFDjwk21ByswNJWOgF5VywcsPbWwytiGOaWMwjZFzT\/HsmtTBJjuVlDcpqv83aRJk4oLeqs9Pv744+I6GXcB8EMPPVQEmQ8++EDuv\/\/+4lqYjRs3ytVXX02A8fdCaSMsN1VpSImFsE2EC5iGbQBS4hBsE+E80yy7NlWAmTdvXpdVlhMnTsiMGTOKbyR1dxrIbd+zZ8\/KD37wg+IUVOUxYsQI+f73vy\/uf7t7sAJDU+kI6FW1fMDSUwurjG2YU8oobFPU\/HMsu7ZEgHEX+7oLdZ999ln5yle+IjfeeKO899578sMf\/lA++eST4sLga665puYKjLvOZuzYscWYQYMGFX94pAu4pjpy5Ihceuml0qdPn\/RCzOwigK3eToEttnoCOpXdPusuwxg6dKj07t1b50kaVLXpA0zIKaQXXnihOL3kThm5AFO5Z8z7779fXDvTq1evIsRcfPHFXTZDZQWm\/V9MmzZNpk6d2qBNZuNp3YqYa6qBAwcW\/jzKE8C2PMvOlbDFVk9Ap7LbZ90fd6aBAKNjXLNq7kW87rqXn\/\/857Jp06biA7P94+mnny7Cy9atW2X48OHdBpg1a9YUX8NmBaacHcB9o+zo0aPFSpa1pipHKL0Ktul2vpnY+oTS\/x7bdLtaM52ru5UIKzA6vt6quV+jdgHG3YV38+bNctlllyUFmO3btxfX2fAoR8DyedlyhNKrYJtu55uJrU8o\/e+xTberNdOya1OcQnIbx93IbsWKFR3uuht6Izs3151CevDBB6ueQrrggguKutV+UoCLeGkqHQG9qpYPWHpqYZWxDXNKGYVtipp\/jmXXpgkw7neQ3D1c9u\/fL+4alGHDhlX9KYHKas2+ffvaTgtVfobgJz\/5ibivW7uvUbvrX5555hn51a9+VfOnCAgw\/gZJGWG5qVI8ypyDbZmaHWthi62egE5ly\/ts0wQYt2k7\/5jj6NGjZeHChTJmzJi2C3OrBRg31\/335557rvjm0YEDB4oLdt1PErgb2HX3DSQ3jwBDU+kI6FW1fMDSUwurjG2YU8oobFPU\/HMsuzZVgPFvqvJHEGDKN3UVLTeVjlh4VWzDrWJHYhsrFj4e23CrmJGWXQkwnj2BABPTKuFjLTdVuILOSGx1XAneeq7Y6tlaPh4QYAgwep1To7LlpmoIaLsnxVZvC2CLrZ6ATmXL+ywBhgCj0zWeqpabqiGgBJi6sLPf6jFjq2Nr2ZUAQ4DR6RoCTENcWYrXZbf8YaAr56+Ord8oZYRlVwIMASalJ7LnWG6qbJzMAthmAnLqUw8Q27rbWj4eEGAIMHVvKFYJdMktH7B05fzVsfUbpY7ANlWu9jzLrgQYAoxO13AKqSGuhENddssfBrpy\/urY+o1SRlh2JcAQYFJ6InuO5abKxsksgG0mIKc59ACxrbut5eMBAYYAU\/eGYpVAl9zyAUtXzl8dW79R6ghsU+U4haQjZ6AqN7LT2YgcrHRcCYd6rthiqyugU93ysZYVGFZgdLrGU9VyUzUEtN2TYqu3BbDFVk9Ap7LlfZYAQ4DR6RoCTENcWSXQZbf8YaAr56+Ord8oZYRlVwIMASalJ7LnWG6qbJzMAthmAtaYji22egI6lS3vswQYAoxO17AC0xBXVmB02S1\/GOjK+atj6zdKGWHZlQBDgEnpiew5lpsqGyezALaZgKzA6AFiW3dby8cDAgwBpu4NxSqBLrnlA5aunL86tn6j1BHYpsrVnmfZlQBDgNHpGk4hNcSVcKjLbvnDQFfOXx1bv1HKCMuuBBgCTEpPZM+x3FTZOJkFsM0E5DSHHiC2dbe1fDwgwBBg6t5QrBLokls+YOnK+atj6zdKHYFtqhynkHTkDFTlTrw6G5GDlY4r4VDPFVtsdQV0qls+1rICwwqMTtd4qlpuqoaAtntSbPW2ALbY6gnoVLa8zxJgCDA6XUOAaYgrqwS67JY\/DHTl\/NWx9RuljLDsSoAhwKT0RPYcy02VjZNZANtMwBrTscVWT0CnsuV9lgBDgNHpGlZgGuLKCowuu+UPA105f3Vs\/UYpIyy7EmAIMCk9kT3HclNl42QWwDYTkBUYPUBs625r+XhAgCHA1L2hWCXQJbd8wNKV81fH1m+UOgLbVLna8yy7EmAIMDpdwymkhrgSDnXZLX8Y6Mr5q2PrN0oZYdmVAEOASemJ7DmWmyobJ7MAtpmAnObQA8S27raWjwcEGAJM3RuKVQJdcssHLF05f3Vs\/UapI7BNleMUko6cgarciVdnI3Kw0nElHOq5YoutroBOdcvHWlZgWIHR6RpPVctN1RDQdk+Krd4WwBZbPQGdypb3WQIMAUanawgwDXFllUCX3fKHga6cvzq2fqOUEZZdCTAEmJSeyJ5juamycTILYJsJWGM6ttjqCehUtrzPEmAIMDpdwwpMQ1xZgdFlt\/xhoCvnr46t3yhlhGVXAgwBJqUnsudYbqpsnMwC2GYCsgKjB4ht3W0tHw8IMASYujcUqwS65JYPWLpy\/urY+o1SR2CbKld7nmVXAgwBRqdrOIXUEFfCoS675Q8DXTl\/dWz9RikjLLsSYAgwKT2RPcdyU2XjZBbANhOQ0xx6gNjW3dby8YAAQ4Cpe0OxSqBLbvmApSvnr46t3yh1BLapcpxC0pEzUJU78epsRA5WOq6EQz1XbLHVFdCpbvlY21QrMCdPnpQtW7bItm3b5MyZMzJu3DhZtGiRjBw5Unr06OHd+r\/+9a9lzZo18tOf\/rQYe+utt8rixYvliiuu6HYuAcbLmjTAclMlgZQ4CdsSMTuVwhZbPQGdypb32aYJMKdPn5aVK1fKnj17ZM6cOTJ48GDZsWOH7N27VzZs2CATJkyoufVdEJk3b55cfvnlMm3aNDl37pw8+uijcsEFF8jGjRvl6quvrjqfAENT6QjoVbV8wNJTC6uMbZhTyihsU9T8cyy7Nk2A2bVrlyxfvlzWr18vEydOLLaaCzXLli2TU6dOybp166Rfv35Vt+aHH35YrNT89re\/LULLZz\/72WLca6+9JjNmzJA777xT5s6dW3UVhwDjb5CUEZabKsWjzDnYlqnZsRa22OoJ6FS2vM82RYBxp4uWLl1anDZau3at9O3bt21L7969WxYsWFCcVho1alTVPeDVV1+Vu+66qwgvlfDjBlbqulWYBx54QC666KIu8wkwNJWOgF5VywcsPbWwytiGOaWMwjZFzT\/HsmtTBJhjx47JzJkzZfz48bJkyZIOW+zgwYMyffp0mT9\/vkyePLnq1nz66aeLgLN582a57LLL\/Fu83QgCTBRX8GDLTRWMoDQQWyVYEcEWWz0BncqW99mmCDCHDh0qQoo7zTNlypQOW7nyd5MmTSpOE1V7PPTQQ\/L666\/LfffdJzt37hQXaNxpp5tvvlm++c1vyrXXXtvtnkOAoal0BPSqWj5g6amFVcY2zCllFLYpav45ll2bKsC4i3A7r7KcOHGiuI7FfSOp8+qM27Rnz54trp05cOBAceGuW4H56le\/Kh988IFs2rSpuI7GXcTb3emnSoBxp6nGjh1b7C2DBg0q\/vBIF3BNdeTIEbn00kulT58+6YWY2UUAW72dAlts9QR0Krt91p3FGDp0qPTu3VvnSRpU1XyAqVzn4i4Cnj17tixcuFB69uxZcLsPUBeKBgwYIA8\/\/HCHa2sq26MSYNpvH\/ctpqlTpzZok9l4WhcsXVMNHDhQevXqZeNNnSfvAlu9DYEttnoCOpXdPuv+jBgxggCjQ1y7auU0UbUVGN8ppI8\/\/lhWrVolL7zwgmzdulWuv\/76Dk\/mvtXkvo795JNPylVXXdXlhVQCjLt\/zJAhQ1iBKWkHcMHy6NGjxUqWtX8VlESUXAbbZDrvRGy9RMkDsE2mqznRuR4\/fpwVGB1ef9Xci3gfeeQRee6554oAM3z48A5P+Mwzz8hjjz1W9e\/cQK6B8W+flBGWz8umeJQ5B9syNTvWwhZbPQGdypb32aY4hZT7Nernn39e7rnnnuKbSO5amcrDXRPjVlZ+9rOfdfsNJQIMTaUjoFfV8gFLTy2sMrZhTimjsE1R88+x7NoUAcZtIncNy4oVKzrcdTf0RnbuWpe7775bhg0bVnwTqXIfmcOHDxf\/ffTo0UXtCy+8sMveQIDxN0jKCMtNleJR5hxsy9RkBUZPE9t62Fo+HjRNgHG\/g+R+t2j\/\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\/ftpSzZ8\/K\/fffLz\/+8Y\/ld7\/7HQHGK6YzwHJT6YiFV8U23Cp2JLaxYuHjsQ23ihlp2bUpAow7XbR06dLitNHatWulb9++bdvPnQpasGBBcSpo1KhR3u26c+fOYsVm6tSpct999xFgvGI6Ayw3lY5YeFVsw61iR2IbKxY+Httwq5iRll2bIsAcO3ZMZs6cKePHj5clS5Z02HYHDx6U6dOny\/z582Xy5Mk1t6sb604\/uetm+vfvL3fccQcBJqYTShxrualKZEoqhW0SW9AkbIOYkgZhm8TmnWTZtSkCzKFDh4o8r0ZnAAAdOUlEQVSQMnfuXJkyZUqHDVb5u0mTJhXBpLtH5Rqa3r17y7333iuvvvoqAca76+sNsNxUemphlbENc0oZhW2KWtgcbMOcYkdZdm2qADNv3rwuqywnTpyQGTNmFBfidl6dqWzoc+fOybPPPiubN2+WjRs3ypVXXimh3y4KHRe7U7X6eMtN1ehti63eFsAWWz0BncqW99mWCDCvvfaa3H333bJq1Sq55ZZbir0kNJhUxrnrbMaOHVvMHTRoUPGHR7qAa6ojR47IpZdeKn369EkvxMwuAtjq7RTYYqsnoFPZ7bPuMoyhQ4eKOwNh6dH0AcZ3CsndO8aFj2HDhhXfYurZs2dSgGm\/0adNm1ZcBMwjXcB9G8w11cCBA6VXr17phZjZRQBbvZ0CW2z1BHQqu33W\/RkxYgQBRoe4dtWci3grKyi1nsHdU2br1q0yfPjwLsMq89esWSNDhgxhBaakHcB9o+zo0aPFSpa1fxWURJRcBttkOu9EbL1EyQOwTaarOdG5Hj9+nBUYHV5\/1ZyvUbsNd+DAgS5P8tZbb8mDDz4o3\/72t+W6666Ta6+9Vi655JJuA4zvhnf+d8GI9gKWz8s2ektjq7cFsMVWT0CnsuV9tilOIbnN6m5kt2LFig533U25kV1lF4m9BoYAU25zWW6qcqXiq2EbbxY6A9tQqfhx2Mabhcyw7No0AcZdy7J48WLZv3+\/uGtQ3DUt1X5KoLJas2\/fvm5PC7mNToAJ2fX1xlhuKj21sMrYhjmljMI2RS1sDrZhTrGjLLs2TYBxG63zjzmOHj1aFi5cKGPGjGn7MUcCTOzu3ZjxlpuqMaL\/\/6zY6m0BbLHVE9CpbHmfbaoAo7N5a1cNXalpxGtr5ue03FSN3i7Y6m0BbLHVE9CpbHmfJcB49hkCDE2lI6BX1fIBS08trDK2YU4po7BNUfPPsexKgCHA+DtAYYTlplLgiiqJbRRX1GBso7iiBmMbxRU82LIrAYYAE9wIZQ603FRlOqXUwjZFLWwOtmFOKaOwTVHzz7HsSoAhwPg7QGGE5aZS4IoqiW0UV9RgbKO4ogZjG8UVPNiyKwGGABPcCGUOtNxUZTql1MI2RS1sDrZhTimjsE1R88+x7EqAIcD4O0BhhOWmUuCKKoltFFfUYGyjuKIGYxvFFTzYsisBhgAT3AhlDrTcVGU6pdTCNkUtbA62YU4po7BNUfPPsexKgCHA+DtAYYTlplLgiiqJbRRX1GBso7iiBmMbxRU82LIrAYYAE9wIZQ603FRlOqXUwjZFLWwOtmFOKaOwTVHzz7HsSoAhwPg7QGGE5aZS4IoqiW0UV9RgbKO4ogZjG8UVPNiyKwGGABPcCGUOtNxUZTql1MI2RS1sDrZhTimjsE1R88+x7EqAIcD4O0BhhOWmUuCKKoltFFfUYGyjuKIGYxvFFTzYsisBhgAT3AhlDrTcVGU6pdTCNkUtbA62YU4po7BNUfPPsexKgCHA+DtAYYTlplLgiiqJbRRX1GBso7iiBmMbxRU82LIrAYYAE9wIZQ603FRlOqXUwjZFLWwOtmFOKaOwTVHzz7HsSoAhwPg7QGGE5aZS4IoqiW0UV9RgbKO4ogZjG8UVPNiyKwGGABPcCGUOtNxUZTql1MI2RS1sDrZhTimjsE1R88+x7EqAIcD4O0BhhOWmUuCKKoltFFfUYGyjuKIGYxvFFTzYsisBhgAT3AhlDrTcVGU6pdTCNkUtbA62YU4po7BNUfPPsexKgCHA+DtAYYTlplLgiiqJbRRX1GBso7iiBmMbxRU82LIrAYYAE9wIZQ603FRlOqXUwjZFLWwOtmFOKaOwTVHzz7HsSoAhwPg7QGGE5aZS4IoqiW0UV9RgbKO4ogZjG8UVPNiyKwGGABPcCGUOtNxUZTql1MI2RS1sDrZhTimjsE1R88+x7EqAIcD4O0BhhOWmUuCKKoltFFfUYGyjuKIGYxvFFTzYsisBhgAT3AhlDrTcVGU6pdTCNkUtbA62YU4po7BNUfPPsexKgCHA+DtAYYTlplLgiiqJbRRX1GBso7iiBmMbxRU82LIrAYYAE9wIZQ603FRlOqXUwjZFLWwOtmFOKaOwTVHzz7HsSoAhwPg7QGGE5aZS4IoqiW0UV9RgbKO4ogZjG8UVPNiyKwGGABPcCGUOtNxUZTql1MI2RS1sDrZhTimjsE1R88+x7EqAIcD4O0BhhOWmUuCKKoltFFfUYGyjuKIGYxvFFTzYsisBhgAT3AhlDrTcVGU6pdTCNkUtbA62YU4po7BNUfPPsexKgCHA+DtAYYTlplLgiiqJbRRX1GBso7iiBmMbxRU82LIrAYYAE9wIZQ603FRlOqXUwjZFLWwOtmFOKaOwTVHzz7HsSoAhwPg7QGGE5aZS4IoqiW0UV9RgbKO4ogZjG8UVPNiyKwGGABPcCGUOtNxUZTql1MI2RS1sDrZhTimjsE1R88+x7EqAIcD4O0BhhOWmUuCKKoltFFfUYGyjuKIGYxvFFTzYsisBhgAT3AhlDrTcVGU6pdTCNkUtbA62YU4po7BNUfPPsexKgCHA+DtAYYTlplLgiiqJbRRX1GBso7iiBmMbxRU82LIrAYYAE9wIZQ603FRlOqXUwjZFLWwOtmFOKaOwTVHzz7HsSoAhwPg7QGGE5aZS4IoqiW0UV9RgbKO4ogZjG8UVPNiya1MFmJMnT8qWLVtk27ZtcubMGRk3bpwsWrRIRo4cKT169Ki5Qc+ePSu7du0q5h84cEAuvvhimThxonzrW9+SK664otu5v\/jFL+SOO+6Q7du3F8\/HoxwBy01VjlB6FWzT7XwzsfUJpf89tul2tWZadm2aAHP69GlZuXKl7NmzR+bMmSODBw+WHTt2yN69e2XDhg0yYcKEbrehCzvf+c53igDjwsiNN94o7733nmzatEk++ugjeeKJJ+T666+vOp8AQ1PpCOhVtXzA0lMLq4xtmFPKKGxT1PxzLLs2TYBx4WP58uWyfv36YuXEPVyoWbZsmZw6dUrWrVsn\/fr1q7o19+3bJ3fddZesWLFCbr\/99rbVmnfffVdmzZol1113naxatUp69erVZT4Bxt8gKSMsN1WKR5lzsC1Ts2MtbLHVE9CpbHmfbYoA41ZQli5dWpw2Wrt2rfTt27dtS+\/evVsWLFhQnFYaNWpU1T3gRz\/6kWzcuLFYqRk+fHjbmHPnzsnq1avllVdekc2bN0v\/\/v0JMDo91KWq5aaqE2G3T4Ot3hbAFls9AZ3KlvfZpggwx44dk5kzZ8r48eNlyZIlHbbywYMHZfr06TJ\/\/nyZPHly1B7w+9\/\/Xu6991556623CDBRcvmDLTdVvk5eBWzz\/GrNxhZbPQGdypb32aYIMIcOHSpCyty5c2XKlCkdtnLl7yZNmlRc0BvzcMFlxowZcttttxXBqGfPnqzAxABmjLXcVBkspUzFthTGqkWwxVZPQKey5X22qQLMvHnzuqyynDhxoggh7htCnVdnau0O77\/\/frFq41Z3Hn\/8cbnmmmuqDq9cA+NOU40dO7YYM2jQoOIPj3QB11RHjhyRSy+9VPr06ZNeiJldBLDV2ymwxVZPQKey22fd59zQoUOld+\/eOk\/SoKotGWDc17HdqSP3jSbfN5gqAab99pk2bZpMnTq1QZvMxtO6r7W7pho4cGDVi6dtvMvGvAts9dyxxVZPQKey22fdnxEjRhBgdIhrV62cJqq2AhN7Csl988h9c2n\/\/v3y6KOPyk033VTzHjKVALNmzRoZMmQIKzAl7QDuguyjR48WK1nW\/lVQElFyGWyT6bwTsfUSJQ\/ANpmu5kTnevz4cVZgdHj9Vcu6iPfNN98svrH0ySefyMMPP9ztt5bavyK+Ru3fPikjLJ+XTfEocw62ZWp2rIUttnoCOpUt77NNcQop92vUbrdw31ZyN8C75JJL5JFHHpHPfe5zQXsLASaIKXqQ5aaKxih5ArYlg7Yrhy22egI6lS3vs00RYNxmdTeyczeia3\/NSuiN7CoX7LqVF3cfGXcxU+iDABMqFTfOclPFSZQ\/GtvyTSsVscVWT0CnsuV9tmkCjLvwdvHixcW1K+4i2mHDhlX9KYHKao27++7WrVuLce6nAtw1LO4nBL74xS922Us+85nPFHfobX+DvMogAgxNpSOgV9XyAUtPLawytmFOKaOwTVHzz7Hs2jQBxm2mzj\/mOHr0aFm4cKGMGTOm7ULczgHG\/WaSu4vv888\/3+2WvuGGG7iRnb8PSh1hualKhUoohm0CWuAUbAOhEoZhm4AWMMWya1MFmIBtVfoQVmBKJy0KWm4qHbHwqtiGW8WOxDZWLHw8tuFWMSMtuxJgPHsCASamVcLHWm6qcAWdkdjquBK89Vyx1bO1fDwgwBBg9DqnRmXLTdUQ0HZPiq3eFsAWWz0BncqW91kCDAFGp2s8VS03VUNACTB1YWe\/1WPGVsfWsisBhgCj0zUEmIa4shSvy275w0BXzl8dW79RygjLrgQYAkxKT2TPsdxU2TiZBbDNBOTUpx4gtnW3tXw8IMAQYOreUKwS6JJbPmDpyvmrY+s3Sh2Bbapc7XmWXQkwBBidruEUUkNcCYe67JY\/DHTl\/NWx9RuljLDsSoAhwKT0RPYcy02VjZNZANtMQE5z6AFiW3dby8cDAgwBpu4NxSqBLrnlA5aunL86tn6j1BHYpspxCklHzkBVbmSnsxE5WOm4Eg71XLHFVldAp7rlYy0rMKzA6HSNp6rlpmoIaLsnxVZvC2CLrZ6ATmXL+ywBhgCj0zUEmIa4skqgy275w0BXzl8dW79RygjLrgQYAkxKT2TPsdxU2TiZBbDNBKwxHVts9QR0KlveZwkwBBidrmEFpiGurMDoslv+MNCV81fH1m+UMsKyKwGGAJPSE9lzLDdVNk5mAWwzAVmB0QPEtu62lo8HBBgCTN0bilUCXXLLByxdOX91bP1GqSOwTZWrPc+yKwGGAKPTNZxCaogr4VCX3fKHga6cvzq2fqOUEZZdCTAEmJSeyJ5juamycTILYJsJyGkOPUBs625r+XhAgCHA1L2hWCXQJbd8wNKV81fH1m+UOgLbVDlOIenIGajKnXh1NiIHKx1XwqGeK7bY6groVLd8rGUFhhUYna7xVLXcVA0Bbfek2OptAWyx1RPQqWx5nyXAEGB0uoYA0xBXVgl02S1\/GOjK+atj6zdKGWHZlQBDgEnpiew5lpsqGyezALaZgDWmY4utnoBOZcv7LAGGAKPTNazANMSVFRhddssfBrpy\/urY+o1SRlh2JcAQYFJ6InuO5abKxsksgG0mICsweoDY1t3W8vGAAEOAqXtDsUqgS275gKUr56+Ord8odQS2qXK151l2JcAQYHS6hlNIDXElHOqyW\/4w0JXzV8fWb5QywrIrAYYAk9IT2XMsN1U2TmYBbDMBOc2hB4ht3W0tHw8IMASYujcUqwS65JYPWLpy\/urY+o1SR2CbKscpJB05A1W5E6\/ORuRgpeNKONRzxRZbXQGd6paPtazAsAKj0zWeqpabqiGg7Z4UW70tgC22egI6lS3vswQYAoxO1xBgGuLKKoEuu+UPA105f3Vs\/UYpIyy7EmAIMCk9kT3HclNl42QWwDYTsMZ0bLHVE9CpbHmfJcAQYHS6hhWYhriyAqPLbvnDQFfOXx1bv1HKCMuuBBgCTEpPZM+x3FTZOJkFsM0EZAVGDxDbuttaPh4QYAgwdW8oVgl0yS0fsHTl\/NWx9RuljsA2Va72PMuuBBgCjE7XcAqpIa6EQ112yx8GunL+6tj6jVJGWHYlwBBgUnoie47lpsrGySyAbSYgpzn0ALGtu63l4wEBhgBT94ZilUCX3PIBS1fOXx1bv1HqCGxT5TiFpCN3HlU9d+6c\/Ou\/\/qs89thj8sorr8iAAQPka1\/7mkydOlUuuuiibl8pd+LV2YgcrHRcCYd6rthiqyugU93ysbZlVmBeeukl+cY3viFjx46VKVOmyBtvvCFPPfWU3HTTTXLfffdJ3759q+49BBiaSkdAr6rlA5aeWlhlbMOcUkZhm6Lmn2PZtSUCzMmTJ4vw4lZdVq9e3RZWXnzxRZk\/f37x3\/7u7\/6OAOPvhdJGWG6q0pASC2GbCBcwDdsApMQh2CbCeaZZdm2JAONWUdypoieeeEJuueWWts19+vRpWbhwofTv31\/uv\/9+ufDCC7vsCqzA6DTVO++8U6yATZs2Ta688kqdJ2nRqtjqbXhssdUT0KlseZ9tiQDz9NNPy\/r162Xr1q0yfPjwtr3EXRfjVl\/cNTGbN28ugkznBwFGp6lw1XF1VbHFVk9ArzL7rY6tZdeWCDAPPfSQ\/PznP5dNmzbJwIEDO+wl7u9eeOEFefLJJ+Wqq64iwOj0EK51ciXA6EJb\/jDQlfNXx9ZvlDLCsmvLBBi3EautsjzzzDPFN5M6r85UdpTKxl+wYEFxATCPcgTeffddWbx4seBajmf7KtiWb1qpiC22egI6lSv77Pbt22XcuHE6T9KgqgQYT4A5fPhw8UH78ssvN2gT8bQIIIAAAgikC7h\/fK9Zs0Yuv\/zy9CLn4cyWDzC+U0hum7kQ4\/7wQAABBBBAoNkEXHCxFl7cNmiJAJNzEW+z7ai8XgQQQAABBFpBoCUCTM7XqFthJ+A9IoAAAggg0GwCLRFgKjeyu+yyyzrcdTfkRnbNtkF5vQgggAACCLSCQEsEGLch\/\/mf\/1m++c1vyoQJE+SrX\/2qvP3220E\/JdAKOwHvEQEEEEAAgWYTaJkAU+3HHP\/hH\/6h+EHHSy65pNm2G68XAQQQQACBlhZomQDT0luZN48AAggggIAxAQKMsQ3K20EAAQQQQKAVBAgwrbCVeY8IIIAAAggYEyDAdLNB3e2X3U8MPPfcc8WIW2+9tbgj7xVXXGFsFyj\/7bhvfW3ZskW2bdsmZ86cKW5fvWjRIhk5cqT06NGj5hOePXtWdu3aVcw\/cOCAXHzxxTJx4kT51re+hb2I5Nh2ht+9e7fMnj272E7WbjEeu1fnuv76178u7nT605\/+lONFJ\/wc2z\/84Q\/ivi26du1aef3112XAgAEyZcqU4trFfv36xW7mlh3vrgH90Y9+VPwm4D\/+4z9Knz59TFgQYKpsxvfff1\/mz58vR48eLf7Xfeg+9dRT8tvf\/lYef\/xxueaaa0xsfI03cfr0aVm5cqXs2bNH5syZI4MHD5YdO3bI3r17ZcOGDcW3wLp7uLDzne98pwgwd9xxh9x4443y3nvvFT\/C+dFHH8kTTzwh119\/vcbLboqaObad3+A777wjs2bNKr6NZ\/E3UmI2aK6ru8\/UvHnzijudTps2TdyHxaOPPioXXHCBbNy4Ua6++uqYl2NqbK7tzp07ZenSpXLzzTcXweWNN94ojsWf\/\/znZd26dYSYgL3l448\/Fvebf\/fff7\/87d\/+rXz3u98lwAS4NeUQd\/BZv359scHdjz+OGDGieB9HjhwpDlIuvNx7773Sq1evpnx\/2i\/ahY\/ly5cXhm7lxD3cQWzZsmVy6tSpmgedffv2yV133SUrVqyQ22+\/vW21xq2GuQ\/b6667TlatWtWy9jm27be7W+VyB7Mf\/\/jH8rvf\/a7lA0yO64cfflisLrp\/3LjQ8tnPfragfu2112TGjBly5513yty5c70rj9p92aj6ObbHjx8v+n7IkCGyevVq6du3b\/E2KiuH3\/ve9+RLX\/pSo97aef+87rPs4MGD8sgjj8i\/\/Mu\/FK\/3y1\/+MgHmvN9yGS+w0jRf+MIXitWAnj17tlVzpzXcv1bdL1cPHTo041lsTnUrKO5fS+5\/3ZJv5YBTOei4X552pytGjRpVFcAtcbp\/sbqVmuHDh7eNcY3oDmCvvPJK1V8Ut6nZ8V3l2rav5v5V64ynTp1a3NixlVdgcl1fffXVInS78FIJ7M66UtetwjzwwANy0UUXtcJu2uE95toeOnRIpk+fXgRAt\/pSeRw7dkxmzpwp48ePlyVLlrSca+gbrvi5MwrOz\/1D0B2TWYEJFWzCcS6xuqZxp44mT57c4R289NJLxUG\/lQ\/4tTZprQNLLVffbvL73\/++WPV66623WjbAlGXrtoM7tedWDfr371+cqmvl\/TnX1f3OmgvlbrXW3embx\/8L5Nr+5je\/KY7FlesPK9fPVT6Y3YeyCzc8qgv8+7\/\/e3H63QVsd3rTrYy7BwHG8B7jzme7g7o7KHW+XqPyd25JbtKkSYYV0t5ad\/9ictUqf+fc3IdnzMMFF7ccf9tttxX\/4mq\/KhZTp5nHlmFbuR6hd+\/eRSB0qwetHmByXd2v2buLS91KllvZcoHGnSp112y4O39fe+21zbzbZb32XFt37YZbKXQr387ZBZkPPvigOP3proVp9euLYjZOZTWMABOj1oRjKyGl2r9KKw3proXpvDrThG+19Jdcy+fEiRNFCHHfdolZ9q1cUO3+NdfKF1Dn2rrTcM8++2yxUuAO\/FdeeaXU2tdL3znO04I5ru5aIvevWvdtOefrVmDcz5S4D1n3L18XGJ11d6dMz1OS0l5Wjm3lRTjjH\/zgB8Up6crDXZf4\/e9\/v+36xNJesOFCBBjDG7f9WyPApG\/oMg5Y7Z\/dff3SrRS4bzT5vsGU\/qqbY2aurbuo9O677y4ugr7llluKN02A+f+VwWr\/KPGF7sqHgrtQ1X0dfeHChW2rg5WL\/t3Xfh9++OEO14M1xx6X\/ypz91kXAN2Fui54f+UrX2n7VuIPf\/hD+eSTT4ovCvCN0LDtRIAJc2r6UbUO6pxCqr15ax2wYk8huQvO3DeX9u\/fX1wgedNNN7XsNzmceo6tC4LuAuphw4YVKwaVU3AEmDxXd4rDBcIXXnihuLC\/81f83Qesu4XAk08+KVdddVXTHxtj30DOPuuey7m6083ulJELMJVrYCqrsu6boM7Y3SuKR20BAkyL7CFcxJu+oXMv2qs885tvvll84Lp\/Zbl\/vbbqEnz7LZFjWwkqtbasu1+P+xBu\/+2v9D2heWbmuLp36a6Hcze7rGbnbsXgbobZiq7OJtfWXffibrzmTscNHDiww07lrjVy4aVVbWM7jAATK9ak4\/kadfqGy\/3apHvmyrdk3C+Euw+Hz33uc+kvyNDMHFu3T7vrNDo\/3MXRDz74oHz7298u7rHjLjhttV9mz3F1ns8\/\/7zcc889Xe5m7K6JcXfm\/dnPftay31DKtXUBxt2Ft9o3vAgwcQc3AkycV9OOrtzI7p\/+6Z86XDTKjezCNqm7HsDdiK79NSuhN7KrLA27lRd30R732ulonmNbbetxCulTlRxXd1xw1xa503Pum0iVex8dPny4+O+jR48u+uHCCy8MayBjo3Js3Vx3CsmF7GqnkNw9dtxxhp8U8O80BBi\/kZkR7v4D7uDjNrq7YZK7CRU\/JRC2ed31Fu43o9y1K+626u7AXu2nBCoN5e6+65aB3Tj3UwHuX63uJwS++MUvdnnCz3zmM8UdetvfIC\/sVdkYlWrb3WkhAsyn+0Wuq1uFcdcWjRkzpsO3kNwHbCt\/cy7XtvK1\/5\/85CfFbSvc16jdP3Lcqblf\/epXLX9hf8xRjQATo2VgbOcfc\/zrv\/7r4l8Drfy7JqGbtfOPt7l\/hbpvaLgDfOVCvM4Bxl2D4e7i6z4MunvccMMNLXsju4pJii0Bxr\/n5ri6Vdtf\/vKXxcXm7maX7j477p5FX\/\/61\/kB0io\/QBpyPKjss+444a4xct88av\/jru4GdnwDyb9fV0YQYMKtGIkAAggggAACCKgK8GvUqrwURwABBBBAAAENAQKMhio1EUAAAQQQQEBVgACjyktxBBBAAAEEENAQIMBoqFITAQQQQAABBFQFCDCqvBRHAAEEEEAAAQ0BAoyGKjURQAABBBBAQFWAAKPKS3EEEEAAAQQQ0BAgwGioUhMBBBBAAAEEVAUIMKq8FEcAAQQQQAABDQECjIYqNRFAAAEEEEBAVYAAo8pLcQQQQAABBBDQECDAaKhSEwEEEEAAAQRUBQgwqrwURwABBBBAAAENgf8D9f1Yg11oVlcAAAAASUVORK5CYII=","height":337.21943048576213,"width":560}}
%---
%[output:1c89a639]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Some output might be missing due to a network interruption. To get the missing output, rerun the script."}}
%---
%[output:05eacf68]
%   data: {"dataType":"warning","outputData":{"text":"Warning: there is a lot of days with less than 50% data coverage"}}
%---
