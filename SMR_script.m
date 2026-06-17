SMR_st.name='SMR'; 
SMR_st.lat=61.84738441;
SMR_st.lon=24.2947798;
SMR_st.SMR=181;
%%

SMR_rd=read_betsy('SMR_2008_2018',SMR_st.name);


%%
datevec(SMR_rd.Time(1))
datevec(SMR_rd.Time(end))
prctile(SMR_rd.U_S,[5 50 95])
%  prctile(SMR_rd.U0_S,[5 50 95])
%  prctile(SMR_rd.U1_S,[5 50 95])
% prctile(SMR_rd.Uu1_S,[5 50 95])
plotFigControl(SMR_rd,SMR_st.name);
% ? size cut, ok
% sc: ok, ? negatives,? clear seasonal cycle
%%
lambdaSC=[450;550;700]*ones(1,4);
lambdaAE=[467;530;660]*ones(1,4);
lambdaAE7=[370 470 520 590 660 880 950];
SMR_cal=compute_exp_SSA(SMR_rd,lambdaSC, lambdaAE);
plotFigControl_cal(SMR_cal, SMR_st.name);

%%
%Questions:
%problems:


%%

SMR_t=timetable(SMR_rd.Time);
% begin at the beginning of a year:
%end: 



%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
