CBW_st.name='CBW'; 
CBW_st.lat=51.9710006714;
CBW_st.lon=4.9270000458;
CBW_st.CBW=0;
%%

CBW_rd=read_betsy('CBW_2008_2018',CBW_st.name);


%%
datevec(CBW_rd.Time(1))
datevec(CBW_rd.Time(end))
prctile(CBW_rd.U_S,[5 50 95])
%  prctile(CBW_rd.U0_S,[5 50 95])
%  prctile(CBW_rd.U1_S,[5 50 95])
% prctile(CBW_rd.Uu1_S,[5 50 95])
plotFigControl(CBW_rd,CBW_st.name);
% ? size cut, ok
% sc: ok, ? negatives,? clear seasonal cycle
%%
lambdaSC=[450;550;700]*ones(1,4);
lambdaAE=[467;530;660]*ones(1,4);
lambdaAE7=[370 470 520 590 660 880 950];
CBW_cal=compute_exp_SSA(CBW_rd,lambdaSC, lambdaAE);
plotFigControl_cal(CBW_cal, CBW_st.name);

%%
%Questions:
%problems:


%%

CBW_t=timetable(CBW_rd.Time);
% begin at the beginning of a year:
%end: 



%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
