BBE_st.name='BBE'; 
BBE_st.lat=29.30;
BBE_st.lon=-103.18;
BBE_st.alt=1066;

%%
BBE_rd=read_betsy('BBE_head',BBE_st.name);
%%
datevec(BBE_rd.Time(1))
datevec(BBE_rd.Time(end))
prctile(BBE_rd.U_S,[5, 50, 95])
plotFigControl(BBE_rd,BBE_st.name);
% 1 size cut
% sc: ok, clear seasonal cycle in summer, 

%%
% lambdaSC=[450;550;700]*ones(1,4);
% lambdaAE3=[467;530;660]*ones(1,4);
% lambdaAE7=[370 470 520 590 660 880 950];
% BBE_cal=compute_exp_SSA(BBE_rd,lambdaSC, lambdaAE3);
% plotFigControl_cal(BBE_cal, BBE_st.name);
%%
%Questions:
%problems:
% MINIMA higher in summer 2011 ?
%modification in 2016 ? rupture ?

%%
BBE_t=timetable(BBE_rd.Time);
% begin at the beginning of a year: 1998
%end 2018

%nothing in Bbsx_: to delete, but Bsx exist, the green canal being the longest. Why no other w in 1983-1984 ?
