ATTO_st.name='ATTO';
ATTO_st.lat=2.1483;
ATTO_st.lon=59.0033;
ATTO_st.alt=130+325;
ATTO_st.env='con';
ATTO_st.footp='f';
%%
ATTO_rd=read_betsy_2026('C:\github_trend\raw_data\ATTO_maap_aeth_neph_house',ATTO_st.name);
%data not ok with clap and psap instead of Ae33 ans maap
names_rd=fieldnames(ATTO_rd);
%%


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
