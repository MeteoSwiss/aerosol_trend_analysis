% This code is used to find the change point in a univariate continuous time series
% using Pettitt Test.
%
%
% The test here assumed is two-tailed test. The hypothesis are as follow:
%  H (Null Hypothesis): There is no change point in the series
%  H(Alternative Hypothesis): There is a change point in the series
%
% Input: univariate data series
% Output:
% The output of the answer in row wise respectively,
% loc: location of the change point in the series, index value in
% the data set
% K: Pettitt Test Statistic for two tail test
% pvalue: p-value of the test
%
%Reference: Pohlert, Thorsten. "Non-Parametric Trend Tests and Change-Point Detection." (2016).
%
function [a, PrctDiff]=pettitt(data, alpha)
[m, ~]=size(data);
%without vectorisation, quite low
% % for t=2:1:m
% %     for j=1:1:m
% %       v(t-1,j)=sign(data(t-1,1)-data(j,1));
% %       V(t-1)=sum(v(t-1,:));
% %     end
% % end

%with vectorisation
t1 = repmat(data',m,1);
t1=t1';
v = sign(t1 - data');
V = sum(v');
%just the same thereafter
U=cumsum(V);
loc=find(abs(U)==max(abs(U)));
if length(loc)>=2
d=length(data)*ones(1,length(loc))-loc;
[~,ind_min]=find(min(d));
loc=loc(ind_min);
end

K=max(abs(U));
pvalue=2*exp((-6*K^2)/(m^3+m^2));



if isempty(loc)
    loc=NaN;
end
if pvalue<=alpha
    a=[loc; K ;pvalue];
    % compute difference after-before the break for min (10%), med (50%) and max (90%) percentiles
PrctDiff=round((prctile(data(loc:end),[10 50 90])-prctile(data(1:loc),[10 50 90])).*100./prctile(data,[10 50 90]));
else
    a=[NaN; K ;NaN];
    PrctDiff=[NaN, NaN, NaN];
end
return