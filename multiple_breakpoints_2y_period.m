function result = multiple_breakpoints_2y_period(data,param, nb_data_min, alpha)
% multiple_breakpoints iteratively search for breakpoints in a time series
% (data) with a defined method (e.g. pettitt). The minimal used period and/or the number of iterations can
% be given.
%   Input: data
%           method
%         nb_data_min
%         nb_iteration
%
% Outpur: results
%
% example:
% mco, March 2026


arguments
    data {mustBeUnderlyingType(data,"timetable")}
    %method {mustBeText}
    param {mustBeText}
    nb_data_min {mustBeNumeric}
    alpha {mustBeNumeric}
end

% arguments
%     result {mustBeUnderlyingType(result,"struct")}
%
% end

% while: the process is done until each segment< nb_data_min data points
% the pre-allocation therefore needs nb data points/(nb_data_min/2)
nbboot=10000;
nb_2yper=year(data.Time(end))-year(data.Time(1));
nb_preallocation=nb_2yper; % the unused data will be removed


%remove potential NAN
%for j=1:length(param{1})
% result.(param{j})=struct("level",NaN(nb_preallocation,1),"pvalue",NaN(nb_preallocation,1),"time",NaT(nb_preallocation,1));
result=struct("level",NaN(nb_preallocation,1),"pvalue",NaN(nb_preallocation,1),"time",NaT(nb_preallocation,1),"pvalue_boot",NaN(nb_preallocation,1),"PrctDiff",NaN(nb_preallocation,3));

data(isnan(data.(param{1})),:)=[];
data.y=year(data.Time);
for i=1:nb_2yper
    data_2yper=data(data.y>=(data.y(1)+(i-1)) & data.y<(data.y(1)+(i-1)+2),:);
    if length( data_2yper.(param{1}))>=nb_data_min

        x=data_2yper.(param{1});
        [a,PrctDiff]=pettitt(x, alpha);
        % result.(param{j}).level(k)=i;
        % result.(param{j}).pvalue(k)=a(3);
        result.pvalue(i)=a(3);
        result.PrctDiff(i,:)=PrctDiff;
        %compute p value from bootstrapping
        test_boot=bootstrp(nbboot,@pettitt,x,alpha);
        indBoot=test_boot(:,2)>=a(2);
        result.pvalue_boot(i)=(1+sum(indBoot))/(nbboot+1);
        if ~isnan(a(1))
            result.time(i)=data_2yper.Time(a(1));
        end
    end
end

%remove empty results
if sum(isnan(result.pvalue))==length(result.pvalue)
    result.pvalue(2:end,:)=[];
    result.pvalue_boot(2:end,:)=[];
    result.level(2:end,:)=[];
    result.time(2:end,:)=[];
    result.PrctDiff(2:end,:)=[];
else
    ind=isnan(result.pvalue);
    result.pvalue(ind,:)=[];
    result.pvalue_boot(ind,:)=[];
    result.level(ind,:)=[];
    result.time(ind,:)=[];
    result.PrctDiff(ind,:)=[];
end


