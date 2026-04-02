function T=make_table_breakpoints(Break_result,param)
% table from break points
nb_param=length(param);
nb_col=8;
nb_row=height(Break_result);
Vartype={'string','datetime', 'double','double','double','double','double','double'};
VarNom={'param','Break_point','p_value','p_value_boot','level','minDiff','medDiff','maxDiff'};

T=table('Size',[nb_row nb_col],'VariableTypes',Vartype,'VariableNames',VarNom);
k=0;
for i=1:height(Break_result)
    if ~isempty(Break_result.station{i})
    for j=1:height(Break_result.results{i}.time)
        k=k+1;
    T.param(k)=Break_result.parameter{i};
     T.Break_point(k)=Break_result.results{i}.time(j);
      T.p_value(k)=Break_result.results{i}.pvalue(j);
      T.p_value_boot(k)=Break_result.results{i}.pvalue_boot(j);
      T.minDiff(k)=Break_result.results{i}.PrctDiff(j,1);
      T.medDiff(k)=Break_result.results{i}.PrctDiff(j,2);
      T.maxDiff(k)=Break_result.results{i}.PrctDiff(j,3);
       T.level(k)=Break_result.results{i}.level(j);

    end
    end
end

T=sortrows(T,'Break_point'); 


