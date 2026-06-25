function data_cal=compute_SSA_D(data,shortnamesSC,shortnamesabs)

%compute the SSA from the given shortnames

%if there is several nephelometers with different wavelengths, lambda is a
%matrix with 4 columns corresponding to the 4 possibility (0, 1, 2 or
%without number)
%takes names
namesOK=data.Properties.VariableNames;
% Cn=startsWith(shortnames,{'Bs'});
% Ca=startsWith(shortnames,{'Ba'});


namesB=namesOK(contains(namesOK,shortnamesSC));
namesBs=namesB(contains(namesB,{'Bs'}));

namesB=namesOK(contains(namesOK,shortnamesabs));
namesBa=namesB(contains(namesB,{'Ba'}));


% initiate timetable
data_cal=timetable(data.Time);

% compute SSAB
if ~isempty(namesBs) && ~isempty(namesBa)
    % compute SSAB
    if length(namesBs)==1 & length(namesBa)==1
data_cal.SSA=data.(namesBs{1})./(data.(namesBs{1})+data.(namesBa{1}));
    elseif length(namesBs)==3 & length(namesBa)==3 %(trhee wavelength + 1 size cuts)
       % 
       % CaB(1,:)=startsWith(namesBa,["Ba2","Bac2","Bax2",BaB]);
       % CsB(1,:)=startsWith(namesBa,'BaB');
data_cal.SSAB=data.(namesBs{startsWith(namesBs,'BsB')})./(data.(namesBs{startsWith(namesBs,'BsB')})+data.(namesBa{startsWith(namesBa,["Ba2","Bac2","Bax2","BaB"])}));
data_cal.SSAG=data.(namesBs{startsWith(namesBs,'BsG')})./(data.(namesBs{startsWith(namesBs,'BsG')})+data.(namesBa{startsWith(namesBa,["Ba3","Bac3","Bax3","BaG"])}));
data_cal.SSAR=data.(namesBs{startsWith(namesBs,'BsR')})./(data.(namesBs{startsWith(namesBs,'BsR')})+data.(namesBa{startsWith(namesBa,["Ba5","Bac5","Bax5","BaR"])}));

   elseif length(namesBs)==6 & length(namesBa)==6 %(trhee wavelength + two size cuts)

data_cal.SSAB0=data.(namesBs{startsWith(namesBs,["BsB_","BsB0"])})./(data.(namesBs{startsWith(namesBs,["BsB_","BsB0"])})+data.(namesBa{startsWith(namesBa,["Ba20","Bac20","Bax20","BaB0","Ba2_","Bac2_","Bax2_","BaB_"])}));
data_cal.SSAG0=data.(namesBs{startsWith(namesBs,["BsG_","BsG0"])})./(data.(namesBs{startsWith(namesBs,["BsG_","BsG0"])})+data.(namesBa{startsWith(namesBa,["Ba30","Bac30","Bax30","BaG0","Ba3_","Bac3_","Bax3_","BaG_"])}));
data_cal.SSAR0=data.(namesBs{startsWith(namesBs,["BsR_","BsR0"])})./(data.(namesBs{startsWith(namesBs,["BsR_","BsR0"])})+data.(namesBa{startsWith(namesBa,["Ba50","Bac50","Bax50","BaR0","Ba5_","Bac5_","Bax5_","BaR_"])}));

data_cal.SSAB1=data.(namesBs{startsWith(namesBs,"BsB1")})./(data.(namesBs{startsWith(namesBs,"BsB1")})+data.(namesBa{startsWith(namesBa,["Ba21","Bac21","Bax21","BaB1"])}));
data_cal.SSAG1=data.(namesBs{startsWith(namesBs,"BsG1")})./(data.(namesBs{startsWith(namesBs,"BsG1")})+data.(namesBa{startsWith(namesBa,["Ba31","Bac31","Bax31","BaG10"])}));
data_cal.SSAR1=data.(namesBs{startsWith(namesBs,"BsR1")})./(data.(namesBs{startsWith(namesBs,"BsR1")})+data.(namesBa{startsWith(namesBa,["Ba51","Bac51","Bax51","BaR1"])}));


    end
end
   

