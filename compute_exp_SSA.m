function data_cal=compute_exp_SSA(data,lambdaSC, lambdaAE )

%compute the pair-wise exponents and the SSA
%if there is several nephelometers with different wavelengths, lambda is a
%matrix with 4 columns corresponding to the 4 possibility (0, 1, 2 or
%without number)
%takes names
namesOK=data.Properties.VariableNames;
data_cal=timetable(data.Time);
CBs=startsWith(namesOK,'Bs');
namesBs=namesOK(CBs);
if ~isempty(namesBs)
    %take the green scat and the red abs
    a=contains(namesBs,'G');
    if sum(a)==0 && ~isempty(namesBs)
        a=1;
    end
    namesBsG=namesBs(a);

    CBbs=startsWith(namesOK,'Bbs');
    namesBbs=namesOK(CBbs);
    if ~isempty(namesBbs)
        a=contains(namesBbs,'G');
        if sum(a)==0 && ~isempty(namesBbs)
            a=1;
        end
        namesBbsG=namesBbs(a);

        %compute the backscattering fraction

        if length(namesBsG)==1 && length(namesBbsG)==1
            %namesBsG=namesBs(a);
            data_cal.BbsFG=data.(namesBbsG{1})./data.(namesBsG{1});
        else
            Cx(1,:)=contains(namesBsG,'0_')& ~contains(namesBsG,'dry') & ~contains(namesBsG,'Q');
            Cx(2,:)=contains(namesBsG,'1_')& ~contains(namesBsG,'dry') & ~contains(namesBsG,'Q');
            Cx(3,:)=contains(namesBsG,'dry');
            Cx(4,:)=contains(namesBsG,'_');
            Cx(4,:)=Cx(4,:) & ~Cx(1,:) & ~Cx(2,:) & ~Cx(3,:);

            Cy(1,:)=contains(namesBbsG,'0_')& ~contains(namesBbsG,'dry') & ~contains(namesBbsG,'Q');
            Cy(2,:)=contains(namesBbsG,'1_')& ~contains(namesBbsG,'dry') & ~contains(namesBbsG,'Q');
            Cy(3,:)=contains(namesBbsG,'dry');
            Cy(4,:)=contains(namesBbsG,'_');
            Cy(4,:)=Cy(4,:) & ~Cy(1,:) & ~Cy(2,:) & ~Cy(3,:);

            e=sum(Cx,2); r=sum(Cy,2);
            for i=1:length(e)
                if e(i)==1 & r(i)==1
                    namesBsx=namesBsG(Cx(i,:));
                    namesBbsx=namesBbsG(Cy(i,:));

                    if i==4
                        namesFG=strcat('BbsFG');
                    else
                        namesFG=strcat('BbsFG', num2str(i-1));
                    end
                    data_cal.(namesFG)=data.(namesBbsx{1})./data.(namesBsx{1});
                end
            end
        end
    else
        namesBsG=namesBs;
    end

    %compute the scattering exponent
    Cs(1,:)=contains(namesBs,'0_')& ~contains(namesBs,'dry') & ~contains(namesBs,'Q');
    Cs(2,:)=contains(namesBs,'1_')& ~contains(namesBs,'dry') & ~contains(namesBs,'Q');
    Cs(3,:)=contains(namesBs,'dry')& ~contains(namesBs,'Q');
    Cs(4,:)=contains(namesBs,'_')& ~contains(namesBs,'Q');
    Cs(4,:)=Cs(4,:) & ~Cs(1,:) & ~Cs(2,:) & ~Cs(3,:);
    for i=1:4
        namesBsx=namesBs(Cs(i,:));
        %ATTENTION: si il existe BsQ ne pas en tenir compte
        if length(namesBsx)>=3
            b=namesBsx(startsWith(namesBsx,'BsB'));
            g=namesBsx(startsWith(namesBsx,'BsG'));
            r=namesBsx(startsWith(namesBsx,'BsR'));
            N1=strcat('expS_bg', num2str(i-1));
            N2=strcat('expS_br', num2str(i-1));
            N3=strcat('expS_gr', num2str(i-1));
            data_cal.(N1)=real(-log(data.(b{1})./data.(g{1}))/log(lambdaSC(1)/lambdaSC(2)));
            data_cal.(N2)=real(-log(data.(b{1})./data.(r{1}))/log(lambdaSC(1)/lambdaSC(3)));
            data_cal.(N3)=real(-log(data.(g{1})./data.(r{1}))/log(lambdaSC(2)/lambdaSC(3)));
            %data_cal.exp_fit=Exp_SC_nan_regr(
        elseif length(namesBsx)==2
            b=namesBsx(startsWith(namesBsx,'BsB'));
            g=namesBsx(startsWith(namesBsx,'BsG'));
            N1=strcat('expS_bg', num2str(i-1));
            data_cal.(N1)=real(-log(data.(b{1})./data.(g{1}))/log(lambdaSC(1)/lambdaSC(2)));
        end
    end
else
    namesBsG=namesBs;
end

% absorption
CBa=startsWith(namesOK,'Ba');
namesBa=namesOK(CBa);
if ~isempty(namesBa)
    a=contains(namesBa,'R') |contains(namesBa,'5');
    if sum(a)==0 && ~isempty(namesBa)
        a=1;
    end
    namesBaR=namesBa(a);

    %compute the absorption exponent
    %ATTENTION: il peut y avoir AE + un a 3lambda
    ax=contains(namesBa,'Ba7') ;namesBa7=namesBa(ax);
    if  ~isempty(namesBa7) & length(namesBa)==7
        b=namesBa(startsWith(namesBa,'Ba2'));
        g=namesBa(startsWith(namesBa,'Ba3'));
        r=namesBa(startsWith(namesBa,'Ba5'));
        data_cal.expA_bg=real(-log(data.(b{1})./data.(g{1}))/log(470/520));
        data_cal.expA_br=real(-log(data.(b{1})./data.(r{1}))/log(470/660));
        data_cal.expA_gr=real(-log(data.(g{1})./data.(r{1}))/log(520/660));
        a1=namesBa(startsWith(namesBa,'Ba1'));
        a4=namesBa(startsWith(namesBa,'Ba4'));
        a6=namesBa(startsWith(namesBa,'Ba6'));
        a7=namesBa(startsWith(namesBa,'Ba7'));
        all_Abs=[data.(a1{1}) data.(b{1}) data.(g{1}) data.(a4{1}) data.(r{1}) data.(a6{1}) data.(a7{1})];
        data_cal.expA_fit=Exp_AE_nan_regr(all_Abs);
    end
    if  length(namesBa)==10 || length(namesBa)==13 || length(namesBa)==16
        warning(' one AE and one/several 3w instrument?');
        a1=namesBa(startsWith(namesBa,'Ba1'));
        b=namesBa(startsWith(namesBa,'Ba2'));
        g=namesBa(startsWith(namesBa,'Ba3'));
        a4=namesBa(startsWith(namesBa,'Ba4'));
        r=namesBa(startsWith(namesBa,'Ba5'));
        a6=namesBa(startsWith(namesBa,'Ba6'));
        a7=namesBa(startsWith(namesBa,'Ba7'));
        data_cal.expA_AE_bg=real(-log(data.(b{1})./data.(g{1}))/log(470/520));
        data_cal.expA_AE_br=real(-log(data.(b{1})./data.(r{1}))/log(470/660));
        data_cal.expA_AE_gr=real(-log(data.(g{1})./data.(r{1}))/log(520/660));
        all_Abs=[data.(a1{1}) data.(b{1}) data.(g{1}) data.(a4{1}) data.(r{1}) data.(a6{1}) data.(a7{1})];
        data_cal.expA_AE_fit=Exp_AE_nan_regr(all_Abs);

        Ca(1,:)=contains(namesBa,'0_');
        Ca(2,:)=contains(namesBa,'1_');
        Ca(3,:)=contains(namesBa,'dry');
        Ca(4,:)=contains(namesBa,'_');
        Ca(4,:)=Ca(4,:) & ~Ca(1,:) & ~Ca(2,:) & ~Ca(3,:);
        for i=1:4
            namesBax=namesBa(Ca(i,:));
            %ATTENTION: changer les noms sinon ils se récrivent l'un sur
            %l'autre
            if length(namesBax)==3
                ba=namesBax(startsWith(namesBax,'BaB'));
                ga=namesBax(startsWith(namesBax,'BaG'));
                ra=namesBax(startsWith(namesBax,'BaR'));
                N1=strcat('expA_bg', num2str(i-1));
                N2=strcat('expA_br', num2str(i-1));
                N3=strcat('expA_gr', num2str(i-1));
                data_cal.(N1)=real(-log(data.(ba{1})./data.(ga{1}))/log(lambdaAE(1)/lambdaAE(2)));
                data_cal.(N2)=real(-log(data.(ba{1})./data.(ra{1}))/log(lambdaAE(1)/lambdaAE(3)));
                data_cal.(N3)=real(-log(data.(ga{1})./data.(ra{1}))/log(lambdaAE(2)/lambdaAE(3)));
            end
        end
    elseif length(namesBa)==12
        %2 x 3w + 2 size cuts
        Ca(1,:)=contains(namesBa,'0_A11');
        Ca(2,:)=contains(namesBa,'1_A11');
        Ca(3,:)=contains(namesBa,'0_A12');
        Ca(4,:)=contains(namesBa,'1_A12');
        % Ca(4,:)=Ca(4,:) & ~Ca(1,:) & ~Ca(2,:) & ~Ca(3,:);
        for i=1:4
            namesBax=namesBa(Ca(i,:));
            %ATTENTION: changer les noms sinon ils se récrivent l'un sur
            %l'autre
            if length(namesBax)==3
                ba=namesBax(startsWith(namesBax,'BaB'));
                ga=namesBax(startsWith(namesBax,'BaG'));
                ra=namesBax(startsWith(namesBax,'BaR'));
                N1=strcat('expA_bg', num2str(i-1));
                N2=strcat('expA_br', num2str(i-1));
                N3=strcat('expA_gr', num2str(i-1));
                data_cal.(N1)=real(-log(data.(ba{1})./data.(ga{1}))/log(lambdaAE(1)/lambdaAE(2)));
                data_cal.(N2)=real(-log(data.(ba{1})./data.(ra{1}))/log(lambdaAE(1)/lambdaAE(3)));
                data_cal.(N3)=real(-log(data.(ga{1})./data.(ra{1}))/log(lambdaAE(2)/lambdaAE(3)));
            end
        end
    else
        Ca(1,:)=contains(namesBa,'0_');
        Ca(2,:)=contains(namesBa,'1_');
        Ca(3,:)=contains(namesBa,'dry');
        Ca(4,:)=contains(namesBa,'_');
        Ca(4,:)=Ca(4,:) & ~Ca(1,:) & ~Ca(2,:) & ~Ca(3,:);
        for i=1:4
            namesBax=namesBa(Ca(i,:));
            %ATTENTION: changer les noms sinon ils se récrivent l'un sur
            %l'autre
            if length(namesBax)==3
                ba=namesBax(startsWith(namesBax,'BaB'));
                ga=namesBax(startsWith(namesBax,'BaG'));
                ra=namesBax(startsWith(namesBax,'BaR'));
                N1=strcat('expA_bg', num2str(i-1));
                N2=strcat('expA_br', num2str(i-1));
                N3=strcat('expA_gr', num2str(i-1));
                data_cal.(N1)=real(-log(data.(ba{1})./data.(ga{1}))/log(lambdaAE(1)/lambdaAE(2)));
                data_cal.(N2)=real(-log(data.(ba{1})./data.(ra{1}))/log(lambdaAE(1)/lambdaAE(3)));
                data_cal.(N3)=real(-log(data.(ga{1})./data.(ra{1}))/log(lambdaAE(2)/lambdaAE(3)));
            end
        end
    end
end

if ~isempty(namesBsG) && ~isempty(namesBa)
    %compute the SSA with, if possible the green scat and the red abs
    %or with any wavelength, since only the trend is important
    CaR(1,:)=startsWith(namesBaR,'BaR0');
    CaR(2,:)=startsWith(namesBaR,'BaR1');
    CaR(3,:)=startsWith(namesBaR,["BacR";"Ba3"]) ;
    CaR(4,:)=startsWith(namesBaR,["BaR_";"Ba5"]);
    CaG=startsWith(namesBa,'BaG_');

    CsG(1,:)=startsWith(namesBsG,'BsG0');
    CsG(2,:)=startsWith(namesBsG,'BsG1');
    CsG(3,:)=startsWith(namesBsG,'BsG') & endsWith(namesBsG,'dry');
    CsG(4,:)=startsWith(namesBsG,'BsG_');
    for i=1:4
        if nansum(CaR(i,:))>0 && nansum(CsG(i,:))>0
            namesBsy=namesBsG(CsG(i,:));
            namesBay=namesBaR(CaR(i,:));
            for j=1:length(namesBsy)
                for k=1:length(namesBay)
                    N_SSA=strcat('SSA', num2str(i-1),num2str(j),num2str(k));
                    data_cal.(N_SSA)=data.(namesBsy{j})./(data.(namesBsy{j})+data.(namesBay{k}));
                end
            end

        elseif nansum(CaG)>0 && nansum(CsG(i,:))>0
            namesBsy=namesBsG(CsG(i,:));
            namesBay=namesBa(CaG);
            for j=1:length(namesBsy)
                for k=1:length(namesBay)
                    N_SSA=strcat('SSA', num2str(i-1),num2str(j),num2str(k));
                    data_cal.(N_SSA)=data.(namesBsy{j})./(data.(namesBsy{j})+data.(namesBay{k}));
                end
            end
        end
    end
end


end