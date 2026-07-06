function data_cal=compute_exp_D(data,shortnames,lambda)


% if scattering data are given, compute the pair-wise expSC and BF
% if absorption data are given, compute the expAE, pair-wise for 3w CLAP or
% PSAP & 2wAE
% and by a fit for  7w AE

% % % %if there is several nephelometers with different wavelengths, lambda is a
% % % %matrix with 4 columns corresponding to the 4 possibility (0, 1, 2 or
% % % %without number)

%takes names
namesOK=data.Properties.VariableNames;

namesB=namesOK(contains(namesOK,shortnames));

namesBs=namesB(contains(namesB,{'Bs'}));
namesBbs=namesB(contains(namesB,{'Bbs'}));
namesBa=namesB(contains(namesB,{'Ba'}));


% initiate timetable
data_cal=timetable(data.Time);

% compute expS and BF
if ~isempty(namesBs)
    if length(namesBs)== 1 & length(namesBbs)==1
        data_cal.BbsFG=data.(namesBbs{1})./data.(namesBs{1});
        if ~isempty(namesBbs) & length(namesBbs)~=1
            error('there is more backscat coef than scat coeff ')
        end
    elseif length(namesBs)== 3
        b=namesBs(startsWith(namesBs,'BsB'));
        g=namesBs(startsWith(namesBs,'BsG'));
        r=namesBs(startsWith(namesBs,'BsR'));
        N1='expS_bg';
        N2='expS_br';
        N3='expS_gr';
        data_cal.(N1)=real(-log(data.(b{1})./data.(g{1}))/log(lambda(1)/lambda(2)));
        data_cal.(N2)=real(-log(data.(b{1})./data.(r{1}))/log(lambda(1)/lambda(3)));
        data_cal.(N3)=real(-log(data.(g{1})./data.(r{1}))/log(lambda(2)/lambda(3)));
        if ~isempty(namesBbs)
            Bb=namesBbs(startsWith(namesBbs,'BbsB'));
            Bg=namesBbs(startsWith(namesBbs,'BbsG'));
            Br=namesBbs(startsWith(namesBbs,'BbsR'));

            data_cal.BsBFb=data.(Bb{1})./data.(b{1});
            data_cal.BsBFg=data.(Bg{1})./data.(g{1});
            data_cal.BsBFr=data.(Br{1})./data.(r{1});
        end
    elseif length(namesBs)== 6
        Cs(1,:)=contains(namesBs,'0_')& ~contains(namesBs,'dry') & ~contains(namesBs,'Q');
        Cs(2,:)=contains(namesBs,'1_')& ~contains(namesBs,'dry') & ~contains(namesBs,'Q');

        for i=1:2
            namesBsx=namesBs(Cs(i,:));
            b=namesBsx(startsWith(namesBsx,'BsB'));
            g=namesBsx(startsWith(namesBsx,'BsG'));
            r=namesBsx(startsWith(namesBsx,'BsR'));
            N1=strcat('expS_bg', num2str(i-1));
            N2=strcat('expS_br', num2str(i-1));
            N3=strcat('expS_gr', num2str(i-1));
            data_cal.(N1)=real(-log(data.(b{1})./data.(g{1}))/log(lambda(1)/lambda(2)));
            data_cal.(N2)=real(-log(data.(b{1})./data.(r{1}))/log(lambda(1)/lambda(3)));
            data_cal.(N3)=real(-log(data.(g{1})./data.(r{1}))/log(lambda(2)/lambda(3)));
        end
    end

    if ~isempty(namesBbs)
        Cbs(1,:)=contains(namesBbs,'0_');
        Cbs(2,:)=contains(namesBbs,'1_');
        for i=1:2
            namesBbsx=namesBbs(Cbs(i,:));
            Bb=namesBbsx(startsWith(namesBbsx,'BbsB'));
            Bg=namesBbsx(startsWith(namesBbsx,'BbsG'));
            Br=namesBbsx(startsWith(namesBbsx,'BbsR'));

            namesBsx=namesBs(Cs(i,:));
            b=namesBsx(startsWith(namesBsx,'BsB'));
            g=namesBsx(startsWith(namesBsx,'BsG'));
            r=namesBsx(startsWith(namesBsx,'BsR'));

            N1=strcat('BbsFb', num2str(i-1));
            N2=strcat('BbsFbg', num2str(i-1));
            N3=strcat('BbsFbr', num2str(i-1));
            data_cal.(N1)=data.(Bb{1})./data.(b{1});
            data_cal.(N2)=data.(Bg{1})./data.(g{1});
            data_cal.(N3)=data.(Br{1})./data.(r{1});
        end

    end
end


% absorption
if ~isempty(namesBa)
    %compute the absorption exponent either for
    % 2w AE
    % 3w CLAP or PSAP
    % 2*3w CLAP or PSAP
    % 7w AE
    if length(namesBa)==2
        if length(lambda)~=2
            error('the nb of wavelengths is different from the number of abs coef')
        end
        N='expA';
        data_cal.(N)=real(-log(data.(namesBa{1})./data.(namesBa{2}))/log(lambda(1)/lambda(2)));

    elseif length(namesBa)==3 % 3w PSAP or CLAP
        if length(lambda)~=3
            error('the nb of wavelengths is different from the number of abs coef')
        end
        ba=namesBa(startsWith(namesBa,'BaB'));
        ga=namesBa(startsWith(namesBa,'BaG'));
        ra=namesBa(startsWith(namesBa,'BaR'));

        data_cal.expA_bg=real(-log(data.(ba{1})./data.(ga{1}))/log(lambda(1)/lambda(2)));
        data_cal.expA_br=real(-log(data.(ba{1})./data.(ra{1}))/log(lambda(1)/lambda(3)));
        data_cal.expA_gr=real(-log(data.(ga{1})./data.(ra{1}))/log(lambda(2)/lambda(3)));


    elseif length(namesBa)==6 %2*3w instruments
        if length(lambda)~=3
            error('the nb of wavelengths is different from the number of abs coef')
        end
        Ca(1,:)=contains(namesBa,'0_');
        Ca(2,:)=contains(namesBa,'1_');
        for i=1:2
            namesBax=namesBa(Ca(i,:));
            ba=namesBax(startsWith(namesBax,'BaB'));
            ga=namesBax(startsWith(namesBax,'BaG'));
            ra=namesBax(startsWith(namesBax,'BaR'));
            N1=strcat('expA_bg', num2str(i-1));
            N2=strcat('expA_br', num2str(i-1));
            N3=strcat('expA_gr', num2str(i-1));
            data_cal.(N1)=real(-log(data.(ba{1})./data.(ga{1}))/log(lambda(1)/lambda(2)));
            data_cal.(N2)=real(-log(data.(ba{1})./data.(ra{1}))/log(lambda(1)/lambda(3)));
            data_cal.(N3)=real(-log(data.(ga{1})./data.(ra{1}))/log(lambda(2)/lambda(3)));
        end


    elseif length(namesBa)==7 % 7w AE
        if length(lambda)~=7
            error('the nb of wavelengths is different from the number of abs coef')
        end

        b=namesBa(startsWith(namesBa,'Ba2'));
        g=namesBa(startsWith(namesBa,'Ba3'));
        r=namesBa(startsWith(namesBa,'Ba5'));
        data_cal.expA_bgAE=real(-log(data.(b{1})./data.(g{1}))/log(470/520));
        data_cal.expA_brAE=real(-log(data.(b{1})./data.(r{1}))/log(470/660));
        data_cal.expA_grAE=real(-log(data.(g{1})./data.(r{1}))/log(520/660));
        a1=namesBa(startsWith(namesBa,'Ba1'));
        a4=namesBa(startsWith(namesBa,'Ba4'));
        a6=namesBa(startsWith(namesBa,'Ba6'));
        a7=namesBa(startsWith(namesBa,'Ba7'));
        all_Abs=[data.(a1{1}) data.(b{1}) data.(g{1}) data.(a4{1}) data.(r{1}) data.(a6{1}) data.(a7{1})];
        data_cal.expA_fit=Exp_AE_nan_regr(all_Abs);
    else
        warning('the number of abs coef does not allow to understand the instrument type')
    end



end

