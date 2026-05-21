function x=Exp_AE_nan_regr (dataAE)

% 1. remplace les lignes des coefficients comportant une valeur negative par
% des NAN, car le fit ne supporte pas les valeurs negatives
% 2. fitte par une loi exponentielle le coefficient de diffusion
% selon la formule data(lambda) = b * lambda^(-a) 

% dataSC == coefficient de diffusion pour les longueurs d'ondes 450 nm 550 nm et 700 nm
%   x (:,1) = exposant a

% Martine Collaud Coen, MeteoSwiss, Aout 2003


lambdaAE = [370 470 520 590 660 880 950];
LongueurFichier = size(dataAE,1);

%variables du nombre de points exclus pour le fit
SCexclus=0;


% si un des coefficients de diffusion  est negatif, il est remplacee par NaN,
% car les fits ne supportent pas les valeurs negatives
for i=1:LongueurFichier    
    for j=1:7       
        if dataAE(i,j)<0 
            dataAE(i,j)=NaN;
            SCexclus=SCexclus+1;
        end
    end
end

%affiche a l'ecran le nombre de valeurs negatives trouvees
SCexclus;

%calcule l'exposant du coefficient d'absorption
x=fitlambda_rob_nan(dataAE,lambdaAE);
x=real(x(:,1));




% Martine Collaud Coen
% MeteoSwiss
% CH-1560 Payerne
% 28.1.2008


