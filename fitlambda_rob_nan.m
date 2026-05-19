function [y,y_err]=fitlambda_rob_nan (data,lambda)

% fit par une loi exponantielle les donnees de data
% selon la formule data(lambda) = b * lambda^(-a) et retourne a et b
% data == les donnees a fitter sans les vecteurs temps
% lambda == vecteur contenant les longueurs d'onde
% y(:,1)== a, y(:,2) == b
%le fit est fait de manière robust et l'erreur sur les paramètres est
%données

LongueurFichier=size(data,1);
nbLOnde=size(lambda,2);
for i=1:LongueurFichier
    if sum(isnan(data(i,:)))==0
        [P,stats]= robustfit(log10(lambda),log10(data(i,1:nbLOnde)));
        y(i,1)= real(-P(2));
        y_err(i,1)=2* stats.robust_s;
        y(i,2)=real(10.^(P(1)));
        %y_err1(i,2)=10.^(P(1)+2*stats.robust_s(1));
        %y_err2(i,2)=10.^(P(1)-2*stats.robust_s(1));
        %y_err(i,2)=max(y_err1(i,2), y_err1(i,2));
%     elseif sum(~isnan(data(i,:)))>= 2
%         ind=~isnan(data(i,:));
%         [P,stats]= robustfit(log10(lambda(:,ind)),log10(data(i,ind)));
%         y(i,1)=-real(P(2));
%         y_err(i,1)=2* stats.robust_s;
%         y(i,2)=10.^(P(1));
        %y_err1(i,2)=10.^(P(1)+2*stats.robust_s(1));
        % y_err2(i,2)=10.^(P(1)-2*stats.robust_s(1));
        %y_err(i,2)=max(y_err1(i,2), y_err1(i,2));
    else
        y(i,1)=NaN;
        y(i,2)=NaN;
        y_err(i,1)=NaN;


    end

end


%Martine Collaud Coen
%MeteoSwiss
% 21.8.2007 / 7.3.2024