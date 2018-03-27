

% donnees du probleme
E = 2e11 ;
nu = 0.3 ;
rho=7800;
S = 4e-5 ;

% effort statique :
Fya = -1000 ;   % amplitude de l'effort
Fyt = @(t) Fya; % fonction temporelle de l'effort

% effort impulsionnel :
% T1 = 0.0012;
% Fya = 1000 ;   % amplitude de l'effort
% Fyt = @(t) Fya.*sin(pi/T1*t).*(t<T1); % fonction temporelle de l'effort


% % effort dynamique sinus :
% Fym = 0 ;   % moyenne de l'effort
% Fya = 1;   % amplitude de l'effort
% fw = 55.8;      % frequence d'excitation
% Fyt = @(t) Fym + Fya*sin(2*pi*fw*t); % fonction temporelle de l'effort


% % effort dynamique sinus balaye ou vobule :
% Fym = 0 ;   % moyenne de l'effort
% Fya = 8;   % amplitude de l'effort
% fw = @(t)20+55.8*t/4;      % frequence d'excitation
% Fyt = @(t) Fym + Fya*sin(2*pi*fw(t)*t); % fonction temporelle de l'effort


element = 'elt_barre2D_lin'; 
element = 'elt_barre2D_nl'; 

npas = 5000;          % nombre de pas de temps
             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coordonnees des noeuds du maillage
vcor = [  0    0
          .5   .01
          1    0 ];


% table des connectivites
kconec = [1 2
          2 3] ;
      

% sollicitations concentrees :
%-----------------------------
% le tableau est organise par numero de noeuds, 
% numero de ddl et valeur imposee
% [noeud ddl norme de Fimp]

sollicitation = { 2 2 Fyt };

                
% conditions aux limites en deplacements :
limites = [ 1 1 0 ; 1 2 0
            3 1 0 ; 3 2 0 ] ;

% conditions initiales
U0 = 0;
Ud0 = 0;

prel = [ E, nu, rho, S ];
clear E nu rho S F

