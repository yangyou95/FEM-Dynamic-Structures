%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%      Analyse d'un probleme dynamique par la methode des EF              %
%      Integration directe explicite
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'\n\n\nDebut du programme d''analyse par la methode des EF \n')
clear all ;  clc; close all; format short g


% chargement du MEF
%------------------
fprintf(1,'lecture des parametres du MEF \n')
%prmaill = uigetfile('*.m','Mesh File : ');  % on demande le maillage a l'utilisateur
prmaill = 'maillage.m';
eval(prmaill(1:(findstr(prmaill,'.'))-1));


% affichage du maillage
figure
patch('Vertices',vcor,'Faces',kconec,'FaceColor','w'); hold on
axis('equal')


% initialisation des variables et matrices necessaires
%-----------------------------------------------------
ndln = feval( element, 0) ; % nombre de ddl/noeuds
nnt = size(vcor,1);     % nb. de noeuds du pb. EF
neq = nnt*ndln ;        % nb. d'equations
nel = size(kconec,1);   % nb. d'elements

U   = zeros(neq,2);    % deplacement
Ud  = zeros(neq,2);    % vitesse
Udd = zeros(neq,2);    % acceleration

Kg  = sparse(neq,neq) ; % matrice raideur globale initialisee a zero
Mg  = sparse(neq,neq) ; % matrice masse globale initialisee a zero
Fg  = sparse(neq,1) ;   % vecteur force global initialise a zero

time = zeros(1,npas);
vsol = zeros(neq,npas);
sigma = zeros(nel,npas);
vfint = zeros(1,npas);
vfext = zeros(1,npas);


% evaluation des matrices globales
%-----------------------------------------
fprintf('assemblage des matrices globales \n')

for (ielt = 1:nel), % boucle sur les elements
    
   kconece = nonzeros(kconec(ielt,:))';                % connectivite de l'element

   ddle = repmat(kconece*ndln,ndln,1) - repmat(ndln-1:-1:0,length(kconece),1)' ;  % localisation des raideurs
   Ke = feval( element, 1, vcor(kconece,:), prel, 0*ddle) ; % evaluation de la matrice raideur elementaire
   Me = feval( element, 2, vcor(kconece,:), prel, 0*ddle) ; 
   Kg(ddle(:),ddle(:)) = Kg(ddle(:),ddle(:)) + Ke ;    % assemblage dans la matrice globale
   Mg(ddle(:),ddle(:)) = Mg(ddle(:),ddle(:)) + Me ;
   
end 



% evaluation du vecteur force (C.L. sollicitation)
%-------------------------------------------------
nsol = 0;
if exist('sollicitation'), 
    fprintf('assemblage des efforts \n')
    nsol = size(sollicitation,1);
    idxf = zeros(1,nsol);
    %Fb = zeros(1,size(sollicitation,1));
    for kf = 1:nsol,

        noe = sollicitation{kf,1};
        plot(vcor(noe,1),vcor(noe,2),'o')

        idxf(kf) = (noe-1)*ndln+sollicitation{kf,2} ;
        %Fb(kf) = sollicitation{kf,3} ;
    end
else
    idxf = [];
end


Cg = 5e-4*abs(Kg);
%Cg = diag(diag(Cg));

R0 = zeros(neq,1);

% formation de la matrice des contraintes cinematiques
%-----------------------------------------------------
fprintf('prise en compte des contraintes cinematiques \n')

ddle = (limites(:,1)-1)*ndln+limites(:,2) ;
kdle = 1:neq;
kdle(ddle) = [];
Mii = Mg(kdle,kdle);
Kii = Kg(kdle,kdle);
Cii = Cg(kdle,kdle);


% Pour le bilan d'energie du schema numerique

T = zeros(1,npas);
V = zeros(1,npas);
dF = zeros(1,npas);
dC = zeros(1,npas);

% resolution du probleme 
%-------------------------------
fprintf('resolution directe (analyse) du pb. EF \n')

dt  = 2/sqrt(eigs(Kii,Mii,1))/pi;
dt2 = dt*dt;

Fg = zeros(neq,1);    % effort initial
Fb = sollicitation{kf,3};
Fg(idxf(kf)) = Fb(0);
F0 = Fg(kdle);
    
U(:,1) = U0;    % d?eplacement
Ud(:,1)= Ud0;    % vitesse
Udd(kdle) = Mii\(F0 - Kii*U(kdle,1) - Cii*Ud(kdle,1)); % acce?le?ration

R0 = Mg*Udd(:,1) + Cg*Ud(:,1) + Kg*U(:,1) - Fg; % bilan d'efforts 
    Fint0 = Kg*U(:,1); % Fint(0)

T(1) = Ud(:,1)'*(Mg*Ud(:,1))*.5;
V(1) = U(:,1)'*(Kg*U(:,1)*.5) ;
    vfext(1) = Fb(0);


for kt=2:npas,
    
    time(kt) = time(kt-1) + dt;

    for kf=1:nsol,
        
        Fb = sollicitation{kf,3};
        Fg(idxf(kf)) = Fb(time(kt));
        
    end
    F = Fg(kdle);

    
    U(:,2) = U(:,1) + Ud(:,1)*dt + .5*Udd(:,1)*dt2;
    
    Fint = zeros(neq,1) ;
    for (ielt = 1:nel), % boucle sur les elements
        
        kconece = nonzeros(kconec(ielt,:))';                % connectivite de l'element
        ddle    = repmat(kconece*ndln,ndln,1) - repmat(ndln-1:-1:0,length(kconece),1)' ;  % localisation des raideurs
        [Fe,Se] = feval( element, 3, vcor(kconece,:), prel, U(ddle,2)) ; % evaluation de la matrice raideur elementaire
        Fint(ddle(:))  = Fint(ddle(:)) + Fe ;
        sigma(ielt,kt) = Se;
    end
    
    Udd(kdle,2) = (Mii + .5*Cii*dt)\(F - Fint(kdle) - Cii*(Ud(kdle,1)+.5*Udd(kdle,1)*dt));
    Ud(:,2)  = Ud(:,1) + .5*(Udd(:,1)+Udd(:,2))*dt;
        
    R = Mg*Udd(:,2) + Cg*Ud(:,2) + Fint - Fg;  % efforts de reactions
    
    T(kt) = .5*Ud(:,2)'*(Mg*Ud(:,2));
    V(kt) = V(kt-1) + .5*(U(:,2)-U(:,1))'*(Fint+Fint0) ;
    
    dF(kt) = dF(kt-1) + .5*(U(:,2)-U(:,1))'*(R+R0) ;
    dF(kt) = dF(kt) + .5*(U(kdle,2)-U(kdle,1))'*(F+F0) ;
        
    dC(kt) = dC(kt-1) + (Ud(:,2)+Ud(:,1))'*(Cg*(Ud(:,2)+Ud(:,1)))*dt*.25;
    
    U(:,1)   = U(:,2);
    Ud(:,1)  = Ud(:,2);
    Udd(:,1) = Udd(:,2);
    F0 = F;
    R0 = R;
    Fint0 = Fint;
    
    vsol(:,kt) = U(:,1);
    vfint(kt) = Fint(kdle(end));
    vfext(kt) = Fb(time(kt));
    
end

T = T - T(1);
V = V - V(1);




% affichage de la deformee  
%-------------------------
fprintf('affichage des resultats et de la deformee \n')

fprintf('  Deplacements maximaux : %g m\n',max(max(abs(vsol))));
fprintf('  Contraintes maximales : %g MPa\n',max(max(abs(sigma)))./1e6);

% Udef = reshape(U(:,1)./max(abs(U(:,1)))*.05*max(max(vcor)),ndln,nnt)' ; % deformee avec amplification
% patch('Vertices',vcor(:,1:2)+Udef(:,1:2),'Faces',kconec,'FaceColor','w'); hold on

amp = .1*(max(vcor(:))-min(vcor(:)))./max(abs(vsol(:))); %amp = 1;
kbound = unique(limites(:,1));
Cdata = colormap;
cdata = linspace(-max(abs(sigma(:))),max(abs(sigma(:))),size(Cdata,1))';

figure(1); clf
Vef = vcor+amp*reshape(vsol(:,1),ndln,nnt)' ; % deformee avec amplification
h2 = plot(vcor(kbound,1),vcor(kbound,2)-0.03,'^','MarkerFaceColor','g','MarkerSize',14); hold on
h1 = patch('Vertices',Vef,'Faces',kconec,'FaceColor','w','LineWidth',3,'EdgeColor','b'); hold on
h3 = plot(vcor(:,1),vcor(:,2),'o','MarkerFaceColor','w','MarkerSize',7); hold on
axis equal; set(gca,'Xlim',[0,1],'Ylim',[-.35,.35],'Ytick',[]); 
caxis([-max(abs(sigma(:))),max(abs(sigma(:)))])
% h = colorbar; axes(h); title('\sigma [Pa]')
drawnow

for kt2 = 2:10:npas
    color = [ interp1(cdata,Cdata(:,1),sigma(1,kt2)), ...
              interp1(cdata,Cdata(:,2),sigma(1,kt2)), ...
              interp1(cdata,Cdata(:,3),sigma(1,kt2)) ];

    Vef = vcor + amp*reshape(vsol(:,kt2),ndln,nnt)';
    set(h1,'Vertices',Vef,'EdgeColor',color);
    set(h3,'XData',Vef(:,1),'YData',Vef(:,2));
    drawnow; 
end


figure
subplot(2,1,1)
plot(time,vsol(4,:))
xlabel('Temps [s]')
ylabel('Deplacement [m]')
title([prmaill,' - DDL #4'],'Interpreter','none')

subplot(2,1,2)
plot(time,sigma(1,:)./1e6)
xlabel('Temps [s]')
ylabel('Contrainte [MPa]')

figure
subplot(2,1,1)
plot(time,vcor(2,2)+vsol(4,:))
xlabel('Temps [s]')
ylabel('Position [m]')
title([prmaill,' - DDL #4'],'Interpreter','none')

subplot(2,1,2)
plot(time,sigma(1,:)./1e6)
xlabel('Temps [s]')
ylabel('Contrainte [MPa]')

figure(4); 
plot(vsol(4,:),vfint,vsol(4,end),vfint(end),'o')
grid on
ylabel('$F_{int}$ [N]','Interpreter','latex')
xlabel('D\''eplacement [m]','Interpreter','latex')
title([prmaill,' - DDL #4'],'Interpreter','none')


figure(5); hold on; clf
plot(time,dF,'y','linewidth',3);  hold on
plot(time,T,time,V,time,dC,time,T+V+dC);
grid on
xlabel('Temps [s]')
ylabel('Energies [Nm]')
legend('dF','T','V','dC','T+V+dC','location','Best')


% evaluation des coefficients de la transformation de Fourier discretisee
% On veut H verifiant Y = H.X
%------------------------------------------------------------------------
k = (1:npas)-1; % partie du signal a transformer
N  = length(k);
T = max(time(N));
% w = (0.21557895 - 0.41663158*cos(2*pi*time/T) + 0.277263158*cos(4*pi*time/T) - ...
%     0.083578947*cos(6*pi*time/T) + 0.006947368*cos(8*pi*time/T))/0.21557895;
w = 1;

xk = w.*vfext(1,k+1);      % extraction des efforts temporels
%yk = w.*vsol(4,k+1);      % extraction des deplacements temporels
yk = w.*sigma(1,k+1);      % extraction des contraintes temporelles

kfin = fix(200*dt*N);  % on coupe a 5000 Hz
X = zeros(1,kfin); Y = zeros(1,kfin); % initialisation des tableaux de coefficients 
f = (0:kfin)/(N*dt); % initialisation du tableau de frequences f
for p=0:kfin,
    X(p+1) = xk*exp(-1j*2*pi*p/N*k')*dt;
    Y(p+1) = yk*exp(-1j*2*pi*p/N*k')*dt;
end
Sxy = X.*conj(Y); % interspectre entre entree et sortie
Syx = Y.*conj(X); % interspectre entre sortie et entree
Sxx = X.*conj(X); % autospectre de l'entree
Syy = Y.*conj(Y); % autospectre de la sortie
Hyx = Syx./Sxx; % FRF (estimateur 1)
%Hxy = Syy./Sxy; % FRF (estimateur 2)

figure(20) % affichage des coefficients X
title('Transform\''ee de Fourier discretis\''ee de x(t)','Interpreter','latex')
semilogy(f,abs(X));  grid on; 
xlabel('f [Hz]'); ylabel('|X| [N/Hz]')
 
% figure(21) % affichage des coefficients X
% title('Transform\''ee de Fourier discretis\''ee de y(t)','Interpreter','latex')
% semilogy(f,abs(Y));  grid on; 
% xlabel('f [Hz]'); ylabel('|Y| [m/Hz]')

figure(22) % affichage des coefficients Y
title('Transform\''ee de Fourier discretis\''ee de h(t)','Interpreter','latex')
semilogy(f,abs(Hyx));  hold on; 
%semilogy(f,abs(Hxy));  grid on; 
xlabel('f [Hz]'); ylabel('|H| [m/N]')


% fin
%----
fprintf('fin de l''analyse du probleme par EF \n\n')

