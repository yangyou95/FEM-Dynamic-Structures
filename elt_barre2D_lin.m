function [res,res2] = elt_l2_lin(code,X,prel,Ue)
%
% element de barre lineaire a deux noeuds en dimension 2
% 

if ( code==0 ),  % nombre de ddl/noeuds
    
    res = 2;  
    return
    
end

% proprietes elementaires
E  = prel(1);
nu = prel(2);
rho= prel(3);
S0 = prel(4);

% longueur de l'element
x1 = X(1,1) ; y1 = X(1,2) ;
x2 = X(2,1) ; y2 = X(2,2) ;

Lx0 = (x2-x1); Ly0 = (y2-y1);
L0 = sqrt( Lx0^2 + Ly0^2 );   

B = [ -1/L0 0 1/L0 0 ];  % relation deplacements-deformation local

Q=[ Lx0  -Ly0 
    Ly0   Lx0 ]/L0; % matrice de changement de repere

vT = blkdiag(Q,Q)'; % 

BT = B*vT;

K = (E*S0*L0)*BT'*BT ; % Matrice de rigidite

if ( code==1 ),  % Matrice de rigidite tangente (Kt)
        
    res = K ;
    
elseif ( code==2 ), % matrice de masse elementaire
    
    %res = rho*S0*L0/6*[ 2 0 1 0  
    %                    0 2 0 1  
    %                    1 0 2 0 
    %                    0 1 0 2 ] ;
                    
    res = .5*rho*S0*L0*eye(4) ; % masses concentrees
    
elseif ( code==3 ), % forces internes et contrainte normale
        
    res = K*Ue;  % forces internes dans le repere global
    
    res2 = E*BT*Ue; % contraintes 
        
end



