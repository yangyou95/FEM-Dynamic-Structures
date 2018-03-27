function [res,res2] = elt_l2(code,X,prel,Ue)
%
% element de barre non-lineaire a deux noeuds en dimension 2
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

% longueur initiale et courante
x1 = X(1,1) ; y1 = X(1,2) ;
x2 = X(2,1) ; y2 = X(2,2) ;

Lx0 = (x2-x1); Ly0 = (y2-y1);
L0 = sqrt( Lx0^2 + Ly0^2 );


Xe = X + [ Ue(1) Ue(2); Ue(3) Ue(4) ]; Xe = Xe'; Xe = Xe(:)';

Lx = Xe(3)-Xe(1); Ly = Xe(4)-Xe(2);

L = sqrt( Lx^2 + Ly^2 ) ;  % Longueur de l'element 

N = E*S0*(L-L0)/L0;  % effort normal dans l'element


if ( code==1 ),  % Matrice de rigidite tangente (Kt)
    
    B = [ -1/L 0 1/L 0 ];  % relation deplacements-deformation local

    Q=[ Lx  -Ly 
        Ly   Lx ]/L; % matrice de changement de repere

    vT = blkdiag(Q,Q)'; % 

    BT = B*vT;

    K = (E*S0*L)*BT'*BT ; % Matrice de rigidite
 
    k = eye(2)*N/L ;  % et stress stiffening 
    
    res = K + [k -k; -k k] ;
        
elseif ( code==2 ), % matrice de masse elementaire
    
    res = rho*S0*L0/6*[ 2 0 1 0  
                        0 2 0 1  
                        1 0 2 0 
                        0 1 0 2 ] ;
    
elseif ( code==3 ), % forces internes, contrainte normale, et energie interne
    
    
    res = [  -Lx
             -Ly 
              Lx 
              Ly ] * N/L;
            
    res2 = N/S0 ;
    
  %  res3 = .5*(L-L0)*N;  % energie interne
    
end



