% Groupe: Najlaa Srifi// Hajar EL fahfouhi // Léa Bagheri

% Séparateurs à Vaste Marge pour régression
% Code pour élèves à développer

clear all
close all

global ikernel  % type de noyau choisi. 1: produit scalaire dans Rp
                                                         % 2: noyau gaussien
global dimX
 
%données X et Y , suivant cas: courbe ou nappe 
%load 'nappe' X Y
load 'courbe' X Y
 
npoints=size(X,2);
dimX=size(X,1);
 
% Les parametres des SVR 
epsilon=0.1; % demi largeur bandea %0.5
ikernel=2;              % choix du type de noyau
marge_souple= true ;  % choix marge souple ou pas
C_souple= 50;            % coef souplesse de la marge souple 
 

 
% dessin des points suivant dimension
if dimX==1
    plot(X,Y, 'o')
    axis equal
    grid on
    hold on
else
    plot3(X(1,:),X(2,:),Y, 'o')
    axis equal
    hold on
    pause
end

% Coeur de la méthode SVR, à programmer
na=2*npoints; % nombre de multiplicateurs alpha
u1=epsilon*ones(na,1);  % vecteurs colonne
u2=[Y'; -Y'];
u3=[ones(npoints,1); -ones(npoints,1)];

% assemblage matrice de la forme quadratique du problème dual A
B=zeros(npoints);
for i=1:npoints
    for j=1:npoints
        B(i,j)=kernel(X(:,i),X(:,j));
    end
end
A=[B -B; -B B]; 

% gradient à pas constant pour problème dual
alph0=2.5*abs(rand(na,1)); % point de départ (0, ou autre): réglable
pasgrad=5e-3;         % pas du gradient : parametre réglable
crit_arret=1;         % initialisation critere d'arret
npas=1;               % comptage du nombre de pas
npas_max=100000;      % garde-fou convergence
epsi=1e-6;                   % seuil convergence
H_alpha = - 0.5 * alph0'* A * alph0 -u1' * alph0 + u2'*alph0;
grad_H = -A * alph0 - u1 + u2;
alpha_k = alph0;
alpha_k_plus_1 = 0;
vit_conv = zeros(1,npas_max);

while and(crit_arret>epsi,npas<npas_max) % boucle gradient
    npas=npas+1;
    alpha_k_plus_1 = alpha_k + pasgrad * grad_H;% gradient pas constant + heuristique approx contrainte
    alpha_k_plus_1 = max(alpha_k_plus_1,0); 
    alpha_k_plus_1 = alpha_k_plus_1 -(u3' * alpha_k_plus_1) * u3 / (norm(u3)^2);
    
   
    %alpha_k_plus_1 = max(alpha_k_plus_1,0);                                % ne pas oublier la projection  projection d' ABORD sur l'hyperplan
    
    if marge_souple == true
        alpha_k_plus_1 = min(alpha_k_plus_1,C_souple);
    end
    
    
    H_alpha = -0.5 * alpha_k_plus_1' * A * alpha_k_plus_1 - u1' * alpha_k_plus_1 +u2'*alpha_k_plus_1 ;
    grad_H = -A * alpha_k_plus_1 - u1 + u2;
    crit_arret = norm(alpha_k_plus_1 - alpha_k);
    vit_conv(npas)=crit_arret; % on stocke pour visualiser le comportement
    alpha_k = alpha_k_plus_1; 
 
end

fprintf(' nb iterations: %d critère: %e \r \r ',npas,crit_arret)
alph = alpha_k_plus_1;
Points_support=find(alph>epsi)';

pause
% on ne se sert plus directement de w, mais de son produit scalaire 
% avec un vecteur, defini par un noyau. Ce produit scalaire dépend des coefs alpha
% et des points d'apprentissage : fonction prodwr 

% Calcul de b (on le calcule pour chaque pt support, puis on moyenne)
% on devrait vérifier que chaque pt support fournit le même b (aux erreurs
% numériques près, c'est pouquoi on prend la moyenne)


% Calcul de b 
nb=0;
for i=Points_support
    if abs(alph(i)-C_souple) > epsi
        nb=nb+1;
        
        if i>npoints
         k=i-npoints;
         p_scal=prodwr(X,alph,X(:,k)); 
         moy_b(nb)=-p_scal+Y(k)+epsilon;
        else
         k=i;
         p_scal=prodwr(X,alph,X(:,k));
         moy_b(nb)=-p_scal+Y(k)-epsilon;
        end
    end
end

b=mean(moy_b);
save 'Solution' alph b   % on sauve dans solution, fichier que va utiliser Dessin Interp pour dessiner



fprintf(' nb iterations: %d critère: %e  \r',npas,crit_arret)

DessinInterp   % on lance le dessin
