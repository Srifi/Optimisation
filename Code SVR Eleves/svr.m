% Séparateurs à Vaste Marge pour régression
% Code pour élèves à développer

clear all
close all

global ikernel  % type de noyau choisi. 1: produit scalaire dans Rp
                                                         % 2: noyau gaussien
 global dimX
 
 %données X et Y , suivant cas: courbe ou nappe (changer le commentaire)
 load 'nappe' X Y
 %load 'courbe' X Y
 
 npoints=size(X,2);
 dimX=size(X,1)
 
% Les parametres des SVR 
epsilon=0.5; % demi largeur bande
 ikernel=2;                   % choix du type de noyau
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
for i=1:npoints
    for j=1:npoints
     % ............
    end
end

% A=.........;  % matrice 

% gradient à pas constant pour problème dual
alph0=2.5*abs(rand(na,1)); % point de départ (0, ou autre): réglable
pasgrad=5e-3;         % pas du gradient : parametre réglable
crit_arret=1;         % initialisation critere d'arret
npas=0;               % comptage du nombre de pas
npas_max=100000;      % garde-fou convergence
epsi=1e-6;                   % seuil convergence


while and(crit_arret>epsi,npas<npas_max) % boucle gradient
    npas=npas+1;
    alph=alph0+ pasgrad*.......; % gradient pas constant + heuristique approx contraintes
                                             % ne pas oublier la projection  projection d' ABORD sur l'hyperplan
    
    if marge_souple
        
    end
    
    
    crit_arret=norm(alph-alph0);
    vit_conv(npas)=crit_arret; % on stocke pour visualiser le comportement
    alph0=alph;
 
end

fprintf(' nb iterations: %d critère: %e \r \r ',npas,crit_arret)
Points_support=find(alph>epsi)';


fprintf('PRESS RETURN TO CONTINUE')
pause
% on ne se sert plus directement de w, mais de son produit scalaire 
% avec un vecteur, defini par un noyau. Ce produit scalaire dépend des coefs alpha
% et des points d'apprentissage : fonction prodwr 

% Calcul de b (on le calcule pour chaque pt support, puis on moyenne)
% on devrait vérifier que chaque pt support fournit le même b (aux erreurs
% numériques près, c'est pouquoi on prend la moyenne)

% Calcul de b 
nb_sup=0;
for i=Points_support
  ..........
end
b=mean(moy_b);   %  on moyenne les b des points supports

save 'Solution' alph b   % on sauve dans solution, fichier que va utiliser Dessin Interp pour dessiner



fprintf(' nb iterations: %d critère: %e  \r',npas,crit_arret)

DessinInterp   % on lance le dessin
