% Groupe: Najlaa Srifi// Hajar EL fahfouhi // Léa Bagheri

% Séparateurs à Vaste Marge pour la Classification
% Cas linéaire, non linéaire, marges "souples"
clear all
close all
global ikernel  % type de noyau choisi.    
load 'data2' X lab      % chargement des données 
ikernel = 2;            % choix du type de noyau: 1 produit scalaire
                                                  %2 noyau gaussien 
                                                  % 4 noyau laplacien
                                                  % 5 noyau polynomial
marge_souple = true;   % choix marge souple ou pas
C_souple = 10;          % coef souplesse de la marge souple
ajout = true;          % Ajouter un nouveau point

% bornes pour le dessin 2D

xmin = min(X(1,:));
ymin = min(X(2,:));
xmax = max(X(1,:));
ymax = max(X(2,:));

na = length(X); % nombre de points (d'apprentissage)

% dessin des points dans R^2
subplot(1,2,1)
for i=1:na
    if lab(i)==1
        plot(X(1,i),X(2,i),'o','linewidth',2)
        hold on
    else
        plot(X(1,i),X(2,i),'x','linewidth',2,'markersize',12)
        hold on
    end
end
axis([xmin xmax ymin ymax])
grid
axis equal
hold on
pause
% Coeur du Programme - Méthode d'UZAWA
% assemblage matrice de la forme quadratique du problème dual

A = zeros(na);
for i = 1:na
    for j=1:na
        A(i,j) = lab(i)*lab(j)* kernel(X(:,i),X(:,j));
    end
end

% Gradient projeté à pas constant pour problème dual

alpha0 = 0.5*ones(na,1); % point de départ (0, ou autre): réglable
pasgrad = 5e-3;          % pas du gradient : parametre réglable
u = ones(na,1);          % vecteur de 1
crit_arret = 1;          % initialisation critere d'arret
npas = 1;                % comptage du nombre de pas
npas_max = 100000;       % garde-fou convergence : on arrete si le nombre d'iterations est trop grand
epsi = 1e-5;             % seuil convergence
% Boucle du gradient projeté à pas constant + UZAWA
H_alpha = - 0.5 * alpha0'* A * alpha0 + u' * alpha0 ;
grad_H = -A * alpha0 + u;
alpha_k = alpha0;
alpha_k_plus_1 = 0;
vit_conv = zeros(1,npas_max);
while crit_arret >= epsi && npas <= npas_max
    
    alpha_k_plus_1 = alpha_k + pasgrad * grad_H;
    alpha_k_plus_1 = alpha_k_plus_1 - (lab' * alpha_k_plus_1) * lab / (norm(lab)^2);
    
    % En cas de contrainte inégalité, il faut projetter notre alpha sur le cadran de (R+)^p
    alpha_k_plus_1 = max(alpha_k_plus_1,0);
    % Dans le cas de marge souple, il faut faire attention dans le calcul
    % de b et prendre alpha et C différents
    if marge_souple == true
        alpha_k_plus_1 = min(alpha_k_plus_1,C_souple);
    end
    
    H_alpha = -0.5 * alpha_k_plus_1' * A * alpha_k_plus_1 + u' * alpha_k_plus_1;
    grad_H = -A * alpha_k_plus_1 + u;
    crit_arret = norm(alpha_k_plus_1 - alpha_k);
    vit_conv(npas) = crit_arret;
    alpha_k = alpha_k_plus_1; 
    npas = npas + 1;
end 
% Recherche des points supports
epsi = 1e-5;
alph = alpha_k_plus_1;
Points_support = find(alph > epsi); 

% on ne se sert plus directement de w, mais de son produit scalaire 
% avec un vecteur, defini par un noyeau : fonction prodw

% Calcul de b (on le calcule pour chaque pt support, puis on moyenne)
% chaque pt support devrait fournit le même b (mais aux erreurs
% numériques près, c'est pouquoi on moyenne)
moy_b = zeros(length(Points_support),1);

for i=1:length(Points_support)
    indice = Points_support(i);
    moy_b(i) = 1/lab(indice) - prodw(X,lab,alph,X(:,indice)); 
end

b = mean(moy_b);      %  on moyenne les b des points supports
grid
%calcul et tracé des isovaleurs

xp = xmin:0.2:xmax;   % création d'une grille pour les besoins de contour
yp = ymin:0.2:ymax;
npx = length(xp);
npy = length(yp);

for i=1:npx
    for j=1:npy
        ps = prodw(X,lab,alph,[xp(i),yp(j)]'); % calcul de <w,x> + b sur une grille,on n'a pas besoin explicitement de w 
        V(i,j) = ps + b;                           
    end
end

hold on
contour(xp,yp,V',[-1 0 1],'linewidth',2,'color','r')
axis([xmin xmax ymin ymax])
title(' Support vecteur Machine avec marge souple ')
grid
pause
subplot(1,2,2)
plot(log(vit_conv),'linewidth',2)
xlabel(' nombre itérations')
ylabel(' log (ecart)')
title('comportement convergence')
grid
% Ajouter un nouveau point

if ajout == true
    X_n = zeros(2,length(X));
    for i=1:length(X_n)
        X_n(:,i) = X(:,i);
    end
    X_n(:,length(X_n))= ginput(1);
    pause
    lab_n= zeros(length(X_n),1);
    for i=1:length(X)
        lab_n(i) = lab(i); 
    end
    rayon = 10;
    compteur_L1 = 0;
    compteur_L2 = 0;
    %  on trace un cercle de centre ce nouveau point, 
    % et on compte le nombre de points de label 1
    % et de label -1 qui l'entourent. S'il y a plus de points de label 1
    % dans ce cercle, on attribut au nouveau point le label 1 sion on lui
    % attribut -1 
    for i=1:na
        if (X(1,i)- X_n(1,length(X_n)))^2 + (X(2,i)-X_n(2,length(X_n)))^2 <= rayon
           if  lab(i) == 1
               compteur_L1 = compteur_L1 + 1;
           else
               compteur_L2 = compteur_L2 + 1;
           end
        end
    end

    if compteur_L1 > compteur_L2
        lab_n(length(X_n)) = 1;
        plot(X_n(1,length(X_n)),X_n(2,length(X_n)),'o','linewidth',2)
    else
        lab_n(length(X_n)) = -1;
        plot(X_n(1,length(X_n)),X_n(2,length(X_n)),'x','linewidth',2,'markersize',12)
    end
end

na = length(X_n);
A = zeros(na);
for i = 1:na
    for j=1:na
        A(i,j) = lab_n(i)*lab_n(j)* kernel(X_n(:,i),X_n(:,j));
    end
end
alpha0 = 0.5*ones(na,1); 
pasgrad = 5e-3;          
u = ones(na,1);         
crit_arret = 1;          
npas = 1;               
npas_max = 100000;       
epsi = 1e-5;            

H_alpha = - 0.5 * alpha0'* A * alpha0 + u' * alpha0 ;
grad_H = -A * alpha0 + u;
alpha_k = alpha0;
alpha_k_plus_1 = 0;
vit_conv = zeros(1,npas_max);

while crit_arret >= epsi && npas <= npas_max
    
    alpha_k_plus_1 = alpha_k + pasgrad * grad_H;
    alpha_k_plus_1 = alpha_k_plus_1 - (lab_n' * alpha_k_plus_1) * lab_n / (norm(lab_n)^2);
    alpha_k_plus_1= max(alpha_k_plus_1,0);  
    if marge_souple == true
        alpha_k_plus_1 = min(alpha_k_plus_1,C_souple);    
    end
    
    H_alpha = -0.5 * alpha_k_plus_1' * A * alpha_k_plus_1 + u' * alpha_k_plus_1;
    grad_H = -A * alpha_k_plus_1 + u;
    crit_arret = norm(alpha_k_plus_1 - alpha_k);
    vit_conv(npas) = crit_arret;
    alpha_k = alpha_k_plus_1; 
    npas = npas + 1;
end

epsi = 1e-5;
alph = alpha_k_plus_1;
Points_support = find(alph > epsi);

for i=1:length(Points_support)
    indice = Points_support(i);
    moy_b(i) = 1/lab_n(indice) - prodw(X_n,lab_n,alph,X_n(:,indice)); 
end

b = mean(moy_b);  
grid

for i=1:npx
    for j=1:npy
        ps = prodw(X_n,lab_n,alph,[xp(i),yp(j)]');
        V_new(i,j) = ps + b;                           
    end
end

hold on
contour(xp,yp,V_new',[-1 0 1],'linewidth',2,'color','b')
axis([xmin xmax ymin ymax])
title('Suport Vector Machine')
grid
hold on 
subplot(1,2,2)
plot(log(vit_conv),'linewidth',2)
xlabel(' nombre itérations')
ylabel(' log (ecart)')
title('comportement convergence')
grid


   
                      




