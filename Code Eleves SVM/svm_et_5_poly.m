% Séparateurs à Vaste Marge
% Cas linéaire, non linéaire, marges "souples"
% structure du code: c'est juste une proposition, à vous de voir
% fonction de décision <w,x> + b 
clear all
close all




color = ['r' 'b' 'g' 'c' 'm'];
 load 'data3' X lab     % chargement des données 
sigma = [  1  ];
cte= 10;
N_supports=zeros(1,length(sigma));
kernel = @(x,y,sigma,cte)((x'*y+ cte)^sigma); 
% bornes pour le dessin 2D
xmin=min(X(1,:));
ymin=min(X(2,:));
xmax=max(X(1,:));
ymax=max(X(2,:));

na=length(X); % nombre de points (d'apprentissage)
C=[ 10000000 ];
for indicesigma = 1 : length(sigma)
for k=1:length(C)

% dessin des points dans R^2

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



% assemblage matrice de la forme quadratique du problème dual
A=zeros(na); %initialisation

% A vous de construire A
for i= 1:na
    for j = 1:na
        A(i,j)= lab(i)*lab(j)*kernel(X(:,i),X(:,j),sigma(indicesigma),cte);
        %A(i,j)= lab(i)*lab(j)*kernel(X(:,i),X(:,j));
    end
end



% gradient à pas constant pour problème dual
alph0=0.5*ones(na,1); % point de départ (0, ou autre): réglable
pasgrad=5e-3;         % pas du gradient : parametre réglable
u=ones(na,1);         % vecteur de 1
crit_arret=1;         % initialisation critere d'arret
npas=0;               % comptage du nombre de pas
npas_max=100000;      % garde-fou convergence (si ça ne converge pas...)
epsi=1e-5;            % seuil convergence




while and(crit_arret>epsi,npas<npas_max)

    % ici la boucle de gradient à pas constant
    
    dH = A*alph0 - u;
    alph = alph0 - pasgrad*dH;    
    alph = alph - ((lab'*alph)*lab)/(norm(lab)^2); %proj sur l'hyperplan 
    alph = max(alph,zeros(na,1)); %proj sur R+
    crit_arret= norm(alph - alph0);
    npas = npas + 1;
    alph0 = alph;                        
end


% Ce n'est pas fini! Il faut construire la fonction de decision <w,x> + b
% tracer ses isovaleurs 0 (la droite de partage) et +1 et -1 : les marges
% donc calculer b d'abord, et <w,x> pour tout x 

% comment calculer b? idée: en se servant des "points supports" x_s pour
% lesquels <w,x_s> +b = + ou - 1
% comment trouver les points support? Qu'est ce qui le caractérise?? 

Points_support=find(alph>epsi)' % points supports = contraintes actives...alph>0
                                % aux erreurs numériques près

% Calcul de b (on le calcule pour chaque pt support, puis on moyenne)
% on devrait vérifier que chaque pt support fournit le même b (aux erreurs
% numériques près, c'est pouquoi on moyenne)
nb_sup=0;

Xs = X(:,Points_support);
labsup= lab(Points_support);

%calcul de w et tracé de la droite et de ses marges

for i = 1 : length(Points_support)  % i décrit le tableau des indices des points supports
  %calcul de b pour chaque point support
  S=0;
  for j=1 : na
   S =  S + alph(j)*lab(j)*kernel(X(:,j), Xs(:,i),sigma(indicesigma),cte);
 
  end
  w0(i) = labsup(i) - S ;
end
if diff(w0)>0.1, error('All b should be the same'), end
b = mean(w0);




% une méthode, qui parait lourde mais qui est la bonne dans le cas de 
% séparation non linéaire, est de construire la fonction de décision
% sur une grille et d'en tracer les trois isovaleurs avec une fonction
% Matlab



xp=xmin:0.2:xmax;   % création d'une grille 
yp=ymin:0.2:ymax;
npx=length(xp);
npy=length(yp);
for i=1:npx
    for j=1:npy
        % calcul de <w,x> + b sur une grille
        x= [xp(i); yp(j)];
        somme =0;
        for l = 1 : na 
            somme = somme + alph(l)*lab(l)*kernel(X(:,l),x,sigma(indicesigma),cte);
            
        end
        V(i,j)=   somme +b;
    end
end
% function s = som(x,p,l,X,na)
% kernel = @(x,y)exp((-(norm(x-y)^2)));
%     s=0;
%     for i=1 : na
%     s = s + p(i)*l(i)*kernel(X(:,i),x);
%     end
% end
hold on
contour(xp,yp,V',[0 0],'linewidth',2,'color','r') % dessin des contours
contour(xp,yp,V',[-1 1],'linewidth',2,'color','r', 'linestyle',':')% dessin des contours
axis([xmin xmax ymin ymax])
title(' Suport Vector Machine')
grid

%N_supports(indicesigma) = length(Points_support);
%pause

end
end

%plot(sigma,N_supports);
%xlabel('sigma') % x-axis label
%ylabel('Nbr points de supports') % y-axis label
