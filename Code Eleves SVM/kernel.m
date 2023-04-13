function kern=kernel(u,v)
%  fonction renvoyant la valeur du noyau K(u,v)
% u et v sont deux vecteurs (colonnes) de Rn
global ikernel % passé en global pour alleger

if ikernel==1
    kern=u'*v;
elseif ikernel==2  % noyau Gaussien
    kern=exp(-0.1*norm(u-v)^2);
elseif ikernel==4
    kern=exp(-0.44*norm(u-v));%noyau laplacien
elseif ikernel==5
    kern=(u'*v+1)^0.5; % noyau polynomail
end

