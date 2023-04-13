function kern=kernel(u,v)

global ikernel
if ikernel==1
    kern=u'*v;
elseif ikernel==2
    kern=exp(-0.1*norm(u-v)^2);
elseif ikernel==3
    kern=(u'*v+100)^0.5;
elseif ikernel==4
    kern=exp(-0.44*norm(u-v));%noyau laplacien
end

