% dessin nappe

global ikernel 
global dimX
 
if dimX==1
    load 'courbe' X Y
else
    load 'nappe' X Y
end
load 'Solution' alph b


%load 'nappe' X Y
%load 'Solution' alph b

ikernel=2;
dimX=size(X,1)

% dessin des points 

if dimX==1
    plot(X,Y,'o')
    grid
else
   plot3(X(1,:),X(2,:),Y, 'o')
end
   axis equal
   hold on

   
   %tracé de la fonction de régression


if dimX==1
xmin=min(X);
xmax=1.05*max(X);
xp=xmin:0.2:xmax;   % grille pour plot 
npx=length(xp);
for i=1:npx
        yp(i)=prodwr(X,alph,xp(i))+b;
end
plot(xp,yp)
grid on


else

ikx=0;
ik=0;
for ix=-12:1:12
     ikx=ikx+1;
     iky=0;
     for iy=-12:1:12
         ik=ik+1;
         iky=iky+1;
         XS(1,ik)=ix;
         XS(2,ik)=iy;
         YS(ik)= prodwr(X,alph,XS(:,ik))+b;
         TXS1(ikx,iky)=XS(1,ik);
         TXS2(ikx,iky)=XS(2,ik);
         TYS(ikx,iky)=YS(ik);
         
     end
end

grid
 surf(TXS1,TXS2,TYS) 
grid
end