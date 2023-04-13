% Génération de points de R3 pour interpolation

 ik=0;
 ikx=0;
 iky=0;
 for ix=-12:4:12
     ikx=ikx+1;
     iky=0;
     for iy=-12:4:12
         ik=ik+1;
         iky=iky+1;
         X(1,ik)=ix+3*rand;
         X(2,ik)=iy +3*rand;
         Y(ik)= ( sin(norm(X(:,ik)/3))  )*norm(X(:,ik)/3)^1.5 +  0.3*(0.5-rand) ;% bruitage
         
     end
 end
 save 'nappe1' X Y