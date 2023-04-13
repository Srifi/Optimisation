% Generation de points pour fonction de R dans R à regresser avec SVR

X=[-3 -1 0 1 2  4  6   7 8 10  12  13 13.3  15 16 16.5 17  18 20 ];
 npoints=length(X);
% choix entre plusieurs valeurs de Y . La première: fonction affine

%Y= 0.5*X+3+ 1.5*(0.5-rand(1,npoints)); %points à approximer
%Y= X/2+(cos(X/3)).*X+  1.5*(0.5-rand(1,npoints));
Y= (sin(X/3)).*X+  1.5*(0.5-rand(1,npoints));

save 'courb1' X Y