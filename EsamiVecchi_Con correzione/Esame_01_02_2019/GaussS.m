function[x,B_Gauss]=GaussS(a,b,x0,nmax)
%Input
% a = matrice
% b = termine noto
% x0 = vettore iniziale
% nmax = numero massimo di iterazioni

%Output
% x = matrice contenente nella i-esima colonna il vettore 
%     soluzione all'i-esima iterazione
% B_Gauss = matrice di iterazione del metodo

%Metodo di Gauss-Seidel
[n,n]=size(a);
iter=0;
x=[];
x_old=x0;
while iter<nmax
    iter=iter+1;
    D=diag(diag(a));
    E=tril(a,-1);
    F=triu(a,1);
    B_Gauss=-inv(D+E)*F;
    x(:,iter)=-inv(D+E)*F*x_old+inv(D+E)*b;
    x_old=x(:,iter);
end