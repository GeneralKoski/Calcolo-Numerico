%Esercizio 3
format long
close all
clear 
clc
%Inserimento dati
A=diag(15*ones(7,1))+diag(1*ones(6,1),1)+diag(1*ones(6,1),-1);
x=[1 0 1 0 1 0 1]';
b=A*x;
%Metodo di Jacobi
[n,n]=size(A);
x0=zeros(n,1);
iter=0;
nmax=40;
while iter<nmax
    iter=iter+1;
    D=diag(diag(A));
    C=A-D;
    B_Jacobi=-inv(D)*C;
    x_J(:,iter)=-inv(D)*C*x0+inv(D)*b;
    x0=x_J(:,iter);
    err_J(iter,1)=norm(x-x_J(:,iter),inf);
end
disp('norma della matrice di Jacobi')
norm(B_Jacobi,inf)
semilogy([1:iter],err_J)
xlabel('numero di iterazioni')
ylabel('errore in norma Inf sulla soluzione esatta')
%Metodo di Gauss-Seidel
x0=zeros(n,1);
iter=0;
x=[];
x_old=x0;
while iter<nmax
    iter=iter+1;
    D=diag(diag(A));
    E=tril(A,-1);
    F=triu(A,1);
    B_Gauss=-inv(D+E)*F;
    x(:,iter)=-inv(D+E)*F*x_old+inv(D+E)*b;
    x_old=x(:,iter);
end

for iter=1:nmax
    err_GS(iter,1)=norm(x_old-x(:,iter),inf);
end
disp('norma della matrice di Gauss-Seidel')
norm(B_Gauss,inf)
hold on
semilogy([1:iter],err_GS,'r')    
legend('Jacobi','Gauss-Seidel')
disp('raggi spettrali Jacobi e Gauss-Seidel')
[max(eig(B_Jacobi)) max(eig(B_Gauss))]

