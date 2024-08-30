clear
clc
close all

A = [8, 2, 3; 2, 4 ,0; 3, 0, 5];
b = [48;2;35];
x = [0 0 0]';
toll = 10^-6;

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

semilogy([1:iter],err_J);
hold on;
semilogy(err_J, err_J, '*');
xlabel('numero di iterazioni')
ylabel('errore in norma Inf sulla soluzione esatta')