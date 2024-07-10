% Esercizio 3
close all
clear
clc
format long e
for i=1:10
    % costruzione della matrice
    n=10*i;
    U=rand(n);
    U=triu(U);
    % costruzione del termine noto
    sol=ones(n,1);
    b=U*sol;
    % risoluzione con / 
    t=tic;
    x=U\b;
    time_x(i)=toc(t);
    err_x(i)=norm(x-sol,inf);
    % risoluzione con sostituzione all'indietro
    t=tic;
    xx=Backward(U,b);
    time_xx(i)=toc(t);
    err_xx(i)=norm(xx-sol,inf);
    condiz(i)=cond(U);
end
format short e
[[1:10]' condiz' err_x' time_x' err_xx' time_xx']


