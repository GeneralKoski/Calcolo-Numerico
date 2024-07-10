% Esercizio 3
close all
clear
clc
format long 

Tab=[];
for n=10:10:100
    %assegnazione della matrice
    A=2*eye(n)-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
    %fattorizzazione di Cholesky
    B=chol(A);
    %verifica della correttezza della fattorizzazione
    prova = norm(A-B'*B,inf);
    %costruzione del termine noto affinchè la soluzione sia un vettore di
    %elementi pari ad 1
    x=ones(n,1);
    b=A*x;
    %risoluzione del sistema lineare con function Matlab
    x_1=A\b;
    % Metodo sostituzione in avanti per risolvere B'y=b
    y = Forward(B', b);
    err_1=norm(x-x_1,inf);
    % Metodo sostituzione all'indietro per risolvere By=b
    x_2=Backward(B,y);
    err_2=norm(x-x_2,inf);
    %confronto degli errori
    Tab=[Tab; err_1 err_2];
end

Tab




