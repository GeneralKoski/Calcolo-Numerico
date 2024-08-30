clear
clc

%assegnazione della matrice
Tab=[];
A=[1 2 3; 0 4 5; 0 0 6];
n = 3;
%costruzione del termine noto affinchè la soluzione sia un vettore di
%elementi pari ad 1
x=ones(n,1);
b=A*x;
%risoluzione del sistema lineare con function Matlab
% Metodo sostituzione all'indietro per risolvere By=b
A=Backward(A,b);
err_2=norm(x-A,inf);
%confronto degli errori
Tab=[Tab; err_2];