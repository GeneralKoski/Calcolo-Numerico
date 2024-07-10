% Esercizio 1
clear all
clc
format long
% matrice di Hilbert di ordine 5
A=hilb(5);
disp(A)
cond(A)
%scambio prima e quinta riga
a=A(1,:);
A(1,:)=A(5,:);
A(5,:)=a;
disp(A)
cond(A)
%scambio seconda e quarta colonna
a=A(:,2);
A(:,2)=A(:,4);
A(:,4)=a;
disp(A)
cond(A)
% matrice di Hilbert di ordine 10
A=hilb(10);
disp(A)
cond(A)