% Esercizio 1
clear
clc
format long
s=5.555*10^9;
% prima espressione
g=s-sqrt(s^2-1)
% seconda espressione
p=1/(s+sqrt(s^2-1))
% le due espressioni sono uguali
syms x real 
simplify(x-sqrt(x^2-1)-(1/(x+sqrt(x^2-1))))
