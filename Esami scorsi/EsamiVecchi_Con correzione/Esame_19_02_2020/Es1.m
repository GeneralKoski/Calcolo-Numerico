% Esercizio 1
clear
clc
format long
x=0.9975; 
% prima espressione
q=(1-x)^4
% sviluppiamo la potenza del binomio
syms z real 
rr=expand((1-z).^4)
% seconda espressione
r=x^4-4*x^3+6*x^2-4*x+1
%double(subs(rr,x))

