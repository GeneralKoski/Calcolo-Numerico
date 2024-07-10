% Esercizio 1
clear
clc
format long
n=10:10:100;
x=2.^n;
% prima espressione
lim_1=4/pi*x.*(sqrt(x.^2+pi)-x);
% seconda espressione
lim_2=4*x./(sqrt(x.^2+pi)+x);
% confronto
[lim_1' lim_2']
