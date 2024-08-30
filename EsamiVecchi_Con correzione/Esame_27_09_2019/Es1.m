% Esercizio 1
clear all
clc
format long

% punto 1a)
eps
N=10^-16;
M=10+N;
M>10
% punto 1b)
S(1)=0;
for i=1:2000
    S(i+1)=S(i)+1/100;
end
S(end);