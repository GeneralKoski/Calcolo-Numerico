%Esercizio 1
clear
clc
%radici esatte
x1=0.5;
x2=0.0001;
%verifica che l'equazione si annulli nelle radici
f=@(x) x.^2-0.5001*x+5*10^-5;
f(x1)
f(x2)
%confronto con i risultati della function "roots"
coeff=[1 -0.5001 5*10^-5];
s=roots(coeff)
abs(s(1)-x1)
abs(s(2)-x2)
%confronto con i risultati della function "fzero"
sol_1=fzero(f,0)
f=@(x) (x.^2-0.5001*x+5*10^-5)./(x-sol_1);
sol_2=fzero(f,0)