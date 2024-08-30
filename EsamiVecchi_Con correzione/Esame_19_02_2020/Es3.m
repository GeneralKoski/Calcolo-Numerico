% Esercizio 3
close all
clear
clc
format long e

% dati
a=sqrt(3);
b=2;
g=@(x) 2*x./sqrt(x.^2+1);

% grafico
ascisse=linspace(a,b);
ordinate=g(ascisse);
plot(ascisse,ordinate)
hold on
pmedio=(a+b)/2;
plot(ascisse,g(pmedio)*ones(size(ascisse)),'r')
legend('g(x)','funzione costante')
xlabel('dominio')
axis square
axis([a b min(ordinate) max(ordinate)])

% formula del punto medio semplice
I_ms=g(pmedio)*(b-a);
% integrazione con "quad"
I_quad=quad(g,a,b);
% confronto
E=abs(I_ms-I_quad);
% % integrazione esatto
% gg=@(x) 2*sqrt(x.^2+1);
% I_esatto=gg(b)-gg(a)

% formula del punto medio composita
err_mc=[];
for m=10:10:100 %numero di sottointervalli della decomposizione
    h=(b-a)/m;
    x=a+h/2:h:b;
    y=g(x);
    int=h*sum(y);
    err_mc = [err_mc; abs(I_quad - int)];
end
format short e
[[10:10:100]; err_mc']


