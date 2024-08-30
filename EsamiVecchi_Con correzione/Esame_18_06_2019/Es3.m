% Esercizio 3
close all
clear
clc
format long e
% costruzione del polinomio interpolatore attraverso la matrice di 
% Vandermonde relativa a 3 punti di interpolazione
x_interp=[-1.4,0,0.1];
y_interp=[1.450850792298967,1,1.110719891073862]';
matrix=vander(x_interp);
coeff=matrix\y_interp;
%grafico
a=-1.5;
b=1.5;
x=linspace(a,b,1000);
y=polyval(coeff,x);
plot(x,y,'b')
hold on
plot(x_interp,y_interp,'ob')
axis square
grid on
xlabel('x')
ylabel('y')
% costruzione del polinomio interpolatore attraverso la matrice di 
% Vandermonde relativa a 4 punti di interpolazione
x_interp=[x_interp,1];
y_interp=[y_interp;5.031038733198447];
matrix=vander(x_interp);
coeff=matrix\y_interp;
y=polyval(coeff,x);
hold on
plot(x,y,'b')
plot(x_interp,y_interp,'ob')
% grafico della funzione
f=@(x) exp(x)./cos(x);
plot(x,f(x),'k')
legend('interpolazione a 3 punti','nodi','interpolazione a 4 punti','nodo aggiunto','funzione')
% costruzione del polinomio interpolatore attraverso la matrice di 
% Vandermonde con un numero crescente di punti di interpolazione
for n=10:10:100
    x_interp=linspace(a,b,n);
    y_interp=f(x_interp)';
    matrix=vander(x_interp);
    condiz(n/10)=cond(matrix);
    coeff=matrix\y_interp;
    y=polyval(coeff,x);
    y_bis=interp1(x_interp,y_interp,x);
    err(n/10)=norm(y-f(x),inf);
    err_bis(n/10)=norm(y_bis-f(x),inf);
end
disp('condizionamento')
[[10:10:100]' condiz']
disp('errore')
[[10:10:100]' err' err_bis']





