% Esercizio 3
close all
clear
clc
format long e
% costruzione del polinomio interpolatore relativa a 3 punti di 
% interpolazione con i polinomi fondamentali di Lagrange
x_interp=[0.05,1.0,1.5];
y_interp=[-1.997500520789933*10,-5.403023058681398*10^-1,-4.715813444513527*10^-2];
a=0.05;
b=1.5;
x=linspace(a,b,1000);
y=zeros(size(x));
for i=1:3
    y=y+y_interp(i)*Lagrange(i,x_interp,x);
end
%grafico
plot(x,y,'b')
hold on
plot(x_interp,y_interp,'ob')
axis square
grid on
xlabel('x')
ylabel('y')
% grafico della funzione
f=@(x) -cos(x)./x;
plot(x,f(x),'k')
% grafico dell'interpolante ottenuto con polyfit
P=polyfit(x_interp,y_interp,length(x_interp)-1);
yy=polyval(P,x);
plot(x,yy,'r')
legend('interpolazione a 3 punti','nodi','funzione','polyfit')
% costruzione del polinomio interpolatore con polyfit e nodi equispaziati
% figure
% plot(x,f(x),'k')
% hold on
for n=5:5:50
    x_interp=linspace(a,b,n);
    y_interp=f(x_interp);
    P=polyfit(x_interp,y_interp,n-1);
    yy=polyval(P,x);
%     plot(x,yy)
    err(n/5)=norm(yy-f(x),inf);
end
% costruzione del polinomio interpolatore Lagrange e nodi equispaziati
% figure
% plot(x,f(x),'k')
% hold on
L=ones(size(x));
for j=1:length(x_interp)
    if j~=i
        L=L.*(x-x_interp(j))/(x_interp(i)-x_interp(j));
    end
end

for n=5:5:50
    y=0;
    x_interp=linspace(a,b,n);
    y_interp=f(x_interp);
    for i=1:n
        y=y+y_interp(i)*L;
    end
%     plot(x,y)
%     plot(x_interp,y_interp,'o')
    err_L(n/5)=norm(y-f(x),inf);
end