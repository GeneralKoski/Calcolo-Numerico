% Esercizio 2
clear
clc
format rat
f=@(x) 3*log(1./x)+sqrt(5*x);
x=[1 3 7];
y_0=f(x);
y=round(y_0);
%costruisco il polinomio interpolatore con 3 punti
p=polyfit(x,y,2)
xx=linspace(min(x),max(x));
yy_0=f(xx);
yy=polyval(p,xx);
plot(xx,yy_0,'k',x,y_0,'k*',xx,yy,'b',x,y,'bo')
pause
close all
x=[1 3 7 12];
y_0=f(x);
y=round(y_0);
%costruisco il polinomio interpolatore con 4 punti
p=polyfit(x,y,3)
xx=linspace(min(x),max(x));
yy_0=f(xx);
yy=polyval(p,xx);
plot(xx,yy_0,'k',x,y_0,'k*',xx,yy,'b',x,y,'bo')