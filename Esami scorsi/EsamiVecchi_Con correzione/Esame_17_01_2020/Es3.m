% Esercizio 3
close all
clear
clc
format long e
% grafico intersezione
r=@(s) 2*s;
z=@(s) sqrt(s.^2+1);
a=0;
b=2;
s=linspace(a,b);
plot(s,r(s),'k',s,z(s),'r')
xlabel('s')
axis square
grid on
hold on
% dati
f=@(s) 2*s-sqrt(s.^2+1);%2*s-sqrt(s.^2+1);
a=0;
b=2;
% intersezione approssimata con il comando fzero
disp('ascissa intersezione calcolata con il comando FZERO')
root=fzero(@(s) f(s),1)
plot(root,r(root),'*r')
legend('r(s)','z(s)','approssimazione intersezione')


