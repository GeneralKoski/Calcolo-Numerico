% Esercizio 3
close all
clear
clc
format long 
f=@(x) 1+exp(-x.^2)/2;
a=-1;
b=1;
% formula trapezi semplice
I_ts=(f(a)+f(b))*(b-a)/2
% formula trapezi composita
m=10;
h=(b-a)/m;
x=a:h:b;
y=f(x);
int=h/2*(y(1)+2*sum(y(2:m))+y(m+1));
x=linspace(a,b,m+1);
I_trapz=trapz(x,f(x)) %formula trapezi composita Matlab
I_trapz-int; %verifica correttezza
I_quad=quad(f,a,b) %quadratura adattiva Matlab
% errore
err=[];
H=[]; 
for m=10:10:100
    h=(b-a)/m;
    x=a:h:b;
    y=f(x);
    int=h/2*(y(1)+2*sum(y(2:m))+y(m+1));
    err=[err; abs(I_quad-int)];
    H=[H;(b-a)/m]%ampiezza sottointervalli
end
format short e
[[10:10:100]' H err err./H.^2]




