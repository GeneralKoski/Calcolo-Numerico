% Esercizio 3
close all
clear
clc
format long 
f=@(x) cos(x);

a=0;
b=pi;
% formula del punto medio semplice
I_ms=f((a+b)/2)*(b-a)
% valore esatto
I_esatto=sin(b)-sin(a)
% errore relativo
Err_rel=abs((I_ms-I_esatto)/I_esatto)

% cambio estremi
a=-pi/2;
b=pi/2;
% formula del punto medio semplice
I_ms=f((a+b)/2)*(b-a)
% valore esatto
I_esatto=sin(b)-sin(a)
% errore relativo
Err_rel=abs((I_ms-I_esatto)/I_esatto)

% formula del punto medio composita

% errore punto medio composita
err_mc=[];
H_mc=[]; 
for m=10:10:100 %numero di sottointervalli della decomposizione
    h=(b-a)/m;
    x=a+h/2:h:b;
    y=f(x);
    int=h*sum(y);
    err_mc=[err_mc; abs(I_esatto-int)];
end
% errore formula trapezi composita Matlab
err_tc=[];
H_tc=[]; 
for m=10:10:100 %numero di sottointervalli della decomposizione
    h=(b-a)/m;
    x=a:h:b;
    y=f(x);
    int=h/2*(y(1)+2*sum(y(2:m))+y(m+1));
    err_tc=[err_tc; abs(I_esatto-int)];
end
format short e
[[10:10:100]' err_mc err_tc]




