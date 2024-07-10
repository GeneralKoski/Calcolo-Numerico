clear
clc

f=@(x) x.^8-2;
a=0.5;
b=2;
nmax=100;
toll=10^-10;
err=abs(b-a);
nit_bis=1;
s_bis=b;
fx=[];
xdif=[];
while nit_bis<nmax && err>toll
    nit_bis=nit_bis+1;
    c=(a+b)/2;
    x=c;
    fc=f(x);
    s_bis=[s_bis;x];
    fx=[fx;fc];
    x=a;
    if fc*f(x)>0
        a=c;
    else
        b=c;
    end
    err=abs(s_bis(nit_bis)-s_bis(nit_bis-1));
    xdif=[xdif;err];
end
x = linspace(1, 100);
plot(x, f(x), 'blue');
hold on;
plot(xdif, xdif, '');