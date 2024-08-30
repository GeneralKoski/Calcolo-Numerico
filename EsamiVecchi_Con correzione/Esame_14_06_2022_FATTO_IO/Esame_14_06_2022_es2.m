clear
clc

f=@(s) s.^8 - 2;
a=0.5;
b=2;
x=linspace(a,b,1000);
plot(x,f(x),'k')
hold on
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

coeff = [1 0 0 0 0 0 0 0 2];
s = roots(coeff);
plot(s, f(s), '*');
abs(s(1)-a);
abs(s(2)-b);
