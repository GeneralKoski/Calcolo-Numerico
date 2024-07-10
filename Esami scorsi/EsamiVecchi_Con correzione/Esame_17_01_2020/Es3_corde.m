clear
clc
f=@(s) 2*s-sqrt(s.^2+1);
a = 0;
b = 2;
uno = 1;

% approssimazione con il metodo delle corde
nmax=100;
toll=10^-10;
root=fzero(@(s) f(s),1)
fa=f(a);
fb=f(b);
r=(fb-fa)/(b-a);
err=b-a;
nit=0;
xvect=[]; 
x=uno;
fx=f(x);
xdif=[];
while nit<nmax && abs(fx)>toll
    nit=nit+1;
    x=x-fx/r;
    xvect=[xvect;x];
    fx=f(x);
end
disp('tabella metodo corde')
[[1:nit]' xvect abs(f(xvect))]
disp('confronto corde-FZERO')
abs(root-xvect(end))
% approssimazione con il metodo delle secanti
fa=f(a);
fb=f(b);
xvect=[]; 
nit=0;
while nit<nmax && abs(fb)>toll
    nit=nit+1;
    x=b-fb*(b-a)/(fb-fa);
    xvect=[xvect;x];
    a=b;
    fa=fb;
    b=x;
    fb=f(b);
end

disp('tabella metodo secanti')
[[1:nit]' xvect abs(f(xvect))]
% confronto con FZERO

disp('confronto secanti-FZERO')
abs(root-xvect(end))