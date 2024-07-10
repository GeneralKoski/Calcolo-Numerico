function[xvect,nit]=chord(a,b,uno,nmax,toll,fun)
%metodo delle corde
%INPUT
%a=estremo di sinistra intervallo di esistenza della radice
%b=estremo di destra intervallo di esistenza della radice
%x0=ascissa di innesco del metodo
%nmax=numero massimo di iterazioni
%toll=tolleranza per il criterio del residuo
%fun=funzione di cui si cerca uno zero
%OUTPUT
%xvect=successione di valori approssimanti la radice
%fx=successione dei residui
%nit=numero di iterazioni
fa=fun(a);
fb=fun(b);
r=(fb-fa)/(b-a);
err=b-a;
nit=0;
xvect=[]; 
x=uno;
fx=fun(x);
xdif=[];
while nit<nmax && abs(fx)>toll
    nit=nit+1;
    x=x-fx/r;
    xvect=[xvect;x];
    fx=fun(x);
end
    