function[xvect,nit]=Newton(x0,nmax,toll,fun,dfun)
%metodo di Newton
%INPUT
%x0=ascissa di innesco del metodo
%nmax=numero massimo di iterazioni
%toll=tolleranza per il criterio dell'incremento
%fun=funzione di cui si cerca uno zero
%dfun=derivata di fun
%OUTPUT
%xvect=successione di valori approssimanti la radice
%nit=numero di iterazioni
err=1;
nit=0;
xvect=x0; 
x=x0;
fx=fun(x);
xdif=[];
while nit<nmax && err>toll
    nit=nit+1;
    x=xvect(nit);
    dfx=dfun(x);
    if dfx==0
        disp('arresto per azzeramento dfun')
    else
    xn=x-fx(nit)/dfx;
    err=abs(xn-x);
    xdif=[xdif;err];
    x=xn;
    xvect=[xvect;x];
    fx=[fx;fun(x)];
    end
end
    