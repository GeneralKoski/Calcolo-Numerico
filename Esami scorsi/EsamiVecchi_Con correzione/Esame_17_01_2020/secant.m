function[xvect,nit]=secant(xm1,x0,nmax,toll,fun)
%metodo delle secanti
%INPUT
%xm1 e x0 = estremi del dominio di definizione della funzion fun
%nmax=numero massimo di iterazioni
%toll=tolleranza per il criterio dell'incremento
%fun=funzione di cui si cerca una radice
%OUTPUT
%xvect=successione di valori approssimanti la radice
%xdif=successione di errori assoluti
%fx=successione dei residui
%nit=numero di iterazioni
fxm1=fun(xm1);
fx0=fun(x0);
xvect=[]; 
nit=0;
while nit<nmax && abs(fx0)>toll
    nit=nit+1;
    x=x0-fx0*(x0-xm1)/(fx0-fxm1);
    xvect=[xvect;x];
    xm1=x0;
    fxm1=fx0;
    x0=x;
    fx0=fun(x0);
end
    