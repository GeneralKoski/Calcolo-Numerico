clear
clc
nodi=linspace(10,20,3);
f=@(x) x-sqrt(2);
y=f(nodi);
Int=trapz(nodi,y)
fint=@(x)(x-sqrt(2))^2/2;
Int_esatto=fint(20)-fint(10)
