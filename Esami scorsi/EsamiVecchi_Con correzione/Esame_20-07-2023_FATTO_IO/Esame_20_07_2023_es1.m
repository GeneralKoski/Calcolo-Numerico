close all
clear
clc

x_interp=[-1 0 1 2];
y_interp=[-0.1 2.03 3.5 6];
a=-1;
b=6;
x=linspace(a,b,1000);
y=zeros(size(x));
for i=1:3
    y=y+y_interp(i)*Lagrange(i,x_interp,x);
end
%grafico
plot(x,y,'b')
hold on
plot(x_interp, y_interp,'ob')
axis square
grid on
xlabel('x')
ylabel('y')