function int=trapezc(a,b,m,fun)
h=(b-a)/m;
x=a:h:b;
y=fun(x);
int=h/2*(y(1)+2*sum(y(2:m))+y(m+1));