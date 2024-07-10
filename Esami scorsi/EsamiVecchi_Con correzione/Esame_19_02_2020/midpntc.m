function int=midpntc(a,b,m,fun)
h=(b-a)/m;
x=a+h/2:h:b;
y=fun(x);
int=h*sum(y);
