clear
clc

n = 3;
%A = rand(n);
A = [2 2 1; 0 1 0; 1 0 2];
x = (1:n)';
b=A*x;

A\b

y = ones(n,1);
D = diag(diag(A));
E = tril(A,-1);
F = triu(A,1);
C = A-D;
toll = 10^(-4);
BGS=-(D+E)^(-1)*F
N1=norm(BGS,1);
N2=norm(BGS);
Ninf=norm(BGS,inf);
aut=eig(BGS);
rspett=max(abs(aut))
pause
count = 0;
while norm(A*y-b, inf) > toll
    
    z = -(D+E)^(-1)*F*y+(D+E)^(-1)*b;
    y = z;
    count = count +1;
    
end
count
y
    