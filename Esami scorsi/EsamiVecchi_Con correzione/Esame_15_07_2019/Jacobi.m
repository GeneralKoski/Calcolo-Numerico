function[x,iter,r]=Jacobi(a,b,x0,nmax,toll)
%Metodo di Jacobi
[n,n]=size(a);
iter=0;
r=b-a*x0;
r0=norm(r);
err=norm(r);
while err>toll && iter<nmax
    iter=iter+1;
    for i=1:n
        s=0;
        for j=1:i-1
            s=s+a(i,j)*x0(j);
        end
        for j=i+1:n
            s=s+a(i,j)*x0(j);
        end
        x(i,1)=(b(i)-s)/a(i,i);
    end
%     D=diag(diag(a));
%     C=a-D;
%     B_Jacobi=-inv(D)*C
%     max(abs(eig(B_Jacobi)))
%     pause
%     x_bis=-inv(D)*C*x0+inv(D)*b;
%     [x_bis,x]
%     pause
    x0=x;
    r(iter)=norm(b-a*x);
    err=r(iter)/r0;
end
    

