function x = Lagrange(i, x_interp, x)

L=ones(size(x));
for j=1:length(x_interp)
    if j~=i
        L=L.*(x-x_interp(j))/(x_interp(i)-x_interp(j));
    end
end