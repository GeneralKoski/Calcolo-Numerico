function L=Lagrange(i,X,x)
%i-esimo polinomio fondamentale di Lagrange valutato nel punto x, con X
%ascisse dei punti di interpolazione
L=ones(size(x));
for j=1:length(X)
    if j~=i
        L=L.*(x-X(j))/(X(i)-X(j));
    end
end
