clear
clc

n = 10000;
A = zeros(n) + eye(n) + diag(ones(n-1, 1), 1) * 2;
disp(norm(A, 1));
disp(norm(A, inf));
cond = cond(A);