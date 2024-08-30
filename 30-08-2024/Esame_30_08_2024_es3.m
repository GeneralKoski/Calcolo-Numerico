clear
clc
format long

M = [4 0 0; 
     2 5 0; 
     0.5 0 3];

[L, U] = lu(M);

Minv = zeros(size(M));

for i = 1:size(M, 1)
    e = zeros(size(M, 1), 1);
    e(i) = 1;
    y = L \ e;
    Minv(:, i) = U \ y;
end

disp('Inversa calcolata con metodo di sostituzione:');
disp(Minv);

Minv_builtin = inv(M);

disp('Inversa calcolata con la funzione built-in:');
disp(Minv_builtin);

difference = Minv - Minv_builtin;
disp('Differenza tra le due inverse:');
disp(difference);
