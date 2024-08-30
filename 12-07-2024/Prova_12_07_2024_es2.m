clear
clc
close all
format long

% Dimensione del sistema
n = 100;

% Costruzione della matrice tridiagonale A
diag_principale = 3 * ones(n, 1); % diagonale principale
diag_secondaria = ones(n-1, 1);    % diagonali secondarie
A = diag(diag_principale) + diag(diag_secondaria, 1) + diag(diag_secondaria, -1);

% Termine noto b
b = [5*ones(1, n-2), 495, 399]';

% Innesco iniziale
x = ones(n, 1);

% Criterio di arresto
toll = 10^-5;

% Decomposizione di A in D e R
D = diag(diag(A));
R = A - D;

% Iterazione del metodo di Jacobi
max_iter = 10000; % numero massimo di iterazioni
for k = 1:max_iter
    x_new = D \ (b - R * x);
    % Controllo del criterio di arresto
    if norm(x_new - x, inf) < toll
        break;
    end
    x = x_new;
end

% Soluzione ottenuta con il metodo di Jacobi
x_jacobi = x;

% Soluzione esatta usando mldivide
x_esatto = A \ b;

% Confronto delle soluzioni
errore_jacobi = norm(x_jacobi - x_esatto, inf);

% Output dei risultati
disp('Errore tra soluzione di Jacobi e soluzione esatta: ');
disp(errore_jacobi);
