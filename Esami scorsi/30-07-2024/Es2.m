clear
clc
close all
format long

x = linspace(-2, 2, 1000); % Genera 1000 punti tra -2 e 2
y = 5.*x.^6 - 2.*x - 7;
plot(x, y);
grid on;
xlabel('x');
ylabel('5x^6 - 2x - 7');
title('Grafico della funzione 5x^6 - 2x - 7');

% Funzione da trovare la radice
f = @(x) 5*x^6 - 2*x - 7;

% Intervallo iniziale
a = 1;
b = 2;

% Tolleranza
tol = 10^-5;

% Metodo di bisezione
while (b - a) / 2 > tol
    c = (a + b) / 2;
    if f(c) == 0
        break;
    elseif f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
end

% Radice approssimata
root = (a + b) / 2;
disp('La radice approssimata è: ');
disp(root);

% Intervallo iniziale
x0 = [1, 2];

% Trovare la radice usando fzero
root_fzero = fzero(f, x0);

% Visualizzare il risultato
disp('La radice approssimata con fzero è: ');
disp(root_fzero);

% Confrontare le due radici
differenza = abs(root - root_fzero);
disp('La differenza tra le due radici è: ');
disp(differenza);
