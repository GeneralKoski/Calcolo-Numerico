clear
clc

% Punti dati
x = [-2, -5, 0, 0.1, pi];
y = [-2, sqrt(2), 7, 7, 7];

% Calcolo delle somme
S_xx = sum(x.^2);
S_xy = sum(x .* y);
S_x = sum(x);
S_y = sum(y);
n = length(x);

% Costruzione del sistema lineare
A = [S_xx, S_x; S_x, n];
b = [S_xy; S_y];

% Risoluzione del sistema lineare per trovare a e b
coeffs = A \ b;
a = coeffs(1);
b = coeffs(2);

% Stampa dei coefficienti
disp(['a = ', num2str(a)]);
disp(['b = ', num2str(b)]);

% Polinomio interpolante usando il metodo di Lagrange

L = 0;
for j = 1:n
    L_j = 1;
    for m = 1:n
        if m ~= j
            L_j = L_j .* (x - x(m)) / (x(j) - x(m));
        end
    end
    L = L + y(j) * L_j;
end


% Generare punti per il grafico
x_values = linspace(min(x)-1, max(x)+1, 500);
y_lagrange = double(subs(p(x), x, x_values));
y_least_squares = a * x_values + b;

% Grafico
figure;
hold on;

% Polinomio interpolante
plot(x_values, y_lagrange, 'b-', 'DisplayName', 'Polinomio Interpolante');

% Retta ai minimi quadrati
plot(x_values, y_least_squares, 'r--', 'DisplayName', 'Retta ai Minimi Quadrati');

% Punti dati
scatter(x, y, 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Punti di Interpolazione');

% Legenda
legend show;
xlabel('x');
ylabel('y');
title('Polinomio Interpolante e Retta ai Minimi Quadrati');
grid on;

hold off;
