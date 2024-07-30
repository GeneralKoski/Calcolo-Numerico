clear
clc
close all
format long

% Punti dati
x = [-1, 1, sqrt(2), 2];
y = [tan(2), 1, 2*sqrt(2), exp(2)];

% a) Polinomio interpolatore semplice
p = polyfit(x, y, length(x) - 1);

% Valutazione del polinomio interpolatore semplice in x = 1.5
y_interp_semplice = polyval(p, 1.5);

% Grafico del polinomio interpolatore semplice
x_plot = linspace(-1, 2, 1000);
y_plot = polyval(p, x_plot);
plot(x_plot, y_plot, '-r', 'DisplayName', 'Polinomio Interpolatore Semplice');
hold on;

% b) Retta ai minimi quadrati
p1 = polyfit(x, y, 1);

% Valutazione della retta ai minimi quadrati in x = 1.5
y_min_quad = polyval(p1, 1.5);

% Grafico della retta ai minimi quadrati
y_plot_p1 = polyval(p1, x_plot);
plot(x_plot, y_plot_p1, '-b', 'DisplayName', 'Retta ai Minimi Quadrati');

% c) Polinomio interpolatore composito di grado 1 (interpolazione lineare)
y_interp_composito = interp1(x, y, x_plot, 'linear');

% Valutazione del polinomio interpolatore composito di grado 1 in x = 1.5
y_interp_comp_1 = interp1(x, y, 1.5, 'linear');

% Grafico del polinomio interpolatore composito di grado 1
plot(x_plot, y_interp_composito, '-g', 'DisplayName', 'Polinomio Interpolatore Composito Grado 1');

% Punti dati
plot(x, y, 'ok', 'MarkerFaceColor', 'k', 'DisplayName', 'Punti dati');

% Aggiungere legende e titoli
xlabel('x');
ylabel('y');
title('Approssimazioni e Interpolazioni');
legend show;
grid on;
hold off;

% d) Valori delle funzioni approssimanti in x = 1.5
disp('Valore del polinomio interpolatore semplice in x = 1.5: ');
disp(y_interp_semplice);
disp('Valore della retta ai minimi quadrati in x = 1.5: ');
disp(y_min_quad);
disp('Valore del polinomio interpolatore composito di grado 1 in x = 1.5: ');
disp(y_interp_comp_1);
