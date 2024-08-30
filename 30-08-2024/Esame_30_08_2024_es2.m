clear
clc
format long

n = 1000;
x = linspace(0, 1, n+1); 
f = sqrt(1 - x.^2); 

integrale = trapz(x, f);

pi_approx = 4 * integrale;

disp('Approssimazione di pi:');
disp(pi_approx);

pi_exact = pi;
disp('Approssimazione di pi:');
disp(pi_exact);

n_values = 2:1000;
errors_trapz = zeros(size(n_values));

for i = 1:length(n_values)
    x = linspace(0, 1, n_values(i)+1);
    f = sqrt(1 - x.^2);
    integral_approx = trapz(x, f);
    pi_approx = 4 * integral_approx;
    errors_trapz(i) = abs(pi_exact - pi_approx);
end

figure;
plot(n_values, errors_trapz);
xlabel('Numero di sottointervalli n');
ylabel('Errore assoluto');
title('Errore nel calcolo di pi con il metodo trapezoidale');
grid on;

max_n = 1000;

errors_leibniz = zeros(max_n, 1);

pi_leibniz = 0;
for n = 0:max_n-1
    pi_leibniz = pi_leibniz + ((-1)^n) / (2*n + 1);
    pi_approx = 4 * pi_leibniz;
    errors_leibniz(n+1) = abs(pi_exact - pi_approx);
end

min_length = min(length(n_values), max_n); 
plot(n_values(1:min_length), errors_trapz(1:min_length), 'b', 'DisplayName', 'Trapezoidale');
hold on;
plot(1:min_length, errors_leibniz(1:min_length), 'r', 'DisplayName', 'Leibniz');

xlabel('Numero di termini/sottointervalli n');
ylabel('Errore assoluto');
title('Confronto tra metodo trapezoidale e formula di Leibniz');
legend;
grid on;