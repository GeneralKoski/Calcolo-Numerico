clear
clc

% Dati forniti
y = [12.51 13.05 11.7 9.26 7.755 6.25 5.34 4.59 5.14 6.36 10.31 13.88];

% Il valore mancante (per Maggio, posizione 5)
missing_index = 5;

% Rimuovere il valore mancante dai dati
x = 1:length(y),100;
x_nuovo = x;
x_nuovo(missing_index) = [];
y_nuovo = y;
y_nuovo(missing_index) = [];

% Creazione del piano cartesiano
figure;
hold on;

% Polinomiale lineare a tratti
linearFit = polyfit(x_nuovo, y_nuovo, 1);
linearCurve = polyval(linearFit, x);
plot(x, linearCurve, '--', 'DisplayName', 'Polinomiale Lineare a Tratti');

% Polinomio di interpolazione (grado n-2 poiché un dato è mancante)
interpolationPoly = polyfit(x_nuovo, y_nuovo, length(x) - 2);
interpolationCurve = polyval(interpolationPoly, x);
plot(x, interpolationCurve, '-', 'DisplayName', 'Polinomio di Interpolazione');

% Parabola approssimante
parabolicFit = polyfit(x_nuovo, y_nuovo, 2);
parabolicCurve = polyval(parabolicFit, x);
plot(x, parabolicCurve, '-.', 'DisplayName', 'Parabola Approssimante');

%Trovo i valori mancanti dei vari dei vari polinomi
mancante_lineare = polyval(linearFit, missing_index);
mancante_interpolante = polyval(interpolationPoly, missing_index);
mancante_parabolica = polyval(parabolicFit, missing_index);

% Aggiungere i valori stimati al grafico
y(missing_index) = mancante_interpolante;
plot(x, y, 'o', 'DisplayName', 'Dati Raccolti');
plot(missing_index, mancante_lineare, 'bx', 'DisplayName', 'Valore Stimato (Lineare)');
plot(missing_index, mancante_interpolante, 'rx', 'DisplayName', 'Valore Stimato (Interpolante)');
plot(missing_index, mancante_parabolica, 'yx', 'DisplayName', 'Valore Stimato (Parabola)');

% Aggiunta delle legende e del titolo
legend('show');
title('Grafico dei Dati e Curve Approssimanti');
xlabel('Mese');
ylabel('Portata (m^3)');
grid on;
hold off;
