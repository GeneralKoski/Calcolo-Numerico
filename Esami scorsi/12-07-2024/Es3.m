clear
clc
close all
format long

% Definizione dei punti di ascissa
x = [-1, -0.8, -0.5, -0.4, -0.3, -0.2, 0.2, 0.3, 0.4, 0.5, 0.8, 1];

% Funzione da integrare
f = @(x) sqrt(1 - x.^2);

% Calcolo della quadratura del punto medio composta
n = length(x) - 1;
integrale_punto_medio = 0;

for i = 1:n
    punto_medio = (x(i) + x(i+1)) / 2;
    integrale_punto_medio = integrale_punto_medio + f(punto_medio) * (x(i+1) - x(i));
end

disp('Integrale con quadratura del punto medio composta: ');
disp(integrale_punto_medio);

% Calcolo dell'integrale usando la formula dei trapezi composta
y = f(x);
integrale_trapezio = trapz(x, y);

% Calcolo dell'errore
errore_medio_trapezio = abs(integrale_punto_medio - integrale_trapezio);

disp('Errore tra la quadratura del punto medio e i trapezi: ');
disp(errore_medio_trapezio);

% Valore esatto dell'integrale
integrale_esatto = pi / 2;

% Calcolo degli errori rispetto al valore esatto
errore_medio_esatto = abs(integrale_punto_medio - integrale_esatto);
errore_trapezio_esatto = abs(integrale_trapezio - integrale_esatto);

disp('Errore quadratura del punto medio rispetto al valore esatto: ');
disp(errore_medio_esatto);
disp('Errore formula dei trapezi rispetto al valore esatto: ');
disp(errore_trapezio_esatto);
