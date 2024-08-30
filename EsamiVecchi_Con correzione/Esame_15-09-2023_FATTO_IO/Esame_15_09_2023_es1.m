clear
clc
hold on

f = @(x) 1./(1+9*x.^2);
a = -1;
b = 1;

i = 1;
for N = 3:2:7
    x_plot = linspace(a, b, 1000);
    y_plot = f(x_plot);
    
    x_equispaced = linspace(a, b, N+1);
    y_equispaced = f(x_equispaced);
    
    %serve solo per lineare semplice
    p_equispaced = polyfit(x_equispaced, y_equispaced, N);
    y_interp_equispaced = polyval(p_equispaced, x_plot);
    %serve solo per lineare composita
    y_interp_equispaced_linear = interp1(x_equispaced, y_equispaced, x_plot, 'linear');

    %grafico lineare semplice
    figure(1);
    subplot(1, 3, i);
    plot(x_plot, y_plot, 'b');
    hold on;
    plot(x_plot, y_interp_equispaced, 'r--');
   
    %grafico lineare composita
    figure(2);
    subplot(1, 3, i);
    plot(x_plot, y_plot, 'b');
    hold on;
    plot(x_plot, y_interp_equispaced_linear, 'r--');
    i = i +1;
end

