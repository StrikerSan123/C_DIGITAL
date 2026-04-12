clc; clear; close all;

% Cargar datos
data = readtable('datamotor1.csv');

% Extraer señales (filtrando ceros si aplica)
idx_valid = find(data.ex_signal ~= 0);
time = data.time(idx_valid);
in = data.ex_signal(idx_valid);
out = data.system_response(idx_valid);

% Ajustar tiempo para que empiece en 0
time = time - time(1);
data_F = table(time, in, out);

% Sacar datos mediante media (estado estacionario)
maximos = data_F.out(data_F.time >= 5.5);
media = mean(maximos);
senal100 = media * ones(size(time));

% Sacar señal base
min_out = min(out);
senal0 = min_out * ones(size(time));

% 1. MÉTODO DE ZIEGLER Y NICHOLS (Tangente)
% Puntos arbitrarios para pendiente elegidos por ti
x1 = 0.05;
x2 = 0.1;

[~, idx1] = min(abs(data_F.time - x1));
[~, idx2] = min(abs(data_F.time - x2));

y1 = data_F.out(idx1);
y2 = data_F.out(idx2);
pendiente = (y2 - y1) / (x2 - x1);

% Intersecciones de la recta tangente
xlow = (min_out - y1) / pendiente + x1; 
if xlow < 0; xlow = 0; end % El tiempo muerto no puede ser negativo

xhigh = (media - y1) / pendiente + x1;

t_tan = linspace(xlow, xhigh, 100);
rtan = pendiente * (t_tan - x1) + y1;

% Parámetros Z-N
K = media / max(in);
theta_zn = xlow;
tau_zn = xhigh - xlow;
G_zn = tf(K, [tau_zn 1], 'InputDelay', theta_zn);
y_zn = lsim(G_zn, in, time);

% PREPARACIÓN PARA MILLER Y ANALÍTICO (Evitar error en interp1)
% Extraemos solo la curva de subida para que los datos no se repitan
idx_subida = find(data_F.out >= media, 1); 
out_subida = data_F.out(1:idx_subida);
time_subida = data_F.time(1:idx_subida);

% 2. MÉTODO DE MILLER
ymiller = 0.6321 * media;

% Interpolación lineal usando solo la curva de subida
t63 = interp1(out_subida, time_subida, ymiller);

taumiller = t63 - theta_zn; % Miller asume el theta calculado por tangente o inspección
G_miller = tf(K, [taumiller 1], 'InputDelay', theta_zn);
y_miller = lsim(G_miller, in, time);

% 3. MÉTODO ANALÍTICO (2 puntos: 28.4% y 63.2%)
yan = 0.284 * media;

% Interpolación de los dos puntos
x63 = interp1(out_subida, time_subida, ymiller);
x28 = interp1(out_subida, time_subida, yan);

tauan = 1.5 * (x63 - x28);
thetan = x63 - tauan;
if thetan < 0; thetan = 0; end % Prevenir tiempo muerto negativo

G_an = tf(K, [tauan 1], 'InputDelay', thetan);
y_an = lsim(G_an, in, time);

% CÁLCULO DEL ERROR CUADRÁTICO MEDIO (RMSE)
% RMSE = sqrt( mean( (ValorReal - ValorSimulado)^2 ) )
rmse_zn = sqrt(mean((out - y_zn).^2));
rmse_miller = sqrt(mean((out - y_miller).^2));
rmse_an = sqrt(mean((out - y_an).^2));

% Imprimir resultados en consola
fprintf('--- RMSE de los Modelos ---\n');
fprintf('Ziegler-Nichols: %.4f\n', rmse_zn);
fprintf('Miller:          %.4f\n', rmse_miller);
fprintf('Analítico:       %.4f\n\n', rmse_an);

% GRÁFICAS
figure('Position', [100, 100, 900, 600]);

% Señales reales
plot(time, in, 'Color', "k", 'LineWidth', 1.5); hold on;
plot(time, out, 'Color', '#FF5733', 'LineWidth', 2);

% Líneas de referencia y tangente
plot(time, senal100, '--k', 'LineWidth', 1);
plot(time, senal0, '--r', 'LineWidth', 1);
plot(t_tan, rtan, '--b', 'LineWidth', 1.5);

% Simulaciones de los modelos
plot(time, y_zn, 'Color', "#9200FF", 'LineWidth', 1.5);
plot(time, y_miller, 'Color', "#EDB000", 'LineWidth', 1.5);
plot(time, y_an, 'Color', "#4DBFEE", 'LineWidth', 1.5);

% Formato del gráfico
title('Ajuste de Modelos de Primer Orden y Comparación RMSE');
xlabel('Tiempo [s]');
ylabel('Amplitud');

% Creación de leyendas dinámicas con los valores de RMSE
str_zn = sprintf('Ziegler y Nichols (RMSE: %.3f)', rmse_zn);
str_miller = sprintf('Miller (RMSE: %.3f)', rmse_miller);
str_an = sprintf('Analítico (RMSE: %.3f)', rmse_an);

legend('Entrada (Escalón)', 'Salida Real', ...
    ['Estado Estacionario: ' num2str(media, '%.2f')], ...
    ['Línea Base: ' num2str(min_out)], ...
    'Recta Tangente', str_zn, str_miller, str_an, 'Location', 'southeast');
grid on;


%% GRAFICA XY

data = readtable('rampa');

figure;
plot(data.V, data.M, 'LineWidth', 1.5);

title('Operación lineal del motor: velocidad vs voltaje');
xlabel('Voltaje [V]');
ylabel('Velocidad [rpm x1000]');

grid on;




