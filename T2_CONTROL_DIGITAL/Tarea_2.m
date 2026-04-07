clear; clc; close all;

%% 1. CARGA DE DATOS
data = readmatrix('valveProp.txt');
V_control = data(:,1); 
posicion  = data(:,2); 
flujo     = data(:,3); 
Ts = 1;
tiempo = (0:length(V_control)-1)' * Ts;

%% 2. IDENTIFICACIÓN SEGUNDO SALTO (Puntos de Operación)
indices_subida = find(diff(V_control) > 2);
if length(indices_subida) < 2
    error('No se detectaron suficientes saltos de voltaje.');
end

idx_s = indices_subida(2); 
t_s = tiempo(idx_s);

du = 4; % Salto de 2V a 6V
y_ini_f = mean(flujo(idx_s-15:idx_s));
y_fin_f = mean(flujo(idx_s+55:idx_s+75)); 
dy_f = y_fin_f - y_ini_f;
K_f = dy_f / du;

%% 3. OBTENCIÓN DE MODELOS (Punto 3)

% --- MÉTODO A: MILLER (28.3% y 63.2%) ---
v28 = y_ini_f + 0.283*dy_f; 
v63 = y_ini_f + 0.632*dy_f;
t28 = tiempo(find(flujo(idx_s:end) >= v28, 1) + idx_s - 1);
t63 = tiempo(find(flujo(idx_s:end) >= v63, 1) + idx_s - 1);

tau_m = 1.5 * (t63 - t28);
L_m = max(0, (t63 - t_s) - tau_m);
G_miller = tf(K_f, [tau_m 1], 'InputDelay', L_m);

% --- MÉTODO B: ANALÍTICO (Smith 63.2%) ---
t_start = tiempo(find(flujo(idx_s:end) >= y_ini_f + 0.02*dy_f, 1) + idx_s - 1);
L_a = max(0, t_start - t_s);
tau_a = max(0.1, (t63 - t_s) - L_a);
G_analitico = tf(K_f, [tau_a 1], 'InputDelay', L_a);

% --- MÉTODO C: ZIEGLER-NICHOLS (Tangente) ---
rango_inf = idx_s:idx_s+40;
df = diff(flujo(rango_inf))./diff(tiempo(rango_inf));
[m_max, idx_rel] = max(df);
idx_inf = idx_rel + idx_s - 1;

t_int_ini = (y_ini_f - flujo(idx_inf))/m_max + tiempo(idx_inf);
t_int_fin = (y_fin_f - flujo(idx_inf))/m_max + tiempo(idx_inf);

L_zn = max(0, t_int_ini - t_s);
tau_zn = max(0.1, t_int_fin - t_int_ini);
G_zn = tf(K_f, [tau_zn 1], 'InputDelay', L_zn);

%% 4. VALIDACIÓN GRÁFICA COMPLETA (Punto 4)
u_sim = V_control - 2; 
y_m = lsim(G_miller, u_sim, tiempo) + y_ini_f;
y_a = lsim(G_analitico, u_sim, tiempo) + y_ini_f;
y_z = lsim(G_zn, u_sim, tiempo) + y_ini_f;

figure('Name', 'Comparativa de Modelos de Flujo');
subplot(2,1,1);
plot(tiempo, V_control, 'b', 'LineWidth', 1.5); grid on;
ylabel('Voltaje [V]'); title('Señal de Entrada');

subplot(2,1,2); hold on;
plot(tiempo, flujo, 'k', 'LineWidth', 1, 'DisplayName', 'Datos Reales');
plot(tiempo, y_m, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Miller');
plot(tiempo, y_a, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Analítico');
plot(tiempo, y_z, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Ziegler-Nichols');
grid on; legend; ylabel('Flujo [l/min]'); xlabel('Tiempo [s]');
title(['Validación Global Flujo - K = ' num2str(K_f, 4)]);

