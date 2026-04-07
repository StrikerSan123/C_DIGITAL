%% 5. CONTROL Y VALIDACIÓN DE POSICIÓN (Punto 5)
% 1. Identificación Miller para posición
yp_ini = mean(posicion(idx_s-15:idx_s));
yp_fin = mean(posicion(idx_s+55:idx_s+75));
dyp = yp_fin - yp_ini; 
Kp = dyp/du;

vp28 = yp_ini + 0.283*dyp; vp63 = yp_ini + 0.632*dyp;
tp28 = tiempo(find(posicion(idx_s:end) >= vp28, 1) + idx_s - 1);
tp63 = tiempo(find(posicion(idx_s:end) >= v63, 1) + idx_s - 1); % Usamos v63 de pos
tau_p = 1.5 * (tp63 - tp28);
L_p = max(0, (tp63 - t_s) - tau_p);

% Forzamos una dinámica mínima si el sistema es demasiado rápido para ver la curva
if tau_p < 0.5, tau_p = 1.2; end 
Gp = tf(Kp, [tau_p 1], 'InputDelay', L_p);

% 2. Diseño de Controlador PI
Ti = tau_p; 
Kc = (0.5 * tau_p) / (Kp * (L_p + 0.5)); 
C = pid(Kc, Kc/Ti);

% 3. Lazo Cerrado (Sistema Realimentado)
T_cl = feedback(C*Gp, 1);

% --- VALIDACIONES ---

% A. Validación Lazo Abierto (El modelo frente a la posición real)
y_open_loop = lsim(Gp, u_sim, tiempo) + yp_ini;

% B. Simulación Lazo Cerrado (Cómo seguiría el controlador esa misma señal)
% Usamos V_control como la REFERENCIA (Setpoint dinámico)
% Restamos el offset para que empiece en 0 y sumamos el nivel inicial al final
y_closed_loop = lsim(T_cl, u_sim, tiempo) + yp_ini;

figure('Name', 'Análisis de Posición: Lazo Abierto vs Lazo Cerrado');

subplot(2,1,1); hold on;
plot(tiempo, posicion, 'k', 'LineWidth', 1, 'DisplayName', 'Posición Real (Datos)');
plot(tiempo, y_open_loop, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Modelo Identificado');
grid on; legend; ylabel('Posición [rad]');
title('Validación de Identificación (Lazo Abierto)');

subplot(2,1,2); hold on;
plot(tiempo, V_control, 'r:', 'LineWidth', 1, 'DisplayName', 'Referencia (V control)');
plot(tiempo, y_closed_loop, 'm', 'LineWidth', 1.5, 'DisplayName', 'Respuesta PI');
grid on; legend; ylabel('Posición [rad]'); xlabel('Tiempo [s]');
title('Simulación de Control: Siguiendo la secuencia del archivo');

fprintf('Modelo Posición: K=%.3f, tau=%.3f, L=%.3f\n', Kp, tau_p, L_p);