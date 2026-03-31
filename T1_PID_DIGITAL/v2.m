clc; clear; close all;
%--------- DATOS DEL MOTOR Y MUESTREO ---------------
Km = 0.5;   % rad/(V*s) - ganancia del motor
Tm = 0.1;   % constante de tiempo del motor
T_sample = 0.01; % tiempo de muestreo (s)
tmax = 1.5; %tiempo  de simulacion
%--------- RESOLUCIONES ADC/DAC ---------------------
NADC = 10; % bits ADC
NDAC = 8;  % bits DAC
Vpp = 10;  % rango total +-5V

resADC = Vpp/2^NADC;   % resolución ADC
resDAC = Vpp/2^NDAC;   % resolución DAC
eqntADC= resADC/2;
eqntDAC= resDAC/2;

fprintf('--- 1. RESULTADOS A/D y D/A ---\n');
fprintf('Res. ADC (10 bits): %.2f mV | Error cuantización: ±%.2f mV\n',  resADC*1000, eqntADC*1000);
fprintf('Res. DAC (8 bits): %.2f mV | Error cuantización: ±%.2f mV\n\n', resDAC*1000, eqntDAC*1000);

%--------- FUNCION TRANSFERENCIA DEL MOTOR ----------
s = tf("s");
G_cont = Km/(s*(Tm*s + 1)); 
G_disc = c2d(G_cont, T_sample, "zoh"); 

%--------- PID (PARÁMETROS DE DISEÑO) ---------------
Mp = 0.05; 
ts = 460e-3; 
zeta = log(1/Mp)/(sqrt(pi^2+log(1/Mp)^2)); 
wn = 4/(ts*zeta);

gobj = tf([0 wn^2],[1  2*zeta*wn wn^2]); 
% polo lejano 
g1 = tf([0 10*wn*zeta],[1 10*wn*zeta]); 
g3orden = series(gobj,g1);

[num, den] = tfdata(g3orden, 'v'); 
ajuste_kd=0.35;
ajuste_ki=2.5;
Kp = (den(3) * Tm) / Km;
Ki =  ajuste_ki*(den(4) * Tm) / Km;
Kd =ajuste_kd*(den(2) * Tm - 1) / Km;
Ti= Kp/Ki;
Td= Kd/Kp;
fprintf('Ganancias analógicas calculadas: Kp=%.2f, Ti=%.2f, Td=%.2f\n', Kp, Ti, Td);

%GRAFICA comparativa 
% figure;
% step(gobj,tmax);hold on;  step(g3orden,tmax);
% legend('2do orden','3er orden')
% title('Respuesta al escalón de sistemas orden 2 y 3');

%--------- SIMULACIÓN DISCRETA (ECUACIONES EN DIFERENCIAS) ----------
%AUMENTAR  A 12V PARA CUMPLIR
u_max = 12; % saturación a +-5V

t = 0:T_sample:tmax; % horizonte de simulación
r =  1*ones(size(t));  % Referencia: escalón unitario

% Inicialización de vectores
u      = zeros(size(t)); % salida del PID ideal
u_dac  = zeros(size(t)); % salida cuantizada real aplicada a la planta
y_real = zeros(size(t)); % salida física del motor
y_adc  = zeros(size(t)); % lectura del microcontrolador
e      = zeros(size(t)); % error

% Extraer coeficientes exactos del modelo discreto (ZOH)
[numG, denG] = tfdata(G_disc, 'v');

% Bucle paso a paso emulando el microcontrolador
for k = 3:length(t)
    % 1. El motor físico evoluciona basándose en voltajes pasados
    y_real(k) = -denG(2)*y_real(k-1) - denG(3)*y_real(k-2) + numG(2)*u_dac(k-1) + numG(3)*u_dac(k-2);

    % 2. El sensor mide, satura en los límites físicos y el ADC cuantiza
    % CORRECCIÓN 3: Saturación FÍSICA de la entrada del ADC
    y_sensor = max(min(y_real(k), u_max), -u_max); 
    y_adc(k) = round(y_sensor/resADC)*resADC;

    % 3. El micro calcula el error actual
    e(k) = r(k) - y_adc(k);

    % 4. Ley de control PID (Forma incremental )
    u(k) = u(k-1) + Kp*((e(k)-e(k-1)) + (T_sample/Ti)*e(k) + (Td/T_sample)*(e(k) - 2*e(k-1) + e(k-2)));

    % 5. Saturación matemática en el código del micro
    u(k) = max(min(u(k), u_max), -u_max);

    % 6. Salida física por el DAC
    u_dac(k) = round(u(k)/resDAC)*resDAC;
end

%--------- GRÁFICAS DE ANÁLISIS DE CUANTIZACIÓN ----------------------
figure('Name', 'Análisis de Cuantización: ADC y DAC', 'Position', [100, 100, 1000, 600]);

% --- 1. ADC: Señales Superpuestas ---
subplot(2,2,1);
plot(t, y_real, 'b', 'LineWidth', 1.5); hold on;
stairs(t, y_adc, 'r', 'LineWidth', 1.2);
yline(max(r), 'k:','Referencia','LineWidth', 1.2 );
step(g3orden,tmax);
title('Sensor: Lectura del ADC (10 bits)');
xlabel('Tiempo (s)'); ylabel('Velocidad (rad/s)');
legend('y_{real} (Física)', 'y_{adc} (Digitalizada)','Respuesta esperada', 'Location', 'Southeast');
grid on;

% --- 2. ADC: Error de Cuantización ---
subplot(2,2,2);
plot(t, y_real - y_adc, 'k', 'LineWidth', 1.2); hold on;
yline(resADC/2, 'b--', 'Límite +\Delta/2'); 
yline(-resADC/2, 'b--', 'Límite -\Delta/2');
title('Error de Cuantización en el ADC');
xlabel('Tiempo (s)'); ylabel('Error (rad/s)');
ylim([-resADC resADC]);
grid on;

% --- 3. DAC: Señales Superpuestas ---
subplot(2,2,3);
stairs(t, u, 'r', 'LineWidth', 1.5); hold on;
stairs(t, u_dac, 'b', 'LineWidth', 1.5);
yline(5, 'k:', '+5V'); yline(-5, 'k:', '-5V');
title('Actuador: Salida del DAC (8 bits)');
xlabel('Tiempo (s)'); ylabel('Voltaje (V)');
legend('u (PID Matemático)', 'u_{dac} (Voltaje Real)', 'Location', 'Southeast');
%ylim([-6 6]);
grid on;

% --- 4. DAC: Error de Cuantización ---
subplot(2,2,4);
stairs(t, u - u_dac, 'k', 'LineWidth', 1.2); hold on;
yline(resDAC/2, 'b--', 'Límite +\Delta/2'); 
yline(-resDAC/2, 'b--', 'Límite -\Delta/2');
title('Error de Cuantización en el DAC');
xlabel('Tiempo (s)'); ylabel('Error (V)');
ylim([-resDAC resDAC]);
grid on;