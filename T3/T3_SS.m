clc; clear; close all;

%--------------modelo intercambiador --------------------------------------
escala1= 1;
escala2= 1;
%modelo continuo
Gs = tf([0 2e-5],[1 0.004299 2e-5]);    % FUN 1
Gw = tf ([0 13e-5],[1 0.004299 2e-5]);  % FUN 2

[A, B, C, D] = ssdata([Gs Gw]);
sys_c = ss(A,B,C,D);

%step(sys_c);
autovalc= eig(A);       % determinarmos estabilidad 

%modelo discreto
Ts= 5; 
sys_d= c2d(sys_c,Ts);
G= sys_d.A;             %in discreta
H = sys_d.B;            % matriz de entrada
autovald= eig(H);       % estable ?

%-------------------- SIMULACION -------------------------------
tsim = 4000;            % tiempo de simulacion

N = round (tsim/Ts);    %n muestras por simulación
x(:,1) = [0;0];
y= zeros(1,N); 
u = [escala1*ones(1,N); escala2*ones(1,N)]; %señales de entrada
for k=1:N
   x(:,k+1) =  G*x(:,k)+H*u(:,k);
   y(k) =  C*x(:,k);
end

t = 0:Ts:(N-1)*Ts;      % vector de tiempo
subplot(2,1,1)

plot(t,y, 'LineWidth',2); grid on;
xlabel('Tiempo (s)'); ylabel('Temperatura (°C)');
title('Respuesta temporal')
subplot(2,1,2)

plot(x(1,:),x(2,:),'k-', 'LineWidth',1.5); grid on;
xlabel('x_1[k]'); ylabel('x_2[k]');
title('Espacio de estados')

%% MATRICES CALCULADAS FORMA CANONICA CONTROLABLE
escala1= 1;
escala2= 1;

A = [0 1;
    -2e-5 -0.004299];

B = [0 0;
     2e-5 13e-5];

C = [1 0];

D = [0 0];

sys = ss(A,B,C,D);
t = 0:1:5000;

u = [escala1*ones(length(t),1) escala2*ones(length(t),1)];

lsim(sys,u,t)

grid on
title('Escalon en valvula de vapor')
