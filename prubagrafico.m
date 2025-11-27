clc
clear
close all

%% Parámetros del sistema Rössler
a = 0.2; b = 0.2; c = 5.7;

% Condiciones iniciales
xi_0 = [0.01, -0.03, 0.04];
xi_1 = [6, -28, 28];
xi_2 = [1 1 1];

% Ventana temporal
t_span = linspace(0, 200, 5000);

% Sistema de Rössler
rossler = @(t,x) [-x(2) - x(3);
                  x(1) + a*x(2);
                  b + x(3)*(x(1) - c)];

% Resolución
[t, x] = ode45(rossler, t_span, xi_0);
[~, y] = ode45(rossler, t_span, xi_1);
[~, z] = ode45(rossler, t_span, xi_2);

%% Figura mejorada
figure('Color','white');
plot3(x(:, 1), x(:, 2), x(:, 3), 'g-', 'LineWidth', 1.5); hold on;
plot3(y(:, 1), y(:, 2), y(:, 3), 'b--', 'LineWidth', 1.5);
plot3(z(:, 1), z(:, 2), z(:, 3), 'r-', 'LineWidth', 1.5);

% Puntos iniciales
plot3(xi_0(1), xi_0(2), xi_0(3), 'r*', 'MarkerSize',10,'MarkerFaceColor', 'g');
plot3(xi_1(1), xi_1(2), xi_1(3), 'go', 'MarkerFaceColor', 'b');
plot3(xi_2(1), xi_2(2), xi_2(3), 'bo', 'MarkerFaceColor', 'r');

xlabel('x'); ylabel('y'); zlabel('z');
title('Diagrama de Fase del Sistema Rössler');
legend('Trayectoria 1', 'Trayectoria 2', 'Trayectoria 3', 'Inicio 1', 'Inicio 2', 'Inicio 3');
grid on; axis tight;
view(3);

% Mejoras visuales
set(gca,'FontSize',12,'LineWidth',1.2);
colormap(parula); % Colores atractivos
camlight headlight; lighting gouraud; % Iluminación suave
hold off;

%% Configuración del video
v = VideoWriter('Rossler3D_Mejorado.mp4','MPEG-4');
v.FrameRate = 30;
open(v);

frames = 360; % Número de pasos
for k = 1:frames
    az = k; % Rotación horizontal
    el = 30 + 15*sin(k*pi/180); % Oscilación vertical suave
    view(az, el);
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
disp('✅ Video guardado como Rossler3D_Mejorado.mp4');